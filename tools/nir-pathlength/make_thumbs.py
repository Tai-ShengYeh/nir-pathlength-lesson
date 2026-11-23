# 課程入口頁的章節縮圖：每章取一個關鍵畫面（該 beat 結尾），存成 portal/img/chNN.jpg
# 需先啟動 python tools/serve.py 8765 site
import asyncio, io, os
from playwright.async_api import async_playwright
from PIL import Image
SHOTS = {'ch01': 'saturate', 'ch02': 'compare', 'ch03': 'shuffle', 'ch04': 'coef', 'ch05': 'four', 'ch06': 'check', 'ch07': 'order', 'ch08': 'agree'}
OUT = os.path.join(os.path.dirname(__file__), '..', '..', 'portal', 'img')
async def main():
    async with async_playwright() as p:
        b = await p.chromium.launch()
        pg = await b.new_page(viewport={'width': 1600, 'height': 960}, reduced_motion='reduce')
        for ch, beat in SHOTS.items():
            await pg.goto(f'http://localhost:8765/nir-pathlength/{ch}/'); await pg.wait_for_function('typeof seek === "function" && typeof BEATS !== "undefined"')
            await pg.evaluate("""(id) => { const i = BEATS.findIndex(b => b.id === id); document.getElementById('cover')?.remove(); seek(i, false); P.t = P.end; evalTo(P.end); }""", beat)
            await pg.wait_for_timeout(400)
            r = await pg.evaluate("(() => { const r = document.getElementById('stage').getBoundingClientRect(); return [r.x, r.y, r.width, r.height]; })()")
            im = Image.open(io.BytesIO(await pg.screenshot(clip={'x': r[0], 'y': r[1], 'width': r[2], 'height': r[3]}))).convert('RGB')
            im = im.resize((800, round(800 * im.height / im.width)), Image.LANCZOS)
            im.save(os.path.join(OUT, ch + '.jpg'), quality=84); print(ch, beat, im.size)
        await b.close()
asyncio.run(main())
