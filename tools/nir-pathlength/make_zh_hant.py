# 把引擎的簡體中文介面字串轉成繁體（台灣），寫成 lib/zh-hant.js（引擎更新後重跑）
import re, json, os
from opencc import OpenCC
here=os.path.dirname(os.path.abspath(__file__)); lib=os.path.join(here,'..','..','site','nir-pathlength','lib')
src=open(os.path.join(lib,'engine.js'),encoding='utf-8').read()
blk=re.search(r"\n  zh: \{(.*?)\n  \},", src, re.S).group(1)
cc=OpenCC('s2twp'); pairs=dict(re.findall(r"(\w+): '([^']*)'", blk))
out={k:cc.convert(v) for k,v in pairs.items()}
js=("// 由 tools/nir-pathlength/make_zh_hant.py 產生：引擎介面改為繁體中文\n"
    "Object.assign(TRANSLATIONS.zh, "+json.dumps(out,ensure_ascii=False,indent=0)+");\n"
    "document.querySelector('#stage')?.setAttribute('aria-label', translate('stageLabel'));\n"
    "{ const pm = document.querySelector('.paused-mark div'); if (pm) pm.textContent = translate('paused'); }\n")
# 引擎在載入時就已畫好封面按鈕、快捷鍵說明等（早於上面的覆寫）：用逐字「簡→繁」對照再轉一次頁面文字，並持續監看之後加入的節點。
cc1 = OpenCC('s2t'); emap = {}
for v in pairs.values():
    for a in v:
        b = cc1.convert(a)
        if len(b) == 1 and b != a: emap[a] = b
emap.update({'钟': '鐘'})
js += ("(() => {\n  const MAP = " + json.dumps(emap, ensure_ascii=False) + ";\n"
       "  const RE = new RegExp('[' + Object.keys(MAP).join('') + ']', 'g'), conv = s => s.replace(RE, c => MAP[c]);\n"
       "  const fix = root => { if (!root) return; if (root.nodeType === 3) { const t = conv(root.nodeValue); if (t !== root.nodeValue) root.nodeValue = t; return; }\n"
       "    const w = document.createTreeWalker(root, NodeFilter.SHOW_TEXT);\n"
       "    for (let n; (n = w.nextNode());) { const t = conv(n.nodeValue); if (t !== n.nodeValue) n.nodeValue = t; }\n"
       "    (root.querySelectorAll ? root.querySelectorAll('[aria-label],[title]') : []).forEach(e => ['aria-label', 'title'].forEach(a => {\n"
       "      const v = e.getAttribute(a); if (v) { const t = conv(v); if (t !== v) e.setAttribute(a, t); } })); };\n"
       "  const run = () => { fix(document.body); new MutationObserver(ms => ms.forEach(m => m.type === 'characterData' ? fix(m.target) : m.addedNodes.forEach(fix))).observe(document.body, { childList: true, subtree: true, characterData: true }); };\n"
       "  document.body ? run() : addEventListener('DOMContentLoaded', run);\n})();\n")
open(os.path.join(lib,'zh-hant.js'),'w',encoding='utf-8').write(js); print(len(out),'strings'); print(out['startLesson'],out['aboutMinutes'],out['quickCheck'])

# 封面／目錄（contents.js 的字串在閉包裡，無法覆寫）：產生「簡→繁」逐字對照，在頁面渲染後轉換文字
csrc=open(os.path.join(lib,'contents.js'),encoding='utf-8').read()
cblk=re.search(r"\n    zh: \{(.*?)\n    \},", csrc, re.S).group(1)
cc2=OpenCC('s2t'); cmap={}
for v in re.findall(r"'([^']*)'", cblk):
    for a in v:   # 逐字轉換（不用詞組轉換，以免把「看」之類兩岸通用字也改掉）
        b=cc2.convert(a)
        if len(b)==1 and b!=a: cmap[a]=b
cmap.update({'钟': '鐘'})   # 分钟 → 分鐘（逐字轉換會給「鍾」）
cjs=("// 由 tools/nir-pathlength/make_zh_hant.py 產生：封面與目錄的介面文字改為繁體\n"
     "(() => {\n  const MAP = "+json.dumps(cmap,ensure_ascii=False)+";\n"
     "  const RE = new RegExp('[' + Object.keys(MAP).join('') + ']', 'g'), conv = s => s.replace(RE, c => MAP[c]);\n"
     "  const fix = root => {\n"
     "    const w = document.createTreeWalker(root, NodeFilter.SHOW_TEXT);\n"
     "    for (let n; (n = w.nextNode());) { const t = conv(n.nodeValue); if (t !== n.nodeValue) n.nodeValue = t; }\n"
     "    (root.querySelectorAll ? root.querySelectorAll('[aria-label],[title]') : []).forEach(e => ['aria-label', 'title'].forEach(a => {\n"
     "      const v = e.getAttribute(a); if (v) { const t = conv(v); if (t !== v) e.setAttribute(a, t); } }));\n"
     "  };\n"
     "  const run = () => { fix(document.body); new MutationObserver(ms => ms.forEach(m => m.type === 'characterData' ? fix(m.target.parentNode || document.body) : m.addedNodes.forEach(n => n.nodeType === 3 ? fix(n.parentNode) : n.nodeType === 1 && fix(n)))).observe(document.body, { childList: true, subtree: true, characterData: true }); };\n"
     "  document.readyState === 'loading' ? addEventListener('DOMContentLoaded', run) : run();\n})();\n")
open(os.path.join(lib,'zh-hant-contents.js'),'w',encoding='utf-8').write(cjs); print(len(cmap),'chars mapped for contents')
