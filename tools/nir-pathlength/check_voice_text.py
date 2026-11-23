# 檢查旁白文字和產生的音檔字幕是否一致（忽略斷字空格）。
# papermorph 的 lint.py 會把 {{量|涼}} 這種中文讀音改寫切成獨立的詞，因此誤報「narration text differs」；用這支確認真的一致。
# 用法：python tools/nir-pathlength/check_voice_text.py SKILL_SCRIPTS_DIR ch01 ch02 ...
import json, os, re, sys
sys.path.insert(0, sys.argv[1]); import lint
root = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..')
bad = 0
for ch in sys.argv[2:]:
    js = open(os.path.join(root, 'site', 'nir-pathlength', ch, 'audio', 'zh', 'timings.js'), encoding='utf-8').read()
    T = json.loads(re.search(r'=\s*(\{.*\})', js, re.S).group(1))
    N = json.load(open(os.path.join(root, 'content', 'nir-pathlength', ch, 'narration.zh.json'), encoding='utf-8'))['beats']
    miss = [b for b, t in N.items() if b not in T or ''.join(lint.words(t)) != ''.join(lint.words(' '.join(c[1] for c in T[b]['cues'])))]
    bad += len(miss); print(ch, 'ok' if not miss else 'differs: ' + ', '.join(miss))
sys.exit(1 if bad else 0)
