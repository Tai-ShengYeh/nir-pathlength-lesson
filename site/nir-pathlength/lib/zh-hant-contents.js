// 由 tools/nir-pathlength/make_zh_hant.py 產生：封面與目錄的介面文字改為繁體
(() => {
  const MAP = {"开": "開", "这": "這", "书": "書", "继": "繼", "续": "續", "单": "單", "钟": "鐘", "个": "個", "种": "種", "风": "風", "参": "參", "献": "獻", "录": "錄", "课": "課", "链": "鏈", "会": "會", "标": "標", "签": "籤", "页": "頁", "该": "該", "图": "圖", "馆": "館", "记": "記", "没": "沒", "线": "線", "来": "來"};
  const RE = new RegExp('[' + Object.keys(MAP).join('') + ']', 'g'), conv = s => s.replace(RE, c => MAP[c]);
  const fix = root => {
    const w = document.createTreeWalker(root, NodeFilter.SHOW_TEXT);
    for (let n; (n = w.nextNode());) { const t = conv(n.nodeValue); if (t !== n.nodeValue) n.nodeValue = t; }
    (root.querySelectorAll ? root.querySelectorAll('[aria-label],[title]') : []).forEach(e => ['aria-label', 'title'].forEach(a => {
      const v = e.getAttribute(a); if (v) { const t = conv(v); if (t !== v) e.setAttribute(a, t); } }));
  };
  const run = () => { fix(document.body); new MutationObserver(ms => ms.forEach(m => m.type === 'characterData' ? fix(m.target.parentNode || document.body) : m.addedNodes.forEach(n => n.nodeType === 3 ? fix(n.parentNode) : n.nodeType === 1 && fix(n)))).observe(document.body, { childList: true, subtree: true, characterData: true }); };
  document.readyState === 'loading' ? addEventListener('DOMContentLoaded', run) : run();
})();
