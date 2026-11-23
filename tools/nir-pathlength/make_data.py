# 由 teaching/data/spectra_long.csv 產生課程動畫用的 lib/nirdata.js（數值與 analysis/pathlength_select.py 同一算法）
import numpy as np, pandas as pd, json, os
from scipy.signal import savgol_filter
here=os.path.dirname(os.path.abspath(__file__)); root=os.path.join(here,'..','..')
d=pd.read_csv(os.path.join(root,'data','spectra_long.csv')); d=d[d.replicate==1]
paths=[0.25,0.5,1,2,5]; conc=[0,3.32,4.87,7.45]
wl=np.sort(d.wavelength_nm.unique())
A={}; D2={}; out={'wl':[round(w,1) for w in wl],'paths':paths,'conc':conc,'A':{},'band':{}}
for p in paths:
    S=np.array([d[(d.path_mm==p)&(d.acid_pct==c)].sort_values('wavelength_nm').absorbance.values for c in conc])
    out['A'][str(p)]=[[round(v,4) for v in s] for s in S]
    d2=savgol_filter(S,9,2,deriv=2,axis=1)
    for name,w0 in [('2258',2258),('1670',1670)]:
        i=int(np.argmin(abs(wl-w0))); sig=d2[:,i]
        sl,ic=np.polyfit(conc,sig,1); pred=np.polyval([sl,ic],conc)
        r2=1-np.sum((sig-pred)**2)/np.sum((sig-sig.mean())**2)
        out['band'].setdefault(name,{})[str(p)]={'i':i,'d2':[float('%.4g'%v) for v in sig],'slope':float('%.4g'%sl),'icpt':float('%.4g'%ic),'r2':round(r2,4)}
    if p==1: out['d2_1mm']=[[float('%.4g'%v) for v in s] for s in d2]
d_all=pd.read_csv(os.path.join(root,'data','spectra_long.csv'))
r2=d_all[(d_all.path_mm==0.25)&(d_all['sample']=='water')&(d_all.replicate==2)].sort_values('wavelength_nm').absorbance.values
out['A_water025_rep2']=[round(v,4) for v in r2]
with open(os.path.join(root,'site','nir-pathlength','lib','nirdata.js'),'w',encoding='utf-8') as f:
    f.write('// 由 tools/nir-pathlength/make_data.py 產生，請勿手改\nconst NIR = '+json.dumps(out,ensure_ascii=False,separators=(',',':'))+';\n')
b=out['band']['2258']; print({p:(b[str(p)]['slope'],round(b[str(p)]['slope']/b['0.25']['slope'],2),b[str(p)]['r2']) for p in paths})
