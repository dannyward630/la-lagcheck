"""Recalculate the local feature measurements used in Figure 3d."""
import numpy as np
from scipy import ndimage
import reconstruction as r
qc, ql = r.qc, r.ql
extent = list(r.ext)
dx = (extent[1]-extent[0])/qc.shape[1]
dy = (extent[2]-extent[3])/qc.shape[0]
x = extent[0]+(np.arange(qc.shape[1])+.5)*dx
y = extent[3]+(np.arange(qc.shape[0])+.5)*dy
anchor = int(np.argmin(abs((r.ys-r.uy.min())/1000-5.2)))
shift = -float(r.disp_um[anchor]/1000)
roi = [13.6,15.0,4.7,6.2]
def measure(img, xx, threshold):
    region=(xx[None,:]>=roi[0])&(xx[None,:]<=roi[1])&(y[:,None]>=roi[2])&(y[:,None]<=roi[3])
    mask=region&np.isfinite(img)&(img>=threshold)
    labels,n=ndimage.label(mask,structure=np.ones((3,3)))
    sizes=np.bincount(labels.ravel()); sizes[0]=0
    lab=int(np.argmax(sizes)); rr,cc=np.where(labels==lab)
    points=np.column_stack([xx[cc],y[rr]])
    center=points.mean(axis=0)
    values,vectors=np.linalg.eigh(np.cov(points.T))
    axis=vectors[:,-1]
    if axis[1]<0: axis=-axis
    angle=float(np.degrees(np.arctan2(axis[0],axis[1])))
    assert points[:,0].min()>roi[0]+dx and points[:,0].max()<roi[1]-dx
    assert points[:,1].min()>roi[2]+dy and points[:,1].max()<roi[3]-dy
    return dict(n_pixels=len(rr),center=center.tolist(),axis=axis.tolist(),angle_deg=angle,
                eigenvalues=values.tolist(),threshold=threshold),labels==lab

comparisons=[]
for threshold in [4.9,5.0,5.1]:
    a,_=measure(qc,x+shift,threshold)
    b,_=measure(ql,x,threshold)
    comparisons.append(dict(threshold=threshold,constant=a,mapped=b,difference_deg=a['angle_deg']-b['angle_deg']))
audit=dict(status='post_hoc_feature_comparison', rate_ppm=r.R_IT*1e6,
           roi_mm=roi, display_window_mm=[13.35,15.25,4.4,6.5],
           constant_translation_mm=shift, anchor_line_zero_based=anchor,
           comparisons=comparisons,
           scope='A feature selected after reconstruction; neither global shear nor an independently measured anatomical direction.')
