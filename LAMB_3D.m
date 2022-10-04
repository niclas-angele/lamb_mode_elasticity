% clear; close all; 
addpath('ffmatlib');
figure
[p,b,t,nv,nbe,nt,labels]=ffreadmesh('exp1.mesh');
[vh]=ffreaddata('exp2.txt');
[u]=ffreaddata('exp4.txt');


% ffpdeplot3D(p,b,t,'VhSeq',vh,'XYZData',u,'Slice',[-1,-1,0.2],[1,-1,0.2],[-1,1,0.2],'Project2D','on','Mesh','off','Colormap','parula');
% ffpdeplot3D(p,b,t,'VhSeq',vh,'XYZData',u,'Mesh','off','Colormap','parula');
ffpdeplot3D(p,b,t,'VhSeq',vh,'XYZData',u,'Slice',[-1,-1,-0.2; -1,-1,-0.2;-1,-1,0.2],[-1,-1,0.2;-1,-1,0.2;-1,1,0.2],[1,-1,-0.2; -1,1,-0.2;1,-1,0.2],'Mesh','off','Colormap','parula');
axis([-1,1,-1,1,-0.2,0.2])
% % colorbar; 
plot3([-1,-1],[-1,1],[0.2,0.2],'k')
plot3([-1,1],[-1,-1],[0.2,0.2],'k')
plot3([-1,-1],[-1,-1],[0.2,-0.2],'k')
view(-30,30)
% caxis([minz,maxz]); 
