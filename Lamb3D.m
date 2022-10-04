clear; close all; 

% parameters 
ct=0.5;
cl=0.9;
mu=ct^2; 
lambda=cl^2-2*ct^2;
w=10; 
h=0.2;
x=linspace(-1,1,100); 
y=linspace(-1,1,100); 
z=linspace(-h,h);
[X,Y]=meshgrid(x,y); 
[t1,t2]=meshgrid(y,z);

%sources 
sig=0.1; 
gaus=@(x,y,z) 1/sig/sig/2/pi*exp(-(x.^2+y.^2)/2/sig/sig).*z; 

%%% Figure 7: 
% fz=@(x,y,z) gaus(x,y,z);
% dfz1=@(x,y,z) -x/sig^2.*gaus(x,y,z);
% dfz2=@(x,y,z) -y/sig^2.*gaus(x,y,z);
% fl1=@(x,y,z) -x.*gaus(x,y,z); 
% fl2=@(x,y,z) -y.*gaus(x,y,z); 
% divfl=@(x,y,z) ((x.^2+y.^2)/(sig^2)-2).*gaus(x,y,z); 
% bsht1=@(x,y) 0*x;
% bsht2=@(x,y) 0*x;

%%%Figure 8:
bsht1=@(x,y) -y/sig/sig/4/pi.*exp(-(x.^2+y.^2)/2/sig/sig);
bsht2=@(x,y) x/sig/sig/4/pi.*exp(-(x.^2+y.^2)/2/sig/sig);
fz=@(x,y,z) x*0;
dfz1=@(x,y,z) x*0;
dfz2=@(x,y,z) x*0;
fl1=@(x,y,z) x*0; 
fl2=@(x,y,z) y*0; 
divfl=@(x,y,z) x*0; 

%%% Figures 7 and 8:
fsh1=@(x,y,z) 0*x; 
fsh2=@(x,y,z) 0*x; 
blt1=@(x,y) 0*x; 
blt2=@(x,y) 0*x; 
blb1=@(x,y) 0*x; 
blb2=@(x,y) 0*x; 
bzt=@(x,y) 0*x; 
bzb=@(x,y) 0*x; 
dbzt1=@(x,y) 0*x; 
dbzt2=@(x,y) 0*x; 
dbzb1=@(x,y) 0*x; 
dbzb2=@(x,y) 0*x; 
divblt=@(x,y) 0*x; 
divblb=@(x,y) 0*x;
bshb1=@(x,y) 0*x; 
bshb2=@(x,y) 0*x; 



fl=@(x,y,z) [fl1(x,y,z);fl2(x,y,z)];
fsh=@(x,y,z) [fsh1(x,y,z);fsh2(x,y,z)];
blt=@(x,y) [blt1(x,y);blt2(x,y)];
blb=@(x,y) [blb1(x,y);blb2(x,y)];
gradfz=@(x,y,z) [dfz1(x,y,z);dfz2(x,y,z)];
gradbzt=@(x,y) [dbzt1(x,y);dbzt2(x,y)];
gradbzb=@(x,y) [dbzb1(x,y);dbzb2(x,y)];
bsht=@(x,y) [bsht1(x,y);bsht2(x,y)];
bshb=@(x,y) [bshb1(x,y);bshb2(x,y)];

%solution obtained with convolution method:
[UX,UY,UZ]=solveLamb3D(w,fz,fl,gradfz,divfl,blt,bzt,gradbzt,divblt,blb,bzb,gradbzb,divblb, fsh,bsht,bshb, x,y,h);


