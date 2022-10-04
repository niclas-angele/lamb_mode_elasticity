clear; close all;

%parameters of the plate
ct=0.5;
cl=0.9;
mu=ct^2; 
lambda=cl^2-2*ct^2;
h=0.1; 

%frequency range
W=linspace(0.1,17,170);

%chosen portion of the waveguide 
x=linspace(-1,5,300);
y=linspace(-h,h,30);
[X,Y]=meshgrid(x,y);

%%% Figure 9, right reconstruction 
% amp=50; 
% supp=3;
% g1=@(x) amp*5/4/4*(x>3.7).*(x<4.2).*(x-3.7).^2.*(x-4.2).^2*2*h;
% g1prim=@(x) amp*5/2/4*(x>3.7).*(x<4.2).*(x-3.7).*(x-4.2).*(2*x-7.9)*2*h;
% g2=@(x) amp*5/4/4*(x>3.4).*(x<4).*(x-3.4).^2.*(x-4).^2*2*h;
% g2prim=@(x) amp*5/2/4*(x>3.4).*(x<4).*(x-3.4).*(x-4).*(2*x-7.4)*2*h;

%%% Figure 9, left reconstruction 
amp=1; %amplitude 
supp=3; %left of the support 
g1=@(x) amp*5/4/4*(x>3.2).*(x<4.2).*(x-3.2).^2.*(x-4.2).^2*2*h;
g1prim=@(x) amp*5/2/4*(x>3.2).*(x<4.2).*(x-3.2).*(x-4.2).*(2*x-7.4)*2*h;
g2=@(x) -7*amp*5/4/4*(x>3.4).*(x<4).*(x-3.4).^2.*(x-4).^2*2*h;
g2prim=@(x) -7*amp*5/2/4*(x>3.4).*(x<4).*(x-3.4).*(x-4).*(2*x-7.4)*2*h;

%%% Figure 10, with changes in amp 
% amp=0.5; 
% supp=2.5; 
% g1=@(x) amp*(x>3).*(x<5).*(x-3).^2.*(x-5).^2*2*h;
% g1prim=@(x) 2*amp*(x>3).*(x<5).*(x-3).*(x-5).*(2*x-8)*2*h;
% g2=@(x) 0*x;
% g2prim=@(x) 0*x;


%%%% résolution FEM
nodes = [-8 -0.1;11 -0.1;11 0.1;-8 0.1];
edges = {{1,2},{2,3},{3,4},{4,1}};
dx = 0.1; 
domain = Domain(nodes,edges);
mesh = Mesh(domain,dx);
mesh=mesh.submesh;
mesh=mesh.submesh;

%cut here for a fast run of the program 
mesh=mesh.submesh;
mesh=mesh.submesh;
mesh=mesh.submesh;
mesh=mesh.submesh;

Xi1=W*0;  % storage of 2k_1 
Xi2=W*0;  %storage of k_1+k_2
% storage of the corresponding fourier transform
F1=W*0; 
F2=W*0; 


for j=1:length(W)
    j
    w=W(j); 
    %computation of the wavenumbers associated to w
    [ks,ka]=dispersionreal(w,5*w,h);
    k1=ks(1); 
    p=@(k) sqrt(w.^2./(cl^2)-k.^2);
    q=@(k) sqrt(w.^2./(ct^2)-k.^2);
    %symetric modes
    jns=@(k) 1i*mu*k.*(q(k).^2+k.^2).*(h*(q(k).^2-k.^2).^2.*sin(q(k)*h).^2+h*4*k.^2.*p(k).^2.*sin(p(k)*h).^2+sin(p(k)*h).*sin(q(k)*h).*(1./p(k).*(q(k).^2-k.^2).*(q(k).^2-k.^2-8*p(k).^2).*sin(q(k)*h).*cos(p(k)*h)+4./q(k).*p(k).^2.*(2*q(k).^2-k.^2).*sin(p(k)*h).*cos(q(k)*h)));
    us=@(k,y) (1i*k*(q(k)^2-k^2)*sin(q(k)*h)*cos(p(k)*y)-2*1i*k*p(k)*q(k)*sin(p(k)*h)*cos(q(k)*y)); 
    vs=@(k,y) (-p(k)*(q(k)^2-k^2)*sin(q(k)*h)*sin(p(k)*y)-2*k^2*p(k)*sin(p(k)*h)*sin(q(k)*y));
    ss=@(k,y) -(q(k)^2-k^2)*(cl^2*k^2+lambda*p(k)^2)*sin(q(k)*h)*cos(p(k)*y)+4*mu*p(k)*q(k)*k^2*sin(p(k)*h)*cos(q(k)*y);
    ts=@(k,y) 2*1i*k*mu*(q(k)^2-k^2)*p(k)*(-sin(q(k)*h)*sin(p(k)*y)+sin(p(k)*h)*sin(q(k)*y));
    %antisymetric modes 
    ua=@(k,y) 1i*k*(q(k)^2-k^2)*cos(q(k)*h)*sin(p(k)*y)-2*1i*k*p(k)*q(k)*cos(p(k)*h)*sin(q(k)*y);
    va=@(k,y) p(k)*(q(k)^2-k^2)*cos(q(k)*h)*cos(p(k)*y)+2*k^2*p(k)*cos(p(k)*h)*cos(q(k)*y);
    sa=@(k,y) -(q(k)^2-k^2)*(cl^2*k^2+lambda*p(k)^2)*cos(q(k)*h)*sin(p(k)*y)+4*mu*p(k)*q(k)*k^2*cos(p(k)*h)*sin(q(k)*y);
    ta=@(k,y) 2*1i*k*mu*(q(k)^2-k^2)*p(k)*(cos(q(k)*h)*cos(p(k)*y)-cos(p(k)*h)*cos(q(k)*y));
    jna=@(k) 1i*mu*k.*(q(k).^2+k.^2).*(h*(q(k).^2-k.^2).^2.*cos(q(k)*h).^2+h*4*k.^2.*p(k).^2.*cos(p(k)*h).^2-cos(p(k)*h).*cos(q(k)*h).*(1./p(k).*(q(k).^2-k.^2).*(q(k).^2-k.^2-8*p(k).^2).*cos(q(k)*h).*sin(p(k)*h)+4./q(k).*p(k).^2.*(2*q(k).^2-k.^2).*cos(p(k)*h).*sin(q(k)*h)));
    %source terms
    f1=@(x,y) 0*x; 
    f2=@(x,y) 0*x;
    b1top=@(x,y) 2*h*g1prim(x).*ss(k1,h).*exp(1i*k1*x); 
    b2top=@(x,y) 0*x; 
    b1bot=@(x,y) -2*h*g2prim(x).*ss(k1,-h).*exp(1i*k1*x); 
    b2bot=@(x,y) 0*x; 
    %boundary conditions
    bc1 = {1,1,0,{b1bot,b2bot}};
    bc3 = {3,1,0,{b1top,b2top}};
    bc2 = {[2,4],0,1,0};
    %PML
    l1=6; 
    l2=-2; 
    at=10; 
    alpha=@(x,y) -1i*w./(-1i*w+at*(x-l1)).*(x>=l1)-1i*w./(-1i*w+at*(-x+l2)).*(x<=l2)+1*(x<l1).*(x>l2);
    C=mesh.P0({@(x,y) (lambda+2*mu)*alpha(x,y), 0, 0 , mu, 0 , lambda, @(x,y) mu*alpha(x,y) , 0 , 0 , mu, @(x,y) lambda*alpha(x,y), 0 , @(x,y) mu*alpha(x,y), 0, 0, lambda+2*mu}); 
    %solution
    u = mesh.solve_vect(C,mesh.P0(@(x,y) -w.^2),{f1,f2},{bc1,bc2,bc3});
    %projection on the observed area 
    temp1=P1togrid(mesh,u(:,1),x,y);
    temp2=P1togrid(mesh,u(:,2),x,y);
    %correction of the solution if we have backward modes 
    for ell=1:length(ks)
        if ks(ell)<0
            [temp1,temp2]=correction_PML_droit(u,mesh,w,h,x,y,ks(ell),l1,l2);
        end
    end
    for ell=1:length(ka)
        if ka(ell)<0
            [temp1,temp2]=correction_PML_droit(u,mesh,w,h,x,y,ka(ell),l1,l2);
        end
    end

    %observation on the left of the defect 
    ind=0; 
    for i=1:length(x)
        if x(i)>supp
            ind=i-1;
            break
        end
    end
    xx=x(1:ind);
    mesu=temp1(length(y),1:ind);
    mesv=temp2(length(y),1:ind);
    %computation of coefficients c_1^1, c_2^1,c_2^1, c_2^2
    [ca,cb]=decompsurface(ks,ka,h,mesu,mesv,xx,w);
    ca=ca/h/ss(k1,h);
    cb=cb/h/ss(k1,h);
    %storage of data 
    Xi1(j)=2*k1; 
    Xi2(j)=k1+ka(1); 
    F1(j)=(ca(1)-cb(1))/2; 
    F2(j)=(ca(length(ks)+1)-cb(length(ks)+1))/2; 
end

%comparison with the true fourier tranform 
xi=linspace(0,50,300);
fourier1=xi*0; 
fourier2=xi*0; 
for i=1:length(xi)
    fourier1(i)=integral(@(x) (g1prim(x)-g2prim(x)).*exp(1i*xi(i)*x),3,5);
    fourier2(i)=integral(@(x) (g1prim(x)+g2prim(x)).*exp(1i*xi(i)*x),3,5);
end

%%%%least square algorithm 
% area of search 
a=3; 
b=5; 
grilx=100; 
x=linspace(a,b,grilx); 

%filtre and least square 
[t1,t2]=filtrefourier(Xi1,F1);
zz=solve1Dcarreb(t2,t1,a,b,grilx,0.08); %moindre carrés gauche
[t1,t2]=filtrefourier(Xi2,F2);
tt=solve1Dcarreb(t2,t1,a,b,grilx,0.08);

%approximations of g1prim and g2prim 
gg1prim=(zz+tt)/2;
gg2prim=(tt-zz)/2;
gg1prim=gg1prim-mean(gg1prim);
gg2prim=gg2prim-mean(gg2prim);

%approimations of g1 and g2
gg1=0*gg1prim;
gg2=0*gg2prim;

for i=2:length(zz)
    gg1(i)=gg1(i-1)+gg1prim(i)*(x(i)-x(i-1));
    gg2(i)=gg2(i-1)+gg2prim(i)*(x(i)-x(i-1));
end

figure 
plot(x,g1(x))
hold on 
plot(x,gg1)

figure
plot(x,g2(x))
hold on 
plot(x,gg2)

%relative error 
error=max(norm(g1(x)-gg1)/norm(g1(x)),norm(g2(x)-gg2)/norm(g2(x))); 

