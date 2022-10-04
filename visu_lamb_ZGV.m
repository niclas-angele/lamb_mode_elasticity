clear; close all;

%parameter of the plate
ct=0.5;
cl=0.9;
mu=ct^2; 
lambda=cl^2-2*ct^2;
w=13.7;
h=@(x) 0.1+x*0; 

%source terms 
xs=0.1;
sig=0.1;
delta=@(x) 1/sig/sqrt(2*pi)*exp(-(x-xs).^2/2/sig/sig);
f1=@(x,y) 100*(x+2*y).*(1.*(((x-0.5)).^2+((y-0.06)/0.015).^2)<1).*(1-(((x-0.5)).^2+((y-0.06)/0.015).^2)); 
f2=@(x,y) 100*(1.*(((x-0.5)).^2+((y-0.06)/0.015).^2)<1).*(1-(((x-0.5)).^2+((y-0.06)/0.015).^2)); 
b1top=@(x,y) delta(x); 
b2top=@(x,y) x.*delta(x); 
b1bot=@(x,y) 20*(x-2).*(x-2.5).*(x>2).*(x<2.5); 
b2bot=@(x,y) 20*sin(x).*(x-2).*(x-2.5).*(x>2).*(x<2.5); 

%area of interest 
x=linspace(-6,6,500);
y=linspace(-h(0),h(0),200);
[X,Y]=meshgrid(x,y);

% % %%% resolution Lamb 

[U,V,SS,T]=solveLamb(w,f1,f2,b1top,b2top,b1bot,b2bot,x,y,h(0));
figure
subplot 221
title('real(u)')
surf(X,Y,real(U),'edgecolor','none')
view(0,90)
colorbar;
subplot 222
title('imag(u)')
surf(X,Y,imag(U),'edgecolor','none')
view(0,90)
colorbar;
subplot 223
title('real(v)')
surf(X,Y,real(V),'edgecolor','none')
view(0,90)
colorbar;
subplot 224
title('imag(v)')
surf(X,Y,imag(V),'edgecolor','none')
view(0,90)
colorbar;


% %%% résolution FEM
nodes = [-8 -0.1;8 -0.1;8 0.1;-8 0.1;];
edges = {{1,2},{2,3},{3,4},{4,1}};
dx = 0.1;
domain = Domain(nodes,edges);
mesh = Mesh(domain,dx);
mesh=mesh.submesh;
mesh=mesh.submesh;
mesh=mesh.submesh;

%cut here to get a bigger mesh 
mesh=mesh.submesh;
mesh=mesh.submesh;
mesh=mesh.submesh;
mesh=mesh.submesh;

%PML
C = mesh.P0({mu,lambda},4); 
bc1 = {1,1,0,{b1bot,b2bot}};
bc2 = {3,1,0,{b1top,b2top}};
l=4; 
at=10; 
alpha=@(x,y) -1i*w./(-1i*w+at*(x-l)).*(x>=l)-1i*w./(-1i*w+at*(-x-l)).*(x<=-l)+1*(x<l).*(x>-l);
C=mesh.P0({@(x,y) (lambda+2*mu)*alpha(x,y), 0, 0 , mu, 0 , lambda, @(x,y) mu*alpha(x,y) , 0 , 0 , mu, @(x,y) lambda*alpha(x,y), 0 , @(x,y) mu*alpha(x,y), 0, 0, lambda+2*mu}); 
%solution
u = mesh.solve_vect(C,mesh.P0(@(x,y) -w.^2),{f1,f2},{bc1,bc2});
temp1=P1togrid(mesh,u(:,1),x,y);
temp2=P1togrid(mesh,u(:,2),x,y);
%correction if backward modes 
[ks,ka]=dispersionreal(w,5*w,h(0));
    for ell=1:length(ks)
        if ks(ell)<0
            [temp1,temp2]=correction_PML_droit(u,mesh,w,h(0),x,y,ks(ell),-l,l);
        end
    end
    for ell=1:length(ka)
        if ka(ell)<0
            [temp1,temp2]=correction_PML_droit(u,mesh,w,h(0),x,y,ka(ell),-l,l);
        end
    end

figure 
subplot 221
surf(X,Y,real(temp1),'edgecolor','none');
view(0,90)
colorbar;

subplot 222
surf(X,Y,imag(temp1),'edgecolor','none');
view(0,90)
colorbar;

subplot 223
surf(X,Y,real(temp2),'edgecolor','none');
view(0,90)
colorbar;

subplot 224
surf(X,Y,imag(temp2),'edgecolor','none');
view(0,90)
colorbar;




