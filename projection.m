function [a,b]=projection(U,mesh,w,H,x,y,k)
%%% projection of an elastic wavefield U defined on a mesh mesh at
%%% frequency w on a grid x,y. H is the width of the waveguide, k is the
%%% wavenumber we want to project on. a=1/jn <U,Xn>, b=1/jn <V,Yn>. 
    ct=0.5;
    cl=0.9;
    mu=ct^2; 
    lambda=cl^2-2*ct^2;
    u=U(:,1); 
    v=U(:,2); 
    grad1=mesh.grad(u);
    grad2=mesh.grad(v); 
    dxu=P1togrid(mesh,P0toP1(mesh,grad1(:,1)),x,y);
    dyu=P1togrid(mesh,P0toP1(mesh,grad1(:,2)),x,y);
    dxv=P1togrid(mesh,P0toP1(mesh,grad2(:,1)),x,y);
    dyv=P1togrid(mesh,P0toP1(mesh,grad2(:,2)),x,y);
    s=lambda*dyv+(lambda+2*mu)*dxu;
    t=mu*(dyu+dxv);
    u=P1togrid(mesh,u,x,y); 
    v=P1togrid(mesh,v,x,y); 
    h=H; 
    a=x*0;
    b=x*0;
    for i=1:length(x)
        p=sqrt(w.^2./(cl^2)-k.^2);
        q=sqrt(w.^2./(ct^2)-k.^2);
        jns=1i*mu*k.*(q.^2+k.^2).*(h*(q.^2-k.^2).^2.*sin(q*h).^2+h*4*k.^2.*p.^2.*sin(p*h).^2+sin(p*h).*sin(q*h).*(1./p.*(q.^2-k.^2).*(q.^2-k.^2-8*p.^2).*sin(q*h).*cos(p*h)+4./q.*p.^2.*(2*q.^2-k.^2).*sin(p*h).*cos(q*h)));
        us= (1i*k*(q^2-k^2)*sin(q*h)*cos(p*y)-2*1i*k*p*q*sin(p*h)*cos(q*y)); 
        vs= (-p*(q^2-k^2)*sin(q*h)*sin(p*y)-2*k^2*p*sin(p*h)*sin(q*y));
        ss= -(q^2-k^2)*(cl^2*k^2+lambda*p^2)*sin(q*h)*cos(p*y)+4*mu*p*q*k^2*sin(p*h)*cos(q*y);
        ts= 2*1i*k*mu*(q^2-k^2)*p*(-sin(q*h)*sin(p*y)+sin(p*h)*sin(q*y));
        scal1=rmmissing(-u(:,i).*(ss.')+t(:,i).*(vs.')); 
        scal2=rmmissing(-(us.').*s(:,i)+(ts.').*v(:,i)); 
        a(i)=sum(scal1(2:length(scal1)))*(y(2)-y(1))/jns; 
        b(i)=sum(scal2(2:length(scal2)))*(y(2)-y(1))/jns; 
    end
end

