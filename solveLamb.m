function [U,V,SS,T]=solveLamb(w,f1,f2,b1top,b2top,b1bot,b2bot,x,y,h)
%%% return the solutions (U,V) of the elastic equation in a waveguide of
%%% width h at frequency w at coordinate x,y. we also return the stress tensor (S,T).
%%% f1,f2,b1top,b2top,b1bot,2bot are the source terms.
    ct=0.5;
    cl=0.9;
    mu=ct^2; 
    lambda=cl^2-2*ct^2;
    p=@(k) sqrt(w.^2./(cl^2)-k.^2);
    q=@(k) sqrt(w.^2./(ct^2)-k.^2);
    %symmetric mode
    jns=@(k) 1i*mu*k.*(q(k).^2+k.^2).*(h*(q(k).^2-k.^2).^2.*sin(q(k)*h).^2+h*4*k.^2.*p(k).^2.*sin(p(k)*h).^2+sin(p(k)*h).*sin(q(k)*h).*(1./p(k).*(q(k).^2-k.^2).*(q(k).^2-k.^2-8*p(k).^2).*sin(q(k)*h).*cos(p(k)*h)+4./q(k).*p(k).^2.*(2*q(k).^2-k.^2).*sin(p(k)*h).*cos(q(k)*h)));
    us=@(k,y) (1i*k*(q(k)^2-k^2)*sin(q(k)*h)*cos(p(k)*y)-2*1i*k*p(k)*q(k)*sin(p(k)*h)*cos(q(k)*y)); 
    vs=@(k,y) (-p(k)*(q(k)^2-k^2)*sin(q(k)*h)*sin(p(k)*y)-2*k^2*p(k)*sin(p(k)*h)*sin(q(k)*y));
    ss=@(k,y) -(q(k)^2-k^2)*(cl^2*k^2+lambda*p(k)^2)*sin(q(k)*h)*cos(p(k)*y)+4*mu*p(k)*q(k)*k^2*sin(p(k)*h)*cos(q(k)*y);
    ts=@(k,y) 2*1i*k*mu*(q(k)^2-k^2)*p(k)*(-sin(q(k)*h)*sin(p(k)*y)+sin(p(k)*h)*sin(q(k)*y));
    %antisymmetric mode 
    ua=@(k,y) 1i*k*(q(k)^2-k^2)*cos(q(k)*h)*sin(p(k)*y)-2*1i*k*p(k)*q(k)*cos(p(k)*h)*sin(q(k)*y);
    va=@(k,y) p(k)*(q(k)^2-k^2)*cos(q(k)*h)*cos(p(k)*y)+2*k^2*p(k)*cos(p(k)*h)*cos(q(k)*y);
    sa=@(k,y) -(q(k)^2-k^2)*(cl^2*k^2+lambda*p(k)^2)*cos(q(k)*h)*sin(p(k)*y)+4*mu*p(k)*q(k)*k^2*cos(p(k)*h)*sin(q(k)*y);
    ta=@(k,y) 2*1i*k*mu*(q(k)^2-k^2)*p(k)*(cos(q(k)*h)*cos(p(k)*y)-cos(p(k)*h)*cos(q(k)*y));
    jna=@(k) 1i*mu*k.*(q(k).^2+k.^2).*(h*(q(k).^2-k.^2).^2.*cos(q(k)*h).^2+h*4*k.^2.*p(k).^2.*cos(p(k)*h).^2-cos(p(k)*h).*cos(q(k)*h).*(1./p(k).*(q(k).^2-k.^2).*(q(k).^2-k.^2-8*p(k).^2).*cos(q(k)*h).*sin(p(k)*h)+4./q(k).*p(k).^2.*(2*q(k).^2-k.^2).*cos(p(k)*h).*sin(q(k)*h)));
    %full dispersion relation 
    [S,A]=dispersionh(w,5*w,h);
    [X,Y]=meshgrid(x,y); 
    %we cut to 20 modes 
    S=S(1:min(length(S),20)); 
    A=A(1:min(length(A),20)); 
    U=X*0; 
    V=X*0; 
    SS=X*0;
    T=X*0;
     %we compute the coefficients on each lamb mode 
    for i=1:length(S)
        k=S(i);
        [a,b]=coefs(x,jns(k),k,@(y) us(k,y),@(y) vs(k,y),f1,f2,b1top,b2top,b1bot,b2bot,h);
        U=U+(ones(length(y),1)*(a.')).*us(k,Y);
        V=V+(ones(length(y),1)*(b.')).*vs(k,Y);
        SS=SS+(ones(length(y),1)*(b.')).*ss(k,Y);
        T=T+(ones(length(y),1)*(a.')).*ts(k,Y);
    end
    for i=1:length(A)
        k=A(i);
        [a,b]=coefs(x,jna(k),k,@(y) ua(k,y),@(y) va(k,y),f1,f2,b1top,b2top,b1bot,b2bot,h);
        U=U+(ones(length(y),1)*(a.')).*ua(k,Y);
        V=V+(ones(length(y),1)*(b.')).*va(k,Y);
        SS=SS+(ones(length(y),1)*(b.')).*sa(k,Y);
        T=T+(ones(length(y),1)*(a.')).*ta(k,Y);
    end
    U=U.*(Y<h).*(Y>-h);
    V=V.*(Y<h).*(Y>-h);
    SS=SS.*(Y<h).*(Y>-h);
    T=T.*(Y<h).*(Y>-h);
end
    
    