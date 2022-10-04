function [ca,cb]=decompsurface(ks,ka,h,mesu,mesv,x,w)
%%% function to compute the coefficients c_1^1, c_1^2, c_2^2, c_2^1. ks and
%%% ka are the wavenumbers, h the width, mesu and mesv the measurements at
%%% the surface of u and v, x the interval of measurement, w the frequency.
%%% We return ca and cb the list of all coefficients c_i^j
    ind=length(ks);
    tot=ind+length(ka);
    %first guess by solving a tot x tot linear system 
    [ca,cb]=guest(ks,ka,x,mesu,mesv,h,w);
    c=[real(ca),real(cb),imag(ca),imag(cb)];
    %least square to get a better approximation of each c_i^j
    gg=@(c) g(c,ks,ka,h,mesu,mesv,x,w);
    opts = optimset('Display','off'); 
    c=fminsearch(gg,c,opts);
    car=c(1:tot); 
    cbr=c((tot+1):2*tot);
    cai=c((2*tot+1):3*tot);
    cbi=c((3*tot+1):4*tot); 
    ca=car+1i*cai; 
    cb=cbr+1i*cbi; 
    ct=0.5;
    cl=0.9;
    mu=ct^2; 
    lambda=cl^2-2*ct^2;
    p=@(k) sqrt(w.^2./(cl^2)-k.^2);
    q=@(k) sqrt(w.^2./(ct^2)-k.^2);
    jns=@(k) 1i*mu*k.*(q(k).^2+k.^2).*(h*(q(k).^2-k.^2).^2.*sin(q(k)*h).^2+h*4*k.^2.*p(k).^2.*sin(p(k)*h).^2+sin(p(k)*h).*sin(q(k)*h).*(1./p(k).*(q(k).^2-k.^2).*(q(k).^2-k.^2-8*p(k).^2).*sin(q(k)*h).*cos(p(k)*h)+4./q(k).*p(k).^2.*(2*q(k).^2-k.^2).*sin(p(k)*h).*cos(q(k)*h)));
    jna=@(k) 1i*mu*k.*(q(k).^2+k.^2).*(h*(q(k).^2-k.^2).^2.*cos(q(k)*h).^2+h*4*k.^2.*p(k).^2.*cos(p(k)*h).^2-cos(p(k)*h).*cos(q(k)*h).*(1./p(k).*(q(k).^2-k.^2).*(q(k).^2-k.^2-8*p(k).^2).*cos(q(k)*h).*sin(p(k)*h)+4./q(k).*p(k).^2.*(2*q(k).^2-k.^2).*cos(p(k)*h).*sin(q(k)*h)));
    us=@(k,y) (1i*k*(q(k)^2-k^2)*sin(q(k)*h)*cos(p(k)*y)-2*1i*k*p(k)*q(k)*sin(p(k)*h)*cos(q(k)*y));
    ua=@(k,y) 1i*k*(q(k)^2-k^2)*cos(q(k)*h)*sin(p(k)*y)-2*1i*k*p(k)*q(k)*cos(p(k)*h)*sin(q(k)*y);
    for i=1:length(ks)
        k=ks(i);
        ca(i)=ca(i)*jns(k)/us(k,h);
        cb(i)=cb(i)*jns(k)/us(k,h); 
    end
    for i=1:length(ka)
        k=ka(i);
        j=i+length(ks);
        ca(j)=ca(j)*jna(k)/ua(k,h);
        cb(j)=cb(j)*jna(k)/ua(k,h); 
    end

end

function [ca,cb]=guest(ks,ka,x,mesu,mesv,h,w)
%%% fonction to have a first guess of the value of coefficient c_i^j using
%%% ks,ka the wavenumbers, mesu,mesv the measurements of u,v at widht h on the interval x and
%%% at frequency w. 
    tot=length(ks)+length(ka);
    %we cut the measurement interval in equal parts 
    in=floor(length(x)/2/tot);
    ind=in*(2*(1:tot)-1);
    %measurements points:
    y=x(ind);
    %Mu and Mv countain the measurements of Lamb modes in y 
    Mu=zeros(tot);
    Mv=zeros(tot);
    ct=0.5;
    cl=0.9;
    p=@(k) sqrt(w.^2./(cl^2)-k.^2);
    q=@(k) sqrt(w.^2./(ct^2)-k.^2);
    us=@(k,y) (1i*k*(q(k)^2-k^2)*sin(q(k)*h)*cos(p(k)*y)-2*1i*k*p(k)*q(k)*sin(p(k)*h)*cos(q(k)*y)); 
    vs=@(k,y) (-p(k)*(q(k)^2-k^2)*sin(q(k)*h)*sin(p(k)*y)-2*k^2*p(k)*sin(p(k)*h)*sin(q(k)*y));
    ua=@(k,y) 1i*k*(q(k)^2-k^2)*cos(q(k)*h)*sin(p(k)*y)-2*1i*k*p(k)*q(k)*cos(p(k)*h)*sin(q(k)*y);
    va=@(k,y) p(k)*(q(k)^2-k^2)*cos(q(k)*h)*cos(p(k)*y)+2*k^2*p(k)*cos(p(k)*h)*cos(q(k)*y);
    for i=1:length(ks)
        k=ks(i);
        Mu(:,i)=(us(k,h)*exp(-1i*k*y)).'; 
        Mv(:,i)=(vs(k,h)*exp(-1i*k*y)).'; 
    end
    for i=1:length(ka)
        k=ka(i);
        j=length(ks)+i;
        Mu(:,j)=(ua(k,h)*exp(-1i*k*y)).'; 
        Mv(:,j)=(va(k,h)*exp(-1i*k*y)).'; 
    end
    %we solve the systems to find each coordinate on each lamb mode 
    ca=Mu\(mesu(ind)).'; 
    cb=Mv\(mesv(ind)).'; 
end


function r=g(c,ks,ka,h,mesu,mesv,x,w)
    ind=length(ks);
    tot=ind+length(ka);
    car=c(1:tot); 
    cbr=c((tot+1):2*tot);
    cai=c((2*tot+1):3*tot);
    cbi=c((3*tot+1):4*tot); 
    r=aminimiser(ks,ka,x,car,cbr, cai,cbi,mesu,mesv,h,w);
end

function r=aminimiser(ks,ka,x,car,cbr, cai,cbi,mesu,mesv,h,w)
    t=0*x; 
    s=0*x; 
    ct=0.5;
    cl=0.9;
    p=@(k) sqrt(w.^2./(cl^2)-k.^2);
    q=@(k) sqrt(w.^2./(ct^2)-k.^2);
    us=@(k,y) (1i*k*(q(k)^2-k^2)*sin(q(k)*h)*cos(p(k)*y)-2*1i*k*p(k)*q(k)*sin(p(k)*h)*cos(q(k)*y)); 
    vs=@(k,y) (-p(k)*(q(k)^2-k^2)*sin(q(k)*h)*sin(p(k)*y)-2*k^2*p(k)*sin(p(k)*h)*sin(q(k)*y));
    ua=@(k,y) 1i*k*(q(k)^2-k^2)*cos(q(k)*h)*sin(p(k)*y)-2*1i*k*p(k)*q(k)*cos(p(k)*h)*sin(q(k)*y);
    va=@(k,y) p(k)*(q(k)^2-k^2)*cos(q(k)*h)*cos(p(k)*y)+2*k^2*p(k)*cos(p(k)*h)*cos(q(k)*y);
    for i=1:length(ks)
        k=ks(i);
        t=t+(car(i)+1i*cai(i))*us(k,h)*exp(-1i*k*x); 
        s=s+(cbr(i)+1i*cbi(i))*vs(k,h)*exp(-1i*k*x); 
    end
    for i=1:length(ka)
        k=ka(i);
        j=length(ks)+i;
        t=t+(car(j)+1i*cai(j))*ua(k,h)*exp(-1i*k*x); 
        s=s+(cbr(j)+1i*cbi(j))*va(k,h)*exp(-1i*k*x); 
    end
    r=norm(t-mesu)+norm(s-mesv); 
end