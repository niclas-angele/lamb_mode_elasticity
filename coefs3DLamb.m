function [alp1,alp2,c]=coefs3DLamb(X,Y,fz,fl,gradfz,divfl,blt,bzt,gradbzt,divblt,blb,bzb,gradbzb,divblb,jn,u,v,k,h)
%%% computation of coefficients in Lamb modes in 3D. alp1 coef u_x, alp2
%%% coef u_y, c coef u_z. X,Y are meshgrid of x,y where we compute
%%% wavefields, and we give all the needed sources. k is the wavenumber we
%%% look at, h the width and mu the 2 Lamé coefficient. jn, u and v are the
%%% scalar product of Xn,Yn and u_n, v_n associated to the mode we are
%%% looking at. 
    x=X(1,:); 
    y=Y(:,1)'; 
    % we extend the grid to a symmetric one 
    [x,y,i1,i2,i3,i4]=nouvxy(x,y); 
    [X,Y]=meshgrid(x,y); 
    xgap=x(length(x))-x(1); 
    ygap=y(length(y))-y(1); 
    xt=[x-xgap,x(2:length(x)-1),x+xgap]; 
    yt=[y-ygap,y(2:length(y)-1),y+ygap]; 
    [Xt,Yt]=meshgrid(xt,yt); 
    
    alp1=X*0; 
    alp2=X*0;
    c=X*0;
    alp=[alp1;alp2];
    %we compute F_n using rectangular integration 
    z=linspace(-h,h,200); 
    for i=1:length(z)-1
        alp=alp-gradfz(X,Y,X*0+z(i))*v(z(i))+1i*k*fl(X,Y,X*0+z(i))*u(z(i));
        c=c-1i*k*fz(X,Y,X*0+z(i))*v(z(i))+divfl(X,Y,X*0+z(i))*u(z(i));
    end 
    alp=(z(2)-z(1))*alp; 
    c=(z(2)-z(1))*c; 

    alp=alp-gradbzt(X,Y)*v(h)-gradbzb(X,Y)*v(-h)+1i*k*blt(X,Y)*u(h)+1i*k*blb(X,Y)*u(-h);
    c=c-1i*k*bzt(X,Y)*v(h)-1i*k*bzb(X,Y)*v(-h)+divblt(X,Y)*u(h)+divblb(X,Y)*u(-h);
    ind=length(alp(:,1));
    alp1=alp(1:ind/2,:); 
    alp2=alp(ind/2+1:ind,:); 
    %we convol with Bessel functions 
    B=-besselh(0,k*sqrt((Xt.^2+Yt.^2))); 
    alp1=1i/4/jn*(x(2)-x(1))*(y(2)-y(1))*conv2(alp1,B,'same'); 
    alp2=1i/4/jn*(x(2)-x(1))*(y(2)-y(1))*conv2(alp2,B,'same'); 
    c=1i/4/jn*(x(2)-x(1))*(y(2)-y(1))*conv2(c,B,'same'); 
    %we come back to the initial grid 
    alp1=alp1(i3:i4,i1:i2); 
    alp2=alp2(i3:i4,i1:i2);
    c=c(i3:i4,i1:i2); 
end


function [p,q,x1,x2,y1,y2]=nouvxy(x,y)
%%% extend lists to symmetric lists
    pasx=x(2)-x(1); 
    pasy=y(2)-y(1); 
    if abs(x(length(x)))>abs(x(1))
        p=-x(length(x)):pasx:x(1)-pasx;
        x1=length(p)+1; 
        x2=length(p)+length(x); 
        p=[p,x]; 
    elseif abs(x(length(x)))<abs(x(1))
        p=x(length(x)):pasx:-x(1)-pasx; 
        x1=1; 
        x2=length(x); 
        p=[x,p]; 
    else
        p=x;
        x1=1; 
        x2=length(x); 
    end 
    if abs(y(length(y)))>abs(y(1))
        q=-y(length(y)):pasy:y(1)-pasy;
        y1=length(q)+1; 
        y2=length(q)+length(y);
        q=[q,y]; 
    elseif abs(y(length(y)))<abs(y(1))
        q=y(length(y)):pasy:-y(1)-pasy; 
        q=[y,q]; 
        y1=1; 
        y2=length(y); 
    else
        q=y;
        y1=1; 
        y2=length(y); 
    end 
end
