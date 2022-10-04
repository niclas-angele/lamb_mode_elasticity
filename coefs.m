function [a,b,F1,F2]=coefs(x,jn,k,u,v,f1,f2,b1top,b2top,b1bot,b2bot,h)
%%% return coefficients a,b and modal source terms F1,F2 in a lamb
%%% decomposition computed on x. jn is the scalar product of Xn and Yn, k
%%% the wavenumber, u,v the fields, h the widht, and all the source terms
    z=linspace(-10,10,10001);
    F1=1/jn*conly(f1,u,z,h) +(1/jn*b1top(z,h)*u(h)+1/jn*b1bot(z,-h)*u(-h)).';
    F2=1/jn*conly(f2,v,z,h) +(1/jn*b2top(z,h)*v(h)+1/jn*b2bot(z,-h)*v(-h)).';
    a=1/2*(convx1(z,F1,k,x)-convx2(z,F2,k,x));
    b=1/2*(convx2(z,F1,k,x)-convx1(z,F2,k,x));
end
    
%fast convolutions functions
function N=convx1(z,F1,k,x)
    M=exp(1i*k*abs(ones(length(x),1)*z-x.'*ones(1,length(z))));
    pas=z(2)-z(1);
    N=pas*M*F1;
end

function N=convx2(z,F2,k,x)
    f=@(z,x) sign(x-z).*exp(1i*k*abs(x-z));
    M=f(ones(length(x),1)*z,x.'*ones(1,length(z)));
    pas=z(2)-z(1);
    N=pas*M*F2;
end

function N=conly(f1,u,x,h) 
    y=linspace(-h,h,1001);
    pas=y(2)-y(1);
    M=f1(x.'*ones(1,length(y)),ones(length(x),1)*y);
    N=pas*M*(u(y).');
end