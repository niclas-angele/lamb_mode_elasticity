function [UX,UY,UZ]=solveLamb3D(w,fz,fl,gradfz,divfl,blt,bzt,gradbzt,divblt,blb,bzb,gradbzb,divblb, fsh,bsht,bshb, x,y,h)
%%% elastic displacement in a plate at frequency w and width h in a grid
%%% [x,y]. Source terms f and b are given decomposed using the HHD and we
%%% also ask for their gradient. 
    ct=0.5;
    cl=0.9;
    mu=ct^2; 
    lambda=cl^2-2*ct^2;
    p=@(k) sqrt(w.^2./(cl^2)-k.^2);
    q=@(k) sqrt(w.^2./(ct^2)-k.^2);
    jns=@(k) 1i*mu*k.*(q(k).^2+k.^2).*(h*(q(k).^2-k.^2).^2.*sin(q(k)*h).^2+h*4*k.^2.*p(k).^2.*sin(p(k)*h).^2+sin(p(k)*h).*sin(q(k)*h).*(1./p(k).*(q(k).^2-k.^2).*(q(k).^2-k.^2-8*p(k).^2).*sin(q(k)*h).*cos(p(k)*h)+4./q(k).*p(k).^2.*(2*q(k).^2-k.^2).*sin(p(k)*h).*cos(q(k)*h)));
    us=@(k,y) (1i*k*(q(k)^2-k^2)*sin(q(k)*h)*cos(p(k)*y)-2*1i*k*p(k)*q(k)*sin(p(k)*h)*cos(q(k)*y)); 
    vs=@(k,y) (-p(k)*(q(k)^2-k^2)*sin(q(k)*h)*sin(p(k)*y)-2*k^2*p(k)*sin(p(k)*h)*sin(q(k)*y));
    ss=@(k,y) -(q(k)^2-k^2)*(cl^2*k^2+lambda*p(k)^2)*sin(q(k)*h)*cos(p(k)*y)+4*mu*p(k)*q(k)*k^2*sin(p(k)*h)*cos(q(k)*y);
    ts=@(k,y) 2*1i*k*mu*(q(k)^2-k^2)*p(k)*(-sin(q(k)*h)*sin(p(k)*y)+sin(p(k)*h)*sin(q(k)*y));
    ua=@(k,y) 1i*k*(q(k)^2-k^2)*cos(q(k)*h)*sin(p(k)*y)-2*1i*k*p(k)*q(k)*cos(p(k)*h)*sin(q(k)*y);
    va=@(k,y) p(k)*(q(k)^2-k^2)*cos(q(k)*h)*cos(p(k)*y)+2*k^2*p(k)*cos(p(k)*h)*cos(q(k)*y);
    sa=@(k,y) -(q(k)^2-k^2)*(cl^2*k^2+lambda*p(k)^2)*cos(q(k)*h)*sin(p(k)*y)+4*mu*p(k)*q(k)*k^2*cos(p(k)*h)*sin(q(k)*y);
    ta=@(k,y) 2*1i*k*mu*(q(k)^2-k^2)*p(k)*(cos(q(k)*h)*cos(p(k)*y)-cos(p(k)*h)*cos(q(k)*y));
    jna=@(k) 1i*mu*k.*(q(k).^2+k.^2).*(h*(q(k).^2-k.^2).^2.*cos(q(k)*h).^2+h*4*k.^2.*p(k).^2.*cos(p(k)*h).^2-cos(p(k)*h).*cos(q(k)*h).*(1./p(k).*(q(k).^2-k.^2).*(q(k).^2-k.^2-8*p(k).^2).*cos(q(k)*h).*sin(p(k)*h)+4./q(k).*p(k).^2.*(2*q(k).^2-k.^2).*cos(p(k)*h).*sin(q(k)*h)));
    [S,A]=dispersionh(w,5*w,h);
    S=S(1:min(length(S),10)); 
    A=A(1:min(length(A),10)); 

%%% for a surface view 
    [X,Y]=meshgrid(x,y); 
    UX=X*0; 
    UY=X*0; 
    UZ=X*0;
    %we compute coefficients and look at the surface 
    for i=1:length(S)
        k=S(i);
        [alp1,alp2,c]=coefs3DLamb(X,Y,fz,fl,gradfz,divfl,blt,bzt,gradbzt,divblt,blb,bzb,gradbzb,divblb,jns(k),@(y) us(k,y),@(y) vs(k,y),k,h);
        UX=UX+alp1*us(k,h); 
        UY=UY+alp2*us(k,h); 
        UZ=UZ+c*vs(k,h); 
    end
    for i=1:length(A)
        k=A(i);
        [alp1,alp2,c]=coefs3DLamb(X,Y,fz,fl,gradfz,divfl,blt,bzt,gradbzt,divblt,blb,bzb,gradbzb,divblb,jna(k),@(y) ua(k,y),@(y) va(k,y),k,h);
        UX=UX+alp1*ua(k,h); 
        UY=UY+alp2*ua(k,h); 
        UZ=UZ+c*va(k,h); 
    end

    SH=sqrt(w^2/mu-(0:10).^2*pi^2/4/h^2);
    for i=0:(length(SH)-1)
        if i==0
            phi=@(y) 1/sqrt(2*h);
        else
            phi=@(y) 1/sqrt(h)*cos(i*pi*(y+h)/2/h);
        end
        k=SH(i+1);
        [alp1,alp2]=coefs3DSH(X,Y,fsh,bsht,bshb,phi,k,h,mu);
        UX=UX+alp1*phi(h); 
        UY=UY+alp2*phi(h); 
    end

% %%% for a left side view 
%     z=linspace(-h,h);
%     [X,Y]=meshgrid(x,y); 
%     [t1,t2]=meshgrid(y,z); 
%     UX=t1*0; 
%     UY=t1*0; 
%     UZ=t1*0;
%     %we compute each mode on the side of the plate
%     for i=1:length(S)
%         k=S(i);
%         [alp1,alp2,c]=coefs3DLamb(X,Y,fz,fl,gradfz,divfl,blt,bzt,gradbzt,divblt,blb,bzb,gradbzb,divblb,jns(k),@(y) us(k,y),@(y) vs(k,y),k,h);
%         alp1=alp1(:,1).';
%         alp2=alp2(:,1).'; 
%         c=c(:,1).';
%         UX=UX+(us(k,z).')*alp1; 
%         UY=UY+(us(k,z).')*alp2; 
%         UZ=UZ+(vs(k,z).')*c; 
%     end
%     for i=1:length(A)
%         k=A(i);
%         [alp1,alp2,c]=coefs3DLamb(X,Y,fz,fl,gradfz,divfl,blt,bzt,gradbzt,divblt,blb,bzb,gradbzb,divblb,jna(k),@(y) ua(k,y),@(y) va(k,y),k,h);
%         alp1=alp1(:,1).';
%         alp2=alp2(:,1).'; 
%         c=c(:,1).';
%         UX=UX+(ua(k,z).')*alp1; 
%         UY=UY+(ua(k,z).')*alp2; 
%         UZ=UZ+(va(k,z).')*c; 
%     end
% 
%     SH=sqrt(w^2/mu-(0:10).^2*pi^2/4/h^2);
%     for i=0:(length(SH)-1)
%         if i==0
%             phi=@(y) 1/sqrt(2*h);
%         else
%             phi=@(y) 1/sqrt(h)*cos(i*pi*(y+h)/2/h);
%         end
%         k=SH(i+1);
%         [alp1,alp2]=coefs3DSH(X,Y,fsh,bsht,bshb,phi,k,h,mu);
%         alp1=alp1(:,1).';
%         alp2=alp2(:,1).'; 
%         UX=UX+(phi(z).')*alp1; 
%         UY=UY+(phi(z).')*alp2; 
%     end


% %%% for a right side view 
%     z=linspace(-h,h);
%     [X,Y]=meshgrid(x,y); 
%     [t1,t2]=meshgrid(y,z); 
%     UX=t1*0; 
%     UY=t1*0; 
%     UZ=t1*0;
%     %we compute each coefficient 
%     for i=1:length(S)
%         k=S(i);
%         [alp1,alp2,c]=coefs3DLamb(X,Y,fz,fl,gradfz,divfl,blt,bzt,gradbzt,divblt,blb,bzb,gradbzb,divblb,jns(k),@(y) us(k,y),@(y) vs(k,y),k,h);
%         alp1=alp1(1,:);
%         alp2=alp2(1,:); 
%         c=c(1,:);
%         UX=UX+(us(k,z).')*alp1; 
%         UY=UY+(us(k,z).')*alp2; 
%         UZ=UZ+(vs(k,z).')*c; 
%     end
%     for i=1:length(A)
%         k=A(i);
%         [alp1,alp2,c]=coefs3DLamb(X,Y,fz,fl,gradfz,divfl,blt,bzt,gradbzt,divblt,blb,bzb,gradbzb,divblb,jna(k),@(y) ua(k,y),@(y) va(k,y),k,h);
%         alp1=alp1(1,:);
%         alp2=alp2(1,:); 
%         c=c(1,:);
%         UX=UX+(ua(k,z).')*alp1; 
%         UY=UY+(ua(k,z).')*alp2; 
%         UZ=UZ+(va(k,z).')*c; 
%     end
% 
%     SH=sqrt(w^2/mu-(0:10).^2*pi^2/4/h^2);
%     for i=0:(length(SH)-1)
%         if i==0
%             phi=@(y) 1/sqrt(2*h);
%         else
%             phi=@(y) 1/sqrt(h)*cos(i*pi*(y+h)/2/h);
%         end
%         k=SH(i+1);
%         [alp1,alp2]=coefs3DSH(X,Y,fsh,bsht,bshb,phi,k,h,mu);
%         alp1=alp1(1,:);
%         alp2=alp2(1,:); 
%         UX=UX+(phi(z).')*alp1; 
%         UY=UY+(phi(z).')*alp2; 
%     end


end