function [U,V]=correction_PML_droit(u,mesh,w,h,x,y,k,l1,l2)
%%% correction of a PML using the method developped in [Bonnet Ben-Dhia
%%% 2014]. u is a waveguide generated on the mesh mesh at the frequency w.
%%% We consider a plate of width 2h and look in the area [x,y]. k is the
%%% backward wavenumber and l1, l2 the right and left coordinate of PML. We
%%% return U,V the corrected displacement fields
    ct=0.5;
    cl=0.9;
    [ar,br]=projection(u,mesh,w,h,l1,y,k);
    [al,bl]=projection(u,mesh,w,h,l2,y,k);
    aplusl=1/2*(al+bl);
    amoinsr=1/2*(-ar+br);
    [X,Y]=meshgrid(x,y); 
    a=(-aplusl*exp(1i*k*(X-l2))+amoinsr*exp(-1i*k*(X-l1))).*(X<l1).*(X>l2);
    b=(-aplusl*exp(1i*k*(X-l2))-amoinsr*exp(-1i*k*(X-l1))).*(X<l1).*(X>l2);
    temp1=P1togrid(mesh,u(:,1),x,y);
    temp2=P1togrid(mesh,u(:,2),x,y);
    p=@(k) sqrt(w.^2./(cl^2)-k.^2);
    q=@(k) sqrt(w.^2./(ct^2)-k.^2);
    us=@(k,y) (1i*k*(q(k)^2-k^2)*sin(q(k)*h)*cos(p(k)*y)-2*1i*k*p(k)*q(k)*sin(p(k)*h)*cos(q(k)*y)); 
    vs=@(k,y) (-p(k)*(q(k)^2-k^2)*sin(q(k)*h)*sin(p(k)*y)-2*k^2*p(k)*sin(p(k)*h)*sin(q(k)*y));
    U=temp1+a.*us(k,Y); 
    V=temp2+b.*vs(k,Y); 
end
