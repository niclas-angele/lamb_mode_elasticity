function [xi,f]=filtrefourier(xi,f)
%%% filter to delete absurd frequencies in the fourier transform. xi is the
%%% frequency data, f the fourier data associated 
     i=3; 
     while i<(length(f))
         %we remove data with too strong variations 
         if abs(abs(f(i))-abs(f(i+1)))>15*var(f,i)
             f(i+1)=[];
             xi(i+1)=[];
         else 
             i=i+1;
         end
     end
end


function v=var(D,n)
    v=abs(abs(D(n-2))-abs(D(n-1)))+abs(abs(D(n-1))-abs(D(n))); 
    v=abs(v)/2; 
end
