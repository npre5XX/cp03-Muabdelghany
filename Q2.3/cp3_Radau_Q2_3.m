% Accuracy Convergence of the Sn with Radau Quadrature
clc;
clear;
close all;
N=0;         % order of the discrete ordinate
L=6;         % thickness of the whole slab
l=162;
del=L/(l/2);

j=1;

while N<8
    
    N=N+2
    
    if N==2
        mu=[1 -1]; 
        w=[1 1];
    elseif N==4
        mu=[0.4472136 1 -1 -0.4472136];
        w=[0.833333 0.16666667 0.16666667 0.833333];
    elseif N==6
        mu=[0.0813570 0.7650553 1 -1 -0.7650553 -0.0813570];
        w=[0.5548584 0.3784750 0.0666667 0.0666667 0.3784750 0.5548584];
    else
       mu=[0.2092992 0.5917002 0.8717401 1 -1 -0.8717401 -0.5917002 -0.2092992];   
       w=[0.4124591 0.3411227 0.2107042 0.0357143 0.0357143 0.2107042 0.3411227 0.4124591];  
    end
    

m=1;
max_ephi=5;
while max_ephi>0.000001         % converging criteria
    
m=m+1;

if m==2;
    
    for i=2:2:l
        
        if i<=l/3        % taking care of the different regions
          s=1;
        else
          s=0;
        end
     Q(i,1)=0.5*s;
     phi(i,m-1)=0;
    end
else
    
    for i=2:2:l
    phi(i,m-1)=0;
    for n=1:N
    phi(i,m-1)=phi(i,m-1)+(w(n)*eps(i,n,m-1));
    end
    ephi(i/2,m-2)=(abs(phi(i,m-1)-phi(i,m-2))/phi(i,m-1));
    
    if i<=l/3
       seg_s=0.5;s=1;
    else
       seg_s=0.3;s=0;
    end
    Q(i,m-1)=0.5*((seg_s*phi(i,m-1))+s);
    end
    max_ephi=max(ephi);
end

for n=1:N
    
if mu(n)>0      % forward sweeping
    
    eps(1,n,m)=0;
    
for i=2:2:l

    if i<=l/3
       seg_t=1;seg_s=0.5;s=1;
    else
       seg_t=1.5;seg_s=0.3;s=0;
    end
    
A(i,n)=((2*mu(n)-seg_t*del)/(2*mu(n)+seg_t*del));
B(i,n)=((2*del)/(2*mu(n)+seg_t*del));

eps(i+1,n,m)=(A(i,n)*eps(i-1,n,m))+(B(i,n)*Q(i,m-1));

end

else                 % backward sweeping
    eps(l+1,n,m)=0;
    
for i=l:-2:2

    if i<=l/3
       seg_t=1;seg_s=0.5;s=1;
    else
       seg_t=1.5;seg_s=0.3;s=0;
    end
A(i,n)=((2*mu(n)+seg_t*del)/(2*mu(n)-seg_t*del));
B(i,n)=((-2*del)/(2*mu(n)-seg_t*del));

eps(i-1,n,m)=(A(i,n)*eps(i+1,n,m))+(B(i,n)*Q(i,m-1));

end

end

end

    for n=1:N
    for i=2:2:l
        eps(i,n,m)=0.5*(eps(i+1,n,m)+eps(i-1,n,m));  % calculating central flux
    end
    end
    
end
m    % number of iterations till we match our criteria

    for i=1:l+1
        flux(i,j)=0;
    for n=1:N
        flux(i,j)=flux(i,j)+eps(i,n,m);
    end
    end
    
    if j>1
    for i=1:l+1
        eflux=(abs(flux(i,j)-flux(i,j-1))/flux(i,j));
    end
        j=j+1
        max_eflux(j-1)=max(eflux);
    else
        j=j+1
        max_eflux(1)=1;
    end 
    
        end

% ploting the error vs number of space meshes
    N=[2 4 6 8];
    plot(N,max_eflux(:),'-- s')
    xlim([0 10])
    grid on
    xlabel('# Radau Quadratures')
    ylabel('Error (\epsilon)')
    title('Error in the Scalar Flux Vs Number of Quadratures')
