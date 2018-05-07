% Spatial Convergence of the Sn with Radau Quadrature
clc;
clear;
close all;
N=4;         % order of the discrete ordinate
L=6;         % thickness of the whole slab
mu=[0.4472136 1 -1 -0.4472136];   % Radau values of the angle
w=[0.833333 0.16666667 0.16666667 0.833333];
l=6;
j=1
max_eflux(1)=0.1;
while max_eflux(j)>0.0002
    
    l=3*l
    o(j)=l/2;
    del=L/(l/2);

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
    for i=2:l/3
        k=(i-1)*3;
        eflux=(abs(flux(k,j)+flux(k+1,j)+flux(k+2,j)-flux(i,j-1))/(2*flux(i,j-1)))-1;
    end
        j=j+1
        max_eflux(j)=max(eflux);
    else
        j=j+1
        max_eflux(2)=0.09;
    end 
    
end
o(j)=o(j-1)+3;

% ploting the error vs number of space meshes
    plot(o(:),max_eflux(:),'-- s')
    xlim([0 o(j)])
    grid on
    xlabel('# of space mesh')
    ylabel('Error (\epsilon)')
    title('Error in the Scalar Flux Vs Number of Space Meshes')
