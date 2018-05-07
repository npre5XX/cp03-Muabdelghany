% Gauss Quadrature
clc;
clear;
N=4;         % order of the discrete ordinate
l=240;       % number of space meshes
L=6;         % thikness of the whole slab
del=L/(l/2);
mu=[0.33998 0.86114 -0.86114 -0.33998];   % Gauss values of the angle
w=[0.65215 0.34785 0.34785 0.65215];

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

% calculating the scalar flux
    for i=1:l+1
        flux(i)=0;
    for n=1:N
        flux(i)=flux(i)+eps(i,n,m);
    end
    end
    
% Calculating the rightward angular flux
for i=1:l+1
    Reps(i)=0;
    Leps(i)=0;
for n=1:2
    Reps(i)=Reps(i)+eps(i,n,m);
end

for n=3:4
    Leps(i)=Leps(i)+eps(i,n,m);
end
end
    
    
m    % number of iterations till we match our criteria
     % plotting eps as a function of position for each ordinate angle
x=linspace(0,L,l+1);

figure(1)
plot(x,Reps(:))
grid on
xlabel('x [Cm]')
ylabel('\psi_+(x)')
title('Rightward Angular Flux Distribution \psi_R(x) Comparison Using Gauss Vs Radau Quadrature')
hold on

figure(2)
plot(x,Leps(:))
grid on
xlabel('x [Cm]')
ylabel('\psi_L(x)')
title('Leftward Angular Flux Distribution \psi_-(x) Comparison Using Gauss Vs Radau Quadrature')
hold on

figure(3)
plot(x,flux(:))
grid on
xlabel('x [Cm]')
ylabel('\phi(x)')
title('Scalar Flux Distribution \phi(x) Comparison Using Gauss Vs Radau Quadrature')
hold on