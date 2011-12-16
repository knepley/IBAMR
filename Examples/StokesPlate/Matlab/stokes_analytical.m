% Analytical Expression 
Nx = 200;
dx = 1.5/Nx;
mu = 1/500;
t  = 5.56189;

for i = 0:Nx-1
    y(i+1) = erfc( (i*dx)/(2*sqrt(mu*t)) );
end


Y = [ fliplr(y) y];

for i = 0:2*Nx-1
    X(i+1) = i*dx;
end

plot(Y,X);
hold on;

load U_x_556189.curve
dis  = U_x_556189(:,1);
simU = U_x_556189(:,2);
plot(simU(1:1:end),dis(1:1:end),'rd')

