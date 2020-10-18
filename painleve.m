% define the system
deq = @(x,y) [y(2); x*y(1)+2*y(1)^3];

% set the integration bounds and output times
x0 = 5;
xn = -8;
sspan = linspace(x0,xn,1000);

% compute the initial values
y0 = [airy(x0); airy(1,x0)];

% set the tolerance and integrate the system
opts = odeset('RelTol',1e-13,'AbsTol',1e-15);
[x,y] = ode45(deq,sspan,y0,opts);

% the function q(x) is the initial entry of y
q=y(:,1);

% set the initial values
dI0=0;
I0=0;
J0=0;


% numerically integrate 
dI=-[0;cumsum((q(1:end-1).^2+q(2:end).^2)/2.*diff(x))]+dI0;
I=-[0;cumsum((dI(1:end-1)+dI(2:end))/2.*diff(x))]+I0;
J=-[0;cumsum((q(1:end-1)+q(2:end))/2.*diff(x))]+J0;

% here we have the distributions
F2=exp(-I);
F1=sqrt(F2.*exp(-J));
F4=sqrt(F2).*(exp(J/2)+exp(-J/2))/2;
x4=x/2^(2/3);

% and here are the probability distributions
f2=gradient(F2,x);
f1=gradient(F1,x);
f4=gradient(F4,x4);

% plotting the probability distributions
plot(x,f1,x,f2,x4,f4,'LineWidth',1)
ylabel('f_\beta(x)',"Rotation",0)
xlabel('x')
legend('\beta = 1','\beta = 2','\beta = 4')
