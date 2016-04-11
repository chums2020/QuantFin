load Factors_July26_July11.dat      % This file, downloaded from French's website, has 5 columns: 
                      %     date, mkt excess return, size factor, book/ME factor, T-bill rate
load Indu10_July26_July11.dat     % The monthly returns on the 10 industry portfolios           
T=1021;
N=10;
bill=Factors_July26_July11(1:T,5) / 100;        
Re=Indu10_July26_July11(1:T,2:11)  / 100;              
ReEX=Re - bill*ones(1,N);                             
gamma=3;               % risk aversion parameter
mu=mean(ReEX);         % estimate the mean and covariance matrix of the 5 excess return 
mu=mu';                % making it a column vector
V=cov(ReEX);
%1)  
a=zeros(N,1);      
b=0.3*ones(N,1);   
gam=0;
c=-mu;
G=gamma*V;
xu=a;
xo=b;
xx=ones(N,1)/N;        % search the optimal portfolio weights 
                        % starting from the equal-weighted 
prt=0;
[x,fct,ier,nsub]=minq(gam,c,G,xu,xo,prt,xx);            % calling the minimization program to get the minimized value and the associated x
w=x;          % an N-vector of the optimal weights subject to the constraints 
w0=(1/gamma)*V^(-1)*mu;          % the optimal unconstrained weights, equation (174)
 
%2)   What are the mean, std and Sharpe ratio of the constrained optimal portfolio?
          
Port=zeros(T,1);      % define this T vector to store the returns on the portfolio
PortS=zeros(T,1);     % to store the accumulative returns   
 
Port(1)=bill(1)+w'*ReEX(1,:)';                   % return in the first period
                                                 % notice that  ReEX(1,:)  is the 1st row, the traspose makes it a column
PortS(1)=1+Port(1);
 
for t=2:T,
   Port(t)=bill(t) + w'*ReEX(t,:)';                 % return at t
   PortS(t)=PortS(t-1)*(1+Port(t));                 % accumulativ return at t
end;
 
PortEX=Port -  bill;                     % the portfolio excess returns 
muP=mean(PortEX);
sig2P=var(PortEX);
SharpeP=muP / sqrt(sig2P);
 
%3)   What are the mean, std and Sharpe ratio of the unconstrained optimal portfolio?
           
Port0=zeros(T,1);      % define this T vector to store the returns on the portfolio
PortS0=zeros(T,1);    % to store the accumulative returns   
 
Port0(1)=bill(1)+w0'*ReEX(1,:)';     % return in the first period
                                     
PortS0(1)=1+Port0(1);
 
for t=2:T,
   Port0(t)=bill(t) + w0'*ReEX(t,:)';                    % return at t
   PortS0(t)=PortS0(t-1)*(1+Port0(t));                   % accumulativ return at t
end;
 
Port0EX=Port0 -  bill;                     % the portfolio excess returns
muP0=mean(Port0EX);
sig2P0=var(Port0EX);
SharpeP0=muP0 / sqrt(sig2P0);
 
%4)   What are the mean, std and Sharpe ratio of the mkt index portfolio?  
 
RmEX=Factors_July26_July11(1:T,2) / 100;          % Get the returns on mkt in excess of T-bill rate, the data are such a way
mu1=mean(RmEX);              % estimate the mean and variance of the mkt excess return   
sig2=var(RmEX);
SharpeRM=mu1/sqrt(sig2);     % Sharpe ratio of the mkt
 
%5)   What are the accumulative returns of investing $1 into the mkt index portfolio, the constrained optimal portfolio, and the unconstrained optimal portfolio respectively for the relevant time period?
 
Rm=RmEX + bill;                           % mkt return, obtained by adding back the riskfree rate
CC_RM=zeros(T,1);
CC_RM(1)=1+Rm(1);                         %  accumulativ return in the first month
for t=2:T,
   CC_RM(t)=CC_RM(t-1)*(1+Rm(t));            % accumulativ return at t
end;
 
fprintf('******  weights on the risky assets: constrained********** \n');
w
fprintf('******  Mean, std abd Sharpe ratios of the constrained portfolios ********** \n');
[muP, sqrt(sig2P),SharpeP]
fprintf('******  Mean, std abd Sharpe ratios of the un-c* portfolios ********** \n');
[muP0, sqrt(sig2P0), SharpeP0]
fprintf('******  Mean, std abd Sharpe ratios of the mkt portfolios ********** \n');
[mu1, sqrt(sig2), SharpeRM]
fprintf('******  Accu return of the mkt, constrained and un-c* the ports  ********** \n');
CC_RM(T) 
PortS(T) 
PortS0(T)

