load Indu10_July26_July11.dat;
load Factors_July26_July11.dat;  
%1)   With the same risk aversion parameter 3, what are the optimal portfolio weights on the 10 industry assets?
gamma = 3;
T = length(Indu10_July26_July11(:,1));
N = 10;
bill=Factors_July26_July11(1:T,5) / 100;      % The T-bill rate or T-bill return
Re=Indu10_July26_July11(1:T,2:11) / 100;    % the 10  returns 
ReEX = Re - bill*ones(1,N);                      % the 10 returns in excess of the T-bill rates
mu=mean(ReEX);         % estimate the mean of the 10 excess return 
mu=mu';                % making it a column vector
V=cov(ReEX);           % estimate the covariance matrix of the 10 excess return 
w=(1/gamma)*V^(-1)*mu;  
fprintf('******  weights on the risky assets ********** \n');
w
 
%2)   What are the mean, std and Sharpe ratio of the optimal portfolio?
Port=zeros(T,1);        % define this 941 vector to store the returns on the portfolio
CC_Port=zeros(T,1);     % to store the accumulative returns   
 
Port(1)=bill(1)+w'*ReEX(1,:)';        % return in the first period
                                      % notice that ReEX(1,:)) is the 1st row, the traspose makes it a column
CC_Port(1)=1+Port(1);
 
for t=2:T,
   Port(t)=bill(t) + w'*ReEX(t,:)';       % return at t
   CC_Port(t)=CC_Port(t-1)*(1+Port(t));            % accumulativ return at t
end;
 
PortEX = Port -  bill;                     % the portfolio excess returns 
muP=mean(PortEX);
sig2P=var(PortEX);
SharpeP=muP / sqrt(sig2P);  
fprintf('******  Mean, std abd Sharpe ratios of the optimal portfolio ********** \n');
[muP, sig2P, SharpeP]
 
%3)   What are the mean, std and Sharpe ratio of the mkt index portfolio?  
RmEX=Factors_July26_July11(1:T,2) / 100;        % Get the returns on mkt in excess of T-bill rate, the data are such a way
mu1=mean(RmEX);              % estimate the mean and variance of the mkt excess return   
sig2=var(RmEX);
SharpeRM=mu1/sqrt(sig2);     % Sharpe ratio of the mkt
fprintf('******  Mean, std abd Sharpe ratios of the mkt ********** \n');
[mu1, sig2, SharpeRM]
 
%4)   What are the accumulative returns of investing $1 into the mkt index portfolio and the optimal portfolio respectively for the relevant time period?
Rm=RmEX + bill;                         % mkt return, obtained by adding back the riskfree rate
CC_RM=zeros(T,1);
CC_RM(1)=1+Rm(1);                         %  accumulativ return in the first month
for t=2:T,
   CC_RM(t)=CC_RM(t-1)*(1+Rm(t));            % accumulativ return at t
end;
 
fprintf('******  Accu return of the mkt and the port   ********** \n');
[CC_RM(T), CC_Port(T)]
 
%5)    Re-do 4) with risk aversion parameter=1 and 5, respectively. 
gamma_1 = 1;
w=(1/gamma_1)*V^(-1)*mu; 
Port=zeros(T,1);        % define this 941 vector to store the returns on the portfolio
CC_Port=zeros(T,1);     % to store the accumulative returns   
Port(1)=bill(1)+w'*ReEX(1,:)';        % return in the first period
                                      % notice that ReEX(1,:)) is the 1st row, the traspose makes it a column
CC_Port(1)=1+Port(1);
for t=2:T,
   Port(t)=bill(t) + w'*ReEX(t,:)';       % return at t
   CC_Port(t)=CC_Port(t-1)*(1+Port(t));            % accumulativ return at t
end;
CC_Port_1 = CC_Port;
 
gamma_5 = 5;
w=(1/gamma_5)*V^(-1)*mu; 
Port=zeros(T,1);        % define this 941 vector to store the returns on the portfolio
CC_Port=zeros(T,1);     % to store the accumulative returns   
Port(1)=bill(1)+w'*ReEX(1,:)';        % return in the first period
                                      % notice that ReEX(1,:)) is the 1st row, the traspose makes it a column
CC_Port(1)=1+Port(1);
for t=2:T,
   Port(t)=bill(t) + w'*ReEX(t,:)';       % return at t
   CC_Port(t)=CC_Port(t-1)*(1+Port(t));            % accumulativ return at t
end;
CC_Port_5 = CC_Port;
 
fprintf('******  Accu return of the mkt and the port with risk aversion parameter=1 ********** \n');
CC_RM(T)
CC_Port_1(T)
fprintf('******  Accu return of the mkt and the port with risk aversion parameter=5 ********** \n');
CC_RM(T)
CC_Port_5(T)
