%Compute alphas of the CAPM regression

load Indu10_July26_July11.dat;
load Factors_July26_July11.dat      % This file has 5 columns: 
                                
T=1021;
N=10;
 
bill=Factors_July26_July11(1:T,5) / 100;      % The T-bill rate or T-bill return
Re=Indu10_July26_July11(1:T,2:11) / 100;      % the 10  returns 
ReEX=Re - bill*ones(1,N);                      % T by N, the 10 returns in excess of the T-bill rates
 
RmEX=Factors_July26_July11(1:T,2) / 100;        % T by 1, mkt returns in excess of T-bill, the data are such a way                     
  
OneT=ones(T,1);           % a vector of ones 
 
%     Computing alpahs and betas 
 
Alpha=zeros(N,1);            %  to store the N=25 alphas
Beta=zeros(N,1);
AlphaTest=zeros(N,1);        % to store the t-stats for the alphas
BetaTest=zeros(N,1);
R2=zeros(N,1);               % to store the regression R^2
 
for i=1:N,         % to run regression of each asset on the mkt
    
    Y=ReEX(:,i);       % the i-th asset return
    X=[OneT RmEX];            %   the constant and mkt
    results=ols(Y,X);   % Run the OLS regression of Y on X, see the
                           %    free program old.m for details
                           
    Alpha(i)=results.beta(1);   % results.beta is the outpur of ols.m, the regression coefficients
    Beta(i)=results.beta(2); 
    AlphaTest(i)=results.tstat(1);    % results.tstat are the t-stats of the regression coefficients
    BetaTest(i)=results.tstat(2); 
    R2(i)=results.rsqr;
end;

AnnualAlpha = 12*Alpha;
 
fprintf('**Annualized  alpahs** \n');
 
AnnualAlpha
