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

%%

%Is the CAPM rejected by the data or not
 
X = [OneT RmEX];            %  T by 2, observations on the regressors, one and the mkt
B=[Alpha'; Beta'];               %  stack the regression coefficients together to 
                               % compute all the model residuals in matrix form quickly
                            % One can of course compute them asset by asset
                            % and then stack them together as E below
E = ReEX - (X*B);                   % matrix of model residuals, T by N
Sigma=(E'*E)/T;                     % Residual Cov estimator
 
           
                  % GRS test of the hypothesis that all alphas are zero
                  
 mu_F=mean(RmEX)';
 var_F=cov(RmEX,1);         % computing the varinace by dividing the sum by T, not (T-1).
 
Theta2=mu_F*mu_F  / var_F;       % The sharo ratio squared of the factor
  
 Q=Alpha'*Sigma^(-1)*Alpha;
 GRS=((T-N-1)/N)*Q/(1+Theta2);
 c=fdis_cdf(GRS,N,T-N-1);     % This function compute the F-distribution which replaces Matlab's 
                                %  function, fcdfH(GRS,N,T-N-1), not available in the basic version
 GRSp=1-c;                    % P-value of the test decided by the F-distribution
 
fprintf('**The GRS test and its P-value** \n');
 
[GRS GRSp]
