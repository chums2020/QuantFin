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
