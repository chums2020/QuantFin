
load Indu10.dat
                           % note: the indutry data is in percentage format
Dates=Indu10(1:947,1);         % single column vector of the dates
Return=Indu10(1:947,2:11)/100;      % 10 columns of the industry returns
                                                 % divide by 100 as the data are in percentage points
                                                 
load Factors_July1926_May05.dat;      % This file has 5 columns:   date, mkt excess return, size factor, book/ME factor, T-bill rate
T=947;
RmEX=Factors_July1926_May05(1:T,2) / 100;      %  Get the returns on mkt in excess of T-bill rate, divide by 100 as the data is in percentage
bill=Factors_July1926_May05(1:T,5) / 100;      % The T-bill return
 
N=10;
ReEX=Return - bill*ones(1,N);  % the excess industry return
Mu=mean(ReEX);               % the sample mean or mu estimator,  
Mu=Mu';                       % put Mu into a column vector
V=cov(ReEX);                % the sample covariance matrix o V estimator       
 
T=947;         % # of period or observations
 
[BB D]=eig(V);         % B are the eigenvectors, and the diagonal elements of D are the egevalues
 
Lam=zeros(N,N);        % define the dimension of Lam to store the egevalues in its diagonal elements
A=zeros(N,N);           % to store the eigenvectors in the notation of the notes
 
for i=1:N,              % Matlab ordered the eigenvalues in increasing order, we reverse the order here
    Lam(i,i)=D(N-i+1,N-i+1);       %  to follow the notes and other authors
    A(:,i)=BB(:,N-i+1);              % Doing the same:  change the last column to the 1st, etc.
end;
 
A1=A(:,1);             % The coefficient for the first PCA, the first column of A, 
                          % i.e.,  the eigenvector corresponding to the largest eigenvalue 
 
A2=A(:,2);             % that for the 2nd PCA
 
RR=ReEX - ones(T,1)*Mu';       % de-meaned excess return as PCA requires zero mean of the variables
 
f=RR*A1;                      %  The realizations of the 1st PCA over time, T by 5
                              %  the PCA analysis is usually applied to de-meaned data, but 
                              % this has no impact on A1, etc, but does on the compuation 
                               % of PCA, i.e., PCA = A1*demeaned variables, so that the factor mean is zero.                           
                                
                               
%     Let us now compare PCA and the mkt index in explaining the returns
 
          %  Analyze the mkt index factor first
 
RM2=zeros(N,1);               % to store the regression R^2
OneT=ones(T,1);           % a vector of ones
 
for i=1:N,         % to run regression of each asset on the mkt
    
    Y=ReEX(:,i);       % the i-th asset return
    X=[OneT RmEX];            %   the constant and mkt
    results=ols(Y,X);   % Run the OLS regression of Y on X, see the
                           %    free program old.m for details                   
    RM2(i)=results.rsqr;
end;
                               
        % Now analyze the PCA factor
 
Rf2=zeros(N,1);               % to store the regression R^2
fEX=f - bill;                  % the excess return of the PC, minics mkt excess return
for i=1:N,         % to run regression of each asset on the mkt
    
    Y=ReEX(:,i);       % the i-th asset return
    X=[OneT fEX];            %   the constant and mkt
    results=ols(Y,X);   % Run the OLS regression of Y on X, see the
                           %    free program old.m for details                   
    Rf2(i)=results.rsqr;
end;
 
fprintf('**eigenvalues and eigenvectors of the sample covariance matrix of the industry data** \n');
Lam
A
 
fprintf('**coefficients for the first PCA** \n');
A1
 
fprintf('**R-squared of the regression of the industry returns on the first PCA** \n');
Rf2(1)
 
fprintf('******   The average R-squared for the market index and PCA  ********** \n');
[mean(RM2)   mean(Rf2)]
