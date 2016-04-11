%I use the daily S&P 500 index from 1980 to 2015 to verify the existence of power law in the tail distribution of stock price movement. 
%The estimation results suggest that the power law exists in both tails. 
%In particular I cannot reject that the power equals to three, as suggested by previous studies which use different or older data sets. (Gabaix, 2011; Gopikrishnan et al., 1999)

clc;
clear all;
 
load SP1980currentDailyReturn.mat;
 
date = SP1980currentDailyReturn(:,1);   %date:YYYYMMDD
sprn = SP1980currentDailyReturn(:,2);   %S&P500 return
spindx = SP1980currentDailyReturn(:,3); %S&P500 index
N = 9017;
 
%define  return
for i = 2:N;
    G(i) = log(spindx(i)) - log(spindx(i-1)); 
end; 
 
%define volatility
v = sqrt(var(G));
 
%define  normalized return
for i = 2:N;
    g(i) = (G(i) - mean(G))/v;
end
 
figure
plot(g)
set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', [])
xlim([0 N])
xlabel('Time')
ylabel('Normalized S&P 500 return, daily')
 
 
%create interval for probability distribution calculation
gmax = max(g);
gmin = min(g);
range = max(g) - min(g);
 
s = 200;
int(1) = gmin;
for i = 2:s;
    int(i) = int(i-1) + range/s;
end;
 
%fraction in each interval
g_sort = sort(g);
count = zeros(s,1);
 
for i = 1:s;
    for j = 1:N;
        if g_sort(j) < int(i);
            count(i) = count(i) + 1;
        end;
    end;
end;
 
%estimate pdf
density(1) = (s/range)*count(1)/N;  %(range/s) is the interval length 
for i = 2:s;
    density(i) =  (s/range)*(count(i) - count(i-1))/N;
end
 
%plot data cdf and normal cdf as a comparison
figure
plot(int, density);
xlabel('Normalized S&P 500 return')
ylabel('Density')
xlim([-20 20])
hold on;
y = -20:0.1:20;
mu = 0;
sigma = 1;
f = exp(-(y-mu).^2./(2*sigma^2))./(sigma*sqrt(2*pi));
plot(y,f)
legend({'Data' 'Standard normal distribution'});
hold off;
 
%estimate cdf
cumden = count/N;
plot(int, cumden);
 
rev_cumden = 1 - cumden;
plot(int, rev_cumden);
 
%Upper Tail
up = find(int> 2, 1); %return that lies 2 std above the mean
plot(int(up:end),rev_cumden(up:end))
 
ones = ones(s-up+1, 1);
Y1 = log(rev_cumden(up:end));
X1 = [ones log(int(up:end)')];
b1 = inv(X1'*X1)*X1'*Y1;
e1 = Y1 - X1*b1; %e is residual
 
[n1,k1] = size(X1);
sig21=(e1'*e1)/(n1-k1); %bias-corrected estimate of error variance
xx1=inv(X1'*X1);
v1=xx1*sig21;
s1=sqrt(diag(v1)); % (Homoskedastic)  standard error of least-sqaures estimator
 
figure
plot(log(int(up:end)),log(rev_cumden(up:end)),'*')
hold on;
plot(log(int(up:end)), X1*b1)
xlabel('log(normalized S&P 500 return)')
ylabel('log(cumulative distribution)')
legend({'Data' 'Model' });
hold off;
 
%Lower Tail
low = find(int>-2, 1); %return that lies 2 std below the mean
plot(abs(int(1:low -1)),cumden(1:low -1))
plot(log(abs(int(1:low -1))),cumden(1:low -1))
plot(log(abs(int(1:low -1))),log(cumden(1:low -1)))
 
clear ones;
m = find(cumden> 0, 1); %find the row index such that cdf is nonzero (to avoid infinite value after taking log)
ones = ones(low -m -1, 1);
Y2 = log(cumden(1+m:low -1));
X2 = [ones log(abs(int(1+m:low -1)'))];  %Since domain of log-function must be positive, take absolute value 
b2 = inv(X2'*X2)*X2'*Y2;
e2 = Y2 - X2*b2; %e is residual
 
[n2,k2] = size(X2);
sig22 = (e2'*e2)/(n2-k2); %bias-corrected estimate of error variance
xx2 =inv(X2'*X2);
v2=xx2*sig22;
s2=sqrt(diag(v2)); % (Homoskedastic) standard error of least-sqaures estimator
 
figure;
plot(log(abs(int(1+m:low-1)')),log(cumden(1+m:low -1)),'o')
hold on;
plot(log(abs(int(1+m:low-1)')), X2*b2)
xlabel('log(normalized S&P 500 return)')
ylabel('log(cumulative distribution)')
legend({'Data' 'Model' });
hold off;
 
%Combine both tails
figure
plot(log(int(up:end)),log(rev_cumden(up:end)),'*')
hold on;
plot(log(int(up:end)), X1*b1)
plot(log(abs(int(1+m:low-1)')),log(cumden(1+m:low -1)),'o')
plot(log(abs(int(1+m:low-1)')), X2*b2, '--')
xlabel('log(normalized S&P 500 return)')
ylabel('log(cumulative distribution)')
legend({'Data - upper tail' 'Model- upper tail' 'Data - lower tail' 'Model- lower tail'});
hold off;
 
%coefficients & confidence interval for coefficient
b1(2)
b2(2)
CI1 = [b1(2)-1.96*s1(2),  b1(2)+1.96*s1(2)]
CI2 = [b2(2)-1.96*s2(2), b2(2)+1.96*s2(2)]
 
%measures of fit
sig2y1 = (Y1 - mean(Y1))'*(Y1 - mean(Y1))/n1;
sig2y2 = (Y2 - mean(Y2))'*(Y2 - mean(Y2))/n2;
R2_1 = 1 - sig21/sig2y1  %R squares for upper tail model
R2_2 = 1 - sig22/sig2y2  %R squares for lower tail model
