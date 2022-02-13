function [teststat pval] = covtestQsvd(X,batch)
% This function performs equal covariance test for high dimensional Gaussian
% data proposed by Srivastava and Yanagihara(2010)
% H0: Sigma1 = Sigma2 =...=Sigmak.
%
% Input: 1) 'X' is (N by p) data matrix (N < p)
%        2) 'batch' is the label of groups(ex. 111222333...)
%
% Output: 1) 'teststat' is the Qk^2 test statistic.
%         2) 'pval' is the corresponding p-value of the 'teststat'.          

%------------------------- s v d -------------------------------------%
[N p] = size(X);   % original dimension of X
if N > p
    error(strcat('Dimension should be greater than sample size'));
end
[U,D,] = svd(X);   % X = UDV, singular value decomposition(svd)
D2 =D(1:N,1:N);    % trim out zero vectors
R = U*D2;          % 'R' is (N x N) data matrix from svd

%------------------------- identify data-------------------------------%
clear X
X=R;
K = length(unique(batch));
Kname = unique(batch);
Ni=zeros(1,K);
xbatch=batch;
for j=1:K
    batch(xbatch==Kname(j))=j;
    Ni(j)=sum(batch==j); % Ni(j)=size(X(batch==j,:),1)
end

ni = Ni-1;  % degree of freedom
n = sum(ni);

% -------- caculate the test statistic -----------------------%
    clear V V1 V2 V3 a1i a2i
    V = zeros(N,N);
    a1i = zeros(1,K);
    a2i = zeros(1,K);
    for j = 1:K
       Vi = (ni(j))*cov(X(batch==j,:));
       getname = genvarname('V',who);
       eval([getname ' = Vi; ' ])
       V = V + Vi;
       a1i(j) = trace(Vi)/(p*ni(j)) ;
       a2i(j) = (trace(Vi^2) - (1/ni(j))*((trace(Vi))^2))/(p*(ni(j)-1)*(ni(j)+2));
    end
            
a1 = trace(V)/(p*n);
a2 = (trace(V^2) - (1/n)*((trace(V))^2))/(p*(n-1)*(n+2));
%a3 = (1/(n*(n^2+3*n+4)))*((1/p)*trace(V^3)-3*n*(n+1)*p*a2*a1 - n*(p^2)*(a1^3));
a3 = (n/((n-1)*(n-2)*(n+2)*(n+4)))*((1/p)*trace(V^3)-3*(n+2)*(n-1)*a2*a1 - n*(p^2)*(a1^3));
c0 = n*(n^3+6*n^2+21*n+18);
c1 = 2*n*(2*(n^2)+6*n+9);
c2 = 2*n*(3*n+2);
c3 = n*(2*(n^2)+5*n+7);
a4 = (1/c0)*((1/p)*trace(V^4)- p*c1*a1 - (p^2)*c2*(a1^2)*a2 ...
      - p*c3*(a2^2) - n*(p^3)*(a1^4));

         % -------cont...calculate the test statistic --------%
        %--------------------------------------------------------------%
        % Q2 by Srivastava(2010)
        %--------------------------------------------------------------%
        gamma = zeros(1,K);
        Xisq = zeros(1,K);
        for j=1:K
            gamma(j) = a2i(j)/(a1i(j)^2);    
    
            Xia = 4/(ni(j)^2);
            Xib = (a2^2)/(a1^4);
            Xic = 2*ni(j)/p;
            Xid = (a2^3)/(a1^6);
            Xie = (2*a2*a3)/(a1^5);
            Xif = a4/(a1^4);    
            Xisq(j) = Xia*(Xib + Xic*(Xid - Xie + Xif));
            if Xisq(j) < 0
                Xisq(j) = 0.00001;
            end
        end

        %Q2 is N(0,1);
        %Q2 = (gamma(1)-gamma(2))/sqrt(sum(Xisq)); 
        %Q2 = real(Q2);
        %Q2sq = Q2^2;

        % generalized version of Qk ~ chi-square(k)

        gambar = sum(gamma./Xisq)/sum(1./Xisq);

        Qksq = sum(((gamma-repmat(gambar,1,K)).^2)./Xisq);


        % ---------- output --------------------%
        teststat = Qksq;    
             
        pval = 1-chi2cdf(Qksq,K-1);

end % end of function