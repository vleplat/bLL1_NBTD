function [A,B,cost] = mu_nmf(X,A0,B0,options)

%---Init

A = A0; B = B0;
cost(1) = 1e10; diff_cost(1) = 1e10;

n = 1;

%---Actual loop

while n<options.nIter && diff_cost(n)>options.kappa
    
    n = n+1;
    
% Update for A

num = ((A*B').^(options.beta-2).*X)*B;
denum = (A*B').^(options.beta-1)*B;
A = A.*((num./denum).^(options.gamma));

% Update for B

num = A'*((A*B').^(options.beta-2).*X);
denum = A'*((A*B').^(options.beta-1));
B = B'.*((num./denum).^(options.gamma));
B = B';

cost(n) = frob(X - A*B','squared');
diff_cost(n) = (cost(n-1) - cost(n))/cost(n-1);



end

