function [A,B,C,cost_hist] = MU_beta_LL1_1L_gpu(Y_1,Y_2,A0,B0,C0,L,P1,P2,P3,options)

%---Init

gpu = gpuDevice();
fprintf('Using a %s .\n', gpu.Name)

A = A0; B = B0; C = C0;
R = size(C0,2);
P = {P1 P2 P3};
rate = 1; 
cost_prev = 1e20;
cost_hist = [];

H1 = tens2mat(Y_1,1,[]); M1 = tens2mat(Y_2,1,[]);
H2 = tens2mat(Y_1,2,[]); M2 = tens2mat(Y_2,2,[]);
H3 = tens2mat(Y_1,3,[]); M3 = tens2mat(Y_2,3,[]);

n = 1;

gpuA = gpuArray(A);
gpuB = gpuArray(B);
gpuC = gpuArray(C);
gpuH1 = gpuArray(H1);
gpuM1 = gpuArray(M1);
gpuH2 = gpuArray(H2);
gpuM2 = gpuArray(M2);
gpuH3 = gpuArray(H3);
gpuM3 = gpuArray(M3);
gpuP1 = gpuArray(P1);
gpuP2 = gpuArray(P2);
gpuP3 = gpuArray(P3);
gpuOnes = gpuArray(ones(1,R));

%---Actual loop

while n<options.nIter && rate>options.kappa
    
    n = n+1;
    
% Update for A
H = pw_kronL(gpuC,gpuP2*gpuB,R,gpuOnes,L)'; 
S = pw_kronL(gpuP3*gpuC,gpuB,R,gpuOnes,L)';
prod1 = gpuP1*gpuA*H;
prod2 = gpuA*S;

num = gpuP1'*(((prod1).^(options.beta-2)).*gpuH1)*H' ...
    + options.lambda*(((prod2).^(options.beta-2)).*gpuM1)*S';
denum = gpuP1'*((prod1).^(options.beta-1))*H' ...
    + options.lambda*((prod2).^(options.beta-1))*S';
gpuA = gpuA.*((num./denum).^(options.gamma));
gpuA = max(gpuA,eps);

% Update for B
H = pw_kronL(gpuC,gpuP1*gpuA,R,gpuOnes,L)'; S = pw_kronL(gpuP3*gpuC,gpuA,R,gpuOnes,L)';
prod1 = gpuP2*gpuB*H;
prod2 = gpuB*S;

num = gpuP2'*(((prod1).^(options.beta-2)).*gpuH2)*H' ...
    + options.lambda*(((prod2).^(options.beta-2)).*gpuM2)*S';
denum = gpuP2'*((prod1).^(options.beta-1))*H' ...
    + options.lambda*((prod2).^(options.beta-1))*S';
gpuB = gpuB.*((num./denum).^(options.gamma));
gpuB = max(gpuB,eps);

% Update for C
H = pw_vecL(gpuA,gpuB,R,L)'; S = pw_vecL(gpuP1*gpuA,gpuP2*gpuB,R,L)';
prod1 = gpuP3*gpuC*H;
prod2 = gpuC*S;

num = options.lambda*gpuP3'*(((prod1).^(options.beta-2)).*gpuM3)*H' ...
    + (((prod2).^(options.beta-2)).*gpuH3)*S';
denum = options.lambda*gpuP3'*((prod1).^(options.beta-1))*H' ...
    + ((prod2).^(options.beta-1))*S';
gpuC = gpuC.*((num./denum).^(options.gamma)); 
gpuC = max(gpuC,eps);

for r=1:R
    gpuC(:,r) = gpuC(:,r)/norm(gpuC(:,r));
end


if mod(n, 100) == 0 && options.verbose
    A = gather(gpuA);
    B = gather(gpuB);
    C = gather(gpuC);
    cost_up = objfun(Y_1,Y_2,options.lambda,A,B,C,P,R,L,options.beta);
    rate = abs(cost_prev-cost_up)/cost_prev;
    fprintf('Loop 1 - iter = %i | obj = %d | err = %d (target is %d) \n',n,cost_up,rate,options.kappa)
    cost_prev = cost_up;
    cost_hist = [cost_hist cost_up];
end


end


end

