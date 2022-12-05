function [A,B,C,cost] = MU_beta_LL1_1L(Y_1,Y_2,A0,B0,C0,L,P1,P2,P3,options)

%---Init

A = A0; B = B0; C = C0;
R = size(C0,2);
cost(1) = 1e10; err = Inf;
P = {P1 P2 P3};

H1 = tens2mat(Y_1,1,[]); M1 = tens2mat(Y_2,1,[]);
H2 = tens2mat(Y_1,2,[]); M2 = tens2mat(Y_2,2,[]);
H3 = tens2mat(Y_1,3,[]); M3 = tens2mat(Y_2,3,[]);

i = 1;

if options.verbose
    fprintf('Init - iter = %i | obj = %d | err = %d (target is %d) \n',i,cost(1),err,options.kappa)
end
%---Actual loop

while i<options.nIter && err(i)>options.kappa
    
    i = i+1;
    
% Update for A
H = pw_kronL(C,P2*B,R,ones(1,R),L)'; 
S = pw_kronL(P3*C,B,R,ones(1,R),L)';
prod1 = P1*A*H;
prod2 = A*S;
num = P1'*(((prod1).^(options.beta-2)).*H1)*H' ...
    + options.lambda*(((prod2).^(options.beta-2)).*M1)*S';
denum = P1'*((prod1).^(options.beta-1))*H' ...
    + options.lambda*((prod2).^(options.beta-1))*S';
A = A.*((num./denum).^(options.gamma));
A = max(A,eps);

% Update for B
H = pw_kronL(C,P1*A,R,ones(1,R),L)'; S = pw_kronL(P3*C,A,R,ones(1,R),L)';
prod1 = P2*B*H;
prod2 = B*S;
num = P2'*(((prod1).^(options.beta-2)).*H2)*H' ...
    + options.lambda*(((prod2).^(options.beta-2)).*M2)*S';
denum = P2'*((prod1).^(options.beta-1))*H' ...
    + options.lambda*((prod2).^(options.beta-1))*S';
B = B.*((num./denum).^(options.gamma));
B = max(B,eps);

% Update for C
H = pw_vecL(A,B,R,L)'; S = pw_vecL(P1*A,P2*B,R,L)';
prod1 = P3*C*H;
prod2 = C*S;
num = options.lambda*P3'*(((prod1).^(options.beta-2)).*M3)*H' ...
    + (((prod2).^(options.beta-2)).*H3)*S';
denum = options.lambda*P3'*((prod1).^(options.beta-1))*H' ...
    + ((prod2).^(options.beta-1))*S';
C = C.*((num./denum).^(options.gamma)); 
C = max(C,eps);

for r=1:R
    C(:,r) = C(:,r)/norm(C(:,r));
end


cost(i) = objfun(Y_1,Y_2,options.lambda,A,B,C,P,R,L,options.beta);
err(i) = (cost(i-1) - cost(i))/cost(i-1);

if options.verbose
    fprintf('Loop 1 - iter = %i | obj = %d | err = %d (target is %d) \n',i,cost(i),err(i),options.kappa)
end

end


end

