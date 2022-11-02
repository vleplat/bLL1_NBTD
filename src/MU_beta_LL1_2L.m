function [A,B,C,P1,P2,cost] = MU_beta_LL1_2L(Y_1,Y_2,A0,B0,C0,L,P10,P20,P3,options)

%---Init

A = A0; B = B0; C = C0; P1 = P10; P2 = P20;
R = size(C0,2);
cost(1) = 1e10; err = Inf;

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

num = P1'*(((P1*A*H).^(options.beta-2)).*H1)*H' ...
    + options.lambda*(((A*S).^(options.beta-2)).*M1)*S';
denum = P1'*((P1*A*H).^(options.beta-1))*H' ...
    + options.lambda*((A*S).^(options.beta-1))*S';
A = A.*((num./denum).^(options.gamma));
A = max(A,eps);

% Update for B
H = pw_kronL(C,P1*A,R,ones(1,R),L)'; S = pw_kronL(P3*C,A,R,ones(1,R),L)';

num = P2'*(((P2*B*H).^(options.beta-2)).*H2)*H' ...
    + options.lambda*(((B*S).^(options.beta-2)).*M2)*S';
denum = P2'*((P2*B*H).^(options.beta-1))*H' ...
    + options.lambda*((B*S).^(options.beta-1))*S';
B = B.*((num./denum).^(options.gamma));
B = max(B,eps);

% Update for C
H = pw_vecL(A,B,R,L)'; S = pw_vecL(P1*A,P2*B,R,L)';

num = options.lambda*P3'*(((P3*C*H).^(options.beta-2)).*M3)*H' ...
    + (((C*S).^(options.beta-2)).*H3)*S';
denum = options.lambda*P3'*((P3*C*H).^(options.beta-1))*H' ...
    + ((C*S).^(options.beta-1))*S';
C = C.*((num./denum).^(options.gamma)); 
C = max(C,eps);

for r=1:R
    C(:,r) = C(:,r)/norm(C(:,r));
end

% Update for P1
V = A*pw_kronL(C,P2*B,R,ones(1,R),L)';
num = (((P1*V).^(options.beta-2)).*H1)*V';
denum = ((P1*V).^(options.beta-1))*V';
P1 = P1.*((num./denum).^(options.gamma));
P1 = max(P1,eps);

% Update for P1
V = B*pw_kronL(C,P1*A,R,ones(1,R),L)';
num = (((P2*V).^(options.beta-2)).*H2)*V';
denum = ((P2*V).^(options.beta-1))*V';
P2 = P2.*((num./denum).^(options.gamma));
P2 = max(P2,eps);

% % Update for P3
% V = C*pw_vecL(A,B,R,L)';
% num = (((P3*V).^(options.beta-2)).*M3)*V';
% denum = ((P3*V).^(options.beta-1))*V';
% P3 = P3.*((num./denum).^(options.gamma));
% P3 = max(P3,eps);


P = {P1 P2 P3};
cost(i) = objfun(Y_1,Y_2,options.lambda,A,B,C,P,R,L,options.beta);
err(i) = (cost(i-1) - cost(i))/cost(i-1);

if options.verbose
    fprintf('Loop 2 - iter = %i | obj = %d | err = %d (target is %d) \n',i,cost(i),err(i),options.kappa)
end


end



end

