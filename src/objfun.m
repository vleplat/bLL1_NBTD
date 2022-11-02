function f = objfun(Y_1,Y_2,lambda,A,B,C,P,R,L,beta)
%objfun Compute the objective function value
Y_1_chap = zeros(size(Y_1));
Y_2_chap = zeros(size(Y_2));
idx = 1;
for r=1:R
    A_r = A(:,idx:idx+L-1);
    B_r = B(:,idx:idx+L-1);
    c_r = C(:,r);
    Y_1_chap = Y_1_chap + outprod(((P{1}*A_r)*(P{2}*B_r)'),c_r);
    Y_2_chap = Y_2_chap + outprod((A_r*B_r'),P{3}*c_r);
    idx = idx + L;
end
%f = betadiv(Y_1,Y_1_chap,beta) + lambda*betadiv(Y_2,Y_2_chap,beta);
f = betaDiv2(Y_1+eps,Y_1_chap+eps,beta) + lambda*betaDiv2(Y_2+eps,Y_2_chap+eps,beta);
end