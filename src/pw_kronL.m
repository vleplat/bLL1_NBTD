function out = pw_kronL(X,Y,R,Lx,Ly)

% Fonction pour réaliser l'opération que je note \odot_p
% entre les facteurs C et A, ou C et B
% avec des L_r différents
% pw_kron veut dire "partition-wise Kronecker product"


out = kron(X(:,1:Lx(1)),Y(:,1:Ly(1)));

for r=2:R
    out = [out kron(X(:,sum(Lx(:,1:r-1))+1:sum(Lx(:,1:r))),Y(:,sum(Ly(:,1:r-1))+1:sum(Ly(:,1:r))))];
end

end

