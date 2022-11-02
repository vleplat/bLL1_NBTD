function out = pw_kron(X,Y,R)

% Fonction pour réaliser l'opération que je note \odot_p
% entre les facteurs C et A, ou C et B
% pw_kron veut dire "partition-wise Kronecker product"


Lx = size(X,2)/R;
Ly = size(Y,2)/R;

out = [];

for r=1:R
    out = [out kron(X(:,(r-1)*Lx+1:r*Lx),Y(:,(r-1)*Ly+1:r*Ly))];
end

end

