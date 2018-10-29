function M = MI(X,Y)
X = double(X);
Y = double(Y);
X_norm = X - min(X(:)) +1; 
Y_norm = Y - min(Y(:)) +1;
matAB(:,1) = X_norm(:);
matAB(:,2) = Y_norm(:);
h = accumarray(int32(matAB+1), 1); % joint histogram
hn = h./sum(h(:)); % normalized joint histogram
y_marg=sum(hn,1); 
x_marg=sum(hn,2);
Hy = - sum(y_marg.*log2(y_marg + (y_marg == 0))); % Entropy of Y
Hx = - sum(x_marg.*log2(x_marg + (x_marg == 0))); % Entropy of X
arg_xy2 = hn.*(log2(hn+(hn==0)));
h_xy = sum(-arg_xy2(:)); % joint entropy
M = Hx + Hy - h_xy; % mutual information
end