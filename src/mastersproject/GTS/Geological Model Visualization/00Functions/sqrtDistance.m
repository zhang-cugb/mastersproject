function D = sqrtDistance(X, Y, Z)
D = (X'-X).^2+(Y'-Y).^2+(Z'-Z).^2;
D = D.^.5;
D = triu(D);