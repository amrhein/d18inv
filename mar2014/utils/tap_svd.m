function [xest] = tap_svd(U,S,V,y,k,gamma)
%function [xest] = tap_svd(U,S,V,y,k,gamma)
% tapered SVD solver. gamma is the tapering parameter. gamma squared gives
% the expected value of the variance at every place and time.

Vk = V(:,1:k);
Uk = U(:,1:k);
Sk = S(1:k,1:k);

% xest =  (Vk*(Sk\Uk'))*y(:);
%xest =  (Vk*Sk*((Sk.^2+gamma^2*eye(k))\Uk'))*y(:);

% weight matrix of gamma and Sk
sm = (Sk.^2+gamma^2*eye(k))\Sk;
%sm = Sk\(Sk.^2+gamma^2*eye(k));
xest =  Vk*(sm * Uk')*y(:);


