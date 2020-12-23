function out = MySmooth(x, N)
% operates on first dimension only

% dd = linspace(-1.5, 1.5, N);
% kern = exp(-(dd.^2));

if N==1
    out = x;
    return;
end

Ncol = size(x, 2);
Nel = size(x, 1);

if mod(N, 2)==0
    N = N+1;
end


kern = ones(1, N);
if N>1
    kern(1:floor(N/2)) = 0; %causal
end

kern = kern./sum(kern);



out = zeros(Nel, Ncol);
for j = 1:Ncol
    
    out(:, j) = conv(x(:, j), kern, 'same');
    
end
