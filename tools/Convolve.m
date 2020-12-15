% Convolve Perform a fast frequency domain convolution.
%
% Author: Glenn Dickins
% Date:   26-10-2006
%
% y = Convolve(h,x)
%
% Perform a fast convolution in the frequency domain.  Break the larger
% signal up into blocks of either 16384 samples (reasonably efficient fft) 
% or the length of the smaller signal if longer.
%
% Not really efficient when both signals are short.
% Can run into problems with memory if both signals are really long.
%
% Handles the multidimensional signals if both inputs have the same number
% of columns, or if one of the vectors has only one column.
%

function y = Convolve(h,x)

if (size(x,2) > size(x,1)) x = x.'; h = h.'; trans = 1; else trans = 0; end;
if (length(h) > length(x)) tmp = x; x = h; h = tmp; end;
if (~(size(x,2)==1 || size(h,2)==1 || size(h,2)==size(x,2))) error('Size Mismatch'); end;
N        = max(size(h,1),8192);     % Block length
M        = size(h,1)+size(x,1)-1;   % Amount of valid data at end
h(2*N,1) = 0;                       % Zero pad the impulse response
H        = fft(h);                  % Get the FFT of impulse
y        = zeros(M,max(size(h,2),size(x,2)));                                    

if (size(H,2)<size(x,2)) H = repmat(H,1,size(x,2)); end;
for (n=1:N:size(x,1)+N+2*N-mod(size(x,1)-1,N)-1-N)
    if     (n==1 && n+N>size(x,1))  X = fft([zeros(N,size(x,2)); x(1:size(x,1),:); zeros(N-(size(x,1)),size(x,2)) ]); 
    elseif (n==1)                   X = fft([zeros(N,size(x,2)); x(1:N,:)]);            
    elseif (n+N>size(x,1))          X = fft([x(n-N:size(x,1),:); zeros(2*N-(size(x,1)-n+N+1),size(x,2))]);
    else                            X = fft(x(n-N:n+N-1,:)); end;        
    if (size(X,2)<size(H,2))        X = repmat(X,1,size(H,2)); end;
    Y  = X.*H;
    yb = ifft(Y);
    y(n:min(n+N-1,M),:) = yb(N+1:min(2*N,N+1+(M-n)),:);
end;
if (trans) y = y.'; end;
    