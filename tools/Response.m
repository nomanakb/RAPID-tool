% Response
%
% Take an impulse response and generate a perceptually banded set of powers
%
% Y = Response(X, Fs, T, Fb, fade)
% X     the data or set of data in colums to calculate the response for
% T     the time to pad to
% Fs    the sample rate
% Fb    the band centers
% fade  the time of X to fade out or [ in out ]
%

function [ Y M ] = Response(x, Fs, T, Fb, fade, M)
Dims  = size(x);
if (length(fade)<1) fade = [0 0]; end;
if (length(fade)<2) fade = [0 fade]; end;
fade = ceil(fade * Fs);

if (~isempty(T) && T*Fs<length(x)) x = x(1:T*Fs,:,:,:,:,:,:); end;
if (sum(fade)>size(x,1)) fade = floor(fade/sum(fade)*size(x,1)); end;
Window = [ sin((0.5:fade(1))/fade(1)*pi/2).^2 ones(1,size(x,1)-sum(fade)) cos((0.5:fade(2))/fade(2)*pi/2).^2]';
x = x .* Window;
if (~isempty(T) && T*Fs>length(x)) x(T*Fs,:,:,:,:,:,:,:,:) = 0; end;
N = size(x,1);

F  = (0:floor(N/2)-1)/N*Fs;

if (nargin<6 || size(M,1)~=length(Fb) || size(M,2) ~= length(F))
    M  = zeros(length(Fb),length(F));
    for (b=1:length(Fb)-1)
        C1 = find(min(abs(Fb(b)-F))  ==abs(Fb(b)-F),1);
        C2 = find(min(abs(Fb(b+1)-F))==abs(Fb(b+1)-F),1);
        if (b==1)            M(1,1:C1  )=1; end;
        if (b==length(Fb)-1) M(b+1,C2:end)=1; end;
        if (C1~=C2)
            w  = cos(((C1:C2)-C1)/(C2-C1)*pi/2).^2;
            M(b,  C1:C2) = w;
            M(b+1,C1:C2) = 1-w;
        else
            M(b,C1) = 1;
        end;
    end;
end;
X      = fft(x(:,:));
X      = X(1:floor(N/2),:);
Y      = reshape((M*abs(X).^2)./sum(M,2),length(Fb),size(x,2),size(x,3),size(x,4),size(x,5),size(x,6),size(x,7));


