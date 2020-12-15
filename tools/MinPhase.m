function h = MinPhase(Fb,H,Fs,N)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MinPhase
%
% Calculate a minimum phase response given a specific set of frequency magnitude pairs
%

if (Fb(1)  >0   ) Fb = [ 0; Fb(:)    ]; H = H ([1 1:end],  :,:,:); end
if (Fb(end)<Fs/2) Fb = [ Fb(:); Fs/2 ]; H = H ([1:end end],:,:,:); end

F = linspace(0,Fs/2,N+1);
for (n=1:size(H,2))
    g = spline(Fb,20*log10(H(:,n)),F)';
    G(:,n) = 10.^([ g; g(end-1:-1:2,:) ]/20);
end;

h = real(ifft(exp(conj(hilbert(log(G))))));
h = h(1:N,:);
