Fs = 48000; 
SweepLength  = 2*Fs;
SweepGap = 0;
SweepSamples = 2*Fs;  
LinearTo = 100;  
OverRun = 141; 
FadeOut = 330;

LinearTime  = floor(((SweepSamples)/(1+log((Fs/2+OverRun)/LinearTo))));               % SweepLength of linear sweep
PhaseDiff   = [linspace(0,LinearTo,LinearTime) logspace(log10(LinearTo),log10(Fs/2+OverRun),SweepSamples-LinearTime)]'/Fs;
Sweep       = sin(2*pi*cumsum(PhaseDiff));
Sweep       = Sweep.*[ones(1,SweepSamples-FadeOut) cos((0.5:FadeOut)/FadeOut/2*pi).^2]';
SweepI      = ifft(1./fft([zeros(50000,1); Sweep; zeros(50000,1)])); 
SweepI      = SweepI(50000+(-1999:SweepSamples+8000));
Offsets     = SweepGap/2 + SweepSamples+2000-round(log([ 1 2:5])*(SweepSamples-LinearTime)/(log(Fs/2+OverRun)-log(LinearTo)));  % Offsets for the non linearities

Sweep       = [ zeros(SweepGap/2,1); Sweep; zeros(SweepGap/2,1) ];
    
    