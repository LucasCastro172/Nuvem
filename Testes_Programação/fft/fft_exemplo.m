Fs = 3500;                   % samples per second
dt = 1/Fs;                   % seconds per sample
StopTime = 0.35;             % seconds
t = (0:dt:StopTime-dt)';     % seconds
%%Sine wave:
Fc = 95;                     % hertz
x = cos(2*pi*Fc*t);
% Plot the signal versus time:
subplot(2,1,1);
plot(t,x);
xlabel('Time(S)');
title('Signal versus Time');
subplot(2,1,2);
f = Fs*linspace(0, 1/2, floor(numel(x)/2));
fr = fft(x);
fr = abs(fr)/numel(x);
fr = fr(1:floor(end/2));
fr(2:end-1) = 2*fr(2:end-1);
plot(f, fr);
xlabel('Time(S)');
title('Frequency Response');