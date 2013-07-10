function fftplot( x, Fs, flim)
%Power Spectral Density Estimates Using FFT
N = length(x);
xdft = fft(x);
xdft = xdft(1:N/2+1,:);
psdx = (1/(Fs*N)).*abs(xdft).^2;
psdx(2:end-1) = 2*psdx(2:end-1);
freq = 0:Fs/length(x):Fs/2;

plot(freq,10*log10(psdx));
grid on; xlim(flim);
title('Power Spectral Density Estimate');
xlabel('Frequency (Hz)'); ylabel('Power/Frequency (dB/Hz)');

% plot(freq,(psdx));
% grid on; xlim(flim);
% title('Power Spectral Density Estimate');
% xlabel('Frequency (Hz)'); ylabel('Power/Frequency (\muV/Hz)');


% %periodogram
% [psdestx,Fxx] = periodogram(x,rectwin(length(x)),length(x),Fs);
% figure,
% plot(Fxx,10*log10(psdestx));
% % plot(Fxx,psdestx);
% grid on; xlim(flim);
% xlabel('Hz'); ylabel('Power/Frequency (dB/Hz)');
% title('Power Spectral Density Estimate');


% L = length(y);                     % Length of signal
% 
% NFFT = 2^nextpow2(L); % Next power of 2 from length of y
% Y = fft(y,NFFT)/L;
% f = Fs/2*linspace(0,1,NFFT/2+1);

% % Plot single-sided amplitude spectrum.
% plot(f,2*abs(Y(1:NFFT/2+1))) 
% title('Single-Sided Amplitude Spectrum')
% xlabel('Frequency (Hz)')
% ylabel('Amplitude (\muV)')


% % Plot power spectrum. 
% plot(f,10*log10((abs(Y(1:NFFT/2+1))).^2)), xlim(flim)
% title('Power Spectrum')
% xlabel('Frequency (Hz)')
% ylabel('Power (dB)')

end

