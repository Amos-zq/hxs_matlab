function fftplot( y, Fs, flim)
L = length(y);                     % Length of signal

NFFT = 2^nextpow2(L); % Next power of 2 from length of y
Y = fft(y,NFFT)/L;
f = Fs/2*linspace(0,1,NFFT/2+1);

% % Plot single-sided amplitude spectrum.
% plot(f,2*abs(Y(1:NFFT/2+1))) 
% title('Single-Sided Amplitude Spectrum')
% xlabel('Frequency (Hz)')
% ylabel('Amplitude (\muV)')


% Plot power spectrum. 
plot(f,((abs(Y(1:NFFT/2+1))))), xlim(flim)
title('Power Spectrum')
xlabel('Frequency (Hz)')
ylabel('Power (dB)')

end

