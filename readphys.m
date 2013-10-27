function [ phi_resp_amp, phi_resp_trig] = readphys( fname, nframes, TR )
%READPHYS Read physiological data file and output ppg and respiratory data
%
% Syntax:  
%     [ phi_resp_amp] = readphys( fname, nframes, TR )
%
% Inputs:
%     fname   : *.physio file name
%     nframes : number of scan frames
%     TR      : TR
%
% Outputs:
%     phi_resp_amp  : ecg triggers's respiratory phase based on amplitude histogram method
%     phi_resp_trig : ecg triggers's respiratory phase based on respiratory triggers
%
% Example:
%     [ phi_resp_amp] = readphys( '*.physio', 295, 2.04 );
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% Author: Xiaoshan Huang, xiaoshanhuang@gmail.com
%
% Versions:
%    v0.1: 2013-09-26 11:32, edited from gary's readphys22.m and catie's retroicor_main.m
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

kern2 = [19:-1:-19];

PRESC = 30.0;		% seconds of prepended data
dtr = 0.040;		% in seconds
dth = 0.010;		% in seconds
Twin = 6;		% integration window

%  open the input file

% fname = input('gimme file = ', 's');
fid = fopen(fname, 'r');
[dat nf] = fscanf(fid, '%g\n');
fclose(fid);

% foo = input('gimme [nframes TR(s)] = ');
% nframes = foo(1);  TR = foo(2);
Tscan = nframes*TR;

% parse the triggers and resp waveform

ntrig = find(dat == -9999)-1;
etrig = dat(1:ntrig)*dth;	% trigger times
ppgloc = find(dat == -8888)-1;
respr = dat(ntrig+2:ppgloc);	% resp waveform
ppgwave = dat(ppgloc+2:end);	% ppg waveform

nresp = length(respr);
ne = length(etrig);
averr = (etrig(ne) - etrig(1))/ne;
avehr = 60/averr;		% BPM
fprintf('average HR = %d BPM\n', fix(avehr));

% check for missing or too many trigs

faverr = .35*averr;

trig(1) = etrig(1);
k = 1;
for j=2:ne
    k = k +1;
    trig(k) = etrig(j);
    if(abs(etrig(j) - etrig(j-1) - averr) < faverr)
        break;
    end
end
i = j;
while(i<=ne)
    dif = etrig(i) - trig(k-1) - averr;
    if(abs(dif) < faverr)
        trig(k) = etrig(i);
        %fprintf('i k dif trigk ok %d %d %f %f\n', i, k, dif, trig(k));
        k = k + 1;
        i = i + 1;
    elseif (dif > faverr)
        trig(k) = trig(k-1) + averr;
        %fprintf('i k dif trigk long %d %d %f %f\n', i, k, dif, trig(k));
        k = k + 1;
    else
        %fprintf('i k dif short %d %d %f\n', i, k, dif);
        i = i + 1;
    end
end
etrig = trig;
ne = k - 1;

ts = nresp*dtr;	% sampled time
time0 = ts - Tscan;		% how much to ignore
fprintf('start time = %.3f\n', time0 - PRESC);
if(time0 < 0)
    fprintf('Eek!  Not enough data- I give up!\n');
    return
end

% % Tout = input('output sample interval, s [2s]) = ');
% if(isempty(Tout))
%   Tout = 2;
% end
% nout = fix(Tscan/Tout);
% time = (0:nout-1)*Tout;
% fprintf('num output samples = %d\n', nout);
% 
% % find ecg intvl for each cardiac cycle and thence hrate
% 
% clear hrate;
% for j=1:nout
%   t = time0 + time(j);
%   i1 = 1; i2 = 1;
%   for i1=1:ntrig
%     if(etrig(i1)>=t-Twin*.5)
%       break;
%     end
%   end
%   for i2=i1:ntrig
%     if(etrig(i2)>t+Twin*.5)
%       break;
%     end
%   end
%   i2 = i2 - 1;
%   if(i2 == i1)              % end of trace
%     i1= i1 - 1;
%   end
%   hrate(j) = (i2-i1)*60/(etrig(i2) - etrig(i1));       % bpm 
% end
% 
% subplot(2,2,1); 
% plot(time, hrate);grid
% ylabel('hrate, BPM');
% xlabel('time, s');
% title(fname);
% 
% % get HRV
% 
% hrv = [diff(hrate) 0];
% subplot(2,2,2); 
% plot(time, hrv);grid
% ylabel('HRV, arb units');
% xlabel('time, s');

%  now do the resp.  find the peaks

respx = max(respr);
respn = min(respr);
resp = 100*(respr - respn)/(respx - respn);
n1 = fix(time0/dtr);
% subplot(2,2,2); 
% plot((1:nresp-n1+1)*dtr, resp(n1:nresp));grid
% ylabel('respiration'); 
% xlabel('time, s');

drdt = conv(resp, kern2);
drdt = drdt(19:nresp+18);
d2rdt2 = conv(drdt, kern2);
d2rdt2 = d2rdt2(19:nresp+18);
rpeak = (d2rdt2 > 0.5e-6);  	% nice threshold

% find the resp trigs

nr = 0; 
for (j=2:nresp)
  if (rpeak(j)==1 & rpeak(j-1)==0)  % first only
    nr = nr + 1;
    rtrig(nr) = j*dtr;
  end
end
averesp = (rtrig(nr) - rtrig(1))/nr;
averrate = 60/averesp;		% breaths/min
fprintf('average resp = %d breaths/min\n', fix(averrate));

% % find resp intvl for each breath
% 
% resptime = diff(rtrig);
% nrespt = nr - 1;
% tresptime = (rtrig(1:nrespt) + rtrig(2:nr))*.5;
% respintvl = spline(tresptime, resptime, time+time0);
% rrate = 60./respintvl;
% subplot(2,2,3); 
% plot(time, rrate);grid
% ylabel('resp rate, BrPM');
% xlabel('time, s');

% % plot rv(t)
% 
% clear rv;
% for j = 1:nout
%   t = time0 + time(j);
%   i1 = fix((t - Twin*.5)/dtr);
%   i2 = min(nresp, fix(t + Twin*.5)/dtr);
%   rv(j) = std(resp(i1:i2));
% end
% subplot(2,2,4); 
% plot(time, rv);  grid
% ylabel('RV')
% xlabel('time, s');

% % save the output
% 
% fnout = input('output file [cr for default, s for none]= ', 's');
% if(isempty(fnout))
%   fnout = sprintf('%s.out',fname);
%   elseif(strcmp(fnout, 's'))
%   return;
% end
% fout = fopen(fnout, 'w');
% fprintf(fout, '%f  %f %f\n', [hrate', rrate', rv']');
% fclose(fout); 
% fprintf('wrote file  %s\n', fnout);


% bin respiration signal into 100 values
[Hb,bins] = hist(resp,100);
% calculate derivative
% first, filter respiratory signal - just in case
f_cutoff = 1; % max allowable freq
fs = 1/dtr;
wn = f_cutoff/(fs/2);
ntaps = 20;
b = fir1(ntaps,wn);
respfilt = filtfilt(b,1,resp);
drdt = diff(respfilt);

phi_resp_trig = [];
phi_resp_amp = [];
phase_times = etrig(etrig>time0);
for i = 1:length(phase_times)
    % respiration: method based on triggers
    prev_trigs = find(rtrig<phase_times(i));
    if isempty(prev_trigs)
        t1 = 0;
    else
        t1 = rtrig(prev_trigs(end));
    end
    next_trigs = find(rtrig>phase_times(i));
    if isempty(next_trigs)
        t2 = nframes*TR;
    else
        t2 = rtrig(next_trigs(1));
    end
    phi_resp_trig(i) = (phase_times(i) - t1)*2*pi/(t2-t1);
    
    % respiration: method based on amplitude histogram
    iphys = max(1,round(phase_times(i)/dtr)); % closest index in resp waveform
    iphys = min(iphys,length(drdt));
    amp = resp(iphys);
    dbins = abs(amp-bins);
    [blah,thisBin] = min(dbins);  %closest resp histo bin
    numer = sum(Hb(1:thisBin));
    phi_resp_amp(i) = pi*sign(drdt(iphys))*(numer/length(respfilt));
end


end

