function [ ppg, resp, etrig ] = loadphys( fname, TR, nframes, framesout, fsout )
%load ppg and resp data from physiological data file
%edit from Gary's readphys22.m
%
% Usage:
%   [ ppg, resp, etrig ] = loadphys( fname, TR, nframes, framesout, fsout )
%
% Inputs:
%   fname:      physio file name
%   TR:         TR
%   nframes:    number of scan frames
%   framesout:  frame range to output, e.g. [4 120]
%   fsout:      output sampling frequency
%
% Outputs:
%   ppg:    ppg
%   resp:   res
%   etrig:  ppg trigger event in s
%
% Author: Huang Xiaoshan, xiaoshanhuang@gmail.com
%
% Versions:
%   v0.1:   19-Jun-2013, orignal
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PRESC = 30.0;		% seconds of prepended data
dtr = 0.040;		% in seconds
dth = 0.010;		% in seconds
Twin = 6;           % integration window

fid = fopen(fname, 'r');
[dat nf] = fscanf(fid, '%g\n');
fclose(fid);
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

Tout = 1/fsout;
time1 = time0+framesout(1)*TR;
time2 = time0+(framesout(2)+1)*TR;
nout = fix(time1*fsout):fix(time2*fsout);
fprintf('num output samples = %d\n', length(nout));

respx = max(respr);
respn = min(respr);
resp = 100*(respr - respn)/(respx - respn);

ppgx = max(ppgwave);
ppgn = min(ppgwave);
ppg = 100*(ppgwave - ppgn)/(ppgx - ppgn);

resp = resample(resp, dtr*1000, Tout*1000);
ppg = resample(ppg, dth*1000, 1);
ppg = decimate(ppg, 1000/fsout);

etrig(etrig<time1 | etrig>time2) = [];
etrig = etrig - time1;
resp = resp(nout);
ppg = ppg(nout);

end

