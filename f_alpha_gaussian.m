function x = f_alpha_gaussian(Nx,seed)
%Generate 1/f noise, reference:  
% https://ccrma.stanford.edu/~jos/sasp/Example_Synthesis_1_F_Noise.html
%
% Syntax:  
%     Nx: number of samples to synthesize
%     seed: rand seed
%
% Inputs:
%     
%
% Outputs:
%     
%
% Example:
%     
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
%
% Versions:Synthesis
%    v0.1: , orignal
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

B = [0.049922035 -0.095993537 0.050612699 -0.004408786];
A = [1 -2.494956002   2.017265875  -0.522189400];
nT60 = round(log(1000)/(1-max(abs(roots(A))))); % T60 est.
rng(seed,'twister');
v = randn(1,Nx+nT60); % Gaussian white noise: N(0,1)
x = filter(B,A,v);    % Apply 1/F roll-off to PSD
x = x(nT60+1:end);    % Skip transient response
x = x./std(x);
