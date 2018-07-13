function [lhsm, sigma, aic, bic]=scmobj(psd,freq,s,lmd)
% create the SCM objective function to utilize opt toolbox
% fisher scoring algorithm to fit single channel spetrum components model
% generate the approximation of information matrix and the score equation
% Input
%         s ---- model parameters [rou, mu, tau, nu], peak amplitude,
%         position, peak half width, peak exponent
%         freq ---- frequency bins, predictors
%         psd ---- power spectrum density for a single channel (in log scale)
%         lmd ---- positive scalar, regularization paramater for L1 norm components
%         amplitudes
% Output
%         U ---- score equations
%      info ---- approximation to fisher information matrix (Demidenko mixed models pp86)

% Shiang Hu, Jul. 2018

nf = length(freq);
ns = 26;  % # of segments
nw = 3.5; % # of slepian windows
s = s(:);
sigma = zeros(nf,1);

% spectrum reconstruction with maximum 15 components
for f=1:nf
    omg = freq(f);
    sigma(f) = stc(s(1:4),omg)+stc(s(5:8),omg)+stc(s(9:12),omg)+stc(s(13:16),omg)+stc(s(17:20),omg)+...
        stc(s(21:24),omg)+stc(s(25:28),omg)+stc(s(29:32),omg)+stc(s(33:36),omg)+stc(s(37:40),omg)+...
        stc(s(41:44),omg)+stc(s(45:48),omg)+stc(s(49:52),omg)+stc(s(53:56),omg)+stc(s(57:60),omg);
end

sigma(sigma==0) = sigma(sigma==0) + eps;
rou = s(1:4:57);

lhsm = sum(log(sigma)+psd./sigma) + lmd*sum(rou);

aic = 2*lhsm+2*sum(s~=0);
bic = 2*lhsm+log(nf*ns*nw*2)*sum(s~=0);
end

function sigma=stc(s,omg)
% student t curve
% s --- the parameter set [rou, mu, tau, nu]
% omg --- frequency bin
% sigma --- spectrum

rou = s(1); mu = s(2); tau = s(3); nu = s(4);
sigma = rou*(1+((omg-mu)/tau)^2)^(-nu);
end