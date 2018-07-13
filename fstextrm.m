function [fma,fmi,ppSpline,p] = fstextrm(freq,psd,maxpsd)
% Return the first fma and fmi in current remained psd curve
% Input
%        freq --- frequency bins
%        psd --- current remained psd curve in the fitting process
%        maxpsd --- maximum value of the original spectrum
% Output
%        fma --- frequency bin of maximum psd
%        fmi --- frequency bin of minimum psd

% Shiang Hu, Jul. 2018

[ppSpline,p] = csaps(freq,psd); % p=1 for simulation
[fma, fmi] = splineMaximaMinima(ppSpline);
if isempty(fmi), return; end
% s_smooth = csaps(freq,psd,p,freq);

fma(ppval(ppSpline,fma)<0.08*maxpsd)=[];
if isempty(fma), fmi=fmi(1); return; else, fma=fma(1); end

a = ppval(ppSpline,fma); % extreme maximum
fmi(a-ppval(ppSpline,fmi)<0.25*a)=[];
fmi = fmi(fmi-fma>0); 
fmi = fmi(1);

% visualization
% figure, plot(freq,[s_smooth,psd]); hold on;
% plot(fma,ppval(ppSpline,fma),'r*','linewidth',2);
% plot(fmi,ppval(ppSpline,fmi),'gs','linewidth',2);
% legend({'Fit','Truth','Pks','Tfs'}); set(gca,'fontsize',12);
% xlabel('Freq'); ylabel('PSD');

end