function [x0, lb, ub, ank] = initialscmopt(psd,freq)
% initialize the starting point and bounds for fmincon
% Input
%        psd --- power spectrum density in natural scale
%        freq --- frequency bins
% ouput
%         x0 --- starting point
%          lb --- lower bounds
%         ub --- upper bounds

% see also initialparascm
% Shiang Hu, Jul. 2018

% check if psd and freq in column
if size(psd,1)==1, psd=psd'; end
if size(freq,1)==1, freq=freq'; end

% paras
nk = 15; % default maximum number of components
x0 = zeros(4*nk,1);

% vis
maxpsd = max(psd);
% figure, subplot(2,(nk+1)/2,1), plot(freq,psd); title('Org');
% xlabel('Frequency (Hz)'), ylabel('PSD (uv^2/Hz)');ylim([-0.2 1]*maxpsd);

% first extremes
[~,fmi] = fstextrm(freq,psd,maxpsd);
if isempty(fmi), lb=0; ub=0; ank=0; return; else, fma = freq(1); end % add Xi


for i=1:nk
    [~,maxloc]=min(abs(freq-fma));
    [~,minloc]=min(abs(freq-fmi));
    rou =psd(maxloc);  % rou
    mu =freq(maxloc);  % mu
    
    [~,phloc] = min(abs(psd-0.5*rou)); % peak half
    tau=abs(freq(phloc)-mu);   % tau
    
    slope = -(psd(maxloc)-psd(minloc))/(freq(maxloc)-freq(minloc));
    nu= abs((log(slope)-log(rou)))/log(1+1/tau); % nu
    
    x0(4*i-3:4*i)=[rou mu tau nu];
    
    psd_det =  rou*(1+((freq-mu)/tau).^2).^(-nu);
    psd = psd - psd_det;
    
    %     subplot(2,(nk+1)/2,i+1), plot(freq,[psd,psd_det]); title(strcat('-',num2str(i),'pk'));
    %     xlabel('Frequency (Hz)'), ylabel('PSD (uv^2/Hz)');ylim([-0.2 1]*maxpsd);
    [fma,fmi] = fstextrm(freq,psd,maxpsd);
    if isempty(fma), break; end
end
% x1 = zeros(60,1);
% x1(1:20)=x0;
% x0=x1;
% rou
idx=1:4:57;lb(idx)=0.25*x0(idx); ub(idx)=0.98*x0(idx);
% mu
idx=2:4:58;lb(idx)=0.999*x0(idx); ub(idx)=1.001*x0(idx);
% tau
idx=3:4:59;lb(idx)=0.06*x0(idx); ub(idx)=100*x0(idx);
% nu
idx=4:4:60;lb(idx)=0.25*x0(idx); ub(idx)=40*x0(idx);

x0 = reshape(x0,4,nk);
[~,odr] = sort(x0(1,:),2,'descend');
[~,loc] = min(x0(2,x0(1,:)~=0));
odr = [loc setdiff(odr, loc, 'stable')];

x0 = x0(:,odr);
lb = reshape(lb,4,nk); lb = lb(:,odr);
ub = reshape(ub,4,nk); ub = ub(:,odr);
ank = sum(x0(1,:)~=0);
end