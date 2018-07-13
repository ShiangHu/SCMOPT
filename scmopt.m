% single channel spctrum components model

clean; dbstop if error; addpath(genpath(cd));
load psdf;
[nf,nc] = size(psdall);
options = optimoptions(@fmincon,'display','notify-detailed','diagnostics','off');
% options.OptimalityTolerance=1e-45;
options.MaxFunctionEvaluations=1e5;
% options.UseParallel=1;

lmd=0;
sigall=zeros(nf,nc) ;
parall=zeros(60,nc);

% parpool(210);
tic;
parfor chn=1:nc
    psd = psdall(:,chn);
    [x0, lb0, ub0, ank] = initialscmopt(psd,freq);
    if ank==0, continue; end
    
    abiclh = zeros(ank,3);
    sigma = zeros(nf,ank);
    xm = zeros(60,ank);
    xs = zeros(size(x0));    lbs = zeros(size(lb0));    ubs = zeros(size(ub0));
    
    for i=1:ank
        xs(:,i) = x0(:,i);
        lbs(:,i) = lb0(:,i); lb = lbs(:);
        ubs(:,i) = ub0(:,i);ub = ubs(:);
        
        xm(:,i) = fmincon(@(x)scmobj(psd,freq,x,lmd),xs(:),[],[],[],[],lb,ub,[],options);
        [abiclh(i,1),sigma(:,i), abiclh(i,2), abiclh(i,3)] = scmobj(psd,freq,xm(:,i),lmd);
        
        xs = reshape(xm(:,i),4,15);
        disp(i);
    end
    
    [~,loc]=min(abiclh(:,2:3)); k = max(loc)+2;
    if k>ank,k=ank;end
    sigall(:,chn)=sigma(:,k) ;
    parall(:,chn)=xm(:,k);
    fprintf('---fitting the %2d channel---\n',chn);
    
end
toc;
save fitrs sigall parall;