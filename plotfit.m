% [~,sigmav] = scmobj(psd,freq,x);


figure,
for i=1:ank
subplot(2,round(ank/2),i),
plot(freq,([psd,sigma(:,i)]),'linewidth',2); 
title(strcat('Fitted cnpmts:',num2str(i)));
xlabel('Freq'); ylabel('PSD');
legend({'Real',strcat('Fit',num2str(i))});
end
% 
% figure,subplot(211),plot(log10(lmdv),'linewidth',2); xlabel('300 LMDs'); ylabel('log10 of LMD values');
% set(gca,'fontsize',12); subplot(212),plot(lmdv,'linewidth',2); xlabel('300 LMDs'); ylabel('LMD values');
% set(gca,'fontsize',12);
%     
% figure, plot(log10(lmdv),[sum(xm(1:4:57,:)~=0); sum(xm(1:4:57,:))],'linewidth',2);
% xlabel('log10 of 300 lmds'); ylabel('sum of rou'); legend({'# of cnpmts','sum of rou'});
% set(gca,'fontsize',12);

% figure,plot(log10(lmdv),[aic,bic],'linewidth',2);xlabel('log10 of 300 lmds'); 
% ylabel('criteria'); legend({'AIC','BIC'});set(gca,'fontsize',12);

% figure, plot(freq,[psd,sigmav],'linewidth',2), legend({'Real','Fit'}); 
% xlabel('Frequency Hz'); ylabel('PSD uv^2/Hz');set(gca,'fontsize',12);
[~,loc1]=min(aic);
[~,loc2]=min(bic);
if loc1-loc2~=0
    warning('Aic ~= Bic');
end
rou = xm(1:4:57,loc2);
% fprintf('lmd = %0.4f\n',lmdv(loc2));
figure, plot(freq,[psd,sigma(:,loc2)],'linewidth',2);title(strcat('Fitted copmts: ',num2str(sum(rou~=0))));
xlabel('Freq (Hz)'); ylabel('PSD (uV^2/Hz)'); set(gca,'fontsize',12);legend({'Real','Fit'});
fprintf('fitting_error = %0.4f\n',rssq(sigma(:,loc2)-psd)./rssq(psd));
% fprintf('min_lh = %4.2f\n',fval);
% figure('Name','Fitted paras vs. Initialized paras'),
% subplot(2,2,1), scatter(x0(1:4:57),x(1:4:57),'filled'), title('Rou');xlabel('Initialized');ylabel('Fitted');
% subplot(2,2,2), scatter(x0(2:4:58),x(2:4:58),'filled'), title('Mu');xlabel('Initialized');ylabel('Fitted');
% subplot(2,2,3), scatter(x0(3:4:59),x(3:4:59),'filled'), title('Tau');xlabel('Initialized');ylabel('Fitted');
% subplot(2,2,4), scatter(x0(4:4:60),x(4:4:60),'filled'), title('Nu');xlabel('Initialized');ylabel('Fitted');


% fprintf('PkAmp = \n')
% fprintf('%4.2f    %4.2f\n',[x0(1:4:57),x(1:4:57)]);
% fprintf('PkLoc = \n')
% fprintf('%4.2f    %4.2f\n',[x0(2:4:58),x(2:4:58)]);
% fprintf('PHBW = \n')
% fprintf('%4.2f    %4.2f\n',[x0(3:4:59),x(3:4:59)]);
% fprintf('PkExp = \n')
% fprintf('%4.2f    %4.2f\n',[x0(4:4:60),x(4:4:60)]);

