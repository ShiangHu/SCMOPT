
clean;
freq = 0.3914*(1:204);
psdall = importdata('./iEEGgtm/results.mat');
idx = importdata('./iEEGgtm/precentralgyrus.mat');

figure, plot(freq,psdall.spect(idx(68),:))

for i=1:12:123
psd = log(psdall.spect(idx(i:i+11),:)');
figure,
for j=1:12
subplot(3,4,j),plot(freq,(psd(:,j)),'linewidth',2); title(num2str(i+j-1)),xlabel('Frequency Hz'), ylabel('PSD uv^2/Hz');
end
end