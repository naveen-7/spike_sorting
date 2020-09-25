%% helper function for getspikes_n 
% written by naveen at JLG on 9/2/19

function F = plotspkgrps_n(WAVES, IND, LIM1, LIM2, Fs)

% WAVES  : matrix of all spike waveforms
% IND    : list of indeces 
% LIM1   : start of spike in ms
% LIM2   : end of spike in ms
% Fs     : sampling frequency (50 kHz)


num = length(unique(IND)); % number of subgroups

F = figure();
clear h;
for i=1:num
    h(i) = subplot(3,3,i);
    hold on;
    plot(-LIM1:1/Fs:LIM2,WAVES(IND==i,:)','color',[0.7 0.7 0.7]);
    plot(-LIM1:1/Fs:LIM2,nanmean(WAVES(IND==i,:)),'color',[1 0 0]);
    xlim([-LIM1 LIM2]);
end
linkaxes(h,'xy');
for i=1:num
    subplot(3,3,i);
    YLIM = ylim;
%     text(LIM2-2,YLIM(2)-0.5,strcat('N = ',num2str(size(WAVES(IND==i,:),1))));
title(strcat('N = ',num2str(size(WAVES(IND==i,:),1))));
text(LIM2-3,YLIM(2)-0.2,num2str(i),'fontsize',20,'color',[0 0 1]);
end












end