%% function to get all spikes from a raw continuous recording
%% written by naveen at JLG on 8/21/19


% ch1 MUST be in 50k Hz
% gets the inital threshold (default=mean+6*std; or user defined)
% checks if the spike amplitude is stable or changes over the recording (if stable, constant threshold; else, adaptive)
% find spikes in a window (default: or user defined: )


function getspikes_n


close all;
clear;
clc;

% Setup directories--------------------------------------------------------
codes_dir = fullfile('E:','NAVEEN_Work','Cerebellum','Codes','CER_codes_NEW','LINEAR');
data_dir  = fullfile('E:\NAVEEN_Work\Cerebellum\Data\RECORDED_CELLS');

disp('!!!!!  getspikes_n has started running  !!!!!')
cd(data_dir);


Fs = 50000/1000;  %%%% frequency of signal MUST BE 50000 Hz


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% LOADING THE FILE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % % disp('*************************************************')
% % % disp('*************************************************')
% % % disp('LOAD THE spk2 CH1 FILE')
% % % [FileName1,PathName1] = uigetfile('*.mat','File to Append');   % Open standard dialog box for retrieving files
% % % SPKfile = cat(2,PathName1,FileName1);
% % % disp(strcat('!!!!!','File you entered is :',FileName1,' !!!!!'));
% % % 
% % % A = load(cat(2,PathName1,FileName1));
% % % names = fieldnames(A);
% % % CH1 = A.(names{1,1});
% % % RAW_S = CH1.values;
% % % RAW_T = CH1.times*1000;



disp('*************************************************')
disp('*************************************************')
disp('LOAD raw Spk2 FILE')
[FileName1,PathName1] = uigetfile('*.smr','File to Append');   % Open standard dialog box for retrieving files
disp(strcat('!!!!!','File you entered is :',FileName1,' !!!!!'));
cd(PathName1);
fid=fopen(FileName1);
[data,header]=SONGetChannel(fid,1);
RAW_S = zscore(double(data));
RAW_T = (linspace(header.start,header.stop,header.npoints)*1000)';

disp('RAW waveforms extracted from Spk2 FILE')





%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% GETTING THE initial THRESHOLD %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



THRESHOLD = nanmean(RAW_S)+6*nanstd(RAW_S);
disp(strcat('!-------------- Default threshold is',{' '},num2str(THRESHOLD),{' '},'---------------!'));
disp('!------ if this is ok, press Y else press N and enter your threshold -------!');

figure;
plot(RAW_T(1:0.05*length(RAW_T)),RAW_S(1:0.05*length(RAW_T)));
hold on;
plot(xlim,[THRESHOLD THRESHOLD],'-r')
grid on;
grid minor;

FLAG=0;
iiter = 0;
while FLAG==0
    iiter=iiter+1;
    if iiter==1
        isok = upper(input('Enter your response here: ','s'));
    end
    if isok == 'Y' | isok == 'N'
        FLAG=1;
    else
        disp('----------Please enter a valid input--------');
        disp('---------- Please enter "Y" or "N" ---------');
        isok = upper(input('Here: ','s'));
    end
end

if strcmp(isok,'N')
    THRESHOLD = input('>>>>>> Enter your THRESHOLD here: ');
end



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% CHECKING STABILITY OF RECORDING %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% algorithm: takes NN = 5 random windows throughout the recording and checks for the spike amplitude (maxpoints). 
% if they are not significantly different (ANOVA), then consider it as a stable recording and use constant threshold.
% else, use adaptive threshold. Here, the threshold is 

for ii=1:5
    NN = 5;
    WINDOW = round(5/100*(length(RAW_S)));
    WIN_STRT = randi([WINDOW+1000 length(RAW_S)],NN,1)-WINDOW-1000;
    WIN_STRT = sort(WIN_STRT,'ascend');
    
    for i=1:NN
        SIG = RAW_S(WIN_STRT(i):WIN_STRT(i)+WINDOW);
        %     plot(SIG);
        PEAKS = findpeaks(SIG,Fs,'MinPeakDistance',3);
        PK{i} = PEAKS;
        LEN(i) = length(PEAKS);
    end
    
    maxmax = nanmax(LEN);
    MAT = NaN(maxmax,NN);
    for i=1:NN
        MAT(1:LEN(i),i) = PK{i};
    end
    
    
    P(ii) = anova1(MAT);
    if P(ii)<0.05 P(ii)=0; else P(ii) = 1; end
    delete(gcf);delete(gcf);
    % figure; subplot(2,2,1); errorbar(1:NN,nanmean(MAT),nanstd(MAT)*2); xlim([0 NN+1]);
end

P = nanmean(P);

if P>3/5
    disp('!----------There is no drift in the signal---------!');
    disp('!-----------continuing CONSTANT THRESHOLD----------!');
    is_adapthr = 0;
else
    disp('!-----------There is a drift in the signal-----------!');
    disp('!-----------switching to ADAPTIVE THRESHOLD----------!');
    is_adapthr = 1;
end







%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% GETTING THE ACTUAL SPIKES %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
LIM=8;
THR = THRESHOLD;
BLOCK = round(5/100*length(RAW_T)); % taking blocks of 5%
count=0;
clear WAVES;
BLK_NUM=0;
for i=1:BLOCK:length(RAW_T)-BLOCK-1
    BLK_NUM=BLK_NUM+1
    clear SIG_S SIG_T SIGG LOCS
    SIG_S = RAW_S(i:i+BLOCK);
    SIG_T = RAW_T(i:i+BLOCK);
    
    if is_adapthr==1
        THR = nanmean(SIG_S)+6*nanstd(SIG_S);
    end
    %     plot(RAW_T(i:i+BLOCK),SIG); hold on; plot(xlim,[THR THR],'-r')
    SIGG = SIG_S; SIGG(find(SIG_S<THR))=NaN;
    [~,LOCS] = findpeaks(SIGG,Fs,'MinPeakDistance',0.5);  %% 1 ms apart
    LOCS=LOCS+SIG_T(1);
    
    for kk=1:length(LOCS)
        XX=LOCS(kk);
        time = find(SIG_T>=XX,1);
        if time>Fs*1 & round(time)+(Fs*(LIM+1))<=length(SIG_S)
            temp_spike = SIG_S(round(time)-(Fs*1):round(time)+(Fs*2));   % getting the temp spike 1 ms before to 2 ms after the temp max point
            [~,zero_point]=nanmax(temp_spike);
            Time = zero_point-(Fs*1)+time;
            
            count=count+1;
            WAVES(count,:)=SIG_S(round(Time)-(Fs*1):round(Time)+(Fs*LIM));
        else
            kk=kk+1;
        end
    end
end



    % Plotting to show the waveform ----------------------------------
    plot(-1+1/Fs:1/Fs:size(WAVES,2)/Fs-1,WAVES','color',[0.8 0.8 0.8]);
    disp(strcat('!---- found', {' '},num2str(length(WAVES)),' spikes ----!'));
    
    
    
    
    
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% DOING BASIC kMEANS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
num = 6;
IND = kmeans(WAVES,num);

% try
%     ccount=0;
%     num = 6;
%     IND = kmeans(WAVES,num);
%     for i=1:num
%         for j=1:num
%             ccount=ccount+1;
%             VAL(i,j)=ccount;
%             if i<j
%                 X = WAVES(IND==i,:);
%                 Y = WAVES(IND==j,:);
%                 for k=1:nanmin([Fs*3 size(WAVES,2)])
%                     temp(1,k) = ttest_NN(X(:,k),Y(:,k));
%                     if temp(1,k)<0.05 temp(1,k)=0; else temp(1,k)=1; end
%                 end
%                 val(i,j) = nanmean(temp); clear temp;
%             end
%         end
%     end
% catch
%     ccount=0;
%     num = 6;
%     IND = kmeans(WAVES,num);
%     for i=1:num
%         for j=1:num
%             ccount=ccount+1;
%             VAL(i,j)=ccount;
%             if i<j
%                 X = WAVES(IND==i,:);
%                 Y = WAVES(IND==j,:);
%                 for k=1:nanmin([Fs*3 size(WAVES,2)])
%                     temp(1,k) = ttest_NN(X(:,k),Y(:,k));
%                     if temp(1,k)<0.05 temp(1,k)=0; else temp(1,k)=1; end
%                 end
%                 val(i,j) = nanmean(temp); clear temp;
%             end
%         end
%     end
%     
% end
% 
% cutoff = 0.4;
% val(val==0)=NaN;
% [VX,VY] = (find(val>cutoff));
% 
% list = [VX VY];
% [~,indd] = sort(list(:,1),'descend');
% list = list(indd,:);
% 
% 
% for i=1:length(list)
%     ch = min(list(i,:));
%     IND(find(IND==list(i,1))) = ch;
%     IND(find(IND==list(i,2))) = ch;
% end
% 
% UNI = unique(IND);
% UNI(:,2)=1:length(UNI);
% for i=1:size(UNI,1)
%     IND(IND==UNI(i,1))=UNI(i,2);
% end




figure();
for i=1:num %size(UNI,1)
    h(i)=subplot(3,3,i);
    hold on;
    plot(-1+1/Fs:1/Fs:size(WAVES,2)/Fs-1,WAVES(IND==i,:)','color',[0.7 0.7 0.7]);
    plot(-1+1/Fs:1/Fs:size(WAVES,2)/Fs-1,nanmean(WAVES(IND==i,:)),'color',[1 0 0]);
%     ylim([-2 2]);
    xlim([-1+1/Fs-1 size(WAVES,2)/Fs-1+1]);
end
linkaxes(h,'xy');

disp('!---if this clasification is ok, press "y" else, press "n"---!');

response = upper(input('input here :','s'));
if strcmp(response,'N')
    disp('enter the itemas that belong to the same group')
    disp('Syntax: {{1,2} {3,4,5}} means 1 and 2 belong to the same group and so do 3,4 and 5');
    disp('remeber to enter ALL the items')
    grp = input('input here :');
    
    
    for i=1:length(grp)
        tt = cell2mat(grp{1,i});
        MIN = nanmin(tt) ;
        tt = sort(tt,'descend');
        
        for ii=1:length(tt)
            IND(find(IND==tt(ii))) = MIN;
        end
        
    end
    
    UNI = unique(IND);
UNI(:,2)=1:length(UNI);
for i=1:size(UNI,1)
    IND(IND==UNI(i,1))=UNI(i,2);
end
    
   F =  figure();
    clear h;
    for i=1:size(UNI,1)
        h = subplot(3,3,i);
        hold on;
        plot(-1+1/Fs:1/Fs:size(WAVES,2)/Fs-1,WAVES(IND==i,:)','color',[0.7 0.7 0.7]);
        plot(-1+1/Fs:1/Fs:size(WAVES,2)/Fs-1,nanmean(WAVES(IND==i,:)),'color',[1 0 0]);
        xlim([-1+1/Fs-1 size(WAVES,2)/Fs-1+1]);
    end
    linkaxes(h,'xy');
    
end









end