%% function to get all spikes from a raw continuous recording
%% written by naveen at JLG on 8/21/19
%% v4 written by naveen at JLG on 9/2/19

% helper functions required:
% plotspkgrps_n
% mergespkgrps_n



% ch1 MUST be in 50k Hz
% gets the inital threshold (default=mean+6*std; or user defined)
% checks if the spike amplitude is stable or changes over the recording (if stable, constant threshold; else, adaptive)
% finds spikes in a window (default: or user defined: )


function getspikes_v4_n


close all;
clear;
clc;
% Setup directories--------------------------------------------------------
codes_dir = 'E:\NAVEEN_Work\Cerebellum\Codes\CER_codes_NEW\LINEAR';
data_dir= 'E:\NAVEEN_Work\Cerebellum\Data\RECORDED_CELLS';
% data_dir= 'E:\NAVEEN_Work\Cerebellum\Data\CS_CELLS\RECORDED_CS';

disp('!!!!!  getspikes_n has started running  !!!!!')
cd(data_dir);

Fs = 50000/1000;  %%%% frequency of signal MUST BE 50000 Hz


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% LOADING THE FILE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% disp('*************************************************')
% disp('*************************************************')
% disp('LOAD THE spk2 CH1 FILE')
% [FileName1,PathName1] = uigetfile('*.mat','File to Append');   % Open standard dialog box for retrieving files
% SPKfile = cat(2,PathName1,FileName1);
% disp(strcat('!!!!!','File you entered is :',FileName1,' !!!!!'));
% 
% A = load(cat(2,PathName1,FileName1));
% names = fieldnames(A);
% CH1 = A.(names{1,1});
% RAW_S = CH1.values;
% RAW_T = CH1.times*1000;

CON_FAC = 6.5534e+03; %(if you open channel through son in matlab)

disp('*************************************************')
disp('*************************************************')
disp('LOAD raw Spk2 FILE')
[FileName1,PathName1] = uigetfile('*.smr','File to Append');   % Open standard dialog box for retrieving files
disp(strcat('!!!!!','File you entered is :',FileName1,' !!!!!'));
cd(PathName1);
fid=fopen(FileName1);
[data,header]=SONGetChannel(fid,1);
% RAW_S = zscore(double(data));
RAW_S = (double(data)/CON_FAC);
RAW_T = (linspace(header.start,header.stop,header.npoints)*1000)';

disp('RAW waveforms extracted from Spk2 FILE')


[data_T,header_T]=SONGetChannel(fid,5);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% GETTING THE initial THRESHOLD %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
THRESHOLD = nanmean(RAW_S)+6*nanstd(RAW_S);
disp(strcat('!-------------- Default threshold is',{' '},num2str(THRESHOLD),{' '},'---------------!'));
disp('!------ if this is ok, press Y else press N and enter your threshold -------!');

f(1) = figure;
plot(RAW_T(1:0.05*length(RAW_T)),RAW_S(1:0.05*length(RAW_T)));
hold on;
plot(xlim,[THRESHOLD THRESHOLD],'-r','linewidth',2)
grid minor;
xlim([1 round(RAW_T(round(0.05*length(RAW_T))))]);

FLAG=0;
iiter = 0;
while FLAG==0
    iiter=iiter+1;
    if iiter==1
        isok = upper(input('Enter your response here: ','s'));
    end
    if isok == 'Y'
        FLAG=1;
    elseif isok == 'N'
        THRESHOLD = input('>>>>>> Enter your THRESHOLD here: ');
        %replot after entering the new threshold---------------------
        f(iiter+1) = figure;
        plot(RAW_T(1:0.05*length(RAW_T)),RAW_S(1:0.05*length(RAW_T)));
        hold on;
        plot(xlim,[THRESHOLD THRESHOLD],'-r','linewidth',2)
        grid minor;
%         xlim([1 round(RAW_T(0.05*length(RAW_T)))])
        %reprompt----------------------------------------------------
        disp('!------ if this is ok, press Y else press N and enter your threshold -------!');
        isok = upper(input('Here: ','s'));
    else
        disp('----------Please enter a valid input--------');
        disp('---------- Please enter "Y" or "N" ---------');
        isok = upper(input('Here: ','s'));
    end
end


close(f(1:size(f,2)-1));
clear f;

if THRESHOLD<0
    RAW_S = RAW_S*-1;
end;
THRESHOLD = abs(THRESHOLD);


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

LIM1=1; LIM2=4;   %%% spk waveform will be from LIM1 to LIM2
LIM3 = 15; %%% for later CS extention 
THR = THRESHOLD;
BLOCK = round(5/100*length(RAW_T)); % taking blocks of 5%
count=0;
clear WAVES;
BLK_NUM=0;
TOTBLK = round((length(RAW_T)-BLOCK-1)/BLOCK);
% disp(strcat('TOTAL NUMBER OF BLOCKS = ',{' '},num2str(TOTBLK)));
disp('>>>>>>>> GETTING SPIKES IN BLOCKS >>>>>>>>>>>');

for i=1:BLOCK:length(RAW_T)-BLOCK-1
    BLK_NUM=BLK_NUM+1;
    
    if round(BLK_NUM/TOTBLK,2)>=0.22 & round(BLK_NUM/TOTBLK,2)<=0.28 disp('~~~~~ 25% complete ~~~~~'); end
    if round(BLK_NUM/TOTBLK,2)>=0.47 & round(BLK_NUM/TOTBLK,2)<=0.52 disp('~~~~~ 50% complete ~~~~~'); end
    if round(BLK_NUM/TOTBLK,2)>=0.72 & round(BLK_NUM/TOTBLK,2)<=0.78 disp('~~~~~ 75% complete ~~~~~'); end
    if round(BLK_NUM/TOTBLK,2)>=0.92 & round(BLK_NUM/TOTBLK,2)<=0.97 disp('~~~~~~ almost done ~~~~~'); end
    
    clear SIG_S SIG_T SIGG LOCS
    SIG_S = RAW_S(i:i+BLOCK);
    SIG_T = RAW_T(i:i+BLOCK);
    
    if is_adapthr==1
        THR = nanmean(SIG_S)+6*nanstd(SIG_S);
    end
    %     plot(RAW_T(i:i+BLOCK),SIG); hold on; plot(xlim,[THR THR],'-r')
    SIGG = SIG_S; SIGG(find(SIG_S<THR))=NaN; 
%     SIGG(find(SIG_S> nanmean(SIG_S)+12*nanstd(SIG_S)))=NaN;  %%%%%% to remove artifact
%     SIGG(find(SIG_S< nanmean(SIG_S)-10*nanstd(SIG_S)))=NaN;  %%%%%% to remove artifact
    [~,LOCS] = findpeaks(SIGG,Fs);  %% 0.5 ms apart
    
    LOCS=LOCS+SIG_T(1);
    
    for kk=1:length(LOCS)
        XX=LOCS(kk);
        time = find(SIG_T>=XX,1);
        if time>Fs*LIM1 & round(time)+(Fs*(LIM3+1))<=length(SIG_S)
            temp_spike = SIG_S(round(time)-(Fs*LIM1):round(time)+(Fs*2));   % getting the temp spike 1 ms before to 2 ms after the temp max point
            [~,zero_point]=nanmax(temp_spike);
            Time = zero_point-(Fs*LIM1)+time;
            
            count=count+1;
            WAVES(count,:) = SIG_S(round(Time)-(Fs*LIM1):round(Time)+(Fs*LIM2));
            WAVES2(count,:)= SIG_S(round(Time)-(Fs*LIM1):round(Time)+(Fs*LIM3));
            TIMES(count,:) = SIG_T(round(Time));
        else
            kk=kk+1;
        end
    end
end



% Plotting to show the waveform ----------------------------------
figure
subplot(2,2,1)
plot(-LIM1:1/Fs:LIM2,WAVES','color',[0.8 0.8 0.8]);
xlim([-LIM1 LIM2]);
subplot(2,2,2)
plot(-LIM1:1/Fs:LIM3,WAVES2','color',[0.8 0.8 0.8]);
xlim([-LIM1 LIM3]);

disp(strcat('!---- found', {' '},num2str(length(WAVES)),' spikes ----!'));
% disp(strcat('!---- removed', {' '},num2str(REMOVED),' spikes ----!'));




%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% DOING BASIC kMEANS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



num = 9;
IND = kmeans(WAVES,num);
plotspkgrps_n(WAVES2,IND,LIM1,LIM3,Fs);


%%%%% MERGE & PLOT ----------------------------------------------
disp('enter the items that belong to the same group')
disp('Syntax: {{1,2} {3,4,5}} means 1 and 2 belong to the same group and so do 3,4 and 5');
disp('remember to enter ALL the items')


IND = mergespkgrps_n(IND);
F(1) = plotspkgrps_n(WAVES2,IND,LIM1,LIM3,Fs);


SATISFIED = 0;
while SATISFIED ==0
    
    disp('!---if this clasification is ok, press "Y" else, press "N"---!');
    response = upper(input('input here :','s'));
    while ~(response == 'Y' |  response == 'N')
        response = upper(input('Please enter a valid input : ','s'));
    end
    
    if response == 'Y'
        SATISFIED=1;
        break;
    end
    
    
    %%%%% RESORT ----------------------------------------------
    if response == 'N'
        resortgrp = input('enter the first group that you want to resort :');
        while ~(isnumeric(resortgrp) & length(resortgrp)==1 & any(IND(:) == resortgrp) )
            disp('!--- PLEASE ENTER ONE GROUP ONLY IN NUMERIC VALUE ---!');
            resortgrp = input('enter the first group that you want to resort :');
        end
        
        
        number= input('>>>>>> How many clusters: ');
        TEMP = WAVES(find(IND==resortgrp),:);
        INDDD=kmeans(TEMP,number);
        
        tempp= find(IND==resortgrp);
        limit=20;
        for i=1:length(INDDD)
            INDDD(i,:)= INDDD(i,:)+limit;
        end
        
        for iii=1:length(INDDD)
            IND(tempp(iii)) = INDDD(iii);
        end
        UNI = unique(IND);
        UNI(:,2)=1:length(UNI);
        for i=1:size(UNI,1)
            IND(IND==UNI(i,1))=UNI(i,2);
        end
        
        F(2) = plotspkgrps_n(WAVES2,IND,LIM1,LIM3,Fs);
        
        disp('enter the items that belong to the same group')
        disp('Syntax: {{1,2} {3,4,5}} means 1 and 2 belong to the same group and so do 3,4 and 5');
        disp('remember to enter ALL the items')
        
        try
            IND = mergespkgrps_n(IND);
        catch
            disp('!!! ERROR DETECTED !!!');
            IND = mergespkgrps_n(IND);
        end
        F(3) = plotspkgrps_n(WAVES2,IND,LIM1,LIM3,Fs);
        
    end
%     delete(gcf); delete(gcf);
    clear F;
end







%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% ORDER ALL GROUPS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


disp('!--------PLEASE ORDER THE CODES SO THAT #1 is SS; #2 is CS and #3-n are others---------------!')
disp('!--------syntax: [1 2 3]---------------!')
order = input('Please enter the order here  : ');

while ~(length(order) == length(unique(IND)))
    disp('!!!!----ERROR ENCOUNTERED----!!!!');
    disp('!--------PLEASE ORDER THE CODES SO THAT #1 is SS; #2 is CS and #3-n are others---------------!')
    disp('!--------syntax: [1 2 3]---------------!')
    order = input('Please enter the order here  : ');
end

xx = (reshape(sort(order,'ascend'),1,[])==reshape(sort(unique(IND),'ascend'),1,[]));

while ~( nanmean(xx(:))==1 )
    
    disp('!!!!----ERROR ENCOUNTERED----!!!!');
    disp('!--------PLEASE ORDER THE CODES SO THAT #1 is SS; #2 is CS and #3-n are others---------------!')
    disp('!--------syntax: [1 2 3]---------------!')
    order = input('Please enter the order here  : ');
    xx = (reshape(sort(order,'ascend'),1,[])==reshape(sort(unique(IND),'ascend'),1,[]));
end

% replace indeces a.k.a order it
clear ind;
ind(:,1)=IND;
for ii=1:length(order)
    %     ind(find(IND==ii),2)=order(ii);
    ind(find(IND==order(ii)),2)=ii;
end
IND=ind(:,2); clear ind;


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SAVING FILES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

F = figure();
s1 = subplot(3,3,1); hold on;
plot(-LIM1:1/Fs:LIM2,WAVES(IND==1,:)','color',[0.7 0.7 0.7]);
plot(-LIM1:1/Fs:LIM2,nanmean(WAVES(IND==1,:)),'color',[1 0 0]);
xlim([-LIM1 LIM2]);
s2 = subplot(3,3,2); hold on;
plot(-LIM1:1/Fs:LIM3,WAVES2(IND==2,:)','color',[0.7 0.7 0.7]);
plot(-LIM1:1/Fs:LIM3,nanmean(WAVES2(IND==2,:)),'color',[1 0 0]);
xlim([-LIM1 LIM3]);
linkaxes([s1 s2],'y');

%%% SIMPLE spikes ----------
S_SPK.times = TIMES(IND==1)/1000;
S_SPK.values = WAVES(IND==1,:);
S_SPK.codes = ones(size(S_SPK.times));
SSfile = strcat(PathName1(1:end-4),FileName1(1:end-4),'_SS_nu');
save(SSfile,'S_SPK');

%%% COMPLEX spikes ----------
C_SPK.times = TIMES(IND==2)/1000;
C_SPK.values = WAVES2(IND==2,:);
C_SPK.codes = ones(size(C_SPK.times));
CSfile = strcat(PathName1(1:end-4),FileName1(1:end-4),'_CS_nu');
save(CSfile,'C_SPK');

%%% TRIGGER --------------------
Trigg.times = data_T;
TRIGGfile = strcat(PathName1(1:end-4),FileName1(1:end-4),'_TRIGG');
save(TRIGGfile,'Trigg');




disp('!---- all spikes are sorted and saved -----!');


end

