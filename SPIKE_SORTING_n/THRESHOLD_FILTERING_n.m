%% Function to get individual spikes from the raw data
%% Written by naveen at cumc on 6/14/16


function THRESHOLD_FILTERING_n



clc;
clear all;
close all;


%% Setup directories ---------------

codes_dir = fullfile('E:','NAVEEN_Work','Cerebellum','Codes','CER_codes_NEW');
data_dir  = fullfile('E:','NAVEEN_Work','Cerebellum','Data','RECORDED_CELLS');

cd(data_dir)
disp('!!! THRESHOLD_FILTERING_n has started running !!!');

THR = 3;
refractory = 3; % dont take spikes for the next 3 ms from the current spike

%% ONE: LOADING THE FILE --------------------------------------------

% LOADING Spike2 - cell CH1 FILE--------------------------------------------------------
disp('*************************************************')
disp('*************************************************')
disp('LOAD a cell FILE')
[FileName1,PathName1] = uigetfile('*.mat','File to Append');   % Open standard dialog box for retrieving files

clear A names;
A = load(cat(2,PathName1,FileName1));
names = fieldnames(A);
CELLfile = A.(names{1,1});
cd(PathName1);

str = FileName1(3:length(FileName1)-4);
[CC,matches] = strsplit(str,{'_'});
CHtemp = cell2mat(CC(3));
% CH_number = str2num(CHtemp(3:length(CHtemp)));
FILE_name = strcat(cell2mat(CC(1)),'-',cell2mat(CC(2)));



Sig_val = CELLfile.values;
Sig_tim = CELLfile.times;


rate=50000/1000;

clear Pre_Spike;

figure
X = Sig_tim*1000;  % in ms
Y = Sig_val;

C = [0 0 0];

DC_OFFSET = nanmean(Sig_val);
DC_DEV = nanstd(Sig_val);
hold on;

temp_spike=[];
spk_count=0;

clear h1 h2 h
subplot(5,3,[1 2 3 4 5 6])
box off;
hold on;
LIM1 = 0.005*length(X);
h2 = plot([0 0],[-1 1],'-R','linewidth',1);
ylim([-1 1]);
h1 = plot(X(1:LIM1),Y(1:LIM1),'-','color',C,'linewidth',0.8);
REFRESH=1;

REF(1:length(X))=0; % Refractory period .... does not take the spike multiple times


subplot(5,3,[7 8 10 11 13 14])
h3 = plot(X(1),Y(1),'-','color',C,'linewidth',0.8);
hold on;
h4 = plot([0 0],ylim,'-r','linewidth',2);
ylim([nanmin(Y) nanmax(Y)])
xlabel('time in ms');
ylabel('Microvolts');
hold on;
grid on;




for k = 1:length(X)
    
    
    % slow update -------------------------------------

        subplot(5,1,[1 2])
        if k==1
            box off;
            hold on;
        end

    
    set(h2,'xdata',[X(k) X(k)],'ydata',[-1 1]);
    drawnow update;
    
    
    if REFRESH==1
        xlim([(length(X(1:LIM1))*(REFRESH-1)) X(length(X(1:LIM1))*REFRESH)])
    end
    
    if REFRESH>1
        xlim([X(length(X(1:LIM1))*(REFRESH-1)) X(length(X(1:LIM1))*REFRESH)])
    end
    
    
    if k==(length(X(1:LIM1))*REFRESH)+1
        set(h1,'xdata',X(LIM1*REFRESH:LIM1*(REFRESH+1)),'ydata',Y(LIM1*REFRESH:LIM1*(REFRESH+1)));
        drawnow;
        REFRESH=REFRESH+1;
    end
    
    
    
    
    % threshold buisness -----
    
    if REF(k)==0
        
        if abs(Y(k))>=DC_OFFSET+(THR*DC_DEV)
            
            clear Waveform
            temp_spike = Y(k-(rate*1):k+(rate*2.2));
            spk_count=spk_count+1;
            
%             if sign(Y(k))==1
                [~,zero_point]=nanmax(temp_spike);
%             elseif sign(Y(k))==-1
%                 [~,zero_point]=nanmin(temp_spike);
%             end
            
            Time_in_ms = X(zero_point-(rate*1)+k-1);
            
            xx = zero_point-(rate*1)+k-1;
            Pre_Spike.values(spk_count,:)= Y(xx-(rate*1):xx+(rate*2.2));
            Pre_Spike.times(spk_count)= Time_in_ms;
            Waveform = Y(xx-(rate*1):xx+(rate*2.2));
            
            Pre_Spike.longvalues(spk_count,:)= Y(xx-(rate*1):xx+(rate*8));
            
            
            subplot(5,3,[9 12])
            hold on;
            plot(Waveform);
            box off;
            axis off;
            
            REF(k+1:k+(rate*refractory)) =  1;
            
        end
        
    end
    
    
    
    
    
    
    
    %% JUST FOR REPRESENTATION PURPOSES ------------
    
    
    subplot(5,3,[7 8 10 11 13 14])
    LIM2 = 1000;
    
    
    if k>350
        set(h3,'xdata',X(k-350:k+350),'ydata',Y(k-350:k+350))  % showing 14 seconds of data
        drawnow;
        xlim([X(k-350) X(k+350)])
        set(h4,'xdata',[X(k) X(k)],'ydata',ylim)  % showing 14 seconds of data
    end
    
    hold on;
    plot(xlim,[DC_DEV*THR+DC_OFFSET DC_DEV*THR+DC_OFFSET],'-b')
    plot(xlim,[DC_OFFSET-(DC_DEV*THR) DC_OFFSET-(DC_DEV*THR)],'-b');
    
    
    
    
end


end


