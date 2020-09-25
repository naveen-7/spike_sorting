%% Function to get individual spikes from the raw data
%% Written by naveen at cumc on 6/14/16


function THRESHOLD_FAST_n


clc;
clear all;
close all;


%% Setup directories ---------------

codes_dir = fullfile('E:','NAVEEN_Work','Cerebellum','Codes','CER_codes_NEW');
data_dir  = fullfile('E:','NAVEEN_Work','Cerebellum','Data','RECORDED_CELLS');

cd(data_dir)
disp('!!! THRESHOLD_FILTERING_n has started running !!!');




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

REFRESH=1;

REF(1:length(X))=0; % Refractory period .... does not take the spike multiple times






for k = 1:length(X)
    
    if k==round(length(X)/4)
        disp('25% complete');
    end
    
     if k==round(length(X)/2)
        disp('50% complete');
     end
    
      if k==round((3*length(X))/4)
        disp('75% complete');
    end
    
    % threshold buisness -----
    
    if REF(k)==0
        
        if abs(Y(k))>=DC_OFFSET+(2*DC_DEV)
            
            
            
            
            
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
            Pre_Spikes.values(spk_count,:)= Y(xx-(rate*1):xx+(rate*2.2));
            Pre_Spikes.times(spk_count)= Time_in_ms;
            Waveform = Y(xx-(rate*1):xx+(rate*2.2));
            
            
            
     
            hold on;
            plot(Waveform);
            pause(0.1)
            box off;
            axis off;
            
            REF(k+1:k+(rate*2.2)) =1;
            
        end
        
    end
    
    
end

save(strcat(FileName1(1:12),'_PreSpikes'),'Pre_Spikes');

end






