
%% To spike sort SS and CS
%% Written by naveen at CUMC on 6/10/16

% SAMPLING RATE HAS TO BE 50kHz

function CER_SPIKE_SORTING_n



clc;
clear all;
close all;


% Setup directories--------------------------------------------------------

codes_dir = fullfile('e:','NAVEEN_Work','Cerebellum','Codes','CER_codes_NEW','SPIKE_SORTING_n');
data_dir  = fullfile('e:','NAVEEN_Work','Cerebellum','Data','RECORDED_CELLS');


cd(data_dir)
disp('!!! Prelim_n has started running !!!');



%% ONE: LOADING THE FILE --------------------------------------------

% LOADING Spike2 - cell CHANNEL FILE--------------------------------------------------------
disp('*************************************************')
disp('*************************************************')
disp('LOAD a cell channel FILE')
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






close all;
figure
hold on;
plot(CELLfile.values(1:size(CELLfile.values,1),:)')




%%%% also do pca
% allspk = CELLfile.values;
% [wcoeff,~,latent,~,explained] = pca(allspk','algorithm','eig');









FLAG=0;
loopcount=0;

while FLAG==0
    
    loopcount=loopcount+1;
    
    if loopcount>1
        delete(gcf); delete(gcf);
    end
    
    
    % get a valid number of clusters --------------------
    repeat=1;
    Rcount=0;
    while repeat==1;
        
        Rcount=Rcount+1;
        if Rcount>1
            disp(' PLEASE RE-ENTER A VALID RESPONSE ');
        end
        
        disp('ENTER THE NUMBER OF CLUSTERS:');
        
        n = input('HERE: ');
        
        if ~isempty(n)
            repeat=0;
        else
            repeat=1;
        end
        
    end
    % -------------------------------------------------------
    
    
    clear idx C;
    
    [idx,C]=kmeans(CELLfile.values,n);
    
    
    samp_rate = 50000;
    length_of_signal = (size(C,2)/samp_rate)*1000; % in ms
    
    
    Temp_max=0;
    for i=1:n
        [~,temp] = nanmax(C(i,:));
        Temp_max = Temp_max+temp;
    end
    
    Temp_max=Temp_max/n;
    Time_of_zero = (Temp_max/samp_rate)*1000; % in ms
    
    x_values = 1:size(C,2);
    x_values = (x_values/samp_rate)*1000;
    x_values=x_values-Time_of_zero;
    
    
    % Plotting mean of each code -----------
    figure
    hold on;
    for i=1:n
        b(i)= plot(x_values,C(i,:));
        legend(b(1:i));
    end
    xlim([x_values(1) x_values(length(x_values))]);
    box off;
    xlabel('time in ms');
    
    
    % Plotting each code separately -----------
    figure
    for k=1:n
        subplot(3,3,k)
        hold on;
        plot(x_values,CELLfile.values(idx==k,:)')
        xlim([x_values(1) x_values(length(x_values))]);
        ylim([-1 1]);
        box off;
        xlabel('time in ms');
        title(strcat('Code',num2str(k)),'fontsize',7);
    end
    suptitle('All the codes');
    
    
    
    
    
    REPEAT=1;
    rcount=0;
    while REPEAT==1;
        
        rcount=rcount+1;
        if rcount>1
            disp(' PLEASE RE-ENTER A VALID RESPONSE ');
        end
        
        disp(' DO YOU WANT TO RESORT? SAY "Y" OR "N"');
        CHOICE = upper(input('HERE: ','s'));
        
        if strcmp(CHOICE,'Y')
            FLAG=0; REPEAT=0;
        elseif strcmp(CHOICE,'N')
            FLAG=1; REPEAT=0;
        else
            REPEAT=1; FLAG=1;
        end   
        
    end

    
    
end





%% manual adjustments----------



REPEAT=1;
rcount=0;
while REPEAT==1;
    
    rcount=rcount+1;
    if rcount>1
        disp(' PLEASE RE-ENTER A VALID RESPONSE ');
    end
    
    TXT = 'ENTER THE NUMBER OF CODES YOU WANT HERE: ';
    No_codes = input(TXT);
    
    if ~isnumeric(No_codes)
        REPEAT=1;
    else
        REPEAT=0;
    end
    
end




clear Spikes;
Spikes.values = CELLfile.values;
Spikes.times = CELLfile.times;

reshape(Spikes.times,[length(Spikes.times),1]);

for kkk = 1:No_codes
    REPEAT=1;
    rcount=0;
    while REPEAT==1;
        
        rcount=rcount+1;
        if rcount>1
            disp(' PLEASE RE-ENTER A VALID RESPONSE ');
        end
        
        TXT = strcat('ENTER THE CODES FOR CODE#',num2str(kkk),' HERE: ');
        CODE{kkk} = input(TXT);
        
        if ~isnumeric(CODE{kkk})
            REPEAT=1;
        else
            REPEAT=0;
        end
        
    end
end

% 
% REPEAT=1;
% rcount=0;
% while REPEAT==1;
%     
%     rcount=rcount+1;
%     if rcount>1
%         disp(' PLEASE RE-ENTER A VALID RESPONSE ');
%     end
%     
%     TXT = 'ENTER THE CODES FOR CODE#2 HERE: ';
%     CODE2 = input(TXT);
%     
%     if ~isnumeric(CODE1)
%         REPEAT=1;
%     else
%         REPEAT=0;
%     end
%     
% end


for kkk=1:No_codes
    for i=1:length(CODE{kkk})
        temp = CODE{kkk};
        Spikes.codes((find(idx==temp(i))),1)=kkk;
    end
end



% 
% for i=1:length(CODE2)
%     Spikes.codes((find(idx==CODE2(i))),1)=2;
% end




%% saving the data

str = FileName1(1:length(FileName1)-4);
[CC,matches] = strsplit(str,{'_'});
FILENAME = strcat(CC{1},'_',CC(2),'_Spikes.mat');



save(fullfile(PathName1,FILENAME{1,1}),'Spikes');

disp('!!! ALL SPIKES SORTED AND STORED !!!');



end