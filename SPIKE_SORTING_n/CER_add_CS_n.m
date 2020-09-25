
%% Part of CER_SPIKE_SORTING_n code series.
%% This function adds a vector of dimentions 1 x  to the datafile CS_Library
%% It also allows to classify CS based on number of spikelets

%% Written by naveen at cumc on 6/13/16

% Recording frequency HAS TO BE 50kHz


function CER_add_CS_n

clc;
clear all;
close all;


disp('!!! CER_add_CS_n has started running !!!');


% Setup directories--------------------------------------------------------

codes_dir = 'e:\NAVEEN_Work\Cerebellum\Codes\CER_codes_NEW\SPIKE_SORTING_n';
Lib_dir  = 'e:\NAVEEN_Work\Cerebellum\Codes\CER_codes_NEW\SPIKE_SORTING_n\Libraries';
data_dir  = fullfile('e:','NAVEEN_Work','Cerebellum','Data','RECORDED_CELLS');


% LOAD LIBRARY ------------------------------

cd(Lib_dir)
FILENAME = strcat('ComplexSpikes_n');
tempwave = load(strcat(Lib_dir,'\',FILENAME));


% LOADING Spike2 - cell CHANNEL FILE--------------------------------------------------------

cd(data_dir)
disp('*************************************************')
disp('*************************************************')
disp('LOAD a cell Ch1 FILE')
[FileName1,PathName1] = uigetfile('*.mat','File to Append');   % Open standard dialog box for retrieving files

clear A names;
A = load(cat(2,PathName1,FileName1));
names = fieldnames(A);
CELLfile = A.(names{1,1});
cd(PathName1);

SIGNAL = CELLfile.values;
TIMES = CELLfile.times*1000;
rate = 50000/1000;

END_FLAG=0;


GLOBAL_MEAN = nanmean(SIGNAL);
GLOBAL_STD = nanstd(SIGNAL);


while END_FLAG==0
    
    
    TXT = strcat('ENTER THE TIME IN MS HERE: ');
    time = find(TIMES>=input(TXT),1);
    
    
    clear Waveform
    temp_spike = SIGNAL(round(time)-(rate*1):round(time)+(rate*10));   % getting the temp spike 1 ms before to 10 ms after the temp max point
    [~,zero_point]=nanmax(temp_spike);
    
    Time = zero_point-(rate*1)+time;
    Time_in_ms = TIMES(Time);
    
    Waveform =SIGNAL(round(Time)-(rate*1):round(Time)+(rate*10));
    
    
    % Plotting to show the waveform ----------------------------------
    
    delete(gcf);
    figure
    plot(Waveform);
    hold on;
    plot(xlim,[GLOBAL_MEAN GLOBAL_MEAN])
    plot(xlim,[GLOBAL_MEAN-2*GLOBAL_STD GLOBAL_MEAN-2*GLOBAL_STD],'-k')
    plot(xlim,[GLOBAL_MEAN+2*GLOBAL_STD GLOBAL_MEAN+2*GLOBAL_STD],'-k')
    box off;
    axis off;
    
    
    
    WAVEFORM = reshape(Waveform,[1,length(Waveform)]);
    
    REPEAT=1;
    rcount=0;
    while REPEAT==1;
        
        rcount=rcount+1;
        if rcount>1
            disp(' PLEASE RE-ENTER A VALID RESPONSE ');
        end
        
        TXT = strcat('IS THIS SPIKE OK? ENTER Y or N: ');
        OK = upper(input(TXT,'s'));
        
        if strcmp(OK,'Y') | strcmp(OK,'N')
            REPEAT=0;
        else
            REPEAT=1;
        end
        
    end
   
    
    if strcmp(OK,'Y')
        TXT = strcat('HOW MANY SPIKELETS?: ');
        Num = input(TXT);
        
        if Num==2
            clear ROW
            ROW = size(tempwave.WAVES.n2,1)+1;
            tempwave.WAVES.n2(ROW,:)=WAVEFORM;
        end
        
        if Num==3
            clear ROW
            ROW = size(tempwave.WAVES.n3,1)+1;
            tempwave.WAVES.n3(ROW,:)=WAVEFORM;
        end
        
        if Num==4
            clear ROW
            ROW = size(tempwave.WAVES.n4,1)+1;
            tempwave.WAVES.n4(ROW,:)=WAVEFORM;
        end
        
        if Num==5
            clear ROW
            ROW = size(tempwave.WAVES.n5,1)+1;
            tempwave.WAVES.n5(ROW,:)=WAVEFORM;
        end
        
        if Num==6
            clear ROW
            ROW = size(tempwave.WAVES.n6,1)+1;
            tempwave.WAVES.n6(ROW,:)=WAVEFORM;
        end
        
        if Num==7
            clear ROW
            ROW = size(tempwave.WAVES.n7,1)+1;
            tempwave.WAVES.n7(ROW,:)=WAVEFORM;
        end
        
        if Num==8
            clear ROW
            ROW = size(tempwave.WAVES.n2,1)+1;
            tempwave.WAVES.n8(ROW,:)=WAVEFORM;
        end
        
        disp('!!! WAVEFORM STORED SUCCESSFULLY !!!');
        
    end
    
    
    
    
    
    
    REPEAT=1;
    rcount=0;
    while REPEAT==1;
        
        rcount=rcount+1;
        if rcount>1
            disp(' PLEASE RE-ENTER A VALID RESPONSE ');
        end
        
        TXT = strcat('DO YOU WANT TO PROCEED? ENTER Y or N: ');
        CONT = upper(input(TXT,'s'));
        
        if strcmp(CONT,'Y') | strcmp(CONT,'N')
            REPEAT=0;
        else
            REPEAT=1;
        end
        
    end
   
    
    
    
    
    if strcmp(CONT,'Y')
        disp('>>> Proceeding >>>');
        END_FLAG=0;
    end
    
    if strcmp(CONT,'N')
        disp('!!! Terminating the process !!!');
        END_FLAG=1;
    end
    
    
    
end




% PLOTTING THE LIBRARY ---------------------------

clear F;
F = figure();


subplot(3,2,1)
hold on;
if ~isempty(tempwave.WAVES.n2)
    TEMPY = mat2gray(tempwave.WAVES.n2);
    plot(tempwave.WAVES.n2');
end
box off;
axis off;
title('n2');

subplot(3,2,2)
hold on;
if ~isempty(tempwave.WAVES.n3)
    TEMPY = mat2gray(tempwave.WAVES.n3);
    plot(tempwave.WAVES.n3');
end
box off;
axis off;
title('n3');

subplot(3,2,3)
hold on;
if ~isempty(tempwave.WAVES.n4)
    TEMPY = mat2gray(tempwave.WAVES.n4);
    plot(tempwave.WAVES.n4');
end
box off;
axis off;
title('n4');

subplot(3,2,4)
hold on;
if ~isempty(tempwave.WAVES.n5)
    TEMPY = mat2gray(tempwave.WAVES.n5);
    plot(tempwave.WAVES.n5');
end
box off;
axis off;
title('n5');

subplot(3,2,5)
hold on;
if ~isempty(tempwave.WAVES.n6)
    TEMPY = mat2gray(tempwave.WAVES.n6);
    plot(tempwave.WAVES.n6');
end
box off;
axis off;
title('n6');

subplot(3,2,6)
hold on;
if ~isempty(tempwave.WAVES.n7)
    TEMPY = mat2gray(tempwave.WAVES.n7);
    plot(tempwave.WAVES.n7');
end
box off;
axis off;
title('n7');

% subplot(3,2,7)
% hold on;
% if ~isempty(tempwave.WAVES.n8)
%     TEMPY = mat2gray(tempwave.WAVES.n8);
%     plot(tempwave.WAVES.n8');
% end
% box off;
% axis off;
% title('n8');

suptitle('CS database library');

cd(Lib_dir)
filename = 'CS_Library';
print(F, '-dpdf', filename, '-r400')


% SAVING THE LIBABRY ------------------------------

WAVES = tempwave.WAVES;
cd(Lib_dir)
save(FILENAME,'WAVES');




end