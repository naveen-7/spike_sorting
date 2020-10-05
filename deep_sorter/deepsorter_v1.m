%% deepsorter_1.0




close all;
clear all;
clc;
% Setup directories--------------------------------------------------------
%%%%% enter your directory here
codes_dir = 'E:\NAVEEN_Work\Cerebellum\Codes\CER_codes_NEW\LINEAR';
data_dir= 'E:\NAVEEN_Work\Cerebellum\Data\RECORDED_CELLS';

disp('!!!!!  getspikes_n has started running  !!!!!')
cd(data_dir);

Fs = 50000/1000;  %%%% frequency of signal MUST BE 50000 Hz


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% LOADING THE FILE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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

 







