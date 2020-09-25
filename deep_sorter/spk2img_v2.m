%% Function to convert raw wave times to images
% written by naveen at JLG on 5/22/20

function IMG = spk2img_v2(signal, img_h, img_w, offset)

FLAG = -1;

Freq = 50; % 50k
% signal = RAW_S; % signal must be 50k Hz
H_IND = -1.5:0.01:1.5;
img_w = 50*Freq; % 50 ms
offset = 0;
count = 0;

mkdir(strcat(pwd,'\spkIMG'))
cd(strcat(pwd,'\spkIMG'));

for i = offset+1:img_w:length(signal)-img_w
    try
        count = count+1;
        clear SIG;
        SIG = signal(i:i+img_w-1);
        
        
        clear F I J;
        fname = strcat('TEMP',num2str(count),'.png');
        F = figure(); plot(SIG); xlim([1 img_w]); ylim([H_IND(1) H_IND(end)]); axis off; saveas(F,fname);
        I = imread(fname);
        I = I(:,:,1);
        J = imresize(I, 0.3);
        %     imagesc(J)
        close(F);
%         delete 'TEMP.png'; 
        
        IMG{i,1} = J;
        
        PERCENT = (count/length(offset+1:img_w:length(signal)))*100;
        if PERCENT>10 & FLAG==-1 disp(strcat(num2str(PERCENT),'% done')); FLAG = 0; end
        if PERCENT>25 & FLAG==0 disp(strcat(num2str(PERCENT),'% done')); FLAG = 1; end
        if PERCENT>50 & FLAG==1 disp(strcat(num2str(PERCENT),'% done')); FLAG = 2; end
        if PERCENT>75 & FLAG==2 disp(strcat(num2str(PERCENT),'% done')); FLAG = 3; end
        if PERCENT>95 & FLAG==3 disp(strcat(num2str(PERCENT),'% done')); FLAG = 4; end
    catch
    end
end

















end