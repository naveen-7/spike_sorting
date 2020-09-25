%% Function to convert raw wave times to images
% written by naveen at JLG on 5/8/20

function IMG = spk2img(signal, img_h, img_w, offset)

FLAG = -1;

Freq = 50; % 50k
% signal = RAW_S; % signal must be 50k Hz
H_IND = -1.5:0.01:1.5;
img_h = length(H_IND);
img_w = 50*Freq; % 50 ms
offset = 0;
count = 0;

for i = offset+1:img_w:length(signal)-img_w
    count = count+1;
    clear SIG; SIG = signal(i:i+img_w-1);
    clear BLANK; BLANK = zeros(img_h, img_w);
    
    for j=1:length(SIG)
        Y_IND(j,1) = find(H_IND>=round(SIG(j),2),1);
    end
    
    y = Y_IND;
    for ii = 1:length(SIG)
        BLANK(y(ii),ii) = 1;
    end
    
    % %     imagesc(IMG{count-1,1});
    IMG{count,1} = logical(BLANK);
    
    PERCENT = (count/length(offset+1:img_w:length(signal)))*100;
    if PERCENT>10 & FLAG==-1 disp(strcat(num2str(PERCENT),'% done')); FLAG = 0; end
    if PERCENT>25 & FLAG==0 disp(strcat(num2str(PERCENT),'% done')); FLAG = 1; end
    if PERCENT>50 & FLAG==1 disp(strcat(num2str(PERCENT),'% done')); FLAG = 2; end
    if PERCENT>75 & FLAG==2 disp(strcat(num2str(PERCENT),'% done')); FLAG = 3; end
    if PERCENT>95 & FLAG==3 disp(strcat(num2str(PERCENT),'% done')); FLAG = 4; end
    
end

















end