
cd('E:\NAVEEN_Work\Cerebellum\Codes\CER_codes_NEW\SPIKE_SORTING_n\Libraries');
load('ComplexSpikes_n.mat')

TEMP  = WAVES;

clear MAT IND
RANDS = randi(10,6,1);

MAT = NaN(1,size(TEMP.n2,2)); IND = NaN(1,1);
for i=2:6
    MAT = [MAT; repmat(eval( strcat('TEMP.n',num2str(i)) ) , RANDS(i),1)];
    IND = [IND; repmat(i, RANDS(i)* size(eval( strcat('TEMP.n',num2str(i))),1),1)];
end
MAT = MAT(2:end,:);
IND = IND(2:end,:);

MAT = MAT+randn(size(MAT))/100;
for i=1:size(MAT,1)
    MAT(i,:) = mat2gray(MAT(i,:));
end



% [e_vec, e_val] = pca_n(MAT);
% score = e_vec';
% F = figure
% scatter3(score(:,1),score(:,2),score(:,3),'.')
% axis equal
% xlabel('PC1')
% ylabel('PC2')
% zlabel('PC3')



NN = 4;
indd = kmeans(MAT,NN);
figure
for i=1:NN
   subplot(2,2,i)
   hold on;
   plot(MAT(find(indd==i),:)','color',[0.7 0.7 0.7]);
   plot(nanmean(MAT(find(indd==i),:))','-r');
end


