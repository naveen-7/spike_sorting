%% helper function for getspikes_n
% written by naveen at JLG on 9/2/19

function IND = mergespkgrps_n(IND)

% IND    : list of indeces
try
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
catch
    disp('!!! ERROR DETECTED !!!');
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
end

end






