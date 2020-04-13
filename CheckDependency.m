[fList,pList] = matlab.codetools.requiredFilesAndProducts('Ua.m');

%%
IsUsed="DesiredEleSizes.m";

if ~any(contains(fList,"\"+IsUsed))
    fprintf("%s is not used by Ua \n",IsUsed)
end


%%

mFiles=dir('*.m');

for I=1:length(mFiles)
    
    if ~any(contains(fList,"\"+mFiles(I).name))
        fprintf("%s is not used by Ua \n",mFiles(I).name)
    end
end
