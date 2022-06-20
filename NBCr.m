function [result,N,time] = NBCr(foldermodelname,modelname)
load(['MLsys/' foldermodelname '/' modelname '_MLsys.mat'])

tic

N = BRS(X,Sigma_c+Sigma_u,transX,Xm);
time = toc;

if isequal(N,X)
    result = true;
else
    result = false;
end

save(['MLsys/' foldermodelname '/' modelname '_NBCr_result.mat'],'result','N')

end