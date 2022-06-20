function [result,N,time] = NBC(foldermodelname,modelname)

load(['MLsys/' foldermodelname '/' modelname '_MLsys.mat'])

tic

Q = FRS(X,Sigma_c+Sigma_u,transX,X0);
[X,transX,~,Xm] = restrict(X,Sigma_c,Sigma_u,transX,X0,Xm,Q);


N = BRS(X,Sigma_c+Sigma_u,transX,Xm);
time = toc;

if isequal(N,X)
    result = true;
else
    result = false;
end

save(['MLsys/' foldermodelname '/' modelname '_NBC_result.mat'],'result','N','Q')

end