function [result,N,time] = TNBCr(foldermodelname,basemodelname,variantmodelname)

% filename_MLsys = ['MLsys/' foldermodelname '/' variantmodelname '_MLsys.mat'];
filename_modeldelta = ['MLsys/' foldermodelname '/' basemodelname '_to_' variantmodelname '_modeldelta.mat'];
filename_NBCresult = ['MLsys/' foldermodelname '/' basemodelname '_NBCr_result.mat'];

% load(filename_MLsys)
load(filename_modeldelta)
load(filename_NBCresult)

if size(N,2)<size(Xprime,2)
    Nprime = [N zeros(1,size(Xprime,2)-size(N,2))];
else
    Nprime = N;
end

tic

Xc = zeros(size(Xprime));
for ti = 1:size(transmin,2)
    td=transmin(:,ti);
    xor=td(1);
%     sigma=td(2);
    xtar=td(3);
    if Nprime(xor)==1 && Nprime(xtar)==1
        Xc(xor)=1;
    end
end

Xa = (Xmmin+Xc)>0;

if sum(Xa)>0
    [Noutput,breakflag] = BRSbreak(Nprime,Sigma_cprime+Sigma_uprime,transprime,Xmprime,Xa);
    %If BRS breaks, breakflag=true, Nprime=Nprime will hold. Otherwise set
    %Nprime to the output of BRS; Noutput
    if ~breakflag
        Nprime=Noutput;
    end
end

check = false;
XsetminN = Xprime-Nprime;

Xd=((Xmplus+XsetminN)==2);
if sum(Xd)>0
    check = true;
else
    for ti = 1:size(transplus,2)
        td=transplus(:,ti);
        xor=td(1);
    %     sigma=td(2);
        xtar=td(3);
        if XsetminN(xor)==1 && Nprime(xtar)==1
            check = true;
            break
        end
    end
end

if check==true
    Nprime = BRS(Xprime,Sigma_cprime+Sigma_uprime,transprime,(Nprime+Xmprime)>0);
end

time=toc;

if isequal(Nprime,Xprime)
    result = true;
else
    result = false;
end
N=Nprime;
save(['MLsys/' foldermodelname '/' variantmodelname '_TNBCr_result.mat'],'result','N')

end
