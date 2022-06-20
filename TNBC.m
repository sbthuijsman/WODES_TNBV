function [result,N,time] = TNBC(foldermodelname,basemodelname,variantmodelname)

filename_modeldelta = ['MLsys/' foldermodelname '/' basemodelname '_to_' variantmodelname '_modeldelta.mat'];
filename_NBCresult = ['MLsys/' foldermodelname '/' basemodelname '_NBC_result.mat'];

load(filename_modeldelta)
load(filename_NBCresult)

if size(N,2)<size(Xprime,2)
    Nprime = [N zeros(1,size(Xprime,2)-size(N,2))];
    Qprime = [Q zeros(1,size(Xprime,2)-size(N,2))];
else
    Nprime = N;
    Qprime = Q;
end

tic

Xbx = zeros(size(Xprime));
for ti = 1:size(transmin,2)
    td=transmin(:,ti);
    xor=td(1);
%     sigma=td(2);
    xtar=td(3);
    if Qprime(xor)==1 && Qprime(xtar)==1
        Xbx(xor)=1;
    end
end

Xb = (X0min+Xbx)>0;

if sum(Xb)>0
    [Qoutput,breakflag] = FRSbreak(Qprime,Sigma_cprime+Sigma_uprime,transprime,X0prime,Xb);
    %If BRS breaks, breakflag=true, Nprime=Nprime will hold. Otherwise set
    %Nprime to the output of BRS; Noutput
    if ~breakflag
        Qprime=Qoutput;
    end
end

check = false;
XsetminQ = Xprime-Qprime;

Xd=((X0plus+XsetminQ)==2);
if sum(Xd)>0
    check = true;
else
    for ti = 1:size(transplus,2)
        td=transplus(:,ti);
        xor=td(1);
    %     sigma=td(2);
        xtar=td(3);
        if Qprime(xor)==1 && XsetminQ(xtar)==1
            check = true;
            break
        end
    end
end

if check==true
    Qprime = FRS(Xprime,Sigma_cprime+Sigma_uprime,transprime,(Qprime+X0prime)>0);
end

% [Xprime,transXprime,~,Xmprime] = restrict(Xprime,Sigma_cprime,Sigma_uprime,transXprime,X0prime,Xmprime,Qprime);

%restrict
Xprime=Xprime.*Qprime;
Xmmin=Xmmin.*Qprime;
Xmplus=Xmplus.*Qprime;

Nprime = Nprime.*Qprime;

transprime = quick_prune_trans(transprime,Xprime,Sigma_cprime+Sigma_uprime,zeros(size(Xprime)));
if ~isempty(transplus)
    transplus = quick_prune_trans(transplus,Xprime,Sigma_cprime+Sigma_uprime,zeros(size(Xprime)));
end
if ~isempty(transmin)
    transmin = quick_prune_trans(transmin,Xprime,Sigma_cprime+Sigma_uprime,zeros(size(Xprime)));
end

%X0 and Sigma do not matter in remaining calculation


%% TNBCr below

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
Q=Qprime;
save(['MLsys/' foldermodelname '/' variantmodelname '_TNBC_result.mat'],'result','N','Q')
end
