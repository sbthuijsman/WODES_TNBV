function trans = quick_prune_trans(trans,present_states,present_sigma,prune_tar_states)
present_states=[present_states 0];
present_sigma=[present_sigma 0];
transxor = ismember(trans(1,:),find(present_states));
transsigma = ismember(trans(2,:),find(present_sigma));
transxtar = ismember(trans(3,:),find(present_states));

removetranstofrom = ((transxor+transsigma+transxtar)~=3);

if sum(prune_tar_states)>0
    transtofound = ismember(trans(3,:),find(prune_tar_states));
    
    removetrans = (removetranstofrom+transtofound~=0);
    trans(:,removetrans)=[];
else
    trans(:,removetranstofrom)=[];
end

end