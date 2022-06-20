function [disc,breakflag] = FRSbreak(X,Sigma,trans,Xinit,Xbreak)
global phi;

breakflag=false;
nX=length(X);
disc = Xinit.*X; %discovered
curlyr = Xinit.*X; %current BFS layer

trans=quick_prune_trans(trans,X,Sigma,Xinit);

while sum(curlyr)~=0
nxtlyr = zeros(1,nX); %next BFS layer
    for u = find(curlyr)
        for edgi = find(trans(1,:)==u)
            phi = phi+1;
            v = trans(3,edgi);
            if disc(v)==0
                disc(v)=1;
                nxtlyr(v)=1;
                if sum(Xbreak-(Xbreak.*disc))==0
                    breakflag=true;
                    return
                end
            end
        end
    end
    curlyr=nxtlyr;
end

end
