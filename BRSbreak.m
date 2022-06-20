function [Xplus,breakflag] = BRSbreak(X,Sigma,trans,Xinit,Xbreak)
%reverse trans for BRS
trans = [trans(3,:)
         trans(2,:)
         trans(1,:)];

%this is same as FRS
[Xplus,breakflag] = FRSbreak(X,Sigma,trans,Xinit,Xbreak);

end