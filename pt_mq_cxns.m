function [X] = pt_mq_cxns(C,mu)
% C is cost of matrix of M rows and N columns
% N is number of partials in frame frame0
% M is number of partials in frame frame1
% c is vector of length N containing the row indices of least cost for each
% column in C. If the cost corresponding to any indices exceeds mu, these
% indices are set to 0. Let there be L such entries.
% X is a matrix containing K=min(M,N)-L entries set to 1, the rest 0.
% X(i,j) corresponds to kth smallest value of c where i=c(j), i != 0
% THIS IS NOT CORRECT YET
[M,N]=size(C);
K=min([M N]);
X=zeros(M,N);
for k=1:K
    [w,wi]=min(C,[],1);
    [w_,wj]=min(w,[],2);
    i_=wi(wj);
    j_=wj;
    if C(i_,j_) < mu
        X(i_,j_)=1;
    end
    C(i_,:)=Inf;
    C(:,j_)=Inf;
end
%    
%[mi,imi]=min(C,[],1);
%P=length(mi);
%% keep only K or P smallest, whatever is smaller
%[s_,i_]=sort(mi,'ascend');
%i_=[i_(1:min([K,P])),zeros(1,P-min([K,P]))];
%imi.*=sum((ones(P,1)*(1:P))==(i_'*ones(1,P)),1);
%imi=imi.*(mi<mu);
%X=(ones(M,1)*imi)==((1:M)'*ones(1,N));
