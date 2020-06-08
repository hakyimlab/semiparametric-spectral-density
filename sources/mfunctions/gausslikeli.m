function y = gausslikeli(data,cova,mu)
%------------------------------------
% gausslikeli(data,cova,mu)
% computes gaussian likelihood 
if nargin<3
    mu=0;
end
cc=chol(cova); invcc=inv(cc); %improve speed by using cc is upper triangular
%% in matlab inverse of sparse(cc) is not faster than inv of cc (tested for n=100, 1000)
xx=invcc'*(data-mu);
y = - sum(log(diag(cc))) - xx'*xx/2;
