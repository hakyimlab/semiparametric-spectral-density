function res=gaussreml(data,cova,M)
% gaussreml(data,cova,M){
% computes gaussian restrited likelihood
% cova=covariance matrix; M = model matrix = [f(x_1) ...f(x_n)]'
    % some reshaping on input
    n = length(data);data=reshape(data,n,1);
    % 
    cc=chol(cova);invcc=inv(cc);
    invcova = invcc * invcc';
    W = M' * invcova * M;
    ww = chol(W); invww=inv(ww);
    invW = invww * invww';
    %                                                                            
    res = - sum(log(diag(cc))) - sum(log(diag(ww))) - ...
      .5 * data' * invcova * ...
    (eye(n) - M * invW * M' * invcova) * data;



