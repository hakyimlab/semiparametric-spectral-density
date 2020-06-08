function [rhon gg] = sieghankel(h_fun,R,varargin)
%% [rhon gg] = sieghankel(@fun,R,varargin)
%%
%% calculates hankel tr of fun
%% g(rho) = int_0^inf r fun(r) J_nu(r rho) dr
%% siegman uses g(rho) = 2pi  int_0^inf r fun(r) J_nu(2pi r rho) dr
%%
%% h_fun is handle of fun
%%
%% R is range of the function. The sampling is exponential on
%% (sqrt(rrho),sqrt(bbeta))*R on spatial domain
%% (sqrt(rrho),sqrt(bbeta))/R on spectral domain
%%
%% varargin has the additional arguments to call h_fun
%%
%% rhon and gg(rhon) is returned
  
%% parameters
  global K1 K2 N
  K1 = 2; %% r0<=1/(K1 b)
  K2 = 2; %% dr_N ~= a b
  N = 2^12; %% number of sample points
  %% b = r0 exp(a N); beta = rho0 exp(a N)
  
  % calc a, bbeta, rrho
  a = fzero(@afun,log(K1/K2));
  bbeta = fzero(@bbetafun,N/K2);
  rrho = a*K2/K1^2;
  
  b = sqrt(bbeta);
  beta = b;
  r0 = sqrt(rrho);
  
  rn = r0*R * exp(a*(0:(N-1)));
  fn = [rn .* h_fun(rn,varargin{:})  zeros(1,N)];
  jm = a*rrho*exp(a*(0:(2*N-1))).*besselj(0,rrho*exp(a*(0:(2*N-1))));

  gm = real(fft(fft(fn).*ifft(jm)));

  rhon = rn/R^2;
  gg = gm(1:N)./rhon;

%% functions
  
function res = afun(x)
  global N K1 K2
  res = x*exp(x*N)-K1/K2;
  
function res = bbetafun(x)
  global N K1 K2
  res = K2*x*log(K1*x) - N;
  
  
