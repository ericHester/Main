function AU_process
% This script is written as a "function" so that the 2 subfunctions 
% below can be included in the same source file.

% **********************************************
% Choose host: ACRO=1, PORI=2, CCA=3, or TURF=4.
host=4;
% **********************************************

load MASTERMATRIX master
load HOST_ISLAND H
% A column of master corresponds to a sample.  A row of master gives
% the counts for an OTU in the various samples.  The row number 
% happens to be (OTU_ID + 1), since the OTUs were numbered starting
% with zero.
% The vector H{host} contains the column numbers of master that
% represent samples of the pure host.

% ***** OTU SCATTERGRAM ***************
% Prepare scattergram.
cols=H{host};
M=master(:,cols);
ss=sum(M);                  % sample sizes (row vector)
num=length(ss);             % number of samples
[whichrows,abun,ubiq]=AU(M);
numrows=length(whichrows);  % total number of OTUs seen on the host
% Plot scattergram.
semilogy(ubiq,abun,'r.')
hold on
% *************************************

% ***** OTU AU-CONFIDENCE LEVELS ******
% Calculate confidance levels for data. Confidence = 1-significance.
% See comments in STATS (in SUBFUNCTIONS below) for connections
% with quantities mentioned in the appendix.
for j=1:numrows
    ab=abun(j);              % relative abundance of the OTU in host pool.
    [m0,p0,E,V]=STATS(ab,ss);
    % Only m0 and p0 are used in this loop.
    k=num*ubiq(j);          % number of samples in which the OTU was seen.
    m0_minus_k=abs(m0-k);           % just n - k
    conf(j,1)=betainc(p0,k+1,m0_minus_k);
end
% ************************************

% ***** AU EXPECTATION CURVE *********
% Prepare expectation curve.
u=linspace(0,.995,501);
log_a=linspace(-5,-.5,501); % logarithmic scale
[U,L]=meshgrid(u,log_a);        
A=10.^L;                % U and A are used below in the contour command.
k=num*u;    % row vector
a=(10.^log_a)'; % column vector
[m0,p0,E,V]=STATS(a,ss); % column vectors
% Only E is used here.
% Plot expectation curve.
u=E/num;
hh=semilogy(u,a);
set(hh,'LineWidth',1)
xlabel('Ubiquity','fontsize',16,'fontweight','b')
ylabel('Relative Abundance ','fontsize',16,'fontweight','b')
hold on
% ************************************

% ***** AU SIGNIFICANCE CONTOUR ******
% Prepare significance function for contour graph.
[K,M0]=meshgrid(k,m0);
[K,P]=meshgrid(k,p0);
M0_minus_K=abs(M0-K);
Z=1-betainc(P,K+1,M0_minus_K);
% Plot confidence contour (significance curve).
Map=[0 0 0];
colormap(Map)
vec=[.01];      % just one contour
[c,h]=contour(U,A,Z,vec);
set(h,'LineWidth',1)
grid on
hold off
% ************************************




% *******************************************************
% ************ SUBFUNCTIONS *****************************
% *******************************************************

function [m0,p0,E,V]=STATS(a,ss)
% ss is a row vector of sample sizes
% a, m0, p0 are column vectors
[SS,A]=meshgrid(ss,a);
PI=1-(1-A).^SS;
E=sum(PI,2);	% expectation E[K] in appendix
s=sum(PI.^2,2);
V=E-s;          % variance Var[K] in appendix
p0=s./E;        % p* in appendix
m0=(E.^2)./s;	% m* in appendix

function  [rows,abun,ubiq]=AU(M)
% Pool the samples.
sM=sum(M,2);
% Save rows present in the pool.
rows=find(sM>0);        % rows where the OTU is seen on the host
% Contract the matrix.
m=M(rows,:);
sm=sum(m,2);
ssm=sum(sm);
abun=sm/ssm;           % pooled relative abundance
m0=m>0;                % logical
k=sum(m0,2);           % number of samples where OTU shows up
n=size(M,2);           % total number of samples
ubiq=k/n;              % ubiquity

