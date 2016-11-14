function sW=strengthcorrect(W,ss,nreps)
% sW=strengthcorrect(W,ss,[nreps])
%
% Iteratively rescale weights in connectivity matrix W to make the strength
% sequence converge onto desired strength sequence ss
%
% Inputs
% W         = weighted connectivity matrix
% ss        = desired strength sequence
% nreps     = number of iterations in the rescaling algorithm 
%             (optional, default=9)
%
% Outputs
% sW        = network with desired strength sequence ss
%
% Note: assumes W is undirected (symmetric)
%
% Reference: Roberts et al. (2016) NeuroImage 124:379-393.
%
% M Breakspear, J Roberts
% QIMR Berghofer, 2015

if nargin<3,
    nreps=9;
end

ss=ss(:)';
N=length(ss);
discnodes=sum(W)==0;
sW=W.*repmat(ss./sum(W),N ,1);
sW(:,discnodes)=0; %ensure disconnected nodes don't introduce NaNs
sW=(sW+sW')/2;
for j=1:nreps,         
    sW=sW.*repmat(ss./sum(sW),N ,1);
    sW(:,discnodes)=0;
    sW=(sW+sW')/2;
end