function [Wwp,Wsp,Wssp]=geombinsurr(W,D,nbins,bintype)
% [Wwp,Wsp,Wssp]=geombinsurr(W,D,nbins,[bintype])
%
% Surrogate random graphs that approximately preserve the relationship
% between weight and distance by first binning the edges by distance and
% then shuffling within the bins
%
% Inputs
% W         = weighted connectivity matrix
% D         = physical distance between nodes
% nbins     = number of bins
% bintype   = 'quantiles' for bins containing equal numbers of edges
%             or 'equalwidth' for bins of equal width
%             (optional, default='quantiles')
%
% Outputs
% Wwp       = random surrogate preserving weights but not strengths
% Wsp       = preserves strength distribution but not sequence
% Wssp      = preserves strength sequence
%
% Note: Wsp and Wssp are generated assuming that W is undirected.
%
% Reference: Roberts et al. (2016) NeuroImage 124:379-393.
%
% J Roberts, M Breakspear
% QIMR Berghofer, 2015

if nargin<4
    bintype='quantiles';
end

%Preliminaries
N = size(W,1);
Wwp = zeros(N,N); 

%Check if directed
drct=1;
if max(max(W-W'))==0, drct=0; end

W(1:N+1:end)=0;         %Ensure there's no self connections
if drct==0,             %Only do one triangle if undirected, then replicate
    W = triu(W);
end

nz= find(W);            %Nonzero entries
w = W(nz);              %Vector of non-zero connections
d = D(nz);              %Corresponding distances

%Set up bins by distance
switch bintype
    case 'quantiles'
        %bin by quantiles (bins contain equal numbers of edges)
        if nbins>2
            binedges=[min(d) quantile(d,nbins-1) max(d)+eps]; 
        elseif nbins==2
            binedges=[min(d) median(d) max(d)+eps];
        else
            binedges=[min(d) max(d)+eps];
        end
    case 'equalwidth'
        %bin into fixed-width bins 
        binedges=linspace(min(d),max(d)+eps,nbins+1);
    otherwise
        error('invalid bintype: %s (must be ''quantiles'' or ''equalwidth'')',bintype)
end

% bin the edges
[counts,inds]=histc(d,binedges);
for j=1:nbins
    % find which edges are in a given bin
    whichw=find(inds==j);  
    if counts(j)>0
        % shuffle all edges within the bin
        w(whichw)=w(whichw(randperm(counts(j))));  
    end
end

Wwp(nz)=w;

if drct==0;
    Wwp=(Wwp+Wwp');
    W=W+W';
end

%Adjust node strengths - NOTE: assumes W is undirected
strengthsW=sum(W);     % original node strengths
strengthsWwp=sum(Wwp); % new node strengths
%Re-order the old strengths to match new random order in Wwp
strengthsWsp=rankreorder(strengthsW,strengthsWwp);
%Adjust strengths to give both original and random strength sequences
Wsp=strengthcorrect(Wwp,strengthsWsp); % orig strengths in new sequence
Wssp=strengthcorrect(Wwp,strengthsW);  % orig strengths in orig sequence

end

function out=rankreorder(x,scaffold)
% reorder vector x to have same rank order as vector scaffold
y(:,1)=scaffold;
y(:,2)=1:length(x);
y=sortrows(y,1);
y(:,1)=sort(x);
y=sortrows(y,2);
out=y(:,1);
end