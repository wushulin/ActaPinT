function [Rr,Rb,mtr,mtb]=RedBlackSubdomains(n,m,nx);
% REDBLACKSUBDOMAINS compute a red-black space-time decomposition
% [Rr,Rb,mtr,mtb]=RedBlackSubdomains(n,m,nx); computes for a 1D
%   space-time all at once matrix in time stepping order with n
%   spatial steps and m time steps representing a second order wave
%   equation for nx red spatial subdomains a red-black subdomain
%   arrangement for the Unmapped Tent Pitching algorithm. The output
%   are mtr red restriction matrices Rr and mtb black restriction
%   matrices Rb in time.

idr=round(n/nx*(0:nx));        % red interface locations
for i=1:nx                
  sxr{i}=idr(i)+1:idr(i+1);    % red subdomain space indices
end
idb=round(n/nx*(0.5:nx-0.5));  % black interface locations
for i=1:nx-1
  sxb{i}=idb(i)+1:idb(i+1) ;   % black subdomain space indices
end
mx=min([diff(idr) diff(idb)]); % smallest subdomain width
% mst=min(floor(mx/2)-1,m);      % maximum tent height, discrete CFL
mst=min(floor(mx/2),m);      % seems ok !?
str{1}=1:mst;                  % first red subdomain half size in time
mtr=1;id=mst;                  % id current end of red subdomain
while id<m
  mtr=mtr+1;
  idn=min(id+2*mst,m);         % end of next red subdomain
  str{mtr}=id+1:idn;           % red subdomain time indices
  id=idn;                      % update current end 
end
mtb=1;id=min(2*mst,m);         % first full size if possible
stb{1}=1:id;               
while id<m
  mtb=mtb+1;
  idn=min(id+2*mst,m);         % end of next black subdomain
  stb{mtb}=id+1:idn;           % black subdomain time indices
  id=idn;                      % update current end 
end
Id=speye(n*m);                 % to extract the R matrices
G=reshape(1:m*n,n,m);          % form unknown ennumeration
for i=1:nx                     % extract the red subdomains Rr
  for j=1:mtr
    id=G(sxr{i},str{j}); Rr{i,j}=Id(id(:),:);
  end
end
for i=1:nx-1                   % extract the black subdomains Rb
  for j=1:mtb
    id=G(sxb{i},stb{j}); Rb{i,j}=Id(id(:),:);
  end
end