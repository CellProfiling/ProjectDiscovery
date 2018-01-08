function ptot = calcOdds(nvotes, nplayers,nclasses,maxchoose)
%This function calculates the odds that a given outcome is seen at least
%that many times at random using a hypergeometric distribution
%
%inputs:
%nvotes - number of votes for a given class
%nplayers - number of total people voting
%nclasses - how many classes there are to choose from
%maxchoose - how many classes it is possible to pick at a time
%
%Updated 6,April,2017 to use the ramanujan approximation for large numbers. This should
%allow us to compute larger numbers with better accuracy.
%https://en.wikipedia.org/wiki/Factorial#Rate_of_growth_and_approximations_for_large_n

% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published
% by the Free Software Foundation; either version 2 of the License,
% or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
% 02110-1301, USA.
%
% For additional information visit proteinatlas.org or
% send email to devin.sullivan@scilifelab.se

if nargin<4
    maxchoose = 5;
end
if nargin<3
    nclasses = 29;
end

ptot = 0;
for n = nvotes:nplayers
%     numerator = (nchoosek(nplayers,n)*nchoosek(nclasses*nplayers-nplayers,maxchoose*nplayers-n));
%    numerator = (bignchoosek(nplayers,n)*bignchoosek(nclasses*nplayers-nplayers,maxchoose*nplayers-n));
%     nclassxnp_p = nclasses*nplayers-nplayers;
%     mch_np = maxchoose*nplayers-n;
%     numerator = exp( gammaln(nplayers+1)-gammaln(n+1)-gammaln(nplayers-n+1) )*...
%         exp( gammaln(nclassxnp_p+1)-gammaln(mch_np+1)-gammaln(nclassxnp_p-mch_np+1) );
%     denominator = nchoosek(nclasses*nplayers,maxchoose*nplayers);
%    denominator = bignchoosek(nclasses*nplayers,maxchoose*nplayers);
%     nclassxnp = nclasses*nplayers;
%     denominator = exp( gammaln(nclassxnp+1)-gammaln(mch_np+1)-gammaln(nclassxnp-mch_np+1) );
%    p = numerator/denominator;
    p = hygepdf(n,nclasses*nplayers,nplayers,maxchoose*nplayers);
    ptot = ptot+p;
end

end



function ram = ramanujan(n)
  ram = n*log(n) - n + log(n*(1 + 4*n*(1+2*n)))/6 + log(pi)/2;
end
% 
% nchoosek <- function(n,k){
%   factorial(n)/(factorial(k)*factorial(n-k))
% } 

function nchk = bignchoosek(n,k)
  if n==k
      nchk=1;
  else
      nchk = exp(ramanujan(n) - ramanujan(k) - ramanujan(n-k));
  end
end

% 
% numerator = 1;
% for n = nvotes:votecount
%     %for i = 1:n
%         %Hypergeometric distribution
%         %numerator = nchoosek(maxchoose,n)*
%         %try
%         numerator = nchoosek(votecount,n)*numerator;
%       %  catch 
%        %     problem = 1
%         %end
%     %end
% 
% end
% 
% numerator = numerator*(nclasses*votecount-votecount*maxchoose)
% 
% denom = nchoosek(nclasses*votecount,maxchoose*votecount);
% 
% p = numerator/denom
% 
% 
