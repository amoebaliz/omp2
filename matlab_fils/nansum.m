function y = nanmean(x)
% NANSUM     Sum of matrix columns, ignoring NaNs
%===================================================================
% NANSUM   1.2   92/04/14  Copyright (C) Phil Morgan 1991
%
% function y = nansum(x)
%
% DESCRIPTION:
%    Sum of matrix columns, ignoring NaNs
% 
% INPUT:
%    x    = vector or matrix 
%
% OUTPUT:
%    y    = column-wise sum of x.  If x a vector then y = sum(x)
%           ignoring all NaNs.  Thus a sum of actual data values.
%
% EXAMPLE:  A = [ 1  2  3;
%                 3 NaN 5];
%           y = sum(x)
%           y = [4 2 8]
%
% CALLER:   general purpose
% CALLEE:   none
%
% AUTHOR:   Phil Morgan 3-09-91
%==================================================================

% @(#)nansum.m   1.2   92/04/14
% 
%--------------------------------------------------------------------

[m,ncols]=size(x);

% IF A ROW VECTOR THEN TRANSPOSE TO COLUMN VECTOR
if m == 1
  x = x';
  ncols = 1;
end

% FOR EACH COLUMN FIND SUM EXCLUDING NaNs
for icol = 1:ncols
   good    = find( ~isnan(x( :,icol)) );
   if length(good)>0
      y(icol) = sum( x(good,icol) );
   else
      y(icol) = NaN;
   end
end


   

