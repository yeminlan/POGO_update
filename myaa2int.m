function  [seq , map] = aa2int(aa, varargin)
%AA2INT converts a string of amino acids from letters to numbers.
%
%   SEQ = AA2INT(AA) converts string AA of amino acids into an array of
%   integers using the following mapping table:
%
%   A R N D C Q E G H I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *  -  ?
%   1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 0
%
%   Where B is D or N (aspartic), Z is E or Q (glutamic), X represents any
%   amino acid, * represents an end terminator, - is a gap, and ? is an
%   unknown amino acid.
%
%   SEQ = AA2INT(AA,...,'UNKNOWN',I) defines the number used to represent
%   unknown amino acid. The default value is 0.
%
%   Example:
%
%   s = aa2int('MATLAB')
%
%   See also AMINOLOOKUP, INT2AA, INT2NT, NT2INT.

%   Copyright 2002-2012 The MathWorks, Inc.

if isempty(aa)
    seq = uint8([]);
    return
end

% If the input is a structure then extract the Sequence data.
if isstruct(aa)
    aa = bioinfoprivate.seqfromstruct(aa);
end

unknown = 0;
origsize = size(aa);
if  nargin > 1
    if rem(nargin,2)== 0
        error(message('bioinfo:aa2int:IncorrectNumberOfArguments', mfilename));
    end
    okargs = {'unknown'};
    for j=1:2:nargin-2
        pname = varargin{j};
        pval = varargin{j+1};
        k = find(strncmpi(pname, okargs,numel(pname)));
        if isempty(k)
            error(message('bioinfo:aa2int:UnknownParameterName', pname));
        elseif length(k)>1
            error(message('bioinfo:aa2int:AmbiguousParameterName', pname));
        else
            unknown = pval;
        end
    end
end

%      A R N D C Q E G H I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  * -  ?
%
%      1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26

map = uint8(...
    [1 21 5 4 7 14 8 9 10 0 12 11 13 3 0 15 6 2 16 17 0 20 18 23 19 22 24 25 0]);
%   'a b  c d e f  g h i  j k  l  m  n o p  q r s  t  u v  w  x  y  z  *  -  ?'

if unknown ~= 0
    map(map == 0) = unknown;
end

aa = lower(aa)-96;  % 96 = 'a'-1

% find gap ('-') and stop ('*') symbols
gaps  = find(aa==-51);
stops = find(aa==-54);

% adjusting values, so we can use the 'map'
aa = min(max(aa,0),27);  %  0 >= aa >= 27
aa(aa==0)  = 29;
aa(aa==27) = 29;
aa(gaps)   = 28;  %#ok -- note that aa has changed
aa(stops)  = 27;  %#ok

seq = map(aa);

seq = reshape(seq,origsize);
