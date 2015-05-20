function  [seq, map] = nt2int(nt,varargin)
%NT2INT convert a string of nucleotides from letters to numbers.
%
%   SEQ = NT2INT(NT) converts string NT of nucleotides into an array of
%   integers. A => 1 , C => 2, G => 3,  T(U) => 4. Unknown characters are
%   mapped to 0.
%
%   SEQ = NT2INT(NT,...,'UNKNOWN',I) defines the number used to represent
%   unknown nucleotides. The default value is 0.
%
%   SEQ = NT2INT(NT,...,'ACGTOnly',TF) forces ambiguous nucleotide characters
%   (R, Y, K, M, S, W, B, D, H, V and N) to be represented by the
%   unknown nucleotide value.
%
%   Example:
%
%      s = nt2int('ACTGCTAGC')
%
%   See also AA2INT, BASELOOKUP, INT2AA, INT2NT.

%   Copyright 2002-2012 The MathWorks, Inc.



% If the input is a structure then extract the Sequence data.
if isstruct(nt)
    nt = bioinfoprivate.seqfromstruct(nt);
end

if isempty(nt)
    seq = uint8([]);
    return
end

unknown = 0;
acgtonly = false;
origsize = size(nt);
nt = nt(:);
nt = nt';

if  nargin > 1
    if rem(nargin,2)== 0
        error(message('bioinfo:nt2int:IncorrectNumberOfArguments', mfilename));
    end
    okargs = {'unknown','acgtonly'};
    for j=1:2:nargin-2
        pname = varargin{j};
        pval = varargin{j+1};
        k = find(strncmpi(pname, okargs,numel(pname)));
        if isempty(k)
            error(message('bioinfo:nt2int:UnknownParameterName', pname));
        elseif length(k)>1
            error(message('bioinfo:nt2int:AmbiguousParameterName', pname));
        else
            switch(k)
                case 1   % unknown
                    unknown = pval;
                case 2  % others forces everything except ACTG to be the unknown value
                    acgtonly = bioinfoprivate.opttf(pval);
                    if isempty(acgtonly)
                        error(message('bioinfo:nt2int:InputOptionNotLogical', upper( char( okargs( k ) ) )));
                    end
            end
        end
    end
end

% * A C G T U R Y K M S  W  B  D  H  V N
% 0 1 2 3 4 4 5 6 7 8 9 10 11 12 13 14 15

%           'a  b c  d e f g  h i j k l m  n o p q r s t u  v  w  x y z  -'
map = uint8([1 11 2 12 0 0 3 13 0 0 7 0 8 15 0 0 0 5 9 4 4 14 10 15 6 0 16]);
if acgtonly
    map(map > 4) = 0;
end
if unknown ~= 0
    map(map == 0) = unknown;
end
nt = lower(nt);
nt = strrep(nt,'{','z');
nt = strrep(nt,'-','{');
nt = strrep(nt,'*','z');
nt = strrep(nt,'?','z');
seq = map((uint8(nt) + 1) - uint8('a'));
seq = reshape(seq,origsize);
