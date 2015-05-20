function dnac = seqcomplement(dna)
%SEQCOMPLEMENT returns the complementary strand of a DNA sequence.
%
%   SEQ = SEQCOMPLEMENT(DNA) calculates the complementary strand (A-->T,
%   C-->G, G-->C, T-->A) of sequence DNA.
%
%   SEQ is returned in the same format as DNA, so if DNA is an integer
%   sequence then so is SEQ.
%
%   Example:
%
%       % Create a random sequence of nucleotides and the complement.
%       seq = randseq(25)
%       compseq = seqcomplement(seq)
%
%   See also SEQRCOMPLEMENT, SEQREVERSE, SEQVIEWER.

%   Copyright 2002-2012 The MathWorks, Inc.


% Full mapping of all nucleotide symbols
% * A C G T(U) R Y K M S W B D H V N -
% * T G C A    Y R M K S W V H D B N -

% Indices used
% * A C G T(U) R Y K M S  W  B  D  H  V  N
% 0 1 2 3 4(4) 5 6 7 8 9 10 11 12 13 14 15

% If the input is a structure then extract the Sequence data.
if isstruct(dna)
    dna = bioinfoprivate.seqfromstruct(dna);
end

origdna = dna;
form = class(dna);
if ischar(dna)
    lowermap = (lower(dna) == dna);
    dna = mynt2int(dna);
end
map =  mynt2int('*TGCAYRMKSWVHDBN-');
dnac = dna;
for count = 1:length(dna)
    dnac(count) = uint8(map(double(dna(count))+1));
end

if strcmp(form,'char')  %#ok (dna now converted to int)
    dnac = myint2nt(dnac);
    % RNA case
    if ~isempty(regexpi(origdna,'u','once'))
        dnac = mydna2rna(dnac);
    end

    dnac(lowermap) = lower(dnac(lowermap));
end





