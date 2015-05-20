function dnar = seqreverse(dna)
%SEQREVERSE returns the reverse strand of a DNA or RNA sequence.
%
%   SEQ = SEQREVERSE(DNA) calculates the reverse strand 3' --> 5' of sequence DNA.
%
%   SEQ is returned in the same format as DNA, so if DNA is an integer sequence then
%   so is SEQ.
%
%   Example:
%
%       % Create a random sequence of nucleotides and reverse it.
%       seq = randseq(25)
%       revseq = seqreverse(seq)
%
%   See also FLIPLR, SEQCOMPLEMENT, SEQRCOMPLEMENT, SEQVIEWER.

%   Copyright 2002-2012 The MathWorks, Inc.


% If the input is a structure then extract the Sequence data.
if isstruct(dna)
    dna = bioinfoprivate.seqfromstruct(dna);
end

dnar = fliplr(dna);




