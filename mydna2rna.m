function rna = dna2rna(dna)
%DNA2RNA converts a DNA sequence into an RNA sequence.
%
%   RNA = DNA2RNA(DNA) converts any thymine nucleotides in a DNA sequence
%   into uracil (T-->U).
%
%   RNA is returned in the same format as DNA, so if DNA is an integer
%   sequence then so is RNA.
%
%   Example:
%
%       dna2rna('ACGATGAGTCATGCTT')
%
%   See also REGEXP, RNA2DNA, STRREP.

%   Copyright 2002-2012 The MathWorks, Inc.


% If the input is a structure then extract the Sequence data.
if isstruct(dna)
    dna = bioinfoprivate.seqfromstruct(dna);
end

if ~ischar(dna)
    rna = dna;
    return
end
if ~isempty(regexpi(dna,'u','once'))
    warning(message('bioinfo:dna2rna:DNAContainsU'));
end

rna = strrep(dna,'T','U');
rna = strrep(rna,'t','u');

if ~isempty(regexpi(rna,'[^ACGURYKMSWBDHVN-]','once'))
    warning(message('bioinfo:dna2rna:UnknownSymbols'));
end




