function dnarc = myseqrcomplement(dna)
%SEQRCOMPLEMENT returns the reverse complementary strand of a DNA sequence.
%
%   SEQ = SEQRCOMPLEMENT(DNA) calculates the reverse complementary strand 
%   3' --> 5' (A-->T, C-->G, G-->C, T-->A) of sequence DNA.
%
%   SEQ is returned in the same format as DNA, so if DNA is an integer
%   sequence then so is SEQ.
%
%   Example:
% 
%       % Create a random sequence of nucleotides and the reverse complement.
%       seq = randseq(25)
%       revcompseq = seqrcomplement(seq)
%
%   See also CODONCOUNT, JOINSEQ, PALINDROMES, SEQCOMPLEMENT, SEQREVERSE,
%   SEQVIEWER.

%   Copyright 2002-2012 The MathWorks, Inc.


dnarc = myseqreverse(myseqcomplement(dna));




