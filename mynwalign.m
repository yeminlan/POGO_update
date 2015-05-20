function [score, alignment, startat, matrices] = mynwalign(seq1,seq2,varargin)
%NWALIGN performs Needleman-Wunsch global alignment of two sequences.
%
%   NWALIGN(SEQ1, SEQ2) returns the score (in bits) for the optimal
%   alignment. Note: The scale factor used to calculate the score is
%   provided by the scoring matrix info (see below). If this is not
%   defined, then NWALIGN returns the raw score.
%
%   [SCORE, ALIGNMENT] = NWALIGN(SEQ1, SEQ2) returns a string showing an
%   optimal global alignment of amino acid (or nucleotide) sequences SEQ1
%   and SEQ2.
%
%   [SCORE, ALIGNMENT, STARTAT] = NWALIGN(SEQ1, SEQ2)  returns a 2x1 vector
%   with the starting point indices indicating the starting point of the
%   alignment in the two sequences. Note: this output is for consistency
%   with SWALIGN and will always be [1;1] because this is a global
%   alignment.
%
%   NWALIGN(..., 'ALPHABET', A) specifies whether the sequences are
%   amino acids ('AA') or nucleotides ('NT'). The default is AA.
%
%   NWALIGN(..., 'SCORINGMATRIX', matrix) defines the scoring matrix to be
%   used for the alignment. The default is BLOSUM50 for AA or NUC44 for NT.
%
%   NWALIGN(..., 'SCALE' ,scale) indicates the scale factor of the scoring
%   matrix to return the score using arbitrary units. If the scoring matrix
%   Info also provides a scale factor, then both are used.
%
%   NWALIGN(..., 'GAPOPEN', penalty) defines the penalty for opening a gap
%   in the alignment. The default gap open penalty is 8.
%
%   NWALIGN(..., 'EXTENDGAP', penalty) defines the penalty for extending a
%   gap in the alignment. If EXTENDGAP is not specified, then extensions to
%   gaps are scored with the same value as GAPOPEN.
%
%   NWALIGN(..., 'GLOCAL', true) performs a semi-global alignment. In a
%   semi-global alignment (also known as GLOCAL) gap penalties at the end
%   of the sequences are null.
%
%   NWALIGN(..., 'SHOWSCORE', true) displays the scoring space and the
%   winning path.
%
%
%   Examples:
%
%       % Return the score in bits and the global alignment using the
%       % default scoring matrix (BLOSUM50).
%       [score, align] = nwalign('VSPAGMASGYD', 'IPGKASYD')
%
%       % Use user-specified scoring matrix and "gap open" penalty.
%       [score, align] = nwalign('IGRHRYHIGG', 'SRYIGRG',...
%                               'scoringmatrix', @pam250, 'gapopen',5)
%
%       % Return the score in nat units (nats).
%       [score, align] = nwalign('HEAGAWGHEE', 'PAWHEAE', 'scale', log(2))
%
%       % Display the scoring space and the winning path.
%       nwalign('VSPAGMASGYD', 'IPGKASYD', 'showscore', true)
%
%   See also ALIGNDEMO, BLOSUM, LOCALALIGN, MULTIALIGN, NT2AA, PAM,
%   PROFALIGN, SEQDOTPLOT, SHOWALIGNMENT, SWALIGN.

%   References:
%   R. Durbin, S. Eddy, A. Krogh, and G. Mitchison. Biological Sequence
%   Analysis. Cambridge UP, 1998.
%   Needleman, S. B., Wunsch, C. D., J. Mol. Biol. (1970) 48:443-453

%   Copyright 2002-2012 The MathWorks, Inc.

glocal = false;
gapopen = -8;
gapextend = -8;
setGapExtend = false;
showscore=false;
isAminoAcid = false;%change!!!
scale=1;

% If the input is a structure then extract the Sequence data.
if isstruct(seq1)
    seq1 = bioinfoprivate.seqfromstruct(seq1);
end
if isstruct(seq2)
    seq2 = bioinfoprivate.seqfromstruct(seq2);
end
if nargin > 2
    if rem(nargin,2) == 1
        error(message('bioinfo:nwalign:IncorrectNumberOfArguments', mfilename));
    end
    okargs = {'scoringmatrix','gapopen','extendgap','alphabet','scale','showscore','glocal'};
    for j=1:2:nargin-2
        pname = varargin{j};
        pval = varargin{j+1};
        k = find(strncmpi(pname, okargs,numel(pname)));
        if isempty(k)
            error(message('bioinfo:nwalign:UnknownParameterName', pname));
        elseif length(k)>1
            error(message('bioinfo:nwalign:AmbiguousParameterName', pname));
        else
            switch(k)
                case 1  % scoring matrix
                    if isnumeric(pval)
                        ScoringMatrix = pval;
                    else
                        if ischar(pval)
                            pval = lower(pval);
                        end
                        try
                            [ScoringMatrix,ScoringMatrixInfo] = feval(pval);
                        catch allExceptions
                            error(message('bioinfo:nwalign:InvalidScoringMatrix'));
                        end
                    end
                case 2 %gap open penalty
                    gapopen = -pval;
                case 3 %gap extend penalty
                    gapextend = -pval;
                    setGapExtend = true;
                case 4 %if sequence is nucleotide
                    isAminoAcid = false;%change!!!  %bioinfoprivate.optAlphabet(pval,okargs{k}, mfilename);
                case 5 % scale
                    scale=pval;
                case 6 % showscore
                    showscore = bioinfoprivate.opttf(pval, okargs{k}, mfilename);
                case 7 % glocal
                    glocal = bioinfoprivate.opttf(pval, okargs{k}, mfilename);
            end
        end
    end
end

% setting the default scoring matrix
if ~exist('ScoringMatrix','var')
    if isAminoAcid
        [ScoringMatrix,ScoringMatrixInfo] = blosum50;
    else
        [ScoringMatrix,ScoringMatrixInfo] = nuc44;
    end
end


% getting the scale from ScoringMatrixInfo, if it exists
if exist('ScoringMatrixInfo','var') && isfield(ScoringMatrixInfo,'Scale')
    scale=scale*ScoringMatrixInfo.Scale;
end

% handle properly "?" characters typically found in pdb files
if isAminoAcid
    if ischar(seq1)
        seq1 = strrep(seq1,'?','X');
    else
        seq1(seq1 == 26) = 23;
    end
    if ischar(seq2)
        seq2 = strrep(seq2,'?','X');
    else
        seq2(seq2 == 26) = 23;
    end
end

% check input sequences
%if isAminoAcid && ~(bioinfoprivate.isaa(seq1) && bioinfoprivate.isaa(seq2))
%    error(message('bioinfo:nwalign:InvalidAminoAcidSequences'));
%elseif ~isAminoAcid && ~(bioinfoprivate.isnt(seq1) && bioinfoprivate.isnt(seq2))
%    error(message('bioinfo:nwalign:InvalidNucleotideSequences'));
%end

% use numerical arrays for easy indexing
if ischar(seq1)
    seq1=upper(seq1); %the output alignment will be all uppercase
    if isAminoAcid
        intseq1 = aa2int(seq1);
    else
        intseq1 = mynt2int(seq1);
    end
else
    intseq1 = uint8(seq1);
    if isAminoAcid
        seq1 = int2aa(intseq1);
    else
        seq1 = int2nt(intseq1);
    end
end
if ischar(seq2)
    seq2 = upper(seq2); %the output alignment will be all uppercase
    if isAminoAcid
        intseq2 = aa2int(seq2);
    else
        intseq2 = mynt2int(seq2);
    end
else
    intseq2 = uint8(seq2);
    if isAminoAcid
        seq2 = int2aa(intseq2);
    else
        seq2 = int2nt(intseq2);
    end
end


m = length(seq1);
n = length(seq2);
if ~n||~m
    error(message('bioinfo:nwalign:InvalidLengthSequences'));
end

% If unknown, ambiguous or gaps appear in the sequence, we need to make
% sure that ScoringMatrix can handle them.

% possible values are
% B  Z  X  *  -  ?
% 21 22 23 24 25 26

scoringMatrixSize = size(ScoringMatrix,1);

highestVal = max([intseq1, intseq2]);
if highestVal > scoringMatrixSize
    % if the matrix contains the 'Any' we map to that
    if isAminoAcid
        anyVal = aa2int('X');
    else
        anyVal = mynt2int('N');
    end
    if scoringMatrixSize >= anyVal
        intseq1(intseq1>scoringMatrixSize) = anyVal;
        intseq2(intseq2>scoringMatrixSize) = anyVal;
    else
        error(message('bioinfo:nwalign:InvalidSymbolsInInputSequences'));
    end
end

if glocal
   algorithm = 3;
else
   algorithm = 1;
end

if setGapExtend 
    % flip order of input sequences for consistency with older versions
    if showscore % return the score matrices
        [score, path(:,2), path(:,1), F(:,:,1), F(:,:,2), F(:,:,3)] = ...
            affinegapmex(intseq2, intseq1, gapopen, gapextend, ScoringMatrix, algorithm); 
        
    elseif nargout ==  4 % return score matrices and pointer matrices
        [score, path(:,2), path(:,1), F(:,:,1), F(:,:,2), F(:,:,3), pointer] = ...
            affinegapmex(intseq2, intseq1, gapopen, gapextend, ScoringMatrix, algorithm); 
        pointer = shiftdim(pointer,1); % for backward compatibility
           
    else % return only score and alignment
        [score, path(:,2), path(:,1)] = affinegapmex(intseq2, intseq1, ...
            gapopen, gapextend, ScoringMatrix, algorithm); 
    end

else
     % flip order of input sequences for consistency with older versions
    if showscore % return the score matrices
        [score, path(:,2), path(:,1), F] = ...
            simplegapmex(intseq2, intseq1, gapopen, ScoringMatrix, algorithm);
        
    elseif nargout == 4
        [score, path(:,2), path(:,1), F, pointer] = ...
            simplegapmex(intseq2, intseq1, gapopen, ScoringMatrix, algorithm);
        
    else
        [score, path(:,2), path(:,1)] = ...
            simplegapmex(intseq2, intseq1, gapopen, ScoringMatrix, algorithm);
    end
end

path = path(sum(path,2)>0,:);
path = flipud(path);

% re-scaling the output score
score = scale * score;

if nargout<=1 && ~showscore
    return
end

% setting the size of the alignment
alignment = repmat(('- -')',1,size(path,1));

% adding sequence to alignment
alignment(1,path(:,1)>0) = seq1;
alignment(3,path(:,2)>0) = seq2;

% find locations where there are no gaps
h=find(all(path>0,2));
if isAminoAcid
    noGaps1=aa2int(alignment(1,h));
    noGaps2=aa2int(alignment(3,h));
else
    noGaps1=mynt2int(alignment(1,h));
    noGaps2=mynt2int(alignment(3,h));
end

% erasing symbols that cannot be scored
htodel=max([noGaps1;noGaps2])>scoringMatrixSize;
h(htodel)=[];
noGaps1(htodel)=[];
noGaps2(htodel)=[];

% score pairs with no gap
value = ScoringMatrix(sub2ind(size(ScoringMatrix),double(noGaps1),double(noGaps2)));

% insert symbols of the match string into the alignment
alignment(2,h(value>=0)) = ':';
alignment(2,h(noGaps1==noGaps2)) = '|';

startat = [1;1];

% undocumented fourth output -- score and pointer matrices
if nargout > 3
    matrices.Scores = F;
    matrices.Pointers = pointer;
end

if showscore
    figure
    F=scale.*max(F(2:end,2:end,:),[],3);
    clim=max(max(max(abs(F(~isinf(F))))),eps);
    imagesc(F,[-clim clim]);
    colormap(privateColorMap(1));
    set(colorbar,'YLim',[min([F(:);-eps]) max([F(:);eps])])
    title('Scoring Space and Winning Path')
    xlabel('Sequence 1')
    ylabel('Sequence 2')
    hold on
    plot(path(all(path>0,2),1),path(all(path>0,2),2),'k.')
end

%=== SIMPLEGAP is now a mex function ===%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [F, pointer] = simplegap(intseq1,m,intseq2,n,ScoringMatrix,gap)
% % Standard Needleman-Wunsch algorithm
% 
% % set up storage for dynamic programming matrix
% F = zeros(n+1,m+1);
% F(2:end,1) = gap * (1:n)';
% F(1,2:end) = gap * (1:m);
% 
% % and for the back tracing matrix
% pointer= repmat(uint8(4),n+1,m+1);
% pointer(:,1) = 2;  % up
% pointer(1,1) = 1;
% 
% 
% % initialize buffers to the first column
% ptr = pointer(:,2); % ptr(1) is always 4
% currentFColumn = F(:,1);
% 
% % main loop runs through the matrix looking for maximal scores
% for outer = 2:m+1
% 
%     % score current column
%     scoredMatchColumn = ScoringMatrix(intseq2,intseq1(outer-1));
%     % grab the data from the matrices and initialize some values
%     lastFColumn    = currentFColumn;
%     currentFColumn = F(:,outer);
%     best = currentFColumn(1);
% 
%     for inner = 2:n+1
%         % score the three options
%         up       = best + gap;
%         left     = lastFColumn(inner) + gap;
%         diagonal = lastFColumn(inner-1) + scoredMatchColumn(inner-1);
% 
%         % max could be used here but it is quicker to use if statements
%         if up > left
%             best = up;
%             pos = 2;
%         else
%             best = left;
%             pos = 4;
%         end
% 
%         if diagonal >= best
%             best = diagonal;
%             ptr(inner) = 1;
%         else
%             ptr(inner) = pos;
%         end
%         currentFColumn(inner) = best;
% 
%     end % inner
%     % put back updated columns
%     F(:,outer)   = currentFColumn;
%     % save columns of pointers
%     pointer(:,outer)  = ptr;
% end % outer
% 


%=== AFFINEGAP is now a mex function ===%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [F,pointer] = affinegap(intseq1,m,intseq2,n,ScoringMatrix,gapopen,gapextend)
% % Needleman-Wunsch algorithm modified to handle affine gaps
% 
% % Set states
% inAlign =   1;
% inGapUp =   2;
% inGapLeft = 3;
% numStates = 3;
% 
% % Set up storage for dynamic programming matrix:
% % for keeping the maximum scores for every state
% 
% F =  zeros(n+1,m+1,numStates);
% F(:,1,:) = -inf;
% F(1,:,:) = -inf;
% F(1,1,inAlign) = 0;
% 
% F(2:end,1,inGapUp)   = gapopen + gapextend * (0:n-1)';
% F(1,2:end,inGapLeft) = gapopen + gapextend * (0:m-1);
% 
% % and for the back tracing pointers
% pointer(n+1,m+1,numStates) = uint8(0);
% pointer(2:end,1,inGapUp)   = 2;  % up
% pointer(1,2:end,inGapLeft) = 4;  % left
% 
% % initialize buffers to the first column
% ptrA = pointer(:,1,inAlign);
% ptrU = pointer(:,1,inGapLeft);
% ptrL = pointer(:,1,inGapUp);
% 
% currentFColumnA = F(:,1,inAlign);
% currentFColumnU = F(:,1,inGapUp);
% currentFColumnL = F(:,1,inGapLeft);
% 
% % main loop runs through the matrix looking for maximal scores
% for outer = 2:m+1
%     % score current column
%     scoredMatchColumn = ScoringMatrix(intseq2,intseq1(outer-1));
%     % grab the data from the matrices and initialize some values for the
%     % first row the most orderly possible
%     lastFColumnA    = currentFColumnA;
%     currentFColumnA = F(:,outer,inAlign);
%     bestA           = currentFColumnA(1);
%     currentinA      = lastFColumnA(1);
% 
%     lastFColumnU    = currentFColumnU;
%     currentFColumnU = F(:,outer,inGapUp);
%     bestU           = currentFColumnU(1);
% 
%     lastFColumnL    = currentFColumnL;
%     currentFColumnL = F(:,outer,inGapLeft);
%     currentinGL     = lastFColumnL(1);
% 
%     for inner = 2:n+1
% 
%         % grab the data from the columns the most orderly possible
%         upOpen      = bestA + gapopen;
%         inA         = currentinA;
%         currentinA  = lastFColumnA(inner);
%         leftOpen    = currentinA + gapopen;
% 
%         inGL        = currentinGL;
%         currentinGL = lastFColumnL(inner);
%         leftExtend  = currentinGL + gapextend;
% 
%         upExtend = bestU + gapextend;
%         inGU     = lastFColumnU(inner-1);
% 
%         % operate state 'inGapUp'
% 
%         if upOpen > upExtend
%             bestU = upOpen; ptr = 1;   % diagonal
%         elseif upOpen < upExtend
%             bestU = upExtend; ptr = 2; % up
%         else % upOpen == upExtend
%             bestU = upOpen; ptr = 3;   % diagonal and up
%         end
%         currentFColumnU(inner)=bestU;
%         ptrU(inner)=ptr;
% 
%         % operate state 'inGapLeft'
% 
%         if leftOpen > leftExtend
%             bestL = leftOpen; ptr = 1;   % diagonal
%         elseif leftOpen < leftExtend
%             bestL = leftExtend; ptr = 4; % left
%         else % leftOpen == leftExtend
%             bestL = leftOpen; ptr = 5;   % diagonal and left
%         end
%         currentFColumnL(inner) = bestL;
%         ptrL(inner) = ptr;
% 
%         % operate state 'inAlign'
% 
%         if  inA > inGU
%             if inA > inGL
%                 bestA = inA; ptr = 1;  % diagonal
%             elseif inGL > inA
%                 bestA = inGL; ptr = 4; % left
%             else
%                 bestA = inA; ptr = 5;  % diagonal and left
%             end
%         elseif inGU > inA
%             if inGU > inGL
%                 bestA = inGU; ptr = 2; % up
%             elseif inGL > inGU
%                 bestA = inGL; ptr = 4; % left
%             else
%                 bestA = inGU; ptr = 6; % up & left
%             end
%         else
%             if inA > inGL
%                 bestA = inA; ptr = 3;  % diagonal & up
%             elseif inGL > inA
%                 bestA = inGL; ptr = 4; % left
%             else
%                 bestA = inA; ptr = 7;  % all
%             end
%         end
% 
%         bestA = bestA + scoredMatchColumn(inner-1);
%         currentFColumnA(inner) = bestA;
%         ptrA(inner) = ptr;
% 
%     end %inner
% 
%     % put back updated columns
%     F(:,outer,inGapLeft) = currentFColumnL;
%     F(:,outer,inGapUp)   = currentFColumnU;
%     F(:,outer,inAlign)   = currentFColumnA;
%     % save columns of pointers
%     pointer(:,outer,inAlign)  = ptrA;
%     pointer(:,outer,inGapUp)  = ptrU;
%     pointer(:,outer,inGapLeft)= ptrL;
% end %outer

function pcmap = privateColorMap(selection)
%PRIVATECOLORMAP returns a custom color map
switch selection
    case 1, pts = [0 0 .3 20;
            0 .1 .8 25;
            0 .9 .5 15;
            .9 1 .9 8;
            1 1 0 26;
            1 0 0 26;
            .4 0 0 0];
    otherwise, pts = [0 0 0 128; 1 1 1 0];
end
xcl=1;
for i=1:size(pts,1)-1
    xcl=[xcl,i+1/pts(i,4):1/pts(i,4):i+1]; %#ok<AGROW>
end
pcmap = interp1(pts(:,1:3),xcl);
