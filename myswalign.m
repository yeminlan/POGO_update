function [score, alignment, startat, matrices] = myswalign(seq1, seq2, varargin)
%SWALIGN performs Smith-Waterman local alignment of two sequences.
%
%   SWALIGN(SEQ1, SEQ2) returns the score (in bits) for the optimal
%   alignment. Note: The scale factor used to calculate the score is
%   provided by the scoring matrix info (see below). If this is not
%   defined, then SWALIGN returns the raw score.
%
%   [SCORE, ALIGNMENT] = SWALIGN(SEQ1, SEQ2) returns a string showing an
%   optimal local alignment of sequences SEQ1 and SEQ2.
%
%   [SCORE, ALIGNMENT, STARTAT] = SWALIGN(SEQ1, SEQ2) returns a 2x1 vector
%   with the starting point indices indicating the starting point of the
%   alignment in the two sequences.
%
%   SWALIGN(..., 'ALPHABET', A) specifies whether the sequences are
%   amino acids ('AA') or nucleotides ('NT'). The default is 'AA'.
%
%   SWALIGN(..., 'SCORINGMATRIX', matrix) specifies the scoring matrix to be
%   used for the alignment. The default is BLOSUM50 for amino acids and
%   NUC44 for nucleotides.
%
%   SWALIGN(..., 'SCALE', scale) indicates the scale factor of the scoring
%   matrix to return the score using arbitrary units. If the scoring matrix
%   info also provides a scale factor, then both are used.
%
%   SWALIGN(..., 'GAPOPEN', penalty) specifies the penalty for opening a gap
%   in the alignment. The default gap open penalty is 8.
%
%   SWALIGN(..., 'EXTENDGAP', penalty) specifies the penalty for extending a
%   gap in the alignment. If EXTENDGAP is not specified, then extensions to
%   gaps are scored with the same value as GAPOPEN.
%
%   SWALIGN(..., 'SHOWSCORE', true) displays the scoring space and the
%   winning path.
%
%
%   Examples:
%
%       % Calculate the score in bits and return the local alignment
%       % using the default scoring matrix (BLOSUM50).
%       [score,align] = swalign('VSPAGMASGYD','IPGKASYD')
%
%       % Use user-specified scoring matrix and gap open penalty.
%       [score,align] = swalign('HEAGAWGHEE','PAWHEAE',...
%                               'scoringmatrix',@pam250,'gapopen',5)
%
%       % Return the score in nat units (nats).
%       [score,align] = swalign('HEAGAWGHEE','PAWHEAE','scale',log(2))
%
%       % Display the scoring space and the winning path.
%       swalign('VSPAGMASGYD','IPGKASYD','showscore',true) 
%
%   See also ALIGNDEMO, BLOSUM, LOCALALIGN, NT2AA, NWALIGN, PAM, SEQDOTPLOT,
%   SHOWALIGNMENT.

%   References:
%   R. Durbin, S. Eddy, A. Krogh, and G. Mitchison. Biological Sequence
%   Analysis. Cambridge UP, 1998.
%   Smith, T. F., Waterman, M. S., J. Mol. Biol. (1981) 147:195-197

%   Copyright 2002-2012 The MathWorks, Inc.

%=== check inputs
% bioinfochecknargin(nargin, 2, mfilename);
[scoringMatrix, gapopen, gapextend, setGapExtend, scale, ...
	isAminoAcid, showscore] = parse_inputs(varargin{:});

%=== the Sequence data from structure
if isstruct(seq1)
    seq1 = bioinfoprivate.seqfromstruct(seq1);
end
if isstruct(seq2)
    seq2 = bioinfoprivate.seqfromstruct(seq2);
end

%=== handle properly "?" characters typically found in pdb files
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

%=== check input sequences
% if isAminoAcid && ~(bioinfoprivate.isaa(seq1) && bioinfoprivate.isaa(seq2))
%     error(message('bioinfo:swalign:InvalidAminoAcidSequences'));
% elseif ~isAminoAcid && ~(bioinfoprivate.isnt(seq1) && bioinfoprivate.isnt(seq2))
%     error(message('bioinfo:swalign:InvalidNucleotideSequences'));
% end

%=== use numerical arrays for easy indexing
if ischar(seq1)
    seq1 = upper(seq1); %the output alignment will be all uppercase
    if isAminoAcid
        intseq1 = myaa2int(seq1);
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
        intseq2 = myaa2int(seq2);
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
    error(message('bioinfo:swalign:InvalidLengthSequences'));
end

%=== set the scoring matrix
if isempty(scoringMatrix)
    if isAminoAcid
        [scoringMatrix,scoringMatrixInfo] = blosum50;
    else
        [scoringMatrix,scoringMatrixInfo] = nuc44;
    end
else
	if ~isnumeric(scoringMatrix)
		try
			[scoringMatrix, scoringMatrixInfo] = feval(scoringMatrix);
		catch allExceptions
			error(message('bioinfo:swalign:InvalidScoringMatrix'));
		end
		
	end
end


%=== get the scale from scoringMatrixInfo, if it exists
if exist('scoringMatrixInfo','var') && isfield(scoringMatrixInfo,'Scale')
    scale = scale * scoringMatrixInfo.Scale;
end

%=== make sure scoringMatrix can handle unknown, ambiguous or gaps

% possible values are
% B  Z  X  * -
% 21 22 23 24 25

scoringMatrixSize = size(scoringMatrix,1);
highestVal = max([intseq1, intseq2]);
if highestVal > scoringMatrixSize
	
    %=== map to the symbol corresponding to 'Any', if the matrix allows it
    if isAminoAcid
        anyVal = myaa2int('X');
    else
        anyVal = mynt2int('N');
    end
    if scoringMatrixSize >= anyVal
        intseq1(intseq1>scoringMatrixSize) = anyVal;
        intseq2(intseq2>scoringMatrixSize) = anyVal;
    else
     error(message('bioinfo:swalign:InvalidSymbolsInInputSequences'));
    end
end

%=== call core functions

if setGapExtend  
% 	[F, pointer] = affinegap(intseq1,m,intseq2,n,scoringMatrix,gapopen,gapextend);
	
	if nargout ==  4 % return score matrices and pointer matrices
		[score, path(:,2), path(:,1), F(:,:,1), F(:,:,2), F(:,:,3), pointer] = ...
			affinegapmex(intseq2, intseq1, gapopen, gapextend, scoringMatrix, 2);
		pointer = shiftdim(pointer,1);
	elseif  showscore % return the score matrices as well
		[score, path(:,2), path(:,1), F(:,:,1), F(:,:,2), F(:,:,3)] = ...
			affinegapmex(intseq2, intseq1, gapopen, gapextend, scoringMatrix, 2);
		
	else % return only score and alignment
		[score, path(:,2), path(:,1)] = affinegapmex(intseq2, intseq1, ...
			gapopen, gapextend, scoringMatrix, 2);
	end

else
    if nargout ==  4 % return score matrices and pointer matrices
     [score, path(:,2), path(:,1), F, pointer] = ...
             simplegapmex(intseq2, intseq1, gapopen, scoringMatrix, 2); 
	
	elseif showscore % return score matrices
        [score, path(:,2), path(:,1), F] = ...
             simplegapmex(intseq2, intseq1, gapopen, scoringMatrix, 2); 
      
    else % return only score and alignment
          [score, path(:,2), path(:,1)] = ...
             simplegapmex(intseq2, intseq1, gapopen, scoringMatrix, 2); 
    end
end

path = path(sum(path,2)>0,:);
path = flipud(path);

%=== re-scale the output score
score = scale * score;

if nargout <= 1 && ~showscore
    return
end

%=== warn if empty alignment
if size(path,1) == 0 
    alignment = char(zeros(3,0));
    startat = zeros(2,1);
    warning(message('bioinfo:swalign:EmptyAlignment'))

else         

    %=== initialize outputs
    alignment = repmat(('- -')',1,size(path,1));

    %=== add sequences to alignment
    alignment(1, path(:,1)>0) = seq1(path(path(:,1)>0, 1));
    alignment(3, path(:,2)>0) = seq2(path(path(:,2)>0, 2));

    %=== find positions where there are no gaps
    h = find(all(path>0,2));
    if isAminoAcid
        noGaps1 = myaa2int(alignment(1,h));
        noGaps2 = myaa2int(alignment(3,h));
    else
        noGaps1 = mynt2int(alignment(1,h));
        noGaps2 = mynt2int(alignment(3,h));
    end

    %=== erase symbols that cannot be scored
    htodel = max([noGaps1;noGaps2]) > scoringMatrixSize;
    h(htodel) = [];
    noGaps1(htodel) = [];
    noGaps2(htodel) = [];

    %=== score pairs with no gap
    value = scoringMatrix(sub2ind(size(scoringMatrix),double(noGaps1),double(noGaps2)));

    %=== insert symbols of the match string into the alignment
    alignment(2,h(value >= 0)) = ':';
    alignment(2,h(noGaps1 == noGaps2)) = '|';

    startat = path(1,:)';
end

% undocumented fourth output -- score and pointer matrices
if nargout > 3
    matrices.Scores = F;
    matrices.Pointers = pointer;
end

if showscore
    figure
    F = scale .* max(F(2:end,2:end,:),[],3);
    clim = max(max(max(abs(F(~isinf(F))))),eps);
    imagesc(F,[-clim clim]);
    colormap(privateColorMap(1));
    set(colorbar,'YLim',[min([F(:);-eps]) max([F(:);eps])])
    title('Scoring Space and Winning Path')
    xlabel('Sequence 1')
    ylabel('Sequence 2')
    hold on
    plot(path(all(path>0,2),1),path(all(path>0,2),2),'k.')
end


%==========================================================================
% SIMPLEGAP is now a mex function --> simplegapmex
%
% function [F, pointer] = simplegap(intseq1,m,intseq2,n,scoringMatrix,gap)
% % Standard Smith-Waterman algorithm
% 
% % set up storage for dynamic programming matrix
% F = zeros(n+1,m+1);
% 
% % and for the back tracing matrix
% pointer(n+1,m+1) = uint8(0);
% 
% % initialize buffers to the first column
% currentFColumn  = zeros(n+1,1);
% ptr(n+1,1)      = uint8(0);
% 
% % main loop runs through the matrix looking for maximal scores
% for outer = 2:m+1
%     
%     % score current column
%     scoredMatchColumn = scoringMatrix(intseq2,intseq1(outer-1));
%     % grab the data from the matrices and initialize some values 
%     lastFColumn    = currentFColumn;
%     currentFColumn = F(:,outer);
%     best = 0;
%            
%     for inner = 2:n+1
%         % score the three options
%         up       = best + gap;
%         left     = lastFColumn(inner) + gap;   
%         diagonal = lastFColumn(inner-1) + scoredMatchColumn(inner-1);
% 
%         % Now figure out which has the highest non negative score
%         best = 0;
%         pos = 0;
% 
%         % max could be used here but it is quicker to use if statements
%         if up > 0 && up > left
%             best = up ;
%             pos = 2;
%         elseif left > 0
%             best = left;
%             pos = 4;
%         end
%         if diagonal >= best
%             best = diagonal;
%             ptr(inner) = 1;
%         else
%             ptr(inner) = pos;
%         end
% 
%         currentFColumn(inner) = best;
%         
%     end % inner
%     % put back updated columns
%     F(:,outer)   = currentFColumn;
%     % save columns of pointers
%     pointer(:,outer)  = ptr;
% end % outer

%=========================================================================
% AFFINEGAP is now a mex function --> affinegapmex
%
% function [F,pointer] = affinegap(intseq1,m,intseq2,n,scoringMatrix,gapopen,gapextend)
% % Smith-Waterman  algorithm modified to handle affine gaps
% 
% % Set states
% inAlign =   1;
% inGapUp =   2;
% inGapLeft = 3;
% numStates = 3;
% 
% % Set up storage for dynamic programming matrix:
% % for keeping the maximum scores for every state
% F = zeros(n+1,m+1,numStates);
% 
% % and for the back tracing pointers
% pointer(n+1,m+1,numStates) = uint8(0);
% 
% % initialize buffers to the first column
% ptrA(n+1) = uint8(0);
% ptrU(n+1) = uint8(0);
% ptrL(n+1) = uint8(0);
% 
% currentFColumnA = F(:,1,inAlign);
% currentFColumnU = F(:,1,inGapUp);
% currentFColumnL = F(:,1,inGapLeft);
% 
% % main loop runs through the matrix looking for maximal scores
% for outer = 2:m+1
%     % score current column
%     scoredMatchColumn = scoringMatrix(intseq2,intseq1(outer-1));
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
%         if bestU > 0
%             currentFColumnU(inner)=bestU;
%             ptrU(inner)=ptr;
%         else
%             bestU = 0;
%             ptrU(inner)=0;
%         end 
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
%         if bestL > 0
%             currentFColumnL(inner) = bestL;
%             ptrL(inner) = ptr;
%         else
%             ptrL(inner) = 0;
%         end 
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
%         if bestA > 0
%             currentFColumnA(inner) = bestA;
%             ptrA(inner) = ptr;
%         else
%             bestA = 0;
%             ptrA(inner) = 0;
%         end
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

%==========================================================================
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
xcl = 1;
for i = 1:size(pts,1)-1
    xcl = [xcl,i+1/pts(i,4):1/pts(i,4):i+1]; %#ok<AGROW>
end
pcmap = interp1(pts(:,1:3),xcl);

%==========================================================================
function [scoringMatrix, gapopen, gapextend, setGapExtend, scale, ...
	isAminoAcid, showscore] = parse_inputs(varargin)
%Parse PV Pairs

%=== defaults
gapopen = -8;
gapextend = -8;
setGapExtend = false;
showscore = false;
isAminoAcid = true;
scale = 1;
scoringMatrix = [];

%=== check for number of inputs
if rem(nargin,2) == 1
	error(message('bioinfo:swalign:IncorrectNumberOfArguments', mfilename));
end

%=== allowed parameters
okargs = {'scoringmatrix','gapopen','extendgap','alphabet','scale','showscore'};
    
for j = 1:2:nargin
	pname = varargin{j};
	pval = varargin{j+1};
	k = find(strncmpi(pname, okargs,numel(pname)));
	if isempty(k)
		error(message('bioinfo:swalign:UnknownParameterName', pname));
	elseif length(k)>1
		error(message('bioinfo:swalign:AmbiguousParameterName', pname));
	else
		switch(k)
			case 1  % scoring matrix
				if isnumeric(pval)
					scoringMatrix = pval;
				else
					scoringMatrix = lower(pval);
				end
				
			case 2 % gap open penalty
				gapopen = -pval;
				
			case 3 % gap extend penalty
				gapextend = -pval;
				setGapExtend = true;
				
			case 4 % alphabet
                isAminoAcid = bioinfoprivate.optAlphabet(pval,okargs{k}, mfilename);
				
			case 5 % scale
				scale = pval;
				
			case 6 % showscore
				showscore = bioinfoprivate.opttf(pval, okargs{k}, mfilename);
		end
	end
end



