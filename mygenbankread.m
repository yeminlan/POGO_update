function data=mygenbankread(gbtext,varargin)
%GENBANKREAD reads GenBank format data files.
%
%   DATA = GENBANKREAD(FILE) reads in a GenBank formatted sequence from
%   FILE and creates a structure DATA containing fields corresponding to
%   the GenBank keywords. If the file contains information about multiple
%   sequences, then the information will be stored in an array of
%   structures.
%
%   FILE can also be a URL or a MATLAB character array that contains the
%   text of a GenBank format file.
%
%   Based on version 179.0 of GenBank
%
%   Examples:
%
%       % Download a GenBank file to your local drive.
%       getgenbank('M10051', 'TOFILE', 'HGENBANKM10051.GBK')
%
%       % Then bring it into a MATLAB sequence.
%       data = genbankread('HGENBANKM10051.GBK')
%
%   See also EMBLREAD, FASTAREAD, GENPEPTREAD, GETGENBANK, SCFREAD, 
%   SEQVIEWER.

% Copyright 2002-2012 The MathWorks, Inc.


getPreambleText = false;
allFeatures = false;
% process input arguments
if  nargin > 1
    if rem(nargin,2) == 0
        error(message('bioinfo:genbankread:IncorrectNumberOfArguments', mfilename));
    end
    okargs = {'preambletext','allfeatures'};
    for j=1:2:nargin-1
        pname = varargin{j};
        pval = varargin{j+1};
        k = find(strncmpi(pname, okargs,numel(pname)));
        if isempty(k)
            error(message('bioinfo:genbankread:UnknownParameterName', pname));
        elseif length(k)>1
            error(message('bioinfo:genbankread:AmbiguousParameterName', pname));
        else
            switch(k)
                case 1  % 'preambletext'
                    getPreambleText = bioinfoprivate.opttf(pval);
                case 2  % 'allFeatures'
                    allFeatures = bioinfoprivate.opttf(pval);
            end
        end
    end
end

if ~ischar(gbtext) && ~iscellstr(gbtext)
    error(message('bioinfo:genbankread:InvalidInput'));
end

if iscellstr(gbtext)
    % do not mess with it, just put it to char and try to decipher in the try-catch trap
    gbtext=char(gbtext);
else %it is char, lets check if it has an url or a file before try to decipher it
    if size(gbtext,1)==1 && ~isempty(strfind(gbtext(1:min(10,end)), '://'))
        % must be a URL
        if (~usejava('jvm'))
            error(message('bioinfo:genbankread:NoJava'))
        end
        try
            gbtext = urlread(gbtext);
        catch allExceptions
            error(message('bioinfo:genbankread:CannotReadURL', gbtext));
        end
        % clean up any &amp s
        gbtext=strrep(gbtext,'&amp;','&');
        % make each line a separate row in a string array
        gbtext = textscan(gbtext,'%s','delimiter','\n','whitespace','');
        gbtext = char(gbtext{1});
    elseif size(gbtext,1)==1
        if (exist(gbtext,'file') || exist(fullfile(pwd,gbtext),'file'))
            fid = fopen(gbtext,'r');
            gbtext = textscan(fid,'%s','delimiter','\n','whitespace','');
            gbtext = char(gbtext{1});
            fclose(fid);   
        else
            gbtext = textscan(gbtext,'%s','delimiter','\n','whitespace','');
            gbtext = char(gbtext{1});
        end
    end
end

% If the input is a string of GenBank data then words LOCUS and DEFINITION must be present
if size(gbtext,1)==1 || isempty(strfind(gbtext(1,:),'LOCUS')) || isempty(strfind(gbtext(2,:),'DEFINITION'));
    error(message('bioinfo:genbankread:NonMinimumRequiredFields'))
end

%line number
ln = 1;

%multiple records possible in one record
record_count=1;

numLines = size(gbtext,1);
try
    while 1,

        %LOCUS - Mandatory keyword/exactly one record.
        data(record_count).LocusName = strtrim(gbtext(ln,13:28));  %#ok<*AGROW>
        data(record_count).LocusSequenceLength =strtrim(gbtext(ln,30:40));
        data(record_count).LocusNumberofStrands = strtrim(gbtext(ln,45:47));
        data(record_count).LocusTopology = strtrim(gbtext(ln,56:63));
        data(record_count).LocusMoleculeType = strtrim(gbtext(ln,48:53));
        data(record_count).LocusGenBankDivision = strtrim(gbtext(ln,65:67));
        data(record_count).LocusModificationDate = strtrim(gbtext(ln,69:79));

        ln=ln+1;

        %DEFINITION - Mandatory keyword/one or more records.
        [~,~,t] = regexp(gbtext(ln,:),'DEFINITION\s+(\w|\W)+'); 
        data(record_count).Definition = strtrim(gbtext(ln,t{1}(1):t{1}(2)));

        ln=ln+1;

        while ~matchstart(gbtext(ln,:),'ACCESSION')
            data(record_count).Definition = [data(record_count).Definition,' ', strtrim(gbtext(ln,:))];
            ln = ln+1;
        end

        %ACCESSION - Mandatory keyword/one or more records.
        [~,~,t] = regexp(gbtext(ln,:),'ACCESSION\s+(\w|\W)+'); 
        data(record_count).Accession = strtrim(gbtext(ln,t{1}(1):t{1}(2)));
        ln=ln+1;

        while ~matchstart(gbtext(ln,:),'VERSION')
            data(record_count).Accession=[data(record_count).Accession ' ' strtrim(gbtext(ln,:))];
            ln=ln+1;
        end

        %VERSION - Mandatory keyword/exactly one record.
        [~,~,t] = regexp(gbtext(ln,:),'VERSION\s+(\w|\W)+'); 
        rmdr = gbtext(ln,t{1}(1):t{1}(2));
        [data(record_count).Version, rmdr] = strtok(rmdr,'GI:'); %#ok<STTOK>
        data(record_count).Version = strtrim(data(record_count).Version);

        %GI - Mandatory (part of version)
        data(record_count).GI = deblank(rmdr(4:end));

        ln=ln+1;

        % NID - Obsolete since 12/1999
        
        %PROJECT - Introduced by NCBI in Feb 2006, treated as optional.
        %          Obsoleted in Release 171.0 April 2009
        data(record_count).Project=[];
        [s,~,t] = regexp(gbtext(ln,:),'PROJECT\s+(\w|\W)+'); 
        if ~isempty(s)
            data(record_count).Project=deblank(gbtext(ln,t{1}(1):t{1}(2)));
            ln=ln+1;
        end
        
        %DBLINK - New field introduced by NCBI in Feb 2009, treated as
        %         Optional keyword/one or more records.
        data(record_count).DBLink = [];
        [s,~,t] = regexp(gbtext(ln,:),'DBLINK\s+(\w|\W)+'); 
        if ~isempty(s)
            data(record_count).DBLink=deblank(gbtext(ln,t{1}(1):t{1}(2)));
            ln=ln+1;
        end
        
        while ~matchstart(gbtext(ln,:),'KEYWORDS')
            data(record_count).DBLink=[data(record_count).DBLink ' ' strtrim(gbtext(ln,:))];
            ln=ln+1;            
        end

        %KEYWORDS - Mandatory keyword in all annotated entries/one or more
        %           records. 
        data(record_count).Keywords=[];
        [s,~,t] = regexp(gbtext(ln,:),'KEYWORDS\s+(\w|\W)+'); 
        if ~isempty(s)
            data(record_count).Keywords=deblank(gbtext(ln,t{1}(1):t{1}(2)));
            ln=ln+1;
        end
        while ~isempty(s) && ~matchstart(gbtext(ln,:),'SEGMENT') && ~matchstart(gbtext(ln,:),'SOURCE') 
            data(record_count).Keywords=strvcat(data(record_count).Keywords, deblank(gbtext(ln,:))); %#ok<*VCAT>
            ln=ln+1;
        end
        if all(~isletter(data(record_count).Keywords)),
            data(record_count).Keywords = [];
        end


        %SEGMENT -  Optional keyword (only in segmented entries)/exactly
        %           one record. 

        data(record_count).Segment=[];
        [s,~,t] = regexp(gbtext(ln,:),'SEGMENT\s+(\w|\W)+'); 
        if ~isempty(s)
            data(record_count).Segment=gbtext(ln,t{1}(1):t{1}(2));
            ln=ln+1;
        end


        %SOURCE - Mandatory keyword in all annotated entries/one or
        %         more records/includes one subkeyword.
        [~,~,t] = regexp(gbtext(ln,:),'SOURCE\s+(\w|\W)+'); 
        data(record_count).Source = deblank(gbtext(ln,t{1}(1):t{1}(2)));
        ln=ln+1;
        while ~matchstart(gbtext(ln,:),'ORGANISM') && ~matchstart(gbtext(ln,:),'FEATURES') && ~matchstart(gbtext(ln,:),'COMMENT') && ~matchstart(gbtext(ln,:),'BASE COUNT')
            data(record_count).Source = [data(record_count).Source ' ' deblank(gbtext(ln,t{1}(1):t{1}(2)))];
            ln=ln+1;
        end


        %ORGANISM - Mandatory for all annotated records.
        data(record_count).SourceOrganism = [];
        [s,~,t] = regexp(gbtext(ln,:),'ORGANISM\s+(\w|\W)+'); 
        if ~isempty(s)
            data(record_count).SourceOrganism = strtrim(gbtext(ln,t{1}(1):t{1}(2)));
            ln=ln+1;
            while ~matchstart(gbtext(ln,:),'REFERENCE') && ~matchstart(gbtext(ln,:),'COMMENT')
                data(record_count).SourceOrganism = strvcat(data(record_count).SourceOrganism, strtrim(gbtext(ln,:))); 
                ln=ln+1;
            end
        end

        %REFERENCE
        [data,gbtext,ln] = referenceparse(data,gbtext,ln,record_count);

        %COMMENT - Optional
        data(record_count).Comment = [];
        [s,~,t] = regexp(gbtext(ln,:),'COMMENT\s+(\w|\W)+'); 
        if ~isempty(s)
            data(record_count).Comment = strtrim(gbtext(ln,t{1}(1):t{1}(2)));
            ln=ln+1;
            while ~matchstart(gbtext(ln,:),'FEATURES') && ~matchstart(gbtext(ln,:),'BASE COUNT')...
                    && ~matchstart(gbtext(ln,:),'ORIGIN')
                data(record_count).Comment=strvcat(data(record_count).Comment, strtrim(gbtext(ln,t{1}(1):t{1}(2)))); 
                ln=ln+1;
            end
        end

        %FEATURES - Optional
        data(record_count).Features = [];
        lnFeatures = inf;
        if matchstart(gbtext(ln,:),'FEATURES')
            lnFeatures = ln; % save this position to get preamble text later
            feats = cell(numLines-ln,1);
            ln=ln+1;
            featCount = 1;
            feats{featCount} = gbtext(ln,1:end);
            ln=ln+1;
            while ~matchstart(gbtext(ln,:),'ORIGIN')
                featCount = featCount+1;
                feats{featCount} = gbtext(ln,1:end);
                ln=ln+1;
            end
            data(record_count).Features = strtrim(strvcat(feats(1:featCount))); 
            % Extract information for the features -- CDS, gene, mRNA,
            % misc_feature and repeat_unit
            try
                data(record_count).CDS = extract_feature(data(record_count).Features,'CDS');
            catch allExceptions %#ok<NASGU>
                if isfield(data(record_count),'CDS')
                    data(record_count).CDS = [];
                end
                warning(message('bioinfo:genbankread:BADCDS'));
            end
            if allFeatures
                try
                    data(record_count).mRNA = extract_feature(data(record_count).Features,'mRNA');
                catch allExceptions %#ok<NASGU>
                    if isfield(data(record_count),'mRNA')
                        data(record_count).mRNA = [];
                    end
                    warning(message('bioinfo:genbankread:BADmRNA'));
                end
                try
                    data(record_count).gene = extract_feature(data(record_count).Features,'gene');
                catch alLExceptions %#ok<NASGU>
                    if isfield(data(record_count),'gene')
                        data(record_count).gene = [];
                    end
                    warning(message('bioinfo:genbankread:BADgene'));
                end
                try
                    data(record_count).misc_feature = extract_feature(data(record_count).Features,'misc_feature');
                catch allExceptions %#ok<NASGU>
                    if isfield(data(record_count),'misc_feature')
                        data(record_count).misc_feature = [];
                    end
                    warning(message('bioinfo:genbankread:BADmisc_feature'));
                end
                try
                    data(record_count).repeat_unit = extract_feature(data(record_count).Features,'repeat_unit');
                catch allExceptions %#ok<NASGU>
                    if isfield(data(record_count),'repeat_unit')
                        data(record_count).repeat_unit = [];
                    end
                    warning(message('bioinfo:genbankread:BADrepeat_unit'));
                end
            end
        end

        %BASECOUNT - obsolete, removed October 2003
        
        %ORIGIN - Mandatory
        if matchstart(gbtext(ln,:),'ORIGIN')
            lnOrigin = ln; % save this position to get preamble text later
            ln=ln+1;
        end

        % SEQUENCE
        data(record_count).Sequence = '';
        startln = ln;
        % the sequence will go up to the start of the next possible record
        ln = find(gbtext(ln:end,1) == '/' | gbtext(ln:end,1) == 'L',1) + ln -1;
        if isempty(ln)
            ln = size(gbtext,1);
        end
        endln = ln-1;
        seq = gbtext(startln:endln,:)';
        seq = seq(:)';
        seq(~isletter(seq)) = '';
        data(record_count).Sequence = seq;

        % PREAMBLETEXT
        if getPreambleText
            data(record_count).PreambleText = gbtext(1:min(lnOrigin,lnFeatures)-1,:);
        end

        if ln < numLines && matchstart(gbtext(ln,:),'//')
            while ln<numLines
                ln=ln+1;
                % another record ?
                if matchstart(gbtext(ln,:),'LOCUS')
                    record_count = record_count+1;
                    break
                end
            end
        end
        if ln == numLines
            return
        end
    end
catch le
    if strcmpi(le.identifier,'matlab:nomem')
        clear data
        rethrow(le)
    else
        warning(message('bioinfo:genbankread:incompleteGenBankData'));
    end
end

%-------------------------------------------------------------------------%
function theStruct = extract_feature(theText,theFeature)
% Extract the feature information
% typically CDS, gene, mRNA
theCellstr = cellstr(theText);

% look for the tag at the start of a line
feat = strmatch(theFeature,strtrim(theCellstr)); 

featCount = numel(feat);

if featCount == 0
    theStruct = [];
    return;
end

% As the items in the CDS fields is unknown, we rely on indentation to tell
% when the CDS field ends. This is not particularly robust.
indent = find(theCellstr{feat(1)} == theFeature(1),1)-1;
if indent>0
    theCellstr = cellstr(theText(:,indent+1:end));
end
%look for lines with first level features.
featureLines = find(cellfun('isempty',regexp(theCellstr,'^\s')));
featureLines(end+1) = numel(theCellstr)+1;

% Sometimes 'CDS' show up in the /note as a single line
theFeat = intersect(featureLines,feat);
featCount= numel(theFeat);
% create empty struct
theStruct(featCount).location = '';
theStruct(featCount).gene = '';
theStruct(featCount).product = '';
theStruct(featCount).codon_start = [];
theStruct(featCount).indices = [];
theStruct(featCount).protein_id = '';
theStruct(featCount).db_xref = '';
theStruct(featCount).note = '';
theStruct(featCount).translation = '';
theStruct(featCount).text = '';

% loop over all of the CDS
for featloop = 1:featCount

    featurePos = find(featureLines == theFeat(featloop));

    endLine = featureLines(featurePos+1)-1;
    textChunk = strtrim(theCellstr(theFeat(featloop):endLine));

    numLines = numel(textChunk);
    theStruct(featloop).text = char(textChunk);

    textChunk{1} = strtrim(strrep(textChunk{1},theFeature,''));
    [featLocation, startLoop] = getFullText(textChunk, 1);
    theStruct(featloop).location = char(strread([featLocation{:}], '%s'));
    theStruct(featloop).indices = featurelocation(theStruct(featloop).location);

    strLoop = startLoop;
    while(strLoop < numLines)
        [fullstr, endpos] = getFullText(textChunk, strLoop);
        lines = size(fullstr, 1);
        [token,rest] = strtok(fullstr{1},'='); 
        rest= strrep(rest,'"','');
        fullstr{1}= rest(2:end);
        if(lines > 1)
            fullstr{lines} = strrep(fullstr{lines}, '"','');
        end
        switch token(2:end)
            case 'gene'
                theStruct(featloop).gene = char(fullstr{:});
            case 'product'
                theStruct(featloop).product =  char(fullstr{:});
            case 'codon_start'
                theStruct(featloop).codon_start =  char(fullstr{:});
            case 'protein_id'
                theStruct(featloop).protein_id = char(fullstr{:});
            case 'db_xref'
                theStruct(featloop).db_xref = char(fullstr{:});
            case 'note'
                theStruct(featloop).note = char(fullstr{:});
            case 'translation'
                theStruct(featloop).translation = char(strread([fullstr{:}],'%s'));
            otherwise % There may be other fields...
                % disp(sprintf('Unknown field: %s',token));
        end
        strLoop = endpos;
    end
end
%-------------------------------------------------------------------------%
function [fullText, endPos] = getFullText(textChunk, i)
% next qualifier (if exists) starts with '/'
nextKey = find(strncmp('/',textChunk(i+1:end),1));
if isempty(nextKey)
    endPos = numel(textChunk)+1;
else
    endPos = i + nextKey(1);
end
fullText = textChunk(i:endPos-1);



