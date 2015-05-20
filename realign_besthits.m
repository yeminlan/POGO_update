function realign_besthits(uidA,uidB)

query1 = myfastaread(['/scratch/yl497/CDS/',uidA,'.faa']);
query2 = myfastaread(['/scratch/yl497/CDS/',uidB,'.faa']);
hits = textread(['/scratch/yl497/ortholog/',uidA,'_vs_',uidB,'.ortholog']);%[seqID_in_uidA,seqID_in_uidB,i_sw,c_sw]

%swalign
for i=1:size(hits,1)
        seq1 = query1(hits(i,1)).Sequence;
        seq2 = query2(hits(i,2)).Sequence;
        seq1(~ismember(seq1,'ARNDCQEGHILKMFPSTWYV')) = 'X';
        seq2(~ismember(seq2,'ARNDCQEGHILKMFPSTWYV')) = 'X';
        [~,t4] = myswalign(seq1,seq2);
        hits(i,3) = sum(t4(2,:)=='|')./size(t4,2);
        hits(i,4) = sum(t4(1,:)~='-')/min(numel(seq1),numel(seq2));
        clear t4 seq1 seq2
end
clear i query1 query2

eval(['save realign/',uidA,'_vs_',uidB,'.mat'])
