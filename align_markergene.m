load ../data.mat AAI
startindex = numel(AAI)+1;
clear AAI

load data.mat

for g=2:numel(gene_name)
  f=myfastaread(['markergene/',gene_name{g},'.fna']);
  for i=startindex:numel(AAI)
    disp(['Processing ',num2str(g),'-',num2str(i)]);
    t1=find(strcmp({f.Header},genome_pair{i,1}));
    t2=find(strcmp({f.Header},genome_pair{i,2}));
    if ~isempty(t1) & ~isempty(t2)
      S1=f(t1(end)).Sequence;
      S2=f(t2(end)).Sequence;
      [~,Align1]=mynwalign(S1,S2,'Alphabet','NT');
      [~,Align2]=mynwalign(S1,myseqrcomplement(S2),'Alphabet','NT');
      t1=length(find(Align1(2,:)=='|'))/min(length(S1),length(S2));
      t2=length(find(Align2(2,:)=='|'))/min(length(S1),length(S2));
      gene_identity(i,g)=max(t1,t2);
    else
      gene_identity(i,g)=NaN;
    end
    clear t1 t2 S1 S2 Align1 Align2
  end
  fclose(f);
  clear i ans f
  save data2.mat
end
clear g

save data.mat %replace


