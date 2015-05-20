load data.mat
load ../codes/ncbi2014.mat

f=fopen('taxonomy.txt');
c=0;
while 1
  line=fgets(f);
  c=c+1;
  if all(line==-1)
    break
  end
  t1 = regexp(line,'\t','split');
  t2 = t1(ismember(t1,Species));
  if ~isempty(t2)
    new_taxonomy.species{c,1} = t2{end};
  else
    new_taxonomy.species{c,1} = 'NA';
  end
  t2 = t1(ismember(t1,Genus));
  if ~isempty(t2)
    new_taxonomy.genus{c,1} = t2{end};
  else
    new_taxonomy.genus{c,1} = 'NA';
  end
  t2 = t1(ismember(t1,Family));
  if ~isempty(t2)
    new_taxonomy.family{c,1} = t2{end};
  else
    new_taxonomy.family{c,1} = 'NA';
  end
  t2 = t1(ismember(t1,Order));
  if ~isempty(t2)
    new_taxonomy.order{c,1} = t2{end};
  else
    new_taxonomy.order{c,1} = 'NA';
  end
  t2 = t1(ismember(t1,Class));
  if ~isempty(t2)
    new_taxonomy.class{c,1} = t2{end};
  else
    new_taxonomy.class{c,1} = 'NA';
  end
  t2 = t1(ismember(t1,Phylum));
  if ~isempty(t2)
    new_taxonomy.phylum{c,1} = t2{end};
  else
    new_taxonomy.phylum{c,1} = 'NA';
  end
  t2 = t1(ismember(t1,Superkingdom));
  if ~isempty(t2)
    new_taxonomy.superkingdom{c,1} = t2{end};
  else
    new_taxonomy.superkingdom{c,1} = 'NA';
  end
  clear t1 t2
end
fclose(f);
clear c line f ans

save data.mat

