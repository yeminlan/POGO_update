function parallel_align(i)

load temp.mat

    for j=1:numel(all_uid)
        disp(['processing pair ',num2str(i),' vs ',num2str(j)])
        if strcmp(new_uid{i},all_uid{j})
                similarity(i,j)=1;
        else
            t1=strcmp({all_rna.Header}',new_uid{i});
            t2=strcmp({all_rna.Header}',all_uid{j});
            t1={all_rna(t1).Sequence}';
            t2={all_rna(t2).Sequence}';
            for m=1:numel(t1)
                for k=1:numel(t2)
                    S1=t1{m};
                    S2=t2{k};
                    [~,Align1] = mynwalign(S1,S2);
                    [~,Align2] = mynwalign(S1,myseqrcomplement(S2));
                    s1=length(find(Align1(2,:)=='|'))/min(length(S1),length(S2));
                    s2=length(find(Align2(2,:)=='|'))/min(length(S1),length(S2));
                    similarity(i,j)=max([similarity(i,j),s1,s2]);
                    clear S1 S2 Align1 Align2 s1 s2
                end
            end
            clear t1 t2 m k
        end
    end

clear j rna all_rna new_rna

eval(['save part',num2str(i,'%.3d'),'.mat similarity'])


