function [liste1, kaynakListe, hedefListe] = komsular(source, target, matris)

    sourceListe(1,1) = source(1,1)-1;
    sourceListe(1,2) = source(1,2)-1;
    sourceListe(2,1) = source(1,1);
    sourceListe(2,2) = source(1,2)-1;
    sourceListe(3,1) = source(1,1)+1 ;
    sourceListe(3,2) = source(1,2)-1;
    sourceListe(4,1) = source(1,1)-1;
    sourceListe(4,2) = source(1,2);
    sourceListe(5,1) = source(1,1)+1;
    sourceListe(5,2) = source(1,2);
    sourceListe(6,1) = source(1,1)-1;
    sourceListe(6,2) = source(1,2)+1;
    sourceListe(7,1) = source(1,1);
    sourceListe(7,2) = source(1,2)+1;
    sourceListe(8,1) = source(1,1)+1;
    sourceListe(8,2) = source(1,2)+1;
    
    targetListe(1,1) = target(1,1)-1;
    targetListe(1,2) = target(1,2)-1;
    targetListe(2,1) = target(1,1);
    targetListe(2,2) = target(1,2)-1;
    targetListe(3,1) = target(1,1)+1 ;
    targetListe(3,2) = target(1,2)-1;
    targetListe(4,1) = target(1,1)-1;
    targetListe(4,2) = target(1,2);
    targetListe(5,1) = target(1,1)+1;
    targetListe(5,2) = target(1,2);
    targetListe(6,1) = target(1,1)-1;
    targetListe(6,2) = target(1,2)+1;
    targetListe(7,1) = target(1,1);
    targetListe(7,2) = target(1,2)+1;
    targetListe(8,1) = target(1,1)+1;
    targetListe(8,2) = target(1,2)+1;
    
    matris1=matris;
    mainSourceList = zeros(8,2);
    mainTargetList = zeros(8,2);
    for i = 1 : 8
        if(matris1(sourceListe(i,1),sourceListe(i,2)) == 1)
            matris1(sourceListe(i,1),sourceListe(i,2)) = 0;
            matris1(targetListe(i,1),targetListe(i,2)) = 0;
            if(~(sourceListe(i,1) == 0 && sourceListe(i,2) == 0))
                mainSourceList(i,:) = sourceListe(i,:);
                mainTargetList(i,:) = targetListe(i,:);
            end
        end
    end
    matris1(source(1,1),source(1,2)) = -1;
    matris1(target(1,1),target(1,2)) = -1;
    listemiz = zeros(8,1);
    for i = 1 : 8
        if(mainSourceList (i,1) == 0 && mainSourceList(i,2) == 0)
            listemiz(i) = i;
        end
    end
    t1 = (find(listemiz>0));
    t = flip(t1);
    for j = 1 : length(t)
        a = t(j);
        mainSourceList(a,:) = [];
        
        mainTargetList(a,:) = [];

    end
    
    liste1 = matris1;
    kaynakListe = mainSourceList;
    hedefListe = mainTargetList;
end