function out = extension(sourcePointArray,targetPointArray,matrix,imageValues)

    findingMatrix = matrix;
    [height, width] = size(matrix);
    sourcePoints = sourcePointArray;
    targetPoints = targetPointArray;
    size1 = length(sourcePoints);
    temp1 = zeros(3); temp2 = zeros(3);
    if ( size1 ~= 0 )
        source.x = sourcePoints(1,1);
        source.y = sourcePoints(1,2);
        target.x = targetPoints(1,1);
        target.y = targetPoints(1,2);
        if(source.x+1 <= height-1 && source.x-1 >= 1 && source.y+1 <= width-1 &&source.y-1 >= 1)
            for i = -1 : 1
                for j = -1 : 1
                    temp1(i+2,j+2) = imageValues(source.x+i,source.y+j);
                end
            end
        end
        if(target.x+1 <= height-1 && target.x-1 >= 1 && target.y+1 <= width-1 && target.y-1 >= 1 )
            for i = -1 : 1
                for j = -1 : 1
                    temp2(i+2,j+2) = imageValues(target.x+i,target.y+j);
                end
            end
        end
        deger = zncc(temp1,temp2);
        clear temp1;
        clear temp2;
        if deger == 1
            [komsuluk,liste1,liste2] = komsular(sourcePoints(1,:),targetPoints(1,:),findingMatrix);
            sourcePointArray(1,:) = [];
            targetPointArray(1,:) = [];
            if(~isempty(liste1))
                denemeListe = [sourcePointArray ; liste1 ];
                denemeListe2 = [targetPointArray ; liste2 ];
                sonuc = extension(denemeListe,denemeListe2,komsuluk,imageValues);
            else
                denemeListe = sourcePointArray ;
                denemeListe2 = targetPointArray ;
                sonuc = extension(denemeListe,denemeListe2,komsuluk,imageValues);
%                sonuc = komsuluk;
            end
            
        else
            denemeListe = sourcePointArray ;
            denemeListe2 = targetPointArray ;
            denemeListe(1,:) = [];
            denemeListe2(1,:) = [];
            sonuc = extension(denemeListe,denemeListe2,findingMatrix,imageValues);
        end
        
    else
        sonuc = findingMatrix;
    end
    
    out = sonuc;
end