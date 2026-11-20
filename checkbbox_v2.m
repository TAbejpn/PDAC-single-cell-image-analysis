function bbox=checkbbox_v2(bbox,imageH,imageW)

for i=1:size(bbox,1)
    %% edit 2022.04.07
    if bbox(i,1) < 1   %Sx is minus or zero
        bbox(i,3) = bbox(i,3) + bbox(i,1);% add this sentence  #################################
        bbox(i,1) = 1;
    end
    
    if bbox(i,1) > imageW   %Sx > imageW
        bbox(i,1) = imageW -10 ;
    end
    
    if ( bbox(i,1) + bbox(i,3) ) > imageW   %Sx+W > imageW
        
        del= ( bbox(i,1) + bbox(i,3) ) - imageW;
        bbox(i,3) = bbox(i,3) - ( del+2);
        
        if bbox(i,3) < 1
            bbox(i,3)=1;
        end
        
    end
    
        %% edit 2022.04.07
    if bbox(i,2) < 1   %Sy is minus or zero
        bbox(i,4) = bbox(i,4) + bbox(i,2);% add this sentence #################################
        bbox(i,2) = 1;
    end
    
    if bbox(i,2) > imageH   %Sy > imageH
        bbox(i,2) = imageH -10 ;
    end
    
    if ( bbox(i,2) + bbox(i,4) ) > imageH   %Sy+H > imageH
        
        del= ( bbox(i,2) + bbox(i,4) ) - imageH;
        bbox(i,4) = bbox(i,4) - ( del+2);
        
        if bbox(i,4) < 1
            bbox(i,4)=1;
        end
        
    end
    
end

end