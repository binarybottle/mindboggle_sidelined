function ZM = reformat_zernike(Z,N)

load coord;

ZM = [];
for c = 1:size(coord,1);
    
    
    if coord(c,1)+1 <= N+1;
        ZM(c,1) = Z(coord(c,1)+1,coord(c,2)+1,coord(c,3)+1);
    else
        break
    end
end
