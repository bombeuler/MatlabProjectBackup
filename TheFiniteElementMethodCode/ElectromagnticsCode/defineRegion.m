% % 
%    |      y_1----------x_2
%    |       |            |
%    |-width-|            |
%    |  PML  |            |
%    |       |            |
%    |      x_1----------y_2
%    |
% gi 是增加的区域

function dl = defineRegion(gi,x1,x2,y1,y2,width,regionStr,regionArray)
    gd = [gi;2 4 x1 x2 x2 x1 y2 y2 y1 y1;...
          2 4 x1-width x2+width x2+width x1-width y2-width y2-width y1+width y1+width]';

        dl = decsg(gd,regionStr,regionArray);
end