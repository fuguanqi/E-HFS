function [is_intersect,x,y]=intersect_at(x1,y1,x2,y2,x3,y3,x4,y4)

if (y3-y4)*(x1-x2)-(x3-x4)*(y1-y2)==0
    is_intersect=0;
    x=-1;
    y=-1;
else
    a1=((y2-y4)*(x3-x4)-(x2-x4)*(y3-y4))/((y3-y4)*(x1-x2)-(x3-x4)*(y1-y2));

    a2=((y2-y4)*(x1-x2)-(x2-x4)*(y1-y2))/((y3-y4)*(x1-x2)-(x3-x4)*(y1-y2));
    if a1>=0 && a1<=1 && a2>=0 && a2<=1
        is_intersect=1;
    else
        is_intersect=0;
    end

    x=a1*x1+(1-a1)*x2;
    y=a1*y1+(1-a1)*y2;

end

end