function out=coupe4(rect0)
x0=rect0.x0;y0=rect0.y0;x1=rect0.x1;y1=rect0.y1;
xm=(x0+x1)/2;
ym=(y0+y1)/2;
rect1.x0=x0;rect1.y0=y0;rect1.x1=xm;rect1.y1=ym;
rect2.x0=xm;rect2.y0=y0;rect2.x1=x1;rect2.y1=ym;
rect3.x0=xm;rect3.y0=ym;rect3.x1=x1;rect3.y1=y1;
rect4.x0=x0;rect4.y0=ym;rect4.x1=xm;rect4.y1=y1;

R(1)=rect1;R(2)=rect2;R(3)=rect3;R(4)=rect4;
out=R;