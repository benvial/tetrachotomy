function trace_rect(rect0,col)
x0=rect0.x0;y0=rect0.y0;x1=rect0.x1;y1=rect0.y1;
rectangle('Position',[x0,y0,x1-x0,y1-y0],'Facecolor',col)
