function out=param(rect0,lieu)
x0=rect0.x0;y0=rect0.y0;x1=rect0.x1;y1=rect0.y1;
l=num2str(lieu);
switch l
    case'1'
        a=0;b=pi/2;
        A=2/pi*(x1-x0);
        B=x0+i*y0;

    case '2'
        a=pi/2;b=pi;
        A=i*2/pi*(y1-y0);
        B=x1+i*(2*y0-y1);

    case '3'
        a=pi;b=3*pi/2;
        A=-2/pi*(x1-x0);
        B=3*x1-2*x0+i*y1;

    case '4'
        a=3*pi/2;b=2*pi;
        A=-i*2/pi*(y1-y0);
        B=x0+i*(4*y1-3*y0);

end

pa.A=A;pa.B=B;
pa.a=a;pa.b=b;

out = pa;
