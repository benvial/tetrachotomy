

function out=ispole(rect0,f,tolre0,tolim0,tolre1,tolim1,tolreR,tolimR,tol_conv)

    s=calc_int(rect0,f,tol_conv);

    I0=s.I0;
    I1=s.I1;
    I2=s.I2;
    R21=I2/I1;
    R10=I1/I0;
    r0=false;r1=false;r2=false;
    P=[];A=[];



    if ((abs(real(I0))<tolre0) && (abs(imag(I0))<tolim0)) && ((abs(real(I1))<tolre1) && (abs(imag(I1))<tolim1))
        message='no poles';
        r0=true;

    else if (abs(real(R21-R10)/real(R10))<tolreR) && (abs(imag(R21-R10)/imag(R10))<tolimR)
            message='one pole';
            P=R21;
            A=I0;
            r1=true;
        else

            message='several poles';
            r2=true;
        end
    end

p.P=P;p.A=A;
p.r0=r0;p.r1=r1;p.r2=r2;
p.message=message;
out=p;
