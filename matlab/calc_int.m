
function out=calc_int(rect0,f,tol_conv)

global nf

for lieu=1:4
    clear J
    pa=param(rect0,lieu);

    a=pa.a;b=pa.b;
    A=pa.A;B=pa.B;


    ind=1;nx=2;

    x=a;
    p=poly_period(a,b,x);
    P=p.P;dP=p.dP;
    Pt=(b-a)*P+a;dPt=(b-a)*dP;
    PT=A*Pt+B;
    fper=A*f(PT).*dPt;nf=nf+1;
    Fper=[fper; Pt.*fper ;Pt.^2.*fper];
    M=Fper;
    x=b;
    p=poly_period(a,b,x);
    P=p.P;dP=p.dP;
    Pt=(b-a)*P+a;dPt=(b-a)*dP;
    PT=A*Pt+B;
    fper=A*f(PT).*dPt;nf=nf+1;
    Fper=[fper; PT.*fper; PT.^2.*fper];
    M=[M,Fper];
    hx=(b-a)/(nx-1);
    Iper=hx*sum(M,2);
    J(:,ind)=Iper;


    ind=2;nx=3;
    x=(a+b)/2;
    p=poly_period(a,b,x);
    P=p.P;dP=p.dP;
    Pt=(b-a)*P+a;dPt=(b-a)*dP;
    PT=A*Pt+B;
    fper=A*f(PT).*dPt;nf=nf+1;
    Fper=[fper; PT.*fper ;PT.^2.*fper];
    M=[M,Fper];
    hx=(b-a)/(nx-1);
    Iper=hx*sum(M,2);
    J(:,ind)=Iper;

    condition=true;
    while condition;
        condition1=  ((abs(real(J(1,ind)-J(1,ind-1)))/(abs(real(J(1,ind)))))<tol_conv ) && ((abs(imag(J(1,ind)-J(1,ind-1)))/(abs(imag(J(1,ind)))))<tol_conv);
        condition2=  ((abs(real(J(2,ind)-J(2,ind-1)))/(abs(real(J(2,ind)))))<tol_conv) && ((abs(imag(J(2,ind)-J(2,ind-1)))/(abs(imag(J(2,ind)))))<tol_conv);
        condition3=  ((abs(real(J(3,ind)-J(3,ind-1)))/(abs(real(J(3,ind)))))<tol_conv ) && ((abs(imag(J(3,ind)-J(3,ind-1)))/(abs(imag(J(3,ind)))))<tol_conv);
        condition=~(condition1 && condition2 && condition3);

        if ind>50
          break
          'break'
        end
        ind=ind+1;
        nx=2^(ind-1)+1;
        for j=1:2^(ind-2)
            x=a+(2*j-1)*(b-a)/2^(ind-1);
            p=poly_period(a,b,x);
            P=p.P;dP=p.dP;
            Pt=(b-a)*P+a;dPt=(b-a)*dP;
            PT=A*Pt+B;
            fper=A*f(PT).*dPt;nf=nf+1;
            Fper=[fper; PT.*fper; PT.^2.*fper];
            M=[M,Fper];
        end
        hx=(b-a)/(nx-1);
        Iper=hx*sum(M,2);
        J(:,ind)=Iper;
    end
    nx_opt=nx;
    Int(:,lieu)=J(:,end);

end

INT=1/(2*i*pi)*sum(Int,2);
s.I0=INT(1);s.I1=INT(2);s.I2=INT(3);

out=s;
