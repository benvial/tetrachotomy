function out=poly_period(a,b,x)
%k=2
s.P=((a - x).^3.*(a.^2 - 5.*a.*b + 3.*a.*x + 10.*b.^2 - 15.*b.*x + 6.*x.^2))/(a - b).^5;
s.dP=-(30.*(a - x).^2.*(b - x).^2)/(a - b).^5;
out=s;