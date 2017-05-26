function npole=pole_2D(f,rect0,tolre0,tolim0,tolre1,tolim1,tolreR,tolimR,tol_conv,trace) ;
global pol res n dist
npole=0 ;
q=ispole(rect0,f,tolre0,tolim0,tolre1,tolim1,tolreR,tolimR,tol_conv);
n=n+1
if trace
    trace_rect(rect0,'none'); 
    hold on;
end
npole=q.r2;
if npole
    R=coupe4(rect0);
    for i=1:4
        npole=pole_2D(f,R(i),tolre0,tolim0,tolre1,tolim1,tolreR,tolimR,tol_conv,trace) ;
    end;
else
    pol=[pol;q.P];res=[res;q.A];
    d_min=min([abs(real(q.P)-rect0.x1),abs(real(q.P)-rect0.x0),abs(imag(q.P)-rect0.y1),abs(imag(q.P)-rect0.y0)]);
    dist=[dist,d_min];
    if q.r1
        fprintf('found : %e + %ei\n',real(q.P),imag(q.P))
    end
end

end
