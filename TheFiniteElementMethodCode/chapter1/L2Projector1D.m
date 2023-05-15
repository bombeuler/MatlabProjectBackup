function L2Projector1D(xmin,xmax,nlist,func)
figure(1);
hold on
for n = nlist
    h = (xmax-xmin)/n;
    xx = xmin:h:xmax;
    M = gallery("tridiag",n+1,h/6,2*h/3,h/6);
    M(1,1) = h/3;
    M(n+1,n+1) = h/3;
    
    b = (func(xx).*h)';
    b(1) = func(xmin).*h./2;
    b(n+1) = func(xmax).*h./2;
    
    Pf = M\b;
    plot(xx,Pf,xx,func(xx));
end
hold off
end