function h=hermite(x,y,dy,t)
    n=length(x);
    m=length(t);
    L=zeros(n,m);
    dL=zeros(n,m);
    [wT,wX]=meshgrid(t,x);
    wx=prod(wT-wX);
    for i=1:n
        xb=x(i)-x(x~=x(i));
        L(i,:)=wx./(t-x(i))./prod(xb);
        dL(i,:)=sum(1./xb);
    end
    A=(1-2.*dL.*(wT-wX)).*(L.^2);
    B=(wT-wX).*(L.^2);
    h=y*A+dy*B;
end