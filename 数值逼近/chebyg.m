function p=chebyg(f,n,t)
   x=cos( (pi/(n)) .*(0:n)); 
  wy=((-1).^(0:n)).*f(x).*[1/2,ones(1,n-1),1/2]; %wk*f(xk)
   [wT,wX]=meshgrid(t,x);
   TX=wT-wX;
   wx=(2^(n-1)/n).*prod(TX);
   WF=1./(TX);
   p=wx.*(wy*WF);
end