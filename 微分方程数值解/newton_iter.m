function X = newton_iter(F,dF,X0,tol)
    Xp=X0;
    Xn=Xp-dF(Xp)\F(Xp);
    while(norm(abs(Xn-Xp),Inf)>=tol)
        Xp=Xn;
        Xn=Xp-dF(Xp)\F(Xp);
    end
    X=Xn;
end