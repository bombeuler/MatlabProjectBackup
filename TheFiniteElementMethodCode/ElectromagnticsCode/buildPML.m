function pml = buildPML(xlr,xrl,ybt,ytb,ll,SIGMA0,n)
    function SIGMA = pmlfunc(x,y)
        SIGMA = 0;
        if (x >= xlr-ll && x <= xlr)
           SIGMA = SIGMA + SIGMA0 * ((xlr - x)/ll)^n;
        elseif (x >= xrl && x <= xrl+ll)
            SIGMA = SIGMA + SIGMA0 * ((x - xrl)/ll)^n;
        end

        if (y >= ybt-ll && y <= ybt)
           SIGMA = SIGMA + SIGMA0 * ((ybt - y)/ll)^n;
        elseif (y >= ytb && y <= ytb+ll)
            SIGMA = SIGMA + SIGMA0 * ((y - ytb)/ll)^n;
        end

    end

pml = @pmlfunc;

end