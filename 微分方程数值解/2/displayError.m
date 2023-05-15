jj_len = 5;
primeN = 4;
for jj = 1:jj_len
    N = primeN*2^jj;
    fprintf("N=%d\n",N);
    uu1 = t2_3(N);
    uu2 = t2_3_remake(N);
    disp([errorNorm(uu1,uu2,1/N,-1,0) errorNorm(uu1,uu2,1/N,0,0)]);
end

disp(["max              ","l2"]);
disp("-------------------------remake--------------------------");
for jj = 1:jj_len
    N = primeN*2^jj;
    fprintf("N=%d\n",N);
    t2_3_remake(N);
end

disp("---------------------------origin------------------------");
for jj = 1:jj_len
    N = primeN*2^jj;
    fprintf("N=%d\n",N);
    t2_3(N);
end