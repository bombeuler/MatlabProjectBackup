jj_len = 5;
primeN = 4;
disp(["max              ","l2"]);

load("u512.mat","u512");

disp("----------------------------------------------------------");
for jj = 1:jj_len
    N = primeN*2^jj;
    fprintf("N=%d\n",N);
    uu = deltaProcess(N);
    errorMax = errorNorm(uu,u512,N,-1);
    errorL2 = errorNorm(uu,u512,N,0);
    disp([errorMax,errorL2]);
end