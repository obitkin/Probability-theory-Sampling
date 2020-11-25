size10 = 10;
size50 = 50;
size1000 = 1000; 
odin = 1/(2^(0.5));
x = -4:0.1:4;

array = [10 100 1000];

for i = 1:3 
    Aver = [];
    Med = [];
    Z_r = [];
    Z_q = [];
    Z_tr = [];
    for j = 1:1000
        vyborka = sort(normrnd(0,1,[1,array(i)]));
        Aver = [Aver average(vyborka)];
        Med = [Med median(vyborka)];
        Z_r = [Z_r halfExtremumSum(vyborka)];
        Z_q = [Z_q halfQuantileSum(vyborka)];
        Z_tr = [Z_tr truncatedMean(vyborka)];
    end
    srendee = [average(Aver) average(Med) average(Z_r) average(Z_q) average(Z_tr)];
    D = [findD(Aver) findD(Med) findD(Z_r) findD(Z_q) findD(Z_tr)];
    disp("Norm: size = " + array(i));
    disp("    Average   Median     Z_R       Z_Q       Z_tr");
    for k = 1:5
        disp(srendee(k))
    end
    for k = 1:5
        disp(D(k))
    end
end

for i = 1:3 
    Aver = [];
    Med = [];
    Z_r = [];
    Z_q = [];
    Z_tr = [];
    for j = 1:1000
        vyborka = sort(trnd(1,1,array(i)));
        Aver = [Aver average(vyborka)];
        Med = [Med median(vyborka)];
        Z_r = [Z_r halfExtremumSum(vyborka)];
        Z_q = [Z_q halfQuantileSum(vyborka)];
        Z_tr = [Z_tr truncatedMean(vyborka)];
    end
    srendee = [average(Aver) average(Med) average(Z_r) average(Z_q) average(Z_tr)];
    D = [findD(Aver) findD(Med) findD(Z_r) findD(Z_q) findD(Z_tr)];
    disp("Caushy: size = " + array(i));
    disp("    Average   Median     Z_R       Z_Q       Z_tr");
    for k = 1:5
        disp(srendee(k))
    end
    disp("---------------------------");
    for k = 1:5
        disp(D(k))
    end
end

for i = 1:3 
    Aver = [];
    Med = [];
    Z_r = [];
    Z_q = [];
    Z_tr = [];
    for j = 1:1000
        vyborka = sort(laprnd(array(i)));
        Aver = [Aver average(vyborka)];
        Med = [Med median(vyborka)];
        Z_r = [Z_r halfExtremumSum(vyborka)];
        Z_q = [Z_q halfQuantileSum(vyborka)];
        Z_tr = [Z_tr truncatedMean(vyborka)];
    end
    srendee = [average(Aver) average(Med) average(Z_r) average(Z_q) average(Z_tr)];
    D = [findD(Aver) findD(Med) findD(Z_r) findD(Z_q) findD(Z_tr)];
    disp("Laplace: size = " + array(i));
    disp("    Average   Median     Z_R       Z_Q       Z_tr");
    for k = 1:5
        disp(srendee(k))
    end
    disp("---------------------------");
    for k = 1:5
        disp(D(k))
    end
end


for i = 1:3 
    Aver = [];
    Med = [];
    Z_r = [];
    Z_q = [];
    Z_tr = [];
    for j = 1:1000
        vyborka = sort(poissrnd(10,1,array(i)));
        Aver = [Aver average(vyborka)];
        Med = [Med median(vyborka)];
        Z_r = [Z_r halfExtremumSum(vyborka)];
        Z_q = [Z_q halfQuantileSum(vyborka)];
        Z_tr = [Z_tr truncatedMean(vyborka)];
    end
    srendee = [average(Aver) average(Med) average(Z_r) average(Z_q) average(Z_tr)];
    D = [findD(Aver) findD(Med) findD(Z_r) findD(Z_q) findD(Z_tr)];
    disp("Poiss: size = " + array(i));
    disp("    Average   Median     Z_R       Z_Q       Z_tr");
    for k = 1:5
        disp(srendee(k))
    end
    disp("---------------------------");
    for k = 1:5
        disp(D(k))
    end
end

for i = 1:3 
    Aver = [];
    Med = [];
    Z_r = [];
    Z_q = [];
    Z_tr = [];
    for j = 1:1000
        vyborka = sort(unifrnd(-sqrt(3),sqrt(3),array(i)));
        Aver = [Aver average(vyborka)];
        Med = [Med median(vyborka)];
        Z_r = [Z_r halfExtremumSum(vyborka)];
        Z_q = [Z_q halfQuantileSum(vyborka)];
        Z_tr = [Z_tr truncatedMean(vyborka)];
    end
    srendee = [average(Aver) average(Med) average(Z_r) average(Z_q) average(Z_tr)];
    D = [findD(Aver) findD(Med) findD(Z_r) findD(Z_q) findD(Z_tr)];
    disp("Unif: size = " + array(i));
    disp("    Average   Median     Z_R       Z_Q       Z_tr");
    for k = 1:5
        disp(srendee(k))
    end
    disp("---------------------------");
    for k = 1:5
        disp(D(k))
    end
end
%{

n1 = normrnd(0,1,[1,size10]);
n2 = normrnd(0,1,[1,size50]);
n3 = normrnd(0,1,[1,size1000]);

figure;
histogram(n1,10,'Normalization','pdf'); grid on; ylabel('density');
hold on;
plot(x,normpdf(x));

figure;
histogram(n2,15,'Normalization','pdf'); grid on; ylabel('density');
hold on;
plot(x,normpdf(x));

figure;
histogram(n3,30,'Normalization','pdf'); grid on; ylabel('density');
hold on;
plot(x,normpdf(x));

x2 = -20:0.1:20;
c1 = trnd(1,1,size10);
c2 = trnd(1,1,size50);
c3 = trnd(1,1,size1000);

figure;
histogram(c1,10,'Normalization','pdf'); grid on; ylabel('density');
hold on;
plot(x2,tpdf(x2,1));

figure;
histogram(c2,15,'Normalization','pdf'); grid on; ylabel('density');
hold on;
plot(x2,tpdf(x2,1));

figure;
histogram(c3,30,'Normalization','pdf'); grid on; ylabel('Log(density)');
hold on;
x33 = min(c3):0.1:max(c3);
plot(x33,tpdf(x33,1));
set(gca,'YScale','log')
x3 = -6:0.1:6;
l1 = laprnd(size10);
l2 = laprnd(size50);
l3 = laprnd(size1000);

figure;
histogram(l1,10,'Normalization','pdf'); grid on; ylabel('density');
hold on;
plot(x3,lappdf(x3,odin,0));

figure;
histogram(l2,20,'Normalization','pdf'); grid on; ylabel('density');
hold on;
plot(x3,lappdf(x3,odin,0));

figure;
histogram(l3,30,'Normalization','pdf'); grid on; ylabel('density');
hold on;
plot(x3,lappdf(x3,odin,0));

x4 = 0:1:20;

p1 = poissrnd(10,1,size10);
p2 = poissrnd(10,1,size50);
p3 = poissrnd(10,1,size1000);

edges = 0:1:20; %for integer bin 

figure;
histogram(p1,edges,'Normalization','pdf'); grid on; ylabel('density');
hold on;
plot(x4,poisspdf(x4,10));

figure;
histogram(p2,edges,'Normalization','pdf'); grid on; ylabel('density');
hold on;
plot(x4,poisspdf(x4,10));

figure;
histogram(p3,edges,'Normalization','pdf'); grid on; ylabel('density');
hold on;
plot(x4,poisspdf(x4,10));

x5 = -3:0.05:3;

r1 = unifrnd(-sqrt(3),sqrt(3),size10);
r2 = unifrnd(-sqrt(3),sqrt(3),size50);
r3 = unifrnd(-sqrt(3),sqrt(3),size1000);

figure;
histogram(r1,10,'Normalization','pdf'); grid on; ylabel('density');
hold on;
plot(x5,unifpdf(x5,-sqrt(3),sqrt(3)));

figure;
histogram(r2,10,'Normalization','pdf'); grid on; ylabel('density');
hold on;
plot(x5,unifpdf(x5,-sqrt(3),sqrt(3)));

figure;
histogram(r3,30,'Normalization','pdf'); grid on; ylabel('density');
hold on;
plot(x5,unifpdf(x5,-sqrt(3),sqrt(3)));
%}

function znach = lappdf(x, a, b)
    znach = (a/2)*exp(-a*abs(x-b));
end

function res = laprnd(size)
odin2 = 1/(2^(0.5));
i = 0;
result = [];
while i < size
    %P = rand * odin2 * 0.5;
    member = (rand * 50)-25; %-log(P*2*sqrt(2))*sqrt(2);
    if (lappdf(member,odin2,0) / odin2 * 0.5) >= rand
        i=i+1;
        if rand > 0.5
            member = -member;
        end
        result = [result member];
    end
end
res = result;
end

function SampleMeanAndCovariance = average(arr)
sz = size(arr);
sz = sz(2);
SampleMeanAndCovariance = sum(arr)/sz;
end

function res1 = median(arr)
sz = size(arr);
sz = sz(2);
if mod(sz, 2) == 0
    resultat = (arr(sz/2) + arr(sz/2+1)) / 2;
else
    resultat = arr((sz+1)/2);
end
res1 = resultat;
end

function result = halfExtremumSum(arr)
sz = size(arr);
sz = sz(2);
result = (arr(1) + arr(sz))/2;
end

function result = halfQuantileSum(arr)
first = 1/4;
second = 3/4;
sz = size(arr);
sz = sz(2);
indexFirst = findIndex(first,sz);
indexSecond = findIndex(second,sz);
result = (arr(indexFirst) + arr(indexSecond))/2;
end

function result = findIndex(quantile,size)
if mod(size*quantile,1) == 0
    result = size*quantile;
else
    result = fix(size*quantile)+1;
end
end

function result = truncatedMean(arr)
sz = size(arr);
sz = sz(2);
r = ceil(sz/4);
summa = 0;
for i = r+1:sz-r
    summa = summa + arr(i);
end
result = summa / (sz-2*r);
end

function result = findD(arr)
sz = size(arr);
sz = sz(2);
srednee = average(arr);
summa = 0;
for i = 1:sz
    summa = summa + (arr(i)-srednee)^2;
end
result = summa / sz;
end