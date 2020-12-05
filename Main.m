size10 = 10;
size50 = 50;
size1000 = 1000; 
odin = 1/(2^(0.5));
x = -4:0.1:4;

norm20 = normrnd(0,1,[1,20]);
norm60 = normrnd(0,1,[1,60]);
norm100 = normrnd(0,1,[1,100]);
x4 = -4:0.01:4;

vyb = norm20;
figure;
cdfplot(vyb);
hold on;
plot(x4,normcdf(x4,0,1));
hold on;
line([max(vyb),4],[1,1]);
line([-4,min(vyb)],[0,0]);
ylim([-0.1 1.1])
xlim([-4 4])
title('Функция нормального распределения N = 20');

vyb = norm60;
figure;
cdfplot(vyb);
hold on;
plot(x4,normcdf(x4,0,1));
hold on;
line([max(vyb),4],[1,1]);
line([-4,min(vyb)],[0,0]);
ylim([-0.1 1.1])
xlim([-4 4])
title('Функция нормального распределения N = 60');

vyb = norm100;
figure;
cdfplot(vyb);
hold on;
plot(x4,normcdf(x4,0,1));
hold on;
line([max(vyb),4],[1,1]);
line([-4,min(vyb)],[0,0]);
ylim([-0.1 1.1])
xlim([-4 4])
title('Функция нормального распределения N = 100');

cau20 = trnd(1,1,20);
cau60 = trnd(1,1,60);
cau100 = trnd(1,1,100);

vyb = cau20;
figure;
cdfplot(vyb);
hold on;
plot(x4,tcdf(x4,1));
hold on;
line([max(vyb),4],[1,1]);
line([-4,min(vyb)],[0,0]);
ylim([-0.1 1.1])
xlim([-4 4])
title('Функция распределения Коши N = 20');

vyb = cau60;
figure;
cdfplot(vyb);
hold on;
plot(x4,tcdf(x4,1));
hold on;
ylim([-0.1 1.1])
xlim([-4 4])
title('Функция распределения Коши N = 100');

vyb = cau100;
figure;
cdfplot(vyb);
hold on;
plot(x4,tcdf(x4,1));
hold on;
ylim([-0.1 1.1])
xlim([-4 4])
title('Функция распределения Коши N = 100');

lap20 = laprnd(20);
lap60 = laprnd(60);
lap100 = laprnd(100);

vyb = lap20;
figure;
cdfplot(vyb);
hold on;
plot(x4,lapcdf(x4,odin,0));
hold on;
line([max(vyb),4],[1,1]);
line([-4,min(vyb)],[0,0]);
ylim([-0.1 1.1])
xlim([-4 4])
title('Функция распределения Лапласа N = 20');

vyb = lap60;
figure;
cdfplot(vyb);
hold on;
plot(x4,lapcdf(x4,odin,0));
hold on;
line([max(vyb),4],[1,1]);
line([-4,min(vyb)],[0,0]);
ylim([-0.1 1.1])
xlim([-4 4])
title('Функция распределения Лапласа N = 60');

vyb = lap100;
figure;
cdfplot(vyb);
hold on;
plot(x4,lapcdf(x4,odin,0));
hold on;
line([max(vyb),4],[1,1]);
line([-4,min(vyb)],[0,0]);
ylim([-0.1 1.1])
xlim([-4 4])
title('Функция распределения Лапласа N = 100');

poiss20 = poissrnd(10,1,20);
poiss60 = poissrnd(10,1,60);
poiss100 = poissrnd(10,1,100);

vyb = poiss20;
figure;
cdfplot(vyb);
hold on;
plot(x4,poisscdf(x4,10));
hold on;
line([max(vyb),4],[1,1]);
line([-4,min(vyb)],[0,0]);
ylim([-0.1 1.1])
xlim([-4 4])
title('Функция распределения Пуассона N = 20');

%{
NormX1 = -2.698;
NormX2 = 2.698;
CaushyX1 = -4;
CaushyX2 = 4;
LaplasQ3 = sqrt(2)*log(2);
LaplasQ1 = -sqrt(2)*log(2);
LaplasX1 = LaplasQ1-1.5*(LaplasQ3-LaplasQ1);
LaplasX2 = LaplasQ3+1.5*(LaplasQ3-LaplasQ1);
PuassX1 = 2;
PuassX2 = 18;
UnifX1 = -3.464;
UnifX2 = 3.464;

percentOutline20 = [];
percentOutline100 = [];
for i = 1:1000
    percentOutline20 = [percentOutline20 percentOfOutLine(normrnd(0,1,[1,20]),NormX1,NormX2)];
    percentOutline100 = [percentOutline100 percentOfOutLine(normrnd(0,1,[1,100]),NormX1,NormX2)];   
end
disp("Нормальное")
disp(average(percentOutline20));
disp(average(percentOutline100));

percentOutline20 = [];
percentOutline100 = [];
for i = 1:1000
    percentOutline20 = [percentOutline20 percentOfOutLine(trnd(1,1,20),CaushyX1,CaushyX2)];
    percentOutline100 = [percentOutline100 percentOfOutLine(trnd(1,1,100),CaushyX1,CaushyX2)];   
end
disp("Коши")
disp(average(percentOutline20));
disp(average(percentOutline100));


percentOutline20 = [];
percentOutline100 = [];
for i = 1:1000
    percentOutline20 = [percentOutline20 percentOfOutLine(laprnd(20),LaplasX1,LaplasX2)];
    percentOutline100 = [percentOutline100 percentOfOutLine(laprnd(100),LaplasX1,LaplasX2)];   
end
disp("Лапласа")
disp(average(percentOutline20));
disp(average(percentOutline100));

percentOutline20 = [];
percentOutline100 = [];
for i = 1:1000
    percentOutline20 = [percentOutline20 percentOfOutLine(poissrnd(10,1,20),PuassX1,PuassX2)];
    percentOutline100 = [percentOutline100 percentOfOutLine(poissrnd(10,1,100),PuassX1,PuassX2)];   
end
disp("Пуассон")
disp(average(percentOutline20));
disp(average(percentOutline100));

percentOutline20 = [];
percentOutline100 = [];
for i = 1:1000
    percentOutline20 = [percentOutline20 percentOfOutLine(unifrnd(-sqrt(3),sqrt(3),1,20),UnifX1,UnifX2)];
    percentOutline100 = [percentOutline100 percentOfOutLine(unifrnd(-sqrt(3),sqrt(3),1,100),UnifX1,UnifX2)];   
end
disp("Равномерное")
disp(average(percentOutline20));
disp(average(percentOutline100));


g1 = repmat({'100'},100,1);
g2 = repmat({'20'},20,1);
grp = [g1; g2];
%rng default  % For reproducibility
n20 = normrnd(0,1,[1,20]);
n100 = normrnd(0,1,[1,100]);
boxNorm = [n100 n20];
figure;
boxplot(boxNorm,grp);
title("Нормальное")

c20 = trnd(1,1,20);
c100 = trnd(1,1,100);
boxCaushy = [c100 c20];
figure;
boxplot(boxCaushy,grp);
title("Коши")

l20 = laprnd(20);
l100 = laprnd(100);
boxLap = [l100 l20];
figure;
boxplot(boxLap,grp);
title("Лаплас")

p20 = poissrnd(10,1,20);
p100 = poissrnd(10,1,100);
boxPoiss = [p100 p20];
figure;
boxplot(boxPoiss,grp);
title("Пуассон")

u20 = unifrnd(-sqrt(3),sqrt(3),1,20);
u100 = unifrnd(-sqrt(3),sqrt(3),1,100);
boxUnif = [u100 u20];
figure;
boxplot(boxUnif,grp);
title("Равномерное")
%}
%{
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
        vyborka = sort(unifrnd(-sqrt(3),sqrt(3),1,array(i)));
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
%}
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

r1 = unifrnd(-sqrt(3),sqrt(3),1,size10);
r2 = unifrnd(-sqrt(3),sqrt(3),1,size50);
r3 = unifrnd(-sqrt(3),sqrt(3),1,size1000);

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

function znach = lapcdf(x, a, b)
sz = size(x);
sz = sz(2);
res = [];
for i = 1:sz
    if (x(i) < b)
        res = [res (1/2)*exp(a*(x(i)-b))];
    else
        res = [res 1-(1/2)*exp((-a)*(x(i)-b))];
    end
end
znach = res;
end

function res = laprnd(size)
odin2 = 1/(2^(0.5));
i = 1;
result = [];
arrOfRand = rand(1,size+1);
while i < (size + 1)
    result = [result (1/odin2)*log((arrOfRand(i)/(arrOfRand(i+1))))];
    i = i + 1;
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

function result = percentOfOutLine(arr,X1,X2)
sz = size(arr);
sz = sz(2);

numberOfOutLine = 0;
for i = 1:sz
    if (arr(i) < X1)
        numberOfOutLine = numberOfOutLine + 1;
    end
    if (arr(i) > X2)
        numberOfOutLine = numberOfOutLine + 1;
    end
end
result = numberOfOutLine / sz;
end