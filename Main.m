size10 = 10;
size50 = 50;
size1000 = 1000; 
odin = 1/(2^(0.5));
x = -4:0.1:4;
pppp = rand * rand - rand + rand(12);
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

p1 = poissrnd(10,size10);
p2 = poissrnd(10,size50);
p3 = poissrnd(10,size1000);

figure;
histogram(p1,10,'Normalization','pdf'); grid on; ylabel('density');
hold on;
plot(x4,poisspdf(x4,10));

figure;
histogram(p2,20,'Normalization','pdf'); grid on; ylabel('density');
hold on;
plot(x4,poisspdf(x4,10));

figure;
histogram(p3,30,'Normalization','pdf'); grid on; ylabel('density');
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
histogram(r2,20,'Normalization','pdf'); grid on; ylabel('density');
hold on;
plot(x5,unifpdf(x5,-sqrt(3),sqrt(3)));

figure;
histogram(r3,30,'Normalization','pdf'); grid on; ylabel('density');
hold on;
plot(x5,unifpdf(x5,-sqrt(3),sqrt(3)));

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
