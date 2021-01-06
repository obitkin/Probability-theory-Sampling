%Lab5
%{
sizes = [20 60 100];
p = [0 0.5 0.9];
mu = [0 0];
disp("Script begin");
%{
R = mvnrnd(mu, [1 0; 0 1], 20);
figure;
plot(R(:,1),R(:,2),'+');
ylim([-4 4])
xlim([-4 4])
hold on;
f=ezplot('x^2 + y^2 =4'); 
set(f, 'LineWidth', 2)
title("");

R = mvnrnd(mu, [1 0.5; 0.5 1], 20);
figure;
plot(R(:,1),R(:,2),'+');
ylim([-4 4])
xlim([-4 4])
hold on;
f=ezplot('x^2 -x*y + y^2 = 4'); 
set(f, 'LineWidth', 2)
title("");

R = mvnrnd(mu, [1 0.9; 0.9 1], 20);
figure;
plot(R(:,1),R(:,2),'+');
ylim([-5 5])
xlim([-5 5])
hold on;
f=ezplot('x^2 - 1.8*x*y + y^2 = 4'); 
set(f, 'LineWidth', 2)
title("");

R = mvnrnd(mu, [1 0; 0 1], 60);
figure;
plot(R(:,1),R(:,2),'+');
ylim([-4 4])
xlim([-4 4])
hold on;
f=ezplot('x^2 + y^2 = 4'); 
set(f, 'LineWidth', 2)
title("");

R = mvnrnd(mu, [1 0.5; 0.5 1], 60);
figure;
plot(R(:,1),R(:,2),'+');
ylim([-4 4])
xlim([-4 4])
hold on;
f=ezplot('x^2 -x*y + y^2 = 4'); 
set(f, 'LineWidth', 2)
title("");

R = mvnrnd(mu, [1 0.9; 0.9 1], 60);
figure;
plot(R(:,1),R(:,2),'+');
ylim([-5 5])
xlim([-5 5])
hold on;
f=ezplot('x^2 - 1.8*x*y + y^2 = 4'); 
set(f, 'LineWidth', 2)
title("");
%}

R = mvnrnd(mu, [1 0; 0 1], 100);
figure;
plot(R(:,1),R(:,2),'+');
ylim([-4 4])
xlim([-4 4])
hold on;
f=ezplot('x^2 + y^2 = 4'); 
set(f, 'LineWidth', 2)
title("");

R = mvnrnd(mu, [1 0.5; 0.5 1], 100);
figure;
plot(R(:,1),R(:,2),'+');
ylim([-4 4])
xlim([-4 4])
hold on;
f=ezplot('x^2 -x*y + y^2 =4'); 
set(f, 'LineWidth', 2)
title("");

R = mvnrnd(mu, [1 0.9; 0.9 1], 100);
figure;
plot(R(:,1),R(:,2),'+');
ylim([-5 5])
xlim([-5 5])
hold on;
f=ezplot('x^2 - 1.8*x*y + y^2 = 4'); 
set(f, 'LineWidth', 2)
title("");
%}
%{
for n = 1:3
    for i = 1:3
        R = [];
        R_q = [];
        R_s = [];
        for j = 1:1000
            v = mvnrnd(mu, [1 p(i); p(i) 1], sizes(n))';
            R = [R r(v)];
            R_q = [R_q r_Q(v)];   
            R_s = [R_s r_S(v)];       
        end
        disp("----------------");
        disp("Size:  " + sizes(n) + "    P:  " +p(i));
        disp("    " + " R " + "       R_s     " + "R_q");
        disp("E(z)" + " & " + E(R) + " & " + E(R_s) + " & " + E(R_q));
        disp("E(z^2)" + " & " + E(R.^2) + " & " + E(R_s.^2) + " & " + E(R_q.^2));
        disp("(z)" + " & " + D(R) + " & " + D(R_s) + " & " + D(R_q));
    end
end


for n = 1:3
        R = [];
        R_q = [];
        R_s = [];
        for j = 1:1000
            v = MvNormRnd(sizes(n));
            R = [R r(v)];
            R_q = [R_q r_Q(v)];   
            R_s = [R_s r_S(v)];       
        end
        disp("----------------");
        disp("Size:  " + sizes(n) + "    P:  " +p(i));
        disp("    " + " R " + "       R_s     " + "R_q");
        disp("E(z)" + " & " + E(R) + " & " + E(R_s) + " & " + E(R_q));
        disp("E(z^2)" + " & " + E(R.^2) + " & " + E(R_s.^2) + " & " + E(R_q.^2));
        disp("(z)" + " & " + D(R) + " & " + D(R_s) + " & " + D(R_q));
end
%}

%Lab6
y_orig = @(c) (2+2*c);
%Выборка без возмущений
x = -1.8:0.2:2;
e = normrnd(0,1,[1,20]);
y = 2 + 2*x + e;
%МНК
xx = x.*x;
yy = y.*y;
xy = x.*y;
B1 = (E(xy)-E(x)*E(y))/(E(xx)-(E(x))^2);
B0 = E(y) - E(x)*B1;
disp("МНК без возмущений");
disp("B0= " + B0 + "  B1= " + B1);
%МНМ
errfcn = @(b) sum(abs(y-b(1)-b(2)*x));
x0 = [0,0];
B = fminsearch(errfcn,x0);
disp("МНМ без возмущений");
disp("B0= " + B(1) + "  B1= " + B(2));
%График
y_LS = @(x) B0 + B1*x; %МНК
y_LAD = @(x) B(1) + B(2)*x; %МНМ

figure;
x = linspace(-1.8,2);
plot(x,y_orig(x),'DisplayName','Модель','LineWidth', 2)
hold on
plot(x,y_LS(x),'DisplayName','МНК','LineWidth', 2)
plot(x,y_LAD(x),'DisplayName','МНМ','LineWidth', 2)
plot(-1.8:0.2:2,y,'o','DisplayName','Выборка','LineWidth', 2)
hold off
xlim([-2, 2.2]);
legend('Location','northwest');

%Выборка c возмущением
x = -1.6:0.2:1.8;
e = normrnd(0,1,[1,18]);
y = 2 + 2*x + e;
x = [-1.8 x 2];
y = [2+2*(-1.8)+normrnd(0,1,[1,1])+10 y 2+2*(2)+normrnd(0,1,[1,1])-10];
%МНК
xx = x.*x;
yy = y.*y;
xy = x.*y;
B1 = (E(xy)-E(x)*E(y))/(E(xx)-(E(x))^2);
B0 = E(y) - E(x)*B1;
disp("МНК с возмущением");
disp("B0= " + B0 + "  B1= " + B1);
%МНМ
errfcn = @(b) sum(abs(y-b(1)-b(2)*x));
x0 = [0,0];
B = fminsearch(errfcn,x0);
disp("МНМ с возмущением");
disp("B0= " + B(1) + "  B1= " + B(2));
%График
y_LS = @(x) B0 + B1*x; %МНК
y_LAD = @(x) B(1) + B(2)*x; %МНМ

figure;
x = linspace(-1.8,2);
plot(x,y_orig(x),'DisplayName','Модель','LineWidth', 2)
hold on
plot(x,y_LS(x),'DisplayName','МНК','LineWidth', 2)
plot(x,y_LAD(x),'DisplayName','МНМ','LineWidth', 2)
plot(-1.8:0.2:2,y,'o','DisplayName','Выборка','LineWidth', 2)
hold off
xlim([-2, 2.2]);
legend('Location','northwest');

function Average = E(arr)
sz = size(arr);
sz = sz(2);
Average = sum(arr)/sz;
end

function dispersya = D(arr)
sz = size(arr);
sz = sz(2);
srednee = E(arr);
summa = 0;
for i = 1:sz
    summa = summa + (arr(i)-srednee)^2;
end
dispersya = summa / sz;
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

function result = r(arr2D) % arr 2*100
n = size(arr2D);
n = n(2);
X = arr2D(1,:);
Y = arr2D(2,:);
x_srendee = E(X);
y_srendee = E(Y);
res = 0;
for i = 1:n
    res = res + (X(i)-x_srendee)*(Y(i)-y_srendee);
end
res = res / n;
result = res / sqrt(D(X)*D(Y));
end

function result = r_Q(arr2D)
N = [0 0 0 0];
n = size(arr2D);
n = n(2);
X = arr2D(1,:);
Y = arr2D(2,:);
medX = median(X);
medY = median(Y);
for i = 1:n
    if (X(i) >= medX) 
        if (Y(i) >= medY)
            N(1) = N(1)+1;
        else
            N(4) = N(4)+1;
        end
    else 
        if (Y(i) >= medY)
            N(2) = N(2)+1;
        else
            N(3) = N(3)+1;
        end   
    end
end
result = ((N(1)+N(3))-(N(2)+N(4)))/n;
end

function result = r_S(arr2D)
n = size(arr2D);
n = n(2);
X = arr2D(1,:);
Y = arr2D(2,:);
u = (n+1)/2;
res = 0;
RangX = zeros(1,n);
RangY = zeros(1,n);
for i = 1:n
    RangX(find(X==min(X))) = i;
    X(find(X==min(X))) = 1000;
    RangY(find(Y==min(Y))) = i;
    Y(find(Y==min(Y))) = 1000;
end
for i = 1:n
    res = res + (RangX(i)-u)*(RangY(i)-u);
end
res = res / n;
result = (res / (sqrt(D(RangX)*D(RangY))));
end

function result = r_qq(x,y)
medY = median(y);
medX = median(x);
end

function vyborka = MvNormRnd(size) %2d 
i = 0;
res = [];
max = 0.3302;
while i < size
    x = 5*[rand rand]-2.5;
    P = NormPdf(x) / max;
    if rand < P
        res = [res x'];
        i = i + 1;
    end   
end
vyborka = res;
%{
P1 = [ 1 0.9; 0.9 1];
P2 = [ sqrt(10) -0.9; -0.9 sqrt(10)];
MU = [0 0];
res = [];
for i = 1:size
    if (rand < 0.9)
        res = [res mvnrnd(MU, P1, 1)'];
    else
        res = [res mvnrnd(MU, P2, 1)'];
    end
end
vyborka = res;
%}
end

function probability = NormPdf(X) %2d 
p1 = [ 1 0.9; 0.9 1];
p2 = [ (10) 0.9; 0.9 (10)];
mu = [0 0];
probability = 0.95*mvnpdf(X,mu,p1)+0.05*mvnpdf(X,mu,p2);
end
