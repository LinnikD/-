%% Планирование испытаний
% Первая выборка X1: НР(t < Tотказ) + ФР(t > Tотказ)
% Вторая выборка X0: НР(Ти)
% Значения обеих выборок изменяются линейно и случайно
% со скоростью 1+k_low в НР

%% Исходные данные

clear;

% Количество невосстанавливаемых ОИ
n = 20;
% Две выборки одинакового объёма
m = n/2;

% Значение выходного параметра по ТУ (выходное напряжение, В)
X = 20;
% Допуск на выходной параметр (+-, В)
dX = 0.01*X;

% СКО
SKO = 0.1; 

% Продолжительность испытаний, ч
T_t = 8640;
% Количетсво замеров выходного параметра
l_t =10;

% Коэффициент ускорения предустановленный
K_set = 3;

%% Экпериментальная часть предварительных испытаний

x_max = X + dX;

% Генерируем начальные значения двух выборок
% в нормальном и форсированном режимах испытаний
X0 = zeros(m, l_t);
X1 = zeros(m, l_t);

% Служебные массивы для отслеживания моментов времени условного отказа
reject_flag = zeros(m, 1);
reject_time = ones(m, 1);
reject = [reject_flag reject_time];

% Коэффициент для НР
k_degradation = dX/l_t;

% Находим среднее значение
x1_av = mean(X1);
% Вычисляем уровень условного отказа
a1_max = x1_av + 0.2*(x_max - x1_av);

% Цикл заполнения массивов испытаний
for j = 1:1:l_t
    for i = 1:1:m
        X0(i, j) = normrnd(X + k_degradation*(j - 1), SKO);
        if (reject(i, 1) == 0)
            X1(i, j) = normrnd(X + k_degradation*(j - 1), SKO);
            if (X1(i, j) > a1_max)
                reject(i, 1) = 1;
                reject(i, 2) = j;
            end
        else
             X1(i, j) = normrnd(X + k_degradation*(reject(i, 2) +...
                 ((j - 1) - reject(i, 2))*K_set), SKO);
        end
    end
end

%% Вычислительная часть

% Задаём начальное значение рассчитываемого коэффициента ускорения
d_K = 2;

K0 = K_set - d_K;
K0_start = K0;
% Шаг коэффициента ускорения
dK = 0.1;

% Количество итераций, которым задаётся максимальное приращение
S = 2*(d_K)/dK;

% Служебный массивчисленных значений статистики T
T_K = zeros(1, S);

for iteration = 1:1:S
    X1_prog = zeros(m, l_t);

    % Заполняем массив прогнозируемых значений
    for i = 1:1:m
        for j = 1:1:l_t
            if j <= reject(i, 2)
                X1_prog(i, j) = X1(i, j);
            else
                X1_prog(i, j) = X1(i, reject(i, 2)) + (X1(i, j) -...
                    X1( i, reject(i, 2)) )/K0;
            end
        end
    end

    Z = cat(2, X0', X1_prog');
    Z_s = Z;
    for i = 1:1:l_t
        Z_s(i, :) = sort(Z(i, :)); 
    end


    d = zeros(1, l_t);
    V = zeros(l_t, l_t);
    
    for j = 1:1:l_t
        d_temp = 0;
        for i = 1:1:m
            if X1_prog(i, j) <= Z_s(j, m)
                d_temp = d_temp + 1;
            end
        end
        d(j) = d_temp / (2 * m) - 0.5;
        
        for i = 1:1:l_t
            if i == j
                V(i, j) = 0.25;
                continue
            end
            v_temp = 0;
            for k = 1:1:2*m
                if Z(i, k) <= Z(j, m)
                    v_temp = v_temp + 1;
                end
                if Z(j, k) <= Z(i, m)
                    v_temp = v_temp + 1;
                end
            end
            V(i, j) = v_temp / (2 * m) - 0.25;
        end
        
    end

    T = 2 * m * d / V * d';
    
    T_K(iteration) = T;
    K0 = K0 + dK;
end


% Определяем минимальное значение статистики
[T_min, i_min] = min(T_K);

% Определяем коэффициент ускорения соответственно минимуму статистики
K_true = K0_start + (i_min - 1)*dK;

% Вычисляем время испытаний
T_t_true = T_t/K_true;

%% Ускоренные испытания
X2 = zeros(n, l_t);

k_degradation = dX/l_t/K_true;

for j = 1:1:l_t
    for i = 1:1:n
         X2(i, j) = normrnd(X + k_degradation*K_set*(j - 1), SKO);
    end
end

% Количество отказавших ОИ
R = 0;

for i = 1:1:n
    if X2(i, l_t) > x_max
        R = R + 1;
    end
end

% ВБР

P = 1 - R/n;

%% Вывод
clc;

figure('Name','Методика испытаний','NumberTitle','off');
subplot(3, 1, 1);
plot(1:1:l_t, X1(:,:), 1:1:l_t, X0(:,:));
title('Значение параметра');
xlabel('Номер измерения');
ylabel('Напряжение');

subplot(3, 1, 2);
plot(1:1:l_t, X1_prog(:,:), 1:1:l_t, X0(:,:));
title('Значение параметра');
xlabel('Номер измерения');
ylabel('Напряжение');

subplot(3, 1, 3);
plot(K0_start:dK:K0_start + dK*(length(T_K) - 1), T_K);
title('Значения статистики T');
xlabel('Коэффициент ускорения');
ylabel('Статистика T');


disp(['Предустановленный коэффициент ускорения равен: ', num2str(K_set)]);
disp(['Вычисленный коэффициент ускорения равен: ',  num2str(K_true)]);
disp(['Длительность испытаний: ', num2str(T_t_true, '%4.2f'),' ч. или ',...
    num2str(T_t_true/24, '%2.2f'),' сут., ', num2str(rem(T_t_true,24),...
    '%2.2f'),' ч.']);
disp(['ВБР: ', num2str(P)]);