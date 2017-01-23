%% ������������ ���������
% ������ ������� X1: ��(t < T�����) + ��(t > T�����)
% ������ ������� X0: ��(��)
% �������� ����� ������� ���������� ������� � ��������
% �� ��������� 1+k_low � ��

%% �������� ������

clear;

% ���������� ������������������� ��
n = 20;
% ��� ������� ����������� ������
m = n/2;

% �������� ��������� ��������� �� �� (�������� ����������, �)
X = 20;
% ������ �� �������� �������� (+-, �)
dX = 0.01*X;

% ���
SKO = 0.1; 

% ����������������� ���������, �
T_t = 8640;
% ���������� ������� ��������� ���������
l_t =10;

% ����������� ��������� �����������������
K_set = 3;

%% ���������������� ����� ��������������� ���������

x_max = X + dX;

% ���������� ��������� �������� ���� �������
% � ���������� � ������������� ������� ���������
X0 = zeros(m, l_t);
X1 = zeros(m, l_t);

% ��������� ������� ��� ������������ �������� ������� ��������� ������
reject_flag = zeros(m, 1);
reject_time = ones(m, 1);
reject = [reject_flag reject_time];

% ����������� ��� ��
k_degradation = dX/l_t;

% ������� ������� ��������
x1_av = mean(X1);
% ��������� ������� ��������� ������
a1_max = x1_av + 0.2*(x_max - x1_av);

% ���� ���������� �������� ���������
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

%% �������������� �����

% ����� ��������� �������� ��������������� ������������ ���������
d_K = 2;

K0 = K_set - d_K;
K0_start = K0;
% ��� ������������ ���������
dK = 0.1;

% ���������� ��������, ������� ������� ������������ ����������
S = 2*(d_K)/dK;

% ��������� ��������������� �������� ���������� T
T_K = zeros(1, S);

for iteration = 1:1:S
    X1_prog = zeros(m, l_t);

    % ��������� ������ �������������� ��������
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


% ���������� ����������� �������� ����������
[T_min, i_min] = min(T_K);

% ���������� ����������� ��������� �������������� �������� ����������
K_true = K0_start + (i_min - 1)*dK;

% ��������� ����� ���������
T_t_true = T_t/K_true;

%% ���������� ���������
X2 = zeros(n, l_t);

k_degradation = dX/l_t/K_true;

for j = 1:1:l_t
    for i = 1:1:n
         X2(i, j) = normrnd(X + k_degradation*K_set*(j - 1), SKO);
    end
end

% ���������� ���������� ��
R = 0;

for i = 1:1:n
    if X2(i, l_t) > x_max
        R = R + 1;
    end
end

% ���

P = 1 - R/n;

%% �����
clc;

figure('Name','�������� ���������','NumberTitle','off');
subplot(3, 1, 1);
plot(1:1:l_t, X1(:,:), 1:1:l_t, X0(:,:));
title('�������� ���������');
xlabel('����� ���������');
ylabel('����������');

subplot(3, 1, 2);
plot(1:1:l_t, X1_prog(:,:), 1:1:l_t, X0(:,:));
title('�������� ���������');
xlabel('����� ���������');
ylabel('����������');

subplot(3, 1, 3);
plot(K0_start:dK:K0_start + dK*(length(T_K) - 1), T_K);
title('�������� ���������� T');
xlabel('����������� ���������');
ylabel('���������� T');


disp(['����������������� ����������� ��������� �����: ', num2str(K_set)]);
disp(['����������� ����������� ��������� �����: ',  num2str(K_true)]);
disp(['������������ ���������: ', num2str(T_t_true, '%4.2f'),' �. ��� ',...
    num2str(T_t_true/24, '%2.2f'),' ���., ', num2str(rem(T_t_true,24),...
    '%2.2f'),' �.']);
disp(['���: ', num2str(P)]);