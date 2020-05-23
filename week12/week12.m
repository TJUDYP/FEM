clc,clear;
i = [0, 0];    % ����i,j,k�������m����λ
j = [0.05, 0.05];
k = [-0.05, 0.05];

E = 1e11; % ��λΪPa
v = 0.3;
t = 0.01; % ��λΪm

K_1 = func_ke_triangle(i,j,k,E,v,t);

%--------��������任����OBC�ĵ�Ԫ�նȾ���---------
lambda_2 = zeros(6,6); % ����3���ⲿ�ڵ�,����lambda_s����Ϊ3
lambda_2_s = [0,-1;-1,0];

i = 1;
j = 1;
while i<= 5
    while j<=5
        lambda_2(i,j) = lambda_2_s(1,1);
        lambda_2(i,j+1) = lambda_2_s(1,2);
        lambda_2(i+1,j) = lambda_2_s(2,1);
        lambda_2(i+1,j+1) = lambda_2_s(2,2);
        
        i = i+2;
        j = j+2;
    end
end

K_2 = lambda_2'*K_1*lambda_2;

%--------��������任����OCD�ĵ�Ԫ�նȾ���---------
lambda_3 = zeros(6,6); % ����3���ⲿ�ڵ�,����lambda_s����Ϊ3
lambda_3_s = [1,0;0,1];

i = 1;
j = 1;
while i<= 5
    while j<=5
        lambda_3(i,j) = lambda_3_s(1,1);
        lambda_3(i,j+1) = lambda_3_s(1,2);
        lambda_3(i+1,j) = lambda_3_s(2,1);
        lambda_3(i+1,j+1) = lambda_3_s(2,2);
        
        i = i+2;
        j = j+2;
    end
end

K_3 = lambda_3'*K_1*lambda_3;

%--------��������任����ODA�ĵ�Ԫ�նȾ���---------
lambda_4 = zeros(6,6); % ����3���ⲿ�ڵ�,����lambda_s����Ϊ3
lambda_4_s = [0,1;1,0];

i = 1;
j = 1;
while i<= 5
    while j<=5
        lambda_4(i,j) = lambda_4_s(1,1);
        lambda_4(i,j+1) = lambda_4_s(1,2);
        lambda_4(i+1,j) = lambda_4_s(2,1);
        lambda_4(i+1,j+1) = lambda_4_s(2,2);
        
        i = i+2;
        j = j+2;
    end
end

K_4 = lambda_4'*K_1*lambda_4;

