close all;clear;clc;

% ��������֮����ܸնȾ�����˴�ʱ�ĵ�Ԫ�ڵ���Ӧ���Ǵ�A��B��C��D��O��˳��

% ���������������ñ�׼��λ

% ���뵯ģ�Ͳ��ɱȺͺ��
E = 1e11; % ��λΪPa
v = 0.3;
t = 0.01; % ��λΪm

% ����ڵ�����ֵ�б�
% ƽ�浥Ԫ������x,y��������
% ��ʽ��node_num,x,y
node=[1,0.05,0.05
    2,-0.05,0.05
    3,-0.05,-0.05
    4,0.05,-0.05
    5,0,0];

% ���뵥Ԫ�б�
% ƽ��3��������ε�Ԫ
% ��ʽ��elem_num,node1,node2,node3
elem=[1,1,2,5;
    2,2,3,5;
    3,3,4,5;
    4,4,1,5];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ���������ε���
elem_num=size(elem,1);
elem_k=zeros(6,6,elem_num); % ���еĵ���

for i=1:elem_num
    % ����ÿ����Ԫ�ĵ���
    node1=elem(i,2);
    node2=elem(i,3);
    node3=elem(i,4);
    point=zeros(3,2);
    for j=1:3
        node_n=elem(i,j+1);
        node_row=find(node(:,1)==node_n);
        point(j,:)=node(node_row,2:3);
    end
    
    elem_k(:,:,i)=func_ke_triangle(point(1,:),point(2,:),point(3,:),E,v,t);
end

% ��������װ���ܸ�
node_num=size(node,1);
total_k=zeros(2*node_num,2*node_num);

for i=1:elem_num
    % ���ɵ�Ԫ�Ķ�λ����
    posi_v=zeros(6,1);
    for j=1:3
        node_n=elem(i,j+1);
        node_row=find(node(:,1)==node_n);    % �ڵ��node_n��һ�����ڸýڵ��ڽڵ�����node�е��к�node_row
        % ���Ϊn�Ľڵ�node_n������λ�Ʊ����ǰ��������к�node_rowȷ���������ǰ��սڵ��node_nȷ��
        posi_v(2*j-1)=node_row*2-1;
        posi_v(2*j)=node_row*2;
    end
    
    % ����װ���ܸ�
    for j=1:6
        for k=1:6
            total_k(posi_v(j),posi_v(k))=total_k(posi_v(j),posi_v(k))+elem_k(j,k,i);
        end
    end
end


%--------��������֮���K_bb_*����--------------
K_bb = total_k(1:8,1:8);
K_bi = total_k(1:8,9:10);
K_ib = total_k(9:10,1:8);
K_ii = total_k(9:10,9:10);

K_bb_n=K_bb-K_bi*(inv(K_ii))*K_ib;

%--------��������֮����ܸնȾ���----------------
lambda = zeros(8,8); % ����4���ⲿ�ڵ�,����lambda_s����Ϊ4

% �ֲ�����ϵ����������ϵ������ת������
lambda_s_1 = [1,0;0,1];
lambda_s_2 = [0,-1;-1,0];
lambda_s_3 = [1,0;0,1];
lambda_s_4 = [0,1;1,0];

lambda_s = [lambda_s_1, lambda_s_2, lambda_s_3, lambda_s_4];


i = 1;
j = 1;
while i<=7
    while j<=7
        lambda(i,j) = lambda_s(1,i);
        lambda(i,j+1) = lambda_s(1,i+1);
        lambda(i+1,j) = lambda_s(2,i);
        lambda(i+1,j+1) = lambda_s(2,i+1);
        i = i+2;
        j = j+2;
    end
end

K = lambda'*K_bb_n*lambda;

