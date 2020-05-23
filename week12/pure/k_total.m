close all;clear;clc;

% 计算凝聚之后的总刚度矩阵，因此此时的单元节点编号应该是从A，B，C，D，O的顺序

% 所有物理量均采用标准单位

% 输入弹模和泊松比和厚度
E = 1e11; % 单位为Pa
v = 0.3;
t = 0.01; % 单位为m

% 输入节点坐标值列表
% 平面单元仅考虑x,y两个方向
% 格式：node_num,x,y
node=[1,0.05,0.05
    2,-0.05,0.05
    3,-0.05,-0.05
    4,0.05,-0.05
    5,0,0];

% 输入单元列表
% 平面3结点三角形单元
% 格式：elem_num,node1,node2,node3
elem=[1,1,2,5;
    2,2,3,5;
    3,3,4,5;
    4,4,1,5];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 计算三角形单刚
elem_num=size(elem,1);
elem_k=zeros(6,6,elem_num); % 所有的单刚

for i=1:elem_num
    % 计算每个单元的单刚
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

% 将单刚组装进总刚
node_num=size(node,1);
total_k=zeros(2*node_num,2*node_num);

for i=1:elem_num
    % 生成单元的定位向量
    posi_v=zeros(6,1);
    for j=1:3
        node_n=elem(i,j+1);
        node_row=find(node(:,1)==node_n);    % 节点号node_n不一定等于该节点在节点数组node中的行号node_row
        % 编号为n的节点node_n其两个位移编码是按照所在行号node_row确定，而不是按照节点号node_n确定
        posi_v(2*j-1)=node_row*2-1;
        posi_v(2*j)=node_row*2;
    end
    
    % 单刚装进总刚
    for j=1:6
        for k=1:6
            total_k(posi_v(j),posi_v(k))=total_k(posi_v(j),posi_v(k))+elem_k(j,k,i);
        end
    end
end


%--------计算凝聚之后的K_bb_*矩阵--------------
K_bb = total_k(1:8,1:8);
K_bi = total_k(1:8,9:10);
K_ib = total_k(9:10,1:8);
K_ii = total_k(9:10,9:10);

K_bb_n=K_bb-K_bi*(inv(K_ii))*K_ib;

%--------计算凝聚之后的总刚度矩阵----------------
lambda = zeros(8,8); % 共有4个外部节点,所以lambda_s个数为4

% 局部坐标系与总体坐标系的坐标转换矩阵
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

