close all;clear;clc;

% 所有物理量均采用标准单位

% 输入弹模和泊松比和厚度
E = 1e11; % 单位为Pa
v = 0.3;
t = 0.01; % 单位为m

% 输入节点坐标值列表
% 平面单元仅考虑x,y两个方向
% 格式：node_num,x,y
node=[1,0,0
    2,0.05,0.05
    3,-0.05,0.05
    4,-0.05,-0.05
    5,0.05,-0.05];

% 输入单元列表
% 平面3结点三角形单元
% 格式：elem_num,node1,node2,node3
elem=[1,1,2,3;
    2,1,3,4;
    3,1,4,5;
    4,1,5,2];

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


