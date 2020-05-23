function K = func_ke_triangle(i,j,k,E,v,t)
   x=[i(1),j(1),k(1)];
   y=[i(2),j(2),k(2)];
%% ������Ϊ�����ɽ̲ĵĵ����γ���ʾ
%% x,y Ϊ3�����x�����y����
%% v0/E0Ϊ��������
%% tΪ�������

%% ��ʼ��  
    K=zeros(6,6);

%% ����ϵ��
    a = zeros(1,3);
    b = zeros(1,3); 
    c = zeros(1,3);
    D = [1 x(1) y(1);
        1 x(2) y(2);
        1 x(3) y(3)];
    A = det(D)/2; %%���
    
    %% p58 eq.2.2.7
    a(1) = x(2)*y(3)-x(3)*y(2);
    a(2) = x(3)*y(1)-x(1)*y(3);
    a(3) = x(1)*y(2)-x(2)*y(1);
    b(1) = y(2) - y(3);
    b(2) = y(3) - y(1);
    b(3) = y(1) - y(2); 
    c(1) = -x(2) + x(3);
    c(2) = -x(3) + x(1);
    c(3) = -x(1) + x(2);

    for i=1:3
        for j=1:3
            Kij = krs(i,j);
            K(2*i-1:2*i,2*j-1:2*j) = Kij;
        end
    end
    
% % K

%% p63 ��ʽ2.2.35-2.2.36
    function Krs = krs(r,s)
        Krs = zeros(2,2);
        Krs(1,1) = b(r)*b(s) + (1-v)/2*c(r)*c(s);   
        Krs(2,1) = v*c(r)*b(s) + (1-v)/2*b(r)*c(s);
        Krs(1,2) = v*b(r)*c(s) + (1-v)/2*c(r)*b(s);
        Krs(2,2) = c(r)*c(s) + (1-v)/2*b(r)*b(s);
        Krs = Krs*E*t/4/(1-v*v)/A;
    end
end 

