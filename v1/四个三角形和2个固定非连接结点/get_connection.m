%%判断是否连通
function is_connec = get_connection(sensor_mat)
global N;
global r;
global r1;

%删除上次的邻接矩阵文件
delete('adjacencyMatrix_dis.mat');
delete('adjacencyMatrix.mat');

%%计算连通性
adjacencyMatrix=zeros(N,N);%定义传感器互联邻接矩阵
adjacencyMatrix_dis=zeros(N,N);%定义传感器互联邻接矩阵
for i=1:1:N
    for j=(i+1):1:N   %因为有多种半径  连通得进行处理
        dis_squ = (sensor_mat(1,i)-sensor_mat(1,j))^2 +...
        (sensor_mat(2,i)-sensor_mat(2,j))^2;
                    %把距离加入
        if i==1 && j==2%都是7.5的
            if dis_squ <= (2*r1)^2;%；两点之间可以感知
                adjacencyMatrix(i,j)=1;%无向图
                adjacencyMatrix(j,i)=1;
                %加入距离
                adjacencyMatrix_dis(i,j)=dis_squ^0.5;%无向图
                adjacencyMatrix_dis(j,i)=dis_squ^0.5;
            end
        elseif i==1 && j~=2%有7.5的也有5的
            if  dis_squ <= (r + r1)^2;%；两点之间可以感知
                adjacencyMatrix(i,j)=1;%无向图
                adjacencyMatrix(j,i)=1;
                
                 %加入距离
                adjacencyMatrix_dis(i,j)=dis_squ^0.5;%无向图
                adjacencyMatrix_dis(j,i)=dis_squ^0.5;
            end
        elseif i==2 && j>2%有个7.5 剩下的都是5
            if dis_squ<=(r + r1)^2;%；两点之间可以感知
                adjacencyMatrix(i,j)=1;%无向图
                adjacencyMatrix(j,i)=1;
                
                %加入距离
                adjacencyMatrix_dis(i,j)=dis_squ^0.5;%无向图
                adjacencyMatrix_dis(j,i)=dis_squ^0.5;
            end
        else%的全是5的
            if dis_squ <= (2*r)^2;%；两点之间可以感知
                adjacencyMatrix(i,j)=1;%无向图
                adjacencyMatrix(j,i)=1;
                
                %加入距离
                adjacencyMatrix_dis(i,j)=dis_squ^0.5;%无向图
                adjacencyMatrix_dis(j,i)=dis_squ^0.5;
            end
        end
    end
end 

%保存最后一次的无线图邻接矩阵
save adjacencyMatrix_dis.mat adjacencyMatrix_dis;
save adjacencyMatrix.mat adjacencyMatrix;
S=zeros(N,N);
for m=1:1:N-1
    S=S+adjacencyMatrix^m;  %全部加到S
end;

%S=M+M^2+M^3+...+M^(N-1)；
%其中N是M的行数或列数
%若S中有元素为零，则不连通；
%S中无零，则连通，
%《基于邻接矩阵图的连通性判定准则》查到的
is_connec = 0;%初始化为不连通
if all(all(S))==1 %判断S向量是不是全为非零元素
    is_connec = 1;%全为非0  则连通
end
end