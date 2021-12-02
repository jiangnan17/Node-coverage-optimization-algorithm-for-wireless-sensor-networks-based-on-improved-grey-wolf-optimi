%%判断是否连通
function is_connec = get_connection(sensor_mat)
global N;
global r;
%%计算连通性
adjacencyMatrix=zeros(N,N);%定义传感器互联邻接矩阵
for i=1:1:N
    for j=(i+1):1:N
        if (sensor_mat(1,i)-sensor_mat(1,j))^2 +...
                (sensor_mat(2,i)-sensor_mat(2,j))^2<=(2*r)^2;%；两点之间可以感知
            adjacencyMatrix(i,j)=1;%无向图
            adjacencyMatrix(j,i)=1;
        end
    end
end 
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