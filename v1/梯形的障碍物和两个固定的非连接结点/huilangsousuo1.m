%%主程序
clc;
clear ;
close all;


global N;
global M;
global L;
global W;
global r;
global r1;
global Grid_cen_x;
global Grid_cen_y;
global Grid_cen_x_and_y;
global ger;%全局的总迭代数
L = 50;%长
W = 50;%宽
%假设1平方米一个网格
M = 2500;%网格总数
r = 5;%感知半径为5
r1 = 7.5;%大圆的感知半径为7.5
N = 28;%r=5的有25个 实际计算差不多23个  r=7.5的来两个，第1  2个
sizepop = 50;%种群规模
dimension = 2;% 空间维数
ger = 50;% 最大迭代次数
pos_limit = [0, 50];            % 设置位置参数限制


pop_new = zeros(dimension,N,sizepop);%新的种群
%%预分配内存
wolf_pop = zeros(dimension,N,sizepop);%种群
wolf_pop_temp =  zeros(dimension,N,sizepop);%临时的一个种群

sensor_fix1 = zeros(2,1);%固定的传感结点1  第1个  坐标为(20,20)
sensor_fix2 = zeros(2,1);%固定的传感结点2  第2个  坐标为(45,20)
sensor_fix1(1,1) = 7.5;  %配合梯形 让面积尽可能的大
sensor_fix1(2,1) = 37;
sensor_fix2(1,1) = 42.5;
sensor_fix2(2,1) = 37;

%func1 = 6*x-80;%左边的三角形的斜线
% 求这条斜线的对称点，所以可以知道该对称点的斜线的k=-1/6 ,可以求得对称点
%对称点 x1_sym = (960-35*x+12*y)/37;   y1_sym = (6*y + x - a)/6;


%func2 = -6*x+220;%右边的三角形的斜线
% 求这条斜线的对称点，所以可以知道该对称点的斜线的k=1/6 ,可以求得对称点
%对称点 x2_sym = (12*y+35*x-2640)/-37;   y2_sym = (6*y - x + a)/6;


%求出新梯形的四个点
syms x;%先定义一个变量
%左上角
k1 = 6;
b1 = -80;

[k1,b1] = get_new_function(k1,b1);
x1 = solve(k1*x+b1==45,x);
%左下角
x2 = solve(k1*x+b1==5,x);

%右上角
k2 = -6;
b2 = 220;
[k2,b2] = get_new_function(k2,b2);
x3 = solve(k2*x+b2==45,x);
%右下角
x4 = solve(k2*x+b2==5,x);
% disp(x1);
% disp(x2);
% disp(x3);
% disp(x4);
%上面测试完成  完全正确
point = zeros(4,2);%用来存四个点
point(1,:) = [x1,45];
point(2,:) = [x2,5];
point(3,:) = [x3,45];
point(4,:) = [x4,5];


%梯形的四边比例  上 1  下 2  左右都为3


%生成初始的种群  已经转换成三维了  种群是三维     2*N*sizepop
for i=1:sizepop
    x = L*rand(1,N);%随机生成节点的横坐标
    y = W*rand(1,N);%随机生成节点的列坐标
    %位置的初始化
    wolf_pop(1,:,i) = x;
    wolf_pop(2,:,i) = y;
    
    %进行固定处理
    
    wolf_pop(:,1,i) = sensor_fix1;
    wolf_pop(:,2,i) = sensor_fix2;
end


%种群每个个体的一些处理   避免布置在障碍物里面  同时处理靠近障碍物的周围的物体  做到极致相切
for i=1:sizepop
    count1 = 0;%用于统计这个个体中到底有多少个传感器节点是进入所在区域的
    count2 = 0;%用于统计这个个体中到底有多少个传感器节点是进入所在区域的
    count3 = 0;%用于统计这个个体中到底有多少个传感器节点是进入所在区域的
    for j=1:N
        if (wolf_pop(1,j,i)>=point(1,1)&&wolf_pop(1,j,i)<=point(3,1)) && (wolf_pop(2,j,i)>=5&&wolf_pop(2,j,i)<=45)%说明在障碍物里面了
            %也是要均匀的放置 均匀放在障碍物的周围 
            %同时这里还需要更加严谨的进行 放在上边还是下边做个布置
            %比例是1:1:1:1
            count1 = count1 + 1;
            if rem(count1,4)==1%做到相切
                wolf_pop(2,j,i) = 45;  %处理上半部分
            elseif rem(count1,4)==2
                wolf_pop(2,j,i) = 5;%处理下部分
            elseif rem(count1,4) == 3 %解决左斜边  为了方便起见，直接把该点的y假设在斜线上求出x
                wolf_pop(1,j,i) = (wolf_pop(2,j,i) - b1)/k1;
            elseif rem(count1,4) == 0%解决右斜边
                wolf_pop(1,j,i) = (wolf_pop(2,j,i) - b2)/k2;
            end
        end
        
        %处理两个斜边
        %先处理左斜边  比例是2:1
         if (wolf_pop(1,j,i)>=point(2,1)&&wolf_pop(1,j,i)<=point(1,1)) && (wolf_pop(2,j,i)>=5&&wolf_pop(2,j,i)<=45)
             if (k1 * wolf_pop(1,j,i) + b1) >= wolf_pop(2,j,i)
                 count2 = count2 + 1;
                 if rem(count2,3) == 1 || rem(count2,3) == 0
                     wolf_pop(2,j,i) = (k1 * wolf_pop(1,j,i) + b1);%让其在平行的那条斜线上 
                 elseif rem(count2,3) == 2%放入下边
                     wolf_pop(2,j,i) = 5;
                 end
             end
         end
         
         
         %处理右斜边  比例是2:1
         if (wolf_pop(1,j,i)>=point(3,1)&&wolf_pop(1,j,i)<=point(4,1)) && (wolf_pop(2,j,i)>=5&&wolf_pop(2,j,i)<=45)
             if (k2 * wolf_pop(1,j,i) + b2) >= wolf_pop(2,j,i)
                 count3 = count3 + 1;
                 if rem(count3,3) == 1 || rem(count2,3) == 0 %左三角形 斜线上
                     wolf_pop(2,j,i) = (k2 * wolf_pop(1,j,i) + b2);
                 elseif rem(count3,3) == 2
                     wolf_pop(2,j,i) = 5;
                     %让其在平行的那条斜线上
                 end
             end
         end
    end
end
                
                
            
wolf_pop_public = wolf_pop;%公共的初始种群

%保存种群到wolf_pop_public.mat文件中
save wolf_pop_public.mat wolf_pop_public;
%%初始的部署后画图  拿第一个粒子拿去初始画图

%求网格中心坐标
X_mat = (0:1:50);%x矩阵
Y_mat = (0:1:50);%y矩阵
Grid_cen_x = zeros(1,L/10);%网格中心点x坐标
Grid_cen_y = zeros(1,W/10);%网格中心点y坐标
%前后两者相加之和除以2
for i=1:L/1
    Grid_cen_x(i) = (X_mat(i)+X_mat(i+1))/2;
    Grid_cen_y(i) = (Y_mat(i)+Y_mat(i+1))/2;
end

%%把横纵坐标丢到一个二维矩阵当中
%用于转坐标  第一行放x轴 第二行放y轴，同一列放一个点坐标
%且先存靠近x坐标的第一行，然后往上存第二行
%网格中心坐标
%从下往上数第一行为第一行
Grid_cen_x_and_y = zeros(L,W,2);%共2500个网格中心，但是每个网格中心有x,y
for i=1:L/1
    for j=1:W/1
        Grid_cen_x_and_y(i,j,1) = Grid_cen_x(j);%1代表x
        Grid_cen_x_and_y(i,j,2) = Grid_cen_y(i);%把y坐标放到第二行
    end
end


x_pos = wolf_pop(1,:,1);%第一个粒子   粒子即是解  的x坐标
y_pos = wolf_pop(2,:,1);%第一个粒子   粒子即是解  的y坐标

iter = 1;%第一代  从1开始
%进行画图
figure(1);
draw_circle(x_pos,y_pos,iter);
title('初始化部署图');

%把随机部署的点的坐标放到矩阵当中
%sensor_mat = zeros(2,N);%预分配内存
sensor_mat = wolf_pop(:,:,1);%第一代中的第一只狼



%%得到联合概率
cover_rate =  get_Grid_cover_unit_and_rate(sensor_mat);
disp('覆盖率：');
disp(cover_rate);



%%计算连通性  不算了
% is_connec = get_connection(sensor_mat);
% if is_connec==1
%     disp('连通');
% else
%     disp('非连通');
% end





%%初始化种群历史值为  无穷小   inf为无穷大
best_all_ger_fitness = -inf;                      % 种群历史最佳适应度  
best_indivi = zeros(2,N);                 % 保存优秀个体

%%在上面已经画图了
%%以上为画出初始化时，50个粒子开始的位置
wolf_one_dis =  zeros(2,N);%第一只狼距离食物的距离
wolf_two_dis = zeros(2,N);%第二只狼距离食物
wolf_three_dis = zeros(2,N);%第三只狼距离食物
wolf_one = zeros(2,N);%第一只狼
wlof_two = zeros(2,N);%第二只狼
wolf_three = zeros(2,N);%第三只狼
wolf_fitness = zeros(sizepop,1);%个体当前适应度
wolf_fitness_temp = zeros(sizepop,1);%临时保存当前适应度
%% 群体更新

record_ger = zeros(ger, 1);          % 记录每次迭代中最好的适应值 
record_pop_ave = zeros(ger,1);       % 记录种群适应值的平均值



while iter <= ger

    %计算该种群中每个个体的适应值
    %factor = 2-2*(iter/ger);%随机因子
    %factor = 1-(iter/ger);%随机因子
    factor = 2-2*(1-(((ger-iter)/ger).^2)).^0.5;
    %factor = sin(((iter*pi)/ger)+pi/2)+1;
    r2 = rand(1);%随机数
    A = factor*r2-factor;
    %C = 2 * r1;
    C1 = 2*rand(dimension,N);%随机数
    C2 = 2*rand(dimension,N);
    C3 = 2*rand(dimension,N);
%     if abs(A)>1
%         disp('扩大狼群搜索范围，全局');
%     else
%         disp('狼群将包围圈缩小，局部搜索');
%     end
    for k=1:sizepop
        sensor_mat(1,:) = wolf_pop(1,:,k);
        sensor_mat(2,:) = wolf_pop(2,:,k);
        wolf_fitness(k,1) = get_Grid_cover_unit_and_rate(sensor_mat) ; % 个体当前适应度
    end
    
    [fitness_order,order_index] = sort(wolf_fitness);%对适应值进行排序
    %让三只狼保存相应的数据
    %得到第一只狼
    wolf_one = wolf_pop(:,:,order_index(sizepop,1));

    %得到第二只狼
    
    wolf_two = wolf_pop(:,:,order_index(sizepop-1,1));
  
    
    %得到第三只狼
    wolf_three = wolf_pop(:,:,order_index(sizepop-2,1));
    
   %更新排序后的种群
    for i=1:sizepop
        wolf_pop_temp(:,:,i) = wolf_pop(:,:,order_index(i,1));
    end
    wolf_pop = wolf_pop_temp;%已经更新了种群
    pop_new = wolf_pop;


%先把5只最好的给到下一代
    for j=1:5
        if rem(j,3)==1
            pop_new(:,:,j) = wolf_one;
        elseif rem(j,3)==2
            pop_new(:,:,j) = wolf_two;
        else
            pop_new(:,:,j) = wolf_three;
        end
    end
    
    [fitness_max_one,max_index_one] = max(wolf_fitness);
    record_ger(iter,1) = fitness_max_one;%保存该代中最优的适应度值
    disp(iter);
%     disp(fitness_max_one);
    disp(wolf_fitness);
    if best_all_ger_fitness < fitness_max_one
        best_all_ger_fitness = fitness_max_one;
        best_indivi = wolf_pop(:,:,max_index_one);
    end
    
    sum_fitness = sum(wolf_fitness);
    record_pop_ave(iter,1) = sum_fitness/sizepop;%算出该种群的平均适应度
    

 

    
   
    
    %只对前sizepop-3个个体进行更新 三只最好的狼依然保存到下一个种群
    for i=1:sizepop-3
        wolf_one_dis = (C1.*wolf_one - wolf_pop(:,:,i));
        wolf_two_dis = (C2.*wolf_two - wolf_pop(:,:,i));
        wolf_three_dis = (C3.*wolf_three - wolf_pop(:,:,i));
        
        X1 = wolf_one -  wolf_one_dis*A;
        X2 = wolf_two - wolf_two_dis*A;
        X3 = wolf_three - wolf_three_dis*A;
        for j=1:N
            %X1的越位更正
            if X1(1,j)>pos_limit(1,2)
                X1(1,j) = wolf_one(1,j);
            end
            if X1(2,j)>pos_limit(1,2)
                X1(2,j) = wolf_one(2,j);
            end
            if X1(1,j)<pos_limit(1,1)
                X1(1,j) = wolf_one(1,j);
            end
            if X1(2,j)<pos_limit(1,1)
                X1(2,j) = wolf_one(2,j);
            end
            %x2位置更正
            if X2(1,j)>pos_limit(1,2)
                X2(1,j) = wolf_two(1,j);
            end
            if X2(2,j)>pos_limit(1,2)
                X2(2,j) = wolf_two(2,j);
            end
            if X2(1,j)<pos_limit(1,1)
                X2(1,j) = wolf_two(1,j);
            end
            if X2(2,j)<pos_limit(1,1)
                X2(2,j) = wolf_two(2,j);
            end
            %x3越位更正
            if X3(1,j)>pos_limit(1,2)
                X3(1,j) = wolf_three(1,j);
            end
            if X3(2,j)>pos_limit(1,2)
                X3(2,j) = wolf_three(2,j);
            end
            if X3(1,j)<pos_limit(1,1)
                X3(1,j) = wolf_three(1,j);
            end
            if X3(2,j)<pos_limit(1,1)
                X3(2,j) = wolf_three(2,j);
            end
        end
        %改进  通过位置占比
%         w1 = X1./(X1+X2+X3);
%         w2 = X2./(X1+X2+X3);
%         w3 = X3./(X1+X2+X3);
%通过适应值占比
%         x1 = get_Grid_cover_unit_and_rate(X1);
%         x2 = get_Grid_cover_unit_and_rate(X2);
%         x3 = get_Grid_cover_unit_and_rate(X3);
%         w1 = x1/(x1+x2+x3);
%         w2 = x2/(x1+x2+x3);
%         w3 = x3/(x1+x2+x3);
        
        p = randperm(N,1);
%         temp = (X1.*w1 + X2.*w2 + X3.*w3);
        temp = (X1+X2+X3)/3;
        %因为这里只有一个  所以用i来进行判断
        for j=1:1
            %障碍物处理  同时把周围的障碍物也处理了
            if (temp(1,p(1,j))>=point(1,1)&&temp(1,p(1,j))<=point(3,1)) && (temp(2,p(1,j))>= 5&&temp(2,p(1,j))<=45)%说明在障碍物里面了
            %也是要均匀的放置 均匀放在障碍物的周围 
                if rem(i,4) == 1
                    temp(2,p(1,j)) = 45;%上部分
                elseif rem(i,4) == 2
                    temp(2,p(1,j)) = 5;%下部分
                elseif rem(i,4) == 3%解决左斜边
                    temp(1,p(1,j)) = (temp(2,p(1,j))-b1)/k1;
                elseif rem(i,4) == 0%右斜边
                    temp(1,p(1,j)) = (temp(2,p(1,j))-b2)/k2;
                end
            end
            
            
            %处理两个斜边
            %先处理左斜边
             if (temp(1,p(1,j))>=point(2,1)&&temp(1,p(1,j))<=point(1,1)) && (temp(2,p(1,j))>=5&&temp(2,p(1,j))<=45)
                 if (k1 * temp(1,p(1,j)) + b1) >= temp(2,p(1,j))
                     if rem(i,2) == 0%左三角 下面的边
                         temp(2,p(1,j)) = 5;
                     else
                        temp(2,p(1,j)) = (k1 * temp(1,p(1,j)) + b1);%让其在平行的那条斜线上
                     end
                 end
             end

             
             %处理右斜边
             if (temp(1,p(1,j))>=point(3,1)&&temp(1,p(1,j))<=point(4,1)) && (temp(2,p(1,j))>=5&&temp(2,p(1,j))<=45)
                 if (k2 * temp(1,p(1,j)) + b2) >= temp(2,p(1,j))
                     if rem(i,2) == 0%右三角形 下面的边
                        temp(2,p(1,j)) = 5;
                     else
                        temp(2,p(1,j)) = (k2 * temp(1,p(1,j)) + b2);%让其在平行的那条斜线上
                     end
                 end
             end
            
            
            %固定处理
            temp(:,1) = sensor_fix1;
            temp(:,2) = sensor_fix2;
            pop_new(:,p(1,j),i) = temp(:,p(1,j)); 
        end
        
        %得到新的种群
%         %对该个体进行越位处理
        for k=1:N
            if pop_new(1,k,i)> pos_limit(1,2)
                pop_new(1,k,i) = pos_limit(1,2);
            end
            if pop_new(2,k,i)>pos_limit(1,2)
                pop_new(2,k,i) = pos_limit(1,2);
            end
            if pop_new(1,k,i)< pos_limit(1,1)
                pop_new(1,k,i) = pos_limit(1,1);
            end
            if pop_new(1,k,i) < pos_limit(1,1)
                pop_new(1,k,i) = pos_limit(1,1);
            end
        end
%         disp(pop_new(:,:,i));
%         pause(5);
    end
    
    %三只最好的狼放到下一代中  把下一代中最小的三个替换掉
%     for i=1:sizepop
%         wolf_fitness(i,1) = get_Grid_cover_unit_and_rate(pop_new(:,:,i));
%     end
%     [wolf_fitness,order] = sort(wolf_fitness);
    pop_new(:,:,sizepop-2) = wolf_three;
    pop_new(:,:,sizepop-1) = wolf_two;
    pop_new(:,:,sizepop) = wolf_one;
    
    wolf_pop = pop_new;%更新种群
    
    iter = iter+1;
end

% %写入到xlsx表中
% disp('数据写入中...');
% xlswrite('huilang1.xlsx',record_ger,1);


%写入到.mat文件中



figure(2);
plot(record_ger);
title('每代中最优的适应值')

figure(3);
plot(record_pop_ave);
title('每个种群的平均适应值')


%保存成mat文件
save best_wolf1.mat best_indivi;

if get_connection(best_indivi)==1
    disp('连通');
else
    disp('非连通');
end

figure(4);
draw_circle(best_indivi(1,:),best_indivi(2,:),iter);
title('优化后的部署图');
disp(['覆盖率：',num2str(max(record_ger))]);

