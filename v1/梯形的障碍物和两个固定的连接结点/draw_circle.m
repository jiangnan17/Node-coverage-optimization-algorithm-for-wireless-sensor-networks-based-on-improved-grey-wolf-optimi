%%画圆
function draw_circle(x_pos,y_pos,iter)
global N;
global L;
global W;
global r;
global r1;%大的半径
global Grid_cen_x;
global Grid_cen_y;
global ger;
angle=0:pi/100:2*pi;%角度
for k=1:N
    plot(x_pos(k),y_pos(k),'.','MarkerSize',10);%画圆心  %先画圆心 再画圆  不然少了一个  不知咋搞的
    %右边
    if rem(k,4)==1
        text(x_pos(k)+0.5,y_pos(k),'N','FontSize',15);%进行标记
        text(x_pos(k)+2,y_pos(k)-0.5,num2str(k),'FontSize',10);%进行标记
    elseif rem(k,4)==2%右边
        text(x_pos(k)-2.5,y_pos(k),'N','FontSize',15);%进行标记
        text(x_pos(k)-1,y_pos(k)-0.5,num2str(k),'FontSize',10);%进行标记
    elseif rem(k,4)==3
        %上边
        text(x_pos(k)- 0.5,y_pos(k)+ 1.5,'N','FontSize',15);%进行标记
        text(x_pos(k)+1,y_pos(k)+0.5,num2str(k),'FontSize',10);%进行标记
    else%下边
        text(x_pos(k)-0.5,y_pos(k)-1.5,'N','FontSize',15);%进行标记
        text(x_pos(k)+1,y_pos(k)-2,num2str(k),'FontSize',10);%进行标记
    end
    hold on;%加上这句  不然少一个圆心
    if k==1||k==2
        plot(r1*cos(angle)+x_pos(k),r1*sin(angle)+y_pos(k),'r','LineWidth',2);
    else
        plot(r*cos(angle)+x_pos(k),r*sin(angle)+y_pos(k));%画圆
    end
    axis([0,50,0,50]);%横纵坐标最小和最大值  
    set(gca,'xtick',(0:1:50));%设置x坐标步长为1
    set(gca,'ytick',(0:1:50));%设置y坐标步长为1   还有个问题 横坐标数字太紧凑了
    axis square;%使横纵比例为1  这样圆才像圆  不然像椭圆
    grid on;%网格化
    hold on;%把图一直保存起来
end
% 把网格中心点的图画出来
% for i=1:L/1
%     for j=1:W/1
%         plot(Grid_cen_x(i),Grid_cen_y(j),'.','MarkerFaceColor','r');
%         hold on;%把图一直保存起来
%     end
% end

%这也是一种画正方形的方法
% line([20,30],[30,30],'linewidth',4,'color','r');
% line([30,30],[30,20],'linewidth',4,'color','r');
% line([20,30],[20,20],'linewidth',4,'color','r');
% line([20,20],[20,30],'linewidth',4,'color','r');

%画矩形
% rectangle('Position',[15,10,20,10],...%左下角坐标和边长
% 'LineWidth',2,'LineStyle','-');
% x1=[15 15 35 35];%坐标分别为  （20,30)  (30,30)..
% y1=[10 20 20 10];
% fill(x1,y1,'b')%想用线条填充，但是没实现  得用到applyhatch函数 

% 画梯形
% line([15,35],[10,10],'linewidth',4,'color','r');
% line([35,30],[10,40],'linewidth',4,'color','r');
% line([30,20],[40,40],'linewidth',4,'color','r');
% line([20,15],[40,10],'linewidth',4,'color','r');
x1 = [15,35,35,30,30,20,20,15];
y1 = [10,10,10,40,40,40,40,10];
line(x1,y1,'LineWidth',1,'LineStyle','-','color','k');
fill(x1,y1,'b');
hold on;

%画三角形
% x = [25,15,35,25];%回去的一个线，所以得四个数
% y = [50,30,30,50];
% line(x,y,'LineWidth',2,'LineStyle','-');
%fill(x,y,'b')%红色填充


if iter-1 == ger
    load adjacencyMatrix_dis.mat;
    load adjacencyMatrix.mat;
%     disp('权重');
%     disp(adjacencyMatrix_dis);
%     disp('是否存在边');
%     disp(adjacencyMatrix);
%     pause(10);
    [weight_sum, span_tree] = kruskal(adjacencyMatrix, adjacencyMatrix_dis);
    [row,colu] = size(span_tree);%得到span_tree的横纵坐标
%     disp('row');
%     disp(row);
%     disp('colu');
%     disp(colu);
%     disp('span_tree');
%     disp(span_tree);
    %注意line的写法
    for i=1:row
        line([x_pos(1,span_tree(i,1)),x_pos(1,span_tree(i,2))],[y_pos(1,span_tree(i,1)),y_pos(1,span_tree(i,2))]...
           ,'LineWidth',2);
    end
    disp('最小生成树的总权值为');
    disp(weight_sum);
end
hold on;
end