%%画圆
function draw_circle(x_pos,y_pos)
global N;
global L;
global W;
global r;
global Grid_cen_x;
global Grid_cen_y;
angle=0:pi/100:2*pi;%角度
for k=1:N
    plot(x_pos(k),y_pos(k),'.');%画圆心  %先画圆心 再话圆  不然少了一个  不知咋搞的
    plot(r*cos(angle)+x_pos(k),r*sin(angle)+y_pos(k));%画圆
    axis([0,50,0,50]);%横纵坐标最小和最大值
    set(gca,'xtick',(0:1:50));%设置x坐标步长为1
    set(gca,'ytick',(0:1:50));%设置y坐标步长为11
    axis square;%使横纵比例为1  这样圆才像圆  不然像椭圆
    grid on;%网格化
    hold on;%把图一直保存起来
end
%把网格中心点的图画出来
% for i=1:L/1
%     for j=1:W/1
%         plot(Grid_cen_x(i),Grid_cen_y(j),'.','MarkerFaceColor','r');
%         hold on;%把图一直保存起来
%     end
% end
end