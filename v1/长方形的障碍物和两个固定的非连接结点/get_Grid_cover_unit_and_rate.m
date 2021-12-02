%%得到覆盖率
function cover_rate = get_Grid_cover_unit_and_rate(sensor_mat1)
global M;
global N;
global Grid_cen_x_and_y;
global r;
global L;
global W;
global r1;

Grid_cover_unit = zeros(L,W);%1*M个用于保存2500个网格的联合概率
Grid_cover_bool = zeros(L,W,N);%每个网格中心监测的概率  只能是0或者
%这里暂时不考虑障碍物的影响，反正这只是过程 关键是联合分布概率
for i=1:L
    for k=1:W
        for j=1:N
            if j==1||j==2%对固定的两个传感器节点进行处理
                if ((Grid_cen_x_and_y(i,k,1)-sensor_mat1(1,j))^2 + (Grid_cen_x_and_y(i,k,2)-sensor_mat1(2,j))^2)...
                        <=r1^2
                    Grid_cover_bool(i,k,j) = 1;%代表(i,k)网格中心点可以被第j个传感器节点覆盖
                end
            else%会移动的传感器节点
                if ((Grid_cen_x_and_y(i,k,1)-sensor_mat1(1,j))^2 + (Grid_cen_x_and_y(i,k,2)-sensor_mat1(2,j))^2)...
                        <=r^2
                    Grid_cover_bool(i,k,j) = 1;%代表(i,k)网格中心点可以被第j个传感器节点覆盖
                end
            end
        end
    end
end

%%计算联合分布概率

%三角形里头的监测点
% case1 = (16:1:35);
% case2 = (17:1:34);
% case3 = (18:1:33);
% case4 = (19:1:32);
% case5 = (20:1:31);
% case6 = (21:1:30);
% case7 = (22:1:29);
% case8 = (23:1:28);
% case9 = (24:1:27);
% case10 = (25:1:26);
for i=1:L
    for k=1:W%注意这里的行列的意思  这是靠近x轴的才是第一行
        if (i>10&&i<=20)&&(k>15&&k<=35)%处理掉矩形的
            continue;
        end
        
        %笨的方法处理三角形
%         if i==31 && (any(case1==k))==1
%             continue;
%         end
%         if (i==32||i==33)&& (any(case2==k))==1
%             continue;
%         end
%         if (i==34||i==35)&& (any(case3==k))==1
%             continue;
%         end
%         if (i==36||i==37)&& (any(case4==k))==1
%             continue;
%         end
%         if (i==38||i==39)&& (any(case5==k))==1
%             continue;
%         end
%         if (i==40||i==41)&& (any(case6==k))==1
%             continue;
%         end
%         if (i==42||i==43)&& (any(case7==k))==1
%             continue;
%         end
%         if (i==44||i==45)&&k== (any(case8==k))==1
%             continue;
%         end
%         if (i==46||i==47)&&k== (any(case9==k))==1
%             continue;
%         end
%         if (i==48||i==49)&&k==(any(case10==k))==1
%             continue;
%         end
          %下面这一句可以不要
%         if (i==50)
%             continue;
%         end
        p = 1;%设置为1
        for j=1:N
            p = p*(1-Grid_cover_bool(i,k,j));%根据公式计算
        end
        Grid_cover_unit(i,k) = 1-p;%联合覆盖概率
    end
end


%计算覆盖率
p_sum = 0;%累加所有网格的覆盖概率
for i=1:L%行
    for j=1:W%列
%         if (i>10&&i<=20)&&(j>15&&j<=35)%处理掉矩形的
%             continue;
%         end
        p_sum = p_sum + Grid_cover_unit(i,j);
    end
end

cover_rate = p_sum*(L*W/M)/((L*W)-200);%覆盖率

end