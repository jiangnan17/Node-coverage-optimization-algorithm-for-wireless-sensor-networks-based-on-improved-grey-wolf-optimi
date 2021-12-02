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

%梯形里头的监测点  共有这三种情况
%对应列的规律
case_colu1 = (16:1:35);
case_colu2 = (17:1:34);
case_colu3 = (18:1:33);
case_colu4 = (19:1:32);
case_colu5 = (20:1:31);
case_colu6 = (21:1:30);

%对应行的规律
case_row1 = (11:1:13);
case_row2 = (14:1:19);
case_row3 = (20:1:25);
case_row4 = (26:1:31);
case_row5 = (32:1:37);
case_row6 = (38:1:40);

for i=1:L
    for k=1:W  
        %笨的方法处理梯形
        if (any(case_row1 == i)) == 1 && (any(case_colu1==k))==1
            continue;
        end
        if (any(case_row2 == i)) == 1 && (any(case_colu2==k))==1
            continue;
        end
        if (any(case_row3 == i)) == 1 && (any(case_colu3==k))==1
            continue;
        end
        if (any(case_row4 == i)) == 1 && (any(case_colu4==k))==1
            continue;
        end
        if (any(case_row5 == i)) == 1 && (any(case_colu5==k))==1
            continue;
        end
        if (any(case_row6 == i)) == 1 && (any(case_colu6==k))==1
            continue;
        end
        p = 1;%设置为1
        for j=1:N
            p = p*(1-Grid_cover_bool(i,k,j));%根据公式计算
        end
        Grid_cover_unit(i,k) = 1-p;%联合覆盖概率
    end
end


%计算覆盖率
p_sum = 0;%累加所有网格的覆盖概率
for i=1:L
    for j=1:W
%         if (any(case_row1 == i)) == 1 && (any(case_colu1==j))==1
%             continue;
%         end
%         if (any(case_row2 == i)) == 1 && (any(case_colu2==j))==1
%             continue;
%         end
%         if (any(case_row3 == i)) == 1 && (any(case_colu3==j))==1
%             continue;
%         end
%         if (any(case_row4 == i)) == 1 && (any(case_colu4==j))==1
%             continue;
%         end
%         if (any(case_row5 == i)) == 1 && (any(case_colu5==j))==1
%             continue;
%         end
%         if (any(case_row6 == i)) == 1 && (any(case_colu6==j))==1
%             continue;
%         end
        p_sum = p_sum + Grid_cover_unit(i,j);
    end
end

cover_rate = p_sum*(L*W/M)/((L*W)-450);%覆盖率

end