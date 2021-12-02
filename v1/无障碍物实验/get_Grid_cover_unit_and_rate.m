%%得到覆盖率
function cover_rate = get_Grid_cover_unit_and_rate(sensor_mat1)
global M;
global N;
global Grid_cen_x_and_y;
global r;
global L;
global W;

Grid_cover_unit = zeros(1,M);%1*M个用于保存2500个网格的联合概率
Grid_cover_bool = zeros(M,N);%每个网格中心监测的概率  只能是0或者
for i=1:M
    for j=1:N
        if ((Grid_cen_x_and_y(1,i)-sensor_mat1(1,j))^2 + (Grid_cen_x_and_y(2,i)-sensor_mat1(2,j))^2)...
                <=r^2
            Grid_cover_bool(i,j) = 1;%代表第i网格中心点可以被第j个传感器节点覆盖
        end
    end
end

%%计算联合分布概率
for i=1:M
    p = 1;%设置为1
    for j=1:N
        p = p*(1-Grid_cover_bool(i,j));%根据公式计算
    end
    Grid_cover_unit(1,i) = 1-p;%联合覆盖概率
end

%计算覆盖率
p_sum = 0;%累加所有网格的覆盖概率
for i=1:M
    p_sum = p_sum + Grid_cover_unit(1,i);
end;
cover_rate = p_sum*(L*W/M)/(L*W);%覆盖率
end