%%先把垂足求出来
function [foot_x,foot_y] =  get_foot(x1,y1,k,b)
syms foot_x foot_y x_t y_t k_t b_t;
x_t = x1;
y_t = y1;
k_t = k;
b_t = b;
[foot_x,foot_y] = solve('k_t*foot_x + b_t = foot_y','(foot_y-y_t)/(foot_x-x_t)=(-1/6)','foot_x','foot_y');
% [foot_x,foot_y] = solve(foot_x + 3 == foot_y,foot_x - 7 == foot_y,foot_x,foot_y);
end
