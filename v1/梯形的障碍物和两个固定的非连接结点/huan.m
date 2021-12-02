function re=huan(A)
    [m n]=size(A);
    while 1
        pre_A=A;
        for i=1:m
            du=0;       %第m个元素的度
            for j=1:n
                if A(i,j)~=0
                    du=du+1;
                end
            end
            if du==1            %元素的度为1时删除这个元素，其相邻元素度减一
               A(i,:)=0;
               A(:,i)=0;
            end
        end
        if pre_A==A     %图中没有度为1的元素则退出
           break; 
        end
    end
    
    if sum(A(:))==0
        re=0;
    else
        re=1;
    end
end