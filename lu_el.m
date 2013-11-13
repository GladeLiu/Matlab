function [L,U,P] =lu_el(A,ptol)
%lu_el函数产生L，U 和P
%3100100050 刘辰昂
%其中L是单位下三角阵，U是上三角阵，P是交换矩阵。
%A为待分解矩阵
%ptol的缺省值为50*eps
[m,n]=size(A);%计算矩阵维数
if nargin<2
    ptol=50*eps;
end
if m~=n
    error('必须是方阵！！！');%进行LU分解的必须是方阵
end
pp=(1:n)';
for i=1:n-1
    [~,p]=max(abs(A(i:n,i)));
    %选主元
    %p为绝对值最大的元素所在的行
    ip=p+i-1;
    if ip~=i
        A([i ip],:)=A([ip i],:);
        pp([i ip])=pp([ip i]);%行交换
    end
    pivot=A(i,i);%取对角元素
    if abs(pivot)<ptol,error('太接近0！');%对角元素不能过小
    end
    for j=i+1:n       %LU分解
            A(j,i)=A(j,i)/pivot;
            A(j,i+1:n)=A(j,i+1:n)-A(j,i)*A(i,i+1:n);
    end
end
L=eye(n)+tril(A,-1);%输出下三角阵
U=triu(A);%输出上三角阵
E=eye(n);
P=E(pp,:);%输出交换阵

end

