function [L,U,P] =lu_el(A,ptol)
%lu_el��������L��U ��P
%3100100050 ������
%����L�ǵ�λ��������U����������P�ǽ�������
%AΪ���ֽ����
%ptol��ȱʡֵΪ50*eps
[m,n]=size(A);%�������ά��
if nargin<2
    ptol=50*eps;
end
if m~=n
    error('�����Ƿ��󣡣���');%����LU�ֽ�ı����Ƿ���
end
pp=(1:n)';
for i=1:n-1
    [~,p]=max(abs(A(i:n,i)));
    %ѡ��Ԫ
    %pΪ����ֵ����Ԫ�����ڵ���
    ip=p+i-1;
    if ip~=i
        A([i ip],:)=A([ip i],:);
        pp([i ip])=pp([ip i]);%�н���
    end
    pivot=A(i,i);%ȡ�Խ�Ԫ��
    if abs(pivot)<ptol,error('̫�ӽ�0��');%�Խ�Ԫ�ز��ܹ�С
    end
    for j=i+1:n       %LU�ֽ�
            A(j,i)=A(j,i)/pivot;
            A(j,i+1:n)=A(j,i+1:n)-A(j,i)*A(i,i+1:n);
    end
end
L=eye(n)+tril(A,-1);%�����������
U=triu(A);%�����������
E=eye(n);
P=E(pp,:);%���������

end

