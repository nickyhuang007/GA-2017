%f(x)=x^a+b  ����Ҫͨ���Ŵ��㷨���a��b
price=[19.89 18.71 18.1 18.61 18.38 19.76 20.16 20.08 20.15 23.28 24.18 23.87 23.8 23.99 22.84];
%ÿ����ƽ�����̼�
sold=[95.3 138 74.9 103 65.7 184 118 77.8 216 252 516 300 312 291 234];
%��ƽ��������

a1=0; a2=1;      %��ڶ�����C���Գ�����һ��
b1=10; b2=25;
m1=14; m2=18;
K=10;          %Kλ��Ⱥ������� 
D=m1+m2;       %2���ƴ���λ��
t=1;
G=randint(K,D);%��ʼ��Ⱥ2���Ʊ���
while t<2000   %�����Ŵ�����
U1=ones(K,m1);
for i=1:1:K
  U1(i,:)=G(i,1:m1);
end  
U2=ones(K,m2);
for i=1:1:K
  U2(i,:)=G(i,(m1+1):D);
end 
X=ones(K,2);

W=1;
 
for i=1:1:K 
 S1=0;
  for j=1:1:m1
    S1=S1+U1(i,j)*power(2,(m1-j));
  end  
  X(i,1)=a1+(S1)*(a2-a1)/(power(2,m1)-1);
end    
for i=1:1:K 
 S2=0;
   for j=1:1:m2
     S2=S2+U2(i,j)*power(2,(m2-j));
   end  
  X(i,2)=b1+(S2)*(b2-b1)/(power(2,m2)-1);  %XΪ��������δ֪�������ֵ
end                                         
 evaluate=ones(K,1);
for i=1:1:K 
  S3=0;
    for j=1:1:15  
      S3=S3+W/abs((power(sold(1,j),X(i,1))+X(i,2))-price(1,j));%�˴�WΪһ��Ȩ�س��������ƻ���֮�������죨Ĭ��Ϊ1��
    end
  evaluate(i,1)=S3;
end                      %��Ӧ�ȼ���  
F=0;
for i=1:1:K
  F=F+evaluate(i,1);    
end                      %������Ӧ���ܺ�
P=ones(K,1);
for i=1:1:K
  P(i,1)=evaluate(i,1)/F;
end                      %������屻ѡ��ĸ���
SP=ones(K,1);           %sum of probility ���ø���ʱ�ĸ���֮�� 
S4=0;
for i=1:1:K
  S4=S4+evaluate(i,1);
  SP(i,1)=S4;
end  
R=rand(K,1);            %����K�������
G1=ones(K,D);          %������һ��(�����ѡ�񣬻�δ���졢����)
l=1;                     %��λ��
for i=1:1:K
    if R(i,:)<=SP(1,:)  
       G1(i,:)=G(1,:);
     else if R(i,:)<=SP(2,:)
             G1(i,:)=G(2,:);
           else if R(i,:)<=SP(3,:)
                   G1(i,:)=G(3,:);
                 else if R(i,:)<=SP(4,:)
                         G1(i,:)=G(4,:);
                       else if R(i,:)<=SP(5,:)
                               G1(i,:)=G(5,:);
                             else if R(i,:)<=SP(6,:)
                                     G1(i,:)=G(6,:);
                                   else if R(i,:)<=SP(7,:)
                                           G1(i,:)=G(7,:);
                                         else if R(i,:)<=SP(8,:)
                                                 G1(i,:)=G(8,:);
                                               else if R(i,:)<=SP(9,:)
                                                       G1(i,:)=G(9,:);
                                                     else if R(i,:)<=SP(10,:)
                                                             G1(i,:)=G(10,:);
                                                         end
                                                   end
                                             end
                                       end
                                 end
                           end
                     end
               end
         end
    end
end
R1=randint(1,1,[1 D]);   %ѡ�񽻻�λ��
P1=0.25;                   %���������ĸ���
R2=rand(K,1);            %����K�������
Count=zeros(K,1);        %ѡ����ѡ�еĸ���
c=0;                      %ͳ�Ʊ�ȡ���ĸ�����                     
for i=1:1:K
  if R2(i,1)<P1  
    Count(i,1)=1;
    c=c+1;
  end
end
Carry=ones(c,D);          %Carry���ر�ѡ�и���
c1=1;                     %ͳ���������
for i=1:1:K 
      if c1<=c
       break   
       else if Count(i,1)==1 
         Carry(c1,:)=G1(i,:);
         c1=c1+1;    
           end
      end
end
tempory1=ones(1,(D-R1));                         %��ʱ����
tempory2=ones(1,(D-R1));  
if c==0
else if c==1
    else if c==10
          for i=2:1:c
           tempory1=Carry(i,R1:D);
           tempory2=Carry((i-1),R1:D);
           Carry(i,R1:D)=tempory2;
           Carry((i-1),R1:D)=tempory1;          
           tempory1=ones(1,(D-R1));                        
           tempory2=ones(1,(D-R1));
          end
         else  for i=2:1:c
                   tempory1=Carry(i,R1:D);
                   tempory2=Carry((i-1),R1:D);
                   Carry(i,R1:D)=tempory2;
                   Carry((i-1),R1:D)=tempory1;          
                   tempory1=ones(1,(D-R1));                        
                   tempory2=ones(1,(D-R1));             %��ɽ���
             end
        end
    end
end
c2=1;                                            %������c1
for i=1:1:K
      if c2==c
       break   
       else if Count(i,1)==1 
         G1(i,:)=Carry(c2,:);
         c2=c2+1;    
           end
      end
end
R3=rand(K,D);                                  %ѡȡ�����
P2=0.05;                                         %������
tempory3=0;
for i=1:1:K
   for j=1:1:D
     if R3(i,j)<P2
       tempory3=~logical(G1(i,j));
       G1(i,j)=tempory3;
       tempory3=0;
     end
   end
end                                                %�������֤��ȷ
G=G1;
G1=ones(K,D);                                    %ע���ѽ�G1��λ
t=t+1;
end
U1=ones(K,m1);
for i=1:1:K
  U1(i,:)=G(i,1:m1);
end  
U2=ones(K,m2);
for i=1:1:K
  U2(i,:)=G(i,(m1+1):D);
end 
X=ones(K,2);
for i=1:1:K ;
 S1=0;
  for j=1:1:m1
    S1=S1+U1(i,j)*power(2,(m1-j));
  end  
  X(i,1)=a1+(S1)*(a2-a1)/(power(2,m1)-1);
end    
for i=1:1:K 
 S2=0;
   for j=1:1:m2
     S2=S2+U2(i,j)*power(2,(m2-j));
   end  
  X(i,2)=b1+(S2)*(b2-b1)/(power(2,m2)-1); 
end                                         
 evaluate=ones(K,1);
for i=1:1:K 
  S3=0;
    for j=1:1:15  
      S3=S3+W/abs((power(sold(1,j),X(i,1))+X(i,2))-price(1,j));
    end
  evaluate(i,1)=S3;
end
evaluate2=evaluate;
for i=2:1:K
   if evaluate2(i,1)>evaluate2((i-1),1)
       i=i+1;
   else evaluate2(i,1)=evaluate2((i-1),1);
       i=i+1;                                 %�ж����Ž�(���Ž�Ϊevaluate2(K,1))
   end
end 
for i=1:1:K
  if evaluate(i,1)==evaluate(K,1)
    a=X(i,1)
    b=X(i,2) 
    break
  end
end  
