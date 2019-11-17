%f(x)=x^a+b  其中要通过遗传算法求出a、b
price=[19.89 18.71 18.1 18.61 18.38 19.76 20.16 20.08 20.15 23.28 24.18 23.87 23.8 23.99 22.84];
%每股周平均收盘价
sold=[95.3 138 74.9 103 65.7 184 118 77.8 216 252 516 300 312 291 234];
%周平均交易量

a1=0; a2=1;      %与第二步中C语言程序定义一致
b1=10; b2=25;
m1=14; m2=18;
K=10;          %K位种群个体个数 
D=m1+m2;       %2进制串总位数
t=1;
G=randint(K,D);%初始种群2进制编码
while t<2000   %控制遗传代数
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
  X(i,2)=b1+(S2)*(b2-b1)/(power(2,m2)-1);  %X为解码后的两未知量的随机值
end                                         
 evaluate=ones(K,1);
for i=1:1:K 
  S3=0;
    for j=1:1:15  
      S3=S3+W/abs((power(sold(1,j),X(i,1))+X(i,2))-price(1,j));%此处W为一个权重常量，控制基因之间个体差异（默认为1）
    end
  evaluate(i,1)=S3;
end                      %适应度计算  
F=0;
for i=1:1:K
  F=F+evaluate(i,1);    
end                      %计算适应度总和
P=ones(K,1);
for i=1:1:K
  P(i,1)=evaluate(i,1)/F;
end                      %计算个体被选择的概率
SP=ones(K,1);           %sum of probility 到该个体时的概率之和 
S4=0;
for i=1:1:K
  S4=S4+evaluate(i,1);
  SP(i,1)=S4;
end  
R=rand(K,1);            %建立K个随机数
G1=ones(K,D);          %定义子一代(刚完成选择，还未变异、交换)
l=1;                     %计位器
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
R1=randint(1,1,[1 D]);   %选择交换位点
P1=0.25;                   %发生交换的概率
R2=rand(K,1);            %产生K个随机数
Count=zeros(K,1);        %选出被选中的个体
c=0;                      %统计被取出的个体数                     
for i=1:1:K
  if R2(i,1)<P1  
    Count(i,1)=1;
    c=c+1;
  end
end
Carry=ones(c,D);          %Carry承载被选中个体
c1=1;                     %统计载入个数
for i=1:1:K 
      if c1<=c
       break   
       else if Count(i,1)==1 
         Carry(c1,:)=G1(i,:);
         c1=c1+1;    
           end
      end
end
tempory1=ones(1,(D-R1));                         %临时变量
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
                   tempory2=ones(1,(D-R1));             %完成交换
             end
        end
    end
end
c2=1;                                            %类似于c1
for i=1:1:K
      if c2==c
       break   
       else if Count(i,1)==1 
         G1(i,:)=Carry(c2,:);
         c2=c2+1;    
           end
      end
end
R3=rand(K,D);                                  %选取变异点
P2=0.05;                                         %变异率
tempory3=0;
for i=1:1:K
   for j=1:1:D
     if R3(i,j)<P2
       tempory3=~logical(G1(i,j));
       G1(i,j)=tempory3;
       tempory3=0;
     end
   end
end                                                %到这儿验证正确
G=G1;
G1=ones(K,D);                                    %注意已将G1复位
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
       i=i+1;                                 %判断最优解(最优解为evaluate2(K,1))
   end
end 
for i=1:1:K
  if evaluate(i,1)==evaluate(K,1)
    a=X(i,1)
    b=X(i,2) 
    break
  end
end  
