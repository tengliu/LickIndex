clear;
clc;
tic;
 readDir = 'E://cas/b6';    %direction of fits








 readPath = [readDir '\*.fit'];
 readList = dir(readPath);
 [m1, n1] = size(readList);  
dim_num=3522;                                                       %total dimensions
x=0; 
 k1=1;%0-5
 k2=1;%5-10
 k3=1;%10-15
 k4=1;%15-20
 k5=1;%20-25
 k6=1;%25-30
 k7=1;%30-40
 k8=1;%40-60
 k9=1;%60-80
 k10=1;
 k11=1;
 k12=1;
 for i1 = 1:m1                                                        %templete
    filename =[readDir '\'  readList(i1, 1).name];    
 %name of file
    %splate1=fitsread(picName);                                      %read data
    %spec=splate1(1,:);                                              %first line
     % filename='E:\cas\b6\spSpec-51691-0342-608.fit'

info = fitsinfo(filename);
l = info.PrimaryData.Keywords;

j = 1;
s = '';
while(j <= length(l))
s = strvcat(s,char(l(j,1)));
j = j + 1;
end
   

               %步长
SN_G   = strmatch('SN_G',s,'exact');   
sn(i1)=l(SN_G,2);   
 %提取信噪比,sn为原始的信噪比数组
 SN=cell2mat(sn(i1));   

if(isnumeric(SN))
    Sn=SN;                                                         %提取数值型的信噪比
else
    Sn=-1;                                                              
end
sn2(1,i1)=Sn(1,1);      
%sn2为最终可用的信噪比,为一行m1列的向量，m1为fit文件个数



%这样，针对一个特定fit文件，我们得到了他的信噪比




NAXIS1 = strmatch('NAXIS1',s,'exact');    %fits文件的维度，注意，很多fits文件的维度是不同的
NAXIS1=l(NAXIS1,2);
NAXIS1=cell2mat(NAXIS1);
COEFF0   = strmatch('COEFF0',s,'exact');  %起始波长,注意1，不同文件起始波长也往往是不一样的  
COEFF0  =l(COEFF0  ,2);                  %注意2，设该值为x,则10进制的起始波长为10^x  
COEFF0  =cell2mat(COEFF0  );
COEFF1   = strmatch('COEFF1',s,'exact');  %步长，注意，也是取了log,设该值为x,则10进制的步长为10^x 
COEFF1  =l(COEFF1  ,2);                  %第i个点的波长为，10^(COEFF0 + COEFF1*i),  具体参考见：http://classic.sdss.org/dr2/products/spectra/read_spSpec.html
COEFF1  =cell2mat(COEFF1);
%--解析头文件里内容--end


splate1=fitsread(filename);              %一般对于dr8及以前的fits文件，用默认的fitsread()函数即可
spec=splate1(1,:);                       %第一行是真正的流量，如果读第x行，就是 splate1(x,:);  

%---波长数组--begin
wave=ones(1,NAXIS1);                     %生成一个长度为NAXIS1的一维数组，用于存储波长
for i=1:NAXIS1   wave(i)=i-1;  end;      %matlab数组的下标从1开始，而不是0；数组的值为0,1,2,3,4........
logwavelength = COEFF0 + wave * COEFF1;   %开始处理波长数组
for i=1:NAXIS1   wave(i)=10^logwavelength(i);  end;  % 转变为10进制的波长，第i个点的10进制波长就是wave(i)
%---波长数组--end






index_bind_L = [4222.250,4281.375,4369.125,4452.125,4514.250,4634.000,4847.875,4977.750,5160.125,5245.650, 5312.125, 5387.500, 5696.625, 5776.625, 5876.875, 4084.750,4321.000,4092.250, 4332.500,4142.125,4142.125,5069.125,5154.125,5936.625,6189.625];
index_bind_R = [4234.750,4316.375,4420.375,4474.625,4559.250,4720.250,4876.625,5054.000,5192.625,5285.650,5352.125,5415.000,5720.375,5796.625,5909.375,4123.500,4364.750,4113.500,4353.500,4177.125,4177.125,5134.125,5196.625,5994.125,6272.125];
blue_L = [4211.000,4266.375,4359.125,4445.875,4504.250,4611.500,4827.875,4946.500,5142.625,5233.150,5304.625,5376.250,5672.875,5765.375,5860.625,4042.850,4284.750,4058.500,4284.750,4080.125,4083.875,4895.125,4895.125,5816.625,6066.625];
blue_R = [4219.750,4282.625,4370.375,4454.625,4514.250,4630.250,4847.875,4977.750,5161.375,5248.150,5315.875,5387.500,5696.625,5775.375,5875.625,4081.000,4321.000,4089.750,4321.000,4117.625,4096.375,4957.625,4957.625,5849.125,6141.625];
red_L = [4241.000,4318.875,4442.875,4477.125,4560.500,4742.750,4876.625,5054.000,5191.375,5285.650,5353.375,5415.000,5722.875,5797.875,5922.125,4129.750,4368.500,4116.000,4356.000,4244.125,4244.125,5301.125,5301.125,6038.625,6372.625];
red_R = [4251.000,4335.125,4455.375,4492.125,4579.250,4756.500,4891.625,5065.250,5206.375,5318.150,5363.375,5425.000,5736.625,5811.625,5948.125,4162.250,4421.000,4138.500,4386.000,4284.125,4284.125,5366.125,5366.125,6103.625,6415.125];

%25个区间的端点
index_L0 = get_index_function(wave,index_bind_L);
index_R0 = get_index_function(wave,index_bind_R);
blue_L0 = get_index_function(wave,blue_L);
blue_R0 = get_index_function(wave,blue_R);
red_L0 = get_index_function(wave,red_L);
red_R0 = get_index_function(wave,red_R);


licks = zeros(1,25);
for i=1:19                                   %在一个带中
wave_Mid=wave(index_L0(i):index_R0(i));
wave_Blue=wave(blue_L0(i):blue_R0(i));
wave_Red=wave(red_L0(i):red_R0(i)); 

Mid_flux=spec(index_L0(i):index_R0(i));
Blue_flux=spec(blue_L0(i):blue_R0(i));
Red_flux=spec(red_L0(i):red_R0(i));

lamda1_b=min(wave_Blue);
lamda2_b=max(wave_Blue);
lamda1_r=min(wave_Red);
lamda2_r=max(wave_Red);
lamda1_mid=min(wave_Mid);
lamda2_mid=max(wave_Mid);
mean_flux_b (1,i)= (lamda1_b+lamda2_b)/2;
mean_flux_b(2,i) = trapz(wave_Blue,Blue_flux)/(lamda2_b-lamda1_b);  
mean_flux_r(1,i)=(lamda1_r+lamda2_r)/2;
mean_flux_r(2,i)=trapz(wave_Red,Red_flux)/(lamda2_r-lamda1_r);

       %flux_mid是一个向量，存储中间带积分值
% red_flux(i) = get_flux (wave_Red(i),wave,spec);
% blue_flux(i) = get_flux(wave_Blue,wave,spec);  %由图中点对应的离散的横坐标（波长）
    temp1 =mean_flux_r(1,i)-mean_flux_b(1,i);%x2-x1
    temp2 = mean_flux_r(2,i)-mean_flux_b(2,i);%y2-y1
    tempk = temp2/temp1;
    
    FClamda = mean_flux_b(2,i)+tempk*(wave_Mid-mean_flux_b(1,i));
    if i<=19  
    licks(1,i) =real(double(trapz(wave_Mid,(1-Mid_flux./FClamda))));
    else
    licks(1,i) =real(double(-2.5*log(trapz(wave_Mid,(Mid_flux./FClamda)/(lamda2_mid-lamda1_mid)))));
    end
    
end

%这样，针对一个fit文件，我们得到了它的25维lick数组licks

% 
% %如果想画散点图，将横坐标snr，纵坐标licks（1，m）（m为lick指数的维度）画出即可
% %以下为画全部fit文件第一维lick数据散点图的代码
% 
% plot(sn2(1,i1),licks(1,1),'*');
% ylim([0 10]);
% xlim([0,80])
% hold on;
%  end;
% 



   if (sn2(1,i1)>=0&&sn2(1,i1)<=5)
       for i=1:19
       licks_with_sn_0_5(k1,i)=licks(1,i);
       
       end
       k1=k1+1;
  
   else if(sn2(1,i1)>=5&&sn2(1,i1)<=10)
           for i=1:19
           licks_with_sn_5_10(k2,i)=licks(1,i);
           
           end
           k2=k2+1;
       else if(sn2(1,i1)>=10&&sn2(1,i1)<=15)
               for i=1:19
               licks_with_sn_10_15(k3,i)=licks(1,i);
               
               end
               k3=k3+1;
           else if(sn2(1,i1)>15&&sn2(1,i1)<20)
                   for i=1:19
                       licks_with_sn_15_20(k4,i)=licks(1,i);
                   end
                   k4=k4+1;
               else if(sn2(1,i1)>20&&sn2(1,i1)<25)
                       for i=1:19
                           licks_with_sn_20_25(k5,i)=licks(1,i);
                       end
                       k5=k5+1;
                    else if(sn2(1,i1)>25&&sn2(1,i1)<30)
                       for i=1:19
                           licks_with_sn_25_30(k6,i)=licks(1,i);
                       end
                       k6=k6+1;
                         else if(sn2(1,i1)>30&&sn2(1,i1)<40)
                       for i=1:19
                           licks_with_sn_30_40(k7,i)=licks(1,i);
                       end
                       k7=k7+1;
                        else if(sn2(1,i1)>40&&sn2(1,i1)<50)
                       for i=1:19
                           licks_with_sn_40_50(k8,i)=licks(1,i);
                       end
                       k8=k8+1;
                        else if(sn2(1,i1)>50&&sn2(1,i1)<60)
                       for i=1:19
                           licks_with_sn_50_60(k9,i)=licks(1,i);
                       end
                       k9=k9+1;
                           else if(sn2(1,i1)>60&&sn2(1,i1)<70)
                       for i=1:19
                           licks_with_sn_60_70(k10,i)=licks(1,i);
                       end
                       k10=k10+1;
                            else if(sn2(1,i1)>70&&sn2(1,i1)<80)
                       for i=1:19
                           licks_with_sn_70_80(k11,i)=licks(1,i);
                       end
                       k11=k11+1;
                       else if(sn2(1,i1)>80)
                       for i=1:19
                           licks_with_sn_80(k12,i)=licks(1,i);
                       end
                       k12=k12+1;
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
end
   end
       

                                                                                     mean_flux 第一维是中点波长，第二维是平均流量
                                                                                     
%  plot(licks);
mean1=mean(licks_with_sn_0_5);
var1=var(licks_with_sn_0_5);
mean2=mean(licks_with_sn_5_10);
var2=var(licks_with_sn_5_10);
mean3=mean(licks_with_sn_10_15);
var3=var(licks_with_sn_10_15);
mean4=mean(licks_with_sn_15_20);
var4=var(licks_with_sn_15_20);
mean5=mean(licks_with_sn_20_25);
var5=var(licks_with_sn_20_25);
mean6=mean(licks_with_sn_25_30);
var6=var(licks_with_sn_25_30);
mean7=mean(licks_with_sn_30_40);
var7=var(licks_with_sn_30_40);
mean8=mean(licks_with_sn_40_50);
var8=var(licks_with_sn_40_50);
mean9=mean(licks_with_sn_50_60);
var9=var(licks_with_sn_50_60);
mean10=mean(licks_with_sn_60_70);
var10=var(licks_with_sn_60_70);
mean11=mean(licks_with_sn_70_80);
var11=var(licks_with_sn_70_80);
mean12=mean(licks_with_sn_80);
var12=var(licks_with_sn_80);








figure(1);
f1=plot(1:19,var7);
set(gca,'ylim',[-10 10 ]);
hold on;
f2=plot(1:19,var8);
f3=plot(1:19,var9);
f4=plot(1:19,var10);
f5=plot(1:19,var11);
f6=plot(1:19,var12);
% f7=plot(1:19,var7);
% f8=plot(1:19,var8);
% f9=plot(1:19,var9);
legend([f1,f2,f3,f4,f5,f6],'30-40','40-50','50-60','60-70','70-80','>80');
title('var');
xlabel('lick_index');
ylabel('lick');
figure(2);
h1=plot(1:19,mean7);
set(gca,'ylim',[-10 10 ]);

hold on;
h2=plot(1:19,mean8);
h3=plot(1:19,mean9);
h4=plot(1:19,mean10);
h5=plot(1:19,mean11);
 h6=plot(1:19,mean12);
% h7=plot(1:19,mean7);
% h8=plot(1:19,mean8);
% h9=plot(1:19,mean9);
legend([h1,h2,h3,h4,h5,h6],'30-40','40-50','50-60','60-70','70-80','>80');
title('mean');
xlabel('lick_index');
ylabel('lick_amount');
toc;
