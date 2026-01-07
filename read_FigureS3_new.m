clear
clc
addpath('F:\data_code\input\', 'F:\data_code\output\', 'F:\data_code\function\');

load('F:\data_code\input\NAO_month_1850_2025.mat')
year_NAO=1850:2025;
idx=find(year_NAO>=1951 & year_NAO<=2024);
NAO=mean(NAO(:,idx),1);

load('F:\data_code\input\nino34_month_1950_2025.mat')
year_nino34=1950:2025;
nino34(nino34<-100)=nan;
idx=find(year_nino34>=1950 & year_nino34<=2024);
nino341=cat(1,nino34(12,1:end-1),nino34(1:2,2:end));
nino34=mean(nino341,1);

load('F:\data_code\input\PDO_month_1854_2025.mat')
year_PDO=1854:2025;
idx=find(year_PDO>=1951 & year_PDO<=2024);
PDO=mean(PDO(:,idx),1);

load('F:\data_code\input\AO_month_1950_2025.mat')
year_AO=1950:2025;
idx=find(year_AO>=1951 & year_AO<=2024);
AO=mean(AO(:,idx),1);

load('F:\data_code\input\ONI_month_1950_2025.mat')
year_ONI=1950:2025;
idx=find(year_ONI>=1951 & year_ONI<=2024);
% ONI=mean(ONI(:,idx),1);
ONI1=cat(1,ONI(12,1:end-1),ONI(1:2,2:end));
ONI=mean(ONI1,1);


load('F:\data_code\input\IOD_month_1854_2024.mat')
year_IOD=1854:2024;
idx=find(year_IOD>=1951 & year_IOD<=2024);
IOD=mean(IOD(:,idx),1);

load('F:\data_code\input\AMO_month_1854_2025.mat')
year_AMO=1854:2025;
idx=find(year_AMO>=1951 & year_AMO<=2024);
for i=1:12
    AMO2(i,:)=detrend(AMO(i,idx));
end
AMO=mean(AMO2(:,:),1);


load('F:\data_code\input\ERA5_region_ecdf_1950_2025.mat')
load('F:\data_code\input\ERA5_global_ecdf_1950_2025.mat')
year_PI=1950:2025;
idx=find(year_PI>=1951 & year_PI<=2024);
PI=series_ERA5(idx);
PI_region=series_ERA5_region(idx,:);

load('F:\data_code\output\PI_global_GSATconstraint_1950_2025.mat')
PI_forced=post_all(2:end-1,1);

year_PI=1951:2024;
idx=find(year_PI>=1981 & year_PI<=2010);
PI=PI-mean(PI(idx));
PI_region=PI_region-mean(PI_region(idx,:),1);
PI_forced=PI_forced-mean(PI_forced(idx));

index(:,1)=(nino34(1:end-1));
index(:,2)=(PDO(1:end))-mean(PDO(idx));
index(:,3)=(IOD(1:end))-mean(IOD(idx));
index(:,4)=(NAO(1:end))-mean(NAO(idx));
index(:,5)=(AO(1:end))-mean(AO(idx));
index(:,6)=(AMO(1:end))-mean(AMO(idx));


for i=1:6
    %     [corr_index(i,1) pvalue(i,1)]=corr(index(:,i),PI(:));
    [corr_index(i,1) pvalue(i,1)]=corr((index(:,i)),PI(:)-PI_forced);
end
PI_1950=PI(:)-PI_forced;
colors = [
    0.00 0.45 0.74;
    0.85 0.33 0.10;
    0.93 0.69 0.13;
    0.49 0.18 0.56;
    0.47 0.67 0.19
    ];



%%
load('F:\data_code\input\ERA5_region_ecdf_1973_2025.mat')
load('F:\data_code\input\ERA5_global_ecdf_1973_2025.mat')
PI=series_ERA5(1:end-1);
PI_region=series_ERA5_region(1:end-1,:);
load('F:\data_code\output\PI_global_GSATconstraint_1973_2025.mat')
PI_forced=post_all(1:end-1,1);

year_PI=1973:2024;
idx=find(year_PI>=1981 & year_PI<=2010);
PI=PI-mean(PI(idx));
PI_region=PI_region-mean(PI_region(idx,:),1);
PI_forced=PI_forced-mean(PI_forced(idx));

year_indice=1951:2024;
idd1=find(year_indice>=1981 & year_indice<=2010);
idd=find(year_indice>=1973 & year_indice<=2024);
index2(:,1)=(nino34(idd));
index2(:,2)=(PDO(idd))-mean(PDO(idd1));
index2(:,3)=(IOD(idd))-mean(IOD(idd1));
index2(:,4)=(NAO(idd))-mean(NAO(idd1));
index2(:,5)=(AO(idd))-mean(AO(idd1));
index2(:,6)=(AMO(idd))-mean(AMO(idd1));


for i=1:6
    %     [corr_index(i,3) pvalue(i,3)]=corr(index2(:,i),PI(:));
    [corr_index(i,2) pvalue(i,2)]=corr((index2(:,i)),PI(:)-PI_forced);
end


index_name={'Nino34' 'PDO' 'IOD' 'NAO' 'AO' 'AMO'};
hold on
subplot(4,2,1)
b=bar(corr_index,'DisplayName','corr_index')
b(1).FaceColor=colors(1,:);
b(1).FaceAlpha=0.7;
b(1).EdgeColor='none';

b(2).FaceColor=colors(2,:);
b(2).FaceAlpha=0.7;
b(2).EdgeColor='none';

set(gca,'Xtick',1:6,'Xticklabel',index_name,'FontName','Times New Roman')
ylabel('Correlation coefficient');

legend(b,{'1951-2024' '1973-2024'})
box on
title('Correlation coefficient')
text(0.1,0.42,'a','FontWeight','bold')
set(gca,'FontName','Arial')



%%
subplot(4,2,2)
b=bar(pvalue,'DisplayName','corr_index')
b(1).FaceColor=colors(1,:);
b(1).FaceAlpha=0.7;
b(1).EdgeColor='none';

b(2).FaceColor=colors(2,:);
b(2).FaceAlpha=0.7;
b(2).EdgeColor='none';

set(gca,'Xtick',1:6,'Xticklabel',index_name,'FontName','Times New Roman')
ylabel('P-value');

legend(b,{'1951-2024' '1973-2024'})
box on
title('(P-value')
text(0.1,1.12,'b','FontWeight','bold')
set(gca,'FontName','Arial')




%%
index_name={'c' 'd' 'e' 'f' 'g' 'h'};
index_name2={'Nino34' 'PDO' 'IOD' 'NAO' 'AO' 'AMO'};

for i=1:6
    subplot(4,2,i+2)

    yyaxis left

    b=bar(1951:2024,index(:,i))
    b.FaceColor=colors(1,:);
    b.FaceAlpha=0.7;
    b.EdgeColor='none';

    hold on
    xlim([1949 2025])
    title(index_name2(i))

    maxvalue=max(max(index(:,i)),abs(min(index(:,i))));
    ylabel(index_name2(i))
    ylim([-1.1*maxvalue 1.1*maxvalue])
    text(1946,1.12*maxvalue,index_name(i),'FontWeight','bold')
    yyaxis right
    h=plot(1951:2024,PI_1950,'color','r','LineWidth',1)
    ylabel('PI residual')

    set(gca,'FontName','Arial')


end


