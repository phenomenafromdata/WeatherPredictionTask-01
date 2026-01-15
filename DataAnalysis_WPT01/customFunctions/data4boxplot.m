function [box,low_whisk,high_whisk,data2plot]=data4boxplot(data)

% computes the stats needed to produce a boxplot from a given data set.

% by  DRL, 2021

% 7-5-2025: added "if" to handle the case of all percentiles being equal.


box = prctile(data,[25, 50, 75]);
iqr=box(3)-box(1);  %inter-quartile range (iqr)
bounds=[box(1)-1.5*iqr box(3)+1.5*iqr]; %criteria for outliers (1.5 times iqr)
data2plot=data>bounds(1) & data<bounds(2); %logical excluding outliers
data2=data(data2plot);
low_whisk=min(data2);
high_whisk=max(data2);


if iqr==0 % case when all three percentiles are the same.
    low_whisk=box(1);
    high_whisk=box(1);
end


%to plot:

% xRange=[1 2];
% 
% plot(xRange,[box(1) box(1)],'linewidth',2,'color','k')
% plot(xRange,[box(2) box(2)],'linewidth',2,'color','r')
% plot(xRange,[box(3) box(3)],'linewidth',2,'color','k')
% 
% plot([xRange(1) xRange(1)],[box(1) box(3)],'linewidth',2,'color','k')
% plot([xRange(2) xRange(2)],[box(1) box(3)],'linewidth',2,'color','k')
% 
% plot([mean(xRange) mean(xRange)],[box(1) low_whisk],'linewidth',2,'linestyle','-','color','k')
% plot([mean(xRange) mean(xRange)],[box(3) high_whisk],'linewidth',2,'linestyle','-','color','k')
