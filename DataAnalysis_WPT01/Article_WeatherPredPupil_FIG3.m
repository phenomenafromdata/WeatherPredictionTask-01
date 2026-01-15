%% calculate strategy score using Excel data
clearvars;
analyze_behav_WeatherPred2018

analyze_strategies_WeatherPred2018
%% PLOT v1 (with Strategy NAMES)


left=0.09;
bott=0.44;
w=0.3;
h=0.45;

posits=[left bott w h;
    left+1.32*w bott 0.68*w h;
    left+2.34*w bott 0.68*w h];

xlimi=[-0.4 6.15];

XTL={'Spat. + NonSpat.' 'Single shape' ...
    'Single posit.' 'NonSpat. + random' ...
    'Spat. + random' 'Totally random'};

fname='Verdana';
fsize=14;
fsizeBig=fsize+7;
ms=20;
fig=figure;
ax=axes('Parent',fig,'TickDir','out', 'position', posits(1,:),...
    'XTick',1:6, 'xticklabel',XTL,...
    'Ydir', 'reverse', ...
    'fontsize',fsize,'fontname',fname);
hold(ax,'on')
ax.XTickLabelRotation=45;

for j=1:size(Scores_ByStrategy,2)
    jitt=randn(1,1)./20;
    plot((1:6)+jitt,Scores_ByStrategy(:,j),...
        'Marker','.','MarkerSize',ms)

end

ylabel('Fit score')

plot(xlimi,[0.1 0.1],':k')

xPos_arrow=0.095;

t1=text(xlimi(1)+0.53,0.08,{'Closer to' 'strategy'},'fontsize',fsize-3);
t1.Rotation=90;
annotation(fig,'arrow',[xPos_arrow xPos_arrow],...
    [0.8 0.9],'headstyle','plain','headwidth',6,'headlength',6);

t2=text(xlimi(1)+0.53,0.28,{'Away from' 'strategy'},'fontsize',fsize-3);
t2.Rotation=90;
annotation(fig,'arrow',[xPos_arrow xPos_arrow],...
    [0.76 0.65],'headstyle','plain','headwidth',6,'headlength',6);


%xlabel('Strategy')

xlim(xlimi)


[~,st]=min(Scores_ByStrategy);

ax=axes('parent',fig,'TickDir','out','position',posits(2,:),...
    'XTick',1:6, 'xticklabel',XTL,...
    'FontSize',fsize,'FontName',fname);
hold(ax,'on')

H=histogram(st);
H.FaceColor=[0.5 0.5 0];

xlabel('Strategy')
ylabel('Number of participants')


ms=12;

ax=axes('parent',fig,'TickDir','out','position',posits(3,:),...
    'XTick',1:6,'xticklabel',XTL,...
    'FontSize',fsize,'FontName',fname);
hold(ax,'on')

Strat_PerfLastBlock=nan(1,2);

for j=1:size(BlocksMean_perf_Strategy,2)-1

    perfdata=BlocksMean_perf_Strategy(:,j);
    perfdata=perfdata(~isnan(perfdata));

    [~,strategyData]=min(Scores_ByStrategy(:,j));
    jitt=0.07*randn(1,1);
    plot(strategyData+jitt,perfdata(end),'marker','.',...
        'MarkerSize',ms,'Color',[0.8 0.8 0.8])

    Strat_PerfLastBlock(j,:)=[strategyData perfdata(end)];

end

strs=1:6;
for j=1:numel(strs)
    
    if strs(j)==5
    else

    logiS=Strat_PerfLastBlock(:,1)==strs(j);
    P=Strat_PerfLastBlock(logiS,2);    

    errorbar(strs(j),mean(P),std(P),'linewidth',2,'Marker','o',...
        'markersize',ms,'MarkerFaceColor','none','color','k')
    end
end




xlim([0.5 6.5])
ylim([0.45 1.01])
%xlabel('Strategy')
ylabel('Performance (fr. correct)')



% plot letters
lettW=0.03;
lettH=0.06;
yPos=0.945;
xLeft=0.01;
xMiddle=0.4;
xRight=0.68;

annotation(fig,'textbox', [xLeft yPos lettW lettH],'String','A','fontweight','bold',...
    'HorizontalAlignment','center','FontSize',fsizeBig,'FitBoxToText','off','EdgeColor','none',...
    'fontname',fname);
annotation(fig,'textbox', [xMiddle yPos lettW lettH],'String','B','fontweight','bold',...
    'HorizontalAlignment','center','FontSize',fsizeBig,'FitBoxToText','off','EdgeColor','none',...
    'fontname',fname);
 annotation(fig,'textbox', [xRight yPos lettW lettH],'String','C','fontweight','bold',...
     'HorizontalAlignment','center','FontSize',fsizeBig,'FitBoxToText','off','EdgeColor','none',...
     'fontname',fname);


set(fig,'PaperUnits','inches')
set(fig, 'PaperPosition', [0 0 9 4.5])
print(fig,[path2Figs filesep 'FIG_03'],'-dpng','-r450')


%% PLOT v2 (with Strategy CODES)


left=0.1;
bott_up=0.54;
bott_bott=0.08;
w=0.44;
h=0.38;

posits=[left bott_up w h;
    left+1.3*w bott_up 0.68*w h;
    1.15*left bott_bott 0.88*w 0.87*h];

xlimiA=[0.7 6.2];


XTL={'1' '2' ...
    '3' '4' ...
    '5' '6'};

fname='Verdana';
fsize=14;
fsizeBig=fsize+13;
ms=20;
fig=figure;
ax1=axes('Parent',fig,'TickDir','out', 'position', posits(1,:),...
    'XTick',1:6, 'xticklabel',XTL,...
    'Ydir', 'reverse', ...
    'fontsize',fsize,'fontname',fname);
hold(ax1,'on')
%ax.XTickLabelRotation=45;

plot(xlimiA,[0.1 0.1],'linewidth',2,'linestyle',':','Color',[0 0 0])

for j=1:size(Scores_ByStrategy,2)
    jitt=randn(1,1)./20;
    plot((1:6)+jitt,Scores_ByStrategy(:,j),...
        'Marker','.','MarkerSize',ms,'linewidth',0.5)

end

ylabel({'Fit score' '(diff expected-actual)'})
xlabel('Strategy')


% xPos_arrow=0.095;
% %TXT1={'Closer to' 'strategy'};
% TXT1='More fit';
% t1=text(xlimi(1)+0.53,0.08,TXT1,'fontsize',fsize-3);
% t1.Rotation=90;
% annotation(fig,'arrow',[xPos_arrow xPos_arrow],...
%     [0.8 0.9],'headstyle','plain','headwidth',6,'headlength',6);
% %TXT2={'Away from' 'strategy'};
% TXT2='Less fit';
% t2=text(xlimi(1)+0.53,0.28,TXT2,'fontsize',fsize-3);
% t2.Rotation=90;
% annotation(fig,'arrow',[xPos_arrow xPos_arrow],...
%     [0.76 0.65],'headstyle','plain','headwidth',6,'headlength',6);


%xlabel('Strategy')

xlim(xlimiA)


[~,st]=min(Scores_ByStrategy);

ax2=axes('parent',fig,'TickDir','out','position',posits(2,:),...
    'XTick',1:6, 'xticklabel',XTL,...
    'FontSize',fsize,'FontName',fname);
hold(ax2,'on')

H=histogram(st);
H.FaceColor=[0.5 0.5 0];

xlabel('Strategy')
ylabel('Number of participants')


ms=14;

ax3=axes('parent',fig,'TickDir','out','position',posits(3,:),...
    'XTick',1:6,'xticklabel',XTL,...
    'FontSize',fsize,'FontName',fname);
hold(ax3,'on')

lw=1.6;
plot([0.5 6.5],[0.94 0.94],'linewidth',lw,'color','k')
plot([0.5 6.5],[0.45 0.45],'linewidth',lw,'color','k','LineStyle',':')
plot([0.5 6.5],[0.5 0.5],'linewidth',lw,'color','k')
plot([0.5 6.5],[0.55 0.55],'linewidth',lw,'color','k','LineStyle',':')




Strat_PerfLastBlock=nan(1,2);
Strat_PerfDiffFirstLast=nan(1,2);
for j=1:size(BlocksMean_perf_Strategy,2)

    perfdata=BlocksMean_perf_Strategy(:,j);
    perfdata=perfdata(~isnan(perfdata));

    [~,strategyData]=min(Scores_ByStrategy(:,j));
    jitt=0.1*randn(1,1);
    plot(strategyData+jitt,perfdata(end)+(jitt/12),'marker','.',...
       'MarkerSize',ms+2,'Color',[0.65 0.65 0.65])

    Strat_PerfLastBlock(j,:)=[strategyData perfdata(end)];
    Strat_PerfDiffFirstLast(j,:)=[strategyData perfdata(end)-perfdata(1)];

end

% strs=1:6;
% for j=1:numel(strs)
% 
%     if strs(j)==5
%     else
% 
%     logiS=PerfStrat(:,1)==strs(j);
%     P=PerfStrat(logiS,2);    
% 
%     errorbar(strs(j),mean(P),std(P),'linewidth',2,'Marker','o',...
%         'markersize',ms,'MarkerFaceColor','none','color','k')
%     end
% end

text(5,0.96,'Criterion','fontweight','bold',...
    'FontName',fname,'FontSize',fsize-2)

text(1.5,0.52,'Chance','fontweight','bold',...
    'FontName',fname,'FontSize',fsize-2)

xlim([0.5 6.5])
ylim([0.41 1.01])
%xlabel('Strategy')
ylabel({'Last block Performance' '(fraction correct)'})
xlabel('Strategy')


tableContent=cell(1,1);
% Use spaces for alignment and 'bold' marker for header
tableContent{1} = 'Strategies'; 
tableContent{2} = '-------------------------------------------'; % Separator
tableContent{2} = '________________________________'; % Separator

% Data Rows
for i = 1:length(Strategy_Code)
    % Use sprintf to ensure fixed width formatting for alignment:
    % %-5d: left-align the number in 5 characters of width
    row = sprintf('%-5d %s', Strategy_Code(i), Strategy_Description{i});
    tableContent{end+1} = row; %#ok<SAGROW>
end


position_norm = [0.58, 0.1, 0.404, 0.26]; 

annotation('textbox', position_norm, ...
           'String', tableContent, ...
           'Interpreter', 'none', ...          % Ensures text is treated literally
           'FontName', fname, ...       
           'FontSize', fsize, ...
           'BackgroundColor', [0.97 0.97 0.97], ...
           'EdgeColor', 'black', ...
           'FitBoxToText', 'off'); % Adjust the box size to fit the content


% plot letters
lettW=0.03;
lettH=0.06;
yPosTop=0.941;
yPosBott=0.43;
xLeft=0.01;
% xMiddle=0.42;
xRight=0.58;
% 
annotation(fig,'textbox', [xLeft yPosTop lettW lettH],'String','A','fontweight','bold',...
    'HorizontalAlignment','center','FontSize',fsizeBig,'FitBoxToText','off','EdgeColor','none',...
    'fontname',fname);
annotation(fig,'textbox', [xRight yPosTop lettW lettH],'String','B','fontweight','bold',...
    'HorizontalAlignment','center','FontSize',fsizeBig,'FitBoxToText','off','EdgeColor','none',...
    'fontname',fname);
 annotation(fig,'textbox', [xLeft yPosBott lettW lettH],'String','C','fontweight','bold',...
     'HorizontalAlignment','center','FontSize',fsizeBig,'FitBoxToText','off','EdgeColor','none',...
     'fontname',fname);
 annotation(fig,'textbox', [xRight 0.9*yPosBott lettW lettH],'String','D','fontweight','bold',...
     'HorizontalAlignment','center','FontSize',fsizeBig,'FitBoxToText','off','EdgeColor','none',...
     'fontname',fname);


set(fig,'PaperUnits','inches')
set(fig, 'PaperPosition', [0 0 10 8.5])
print(fig,'FIG_03','-dpng','-r450')