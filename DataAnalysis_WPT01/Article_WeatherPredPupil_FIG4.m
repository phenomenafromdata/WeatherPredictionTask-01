%% LOAD DATA, CALCULATE BEHAVIORAL VARIABLES
clearvars;clc
analyze_behav_WeatherPred2018


load('DATA_pupil_BPcenter_zscore_WP2018.mat') %mat file produced by script 'analyze_Pupil_WeatherPred'
%load('DATA_pupil_BPcenter_subtract_WP2018.mat')
PUPILDATA_BP=PUPILDATA;
clear PUPILDATA
load('DATA_pupil_STIMcenter_zscore_WP2018.mat')
%load('DATA_pupil_STIMcenter_subtract_WP2018.mat')
PUPILDATA_STIM=PUPILDATA;
clear PUPILDATA
load('DATA_pupil_FEEDBcenter_zscore_WP2018.mat')
PUPILDATA_FEEDB=PUPILDATA;

disp('***** done loading pupil data')
%% PRODUCE PERI-EVENT PUPIL EPOCHS 
% (requires having loaded PUPILDATA structure array)

baseline_s=0.5;  %duration of baseline, in seconds

pnames_pupil=fieldnames(PUPILDATA_BP.data);
%preallocate
pupilEpochs_BP=nan(1,numel(particip_IDs));
performance_lastBlock=nan(1,numel(particip_IDs));
% counter
c=0;
for thisParticip = 1 : numel(pnames_pupil)
    %check if participant is in behavior database
    id=pnames_pupil{thisParticip}(1:2);
    F=strcmp(id,particip_IDs);
    if sum(F)>0
        c=c+1;
        
        pupilVector=median(PUPILDATA_BP.data.(pnames_pupil{thisParticip}).pupilEpochs_blinkinterp,2,'omitmissing');
        %check Fs and correct if needed
        if PUPILDATA_BP.data.(pnames_pupil{thisParticip}).Fs==1000
            pupilVector=downsample(pupilVector,2);
        end
        pupilEpochs_BP(1:numel(pupilVector),c)=pupilVector;
        performance_lastBlock(c)=mean(PUPILDATA_BP.data.(pnames_pupil{thisParticip}).performance_nomisses(end-49:end));
    end
end

Fs=500;

%tempo=(0:size(pupilEpochs,1)-1)./Fs;
%tempo_EvRelat=tempo+PUPILDATA.info.epoch_before_after_event_s(1);

% calculate pupil maxima times
pnames_pupil=fieldnames(PUPILDATA_BP.data);

pnames_pupil=pnames_pupil(PUPILDATA_BP.info.filter_min_n_trials);

AUC_pattern=nan(numel(pnames_pupil),1);  % area under the curve
AUC_feedback=nan(numel(pnames_pupil),1);
TAB_pattern=nan(numel(pnames_pupil),1); % fraction of time above baseline
TAB_feedback=nan(numel(pnames_pupil),1);

AUC_PATTERN=nan(numel(pnames_pupil),1);  % area under the curve
AUC_FEEDB=nan(numel(pnames_pupil),1);

Maxima_feedb=nan(300,numel(pnames_pupil));
Maxima_patt=nan(300,numel(pnames_pupil));


for thisP = 1:numel(pnames_pupil)

    PM=PUPILDATA_BP.data.(pnames_pupil{thisP}).pupilEpochs_blinkinterp;
    tempo=(0:size(PM,1)-1)./PUPILDATA_BP.data.(pnames_pupil{thisP}).Fs;
    tempo_EvRelat_BP=tempo+PUPILDATA_BP.info.epoch_before_after_event_s(1);
    auc_patt=nan(size(PM,2),1);
    auc_feedb=nan(size(PM,2),1);
    fract_time_aboveBaseline_patt=nan(size(PM,2),1);
    fract_time_aboveBaseline_feedb=nan(size(PM,2),1);
 
    for j=1:size(PM,2) %loop through trials
        pupil=PM(:,j);
        
        pupil=MovMean_DRL(pupil,20);
        
        rt=PUPILDATA_BP.data.(pnames_pupil{thisP}).RT_nomisses(j);
        logi_pattern=tempo_EvRelat_BP<=1 & tempo_EvRelat_BP>-rt;
        logi_feedback=tempo_EvRelat_BP>1;

        logi_beforePattern=tempo_EvRelat_BP<-rt;
        
        %find the index where pattern starts
        t_patternStart=-rt;
        distances=abs(tempo_EvRelat_BP - t_patternStart);
        distances(tempo_EvRelat_BP>=0)=nan; %exclude positive values
        [~, idx_PatternStart] = min(distances);
       
        nSamples_baseline=round(Fs*baseline_s);
        if sum(logi_beforePattern)<nSamples_baseline
            warning('there are not enough samples before pattern to form the baseline')
            disp(['trial = ' num2str(j)])
            disp(pnames_pupil{thisP})
        end
        
        % now build baseline: right before pattern start
        logi_baseline=false(size(tempo_EvRelat_BP));
        logi_baseline(idx_PatternStart-nSamples_baseline:idx_PatternStart)=true;


        auc_patt(j)=trapz(pupil(logi_pattern),tempo(logi_pattern));
        auc_feedb(j)=trapz(pupil(logi_feedback),tempo(logi_feedback));
        
        BL=mean(pupil(logi_baseline));
        
        aboveBL_pattern=sum(pupil(logi_pattern)>BL);
        aboveBL_feedb=sum(pupil(logi_feedback)>BL);

        fract_time_aboveBaseline_patt(j)=aboveBL_pattern./sum(logi_pattern);
        fract_time_aboveBaseline_feedb(j)=aboveBL_feedb./sum(logi_feedback);

        p1=pupil(logi_pattern);
        [max_val1, max_idx1] = max(p1);
        Maxima_patt(j,thisP)=max_val1;
        
        p2=pupil(logi_feedback);
        [max_val2, max_idx2] = max(p2);
        Maxima_feedb(j,thisP)=max_val2;
        
    end
    AUC_pattern(thisP)=mean(auc_patt,'omitmissing');
    AUC_feedback(thisP)=mean(auc_feedb,'omitmissing');
    TAB_pattern(thisP)=mean(fract_time_aboveBaseline_patt);
    TAB_feedback(thisP)=mean(fract_time_aboveBaseline_feedb);
    
    %get block (50-trial) values
    N=50;
    L=numel(auc_patt);
    numBlocks=floor(L/N);
    %auc_temp=reshape(auc_patt,N,numBlocks);
    
    AUC_PATTERN(1:numel(auc_patt),thisP)=auc_patt;
    AUC_FEEDB(1:numel(auc_patt),thisP)=auc_feedb;
    
end

MAX_pattern=mean(Maxima_patt,1,'omitnan');
MAX_feedback=mean(Maxima_feedb,1,'omitnan');

% now for Stim-related 

pnames_pupil=fieldnames(PUPILDATA_STIM.data);
%preallocate
pupilEpochs_STIM=nan(1,numel(particip_IDs));
performance_lastBlock=nan(1,numel(particip_IDs));
% counter
c=0;
for thisParticip = 1 : numel(pnames_pupil)
    %check if participant is in behavior database
    id=pnames_pupil{thisParticip}(1:2);
    F=strcmp(id,particip_IDs);
    if sum(F)>0
        c=c+1;
        
        pupilVector=median(PUPILDATA_STIM.data.(pnames_pupil{thisParticip}).pupilEpochs_blinkinterp,2,'omitmissing');
        %check Fs and correct if needed
        if PUPILDATA_STIM.data.(pnames_pupil{thisParticip}).Fs==1000
            pupilVector=downsample(pupilVector,2);
        end
        pupilEpochs_STIM(1:numel(pupilVector),c)=pupilVector;
        performance_lastBlock(c)=mean(PUPILDATA_STIM.data.(pnames_pupil{thisParticip}).performance_nomisses(end-49:end));
    end
end

tempo=(0 : size(PUPILDATA_STIM.data.aa_54_24.pupilEpochs_blinkinterp)-1)./500;
tempo_EvRelat_STIM=tempo+PUPILDATA_STIM.info.epoch_before_after_event_s(1);



%now for feedback
pnames_pupil=fieldnames(PUPILDATA_FEEDB.data);
%preallocate
pupilEpochs_FEEDB=nan(1,numel(particip_IDs));
performance_lastBlock=nan(1,numel(particip_IDs));
% counter
c=0;
for thisParticip = 1 : numel(pnames_pupil)
    %check if participant is in behavior database
    id=pnames_pupil{thisParticip}(1:2);
    F=strcmp(id,particip_IDs);
    if sum(F)>0
        c=c+1;
        
        pupilVector=median(PUPILDATA_FEEDB.data.(pnames_pupil{thisParticip}).pupilEpochs_blinkinterp,2,'omitmissing');
        %check Fs and correct if needed
        if PUPILDATA_FEEDB.data.(pnames_pupil{thisParticip}).Fs==1000
            pupilVector=downsample(pupilVector,2);
        end
        pupilEpochs_FEEDB(1:numel(pupilVector),c)=pupilVector;
        performance_lastBlock(c)=mean(PUPILDATA_FEEDB.data.(pnames_pupil{thisParticip}).performance_nomisses(end-49:end));
    end
end

tempo=(0 : size(PUPILDATA_FEEDB.data.aa_54_24.pupilEpochs_blinkinterp)-1)./500;
tempo_EvRelat_FEEDB=tempo+PUPILDATA_FEEDB.info.epoch_before_after_event_s(1);


disp('done calculating pupil Epochs')

%% plot - event-related pupil response (BUTTON PRESS)


fsize=13;
fsizeBig=fsize+10;
fname='verdana';

w=0.39;
h=0.42;
left=0.08;
bott_up=0.53;
bott_bott=0.05;


positions=[left bott_up w h;
    left+1.24*w bott_up w h;
    left*2 bott_bott w*0.7 h*0.73;
    left+1.4*w bott_bott w*0.7 h*0.73];


ms=18;

yLAB='Pupil size (z-score)';
ylimi1=[-1.27 1.45];
%ylimi1=[-200 800];
xlimi=[-1.8 3.02];
% incorrText=[-0.76,-0.7];
% corrText=[-0.8,0.6];
% signifposit=-0.8;
% 
% redCol=[0.7 0.2 0.2];
% blueCol=[0.2 0.2 0.7];


boxPlotYwidth=[1 1.18];

fig=figure('color','w');

lw=1;

%%%%%%%%%%% A

ax1=axes('parent',fig, 'position', positions(1,:),'tickdir','out','fontsize',fsize,...
    'fontname', fname,'xtick',-2:3,'xgrid','on');
hold(ax1,'on')


greycol_patt=[0.93 0.93 0.93];
greycol_feedb=[0.82 0.82 0.82];

feedb_area=true(size(tempo_EvRelat_BP));
feedb_area(tempo_EvRelat_BP<1)=false;
feedb_area(tempo_EvRelat_BP>3)=false;

area(tempo_EvRelat_BP,feedb_area*ylimi1(2),'EdgeColor','none','facecolor',greycol_feedb)
area(tempo_EvRelat_BP,feedb_area*(ylimi1(1)*0.99),'EdgeColor','none','facecolor',greycol_feedb)

patt_area=true(size(tempo_EvRelat_BP));
patt_area(tempo_EvRelat_BP<-1.05)=false;
patt_area(tempo_EvRelat_BP>1)=false;

area(tempo_EvRelat_BP,patt_area*ylimi1(2),'EdgeColor','none','facecolor',greycol_patt)
area(tempo_EvRelat_BP,patt_area*(ylimi1(1)*0.99),'EdgeColor','none','facecolor',greycol_patt)

[box,low_whisk,high_whisk]=data4boxplot(-median(RTMat,'omitmissing'));

Y=[zeros(size(box))+boxPlotYwidth(1); zeros(size(box))+boxPlotYwidth(2)];
%plot distribution of onset times
plot([box; box],Y,'linewidth',lw,'color','k')
plot([box(2) box(2)],Y(:,1),'linewidth',lw,'color','r')
plot([box(1) box(3)],[min(min(Y)) min(min(Y))],'linewidth',lw,'color','k')
plot([box(1) box(3)],[max(max(Y)) max(max(Y))],'linewidth',lw,'color','k')
plot([low_whisk box(1)],[mean(mean(Y)) mean(mean(Y))],'linewidth',lw,'color','k')
plot([box(3) high_whisk],[mean(mean(Y)) mean(mean(Y))],'linewidth',lw,'color','k')

annotation(fig,'textbox', [0.04 0.81 0.15 0.15],'String',{'pattern' 'onset'},...
    'HorizontalAlignment','center','FontSize',fsize-4,'FitBoxToText','off','EdgeColor','none',...
    'fontname',fname);

grey_cols=[0.7 0.7 0.7;
    0.6 0.6 0.6;
    0.5 0.5 0.5];

%plot individual participants
for i=1:size(pupilEpochs_BP,2)
    yy=pupilEpochs_BP(:,i);
    yy=MovMean_DRL(yy,35);
    %plot(tempo, yy,'color',[0.35 0.35 0.35])
    
    R=randi([1 3],1,1);
    c=grey_cols(R,:);

    plot(tempo_EvRelat_BP, yy,'color',c)

end
lw=4;


m=mean(pupilEpochs_BP,2,'omitnan');
%cosmetics
win=25;
m=movmean(m,win);

plot(tempo_EvRelat_BP,m,'linewidth',lw+2,'color',[0.99 0.89 0])

plot(xlimi,[0 0], 'k')


xlabel('Time to button press (s)')
ylabel(yLAB)


annotation(fig,'textbox', [0.15 0.84 0.15 0.15],'String','pattern',...
    'HorizontalAlignment','center','FontSize',fsize,'FitBoxToText','off','EdgeColor','none',...
    'fontname',fname);
annotation(fig,'textbox', [0.32 0.84 0.15 0.15],'String','feedback',...
    'HorizontalAlignment','center','FontSize',fsize,'FitBoxToText','off','EdgeColor','none',...
    'fontname',fname);

xlim(xlimi)
ylim(ylimi1)


%%%%%%%%%%%%%%%%%%%  B

ax2=axes('parent',fig,'TickDir','out',...
    'position',positions(2,:),...
    'FontName',fname,'FontSize',fsize);
hold(ax2,'on')

imagesc(tempo_EvRelat_BP,1:size(pupilEpochs_BP,2),pupilEpochs_BP')

plot([1 1],[1 size(pupilEpochs_BP,2)],'linewidth',3,'color','w')

xlabel('Time to button press (s)')
ylabel('Participant')
cb=colorbar;
colormap('jet')
ylabel(cb,'Pupil size (z-score)')

xlim(xlimi)
ylim([1 size(pupilEpochs_BP,2)])


annotation(fig,'textbox', [0.58 0.84 0.15 0.15],'String','pattern',...
    'HorizontalAlignment','center','FontSize',fsize,'FitBoxToText','off','EdgeColor','none',...
    'fontname',fname);
annotation(fig,'textbox', [0.73 0.84 0.15 0.15],'String','feedback',...
    'HorizontalAlignment','center','FontSize',fsize,'FitBoxToText','off','EdgeColor','none',...
    'fontname',fname);



xlimi_CD=[-0.7 1.7];


%%%%%%%%%%%%%  C

ax3=axes('Parent',fig,'TickDir','out', 'position', positions(3,:),...
    'XTick',[0 1],'XTickLabel',{'pattern' 'feedback'},...
    'FontSize',fsize,'FontName',fname);
hold(ax3,'on')
X=randn(size(MAX_pattern))./15;
X2=randn(size(MAX_feedback))./15;

%plot(xlimi_CD,[0 0], ':k')

% plot the pattern period data

plot(X,MAX_pattern,'linestyle','none','marker','.','markersize',ms,'Color',[0.75 0.75 0.75])

[box,low_whisk,high_whisk]=data4boxplot(MAX_pattern);

xRange=[-0.2 0.2];
% 
plot(xRange,[box(1) box(1)],'linewidth',2,'color','k')
plot(xRange,[box(2) box(2)],'linewidth',2,'color','r')
plot(xRange,[box(3) box(3)],'linewidth',2,'color','k')
 
plot([xRange(1) xRange(1)]*0.95,[box(1) box(3)],'linewidth',2,'color','k')
plot([xRange(2) xRange(2)]*0.95,[box(1) box(3)],'linewidth',2,'color','k')
 
plot([mean(xRange) mean(xRange)],[box(1) low_whisk],'linewidth',2,'linestyle','-','color','k')
plot([mean(xRange) mean(xRange)],[box(3) high_whisk],'linewidth',2,'linestyle','-','color','k')

% now the feedback period data

plot((X2+1),MAX_feedback,'linestyle','none','marker','.','markersize',ms,'Color',[0.75 0.75 0.75])

[box,low_whisk,high_whisk]=data4boxplot(MAX_feedback);
xRange=[-0.2 0.2]+1;
% 
plot(xRange,[box(1) box(1)],'linewidth',2,'color','k')
plot(xRange,[box(2) box(2)],'linewidth',2,'color','r')
plot(xRange,[box(3) box(3)],'linewidth',2,'color','k')
 
plot([xRange(1) xRange(1)]*1.015,[box(1) box(3)],'linewidth',2,'color','k')
plot([xRange(2) xRange(2)]*0.9895,[box(1) box(3)],'linewidth',2,'color','k')
 
plot([mean(xRange) mean(xRange)],[box(1) low_whisk],'linewidth',2,'linestyle','-','color','k')
plot([mean(xRange) mean(xRange)],[box(3) high_whisk],'linewidth',2,'linestyle','-','color','k')

ylabel({'Mean pupil maxima' '(z-score)'})

xlim(xlimi_CD)

temp=vertcat(MAX_feedback',MAX_pattern');
ylim([-0.1 1.1*max(temp)])

%%%%%%%% D


ax=axes('Parent',fig,'TickDir','out', 'position', positions(4,:),...
    'XTick',[0 1],'XTickLabel',{'pattern' 'feedback'},...
    'YTick',[0 50 100],...
    'FontSize',fsize,'FontName',fname);
hold(ax,'on')
X=randn(size(TAB_pattern))./15;
X2=randn(size(TAB_feedback))./15;

plot(xlimi_CD,[50 50], ':k')

plot(X,TAB_pattern.*100,'linestyle','none','marker','.','markersize',ms,'Color',[0.75 0.75 0.75])

[box,low_whisk,high_whisk]=data4boxplot(TAB_pattern.*100);
xRange=[-0.2 0.2];
% 
plot(xRange,[box(1) box(1)],'linewidth',2,'color','k')
plot(xRange,[box(2) box(2)],'linewidth',2,'color','r')
plot(xRange,[box(3) box(3)],'linewidth',2,'color','k')
 
plot([xRange(1) xRange(1)]*0.95,[box(1) box(3)],'linewidth',2,'color','k')
plot([xRange(2) xRange(2)]*0.95,[box(1) box(3)],'linewidth',2,'color','k')
 
plot([mean(xRange) mean(xRange)],[box(1) low_whisk],'linewidth',2,'linestyle','-','color','k')
plot([mean(xRange) mean(xRange)],[box(3) high_whisk],'linewidth',2,'linestyle','-','color','k')


plot((X2+1),TAB_feedback.*100,'linestyle','none','marker','.','markersize',ms,'Color',[0.75 0.75 0.75])

[box,low_whisk,high_whisk]=data4boxplot(TAB_feedback.*100);
xRange=[-0.2 0.2]+1;
% 
plot(xRange,[box(1) box(1)],'linewidth',2,'color','k')
plot(xRange,[box(2) box(2)],'linewidth',2,'color','r')
plot(xRange,[box(3) box(3)],'linewidth',2,'color','k')
 
plot([xRange(1) xRange(1)]*1.015,[box(1) box(3)],'linewidth',2,'color','k')
plot([xRange(2) xRange(2)]*0.9895,[box(1) box(3)],'linewidth',2,'color','k')
 
plot([mean(xRange) mean(xRange)],[box(1) low_whisk],'linewidth',2,'linestyle','-','color','k')
plot([mean(xRange) mean(xRange)],[box(3) high_whisk],'linewidth',2,'linestyle','-','color','k')

ylabel({'Mean time pupil spent' 'above baseline (%)'})

xlim(xlimi_CD)

ylim([0 101])

xLeft=0.01;
xRight=0.49;
yPosBottAB=0.94;
yposBottCD=0.39;
lettW=0.03;
lettH=0.06;

annotation(fig,'textbox', [xLeft yPosBottAB lettW lettH],'String','A','fontweight','bold',...
    'HorizontalAlignment','center','FontSize',fsizeBig,'FitBoxToText','off','EdgeColor','none',...
    'fontname',fname);
annotation(fig,'textbox', [xRight yPosBottAB lettW lettH],'String','B','fontweight','bold',...
    'HorizontalAlignment','center','FontSize',fsizeBig,'FitBoxToText','off','EdgeColor','none',...
    'fontname',fname);
annotation(fig,'textbox', [xLeft+0.05 yposBottCD lettW lettH],'String','C','fontweight','bold',...
    'HorizontalAlignment','center','FontSize',fsizeBig,'FitBoxToText','off','EdgeColor','none',...
    'fontname',fname);
annotation(fig,'textbox', [xRight+0.01 yposBottCD lettW lettH],'String','D','fontweight','bold',...
    'HorizontalAlignment','center','FontSize',fsizeBig,'FitBoxToText','off','EdgeColor','none',...
    'fontname',fname);

set(fig,'PaperUnits','inches')
set(fig, 'PaperPosition', [0 0 10.5 7.5])
print(fig,'FIG_04','-dpng','-r450')

%% event-related (PATTERN ONSET)

%compute test (must be set to 'true' only first time)
computePermTest=true;

% type of conf interv to use
ci_type='tdist';  %based on t-distribution
% ci_type='trim';  % based on resampling and trimming data distribution

fsize=13;
fsizeBig=fsize+10;
fname='verdana';

grey_cols=[0.7 0.7 0.7;
    0.6 0.6 0.6;
    0.5 0.5 0.5];


w=0.39;
h=0.72;
left=0.08;
bott_bott=0.15;


positions=[left bott_bott w h;
    left+1.24*w bott_bott w h];


ms=18;

yLAB='Pupil size (z-score)';
ylimi1=[-1.21 1.3];
%xlimi=PUPILDATA_STIM.info.epoch_before_after_event_s;
xlimi=[-1 2.5];
% incorrText=[-0.76,-0.7];
% corrText=[-0.8,0.6];
% signifposit=-0.8;
% 
% redCol=[0.7 0.2 0.2];
% blueCol=[0.2 0.2 0.7];

fig=figure('color','w');

lw=1;

%%%%%%%%%%% A

ax1=axes('parent',fig, 'position', positions(1,:),'tickdir','out','fontsize',fsize,...
    'fontname', fname,'xtick',-2:3,'xgrid','on');
hold(ax1,'on')


%plot individual participants
for i=1:size(pupilEpochs_STIM,2)
    yy=pupilEpochs_STIM(:,i);
    yy=MovMean_DRL(yy,30);
    %plot(tempo, yy,'color',[0.35 0.35 0.35])
    
    R=randi([1 3],1,1);
    c=grey_cols(R,:);

    plot(tempo_EvRelat_STIM, yy,'color',c)

end
lw=4;

x=tempo_EvRelat_STIM;
if strcmpi(ci_type,'tdist')
    ci=GetCIs(pupilEpochs_STIM,0.975,2);
    ci=std(pupilEpochs_STIM,[],2);
    y1=(mean(pupilEpochs_STIM,2,'omitnan')-ci);
    y2=(mean(pupilEpochs_STIM,2,'omitnan')+ci);
    m=mean(pupilEpochs_STIM,2,'omitnan');
    m2=median(pupilEpochs_STIM,2,'omitnan');
elseif strcmpi(ci_type,'trim')
    
    data=reshape(pupilEpochs_STIM,[1 size(pupilEpochs_STIM)]);
    percent=10;
    alphav=0.05;
    
    [~,tmdata,ci]=limo_trimci(data, percent, alphav);
    ci=squeeze(ci);
    y1=ci(:,1);
    y2=ci(:,2);
    m=tmdata;
end



%cosmetics
win=25;
y1=MovMean_DRL(y1,win);
y2=MovMean_DRL(y2,win);
m=MovMean_DRL(m,win);
m2=MovMean_DRL(m2,win);

%fill([x, fliplr(x)] ,[y1; flipud(y2)], [0.88 0.3 0.1],'facealpha',0.8,'linestyle','none')
% plot(tempo,m,'linewidth',lw,'color',[0.5 0 0])
plot(tempo_EvRelat_STIM,m,'linewidth',lw+2,'color',[0.99 0.59 0])
%plot(tempo,m2,'linewidth',lw+2,'color',[0.9 0.8 0])

plot(xlimi,[0 0], 'k')


xlabel('Time to pattern onset (s)')
ylabel(yLAB)


%IM=imread('ButtonPress.png');
% insert button press image
% ax2=axes('parent',fig1, 'position', positions(2,:),'tickdir','out','fontsize',fsize,...
%     'fontname', fname,'ycolor',greycol_patt,'xcolor',greycol_patt,...
%     'yticklabel','','xticklabel','','xgrid','on');
% hold(ax2,'all')
% imagesc(flipud(IM))
% xlim([1 size(IM,2)])
% ylim([1 size(IM,1)])
% axis xy


xlim(xlimi)
ylim(ylimi1)


% [performance_lastBlock_sorted,idcs]=sort(performance_lastBlock);
% pupilEpochs=pupilEpochs(:,idcs);


%%%%%%%%%%%%%%%%%%%  B

ax2=axes('parent',fig,'TickDir','out',...
    'position',positions(2,:),...
    'FontName',fname,'FontSize',fsize);
hold(ax2,'on')

imagesc(tempo_EvRelat_STIM,1:size(pupilEpochs_STIM,2),pupilEpochs_STIM')

plot([0 0],[1 size(pupilEpochs_STIM,2)],'linewidth',3,'color','w')

xlabel('Time to pattern onset (s)')
ylabel('Participant')
cb=colorbar;
colormap('jet')
ylabel(cb,'Pupil size (z-score)')

xlim(xlimi)
ylim([1 size(pupilEpochs_STIM,2)])

set(fig,'PaperUnits','inches')
set(fig, 'PaperPosition', [0 0 10.5 5.5])
print(fig,[path2Figs filesep 'FIG_04part2'],'-dpng','-r450')

%% event-related  FEEDBACK ONSET


fsize=13;
fsizeBig=fsize+10;
fname='verdana';

grey_cols=[0.7 0.7 0.7;
    0.6 0.6 0.6;
    0.5 0.5 0.5];


w=0.39;
h=0.72;
left=0.08;
bott_bott=0.15;


positions=[left bott_bott w h;
    left+1.24*w bott_bott w h];


ms=18;

yLAB='Pupil size (z-score)';
ylimi1=[-1.21 1.3];
%xlimi=PUPILDATA_STIM.info.epoch_before_after_event_s;
xlimi=[-1 2.05];
% incorrText=[-0.76,-0.7];
% corrText=[-0.8,0.6];
% signifposit=-0.8;
% 
% redCol=[0.7 0.2 0.2];
% blueCol=[0.2 0.2 0.7];

fig=figure('color','w');

lw=1;

%%%%%%%%%%% A

ax1=axes('parent',fig, 'position', positions(1,:),'tickdir','out','fontsize',fsize,...
    'fontname', fname,'xtick',-2:3,'xgrid','on');
hold(ax1,'on')


%plot individual participants
for i=1:size(pupilEpochs_FEEDB,2)
    yy=pupilEpochs_FEEDB(:,i);
    yy=MovMean_DRL(yy,30);
    %plot(tempo, yy,'color',[0.35 0.35 0.35])
    
    R=randi([1 3],1,1);
    c=grey_cols(R,:);

    plot(tempo_EvRelat_FEEDB, yy,'color',c)

end
lw=4;

x=tempo_EvRelat_FEEDB;
m=mean(pupilEpochs_FEEDB,2,'omitnan');



%cosmetics
win=25;
m=MovMean_DRL(m,win);

%fill([x, fliplr(x)] ,[y1; flipud(y2)], [0.88 0.3 0.1],'facealpha',0.8,'linestyle','none')
% plot(tempo,m,'linewidth',lw,'color',[0.5 0 0])
plot(tempo_EvRelat_FEEDB,m,'linewidth',lw+2,'color',[0.99 0.59 0])
%plot(tempo,m2,'linewidth',lw+2,'color',[0.9 0.8 0])

plot(xlimi,[0 0], 'k')


xlabel('Time to feedback onset (s)')
ylabel(yLAB)


%IM=imread('ButtonPress.png');
% insert button press image
% ax2=axes('parent',fig1, 'position', positions(2,:),'tickdir','out','fontsize',fsize,...
%     'fontname', fname,'ycolor',greycol_patt,'xcolor',greycol_patt,...
%     'yticklabel','','xticklabel','','xgrid','on');
% hold(ax2,'all')
% imagesc(flipud(IM))
% xlim([1 size(IM,2)])
% ylim([1 size(IM,1)])
% axis xy


xlim(xlimi)
ylim(ylimi1)


% [performance_lastBlock_sorted,idcs]=sort(performance_lastBlock);
% pupilEpochs=pupilEpochs(:,idcs);


%%%%%%%%%%%%%%%%%%%  B

ax2=axes('parent',fig,'TickDir','out',...
    'position',positions(2,:),...
    'FontName',fname,'FontSize',fsize);
hold(ax2,'on')

imagesc(tempo_EvRelat_FEEDB,1:size(pupilEpochs_FEEDB,2),pupilEpochs_FEEDB')

plot([0 0],[1 size(pupilEpochs_FEEDB,2)],'linewidth',3,'color','w')

xlabel('Time to feedback onset (s)')
ylabel('Participant')
cb=colorbar;
colormap('jet')
ylabel(cb,'Pupil size (z-score)')

xlim(xlimi)
ylim([1 size(pupilEpochs_FEEDB,2)])

set(fig,'PaperUnits','inches')
set(fig, 'PaperPosition', [0 0 10.5 5.5])
print(fig,[path2Figs filesep 'FIG_04part3'],'-dpng','-r450')



%%

figure;
hold on
for j=1:size(AUC_FEEDB,2)
    y=AUC_FEEDB(:,j);
    y=y(~isnan(y));
    x=1:numel(y);

    [B,B0]=LeastSquaresLine(x,y);

    plot(x,y,'linestyle','none','marker','.',...
        'color',[0.65 0.65 0.65]);
    plot([min(x) max(x)],B0+B.*[min(x) max(x)],'linewidth',2)
end

title('feedback')


figure;
hold on
for j=1:size(AUC_PATTERN,2)
    y=AUC_PATTERN(:,j);
    y=y(~isnan(y));
    x=1:numel(y);

    [B,B0]=LeastSquaresLine(x,y);

    plot(x,y,'linestyle','none','marker','.',...
        'color',[0.65 0.65 0.65]);
    plot([min(x) max(x)],B0+B.*[min(x) max(x)],'linewidth',2)
end

title('pattern')
