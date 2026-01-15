%% LOAD DATA
%clearvars
% requires having run "analyze_LinModels_eye_Behavior_WeatherPred.m"

%load data
load('LinMod_PUPIL_WP_BPcenter_Xtrial_Ypupil.mat')
load('LinMod_PUPIL_NULL_WP_BPcenter_Xtrial_Ypupil.mat')

%load image


tempo=(0:size(LM_PUPIL.data.betaMat,1)-1)./LM_PUPIL.data.Fs;
tempo=tempo + LM_PUPIL.info.epoch_before_after_event_s(1);

I2_LM=imread('LM2_LM_v2.png');

%crop image (otherwise shows grey borders)
border_size = 10; % in pixels

start_index = border_size + 1;
end_index = size(I2_LM, 1) - border_size; % size(I2_LM, 1) is the original height (rows)
end_index_cols = size(I2_LM, 2) - border_size; % size(I2_LM, 2) is the original width (columns)

I2_LM = I2_LM(start_index:end_index, start_index:end_index_cols, :);

I2_LM(:,end-70:end,:)=[];


%% plot figure

left=-0.04;
bott_mid=0.55;
bott_bott1=0.18;
bott_bott2=0.1;
w=0.64;
h=0.64;
w2=0.33;
h2=0.39;

left2=0.64;

positions=[left bott_bott1 w h;
    left2 bott_mid w2 h2;
    left2 bott_bott2 w2 h2];

ylimiBeta=[-0.35e-2 0.17e-2];

ylimiR2=[0 0.06];
xlimi=[-1.5 3];

ms=4;

Betas_area_col=[0.5 0.88 0.5];
Betas_line_col=[0.15 0.6 0.15];

R2_area_col=Betas_area_col;
R2_line_col=Betas_line_col;

yval_significanceBar=-0.0033;

Null_area_col=[0.88 0.88 0.5];
Null_line_col=[0.4 0.4 0];

fname='Verdana';
fsize=14;
fsizeBig=fsize+11;

fig=figure;

ax1=axes('Parent',fig,'Position',positions(1,:),'Visible','off');
hold(ax1,'on')

imshow(I2_LM)

ax2=axes('Parent',fig,'TickDir','out','position',positions(2,:),...
    'FontName',fname,'FontSize',fsize);
hold(ax2,'on')
greycol_patt=[0.93 0.93 0.93];
greycol_feedb=[0.82 0.82 0.82];

feedb_area=true(size(tempo));
feedb_area(tempo<1)=false;
feedb_area(tempo>3)=false;

area(tempo,feedb_area*ylimiBeta(2),'EdgeColor','none','facecolor',greycol_feedb)
area(tempo,feedb_area*(ylimiBeta(1)*0.99),'EdgeColor','none','facecolor',greycol_feedb)

patt_area=true(size(tempo));
patt_area(tempo<-1.05)=false;
patt_area(tempo>1)=false;

area(tempo,patt_area*ylimiBeta(2),'EdgeColor','none','facecolor',greycol_patt)
area(tempo,patt_area*(ylimiBeta(1)*0.99),'EdgeColor','none','facecolor',greycol_patt)

grey_cols=[0.7 0.7 0.7;
    0.6 0.6 0.6;
    0.5 0.5 0.5];




factor4visualiz=1;

% null
Mat2plot=LM_PUPIL_NULL.data.betaMat(:,LM_PUPIL.info.filter_min_n_trials).*factor4visualiz;

ci=GetCIs(Mat2plot,0.975,2).*5;
y1=(mean(Mat2plot,2,'omitnan')-ci);
y2=(mean(Mat2plot,2,'omitnan')+ci);

m=mean(Mat2plot,2,'omitnan');


win=25;
y1=MovMean_DRL(y1,win);
y2=MovMean_DRL(y2,win);
m=MovMean_DRL(m,win);

fill([tempo, fliplr(tempo)] ,[y1; flipud(y2)], Null_area_col,'facealpha',0.6,'linestyle','none')
plot(tempo, m,'linewidth',4,...
    'Color',Null_line_col)


% data
Mat2plot=LM_PUPIL.data.betaMat(:,LM_PUPIL.info.filter_min_n_trials).*factor4visualiz;

ci=GetCIs(Mat2plot,0.975,2);
y1=(mean(Mat2plot,2,'omitnan')-ci);
y2=(mean(Mat2plot,2,'omitnan')+ci);

m=mean(Mat2plot,2,'omitnan');


win=25;
y1=MovMean_DRL(y1,win);
y2=MovMean_DRL(y2,win);
m=MovMean_DRL(m,win);

fill([tempo, fliplr(tempo)] ,[y1; flipud(y2)], Betas_area_col,'facealpha',0.6,'linestyle','none')
plot(tempo, m,'linewidth',4,...
    'Color',Betas_line_col)


ax2.YAxis.ExponentMode='manual';  %so it does not use scientific notation

ylim(ylimiBeta)
xlim(xlimi)

%xlabel('Time to button press (s)')
ylabel('\beta values (z-unit/trial)')

annotation(fig,'textbox', [0.62 0.88 0.2 0.1],'String','pattern',...
    'HorizontalAlignment','center','FontSize',fsize,'FitBoxToText','off','EdgeColor','none',...
    'fontname',fname);
annotation(fig,'textbox', [0.8 0.88 0.2 0.1],'String','feedback',...
    'HorizontalAlignment','center','FontSize',fsize,'FitBoxToText','off','EdgeColor','none',...
    'fontname',fname);

text(1.7,0.05e-2, 'null','FontSize',fsize,'FontName',fname,...
    'Color',[0.5 0.5 0],'FontWeight','bold');
text(1.2,-0.24e-2, 'data','FontSize',fsize,'FontName',fname,...
    'Color',[0 0.5 0],'FontWeight','bold');

% perform ttest

Mat_NULL=LM_PUPIL_NULL.data.betaMat;
Mat_DATA=LM_PUPIL.data.betaMat;

if exist('final_corrected_mask','var')

else % compute cluster-based permutations test 

    final_alpha = 0.05;

    n_perms=1000;
    
    % 1. Get the T-stats from your original test
    disp('computing permutations ttest for beta values...')
    pvals_permtest=nan(size(Mat_DATA,1),1);
    tstats_permtest=nan(size(Mat_DATA,1),1);
    hedgesg=nan(size(Mat_DATA,1),1);
    for thisTime=1:size(Mat_DATA,1) % loop through timepoints

        [pv, tv] = PermutePairedTest(Mat_DATA(thisTime,:),Mat_NULL(thisTime,:),n_perms,...
            'tstat','both');

        pvals_permtest(thisTime)=pv;

        tstats_permtest(thisTime)=tv;
        
        temp=mes(Mat_DATA(thisTime,:)',Mat_NULL(thisTime,:)','hedgesg');

        hedgesg(thisTime,1) = temp.hedgesg;
        hedgesg(thisTime,2:3) = temp.hedgesgCi'; 

    end
    disp('permutations done.')
    

    % 2. Identify streaks of significance
    signif_mask = pvals_permtest < final_alpha;
    labels = cumsum([signif_mask(1); diff(signif_mask) == 1] .* signif_mask);
    % This creates a vector like [0 0 1 1 1 0 0 2 2 0] where each cluster has a unique ID

    % 3. Calculate the Mass for each cluster
    num_clusters = max(labels);
    cluster_masses = zeros(num_clusters, 1);

    for c = 1:num_clusters
        cluster_indices = find(labels == c);
        % Sum the T-statistics for this cluster (the "Mass")
        cluster_masses(c) = sum(abs(tstats_permtest(cluster_indices)));
    end


    % --- 4. Build Null Distribution for Cluster Correction ---
    n_null_perms = 1000; % 1000 is standard for the cluster-level test
    max_null_mass = zeros(n_null_perms, 1);
    diff_mat = Mat_DATA - Mat_NULL; % The pairwise difference

    disp('Building null distribution of cluster masses...')
    for p = 1:n_null_perms
        % Shuffle: flip the sign of the difference for each participant randomly
        flip_mask = sign(randn(1, size(Mat_DATA, 2)));
        perm_diff = diff_mat .* flip_mask;

        % Fast T-test across all timepoints for this shuffle
        t_perm = mean(perm_diff, 2) ./ (std(perm_diff, 0, 2) ./ sqrt(size(Mat_DATA, 2)));

        % Find clusters in this random shuffle (using your 0.01 threshold)
        perm_mask = abs(t_perm) > tinv(1 - 0.01/2, size(Mat_DATA, 2)-1);

        % Manual clumping for the shuffle
        perm_labels = cumsum([perm_mask(1); diff(perm_mask) == 1] .* perm_mask);
        num_perm_clusters = max(perm_labels);

        if num_perm_clusters > 0
            % Calculate masses for all clusters in this shuffle and take the MAX
            temp_masses = zeros(num_perm_clusters, 1);
            for pc = 1:num_perm_clusters
                temp_masses(pc) = sum(abs(t_perm(perm_labels == pc)));
            end
            max_null_mass(p) = max(temp_masses);
        else
            max_null_mass(p) = 0;
        end
    end

    % --- 5. Calculate Cluster P-values and Create Final Mask ---
    cluster_pvals = zeros(num_clusters, 1);
    for c = 1:num_clusters
        % P-value is the proportion of null masses >= the observed mass
        cluster_pvals(c) = mean(max_null_mass >= cluster_masses(c));
    end

    % Set your final alpha (usually 0.05)

    sig_cluster_ids = find(cluster_pvals < final_alpha);

    % Create the final logical vector
    final_corrected_mask = ismember(labels, sig_cluster_ids);

    disp('DONE!!')
end


dumm=ones(size(tempo)).*yval_significanceBar;

plot(tempo(final_corrected_mask), dumm(final_corrected_mask),'Marker','square',...
    'linestyle','none',...
    'markersize',ms,'markerfacecolor',[0 0 0],'color',[0 0 0])

text(0,yval_significanceBar*0.95,'*','fontname',fname,'fontsize',fsize)
text(2.2,yval_significanceBar*0.95,'*','fontname',fname,'fontsize',fsize)

ylim(ylimiBeta)

% plot R squared
ax3=axes('Parent',fig,'TickDir','out','position',positions(3,:),...
    'FontName',fname,'FontSize',fsize);
hold(ax3,'on')
greycol_patt=[0.93 0.93 0.93];
greycol_feedb=[0.82 0.82 0.82];

feedb_area=true(size(tempo));
feedb_area(tempo<1)=false;
feedb_area(tempo>3)=false;

area(tempo,feedb_area*ylimiR2(2),'EdgeColor','none','facecolor',greycol_feedb)
area(tempo,feedb_area*(ylimiR2(1)*0.99),'EdgeColor','none','facecolor',greycol_feedb)

patt_area=true(size(tempo));
patt_area(tempo<-1.05)=false;
patt_area(tempo>1)=false;

area(tempo,patt_area*ylimiR2(2),'EdgeColor','none','facecolor',greycol_patt)
area(tempo,patt_area*(ylimiR2(1)*0.99),'EdgeColor','none','facecolor',greycol_patt)

factor4visualiz=1;

% null
Mat2plot=LM_PUPIL_NULL.data.r2Mat(:,LM_PUPIL.info.filter_min_n_trials).*factor4visualiz;

ci=GetCIs(Mat2plot,0.975,2).*2;
y1=(mean(Mat2plot,2,'omitnan')-ci);
y2=(mean(Mat2plot,2,'omitnan')+ci);

m=mean(Mat2plot,2,'omitnan');
 
 
win=25;
y1=movmean(y1,win);
y2=movmean(y2,win);
m=movmean(m,win);
% 
fill([tempo, fliplr(tempo)] ,[y1; flipud(y2)], Null_area_col,'facealpha',0.6,'linestyle','none')
plot(tempo, m,'linewidth',4,...
     'Color',Null_line_col)


%data
Mat2plot=LM_PUPIL.data.r2Mat(:,LM_PUPIL.info.filter_min_n_trials).*factor4visualiz;

ci=GetCIs(Mat2plot,0.975,2);
y1=(mean(Mat2plot,2,'omitnan')-ci);
y2=(mean(Mat2plot,2,'omitnan')+ci);

m=mean(Mat2plot,2,'omitnan');


win=25;
y1=movmean(y1,win);
y2=movmean(y2,win);
m=movmean(m,win);

fill([tempo, fliplr(tempo)] ,[y1; flipud(y2)], R2_area_col,'facealpha',0.6,'linestyle','none')


plot(tempo, m,'linewidth',4,...
    'Color',R2_line_col)

plot(xlimi,[0 0], 'k')


xlim(xlimi)

ylabel('R^2')
xlabel('Time to button press (s)')

ylim(ylimiR2)

text(1.7,0.01, 'null','FontSize',fsize,'FontName',fname,...
    'Color',[0.5 0.5 0],'FontWeight','bold');
text(-0.1,0.04, 'data','FontSize',fsize,'FontName',fname,...
    'Color',[0 0.5 0],'FontWeight','bold');




xLeft=0.01;
xRight=0.53;
yPosBott=0.92;
lettW=0.03;
lettH=0.06;



annotation(fig,'textbox', [xLeft 0.95*yPosBott lettW lettH],'String','A','fontweight','bold',...
    'HorizontalAlignment','center','FontSize',fsizeBig,'FitBoxToText','off','EdgeColor','none',...
    'fontname',fname);
annotation(fig,'textbox', [0.98*xRight yPosBott lettW lettH],'String','B','fontweight','bold',...
    'HorizontalAlignment','center','FontSize',fsizeBig,'FitBoxToText','off','EdgeColor','none',...
    'fontname',fname);
annotation(fig,'textbox', [xRight 0.5*yPosBott lettW lettH],'String','C','fontweight','bold',...
    'HorizontalAlignment','center','FontSize',fsizeBig,'FitBoxToText','off','EdgeColor','none',...
    'fontname',fname);


set(fig,'PaperUnits','inches')
set(fig, 'PaperPosition', [0 0 10.5 7.5])
print(fig,'FIG_05','-dpng','-r450')
