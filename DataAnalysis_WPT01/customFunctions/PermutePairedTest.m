function [pval, stat_no_perm, dof, stat_perm]=PermutePairedTest(data1, data2, n_perms, statistic, tails)

% inputs:
% data1 and data2 must be column vectors with non- paired data (eg two exp conditions)
% n_perms is the number of permutations
% statistic is the test to use: 'tstat', 'dCohen'
% tails is a string can be 'both', 'left' or 'right'.

%outputs:
%pval: p-value of test
%stat_no_perm: value of test statistic for the actual data comparison
%stat_perm: collection of values for test statistic from all permutations


% test is based on Maris & Oostenveld 2007 (doi:10.1016/j.jneumeth.2007.03.024)
%summary of procedure:
% (1) Collect the trials of the two experimental conditions in a
% single set.
% (2) Randomly draw as many trials from this combined data set
% as there were trials in condition 1 and place those trials into
% subset 1. Place the remaining trials in subset 2. The result
% of this procedure is called a random partition.
% (3) Calculate the test statistic on this random partition.
% (4) Repeat steps 2 and 3 a large number of times and construct
% a histogram of the test statistics.
% (5) From the test statistic that was actually observed and the
% histogram in step 4, calculate the proportion of random partitions
% that resulted in a larger test statistic than the observed
% one. This proportion is called the p-value.


% by Daniel Rojas Libano, July 2019

% 04/07/2023 Changed function to match 'PermuteNonPairedTest_v2'
% March 2024: Added degrees of freedom (dof) as output

%%
if ~isequal(numel(data1), numel(data2))
    error('data1 and data2 must have same dimensions')
end

% 1) calculate non-permuted statistic

if strcmp(statistic, 'tstat')
%     stat_no_perm=TTEST_DRL(data1,data2); %t-value
    [~,~,~,stats]=ttest(data1,data2); %t-value
    stat_no_perm=stats.tstat;
    dof=stats.df;
elseif strcmp(statistic,'dCohen')
    d=nanmean(data1)-nanmean(data2);
%     s=nanstd(data1);
    s=PooledStd(data1,data2);
    stat_no_perm=d/s;
end

% 2) Perform permutations 

stat_perm=nan(n_perms,1); %preallocate
for thisPermutation=1:n_perms
    %construct random partition
    both=[data1; data2];
    r=randperm(numel(both));
    both=both(r);
    
    if strcmp(statistic, 'tstat')
        [~,~,~,d]=ttest(both(1:numel(data1)),both(numel(data1)+1:end));
        stat_perm(thisPermutation)=d.tstat;
    elseif strcmp(statistic,'dCohen')
        d=nanmean(both(1:numel(data1))-nanmean(both(numel(data1)+1:end)));
%         s=nanstd(data1);
        s=PooledStd(data1,data2);
        stat_perm(thisPermutation)=d/s;
    end
end

% 3) calculate p-value

% proportion of random partitions that resulted in a larger
% test statistic than the observed one. 

if strcmp(tails,'right')
    logi4pvalue=stat_perm>stat_no_perm;
elseif strcmp(tails,'left')
    logi4pvalue=stat_perm<stat_no_perm;
elseif strcmp(tails,'both')
    if stat_no_perm<0
        logi4pvalue=stat_perm<stat_no_perm | stat_perm>-stat_no_perm;
    else
        logi4pvalue=stat_perm<-stat_no_perm | stat_perm>stat_no_perm;
    end
end

pval=sum(logi4pvalue)./n_perms;