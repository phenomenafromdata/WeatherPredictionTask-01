function [ERR,t_crit]=GetCIs(datamatrix,TwoTailedSignifLevel,dimension)

% this function calculates confidence limits around the mean of a sample.
% confidence limits are based on a t-distribution, with significance level
% determined by 'TwoTailedSignifLevel', which should be a number between 0 and 1.
% e.g., alpha = 2*(1-TwoTailedSignifLevel)

% since the input data is a matrix, 'dimension' specifies the dimension
% along which the mean is calculates (dimension =1 produces a row vector
% whose entries are the mean of each column; dimension =2 produces a column
% vector whose entries are the mean of each row).

%the output 'ERR' is a vector containing the values t[alpha,df]*(s/sqrt(n)) for each
%row or column, where t = t units corresponding to alpha and degrees of
%freedom df, s= standard deviation of the row/column and n = sample size of
%the row/column.

% by DRL, 2009

dim1=size(datamatrix,1);
dim2=size(datamatrix,2);

if dim1==1 && dimension==1
    ERR=zeros(size(datamatrix)); t_crit=NaN;
    disp('Can''t compute conf interv b/c only one sample')
elseif dim2==1 && dimension==2
    ERR=zeros(size(datamatrix)); t_crit=NaN;
    disp('Can''t compute conf interv b/c only one sample')

else
   

% first, loop to figure out how many of the observations are NaNs (if there are some)

A=ones(dim1,dim2).*-1;
for i=1:dim2
    xm=isnan(datamatrix(:,i));
    A(:,i)=logical(xm-1);
end

%A is a logical matrix containing ones marking non-nan in the corresponding
%entries of 'datamatrix'.




if dimension==1
    Ns=sum(A,1);  %efective sample size for each column
    t_crit=zeros(1,dim2);    %preallocating
    ERR=zeros(1,dim2);
    forloopcycle=dim2;
    sdatamatrix=nanstd(datamatrix);
    
    for j=1:forloopcycle
        t=tinv(TwoTailedSignifLevel, Ns(1,j)-1);
        t_crit(1,j)=t;
    end
    
    for i=1:forloopcycle
            ERR(1,i)=(t_crit(1,i).*sdatamatrix(1,i)./(sqrt(Ns(1,i))));
    end
    
    
    
elseif dimension==2
    Ns=sum(A,2);   %efective sample size for each row
    t_crit=zeros(dim1,1);
    ERR=zeros(dim1,1);
    forloopcycle=dim1;
    sdatamatrix=nanstd(datamatrix');sdatamatrix=sdatamatrix';
    
    for j=1:forloopcycle
        t=tinv(TwoTailedSignifLevel, Ns(j,1)-1);
        %disp(Ns(j,1))
        t_crit(j,1)=t;
    end
    
    
    for i=1:forloopcycle
        ERR(i,1)=(t_crit(i,1).*sdatamatrix(i,1)./(sqrt(Ns(i,1))));
    end
    
end

end
