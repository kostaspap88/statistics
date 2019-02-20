%Description: Higher order univariate t-test computation using formula from 
%"Leakage Assessment Methodology, Schneider, Moradi". See formulas (5,6,7),
%section 4.2. Computing the t-test threshold value uses "Towards Sound and
%Optimal Leakage Detection Procedure, Zhang et al." which is based on Sidak
%correction


%Inputs: the data matrix, the order of the t-test, a function to compute
%central moments and the significance level alpha
%Output: the t test values, the threshold values (for normal and student
%distribution) and the central moments

function [t, threshold_normal, threshold_student, cm1, mu1, cm2, mu2] = ttest_ho_univariate(data1, data2, tt_order,cm_method, alpha)


%number of traces in data matrix 1
nA=size(data1,1);
%number of traces in data matrix 2
nB=size(data2,1);

%number of samples in data matrix 1 and data matrix 2 (same number)
no_samples = size(data1,2);

%init
mu1=zeros(1,no_samples);
mu2=zeros(1,no_samples);
cm1=zeros(no_samples, 2*tt_order);
cm2=zeros(no_samples, 2*tt_order);

%compute the moments for every sample in the traceset
for i=1:no_samples
    [cm_d1, mu_d1]=cm_method(data1(:,i),2*tt_order);
    mu1(i)=mu_d1;
    cm1(i,:)=cell2mat(cm_d1);
    [cm_d2, mu_d2]=cm_method(data2(:,i),2*tt_order);
    mu2(i)=mu_d2;
    cm2(i,:)=cell2mat(cm_d2);
end


switch tt_order
    case 1
        mA=mu1;
        mB=mu2;
        varA=cm1(:,2);
        varB=cm2(:,2);
    case 2      
        mA=cm1(:,2);
        mB=cm2(:,2);
        varA=cm1(:,4)-cm1(:,2).^2;
        varB=cm2(:,4)-cm2(:,2).^2;
    otherwise
        mA= cm1(:, tt_order)/(sqrt(cm1(:, 2)).^tt_order);
        mB= cm2(:, tt_order)/(sqrt(cm2(:, 2)).^tt_order);
        varA=(cm1(:,2*tt_order)-cm1(:,tt_order).^2)/(cm1(:,2).^tt_order);
        varB=(cm2(:,2*tt_order)-cm2(:,tt_order).^2)/(cm2(:,2).^tt_order);
end

%t-test value computation
t = ( mA - mB )./ sqrt(varA'/nA+varB'/nB); 

%compute the threhold for the test

%sA, sB is the sampled std deviation, needs to be estimated from dataset
sA = sqrt(cm1(:,2)*nA/(nA-1)); 
sB = sqrt(cm2(:,2)*nB/(nB-1));

%family-wise error rate
fwer = 1-(1-alpha)^no_samples ;
%Sidak correction
sidak_a = 1-(1-alpha)^(1/no_samples) ;

%compute the degrees of freedom
df=zeros(1,no_samples);
threshold_student=zeros(1,no_samples);
for i=1:no_samples
    df(i)=( (sA(i)^2)/nA + (sB(i)^2)/nB )^2 /  ( ((sA(i)^2)/nA)^2 /(nA-1) + ((sB(i)^2)/nB)^2 /(nB-1) );
    %threshold value for student t distribution
    threshold_student(i)=tinv(1-sidak_a/2,df(i));
end

%threshold value for normal(0,1)
threshold_normal=norminv(1-sidak_a/2);


end