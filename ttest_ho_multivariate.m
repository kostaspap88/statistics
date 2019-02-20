%Description: Higher order multivariate t-test computation using formula from 
%"Leakage Assessment Methodology, Schneider, Moradi". See formulas (10,12,13)
%Computing the t-test threshold value uses "Towards Sound and
%Optimal Leakage Detection Procedure, Zhang et al." which is based on Sidak
%correction


%Inputs: the data matrix, the order of the t-test, a function to compute
%multivariate central moments and the significance level alpha
%Output: the t test values, the threshold values (normal and student), the central moments

function [t, threshold_normal, threshold_student, cm_d1, mu_d1, cm_d2, mu_d2] = ttest_ho_multivariate(data1, data2, tt_order,cm_method, alpha)
    

%number of traces in data matrix 1
nA=size(data1,1);
%number of traces in data matrix 2
nB=size(data2,1);

%number of samples in data matrix 1 and data matrix 2 (same number)
no_samples = size(data1,2);

%compute the multivariate moments 
[cm_d1, mu_d1]=cm_method(data1,2*tt_order);
[cm_d2, mu_d2]=cm_method(data2,2*tt_order);

denominator_index1 = 1;
denominator_index2 = 1;
numerator_index1 = 2;
numerator_index2 = 2;


%t-test computation

%find entries with [1,2,3,...,no_poi]
mA=cm_d1{tt_order,numerator_index1};
mB=cm_d2{tt_order,numerator_index2};

%find entries with [1,2,3,...,no_poi, 1,2,3,...no_poi]
varA=cm_d1{2*tt_order,denominator_index1} - mA^2;
varB=cm_d2{2*tt_order,denominator_index2} - mB^2;

%t-test value computation
t = ( mA - mB )/ sqrt(varA/nA+varB/nB); 

%compute the threhold for the test

%sA, sB is the sampled std deviation, needs to be estimated from dataset
for i=1:size(cm_d1(2,:),2)
    temp1(i)=cm_d1{2,i};
    temp2(i)=cm_d2{2,i};
end
sA = sqrt(temp1*nA/(nA-1)); 
sB = sqrt(temp2*nB/(nB-1));

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