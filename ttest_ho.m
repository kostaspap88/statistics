%Description: Higher order t-test computation using formula from 
%"Leakage Assessment Methodology, Schneider, Moradi". See formulas (5,6,7),
%section 4.2



%Inputs: the data matrix, the order of the t-test and a function to compute
%central moments
%Output: the t test values and the threshold value

function [t] = ttest_ho(data1, data2, tt_order,cm_method)

%number of traces in data matrix 1
nA=size(data1,1);
%number of traces in data matrix 2
nB=size(data2,1);


%number of samples in data matrix 1 and data matrix 2
no_samples = size(data1,2);

mu1=zeros(1,no_samples);
mu2=zeros(1,no_samples);
cm1=zeros(no_samples, 2*tt_order);
cm2=zeros(no_samples, 2*tt_order);

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

%t-test computation
t = ( mA - mB )./ power((varA'/nA+varB'/nB),0.5); 



end