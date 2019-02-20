%main file of High Performance Statistics 
%author: Kostas - kostaspap88@gmail.com - kpcrypto.net

%CAROLINE DELETED
%GOODBYE...CAROLINE

clear all
close all

%the main function tests the equivalence and numerical stability of
%different computational variants for central moments and other statistics

format longG



%TESTING UNIVARIATE CENTRAL MOMENTS----------------------------------------

%univariate dataset simulation
data_univ=normrnd(9,2,1000,1);

%specify central moments order
order_univ=4;

%test the standard matlab formula 
%(no clue what the official MATLAB implementation is)
cm_matlab=cell(1,order_univ);
for i=1:order_univ
    cm_matlab{i}=moment(data_univ,i);   
end
mu_matlab = mean(data_univ);

tic
%test the corrected 2 pass formula (serialized)
[cm_2pcorr_serial, mu_2pcorr_serial]=cm_2passcorrected_serial(data_univ,order_univ);
toc

%test the corrected 2 pass formula (vectorized)
[cm_2pcorr_vectorized, mu_2pcorr_vectorized]=cm_2passcorrected_vectorized(data_univ,order_univ);

%test the naive 2 pass formla (serialized) - just for reference
[cm_2pnaive_serial, mu_2pnaive_serial]=cm_2passnaive_serial(data_univ,order_univ);

%test the univariate singleton 1 pass formula (serialized by construction)
[cm_1psingle, mu_1psingle]=cm_1pass_univariate_singleton(data_univ,order_univ);

%--------------------------------------------------------------------------



%TESTING MULTIVARIATE CENTRAL MOMENTS--------------------------------------
    
%multivariate dataset simulation with POIs
data_mv=normrnd(9,2,50,2);

%specify central moments order
order_mv=4;

%test the multivariate singleton 1 pass formula (serialized by construction)
[cm_mv_1psingle, mu_mv_1psingle]=cm_1pass_multivariate_singleton(data_mv,order_mv);

%test the naive multivariate 2 pass formula - just for reference
order_mv=2;
[cm_mv_2pnaive, mu_mv_2pnaive]=cm_2pass_multivariate_naive(data_mv,order_mv);


%--------------------------------------------------------------------------





%TESTING UNIVARIATE T-TEST-------------------------------------------------

%dataset simulation
data1=normrnd(9,2,100,3);
data2=normrnd(9,2,100,3);

%specify the order of the t-test
tt_order_univ = 1;
%specify the significance level alpha of the t-test
alpha_univ = 0.00001;

%test the standard matlab t-test formula 
%(no clue what the official MATLAB implementation is)
[h,p,ci,t_matlab]= ttest(data1,data2,alpha_univ);

%test the higher order univariate t-test formula.
%first choose a method to compute the central moments
chosen_univ_cm_method = @cm_2passcorrected_serial;
%then compute the higher order t test
[tval_with2pcorr, tnorm_2pcorr, tstu_2pcorr, cm1, mu1, cm2, mu2]=ttest_ho_univariate(data1, data2, tt_order_univ, chosen_univ_cm_method,alpha_univ);
%same with another method 
chosen_univ_cm_method = @cm_1pass_univariate_singleton;
[tval_1p_univ, thnorm_1psingle, thstu_1psingle]=ttest_ho_univariate(data1, data2, tt_order_univ, chosen_univ_cm_method,alpha_univ);

%--------------------------------------------------------------------------


%TESTING MULTIVARIATE T-TEST ----------------------------------------------

no_traces=100;
no_values=4;

%value simulation
r1=randi(no_values,1,no_traces)-1; 
r2=randi(no_values,1,no_traces)-1; 
r3=randi(no_values,1,no_traces)-1; 
r4=randi(no_values,1,no_traces)-1; 

in_rand=randi(no_values,1,no_traces)-1;
in_fixed=zeros(1,no_traces); 

masked_value_inrand=bitxor(in_rand,r1); 
masked_value_infixed=bitxor(in_fixed,r2);

%dataset simulation
data3=[r1' masked_value_inrand']+normrnd(0,1,size(in_rand,1),2);
data4=[r2' masked_value_infixed']+normrnd(0,1,size(in_fixed,1),2);

%specify the order of the t-test
tt_order_mv = 2;
%specify the significance level alpha of the t-test
alpha_mv = 0.00001;

%test the higher order multivariate t-test formula.
%first choose a method to compute the central moments
chosen_mv_cm_method = @cm_1pass_multivariate_singleton;
%then compute the higher order t test
[tval_1p_mv, th_norm_1pmv, th_stu_1pmv, cm3, mu3, cm4, mu4]=ttest_ho_multivariate(data3, data4, tt_order_mv, chosen_mv_cm_method,alpha_mv);

cvm3=cov(data3);
cvm4=cov(data4);


%--------------------------------------------------------------------------

