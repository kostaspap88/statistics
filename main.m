%main file of High Performance Statistics 
%author: Kostas - kostaspap88@gmail.com - kpcrypto.net

%CAROLINE DELETED
%GOODBYE...CAROLINE

clear all
close all

%the main function tests the equivalence and numerical stability of
%different computational variants for central moments

format longG

%dataset simulation
data=normrnd(9,2,1,2^10);

%dataset length
n=length(data);
%specify order
order=6;

%test the standard matlab formula (no clue what the official MATLAB implementation is)
cm_matlab=cell(1,order);
for i=1:order
    cm_matlab{i}=moment(data,i);   
end
mu_matlab = mean(data);

%test the corrected 2 pass formula (serialized)
[cm_2pcorr_serial mu_2pcorr_serial]=cm_2passcorrected_serial(data,order);

%test the corrected 2 pass formula (vectorized)
[cm_2pcorr_vectorized mu_2pcorr_vectorized]=cm_2passcorrected_vectorized(data,order);

%test the naive 2 pass formla (serialized) - just for reference
[cm_2pnaive_serial mu_2pnaive_serial]=cm_2passnaive_serial(data,order);

%test the singleton 1 pass formula (serialized by construction)
[cm_1psingle mu_1psingle]=cm_1pass_singleton(data,order);

%-------------------------------------------------------------------------

%dataset simulation
data1=normrnd(9,2,50,10);
data2=normrnd(9,2,50,10);

%specify the order of the t-test
tt_order = 1;

%test the standard matlab t-test formula (no clue what the official MATLAB implementation is)
t_matlab= ttest(data1,data2,0.0001);

%test the higher order t-test formula.
%first choose a method to compute the central moments
chosen_cm_method = @cm_1pass_singleton;
t_ho=ttest_ho(data1, data2, tt_order, chosen_cm_method);



%display results    
% cm_matlab
% cm_2pcorr_serial
% cm_2pcorr_vectorized
% cm_2pnaive_serial
% cm_1psingle
%cm_1ppair
% cm_AB{1}