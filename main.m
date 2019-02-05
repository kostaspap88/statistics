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

% n=legnth(data);
% current_length=n;
% 
% 
% current_length = 1;
% index_counter = 1;
% for i=1:2:n
%     %pairwise 
%     
%     
% end
% 
% 
% current_length = 1;
% index_counter = 1;
% for i=1:index_counter
%     %pairwise
%     
%         %compute a pairwise update
%     [cm_AB{index}, mu_AB{index}, lenAB{index}] = cm_1pass_pairwise(order,cm_setA, cm_setB, mu_setA, mu_setB, lenA, lenB);
%     
%     
% end


%  - we need data wrapper around this
% %first compute the central moments of set A and set B (using e.g. corrected 2 pass)
% %set A (one third of data set)
% indexA = 1:floor(length(data)/3);
% lenA = length(indexA);
% [cm_setA mu_setA] = cm_2passcorrected_serial(data(indexA),order);
% %set B (two thirds of data set)
% indexB = floor(length(data)/3)+1 : length(data);
% lenB = length(indexB);
% [cm_setB mu_setB]= cm_2passcorrected_serial(data(indexB),order);
% %combine the central moments of sets A and B
% [cm_1ppair,mu_1ppair,len_1ppair] = cm_1pass_pairwise(order,cm_setA, cm_setB, mu_setA, mu_setB, lenA, lenB);


%test the pairwise 1 pass formula - we need data wrapper around this
%start form singletons and combine them in a pairwise manner

% index1 = 1:2:n;
% index2 = 2:2:n;
% 
% index = 1;
% for i=1:2:n
%     
%     %initialize length, means and central moments pairwise
%     lenA=1;
%     lenB=1;
%     mu_setA = data(i);
%     mu_setB = data(i+1);
%     for p=1:order
%         cm_setA{p} = moment(data(i),i);
%         cm_setB{p} = moment(data(i+1),i);
%     end
%     %compute a pairwise update
%     [cm_AB{index}, mu_AB{index}, lenAB{index}] = cm_1pass_pairwise(order,cm_setA, cm_setB, mu_setA, mu_setB, lenA, lenB);
%     index=index+1;
% end
% 
% while (index>1)
% for i=1:2:(index-1)
%     
%      cmA= cm_AB{i};
%      cmB=cm_AB{i+1};
%      mA= mu_AB{i};
%      mB= mu_AB{i+1};
%      lA=lenAB{i};
%      lB=lenAB{i+1};
%      [cm_AB{i}, mu_AB{i}, lenAB{i}] =cm_1pass_pairwise(order,cmA,cmB,mA,mB,lA,lB);
%     
% end
% index = length(1:2:(index-1));
% end


%display results    
cm_matlab
cm_2pcorr_serial
cm_2pcorr_vectorized
cm_2pnaive_serial
cm_1psingle
%cm_1ppair
% cm_AB{1}