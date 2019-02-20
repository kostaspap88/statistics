%Initialize the univariate 1 pass singleton variables. See Pebay formula (3.7)

%Input: first trace, order
%Outpt: the first central moments, first mu, first length

function [first_cm, first_mu, first_length] = cm_1pass_univariate_singleton_initialize(first_singleton,order)

%init central moments
first_cm=cell(1,order);
first_cm(1,:)={0};

%compute the first mean
first_mu=first_singleton;
%compute the first central moments
for i=1:order
    first_cm{i}=moment(first_singleton,i);      
end

%compute the first set length
first_length=1;


end