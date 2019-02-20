%Data wrapper around the multivariate singleton 1 pass formula

function [cm, mu] = cm_1pass_multivariate_singleton(data,order)

%number of POIs in dataset
no_poi=size(data,2);

%if order is equal to the number of POIs then we use the standard set
%[1,2,...,noPOIs]. If it is equal to double the number of POIs we use the
%multiset [1,2,...,noPOIs,1,2,...,noPOIs]. If neither holds then the
%function stops.

switch no_poi
    case order
        set_choice=0;
    case order/2
        set_choice=1;
    otherwise
        fprintf ('POIs do not correspond to order!\n');
        cm=[];
        mu=[];
        return
end

%first, initialize using the first data element data(1)
[cm_1psingle,mu_1psingle,len_1psingle, powerset,ps_lut] = cm_1pass_multivariate_singleton_initialize(data(1,:),order,set_choice);

%then create the formula and print it
cm_1pass_multivariate_singleton_formula(order,powerset);

%then iterate on dataset
%note that the data length does not need to be known in advance 
%note also that we can keep updating with more data, unlike the 2pass
%formulas where any new incoming data requires us to pass again

n=size(data,1);

for i=2:n
    
    %perform incremental update on cm,mu,length using every singleton
    [cm_1psingle,mu_1psingle,len_1psingle] = cm_1pass_multivariate_singleton_update(data(i,:),order,...
                                                        cm_1psingle,mu_1psingle,len_1psingle,powerset,ps_lut);
                                                    
end

%finally, compute moment instead of sum
for i=2:order
    for j=1:size(powerset{i},1)
        cm_1psingle{i,j}=cm_1psingle{i,j}/n;
    end
end

%return moments and mean
cm=cm_1psingle;
mu=mu_1psingle;


end
