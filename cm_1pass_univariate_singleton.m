%Data wrapper around the univariate singleton 1 pass formula

function [cm, mu] = cm_1pass_univariate_singleton(data,order)

%first, initialize using the first data element data(1)
[cm_1psingle,mu_1psingle,len_1psingle] = cm_1pass_univariate_singleton_initialize(data(1),order);

n=length(data);

%then iterate on dataset
%note that the data length does not need to be known in advance 
%note also that we can keep updating with more data, unlike the 2pass
%formulas where any new incoming data requires us to pass again
for i=2:n
    %perform incremental update on cm,mu,len using every singleton
    [cm_1psingle,mu_1psingle,len_1psingle] = cm_1pass_univariate_singleton_update(data(i),order,cm_1psingle,mu_1psingle,len_1psingle);
end

%finally, compute moment instead of sum
for i=1:order
    cm_1psingle{i}=cm_1psingle{i}/n;
end

%return moments and mean
cm=cm_1psingle;
mu=mu_1psingle;


end
