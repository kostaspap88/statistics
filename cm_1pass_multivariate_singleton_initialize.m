%Initialize the multivariate 1 pass singleton variables. See Schneider formula (12)

%Input: first trace, moment order, standard set or multiset choice
%Outpt: the first central moments, first mu, first length,
%the powerset of standard set [1,2,...,noPOIs] or the powerset of multiset
%[1,2,...,noPOIs,1,2,...,noPOIs], the powerset lookuptable for quicker
%indexing by cm_1pass_multivariate_singleton_update function

function [first_cm, first_mu, first_length, ps, ps_lut] = cm_1pass_multivariate_singleton_initialize(first_singleton,order,set_choice)

%number of POIs
no_poi = length(first_singleton);

%compute the first mean
first_mu=first_singleton;

%init the LUT
ps_lut = containers.Map();

%init the powerset
ps=cell(1,order);

%generate 1st-order powerset until current-order powerset
for current_order=2:order
    
    if (set_choice==0)
        set=unique(sort(nchoosek(1:no_poi,current_order),2),'rows');
    end
    
    if (set_choice==1)
        set=unique(sort(nchoosek([1:no_poi 1:no_poi],current_order),2),'rows');
    end
    
    ps{current_order}=set;

    for index=1:size(set,1)
        if (isequal(set(index,:),1:no_poi))
            numerator_index=index;
        end
        
        if (isequal(set(index,:),sort([1:no_poi 1:no_poi])))
            denominator_index=index;
        end
        
        first_cm{current_order,index}=0;   
        ps_lut(mat2str(set(index,:)))=index;
        
    end
end


%compute the first set length
first_length=1;


end