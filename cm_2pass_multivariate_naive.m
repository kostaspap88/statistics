

function [cm, mu] = cm_2pass_multivariate_naive(data,order)

no_poi = size(data,2);

%if order is equal to the number of POIs then we use the standard set
%[1,2,...,noPOIs]. If it is equal to double the number of POIs we use the
%multiset [1,2,...,noPOIs,1,2,...,noPOIs]. If neither holds then the
%function stops.



switch no_poi
    case 1
        fprintf ('Try the univariate formula!\n');
        cm=[];
        mu=[];
        return   
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





mu=mean(data);


%naive loop
for i=1:size(data,1)        
    for current_order=2:order
        
    if(set_choice==0)
        powerset=unique(sort(nchoosek(1:no_poi,current_order),2),'rows'); 
    end
    if(set_choice==1)
        powerset=unique(sort(nchoosek([1:no_poi 1:no_poi],current_order),2),'rows'); 
    end
    
        for set_index=1:size(powerset,1)
            set=powerset(set_index,:);
            if (i==1)
                cm{current_order,set_index}=prod(data(i,set)-mu(set));
            else
                cm{current_order,set_index}=cm{current_order,set_index} + prod(data(i,set)-mu(set));
            end
        end

    end
end

n=size(data,1);

for i=2:order
    for j=1:size(cm(i,:),2)
        cm{i,j}=cm{i,j}/n;
    end
end



end