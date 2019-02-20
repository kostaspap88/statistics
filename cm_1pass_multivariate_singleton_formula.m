%Description: Generating the multivariate formula in unrolled fashion,
%based on formula (12) of Schneider and Appendix C.

%Input: order, powerset of multiset [1,2,...,noPOIs,1,2,...,noPOIs]

function [] = cm_1pass_multivariate_singleton_formula(order,powerset)


%loop from 1 until specified order
for current_order=2:order 
    
    %get all elements of the current order powerset
    current_powerset = powerset{current_order};
    %POIs in the current order powerset
    current_no_poi = size(current_powerset,2);
    
    
    %loop through all elements (i.e. sets) of the current order powerset
    for set_index=1:size(current_powerset,1)
             
        current_set = current_powerset(set_index,:);
        
        fprintf ('\n');
        fprintf ('C_%d,Q'',',current_order);
        fprintf ('%s',mat2str(current_set));
        fprintf ('=');
        fprintf ('C_%d,Q,',current_order);
        fprintf ('%s',mat2str(current_set));
        fprintf (' + ');
        

        %loop from 2 until number of POIs -1
        for k=2:(current_no_poi-1)
       
            %get the k-th order powerset of current set
            inner_current_powerset = nchoosek(current_set,k);
            %loop over the elements (i.e. sets) of the k-th order powerset
            for inner_set_index=1:size(inner_current_powerset,1)
                
                inner_current_set = inner_current_powerset(inner_set_index,:);
                  
                %create set JexceptS
                uA  = unique(current_set);
                hca = histc(current_set,uA); 
                hcb = histc(inner_current_set,uA);
                JexceptS = repelem(uA,hca-hcb);
                
                %print formula
                if (~isempty(JexceptS))
                    fprintf ('C_%d,',k);
                    fprintf ('%s',mat2str(inner_current_set));
                    fprintf ('*PROD(Delta');
                    fprintf ('%s',mat2str(JexceptS));
                    fprintf ('/(-n))+');
                end                
                
            end         
        end
        
        %print formula
        fprintf ('PROD(Delta');
        fprintf ('%s',mat2str(current_set));
        fprintf (')');
        fprintf (' * [(-1)^%d * (n-1) + (n-1)^%d]/[ n^%d ]',current_order,current_order,current_order);
        fprintf ('\n');

    end

end

end