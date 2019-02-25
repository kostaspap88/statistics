
%Description: Generating the multivariate formula in unrolled fashion,
%based on formula (12) of Schneider and Appendix C. Used for python
%functions

%Input: order, powerset of multiset [1,2,...,noPOIs,1,2,...,noPOIs]

function [] = cm_1pass_multivariate_singleton_formula_python(order,powerset,flat_lut,no_poi)

for i=1:no_poi
    fprintf('delta%d = delta[%d]',i-1,i-1);
    fprintf('\n');
end

for current_order=2:order 
    current_powerset = powerset{current_order};
    current_no_poi = size(current_powerset,2);
      for set_index=1:size(current_powerset,1)
         
        current_set = current_powerset(set_index,:);
        t4 = current_set-1;
        t4=mat2str(t4);
        t4(t4=='[') = [];
        t4(t4==']') = [];
        t4= t4(find(~isspace(t4)));
        fprintf('delta%s = np.prod([ ',t4);
        for j=1:(length(t4)-1)
            fprintf('delta[%d], ',current_set(j)-1);
        end
        
        fprintf('delta[%d] ])\n',current_set(length(t4))-1);
      end
end 



full_index = 0;

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
        fprintf ('upd_cm[%d] = current_cm[%d]',full_index,full_index);
        fprintf (' + ');
        full_index = full_index+1;

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
                    
                    
                    t3=flat_lut(mat2str(inner_current_set))-1;
                    fprintf('current_cm[%d]',t3);
                    fprintf (' * delta');
                    
                    t2=JexceptS-1;
                    t2=mat2str(t2);
                    t2(t2=='[') = [];
                    t2(t2==']') = [];
                    t2= t2(find(~isspace(t2)));
                    
                    fprintf ('%s',t2);
                    
                    fprintf ('/((-cl)**%d) + ', length(t2));
                end                
                
            end         
        end
        
        %print formula
        fprintf ('delta');
        t1=current_set-1;
        t1=mat2str(t1);
        t1(t1=='[') = [];
        t1(t1==']') = [];
        t1= t1(find(~isspace(t1)));
        fprintf ('%s',t1);

        fprintf (' * final_term[%d]',current_order-2);
        fprintf ('\n');

    end

end

end