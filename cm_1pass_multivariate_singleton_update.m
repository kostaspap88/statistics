%Description: Central moment computation using 1-pass formula from "Leakage Assessment Methodology
%a clear roadmap for side-channel evaluations, Schneider, Moradi". See formula (12).

%The formula is sigleton incremental so it has minimal memory requirements.
%Note that the formula does not need to know the data size a priori.

%Input: the singleton value, order, the current central moments, the
%current data length and the current mean
%Output: the updated central moments, updated length, updated mean

function [upd_cm,upd_mu,upd_length] = cm_1pass_multivariate_singleton_update(singleton,order,current_cm,current_mu,...
                                                      current_length,powerset,ps_lut)

%init updated central moments
size_upd_cm=size(current_cm);
upd_cm=cell(size_upd_cm);

%update length, i.e. increase length by 1 since we use singleton
upd_length=current_length+1;
cl=upd_length;

%compute delta
delta=singleton-current_mu;
%update the mean
upd_mu = current_mu + delta./cl; 


%loop from 1 until specified order
for current_order=2:order 
    
    %get all elements of the current order powerset
    current_powerset = powerset{current_order};
    current_no_poi = size(current_powerset,2);
    %loop through all elements (i.e. sets) of the current order powerset
    for set_index=1:size(current_powerset,1)
             
        current_set = current_powerset(set_index,:);
        %set_entry = table_powerset(mat2str(current_set));
        %term1 = current_cm{order,set_entry};
        term1 = current_cm{current_order,set_index};
        

        %loop from 2 until number of POIs -1
        term2 = 0;
        for k=2:(current_no_poi-1)
       
            %get the k-th order powerset of the current set
            inner_current_powerset = nchoosek(current_set,k);
            %loop over the elements (i.e. sets) of the k-th order powerset
            for inner_set_index=1:size(inner_current_powerset,1)
                
                inner_current_set = inner_current_powerset(inner_set_index,:);
                %inner_set_entry = table_powerset(mat2str(inner_current_set));
                         
                %create set JexceptS 
                uA  = unique(current_set);
                hca = histc(current_set,uA); 
                hcb = histc(inner_current_set,uA);
                JexceptS = repelem(uA,hca-hcb);

                
                if (~isempty(JexceptS))
                    prod_term = prod(delta(JexceptS)./(-cl));
                    %term2 = term2 + current_cm{inner_set_entry} * prod_term;
                    term2 = term2 + current_cm{k,ps_lut(mat2str(inner_current_set))} * prod_term;
                end
            end         
        end
        
        term3 = prod(delta(current_set)) * ( (-1)^current_order * (cl-1) + (cl-1)^current_order )/(cl^current_order) ;  
        %upd_cm{set_entry} = term1 + term2 + term3;
        upd_cm{current_order,set_index} = term1 + term2 + term3;
    end

end

end