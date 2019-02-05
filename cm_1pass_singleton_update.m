%Description: Central moment computation using corrected 1-pass formula from "Numerically stable, scalable
%formulas for parallel and online computation of ho multivariate central moments with arbitary weights,
%Pebay, Terriberry, Kolla, Bennet". See formula (3.7).

%The formula is sigleton incremental so it has minimal memory requirements.
%Note that the formula does not need to know the data size a priori.

%Input: the singleton value, order, the current central moments, the
%current data length and the current mean
%Output: the updated central moments, updated length, updated mean

function [upd_cm,upd_mu,upd_length] = cm_1pass_singleton_update(singleton,order,current_cm,current_mu,current_length)


%init updated central moments
upd_cm=cell(1,order);
upd_cm(1,:)={0};

%update length, i.e. increase length by 1 since we use singleton
upd_length=current_length+1;
cl=upd_length;

%compute delta
delta=singleton-current_mu;
%updated the mean
upd_mu = current_mu +  delta/cl;

%compute all central moments from order 1 until specified order
for i=1:order

    term1= ( (cl-1)/((-cl)^i) + ((cl-1)/cl)^i ) * delta^i;

    term2=0;
    for k=1:(i-2)
        term2=term2+nchoosek(i,k)*current_cm{i-k}*(-delta/cl)^k;
    end
    
    upd_cm{i}=current_cm{i}+term2+term1;
    
end



end