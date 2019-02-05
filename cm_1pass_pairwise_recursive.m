%Description: Central moment computation using corrected 1-pass formula from "Numerically stable, scalable
%formulas for parallel and online computation of ho multivariate central moments with arbitary weights,
%Pebay, Terriberry, Kolla, Bennet". See formula (3.1).

%The formula combines the central moments of setA and setB
%Note that the formula does not need to know the data size a priori.

%Input: mean, central moments and length of set A and set B
%Output current mean, central moments and length of A union B n


function [upd_cm,upd_mu,len] = cm_1pass_pairwise_recursive(order, dataA, dataB, cm_setA, cm_setB, mu_setA, mu_setB, lenA, lenB)

if ((lenA==1)&&(lenB==1))
    mu_setA = dataA;
    mu_setB = dataB;
    for i=1:order
        cm_setA{i} = moments(dataA,i);
        cm_setB{i} = moments(dataB,i);
    end
else
    
    cm_1pass_pairwise_recursive(order, dataA(1:lenA/2), dataB, cm_setA, cm_setB, mu_setA, mu_setB, lenA, lenB)

%init updated central moments
upd_cm=cell(1,order);
upd_cm(1,:)={0};

%convert to sum 
for i=1:order
    cm_setA{i}=cm_setA{i}*lenA;
    cm_setB{i}=cm_setB{i}*lenB;
end

%update the length
len=lenA+lenB;

%compute delta
delta = mu_setA - mu_setB;

%update mean
upd_mu = mu_setA + lenB/len * delta;

for i=1:order
    
    term1 = lenA*(-lenB/len * delta)^i + lenB*(lenA/len * delta)^i;
    term2 = 0;
    for k=1:i-2
        term2 = term2 + nchoosek(i,k)*delta^k * ( cm_setA{i-k}*(-lenB/len)^k + (cm_setB{i-k}*(lenA/len)^k) ) ;
    end
    upd_cm{i} = cm_setA{i} + cm_setB{i} + term1 + term2;
    
end

%convert to central moment 
for i=1:order
    upd_cm{i}=upd_cm{i}/len;
end


end