%Description: Central moment computation using corrected 2-pass formula from "Numerically stable, scalable
%formulas for parallel and online computation of ho multivariate central moments with arbitary weights,
%Pebay, Terriberry, Kolla, Bennet". See formula (3.40).


%This standard implementation does not vectorize operations s.t. the memory
%requirements stem only from data. I.e. the extra variables needed (delta,
%thetas) have minimal impact.


%Inputs: the data vector and the statistical order of the central moment 
%Output: the central moments from 1 until the specified order and the mean

function [cm,mu]=cm_2passcorrected_serial(data,order)

%init central moments
cm=cell(1,order);
cm(1,:)={0};

%compute the length of the data
n=length(data);

%1st pass computes the mean
mu=mean(data);


%2nd pass computes the sums required

%init thetas (exponents of delta) and their sums
theta=cell(1,order+1);
sum_theta=cell(1,order+1);
sum_theta(1,:)={0};

%2nd pass
for i=1:n
    
    %delta computation
    delta=data(i)-mu; 
    
    for k=0:order
        %theta computation (exponents of delta)
        theta{k+1}=delta^k;
        %sum of thetas
        sum_theta{k+1}=sum_theta{k+1}+theta{k+1};
    end       
end


%Compute all central moments from 1 until specified order
for i=1:order
    current_order=i;
    %combine sums to compute cm{i}
    for k=0:current_order
        t=current_order-k;
        %note that sum_theta{2} is equal to delta
        cm{i}=cm{i}+nchoosek(current_order,k)*sum_theta{t+1}*((-1/n)*sum_theta{2})^k;   
    end
    cm{i}=(1/n)*cm{i};
end


end