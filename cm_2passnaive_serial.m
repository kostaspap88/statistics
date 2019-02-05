%A naive serialized implementation of central moments. The formula is the
%standard statistical definition for moment: SUM(x-mu)^p 
%Implemented just for reference

%Input: data and specified order
%Output: central moments and mean

function [cm,mu] = cm_2passnaive_serial(data,order)

%init central moments
cm=cell(1,order);
cm(1,:)={0};

%compute the length of the data
n=length(data);

%1st pass computes the mean
mu=sum(data)/n;


%2nd pass computes the sums required
for j=1:order
    for i=1:n
        delta=data(i)-mu;
        cm{j}=cm{j}+delta^j;
    end
    cm{j}=cm{j}/n;
end


end