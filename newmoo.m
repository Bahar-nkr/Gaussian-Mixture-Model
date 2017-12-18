function [ newmoo ] = newmoo(r,landa)
for i=1:size(r,1)
    prod(i,:)=landa(i).*r(i,:);
end
newmoo=sum(prod)./sum(landa);
end

