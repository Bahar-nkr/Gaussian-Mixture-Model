function [ newsigma ] = newsigma( landa,r,moo )
pro=zeros(size(r,2),size(r,2));
for i=1:size(r,1)
    pr=(landa(i)*(r(i,:)-moo)'*(r(i,:)-moo));
    pro=pro+pr;
end
newsigma=pro./sum(landa);

end

