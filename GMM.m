clc;clear all;close all
T=100e-6;
w=5e6;
P=200;
a=-4;
samples=T.*w;
etha=10.^(-17.4);
Cn_x=[0 1].*1000;
Cn_y=[0 0].*1000;
plot(Cn_x,Cn_y,'o');axis([-500 1500 -1000 500]);
Cm_x=[0.5 0].*1000;
Cm_y=[0 -0.5].*1000;
hold on;plot(Cm_x,Cm_y,'^'),grid on;legend('SU','PU')
m=2;n=2;
for i=1:m
    for j=1:n
        g(i,j)=sqrt(((Cm_x(i)-Cn_x(j)).^2)+((Cm_y(i)-Cn_y(j)).^2));
    end
end
d=g.^a;
prod=(2.*T.*P.*10.^-3)./etha;
mean0=(ones(2,1).*2.*samples)';
mean1=(2.*samples)+(d(1,:).*prod);
mean2=(2.*samples)+(d(2,:).*prod);
mean3=(2.*samples)+((d(1,:)+d(2,:)).*prod);
sigma0=diag(ones(2,1).*4.*samples);
sigma1=diag((4.*samples)+(4.*(d(1,:).*prod)));
sigma2=diag((4.*samples)+(4.*(d(2,:).*prod)));
sigma3=diag((4.*samples)+(4.*((d(1,:)+d(2,:)).*prod)));
r0=mvnrnd(mean0,sigma0,120);
r1=mvnrnd(mean1,sigma1,120);
r2=mvnrnd(mean2,sigma2,120);
r3=mvnrnd(mean3,sigma3,120);
r=[r0;r1; r2 ;r3];
label=[ones(120,1);ones(360,1).*2];
c1=r(label==1,:);
c2=r(label==2,:);figure;
plot(c1(:,1),c1(:,2),'c*'),hold on
plot(c2(:,1),c2(:,2),'ms');title('original data');axis([800 1400 800 1400])
moo_g1=mean0;sigma_g1=sigma0;v1=0.25;v2=1-v1;
moo_g2=[1250,1250];sigma_g2=[3000,0;0,3000];
for t=1:size(r,1)
    phi1(t)=(1./(((2.*pi).^(n/2)).*sqrt(det(sigma_g1)))).*exp(-0.5*(r(t,:)-moo_g1)*((sigma_g1)^(-1))...
        *(r(t,:)-moo_g1)');
end
for iteration=1:100
    for t=1:size(r,1)
        phi2(t)=(1./(((2.*pi).^(n/2)).*sqrt(det(sigma_g2)))).*exp(-0.5*(r(t,:)-moo_g2)*((sigma_g2)^(-1))...
            *(r(t,:)-moo_g2)');
    end
    %%%%%%%%E Step%%%%%%%%%
    landa1=(v1.*phi1)./((v1.*phi1)+(v2.*phi2));
    landa2=(v2.*phi2)./((v1.*phi1)+(v2.*phi2));
    %%%%%%%M Step%%%%%%%%%%
    moo_g2=newmoo(r,landa2);
    sigma_g2=newsigma(landa2,r,moo_g2);
    v1=newv(r,landa1);
    v2=newv(r,landa2);
end
ggg = gmdistribution([moo_g1;moo_g2], cat(3,sigma_g1,sigma_g2), [v2;v1]);
figure;
for x=1:size(r,1)
    pos(x,:)=posterior(ggg,r(x,:));
end
for y=1:size(pos,1)
    if pos(y,1)>pos(y,2)
        ll(y)=1;
    else
        ll(y)=2;
    end
end
c1=r(ll==1,:);
c2=r(ll==2,:);
plot(c1(:,1),c1(:,2),'c*'),hold on
plot(c2(:,1),c2(:,2),'ms');axis([800 1400 800 1400]);title('clustered data');
conmat=confusionmat(label,ll)
error=((conmat(1,2)+conmat(2,1))./sum(sum(conmat))).*100 %%%%IN PERCENT