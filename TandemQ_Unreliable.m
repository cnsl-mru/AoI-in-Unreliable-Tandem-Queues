    %% Initialization
clear ; close all; clc
%% ==================== Part 1: Basic Function ====================
% Complete warmUpExercise.m
fprintf('Average AoI Running ... \n');

fprintf('Average AoI: \n');
fprintf('Program paused. Press enter to continue.\n');
pause;
%% ======================= Part 2: Base Line Model=======================
fprintf('Plotting Data ...for M/M/1 tandem queue\n')
rho=[];
la1=0;
mu1=1;mu2=1;
%mu1=0.7;mu2=1;
%mu1=0.9;mu2=1;
m1=1;
m2=1;
AoI1=[]
apeak=[]
mu3=0.5:0.2:1.1;
for k=0.1:0.01:0.9
     la1=k;
           for i=1:2
                 rho(i)=la1/mu1;
           end
            AoI11=(1/mu1)*(1+1/rho(1)+rho(1)^2/(1-rho(1)))+(1/mu2)+la1^2/(mu2^2*(mu2-la1))+la1^2/(mu1*mu2*(mu1+mu2-la1))
            AoI1(m1)=AoI11;
            s1=1/la1;
            for j=1:m2
               s1=s1+1/(mu3(i)-la1)
            end
            apeak(m1)=s1;
            m1=m1+1
end
x=(0.1:0.01:0.9)
plot(x,AoI1,'*--')
hold on;
%% %%====================Part-3: TandemQ with unreliable===================
fprintf('Plotting Data ...for M/M/1 tandem queue\n')
rho=[];
la1=0;
m2=4;
ga1=1;
aaoi1=[]
apeak1=[]
%mu1=1;mu2=0.7;
mu1=1;mu2=0.1;
%mu1=0.7;mu2=0.9;
mu3=0.3:0.2:1.1;
AoI1=[];
al1=0.1:0.2:0.5;
alpa=1/0.5;
m1=1;
for k=0.1:0.01:0.8
    la1=k;
  p0=ga1/(ga1+alpa)-la1*(1/mu1+1/mu2);      
  %  p0=ga1/(ga1+alpa)-la1*(1/mu1+1/mu2); 
     rho(1)=la1/mu1; rho(2)=la1/mu2;
     f1=@(z2)  mu1*(z2^2)/(mu1+mu2*z2)
     D1=@(z1,z2) z1*z2*(la1*(1-z1)+ga1)+mu1*z2*(z1-z2)+mu2*z1*(z2-1);
     C1=@(z1,z2) mu1*z2*(z2-z1)+mu2*z1*(z1-z2);
     E1=@(z1) alpa/(la1*(1-z1)+ga1);
     w1=@(s1) E1(f1(1-s1/la1))*(C1(f1(1-s1/la1),(1-s1/la1))*(p0))/(E1(f1(1-s1/la1))*(D1(f1(1-s1/la1),1-s1/la1))-la1*ga1*f1(1-s1/la1)*(1-s1/la1));
     h1=@(s1) (mu1/(s1+mu1*(ga1+alpa)/ga1))*(mu2/(s1+mu2*(ga1+alpa)/ga1))*((ga1+alpa)/ga1)^2;
     %h1=@(s1) (mu1/(s1+mu1))*(mu2/(s1+mu2));
     syms apeak(s1)
        apeak(s1)=w1(s1)*h1(s1)-w1(s1+la1)*s1/(s1+la1);
        Apeak=-diff(apeak(s1),s1);
    
    syms aoi1(s1)
        aoi1(s1)=w1(s1)-(p0)*s1*h1(s1)/(s1+la1*h1(s1+la1)); 
        %aoi1(s1)=la1*(w1(s1)-apeak(s1))/s1;
         AAoi=-diff(aoi1(s1),s1);


    aaoi1(m1)=double(limit(AAoi,s1,0));
    apeak1(m1)=double(limit(Apeak,s1,0));
    m1=m1+1;
end
x=0.1:0.01:0.8
plot(x,aaoi1)
hold on
%% ======================= Part 4: Bufferless Breakdown Model with respect to lambda1 =======================
fprintf('Plotting Data ...\n')
for M=1:2% number of servers
lambda=[];
rho=[];
AoI=[];
AoI1=[];
alpa=0.9;
% % %Exponential
% gamma1=1;%expected repair time
% gamma2=2*gamma1^2;% second order repair time
% beta11=1; %expected service time
% beta12=2*beta11^2 % second order service time
% ES=beta11*(1+alpa*gamma1)% expected service time including breakdowns

%Erlang-2
gamma1=1; %expected repair time
gamma2=2*gamma1^2; % second order moment of repair time
beta11=1; %expected service time
beta12=2*beta11^2 % second order moment of service time
ES=beta11*(1+alpa*gamma1)% expected service time including breakdowns


% %HyperExponential with p
% p1=0.9;
% gamma1=p1*1/2+(1-p1)*1/(2^2); % expected repair time
% gamma2=p1*2/(2^2)+(1-p1)*2/(2^4); % second order moment of repair time
% beta11=p1*0.8+(1-p1)*0.8^2;   % expected service time
% beta12=p1*2*0.8^2+(1-p1)*2*0.8^4; % second order moment of service time
% ES=beta11*(1+alpa*gamma1) % expected service time including breakdowns
m1=1;
%AAoI manipulation for different sources
for k=0.1:0.01:0.4
    lambda(1)=k;
    for i=2:M
    lambda(i)=0.1;
    end
     l1=sum(lambda);
    for i=1:M
    rho(i)=lambda(i)*ES;
    end
    rho(1)
    r1=sum(rho);
ES2=l1*(beta11*alpa*gamma2+r1*beta12)  %2nd order moment of service time
EW=l1*ES2/(2*(1-r1)); % Expected waiting Time
LS1=(1/ES)/((1/ES)+lambda(1));    % S*(lambda(1)))
LS2=-(1/ES)/((1/ES)+lambda(1))^2; % S*'(lambda(1))) 
LS3=(2/ES)/((1/ES)+lambda(1))^3;  % S*''(lambda(1))) 
WS1=((1-r1)*lambda(1)*LS1)/((lambda(1)-l1*(1-LS1)))
WS2=(1-r1)*(l1*LS1^2+(lambda(1)^2-lambda(1)*l1)*LS2-l1*LS1)/(lambda(1)-l1*(1-LS1))^2
WS3=0;
for k1=2:M
LS11=(1/ES)/((1/ES)+lambda(k1));    % S*(lambda(1)))
LS21=-(1/ES)/((1/ES)+lambda(k1))^2; % S*'(lambda(1))) 
LS31=(2/ES)/((1/ES)+lambda(k1))^3;  % S*''(lambda(1))) 
WS11=((1-r1)*lambda(1)*LS11)/((lambda(1)-l1*(1-LS11)))
WS21=(1-r1)*(l1*LS11^2+(lambda(1)^2-lambda(1)*l1)*LS21-l1*LS11)/(lambda(1)-l1*(1-LS11))^2;
WS3=WS3+lambda(k1)*((2/lambda(1))+WS21-(2.0*WS11/lambda(1))) 
end;
AoI(m1)=EW+2.0*ES+(2.0*WS1/lambda(1))-WS2-(1/lambda(1))+ES*WS3
AoI1(m1)=EW+2.0*ES+1/lambda(1)
m1=m1+1;
end;
x=(0.1:0.01:0.4)
p(M)=plot(x,AoI,'--')
hold on;
end;
%legend([p(2),p(3),p(4),p(5)],'k=2','k=3','k=4','k=5')

%% ======================= Part 5: Bufferless Breakdown Model with respect to lambda1 =======================
fprintf('Plotting Data ...\n')
for M=1:1% number of servers
lambda=[];
rho=[];
AoI=[];
AoI1=[];
beta11=[]
SS=[];
SS1=[];
alpa=1;l1=0.4
% %Exponential
gamma1=0.5;%expected repair time
gamma2=2*gamma1^2;% second order repair time
beta11(1)=1; %expected service time
beta12(1)=2*beta11(1)^2 % second order service time
SS(1)=beta11(1)*(1+alpa*gamma1)% expected service time including breakdowns
ES=0;ES2=0;
for i=2:M
    beta11(i)=i*0.5;
    beta12(i)=2*beta11(1)^2; % second order service time
    SS(i)=beta11(i)*(1+alpa*gamma1);% expected service time including breakdowns
    SS1(i)=(beta11(i)*alpa*gamma2+l1*SS(i)*beta12(i));
    ES=ES+SS(i);
    ES2=ES2+SS1(i)+2*SS(i)*SS(i-1);
end
% %Erlang-2
% gamma1=1; %expected repair time
% gamma2=1*gamma1^2; % second order moment of repair time
% beta11=0.89; %expected service time
% beta12=1*beta11^2 % second order moment of service time
% ES=beta11*(1+alpa*gamma1)% expected service time including breakdowns


% %HyperExponential with p
% p1=0.9;
% gamma1=p1*1/2+(1-p1)*1/(2^2); % expected repair time
% gamma2=p1*2/(2^2)+(1-p1)*2/(2^4); % second order moment of repair time
% beta11=p1*0.8+(1-p1)*0.8^2;   % expected service time
% beta12=p1*2*0.8^2+(1-p1)*2*0.8^4; % second order moment of service time
% ES=beta11*(1+alpa*gamma1) % expected service time including breakdowns
m1=1;
%AAoI manipulation for different sources
for k=0.1:0.01:0.55
    lambda(1)=k;
    r1=k*1/beta11(M);
    st=[];
    rt1=[];
    s2=1;
    for j=1:M
        rt1(j)=gamma1/(gamma1+k);
        st(j)=(k+alpa*(1-rt1(j)))/(k+alpa*(1-rt1(j))+1/ES);
        s2=s2*st(j)
    end
AoI(m1)=ES+k*ES2/(2*(1-r1))+(1-r1)/(k*s2)
%AoI1(m1)=EW+2.0*ES+1/lambda(1)
m1=m1+1;
end;
x=(0.1:0.01:0.55)
p(M)=plot(x,AoI,'--')
hold on;
end;
%legend([p(2),p(3),p(4),p(5)],'k=2','k=3','k=4','k=5')

%% %% ======================= Part 6: Bufferless Breakdown Model AAoI with respect to CV^2=======================
fprintf('Plotting Data ...\n')
for M=1:4% number of servers
lambda=[];
rho=[];
AoI=[];
AoI1=[];
alpa=0.6;
% % %Exponential
% gamma1=1;%expected repair time
% gamma2=2*gamma1^2;% second order repair time
% beta11=1; %expected service time
% beta12=2*beta11^2 % second order service time
% ES=beta11*(1+alpa*gamma1)% expected service time including breakdowns

%Erlang-m2
% m2=2;
% gamma1=1; %expected repair time
% gamma2=(1/m2)*gamma1^2; % second order moment of repair time
% beta11=1; %expected service time
% beta12=(1/m2)*beta11^2 % second order moment of service time
% cv2=1/m2
% ES=beta11*(1+alpa*gamma1)% expected service time including breakdowns


% %HyperExponential with p
cv2=1
m2=1
p1=(1+sqrt((cv2-1)/(cv2+1)))/2;
j1=2*p1/m2;j2=2*(1-p1)/m2;
gamma1=p1*1/j1+(1-p1)*1/j2; % expected repair time
gamma2=p1*2/(j1^2)+(1-p1)*2/(j2^2); % second order moment of repair time
beta11=p1*1/j1+(1-p1)*1/j2;   % expected service time
beta12=p1*2/j1^2+(1-p1)*2/j2^2; % second order moment of service time
ES=beta11*(1+alpa*gamma1) % expected service time including breakdowns
m1=1;
%AAoI manipulation for different sources
for k=0.1:0.01:0.35
    lambda(1)=k;
    for i=2:M
    lambda(i)=0.1;
    end
     l1=sum(lambda);
    for i=1:M
    rho(i)=lambda(i)*ES;
    end
    rho(1)
    r1=sum(rho);
ES2=l1*(beta11*alpa*gamma2+r1*beta12)  %2nd order moment of service time
EW=l1*ES2/(2*(1-r1)); % Expected waiting Time
LS1=(1/ES)/((1/ES)+lambda(1));    % S*(lambda(1)))
LS2=-(1/ES)/((1/ES)+lambda(1))^2; % S*'(lambda(1))) 
LS3=(2/ES)/((1/ES)+lambda(1))^3;  % S*''(lambda(1))) 
WS1=((1-r1)*lambda(1)*LS1)/((lambda(1)-l1*(1-LS1)))
WS2=(1-r1)*(l1*LS1^2+(lambda(1)^2-lambda(1)*l1)*LS2-l1*LS1)/(lambda(1)-l1*(1-LS1))^2
WS3=0;
for k1=2:M
LS11=(1/ES)/((1/ES)+lambda(k1));    % S*(lambda(1)))
LS21=-(1/ES)/((1/ES)+lambda(k1))^2; % S*'(lambda(1))) 
LS31=(2/ES)/((1/ES)+lambda(k1))^3;  % S*''(lambda(1))) 
WS11=((1-r1)*lambda(1)*LS11)/((lambda(1)-l1*(1-LS11)))
WS21=(1-r1)*(l1*LS11^2+(lambda(1)^2-lambda(1)*l1)*LS21-l1*LS11)/(lambda(1)-l1*(1-LS11))^2;
WS3=WS3+lambda(k1)*((2/lambda(1))+WS21-(2.0*WS11/lambda(1))) 
end;
AoI(m1)=EW+2.0*ES+(2.0*WS1/lambda(1))-WS2-(1/lambda(1))+ES*WS3
EW=l1*ES2/(2*(1-r1));
AoI1(m1)=EW+2.0*ES+1/lambda(1)
m1=m1+1;
end;
x=(0.1:0.01:0.35)
p(M)=plot(x,AoI,'--g')
hold on;
end;
%legend([p(2),p(3),p(4),p(5)],'k=2','k=3','k=4','k=5')

