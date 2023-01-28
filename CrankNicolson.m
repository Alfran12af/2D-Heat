clear
clc

%Data
P1=[0.5,0.4];
P2=[0.5;0.7];
P3=[1.1;0.8];

rho=[1500 1600 1900 2500];

cp=[750 770 810 930];

lambda=[170 140 200 140];

time_span=5000;    %time laps studied in seconds

T0=8+273;
Tbot=23+273;  %in kelvin
Text=33+273;
alfa_ext=9;   %convection constant
qw=54.55;     %Thermic flow in top



%Simulation parameters
deltat=1;
deltax=0.01;
deltay=0.01;
L=1.1;
W=0.8;
n=L/deltax;
m=W/deltay;
N=n+1;
M=m+1;



tol=0.1;      %Tolerance


%%
%Mesh properties
lambdaP=zeros(N,M);
cpP=zeros(N,M);
rhoP=zeros(N,M);

for i=1:N
    for j=1:M
        x=i*deltax;
        y=j*deltax;
        if x<P1(1) && y<P1(2)  %M1
            lambdaP(i,j)=lambda(1);
            cpP(i,j)=cp(1);
            rhoP(i,j)=rho(1);
        elseif x>P1(1) && y<P2(2) %M2
            lambdaP(i,j)=lambda(2);
            cpP(i,j)=cp(2);
            rhoP(i,j)=rho(2);
        elseif x<P1(1) && y>P1(2) %M3
            lambdaP(i,j)=lambda(3);
            cpP(i,j)=cp(3);
            rhoP(i,j)=rho(3);
        elseif x>P1(1) && y>P2(2) %M4
            lambdaP(i,j)=lambda(4);
            cpP(i,j)=cp(4);
            rhoP(i,j)=rho(4);
        elseif x==P1(1) && y<P1(2) %border M1-M2
            lambdaP(i,j)=0.5*(lambda(1)+lambda(2));
            cpP(i,j)=0.5*(cp(1)+cp(2));
            rhoP(i,j)=0.5*(rho(1)+rho(2));
        elseif x==P1(1) && y>P1(2) && y<P2(2) %Border M3-M2
            lambdaP(i,j)=0.5*(lambda(3)+lambda(2));
            cpP(i,j)=0.5*(cp(3)+cp(2));
            rhoP(i,j)=0.5*(rho(3)+rho(2));
        elseif x==P1(1) && y>P2(2) %Border M3-M4
            lambdaP(i,j)=0.5*(lambda(3)+lambda(4));
            cpP(i,j)=0.5*(cp(3)+cp(4));
            rhoP(i,j)=0.5*(rho(3)+rho(4));
        elseif x<P1(1) && y==P1(2) %Border M1-M3
            lambdaP(i,j)=0.5*(lambda(1)+lambda(3));
            cpP(i,j)=0.5*(cp(1)+cp(3));
            rhoP(i,j)=0.5*(rho(1)+rho(3));
        elseif x>P1(1) && y==P2(2) %Border M2-M4
            lambdaP(i,j)=0.5*(lambda(2)+lambda(4));
            cpP(i,j)=0.5*(cp(4)+cp(2));
            rhoP(i,j)=0.5*(rho(4)+rho(2));
        elseif x==P1(1) && y==P1(2) %Point P1 considering it to be in M1-M2
            lambdaP(i,j)=0.5*(lambda(1)+lambda(2));
            cpP(i,j)=0.5*(cp(1)+cp(2));
            rhoP(i,j)=0.5*(rho(1)+rho(2));
        elseif x==P2(1) && y==P2(2) %Point P2 considering it to be in M3-M4
            lambdaP(i,j)=0.5*(lambda(3)+lambda(4));
            cpP(i,j)=0.5*(cp(3)+cp(4));
            rhoP(i,j)=0.5*(rho(3)+rho(4));
        end
    end
end


%Vector definition
T=zeros(N,M,time_span/deltat+1);
T(:,:,1)=T0;  
ae=zeros(N,M);
aw=zeros(N,M);
an=zeros(N,M);
as=zeros(N,M);
ap=zeros(N,M);
b=zeros(N,M);


%%
%Coeficient determination (non-time dependent)
%[ae,aw,an,as,ap] = coeficients(ae,aw,an,as,ap);
%Left wall is convection so it depends on wall Temperature which is a
%function of time so it will be calculated later
    
for i=1:N
    for j=1:M
        if i==1 && (j~=1 && j~=M) %Left wall
            ae(i,j)=lambdaP(i+1,j);
            aw(i,j)=0;
            an(i,j)=0;
            as(i,j)=0;
            ap(i,j)=0;%Convection parametre, function of time
        elseif i==N && (j~=0 && j~=M) %Right wall
            ae(i,j)=0;
            aw(i,j)=0;
            an(i,j)=0;
            as(i,j)=0;
            ap(i,j)=1;
        elseif j==1 && (i~=1 && i~=N) %Bottom wall
            ae(i,j)=0;
            aw(i,j)=0;
            an(i,j)=0;
            as(i,j)=0;
            ap(i,j)=1;
        elseif j==M && (i~=1 && i~=N) %Top wall
            ae(i,j)=0;
            aw(i,j)=0;
            an(i,j)=0;
            as(i,j)=lambdaP(i,j-1);
            ap(i,j)=as(i,j);
        elseif i==1 && j==1 %Corner SW
            ae(1,1)=0;
            aw(1,1)=0;
            an(1,1)=0;
            as(1,1)=0;
            ap(1,1)=1;
        elseif i==N && j==1 %Corner SE
            ae(N,1)=0;
            aw(N,1)=0;
            an(N,1)=0;
            as(N,1)=0;
            ap(N,1)=1;
        elseif i==1 && j==M %Corner NW
            ae(1,M)=0;
            aw(1,M)=0;
            an(1,M)=0;
            as(1,M)=0;
            ap(1,M)=0;
        elseif i==N && j==M %Corner NE
            ae(N,M)=0;
            aw(N,M)=0;
            an(N,M)=0;
            as(N,M)=0;
            ap(N,M)=1;
        else %General CV
            ae(i,j)=lambdaP(i+1,j);
            aw(i,j)=lambdaP(i-1,j);
            an(i,j)=lambdaP(i,j+1);
            as(i,j)=lambdaP(i,j-1);
            ap(i,j)=ae(i,j)+aw(i,j)+an(i,j)+as(i,j)+2*rhoP(i,j)*deltax^2*cpP(i,j)/deltat;
        end
    end
end




%%

Taux=T(:,:,1);
for t=1:(time_span/deltat+1)

    error=inf;
    %Bot and left BC
    T(:,1 ,t)=Tbot+273;
    T(N,:,t)=8+0.005*(t*deltat-1)+273;

    
    %Coeficient b determination (time dependent)
    %[b] = bp(T(:,:,t));


    for i=1:N
        for j=1:M

            if i==1 && (j~=1 && j~=M) %Left wall ap calculation
                ap(i,j)=ae(i,j)+alfa_ext*(T(i,j,t)-Text);
            elseif i==1 && j==M %SW Corner
                ap(1,M)=alfa_ext*(T(i,j-1,t)-Text);
            end


            if i==1 && (j~=1 && j~=M) %Left wall
                b(i,j)=alfa_ext*deltay*Text;%((rhoP(i,j)*cpP(i,j)*deltax*deltay)/deltat)*Taux(i,j)+alfa_ext*deltay*Text;
            elseif i==N && (j~=0 && j~=M) %Right wall
                b(i,j)=8+0.005*(t*deltat-1)+273;
            elseif j==1 && (i~=1 && i~=N) %Bottom wall
                b(i,j)=Tbot;
            elseif j==M && (i~=1 && i~=N) %Top wall
                b(i,j)=qw*deltax;
            elseif i==1 && j==1 %Corner SW
                b(i,j)=Tbot;
            elseif i==N && j==1 %Corner SE
                b(i,j)=Tbot;
            elseif i==1 && j==M %Corner NW
                b(i,j)=qw*deltax;
            elseif i==N && j==M %Corner NE
                b(i,j)=8+0.005*(t*deltat-1)+273;
            else %general nodes
                b(i,j)=((rhoP(i,j)*cpP(i,j)-deltax*deltay)/deltat)*Taux(i,j);
            end

        end
    end



    

    while error>tol
        for i=1:N
            for j=1:M
                
                %[T(i,j,t+1)]= CN(ap(i,j),aw(i,j),ae(i,j),an(i,j),as(i,j),b(i,j),Taux(i-1,j),Taux(i+1,j),Taux(i,j+1),Taux(i,j-1));
                
                
                if j==1 && i~=1 && i~=N
                    T(i,j,t+1)=(an(i,j)*Taux(i,j+1)+aw(i,j)*Taux(i-1,j)+ae(i,j)*Taux(i+1,j)+b(i,j))/ap(i,j);
                elseif j==M && i~=1 && i~=N
                    T(i,j,t+1)=(as(i,j)*Taux(i,j-1)+aw(i,j)*Taux(i-1,j)+ae(i,j)*Taux(i+1,j)+b(i,j))/ap(i,j);
                elseif i==1 && j~=1 && j~=M
                    T(i,j,t+1)=(an(i,j)*Taux(i,j+1)+as(i,j)*Taux(i,j-1)+ae(i,j)*Taux(i+1,j)+b(i,j))/ap(i,j);
                elseif i==N && j~=1 && j~=M
                    T(i,j,t+1)=(an(i,j)*Taux(i,j+1)+as(i,j)*Taux(i,j-1)+aw(i,j)*Taux(i-1,j)+b(i,j))/ap(i,j);
                elseif (i==1 && j==1) || (i==1 && j==M) || (i==N && j==1) || (i==N && j==M)
                    T(i,j,t+1)=0;
                else
                    T(i,j,t+1)=(an(i,j)*Taux(i,j+1)+as(i,j)*Taux(i,j-1)+ae(i,j)*Taux(i+1,j)+aw(i,j)*Taux(i-1,j)+b(i,j))/ap(i,j);
                end

            end
        end
        
        T(1,1,t+1)=0.5*(T(1,2,t+1)+T(2,1,t+1));
        T(N,1,t+1)=0.5*(T(N-1,1,t+1)+T(N,2,t+1));
        T(1,M,t+1)=0.5*(T(2,M,t+1)+T(1,M-1,t+1));
        T(N,M,t+1)=0.5*(T(N-1,M,t+1)+T(N,M-1,t+1));


        error=max(abs(T(:,:,t+1)-Taux(:,:)));
        Taux=T(:,:,t+1);
 
        
    end
    

t
end




