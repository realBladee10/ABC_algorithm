function [freq]=oct_fraction(fp,fk,fraction)


% f=linspace(0,10000,length(TG_J));


N=fraction; %1/N oktawy
n=-10*N:1:10*N;

    fc(1,:)=1000*(2.^(n/N)); %czêstotliwoœci srodkowe
    fd=2^(1/(2*N));
    fc(2,:)=fc(1,:)./fd; %czestotliwoœæ dolna pasma
    fc(3,:)=fc(1,:)*fd; %czestotliwoœæ górna pasma


freq1=fc(1,:);
freq=freq1(and(freq1>=fp/(2.^(1/(N*9))),freq1<=fk*(2.^(1/(N*9)))));



