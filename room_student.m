
function [stdLp,sLp3,sLp12]=room_student(dane)


r0=dane(1:3);
r=dane(4:6);
GF=GreenFunction_OK(r0,r);


Lp=(20*log10(abs(1.21*GF)/2e-5));
sLp3=smooth(Lp,3);
sLp12=smooth(Lp,12);
stdLp=std(sLp3)+std(smooth(Lp,12))*3;

