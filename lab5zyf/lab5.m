clear all
format long
L = 0.01;
R = 99.94;
C = 233.02*10^(-9);
f = [0.316968
   
   0.476065
   
   0.696943
   
   0.757990
   
   1.170060
   
   1.337000
   
   1.487060
   
   1.757000
   
    1.947010
   
   2.017010
   
   2.317020
   
   2.377020
   
   2.467020
   
   2.607010
   
    2.847000
   
   3.097020
   
   3.307050
   
   3.587020
   
   3.967080
   
    7.027010
   
   10.027000
   
   ]';
f = f*1000;
v=[0.40 
0.64 
0.96 
1.12 
1.76 
2.08 
2.40 
3.04 
3.36 
3.52 
3.84 
3.76 
3.68 
3.52 
3.20 
2.88 
2.72 
2.40 
2.08 
1.12 
0.80 
]';
f1=1.75700*1000;
f2=3.30705*1000;
f0=2.31702*1000;
muR=0.01;
muf=0.001;
muepsilon=0.001;
muC=0.01*10^(-9);
muU=0.001;

%--------------------------------------------------------------

[~,i]=max(v);

ff0=f./f(i);
vvm=v./max(v);
fitheo = atan((2*pi.*f*L - 1./(2*pi.*f*C))/R);
fiex=acos(vvm);
Qex=f0/(f2-f1);
Qth=sqrt(L*C)/(R*C);
Qth2=2*pi*f0*L/R;
Qth3=1/(2*pi*f0*R*C);
um=max(v);


muIIm=sqrt((muU/um)^2+(-(v./um.^2).*muU).^2);
muff0=sqrt((muf/f0)^2+(-f./f0^2.*muf).^2);
mufiex=sqrt((-muU./sqrt(um^2-v.^2)).^2+(v./(um.*sqrt(um^2-v.^2)).*muU).^2);
mufiex1=R*(2*pi*L+1./(2*pi.*f.^2.*C))./(R^2+(2*pi.*f.*L-1./(2*pi.*f.*C)).^2).*muf;
mufiex2=R./((2*pi.*f.*L-1./(2*pi.*f.*C)).^2+R^2).*muC;
mufiex3=(1./(2*pi.*f.*C)-2*pi.*f.*L)./((2*pi.*f.*L-1./(2*pi.*f.*C)).^2+R^2).*muR;
mufitheo=sqrt(mufiex1.^2+mufiex2.^2+mufiex3.^2);
muqex=sqrt((muf/(f2-f1))^2+(f0/(f2-f1)^2*muf)^2+(-f0/(f2-f1)^2*muf)^2);
muqtheo=sqrt((-sqrt(L*C)/(R^2*C)*muR)^2+(-sqrt(L)/(2*R*C^(-3/2))*muC)^2);
final=[f',v',vvm',muIIm',ff0',muff0',fiex',mufiex',fitheo',mufitheo'];
csvwrite('1.csv',final);
Qex
muqex
Qth
muqtheo
