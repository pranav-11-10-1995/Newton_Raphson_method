clc;   
      %s.no from   to    admit         admit(1/2 line)
data=[ 1    1      2     2.035-8.61i      0.02i;
       2    2      3     2.035-8.61i      0.02i;
       3    3      1     2.035-8.61i      0.02i;   ];  
   y=zeros(3,3);% returns n*n matrix of zeros
   len=length(data(:,1)); % length finds no. of elements in a list
   
for i=1:len
    for j=1:len
        m=data(i,2);
        n=data(j,3);
        if(m==n)% m=n -> compares elements in 2nd & 3rd columns
        y(m,n)= data(i,4)+data(j,4)+data(i,5)+data(j,5);
        elseif(i==j)% i=j -> compares the elements in same row
        y(m,n)=-data(i,4);
        else
        y(m,n)=y(n,m);% since y(1,2),y(2,3) occur before y(2,1) & y(3,2)
        end
    end
 end
y(n,m)=y(m,n);% since y(1,3) occurs before y(3,1) ie., y(3,1) -> y(1,3)
disp('          FORMING Y BUS')
y  
disp('    SEPARATING REAL & IMAGINARY PARTS OF Y BUS')
disp('          REAL PART OF Y BUS')
g=real(y)
disp('          IMAGINARY PART OF Y BUS')
b=imag(y)
for i=1:len
    for j=1:len
 mag(i,j)=sqrt(real(y(i,j))*real(y(i,j))+imag(y(i,j))*imag(y(i,j)));
    end 
end
disp('          MAGNITUDE MATRIX OF Y BUS')
mag
for i=1:len
    for j=1:len
        if i==j
        del(i,j)=atand(imag(y(i,j))/real(y(i,j)));
        else
        del(i,j)=180+atand(imag(y(i,j))/real(y(i,j)));
        end
     end 
end
disp('          DELTA MATRIX OF Y BUS')
del
        %s       %pg        qg         Pl       Ql      v          del                
data1=[  1      1.00         0.5       0        0       1.03       0      ;   %slack bus 
         2      1.5          0         0        0       1.03       0      ;   %PV bus
         3      0            0        -1.2     -0.5     1          0      ]; %PQ bus
delp2=1;
while(delp2>0.01)
p2=0;
q=0;
i=2;
disp('          FINDING P2,Q2,P3,Q3')
disp('          FOR PV BUS-> P2,Q2')
for j=1:3
    p2=p2+data1(i,6)*data1(j,6)*mag(i,j)*cosd(del(i,j)+data1(j,7)-data1(i,7));
    q=q+data1(i,6)*data1(j,6)*mag(i,j)*sind(del(i,j)+data1(j,7)-data1(i,7));
end
  p2
  q2=-q
 p3=0;
 q=0;
 i=3;
 disp('          FOR PQ BUS-> P3,Q3')
for j=1:3
    p3=p3+data1(i,6)*data1(j,6)*mag(i,j)*cosd(del(i,j)+data1(j,7)-data1(i,7));
    q=q+data1(i,6)*data1(j,6)*mag(i,j)*sind(del(i,j)+data1(j,7)-data1(i,7));
end
p3
q3=-q

disp('          FORMING JACOBIAN')
h22=-q2-data1(2,6)^2*b(2,2);%H22=-Q2-|V2|^2*B22
h33=-q3-data1(3,6)^2*b(3,3);%H33=-Q3-|V3|^2*B33
l33=+q3-data1(3,6)^2*b(3,3);%L33=Q3-|V3|^2*B33
j33=+p3-data1(3,6)^2*g(3,3);%J33=P3-|V3|^2*G33
n33=+p3+data1(3,6)^2*g(3,3);%N33=P3+|V3|^2*G33
i=2;
j=3;
r0=data1(i,6)*data1(j,6)*mag(i,j);%R0=|V2|*|V3|*|Y23|
s0=data1(j,6)*data1(i,6)*mag(j,i);%S0=|V3|*|V2|*|Y32|
%all elements with suffixes 23 or 32 have -ve sign
h23=-r0*sind(del(i,j)+data1(j,7)-data1(i,7));% DEL-> YBUS DELTA
h32=-s0*sind(del(j,i)+data1(i,7)-data1(j,7));
j23=-r0*cosd(del(i,j)+data1(j,7)-data1(i,7));
j32=-s0*cosd(del(j,i)+data1(i,7)-data1(j,7));
n23=-j23;
jac=[h22 h23 n23
     h32 h33 n33
     j32 j33 l33]
delp2=data1(2,2)-p2;
delp3=data1(3,4)-p3;
delq3=data1(3,5)-q3;
disp('          CALCULATING RESIDUE')
disp('     TARGET: MAKING THE RESIDUE ZERO')
residue=[delp2
         delp3
         delq3]
change=inv(jac)*residue;
data1(2,7)=data1(2,7)+change(1,1);
data1(3,7)=data1(3,7)+change(2,1);
data1(3,6)=data1(3,6)*(1+change(3,1));
end

disp('          GETTING V FOR PQ BUS & DELTA FOR PV & PQ BUSES')
del2=data1(2,7)+change(1,1)
del3=data1(3,7)+change(2,1)
v3=data1(3,6)*(1+change(3,1))










