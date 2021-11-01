clc
clear all
H00=0;H01=1;dt=0.0001;E=0.1;gr=inv(E+j*dt-H00);
for i=1:10
    C(1,i)=-imag(gr);
    C(2,i)=real(gr);
    gr=inv(E+j*dt-H00-gr);
end
hold on;
plot(C(1,:),'r');plot(C(2,:),'b')
