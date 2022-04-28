function ATauc = Tau_c(ASTF,delta)
% This function is coded to calculate the apparent cheracteristic duration
% ASTF is the apparent source time function, delta is the sampling rate

ASTF_int=zeros(length(ASTF),1);
tt = 0:delta:delta*(length(ASTF)-1);
for i=2:length(ASTF)
    ASTF_int(i)=ASTF_int(i-1)+(ASTF(i-1)+ASTF(i))/2;
end
ASTF_int=ASTF_int*delta;

ASTF=ASTF/ASTF_int(end);
for i=1:length(ASTF)
    temp1(i)=tt(i)*ASTF(i);
    temp2(i)=tt(i)*tt(i)*ASTF(i);
end

miu=zeros(length(ASTF),1);
miu2=zeros(length(ASTF),1);
for i=2:length(ASTF)
    miu(i)=miu(i-1)+(temp1(i-1)+temp1(i))/2;
    miu2(i)=miu2(i-1)+(temp2(i-1)+temp2(i))/2;
end
miu=miu*delta;miu2=miu2*delta;
ATauc=2*sqrt(miu2(end)-miu(end)*miu(end));