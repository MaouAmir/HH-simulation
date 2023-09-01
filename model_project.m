ro = 28.09;
gglia = 66.67;
gamma = 1.86;
beta = 7;
gcl = 0.02;
gnaf = 24;
gca  = 0.08;
gh = 0.05;
gkdr = 3;
gkl = 0.02;
gka = 0.25;
gkm = 1;
gkca = 0.55;
gnal = 0.07;
cao=1.2;
epsilonk=1.4;
kbath = 4;
timestep=0.01;
vm = zeros(1,1000000);
vm(1:150000)=-78;
n(1)=0;
z(1)=0;
b(1)=0;
h(1)=0;
s(1)=0;
r(1)=0;
c(1)=0;
cai(1)=0.001;
nai(1)=15;
ko(1)=5;
ta = 1000;
cainf = 50;
tas = 1;
ik(1) = 0;
iglia(1) = 0 ;
idiff(1)=0;
ipump(1)=0;
ina(1)=0;
I_stim = zeros(1000000,1);
I_stim(150000:850000)= -50;
taca = 1;% not mentioned in article
%
for i=2:1000000
sinf(i) = 1 / (1 + exp((-20-vm(i))/10));
s(i) = ((sinf(i) - s(i-1)) /tas)*timestep + s(i-1);
vca(i) = 26.64*log(cao/cai(i-1));%loop 1
ica(i) = gca * s(i)^2 * (vm(i) - vca(i));
cai(i) = (1/ta * (-1*gamma*ica(i) + (cainf-cai(i-1))/taca))*timestep+cai(i-1);
ninf(i) = 1 / (1 + exp((-30-vm(i))/9.5));
zinf(i) = 1 / (1 + exp((-39-vm(i))/5));
binf(i) = 1 / (1 + exp((-80-vm(i))/6));
ainf(i) = 1 / (1 + exp((-50-vm(i))/20));
minf(i) = 1 / (1 + exp((-30-vm(i))/9.5));
hinf(i) = 1 / (1 + exp((-53-vm(i))/-7));
rinf(i) = 1 / (1 + exp((-84-vm(i))/10.2));

cinf(i) = 1 / (1 + 0.03/48*cai(i)^2);


%
tan(i) = 0.37 + 1.85 * (1/(1+ exp((vm(i)+27)/15)));
tac(i) = 0.2148 / (48*c(i-1)^2 + 0.03);%loop 2
taz = 75;
tab = 15;
tah(i) = 0.37 + 2.78 * (1/(1+ exp((vm(i)+40.5)/6)));

tar(i) = 1/(exp(-14.59-0.086*vm(i))+exp(-1.87 + 0.0701*vm(i)));
%
n(i) = ((ninf(i) - n(i-1)) /tan(i))*timestep + n(i-1);
z(i) = ((zinf(i) - z(i-1)) /taz)*timestep + z(i-1);
b(i) = ((binf(i) - b(i-1)) /tab)*timestep + b(i-1);
h(i) = ((hinf(i) - h(i-1)) /tah(i))*timestep + h(i-1);

r(i) = ((rinf(i) - r(i-1)) /tar(i))*timestep + r(i-1);
c(i) = ((cinf(i) - c(i-1)) /tac(i))*timestep + c(i-1);
%
ko(i) = (1/ta * (gamma*beta*ik(i-1)-2*beta*gamma*ipump(i-1)-iglia(i-1)-idiff(i-1)))*timestep + ko(i-1);%third and forth anf fith and sixth loop
nai(i) = (1/ta * (-1*gamma*ina(i-1)-3-gamma*ipump(i-1)))*timestep + nai(i-1);%seventh loop

ki(i) = 140 + (18 - nai(i));
cli(i) = nai(i) + ki(i) + 2*cai(i) -150;
nao(i) = 144 + beta*(nai(i)-18);
clo(i) = nao(i) + ko(i) + 2*cao;




%
vna(i) = 26.64*log(nao(i)/nai(i));
vk(i) = 26.64*log(ko(i)/ki(i));
vh(i) = 26.64*log((0.2*nao(i)+ ko(i))/(0.2*nai(i) + ki(i)));
vcl(i) = 26.64*log(clo(i)/cli(i));


% 
idiff(i) = epsilonk*(ko(i)-kbath);
iglia(i) = gglia/(1+exp((18-ko(i))/2.5));
ipump(i) = ro/((1+exp(25-nai(i)/3))*(1+exp(5.5-ko(i))));




ikdr(i) = gkdr * n(i)^4 * (vm(i) - vk(i));
ikm(i) = gkm * z(i)* (vm(i) - vk(i));
ika(i) = gka * ainf(i)^3 * (vm(i) - vk(i));
ikca(i) = gkca * c(i)^2 * (vm(i) - vk(i));
ikl(i) = gkl * (vm(i) - vk(i));
inaf(i) = gnaf * minf(i)^3 * h(i)  * (vm(i) - vna(i));
inal(i) = gnal * (vm(i) - vna(i));
ih(i) = gh * r(i) * (vm(i) - vh(i));
icl(i) = gcl * (vm(i) - vcl(i));


ik(i)= ikdr(i)+ ikm(i) + ika(i) + ikl(i) + ikca(i) ; 
ina(i) = inaf(i)+ inal(i);
I_m(i) = -(ik(i)+ina(i)+ih(i)+icl(i)+ica(i));
%

vm(i+1) = (I_m(i) + I_stim(i) + ipump(i)/1.86) * timestep / gamma +vm(i) ; 
end
%%
plot (vm);