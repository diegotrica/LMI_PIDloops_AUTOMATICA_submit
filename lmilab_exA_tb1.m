clear

myfile{1}='github\caseA\lmilab_exA_SOFS_L2N_Nl=20_PhiD=100_gL2=1000.mat';
myfile{2}='github\caseA\lmilab_exA_SOFS_L2N_Nl=20_PhiD=100_gL2=100.mat';
myfile{3}='github\caseA\lmilab_exA_SOFS_L2N_Nl=20_PhiD=100_gL2=10.mat';
myfile{4}='github\caseA\lmilab_exA_SOFS_L2N_Nl=20_PhiD=100_gL2=1.mat';

myfile{5}='github\caseA\lmilab_exA_SOFS_L2N_Nl=20_PhiD=10_gL2=1000.mat';
myfile{6}='github\caseA\lmilab_exA_SOFS_L2N_Nl=20_PhiD=10_gL2=100.mat';
myfile{7}='github\caseA\lmilab_exA_SOFS_L2N_Nl=20_PhiD=10_gL2=10.mat';
myfile{8}='github\caseA\lmilab_exA_SOFS_L2N_Nl=20_PhiD=10_gL2=1.mat';

tstart=0; tend=300; tstep=0.1; t=(tstart:tstep:tend)'; %time domain
yf=[0.1;5]; dw=[0.1;5];

Nrlt=0;
save('github\caseA\lmilab_exA_tb1_rlt.mat');
while Nrlt<length(myfile);
rlt(Nrlt+1).data=myfile{Nrlt+1};
load(rlt(Nrlt+1).data);
rlt(Nrlt+1).gL2=gL2;
rlt(Nrlt+1).PhiD=nthroot(det(PhiD),Ny);

%This work proposed technique PI control
rlt(Nrlt+1).TW.PI.KP=KPpi;
rlt(Nrlt+1).TW.PI.KI=KIpi;
rlt(Nrlt+1).TW.PI.TI=TIpi;
rlt(Nrlt+1).TW.PI.alp=alppi;
rlt(Nrlt+1).TW.PI.ElapTime=elapsedTime_pi;

K_=[KPpi,KIpi];
rlt(Nrlt+1).TW.PI.Gload=[];
rlt(Nrlt+1).TW.PI.Gsp=[];
rlt(Nrlt+1).TW.PI.eig=[];
rlt(Nrlt+1).TW.PI.hinf=[];
rlt(Nrlt+1).TW.PI.info=[];
if ~isempty(alppi)
rlt(Nrlt+1).TW.PI.Gload=ss(A_pi-Bu_pi*K_*Cy_pi, Bw_pi, Cy_pi(1:Ny,:), 0);
Br_pi=[Bu*KPpi;-eye(Ny)];
rlt(Nrlt+1).TW.PI.Gsp=ss(A_pi-Bu_pi*K_*Cy_pi, Br_pi, Cy_pi(1:Ny,:), 0);
rlt(Nrlt+1).TW.PI.eig=eig(rlt(Nrlt+1).TW.PI.Gload);
rlt(Nrlt+1).TW.PI.hinf=hinfnorm(rlt(Nrlt+1).TW.PI.Gload);
rlt(Nrlt+1).TW.PI.info=stepinfo_caseA(KPpi,TIpi,[],[],par,t,x0,u0,[yf;dw],'PI');
end

%This work proposed technique PID control
rlt(Nrlt+1).TW.PID.KP=KPpid;
rlt(Nrlt+1).TW.PID.KI=KIpid;
rlt(Nrlt+1).TW.PID.TI=TIpid;
rlt(Nrlt+1).TW.PID.KD=KDpid;
rlt(Nrlt+1).TW.PID.TD=TDpid;
rlt(Nrlt+1).TW.PID.alp=alppid;
rlt(Nrlt+1).TW.PID.ElapTime=elapsedTime_pid;

K_=[KPpid,KIpid,KDpid];
rlt(Nrlt+1).TW.PID.Gload=[];
rlt(Nrlt+1).TW.PID.Gsp=[];
rlt(Nrlt+1).TW.PID.eig=[];
rlt(Nrlt+1).TW.PID.hinf=[];
rlt(Nrlt+1).TW.PID.info=[];
if ~isempty(alppid)
rlt(Nrlt+1).TW.PID.Gload=ss(A_pid-Bu_pid*K_*Cy_pid, Bw_pid, Cy_pid(1:Ny,:), 0);
Br_pid=[Bu*KPpid;-eye(Ny);-PhiD*Cy*Bu*KPpid];
rlt(Nrlt+1).TW.PID.Gsp=ss(A_pid-Bu_pid*K_*Cy_pid, Br_pid, Cy_pid(1:Ny,:), 0);
rlt(Nrlt+1).TW.PID.eig=eig(rlt(Nrlt+1).TW.PID.Gload);
rlt(Nrlt+1).TW.PID.hinf=hinfnorm(rlt(Nrlt+1).TW.PID.Gload);
rlt(Nrlt+1).TW.PID.info=stepinfo_caseA(KPpid,TIpid,TDpid,PhiD,par,t,x0,u0,[yf;dw],'PID');
end

%ILMI technique PI control
rlt(Nrlt+1).ILMI.PI.KP=KPpi_ILMI;
rlt(Nrlt+1).ILMI.PI.KI=KIpi_ILMI;
rlt(Nrlt+1).ILMI.PI.TI=TIpi_ILMI;
rlt(Nrlt+1).ILMI.PI.alp=alppi_ILMI;
rlt(Nrlt+1).ILMI.PI.iter=iterpi_ILMI;
rlt(Nrlt+1).ILMI.PI.ElapTime=elapsedTime_pi_ILMI;

K_=[KPpi_ILMI,KIpi_ILMI];
rlt(Nrlt+1).ILMI.PI.Gload=[];
rlt(Nrlt+1).ILMI.PI.Gsp=[];
rlt(Nrlt+1).ILMI.PI.eig=[];
rlt(Nrlt+1).ILMI.PI.hinf=[];
rlt(Nrlt+1).ILMI.PI.info=[];
if ~isempty(alppi_ILMI)
rlt(Nrlt+1).ILMI.PI.Gload=ss(A_pi-Bu_pi*K_*Cy_pi, Bw_pi, Cy_pi(1:Ny,:), 0);
Br_pi_ILMI=[Bu*KPpi_ILMI;-eye(Ny)];
rlt(Nrlt+1).ILMI.PI.Gsp=ss(A_pi-Bu_pi*K_*Cy_pi, Br_pi_ILMI, Cy_pi(1:Ny,:), 0);
rlt(Nrlt+1).ILMI.PI.eig=eig(rlt(Nrlt+1).ILMI.PI.Gload);
rlt(Nrlt+1).ILMI.PI.hinf=hinfnorm(rlt(Nrlt+1).ILMI.PI.Gload);
rlt(Nrlt+1).ILMI.PI.info=stepinfo_caseA(KPpi_ILMI,TIpi_ILMI,[],[],par,t,x0,u0,[yf;dw],'PI');
end

%ILMI technique PID control
rlt(Nrlt+1).ILMI.PID.KP=KPpid_ILMI;
rlt(Nrlt+1).ILMI.PID.KI=KIpid_ILMI;
rlt(Nrlt+1).ILMI.PID.TI=TIpid_ILMI;
rlt(Nrlt+1).ILMI.PID.KD=KDpid_ILMI;
rlt(Nrlt+1).ILMI.PID.TD=TDpid_ILMI;
rlt(Nrlt+1).ILMI.PID.alp=alppid_ILMI;
rlt(Nrlt+1).ILMI.PID.iter=iterpid_ILMI;
rlt(Nrlt+1).ILMI.PID.ElapTime=elapsedTime_pid_ILMI;

K_=[KPpid_ILMI,KIpid_ILMI,KDpid_ILMI];
rlt(Nrlt+1).ILMI.PID.Gload=[];
rlt(Nrlt+1).ILMI.PID.Gsp=[];
rlt(Nrlt+1).ILMI.PID.eig=[];
rlt(Nrlt+1).ILMI.PID.hinf=[];
rlt(Nrlt+1).ILMI.PID.info=[];
if ~isempty(alppid_ILMI)
rlt(Nrlt+1).ILMI.PID.Gload=ss(A_pid-Bu_pid*K_*Cy_pid, Bw_pid, Cy_pid(1:Ny,:), 0);
Br_pid_ILMI=[Bu*KPpid_ILMI;-eye(Ny);-PhiD*Cy*Bu*KPpid_ILMI];
rlt(Nrlt+1).ILMI.PID.Gsp=ss(A_pid-Bu_pid*K_*Cy_pid, Br_pid_ILMI, Cy_pid(1:Ny,:), 0);
rlt(Nrlt+1).ILMI.PID.eig=eig(rlt(Nrlt+1).ILMI.PID.Gload);
rlt(Nrlt+1).ILMI.PID.hinf=hinfnorm(rlt(Nrlt+1).ILMI.PID.Gload);
rlt(Nrlt+1).ILMI.PID.info=stepinfo_caseA(KPpid_ILMI,TIpid_ILMI,TDpid_ILMI,PhiD,par,t,x0,u0,[yf;dw],'PID');
end
Nrlt=Nrlt+1;
save('github\caseA\lmilab_exA_tb1_rlt.mat','rlt','-append');
load('github\caseA\lmilab_exA_tb1_rlt.mat','myfile');
end

%Reload results
clear
load('github\caseA\lmilab_exA_tb1_rlt.mat','rlt');