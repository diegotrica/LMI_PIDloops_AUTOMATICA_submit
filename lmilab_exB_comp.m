clear, close all
load('github\caseB\lmilab_exB_SOFS_L2N_Nd=5_Nl=7_PhiD=100_gL2=1000.mat');

tstart=0; tend=300; tstep=0.1; t=(tstart:tstep:tend)'; %time domain
yf=[0.1;0.1]; dw=[0.245];

%Tuning from this work
KP_TW=KPpid; KI_TW=KIpid; KD_TW=KDpid;

%Tuning from Boyd, S., Hast, M., Astrom, J. (2016)
KP_BHA16=...
	[0.1535 	  0
	 0      -0.0692];
	 
KI_BHA16=...
	[0.0210 	  0
	 0      -0.0136];

KD_BHA16=...
	[0.1714 	  0
	 0      -0.1725];
	 
%Tuning from Chen, D., Seaborg, D. E. (2003)
KP_CS03=...
	[0.436 	  0
	 0      -0.0945];
	 
KI_CS03=...
	[11\0.436 	  0
	 0      -15.5\0.0945];
	 
%BLT tuning Luyben, W. L. (1986)
KP_BLT=...
	[0.375 	  0
	 0      -0.075];
	 
KI_BLT=...
	[8.29\0.375 	  0
	 0      -23.6\0.075];	 

Br_=[Bu*KP_TW;
	 -eye(Ny);
	 -PhiD*Cy*Bu*KP_TW];
Bw_=[Bw
	 zeros(Ny,Nw)
	 -PhiD*Cy*Bw];
sys_TW=ss(A_pid-Bu_pid*[KP_TW,KI_TW,KD_TW]*Cy_pid, [Br_,Bw_], Cy_pid(1:Ny,:), []);
set(sys_TW,'Name','WB model TW PID tuning');
set(sys_TW,'InputName',{'w_{D,r}','w_{B,r}','f'});
set(sys_TW,'OutputName',{'w_D','w_B'});
set(sys_TW,'TimeUnit','minutes');

Br_=[Bu*KP_BHA16;
	 -eye(Ny);
	 -PhiD*Cy*Bu*KP_BHA16];
Bw_=[Bw
	 zeros(Ny,Nw)
	 -PhiD*Cy*Bw];
sys_BHA16=ss(A_pid-Bu_pid*[KP_BHA16,KI_BHA16,KD_BHA16]*Cy_pid, [Br_,Bw_], Cy_pid(1:Ny,:), []);
set(sys_BHA16,'Name','WB model BHA16 PID tuning');
set(sys_BHA16,'InputName',{'w_{D,r}','w_{B,r}','f'});
set(sys_BHA16,'OutputName',{'w_D','w_B'});
set(sys_BHA16,'TimeUnit','minutes');

Br_=[Bu*KP_CS03;
	 -eye(Ny)];
Bw_=[Bw
	 zeros(Ny,Nw)];
sys_CS03=ss(A_pi-Bu_pi*[KP_CS03,KI_CS03]*Cy_pi,[Br_,Bw_], Cy_pi(1:Ny,:), []);
set(sys_CS03,'Name','WB model CS03 PID tuning');
set(sys_CS03,'InputName',{'w_{D,r}','w_{B,r}','f'});
set(sys_CS03,'OutputName',{'w_D','w_B'});
set(sys_CS03,'TimeUnit','minutes');

Br_=[Bu*KP_BLT;
	 -eye(Ny)];
Bw_=[Bw
	 zeros(Ny,Nw)];
sys_BLT=ss(A_pi-Bu_pi*[KP_BLT,KI_BLT]*Cy_pi,[Br_,Bw_], Cy_pi(1:Ny,:), []);
set(sys_BLT,'Name','WB model BLT PID tuning');
set(sys_BLT,'InputName',{'w_{D,r}','w_{B,r}','f'});
set(sys_BLT,'OutputName',{'w_D','w_B'});
set(sys_BLT,'TimeUnit','minutes');

%Plot step response
figh=figure(1);
dataopt=stepDataOptions('StepAmplitude',[yf;dw]);
plotopt=timeoptions;
plotopt.Title.String='';
plotopt.XLabel.String='time';
plotopt.XLabel.FontSize=12;
plotopt.TimeUnits='minutes';
plotopt.YLabel.String='';
plotopt.OutputLabels.FontSize=10;
plotopt.OutputLabels.Color=[0 0 0];
plotopt.InputLabels.FontSize=10;
plotopt.InputLabels.Color=[0 0 0];
plotopt.TickLabel.FontSize=8;

%Rendering figure
figh.Units='centimeters';
figh.Color='w';
figh.OuterPosition=[0 0 16 12];

stepplot(sys_TW,t,plotopt,dataopt);
hold on; stepplot(sys_BHA16,t,plotopt,dataopt);
hold on; stepplot(sys_CS03,t,plotopt,dataopt);
hold on; stepplot(sys_BLT,t,plotopt,dataopt);

%Set legend
h=legend('This work','BHA16','CS03','BLT');
%h.Position=[0.9,0.820,0.1,0.1];
h.Position=[0.35,0.925,0.3,0.1];
h.Orientation='horizontal';

%Export figure
print(figh,'C:\Users\HP\Dropbox\Livro_LMIChE\paper\fig4','-dsvg');
print(figh,'C:\Users\HP\Dropbox\Livro_LMIChE\paper\fig4','-deps');

%ITAE comparison
info_TW=stepinfo_caseB(KP_TW,KI_TW,KD_TW,PhiD,ss(A_pid,[Bu_pid,Bw_pid],Cy_pid,0),t,[yf;dw],'PID');
info_BHA16=stepinfo_caseB(KP_BHA16,KI_BHA16,KD_BHA16,PhiD,ss(A_pid,[Bu_pid,Bw_pid],Cy_pid,0),t,[yf;dw],'PID');
info_CS03=stepinfo_caseB(KP_CS03,KI_CS03,[],[],ss(A_pi,[Bu_pi,Bw_pi],Cy_pi,0),t,[yf;dw],'PI');
info_BLT=stepinfo_caseB(KP_BLT,KI_BLT,[],[],ss(A_pi,[Bu_pi,Bw_pi],Cy_pi,0),t,[yf;dw],'PI');