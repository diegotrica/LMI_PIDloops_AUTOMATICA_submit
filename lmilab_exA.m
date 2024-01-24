clear
%--------------------------------------------------------------------
%VAN DE VUSSE NON-ISOTHERMAL CSTR MODEL (1993) ----------------------
%Parameters:
par(1)	=  60\1.287e+12;%K1o: Arrhenius spec. reaction rate  	[1/min]
par(2)	=  60\1.287e+12;%K2o: Arrhenius spec. reaction rate  	[1/min]
par(3)	=  60\9.043e+09;%K3o: Arrhenius spec. reaction rate  	[m^3/kmol min]
par(4)	=  9758.3;      %E1R: Activion energy in terms of temp. [K]
par(5)	=  9758.3;      %E2R: Activion energy in terms of temp. [K]
par(6)	=  8560;        %E3R: Activion energy in terms of temp. [K]
par(7)	=  4200;        %DH1: Molar reaction heat 			  	[kJ/kmol]
par(8)	=  -11000;      %DH2: Molar reaction heat 			  	[kJ/kmol]
par(9)	=  -41850;      %DH3: Molar reaction heat 			  	[kJ/kmol]
par(10)	=  934.2;       %rho: liquid mass density 			  	[kg/m^3]
par(11)	=  3.01;        %cp:  liquid specific heat capacity 	[kJ/kg K]
par(12)	=  60\4032;     %Uf:  Global heat transfer coeff.    	[kJ/min m^2 K]
par(13)	=  0.215;       %A:	  Heat transfer area      		  	[m^2]
par(14)	=  0.01;        %V:	  liquid reaction volume 			[m^3]
par(15)	=  5.00;        %mc:  cooling medium mass hold-up 		[kg]
par(16)	=  2.0;         %cpc: cooling medium spec. heat cap.	[kJ/kg K]

%----- Rigorous steady-state computation:
%States initial guess
CA0		= 1.235;		%component A molar concentration 			 [kmol/m3]
CB0		= 0.900;		%component B molar concentration 			 [kmol/m3]
T0 		= 134.14+273.15;%reactor temperature						 [K]
Tc0 	= 128.95+273.15;%temperature of cooling medium				 [K]

%Steady-state states allocation
x0=[CA0;CB0;T0;Tc0];

%Input initial guess
th0 	= 60\18.83;     %inlet vol. flow rate per liq. reaction vol. [1/min]
Qc0 	= 60\4496;		%cooling medium heat removal rate			 [kJ/min]
CAin0 	= 5.100;        %component A inlet molar concentration 		 [kmol/m3]
Tin0 	= 130+273.15;   %inlet temperature 					    	 [K]

%Steady-state input allocation
%u0=[th0;Tc0;CAin0;Tin0];
u0=[th0;Qc0;CAin0;Tin0];

%Fixed and free states:
idx 	= [0;1;1;0];
fixed_x	= x0(idx==true);
free_x 	= x0(idx==false);

%Fixed and free inputs:
idu 	= [0;0;1;1];
fixed_u = u0(idu==true);
free_u 	= u0(idu==false);

%Steady-state computation
fh=@(z)sfunc_ss(@vandevusse,z,fixed_x,idx,fixed_u,idu,[],[],par);
opt = optimoptions('fsolve','Display','off');
z = fsolve(fh,[free_x;free_u],opt);
free_x = z(1:size(find(idx==false),1));
free_u = z(size(find(idx==false),1)+1:end);
x0(idx==false)=free_x;
u0(idu==false)=free_u;
%--------------------------------------------------------------------
%OPEN-LOOP SYSTEM MATRICES ------------------------------------------
sys=sfunc2ss(@vandevusse,x0,u0,par);

A=sys.A;
Bu=sys.B(:,1:2);
Bw=sys.B(:,3:4);
Cy=sys.C;

Nx=size(A,1);   %State variables size
Ny=size(Cy,1);  %Controlled output variables size
Nu=size(Bu,2);  %Input variables size
Nw=size(Bw,2);  %Exogenous input variables size

PhiD=1e2*eye(Ny);	%Derivative kick filter
%PhiD=1e1*eye(Ny);	%Derivative kick filter
%--------------------------------------------------------------------

%--------------------------------------------------------------------
%POLYTOPIC LDI SYSTEM -----------------------------------------------
%Nominal LTI system
Poly.A0=A;
Poly.Bu0=Bu;
Poly.Bw0=Bw;
Poly.Cy0=Cy;
Poly.deig0=eig(A);

%Parameters uncertainty
dpar(1)	=  60\0.04e+12;	%K1o: Arrhenius spec. reaction rate     [1/min]
dpar(2) =  60\0.04e+12; %K2o: Arrhenius spec. reaction rate     [1/min]
dpar(3)	=  60\0.27e+09; %K3o: Arrhenius spec. reaction rate     [m^3/kmol min]
dpar(4) =  false;       %E1R: Activion energy in terms of temp. [K]
dpar(5) =  false;       %E2R: Activion energy in terms of temp. [K]
dpar(6) =  false;       %E3R: Activion energy in terms of temp. [K]
dpar(7)	=  2360;        %DH1: Molar reaction heat 			  	[kJ/kmol]
dpar(8)	=  1920;        %DH2: Molar reaction heat 			  	[kJ/kmol]
dpar(9)	=  1410;        %DH3: Molar reaction heat 			  	[kJ/kmol]
dpar(10)=  0.4;         %rho: liquid mass density 			  	[kg/m^3]
dpar(11)=  0.04;        %cp:  liquid specific heat capacity 	[kJ/kg K]
dpar(12)=  60\120;      %Uf:  Global heat transfer coeff.       [kJ/min m^2 K]
dpar(13)=  false;       %A:	  Heat transfer area      		  	[m^2]
dpar(14)=  false;       %V:	  liquid reaction volume 			[m^3]
dpar(15)=  false;       %mc:  cooling medium mass hold-up 		[kg]
dpar(16)=  0.05;        %cpc: cooling medium spec. heat cap.	[kJ/kg K]

Poly.par={}; Poly.Nl=0;
for i=1:length(dpar)
 if dpar(i)~=false
	Poly.par{Poly.Nl+1}=zeros(size(par,1),size(par,2));
	Poly.par{Poly.Nl+1}(i)=dpar(i);
	Poly.Nl=Poly.Nl+1;
	Poly.par{Poly.Nl+1}=zeros(size(par,1),size(par,2));
	Poly.par{Poly.Nl+1}(i)=-dpar(i);
	Poly.Nl=Poly.Nl+1;
 end
end
for i=1:Poly.Nl
	Poly.par{i}=par+Poly.par{i};
end

Poly.Nx=Nx; Poly.Ny=Ny; Poly.Nu=Nu; Poly.Nw=Nw;
Poly.A=cell(Poly.Nl,1);
Poly.Bu=cell(Poly.Nl,1);
Poly.Bw=cell(Poly.Nl,1);
Poly.Cy=cell(Poly.Nl,1);
Poly.deig=cell(Poly.Nl,1);
Poly.x0=cell(Poly.Nl,1);
Poly.u0=cell(Poly.Nl,1);

for i=1:Poly.Nl
%Steady-state for polytopic vertex
fh=@(z)sfunc_ss(@vandevusse,z,fixed_x,idx,fixed_u,idu,[],[],Poly.par{i});
opt = optimoptions('fsolve','Display','off');
z = fsolve(fh,[free_x;free_u],opt);
free_x = z(1:size(find(idx==false),1));
free_u = z(size(find(idx==false),1)+1:end);
Poly.x0{i}=x0; Poly.x0{i}(idx==false)=free_x;
Poly.u0{i}=u0; Poly.u0{i}(idu==false)=free_u;

%Polytopic vertex LTI model
sys=sfunc2ss(@vandevusse,Poly.x0{i},Poly.u0{i},Poly.par{i});

%Polytopic vertex state-space realization
Poly.A{i}=sys.A;
Poly.Bu{i}=sys.B(:,1:2);
Poly.Bw{i}=sys.B(:,3:4);
Poly.Cy{i}=sys.C;
Poly.deig{i}=eig(sys.A);
end

%{
%Uncommenting this section and for LTI system instead of PLDI
clear Poly
Poly.A0=A;
Poly.Bu0=Bu;
Poly.Bw0=Bw;
Poly.Cy0=Cy;
Poly.Nl=1; Poly.Nx=Nx; Poly.Ny=Ny; Poly.Nu=Nu; Poly.Nw=Nw;
Poly.A{1}=A;
Poly.Bu{1}=Bu;
Poly.Bw{1}=Bw;
Poly.Cy{1}=Cy;
Poly.deig{1}=eig(A);
Poly.x0=x0;
Poly.u0=u0;
%}
%-------------------------------------------------------------------- 
%LMI PROBLEM STATEMENT ----------------------------------------------
%control action structure matrix
Tk=...
	[1     0
     0     -1];

%Control parameters sector constratins
KPlb=...                %Proportional gain lower bound
    [0.1 	0
     0  	-0.1];
KPub=...                %Proportional gain upper bound
    [1 		0
     0   	-2];
TIlb=...                %Integral time lower bound
    [1  	0
     0    	1];
TIub=...                %Integral time upper bound
    [10  	0
     0    	10];
TDlb=...                %Derivative time lower bound
    [0  	0
     0    	0];
TDub=...                %Derivative time upper bound
    [10/60  0
     0    	10/60];

%Higher KPub values for re-simulate #6 ILMI,
%which converged on the boundaries
%{
KPlb=...                %Proportional gain lower bound
    [0.1 	0
     0  	-0.1];
KPub=...                %Proportional gain upper bound
    [1 		0
     0   	-4];
TIlb=...                %Integral time lower bound
    [1  	0
     0    	2];
TIub=...                %Integral time upper bound
    [10  	0
     0    	10];
TDlb=...                %Derivative time lower bound
    [0  	0
     0    	0];
TDub=...                %Derivative time upper bound
    [10/60  0
     0    	10/60];
%}

%Upper bound on L2 norm
%gL2=1e3;
%gL2=1e2;
gL2=1e1;
%gL2=1e0;

myfile=strjoin({'lmilab_exA_SOFS_L2N_Nl=',num2str(Poly.Nl),'_PhiD=',num2str(nthroot(det(PhiD),Ny)),'_gL2=',num2str(gL2)},'');
diary(myfile);
%--------------------------------------------------------------------
%AUGMENTED PI OPEN-LOOP SYSTEM MATRICES -----------------------------
%Open-loop autonomous system matrix (Nx+Ny,Nx+Ny)
A_pi= [A              zeros(Nx,Ny)
     Cy             zeros(Ny,Ny)];

%Open-loop system input gain matrix (Nx+Ny,Nu)
Bu_pi=[Bu
	 zeros(Ny,Nu)];
	 
%Open-loop system exogenous input gain matrix (Nx+Ny,Nw)
Bw_pi=[Bw
	 zeros(Ny,Nw)];

%Controlled Output-State matrix relationship (Ny+Ny,Ny+Ny)
Cy_pi=[Cy             zeros(Ny,Ny)
	 zeros(Ny,Nx)   eye(Ny,Ny)];

%--------------------------------------------------------------------
%PI RESULTS ---------------------------------------------------------
%PI tuning by SOF stabilization
disp('----- start of PI calculation ---------------------------------------');
timerVal = tic;
[KPpi,KIpi,~,TIpi,~,alppi,P_pi]=...
    PID_SOFS_L2N(Poly,gL2,Tk,[],KPlb,KPub,TIlb,TIub,[],[],'PI');
elapsedTime_pi = toc(timerVal);

timerVal = tic;
[KPpi_ILMI,KIpi_ILMI,~,TIpi_ILMI,~,alppi_ILMI,P_pi_ILMI,iterpi_ILMI]=...
    PID_SOFS_L2N_ILMI(Poly,gL2,Tk,[],KPlb,KPub,TIlb,TIub,[],[],'PI');
elapsedTime_pi_ILMI = toc(timerVal);
disp('----- end of PI calculation -----------------------------------------');
%--------------------------------------------------------------------
%AUGMENTED PID OPEN-LOOP SYSTEM MATRICES ----------------------------
%Open-loop autonomous system matrix (Nx+2*Ny,Nx+2*Ny)
A_pid= [A              zeros(Nx,Ny)	zeros(Nx,Ny)
     Cy             zeros(Ny)		zeros(Ny)
	 -PhiD*Cy*A 	zeros(Ny)		-PhiD*eye(Ny)];

%Open-loop system input gain matrix (Nx+2*Ny,Nu)
Bu_pid=[Bu
	 zeros(Ny,Nu)
	 -PhiD*Cy*Bu];
	 
%Open-loop system exogenous input gain matrix (Nx+2*Ny,Nw)
Bw_pid=[Bw
	 zeros(Ny,Nw)
	 -PhiD*Cy*Bw];

%Controlled Output-State matrix relationship (Ny+2*Ny,Ny+2*Ny)
Cy_pid=[Cy             zeros(Ny)       zeros(Ny)
	 zeros(Ny,Nx)   eye(Ny)			zeros(Ny)
	 zeros(Ny,Nx)	zeros(Ny)		eye(Ny)];
%-------------------------------------------------------------------- 
%PID RESULTS ---------------------------------------------------------
%PID tuning by SOF stabilization
disp('----- start of PID calculation --------------------------------------');
timerVal = tic;
[KPpid,KIpid,KDpid,TIpid,TDpid,alppid,P_pid]=...
    PID_SOFS_L2N(Poly,gL2,Tk,PhiD,KPlb,KPub,TIlb,TIub,TDlb,TDub,'PID');
elapsedTime_pid = toc(timerVal);

timerVal = tic;
[KPpid_ILMI,KIpid_ILMI,KDpid_ILMI,TIpid_ILMI,TDpid_ILMI,alppid_ILMI,P_pid_ILMI,iterpid_ILMI]=...
    PID_SOFS_L2N_ILMI(Poly,gL2,Tk,PhiD,KPlb,KPub,TIlb,TIub,TDlb,TDub,'PID');
elapsedTime_pid_ILMI = toc(timerVal);
disp('----- end of PID calculation ----------------------------------------');
%--------------------------------------------------------------------
diary off;