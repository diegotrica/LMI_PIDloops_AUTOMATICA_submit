function [sys,x0,str,ts,simStateCompliance]=vandevusse_PIctrl(t,x,u,flag,x0,par)
%vandevusse - Van de Vusse non-isothermal CSTR non-linear model 
%S-Function under PI control
%
%This model describes the dynamic behavior of the Van de Vusse non-
%isothermal CSTR model with the following chemical reaction scheme:
%
%	A --> B --> C
%	2A --> D
%
%Chemical components:
%(A) cyclopentadiene
%(B) cyclopentenol  
%(C) cyclopentanediol
%(D) di-cyclopentadiene
%
%[sys,x0,str,ts,simStateCompliance]=vandevusse_PIctrl(t,x,u,flag,x0,par)
%
%sys is the S-function output, which depends on which flag is 
%being inputted by simulink. x0 are the states initial condi-
%tions. See the S-function Level-1 MATLAB reference page for more
%details about str, ts and simStateCompliance S-Function ouputs.
%
%t is the independent variable time. x are the state variables.
%u are the input variables, combined into one vector of manipula-
%ted and disturbs. flag is the Simulink variable inputted to get a
%specific output from the S-function. x0 are the variable set for 
%state initial condition computation. par are the model parameters.
%
%----- Model parameters ---------------------------------------------
%K1o: 		Arrhenius spec. reaction rate    			[1/min]
%K2o: 		Arrhenius spec. reaction rate    			[1/min]
%K3o: 		Arrhenius spec. reaction rate    			[m^3/kmol min]
%E1R:		Activion energy in terms of temp.			[K]
%E2R:		Activion energy in terms of temp.			[K]
%E3R:		Activion energy in terms of temp.			[K]
%DH1:		Molar reaction heat 			 			[kJ/kmol]
%DH2:		Molar reaction heat 			 			[kJ/kmol]
%DH3:		Molar reaction heat 			 			[kJ/kmol]
%rho:		liquid mass density 			 			[kg/m^3]
%cp: 		liquid specific heat capacity 	 			[kJ/kg K]
%Uf: 		Global heat transfer coeff.      			[kJ/min m^2 K]
%A:			Heat transfer area      		 			[m^2]
%V:			liquid reaction volume 			 			[m^3]
%mc:		cooling medium mass hold-up 			 	[kg]
%cpc:		cooling medium specific heat capacity 	 	[kJ/kg K]
%KP1:		loop #1 proportional gain 					[]
%KP2:		loop #2 proportional gain 					[]
%TI1:		loop #1 integral time 						[min]
%TI2:		loop #2 integral time 						[min]
%th0:		bias for inlet vol. flow rate per liq. vol.	[1/min]
%Qc0:		bias for cooling medium heat removal rate	[kJ/min]
%--------------------------------------------------------------------
%----- States -------------------------------------------------------
%x(1): CA	component A molar concentration 			[kmol/m3]
%x(2): CB	component B molar concentration 			[kmol/m3]
%x(3): T 	reactor temperature							[K]
%x(4): Tc 	cooling medium temperature					[K]
%x(5): Int1 loop #1 integral action						[]
%x(6): Int2 loop #2 integral action						[]
%--------------------------------------------------------------------
%----- Inputs -------------------------------------------------------
%u(1): CBsp	component B molar concentration set-point 	[kmol/m3]
%u(2): Tsp	reactor temperature set-point				[K]
%u(3): CAin	component A inlet molar concentration 		[kmol/m3]
%u(4): Tin	inlet temperature 							[K]
%------ Outputs -----------------------------------------------------
%y(1): CB	component B molar concentration 			[kmol/m3]
%y(2): T	reactor temperature							[K]
%--------------------------------------------------------------------

%----- Parameters:
K1o = par(1);
K2o = par(2);
K3o = par(3);
E1R = par(4);
E2R = par(5);
E3R = par(6);
DH1 = par(7);
DH2 = par(8); 
DH3 = par(9);
rho = par(10);
cp 	= par(11);
Uf 	= par(12);
A 	= par(13);
V 	= par(14);
mc 	= par(15);
cpc = par(16);
KP1 = par(17);
KP2 = par(18);
TI1 = par(19);
TI2 = par(20);
th0 = par(21);
Qc0 = par(22);

switch flag,
%----- S-FUNCTION BUILDING PARAMETERS AND INITIAL CONDITIONS --------
    case 0,        
sys(1)=6;	%Number of continuous states
sys(2)=0;	%Number of discrete states
sys(3)=2;	%Number of outputs
sys(4)=4;	%Number of inputs
sys(5)=0;	%Reserved for root finding. Must be zero.
sys(6)=0;	%Direct Feedthrough S-function parameter
sys(7)=1;	%Number of sample times (at least one sample time is needed)

%----- States initial condition:
x0=x0(:);

%----- Other building parameters:
str=[];     %str is always an empty matrix
ts=[0 0];   %Sample time array [0 0] for continuous system
simStateCompliance='UnknownSimState';
%--------------------------------------------------------------------
%----- STATE DERIVATIVES --------------------------------------------
    case 1,   
%----- States:
CA 		= x(1);
CB 		= x(2);
T 		= x(3);
Tc 		= x(4);
Int1 	= x(5);
Int2 	= x(6);

%----- Inputs:
CBsp 	= u(1);
Tsp 	= u(2);
CAin 	= u(3);
Tin 	= u(4);   

%----- State derivatives:
dInt1dt = CBsp-CB;

dInt2dt = Tsp-T;

%u(1): th	inlet vol. flow rate per liq. reaction vol.	[1/min]
th = th0 + KP1*(dInt1dt + Int1/TI1);

%u(2): Qc	cooling medium heat removal rate			[kJ/min]
Qc = Qc0 + KP2*(dInt2dt + Int2/TI2);

dCAdt = th*(CAin-CA) - K1o*exp(-E1R/T)*CA - K3o*exp(-E3R/T)*CA^2;
            
dCBdt = -th*CB + K1o*exp(-E1R/T)*CA - K2o*exp(-E2R/T)*CB;
            
dTdt = th*(Tin-T) - Uf*A/(rho*cp*V)*(T-Tc) ...
	+ K1o*exp(-E1R/T)*CA*(-DH1)/(rho*cp)... 
	+ K2o*exp(-E2R/T)*CB*(-DH2)/(rho*cp) ...
	+ K3o*exp(-E3R/T)*CA^2*(-DH3)/(rho*cp);

dTcdt = (mc*cpc)\(-Qc + Uf*A*(T-Tc));

%dxdt allocation:
dxdt=[dCAdt; dCBdt; dTdt; dTcdt; dInt1dt; dInt2dt];

%----- State derivatives for sys output flag=1
sys=dxdt;
sys=sys(:);
%--------------------------------------------------------------------
    case 2,
sys=[];
%----- OUPUTS -------------------------------------------------------
    case 3,       
%----- States:
%CA 	= x(1);
CB 		= x(2);
T 		= x(3);
%Tc 	= x(4);

%y allocation
y=[CB; T];

%----- Outputs for sys output flag=3:
sys=y;
sys=sys(:);
%--------------------------------------------------------------------
%----- S-FUNCTION OTHER FLAGS ---------------------------------------
    case 4,
sys=[];
    case 9,
sys=[];
    otherwise   % Unexpected flags
DAStudio.error('Simulink:blocks:unhandledFlag', num2str(flag));
end
%--------------------------------------------------------------------    
end