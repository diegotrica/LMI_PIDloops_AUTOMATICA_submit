function [KP,KI,KD,TI,TD,alp,P_]=...
	PID_SOFS(Poly,Tk,PhiD,KPlb,KPub,TIlb,TIub,TDlb,TDub,law)
%	Semidefinite programming (SDP) problem to obtain feedback PID 
%	controller gain which maximizes system decay rate.
%
%	The PID tuning is cast as a static output feedback (SOF) stabili-
%	zation problem. SDP seeks to maximize system decay rate solving 
%	the generalized eigenvalue problem (GEVP) using MATLAB(C) LMI 
%	Toolbox.
%
% 	[KP,KI,KD,TI,TD,alp,P_]=...
%		PID_SOFS(Poly,Tk,PhiD,KPlb,KPub,TIlb,TIub,TDlb,TDub,law)
%
%---- Outputs -------------------------------------------------------
%   KP: 	Proportional action gain matrix
%   KI: 	Integral action gain matrix
%   KD: 	Derivative action gain matrix
%   TI: 	Integral time matrix
%   TD: 	Derivative time matrix
%   alp: 	Decay rate
%   P_: 	Augmented system Lyapunov matrix
%--------------------------------------------------------------------
%---- Inputs --------------------------------------------------------
%	Poly: 	Structure storing LDI (A,Bu,Bw,Cy) polytope verteces
%		.Nx:		Number of states
%		.Ny:		Number of outputs
%		.Nw:		Number of exogenous inputs
%		.Nl:		Number of polytope vertices
%		.A{1:Nl}:	Autonomous system matrix (Nx,Nx)
%   	.Bu{1:Nl}: 	Maniputaled inputs gain matrix (Nx,Ny)
%   	.Cy{1:Nl}: 	Output-state matrix relationship (Ny,Nx)
%	PhiD:	Derivative kick preventor constant matrix (diagonal)
%   Tk: 	Control action structure matrix
%   KPlb: 	Proportional gain matrix lower bound
%   KPub: 	Proportional gain matrix upper bound
%   TIlb: 	Integral time matrix lower bound
%   TIub: 	Integral time matrix upper bound
%   TDlb: 	Derivative time matrix lower bound
%   TDub: 	Derivative time matrix upper bound
%   law: 	PID law (default values: 'P', 'PI', 'PID')
%--------------------------------------------------------------------

%	Developed for the paper:
%	Trica, D. J. 2024. Robust tuning of multiple single-paired PID 
%	loops by non-iterative LMI approach.
%
%   Diego Trica (diegotrica@gmail.com)
%--------------------------------------------------------------------
warning('off','MATLAB:rankDeficientMatrix');

%Initialization parameters ------------------------------------------
rtol=1e-6;		%SDP relative tolerance
Nx=Poly.Nx; 	%State variables size
Ny=Poly.Ny;		%Output variables size
Nl=Poly.Nl;		%Number of polytope LTI system verteces
Tk=sign(Tk); 	%Control action structure normalization
Uk=abs(Tk); 	%Diagonalization matrix
Bu=Poly.Bu;		%System maniputaled inputs polytope
Cy=Poly.Cy;		%Controlled output-state relationship polytope
Tl=cell(Nl,1);	%Left-multiplier matrix of equivalence relation
Tr=cell(Nl,1);	%Right-multiplier matrix of equivalence relation
R=cell(Nl,1);	%Weighting matrix of equivalence relation
for i=1:Nl
	R{i}=Bu{i}'*Bu{i};
	Tl{i}=[Bu{i} null(Bu{i}')];
	Tr{i}=[Cy{i};null(Cy{i})']';
end

%Storing into Poly
Poly.Tl=Tl; Poly.Tr=Tr; Poly.R=R;

%K_ variable structure
switch law
    case 'P'
fixed_K_=find(Tk~=0);
sK_=zeros(Ny,Ny); sK_(fixed_K_)=1:length(fixed_K_);

    case 'PI'
fixed_K_=find([Tk,Tk]~=0);
sK_=zeros(Ny,2*Ny); sK_(fixed_K_)=1:length(fixed_K_);

    case 'PID'
fixed_K_=find([Tk,Tk,Tk]~=0);
sK_=zeros(Ny,3*Ny); sK_(fixed_K_)=1:length(fixed_K_);
end   
%--------------------------------------------------------------------
%Set LMI constraints ------------------------------------------------
switch law
    case 'P'
LMISYS=LMI_P(Poly,sK_,Tk,KPlb,KPub);	
    
	case 'PI'
LMISYS=LMI_PI(Poly,sK_,Tk,KPlb,KPub,TIlb,TIub);

    case 'PID'
LMISYS=LMI_PID(Poly,sK_,Tk,PhiD,KPlb,KPub,TIlb,TIub,TDlb,TDub);
end 

%Solve SDP problem
LMIOPT=[rtol,1000,-1,0,0]; %LMI solver options
[alp,z]=gevp(LMISYS,Nl,LMIOPT,[],[],-Inf);
%--------------------------------------------------------------------
%Output solution ----------------------------------------------------
if isempty(alp)
	warning('Unfeasible problem. Check or modify constraints');
	KP=[];
	KI=[];
	KD=[];
	TI=[];
	TD=[];
	alp=[];
	P_=[];
else
%Augmented control gain matrix
K_=dec2mat(LMISYS,z,1);

%Augmented proportional gain matrix
if Ny<Nx
	Kprp=dec2mat(LMISYS,z,6);
else
	Kprp=dec2mat(LMISYS,z,3);
end

%KP, KI, KD and P_
switch law
    case 'P'
KP=K_(:,1:Ny);
KI=[]; TI=[];
KD=[]; TD=[];
Pprp=cell(Nl,1);
for i=1:Nl
	Pprp{i}=0.5*(Tl{i}*Kprp*Tr{i}' + (Tl{i}*Kprp*Tr{i}')');
end
P_=Pprp;

    case 'PI'
KP=K_(:,1:Ny);
KI=K_(:,1+Ny:2*Ny); TI=(KP*Tk')/(KI*Tk')*Uk;
KD=[]; TD=[];
Pprp=cell(Nl,1); P_=cell(Nl,1);
for i=1:Nl
	Pprp{i}=0.5*(Tl{i}*Kprp*Tr{i}' + (Tl{i}*Kprp*Tr{i}')');
	P_{i}=[Pprp{i}        Bu{i}*Tk
		   Tk'*Bu{i}'     eye(Ny)];
end

    case 'PID'
KP=K_(:,1:Ny);
KI=K_(:,1+Ny:2*Ny); TI=(KP*Tk')/(KI*Tk')*Uk;
KD=K_(:,1+2*Ny:3*Ny); TD=(KD*Tk')/(KP*Tk')*Uk;
Pprp=cell(Nl,1); SD=cell(Nl,1); P_=cell(Nl,1);
for i=1:Nl
	SD{i}=Bu{i}*KD+Cy{i}'*(KD'*KD);
	Pprp{i}=0.5*(Tl{i}*Kprp*Tr{i}' + (Tl{i}*Kprp*Tr{i}')' ...
					+ (Cy{i}'*PhiD*SD{i}') + (Cy{i}'*PhiD*SD{i}')');
	P_{i}=[Pprp{i}            	 Bu{i}*Tk   (Bu{i}+Cy{i}'*KD')*Tk
		   Tk'*Bu{i}'         	 eye(Ny)     zeros(Ny)
		   Tk'*(Bu{i}'+KD*Cy{i}) zeros(Ny)	 inv(PhiD)];
end
end
end
%--------------------------------------------------------------------
warning('on','MATLAB:rankDeficientMatrix');
end

function LMISYS=LMI_P(Poly,sK_,Tk,KPlb,KPub)
%Parameters
%Uk=abs(Tk); 	%Diagonalization matrix
Nx=Poly.Nx; 	%State variables size
Ny=Poly.Ny;		%Output variables size
Nl=Poly.Nl;		%Number of polytope LTI system verteces
A=Poly.A;		%Autonomous system matrix polytope
%Bu=Poly.Bu;	%System maniputaled inputs polytope
%Cy=Poly.Cy;	%Controlled output-state relationship polytope
Tl=Poly.Tl;     %Left-multiplier matrix of equivalence relation
Tr=Poly.Tr;     %Right-multiplier matrix of equivalence relation
%R=Poly.R;      %Weighting matrix of equivalence relation

%Initialize LMI system setting
setlmis([]);
	
%LMI system variables
lmivar(3,sK_);  					%Augmented control gain
[KP,~,sKP]=lmivar(3,sK_(:,1:Ny));	%Proportional gain

%Structure for augmented propotional gain
if Ny<Nx
	[~,~,sXP]=lmivar(2,[Ny,Nx-Ny]);
	[~,~,sYP]=lmivar(2,[Nx-Ny,Ny]);
	[~,~,sZP]=lmivar(2,[Nx-Ny,Nx-Ny]);
else
	sXP=[]; sYP=[]; sZP=[];
end
	
Kprp=lmivar(3,[[sKP,sXP];[sYP,sZP]]);%Augmented prop. gain

%LMI system constraints
neq=0;

for i=1:nl
%P_{i} > 0
	neq=neq+1;
	lmiterm([-neq 1 1 Kprp],0.5*Tl{i},Tr{i}','s');
end  	
	
%|KP| > |KPlb|
neq=neq+1;
lmiterm([-neq 1 1 KP],1,Tk','s');
lmiterm([neq 1 1 0],KPlb*Tk'+(KPlb*Tk')');
	
%|KP| < |KPub|
neq=neq+1;
lmiterm([neq 1 1 KP],1,Tk','s');
lmiterm([-neq 1 1 0],KPub*Tk'+(KPub*Tk')');

for i=1:nl
%P_{i} > 0
	neq=neq+1;
	lmiterm([-neq 1 1 Kprp],0.5*Tl{i},Tr{i}','s');

%(LFC) F{i} < 2*alp*P_{i}
	neq=neq+1;
	lmiterm([neq 1 1 Kprp],0.5*A{i}'*Tl{i},Tr{i}','s');
	lmiterm([neq 1 1 Kprp],Tl{i},0.5*Tr{i}'*A{i},'s');
	lmiterm([-neq 1 1 Kprp],Tl{i},Tr{i}','s');
end    
	
%End and store LMI system
LMISYS=getlmis;
end

function LMISYS=LMI_PI(Poly,sK_,Tk,KPlb,KPub,TIlb,TIub)
%Parameters
Uk=abs(Tk); 	%Diagonalization matrix
Nx=Poly.Nx; 	%State variables size
Ny=Poly.Ny;		%Output variables size
Nl=Poly.Nl;		%Number of polytope LTI system verteces
A=Poly.A;		%Autonomous system matrix polytope
Bu=Poly.Bu;     %System maniputaled inputs polytope
Cy=Poly.Cy;     %Controlled output-state relationship polytope
Tl=Poly.Tl;     %Left-multiplier matrix of equivalence relation
Tr=Poly.Tr;     %Right-multiplier matrix of equivalence relation
R=Poly.R;       %Weighting matrix of equivalence relation

%Initialize LMI system setting
setlmis([]);
	
%LMI system variables
lmivar(3,sK_);  					%Augmented control gain
[KP,~,sKP]=lmivar(3,sK_(:,1:Ny));	%Proportional gain

%Structure for augmented propotional gain
if Ny<Nx
	[~,~,sXP]=lmivar(2,[Ny,Nx-Ny]);
	[~,~,sYP]=lmivar(2,[Nx-Ny,Ny]);
	[~,~,sZP]=lmivar(2,[Nx-Ny,Nx-Ny]);
else
	sXP=[]; sYP=[]; sZP=[];
end
	
Kprp=lmivar(3,[[sKP,sXP];[sYP,sZP]]);%Augmented prop. gain
KI=lmivar(3,sK_(:,1+Ny:2*Ny));		 %Integral gain

%LMI system constraints
neq=0;

for i=1:Nl
%P_{i} > 0
	neq=neq+1;
	lmiterm([-neq 1 1 Kprp],0.5*Tl{i},Tr{i}','s');
	lmiterm([-neq 2 1 0],Tk'*Bu{i}');
	lmiterm([-neq 2 2 0],eye(Ny));
end

%|KP| > |KPlb|
neq=neq+1;
lmiterm([-neq 1 1 KP],1,Tk','s');
lmiterm([neq 1 1 0],KPlb*Tk'+(KPlb*Tk')');
	
%|KP| < |KPub|
neq=neq+1;
lmiterm([neq 1 1 KP],1,Tk','s');
lmiterm([-neq 1 1 0],KPub*Tk'+(KPub*Tk')');

%|KI| > |TIub\KP|
neq=neq+1;
lmiterm([-neq 1 1 KI],1,Tk','s');
lmiterm([neq 1 1 KP],inv(TIub*Uk'),Tk','s');
	
%|KI| < |TIlb\KP|
neq=neq+1;
lmiterm([neq 1 1 KI],1,Tk','s');
lmiterm([-neq 1 1 KP],inv(TIlb*Uk'),Tk','s');

for i=1:Nl
%(LFC) F{i} < 2*alp*P_{i}
	neq=neq+1;
	lmiterm([neq 1 1 Kprp],0.5*A{i}'*Tl{i},Tr{i}','s');
	lmiterm([neq 1 1 Kprp],Tl{i},0.5*Tr{i}'*A{i},'s');
	lmiterm([neq 1 1 KI],Bu{i},Cy{i},'s');
	lmiterm([neq 2 1 0],Tk'*Bu{i}'*A{i});
	lmiterm([neq 2 1 KP],-Tk'*(2*R{i}),Cy{i});
	lmiterm([neq 2 1 KI],Tk',Cy{i});
	lmiterm([neq 2 2 0],-Tk'*(2*R{i})*Tk);
	lmiterm([-neq 1 1 Kprp],Tl{i},Tr{i}','s');
	lmiterm([-neq 2 1 0],2*Tk'*Bu{i}');
	lmiterm([-neq 2 2 0],2*eye(Ny));
end

%End and store LMI system
LMISYS=getlmis;
end

function LMISYS=...
	LMI_PID(Poly,sK_,Tk,PhiD,KPlb,KPub,TIlb,TIub,TDlb,TDub)
%Parameters
Uk=abs(Tk); 	%Diagonalization matrix
Nx=Poly.Nx; 	%State variables size
Ny=Poly.Ny;		%Output variables size
Nl=Poly.Nl;		%Number of polytope LTI system verteces
A=Poly.A;		%Autonomous system matrix polytope
Bu=Poly.Bu;     %System maniputaled inputs polytope
Cy=Poly.Cy;     %Controlled output-state relationship polytope
Tl=Poly.Tl;     %Left-multiplier matrix of equivalence relation
Tr=Poly.Tr;     %Right-multiplier matrix of equivalence relation
R=Poly.R;       %Weighting matrix of equivalence relation

%Initialize LMI system setting
setlmis([]);
	
%LMI system variables
lmivar(3,sK_);  					%Augmented control gain
[KP,~,sKP]=lmivar(3,sK_(:,1:Ny));	%Proportional gain

if Ny<Nx
	[~,~,sXP]=lmivar(2,[Ny,Nx-Ny]);
	[~,~,sYP]=lmivar(2,[Nx-Ny,Ny]);
	[~,~,sZP]=lmivar(2,[Nx-Ny,Nx-Ny]);
else
	sXP=[]; sYP=[]; sZP=[];
end
	
Kprp=lmivar(3,[[sKP,sXP];[sYP,sZP]]);%Augmented prop. gain
KI=lmivar(3,sK_(:,1+Ny:2*Ny));		 %Integral gain
KD=lmivar(3,sK_(:,1+2*Ny:3*Ny));	 %Derivative gain

%LMI system constraints
neq=0;

for i=1:Nl
%P_{i} > 0
	neq=neq+1;
	lmiterm([-neq 1 1 Kprp],0.5*Tl{i},Tr{i}','s');
	lmiterm([-neq 1 1 KD],0.5*Bu{i},PhiD*Cy{i},'s');
	lmiterm([-neq 2 1 0],Tk'*Bu{i}');
	lmiterm([-neq 3 1 0],Tk'*Bu{i}');
	lmiterm([-neq 3 1 KD],Tk',Cy{i});
	lmiterm([-neq 2 2 0],eye(Ny));
	lmiterm([-neq 3 3 0],eye(Ny)/PhiD);
end
	
%|KP| > |KPlb|
neq=neq+1;
lmiterm([-neq 1 1 KP],1,Tk','s');
lmiterm([neq 1 1 0],KPlb*Tk'+(KPlb*Tk')');
	
%|KP| < |KPub|
neq=neq+1;
lmiterm([neq 1 1 KP],1,Tk','s');
lmiterm([-neq 1 1 0],KPub*Tk'+(KPub*Tk')');

%|KI| > |TIub\KP|
neq=neq+1;
lmiterm([-neq 1 1 KI],1,Tk','s');
lmiterm([neq 1 1 KP],inv(TIub*Uk'),Tk','s');
	
%|KI| < |TIlb\KP|
neq=neq+1;
lmiterm([neq 1 1 KI],1,Tk','s');
lmiterm([-neq 1 1 KP],inv(TIlb*Uk'),Tk','s');

%|KD| > |TDlb*KP|
neq=neq+1;
lmiterm([-neq 1 1 KD],1,Tk','s');
lmiterm([neq 1 1 KP],(TDlb*Uk'),Tk','s');
	
%|KD| < |TDub*KP|
neq=neq+1;
lmiterm([neq 1 1 KD],1,Tk','s');
lmiterm([-neq 1 1 KP],(TDub*Uk'),Tk','s');

for i=1:Nl
%(LFC) F{i} < 2*alp*P_{i}
	neq=neq+1;
	lmiterm([neq 1 1 Kprp],0.5*A{i}'*Tl{i},Tr{i}','s');
	lmiterm([neq 1 1 Kprp],Tl{i},0.5*Tr{i}'*A{i},'s');
	lmiterm([neq 1 1 KD],0.5*A{i}'*Bu{i},PhiD*Cy{i},'s');
	lmiterm([neq 1 1 KD],-0.5*Bu{i},PhiD*Cy{i}*A{i},'s');
	lmiterm([neq 1 1 KI],Bu{i},Cy{i},'s');
	lmiterm([neq 2 1 0],Tk'*Bu{i}'*A{i});
	lmiterm([neq 2 1 KP],-Tk'*(2*R{i}),Cy{i});
	lmiterm([neq 2 1 KI],Tk',Cy{i});
	lmiterm([neq 3 1 0],Tk'*Bu{i}'*A{i});
	lmiterm([neq 3 1 0],-PhiD*Tk'*Bu{i}');
	lmiterm([neq 3 1 KP],-Tk'*(2*R{i}),Cy{i});
	lmiterm([neq 3 1 KD],-PhiD*Tk',Cy{i});
	lmiterm([neq 2 2 0],-Tk'*(2*R{i})*Tk);
	lmiterm([neq 3 2 0],-Tk'*(2*R{i})*Tk);
	lmiterm([neq 3 3 0],-Tk'*2*(R{i}+eye(Ny))*Tk);
	lmiterm([-neq 1 1 Kprp],Tl{i},Tr{i}','s');
	lmiterm([-neq 1 1 KD],Bu{i},PhiD*Cy{i},'s');
	lmiterm([-neq 2 1 0],2*Tk'*Bu{i}');
	lmiterm([-neq 3 1 0],2*Tk'*Bu{i}');
	lmiterm([-neq 3 1 KD],2*Tk',Cy{i});
	lmiterm([-neq 2 2 0],2*eye(Ny));
	lmiterm([-neq 3 3 0],2*eye(Ny)/PhiD);
end

%End and store LMI system
LMISYS=getlmis;
end