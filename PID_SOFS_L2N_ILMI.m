function [KP,KI,KD,TI,TD,alp,P_,iter]=...
	PID_SOFS_L2N_ILMI(Poly,gL2,Tk,PhiD,KPlb,KPub,TIlb,TIub,TDlb,TDub,law)
%	Semidefinite programming (SDP) problem to obtain feedback PID 
%	controller gain which maximizes system decay rate with upper
%	bound on Hinf norm by L2 function norms. The iterative LMI 
%	approach is given by Wang et. al. (2008) PID Control for Multi-
%	variable Processes, Section 7.5.
%
%	The PID tuning is cast as a static output feedback (SOF) stabili-
%	zation problem. SDP seeks to maximize system decay rate solving 
%	the generalized eigenvalue problem (GEVP) using MATLAB(C) LMI 
%	Toolbox.
%
% 	[KP,KI,KD,TI,TD,alp,P_,iter]=...
%		PID_SOFS_L2N_ILMI(Poly,gL2,PhiD,KPlb,KPub,TIlb,TIub,TDlb,TDub,law)
%
%---- Outputs -------------------------------------------------------
%   KP: 	Proportional action gain matrix
%   KI: 	Integral action gain matrix
%   KD: 	Derivative action gain matrix
%   TI: 	Integral time matrix
%   TD: 	Derivative time matrix
%   alp: 	Decay rate
%   P_: 	Augmented system Lyapunov matrix
%	iter:	Number of LMI iterations
%--------------------------------------------------------------------
%---- Inputs --------------------------------------------------------
%	Poly: 	Structure storing LDI (A,Bu,Bw,Cy) polytope verteces
%		.Nx:		Number of states
%		.Ny:		Number of outputs
%		.Nw:		Number of exogenous inputs
%		.Nl:		Number of polytope vertices
%		.A{1:Nl}:	Autonomous system matrix (Nx,Nx)
%   	.Bu{1:Nl}: 	Maniputaled inputs gain matrix (Nx,Ny)
%   	.Bw{1:Nl}: 	Exogenous inputs gain matrix (Nx,Nw)
%   	.Cy{1:Nl}: 	Output-state matrix relationship (Ny,Nx)
%	gL2:	Upper bound on L2 norm
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
Nu=Ny;
Nl=Poly.Nl;		%Number of polytope LTI system verteces
Tk=sign(Tk); 	%Control action structure normalization
Uk=abs(Tk); 	%Diagonalization matrix

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
%Solve SDP problem by iterative LMI approach ------------------------
switch law
    case 'P'
A_=Poly.A0;
Bu_=Poly.Bu0;
    
	case 'PI'
A_ = [Poly.A0 				zeros(Nx,Ny)
	  Poly.Cy0  			zeros(Ny,Ny)];

Bu_= [Poly.Bu0
	  zeros(Ny,Nu)];

    case 'PID'
A_ = [Poly.A0               zeros(Nx,Ny)	zeros(Nx,Ny)
      Poly.Cy0              zeros(Ny)		zeros(Ny)
	 -PhiD*Poly.Cy0*Poly.A0 zeros(Ny)		-PhiD*eye(Ny)];

Bu_= [Poly.Bu0
	 zeros(Ny,Nu)
	 -PhiD*Poly.Cy0*Poly.Bu0];
	 
end
Q0_=eye(size(A_,1));
[P_,~,~]=care(A_,Bu_,Q0_);	%Step 1

iter=0; maxiter=100; alp=[]; zOP1=[];
LMIOPT=[rtol,500,-1,0,0]; %LMI solver options
while iter<maxiter
iter=iter+1;
X_=P_;
disp(strjoin({'----- Iteration #',num2str(iter),' --------------------------------------------------'},''));
disp('----- solving gevp for SOFHinf --------------------------------------');
LMISYS=LMIOP1(Poly,gL2,sK_,Tk,PhiD,KPlb,KPub,TIlb,TIub,TDlb,TDub,X_,law);  		%Step 2
[alp,zOP1]=gevp(LMISYS,Nl,LMIOPT,alp,zOP1,-Inf);
if isempty(alp)
	warning('Unfeasible problem. Check or modify constraints');
	break;
end
K_=dec2mat(LMISYS,zOP1,1);
P_=dec2mat(LMISYS,zOP1,2);
if alp<0	%Step 3
	break;
end
disp('----- alp > 0. Minimizing tr(P) for next iteration ------------------');
LMISYS=LMIOP2(Poly,gL2,sK_,Tk,PhiD,KPlb,KPub,TIlb,TIub,TDlb,TDub,X_,alp,law);	%Step 4
c=zeros(decnbr(LMISYS),1); idtrP_=diag(decinfo(LMISYS,2)); c(idtrP_)=1;
zOP2=zOP1(1:size(K_(K_~=0),1)+size(P_(tril(P_)~=0),1));
[trP,zOP2]=mincx(LMISYS,c,[LMIOPT(1:end-1),1],zOP2,0);
%{
Remark: The step 4 SOFHinf algorithm was modified as below. 
P_ is updated by P_* only if the min tr(P) problem is feasible.
It was not clear which value of alp+dalp must be given to achieve
feasilibty of min tr(P). Issues about feasbility of min tr(P) problem
is mentioned in the paper where originally the algorithm was
developed (Cao et al, 1998).
%}
if ~isempty(zOP2)
	K_=dec2mat(LMISYS,zOP2,1);
	P_=dec2mat(LMISYS,zOP2,2);
	zOP1(1:size(K_(K_~=0),1)+size(P_(tril(P_)~=0),1))=zOP2;
	disp(strjoin({'----- P_ updated by P_*. tr(P) = ',num2str(trP)},''));
else
	disp('----- P_ not updated by P_* because tr(P) problem is unfeasible -----');
end
if norm(X_*Bu_-P_*Bu_)<rtol	%Step 5
	warning('It cannot be decided by this algorithm whether SOFL2 problem is solvable');
    break;
end
end
if iter==maxiter
    warning('Maximum number of LMI iterations achieved');
end
%--------------------------------------------------------------------
%Output solution ----------------------------------------------------
if isempty(alp)
	KP=[];
	KI=[];
	KD=[];
	TI=[];
	TD=[];
	alp=[];
	P_=[];
else
%KP, KI, KD and P_
switch law
    case 'P'
KP=K_(:,1:Ny);
KI=[]; TI=[];
KD=[]; TD=[];

    case 'PI'
KP=K_(:,1:Ny);
KI=K_(:,1+Ny:2*Ny); TI=(KP*Tk')/(KI*Tk')*Uk;
KD=[]; TD=[];

    case 'PID'
KP=K_(:,1:Ny);
KI=K_(:,1+Ny:2*Ny); TI=(KP*Tk')/(KI*Tk')*Uk;
KD=K_(:,1+2*Ny:3*Ny); TD=(KD*Tk')/(KP*Tk')*Uk;
end
end
%--------------------------------------------------------------------
warning('on','MATLAB:rankDeficientMatrix');
end

function LMISYS=...
	LMIOP1(Poly,gL2,sK_,Tk,PhiD,KPlb,KPub,TIlb,TIub,TDlb,TDub,X_,law)
%Parameters
Uk=abs(Tk); 	%Diagonalization matrix
Nx=Poly.Nx; 	%State variables size
Ny=Poly.Ny;		%Output variables size
Nw=Poly.Nw;		%Exogenous input variables size
Nl=Poly.Nl;		%Number of polytope LTI system verteces

A_=cell(Nl,1); Bu_=cell(Nl,1); Cy_=cell(Nl,1);
for i=1:Nl
switch law
    case 'P'
A_{i}=Poly.A{i};
Bu_=Poly.Bu{i};
Bw_=Poly.Bw{i};
Cy_{i}=Poly.Cy{i};
    
	case 'PI'
A_{i} = [Poly.A{i} 			zeros(Nx,Ny)
	  Poly.Cy{i}  			zeros(Ny,Ny)];

Bu_{i}= [Poly.Bu{i}
	  zeros(Ny)];
	  
Bw_{i}= [Poly.Bw{i}
	  zeros(Ny,Nw)];	  
	 
Cy_{i} = [Poly.Cy{i}     	zeros(Ny,Ny)
	   zeros(Ny,Nx) 		eye(Ny,Ny)];	

    case 'PID'
A_{i} = [Poly.A{i}              zeros(Nx,Ny)	zeros(Nx,Ny)
      Poly.Cy{i}              	zeros(Ny)		zeros(Ny)
	 -PhiD*Poly.Cy{i}*Poly.A{i} zeros(Ny)		-PhiD*eye(Ny)];

Bu_{i}= [Poly.Bu{i}
	 zeros(Ny)
	 -PhiD*Poly.Cy{i}*Poly.Bu{i}];
	 
Bw_{i}= [Poly.Bw{i}
	 zeros(Ny,Nw)
	 -PhiD*Poly.Cy{i}*Poly.Bw{i}];	 

Cy_{i}= [Poly.Cy{i}        	zeros(Ny)       zeros(Ny)
	  zeros(Ny,Nx)   		eye(Ny)			zeros(Ny)
	  zeros(Ny,Nx)			zeros(Ny)		eye(Ny)];
	 
end
end

%Initialize LMI system setting
setlmis([]);
	
%LMI system variables
K_=lmivar(3,sK_);  					%Augmented control gain

switch law
    case 'P'                        %Lyapunov matrix for P control
P_=lmivar(1,[Nx 1]);
Sgm_=lmivar(1,[Nx 1]);
    
    case 'PI'                       %Lyapunov matrix for PI control
P_=lmivar(1,[Nx+Ny 1]);
Sgm_=lmivar(1,[Nx+Ny 1]);

    case 'PID'                      %Lyapunov matrix for PID control
P_=lmivar(1,[Nx+Ny+Ny 1]);
Sgm_=lmivar(1,[Nx+Ny+Ny 1]);
end

KP=lmivar(3,sK_(:,1:Ny));			%Proportional gain

%LMI system constraints
neq=0;

%P_ > 0
neq=neq+1;
lmiterm([-neq 1 1 P_],1,1);

%Sgm_ < 0
neq=neq+1;
lmiterm([neq 1 1 Sgm_],1,1);

%|KP| > |KPlb|
neq=neq+1;
lmiterm([-neq 1 1 KP],1,Tk','s');
lmiterm([neq 1 1 0],KPlb*Tk'+(KPlb*Tk')');

%|KP| < |KPub|
neq=neq+1;
lmiterm([neq 1 1 KP],1,Tk','s');
lmiterm([-neq 1 1 0],KPub*Tk'+(KPub*Tk')');

switch law
	case 'PI'
KI=lmivar(3,sK_(:,1+Ny:2*Ny));		%Integral gain

%|KI| > |TIub\KP|
neq=neq+1;
lmiterm([-neq 1 1 KI],1,Tk','s');
lmiterm([neq 1 1 KP],inv(TIub*Uk'),Tk','s');
	
%|KI| < |TIlb\KP|
neq=neq+1;
lmiterm([neq 1 1 KI],1,Tk','s');
lmiterm([-neq 1 1 KP],inv(TIlb*Uk'),Tk','s');

    case 'PID'
KI=lmivar(3,sK_(:,1+Ny:2*Ny));		%Integral gain
KD=lmivar(3,sK_(:,1+2*Ny:3*Ny));	%Derivative gain

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
end

for i=1:Nl
%{
H{i} < 0
H{i}(P_,K_,alp) = 
	A_'*P_ + P_*A_ + Cy'*Cy + P_*((Bw_/gL2^2)*Bw_')*P_ - P_*Bu_*K_*Cy_ - Cy_'*K_'*Bu_'*P_ < 2*alp*P_
	A_'*P_ + P_*A_ + Cy'*Cy + P_*((Bw_/gL2^2)*Bw_')*P_ - P_*(Bu_*Bu_')*P_ + (Bu_'*P_-K_*Cy_)'*(Bu_'*P_-K_*Cy_) - (K_*Cy_)'*(K_*Cy_) < 2*alp*P_
dropping - (K_*Cy_)'*(K_*Cy_) yields to:
	A_'*P_ + P_*A_ + Cy'*Cy - Psi_ + P_*((Bw_/gL2^2)*Bw_')*P_ + (Bu_'*P_-K_*Cy_)'*(Bu_'*P_-K_*Cy_) < 2*alp*P_
By schur complements
[A_'*P_ + P_*A_ - Psi_	Cy'		P_*Bw_		(Bu_'*P_-K_*Cy_)'] < [2*alp*P_	0	0   0  ]
[Cy 					-I		0			0				 ]   [0     	eps	0	0  ]
[Bw_'*P_ 				0		-gL2^2		0 				 ]   [0        	0	eps 0  ]
[(Bu_'*P_-K_*Cy_)		0		0			-I 				 ]   [0        	0   0	eps]
Psi_ = X_'*(Bu_*Bu_')*P_ + P_*(Bu_*Bu_')*X_ - X_'*(Bu_*Bu_')*X_
%}
	neq=neq+1;
	lmiterm([neq 1 1 Sgm_],1,1);
	lmiterm([neq 1 1 P_],A_{i}',1,'s');
	lmiterm([neq 1 1 P_],-X_'*(Bu_{i}*Bu_{i}'),1,'s');
	lmiterm([neq 1 1 0],X_'*(Bu_{i}*Bu_{i}')*X_);	
	lmiterm([neq 2 1 0],Cy_{i});
    lmiterm([neq 3 1 P_],Bw_{i}',1);
	lmiterm([neq 4 1 P_],Bu_{i}',1);
	lmiterm([neq 4 1 K_],1,-Cy_{i});
	lmiterm([neq 2 2 0],-1);
    lmiterm([neq 3 3 0],-gL2^2);
	lmiterm([neq 4 4 0],-1);

end

for i=1:Nl
%(LFC) A_'*P_ + P_*A_ - Psi_ - Sgm_ < 2*alp*P_	
	neq=neq+1;
	lmiterm([neq 1 1 P_],A_{i}',1,'s');
	lmiterm([neq 1 1 P_],-X_'*(Bu_{i}*Bu_{i}'),1,'s');
	lmiterm([neq 1 1 0],X_'*(Bu_{i}*Bu_{i}')*X_);
	lmiterm([neq 1 1 Sgm_],-1,1);
	lmiterm([-neq 1 1 P_],1,1,'s');
end

%End and store LMI system
LMISYS=getlmis;
end

function LMISYS=...
	LMIOP2(Poly,gL2,sK_,Tk,PhiD,KPlb,KPub,TIlb,TIub,TDlb,TDub,X_,alp,law)
%Parameters
Uk=abs(Tk); 	%Diagonalization matrix
Nx=Poly.Nx; 	%State variables size
Ny=Poly.Ny;		%Output variables size
Nw=Poly.Nw;		%Exogenous input variables size
Nl=Poly.Nl;		%Number of polytope LTI system verteces

A_=cell(Nl,1); Bu_=cell(Nl,1); Cy_=cell(Nl,1);
for i=1:Nl
switch law
    case 'P'
A_{i}=Poly.A{i};
Bu_=Poly.Bu{i};
Bw_=Poly.Bw{i};
Cy_{i}=Poly.Cy{i};
    
	case 'PI'
A_{i} = [Poly.A{i} 			zeros(Nx,Ny)
	  Poly.Cy{i}  			zeros(Ny,Ny)];

Bu_{i}= [Poly.Bu{i}
	  zeros(Ny)];

Bw_{i}= [Poly.Bw{i}
	  zeros(Ny,Nw)];
	 
Cy_{i} = [Poly.Cy{i}     	zeros(Ny,Ny)
	   zeros(Ny,Nx) 		eye(Ny,Ny)];	

    case 'PID'
A_{i} = [Poly.A{i}              zeros(Nx,Ny)	zeros(Nx,Ny)
      Poly.Cy{i}              	zeros(Ny)		zeros(Ny)
	 -PhiD*Poly.Cy{i}*Poly.A{i} zeros(Ny)		-PhiD*eye(Ny)];

Bu_{i}= [Poly.Bu{i}
	 zeros(Ny)
	 -PhiD*Poly.Cy{i}*Poly.Bu{i}];

Bw_{i}= [Poly.Bw{i}
	 zeros(Ny,Nw)
	 -PhiD*Poly.Cy{i}*Poly.Bw{i}];	 

Cy_{i}= [Poly.Cy{i}        	zeros(Ny)       zeros(Ny)
	  zeros(Ny,Nx)   		eye(Ny)			zeros(Ny)
	  zeros(Ny,Nx)			zeros(Ny)		eye(Ny)];
	 
end
end

%Initialize LMI system setting
setlmis([]);
	
%LMI system variables
K_=lmivar(3,sK_);  					%Augmented control gain

switch law
    case 'P'                        %Lyapunov matrix for P control
P_=lmivar(1,[Nx 1]);
    
    case 'PI'                       %Lyapunov matrix for PI control
P_=lmivar(1,[Nx+Ny 1]);

    case 'PID'                      %Lyapunov matrix for PID control
P_=lmivar(1,[Nx+Ny+Ny 1]);
end

KP=lmivar(3,sK_(:,1:Ny));			%Proportional gain

%LMI system constraints
neq=0;

%P_ > 0
neq=neq+1;
lmiterm([-neq 1 1 P_],1,1);

%|KP| > |KPlb|
neq=neq+1;
lmiterm([-neq 1 1 KP],1,Tk','s');
lmiterm([neq 1 1 0],KPlb*Tk'+(KPlb*Tk')');

%|KP| < |KPub|
neq=neq+1;
lmiterm([neq 1 1 KP],1,Tk','s');
lmiterm([-neq 1 1 0],KPub*Tk'+(KPub*Tk')');

switch law
	case 'PI'
KI=lmivar(3,sK_(:,1+Ny:2*Ny));		%Integral gain

%|KI| > |TIub\KP|
neq=neq+1;
lmiterm([-neq 1 1 KI],1,Tk','s');
lmiterm([neq 1 1 KP],inv(TIub*Uk'),Tk','s');
	
%|KI| < |TIlb\KP|
neq=neq+1;
lmiterm([neq 1 1 KI],1,Tk','s');
lmiterm([-neq 1 1 KP],inv(TIlb*Uk'),Tk','s');

    case 'PID'
KI=lmivar(3,sK_(:,1+Ny:2*Ny));		%Integral gain
KD=lmivar(3,sK_(:,1+2*Ny:3*Ny));	%Derivative gain

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
end

for i=1:Nl
%{
H{i} < 0
H{i}(P_,K_) = 
	A_'*P_ + P_*A_ + Cy'*Cy + P_*((Bw_/gL2^2)*Bw_')*P_ - P_*Bu_*K_*Cy_ - Cy_'*K_'*Bu_'*P_ < 2*alp*P_
	A_'*P_ + P_*A_ + Cy'*Cy + P_*((Bw_/gL2^2)*Bw_')*P_ - P_*(Bu_*Bu_')*P_ + (Bu_'*P_-K_*Cy_)'*(Bu_'*P_-K_*Cy_) - (K_*Cy_)'*(K_*Cy_) < 2*alp*P_
dropping - (K_*Cy_)'*(K_*Cy_) yields to:
	A_'*P_ + P_*A_ + Cy'*Cy - Psi_ -2*alp*P_ + P_*((Bw_/gL2^2)*Bw_')*P_ + (Bu_'*P_-K_*Cy_)'*(Bu_'*P_-K_*Cy_) < 0
By schur complements
[A_'*P_ + P_*A_ - Psi_ - 2*alp*P_	Cy'		P_*Bw_		(Bu_'*P_-K_*Cy_)'] < 0
[Cy 								-I		0			0				 ]
[Bw_'*P_ 							0		-gL2^2		0 				 ]
[(Bu_'*P_-K_*Cy_)					0		0			-I 				 ]
Psi_ = X_'*(Bu_*Bu_')*P_ + P_*(Bu_*Bu_')*X_ - X_'*(Bu_*Bu_')*X_
%}
	neq=neq+1;
	lmiterm([neq 1 1 P_],A_{i}',1,'s');
	lmiterm([neq 1 1 P_],-alp,1,'s');
	lmiterm([neq 1 1 P_],-X_'*(Bu_{i}*Bu_{i}'),1,'s');
	lmiterm([neq 1 1 0],X_'*(Bu_{i}*Bu_{i}')*X_);
	lmiterm([neq 2 1 0],Cy_{i});
    lmiterm([neq 3 1 P_],Bw_{i}',1);
	lmiterm([neq 4 1 P_],Bu_{i}',1);
	lmiterm([neq 4 1 K_],1,-Cy_{i});
	lmiterm([neq 2 2 0],-1);
    lmiterm([neq 3 3 0],-gL2^2);
	lmiterm([neq 4 4 0],-1);
end

%End and store LMI system
LMISYS=getlmis;
end