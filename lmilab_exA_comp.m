clear, close all

Nrlt=0;
load('github\caseA\lmilab_exA_tb_rlt.mat','rlt');

tstart=0; tend=120; tstep=0.1; t=(tstart:tstep:tend)'; %time domain
yf=[0.1;5]; dw=[0.1;5];
InputName={'C_{B,r} [kmol/m^3]','T_r [K]','C_{A,in} [kmol/m^3]','T_{in} [K]'};
OutputName={'C_B [kmol/m^3]','T [K]'};
Ny=size(yf,1); Nw=size(dw,1);
%--------------------------------------------------------------------
%{
Figure 1. Output step response of Van de Vusse non-isothermal
CSTR non-linear model for TW and ILMI PI tuning for gL2=1000
%}
load(rlt(Nrlt+1).data,'u0','x0','par');
y1_nlin=zeros(size(t,1),Ny,Ny);
u=[x0(2:3);u0(3:4)]*ones(1,size(u0,1))+diag([yf;dw]);
par=[par';diag(rlt(Nrlt+1).TW.PI.KP);diag(rlt(Nrlt+1).TW.PI.TI);u0(1);u0(2)];
for i=1:Ny+Nw
[~,x1]=ode15s(@(t,x)vandevusse_PIctrl(t,x,u(:,i),1,[x0;0;0],par),t,[x0;0;0]);
y1_nlin(:,:,i)=x1(:,2:3)-ones(size(t,1),1)*x0(2:3)';
end
sys=ss(rlt(Nrlt+1).TW.PI.Gsp.A,[rlt(Nrlt+1).TW.PI.Gsp.B,rlt(Nrlt+1).TW.PI.Gload.B],rlt(Nrlt+1).TW.PI.Gsp.C,[]);
y1_lin=step(sys,t,stepDataOptions('StepAmplitude',[yf;dw]));

load(rlt(Nrlt+1).data,'u0','x0','par');
y2_nlin=zeros(size(t,1),Ny,Ny);
u=[x0(2:3);u0(3:4)]*ones(1,size(u0,1))+diag([yf;dw]);
par=[par';diag(rlt(Nrlt+1).ILMI.PI.KP);diag(rlt(Nrlt+1).ILMI.PI.TI);u0(1);u0(2)];
for i=1:Ny+Nw
load(rlt(Nrlt+1).data,'u0','x0');
[~,x2]=ode15s(@(t,x)vandevusse_PIctrl(t,x,u(:,i),1,[x0;0;0],par),t,[x0;0;0]);
y2_nlin(:,:,i)=x2(:,2:3)-ones(size(t,1),1)*x0(2:3)';
end
sys=ss(rlt(Nrlt+1).ILMI.PI.Gsp.A,[rlt(Nrlt+1).ILMI.PI.Gsp.B,rlt(Nrlt+1).ILMI.PI.Gload.B],rlt(Nrlt+1).ILMI.PI.Gsp.C,[]);
y2_lin=step(sys,t,stepDataOptions('StepAmplitude',[yf;dw]));

%Rendering figure
figh=figure(Nrlt+1); pos=0;
ytick{1}=(-0.6*yf(1):0.3*yf(1):1.2*yf(1))';
ytick{2}=(-0.2*yf(2):0.3*yf(2):1.6*yf(2))';

%Ploting results
for i=1:Ny
	for j=1:Ny+Nw
		pos=pos+1;
		ax=subplot(Ny,Ny+Nw,pos);
		plot(ax,t,y1_nlin(:,i,j),'b-',t,y1_lin(:,i,j),'b--',t,y2_nlin(:,i,j),'r-',t,y2_lin(:,i,j),'r--');
        if i==1
			title(ax,InputName{j},'FontName','Helvetica','FontSize',8);
		end
		xlabel(ax,'t [min]','FontName','Helvetica','FontSize',8,'HorizontalAlignment','center','Rotation',0);
		set(ax,'XLim',[0 t(end)],'XTick',(0:20:t(end))','YLim',[ytick{i}(1) ytick{i}(end)],'YTick',ytick{i},'FontName','Helvetica','FontSize',8);
        hold on;
	end
end

ax=subplot(Ny,Ny+Nw,1);
ylabel(ax,OutputName{1},'FontName','Helvetica','FontSize',8);

ax=subplot(Ny,Ny+Nw,1+Ny+Nw);
ylabel(ax,OutputName{2},'FontName','Helvetica','FontSize',8);

figh.Units='centimeters';
figh.Color='none';
figh.OuterPosition=[0 0 19 10.5];
figh.Position=[0 0 19 10.5];

%Set legend
h=legend('TW','TW linear','ILMI','ILMI linear');
h.Position=[0.40,0.955,0.2,0.05];
h.Orientation='horizontal';
h.FontName='Helvetica';
h.FontSize=8;

%Saving figure
print(figh,rlt(Nrlt+1).data(1:end-4),'-dsvg');
print(figh,rlt(Nrlt+1).data(1:end-4),'-depsc');
%--------------------------------------------------------------------
%{
Figure 2. Output step response of Van de Vusse non-isothermal
CSTR non-linear model for TW and ILMI PI tuning for gL2=100
%}
Nrlt=Nrlt+1;
load(rlt(Nrlt+1).data,'u0','x0','par');
y1_nlin=zeros(size(t,1),Ny,Ny);
u=[x0(2:3);u0(3:4)]*ones(1,size(u0,1))+diag([yf;dw]);
par=[par';diag(rlt(Nrlt+1).TW.PI.KP);diag(rlt(Nrlt+1).TW.PI.TI);u0(1);u0(2)];
for i=1:Ny+Nw
[~,x1]=ode15s(@(t,x)vandevusse_PIctrl(t,x,u(:,i),1,[x0;0;0],par),t,[x0;0;0]);
y1_nlin(:,:,i)=x1(:,2:3)-ones(size(t,1),1)*x0(2:3)';
end
sys=ss(rlt(Nrlt+1).TW.PI.Gsp.A,[rlt(Nrlt+1).TW.PI.Gsp.B,rlt(Nrlt+1).TW.PI.Gload.B],rlt(Nrlt+1).TW.PI.Gsp.C,[]);
y1_lin=step(sys,t,stepDataOptions('StepAmplitude',[yf;dw]));

load(rlt(Nrlt+1).data,'u0','x0','par');
y2_nlin=zeros(size(t,1),Ny,Ny);
u=[x0(2:3);u0(3:4)]*ones(1,size(u0,1))+diag([yf;dw]);
par=[par';diag(rlt(Nrlt+1).ILMI.PI.KP);diag(rlt(Nrlt+1).ILMI.PI.TI);u0(1);u0(2)];
for i=1:Ny+Nw
load(rlt(Nrlt+1).data,'u0','x0');
[~,x2]=ode15s(@(t,x)vandevusse_PIctrl(t,x,u(:,i),1,[x0;0;0],par),t,[x0;0;0]);
y2_nlin(:,:,i)=x2(:,2:3)-ones(size(t,1),1)*x0(2:3)';
end
sys=ss(rlt(Nrlt+1).ILMI.PI.Gsp.A,[rlt(Nrlt+1).ILMI.PI.Gsp.B,rlt(Nrlt+1).ILMI.PI.Gload.B],rlt(Nrlt+1).ILMI.PI.Gsp.C,[]);
y2_lin=step(sys,t,stepDataOptions('StepAmplitude',[yf;dw]));

%Rendering figure
figh=figure(Nrlt+1); pos=0;
ytick{1}=(-0.6*yf(1):0.3*yf(1):1.2*yf(1))';
ytick{2}=(-0.2*yf(2):0.3*yf(2):1.6*yf(2))';

%Ploting results
for i=1:Ny
	for j=1:Ny+Nw
		pos=pos+1;
		ax=subplot(Ny,Ny+Nw,pos);
		plot(ax,t,y1_nlin(:,i,j),'b-',t,y1_lin(:,i,j),'b--',t,y2_nlin(:,i,j),'r-',t,y2_lin(:,i,j),'r--');
        if i==1
			title(ax,InputName{j},'FontName','Helvetica','FontSize',8);
		end
		xlabel(ax,'t [min]','FontName','Helvetica','FontSize',8,'HorizontalAlignment','center','Rotation',0);
		set(ax,'XLim',[0 t(end)],'XTick',(0:20:t(end))','YLim',[ytick{i}(1) ytick{i}(end)],'YTick',ytick{i},'FontName','Helvetica','FontSize',8);
        hold on;
	end
end

ax=subplot(Ny,Ny+Nw,1);
ylabel(ax,OutputName{1},'FontName','Helvetica','FontSize',8);

ax=subplot(Ny,Ny+Nw,1+Ny+Nw);
ylabel(ax,OutputName{2},'FontName','Helvetica','FontSize',8);

figh.Units='centimeters';
figh.Color='none';
figh.OuterPosition=[0 0 19 10.5];
figh.Position=[0 0 19 10.5];

%Set legend
h=legend('TW','TW linear','ILMI','ILMI linear');
h.Position=[0.40,0.955,0.2,0.05];
h.Orientation='horizontal';
h.FontName='Helvetica';
h.FontSize=8;

%Saving figure
print(figh,rlt(Nrlt+1).data(1:end-4),'-dsvg');
print(figh,rlt(Nrlt+1).data(1:end-4),'-depsc');
%--------------------------------------------------------------------
%{
Figure 3. Output step response of Van de Vusse non-isothermal
CSTR non-linear model for TW and ILMI PI tuning for gL2=10
%}
Nrlt=Nrlt+1;
load(rlt(Nrlt+1).data,'u0','x0','par');
y1_nlin=zeros(size(t,1),Ny,Ny);
u=[x0(2:3);u0(3:4)]*ones(1,size(u0,1))+diag([yf;dw]);
par=[par';diag(rlt(Nrlt+1).TW.PI.KP);diag(rlt(Nrlt+1).TW.PI.TI);u0(1);u0(2)];
for i=1:Ny+Nw
[~,x1]=ode15s(@(t,x)vandevusse_PIctrl(t,x,u(:,i),1,[x0;0;0],par),t,[x0;0;0]);
y1_nlin(:,:,i)=x1(:,2:3)-ones(size(t,1),1)*x0(2:3)';
end
sys=ss(rlt(Nrlt+1).TW.PI.Gsp.A,[rlt(Nrlt+1).TW.PI.Gsp.B,rlt(Nrlt+1).TW.PI.Gload.B],rlt(Nrlt+1).TW.PI.Gsp.C,[]);
y1_lin=step(sys,t,stepDataOptions('StepAmplitude',[yf;dw]));

load(rlt(Nrlt+1).data,'u0','x0','par');
y2_nlin=zeros(size(t,1),Ny,Ny);
u=[x0(2:3);u0(3:4)]*ones(1,size(u0,1))+diag([yf;dw]);
par=[par';diag(rlt(Nrlt+1).ILMI.PI.KP);diag(rlt(Nrlt+1).ILMI.PI.TI);u0(1);u0(2)];
for i=1:Ny+Nw
load(rlt(Nrlt+1).data,'u0','x0');
[~,x2]=ode15s(@(t,x)vandevusse_PIctrl(t,x,u(:,i),1,[x0;0;0],par),t,[x0;0;0]);
y2_nlin(:,:,i)=x2(:,2:3)-ones(size(t,1),1)*x0(2:3)';
end
sys=ss(rlt(Nrlt+1).ILMI.PI.Gsp.A,[rlt(Nrlt+1).ILMI.PI.Gsp.B,rlt(Nrlt+1).ILMI.PI.Gload.B],rlt(Nrlt+1).ILMI.PI.Gsp.C,[]);
y2_lin=step(sys,t,stepDataOptions('StepAmplitude',[yf;dw]));

%Rendering figure
figh=figure(Nrlt+1); pos=0;
ytick{1}=(-0.6*yf(1):0.3*yf(1):1.2*yf(1))';
ytick{2}=(-0.2*yf(2):0.3*yf(2):1.6*yf(2))';

%Ploting results
for i=1:Ny
	for j=1:Ny+Nw
		pos=pos+1;
		ax=subplot(Ny,Ny+Nw,pos);
		plot(ax,t,y1_nlin(:,i,j),'b-',t,y1_lin(:,i,j),'b--',t,y2_nlin(:,i,j),'r-',t,y2_lin(:,i,j),'r--');
        if i==1
			title(ax,InputName{j},'FontName','Helvetica','FontSize',8);
		end
		xlabel(ax,'t [min]','FontName','Helvetica','FontSize',8,'HorizontalAlignment','center','Rotation',0);
		set(ax,'XLim',[0 t(end)],'XTick',(0:20:t(end))','YLim',[ytick{i}(1) ytick{i}(end)],'YTick',ytick{i},'FontName','Helvetica','FontSize',8);
        hold on;
	end
end

ax=subplot(Ny,Ny+Nw,1);
ylabel(ax,OutputName{1},'FontName','Helvetica','FontSize',8);

ax=subplot(Ny,Ny+Nw,1+Ny+Nw);
ylabel(ax,OutputName{2},'FontName','Helvetica','FontSize',8);

figh.Units='centimeters';
figh.Color='none';
figh.OuterPosition=[0 0 19 10.5];
figh.Position=[0 0 19 10.5];

%Set legend
h=legend('TW','TW linear','ILMI','ILMI linear');
h.Position=[0.40,0.955,0.2,0.05];
h.Orientation='horizontal';
h.FontName='Helvetica';
h.FontSize=8;

%Saving figure
print(figh,rlt(Nrlt+1).data(1:end-4),'-dsvg');
print(figh,rlt(Nrlt+1).data(1:end-4),'-depsc');
%--------------------------------------------------------------------