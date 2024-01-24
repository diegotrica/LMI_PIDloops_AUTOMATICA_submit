function info=stepinfo_caseA(KP,TI,TD,PhiD,par,t,x0,u0,du,law)
Ny=size(KP,1); Nw=size(u0,1)-Ny;
u=[x0(2:3);u0(3:4)]*ones(1,size(u0,1))+diag(du);
y=zeros(size(t,1),Ny,Ny);
switch law
	case 'PI'
for i=1:Ny+Nw
parPI=[par';diag(KP);diag(TI);u0(1);u0(2)];
[~,x]=ode15s(@(t,x)vandevusse_PIctrl(t,x,u(:,i),1,[x0;0;0],parPI),t,[x0;0;0]);
y(:,:,i)=x(:,2:3)-ones(size(t,1),1)*x0(2:3)';
end

	case 'PID'
for i=1:Ny+Nw
parPID=[par';diag(KP);diag(TI);diag(TD);diag(PhiD);u0(1);u0(2)];
[~,x]=ode15s(@(t,x)vandevusse_PIDctrl(t,x,u(:,i),1,[x0;0;0;0;0],parPID),t,[x0;0;0;0;0]);
y(:,:,i)=x(:,2:3)-ones(size(t,1),1)*x0(2:3)';
end
end

info=struct(...
    'ITAE',cell(Ny,Ny+Nw),...
    'RiseTime',cell(Ny,Ny+Nw),...
    'SettlingTime',cell(Ny,Ny+Nw),...
    'Peak',cell(Ny,Ny+Nw),...
    'PeakTime',cell(Ny,Ny+Nw));
dt=t(2:end)-t(1:end-1);
for i=1:Ny
    for j=1:Ny
        if i==j
            info(i,j).ITAE=sum(t(1:end-1).*abs(du(i)-y(1:end-1,i,j)).*dt);
			if max(y(:,i,j))>du(i)
                info(i,j).RiseTime=t(min(find(y(:,i,j)>du(i))));
            end
        else
            info(i,j).ITAE=sum(t(1:end-1).*abs(y(1:end-1,i,j)).*dt);
        end
		info(i,j).SettlingTime=t(max(find(abs(y(2:end,i,j)-y(1:end-1,i,j))>=1e-6)));
		info(i,j).Peak=max(y(:,i,j));
		info(i,j).PeakTime=t(find(max(y(:,i,j))));
    end
    for j=1+Ny:Ny+Nw
		info(i,j).ITAE=sum(t(1:end-1).*abs(y(1:end-1,i,j)).*dt);
		info(i,j).SettlingTime=t(max(find(abs(y(2:end,i,j)-y(1:end-1,i,j))>=1e-6)));
		info(i,j).Peak=max(y(:,i,j));
		info(i,j).PeakTime=t(find(max(y(:,i,j))));
    end
end
end