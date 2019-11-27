% load data file from simulation, and plot the results. This code is
% specifically made to create figure 5(?) in the paper, using the files
% 'FD_along_p_long.mat' and 'FD_along_p_short.mat'. These files can be
% replaced with any other saved workspace from simulation.
load('FD_along_p_short.mat'); % first data file to be loaded.
pulse_width=zeros(nm,np);
btot=zeros(nm,np);
thress=0.1; % detection thresshold for biomass, below that we assume b_i=0.
for i=1:nm
    for j=1:np
        xdat(i,j,xdat(i,j,:)<thress)=0; % assign zero biomass below detection level
        bmax(i,j)=max(xdat(i,j,:)); % find maximal species
        btot(i,j)=sum(xdat(i,j,1:prm.n))/prm.n;
        for w=1:prm.n % calculate the pulse width for the current m,p parameters.
            if (xdat(i,j,w)>thress)
                pulse_width(i,j)=pulse_width(i,j)+1/prm.n;
            end
  
        end
    end
end

figure(1);
set(gcf,'position',[100,100,800,500]);
subplot(2,2,3);
hold all;
show=[1 2];
str=['xb';'+g'];
for i=show
    [val,ind]=min(pulse_width(i,pulse_width(i,:)~=0)); % find the first nonzero result to remove points on the X axis.
    plot(prm.p_range(1:ind+1),pulse_width(i,1:ind+1),str(i),'markersize',8,'linewidth',1.2)
end
axis([0,400,0,.05]);
leg={'Non grazed','Grazed'};
legend(leg,'location','northwest')
xlabel('Precipitation','fontsize',14);
ylabel('Functional Diversity (pulse width)','fontsize',14);
set(gca,'fontsize',14);
text(0.9,0.1,'(c)','units','normalized','fontsize',16);

subplot(2,2,1);
hold all;
title('Short grazing history','fontsize',14);
show=[1 2];
str=['xb';'+g'];
for i=show
    [val,ind]=min(btot(i,btot(i,:)~=0));
    plot(prm.p_range(1:ind+1),btot(i,1:ind+1),str(i),'markersize',8,'linewidth',1.2)
end
axis([0,400,0,0.015]);
leg={'Non grazed','Grazed'};
legend(leg,'location','northwest')
xlabel('','fontsize',14);
ylabel('Total Biomass','fontsize',14);
set(gca,'fontsize',14);
text(0.9,0.1,'(a)','units','normalized','fontsize',16);


% repeat the process to displat another file
load('FD_along_p_long.mat');  % second data file to be loaded.
pulse_width=zeros(nm,np);
btot=zeros(nm,np);
thress=0.1; % detection thresshold for biomass, below that we assume b_i=0.
for i=1:nm
    for j=1:np
        xdat(i,j,xdat(i,j,:)<thress)=0; % assign zero biomass below detection level
        bmax(i,j)=max(xdat(i,j,:)); % find maximal species
        btot(i,j)=sum(xdat(i,j,1:prm.n))/prm.n;
        for w=1:prm.n % calculate the pulse width for the current m,p parameters.
            if (xdat(i,j,w)>thress)
                pulse_width(i,j)=pulse_width(i,j)+1/prm.n;
            end
  
        end
    end
end

figure(1);
set(gcf,'position',[100,100,800,500]);
subplot(2,2,4);
hold all;
show=[1 2];
str=['xb';'+g'];
for i=show
    [val,ind]=min(pulse_width(i,pulse_width(i,:)~=0)); % find the first nonzero result to remove points on the X axis.
    plot(prm.p_range(1:ind+1),pulse_width(i,1:ind+1),str(i),'markersize',8,'linewidth',1.2)
end
axis([0,400,0,.05]);
leg={'Non grazed','Grazed'};
legend(leg,'location','northwest')
xlabel('Precipitation','fontsize',14);

ylabel('Functional Diversity (pulse width)','fontsize',14);
set(gca,'fontsize',14);
text(0.9,0.1,'(d)','units','normalized','fontsize',16);

subplot(2,2,2);
hold all;
title('Long grazing history','fontsize',14);
show=[1 2];
str=['xb';'+g'];
for i=show
    [val,ind]=min(btot(i,btot(i,:)~=0));
    plot(prm.p_range(1:ind+1),btot(i,1:ind+1),str(i),'markersize',8,'linewidth',1.2)
end
axis([0,400,0,0.015]);
leg={'Non grazed','Grazed'};
legend(leg,'location','northwest')
xlabel('','fontsize',14);
ylabel('Total Biomass','fontsize',14);
set(gca,'fontsize',14);
text(0.9,0.1,'(b)','units','normalized','fontsize',16);
