%clear up
% clear all;
% close all;

% exponent: 2.6
% do allow partial abatement 
% fix / dont fix 2015 abatement
% timescale 20 and 40 
% 1.5 * damage 
% other rho


% % % %set choices: 
tmax = 60; %nr of time steps
dt=5; %
t0=2015;
sr_flag=1; %=0: saving rate fixed; =1: freely optimised
sr_0=0.25; % saving rate, if not freely optimised (if using ssp data, use that saving rate sr_in!!)
plotyn=1; %=1: make plots at the end
pliavec=[0 0.5 1]; plian=length(pliavec); %pliabilities to consider
% if pliavec gets another dimension than 3, then adjust the legends in the
% plots by adding or removing entries, and in case pliavec gets >3
% elements, add entries to dashvec. 
%otherwise, the plotting will give errors. 



%prepare some input data 
 sr_in=sr_0; %saving rate in case it's not optimised


for pli=1:plian
    plia=pliavec(pli); 

 %here we put the inputs defined above into a file to be used in
 %evalDICEPACE
 
 saveyn=0; %do not save results within evalDICEPACE while optimising... 
  save('DICEdataDump.mat','sr_in', 't0','dt','plia','saveyn');
    
% set variable's initial condition and bounds

x0=[]; x0_lo=[]; x0_up=[]; 
mu=0.0*ones(1,tmax);       %abatement initial guess
  mu_lo=zeros(1,tmax);     %abatement lower bound
  mu_up=ones(1,tmax);      %abatement upper bound 
sr = 0.25*ones(1,tmax);    %saving rate
  sr_lo=zeros(1,tmax); 
  sr_up=ones(1,tmax); 
  
%if you want to fix the first time step's abatement (As Nordhaus sometime does),
% set the initial condition and upper and lower bound to that value... 
%   mu(1,1)=0.03; mu_lo(1,1)=0.03; mu_up(1,1)=0.03; 


x0=mu; x0_lo=mu_lo; x0_up=mu_up; 
if sr_flag==1;   
    x0=[x0; sr]; x0_lo=[x0_lo; sr_lo]; x0_up=[x0_up; sr_up];
end



%optimize
opt = optimset('MaxFunEvals',1000000, 'MaxIter',10000);
[x,y0,exitflag] = fmincon(@evalDICEPACE,x0,[],[],[],[],x0_lo,x0_up,[],opt);
exitflag



%now we call the function again for the optimal x, but this time with saving outputs (saveyn=1)
 saveyn=1; 
 nameadd=[]; %don't add anything to the file name when saving results in evalDICE... 
 
  save('DICEdataDump.mat','sr_in','t0','dt','plia','saveyn','nameadd')


[J]=evalDICEPACE(x);

if pli==1;
    nameadd='_no_abate';
    nameadd0=nameadd; %for later use when plotting
    x0=x; x0(1,:)=0; 
    save('DICEdataDump.mat','sr_in', 't0','dt','plia','saveyn','nameadd')
    [J]=evalDICEPACE(x0);
end

end %pliavec


%%%%%%%%%%%%%%%%%%%%%% plotting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



if plotyn==1

fonty=18;
colvec1=['r','m','b','c','g','k'];
dashvec=['- ';'--';': ']; 


figure; 
set(gca,'fontsize',fonty)
hold on
title('Losses from climate damage')
for pli=1:plian
    plia=pliavec(pli);
    if pli==1
            load(['DICE_plia_' num2str(plia*100)  nameadd0 '.mat'])
            if damflag==3
             plot(tvec,(1-1./(1+D)),['r' dashvec(pli,:)],'linewidth',2) %Weitzman
            else
             plot(tvec,D,['r' dashvec(pli,:)],'linewidth',2)
            end  
    end
    load(['DICE_plia_' num2str(plia*100)  '.mat'])
    if damflag==3
      plot(tvec,(1-1./(1+D)),['b' dashvec(pli,:)],'linewidth',2) %Weitzman
    else
      plot(tvec,D,['b' dashvec(pli,:)],'linewidth',2)
    end   
end
legend('baseline',['abate p = ' num2str(pliavec(1))],['abate p = ' num2str(pliavec(2))],...
    ['abate p = ' num2str(pliavec(3))])
xlabel('time (yr)')
ylabel('expenditure [fraction of GDP]')

figure; 
set(gca,'fontsize',fonty)
hold on
title('Abatement expenditure')
for pli=1:plian
    plia=pliavec(pli);
    if pli==1
     load(['DICE_plia_' num2str(plia*100)  nameadd0 '.mat'])
     plot(tvec,abcost,['r' dashvec(pli,:)],'linewidth',2)          
    end
    load(['DICE_plia_' num2str(plia*100) '.mat'])
    plot(tvec,abcost,['b' dashvec(pli,:)],'linewidth',2)    
end
legend('baseline',['abate p = ' num2str(pliavec(1))],['abate p = ' num2str(pliavec(2))],...
    ['abate p = ' num2str(pliavec(3))])
xlabel('time (yr)')
ylabel('GDP loss  [fraction of GDP]')

figure; 
set(gca,'fontsize',fonty)
hold on
title('Industrial Emissions in GtC/yr')
for pli=1:plian
    plia=pliavec(pli);
    if pli==1
     load(['DICE_plia_' num2str(plia*100) nameadd0 '.mat'])
     plot(tvec,E,['r' dashvec(pli,:)],'linewidth',2)          
    end
    load(['DICE_plia_' num2str(plia*100)  '.mat'])
    plot(tvec,E,['b' dashvec(pli,:)],'linewidth',2)
end
legend('baseline',['abate p = ' num2str(pliavec(1))],['abate p = ' num2str(pliavec(2))],...
    ['abate p = ' num2str(pliavec(3))])
xlabel('time (yr)')
ylabel('CO2 emissions (GtC/yt)')

figure; 
set(gca,'fontsize',fonty)
hold on
title('Carbon content above Pre-Industry in GtC')
for pli=1:plian
    plia=pliavec(pli);
    if pli==1
      load(['DICE_plia_' num2str(plia*100) nameadd0 '.mat'])
      plot(tvec,M_AT-M_AT_PI,['r' dashvec(pli,:)],'linewidth',2)
%     plot(tvec,EC,['m' dashvec(pli,:)],'linewidth',2)
    end
    load(['DICE_plia_' num2str(plia*100)  '.mat'])
    plot(tvec,M_AT-M_AT_PI,['b' dashvec(pli,:)],'linewidth',2)
%     plot(tvec,EC,['c' dashvec(pli,:)],'linewidth',2)
end
legend('baseline',['abate p = ' num2str(pliavec(1))],['abate p = ' num2str(pliavec(2))],...
    ['abate p = ' num2str(pliavec(3))])
xlabel('time (yr)')
ylabel('CO2 content above PI (GtC)')


figure; 
set(gca,'fontsize',fonty)
hold on
title('Global Mean Surface Temperature change since PI  ')
for pli=1:plian
    plia=pliavec(pli);
    if pli==1
    load(['DICE_plia_' num2str(plia*100)  nameadd0 '.mat'])
    plot(tvec,sum(T,1),['r' dashvec(pli,:)],'linewidth',2)
    end 
    load(['DICE_plia_' num2str(plia*100) '.mat'])
    plot(tvec,sum(T,1),['b' dashvec(pli,:)],'linewidth',2)
end
legend('baseline',['abate p = ' num2str(pliavec(1))],['abate p = ' num2str(pliavec(2))],...
    ['abate p = ' num2str(pliavec(3))])
xlabel('time (yr)')
ylabel('temp. change (K)')


end %plotyn
