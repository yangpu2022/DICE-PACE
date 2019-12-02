function [J] = evalDICEPACE(x); 
%--------------------------------------------------------------------------
% this function uses as input x, which contains the saving rate (sr) and
% the abatement fraction (mu). 
% it returns J which is (-1) times the welfare. The welfare optimisation is
% thus solved by minimising J using matlab's fmincon software.

%The function can be called from optimiseDetDice.m (for optimisation)
%     (in this case, make sure save_run=0 and plotyn=0!) 
%or from run_plot_evalDICE.m which serves to make plots 
%     (set plotyn=1). 


%--------------------------------------------------------------------------

%%% TO DO



%%%% prepare

% ---- SIMULATION SETTINGS ---- %
 

plotyn=0 ; % =1 : do plots; 0: don't plot


%USER-DEFINED OPTIONS 
prodflag=1; %total factor productivity
 %1 = DICE2016; 2 = RCP8.5-based (slower increase)
carbintflag=1; %carbon intensity (business as usual CO2 emis. per $ GDP)
 %1 = DICE2016; 2 = RCP8.5-based (slower decrease)
exoglearnflag=1; % exogenous reduction in abatement cost
 %1=DICE2016; 2= user-defined (currently just a constant)
abcostg=0.01/2; %decline in abatement cost 
 %DICE has 0.025/5=0.01/2 (per year); 
diceCparflag=13; %which carbon cycle parameters to use in DICE
 %13: carbon cycle from 2013; 16: from 2016
 %%% remark: DICE2016 is tuned for very long (>centuries) time scales; not needed here. 
damflag=1; %damage function
 %1=DICE2016 (temp. dependent); 2=dependent on cumulative emissions 3=Weitzmann
 %4=DICE2016+Paris
otherforceflag=1; %non-CO2 forcing 
 %0 = no non-CO2 forcing 
 %1 = DICE2016 (but maybe with partial abatement, see abfrac_other);
 
% USER-DEFINED PARAMETERS (i.e. those often adjusted)  
abcostexp=2.6; %the exponent in the abatement cost function 
 %DICE2016 has 2.6; Grubb et al had 2.0 
transtime=30; %typical technological transition time in years (relevant if plia>0) 
climsens=3.1; %climate sensitivity (warming [K] / doubling CO2 ; DICE2016 standard is =3.1
abfrac_other=0.6; %fraction of non-CO2 forcing that can be abated 
 %(along with CO2, at no extra cost)
 %DICE2016 has 0; Helwegen et al 2019 used 0.5
damthreshold=2; %temp threshold in case of damflag==4 [K] 
damfactor=1.0; %multiplying the original damage by a factor (only if damflag==1)
rhofactor=1; %multiplying rate of pure time preference 

%%%%%%%%%%%% taking in the input from optimizeDetDICE or run_plot_evalDICE
% note: the last row of x contains settings 

load('DICEdataDump.mat');


if size(x,1)==2  %i.e. if 2 time series are given, they are abatement and saving rate
 mu=x(1,:); %abatement fraction
 nt=length(mu);
 sr=x(2,:); %saving rate
elseif size(x,1)==1 %if only one time series is given, it is abatement
 mu=x(1,:); 
 nt=length(mu);
 sr=sr_in*ones(1,nt); %then we need to prescribe saving rate
end

tvec=(0:nt-1)*dt+t0; 



%%%%%%%%%%% OTHER MODEL PARAMETERS AND EXOGENOUS FUNCTIONS %%%%%%%%%%%%%%%%

%%%%% population %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

popadj=0.134/5*dt; %population adjustment rate
popasym=11500; %asympt. value of population (in million)
L     = zeros(1,nt); %the population in million
L(1)  = 7403; 
for t=2:nt
    L(t) = L(t-1)*(popasym/L(t-1))^popadj; 
end


%%%%% capital depreciation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

delK=(0.0633+0.0012*dt+0.00025*dt^2)/0.756; %effectively 0.1 for time step 5 yr
 %one must adjust delK depending on time step such as to get the same
 %effective capital growth for the same saving rate.
 

%%%%% total factor productivity and output %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

gamma = 0.3; %elasticity capital vs labour in production function 


if prodflag==1 %DICE2016 original
  A=zeros(1,nt);
  g_A0 = 0.076/5; % adjusted general time step
  delta_A = 0.005; 
  A(1)     = 5.115/(1000^(1-gamma)); 
   %note: I express output in the 10^12 $, and labour in 10^6 persons;
   %hence normalisation by 1000^(1-gamma) w.r.t. original DICE
  for t=2:nt
    A(t) = A(t-1)/(1 - dt*g_A0*exp(-delta_A*dt*(t-1))); % 
  end
  
elseif prodflag==2 %inspired by RCP8.5 (but same init. cond. as DICE2016 
    A=zeros(1,nt);
    A(1)     = 5.115/(1000^(1-gamma)); 
    for t=2:nt
      A(t)=A(1)*exp(3.66/pi*atan((t-1)*dt/100));
    end
end


    

%%%%% carbon intensity %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


sigma=zeros(1,nt);
if carbintflag==1 %DICE original
  gsigma0= -0.0152;
  dsigma=-0.001;
    E0=35.85*(12/44); %ref. init. emission (Gt(C)/yr) - nordhaus had GtCO2, -> *(12/44)
    Y0=105.5; %ref. init. output (trillion $2010)
    mu0=0.03; %ref. init. abatement.
  sigma(1) = E0/Y0/(1-mu0);
  for t=2:nt
    gsigma=gsigma0*(1+dsigma)^((t-1)*dt);
    sigma(t)=sigma(t-1)*exp(gsigma*dt);   
  end

elseif carbintflag==2 %inspired by RCP8.5: slower decay carb.int.    
    E0=35.85*(12/44); %ref. init. emission (Gt(C)/yr) - nordhaus had GtCO2, -> *(12/44)
    Y0=105.5; %ref. init. output (trillion $2010)
    mu0=0.03; %ref. init. abatement.
  sigma(1) = E0/Y0/(1-mu0);   
  for t=2:nt
    sigma(t)=sigma(1)*exp(-0.0057*(t-1)*dt);  %in Gt(C)/10^12$   (or kgC/$)
  end
end 


%%%%% abatement cost %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%pliability (plia) was set above

if exoglearnflag==1 %DICE-like exogenous learning for abatement cost
  abcostfac0=550/2.6/1000*44/12; %initial ab.cost factor, in $/kgC or 10^12$/GtC 
    %(DICE has 550/2.6 $/tCO2, hence I normalise) 
  abcostfac=zeros(1,nt);
  for t=1:nt
  abcostfac(t)=abcostfac0*(1-abcostg)^((t-1)*dt); 
  end

elseif exoglearnflag==2 %user-defined; here: no exogenous learning at all 
  abcostfac0=550/2.6/1000*44/12; %initial ab.cost factor, in $/kgC or 10^12$/GtC 
    %(DICE has 550/2.6 $/tCO2, hence I normalise) 
          %%%%% REMARK: Grubb working paper seems to have abcostfac0=0.024 
          %%%%%            i.e. about 1/3 of DICE2016
  abcostfac=zeros(1,nt);
  for t=1:nt
  abcostfac(t)=abcostfac0;
  end
end

%%%%% damage function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% grubbs paper suggests extra 2% GDP loss per year for 500GtC cum emiss
% above current
% (=1K -> 2K by 2050), but in general a quadratic form. so i guess
%D = damEC *EC^2  (EC = cum emiss, damCE constant) 
% damEC such that damEC *(current EC + 500GtC) - damEC *(current EC) = 2% GDP 

if damflag==1 %DICE-like (D = damT * T^2)
    damT  = 0.00236*damfactor;  %quadratic term 
elseif damflag==2 %depends on cum. emission
    EC0=500; % GtC emitted til now (roughly) 
    EC1=1000; %GtC emitted extra (reference point) 
    damEC = 0.02 / (EC1^2-EC0^2);
elseif damflag==3 
    damweitz1=1/20.2; 
    damweitz2=1/6.08; 
    damexweitz1=2;
    damexweitz2=6.76;
elseif damflag==4 %DICE-like (D = damT * T^2) + threshold ("Paris")
    damT  = 0.00236;  %quadratic term 
    damsharpness=0.02; %how sharp the arctan is (heaviside doesn't work with optimising)
    Parisfactor=0.1; %factor on the "Paris" damage term
end

%%%%% climate system %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
 if diceCparflag==16    %carbon model parameters from 2016 
     %%% -> high atm. cab. conc. due to small ocean reservoir. Nordhaus
     %%% made this change to represent long-term (slow) carbon uptake
    
  M_AT_PI = 588; %atmosphere  \
  M_UP_PI = 360; %upper ocean  } pre-industrial carbon content (GtC)
  M_LO_PI =1720; %lower ocean /

  %Carbon transition Matrix 
  phi12=0.12/5*dt; 
  phi23=0.007/5*dt; 
  
 elseif diceCparflag==13
     %%% -> faster uptake, probably more relaistic on scales of a few
     %%% decades. more in line with C-ROADS, Joost 
     
  M_AT_PI =  588; %atmosphere  \
  M_UP_PI = 1350; %upper ocean  } pre-industrial carbon content (GtC)
  M_LO_PI =10000; %lower ocean /

  %Carbon transition Matrix 
  phi12=0.088/5*dt; 
  phi23=0.0025/5*dt; 
  
 end %diceCparflag
  
  phi = [ 1 - phi12            phi12                               0; ...
       phi12*M_AT_PI/M_UP_PI   1 - phi12*M_AT_PI/M_UP_PI - phi23   phi23; ...
       0                       phi23*M_UP_PI/M_LO_PI               1-phi23*M_UP_PI/M_LO_PI]; 
    
  nu = 3.7; %increase in radiative forcing due to doubling of CO2 concenrations
     %nu/lambda is the climate sensitivity (climsens is set in beginning)
  lambda=nu/climsens;

   %temperature dynamics parameters
  sigma1 = .1105/5*dt;
  sigma2 = .088;
  sigma3 = .025/5*dt; 
   

%%%%% land-use CO2 emissions and other forcing %%%%%%%%%%%%%%%%%%%%%%%%%%%%

%DICE2016 assumes an additional land use contribution of CO2 that does not
%fall under abatement. it is so small that it doesn't matter much. 
%DICE2016 also assumes a term for "other forcing agents" (CH4, aerosol,...)
%that are treated as exogenous and are a positive forcing. 

%In Helwegen et al. 2019 we assumed that this part can be partly abated
%(abfrac_other=0.5 i.e. 50%), and assumed 
%Other_forcing = Other_forcing_DICE * (1-abfrac_other * mu)
%with no extra costs involved for abatement other forcing 
%DICE would have abfrac_other = 0

%Land-Use CO2 emissions 

delLU=0.115;
LU=zeros(1,nt);
LU(1)= 2.6*(12/44); %annual landuse em. in GtC - nordhaus had GtCO2, -> *(12/44)
for t=2:nt
    LU(t) = (1-delLU)*LU(t-1);
end

    
% non-CO2 forcing 

Foth0=zeros(1,nt);

if otherforceflag==1  %non-CO2 forcing
  for t=1:nt 
     Foth0(t)  = min(0.5+0.5/85*dt*(t-1),1); 
  end
elseif otherforceflag==0
  Foth0(1:nt)=0;
end 

    
%%%%% time preference and utility %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

elasmu=1.45; %inequality aversion parameter
rho0=0.015*rhofactor; %pure rate of time preferece (assumed constant) 

R  = zeros(1,nt); %impatience (pure time preference) weight factor
for t=1:nt
    R(t)=1/(1+rho0)^(dt*(t-1));
end

welfare_multiplicativescaling =  1/333.508503991781;
welfare_additivescaling       =  1493.7692945736;

 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%% define arrays and set initial conditions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%% climate system incl CO2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
  M_AT  = zeros(1,nt); %CO2 in atmosphere (GtC) 
  M_UP  = zeros(1,nt); %CO2 in upper ocean + biosphere 
  M_LO  = zeros(1,nt); %CO2 in lower ocean
  FCO2  = zeros(1,nt); %CO2-induced forcing W/m^2
  F     = zeros(1,nt); %total forcing (FCO2 + Foth; the latter is defined above)
  T     = zeros(1,nt); %surface temp. change w.r.t. pre-ind. in K 
  T_LO  = zeros(1,nt); %deep ocean temp. change w.r.t pre-ind. 
  
  
 if diceCparflag==16  
  M_AT(1)  = 851;
  M_UP(1)  = 460;
  M_LO(1)  = 1740;
 elseif diceCparflag==13  %obtained running DICE2013 CO2model 5yr
    %starting from Nordhaus' init. cond. for 2010
    % remark: 864GtC (atm) is actually more realistic for 2015 than is
    % 851GtC
  M_AT(1)  = 864;
  M_UP(1)  = 1541;
  M_LO(1)  = 10010;
 end
   EC0=M_AT(1)+M_UP(1)+M_LO(1)-M_AT_PI-M_UP_PI-M_LO_PI; %initial cum.Emiss. 
   %remark: init. cum. emiss. do not fit well with obs. 
  T_LO(1)  = 0.0086;
  T(1)     = 0.85;
    

%%%%% economy (capital, output, utility, emission etc) %%%%%%%%%%%%%%%%%%%%

Y     = zeros(1,nt); %economic output in 10^12 $ / yr
Ygross= zeros(1,nt);  %gross putput (no abatement cost or damage)
I     = zeros(1,nt); %investment in  10^12 $ / yr
K     = zeros(1,nt); %capital in  10^12 $ 
  K(1)=223; %initial capital 

c     = zeros(1,nt); %consumption per capita in 10^12$ / 10^6 people per year 
D     = zeros(1,nt); %damage 
abcost= zeros(1,nt); %abatement costs

E     = zeros(1,nt); %industrial carbon emissions GtC/yr
ET    = zeros(1,nt); %total (ind+ landuse) carbon emissions GtC/yr
EC    = zeros(1,nt); %cumulative carbon emissions in GtC        
  EC(1) = EC0; %cum. emiss. in first year (20150
U     = zeros(1,nt); %utility 
W     = zeros(1,nt); %welfare (integrated utility) 




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%    RUN THE MODEL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%% RUN FOR t=1 TO GET THOSE VARIABLES THAT ARE NOT INITIALISED ABOVE

t = 1;

%%%%% climate %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% is already initiated, so no need to compute that in timestep 1
    Foth(t)=Foth0(t)*(1-abfrac_other*mu(t));
    FCO2(t) = nu * log( M_AT(t) / M_AT_PI ) / log(2) ;
    F(t) = FCO2(t) + Foth(t); 
%DICE uses last step's forcing for this year's temp., need to provide
%F(t=1) for next time step


%%%%% economy %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% gross output and emission

Ygross(t) = A(t)*K(t)^gamma*(L(t))^(1-gamma);  %gross: no climate-related costs

E(t) = (1-mu(t))*sigma(t)*Ygross(t); %industrial
ET(t) = E(t) + LU(t); %total
%cumulative at t=1: earlier we put EC(1) to some init.cond ("BEGIN of
%period t=1"); ET from t=1 will be added to EC(2) later. 

%%%%% abatement cost

abcost(t)= abcostfac(t) *sigma(t) * ...
     ((1-plia) * (mu(t))^abcostexp + ...
       plia * transtime^abcostexp /(abcostexp +1)* (abs((mu(t)-0)/dt))^abcostexp);


%%%%% damage
if damflag==1 %DICE2016  
   D(t) = damT * T(t)^2;
elseif damflag==2 %cumulated emission
  D(t) = damEC * EC(t)^2; 
elseif damflag==3 %Weitzmann
   D(t) = (damweitz1 * T(t))^damexweitz1 + (damweitz2 * T(t))^damexweitz2;
elseif damflag==4 %Paris agreement 
   D(t) = damT * T(t)^2 + Parisfactor*(1 + 2/pi * atan((T(t)-damthreshold)/damsharpness))^2;
end
    
%%%%% net output

if damflag==3
Y(t) =  Ygross(t)/(1+D(t)) - Ygross(t) * abcost(t) ;
else
Y(t) =  Ygross(t)*(1-D(t)) - Ygross(t) * abcost(t) ;    
end

%%%%% new investment and new capital

I(t) = sr(t)*Y(t);
if t<nt
K(t+1) = (1-delK)^dt*K(t)+dt*I(t);
end


%utility
c(t)   = (Y(t) - I(t))/L(t);
U(t)   = L(t)*(c(t))^(1-elasmu)/(1-elasmu);
W(t)   = U(t)*R(t);


%%%%% RUN FOR t>1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%main loop
for t = 2:nt

%%%%% climate  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%influence abatement on non-CO2 forcing (no influence if abfrac_other==0)
Foth(t)=Foth0(t)*(1-abfrac_other*mu(t));

    %atmospheric CO2 concentration
    M_AT(t) =   dt * ET(t-1) ... %emissions
              + phi(1,1)*M_AT(t-1) ... %decay
              + phi(2,1)*M_UP(t-1); %exchange with upper oceans/biosphere

    %CO2 concentration upper oceans/biosphere
    M_UP(t) =   phi(2,2) * M_UP(t-1) ... %from own
              + phi(1,2) * M_AT(t-1) +  ... %exchange atmosphere
              + phi(3,2) * M_LO(t-1); %exchange lower oceans

    M_LO(t) = phi(3,3)*M_LO(t-1) + phi(2,3)*M_UP(t-1);      


    FCO2(t) = nu * log( M_AT(t) / M_AT_PI ) / log(2) ;
    F(t) = FCO2(t) + Foth(t); 
    
    %remark: for current temp., corfing from previous period is used (i.e.
    
    %temperature
    T(t) =      T(t-1) ...
              + sigma1*(F(t-1)- lambda*T(t-1) - sigma2*(T(t-1) - T_LO(t-1))); 

    %Temperature in lower oceans
    T_LO(t) =   T_LO(t-1) + sigma3*(T(t-1) - T_LO(t-1));
 

%%%%% economy  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% gross output and emission

Ygross(t) = A(t)*K(t)^gamma*(L(t))^(1-gamma);  %gross: no climate-related costs

E(t) = (1-mu(t))*sigma(t)*Ygross(t); %industrial
ET(t) = E(t) + LU(t); %total
EC(t) = EC(t-1)+ET(t-1)*dt; %cumulative emiss. ("beginning of period t") 

%%%%% abatement cost

abcost(t)= abcostfac(t) *sigma(t) * ...
     ((1-plia) * (mu(t))^abcostexp + ...
       plia * transtime^abcostexp /(abcostexp +1)* (abs((mu(t)-mu(t-1)))/dt)^abcostexp);


%%%%% damage
if damflag==1 %DICE2016
   D(t) = damT * T(t)^2;
elseif damflag==2 %cumulated emission
  D(t) = damEC * EC(t)^2; 
elseif damflag==3 %Weitzmann
    D(t) = (damweitz1 * T(t))^damexweitz1 + (damweitz2 * T(t))^damexweitz2;
elseif damflag==4 %Paris agreement
   D(t) = damT * T(t)^2 + Parisfactor*(1 + 2/pi * atan((T(t)-damthreshold)/damsharpness))^2;
end
    
%%%%% net output

if damflag==3
Y(t) =  Ygross(t)/(1+D(t)) - Ygross(t) * abcost(t) ;
else
Y(t) =  Ygross(t)*(1-D(t)) - Ygross(t) * abcost(t) ;    
end

%%%%% new investment and new capital

I(t) = sr(t)*Y(t);
if t<nt
K(t+1) = (1-delK)^dt*K(t)+dt*I(t);
end


    
    
%%%%% consumption, utility, welfare 
c(t)   = (Y(t) - I(t))/L(t);
U(t)   = L(t)*(c(t))^(1-elasmu)/(1-elasmu);
W(t)   = U(t)*R(t);
    
end  %time step

% figure; plot(D); title(num2str(damEC)); hekla
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Wrap up %%%%%%%%%%%%%%%%%%%%%%%%%%%%

%calculate & scale total welfare
Wt = sum(W)*welfare_multiplicativescaling + welfare_additivescaling;
%convert utility to loss for use in optimization
J = -Wt;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% saving  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if saveyn==1
save(['DICE_plia_' num2str(plia*100) nameadd '.mat']); 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% plotting %%%%%%%%%%%%%%%%%%%%%%%%%%%%

% if plotyn==1
% 
% fonty=20;
% colvec1=['r','m','b','c','g','k'];
% 
%  figure; 
% set(gca,'fontsize',fonty)
% hold on
% title('abatement fraction')
% plot(tvec,x(1,:),'b:','linewidth',2)
% ylabel('fraction of CO2 emiss. avoided')
% xlabel('time (yr)')
% legend('p=0.0','p=0.5','p=1.0')
% % legend('0.0','0.5','1.0','base 0.0','base 1.0')
% 
% if size(x,1)==2
% figure; 
% set(gca,'fontsize',fonty)
% hold on
% title('saving rate')
% plot(tvec,x(2,:),'b','linewidth',2)
% ylabel('fraction GDP re-invested')
% xlabel('time (yr)')
% end 
% 
% figure; 
% set(gca,'fontsize',fonty)
% hold on
% title('econ. output in trillion $')
% plot(tvec,Y,'b','linewidth',2)
% plot(tvec,Ygross,'r','linewidth',2)
% legend('net output','gross output')
% xlabel('time (yr)')
% ylabel('output (trillion $)')
% 
% figure; 
% set(gca,'fontsize',fonty)
% hold on
% title('losses (fraction GDP) for abatement and damage')
% if damflag==3
%  plot(tvec,1-1./(1+D),'b','linewidth',2) %Weitzman
% else
%  plot(tvec,D,'b','linewidth',2)
% end
% plot(tvec,abcost,'r','linewidth',2)
% legend('damage','abatement')
% xlabel('time (yr)')
% ylabel('costs  [fraction GDP]')
% 
% % figure; 
% % set(gca,'fontsize',fonty)
% % hold on
% % title('econ. growth')
% % plot(tvec(2:end),(Y(2:end)./Y(1:end-1)).^(1/dt),'b','linewidth',2)
% % legend('growth')
% % xlabel('time (yr)')
% % ylabel('econ. growth (current Y / previous Y')
% 
% figure; 
% set(gca,'fontsize',fonty)
% hold on
% title('Emissions in GtC/yr')
% plot(tvec,LU,'g','linewidth',2)
% plot(tvec,E,'b','linewidth',2)
% plot(tvec,ET,'r','linewidth',2)
% legend('landuse','industrial','total')
% xlabel('time (yr)')
% ylabel('CO2 emissions (GtC/yt)')

% elseif climmodflag ==1 
% figure; 
% set(gca,'fontsize',fonty)
% hold on
% title('Carbon content above Pre-Industry in GtC')
% plot(tvec,M_AT-M_AT_PI,'g','linewidth',2)
% plot(tvec,M_UP-M_UP_PI,'c','linewidth',2)
% plot(tvec,M_LO-M_LO_PI,'b','linewidth',2)
% plot(tvec,EC,'m','linewidth',2)
% legend('atm','upper ocean','lower ocean','cum Emis')
% xlabel('time (yr)')
% ylabel('CO2 content above PI (GtC)')
% 
% figure; 
% set(gca,'fontsize',fonty)
% hold on
% title('Temperature change since PI (K)')
% plot(tvec,T,'g','linewidth',2)
% plot(tvec,T_LO,'b','linewidth',2)
% legend('surface','lower ocean')
% xlabel('time (yr)')
% ylabel('temp. change (K)')    
%     
% 
% 
%    hekla %just a way to make the code stop after plotting. 
%    %otherwise if you forget switching off the ploting if-statement, it
%    %tried to plot during optimisation, which is a nuissance. on the other
%    %hand, if you run it in order to plot, a little error at the end won't
%    %harm anybody. 

% %plotting if




end