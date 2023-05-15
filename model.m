function  [Y_real,lock_start2,lock_end,defGDP,debGDP,labforce1,labforce2,Lev_a,Lev_ab,Lev_al,Lev_af,Lev_ak,defaults,defaults_k,Hdemand,tot_healthcare,hospital_capacity,CPI_b,CPI_l,dprob,lock_duration,detect_new,inf_detected,distance,inf_share2,sus_share2,rec_share2,healthcare,serious_share,inf_cum,inf_new,inf_cum_v,inf_new_v,sus_share,inf_share,rec_share,death_total,gdp_deflator, Investment,EXP, consumption, Un,stock_bonds,wages_t, DC, totK,I,totE,bankruptcy_rate,et,redcons,dis_share,dis_prob_tt,vacc_share,epidemic_end,epidemic_end_tt,epidemic_start,death_new,vax_stats,dead_share,inf_current,detect_current,vax_stats2,dead_age1,dead_age2,dead_age3,detect_new_v] = model(seed, T, novax_share,healthcare_exp,vac_immunity,nat_immunity)

%dummies for scenarios
epidemic=1;
lockdown_exp=1;
social_distancing=1;
CIGS=0;
liquidity_transfers1=0;
liquidity_transfers2=0;
liquidity_transfers3=0;
credit_guarantees=0;
equity_support=0;
income_support=0;
vaccine=0;
vaccine_strategy_1=0;   %random vaccine
vaccine_strategy_2=0;   %priority to old
vaccine_strategy_3=0;   %priority to workers
mutation=0;
vaccine_improved=0;
mutation2=0;
no_vax_exp=0;

%numbers of agents
Bk=1;                               %no. of banks
F=600;                              %no. of C-firms
F_l=360;                            %no. of "luxury" firms
F_b=240;                            %no. of basic firms
W=30000;                            %no. of workers
N = 200;                            %no. of capital producing firms

%%main seed
rng(seed,'twister')

%%creation of random numbers
shock_pk = rand(T,N);
shock_p  = rand(T,F);
prob_k = rand(T,F);

shocks_healthcare2=rand(T*4,W);
partners_b=zeros(W+F+N,2);
partners_l=zeros(W+F+N,2);
partners_k=zeros(F,2);
partners_b(:,1)=randi(F_b,W+F+N,1);
partners_l(:,1)=F_b+randi(F_l,W+F+N,1);
partners_k(:,1)=randi(N,F,1);
rage=rand(W,1);

%age structure of population
age=zeros(1,W);
age(rage<0.15)=1;       %young
age(rage>0.15)=2;       %middle-aged
age(rage>0.8)=3;        %old 

%economic model parameters
r_f=0.01;                          %general refinancing rate
r_f=r_f/3;                         % dividing by 3 to get a rough monthly rate
r_d=r_f/2;                          % deposit rate (markdown on risk-free rate)
z_c=2;                              %no. of aplications in consumption good market
z_k = 2;                         %no. of aplications in capiatal good market
z_e = 5;                            %number of job applications  
xi = 0.55;%0.5,0.05                          %memory parameter human wealth
chi = 0.00835;%0.00825,0.05                            %fraction of wealth devoted to consumption	
q_adj = 0.2;%par(6);%0.2                            %quantity adjustment parameter	
p_adj = 0.08;%par(7);0.07                           %price adjustment parameter    
mu =  1.0070;                           %bank's gross mark-up
eta = 0.0100;%0.01                         %capital depreciation
Iprob=0.4;%par(10)                      %probability of investing
phi =  0.0025;                        %bank's leverage parameter
theta=0.005;%par(12);                         %rate of debt reimbursment
delta = 0.2; %par(13);%0.2                        %memory parameter in the capital utilization rate
%k = 1/9;
k = 1/6;%8/9; %1                           %capital productivity
alpha=2/3;%2/9;                              %labour productivity
div =0.2500;                            %share of dividends
div_B=0.2500;                                %bank dividends
barX=0.8500;%0.85;                          %desired capital utilization
inventory_depreciation = 0.0800;%0.08;0.3           %rate at which capital firms' inventories depreciate
b1 = -10;   %-15;
b2 = 10;   %13;                      %Parameters for risk evaluation by banks
b_k1 = -15; %-5;
b_k2 = 15; %5 ;
interest_rate = 0.003333333;
subsidy = 0.75;%par(24);0.75,0.40          %unemployment subsidy --> now roughly in line with net replacement rate for Italy acc. to OECD
tax_rate = 0.2750;  %0.15                %tax rate --> needs to be changed in line with unemployment & inactive subsidy & healthcare spending
wage_update_up = 0.0330;
wage_update_down = 0.0033;
u_target = 0.1000;%0.10, 0.05
bond_interest_rate = interest_rate;   
unemployment_subsidy_init = subsidy;
inactive_sub=1.2;                     %subsidy received by inactives (as fraction of unemployment benefit) --> in line with pensions for Italy
tax_rate_f = 0.3000;        %taxes on dividends --> needs to be calibrated in line with gov. spending on subsidies & healthcare
tax_rate_k = 0.3000;
tax_rate_b = 0.3000;
b = [b1;b2];
b_k = [b_k1;b_k2];
targetLev=0.2;
entry=-40;
price_sensitivity=40;
mu_0=[0.15 0.15]; %0.15,0.5
wb=1/3;                                 % initial wage rate

%epidemic model parameters
ncons=((W*(W-1)/2));   %maximum number of possible connections
n_perm_cons=round(ncons/15000); %number of permanent connections 
shock_l=0.00165;
shock_b=0.00055;
reduce_cons_lockdown=0.25;
reduce_workcons_lockdown=0.375;
tr_annual=zeros(1,12);
tr_annual(1:4)=0.07;%0.0825,0.085,0.0875
tr_annual(5:9)=0.04;%0.0575,0.05,0.0525
tr_annual(10:12)=0.07;%0.0825,0.085,0.0875
transmission_rate=repmat(tr_annual,1,73);
%transmission_rate=0.075;
distancing_effect=1/3;
lockdown_threshold=30;%25
lock_min_duration=3;
lock_constraint=1;
lock_end_threshold=12.5;
post_lock_adjustment=0.0775;%0.08;%1/3;
distancing_persistence=0.725;%0.8;%0.8;%0.8
intensity=[0.05 0.5 1];
dis_cost=6;
dis_cost_lockdown=-6;
dis_threshold=5;
detectprob=0.02; 
max_detectprob=0.125;
dieprob=0.0075;%0.01
detect_adjustment=0.0005;%0.00005
seriousprob=[0.01 0.02 0.525];%0.02 0.06 0.3
h_1=2;%2.2
h_2=0.1;
h_3=0.0375;%0.035
n_shopcons=3;
susprob=0.0025;
%"normal" disease
infprob2=0.0012;
susprob2=0.1;
max_dur2=4*ones(1,W);

%vaccine parameters
coverage_rate=0.01;                                 %coverage rate
vaccine_efficacy1=0.8;%0.7                              %vaccine efficacy wrt susceptibility
vaccine_efficacy2=0.95;%0.9                              %vaccine efficacy wrt pathology
%average_immunity_weeks=52;                          %length of vaccine-induced immunity
average_immunity_weeks=vac_immunity;
%natural_immunity=52;
natural_immunity=nat_immunity;

%variant parameters
mutation_start=2348;
%transmission_rate_v=transmission_rate.*1.5;
tr_annual_v=zeros(1,12);
tr_annual_v(1:4)=0.105;
tr_annual_v(5:9)=0.095;
tr_annual_v(10:12)=0.105;
transmission_rate_v=repmat(tr_annual_v,1,73);
dieprob_v=1;    
distancing_effect_v=distancing_effect*0.25;
vaccine_efficacy1_v=0.8*vaccine_efficacy1*ones(1,W);        %vaccine efficacy 1 with variant
if mutation2==1
    vaccine_efficacy2_v=0.8*vaccine_efficacy2*ones(1,W);    %vaccine efficacy 2 with variant
else
    vaccine_efficacy2_v=vaccine_efficacy2*ones(1,W);
end
%initialise agents' states
active_f=ones(F,1);	
active_k=ones(N,1);
active=ones(1,W);       %which agents are active?
active(age==3)=0;       %Note that here (unless there is an epidemic) old and inactive codincide! i.e. there are no younger inactive! Could be changed if desired
sickpay=zeros(1,W);     %to indicate which agent receives sickpay
dead=zeros(1,W);        %to record dead agents
susceptible2=ones(1,W);
infected2=zeros(1,W);
recovered2=zeros(1,W);
susceptible=ones(1,W);
infected=zeros(1,W);
detected=zeros(1,W);
recovered=zeros(1,W);
duration=zeros(1,W);
duration2=zeros(1,W);
recovered_duration=zeros(1,W);
recovered_max=zeros(1,W);
serious=zeros(1,W);
bankrupt_f=zeros(1,F);
locked_f=zeros(F,1);
bankrupt_k=zeros(1,N);
homeoffice=zeros(1,F);
homeoffice_k=zeros(1,N);
distancing=zeros(1,W);
health_demanded=zeros(1,W);
health_demanded2=zeros(1,W);
health_supplied=zeros(1,W);
health_supplied2=zeros(1,W);
work_cons=[];
labforce=sum(active);   %size of labour force
vaccinated=zeros(1,W);
immunity_duration=zeros(1,W);
variant=zeros(1,W);
immunity_weeks=zeros(1,W);

%epidemic time series
inf_share2=zeros(1,T*4);   %share of infected (normal disease)
sus_share2=zeros(1,T*4);   %share of susceptible (normal disease)
rec_share2=zeros(1,T*4);   %share of recovered (normal disease)
Hdemand=zeros(1,T*4);
inf_current=zeros(1,T*4);
detect_current=zeros(1,T*4);
sus_share=zeros(1,T*4);   %share of susceptible
inf_share=zeros(1,T*4);   %share of infected
inf_detected=zeros(1,T*4);   %share of infected & detected
rec_share=zeros(1,T*4);   %share of recovered
dead_share=zeros(1,T*4);  %share dead
inf_cum=zeros(1,T*4);     %cumulative cases
inf_new=zeros(1,T*4);     %new cases
inf_newp=zeros(1,T*4);
detect_new=(zeros(1,T*4));
serious_share=zeros(1,T*4);
distance=zeros(1,T*4);
dprob=detectprob*ones(1,T*4);
redcons=zeros(1,T*4);
dis_share=zeros(1,T*4);
health_rationing=zeros(1,T*4);
health_rationing2=zeros(1,T*4);
dis_prob_tt=zeros(T*4,W);
vacc_share=zeros(1,T*4);
inf_detected_v=zeros(1,T*4); 
inf_cum_v=zeros(1,T*4);     
inf_new_v=zeros(1,T*4);     
detect_new_v=(zeros(1,T*4));
death_new=(zeros(1,T*4));
death_total=(zeros(1,T*4));
dead_age1=zeros(1,T*4); 
dead_age2=zeros(1,T*4); 
dead_age3=zeros(1,T*4); 

%Policy
labour_transfer=zeros(1,F);
labour_transfer_k=zeros(1,N);
liquidity_transfer=zeros(1,F);
liquidity_transfer_k=zeros(1,N);
loans_g=0;
deb_g=zeros(1,F);
deb_g_k=zeros(1,N);
equity_injection=zeros(1,F);
equity_injection_k=zeros(1,N);
gov_ownership=zeros(1,F);
gov_ownership_k=zeros(1,N);
div_g=zeros(1,F);
div_g_k=zeros(1,N);
defaults_gov=0;

%Economic model time series
Y_tot_0=(labforce-u_target*labforce)*alpha;
Y_firm_0=Y_tot_0/(F+N);

price_k=zeros(1,(T+1));
G=zeros(1,T);                       %government expenditures on subsidies
Healthcare=zeros(1,T);              %nominal healthcare spending
TA=zeros(1,T);                      %government income
GB=zeros(1,T);                      %governament budget  GB = TA - subsidies - EXP - Healthcare - interest
DC=zeros(1,T);
Lev_a=zeros(1,T);
Lev_af=zeros(1,T);
Lev_ab=zeros(1,T);
Lev_al=zeros(1,T);
Lev_ak=zeros(1,T);
labforce1=zeros(1,T);
labforce2=zeros(1,T);
defGDP=zeros(1,T);
debGDP=zeros(1,T);
average_interest_rate=zeros(1,T);
average_interest_rate_k=zeros(1,T);
pub_exp_c= zeros(1,T); 
INT=zeros(1,T);
stock_bonds=zeros(1,T);
Y_nominal_k = NaN(1,T);
Y_nominal_c = NaN(1,T);
Y_nominal_tot= NaN(1,T);
Y_real =NaN(1,T);
%Y_real(1)=Y_tot_0;
gdp_deflator = NaN(1,T);
I = NaN(1,T);
bonds = zeros(1,T);
healthcare=zeros(1,T);                     %healthcare spending time series
totE=zeros(1,T);  
Un=zeros(1,T);                      %unemployment
wages_t = NaN(1,T);                 %wages    
dividends=zeros(1,T);               %total dividends
consumption=zeros(1,T);             %total cunsumption                     
price=zeros(1,(T+1));                %consumer price index time series
defaults=zeros(1,T);                %number of bankruptcies
profitsB=zeros(T,Bk);               %bank profit time series
dividendsB=zeros(T,Bk);
unsatisfiedDemand = zeros(1,T);
totK = zeros(1,T);
Investment = NaN(1,T);
defaults_k=zeros(1,T); 
bankruptcy_rate=zeros(1,T);
CPI_l=zeros(1,T);
CPI_b=zeros(1,T);
active_f_ratio=zeros(1,T);
active_k_ratio=zeros(1,T);
active_tot_ratio=zeros(1,T);
times=zeros(1,T);
exp_c=zeros(F,T);         % valore EROGATO individualmente per updating liquidity
quota_exp_c=zeros(F,T);   % quota singola erogabile impresa per totale
public_dem_c=zeros(F,T);  % domanda pubblica per le consumption
quota_health=zeros(F+N,T); %share of healthcare spending going to each firm

Inventories=zeros(1,T);
gross_investment_share=zeros(1,T);
net_investment_share=zeros(1,T);
inventories_share=zeros(1,T);
consumption_share=zeros(1,T);

%Initialisation Economic model
healthcare(:)=0.04*labforce*alpha;         %healthcare spending is 4% of full employment output (needs to be calibrated but note that this shouldn't be calibrated to total healthcare spending over GDP!)
tot_healthcare=healthcare(1)*ones(1,T*4);
if healthcare_exp==1
    healthcare(578:length(healthcare))=0.05*labforce*alpha;
    tot_healthcare(2309:length(tot_healthcare))=0.05*labforce*alpha;
end
hospital_capacity=zeros(1,T*4);
Gov = 0;
quota_exp_c(:,1)=1/F;     % time 1: occhio se cambi iniziando con t=2
quota_health(:,1)=1/(F+N);
RIC=ones(1,F);
RIC_k=ones(1,N);

%capital firm
Leff_k = zeros(1,N);
Y_k =zeros(1,N);
Y_prev_k=Y_firm_0*ones(1,N);
Y_kd=Y_prev_k;
P_k=(wb/alpha)*(1+mu_0(1))*ones(1,N); %mu_0=1 -> P_k=3
A_k =30*ones(1,N);
A_control=A_k;
liquidity_k = A_k;
De_k=Y_prev_k;
deb_k=zeros(1,N);
price_k(1:2)=mean(P_k);                %capital price index
Q_k=zeros(1,N);                     %sales 
Ftot_k=zeros(1,N);
interest_r_k=zeros(1,N);

%firms
K_0=(Y_firm_0/k)/barX;
K=K_0*ones(1,F);
A=30*ones(1,F)+ K*price_k(1); 
liquidity=A-K*price_k(1);
capital_value = K*price_k(1);
value_investments=zeros(1,F);
investment=zeros(1,F);
P=(wb/alpha)*(1+mu_0(2))*ones(1,F);                    %past production
Y_prev=Y_firm_0*ones(1,F);
Yd=Y_prev;
Q=zeros(1,F);                       %sales   
Leff=zeros(1,F);                    %current employees
De=Y_prev;                       %expected demand
deb=zeros(1,F);                     %firm debts
barK=K;
barYK=Y_prev/k;
x = barYK./barK;
interest_r=zeros(1,F);

%households
w=zeros(1,W);                       %earned wages
PA=1/3*ones(1,W+F+N);                 %household personal assets (saving stock)
Oc=zeros(1,W);                      %employment status: Oc(i)=j --> worker j employed by firm i; if i=0, j is unemployed

%total deposits of households and firms
deposits=sum(liquidity)+sum(liquidity_k)+sum(PA);
loans=sum(deb);                            

%government
GB(1)=-(1000+deposits);
EXP = Gov*ones (1,T);       % spesa pubblica EROGABILE, update erogata in fondo
EXP(1,1)=0;               % inizializzazione, vedi sotto
bonds(1)=-GB(1);
stock_bonds(1)=bonds(1);
A_init=A;
A_init_k=A_k;

%bank
E=loans+bonds(1)-deposits;
totE(1:2)=E;


price(1:2)=P(1);
Ftot=zeros(1,F);                    %borrowings
CPI_l(1)=mean(P);                         %CPI for luxury goods
CPI_b(1)=mean(P);                         %CPI for basic goods
dividends_income_k = zeros(1,N);
dividends_income   = zeros(1,F);
permanent_income = 1/3*ones(F+W+N,1)./price(1);
permanent_income(W+1:W+F)=0.1/price(1);

Y=zeros(1,F);	
stock=zeros(1,F);	
stock_k=zeros(1,N);	
pi=zeros(1,F);	
pi_k=zeros(1,N);


%Initialisation epidemic model
dis_index=-3*ones(1,W);
dis_prob=zeros(1,W);
constraint_k=ones(1,N);         %degree to which K-production is constrained during lockdown
constraint_f=ones(1,F);         %degree to which L-firms that implement smart working reduce production
lockdown=0;
lf_share=round(F_l/3);
deactivate_f=F_b+1:F_b+lf_share;
lock_start=0;
reduce_cons=1;          
lock_duration=0;
lock_ended=0;
c_shock_l=ones(1,W);
c_shock_b=ones(1,W);
dead_assets=0;
epidemic_end=T+1;
epidemic_end_tt=T*4;
a_v1=0;
a_v2=0;
a_v3=0;
tt=5;
lock_start2=0;
lock_end=0;
epidemic_start=T;

tic()
 for t=2:T 
      %disp(t)
      
     if t==epidemic_end
        clear('connections','con','cons','colleagues','customers','fixed_cons','perm_cons','shop_cons','work_cons','shop_encounters','shop_weights','perm_weights','work_weights')
     end
     
     if epidemic==1 && t>577
        deadA=find(dead==1);
        for i=1:length(deadA)
            if rand<0.0125
                dead(deadA(i))=0;
                active(deadA(i))=1;
                Oc(deadA(i))=0;
                age(deadA(i))=1;
                susceptible(deadA(i))=1;
                infected(deadA(i))=0;
                detected(deadA(i))=0;
                vaccinated(deadA(i))=0;
                %if rand<seriousprob(age(deadA(i)))
                if rand<seriousprob(1)
                  serious(deadA(i))=1;
                else
                  serious(deadA(i))=0;
                end
            end
        end    
     end
     
     if epidemic==1 && t==577
        shocks_healthcare=rand(T*4,W);
        max_dur=randi([3,5],1,W);              % maximum duration of the disease
        for i=1:length(serious)
            if age(i)==1
               if rand<seriousprob(1)     %probabilities of getting serious symptoms
                  serious(i)=1;
               end
            end
            if age(i)==2
               if rand<seriousprob(2)     
                  serious(i)=1;
               end
            end
            if age(i)==3
               if rand<seriousprob(3)     
                  serious(i)=1;
               end
            end
        end
        
        perm_cons=zeros(n_perm_cons,2);
        for i=1:n_perm_cons
            con=randperm(W,2);
            perm_cons(i,:)=con;
        end
        perm_weights=3*ones(length(perm_cons),1);
        perm_cons=[perm_cons perm_weights];
        
        no_vax=zeros(1,W);
        if no_vax_exp==1
            no_vax_share=novax_share;
            W_input=round(W*no_vax_share);
            id_nv=randperm(W,W_input);          %random selection of novax agents
            no_vax(id_nv)=1;
        end

     end
    
     %% Lockdown ends
     if lockdown==1 && t-lock_start>=lock_min_duration && mean(detect_new(tt-1:tt))<=lock_end_threshold
        lockdown=0;                 
        sumY=sum(Yd(F_b+1:F_b+F_l));
        %divisor=F_l-sum(active_f==0);
        meanY=sumY/F_l;
        indexf=1:F;
        P(active_f==0 & bankrupt_f'==0 & indexf'>F_b)=mean(P(active_f==1 & bankrupt_f'==0 & indexf'>F_b));
        RIC(active_f==0 & bankrupt_f'==0 & indexf'>F_b)=mean(RIC(active_f==1 & bankrupt_f'==0 & indexf'>F_b));
        YP=Y_prev;
        YP(YP>meanY)=meanY;
        Y_prev(F_b+1:F_b+F_l)=YP(F_b+1:F_b+F_l);
        Y_prev(active_f==0 & bankrupt_f'==0 & indexf'>F_b)=meanY;
        Yd(F_b+1:F_b+F_l)=Y_prev(F_b+1:F_b+F_l);
        barYK(F_b+1:F_b+F_l) = Y_prev(F_b+1:F_b+F_l)/k;
        Y_prev(bankrupt_f==1)=0;
        Yd(bankrupt_f==1)=0;
        active_f(bankrupt_f==0)=1;
        locked_f(:)=0;
        lock_duration=tt-lock_start2;
        lock_end=tt;
        lock_ended=1;
    end
    
    if lock_ended==1
       for i=1:F
           if homeoffice(i)==1 && rand<post_lock_adjustment
              homeoffice(i)=0;
              constraint_f(i)=1;
           end
       end
       for i=1:N
           if homeoffice_k(i)==1 && rand<post_lock_adjustment
              homeoffice_k(i)=0;
              constraint_k(i)=1;
           end
       end  
    end
    
%     if lockdown==1
%         disp('lock')
%     end
    
    inactive_firms=find(active_f==0 & locked_f==0);
     for j=1:length(inactive_firms)	
        i=inactive_firms(j);	
        if i<=F_b	
           prob_e=1/(1+exp(entry*profit_rate_b+2));	
        else	
           prob_e=1/(1+exp(entry*profit_rate_l+2));	
        end	
        if rand<=prob_e	
        A(i)=PA(W+i)+K(i)*price_k(t);                                                   %initialize new firm 	
        capital_value(i)=K(i)*price_k(t);	
        PA(W+i)=0;	
        liquidity(i)=A(i)-K(i)*price_k(t);	
        deb(i)=0;	
        P(i)=mean(P);	
        mxY=((A(i)+targetLev*A(i)/(1-targetLev))*1/wb)*alpha;	
        Y_prev(i)=min(trimmean(Y_prevp,10),mxY);	
        Yd(i)=Y_prev(i);	
        x(i)=Y_prev(i)/k/K(i);	
        barK(i)=K(i);	
        barYK(i)=Y_prev(i)/k;	
        Y(i)=0;	
        stock(i)=0;	
        interest_r(i)=r_f;	
        active_f(i)=1;	
        bankrupt_f(i)=0;
        end	
     end	
     	
     inactive_kfirms=find(active_k==0);	
     for j=1:length(inactive_kfirms)	
        if rand<1/(1+exp(entry*profit_rate_k+2))	
        i=inactive_kfirms(j);	
        A_k(i)=PA(W+F+i);                                                   %initialize new firm 	
        A_control(i)=A_k(i);	
        PA(W+F+i)=0;	
        liquidity_k(i)=A_k(i);	
        P_k(i)=mean(P_k);	
        mxY=((A_k(i)+targetLev*A_k(i)/(1-targetLev))*1/wb)*alpha;	
        Y_prev_k(i)=min(trimmean(Y_prevp_k,10),mxY);	
        Y_kd(i)=Y_prev_k(i);	
        Y_k(i)=0;	
        stock_k(i)=0;	
        interest_r_k(i)=r_f;	
        active_k(i)=1;	
        bankrupt_k(i)=0;
        end	
     end
     
     active_f_ratio(t)=sum(active_f==1)/(F);
     active_k_ratio(t)=sum(active_k==1)/(N);
     active_tot_ratio(t)=(sum(active_f==1)+sum(active_k==1))/(N+F);
     
     %define shares by which gov expenditure (discretionary and health) is
     %distributed between firms
     if lockdown==0
     quota_exp_c(:,t)= (RIC./sum(RIC));  % quota di mercato calcolata come share entrate
     quota_health(1:F,t)=(RIC./(sum(RIC)+sum(RIC_k)));
     quota_health(F+1:F+N,t)=(RIC_k./(sum(RIC)+sum(RIC_k)));
     else
     rev=RIC;
     revk=RIC_k;
     rev(deactivate_f)=0;
     revk(active_k==0)=0;
     quota_exp_c(:,t)= (rev./sum(rev));  % quota di mercato calcolata come share entrate
     quota_health(1:F,t)=(rev./(sum(rev)+sum(revk)));
     quota_health(F+1:F+N,t)=(revk./(sum(rev)+sum(revk)));
     end
     
    Oc(dead==1)=-1;             %make sure that dead are not counted as employed or unemployed
    
    if A<=0                                                             %if all the firms are defaulted, then stop simulation 
        %break
    end
    totK(t)=sum(K);
    
    %firms define expected demands and bid prices
    %note that given search and matching in the demand, the demand
    %perceived by the firms is overestimating actual demand. 
    stock = Y_prev-Yd;
    De(active_f==1)=Y_prev(active_f==1)-stock(active_f==1)*q_adj;
    De(active_f==0)=0;
    De(active_f'==1 & De<alpha)=alpha;
    P(active_f'==1 & stock<=0 & P<price(t))=P(active_f'==1 & stock<=0 & P<price(t)).*(1+shock_p(t,(active_f'==1 & stock<=0 & P<price(t)))*p_adj);
    P(active_f'==1 & stock>0 & P>price(t))=P(active_f'==1 & stock>0 & P>price(t)).*(1-shock_p(t,(active_f'==1 & stock>0 & P>price(t)))*p_adj);   
    
    %CAPITAL PRODUCTION DECISION
    %%capital good is durable, inventories depreciates
    inv_dep=inventory_depreciation*Y_k;
    inv_dep=inv_dep.*P_k;
    inventory_k=(1-inventory_depreciation)*Y_k;
    invent_start=inventory_k.*P_k;
    
    %for K-firms everything is as above for C-firms!
    stock_k = Y_prev_k-Y_kd; %virtual variation of inventories
    
    De_k(active_k==1)=Y_prev_k(active_k==1)-stock_k(active_k==1)*q_adj-inventory_k(active_k==1);
    De_k(active_k==0)=0;
    De_k(active_k'==1 & De_k<alpha)=alpha;
    P_k(active_k'==1 & stock_k<=0 & P_k<price_k(t))=P_k(active_k'==1 & stock_k<=0 & P_k<price_k(t)).*(1+shock_pk(t,(active_k'==1 & stock_k<=0 & P_k<price_k(t)))*p_adj);
    P_k(active_k'==1 & stock_k>0 & P_k>price_k(t))=P_k(active_k'==1 & stock_k>0 & P_k>price_k(t)).*(1-shock_pk(t,(active_k'==1 & stock_k>0 & P_k>price_k(t)))*p_adj);
   
    
    %% investments
    %check who is allowed to invest and formulate investment demand
    prob=prob_k(t,:);
    K_dem=zeros(1,F);
    K_des=zeros(1,F);
    K_des(prob<Iprob) = barYK(prob<Iprob)/barX;
    %only used capital depreciates
    depreciation = barX*K_des*eta*1/(Iprob);
    K_dem(prob<Iprob)= K_des(prob<Iprob) - K(prob<Iprob)+depreciation(prob<Iprob);

    K_dem(K_dem<0)=0;
    
    %labour requirement (consumption good)
    Ld=min(ceil(De./alpha), ceil(K.*k/alpha));    %%demand for labour is the minimum between the actual production needs and the employable workers given the capital                                                 %labour demand approximated by excess
    wages=wb*Ld;   
   
    %the consumption good firms looks at the prices on the capital market
    %have an idea of the investment costs
    p_k=price_k(t);
    
   
    %labour requirement (capital good)
    Ld_k=ceil(De_k./alpha);
    wages_k = wb*Ld_k;
	
	if CIGS==1 && t>=lock_start && lock_start>0 && t<=lock_start+11
    labour_transfer(:)=0;
    labour_transfer_k(:)=0;
       for i=1:F
           if Leff(i)>Ld(i)
              labour_transfer(i)=(Leff(i)-Ld(i))*wb;
              liquidity(i)=liquidity(i)+labour_transfer(i);
              A(i)=A(i)+labour_transfer(i);
           end
       end
       for j=1:N
           if Leff_k(j)>Ld_k(j)
              labour_transfer_k(j)=(Leff_k(j)-Ld_k(j))*wb;
              liquidity_k(j)=liquidity_k(j)+labour_transfer_k(j);
              A_k(j)=A_k(j)+labour_transfer_k(j);
           end
       end
    elseif CIGS==1 && t>lock_start+11 && lock_start>0
        labour_transfer(:)=0;
        labour_transfer_k(:)=0;
    end





%% CREDIT MARKET OPENS
             
    
    Ftot(:)=0;
    Ftot_k(:)=0;
   
        
    %%THE FIRMS ASK FOR LIQUIDITY TO BUY THE CAPITAL AT A PRICE EQUAL TO
    %%THE ACTUAL AVERAGE CAPITAL PRICE
    
    %compute financial gap
    B=wages+K_dem*p_k-liquidity;                                        %financial gap          
    B_k =wages_k-liquidity_k;                                           %financial gap of capital producers
    B(B<0)=0;
    B_k(B_k<0)=0;
    
    %compute leverage and use it to define bankruptcy probability and
    %interest rates
    lev=(deb+B)./(A+B+deb);                                             %leverage
    lev_k=(deb_k+B_k)./(A_k+B_k+deb_k);
    
	loan_applications=find(B>0 & bankrupt_f==0);                                        %only firms with positive credit demand will apply for an additional loan
    loan_applications_k=find(B_k>0 & bankrupt_k==0); 
    %evaluate bankruptcy probability and expected survival
    pr(:,1)= exp(b(1)+b(2)*lev)./(1+exp(b(1)+b(2)*lev));                           %banks evaluate bankruptcy probability of each firm
    pr_k(:,1)=exp(b_k(1)+b_k(2)*lev_k)./(1+exp(b_k(1)+b_k(2)*lev_k));                                                                                              
                                                           
    Xi=(1-(1-theta).^(1+1./pr))./(1-(1-theta));                         
    Xi_k=(1-(1-theta).^(1+1./pr_k))./(1-(1-theta));
    %proposed rate depends on the estimated bankruptcy probability
    proposed_rate=mu*((1+r_f/theta)./Xi - theta)';
    proposed_rate_k=mu*((1+r_f/theta)./Xi_k - theta)';
    
    

    %for each firm the bank computes the maximum loan and gives loans up to the maximum amount 
    for i=loan_applications
        credit=B(i);
        %the bank gives a maximum credit depending
        %on  maximum expected loss
        
        maxL = (phi*totE(t)-pr(i)*deb(i))/pr(i);
        maxL=max(0,maxL); %maxL never negative        
        credit=min(credit,maxL); %%credit given to firm i           

        
        deb_0=deb(i);
        deb(i)=deb(i)+credit;                                           %update firm's debt stock
        loans=loans+credit;                                             %update bank's credit stock
        Ftot(i)=credit;                                                 %record flow of new credit for firm i
        liquidity(i)=liquidity(i)+Ftot(i);
        %compute new average interest rate
        if deb(i)>0
            interest_r(i)=(deb_0*interest_r(i)+proposed_rate(i)*credit)/deb(i);
        end
    end       
    %weighted average interest rate
     average_interest_rate(t)=sum(interest_r.*deb)/sum(deb);

     
    %mutatis mutandis for capital firms
    for i=loan_applications_k
        credit=B_k(i);
        maxL = (phi*totE(t)-pr_k(i)*deb_k(i))/pr_k(i);
        maxL=max(0,maxL);
        credit=min(credit,maxL);
         
        deb_0=deb_k(i);
        deb_k(i)=deb_k(i)+credit;                                       %update firm's debt stock
        loans=loans+credit;                                             %update bank's credit stock
        Ftot_k(i)=credit;                                               %record flow of new credit for firm i
        liquidity_k(i)=liquidity_k(i)+Ftot_k(i);
        if deb_k(i)>0
            interest_r_k(i)=(deb_0*interest_r_k(i)+proposed_rate_k(i)*credit)/deb_k(i);
        end
    end 
    average_interest_rate_k(t)=sum(interest_r_k.*deb_k)/sum(deb_k);

    %%CREDIT MARKET CLOSES    
    
    if liquidity_transfers1==1 && t>=lock_start && lock_start>0 && t<=lock_start+11
		liquidity_transfer(:)=0;
		liquidity_transfer_k(:)=0;
		B=wages+K_dem*p_k-liquidity;                                        %financial gap          
		B_k =wages_k-liquidity_k;                                           %financial gap of capital producers
		gap=find(B>0);
		gap_k=find(B_k>0);
		for i=gap
			liquidity_transfer(i)=liquidity_transfer(i)+B(i);
			liquidity(i)=liquidity(i)+liquidity_transfer(i);
			A(i)=liquidity(i)+capital_value(i)-deb(i)-deb_g(i);
		end
		for i=gap_k
			liquidity_transfer_k(i)=liquidity_transfer_k(i)+B_k(i);
			liquidity_k(i)=liquidity_k(i)+liquidity_transfer_k(i);
			A_k(i)=liquidity_k(i)+Y_k(i)*P_k(i)-deb_k(i)-deb_g_k(i);
		end
    elseif liquidity_transfers1==1 && t>lock_start+11 && lock_start>0
		liquidity_transfer(:)=0;
		liquidity_transfer_k(:)=0;
	end
	
	
    if liquidity_transfers2==1 && t>=lock_start && lock_start>0 && t<=lock_start+1
		liquidity_transfer(:)=0;
		liquidity_transfer_k(:)=0;
		liquidity(bankrupt_f==0)=liquidity(bankrupt_f==0)+10;
		liquidity_k(bankrupt_k==0)=liquidity_k(bankrupt_k==0)+10;
		liquidity_transfer_k(bankrupt_k==0)=liquidity_transfer_k(bankrupt_k==0)+10;
		liquidity_transfer(bankrupt_f==0)=liquidity_transfer(bankrupt_f==0)+10;
		A_k=liquidity_k+Y_k.*P_k-deb_k-deb_g_k;
		A=liquidity+capital_value-deb-deb_g;
    elseif liquidity_transfers2==1 && t>lock_start+1 && lock_start>0
		liquidity_transfer(:)=0;
		liquidity_transfer_k(:)=0;
    end
    
    if liquidity_transfers3==1 && t>=lock_start && lock_start>0 && t<=lock_start+1
    %if liquidity_transfers3==1 && t>=lock_start && lock_start>0 && t<=lock_start+11
		liquidity_transfer(:)=0;
		liquidity_transfer_k(:)=0;
		growth=(De-Y_prev)./Y_prev;
        growth(isnan(growth))=0;
        growth_k=(De_k-Y_prev_k)./Y_prev_k;
        growth_k(isnan(growth_k))=0;
        rel_size=Y_prev./mean(Y_prev);
        rel_size_k=Y_prev_k./mean(Y_prev_k);
        rel_profits=profit_rates./mean(profit_rates);
        rel_profits_k=profit_rates_k./mean(profit_rates_k);
%       B=wages+K_dem*p_k-liquidity;                                        %financial gap          
% 		B_k =wages_k-liquidity_k;                                           %financial gap of capital producers
% 		gap=find(B>0 & growth>0);
% 		gap_k=find(B_k>0 & growth_k>0);
%       gap=find(B>0 & rel_size>1);
%       gap_k=find(B_k>0 & rel_size_k>1);
%       gap=find(B>0 & rel_profits>1);
%       gap_k=find(B_k>0 & rel_profits_k>1);
% 		for i=gap
% 			liquidity_transfer(i)=liquidity_transfer(i)+B(i);
% 			liquidity(i)=liquidity(i)+liquidity_transfer(i);
% 			A(i)=liquidity(i)+capital_value(i)-deb(i)-deb_g(i);
% 		end
% 		for i=gap_k
% 			liquidity_transfer_k(i)=liquidity_transfer_k(i)+B_k(i);
% 			liquidity_k(i)=liquidity_k(i)+liquidity_transfer_k(i);
% 			A_k(i)=liquidity_k(i)+Y_k(i)*P_k(i)-deb_k(i)-deb_g_k(i);
% 		end

        liquidity(bankrupt_f==0 & growth>0)=liquidity(bankrupt_f==0 & growth>0)+10;
		liquidity_k(bankrupt_k==0& growth_k>0)=liquidity_k(bankrupt_k==0& growth_k>0)+10;
		liquidity_transfer_k(bankrupt_k==0& growth_k>0)=liquidity_transfer_k(bankrupt_k==0& growth_k>0)+10;
		liquidity_transfer(bankrupt_f==0& growth>0)=liquidity_transfer(bankrupt_f==0& growth>0)+10;
		
%       liquidity(bankrupt_f==0 & rel_size>1)=liquidity(bankrupt_f==0 & rel_size>1)+10;
% 		liquidity_k(bankrupt_k==0& rel_size_k>1)=liquidity_k(bankrupt_k==0& rel_size_k>1)+10;
% 		liquidity_transfer_k(bankrupt_k==0& rel_size_k>1)=liquidity_transfer_k(bankrupt_k==0& rel_size_k>1)+10;
% 		liquidity_transfer(bankrupt_f==0& rel_size>1)=liquidity_transfer(bankrupt_f==0& rel_size>1)+10;
% 		
%         liquidity(bankrupt_f==0 & rel_profits>1)=liquidity(bankrupt_f==0 & rel_profits>1)+10;
% 		liquidity_k(bankrupt_k==0& rel_profits_k>1)=liquidity_k(bankrupt_k==0& rel_profits_k>1)+10;
% 		liquidity_transfer_k(bankrupt_k==0& rel_profits_k>1)=liquidity_transfer_k(bankrupt_k==0& rel_profits_k>1)+10;
% 		liquidity_transfer(bankrupt_f==0& rel_profits>1)=liquidity_transfer(bankrupt_f==0& rel_profits>1)+10;
% 		
        
        A_k=liquidity_k+Y_k.*P_k-deb_k-deb_g_k;
		A=liquidity+capital_value-deb-deb_g;
    elseif liquidity_transfers3==1 && t>lock_start+1 && lock_start>0
    %elseif liquidity_transfers3==1 && t>lock_start+11 && lock_start>0    
		liquidity_transfer(:)=0;
		liquidity_transfer_k(:)=0;
	end
	
	if credit_guarantees==1 && t>=lock_start && lock_start>0 && t<=lock_start+11
		defaults_gov=0;
		B=wages+K_dem*p_k-liquidity;                                        %financial gap          
		B_k =wages_k-liquidity_k;                                           %financial gap of capital producers
		gap=find(B>0);
		gap_k=find(B_k>0);
    
		for i=gap
			credit_g=B(i);       
			deb_g(i)=deb_g(i)+credit_g;                                           %update firm's debt stock
			loans_g=loans_g+credit_g;                                             %update bank's credit stock
			liquidity(i)=liquidity(i)+credit_g;
		end
		for i=gap_k
			credit_g_k=B_k(i);       
			deb_g_k(i)=deb_g_k(i)+credit_g_k;                                           %update firm's debt stock
			loans_g=loans_g+credit_g_k;                                             %update bank's credit stock
			liquidity_k(i)=liquidity_k(i)+credit_g_k;
		end
    elseif credit_guarantees==1 && t>lock_start+11 && lock_start>0
		defaults_gov=0;
	end

    %% JOB MARKET OPENS
     
    %determine desired labour and vacancies given available liquidity
    Ld_k=min(Ld_k, (liquidity_k)/wb);
    Ld_k=floor(Ld_k);
    Ld_k(Ld_k<1)=1;
    vacancies_k=Ld_k-Leff_k;
    surplus_k=find(vacancies_k<0);                                          %firms with too many workers
    
    %%CONSUPTION GOOD

    %%re-define labour demand given available liquidity
    Ld=min(Ld,(liquidity)/wb);                                     %demand (stock)     
    Ld=floor(Ld);
    Ld(Ld<1)=1; %%since there is the capital the firms can have positive equity and negative liquidity, in this latter case Ld would be negative, which is impossible
    vacancies=Ld-Leff;                                                  %definitive labour demand (flow)    
    
    %%JOB MARKET OPENS
    surplus=find(vacancies<0);                                          %firms with too many workers
    
    if CIGS==1 && t>=lock_start && lock_start>0 && t<=lock_start+11
    else
    for i=surplus_k
        workforce_k=find(Oc==F+i);
        f_k=randperm(length(workforce_k));
        f_k=f_k(1:-vacancies_k(i));                                     %take randomly "-vacancies(i)" workers and fire them
        fired_k=workforce_k(f_k);  
        Oc(fired_k)=0;
        w(fired_k)=0;
        Leff_k(i)=Ld_k(i);                                              %update no. of workers
    end
	end

 %firms with excess workforce fire
	if CIGS==1 && t>=lock_start && lock_start>0 && t<=lock_start+11
    else
    for i=surplus
        workforce=find(Oc==i);
        f=randperm(length(workforce));
        f=f(1:-vacancies(i));                                           %take randomly "-vacancies(i)" workers and fire them
        fired=workforce(f);  
        Oc(fired)=0;
        w(fired)=0;
        Leff(i)=Ld(i);                                                  %update no. of workers
    end
    end
    
%% UNEMPLOYED WORKERS LOOK FOR A JOB
    
    %only active unemployed participate in labour market
    unemployed=find(Oc==0 & active==1);
    vec = randperm(length(unemployed));
    active_cfirms=find(active_f==1);
    active_kfirms=find(active_k==1);
    active_firms=[active_cfirms; F+active_kfirms];
    for un=vec      
        j=unemployed(un);                                               %randomly pick an unemployed worker
        Z_e = randperm(length(active_firms),z_e);
        flag=1;        
        while (Oc(j)==0 && flag<=z_e)                              %continue searching until you are unemployed and didn't search at all available firms
            f=active_firms(Z_e(flag));                                        %pick first best firm; with flag increasing, pick the second best, the third...   
            if f>F %selected firm is a capital firm
                if vacancies_k(f-F)>0                                   %if the selected firm has an open vacancy, take the job
                Oc(j)=f;                                                 %update employed status
                w(j)=wb;                                                 %salary   
                Leff_k(f-F)=Leff_k(f-F)+1;                               %firm's workforce   
                vacancies_k(f-F)=vacancies_k(f-F)-1;                     %firm's vacancies   
                end
            else %selected firm is a consuption firm
                if vacancies(f)>0                                     %if the selected firm has an open vacancy, take the job
                Oc(j)=f;
                w(j)=wb;
                Leff(f)=Leff(f)+1;
                vacancies(f)=vacancies(f)-1;
                end                
            end
            flag=flag+1;                                                %increase counter
        end
        
    end  	
    
    %%JOB MARKET CLOSES
    
    if lockdown==1 && t==lock_start+1
        active_f(deactivate_f)=0;
        locked_f(deactivate_f)=1;
        constraint_k(:)=lock_constraint;         
        constraint_f(:)=lock_constraint;
    end
    
    if epidemic==1 && t>=577 && t<epidemic_end
    dieprob=max(0.0075/3,dieprob*(1-post_lock_adjustment));%post_lock_adjustment
    h_3=max(0.0375/3,h_3*(1-post_lock_adjustment));
    work_cons=[];
    for j=1:(F)
        if active_f(j)==1
        colleagues=find(Oc==j);
        if homeoffice(j)==1
              colleagues=randsample(colleagues,round(reduce_workcons_lockdown*length(colleagues)));
        end
        if length(colleagues)>1
           work_cons=[work_cons; nchoosek(colleagues,2)];
        end
        end
    end
    for j=1:(N)
        if active_k(j)==1
        colleagues=find(Oc==F+j);
        if homeoffice_k(j)==1
              colleagues=randsample(colleagues,round(reduce_workcons_lockdown*length(colleagues)));
        end
        if length(colleagues)>1
            work_cons=[work_cons; nchoosek(colleagues,2)];
        end
        end
    end
    work_weights=2*ones(length(work_cons),1);
    work_cons=[work_cons work_weights];
    end
    
    %% production
    %produce capital
    Y_k=min(De_k,constraint_k.*Leff_k*alpha);                              
    Y_prev_k=Y_k;
    Y_k = Y_k+inventory_k;                                   %Y_k is increased by the inventories (capital good is durable)   
    
    %produce consuption
    Yp = min(constraint_f.*Leff*alpha, K*k);                               %production frontier given available resources (Leontieff)
    Y=min(De,Yp);                                            %actual production
    Y_prev=Y;
    
    %here we compute the average production cost to impose the minimum
    %price
    %compute interests to be paid at the end of the period
    interests=interest_r.*deb+r_f.*deb_g; 
    interests_k=interest_r_k.*deb_k+r_f.*deb_g_k;
    %total wages paid by firms
    wages=wb*Leff;
    wages_k=wb*Leff_k;
    liquidity=liquidity-wages;
    liquidity_k=liquidity_k-wages_k;
    
    %minimum prices
    
    Pl = 1.01*(wages-labour_transfer+interests+Y/k*eta*p_k+theta*deb+theta*deb_g)./Y_prev;
    Pl(isinf(Pl)) = 0;
    P(P<Pl) = Pl(P<Pl);

    
    %minimum prices capital
    Pl_k = 1.01*(wages_k-labour_transfer_k+interests_k+theta*deb_k+theta*deb_g_k)./(Y_prev_k);
    Pl_k(isinf(Pl_k)) = 0;
    P_k(P_k<Pl_k) = Pl_k(P_k<Pl_k);
    

  
    %%CAPITAL GOODS MARKET OPENS
    %%CONSUMPTION FIRMS BUY THE CAPITAL THEY NEED
    % The capital is bought now and will be used for production in the next
    % period. The capital market is a search and matching as the consumption
    % market

    capital_budget=liquidity;                              %amount of liquidity available to buy capital goods
    capital_demanders=find(K_dem>0);                                    
    investment(:)=0;
    value_investments(:)=0;
    Q_k(:)=0;
    Y_kd=zeros(1,N);  %%record demand
    
    %Government buys goods from K-sector for healthcare; share based on
    %previous revenue (see below)
    health_demand_k=zeros(1,N);
    hdk=quota_health(F+1:F+N,t)*healthcare(t);
    while sum(hdk)>0 && sum(Y_k)>0
    for j=1:N
        Y_kd(j) = Y_kd(j) + hdk(j);
        if Y_k(j)>= hdk(j)
           Y_k(j)=Y_k(j)-hdk(j);   
           Q_k(j)= Q_k(j)+hdk(j);
           health_demand_k(j)=health_demand_k(j)+hdk(j);
           hdk(j)=0;
        else 
           Q_k(j)=Q_k(j)+Y_k(j);
           health_demand_k(j)=health_demand_k(j)+Y_k(j);              
           hdk(j)=hdk(j)-Y_k(j);
           Y_k(j)= 0;
        end
    end
    hresid=sum(hdk);
    quota_new=Y_k/sum(Y_k);
    hdk=quota_new*hresid;
    end
    Healthcare_k= sum(health_demand_k.*P_k);
    healthcare_k=sum(health_demand_k);
    
    
    %the capital market works exactly as the consumption market (below).
    %Search and matching!
    partners_k(:,2)=randi(length(active_kfirms),F,1);	
    partners_k(:,2)=active_kfirms(partners_k(:,2));
    for firm=randperm(length(capital_demanders))
        j=capital_demanders(firm);
        %randomly chosen firms
        Z = partners_k(j,:);	
        while Z(1)==Z(2)	
            Z(2)=randi(randi(length(active_kfirms)));	
            Z(2)=active_kfirms(Z(2));	
        end	
        if rand<(1-exp(price_sensitivity*(P_k(Z(2))-P_k(Z(1)))/P_k(Z(1)))) 	
           partners_k(j,1)=Z(2);	
           partners_k(j,2)=Z(1);	
           Z=partners_k(j,:);	
        end
        flag=1;
        while (capital_budget(j)>0 && K_dem(j)>0 && flag<=z_k) 
            best=Z(flag);
            Y_kd(best)=Y_kd(best)+min(capital_budget(j)/P_k(best), K_dem(j));
                if Y_k(best)>0                                                %buy if 'best' firm has still positive stocks 
                    pk=P_k(best);
                    budget=min(capital_budget(j), K_dem(j)*pk);
                    if Y_k(best) > budget/pk           
                        Y_k(best)=Y_k(best)-budget/pk;                   %reduce stocks
                        Q_k(best)=Q_k(best)+budget/pk;                   %update sales
                        K_dem(j)=K_dem(j)-budget/pk;
                        capital_budget(j) = capital_budget(j)-budget;
                        liquidity(j)=liquidity(j)-budget;
                        investment(j)=investment(j)+budget/pk;
                        value_investments(j)= value_investments(j)+budget;                                 %j spends all its budget  
                    elseif  Y_k(best) <= budget/pk  
                        K_dem(j)=K_dem(j)-Y_k(best);
                        capital_budget(j)=capital_budget(j) - Y_k(best)*pk;
                        liquidity(j)=liquidity(j)- Y_k(best)*pk;
                        Q_k(best)=Q_k(best)+Y_k(best);    
                        investment(j)=investment(j)+Y_k(best);
                        value_investments(j)= value_investments(j)+Y_k(best)*pk;
                        Y_k(best)=0;    
                    end    
                end
                flag=flag+1;                                                %increase counter
        end
    end
    


    %%CAPITAL GOOD MARKET CLOSES
 
    %%CONSUMPTION GOOD MARKET OPENS
    w(w>0) = wb;
    w(Oc==0)=0;
    w(dead==1)=0;
    PA(dead==1)=0;
    %taxes on wages
    wn=w*(1-tax_rate);  
    TA(t) = TA(t)+sum(w*tax_rate);
    %workers receive wages
    interest_workers=r_d*PA(1:W);
    interest_deposits=sum(r_d*PA(1:W));
    PA(1:W)=PA(1:W)+r_d*PA(1:W);
    PA(1:W)=PA(1:W)+wn;
   

    unemployment_subsidy = unemployment_subsidy_init;    
    
    %pay unemployement subsidy to unemployed, sickpay to sick and inactive income to
    %inactive
    oPA=sum(PA(Oc==0 & active==1 & dead==0));
    nPA=sum(PA(Oc==0 & active==1 & dead==0)+unemployment_subsidy*wb*(1-tax_rate));
    dPA=nPA-oPA;
    G(t)=dPA;
    PA(Oc==0 & active==1 & dead==0)=PA(Oc==0 & active==1 & dead==0)+unemployment_subsidy*wb*(1-tax_rate);
    oPA=sum(PA(active==0 & sickpay==0 & dead==0));
    nPA=sum(PA(active==0 & sickpay==0 & dead==0)+inactive_sub*unemployment_subsidy*wb*(1-tax_rate));
    dPA=nPA-oPA;
    G(t)=G(t)+dPA;
    PA(active==0 & sickpay==0 & dead==0)=PA(active==0 & sickpay==0 & dead==0)+inactive_sub*unemployment_subsidy*wb*(1-tax_rate);
    oPA=sum(PA(active==0 & sickpay==1 & dead==0));
    nPA=sum(PA(active==0 & sickpay==1 & dead==0)+unemployment_subsidy*wb*(1-tax_rate));      %assume sickpay=unemployment benefit
    PA(active==0 & sickpay==1 & dead==0)=PA(active==0 & sickpay==1 & dead==0)+unemployment_subsidy*wb*(1-tax_rate); 
    dPA=nPA-oPA;
    G(t)=G(t)+dPA;
	if income_support==1 && t>=lock_start && lock_start>0 && t<=lock_start+5
       PA(Oc>=0)=PA(Oc>=0)+0.5*unemployment_subsidy*wb*(1-tax_rate);
       G(t)=G(t)+0.5*unemployment_subsidy*wb*(1-tax_rate)*length(PA(Oc>=0));
    end
    
    workers_income = wn;
    workers_income(Oc==0 & active==1) = unemployment_subsidy*wb*(1-tax_rate);
    workers_income(active==0 & sickpay==0) = inactive_sub*unemployment_subsidy*wb*(1-tax_rate);
    workers_income(active==0 & sickpay==1) = unemployment_subsidy*wb*(1-tax_rate);
    if income_support==1 && t>=lock_start && lock_start>0 && t<=lock_start+5
    workers_income = workers_income+0.5*unemployment_subsidy*wb*(1-tax_rate);
    end
	workers_income(dead==1)=0;
    workers_income=workers_income+interest_workers;
    
    %consumers compute their consumption budgets
    income =  ([workers_income,dividends_income,dividends_income_k]')./price(t); 

    permanent_income = permanent_income*xi + (1-xi)*income;

    target = 1*permanent_income' + chi*PA./price(t) ; %0.05
    target_b=F_b/F*target.*CPI_l(t-1)/CPI_b(t-1);
    target_l=target-target_b;
    %compute a budget for basic goods as fraction of total cons. budget.
    %Here this is also a function of relative prices CPI_l and CPI_b
    %lexicographic ordering (only spend on luxury if sufficient money
    %available)
    cons_budget_b=min(PA,target_b.*CPI_b(t-1));
    cons_budget_l=max(0,min(PA-cons_budget_b,target_l.*CPI_l(t-1)));
    if t>=577 && epidemic==1 && t<epidemic_end
    for i=1:W
    if distancing(i)==1
	c_shock_b(i)=min(1.5,1+shock_b*sum(detected==1 & infected==1));%*exp(-0.1*(t-649));
    c_shock_l(i)=max(0.5,1-shock_l*sum(detected==1 & infected==1));%*exp(-0.1*(t-649));
    cons_budget_b(i)=min(PA(i),c_shock_b(i)*target_b(i)*CPI_b(t-1));
    cons_budget_l(i)=max(0,min(PA(i)-cons_budget_b(i),target_l(i)*CPI_l(t-1)));   
    cons_budget_l(i)=c_shock_l(i)*cons_budget_l(i);
    end
    end
    end
    
    PA=PA-cons_budget_b-cons_budget_l;        
    
    DC(t) = sum(cons_budget_b)/CPI_b(t-1)+sum(cons_budget_l)/CPI_l(t-1);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%INTRO SPESA PUBBLICA
    Q(:)=0;
    Yd=zeros(1,F);
    
    %government buys from all C-firms (including "luxury" --> can view this
    %as consumption by doctors??)
    health_demand_f=zeros(1,F);
    hd=quota_health(1:F,t)*healthcare(t);
    while sum(hd)>0 && sum(Y)>0
    for i=1:F
        Yd(i) = Yd(i) + hd(i);
        if Y(i)>= hd(i)
           Y(i)=Y(i)-hd(i);   
           Q(i)= Q(i)+hd(i);
           health_demand_f(i)=health_demand_f(i)+hd(i);
           hd(i)=0;
        else 
           Q(i)=Q(i)+Y(i);
           health_demand_f(i)=health_demand_f(i)+Y(i);
           hd(i)=hd(i)-Y(i);
           Y(i)= 0;      
        end
    end
    hresid=sum(hd);
    quota_new=Y/sum(Y);
    hd=quota_new*hresid;
    end
    Healthcare_f= sum(health_demand_f.*P);
    healthcare_f=sum(health_demand_f);
    
    Healthcare(t)=Healthcare_k+Healthcare_f;
    healthcare(t)=healthcare_k+healthcare_f;
   %search and matching starts
    
    C = zeros(1,W+F+N);
    
    %same as with capital goods market
    vec = randperm(W+F+N);
    active_bfirms=find(active_f(1:F_b)==1);	
    active_lfirms=F_b+find(active_f(F_b+1:F)==1);	
    partners_b(:,2)=randi(length(active_bfirms),W+F+N,1);	
    partners_b(:,2)=active_bfirms(partners_b(:,2));	
    partners_l(:,2)=randi(length(active_lfirms),W+F+N,1);	
    partners_l(:,2)=active_lfirms(partners_l(:,2));
    
    for j=vec     
        if j<=W && dead(j)==1
        else
        Z = partners_b(j,:);	
        while Z(1)==Z(2)	
            Z(2)=randi(length(active_bfirms));	
            Z(2)=active_bfirms(Z(2));	
        end	
        if rand<(1-exp(price_sensitivity*(P(Z(2))-P(Z(1)))/P(Z(1)))) 	
           partners_b(j,1)=Z(2);	
           partners_b(j,2)=Z(1);	
           Z=partners_b(j,:);
        end
        flag=1;
        %first shop with basic firms, luxury below
        while (cons_budget_b(j)>0 && flag<=z_c)                              %continue buying till budget is positive and there are firms available
            best=Z(flag);                                       %pick first best firm; with flag increasing, pick the second best, the third...   
            Yd(best)=Yd(best)+cons_budget_b(j)/P(best);
            if Y(best)>0                                                %buy if 'best' firm has still positive stocks 
                p=P(best);
                if Y(best) > cons_budget_b(j)/p           
                    Y(best)=Y(best)-cons_budget_b(j)/p;                   %reduce stocks
                    Q(best)=Q(best)+cons_budget_b(j)/p; 
                    C(j) = C(j)+cons_budget_b(j)/p; %update sales
                    consumption(t)=consumption(t)+cons_budget_b(j)/p;     
                    cons_budget_b(j)=0;                                   %j spends all its budget  
                elseif  Y(best) <= cons_budget_b(j)/p  
                    cons_budget_b(j)=cons_budget_b(j)- Y(best)*p;   
                    Q(best)=Q(best)+Y(best);
                    C(j) = C(j)+Y(best); 
                    consumption(t)=consumption(t)+Y(best);
                    Y(best)=0;
                end    
            end
            flag=flag+1;                                                %increase counter
        end 
        Z = partners_l(j,:);	
        while Z(1)==Z(2)	
            Z(2)=randi(length(active_lfirms));	
            Z(2)=active_lfirms(Z(2));	
        end	
        if rand<(1-exp(price_sensitivity*(P(Z(2))-P(Z(1)))/P(Z(1)))) 	
           partners_l(j,1)=Z(2);	
           partners_l(j,2)=Z(1);	
           Z=partners_l(j,:);	
        end                                           %sort prices in ascending order
        flag=1;
        while (cons_budget_l(j)>0 && flag<=z_c)                              %continue buying till budget is positive and there are firms available
            best=Z(flag);                                        %pick first best firm; with flag increasing, pick the second best, the third...   
            Yd(best)=Yd(best)+cons_budget_l(j)/P(best);
            if Y(best)>0                                                %buy if 'best' firm has still positive stocks 
                p=P(best);
                if Y(best) > cons_budget_l(j)/p           
                    Y(best)=Y(best)-cons_budget_l(j)/p;                   %reduce stocks
                    Q(best)=Q(best)+cons_budget_l(j)/p; 
                    C(j) = C(j)+cons_budget_l(j)/p; %update sales
                    consumption(t)=consumption(t)+cons_budget_l(j)/p;     
                    cons_budget_l(j)=0;                                   %j spends all its budget  
                elseif  Y(best) <= cons_budget_l(j)/p  
                    cons_budget_l(j)=cons_budget_l(j)- Y(best)*p;   
                    Q(best)=Q(best)+Y(best);
                    C(j) = C(j)+Y(best); 
                    consumption(t)=consumption(t)+Y(best);
                    Y(best)=0;
                end    
            end
            flag=flag+1;                                                %increase counter
        end
       end
    end    
    
    %search and matching ends
    %unvoluntary savings are added to the vuluntary ones
    unsatisfiedDemand(t) = sum(cons_budget_b)+sum(cons_budget_l);            
    PA=PA+cons_budget_b+cons_budget_l;
    
    % dopo vende al settore pubblico, se ha prodotto abbastanza
    %this is discretionary gov. expenditure
    pub_exp_c(t)= 1*EXP(t);
    public_dem_c(:,t)=quota_exp_c(:,t).*(pub_exp_c(t)./P');  % dove pub_exp è quella EROGABILE
    for i=1:F
        Yd(i) = Yd(i) + public_dem_c(i,t);
        if Y(i)>= public_dem_c(i,t)
           Y(i)=Y(i)-public_dem_c(i,t);     % riduco lo stock di offerta
           Q(i)= Q(i)+public_dem_c(i,t);    % ho venduto l'ammontare "public dem"
        else 
           Q(i)=Q(i)+Y(i);
           public_dem_c(i,t)=Y(i);              % così riflette il fatto che l'impresa potrebbe non avere
           Y(i)= 0;      
        end                                     % abbastanza offerta per soddisfare la domanda pubblica
    end
    exp_c(:,t)= public_dem_c(:,t).*P';     % per aggiornamento EXP effettivamente erogate complessivam
    pub_exp_c(t)=sum(exp_c(:,t));
    EXP(t)=pub_exp_c(t); % qui sempre per la storia che se faccio come lo volevo fare mi si impalla
    
%%CONSUMPTION GOOD MARKET CLOSES
    %Capital price index
    if sum(Q_k)>0 
        share_k=(Q_k)/sum(Q_k);
    else
        share_k=ones(1,N)/N;
    end
    
    price_k(t+1)=P_k*share_k';

%%ACCOUNTING

    %capital capacity update
    barYK = delta*barYK + (1-delta)*Y_prev/k;
    %barYK = delta*barYK + (1-delta)*Yd/k;
    dep = eta*Y_prev/k; 
    
    %capital reduced by depreciation
    K = K - dep;
    %capital increased by bought capital 
    K = K + investment;
    I(t) = sum(investment);

    %update capital value in the book
    depreciation_value =  dep.*capital_value./(K+dep-investment);
    capital_value = capital_value - depreciation_value + value_investments;
    
    %firm revenues
    RIC=P.*Q;   
    %firm uodate their liquidity, pay wages, interests and installments
    %consumption firm profits
    pi = RIC-wages-interests - depreciation_value;
    liquidity=liquidity+RIC-interests-theta*deb-theta*deb_g; %%investments already paid
    %loans are updated (the installments have been paid)
    loans=loans-sum(deb)*theta;
    loans_g=loans_g-sum(deb_g)*theta;
    deb=(1-theta)*deb;
    deb_g=(1-theta)*deb_g;
    
    %equity law of motion!
    A = A+pi;
    
    %Capital producer accounting
    RIC_k = P_k.*Q_k;
    %update liquidity
    invent_end=Y_k.*P_k;
    pi_k=RIC_k+(invent_end-invent_start)-wages_k-interests_k-inv_dep;
    liquidity_k = liquidity_k + RIC_k-interests_k-theta*deb_k-theta*deb_g_k;
    %update loand
    loans=loans-sum(deb_k)*theta;
    loans_g=loans_g-sum(deb_g_k)*theta;
    deb_k=(1-theta)*deb_k;
    deb_g_k=(1-theta)*deb_g_k;
    
    %pre-dividends average profit rate
    profit_rate_b=mean(pi(active_f(1:F_b)==1))*12/mean(A(active_f(1:F_b)==1));	
    profit_rate_l=mean(pi(active_f(F_b+1:F)==1))*12/mean(A(active_f(F_b+1:F)==1));	
    profit_rate_k=mean(pi_k(active_k==1))*12/mean(A_k(active_k==1));	
    profit_rates=pi*12./A;
    profit_rates_k=pi_k*12./A_k;
    
    
    %dividends
    pospi=find(pi>0 & liquidity>0);                                                   %pick firms with positive profits
    dividends_income(:)=0;
	
	if t>577
	div_g(:)=0;
    div_g_k(:)=0;
	end
	
    for i=pospi
        TA(t)=TA(t)+pi(i)*tax_rate_f;
        tax_pr=pi(i)*tax_rate_f;
        pi(i)=(1-tax_rate_f)*pi(i);
        liquidity(i)=liquidity(i)-tax_pr;
        A(i)=A(i)-tax_pr;
        di=max(0,min(div*pi(i),liquidity(i)));  %%dividends                                              %compute dividends paid by firm i
        div_g(i)=gov_ownership(i)/(gov_ownership(i)+A_init(i))*di;
        divi=di-div_g(i); %dividends after taxes
        PA(W+i)=PA(W+i)+divi;                                          %dividends paid to firm i owner
        dividends_income(i) = divi;
        liquidity(i)=liquidity(i)-di;
        A(i)=A(i)-di;
        dividends(t)=dividends(t)+di; %lordi
        pi(i)=pi(i)-di;
    end
    dividends_income=dividends_income+r_d*PA(W+1:W+F);
    interest_deposits=interest_deposits+sum(r_d*PA(W+1:W+F));
    PA(W+1:W+F)=PA(W+1:W+F)+r_d*PA(W+1:W+F);
    
    A_control=A_control+pi_k;
    pospi_k=find(pi_k>0 & liquidity_k>0); 
    dividends_income_k(:)=0;%pick firms with positive profits
    for i=pospi_k
        TA(t)=TA(t)+pi_k(i)*tax_rate_k;
        tax_pr_k=pi_k(i)*tax_rate_k;
        pi_k(i)=(1-tax_rate_k)*pi_k(i);
        liquidity_k(i)=liquidity_k(i)-tax_pr_k;
        A_control(i)=A_control(i)-tax_pr_k;
        di=max(0,min((div)*pi_k(i),liquidity_k(i))); 
        div_g_k(i)=gov_ownership_k(i)/(gov_ownership_k(i)+A_init_k(i))*di;
        divi=di-div_g_k(i);
        PA(W+F+i)=PA(W+F+i)+divi;                                          %dividends paid to firm i owner
        dividends_income_k(i)=divi;
        A_control(i)=A_control(i)-di;
        liquidity_k(i)=liquidity_k(i)-di;
        dividends(t)=dividends(t)+di;
        pi_k(i)=pi_k(i)-di;
       
    end
    dividends_income_k=dividends_income_k+r_d*PA(W+F+1:W+F+N);
    interest_deposits=interest_deposits+sum(r_d*PA(W+F+1:W+F+N));
    PA(W+F+1:W+F+N)=PA(W+F+1:W+F+N)+r_d*PA(W+F+1:W+F+N);
   
    A_k=liquidity_k+Y_k.*P_k-deb_k;

    
    %replacement of bankrupted consumption firms 
    piB=0;                                                              %reset bank's profits
       
    
    %% time series (before bankruptcies)
    %inflation rate
    if sum(Q)>0 
        share=(Q)/sum(Q);
    else
        share=ones(1,F)/F;
    end
    RPI=P*share';                                                       %retail price index
    price(t+1)=RPI;
    P_b=P(1:F_b);
    Q_b=Q(1:F_b);
    if sum(Q_b)>0
    CPI_b(t)=sum(P_b.*(Q_b./(sum(Q_b))));
    else
    CPI_b(t)=sum(P_b.*(ones(1,F_b)./F_b));
    end
    P_l=P(F_b+1:F);
    Q_l=Q(F_b+1:F);
    if sum(Q_l)>0
    CPI_l(t)=sum(P_l.*(Q_l./(sum(Q_l))));
    else
    CPI_l(t)=sum(P_l.*(ones(1,F_l)./F_l));   
    end
    %unemployment rate
    Un(t)=(sum(active)-length(Oc(Oc>0)))/sum(active);
    
    Y_nominal_k(t) = sum(Y_prev_k)*price_k(t);
    Y_nominal_c(t) = sum(Y_prev)*price(t);
    Y_nominal_tot(t)= Y_nominal_k(t)+Y_nominal_c(t);
    Y_real(t) = sum(Y_prev)*price(1) + sum(Y_prev_k)*price_k(1);
    
    gdp_deflator(t) = Y_nominal_tot(t)/ Y_real(t);
    consumption(t)=consumption(t)*price(1);
    
    Investment(t)=I(t)*price_k(1)+sum(Y_k-inventory_k)*price_k(1)+price(1)*sum(Y); %total investment is investment plus inventory variation
    Inventories(t)=sum(Y_k-inventory_k)*price_k(1)+price(1)*sum(Y);
    gross_investment_share(t)=Investment(t)/Y_real(t);
    net_investment_share(t)=I(t)*price_k(1)/Y_real(t);
    inventories_share(t)=Inventories(t)/Y_real(t);
    consumption_share(t)=consumption(t)/Y_real(t);
    
    
    Lev_af(t)=sum(deb./(A+deb).*(A/sum(A)));
    Lev_ab(t)=sum(deb(1:F_b)./(A(1:F_b)+deb(1:F_b)).*(A(1:F_b)/sum(A(1:F_b))));
    Lev_al(t)=sum(deb((F_b+1):F_l)./(A((F_b+1):F_l)+deb((F_b+1):F_l)).*(A((F_b+1):F_l)/sum(A((F_b+1):F_l))));
    Lev_ak(t)=sum(deb_k./(A_k+deb_k).*(A_k/sum(A_k)));
    debtot=[deb,deb_k];
    atot=[A,A_k];
    Lev_a(t)=sum(debtot./(atot+debtot).*atot/(sum(atot)));
    
    negcash_k=find(A_k<=0|liquidity_k<0);
    negcash=find(A<=0|liquidity<0);
    
    NetEq = liquidity-deb;
    Y_prevp=Y_prev(NetEq>0);  
    
    if t>=577
       pb_prev=partners_b;
       pl_prev=partners_l;
	   equity_injection(:)=0;
	   equity_injection_k(:)=0;
       defaults_gov=0;
	   meanA=mean(A(A>0));
	   meanA_k=mean(A_k(A_k>0));
    end
    
    %update bankrupted firms!
    for i=negcash                                                       %pick sequentially failed firms
        %if active_f(i)==1
        if bankrupt_f(i)==0
            if liquidity(i)<0 && A(i)>0 && PA(W+i)>=(-liquidity(i))
               liquidity(i)=liquidity(i)+PA(W+i);
               PA(W+i)=0;
               A(i)=liquidity(i)+capital_value(i)-deb(i)-deb_g(i);
            else
                if equity_support==1 && t>=lock_start && lock_start>0 && t<=lock_start+11
                    if -A(i)+meanA>=(-liquidity(i))
                        equity_injection(i)=-A(i)+meanA;
                    else
                        equity_injection(i)=-liquidity(i);    
                    end
                    liquidity(i)=liquidity(i)+equity_injection(i);
                    gov_ownership(i)=gov_ownership(i)+equity_injection(i);
                    A(i)=A(i)+equity_injection(i);
                else
                    defaults(t)=defaults(t)+1;
                    zzz=deb(i);                                                     %take residual debts     
                    xxx=deb_g(i);
                    piB=piB+(liquidity(i)-deb(i)-deb_g(i));                                               %account for bad debts
                    loans=loans-zzz;
                    loans_g=loans_g-xxx;
                    defaults_gov=defaults_gov+xxx;
                    deb(i)=0;
                    deb_g(i)=0;
                    active_f(i)=0;
                    bankrupt_f(i)=1;
                    A(i)=K(i)*price_k(t+1);                                                   %initialize new firm 
                    capital_value(i)=K(i)*price_k(t+1);
                    liquidity(i)=0;
                    %%fire workers
                    workforce=find(Oc==i);                                          %pick all firm i's workers
                    fired=workforce;  
                    Oc(fired)=0;
                    Leff(i)=0;
                    if i<=F_b
                       active_bfirms(active_bfirms==i)=[];
                       partners=find(partners_b(:,1)==i);
                       if ~isempty(partners)
                       partners_b(partners,1)=randi(length(active_bfirms),length(partners),1);
                       partners_b(partners,1)=active_bfirms(partners_b(partners,1));
                       end
                    else
                       active_lfirms(active_lfirms==i)=[];
                       partners=find(partners_l(:,1)==i);
                       if ~isempty(partners)
                       partners_l(partners,1)=randi(length(active_lfirms),length(partners),1);
                       partners_l(partners,1)=active_lfirms(partners_l(partners,1));
                       end
                    end
                end
            end
        end
    end
    
    NetEq_k = A_k;
    if isempty(find(NetEq_k>0, 1))
          warning('all capital firms are bankrupted')
        %break
        %%keep the same variable of last year        
    else
        Y_prevp_k=Y_prev_k(NetEq_k>0);  
    end
                  
    for i=negcash_k                                                       %pick sequentially failed firms
        %if active_k(i)==1
        if bankrupt_k(i)==0    
            if liquidity_k(i)<0 && A_k(i)>0 && PA(W+F+i)>=(-liquidity_k(i))
               liquidity_k(i)=liquidity_k(i)+PA(W+F+i);
               PA(W+F+i)=0;
               A_k(i)=liquidity_k(i)+Y_k(i)*P_k(i)-deb_k(i)-deb_g_k(i);
               A_control(i)=A_k(i);
            else
                if equity_support==1 && t>=lock_start && lock_start>0 && t<=lock_start+11
                    if -A_k(i)+meanA_k>=(-liquidity_k(i))
                        equity_injection_k(i)=-A_k(i)+meanA_k;
                    else
                        equity_injection_k(i)=-liquidity_k(i);    
                    end
                    liquidity_k(i)=liquidity_k(i)+equity_injection_k(i);
                    gov_ownership_k(i)=gov_ownership_k(i)+equity_injection_k(i);
                    A_k(i)=liquidity_k(i)+Y_k(i)*price_k(t)-deb_k(i)-deb_g_k(i);
                else
                    defaults_k(t)=defaults_k(t)+1;
                    zzz=deb_k(i);
                    xxx=deb_g_k(i);
                    piB=piB+(liquidity_k(i)-deb_k(i)-deb_g_k(i));                                               %account for bad debts
                    loans=loans-zzz;
                    loans_g=loans_g-xxx;
                    defaults_gov=defaults_gov+xxx;
                    liquidity_k(i)=0;
                    deb_k(i)=0;
                    deb_g_k(i)=0;
                    active_k(i)=0;
                    bankrupt_k(i)=1;
                    A_k(i)=0;                                                   %initialize new firm 
                    A_control(i)=A_k(i);
                    workforce_k=find(Oc==F+i);                 %pick all firm i's workers
                    fired_k=workforce_k;  
                    Oc(fired_k)=0;
                    Leff_k(i)=0; 
                    stock_k(i)=0;
                    Y_k(i)=0;
                    active_kfirms(active_kfirms==i)=[];
                    partners=find(partners_k(:,1)==i);
                    if ~isempty(partners)
                        partners_k(partners,1)=randi(length(active_kfirms),length(partners),1);
                        partners_k(partners,1)=active_kfirms(partners_k(partners,1));
                    end
                end
            end
        end
    end
    
    bankruptcy_rate(t) = (length(negcash)+length(negcash_k))/(F+N);
    
    %% bank accounting
    INT(t)=bond_interest_rate*bonds(t-1);
    piB=piB+defaults_gov+sum(interests)+sum(interests_k) + INT(t)-interest_deposits+dead_assets;                                             %bank profits  

    if piB>0 && totE(t)>0
        TA(t)=TA(t)+piB*tax_rate_b;
        piB=piB*(1-tax_rate_b);
        dividendsB(t)=div_B*piB;
        piB=(1-div_B)*piB;
        PA(W+1:W+F+N)=PA(W+1:W+F+N)+dividendsB(t)/(F+N);
    else 
        dividendsB(t) = 0;
    end    
    E=E+piB;                                                            %update bank capital

    %%add bank's dividends to income of capitalists
    dividends_income_k = dividends_income_k + dividendsB(t)/(F+N);
    dividends_income = dividends_income + dividendsB(t)/(F+N);

    profitsB(t)=piB;
    
    %Bank bail-in
    
    if E<0
    bank_loss=-E;
    E=(loans+bonds(t-1))*0.06;                                     %re-capitalize bank at minimun CAR
    deposit_share=PA/sum(PA);                                       %compute households' deposit share
    if (bank_loss + E) > sum(PA)
        E=sum(PA)-bank_loss;
        PA=zeros(1,W+F+N); 
    else 
        PA=PA-(bank_loss+E)*deposit_share;                             %depositors bear the loss incurred by the bank proportionally to their deposit share
    end
    %BT(t)=bank_loss+Eb;
    end
    %}
    totE(t+1)=E;
    %PA_t(t,:)=PA;

  
  
 GB(t)= TA(t)+sum(div_g)+sum(div_g_k) - G(t)-defaults_gov- EXP(t)- bond_interest_rate*bonds(t-1)-Healthcare(t)-sum(labour_transfer)-sum(labour_transfer_k)-sum(liquidity_transfer)-sum(liquidity_transfer_k)-sum(equity_injection)-sum(equity_injection_k);
 
 stock_bonds(t) =sum(-GB(1:t));
 bonds(t) = max(0,stock_bonds(t)); 

 debGDP(t)=stock_bonds(t)/Y_nominal_tot(t);
 defGDP(t)=GB(t)/Y_nominal_tot(t);
 
%disp(sum(deb)+sum(deb_k)+sum(deb_g)+sum(deb_g_k)+stock_bonds(t)-sum(PA)-sum(liquidity)-sum(liquidity_k)-E)
%disp(sum(PA)+sum(A)+sum(A_k)+E-stock_bonds(t)-sum(capital_value)-sum(Y_k.*P_k))

%% wage rule

if u_target - Un(t)>0
    wb = wb * (1+ wage_update_up * (u_target - Un(t))); 
else
    wb = wb * (1+ wage_update_down * (u_target - Un(t))); 
end

wages_t(t) = wb;


%% Epidemic sub-model
dead_assets=0;

for tt=((t-1)*4+1):((t-1)*4+4)
    Hdemand(tt)=sum(health_demanded)+sum(health_demanded2);
    
    if tt==2309
       patient0=randperm(W,5);  %choose a few agents (5 in this case) to be the initial cases of the disease
       infected(patient0)=1;
       susceptible(patient0)=0;
       inf_new(tt)=inf_new(tt)+length(patient0);
       epidemic_start=t;
       epidemic_start_tt=tt;
    end
    
	if social_distancing==1 && lockdown_exp==0
		lock_end=2320;
	end
	
    inf_current(tt)=sum(infected);
    detect_current(tt)=sum(detected==1 & infected==1);

    health_supply=max(0,healthcare(t)-sum(health_supplied)-sum(health_supplied2));
    if sum(health_demanded)+sum(health_demanded2)>healthcare(t)
        %disp('over')
    end
    hospital_capacity(tt)=min(1,Hdemand(tt)/tot_healthcare(tt));
    
    if social_distancing==1
        if lockdown==1
           dis_cost=dis_cost_lockdown;
        end
        if lock_ended==1 || lockdown_exp==0 && t>=epidemic_start+3
           dis_cost=min(6,dis_cost+post_lock_adjustment*(tt-lock_end));
           shock_b=max(0.00055/10,0.00055*exp(-post_lock_adjustment*(tt-lock_end)));	
           shock_l=max(0.00165/10,0.00165*exp(-post_lock_adjustment*(tt-lock_end)));
        end
        
        dis_vec=[sum(infected==1 & detected==1)-dis_threshold; sum(distancing==1)/W-sum(distancing==0)/W; -dis_cost];
    
         if epidemic==1 && t>=577 && sum(infected)>0
            for i=1:W
                dis_index(i)=distancing_persistence*dis_index(i)+(1-distancing_persistence)*intensity*dis_vec; 
                dis_prob(i)=1/(1+exp(-dis_index(i)));
                if rand<dis_prob(i) || detected(i)==1
                   distancing(i)=1;
                else
                   distancing(i)=0;
                end
            end
        end
    end
    
    if lockdown==0 && lock_ended==0
       reduce_cons=1;
    elseif lockdown==1
       reduce_cons=reduce_cons_lockdown;
    end
    if lock_ended==1 || lockdown_exp==0 && t>=epidemic_start+3 && social_distancing==1
       reduce_cons=min(1,max(reduce_cons_lockdown,(1-sum(distancing)/W)));
    end
    if lock_ended==1 || lockdown_exp==0 && t>=epidemic_start+3 && social_distancing==1
       detectprob=min(max_detectprob,0.02+detect_adjustment*(detect_new(tt-1)));
    end
    redcons(tt)=reduce_cons;
    dis_share(tt)=sum(distancing)/W;
    
    %% Vaccination
    
    if vaccine==1 && lock_ended==1 && t>=epidemic_start+11
        coverage_rate=min(0.05,coverage_rate+0.001);
        vaccine_capacity=round(coverage_rate*sum(dead==0));
        %vaccine_capacity=W;
        if mutation==0
            wait_list=find(susceptible==1 & vaccinated==0 & dead==0 & no_vax==0 | recovered==1 & vaccinated==0 & dead==0 & no_vax==0 | infected==1 & detected==0 & vaccinated==0 & dead==0 & no_vax==0);
            if length(wait_list)>=vaccine_capacity
                new_vaccine=vaccine_capacity;
            else
                new_vaccine=length(wait_list);
            end
        else
            wait_list=find(susceptible==1 & vaccinated==0 & dead==0 & no_vax==0 | recovered==1 & vaccinated==0 & dead==0 & no_vax==0 | infected==1 & detected==0 & vaccinated==0 & dead==0 & no_vax==0);
            if length(wait_list)>=vaccine_capacity
                new_vaccine=vaccine_capacity;
            else
                new_vaccine=length(wait_list);
            end
            if length(wait_list)<vaccine_capacity && tt>=mutation_start+24 && vaccine_improved==1     %new vaccine available
                    wait_list2=find(vaccine_efficacy1_v<vaccine_efficacy1 & vaccinated==1 & dead==0 & detected==0);                    
                    vaccine_capacity_res=vaccine_capacity-length(wait_list);
                    %wait_list=[wait_list wait_list2];
                    if length(wait_list2)>=vaccine_capacity_res
                        new_vaccine2=vaccine_capacity_res;
                    else
                        new_vaccine2=length(wait_list2);
                    end
            end
        end
       
       if vaccine_strategy_1==1 %random vaccination
            vax=datasample(wait_list,new_vaccine,'Replace',false);
            if length(wait_list)<vaccine_capacity && tt>=mutation_start+24 && vaccine_improved==1 
            vax2=datasample(wait_list2,new_vaccine2,'Replace',false);
            end
        elseif vaccine_strategy_2==1 %priority to old
            priority=exp(age(wait_list));
            vax=datasample(wait_list,new_vaccine,'Replace',false,'Weights',priority);
            if length(wait_list)<vaccine_capacity && tt>=mutation_start+24 && vaccine_improved==1 
            priority2=exp(age(wait_list2));
            vax2=datasample(wait_list2,new_vaccine2,'Replace',false,'Weights',priority2);
            end
        elseif vaccine_strategy_3==1 %priority to workers --> equal priority to age 1&2!
            priority=age(wait_list);
            priority(priority==1)=2;
            priority(priority==3)=1;
            priority(priority==2)=3;
            priority=exp(priority);
            vax=datasample(wait_list,new_vaccine,'Replace',false,'Weights',priority);
            if length(wait_list)<vaccine_capacity && tt>=mutation_start+24 && vaccine_improved==1 
            priority2=exp(age(wait_list2));
            vax2=datasample(wait_list2,new_vaccine2,'Replace',false,'Weights',priority2);
            end
       end
       
        vaccinated(vax)=1;
        immunity_duration(vax)=1;
        immunity_weeks(vax)=round(normrnd(average_immunity_weeks,2,new_vaccine,1));%sd=2
        if length(wait_list)<vaccine_capacity && tt>=mutation_start+24 && vaccine_improved==1 
        vaccinated(vax2)=1;
        immunity_duration(vax2)=1;
        immunity_weeks(vax2)=round(normrnd(average_immunity_weeks,2,new_vaccine2,1));%sd=2
        end
        
        vacc_share(tt)=sum(vaccinated)/sum(dead==0);
        
        
       if mutation==1 && vaccine_improved==1 && tt>=mutation_start+24  %new vaccine available with higher efficacy against new variant
           vaccine_new=vax(vaccine_efficacy1_v(vax)<vaccine_efficacy1);
           vaccine_efficacy1_v(vaccine_new)=vaccine_efficacy1;
           vaccine_efficacy2_v(vaccine_new)=vaccine_efficacy2;
           if length(wait_list)<vaccine_capacity
           vaccine_new2=vax2(vaccine_efficacy1_v(vax2)<vaccine_efficacy1);
           vaccine_efficacy1_v(vaccine_new2)=vaccine_efficacy1;
           vaccine_efficacy2_v(vaccine_new2)=vaccine_efficacy2;
           end
       end
   end
    
    if epidemic==1 && t>=577 && sum(infected)>0
        duration(infected==1)=duration(infected==1)+1; %update duration of disease
        
        inf_cum(tt)=inf_cum(tt-1);
        inf_detected(tt)=inf_detected(tt-1);
        inf_cum_v(tt)=inf_cum_v(tt-1);
        inf_detected_v(tt)=inf_detected_v(tt-1);
        death_total(tt)=death_total(tt-1);
        
        shop_cons=[];
        for j=1:(F)
            if j<=F_b
            customers=find(pb_prev(1:W,1)==j);
            customers=[customers; find(pb_prev(1:W,2)==j)];
            if ~isempty(customers)
            shop_encounters=randi(length(customers),round(length(customers)/(n_shopcons/reduce_cons)),2);
            shop_encounters(:,1)=customers(shop_encounters(:,1));
            shop_encounters(:,2)=customers(shop_encounters(:,2));
            if length(shop_encounters)>1
            shop_cons=[shop_cons; shop_encounters];
            end
            end
            else
            customers=find(pl_prev(1:W,1)==j);
            customers=[customers; find(pl_prev(1:W,2)==j)];
            if ~isempty(customers)
            shop_encounters=randi(length(customers),round(length(customers)/(n_shopcons/reduce_cons)),2);
            shop_encounters(:,1)=customers(shop_encounters(:,1));
            shop_encounters(:,2)=customers(shop_encounters(:,2));
            if length(shop_encounters)>1
            shop_cons=[shop_cons; shop_encounters];
            end
            end
            end
        end


        shop_weights=ones(length(shop_cons),1);
        shop_cons=[shop_cons shop_weights];
        if reduce_cons==1
            fixed_cons=perm_cons;
        else
            fixed_cons=datasample(perm_cons,round(reduce_cons*size(perm_cons,1)),'Replace',false);
        end

        connections=[fixed_cons; work_cons; shop_cons];
        
        %% Mutation
        if mutation==1 && tt==mutation_start
        current_infected=find(infected==1 & detected==0 & duration<3);
        degrees=zeros(1,length(current_infected));
        for i=1:length(current_infected)
           patientcons=find(connections(:,1)==current_infected(i) | connections(:,2)==current_infected(i)); 
           patientcons=connections(patientcons,:);
           patientcons=find(susceptible(patientcons(:,1))==1 & infected(patientcons(:,2))==1 | susceptible(patientcons(:,2))==1 & infected(patientcons(:,1))==1);
           degrees(i)=length(patientcons);
        end
        [~, deg_order] = sort(degrees);
        newCurrent = current_infected(deg_order);
        newCurrent=flip(newCurrent);
        patient00=newCurrent(1:5);
        infected(patient00)=1;
        variant(patient00)=1;
        end
        
        con=find(susceptible(connections(:,1))==1 & infected(connections(:,2))==1 & variant(connections(:,2))==0 & detected(connections(:,2))==0 | susceptible(connections(:,2))==1 & infected(connections(:,1))==1 & variant(connections(:,1))==0  & detected(connections(:,1))==0);
        cons=connections(con,1:2);
        weights=connections(con,3);
        
        con2=find(susceptible(connections(:,1))==1 & infected(connections(:,2))==1 & variant(connections(:,2))==1 & detected(connections(:,2))==0 | susceptible(connections(:,2))==1 & infected(connections(:,1))==1 & variant(connections(:,1))==1  & detected(connections(:,1))==0);
        cons2=connections(con2,1:2);
        weights2=connections(con2,3);
        
        new_infections=max(1,round(transmission_rate(t)*length(con)));
        new_infections2=max(1,round(transmission_rate_v(t)*length(con2)));
        
       if length(con)<new_infections
          new_infections=length(con);
       end
       if length(con2)<new_infections2
          new_infections2=length(con2);
       end
       
       tosample=[1:size(cons,1)];
       con_sample=datasample(tosample,new_infections,'Replace',false,'Weights',weights);
       tosample2=[1:size(cons2,1)];
       con_sample2=datasample(tosample2,new_infections2,'Replace',false,'Weights',weights2);
       con_sample=con_sample+size(cons2,1);
       con_sample=[con_sample2 con_sample];
       new_infections=new_infections+new_infections2;
       cons=[cons2; cons];
       weights=[weights2; weights];
       inf_newp(tt)=new_infections;
       
       %% Contagion      
       for i = 1:new_infections
          con=cons(con_sample(i),:);
          if variant(con(1))+variant(con(2))==0
              %without variant with vaccine
              if infected(con(2))==1 && vaccinated(con(1))==1 || infected(con(1))==1 && vaccinated(con(2))==1
                  if rand<(1-vaccine_efficacy1)*(1-distancing_effect*distancing(con(1))-distancing_effect*distancing(con(2))) && susceptible(con(1))==1 && infected(con(2))==1 && detected(con(2))==0 && vaccinated(con(1))==1 || ...
                     rand<(1-vaccine_efficacy1)*(1-distancing_effect*distancing(con(1))-distancing_effect*distancing(con(2))) && susceptible(con(2))==1 && infected(con(1))==1 && detected(con(1))==0 && vaccinated(con(2))==1
                      newInf=con(susceptible(con)==1);
                      infected(con(1))=1;
                      infected(con(2))=1;
                      inf_cum(tt)=inf_cum(tt)+1;
                      inf_new(tt)=inf_new(tt)+1;
                      susceptible(con(1))=0;
                      susceptible(con(2))=0;
                       if serious(newInf)==1 && vaccinated(newInf)==1 && rand<vaccine_efficacy2
                           serious(newInf)=0;
                       end
                  end
              %without variant without vaccine
              else
                  if rand<(1-distancing_effect*distancing(con(1))-distancing_effect*distancing(con(2))) && susceptible(con(1))==1 && infected(con(2))==1 &&  detected(con(2))==0 && vaccinated(con(1))==0 || ...
                     rand<(1-distancing_effect*distancing(con(1))-distancing_effect*distancing(con(2))) && susceptible(con(2))==1 && infected(con(1))==1 &&  detected(con(1))==0 && vaccinated(con(2))==0  %|| ...
                      infected(con(1))=1;
                      infected(con(2))=1;
                      inf_cum(tt)=inf_cum(tt)+1;
                      inf_new(tt)=inf_new(tt)+1;
                      susceptible(con(1))=0;
                      susceptible(con(2))=0;
                  end
              end
          else
              %with variant with vaccine
             if infected(con(2))==1 && vaccinated(con(1))==1 || infected(con(1))==1 && vaccinated(con(2))==1
                  if rand<(1-vaccine_efficacy1_v(con(1)))*(1-distancing_effect_v*distancing(con(1))-distancing_effect_v*distancing(con(2))) && susceptible(con(1))==1 && infected(con(2))==1 && detected(con(2))==0 && vaccinated(con(1))==1 || ...
                     rand<(1-vaccine_efficacy1_v(con(2)))*(1-distancing_effect_v*distancing(con(1))-distancing_effect_v*distancing(con(2))) && susceptible(con(2))==1 && infected(con(1))==1 && detected(con(1))==0 && vaccinated(con(2))==1
                      newInf=con(susceptible(con)==1);
                      infected(con(1))=1;
                      infected(con(2))=1;
                      variant(con(1))=1;
                      variant(con(2))=1;
                      inf_cum(tt)=inf_cum(tt)+1;
                      inf_new(tt)=inf_new(tt)+1;
                      inf_cum_v(tt)=inf_cum_v(tt)+1;
                      inf_new_v(tt)=inf_new_v(tt)+1;
                      susceptible(con(1))=0;
                      susceptible(con(2))=0;
                      if serious(newInf)==1 && vaccinated(newInf)==1 && rand<vaccine_efficacy2_v(newInf)
                         serious(newInf)=0;
                      end
                  end
             else
                 %with variant without vaccine
                  if rand<(1-distancing_effect_v*distancing(con(1))-distancing_effect_v*distancing(con(2))) && susceptible(con(1))==1 && infected(con(2))==1 && detected(con(2))==0 && vaccinated(con(1))==0 || ...
                     rand<(1-distancing_effect_v*distancing(con(1))-distancing_effect_v*distancing(con(2))) && susceptible(con(2))==1 && infected(con(1))==1 && detected(con(1))==0  && vaccinated(con(2))==0
                      infected(con(1))=1;
                      infected(con(2))=1;
                      variant(con(1))=1;
                      variant(con(2))=1;
                      inf_cum(tt)=inf_cum(tt)+1;
                      inf_new(tt)=inf_new(tt)+1;
                      inf_cum_v(tt)=inf_cum_v(tt)+1;
                      inf_new_v(tt)=inf_new_v(tt)+1;
                      susceptible(con(1))=0;
                      susceptible(con(2))=0;
                  end  
             end
          end
       end
       
       %% Detection and quarantine
       for i=randperm(W)
           %if currently infected
           if infected(i)==1 && duration(i)<=max_dur(i)
               %random detection of disease
              if rand<detectprob && detected(i)==0 && serious(i)==0
                 detected(i)=1;
                 detect_new(tt)=detect_new(tt)+1;
                 inf_detected(tt)=inf_detected(tt)+1;
                 if variant(i)==1
                     detect_new_v(tt)=detect_new_v(tt)+1;
                     inf_detected_v(tt)=inf_detected_v(tt)+1;
                 end
                 %detected cases become inactive; firm loses worker
                 if active(i)==1
                 active(i)=0;
                 sickpay(i)=1;  %only previously active agents receive sick pay; others just continue to receive pension
                 end
                 Firm=Oc(i);
                   if Firm>0
                      if Firm<=F
                         Leff(Firm)=Leff(Firm)-1;
                      end
                      if Firm>F
                         Leff_k(Firm-F)=Leff_k(Firm-F)-1;
                      end
                   end
                 Oc(i)=0;
              end
              %all serious cases are detected & leave work
              if serious(i)==1 && detected(i)==0
                  detected(i)=1;
                  detect_new(tt)=detect_new(tt)+1;
                  inf_detected(tt)=inf_detected(tt)+1;
                  if variant(i)==1
                     detect_new_v(tt)=detect_new_v(tt)+1;
                     inf_detected_v(tt)=inf_detected_v(tt)+1;
                  end
                  if active(i)==1
                  active(i)=0;
                  sickpay(i)=1;  %only previously active agents receive sick pay; others just continue to receive pension
                  end
                  Firm=Oc(i);
                   if Firm>0
                      if Firm<=F
                         Leff(Firm)=Leff(Firm)-1;
                      end
                      if Firm>F
                         Leff_k(Firm-F)=Leff_k(Firm-F)-1;
                      end
                   end
                  Oc(i)=0;
              end
              
              %% Hospitalization
              if serious(i)==1 && dead(i)==0
                   if health_demanded(i)==0
                      health_demanded(i)=h_1*age(i)+shocks_healthcare(tt,i)*h_2;
                   end
                   if health_demanded(i)>0 && health_supplied(i)<health_demanded(i)
                      hdem=health_demanded(i)-health_supplied(i);
                      health_supplied(i)=health_supplied(i)+min(health_supply,hdem);
                      health_supply=max(0,health_supply-hdem);
                      health_rationing(tt)=health_rationing(tt)+max(0,hdem-health_supply);
                   end
                   if variant(i)==1
                      die=dieprob*dieprob_v;
                   else
                      die=dieprob;
                   end
                   if rand<die*age(i)^3+h_3*(health_demanded(i)-health_supplied(i))
                       infected(i)=0;
                       dead(i)=1;
                       variant(i)=0;
                       if age(i)==1
                          a_v1=a_v1+1;
                          dead_age1(tt)=dead_age1(tt)+1;
                       elseif age(i)==2
                          a_v2=a_v2+1;
                          dead_age2(tt)=dead_age2(tt)+1;
                       elseif age(i)==3
                          a_v3=a_v3+1;
                          dead_age3(tt)=dead_age3(tt)+1;
                       end
                       death_new(tt)=death_new(tt)+1;
                       death_total(tt)=death_total(tt)+1;
                       if detected(i)==0
                       detected(i)=1;
                       detect_new(tt)=detect_new(tt)+1;
                       inf_detected(tt)=inf_detected(tt)+1;
                           if variant(i)==1
                           detect_new_v(tt)=detect_new_v(tt)+1;
                           inf_detected_v(tt)=inf_detected_v(tt)+1;
                           end
                       end
                       Firm=Oc(i);
                       if Firm>0
                          if Firm<=F
                             Leff(Firm)=Leff(Firm)-1;
                          end
                          if Firm>F
                             Leff_k(Firm-F)=Leff_k(Firm-F)-1;
                          end
                       end
                       Oc(i)=-1;
                       health_demanded(i)=0;
                       health_supplied(i)=0;
                       health_supply=max(0,healthcare(t)-sum(health_supplied)-sum(health_supplied2));
                       dead_assets=dead_assets+PA(i);
                       PA(i)=0;
                       active(i)=0;
                       sickpay(i)=0;
                   end             
              end
           end
           
          if infected(i)==1 && duration(i)>max_dur(i) && dead(i)==0
          infected(i)=0;
          recovered(i)=1;
          variant(i)=0;
          health_demanded(i)=0;
          health_supplied(i)=0;
          detected(i)=0;
          health_supply=max(0,healthcare(t)-sum(health_supplied)-sum(health_supplied2));
          if age(i)<3
             active(i)=1;
          end
          sickpay(i)=0;
          end
       end
       else
        health_demanded(:)=0;
        health_supplied(:)=0;
    end
    
    %% Update immunity duration
    %update immunity duration
    if vaccine==1 && sum(vaccinated)>0
         vac=find(vaccinated==1); 
         for i=vac
             if immunity_duration(i)>0
                 immunity_duration(i)=max(0,immunity_duration(i)-1/immunity_weeks(i));
             else
                 vaccinated(i)=0;               %if immunity is over, the vaccinated status is lost
                 if rand<seriousprob(age(i))    %and individuals become exposed again to serious symptoms with a certain probability depending on age
                    serious(i)=1;
                 end
             end
        end
    end
    
    %calculate aggregate stats    
    sus_share(tt)=sum(susceptible)/sum(dead==0);
    inf_share(tt)=sum(infected)/sum(dead==0);
    rec_share(tt)=sum(recovered)/sum(dead==0);
    dead_share(tt)=sum(dead)/W;
    serious_share(tt)=sum(infected==1 & serious==1)/sum(dead==0);
    serious_share2(tt)=sum(dead==0 & serious==1)/sum(dead==0);

    rec=find(recovered==1 & dead==0);
    newrec=find(recovered==1 & dead==0 & recovered_duration==0);
    immu=normrnd(natural_immunity,2,[1,length(newrec)]);
    recovered_max(newrec)=immu;
    recovered_duration(rec)=recovered_duration(rec)+1;
    new_sus=find(recovered_duration>recovered_max);
    susceptible(new_sus)=1;
    detected(new_sus)=0;
    recovered(new_sus)=0;
    duration(new_sus)=0;
    recovered_duration(new_sus)=0;
    for i=1:length(new_sus)
        if rand<seriousprob(age(new_sus(i)))
        serious(new_sus(i))=1;
        else
        serious(new_sus(i))=0;
        end
    end
    
    
    %% "Normal" disease
    %increment duration for infected
    duration2(infected2==1)=duration2(infected2==1)+1;
    %agents who have exceeded the max. duration of the disease recover
    infected2(duration2>=max_dur2)=0;
    recovered2(duration2>=max_dur2)=1;
    duration2(recovered2==1)=0;
    health_demanded2(recovered2==1)=0;
    health_supplied2(recovered2==1)=0;
    health_supply=max(0,healthcare(t)-sum(health_supplied)-sum(health_supplied2));
    %recovered agents become active again
    active(recovered2==1 & age<3)=1;
    sickpay(recovered2==1)=0;
    %find susceptible agents
    sus=find(susceptible2==1 & dead==0);
    %randomly infect some susceptible agents
    newinf2=round(infprob2*sum(dead==0));
    infect=randsample(sus,newinf2,false);
    infected2(infect)=1;
    susceptible2(infect)=0;
    sick=find(infected2==1);
    svec=randperm(length(sick));
    for si=svec
           i=sick(si);
           if health_demanded2(i)==0
              health_demanded2(i)=h_1*age(i)+shocks_healthcare2(tt,i)*h_2;
           end
           if health_demanded2(i)>0 && health_supplied2(i)<health_demanded2(i)
              hdem=health_demanded2(i)-health_supplied2(i);
              health_supplied2(i)=health_supplied2(i)+min(health_supply,hdem);
              health_supply=max(0,health_supply-hdem);
              health_rationing2(tt)=health_rationing2(tt)+max(0,hdem-health_supply);
           end
           if active(i)==1
           active(i)=0;
           sickpay(i)=1;  %only previously active agents receive sick pay; others just continue to receive pension
           Firm=Oc(i);
            if Firm>0
              if Firm<=F
                 Leff(Firm)=Leff(Firm)-1;
              end
              if Firm>F
                 Leff_k(Firm-F)=Leff_k(Firm-F)-1;
              end
            end
           end
           Oc(i)=0;       
    end

    %a share of recovered become susceptible again
    rec=find(recovered2==1 & dead==0);
    new_sus=rand(1,length(rec));
    susceptible2(rec(new_sus<susprob2))=1;
    recovered2(susceptible2==1)=0;

    inf_share2(tt)=sum(infected2)/W;
    sus_share2(tt)=sum(susceptible2)/W;
    rec_share2(tt)=sum(recovered2)/W;
    distance(tt)=sum(distancing)/W;

%% Lockdown starts
if detect_new(tt)>=lockdown_threshold && t>=578 && lockdown_exp==1 && lockdown==0 && lock_ended==0
    lockdown=1;
    homeoffice(:)=1;
    homeoffice_k(:)=1;
    lock_start=t;
    lock_start2=tt;
    lock_ended=0;
    activel=active_lfirms;
    for i=1:length(deactivate_f)
            j=deactivate_f(i);
            activel(activel==j)=[];
            partners=find(partners_l(:,1)==j);
            if ~isempty(partners)
            partners_l(partners,1)=randi(length(activel),length(partners),1);
            partners_l(partners,1)=activel(partners_l(partners,1));
            end
    end
    work_cons=[];
    for j=1:(F)
        if active_f(j)==1
        colleagues=find(Oc==j);
        if homeoffice(j)==1
              colleagues=randsample(colleagues,round(reduce_workcons_lockdown*length(colleagues)));
        end
        if length(colleagues)>1
           work_cons=[work_cons; nchoosek(colleagues,2)];
        end
        end
    end
    for j=1:(N)
        if active_k(j)==1
        colleagues=find(Oc==F+j);
        if homeoffice_k(j)==1
              colleagues=randsample(colleagues,round(reduce_workcons_lockdown*length(colleagues)));
        end
        if length(colleagues)>1
            work_cons=[work_cons; nchoosek(colleagues,2)];
        end
        end
    end
    work_weights=2*ones(length(work_cons),1);
    work_cons=[work_cons work_weights];
end

% if death_new(tt)>0
% count=['new deaths ', num2str(death_new(tt))];
% disp(count);
% end


%% Epidmemic ends
if inf_detected(tt)>0 && sum(infected)==0 && epidemic_end==T+1
   epidemic_end=t+1;
   epidemic_end_tt=tt+1;
end

end
% 
% if t==590
%    infected(:)=0;
%    epidemic_end=t+1;
%    epidemic_end_tt=tt;
% end


labforce1(t)=sum(active==1)/W;
labforce2(t)=sum(active==1)/sum(age<3);

 end
 
%vaccination statistics
vax_start=(epidemic_start+11)*4-4;
%mutation_start_tt=mutation_start;
dead_after_vax=sum(death_new(vax_start:epidemic_end_tt));
total_deaths=sum(death_new);
total_det=sum(detect_new);
length_vax=epidemic_end_tt-vax_start;
vax_stats2=[total_deaths;total_det;epidemic_end_tt;length_vax];

a_v4=sum(dead_after_vax);
vax_stats=[a_v1; a_v2; a_v3; a_v4];


et=toc();

 