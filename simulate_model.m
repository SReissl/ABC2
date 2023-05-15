clear

%number of runs
S = 1;
%periods
T = 756;
%parfor s=1:S
%disp(s)
for s = 1:S
    feature('jit','off')
    feature('accel', 'off')
    %choose which series to output (also change in model script)
    [Y_real,lock_start2,lock_end,defGDP,debGDP,labforce1,labforce2,Lev_a,Lev_ab,Lev_al,Lev_af,Lev_ak,defaults,defaults_k,Hdemand,tot_healthcare,hospital_capacity,CPI_b,CPI_l,dprob,lock_duration,detect_new,inf_detected,distance,inf_share2,sus_share2,rec_share2,healthcare,serious_share,inf_cum,inf_new,inf_cum_v,inf_new_v,sus_share,inf_share,rec_share,death_total,gdp_deflator, Investment,EXP, consumption, Un,stock_bonds,wages_t, DC, totK,I,totE,bankruptcy_rate,et,redcons,dis_share,dis_prob_tt,vacc_share,epidemic_end,epidemic_end_tt,epidemic_start,death_new,vax_stats,dead_share,inf_current,detect_current,vax_stats2,dead_age1,dead_age2,dead_age3,detect_new_v] = model(s,T,0,0,52,52);
    %this is the data which will be saved for multiple runs
    Y(:,s)  = Y_real;
    P(:,s)  = gdp_deflator;
    Ig(:,s) = I;
    Inv(:,s)  = Investment;
    C(:,s)  = consumption;
    %DC(:,s) = DC;
    U(:,s) = 1-Un;
    W(:,s) = wages_t;
    elapsed(:,s) = et;
    G(:,s)=EXP;
    B(:,s)=stock_bonds;
    Infected(:,s)=inf_cum;
    Infected_new(:,s)=inf_new;
    Infected_c(:,s)=inf_current;
    Dead_s(:,s)=dead_share;
    Dead(:,s)=death_total;
    Dead_new(:,s)=death_new;
    Serious(:,s)=serious_share;
    Distance(:,s)=distance;
    Detected(:,s)=inf_detected;
    Detected_c(:,s)=detect_current;
    Detected_new(:,s)=detect_new;
    Lock_duration(s)=lock_duration;
    Detectprob(:,s)=dprob;
    bankrupt(:,s)=bankruptcy_rate;
    bankrupt2(:,s)=defaults+defaults_k;
    P_l(:,s)=CPI_l;
    P_b(:,s)=CPI_b;
    HealthD(:,s)=Hdemand;
    Lev_tot(:,s)=Lev_a;
    Lev_f(:,s)=Lev_af;
    Lev_k(:,s)=Lev_ak;
    Lev_b(:,s)=Lev_ab;
    Lev_l(:,s)=Lev_al;
    lab1(:,s)=labforce1;
    lab2(:,s)=labforce2;
    DefGDP(:,s)=defGDP;
    DebGDP(:,s)=debGDP;
    lstart(s)=lock_start2;
    lend(s)=lock_end;
    vax_deaths(:,s)=vax_stats;
    vax_data(:,s)=vax_stats2;
    dead_a1(:,s)=dead_age1;
    dead_a2(:,s)=dead_age2;
    dead_a3(:,s)=dead_age3;
    Infected_v(:,s)=inf_cum_v;
    Infected_v_new(:,s)=inf_new_v;
    Detected_v_new(:,s)=detect_new_v;
end

