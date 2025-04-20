function [x,y]=drawCO2(i,data,new)

TPumpIn = data.data.TPumpIn(end)+273.15;                   % K
TPumpOut = data.data.TPumpOut(end)+273.15;                 % K
TThrottleIn = data.data.TThrottleIn(end)+273.15;           % K
TThrottleOut = data.data.TThrottleOut(end)+273.15;         % K
pThrottleIn = (data.data.pThrottleIn(end)+1)/10;           % MPa
pThrottleOut = (data.data.pThrottleOut(end)+1)/10;         % MPa


n=200;
T=NaN(4,n);
s=NaN(4,n);

T(1,1:2)=[TPumpIn,TPumpOut];                                                     %compression
T(2,:)=linspace(TPumpOut,TThrottleIn,n);                                         %heating
T(3,1:2)=[TThrottleIn,TThrottleOut];                                             %expansion
T(4,:)=linspace(TPumpIn,TThrottleOut,n);                                         %isobar low

s(1,1:2)=CO2.s_pT([pThrottleOut,pThrottleIn],T(1,1:2));        %compression
s(2,:)=CO2.s_pT(pThrottleIn,T(2,:));                           %heating
s(3,1:2)=CO2.s_pT([pThrottleIn,pThrottleOut],T(3,1:2));        %expansion
s(4,:)=CO2.s_pT(pThrottleOut,T(4,:));                          %isobar low


T=(T-273.15)';
s=s';


x=reshape(s,numel(s),1);
[x,x_ind]=sort(x);


ind=repmat((1:4),n,1);
ind=ind(x_ind);
ind=sub2ind([numel(x),4],(1:numel(x))',ind);


T=reshape(T,numel(T),1);
T=T(x_ind);

y=NaN(numel(x),4);
y(ind)=T(:);
y=fillmissing(y,'linear',1,'EndValues','none');

% miss=isnan(x);
% y(miss,:)=[];
% x(miss)=[];



if new
    figure(i);
    
    % two phase region
    T_2p = linspace(CO2.Tt, CO2.Tc,100)';
    rho_l = CO2.SaturatedLiquidDensity(T_2p);
    rho_v = CO2.SaturatedVaporDensity(T_2p);
    s_l = CO2.s_rhoT(rho_l, T_2p);
    s_v = CO2.s_rhoT(rho_v, T_2p);
    plot([s_l,s_v], T_2p-273.15, 'k');
    
    % Critical isobar and isotherm
    T = linspace(-40, 150, 500) + 273.15;
    hold on
    plot(CO2.s_pT(CO2.pc, T), T-273.15);
    plot([min(s_l) max(s_v)], [CO2.Tc, CO2.Tc]-273.15);
    
    %remove visibility of children
    lines=get(gca,'Children');
    set(lines,'HandleVisibility','off');
    
    
    plot(x,y);
    hold off
    
    
    xlabel('Entropy in kJ/kgK');
    ylabel('Temperature in °C');
    
end






end
























% %% Base diagram
% figure(1); clf; hold on
% 
% % two phase region
% T_2p = linspace(CO2.Tt, CO2.Tc,100)';
% T_2p(1) = [];
% rho_l = CO2.SaturatedLiquidDensity(T_2p);
% rho_v = CO2.SaturatedVaporDensity(T_2p);
% s_l = CO2.s_rhoT(rho_l, T_2p);
% s_v = CO2.s_rhoT(rho_v, T_2p);
% plot(s_l, T_2p-273.15, 'k');
% plot(s_v, T_2p-273.15, 'k');
% 
% % Critical isobar and isotherm
% T = linspace(-40, 400, 500) + 273.15;
% plot(CO2.s_pT(CO2.pc, T), T-273.15);
% plot([min(s_l) max(s_v)], [CO2.Tc, CO2.Tc]-273.15);
% 
% Tc1 = 10+273.15; % K
% Tc2 = 15+273.15; % K
% Te1 = 200+273.15; % K
% Te2 = 170+273.15; % K
% ph = (80+1)/10; % MPa
% pl = (65+1)/10; % MPa
% h = gobjects(4,1);
% % compression
% h(1) = plot(CO2.s_pT([pl ph], [Tc1 Tc2]), [Tc1 Tc2]-273.15);
% % heating
% T = linspace(Tc2, Te1, 200);
% h(2) = plot(CO2.s_pT(ph, T), T-273.15);
% % expansion
% h(1) = plot(CO2.s_pT([ph pl], [Te1 Te2]), [Te1 Te2]-273.15);
% % isobar low
% T = linspace(Te2, Tc1, 200);
% h(2) = plot(CO2.s_pT(pl, T), T-273.15);
% 
% % xlim([-2 -0.4]);
% % ylim([-10 200]);
% xlabel('Entropy in kJ/kgK');
% ylabel('Temperature in °C');
