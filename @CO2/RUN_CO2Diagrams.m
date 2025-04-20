

Ttc = [21 380]';
ptc = [CO2.VaporPressure(Ttc(1)+273.15)*10+1e-2 220]';

Tsc = [35 380]';
psc = [80 220]';

Tlim = [-60 400];
plim = [5 1e4];
slim = [-2.5 0];

% rho and T along the cycle
[rhoCycleTc, TCycleTc] = cycle(ptc/10, Ttc+273.15, 1.1, 1.1);
pCycleTc = CO2.p_rhoT( rhoCycleTc, TCycleTc )*10;
% TCycleTc = TCycleTc - 273.15;
[rhoCycleSc, TCycleSc] = cycle(psc/10, Tsc+273.15, 1.2, 1.1);
pCycleSc = CO2.p_rhoT( rhoCycleSc, TCycleSc )*10;

% pT Diagram
f = figure(1); clf; hold on
ax = gca;

T = linspace(CO2.Tt-273.15, CO2.Tc-273.15, 100)';
h(1) = plot(T, CO2.VaporPressure(T+273.15)*10, 'k', 'DisplayName', 'phase boundaries');
T = linspace(CO2.Tt-273.15, Tsc(2), 100)';
plot(T, CO2.MeltingPressure(T+273.15)*10, 'k', 'DisplayName', 'melting pressure');
% plot(CO2.Tc-273.15, CO2.pc*10, '*', 'DisplayName', 'Critical point')

MT = [1 1 0 0 1; 0 0 1 1 0]';
Mp = [0 1 1 0 0; 1 0 0 1 1]';

% h(2) = plot(MT*Ttc, Mp*ptc, 'DisplayName', 'Trans critical cycle', 'LineWidth', 2);
% h(3) = plot(MT*Tsc, Mp*psc, 'DisplayName', 'Super critical cycle', 'LineWidth', 2);
h(2) = plot(TCycleTc-273.15, pCycleTc, 'DisplayName', 'trans critical cycle', 'LineWidth', 1.5);
h(3) = plot(TCycleSc-273.15, pCycleSc, 'DisplayName', 'super critical cycle', 'LineWidth', 1.5);

plot([CO2.Tc CO2.Tc]-273.15, ([CO2.pc CO2.MeltingPressure(CO2.Tc)])*10, 'k--',  'DisplayName', 'fluid-Super critical');
plot([CO2.Tc-273.15 Tsc(2)], ([CO2.pc CO2.pc])*10, 'k--', 'DisplayName', 'gas-Super critical');

ylabel('pressure in bar');
xlabel('temperature °C');
ax.YScale = 'log';
xlim(Tlim);
ylim(plim);
legend(h);
f.Units = 'centimeters';
f.Position(3:4) = [12 6];

% Ts diagram
stc = CO2.s_pT(ptc/10,Ttc+273.15);
ssc = CO2.s_pT(psc/10,Tsc+273.15);

f = figure(2); clf; hold on
ax = gca;

% two phase region
T_2p = linspace(CO2.Tt, CO2.Tc, 100)';
rho_l = CO2.SaturatedLiquidDensity(T_2p);
rho_v = CO2.SaturatedVaporDensity(T_2p);
s_l = CO2.s_rhoT(rho_l, T_2p);
s_v = CO2.s_rhoT(rho_v, T_2p);
h = plot([s_l s_v], T_2p-273.15, 'k', 'DisplayName', 'phase boundaries');

% Critical isobar and isotherm
T = linspace(Tlim(1), Tlim(2), 500) + 273.15;
hold on
plot(CO2.s_pT(CO2.pc, T), T-273.15, 'k--', 'DisplayName', 'Critical isobar');
plot([min(s_l) max(s_v)], [CO2.Tc, CO2.Tc]-273.15, 'k--', 'DisplayName', 'critical isotherm');

% T = linspace(Ttc(1), Ttc(2), 100);
% s = [CO2.s_pT(ptc(1)/10, T+273.15); CO2.s_pT(ptc(2)/10, fliplr(T)+273.15)];
% s = [s; s(1)];
% h(2) = plot(s, [T fliplr(T) T(1)], 'DisplayName', 'Trans critical cycle', 'LineWidth', 2);
h(2) = plot(CO2.s_rhoT(rhoCycleTc, TCycleTc), TCycleTc-273.15, 'DisplayName', 'trans critical cycle', 'LineWidth', 1.5);


% T = linspace(Tsc(1), Tsc(2), 100);
% s = [CO2.s_pT(psc(1)/10, T+273.15); CO2.s_pT(psc(2)/10, fliplr(T)+273.15)];
% s = [s; s(1)];
% h(3) = plot(s, [T fliplr(T) T(1)], 'DisplayName', 'Super critical cycle', 'LineWidth', 2);
h(3) = plot(CO2.s_rhoT(rhoCycleSc, TCycleSc), TCycleSc-273.15, 'DisplayName', 'super critical cycle', 'LineWidth', 1.5);

xlabel('specific entropy in kJ/kgK');
ylabel('temperature in °C');
xlim(slim);
ylim(Tlim);
l = legend(h);
l.Location = 'northwest';
f.Units = 'centimeters';
f.Position(3:4) = [12 6];

%%
% T = linspace(Tlim(1), Tlim(2), 200);
% p = logspace(log10(plim(1)), log10(plim(2)), 200);
% [X, Y] = meshgrid(T+273.15, p/10);
% Z = CO2.rho_pT(Y,X); Z = reshape(Z, size(X));
% 
% figure(1)
% 
% contour(X+273.15, Y*10, Z)

%%

function [rho, T] = cycle(pIn, TIn, n1, n2)
    T = [];
    rho = [];

    % polytropic compression from T1 p1 to p2
    n = n1;
    p_ = logspace(log10(pIn(1)), log10(pIn(2)), 500)';
    T_ = TIn(1)*(p_/pIn(1)).^((n-1)/n);
%     for i = 2:numel(p_)
%         n = CO2.cp_pT(p_(i-1),T_(i-1))/CO2.cv_pT(p_(i-1),T_(i-1));
%         T_(i) = T_(i-1)*(p_(i)/p_(i-1)).^((n-1)/n);
%     end
    rho_ = CO2.rho_pT(p_, T_);
    
    T = [T; T_];
    rho = [rho; rho_];
    
    % isobaric heating from T(end) p2 to T2
    T_ = linspace(T_(end), TIn(2), 500)';
    rho_ = CO2.rho_pT(pIn(2), T_);

    T = [T; T_];
    rho = [rho; rho_];

    % adiabatic expansion from T2 p2 to p1
    n = n2;
    p_ = logspace(log10(pIn(2)), log10(pIn(1)), 500)';
    T_ = TIn(2)*(p_/pIn(2)).^((n-1)/n);
    rho_ = CO2.rho_pT(p_, T_);
    
    T = [T; T_];
    rho = [rho; rho_];

    % isobaric cooling from T(end) p1 to T1
    T_ = linspace(T_(end), TIn(1), 500)';
    rho_ = CO2.rho_pT(pIn(1), T_);

    T = [T; T_];
    rho = [rho; rho_];

%     T = [T; T(1)];
%     rho = [rho; rho(1)];
end