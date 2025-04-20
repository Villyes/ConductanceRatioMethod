%% Phase diagram

figure(2); clf; hold on
CO2.plotPhaseDiagram()
grid on

%% T s diagram

T_2p = linspace(CO2.Tt, CO2.Tc,100)';
T_2p(1) = [];
rho_l = CO2.rhoLiqSat(T_2p);
rho_v = CO2.rhoVapSat(T_2p);

p_l = CO2.p_rhoT(rho_l(:), T_2p(:));
p_v = CO2.p_rhoT(rho_v(:), T_2p(:));

s_l = CO2.s_rhoT(rho_l, T_2p);
s_v = CO2.s_rhoT(rho_v, T_2p);

T = linspace(-40, 400, 500) + 273.15;
figure(3); clf; hold on
plot(s_l, T_2p-273.15, 'k');
plot(s_v, T_2p-273.15, 'k');
h(1) = plot(CO2.s_pT(5, T), T-273.15, 'DisplayName', '50 bar');
h(2) = plot(CO2.s_pT(CO2.pc, T), T-273.15, 'DisplayName', 'critical isobar');
h(3) = plot(CO2.s_pT(12, T), T-273.15, 'DisplayName', '120 bar');
h(4) = plot(CO2.s_pT(22, T), T-273.15, 'DisplayName', '220 bar');

xlim([-2 0]);
ylim([-10 400]);
xlabel('entropy in kJ/kg K');
ylabel('temperature in °C');
grid on
legend(h(end:-1:1));



















