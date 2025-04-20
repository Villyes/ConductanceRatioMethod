classdef CO2
    %CO2 Calculate thermodynamic properties of CO2 in the fluid region
    %   Based on R. Span, W. Wagner 1994: A New Equation of State for
    %   Carbon Dioxide Covering the Fluid Region from the Tripple-Point
    %   Temperature to 1100 K at Pressures up to 800 MPa
    
    properties(Constant)
        M = 44.0098; % g/mol Molar mass
        RM = 8.314510; % J/(mol K) Molar gas constant
        R = 0.1889241; % kJ/(kg K) Specific gas constant
        
        % Triple point
        Tt = 216.592; % K
        pt = 0.51795; % MPa
        rhotLiquid = 1178.53; % kg/m3
        rhotVapor = 13.7614; % kg/m3
        
        % Critical point
        Tc = 304.1282; % K
        pc = 7.3773; % MPa
        rhoc = 467.6; % kg/m3
        
        T0 = 298.15; % K Reference temperature
        p0 = 0.1013258 % MPa Reference pressure
        h0 = 0; % kJ/kg Reference enthalpy in the ideal gas state at T0
        s0 = 0; % kJ/kg K Reference entropy in the ideal gas state at T0, p0
    end
    
    methods(Static) % Phase boundaries
        function pMelting = MeltingPressure(T)
            % According to Eq. (3.10)
            a1 = 1955.5390;
            a2 = 2055.4593;
            
            pMelting = nan(size(T));
            idx = T >= CO2.Tt;
            
            T_ = T(idx)/CO2.Tt - 1;
            pMelting(idx) = CO2.pt * ( 1 + a1*T_ + a2*T_.^2 );
        end
        
        function pSublimation = SublimationPressure(T)
            % According to Eq. (3.12)
            a1 = -14.740846;
            a2 = 2.4327015;
            a3 = -5.3061778;
            
            pSublimation = nan(size(T));
            idx = T <= CO2.Tt;
            
            T_ = 1 - T(idx)/CO2.Tt;
            pSublimation(idx) = CO2.pt * exp( CO2.Tt ./ T(idx) .*  ...
                ( a1*T_ + a2*T_.^1.9 + a3*T_.^2.9 ) );
        end
        
        function pSaturation = VaporPressure(T)
            % According to Eq. (3.13)
            a1 = -7.0602087;
            a2 = 1.9391218;
            a3 = -1.6463597;
            a4 = -3.2995634;
            t1 = 1;
            t2 = 1.5;
            t3 = 2;
            t4 = 4;
            
            pSaturation = nan(size(T));
            idx = CO2.Tt <= T & T <= CO2.Tc;
            
            T_ = 1 - T(idx)/CO2.Tc;
            pSaturation(idx) = CO2.pc * exp( CO2.Tc ./ T(idx) .* ...
                (a1*T_.^t1 + a2*T_.^t2 + a3*T_.^t3 + a4*T_.^t4 ) ); 
        end
        
        function rhoLiquid = SaturatedLiquidDensity(T)
            % According to Eq. (3.14)
            a1 = 1.9245108;
            a2 = -0.62385555;
            a3 = -0.32731127;
            a4 = 0.39245142;
            t1 = 0.34;
            t2 = 0.5;
            t3 = 10/6;
            t4 = 11/6;

            rhoLiquid = nan(size(T));
            idx = CO2.Tt <= T & T <= CO2.Tc;
            
            T_ = 1 - T(idx)/CO2.Tc;
            rhoLiquid(idx) = CO2.rhoc * exp( ...
                (a1*T_.^t1 + a2*T_.^t2 + a3*T_.^t3 + a4*T_.^t4 ) ); 
        end
        
        function rhoVapor = SaturatedVaporDensity(T)
            % According to Eq. (3.15)
            a1 = -1.7074879;
            a2 = -0.82274670;
            a3 = -4.6008549;
            a4 = -10.111178;
            a5 = -29.742252;
            t1 = 0.34;
            t2 = 0.5;
            t3 = 1;
            t4 = 7/3;
            t5 = 14/3;

            rhoVapor = nan(size(T));
            idx = CO2.Tt <= T & T <= CO2.Tc;
            
            T_ = 1 - T(idx)/CO2.Tc;
            rhoVapor(idx) = CO2.rhoc * exp( ...
                (a1*T_.^t1 + a2*T_.^t2 + a3*T_.^t3 + a4*T_.^t4 + a5*T_.^t5 ) ); 
        end
        
        function state = state_rhoT(rho, T)
            % States
            %   solid          -1
            %   liquid          0
            %   two-phase       1
            %   gas             2
            %   super critical  3 
            state = -ones(max(numel(rho), numel(T)), 1);
            rhoLiquid = SaturatedLiquidDensity(T);
            rhoVapor = SaturatedVaporDensity(T);
            
            
            
        end
        
        function plotPhaseDiagram()
            newplot(); hold on
            
            T = [linspace(CO2.Tt-50, CO2.Tt, 50) linspace(CO2.Tt, 1100, 500)];
            TC = T - 273.15;
            plot(TC, CO2.SublimationPressure(T)*10, 'DisplayName', 'Sublimation');
            plot(TC, CO2.MeltingPressure(T)*10, 'DisplayName', 'Melting');
            plot(TC, CO2.VaporPressure(T)*10, 'DisplayName', 'Vapor');
            plot(CO2.Tt-273.15, CO2.pt*10, '*', 'DisplayName', 'Triple point')
            plot(CO2.Tc-273.15, CO2.pc*10, '*', 'DisplayName', 'Critical point')
            legend('show')
            plot([CO2.Tc CO2.Tc]-273.15, [CO2.pc CO2.MeltingPressure(CO2.Tc)]*10, 'k--',  'DisplayName', 'Fluid-Super critical');
            plot([CO2.Tc T(end)]-273.15, [CO2.pc CO2.pc]*10, 'k--', 'DisplayName', 'Gas-Super critical');
            ylabel('p in bar');
            xlabel('T in °C');
            set(gca, 'yscale', 'log');
        end
        
        function plotPhaseDiagramLin()
            newplot(); hold on
            
            T = [linspace(CO2.Tt-50, CO2.Tt, 50) linspace(CO2.Tt, 1100, 500)];
            plot(T-273.15, (CO2.SublimationPressure(T)), 'DisplayName', 'Sublimation');
            plot(T-273.15, (CO2.MeltingPressure(T)), 'DisplayName', 'Melting');
            plot(T-273.15, (CO2.VaporPressure(T)), 'DisplayName', 'Vapor');
            plot(CO2.Tt-273.15, (CO2.pt), '*', 'DisplayName', 'Triple point')
            plot(CO2.Tc-273.15, (CO2.pc), '*', 'DisplayName', 'Critical point')
            legend('show')
            plot([CO2.Tc CO2.Tc]-273.15, ([CO2.pc CO2.MeltingPressure(CO2.Tc)]), 'k--',  'DisplayName', 'Fluid-Super critical');
            plot([CO2.Tc T(end)]-273.15, ([CO2.pc CO2.pc]), 'k--', 'DisplayName', 'Gas-Super critical');
            ylabel('p / MPa');
            xlabel('T / °C');
        end
    end
    
    methods(Static) % Properties in (rho, T) 
        % Thermodynamic properties as functions of density in kg/m³ and
        % temperature K according to Table 3.
        % A validity check is performed. For points in the two-phase region
        % nan is returned. If the validity check is not required the
        % internal functions with the suffix rhoT_i can be used instead.
        
        function p = p_rhoT(rho, T)
            % Pressure in MPa as function of density in kg/m3 and temperature in K
            [v, rho, T] = CO2.isValid_rhoT(rho, T);
            p = nan(size(v));
            p(v) = CO2.p_rhoT_i(rho, T);
        end
        
        function s = s_rhoT(rho, T)
            % Entropy in kJ/(kg K) as function of density in kg/m3 and temperature in K
            [v, rho, T] = CO2.isValid_rhoT(rho, T);
            s = nan(size(v));
            s(v) = CO2.s_rhoT_i(rho, T);
        end
        
        function u = u_rhoT(rho, T)
            % Internal energy in kJ/kg as function of density in kg/m3 and temperature in K
            [v, rho, T] = CO2.isValid_rhoT(rho, T);
            u = nan(size(v));
            u(v) = CO2.u_rhoT_i(rho, T);
        end
        
        function cv = cv_rhoT(rho, T)
            % Isochoric heat capacity in kJ/(kg K) as function of density in kg/m3 and temperature in K
            [v, rho, T] = CO2.isValid_rhoT(rho, T);
            cv = nan(size(v));
            cv(v) = CO2.cv_rhoT_i(rho, T);
        end        
 
        function h = h_rhoT(rho, T)
            % Enthalpy in kJ/kg as function of density in kg/m3 and temperature in K
            [v, rho, T] = CO2.isValid_rhoT(rho, T);
            h = nan(size(v));
            h(v) = CO2.h_rhoT_i(rho, T);
        end

        function cp = cp_rhoT(rho, T)
            % Isobaric heat capacity in kJ/(kg K) as function of density in kg/m3 and temperature in K
            [v, rho, T] = CO2.isValid_rhoT(rho, T);
            cp = nan(size(v));
            cp(v) = CO2.cp_rhoT_i(rho, T);
        end        

        function cs = cs_rhoT(rho, T)
            % Saturated liquid heat capacity in kJ/(kg K) as function of density in kg/m3 and temperature in K
            [v, rho, T] = CO2.isValid_rhoT(rho, T);
            cs = nan(size(v));
            cs(v) = CO2.cs_rhoT_i(rho, T);
        end
        
        function w = w_rhoT(rho, T)
            % Speed of sound in m/s as function of density in kg/m3 and temperature in K
            [v, rho, T] = CO2.isValid_rhoT(rho, T);
            w = nan(size(v));
            w(v) = CO2.w_rhoT_i(rho, T);
        end
        
        function mu = mu_rhoT(rho, T)
            % Joule-Thomson coefficient in K/MPa as function of density in kg/m3 and temperature in K
            [v, rho, T] = CO2.isValid_rhoT(rho, T);
            mu = nan(size(v));
            mu(v) = CO2.mu_rhoT_i(rho, T);
        end
        
        function f = f_rhoT(rho, T)
            % Fugacity as function of density in kg/m3 and temperature in K
            [v, rho, T] = CO2.isValid_rhoT(rho, T);
            f = nan(size(v));
            f(v) = CO2.f_rhoT_i(rho, T);
        end
    end
    
    methods(Static) % Virial coefficients in T 
        % Virial coefficients as functions of temperature according to Table 3.

        function B = B_T(T)
            % Second virial coefficient in m^3/kg as function of temperature in K
            t = CO2.Tc ./ T(:);
            B = CO2.dphir_dd(1e-18,t) / CO2.rhoc;
        end

        function C = C_T(T)
            % Third virial coefficient in (m^3/kg)^2 as function of temperature in K
            t = CO2.Tc ./ T(:);
            C = CO2.ddphir_dddd(1e-18,t) / CO2.rhoc.^2;
        end
    end
    
    methods(Static) % Properties in (p, T)
        % Thermodynamic properties as functions of pressure in MPa and temperature in K
        
        function rho = rho_pT(p, T)
            % Density in kg/m3 as function of pressure in MPa and temperature in K
            p = p(:); T = T(:);
            % This function is used in all p,T functions to compute the
            % density. Then the requested property is computed with the
            % rho,T functions.
            if isscalar(T)
                T = repmat(T, size(p));
            end
            % Reduced density of the ideal gas
            d_ig = p(:).*1e3/CO2.R./T(:) / CO2.rhoc;
            % Reduced density of liquid by linear interpolation between p_v and pt
            p_v = CO2.VaporPressure(T(:));
%             p_m = CO2.MeltingPressure(T(:));
%             d_l = ( CO2.SaturatedLiquidDensity(T(:)).*( p_m-p(:) ) + CO2.rhotLiquid*( p(:)-p_v ) )./( p_m-p_v )/CO2.rhoc;
            d_l = CO2.SaturatedLiquidDensity(T(:))/CO2.rhoc;
            d_v = CO2.SaturatedVaporDensity(T(:))/CO2.rhoc;
            isLiquid = T < CO2.Tc & p > p_v;

            % Compute reduced density with inverse function
            t = CO2.Tc ./ T(:);
            d = zeros(size(d_ig));
            for i = 1:numel(d)
                % Reset starting value if there was a gas/liquid change or
                % if the last result was invalid
                if i == 1 || (T(i-1) < CO2.Tc && T(i) < CO2.Tc && isLiquid(i)~=isLiquid(i-1)) || ~isfinite(d(i-1))
                    % Use new starting value
                    if T(i)>CO2.Tc
                        d0 = d_ig(i);
                    elseif T(i) <= CO2.Tt
                        d0 = CO2.rhotLiquid;
                    elseif isLiquid(i)
                        d0 = d_l(i);
                    else
                        d0 = d_v(i);
                    end
                else
                    % Use result of last step as starting value
                    d0 = d(i-1);
                end
                
%                 d(i) = CO2.d_fzero(t(i), d_ig(i), d0);
                d(i) = CO2.d_newton(t(i), d_ig(i), d0);
            end
            
            rho = d.*CO2.rhoc;
        end
        
        function s = s_pT(p, T)
            % Entropy in kJ/(kg K) as function of pressure in MPa and temperature in K
            rho = CO2.rho_pT(p, T);
            s = CO2.s_rhoT_i(rho, T);
        end
        
        function u = u_pT(p, T)
            % Internal energy in kJ/kg as function of pressure in MPa and temperature in K
            rho = CO2.rho_pT(p, T);
            u = CO2.u_rhoT_i(rho, T);
        end
        
        function cv = cv_pT(p, T)
            % Isochoric heat capacity in kJ/(kg K) as function of pressure in MPa and temperature in K
            rho = CO2.rho_pT(p, T);
            cv = CO2.cv_rhoT_i(rho, T);
        end        
 
        function h = h_pT(p, T)
            % Enthalpy in kJ/kg as function of pressure in MPa and temperature in K
            rho = CO2.rho_pT(p, T);
            h = CO2.h_rhoT_i(rho, T);
        end

        function cp = cp_pT(p, T)
            % Isobaric heat capacity in kJ/(kg K) as function of pressure in MPa and temperature in K
            rho = CO2.rho_pT(p, T);
            cp = CO2.cp_rhoT_i(rho, T);
        end        

        function cs = cs_pT(p, T)
            % Saturated liquid heat capacity in kJ/(kg K) as function of pressure in MPa and temperature in K
            rho = CO2.rho_pT(p, T);
            cs = CO2.cs_rhoT_i(rho, T);
        end
        
        function w = w_pT(p, T)
            % Speed of sound in m/s as function of pressure in MPa and temperature in K
            rho = CO2.rho_pT(p, T);
            w = CO2.w_rhoT_i(rho, T);
        end
        
        function mu = mu_pT(p, T)
            % Joule-Thomson coefficient in K/MPa as function of pressure in MPa and temperature in K
            rho = CO2.rho_pT(p, T);
            mu = CO2.mu_rhoT_i(rho, T);
        end
        
        function f = f_pT(p, T)
            % Fugacity as function of pressure in MPa and temperature in K
            rho = CO2.rho_pT(p, T);
            f = CO2.f_rhoT_i(rho, T);
        end
    end
    
    methods(Static, Hidden) % Properties in (rho, T) internal 
        % Thermodynamic properties as functions of density and temperature
        % according to Table 3.
        % No validity check is performed for rho and T.

        function p = p_rhoT_i(rho, T)
            % Pressure in MPa as function of density in kg/m3 and temperature in K
            d = rho(:) ./ CO2.rhoc;
            t = CO2.Tc ./ T(:);
            p = rho.*CO2.R.*T(:).*( 1 + d.*CO2.dphir_dd(d, t) ) * 1e-3;
        end

        function s = s_rhoT_i(rho, T)
            % Entropy in kJ/(kg K) as function of density in kg/m3 and temperature in K
            d = rho(:) ./ CO2.rhoc;
            t = CO2.Tc ./ T(:);
            s = CO2.R.*( t.*( CO2.dphi0_dt(d,t)+CO2.dphir_dt(d,t) ) - CO2.phi0(d,t) - CO2.phir(d,t) );
        end

        function u = u_rhoT_i(rho, T)
            % Internal energy in kJ/kg as function of density in kg/m3 and temperature in K
            d = rho(:) ./ CO2.rhoc;
            t = CO2.Tc ./ T(:);
            u = CO2.R.*T(:).*( t.*( CO2.dphi0_dt(d,t)+CO2.dphir_dt(d,t) ) );
        end

        function cv = cv_rhoT_i(rho, T)
            % Isochoric heat capacity in kJ/(kg K) as function of density in kg/m3 and temperature in K
            d = rho(:) ./ CO2.rhoc;
            t = CO2.Tc ./ T(:);
            cv = CO2.R.*( -t.^2.*( CO2.ddphi0_dtdt(d,t)+CO2.ddphir_dtdt(d,t) ) );
        end        

        function h = h_rhoT_i(rho, T)
            % Enthalpy in kJ/kg as function of density in kg/m3 and temperature in K
            d = rho(:) ./ CO2.rhoc;
            t = CO2.Tc ./ T(:);
            h = CO2.R.*T(:).*( 1 + t.*( CO2.dphi0_dt(d,t)+CO2.dphir_dt(d,t) ) + d.*CO2.dphir_dd(d,t) );
        end

        function cp = cp_rhoT_i(rho, T)
            % Isobaric heat capacity in kJ/(kg K) as function of density in kg/m3 and temperature in K
            d = rho(:) ./ CO2.rhoc;
            t = CO2.Tc ./ T(:);
            cp = CO2.R.*( -t.^2.*( CO2.ddphi0_dtdt(d,t)+CO2.ddphir_dtdt(d,t) ) + ...
                ( 1 + d.*CO2.dphir_dd(d,t) - d.*t.*CO2.ddphir_dddt(d,t) ).^2 ./ ...
                ( 1 + 2*d.*CO2.dphir_dd(d,t) + d.^2.*CO2.ddphir_dddd(d,t) ) );
        end        

        function cs = cs_rhoT_i(rho, T)
            % Saturated liquid heat capacity in kJ/(kg K) as function of density in kg/m3 and temperature in K
            d = rho(:) ./ CO2.rhoc;
            t = CO2.Tc ./ T(:);
            cs = CO2.R.*( -t.^2.*( CO2.ddphi0_dtdt(d,t)+CO2.ddphir_dtdt(d,t) ) + ...
                ( 1 + d.*CO2.dphir_dd(d,t) - d.*t.*CO2.ddphir_dddt(d,t) ) ./ ...
                ( 1 + 2*d.*CO2.dphir_dd(d,t) + d.^2.*CO2.ddphir_dddd(d,t) ) .* ...
                ( ( 1 + d.*CO2.dphir_dd(d,t) - d.*t.*CO2.ddphir_dddt(d,t) ) - ...
                CO2.rhoc/CO2.R./d .* CO2.dps_dT(T) ));
        end

        function w = w_rhoT_i(rho, T)
            % Speed of sound in m/s as function of density in kg/m3 and temperature in K
            d = rho(:) ./ CO2.rhoc;
            t = CO2.Tc ./ T(:);
            w = sqrt( CO2.R.*1e3.*T(:).*( ...
                1 + 2*d.*CO2.dphir_dd(d,t) + d.^2.*CO2.ddphir_dddd(d,t) - ...
                ( 1 + d.*CO2.dphir_dd(d,t) - d.*t.*CO2.ddphir_dddt(d,t) ).^2 ./ ...
                ( t.^2.*( CO2.ddphi0_dtdt(d,t) + CO2.ddphir_dtdt(d,t) ) ) ...
                ));
        end

        function mu = mu_rhoT_i(rho, T)
            % Joule-Thomson coefficient in K/MPa as function of density in kg/m3 and temperature in K
            d = rho(:) ./ CO2.rhoc;
            t = CO2.Tc ./ T(:);
            mu = 1e3./( CO2.R*rho ) .* ... 
                ( -d.*CO2.dphir_dd(d,t) - d.^2.*CO2.ddphir_dddd(d,t) - d.*t.*CO2.ddphir_dddt(d,t) ) ./ ...
                ( ( 1 + d.*CO2.dphir_dd(d,t) - d.*t.*CO2.ddphir_dddt(d,t) ).^2 - ...
                t.^2.*( CO2.ddphi0_dtdt(d,t) + CO2.ddphir_dtdt(d,t) ).*( 1 + 2.*d.*CO2.dphir_dd(d,t) + d.^2.*CO2.ddphir_dddd(d,t) ) );
        end

        function f = f_rhoT_i(rho, T)
            % Fugacity as function of density in kg/m3 and temperature in K
            d = rho(:) ./ CO2.rhoc;
            t = CO2.Tc ./ T(:);
            f = exp( CO2.phir(d,t) + d.*CO2.dphir_dd(d,t) - log( 1 + d.*CO2.dphir_dd(d,t) ) );
        end
    end

    methods(Static, Hidden) % Auxilary functions
        function d = d_fzero(t, d_ig, d0)
            % Find the reduced density that satisfies p(d,t) = p
            
            % The corresponding problem is finding the zero point of
            % f(d) = d + d^2*dphir_dd - d_ig = 0
            % t ... reduced temperature
            % d_ig ... reduced density of the ideal gas
            % d0 ... starting point
            
            % Two step algorithm:
            % 1) Use newtons' algorithm to find an interval that contains
            %   the zero point. i.e. sign(f(a)) ~= sign(f(b))
            % 2) Use fzero type algorithm to find zero point
            %   Based on Brent 1973
            %   Algorithms for minimization without derivatives

            reltol = 1e-6;
            abstol = 1e-4;
            maxiter = 1e3;
            k = 0;
            
            % Use newtons' algorithm to find an interval that contains
            % the zero point. i.e. sign(f(a)) ~= sign(f(b))            
            
            % b is the best approximation to the zero point so far
            % a is the value of the last step
            a = d0;
            b = d0;
            dphir = CO2.dphir_dd(d0,t);
            fb = d0 + d0^2*dphir - d_ig;
            fa = fb;
            while sign(fb) == sign(fa) && k <= maxiter
                % Move b to a
                a = b;
                fa = fb;
                
                % Newton step f(x) / f'(x)
                delta = fa / ( 1+2*a*dphir + a.^2*CO2.ddphir_dddd(a,t) );
                % Limit step
                if abs(delta)/a > 0.5
                    delta = 0.5*a*sign(delta);
                end
                b = a - delta;
                
                % Compute new function value
                dphir = CO2.dphir_dd(b,t);
                fb = b + b^2*dphir - d_ig;
                
                k = k+1;
            end

            % Use fzero type algorithm to find zero point
            %   Based on Brent 1973
            %   Algorithms for minimization without derivatives
            
            % Target function
            f = @(d) d+d.^2.*CO2.dphir_dd(d,t) - d_ig;
            
            % b is the best approximation to the zero point so far
            % a is the value of the last step
            % c is the best approximation where f(c) opposite sign of f(b)
                        
            % Auxilary variables
            c = a;
            fc = fa;
            d = b - a;
            e = d;
            
            tol = ( 2*eps + reltol )*abs(b);
            m = 0.5 * ( c-b );
            
            while ( abs(m) > tol || abs(fb) > abstol ) && k < maxiter
                % Make sure that f(b) is the value closest to zero
                if abs(fc) < abs(fb)
                    a = b; b = c; c = a;
                    fa = fb; fb = fc; fc = fa;
                end
                
                tol = ( 2*eps + reltol )*abs(b);
                m = 0.5 * ( c-b );
                
                if abs(e) < tol && abs(fa) <= abs(fb)
                    % Bisection
                    d = m;
                    e = m;
                else
                    % Interpolation
                    s = fb/fa;
                    if a == c
                        % Linear interpolation
                        p = 2*m*s;
                        q = 1-s;
                    else
                        % Inverse quadratic interpolation
                        q = fa/fc;
                        r = fb/fc;
                        p = s*( 2*m*q*( q-r ) - ( b-a )*( r-1) );
                        q = ( q-1 )*( r-1 )*( s-1 );
                    end
                    
                    if p > 0
                        q = -q;
                    else
                        p = -p;
                    end
                    
                    s = e;
                    e = d;
                    
                    if 2*p <  3*m*q-abs(tol*q) || p < abs(0.5*s*q)
                        d = p/q;
                    else
                        d = m;
                        e = m;
                    end
                end
                
                % Move b to a
                a = b;
                fa = fb;
                
                % Calculate new b
                if abs(d) > tol
                    b = b + d;
                elseif m > 0
                    b = b + tol;
                else
                    b = b - tol;
                end
                fb = f(b);
                
                if ( fb>0 ) == ( fc>0 )
                    c = a;
                    fc = fa;
                    d = b - a;
                    e = d;
                end
                
                k = k + 1;
            end % while
            
            % Return reduced density value
            d = b;
        end
        
        function d = d_newton(t, d_ig, d0)
            % Find the reduced density that satisfies p(d,t) = p
            
            % The corresponding problem is finding the zero point of
            % f(d) = d + d^2*dphir_dd - d_ig = 0
            % t ... reduced temperature
            % d_ig ... reduced density of the ideal gas
            % d0 ... starting point
            
            % Uses newtons' algorithm
            
            reltol = 1e-6;
            abstol = 1e-4;
            maxiter = 1e3;
            k = 0;
            
            % b is the best approximation to the zero point so far
            err = inf;
            b = d0;
            dphir = CO2.dphir_dd(d0,t);
            fb = d0 + d0^2*dphir - d_ig;
            while (abs(fb) > abstol || err > reltol) && k <= maxiter                
                % Newton step f(x) / f'(x)
                delta = fb / ( 1+2*b*dphir + b.^2*CO2.ddphir_dddd(b,t) );
                % Limit step
                err = abs(delta)/b;
                if err > 0.5
                    delta = 0.5*b*sign(delta);
                end
                b = b - delta;
                
                % Compute new function value
                dphir = CO2.dphir_dd(b,t);
                fb = b + b^2*dphir - d_ig;
                
                k = k+1;
            end                        
            
            % Return reduced density value
            d = b;
        end
        
        function val = dps_dT(T)
            % Derivative of saturation pressure with respect to temperature
            % According to Eq. (3.13)
            a1 = -7.0602087;
            a2 = 1.9391218;
            a3 = -1.6463597;
            a4 = -3.2995634;
            t1 = 1;
            t2 = 1.5;
            t3 = 2;
            t4 = 4;
            
            idx = CO2.Tt <= T & T <= CO2.Tc;
            p_s = CO2.VaporPressure(T(idx));
            T_ = 1 - T(idx)/CO2.Tc;
            
            val = nan(size(T));
            val(idx) = -p_s./T(idx).*( log(p_s/CO2.pc) + ...
                ( a1*t1*T_.^(t1-1) + a2*t2*T_.^(t2-1) + a3*t3*T_.^(t3-1) + a4*t4*T_.^(t4-1) ) ); 
        end
        
        function [valid, rho, T] = isValid_rhoT(rho, T)
            % Check if the (rho,T) couples are valid and return valid
            % couples.
            % (rho,T) couples are invalid if they are within the
            % two phase region. 
            % Solid state would also be invalid, but it is not checked.
            
            idx = CO2.Tt < T & T < CO2.Tc;
            rhoLiquid = CO2.SaturatedLiquidDensity(T(idx));
            rhoVapor = CO2.SaturatedVaporDensity(T(idx));
            
            if isscalar(rho)
                % rho is scalar; T may be of any size.
                valid = true(size(T));
                valid(idx) = ~(rhoVapor < rho & rho < rhoLiquid);
                rho = rho(any(valid));
                T = T(valid);
            else
                if isscalar(T)
                    if idx
                        % rho is of any size; T is a scalar and T < Tc
                        valid = ~(rhoVapor < rho & rho < rhoLiquid);
                        rho = rho(valid);
                        T = T(any(valid));
                    else
                        % rho is of any size; T is a scalar and T >= Tc
                        valid = true(size(rho));
                    end
                else
                    % rho and T must be of the same size
                    valid = true(size(rho));
                    valid(idx) = ~(rhoVapor < rho(idx) & rho(idx) < rhoLiquid);
                    rho = rho(valid);
                    T = T(valid);
                end
            end
        end
    end
    
    % Helmholz function
    properties(Constant, Hidden) 
        % Coefficients for the calculation of phi0 from Table 27
        phi0_a = [8.37304456 -3.70454304 2.5 1.99427042 0.62105248 0.41195293 1.04028922 0.08327678];
        phi0_theta = [0 0 0 3.15163 6.11190 6.77708 11.32384 27.08792];
        
        % Coefficients for the calculation of phir from Table 31
        phir_n = [...
             0.38856823203161e0   0.29385475942740e1  -0.55867188534934e1  -0.76753199592477e0   0.31729005580416e0   0.54803315897767e0   0.12279411220335e0 ... 1-7
             0.21658961543220e1   0.15841735109724e1  -0.23132705405503e0   0.58116916431436e-1 -0.55369137205382e0   0.48946615909422e0  -0.24275739843501e-1 ... 8-14
             0.62494790501678e-1 -0.12175860225246e0  -0.37055685270086e0  -0.16775879700426e-1 -0.11960736637987e0  -0.45619362508778e-1  0.35612789270346e-1 ... 15-21
            -0.74427727132052e-2 -0.17395704902432e-2 -0.21810121289527e-1  0.24332166559236e-1 -0.37440133423463e-1  0.14338715756878e0  -0.13491969083286e0 ... 22-28
            -0.23151225053480e-1  0.12363125492901e-1  0.21058321972940e-2 -0.33958519026368e-3  0.55993651771592e-2 -0.30335118055646e-3 ... 29-34
            -0.21365488688320e3   0.26641569149272e5  -0.24027212204557e5  -0.28341603423999e3   0.21247284400179e3... 35-39
            -0.66642276540751e0   0.72608632349897e0   0.55068668612842e-1 ... 40-42
            ];
        phir_d = [...
            1 1 1 1 2 2 3 ... 1-7
            1 2 4 5 5 5 6 ... 8-14
            6 6 1 1 4 4 4 ... 15-21
            7 8 2 3 3 5 5 ... 22-28
            6 7 8 10 4 8 ... 29-34
            2 2 2 3 3 ... 35-39
            ];
        phir_t = [...
            0.00 0.75 1.00 2.00 0.75 2.00 0.75 ... 1-7
            1.50 1.50 2.50 0.00 1.50 2.00 0.00 ... 8-14
            1.00 2.00 3.00 6.00 3.00 6.00 8.00 ... 15-21
            6.00 0.00 7.00 12.0 16.0 22.0 24.0 ... 22-28
            16.0 24.0 8.00 2.00 28.0 14.0 ... 29-34
            1.00 0.00 1.00 3.00 3.00 ... 35-39
            ];
        phir_c = [...
            1 1 1 1 1 1 1 ... 8-14
            1 1 2 2 2 2 2 ... 15-21
            2 2 3 3 3 4 4 ... 22-28
            4 4 4 4 5 6 ... 29-34
            ];
        
        phir_alpha = [25 25 25 15 20]; % 35-39
        phir_beta = [325 300 300 275 275]; % 35-39
        phir_gamma = [1.16 1.19 1.19 1.25 1.22]; % 35-39
        phir_epsilon = [1 1 1 1 1]; % 35-39
        
        phir_a = [3.5 3.5 3.0]; % 40-42
        phir_b = [0.875 0.925 0.875]; % 40-42
        phir_bb = [0.3 0.3 0.3];  % 40-42
        phir_A = [0.7 0.7 0.7]; % 40-42
        phir_B = [0.3 0.3 1]; % 40-42
        phir_C = [10 10 12.5]; % 40-42
        phir_D = [275 275 275]; % 40-42
    end

    methods(Static, Hidden) % phi0
        % Dimensionless helmholz function phi0 and its derivatives
        % according to Table 28
        % phi0 ... dimensionless helmholz function
        % d ... reduced density delta = rho / rhoc
        % t ... inverse reduced temperature tau = Tc / T
        function val = phi0(d, t)
            val = log(d) + CO2.phi0_a(1) + CO2.phi0_a(2)*t + CO2.phi0_a(3)*log(t) + ...
                sum( CO2.phi0_a(4:8) .* log( 1-exp(-CO2.phi0_theta(4:8).*t) ) ,2);
        end
        
        function val = dphi0_dd(d, ~)
            val = d.^-1;
        end
        
        function val = ddphi0_dddd(d, ~)
            val = -d.^-2;
        end
        
        function val = ddphi0_dddt(d, ~)
            val = zeros(size(d));
        end
        
        function val = dphi0_dt(~, t)
            e = exp(-CO2.phi0_theta(4:8).*t);
            val = CO2.phi0_a(2) + CO2.phi0_a(3)./t + ...
                sum( CO2.phi0_a(4:8).*CO2.phi0_theta(4:8).*( (1-e).^-1 - 1 ) ,2);
        end
        
        function val = ddphi0_dtdt(~, t)
            e = exp(-CO2.phi0_theta(4:8).*t);
            val = -CO2.phi0_a(3)./t.^2 - ...
                sum( CO2.phi0_a(4:8).*CO2.phi0_theta(4:8).^2.*e.*( 1-e ).^-2 ,2);
        end
    end
    
    methods(Static, Hidden) % phir
        % Residual part of the dimensionless helmholz function phir and its derivatives
        % according to Table 32
        % phir ... residual part of the dimensionless helmholz function
        % d ... reduced density delta = rho / rhoc
        % t ... inverse reduced temperature tau = Tc / T
        % Delta ... Distance function
        % Theta ... Some other function
        % Psi ... Exponential function
        
        function val = phir(d, t)
            val = sum( CO2.phir_n(1:7) .* d.^CO2.phir_d(1:7) .* t.^CO2.phir_t(1:7) ,2) + ...
                sum( CO2.phir_n(8:34) .* d.^CO2.phir_d(8:34) .* t.^CO2.phir_t(8:34) .* exp(-d.^CO2.phir_c) ,2) + ...
                sum( CO2.phir_n(35:39) .* d.^CO2.phir_d(35:39) .* t.^CO2.phir_t(35:39) .* ...
                    exp( -CO2.phir_alpha.*(d-CO2.phir_epsilon).^2 - CO2.phir_beta.*(t-CO2.phir_gamma).^2 ) ,2) + ...
                sum( CO2.phir_n(40:42) .* CO2.Delta(d, t).^CO2.phir_b .* d .* CO2.Psi(d, t), 2);
        end
        
        function val = dphir_dd(d, t)
            val = sum( CO2.phir_n(1:7) .* CO2.phir_d(1:7) .* d.^(CO2.phir_d(1:7)-1) .* t.^CO2.phir_t(1:7) ,2) + ...
                sum( CO2.phir_n(8:34).*exp(-d.^CO2.phir_c).*( d.^(CO2.phir_d(8:34)-1) .* t.^CO2.phir_t(8:34) .* ( CO2.phir_d(8:34) - CO2.phir_c.*d.^CO2.phir_c ) ) ,2) + ...
                sum( CO2.phir_n(35:39) .* d.^CO2.phir_d(35:39) .* t.^CO2.phir_t(35:39) .* ...
                    exp( -CO2.phir_alpha.*(d-CO2.phir_epsilon).^2 - CO2.phir_beta.*(t-CO2.phir_gamma).^2 ) .* ...
                    ( CO2.phir_d(35:39)./d - 2*CO2.phir_alpha.*(d-CO2.phir_epsilon) ) ,2) + ...
                sum( CO2.phir_n(40:42) .* ( CO2.Delta(d,t).^CO2.phir_b.*( CO2.Psi(d,t) + d.*CO2.dPsi_dd(d,t) ) + ...
                        CO2.dDeltab_dd(d,t) .* d .* CO2.Psi(d,t) ) ,2);
        end
        
        function val = ddphir_dddd(d, t)
            val = sum( CO2.phir_n(1:7) .* CO2.phir_d(1:7) .* (CO2.phir_d(1:7)-1) .* d.^(CO2.phir_d(1:7)-2) .* t.^CO2.phir_t(1:7) ,2) + ...
                sum( CO2.phir_n(8:34).*exp(-d.^CO2.phir_c).*( d.^(CO2.phir_d(8:34)-2) .* t.^CO2.phir_t(8:34) .* ( ...
                    ( CO2.phir_d(8:34) - CO2.phir_c.*d.^CO2.phir_c ) .* ( CO2.phir_d(8:34) - 1 - CO2.phir_c.*d.^CO2.phir_c ) - ...
                    CO2.phir_c.^2.*d.^CO2.phir_c )) ,2) + ...
                sum( CO2.phir_n(35:39) .* t.^CO2.phir_t(35:39) .* exp( -CO2.phir_alpha.*(d-CO2.phir_epsilon).^2 - CO2.phir_beta.*(t-CO2.phir_gamma).^2 ) .* ...
                    ( -2*CO2.phir_alpha.*d.^CO2.phir_d(35:39) + 4*CO2.phir_alpha.^2.*d.^CO2.phir_d(35:39).*( d-CO2.phir_epsilon ).^2 - ...
                    4*CO2.phir_d(35:39).*CO2.phir_alpha.*d.^( CO2.phir_d(35:39)-1 ).*( d-CO2.phir_epsilon ) + CO2.phir_d(35:39).*( CO2.phir_d(35:39)-1 ).*d.^( CO2.phir_d(35:39)-2 ) ) ,2) + ...
                sum( CO2.phir_n(40:42) .* ( CO2.Delta(d,t).^CO2.phir_b.*( 2*CO2.dPsi_dd(d,t) + d.*CO2.ddPsi_dddd(d,t) ) + ...
                    2*CO2.dDeltab_dd(d,t).*( CO2.Psi(d,t) + d.*CO2.dPsi_dd(d,t) ) + CO2.ddDeltab_dddd(d,t).*d.*CO2.Psi(d,t) ) ,2);
        end
        
        function val = dphir_dt(d, t)
            val = sum( CO2.phir_n(1:7) .* CO2.phir_t(1:7) .* d.^CO2.phir_d(1:7) .* t.^( CO2.phir_t(1:7)-1 ) ,2) + ...
                sum( CO2.phir_n(8:34) .* CO2.phir_t(8:34) .* d.^CO2.phir_d(8:34) .* t.^( CO2.phir_t(8:34)-1 ) .* exp( -d.^CO2.phir_c ) ,2) + ...
                sum( CO2.phir_n(35:39) .* d.^CO2.phir_d(35:39) .* t.^CO2.phir_t(35:39) .* ...
                    exp( -CO2.phir_alpha.*(d-CO2.phir_epsilon).^2 - CO2.phir_beta.*(t-CO2.phir_gamma).^2 ) .* ...
                    ( CO2.phir_t(35:39)./t - 2*CO2.phir_beta.*(t-CO2.phir_gamma) ) ,2) + ...
                sum( CO2.phir_n(40:42).*d.*( CO2.dDeltab_dt(d,t).*CO2.Psi(d,t) + CO2.Delta(d,t).^CO2.phir_b.*CO2.dPsi_dt(d,t) ) ,2);
        end
        
        function val = ddphir_dtdt(d, t)
            val = sum( CO2.phir_n(1:7) .* CO2.phir_t(1:7) .* ( CO2.phir_t(1:7)-1 ) .* d.^CO2.phir_d(1:7) .* t.^( CO2.phir_t(1:7)-2 ) ,2) + ...
                sum( CO2.phir_n(8:34) .* CO2.phir_t(8:34) .* ( CO2.phir_t(8:34)-1 ) .* d.^CO2.phir_d(8:34) .* t.^( CO2.phir_t(8:34)-2 ) .* exp(-d.^CO2.phir_c) ,2) + ...
                sum( CO2.phir_n(35:39) .* d.^CO2.phir_d(35:39) .* t.^CO2.phir_t(35:39) .* ...
                    exp( -CO2.phir_alpha.*(d-CO2.phir_epsilon).^2 - CO2.phir_beta.*(t-CO2.phir_gamma).^2 ) .* ...
                    ( ( CO2.phir_t(35:39)./t - 2*CO2.phir_beta.*(t-CO2.phir_gamma) ).^2 - CO2.phir_t(35:39)./t.^2 - 2*CO2.phir_beta ) ,2)  + ...
                sum( CO2.phir_n(40:42).*d.*( CO2.ddDeltab_dtdt(d,t).*CO2.Psi(d,t) + 2*CO2.dDeltab_dt(d,t).*CO2.dPsi_dt(d,t) + CO2.Delta(d,t).^CO2.phir_b.*CO2.ddPsi_dtdt(d,t) ) ,2);
        end
        
        function val = ddphir_dddt(d, t)
            val = sum( CO2.phir_n(1:7) .* CO2.phir_d(1:7) .* CO2.phir_t(1:7) .* d.^( CO2.phir_d(1:7)-1 ) .* t.^( CO2.phir_t(1:7)-1 ) ,2) + ...
                sum( CO2.phir_n(8:34).*exp(-d.^CO2.phir_c).*d.^( CO2.phir_d(8:34)-1 ).*CO2.phir_t(8:34).* t.^( CO2.phir_t(8:34)-1 ).*( CO2.phir_d(8:34) - CO2.phir_c.*d.^CO2.phir_c ) ,2) + ...
                sum( CO2.phir_n(35:39) .* d.^CO2.phir_d(35:39) .* t.^CO2.phir_t(35:39) .* ...
                    exp( -CO2.phir_alpha.*(d-CO2.phir_epsilon).^2 - CO2.phir_beta.*(t-CO2.phir_gamma).^2 ) .* ...
                    ( CO2.phir_d(35:39)./d - 2*CO2.phir_alpha.*(d-CO2.phir_epsilon) ) .* ...
                    ( CO2.phir_t(35:39)./t - 2*CO2.phir_beta.*(t-CO2.phir_gamma) ) ,2) + ...
                sum( CO2.phir_n(40:42) .* (  CO2.Delta(d,t).^CO2.phir_b.*( CO2.dPsi_dt(d,t) + d.*CO2.ddPsi_dddt(d,t) ) + ...
                    d.*CO2.dDeltab_dd(d,t).*CO2.dPsi_dt(d,t) + CO2.dDeltab_dt(d,t).*( CO2.Psi(d,t) + d.*CO2.dPsi_dd(d,t) ) + ...
                    CO2.ddDeltab_dddt(d,t).*d.*CO2.Psi(d,t) ) ,2);
        end
        
        function val = Delta(d, t)
            val = CO2.Theta(d, t).^2 + CO2.phir_B .* ( (d-1).^2 ).^CO2.phir_a;
        end
        
        function val = Theta(d, t)
            val = (1-t) + CO2.phir_A .* ( (d-1).^2 ).^( 1./(2*CO2.phir_bb) );
        end
        
        function val = dDeltab_dd(d, t)
            val = CO2.phir_b .* CO2.Delta(d, t).^(CO2.phir_b-1) .* CO2.dDelta_dd(d, t);
            val(isnan(val)) = 0;
        end
        
        function val = ddDeltab_dddd(d, t)
            val = CO2.phir_b .* ( ...
                CO2.Delta(d, t).^(CO2.phir_b-1) .* CO2.ddDelta_dddd(d, t) + ...
                (CO2.phir_b-1) .* CO2.Delta(d, t).^(CO2.phir_b-2) .* CO2.dDelta_dd(d, t).^2 ...
                );
            val(isnan(val)) = 0;
        end
        
        function val = dDeltab_dt(d, t)
            val = -2 .* CO2.Theta(d,t) .* CO2.phir_b .* CO2.Delta(d,t).^( CO2.phir_b-1 );
            val(isnan(val)) = 0;
        end
        
        function val = ddDeltab_dtdt(d, t)
            val = 2 .* CO2.phir_b .* CO2.Delta(d,t).^( CO2.phir_b-1 ) + ...
                4 .* CO2.Theta(d, t).^2 .* CO2.phir_b .* (CO2.phir_b-1) .* CO2.Delta(d, t).^(CO2.phir_b-2);
            val(isnan(val)) = 0;
        end
        
        function val = ddDeltab_dddt(d, t)
            val = -CO2.phir_A .* CO2.phir_b .* 2./CO2.phir_bb .* CO2.Delta(d, t).^(CO2.phir_b-1) .* (d-1) .* ( (d-1).^2 ).^( 1./(2*CO2.phir_bb) - 1 ) - ...
                2 * CO2.Theta(d, t) .* CO2.phir_b .* (CO2.phir_b-1) .* CO2.Delta(d, t).^(CO2.phir_b-2) .* CO2.dDelta_dd(d, t);
            val(isnan(val)) = 0;
        end
        
        function val = dDelta_dd(d, t)
            val = (d-1) .* ( ...
                CO2.phir_A .* CO2.Theta(d, t) .* 2./CO2.phir_bb .* ( (d-1).^2 ).^( 1./(2*CO2.phir_bb) - 1) + ...
                2 * CO2.phir_B .* CO2.phir_a .* ( (d-1).^2 ).^(CO2.phir_a-1) ...
                );
        end
        
        function val = ddDelta_dddd(d, t)
            val = 1./(d-1) .* CO2.dDelta_dd(d, t) + (d-1).^2 .* ( ...
                4 * CO2.phir_B .* CO2.phir_a .* ( CO2.phir_a-1 ) .* ( (d-1).^2 ).^( CO2.phir_a-2 ) + ...
                2 * CO2.phir_A.^2 .* CO2.phir_bb.^-2 .* (( (d-1).^2 ).^( 1./(2*CO2.phir_bb) - 1 )).^2 + ...
                CO2.phir_A .* CO2.Theta(d, t) .* 4 ./ CO2.phir_bb .* (1./(2.*CO2.phir_bb) -1) .* ( (d-1).^2 ).^( 1./(2*CO2.phir_bb) - 2 ) ...
                );
        end
        
        function val = Psi(d, t)
            val = exp( -CO2.phir_C.*(d-1).^2 - CO2.phir_D.*(t-1).^2 );
        end
        
        function val = dPsi_dd(d, t)
            val = -2 .* CO2.phir_C .* (d-1) .* CO2.Psi(d, t);
        end
                
        function val = ddPsi_dddd(d, t)
            val = ( 2.*CO2.phir_C.*(d-1).^2 - 1).*2.*CO2.phir_C.*CO2.Psi(d, t);
        end
                
        function val = dPsi_dt(d, t)
            val = -2 .* CO2.phir_D .* (t-1) .* CO2.Psi(d, t);
        end
                
        function val = ddPsi_dtdt(d, t)
            val = ( 2.*CO2.phir_D.*(t-1).^2 - 1).*2.*CO2.phir_D.*CO2.Psi(d, t);
        end
                
        function val = ddPsi_dddt(d, t)
            val = 4.*CO2.phir_C.*CO2.phir_D.*(d-1).*(t-1).*CO2.Psi(d, t);
        end
    end
end

