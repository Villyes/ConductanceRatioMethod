classdef CO2
    % CO2 Calculate thermodynamic equlibrium and transport properties of
    %   carbon dioxide in the fluid region.
    %
    % The functions in this class are based on scientific publications,
    % which are listed a bit further down. If you use this code, make sure
    % to cite these publications accordingly. These papers also contain
    % further information on equations, such as the validity range. 
    %
    %
    % Phase boundaries and virial coefficients as a function of temperature
    % in K
    %   melting pressure pMelt in MPa
    %   sublimation pressure pSub in MPa
    %   vapor pressure pVap in MPa
    %   saturated liquid density rhoLiqSat in kg/m^3
    %   saturated vapor density rhoVapSat in kg/m^3
    %   second virial coefficient B in m^3/kg
    %   third virial coefficient C in (m^3/kg)^2 
    % Phase boundary functions return nan outside their validity range. The
    % validity range is not checked for the virial coefficients.
    %
    % Equilibrium properties can be calculated in two ways: 
    % 1) as a function of density in kg/m^3 and temperature in K
    %   e.g. entropy = CO2.s_rhoT(density, temperature)
    % 2) as a function of pressure in MPa and temperature in K
    %   e.g. entropy = CO2.s_pT(pressure, temperature)
    % 3) as a function of pressure in Pa and temperature in K
    %   e.g. entropy = CO2.s(pressure, temperature)
    %   to use with hex() 
    % The inputs must have the same size or one of them has to be a scalar.
    %
    % Equilibrium properties:
    %   pressure p in MPa
    %   density rho in kg/m^3
    %   entropy s in kJ/(kg K) 
    %   internal energy u in kJ/kg
    %   isochoric heat capacity cv in kJ/(kg K)
    %   enthalpy h in kJ/kg 
    %   isobaric heat capacity cp in kJ/(kg K) 
    %   saturated liquid heat capacity cs in kJ/(kg K) 
    %   speed of sound w in m/s 
    %   Joule-Thomson coefficient mu in K/MPa 
    %   fugacity coefficient f
    %
    % Transport properties have to be calculated both at once, since they
    % depend on each other in the critical region. Just as with equilibrium
    % properties, they can be calculated in two ways: 
    % 1) as function of density in  kg/m3 and temperature in K
    %   [mu, lambda] = CO2.transport_rhoT(density, temperature)
    % 2) as a function of pressure in MPa and temperature in K
    %   [mu, lambda] = CO2.transport_pT(pressure, temperature)
    %
    % Transport properties:
    %   viscosity mu in mPa s
    %   thermal conductivity lambda in mW/(m K)
    %
    %
    % 
    % +++ Implementation notes +++
    % 
    % The validity range of the equations is not checked automatically.
    % At low temperatures/high pressures, you might have to check for the
    % solid region.
    % Only the .*_rhoT functions identify values in the two-phase region and
    % return nan. If this check is not required the internal functions
    % .*_rhoT_i can be used. The internal functions will return arbitrary
    % values in the two-phase region and they always return column vectors
    % regarless of the input dimensions.
    %
    % The .*_pT functions first compute the density with an iterative
    % algorithm and then the requested property with the corresponding
    % .*_rhoT function. If more than one property is computed at the same
    % p-T conditions, it is more efficient to compute the density
    % explicitly with the rho_pT function and then use the .*_rhoT 
    % functions.
    %
    % The calculation of the transport properties uses values from the
    % equation of state by Span et al. (1994). Unfortunately, this EOS
    % shows some unphysical behavior in close vicinty of the critical
    % point. Specifically, there is a "dent" in the critical peak (see 
    % Figure 15 in Span et al. (1994) and also Figure 14 and 16 in Huber et
    % al. (2016)). 
    % Huber et al. propose to use the scaled equation by Albright et al. in
    % the region 303.1282 < T/K < 305.1282 and 350 < rho/(kg/m) < 530 if 
    % highly accurate results are required.
    %
    % For the transport properties two variants are implemented. The 
    % up-to-date one, based on the work of Huber et al. (2016), Laesecke et
    % al. (2017) and Luettmer-Strahtmann (1995), which should be equivalent
    % to the implementation in Refprop 10. The depreciated one, based on 
    % Vesovic et al. (1990) and Fenghour et. al (1997), which should be 
    % equivalent to Refprop 7.
    % More information is included in the comments of the transport_rhoT_i
    % function.
    %
    %
    %
    % +++ References +++
    %
    % This class is based on the following publications. If you use it, 
    % make sure that you cite these papers.
    %
    % Equation of State/Equilibrium properties
    %   R. Span and W. Wagner (1994): A New Equation of State for Carbon
    %   Dioxide Covering the Fluid Region from the Triple-Point Temperature
    %   to 1100 K at Pressures up to 800 MPa. In J. Phys. Chem. Ref. Data,
    %   Vol. 25, No. 6
    %   DOI: https://doi.org/10.1063/1.555991
    %
    % Transport properties
    %   M. L. Huber, E. A. Sykioti, M. J. Assael and R. A. Perkins (2016): 
    %   Reference Correlation of the Thermal Conductivity of Carbon Dioxide
    %   from the Triple Point to 1100 K and up to 200 MPa. In J. Phys.
    %   Chem. Ref. Data, Vol. 45, No. 1
    %   DOI: https://doi.org/10.1063/1.4940892
    %
    %   A. Laesecke and C. D. Muzny (2017): Reference Correlation for the
    %   Viscosity of Carbon Dioxide. In J. Phys. Chem. Ref. Data, Vol. 46,
    %   No. 1
    %   DOI: https://doi.org/10.1063/1.4977429
    %
    %   J. Luettmer-Strathmann, J. V. Sengers and G. A. Olchowy (1995):
    %   Non-asymptotic critical behavior of the transport properties of 
    %   fluids. In J. Chem. Phys. 103 (17)
    %   DOI: https://doi.org/10.1063/1.470718
    %
    % Transport properties (depreciated)
    %   V. Vesovic, W. A. Wakeham, G. A. Olchowy, J. V. Sengers, J. T. R.
    %   Watson and J. Millat  (1990): The Transport Properties of Carbon
    %   Dioxide. In J. Phys. Chem. Ref. Data, Vol. 19, No. 3
    %   DOI: https://doi.org/10.1063/1.555875
    %
    %   A. Fenghour, W. A. Wakeham and V. Vesovic (1997): The Viscosity of 
    %   Carbon Dioxide. In J. Phys. Chem. Ref. Data, Vol. 27, No. 1
    %   DOI: https://doi.org/10.1063/1.556013
    %
    %
    %
    % +++ Disclaimer +++
    %
    % I did my best to test and validate this code, but I cannot guarantee
    % that there are no errors in here.
    %
    % THE SOFTWARE IS PROVIDED AS IS, WITHOUT WARRANTY OF ANY KIND, 
    % EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF 
    % MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND 
    % NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR ANY 
    % CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, 
    % TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE 
    % SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.    
    %
    %
    %   (C) Felix Birkelbach 2020, CC BY-SA 4.0
    %   You are free to
    %       + copy and redistribute the material in any medium or format 
    %       + remix, transform, and build upon the material for any
    %         purpose, even commercially
    %   under the following terms:
    %       + You must give appropriate credit, provide a link to the 
    %         license, and indicate if changes were made. You may do so in 
    %         any reasonable manner, but not in any way that suggests the 
    %         licensor endorses you or your use.
    %       + If you remix, transform, or build upon the material, you must
    %         distribute your contributions under the same license as the 
    %         original. 
    %   https://creativecommons.org/licenses/by-sa/4.0/
    %
    %   Contributions by Theresa Brunauer
    %   Contributions by Viktoria Illyes

    properties(Constant)
        M = 44.0098; % g/mol molar mass
        RM = 8.314510; % J/(mol K) molar gas constant
        R = 0.1889241; % kJ/(kg K) specific gas constant
        
        % Triple point
        Tt = 216.592; % K triple point temperature
        pt = 0.51795; % MPa triple point pressure
        rhotLiquid = 1178.53; % kg/m3 liquid density at triple point
        rhotVapor = 13.7614; % kg/m3 vapor density at triple point
        
        % Critical point
        Tc = 304.1282; % K critical temperature
        pc = 7.3773; % MPa critical pressure
        rhoc = 467.6; % kg/m3 density at critical point
        
        T0 = 298.15; % K reference temperature
        p0 = 0.1013258 % MPa reference pressure
        h0 = 0; % kJ/kg reference enthalpy in the ideal gas state at T0
        s0 = 0; % kJ/kg K reference entropy in the ideal gas state at T0, p0
        
        % In the NIST webbook an other reference state, the IIR convention,
        % is used: 
        %   h = 200 kJ/kg at 0C for saturated liquid
        %   s = 1 J/gK at 0C for saturated liquid
        % In most application it makes no difference, because only 
        % differences are of interest.
        % To convert to the IRR convention use
        %   h_irr = h + CO2.delta_h
        %   s_irr = s + CO2.delta_s
        %   u_irr = u + CO2.delta_h
        delta_h = 506.7770; % kJ/kg offset to IRR reference state
        delta_s = 2.7390; % kJ/kg K offset to IRR reference state
    end


    properties(Constant)
        cf = 1.0001;    %properties of the fluid                                                    
        q0 = 150000;  %substance-specific values
        d0 = 10^-2 ;
        Ra0 = 10^-6;
        alpha_0 = 18890; %W/m^2K
        p_c = 73.8 * 10^5 %Pa
    end

    methods(Static)
        function h = h_rhoT(rho,T)
            h = CO2.h_rho_T(rho, T).*10^3;
        end

        function Pr=Pr(rho,T)   %Prandtl-number
            %rho in kg/m³, T in K
            c_p = CO2.cp_rhoT(rho, T).*10^3;
            my = CO2.transport_rhoT(rho, T).*10^-3;
            lambda = CO2.lambda_in(rho,T).*10^-3;
            Pr=c_p.*my./lambda;
        end

        function lambda = CO2.lambda(rho,T)
            lambda = CO2.lambda_in(rho,T).*10^-3;
        end

        function is2phase=is2phase(rho,T)
            is2phase = ~CO2.isSinglePhase_rhoT(rho, T);
        end 

        function rho=rho_satL(T)
            rho = CO2.rhoLiqSat(T);
        end

        function rho = rho_satV(T)
            rho = CO2.rhoVapSat(T);
        end

        function x = x(rho,T)
            rho0 = CO2.rho_satL(T);
            rho1 = CO2.rho_satV(T);
            x = CO2.x_intern(rho,rho0,rho1);
        end

        function x=x_intern(rho,rho0,rho1)
            x=(rho.^-1-rho0.^-1)./(rho1.^-1-rho0.^-1);
        end

    end

   methods(Static) % Properties in (p, T) with (Pa, K) for hex()-class
   % Thermodynamic properties as functions of pressure in Pa and 
   % temperature in K.    

        function h = h(p,T)
            % Enthalpy in J/kg as function of pressure in Pa and temperature in K
            p = p .* 10^-6;
            h = CO2.h_pT(p, T).*10^3;
        end

        function rho = rho(p, T)
            % Density in kg/m3 as function of pressure in Pa and temperature in K
            % This function is used in all p-T functions to compute the
            % density. Then the requested property is computed with the
            % rho,T functions.
            p=p.*1e-6;
            [p, T, sz] = CO2.checkInput(p, T);
            if isscalar(T)
                T = repmat(T, size(p));
            end
            
            isLiquid = p > CO2.pVap(T(:)); % no need to check T < CO2.Tc 
            d_ig = p.*1e3/CO2.R./T/CO2.rhoc; % reduced density of ideal gas
            
            % initial values for iteration
            d0 = nan(size(T));
            % estimate density in liquid with saturated liquid density at T
            d0(isLiquid) = CO2.rhoLiqSat(T(isLiquid)) /CO2.rhoc;
            % estimate density in gas/sc-fluid with density of ideal gas
            d0(~isLiquid) = d_ig(~isLiquid);
            
            % Compute reduced density with inverse function
            t = CO2.Tc ./ T(:);
            d = CO2.d_newton(t, d_ig, d0);

            rho = reshape(d.*CO2.rhoc, sz);
        end  

        function v = v(p, T)
            % Specific volume in m3/kg as function of pressure in Pa and temperature in K
            p=p.*1e6;
            [p, T, sz] = CO2.checkInput(p, T);
            if isscalar(T)
                T = repmat(T, size(p));
            end
            
            isLiquid = p > CO2.pVap(T(:)); % no need to check T < CO2.Tc 
            d_ig = p.*1e3/CO2.R./T/CO2.rhoc; % reduced density of ideal gas
            
            % initial values for iteration
            d0 = nan(size(T));
            % estimate density in liquid with saturated liquid density at T
            d0(isLiquid) = CO2.rhoLiqSat(T(isLiquid)) /CO2.rhoc;
            % estimate density in gas/sc-fluid with density of ideal gas
            d0(~isLiquid) = d_ig(~isLiquid);
            
            % Compute reduced density with inverse function
            t = CO2.Tc ./ T(:);
            d = CO2.d_newton(t, d_ig, d0);

            rho = reshape(d.*CO2.rhoc, sz);
            v=1./rho;
        end

        function T=T_ph(p,h)
            % Temperature in K as function of pressure in Pa and enthalpy in J/kg
            p = p.*10^-6;
            h = h .* 10^-3;
        
            T = NaN(1,size(h,2));
            for i=1:size(h,2)
                T(i) = fzero(@(T) CO2.h_pT(p(i),T)-h(i),[218, 1500]);
            end
        end

        function w=w(p,T,x)
            if nargin<3 || isempty(x)
                x=NaN;
            end
            p = p.*10^-6;
            rho = CO2.rho_pT(p,T);
            w = CO2.w_rhoT(rho,T);
        end

        function p=p(rho,T) 
        % Pressure in Pa as function of temperature in K and density in
        % kg/m^3
            p = CO2.p_rhoT(rho, T).*10^6; 
        end

        function my=my(rho,T) 
            % Viscosity eta in Pa s
            [MY, ~] = CO2.transport_rhoT(rho, T);
            my = MY.*10^-3;
        end
        
        function lambda = lambda(rho,T)
            % Thermal conductivity lambda in W/m K
            lambda = CO2.lambda_in(rho,T).*10^-3;
        end

        function cp = cp(p, T)
            % Isobaric heat capacity in J/(kg K) as function of pressure in Pa and temperature in K
            p=p./1e6;
            cp=CO2.cp_pT(p, T).*1000;
        end

        function Tpc=Tpc(p)
            % pseudocritical temperature in K as a function of pressure in Pa
            Tpc = 3.700621431891665e+02.*exp(-2.684598659785788e-10.*p) -1.345178030313569e+02.*exp(-9.751799518360223e-08.*p);

            % n=1000;
            % 
            % p1=repmat((linspace(75,300,n)/10)',1,n); %MPa
            % T1=repmat(linspace(31,120,n)+273.15,n,1);
            % cp1=CO2.cp_pT(p1,T1);
            % 
            % [~,idx]=max(cp1,[],2);
            % 
            % figure(2)
            % p=p1(:,1)*10^6;
            % T=T1(1,idx);
            % plot(p,T)
            % Fit Name: untitled fit 1
            % 
            % Exponential Curve Fit (exp2)
            % f(x) = a*exp(b*x) + c*exp(d*x)
            % 
            % Coefficients and 95% Confidence Bounds
            %        Value        Lower        Upper   
            % a    370.0621     369.4163     370.7080 
            % b    -0.0000      -0.0000      -0.0000  
            % c    -134.5178    -134.9099    -134.1257
            % d    -0.0000      -0.0000      -0.0000  
            % 
            % Goodness of Fit
            %              Value   
            % SSE         7.1792  
            % R-square    1.0000  
            % DFE         996.0000
            % Adj R-sq    1.0000  
            % RMSE        0.0849  
            %
            %coeffvalues(fittedmodel)
        end

   end

    methods(Static) % Phase boundaries and properties in T
        function pM = pMelt(T)
            % Melting pressure as funktion of temperature in K
            % according to Eq. (3.10)
            a1 = 1955.5390;
            a2 = 2055.4593;
            
            pM = nan(size(T));
            idx = T >= CO2.Tt;
            
            T_ = T(idx)/CO2.Tt - 1;
            pM(idx) = CO2.pt .* ( 1 + a1.*T_ + a2.*T_.^2 );
        end
        
        function pS = pSub(T)
            % Sublimation pressure as funktion of temperature in K
            % according to Eq. (3.12)
            a1 = -14.740846;
            a2 = 2.4327015;
            a3 = -5.3061778;
            
            pS = nan(size(T));
            idx = T <= CO2.Tt;
            
            T_ = 1 - T(idx)/CO2.Tt;
            pS(idx) = CO2.pt .* exp( CO2.Tt ./ T(idx) .*  ...
                ( a1.*T_ + a2.*T_.^1.9 + a3.*T_.^2.9 ) );
        end
        
        function pV = pVap(T)
            % Vapor pressure as function of temperature in K
            % according to Eq. (3.13)
            a1 = -7.0602087;
            a2 = 1.9391218;
            a3 = -1.6463597;
            a4 = -3.2995634;
            t1 = 1;
            t2 = 1.5;
            t3 = 2;
            t4 = 4;
            
            pV = nan(size(T));
            idx = CO2.Tt <= T & T <= CO2.Tc;
            
            T_ = 1 - T(idx)/CO2.Tc;
            pV(idx) = CO2.pc .* exp( CO2.Tc ./ T(idx) .* ...
                (a1.*T_.^t1 + a2.*T_.^t2 + a3.*T_.^t3 + a4.*T_.^t4 ) ); 
        end
        
        function rhoLiquid = rhoLiqSat(T)
            % Density of saturated liquid as function of temperature in K
            % according to Eq. (3.14)
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
            rhoLiquid(idx) = CO2.rhoc .* exp( ...
                (a1.*T_.^t1 + a2.*T_.^t2 + a3.*T_.^t3 + a4.*T_.^t4 ) ); 
        end
        
        function rhoVapor = rhoVapSat(T)
            % Density of saturated vapor as function of temperature in K
            % according to Eq. (3.15)
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
            rhoVapor(idx) = CO2.rhoc .* exp( ...
                (a1.*T_.^t1 + a2.*T_.^t2 + a3.*T_.^t3 + a4.*T_.^t4 + a5.*T_.^t5 ) ); 
        end
                
        function plotPhaseDiagram()
            newplot(); hold on
            
            T = [linspace(CO2.Tt-50, CO2.Tt, 50) linspace(CO2.Tt, 1100, 500)];
            TC = T - 273.15;
            plot(TC, CO2.pSub(T).*10, 'DisplayName', 'sublimation');
            plot(TC, CO2.pMelt(T).*10, 'DisplayName', 'melting');
            plot(TC, CO2.pVap(T).*10, 'DisplayName', 'vapor');
            plot(CO2.Tt-273.15, CO2.pt.*10, '.*', 'DisplayName', 'triple point')
            plot(CO2.Tc-273.15, CO2.pc.*10, '.*', 'DisplayName', 'critical point')
            legend('show')
            plot([CO2.Tc CO2.Tc]-273.15, [CO2.pc CO2.pMelt(CO2.Tc)].*10, 'k--',  'DisplayName', 'liquid-super critical');
            plot([CO2.Tc T(end)]-273.15, [CO2.pc CO2.pc].*10, 'k--', 'DisplayName', 'gas-super critical');
            ylabel('p in bar');
            xlabel('T in C');
            set(gca, 'yscale', 'log');
        end
        
        function b = B(T)
            % Second virial coefficient in m^3/kg as function of temperature in K
            % according to Table 3.
            t = CO2.Tc ./ T(:);
            b = CO2.dphir_dd(1e-18,t) / CO2.rhoc;
        end

        function c = C(T)
            % Third virial coefficient in (m^3/kg)^2 as function of temperature in K
            % according to Table 3.
            t = CO2.Tc ./ T(:);
            c = CO2.ddphir_dddd(1e-18,t) / CO2.rhoc.^2;
        end
    end
    
    methods(Static) % Properties in (rho, T) 
        % Equilibrium properties as functions of density in kg/m and
        % temperature in K according to Table 3 in Span et al. (1994).
        % A validity check is performed. For points in the two-phase region
        % nan is returned. If the validity check is not required the
        % internal functions with the suffix rhoT_i can be used instead.
        % The solid region is not checked.
        
        function p = p_rhoT(rho, T)
            % Pressure in MPa as function of density in kg/m3 and temperature in K
            [v, rho, T] = CO2.isSinglePhase_rhoT(rho, T);
            p = nan(size(v));
            p(v) = CO2.p_rhoT_i(rho, T);
        end
        
        function s = s_rhoT(rho, T)
            % Entropy in kJ/(kg K) as function of density in kg/m3 and temperature in K
            [v, rho, T] = CO2.isSinglePhase_rhoT(rho, T);
            s = nan(size(v));
            s(v) = CO2.s_rhoT_i(rho, T);
        end
        
        function u = u_rhoT(rho, T)
            % Internal energy in kJ/kg as function of density in kg/m3 and temperature in K
            [v, rho, T] = CO2.isSinglePhase_rhoT(rho, T);
            u = nan(size(v));
            u(v) = CO2.u_rhoT_i(rho, T);
        end
        
        function cv = cv_rhoT(rho, T)
            % Isochoric heat capacity in kJ/(kg K) as function of density in kg/m3 and temperature in K
            [v, rho, T] = CO2.isSinglePhase_rhoT(rho, T);
            cv = nan(size(v));
            cv(v) = CO2.cv_rhoT_i(rho, T);
        end        
 
        function h = h_rho_T(rho, T)
            % Enthalpy in kJ/kg as function of density in kg/m3 and temperature in K
            [v, rho, T] = CO2.isSinglePhase_rhoT(rho, T);
            h = nan(size(v));
            h(v) = CO2.h_rhoT_i(rho, T);
        end

        function cp = cp_rhoT(rho, T)
            % Isobaric heat capacity in kJ/(kg K) as function of density in kg/m3 and temperature in K
            [v, rho, T] = CO2.isSinglePhase_rhoT(rho, T);
            cp = nan(size(v));
            cp(v) = CO2.cp_rhoT_i(rho, T);
        end        

        function cs = cs_rhoT(rho, T)
            % Saturated liquid heat capacity in kJ/(kg K) as function of density in kg/m3 and temperature in K
            [v, rho, T] = CO2.isSinglePhase_rhoT(rho, T);
            cs = nan(size(v));
            cs(v) = CO2.cs_rhoT_i(rho, T);
        end
        
        function w = w_rhoT(rho, T)
            % Speed of sound in m/s as function of density in kg/m3 and temperature in K
            [v, rho, T] = CO2.isSinglePhase_rhoT(rho, T);
            w = nan(size(v));
            w(v) = CO2.w_rhoT_i(rho, T);
        end
        
        function mu = mu_rhoT(rho, T)
            % Joule-Thomson coefficient in K/MPa as function of density in kg/m3 and temperature in K
            [v, rho, T] = CO2.isSinglePhase_rhoT(rho, T);
            mu = nan(size(v));
            mu(v) = CO2.mu_rhoT_i(rho, T);
        end
        
        function f = f_rhoT(rho, T)
            % Fugacity coefficient as function of density in kg/m3 and temperature in K
            [v, rho, T] = CO2.isSinglePhase_rhoT(rho, T);
            f = nan(size(v));
            f(v) = CO2.f_rhoT_i(rho, T);
        end
    end
    
    methods(Static) % Properties in (p, T)
        % Thermodynamic properties as functions of pressure in MPa and 
        % temperature in K.

        function rho = rho_pT(p, T)
            % Density in kg/m3 as function of pressure in MPa and temperature in K
            % This function is used in all p-T functions to compute the
            % density. Then the requested property is computed with the
            % rho,T functions.
            [p, T, sz] = CO2.checkInput(p, T);
            if isscalar(T)
                T = repmat(T, size(p));
            end
            
            isLiquid = p > CO2.pVap(T(:)); % no need to check T < CO2.Tc 
            d_ig = p.*1e3/CO2.R./T/CO2.rhoc; % reduced density of ideal gas
            
            % initial values for iteration
            d0 = nan(size(T));
            % estimate density in liquid with saturated liquid density at T
            d0(isLiquid) = CO2.rhoLiqSat(T(isLiquid)) /CO2.rhoc;
            % estimate density in gas/sc-fluid with density of ideal gas
            d0(~isLiquid) = d_ig(~isLiquid);
            
            % Compute reduced density with inverse function
            t = CO2.Tc ./ T(:);
            d = CO2.d_newton(t, d_ig, d0);

            rho = reshape(d.*CO2.rhoc, sz);
        end        
        
        function rho = rho_pT_sequential(p, T)
            % Density in kg/m3 as function of pressure in MPa and temperature in K
            % It assumes that the (p,T) vectors are in order - i.e. that
            % the change of density from one point to the other will be
            % small (except for phase changes) so that each value can be
            % used as initial value for the subsequent iteration.
            [p, T, sz] = CO2.checkInput(p, T);
            if isscalar(T)
                T = repmat(T, size(p));
            end
            % Reduced density of the ideal gas
            d_ig = p(:).*1e3/CO2.R./T(:) / CO2.rhoc;
            p_v = CO2.pVap(T(:));
%             d_l = CO2.rhoLiqSat(T(:))/CO2.rhoc;
%             d_v = CO2.rhoVapSat(T(:))/CO2.rhoc;
            isLiquid = p > p_v; % no need to check T < CO2.Tc

            % Compute reduced density with inverse function
            t = CO2.Tc ./ T(:);
            d = zeros(size(d_ig));
            for i = 1:numel(d)
                % Reset starting value if there was a gas/liquid change or
                % if the last result was invalid
                if i == 1 || (T(i-1) < CO2.Tc && T(i) < CO2.Tc && isLiquid(i)~=isLiquid(i-1)) || ~isfinite(d(i-1))
                    % Use new starting value
                    if isLiquid(i)
                        d0 = CO2.rhoLiqSat(T(i))/CO2.rhoc;
                    else
                        d0 = d_ig(i);
                    end
                else
                    % Use result of last step as starting value
                    d0 = d(i-1);
                end
                
%                 d(i) = CO2.d_fzero(t(i), d_ig(i), d0);
                d(i) = CO2.d_newton(t(i), d_ig(i), d0);
            end
            
            rho = reshape(d.*CO2.rhoc, sz);
        end
        
        function s = s_pT(p, T)
            % Entropy in kJ/(kg K) as function of pressure in MPa and temperature in K
            rho = CO2.rho_pT(p, T);
            s = CO2.s_rhoT_i(rho, T);
            s = reshape(s, size(rho));
        end
        
        function u = u_pT(p, T)
            % Internal energy in kJ/kg as function of pressure in MPa and temperature in K
            rho = CO2.rho_pT(p, T);
            u = CO2.u_rhoT_i(rho, T);
            u = reshape(u, size(rho));
        end
        
        function cv = cv_pT(p, T)
            % Isochoric heat capacity in kJ/(kg K) as function of pressure in MPa and temperature in K
            rho = CO2.rho_pT(p, T);
            cv = CO2.cv_rhoT_i(rho, T);
            cv = reshape(cv, size(rho));
        end        
 
        function h = h_pT(p, T)
            % Enthalpy in kJ/kg as function of pressure in MPa and temperature in K
            rho = CO2.rho_pT(p, T);
            h = CO2.h_rhoT_i(rho, T);
            h = reshape(h, size(rho));
        end

        function cp = cp_pT(p, T)
            % Isobaric heat capacity in kJ/(kg K) as function of pressure in MPa and temperature in K
            rho = CO2.rho_pT(p, T);
            cp = CO2.cp_rhoT_i(rho, T);
            cp = reshape(cp, size(rho));
        end        

        function cs = cs_pT(p, T)
            % Saturated liquid heat capacity in kJ/(kg K) as function of pressure in MPa and temperature in K
            rho = CO2.rho_pT(p, T);
            cs = CO2.cs_rhoT_i(rho, T);
            cs = reshape(cs, size(rho));
        end
        
        function w = w_pT(p, T)
            % Speed of sound in m/s as function of pressure in MPa and temperature in K
            rho = CO2.rho_pT(p, T);
            w = CO2.w_rhoT_i(rho, T);
            w = reshape(w, size(rho));
        end
        
        function mu = mu_pT(p, T)
            % Joule-Thomson coefficient in K/MPa as function of pressure in MPa and temperature in K
            rho = CO2.rho_pT(p, T);
            mu = CO2.mu_rhoT_i(rho, T);
            mu = reshape(mu, size(rho));
        end
        
        function f = f_pT(p, T)
            % Fugacity coefficient as function of pressure in MPa and temperature in K
            rho = CO2.rho_pT(p, T);
            f = CO2.f_rhoT_i(rho, T);
            f = reshape(f, size(rho));
        end
    end
    
    methods(Static, Hidden) % Properties in (rho, T) internal 
        % Thermodynamic properties as functions of density and temperature
        % according to Table 3.
        % No validity check is performed for rho and T.
        % The result will always be returned as a column vector, regardless
        % of the input dimensions.

        function p = p_rhoT_i(rho, T)
            % Pressure in MPa as function of density in kg/m3 and temperature in K
            d = rho(:) ./ CO2.rhoc;
            t = CO2.Tc ./ T(:);
            phir = CO2.phir_master(d,t,1,0);
            p = rho.*CO2.R.*T(:).*( 1 + d.*phir{2,1} ) .* 1e-3;
        end

        function s = s_rhoT_i(rho, T)
            % Entropy in kJ/(kg K) as function of density in kg/m3 and temperature in K
            d = rho(:) ./ CO2.rhoc;
            t = CO2.Tc ./ T(:);
            phi0 = CO2.phi0_master(d,t,0,1);
            phir = CO2.phir_master(d,t,0,1);
            s = CO2.R.*( t.*( phi0{1,2}+phir{1,2} ) - phi0{1,1} - phir{1,1} );
        end

        function u = u_rhoT_i(rho, T)
            % Internal energy in kJ/kg as function of density in kg/m3 and temperature in K
            d = rho(:) ./ CO2.rhoc;
            t = CO2.Tc ./ T(:);
            phi0 = CO2.phi0_master(d,t,0,1);
            phir = CO2.phir_master(d,t,0,1);
            u = CO2.R.*T(:).*( t.*( phi0{1,2}+phir{1,2} ) );
        end

        function cv = cv_rhoT_i(rho, T)
            % Isochoric heat capacity in kJ/(kg K) as function of density in kg/m3 and temperature in K
            d = rho(:) ./ CO2.rhoc;
            t = CO2.Tc ./ T(:);
            phi0 = CO2.phi0_master(d,t,0,2);
            phir = CO2.phir_master(d,t,0,2);
            cv = CO2.R.*( -t.^2.*( phi0{1,3}+phir{1,3} ) );
        end        

        function h = h_rhoT_i(rho, T)
            % Enthalpy in kJ/kg as function of density in kg/m3 and temperature in K
            d = rho(:) ./ CO2.rhoc;
            t = CO2.Tc ./ T(:);
            phi0 = CO2.phi0_master(d,t,0,1);
            phir = CO2.phir_master(d,t,1,1);
            h = CO2.R.*T(:).*( 1 + t.*( phi0{1,2}+phir{1,2}  ) + d.*phir{2,1} );
        end

        function cp = cp_rhoT_i(rho, T)
            % Isobaric heat capacity in kJ/(kg K) as function of density in kg/m3 and temperature in K
            d = rho(:) ./ CO2.rhoc;
            t = CO2.Tc ./ T(:);
            phi0 = CO2.phi0_master(d,t,0,2);
            phir = CO2.phir_master(d,t,2,2);
            cp = CO2.R.*( -t.^2.*( phi0{1,3}+phir{1,3} ) + ...
                ( 1 + d.*phir{2,1} - d.*t.*phir{2,2} ).^2 ./ ...
                ( 1 + 2.*d.*phir{2,1} + d.^2.*phir{3,1} ) );
        end        

        function cs = cs_rhoT_i(rho, T)
            % Saturated liquid heat capacity in kJ/(kg K) as function of density in kg/m3 and temperature in K
            d = rho(:) ./ CO2.rhoc;
            t = CO2.Tc ./ T(:);
            phi0 = CO2.phi0_master(d,t,0,2);
            phir = CO2.phir_master(d,t,2,2);
            cs = CO2.R.*( -t.^2.*( phi0{1,3}+phir{1,3} ) + ...
                ( 1 + d.*phir{2,1} - d.*t.*phir{2,2} ) ./ ...
                ( 1 + 2.*d.*phir{2,1} + d.^2.*phir{3,1} ) .* ...
                ( ( 1 + d.*phir{2,1} - d.*t.*phir{2,2} ) - ...
                CO2.rhoc/CO2.R./d .* CO2.dps_dT(T(:)) ));
        end

        function w = w_rhoT_i(rho, T)
            % Speed of sound in m/s as function of density in kg/m3 and temperature in K
            d = rho(:) ./ CO2.rhoc;
            t = CO2.Tc ./ T(:);
            phi0 = CO2.phi0_master(d,t,0,2);
            phir = CO2.phir_master(d,t,2,2);
            w = sqrt( CO2.R.*1e3.*T(:).*( ...
                1 + 2.*d.*phir{2,1} + d.^2.*phir{3,1} - ...
                ( 1 + d.*phir{2,1} - d.*t.*phir{2,2} ).^2 ./ ...
                ( t.^2.*( phi0{1,3}+phir{1,3} ) ) ...
                ));
        end

        function mu = mu_rhoT_i(rho, T)
            % Joule-Thomson coefficient in K/MPa as function of density in kg/m3 and temperature in K
            d = rho(:) ./ CO2.rhoc;
            t = CO2.Tc ./ T(:);
            phi0 = CO2.phi0_master(d,t,0,2);
            phir = CO2.phir_master(d,t,2,2);
            mu = 1e3./( CO2.R.*rho(:) ) .* ... 
                ( -d.*phir{2,1} - d.^2.*phir{3,1} - d.*t.*phir{2,2} ) ./ ...
                ( ( 1 + d.*phir{2,1} - d.*t.*phir{2,2} ).^2 - ...
                t.^2.*( phi0{1,3}+phir{1,3} ).*( 1 + 2.*d.*phir{2,1} + d.^2.*phir{3,1} ) );
        end

        function f = f_rhoT_i(rho, T)
            % Fugacity coefficient as function of density in kg/m3 and temperature in K
            d = rho(:) ./ CO2.rhoc;
            t = CO2.Tc ./ T(:);
            phir = CO2.phir_master(d,t,1,0);
            f = exp( phir{1,1} + d.*phir{2,1} - log( 1 + d.*phir{2,1} ) );
        end
    end

    methods(Static, Hidden) % Auxilary functions 
        function [a, b, sz] = checkInput(a, b)
            % Check if inputs are of same size or if one of them is scalar.
            % Return the row vectors of each input and the output size.
            
            if ~isscalar(a) 
                if ~isscalar(b)
                    if ~all(size(a) == size(b))
                        ME = MException('CO2:inputError','Input must be of same size.');
                        throwAsCaller(ME)
                    end
                end
                sz = size(a);
            else
                sz = size(b);
            end
            
            a = a(:);
            b = b(:);
        end
         
        function d = d_newton(t, d_ig, d0)
            % Find the reduced density that satisfies p(d,t) = p
            
            % The corresponding problem is finding the zero point of
            % f(d) = d + d^2.*dphir_dd - d_ig = 0
            % t ... reduced temperature
            % d_ig ... reduced density of the ideal gas
            % d0 ... starting point
            
            % Uses newtons' algorithm
            
            reltol = 1e-6;
            abstol = 1e-4;
            maxiter = 1e2;
            k = 0;
            
            
            d = d0; % best approximation to the zero point so far
            fd = zeros(size(d0)); % function value at d
            todo = true(size(d0)); % points that are not yet within tol
            err = zeros(size(d0)); % relative error
            
            while any(todo) && k <= maxiter                
                % Compute new function value
                d_ = d(todo);
                phir = CO2.phir_master(d_,t(todo),2,0);
                fd(todo) = d_ + d_.^2.*phir{2,1} - d_ig(todo);
                
                % Newton step f(x) / f'(x)
                delta = fd(todo) ./ ( 1+2.*d_.*phir{2,1} + d_.^2.*phir{3,1} );
                err(todo) = abs(delta)./d_;
                % Limit step to 0.5.*d to avoid negative densities
                idx = ( err(todo) > 0.5 );
                delta(idx)= 0.5.*d_(idx).*sign(delta(idx));

                d(todo) = d_ - delta;
                
                todo = abs(fd) > abstol | err > reltol;
                k = k+1;
            end                        
        end % function
        
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
            p_s = CO2.pVap(T(idx));
            T_ = 1 - T(idx)/CO2.Tc;
            
            val = nan(size(T));
            val(idx) = -p_s./T(idx).*( log(p_s/CO2.pc) + ...
                ( a1.*t1.*T_.^(t1-1) + a2.*t2.*T_.^(t2-1) + a3.*t3.*T_.^(t3-1) + a4.*t4.*T_.^(t4-1) ) ); 
        end
        
        function [valid, rho, T] = isSinglePhase_rhoT(rho, T)
            % Exclude (rho,T) couples within the two phase region. 
            [rho, T, sz] = CO2.checkInput(rho, T);
            
            idx = CO2.Tt < T & T < CO2.Tc;
            rhoLiquid = CO2.rhoLiqSat(T(idx));
            rhoVapor = CO2.rhoVapSat(T(idx));
            
            if isscalar(rho)
                % rho is scalar; T may be of any size.
                valid = true(size(T));
                valid(idx) = ~(rhoVapor < rho & rho < rhoLiquid); 
                
                rho = rho(any(valid));
                T = T(valid);
            elseif isscalar(T)
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
            
            % make size of valid match the input size 
            valid = reshape(valid, sz);
        end % function isSinglePhase_rhoT
    end % methods auxilary
    
    methods(Static, Hidden) % Helmholz function 
        function val = phi0_master(d, t, nd, nt)
            % Dimensionless helmholz function phi0 and its derivatives
            % according to Table 28 in Span et al. (1994)
            % phi0 ... dimensionless helmholz function
            % d ... reduced density delta = rho / rhoc
            % t ... inverse reduced temperature tau = Tc / T
            % nd ... number of derivatives w.r.t. d
            % nt ... number of derivatives w.r.t. t
            % Returns a cell array with (nd+1 x nt+1) cells. The (1,1) cell
            % contains the values of phi0, in column direction the
            % derivatives w.r.t. d are places; in row direction the
            % derivatives w.r.t. t. Consequently, the second derivative
            % w.r.t. t is stored in cell (1,3); the first mixed derivative
            % in (2,2).
            
            val = cell(nd+1, nt+1);
            
            % Coefficients for the calculation of phi0 from Table 27
            phi0_a = [8.37304456 -3.70454304 2.5 1.99427042 0.62105248 0.41195293 1.04028922 0.08327678];
            phi0_theta = [0 0 0 3.15163 6.11190 6.77708 11.32384 27.08792];

            
            phi0 = log(d) + phi0_a(1) + phi0_a(2).*t + phi0_a(3).*log(t) + ...
                sum( phi0_a(4:8) .* log( 1-exp(-phi0_theta(4:8).*t) ) ,2);
            val{1,1} = phi0;
            
            if nd > 0
                dphi0_dd = d.^-1;
                val{2,1} = dphi0_dd;
            end
            
            if nd > 1
                ddphi0_dddd = -d.^-2;
                val{3,1} = ddphi0_dddd;
            end
            
            if nt > 0
                e = exp(-phi0_theta(4:8).*t);
                dphi0_dt = phi0_a(2) + phi0_a(3)./t + ...
                    sum( phi0_a(4:8).*phi0_theta(4:8).*( (1-e).^-1 - 1 ) ,2);
                val{1,2} = dphi0_dt;
            end
            
            if nt > 1
                ddphi0_dtdt = -phi0_a(3)./t.^2 - ...
                    sum( phi0_a(4:8).*phi0_theta(4:8).^2.*e.*( 1-e ).^-2 ,2);
                val{1,3} = ddphi0_dtdt;
            end
            
            if nd > 0 && nt > 0
                ddphi0_dddt = zeros(size(d));
                val{2,2} = ddphi0_dddt;
            end
        end % function
        
        function val = phir_master(d, t, nd, nt)
            % Residual part of the dimensionless helmholz function phir and its derivatives
            % according to Table 32  in Span et al. (1994)
            % phir ... residual part of the dimensionless helmholz function
            % d ... reduced density delta = rho / rhoc
            % t ... inverse reduced temperature tau = Tc / T
            % nd ... number of derivatives w.r.t. d
            % nt ... number of derivatives w.r.t. t
            % Returns a cell array with (nd+1 x nt+1) cells. The (1,1) cell
            % contains the values of phir, in column direction the
            % derivatives w.r.t. d are places; in row direction the
            % derivatives w.r.t. t. Consequently, the second derivative
            % w.r.t. t is stored in cell (1,3); the first mixed derivative
            % in (2,2).
            %
            % Delta ... Distance function
            % Theta ... Some other function
            % Psi ... Exponential function
            
            val = cell(nd+1, nt+1);

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

            
            Theta = (1-t) + phir_A .* ( (d-1).^2 ).^( 1./(2.*phir_bb) );
            Delta = Theta.^2 + phir_B .* ( (d-1).^2 ).^phir_a;
            Psi = exp( -phir_C.*(d-1).^2 - phir_D.*(t-1).^2 );
        
            phir = sum( phir_n(1:7) .* d.^phir_d(1:7) .* t.^phir_t(1:7) ,2) + ...
                sum( phir_n(8:34) .* d.^phir_d(8:34) .* t.^phir_t(8:34) .* exp(-d.^phir_c) ,2) + ...
                sum( phir_n(35:39) .* d.^phir_d(35:39) .* t.^phir_t(35:39) .* ...
                    exp( -phir_alpha.*(d-phir_epsilon).^2 - phir_beta.*(t-phir_gamma).^2 ) ,2) + ...
                sum( phir_n(40:42) .* Delta.^phir_b .* d .* Psi, 2);
            val{1,1} = phir;
            
            if nd > 0
                dDelta_dd = (d-1) .* ( ...
                    phir_A .* Theta .* 2./phir_bb .* ( (d-1).^2 ).^( 1./(2.*phir_bb) - 1) + ...
                    2 .* phir_B .* phir_a .* ( (d-1).^2 ).^(phir_a-1) ...
                    );
                dDeltab_dd = phir_b .* Delta.^(phir_b-1) .* dDelta_dd;
                dDeltab_dd(isnan(dDeltab_dd)) = 0;
                dPsi_dd = -2 .* phir_C .* (d-1) .* Psi;

                
                dphir_dd = sum( phir_n(1:7) .* phir_d(1:7) .* d.^(phir_d(1:7)-1) .* t.^phir_t(1:7) ,2) + ...
                    sum( phir_n(8:34).*exp(-d.^phir_c).*( d.^(phir_d(8:34)-1) .* t.^phir_t(8:34) .* ( phir_d(8:34) - phir_c.*d.^phir_c ) ) ,2) + ...
                    sum( phir_n(35:39) .* d.^phir_d(35:39) .* t.^phir_t(35:39) .* ...
                        exp( -phir_alpha.*(d-phir_epsilon).^2 - phir_beta.*(t-phir_gamma).^2 ) .* ...
                        ( phir_d(35:39)./d - 2.*phir_alpha.*(d-phir_epsilon) ) ,2) + ...
                    sum( phir_n(40:42) .* ( Delta.^phir_b.*( Psi + d.*dPsi_dd ) + ...
                            dDeltab_dd .* d .* Psi ) ,2);
                val{2,1} = dphir_dd;
            end
            
            if nd > 1
                ddDelta_dddd = 1./(d-1) .* dDelta_dd + (d-1).^2 .* ( ...
                    4 .* phir_B .* phir_a .* ( phir_a-1 ) .* ( (d-1).^2 ).^( phir_a-2 ) + ...
                    2 .* phir_A.^2 .* phir_bb.^-2 .* (( (d-1).^2 ).^( 1./(2.*phir_bb) - 1 )).^2 + ...
                    phir_A .* Theta .* 4 ./ phir_bb .* (1./(2.*phir_bb) -1) .* ( (d-1).^2 ).^( 1./(2.*phir_bb) - 2 ) ...
                    );
                ddDeltab_dddd = phir_b .* ( ...
                    Delta.^(phir_b-1) .* ddDelta_dddd + ...
                    (phir_b-1) .* Delta.^(phir_b-2) .* dDelta_dd.^2 ...
                    );
                ddDeltab_dddd(isnan(ddDeltab_dddd)) = 0;
                ddPsi_dddd = ( 2.*phir_C.*(d-1).^2 - 1).*2.*phir_C.*Psi;
                
                ddphir_dddd = sum( phir_n(1:7) .* phir_d(1:7) .* (phir_d(1:7)-1) .* d.^(phir_d(1:7)-2) .* t.^phir_t(1:7) ,2) + ...
                    sum( phir_n(8:34).*exp(-d.^phir_c).*( d.^(phir_d(8:34)-2) .* t.^phir_t(8:34) .* ( ...
                        ( phir_d(8:34) - phir_c.*d.^phir_c ) .* ( phir_d(8:34) - 1 - phir_c.*d.^phir_c ) - ...
                        phir_c.^2.*d.^phir_c )) ,2) + ...
                    sum( phir_n(35:39) .* t.^phir_t(35:39) .* exp( -phir_alpha.*(d-phir_epsilon).^2 - phir_beta.*(t-phir_gamma).^2 ) .* ...
                        ( -2.*phir_alpha.*d.^phir_d(35:39) + 4.*phir_alpha.^2.*d.^phir_d(35:39).*( d-phir_epsilon ).^2 - ...
                        4.*phir_d(35:39).*phir_alpha.*d.^( phir_d(35:39)-1 ).*( d-phir_epsilon ) + phir_d(35:39).*( phir_d(35:39)-1 ).*d.^( phir_d(35:39)-2 ) ) ,2) + ...
                    sum( phir_n(40:42) .* ( Delta.^phir_b.*( 2.*dPsi_dd + d.*ddPsi_dddd ) + ...
                        2.*dDeltab_dd.*( Psi + d.*dPsi_dd ) + ddDeltab_dddd.*d.*Psi ) ,2);
                val{3,1} = ddphir_dddd;
            end
            
            if nt > 0
                dDeltab_dt = -2 .* Theta .* phir_b .* Delta.^( phir_b-1 );
                dDeltab_dt(isnan(dDeltab_dt)) = 0;

                dPsi_dt = -2 .* phir_D .* (t-1) .* Psi;

                dphir_dt = sum( phir_n(1:7) .* phir_t(1:7) .* d.^phir_d(1:7) .* t.^( phir_t(1:7)-1 ) ,2) + ...
                    sum( phir_n(8:34) .* phir_t(8:34) .* d.^phir_d(8:34) .* t.^( phir_t(8:34)-1 ) .* exp( -d.^phir_c ) ,2) + ...
                    sum( phir_n(35:39) .* d.^phir_d(35:39) .* t.^phir_t(35:39) .* ...
                        exp( -phir_alpha.*(d-phir_epsilon).^2 - phir_beta.*(t-phir_gamma).^2 ) .* ...
                        ( phir_t(35:39)./t - 2.*phir_beta.*(t-phir_gamma) ) ,2) + ...
                    sum( phir_n(40:42).*d.*( dDeltab_dt.*Psi + Delta.^phir_b.*dPsi_dt ) ,2);
                val{1,2} = dphir_dt;
            end
            
            if nt > 1
                ddDeltab_dtdt = 2 .* phir_b .* Delta.^( phir_b-1 ) + ...
                    4 .* Theta.^2 .* phir_b .* (phir_b-1) .* Delta.^(phir_b-2);
                ddDeltab_dtdt(isnan(ddDeltab_dtdt)) = 0;
                ddPsi_dtdt = ( 2.*phir_D.*(t-1).^2 - 1).*2.*phir_D.*Psi;

                ddphir_dtdt = sum( phir_n(1:7) .* phir_t(1:7) .* ( phir_t(1:7)-1 ) .* d.^phir_d(1:7) .* t.^( phir_t(1:7)-2 ) ,2) + ...
                    sum( phir_n(8:34) .* phir_t(8:34) .* ( phir_t(8:34)-1 ) .* d.^phir_d(8:34) .* t.^( phir_t(8:34)-2 ) .* exp(-d.^phir_c) ,2) + ...
                    sum( phir_n(35:39) .* d.^phir_d(35:39) .* t.^phir_t(35:39) .* ...
                        exp( -phir_alpha.*(d-phir_epsilon).^2 - phir_beta.*(t-phir_gamma).^2 ) .* ...
                        ( ( phir_t(35:39)./t - 2.*phir_beta.*(t-phir_gamma) ).^2 - phir_t(35:39)./t.^2 - 2.*phir_beta ) ,2)  + ...
                    sum( phir_n(40:42).*d.*( ddDeltab_dtdt.*Psi + 2.*dDeltab_dt.*dPsi_dt + Delta.^phir_b.*ddPsi_dtdt ) ,2);
                val{1,3} = ddphir_dtdt;
            end
            
            if nd > 0 && nt > 0
                ddDeltab_dddt = -phir_A .* phir_b .* 2./phir_bb .* Delta.^(phir_b-1) .* (d-1) .* ( (d-1).^2 ).^( 1./(2.*phir_bb) - 1 ) - ...
                    2 .* Theta .* phir_b .* (phir_b-1) .* Delta.^(phir_b-2) .* dDelta_dd;
                ddDeltab_dddt(isnan(ddDeltab_dddt)) = 0;
                ddPsi_dddt = 4.*phir_C.*phir_D.*(d-1).*(t-1).*Psi;
                
                ddphir_dddt = sum( phir_n(1:7) .* phir_d(1:7) .* phir_t(1:7) .* d.^( phir_d(1:7)-1 ) .* t.^( phir_t(1:7)-1 ) ,2) + ...
                    sum( phir_n(8:34).*exp(-d.^phir_c).*d.^( phir_d(8:34)-1 ).*phir_t(8:34).* t.^( phir_t(8:34)-1 ).*( phir_d(8:34) - phir_c.*d.^phir_c ) ,2) + ...
                    sum( phir_n(35:39) .* d.^phir_d(35:39) .* t.^phir_t(35:39) .* ...
                        exp( -phir_alpha.*(d-phir_epsilon).^2 - phir_beta.*(t-phir_gamma).^2 ) .* ...
                        ( phir_d(35:39)./d - 2.*phir_alpha.*(d-phir_epsilon) ) .* ...
                        ( phir_t(35:39)./t - 2.*phir_beta.*(t-phir_gamma) ) ,2) + ...
                    sum( phir_n(40:42) .* (  Delta.^phir_b.*( dPsi_dt + d.*ddPsi_dddt ) + ...
                        d.*dDeltab_dd.*dPsi_dt + dDeltab_dt.*( Psi + d.*dPsi_dd ) + ...
                        ddDeltab_dddt.*d.*Psi ) ,2);
                val{2,2} = ddphir_dddt;
            end

        end % function
    end % methods

    methods(Static) % Transport properties 
        function [eta, lambda] = transport_rhoT(rho, T)
            % Transport properties as functions of density in kg/m and
            % temperature in K. The model can be selected in
            % CO2.transport_rhoT_i(rho, T).
            % A validity check is performed. For points in the two-phase region
            % nan is returned. If the validity check is not required the
            % internal functions with the suffix rhoT_i can be used instead.
            % The solid region is not checked.
            % Viscosity eta in mPa s
            % Thermal conductivity lambda in mW/m K
            
            [v, rho, T] = CO2.isSinglePhase_rhoT(rho, T);
            eta = nan(size(v));
            lambda = nan(size(v));
            [eta(v), lambda(v)] = CO2.transport_rhoT_i(rho, T);
        end
        
        function lambda=lambda_in(rho,T)
            [~,lambda]=CO2.transport_rhoT(rho,T);
        end
        
        function [eta, lambda] = transport_pT(p, T)
            % eta in mPa s
            % lambda in mW/m K
            rho = CO2.rho_pT(p, T);
            [eta, lambda] = CO2.transport_rhoT_i(rho, T);
            eta = reshape(eta, size(rho));
            lambda = reshape(lambda, size(rho));
        end
        
        function [eta, lambda] = transport_rhoT_i(rho, T)
            % Transport properties as functions of density in kg/m and
            % temperature in K. Two models are available:
            % 1) The modern one based on Laesecke et al. (2017),
            %   Luettmer-Strathmann et al. (1995) and Huber et al. (2016),
            %   which should be equivalent to REFPROP version 10.
            % 2) The depreciated one based on Vesovic (1990) and Fenghour
            %   (1997), which should be equivalent to REFPROP version 7.
            % The model can be selected by commenting the corresponding
            % lines of code.
            % Viscosity eta in mPa s
            % Thermal conductivity lambda in mW/(m K)
            
            T = T(:); rho = rho(:);
            
            % viscosity from Laesecke (2017) with critical enhancement from
            % Luettmer-Strathmann (1995)
            % thermal conductivity and critical enhancement from Huber
            % (2016). This should be equivalent to REFPROP version 10
            eta = CO2.eta_rhoT_Laesecke(rho, T); % mPa
            lambda = CO2.lambda_rhoT_Huber(rho, T); % mW/m K
            [eta, lambda] = CO2.criticalEnhancement_LS_Huber(eta, lambda, rho, T);
            % The critical enhancement of the viscosity affects only a very
            % narrow region around the critical point and also it is not
            % very pronounced. Often it can be neglected. At the same time,
            % it requires a comparably large computational effort. If
            % performance is an issue, and the critical enhancement of
            % viscosity can be neglected the following line can be used
            % instead of the one above. See Olehowy and Sengers (1989): A 
            % Simplified Representation for the Thermal Conductivity of 
            % Fluids in the Critical Region. International Journal of 
            % Thermophysics, Vol. 10, No. 2. https://doi.org/10.1007/BF01133538
%             lambda = CO2.criticalEnhancement_Huber(eta, lambda, rho, T);
 
            % viscosiy from Vesovic (1990) and Fenghour (1997)
            % thermal conductivity from Vesovic (1990)
            % accomodation function for the vicintiy of the critical point
            % according to Vesovic (1990) (with crossover function for the
            % viscosity and the simplified approach for the thermal
            % conductivity)
            % this should be equivalent to REFPROP version 7
%             lambda = CO2.lambda_rhoT_Vesovic(rho, T); % mW/(m K)
%             eta = CO2.eta_rhoT_Vesovic(rho, T).*1e-3; % mPa
%             [eta, lambda] = CO2.criticalEnhancement_Vesovic(eta, lambda, rho, T);
       end
    end % methods Transport properties
    
    methods(Static, Hidden) % Transport auxiliary 
        function [eta, lambda] = criticalEnhancement_Vesovic(eta, lambda, rho, T)
            % according to Vesovic et. al (1990)
            % eta in muPa s
            % lambda in mW/(m K)
            % rho in kg/m
            % T in K
            
            % constants according to Table 5
            k = 1.38064852e-23; % J/K
            
            % prepare input for EOS to compute isothermal compressibility
            d = rho./CO2.rhoc;
            t = CO2.Tc./T;
            x = max(numel(d), numel(t));
            % correction according to Eq. (49)
            idxCorr = T > 445;
            t(idxCorr) = CO2.Tc./445;
            % add values for reference temperature
            if isscalar(t)
                t = repmat(t, size(d));
            end
            t = [t; repmat(CO2.Tc/450, size(d))];
            if ~isscalar(d)
                d = [d; d];
            end
            % isothermal compressibility from EOS
            phir = CO2.phir_master(d,t,2,0);
            drho_dp = 1./( CO2.Tc./t.*CO2.R.*1e-3.*( 1+ 2.*d.*phir{2,1} + d.^2.*phir{3,1}) );
            
            % reduced symetrized compressibility difference according to
            % Eq. (47) and (40)
            dX = CO2.pc/CO2.rhoc^2/CO2.Tc .* rho .* ...
                ( T.*drho_dp(1:x) - 450.^2./T.*drho_dp(x+1:end) );
            % xi is bounded by 0; for dX <= 0, xi = 0
            % A Note on the Critical Enhancement of Transport Properties 
            % and Correlation Length of Fluids DOI 10.1007/s10765-013-1519-7
            % For xi = 0, critical enhancement of the transport properties
            % is zero. Consequently values of dX <= 0 can be omitted.
            idxC = dX > 0;
            
            % correlation length according to Eq. (46)
            xi = zeros(size(dX));
            xi(idxC) = 1.5e-10.*(dX(idxC)/0.052).^(0.630/1.2415);
            % correction according to Eq. (49)
            xi(idxCorr) = xi(idxCorr).*exp(-(T(idxCorr) - 445)/10);

            
            % equilibrium properties from EOS
            if isscalar(rho)
                rhoC = rho(any(idxC));
            else
                rhoC = rho(idxC);
            end
            if isscalar(T)
                TC = T(any(idxC));
            else
                TC = T(idxC);
            end
            t = CO2.Tc./TC;
            d = rhoC./CO2.rhoc;
            phi0 = CO2.phi0_master(d,t,0,2);
            phir = CO2.phir_master(d,t,2,2);
            cp = CO2.R.*( -t.^2.*( phi0{1,3}+phir{1,3} ) + ...
                ( 1 + d.*phir{2,1} - d.*t.*phir{2,2} ).^2 ./ ...
                ( 1 + 2.*d.*phir{2,1} + d.^2.*phir{3,1} ) );
            cv = CO2.R.*( -t.^2.*( phi0{1,3}+phir{1,3} ) );

            % parameters according to Table 3 and 4
            qD = (2.3e-10)^-1; % 1/m
            yD = atan(qD .* xi(idxC));
            foo = (1+(qD.*xi(idxC)).^2).^(1/2);
            yr = ( atan(qD.*xi(idxC)./foo) -yD)./foo;
            ya = rhoC.*k.*TC./( 8.*pi.*eta(idxC).^2.*xi(idxC) ) .*1e6;
            yb = lambda(idxC)./( eta(idxC).*( cp - cv ) ) .*1e-3;
            yg = cv./( cp - cv ); 
            
            H = CO2.eta_crossover_Vesovic(yD, yr, ya, yb, yg);

            % I did not manage to implement the crossover function for
            % the thermal conductivity. For this reason, only the
            % simplified model for critical enhancement according to
            % section 4.5 in Vesovic 1990 is included. At the NIST webpage,
            % they also only use the simplified model, so I guess its fine.
            qD_ = (4e-10)^-1;
            O = 2/pi.*( (cp-cv)./cp.*atan(qD_.*xi(idxC)) + cv./cp.*qD_.*xi(idxC) );
            O0 = 2/pi.*( 1-exp( -1./( (qD_.*xi(idxC)).^-1 + 1/3.*( qD_.*xi(idxC).*CO2.rhoc./rhoC ).^2 ) ) );
            lambda_c = rhoC.*cp.*1.01.*k.*TC ./ (6.*pi.*eta(idxC).*xi(idxC)) .* (O - O0) .*1e9;
            
            eta(idxC) = eta(idxC).*exp(0.06.*H);
            lambda(idxC) = lambda(idxC) + lambda_c;
        end % function criticalEnhancement_Vesovic
        
        function H = eta_crossover_Vesovic(yD, yr, ya, yb, yg)
            % viscosity crossover function for the critical region
            % according to Vesovic 1990 Table 4.
            ye = (yr + yb./ya)./yD;
            yn = yg.*yr./yD;
            
            % according to part IV in Table 4
            nu = zeros(numel(yD),3);
            for i = 1:numel(yD)
                nu(i,:) = -roots([1 ye(i) yg(i) yn(i)]);
            end

            % according to part III in Table 4
            M = cat(3, ones(numel(yD),3), -nu+sum(nu, 2), prod(nu,2)./nu);

            % according to part II in Table 4
            p = [...
                ye.^4 - 4.*ye.^2.*yg - 2.*ye.^2 + 2.*ye.*yn + 3.*yg.^2 + 4.*yg + 1, ...
                (ye.*yg-yn).*(ye.^2-3.*yg-2) + yg.^2.*(1+yg).^2./(ye.*yg-yn),...
                ye.*yn.*(ye.^2-3.*yg-2) + yn.^2 + yn.*yg.*(1+yg).^2./(ye.*yg-yn)...
                ];
            
            % solve equation systems M.*N = p
            N = zeros(size(nu));
            M = permute(M, [3 2 1]);
            p = permute(p, [2 3 1]);
            for i = 1:numel(yD)
                N(i,:) = M(:,:,i) \ p(:,:,i);
            end
            
            foo = (1-nu.^2).^(1/2);
            F = 1./foo .*log( (1 + nu + foo.*tan(yD/2))./(1 + nu - foo.*tan(yD/2)) ) ;
            
            H = (3.*yg.*ye + 3.*ye/2 - ye.^3 - yn).*yD ...
                + (ye.^2 - 2.*yg - 5/4).*sin(yD) - ye.*sin(2.*yD)/4 + sin(3.*yD)/12 ...
                + (yg.*(1+yg)).^(3/2)./(yn-yg.*ye).*atan((yg./(1+yg)).^(1/2).*tan(yD)) ...
                + sum(N.*F,2);
            H = real(H); % some complex part might remain due to numeric errors.
        end % function eta_crossover_Vesovic
        
        function eta = eta_rhoT_Vesovic(rho, T)
            % Viscosity of CO2 in muPa s as a function of density in kg/m
            % and temperature in K.
            % Valid from 200 K to 1500 K and densities up to 1400 kg/m. In
            % terms of pressure it is valid up to 300 MPa for temperatures
            % below 1000 K. Above 1000 K the limit is 30 MPa.
            % based on Vesovic et al. (1990) "The transport properties of
            % carbon dioxide" in J. Phys. Chem. Ref. Data, Vol. 19, No. 3
            % and Fenghour et al. (1998) "The viscosity of Carbon Dioxide"
            % in J. Phys. Chem. Ref. Data, Vol. 27, No. 1.
            
            
            % reduced temperature according to Eq. (5) and (6)
            T_ = T./251.196;
            
            % constants from Table 1
            a = [0.235156 -0.491266 5.211155e-2 5.347906e-2 -1.537102e-2];            
            % reduced effective crossection according to Eq. (4)
            G = exp(sum( a.*log(T_).^(0:4) ,2));
            % zero density eta according to Eq. (3)
            eta_0 = 1.00697.*T.^(1/2)./G;
            
            % constants from Table 3 in Fenghour 1998
            eta_d = [0.4071119e-2 0.7198037e-4 0.2411697e-16 0.2971072e-22 -0.1627888e-22];
            % excess eta according to Eq. (8) in Fenghour 1998
            eta_ex = sum(eta_d .* [rho rho.^2 rho.^6./T_.^3 rho.^8 rho.^8./T_] ,2);
            
%             % constants from Tabel II in Luettmer-Strathmann 1995
%             eta_c = [5.1512e-3 6.2634e-5 -3.9414e-9 3.1447e-11];
%             % excess eta according to Eq. (12) in Luettmer-Strathmann 1995
%             eta_ex = sum(eta_c.*rho.^(1:4), 2);
            
            eta = eta_0 + eta_ex;
        end % function eta_rhoT_Vesovic
        
        function lambda = lambda_rhoT_Vesovic(rho, T)
            % Thermal conductivity of CO2 in mW/m K as a function of
            % density in kg/m and temperature in K.
            % Valid from 200 K to 1000 K and pressures up to 100 MPa
            % based on Vesovic et al. (1990) "The transport properties of
            % carbon dioxide" in J. Phys. Chem. Ref. Data, Vol. 19, No. 3.
            % reduced temperature according to Eq. (5) and (6)
            T_ = T./251.196;
            
            % constants from Table 1
            b = [0.4226159 0.6280115 -0.5387661 0.6735941 0 0 -0.4362677 0.2255388];
            c = [2.387869e-2 4.350794 -10.33404 7.981590 -1.940558];
            
            % c_int/k according to Eq. (31)
            cintk = 1 + exp(-183.5./T) .* sum( c.*(T/100).^(1:-1:-3), 2);
            % reduced effective crossection according to Eq. (30)
            G = sum(b./T_.^(0:7), 2);
            % zero density lambda according to Eq. (29) and (12)
            lambda_0 = 0.475598 .* T.^(1/2) .* ( 1 + 2/5.*cintk ) ./ G;
            
            % constants from Table 8
            lambda_d = [2.447164e-2 8.705605e-5 -6.547950e-8 6.594919e-11];
            % excess lambda according to Eq. (63)
            lambda_ex = sum(lambda_d .* rho.^(1:4),2); 
            
            lambda = lambda_0 + lambda_ex;
        end % function lambda_rhoT_Vesovic
        
        function eta = eta_rhoT_Laesecke(rho, T)
            % Viscosity of CO2 in mPa s as a function of density in kg/m
            % and temperature in K.
            % Valid from 100 to 2000 K for gaseous CO2 and from 220 to 700
            % K with pressures along the melting line up to 8000 MPa for
            % compressed and supercritical liquid states.
            % based on Laesecke and Muzny (2017) "Reference Correlation for
            % the Viscosity of Carbon Dioxide" in J. Phys. Chem. Ref. Data
            % Vol. 46, No. 1
            T = T(:); rho = rho(:);
            
            % constants from Table 2
            a = [1749.354893188350 -369.069300007128 5423856.34887691 -2.21283852168356 -269503.247933569 73145.021531826 5.34368649509278];
            % zero density eta according to Eq. (4)
            eta_0 = 1.0055.*sqrt(T) ./ ( ...
                a(1) + a(2).*T.^(1/6) + a(3).*exp(a(4).*T.^(1/3)) ...
                + ( a(5)+a(6).*T.^(1/3) )./( exp(T.^(1/3)) ) + a(7).*sqrt(T) ...
                );
            
            %constants from Table 3
            b0 = -19.572881;
            b = [219.73999 -1015.3226 2471.0125 -3375.1717 2491.6597 -787.26086 14.085455 -0.34664158];
            t = [0.25 0.5 0.75 1 1.25 1.5 2.5 5.5];
            % other constants
            e_kb = 200.760; % K - energy scaling parameter            
            Na = 6.02214085774e23; % 1/mol - Avogadro constant
            M = 44.0095e-3; % kg/mol - molar mass of CO2
            s = 0.378421; % nm - length scaling parameter
            % linear-in-density viscosity coefficient according to Eq. (5)
            % and (6)
            B = b0 + sum(b./(T./e_kb).^t, 2);
            eta_l = eta_0.*B.*s^3.*Na/M.*1e-27;
            
            % constants from Table 8
            g = 8.06282737481277;
            c1 = 0.360603235428487;
            c2 = 0.121550806591497;
            % dimensioning factor according to Eq. (9)
            eta_tl = 1e3.*CO2.rhotLiquid^(2/3).*sqrt(CO2.RM.*CO2.Tt)/M^(1/6)/Na^(1/3);
            % residual viscosity according to Eq. (8)
            Tr = T/CO2.Tt;
            rhor = rho/CO2.rhotLiquid;
            eta_r = eta_tl.*( ...
                c1.*Tr.*rhor.^3 ...
                + (rhor.^2 + rhor.^g)./(Tr - c2)...
                );
            
            % background viscosity according to Eq. (10)
            eta = eta_0 + rho.*eta_l + eta_r;
        end % function eta_rhoT_Laesecke
        
        function [eta, lambda] = criticalEnhancement_LS(eta, lambda, rho, T)
            % according to Luettmer-Strathmann et al. (1995)
            % eta in mPa s
            % lambda in mW/m K
            % rho in kg/m
            % T in K
            
            % constants
            k = 1.3806485279e-23; % J/K - Bolzmann constant
            Tref = 3/2.*CO2.Tc; % K
            
            % prepare input for EOS to compute isothermal compressibility
            d = rho./CO2.rhoc;
            t = CO2.Tc./T;
            x = max(numel(d), numel(t));
            % add values for reference temperature
            if isscalar(t)
                t = repmat(t, size(d));
            end
            t = [t; repmat(CO2.Tc/Tref, size(d))];
            if ~isscalar(d)
                d = [d; d];
            end
            % isothermal compressibility from EOS
            phir = CO2.phir_master(d,t,2,0);
            drho_dp = 1./( CO2.Tc./t.*CO2.R.*1e-3.*( 1+ 2.*d.*phir{2,1} + d.^2.*phir{3,1}) );
            
            % reduced symetrized compressibility difference according to
            % Eq. (6.10)
            dX = CO2.pc/CO2.rhoc^2 .* rho .* ...
                ( drho_dp(1:x) - Tref./T.*drho_dp(x+1:end) );
            % xi is bounded by 0; for dX <= 0, xi = 0
            % A Note on the Critical Enhancement of Transport Properties 
            % and Correlation Length of Fluids DOI 10.1007/s10765-013-1519-7
            % For xi = 0, critical enhancement of the transport properties
            % is zero. Consequently values of dX <= 0 can be omitted.
            idxC = dX > 0;
            
            % correlation length according to Eq. (6.10)
            % constants from Table II
            xi = zeros(size(dX));
            xi(idxC) = 1.5e-10.*(dX(idxC)/0.0481).^(0.630/1.239);

            % equilibrium properties from EOS
            if isscalar(rho)
                rhoC = rho(any(idxC));
            else
                rhoC = rho(idxC);
            end
            if isscalar(T)
                TC = T(any(idxC));
            else
                TC = T(idxC);
            end
            t = CO2.Tc./TC;
            d = rhoC./CO2.rhoc;
            phi0 = CO2.phi0_master(d,t,0,2);
            phir = CO2.phir_master(d,t,2,2);
            cp = CO2.R.*( -t.^2.*( phi0{1,3}+phir{1,3} ) + ...
                ( 1 + d.*phir{2,1} - d.*t.*phir{2,2} ).^2 ./ ...
                ( 1 + 2.*d.*phir{2,1} + d.^2.*phir{3,1} ) );
            cv = CO2.R.*( -t.^2.*( phi0{1,3}+phir{1,3} ) );

            % parameters according to Table I
            qD = (2.3e-10)^-1; % 1/m
            yD = atan(qD .* xi(idxC));
            foo = (1+(qD.*xi(idxC)).^2).^(1/2);
            yd = ( atan(qD.*xi(idxC)./foo) - yD )./foo;
            ya = rhoC.*k.*TC./( 8.*pi.*eta(idxC).^2.*xi(idxC) ) .*1e6;
            yb = lambda(idxC)./( eta(idxC).*( cp - cv ) ) .* 1e-3;
            yg = cv./( cp - cv ); 
            % adjustments for asymptotic behavior
            eta_ = 2-1.239/0.630;
            alpha = 2-3.*0.630;
            yb = yb.^(1/(1-eta_/2));
            yg = yg.^(1/(1-eta_/2-alpha/2/0.630));

            H = CO2.eta_crossover_LS(yD, yd, ya, yb, yg);
            eta(idxC) = eta(idxC).*exp(0.063.*H);

            if nargout > 1
                O = CO2.lambda_crossover_LS(yD, yd, ya, yb, yg);
                O0 = 2/pi.*( 1-exp( -1./( 1./(qD.*xi) + (qD.*xi).^3.*(CO2.rhoc/3./rho).^2 ))) ...
                    ./ (1 + ya.*(yD+yd) + yb.*(1+yg).^-1);
                lambda_c = rhoC.*cp.*1.05.*k.*TC ./ (6.*pi.*eta(idxC).*xi(idxC)) .* (O - O0) .*1e9;
                lambda(idxC) = lambda(idxC) + lambda_c;
            end
        end % function criticalEnhancement_LS
        
        function H = eta_crossover_LS(yD, yd, ya, yb, yg)
            % viscosity crossover function for the critical region
            % according to Luettmer-Strathmann (1995) Table I
            % this crossover function is also recommended by Laesecke (2017)
            ye = (yd + yb./ya)./yD;
            yn = yg.*yd./yD;

            c0 = ye.*yn.*( ye.^2-3.*yg-2 ) + yn.^2 - yg.*yn.*(1+yg).^2./( yn-yg.*ye );
            c1 = -1.*( yn-yg.*ye ).*( ye.^2-3.*yg-2 ) - ( yg.^2.*(1+yg).^2./(yn-yg.*ye) );
            c2 = ye.^4 - 4.*yg.*ye.^2 - 2.*ye.^2 + 2.*ye.*yn + 3.*yg.^2 + 4.*yg + 1;
            
            % roots of polynomial
            z = zeros(numel(yD),3);
            for i = 1:numel(yD)
                z(i,:) = roots([1 ye(i) yg(i) yn(i)]);
            end
            
            % F-function
            foo = (1-z.^2).^(1/2);
            F = 1./foo .*log( (1 - z + foo.*tan(yD/2))./(1 - z - foo.*tan(yD/2)) );
            
            % denominator of the term in the sum in H-function
            M = [1 -1 0; 1 0 -1];
            zz = zeros(size(z));
            for i = 1:3
                zz(:,i) = prod( z * M', 2);
                M = circshift(M, 1, 2);
            end
            
            % value of crossover function
            H = (3.*yg.*ye + 3.*ye/2 - ye.^3 - yn).*yD ...
                + (ye.^2 - 2.*yg - 5/4).*sin(yD) - ye.*sin(2.*yD)/4 + sin(3.*yD)/12 ...
                + (yg.*(1+yg)).^(3/2) ./ (yn-yg.*ye) .* atan((yg./(1+yg)).^(1/2).*tan(yD)) ...
                + sum((c2.*z.^2+c1.*z+c0)./zz.*F,2);
            H = real(H); % some complex part might remain due to numeric errors.
        end % function eta_crossover_LS
        
        function O = lambda_crossover_LS(yD, yd, ya, yb, yg)
            % viscosity crossover function for the critical region
            % according to Luettmer-Strathmann (1995) Table I
            
            a0 = yg.^2 - ya.*yg.*yd;
            a1 = -ya.*yg.*yD;
            a2 = yg - yb - ya.*yd;
            a3 = -ya.*yD;
            
            b0 = ya.*yg.*yd;
            b1 = ya.*yg.*yD;
            b2 = yg + yb + ya.*yd;
            b3 = ya.*yD;
            
            % roots of polynomial
            z = zeros(numel(yD),4);
            for i = 1:numel(yD)
                z(i,:) = roots([1 b3(i) b2(i) b1(i) b0(i)]);
            end
            
            % F-function
            foo = (1-z.^2).^(1/2);
            F = 1./foo .*log( (1 - z + foo.*tan(yD/2))./(1 - z - foo.*tan(yD/2)) ) ;
            
            % denominator of the term in the sum in O-function
            M = [1 -1 0 0; 1 0 -1 0; 1 0 0 -1];
            zz = zeros(size(z));
            for i = 1:4
                zz(:,i) = prod( z .* M', 2);
                M = circshift(M, 1, 2);
            end
            
            % value of crossover function
            O = 2/pi./(1+yg) .* ( yD + sum( (a3.*z.^3 + a2.*z.^2 + a1.*z + a0)./zz .*F ,2) );
            O = real(O); % some complex part might remain due to numeric errors.
        end % function lambda_crossover_LS
        
        function lambda = lambda_rhoT_Huber(rho, T)
            % Reference Correlation of the Thermal Conductivity of Carbon Dioxide from the Triple Point to 1100 K and up to 200 MPa
            Tr = reshape(T./ CO2.Tc,1,numel(T));
            rho=reshape(rho,1,numel(rho));
            
            % thermal conductivity at zero density according to Eq. (3)
            % constants from Table 3
            L = [1.51874307e-2;2.80674040e-2;2.28564190e-2;-7.41624210e-3];
            lambda_0 = sqrt(Tr)./sum( L ./ Tr.^((0:3)') ); % mW/m K
            
            % residual thermal conductivity according to Eq. (4)
            % constants from Table (4)
            B1 = [1.00128e-2;5.60488e-2;-8.11620e-2;6.24337e-2;-2.06336e-2;2.53248e-3];
            B2 = [4.30829e-3;-3.58563e-2;6.71480e-2;-5.22855e-2;1.74571e-2;-1.96414e-3];
            lambda_r = sum(( B1+B2.*Tr ).*( rho./CO2.rhoc ).^((1:6)') ,1).*1e3; % mW/m K
            
            lambda = (lambda_0 + lambda_r)';
        end % function lambda_rhoT_Huber
        
        function lambda = criticalEnhancement_Huber(eta, lambda, rho, T)
            % eta in mPa s
            % lambda in mW/m K
            
            % constants
            k = 1.3806485279e-23; % J/K - Bolzmann constant
            Tref = 3/2.*CO2.Tc; % K
            
            % prepare input for EOS to compute isothermal compressibility
            d = rho./CO2.rhoc;
            t = CO2.Tc./T;
            x = max(numel(d), numel(t));
            % add values for reference temperature
            if isscalar(t)
                t = repmat(t, size(d));
            end
            t = [t; repmat(CO2.Tc/Tref, size(d))];
            if ~isscalar(d)
                d = [d; d];
            end
            % isothermal compressibility from EOS
            phir = CO2.phir_master(d,t,2,0);
            drho_dp = 1./( CO2.Tc./t.*CO2.R.*1e-3.*( 1+ 2.*d.*phir{2,1} + d.^2.*phir{3,1}) );
            
            % reduced symetrized compressibility difference according to
            % Eq. (8)
            dX = CO2.pc/CO2.rhoc^2 .* rho .* ...
                ( drho_dp(1:x) - Tref./T.*drho_dp(x+1:end) );
            % xi is bounded by 0; for dX <= 0, xi = 0
            % A Note on the Critical Enhancement of Transport Properties 
            % and Correlation Length of Fluids DOI 10.1007/s10765-013-1519-7
            % For xi = 0, critical enhancement of the transport properties
            % is zero. Consequently values of dX <= 0 can be omitted.
            idxC = dX > 0;
            
            % correlation length according to Eq. (8)
            % constants from the text below Eq. (8)
            xi = zeros(size(dX));
            xi(idxC) = 1.5e-10.*(dX(idxC)/0.052).^(0.630/1.239);

            % equilibrium properties from EOS
            if isscalar(rho)
                rhoC = rho(any(idxC));
            else
                rhoC = rho(idxC);
            end
            if isscalar(T)
                TC = T(any(idxC));
            else
                TC = T(idxC);
            end
            t = CO2.Tc./TC;
            d = rhoC./CO2.rhoc;
            phi0 = CO2.phi0_master(d,t,0,2);
            phir = CO2.phir_master(d,t,2,2);
            cp = CO2.R.*( -t.^2.*( phi0{1,3}+phir{1,3} ) + ...
                ( 1 + d.*phir{2,1} - d.*t.*phir{2,2} ).^2 ./ ...
                ( 1 + 2.*d.*phir{2,1} + d.^2.*phir{3,1} ) );
            cv = CO2.R.*( -t.^2.*( phi0{1,3}+phir{1,3} ) );

            % critical enhancement according to Eq. (5) - (7)
            qD_ = (4e-10)^-1; % 1/m
            O = 2/pi.*( (cp-cv)./cp.*atan(qD_.*xi(idxC)) + cv./cp.*qD_.*xi(idxC) );
            O0 = 2/pi.*( 1-exp( -1./( (qD_.*xi(idxC)).^-1 + 1/3.*( qD_.*xi(idxC).*CO2.rhoc./rhoC ).^2 ) ) );
            lambda_c = rhoC.*cp.*1.02.*k.*TC ./ (6.*pi.*eta(idxC).*xi(idxC)) .* (O - O0) .*1e9;
            
            lambda(idxC) = lambda(idxC) + lambda_c;
        end % function criticalEnhancement_Huber
                
        function [eta, lambda] = criticalEnhancement_LS_Huber(eta, lambda, rho, T) % function criticalEnhancement_LS_Huber
            % according to Luettmer-Strathmann et al. (1995) for the 
            % viscosity and Huber et al. (2016) for the thermal
            % conductivity. Most of the preliminary calculation is identical.
            % Only the definition of the correlation length xi differs
            % slightly.
            % eta in mPa s
            % lambda in mW/m K
            % rho in kg/m
            % T in K
            
            % constants
            k = 1.3806485279e-23; % J/K - Bolzmann constant
            Tref = 3/2.*CO2.Tc; % K
            
            % prepare input for EOS to compute isothermal compressibility
            d = rho./CO2.rhoc;
            t = CO2.Tc./T;
            x = max(numel(d), numel(t));
            % add values for reference temperature
            if isscalar(t)
                t = repmat(t, size(d));
            end
            t = [t; repmat(CO2.Tc/Tref, size(d))];
            if ~isscalar(d)
                d = [d; d];
            end
            % isothermal compressibility from EOS
            phir = CO2.phir_master(d,t,2,0);
            drho_dp = 1./( CO2.Tc./t.*CO2.R.*1e-3.*( 1+ 2.*d.*phir{2,1} + d.^2.*phir{3,1}) );
            
            % reduced symetrized compressibility difference according to
            % Eq. (6.10) in Luettmer-Strathmann (1995) and Eq. (8) in Huber
            % (2016)
            dX = CO2.pc/CO2.rhoc^2 .* rho .* ...
                ( drho_dp(1:x) - Tref./T.*drho_dp(x+1:end) );
            % xi is bounded by 0; for dX <= 0, xi = 0
            % A Note on the Critical Enhancement of Transport Properties 
            % and Correlation Length of Fluids DOI 10.1007/s10765-013-1519-7
            % For xi = 0, critical enhancement of the transport properties
            % is zero. Consequently values of dX <= 0 can be omitted.
            idxC = dX > 0;
            
            % correlation length according to Eq. (6.10) in
            % Luettmer-Strathmann (1995); constants from Table II.
            % and according to Eq. (8) in Huber (2016); constants from the
            % text below Eq. (8).
            gn = 0.630/1.239; % critical exponent ratio
            xi = zeros(size(dX));
            xi(idxC) = 1.5e-10.*dX(idxC).^gn;            
            xi_LS = xi.*0.0481^-gn;
            xi_Huber = xi.*0.052^-gn;

            % equilibrium properties from EOS
            if isscalar(rho)
                rhoC = rho(any(idxC));
            else
                rhoC = rho(idxC);
            end
            if isscalar(T)
                TC = T(any(idxC));
            else
                TC = T(idxC);
            end
            t = CO2.Tc./TC;
            d = rhoC./CO2.rhoc;
            phi0 = CO2.phi0_master(d,t,0,2);
            phir = CO2.phir_master(d,t,2,2);
            cp = CO2.R.*( -t.^2.*( phi0{1,3}+phir{1,3} ) + ...
                ( 1 + d.*phir{2,1} - d.*t.*phir{2,2} ).^2 ./ ...
                ( 1 + 2.*d.*phir{2,1} + d.^2.*phir{3,1} ) );
            cv = CO2.R.*( -t.^2.*( phi0{1,3}+phir{1,3} ) );

            % parameters according to Table I in Luettmer-Strathmann (1995)
            qD = (2.3e-10)^-1; % 1/m
            yD = atan(qD .* xi_LS(idxC));
            foo = (1+(qD.*xi_LS(idxC)).^2).^(1/2);
            yd = ( atan(qD.*xi_LS(idxC)./foo) - yD )./foo;
            ya = rhoC.*k.*TC./( 8.*pi.*eta(idxC).^2.*xi_LS(idxC) ) .*1e6;
            yb = lambda(idxC)./( eta(idxC).*( cp - cv ) ) .* 1e-3;
            yg = cv./( cp - cv ); 
            % adjustments for asymptotic behavior
            eta_ = 2-1.239/0.630;
            alpha = 2-3.*0.630;
            yb = yb.^(1/(1-eta_/2));
            yg = yg.^(1/(1-eta_/2-alpha/2/0.630));
            % viscosity crossover function in Luettmer-Strathmann (1995)
            H = CO2.eta_crossover_LS(yD, yd, ya, yb, yg);
            
            % critical enhancement according to Eq. (5) - (7) in Huber
            % (2016)
            qD_ = (4e-10)^-1; % 1/m
            O = 2/pi.*( (cp-cv)./cp.*atan(qD_.*xi_Huber(idxC)) + cv./cp.*qD_.*xi_Huber(idxC) );
            O0 = 2/pi.*( 1-exp( -1./( (qD_.*xi_Huber(idxC)).^-1 + 1/3.*( qD_.*xi_Huber(idxC).*CO2.rhoc./rhoC ).^2 ) ) );
            lambda_c = rhoC.*cp.*1.02.*k.*TC ./ (6.*pi.*eta(idxC).*xi_Huber(idxC)) .* (O - O0) .*1e9;
            
            % critical enhancements
            eta(idxC) = eta(idxC).*exp(0.063.*H);
            lambda(idxC) = lambda(idxC) + lambda_c;
        end
        
    end % methods Transport auxiliary
    
end % classdef
    
%     properties(Constant, Hidden) 
%         % Coefficients for the calculation of phi0 from Table 27
%         phi0_a = [8.37304456 -3.70454304 2.5 1.99427042 0.62105248 0.41195293 1.04028922 0.08327678];
%         phi0_theta = [0 0 0 3.15163 6.11190 6.77708 11.32384 27.08792];
%         
%         % Coefficients for the calculation of phir from Table 31
%         phir_n = [...
%              0.38856823203161e0   0.29385475942740e1  -0.55867188534934e1  -0.76753199592477e0   0.31729005580416e0   0.54803315897767e0   0.12279411220335e0 ... 1-7
%              0.21658961543220e1   0.15841735109724e1  -0.23132705405503e0   0.58116916431436e-1 -0.55369137205382e0   0.48946615909422e0  -0.24275739843501e-1 ... 8-14
%              0.62494790501678e-1 -0.12175860225246e0  -0.37055685270086e0  -0.16775879700426e-1 -0.11960736637987e0  -0.45619362508778e-1  0.35612789270346e-1 ... 15-21
%             -0.74427727132052e-2 -0.17395704902432e-2 -0.21810121289527e-1  0.24332166559236e-1 -0.37440133423463e-1  0.14338715756878e0  -0.13491969083286e0 ... 22-28
%             -0.23151225053480e-1  0.12363125492901e-1  0.21058321972940e-2 -0.33958519026368e-3  0.55993651771592e-2 -0.30335118055646e-3 ... 29-34
%             -0.21365488688320e3   0.26641569149272e5  -0.24027212204557e5  -0.28341603423999e3   0.21247284400179e3... 35-39
%             -0.66642276540751e0   0.72608632349897e0   0.55068668612842e-1 ... 40-42
%             ];
%         phir_d = [...
%             1 1 1 1 2 2 3 ... 1-7
%             1 2 4 5 5 5 6 ... 8-14
%             6 6 1 1 4 4 4 ... 15-21
%             7 8 2 3 3 5 5 ... 22-28
%             6 7 8 10 4 8 ... 29-34
%             2 2 2 3 3 ... 35-39
%             ];
%         phir_t = [...
%             0.00 0.75 1.00 2.00 0.75 2.00 0.75 ... 1-7
%             1.50 1.50 2.50 0.00 1.50 2.00 0.00 ... 8-14
%             1.00 2.00 3.00 6.00 3.00 6.00 8.00 ... 15-21
%             6.00 0.00 7.00 12.0 16.0 22.0 24.0 ... 22-28
%             16.0 24.0 8.00 2.00 28.0 14.0 ... 29-34
%             1.00 0.00 1.00 3.00 3.00 ... 35-39
%             ];
%         phir_c = [...
%             1 1 1 1 1 1 1 ... 8-14
%             1 1 2 2 2 2 2 ... 15-21
%             2 2 3 3 3 4 4 ... 22-28
%             4 4 4 4 5 6 ... 29-34
%             ];
%         
%         phir_alpha = [25 25 25 15 20]; % 35-39
%         phir_beta = [325 300 300 275 275]; % 35-39
%         phir_gamma = [1.16 1.19 1.19 1.25 1.22]; % 35-39
%         phir_epsilon = [1 1 1 1 1]; % 35-39
%         
%         phir_a = [3.5 3.5 3.0]; % 40-42
%         phir_b = [0.875 0.925 0.875]; % 40-42
%         phir_bb = [0.3 0.3 0.3];  % 40-42
%         phir_A = [0.7 0.7 0.7]; % 40-42
%         phir_B = [0.3 0.3 1]; % 40-42
%         phir_C = [10 10 12.5]; % 40-42
%         phir_D = [275 275 275]; % 40-42
%     end

%     methods(Static, Hidden) % phi0
%         % Dimensionless helmholz function phi0 and its derivatives
%         % according to Table 28
%         % phi0 ... dimensionless helmholz function
%         % d ... reduced density delta = rho / rhoc
%         % t ... inverse reduced temperature tau = Tc / T
%         function val = phi0(d, t)
%             val = log(d) + CO2.phi0_a(1) + CO2.phi0_a(2).*t + CO2.phi0_a(3).*log(t) + ...
%                 sum( CO2.phi0_a(4:8) .* log( 1-exp(-CO2.phi0_theta(4:8).*t) ) ,2);
%         end
%         
%         function val = dphi0_dd(d, ~)
%             val = d.^-1;
%         end
%         
%         function val = ddphi0_dddd(d, ~)
%             val = -d.^-2;
%         end
%         
%         function val = ddphi0_dddt(d, ~)
%             val = zeros(size(d));
%         end
%         
%         function val = dphi0_dt(~, t)
%             e = exp(-CO2.phi0_theta(4:8).*t);
%             val = CO2.phi0_a(2) + CO2.phi0_a(3)./t + ...
%                 sum( CO2.phi0_a(4:8).*CO2.phi0_theta(4:8).*( (1-e).^-1 - 1 ) ,2);
%         end
%         
%         function val = ddphi0_dtdt(~, t)
%             e = exp(-CO2.phi0_theta(4:8).*t);
%             val = -CO2.phi0_a(3)./t.^2 - ...
%                 sum( CO2.phi0_a(4:8).*CO2.phi0_theta(4:8).^2.*e.*( 1-e ).^-2 ,2);
%         end
%     end
%     
%     methods(Static, Hidden) % phir
%         % Residual part of the dimensionless helmholz function phir and its derivatives
%         % according to Table 32
%         % phir ... residual part of the dimensionless helmholz function
%         % d ... reduced density delta = rho / rhoc
%         % t ... inverse reduced temperature tau = Tc / T
%         % Delta ... Distance function
%         % Theta ... Some other function
%         % Psi ... Exponential function
%         
%         function val = phir(d, t)
%             val = sum( CO2.phir_n(1:7) .* d.^CO2.phir_d(1:7) .* t.^CO2.phir_t(1:7) ,2) + ...
%                 sum( CO2.phir_n(8:34) .* d.^CO2.phir_d(8:34) .* t.^CO2.phir_t(8:34) .* exp(-d.^CO2.phir_c) ,2) + ...
%                 sum( CO2.phir_n(35:39) .* d.^CO2.phir_d(35:39) .* t.^CO2.phir_t(35:39) .* ...
%                     exp( -CO2.phir_alpha.*(d-CO2.phir_epsilon).^2 - CO2.phir_beta.*(t-CO2.phir_gamma).^2 ) ,2) + ...
%                 sum( CO2.phir_n(40:42) .* CO2.Delta(d, t).^CO2.phir_b .* d .* CO2.Psi(d, t), 2);
%         end
%         
%         function val = dphir_dd(d, t)
%             val = sum( CO2.phir_n(1:7) .* CO2.phir_d(1:7) .* d.^(CO2.phir_d(1:7)-1) .* t.^CO2.phir_t(1:7) ,2) + ...
%                 sum( CO2.phir_n(8:34).*exp(-d.^CO2.phir_c).*( d.^(CO2.phir_d(8:34)-1) .* t.^CO2.phir_t(8:34) .* ( CO2.phir_d(8:34) - CO2.phir_c.*d.^CO2.phir_c ) ) ,2) + ...
%                 sum( CO2.phir_n(35:39) .* d.^CO2.phir_d(35:39) .* t.^CO2.phir_t(35:39) .* ...
%                     exp( -CO2.phir_alpha.*(d-CO2.phir_epsilon).^2 - CO2.phir_beta.*(t-CO2.phir_gamma).^2 ) .* ...
%                     ( CO2.phir_d(35:39)./d - 2.*CO2.phir_alpha.*(d-CO2.phir_epsilon) ) ,2) + ...
%                 sum( CO2.phir_n(40:42) .* ( CO2.Delta(d,t).^CO2.phir_b.*( CO2.Psi(d,t) + d.*CO2.dPsi_dd(d,t) ) + ...
%                         CO2.dDeltab_dd(d,t) .* d .* CO2.Psi(d,t) ) ,2);
%         end
%         
%         function val = ddphir_dddd(d, t)
%             val = sum( CO2.phir_n(1:7) .* CO2.phir_d(1:7) .* (CO2.phir_d(1:7)-1) .* d.^(CO2.phir_d(1:7)-2) .* t.^CO2.phir_t(1:7) ,2) + ...
%                 sum( CO2.phir_n(8:34).*exp(-d.^CO2.phir_c).*( d.^(CO2.phir_d(8:34)-2) .* t.^CO2.phir_t(8:34) .* ( ...
%                     ( CO2.phir_d(8:34) - CO2.phir_c.*d.^CO2.phir_c ) .* ( CO2.phir_d(8:34) - 1 - CO2.phir_c.*d.^CO2.phir_c ) - ...
%                     CO2.phir_c.^2.*d.^CO2.phir_c )) ,2) + ...
%                 sum( CO2.phir_n(35:39) .* t.^CO2.phir_t(35:39) .* exp( -CO2.phir_alpha.*(d-CO2.phir_epsilon).^2 - CO2.phir_beta.*(t-CO2.phir_gamma).^2 ) .* ...
%                     ( -2.*CO2.phir_alpha.*d.^CO2.phir_d(35:39) + 4.*CO2.phir_alpha.^2.*d.^CO2.phir_d(35:39).*( d-CO2.phir_epsilon ).^2 - ...
%                     4.*CO2.phir_d(35:39).*CO2.phir_alpha.*d.^( CO2.phir_d(35:39)-1 ).*( d-CO2.phir_epsilon ) + CO2.phir_d(35:39).*( CO2.phir_d(35:39)-1 ).*d.^( CO2.phir_d(35:39)-2 ) ) ,2) + ...
%                 sum( CO2.phir_n(40:42) .* ( CO2.Delta(d,t).^CO2.phir_b.*( 2.*CO2.dPsi_dd(d,t) + d.*CO2.ddPsi_dddd(d,t) ) + ...
%                     2.*CO2.dDeltab_dd(d,t).*( CO2.Psi(d,t) + d.*CO2.dPsi_dd(d,t) ) + CO2.ddDeltab_dddd(d,t).*d.*CO2.Psi(d,t) ) ,2);
%         end
%         
%         function val = dphir_dt(d, t)
%             val = sum( CO2.phir_n(1:7) .* CO2.phir_t(1:7) .* d.^CO2.phir_d(1:7) .* t.^( CO2.phir_t(1:7)-1 ) ,2) + ...
%                 sum( CO2.phir_n(8:34) .* CO2.phir_t(8:34) .* d.^CO2.phir_d(8:34) .* t.^( CO2.phir_t(8:34)-1 ) .* exp( -d.^CO2.phir_c ) ,2) + ...
%                 sum( CO2.phir_n(35:39) .* d.^CO2.phir_d(35:39) .* t.^CO2.phir_t(35:39) .* ...
%                     exp( -CO2.phir_alpha.*(d-CO2.phir_epsilon).^2 - CO2.phir_beta.*(t-CO2.phir_gamma).^2 ) .* ...
%                     ( CO2.phir_t(35:39)./t - 2.*CO2.phir_beta.*(t-CO2.phir_gamma) ) ,2) + ...
%                 sum( CO2.phir_n(40:42).*d.*( CO2.dDeltab_dt(d,t).*CO2.Psi(d,t) + CO2.Delta(d,t).^CO2.phir_b.*CO2.dPsi_dt(d,t) ) ,2);
%         end
%         
%         function val = ddphir_dtdt(d, t)
%             val = sum( CO2.phir_n(1:7) .* CO2.phir_t(1:7) .* ( CO2.phir_t(1:7)-1 ) .* d.^CO2.phir_d(1:7) .* t.^( CO2.phir_t(1:7)-2 ) ,2) + ...
%                 sum( CO2.phir_n(8:34) .* CO2.phir_t(8:34) .* ( CO2.phir_t(8:34)-1 ) .* d.^CO2.phir_d(8:34) .* t.^( CO2.phir_t(8:34)-2 ) .* exp(-d.^CO2.phir_c) ,2) + ...
%                 sum( CO2.phir_n(35:39) .* d.^CO2.phir_d(35:39) .* t.^CO2.phir_t(35:39) .* ...
%                     exp( -CO2.phir_alpha.*(d-CO2.phir_epsilon).^2 - CO2.phir_beta.*(t-CO2.phir_gamma).^2 ) .* ...
%                     ( ( CO2.phir_t(35:39)./t - 2.*CO2.phir_beta.*(t-CO2.phir_gamma) ).^2 - CO2.phir_t(35:39)./t.^2 - 2.*CO2.phir_beta ) ,2)  + ...
%                 sum( CO2.phir_n(40:42).*d.*( CO2.ddDeltab_dtdt(d,t).*CO2.Psi(d,t) + 2.*CO2.dDeltab_dt(d,t).*CO2.dPsi_dt(d,t) + CO2.Delta(d,t).^CO2.phir_b.*CO2.ddPsi_dtdt(d,t) ) ,2);
%         end
%         
%         function val = ddphir_dddt(d, t)
%             val = sum( CO2.phir_n(1:7) .* CO2.phir_d(1:7) .* CO2.phir_t(1:7) .* d.^( CO2.phir_d(1:7)-1 ) .* t.^( CO2.phir_t(1:7)-1 ) ,2) + ...
%                 sum( CO2.phir_n(8:34).*exp(-d.^CO2.phir_c).*d.^( CO2.phir_d(8:34)-1 ).*CO2.phir_t(8:34).* t.^( CO2.phir_t(8:34)-1 ).*( CO2.phir_d(8:34) - CO2.phir_c.*d.^CO2.phir_c ) ,2) + ...
%                 sum( CO2.phir_n(35:39) .* d.^CO2.phir_d(35:39) .* t.^CO2.phir_t(35:39) .* ...
%                     exp( -CO2.phir_alpha.*(d-CO2.phir_epsilon).^2 - CO2.phir_beta.*(t-CO2.phir_gamma).^2 ) .* ...
%                     ( CO2.phir_d(35:39)./d - 2.*CO2.phir_alpha.*(d-CO2.phir_epsilon) ) .* ...
%                     ( CO2.phir_t(35:39)./t - 2.*CO2.phir_beta.*(t-CO2.phir_gamma) ) ,2) + ...
%                 sum( CO2.phir_n(40:42) .* (  CO2.Delta(d,t).^CO2.phir_b.*( CO2.dPsi_dt(d,t) + d.*CO2.ddPsi_dddt(d,t) ) + ...
%                     d.*CO2.dDeltab_dd(d,t).*CO2.dPsi_dt(d,t) + CO2.dDeltab_dt(d,t).*( CO2.Psi(d,t) + d.*CO2.dPsi_dd(d,t) ) + ...
%                     CO2.ddDeltab_dddt(d,t).*d.*CO2.Psi(d,t) ) ,2);
%         end
%         
%         function val = Delta(d, t)
%             val = CO2.Theta(d, t).^2 + CO2.phir_B .* ( (d-1).^2 ).^CO2.phir_a;
%         end
%         
%         function val = Theta(d, t)
%             val = (1-t) + CO2.phir_A .* ( (d-1).^2 ).^( 1./(2.*CO2.phir_bb) );
%         end
%         
%         function val = dDeltab_dd(d, t)
%             val = CO2.phir_b .* CO2.Delta(d, t).^(CO2.phir_b-1) .* CO2.dDelta_dd(d, t);
%             val(isnan(val)) = 0;
%         end
%         
%         function val = ddDeltab_dddd(d, t)
%             val = CO2.phir_b .* ( ...
%                 CO2.Delta(d, t).^(CO2.phir_b-1) .* CO2.ddDelta_dddd(d, t) + ...
%                 (CO2.phir_b-1) .* CO2.Delta(d, t).^(CO2.phir_b-2) .* CO2.dDelta_dd(d, t).^2 ...
%                 );
%             val(isnan(val)) = 0;
%         end
%         
%         function val = dDeltab_dt(d, t)
%             val = -2 .* CO2.Theta(d,t) .* CO2.phir_b .* CO2.Delta(d,t).^( CO2.phir_b-1 );
%             val(isnan(val)) = 0;
%         end
%         
%         function val = ddDeltab_dtdt(d, t)
%             val = 2 .* CO2.phir_b .* CO2.Delta(d,t).^( CO2.phir_b-1 ) + ...
%                 4 .* CO2.Theta(d, t).^2 .* CO2.phir_b .* (CO2.phir_b-1) .* CO2.Delta(d, t).^(CO2.phir_b-2);
%             val(isnan(val)) = 0;
%         end
%         
%         function val = ddDeltab_dddt(d, t)
%             val = -CO2.phir_A .* CO2.phir_b .* 2./CO2.phir_bb .* CO2.Delta(d, t).^(CO2.phir_b-1) .* (d-1) .* ( (d-1).^2 ).^( 1./(2.*CO2.phir_bb) - 1 ) - ...
%                 2 .* CO2.Theta(d, t) .* CO2.phir_b .* (CO2.phir_b-1) .* CO2.Delta(d, t).^(CO2.phir_b-2) .* CO2.dDelta_dd(d, t);
%             val(isnan(val)) = 0;
%         end
%         
%         function val = dDelta_dd(d, t)
%             val = (d-1) .* ( ...
%                 CO2.phir_A .* CO2.Theta(d, t) .* 2./CO2.phir_bb .* ( (d-1).^2 ).^( 1./(2.*CO2.phir_bb) - 1) + ...
%                 2 .* CO2.phir_B .* CO2.phir_a .* ( (d-1).^2 ).^(CO2.phir_a-1) ...
%                 );
%         end
%         
%         function val = ddDelta_dddd(d, t)
%             val = 1./(d-1) .* CO2.dDelta_dd(d, t) + (d-1).^2 .* ( ...
%                 4 .* CO2.phir_B .* CO2.phir_a .* ( CO2.phir_a-1 ) .* ( (d-1).^2 ).^( CO2.phir_a-2 ) + ...
%                 2 .* CO2.phir_A.^2 .* CO2.phir_bb.^-2 .* (( (d-1).^2 ).^( 1./(2.*CO2.phir_bb) - 1 )).^2 + ...
%                 CO2.phir_A .* CO2.Theta(d, t) .* 4 ./ CO2.phir_bb .* (1./(2.*CO2.phir_bb) -1) .* ( (d-1).^2 ).^( 1./(2.*CO2.phir_bb) - 2 ) ...
%                 );
%         end
%         
%         function val = Psi(d, t)
%             val = exp( -CO2.phir_C.*(d-1).^2 - CO2.phir_D.*(t-1).^2 );
%         end
%         
%         function val = dPsi_dd(d, t)
%             val = -2 .* CO2.phir_C .* (d-1) .* CO2.Psi(d, t);
%         end
%                 
%         function val = ddPsi_dddd(d, t)
%             val = ( 2.*CO2.phir_C.*(d-1).^2 - 1).*2.*CO2.phir_C.*CO2.Psi(d, t);
%         end
%                 
%         function val = dPsi_dt(d, t)
%             val = -2 .* CO2.phir_D .* (t-1) .* CO2.Psi(d, t);
%         end
%                 
%         function val = ddPsi_dtdt(d, t)
%             val = ( 2.*CO2.phir_D.*(t-1).^2 - 1).*2.*CO2.phir_D.*CO2.Psi(d, t);
%         end
%                 
%         function val = ddPsi_dddt(d, t)
%             val = 4.*CO2.phir_C.*CO2.phir_D.*(d-1).*(t-1).*CO2.Psi(d, t);
%         end
%     end
             
%         function d = d_fzero(t, d_ig, d0)
%             % Find the reduced density that satisfies p(d,t) = p
%             
%             % The corresponding problem is finding the zero point of
%             % f(d) = d + d^2.*dphir_dd - d_ig = 0
%             % t ... reduced temperature
%             % d_ig ... reduced density of the ideal gas
%             % d0 ... starting point
%             
%             % Two step algorithm:
%             % 1) Use newtons' algorithm to find an interval that contains
%             %   the zero point. i.e. sign(f(a)) ~= sign(f(b))
%             % 2) Use fzero type algorithm to find zero point
%             %   Based on Brent 1973
%             %   Algorithms for minimization without derivatives
% 
%             reltol = 1e-6;
%             abstol = 1e-4;
%             maxiter = 1e3;
%             k = 0;
%             
%             % Use newtons' algorithm to find an interval that contains
%             % the zero point. i.e. sign(f(a)) ~= sign(f(b))            
%             
%             % b is the best approximation to the zero point so far
%             % a is the value of the last step
%             a = d0;
%             b = d0;
%             dphir = CO2.dphir_dd(d0,t);
%             fb = d0 + d0^2.*dphir - d_ig;
%             fa = fb;
%             while sign(fb) == sign(fa) && k <= maxiter
%                 % Move b to a
%                 a = b;
%                 fa = fb;
%                 
%                 % Newton step f(x) / f'(x)
%                 delta = fa / ( 1+2.*a.*dphir + a.^2.*CO2.ddphir_dddd(a,t) );
%                 % Limit step
%                 if abs(delta)/a > 0.5
%                     delta = 0.5.*a.*sign(delta);
%                 end
%                 b = a - delta;
%                 
%                 % Compute new function value
%                 dphir = CO2.dphir_dd(b,t);
%                 fb = b + b^2.*dphir - d_ig;
%                 
%                 k = k+1;
%             end
% 
%             % Use fzero type algorithm to find zero point
%             %   Based on Brent 1973
%             %   Algorithms for minimization without derivatives
%             
%             % Target function
%             f = @(d) d+d.^2.*CO2.dphir_dd(d,t) - d_ig;
%             
%             % b is the best approximation to the zero point so far
%             % a is the value of the last step
%             % c is the best approximation where f(c) opposite sign of f(b)
%                         
%             % Auxilary variables
%             c = a;
%             fc = fa;
%             d = b - a;
%             e = d;
%             
%             tol = ( 2.*eps + reltol ).*abs(b);
%             m = 0.5 .* ( c-b );
%             
%             while ( abs(m) > tol || abs(fb) > abstol ) && k < maxiter
%                 % Make sure that f(b) is the value closest to zero
%                 if abs(fc) < abs(fb)
%                     a = b; b = c; c = a;
%                     fa = fb; fb = fc; fc = fa;
%                 end
%                 
%                 tol = ( 2.*eps + reltol ).*abs(b);
%                 m = 0.5 .* ( c-b );
%                 
%                 if abs(e) < tol && abs(fa) <= abs(fb)
%                     % Bisection
%                     d = m;
%                     e = m;
%                 else
%                     % Interpolation
%                     s = fb/fa;
%                     if a == c
%                         % Linear interpolation
%                         p = 2.*m.*s;
%                         q = 1-s;
%                     else
%                         % Inverse quadratic interpolation
%                         q = fa/fc;
%                         r = fb/fc;
%                         p = s.*( 2.*m.*q.*( q-r ) - ( b-a ).*( r-1) );
%                         q = ( q-1 ).*( r-1 ).*( s-1 );
%                     end
%                     
%                     if p > 0
%                         q = -q;
%                     else
%                         p = -p;
%                     end
%                     
%                     s = e;
%                     e = d;
%                     
%                     if 2.*p <  3.*m.*q-abs(tol.*q) || p < abs(0.5.*s.*q)
%                         d = p/q;
%                     else
%                         d = m;
%                         e = m;
%                     end
%                 end
%                 
%                 % Move b to a
%                 a = b;
%                 fa = fb;
%                 
%                 % Calculate new b
%                 if abs(d) > tol
%                     b = b + d;
%                 elseif m > 0
%                     b = b + tol;
%                 else
%                     b = b - tol;
%                 end
%                 fb = f(b);
%                 
%                 if ( fb>0 ) == ( fc>0 )
%                     c = a;
%                     fc = fa;
%                     d = b - a;
%                     e = d;
%                 end
%                 
%                 k = k + 1;
%             end % while
%             
%             % Return reduced density value
%             d = b;
%         end
