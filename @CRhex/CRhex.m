classdef CRhex < matlab.mixin.Copyable
    % method based on:
    % Crespi et al. doi: 10.1115/GT2017-64908 

    properties
        name
        k                   % devisions (number of nodes)

        process_medium_H    % available: "CO2"
        process_medium_C    % available: "CO2"

        HTcorr="dittusboelter";      % heat transfer correlation; available: "dittusboelter", "meshram", "kim", "ngo_zigzag", "ngo_sshapedfins", "zhao", "saeedkim", "cheng"
                % references
                % Meshram et al. http://dx.doi.org/10.1016/j.applthermaleng.2016.05.033
                % Kim et al. http://dx.doi.org/10.1016/j.anucene.2016.01.019
                % Saeed Kim https://doi.org/10.1016/j.enconman.2019.04.058
                % Cheng et al. https://doi.org/10.1016/j.applthermaleng.2021.116882
                % Jackson http://dx.doi.org/10.1016/j.nucengdes.2012.09.040; 
                % Pitla, Groll, Ramadhyani 2002 New correlation to predict the heat transfer coefficient during in-tube cooling of turbulent supercritical CO2, Int. Journal of Refrigeration
                % Gnielinski VDI Heat Atlas 2013, chapter G1
        hAratio             % conductance ratio

        design_case         % use function ondesigninput or format: table 
           % .Qtot  / W        transferred heat
           % .mdot_H  / kg/s   hot side mass flow rate
           % .T_Hin  / °C      hot side inlet temperature
           % .T_Hout  / °C     hot side outlet temperature
           % .P_Hin  / Pa      hot side inlet pressure
           % .P_Hout  / Pa     hot side outlet pressure                 
           % .mdot_C  / kg/s   cold side mass flow rate
           % .T_Cin  / °C      cold side inlet temperature
           % .T_Cout  / °C     cold side outlet temperature
           % .P_Cin  / Pa      cold side inlet pressure
           % .P_Cout  / Pa     cold side outlet pressure
        design              % storage table for solved design case, between nodes
           % .hA_H  / W/K           hot side thermal conductance
           % .hA_C  / W/K           cold side thermal conductance
           % .lambda_H  / W/mK      hot side thermal conductivity
           % .lambda_C  / W/mK      cold side thermal conductivity
           % .Re_H  /               hot side "Reynolds number" without area
           % .Re_C  /               cold side "Reynolds number" without area
           % .Pr_H  / -             hot side Prandtl number
           % .Pr_C  / -             cold side Prandtl number
           % .T_H  / °C             hot side temperature
           % .T_C  / °C             cold side temperature
           % .P_H   / Pa            hot side pressure
           % .P_C   / Pa            cold side pressure
           % .TD  / K               mean temperature difference
           % .UA  / W/m             overall heat transfer coefficient
           % .dQ  / W               transferred heat
        offdesign_cases     % use function ondesigninput or format: table
                            % see design_case for variables
        offdesign_solv      % storage for solved offdesign cases

        T_H  % °C           solved design case: hot side temperature, per node
        T_C  % °C           solved design case: cold side temperature, per node
        P_H  % Pa           solved design case: hot side pressure, per node               
        P_C  % Pa           solved design case: cold side pressure, per node  
        h_H  % J/kgK        solved design case: hot side enthalpy, per node  
        h_C  % J/kgK        solved design case: cold side enthalpy, per node 

        ispinch
    end

    properties (Hidden)
        % calculation parameters
        tol        % tolerance, stop criterion iterative offdesign calculation


        db_x=0.8;           % Dittus-Boelter coefficients
        db_yH=0.3;
        db_yC=0.4;
        k_m1_H=0.8138;      % Kim coefficients
        k_m1_C=0.8742;
        m_m1_H=0.893;       % Meshram coefficients
        m_m1_H2=0.869;
        m_m1_C=0.871;
        m_m1_C2=0.876;
        m_m2=0.7;
        z_m1_H=0.6051;       % Zhao coefficients
        z_m1_C=0.6637;
        z_m2_H=0.2127;
        z_m2_C=0.3536;
        n_z_m1=0.629;        % Ngo coefficients (zigzag channels)
        n_z_m2=0.317;
        n_s_m1=0.593;        % Ngo coefficients (s-shaped fins)
        n_s_m2=0.43;
        s_m1=0.83;           % Saeed & Kim coefficients
        s_m2=0.95;
        c_m1_H=0.76214;      % Cheng coefficients
        c_m1_C=0.7678;
        cu_m1=0.83;          % custom coefficients
        cu_m2=0.95;

        
        prop_H              % references to property class, hot fluid
        prop_C              % references to property class, cold fluid
        
        % performance indicators
        calctime      % elapsed time
        count_it      % count of iterations
    end

    methods
        function obj = CRhex(process_medium_H, process_medium_C)
            if nargin==0
                return
            end
            
            % assigns property class automatically
            if strcmp(process_medium_H, "CO2")
                obj.process_medium_H = "CO2";
                obj.prop_H = CO2();
            end

            if strcmp(process_medium_C, "CO2")
                obj.process_medium_C = "CO2";
                obj.prop_C = CO2();
            end
        end

        function ondesigninput(obj, mdot_H, T_Hin, T_Hout, P_Hin, P_Hout, mdot_C, T_Cin, T_Cout, P_Cin, P_Cout)
             obj.design_case = table('Size',[1 11],'VariableTypes',repmat({'double'},1,11),'VariableNames',{'Qtot','mdot_H','T_Hin','T_Hout','P_Hin','P_Hout','mdot_C','T_Cin','T_Cout','P_Cin','P_Cout'});
             obj.design_case.mdot_H = mdot_H;
             obj.design_case.T_Hin = T_Hin;
             obj.design_case.T_Hout = T_Hout;
             obj.design_case.P_Hin = P_Hin;
             obj.design_case.P_Hout = P_Hout;
             obj.design_case.mdot_C = mdot_C;
             obj.design_case.T_Cin = T_Cin;
             obj.design_case.T_Cout = T_Cout;
             obj.design_case.P_Cin = P_Cin;
             obj.design_case.P_Cout = P_Cout;
        end

        function obj = ondesign(obj,hAratio,k)
             obj.hAratio=hAratio;
             obj.k=k;
             % storage container for saved design case
             obj.design= table('Size',[obj.k-1 16],'VariableTypes',repmat({'double'},1,16),'VariableNames',{'hA_H', 'hA_C','lambda_H','lambda_C', 'rho_H', 'rho_C', 'Re_H', 'Re_C', 'Pr_H', 'Pr_C', 'T_H', 'T_C', 'P_H', 'P_C','TD','UA'});
             
             % Qtot calculation
             h_Hin=obj.prop_H.h(obj.design_case.P_Hin,obj.design_case.T_Hin+273.15);
             h_Hout=obj.prop_H.h(obj.design_case.P_Hout,obj.design_case.T_Hout+273.15);
             
             obj.design_case.Qtot=obj.design_case.mdot_H*(h_Hin-h_Hout);
             dQ= obj.design_case.Qtot/(obj.k-1);
             obj.design.dQ=dQ*ones(obj.k-1,1); 

             % properties along nodes
                % assumption: linear pressure drop
             obj.P_H=linspace(obj.design_case.P_Hin, obj.design_case.P_Hout, obj.k);
             obj.P_C=linspace(obj.design_case.P_Cout, obj.design_case.P_Cin, obj.k);

             % spacing: linear enthalpy distribution
             obj.h_H=linspace(h_Hin, h_Hout, obj.k);
            
             h_Cin=obj.prop_C.h(obj.design_case.P_Cin,obj.design_case.T_Cin+273.15);
             h_Cout=obj.prop_C.h(obj.design_case.P_Cout,obj.design_case.T_Cout+273.15);
             
             obj.h_C=linspace(h_Cout, h_Cin, obj.k);
                
             % T=f(p,h)
             obj.T_H=obj.prop_H.T_ph(obj.P_H,obj.h_H)-273.15;
             obj.T_C=obj.prop_C.T_ph(obj.P_C,obj.h_C)-273.15;
                % Pinch Point Check
             obj.ispinch=min(obj.T_H-obj.T_C)<5;

             % mean values in each division (k-1)
                % temperature difference, pressure   
             for i=1:obj.k-1
                obj.design.T_H(i)=(obj.T_H(i)+obj.T_H(i+1))/2;
                obj.design.T_C(i)=(obj.T_C(i)+obj.T_C(i+1))/2;
                obj.design.P_H(i)=(obj.P_H(i)+obj.P_H(i+1))/2;
                obj.design.P_C(i)=(obj.P_C(i)+obj.P_C(i+1))/2;
             end
             obj.design.TD= obj.design.T_H-obj.design.T_C;
             obj.design.UA= dQ./obj.design.TD;
                % use hAratio to obtain 
             obj.design.hA_H= obj.design.UA.*(1+obj.hAratio);
             obj.design.hA_C= obj.design.UA./obj.hAratio.*(1+obj.hAratio);
                
             % calculating all properties
             [obj.design.rho_H,obj.design.rho_C,obj.design.lambda_H,obj.design.lambda_C,obj.design.Pr_H,obj.design.Pr_C,obj.design.Re_H,obj.design.Re_C]=CRhex.calctherm(obj,obj.design);
        end

        function offdesigninput(obj, mdot_H, T_Hin, P_Hin, mdot_C, T_Cin, P_Cin)
            obj.offdesign_cases = table('Size',[length(mdot_H) 11],'VariableTypes',repmat({'double'},1,11),'VariableNames',{'Qtot','mdot_H','T_Hin','T_Hout','P_Hin','P_Hout','mdot_C','T_Cin','T_Cout','P_Cin','P_Cout'});
             obj.offdesign_cases.mdot_H = mdot_H;
             obj.offdesign_cases.T_Hin = T_Hin;
             obj.offdesign_cases.P_Hin = P_Hin;
             obj.offdesign_cases.mdot_C = mdot_C;
             obj.offdesign_cases.T_Cin = T_Cin;
             obj.offdesign_cases.P_Cin = P_Cin;
             obj.offdesign_solv=cell(length(mdot_H),1);
             obj.calctime=NaN(length(mdot_H),1);
        end
        
        function obj = offdesign(obj,tol,n_case)
            tic
            obj.tol=tol;
            obj.count_it=1;
            % initialization of variables before while-loop
             dQ=100000;
             %dQrel=10;
             Qx=ones(1,obj.k-1); %between nodes
             T_Hx=NaN(1,obj.k-1);
             T_Cx=NaN(1,obj.k-1);
             dP_Hx=NaN(1,obj.k-1);
             dP_Cx=NaN(1,obj.k-1); 
             P_Hx=NaN(1,obj.k-1);
             P_Cx=NaN(1,obj.k-1);
             Re_Hx=NaN(1,obj.k-1);
             Re_Cx=NaN(1,obj.k-1);
             T_H=NaN(1,obj.k); %per node
             T_C=NaN(1,obj.k);               
             P_H=NaN(1,obj.k);
             P_C=NaN(1,obj.k);               
             h_H=NaN(1,obj.k);
             h_C=NaN(1,obj.k);

             for i=1:obj.k-1
                dP_H(i)=obj.P_H(i)-obj.P_H(i+1);
                dP_C(i)=-obj.P_C(i)+obj.P_C(i+1);
             end       
             P_H(1)=obj.offdesign_cases.P_Hin(n_case);
             P_C(end)=obj.offdesign_cases.P_Cin(n_case);

             P_H(2:end)=P_H(1)-cumsum(dP_H); 
             P_C(1:end-1)=P_C(end)-flip(cumsum(flip(dP_C))); 

             %initial conditions: TQ-graph. Assume an approach temperature
             %and check if temperature curves are crossing in the TQ-graph.
             %If so, increase the approach temperature and try again.
             h_H(1)=obj.prop_H.h(P_H(1),obj.offdesign_cases.T_Hin(n_case)+273.15); 
             h_C(end)=obj.prop_C.h(P_C(end),obj.offdesign_cases.T_Cin(n_case)+273.15);

             dTstart=10;
             TQok=false;
             while ~TQok
                 T_H(end)=obj.offdesign_cases.T_Cin(n_case)+dTstart; 
                 h_H(end)=obj.prop_H.h(P_H(end),T_H(end)+273.15);
                 h_H=linspace(h_H(1),h_H(end),obj.k);

                 Qdot=obj.offdesign_cases.mdot_H(n_case)*(h_H(1)-h_H(end));

                 h_C(1)=Qdot./obj.offdesign_cases.mdot_C(n_case)+h_C(end);
                 h_C=linspace(h_C(1),h_C(end),obj.k);

                 T_H=obj.prop_H.T_ph(P_H,h_H)-273.15;
                 T_C=obj.prop_C.T_ph(P_C,h_C)-273.15;
                
                 if all(T_H>T_C)
                     TQok=true;
                 else
                     dTstart=dTstart+10;
                 end
             end


                
             % iterative calculation
             while dQ>tol %|| dQrel>0.005
                
               % per divisions
                 for i=1:obj.k-1
                    T_Hx(i)=(T_H(i)+T_H(i+1))/2;
                    T_Cx(i)=(T_C(i)+T_C(i+1))/2;
                    P_Hx(i)=(P_H(i)+P_H(i+1))/2;
                    P_Cx(i)=(P_C(i)+P_C(i+1))/2;
                 end
                             
                 % other properties depending on T, p
                 rho_Hx=obj.prop_H.rho(P_Hx,T_Hx+273.15);
                 rho_Cx=obj.prop_C.rho(P_Cx,T_Cx+273.15);
                 lambda_Hx=obj.prop_H.lambda(rho_Hx,T_Hx+273.15);
                 lambda_Cx=obj.prop_C.lambda(rho_Cx,T_Cx+273.15);
                 Pr_Hx=obj.prop_H.Pr(rho_Hx,T_Hx+273.15);
                 Pr_Cx=obj.prop_C.Pr(rho_Cx,T_Cx+273.15);
                 my_Hx=obj.prop_H.my(rho_Hx,T_Hx+273.15);
                 my_Cx=obj.prop_C.my(rho_Cx,T_Cx+273.15);
                 Re_Hx=obj.offdesign_cases.mdot_H(n_case)./my_Hx;
                 Re_Cx=obj.offdesign_cases.mdot_C(n_case)./my_Cx;
                 
             switch obj.HTcorr
                 case "dittusboelter"
                     %acc. to Hoopes, Crespi
                 hA_Hx=obj.design.hA_H'.*(lambda_Hx./obj.design.lambda_H').*(Re_Hx./obj.design.Re_H').^obj.db_x.*(Pr_Hx./obj.design.Pr_H').^obj.db_yH;
                 hA_Cx=obj.design.hA_C'.*(lambda_Cx./obj.design.lambda_C').*(Re_Cx./obj.design.Re_C').^obj.db_x.*(Pr_Cx./obj.design.Pr_C').^obj.db_yC;
                 
                 case "kim"
                       % acc. to Kim
                 hA_Hx=obj.design.hA_H'.*(lambda_Hx./obj.design.lambda_H').*(Re_Hx./obj.design.Re_H').^obj.k_m1_H;
                 hA_Cx=obj.design.hA_C'.*(lambda_Cx./obj.design.lambda_C').*(Re_Cx./obj.design.Re_C').^obj.k_m1_C;

                 case "meshram"
                     %acc. to Meshram
                 isTlow= T_Hx<(580-273.15);
                 isTmiddle= T_Hx>=(580-273.15) & T_Hx<=(630-273.15);
                 isThigh= T_Hx>(630-273.15);
                 hA_Hx_low=obj.design.hA_H'.*(lambda_Hx./obj.design.lambda_H').*(Re_Hx./obj.design.Re_H').^obj.m_m1_H.*(Pr_Hx./obj.design.Pr_H').^obj.m_m2;
                 hA_Hx_high=obj.design.hA_H'.*(lambda_Hx./obj.design.lambda_H').*(Re_Hx./obj.design.Re_H').^obj.m_m1_H2.*(Pr_Hx./obj.design.Pr_H').^obj.m_m2;
                 hA_Hx(isTlow)=hA_Hx_low(isTlow);
                 hA_Hx(isTmiddle)=(hA_Hx_low(isTmiddle)+hA_Hx_high(isTmiddle))./2;
                 hA_Hx(isThigh)=hA_Hx_high(isThigh);

                 isTlow= T_Hx<(500-273.15);
                 isTmiddle= T_Hx>=(500-273.15) & T_Hx<=(520-273.15);
                 isThigh= T_Hx>(520-273.15);
                 hA_Cx_low=obj.design.hA_C'.*(lambda_Cx./obj.design.lambda_C').*(Re_Cx./obj.design.Re_C').^obj.m_m1_C.*(Pr_Cx./obj.design.Pr_C').^obj.m_m2;
                 hA_Cx_high=obj.design.hA_C'.*(lambda_Cx./obj.design.lambda_C').*(Re_Cx./obj.design.Re_C').^obj.m_m1_C2.*(Pr_Cx./obj.design.Pr_C').^obj.m_m2;
                 hA_Cx(isTlow)=hA_Cx_low(isTlow);
                 hA_Cx(isTmiddle)=(hA_Cx_low(isTmiddle)+hA_Cx_high(isTmiddle))./2;
                 hA_Cx(isThigh)=hA_Cx_high(isThigh);

                 case "ngo_zigzag"
                     %acc. to Ngo
                 hA_Hx=obj.design.hA_H'.*(lambda_Hx./obj.design.lambda_H').*(Re_Hx./obj.design.Re_H').^obj.n_z_m1.*(Pr_Hx./obj.design.Pr_H').^obj.n_z_m2;
                 hA_Cx=obj.design.hA_C'.*(lambda_Cx./obj.design.lambda_C').*(Re_Cx./obj.design.Re_C').^obj.n_z_m1.*(Pr_Cx./obj.design.Pr_C').^obj.n_z_m2;

                 case "ngo_sshapedfins"
                     %acc. to Ngo
                 hA_Hx=obj.design.hA_H'.*(lambda_Hx./obj.design.lambda_H').*(Re_Hx./obj.design.Re_H').^obj.n_s_m1.*(Pr_Hx./obj.design.Pr_H').^obj.n_s_m2;
                 hA_Cx=obj.design.hA_C'.*(lambda_Cx./obj.design.lambda_C').*(Re_Cx./obj.design.Re_C').^obj.n_s_m1.*(Pr_Cx./obj.design.Pr_C').^obj.n_s_m2;

                 case "zhao"
                     %acc. to Zhao
                 hA_Hx=obj.design.hA_H'.*(lambda_Hx./obj.design.lambda_H').*(Re_Hx./obj.design.Re_H').^obj.z_m1_H.*(Pr_Hx./obj.design.Pr_H').^obj.z_m2_H;
                 hA_Cx=obj.design.hA_C'.*(lambda_Cx./obj.design.lambda_C').*(Re_Cx./obj.design.Re_C').^obj.z_m1_C.*(Pr_Cx./obj.design.Pr_C').^obj.z_m2_C;

                 case "saeedkim"
                 hA_Hx=obj.design.hA_H'.*(lambda_Hx./obj.design.lambda_H').*(Re_Hx./obj.design.Re_H').^obj.s_m1.*(Pr_Hx./obj.design.Pr_H').^obj.s_m2;
                 hA_Cx=obj.design.hA_C'.*(lambda_Cx./obj.design.lambda_C').*(Re_Cx./obj.design.Re_C').^obj.s_m1.*(Pr_Cx./obj.design.Pr_C').^obj.s_m2;

                 case "cheng"
                 hA_Hx=obj.design.hA_H'.*(lambda_Hx./obj.design.lambda_H').*(Re_Hx./obj.design.Re_H').^obj.c_m1_H;
                 hA_Cx=obj.design.hA_C'.*(lambda_Cx./obj.design.lambda_C').*(Re_Cx./obj.design.Re_C').^obj.c_m1_C;

                 case "custom"
                 hA_Hx=obj.design.hA_H'.*(lambda_Hx./obj.design.lambda_H').*(Re_Hx./obj.design.Re_H').^obj.cu_m1.*(Pr_Hx./obj.design.Pr_H').^obj.cu_m2;
                 hA_Cx=obj.design.hA_C'.*(lambda_Cx./obj.design.lambda_C').*(Re_Cx./obj.design.Re_C').^obj.cu_m1.*(Pr_Cx./obj.design.Pr_C').^obj.cu_m2;

                 otherwise
                     disp("wrong input for heat transfer correlavtion HTcorr")
             end

                 UAx=1./(1./hA_Hx+1./hA_Cx);
                 TDx=T_Hx-T_Cx;
                 Qxnew=UAx.*TDx;

                 % stop criterion
                 dQnew=max(abs(Qx-Qxnew));
                 Qx=(Qxnew+Qx)./2;
                 dQ=dQnew;


                 h_H(2:end)=h_H(1:end-1)-Qx./obj.offdesign_cases.mdot_H(n_case);
                 h_C(1:end-1)=h_C(2:end)+Qx./obj.offdesign_cases.mdot_C(n_case);
                 
                T_Hnew=(obj.prop_H.T_ph(P_H,h_H)-273.15);
                T_Cnew=(obj.prop_C.T_ph(P_C,h_C)-273.15);
                %old stop criterion: dTH=max(max(abs(T_Hnew-T_H)),max(abs(T_Cnew-T_C)));

                isneg=(T_Hnew-T_Cnew)<0;

                if any(isneg)
                    warning('T-profiles cross')
                end
                
                % pressure drop scaling formulas
                dP_Hx=dP_H.*((obj.offdesign_cases.mdot_H(n_case)^2./rho_Hx)./(obj.design_case.mdot_H^2./obj.design.rho_H'));
                dP_Cx=dP_C.*((obj.offdesign_cases.mdot_C(n_case)^2./rho_Cx)./(obj.design_case.mdot_C^2./obj.design.rho_C'));

                % per node
                P_H(2:end)=P_H(1)-cumsum(dP_Hx); 
                P_C(1:end-1)=P_C(end)-flip(cumsum(flip(dP_Cx))); 
                T_H=T_Hnew; 
                T_C=T_Cnew;  
                
                
                obj.count_it=obj.count_it+1;

                % %%use to check convergence
                % a=400;
                % mapcolors=jet(a);
                % mapcolors_w=colormap(winter(a));
                % mapcolors_a=colormap(autumn(a));

                if mod(obj.count_it,10)==0
                    % %%use to check convergence
                    % figure(1)
                    % hold on
                    % ylabel('T / °C')
                    % plot(T_H,'LineStyle','-','Color', mapcolors_a(obj.count_it,:))
                    % plot(T_C,'LineStyle','-','Color', mapcolors_w(obj.count_it,:))
                    % plot(obj.T_H,'k-')
                    % plot(obj.T_C,'k-')
                    % 
                    % figure(2)
                    % hold on
                    % ylabel('dQ / W')
                    % plot(obj.design.dQ,'k-')
                    % plot(Qx,'LineStyle','-','Color',mapcolors(obj.count_it,:))

                    disp("case "+n_case+": iteration "+obj.count_it+", error "+dQ+"W")
                end

                if obj.count_it>10000
                    warning('Calculation terminated.')
                    break
                end

             end

             obj.calctime(n_case)=toc;
             disp("finished off-design case no. "+n_case+": "+obj.count_it+" iterations in "+obj.calctime(n_case)+"s, error "+dQ+"W") %disp("finished off-design case no. "+n_case+": "+obj.count_it+" iterations in "+obj.calctime(n_case)+"s, error "+dQ+"W"+", relative error "+dQrel)
             obj.offdesign_cases.Qtot(n_case)=sum(Qx);
             obj.offdesign_cases.T_Hout(n_case)=T_H(end);
             obj.offdesign_cases.P_Hout(n_case)=P_H(end);
             obj.offdesign_cases.T_Cout(n_case)=T_C(1);
             obj.offdesign_cases.P_Cout(n_case)=P_C(1);
             obj.offdesign_solv(n_case,1)={table(T_H',T_C',P_H',P_C',h_H',h_C','VariableNames',["T_H","T_C","P_H","P_C","h_H","h_C"])};
             obj.offdesign_solv(n_case,2)={table(hA_Hx',hA_Cx',UAx',TDx',Qx',dP_Hx',dP_Cx','VariableNames',["hA_Hx","hA_Cx","UAx","TDx","Qx","dP_Hx","dP_Cx"])};

        end
    end

    methods (Static)
        function [rho_H,rho_C,lambda_H,lambda_C,Pr_H,Pr_C,Re_H,Re_C]=calctherm(obj,objtable)
            rho_H=obj.prop_H.rho(objtable.P_H,objtable.T_H+273.15);
            rho_C=obj.prop_C.rho(objtable.P_C,objtable.T_C+273.15);
            my_H=obj.prop_H.my(rho_H,objtable.T_H+273.15);
            my_C=obj.prop_C.my(rho_C,objtable.T_C+273.15);
            lambda_H=obj.prop_H.lambda(rho_H,objtable.T_H+273.15);
            lambda_C=obj.prop_C.lambda(rho_C,objtable.T_C+273.15);
            Pr_H=obj.prop_H.Pr(rho_H,objtable.T_H+273.15);
            Pr_C=obj.prop_C.Pr(rho_C,objtable.T_C+273.15);
            
            % not real dimensionless Reynoldsnumber but without
            % characteristic length L and A for velocity calculation!
            Re_H=obj.design_case.mdot_H./my_H;
            Re_C=obj.design_case.mdot_C./my_C;

        end
    end

end