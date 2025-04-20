 % villyes, 21.11.24
 % MATLAB R2024b

% Example call for CRhex-class-call; for PCHE-type heat exchangers
% make sure to install Parallel Computing Toolbox

%% INPUTS: conductance ratio method specific
k_devisions=200;
dQ=0.05;
hAratio=[1,0.8];
HTcorr=["dittusboelter", "meshram", "kim", "ngo_zigzag"];
ref_case=1;

%% INPUTS: example data for recalculation, SI-units
% if using experimental data and later comparing, make sure to make the
% data fit the energy conversion on both sides exactly.

tab=table;
tab.id=["dataset_1"; "dataset_2"];
tab.mCO2=[0.6;0.4]; %kg/s (only given once, since recuperator.)
tab.THotIn=[360;450]+273.15; %K
tab.pHotIn=[61*10^5;60*10^5]; %Pa (absolute pressure, if using experimental data make sure to add ambient p)
tab.TColdIn=[20;20]; %°C
tab.pColdIn=[150*10^5;200*10^5]; %Pa

% for the reference case, make sure the data fits together (both sides produce the same transferred heat)
tab.THotOut(ref_case)=88.6621;
tab.pHotOut(ref_case)=61*10^5;
tab.TColdOut(ref_case)=147.5767;
tab.pColdOut(ref_case)=150*10^5;


%% PREPARATIONS
l_r=length(hAratio);
l_c=length(HTcorr);

pche1=CRhex("CO2","CO2");
pche1.name="zigzag";
pche1.ondesigninput(tab.mCO2(ref_case), tab.THotIn(ref_case), tab.THotOut(ref_case), tab.pHotIn(ref_case), tab.pHotOut(ref_case), tab.mCO2(ref_case), tab.TColdIn(ref_case), tab.TColdOut(ref_case), tab.pColdIn(ref_case), tab.pColdOut(ref_case));
pche1.offdesigninput(tab.mCO2, tab.THotIn, tab.pHotIn, tab.mCO2, tab.TColdIn, tab.pColdIn);

pche=repmat(CRhex,[l_r,l_c]);
for i=1:l_c*l_r
    pche(i)=copy(pche1);
end

for i=2:l_c
    for j=1:l_r
        pche(j,i).HTcorr=HTcorr(i);
    end
end

clear pche1

%% ON-DESIGN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dsolv=cell(l_r,1);
dsolv2=cell(l_r,1);
TH=cell(l_r,1);
TC=cell(l_r,1);
pH=cell(l_r,1);
pC=cell(l_r,1);
hH=cell(l_r,1);
hC=cell(l_r,1);
for j=1:l_c
    for i=1:l_r
        local=pche(i,j);
        local.ondesign(hAratio(i), k_devisions);
        dsolv{i}=local.design;
        dsolv2{i}=local.design_case;
        TH{i}=local.T_H;
        TC{i}=local.T_C;
        pH{i}=local.P_H;
        pC{i}=local.P_C;
        hH{i}=local.h_H;
        hC{i}=local.h_C;
    end
    
    for i=1:l_r
        pche(i,j).hAratio=hAratio(i);
        pche(i,j).k=k_devisions;
    
        pche(i,j).T_H=TH{i};
        pche(i,j).T_C=TC{i};
        pche(i,j).P_H=pH{i};
        pche(i,j).P_C=pC{i};
        pche(i,j).h_H=hH{i};
        pche(i,j).h_C=hC{i};
    
        pche(i,j).design=dsolv{i};
        pche(i,j).design_case=dsolv2{i};
    end
end

disp('Design case calculation finished.')
% disp(datetime)

clear TH TC pH pC hH hC dsolv dsolv2

%%
% off-design
htab=height(tab);

%% OFF-DESIGN
for i=1:l_c*l_r
    tab_offdesign=cell(htab,1);
    solv=cell(htab,2);
    calctime=NaN(htab,1);
    count=NaN(htab,1);

    parfor j=1:htab
        local=pche(i);
        local.offdesign(dQ,j);
        tab_offdesign{j,:}=local.offdesign_cases(j,:);
        solv(j,:)=local.offdesign_solv(j,:);
        calctime(j)=local.calctime(j);
        count(j)=local.count_it;
    end

    for j=1:htab
        pche(i).offdesign_cases(j,:)=tab_offdesign{j,:};
        pche(i).calctime=calctime;
        pche(i).count_it=count;
        pche(i).offdesign_solv=solv;
    end
    
    disp("DONE: pche calculation "+i)
    % save('savePCHE')
end

