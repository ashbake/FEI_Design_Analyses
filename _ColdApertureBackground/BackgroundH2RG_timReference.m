clear all;

h=6.62*10^-34 ; % J.s
c=299*10^6; % m.s-1
kb=1.38*10^-23; % J.K-1

param_pixel_size=18*18*10^-12; % m^2
param_size_H2RG=0.036; % m
param_distance_sphere_sensor=0.4;% m

param_H2RG_A_CutOff=2.5*10^-6; % cut off wavelenght of the H2RG
param_H2RG_SNAP_CutOff=1.7*10^-6; % cut off wavelenght of the H2RG SNAP

lambda=(100:1:10000) *10^-9; % m

%% geometry



%% Black Body radiation
T_BB=290;
param_BB_size=(0.000018)^2; % m-2
BB=(2*c./lambda.^4)./(exp(h*c./(lambda*kb*T_BB))-1)*param_BB_size; %ph.s-1.m-1.sr-1


%% Transmittance Filter


param_Taverage=0.95; % percentage
param_pass=[830 940] *10^-9; % bandwidth below 1percent
OD4=0.0001;
T_Filter=ones(size(lambda,2),1)*OD4;
for i=1:size(lambda,2);
    if lambda(i)< param_pass(1);
         T_Filter(i)=0.0001;
    end
    if lambda(i)>param_pass(1) && lambda(i)<param_pass(2) ;
        T_Filter(i)=param_Taverage;
    end
    if lambda(i)>param_pass(2);
        T_Filter(i)=0.0001;    
    end
end

T_Filter_Fast=ones(size(lambda,2),1)*OD4;
for i=1:size(lambda,2);
    if lambda(i)< param_pass(1);
         T_Filter_Fast(i)=OD4;
    end
    if lambda(i)>param_pass(1) && lambda(i)<param_pass(2) ;
        T_Filter_Fast(i)=param_Taverage;
    end
    if lambda(i)>param_pass(2);
        T_Filter_Fast(i)=OD4;    
    end
end

%% photon current over one pixel of the H2RG

% CASE A: Filter F/0.5 and Cold Baffle F/1.2
FN_ColdBaffle=1.2;
SolidAngle_ColdBaffle=atan(1/2/FN_ColdBaffle)*2;
FN_Filter=3;
SolidAngle_Filter=atan(1/2/FN_Filter)*2;
R=BB*(SolidAngle_ColdBaffle-SolidAngle_Filter).*T_Filter_Fast'+BB*SolidAngle_Filter.*T_Filter';
for i=1:size(lambda,2)
Background(i)=sum(R(1:i),2)*10^-9;
end
figure; semilogy(lambda*10^6,Background,'-.m');ylim([1E-9 1E5]);xlim([0.5 3.00]);
hold on;
% 
% CASE B: Filter F/0.5 and Cold Baffle F/8
FN_ColdBaffle=8;
SolidAngle_ColdBaffle=atan(1/2/FN_ColdBaffle)*2;
FN_Filter=0.5;
SolidAngle_Filter=atan(1/2/FN_Filter)*2;
R=BB*(SolidAngle_ColdBaffle-SolidAngle_Filter).*T_Filter_Fast'+BB*SolidAngle_Filter.*T_Filter';
for i=1:size(lambda,2)
Background(i)=sum(R(1:i),2)*10^-9;
end
semilogy(lambda*10^6,Background,'-.r');
hold on;
% % CASE C: Filter F/1.4 and Cold Baffle F/1.4
% FN_ColdBaffle=1.4;
% SolidAngle_ColdBaffle=atan(1/2/FN_ColdBaffle)*2;
% FN_Filter=1.4;
% SolidAngle_Filter=atan(1/2/FN_Filter)*2;
% R=BB*(SolidAngle_ColdBaffle-SolidAngle_Filter).*T_Filter_Fast'+BB*SolidAngle_Filter.*T_Filter';
% for i=1:size(lambda,2)
% Background(i)=sum(R(1:i),2)*10^-9;
% end
% semilogy(lambda*10^6,Background,'-.g');
% hold on;
% 
% % CASE D: Pinhole Filter F/2.8 and Cold Baffle F/2.8
% FN_ColdBaffle=2.8;
% SolidAngle_ColdBaffle=atan(1/2/FN_ColdBaffle)*2;
% FN_Filter=2.8;
% SolidAngle_Filter=atan(1/2/FN_Filter)*2;
% R=BB*(SolidAngle_ColdBaffle-SolidAngle_Filter).*T_Filter_Fast'+BB*SolidAngle_Filter.*T_Filter';
% for i=1:size(lambda,2)
% Background(i)=sum(R(1:i),2)*10^-9;
% end
% semilogy(lambda*10^6,Background,'-.b');
% hold on;

vline=line([param_H2RG_A_CutOff*10^6 param_H2RG_A_CutOff*10^6],ylim,'LineWidth',4,'Color','m');
vline=line([param_H2RG_SNAP_CutOff*10^6 param_H2RG_SNAP_CutOff*10^6],ylim,'LineWidth',4,'Color','b');
ylabel('Number of photons');
hold on;
yyaxis right;
ylabel('Transmittance (%)');
semilogy(lambda*10^6,T_Filter,'black');
ylim([0.001 1]);
hold on;
legend( 'CASE A: Filter F/0.5 and Cold Baffle F/1.2',...
    'CASE B: Filter F/0.5 and Cold Baffle F/8',...
         'H2RG 2.5um Cutoff','H2RG 1.7um Cutoff',...
        'Filter Transmission');



title('Integrated photons over one pixel during one second VS cuting wavelength');
xlabel('wavelength (um)'); 


