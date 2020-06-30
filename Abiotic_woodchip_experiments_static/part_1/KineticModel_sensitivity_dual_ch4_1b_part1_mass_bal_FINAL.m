clear all
close all

%first-type BC;
%kinetic model;
%see Van Genuchten pg. 16 and 28:



%data = csvread('ObjC_1_gas.csv',1,0);
data = csvread('Obj1b_part1_gas.csv',1,0);
data_time = data(:,1).*60; %time in sec
ethane = data(:,2);
helium = data(:,3);
n2o = data(:,4);
sf6 = data(:,5);
ch4 = data(:,6);



data_br = csvread('Obj1b_part1_br.csv',1,0);
%data_br = csvread('Obj1b_br.csv',1,0);
Br_time = data_br(:,1).*60;
Br = data_br(:,2);

measured = Br; %define which solute you are modeling

%If Vg/Vw = .01, then
%R_He = 2
%R_SF6 = 2.5
%R_ethane = 1.2
%R_N2O = 1.015

%time units are in seconds;
%length units are in m;

C_o = 1; %dimensionless
t = 0:50:49200;   %default delta t is 25
%D = [2.125e-6:0.0625e-6:2.375e-6]; 
D = 2.75e-6;

%2.69e-6, dispersion coefficients, m2/s, default 5e-7; determined from bromide only.
    %%%%%%3.1e-6 for bromide
%v_m = 7.75e-5:0.0625e-5:8e-5;  
v_m = 7.125e-5;


%velocity, 7.09e-5,  m/s, default 1e-4 m/s
    %%%%%%%3e-5 for bromide
    %try 1.19e-4 to 1.43e-4, expected values based on reactor geometry +
    %porosity
L = .3;
theta_m = 0.586; %mobile phase volume (dynamic, macropores)
theta_im = 0.31; %immobile phase volume (stagnant, micropores)
phi_m = theta_m/(theta_m + theta_im)  %fraction of liquid phase that is mobile


data_PV = v_m(1).*data_time.*phi_m./L

data_PV_br_plot = v_m(1).*Br_time.*phi_m./L


figure(11)
clf
h1 = plot(data_PV_br_plot, Br,'ok')
hold on
h2 = plot(data_PV, ethane, 'sb')
hold on
h3 = plot(data_PV, helium, '^r')
hold on
h4 = plot(data_PV, sf6, 'vg')
hold on
h5 = plot(data_PV, n2o, 'dm')
hold on
h6 = plot(data_PV, ch4, '+c')
hold on


R = 1.0001; %ethane 2.7
%R = [2.0:.5:6.0];
R_m = 1.0001;
%R_m = [1.0:.5:3.0];
alpha = 0.0001;  %ethane 1.3e-4
%alpha = [.9e-4:.1e-4:2.5e-4];
loopcnt = 0;

results_matrix = zeros(length(v_m)*length(D),4);
%column 1 of results_matrix is loopcnt
%column 2 of results_matrix is R
%column 3 of results_matrix is alpha
%column 4 of results_matrix is sum of squared residuals



beta = phi_m*R_m/R;


for j = 1:length(v_m)
    
    for k = 1:length(D)
        
        
w = alpha*L*(R-1)/(v_m(j)*theta_m);
P = v_m(j)*L/D(k); %Peclet number
PV = v_m(j).*t.*phi_m./L;     %pore volumes, tau in equations in Table 2 in Van Genuchten
data_PV_br = v_m(j).*Br_time.*phi_m./L;


%G term in Van Genuchten, equation SI-1 in Table 2:
G_T = .5 *  erfc( sqrt(P./(4.*beta.*R.*PV)).*(beta*R-PV)) + .5*exp(P).*erfc( sqrt (P./(4.*beta.*R.*PV) ).*(beta*R + PV) );
%IMPORTANT NOTE: Equation SI-1 in Table 2 has the final term as beta*R-tau.
%I think this is an error and that it should be beta*R+tau.  At the top of page 18 they note: 
%"Expressions for G(T) in Table 2 follow from those given in Table 1 by replacing T by
%tau and R by (beta*R).  In Table 1, the final term is R+T, so based on
%their statement the final term should be B*R+tau.
%This is how it is coded above, and this equation gives an effluent C
%ranging from 0 to 1, while the original equation from Van Genuchten has an
%effluent C ranging from 0 to a number much greater than 1.

%exponential term multiplied by G
term = exp(-w.*PV/(beta*R));


integrand_vec = zeros(1,length(PV));

delta_tau = 0.001;   %tau time increment, also d_tau dummy variable
                    %default .001

for i = 2:length(PV);  %2:length(PV)
    
    tau = [0:delta_tau:PV(i)];
    
G_tau = .5 *  erfc( sqrt(P./(4.*beta.*R.*tau)).*(beta*R-tau)) + .5*exp(P).*erfc( sqrt (P./(4.*beta.*R.*tau) ).*(beta*R + tau) );
    
    a = w.*tau;
    b = w.*(PV(i) - tau)./( (1 - beta)*R );
    
    squig = 2*sqrt(a.*b);
    
%In equation 36, the I_0 and I_1 are modified bessel functions, see
%notation on pg. 50-51
I_0 = besseli(0,squig);  %modified bessel function of the first kind, order 0
I_1 = besseli(1,squig);  %modified bessel function of the first kind, order 1
    
%H term defined by Equation 36:
H = exp(-a-b).* ( I_0/beta + (squig.*I_1)./(2.*b.*(1-beta)) );


prod = G_tau.*H;

integrand = delta_tau*trapz(prod(1:length(prod)-1));   %last element of product is NaN, since the last element of H is NaN ...
                                                    %due to division by 0
                                                    %in b term
integrand_vec(i) = integrand;


end


BrC = G_T.*term + w/R .* integrand_vec;


%find indices in PV that are closest PV to data: 
edges = [-Inf, mean([PV(2:end); PV(1:end-1)]), +Inf];
%I = discretize(data_PV_br, edges);
I = discretize(data_PV_br, edges);

RSS = sum((BrC(I) - measured').^2)


loopcnt = loopcnt + 1

     results_matrix(loopcnt,1) = loopcnt;
     results_matrix(loopcnt,2) = v_m(j);
     results_matrix(loopcnt,3) = D(k);
     results_matrix(loopcnt,4) = RSS;

    end
    
end




figure(111)
clf
plot(PV, BrC,'-')
ylim([-0.1 1.1])
xlim([0 2])
ylabel('C/C_o')
xlabel('PV')

%column 1 of results_matrix is loopcnt
%column 2 of results_matrix is R
%column 3 of results_matrix is alpha
%column 4 of results_matrix is sum of squared residuals

y = [results_matrix(:,1)'; results_matrix(:,2)'; results_matrix(:,3)'; results_matrix(:,4)']
fileID = fopen('results_matrix.txt','w');
%fprintf(fileID, 'Exponential Function\n\n');
fprintf(fileID,'%0.3f %0.3e %0.3e %0.4f\n',y);
fclose(fileID);




y = [PV; BrC];
fileID = fopen('exptable.txt','w');
%fprintf(fileID, 'Exponential Function\n\n');
fprintf(fileID,'%f %f\n',y);
fclose(fileID);


%find indices in PV that are closest PV to data: 
edges = [-Inf, mean([PV(2:end); PV(1:end-1)]), +Inf];
%I = discretize(data_PV_br, edges);
I = discretize(data_PV_br, edges);

RSS = sum((BrC(I) - measured').^2)
TSS = sum((measured - mean(measured)).^2)
R2 = 1 - (RSS/TSS)



figure(11)
plot(PV, BrC,'-r')
ylim([-0.1 1.1])
xlim([0 9])
ylabel('C/C_o')
xlabel('PV')
text(.4, .8, sprintf('R = %0.2f', R))
text(.4, .9, sprintf('alpha = %0.2e s^{-1}', alpha))
text(.4, .7, sprintf('R^2 = %0.3f', R2))
saveas(gcf,'output.eps','epsc')


figure (101)
brarea = cumtrapz(PV,BrC);
totbrarea = trapz(PV,BrC);
totbrmass = data_PV(end)*1;
normalizedbr = brarea/totbrmass;
plot(PV, normalizedbr)



measured = ch4; %define which solute you are modeling

%fit parameters:
%R = 2;      %MUST BE > 1
%alpha = 1e-3;    %rate constant, try starting in the range 5e-4 to 5e-3

                    %1/s, 100 d^-1 from Geistlinger (2005) = 1e-3 1/s
                    %vulava has 5e-4 - 1e-3  1/s for alpha for SF6, so
                    %probably faster for Helium, N2O
                    %velocity in bioreactor on the order of .003 m/s, much slower than in Vulava, so 
                    %expect lower alpha
%for determination of R-squared:



%dimensionless parameters
%defined in Van Genuchten pg. 15, with the exception of w, beta, and R
%For w and beta, defined on pg. 28 for one-site kinetic
%non-equilibrium adsorption.  For R, comes from Fry et al. (1994) and Vulava et al. (2002): 
%R = 1 + KH*(Vg/Vw)   %use this with specific K_H values to determine Vg/Vw
%based on fit R parameter




R = 3.9; %ethane 2.7
%R = [1.1:0.02:1.2];
%R_m = [1.1:0.1:2.8];
R_m = 1.3;
alpha = 5.5e-6;  %ethane 1.3e-4
%alpha = [1e-5 1e-6 1e-7 1e-8 1e-9];%8e-6:4e-6:8e-5];

loopcnt = 0;

results_matrix = zeros(length(R)*length(alpha),5);
%column 1 of results_matrix is loopcnt
%column 2 of results_matrix is R
%column 3 of results_matrix is R_m
%column 4 of results_matrix is alpha
%column 5 of results_matrix is sum of squared residuals




for j = 1:length(R)
    
    for k = 1:length(alpha)
        
        for l = 1:length(R_m)
        
        
w = alpha(k)*L*(R(j)-1)/(v_m*theta_m);
beta = phi_m*R_m(l)/R(j);



%G term in Van Genuchten, equation SI-1 in Table 2:
G_T = .5 *  erfc( sqrt(P./(4.*beta.*R(j).*PV)).*(beta*R(j)-PV)) + .5*exp(P).*erfc( sqrt (P./(4.*beta.*R(j).*PV) ).*(beta*R(j) + PV) );
%IMPORTANT NOTE: Equation SI-1 in Table 2 has the final term as beta*R-tau.
%I think this is an error and that it should be beta*R+tau.  At the top of page 18 they note: 
%"Expressions for G(T) in Table 2 follow from those given in Table 1 by replacing T by
%tau and R by (beta*R).  In Table 1, the final term is R+T, so based on
%their statement the final term should be B*R+tau.
%This is how it is coded above, and this equation gives an effluent C
%ranging from 0 to 1, while the original equation from Van Genuchten has an
%effluent C ranging from 0 to a number much greater than 1.

%exponential term multiplied by G
term = exp(-w.*PV/(beta*R(j)));

%some figures to visualize how each of these terms are affected by
%different fit parameters:
%use these to troubleshoot problems if the model acts weird with certain
%paremeters.
%figure(5)
%plot(PV,term)
%title('Term')

integrand_vec = zeros(1,length(PV));

delta_tau = 0.001;   %tau time increment, also d_tau dummy variable
                    %default .001

for i = 2:length(PV);  %2:length(PV)
    
    tau = [0:delta_tau:PV(i)];
    
G_tau = .5 *  erfc( sqrt(P./(4.*beta.*R(j).*tau)).*(beta*R(j)-tau)) + .5*exp(P).*erfc( sqrt (P./(4.*beta.*R(j).*tau) ).*(beta*R(j) + tau) );
    
    a = w.*tau;
    b = w.*(PV(i) - tau)./( (1 - beta)*R(j) );
    
    squig = 2*sqrt(a.*b);
    
%In equation 36, the I_0 and I_1 are modified bessel functions, see
%notation on pg. 50-51
I_0 = besseli(0,squig);  %modified bessel function of the first kind, order 0
I_1 = besseli(1,squig);  %modified bessel function of the first kind, order 1
    
%H term defined by Equation 36:
H = exp(-a-b).* ( I_0/beta + (squig.*I_1)./(2.*b.*(1-beta)) );


prod = G_tau.*H;

integrand = delta_tau*trapz(prod(1:length(prod)-1));   %last element of product is NaN, since the last element of H is NaN ...
                                                    %due to division by 0
                                                    %in b term
integrand_vec(i) = integrand;


end


%first_term = G_T.*term;
%figure(14)
%plot(PV,first_term)
%title('G*exp term')

%C = G(1:length(G)-1).*term(1:length(term)-1) + (w/R).*ConvInt;   %equation 35 in van genuchten
C = G_T.*term + w/R(j) .* integrand_vec;


%find indices in PV that are closest PV to data: 
edges = [-Inf, mean([PV(2:end); PV(1:end-1)]), +Inf];
%I = discretize(data_PV_br, edges);
I = discretize(data_PV, edges);

RSS = sum((C(I) - measured').^2)


loopcnt = loopcnt + 1

     results_matrix(loopcnt,1) = loopcnt;
     results_matrix(loopcnt,2) = R(j);
     results_matrix(loopcnt,3) = R_m(l);
     results_matrix(loopcnt,4) = alpha(k);
     results_matrix(loopcnt,5) = RSS;

    end
    
    end

end

C(C > 1) = 1;

figure(11)
plot(PV, C,'-m')
hold on

%Effluent Mass Balance
figure (101)
hold on
ch4area = cumtrapz(PV,C);
totch4area = trapz(PV,C);
totch4mass = (data_PV(end))*1;
normalizedch4 = ch4area/totch4mass;
plot(PV, normalizedch4)
% massbaln2o = normalizedn2o./normalizedbr;
% plot(PV, massbaln2o)


%Headspace Mass Balance
tracersoldata = csvread('Obj1b_part1_tracer_sol.csv',1,0);
ch4tracer = tracersoldata(:,6);
tracertime = tracersoldata(:,1).*60;
tracerPV = v_m.*tracertime.*phi_m./L;
interptracer = interp1(tracerPV, ch4tracer, PV);
aveinterptracer = mean([interptracer(1:end-1);interptracer(2:end)]);
avetracermass = PV(2).*aveinterptracer;
avetracermass(isnan(avetracermass))=0;
tottracermass = sum(avetracermass);
headspacedata = csvread('Obj1b_part1_headspace_ch4.csv',1,0);
seg1time = headspacedata(:,1).*60;
seg1data = headspacedata(:,2);
seg1PV = v_m.*seg1time.*phi_m./L;
seg2time = headspacedata(:,3).*60;
seg2data = headspacedata(:,4);
seg2PV = v_m.*seg2time.*phi_m./L;
seg3time = headspacedata(:,5).*60;
seg3data = headspacedata(:,6);
seg3PV = v_m.*seg3time.*phi_m./L;
seg1interp = interp1(seg1PV, seg1data, PV);
seg2interp = interp1(seg2PV, seg2data, PV);
seg3interp = interp1(seg3PV, seg3data, PV);
headspaceinterp = seg1interp+seg2interp+seg3interp;
normalizedheadspace = headspaceinterp/tottracermass;
combinedheadspace = csvread('Obj1b_part1_headspace_ch4_combined.csv',1,0);
combinedtime = combinedheadspace(:,1).*60;
combineddata = combinedheadspace(:,2);
combinedPV = v_m.*combinedtime.*phi_m./L;
normalizedcombined = combineddata/tottracermass;

%Bubble Mass Balance
bubblech4 = normalizedbr-(normalizedch4);%+normalizedheadspace);
figure(22)
plot(PV,bubblech4)
bubblech4change = diff(bubblech4);
find(bubblech4change == max(bubblech4change(:)))


%truncate data to 6.64 pore volumes
shortnormalizedch4 = normalizedch4(1:856);
shortbubblech4 = bubblech4(1:856);
shortnormalizedheadspace = normalizedheadspace(1:856);
shortPV = PV(1:856);

figure(222)
allch4 = [shortnormalizedch4' shortnormalizedheadspace' shortbubblech4'];
massbalpv = area(shortPV,allch4)
massbalpv(1).FaceColor = [0/255 128/255 255/255];
massbalpv(2).FaceColor = [252/255 226/255 5/255];
massbalpv(3).FaceColor = [255/255 40/255 0/255];
ylabel('M/M_T', 'fontsize', 25);
xlabel('Pore Volumes', 'fontsize', 25);
ylim([0 1]);
xlim([0 8])
%title('CH_4 Distribution Throughout Reactor', 'fontsize', 36);
%legend('Effluent', 'Headspace', 'Bubble Entrapped', 'Location', 'northwest');
set(gca, 'FontSize', 20);

figure(333)
finaleffluent = normalizedch4(end);
finalbubble = 0;
finalheadspace = normalizedcombined(end);
remainingheadspace = shortnormalizedheadspace(end);
liberatedbubble = finalheadspace-shortnormalizedheadspace(end);
predbubble = shortbubblech4(end);
if liberatedbubble-predbubble<0;
    unantbubble = 0;
else
    unantbubble = liberatedbubble - predbubble;
end
barallch4 = [finaleffluent, remainingheadspace, predbubble, unantbubble; nan(1,4)];
massbalbar = bar(barallch4, 'stacked')
ylim([0 1]);
set(massbalbar,{'FaceColor'},{[0/255 128/255 255/255];[252/255 226/255 5/255];[255/255 40/255 0/255];'k';});
xlim=get(gca,'xlim');
hold on
%For making legend with dashed line
%plot(xlim,[0.5 0.5], '--k', 'LineWidth', 2)
%legend('Effluent', 'Evolved Headspace', 'Predicted Bubble Entrapment/Liberation', 'Unanticipated Bubble Liberation' , 'Reactor Drainage');
set(gca, 'FontSize', 20);

recovery = finaleffluent + remainingheadspace + predbubble + unantbubble


%Extend data to 7.64 PVs (Run for 6.64 PVs and add 1 PV for "effluent"
%remaining in reactor at experiment end
integralPV = PV(1:985);
integralBrC = BrC(1:985);
integraleff = C(1:985);
Br_area = trapz(integralPV, integralBrC);
Eff_area = trapz(integralPV, integraleff);
bubblearea = Br_area-Eff_area;
btestimate_effluent = Eff_area/Br_area;
btestimate_retainedinreactor = bubblearea/Br_area;