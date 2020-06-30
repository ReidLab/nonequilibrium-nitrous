clear

%first-type BC;
%kinetic model;
%see Van Genuchten pg. 16 and 28:



%data = csvread('ObjC_1_gas.csv',1,0);
data = csvread('Obj2b_part1_gas.csv',1,0);
data_time = data(:,1).*60; %time in sec
ethane = data(:,2);
helium = data(:,3);
n2o = data(:,4);
sf6 = data(:,5);
ch4 = data(:,6);



data_br = csvread('Obj2b_part1_br.csv',1,0);
%data_br = csvread('Obj1b_br.csv',1,0);
Br_time = data_br(:,1).*60;
Br = data_br(:,2);



%If Vg/Vw = .01, then
%R_He = 2
%R_SF6 = 2.5
%R_ethane = 1.2
%R_N2O = 1.015

%time units are in seconds;
%length units are in m;

C_o = 1; %dimensionless
t = 0:100:60000; %0:100:200000; %default delta t is 25 %secondary t is 50
D = 1.025e-06; 


%2.69e-6, dispersion coefficients, m2/s, default 5e-7; determined from bromide only.
    %%%%%%3.1e-6 for bromide
v_m = 3.625e-5;  

%velocity, 7.09e-5,  m/s, default 1e-4 m/s
    %%%%%%%3e-5 for bromide
    %try 1.19e-4 to 1.43e-4, expected values based on reactor geometry +
    %porosity
L = .3;
theta_m = 0.586; %mobile phase volume (dynamic, macropores)
theta_im = 0.31;  %immobile phase volume (stagnant, micropores)
phi_m = theta_m/(theta_m + theta_im);  %fraction of liquid phase that is mobile

P = v_m*L/D; %Peclet number
PV = v_m.*t.*phi_m./L;     %pore volumes, tau in equations in Table 2 in Van Genuchten
data_PV_br = v_m.*Br_time.*phi_m./L;

data_PV = v_m.*data_time.*phi_m./L;

data_PV_br_plot = v_m.*Br_time.*phi_m./L;


figure(11)
clf
h1 = plot(data_PV_br_plot, Br, '*k','MarkerFaceColor', 'k', 'MarkerSize', 18)
hold on
h2 = plot(data_PV, n2o, 'dm','MarkerFaceColor', 'm', 'MarkerSize', 18)
hold on
h3 = plot(data_PV, sf6, 'vg','MarkerFaceColor', 'g', 'MarkerSize', 18)
hold on
h4 = plot(data_PV, ethane, 'sb','MarkerFaceColor', 'b', 'MarkerSize', 18)
hold on
h5 = plot(data_PV, ch4, 'oc','MarkerFaceColor', 'c', 'MarkerSize', 18)
hold on
h6 = plot(data_PV, helium, '^r','MarkerFaceColor', 'r', 'MarkerSize', 18)
hold on






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






measured = Br; %define which solute you are modeling




R = 1.000000001; %ethane 2.7
%R = [50:2.5:70];
%R_m = [1.01:0.05:1.41];
R_m = 1.000000001;
alpha = 1.0e-4;  %ethane 1.3e-4
%alpha = [1e-7:1e-7:20e-7];

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

%RSS = sum((C(I) - measured').^2)


loopcnt = loopcnt + 1

     results_matrix(loopcnt,1) = loopcnt;
     results_matrix(loopcnt,2) = R(j);
     results_matrix(loopcnt,3) = R_m(l);
     results_matrix(loopcnt,4) = alpha(k);
     %results_matrix(loopcnt,5) = RSS;

    end
    
    end

end






figure(11)
plot(PV, C,'-k', 'LineWidth',1.5)
ylim([-0.1 1.1])
xlim([0 2])
ylabel('C/C_o')
xlabel('PV')
hold on

%column 1 of results_matrix is loopcnt
%column 2 of results_matrix is R
%column 3 of results_matrix is R_m
%column 4 of results_matrix is alpha
%column 5 of results_matrix is sum of squared residuals

% y = [results_matrix(:,1)'; results_matrix(:,2)'; results_matrix(:,3)'; results_matrix(:,4)'; results_matrix(:,5)']
% fileID = fopen('results_matrix.txt','w');
% %fprintf(fileID, 'Exponential Function\n\n');
% fprintf(fileID,'%0.3f %0.3e %0.3e %0.3e %0.4f\n',y);
% fclose(fileID);




% y = [PV; C];
% fileID = fopen('exptable.txt','w');
% %fprintf(fileID, 'Exponential Function\n\n');
% fprintf(fileID,'%f %f\n',y);
% fclose(fileID);



%find indices in PV that are closest PV to data: 
%edges = [-Inf, mean([PV(2:end); PV(1:end-1)]), +Inf];
%I = discretize(data_PV, edges);

%edges = [-Inf, mean([Br_time(2:end); Br_time(1:end-1)]), +Inf];
%I = discretize(Br_time, edges);



% %find indices in PV that are closest PV to data: 
% edges = [-Inf, mean([PV(2:end); PV(1:end-1)]), +Inf];
% %I = discretize(data_PV_br, edges);
% I = discretize(data_PV, edges);
% 
% RSS = sum((C(I) - measured').^2)
% TSS = sum((measured - mean(measured)).^2)
% R2 = 1 - (RSS/TSS)


measured = sf6; %define which solute you are modeling




R = 57.5; %ethane 2.7
%R = [50:2.5:70];
%R_m = [1.01:0.05:1.41];
R_m = 1.01;
alpha = 1.4e-6;  %ethane 1.3e-4
%alpha = [1e-7:1e-7:20e-7];

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






figure(11)
plot(PV, C,'-k', 'LineWidth',1.5)
ylim([-0.1 1.1])
xlim([0 2])
ylabel('C/C_o')
xlabel('PV')
hold on 





measured = helium; %define which solute you are modeling




R = 39; %ethane 2.7
%R = [50:2.5:70];
%R_m = [1.01:0.05:1.41];
R_m = 6.1;
alpha = 1.05e-6;  %ethane 1.3e-4
%alpha = [1e-7:1e-7:20e-7];

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






figure(11)
plot(PV, C,'-r', 'LineWidth',1.5)
ylim([-0.1 1.1])
xlim([0 2])
ylabel('C/C_o')
xlabel('PV')
hold on

%column 1 of results_matrix is loopcnt
%column 2 of results_matrix is R
%column 3 of results_matrix is R_m
%column 4 of results_matrix is alpha
%column 5 of results_matrix is sum of squared residuals

% y = [results_matrix(:,1)'; results_matrix(:,2)'; results_matrix(:,3)'; results_matrix(:,4)'; results_matrix(:,5)']
% fileID = fopen('results_matrix.txt','w');
% %fprintf(fileID, 'Exponential Function\n\n');
% fprintf(fileID,'%0.3f %0.3e %0.3e %0.3e %0.4f\n',y);
% fclose(fileID);




% y = [PV; C];
% fileID = fopen('exptable.txt','w');
% %fprintf(fileID, 'Exponential Function\n\n');
% fprintf(fileID,'%f %f\n',y);
% fclose(fileID);



%find indices in PV that are closest PV to data: 
%edges = [-Inf, mean([PV(2:end); PV(1:end-1)]), +Inf];
%I = discretize(data_PV, edges);

%edges = [-Inf, mean([Br_time(2:end); Br_time(1:end-1)]), +Inf];
%I = discretize(Br_time, edges);



% %find indices in PV that are closest PV to data: 
% edges = [-Inf, mean([PV(2:end); PV(1:end-1)]), +Inf];
% %I = discretize(data_PV_br, edges);
% I = discretize(data_PV, edges);
% 
% RSS = sum((C(I) - measured').^2)
% TSS = sum((measured - mean(measured)).^2)
% R2 = 1 - (RSS/TSS)


measured = sf6; %define which solute you are modeling




R = 57.5; %ethane 2.7
%R = [50:2.5:70];
%R_m = [1.01:0.05:1.41];
R_m = 1.01;
alpha = 1.4e-6;  %ethane 1.3e-4
%alpha = [1e-7:1e-7:20e-7];

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






figure(11)
plot(PV, C,'-g', 'LineWidth',1.5)
ylim([-0.1 1.1])
xlim([0 2])
ylabel('C/C_o')
xlabel('PV')
hold on 





measured = ethane; %define which solute you are modeling




R = 8.5; %ethane 2.7
%R = [50:2.5:70];
%R_m = [1.01:0.05:1.41];
R_m = 1.01;
alpha = 1.2e-5;  %ethane 1.3e-4
%alpha = [1e-7:1e-7:20e-7];

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






figure(11)
plot(PV, C,'-b', 'LineWidth',1.5)
ylim([-0.1 1.1])
xlim([0 3.5])
ylabel('C/C_o')
xlabel('PV')
hold on

measured = n2o; %define which solute you are modeling




R = 1.6; %ethane 2.7
%R = [50:2.5:70];
%R_m = [1.01:0.05:1.41];
R_m = 1.2;
alpha = 3e-5;  %ethane 1.3e-4
%alpha = [1e-7:1e-7:20e-7];

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






figure(11)
plot(PV, C,'--m', 'LineWidth',1.5)
ylim([-0.1 1.1])
xlim([0 3.5])
ylabel('C/C_o')
xlabel('PV')






measured = ch4; %define which solute you are modeling




R = 6.9; %ethane 2.7
%R = [50:2.5:70];
%R_m = [1.01:0.05:1.41];
R_m = 1.15;
alpha = 2.25e-5;  %ethane 1.3e-4
%alpha = [1e-7:1e-7:20e-7];

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






figure(11)
plot(PV, C,'--c', 'LineWidth',1.5)
ylim([-0.1 1.1])
xlim([0 3.5])
ylabel('C/C_o')
xlabel('PV')



figure(111)
plot(PV, C,'-r')
ylim([-0.1 1.2])
xlim([0 4])
ylabel('C/C_o')
xlabel('PV')
%text(.4, .8, sprintf('R = %0.2f', R))
%text(.4, .9, sprintf('alpha = %0.2e s^{-1}', alpha))
%text(.4, .7, sprintf('R^2 = %0.3f', R2))
%saveas(gcf,'output.eps','epsc')
figure(11)
legend('Br^-','N_2O', 'SF_6', 'C_2H_6', 'CH_4', 'He', 'Location', 'northwest'); %
title('Biotic, Transient, Woodchip')
set(gca, 'FontSize', 34);








