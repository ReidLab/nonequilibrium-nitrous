clear



%first-type BC;
%kinetic model;
%see Van Genuchten pg. 16 and 28:

%data = csvread('Obj1b_part1_data.csv',1,0);
gas_min = xlsread('GB1_gas_data.csv','A2:A21'); %time in min
data_time = gas_min.*60;
ethane = xlsread('GB1_gas_data.csv','B2:B21');
helium = xlsread('GB1_gas_data.csv','C2:C21');
n2o = xlsread('GB1_gas_data.csv','D2:D21');
sf6 = xlsread('GB1_gas_data.csv','E2:E21');
ch4 = xlsread('GB1_gas_data.csv','F2:F21');
br_min = xlsread('GB1_br_data.csv','A2:A28');
br_time = br_min.*60;
Br = xlsread('GB1_br_data.csv','B2:B28');


%time units are in seconds;
%length units are in m;

C_o = 1; %dimensionless
t = 0:50:40000;
D = 1.8e-7;    %dispersion coefficients, m2/s, default 5e-7; determined from bromide only.
    %%%%%%1.3e-6 for bromide
v = 7.45e-5;   %velocity, m/s, default 1e-4 m/s
    %%%%%%%4.65e-5 for bromide
    %try 1.19e-4 to 1.43e-4, expected values based on reactor geometry +
    %porosity
L = .3048;

data_PV = v.*data_time./L;
brdata_PV = v.*br_time./L;

figure(11)
clf
h1 = plot(brdata_PV, Br, '*k', 'MarkerSize', 18);
hold on
h2 = plot(data_PV, n2o, 'dm', 'MarkerFaceColor', 'm', 'MarkerSize', 18);
hold on
h3 = plot(data_PV, sf6, 'vg','MarkerFaceColor', 'g', 'MarkerSize', 18);
hold on
h4 = plot(data_PV, ethane, 'sb','MarkerFaceColor', 'b', 'MarkerSize', 18);
hold on
h5 = plot(data_PV, ch4, 'oc','MarkerFaceColor', 'c', 'MarkerSize', 18);
hold on
h6 = plot(data_PV, helium, '^r','MarkerFaceColor', 'r', 'MarkerSize', 18);
hold on

%%%Bromide
measured = Br;

D = 1.8e-7;%1.6e-7:0.025e-7:2e-7;    %dispersion coefficients, m2/s, default 5e-7; determined from bromide only.
    %%%%%%1.3e-6 for bromide
v = 7.45e-5;%7.4e-5:0.025e-5:7.6e-5;   %velocity, m/s, default 1e-4 m/s

R = 1.00000001;%ethR0;
%heR = 5;        %MUST BE > 1
alpha = .0001;    %rate constant, try starting in the range 5e-4 to 5e-3

                    %1/s, 100 d^-1 from Geistlinger (2005) = 1e-3 1/s
                    %vulava has 5e-4 - 1e-3  1/s for alpha for SF6, so
                    %probably faster for Helium, N2O
                    %velocity in bioreactor on the order of .003 m/s, much slower than in Vulava, so 
                    %expect lower alpha

%dimensionless parameters
%defined in Van Genuchten pg. 15, with the exception of w, beta, and R
%For w and beta, defined on pg. 28 for one-site kinetic
%non-equilibrium adsorption.  For R, comes from Fry et al. (1994) and Vulava et al. (2002): 
%R = 1 + KH*(Vg/Vw)   %use this with specific K_H values to determine Vg/Vw
%based on fit R parameter

loopcnt = 0;

results_matrix = zeros(length(D)*length(v),4);



for j = 1:length(v)
    
    for k = 1:length(D)

        
P = v(j)*L/D(k); %Peclet number
PV = v(j).*t./L;     %pore volumes, tau in equations in Table 2 in Van Genuchten
data_PV_br = v(j).*br_time./L;
w = alpha*L*(R-1)/v(j);
beta = 1/R;

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

%some figures to visualize how each of these terms are affected by
%different fit parameters:
%use these to troubleshoot problems if the model acts weird with certain
%paremeters.
%figure(5)
%plot(PV,term)
%title('Term')

integrand_vec = zeros(1,length(PV));

delta_tau = 0.001;   %tau time increment, also d_tau dummy variable

for i = 2:length(PV)  %2:length(PV)
    
    tau = (0:delta_tau:PV(i));
    
G_tau = .5 *  erfc( sqrt(P./(4.*beta.*R.*tau)).*(beta*R-tau)) + .5*exp(P).*erfc( sqrt (P./(4.*beta.*R.*tau) ).*(beta*R + tau) );
    
    a = w.*tau;
    b = w.*(PV(i) - tau)./( (1 - beta)*R );
    
    squig = 2*sqrt(a.*b);
    
%In equation 36, the I_0 and I_1 are modified bessel functions, see
%notation on pg. 50-51
I_0 = besseli(0,squig);  %modified bessel function of the first kind, order 0
I_1 = besseli(1,squig);  %modified bessel function of the first kind, order 1
    
%H term defined by Equation 36:
H = exp(-a-b).* (I_0/beta + (squig.*I_1)./(2.*b.*(1-beta)) );


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
C = G_T.*term + w/R .* integrand_vec;


%find indices in PV that are closest PV to data: 
edges = [-Inf, mean([PV(2:end); PV(1:end-1)]), +Inf];
%I = discretize(data_PV_br, edges);
I = discretize(data_PV_br, edges);

RSS = sum((C(I) - measured').^2)


loopcnt = loopcnt + 1

     results_matrix(loopcnt,1) = loopcnt;
     results_matrix(loopcnt,2) = v(j);
     results_matrix(loopcnt,4) = D(k);
     results_matrix(loopcnt,5) = RSS;
     
     
    end
end

plot(PV, C, '-k', 'LineWidth',1.5)
hold on



%%%N2O
measured = n2o;

R = 1.09;%1.1:0.05:1.7;%1.7;%ethR0;
%heR = 5;        %MUST BE > 1
alpha = 2e-4;%10e-5:3e-5:40e-5; %rate constant, try starting in the range 5e-4 to 5e-3

                    %1/s, 100 d^-1 from Geistlinger (2005) = 1e-3 1/s
                    %vulava has 5e-4 - 1e-3  1/s for alpha for SF6, so
                    %probably faster for Helium, N2O
                    %velocity in bioreactor on the order of .003 m/s, much slower than in Vulava, so 
                    %expect lower alpha

%dimensionless parameters
%defined in Van Genuchten pg. 15, with the exception of w, beta, and R
%For w and beta, defined on pg. 28 for one-site kinetic
%non-equilibrium adsorption.  For R, comes from Fry et al. (1994) and Vulava et al. (2002): 
%R = 1 + KH*(Vg/Vw)   %use this with specific K_H values to determine Vg/Vw
%based on fit R parameter

loopcnt = 0;

results_matrix = zeros(length(R)*length(alpha),4);



for j = 1:length(R)
    
    for k = 1:length(alpha)

        
P = v*L/D; %Peclet number
PV = v.*t./L;     %pore volumes, tau in equations in Table 2 in Van Genuchten
w = alpha(k)*L*(R(j)-1)/v;
beta = 1/R(j);

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

for i = 2:length(PV)  %2:length(PV)
    
    tau = (0:delta_tau:PV(i));
    
G_tau = .5 *  erfc( sqrt(P./(4.*beta.*R(j).*tau)).*(beta*R(j)-tau)) + .5*exp(P).*erfc( sqrt (P./(4.*beta.*R(j).*tau) ).*(beta*R(j) + tau) );
    
    a = w.*tau;
    b = w.*(PV(i) - tau)./( (1 - beta)*R(j) );
    
    squig = 2*sqrt(a.*b);
    
%In equation 36, the I_0 and I_1 are modified bessel functions, see
%notation on pg. 50-51
I_0 = besseli(0,squig);  %modified bessel function of the first kind, order 0
I_1 = besseli(1,squig);  %modified bessel function of the first kind, order 1
    
%H term defined by Equation 36:
H = exp(-a-b).* (I_0/beta + (squig.*I_1)./(2.*b.*(1-beta)) );


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
     results_matrix(loopcnt,4) = alpha(k);
     results_matrix(loopcnt,5) = RSS;
     
     
    end
end

plot(PV, C,'-m', 'LineWidth',1.5)
hold on






%%%SF6
measured = sf6;

R = 2.1;%1.1:0.05:1.7;%1.7;%ethR0;
%heR = 5;        %MUST BE > 1
alpha = 1.2e-4;%10e-5:3e-5:40e-5; %rate constant, try starting in the range 5e-4 to 5e-3

                    %1/s, 100 d^-1 from Geistlinger (2005) = 1e-3 1/s
                    %vulava has 5e-4 - 1e-3  1/s for alpha for SF6, so
                    %probably faster for Helium, N2O
                    %velocity in bioreactor on the order of .003 m/s, much slower than in Vulava, so 
                    %expect lower alpha

%dimensionless parameters
%defined in Van Genuchten pg. 15, with the exception of w, beta, and R
%For w and beta, defined on pg. 28 for one-site kinetic
%non-equilibrium adsorption.  For R, comes from Fry et al. (1994) and Vulava et al. (2002): 
%R = 1 + KH*(Vg/Vw)   %use this with specific K_H values to determine Vg/Vw
%based on fit R parameter

loopcnt = 0;

results_matrix = zeros(length(R)*length(alpha),4);



for j = 1:length(R)
    
    for k = 1:length(alpha)

        
P = v*L/D; %Peclet number
PV = v.*t./L;     %pore volumes, tau in equations in Table 2 in Van Genuchten
w = alpha(k)*L*(R(j)-1)/v;
beta = 1/R(j);

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

for i = 2:length(PV)  %2:length(PV)
    
    tau = (0:delta_tau:PV(i));
    
G_tau = .5 *  erfc( sqrt(P./(4.*beta.*R(j).*tau)).*(beta*R(j)-tau)) + .5*exp(P).*erfc( sqrt (P./(4.*beta.*R(j).*tau) ).*(beta*R(j) + tau) );
    
    a = w.*tau;
    b = w.*(PV(i) - tau)./( (1 - beta)*R(j) );
    
    squig = 2*sqrt(a.*b);
    
%In equation 36, the I_0 and I_1 are modified bessel functions, see
%notation on pg. 50-51
I_0 = besseli(0,squig);  %modified bessel function of the first kind, order 0
I_1 = besseli(1,squig);  %modified bessel function of the first kind, order 1
    
%H term defined by Equation 36:
H = exp(-a-b).* (I_0/beta + (squig.*I_1)./(2.*b.*(1-beta)) );


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
     results_matrix(loopcnt,4) = alpha(k);
     results_matrix(loopcnt,5) = RSS;
     
     
    end
end

plot(PV, C, '-g', 'LineWidth',1.5)
hold on







%%%Ethane
measured = ethane;

R = 1.61;%1.1:0.05:1.7;%1.7;%ethR0;
%heR = 5;        %MUST BE > 1
alpha = 1.7e-4;%10e-5:3e-5:40e-5; %rate constant, try starting in the range 5e-4 to 5e-3

                    %1/s, 100 d^-1 from Geistlinger (2005) = 1e-3 1/s
                    %vulava has 5e-4 - 1e-3  1/s for alpha for SF6, so
                    %probably faster for Helium, N2O
                    %velocity in bioreactor on the order of .003 m/s, much slower than in Vulava, so 
                    %expect lower alpha

%dimensionless parameters
%defined in Van Genuchten pg. 15, with the exception of w, beta, and R
%For w and beta, defined on pg. 28 for one-site kinetic
%non-equilibrium adsorption.  For R, comes from Fry et al. (1994) and Vulava et al. (2002): 
%R = 1 + KH*(Vg/Vw)   %use this with specific K_H values to determine Vg/Vw
%based on fit R parameter

loopcnt = 0;

results_matrix = zeros(length(R)*length(alpha),4);



for j = 1:length(R)
    
    for k = 1:length(alpha)

        
P = v*L/D; %Peclet number
PV = v.*t./L;     %pore volumes, tau in equations in Table 2 in Van Genuchten
w = alpha(k)*L*(R(j)-1)/v;
beta = 1/R(j);

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

for i = 2:length(PV)  %2:length(PV)
    
    tau = (0:delta_tau:PV(i));
    
G_tau = .5 *  erfc( sqrt(P./(4.*beta.*R(j).*tau)).*(beta*R(j)-tau)) + .5*exp(P).*erfc( sqrt (P./(4.*beta.*R(j).*tau) ).*(beta*R(j) + tau) );
    
    a = w.*tau;
    b = w.*(PV(i) - tau)./( (1 - beta)*R(j) );
    
    squig = 2*sqrt(a.*b);
    
%In equation 36, the I_0 and I_1 are modified bessel functions, see
%notation on pg. 50-51
I_0 = besseli(0,squig);  %modified bessel function of the first kind, order 0
I_1 = besseli(1,squig);  %modified bessel function of the first kind, order 1
    
%H term defined by Equation 36:
H = exp(-a-b).* (I_0/beta + (squig.*I_1)./(2.*b.*(1-beta)) );


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
     results_matrix(loopcnt,4) = alpha(k);
     results_matrix(loopcnt,5) = RSS;
     
     
    end
end

plot(PV, C, '-b', 'LineWidth',1.5)
hold on









%%%Methane
measured = ch4;

R = 1.55;%1.1:0.05:1.7;%1.7;%ethR0;
%heR = 5;        %MUST BE > 1
alpha = 2.2e-4;%10e-5:3e-5:40e-5; %rate constant, try starting in the range 5e-4 to 5e-3

                    %1/s, 100 d^-1 from Geistlinger (2005) = 1e-3 1/s
                    %vulava has 5e-4 - 1e-3  1/s for alpha for SF6, so
                    %probably faster for Helium, N2O
                    %velocity in bioreactor on the order of .003 m/s, much slower than in Vulava, so 
                    %expect lower alpha

%dimensionless parameters
%defined in Van Genuchten pg. 15, with the exception of w, beta, and R
%For w and beta, defined on pg. 28 for one-site kinetic
%non-equilibrium adsorption.  For R, comes from Fry et al. (1994) and Vulava et al. (2002): 
%R = 1 + KH*(Vg/Vw)   %use this with specific K_H values to determine Vg/Vw
%based on fit R parameter

loopcnt = 0;

results_matrix = zeros(length(R)*length(alpha),4);



for j = 1:length(R)
    
    for k = 1:length(alpha)

        
P = v*L/D; %Peclet number
PV = v.*t./L;     %pore volumes, tau in equations in Table 2 in Van Genuchten
w = alpha(k)*L*(R(j)-1)/v;
beta = 1/R(j);

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

for i = 2:length(PV)  %2:length(PV)
    
    tau = (0:delta_tau:PV(i));
    
G_tau = .5 *  erfc( sqrt(P./(4.*beta.*R(j).*tau)).*(beta*R(j)-tau)) + .5*exp(P).*erfc( sqrt (P./(4.*beta.*R(j).*tau) ).*(beta*R(j) + tau) );
    
    a = w.*tau;
    b = w.*(PV(i) - tau)./( (1 - beta)*R(j) );
    
    squig = 2*sqrt(a.*b);
    
%In equation 36, the I_0 and I_1 are modified bessel functions, see
%notation on pg. 50-51
I_0 = besseli(0,squig);  %modified bessel function of the first kind, order 0
I_1 = besseli(1,squig);  %modified bessel function of the first kind, order 1
    
%H term defined by Equation 36:
H = exp(-a-b).* (I_0/beta + (squig.*I_1)./(2.*b.*(1-beta)) );


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
     results_matrix(loopcnt,4) = alpha(k);
     results_matrix(loopcnt,5) = RSS;
     
     
    end
end

plot(PV, C, '-c', 'LineWidth',1.5)
hold on









%%%Helium
measured = helium;

R = 3;%1.1:0.05:1.7;%1.7;%ethR0;
%heR = 5;        %MUST BE > 1
alpha = 1.5e-4;%10e-5:3e-5:40e-5; %rate constant, try starting in the range 5e-4 to 5e-3

                    %1/s, 100 d^-1 from Geistlinger (2005) = 1e-3 1/s
                    %vulava has 5e-4 - 1e-3  1/s for alpha for SF6, so
                    %probably faster for Helium, N2O
                    %velocity in bioreactor on the order of .003 m/s, much slower than in Vulava, so 
                    %expect lower alpha

%dimensionless parameters
%defined in Van Genuchten pg. 15, with the exception of w, beta, and R
%For w and beta, defined on pg. 28 for one-site kinetic
%non-equilibrium adsorption.  For R, comes from Fry et al. (1994) and Vulava et al. (2002): 
%R = 1 + KH*(Vg/Vw)   %use this with specific K_H values to determine Vg/Vw
%based on fit R parameter

loopcnt = 0;

results_matrix = zeros(length(R)*length(alpha),4);



for j = 1:length(R)
    
    for k = 1:length(alpha)

        
P = v*L/D; %Peclet number
PV = v.*t./L;     %pore volumes, tau in equations in Table 2 in Van Genuchten
w = alpha(k)*L*(R(j)-1)/v;
beta = 1/R(j);

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

for i = 2:length(PV)  %2:length(PV)
    
    tau = (0:delta_tau:PV(i));
    
G_tau = .5 *  erfc( sqrt(P./(4.*beta.*R(j).*tau)).*(beta*R(j)-tau)) + .5*exp(P).*erfc( sqrt (P./(4.*beta.*R(j).*tau) ).*(beta*R(j) + tau) );
    
    a = w.*tau;
    b = w.*(PV(i) - tau)./( (1 - beta)*R(j) );
    
    squig = 2*sqrt(a.*b);
    
%In equation 36, the I_0 and I_1 are modified bessel functions, see
%notation on pg. 50-51
I_0 = besseli(0,squig);  %modified bessel function of the first kind, order 0
I_1 = besseli(1,squig);  %modified bessel function of the first kind, order 1
    
%H term defined by Equation 36:
H = exp(-a-b).* (I_0/beta + (squig.*I_1)./(2.*b.*(1-beta)) );


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
     results_matrix(loopcnt,4) = alpha(k);
     results_matrix(loopcnt,5) = RSS;
     
     
    end
end

plot(PV, C, '-r', 'LineWidth',1.5)
hold on



ylim([-0.1 1.2])
xlim([0 10])
ylabel('C/C_o [-]')
xlabel('PV')
%legend('Br^-','N_2O', 'SF_6', 'C_2H_6', 'CH_4', 'He', 'Location', 'southeast'); %
legend('Bromide','N_2O', 'SF_6', 'C_2H_6', 'CH_4', 'He', 'Location', 'southeast'); %
set(gca, 'FontSize', 34);




