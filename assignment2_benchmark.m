clc
clearvars
close all
%% Assignment 1:
% Material: Ti-6Al-4V
%% Things to do:
% 1. Make the sections actually sensible
% 2. Optimise a bit


%% Geometry
L = 0.03; % m

rho_As = [7840 0.44];
k_As = [13.1947 0.0126919];
Cp_As = [490 0.0733333];

%% Material Properties and constants
rho = @(T) rho_As(1) + rho_As(2); %kg/m^3
k = @(T) k_As(1) + k_As(2).*T; % W/m/K
Cp = @(T) Cp_As(1) + Cp_As(2).*T; % J/kg/K
alpha = @(T) k(T) ./ (rho(T) .* Cp(T)); %Unitless (I think)

eps = 1; %emissivity 
sb = 5.67e-8; %steffan-boltzmann constant

%% Parameter of Interest
hs = [50 50 50]; % W m^-2 K^-1
h_top = hs(1);
h_side = hs(2);
h_bot = hs(3);

%% Initial Values
T0 = 1000; %Furnace Temperature
T_inf = 293; % K, far-field temperature 

% Simulation Controls
n = 50;
m = n;
time = 0;
t_max = 18000;

%% 1f. Discretization
dx = L / (m-1);
dy = L / (n-1);
Xvec = 0:dx:L;          % spatial grid (m)
Yvec = L:-dx:0;          % assuming dx = dy
CFL = 0.1; % Courant-Friedrichs-Lewy condition

% %% Plotting parameters
% fsize = 14;
% mymu = '\mu';


%% 1g. Initialise Temperature distribution and material properties
T = ones(m, n) * T0;
dt = zeros(m, n);
b = zeros(m,n); %b matrix 

%% Interpolation points
xs = [0.015 0.027 0.027]; % X Points
ys = [0.015 0.015 0.027]; % Y Points

%% Interpolation points and save criteria
T_interp = interp2(Xvec,Yvec,T,xs,ys);
save_num = 1;
save_freq = 10;
save_times = time;
save_data_mat = T_interp;

it = 0; % iteration counter 

Bi_mat = zeros(n,n);
%% Time Evolution
while max(T(:)) > 1.1*T_inf
    %update/initialise alpha, Fo, k for this step
    k_mat = k(T);
    al_mat = alpha(T);
    Br_mat = (eps*sb*dx) ./ k_mat; %radiation term
    
    %% Update Biot matrix:
    % Bi differs depending on position on the grid
    % Edges are based on just one of h_top, h_side, or h_bot
    % Corners are based on a sum of 2 of h_top, h_side, and h_bot
    % This is mostly so I don't have to change too much of the code
    %Edges:
    Bi_mat(1,2:n-1) = (dx * h_top) ./ k_mat(1,2:n-1);%top
    Bi_mat(n,2:n-1) = (dx * h_bot) ./ k_mat(n,2:n-1);%bottom
    Bi_mat(2:n-1,1) = (dx * h_side) ./ k_mat(2:n-1,1);%left
    Bi_mat(2:n-1,m) = (dx * h_side) ./ k_mat(2:n-1,m); %right
    
    %Corners:
    Bi_mat(1,1) = (dx * (h_top + h_side)) / (2*k_mat(1,1));
    Bi_mat(1,m) = (dx * (h_top + h_side)) / (2*k_mat(1,m));
    Bi_mat(n,1) = (dx * (h_bot + h_side)) / (2*k_mat(n,1));
    Bi_mat(n, m) = (dx * (h_top + h_side)) / (2*k_mat(n,m));
    % Divided by two so that Bi(corner) = 0.5 (B(side) + B(T/B))
    % It just means the rest of the maths can stay as if only one value of
    % h existed
    
    %% Time Step Calculation
    % EDGES:
    
    %TODO: Make this look cleaner
    dt(2:n-1,1) = (CFL*dx^2) ./ (2*al_mat(2:n-1,1).*(2+Bi_mat(2:n-1,1)+Br_mat(2:n-1,1).*T(2:n-1,1).^3)); %left
    dt(2:n-1,m) = CFL*(dx^2) ./ (2*al_mat(2:n-1,m).*(2+Bi_mat(2:n-1,m)+Br_mat(2:n-1,m).*T(2:n-1,m).^3)); %right
    dt(1,2:m-1) = CFL*(dx^2) ./ (2*al_mat(1,2:n-1).*(2+Bi_mat(1,2:n-1)+Br_mat(1,2:m-1).*T(1,2:n-1).^3)); %top
    dt(n,2:m-1) = CFL*(dx^2) ./ (2*al_mat(n,2:n-1).*(2+Bi_mat(n,2:n-1)+Br_mat(n,2:m-1).*T(n,2:n-1).^3)); %bot
    % CORNERS:
    dt(1,1) = CFL*(dx^2) /...
        ((4*al_mat(1,1)) * (1 + Bi_mat(1,1) + Br_mat(1,1)*T(1,1)^3));
    dt(1,m) = CFL*(dx^2) /...
        ((4*al_mat(1,m)) * (1 + Bi_mat(1,m) + Br_mat(1,m)*T(1,m)^3));
    dt(n,1) = CFL*(dx^2) /...
        ((4*al_mat(n,1)) * (1 + Bi_mat(n,1) + Br_mat(n,1)*T(n,1)^3));
    dt(n,m) = CFL*(dx^2) /...
        ((4*al_mat(n,m)) * (1 + Bi_mat(n,m) + Br_mat(n,m)*T(n,m)^3));
    
    
    % INTERIOR:
    
    dt(2:n-1,2:m-1) = CFL * (dx^2) ./ (4*al_mat(2:n-1,2:m-1));
    dt_min = min(min(dt));
    
    %update Fo
    Fo_mat = (al_mat .* dt_min) ./ (dx^2);
    
    %% update b matrix
    
    %edges and corners:
    b(1,1:n) = 2*Fo_mat(1,1:n).*(Bi_mat(1,1:n)+Br_mat(1,1:n).*T_inf^3)*T_inf;
    b(n,1:n) = 2*Fo_mat(n,1:n).*(Bi_mat(n,1:n)+Br_mat(n,1:n).*T_inf^3)*T_inf;
    b(1:n,1) = 2*Fo_mat(1:n,1).*(Bi_mat(1:n,1)+Br_mat(1:n,1).*T_inf^3)*T_inf;
    b(1:n,n) = 2*Fo_mat(1:n,n).*(Bi_mat(1:n,n)+Br_mat(1:n,n).*T_inf^3)*T_inf;
    
    %For corners: just multiply all terms but temperature by 2
    b(1,1) = b(1,1) * 2;
    b(1,n) = b(1,n) * 2;
    b(n,1) = b(n,1) * 2;
    b(n,n) = b(n,n) * 2;
    
    % add T to all of it
    b(2:n-1,2:n-1) = 0; %This is a bodge to make sure old temps don't stack
    b = b+T;
    
    %% Gauss-Siedel Iterative Method
    T_new = T;
    T_last = T;
    
    % Error conditions
    err = 1;    
    err_mat = zeros(n, n);
    err_max = 1e-4;
    tic
    while (err > err_max)
        
        %top-left corner

        crnr_part(1) = 2*T_new(2,1);
        crnr_part(2) = 2*T_new(1,2);
        A_ii = 1+4*Fo_mat(1,1)*(1+Bi_mat(1,1)+Br_mat(1,1)*T_new(1,1)^3);
        
        T_new(1,1) = (b(1,1) + Fo_mat(1,1) * sum(crnr_part)) / A_ii;
        
        %Top edge
        for j = 2:n-1
            edge_part(1) = 2*T_new(2,j);
            edge_part(2) = T_new(1,j-1);
            edge_part(3) = T_new(1, j+1);
            A_ii = 1+2*Fo_mat(1,j)*(2+Bi_mat(1,j)+Br_mat(1,j)*T_new(1,j)^3);
            
            T_new(1,j) = (b(1,j) + Fo_mat(1,j) * sum(edge_part)) / A_ii;
        end
        
        %top-right corner (1,m)
        crnr_part(1) = 2*T_new(2,m);
        crnr_part(2) = 2*T_new(1,m-1);
        A_ii = 1+4*Fo_mat(1,m)*(1+Bi_mat(1,m)+Br_mat(1,m)*T_new(1,m)^3);
        
        T_new(1,m) = (b(1,m) + Fo_mat(1,m) * sum(crnr_part)) / A_ii;
        
        %Left Edge
        for i = 2:n-1
            edge_part(1) = T_new(i+1,1);
            edge_part(2) = T_new(i-1,1);
            edge_part(3) = 2 * T_new(i, 2);
            A_ii = 1+2*Fo_mat(i,1)*(2+Bi_mat(i,1)+Br_mat(i,1)*T_new(i,1)^3);
            
            T_new(i,1) = (b(i,1) + Fo_mat(i,1) *sum(edge_part)) / A_ii;
        end
        %Right Edge
        for i =  2:n-1
            edge_part(1) = T_new(i+1,m);
            edge_part(2) = 2*T_new(i,m-1);
            edge_part(3) = T_new(i-1, m);
            A_ii = 1+2*Fo_mat(i,m)*(2+Bi_mat(i,m)+Br_mat(i,m)*T_new(i,m)^3);
            
            T_new(i,m) = (b(i,m) + Fo_mat(i,m) * sum(edge_part)) / A_ii;
        end
        
        %% interior calculation
        for i = 2:n-1
            for j = 2:n-1
                %calculate new beta values (THESE MIGHT BE WRONG)
                beta_x = (k_mat(i,j + 1) - k_mat(i,j - 1)) / (4*k_mat(i, j));
                beta_y = (k_mat(i + 1, j) - k_mat(i - 1, j)) /(4*k_mat(i, j));
                %calculate new temperatures using old Fo
                Fo_ij = Fo_mat(i, j);
                A_ii = 1 + 4*Fo_ij;
                
                soln_part(1) = -(beta_x + 1) * T_new(i,j+1);
                soln_part(2) = (beta_x - 1) * T_new(i,j-1);
                soln_part(3) = -(beta_y + 1) * T_new(i+1, j);
                soln_part(4) = (beta_y - 1) * T_new(i-1,j);
                
                %Don't need to use b, but doing so for consistency
                T_new(i, j) = (b(i,j) - Fo_ij * sum(soln_part)) / A_ii;

            end 
        end
        

        %bottom-left corner (n, 1)
        crnr_part(1) = 2*T_new(n-1,1);
        crnr_part(2) = 2*T_new(n,2);
        A_ii = 1+4*Fo_mat(n,1)*(1+Bi_mat(n,1)+Br_mat(n,1)*T_new(n,1)^3);
        
        T_new(n,1) = (b(n,1) + Fo_mat(n,1) * sum(crnr_part)) / A_ii;
        %Bottom Edge
        for j = 2:n-1
            edge_part(1) = T_new(n,j-1);
            edge_part(2) = 2* T_new(n-1,j);
            edge_part(3) = T_new(n, j+1);
            A_ii = 1+2*Fo_mat(n,j)*(2+Bi_mat(n,j)+Br_mat(n,j)*T_new(n,j)^3);
            
            T_new(n,j) = (b(n,j) + Fo_mat(n,j) * sum(edge_part)) / A_ii;
        end

        %bottom-right corner (n,m)
        crnr_part(1) = 2*T_new(n-1,m);
        crnr_part(2) = 2*T_new(n,m-1);
        A_ii = 1+4*Fo_mat(n,m)*(1+Bi_mat(n,m)+Br_mat(n,m)*T_new(n,m)^3);
        
        T_new(n,m) = (b(n,m) + Fo_mat(n,m) * sum(crnr_part)) / A_ii;
        
        %% error Calculation
        for i = 1: n
            for j = 1: n
                if T_new(i,j) > 1E-10
                    err_mat(i,j) = (T_new(i,j)-T_last(i,j))/T_new(i,j);
                end
            end
        end

        err = max(max(abs(err_mat)));
        T_last = T_new; 
    end
    % Update solution
    T = T_new;
    time = time + dt_min;
    it = it+1;
    if it >= save_freq
        %% Interpolate points and save data
        T_interp = interp2(Xvec,Yvec,T,xs,ys); %Interpolation
        save_num = save_num + 1;
        save_times(save_num) = time;
        save_data_mat(save_num, 1:length(xs)) = T_interp;

%         %% Plot graphs
%         figure(1)
%         contourf(Xvec,Yvec,T,100,'LineStyle','none')
%         cb = colorbar;
%         pos=get(cb,'Position');
%         set(cb,'Position',pos+[0.11,0,0,0]); 
%         xlabel(['x distance (m)'],'fontsize',fsize)
%         ylabel(['y distance (m)'],'fontsize',fsize)
%         
%         axis equal
%         title(['Time: ',num2str(time),'[s]'],'fontsize',fsize)
%         set(gca,'fontsize',fsize)
%         hold on
%         scatter(xs,ys,'kx') %Identify points of interest
%         c = strsplit(num2str(T_interp));
%         dx0 = 0.001;
%         dy0 = dx0;
%         text(xs+dx0, ys+dy0, c);
%         hold off
        it = 0;
    end
end

%% Benchmark against given data
bench_data = load('Benchmark_data.txt');
bench_time = bench_data(:,1);
thermo_data = bench_data(:, 2:end);

 % interpolate
T_intp = interp1(save_times,save_data_mat,bench_time);
pred_err = abs(T_intp - thermo_data) ./ thermo_data(i);

errsum = 0;
[rows, cols] = size(thermo_data);
for j = 1:cols  

   err_intgrt = trapz(bench_time, pred_err(:,j)); %summed error
   errsum = errsum + err_intgrt; 
end

disp(errsum / bench_time(end));