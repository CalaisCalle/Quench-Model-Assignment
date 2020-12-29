clc
clearvars
close all
%% Assignment 1:
% Material: Ti-6Al-4V
%% Things to do:
% 1. Change to using n and m for indexing (not just n)
% 2. Fix the interior iterations
% 3. add interpolated nodes
% 4. Add nodes from 3. to graph
% 5. save the results
% 6. Fix the indexing for the extended array
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. Pre-Processor
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 1a. Geometry
L = 0.25; % m

%% 1b. Material Properties and constants
rho = 4430; %kg/m^3
k = @(T) 1.116 + 0.0174.*T; % W/m/K
Cp = @(T) 546.31 + 0.219.*T; % J/kg/K
alpha = @(T) k(T) ./ (rho .* Cp(T));
h = 150; % W m^-2 K^-1
eps = 0.279; %emissivity 
sb = 5.67e-8; %steffan-boltzmann constant

%% 1c. Initial Values
T0 = 960; %Furnace Temperature
T_inf = 20 + 273.15; % K, far-field temperature 

%% 1d. Boundary Conditions
% 1e. Simulation Controls
m = 100;
n = m;
time = 0;
t_max = 200;
 
%% 1f. Discretization
dx = L / (m-1);
dy = L / (n-1);
Xvec = 0:dx:L;          % spatial grid (m)
Yvec = 0:dx:L;          % assuming dx = dy
CFL = 0.1; % Courant-Friedrichs-Lewy condition

%% Plotting parameters
fsize = 14;
mymu = '\mu';


%% 1g. Initialise Temperature distribution and material properties
T = ones(m, n) * T0;
dt = zeros(m, n);
k_mat = k(T);
alpha_mat = alpha(T);
Cp_mat = Cp(T);
Fo_mat = alpha(T) .* (dt / (dx^2));
Bi_mat = dx * h ./ k_mat;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2. Solver
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 2a. Solve the partial differential equation

%% Time Step and Model Initialisation


%% Time Evolution
while time < t_max
    T_new = T;
    T_last = T;
    %update alpha, Fo, k for this step
    k_mat = k(T);
    alpha_mat = k_mat ./ (rho .* Cp_mat);

    Bi_mat = dx * h ./ k_mat;
    
    %% Time Step Calculation
    % EDGES:
    dt(2:n-1,1) = CFL*(dx^2) ./ (alpha_mat(2:n-1,1).*(4+2*Bi_mat(2:n-1,1))); %left
    dt(2:n-1,m) = CFL*(dx^2) ./ (alpha_mat(2:n-1,1).*(4+2*Bi_mat(2:n-1,1))); %right
    dt(1,2:m-1) = CFL*(dx^2) ./ (alpha_mat(2:n-1,1).*(4+2*Bi_mat(2:n-1,1))); %top
    dt(n,2:m-1,1) = CFL*(dx^2) ./ (alpha_mat(2:n-1,1).*(4+2*Bi_mat(2:n-1,1))); %bot
    
    % CORNERS:
    dt(1,1) = CFL*(dx^2) /...
        ((4*alpha_mat(1,1)) * (1 + Bi_mat(1,1)*(1 + (eps*sb/h)*(T(1,1)^3))));
    dt(1,m) = CFL*(dx^2) /...
        ((4*alpha_mat(1,m)) * (1 + Bi_mat(1,m)*(1 + (eps*sb/h)*(T(1,m)^3))));
    dt(n,1) = CFL*(dx^2) /...
        ((4*alpha_mat(n,1)) * (1 + Bi_mat(n,1)*(1 + (eps*sb/h)*(T(n,1)^3))));
    dt(n,m) = CFL*(dx^2) /...
        ((4*alpha_mat(n,m)) * (1 + Bi_mat(n,m)*(1 + (eps*sb/h)*(T(n,m)^3))));
    
    
    % INTERIOR:
    
    dt(2:n-1,2:m-1) = CFL * (dx^2) ./ (4*alpha_mat(2:n-1,2:m-1));
    dt_min = min(min(dt));
    Fo_mat = (alpha_mat .* dt_min) / (dx^2);
    %% Gauss-Siedel Iterative Method
    err = 1;    
    err_mat = zeros(n, n);
    err_max = 1e-4;
    
    while (err > err_max)
        
        
        %% Corners
        %top-left
        crnr_part(1) = 4*sb*eps*Bi_mat(1,1)*(T_inf^4 - T_new(1,1)^4)/h;
        crnr_part(2) = 4 * Bi_mat(1,1)*T_inf;
        crnr_part(3) = 2*T_new(2,1);
        crnr_part(4) = 2*T_new(1,2);
        A_ij = 1 + 2*Fo_mat(1,1) + 4*Bi_mat(1,1)*Fo_mat(1,1);
        
        T_new(1,1) = (T(1,1) + Fo_mat(1,1) * sum(crnr_part)) / A_ij;
        
        %top-right
        crnr_part(1) = 4*sb*eps*Bi_mat(1,m)*(T_inf^4 - T_new(1,m)^4)/h;
        crnr_part(2) = 4 * Bi_mat(1,m)*T_inf;
        crnr_part(3) = 2*T_new(2,m);
        crnr_part(4) = 2*T_new(1,m-1);
        A_ij = 1 + 2*Fo_mat(1,m) + 4*Bi_mat(1,m)*Fo_mat(1,m);
        
        T_new(1,m) = (T(1,m) + Fo_mat(1,m) * sum(crnr_part)) / A_ij;
        
        %bottom-right (n,m)
        crnr_part(1) = 4*sb*eps*Bi_mat(n,m)*(T_inf^4 - T_new(n,m)^4)/h;
        crnr_part(2) = 4 * Bi_mat(n,m)*T_inf;
        crnr_part(3) = 2*T_new(n-1,m);
        crnr_part(4) = 2*T_new(n,m-1);
        A_ij = 1 + 2*Fo_mat(n,m) + 4*Bi_mat(n,m)*Fo_mat(n,m);
        
        T_new(n,m) = (T(n,m) + Fo_mat(n,m) * sum(crnr_part)) / A_ij;
        
        %bottom-left (n, 1)
        crnr_part(1) = 4*sb*eps*Bi_mat(n,1)*(T_inf^4 - T_new(n,1)^4)/h;
        crnr_part(2) = 4 * Bi_mat(n,1)*T_inf;
        crnr_part(3) = 2*T_new(n-1,1);
        crnr_part(4) = 2*T_new(n,2);
        A_ij = 1 + 2*Fo_mat(n,1) + 4*Bi_mat(n,1)*Fo_mat(n,1);
        
        T_new(n,1) = (T(n,1) + Fo_mat(n,1) * sum(crnr_part)) / A_ij;
        
        %% Edges
        % Had to increment/decrement some iterators thanks to funky
        % stuff happening with the extended matrix
        %Left
        for i = 2:n-1
            edge_part(1) = 2*Bi_mat(i, 1) * T_inf;
            edge_part(2) = 2*sb*eps*Bi_mat(i,1)*(T_inf^4 - T_new(i,1)^4) / h;
            edge_part(3) = T_new(i+1,1);
            edge_part(4) = T_new(i-1,1);
            edge_part(5) = 2 * T_new(i, 2);
            A_ij = 1 + 4*Fo_mat(i, 1) + 2*Bi_mat(i,1);
            
            T_new(i,1) = (T(i,1) + Fo_mat(i,1) *sum(edge_part)) / A_ij;
        end
        
        %Right
        for i = 2:n-1
            edge_part(1) = 2*Bi_mat(i, m) * T_inf;
            edge_part(2) = 2*sb*eps*Bi_mat(i,m)*(T_inf^4 - T_new(i,m)^4) / h;
            edge_part(3) = T_new(i+1,m);
            edge_part(4) = 2*T_new(i,m-1);
            edge_part(5) = T_new(i-1, m);
            A_ij = 1 + 4*Fo_mat(i, m) + 2*Bi_mat(i,m);
            
            T_new(i,m) = (T(i,m) + Fo_mat(i,m) * sum(edge_part)) / A_ij;
        end
        %Top
        for j = 2:n-1
            edge_part(1) = 2*Bi_mat(1, j) * T_inf;
            edge_part(2) = 2*sb*eps*Bi_mat(1,j)*(T_inf^4 - T_new(1,j)^4) / h;
            edge_part(3) = 2*T_new(2,j);
            edge_part(4) = T_new(1,j-1);
            edge_part(5) = T_new(1, j+1);
            A_ij = 1 + 4*Fo_mat(1, j) + 2*Bi_mat(1,j);
            
            T_new(1,j) = (T(1,j) + Fo_mat(1,j) * sum(edge_part)) / A_ij;
        end
        %Bottom
        for j = 2:n-1
            edge_part(1) = 2*Bi_mat(n, j) * T_inf;
            edge_part(2) = 2*sb*eps*Bi_mat(n,j)*(T_inf^4 - T_new(n,j)^4) / h;
            edge_part(3) = T_new(n,j-1);
            edge_part(4) = 2* T_new(n-1,j);
            edge_part(5) = T_new(n, j+1);
            A_ij = 1 + 4*Fo_mat(n, j) + 2*Bi_mat(n,j);
            
            T_new(n-1,j) = (T(n,j) + Fo_mat(n,j) * sum(edge_part)) / A_ij;
        end
        
        %% interior calculation
        % double letters for original matrices
        for i = 2:n-1
            for j = 2:n-1
                %calculate new beta values (THESE MIGHT BE WRONG)
                beta_x = (k_mat(i, j + 1) - k_mat(i, j - 1)) / k_mat(i, j);
                beta_y = (k_mat(i + 1, j) - k_mat(i - 1, j)) /k_mat(i, j);
                %calculate new temperatures using old Fo
                Fo_ij = Fo_mat(i, j);
                A_ij = 1 + 4*Fo_ij;
                
                soln_part(1) = beta_x*(T_new(i, j-1) - T_new(i, j+1)) / 4;
                soln_part(2) = beta_y*(T_new(i-1, j) - T_new(i+1, j)) / 4;
                soln_part(3) = T_new(i, j+1) + T_new(i, j-1);
                soln_part(4) = T_new(i+1, j) + T_new(i-1, j);
                
                T_new(i, j) = (T(i, j) - Fo_ij * sum(soln_part)) / A_ij;

            end 
        end
        %% error Calculation
        for i = 1: n
            for j = 1: n
                if T_new(i,j) > 1E-10
                    err_mat(i,j) = (T_new(i,j)-T_last(i,j))/T_new(i,j);
                end
            end
        end
        %% 2b. Update Solution
        err = max(max(abs(err_mat)));
        T_last = T_new; 
    end
    
    T = T_new;
    time = time + dt_min;
    
    %% Plot graphs
    figure(1)
    contourf(Xvec,Yvec,T)
    colorbar
    xlabel(['x distance (',mymu,'m)'],'fontsize',fsize)
    ylabel(['y distance (',mymu,'m)'],'fontsize',fsize)
    axis equal
    title('Numerical Solution (K)','fontsize',fsize)
    set(gca,'fontsize',fsize)


    %Save graphs to file
    %Use interpolation to determine thermistor values?
end

%% 2c. Save Variables

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3. Post-Processor
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 3a. Display Results

%% 3b. Solution