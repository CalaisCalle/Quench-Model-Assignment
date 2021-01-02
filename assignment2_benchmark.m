clc
clearvars
close all
%% Assignment 1:
% Material: Ti-6Al-4V
%% Things to do:
% 1. Fix the graph
% 2. Save the data to a file 
% 3. Optimise a bit


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. Pre-Processor
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 1a. Geometry
L = 0.03; % m

%% 1b. Material Properties and constants
rho_f = @(T) 7840 + 0.44*T; %kg/m^3
k = @(T) 13.1947 + 0.0126919.*T; % W/m/K
Cp = @(T) 490 + 0.0733333.*T; % J/kg/K
alpha = @(T) k(T) ./ (rho_f(T) .* Cp(T));
h = 50; % W m^-2 K^-1
eps = 1; %emissivity 
sb = 5.67e-8; %steffan-boltzmann constant

%% 1c. Initial Values
T0 = 1000; %Furnace Temperature
T_inf = 20 + 273; % K, far-field temperature 

%% 1d. Boundary Conditions
% 1e. Simulation Controls
m = 100;
n = m;
time = 0;
t_max = 1360;

%% 1f. Discretization
dx = L / (m-1);
dy = L / (n-1);
Xvec = 0:dx:L;          % spatial grid (m)
Yvec = 0:dx:L;          % assuming dx = dy
CFL = 0.5; % Courant-Friedrichs-Lewy condition

%% Plotting parameters
fsize = 14;
mymu = '\mu';


%% 1g. Initialise Temperature distribution and material properties
T = ones(m, n) * T0;
dt = zeros(m, n);

%% Savepoints and saving
xs = [0.5, 0.9, 0.9]*L; % X Points
ys = [0.5, 0.5, 0.9]*L; % Y Points

T_interp = interp2(Xvec,Yvec,T,xs,ys);
save_num = 1;
save_time(save_num) = time;
save_freq = 100;

it = 0; % iteration counter 
save_data_mat = T_interp;
k_mat = zeros(n,m);
Cp_mat = k_mat;
al_mat = k_mat;
Bi_mat = k_mat;
Br_mat = k_mat;
Fo_mat = k_mat;
%% Time Evolution
while time < t_max
    tic
    %update/initialise alpha, Fo, k for this step
    
    for i = 1:n
        for j = 1:m
            k_mat(i,j) = k(T(i,j));
            Cp_mat(i,j) = Cp(T(i,j));
            al_mat(i,j) = k_mat(i,j) ./ (rho_f(T(i,j)) .* Cp_mat(i,j));
            Bi_mat(i,j) = (dx * h) ./ k_mat(i,j);
            Br_mat(i,j) = (eps*sb*dx) ./ k_mat(i,j); %radiation term
        end
    end
    %% Time Step Calculation
    % EDGES:
    for i = 2:n-1
        dt(i,1) = (CFL*dx^2) ./ (2*al_mat(i,1).*(2+Bi_mat(i,1)+Br_mat(i,1).*T(i,1).^3)); %left
        dt(i:n-1,m) = CFL*(dx^2) ./ (2*al_mat(i,m).*(2+Bi_mat(i,m)+Br_mat(i,m).*T(i,m).^3)); %right
        dt(1,i) = CFL*(dx^2) ./ (2*al_mat(1,i).*(2+Bi_mat(1,i)+Br_mat(1,i).*T(1,i).^3)); %top
        dt(n,i) = CFL*(dx^2) ./ (2*al_mat(n,i).*(2+Bi_mat(n,i)+Br_mat(n,i).*T(n,i).^3)); %bot
    end
    
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
    for i = 2:n-1
        for j = 2:n-1
        dt(i,j) = CFL * (dx^2) ./ (4*al_mat(i,j));
        end
    end
    dt_min = min(min(dt));
    
    %update Fo
    for i = 1:n
        for j = 1:m
            Fo_mat(i,j) = (al_mat(i,j) * dt_min) / (dx^2);
        end
    end        
    %% Gauss-Siedel Iterative Method
    T_new = T;
    T_last = T;
    
    err = 1;    
    err_mat = zeros(n, n);
    err_max = 1e-4;
    
    while (err > err_max)
        
        %% Edges

        %Left
        for i = 2:n-1
            edge_part(1) = 2*Bi_mat(i, 1) * T_inf;
            edge_part(2) = 2*Br_mat(i, 1) * T_inf^4;
            edge_part(3) = T_new(i+1,1);
            edge_part(4) = T_new(i-1,1);
            edge_part(5) = 2 * T_new(i, 2);
            A_ii = 1+2*Fo_mat(i,1)*(2+Bi_mat(i,1)+Br_mat(i,1)*T_new(i,1)^3);
            
            T_new(i,1) = (T(i,1) + Fo_mat(i,1) *sum(edge_part)) / A_ii;
        end
        
        %Right
        for i = 2:n-1
            edge_part(1) = 2*Bi_mat(i, m) * T_inf;
            edge_part(2) = 2*Br_mat(i, m) * T_inf^4;
            edge_part(3) = T_new(i+1,m);
            edge_part(4) = 2*T_new(i,m-1);
            edge_part(5) = T_new(i-1, m);
            A_ii = 1+2*Fo_mat(i,m)*(2+Bi_mat(i,m)+Br_mat(i,m)*T_new(i,m)^3);
            
            T_new(i,m) = (T(i,m) + Fo_mat(i,m) * sum(edge_part)) / A_ii;
        end
        %Top
        for j = 2:n-1
            edge_part(1) = 2*Bi_mat(1, j) * T_inf;
            edge_part(2) = 2*Br_mat(1, j) * T_inf^4;
            edge_part(3) = 2*T_new(2,j);
            edge_part(4) = T_new(1,j-1);
            edge_part(5) = T_new(1, j+1);
            A_ii = 1+2*Fo_mat(1,j)*(2+Bi_mat(1,j)+Br_mat(1,j)*T_new(1,j)^3);
            
            T_new(1,j) = (T(1,j) + Fo_mat(1,j) * sum(edge_part)) / A_ii;
        end
        %Bottom
        for j = 2:n-1
            edge_part(1) = 2*Bi_mat(n, j) * T_inf;
            edge_part(2) = 2*Br_mat(n, j) * T_inf^4;
            edge_part(3) = T_new(n,j-1);
            edge_part(4) = 2* T_new(n-1,j);
            edge_part(5) = T_new(n, j+1);
            A_ii = 1+2*Fo_mat(n,j)*(2+Bi_mat(n,j)+Br_mat(n,j)*T_new(n,j)^3);
            
            T_new(n,j) = (T(n,j) + Fo_mat(n,j) * sum(edge_part)) / A_ii;
        end
        
        %% Corners
        %top-left
        crnr_part(1) = 4 * Br_mat(1,1) * T_inf^4;
        crnr_part(2) = 4 * Bi_mat(1,1)*T_inf;
        crnr_part(3) = 2*T_new(2,1);
        crnr_part(4) = 2*T_new(1,2);
        A_ii = 1+4*Fo_mat(1,1)*(1+Bi_mat(1,1)+Br_mat(1,1)*T_new(1,1)^3);
        
        T_new(1,1) = (T(1,1) + Fo_mat(1,1) * sum(crnr_part)) / A_ii;
        
        %top-right
        crnr_part(1) = 4 * Br_mat(1,m) *T_inf^4;
        crnr_part(2) = 4 * Bi_mat(1,m)* T_inf;
        crnr_part(3) = 2*T_new(2,m);
        crnr_part(4) = 2*T_new(1,m-1);
        A_ii = 1+4*Fo_mat(1,m)*(1+Bi_mat(1,m)+Br_mat(1,m)*T_new(1,m)^3);
        
        T_new(1,m) = (T(1,m) + Fo_mat(1,m) * sum(crnr_part)) / A_ii;
        
        %bottom-right (n,m)
        crnr_part(1) = 4 * Br_mat(n,m) *T_inf^4;
        crnr_part(2) = 4 * Bi_mat(n,m)*T_inf;
        crnr_part(3) = 2*T_new(n-1,m);
        crnr_part(4) = 2*T_new(n,m-1);
        A_ii = 1+4*Fo_mat(n,m)*(1+Bi_mat(n,m)+Br_mat(n,m)*T_new(n,m)^3);
        
        T_new(n,m) = (T(n,m) + Fo_mat(n,m) * sum(crnr_part)) / A_ii;
        
        %bottom-left (n, 1)
        crnr_part(1) = 4 * Br_mat(n,1) * T_inf^4;
        crnr_part(2) = 4 * Bi_mat(n,1)*T_inf;
        crnr_part(3) = 2*T_new(n-1,1);
        crnr_part(4) = 2*T_new(n,2);
        A_ii = 1+4*Fo_mat(n,1)*(1+Bi_mat(n,1)+Br_mat(n,1)*T_new(n,1)^3);
        
        T_new(n,1) = (T(n,1) + Fo_mat(n,1) * sum(crnr_part)) / A_ii;
        
        %% interior calculation
        % double letters for original matrices
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
                
                T_new(i, j) = (T(i, j) - Fo_ij * sum(soln_part)) / A_ii;

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
    it = it+1;
    if it >= save_freq
        
        %% Interpolate points and save data
        T_interp = interp2(Xvec,Yvec,T,xs,ys); %Interpolation
        save_num = save_num + 1;
        save_time(save_num) = time;
        save_data_mat(save_num,1:3) = T_interp;
        
    end
    disp(toc)
end

%% 2c. Save Variables

for i = 1:3
    fname = ['ai',num2str(i),'.txt']; %Create File Name
    fileID=fopen(fname,'w');
    for j = 1: SaveNo
        fprintf(fileID,'%f,%f\n',Save_time(j),Save_TempMat(j,i));   % Comma Separated
    end
    fclose(fileID);
end
