clc
clearvars
close all
%% Assignment 1:
% Material: Ti-6Al-4V

%% Geometry
L = 0.25; % m

% Coefficients for material properties
% Can be turned into input parameters for alternate materials
rho_As = [4430 0];
k_As = [1.117 0.0174];
Cp_As = [546.31 0.219];

%% Material Properties Calculations 
% Anonymous functions where temperature is the only variable.

rho = @(T) rho_As(1) + rho_As(2); %kg/m^3
k = @(T) k_As(1) + k_As(2).*T; % W/m/K
Cp = @(T) Cp_As(1) + Cp_As(2).*T; % J/kg/K
alpha = @(T) k(T) ./ (rho(T) .* Cp(T));

% Radiation constants
eps = 0.279; %emissivity 
sb = 5.67e-8; %steffan-boltzmann constant

%% Heat Transfer Coefficients
hs = [150 150 150]; % W m^-2 K^-1
h_top = hs(1);
h_side = hs(2);
h_bot = hs(3);

%% Grid spacing 
n = 50;
m = n;

%% Discretization and grid creation
dx = L / (m-1); % grid spacing
dy = L / (n-1);
Xvec = 0:dx:L;          % spatial grid (m)
Yvec = L:-dx:0;          % assuming dx = dy
CFL = 0.1; % Courant-Friedrichs-Lewy condition
dt = zeros(m, n); %time step matrix

%% Initial temperature values
T0 = 960 + 273.15; %Furnace Temperature
T_inf = 20 + 273.15; % K, far-field temperature 
T = ones(m, n) * T0; % Initial grid

time = 0; % Set Initial Time
t_max = 18000; % Set Maximum time (Could take from file in future)

%% Plotting parameters
fsize = 14;
mymu = '\mu';

%% Interpolation points: Virtual Thermistor Locations for Plottings
xs = [0.5, 0.05, 0.5, 0.95, 0.95, 0.95, 0.5, 0.05, 0.05]*L; % X Points
ys = [0.5, 0.05, 0.05, 0.05, 0.5, 0.95, 0.95, 0.95, 0.5]*L; % Y Points

%% Save Criteria and initial interpolation
T_interp = interp2(Xvec,Yvec,T,xs,ys); %2D interpolation
save_num = 1; % Indexing for save data
save_freq = 10; % Define save points
save_time = time; % Save initial time
save_data_mat = T_interp; % Save initial thermistor temperatures

it = 0; % iteration counter 

%% Matrices for parameters used in calculation
% Attempt at minor optimisation: create matrices only once and overwrite
Bi_mat = zeros(n,n); % Biot Numbers
k_mat = zeros(n,n); % Thermal conductivities matrix
Br_mat = zeros(n,n); % Radiation Terms
al_mat = zeros(n,n); % Alpha terms
Fo_mat = zeros(n,n); % Fourier Numbers 
b = zeros(n,n); % b matrix in Ax = b calculation

    
%% Error conditions for Gauss-Seidel method

err_mat = zeros(n, n);
err_max = 1e-4; % Maximum error allowed
%% Simulation
while max(T(:)) > 1.1*T_inf
    %update/initialise alpha, Fo, k for this step
    k_mat(:,:) = k(T);
    al_mat(:,:) = alpha(T);
    Br_mat(:,:) = (eps*sb*dx) ./ k_mat; 
    
    %% Update Biot matrix using separated HTC's:
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

    %Subsequent corner Calculation assumes Bi total is sum of Bi on two edges 
    %% Time Step Calculation
    
    % EDGES:
        
    %TODO: Make this look cleaner
    dt(2:n-1,1) = (CFL*dx^2) ./ (2*al_mat(2:n-1,1).*(2+Bi_mat(2:n-1,1)+Br_mat(2:n-1,1).*T(2:n-1,1).^3)); %left
    dt(2:n-1,m) = CFL*(dx^2) ./ (2*al_mat(2:n-1,m).*(2+Bi_mat(2:n-1,m)+Br_mat(2:n-1,m).*T(2:n-1,m).^3)); %right
    dt(1,2:m-1) = CFL*(dx^2) ./ (2*al_mat(1,2:n-1).*(2+Bi_mat(1,2:n-1)+Br_mat(1,2:m-1).*T(1,2:n-1).^3)); %top
    dt(n,2:m-1) = CFL*(dx^2) ./ (2*al_mat(n,2:n-1).*(2+Bi_mat(n,2:n-1)+Br_mat(n,2:m-1).*T(n,2:n-1).^3)); %bot
    
    % CORNERS:
    dt(1,1) = CFL*(dx^2) /...
        ((2*al_mat(1,1)) * (2 + Bi_mat(1,1) + 2*Br_mat(1,1)*T(1,1)^3));
    dt(1,m) = CFL*(dx^2) /...
        ((2*al_mat(1,m)) * (2 + Bi_mat(1,m) + 2*Br_mat(1,m)*T(1,m)^3));
    dt(n,1) = CFL*(dx^2) /...
        ((2*al_mat(n,1)) * (2 + Bi_mat(n,1) + 2*Br_mat(n,1)*T(n,1)^3));
    dt(n,m) = CFL*(dx^2) /...
        ((2*al_mat(n,m)) * (2 + Bi_mat(n,m) + 2*Br_mat(n,m)*T(n,m)^3));
    
    
    % INTERIOR:
    dt(2:n-1,2:m-1) = CFL * (dx^2) ./ (4*al_mat(2:n-1,2:m-1));
    
    % Take the minimum time step
    dt_min = min(min(dt));
    
    %update Fo using new time step
    Fo_mat(:,:) = (al_mat .* dt_min) ./ (dx^2);
    
    %% update b matrix
    
    % Edges
    for i = 2:n-1
       b(1,i) = 2*Fo_mat(1,i) .* (Br_mat(1,i) * T_inf^3 + Bi_mat(1,i))*T_inf;  
       b(n,i) = 2*Fo_mat(n,i) .* (Br_mat(n,i) * T_inf^3 + Bi_mat(n,i))*T_inf;  
       b(i,1) = 2*Fo_mat(n,i) .* (Br_mat(n,i) * T_inf^3 + Bi_mat(n,i))*T_inf;  
       b(i,n) = 2*Fo_mat(n,i) .* (Br_mat(n,i) * T_inf^3 + Bi_mat(n,i))*T_inf;  
    end
    
    %corners: 
    b(1,1) = 2*Fo_mat(1,1) * (2*Br_mat(1,1) * T_inf^3 + Bi_mat(1,1)) * T_inf;
    b(1,n) = 2*Fo_mat(1,n) * (2*Br_mat(1,n) * T_inf^3 + Bi_mat(1,n)) * T_inf;
    b(n,1) = 2*Fo_mat(n,1) * (2*Br_mat(n,1) * T_inf^3 + Bi_mat(n,1)) * T_inf;
    b(n,n) = 2*Fo_mat(n,n) * (2*Br_mat(n,n) * T_inf^3 + Bi_mat(n,n)) * T_inf;
    
    b(2:n-1,2:n-1) = 0; %remove old temperatures in interior
    b = b+T; % Each term in b has a Temperature term
    
    %% Gauss-Siedel Iterative Method
    
    % update matrices of old and new temperature for comparison
    T_new = T;
    T_last = T;
    
    err = 1; % Reset Err
    while (err > err_max)
        
        % Roughly: Update values as close to row-by-row as possible
        
        %top-left corner

        crnr_part(1) = 2*T_new(2,1);
        crnr_part(2) = 2*T_new(1,2);
        A_ii = 1+2*Fo_mat(1,1)*(2+Bi_mat(1,1)+2*Br_mat(1,1)*T_new(1,1)^3);
        
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
        A_ii = 1+2*Fo_mat(1,m)*(2+Bi_mat(1,m)+2*Br_mat(1,m)*T_new(1,m)^3);
        
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
        A_ii = 1+2*Fo_mat(n,1)*(2+Bi_mat(n,1)+2*Br_mat(n,1)*T_new(n,1)^3);
        
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
        A_ii = 1+2*Fo_mat(n,m)*(2+Bi_mat(n,m)+2*Br_mat(n,m)*T_new(n,m)^3);
        
        T_new(n,m) = (b(n,m) + Fo_mat(n,m) * sum(crnr_part)) / A_ii;
        
        %% error Calculation
        for i = 1: n
            for j = 1: n
                if T_new(i,j) > 1E-10
                    err_mat(i,j) = (T_new(i,j)-T_last(i,j))/T_new(i,j);
                end
            end
        end
        
        % Take maximum error of anywhere
        err = max(max(abs(err_mat)));
        T_last = T_new; % For next iteration: save current iteration
    end
    
    % Update solution
    T = T_new; % Set T to T_new for plotting
    time = time + dt_min;
    it = it+1; % Increment completed iterations
    
    %% Save Data
    if it >= save_freq
        
        %% Interpolate thermistor temperatures from T
        T_interp = interp2(Xvec,Yvec,T,xs,ys); %Interpolation
        save_num = save_num + 1; % Increment save_num for indexing
        save_time(save_num) = time; % Save time
        save_data_mat(save_num, 1:length(xs)) = T_interp; % Save Temps

        % Plot graphs
%         figure(1)
%         contourf(Xvec,Yvec,T,100,'LineStyle','none')
%         cb = colorbar;
%         pos=get(cb,'Position');
%         set(cb,'Position',pos+[0.1,0,0.01,0.01]); 
%         xlabel(['x distance (m)'],'fontsize',fsize)
%         ylabel(['y distance (m)'],'fontsize',fsize)
%         
%         % Add thermistor data to graph
%         axis equal
%         title(['Time: ',num2str(time),'[s]'],'fontsize',fsize)
%         set(gca,'fontsize',fsize)
%         hold on
%         %add thermistor positions and temperatures as text
%         scatter(xs,ys,'kx') 
%         c = strsplit(num2str(T_interp));
%         text(xs+0.001, ys+0.001, c);
%         hold off
%         
        it = 0; % Reset iteration count to prevent constant saving
    end
end

%% Save Generated Data
for i = 1:9
    % Create numbered file names:
    fname = ['ai',num2str(i),'.txt']; 
    fileID=fopen(fname,'w'); % open file in write mode
    fprintf(fileID, 'Time [s], Temperature [K],\n');
    for j = 1: save_num
        %save each column of data as well as time.
        fprintf(fileID,'%.3f,%.3f\n',save_time(j),save_data_mat(j,i));   
    end
    fclose(fileID);
end
