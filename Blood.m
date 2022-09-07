clear; close all;

%% Model parameters
%Initialize values with dimensions
omega = 90/60*2*pi;           %angular velocity
amplitude = 1;                    %amplitude
P = -21;                         %pressure drop of the entire tube length
L = 10^(-2);                      %tube length
R0 = 2*10^(-4);                   %tube Radius
Hct = 0.42;                       %Blood Hematocrit
np0 = 1.67*10^(-3);               %plasma viscosity
T0 = 296.16;                      %Reference temperature
Te = T0;                          %Temperature

%according to Casson fit
eta = 0.0031;%np0*(1+2.0703*Hct+3.7222*Hct^2)*exp(-7.0276*(1-T0/Te)); %Viscosity
% eta = 1.25*10^(-3);             %Viscosity
rho = 1200;                       %Density
t = linspace(0,1, 1000)';
tm = 0.1/omega;                         %Characteristic material time
G = 0.019;                        %elastic modulus Stephanou model
tr = 0.14;                       %Tule parameters Staphanou model
tl =0.14;
% tlambda = 1960;                   %Time parameters structural equation
% ta = 483;
% tb = 10;
a1 = 0.1873;
a2 = 2.4654;
b = 0.2072;
eta_0 = 0.0075;
eta_R = 0.0215;
eta_inf = 0.003067;
tau_c = 0.03357;
G = 3.676e-04;
tau_R = 0.0112;

tlambda = 10/omega;
ta = 0.01/omega;
tb = 0.1/omega;
tr = 0.14/omega;

%non-linear models
nCycles = 5;
t = t';                           %dimensionless time
t = nCycles*2*pi*t;                       %dimensionless period
tm = tm*omega;                    %dimensionless Maxwell time
pbar = (P/L)/(eta_inf*omega/R0);    %dimensionless pressure gradient
% pbar = -15;
G = G/(eta_inf*omega);
Wo = sqrt(rho*omega*R0^2/eta_inf);    %Womersley number
p = 1;                            %exponent Carreau model
eta0 = 10;                        %Carreau model zero shear viscosity
params.tr = tr*omega;
params.tl = tl*omega;
params.tlambda = tlambda*omega;
params.ta = ta*omega;
params.tb = tb*omega;
params.kc = 1;
params.kn = 0.5;
params.Wo = Wo;
params.tm = tm;
params.G = G;
params.eta = eta;
params.a1 = a1*omega;
params.a2 = a2*omega;
params.b = b*omega;
params.eta_0 = eta_0/eta_inf;
params.eta_R = eta_R/eta_R;
params.eta_inf = eta_inf/eta_inf;
params.tau_c = tau_c*omega;
params.tau_R = tau_R*omega;
params.phi_hct0 = 0.42;
psi_init = sqrt(params.phi_hct0/(1-params.phi_hct0));

%Add more specific dimensionless model parameters as needed
params_Newtonian = params;
params_Maxwell = params;
params_Stephanou = params;

%--------------------------------------------------------------------------
%% Model and solver conditions
% Number of Chebyshev nodes
m = 12;                            

% Pressure forcing (CHANGE AS NEEDED)
% dpdz = @(t) pbar*amplitude*(1+sin(t));     % Sinusoidal unidirectional pressure field
% dpdz = @(t) pbar*amplitude*(sin(t));     % Sinusoidal pressure field
% dpdz = @(t) (t<pi/2)*pbar*amplitude*(sin(t)) + (t>pi/2)*pbar*amplitude;
% dpdz = @(t) pbar*amplitude*(t- floor(t/10));  % Sawtooth pressure field 
dpdz = @(t) pbar; % Constant profile

% Fluid constitutive models (CHANGE AS NEEDED)
RBCMigrationNeumann = @(t,AA,m,dpdz,B, B1, B2, Bz, B1z, B2z, B1zor, Br, B1r, B2r,B1ror, Bor, B1or, B2or,...
    P, P1, P2, Por, P1or, P2or,r,params) ...
        fRBCMigrationNeumann(t,AA,m,dpdz,B,B1, B2, Bz, B1z, B2z, B1zor, Br, B1r, B2r,B1ror, Bor, B1or, B2or,...
    P, P1, P2, Por, P1or, P2or,r,params);   
    
% Boundary layer thickness (for refinement)
rb = sqrt(0.1);

% Initial condition Chebyshev coefficients 
% (Use fast Chebyshev transform for a specific intial condition)
newtonian_init = omega*R0*[1 zeros(1, m-1)];

%--------------------------------------------------------------------------
%% Initialize Chebyshev domain
r = chebyshevCoefficients(0,1,m)';

% Quadratic bilinear transform
c = (0.25-rb^2)/(0.25*rb^2); % Transform coefficient
r = sqrt((c+1).*r.^2./(c.*r.^2+1)); % Transformed discretized radius

% Generate Chebyshev basis using function values discretized at nodes
[B, B1, B2, Bz, B1z, B2z, B1zor, Br, B1r, B2r,B1ror, Bor, B1or, B2or,...
    P, P1, P2, Por, P1or, P2or] = chebyshevBasis(m,m,r);

%% Solve the time dependent ODE
opts = odeset('RelTol',1e-6,'AbsTol',1/2^(m+1));

RBCmigration_init = [-0.0865/2 zeros(1,m-1), ... %omega*R0*newtonian(end,:),...
                    [0 zeros(1,m-1)],...
                    [params.phi_hct0 zeros(1,m-1)]];

RBCMigration = @(t, AA) RBCMigrationNeumann(t,AA,m,dpdz,B,B1, B2, Bz, B1z, B2z, B1zor, Br, B1r, B2r,B1ror, Bor, B1or, B2or,...
    P, P1, P2, Por, P1or, P2or,r,params);
sol4 = ode15s(RBCMigration, [t(1), t(end)], RBCmigration_init, opts);

%% Obtain solution values at additional radial points
rpoints = 100;
r = linspace(-1, 1, rpoints);

[B, B1, B2, Bz, B1z, B2z, B1zor, Br, B1r, B2r,B1ror, Bor, B1or, B2or,...
    P, P1, P2, Por, P1or, P2or] = chebyshevBasis(m,rpoints, r);

%% Final output variable(s)
RBCmigration = deval(sol4,t)';
vz_stressmigration = omega*R0*RBCmigration(:,1:m)*Bz;
gamma_stressmigration = RBCmigration(:,1:m)*B1z;
phi_stressmigration = (RBCmigration(:,2*m+1:3*m)*P);
phi = phi_stressmigration;
vz = vz_stressmigration;
