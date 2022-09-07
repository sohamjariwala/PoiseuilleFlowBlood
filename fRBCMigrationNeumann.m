function dAA = fRBCMigrationNeumann(t,AA,m,dpdz,B,B1, B2, Bz, B1z, B2z, B1zor, Br, B1r, B2r,B1ror, Bor, B1or, B2or,...
    P, P1, P2, Por, P1or, P2or,r,params)
% Modified function to account for Neumann boundary condition at the walls
% to indicate no flux for volume fraction.
dAA = zeros(size(AA))';

% The solution/idependent variables include Velocity, conformation tensor &
% Hematocrit, respectively along with their bases.
v_c = AA(1:m)'; %Bz
c_c = AA(m+1:2*m)'; %Br
phi_c = AA(2*m+1:3*m)'; %P

% Dependent variables can be expressed in terms of the above
v_z = v_c*Bz; % z direction velocity
gamma_dot = (v_c*B1z); % Shear rate
C_rz =  c_c*Br; % Conformation tensor rz component
phi = phi_c*P; % Hematocrit

% Owen's model for structure
F = params.b*(gamma_dot).^2./(params.a1 + params.a2.*abs(gamma_dot));

M = F/2.*(-1 + realsqrt(1+4.*phi./(F)));
lambda = (phi-M)./phi;

eta_0 = params.eta_0.*(1-phi).^2;
eta_inf = params.eta_inf.*(1-phi).^2;
eta_R = params.eta_R.*(1-phi).^2;

% Stephanou and Geourgiou model for stress
sigma_rz = (params.G.*lambda.*C_rz);

% Cauchy momentum equation
dp = dpdz(t);

solvent_viscosity = eta_R.*lambda.^1.5 ...
          +((eta_0-eta_inf)...
         ./(1+params.tau_c.*abs(gamma_dot))+eta_inf);

beta_rinvSigma_rz = solvent_viscosity.*v_c*B1zor;
beta_dSigma_rz_dr = solvent_viscosity.*v_c*B2z;

rinvSigma_rz = (sigma_rz/Br*B);
dsigma_rzdr = (sigma_rz/B)*B1;
dVzdt = (params.Wo^(-2))*(- dp +  rinvSigma_rz + dsigma_rzdr +  beta_rinvSigma_rz + beta_dSigma_rz_dr)/Bz;

% Conformation tensor evolution equation
dC_rzdt = (-(1-lambda)/params.tau_R.*(C_rz)+ gamma_dot)/Br;

% Hematocrit diffusion equation
omega = 90/60*2*pi;           %angular velocity
R0 = 2*10^(-4);                   %tube Radius
a = 7e-6/R0;
params.kc = params.kc*a^2/4./solvent_viscosity;
D = 5e-10*omega^2/(2e-4);

phi = phi_c*P;
phior = phi_c*Por;
dphidr = phi_c*P1;
d2phidr2 = phi_c*P2;
dphidror = phi_c*P1or;

gamma_DDot = v_c*Bz/B*B2;

sigma_rz = sigma_rz;
                
% Phillips' model
dsigma_rzdr = (sigma_rz/B)*(B1);
dsigma_rzdror = (sigma_rz/B)*B1or;
d2sigma_rzdr2 = (sigma_rz/B)*B2;

% dpsi_dt = (params.kc.*(phior+dphidr).*(dphidr.*sigma_rz+phi.*dsigma_rzdr) ...
%             + params.kc.*phi.*(phi.*(dsigma_rzdror+d2sigma_rzdr2) ...
%             + 2*dphidr.*dsigma_rzdr ...
%             + sigma_rz.*(dphidror+d2phidr2)))/P;

%% Modified Phillips model as per Horner, J.S. Thesis.
% Nc = -params.kc.*phi.*a^2.*(gamma_dot.*phi)/B*B1;
% N_n = -params.kn.*a^2.*gamma_dot.*phi.^2./B*B1;
f_2 = 0.313;

% Malipeddi & Sarkar model
dpsi_dt = (((gamma_dot.*phi.*a^2.*r.*f_2.*dphidr)/P)*P1or)/P;

dpsi_dt = a^2.*f_2*((gamma_dot.*phior.*dphidr + phi.*gamma_DDot.*dphidr ...
               + gamma_dot.*dphidr.^2 + gamma_dot.*phi.*d2phidr2) + D.*(dphidror+d2phidr2))/P;

dAA(1:m) = dVzdt;
dAA(m+1:2*m) = dC_rzdt;
dAA(2*m+1:3*m) = 0*dpsi_dt;
dAA = dAA';

fprintf("\n\nTime %f :", t);
fprintf("\n\nStress: \n");
fprintf("%1.4f ", sigma_rz);
fprintf("\n\nGamma_dot: \n");
fprintf("%1.4f ", gamma_dot);
fprintf("\n\nv_z: \n");
fprintf("%1.4f ", v_z);
fprintf("\n\ndv_zdt: \n");
fprintf("%1.4f ", dVzdt*Bz);
fprintf("\n\ndp: \n");
fprintf("%1.4f ", -dp);
fprintf("\n\nrinvSigma_rz: \n");
fprintf("%1.4f ", rinvSigma_rz);
fprintf("\n\ndsigma_rzdrz: \n");
fprintf("%1.4f ", dsigma_rzdr);
fprintf("\n\n F: \n");
fprintf("%1.4f ", F);
fprintf("\n\nM: \n");
fprintf("%1.4f ", M);
fprintf("\n\nLambda: \n");
fprintf("%1.4f ", lambda);
fprintf("\n\nC_rz: \n");
fprintf("%1.4f ", C_rz);
fprintf("\n\nPhi: \n");
fprintf("%1.8f ", phi);
fprintf("\n\ndphi_dt: \n");
fprintf("%1.8f ", dpsi_dt*P);

% figure(1); plot(r, lambda); hold on
% figure(2); plot(r, sigma_rz); hold on
% figure(4); plot(r, phi); hold on
end