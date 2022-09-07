%% Velocity contour plot
figure('Name', 'Velocity Contour Plot')
contourf(t/omega, r(2:end),omega*R0*vz(:,2:end)')
xlabel('Time [s]')
ylabel('Radial distance from the center [m]')
zlabel('Velocity [m/s]')
set(gca,'FontSize',20,'linewidth',2,'DataAspectRatio',[2*R0, L/4, 1])