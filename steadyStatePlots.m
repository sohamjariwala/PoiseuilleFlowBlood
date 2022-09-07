figure('Name','Steady state graphs')

xlabel('Radius, $r$ (m)','Interpreter','latex'); hold on;

plot(r*R0,lambda(end,:),...
    r*R0,M_stressmigration(end,:), ...
    r*R0,C_stressmigration(end,:), ...
    'LineWidth', 2);

box on
xline(R0*rCheb)
legend('\lambda','M','C_{xy}','\gamma')
set(gca,'FontSize',20,'FontWeight','bold','linewidth',2,'TickLabelInterpreter','latex');

figure('Name','Steady state velocity profile')
xlabel('Radius, $r$ (m)','Interpreter','latex'); hold on;
ylabel('$v_z$ (m/s)','Interpreter','latex');
plot(r*R0,vz(end,:),...
    'LineWidth', 2);
box on
xline(R0*rCheb)
set(gca,'FontSize',20,'FontWeight','bold','linewidth',2,'TickLabelInterpreter','latex');

% Contour plot
figure('Name', 'Velocity Contour Plot')
patch_waterfall = contourf(t/omega, r(1:end),lambda(:,1:end)');
xlabel('Time [s]')
ylabel('Radial distance from the center [m]')
zlabel('Velocity [m/s]')
% set(gca,'FontSize',20,'linewidth',2,'DataAspectRatio',[1, L/(2*R0), 1],'Xtick',[])
