close;
v = VideoWriter('./migrationVz','MPEG-4');
vFrameRate = 60;
open(v);
tm = t/omega;
velocity = vz(1:end,:);

figure('Position',[1 0 720 720], ...
    'MenuBar','none','ToolBar','none','resize','off');

set(gca,'FontSize',20,'FontWeight','bold','linewidth',2,'TickLabelInterpreter','latex');

axis([-R0, R0,min(velocity,[],'all'), max(velocity,[],'all')]); box on;
% xline(R0*rCheb)
ylabel('Velocity $v_z$ (m/s)','Interpreter','latex');
xlabel('Radius, $r$ (m)','Interpreter','latex'); hold on;

for i = 1:length(tm)
    title(['Time  = ' num2str(tm(i),'%.2f') ' (s)'],'Interpreter','latex')

    velplot1 = semilogx(r*R0,velocity(i,:), 'b', 'LineWidth', 2);

    frame = getframe(gcf);
    writeVideo(v,frame);

    delete(velplot1)
end

close(v)