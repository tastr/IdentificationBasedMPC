tmp = load('./simulations/sols4paperPlot');

sol_dd = tmp.sols4paperPlot.sol_dd;
sol_lsq = tmp.sols4paperPlot.sol_lsq;
sol_mdl =  tmp.sols4paperPlot.sol_mld;



%% calculate the values of u 

u_clMDL = zeros(1,length(sol_mdl.u_cl));
u_clMDL(1) =sol_dd.u_cl(1);
for i = 2:length(u_clMDL)
    u_clMDL(i) = u_clMDL(i-1) + sol_mdl.u_cl(i);
end

u_clLSQ = zeros(1,length(sol_lsq.u_cl));
u_clLSQ(1) =sol_dd.u_cl(1);
for i = 2:length(u_clLSQ)
    u_clLSQ(i) = u_clLSQ(i-1) + sol_lsq.u_cl(i);
end

%%
close all

lgnd = {};
subplot(2,1,1)

hold on 
plot(sol_dd.x_cl(2,:));
lgnd{end+1} = "Data-driven MPC with $N = " + num2str(sol_dd.params.N) + "$"
plot(sol_lsq.x_cl(2,:));
lgnd{end+1} = "Identification-based MPC with $N = " + num2str(sol_lsq.params.N) + "$"

plot(sol_mdl.x_cl(2,:))
lgnd{end+1} = "Model-based MPC"

plot(1:length(sol_dd.x_cl(2,:)),sol_dd.y_T(1,:)*ones(1,length(sol_dd.x_cl(2,:))), '--');
lgnd{end+1} = "Setpoint $y^{\mathrm{r}}$"
hold off


xlim([0,2500])
xlabel('Time-step $t$')
ylabel('Solution $y(t)$')
legend(lgnd,'Location','southeast' )
grid on 

subplot(2,1,2)
plot(sol_dd.u_cl(1,:));
hold on 
plot(u_clLSQ)
plot(u_clMDL)
hold off

xlim([0,2500])
xlabel('Time-step $t$')
ylabel('Input $u(t)$')
legend(lgnd{1:end-1},'Location','southeast' )
grid on 
matlab2tikz('solExample.tex')
%exportgraphics(gcf,"./sols4paper/solExample.jpg")