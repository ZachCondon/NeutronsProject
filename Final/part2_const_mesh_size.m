M = 5;
flux_t_M5 = readmatrix('Part2_flux_from_transport_norm_convergence_M5.csv');
x_t_M5 = linspace(-M,M,100*M/5);
flux_d_M5 = readmatrix('Part2_flux_from_diffusion_norm2_convergence_M5.csv');
x_d_M5 = linspace(-(M+2/3),(M+2/3),100*M/5);
flux_d_interp_M5 = zeros(1,100*M/5);
for i = 1:length(x_t_M5)
    x = zeros(1,2);
    y = zeros(1,2);
    x(1) = max(x_d_M5(x_d_M5 < x_t_M5(i)));
    x(2) = min(x_d_M5(x_d_M5(:) > x_t_M5(i)));
    y1 = find(x_d_M5(:)==x(1));
    y(1) = flux_d_M5(y1);
    y2 = find(x_d_M5(:)==x(2));
    y(2) = flux_d_M5(y2);
    flux_d_interp_M5(i) = interp1(x,y,x_t_M5(i));
end

figure()
plot(x_d_M5,flux_d_M5,x_t_M5,flux_t_M5,'linewidth',2)
title({"Part 2 - Final Project","Comparing Diffusion to Transport","M="+M})
ylabel('Flux $\phi$ ($n/(cm^2 s)$)','Interpreter','latex')
xlabel('Position $(cm)$','Interpreter','latex')
legend('diffusion','transport','location','best')
xlim([-(M+2/3) (M+2/3)])
ylim([0 2])
norm2_flux(M/5) = norm(abs(flux_d_M5-flux_t_M5),2)/norm(flux_t_M5,2);
fprintf('M5 norm: %4.2f\n',norm2_flux)
filename = ['Flux_comparison_M_' num2str(M) '.jpg'];
saveas(gcf,filename)

%%
M = 10;
flux_t_M10 = readmatrix('Part2_flux_from_transport_norm_convergence_M10.csv');
x_t_M10 = linspace(-M,M,100*M/5);
flux_d_M10 = readmatrix('Part2_flux_from_diffusion_norm2_convergence_M10.csv');
x_d_M10 = linspace(-(M+2/3),(M+2/3),100*M/5);
flux_d_interp_M10 = zeros(1,100*M/5);
for i = 1:length(x_t_M10)
    x = zeros(1,2);
    y = zeros(1,2);
    x(1) = max(x_d_M10(x_d_M10 < x_t_M10(i)));
    x(2) = min(x_d_M10(x_d_M10(:) > x_t_M10(i)));
    y1 = find(x_d_M10(:)==x(1));
    y(1) = flux_d_M10(y1);
    y2 = find(x_d_M10(:)==x(2));
    y(2) = flux_d_M10(y2);
    flux_d_interp_M10(i) = interp1(x,y,x_t_M10(i));
end

figure()
plot(x_d_M10,flux_d_M10,x_t_M10,flux_t_M10,'linewidth',2)
title({"Part 2 - Final Project","Comparing Diffusion to Transport","M="+M})
ylabel('Flux $\phi$ ($n/(cm^2 s)$)','Interpreter','latex')
xlabel('Position $(cm)$','Interpreter','latex')
legend('diffusion','transport','location','best')
xlim([-(M+2/3) (M+2/3)])
ylim([0 2])
norm2_flux(M/5) = norm(abs(flux_d_M10-flux_t_M10),2)/norm(flux_t_M10,2);
fprintf('M10 norm: %4.2f\n',norm2_flux)
filename = ['Flux_comparison_M_' num2str(M) '.jpg'];
saveas(gcf,filename)

%%
M = 15;
flux_t_M15 = readmatrix('Part2_flux_from_transport_norm_convergence_M15.csv');
x_t_M15 = linspace(-M,M,100*M/5);
flux_d_M15 = readmatrix('Part2_flux_from_diffusion_norm2_convergence_M15.csv');
x_d_M15 = linspace(-(M+2/3),(M+2/3),100*M/5);
flux_d_interp_M15 = zeros(1,100*M/5);
for i = 1:length(x_t_M15)
    x = zeros(1,2);
    y = zeros(1,2);
    x(1) = max(x_d_M15(x_d_M15 < x_t_M15(i)));
    x(2) = min(x_d_M15(x_d_M15(:) > x_t_M15(i)));
    y1 = find(x_d_M15(:)==x(1));
    y(1) = flux_d_M15(y1);
    y2 = find(x_d_M15(:)==x(2));
    y(2) = flux_d_M15(y2);
    flux_d_interp_M15(i) = interp1(x,y,x_t_M15(i));
end

figure()
plot(x_d_M15,flux_d_M15,x_t_M15,flux_t_M15,'linewidth',2)
title({"Part 2 - Final Project","Comparing Diffusion to Transport","M="+M})
ylabel('Flux $\phi$ ($n/(cm^2 s)$)','Interpreter','latex')
xlabel('Position $(cm)$','Interpreter','latex')
legend('diffusion','transport','location','best')
xlim([-(M+2/3) (M+2/3)])
ylim([0 2])
norm2_flux(M/5) = norm(abs(flux_d_M15-flux_t_M15),2)/norm(flux_t_M15,2);
fprintf('M15 norm: %4.2f\n',norm2_flux)
filename = ['Flux_comparison_M_' num2str(M) '.jpg'];
saveas(gcf,filename)

%%
M = 20;
flux_t_M20 = readmatrix('Part2_flux_from_transport_norm_convergence_M20.csv');
x_t_M20 = linspace(-M,M,100*M/5);
flux_d_M20 = readmatrix('Part2_flux_from_diffusion_norm2_convergence_M20.csv');
x_d_M20 = linspace(-(M+2/3),(M+2/3),100*M/5);
flux_d_interp_M20 = zeros(1,100*M/5);
for i = 1:length(x_t_M20)
    x = zeros(1,2);
    y = zeros(1,2);
    x(1) = max(x_d_M20(x_d_M20 < x_t_M20(i)));
    x(2) = min(x_d_M20(x_d_M20(:) > x_t_M20(i)));
    y1 = find(x_d_M20(:)==x(1));
    y(1) = flux_d_M20(y1);
    y2 = find(x_d_M20(:)==x(2));
    y(2) = flux_d_M20(y2);
    flux_d_interp_M20(i) = interp1(x,y,x_t_M20(i));
end

figure()
plot(x_d_M20,flux_d_M20,x_t_M20,flux_t_M20,'linewidth',2)
title({"Part 2 - Final Project","Comparing Diffusion to Transport","M="+M})
ylabel('Flux $\phi$ ($n/(cm^2 s)$)','Interpreter','latex')
xlabel('Position $(cm)$','Interpreter','latex')
legend('diffusion','transport','location','best')
xlim([-(M+2/3) (M+2/3)])
ylim([0 2])
norm2_flux(M/5) = norm(abs(flux_d_M20-flux_t_M20),2)/norm(flux_t_M20,2);
fprintf('M20 norm: %4.2f\n',norm2_flux)
filename = ['Flux_comparison_M_' num2str(M) '.jpg'];
saveas(gcf,filename)

%%
M = 25;
flux_t_M25 = readmatrix('Part2_flux_from_transport_norm_convergence_M25.csv');
x_t_M25 = linspace(-M,M,100*M/5);
flux_d_M25 = readmatrix('Part2_flux_from_diffusion_norm2_convergence_M25.csv');
x_d_M25 = linspace(-(M+2/3),(M+2/3),100*M/5);
flux_d_interp_M25 = zeros(1,100*M/5);
for i = 1:length(x_t_M25)
    x = zeros(1,2);
    y = zeros(1,2);
    x(1) = max(x_d_M25(x_d_M25 < x_t_M25(i)));
    x(2) = min(x_d_M25(x_d_M25(:) > x_t_M25(i)));
    y1 = find(x_d_M25(:)==x(1));
    y(1) = flux_d_M25(y1);
    y2 = find(x_d_M25(:)==x(2));
    y(2) = flux_d_M25(y2);
    flux_d_interp_M25(i) = interp1(x,y,x_t_M25(i));
end

figure()
plot(x_d_M25,flux_d_M25,x_t_M25,flux_t_M25,'linewidth',2)
title({"Part 2 - Final Project","Comparing Diffusion to Transport","M="+M})
ylabel('Flux $\phi$ ($n/(cm^2 s)$)','Interpreter','latex')
xlabel('Position $(cm)$','Interpreter','latex')
legend('diffusion','transport','location','best')
xlim([-(M+2/3) (M+2/3)])
ylim([0 2])
norm2_flux(M/5) = norm(abs(flux_d_M25-flux_t_M25),2)/norm(flux_t_M25,2);
fprintf('M25 norm: %4.2f\n',norm2_flux)
filename = ['Flux_comparison_M_' num2str(M) '.jpg'];
saveas(gcf,filename)

%%
M = 30;
flux_t_M30 = readmatrix('Part2_flux_from_transport_norm_convergence_M30.csv');
x_t_M30 = linspace(-M,M,100*M/5);
flux_d_M30 = readmatrix('Part2_flux_from_diffusion_norm2_convergence_M30.csv');
x_d_M30 = linspace(-(M+2/3),(M+2/3),100*M/5);
flux_d_interp_M30 = zeros(1,100*M/5);
for i = 1:length(x_t_M30)
    x = zeros(1,2);
    y = zeros(1,2);
    x(1) = max(x_d_M30(x_d_M30 < x_t_M30(i)));
    x(2) = min(x_d_M30(x_d_M30(:) > x_t_M30(i)));
    y1 = find(x_d_M30(:)==x(1));
    y(1) = flux_d_M30(y1);
    y2 = find(x_d_M30(:)==x(2));
    y(2) = flux_d_M30(y2);
    flux_d_interp_M30(i) = interp1(x,y,x_t_M30(i));
end

figure()
plot(x_d_M30,flux_d_M30,x_t_M30,flux_t_M30,'linewidth',2)
title({"Part 2 - Final Project","Comparing Diffusion to Transport","M="+M})
ylabel('Flux $\phi$ ($n/(cm^2 s)$)','Interpreter','latex')
xlabel('Position $(cm)$','Interpreter','latex')
legend('diffusion','transport','location','best')
xlim([-(M+2/3) (M+2/3)])
ylim([0 2])
norm2_flux(M/5) = norm(abs(flux_d_M30-flux_t_M30),2)/norm(flux_t_M30,2);
fprintf('M30 norm: %4.2f\n',norm2_flux)
filename = ['Flux_comparison_M_' num2str(M) '.jpg'];
saveas(gcf,filename)

%%
M = 35;
flux_t_M35 = readmatrix('Part2_flux_from_transport_norm_convergence_M35.csv');
x_t_M35 = linspace(-M,M,100*M/5);
flux_d_M35 = readmatrix('Part2_flux_from_diffusion_norm2_convergence_M35.csv');
x_d_M35 = linspace(-(M+2/3),(M+2/3),100*M/5);
flux_d_interp_M35 = zeros(1,100*M/5);
for i = 1:length(x_t_M35)
    x = zeros(1,2);
    y = zeros(1,2);
    x(1) = max(x_d_M35(x_d_M35 < x_t_M35(i)));
    x(2) = min(x_d_M35(x_d_M35(:) > x_t_M35(i)));
    y1 = find(x_d_M35(:)==x(1));
    y(1) = flux_d_M35(y1);
    y2 = find(x_d_M35(:)==x(2));
    y(2) = flux_d_M35(y2);
    flux_d_interp_M35(i) = interp1(x,y,x_t_M35(i));
end

figure()
plot(x_d_M35,flux_d_M35,x_t_M35,flux_t_M35,'linewidth',2)
title({"Part 2 - Final Project","Comparing Diffusion to Transport","M="+M})
ylabel('Flux $\phi$ ($n/(cm^2 s)$)','Interpreter','latex')
xlabel('Position $(cm)$','Interpreter','latex')
legend('diffusion','transport','location','best')
xlim([-(M+2/3) (M+2/3)])
ylim([0 2])
norm2_flux(M/5) = norm(abs(flux_d_M35-flux_t_M35),2)/norm(flux_t_M35,2);
fprintf('M35 norm: %4.2f\n',norm2_flux)
filename = ['Flux_comparison_M_' num2str(M) '.jpg'];
saveas(gcf,filename)

%%
M = 40;
flux_t_M40 = readmatrix('Part2_flux_from_transport_norm_convergence_M40.csv');
x_t_M40 = linspace(-M,M,100*M/5);
flux_d_M40 = readmatrix('Part2_flux_from_diffusion_norm2_convergence_M40.csv');
x_d_M40 = linspace(-(M+2/3),(M+2/3),100*M/5);
flux_d_interp_M40 = zeros(1,100*M/5);
for i = 1:length(x_t_M40)
    x = zeros(1,2);
    y = zeros(1,2);
    x(1) = max(x_d_M40(x_d_M40 < x_t_M40(i)));
    x(2) = min(x_d_M40(x_d_M40(:) > x_t_M40(i)));
    y1 = find(x_d_M40(:)==x(1));
    y(1) = flux_d_M40(y1);
    y2 = find(x_d_M40(:)==x(2));
    y(2) = flux_d_M40(y2);
    flux_d_interp_M40(i) = interp1(x,y,x_t_M40(i));
end

figure()
plot(x_d_M40,flux_d_M40,x_t_M40,flux_t_M40,'linewidth',2)
title({"Part 2 - Final Project","Comparing Diffusion to Transport","M="+M})
ylabel('Flux $\phi$ ($n/(cm^2 s)$)','Interpreter','latex')
xlabel('Position $(cm)$','Interpreter','latex')
legend('diffusion','transport','location','best')
xlim([-(M+2/3) (M+2/3)])
ylim([0 2])
norm2_flux(M/5) = norm(abs(flux_d_M40-flux_t_M40),2)/norm(flux_t_M40,2);
fprintf('M40 norm: %4.2f\n',norm2_flux)
filename = ['Flux_comparison_M_' num2str(M) '.jpg'];
saveas(gcf,filename)

%%
M = 45;
flux_t_M45 = readmatrix('Part2_flux_from_transport_norm_convergence_M45.csv');
x_t_M45 = linspace(-M,M,100*M/5);
flux_d_M45 = readmatrix('Part2_flux_from_diffusion_norm2_convergence_M45.csv');
x_d_M45 = linspace(-(M+2/3),(M+2/3),100*M/5);
flux_d_interp_M45 = zeros(1,100*M/5);
for i = 1:length(x_t_M45)
    x = zeros(1,2);
    y = zeros(1,2);
    x(1) = max(x_d_M45(x_d_M45 < x_t_M45(i)));
    x(2) = min(x_d_M45(x_d_M45(:) > x_t_M45(i)));
    y1 = find(x_d_M45(:)==x(1));
    y(1) = flux_d_M45(y1);
    y2 = find(x_d_M45(:)==x(2));
    y(2) = flux_d_M45(y2);
    flux_d_interp_M45(i) = interp1(x,y,x_t_M45(i));
end

figure()
plot(x_d_M45,flux_d_M45,x_t_M45,flux_t_M45,'linewidth',2)
title({"Part 2 - Final Project","Comparing Diffusion to Transport","M="+M})
ylabel('Flux $\phi$ ($n/(cm^2 s)$)','Interpreter','latex')
xlabel('Position $(cm)$','Interpreter','latex')
legend('diffusion','transport','location','best')
xlim([-(M+2/3) (M+2/3)])
ylim([0 2])
norm2_flux(M/5) = norm(abs(flux_d_M45-flux_t_M45),2)/norm(flux_t_M45,2);
fprintf('M45 norm: %4.2f\n',norm2_flux)
filename = ['Flux_comparison_M_' num2str(M) '.jpg'];
saveas(gcf,filename)

%%
M = 50;
flux_t_M50 = readmatrix('Part2_flux_from_transport_norm_convergence_M50.csv');
x_t_M50 = linspace(-M,M,100*M/5);
flux_d_M50 = readmatrix('Part2_flux_from_diffusion_norm2_convergence_M50.csv');
x_d_M50 = linspace(-(M+2/3),(M+2/3),100*M/5);
flux_d_interp_M50 = zeros(1,100*M/5);
for i = 1:length(x_t_M50)
    x = zeros(1,2);
    y = zeros(1,2);
    x(1) = max(x_d_M50(x_d_M50 < x_t_M50(i)));
    x(2) = min(x_d_M50(x_d_M50(:) > x_t_M50(i)));
    y1 = find(x_d_M50(:)==x(1));
    y(1) = flux_d_M50(y1);
    y2 = find(x_d_M50(:)==x(2));
    y(2) = flux_d_M50(y2);
    flux_d_interp_M50(i) = interp1(x,y,x_t_M50(i));
end

figure()
plot(x_d_M50,flux_d_M50,x_t_M50,flux_t_M50,'linewidth',2)
title({"Part 2 - Final Project","Comparing Diffusion to Transport","M="+M})
ylabel('Flux $\phi$ ($n/(cm^2 s)$)','Interpreter','latex')
xlabel('Position $(cm)$','Interpreter','latex')
legend('diffusion','transport','location','best')
xlim([-(M+2/3) (M+2/3)])
ylim([0 2])
norm2_flux(M/5) = norm(abs(flux_d_M50-flux_t_M50),2)/norm(flux_t_M50,2);
fprintf('M50 norm: %4.2f\n',norm2_flux)
filename = ['Flux_comparison_M_' num2str(M) '.jpg'];
saveas(gcf,filename)

%%
M = [5,10,15,20,25,30,35,40,45,50];
figure()
plot(M,norm2_flux,'linewidth',2)
xlabel('M')
ylabel('Euclidean distance')
title('Euclidean distance between the transport flux and diffusion flux')
filename = ['Euclidean difference.jpg'];
saveas(gcf,filename)