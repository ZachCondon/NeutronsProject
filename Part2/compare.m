i = 1;
for M = 5:5:50
    figure()
    hold on
    x = -(M-M/100):M/100:M;
    plot(x,Part2fluxfromtransport1(i,:))
    plot(xfromdiffusion1(i,:),Part2fluxfromdiffusion1(i,:))
    title('M = %d',M)
    legend('transport','diffusion')
    i = i + 1;
    hold off
end
