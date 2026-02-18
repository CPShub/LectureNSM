function plot_result(x_axis,Q,M,w,w_analytical)
    subplot(3,1,1)
    plot(x_axis,Q)
    title('Querkraft')
    ylabel('Q(x)')
    xlabel('x')
    grid on
    
    subplot(3,1,2)
    plot(x_axis,M)
    title('Biegemoment')
    ylabel('M(x)')
    xlabel('x')
    grid on
    
    subplot(3,1,3)
    plot(x_axis,w)
    hold on
    plot(x_axis,w_analytical(x_axis))
    hold off
    set(gca, 'YDir','reverse')
    title('Durchbiegung')
    ylabel('w(x)')
    xlabel('x')
    grid on
end

