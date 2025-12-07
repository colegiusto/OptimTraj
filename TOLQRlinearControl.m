function u = TOLQRlinearControl(t,x,c)

    % Initialize u
    u = 0;

    % Determine most recent knot point
    index = find(c.t_k > t, 1)-1;
    if t >= c.t_k(end)
        index = length(t_k);
    end

    % Interpolate trajectory reference state
    x_star = c.x_star(:,index-1)+(c.x_star(:,index) ...
        -c.x_star(:,index-1))-(t-t_k(index-1)) ...
        /(t_k(index)-t_k(index-1));

    % Interpolate trajectory reference input
    u_star = c.u_star(index-1)+(c.u_star(index)-u(index-1)) ...
        *(t-t_k(index-1))/(t_k(index)-t_k(index-1));

    % Pick K linearized around knot point
    K = c.K(:,:,index);

    % Calculate input so system converges to trajectory
    u(1) = u_star - K*(x-x_star);

end