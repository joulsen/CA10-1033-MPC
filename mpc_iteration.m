function deltaU = mpc_iteration(x, u, r, d, Hp_bar, Hu_bar, param, options)
    %MPC_ITERATION_BASE Performs a single iteration of the MPC
    if options.MPC_Enabled
        % Cropping to Hp_bar
        dim = struct("p", size(param.msys.Cz, 1), ...
                     "m", size(param.msys.A, 1), ...
                     "n", size(param.msys.B, 2), ...
                     "j", size(param.msys.Bd, 2));
        Psi = param.Psi(1:dim.p*Hp_bar, :);
        Ups = param.Ups(1:dim.p*Hp_bar, :);
        Theta = param.Theta(1:dim.p*Hp_bar, 1:dim.n*Hu_bar);
        Xi = param.Xi(1:dim.p*Hp_bar, 1:dim.j*Hp_bar);
        Q = param.Q(1:Hp_bar, 1:Hp_bar);
        R = param.R(1:Hu_bar, 1:Hu_bar);
        H = Theta' * Q * Theta + R;
        H = (H + H')/2;
        Eps = r - Psi * x - Ups * u - Xi * d;
        G = 2*Theta' * Q * Eps;
        [A_con, b_con] = get_lifted_constraints(Hu_bar, u, d);
        qp_opt = optimoptions('quadprog', 'Algorithm', 'active-set', ...
                              'Display', 'off');
        deltaU = quadprog(2*H, -G', A_con, b_con, ...
                          [],[],[],[],zeros(Hu_bar, 1),qp_opt);
    else
        deltaU = 0;
    end
    function [A_con, b_con] = get_lifted_constraints(Hu, u0, d)
        nu = 1;
        nf = 3;
        % Input constraint F
        F = kron(eye(Hu, Hu), [-1; 1; -1]);
        f1 = kron(ones(Hu, 1), [param.mu1; -param.mu2; param.mu3]);
        f2 = kron(eye(Hu), [0; 0; 1]) * d(1:Hu);
        f = f1+f2;
        F = [F, f];
        % Input rate constrain E = [W w]
        W = [eye(Hu); -eye(Hu)];
        w = [ones(Hu, 1); ones(Hu, 1)] * param.mu4 * param.kappa;
        % Lifting input constraint F
        calF = zeros(Hu*nf, Hu*nu);
        for i=1:Hu
            calFi = zeros(Hu*nf, nu);
            for j=i:Hu
                Fi = F(:,nu*j);
                calFi = calFi + Fi;
            end
            calF(:,i) = calFi;
        end
        A_con = [calF; W];
        b_con = [-calF(:,1:nu) * u0 - f; w];
    end
end