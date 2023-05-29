function [u, x] = mpc_iteration_base(r, d, y, param, options)
    %MPC_ITERATION_BASE Performs a single iteration of the MPC
    persistent sys x0 u0 Hu Hp;
    if isempty(sys)
        sys = param.sys;
        x0 = param.x0;
        u0 = param.u0;
        Hu = param.Hu;
        Hp = param.Hp;
    end
    y_est = sys.Cy * x0 + sys.Dy * u0 + sys.Ddy * d(1);
    x = sys.A * x0 + sys.B * u0 + sys.Bd * d(1) + param.L*(y - y_est);
    if options.MPC_Enabled
        r = match_Hp(r, Hp);
        d = match_Hp(d, Hp);
        Eps = r - param.Psi * x - param.Ups * u0 - param.Xi * d;
        G = 2*param.Theta' * param.Q * Eps;
        [A_con, b_con] = get_lifted_constraints(Hu, u0, d);
        qp_opt = optimoptions('quadprog', 'Algorithm', 'active-set', ...
                              'Display', 'off');
        deltaU = quadprog(2*param.H, -G', A_con, b_con, ...
                          [],[],[],[],zeros(param.Hu, 1),qp_opt);
        u = u0 + deltaU(1);
    else
        u = options.Input_Default;
    end
    x0 = x;
    u0 = u;
    x = [x; u0(1) - x(3)]';
    u = u';
    function mat = match_Hp(mat, Hp)
        size_diff = size(mat, 1) - Hp;
        if size_diff < 0
            mat = [mat; kron(ones(abs(size_diff), 1), mat(end,:))];
        elseif size_diff > 0
            mat = mat(1:Hp);
        end
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
        w = [ones(Hu, 1); ones(Hu, 1)] * param.mu4;
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

