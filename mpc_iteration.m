function [u, x] = mpc_iteration(r, d, y, options)
    persistent sys L x0 u0 opt Psi Theta Ups Xi Hp Hu H Q R QP_opt mu1 mu2;
    if isempty(opt)
        if ~exist("options", "var")
            opt = struct("MPC_Enabled", true, "Input_Default", 0, ...
                         "Hp", 1000, "Hu", 20, ...
                         "mu1", 1, "mu2", 60, "mu3", 2, "mu4", 1);
        else
            opt = options;
        end
    end
    if isempty(sys)
        Hu = opt.Hu;
        Hp = opt.Hp;
        % Plant system without hot water production
        sys = get_plant("../../Fitting/models/hybrid_model_2seg_2_best.mat");
%         x0 = zeros(size(sys.A, 1), 1);
        x0 = [35; 30; 25; -2; 20];
        u0 = d(1) + opt.mu3;
        % Luenberger observer
        P = [0.01, 0.02, 0.03, 0.04, 0.05];
        L = place(sys.A', sys.Cy', P)';
        % Construct lifted input matrices
        [Psi, Theta, Ups, Xi] = get_lifted(sys, Hp, Hu);
        Q = eye(Hp);
        R = zeros(Hu);
        H = Theta' * Q * Theta + R;
        QP_opt =  optimoptions('quadprog', 'Display','off');
    end
    % Estimate states using original model (used for room temperature)
    y_est = sys.Cy * x0 + sys.Dy * u0 + sys.Ddy * d(1);
    x = sys.A * x0 + sys.B * u0 + sys.Bd * d(1) + L*(y - y_est);
    if opt.MPC_Enabled
        r = match_Hp(r, Hp);
        d = match_Hp(d, Hp);
        Eps = r - Psi * x - Ups * u0 - Xi * d;
        G = 2*Theta' * Q * Eps;
        [A_con, b_con] = get_lifted_constraints(Hu, u0, d);
        deltaU = quadprog(2*H, -G', A_con, b_con,[],[],[],[],[], QP_opt);
        u = u0 + deltaU(1);
    else
        u = opt.Input_Default;
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
        f1 = kron(ones(Hu, 1), [opt.mu1; -opt.mu2; opt.mu3]);
        f2 = kron(eye(Hu), [0; 0; 1]) * d(1:Hu);
        f = f1+f2;
        F = [F, f];
        % Input rate constrain E = [W w]
        W = [eye(Hu); -eye(Hu)];
        w = [ones(Hu, 1); ones(Hu, 1)] * opt.mu4;
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
    function [Psi, Theta, Ups, Xi] = get_lifted(sys, Hp, Hu)
        calA = cell([Hp, 1]);
        for i=1:Hp
            calA{i} = sys.A^i;
        end
        calB = cell([Hp, 1]);
        sum = zeros(size(sys.B));
        for i=1:Hp
            sum = sum + sys.A^(i-1)*sys.B;
            calB{i} = sum;
        end
        calBu = cell([Hp, Hp-Hu]);
        column = calB;
        for i=1:Hu
            calBu(:,i) = column;
            column{end} = zeros(size(sys.B));
            column = circshift(column, 1);
        end
        % Construct lifted output matrices
        calC = kron(eye(Hp), sys.Cz);
        Psi = calC * cell2mat(calA);
        Ups = calC * cell2mat(calB);
        Theta = calC * cell2mat(calBu);
        Xi = cell([Hp, Hp]);
        column = cell([Hp, 1]);
        for i=1:Hp
            column{i} = sys.Cz*sys.A^(i-1)*sys.Bd;
        end
        for i=1:Hp
            Xi(:,i) = column;
            column{end} = zeros(size(sys.Cz*sys.Bd));
            column = circshift(column, 1);
        end
        Xi = cell2mat(Xi);
    end
end