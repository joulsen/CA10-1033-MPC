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
        Lambda = param.Lambda(1:dim.p*Hp_bar, :);
        LambdaD = param.LambdaD(1:dim.p*Hp_bar, 1:dim.n*Hu_bar);
        TE = Ups + Lambda;
        TL = Theta + LambdaD;
        Q = param.Q(1:dim.p*Hp_bar, 1:dim.p*Hp_bar);
        R = param.R(1:Hu_bar, 1:Hu_bar);
        H = TL' * Q * TL + R;
        H = (H + H')/2;
        Z = Psi * x + TE * u + Xi * d;
        Eps = r - Z;
        G = 2*TL' * Q * Eps;
        rm = movmean(r(1:2:end), [param.mu5, 0]);
        [A_con, b_con] = get_lifted_constraints(Hu_bar, Hp_bar, u, d, rm, Z, TL);
        [deltaU, ~, exitflag] = quadprog(2*H, -G', A_con, b_con, ...
                                [],[],[],[],zeros(Hu_bar, 1),options.qp_opt);
        if exitflag ~= 1
            if exitflag == -2
                warning("Constraints have been violated. Default input is set")
            else
                warning("Quadratic problem does not produce solution. Default input is set")
            end
            deltaU = param.u_fallback - u;
        end
    else
        deltaU = 0;
    end
    function [A_con, b_con] = get_lifted_constraints(Hu, Hp, u0, d, rm, Z, TL)
        nu = 1;
        nf = 3;
        % Input constraint F
        F = kron(eye(Hu, Hu), [-1; 1; -1]);
        f1 = kron(ones(Hu, 1), [param.mu1; -param.mu2; param.mu3]);
        f2 = kron(eye(Hu), [0; 0; 1]) * d(1:Hu);
        f = f1+f2;
        F = [F, f];
        % Output constraint G = [Gf, gff]
        gff = nan(2*size(rm, 1), 1);
        gff(1:2:end) = rm;
        gff(2:2:end) = -rm;
        gff = gff - param.mu6;
        Gf = kron(eye(Hp), [-1, 0; 1, 0]);
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
        A_con = [calF; Gf*TL; W];
        b_con = [-calF(:,1:nu) * u0 - f; -Gf*Z-gff; w];
%         A_con = [calF; W];
%         b_con = [-calF(:,1:nu) * u0 - f; w];
    end
end