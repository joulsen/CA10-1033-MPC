function param = mpc_iteration_initialize(param)
    % MPC_ITERATION_INITIALIZE Constructs the param struct for use in
    % mpc_iteration_base.m
    % Plant system without hot water production
    param.sys = get_plant("../../Fitting/models/hybrid_model_2seg_2_best.mat", 1);
    param.msys = get_plant("../../Fitting/models/hybrid_model_2seg_2_best.mat", param.kappa);
%     deleteproperty(sys)
    % Luenberger observer
    param.L = place(param.sys.A', param.sys.Cy', param.P)';
    % Resample 
    % Construct lifted input matrices
    [param.Psi, param.Theta, param.Ups, param.Xi] = get_lifted(param.msys, param.Hp, param.Hu);
    param.H = param.Theta' * param.Q * param.Theta + param.R;
    param.H = (param.H + param.H')/2;
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

