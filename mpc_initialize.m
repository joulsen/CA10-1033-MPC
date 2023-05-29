function param = mpc_initialize(param)
    % MPC_ITERATION_INITIALIZE Constructs the param struct for use in
    % mpc_iteration_base.m
    % Plant system without hot water production
    param.sys = get_plant("model.mat", 1);
    param.psys = get_plant("model.mat", 1);
%     param.psys.A(2,3) = param.psys.A(2,3) - 1e-4;
%     param.sys.A(2,3) = param.sys.A(2,3) - 1e-4;
    param.msys = get_plant("model.mat", param.kappa);
    % Luenberger observer
    param.L = place(param.sys.A', param.sys.Cy', param.P)';
    % Resample 
    % Construct lifted input matrices
    [param.Psi, param.Theta, param.Ups, ...
     param.Xi, param.Lambda, param.LambdaD] = get_lifted(param.msys, param.Hp, param.Hu);
    function [Psi, Theta, Ups, Xi, Lambda, LambdaD] = get_lifted(sys, Hp, Hu)
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
        Lambda = kron(ones(Hp, 1), sys.Dz);
        LambdaD = kron(tril(ones(Hp, Hu)), sys.Dz);
    end
end

