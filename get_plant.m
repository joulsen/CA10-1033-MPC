function sys_plant = get_plant(fname, kappa)
%GET_PLANT Constructs the plant model using the hybrid_2seg_2_best model
    sys = load(fname).sys;
    if kappa > 1
        sys = d2d(sys, sys.Ts * kappa, "zoh");
    end
%     subplot(1,2,1);
%     [p1, z1] = pzmap(sys)
%     pzmap(sys)
%     sys = d2d(sys, sys.Ts * 20, "zoh");
%     subplot(1,2,2)
%     [p2, z2] = pzmap(sys)
%     pzmap(sys)
    % Initial variables
    nx = size(sys.A, 1) - 1;
    nu = size(sys.B, 2) - 2;
    ny = size(sys.C, 1);
    % Plant system withsys_plant hot water production
%     sys_plant = ss(sys.A(1:nx, 1:nx), sys.B(1:nx, 1:nu), ...
%                    sys.C(1:ny, 1:nx), sys.D(1:ny, 1:nu), ...
%                    sys.Ts);
    sys_plant = struct("A", sys.A(1:nx, 1:nx), ...
                       "B", sys.B(1:nx, 1), "Bd", sys.B(1:nx, 2), ...
                       "Cy", sys.C(1,1:nx), "Cz", sys.C(2,1:nx), ...
                       "Dy", sys.D(1,1), "Ddy", sys.D(1,2), ...
                       "Dz", sys.D(2:ny,1), "Ddz", sys.D(2:ny,2), ...
                       "Ts", sys.Ts, ...
                       "nx", nx, "nu", nu, "ny", ny, ...
                       "hotel", "Trivago");
%     sys_plant.InputName = ["Tin", "Tamb"];
%     sys_plant.OutputName = ["T_w_r", "T_r"];
%     sys_plant.StateName = ["T_w1", "T_w2", "T_w3", "T_f", "T_r"];
end
