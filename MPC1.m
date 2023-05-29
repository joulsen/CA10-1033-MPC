addpath("../../commoncode/")

sys = load("../../Fitting/models/hybrid_model_2seg_2_best.mat").sys;
x0 = sys.x0;
A = sys.A(1:5, 1:5);
B = sys.B(1:5, 1:2);
C = sys.C(:, 1:5);
D = sys.C(:, 1:2);
sys = ss(A, B, C, D, sys.Ts);

sys.InputName = {'ST_in', 'T_a'};
sys.OutputName = {'T_w,r', 'T_r'};
sys.StateName = {'T_w1', 'T_w2', 'T_w3', 'T_f', 'T_r'};

sys = setmpcsignals(sys,"ManipulatedVariables", 1, ...
                        "MeasuredDisturbances", 2, ...
                        "MeasuredOutputs", 1, ...
                        "UnmeasuredOutputs", 2);

controller_ts = sys.Ts*10;
horz_p = 100;
horz_c = 10;
mObj = struct("Units", "°C", "Min", 20, "Max", 60);
oObj = [struct("Units", "°C"),
        struct("Units", "°C")];
dObj = struct("Units", "°C");
wObj = struct("ManipulatedVariables", 0, ...
              "ManipulatedVariablesRate", 0.1, ...
              "OutputVariables", [0, 1]);

simu_time = 5e5;

mpcobj = mpc(sys, controller_ts, horz_p, horz_c, wObj, mObj, oObj, dObj);

data = read_bitzer_data("../../BitzerData/LogData-d-60/HPLog A023.csv", [1, -1]);
simu_time = 0:controller_ts:simu_time;
Tamb_time = (0:size(data, 1)-1)*60;
Tamb = interp1(Tamb_time, data.Tamb, simu_time);

sim(mpcobj, 1000, [0 0; 0, 25], Tamb');

% fig = default_fig();
% 
% for i=1:10
%     mpcobj = mpc(sys, sys.Ts*2^(i-1));
%     mpcobj.PredictionHorizon = 1000;
%     mpcobj.ControlHorizon = 10;
%     mpcobj.Weights.OutputVariables = [0, 1]; % Only weight T_r
%     mpcobj.MV.min = 0;
%     mpcobj.MV.max = 60;
%     % review(mpcobj)
%     opt = mpcsimopt(mpcobj);
%     sys2 = sys;
%     sys2.A(1,1) = -0.0106;
%     opt.Model = sys2;
%     [y, t, u, xp, xmpc, opt] = sim(mpcobj, 1000, [0 20; 0 25], 20, opt);
%     subplot(3,1,1)
%     hold on
%     plot(t,y(:,1))
%     subplot(3,1,2)
%     hold on
%     plot(t,y(:,2))
%     subplot(3,1,3)
%     hold on
%     plot(t,u(:,1))
% end
% 
% for i=1:3
%     subplot(3,1,i)
%     xlim([0, 10e5])
%     if i == 2
%         legend(string(2.^(1:10)), "Location","bestoutside")
%     end
% end