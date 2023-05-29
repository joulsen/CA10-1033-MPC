addpath("..\..\commoncode\")
clear mpc_iteration
sys = get_plant("../../Fitting/models/hybrid_model_2seg_2_best.mat");

data = read_bitzer_data("../../BitzerData/Logdata-d-60/HPLog 026.csv", [1, -1]);
steps = size(data, 1);
u = zeros(steps, 1);
r = ones(steps, 1)*25;
d = data.Tamb;
x = zeros(steps, size(sys.A, 1));
xest = zeros(steps, size(sys.A, 1) + 1);
y_C = zeros(steps, size(sys.Cy, 1));
wait = waitbar(0, "Computing MPC");
tic
for i=2:steps
    ri = r(i:end,:);
    di = d(i:end,:);
    if i == 3116
        disp("here");
    end
    [u(i,:), xest(i,:)] = mpc_iteration(ri, di, y_C(i-1,:));
    x(i,:) = sys.A * x(i-1,:)' + sys.B * u(i-1,:)' + sys.Bd * d(i-1,:)';
    y_C(i,:) = sys.Cy * x(i,:)' + sys.Dy * u(i,:)' + sys.Ddy * d(i-1,:)';
    if mod(i, 100) == 0
        waitbar(i / steps, wait, "Computing MPC")
    end
end
toc
close(wait)

