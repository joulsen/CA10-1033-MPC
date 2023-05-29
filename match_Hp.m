function [mat, Hp_bar] = match_Hp(mat, Hp, kappa)
    v = size(mat, 1);
    Hp_bar = min(Hp, floor(v / kappa));
    mat = mat(1:kappa:Hp_bar*kappa);
end
