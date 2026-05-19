%Alignment checking
fTt1 = tr.Transform.fTt(:,:,1);

% Pull first-row flexion (robust to table/array + tibiofemoral/tibiafemoral naming)
K = tr.Kinematics;
if isfield(K,'tibiofemoral') && ~isempty(K.tibiofemoral)
    TF = K.tibiofemoral;
elseif isfield(K,'tibiafemoral') && ~isempty(K.tibiafemoral)
    TF = K.tibiafemoral;
else
    error('No tibiofemoral kinematics found');
end

if istable(TF)
    flex0 = TF{1,'flexion'};
else
    flex0 = TF(1,1); % assumes flexion is column 1 (matches your printout)
end

disp('fTt(:,:,1) = ');
disp(fTt1);
fprintf('flexion(1) = %.6f (deg)\n', flex0);
disp('===================================');
