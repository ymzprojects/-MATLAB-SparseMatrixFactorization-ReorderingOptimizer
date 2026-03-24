%% Nettoyage de l'environnement
close all;
clear all;

%% Chargement de la matrice
% Pour tester, décommentez l'une des lignes ci-dessous selon la matrice choisie.
load A;   % La matrice A est chargée
% load pde225_5e-1; % La matrice A est chargée
% load piston;     % La matrice A est chargée
n = size(A,1);
b = (1:n)';   % Vecteur b = [1; 2; 3; ...; n]

%% 0. Résolution directe (sans factorisation explicite) pour référence
x_direct = A\b;
residual_direct = norm(b - A*x_direct) / (norm(full(A))*norm(x_direct) + norm(b));
fprintf('=== Solution directe (A\\b) ===\n');
fprintf('Residual error: %e\n\n', residual_direct);

%% 1. Factorisation LU sans permutation externe (référence)
[L, U, P] = lu(A);
x_lu = U \ (L \ (P*b));
residual_lu = norm(b - A*x_lu) / (norm(full(A))*norm(x_lu) + norm(b));

% Comptage des non-zéros et calcul des flops
nnz_L = nnz(L);
nnz_U = nnz(U);
row_counts_L = full(sum(spones(L), 2));
flops_forward = nnz(L) + 2 * sum(row_counts_L(2:end)) - 2*(n-1);
row_counts_U = full(sum(spones(U), 2));
flops_backward = nnz(U) + 2 * sum(row_counts_U(1:end-1)) - 2*(n-1);
flops_total = flops_forward + flops_backward;

fprintf('=== LU sans permutation externe ===\n');
fprintf('Nombre de non-zéros dans L : %d\n', nnz_L);
fprintf('Nombre de non-zéros dans U : %d\n', nnz_U);
fprintf('Opérations flottantes : %d\n\n', flops_total);

%% 2. Comparaison des méthodes de réordonnancement
methods = {'amd', 'colamd', 'symamd', 'symrcm', 'colperm'};
num_methods = length(methods);
results = struct('method',[],'nnz_L',[],'nnz_U',[],'flops',[],'residual',[],'error',[],...
                 'A_perm', [], 'L_perm', [], 'U_perm', []);

for i = 1:num_methods
    method = methods{i};
    try
        switch method
            case 'colamd', Q = colamd(A);
            case 'colperm', Q = colperm(A);
            case 'amd', Q = amd(A);
            case 'symamd', Q = symamd(A);
            case 'symrcm', Q = symrcm(A);
        end
        
        A_perm = A(:, Q);
        [L_perm, U_perm, P2] = lu(A_perm);
        y_perm = U_perm \ (L_perm \ (P2 * b));
        x_perm_full = zeros(n,1);
        x_perm_full(Q) = y_perm;
        residual_perm = norm(b - A*x_perm_full) / (norm(full(A))*norm(x_perm_full) + norm(b));
        
        nnz_L_c = nnz(L_perm);
        nnz_U_c = nnz(U_perm);
        row_counts_L_c = full(sum(spones(L_perm),2));
        flops_f = nnz(L_perm) + 2 * sum(row_counts_L_c(2:end)) - 2*(n-1);
        row_counts_U_c = full(sum(spones(U_perm),2));
        flops_b = nnz(U_perm) + 2 * sum(row_counts_U_c(1:end-1)) - 2*(n-1);
        
        results(i).method   = method;
        results(i).nnz_L    = nnz_L_c;
        results(i).nnz_U    = nnz_U_c;
        results(i).flops    = flops_f + flops_b;
        results(i).residual = residual_perm;
        results(i).error    = '';
        results(i).A_perm   = A_perm;
        results(i).L_perm   = L_perm;
        results(i).U_perm   = U_perm;
    catch ME
        results(i).error    = ME.message;
    end
end

% Affichage du tableau récapitulatif
T = table({results.method}', [results.nnz_L]', [results.nnz_U]', [results.flops]', [results.residual]', ...
    'VariableNames', {'Methode', 'nnz_L', 'nnz_U', 'Flops', 'Residual'});
disp(T);

%% 4. Visualisation globale du fill-in
nrows = 1 + num_methods;
ncols = 3; 
hFig = figure('Name', 'Analyse Fill-in LU', 'NumberTitle', 'off', 'Position', [100, 100, 900, 1200]);

% Ligne 1 : cas sans permutation
subplot(nrows, ncols, 1); spy(A); title('A originale');
subplot(nrows, ncols, 2); spy(L); title('L sans réord.');
subplot(nrows, ncols, 3); spy(U); title('U sans réord.');

% Lignes suivantes : méthodes de permutation
for i = 1:num_methods
    row = i + 1;
    if isempty(results(i).error)
        subplot(nrows, ncols, (row-1)*ncols + 1); spy(results(i).A_perm);
        title(sprintf('A (%s)', results(i).method));
        
        subplot(nrows, ncols, (row-1)*ncols + 2); spy(results(i).L_perm);
        title(sprintf('L (%s)', results(i).method));
        
        subplot(nrows, ncols, (row-1)*ncols + 3); spy(results(i).U_perm);
        title(sprintf('U (%s)', results(i).method));
    end
end

%% 5. Sauvegarde de la figure en PNG
% exportgraphics est préférable à saveas pour la qualité des graphiques de sparsité
filename = 'fillin_analysis_LU.png';
exportgraphics(hFig, filename, 'Resolution', 300);

fprintf('\nFigure sauvegardée avec succès : %s\n', filename);