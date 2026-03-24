%% Nettoyage de l'environnement
close all;
clear all;

%% Chargement de la matrice SPD
% La matrice A doit être symétrique définie positive (exemple : mat3)
load mat3;  % La matrice A est chargée
n = size(A,1);
b = (1:n)';   % Vecteur b = [1; 2; ...; n]

%% Vérification que A est SPD
try 
    L_test = chol(A, 'lower');
catch ME
    error('La matrice A n''est pas symétrique définie positive.');
end

%% 0. Résolution directe (sans factorisation explicite) pour référence
x_direct = A\b;
residual_direct = norm(b - A*x_direct) / (norm(full(A))*norm(x_direct) + norm(b));
fprintf('=== Solution directe (A\\b) ===\n');
fprintf('Residual error: %e\n\n', residual_direct);

%% 1. Factorisation Cholesky sans réordonnancement (référence)
L_orig = chol(A, 'lower');  % Facteur de Cholesky (inférieur)
% Résolution par substitution (forward puis backward)
y_orig = L_orig \ b;
x_orig = L_orig' \ y_orig;
residual_orig = norm(b - A*x_orig) / (norm(full(A))*norm(x_orig) + norm(b));

% Calcul du nombre de flops pour la résolution triangulaire :
row_counts_orig = full(sum(spones(L_orig), 2));
flops_forward_orig = nnz(L_orig) + 2 * sum(row_counts_orig(2:end)) - 2*(n-1); 
flops_chol_orig = 2 * flops_forward_orig; 
nnz_L_orig = nnz(L_orig);

fprintf('=== Cholesky sans réordonnancement ===\n');
fprintf('Nombre de non-zéros dans L : %d\n', nnz_L_orig);
fprintf('Opérations flottantes (forward+backward) : %d\n', flops_chol_orig);
fprintf('Residual error: %e\n\n', residual_orig);

%% 2. Comparaison de plusieurs méthodes de réordonnancement
methods = {'amd', 'symamd', 'symrcm'};
num_methods = length(methods);
results = struct('method',[],'nnz_L',[],'flops',[],'residual',[],'error',[],...
                 'A_perm', [], 'L_perm', []);

for i = 1:num_methods
    method = methods{i};
    try
        switch method
            case 'amd'
                Q = amd(A);
            case 'symamd'
                Q = symamd(A);
            case 'symrcm'
                Q = symrcm(A);
        end
        
        A_perm = A(Q, Q);
        b_perm = b(Q);
        L_perm = chol(A_perm, 'lower');
        
        y_perm = L_perm \ b_perm;
        x_perm = L_perm' \ y_perm;
        x_perm_full = zeros(n,1);
        x_perm_full(Q) = x_perm;
        residual_perm = norm(b - A*x_perm_full) / (norm(full(A))*norm(x_perm_full) + norm(b));
        
        row_counts_perm = full(sum(spones(L_perm), 2));
        flops_forward_perm = nnz(L_perm) + 2 * sum(row_counts_perm(2:end)) - 2*(n-1);
        flops_chol_perm = 2 * flops_forward_perm;
        
        results(i).method   = method;
        results(i).nnz_L    = nnz(L_perm);
        results(i).flops    = flops_chol_perm;
        results(i).residual = residual_perm;
        results(i).error    = '';
        results(i).A_perm   = A_perm;
        results(i).L_perm   = L_perm;
    catch ME
        results(i).error    = ME.message;
    end
end

%% 3. Comparaison et affichage
best_flops = Inf;
best_method = '';
for i = 1:num_methods
    if isempty(results(i).error) && results(i).flops < best_flops
         best_flops = results(i).flops;
         best_method = results(i).method;
    end
end

% Affichage du tableau
T = table({results.method}', [results.nnz_L]', [results.flops]', [results.residual]', ...
    'VariableNames', {'Methode', 'nnz_L', 'Flops', 'Residual'});
disp(T);

%% 4. Visualisation globale des structures
nrows = 1 + num_methods;
ncols = 3;

% On stocke le handle de la figure pour la sauvegarde
hFig = figure('Name', 'Comparaison Reordonnancement', 'NumberTitle', 'off');

% Ligne 1 : cas sans permutation
subplot(nrows, ncols, 1); spy(A); title('A originale');
subplot(nrows, ncols, 2); spy(L_orig); title('Cholesky L (original)');
subplot(nrows, ncols, 3); spy(spones(L_orig+L_orig')); title('Sparsité L+L'' (original)');

% Lignes suivantes
for i = 1:num_methods
    row = i + 1;
    if isempty(results(i).error)
        subplot(nrows, ncols, (row-1)*ncols + 1); spy(results(i).A_perm);
        title(sprintf('A permutée (%s)', results(i).method));
        
        subplot(nrows, ncols, (row-1)*ncols + 2); spy(results(i).L_perm);
        title(sprintf('L (%s)', results(i).method));
        
        subplot(nrows, ncols, (row-1)*ncols + 3); spy(spones(results(i).L_perm + results(i).L_perm'));
        title(sprintf('Sparsité L+L'' (%s)', results(i).method));
    end
end

%% 5. Sauvegarde de la figure en PNG
% exportgraphics est idéal car il recadre automatiquement l'image
filename = 'comparaison_structures.png';
exportgraphics(hFig, filename, 'Resolution', 300);

fprintf('\nFigure sauvegardée avec succès sous le nom : %s\n', filename);