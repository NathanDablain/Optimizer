function compare_solvers()
    % Define the system size (keep small as the textbook solver is slow)
    n = 250;

    % Generate a random dense matrix and right-hand side vector
    rng(42); % Seed for reproducibility
    A = rand(n, n) + n * eye(n); % Diagonally dominant to ensure stability
    b = rand(n, 1);

    % ==========================================
    % 1. Textbook Gaussian Elimination Solver
    % ==========================================
    tic;
    x_textbook = textbook_gaussian_elimination(A, b);
    time_textbook = toc;

    tic;
    x_modern = modern_direct_solver(A, b);
    time_modern = toc;

    % ==========================================
    % 2. Modern Direct Solver (MATLAB Built-in)
    % ==========================================
    tic;
    x_backslash = A \ b;
    time_backslash = toc;

    % ==========================================
    % Verification and Performance Display
    % ==========================================
    % Compute the norm of the difference to ensure both got the same answer
    accuracy_check = norm(x_textbook - x_modern);

    fprintf('--- Performance Results (Size: %d x %d) ---\n', n, n);
    fprintf('Textbook Solver Time : %.6f seconds\n', time_textbook);
    fprintf('Modern Direct Solver : %.6f seconds\n', time_modern);
    fprintf('Backslash Direct Solver : %.6f seconds\n', time_backslash);
    fprintf('Speedup Factor       : %.1fx faster\n', time_textbook / time_modern);
    fprintf('Solution Difference  : %e\n', accuracy_check);
end

%% --- THE HIGH-PERFORMANCE MODERN SOLVER ---
function x = modern_direct_solver(A, b)
    % Solves Ax = b by implementing a vectorized LU factorization
    % with partial pivoting, followed by vectorized forward/back substitution.
    n = size(A, 1);
    piv = (1:n)'; % Vector to track row permutations (P matrix)

    % --- Phase 1: Vectorized LU Factorization (PA = LU) ---
    for i = 1:n-1
        % 1. Vectorized Pivot Selection: Find max element in column i (from row i down)
        [~, max_idx] = max(abs(A(i:end, i)));
        pivot_row = max_idx + i - 1; % Adjust index relative to the full matrix

        % 2. Vectorized Row Swaps (Permutation)
        if pivot_row ~= i
            % Swap rows in A
            temp_row = A(i, :);
            A(i, :) = A(pivot_row, :);
            A(pivot_row, :) = temp_row;

            % Record the permutation in the pivoting vector
            temp_piv = piv(i);
            piv(i) = piv(pivot_row);
            piv(pivot_row) = temp_piv;
        end

        % 3. Vectorized Elimination (No Loops)
        % Compute multipliers for the column all at once
        A(i+1:end, i) = A(i+1:end, i) / A(i, i);

        % Rank-1 outer product update (Processes the entire remaining submatrix at once)
        % This line utilizes BLAS-like vectorized memory handling
        A(i+1:end, i+1:end) = A(i+1:end, i+1:end) - A(i+1:end, i) * A(i, i+1:end);
    end

    % Apply the permutation vector to the right-hand side vector b (Pb)
    b = b(piv);

    % --- Phase 2: Vectorized Forward Substitution (Ly = Pb) ---
    % Since L has 1s on the diagonal, we compute row-by-row but vector-multiply columns
    y = zeros(n, 1);
    for i = 1:n
        y(i) = b(i) - A(i, 1:i-1) * y(1:i-1);
    end

    % --- Phase 3: Vectorized Backward Substitution (Ux = y) ---
    x = zeros(n, 1);
    for i = n:-1:1
        x(i) = (y(i) - A(i, i+1:end) * x(i+1:end)) / A(i, i);
    end
end

%% --- Helper Function: Textbook Gaussian Elimination ---
function x = textbook_gaussian_elimination(A, b)
    n = length(b);
    % Create augmented matrix [A | b]
    M = [A, b];

    % Forward Elimination Phase
    for i = 1:n
        % 1. Partial Pivoting: Find the largest element in the current column
        pivot_row = i;
        for r = i+1:n
            if abs(M(r,i)) > abs(M(pivot_row,i))
                pivot_row = r;
            end
        end

        % Swap current row with pivot row if necessary
        if pivot_row ~= i
            temp = M(i,:);
            M(i,:) = M(pivot_row,:);
            M(pivot_row,:) = temp;
        end

        % 2. Eliminate entries below the pivot using loops and checks
        for r = i+1:n
            % Explicitly check if the entry is already zero to "skip" math
            if M(r,i) ~= 0
                factor = M(r,i) / M(i,i);
                for c = i:n+1
                    M(r,c) = M(r,c) - factor * M(i,c);
                end
            end
        end
    end

    % Back Substitution Phase
    x = zeros(n, 1);
    for i = n:-1:1
        sum_terms = 0;
        for j = i+1:n
            sum_terms = sum_terms + M(i,j) * x(j);
        end
        x(i) = (M(i,n+1) - sum_terms) / M(i,i);
    end
end
