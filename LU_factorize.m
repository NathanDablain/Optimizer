function [A, A_rows, A_cols, L_rows, L_cols, U_rows, U_cols, p_perm, q_perm] = LU_factorize(cons)
  H_contribution = (cons.H_block * cons.N_knots) + cons.N_sv;
  J_contribution = 2*(cons.J_block * (cons.N_knots - 1) + cons.N_ic + 2*cons.N_sv);
  A_nz = H_contribution + J_contribution;

  N = cons.N_con;
  M = cons.N_dec;
  % The hessian will always be diagonal, the nonlinear defects sit in blocks on the diagonal
  % The quadratic final state from the cost function will sit on the diagonal
  % The slack vars sit on the diagonal
  H_test = evaluate_hessian(zeros(M,1), zeros(N,1), cons, true);
  J_test = evaluate_jacobian(zeros(M,1), cons, true);

  A = [H_test J_test';...
       J_test zeros(N,N)];

  [A_rows, A_cols] = find(A);

  q_perm = colamd(A);
  A_colperm = A(:, q_perm);

  [L_struct, U_struct, p_perm] = lu(A_colperm);

  [L_rows, L_cols] = find(L_struct);
  [U_rows, U_cols] = find(U_struct);

  %% Begin copied code
  function x = my_phase2_collocation_solver(new_A_vals, rows_A, cols_A, ...
                L_rows, L_cols, num_L_nz, U_rows, U_cols, num_U_nz, ...
                b, p_perm, q_perm, matrix_dim)

    % --- 1. PRESERVE THE LAYOUT & INJECT VALUES ---
    % Reconstruct A using the exact static mapping to preserve structure
    A_updated = sparse(rows_A, cols_A, new_A_vals, matrix_dim, matrix_dim);

    % Apply the pre-calculated static column permutation
    A_colperm = A_updated(:, q_perm);

    % --- 2. NUMERICAL FACTORIZATION (Bypassing Reordering/Symbolic Steps) ---
    % Allocate static memory arrays
    L_vals = zeros(num_L_nz, 1);
    U_vals = zeros(num_U_nz, 1);

    % loop through each column of the collocation grid
    for j = 1:matrix_dim
        % Dense workspace vector for column j (Scatter step)
        w = zeros(matrix_dim, 1);

        % Gather A's column values into the dense workspace based on the static row mapping
        % (Applying the pre-calculated row permutation 'p_perm')
        A_col_idx = find(A_colperm.cols == j); % Conceptually accessing column j
        % In MATLAB matrix terms, we map the row-permuted column:
        w = A_colperm(p_perm, j);

        % LEFT-LOOKING UPDATE: Eliminate entries using prior columns
        for k = 1:(j-1)
            if w(k) ~= 0
                % Find structural non-zeros in column k of L below the diagonal
                l_idx = find(L_cols == k & L_rows > k);
                for idx = l_idx'
                    row = L_rows(idx);
                    w(row) = w(row) - L_vals(idx) * w(k);
                end
            end
        end

        % STATIC PIVOT GUARD (Crucial for Trapezoidal KKT Saddle-Points)
        % If a diagonal pivot drops near zero due to an indefinite step,
        % force it to a small threshold value to preserve the static structure.
        diag_val = w(j);
        if abs(diag_val) < 1e-12
            diag_val = sign(diag_val) * 1e-12 + (diag_val == 0) * 1e-12;
        end

        % Store values into the static U structure (Upper triangle)
        u_idx = find(U_cols == j);
        for idx = u_idx'
            U_vals(idx) = w(U_rows(idx));
        end

        % Store values into the static L structure (Lower triangle, unit diagonal)
        l_idx = find(L_cols == j);
        for idx = l_idx'
            row = L_rows(idx);
            if row == j
                L_vals(idx) = 1.0;
            else
                L_vals(idx) = w(row) / diag_val;
            end
        end
    end

    % --- 3. FORWARD SUBSTITUTION (L * y = P * b) ---
    b_perm = b(p_perm);
    y = zeros(matrix_dim, 1);
    for i = 1:matrix_dim
        idx_row_i = find(L_rows == i & L_cols < i);
        sum_ly = 0;
        for idx = idx_row_i'
            sum_ly = sum_ly + L_vals(idx) * y(L_cols(idx));
        end
        y(i) = b_perm(i) - sum_ly; % L has a unit diagonal of 1.0
    end

    % --- 4. BACKWARD SUBSTITUTION (U * x_perm = y) ---
    x_perm = zeros(matrix_dim, 1);
    for i = matrix_dim:-1:1
        idx_row_i = find(U_rows == i & U_cols > i);
        sum_ux = 0;
        for idx = idx_row_i'
            sum_ux = sum_ux + U_vals(idx) * x_perm(U_cols(idx));
        end
        diag_idx = find(U_rows == i & U_cols == i);
        x_perm(i) = (y(i) - sum_ux) / U_vals(diag_idx);
    end

    % --- 5. INVERSE PERMUTATION ---
    % Revert the column permutation to get back original variable layout
    x = zeros(matrix_dim, 1);
    x(q_perm) = x_perm;
end
  %% end copied code
##  num_L_nz = length(L_rows);
##  num_U_nz = length(U_rows);
end
