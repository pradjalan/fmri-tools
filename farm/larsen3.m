% Change Log
% 1. The max is to be take over inactive set. I inserted in the following
%    line: temp = [(c(I) - cmax) ./ (cd(I) - cmax); (c(I) + cmax) ./ (cd(I) + cmax)];
%
% 2. Exact comutation of G using G_array = X(:,A)' * (y - mu) ./ sb / a2;
%
% 3. Sanity check for gradients using if ((max(G_array) - min(G_array)) / max(G_array) > 1e-5)
%
% 4. If LassoCond becomes true, then active set should be updated only if we
%    go into the next iteration, else it will interfere with the trimming at the end of the interation.
%    The active set update moved to the beginning of the loop.
% 
% 5. Final bug: In the last step while trimming for exact g value, the sg
%    computed is incorrect if the lassoCond became true in the last
%    iteration. sg computation should take care of lassoCond and get
%    previous sign for the beta (b) which has just become zero. Fixed this
%    and the code now works correctly. 
%
% 6. Another bug found. In the following for computing gamma_tilde replaced end-1 with end. The correct one is given below
%    gamma_tilde = b(A(1: end)) ./ (b(A(1: end)) - b_OLS(1: end, 1));
%
% Usage:  [b, g_out, steps] = larsen(flat_mri, model_index, g, max_l1, min_l2, max_ss, max_steps)
%
% 7. Another bug (actually not a bug but improvement) in the original lasso algorithm: If the matrix of active set becomes close to singular
%    then the algorithm becomes numerically unstable and sometimes give incorrect results. Fixed this by adding the following condition  
%    if (cond(Gram) > 1e2) then skip the current variable and remove it from inactive set. This is a topic for further research.
%    Now the algorithm is robust and gives good results and is also fast.
%
%


function [b, g_out, steps] = larsen(flat_mri, model_index, g, max_l1, min_l2, max_ss, max_steps)


X = double(flat_mri);
X(:, model_index) = [];
X = X(1: end - 1, :);
y = flat_mri(2: end, model_index);
X = zscore(X) / sqrt(size(X, 1) - 1);
y = zscore(y) / sqrt(size(y, 1) - 1);
% X = zscore(X);
% y = zscore(y);
[n p] = size(X);

if (g <= 0) g = 0.001; end;
if (max_ss <= 0) max_ss = n; end;
maxVariables = min([n  p max_ss]);
if (max_steps <= 0) max_steps = 8 * maxVariables; end;
max_steps = min(8*maxVariables, max_steps);
err = 1;
l1 = 0;
if max_l1 <= 0 max_l1 = 1000; end;
b = double(zeros(p, 1));
b_prev = b;
mu = double(zeros(n, 1));
I = 1:p;
A = [];
lassoCond = 0;
step = 0;
deltar = 0;
deltax = 0;
G = g+1;

while (length(A) < maxVariables) && (step < max_steps) && (G > g) && (err > min_l2) && (l1 < max_l1)

    if lassoCond
        I = [I A(dropIdx)];
        A(dropIdx) = [];
    end

    r = y - mu;
    c = X' * r;
    [cmax cidxI] = max(abs(c(I)));
    % disp(sprintf('c(A) max, min; c(I) %g %g %g %g\n', max(c(A)), min(c(A)), max(c(I)), min(c(I))));
    cidx = I(cidxI);
    if ~lassoCond 
        A = [A cidx];
        I(cidxI) = [];
    else
        lassoCond = 0;
    end

    Gram = X(:, A)' * X(:, A);
    % For a variant of algorithm that works on well conditioned matrices
    if (cond(Gram) > 1e2)
      disp(sprintf('Index (%d %d) makes the matrix ill conditioned. Dropping', cidx, find(A == cidx)));
      A(find(A == cidx)) = [];
      continue
    end

    b_OLS = Gram \ (X(: , A)' * y);
    b_2 = inv(Gram) * (X(: , A)' * y);
    sb = sign(b(A));
    sb(length(b(A))) = sign(c(cidx));
    b_d = inv(Gram) * sb;
    % display(b_d ./ (b(A) - b_OLS));
    % disp(sprintf('Step %d, condition %g err in b_OLS %g\n', step, cond(Gram), norm(b_OLS - b_2)/norm(b_OLS)));
    
    d = X(: , A) * b_OLS - mu;

    gamma_tilde = b(A(1: end)) ./ (b(A(1: end)) - b_OLS(1: end, 1));
    gamma_tilde(gamma_tilde <= 0) = inf;
    [gamma_tilde dropIdx] = min(gamma_tilde);

    if isempty(I)
        gamma = 1;
    else
        cd = X'*d;
        temp = [(c(I) - cmax) ./ (cd(I) - cmax); (c(I) + cmax) ./ (cd(I) + cmax)];
        inds1 = find(temp>0);
        [temp2 inds2] = sort(temp(temp > 0));
        % display('temp2(1:10)'); display(temp(A)); display(temp2(1:10)); display(inds2(1:10));
        if isempty(temp2)
            printf('Error no +ve direction\n');
            return;
        end
        gamma = temp2(1);
        gamma2 = gamma;
    end

    if gamma_tilde < gamma
        lassoCond = 1;
        gamma = gamma_tilde;
    end
    b_prev = b;
    sb_prev = sign(b_prev);
    b(A) = b(A) + gamma*(b_OLS - b(A));
    sb = sign(b);
    % delta = find(sb - sb_prev);
    % if (length(delta) > 1) 
    %    disp(sprintf('Warning: step %d sign of b changed at more than one places\n', step));
    %    display(delta);
    %    display([b_prev(A) b(A)]);  
    % end

    mu = mu + gamma * d;
    step = step + 1;
    l1 = norm(b, 1);
    err = norm((y - mu), 2);

    if lassoCond
        sb(A(dropIdx)) = sb_prev(A(dropIdx));
    end
    G_array = X(:,A)' * (y - mu) ./ sb(A) / err;
    if ((max(G_array) - min(G_array)) / max(G_array) > 1e-5)
        disp(sprintf('Warning: Active sets do not seem to have same derivatives (possible numerical unstability)\nMax_g min_g %g %g %g\n', max(G_array), min(G_array), (max(G_array)-min(G_array))/max(G_array)));
    end
    G = min(G_array);
    % display(A);
    disp(sprintf('Step %d, size(A) %d, err %g, L1 %g, G %g\n', step, length(A), err, l1, G));
    % display(b(A)); 
    % display(A);
    check = (X(:, A)' * X(:, A) * b(A) - X(:, A)' * y) ./ (err * sign(b(A))); % = -g
    mx = max(check);
    mn = min(check);
    ch = (mx - mn) * 200 / abs(mx + mn);
    % disp(sprintf('%d, %d, %f, %f, %f %g %g %g\n', step, size(A, 2), mx, mn, ch, gamma, gamma_tilde, gamma2));
    % if step > 1
            % G = -((upper1 - a2) / (normb - a1));
    %        if G < g
    %            break;
    %        end
    % end
    
    %    upper1 = a2;
    %    normb = a1;
end

if (err < min_l2) 
    disp(sprintf('Applying error bound correction: quadratic equation here. Err %g min_l2 %g\n', err, min_l2));
    % Apply L2 error correction here. Update G and L1 norms.
end
if (l1 > max_l1)
    disp(sprintf('Applying L1 bound correction: L1 %g max_l1 %g\n', l1, max_l1));
    % Apply L1 error correction here. Update G.
end
if (G < g)
    disp(sprintf('Applying G  correction:G %g g %g\n', G, g));
    % Apply G correction here (already done below. Need to just move.
end

g_out = G;

tmp = b;
display('beta before');
% display(b(A));
c = X' * (y - X * b);
disp(sprintf('c(A) max, min; c(I) %g %g %g %g\n', max(abs(c(A))), min(abs(c(A))), max(c(I)), min(c(I))));

sg = sign(b(A));
if lassoCond
    sg(dropIdx) = sb_prev(A(dropIdx));
end
XA = X(:, A);

Yh = y - XA * inv(XA' * XA) * XA' * y;
Z = g * XA * inv(XA' * XA) * sg;

p = 1 - Z' * Z;
q = Yh' * Z + Z' * Yh;                                                     % Yh' * Z + Z' * Yh = 2 * Yh' * Z
r = -Yh' * Yh;

a21 = (-q + sqrt(q * q - 4 * p * r)) / (2 * p);
a22 = (-q - sqrt(q * q - 4 * p * r)) / (2 * p);
if a21 > 0 && a22 > 0
display('controversy');
elseif a21 > 0
a2 = a21;
elseif a22 > 0
a2 = a22;
end

display(a2);

b(A) = inv(XA' * XA) * (XA' * y - g * a2 * sg);

display('beta after');
% display(b(A));
c = X' * (y - X * b);
disp(sprintf('c(A) max, min; c(I) %g %g %g %g\n', max(abs(c(A))), min(abs(c(A))), max(c(I)), min(c(I))));
% verify = sum(sign(tmp) == sign(b(A)));
% verify = verify == size(A, 2);
% display('signs match?');
% display(verify);

verify = sum(sg == sign(b(A)));
verify = verify == size(A, 2);
display('signs match?');
display(verify);
% 306 319
end
