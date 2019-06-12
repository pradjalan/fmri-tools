% Change Log
% L1 norm(b) changed everywhere to L1 norm(b(A))
% Implemented correction equations
% Bug in Line 56 Fixed. When the matrix became ill-conditioned the index
% was dropped and next iteration of loop was executed.
% Now after dropping the index, the Gram matrix is recomputed and the same 
% loop iteration is executed with the changed inactive set and updating 
% b(A) for same active set as in the previous iteration.

% Added functionality for cross-validation. k1 and k2 are start and end
% indices of train set respectively, k3 and k4 are start and end indices 
% for cross-validation set respectively.
function [b, l1_out, err_out, g_out, step,test_err,train_err] = larsen(flat_mri, model_index, max_l1, min_l2, g, max_ss, max_steps,k1,k2,k3,k4)

if nargin <= 7 %if set indices not provides,run without cross-validation
    k1=0;
    k2=0;
    k3=0;
    k4=0;
end
if k1>0 || k2>0 || k3>0 || k4>0  
    storepath=true;
else 
    storepath=false;
    test_err=inf;
    train_err=inf;
end
X = double(flat_mri);
X(:, model_index) = [];
% Set default values 
if storepath
    if k1 <= 0 k1 = 1; end;
    if k2 <= 0 k2 = size(X,1); end;
    if k3 <= 0 k3 = 1; end;
    if k4 <= 0 k4 = size(X,1); end;

    X = X(k1: k2 - 1, :);
    y = flat_mri(k1+1: k2, model_index);
    X = zscore(X) / sqrt(size(X, 1) - 1);
    y = zscore(y) / sqrt(size(y, 1) - 1);
    Xtot = double(flat_mri);
    Xtot(:, model_index) = [];
    Xtot = Xtot(1: end - 1, :);
    ytot = flat_mri(2: end, model_index);
    Xtot = zscore(Xtot) / sqrt(size(Xtot, 1) - 1);
    ytot = zscore(ytot) / sqrt(size(ytot, 1) - 1);
    Xtest = double(flat_mri);
    Xtest(:, model_index) = [];
    Xtest = Xtest(k3: k4 - 1, :);
    ytest = flat_mri(k3+1: k4, model_index);
    Xtest = zscore(Xtest) / sqrt(size(Xtest, 1) - 1);
    ytest = zscore(ytest) / sqrt(size(ytest, 1) - 1);
else
    X = X(1: end - 1, :);
    y = flat_mri(2: end, model_index);
    X = zscore(X) / sqrt(size(X, 1) - 1);
    y = zscore(y) / sqrt(size(y, 1) - 1);
end
[n p] = size(X);

if max_l1 <= 0 max_l1 = 1000; end;
if (g <= 0) g = 0.001; end;
if (max_ss <= 0) max_ss = n; end;
maxVariables = min([n  p max_ss]);
if (max_steps <= 0) max_steps = 8 * maxVariables; end;
max_steps = min(8 * maxVariables, max_steps);

err = 1;
l1 = 0;

if storepath
    active_set_size=min(8 * maxVariables, max_steps);
    active_mat=zeros(max_steps,active_set_size);  % to store active set index for path
    b = zeros(p,100);
    mu = double(zeros(k2-k1, 1));
else
    b = double(zeros(p, 1));
    b_prev=b;
    mu = double(zeros(n, 1));
end
I = 1:p;
A = [];
lassoCond = 0;
step = 1;
G = g + 1;

while (length(A) < maxVariables) && (step < max_steps) && (G > g) && (err > min_l2) && (l1 < max_l1)

	if lassoCond
		I = [I A(dropIdx)];
		A(dropIdx) = [];
	end

	r = y - mu;
	c = X' * r;
	[cmax cidxI] = max(abs(c(I)));
	cidx = I(cidxI);
	if ~lassoCond 
		A = [A cidx];
		I(cidxI) = [];
	else
		lassoCond = 0;
	end

	Gram = X(:, A)' * X(:, A);
	if (cond(Gram) > 1e2)
		disp(sprintf('Index (%d %d) makes the matrix ill conditioned. Dropping', cidx, find(A == cidx)));
		A(find(A == cidx)) = [];
        Gram = X(:, A)' * X(:, A);
        %continue
    end
    
	b_OLS = Gram \ (X(: , A)' * y);
	d = X(: , A) * b_OLS - mu;
    if storepath
        active_mat(step+1, : )=[A zeros( 1,( size(active_mat,2) - size(A,2) )  ) ]; %storing active set for current beta's
        gamma_tilde = b(A(1:end-1),step)./(b(A(1:end-1),step) - b_OLS(1:end-1,1));
    else
        gamma_tilde = b(A) ./ (b(A) - b_OLS);
    end
	gamma_tilde(gamma_tilde <= 0) = inf;
	[gamma_tilde dropIdx] = min(gamma_tilde);

	if isempty(I)
		gamma = 1;
	else
		cd = X'*d;
		temp = [(c(I) - cmax) ./ (cd(I) - cmax); (c(I) + cmax) ./ (cd(I) + cmax)];
		[temp2 inds2] = sort(temp(temp > 0));
		if isempty(temp2)
			printf('Error no +ve direction\n');
			return;
		end
		gamma = temp2(1);
	end

	if gamma_tilde < gamma
		lassoCond = 1;
		gamma = gamma_tilde;
    end
    
    if storepath
    % check if beta must grow
    if size(b,2) < step + 1
      b = [b zeros(p, 100)];
    end
    b(A,step + 1) = b(A,step) + gamma*(b_OLS - b(A,step)); % update beta
    else
        b_prev = b;
        b(A) = b(A) + gamma*(b_OLS - b(A));
    end
    
	mu = mu + gamma * d;
	step = step + 1;
	l1 = norm(b(A), 1);
	err = norm(y - mu, 2);

	G_array = abs(X(:,A)' * (y - mu) / err);
	if ((max(G_array) - min(G_array)) / max(G_array) > 1e-5)
		disp(sprintf('Warning: Active sets do not seem to have same derivatives (possible numerical unstability)\nMax_g min_g %g %g %g\n', max(G_array), min(G_array), (max(G_array)-min(G_array))/max(G_array)));
	end
	G = min(G_array);
    
end

if storepath && size(b,2) > step
  b(:,step + 1:end) = [];
end
if storepath && size(active_mat,1) > step
  active_mat(step + 1:end,:) = [];
end


if storepath
    sb = sign(b(A,step)); 
    if lassoCond
        sb(dropIdx) = sign(b(A(dropIdx),step-1));
    end
else
    sb = sign(b(A));
    if lassoCond
        sb(dropIdx) = sign(b_prev(A(dropIdx)));
    end
end

if storepath
    muu=Xtest*b; % selecting beta's with min_err with cross validation set
    min_err=inf;
    AIdx=0;
    for i=1:size(b,2)
        if norm(ytest-muu(:,i),2) < min_err
            min_err=norm(ytest-muu(:,i),2);
            AIdx=i;
        end
    end
    b=b(:,AIdx);
    A=active_mat(AIdx,:); % active set for corresponding beta's
    A=A(A~=0);
    muu=X(:,A)*b(A);
    train_err=norm(y-X(:,A)*b(A),2); % recomputation of error and G
    G_array = abs(X(:,A)' * (y - muu) / train_err);
	G = min(G_array);
    test_err=norm(ytest-Xtest(:,A)*b(A),2);
    err=norm(ytot-Xtot(:,A)*b(A),2);
    l1 = norm(b(A), 1);
end

if ~storepath
    
        XA = X(:, A);
        sb = sign(b(A));
        if (l1 > max_l1)
            disp(sprintf('Applying L1 bound correction: L1 %g max_l1 %g\n', l1, max_l1));
            l1 = norm(b_prev(A), 1);
            delta = sum(sign(b_prev(A)) .* (b_OLS - b_prev(A)));
            gamma = (max_l1 - l1) / delta;
            b(A) = b_prev(A) + gamma * (b_OLS - b_prev(A));
            l1 = max_l1;
            mu = X(:, A) * b(A);
            err = norm(y - mu, 2);
            G_array = abs(X(:,A)' * (y - mu) / err);
            G = min(G_array);
        end
        if (err < min_l2)
            disp(sprintf('Applying error bound correction: quadratic equation here. Err %g min_l2 %g\n', err, min_l2));
            yhyh = y' * (y - XA * inv(XA' * XA) * XA' * y);
            zz = sb' * inv(XA' * XA) * sb;
            G = sqrt((min_l2 * min_l2 - yhyh) / (min_l2 * min_l2 * zz));
            err = min_l2;
            b(A) = inv(XA' * XA) * (XA' * y - G * err * sb);
            l1 = norm(b(A), 1);
        end
        if (G < g)
            disp(sprintf('Applying G  correction:G %g g %g\n', G, g));
            yhyh = y' * (y - XA * inv(XA' * XA) * XA' * y);
            zz = sb' * inv(XA' * XA) * sb;
            err = sqrt(yhyh / (1 - g * g * zz));
            G = g;
            b(A) = inv(XA' * XA) * (XA' * y - G * err * sb);
            l1 = norm(b(A), 1);
        end
        verify = sum(sign(b(A)) == sb);
        if verify ~= size(A, 2)
            %display([sb sign(b(A))]);
            disp(sprintf('Signs dont match!\n'));
        end
end

step=step-1;
l1_out = l1;
err_out = err;
g_out = G;
end
