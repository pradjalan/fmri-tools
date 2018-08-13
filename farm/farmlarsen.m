function [b,steps,G,a2,error,drop] = farmlarsen(X, y, delta, stop, Gram, storepath, verbose, g)
%LARSEN The LARS-EN algorithm for estimating Elastic Net solutions.
%
%   BETA = LARSEN(X, Y, DELTA, STOP, GRAM, STOREPATH, VERBOSE) evaluates
%   the LARS-EN algmorithm [1] using the variables in X to approximate the
%   response y given regularization parameters DELTA and STOP (LAMBDA).
%   See the function ELASTICNET for a description of parameters X, Y,
%   DELTA, STOP, STOREPATH and VERBOSE. GRAM represents an optional
%   precomputed Gram matrix of size p by p (X is n by p). The number of
%   iterations performed is returned as a second output.
%   
%   Note: The main purpose of this function is to act as an inner function
%   for the more user-friendly functions ELASTICNET and LASSO. Direct use
%   of this function requires good understanding of the algorithm and its
%   implementation.
%
%   The algorithm is a variant of the LARS-EN algorithm [1].
%
%   References
%   -------
%   [1] H. Zou and T. Hastie. Regularization and variable selection via the
%   elastic net. J. Royal Stat. Soc. B. 67(2):301-320, 2005. 
%


fprintf('farmlarsen called with delta=%g, stop=%d, storepath=%d, verbose=%d, g=%g\n', delta, stop, storepath, verbose, g);
G=0;
%% algorithm setup
[n p] = size(X);
% verbose =1;
% Determine maximum number of active variables
if delta < eps
  maxVariables = min(n,p); %LASSO
else
  maxVariables = p; % Elastic net
end

maxSteps = 8*maxVariables; % Maximum number of algorithm steps

% set up the LASSO coefficient vector
if storepath
  b = zeros(p, 2*p);
else
  b = zeros(p, 1);
  b_prev = b;
end

% current "position" as LARS travels towards lsq solution
mu = zeros(n, 1);

% Is a precomputed Gram matrix supplied?
useGram = ~isempty(Gram);

I = 1:p; % inactive set
A = []; % active set
if ~useGram
  R = []; % Cholesky factorization R'R = X'X where R is upper triangular
end

% correction of stopping criterion to fit na�ve Elastic Net
if delta > 0 && stop > 0,
  stop = stop/(1 + delta);
end

lassoCond = 0; % LASSO condition boolean
stopCond = 0; % Early stopping condition boolean
step = 1; % step count
deltar = 0; % change in residual
deltax = 0; % change in norm x
G = 0;
if verbose
  fprintf('Step\tAdded\tDropped\t\tActive set size\n');
end

%% LARS main loop
% while not at OLS solution, early stopping criterion is met, or too many
% steps have passed 
while length(A) < maxVariables && ~stopCond && step < maxSteps
  r = y - mu;

  % find max correlation
  c = X'*r;
  [cmax cidxI] = max(abs(c(I)));
  cidx = I(cidxI); % index of next active variable
  
  if ~lassoCond 
    % add variable
    if ~useGram
      [R] = cholinsert1(R, X(:,cidx), X(:,A), delta);
    end
    if verbose
      fprintf('%d\t\t%d\t\t\t\t\t%d\n', step, cidx, length(A) + 1);
    end
    A = [A cidx]; % add to active set
    I(cidxI) = []; % ...and drop from inactive set
  else
    % if a variable has been dropped, do one step with this
    % configuration (don't add new one right away) 
    lassoCond = 0;
  end

  % partial OLS solution and direction from current position to the OLS
  % solution of X_A
  if useGram
    b_OLS = Gram(A,A)\(X(:,A)'*y); % same as X(:,A)\y, but faster
  else
    b_OLS = R\(R'\(X(:,A)'*y)); % same as X(:,A)\y, but faster
  end
  d = X(:,A)*b_OLS - mu;
  
  % compute length of walk along equiangular direction
  if storepath
    gamma_tilde = b(A(1:end-1),step)./(b(A(1:end-1),step) - b_OLS(1:end-1,1));
  else
    gamma_tilde = b(A(1:end-1))./(b(A(1:end-1)) - b_OLS(1:end-1,1));
  end
  gamma_tilde(gamma_tilde <= 0) = inf;
  [gamma_tilde dropIdx] = min(gamma_tilde);

  if isempty(I)
    % if all variables active, go all the way to the OLS solution
    gamma = 1;
  else
    cd = X'*d;
    temp = [ (c(I) - cmax)./(cd(I) - cmax); (c(I) + cmax)./(cd(I) + cmax) ];
    temp = sort(temp(temp > 0)); % faster than min(temp(temp > 0)) (!)
    if isempty(temp)
      %error('SpaSM:larsen', 'Could not find a positive direction towards the next event.');
      continue
    end
    gamma = temp(1);
  end
  
  % check if variable should be dropped
  if gamma_tilde < gamma,
    lassoCond = 1;
    gamma = gamma_tilde;
  end
    
  % update beta
  if storepath
    % check if beta must grow
    if size(b,2) < step + 1
      b = [b zeros(p, size(b,2))];
    end
    b(A,step + 1) = b(A,step) + gamma*(b_OLS - b(A,step)); % update beta
  else
    b_prev = b;
    b(A) = b(A) + gamma*(b_OLS - b(A)); % update beta
  end

  % update position
  mu = mu + gamma*d;
  
  % increment step counter
  step = step + 1;

  % Early stopping at specified bound on L1 norm of beta
  if stop > 0
    if storepath
      t2 = sum(abs(b(:,step)));
      if t2 >= stop
        t1 = sum(abs(b(:,step - 1)));
        s = (stop - t1)/(t2 - t1); % interpolation factor 0 < s < 1
        b(:,step) = b(:,step - 1) + s*(b(:,step) - b(:,step - 1));
        stopCond = 1;
      end
    else
      t2 = sum(abs(b));
      if t2 >= stop
        t1 = sum(abs(b_prev));
        s = (stop - t1)/(t2 - t1); % interpolation factor 0 < s < 1
        b = b_prev + s*(b - b_prev);
        stopCond = 1;
      end
    end
  end
  
  % If LASSO condition satisfied, drop variable from active set
  if lassoCond
    if verbose
      %fprintf('%d\t\t\t\t%d\t\t\t%d\n', step, A(dropIdx), length(A)-1);
    end
    if ~useGram
      [R ] = choldelete(R, dropIdx);
    end
    I = [I A(dropIdx)]; % add dropped variable to inactive set
    A(dropIdx) = []; % ...and remove from active set
  end
  
  % Early stopping at specified number of variables
  if stop < 0
    stopCond = length(A) >= -stop;
  end
  Y = y;
  Y_new = X*b;
  a1 = norm(b,1);
  a2 = norm((Y-Y_new),2);
  if step > 2
      Gprev = G;
      G = -((upper1-a2)/(normb- a1));
      error = a2/norm(Y_new,2);
      if G < g
	  fprintf('Breaking because G=%g, Gprev=%g, active set size %d\n', G, Gprev, length(A));
          break
      end
  end
    
    upper1 = a2;

    normb = a1;

  
end

drop = length(A);
% trim beta
if storepath && size(b,2) > step
  b(:,step + 1:end) = [];
end

% return number of iterations
steps = step - 1;

% issue warning if algorithm did not converge
if step == maxSteps
  warning('SpaSM:larsen', 'Forced exit. Maximum number of steps reached.');
end

end
		  
%% Fast Cholesky insert and remove functions
% Updates R in a Cholesky factorization R'R = X'X of a data matrix X. R is
% the current R matrix to be updated. x is a column vector representing the
% variable to be added and X is the data matrix containing the currently
% active variables (not including x).
function R = cholinsert1(R, x, X, lambda)
diag_k = (x'*x + lambda)/(1 + lambda); % diagonal element k in X'X matrix
if isempty(R)
  R = sqrt(diag_k);
else
  col_k = x'*X/(1 + lambda); % elements of column k in X'X matrix
  R_k = R'\col_k'; % R'R_k = (X'X)_k, solve for R_k
  R_kk = sqrt(diag_k - R_k'*R_k); % norm(x'x) = norm(R'*R), find last element by exclusion
  R = [R R_k; [zeros(1,size(R,2)) R_kk]]; % update R
end
end
% Deletes a variable from the X'X matrix in a Cholesky factorisation R'R =
% X'X. Returns the downdated R. This function is just a stripped version of
% Matlab's qrdelete.
function R = choldelete(R,j)
R(:,j) = []; % remove column j
n = size(R,2);
for k = j:n
  p = k:k+1;
  [G,R(p,k)] = planerot(R(p,k)); % remove extra element in column
  if k < n
    R(p,k+1:n) = G*R(p,k+1:n); % adjust rest of row
  end
end
R(end,:) = []; % remove zero'ed out row
    
end