% INPUT:
% x: [Real(lambdas) ; Imag(lambdas)]
% h: output measurements along trajectories vectorized trajectory per trajectory

% OUTPUT:
% Returns the function value and gradient of the function p = -h'*L*(L'*L)^{-1}L'*h
function [f,grad,L,A] = eigOptim_grad(x,h,Ntraj,TrajLen)
N = numel(x) / 2;
lams_R = x(1:N);
lams_I = x(N+1:end);
lams = lams_R + 1i*lams_I;

Lam = complex(zeros(TrajLen,N));
for i = 1:N
    Lam(:,i) = lams(i).^([0:TrajLen-1]');
end

% Build derivatives of L wrt the real part of lambda. The derivative wrt
% the imaginary part is 1i * the derivative wrt the real part.

Ldiff = cell(1,N);
vals_L = []; indr_L = []; indc_L = [];
for i = 1:N
    indr = [1:TrajLen*Ntraj]';
    indc = kron( Ntraj*(i-1) + (1:Ntraj)',ones(TrajLen,1));
    vals = kron(ones(Ntraj,1), [0;Lam(1:end-1,i)].*([0:TrajLen-1]') )  ;
    
    vals_L = [vals_L ; kron(ones(Ntraj,1), [Lam(:,i)])];
    indr_L = [indr_L ; indr];
    indc_L = [indc_L ; indc];
    
    Ldiff{i} = sparse(indr,indc,vals, Ntraj*TrajLen , Ntraj*N);
end
L = sparse(indr_L,indc_L,vals_L);



% Compute gradients
grad_R = complex(zeros(N,1));
grad_I = complex(zeros(N,1));


%A = inv(L'*L)*L'*h;
A = (L'*L) \ (L'*h);


% Function value
f = norm(h)^2 - real(h'*L*A);

if(nargout > 1)    
    for i = 1:N
        
        grad_R(i) = -h'*Ldiff{i}*A;
        grad_R(i) = grad_R(i) + grad_R(i)';
        tmp = Ldiff{i}'*L;
        grad_R(i) = grad_R(i) + A'*(tmp+tmp')*A;
        
        
        grad_I(i) = -h'*(1i*Ldiff{i})*A;
        grad_I(i) = grad_I(i) + grad_I(i)';
        tmp = (1i*Ldiff{i})'*L;
        grad_I(i) = grad_I(i) + A'*(tmp+tmp')*A;
    end   
    grad = real([grad_R ; grad_I]);
end




end

