function [  U,C ] = iterationNew( A, S, U, C,belta, alpha )
% iteratively updates U with C fixed and then C with U fixed
    lossPre = realmax;
    lossCur = lossPre - (lossPre * 1e-6);
    SExt = [S;zeros(1,size(S,2))];
    i = 1;
    while (abs(lossPre - lossCur) / lossPre) > 1e-7 || i < 500
        % update C with U fixed
        UExt = [U;sqrt(alpha) * ones(1,size(U,2))];
        C = C .* ((UExt' * SExt) ./ max(realmin,(UExt' * UExt * C)));
        
        % update U with C fixed
         U = U .* power(  max(( (S * C'- U * (C * C')) + 2 * belta * A * U ),realmin) ./ max(realmin,(2 * belta * (U * U') * U)), (1/4));
        
        %calculate loss function
        lossPre = lossCur;
%       lossCur = lost( A,S,U,C,alpha,belta);
        L1 = trace((S - U*C) * (S - U*C)');
        L2 = belta * trace((A - U*U') * (A - U*U')');
        L3 = alpha * sum(sum(C).^2);
        lossCur = L1 + L2 + L3;
        %disp([num2str(i),' th iteration with loss ',num2str(lossCur)]);
        i = i + 1;
        if( i > 500)
            break;
        end
    end  

end

