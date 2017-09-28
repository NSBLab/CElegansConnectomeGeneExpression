function p = probMMatches(N,n1,n2,m,doExact)
%-------------------------------------------------------------------------------
% Computes the probability of having m matches from two binary vectors of
% length N given n1 1s in vector 1, and n2 1s in vector 2
%-------------------------------------------------------------------------------

if nargin < 5
    doExact = true;
end

% Not possible to get any other matches than what was observed:
% (simple case):
if (N-n2) < (n1-m)
    p = 0;
    return
end

%===============================================================================
if m > n1 || m > n2
    error('More matches than ones');
    p = NaN;
end
if doExact
    % Evaluate large factorials using symbolic math toolbox (in subfunctions):
    syms N_s n1_s n2_s m_s

    if (m==0) || (n1-m==0)
        p = subs(factorial(n1_s)*factorial(n2_s)*factorial(N_s-n1_s)*factorial(N_s-n2_s)/...
                (factorial(N_s)*factorial(n2_s-m_s)*factorial(n1_s-m_s)*factorial(m_s)*...
                    factorial(N_s-n1_s-n2_s+m_s)),[N_s,n1_s,n2_s,m_s],[N,n1,n2,m]);
    else
        p = subs(nchoosek(n2_s,m_s)*nchoosek(N_s-n2_s,n1_s-m_s)/nchoosek(N_s,n1_s),...
                        [N_s,n1_s,n2_s,m_s],[N,n1,n2,m]);
    end
else
    % Smaller numbers can be done numerically:
    warning('off')
    if (m==0) || (n1-m==0)
        p = factorial(n1)*factorial(n2)*factorial(N-n1)*factorial(N-n2)/...
                (factorial(N)*factorial(n2-m)*factorial(n1-m)*factorial(m)*...
                    factorial(N-n1-n2+m));
    else
        p = nchoosek(n2,m)*nchoosek(N-n2,n1-m)/nchoosek(N,n1);
    end
    warning('on')
end

end
