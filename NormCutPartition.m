function [Seg  Ncut] = NormCutPartition(I, W, segNcut, segArea)
N = length(W);
d = sum(W, 2);
D = spdiags(d, 0, N, N); 
[U,S] = eigs(D-W, D, 2, 'sm');
U2 = U(:, 2);
t = mean(U2);
t = fminsearch('NcutValue', t, [], U2, W, D);
A = find(U2 > t);
B = find(U2 <= t);

ncut = NcutValue(t, U2, W, D);
if (length(A) < segArea || length(B) < segArea) || ncut > segNcut
    Seg{1}   = I;
%    Id{1}   = id;          % for debugging
    Ncut{1} = ncut;        % for duebuggin
    return;
end

[SegA  NcutA] = NormCutPartition(I(A), W(A, A), segNcut, segArea);
[SegB  NcutB] = NormCutPartition(I(B), W(B, B), segNcut, segArea);

Seg   = [SegA SegB];

Ncut = [NcutA NcutB];
end
