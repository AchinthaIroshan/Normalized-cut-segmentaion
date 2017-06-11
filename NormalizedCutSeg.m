clear all;
close all;
Im = imread('eye.jpg');
r = 1.5; % Spatial threshold (less than r pixels apart)

ImS   = 5;                  % Color similarity
SpS   = 6;                  % spatical similarity  

segNcut = 1;
segArea = 100;
imshow(Im);

[numRow, numCol,c] = size(Im);   
N = numRow * numCol;
V = reshape(Im, N, c); % Vertices of the graph 

%% Computing W matrix (affine matrix) 

W = sparse(N,N);        % Computing W 
WI = reshape(Im, N, 1, c); %Intnsity value vector
Xs = cat(3, repmat((1:numRow)', 1, numCol), repmat((1:numCol), numRow, 1));
Xs = reshape(Xs, N, 1, 2);  %column vector for spatial locaiton

for icol=1:numCol  
    for irow=1:numRow 
        %here a block of image is extracted where distance between 2 vertice <= r 
        %block (8 neighbour vertices when r=1)
        jcol = (icol - floor(r)) : (icol + floor(r));
        jrow = ((irow - floor(r)) :(irow + floor(r)))';
        jcol = jcol(jcol >= 1 & jcol <= numCol);
        jrow = jrow(jrow >= 1 & jrow <= numRow);
        jN = length(jcol) * length(jrow);
        %Index at vertex(V(a))
        a = irow + (icol - 1) * numRow;
        b = repmat(jrow, 1, length(jcol)) + repmat((jcol -1) * numRow, length(jrow), 1);
        %Index if neighbor vertices within the block 
        b = reshape(b, length(jcol) * length(jrow), 1);
        %Difference of distance between spatial laction of a n b
        %b is a set of points
        XB = Xs(b, 1, :); 
        XA = repmat(Xs(a, 1, :), length(b), 1);
        DXab = XA - XB;
        DXab = sum(DXab .* DXab, 3);
        constraint = find(sqrt(DXab) <= r);
        b = b(constraint);
        DXab = DXab(constraint);
        %Intensity differences (Intensity dissimilarity)
        IB = WI(b, 1, :);
        IA = repmat(WI(a, 1, :), length(b), 1);
        DifIab = IA - IB;
        DifIab = sum(DifIab .* DifIab, 3);
        W(a, b) = exp(-DifIab / (ImS*ImS)) .* exp(-DXab / (SpS*SpS));
        
    end    
end

%% ncutPartition

Seg = (1:N)';
id = 'ROOT';

% Computing D
N = length(W);
d=sum(W,2);
D = spdiags(d, 0, N, N);




[U,S] = eigs(D-W, D, 2, 'sm'); 

%returns diagonal matrix S of generalized eigenvalues and 
%full matrix V whose columns are the corresponding right eigenvectors, 
%so that (D-W)*U = D*U*S.

U2 = U(:,2); % Second smallest eigen vector

%Bipartiioning the graph where the normalized cut is mininized 
t=mean(U2);
t=fminsearch('NormcutValue',t,[],U2,W,D);
segA = find(U2 > t);
segB = find(U2 <= t);

x=(U2 > t);
x=(2*x)-1;
d=diag(D);
k=sum(d(x>0))/sum(d);

x = (U2 > t);
x = (2 * x) - 1;
d = diag(D);
k = sum(d(x > 0)) / sum(d);
b = k / (1 - k);
y = (1 + x) - b * (1 - x);
ncut = (y' * (D - W) * y) / ( y' * D * y );

[SegA  NcutA] = NormCutPartition(Seg(segA), W(segA, segA), segNcut, segArea );
[SegB  NcutB] = NormCutPartition(Seg(segB), W(segB, segB), segNcut, segArea );
Seg  = [SegA SegB];

Ncut = [NcutA NcutB];
NcutImage  = zeros(size(Im),'uint8');
for k=1:length(Seg)
 [r, c] = ind2sub(size(Im),Seg{k});
 for a=1:length(r)
 NcutImage(r(a),c(a),1:3) = uint8(round(mean(V(Seg{k}, :))));
 end
end

figure;
imshow(NcutImage)





