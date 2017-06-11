function normcut = NormcutValue(t, U2, W, D)

x = (U2 > t);
x = (2 * x) - 1;
d = diag(D);
k = sum(d(x > 0)) / sum(d);
b = k / (1 - k);
y = (1 + x) - b * (1 - x);
normcut = (y' * (D - W) * y) / ( y' * D * y );

end