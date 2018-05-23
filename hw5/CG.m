a = [4,3,0;3,4,-1;0,-1,4]
b = [3;5;-5]
x0 = [0;0;0]
e = 0.00000000000001

x = cg(a, b, x0, e)


function x = cg(A, b, x0, e)
n = size(A, 1);
x = x0;
r = b - A * x
d = r
for k = 0: (n - 1)
    alpha = (r'*r) / (d' * A * d)
    x = x + alpha * d
    r2 = b - A * x
    if ((norm(r2) <= e) | (k == n - 1))
        x;
        break;
    end
    beta = norm(r2) ^ 2 / norm(r) ^ 2
    d = r2 + beta * d
    r = r2;
end
end
