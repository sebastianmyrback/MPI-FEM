base_path = "/home/s/m/smyrback/Documents/courses/SF2568-Parallel/build/";

A = load(base_path + "matrix.dat");
B = spconvert(A);

b = load(base_path + "rhs.dat");
x = load(base_path + "x.dat");

ucg = load(base_path + "uh.dat");


uh = B \ b;

u = @(x) 2*sin(2*pi*x);

u_vals = u(x);

plot(x, uh, 'o');
hold on;
plot(x, u_vals);
hold on;
plot(x, ucg, '^')

norm(u_vals - ucg)


syms x;
assume(x,'real')

u = 2*sin(2*pi*x)*sin(pi*x/10);
%u = 2*sin(2*pi*x)

f = simplify(-diff(u, x, 2))