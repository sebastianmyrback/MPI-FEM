%base_path = "/home/s/m/smyrback/Documents/courses/SF2568-Parallel/build/";
base_path = "/Users/sebastianmyrback/Documents/doktorandprojekt/kurser/SF2568_Parallel_Computations_for_Large-Scale_Problems/project-git/build/";

A = load(base_path + "matrix.dat");
B = spconvert(A);

b = load(base_path + "rhs.dat");
x = load(base_path + "x.dat");

ucg = load(base_path + "uh.dat");


uh = B \ b;

C = full(B)

b


% correct b
%  50.0000
%    99.5077
%     0.2292
%     0.9514
%     1.2169
%     0.7798
%    -0.2443
%    -1.4041
%    -2.1085
%    98.0886
%   100.0000


%u = @(x) 2*sin(2*pi*x);
u = @(x) 2*sin(2*pi*x).*sin(pi*x/10) + 10;

u_vals = u(x);

%plot(x, uh, 'o');
%hold on;
plot(x, u_vals);
hold on;
plot(x, ucg, '^')

norm(u_vals - ucg)


syms x;
assume(x,'real')

u = 2*sin(2*pi*x)*sin(pi*x/10);
%u = 2*sin(2*pi*x)

f = simplify(-diff(u, x, 2))