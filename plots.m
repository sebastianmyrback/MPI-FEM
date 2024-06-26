%base_path = "/home/s/m/smyrback/Documents/courses/SF2568-Parallel/build/";
base_path = "/Users/sebastianmyrback/Documents/doktorandprojekt/kurser/SF2568_Parallel_Computations_for_Large-Scale_Problems/project-git/build/";

% A = load(base_path + "matrix.dat");
% B = spconvert(A);
% 
% b = load(base_path + "rhs.dat");
% x = load(base_path + "x.dat");
% 
% ucg = load(base_path + "uh.dat");
% 
% 
% uh = B \ b;
% 
% C = full(B)
% 
% b


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

ucgs = cell(5,1);
ucg = [];
xs = cell(5,1);
x = [];

for i=1:5
    ucgi = load(base_path + "solution_" + num2str(i-1) + ".dat")';
    ucgs{i} = [ucgs{i}, ucgi];
    ucg = [ucg, ucgi];
    %xs{i} = [xs{i}, linspace(0,1,length(ucgi))];
    %x = [x, xs{i}];
end
x = linspace(0,1,length(ucg));

% xs{1} = x(1:21);
% xs{2} = x(21:21+19);
% xs{3} = x(21+)

lenold = 1;
for i=1:5
    lennew = length(ucgs{i});
    xs{i} = x(lenold:lenold + lennew - 1);
    lenold = lenold + lennew;
end

%u = @(x) 2*sin(2*pi*x);
u = @(x) 2*sin(2*pi*x).*sin(pi*x/10) + 10;

%x = linspace(0,1,length(ucg));

u_vals = u(x);

%plot(x, uh, 'o');
%hold on;
plot(x, u_vals, 'Color', 'black');
colormap(parula(5));
markers = ['^', 'o', 's', '*', 'p'];
legend('Exact solution', 'FontSize', 15)
for i=1:5
hold on;

plot(xs{i}, ucgs{i}, markers(i), 'MarkerSize', 12, 'DisplayName',strcat('Rank: ',num2str(i-1)))
%legend('Numerical solution, rank' + num2str(i), 'FontSize', 15)
end

xlabel('x', 'FontSize', 20)
ylabel('Solution', 'Fontsize', 20)
%norm(u_vals' - ucg)


syms x;
assume(x,'real')

u = 2*sin(2*pi*x)*sin(pi*x/10);
%u = 2*sin(2*pi*x)

f = simplify(-diff(u, x, 2))