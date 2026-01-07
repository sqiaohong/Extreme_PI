function out = ncs_random_smooth(x, y, varargin)
p = inputParser;
p.addParameter('NRepeat',   10,    @(v)isnumeric(v)&&isscalar(v)&&v>=1);
p.addParameter('SpanYears', [15 30], @(v)isnumeric(v)&&numel(v)==2&&v(1)>0&&v(2)>=v(1));
p.addParameter('ReturnAll', false,  @(v)islogical(v)||ismember(v,[0 1]));
p.addParameter('Seed',      [],     @(v)isempty(v)||(isscalar(v)&&isnumeric(v)));
p.parse(varargin{:});
NRepeat   = p.Results.NRepeat;
SpanYears = p.Results.SpanYears(:).';
ReturnAll = p.Results.ReturnAll;
Seed      = p.Results.Seed;

if ~isempty(Seed)
    rng(Seed);
end

x = x(:); y = y(:);
ok = isfinite(x) & isfinite(y);
x = x(ok); y = y(ok);
[x, idx] = sort(x); y = y(idx);

n = numel(x);
assert(n >= 8, 'sample size is too small');

y_each = zeros(n, NRepeat);
knots_all = cell(NRepeat,1);

for k = 1:NRepeat
    
    knots = gen_random_knots(x(1), x(end), SpanYears);
    knots_all{k} = knots;

 
    Phi = natural_cubic_basis(x, knots);   
    Xls = [ones(n,1), Phi];                
    beta = Xls \ y;                       

    y_each(:,k) = Xls * beta;
end


y_smooth = mean(y_each, 2);


out = struct('x', x, 'y_smooth', y_smooth);
if ReturnAll
    out.y_each    = y_each;
    out.knots_all = knots_all;
end


function knots = gen_random_knots(xmin, xmax, span)

a = span(1); b = span(2);
L = xmax - xmin;
if L < a*2
    
    knots = (xmin + xmax)/2;
    return;
end


pos = xmin;
knots = [];
while true
    gap = a + (b-a)*rand;
    pos = pos + gap;
    if pos < xmax
        knots(end+1) = pos; 
    else
        break;
    end
end


if numel(knots) < 1

    knots = (xmin + xmax)/2;
end


knots = unique(knots);
end

end


function Phi = natural_cubic_basis(X, knots)

X = X(:);
K = numel(knots);
if K < 1

    Phi = X;
    return;
end
t = sort(knots(:));
n = numel(X);

Phi = zeros(n, K);         
Phi(:,1) = X;

tK = t(end);
dK = ppart_cube(X - tK);  


for i = 1:K-1
    ti = t(i);
    numer = ppart_cube(X - ti) - dK;
    denom = (tK - ti);
    Phi(:, i+1) = numer ./ denom;
end


if K >= 2
    dK_1 = Phi(:, K);      
    for j = 2:K-1
        Phi(:, j) = Phi(:, j) - dK_1;
    end
    Phi(:, K) = [];       
end
end


function y = ppart_cube(t)
y = max(t, 0).^3;
end
