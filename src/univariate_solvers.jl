bexport golden_search, broyden, bisection, newton

function golden_search(f::Function, a::F, b::F, args...; tol::F = 1e-10, max_iters = 100) where F
    φ = (1+sqrt(5))/2
    for iter in 1:max_iters
        x = (b+φ*a)/(φ+1)
        y = a + (b-x)
        fx = f(x,args...)
        fy = f(y,args...)
        boolf = fx ≤ fy
        bool = x < y
        a = boolf*bool*a + boolf*(1-bool)*y + (1-boolf)*bool*x + (1-boolf)*(1-bool)*a
        b = boolf*bool*y + boolf*(1-bool)*b + (1-boolf)*bool*b + (1-boolf)*(1-bool)*x
        err = abs(x-y)
        iter > 1 && err < tol && begin
            _x = (x+y)/2
            return _x, f(_x,args...)
        end
    end
end
function broyden()
end
function bisection()
end
function newton()
end
