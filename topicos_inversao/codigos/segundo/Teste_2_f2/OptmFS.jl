module OptmFS

using LinearAlgebra
using ForwardDiff
using Plots

##*********************************************
mutable struct Arr
    arr::Array{Any,1}
end
init_arr!(x::AbstractArray) = Arr(x)
struct VecStrc
    vec::Array{Arr,1}
end
VecStrc!(n) = VecStrc(Array{Arr,1}(undef,n))

function select(v::VecStrc,k::Integer;fr1::Int64=1,fr2::Int64=2)
x=zeros(k)
y=zeros(k)

for ic=1:k
    x[ic]=v.vec[ic].arr[fr1]
    y[ic]=v.vec[ic].arr[fr2]
end

return x,y
end
###************************

function bracket_minimum(f::Function,x::Real=0.0;
    s::Real=1e-2,k::Integer=2)

    a, ya = x, f(x)
    b, yb = a + s, f(a + s)
    if yb > ya
        a, b = b, a
        ya, yb = yb, ya
        s = -s
    end

    while true
        c, yc = b + s, f(b + s)
        if yc > yb
            return a < c ? (a, c) : (c, a)
        end
        a, ya, b, yb = b, yb, c, yc
        s *= k
    end
end

function brent_dekker(f::Function, 
    x0::Real,x1::Real;ϵ::Real=1e-7,kmax::Integer=50)

    xtol=ϵ 
    ytol=2ϵ

    y0 = f(x0)
    y1 = f(x1)
    if abs(y0) < abs(y1)
        # swap lower and upper bounds #
        x0, x1 = x1, x0
        y0, y1 = y1, y0        
    end
    x2 = x0
    y2 = y0
    x3 = x2
    mflag = true
    k=0
    while k < kmax
        k+=1
        # x-tolerance #
        if abs(x1-x0) < xtol
            return x1            
        end

        # inverse quadratic interpolation #        
        if abs(y0-y2) > ytol && abs(y1-y2) > ytol

            aux1 = x0*((y1*y2)/((y0-y1)*(y0-y2)))
            aux2 = x1*((y0*y2)/((y1-y0)*(y1-y2)))
            aux3 = x2*((y0*y1)/((y2-y0)*(y2-y1)))
            x = aux1 + aux2 + aux3               
        else
        # secant method #
            x = x1 - y1 * (x1-x0)/(y1-y0)            
        end

        # bisection method #
        δ = abs(2ϵ*abs(x1))
        aux1 = abs(x-x1)
        aux2 = abs(x1-x2)
        aux3 = abs(x2-x3)
        if (x < (3x0+x1)/4 && x > x1) ||
            (mflag && aux1 >= aux2/2) ||
            (!mflag && aux1 >= aux3/2) ||
            (mflag && aux2 < δ) ||
            (!mflag && aux3 < δ)

            x = (x0+x1)/2            
            mflag = true
        else
            mflag = false
        end
        y = f(x)
        # y-tolerance #
        if abs(y) < ytol
            return x            
        end
        x3 = x2
        x2 = x1
        if sign(y0) != sign(y)
            x1 = x
            y1 = y            
        else
            x0 = x
            y0 = y        
        end
        if abs(y0) < abs(y1)
            # swap lower and upper bounds #
            x0, x1 = x1, x0
            y0, y1 = y1, y0
        end        
    end    
end

function line_search(f::Function,
    x::AbstractArray,d::AbstractArray;ϵ::Real=1e-7)
    
    ϕ = α -> f(x + α*d)
    a, b = bracket_minimum(ϕ)
    α = brent_dekker(ϕ, a, b)
    return α
end

function backtracting_line_search(f::Function,
    x::AbstractArray,d::AbstractArray,α::Real; 
    p::Real=0.5, β::Real=1e-4)

    y, g = f(x), gradient(f,x)
    while f(x + α*d) > y + β*α*dot(g,d)
        α *= p
    end
    return α
end

function strong_backtracking_to_test(f::Function,
    x::AbstractArray, d::AbstractArray; 
    α::Real=1.0, β::Real=1.0e-4, σ::Real=0.1, αmax::Real=30.0)

    y0, g0, y_prev, α_prev = f(x), dot(gradient(f,x) ,d), NaN, 0.0
    αlo,  αhi = NaN, NaN
    ϵ = 1.0e-8

    # bracket phase
    while true
        y = f(x + α*d)
        
        if y > y0 + β*α*g0 || (!isnan(y_prev) && y >= y_prev)
            αlo, αhi = α_prev, α
            break
        end
        g = dot(gradient(f,x + α*d),d)
        if  abs(g) <= -σ*g0 #g >= σ*g0
            return α
        elseif g >= 0
            αlo, αhi = α, α_prev
            break
        end

        y_prev, α_prev, α = y, α, 2.0α
    end

    # zoom phase
    ylo = f(x + αlo*d)
    while true
        α = (αlo + αhi)/2.0
        y = f(x + α*d)
        
        println("α = ",α ,"  y = ", y)
        if y > y0 + β*α*g0 || y >= ylo
            αhi = α
        else
            g = dot(gradient(f,x + α*d),d)
            if  abs(g) <= -σ*g0 #g >= σ*g0
                println("success = ",α)
                return α
            elseif g*(αhi - αlo) >= 0
               αhi = αlo
            end
            αlo = α
        end        
        
        if y - f(x + αlo*d) < abs(y)*ϵ
            println("failed = ", α)
            return α                        
        end
    end
end


function zoom(f::Function,x::AbstractArray,d::AbstractArray,
    β::Real,σ::Real;αlo::Real,αhi::Real,kmax::Integer=10)

    α, y0, g0 = αlo, f(x), dot(gradient(f,x),d)
    ylo = f(x + αlo*d)

    k=0
    while k < kmax
        k+=1
        
        α = (αlo + αhi)/2
        y = f(x + α*d)
        
        if y > y0 + β*α*g0 || y >= ylo
            αhi =  α
        else
            g = dot(gradient(f,x + α*d),d)
            if abs(g) <= -σ*g0              
                return α
            end
            if g*(αhi - αlo) >= 0.0
               αhi = αlo      
            end
            αlo = α
        end
    end    
    return α
end

function strong_backtracking(f::Function,x::AbstractArray,d::AbstractArray;
    α::Real=1.0, αmax::Real=10.0, β::Real=1.0e-4, σ::Real=0.1, kmax::Integer=50) 

    if β < 0.0 || σ < β || σ > 1.0
        exit()
    end

    α0 = 0.0
    if α < α0 || α > αmax
        α = 1.0
    end
    y0, g0 = f(x), dot(gradient(f,x),d)
    y_prev = y0

    δ = 2.0
    k=0
    while k < kmax
        k+=1
        
        y = f(x + α*d)
        if ( y > y0 + β*α*g0 ) || 
            ( y >= y_prev && k > 1 )
           return zoom(f,x,d,β,σ,αlo=α0,αhi=α)            
        end

        g = dot(gradient(f,x + α*d),d)
        if abs(g) <= -σ*g0            
            return α
        end
        
        if g >= 0.0            
            return zoom(f,x,d,β,σ,αlo=α,αhi=α0)
        end

        α0, y_prev = α, y

        if α0 <= α && α <= αmax
            α = α0 * δ        
        else
            α = 1.0            
        end
    end
end


abstract type DescentMethod end
mutable struct GradientDescent <: DescentMethod
	α::Real
    g::AbstractArray
end
function init_sd!(α::Real,x::AbstractArray)
   g::Array{Real,1} = zeros(Real,length(x)) 
   return GradientDescent(α,g)
end
function step!(M::GradientDescent,f::Function,
    x::AbstractArray)
	α = M.α

    M.g = gradient(f,x)
    
    # descent direction #
    d = -M.g

    # descent direction normalized #
    d = d/norm(d,Inf)

    #α = backtracting_line_search(f, x, d, α)
    α = strong_backtracking_to_test(f, x, d, α=α)

    M.α = α
    return x + α*d
end

function steepest_descent!(M::GradientDescent, 
    f::Function, x::AbstractArray, kmax::Integer)  

   xall=VecStrc!(kmax)   
   ϵ = 1e-5

   k=0
   while k < kmax
     k += 1

     x0 = x
     xall.vec[k] = init_arr!(x0) 

     y = f(x0)
     x = step!(M,f,x0)
     
     if norm(M.g,2) < ϵ && (abs(y - f(x)) <= abs(y)*ϵ)
        println("Optimal")
        break
     end
          
     #println("k = ",k," and f(x) = ", y)     
   end

   return x, xall, k
end


mutable struct ConjugateGradientDescent <: DescentMethod
    α::Real
    d::AbstractArray
    g::AbstractArray
end
function init_cg!(f::Function,x::AbstractArray, α::Real)
    g::Array{Real,1} = gradient(f,x)  
    
    # descent direction #
    d::Array{Real,1} = -g    

    # descent direction normalized #
    d = d/norm(d,Inf) 
    return ConjugateGradientDescent(α, d, g)
end
function step!(M::ConjugateGradientDescent, 
    f::Function, x::AbstractArray)
    α, d, g = M.α, M.d, M.g

    #α = backtracting_line_search(f, x, d, α)
    α = strong_backtracking_to_test(f, x, d, α=α)              

    xc = x + α*d
    gc = gradient(f,xc)

    β = dot(gc, gc-g)/dot(g,g)
    β = max(0, β)
    dc = -gc + β*d

    M.α, M.g, M.d = α, gc, dc
    return xc
end

function conjugate_gradient_descent!(M::ConjugateGradientDescent, 
    f::Function, x::AbstractArray, kmax::Integer)  

    xall=VecStrc!(kmax)  
    ϵ = 1e-3

    k=0
    while k < kmax
      k += 1   
      
      x0 = x
      xall.vec[k] = init_arr!(x0) 

      y = f(x0)
      x= step!(M,f,x0)

      if norm(M.g,2) < ϵ && (abs(y - f(x)) <= abs(y)*ϵ)
         println("Optimal")
         break
      end
      
      #println("k = ",k," and f(x) = ", y)
    end

    return x, xall, k
end


mutable struct Momentum <: DescentMethod
    α::Real # learning rate
    β::Real # momentum decay
    v::AbstractArray # momentum
    g::AbstractArray # gradient
end
function init_momentum!(x::AbstractArray,α::Real,β::Real)
    n=length(x)
    v::Array{Real,1} = zeros(Real,n)
    g::Array{Real,1} = zeros(Real,n)
    return Momentum(α,β,v,g)
end
function step!(M::Momentum, f::Function, x::AbstractArray)
    α, β, v = M.α, M.β, M.v

    M.g = gradient(f,x)

    v = β*v - α*M.g
    M.v = v
    return x + v
end

function momentum!(M::Momentum, 
    f::Function, x::AbstractArray, kmax::Integer)  
      
    xall=VecStrc!(kmax)    
    ϵ = 1e-4

    k=0
    while k < kmax
      k += 1   

      x0 = x
      xall.vec[k] = init_arr!(x0) 
 
      y= f(x0)
      x = step!(M,f,x0)
      
      if norm(M.g,2) < ϵ && (abs(y - f(x)) <= abs(y)*ϵ)
         println("Optimal")
         break
      end

      #println("k = ",k," and f(x) = ", y)     
    end

    return x, xall, k
end


mutable struct NesterovMomentum <: DescentMethod
    α::Real # learning rate
    β::Real # momentum decay
    v::AbstractArray # momentum
    g::AbstractArray # gradient
end
function init_nesterov!(x::AbstractArray,α::Real,β::Real)
    n=length(x)
    v::Array{Real,1} = zeros(Real,n)
    g::Array{Real,1} = zeros(Real,n)
    return NesterovMomentum(α,β,v,g)
end
function step!(M::NesterovMomentum, f::Function, x::AbstractArray)
    α, β, v = M.α, M.β, M.v
    
    M.g = gradient(f,x + β*v)

    v = β*v - α*M.g
    M.v = v
    return x + v
end

function nesterov_momentum!(M::NesterovMomentum, 
    f::Function, x::AbstractArray, kmax::Integer)  
      
    xall=VecStrc!(kmax)
    ϵ = 1e-4

    k=0
    while k < kmax
      k += 1   
      
      x0 = x
      xall.vec[k] = init_arr!(x0) 

      y = f(x0)
      x = step!(M,f,x0)
      
      if norm(M.g,2) < ϵ && (abs(y - f(x)) <= abs(y)*ϵ)
         println("Optimal")
         break
      end
      
      #println("k = ",k," and f(x) = ", y)     
    end

    return x, xall, k
end


mutable struct RMSProp <: DescentMethod
    α # learning rate
    γ # decay
    ε # small value
    s # sum of squared gradient
    g # gradient
end
function init_rmsprop!(x::AbstractArray;
    α::Real=0.1, γ::Real=0.9)

    ε::Real = 1.0e-8
    n=length(x)
    s::Array{Real,1} = zeros(Real,n)
    g::Array{Real,1} = zeros(Real,n)
    return RMSProp(α,γ,ε,s,g)
end
function step!(M::RMSProp, f::Function, x::AbstractArray)
    α, γ, ε, s = M.α, M.γ, M.ε, M.s
    
    M.g = gradient(f,x)

    s = γ*s + (1.0-γ)*(M.g.*M.g)

    M.s = s
    return x - α*M.g ./ (sqrt.(s) .+ ε)
end

function rmsprop!(M::RMSProp, 
    f::Function, x::AbstractArray, kmax::Integer)  
      
    xall=VecStrc!(kmax)
    ϵ = 1e-5

    k=0
    while k < kmax
      k += 1   

      x0 = x
      xall.vec[k] = init_arr!(x0) 
      
      y = f(x0)
      x = step!(M,f,x0)
      
      if norm(M.g,2) < ϵ && (abs(y - f(x)) <= abs(y)*ϵ)
         println("Optimal")
         break
      end
      
      #println("k = ",k," and f(x) = ", y)     
    end

    return x, xall, k
end


mutable struct Adam <: DescentMethod
    α::Real  # learning rate
    γv::Real # decay
    γs::Real # decay
    ε::Real  # small value
    g::AbstractArray # gradient
    v::AbstractArray # 1st moment estimate
    s::AbstractArray # 2nd moment estimate    
end
function init_adam!(x::AbstractArray;α::Real=0.001,
    γv::Real=0.9,γs::Real=0.999)
    ε::Real=1e-8
    
    n=length(x)
    g::Array{Real,1} = zeros(Real,n)
    v::Array{Real,1} = zeros(Real,n)
    s::Array{Real,1} = zeros(Real,n)
    return Adam(α,γv,γs,ε,g,v,s)
end
function step!(M::Adam, f::Function, x::AbstractArray, k::Integer)
    α, γv, γs, ε = M.α, M.γv, M.γs, M.ε
    s, v, M.g = M.s, M.v, gradient(f,x)

    v = γv*v + (1.0-γv)*M.g
    s = γs*s + (1.0-γs)*M.g.*M.g

    v_hat = v ./ (1.0 - γv^k)
    s_hat = s ./ (1.0 - γs^k)

    M.v, M.s = v, s
    return x - α*v_hat ./ (sqrt.(s_hat) .+ ε)
end

function adam!(M::Adam, 
    f::Function, x::AbstractArray, kmax::Integer)  
      
    xall=VecStrc!(kmax)
    ϵ = 1e-3

    k=0
    while k < kmax
      k += 1   

      x0 = x
      xall.vec[k] = init_arr!(x0) 
      
      y = f(x0)
      x = step!(M,f,x0,k)
      
      if norm(M.g,2) < ϵ && (abs(y - f(x)) <= abs(y)*ϵ)
         println("Optimal")
         break
      end
      
      #println("k = ",k," and f(x) = ", y)     
    end

    return x, xall, k
end


mutable struct Newton <: DescentMethod
    α::Real
    g::AbstractArray
    d::AbstractArray
end
function init_newton!(x::AbstractArray,α::Real=1.0)
    n=length(x)
    g = Vector{Real}(undef,n)
    d = Vector{Real}(undef,n)
    return Newton(α,g,d)
end
function step!(M::Newton,f::Function,x::AbstractArray)
    α0 = M.α
    ϵ = 1e-7

    H = hessian(f,x) 
    M.g = gradient(f,x)
    
    # descent direction #
    d = -H\M.g

    # descent direction normalized #
    #d = M.d/norm(M.d,Inf)

    #α = backtracting_line_search(f,x,d,α0)
    α = strong_backtracking_to_test(f, x, d, α=α0) 
    if (α0 - α) < ϵ
        α = 1.0        
    end

    M.α = α
    return x + α*d
end

function newtons_method!(M::Newton,f::Function,x::AbstractArray,kmax::Integer)
    
    xall=VecStrc!(kmax)
    ϵ = 1e-7

    k = 0
    while k < kmax
        k=k+1        
        
        x0 = x
        xall.vec[k] = init_arr!(x0) 

        y = f(x0)        
        x = step!(M,f,x0)        

        if norm(M.g,2) < ϵ && (abs(y - f(x)) <= abs(y)*ϵ)
            println("Optimal")
            break
        end

        #println("k = ",k," and f(x) = ", y)
    end
    return x, xall, k
end


mutable struct BFGS <: DescentMethod
    α::Real
    g::AbstractArray    
    Q::AbstractArray
end
function init_bfgs!(x::AbstractArray,α::Real=1.0)
    m = length(x)
    
    g = Vector{Real}(undef,m)
    Q = Matrix{Real}(1.0I, m, m)
    # scale #
    #δ = 0.1
    #g = gradient(f,x)
    #Q = (δ/norm(g))*Q

    return BFGS(α,g,Q)
end
function step!(M::BFGS, f::Function, x_prev::AbstractArray,k::Int64)
    Q, g_prev = M.Q, gradient(f,x_prev)   

    α0 = M.α
    # descent direction #
    d = -Q*g_prev
    
    # descent direction normalized #
    #d = d/norm(d,Inf)

    #α = backtracting_line_search(f, x_prev, d, α0)
    α = strong_backtracking(f, x_prev, d, α=α0)
    M.α = α

    x = x_prev + α*d

    g = gradient(f,x)
    δ = x - x_prev
    γ = g - g_prev
    
    # scale #
    if k == 1         
        Q = ((γ'*δ)/(γ'*γ))*Q        
    end
    Q = Q - (δ*γ'*Q + Q*γ*δ')/(δ'*γ) +
            (1.0 + (γ'*Q*γ)/(δ'*γ)) * (δ*δ')/(δ'*γ)         

    M.g, M.Q = g, Q            
    return x
end


function bfgs!(M::BFGS,f::Function,x::AbstractArray,kmax::Integer)
    
    xall=VecStrc!(kmax)
    ϵ = 1e-5

    k = 0
    while k < kmax
        k=k+1
        
        x0 = x
        xall.vec[k] = init_arr!(x0) 

        y = f(x0)        
        x = step!(M,f,x0,k)
        
        if norm(M.g,2) < ϵ && (abs(y - f(x)) <= abs(y)*ϵ)
            println("Optimal")
            break
        end

        #println("k = ",k," and f(x) = ", y)
    end
    return x, xall, k
end


mutable struct LimitedMemoryBFGS <: DescentMethod
    m::Integer
    α::Real
    g::AbstractArray
    δs::AbstractArray
    γs::AbstractArray
    qs::AbstractArray
end
function init_lbfgs!(x::AbstractArray,α::Real=1.0,m::Integer=3)
    n=length(x)
    g = Vector{Real}(undef,n)
    δs = []
    γs = []
    qs = []    
    return LimitedMemoryBFGS(m,α,g,δs,γs,qs)
end
function step!(M::LimitedMemoryBFGS, f::Function , x_prev::AbstractArray)
    δs, γs, qs, g_prev = M.δs, M.γs, M.qs, gradient(f,x_prev)
    m = length(δs)
    
    α0 = M.α

    if m > 0
        q = g_prev
        for i in m : -1 : 1
            qs[i] = copy(q)            
            q -= (δs[i]⋅q)/(γs[i]⋅δs[i])*γs[i]
        end

        # scale #
        z = (γs[m] .* δs[m] .* q) / (γs[m]⋅γs[m])
        for i in 1 : m
            z += δs[i]*(δs[i]⋅qs[i] - γs[i]⋅z)/(γs[i]⋅δs[i])
        end
        # descent direction #
        d = -z

        # descent direction normalized # 
        #d = d/norm(d,Inf)
        
        #α = backtracting_line_search(f, x_prev, d, α0)
        α = strong_backtracking(f, x_prev, d, α=α0)
        M.α = α

        x = x_prev + α*d
    else
        # descent direction #
        d = -g_prev

        # descent direction normalized #
        d = d/norm(d,Inf)

        #α = backtracting_line_search(f, x_prev, d, α0)
        α = strong_backtracking(f, x_prev, d, α=α0)
        M.α = α

        x = x_prev + α*d
    end

    g = gradient(f,x)
    M.g = g

    push!(δs, x - x_prev) 
    push!(γs, g - g_prev)
    push!(qs, zeros(length(x_prev)))
    while length(δs) > M.m
        popfirst!(δs) 
        popfirst!(γs)
        popfirst!(qs)
    end
    
    return x
end

function lbfgs!(M::LimitedMemoryBFGS,f::Function,x::AbstractArray,kmax::Integer)
    
    xall=VecStrc!(kmax)
    ϵ = 1e-6

    k = 0
    while k < kmax
        k=k+1
        
        x0 = x
        xall.vec[k] = init_arr!(x0) 

        y = f(x0)        
        x = step!(M,f,x0)
        
        if norm(M.g,2) < ϵ && (abs(y - f(x)) <= abs(y)*ϵ)
            println("Optimal")
            break
        end

        #println("k = ",k," and f(x) = ", y)
    end
    return x, xall, k
end

function graph(x::Array{Float64,1},y::Array{Float64,1},
    z::Array{Float64,2};  
    p1::Array{Float64,1}=[ ],p2::Array{Float64,1}=[ ], 
    tp::Bool=true, lx::String=" ", ly::String=" ",
    tn::String="",sav::Bool=false)
    
    xmin = minimum(x)
    xmax = maximum(x)
    ymin = minimum(y)
    ymax = maximum(y)

    q =contour(x,y,z,    
    title = tn, titlefontsize = 14,xlabel = lx, ylabel = ly, labelfont = font(10, "Courier"),
    xlims = (xmin,xmax), ylims = (ymin,ymax), 
    xticks = xmin:2:xmax,yticks = ymin:2:ymax,
    xtickfont = font(10, "Courier"),
    ytickfont = font(10, "Courier"),      
    levels=[1,2,3,5,10,20,50,100],
    nx=length(x), ny=length(y),
    colorbar = false,
    linecolor = [:yellow, :green, :orange],
    legend = :none, 
    linestyle = :solid, 
    linewidth = 2,    
    framestyle = :box,
    grid = false,
    reuse = true,
    size = (500,500)
    )

    if tp == true
    plot!(p1,p2,
    color = :black, 
    linestyle = :solid,
    linewidth = 2,
    legend =:none,
    xlims = (xmin,xmax),
    ylims = (ymin,ymax),
    xticks = xmin:2:xmax,
    yticks = floor(ymin):2:floor(ymax),    
    framestyle = :box,
    grid = false
    )
          
    else
    scatter!([p1[1]], [p1[2]],
    markershape = :circle,
    markersize = 5,
    markercolor = :red,
    markerstrokecolor = :black,
    legend = :none,
    label = "Optimal"
    )

    scatter!([p2[1]], [p2[2]],
    markershape = :circle,
    markersize = 5,
    markercolor = :blue,
    markerstrokecolor = :black,
    legend = :none,
    label = "Initial"
    )
    end

    if sav == true
        savefig(q,tn*".pdf")
    end
    return q
end
#*****************************************************************

rosenbrock2d(x; a=3.0, b=24.0) = a*x[1]^2 + b*(x[2])^2
gener_rosenbrock(x; a=1.0, b=5.0) = b*(x[2] - x[1]^2)^2 + (a-x[1])^2 + b*(x[3] - x[2]^2)^2 + (a-x[2])^2 + b*(x[4] - x[3]^2)^2 + (a-x[3])^2 + b*(x[5] - x[4]^2)^2 + (a-x[4])^2 + b*(x[6] - x[5]^2)^2 + (a-x[5])^2

gradient(f,x0) = ForwardDiff.gradient(f,x0)
hessian(f,x0) = ForwardDiff.hessian(f,x0)

#**************************************************************#

function Main()

# two parameters #
#
x1::Array{Float64,1} = range(-2,stop=2,length=200)
x2::Array{Float64,1} = range(-2,stop=2,length=200)
nd1=length(x1)
nd2=length(x2)
frame1::Int64=1
frame2::Int64=2

rb::Array{Float64,2} = zeros(Float64,nd1,nd2) 

rb = [rosenbrock2d([x,y]) for y in x2, x in x1]
xo = [0.0, 0.0] # minimum point #
x0 = [-5.75, 1.75] # start point #

q = graph(x1,x2,rb,p1=xo,p2=x0,tp=false,tn="Optimal and initial design points",sav=true)
display(q)
#
#=
# six parameters #
#
x1::Array{Float64,1} = range(-2,stop=2,length=10)
x2::Array{Float64,1} = range(-2,stop=2,length=10)
x3::Array{Float64,1} = x1
x4::Array{Float64,1} = x1
x5::Array{Float64,1} = x1
x6::Array{Float64,1} = x1
nd1=length(x1)
nd2=length(x2)
frame1::Int64=5
frame2::Int64=6

xo = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0] # minimum point #
x0 = [-1.7, -1.7, -1.7, -1.7, -1.7, -1.7] # start point #

grb::Array{Float64,6} = zeros(Float64,nd1,nd1,nd1,nd1,nd1,nd1)
grb = [gener_rosenbrock([a,b,c,d,e,f]) for f in x6, e in x5, d in x4, c in x3, b in x2, a in x1]

rb = grb[7,7,7,7,:,:]
q = graph(x1,x2,rb,p1=[xo[frame1],xo[frame2]],p2=[x0[frame1],x0[frame2]],tp=false,tn="Optimal and initial design points",sav=true)
display(q)
#
=#
#*************#
#  Inversion  #
#*************#
# 
fobj = rosenbrock2d
#fobj = gener_rosenbrock

########Steepest Descent###########
println("********************")
println("* Steepest Descent *")
println("********************")

kmax = 5000
α0 = 6.0
x0 = [-5.75, 1.75]
#x0 = [-1.7, -1.7, -1.7, -1.7, -1.7, -1.7]

xvec=VecStrc!(kmax)
sd=init_sd!(α0,x0)
x,xvec,k = steepest_descent!(sd, fobj, x0, kmax)

println("iter = ",k, "  f(x) = ",fobj(x))
println("x st = ",x0)
println("x op = ",xo)
println("x pr = ",x)

cx,cy=select(xvec,k,fr1=frame1,fr2=frame2)
q1 = graph(x1,x2,rb,p1=cx,p2=cy,tn="Steepest descent",sav=true)
#
#
############Conjugate Gradient Descent#################
println("**********************")
println("* Conjugate Gradient *")
println("**********************")

kmax = 5000
α0 = 6.0
x0 = [-5.75, 1.75]
#x0 = [-1.7, -1.7, -1.7, -1.7, -1.7, -1.7]

xvec=VecStrc!(kmax)
cg=init_cg!(fobj,x0,α0)
x,xvec,k = conjugate_gradient_descent!(cg, fobj, x0, kmax)

println("iter = ",k, "  f(x) = ",fobj(x))
println("x st = ",x0)
println("x op = ",xo)
println("x pr = ",x)

cx,cy=select(xvec,k,fr1=frame1,fr2=frame2)
q2 = graph(x1,x2,rb,p1=cx,p2=cy,tn="Conjugate gradient",sav=true)
#
#=
############Momentum#################
println("************")
println("* Momentum *")
println("************")

kmax = 5000
α = 0.01 
β = 0.85 # β-hyperparameter is defined in the range [0.0, 1.0]
#x0 = [-1.7, -1.7]
x0 = [-1.7, -1.7, -1.7, -1.7, -1.7, -1.7]

xvec=VecStrc!(kmax)
m=init_momentum!(x0,α,β)
x,xvec,k = momentum!(m, fobj, x0, kmax)

println("iter = ",k, "  f(x) = ",fobj(x))
println("x st = ",x0)
println("x op = ",xo)
println("x pr = ",x)

cx,cy=select(xvec,k,fr1=frame1,fr2=frame2)
q3 = graph(x1,x2,rb,p1=cx,p2=cy,tn="Momentum",sav=true)
=#
#
#=
println("*******************")
println("* RMS Propagation *")
println("*******************")

kmax = 5000
# hyperparameters #
α = 0.03
γ = 0.999
#x0 = [-1.7, -1.7]
x0 = [-1.7, -1.7, -1.7, -1.7, -1.7, -1.7]

xvec=VecStrc!(kmax)
rmsp=init_rmsprop!(x0,α=α,γ=γ)
x,xvec,k = rmsprop!(rmsp, fobj, x0, kmax)

println("iter = ",k, "  f(x) = ",fobj(x))
println("x st = ",x0)
println("x op = ",xo)
println("x pr = ",x)

cx,cy=select(xvec,k,fr1=frame1,fr2=frame2)
q3 = graph(x1,x2,rb,p1=cx,p2=cy,tn="RMSProp",sav=true)
#
println("*********************")
println("* Nesterov Momentum *")
println("*********************")

kmax = 10000
α = 0.001
β = 0.9 # β-hyperparameter is defined in the range [0.0, 1.0] #
#x0 = [-1.7, -1.7]
x0 = [-1.7, -1.7, -1.7, -1.7, -1.7, -1.7]

xvec=VecStrc!(kmax)
nm=init_nesterov!(x0,α,β)
x,xvec,k = nesterov_momentum!(nm, fobj, x0, kmax)

println("iter = ",k, "  f(x) = ",fobj(x))
println("x st = ",x0)
println("x op = ",xo)
println("x pr = ",x)

cx,cy=select(xvec,k,fr1=frame1,fr2=frame2)
q4 = graph(x1,x2,rb,p1=cx,p2=cy,tn="Nesterov momentum",sav=true)
#
#
println("********************************")
println("* Adaptive Momentum Estimation *")
println("********************************")

kmax = 5000
# hyperparameters #
α = 0.1
γv = 0.9
γs = 0.99
#x0 = [-1.7, -1.7]
x0 = [-1.7, -1.7, -1.7, -1.7, -1.7, -1.7]

xvec=VecStrc!(kmax)
ad=init_adam!(x0,α=α,γv=γv,γs=γs)
x,xvec,k = adam!(ad, fobj, x0, kmax)

println("iter = ",k, "  f(x) = ",fobj(x))
println("x st = ",x0)
println("x op = ",xo)
println("x pr = ",x)

cx,cy=select(xvec,k,fr1=frame1,fr2=frame2)
q5 = graph(x1,x2,rb,p1=cx,p2=cy,tn="Adam",sav=true)
=#
#
println("*******************")
println("* Newton's Method *")
println("*******************")

kmax = 1000
#α0 = 6.0 # Use this value to backtracting_line_search #
α0 = 1.0
x0 = [-5.75, 1.75]
#x0 = [-1.7, -1.7, -1.7, -1.7, -1.7, -1.7]

xvec=VecStrc!(kmax)
nm=init_newton!(x0,α0)
x,xvec,k = newtons_method!(nm, fobj, x0, kmax)

println("iter = ",k, "  f(x) = ",fobj(x))
println("x st = ",x0)
println("x op = ",xo)
println("x pr = ",x)

cx,cy=select(xvec,k,fr1=frame1,fr2=frame2)
q6 = graph(x1,x2,rb,p1=cx,p2=cy,tn="Newton",sav=true)
#
#=
println("***************")
println("* Method BFGS *")
println("***************")

kmax = 1000
α0 = 1.0 
#x0 = [-1.7, -1.7]
x0 = [-1.7, -1.7, -1.7, -1.7, -1.7, -1.7]

xvec=VecStrc!(kmax)
bf=init_bfgs!(x0,α0)
x,xvec,k = bfgs!(bf, fobj, x0, kmax)

println("iter = ",k, "  f(x) = ",fobj(x))
println("x st = ",x0)
println("x op = ",xo)
println("x pr = ",x)

cx,cy=select(xvec,k,fr1=frame1,fr2=frame2)
q7 = graph(x1,x2,rb,p1=cx,p2=cy,tn="BFGS",sav=true)
=#
#=
println("*****************")
println("* Method L-BFGS *")
println("*****************")

kmax = 1000
α0 = 1.0 
m = 5
#x0 = [-1.7, -1.7]
x0 = [-1.7, -1.7, -1.7, -1.7, -1.7, -1.7]

xvec=VecStrc!(kmax)
lbf=init_lbfgs!(x0,α0,m)
x,xvec,k = lbfgs!(lbf, fobj, x0, kmax)

println("iter = ",k, "  f(x) = ",fobj(x))
println("x st = ",x0)
println("x op = ",xo)
println("x pr = ",x)

cx,cy=select(xvec,k,fr1=frame1,fr2=frame2)
q8 = graph(x1,x2,rb,p1=cx,p2=cy,tn="L-BFGS",sav=true)
#
display(plot(q1, q2, q3, q4, q5, q6, q7, q8, layout = (2, 4), legend = false, size = (1130,590), reuse = true) )
=#

end#main

end#module
