module STomog

using LinearAlgebra
using Random
using Statistics
using PyCall
#using Plots
#pyplot()

plt = pyimport("matplotlib.pyplot")
color = pyimport("matplotlib.colors")
plt.rcParams["text.usetex"] = "True"
 
#Random.seed!(123) # Setting the seed

const Nv=4
const ndim=2
const eps=1.e-6
const TOL=1.0e-6
const ABSTOL=1.e-5
const RELTOL=1.e-3

mutable struct Tgrd
n1::Int64
n2::Int64
o1::Float64
o2::Float64
d1::Float64
d2::Float64
end


mutable struct Topr
nr::Int64              # numero de linhas
nc::Int64              # numero de colunas
nsize::Int64           # dimensao de col[] e val[]
nzero::Int64           # numero de elementos nao nulos <= nsize
row::Array{Int64,1}    # row[1:nr+1] = 
col::Array{Int64,1}    # col[nsize]  = 
val::Array{Float64,1}  # val
norm::Float64
#
function Topr(nrows::Int64,ncols::Int64,size::Int64)
nr=nrows
nc=ncols
nsize=size
nzero=0
row=zeros(Int64,nr+1)
col=zeros(Int64,nsize)
val=zeros(Float64,nsize)
norm=0.0
new(nr,nc,nsize,nzero,row,col,val,norm)
end
end


mutable struct Tdat
nd::Int64               # number of data values
xs::Array{Float64,ndim}    # 
xr::Array{Float64,ndim}
od::Array{Float64,1}    # observed data
function Tdat(nd::Int64)
xs=zeros(Float64,ndim,nd)
xr=zeros(Float64,ndim,nd)
od=zeros(Float64,nd)
new(nd,xs,xr,od)
end
end


function saxpy!(adj::Bool,add::Bool,oper::Topr,x::Array{Float64,1},y::Array{Float64,1})

if size(x,1) != oper.nc 
   return false
end 

if size(y,1) != oper.nr 
   return false
end 

if adj == false
   if add == false
      for ir=1:oper.nr
          y[ir]= 0.0
      end
   end

# Gather:
   @inbounds for ir=1:oper.nr
       for ie=oper.row[ir]:oper.row[ir+1]-1
               y[ir] += oper.val[ie] * x[oper.col[ie]]
       end
   end

else

   if add == false
      for ic=1:oper.nc
          x[ic]= 0.0
      end
   end

# Scatter:
   @inbounds for ir=1:oper.nr
       for ie=oper.row[ir]:oper.row[ir+1]-1
               x[oper.col[ie]] += oper.val[ie] * y[ir]
       end
   end

end

return true

end


function cgls!(oper::Topr,lambda::Float64,sig2::Float64,
    x::Array{Float64,1},b::Array{Float64,1},
    x0::Array{Float64,1},saxpy::Function,maxit::Int64)

@assert( size(x,1)  == oper.nc )
@assert( size(x0,1) == oper.nc )
@assert( size(b,1) == oper.nr )

p::Array{Float64,1}=zeros(Float64,oper.nc)
q::Array{Float64,1}=zeros(Float64,oper.nr)
r::Array{Float64,1}=zeros(Float64,oper.nr)
s::Array{Float64,1}=zeros(Float64,oper.nc)

@inbounds for ic=1:oper.nc
x[ic]=x0[ic]
#x[ic]=0.0
end

saxpy(false,false,oper,x,q)
@inbounds for ir=1:oper.nr
r[ir]=b[ir]-q[ir]
end

saxpy(true,false,oper,p,r)
@inbounds for ic=1:oper.nc
s[ic]=p[ic]-lambda*x[ic]
end

# Initialize
@inbounds for ic=1:oper.nc
p[ic]=s[ic]
end

norms0::Float64=norm(s,2)
gamma::Float64=norms0*norms0
normx::Float64=norm(x,2)
xmax::Float64=normx
alpha::Float64=0.0
delta::Float64=0.0
shift::Float64=lambda #*oper.norm

k::Int64=0
flag::Int32=0

while ( k < maxit ) && ( flag == 0 )

k=k+1
saxpy(false,false,oper,p,q)
delta=norm(q,2)^2+shift*norm(p,2)^2
if delta <= 0.0
flag=3
end

if delta == 0.0
delta=TOL
end
alpha=gamma/delta

@inbounds for ic=1:oper.nc
x[ic] = x[ic] + alpha * p[ic]
end

@inbounds for ir=1:oper.nr
r[ir] = r[ir] - alpha * q[ir]
end

saxpy(true,false,oper,s,r)
@inbounds for ic=1:oper.nc
s[ic] = s[ic] - shift * (x[ic]-x0[ic])
end
norms = norm(s,2)
gamma1= gamma
gamma = norms*norms
beta  = gamma/gamma1

@inbounds for ic=1:oper.nc
p[ic] = s[ic] + beta * p[ic]
end

normx=norm(x,2)
xmax=max(xmax,normx)

if norms < TOL*norms0
flag=1
end

end #while

println("norms0 = ",norms0);
println("norms = ",norm(s,2));
println("normr = ",norm(r,2));
println("norm1(x) = ",norm(x,1))
println("normx = ",normx);

if k == maxit
flag=2
end

if abs(alpha)*norm(p,2) <= sqrt(TOL)*xmax
flag=3
end


return flag

end # cgls


function cross2d(x::Array{Float64,1},y::Array{Float64,1})
        return x[1]*y[2]-x[2]*y[1]
end


function inside_cell(x::Array{Float64,1},vrt::Array{Float64,2})
for j=1:Nv
    jp=mod(j,Nv)+1;
    if  cross2d(vrt[1:ndim,jp]-vrt[1:ndim,j],x-vrt[1:ndim,jp]) < 0.0
        return false;
    end
end
return true
end


function ray_cell!(grd::Tgrd, 
                   x0::Array{Float64,1},
                   nu::Array{Float64,1},
                   i1::Int64,i2::Int64,
                   x1::Array{Float64,1})
#=
%
% raytracing in a homogeneous cell
%
% input
%
% x0(1:2) - origin at a cell boundary
% n(1:2)  - unit vector at ray direction
%
% output
%
% x(1:2)  - ray exit point
% s       - ray length
%
%
=#
EPS::Float64=1.e-3*min(grd.d1,grd.d2)
vrt::Array{Float64,2}=zeros(ndim,Nv)
alpha::Array{Float64,1}=zeros(Float64,Nv)
vrt[1,1]=grd.o1+Float64(i1-1)*grd.d1
vrt[2,1]=grd.o2+Float64(i2-1)*grd.d2
vrt[1,2]=vrt[1,1]+grd.d1
vrt[2,2]=vrt[2,1]
vrt[1,3]=vrt[1,2]
vrt[2,3]=vrt[2,2]+grd.d2
vrt[1,4]=vrt[1,3]-grd.d1
vrt[2,4]=vrt[2,3]

if inside_cell(x0+EPS*nu,vrt) == false
   s=0.0
   x .= x0
   return s
end

nu .= nu / sqrt(dot(nu,nu));

for j=1:Nv
    jp=mod(j,Nv)+1;
    a = cross2d(nu,vrt[1:ndim,jp]-vrt[1:ndim,j]); 
    if abs(a) > 1.0e-9*EPS
            b = cross2d(vrt[1:ndim,jp]-vrt[1:ndim,j],x0[1:ndim]-vrt[1:ndim,j]);
            alpha[j] = b / a;
            #println("alpha[",j,"]= ",alpha[j])
    end
end
s=grd.d1+grd.d2
for j=1:Nv
    if alpha[j] > 0.0 
       s = min(s,alpha[j])
    end
end

for j=1:ndim
x1[j]= x0[j] + s*nu[j];
end

return s

end

function raytracing!(grd::Tgrd,dat::Tdat,oper::Topr)

EPS::Float64=1.e-6*min(grd.d1,grd.d2)
TOL::Float64=10.0*EPS

xs::Array{Float64,1}=zeros(Float64,ndim)
xr::Array{Float64,1}=zeros(Float64,ndim)
x0::Array{Float64,1}=zeros(Float64,ndim)
x1::Array{Float64,1}=zeros(Float64,ndim)
nu::Array{Float64,1}=zeros(Float64,ndim)
nzero::Int64=0

i1::Int64=0
i2::Int64=0
length::Float64=0.0
s::Float64=0.0


oper.norm=0.0
for id=1:dat.nd
        xs .= dat.xs[1:ndim,id]
        xr .= dat.xr[1:ndim,id]
        nu .= xr-xs

        length=sqrt(dot(nu,nu))
        nu .= nu /length

s=0.0
i1=0
i2=0

x1 .= xs
oper.row[id]=nzero+1
#
# Raytracing loop
while abs(s-length) > TOL 

       x0 .= x1

       i1=Int64(floor((x0[1]+EPS*nu[1]-grd.o1)/grd.d1))+1
       i2=Int64(floor((x0[2]+EPS*nu[2]-grd.o2)/grd.d2))+1
        
       ds=ray_cell!(grd,x0,nu,i1,i2,x1)

       if  ds > 0.0
           nzero += 1
           oper.col[nzero]=i1+(i2-1)*grd.n1
           oper.val[nzero]=ds
           oper.norm += ds*ds
           s += ds
       end
 
end # ray
dat.od[id]=s
end # data
oper.nzero=nzero
oper.row[oper.nr+1]=nzero+1
oper.norm=sqrt(oper.norm)
return nothing
end

function test_model!(n1::Int64,n2::Int64,mdl::Array{Float64,1})
@assert ( size(mdl,1) == n1*n2 )

@inbounds for i2=1:n2
for i1=1:n1
    mdl[i1+(i2-1)*n1]=1.0
end
end

# slowness model
i1beg=6
i1=i1beg
i2=6
mdl[i1+(i2-1)*n1]=0.5
i1=i1beg+1
for i2=5:6
mdl[i1+(i2-1)*n1]=0.5
end
i1=i1beg+2
for i2=3:8
mdl[i1+(i2-1)*n1]=0.5
end
i1=i1beg+3
for i2=4:7
mdl[i1+(i2-1)*n1]=0.5
end
i1=i1beg+4
for i2=4:6
mdl[i1+(i2-1)*n1]=0.5
end
i1=i1beg+5
for i2=3:4
mdl[i1+(i2-1)*n1]=0.5
end
i1=i1beg+6
for i2=2:3
mdl[i1+(i2-1)*n1]=0.5
end
i1=i1beg+7
i2=2
mdl[i1+(i2-1)*n1]=0.5
for i2=5:6
mdl[i1+(i2-1)*n1]=2.0
end
i1=i1beg+8
for i2=4:7
mdl[i1+(i2-1)*n1]=2.0
end
#
end

function denoiser!(grd::Tgrd,x::Array{Float64,1},xf::Array{Float64,1})
tile::Array{Float64,1}=zeros(Float64,9)
@inbounds for i2=2:grd.n2-1
    for i1=2:grd.n1-1
        # moving median
        it=0
        for l2=-1:1
            for l1=-1:1
                it = it + 1
                tile[it]=x[(i1+l1)+(i2+l2-1)*grd.n1]
            end
        end
        xf[i1+(i2-1)*grd.n1]=median(tile)
    end
end
return nothing
end


function addnoise!(dat::Tdat,perc::Float64)
sig::Float64=0.0
dmin::Float64=1e32;
dmax::Float64=0.0
@inbounds for id=1:dat.nd
    dmax=max(dmax,abs(dat.od[id]))
    dmin=min(dmin,abs(dat.od[id]))
end

#sig=sqrt(3.0/2.0)*perc*sig
sig=perc*(dmax-dmin)

@inbounds for id=1:dat.nd
        dat.od[id] += sig*(2.0*rand(Float64)-1.0)
end
return sig
end


function shrink(a::Float64,kappa::Float64)
    return max(0.0,a-kappa)-max(0.0,-a-kappa)
end


function shrink(a::Array{Float64,1},kappa::Float64)
    n=size(a,1)
    soft::Array{Float64,1}=zeros(Float64,n)
    for i=1:n
        soft[i]=max(0.0,a[i]-kappa)-max(0.0,-a[i]-kappa)
    end
return soft
end


function normest(oper::Topr,saxpy::Function;tol::Float64=1.0e-06,maxit::Int64=500)

    p::Array{Float64,1}=zeros(Float64,oper.nc)
    q::Array{Float64,1}=zeros(Float64,oper.nr)
    s::Array{Float64,1}=zeros(Float64,oper.nc)    
    x::Array{Float64,1}=zeros(Float64,oper.nc)    
    normx::Float64=0.0
    tau::Float64=0.0
    tautil::Float64=0.0
    
    @inbounds for ic=1:oper.nc
        x[ic] = rand(Float64)
    end

    k=0
    while k < maxit

        k = k + 1

        tautil = tau

        saxpy(false,false,oper,x,q)
        saxpy(true,false,oper,s,q)
        normx = norm(x,2) 
        @inbounds for ic=1:oper.nc
            p[ic] = s[ic]/normx
        end  

        tau = sqrt(norm(p,2))

        if  abs(tau - tautil) < tol*tau
            println("optimal")
            break
        end

        @inbounds for ic=1:oper.nc
            x[ic] = p[ic]
        end

        println("t = ",tau)
        println("|t - t0| = ", abs(tau - tautil))
        
    end

return tau
end


function normestfb(oper::Topr,saxpy::Function;tol::Float64=1.0e-06,maxit::Int64=500)

    q::Array{Float64,1}=zeros(Float64,oper.nr)
    s::Array{Float64,1}=zeros(Float64,oper.nc)    
    x::Array{Float64,1}=zeros(Float64,oper.nc)    
    normx::Float64=0.0
    tau::Float64=0.0
    tautil::Float64=0.0
    
    @inbounds for ic=1:oper.nc
        x[ic] = rand(Float64)
    end

    k=0
    while k < maxit

        k = k + 1

        tautil = tau

        saxpy(false,false,oper,x,q)
        saxpy(true,false,oper,s,q)
        normx = norm(x,2)         

        tau = norm(s,2)/normx

        if  abs(tau - tautil) < tol*tau
            println("optimal")
            break
        end

        @inbounds for ic=1:oper.nc
            x[ic] = s[ic]
        end

        println("t = ",tau)
        println("|t - t0| = ", abs(tau - tautil))
        
    end

return tau
end


function red_sd!(grd::Tgrd,oper::Topr,lambda::Float64,sig2::Float64,
              x::Array{Float64,1},b::Array{Float64,1},
              x0::Array{Float64,1},saxpy::Function,denoiser::Function,maxit)
# Initialize
@assert( size(x,1)  == oper.nc )
@assert( size(x0,1) == oper.nc )
@assert( size(b,1) == oper.nr )

p::Array{Float64,1}=zeros(Float64,oper.nc)
q::Array{Float64,1}=zeros(Float64,oper.nr)
r::Array{Float64,1}=zeros(Float64,oper.nr)
s::Array{Float64,1}=zeros(Float64,oper.nc)

alpha::Float64=0.0
delta::Float64=0.0
shift::Float64=lambda #*oper.norm

objls::Array{Float64,1}=zeros(Float64,maxit)
objr::Array{Float64,1}=zeros(Float64,maxit)
obj::Array{Float64,1}=zeros(Float64,maxit)
iter::Array{Float64,1}=zeros(Float64,maxit)

@inbounds for ic=1:oper.nc
    x[ic]=x0[ic]
end

saxpy(false,false,oper,x,q)
@inbounds for ir=1:oper.nr
    r[ir]=b[ir]-q[ir]
end

saxpy(true,false,oper,s,r)
@inbounds for ic=1:oper.nc
   s[ic]=s[ic]-shift*(x[ic]-x0[ic])
end


norms0::Float64=norm(s,2)
gamma::Float64=norms0*norms0
normx::Float64=norm(x,2)
xmax::Float64=normx
rnorm::Float64=0.0

alpha=2.0e-2/(((1.0/sig2)+lambda)*oper.norm)

# RED fixed point iterations with median denoiser
k=0
while k < maxit
        
k=k+1
iter[k]=k

saxpy(false,false,oper,x,q)
@inbounds for ir=1:oper.nr
    r[ir]=b[ir]-q[ir]
end
saxpy(true,false,oper,s,r)
@inbounds for ic=1:oper.nc
    s[ic] = s[ic] - shift * (x[ic]-x0[ic])
end

@inbounds for ic=1:oper.nc
    x[ic] = x[ic] + alpha * s[ic]
end

# projection
for ic=1:oper.nc
         x[ic]=min(max(x[ic],0.2),5.0)
end 

# median denoise
denoiser(grd,x,x0)

# objective function
rnorm=norm(r,2)
objls[k]=0.5*rnorm*rnorm
objr[k]=shift*dot(x,x .- x0)
 
obj[k]=0.5*rnorm*rnorm+shift*dot(x,x .- x0)

if norm(s,2) < RELTOL * norms0
   println("dual")
   break
end

println("norm(s)= ",norm(s,2))

end #while

println("ITER= ",k," norm(s)= ",norm(s,2)," norm0= ",norms0)

graph(iter,objls,nt="Least Square ",nx="k",ny="1/2|| Ax -b ||_2^2",sav=true)
graph(iter,objr,nt="Regularization",nx="k",ny="ld x^T(x - f(x))",sav=true)
graph(iter,obj,nt="Objective function",nx="k",ny="1/2|| Ax -b ||_2^2 + ld x^T(x - f(x))",sav=true)

return nothing
end # denoiser!


function red_fp!(grd::Tgrd,oper::Topr,lambda::Float64,sig2::Float64,
              x::Array{Float64,1},b::Array{Float64,1},
              x0::Array{Float64,1},saxpy::Function,denoiser::Function,
              maxit1::Int64,maxit2::Int64)

# Initialize
@inbounds for ic=1:oper.nc
    x[ic]=x0[ic]
end

# RED fixed point iterations with median denoiser
k=0
while k < maxit1
        
k=k+1

flag=cgls!(oper,lambda,sig2,
          x,b,
          x0,saxpy,maxit2)

# projection
        for ic=1:oper.nc
                x[ic]=min(max(x[ic],0.2),5.0)
        end 
        
# median denoise
denoiser(grd,x,x0)

end #while

return nothing
end # denoiser!


function red_admm!(grd::Tgrd,oper::Topr,
                   lambda::Float64,sig2::Float64,beta::Float64,
                   x::Array{Float64,1},b::Array{Float64,1},
                   x0::Array{Float64,1},saxpy::Function,denoiser::Function,
                   maxit1::Int64,maxit2::Int64,maxit3::Int64)

u::Array{Float64,1}=zeros(Float64,oper.nc)
v::Array{Float64,1}=zeros(Float64,oper.nc)
xhat::Array{Float64,1}=zeros(Float64,oper.nc)

p::Array{Float64,1}=zeros(Float64,oper.nc)
q::Array{Float64,1}=zeros(Float64,oper.nr)
r::Array{Float64,1}=zeros(Float64,oper.nc)
s::Array{Float64,1}=zeros(Float64,oper.nr)
norms::Float64=oper.norm
lshift::Float64=lambda*norms
bshift::Float64=beta*norms
rnorm::Float64=0.0
tau::Float64=2.0
murho::Float64=10.0
alpha::Float64=1.0
nroot::Float64=sqrt(Float64(oper.nc))

# Initialize
@inbounds for ic=1:oper.nc
    x[ic]=x0[ic]
end

@inbounds for ic=1:oper.nc
    v[ic]=x0[ic]
end

@inbounds for ic=1:oper.nc
    u[ic]=0.0
end

obj::Array{Float64,1}=zeros(Float64,maxit1)
rpri::Array{Float64,1}=zeros(Float64,maxit1)
sdual::Array{Float64,1}=zeros(Float64,maxit1)
epspri::Array{Float64,1}=zeros(Float64,maxit1)
epsdual::Array{Float64,1}=zeros(Float64,maxit1)
rhotil::Array{Float64,1}=zeros(Float64,maxit1)
iter::Array{Float64,1}=zeros(Float64,maxit1)

# RED fixed point iterations with median denoiser
k=0
while k < maxit1
        
k=k+1
iter[k]=k

@inbounds for ic=1:oper.nc
    p[ic]=x[ic]
end
saxpy(false,false,oper,p,q)
for ir=1:oper.nr
    s[ir]=b[ir]-q[ir]
end
saxpy(true,false,oper,r,s)

for j=1:maxit2
        for ic=1:oper.nc
            r[ic]=r[ic]+bshift*(p[ic]-v[ic]+u[ic])
        end 
        saxpy(false,false,oper,r,q)
        mu=dot(r,r)/(dot(q,q) + bshift*dot(r,r))
        for ic=1:oper.nc
            p[ic]=p[ic]+mu*r[ic]
        end 
        for ir=1:oper.nr
            s[ir]=s[ir]-mu*q[ir]
        end
        saxpy(true,false,oper,r,s)
        # projection
        for ic=1:oper.nc
                p[ic]=min(max(p[ic],0.2),5.0)
        end 
end

@inbounds for ic=1:oper.nc
        x[ic]=p[ic]
end

@inbounds for ic=1:oper.nc
        p[ic]=v[ic]
end

# over-relaxation
@inbounds for ic=1:oper.nc
    xhat[ic]=alpha*x[ic] - (1.0 - alpha)*p[ic]
end

# objective function
rnorm=norm(r,2)
denoiser(grd,p,r)
obj[k]=0.5*rnorm*rnorm + 0.5*lambda*(dot(p,p .- r))

# v-update 
for j=1:maxit3
    denoiser(grd,p,r)
    @inbounds for ic=1:oper.nc
        p[ic]=(lshift*r[ic]+bshift*(x[ic]+u[ic]))/(bshift+lshift)
    end
end
# primal and dual residue
println("norm(p-v) = ",norm(p .- v,2))
println("norm(x-v) = ",norm(x .- v,2))
rpri[k]=norm(x .- v,2)
sdual[k]=norm(beta.*(p .- v,2))

# varying penalty parameter
if mod(k,30) == 0
if rpri[k] > murho*sdual[k]
    beta = tau*beta
    u .= u ./ tau
end
    
if sdual[k] > murho*rpri[k]
    beta = beta/tau
    u .= u .* tau
end
bshift=beta*norms
end
rhotil[k]=beta

epsdual[k] = ABSTOL*nroot + RELTOL * norm(beta.*u,2)
epspri[k] = ABSTOL*nroot + RELTOL * max(norm(x,2),norm(v,2))

if (rpri[k] <= epspri[k]) && (sdual[k] <= epsdual[k])
   println("dual")
   break
end

@inbounds for ic=1:oper.nc
        v[ic]=p[ic]
end

@inbounds for ic=1:oper.nc
    u[ic] += x[ic]-v[ic]
end

end # while

println("k = ",k)
println("norm(b) = ",norm(b,2))
println("norm(s) = ",norm(s,2))
println("norm(r) = ",rnorm)
println("norm(x-v) = ",norm( x .- v,2))
println("norm(x) = ",norm(x,2))
println("norm(v) = ",norm(v,2))

graph(iter,obj,nt="Objective function",nx="k",ny="1/2|| Ax -b ||_2^2 + ld/2 < x, x-f(x)> ",sav=true)
graph(iter,rpri,nt="Primal residual",nx="k",ny="|| r ||_2",sav=true)
graph(iter,sdual,pl=false,nt="Dual residual",nx="k",ny="|| s ||_2",sav=true)
graph(iter,epspri,nt="Primal tolerance",nx="k",ny="eps^primal",sav=true)
graph(iter,epsdual,pl=false,nt="Dual tolerance",nx="k",ny="eps^dual",sav=true)
graph(iter,rhotil,nt="Penality term",nx="k",ny="rho",sav=true)

return nothing
end # red_admm


function red_admm!(grd::Tgrd,oper::Topr,
    lambda::Float64,sig2::Float64,beta::Float64,
    x::Array{Float64,1},b::Array{Float64,1},
    x0::Array{Float64,1},saxpy::Function,denoiser::Function,
    maxit1::Int64,maxit2::Int64,maxit3::Int64)

u::Array{Float64,1}=zeros(Float64,oper.nc)
v::Array{Float64,1}=zeros(Float64,oper.nc)
xhat::Array{Float64,1}=zeros(Float64,oper.nc)

p::Array{Float64,1}=zeros(Float64,oper.nc)
q::Array{Float64,1}=zeros(Float64,oper.nr)
r::Array{Float64,1}=zeros(Float64,oper.nc)
s::Array{Float64,1}=zeros(Float64,oper.nr)
normb::Float64=norm(b,2)
lshift::Float64=lambda*normb
bshift::Float64=beta*normb
rnorm::Float64=0.0
tau::Float64=2.0
murho::Float64=10.0
alpha::Float64=1.0
nroot::Float64=sqrt(Float64(oper.nc))

# Initialize
@inbounds for ic=1:oper.nc
    x[ic]=x0[ic]
end

@inbounds for ic=1:oper.nc
    v[ic]=x0[ic]
end

@inbounds for ic=1:oper.nc
    u[ic]=0.0
end

obj::Array{Float64,1}=zeros(Float64,maxit1)
rpri::Array{Float64,1}=zeros(Float64,maxit1)
sdual::Array{Float64,1}=zeros(Float64,maxit1)
epspri::Array{Float64,1}=zeros(Float64,maxit1)
epsdual::Array{Float64,1}=zeros(Float64,maxit1)
rhotil::Array{Float64,1}=zeros(Float64,maxit1)
iter::Array{Float64,1}=zeros(Float64,maxit1)

# RED fixed point iterations with median denoiser
k=0
while k < maxit1

k=k+1
iter[k]=k

@inbounds for ic=1:oper.nc
    p[ic]=x[ic]
end
saxpy(false,false,oper,p,q)
for ir=1:oper.nr
    s[ir]=b[ir]-q[ir]
end
saxpy(true,false,oper,r,s)

for j=1:maxit2
for ic=1:oper.nc
    r[ic]=r[ic]-bshift*(p[ic]-v[ic]+u[ic])
end 
saxpy(false,false,oper,r,q)
mu=dot(r,r)/(dot(q,q) + bshift*dot(r,r))
for ic=1:oper.nc
    p[ic]=p[ic]+mu*r[ic]
end 
for ir=1:oper.nr
    s[ir]=s[ir]-mu*q[ir]
end
saxpy(true,false,oper,r,s)
# projection
for ic=1:oper.nc
    p[ic]=min(max(p[ic],0.2),5.0)
end 
end

@inbounds for ic=1:oper.nc
    x[ic]=p[ic]
end

@inbounds for ic=1:oper.nc
    p[ic]=v[ic]
end

# over-relaxation
@inbounds for ic=1:oper.nc
    xhat[ic]=alpha*x[ic] - (1.0 - alpha)*p[ic]
end

# objective function
rnorm=norm(r,2)
denoiser(grd,p,r)
obj[k]=0.5*rnorm*rnorm + 0.5*lambda*(dot(p,p .- r))

# v-update 
for j=1:maxit3
denoiser(grd,p,r)
@inbounds for ic=1:oper.nc
    p[ic]=(lshift*r[ic]+bshift*(xhat[ic]+u[ic]))/(bshift+lshift)
end
end

# primal and dual residue 
println("norm(p-v) = ",norm(p .- v,2))
println("norm(x-v) = ",norm(x .- v,2))
rpri[k]=norm(x .- v,2)
sdual[k]=norm(beta.*(p .- v,2))

# varying penalty parameter
if mod(k,30) == 0
if rpri[k] > murho*sdual[k]
beta = tau*beta
u .= u ./ tau
end

if sdual[k] > murho*rpri[k]
beta = beta/tau
u .= u .* tau
end
bshift=beta*normb
end
rhotil[k]=beta

# u-update
@inbounds for ic=1:oper.nc
u[ic] += x[ic]-p[ic]
end

@inbounds for ic=1:oper.nc
    v[ic]=p[ic]
end

epsdual[k] = ABSTOL*nroot + RELTOL * norm(beta.*u,2)
epspri[k] = ABSTOL*nroot + RELTOL * max(norm(x,2),norm(v,2))

if (rpri[k] <= epspri[k]) && (sdual[k] <= epsdual[k])
   println("dual")
   break
end

end # while

println("k = ",k)
println("norm(b) = ",norm(b,2))
println("norm(s) = ",norm(s,2))
println("norm(r) = ",rnorm)
println("norm(x-v) = ",norm( x .- v,2))
println("norm(x) = ",norm(x,2))
println("norm(v) = ",norm(v,2))

graph(iter,obj,nt="Objective function",nx="k",ny="1/2|| Ax -b ||_2^2 + ld/2 < x, x-f(x)> ",sav=true)
graph(iter,rpri,nt="Primal residual",nx="k",ny="|| r ||_2",sav=true)
graph(iter,sdual,pl=false,nt="Dual residual",nx="k",ny="|| s ||_2",sav=true)
graph(iter,epspri,nt="Primal tolerance",nx="k",ny="eps^primal",sav=true)
graph(iter,epsdual,pl=false,nt="Dual tolerance",nx="k",ny="eps^dual",sav=true)
graph(iter,rhotil,nt="Penality term",nx="k",ny="rho",sav=true)

return nothing
end # red_admm



function cglsp!(oper::Topr,lambda::Float64,sig2::Float64,
    x::Array{Float64,1},b::Array{Float64,1},
    x0::Array{Float64,1},xp::Array{Float64,1},saxpy::Function,maxit::Int64)

@assert( size(x,1)  == oper.nc )
@assert( size(x0,1) == oper.nc )
@assert( size(b,1) == oper.nr )

p::Array{Float64,1}=zeros(Float64,oper.nc)
q::Array{Float64,1}=zeros(Float64,oper.nr)
r::Array{Float64,1}=zeros(Float64,oper.nr)
s::Array{Float64,1}=zeros(Float64,oper.nc)

@inbounds for ic=1:oper.nc
x[ic]=x0[ic]
#x[ic]=0.0
end

saxpy(false,false,oper,x,q)
@inbounds for ir=1:oper.nr
r[ir]=b[ir]-q[ir]
end

saxpy(true,false,oper,p,r)
@inbounds for ic=1:oper.nc
s[ic]=p[ic]-lambda*x[ic]
end

# Initialize
@inbounds for ic=1:oper.nc
p[ic]=s[ic]
end

norms0::Float64=norm(s,2)
gamma::Float64=norms0*norms0
normx::Float64=norm(x,2)
xmax::Float64=normx
alpha::Float64=0.0
delta::Float64=0.0

k::Int64=0
flag::Int32=0

while ( k < maxit ) && ( flag == 0 )

k=k+1
saxpy(false,false,oper,p,q)
delta=norm(q,2)^2+lambda*norm(p,2)^2
if delta <= 0.0
flag=3
end

if delta == 0.0
delta=TOL
end
alpha=gamma/delta

@inbounds for ic=1:oper.nc
x[ic] = x[ic] + alpha * p[ic]
end

@inbounds for ir=1:oper.nr
r[ir] = r[ir] - alpha * q[ir]
end

saxpy(true,false,oper,s,r)
@inbounds for ic=1:oper.nc
s[ic] = s[ic] - lambda * (x[ic]-xp[ic])
end
norms = norm(s,2)
gamma1= gamma
gamma = norms*norms
beta  = gamma/gamma1

@inbounds for ic=1:oper.nc
p[ic] = s[ic] + beta * p[ic]
end

normx=norm(x,2)
xmax=max(xmax,normx)

if norms < TOL*norms0
flag=1
end

#println("norms0 = ",norms0);
#println("norms = ",norms);
#println("normx = ",normx);
#println("normp = ",alpha*norm(p,2));
end #while

if k == maxit
flag=2
end

if abs(alpha)*norm(p,2) <= sqrt(TOL)*xmax
flag=3
end

println("flag = ",flag);

return flag

end # cglsp



function lasso_admm!(grd::Tgrd,oper::Topr,
    lambda::Float64,sig2::Float64,beta::Float64,
    x::Array{Float64,1},b::Array{Float64,1},
    x0::Array{Float64,1},saxpy::Function,shrink::Function,
    maxit1::Int64,maxit2::Int64)

u::Array{Float64,1}=zeros(Float64,oper.nc)
v::Array{Float64,1}=zeros(Float64,oper.nc)
vold::Array{Float64,1}=zeros(Float64,oper.nc)
xhat::Array{Float64,1}=zeros(Float64,oper.nc)
xprox::Array{Float64,1}=zeros(Float64,oper.nc)

p::Array{Float64,1}=zeros(Float64,oper.nc)
q::Array{Float64,1}=zeros(Float64,oper.nr)
r::Array{Float64,1}=zeros(Float64,oper.nc)
s::Array{Float64,1}=zeros(Float64,oper.nr)
rnorm::Float64=0.0
alpha::Float64=1.0
rho::Float64=beta
tau::Float64=2.0
murho::Float64=10.0
lmbd::Float64=1.0
nroot::Float64=sqrt(Float64(oper.nc))

rho=beta
lmbd=lambda;

objls::Array{Float64,1}=zeros(Float64,maxit1)
objr::Array{Float64,1}=zeros(Float64,maxit1)
obj::Array{Float64,1}=zeros(Float64,maxit1)
rpri::Array{Float64,1}=zeros(Float64,maxit1)
sdual::Array{Float64,1}=zeros(Float64,maxit1)
epspri::Array{Float64,1}=zeros(Float64,maxit1)
epsdual::Array{Float64,1}=zeros(Float64,maxit1)
rhotil::Array{Float64,1}=zeros(Float64,maxit1)
iter::Array{Float64,1}=zeros(Float64,maxit1)

# Initialize
@inbounds for ic=1:oper.nc
x[ic]=x0[ic]
end

@inbounds for ic=1:oper.nc
v[ic]=0.0
end

@inbounds for ic=1:oper.nc
u[ic]=0.0
end


# Smooth objective function minimization
k=0
while k < maxit1

k=k+1
iter[k]=k

@inbounds for ic=1:oper.nc
p[ic]=x[ic]
end

@inbounds for ic=1:oper.nc
    xprox[ic]=v[ic]-u[ic]
end
cglsp!(oper,rho,sig2,p,b,x0,xprox,saxpy!,maxit2)

# projection on convex set
for ic=1:oper.nc
 p[ic]=min(max(p[ic],0.2),5.0)
end

@inbounds for ic=1:oper.nc
x[ic]=p[ic]
end

saxpy(false,false,oper,p,q)
for ir=1:oper.nr
s[ir]=b[ir]-q[ir]
end

# over-relaxation
@inbounds for ic=1:oper.nc
vold[ic]=v[ic]
xhat[ic]=alpha*x[ic]-(1.0-alpha)*vold[ic]
end

# L1 minimization
@inbounds for ic=1:oper.nc
p[ic]=xhat[ic]+u[ic]
end
v = shrink(p,lmbd/rho)


@inbounds for ic=1:oper.nc
u[ic] += x[ic]-v[ic]
end


rnorm=norm(s,2)

objls[k]=0.5*rnorm*rnorm
objr[k]=lmbd*norm(v,1)

obj[k]=0.5*rnorm*rnorm+lmbd*norm(v,1)
rpri[k]=norm(x .- v,2)
sdual[k]=norm(rho.*(v .- vold,2))

#
# varying penalty parameter
#
if mod(k,50) == 0
if rpri[k] > murho*sdual[k]
rho = tau*rho
u .= u ./ tau
end

if sdual[k] > murho*rpri[k]
rho = rho/tau
u .= u .* tau
end
end
#
rhotil[k]=rho

println("lsq = ",0.5*rnorm*rnorm)
println("obj = ",obj[k])
println("rpri = ",rpri[k])
println("sdual = ",sdual[k])

# stopping
epsdual[k] = ABSTOL*nroot + RELTOL * norm(rho.*u,2)
epspri[k] = ABSTOL*nroot + RELTOL * max(norm(x,2),norm(v,2))
if (rpri[k] <= epspri[k]) && (sdual[k] <= epsdual[k])
   println("dual")
   break
end

end # while

println("k = ",k)
println("norm(b) = ",norm(b,2))
println("norm(s) = ",norm(s,2))
println("norm(r) = ",rnorm)
println("norm(x-v) = ",norm( x .- v,2))
println("norm1(x) = ",norm(x,1))
println("norm(x) = ",norm(x,2))
println("norm(v) = ",norm(v,2))
println("lmbd = ",lmbd)
println("rho  = ",rho)
println("oper.norm  = ",oper.norm)

graph(iter,objls,nt="Least Square ",nx="k",ny="1/2|| Ax -b ||_2^2",sav=true)
graph(iter,objr,nt="Regularization",nx="k",ny="ld|| x ||_1",sav=true)
graph(iter,obj,nt="Objective function",nx="k",ny="1/2|| Ax -b ||_2^2 + ld|| x ||_1",sav=true)
graph(iter,rpri,nt="Primal residual",nx="k",ny="|| r ||_2",sav=true)
graph(iter,sdual,nt="Dual residual",nx="k",ny="|| s ||_2",sav=true)
graph(iter,epspri,pl=false,nt="Primal tolerance",nx="k",ny="eps^primal",sav=true)
graph(iter,epsdual,pl=false,nt="Dual tolerance",nx="k",ny="eps^dual",sav=true)
graph(iter,rhotil,nt="Penality term",nx="k",ny="rho",sav=true)


return nothing
end # lasso_admm


function lasso_ista!(oper::Topr,
    sig2::Float64,       # level of noise;
    Lf::Float64,         # Lipschitz constant;
    lambda::Float64,     # penalty parameter;
    x::Array{Float64,1}, # the output model;
    b::Array{Float64,1}, # the observed data;
    x0::Array{Float64,1},# the input model;  
    shrink::Function,    # soft thresholding
    saxpy::Function,     # forward and backward solves;
    maxit::Int64)        # number of outer iterations;

    @assert( length(x)  == oper.nc )
    @assert( length(x0) == oper.nc )
    @assert( length(b) == oper.nr )
    
    p::Array{Float64,1}=zeros(Float64,oper.nc)
    q::Array{Float64,1}=zeros(Float64,oper.nr)
    r::Array{Float64,1}=zeros(Float64,oper.nr)
    s::Array{Float64,1}=zeros(Float64,oper.nc)
    
    Tol::Float64 = 1.0e-05
    lmbd::Float64=lambda
    alpha::Float64=0.0

    objls::Array{Float64,1}=zeros(Float64,maxit)
    objr::Array{Float64,1}=zeros(Float64,maxit)
    obj::Array{Float64,1}=zeros(Float64,maxit)
    iter::Array{Float64,1}=zeros(Float64,maxit)

    # intialize #
    for ic=1:oper.nc
        x[ic] = x0[ic]
    end
    
    # step-size #
    alpha=1.0e-02/(2.0*Lf)
    
    k = 0
    while k < maxit
       
        k = k + 1
        iter[k]=k

        @inbounds for ic=1:oper.nc
            x0[ic] = x[ic]
        end

        # data residual #
        saxpy(false,false,oper,x,q)
        @inbounds for ir=1:oper.nr
            r[ir] = b[ir]-q[ir]
        end

        # steepest descent #
        saxpy(true,false,oper,s,r)        
        @inbounds for ic = 1:oper.nc
            x[ic] = x[ic] + alpha*s[ic]
        end
    
        # soft thresholding #
        p = shrink(x,lmbd*alpha)

        # objective function #
        objls[k]=0.5*norm(r,2)*norm(r,2)
        objr[k]=lmbd*norm(p,1)

        obj[k]=0.5*norm(r,2)*norm(r,2) + lmbd*norm(p,1)

        # stopping criterion #
        if norm(p - x0,2)/(1.0 + norm(x0,2))  < Tol
            println("optimal")
            break
        end

        @inbounds for ic=1:oper.nc
            x[ic]=p[ic]
        end
        println("lsq = ",0.5*norm(r,2)*norm(r,2) )
        println("norm(|x - x0|) = ", norm(x - x0,2) ) 
       
    end

    println("k = ",k)
    println("lmbd = ",lmbd)

    graph(iter,objls,nt="Least Square ",nx="k",ny="1/2|| Ax -b ||_2^2",sav=true)
    graph(iter,objr,nt="Regularization",nx="k",ny="ld|| x ||_1",sav=true)
    graph(iter,obj,nt="Objective function",nx="k",ny="1/2|| Ax -b ||_2^2 + ld|| x ||_1",sav=true)

return nothing
end#lasso_ista


function lasso_fista!(oper::Topr,
    sig2::Float64,       # level of noise;
    Lf::Float64,         # Lipschitz constant;
    lambda::Float64,     # penalty parameter;
    x::Array{Float64,1}, # the output model;
    b::Array{Float64,1}, # the observed data;
    x0::Array{Float64,1},# the input model;  
    shrink::Function,    # soft thresholding
    saxpy::Function,     # forward and backward solves;
    maxit::Int64)        # number of outer iterations;

    @assert( length(x)  == oper.nc )
    @assert( length(x0) == oper.nc )
    @assert( length(b) == oper.nr )
    
    p::Array{Float64,1}=zeros(Float64,oper.nc)
    q::Array{Float64,1}=zeros(Float64,oper.nr)
    r::Array{Float64,1}=zeros(Float64,oper.nr)
    s::Array{Float64,1}=zeros(Float64,oper.nc)
    y::Array{Float64,1}=zeros(Float64,oper.nc)
    
    Tol::Float64 = 5.0e-04
    lmbd::Float64=lambda
    alpha::Float64=0.0
    theta::Float64=0.0
    thetaold::Float64=0.0

    objls::Array{Float64,1}=zeros(Float64,maxit)
    objr::Array{Float64,1}=zeros(Float64,maxit)
    obj::Array{Float64,1}=zeros(Float64,maxit)
    iter::Array{Float64,1}=zeros(Float64,maxit)

    # intialize #
    for ic=1:oper.nc
        x[ic] = x0[ic]
    end

    for ic=1:oper.nc
        y[ic] = 0.0 
    end
    
    # step-size #
    alpha=1.0e-02/(2.0*Lf)
    
    k = 0
    while k < maxit
       
        k = k + 1
        iter[k]=k

        @inbounds for ic=1:oper.nc
            x0[ic] = x[ic]
        end

        # data residual #
        saxpy(false,false,oper,y,q)
        @inbounds for ir=1:oper.nr
            r[ir] = b[ir]-q[ir]
        end

        # steepest descent #
        saxpy(true,false,oper,s,r)        
        @inbounds for ic = 1:oper.nc
            x[ic] = y[ic] + alpha*s[ic]
        end
    
        # soft thresholding #
        p = shrink(x,lmbd*alpha)

        thetaold = theta
        theta = (1.0 + sqrt(1.0+4.0*theta*theta) )/2.0

        # update the extrapolation point #
        @inbounds for ic = 1:oper.nc
            y[ic] = p[ic] + ( (thetaold -1.0)/theta ) * (p[ic] - x0[ic])
        end

        # objective function #
        objls[k]=0.5*norm(r,2)*norm(r,2)
        objr[k]=lmbd*norm(p,1)

        obj[k]=0.5*norm(r,2)*norm(r,2) + lmbd*norm(p,1)

        # stopping criterion #
        if norm(p - x0,2)/(1.0 + norm(x0,2))  < Tol
            println("optimal")
            break
        end

        @inbounds for ic=1:oper.nc
            x[ic]=p[ic]
        end
        println("lsq = ",0.5*norm(r,2)*norm(r,2) )
        println("norm(|x - x0|) = ", norm(x - x0,2) ) 
       
    end

    println("k = ",k)
    println("lmbd = ",lmbd)

    graph(iter,objls,nt="Least Square ",nx="k",ny="1/2|| Ax -b ||_2^2",sav=true)
    graph(iter,objr,nt="Regularization",nx="k",ny="ld|| x ||_1",sav=true)
    graph(iter,obj,nt="Objective function",nx="k",ny="1/2|| Ax -b ||_2^2 + ld|| x ||_1",sav=true)

return nothing
end#lasso_fista


function display(image::Array{Float64,1},
    n1::Int64,n2::Int64,
    o1::Float64,o2::Float64,
    d1::Float64,d2::Float64, 
    name::String="",sav::Bool=false)

bmxm = maximum(reshape(image, (n1, n2)))
bmnm = minimum(reshape(image, (n1, n2)))
bmed = (bmnm + bmxm) / 2
barticks = [bmnm, (bmnm + bmed) / 2, bmed, (bmxm + bmed) / 2, bmxm]

e1::Int64=mod(n1,2)
e2::Int64=mod(n2,2)

fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(6, 7))
im = ax.imshow(reshape(image, (n1, n2)),
cmap=plt.cm.jet,
origin="upper",
interpolation="None",
extent=[o2,(n2 - e2) * d2, (n1 - e1) * d1, o1])

ax.set_title(label= name ,loc="center",fontsize=18,fontweight="bold")
ax.set_xlabel("Distance (m)", fontsize=13, fontweight="bold")
ax.set_ylabel("Depth (m)", fontsize=13, fontweight="bold")  
ax.locator_params(axis="x", nbins=5)
ax.axis("image")

cbar = fig.colorbar(im,
ax=ax,
orientation="vertical",
spacing="uniform",
shrink=0.94,
ticks=barticks,
fraction=0.046*(n1/n2), pad=0.04)
cbar.set_label("Slowless (s/m)",fontsize=13,fontweight="bold")

fig.show()
if sav == true
    fig.savefig(join([name,".pdf"]))
end

return nothing    
end#display


function graph(x::Array{Float64,1},y::Array{Float64,1};pl::Bool=false,
    nt::String="",nx::String="x",ny::String="y",sav::Bool=false)

k=0
for ic=1:length(x)
    if x[ic] != 0.0
        k+=1
    end
end

fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(9, 5))
if pl==false
ax.plot(x[1:k],y[1:k],linewidth=2, color="b")
else
#yn::Array{Float64,1}=zeros(Float64,k)
#for i=1:k
#    yn[i]=max(1.e-8,y[i])
#end
ax.semilogy(x[1:k],y[1:k],linewidth=2, color="b")
#ax.loglog(x[1:k],y[1:k],linewidth=2, color="b")
end

ax.set_title(label=nt ,loc="center",fontsize=18,fontweight="bold")
ax.set_xlabel(nx, fontsize=13, fontweight="bold")
ax.set_ylabel(ny, fontsize=13, fontweight="bold")  

ax.set_xlim(minimum(x[1:k]), maximum(x[1:k]))
#ax.set_ylim(minimum(y[1:k]), maximum(y[1:k]))    

fig.show()
if sav == true
    fig.savefig(join([nt,".pdf"]))
end

return nothing
end#graph



function display(image::Array{Float64,1},
    xs::Array{Float64,1},
    ys::Array{Float64,1},
    xr::Array{Float64,1},
    yr::Array{Float64,1},
    n1::Int64,n2::Int64,
    o1::Float64,o2::Float64,
    d1::Float64,d2::Float64, 
    name::String="",sav::Bool=false)

bmxm = maximum(reshape(image, (n1, n2)))
bmnm = minimum(reshape(image, (n1, n2)))
bmed = (bmnm + bmxm) / 2
barticks = [bmnm, (bmnm + bmed) / 2, bmed, (bmxm + bmed) / 2, bmxm]
 
e1::Int64=mod(n1,2)
e2::Int64=mod(n2,2)

fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(6,7))
im = ax.imshow(reshape(image, (n1, n2)),
cmap=plt.cm.jet,
origin="upper",
interpolation="None",
extent=[o2,(n2 - e2) * d2, (n1 - e1) * d1, o1])

ax.set_title(label= name ,loc="center",fontsize=18,fontweight="bold")
ax.set_xlabel("Distance (m)", fontsize=13, fontweight="bold")
ax.set_ylabel("Depth (m)", fontsize=13, fontweight="bold")  
ax.plot(xs,ys,"ro")
ax.plot(xr,yr,"wo")

cbar = fig.colorbar(im,
ax=ax,
orientation="vertical",
spacing="uniform",
shrink=0.94,
ticks=barticks,
fraction=0.046*(n1/n2), pad=0.04)
cbar.set_label("Slowless (s/m)",fontsize=13,fontweight="bold")

fig.show()
if sav == true
    fig.savefig(join([name,".pdf"]))
end

return nothing 
end#display


function acqg!(grd::Tgrd,dat::Tdat,
    ns::Int64,nr::Int64;tg::String="cross-well")

xs::Array{Float64,1} = zeros(Float64, ndim)
xr::Array{Float64,1} = zeros(Float64, ndim)
nde::Int64=0
nse::Int64=0

if tg == "cross-well"

    id = 0
    @inbounds for is = 1:ns
        xs[1] = grd.o1 + (Float64(is) - 0.5) * grd.d1
        xs[2] = grd.o2
        for ir = 1:nr
            xr[1] = grd.o1 + (Float64(ir) - 0.5) * grd.d1
            xr[2] = grd.o2 + grd.n2 * grd.d2
            id += 1
            dat.xs[1,id] = xs[1]
            dat.xs[2,id] = xs[2]
            dat.xr[1,id] = xr[1]
            dat.xr[2,id] = xr[2]
        end
    end

elseif tg == "cw-top" 
    
    nde=grd.n1*grd.n2
    id=0
    @inbounds for is = 1:ns
        xs[1] = grd.o1 + (Float64(1) - 0.5) * grd.d1
        xs[2] = grd.o2
        for ir = 1:grd.n2 
            xr[1] = grd.o1
            xr[2] = grd.o2 + (Float64(ir) - 0.5) * grd.d2
            id += 1
            dat.xs[1,id] = xs[1]
            dat.xs[2,id] = xs[2]
            dat.xr[1,id] = xr[1]
            dat.xr[2,id] = xr[2]
        end
    end
    
    id = 0
    @inbounds for is = 1:ns    
        xs[1] = grd.o1 + (Float64(is) - 0.5) * grd.d1 
        xs[2] = grd.o2
        for ir = 1:nr
            xr[1] = grd.o1 + (Float64(ir) - 0.5) * grd.d1
            xr[2] = grd.o2 + grd.n2 * grd.d2
            id += 1
            dat.xs[1,id+nde] = xs[1]
            dat.xs[2,id+nde] = xs[2]
            dat.xr[1,id+nde] = xr[1]
            dat.xr[2,id+nde] = xr[2]
        end
    end    

elseif tg == "cw-top-botton"    

    nde=grd.n1*grd.n2
    nse=grd.n1*grd.n1
    id=0
    @inbounds for is = 1:ns
        xs[1] = grd.o1 + (Float64(1) - 0.5) * grd.d1
        xs[2] = grd.o2
        for ir = 1:grd.n2 
            xr[1] = grd.o1
            xr[2] = grd.o2 + (Float64(ir) - 0.5) * grd.d2
            id += 1
            dat.xs[1,id] = xs[1]
            dat.xs[2,id] = xs[2]
            dat.xr[1,id] = xr[1]
            dat.xr[2,id] = xr[2]
        end
    end
    
    id = 0
    @inbounds for is = 1:ns    
        xs[1] = grd.o1 + (Float64(is) - 0.5) * grd.d1 
        xs[2] = grd.o2
        for ir = 1:nr
            xr[1] = grd.o1 + (Float64(ir) - 0.5) * grd.d1
            xr[2] = grd.o2 + grd.n2 * grd.d2
            id += 1
            dat.xs[1,id+nde] = xs[1]
            dat.xs[2,id+nde] = xs[2]
            dat.xr[1,id+nde] = xr[1]
            dat.xr[2,id+nde] = xr[2]
        end
    end
    
    id = 0
    @inbounds for is = 1:ns    
        xs[1] = (grd.o1 + grd.n1 * grd.d1) - (Float64(1) - 0.5) * grd.d1
        xs[2] = grd.o2
        for ir = 1:grd.n2
            xr[1] = grd.o1 + grd.n1 * grd.d1
            xr[2] = (grd.o2+grd.n2*grd.d2)  - (Float64(ir) - 0.5) * grd.d2
            id += 1
            dat.xs[1,id+nde+nse] = xs[1]
            dat.xs[2,id+nde+nse] = xs[2]
            dat.xr[1,id+nde+nse] = xr[1]
            dat.xr[2,id+nde+nse] = xr[2]
        end
    end

elseif tg == "full"     

nde=grd.n1*grd.n2
nse=grd.n1*grd.n1
nce=grd.n2*grd.n2
nre=nse + 2*nde
nxe=nre + nse + nce
npe=nxe + 2*nse

#source left
id=0
@inbounds for is = 1:grd.n1
    xs[1] = (grd.o1+grd.n1*grd.d1) - (Float64(1) - 0.5) * grd.d1
    xs[2] = grd.o2
    for ir = 1:grd.n2 
        xr[1] = grd.o1+grd.n1*grd.d1
        xr[2] = grd.o2 + (Float64(ir) - 0.5) * grd.d2
        id += 1
        dat.xs[1,id] = xs[1]
        dat.xs[2,id] = xs[2]
        dat.xr[1,id] = xr[1]
        dat.xr[2,id] = xr[2]
    end
end

id = 0
@inbounds for is = 1:grd.n1    
    xs[1] = (grd.o1+grd.n1*grd.d1) - (Float64(is) - 0.5) * grd.d1 
    xs[2] = grd.o2
    for ir = 1:grd.n1
        xr[1] = (grd.o1+grd.n1*grd.d1) - (Float64(ir) - 0.5) * grd.d1
        xr[2] = grd.o2 + grd.n2 * grd.d2
        id += 1
        dat.xs[1,id+nde] = xs[1]
        dat.xs[2,id+nde] = xs[2]
        dat.xr[1,id+nde] = xr[1]
        dat.xr[2,id+nde] = xr[2]
    end
end

id = 0
@inbounds for is = 1:grd.n1  
    xs[1] = grd.o1 + (Float64(1) - 0.5) * grd.d1
    xs[2] = grd.o2
    for ir = 1:grd.n2
        xr[1] = grd.o1
        xr[2] = (grd.o2+grd.n2*grd.d2) - (Float64(ir) - 0.5) * grd.d2
        id += 1
        dat.xs[1,id+nde+nse] = xs[1]
        dat.xs[2,id+nde+nse] = xs[2]
        dat.xr[1,id+nde+nse] = xr[1]
        dat.xr[2,id+nde+nse] = xr[2]
    end
end

#source botton 
id=0
@inbounds for is = 1:grd.n2
    xs[1] = grd.o1+grd.n1*grd.d1
    xs[2] = (grd.o2+grd.n2*grd.d2) -  (Float64(1) - 0.5) * grd.d2
    for ir = 1:grd.n1
        xr[1] = (grd.o1+grd.n1*grd.d1)  - (Float64(ir) - 0.5) * grd.d1
        xr[2] = grd.o2+grd.n2*grd.d2
        id += 1
        dat.xs[1,id+nre] = xs[1]
        dat.xs[2,id+nre] = xs[2]
        dat.xr[1,id+nre] = xr[1]
        dat.xr[2,id+nre] = xr[2]
    end
end

id = 0
@inbounds for is = 1:grd.n2
    xs[1] = grd.o1+grd.n1*grd.d1
    xs[2] = (grd.o2+grd.n2*grd.d2) - (Float64(is) - 0.5) * grd.d2
    for ir = 1:grd.n2
        xr[1] = grd.o1
        xr[2] = (grd.o2+grd.n2*grd.d2) - (Float64(ir) - 0.5) * grd.d2
        id += 1
        dat.xs[1,id+nre+nde] = xs[1]
        dat.xs[2,id+nre+nde] = xs[2]
        dat.xr[1,id+nre+nde] = xr[1]
        dat.xr[2,id+nre+nde] = xr[2]
    end
end

id = 0
@inbounds for is = 1:grd.n2 
    xs[1] = grd.o1+grd.n1*grd.d1
    xs[2] = grd.o2 + (Float64(1) - 0.5) * grd.d2
    for ir = 1:grd.n1
        xr[1] = grd.o1 + (Float64(ir) - 0.5) * grd.d1
        xr[2] = grd.o2
        id += 1
        dat.xs[1,id+nre+nde+nce] = xs[1]
        dat.xs[2,id+nre+nde+nce] = xs[2]
        dat.xr[1,id+nre+nde+nce] = xr[1]
        dat.xr[2,id+nre+nde+nce] = xr[2]        
    end
end    

# source right
id=0
@inbounds for is = 1:grd.n1
    xs[1] = grd.o1 + (Float64(1) - 0.5) * grd.d1
    xs[2] = grd.o2 + grd.n2 * grd.d2
    for ir = 1:grd.n2 
        xr[1] = grd.o1
        xr[2] = (grd.o2 + grd.n2 * grd.d2) - (Float64(ir) - 0.5) * grd.d2
        id += 1
        dat.xs[1,id+nxe] = xs[1]
        dat.xs[2,id+nxe] = xs[2]
        dat.xr[1,id+nxe] = xr[1]
        dat.xr[2,id+nxe] = xr[2]
    end
end

id = 0
@inbounds for is = 1:grd.n1
    xs[1] = grd.o1 + (Float64(is) - 0.5) * grd.d1 
    xs[2] = grd.o2 + grd.n2 * grd.d2
    for ir = 1:grd.n1
        xr[1] = grd.o1 + (Float64(ir) - 0.5) * grd.d1
        xr[2] = grd.o2
        id += 1
        dat.xs[1,id+nxe+nde] = xs[1]
        dat.xs[2,id+nxe+nde] = xs[2]
        dat.xr[1,id+nxe+nde] = xr[1]
        dat.xr[2,id+nxe+nde] = xr[2]
    end
end

id = 0
@inbounds for is = 1:grd.n1   
    xs[1] = (grd.o1 + grd.n1 * grd.d1) - (Float64(1) - 0.5) * grd.d1
    xs[2] = grd.o2+grd.n2*grd.d2
    for ir = 1:grd.n2
        xr[1] = grd.o1+grd.n1*grd.d1
        xr[2] = grd.o2 + (Float64(ir) - 0.5) * grd.d2
        id += 1
        dat.xs[1,id+nxe+nde+nse] = xs[1]
        dat.xs[2,id+nxe+nde+nse] = xs[2]
        dat.xr[1,id+nxe+nde+nse] = xr[1]
        dat.xr[2,id+nxe+nde+nse] = xr[2]
    end
end

#source top
id=0
@inbounds for is = 1:grd.n2
    xs[1] = grd.o1 
    xs[2] = grd.o2 + (Float64(1) - 0.5) * grd.d2
    for ir = 1:grd.n1
        xr[1] = grd.o1  + (Float64(ir) - 0.5) * grd.d1
        xr[2] = grd.o2
        id += 1
        dat.xs[1,id+npe] = xs[1]
        dat.xs[2,id+npe] = xs[2]
        dat.xr[1,id+npe] = xr[1]
        dat.xr[2,id+npe] = xr[2]
    end
end

id = 0
@inbounds for is = 1:grd.n2
    xs[1] = grd.o1
    xs[2] = grd.o2 +  (Float64(is) - 0.5) * grd.d2
    for ir = 1:grd.n2
        xr[1] = grd.o1 + grd.n1 * grd.d1
        xr[2] = grd.o2 + (Float64(ir) - 0.5) * grd.d2
        id += 1
        dat.xs[1,id+npe+nde] = xs[1]
        dat.xs[2,id+npe+nde] = xs[2]
        dat.xr[1,id+npe+nde] = xr[1]
        dat.xr[2,id+npe+nde] = xr[2]
    end
end

id = 0
@inbounds for is = 1:grd.n2  
    xs[1] = grd.o1
    xs[2] = (grd.o2 + grd.n2 * grd.d2) - (Float64(1) - 0.5) * grd.d2
    for ir = 1:grd.n1
        xr[1] = (grd.o1+grd.n1*grd.d1)  - (Float64(ir) - 0.5) * grd.d1
        xr[2] = grd.o2 + grd.n2 * grd.d2 
        id += 1
        dat.xs[1,id+npe+nde+nce] = xs[1]
        dat.xs[2,id+npe+nde+nce] = xs[2]
        dat.xr[1,id+npe+nde+nce] = xr[1]
        dat.xr[2,id+npe+nde+nce] = xr[2]        
    end
end

end#

return nothing
end#acqg



function Main()
#=
Test REgularization by Denoise using straight ray tomography 
=#

n1=20
n2=10
o1=0.0
o2=0.0
d1=5.0
d2=5.0

#cross well#
ns=n1
nr=n1
nd::Int64=ns*nr

#cross well and top#
#nd=nd+n1*n2

#cross well and, top and botton#
#nd=nd+2*(n1*n2)

#full
nd=nd+2*n1*n2     #until right
nd=nd+n1*n1+n2*n2 #until botton
nd=nd+2*n1*n1     #until left
nd=nd+n1*n1+n2*n2 #until top

grid::Tgrd=Tgrd(n1,n2,o1,o2,d1,d2)
dat::Tdat=Tdat(nd)
nrows::Int64=dat.nd
ncols::Int64=grid.n1*grid.n2
size::Int64=dat.nd*(grid.n1+grid.n2)

oper::Topr=Topr(nrows,ncols,size)

xs::Array{Float64,1}=zeros(Float64,ndim)
xr::Array{Float64,1}=zeros(Float64,ndim)

#***********************************************

acqg!(grid,dat,ns,nr,tg="full")

#***********************************************


raytracing!(grid,dat,oper)

mdl::Array{Float64,1}=zeros(Float64,n1*n2)
dfit::Array{Float64,1}=zeros(Float64,dat.nd)
rsd::Array{Float64,1}=zeros(Float64,dat.nd)

test_model!(grid.n1,grid.n2,mdl)

saxpy!(false,false,oper,mdl,dat.od)

perc=0.01
addnoise!(dat,perc)


#ext=heatmap(reshape(mdl,(n1,n2)),yflip=true,color=:jet1,clims=(0,2))
plt.close("all")
display(mdl, n1, n2, o1, o2, d1, d2, "true-model", true)


display(mdl, dat.xs[2,:], dat.xs[1,:],
                dat.xr[2,:], dat.xr[1,:],
grid.n1, grid.n2, grid.o1, grid.o2, grid.d1, grid.d2, "acq-geom", true)


# RED inversion
xls::Array{Float64,1}=zeros(Float64,n1*n2)
x0::Array{Float64,1}=zeros(Float64,n1*n2)

avg=0.0
for id=1:dat.nd
    xs .= dat.xs[1:2,id]
    xr .= dat.xr[1:2,id]
    avg += dat.od[id]/norm(xr .- xs,2) 
end
avg=avg/dat.nd

# Least-squares solution #
lambda=0.01
sig2=1.0
maxit1=oper.nc
x0 .= avg

cgls!(oper,lambda*oper.norm,sig2,xls,dat.od,x0,saxpy!,maxit1)

#lsq=heatmap(reshape(xls,(n1,n2)),yflip=true,color=:jet1,clims=(0,2))
display(xls, n1, n2, o1, o2, d1, d2, "cgls-model", true)

saxpy!(false, false, oper, xls, dfit)
for id = 1:dat.nd
    rsd[id] = dat.od[id] - dfit[id]
end
println("data norm: ", norm(dat.od, 2), " residual norm: ", norm(rsd, 2))

x0 .= mdl .- xls
mfit = 100.0*( 1.0 - norm(x0,2)/norm(mdl,2))
println("model fit =", mfit)

display(x0, n1, n2, o1, o2, d1, d2, "cgls-model-res", true)
println("mod res norm: ", norm(x0,2))
x0 .= 0.0

#---------------
# RED inversion
#---------------

# choose redtype: "sd"; "fp"; "admm"; "ista"; "fista"; "lasso"
#
redtype = "fp"
if redtype == "sd"

    println(">> RED-SD")
    lambda=0.0001 #0.1
    sig2=1.0
    maxit=5000
    x0 .= avg

    lambda_scl=0.0
    tol=1.0e-06
    maxitn=1000
    lambda_scl = normestfb(oper,saxpy!,tol=tol,maxit=maxitn)
    lambda=lambda*lambda_scl

    red_sd!(grid,oper,lambda,sig2,
        xls,dat.od,x0,
        saxpy!,denoiser!,maxit)

elseif redtype == "fp"

    println(">> RED-FP")
    lambda=0.000001 #0.1
    sig2=1.0
    maxit1=1000
    maxit2=20
    x0 .= avg

    lambda_scl=0.0
    tol=1.0e-06
    maxitn=1000
    lambda_scl = normestfb(oper,saxpy!,tol=tol,maxit=maxitn)
    lambda=lambda*lambda_scl

    red_fp!(grid,oper,lambda,sig2,
        xls,dat.od,x0,
        saxpy!,denoiser!,
        maxit1,maxit2) 

elseif redtype == "ista"
    
    println(">> LASSO-ISTA")              
    sig2=1.0
    maxit1=oper.nc
    maxit2=oper.nc
    maxit3=10
    lambda=0.01
    beta=0.001  
    x0 .= avg

    saxpy!(false,false,oper,x0,dfit)

    dfit .= dat.od .- dfit 
    rsdnorm=norm(dfit,2)^2
    x0norm1=norm(x0,1)
    lambda_scl=rsdnorm/x0norm1
    lambda=lambda*lambda_scl

    Lip=0.0
    tol=1.0e-06
    maxitn=1000
    Lip = normest(oper,saxpy!,tol=tol,maxit=maxitn)

    lasso_ista!(oper,sig2,Lip,lambda, 
            xls,dat.od,x0,shrink,
            saxpy!,20*maxit1)


elseif redtype == "fista"

    
    println(">> LASSO-FISTA")                  
    sig2=1.0
    maxit1=oper.nc
    maxit2=oper.nc
    maxit3=10
    lambda=0.1
    beta=0.1
    x0 .= avg

    saxpy!(false,false,oper,x0,dfit)
    
    dfit .= dat.od .- dfit 
    rsdnorm=norm(dfit,2)^2
    x0norm1=norm(x0,1)
    lambda_scl=rsdnorm/x0norm1
    lambda=lambda*lambda_scl
    
    Lip=0.0
    tol=1.0e-06
    maxitn=1000
    Lip = normest(oper,saxpy!,tol=tol,maxit=maxitn)
    
    lasso_fista!(oper,sig2,Lip,lambda, 
            xls,dat.od,x0,shrink,
            saxpy!,20*maxit1)


elseif redtype == "lasso"

    println(">> LASSO-ADMM")              
    sig2=1.0
    maxit1=1000
    maxit2=20
    lambda=0.01
    beta=0.005
    x0 .= avg

    saxpy!(false,false,oper,x0,dfit)

    dfit .= dat.od .- dfit 
    rsdnorm=norm(dfit,2)^2
    x0norm1=norm(x0,1)
    lambda_scl=rsdnorm/x0norm1
    lambda=lambda*lambda_scl

    beta_scl=0.0
    tol=1.0e-06
    maxitn=1000
    beta_scl = normestfb(oper,saxpy!,tol=tol,maxit=maxitn)
    beta=beta*beta_scl

    lasso_admm!(grid,oper,lambda,sig2,beta,
                   xls,dat.od,x0,saxpy!,shrink,
                   maxit1,maxit2)

else

    println(">> RED-ADMM")              
    sig2=1.0
    maxit1=3000
    maxit2=24
    maxit3=24
    lambda=0.00001
    beta=0.000005
    x0 .= avg

    red_admm!(grid,oper,lambda,sig2,beta,
            xls,dat.od,x0,
            saxpy!,denoiser!,
            maxit1,maxit2,maxit3)

end#

println("norm(x) = ",norm(xls,2))

#fit=heatmap(reshape(xls,(n1,n2)),yflip=true,color=:jet1,clims=(0,2))
#display(plot(ext,fit,lsq))
#display(plot(fit,lsq,ext))

display(xls, n1, n2, o1, o2, d1, d2, join([redtype,"-model"]), true)

saxpy!(false,false,oper,xls,dfit)
for id=1:dat.nd
        rsd[id]=dat.od[id]-dfit[id]
end
println("data norm: ",norm(dat.od,2)," residual norm: ",norm(rsd,2))

x0 .= mdl .- xls
mfit = 100.0*( 1.0 - norm(x0,2)/norm(mdl,2))
println("model fit =", mfit)

display(x0, n1, n2, o1, o2, d1, d2, join([redtype,"-model-res"]), true)
println("mod res norm: ", norm(x0,2))

end # Main

end # module
