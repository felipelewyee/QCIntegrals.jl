module ERIs

using HypergeometricFunctions
using LinearAlgebra

greet() = print("Hello World!")

function doublefactorial(n)
    if (n <= 1)
        return 1
    end
    return n * doublefactorial(n - 2)
end

function F_n(n,x)
    return HypergeometricFunctions.M(n+1/2,n+3/2,-x)/(2*n+1)
end


struct shell
    l::Int8
    exp::Array{Float64}
    coef::Array{Float64}
    coord::Array{Float64}
    N::Array{Float64}
    orientaciones::Array{Array{Int}}
end

shell(l, exp, coef, coord) = shell(l, exp, coef, coord, ones(Int((l+1)*(l+2)/2)),[[lx,ly,l-lx-ly] for lx in l:-1:0 for ly in l-lx:-1:0])

function normalized_shell(l,exp,coef,coord)

    #normalizacion
    for i in eachindex(exp)
        coef[i] = coef[i]*(2*exp[i]/pi)^(3/4)*(4*exp[i])^(l/2)
    end
    
    norm = 0
    for j in eachindex(exp)
        for k in eachindex(exp)
            norm += coef[j]*coef[k]*(pi/(exp[j]+exp[k]))^(3/2)/(2*(exp[j]+exp[k]))^(l)
        end
    end

    for i in eachindex(exp)
        coef[i] = coef[i]/sqrt(norm)
    end
        
    orientaciones = []
    for lx in l:-1:0
        for ly in l-lx:-1:0
            lz = l-lx-ly
            append!(orientaciones,[[lx,ly,lz]])      
        end
    end        

    N = ones(Int((l+1)*(l+2)/2))
    for (i,(lx,ly,lz)) in enumerate(orientaciones)
        N[i] = sqrt(1/(doublefactorial(2*lx-1)*doublefactorial(2*ly-1)*doublefactorial(2*lz-1)))
    end
    
    return shell(l, exp, coef, coord, N, orientaciones)
    
end


end # module
