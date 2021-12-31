module QCInts

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

struct shell_pair
    
    g_a::shell
    g_b::shell
    
    exp::Array{Float64}
    coef::Array{Float64}
    coord::Array{Array{Float64}}
    alpha::Array{Float64}
    
end
    
function build_shell_pair(g_a, g_b)
    exp = vec((g_a.exp .+ g_b.exp')')
    coef = vec((g_a.coef .* g_b.coef')')

    coord = []
    for i in eachindex(g_a.exp)
        for j in eachindex(g_b.exp)
            append!(coord,[(g_a.exp[i]*g_a.coord+g_b.exp[j]*g_b.coord)/(g_a.exp[i]+g_b.exp[j])])
        end
    end
            
    alpha = vec(((g_a.exp .* g_b.exp')./(g_a.exp .+ g_b.exp'))')
    
    return shell_pair(g_a, g_b, exp, coef, coord, alpha)

end

function get_set_of_shells(element,coord,basis_name)

    element = uppercase(element)
    shells = []
    
    # Basis file
    f = open("basis/"*basis_name*".bas","r")
    lines = readlines(f)
    close(f)
    
    # Find_atom
    atom_line = 0
    for (i,line) in enumerate(lines)
        splitted_line = split(line)
        if(size(splitted_line)[1]>0)
            if(occursin("O-",splitted_line[1]))
                if(splitted_line[2] == element)
                    atom_line = i
                    break
                end
            end
        end
    end

    # Pass the scheme contraction info line
    num_of_shells = parse(Int,lines[atom_line + 2])
    
    # Read shells
    new_shell = true
    idx_shell = 0
    idx_prim = 0
    exp = []
    coef = []
    idx,l,num_of_primitives = 0,0,0
    for line in lines[atom_line + 3:size(lines)[1]]
        splitted_line = split(line)
        #If new shell, read angular momentum and number of primitives
        if new_shell
            idx_shell += 1
            idx_prim = 0        
            idx,l,num_of_primitives = parse(Int,splitted_line[1]),parse(Int,splitted_line[2]),parse(Int,splitted_line[3])
            exp = []
            coef = []
            new_shell = false
            continue
        end
        #If not a new shell, read exponents and coefficients
        if !new_shell
            idx_prim += 1
            append!(exp,parse(Float64,splitted_line[1]))           
            append!(coef,parse(Float64,splitted_line[2]))
            if(idx_prim == num_of_primitives)
                push!(shells,normalized_shell(l,exp,coef,coord))
                #If all shells readed, break, else, read new shell
                if(idx_shell == num_of_shells)
                    break
                else
                    new_shell = true
                end
            end
            continue
        end
    end
            
    return shells
                        
end

function get_all_shells_from_xyz(molecule,basis_name)
    molecule = split(molecule,"\n")
    natoms = size(molecule)[1]
    shells = []
    for iatom in 1:natoms
        atom = split(molecule[iatom])
        if size(atom)[1] >= 1
            shell = get_set_of_shells(atom[1],[parse(Float64,atom[2]),parse(Float64,atom[3]),parse(Float64,atom[4])],basis_name)
            push!(shells,shell)
        end
    end
    return shells
end

function build_R(lt,alpha,Rpq,F,pairAB,pairCD)
    R = [[[[[[0.0 for lz in 0:lt-n-lx-ly] for ly in 0:lt-n-lx] for lx in 0:lt-n] for n in 0:lt] for jj in eachindex(pairCD.exp)] for ii in eachindex(pairAB.exp)]
    for ii in eachindex(pairAB.exp)
        for jj in eachindex(pairCD.exp)
            for nini in 0:lt
                R[ii][jj][nini+1][1][1][1] = (-2*alpha[ii][jj])^nini*2*pi^(5/2)/(pairAB.exp[ii]*pairCD.exp[jj]*sqrt(pairAB.exp[ii]+pairCD.exp[jj]))*F[ii][jj][nini+1]
                l = 0
                for n in nini-1:-1:0
                    l += 1
                    for lx in 0:l
                        for ly in 0:l-lx
                            lz = l-lx-ly
                            if(lz>=2)
                                R[ii][jj][n+1][lx+1][ly+1][lz+1] += (lz-1)*R[ii][jj][n+2][lx+1][ly+1][lz-1]
                            end
                            if(lz>=1)
                                R[ii][jj][n+1][lx+1][ly+1][lz+1] += Rpq[ii][jj][3]*R[ii][jj][n+2][lx+1][ly+1][lz]
                            end
                        end
                        ly = l-lx
                        if(ly>=2)
                            R[ii][jj][n+1][lx+1][l-lx+1][1] += (ly-1)*R[ii][jj][n+2][lx+1][l-lx-1][1]
                        end
                        if(ly>=1)
                            R[ii][jj][n+1][lx+1][l-lx+1][1] += Rpq[ii][jj][2]*R[ii][jj][n+2][lx+1][l-lx][1]
                        end
                    end
                    lx = l
                    if(lx>=2)
                        R[ii][jj][n+1][l+1][1][1] += (lx-1)*R[ii][jj][n+2][l-1][1][1]
                    end
                    if(lx>=1)
                        R[ii][jj][n+1][l+1][1][1] += Rpq[ii][jj][1]*R[ii][jj][n+2][l][1][1]                
                    end
                end
            end
        end
    end
             
    return R

end

function build_E(l1,l2,pair)
    E = [[[[[0.0 for t in 0:la+lb] for lb in 0:l2] for la in 0:l1] for xyz in 1:3] for ii in eachindex(pair.exp)]

    for ii in eachindex(pair.exp)
        for xyz in 1:3
            E[ii][xyz][1][1][1] = exp(-pair.alpha[ii]*(pair.g_b.coord[xyz]-pair.g_a.coord[xyz])^2)
            for la in 0:l1
                for lb in 0:l2
                    for t in 0:la+lb
                        if(lb>0)
                            if(t>0)
                                E[ii][xyz][la+1][lb+1][t+1] += 1/(2*pair.exp[ii])*E[ii][xyz][la+1][lb][t]
                            end
                            if(t<=la+lb-1)
                                E[ii][xyz][la+1][lb+1][t+1] += (pair.coord[ii][xyz]-pair.g_b.coord[xyz])*E[ii][xyz][la+1][lb][t+1]
                            end
                            if(t+1<=la+lb-1)
                                E[ii][xyz][la+1][lb+1][t+1] += (t+1)*E[ii][xyz][la+1][lb][t+2]
                            end
                        elseif(la>0)
                            if(t>0)
                                E[ii][xyz][la+1][lb+1][t+1] += 1/(2*pair.exp[ii])*E[ii][xyz][la][lb+1][t]
                            end
                            if(t<=la-1+lb)
                                E[ii][xyz][la+1][lb+1][t+1] += (pair.coord[ii][xyz]-pair.g_a.coord[xyz])*E[ii][xyz][la][lb+1][t+1]
                            end
                            if(t+1<=la-1+lb)
                                E[ii][xyz][la+1][lb+1][t+1] += (t+1)*E[ii][xyz][la][lb+1][t+2]
                            end
                        end
                    end
                end
            end
        end
    end

    return E

end




end # module
