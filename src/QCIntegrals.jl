module QCIntegrals

using HypergeometricFunctions
using SpecialFunctions
using LinearAlgebra

function doublefactorial(n)
    if (n <= 1)
        return 1
    end
    return n * doublefactorial(n - 2)
end

function F_n(n,x)
    if(x!=0)
        return gamma(0.5 + n) * gamma_inc(0.5 + n, x, 0)[1] / (2*x^(0.5 + n))
    else
        return HypergeometricFunctions.M(n+1/2,n+3/2,-x)/(2*n+1)
    end
end

struct shell
    l::Int8
    exp::Vector{Float64}
    coef::Vector{Float64}
    coord::Vector{Float64}
    N::Vector{Float64}
    orientaciones::Vector{Vector{Int}}
end

function build_shell(l,exp,coef,coord;normalized=true)

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
    if normalized
        for (i,(lx,ly,lz)) in enumerate(orientaciones)
            N[i] = sqrt(1/(doublefactorial(2*lx-1)*doublefactorial(2*ly-1)*doublefactorial(2*lz-1)))
        end
    end
    
    return shell(l, exp, coef, coord, N, orientaciones)
    
end

function build_zero_shell()

    #normalizacion
    l = 0
    exp = [0.0]
    coef = [1.0]
    coord = [0.0,0.0,0.0]
    N = [1.0]
    orientaciones = [[0,0,0]]
    
    return shell(l, exp, coef, coord, N, orientaciones)
    
end

g_a = build_shell(0,[13.01, 1.962, 0.4446, 0.122],[0.01968500,0.1379770,0.4781480,0.5012400],[0.0,0.0,0.0])
g_b = build_shell(1,[0.07896],[1.0],[0,0,0])
g_c = build_shell(0,[0.141],[1.0],[0.6068,-0.2383,-0.7169])
g_d = build_shell(0,[0.07896],[1.0],[0,0,0])

struct shell_pair
    
    g_a::shell
    g_b::shell
    
    exp::Vector{Float64}
    coef::Vector{Float64}
    coord::Vector{Vector{Float64}}
    alpha::Vector{Float64}
    
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

function get_set_of_shells(element,coord,basis_name;normalized=true,auxiliar=false)

    key_to_l = Dict("S" => 0, "P" => 1, "D" => 2, "F" => 3, "G" => 4, "H" => 5, "I" => 6)

    element = uppercase(element)
    shells = []
    
    # Basis file
    f = open("basis/"*basis_name*".gbs","r")
    lines = readlines(f)
    close(f)
    
    # Find_atom
    atom_line = 0
    for (i,line) in enumerate(lines)
        splitted_line = split(line)
        if(size(splitted_line)[1]>0)
            if(occursin("****",splitted_line[1]))
                next_line = split(lines[i+1])
                if(size(next_line)[1]>0)        
                    if(next_line[1] == element)
                        atom_line = i+1
                        break
                    end
                end
            end
        end
    end

    # Read shells
    num_of_primitives = 0
    l = 0
    exp = []
    coef = []
    for line in lines[atom_line+1:size(lines)[1]]
        splitted_line = split(line)
        if(size(splitted_line)[1]==3)
            exp = []
            coef = []
            try
                l = parse(Int,splitted_line[1])
            catch
                l = key_to_l[splitted_line[1]]
            end
            num_of_primitives = parse(Int,splitted_line[2])
        elseif(size(splitted_line)[1]==2)
            append!(exp,parse(Float64,splitted_line[1]))           
            append!(coef,parse(Float64,splitted_line[2]))
            if(size(exp)[1] == num_of_primitives)
                if(auxiliar)
                    shell = build_aux_shell(l,exp,coef,coord)
                else
                    shell = build_shell(l,exp,coef,coord,normalized=normalized)
                end
                push!(shells,shell)
            end
        end
        if(occursin("****",splitted_line[1]))
            break
        end
    end
            
    return shells
                        
end

function get_all_shells_from_xyz(molecule,basis_name;normalized=true,auxiliar=false)
    molecule = split(molecule,"\n")    
    natoms = size(molecule)[1]
    shells = []
    for iatom in 1:natoms
        atom = split(molecule[iatom])
        if size(atom)[1] >= 1
            shell = get_set_of_shells(atom[1],[parse(Float64,atom[2]),parse(Float64,atom[3]),parse(Float64,atom[4])],basis_name,normalized=normalized,auxiliar=auxiliar)
            append!(shells,shell)
        end
    end
    return shells        
end

function get_Z_xyz(molecule)
    symbol_to_Z = Dict("H" => 1, "He" => 2, "Li" => 3, "Be" => 4, "B" => 5, "C" => 6, "N" => 7, "O" => 8, "F" => 9, "Ne" => 10) 
    molecule = split(molecule,"\n")    
    natoms = size(molecule)[1]

    Zs = []
    coords = []
    for iatom in 1:natoms
        atom = split(molecule[iatom])
        if size(atom)[1] >= 1
            append!(Zs,symbol_to_Z[atom[1]])
            append!(coords,[[parse(Float64,atom[2]),parse(Float64,atom[3]),parse(Float64,atom[4])]])
        end
    end

    return Zs,coords
    
end

function get_nbf(shells)
    
    nbf = 0
    
    for shell in shells
        nbf += size(shell.orientaciones)[1]
    end
    
    return nbf
    
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

function overlap(pairAB)
    
    g_a = pairAB.g_a
    g_b = pairAB.g_b

    Eab = build_E(g_a.l,g_b.l,pairAB)
    result = zeros(size(g_a.orientaciones)[1],size(g_b.orientaciones)[1])

    for (a,(lax,lay,laz)) in enumerate(g_a.orientaciones)
        for (b,(lbx,lby,lbz)) in enumerate(g_b.orientaciones)
        
            integral = 0.0
            for ii in eachindex(pairAB.exp)
                t_z = 0
                t_y = 0
                t_x = 0
                integral += Eab[ii][1][lax+1][lbx+1][t_x+1]*
                            Eab[ii][2][lay+1][lby+1][t_y+1]*
                            Eab[ii][3][laz+1][lbz+1][t_z+1]*
                            (pi/pairAB.exp[ii])^1.5*pairAB.coef[ii]
            end           
            result[a,b] = g_a.N[a]*g_b.N[b]*integral
        end
    end
    
    return result
    
end

function get_S(shells)
    
    nbf = get_nbf(shells)
    
    S = zeros(nbf,nbf)
    i = 1
    for g_a in shells
        j = 1
        for g_b in shells
            pair = build_shell_pair(g_a,g_b)
            results = overlap(pair)
            S[i:i+size(g_a.orientaciones)[1]-1,j:j+size(g_b.orientaciones)[1]-1] = results
            j += size(g_b.orientaciones)[1]
        end
        i += size(g_a.orientaciones)[1]
    end
    
    return S
    
end

function Tij(la,lb,xyz,pairAB)
    
    Eab = build_E(pairAB.g_a.l,pairAB.g_b.l,pairAB)
    
    S = [[[0.0 for j in 0:lb] for i in 0:la] for ii in eachindex(pairAB.exp)]
    T = [[[0.0 for j in 0:lb] for i in 0:la] for ii in eachindex(pairAB.exp)]
    
    for ii in eachindex(pairAB.exp)
        i = 0
        t = 0
        for j in 0:lb
            S[ii][i+1][j+1] = Eab[ii][xyz][i+1][j+1][t+1]*(pi/pairAB.exp[ii])^0.5            
        end
        for i in 1:la
            for j in max(0,lb-la+i):lb
                S[ii][i+1][j+1] = Eab[ii][xyz][i+1][j+1][t+1]*(pi/pairAB.exp[ii])^0.5
            end
        end
    end

    ii = 0
    for jprim in 1:size(pairAB.g_b.exp)[1], iprim in 1:size(pairAB.g_a.exp)[1]
        ii += 1
        i=0
        j=0
        T[ii][i+1][j+1] = (pairAB.g_a.exp[iprim] - 2*pairAB.g_a.exp[iprim]^2 * ((pairAB.coord[ii][xyz]-pairAB.g_a.coord[xyz])^2 + 1/(2*pairAB.exp[ii]))) * S[ii][i+1][j+1]
        for j in 1:lb
            T[ii][i+1][j+1] = (pairAB.coord[ii][xyz]-pairAB.g_b.coord[xyz])*T[ii][i+1][j+1] + pairAB.g_a.exp[iprim]/pairAB.exp[ii]*2*pairAB.g_b.exp[jprim]*S[ii][i+1][j+1] 
            if(j>1)
                T[ii][i+1][j+1] += (j-1)/(2*pairAB.exp[ii])*T[ii][i+1][j-1] - pairAB.g_a.exp[iprim]/pairAB.exp[ii]*(j)*S[ii][i+1][j-1]
            end
        end
        for i in 1:la
            for j in max(0,lb-la+i):lb
                T[ii][i+1][j+1] = (pairAB.coord[ii][xyz]-pairAB.g_a.coord[xyz])*T[ii][i][j+1] + pairAB.g_b.exp[jprim]/pairAB.exp[ii]*2*pairAB.g_a.exp[iprim]*S[ii][i+1][j+1] 
                if(i>1)
                    T[ii][i+1][j+1] += (i-1)/(2*pairAB.exp[ii])*T[ii][i-1][j+1] - pairAB.g_b.exp[jprim]/pairAB.exp[ii]*(i-1)*S[ii][i-1][j+1]
                end
                if(j>0)
                    T[ii][i+1][j+1] += (j)/(2*pairAB.exp[ii])*T[ii][i][j]# - pairAB.g_b.exp[jprim]/pairAB.exp[ii]*(i-1)*S[ii][i-1][j+1]
                end
            end
        end
    end       
        
    return [T[ii][la+1][lb+1] for ii in eachindex(pairAB.exp)]
        
end

function kinetic(pairAB)
    
    g_a = pairAB.g_a
    g_b = pairAB.g_b

    Eab = build_E(g_a.l,g_b.l,pairAB)

    result = zeros((size(g_a.orientaciones)[1],size(g_b.orientaciones)[1]))

    for (a,(lax,lay,laz)) in enumerate(g_a.orientaciones)
        for (b,(lbx,lby,lbz)) in enumerate(g_b.orientaciones)
        
            Tx = Tij(lax,lbx,1,pairAB)
            Ty = Tij(lay,lby,2,pairAB)
            Tz = Tij(laz,lbz,3,pairAB)
        
            integral = 0.0
            for ii in eachindex(pairAB.exp)
                tx,ty,tz = 0,0,0
                Sx = Eab[ii][1][lax+1][lbx+1][tx+1]*(pi/pairAB.exp[ii])^0.5
                Sy = Eab[ii][2][lay+1][lby+1][ty+1]*(pi/pairAB.exp[ii])^0.5
                Sz = Eab[ii][3][laz+1][lbz+1][tz+1]*(pi/pairAB.exp[ii])^0.5
                integral += pairAB.coef[ii]*(Tx[ii]*Sy*Sz + Ty[ii]*Sx*Sz + Tz[ii]*Sx*Sy)
            end           
            result[a,b] = g_a.N[a]*g_b.N[b]*integral
        end
    end
    
    return result
    
end

function get_T(shells)
    
    nbf = get_nbf(shells)
    
    T = zeros(nbf,nbf)
    i = 1
    for g_a in shells
        j = 1
        for g_b in shells
            pair = build_shell_pair(g_a,g_b)
            results = kinetic(pair)
            T[i:i+size(g_a.orientaciones)[1]-1,j:j+size(g_b.orientaciones)[1]-1] = results
            j += size(g_b.orientaciones)[1]
        end
        i += size(g_a.orientaciones)[1]
    end
    
    return T
    
end

function build_R_one(coord,pairAB)
    
    p = pairAB.exp
    Rpc = [pairAB.coord[ii]-coord for ii in eachindex(pairAB.exp)]
    F = [[F_n(n,p[ii]*norm(Rpc[ii])^2) for n in 0:pairAB.g_a.l + pairAB.g_b.l] for ii in eachindex(pairAB.exp)]
    lt = pairAB.g_a.l + pairAB.g_b.l
    R = [[[[[0.0 for lz in 0:lt-n-lx-ly] for ly in 0:lt-n-lx] for lx in 0:lt-n] for n in 0:lt] for ii in eachindex(pairAB.exp)]

    for ii in eachindex(pairAB.exp)
        for nini in 0:lt
            R[ii][nini+1][1][1][1] = (-2*p[ii])^nini*2*pi/p[ii]*F[ii][nini+1]
            l = 0
            for n in nini-1:-1:0
                l += 1
                for lx in 0:l
                    for ly in 0:l-lx
                        lz = l-lx-ly
                        if(lz>=2)
                            R[ii][n+1][lx+1][ly+1][lz+1] += (lz-1)*R[ii][n+2][lx+1][ly+1][lz-1]
                        end
                        if(lz>=1)
                            R[ii][n+1][lx+1][ly+1][lz+1] += Rpc[ii][3]*R[ii][n+2][lx+1][ly+1][lz]
                        end
                    end
                    ly = l-lx
                    if(ly>=2)
                        R[ii][n+1][lx+1][l-lx+1][1] += (ly-1)*R[ii][n+2][lx+1][l-lx-1][1]
                    end
                    if(ly>=1)
                        R[ii][n+1][lx+1][l-lx+1][1] += Rpc[ii][2]*R[ii][n+2][lx+1][l-lx][1]
                    end
                end
                lx = l
                if(lx>=2)
                    R[ii][n+1][l+1][1][1] += (lx-1)*R[ii][n+2][l-1][1][1]
                end
                if(lx>=1)
                    R[ii][n+1][l+1][1][1] += Rpc[ii][1]*R[ii][n+2][l][1][1]                
                end
            end
        end
    end
             
    return R

end

function potential(pairAB,Z,coord)

    g_a = pairAB.g_a
    g_b = pairAB.g_b
    
    R = build_R_one(coord,pairAB)
    Eab = build_E(pairAB.g_a.l,pairAB.g_b.l,pairAB)
    
    result = zeros(size(g_a.orientaciones)[1],size(g_b.orientaciones)[1])

    for (a,(lax,lay,laz)) in enumerate(g_a.orientaciones)
        for (b,(lbx,lby,lbz)) in enumerate(g_b.orientaciones)
        
            integral = 0.0
            for ii in eachindex(pairAB.exp)
                aux3 = 0
                for t_z in 0:laz+lbz
                    aux2 = 0    
                    for t_y in 0:lay+lby
                        aux1 = 0       
                        for t_x in 0:lax+lbx
                            aux1 += Eab[ii][1][lax+1][lbx+1][t_x+1]*R[ii][1][t_x+1][t_y+1][t_z+1]
                        end
                        aux2 += Eab[ii][2][lay+1][lby+1][t_y+1]*aux1
                    end
                    aux3 += Eab[ii][3][laz+1][lbz+1][t_z+1]*aux2
                end
                integral += aux3*pairAB.coef[ii]
            end           
            result[a,b] = g_a.N[a]*g_b.N[b]*integral
        end
    end
    
    return Z*result
    
end

function get_V(shells,Zs,coords)
    
    nbf = get_nbf(shells)
    
    V = zeros(nbf,nbf)
    i = 1
    for g_a in shells
        j = 1
        for g_b in shells
            pair = build_shell_pair(g_a,g_b)
            for (Z,coord) in zip(Zs,coords)
                results = potential(pair,Z,coord)
                V[i:i+size(g_a.orientaciones)[1]-1,j:j+size(g_b.orientaciones)[1]-1] += results
            end
            j += size(g_b.orientaciones)[1]
        end
        i += size(g_a.orientaciones)[1]
    end
    
    return V
    
end

function ERI(pairAB,pairCD)

    g_a = pairAB.g_a
    g_b = pairAB.g_b
    g_c = pairCD.g_a
    g_d = pairCD.g_b
    
    lt = pairAB.g_a.l + pairAB.g_b.l + pairCD.g_a.l + pairCD.g_b.l

    # Build boys
    alpha = [[pairAB.exp[ii]*pairCD.exp[jj]/(pairAB.exp[ii]+pairCD.exp[jj]) for jj in eachindex(pairCD.exp)] for ii in eachindex(pairAB.exp)]
    Rpq = [[pairAB.coord[ii]-pairCD.coord[jj] for jj in eachindex(pairCD.exp)] for ii in eachindex(pairAB.exp)]
    F = [[[F_n(n,alpha[ii][jj]*norm(Rpq[ii][jj])^2) for n in 0:lt] for jj in eachindex(pairCD.exp)] for ii in eachindex(pairAB.exp)]

    # Build R
    R = build_R(lt,alpha,Rpq,F,pairAB,pairCD)

    # Build E
    Eab = build_E(pairAB.g_a.l,pairAB.g_b.l,pairAB)
    Ecd = build_E(pairCD.g_a.l,pairCD.g_b.l,pairCD)

    result = zeros((size(g_a.orientaciones)[1],size(g_b.orientaciones)[1],size(g_c.orientaciones)[1],size(g_d.orientaciones)[1]))
    for (a,(lax,lay,laz)) in enumerate(g_a.orientaciones)
        for (b,(lbx,lby,lbz)) in enumerate(g_b.orientaciones)
            for (c,(lcx,lcy,lcz)) in enumerate(g_c.orientaciones)
                for (d,(ldx,ldy,ldz)) in enumerate(g_d.orientaciones)

                    integral = 0.0
                    for ii in eachindex(pairAB.exp)
                        aux7 = 0
                        for t_z in 0:laz+lbz
                            aux6 = 0    
                            for t_y in 0:lay+lby
                                aux5 = 0        
                                for t_x in 0:lax+lbx
                                    aux4 = 0
                                    for jj in eachindex(pairCD.exp)
                                        aux3 = 0
                                        for u_z in 0:lcz+ldz
                                            aux2 = 0
                                            for u_y in 0:lcy+ldy
                                                aux1 = 0
                                                for u_x in 0:lcx+ldx
                                                    aux1 += (-1)^(u_x+u_y+u_z)*Ecd[jj][1][lcx+1][ldx+1][u_x+1]*R[ii][jj][1][t_x+u_x+1][t_y+u_y+1][t_z+u_z+1]
                                                end
                                                aux2 += Ecd[jj][2][lcy+1][ldy+1][u_y+1]*aux1
                                            end
                                            aux3 += Ecd[jj][3][lcz+1][ldz+1][u_z+1]*aux2
                                        end
                                        aux4 += aux3*pairCD.coef[jj]
                                    end
                                    aux5 += Eab[ii][1][lax+1][lbx+1][t_x+1]*aux4
                                end
                                aux6 += Eab[ii][2][lay+1][lby+1][t_y+1]*aux5
                            end
                            aux7 += Eab[ii][3][laz+1][lbz+1][t_z+1]*aux6
                        end
                        integral += aux7*pairAB.coef[ii]
                    end           
                    result[a,b,c,d] = g_a.N[a]*g_b.N[b]*g_c.N[c]*g_d.N[d]*integral
                end
            end
        end
    end
    
    return result
    
end      

function get_I4(shells)
    
    nbf = get_nbf(shells)
    
    I = zeros(nbf,nbf,nbf,nbf)
    i = 1
    for g_a in shells
        j = 1
        for g_b in shells
            pairAB = build_shell_pair(g_a,g_b)
            k = 1
            for g_c in shells
                l = 1
                for g_d in shells
                    pairCD = build_shell_pair(g_c,g_d)
                    results = ERI(pairAB,pairCD)
                    I[i:i+size(g_a.orientaciones)[1]-1,
                      j:j+size(g_b.orientaciones)[1]-1,
                      k:k+size(g_c.orientaciones)[1]-1,
                      l:l+size(g_d.orientaciones)[1]-1] = results
                    l += size(g_d.orientaciones)[1]
                end
                k += size(g_c.orientaciones)[1]
            end
            j += size(g_b.orientaciones)[1]
        end
        i += size(g_a.orientaciones)[1]
    end
    
    return I
    
end

function normalize_aux_shell(g_a)
    
    g_1 = build_zero_shell()
    
    pairAB = build_shell_pair(g_a,g_1)
    results = ERI(pairAB,pairAB)

    N = diag(results[1:size(g_a.orientaciones)[1],1,1:size(g_a.orientaciones)[1],1])
    N = 1.0 ./ (sqrt.(N))
    
    return N
    
    l = g_a.l
    orientaciones = []
    for lx in l:-1:0
        for ly in l-lx:-1:0
            lz = l-lx-ly
            append!(orientaciones,[[lx,ly,lz]])      
        end
    end        
    
    aux_normalized_shell = shell(g_a.l,g_a.exp,g_a.coef,g_a.coord,N,orientaciones)
    
    return aux_normalized_shell
    
end

function build_aux_shell(l,exp,coef,coord)

    orientaciones = []
    for lx in l:-1:0
        for ly in l-lx:-1:0
            lz = l-lx-ly
            append!(orientaciones,[[lx,ly,lz]])      
        end
    end        

    N = ones(Int((l+1)*(l+2)/2))

    tmp_shell = shell(l, exp, coef, coord, N, orientaciones)

    N = normalize_aux_shell(tmp_shell)    
    
    return shell(l, exp, coef, coord, N, orientaciones)
    
end

function get_I2(shells_aux)
    
    g_1 = build_zero_shell()
    
    nbfaux = get_nbf(shells_aux)
    
    I = zeros(nbfaux,nbfaux)
    i = 1
    for g_a in shells_aux
        pairAB = build_shell_pair(g_a,g_1)
        j = 1
        for g_b in shells_aux
            pairCD = build_shell_pair(g_b,g_1)
            results = ERI(pairAB,pairCD)
            I[i:i+size(g_a.orientaciones)[1]-1,
              j:j+size(g_b.orientaciones)[1]-1] = results[1:size(g_a.orientaciones)[1],1,1:size(g_b.orientaciones)[1],1]
            j += size(g_b.orientaciones)[1]
        end
        i += size(g_a.orientaciones)[1]
    end
    
    return I
    
end

function get_I3(shells,shells_aux)
    
    g_1 = build_zero_shell()
    
    nbf = get_nbf(shells)
    nbfaux = get_nbf(shells_aux)
    
    I = zeros(nbf,nbf,nbfaux)
    i = 1
    for g_a in shells
        j = 1
        for g_b in shells
            pairAB = build_shell_pair(g_a,g_b)
            k = 1
            for g_c in shells_aux
                pairCD = build_shell_pair(g_c,g_1)
                results = ERI(pairAB,pairCD)
                I[i:i+size(g_a.orientaciones)[1]-1,
                  j:j+size(g_b.orientaciones)[1]-1,
                  k:k+size(g_c.orientaciones)[1]-1] = results[1:size(g_a.orientaciones)[1],1:size(g_b.orientaciones)[1],1:size(g_c.orientaciones)[1],1]
                k += size(g_c.orientaciones)[1]
            end
            j += size(g_b.orientaciones)[1]
        end
        i += size(g_a.orientaciones)[1]
    end
    
    return I
    
end

@btime I = get_I3(shells,shells_aux)

end # module
