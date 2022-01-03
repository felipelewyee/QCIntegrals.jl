# This file download basis set from Basis Set Exchange (https://www.basissetexchange.org)
#
using Downloads

basis_list = [
	      "STO-2G"
	      "STO-3G"
	      "STO-4G"
	      "STO-5G"
	      "STO-6G"
	      "6-31G"
	      "6-31+G"
	      "6-31++G"
	      "6-31G*"
	      "6-31+G*"
	      "6-31++G*"
	      "6-31G**"
	      "6-31+G**"
	      "6-31++G**"
	      "6-311G"
	      "6-311+G"
	      "6-311++G"
	      "6-311G*"
	      "6-311+G*"
	      "6-311++G*"
	      "6-311G**"
	      "6-311+G**"
	      "6-311++G**"
	      "cc-pVDZ"
	      "cc-pVTZ"
	      "cc-pVQZ"
	      "cc-pV5Z"
	      "aug-cc-pVDZ"
	      "aug-cc-pVTZ"
	      "aug-cc-pVQZ"
	      "aug-cc-pV5Z"
	      "def2-SVP"
	      "def2-SVPD"
	      "def2-SVP"
	      "def2-TZVP"
	      "def2-TZVPD"
	      "def2-TZVPP"
	      "def2-TZVPPD"
	      "def2-QZVP"
	      "def2-QZVPD"
	      "def2-QZVPP"
	      "def2-QZVPPD"
	      "pc-0"
	      "pc-1"
	      "pc-2"
	      "pc-3"
	      "pc-4"
	      "pcseg-0"
	      "pcseg-1"
	      "pcseg-2"
	      "pcseg-3"
	      "pcseg-4"
	      "aug-pc-0"
	      "aug-pc-1"
	      "aug-pc-2"
	      "aug-pc-3"
	      "aug-pc-4"
	      "aug-pcseg-0"
	      "aug-pcseg-1"
	      "aug-pcseg-2"
	      "aug-pcseg-3"
	      "aug-pcseg-4"
	      ]

for basis in basis_list

    println(basis)
    url_basis_name = replace(lowercase(basis), "*" => "_st_")

    result = String(take!(Downloads.download("https://www.basissetexchange.org/download_basis/basis/"*url_basis_name*"/format/psi4/?uncontract_spdf=true",IOBuffer())))

    file = open(lowercase(basis*".gbs"),"w")
    write(file,replace(replace(result, "D+" => "E+"), "D-" => "E-"))
    close(file)
end
