# A GF library for Julia (only GF256 is implemented)

```
(@v1.5) pkg> add https://github.com/tavimori/gf
    Cloning git-repo `https://github.com/tavimori/gf`
   Updating git-repo `https://github.com/tavimori/gf`
   Updating registry at `~/.julia/registries/General`
######################################################################## 100.0%
julia> using gf

julia> a = rand(0:255,3,3)
3×3 Array{Int64,2}:
  78   33   48
 246  119  122
   8    6  116

julia> b = GF.(a)
3×3 Array{GF,2}:
 GF(78)   GF(33)   GF(48)
 GF(246)  GF(119)  GF(122)
 GF(8)    GF(6)    GF(116)

julia> r = rank(b)
3

```