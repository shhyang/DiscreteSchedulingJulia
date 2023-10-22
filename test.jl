using Revise
includet("./DiscreteScheduling.jl")

using .DiscreteScheduling
using Graphs
using BenchmarkTools
using CDDLib
using Polyhedra

wg = linenet(5,2)
rrmsg = @time RateRegion(maxschedgraph(wg,1),4);
rrvf = @time RateRegion(genVF(maxschedgraph(wg,1)),4);
rrsg = @time RateRegion(schedgraph(wg,1,1),4);

p = Array{Any}(undef,4)
for i = 1:4
    p[i] = polyhedron(vrep(rrmsg[i]),CDDLib.Library())
    removevredundancy!(p[i])
end       

pp = polyhedron(vrep(rrvf),CDDLib.Library())
removevredundancy!(pp)

vf = genVF(maxschedgraph(wg,1));
sg = schedgraph(wg,1,1);
@btime simplecycles(sg.F);
@btime simplecycles(vf.F);
