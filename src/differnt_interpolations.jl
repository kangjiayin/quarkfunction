# Differnt interpolations
# 2021/09/01
using JLD2
using Interpolations
using Dierckx
using Plots

# A, B,k=load("/Users/kjy/Desktop/program/julia/Gamma5/data/ABk1024.jld2","A","B","k")
# 约在688处为361
AD=Spline1D(k,A)
BD=Spline1D(k,B)

AI1=LinearInterpolation(k,A)
BI1=LinearInterpolation(k,B)

AI2=LinearInterpolation(k,A,extrapolation_bc=Line())
BI2=LinearInterpolation(k,B,extrapolation_bc=Line())