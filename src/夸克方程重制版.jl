#Julia
#2021.08.10
using Gaussquad
using JLD2
using Plots
###参数列表
const mt = 0.5;
const τ = ℯ^2-1;
const Λ = 0.234;
const ω = 0.5;
const dd = (0.82)^3/ω;
const Nf = 4;
const rm = 12/(33 - 2*Nf);
const m = 0.003;

upper = 4
# step = 5;
intstep = 2^8;
cutup = 10. ^upper;
cutdown = 10. ^(-upper);

# 分配需要变量的内存
# z1 = 1.::Float64;
z2 = 1.::Float64;
z4 = 1.::Float64;
A = Array{Float64}(undef,intstep,1)
B = Array{Float64}(undef,intstep,1)
# 点与权重
k,w = gausslegendremesh(cutdown,cutup,intstep,2)
# w=w .*4/3
# 所需的函数
stepfunction(x)=x>0
delta(x,y)= ==(x,y)
F(x)=(1-exp(-x/(4*mt)^2))/x;
# 这里考虑下ω的几次方，这里取4
D(t)=8*pi^2*(dd*exp(-t/(ω^2))/ω^4+rm*F(t)/log(τ+(1+t/Λ^2)^2));

# 先做角度积分
# 采用第二类切比雪夫
# i代表外动量，j代表内动量
IntA=Array{Float64}(undef,intstep+1,intstep)
IntB=Array{Float64}(undef,intstep+1,intstep)
IntBz4=Array{Float64}(undef,intstep)
Threads.@threads for i= 1:intstep+1
    for j = 1:intstep
        q2=k[j]
        if i==intstep+1
            IntBz4[j]=(4/3)*gausschebyshevint64(z->D(19. ^2+q2+19. *sqrt(q2)*z))/(2*pi)^3
            k2=19. ^2
        else
            k2=k[i]
        end
        kdotq(z)=sqrt(k2*q2)*z
        IntA[i,j]=(4/3)*gausschebyshevint64(z->D(k2+q2-2*kdotq(z))*(3*k2*kdotq(z)+3*q2*kdotq(z)-2*k2*q2-4*(kdotq(z))^2)/k2/(k2+q2-2*kdotq(z)))/(2*pi)^3
        IntB[i,j]=(4/3)*gausschebyshevint64(z->(D(k2+q2-2*kdotq(z))-D(19. ^2+q2+19. *sqrt(q2)*z)))/(2*pi)^3
    end
end
#########角度积分完结####################

function getz4()
    sumB=0.::Float64
    for j=1:intstep
        sumB+=3*k[j]*B[j]/(k[j]*A[j]^2+B[j]^2)*IntBz4[j]*w[j]
    end
    return sumB
end 

function FAf(i)
    sum=0.::Float64
    for j=1:intstep
        sum+=k[j]*A[j]/(k[j]*A[j]^2+B[j]^2)*IntA[i,j]*w[j]
    end
    result=A[i]-z2-z2^2*sum
    return result
end

function FBf(i)
    sum=0.::Float64
    for j=1:intstep
        sum+=3*k[j]*B[j]/(k[j]*A[j]^2+B[j]^2)*IntB[i,j]*w[j]
    end
    result=B[i]-m-z2^2*sum
    return result
end

function renormalpoint()
    sumA=0.::Float64
    sumB=0.::Float64
    for j=1:intstep
        sumA+=k[j]*A[j]/(k[j]*A[j]^2+B[j]^2)*IntA[intstep+1,j]*w[j]
        sumB+=3*k[j]*B[j]/(k[j]*A[j]^2+B[j]^2)*IntB[intstep+1,j]*w[j]
    end
    return sumA,sumB
end

jacobifAA(i,j)=delta(i,j)-z2^2*w[j]*k[j]*(B[j]^2-k[j]*A[j]^2)/(k[j]*A[j]^2+B[j]^2)^2*IntA[i,j]
jacobifAB(i,j)=-z2^2*w[j]*k[j]*(-2*A[j]*B[j])/(k[j]*A[j]^2+B[j]^2)^2*IntA[i,j]
jacobifBA(i,j)=-3*z2^2*w[j]*k[j]*(-2*A[j]*B[j]*k[j])/(k[j]*A[j]^2+B[j]^2)^2*IntB[i,j]
jacobifBB(i,j)=delta(i,j)-3*z2^2*w[j]*k[j]*(k[j]*A[j]^2-B[j]^2)/(k[j]*A[j]^2+B[j]^2)^2*IntB[i,j]
# 给定a和b初值


A=fill(1.4,intstep);
B=fill(1.,intstep);
Δ=fill(1.,2*intstep)
st=0::Int
z2old=0.
# z4old=0.

# k19sumA, k19sumB=renormalpoint()
while maximum(abs.(Δ))>10^-6 || abs(1- z2/z2old) > 10^-6 
    global A, B, Δ, st, z2, z4, z2old#, z4old 
    st+=1
    k19sumA, k19sumB=renormalpoint()
    z2old=z2
    z2=(-1+sqrt(1+4*k19sumA))/(2*k19sumA)
    FA=[FAf(i) for i=1:(intstep)]
    FB=[FBf(i) for i=1:(intstep)]
    jacobiAA=[jacobifAA(i,j) for i=1:(intstep),j=1:(intstep)]
    jacobiAB=[jacobifAB(i,j) for i=1:(intstep),j=1:(intstep)]
    jacobiBA=[jacobifBA(i,j) for i=1:(intstep),j=1:(intstep)]
    jacobiBB=[jacobifBB(i,j) for i=1:(intstep),j=1:(intstep)]
    jacobi=[jacobiAA jacobiAB; jacobiBA jacobiBB]
    Δ=jacobi\[FA; FB]
    A-=Δ[1:(intstep)]
    B-=Δ[(intstep+1):(2intstep)]
end

# jldsave("/Users/kjy/Desktop/program/julia/Gamma5/data/ABk$upper.jld2";A, B, k)
# A4,B4,k4=load("./data/ABk4.jld2","A","B", "k")
# plot!(k6,A6,scale=:log10)
# scatter!(k6,A6,scale=:log10)
# plot(k,A,scale=:log10,title="A")
# plot(k,B,scale=:log10,title="B")

z4=(m-z2^2*getz4())/m