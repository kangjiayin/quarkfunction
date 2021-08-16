#Julia
#2021.08.10
include("/Users/kjy/Desktop/program/julia/module/GaussQuad/src/gaussQuad.jl")
include("/Users/kjy/Desktop/program/julia/module/GaussQuad/src/gaussmesh.jl")
###参数列表
const mt=0.5;
const τ=ℯ^2-1;
const Λ=0.234;
const ω=0.5;
const dd=(0.82)^3/ω;
const Nf=4;
const rm=12/(33 - 2*Nf);
const m = 0.003;


step = 5;
intstep = 2^5;
cutup = 10^4+0.;
cutdown = 10^(-4);

# 分配需要变量的内存
z1 = 1.::Float64;
z2 = 1.::Float64;
z4 = 1.::Float64;
A=Array{Float64}(undef,intstep,1)
B=Array{Float64}(undef,intstep,1)
# 点与权重
k,w= gausslegendremesh(cutdown,cutup,intstep,2)

stepfunction(x)= x>0
delta(x,y)= ==(x,y)
F(x)=(1-exp(-x/(4*mt)^2))/x;
D(t)=8*pi^2*(dd*exp(-t/(ω^2))/ω^4+rm*F(t)/log(τ+(1+t/Λ^2)^2));


# 先做角度积分
# 采用第二类切比雪夫
# i代表外动量，j代表内动量
IntA=Array{Float64}(undef,intstep+1,intstep)
IntB=Array{Float64}(undef,intstep+1,intstep)
for i= 1:intstep+1, j= 1:intstep
    if i==intstep+1
        k2=19. ^2
    else
        k2=k[i]
    end
    q2=k[j]
        kdotq(z)=sqrt(k2*q2)*z
        IntA[i,j]=gausschebyshevint64(z->D(k2+q2-2kdotq(z))*(3*k2*kdotq(z)+3*q2*kdotq(z)-2*k2*q2-4*(kdotq(z))^2)/k2/(k2+q2-2kdotq(z)))
        IntB[i,j]=gausschebyshevint64(z->D(k2+q2-sqrt(k2*q2)z))
end

function FAf(i)
    sum=0.::Float64
    for j=1:intstep
        sum+=k[j]*A[j]/(k[j]*A[j]^2+B[j]^2)*IntA[i,j]*w[j]
    end
    sum*=1/(2*pi)^3
    result=A[i]-z2-z2^2*sum
    return result
end

function FBf(i)
    sum=0.::Float64
    for j=1:intstep
        sum+=3*k[j]*B[j]/(k[j]*A[j]^2+B[j]^2)*IntB[i,j]*w[j]
    end
    sum*=1/(2*pi)^3
    result=B[i]-m*z4-z2^2*sum
    return result
end

jacobifAA(i,j)=delta(i,j)-z2^2*(2*pi)^(-3)*w[j]*k[j]*(B[j]^2-k[j]*A[j]^2)/(k[j]*A[j]^2+B[j]^2)^2*IntA[i,j]
jacobifAB(i,j)=-z2^2*(2*pi)^(-3)*w[j]*k[j]*(-2*A[j]*B[j])/(k[j]*A[j]^2+B[j]^2)^2*IntA[i,j]
jacobifBA(i,j)=-3*z2^2*(2*pi)^(-3)*w[j]*k[j]*(-2*A[j]*B[j]*k[j])/(k[j]*A[j]^2+B[j]^2)^2*IntB[i,j]
jacobifBB(i,j)=delta(i,j)-3*z2^2*(2*pi)^(-3)*w[j]*k[j]*(k[j]*A[j]^2-B[j]^2)/(k[j]*A[j]^2+B[j]^2)^2*IntB[i,j]


A=fill(2,intstep);
B=fill(1,intstep);
Δ=Array{Float64}(undef,2*intstep,1)


        # global A
        # global B
        # global z2
        # global z4
        # global Δ
        # realSumA19=realSumA(intstep+2)
        # z2=(-1+sqrt(1+4realSumA19))/(2*realSumA19)
        # z4=1-z2^2*realSumB(intstep+2)/m
        FA=[FAf(i) for i=1:(intstep)]
        FB=[FBf(i) for i=1:(intstep)]
        jacobiAA_M=[jacobifAA(i,j) for i=1:(intstep),j=1:(intstep)]
        jacobiAB_M=[jacobifAB(i,j) for i=1:(intstep),j=1:(intstep)]
        jacobiBA_M=[jacobifBA(i,j) for i=1:(intstep),j=1:(intstep)]
        jacobiBB_M=[jacobifBB(i,j) for i=1:(intstep),j=1:(intstep)]
        jacobi=[jacobiAA_M jacobiAB_M ;jacobiBA_M jacobiBB_M]
        Δ=jacobi\[FA; FB]
        A-=Δ[1:(intstep)]
        B-=Δ[(intstep+1):(2intstep)]




#plot(X1,A,scale=:log10,title="A")
#plot(X1,B,scale=:log10,title="B")