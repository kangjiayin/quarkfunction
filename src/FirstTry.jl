#First program of Julia
#2021.04.19




import QuadGK
using BenchmarkTools
using Plots


global A
global B
global z2
global z4
###自定义函数
stepfunction(x)= x>0
delta(x,y)= ==(x,y)
const mt=0.5;
const τ=ℯ^2-1;
const Λ=0.234;
const dd=0.93;
const Nf=4;
const rm=12/(33 - 2*Nf);
const ω=0.5;
z1 = 1.0;
z2 = 1;
z4 = 1;
const m = 0.003;
const step = 5;
const intstep = 2^7;
A=fill(2,intstep+1);
B=fill(2,intstep+1);
######################
const cutup = 1024;
const cutdown = 1/64;
X=[10^(log10(cutdown)+i/intstep*(-log10(cutdown) + log10(cutup))) for i= 0:intstep];
X1=[10^(log10(cutdown)+i/intstep*(-log10(cutdown) + log10(cutup))) for i= 0:intstep];
push!(X,361.);
F(x)=(1-exp(-x/(4*mt)^2))/x;
G(t)=8*pi^2*(dd*exp(-t/(ω^2))/ω^4+rm*F(t)/log(τ+(1+t/Λ^2)^2));
Z(i,j,k)=X[i]+X[j]-2*cos(k)*sqrt(X[i]*X[j]);
Aback(i,j,k)=sqrt(X[i]*X[j])cos(k)/X[i]+2(X[i]-sqrt(X[i]*X[j])cos(k))*(sqrt(X[i]*X[j])cos(k)-X[j])/(Z(i,j,k)*X[i]);
GetIntA_bef(i,j,θ)=sin(θ)^2*G(Z(i,j,θ))*Aback(i,j,θ);
GetIntB_bef(i,j,θ)=sin(θ)^2*G(Z(i,j,θ));



function getRealSumA(i,j)
    if j==1
        X[j]^2*A[j]/(X[j]*A[j]^2+B[j]^2)*IntA[i,j]
    else
        X[j]*(X[j]-X[j-1])*A[j]/(X[j]*A[j]^2+B[j]^2)*IntA[i,j]
    end
end
function realSumA(i)
    (2*pi)^(-3)*sum(j->getRealSumA(i,j),1:(intstep+1))
end
function getRealSumB(i,j)
    if j==1
        X[j]^2*B[j]/(X[j]*A[j]^2+B[j]^2)*IntB[i,j]
    else
        X[j]*(X[j]-X[j-1])*B[j]/(X[j]*A[j]^2+B[j]^2)*IntB[i,j]
    end
end
function realSumB(i)
    3*(2*pi)^(-3)*sum(j->getRealSumB(i,j),1:(intstep+1))
end
sumFA(t)=A[t]-z2-z2^2*realSumA(t);
sumFB(t)=B[t]-z4*m-z2^2*realSumB(t);
IntA=[QuadGK.quadgk(x->GetIntA_bef(i,j,x),0,pi)[1] for i=1:(intstep+2),j=1:(intstep+1)];
IntB=[QuadGK.quadgk(x->GetIntB_bef(i,j,x),0,pi)[1] for i=1:(intstep+2),j=1:(intstep+1)];
function jacobiAA(i,j)
    if j==1
        delta(i,j)-z2^2*(2*pi)^(-3)*X[j]^2*(B[j]^2-X[j]*A[j]^2)/(X[j]*A[j]^2+B[j]^2)^2*IntA[i,j]
    else
        delta(i,j)-z2^2*(2*pi)^(-3)*X[j]*(X[j]-X[j-1])*(B[j]^2-X[j]*A[j]^2)/(X[j]*A[j]^2+B[j]^2)^2*IntA[i,j]
    end
end
function jacobiAB(i,j)
    if j==1
        z2^2*(2*pi)^(-3)*X[j]^2*(B[j]*A[j]*2)/(X[j]*A[j]^2+B[j]^2)^2*IntA[i,j]
    else
        z2^2*(2*pi)^(-3)*X[j]*(X[j]-X[j-1])*(B[j]*A[j]*2)/(X[j]*A[j]^2+B[j]^2)^2*IntA[i,j]
    end
end
function jacobiBA(i,j)
    if j==1
        3*z2^2*(2*pi)^(-3)*X[j]^3*(B[j]*A[j]*2)/(X[j]*A[j]^2+B[j]^2)^2*IntB[i,j]
    else
        3*z2^2*(2*pi)^(-3)*X[j]^2*(X[j]-X[j-1])*(B[j]*A[j]*2)/(X[j]*A[j]^2+B[j]^2)^2*IntB[i,j]
    end
end
function jacobiBB(i,j)
    if j==1
        delta(i,j)+3*z2^2*(2*pi)^(-3)*X[j]^2*(B[j]^2-X[j]*A[j]^2)/(X[j]*A[j]^2+B[j]^2)^2*IntB[i,j]
    else
        delta(i,j)+3*z2^2*(2*pi)^(-3)*X[j]*(X[j]-X[j-1])*(B[j]^2-X[j]*A[j]^2)/(X[j]*A[j]^2+B[j]^2)^2*IntB[i,j]
    end
end
Δ=A=fill(1,intstep+1)
# while maximum(Δ)>0.01 
#     map(1:step) do st
#         global A
#         global B
#         global z2
#         global z4
#         global Δ
#         realSumA19=realSumA(intstep+2)
#         z2=(-1+sqrt(1+4realSumA19))/(2*realSumA19)
#         z4=1-z2^2*realSumB(intstep+2)/m
#         FA=[sumFA(i) for i=1:(intstep+1)]
#         FB=[sumFB(i) for i=1:(intstep+1)]
#         jacobiAA_M=[jacobiAA(i,j) for i=1:(intstep+1),j=1:(intstep+1)]
#         jacobiAB_M=[jacobiAB(i,j) for i=1:(intstep+1),j=1:(intstep+1)]
#         jacobiBA_M=[jacobiBA(i,j) for i=1:(intstep+1),j=1:(intstep+1)]
#         jacobiBB_M=[jacobiBB(i,j) for i=1:(intstep+1),j=1:(intstep+1)]
#         jacobi=[[jacobiAA_M jacobiAB_M] ;[jacobiBA_M jacobiBB_M]]
#         Δ=inv(jacobi)*[FA;FB]
#         A-=Δ[1:(intstep+1)]
#         B-=Δ[(intstep+2):(2intstep+2)]
#     end 
# end
realSumA19=realSumA(intstep+2)



#plot(X1,A,scale=:log10,title="A")
#plot(X1,B,scale=:log10,title="B")