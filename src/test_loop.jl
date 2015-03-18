##################################################################
##################################################################
loop = true
while loop
  tic()
ind +=1

println("up_y")
# update of y
yTMP = copy(yp)
y_est = yc - tau_y/rho_y
yp = prox_visi(y_est)

# Consensus

yc =  yp #(yp +  Fx +(-tau_y/rho_y + tau_yc  ))./2

println("up_x")
# update of x
#function estim_x(x)
xTMP = copy(x)
x,Fx = estime_x(x)
println("up s")
# # update of S
 S,sH = estime_S(S)

println("up_p")
#  # update of P
 Pti = x-tau_P/rho_P
 P = max(0,Pti)

println("up_t")
#  # update of T
 Tti = Hx - tau_T/rho_T
 T = prox_u(Tti,μ__tv)

println("up_v")
#  # update of V
 Vti = sH - tau_V/rho_V
 V = prox_u(Vti,μ__ss)

println("up_lagr")
 # update of Lagrange multipliers

 #tau_yc= tau_yc+rho_y *(y_est - yc)
 tau_y = tau_y + rho_y*(yc-Fx)
 tau_P = tau_P + rho_P*(P-x)
 tau_T = tau_T + rho_T*(T-Hx)
 tau_V = tau_V + rho_V*(V - sH)
 tau_S = tau_S + rho_S*(S-x)

n1 = norm(vec(x -xTMP))
n2 = norm(yc-yTMP)


  if (ind >= nbitermax)||( (n1<1E-8)&&(n2<1E-8) )
   loop = false
  end
  println(ind)
 toc()
  println(n1,"    ",n2)
end
