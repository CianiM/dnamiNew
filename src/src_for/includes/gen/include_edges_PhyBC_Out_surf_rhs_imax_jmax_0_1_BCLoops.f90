

!***********************************************************
!                                                           
! Start building layers for BC : imax jmax None ************
!                                                           
!***********************************************************




!***********************************************************
!                                                           
! BC layer: 0 1 None ***************************************
!                                                           
!***********************************************************




!***********************************************************
!                                                           
! building source terms in RHS for layer 0 1 None d(rho)/dt 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! deltaxI*([rho*u]_1x)+deltayI*([rho*v]_1y)+(1.0_wp/c**2*((u*(-c**2*[rho]_1x+1.0_wp*[p]_1x))*deltaxI+0.5_wp*(((c+u)*(c*rho*[u]_1x+1.0_wp*[p]_1x))*deltaxI+-((c+u)*(c*rho*[u]_1x+1.0_wp*[p]_1x))*deltaxI)))
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_rhs_rho_dx_0_nxp2p0p0nyp2m1k = q(nx+2+0+0,ny+2-1,indvars(1))*q(nx+2+0+0,ny+2-1,indvars(2))

d1_rhs_rho_dx_0_nxp2p0m1nyp2m1k = q(nx+2+0-1,ny+2-1,indvars(1))*q(nx+2+0-1,ny+2-1,indvars(2))

d1_rhs_rho_dx_0_nxp2p0m2nyp2m1k = q(nx+2+0-2,ny+2-1,indvars(1))*q(nx+2+0-2,ny+2-1,indvars(2))

d1_rhs_rho_dx_0_nxp2p0nyp2m1k = 1.5_wp*d1_rhs_rho_dx_0_nxp2p0p0nyp2m1k-&
          2.0_wp*d1_rhs_rho_dx_0_nxp2p0m1nyp2m1k+&
          0.5_wp*d1_rhs_rho_dx_0_nxp2p0m2nyp2m1k

d1_rhs_rho_dx_0_nxp2p0nyp2m1k = d1_rhs_rho_dx_0_nxp2p0nyp2m1k*param_float(1)

d1_rhs_rho_dx_1_nxp2p0p0nyp2m1k = q(nx+2+0+0,ny+2-1,indvars(1))

d1_rhs_rho_dx_1_nxp2p0m1nyp2m1k = q(nx+2+0-1,ny+2-1,indvars(1))

d1_rhs_rho_dx_1_nxp2p0m2nyp2m1k = q(nx+2+0-2,ny+2-1,indvars(1))

d1_rhs_rho_dx_1_nxp2p0nyp2m1k = 1.5_wp*d1_rhs_rho_dx_1_nxp2p0p0nyp2m1k-&
          2.0_wp*d1_rhs_rho_dx_1_nxp2p0m1nyp2m1k+&
          0.5_wp*d1_rhs_rho_dx_1_nxp2p0m2nyp2m1k

d1_rhs_rho_dx_1_nxp2p0nyp2m1k = d1_rhs_rho_dx_1_nxp2p0nyp2m1k*param_float(1)

d1_rhs_rho_dx_2_nxp2p0p0nyp2m1k = (param_float(3 + 5))*q(nx+2+0+0,ny+2-1,indvars(1))*((q(nx+2+0+0,ny+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0+0,ny+2-1,indvars(2))*q(nx+2+0+0,ny+2-1,indvars(2))+&
                    q(nx+2+0+0,ny+2-1,indvars(3))*q(nx+2+0+0,ny+2-1,indvars(3)))))

d1_rhs_rho_dx_2_nxp2p0m1nyp2m1k = (param_float(3 + 5))*q(nx+2+0-1,ny+2-1,indvars(1))*((q(nx+2+0-1,ny+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0-1,ny+2-1,indvars(2))*q(nx+2+0-1,ny+2-1,indvars(2))+&
                    q(nx+2+0-1,ny+2-1,indvars(3))*q(nx+2+0-1,ny+2-1,indvars(3)))))

d1_rhs_rho_dx_2_nxp2p0m2nyp2m1k = (param_float(3 + 5))*q(nx+2+0-2,ny+2-1,indvars(1))*((q(nx+2+0-2,ny+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0-2,ny+2-1,indvars(2))*q(nx+2+0-2,ny+2-1,indvars(2))+&
                    q(nx+2+0-2,ny+2-1,indvars(3))*q(nx+2+0-2,ny+2-1,indvars(3)))))

d1_rhs_rho_dx_2_nxp2p0nyp2m1k = 1.5_wp*d1_rhs_rho_dx_2_nxp2p0p0nyp2m1k-&
          2.0_wp*d1_rhs_rho_dx_2_nxp2p0m1nyp2m1k+&
          0.5_wp*d1_rhs_rho_dx_2_nxp2p0m2nyp2m1k

d1_rhs_rho_dx_2_nxp2p0nyp2m1k = d1_rhs_rho_dx_2_nxp2p0nyp2m1k*param_float(1)

d1_rhs_rho_dx_3_nxp2p0p0nyp2m1k = q(nx+2+0+0,ny+2-1,indvars(2))

d1_rhs_rho_dx_3_nxp2p0m1nyp2m1k = q(nx+2+0-1,ny+2-1,indvars(2))

d1_rhs_rho_dx_3_nxp2p0m2nyp2m1k = q(nx+2+0-2,ny+2-1,indvars(2))

d1_rhs_rho_dx_3_nxp2p0nyp2m1k = 1.5_wp*d1_rhs_rho_dx_3_nxp2p0p0nyp2m1k-&
          2.0_wp*d1_rhs_rho_dx_3_nxp2p0m1nyp2m1k+&
          0.5_wp*d1_rhs_rho_dx_3_nxp2p0m2nyp2m1k

d1_rhs_rho_dx_3_nxp2p0nyp2m1k = d1_rhs_rho_dx_3_nxp2p0nyp2m1k*param_float(1)

d1_rhs_rho_dy_0_nxp2p0nyp2m1m1k = q(nx+2+0,ny+2-1-1,indvars(1))*q(nx+2+0,ny+2-1-1,indvars(3))

d1_rhs_rho_dy_0_nxp2p0nyp2m1p1k = q(nx+2+0,ny+2-1+1,indvars(1))*q(nx+2+0,ny+2-1+1,indvars(3))

d1_rhs_rho_dy_0_nxp2p0nyp2m1k = -&
          0.5_wp*d1_rhs_rho_dy_0_nxp2p0nyp2m1m1k+&
          0.5_wp*d1_rhs_rho_dy_0_nxp2p0nyp2m1p1k

d1_rhs_rho_dy_0_nxp2p0nyp2m1k = d1_rhs_rho_dy_0_nxp2p0nyp2m1k*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer 0 1 None d(rho)/dt *************
!                                                           
!***********************************************************


rhs(nx+2+0,ny+2-1,indvars(1)) =   -  ( qst(nx+2+0,ny+2-1,indvarsst(10))*(d1_rhs_rho_dx_0_nxp2p0nyp2m1k)+&
                    qst(nx+2+0,ny+2-1,indvarsst(11))*(d1_rhs_rho_dy_0_nxp2p0nyp2m1k)+&
                    (1.0_wp/(param_float(23 + 5)*(param_float(3 + 5))*q(nx+2+0,ny+2-1,indvars(1))*((q(nx+2+0,ny+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2-1,indvars(2))*q(nx+2+0,ny+2-1,indvars(2))+&
                    q(nx+2+0,ny+2-1,indvars(3))*q(nx+2+0,ny+2-1,indvars(3)))))/q(nx+2+0,ny+2-1,indvars(1)))**2*((q(nx+2+0,ny+2-1,indvars(2))*(-&
                    (param_float(23 + 5)*(param_float(3 + 5))*q(nx+2+0,ny+2-1,indvars(1))*((q(nx+2+0,ny+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2-1,indvars(2))*q(nx+2+0,ny+2-1,indvars(2))+&
                    q(nx+2+0,ny+2-1,indvars(3))*q(nx+2+0,ny+2-1,indvars(3)))))/q(nx+2+0,ny+2-1,indvars(1)))**2*d1_rhs_rho_dx_1_nxp2p0nyp2m1k+&
                    1.0_wp*d1_rhs_rho_dx_2_nxp2p0nyp2m1k))*qst(nx+2+0,ny+2-1,indvarsst(10))+&
                    0.5_wp*((((param_float(23 + 5)*(param_float(3 + 5))*q(nx+2+0,ny+2-1,indvars(1))*((q(nx+2+0,ny+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2-1,indvars(2))*q(nx+2+0,ny+2-1,indvars(2))+&
                    q(nx+2+0,ny+2-1,indvars(3))*q(nx+2+0,ny+2-1,indvars(3)))))/q(nx+2+0,ny+2-1,indvars(1)))+&
                    q(nx+2+0,ny+2-1,indvars(2)))*((param_float(23 + 5)*(param_float(3 + 5))*q(nx+2+0,ny+2-1,indvars(1))*((q(nx+2+0,ny+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2-1,indvars(2))*q(nx+2+0,ny+2-1,indvars(2))+&
                    q(nx+2+0,ny+2-1,indvars(3))*q(nx+2+0,ny+2-1,indvars(3)))))/q(nx+2+0,ny+2-1,indvars(1)))*q(nx+2+0,ny+2-1,indvars(1))*d1_rhs_rho_dx_3_nxp2p0nyp2m1k+&
                    1.0_wp*d1_rhs_rho_dx_2_nxp2p0nyp2m1k))*qst(nx+2+0,ny+2-1,indvarsst(10))+&
                    -&
                    (((param_float(23 + 5)*(param_float(3 + 5))*q(nx+2+0,ny+2-1,indvars(1))*((q(nx+2+0,ny+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2-1,indvars(2))*q(nx+2+0,ny+2-1,indvars(2))+&
                    q(nx+2+0,ny+2-1,indvars(3))*q(nx+2+0,ny+2-1,indvars(3)))))/q(nx+2+0,ny+2-1,indvars(1)))+&
                    q(nx+2+0,ny+2-1,indvars(2)))*((param_float(23 + 5)*(param_float(3 + 5))*q(nx+2+0,ny+2-1,indvars(1))*((q(nx+2+0,ny+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2-1,indvars(2))*q(nx+2+0,ny+2-1,indvars(2))+&
                    q(nx+2+0,ny+2-1,indvars(3))*q(nx+2+0,ny+2-1,indvars(3)))))/q(nx+2+0,ny+2-1,indvars(1)))*q(nx+2+0,ny+2-1,indvars(1))*d1_rhs_rho_dx_3_nxp2p0nyp2m1k+&
                    1.0_wp*d1_rhs_rho_dx_2_nxp2p0nyp2m1k))*qst(nx+2+0,ny+2-1,indvarsst(10))))) ) 



!***********************************************************
!                                                           
! building source terms in RHS for layer 0 1 None d(rho u)/dt 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! u*(1.0_wp/c**2*((u*(-c**2*[rho]_1x+1.0_wp*[p]_1x))*deltaxI+0.5_wp*(((c+u)*(c*rho*[u]_1x+1.0_wp*[p]_1x))*deltaxI+-((c+u)*(c*rho*[u]_1x+1.0_wp*[p]_1x))*deltaxI)))+rho*(1.0_wp/(2*rho*c)*(((c+u)*(c*rho*[u]_1x+1.0_wp*[p]_1x))*deltaxI--((c+u)*(c*rho*[u]_1x+1.0_wp*[p]_1x))*deltaxI))+[rho*u*v]_1y+deltaxI*([-visc_t*(2.0_wp*deltaxI*({u}_1x)-2.0_wp/3.0_wp*(deltaxI*({u}_1x)+deltayI*({v}_1y)))]_1x)+deltayI*([-visc_t*(deltayI*({u}_1y)+deltaxI*({v}_1x))]_1y)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d2_rhs_u_dxdx_10_0_nxp2p0p0nyp2m1k_nxp2p0p0p0nyp2m1k = q(nx+2+0+0+0,ny+2-1,indvars(2))

d2_rhs_u_dxdx_10_0_nxp2p0p0nyp2m1k_nxp2p0p0m1nyp2m1k = q(nx+2+0+0-1,ny+2-1,indvars(2))

d2_rhs_u_dxdx_10_0_nxp2p0p0nyp2m1k_nxp2p0p0m2nyp2m1k = q(nx+2+0+0-2,ny+2-1,indvars(2))

d2_rhs_u_dxdx_10_0_nxp2p0p0nyp2m1k = 1.5_wp*d2_rhs_u_dxdx_10_0_nxp2p0p0nyp2m1k_nxp2p0p0p0nyp2m1k-&
          2.0_wp*d2_rhs_u_dxdx_10_0_nxp2p0p0nyp2m1k_nxp2p0p0m1nyp2m1k+&
          0.5_wp*d2_rhs_u_dxdx_10_0_nxp2p0p0nyp2m1k_nxp2p0p0m2nyp2m1k

d2_rhs_u_dxdx_10_0_nxp2p0p0nyp2m1k = d2_rhs_u_dxdx_10_0_nxp2p0p0nyp2m1k*param_float(1)

d2_rhs_u_dxdx_10_0_nxp2p0m1nyp2m1k_nxp2p0m1m1nyp2m1k = q(nx+2+0-1-1,ny+2-1,indvars(2))

d2_rhs_u_dxdx_10_0_nxp2p0m1nyp2m1k_nxp2p0m1p1nyp2m1k = q(nx+2+0-1+1,ny+2-1,indvars(2))

d2_rhs_u_dxdx_10_0_nxp2p0m1nyp2m1k = -&
          0.5_wp*d2_rhs_u_dxdx_10_0_nxp2p0m1nyp2m1k_nxp2p0m1m1nyp2m1k+&
          0.5_wp*d2_rhs_u_dxdx_10_0_nxp2p0m1nyp2m1k_nxp2p0m1p1nyp2m1k

d2_rhs_u_dxdx_10_0_nxp2p0m1nyp2m1k = d2_rhs_u_dxdx_10_0_nxp2p0m1nyp2m1k*param_float(1)

d2_rhs_u_dxdx_10_0_nxp2p0m2nyp2m1k_nxp2p0m2m1nyp2m1k = q(nx+2+0-2-1,ny+2-1,indvars(2))

d2_rhs_u_dxdx_10_0_nxp2p0m2nyp2m1k_nxp2p0m2p1nyp2m1k = q(nx+2+0-2+1,ny+2-1,indvars(2))

d2_rhs_u_dxdx_10_0_nxp2p0m2nyp2m1k = -&
          0.5_wp*d2_rhs_u_dxdx_10_0_nxp2p0m2nyp2m1k_nxp2p0m2m1nyp2m1k+&
          0.5_wp*d2_rhs_u_dxdx_10_0_nxp2p0m2nyp2m1k_nxp2p0m2p1nyp2m1k

d2_rhs_u_dxdx_10_0_nxp2p0m2nyp2m1k = d2_rhs_u_dxdx_10_0_nxp2p0m2nyp2m1k*param_float(1)

d2_rhs_u_dxdy_10_0_nxp2p0p0nyp2m1k_nxp2p0p0nyp2m1m1k = q(nx+2+0+0,ny+2-1-1,indvars(3))

d2_rhs_u_dxdy_10_0_nxp2p0p0nyp2m1k_nxp2p0p0nyp2m1p1k = q(nx+2+0+0,ny+2-1+1,indvars(3))

d2_rhs_u_dxdy_10_0_nxp2p0p0nyp2m1k = -&
          0.5_wp*d2_rhs_u_dxdy_10_0_nxp2p0p0nyp2m1k_nxp2p0p0nyp2m1m1k+&
          0.5_wp*d2_rhs_u_dxdy_10_0_nxp2p0p0nyp2m1k_nxp2p0p0nyp2m1p1k

d2_rhs_u_dxdy_10_0_nxp2p0p0nyp2m1k = d2_rhs_u_dxdy_10_0_nxp2p0p0nyp2m1k*param_float(2)

d2_rhs_u_dxdy_10_0_nxp2p0m1nyp2m1k_nxp2p0m1nyp2m1m1k = q(nx+2+0-1,ny+2-1-1,indvars(3))

d2_rhs_u_dxdy_10_0_nxp2p0m1nyp2m1k_nxp2p0m1nyp2m1p1k = q(nx+2+0-1,ny+2-1+1,indvars(3))

d2_rhs_u_dxdy_10_0_nxp2p0m1nyp2m1k = -&
          0.5_wp*d2_rhs_u_dxdy_10_0_nxp2p0m1nyp2m1k_nxp2p0m1nyp2m1m1k+&
          0.5_wp*d2_rhs_u_dxdy_10_0_nxp2p0m1nyp2m1k_nxp2p0m1nyp2m1p1k

d2_rhs_u_dxdy_10_0_nxp2p0m1nyp2m1k = d2_rhs_u_dxdy_10_0_nxp2p0m1nyp2m1k*param_float(2)

d2_rhs_u_dxdy_10_0_nxp2p0m2nyp2m1k_nxp2p0m2nyp2m1m1k = q(nx+2+0-2,ny+2-1-1,indvars(3))

d2_rhs_u_dxdy_10_0_nxp2p0m2nyp2m1k_nxp2p0m2nyp2m1p1k = q(nx+2+0-2,ny+2-1+1,indvars(3))

d2_rhs_u_dxdy_10_0_nxp2p0m2nyp2m1k = -&
          0.5_wp*d2_rhs_u_dxdy_10_0_nxp2p0m2nyp2m1k_nxp2p0m2nyp2m1m1k+&
          0.5_wp*d2_rhs_u_dxdy_10_0_nxp2p0m2nyp2m1k_nxp2p0m2nyp2m1p1k

d2_rhs_u_dxdy_10_0_nxp2p0m2nyp2m1k = d2_rhs_u_dxdy_10_0_nxp2p0m2nyp2m1k*param_float(2)

d1_rhs_u_dx_10_nxp2p0p0nyp2m1k = -((1+&
                    param_float(21 + 5))/(((q(nx+2+0+0,ny+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0+0,ny+2-1,indvars(2))*q(nx+2+0+0,ny+2-1,indvars(2))+&
                    q(nx+2+0+0,ny+2-1,indvars(3))*q(nx+2+0+0,ny+2-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2+0+0,ny+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0+0,ny+2-1,indvars(2))*q(nx+2+0+0,ny+2-1,indvars(2))+&
                    q(nx+2+0+0,ny+2-1,indvars(3))*q(nx+2+0+0,ny+2-1,indvars(3)))))/param_float(4 + 5)**1.5)*(1+&
                    ((q(nx+2+0+0,ny+2-1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+2+0+0,ny+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0+0,ny+2-1,indvars(2))*q(nx+2+0+0,ny+2-1,indvars(2))+&
                    q(nx+2+0+0,ny+2-1,indvars(3))*q(nx+2+0+0,ny+2-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2+0+0,ny+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0+0,ny+2-1,indvars(2))*q(nx+2+0+0,ny+2-1,indvars(2))+&
                    q(nx+2+0+0,ny+2-1,indvars(3))*q(nx+2+0+0,ny+2-1,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+2+0+0,ny+2-1,indvars(1)))**3/((q(nx+2+0+0,ny+2-1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+2+0+0,ny+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0+0,ny+2-1,indvars(2))*q(nx+2+0+0,ny+2-1,indvars(2))+&
                    q(nx+2+0+0,ny+2-1,indvars(3))*q(nx+2+0+0,ny+2-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2+0+0,ny+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0+0,ny+2-1,indvars(2))*q(nx+2+0+0,ny+2-1,indvars(2))+&
                    q(nx+2+0+0,ny+2-1,indvars(3))*q(nx+2+0+0,ny+2-1,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+2+0+0,ny+2-1,indvars(1)))**3+&
                    param_float(13 + 5)**3))*(q(nx+2+0+0,ny+2-1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+2+0+0,ny+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0+0,ny+2-1,indvars(2))*q(nx+2+0+0,ny+2-1,indvars(2))+&
                    q(nx+2+0+0,ny+2-1,indvars(3))*q(nx+2+0+0,ny+2-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2+0+0,ny+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0+0,ny+2-1,indvars(2))*q(nx+2+0+0,ny+2-1,indvars(2))+&
                    q(nx+2+0+0,ny+2-1,indvars(3))*q(nx+2+0+0,ny+2-1,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+2+0+0,ny+2-1,indvars(1))))*param_float(1 + 5)*(2.0_wp*qst(nx+2+0+0,ny+2-1,indvarsst(10))*(d2_rhs_u_dxdx_10_0_nxp2p0p0nyp2m1k)-&
                    2.0_wp/3.0_wp*(qst(nx+2+0+0,ny+2-1,indvarsst(10))*(d2_rhs_u_dxdx_10_0_nxp2p0p0nyp2m1k)+&
                    qst(nx+2+0+0,ny+2-1,indvarsst(11))*(d2_rhs_u_dxdy_10_0_nxp2p0p0nyp2m1k)))

d1_rhs_u_dx_10_nxp2p0m1nyp2m1k = -((1+&
                    param_float(21 + 5))/(((q(nx+2+0-1,ny+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0-1,ny+2-1,indvars(2))*q(nx+2+0-1,ny+2-1,indvars(2))+&
                    q(nx+2+0-1,ny+2-1,indvars(3))*q(nx+2+0-1,ny+2-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2+0-1,ny+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0-1,ny+2-1,indvars(2))*q(nx+2+0-1,ny+2-1,indvars(2))+&
                    q(nx+2+0-1,ny+2-1,indvars(3))*q(nx+2+0-1,ny+2-1,indvars(3)))))/param_float(4 + 5)**1.5)*(1+&
                    ((q(nx+2+0-1,ny+2-1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+2+0-1,ny+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0-1,ny+2-1,indvars(2))*q(nx+2+0-1,ny+2-1,indvars(2))+&
                    q(nx+2+0-1,ny+2-1,indvars(3))*q(nx+2+0-1,ny+2-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2+0-1,ny+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0-1,ny+2-1,indvars(2))*q(nx+2+0-1,ny+2-1,indvars(2))+&
                    q(nx+2+0-1,ny+2-1,indvars(3))*q(nx+2+0-1,ny+2-1,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+2+0-1,ny+2-1,indvars(1)))**3/((q(nx+2+0-1,ny+2-1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+2+0-1,ny+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0-1,ny+2-1,indvars(2))*q(nx+2+0-1,ny+2-1,indvars(2))+&
                    q(nx+2+0-1,ny+2-1,indvars(3))*q(nx+2+0-1,ny+2-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2+0-1,ny+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0-1,ny+2-1,indvars(2))*q(nx+2+0-1,ny+2-1,indvars(2))+&
                    q(nx+2+0-1,ny+2-1,indvars(3))*q(nx+2+0-1,ny+2-1,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+2+0-1,ny+2-1,indvars(1)))**3+&
                    param_float(13 + 5)**3))*(q(nx+2+0-1,ny+2-1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+2+0-1,ny+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0-1,ny+2-1,indvars(2))*q(nx+2+0-1,ny+2-1,indvars(2))+&
                    q(nx+2+0-1,ny+2-1,indvars(3))*q(nx+2+0-1,ny+2-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2+0-1,ny+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0-1,ny+2-1,indvars(2))*q(nx+2+0-1,ny+2-1,indvars(2))+&
                    q(nx+2+0-1,ny+2-1,indvars(3))*q(nx+2+0-1,ny+2-1,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+2+0-1,ny+2-1,indvars(1))))*param_float(1 + 5)*(2.0_wp*qst(nx+2+0-1,ny+2-1,indvarsst(10))*(d2_rhs_u_dxdx_10_0_nxp2p0m1nyp2m1k)-&
                    2.0_wp/3.0_wp*(qst(nx+2+0-1,ny+2-1,indvarsst(10))*(d2_rhs_u_dxdx_10_0_nxp2p0m1nyp2m1k)+&
                    qst(nx+2+0-1,ny+2-1,indvarsst(11))*(d2_rhs_u_dxdy_10_0_nxp2p0m1nyp2m1k)))

d1_rhs_u_dx_10_nxp2p0m2nyp2m1k = -((1+&
                    param_float(21 + 5))/(((q(nx+2+0-2,ny+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0-2,ny+2-1,indvars(2))*q(nx+2+0-2,ny+2-1,indvars(2))+&
                    q(nx+2+0-2,ny+2-1,indvars(3))*q(nx+2+0-2,ny+2-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2+0-2,ny+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0-2,ny+2-1,indvars(2))*q(nx+2+0-2,ny+2-1,indvars(2))+&
                    q(nx+2+0-2,ny+2-1,indvars(3))*q(nx+2+0-2,ny+2-1,indvars(3)))))/param_float(4 + 5)**1.5)*(1+&
                    ((q(nx+2+0-2,ny+2-1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+2+0-2,ny+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0-2,ny+2-1,indvars(2))*q(nx+2+0-2,ny+2-1,indvars(2))+&
                    q(nx+2+0-2,ny+2-1,indvars(3))*q(nx+2+0-2,ny+2-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2+0-2,ny+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0-2,ny+2-1,indvars(2))*q(nx+2+0-2,ny+2-1,indvars(2))+&
                    q(nx+2+0-2,ny+2-1,indvars(3))*q(nx+2+0-2,ny+2-1,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+2+0-2,ny+2-1,indvars(1)))**3/((q(nx+2+0-2,ny+2-1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+2+0-2,ny+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0-2,ny+2-1,indvars(2))*q(nx+2+0-2,ny+2-1,indvars(2))+&
                    q(nx+2+0-2,ny+2-1,indvars(3))*q(nx+2+0-2,ny+2-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2+0-2,ny+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0-2,ny+2-1,indvars(2))*q(nx+2+0-2,ny+2-1,indvars(2))+&
                    q(nx+2+0-2,ny+2-1,indvars(3))*q(nx+2+0-2,ny+2-1,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+2+0-2,ny+2-1,indvars(1)))**3+&
                    param_float(13 + 5)**3))*(q(nx+2+0-2,ny+2-1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+2+0-2,ny+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0-2,ny+2-1,indvars(2))*q(nx+2+0-2,ny+2-1,indvars(2))+&
                    q(nx+2+0-2,ny+2-1,indvars(3))*q(nx+2+0-2,ny+2-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2+0-2,ny+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0-2,ny+2-1,indvars(2))*q(nx+2+0-2,ny+2-1,indvars(2))+&
                    q(nx+2+0-2,ny+2-1,indvars(3))*q(nx+2+0-2,ny+2-1,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+2+0-2,ny+2-1,indvars(1))))*param_float(1 + 5)*(2.0_wp*qst(nx+2+0-2,ny+2-1,indvarsst(10))*(d2_rhs_u_dxdx_10_0_nxp2p0m2nyp2m1k)-&
                    2.0_wp/3.0_wp*(qst(nx+2+0-2,ny+2-1,indvarsst(10))*(d2_rhs_u_dxdx_10_0_nxp2p0m2nyp2m1k)+&
                    qst(nx+2+0-2,ny+2-1,indvarsst(11))*(d2_rhs_u_dxdy_10_0_nxp2p0m2nyp2m1k)))

d1_rhs_u_dx_10_nxp2p0nyp2m1k = 1.5_wp*d1_rhs_u_dx_10_nxp2p0p0nyp2m1k-&
          2.0_wp*d1_rhs_u_dx_10_nxp2p0m1nyp2m1k+&
          0.5_wp*d1_rhs_u_dx_10_nxp2p0m2nyp2m1k

d1_rhs_u_dx_10_nxp2p0nyp2m1k = d1_rhs_u_dx_10_nxp2p0nyp2m1k*param_float(1)

d1_rhs_u_dy_0_nxp2p0nyp2m1m1k = q(nx+2+0,ny+2-1-1,indvars(1))*q(nx+2+0,ny+2-1-1,indvars(2))*q(nx+2+0,ny+2-1-1,indvars(3))

d1_rhs_u_dy_0_nxp2p0nyp2m1p1k = q(nx+2+0,ny+2-1+1,indvars(1))*q(nx+2+0,ny+2-1+1,indvars(2))*q(nx+2+0,ny+2-1+1,indvars(3))

d1_rhs_u_dy_0_nxp2p0nyp2m1k = -&
          0.5_wp*d1_rhs_u_dy_0_nxp2p0nyp2m1m1k+&
          0.5_wp*d1_rhs_u_dy_0_nxp2p0nyp2m1p1k

d1_rhs_u_dy_0_nxp2p0nyp2m1k = d1_rhs_u_dy_0_nxp2p0nyp2m1k*param_float(2)

d2_rhs_u_dydx_1_0_nxp2p0nyp2m1m1k_nxp2p0p0nyp2m1m1k = q(nx+2+0+0,ny+2-1-1,indvars(3))

d2_rhs_u_dydx_1_0_nxp2p0nyp2m1m1k_nxp2p0m1nyp2m1m1k = q(nx+2+0-1,ny+2-1-1,indvars(3))

d2_rhs_u_dydx_1_0_nxp2p0nyp2m1m1k_nxp2p0m2nyp2m1m1k = q(nx+2+0-2,ny+2-1-1,indvars(3))

d2_rhs_u_dydx_1_0_nxp2p0nyp2m1m1k = 1.5_wp*d2_rhs_u_dydx_1_0_nxp2p0nyp2m1m1k_nxp2p0p0nyp2m1m1k-&
          2.0_wp*d2_rhs_u_dydx_1_0_nxp2p0nyp2m1m1k_nxp2p0m1nyp2m1m1k+&
          0.5_wp*d2_rhs_u_dydx_1_0_nxp2p0nyp2m1m1k_nxp2p0m2nyp2m1m1k

d2_rhs_u_dydx_1_0_nxp2p0nyp2m1m1k = d2_rhs_u_dydx_1_0_nxp2p0nyp2m1m1k*param_float(1)

d2_rhs_u_dydx_1_0_nxp2p0nyp2m1p1k_nxp2p0p0nyp2m1p1k = q(nx+2+0+0,ny+2-1+1,indvars(3))

d2_rhs_u_dydx_1_0_nxp2p0nyp2m1p1k_nxp2p0m1nyp2m1p1k = q(nx+2+0-1,ny+2-1+1,indvars(3))

d2_rhs_u_dydx_1_0_nxp2p0nyp2m1p1k_nxp2p0m2nyp2m1p1k = q(nx+2+0-2,ny+2-1+1,indvars(3))

d2_rhs_u_dydx_1_0_nxp2p0nyp2m1p1k = 1.5_wp*d2_rhs_u_dydx_1_0_nxp2p0nyp2m1p1k_nxp2p0p0nyp2m1p1k-&
          2.0_wp*d2_rhs_u_dydx_1_0_nxp2p0nyp2m1p1k_nxp2p0m1nyp2m1p1k+&
          0.5_wp*d2_rhs_u_dydx_1_0_nxp2p0nyp2m1p1k_nxp2p0m2nyp2m1p1k

d2_rhs_u_dydx_1_0_nxp2p0nyp2m1p1k = d2_rhs_u_dydx_1_0_nxp2p0nyp2m1p1k*param_float(1)

d2_rhs_u_dydy_1_0_nxp2p0nyp2m1m1k_nxp2p0nyp2m1m1m1k = q(nx+2+0,ny+2-1-1-1,indvars(2))

d2_rhs_u_dydy_1_0_nxp2p0nyp2m1m1k_nxp2p0nyp2m1m1p1k = q(nx+2+0,ny+2-1-1+1,indvars(2))

d2_rhs_u_dydy_1_0_nxp2p0nyp2m1m1k = -&
          0.5_wp*d2_rhs_u_dydy_1_0_nxp2p0nyp2m1m1k_nxp2p0nyp2m1m1m1k+&
          0.5_wp*d2_rhs_u_dydy_1_0_nxp2p0nyp2m1m1k_nxp2p0nyp2m1m1p1k

d2_rhs_u_dydy_1_0_nxp2p0nyp2m1m1k = d2_rhs_u_dydy_1_0_nxp2p0nyp2m1m1k*param_float(2)

d2_rhs_u_dydy_1_0_nxp2p0nyp2m1p1k_nxp2p0nyp2m1p1p0k = q(nx+2+0,ny+2-1+1+0,indvars(2))

d2_rhs_u_dydy_1_0_nxp2p0nyp2m1p1k_nxp2p0nyp2m1p1m1k = q(nx+2+0,ny+2-1+1-1,indvars(2))

d2_rhs_u_dydy_1_0_nxp2p0nyp2m1p1k_nxp2p0nyp2m1p1m2k = q(nx+2+0,ny+2-1+1-2,indvars(2))

d2_rhs_u_dydy_1_0_nxp2p0nyp2m1p1k = 1.5_wp*d2_rhs_u_dydy_1_0_nxp2p0nyp2m1p1k_nxp2p0nyp2m1p1p0k-&
          2.0_wp*d2_rhs_u_dydy_1_0_nxp2p0nyp2m1p1k_nxp2p0nyp2m1p1m1k+&
          0.5_wp*d2_rhs_u_dydy_1_0_nxp2p0nyp2m1p1k_nxp2p0nyp2m1p1m2k

d2_rhs_u_dydy_1_0_nxp2p0nyp2m1p1k = d2_rhs_u_dydy_1_0_nxp2p0nyp2m1p1k*param_float(2)

d1_rhs_u_dy_1_nxp2p0nyp2m1m1k = -((1+&
                    param_float(21 + 5))/(((q(nx+2+0,ny+2-1-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2-1-1,indvars(2))*q(nx+2+0,ny+2-1-1,indvars(2))+&
                    q(nx+2+0,ny+2-1-1,indvars(3))*q(nx+2+0,ny+2-1-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2+0,ny+2-1-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2-1-1,indvars(2))*q(nx+2+0,ny+2-1-1,indvars(2))+&
                    q(nx+2+0,ny+2-1-1,indvars(3))*q(nx+2+0,ny+2-1-1,indvars(3)))))/param_float(4 + 5)**1.5)*(1+&
                    ((q(nx+2+0,ny+2-1-1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+2+0,ny+2-1-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2-1-1,indvars(2))*q(nx+2+0,ny+2-1-1,indvars(2))+&
                    q(nx+2+0,ny+2-1-1,indvars(3))*q(nx+2+0,ny+2-1-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2+0,ny+2-1-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2-1-1,indvars(2))*q(nx+2+0,ny+2-1-1,indvars(2))+&
                    q(nx+2+0,ny+2-1-1,indvars(3))*q(nx+2+0,ny+2-1-1,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+2+0,ny+2-1-1,indvars(1)))**3/((q(nx+2+0,ny+2-1-1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+2+0,ny+2-1-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2-1-1,indvars(2))*q(nx+2+0,ny+2-1-1,indvars(2))+&
                    q(nx+2+0,ny+2-1-1,indvars(3))*q(nx+2+0,ny+2-1-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2+0,ny+2-1-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2-1-1,indvars(2))*q(nx+2+0,ny+2-1-1,indvars(2))+&
                    q(nx+2+0,ny+2-1-1,indvars(3))*q(nx+2+0,ny+2-1-1,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+2+0,ny+2-1-1,indvars(1)))**3+&
                    param_float(13 + 5)**3))*(q(nx+2+0,ny+2-1-1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+2+0,ny+2-1-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2-1-1,indvars(2))*q(nx+2+0,ny+2-1-1,indvars(2))+&
                    q(nx+2+0,ny+2-1-1,indvars(3))*q(nx+2+0,ny+2-1-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2+0,ny+2-1-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2-1-1,indvars(2))*q(nx+2+0,ny+2-1-1,indvars(2))+&
                    q(nx+2+0,ny+2-1-1,indvars(3))*q(nx+2+0,ny+2-1-1,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+2+0,ny+2-1-1,indvars(1))))*param_float(1 + 5)*(qst(nx+2+0,ny+2-1-1,indvarsst(11))*(d2_rhs_u_dydy_1_0_nxp2p0nyp2m1m1k)+&
                    qst(nx+2+0,ny+2-1-1,indvarsst(10))*(d2_rhs_u_dydx_1_0_nxp2p0nyp2m1m1k))

d1_rhs_u_dy_1_nxp2p0nyp2m1p1k = -((1+&
                    param_float(21 + 5))/(((q(nx+2+0,ny+2-1+1,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2-1+1,indvars(2))*q(nx+2+0,ny+2-1+1,indvars(2))+&
                    q(nx+2+0,ny+2-1+1,indvars(3))*q(nx+2+0,ny+2-1+1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2+0,ny+2-1+1,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2-1+1,indvars(2))*q(nx+2+0,ny+2-1+1,indvars(2))+&
                    q(nx+2+0,ny+2-1+1,indvars(3))*q(nx+2+0,ny+2-1+1,indvars(3)))))/param_float(4 + 5)**1.5)*(1+&
                    ((q(nx+2+0,ny+2-1+1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+2+0,ny+2-1+1,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2-1+1,indvars(2))*q(nx+2+0,ny+2-1+1,indvars(2))+&
                    q(nx+2+0,ny+2-1+1,indvars(3))*q(nx+2+0,ny+2-1+1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2+0,ny+2-1+1,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2-1+1,indvars(2))*q(nx+2+0,ny+2-1+1,indvars(2))+&
                    q(nx+2+0,ny+2-1+1,indvars(3))*q(nx+2+0,ny+2-1+1,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+2+0,ny+2-1+1,indvars(1)))**3/((q(nx+2+0,ny+2-1+1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+2+0,ny+2-1+1,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2-1+1,indvars(2))*q(nx+2+0,ny+2-1+1,indvars(2))+&
                    q(nx+2+0,ny+2-1+1,indvars(3))*q(nx+2+0,ny+2-1+1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2+0,ny+2-1+1,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2-1+1,indvars(2))*q(nx+2+0,ny+2-1+1,indvars(2))+&
                    q(nx+2+0,ny+2-1+1,indvars(3))*q(nx+2+0,ny+2-1+1,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+2+0,ny+2-1+1,indvars(1)))**3+&
                    param_float(13 + 5)**3))*(q(nx+2+0,ny+2-1+1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+2+0,ny+2-1+1,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2-1+1,indvars(2))*q(nx+2+0,ny+2-1+1,indvars(2))+&
                    q(nx+2+0,ny+2-1+1,indvars(3))*q(nx+2+0,ny+2-1+1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2+0,ny+2-1+1,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2-1+1,indvars(2))*q(nx+2+0,ny+2-1+1,indvars(2))+&
                    q(nx+2+0,ny+2-1+1,indvars(3))*q(nx+2+0,ny+2-1+1,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+2+0,ny+2-1+1,indvars(1))))*param_float(1 + 5)*(qst(nx+2+0,ny+2-1+1,indvarsst(11))*(d2_rhs_u_dydy_1_0_nxp2p0nyp2m1p1k)+&
                    qst(nx+2+0,ny+2-1+1,indvarsst(10))*(d2_rhs_u_dydx_1_0_nxp2p0nyp2m1p1k))

d1_rhs_u_dy_1_nxp2p0nyp2m1k = -&
          0.5_wp*d1_rhs_u_dy_1_nxp2p0nyp2m1m1k+&
          0.5_wp*d1_rhs_u_dy_1_nxp2p0nyp2m1p1k

d1_rhs_u_dy_1_nxp2p0nyp2m1k = d1_rhs_u_dy_1_nxp2p0nyp2m1k*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer 0 1 None d(rho u)/dt ***********
!                                                           
!***********************************************************


rhs(nx+2+0,ny+2-1,indvars(2)) =   -  ( q(nx+2+0,ny+2-1,indvars(2))*(1.0_wp/(param_float(23 + 5)*(param_float(3 + 5))*q(nx+2+0,ny+2-1,indvars(1))*((q(nx+2+0,ny+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2-1,indvars(2))*q(nx+2+0,ny+2-1,indvars(2))+&
                    q(nx+2+0,ny+2-1,indvars(3))*q(nx+2+0,ny+2-1,indvars(3)))))/q(nx+2+0,ny+2-1,indvars(1)))**2*((q(nx+2+0,ny+2-1,indvars(2))*(-&
                    (param_float(23 + 5)*(param_float(3 + 5))*q(nx+2+0,ny+2-1,indvars(1))*((q(nx+2+0,ny+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2-1,indvars(2))*q(nx+2+0,ny+2-1,indvars(2))+&
                    q(nx+2+0,ny+2-1,indvars(3))*q(nx+2+0,ny+2-1,indvars(3)))))/q(nx+2+0,ny+2-1,indvars(1)))**2*d1_rhs_rho_dx_1_nxp2p0nyp2m1k+&
                    1.0_wp*d1_rhs_rho_dx_2_nxp2p0nyp2m1k))*qst(nx+2+0,ny+2-1,indvarsst(10))+&
                    0.5_wp*((((param_float(23 + 5)*(param_float(3 + 5))*q(nx+2+0,ny+2-1,indvars(1))*((q(nx+2+0,ny+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2-1,indvars(2))*q(nx+2+0,ny+2-1,indvars(2))+&
                    q(nx+2+0,ny+2-1,indvars(3))*q(nx+2+0,ny+2-1,indvars(3)))))/q(nx+2+0,ny+2-1,indvars(1)))+&
                    q(nx+2+0,ny+2-1,indvars(2)))*((param_float(23 + 5)*(param_float(3 + 5))*q(nx+2+0,ny+2-1,indvars(1))*((q(nx+2+0,ny+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2-1,indvars(2))*q(nx+2+0,ny+2-1,indvars(2))+&
                    q(nx+2+0,ny+2-1,indvars(3))*q(nx+2+0,ny+2-1,indvars(3)))))/q(nx+2+0,ny+2-1,indvars(1)))*q(nx+2+0,ny+2-1,indvars(1))*d1_rhs_rho_dx_3_nxp2p0nyp2m1k+&
                    1.0_wp*d1_rhs_rho_dx_2_nxp2p0nyp2m1k))*qst(nx+2+0,ny+2-1,indvarsst(10))+&
                    -&
                    (((param_float(23 + 5)*(param_float(3 + 5))*q(nx+2+0,ny+2-1,indvars(1))*((q(nx+2+0,ny+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2-1,indvars(2))*q(nx+2+0,ny+2-1,indvars(2))+&
                    q(nx+2+0,ny+2-1,indvars(3))*q(nx+2+0,ny+2-1,indvars(3)))))/q(nx+2+0,ny+2-1,indvars(1)))+&
                    q(nx+2+0,ny+2-1,indvars(2)))*((param_float(23 + 5)*(param_float(3 + 5))*q(nx+2+0,ny+2-1,indvars(1))*((q(nx+2+0,ny+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2-1,indvars(2))*q(nx+2+0,ny+2-1,indvars(2))+&
                    q(nx+2+0,ny+2-1,indvars(3))*q(nx+2+0,ny+2-1,indvars(3)))))/q(nx+2+0,ny+2-1,indvars(1)))*q(nx+2+0,ny+2-1,indvars(1))*d1_rhs_rho_dx_3_nxp2p0nyp2m1k+&
                    1.0_wp*d1_rhs_rho_dx_2_nxp2p0nyp2m1k))*qst(nx+2+0,ny+2-1,indvarsst(10)))))+&
                    q(nx+2+0,ny+2-1,indvars(1))*(1.0_wp/(2*q(nx+2+0,ny+2-1,indvars(1))*(param_float(23 + 5)*(param_float(3 + 5))*q(nx+2+0,ny+2-1,indvars(1))*((q(nx+2+0,ny+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2-1,indvars(2))*q(nx+2+0,ny+2-1,indvars(2))+&
                    q(nx+2+0,ny+2-1,indvars(3))*q(nx+2+0,ny+2-1,indvars(3)))))/q(nx+2+0,ny+2-1,indvars(1))))*((((param_float(23 + 5)*(param_float(3 + 5))*q(nx+2+0,ny+2-1,indvars(1))*((q(nx+2+0,ny+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2-1,indvars(2))*q(nx+2+0,ny+2-1,indvars(2))+&
                    q(nx+2+0,ny+2-1,indvars(3))*q(nx+2+0,ny+2-1,indvars(3)))))/q(nx+2+0,ny+2-1,indvars(1)))+&
                    q(nx+2+0,ny+2-1,indvars(2)))*((param_float(23 + 5)*(param_float(3 + 5))*q(nx+2+0,ny+2-1,indvars(1))*((q(nx+2+0,ny+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2-1,indvars(2))*q(nx+2+0,ny+2-1,indvars(2))+&
                    q(nx+2+0,ny+2-1,indvars(3))*q(nx+2+0,ny+2-1,indvars(3)))))/q(nx+2+0,ny+2-1,indvars(1)))*q(nx+2+0,ny+2-1,indvars(1))*d1_rhs_rho_dx_3_nxp2p0nyp2m1k+&
                    1.0_wp*d1_rhs_rho_dx_2_nxp2p0nyp2m1k))*qst(nx+2+0,ny+2-1,indvarsst(10))-&
                    -&
                    (((param_float(23 + 5)*(param_float(3 + 5))*q(nx+2+0,ny+2-1,indvars(1))*((q(nx+2+0,ny+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2-1,indvars(2))*q(nx+2+0,ny+2-1,indvars(2))+&
                    q(nx+2+0,ny+2-1,indvars(3))*q(nx+2+0,ny+2-1,indvars(3)))))/q(nx+2+0,ny+2-1,indvars(1)))+&
                    q(nx+2+0,ny+2-1,indvars(2)))*((param_float(23 + 5)*(param_float(3 + 5))*q(nx+2+0,ny+2-1,indvars(1))*((q(nx+2+0,ny+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2-1,indvars(2))*q(nx+2+0,ny+2-1,indvars(2))+&
                    q(nx+2+0,ny+2-1,indvars(3))*q(nx+2+0,ny+2-1,indvars(3)))))/q(nx+2+0,ny+2-1,indvars(1)))*q(nx+2+0,ny+2-1,indvars(1))*d1_rhs_rho_dx_3_nxp2p0nyp2m1k+&
                    1.0_wp*d1_rhs_rho_dx_2_nxp2p0nyp2m1k))*qst(nx+2+0,ny+2-1,indvarsst(10))))+&
                    d1_rhs_u_dy_0_nxp2p0nyp2m1k+&
                    qst(nx+2+0,ny+2-1,indvarsst(10))*(d1_rhs_u_dx_10_nxp2p0nyp2m1k)+&
                    qst(nx+2+0,ny+2-1,indvarsst(11))*(d1_rhs_u_dy_1_nxp2p0nyp2m1k) ) 



!***********************************************************
!                                                           
! building source terms in RHS for layer 0 1 None d(rho v)/dt 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! v*(1.0_wp/c**2*((u*(-c**2*[rho]_1x+1.0_wp*[p]_1x))*deltaxI+0.5_wp*(((c+u)*(c*rho*[u]_1x+1.0_wp*[p]_1x))*deltaxI+-((c+u)*(c*rho*[u]_1x+1.0_wp*[p]_1x))*deltaxI)))+rho*((1.0_wp*u*[v]_1x)*deltaxI)+[rho*v*v]_1y+deltaxI*([-visc_t*(deltayI*({u}_1y)+deltaxI*({v}_1x))]_1x)+deltayI*([-visc_t*(2.0_wp*deltayI*({v}_1y)-2.0_wp/3.0_wp*(deltaxI*({u}_1x)+deltayI*({v}_1y)))]_1y)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_rhs_v_dx_6_nxp2p0p0nyp2m1k = q(nx+2+0+0,ny+2-1,indvars(3))

d1_rhs_v_dx_6_nxp2p0m1nyp2m1k = q(nx+2+0-1,ny+2-1,indvars(3))

d1_rhs_v_dx_6_nxp2p0m2nyp2m1k = q(nx+2+0-2,ny+2-1,indvars(3))

d1_rhs_v_dx_6_nxp2p0nyp2m1k = 1.5_wp*d1_rhs_v_dx_6_nxp2p0p0nyp2m1k-&
          2.0_wp*d1_rhs_v_dx_6_nxp2p0m1nyp2m1k+&
          0.5_wp*d1_rhs_v_dx_6_nxp2p0m2nyp2m1k

d1_rhs_v_dx_6_nxp2p0nyp2m1k = d1_rhs_v_dx_6_nxp2p0nyp2m1k*param_float(1)

d2_rhs_v_dxdx_7_0_nxp2p0p0nyp2m1k_nxp2p0p0p0nyp2m1k = q(nx+2+0+0+0,ny+2-1,indvars(3))

d2_rhs_v_dxdx_7_0_nxp2p0p0nyp2m1k_nxp2p0p0m1nyp2m1k = q(nx+2+0+0-1,ny+2-1,indvars(3))

d2_rhs_v_dxdx_7_0_nxp2p0p0nyp2m1k_nxp2p0p0m2nyp2m1k = q(nx+2+0+0-2,ny+2-1,indvars(3))

d2_rhs_v_dxdx_7_0_nxp2p0p0nyp2m1k = 1.5_wp*d2_rhs_v_dxdx_7_0_nxp2p0p0nyp2m1k_nxp2p0p0p0nyp2m1k-&
          2.0_wp*d2_rhs_v_dxdx_7_0_nxp2p0p0nyp2m1k_nxp2p0p0m1nyp2m1k+&
          0.5_wp*d2_rhs_v_dxdx_7_0_nxp2p0p0nyp2m1k_nxp2p0p0m2nyp2m1k

d2_rhs_v_dxdx_7_0_nxp2p0p0nyp2m1k = d2_rhs_v_dxdx_7_0_nxp2p0p0nyp2m1k*param_float(1)

d2_rhs_v_dxdx_7_0_nxp2p0m1nyp2m1k_nxp2p0m1m1nyp2m1k = q(nx+2+0-1-1,ny+2-1,indvars(3))

d2_rhs_v_dxdx_7_0_nxp2p0m1nyp2m1k_nxp2p0m1p1nyp2m1k = q(nx+2+0-1+1,ny+2-1,indvars(3))

d2_rhs_v_dxdx_7_0_nxp2p0m1nyp2m1k = -&
          0.5_wp*d2_rhs_v_dxdx_7_0_nxp2p0m1nyp2m1k_nxp2p0m1m1nyp2m1k+&
          0.5_wp*d2_rhs_v_dxdx_7_0_nxp2p0m1nyp2m1k_nxp2p0m1p1nyp2m1k

d2_rhs_v_dxdx_7_0_nxp2p0m1nyp2m1k = d2_rhs_v_dxdx_7_0_nxp2p0m1nyp2m1k*param_float(1)

d2_rhs_v_dxdx_7_0_nxp2p0m2nyp2m1k_nxp2p0m2m1nyp2m1k = q(nx+2+0-2-1,ny+2-1,indvars(3))

d2_rhs_v_dxdx_7_0_nxp2p0m2nyp2m1k_nxp2p0m2p1nyp2m1k = q(nx+2+0-2+1,ny+2-1,indvars(3))

d2_rhs_v_dxdx_7_0_nxp2p0m2nyp2m1k = -&
          0.5_wp*d2_rhs_v_dxdx_7_0_nxp2p0m2nyp2m1k_nxp2p0m2m1nyp2m1k+&
          0.5_wp*d2_rhs_v_dxdx_7_0_nxp2p0m2nyp2m1k_nxp2p0m2p1nyp2m1k

d2_rhs_v_dxdx_7_0_nxp2p0m2nyp2m1k = d2_rhs_v_dxdx_7_0_nxp2p0m2nyp2m1k*param_float(1)

d2_rhs_v_dxdy_7_0_nxp2p0p0nyp2m1k_nxp2p0p0nyp2m1m1k = q(nx+2+0+0,ny+2-1-1,indvars(2))

d2_rhs_v_dxdy_7_0_nxp2p0p0nyp2m1k_nxp2p0p0nyp2m1p1k = q(nx+2+0+0,ny+2-1+1,indvars(2))

d2_rhs_v_dxdy_7_0_nxp2p0p0nyp2m1k = -&
          0.5_wp*d2_rhs_v_dxdy_7_0_nxp2p0p0nyp2m1k_nxp2p0p0nyp2m1m1k+&
          0.5_wp*d2_rhs_v_dxdy_7_0_nxp2p0p0nyp2m1k_nxp2p0p0nyp2m1p1k

d2_rhs_v_dxdy_7_0_nxp2p0p0nyp2m1k = d2_rhs_v_dxdy_7_0_nxp2p0p0nyp2m1k*param_float(2)

d2_rhs_v_dxdy_7_0_nxp2p0m1nyp2m1k_nxp2p0m1nyp2m1m1k = q(nx+2+0-1,ny+2-1-1,indvars(2))

d2_rhs_v_dxdy_7_0_nxp2p0m1nyp2m1k_nxp2p0m1nyp2m1p1k = q(nx+2+0-1,ny+2-1+1,indvars(2))

d2_rhs_v_dxdy_7_0_nxp2p0m1nyp2m1k = -&
          0.5_wp*d2_rhs_v_dxdy_7_0_nxp2p0m1nyp2m1k_nxp2p0m1nyp2m1m1k+&
          0.5_wp*d2_rhs_v_dxdy_7_0_nxp2p0m1nyp2m1k_nxp2p0m1nyp2m1p1k

d2_rhs_v_dxdy_7_0_nxp2p0m1nyp2m1k = d2_rhs_v_dxdy_7_0_nxp2p0m1nyp2m1k*param_float(2)

d2_rhs_v_dxdy_7_0_nxp2p0m2nyp2m1k_nxp2p0m2nyp2m1m1k = q(nx+2+0-2,ny+2-1-1,indvars(2))

d2_rhs_v_dxdy_7_0_nxp2p0m2nyp2m1k_nxp2p0m2nyp2m1p1k = q(nx+2+0-2,ny+2-1+1,indvars(2))

d2_rhs_v_dxdy_7_0_nxp2p0m2nyp2m1k = -&
          0.5_wp*d2_rhs_v_dxdy_7_0_nxp2p0m2nyp2m1k_nxp2p0m2nyp2m1m1k+&
          0.5_wp*d2_rhs_v_dxdy_7_0_nxp2p0m2nyp2m1k_nxp2p0m2nyp2m1p1k

d2_rhs_v_dxdy_7_0_nxp2p0m2nyp2m1k = d2_rhs_v_dxdy_7_0_nxp2p0m2nyp2m1k*param_float(2)

d1_rhs_v_dx_7_nxp2p0p0nyp2m1k = -((1+&
                    param_float(21 + 5))/(((q(nx+2+0+0,ny+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0+0,ny+2-1,indvars(2))*q(nx+2+0+0,ny+2-1,indvars(2))+&
                    q(nx+2+0+0,ny+2-1,indvars(3))*q(nx+2+0+0,ny+2-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2+0+0,ny+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0+0,ny+2-1,indvars(2))*q(nx+2+0+0,ny+2-1,indvars(2))+&
                    q(nx+2+0+0,ny+2-1,indvars(3))*q(nx+2+0+0,ny+2-1,indvars(3)))))/param_float(4 + 5)**1.5)*(1+&
                    ((q(nx+2+0+0,ny+2-1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+2+0+0,ny+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0+0,ny+2-1,indvars(2))*q(nx+2+0+0,ny+2-1,indvars(2))+&
                    q(nx+2+0+0,ny+2-1,indvars(3))*q(nx+2+0+0,ny+2-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2+0+0,ny+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0+0,ny+2-1,indvars(2))*q(nx+2+0+0,ny+2-1,indvars(2))+&
                    q(nx+2+0+0,ny+2-1,indvars(3))*q(nx+2+0+0,ny+2-1,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+2+0+0,ny+2-1,indvars(1)))**3/((q(nx+2+0+0,ny+2-1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+2+0+0,ny+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0+0,ny+2-1,indvars(2))*q(nx+2+0+0,ny+2-1,indvars(2))+&
                    q(nx+2+0+0,ny+2-1,indvars(3))*q(nx+2+0+0,ny+2-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2+0+0,ny+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0+0,ny+2-1,indvars(2))*q(nx+2+0+0,ny+2-1,indvars(2))+&
                    q(nx+2+0+0,ny+2-1,indvars(3))*q(nx+2+0+0,ny+2-1,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+2+0+0,ny+2-1,indvars(1)))**3+&
                    param_float(13 + 5)**3))*(q(nx+2+0+0,ny+2-1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+2+0+0,ny+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0+0,ny+2-1,indvars(2))*q(nx+2+0+0,ny+2-1,indvars(2))+&
                    q(nx+2+0+0,ny+2-1,indvars(3))*q(nx+2+0+0,ny+2-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2+0+0,ny+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0+0,ny+2-1,indvars(2))*q(nx+2+0+0,ny+2-1,indvars(2))+&
                    q(nx+2+0+0,ny+2-1,indvars(3))*q(nx+2+0+0,ny+2-1,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+2+0+0,ny+2-1,indvars(1))))*param_float(1 + 5)*(qst(nx+2+0+0,ny+2-1,indvarsst(11))*(d2_rhs_v_dxdy_7_0_nxp2p0p0nyp2m1k)+&
                    qst(nx+2+0+0,ny+2-1,indvarsst(10))*(d2_rhs_v_dxdx_7_0_nxp2p0p0nyp2m1k))

d1_rhs_v_dx_7_nxp2p0m1nyp2m1k = -((1+&
                    param_float(21 + 5))/(((q(nx+2+0-1,ny+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0-1,ny+2-1,indvars(2))*q(nx+2+0-1,ny+2-1,indvars(2))+&
                    q(nx+2+0-1,ny+2-1,indvars(3))*q(nx+2+0-1,ny+2-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2+0-1,ny+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0-1,ny+2-1,indvars(2))*q(nx+2+0-1,ny+2-1,indvars(2))+&
                    q(nx+2+0-1,ny+2-1,indvars(3))*q(nx+2+0-1,ny+2-1,indvars(3)))))/param_float(4 + 5)**1.5)*(1+&
                    ((q(nx+2+0-1,ny+2-1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+2+0-1,ny+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0-1,ny+2-1,indvars(2))*q(nx+2+0-1,ny+2-1,indvars(2))+&
                    q(nx+2+0-1,ny+2-1,indvars(3))*q(nx+2+0-1,ny+2-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2+0-1,ny+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0-1,ny+2-1,indvars(2))*q(nx+2+0-1,ny+2-1,indvars(2))+&
                    q(nx+2+0-1,ny+2-1,indvars(3))*q(nx+2+0-1,ny+2-1,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+2+0-1,ny+2-1,indvars(1)))**3/((q(nx+2+0-1,ny+2-1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+2+0-1,ny+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0-1,ny+2-1,indvars(2))*q(nx+2+0-1,ny+2-1,indvars(2))+&
                    q(nx+2+0-1,ny+2-1,indvars(3))*q(nx+2+0-1,ny+2-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2+0-1,ny+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0-1,ny+2-1,indvars(2))*q(nx+2+0-1,ny+2-1,indvars(2))+&
                    q(nx+2+0-1,ny+2-1,indvars(3))*q(nx+2+0-1,ny+2-1,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+2+0-1,ny+2-1,indvars(1)))**3+&
                    param_float(13 + 5)**3))*(q(nx+2+0-1,ny+2-1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+2+0-1,ny+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0-1,ny+2-1,indvars(2))*q(nx+2+0-1,ny+2-1,indvars(2))+&
                    q(nx+2+0-1,ny+2-1,indvars(3))*q(nx+2+0-1,ny+2-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2+0-1,ny+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0-1,ny+2-1,indvars(2))*q(nx+2+0-1,ny+2-1,indvars(2))+&
                    q(nx+2+0-1,ny+2-1,indvars(3))*q(nx+2+0-1,ny+2-1,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+2+0-1,ny+2-1,indvars(1))))*param_float(1 + 5)*(qst(nx+2+0-1,ny+2-1,indvarsst(11))*(d2_rhs_v_dxdy_7_0_nxp2p0m1nyp2m1k)+&
                    qst(nx+2+0-1,ny+2-1,indvarsst(10))*(d2_rhs_v_dxdx_7_0_nxp2p0m1nyp2m1k))

d1_rhs_v_dx_7_nxp2p0m2nyp2m1k = -((1+&
                    param_float(21 + 5))/(((q(nx+2+0-2,ny+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0-2,ny+2-1,indvars(2))*q(nx+2+0-2,ny+2-1,indvars(2))+&
                    q(nx+2+0-2,ny+2-1,indvars(3))*q(nx+2+0-2,ny+2-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2+0-2,ny+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0-2,ny+2-1,indvars(2))*q(nx+2+0-2,ny+2-1,indvars(2))+&
                    q(nx+2+0-2,ny+2-1,indvars(3))*q(nx+2+0-2,ny+2-1,indvars(3)))))/param_float(4 + 5)**1.5)*(1+&
                    ((q(nx+2+0-2,ny+2-1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+2+0-2,ny+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0-2,ny+2-1,indvars(2))*q(nx+2+0-2,ny+2-1,indvars(2))+&
                    q(nx+2+0-2,ny+2-1,indvars(3))*q(nx+2+0-2,ny+2-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2+0-2,ny+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0-2,ny+2-1,indvars(2))*q(nx+2+0-2,ny+2-1,indvars(2))+&
                    q(nx+2+0-2,ny+2-1,indvars(3))*q(nx+2+0-2,ny+2-1,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+2+0-2,ny+2-1,indvars(1)))**3/((q(nx+2+0-2,ny+2-1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+2+0-2,ny+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0-2,ny+2-1,indvars(2))*q(nx+2+0-2,ny+2-1,indvars(2))+&
                    q(nx+2+0-2,ny+2-1,indvars(3))*q(nx+2+0-2,ny+2-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2+0-2,ny+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0-2,ny+2-1,indvars(2))*q(nx+2+0-2,ny+2-1,indvars(2))+&
                    q(nx+2+0-2,ny+2-1,indvars(3))*q(nx+2+0-2,ny+2-1,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+2+0-2,ny+2-1,indvars(1)))**3+&
                    param_float(13 + 5)**3))*(q(nx+2+0-2,ny+2-1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+2+0-2,ny+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0-2,ny+2-1,indvars(2))*q(nx+2+0-2,ny+2-1,indvars(2))+&
                    q(nx+2+0-2,ny+2-1,indvars(3))*q(nx+2+0-2,ny+2-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2+0-2,ny+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0-2,ny+2-1,indvars(2))*q(nx+2+0-2,ny+2-1,indvars(2))+&
                    q(nx+2+0-2,ny+2-1,indvars(3))*q(nx+2+0-2,ny+2-1,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+2+0-2,ny+2-1,indvars(1))))*param_float(1 + 5)*(qst(nx+2+0-2,ny+2-1,indvarsst(11))*(d2_rhs_v_dxdy_7_0_nxp2p0m2nyp2m1k)+&
                    qst(nx+2+0-2,ny+2-1,indvarsst(10))*(d2_rhs_v_dxdx_7_0_nxp2p0m2nyp2m1k))

d1_rhs_v_dx_7_nxp2p0nyp2m1k = 1.5_wp*d1_rhs_v_dx_7_nxp2p0p0nyp2m1k-&
          2.0_wp*d1_rhs_v_dx_7_nxp2p0m1nyp2m1k+&
          0.5_wp*d1_rhs_v_dx_7_nxp2p0m2nyp2m1k

d1_rhs_v_dx_7_nxp2p0nyp2m1k = d1_rhs_v_dx_7_nxp2p0nyp2m1k*param_float(1)

d1_rhs_v_dy_0_nxp2p0nyp2m1m1k = q(nx+2+0,ny+2-1-1,indvars(1))*q(nx+2+0,ny+2-1-1,indvars(3))*q(nx+2+0,ny+2-1-1,indvars(3))

d1_rhs_v_dy_0_nxp2p0nyp2m1p1k = q(nx+2+0,ny+2-1+1,indvars(1))*q(nx+2+0,ny+2-1+1,indvars(3))*q(nx+2+0,ny+2-1+1,indvars(3))

d1_rhs_v_dy_0_nxp2p0nyp2m1k = -&
          0.5_wp*d1_rhs_v_dy_0_nxp2p0nyp2m1m1k+&
          0.5_wp*d1_rhs_v_dy_0_nxp2p0nyp2m1p1k

d1_rhs_v_dy_0_nxp2p0nyp2m1k = d1_rhs_v_dy_0_nxp2p0nyp2m1k*param_float(2)

d2_rhs_v_dydx_1_0_nxp2p0nyp2m1m1k_nxp2p0p0nyp2m1m1k = q(nx+2+0+0,ny+2-1-1,indvars(2))

d2_rhs_v_dydx_1_0_nxp2p0nyp2m1m1k_nxp2p0m1nyp2m1m1k = q(nx+2+0-1,ny+2-1-1,indvars(2))

d2_rhs_v_dydx_1_0_nxp2p0nyp2m1m1k_nxp2p0m2nyp2m1m1k = q(nx+2+0-2,ny+2-1-1,indvars(2))

d2_rhs_v_dydx_1_0_nxp2p0nyp2m1m1k = 1.5_wp*d2_rhs_v_dydx_1_0_nxp2p0nyp2m1m1k_nxp2p0p0nyp2m1m1k-&
          2.0_wp*d2_rhs_v_dydx_1_0_nxp2p0nyp2m1m1k_nxp2p0m1nyp2m1m1k+&
          0.5_wp*d2_rhs_v_dydx_1_0_nxp2p0nyp2m1m1k_nxp2p0m2nyp2m1m1k

d2_rhs_v_dydx_1_0_nxp2p0nyp2m1m1k = d2_rhs_v_dydx_1_0_nxp2p0nyp2m1m1k*param_float(1)

d2_rhs_v_dydx_1_0_nxp2p0nyp2m1p1k_nxp2p0p0nyp2m1p1k = q(nx+2+0+0,ny+2-1+1,indvars(2))

d2_rhs_v_dydx_1_0_nxp2p0nyp2m1p1k_nxp2p0m1nyp2m1p1k = q(nx+2+0-1,ny+2-1+1,indvars(2))

d2_rhs_v_dydx_1_0_nxp2p0nyp2m1p1k_nxp2p0m2nyp2m1p1k = q(nx+2+0-2,ny+2-1+1,indvars(2))

d2_rhs_v_dydx_1_0_nxp2p0nyp2m1p1k = 1.5_wp*d2_rhs_v_dydx_1_0_nxp2p0nyp2m1p1k_nxp2p0p0nyp2m1p1k-&
          2.0_wp*d2_rhs_v_dydx_1_0_nxp2p0nyp2m1p1k_nxp2p0m1nyp2m1p1k+&
          0.5_wp*d2_rhs_v_dydx_1_0_nxp2p0nyp2m1p1k_nxp2p0m2nyp2m1p1k

d2_rhs_v_dydx_1_0_nxp2p0nyp2m1p1k = d2_rhs_v_dydx_1_0_nxp2p0nyp2m1p1k*param_float(1)

d2_rhs_v_dydy_1_0_nxp2p0nyp2m1m1k_nxp2p0nyp2m1m1m1k = q(nx+2+0,ny+2-1-1-1,indvars(3))

d2_rhs_v_dydy_1_0_nxp2p0nyp2m1m1k_nxp2p0nyp2m1m1p1k = q(nx+2+0,ny+2-1-1+1,indvars(3))

d2_rhs_v_dydy_1_0_nxp2p0nyp2m1m1k = -&
          0.5_wp*d2_rhs_v_dydy_1_0_nxp2p0nyp2m1m1k_nxp2p0nyp2m1m1m1k+&
          0.5_wp*d2_rhs_v_dydy_1_0_nxp2p0nyp2m1m1k_nxp2p0nyp2m1m1p1k

d2_rhs_v_dydy_1_0_nxp2p0nyp2m1m1k = d2_rhs_v_dydy_1_0_nxp2p0nyp2m1m1k*param_float(2)

d2_rhs_v_dydy_1_0_nxp2p0nyp2m1p1k_nxp2p0nyp2m1p1p0k = q(nx+2+0,ny+2-1+1+0,indvars(3))

d2_rhs_v_dydy_1_0_nxp2p0nyp2m1p1k_nxp2p0nyp2m1p1m1k = q(nx+2+0,ny+2-1+1-1,indvars(3))

d2_rhs_v_dydy_1_0_nxp2p0nyp2m1p1k_nxp2p0nyp2m1p1m2k = q(nx+2+0,ny+2-1+1-2,indvars(3))

d2_rhs_v_dydy_1_0_nxp2p0nyp2m1p1k = 1.5_wp*d2_rhs_v_dydy_1_0_nxp2p0nyp2m1p1k_nxp2p0nyp2m1p1p0k-&
          2.0_wp*d2_rhs_v_dydy_1_0_nxp2p0nyp2m1p1k_nxp2p0nyp2m1p1m1k+&
          0.5_wp*d2_rhs_v_dydy_1_0_nxp2p0nyp2m1p1k_nxp2p0nyp2m1p1m2k

d2_rhs_v_dydy_1_0_nxp2p0nyp2m1p1k = d2_rhs_v_dydy_1_0_nxp2p0nyp2m1p1k*param_float(2)

d1_rhs_v_dy_1_nxp2p0nyp2m1m1k = -((1+&
                    param_float(21 + 5))/(((q(nx+2+0,ny+2-1-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2-1-1,indvars(2))*q(nx+2+0,ny+2-1-1,indvars(2))+&
                    q(nx+2+0,ny+2-1-1,indvars(3))*q(nx+2+0,ny+2-1-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2+0,ny+2-1-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2-1-1,indvars(2))*q(nx+2+0,ny+2-1-1,indvars(2))+&
                    q(nx+2+0,ny+2-1-1,indvars(3))*q(nx+2+0,ny+2-1-1,indvars(3)))))/param_float(4 + 5)**1.5)*(1+&
                    ((q(nx+2+0,ny+2-1-1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+2+0,ny+2-1-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2-1-1,indvars(2))*q(nx+2+0,ny+2-1-1,indvars(2))+&
                    q(nx+2+0,ny+2-1-1,indvars(3))*q(nx+2+0,ny+2-1-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2+0,ny+2-1-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2-1-1,indvars(2))*q(nx+2+0,ny+2-1-1,indvars(2))+&
                    q(nx+2+0,ny+2-1-1,indvars(3))*q(nx+2+0,ny+2-1-1,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+2+0,ny+2-1-1,indvars(1)))**3/((q(nx+2+0,ny+2-1-1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+2+0,ny+2-1-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2-1-1,indvars(2))*q(nx+2+0,ny+2-1-1,indvars(2))+&
                    q(nx+2+0,ny+2-1-1,indvars(3))*q(nx+2+0,ny+2-1-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2+0,ny+2-1-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2-1-1,indvars(2))*q(nx+2+0,ny+2-1-1,indvars(2))+&
                    q(nx+2+0,ny+2-1-1,indvars(3))*q(nx+2+0,ny+2-1-1,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+2+0,ny+2-1-1,indvars(1)))**3+&
                    param_float(13 + 5)**3))*(q(nx+2+0,ny+2-1-1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+2+0,ny+2-1-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2-1-1,indvars(2))*q(nx+2+0,ny+2-1-1,indvars(2))+&
                    q(nx+2+0,ny+2-1-1,indvars(3))*q(nx+2+0,ny+2-1-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2+0,ny+2-1-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2-1-1,indvars(2))*q(nx+2+0,ny+2-1-1,indvars(2))+&
                    q(nx+2+0,ny+2-1-1,indvars(3))*q(nx+2+0,ny+2-1-1,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+2+0,ny+2-1-1,indvars(1))))*param_float(1 + 5)*(2.0_wp*qst(nx+2+0,ny+2-1-1,indvarsst(11))*(d2_rhs_v_dydy_1_0_nxp2p0nyp2m1m1k)-&
                    2.0_wp/3.0_wp*(qst(nx+2+0,ny+2-1-1,indvarsst(10))*(d2_rhs_v_dydx_1_0_nxp2p0nyp2m1m1k)+&
                    qst(nx+2+0,ny+2-1-1,indvarsst(11))*(d2_rhs_v_dydy_1_0_nxp2p0nyp2m1m1k)))

d1_rhs_v_dy_1_nxp2p0nyp2m1p1k = -((1+&
                    param_float(21 + 5))/(((q(nx+2+0,ny+2-1+1,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2-1+1,indvars(2))*q(nx+2+0,ny+2-1+1,indvars(2))+&
                    q(nx+2+0,ny+2-1+1,indvars(3))*q(nx+2+0,ny+2-1+1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2+0,ny+2-1+1,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2-1+1,indvars(2))*q(nx+2+0,ny+2-1+1,indvars(2))+&
                    q(nx+2+0,ny+2-1+1,indvars(3))*q(nx+2+0,ny+2-1+1,indvars(3)))))/param_float(4 + 5)**1.5)*(1+&
                    ((q(nx+2+0,ny+2-1+1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+2+0,ny+2-1+1,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2-1+1,indvars(2))*q(nx+2+0,ny+2-1+1,indvars(2))+&
                    q(nx+2+0,ny+2-1+1,indvars(3))*q(nx+2+0,ny+2-1+1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2+0,ny+2-1+1,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2-1+1,indvars(2))*q(nx+2+0,ny+2-1+1,indvars(2))+&
                    q(nx+2+0,ny+2-1+1,indvars(3))*q(nx+2+0,ny+2-1+1,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+2+0,ny+2-1+1,indvars(1)))**3/((q(nx+2+0,ny+2-1+1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+2+0,ny+2-1+1,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2-1+1,indvars(2))*q(nx+2+0,ny+2-1+1,indvars(2))+&
                    q(nx+2+0,ny+2-1+1,indvars(3))*q(nx+2+0,ny+2-1+1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2+0,ny+2-1+1,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2-1+1,indvars(2))*q(nx+2+0,ny+2-1+1,indvars(2))+&
                    q(nx+2+0,ny+2-1+1,indvars(3))*q(nx+2+0,ny+2-1+1,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+2+0,ny+2-1+1,indvars(1)))**3+&
                    param_float(13 + 5)**3))*(q(nx+2+0,ny+2-1+1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+2+0,ny+2-1+1,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2-1+1,indvars(2))*q(nx+2+0,ny+2-1+1,indvars(2))+&
                    q(nx+2+0,ny+2-1+1,indvars(3))*q(nx+2+0,ny+2-1+1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2+0,ny+2-1+1,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2-1+1,indvars(2))*q(nx+2+0,ny+2-1+1,indvars(2))+&
                    q(nx+2+0,ny+2-1+1,indvars(3))*q(nx+2+0,ny+2-1+1,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+2+0,ny+2-1+1,indvars(1))))*param_float(1 + 5)*(2.0_wp*qst(nx+2+0,ny+2-1+1,indvarsst(11))*(d2_rhs_v_dydy_1_0_nxp2p0nyp2m1p1k)-&
                    2.0_wp/3.0_wp*(qst(nx+2+0,ny+2-1+1,indvarsst(10))*(d2_rhs_v_dydx_1_0_nxp2p0nyp2m1p1k)+&
                    qst(nx+2+0,ny+2-1+1,indvarsst(11))*(d2_rhs_v_dydy_1_0_nxp2p0nyp2m1p1k)))

d1_rhs_v_dy_1_nxp2p0nyp2m1k = -&
          0.5_wp*d1_rhs_v_dy_1_nxp2p0nyp2m1m1k+&
          0.5_wp*d1_rhs_v_dy_1_nxp2p0nyp2m1p1k

d1_rhs_v_dy_1_nxp2p0nyp2m1k = d1_rhs_v_dy_1_nxp2p0nyp2m1k*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer 0 1 None d(rho v)/dt ***********
!                                                           
!***********************************************************


rhs(nx+2+0,ny+2-1,indvars(3)) =   -  ( q(nx+2+0,ny+2-1,indvars(3))*(1.0_wp/(param_float(23 + 5)*(param_float(3 + 5))*q(nx+2+0,ny+2-1,indvars(1))*((q(nx+2+0,ny+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2-1,indvars(2))*q(nx+2+0,ny+2-1,indvars(2))+&
                    q(nx+2+0,ny+2-1,indvars(3))*q(nx+2+0,ny+2-1,indvars(3)))))/q(nx+2+0,ny+2-1,indvars(1)))**2*((q(nx+2+0,ny+2-1,indvars(2))*(-&
                    (param_float(23 + 5)*(param_float(3 + 5))*q(nx+2+0,ny+2-1,indvars(1))*((q(nx+2+0,ny+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2-1,indvars(2))*q(nx+2+0,ny+2-1,indvars(2))+&
                    q(nx+2+0,ny+2-1,indvars(3))*q(nx+2+0,ny+2-1,indvars(3)))))/q(nx+2+0,ny+2-1,indvars(1)))**2*d1_rhs_rho_dx_1_nxp2p0nyp2m1k+&
                    1.0_wp*d1_rhs_rho_dx_2_nxp2p0nyp2m1k))*qst(nx+2+0,ny+2-1,indvarsst(10))+&
                    0.5_wp*((((param_float(23 + 5)*(param_float(3 + 5))*q(nx+2+0,ny+2-1,indvars(1))*((q(nx+2+0,ny+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2-1,indvars(2))*q(nx+2+0,ny+2-1,indvars(2))+&
                    q(nx+2+0,ny+2-1,indvars(3))*q(nx+2+0,ny+2-1,indvars(3)))))/q(nx+2+0,ny+2-1,indvars(1)))+&
                    q(nx+2+0,ny+2-1,indvars(2)))*((param_float(23 + 5)*(param_float(3 + 5))*q(nx+2+0,ny+2-1,indvars(1))*((q(nx+2+0,ny+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2-1,indvars(2))*q(nx+2+0,ny+2-1,indvars(2))+&
                    q(nx+2+0,ny+2-1,indvars(3))*q(nx+2+0,ny+2-1,indvars(3)))))/q(nx+2+0,ny+2-1,indvars(1)))*q(nx+2+0,ny+2-1,indvars(1))*d1_rhs_rho_dx_3_nxp2p0nyp2m1k+&
                    1.0_wp*d1_rhs_rho_dx_2_nxp2p0nyp2m1k))*qst(nx+2+0,ny+2-1,indvarsst(10))+&
                    -&
                    (((param_float(23 + 5)*(param_float(3 + 5))*q(nx+2+0,ny+2-1,indvars(1))*((q(nx+2+0,ny+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2-1,indvars(2))*q(nx+2+0,ny+2-1,indvars(2))+&
                    q(nx+2+0,ny+2-1,indvars(3))*q(nx+2+0,ny+2-1,indvars(3)))))/q(nx+2+0,ny+2-1,indvars(1)))+&
                    q(nx+2+0,ny+2-1,indvars(2)))*((param_float(23 + 5)*(param_float(3 + 5))*q(nx+2+0,ny+2-1,indvars(1))*((q(nx+2+0,ny+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2-1,indvars(2))*q(nx+2+0,ny+2-1,indvars(2))+&
                    q(nx+2+0,ny+2-1,indvars(3))*q(nx+2+0,ny+2-1,indvars(3)))))/q(nx+2+0,ny+2-1,indvars(1)))*q(nx+2+0,ny+2-1,indvars(1))*d1_rhs_rho_dx_3_nxp2p0nyp2m1k+&
                    1.0_wp*d1_rhs_rho_dx_2_nxp2p0nyp2m1k))*qst(nx+2+0,ny+2-1,indvarsst(10)))))+&
                    q(nx+2+0,ny+2-1,indvars(1))*((1.0_wp*q(nx+2+0,ny+2-1,indvars(2))*d1_rhs_v_dx_6_nxp2p0nyp2m1k)*qst(nx+2+0,ny+2-1,indvarsst(10)))+&
                    d1_rhs_v_dy_0_nxp2p0nyp2m1k+&
                    qst(nx+2+0,ny+2-1,indvarsst(10))*(d1_rhs_v_dx_7_nxp2p0nyp2m1k)+&
                    qst(nx+2+0,ny+2-1,indvarsst(11))*(d1_rhs_v_dy_1_nxp2p0nyp2m1k) ) 



!***********************************************************
!                                                           
! building source terms in RHS for layer 0 1 None d(rho et)/dt 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! 0.5_wp*(u**2+v**2)*(1.0_wp/c**2*((u*(-c**2*[rho]_1x+1.0_wp*[p]_1x))*deltaxI+0.5_wp*(((c+u)*(c*rho*[u]_1x+1.0_wp*[p]_1x))*deltaxI+-((c+u)*(c*rho*[u]_1x+1.0_wp*[p]_1x))*deltaxI)))+1.0_wp/gamma_m1*(0.5_wp*(((c+u)*(c*rho*[u]_1x+1.0_wp*[p]_1x))*deltaxI+-((c+u)*(c*rho*[u]_1x+1.0_wp*[p]_1x))*deltaxI))+rho*u*(1.0_wp/(2*rho*c)*(((c+u)*(c*rho*[u]_1x+1.0_wp*[p]_1x))*deltaxI--((c+u)*(c*rho*[u]_1x+1.0_wp*[p]_1x))*deltaxI))+rho*v*((1.0_wp*u*[v]_1x)*deltaxI)+[(rho*et+p)*v]_1y+deltaxI*([-kappa*deltaxI*({T}_1x)-u*(visc_t*(2.0_wp*deltaxI*({u}_1x)-2.0_wp/3.0_wp*(deltaxI*({u}_1x)+deltayI*({v}_1y))))-v*(visc_t*(deltayI*({u}_1y)+deltaxI*({v}_1x)))]_1x)+deltayI*([-kappa*deltayI*({T}_1y)-u*(visc_t*(deltayI*({u}_1y)+deltaxI*({v}_1x)))-v*(visc_t*(2.0_wp*deltayI*({v}_1y)-2.0_wp/3.0_wp*(deltaxI*({u}_1x)+deltayI*({v}_1y))))]_1y)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d2_rhs_et_dxdx_15_0_nxp2p0p0nyp2m1k_nxp2p0p0p0nyp2m1k = ((q(nx+2+0+0+0,ny+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0+0+0,ny+2-1,indvars(2))*q(nx+2+0+0+0,ny+2-1,indvars(2))+&
                    q(nx+2+0+0+0,ny+2-1,indvars(3))*q(nx+2+0+0+0,ny+2-1,indvars(3)))))/param_float(4 + 5)

d2_rhs_et_dxdx_15_0_nxp2p0p0nyp2m1k_nxp2p0p0m1nyp2m1k = ((q(nx+2+0+0-1,ny+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0+0-1,ny+2-1,indvars(2))*q(nx+2+0+0-1,ny+2-1,indvars(2))+&
                    q(nx+2+0+0-1,ny+2-1,indvars(3))*q(nx+2+0+0-1,ny+2-1,indvars(3)))))/param_float(4 + 5)

d2_rhs_et_dxdx_15_0_nxp2p0p0nyp2m1k_nxp2p0p0m2nyp2m1k = ((q(nx+2+0+0-2,ny+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0+0-2,ny+2-1,indvars(2))*q(nx+2+0+0-2,ny+2-1,indvars(2))+&
                    q(nx+2+0+0-2,ny+2-1,indvars(3))*q(nx+2+0+0-2,ny+2-1,indvars(3)))))/param_float(4 + 5)

d2_rhs_et_dxdx_15_0_nxp2p0p0nyp2m1k = 1.5_wp*d2_rhs_et_dxdx_15_0_nxp2p0p0nyp2m1k_nxp2p0p0p0nyp2m1k-&
          2.0_wp*d2_rhs_et_dxdx_15_0_nxp2p0p0nyp2m1k_nxp2p0p0m1nyp2m1k+&
          0.5_wp*d2_rhs_et_dxdx_15_0_nxp2p0p0nyp2m1k_nxp2p0p0m2nyp2m1k

d2_rhs_et_dxdx_15_0_nxp2p0p0nyp2m1k = d2_rhs_et_dxdx_15_0_nxp2p0p0nyp2m1k*param_float(1)

d2_rhs_et_dxdx_15_0_nxp2p0m1nyp2m1k_nxp2p0m1m1nyp2m1k = ((q(nx+2+0-1-1,ny+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0-1-1,ny+2-1,indvars(2))*q(nx+2+0-1-1,ny+2-1,indvars(2))+&
                    q(nx+2+0-1-1,ny+2-1,indvars(3))*q(nx+2+0-1-1,ny+2-1,indvars(3)))))/param_float(4 + 5)

d2_rhs_et_dxdx_15_0_nxp2p0m1nyp2m1k_nxp2p0m1p1nyp2m1k = ((q(nx+2+0-1+1,ny+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0-1+1,ny+2-1,indvars(2))*q(nx+2+0-1+1,ny+2-1,indvars(2))+&
                    q(nx+2+0-1+1,ny+2-1,indvars(3))*q(nx+2+0-1+1,ny+2-1,indvars(3)))))/param_float(4 + 5)

d2_rhs_et_dxdx_15_0_nxp2p0m1nyp2m1k = -&
          0.5_wp*d2_rhs_et_dxdx_15_0_nxp2p0m1nyp2m1k_nxp2p0m1m1nyp2m1k+&
          0.5_wp*d2_rhs_et_dxdx_15_0_nxp2p0m1nyp2m1k_nxp2p0m1p1nyp2m1k

d2_rhs_et_dxdx_15_0_nxp2p0m1nyp2m1k = d2_rhs_et_dxdx_15_0_nxp2p0m1nyp2m1k*param_float(1)

d2_rhs_et_dxdx_15_0_nxp2p0m2nyp2m1k_nxp2p0m2m1nyp2m1k = ((q(nx+2+0-2-1,ny+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0-2-1,ny+2-1,indvars(2))*q(nx+2+0-2-1,ny+2-1,indvars(2))+&
                    q(nx+2+0-2-1,ny+2-1,indvars(3))*q(nx+2+0-2-1,ny+2-1,indvars(3)))))/param_float(4 + 5)

d2_rhs_et_dxdx_15_0_nxp2p0m2nyp2m1k_nxp2p0m2p1nyp2m1k = ((q(nx+2+0-2+1,ny+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0-2+1,ny+2-1,indvars(2))*q(nx+2+0-2+1,ny+2-1,indvars(2))+&
                    q(nx+2+0-2+1,ny+2-1,indvars(3))*q(nx+2+0-2+1,ny+2-1,indvars(3)))))/param_float(4 + 5)

d2_rhs_et_dxdx_15_0_nxp2p0m2nyp2m1k = -&
          0.5_wp*d2_rhs_et_dxdx_15_0_nxp2p0m2nyp2m1k_nxp2p0m2m1nyp2m1k+&
          0.5_wp*d2_rhs_et_dxdx_15_0_nxp2p0m2nyp2m1k_nxp2p0m2p1nyp2m1k

d2_rhs_et_dxdx_15_0_nxp2p0m2nyp2m1k = d2_rhs_et_dxdx_15_0_nxp2p0m2nyp2m1k*param_float(1)

d1_rhs_et_dx_15_nxp2p0p0nyp2m1k = -param_float(2 + 5)*qst(nx+2+0+0,ny+2-1,indvarsst(10))*(d2_rhs_et_dxdx_15_0_nxp2p0p0nyp2m1k)-&
                    q(nx+2+0+0,ny+2-1,indvars(2))*(((1+&
                    param_float(21 + 5))/(((q(nx+2+0+0,ny+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0+0,ny+2-1,indvars(2))*q(nx+2+0+0,ny+2-1,indvars(2))+&
                    q(nx+2+0+0,ny+2-1,indvars(3))*q(nx+2+0+0,ny+2-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2+0+0,ny+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0+0,ny+2-1,indvars(2))*q(nx+2+0+0,ny+2-1,indvars(2))+&
                    q(nx+2+0+0,ny+2-1,indvars(3))*q(nx+2+0+0,ny+2-1,indvars(3)))))/param_float(4 + 5)**1.5)*(1+&
                    ((q(nx+2+0+0,ny+2-1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+2+0+0,ny+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0+0,ny+2-1,indvars(2))*q(nx+2+0+0,ny+2-1,indvars(2))+&
                    q(nx+2+0+0,ny+2-1,indvars(3))*q(nx+2+0+0,ny+2-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2+0+0,ny+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0+0,ny+2-1,indvars(2))*q(nx+2+0+0,ny+2-1,indvars(2))+&
                    q(nx+2+0+0,ny+2-1,indvars(3))*q(nx+2+0+0,ny+2-1,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+2+0+0,ny+2-1,indvars(1)))**3/((q(nx+2+0+0,ny+2-1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+2+0+0,ny+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0+0,ny+2-1,indvars(2))*q(nx+2+0+0,ny+2-1,indvars(2))+&
                    q(nx+2+0+0,ny+2-1,indvars(3))*q(nx+2+0+0,ny+2-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2+0+0,ny+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0+0,ny+2-1,indvars(2))*q(nx+2+0+0,ny+2-1,indvars(2))+&
                    q(nx+2+0+0,ny+2-1,indvars(3))*q(nx+2+0+0,ny+2-1,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+2+0+0,ny+2-1,indvars(1)))**3+&
                    param_float(13 + 5)**3))*(q(nx+2+0+0,ny+2-1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+2+0+0,ny+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0+0,ny+2-1,indvars(2))*q(nx+2+0+0,ny+2-1,indvars(2))+&
                    q(nx+2+0+0,ny+2-1,indvars(3))*q(nx+2+0+0,ny+2-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2+0+0,ny+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0+0,ny+2-1,indvars(2))*q(nx+2+0+0,ny+2-1,indvars(2))+&
                    q(nx+2+0+0,ny+2-1,indvars(3))*q(nx+2+0+0,ny+2-1,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+2+0+0,ny+2-1,indvars(1))))*param_float(1 + 5)*(2.0_wp*qst(nx+2+0+0,ny+2-1,indvarsst(10))*(d2_rhs_u_dxdx_10_0_nxp2p0p0nyp2m1k)-&
                    2.0_wp/3.0_wp*(qst(nx+2+0+0,ny+2-1,indvarsst(10))*(d2_rhs_u_dxdx_10_0_nxp2p0p0nyp2m1k)+&
                    qst(nx+2+0+0,ny+2-1,indvarsst(11))*(d2_rhs_u_dxdy_10_0_nxp2p0p0nyp2m1k))))-&
                    q(nx+2+0+0,ny+2-1,indvars(3))*(((1+&
                    param_float(21 + 5))/(((q(nx+2+0+0,ny+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0+0,ny+2-1,indvars(2))*q(nx+2+0+0,ny+2-1,indvars(2))+&
                    q(nx+2+0+0,ny+2-1,indvars(3))*q(nx+2+0+0,ny+2-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2+0+0,ny+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0+0,ny+2-1,indvars(2))*q(nx+2+0+0,ny+2-1,indvars(2))+&
                    q(nx+2+0+0,ny+2-1,indvars(3))*q(nx+2+0+0,ny+2-1,indvars(3)))))/param_float(4 + 5)**1.5)*(1+&
                    ((q(nx+2+0+0,ny+2-1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+2+0+0,ny+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0+0,ny+2-1,indvars(2))*q(nx+2+0+0,ny+2-1,indvars(2))+&
                    q(nx+2+0+0,ny+2-1,indvars(3))*q(nx+2+0+0,ny+2-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2+0+0,ny+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0+0,ny+2-1,indvars(2))*q(nx+2+0+0,ny+2-1,indvars(2))+&
                    q(nx+2+0+0,ny+2-1,indvars(3))*q(nx+2+0+0,ny+2-1,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+2+0+0,ny+2-1,indvars(1)))**3/((q(nx+2+0+0,ny+2-1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+2+0+0,ny+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0+0,ny+2-1,indvars(2))*q(nx+2+0+0,ny+2-1,indvars(2))+&
                    q(nx+2+0+0,ny+2-1,indvars(3))*q(nx+2+0+0,ny+2-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2+0+0,ny+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0+0,ny+2-1,indvars(2))*q(nx+2+0+0,ny+2-1,indvars(2))+&
                    q(nx+2+0+0,ny+2-1,indvars(3))*q(nx+2+0+0,ny+2-1,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+2+0+0,ny+2-1,indvars(1)))**3+&
                    param_float(13 + 5)**3))*(q(nx+2+0+0,ny+2-1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+2+0+0,ny+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0+0,ny+2-1,indvars(2))*q(nx+2+0+0,ny+2-1,indvars(2))+&
                    q(nx+2+0+0,ny+2-1,indvars(3))*q(nx+2+0+0,ny+2-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2+0+0,ny+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0+0,ny+2-1,indvars(2))*q(nx+2+0+0,ny+2-1,indvars(2))+&
                    q(nx+2+0+0,ny+2-1,indvars(3))*q(nx+2+0+0,ny+2-1,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+2+0+0,ny+2-1,indvars(1))))*param_float(1 + 5)*(qst(nx+2+0+0,ny+2-1,indvarsst(11))*(d2_rhs_v_dxdy_7_0_nxp2p0p0nyp2m1k)+&
                    qst(nx+2+0+0,ny+2-1,indvarsst(10))*(d2_rhs_v_dxdx_7_0_nxp2p0p0nyp2m1k)))

d1_rhs_et_dx_15_nxp2p0m1nyp2m1k = -param_float(2 + 5)*qst(nx+2+0-1,ny+2-1,indvarsst(10))*(d2_rhs_et_dxdx_15_0_nxp2p0m1nyp2m1k)-&
                    q(nx+2+0-1,ny+2-1,indvars(2))*(((1+&
                    param_float(21 + 5))/(((q(nx+2+0-1,ny+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0-1,ny+2-1,indvars(2))*q(nx+2+0-1,ny+2-1,indvars(2))+&
                    q(nx+2+0-1,ny+2-1,indvars(3))*q(nx+2+0-1,ny+2-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2+0-1,ny+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0-1,ny+2-1,indvars(2))*q(nx+2+0-1,ny+2-1,indvars(2))+&
                    q(nx+2+0-1,ny+2-1,indvars(3))*q(nx+2+0-1,ny+2-1,indvars(3)))))/param_float(4 + 5)**1.5)*(1+&
                    ((q(nx+2+0-1,ny+2-1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+2+0-1,ny+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0-1,ny+2-1,indvars(2))*q(nx+2+0-1,ny+2-1,indvars(2))+&
                    q(nx+2+0-1,ny+2-1,indvars(3))*q(nx+2+0-1,ny+2-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2+0-1,ny+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0-1,ny+2-1,indvars(2))*q(nx+2+0-1,ny+2-1,indvars(2))+&
                    q(nx+2+0-1,ny+2-1,indvars(3))*q(nx+2+0-1,ny+2-1,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+2+0-1,ny+2-1,indvars(1)))**3/((q(nx+2+0-1,ny+2-1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+2+0-1,ny+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0-1,ny+2-1,indvars(2))*q(nx+2+0-1,ny+2-1,indvars(2))+&
                    q(nx+2+0-1,ny+2-1,indvars(3))*q(nx+2+0-1,ny+2-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2+0-1,ny+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0-1,ny+2-1,indvars(2))*q(nx+2+0-1,ny+2-1,indvars(2))+&
                    q(nx+2+0-1,ny+2-1,indvars(3))*q(nx+2+0-1,ny+2-1,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+2+0-1,ny+2-1,indvars(1)))**3+&
                    param_float(13 + 5)**3))*(q(nx+2+0-1,ny+2-1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+2+0-1,ny+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0-1,ny+2-1,indvars(2))*q(nx+2+0-1,ny+2-1,indvars(2))+&
                    q(nx+2+0-1,ny+2-1,indvars(3))*q(nx+2+0-1,ny+2-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2+0-1,ny+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0-1,ny+2-1,indvars(2))*q(nx+2+0-1,ny+2-1,indvars(2))+&
                    q(nx+2+0-1,ny+2-1,indvars(3))*q(nx+2+0-1,ny+2-1,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+2+0-1,ny+2-1,indvars(1))))*param_float(1 + 5)*(2.0_wp*qst(nx+2+0-1,ny+2-1,indvarsst(10))*(d2_rhs_u_dxdx_10_0_nxp2p0m1nyp2m1k)-&
                    2.0_wp/3.0_wp*(qst(nx+2+0-1,ny+2-1,indvarsst(10))*(d2_rhs_u_dxdx_10_0_nxp2p0m1nyp2m1k)+&
                    qst(nx+2+0-1,ny+2-1,indvarsst(11))*(d2_rhs_u_dxdy_10_0_nxp2p0m1nyp2m1k))))-&
                    q(nx+2+0-1,ny+2-1,indvars(3))*(((1+&
                    param_float(21 + 5))/(((q(nx+2+0-1,ny+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0-1,ny+2-1,indvars(2))*q(nx+2+0-1,ny+2-1,indvars(2))+&
                    q(nx+2+0-1,ny+2-1,indvars(3))*q(nx+2+0-1,ny+2-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2+0-1,ny+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0-1,ny+2-1,indvars(2))*q(nx+2+0-1,ny+2-1,indvars(2))+&
                    q(nx+2+0-1,ny+2-1,indvars(3))*q(nx+2+0-1,ny+2-1,indvars(3)))))/param_float(4 + 5)**1.5)*(1+&
                    ((q(nx+2+0-1,ny+2-1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+2+0-1,ny+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0-1,ny+2-1,indvars(2))*q(nx+2+0-1,ny+2-1,indvars(2))+&
                    q(nx+2+0-1,ny+2-1,indvars(3))*q(nx+2+0-1,ny+2-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2+0-1,ny+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0-1,ny+2-1,indvars(2))*q(nx+2+0-1,ny+2-1,indvars(2))+&
                    q(nx+2+0-1,ny+2-1,indvars(3))*q(nx+2+0-1,ny+2-1,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+2+0-1,ny+2-1,indvars(1)))**3/((q(nx+2+0-1,ny+2-1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+2+0-1,ny+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0-1,ny+2-1,indvars(2))*q(nx+2+0-1,ny+2-1,indvars(2))+&
                    q(nx+2+0-1,ny+2-1,indvars(3))*q(nx+2+0-1,ny+2-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2+0-1,ny+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0-1,ny+2-1,indvars(2))*q(nx+2+0-1,ny+2-1,indvars(2))+&
                    q(nx+2+0-1,ny+2-1,indvars(3))*q(nx+2+0-1,ny+2-1,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+2+0-1,ny+2-1,indvars(1)))**3+&
                    param_float(13 + 5)**3))*(q(nx+2+0-1,ny+2-1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+2+0-1,ny+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0-1,ny+2-1,indvars(2))*q(nx+2+0-1,ny+2-1,indvars(2))+&
                    q(nx+2+0-1,ny+2-1,indvars(3))*q(nx+2+0-1,ny+2-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2+0-1,ny+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0-1,ny+2-1,indvars(2))*q(nx+2+0-1,ny+2-1,indvars(2))+&
                    q(nx+2+0-1,ny+2-1,indvars(3))*q(nx+2+0-1,ny+2-1,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+2+0-1,ny+2-1,indvars(1))))*param_float(1 + 5)*(qst(nx+2+0-1,ny+2-1,indvarsst(11))*(d2_rhs_v_dxdy_7_0_nxp2p0m1nyp2m1k)+&
                    qst(nx+2+0-1,ny+2-1,indvarsst(10))*(d2_rhs_v_dxdx_7_0_nxp2p0m1nyp2m1k)))

d1_rhs_et_dx_15_nxp2p0m2nyp2m1k = -param_float(2 + 5)*qst(nx+2+0-2,ny+2-1,indvarsst(10))*(d2_rhs_et_dxdx_15_0_nxp2p0m2nyp2m1k)-&
                    q(nx+2+0-2,ny+2-1,indvars(2))*(((1+&
                    param_float(21 + 5))/(((q(nx+2+0-2,ny+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0-2,ny+2-1,indvars(2))*q(nx+2+0-2,ny+2-1,indvars(2))+&
                    q(nx+2+0-2,ny+2-1,indvars(3))*q(nx+2+0-2,ny+2-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2+0-2,ny+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0-2,ny+2-1,indvars(2))*q(nx+2+0-2,ny+2-1,indvars(2))+&
                    q(nx+2+0-2,ny+2-1,indvars(3))*q(nx+2+0-2,ny+2-1,indvars(3)))))/param_float(4 + 5)**1.5)*(1+&
                    ((q(nx+2+0-2,ny+2-1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+2+0-2,ny+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0-2,ny+2-1,indvars(2))*q(nx+2+0-2,ny+2-1,indvars(2))+&
                    q(nx+2+0-2,ny+2-1,indvars(3))*q(nx+2+0-2,ny+2-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2+0-2,ny+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0-2,ny+2-1,indvars(2))*q(nx+2+0-2,ny+2-1,indvars(2))+&
                    q(nx+2+0-2,ny+2-1,indvars(3))*q(nx+2+0-2,ny+2-1,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+2+0-2,ny+2-1,indvars(1)))**3/((q(nx+2+0-2,ny+2-1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+2+0-2,ny+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0-2,ny+2-1,indvars(2))*q(nx+2+0-2,ny+2-1,indvars(2))+&
                    q(nx+2+0-2,ny+2-1,indvars(3))*q(nx+2+0-2,ny+2-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2+0-2,ny+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0-2,ny+2-1,indvars(2))*q(nx+2+0-2,ny+2-1,indvars(2))+&
                    q(nx+2+0-2,ny+2-1,indvars(3))*q(nx+2+0-2,ny+2-1,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+2+0-2,ny+2-1,indvars(1)))**3+&
                    param_float(13 + 5)**3))*(q(nx+2+0-2,ny+2-1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+2+0-2,ny+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0-2,ny+2-1,indvars(2))*q(nx+2+0-2,ny+2-1,indvars(2))+&
                    q(nx+2+0-2,ny+2-1,indvars(3))*q(nx+2+0-2,ny+2-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2+0-2,ny+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0-2,ny+2-1,indvars(2))*q(nx+2+0-2,ny+2-1,indvars(2))+&
                    q(nx+2+0-2,ny+2-1,indvars(3))*q(nx+2+0-2,ny+2-1,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+2+0-2,ny+2-1,indvars(1))))*param_float(1 + 5)*(2.0_wp*qst(nx+2+0-2,ny+2-1,indvarsst(10))*(d2_rhs_u_dxdx_10_0_nxp2p0m2nyp2m1k)-&
                    2.0_wp/3.0_wp*(qst(nx+2+0-2,ny+2-1,indvarsst(10))*(d2_rhs_u_dxdx_10_0_nxp2p0m2nyp2m1k)+&
                    qst(nx+2+0-2,ny+2-1,indvarsst(11))*(d2_rhs_u_dxdy_10_0_nxp2p0m2nyp2m1k))))-&
                    q(nx+2+0-2,ny+2-1,indvars(3))*(((1+&
                    param_float(21 + 5))/(((q(nx+2+0-2,ny+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0-2,ny+2-1,indvars(2))*q(nx+2+0-2,ny+2-1,indvars(2))+&
                    q(nx+2+0-2,ny+2-1,indvars(3))*q(nx+2+0-2,ny+2-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2+0-2,ny+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0-2,ny+2-1,indvars(2))*q(nx+2+0-2,ny+2-1,indvars(2))+&
                    q(nx+2+0-2,ny+2-1,indvars(3))*q(nx+2+0-2,ny+2-1,indvars(3)))))/param_float(4 + 5)**1.5)*(1+&
                    ((q(nx+2+0-2,ny+2-1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+2+0-2,ny+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0-2,ny+2-1,indvars(2))*q(nx+2+0-2,ny+2-1,indvars(2))+&
                    q(nx+2+0-2,ny+2-1,indvars(3))*q(nx+2+0-2,ny+2-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2+0-2,ny+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0-2,ny+2-1,indvars(2))*q(nx+2+0-2,ny+2-1,indvars(2))+&
                    q(nx+2+0-2,ny+2-1,indvars(3))*q(nx+2+0-2,ny+2-1,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+2+0-2,ny+2-1,indvars(1)))**3/((q(nx+2+0-2,ny+2-1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+2+0-2,ny+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0-2,ny+2-1,indvars(2))*q(nx+2+0-2,ny+2-1,indvars(2))+&
                    q(nx+2+0-2,ny+2-1,indvars(3))*q(nx+2+0-2,ny+2-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2+0-2,ny+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0-2,ny+2-1,indvars(2))*q(nx+2+0-2,ny+2-1,indvars(2))+&
                    q(nx+2+0-2,ny+2-1,indvars(3))*q(nx+2+0-2,ny+2-1,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+2+0-2,ny+2-1,indvars(1)))**3+&
                    param_float(13 + 5)**3))*(q(nx+2+0-2,ny+2-1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+2+0-2,ny+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0-2,ny+2-1,indvars(2))*q(nx+2+0-2,ny+2-1,indvars(2))+&
                    q(nx+2+0-2,ny+2-1,indvars(3))*q(nx+2+0-2,ny+2-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2+0-2,ny+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0-2,ny+2-1,indvars(2))*q(nx+2+0-2,ny+2-1,indvars(2))+&
                    q(nx+2+0-2,ny+2-1,indvars(3))*q(nx+2+0-2,ny+2-1,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+2+0-2,ny+2-1,indvars(1))))*param_float(1 + 5)*(qst(nx+2+0-2,ny+2-1,indvarsst(11))*(d2_rhs_v_dxdy_7_0_nxp2p0m2nyp2m1k)+&
                    qst(nx+2+0-2,ny+2-1,indvarsst(10))*(d2_rhs_v_dxdx_7_0_nxp2p0m2nyp2m1k)))

d1_rhs_et_dx_15_nxp2p0nyp2m1k = 1.5_wp*d1_rhs_et_dx_15_nxp2p0p0nyp2m1k-&
          2.0_wp*d1_rhs_et_dx_15_nxp2p0m1nyp2m1k+&
          0.5_wp*d1_rhs_et_dx_15_nxp2p0m2nyp2m1k

d1_rhs_et_dx_15_nxp2p0nyp2m1k = d1_rhs_et_dx_15_nxp2p0nyp2m1k*param_float(1)

d1_rhs_et_dy_0_nxp2p0nyp2m1m1k = (q(nx+2+0,ny+2-1-1,indvars(1))*q(nx+2+0,ny+2-1-1,indvars(4))+&
                    (param_float(3 + 5))*q(nx+2+0,ny+2-1-1,indvars(1))*((q(nx+2+0,ny+2-1-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2-1-1,indvars(2))*q(nx+2+0,ny+2-1-1,indvars(2))+&
                    q(nx+2+0,ny+2-1-1,indvars(3))*q(nx+2+0,ny+2-1-1,indvars(3))))))*q(nx+2+0,ny+2-1-1,indvars(3))

d1_rhs_et_dy_0_nxp2p0nyp2m1p1k = (q(nx+2+0,ny+2-1+1,indvars(1))*q(nx+2+0,ny+2-1+1,indvars(4))+&
                    (param_float(3 + 5))*q(nx+2+0,ny+2-1+1,indvars(1))*((q(nx+2+0,ny+2-1+1,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2-1+1,indvars(2))*q(nx+2+0,ny+2-1+1,indvars(2))+&
                    q(nx+2+0,ny+2-1+1,indvars(3))*q(nx+2+0,ny+2-1+1,indvars(3))))))*q(nx+2+0,ny+2-1+1,indvars(3))

d1_rhs_et_dy_0_nxp2p0nyp2m1k = -&
          0.5_wp*d1_rhs_et_dy_0_nxp2p0nyp2m1m1k+&
          0.5_wp*d1_rhs_et_dy_0_nxp2p0nyp2m1p1k

d1_rhs_et_dy_0_nxp2p0nyp2m1k = d1_rhs_et_dy_0_nxp2p0nyp2m1k*param_float(2)

d2_rhs_et_dydy_1_0_nxp2p0nyp2m1m1k_nxp2p0nyp2m1m1m1k = ((q(nx+2+0,ny+2-1-1-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2-1-1-1,indvars(2))*q(nx+2+0,ny+2-1-1-1,indvars(2))+&
                    q(nx+2+0,ny+2-1-1-1,indvars(3))*q(nx+2+0,ny+2-1-1-1,indvars(3)))))/param_float(4 + 5)

d2_rhs_et_dydy_1_0_nxp2p0nyp2m1m1k_nxp2p0nyp2m1m1p1k = ((q(nx+2+0,ny+2-1-1+1,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2-1-1+1,indvars(2))*q(nx+2+0,ny+2-1-1+1,indvars(2))+&
                    q(nx+2+0,ny+2-1-1+1,indvars(3))*q(nx+2+0,ny+2-1-1+1,indvars(3)))))/param_float(4 + 5)

d2_rhs_et_dydy_1_0_nxp2p0nyp2m1m1k = -&
          0.5_wp*d2_rhs_et_dydy_1_0_nxp2p0nyp2m1m1k_nxp2p0nyp2m1m1m1k+&
          0.5_wp*d2_rhs_et_dydy_1_0_nxp2p0nyp2m1m1k_nxp2p0nyp2m1m1p1k

d2_rhs_et_dydy_1_0_nxp2p0nyp2m1m1k = d2_rhs_et_dydy_1_0_nxp2p0nyp2m1m1k*param_float(2)

d2_rhs_et_dydy_1_0_nxp2p0nyp2m1p1k_nxp2p0nyp2m1p1p0k = ((q(nx+2+0,ny+2-1+1+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2-1+1+0,indvars(2))*q(nx+2+0,ny+2-1+1+0,indvars(2))+&
                    q(nx+2+0,ny+2-1+1+0,indvars(3))*q(nx+2+0,ny+2-1+1+0,indvars(3)))))/param_float(4 + 5)

d2_rhs_et_dydy_1_0_nxp2p0nyp2m1p1k_nxp2p0nyp2m1p1m1k = ((q(nx+2+0,ny+2-1+1-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2-1+1-1,indvars(2))*q(nx+2+0,ny+2-1+1-1,indvars(2))+&
                    q(nx+2+0,ny+2-1+1-1,indvars(3))*q(nx+2+0,ny+2-1+1-1,indvars(3)))))/param_float(4 + 5)

d2_rhs_et_dydy_1_0_nxp2p0nyp2m1p1k_nxp2p0nyp2m1p1m2k = ((q(nx+2+0,ny+2-1+1-2,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2-1+1-2,indvars(2))*q(nx+2+0,ny+2-1+1-2,indvars(2))+&
                    q(nx+2+0,ny+2-1+1-2,indvars(3))*q(nx+2+0,ny+2-1+1-2,indvars(3)))))/param_float(4 + 5)

d2_rhs_et_dydy_1_0_nxp2p0nyp2m1p1k = 1.5_wp*d2_rhs_et_dydy_1_0_nxp2p0nyp2m1p1k_nxp2p0nyp2m1p1p0k-&
          2.0_wp*d2_rhs_et_dydy_1_0_nxp2p0nyp2m1p1k_nxp2p0nyp2m1p1m1k+&
          0.5_wp*d2_rhs_et_dydy_1_0_nxp2p0nyp2m1p1k_nxp2p0nyp2m1p1m2k

d2_rhs_et_dydy_1_0_nxp2p0nyp2m1p1k = d2_rhs_et_dydy_1_0_nxp2p0nyp2m1p1k*param_float(2)

d1_rhs_et_dy_1_nxp2p0nyp2m1m1k = -param_float(2 + 5)*qst(nx+2+0,ny+2-1-1,indvarsst(11))*(d2_rhs_et_dydy_1_0_nxp2p0nyp2m1m1k)-&
                    q(nx+2+0,ny+2-1-1,indvars(2))*(((1+&
                    param_float(21 + 5))/(((q(nx+2+0,ny+2-1-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2-1-1,indvars(2))*q(nx+2+0,ny+2-1-1,indvars(2))+&
                    q(nx+2+0,ny+2-1-1,indvars(3))*q(nx+2+0,ny+2-1-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2+0,ny+2-1-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2-1-1,indvars(2))*q(nx+2+0,ny+2-1-1,indvars(2))+&
                    q(nx+2+0,ny+2-1-1,indvars(3))*q(nx+2+0,ny+2-1-1,indvars(3)))))/param_float(4 + 5)**1.5)*(1+&
                    ((q(nx+2+0,ny+2-1-1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+2+0,ny+2-1-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2-1-1,indvars(2))*q(nx+2+0,ny+2-1-1,indvars(2))+&
                    q(nx+2+0,ny+2-1-1,indvars(3))*q(nx+2+0,ny+2-1-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2+0,ny+2-1-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2-1-1,indvars(2))*q(nx+2+0,ny+2-1-1,indvars(2))+&
                    q(nx+2+0,ny+2-1-1,indvars(3))*q(nx+2+0,ny+2-1-1,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+2+0,ny+2-1-1,indvars(1)))**3/((q(nx+2+0,ny+2-1-1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+2+0,ny+2-1-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2-1-1,indvars(2))*q(nx+2+0,ny+2-1-1,indvars(2))+&
                    q(nx+2+0,ny+2-1-1,indvars(3))*q(nx+2+0,ny+2-1-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2+0,ny+2-1-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2-1-1,indvars(2))*q(nx+2+0,ny+2-1-1,indvars(2))+&
                    q(nx+2+0,ny+2-1-1,indvars(3))*q(nx+2+0,ny+2-1-1,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+2+0,ny+2-1-1,indvars(1)))**3+&
                    param_float(13 + 5)**3))*(q(nx+2+0,ny+2-1-1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+2+0,ny+2-1-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2-1-1,indvars(2))*q(nx+2+0,ny+2-1-1,indvars(2))+&
                    q(nx+2+0,ny+2-1-1,indvars(3))*q(nx+2+0,ny+2-1-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2+0,ny+2-1-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2-1-1,indvars(2))*q(nx+2+0,ny+2-1-1,indvars(2))+&
                    q(nx+2+0,ny+2-1-1,indvars(3))*q(nx+2+0,ny+2-1-1,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+2+0,ny+2-1-1,indvars(1))))*param_float(1 + 5)*(qst(nx+2+0,ny+2-1-1,indvarsst(11))*(d2_rhs_u_dydy_1_0_nxp2p0nyp2m1m1k)+&
                    qst(nx+2+0,ny+2-1-1,indvarsst(10))*(d2_rhs_u_dydx_1_0_nxp2p0nyp2m1m1k)))-&
                    q(nx+2+0,ny+2-1-1,indvars(3))*(((1+&
                    param_float(21 + 5))/(((q(nx+2+0,ny+2-1-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2-1-1,indvars(2))*q(nx+2+0,ny+2-1-1,indvars(2))+&
                    q(nx+2+0,ny+2-1-1,indvars(3))*q(nx+2+0,ny+2-1-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2+0,ny+2-1-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2-1-1,indvars(2))*q(nx+2+0,ny+2-1-1,indvars(2))+&
                    q(nx+2+0,ny+2-1-1,indvars(3))*q(nx+2+0,ny+2-1-1,indvars(3)))))/param_float(4 + 5)**1.5)*(1+&
                    ((q(nx+2+0,ny+2-1-1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+2+0,ny+2-1-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2-1-1,indvars(2))*q(nx+2+0,ny+2-1-1,indvars(2))+&
                    q(nx+2+0,ny+2-1-1,indvars(3))*q(nx+2+0,ny+2-1-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2+0,ny+2-1-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2-1-1,indvars(2))*q(nx+2+0,ny+2-1-1,indvars(2))+&
                    q(nx+2+0,ny+2-1-1,indvars(3))*q(nx+2+0,ny+2-1-1,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+2+0,ny+2-1-1,indvars(1)))**3/((q(nx+2+0,ny+2-1-1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+2+0,ny+2-1-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2-1-1,indvars(2))*q(nx+2+0,ny+2-1-1,indvars(2))+&
                    q(nx+2+0,ny+2-1-1,indvars(3))*q(nx+2+0,ny+2-1-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2+0,ny+2-1-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2-1-1,indvars(2))*q(nx+2+0,ny+2-1-1,indvars(2))+&
                    q(nx+2+0,ny+2-1-1,indvars(3))*q(nx+2+0,ny+2-1-1,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+2+0,ny+2-1-1,indvars(1)))**3+&
                    param_float(13 + 5)**3))*(q(nx+2+0,ny+2-1-1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+2+0,ny+2-1-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2-1-1,indvars(2))*q(nx+2+0,ny+2-1-1,indvars(2))+&
                    q(nx+2+0,ny+2-1-1,indvars(3))*q(nx+2+0,ny+2-1-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2+0,ny+2-1-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2-1-1,indvars(2))*q(nx+2+0,ny+2-1-1,indvars(2))+&
                    q(nx+2+0,ny+2-1-1,indvars(3))*q(nx+2+0,ny+2-1-1,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+2+0,ny+2-1-1,indvars(1))))*param_float(1 + 5)*(2.0_wp*qst(nx+2+0,ny+2-1-1,indvarsst(11))*(d2_rhs_v_dydy_1_0_nxp2p0nyp2m1m1k)-&
                    2.0_wp/3.0_wp*(qst(nx+2+0,ny+2-1-1,indvarsst(10))*(d2_rhs_v_dydx_1_0_nxp2p0nyp2m1m1k)+&
                    qst(nx+2+0,ny+2-1-1,indvarsst(11))*(d2_rhs_v_dydy_1_0_nxp2p0nyp2m1m1k))))

d1_rhs_et_dy_1_nxp2p0nyp2m1p1k = -param_float(2 + 5)*qst(nx+2+0,ny+2-1+1,indvarsst(11))*(d2_rhs_et_dydy_1_0_nxp2p0nyp2m1p1k)-&
                    q(nx+2+0,ny+2-1+1,indvars(2))*(((1+&
                    param_float(21 + 5))/(((q(nx+2+0,ny+2-1+1,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2-1+1,indvars(2))*q(nx+2+0,ny+2-1+1,indvars(2))+&
                    q(nx+2+0,ny+2-1+1,indvars(3))*q(nx+2+0,ny+2-1+1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2+0,ny+2-1+1,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2-1+1,indvars(2))*q(nx+2+0,ny+2-1+1,indvars(2))+&
                    q(nx+2+0,ny+2-1+1,indvars(3))*q(nx+2+0,ny+2-1+1,indvars(3)))))/param_float(4 + 5)**1.5)*(1+&
                    ((q(nx+2+0,ny+2-1+1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+2+0,ny+2-1+1,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2-1+1,indvars(2))*q(nx+2+0,ny+2-1+1,indvars(2))+&
                    q(nx+2+0,ny+2-1+1,indvars(3))*q(nx+2+0,ny+2-1+1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2+0,ny+2-1+1,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2-1+1,indvars(2))*q(nx+2+0,ny+2-1+1,indvars(2))+&
                    q(nx+2+0,ny+2-1+1,indvars(3))*q(nx+2+0,ny+2-1+1,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+2+0,ny+2-1+1,indvars(1)))**3/((q(nx+2+0,ny+2-1+1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+2+0,ny+2-1+1,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2-1+1,indvars(2))*q(nx+2+0,ny+2-1+1,indvars(2))+&
                    q(nx+2+0,ny+2-1+1,indvars(3))*q(nx+2+0,ny+2-1+1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2+0,ny+2-1+1,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2-1+1,indvars(2))*q(nx+2+0,ny+2-1+1,indvars(2))+&
                    q(nx+2+0,ny+2-1+1,indvars(3))*q(nx+2+0,ny+2-1+1,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+2+0,ny+2-1+1,indvars(1)))**3+&
                    param_float(13 + 5)**3))*(q(nx+2+0,ny+2-1+1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+2+0,ny+2-1+1,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2-1+1,indvars(2))*q(nx+2+0,ny+2-1+1,indvars(2))+&
                    q(nx+2+0,ny+2-1+1,indvars(3))*q(nx+2+0,ny+2-1+1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2+0,ny+2-1+1,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2-1+1,indvars(2))*q(nx+2+0,ny+2-1+1,indvars(2))+&
                    q(nx+2+0,ny+2-1+1,indvars(3))*q(nx+2+0,ny+2-1+1,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+2+0,ny+2-1+1,indvars(1))))*param_float(1 + 5)*(qst(nx+2+0,ny+2-1+1,indvarsst(11))*(d2_rhs_u_dydy_1_0_nxp2p0nyp2m1p1k)+&
                    qst(nx+2+0,ny+2-1+1,indvarsst(10))*(d2_rhs_u_dydx_1_0_nxp2p0nyp2m1p1k)))-&
                    q(nx+2+0,ny+2-1+1,indvars(3))*(((1+&
                    param_float(21 + 5))/(((q(nx+2+0,ny+2-1+1,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2-1+1,indvars(2))*q(nx+2+0,ny+2-1+1,indvars(2))+&
                    q(nx+2+0,ny+2-1+1,indvars(3))*q(nx+2+0,ny+2-1+1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2+0,ny+2-1+1,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2-1+1,indvars(2))*q(nx+2+0,ny+2-1+1,indvars(2))+&
                    q(nx+2+0,ny+2-1+1,indvars(3))*q(nx+2+0,ny+2-1+1,indvars(3)))))/param_float(4 + 5)**1.5)*(1+&
                    ((q(nx+2+0,ny+2-1+1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+2+0,ny+2-1+1,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2-1+1,indvars(2))*q(nx+2+0,ny+2-1+1,indvars(2))+&
                    q(nx+2+0,ny+2-1+1,indvars(3))*q(nx+2+0,ny+2-1+1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2+0,ny+2-1+1,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2-1+1,indvars(2))*q(nx+2+0,ny+2-1+1,indvars(2))+&
                    q(nx+2+0,ny+2-1+1,indvars(3))*q(nx+2+0,ny+2-1+1,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+2+0,ny+2-1+1,indvars(1)))**3/((q(nx+2+0,ny+2-1+1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+2+0,ny+2-1+1,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2-1+1,indvars(2))*q(nx+2+0,ny+2-1+1,indvars(2))+&
                    q(nx+2+0,ny+2-1+1,indvars(3))*q(nx+2+0,ny+2-1+1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2+0,ny+2-1+1,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2-1+1,indvars(2))*q(nx+2+0,ny+2-1+1,indvars(2))+&
                    q(nx+2+0,ny+2-1+1,indvars(3))*q(nx+2+0,ny+2-1+1,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+2+0,ny+2-1+1,indvars(1)))**3+&
                    param_float(13 + 5)**3))*(q(nx+2+0,ny+2-1+1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+2+0,ny+2-1+1,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2-1+1,indvars(2))*q(nx+2+0,ny+2-1+1,indvars(2))+&
                    q(nx+2+0,ny+2-1+1,indvars(3))*q(nx+2+0,ny+2-1+1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2+0,ny+2-1+1,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2-1+1,indvars(2))*q(nx+2+0,ny+2-1+1,indvars(2))+&
                    q(nx+2+0,ny+2-1+1,indvars(3))*q(nx+2+0,ny+2-1+1,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+2+0,ny+2-1+1,indvars(1))))*param_float(1 + 5)*(2.0_wp*qst(nx+2+0,ny+2-1+1,indvarsst(11))*(d2_rhs_v_dydy_1_0_nxp2p0nyp2m1p1k)-&
                    2.0_wp/3.0_wp*(qst(nx+2+0,ny+2-1+1,indvarsst(10))*(d2_rhs_v_dydx_1_0_nxp2p0nyp2m1p1k)+&
                    qst(nx+2+0,ny+2-1+1,indvarsst(11))*(d2_rhs_v_dydy_1_0_nxp2p0nyp2m1p1k))))

d1_rhs_et_dy_1_nxp2p0nyp2m1k = -&
          0.5_wp*d1_rhs_et_dy_1_nxp2p0nyp2m1m1k+&
          0.5_wp*d1_rhs_et_dy_1_nxp2p0nyp2m1p1k

d1_rhs_et_dy_1_nxp2p0nyp2m1k = d1_rhs_et_dy_1_nxp2p0nyp2m1k*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer 0 1 None d(rho et)/dt **********
!                                                           
!***********************************************************


rhs(nx+2+0,ny+2-1,indvars(4)) =   -  ( 0.5_wp*(q(nx+2+0,ny+2-1,indvars(2))**2+&
                    q(nx+2+0,ny+2-1,indvars(3))**2)*(1.0_wp/(param_float(23 + 5)*(param_float(3 + 5))*q(nx+2+0,ny+2-1,indvars(1))*((q(nx+2+0,ny+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2-1,indvars(2))*q(nx+2+0,ny+2-1,indvars(2))+&
                    q(nx+2+0,ny+2-1,indvars(3))*q(nx+2+0,ny+2-1,indvars(3)))))/q(nx+2+0,ny+2-1,indvars(1)))**2*((q(nx+2+0,ny+2-1,indvars(2))*(-&
                    (param_float(23 + 5)*(param_float(3 + 5))*q(nx+2+0,ny+2-1,indvars(1))*((q(nx+2+0,ny+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2-1,indvars(2))*q(nx+2+0,ny+2-1,indvars(2))+&
                    q(nx+2+0,ny+2-1,indvars(3))*q(nx+2+0,ny+2-1,indvars(3)))))/q(nx+2+0,ny+2-1,indvars(1)))**2*d1_rhs_rho_dx_1_nxp2p0nyp2m1k+&
                    1.0_wp*d1_rhs_rho_dx_2_nxp2p0nyp2m1k))*qst(nx+2+0,ny+2-1,indvarsst(10))+&
                    0.5_wp*((((param_float(23 + 5)*(param_float(3 + 5))*q(nx+2+0,ny+2-1,indvars(1))*((q(nx+2+0,ny+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2-1,indvars(2))*q(nx+2+0,ny+2-1,indvars(2))+&
                    q(nx+2+0,ny+2-1,indvars(3))*q(nx+2+0,ny+2-1,indvars(3)))))/q(nx+2+0,ny+2-1,indvars(1)))+&
                    q(nx+2+0,ny+2-1,indvars(2)))*((param_float(23 + 5)*(param_float(3 + 5))*q(nx+2+0,ny+2-1,indvars(1))*((q(nx+2+0,ny+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2-1,indvars(2))*q(nx+2+0,ny+2-1,indvars(2))+&
                    q(nx+2+0,ny+2-1,indvars(3))*q(nx+2+0,ny+2-1,indvars(3)))))/q(nx+2+0,ny+2-1,indvars(1)))*q(nx+2+0,ny+2-1,indvars(1))*d1_rhs_rho_dx_3_nxp2p0nyp2m1k+&
                    1.0_wp*d1_rhs_rho_dx_2_nxp2p0nyp2m1k))*qst(nx+2+0,ny+2-1,indvarsst(10))+&
                    -&
                    (((param_float(23 + 5)*(param_float(3 + 5))*q(nx+2+0,ny+2-1,indvars(1))*((q(nx+2+0,ny+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2-1,indvars(2))*q(nx+2+0,ny+2-1,indvars(2))+&
                    q(nx+2+0,ny+2-1,indvars(3))*q(nx+2+0,ny+2-1,indvars(3)))))/q(nx+2+0,ny+2-1,indvars(1)))+&
                    q(nx+2+0,ny+2-1,indvars(2)))*((param_float(23 + 5)*(param_float(3 + 5))*q(nx+2+0,ny+2-1,indvars(1))*((q(nx+2+0,ny+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2-1,indvars(2))*q(nx+2+0,ny+2-1,indvars(2))+&
                    q(nx+2+0,ny+2-1,indvars(3))*q(nx+2+0,ny+2-1,indvars(3)))))/q(nx+2+0,ny+2-1,indvars(1)))*q(nx+2+0,ny+2-1,indvars(1))*d1_rhs_rho_dx_3_nxp2p0nyp2m1k+&
                    1.0_wp*d1_rhs_rho_dx_2_nxp2p0nyp2m1k))*qst(nx+2+0,ny+2-1,indvarsst(10)))))+&
                    1.0_wp/param_float(3 + 5)*(0.5_wp*((((param_float(23 + 5)*(param_float(3 + 5))*q(nx+2+0,ny+2-1,indvars(1))*((q(nx+2+0,ny+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2-1,indvars(2))*q(nx+2+0,ny+2-1,indvars(2))+&
                    q(nx+2+0,ny+2-1,indvars(3))*q(nx+2+0,ny+2-1,indvars(3)))))/q(nx+2+0,ny+2-1,indvars(1)))+&
                    q(nx+2+0,ny+2-1,indvars(2)))*((param_float(23 + 5)*(param_float(3 + 5))*q(nx+2+0,ny+2-1,indvars(1))*((q(nx+2+0,ny+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2-1,indvars(2))*q(nx+2+0,ny+2-1,indvars(2))+&
                    q(nx+2+0,ny+2-1,indvars(3))*q(nx+2+0,ny+2-1,indvars(3)))))/q(nx+2+0,ny+2-1,indvars(1)))*q(nx+2+0,ny+2-1,indvars(1))*d1_rhs_rho_dx_3_nxp2p0nyp2m1k+&
                    1.0_wp*d1_rhs_rho_dx_2_nxp2p0nyp2m1k))*qst(nx+2+0,ny+2-1,indvarsst(10))+&
                    -&
                    (((param_float(23 + 5)*(param_float(3 + 5))*q(nx+2+0,ny+2-1,indvars(1))*((q(nx+2+0,ny+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2-1,indvars(2))*q(nx+2+0,ny+2-1,indvars(2))+&
                    q(nx+2+0,ny+2-1,indvars(3))*q(nx+2+0,ny+2-1,indvars(3)))))/q(nx+2+0,ny+2-1,indvars(1)))+&
                    q(nx+2+0,ny+2-1,indvars(2)))*((param_float(23 + 5)*(param_float(3 + 5))*q(nx+2+0,ny+2-1,indvars(1))*((q(nx+2+0,ny+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2-1,indvars(2))*q(nx+2+0,ny+2-1,indvars(2))+&
                    q(nx+2+0,ny+2-1,indvars(3))*q(nx+2+0,ny+2-1,indvars(3)))))/q(nx+2+0,ny+2-1,indvars(1)))*q(nx+2+0,ny+2-1,indvars(1))*d1_rhs_rho_dx_3_nxp2p0nyp2m1k+&
                    1.0_wp*d1_rhs_rho_dx_2_nxp2p0nyp2m1k))*qst(nx+2+0,ny+2-1,indvarsst(10))))+&
                    q(nx+2+0,ny+2-1,indvars(1))*q(nx+2+0,ny+2-1,indvars(2))*(1.0_wp/(2*q(nx+2+0,ny+2-1,indvars(1))*(param_float(23 + 5)*(param_float(3 + 5))*q(nx+2+0,ny+2-1,indvars(1))*((q(nx+2+0,ny+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2-1,indvars(2))*q(nx+2+0,ny+2-1,indvars(2))+&
                    q(nx+2+0,ny+2-1,indvars(3))*q(nx+2+0,ny+2-1,indvars(3)))))/q(nx+2+0,ny+2-1,indvars(1))))*((((param_float(23 + 5)*(param_float(3 + 5))*q(nx+2+0,ny+2-1,indvars(1))*((q(nx+2+0,ny+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2-1,indvars(2))*q(nx+2+0,ny+2-1,indvars(2))+&
                    q(nx+2+0,ny+2-1,indvars(3))*q(nx+2+0,ny+2-1,indvars(3)))))/q(nx+2+0,ny+2-1,indvars(1)))+&
                    q(nx+2+0,ny+2-1,indvars(2)))*((param_float(23 + 5)*(param_float(3 + 5))*q(nx+2+0,ny+2-1,indvars(1))*((q(nx+2+0,ny+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2-1,indvars(2))*q(nx+2+0,ny+2-1,indvars(2))+&
                    q(nx+2+0,ny+2-1,indvars(3))*q(nx+2+0,ny+2-1,indvars(3)))))/q(nx+2+0,ny+2-1,indvars(1)))*q(nx+2+0,ny+2-1,indvars(1))*d1_rhs_rho_dx_3_nxp2p0nyp2m1k+&
                    1.0_wp*d1_rhs_rho_dx_2_nxp2p0nyp2m1k))*qst(nx+2+0,ny+2-1,indvarsst(10))-&
                    -&
                    (((param_float(23 + 5)*(param_float(3 + 5))*q(nx+2+0,ny+2-1,indvars(1))*((q(nx+2+0,ny+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2-1,indvars(2))*q(nx+2+0,ny+2-1,indvars(2))+&
                    q(nx+2+0,ny+2-1,indvars(3))*q(nx+2+0,ny+2-1,indvars(3)))))/q(nx+2+0,ny+2-1,indvars(1)))+&
                    q(nx+2+0,ny+2-1,indvars(2)))*((param_float(23 + 5)*(param_float(3 + 5))*q(nx+2+0,ny+2-1,indvars(1))*((q(nx+2+0,ny+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2-1,indvars(2))*q(nx+2+0,ny+2-1,indvars(2))+&
                    q(nx+2+0,ny+2-1,indvars(3))*q(nx+2+0,ny+2-1,indvars(3)))))/q(nx+2+0,ny+2-1,indvars(1)))*q(nx+2+0,ny+2-1,indvars(1))*d1_rhs_rho_dx_3_nxp2p0nyp2m1k+&
                    1.0_wp*d1_rhs_rho_dx_2_nxp2p0nyp2m1k))*qst(nx+2+0,ny+2-1,indvarsst(10))))+&
                    q(nx+2+0,ny+2-1,indvars(1))*q(nx+2+0,ny+2-1,indvars(3))*((1.0_wp*q(nx+2+0,ny+2-1,indvars(2))*d1_rhs_v_dx_6_nxp2p0nyp2m1k)*qst(nx+2+0,ny+2-1,indvarsst(10)))+&
                    d1_rhs_et_dy_0_nxp2p0nyp2m1k+&
                    qst(nx+2+0,ny+2-1,indvarsst(10))*(d1_rhs_et_dx_15_nxp2p0nyp2m1k)+&
                    qst(nx+2+0,ny+2-1,indvarsst(11))*(d1_rhs_et_dy_1_nxp2p0nyp2m1k) ) 

