

!***********************************************************
!                                                           
! Start building layers for BC : i1 jmax None **************
!                                                           
!***********************************************************




!***********************************************************
!                                                           
! BC layer: 1 0 None ***************************************
!                                                           
!***********************************************************




!***********************************************************
!                                                           
! building source terms in RHS for layer 1 0 None d(rho)/dt 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! deltaxI*([rho*u]_1x)+deltayI*([rho*v]_1y)+(1.0_wp/c**2*((v*(-c**2*[rho]_1y+1.0_wp*[p]_1y))*deltayI+0.5_wp*(((c+v)*(c*rho*[v]_1y+1.0_wp*[p]_1y))*deltayI+esse*c*(1-M_jmax*M_jmax)/L_ref*(p-P0))))
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_rhs_rho_dx_0_1m2p1m1nyp2p0k = q(1-2+1-1,ny+2+0,indvars(1))*q(1-2+1-1,ny+2+0,indvars(2))

d1_rhs_rho_dx_0_1m2p1p1nyp2p0k = q(1-2+1+1,ny+2+0,indvars(1))*q(1-2+1+1,ny+2+0,indvars(2))

d1_rhs_rho_dx_0_1m2p1nyp2p0k = -&
          0.5_wp*d1_rhs_rho_dx_0_1m2p1m1nyp2p0k+&
          0.5_wp*d1_rhs_rho_dx_0_1m2p1p1nyp2p0k

d1_rhs_rho_dx_0_1m2p1nyp2p0k = d1_rhs_rho_dx_0_1m2p1nyp2p0k*param_float(1)

d1_rhs_rho_dy_0_1m2p1nyp2p0p0k = q(1-2+1,ny+2+0+0,indvars(1))*q(1-2+1,ny+2+0+0,indvars(3))

d1_rhs_rho_dy_0_1m2p1nyp2p0m1k = q(1-2+1,ny+2+0-1,indvars(1))*q(1-2+1,ny+2+0-1,indvars(3))

d1_rhs_rho_dy_0_1m2p1nyp2p0m2k = q(1-2+1,ny+2+0-2,indvars(1))*q(1-2+1,ny+2+0-2,indvars(3))

d1_rhs_rho_dy_0_1m2p1nyp2p0k = 1.5_wp*d1_rhs_rho_dy_0_1m2p1nyp2p0p0k-&
          2.0_wp*d1_rhs_rho_dy_0_1m2p1nyp2p0m1k+&
          0.5_wp*d1_rhs_rho_dy_0_1m2p1nyp2p0m2k

d1_rhs_rho_dy_0_1m2p1nyp2p0k = d1_rhs_rho_dy_0_1m2p1nyp2p0k*param_float(2)

d1_rhs_rho_dy_1_1m2p1nyp2p0p0k = q(1-2+1,ny+2+0+0,indvars(1))

d1_rhs_rho_dy_1_1m2p1nyp2p0m1k = q(1-2+1,ny+2+0-1,indvars(1))

d1_rhs_rho_dy_1_1m2p1nyp2p0m2k = q(1-2+1,ny+2+0-2,indvars(1))

d1_rhs_rho_dy_1_1m2p1nyp2p0k = 1.5_wp*d1_rhs_rho_dy_1_1m2p1nyp2p0p0k-&
          2.0_wp*d1_rhs_rho_dy_1_1m2p1nyp2p0m1k+&
          0.5_wp*d1_rhs_rho_dy_1_1m2p1nyp2p0m2k

d1_rhs_rho_dy_1_1m2p1nyp2p0k = d1_rhs_rho_dy_1_1m2p1nyp2p0k*param_float(2)

d1_rhs_rho_dy_2_1m2p1nyp2p0p0k = (param_float(3 + 5))*q(1-2+1,ny+2+0+0,indvars(1))*((q(1-2+1,ny+2+0+0,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2+0+0,indvars(2))*q(1-2+1,ny+2+0+0,indvars(2))+&
                    q(1-2+1,ny+2+0+0,indvars(3))*q(1-2+1,ny+2+0+0,indvars(3)))))

d1_rhs_rho_dy_2_1m2p1nyp2p0m1k = (param_float(3 + 5))*q(1-2+1,ny+2+0-1,indvars(1))*((q(1-2+1,ny+2+0-1,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2+0-1,indvars(2))*q(1-2+1,ny+2+0-1,indvars(2))+&
                    q(1-2+1,ny+2+0-1,indvars(3))*q(1-2+1,ny+2+0-1,indvars(3)))))

d1_rhs_rho_dy_2_1m2p1nyp2p0m2k = (param_float(3 + 5))*q(1-2+1,ny+2+0-2,indvars(1))*((q(1-2+1,ny+2+0-2,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2+0-2,indvars(2))*q(1-2+1,ny+2+0-2,indvars(2))+&
                    q(1-2+1,ny+2+0-2,indvars(3))*q(1-2+1,ny+2+0-2,indvars(3)))))

d1_rhs_rho_dy_2_1m2p1nyp2p0k = 1.5_wp*d1_rhs_rho_dy_2_1m2p1nyp2p0p0k-&
          2.0_wp*d1_rhs_rho_dy_2_1m2p1nyp2p0m1k+&
          0.5_wp*d1_rhs_rho_dy_2_1m2p1nyp2p0m2k

d1_rhs_rho_dy_2_1m2p1nyp2p0k = d1_rhs_rho_dy_2_1m2p1nyp2p0k*param_float(2)

d1_rhs_rho_dy_3_1m2p1nyp2p0p0k = q(1-2+1,ny+2+0+0,indvars(3))

d1_rhs_rho_dy_3_1m2p1nyp2p0m1k = q(1-2+1,ny+2+0-1,indvars(3))

d1_rhs_rho_dy_3_1m2p1nyp2p0m2k = q(1-2+1,ny+2+0-2,indvars(3))

d1_rhs_rho_dy_3_1m2p1nyp2p0k = 1.5_wp*d1_rhs_rho_dy_3_1m2p1nyp2p0p0k-&
          2.0_wp*d1_rhs_rho_dy_3_1m2p1nyp2p0m1k+&
          0.5_wp*d1_rhs_rho_dy_3_1m2p1nyp2p0m2k

d1_rhs_rho_dy_3_1m2p1nyp2p0k = d1_rhs_rho_dy_3_1m2p1nyp2p0k*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer 1 0 None d(rho)/dt *************
!                                                           
!***********************************************************


rhs(1-2+1,ny+2+0,indvars(1)) =   -  ( qst(1-2+1,ny+2+0,indvarsst(10))*(d1_rhs_rho_dx_0_1m2p1nyp2p0k)+&
                    qst(1-2+1,ny+2+0,indvarsst(11))*(d1_rhs_rho_dy_0_1m2p1nyp2p0k)+&
                    (1.0_wp/(param_float(23 + 5)*(param_float(3 + 5))*q(1-2+1,ny+2+0,indvars(1))*((q(1-2+1,ny+2+0,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2+0,indvars(2))*q(1-2+1,ny+2+0,indvars(2))+&
                    q(1-2+1,ny+2+0,indvars(3))*q(1-2+1,ny+2+0,indvars(3)))))/q(1-2+1,ny+2+0,indvars(1)))**2*((q(1-2+1,ny+2+0,indvars(3))*(-&
                    (param_float(23 + 5)*(param_float(3 + 5))*q(1-2+1,ny+2+0,indvars(1))*((q(1-2+1,ny+2+0,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2+0,indvars(2))*q(1-2+1,ny+2+0,indvars(2))+&
                    q(1-2+1,ny+2+0,indvars(3))*q(1-2+1,ny+2+0,indvars(3)))))/q(1-2+1,ny+2+0,indvars(1)))**2*d1_rhs_rho_dy_1_1m2p1nyp2p0k+&
                    1.0_wp*d1_rhs_rho_dy_2_1m2p1nyp2p0k))*qst(1-2+1,ny+2+0,indvarsst(11))+&
                    0.5_wp*((((param_float(23 + 5)*(param_float(3 + 5))*q(1-2+1,ny+2+0,indvars(1))*((q(1-2+1,ny+2+0,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2+0,indvars(2))*q(1-2+1,ny+2+0,indvars(2))+&
                    q(1-2+1,ny+2+0,indvars(3))*q(1-2+1,ny+2+0,indvars(3)))))/q(1-2+1,ny+2+0,indvars(1)))+&
                    q(1-2+1,ny+2+0,indvars(3)))*((param_float(23 + 5)*(param_float(3 + 5))*q(1-2+1,ny+2+0,indvars(1))*((q(1-2+1,ny+2+0,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2+0,indvars(2))*q(1-2+1,ny+2+0,indvars(2))+&
                    q(1-2+1,ny+2+0,indvars(3))*q(1-2+1,ny+2+0,indvars(3)))))/q(1-2+1,ny+2+0,indvars(1)))*q(1-2+1,ny+2+0,indvars(1))*d1_rhs_rho_dy_3_1m2p1nyp2p0k+&
                    1.0_wp*d1_rhs_rho_dy_2_1m2p1nyp2p0k))*qst(1-2+1,ny+2+0,indvarsst(11))+&
                    param_float(19 + 5)*(param_float(23 + 5)*(param_float(3 + 5))*q(1-2+1,ny+2+0,indvars(1))*((q(1-2+1,ny+2+0,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2+0,indvars(2))*q(1-2+1,ny+2+0,indvars(2))+&
                    q(1-2+1,ny+2+0,indvars(3))*q(1-2+1,ny+2+0,indvars(3)))))/q(1-2+1,ny+2+0,indvars(1)))*(1-&
                    qface_j(1-2+1,1)*qface_j(1-2+1,1))/param_float(20 + 5)*((param_float(3 + 5))*q(1-2+1,ny+2+0,indvars(1))*((q(1-2+1,ny+2+0,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2+0,indvars(2))*q(1-2+1,ny+2+0,indvars(2))+&
                    q(1-2+1,ny+2+0,indvars(3))*q(1-2+1,ny+2+0,indvars(3)))))-&
                    param_float(26 + 5))))) ) 



!***********************************************************
!                                                           
! building source terms in RHS for layer 1 0 None d(rho u)/dt 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! u*(1.0_wp/c**2*((v*(-c**2*[rho]_1y+1.0_wp*[p]_1y))*deltayI+0.5_wp*(((c+v)*(c*rho*[v]_1y+1.0_wp*[p]_1y))*deltayI+esse*c*(1-M_jmax*M_jmax)/L_ref*(p-P0))))+rho*(1.0_wp/(2*rho*c)*(((c+v)*(c*rho*[v]_1y+1.0_wp*[p]_1y))*deltayI-esse*c*(1-M_jmax*M_jmax)/L_ref*(p-P0)))+[rho*u*v]_1y+deltaxI*([-visc_t*(2.0_wp*deltaxI*({u}_1x)-2.0_wp/3.0_wp*(deltaxI*({u}_1x)+deltayI*({v}_1y)))]_1x)+deltayI*([-visc_t*(deltayI*({u}_1y)+deltaxI*({v}_1x))]_1y)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d2_rhs_rhou_dxdx_0_0_1m2p1m1nyp2p0k_1m2p1m1p0nyp2p0k = q(1-2+1-1+0,ny+2+0,indvars(2))

d2_rhs_rhou_dxdx_0_0_1m2p1m1nyp2p0k_1m2p1m1p1nyp2p0k = q(1-2+1-1+1,ny+2+0,indvars(2))

d2_rhs_rhou_dxdx_0_0_1m2p1m1nyp2p0k_1m2p1m1p2nyp2p0k = q(1-2+1-1+2,ny+2+0,indvars(2))

d2_rhs_rhou_dxdx_0_0_1m2p1m1nyp2p0k = -&
          1.5_wp*d2_rhs_rhou_dxdx_0_0_1m2p1m1nyp2p0k_1m2p1m1p0nyp2p0k+&
          2.0_wp*d2_rhs_rhou_dxdx_0_0_1m2p1m1nyp2p0k_1m2p1m1p1nyp2p0k-&
          0.5_wp*d2_rhs_rhou_dxdx_0_0_1m2p1m1nyp2p0k_1m2p1m1p2nyp2p0k

d2_rhs_rhou_dxdx_0_0_1m2p1m1nyp2p0k = d2_rhs_rhou_dxdx_0_0_1m2p1m1nyp2p0k*param_float(1)

d2_rhs_rhou_dxdx_0_0_1m2p1p1nyp2p0k_1m2p1p1m1nyp2p0k = q(1-2+1+1-1,ny+2+0,indvars(2))

d2_rhs_rhou_dxdx_0_0_1m2p1p1nyp2p0k_1m2p1p1p1nyp2p0k = q(1-2+1+1+1,ny+2+0,indvars(2))

d2_rhs_rhou_dxdx_0_0_1m2p1p1nyp2p0k = -&
          0.5_wp*d2_rhs_rhou_dxdx_0_0_1m2p1p1nyp2p0k_1m2p1p1m1nyp2p0k+&
          0.5_wp*d2_rhs_rhou_dxdx_0_0_1m2p1p1nyp2p0k_1m2p1p1p1nyp2p0k

d2_rhs_rhou_dxdx_0_0_1m2p1p1nyp2p0k = d2_rhs_rhou_dxdx_0_0_1m2p1p1nyp2p0k*param_float(1)

d2_rhs_rhou_dxdy_0_0_1m2p1m1nyp2p0k_1m2p1m1nyp2p0p0k = q(1-2+1-1,ny+2+0+0,indvars(3))

d2_rhs_rhou_dxdy_0_0_1m2p1m1nyp2p0k_1m2p1m1nyp2p0m1k = q(1-2+1-1,ny+2+0-1,indvars(3))

d2_rhs_rhou_dxdy_0_0_1m2p1m1nyp2p0k_1m2p1m1nyp2p0m2k = q(1-2+1-1,ny+2+0-2,indvars(3))

d2_rhs_rhou_dxdy_0_0_1m2p1m1nyp2p0k = 1.5_wp*d2_rhs_rhou_dxdy_0_0_1m2p1m1nyp2p0k_1m2p1m1nyp2p0p0k-&
          2.0_wp*d2_rhs_rhou_dxdy_0_0_1m2p1m1nyp2p0k_1m2p1m1nyp2p0m1k+&
          0.5_wp*d2_rhs_rhou_dxdy_0_0_1m2p1m1nyp2p0k_1m2p1m1nyp2p0m2k

d2_rhs_rhou_dxdy_0_0_1m2p1m1nyp2p0k = d2_rhs_rhou_dxdy_0_0_1m2p1m1nyp2p0k*param_float(2)

d2_rhs_rhou_dxdy_0_0_1m2p1p1nyp2p0k_1m2p1p1nyp2p0p0k = q(1-2+1+1,ny+2+0+0,indvars(3))

d2_rhs_rhou_dxdy_0_0_1m2p1p1nyp2p0k_1m2p1p1nyp2p0m1k = q(1-2+1+1,ny+2+0-1,indvars(3))

d2_rhs_rhou_dxdy_0_0_1m2p1p1nyp2p0k_1m2p1p1nyp2p0m2k = q(1-2+1+1,ny+2+0-2,indvars(3))

d2_rhs_rhou_dxdy_0_0_1m2p1p1nyp2p0k = 1.5_wp*d2_rhs_rhou_dxdy_0_0_1m2p1p1nyp2p0k_1m2p1p1nyp2p0p0k-&
          2.0_wp*d2_rhs_rhou_dxdy_0_0_1m2p1p1nyp2p0k_1m2p1p1nyp2p0m1k+&
          0.5_wp*d2_rhs_rhou_dxdy_0_0_1m2p1p1nyp2p0k_1m2p1p1nyp2p0m2k

d2_rhs_rhou_dxdy_0_0_1m2p1p1nyp2p0k = d2_rhs_rhou_dxdy_0_0_1m2p1p1nyp2p0k*param_float(2)

d1_rhs_rhou_dx_0_1m2p1m1nyp2p0k = -((1+&
                    param_float(21 + 5))/(((q(1-2+1-1,ny+2+0,indvars(4))-&
                    0.5_wp*(q(1-2+1-1,ny+2+0,indvars(2))*q(1-2+1-1,ny+2+0,indvars(2))+&
                    q(1-2+1-1,ny+2+0,indvars(3))*q(1-2+1-1,ny+2+0,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(1-2+1-1,ny+2+0,indvars(4))-&
                    0.5_wp*(q(1-2+1-1,ny+2+0,indvars(2))*q(1-2+1-1,ny+2+0,indvars(2))+&
                    q(1-2+1-1,ny+2+0,indvars(3))*q(1-2+1-1,ny+2+0,indvars(3)))))/param_float(4 + 5)**1.5)*(1+&
                    ((q(1-2+1-1,ny+2+0,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(1-2+1-1,ny+2+0,indvars(4))-&
                    0.5_wp*(q(1-2+1-1,ny+2+0,indvars(2))*q(1-2+1-1,ny+2+0,indvars(2))+&
                    q(1-2+1-1,ny+2+0,indvars(3))*q(1-2+1-1,ny+2+0,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(1-2+1-1,ny+2+0,indvars(4))-&
                    0.5_wp*(q(1-2+1-1,ny+2+0,indvars(2))*q(1-2+1-1,ny+2+0,indvars(2))+&
                    q(1-2+1-1,ny+2+0,indvars(3))*q(1-2+1-1,ny+2+0,indvars(3)))))/param_float(4 + 5)**1.5*q(1-2+1-1,ny+2+0,indvars(1)))**3/((q(1-2+1-1,ny+2+0,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(1-2+1-1,ny+2+0,indvars(4))-&
                    0.5_wp*(q(1-2+1-1,ny+2+0,indvars(2))*q(1-2+1-1,ny+2+0,indvars(2))+&
                    q(1-2+1-1,ny+2+0,indvars(3))*q(1-2+1-1,ny+2+0,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(1-2+1-1,ny+2+0,indvars(4))-&
                    0.5_wp*(q(1-2+1-1,ny+2+0,indvars(2))*q(1-2+1-1,ny+2+0,indvars(2))+&
                    q(1-2+1-1,ny+2+0,indvars(3))*q(1-2+1-1,ny+2+0,indvars(3)))))/param_float(4 + 5)**1.5*q(1-2+1-1,ny+2+0,indvars(1)))**3+&
                    param_float(13 + 5)**3))*(q(1-2+1-1,ny+2+0,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(1-2+1-1,ny+2+0,indvars(4))-&
                    0.5_wp*(q(1-2+1-1,ny+2+0,indvars(2))*q(1-2+1-1,ny+2+0,indvars(2))+&
                    q(1-2+1-1,ny+2+0,indvars(3))*q(1-2+1-1,ny+2+0,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(1-2+1-1,ny+2+0,indvars(4))-&
                    0.5_wp*(q(1-2+1-1,ny+2+0,indvars(2))*q(1-2+1-1,ny+2+0,indvars(2))+&
                    q(1-2+1-1,ny+2+0,indvars(3))*q(1-2+1-1,ny+2+0,indvars(3)))))/param_float(4 + 5)**1.5*q(1-2+1-1,ny+2+0,indvars(1))))*param_float(1 + 5)*(2.0_wp*qst(1-2+1-1,ny+2+0,indvarsst(10))*(d2_rhs_rhou_dxdx_0_0_1m2p1m1nyp2p0k)-&
                    2.0_wp/3.0_wp*(qst(1-2+1-1,ny+2+0,indvarsst(10))*(d2_rhs_rhou_dxdx_0_0_1m2p1m1nyp2p0k)+&
                    qst(1-2+1-1,ny+2+0,indvarsst(11))*(d2_rhs_rhou_dxdy_0_0_1m2p1m1nyp2p0k)))

d1_rhs_rhou_dx_0_1m2p1p1nyp2p0k = -((1+&
                    param_float(21 + 5))/(((q(1-2+1+1,ny+2+0,indvars(4))-&
                    0.5_wp*(q(1-2+1+1,ny+2+0,indvars(2))*q(1-2+1+1,ny+2+0,indvars(2))+&
                    q(1-2+1+1,ny+2+0,indvars(3))*q(1-2+1+1,ny+2+0,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(1-2+1+1,ny+2+0,indvars(4))-&
                    0.5_wp*(q(1-2+1+1,ny+2+0,indvars(2))*q(1-2+1+1,ny+2+0,indvars(2))+&
                    q(1-2+1+1,ny+2+0,indvars(3))*q(1-2+1+1,ny+2+0,indvars(3)))))/param_float(4 + 5)**1.5)*(1+&
                    ((q(1-2+1+1,ny+2+0,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(1-2+1+1,ny+2+0,indvars(4))-&
                    0.5_wp*(q(1-2+1+1,ny+2+0,indvars(2))*q(1-2+1+1,ny+2+0,indvars(2))+&
                    q(1-2+1+1,ny+2+0,indvars(3))*q(1-2+1+1,ny+2+0,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(1-2+1+1,ny+2+0,indvars(4))-&
                    0.5_wp*(q(1-2+1+1,ny+2+0,indvars(2))*q(1-2+1+1,ny+2+0,indvars(2))+&
                    q(1-2+1+1,ny+2+0,indvars(3))*q(1-2+1+1,ny+2+0,indvars(3)))))/param_float(4 + 5)**1.5*q(1-2+1+1,ny+2+0,indvars(1)))**3/((q(1-2+1+1,ny+2+0,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(1-2+1+1,ny+2+0,indvars(4))-&
                    0.5_wp*(q(1-2+1+1,ny+2+0,indvars(2))*q(1-2+1+1,ny+2+0,indvars(2))+&
                    q(1-2+1+1,ny+2+0,indvars(3))*q(1-2+1+1,ny+2+0,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(1-2+1+1,ny+2+0,indvars(4))-&
                    0.5_wp*(q(1-2+1+1,ny+2+0,indvars(2))*q(1-2+1+1,ny+2+0,indvars(2))+&
                    q(1-2+1+1,ny+2+0,indvars(3))*q(1-2+1+1,ny+2+0,indvars(3)))))/param_float(4 + 5)**1.5*q(1-2+1+1,ny+2+0,indvars(1)))**3+&
                    param_float(13 + 5)**3))*(q(1-2+1+1,ny+2+0,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(1-2+1+1,ny+2+0,indvars(4))-&
                    0.5_wp*(q(1-2+1+1,ny+2+0,indvars(2))*q(1-2+1+1,ny+2+0,indvars(2))+&
                    q(1-2+1+1,ny+2+0,indvars(3))*q(1-2+1+1,ny+2+0,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(1-2+1+1,ny+2+0,indvars(4))-&
                    0.5_wp*(q(1-2+1+1,ny+2+0,indvars(2))*q(1-2+1+1,ny+2+0,indvars(2))+&
                    q(1-2+1+1,ny+2+0,indvars(3))*q(1-2+1+1,ny+2+0,indvars(3)))))/param_float(4 + 5)**1.5*q(1-2+1+1,ny+2+0,indvars(1))))*param_float(1 + 5)*(2.0_wp*qst(1-2+1+1,ny+2+0,indvarsst(10))*(d2_rhs_rhou_dxdx_0_0_1m2p1p1nyp2p0k)-&
                    2.0_wp/3.0_wp*(qst(1-2+1+1,ny+2+0,indvarsst(10))*(d2_rhs_rhou_dxdx_0_0_1m2p1p1nyp2p0k)+&
                    qst(1-2+1+1,ny+2+0,indvarsst(11))*(d2_rhs_rhou_dxdy_0_0_1m2p1p1nyp2p0k)))

d1_rhs_rhou_dx_0_1m2p1nyp2p0k = -&
          0.5_wp*d1_rhs_rhou_dx_0_1m2p1m1nyp2p0k+&
          0.5_wp*d1_rhs_rhou_dx_0_1m2p1p1nyp2p0k

d1_rhs_rhou_dx_0_1m2p1nyp2p0k = d1_rhs_rhou_dx_0_1m2p1nyp2p0k*param_float(1)

d1_rhs_rhou_dy_6_1m2p1nyp2p0p0k = q(1-2+1,ny+2+0+0,indvars(1))*q(1-2+1,ny+2+0+0,indvars(2))*q(1-2+1,ny+2+0+0,indvars(3))

d1_rhs_rhou_dy_6_1m2p1nyp2p0m1k = q(1-2+1,ny+2+0-1,indvars(1))*q(1-2+1,ny+2+0-1,indvars(2))*q(1-2+1,ny+2+0-1,indvars(3))

d1_rhs_rhou_dy_6_1m2p1nyp2p0m2k = q(1-2+1,ny+2+0-2,indvars(1))*q(1-2+1,ny+2+0-2,indvars(2))*q(1-2+1,ny+2+0-2,indvars(3))

d1_rhs_rhou_dy_6_1m2p1nyp2p0k = 1.5_wp*d1_rhs_rhou_dy_6_1m2p1nyp2p0p0k-&
          2.0_wp*d1_rhs_rhou_dy_6_1m2p1nyp2p0m1k+&
          0.5_wp*d1_rhs_rhou_dy_6_1m2p1nyp2p0m2k

d1_rhs_rhou_dy_6_1m2p1nyp2p0k = d1_rhs_rhou_dy_6_1m2p1nyp2p0k*param_float(2)

d2_rhs_rhou_dydx_7_0_1m2p1nyp2p0p0k_1m2p1m1nyp2p0p0k = q(1-2+1-1,ny+2+0+0,indvars(3))

d2_rhs_rhou_dydx_7_0_1m2p1nyp2p0p0k_1m2p1p1nyp2p0p0k = q(1-2+1+1,ny+2+0+0,indvars(3))

d2_rhs_rhou_dydx_7_0_1m2p1nyp2p0p0k = -&
          0.5_wp*d2_rhs_rhou_dydx_7_0_1m2p1nyp2p0p0k_1m2p1m1nyp2p0p0k+&
          0.5_wp*d2_rhs_rhou_dydx_7_0_1m2p1nyp2p0p0k_1m2p1p1nyp2p0p0k

d2_rhs_rhou_dydx_7_0_1m2p1nyp2p0p0k = d2_rhs_rhou_dydx_7_0_1m2p1nyp2p0p0k*param_float(1)

d2_rhs_rhou_dydx_7_0_1m2p1nyp2p0m1k_1m2p1m1nyp2p0m1k = q(1-2+1-1,ny+2+0-1,indvars(3))

d2_rhs_rhou_dydx_7_0_1m2p1nyp2p0m1k_1m2p1p1nyp2p0m1k = q(1-2+1+1,ny+2+0-1,indvars(3))

d2_rhs_rhou_dydx_7_0_1m2p1nyp2p0m1k = -&
          0.5_wp*d2_rhs_rhou_dydx_7_0_1m2p1nyp2p0m1k_1m2p1m1nyp2p0m1k+&
          0.5_wp*d2_rhs_rhou_dydx_7_0_1m2p1nyp2p0m1k_1m2p1p1nyp2p0m1k

d2_rhs_rhou_dydx_7_0_1m2p1nyp2p0m1k = d2_rhs_rhou_dydx_7_0_1m2p1nyp2p0m1k*param_float(1)

d2_rhs_rhou_dydx_7_0_1m2p1nyp2p0m2k_1m2p1m1nyp2p0m2k = q(1-2+1-1,ny+2+0-2,indvars(3))

d2_rhs_rhou_dydx_7_0_1m2p1nyp2p0m2k_1m2p1p1nyp2p0m2k = q(1-2+1+1,ny+2+0-2,indvars(3))

d2_rhs_rhou_dydx_7_0_1m2p1nyp2p0m2k = -&
          0.5_wp*d2_rhs_rhou_dydx_7_0_1m2p1nyp2p0m2k_1m2p1m1nyp2p0m2k+&
          0.5_wp*d2_rhs_rhou_dydx_7_0_1m2p1nyp2p0m2k_1m2p1p1nyp2p0m2k

d2_rhs_rhou_dydx_7_0_1m2p1nyp2p0m2k = d2_rhs_rhou_dydx_7_0_1m2p1nyp2p0m2k*param_float(1)

d2_rhs_rhou_dydy_7_0_1m2p1nyp2p0p0k_1m2p1nyp2p0p0p0k = q(1-2+1,ny+2+0+0+0,indvars(2))

d2_rhs_rhou_dydy_7_0_1m2p1nyp2p0p0k_1m2p1nyp2p0p0m1k = q(1-2+1,ny+2+0+0-1,indvars(2))

d2_rhs_rhou_dydy_7_0_1m2p1nyp2p0p0k_1m2p1nyp2p0p0m2k = q(1-2+1,ny+2+0+0-2,indvars(2))

d2_rhs_rhou_dydy_7_0_1m2p1nyp2p0p0k = 1.5_wp*d2_rhs_rhou_dydy_7_0_1m2p1nyp2p0p0k_1m2p1nyp2p0p0p0k-&
          2.0_wp*d2_rhs_rhou_dydy_7_0_1m2p1nyp2p0p0k_1m2p1nyp2p0p0m1k+&
          0.5_wp*d2_rhs_rhou_dydy_7_0_1m2p1nyp2p0p0k_1m2p1nyp2p0p0m2k

d2_rhs_rhou_dydy_7_0_1m2p1nyp2p0p0k = d2_rhs_rhou_dydy_7_0_1m2p1nyp2p0p0k*param_float(2)

d2_rhs_rhou_dydy_7_0_1m2p1nyp2p0m1k_1m2p1nyp2p0m1m1k = q(1-2+1,ny+2+0-1-1,indvars(2))

d2_rhs_rhou_dydy_7_0_1m2p1nyp2p0m1k_1m2p1nyp2p0m1p1k = q(1-2+1,ny+2+0-1+1,indvars(2))

d2_rhs_rhou_dydy_7_0_1m2p1nyp2p0m1k = -&
          0.5_wp*d2_rhs_rhou_dydy_7_0_1m2p1nyp2p0m1k_1m2p1nyp2p0m1m1k+&
          0.5_wp*d2_rhs_rhou_dydy_7_0_1m2p1nyp2p0m1k_1m2p1nyp2p0m1p1k

d2_rhs_rhou_dydy_7_0_1m2p1nyp2p0m1k = d2_rhs_rhou_dydy_7_0_1m2p1nyp2p0m1k*param_float(2)

d2_rhs_rhou_dydy_7_0_1m2p1nyp2p0m2k_1m2p1nyp2p0m2m1k = q(1-2+1,ny+2+0-2-1,indvars(2))

d2_rhs_rhou_dydy_7_0_1m2p1nyp2p0m2k_1m2p1nyp2p0m2p1k = q(1-2+1,ny+2+0-2+1,indvars(2))

d2_rhs_rhou_dydy_7_0_1m2p1nyp2p0m2k = -&
          0.5_wp*d2_rhs_rhou_dydy_7_0_1m2p1nyp2p0m2k_1m2p1nyp2p0m2m1k+&
          0.5_wp*d2_rhs_rhou_dydy_7_0_1m2p1nyp2p0m2k_1m2p1nyp2p0m2p1k

d2_rhs_rhou_dydy_7_0_1m2p1nyp2p0m2k = d2_rhs_rhou_dydy_7_0_1m2p1nyp2p0m2k*param_float(2)

d1_rhs_rhou_dy_7_1m2p1nyp2p0p0k = -((1+&
                    param_float(21 + 5))/(((q(1-2+1,ny+2+0+0,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2+0+0,indvars(2))*q(1-2+1,ny+2+0+0,indvars(2))+&
                    q(1-2+1,ny+2+0+0,indvars(3))*q(1-2+1,ny+2+0+0,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(1-2+1,ny+2+0+0,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2+0+0,indvars(2))*q(1-2+1,ny+2+0+0,indvars(2))+&
                    q(1-2+1,ny+2+0+0,indvars(3))*q(1-2+1,ny+2+0+0,indvars(3)))))/param_float(4 + 5)**1.5)*(1+&
                    ((q(1-2+1,ny+2+0+0,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(1-2+1,ny+2+0+0,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2+0+0,indvars(2))*q(1-2+1,ny+2+0+0,indvars(2))+&
                    q(1-2+1,ny+2+0+0,indvars(3))*q(1-2+1,ny+2+0+0,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(1-2+1,ny+2+0+0,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2+0+0,indvars(2))*q(1-2+1,ny+2+0+0,indvars(2))+&
                    q(1-2+1,ny+2+0+0,indvars(3))*q(1-2+1,ny+2+0+0,indvars(3)))))/param_float(4 + 5)**1.5*q(1-2+1,ny+2+0+0,indvars(1)))**3/((q(1-2+1,ny+2+0+0,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(1-2+1,ny+2+0+0,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2+0+0,indvars(2))*q(1-2+1,ny+2+0+0,indvars(2))+&
                    q(1-2+1,ny+2+0+0,indvars(3))*q(1-2+1,ny+2+0+0,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(1-2+1,ny+2+0+0,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2+0+0,indvars(2))*q(1-2+1,ny+2+0+0,indvars(2))+&
                    q(1-2+1,ny+2+0+0,indvars(3))*q(1-2+1,ny+2+0+0,indvars(3)))))/param_float(4 + 5)**1.5*q(1-2+1,ny+2+0+0,indvars(1)))**3+&
                    param_float(13 + 5)**3))*(q(1-2+1,ny+2+0+0,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(1-2+1,ny+2+0+0,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2+0+0,indvars(2))*q(1-2+1,ny+2+0+0,indvars(2))+&
                    q(1-2+1,ny+2+0+0,indvars(3))*q(1-2+1,ny+2+0+0,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(1-2+1,ny+2+0+0,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2+0+0,indvars(2))*q(1-2+1,ny+2+0+0,indvars(2))+&
                    q(1-2+1,ny+2+0+0,indvars(3))*q(1-2+1,ny+2+0+0,indvars(3)))))/param_float(4 + 5)**1.5*q(1-2+1,ny+2+0+0,indvars(1))))*param_float(1 + 5)*(qst(1-2+1,ny+2+0+0,indvarsst(11))*(d2_rhs_rhou_dydy_7_0_1m2p1nyp2p0p0k)+&
                    qst(1-2+1,ny+2+0+0,indvarsst(10))*(d2_rhs_rhou_dydx_7_0_1m2p1nyp2p0p0k))

d1_rhs_rhou_dy_7_1m2p1nyp2p0m1k = -((1+&
                    param_float(21 + 5))/(((q(1-2+1,ny+2+0-1,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2+0-1,indvars(2))*q(1-2+1,ny+2+0-1,indvars(2))+&
                    q(1-2+1,ny+2+0-1,indvars(3))*q(1-2+1,ny+2+0-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(1-2+1,ny+2+0-1,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2+0-1,indvars(2))*q(1-2+1,ny+2+0-1,indvars(2))+&
                    q(1-2+1,ny+2+0-1,indvars(3))*q(1-2+1,ny+2+0-1,indvars(3)))))/param_float(4 + 5)**1.5)*(1+&
                    ((q(1-2+1,ny+2+0-1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(1-2+1,ny+2+0-1,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2+0-1,indvars(2))*q(1-2+1,ny+2+0-1,indvars(2))+&
                    q(1-2+1,ny+2+0-1,indvars(3))*q(1-2+1,ny+2+0-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(1-2+1,ny+2+0-1,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2+0-1,indvars(2))*q(1-2+1,ny+2+0-1,indvars(2))+&
                    q(1-2+1,ny+2+0-1,indvars(3))*q(1-2+1,ny+2+0-1,indvars(3)))))/param_float(4 + 5)**1.5*q(1-2+1,ny+2+0-1,indvars(1)))**3/((q(1-2+1,ny+2+0-1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(1-2+1,ny+2+0-1,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2+0-1,indvars(2))*q(1-2+1,ny+2+0-1,indvars(2))+&
                    q(1-2+1,ny+2+0-1,indvars(3))*q(1-2+1,ny+2+0-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(1-2+1,ny+2+0-1,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2+0-1,indvars(2))*q(1-2+1,ny+2+0-1,indvars(2))+&
                    q(1-2+1,ny+2+0-1,indvars(3))*q(1-2+1,ny+2+0-1,indvars(3)))))/param_float(4 + 5)**1.5*q(1-2+1,ny+2+0-1,indvars(1)))**3+&
                    param_float(13 + 5)**3))*(q(1-2+1,ny+2+0-1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(1-2+1,ny+2+0-1,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2+0-1,indvars(2))*q(1-2+1,ny+2+0-1,indvars(2))+&
                    q(1-2+1,ny+2+0-1,indvars(3))*q(1-2+1,ny+2+0-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(1-2+1,ny+2+0-1,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2+0-1,indvars(2))*q(1-2+1,ny+2+0-1,indvars(2))+&
                    q(1-2+1,ny+2+0-1,indvars(3))*q(1-2+1,ny+2+0-1,indvars(3)))))/param_float(4 + 5)**1.5*q(1-2+1,ny+2+0-1,indvars(1))))*param_float(1 + 5)*(qst(1-2+1,ny+2+0-1,indvarsst(11))*(d2_rhs_rhou_dydy_7_0_1m2p1nyp2p0m1k)+&
                    qst(1-2+1,ny+2+0-1,indvarsst(10))*(d2_rhs_rhou_dydx_7_0_1m2p1nyp2p0m1k))

d1_rhs_rhou_dy_7_1m2p1nyp2p0m2k = -((1+&
                    param_float(21 + 5))/(((q(1-2+1,ny+2+0-2,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2+0-2,indvars(2))*q(1-2+1,ny+2+0-2,indvars(2))+&
                    q(1-2+1,ny+2+0-2,indvars(3))*q(1-2+1,ny+2+0-2,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(1-2+1,ny+2+0-2,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2+0-2,indvars(2))*q(1-2+1,ny+2+0-2,indvars(2))+&
                    q(1-2+1,ny+2+0-2,indvars(3))*q(1-2+1,ny+2+0-2,indvars(3)))))/param_float(4 + 5)**1.5)*(1+&
                    ((q(1-2+1,ny+2+0-2,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(1-2+1,ny+2+0-2,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2+0-2,indvars(2))*q(1-2+1,ny+2+0-2,indvars(2))+&
                    q(1-2+1,ny+2+0-2,indvars(3))*q(1-2+1,ny+2+0-2,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(1-2+1,ny+2+0-2,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2+0-2,indvars(2))*q(1-2+1,ny+2+0-2,indvars(2))+&
                    q(1-2+1,ny+2+0-2,indvars(3))*q(1-2+1,ny+2+0-2,indvars(3)))))/param_float(4 + 5)**1.5*q(1-2+1,ny+2+0-2,indvars(1)))**3/((q(1-2+1,ny+2+0-2,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(1-2+1,ny+2+0-2,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2+0-2,indvars(2))*q(1-2+1,ny+2+0-2,indvars(2))+&
                    q(1-2+1,ny+2+0-2,indvars(3))*q(1-2+1,ny+2+0-2,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(1-2+1,ny+2+0-2,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2+0-2,indvars(2))*q(1-2+1,ny+2+0-2,indvars(2))+&
                    q(1-2+1,ny+2+0-2,indvars(3))*q(1-2+1,ny+2+0-2,indvars(3)))))/param_float(4 + 5)**1.5*q(1-2+1,ny+2+0-2,indvars(1)))**3+&
                    param_float(13 + 5)**3))*(q(1-2+1,ny+2+0-2,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(1-2+1,ny+2+0-2,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2+0-2,indvars(2))*q(1-2+1,ny+2+0-2,indvars(2))+&
                    q(1-2+1,ny+2+0-2,indvars(3))*q(1-2+1,ny+2+0-2,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(1-2+1,ny+2+0-2,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2+0-2,indvars(2))*q(1-2+1,ny+2+0-2,indvars(2))+&
                    q(1-2+1,ny+2+0-2,indvars(3))*q(1-2+1,ny+2+0-2,indvars(3)))))/param_float(4 + 5)**1.5*q(1-2+1,ny+2+0-2,indvars(1))))*param_float(1 + 5)*(qst(1-2+1,ny+2+0-2,indvarsst(11))*(d2_rhs_rhou_dydy_7_0_1m2p1nyp2p0m2k)+&
                    qst(1-2+1,ny+2+0-2,indvarsst(10))*(d2_rhs_rhou_dydx_7_0_1m2p1nyp2p0m2k))

d1_rhs_rhou_dy_7_1m2p1nyp2p0k = 1.5_wp*d1_rhs_rhou_dy_7_1m2p1nyp2p0p0k-&
          2.0_wp*d1_rhs_rhou_dy_7_1m2p1nyp2p0m1k+&
          0.5_wp*d1_rhs_rhou_dy_7_1m2p1nyp2p0m2k

d1_rhs_rhou_dy_7_1m2p1nyp2p0k = d1_rhs_rhou_dy_7_1m2p1nyp2p0k*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer 1 0 None d(rho u)/dt ***********
!                                                           
!***********************************************************


rhs(1-2+1,ny+2+0,indvars(2)) =   -  ( q(1-2+1,ny+2+0,indvars(2))*(1.0_wp/(param_float(23 + 5)*(param_float(3 + 5))*q(1-2+1,ny+2+0,indvars(1))*((q(1-2+1,ny+2+0,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2+0,indvars(2))*q(1-2+1,ny+2+0,indvars(2))+&
                    q(1-2+1,ny+2+0,indvars(3))*q(1-2+1,ny+2+0,indvars(3)))))/q(1-2+1,ny+2+0,indvars(1)))**2*((q(1-2+1,ny+2+0,indvars(3))*(-&
                    (param_float(23 + 5)*(param_float(3 + 5))*q(1-2+1,ny+2+0,indvars(1))*((q(1-2+1,ny+2+0,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2+0,indvars(2))*q(1-2+1,ny+2+0,indvars(2))+&
                    q(1-2+1,ny+2+0,indvars(3))*q(1-2+1,ny+2+0,indvars(3)))))/q(1-2+1,ny+2+0,indvars(1)))**2*d1_rhs_rho_dy_1_1m2p1nyp2p0k+&
                    1.0_wp*d1_rhs_rho_dy_2_1m2p1nyp2p0k))*qst(1-2+1,ny+2+0,indvarsst(11))+&
                    0.5_wp*((((param_float(23 + 5)*(param_float(3 + 5))*q(1-2+1,ny+2+0,indvars(1))*((q(1-2+1,ny+2+0,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2+0,indvars(2))*q(1-2+1,ny+2+0,indvars(2))+&
                    q(1-2+1,ny+2+0,indvars(3))*q(1-2+1,ny+2+0,indvars(3)))))/q(1-2+1,ny+2+0,indvars(1)))+&
                    q(1-2+1,ny+2+0,indvars(3)))*((param_float(23 + 5)*(param_float(3 + 5))*q(1-2+1,ny+2+0,indvars(1))*((q(1-2+1,ny+2+0,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2+0,indvars(2))*q(1-2+1,ny+2+0,indvars(2))+&
                    q(1-2+1,ny+2+0,indvars(3))*q(1-2+1,ny+2+0,indvars(3)))))/q(1-2+1,ny+2+0,indvars(1)))*q(1-2+1,ny+2+0,indvars(1))*d1_rhs_rho_dy_3_1m2p1nyp2p0k+&
                    1.0_wp*d1_rhs_rho_dy_2_1m2p1nyp2p0k))*qst(1-2+1,ny+2+0,indvarsst(11))+&
                    param_float(19 + 5)*(param_float(23 + 5)*(param_float(3 + 5))*q(1-2+1,ny+2+0,indvars(1))*((q(1-2+1,ny+2+0,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2+0,indvars(2))*q(1-2+1,ny+2+0,indvars(2))+&
                    q(1-2+1,ny+2+0,indvars(3))*q(1-2+1,ny+2+0,indvars(3)))))/q(1-2+1,ny+2+0,indvars(1)))*(1-&
                    qface_j(1-2+1,1)*qface_j(1-2+1,1))/param_float(20 + 5)*((param_float(3 + 5))*q(1-2+1,ny+2+0,indvars(1))*((q(1-2+1,ny+2+0,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2+0,indvars(2))*q(1-2+1,ny+2+0,indvars(2))+&
                    q(1-2+1,ny+2+0,indvars(3))*q(1-2+1,ny+2+0,indvars(3)))))-&
                    param_float(26 + 5)))))+&
                    q(1-2+1,ny+2+0,indvars(1))*(1.0_wp/(2*q(1-2+1,ny+2+0,indvars(1))*(param_float(23 + 5)*(param_float(3 + 5))*q(1-2+1,ny+2+0,indvars(1))*((q(1-2+1,ny+2+0,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2+0,indvars(2))*q(1-2+1,ny+2+0,indvars(2))+&
                    q(1-2+1,ny+2+0,indvars(3))*q(1-2+1,ny+2+0,indvars(3)))))/q(1-2+1,ny+2+0,indvars(1))))*((((param_float(23 + 5)*(param_float(3 + 5))*q(1-2+1,ny+2+0,indvars(1))*((q(1-2+1,ny+2+0,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2+0,indvars(2))*q(1-2+1,ny+2+0,indvars(2))+&
                    q(1-2+1,ny+2+0,indvars(3))*q(1-2+1,ny+2+0,indvars(3)))))/q(1-2+1,ny+2+0,indvars(1)))+&
                    q(1-2+1,ny+2+0,indvars(3)))*((param_float(23 + 5)*(param_float(3 + 5))*q(1-2+1,ny+2+0,indvars(1))*((q(1-2+1,ny+2+0,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2+0,indvars(2))*q(1-2+1,ny+2+0,indvars(2))+&
                    q(1-2+1,ny+2+0,indvars(3))*q(1-2+1,ny+2+0,indvars(3)))))/q(1-2+1,ny+2+0,indvars(1)))*q(1-2+1,ny+2+0,indvars(1))*d1_rhs_rho_dy_3_1m2p1nyp2p0k+&
                    1.0_wp*d1_rhs_rho_dy_2_1m2p1nyp2p0k))*qst(1-2+1,ny+2+0,indvarsst(11))-&
                    param_float(19 + 5)*(param_float(23 + 5)*(param_float(3 + 5))*q(1-2+1,ny+2+0,indvars(1))*((q(1-2+1,ny+2+0,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2+0,indvars(2))*q(1-2+1,ny+2+0,indvars(2))+&
                    q(1-2+1,ny+2+0,indvars(3))*q(1-2+1,ny+2+0,indvars(3)))))/q(1-2+1,ny+2+0,indvars(1)))*(1-&
                    qface_j(1-2+1,1)*qface_j(1-2+1,1))/param_float(20 + 5)*((param_float(3 + 5))*q(1-2+1,ny+2+0,indvars(1))*((q(1-2+1,ny+2+0,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2+0,indvars(2))*q(1-2+1,ny+2+0,indvars(2))+&
                    q(1-2+1,ny+2+0,indvars(3))*q(1-2+1,ny+2+0,indvars(3)))))-&
                    param_float(26 + 5))))+&
                    d1_rhs_rhou_dy_6_1m2p1nyp2p0k+&
                    qst(1-2+1,ny+2+0,indvarsst(10))*(d1_rhs_rhou_dx_0_1m2p1nyp2p0k)+&
                    qst(1-2+1,ny+2+0,indvarsst(11))*(d1_rhs_rhou_dy_7_1m2p1nyp2p0k) ) 



!***********************************************************
!                                                           
! building source terms in RHS for layer 1 0 None d(rho v)/dt 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! v*(1.0_wp/c**2*((v*(-c**2*[rho]_1y+1.0_wp*[p]_1y))*deltayI+0.5_wp*(((c+v)*(c*rho*[v]_1y+1.0_wp*[p]_1y))*deltayI+esse*c*(1-M_jmax*M_jmax)/L_ref*(p-P0))))+rho*((1.0_wp*v*[u]_1y)*deltayI)+[rho*v*v]_1y+deltaxI*([-visc_t*(deltayI*({u}_1y)+deltaxI*({v}_1x))]_1x)+deltayI*([-visc_t*(2.0_wp*deltayI*({v}_1y)-2.0_wp/3.0_wp*(deltaxI*({u}_1x)+deltayI*({v}_1y)))]_1y)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d2_rhs_rhov_dxdx_0_0_1m2p1m1nyp2p0k_1m2p1m1p0nyp2p0k = q(1-2+1-1+0,ny+2+0,indvars(3))

d2_rhs_rhov_dxdx_0_0_1m2p1m1nyp2p0k_1m2p1m1p1nyp2p0k = q(1-2+1-1+1,ny+2+0,indvars(3))

d2_rhs_rhov_dxdx_0_0_1m2p1m1nyp2p0k_1m2p1m1p2nyp2p0k = q(1-2+1-1+2,ny+2+0,indvars(3))

d2_rhs_rhov_dxdx_0_0_1m2p1m1nyp2p0k = -&
          1.5_wp*d2_rhs_rhov_dxdx_0_0_1m2p1m1nyp2p0k_1m2p1m1p0nyp2p0k+&
          2.0_wp*d2_rhs_rhov_dxdx_0_0_1m2p1m1nyp2p0k_1m2p1m1p1nyp2p0k-&
          0.5_wp*d2_rhs_rhov_dxdx_0_0_1m2p1m1nyp2p0k_1m2p1m1p2nyp2p0k

d2_rhs_rhov_dxdx_0_0_1m2p1m1nyp2p0k = d2_rhs_rhov_dxdx_0_0_1m2p1m1nyp2p0k*param_float(1)

d2_rhs_rhov_dxdx_0_0_1m2p1p1nyp2p0k_1m2p1p1m1nyp2p0k = q(1-2+1+1-1,ny+2+0,indvars(3))

d2_rhs_rhov_dxdx_0_0_1m2p1p1nyp2p0k_1m2p1p1p1nyp2p0k = q(1-2+1+1+1,ny+2+0,indvars(3))

d2_rhs_rhov_dxdx_0_0_1m2p1p1nyp2p0k = -&
          0.5_wp*d2_rhs_rhov_dxdx_0_0_1m2p1p1nyp2p0k_1m2p1p1m1nyp2p0k+&
          0.5_wp*d2_rhs_rhov_dxdx_0_0_1m2p1p1nyp2p0k_1m2p1p1p1nyp2p0k

d2_rhs_rhov_dxdx_0_0_1m2p1p1nyp2p0k = d2_rhs_rhov_dxdx_0_0_1m2p1p1nyp2p0k*param_float(1)

d2_rhs_rhov_dxdy_0_0_1m2p1m1nyp2p0k_1m2p1m1nyp2p0p0k = q(1-2+1-1,ny+2+0+0,indvars(2))

d2_rhs_rhov_dxdy_0_0_1m2p1m1nyp2p0k_1m2p1m1nyp2p0m1k = q(1-2+1-1,ny+2+0-1,indvars(2))

d2_rhs_rhov_dxdy_0_0_1m2p1m1nyp2p0k_1m2p1m1nyp2p0m2k = q(1-2+1-1,ny+2+0-2,indvars(2))

d2_rhs_rhov_dxdy_0_0_1m2p1m1nyp2p0k = 1.5_wp*d2_rhs_rhov_dxdy_0_0_1m2p1m1nyp2p0k_1m2p1m1nyp2p0p0k-&
          2.0_wp*d2_rhs_rhov_dxdy_0_0_1m2p1m1nyp2p0k_1m2p1m1nyp2p0m1k+&
          0.5_wp*d2_rhs_rhov_dxdy_0_0_1m2p1m1nyp2p0k_1m2p1m1nyp2p0m2k

d2_rhs_rhov_dxdy_0_0_1m2p1m1nyp2p0k = d2_rhs_rhov_dxdy_0_0_1m2p1m1nyp2p0k*param_float(2)

d2_rhs_rhov_dxdy_0_0_1m2p1p1nyp2p0k_1m2p1p1nyp2p0p0k = q(1-2+1+1,ny+2+0+0,indvars(2))

d2_rhs_rhov_dxdy_0_0_1m2p1p1nyp2p0k_1m2p1p1nyp2p0m1k = q(1-2+1+1,ny+2+0-1,indvars(2))

d2_rhs_rhov_dxdy_0_0_1m2p1p1nyp2p0k_1m2p1p1nyp2p0m2k = q(1-2+1+1,ny+2+0-2,indvars(2))

d2_rhs_rhov_dxdy_0_0_1m2p1p1nyp2p0k = 1.5_wp*d2_rhs_rhov_dxdy_0_0_1m2p1p1nyp2p0k_1m2p1p1nyp2p0p0k-&
          2.0_wp*d2_rhs_rhov_dxdy_0_0_1m2p1p1nyp2p0k_1m2p1p1nyp2p0m1k+&
          0.5_wp*d2_rhs_rhov_dxdy_0_0_1m2p1p1nyp2p0k_1m2p1p1nyp2p0m2k

d2_rhs_rhov_dxdy_0_0_1m2p1p1nyp2p0k = d2_rhs_rhov_dxdy_0_0_1m2p1p1nyp2p0k*param_float(2)

d1_rhs_rhov_dx_0_1m2p1m1nyp2p0k = -((1+&
                    param_float(21 + 5))/(((q(1-2+1-1,ny+2+0,indvars(4))-&
                    0.5_wp*(q(1-2+1-1,ny+2+0,indvars(2))*q(1-2+1-1,ny+2+0,indvars(2))+&
                    q(1-2+1-1,ny+2+0,indvars(3))*q(1-2+1-1,ny+2+0,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(1-2+1-1,ny+2+0,indvars(4))-&
                    0.5_wp*(q(1-2+1-1,ny+2+0,indvars(2))*q(1-2+1-1,ny+2+0,indvars(2))+&
                    q(1-2+1-1,ny+2+0,indvars(3))*q(1-2+1-1,ny+2+0,indvars(3)))))/param_float(4 + 5)**1.5)*(1+&
                    ((q(1-2+1-1,ny+2+0,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(1-2+1-1,ny+2+0,indvars(4))-&
                    0.5_wp*(q(1-2+1-1,ny+2+0,indvars(2))*q(1-2+1-1,ny+2+0,indvars(2))+&
                    q(1-2+1-1,ny+2+0,indvars(3))*q(1-2+1-1,ny+2+0,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(1-2+1-1,ny+2+0,indvars(4))-&
                    0.5_wp*(q(1-2+1-1,ny+2+0,indvars(2))*q(1-2+1-1,ny+2+0,indvars(2))+&
                    q(1-2+1-1,ny+2+0,indvars(3))*q(1-2+1-1,ny+2+0,indvars(3)))))/param_float(4 + 5)**1.5*q(1-2+1-1,ny+2+0,indvars(1)))**3/((q(1-2+1-1,ny+2+0,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(1-2+1-1,ny+2+0,indvars(4))-&
                    0.5_wp*(q(1-2+1-1,ny+2+0,indvars(2))*q(1-2+1-1,ny+2+0,indvars(2))+&
                    q(1-2+1-1,ny+2+0,indvars(3))*q(1-2+1-1,ny+2+0,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(1-2+1-1,ny+2+0,indvars(4))-&
                    0.5_wp*(q(1-2+1-1,ny+2+0,indvars(2))*q(1-2+1-1,ny+2+0,indvars(2))+&
                    q(1-2+1-1,ny+2+0,indvars(3))*q(1-2+1-1,ny+2+0,indvars(3)))))/param_float(4 + 5)**1.5*q(1-2+1-1,ny+2+0,indvars(1)))**3+&
                    param_float(13 + 5)**3))*(q(1-2+1-1,ny+2+0,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(1-2+1-1,ny+2+0,indvars(4))-&
                    0.5_wp*(q(1-2+1-1,ny+2+0,indvars(2))*q(1-2+1-1,ny+2+0,indvars(2))+&
                    q(1-2+1-1,ny+2+0,indvars(3))*q(1-2+1-1,ny+2+0,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(1-2+1-1,ny+2+0,indvars(4))-&
                    0.5_wp*(q(1-2+1-1,ny+2+0,indvars(2))*q(1-2+1-1,ny+2+0,indvars(2))+&
                    q(1-2+1-1,ny+2+0,indvars(3))*q(1-2+1-1,ny+2+0,indvars(3)))))/param_float(4 + 5)**1.5*q(1-2+1-1,ny+2+0,indvars(1))))*param_float(1 + 5)*(qst(1-2+1-1,ny+2+0,indvarsst(11))*(d2_rhs_rhov_dxdy_0_0_1m2p1m1nyp2p0k)+&
                    qst(1-2+1-1,ny+2+0,indvarsst(10))*(d2_rhs_rhov_dxdx_0_0_1m2p1m1nyp2p0k))

d1_rhs_rhov_dx_0_1m2p1p1nyp2p0k = -((1+&
                    param_float(21 + 5))/(((q(1-2+1+1,ny+2+0,indvars(4))-&
                    0.5_wp*(q(1-2+1+1,ny+2+0,indvars(2))*q(1-2+1+1,ny+2+0,indvars(2))+&
                    q(1-2+1+1,ny+2+0,indvars(3))*q(1-2+1+1,ny+2+0,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(1-2+1+1,ny+2+0,indvars(4))-&
                    0.5_wp*(q(1-2+1+1,ny+2+0,indvars(2))*q(1-2+1+1,ny+2+0,indvars(2))+&
                    q(1-2+1+1,ny+2+0,indvars(3))*q(1-2+1+1,ny+2+0,indvars(3)))))/param_float(4 + 5)**1.5)*(1+&
                    ((q(1-2+1+1,ny+2+0,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(1-2+1+1,ny+2+0,indvars(4))-&
                    0.5_wp*(q(1-2+1+1,ny+2+0,indvars(2))*q(1-2+1+1,ny+2+0,indvars(2))+&
                    q(1-2+1+1,ny+2+0,indvars(3))*q(1-2+1+1,ny+2+0,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(1-2+1+1,ny+2+0,indvars(4))-&
                    0.5_wp*(q(1-2+1+1,ny+2+0,indvars(2))*q(1-2+1+1,ny+2+0,indvars(2))+&
                    q(1-2+1+1,ny+2+0,indvars(3))*q(1-2+1+1,ny+2+0,indvars(3)))))/param_float(4 + 5)**1.5*q(1-2+1+1,ny+2+0,indvars(1)))**3/((q(1-2+1+1,ny+2+0,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(1-2+1+1,ny+2+0,indvars(4))-&
                    0.5_wp*(q(1-2+1+1,ny+2+0,indvars(2))*q(1-2+1+1,ny+2+0,indvars(2))+&
                    q(1-2+1+1,ny+2+0,indvars(3))*q(1-2+1+1,ny+2+0,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(1-2+1+1,ny+2+0,indvars(4))-&
                    0.5_wp*(q(1-2+1+1,ny+2+0,indvars(2))*q(1-2+1+1,ny+2+0,indvars(2))+&
                    q(1-2+1+1,ny+2+0,indvars(3))*q(1-2+1+1,ny+2+0,indvars(3)))))/param_float(4 + 5)**1.5*q(1-2+1+1,ny+2+0,indvars(1)))**3+&
                    param_float(13 + 5)**3))*(q(1-2+1+1,ny+2+0,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(1-2+1+1,ny+2+0,indvars(4))-&
                    0.5_wp*(q(1-2+1+1,ny+2+0,indvars(2))*q(1-2+1+1,ny+2+0,indvars(2))+&
                    q(1-2+1+1,ny+2+0,indvars(3))*q(1-2+1+1,ny+2+0,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(1-2+1+1,ny+2+0,indvars(4))-&
                    0.5_wp*(q(1-2+1+1,ny+2+0,indvars(2))*q(1-2+1+1,ny+2+0,indvars(2))+&
                    q(1-2+1+1,ny+2+0,indvars(3))*q(1-2+1+1,ny+2+0,indvars(3)))))/param_float(4 + 5)**1.5*q(1-2+1+1,ny+2+0,indvars(1))))*param_float(1 + 5)*(qst(1-2+1+1,ny+2+0,indvarsst(11))*(d2_rhs_rhov_dxdy_0_0_1m2p1p1nyp2p0k)+&
                    qst(1-2+1+1,ny+2+0,indvarsst(10))*(d2_rhs_rhov_dxdx_0_0_1m2p1p1nyp2p0k))

d1_rhs_rhov_dx_0_1m2p1nyp2p0k = -&
          0.5_wp*d1_rhs_rhov_dx_0_1m2p1m1nyp2p0k+&
          0.5_wp*d1_rhs_rhov_dx_0_1m2p1p1nyp2p0k

d1_rhs_rhov_dx_0_1m2p1nyp2p0k = d1_rhs_rhov_dx_0_1m2p1nyp2p0k*param_float(1)

d1_rhs_rhov_dy_4_1m2p1nyp2p0p0k = q(1-2+1,ny+2+0+0,indvars(2))

d1_rhs_rhov_dy_4_1m2p1nyp2p0m1k = q(1-2+1,ny+2+0-1,indvars(2))

d1_rhs_rhov_dy_4_1m2p1nyp2p0m2k = q(1-2+1,ny+2+0-2,indvars(2))

d1_rhs_rhov_dy_4_1m2p1nyp2p0k = 1.5_wp*d1_rhs_rhov_dy_4_1m2p1nyp2p0p0k-&
          2.0_wp*d1_rhs_rhov_dy_4_1m2p1nyp2p0m1k+&
          0.5_wp*d1_rhs_rhov_dy_4_1m2p1nyp2p0m2k

d1_rhs_rhov_dy_4_1m2p1nyp2p0k = d1_rhs_rhov_dy_4_1m2p1nyp2p0k*param_float(2)

d1_rhs_rhov_dy_5_1m2p1nyp2p0p0k = q(1-2+1,ny+2+0+0,indvars(1))*q(1-2+1,ny+2+0+0,indvars(3))*q(1-2+1,ny+2+0+0,indvars(3))

d1_rhs_rhov_dy_5_1m2p1nyp2p0m1k = q(1-2+1,ny+2+0-1,indvars(1))*q(1-2+1,ny+2+0-1,indvars(3))*q(1-2+1,ny+2+0-1,indvars(3))

d1_rhs_rhov_dy_5_1m2p1nyp2p0m2k = q(1-2+1,ny+2+0-2,indvars(1))*q(1-2+1,ny+2+0-2,indvars(3))*q(1-2+1,ny+2+0-2,indvars(3))

d1_rhs_rhov_dy_5_1m2p1nyp2p0k = 1.5_wp*d1_rhs_rhov_dy_5_1m2p1nyp2p0p0k-&
          2.0_wp*d1_rhs_rhov_dy_5_1m2p1nyp2p0m1k+&
          0.5_wp*d1_rhs_rhov_dy_5_1m2p1nyp2p0m2k

d1_rhs_rhov_dy_5_1m2p1nyp2p0k = d1_rhs_rhov_dy_5_1m2p1nyp2p0k*param_float(2)

d2_rhs_rhov_dydx_6_0_1m2p1nyp2p0p0k_1m2p1m1nyp2p0p0k = q(1-2+1-1,ny+2+0+0,indvars(2))

d2_rhs_rhov_dydx_6_0_1m2p1nyp2p0p0k_1m2p1p1nyp2p0p0k = q(1-2+1+1,ny+2+0+0,indvars(2))

d2_rhs_rhov_dydx_6_0_1m2p1nyp2p0p0k = -&
          0.5_wp*d2_rhs_rhov_dydx_6_0_1m2p1nyp2p0p0k_1m2p1m1nyp2p0p0k+&
          0.5_wp*d2_rhs_rhov_dydx_6_0_1m2p1nyp2p0p0k_1m2p1p1nyp2p0p0k

d2_rhs_rhov_dydx_6_0_1m2p1nyp2p0p0k = d2_rhs_rhov_dydx_6_0_1m2p1nyp2p0p0k*param_float(1)

d2_rhs_rhov_dydx_6_0_1m2p1nyp2p0m1k_1m2p1m1nyp2p0m1k = q(1-2+1-1,ny+2+0-1,indvars(2))

d2_rhs_rhov_dydx_6_0_1m2p1nyp2p0m1k_1m2p1p1nyp2p0m1k = q(1-2+1+1,ny+2+0-1,indvars(2))

d2_rhs_rhov_dydx_6_0_1m2p1nyp2p0m1k = -&
          0.5_wp*d2_rhs_rhov_dydx_6_0_1m2p1nyp2p0m1k_1m2p1m1nyp2p0m1k+&
          0.5_wp*d2_rhs_rhov_dydx_6_0_1m2p1nyp2p0m1k_1m2p1p1nyp2p0m1k

d2_rhs_rhov_dydx_6_0_1m2p1nyp2p0m1k = d2_rhs_rhov_dydx_6_0_1m2p1nyp2p0m1k*param_float(1)

d2_rhs_rhov_dydx_6_0_1m2p1nyp2p0m2k_1m2p1m1nyp2p0m2k = q(1-2+1-1,ny+2+0-2,indvars(2))

d2_rhs_rhov_dydx_6_0_1m2p1nyp2p0m2k_1m2p1p1nyp2p0m2k = q(1-2+1+1,ny+2+0-2,indvars(2))

d2_rhs_rhov_dydx_6_0_1m2p1nyp2p0m2k = -&
          0.5_wp*d2_rhs_rhov_dydx_6_0_1m2p1nyp2p0m2k_1m2p1m1nyp2p0m2k+&
          0.5_wp*d2_rhs_rhov_dydx_6_0_1m2p1nyp2p0m2k_1m2p1p1nyp2p0m2k

d2_rhs_rhov_dydx_6_0_1m2p1nyp2p0m2k = d2_rhs_rhov_dydx_6_0_1m2p1nyp2p0m2k*param_float(1)

d2_rhs_rhov_dydy_6_0_1m2p1nyp2p0p0k_1m2p1nyp2p0p0p0k = q(1-2+1,ny+2+0+0+0,indvars(3))

d2_rhs_rhov_dydy_6_0_1m2p1nyp2p0p0k_1m2p1nyp2p0p0m1k = q(1-2+1,ny+2+0+0-1,indvars(3))

d2_rhs_rhov_dydy_6_0_1m2p1nyp2p0p0k_1m2p1nyp2p0p0m2k = q(1-2+1,ny+2+0+0-2,indvars(3))

d2_rhs_rhov_dydy_6_0_1m2p1nyp2p0p0k = 1.5_wp*d2_rhs_rhov_dydy_6_0_1m2p1nyp2p0p0k_1m2p1nyp2p0p0p0k-&
          2.0_wp*d2_rhs_rhov_dydy_6_0_1m2p1nyp2p0p0k_1m2p1nyp2p0p0m1k+&
          0.5_wp*d2_rhs_rhov_dydy_6_0_1m2p1nyp2p0p0k_1m2p1nyp2p0p0m2k

d2_rhs_rhov_dydy_6_0_1m2p1nyp2p0p0k = d2_rhs_rhov_dydy_6_0_1m2p1nyp2p0p0k*param_float(2)

d2_rhs_rhov_dydy_6_0_1m2p1nyp2p0m1k_1m2p1nyp2p0m1m1k = q(1-2+1,ny+2+0-1-1,indvars(3))

d2_rhs_rhov_dydy_6_0_1m2p1nyp2p0m1k_1m2p1nyp2p0m1p1k = q(1-2+1,ny+2+0-1+1,indvars(3))

d2_rhs_rhov_dydy_6_0_1m2p1nyp2p0m1k = -&
          0.5_wp*d2_rhs_rhov_dydy_6_0_1m2p1nyp2p0m1k_1m2p1nyp2p0m1m1k+&
          0.5_wp*d2_rhs_rhov_dydy_6_0_1m2p1nyp2p0m1k_1m2p1nyp2p0m1p1k

d2_rhs_rhov_dydy_6_0_1m2p1nyp2p0m1k = d2_rhs_rhov_dydy_6_0_1m2p1nyp2p0m1k*param_float(2)

d2_rhs_rhov_dydy_6_0_1m2p1nyp2p0m2k_1m2p1nyp2p0m2m1k = q(1-2+1,ny+2+0-2-1,indvars(3))

d2_rhs_rhov_dydy_6_0_1m2p1nyp2p0m2k_1m2p1nyp2p0m2p1k = q(1-2+1,ny+2+0-2+1,indvars(3))

d2_rhs_rhov_dydy_6_0_1m2p1nyp2p0m2k = -&
          0.5_wp*d2_rhs_rhov_dydy_6_0_1m2p1nyp2p0m2k_1m2p1nyp2p0m2m1k+&
          0.5_wp*d2_rhs_rhov_dydy_6_0_1m2p1nyp2p0m2k_1m2p1nyp2p0m2p1k

d2_rhs_rhov_dydy_6_0_1m2p1nyp2p0m2k = d2_rhs_rhov_dydy_6_0_1m2p1nyp2p0m2k*param_float(2)

d1_rhs_rhov_dy_6_1m2p1nyp2p0p0k = -((1+&
                    param_float(21 + 5))/(((q(1-2+1,ny+2+0+0,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2+0+0,indvars(2))*q(1-2+1,ny+2+0+0,indvars(2))+&
                    q(1-2+1,ny+2+0+0,indvars(3))*q(1-2+1,ny+2+0+0,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(1-2+1,ny+2+0+0,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2+0+0,indvars(2))*q(1-2+1,ny+2+0+0,indvars(2))+&
                    q(1-2+1,ny+2+0+0,indvars(3))*q(1-2+1,ny+2+0+0,indvars(3)))))/param_float(4 + 5)**1.5)*(1+&
                    ((q(1-2+1,ny+2+0+0,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(1-2+1,ny+2+0+0,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2+0+0,indvars(2))*q(1-2+1,ny+2+0+0,indvars(2))+&
                    q(1-2+1,ny+2+0+0,indvars(3))*q(1-2+1,ny+2+0+0,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(1-2+1,ny+2+0+0,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2+0+0,indvars(2))*q(1-2+1,ny+2+0+0,indvars(2))+&
                    q(1-2+1,ny+2+0+0,indvars(3))*q(1-2+1,ny+2+0+0,indvars(3)))))/param_float(4 + 5)**1.5*q(1-2+1,ny+2+0+0,indvars(1)))**3/((q(1-2+1,ny+2+0+0,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(1-2+1,ny+2+0+0,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2+0+0,indvars(2))*q(1-2+1,ny+2+0+0,indvars(2))+&
                    q(1-2+1,ny+2+0+0,indvars(3))*q(1-2+1,ny+2+0+0,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(1-2+1,ny+2+0+0,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2+0+0,indvars(2))*q(1-2+1,ny+2+0+0,indvars(2))+&
                    q(1-2+1,ny+2+0+0,indvars(3))*q(1-2+1,ny+2+0+0,indvars(3)))))/param_float(4 + 5)**1.5*q(1-2+1,ny+2+0+0,indvars(1)))**3+&
                    param_float(13 + 5)**3))*(q(1-2+1,ny+2+0+0,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(1-2+1,ny+2+0+0,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2+0+0,indvars(2))*q(1-2+1,ny+2+0+0,indvars(2))+&
                    q(1-2+1,ny+2+0+0,indvars(3))*q(1-2+1,ny+2+0+0,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(1-2+1,ny+2+0+0,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2+0+0,indvars(2))*q(1-2+1,ny+2+0+0,indvars(2))+&
                    q(1-2+1,ny+2+0+0,indvars(3))*q(1-2+1,ny+2+0+0,indvars(3)))))/param_float(4 + 5)**1.5*q(1-2+1,ny+2+0+0,indvars(1))))*param_float(1 + 5)*(2.0_wp*qst(1-2+1,ny+2+0+0,indvarsst(11))*(d2_rhs_rhov_dydy_6_0_1m2p1nyp2p0p0k)-&
                    2.0_wp/3.0_wp*(qst(1-2+1,ny+2+0+0,indvarsst(10))*(d2_rhs_rhov_dydx_6_0_1m2p1nyp2p0p0k)+&
                    qst(1-2+1,ny+2+0+0,indvarsst(11))*(d2_rhs_rhov_dydy_6_0_1m2p1nyp2p0p0k)))

d1_rhs_rhov_dy_6_1m2p1nyp2p0m1k = -((1+&
                    param_float(21 + 5))/(((q(1-2+1,ny+2+0-1,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2+0-1,indvars(2))*q(1-2+1,ny+2+0-1,indvars(2))+&
                    q(1-2+1,ny+2+0-1,indvars(3))*q(1-2+1,ny+2+0-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(1-2+1,ny+2+0-1,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2+0-1,indvars(2))*q(1-2+1,ny+2+0-1,indvars(2))+&
                    q(1-2+1,ny+2+0-1,indvars(3))*q(1-2+1,ny+2+0-1,indvars(3)))))/param_float(4 + 5)**1.5)*(1+&
                    ((q(1-2+1,ny+2+0-1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(1-2+1,ny+2+0-1,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2+0-1,indvars(2))*q(1-2+1,ny+2+0-1,indvars(2))+&
                    q(1-2+1,ny+2+0-1,indvars(3))*q(1-2+1,ny+2+0-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(1-2+1,ny+2+0-1,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2+0-1,indvars(2))*q(1-2+1,ny+2+0-1,indvars(2))+&
                    q(1-2+1,ny+2+0-1,indvars(3))*q(1-2+1,ny+2+0-1,indvars(3)))))/param_float(4 + 5)**1.5*q(1-2+1,ny+2+0-1,indvars(1)))**3/((q(1-2+1,ny+2+0-1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(1-2+1,ny+2+0-1,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2+0-1,indvars(2))*q(1-2+1,ny+2+0-1,indvars(2))+&
                    q(1-2+1,ny+2+0-1,indvars(3))*q(1-2+1,ny+2+0-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(1-2+1,ny+2+0-1,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2+0-1,indvars(2))*q(1-2+1,ny+2+0-1,indvars(2))+&
                    q(1-2+1,ny+2+0-1,indvars(3))*q(1-2+1,ny+2+0-1,indvars(3)))))/param_float(4 + 5)**1.5*q(1-2+1,ny+2+0-1,indvars(1)))**3+&
                    param_float(13 + 5)**3))*(q(1-2+1,ny+2+0-1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(1-2+1,ny+2+0-1,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2+0-1,indvars(2))*q(1-2+1,ny+2+0-1,indvars(2))+&
                    q(1-2+1,ny+2+0-1,indvars(3))*q(1-2+1,ny+2+0-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(1-2+1,ny+2+0-1,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2+0-1,indvars(2))*q(1-2+1,ny+2+0-1,indvars(2))+&
                    q(1-2+1,ny+2+0-1,indvars(3))*q(1-2+1,ny+2+0-1,indvars(3)))))/param_float(4 + 5)**1.5*q(1-2+1,ny+2+0-1,indvars(1))))*param_float(1 + 5)*(2.0_wp*qst(1-2+1,ny+2+0-1,indvarsst(11))*(d2_rhs_rhov_dydy_6_0_1m2p1nyp2p0m1k)-&
                    2.0_wp/3.0_wp*(qst(1-2+1,ny+2+0-1,indvarsst(10))*(d2_rhs_rhov_dydx_6_0_1m2p1nyp2p0m1k)+&
                    qst(1-2+1,ny+2+0-1,indvarsst(11))*(d2_rhs_rhov_dydy_6_0_1m2p1nyp2p0m1k)))

d1_rhs_rhov_dy_6_1m2p1nyp2p0m2k = -((1+&
                    param_float(21 + 5))/(((q(1-2+1,ny+2+0-2,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2+0-2,indvars(2))*q(1-2+1,ny+2+0-2,indvars(2))+&
                    q(1-2+1,ny+2+0-2,indvars(3))*q(1-2+1,ny+2+0-2,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(1-2+1,ny+2+0-2,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2+0-2,indvars(2))*q(1-2+1,ny+2+0-2,indvars(2))+&
                    q(1-2+1,ny+2+0-2,indvars(3))*q(1-2+1,ny+2+0-2,indvars(3)))))/param_float(4 + 5)**1.5)*(1+&
                    ((q(1-2+1,ny+2+0-2,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(1-2+1,ny+2+0-2,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2+0-2,indvars(2))*q(1-2+1,ny+2+0-2,indvars(2))+&
                    q(1-2+1,ny+2+0-2,indvars(3))*q(1-2+1,ny+2+0-2,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(1-2+1,ny+2+0-2,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2+0-2,indvars(2))*q(1-2+1,ny+2+0-2,indvars(2))+&
                    q(1-2+1,ny+2+0-2,indvars(3))*q(1-2+1,ny+2+0-2,indvars(3)))))/param_float(4 + 5)**1.5*q(1-2+1,ny+2+0-2,indvars(1)))**3/((q(1-2+1,ny+2+0-2,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(1-2+1,ny+2+0-2,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2+0-2,indvars(2))*q(1-2+1,ny+2+0-2,indvars(2))+&
                    q(1-2+1,ny+2+0-2,indvars(3))*q(1-2+1,ny+2+0-2,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(1-2+1,ny+2+0-2,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2+0-2,indvars(2))*q(1-2+1,ny+2+0-2,indvars(2))+&
                    q(1-2+1,ny+2+0-2,indvars(3))*q(1-2+1,ny+2+0-2,indvars(3)))))/param_float(4 + 5)**1.5*q(1-2+1,ny+2+0-2,indvars(1)))**3+&
                    param_float(13 + 5)**3))*(q(1-2+1,ny+2+0-2,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(1-2+1,ny+2+0-2,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2+0-2,indvars(2))*q(1-2+1,ny+2+0-2,indvars(2))+&
                    q(1-2+1,ny+2+0-2,indvars(3))*q(1-2+1,ny+2+0-2,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(1-2+1,ny+2+0-2,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2+0-2,indvars(2))*q(1-2+1,ny+2+0-2,indvars(2))+&
                    q(1-2+1,ny+2+0-2,indvars(3))*q(1-2+1,ny+2+0-2,indvars(3)))))/param_float(4 + 5)**1.5*q(1-2+1,ny+2+0-2,indvars(1))))*param_float(1 + 5)*(2.0_wp*qst(1-2+1,ny+2+0-2,indvarsst(11))*(d2_rhs_rhov_dydy_6_0_1m2p1nyp2p0m2k)-&
                    2.0_wp/3.0_wp*(qst(1-2+1,ny+2+0-2,indvarsst(10))*(d2_rhs_rhov_dydx_6_0_1m2p1nyp2p0m2k)+&
                    qst(1-2+1,ny+2+0-2,indvarsst(11))*(d2_rhs_rhov_dydy_6_0_1m2p1nyp2p0m2k)))

d1_rhs_rhov_dy_6_1m2p1nyp2p0k = 1.5_wp*d1_rhs_rhov_dy_6_1m2p1nyp2p0p0k-&
          2.0_wp*d1_rhs_rhov_dy_6_1m2p1nyp2p0m1k+&
          0.5_wp*d1_rhs_rhov_dy_6_1m2p1nyp2p0m2k

d1_rhs_rhov_dy_6_1m2p1nyp2p0k = d1_rhs_rhov_dy_6_1m2p1nyp2p0k*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer 1 0 None d(rho v)/dt ***********
!                                                           
!***********************************************************


rhs(1-2+1,ny+2+0,indvars(3)) =   -  ( q(1-2+1,ny+2+0,indvars(3))*(1.0_wp/(param_float(23 + 5)*(param_float(3 + 5))*q(1-2+1,ny+2+0,indvars(1))*((q(1-2+1,ny+2+0,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2+0,indvars(2))*q(1-2+1,ny+2+0,indvars(2))+&
                    q(1-2+1,ny+2+0,indvars(3))*q(1-2+1,ny+2+0,indvars(3)))))/q(1-2+1,ny+2+0,indvars(1)))**2*((q(1-2+1,ny+2+0,indvars(3))*(-&
                    (param_float(23 + 5)*(param_float(3 + 5))*q(1-2+1,ny+2+0,indvars(1))*((q(1-2+1,ny+2+0,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2+0,indvars(2))*q(1-2+1,ny+2+0,indvars(2))+&
                    q(1-2+1,ny+2+0,indvars(3))*q(1-2+1,ny+2+0,indvars(3)))))/q(1-2+1,ny+2+0,indvars(1)))**2*d1_rhs_rho_dy_1_1m2p1nyp2p0k+&
                    1.0_wp*d1_rhs_rho_dy_2_1m2p1nyp2p0k))*qst(1-2+1,ny+2+0,indvarsst(11))+&
                    0.5_wp*((((param_float(23 + 5)*(param_float(3 + 5))*q(1-2+1,ny+2+0,indvars(1))*((q(1-2+1,ny+2+0,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2+0,indvars(2))*q(1-2+1,ny+2+0,indvars(2))+&
                    q(1-2+1,ny+2+0,indvars(3))*q(1-2+1,ny+2+0,indvars(3)))))/q(1-2+1,ny+2+0,indvars(1)))+&
                    q(1-2+1,ny+2+0,indvars(3)))*((param_float(23 + 5)*(param_float(3 + 5))*q(1-2+1,ny+2+0,indvars(1))*((q(1-2+1,ny+2+0,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2+0,indvars(2))*q(1-2+1,ny+2+0,indvars(2))+&
                    q(1-2+1,ny+2+0,indvars(3))*q(1-2+1,ny+2+0,indvars(3)))))/q(1-2+1,ny+2+0,indvars(1)))*q(1-2+1,ny+2+0,indvars(1))*d1_rhs_rho_dy_3_1m2p1nyp2p0k+&
                    1.0_wp*d1_rhs_rho_dy_2_1m2p1nyp2p0k))*qst(1-2+1,ny+2+0,indvarsst(11))+&
                    param_float(19 + 5)*(param_float(23 + 5)*(param_float(3 + 5))*q(1-2+1,ny+2+0,indvars(1))*((q(1-2+1,ny+2+0,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2+0,indvars(2))*q(1-2+1,ny+2+0,indvars(2))+&
                    q(1-2+1,ny+2+0,indvars(3))*q(1-2+1,ny+2+0,indvars(3)))))/q(1-2+1,ny+2+0,indvars(1)))*(1-&
                    qface_j(1-2+1,1)*qface_j(1-2+1,1))/param_float(20 + 5)*((param_float(3 + 5))*q(1-2+1,ny+2+0,indvars(1))*((q(1-2+1,ny+2+0,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2+0,indvars(2))*q(1-2+1,ny+2+0,indvars(2))+&
                    q(1-2+1,ny+2+0,indvars(3))*q(1-2+1,ny+2+0,indvars(3)))))-&
                    param_float(26 + 5)))))+&
                    q(1-2+1,ny+2+0,indvars(1))*((1.0_wp*q(1-2+1,ny+2+0,indvars(3))*d1_rhs_rhov_dy_4_1m2p1nyp2p0k)*qst(1-2+1,ny+2+0,indvarsst(11)))+&
                    d1_rhs_rhov_dy_5_1m2p1nyp2p0k+&
                    qst(1-2+1,ny+2+0,indvarsst(10))*(d1_rhs_rhov_dx_0_1m2p1nyp2p0k)+&
                    qst(1-2+1,ny+2+0,indvarsst(11))*(d1_rhs_rhov_dy_6_1m2p1nyp2p0k) ) 



!***********************************************************
!                                                           
! building source terms in RHS for layer 1 0 None d(rho et)/dt 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! 0.5_wp*(u**2+v**2)*(1.0_wp/c**2*((v*(-c**2*[rho]_1y+1.0_wp*[p]_1y))*deltayI+0.5_wp*(((c+v)*(c*rho*[v]_1y+1.0_wp*[p]_1y))*deltayI+esse*c*(1-M_jmax*M_jmax)/L_ref*(p-P0))))+1.0_wp/gamma_m1*(0.5_wp*(((c+v)*(c*rho*[v]_1y+1.0_wp*[p]_1y))*deltayI+esse*c*(1-M_jmax*M_jmax)/L_ref*(p-P0)))+rho*u*(1.0_wp/(2*rho*c)*(((c+v)*(c*rho*[v]_1y+1.0_wp*[p]_1y))*deltayI-esse*c*(1-M_jmax*M_jmax)/L_ref*(p-P0)))+rho*v*((1.0_wp*v*[u]_1y)*deltayI)+[(rho*et+p)*v]_1y+deltaxI*([-kappa*deltaxI*({T}_1x)-u*(visc_t*(2.0_wp*deltaxI*({u}_1x)-2.0_wp/3.0_wp*(deltaxI*({u}_1x)+deltayI*({v}_1y))))-v*(visc_t*(deltayI*({u}_1y)+deltaxI*({v}_1x)))]_1x)+deltayI*([-kappa*deltayI*({T}_1y)-u*(visc_t*(deltayI*({u}_1y)+deltaxI*({v}_1x)))-v*(visc_t*(2.0_wp*deltayI*({v}_1y)-2.0_wp/3.0_wp*(deltaxI*({u}_1x)+deltayI*({v}_1y))))]_1y)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d2_rhs_et_dxdx_0_0_1m2p1m1nyp2p0k_1m2p1m1p0nyp2p0k = ((q(1-2+1-1+0,ny+2+0,indvars(4))-&
                    0.5_wp*(q(1-2+1-1+0,ny+2+0,indvars(2))*q(1-2+1-1+0,ny+2+0,indvars(2))+&
                    q(1-2+1-1+0,ny+2+0,indvars(3))*q(1-2+1-1+0,ny+2+0,indvars(3)))))/param_float(4 + 5)

d2_rhs_et_dxdx_0_0_1m2p1m1nyp2p0k_1m2p1m1p1nyp2p0k = ((q(1-2+1-1+1,ny+2+0,indvars(4))-&
                    0.5_wp*(q(1-2+1-1+1,ny+2+0,indvars(2))*q(1-2+1-1+1,ny+2+0,indvars(2))+&
                    q(1-2+1-1+1,ny+2+0,indvars(3))*q(1-2+1-1+1,ny+2+0,indvars(3)))))/param_float(4 + 5)

d2_rhs_et_dxdx_0_0_1m2p1m1nyp2p0k_1m2p1m1p2nyp2p0k = ((q(1-2+1-1+2,ny+2+0,indvars(4))-&
                    0.5_wp*(q(1-2+1-1+2,ny+2+0,indvars(2))*q(1-2+1-1+2,ny+2+0,indvars(2))+&
                    q(1-2+1-1+2,ny+2+0,indvars(3))*q(1-2+1-1+2,ny+2+0,indvars(3)))))/param_float(4 + 5)

d2_rhs_et_dxdx_0_0_1m2p1m1nyp2p0k = -&
          1.5_wp*d2_rhs_et_dxdx_0_0_1m2p1m1nyp2p0k_1m2p1m1p0nyp2p0k+&
          2.0_wp*d2_rhs_et_dxdx_0_0_1m2p1m1nyp2p0k_1m2p1m1p1nyp2p0k-&
          0.5_wp*d2_rhs_et_dxdx_0_0_1m2p1m1nyp2p0k_1m2p1m1p2nyp2p0k

d2_rhs_et_dxdx_0_0_1m2p1m1nyp2p0k = d2_rhs_et_dxdx_0_0_1m2p1m1nyp2p0k*param_float(1)

d2_rhs_et_dxdx_0_0_1m2p1p1nyp2p0k_1m2p1p1m1nyp2p0k = ((q(1-2+1+1-1,ny+2+0,indvars(4))-&
                    0.5_wp*(q(1-2+1+1-1,ny+2+0,indvars(2))*q(1-2+1+1-1,ny+2+0,indvars(2))+&
                    q(1-2+1+1-1,ny+2+0,indvars(3))*q(1-2+1+1-1,ny+2+0,indvars(3)))))/param_float(4 + 5)

d2_rhs_et_dxdx_0_0_1m2p1p1nyp2p0k_1m2p1p1p1nyp2p0k = ((q(1-2+1+1+1,ny+2+0,indvars(4))-&
                    0.5_wp*(q(1-2+1+1+1,ny+2+0,indvars(2))*q(1-2+1+1+1,ny+2+0,indvars(2))+&
                    q(1-2+1+1+1,ny+2+0,indvars(3))*q(1-2+1+1+1,ny+2+0,indvars(3)))))/param_float(4 + 5)

d2_rhs_et_dxdx_0_0_1m2p1p1nyp2p0k = -&
          0.5_wp*d2_rhs_et_dxdx_0_0_1m2p1p1nyp2p0k_1m2p1p1m1nyp2p0k+&
          0.5_wp*d2_rhs_et_dxdx_0_0_1m2p1p1nyp2p0k_1m2p1p1p1nyp2p0k

d2_rhs_et_dxdx_0_0_1m2p1p1nyp2p0k = d2_rhs_et_dxdx_0_0_1m2p1p1nyp2p0k*param_float(1)

d1_rhs_et_dx_0_1m2p1m1nyp2p0k = -param_float(2 + 5)*qst(1-2+1-1,ny+2+0,indvarsst(10))*(d2_rhs_et_dxdx_0_0_1m2p1m1nyp2p0k)-&
                    q(1-2+1-1,ny+2+0,indvars(2))*(((1+&
                    param_float(21 + 5))/(((q(1-2+1-1,ny+2+0,indvars(4))-&
                    0.5_wp*(q(1-2+1-1,ny+2+0,indvars(2))*q(1-2+1-1,ny+2+0,indvars(2))+&
                    q(1-2+1-1,ny+2+0,indvars(3))*q(1-2+1-1,ny+2+0,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(1-2+1-1,ny+2+0,indvars(4))-&
                    0.5_wp*(q(1-2+1-1,ny+2+0,indvars(2))*q(1-2+1-1,ny+2+0,indvars(2))+&
                    q(1-2+1-1,ny+2+0,indvars(3))*q(1-2+1-1,ny+2+0,indvars(3)))))/param_float(4 + 5)**1.5)*(1+&
                    ((q(1-2+1-1,ny+2+0,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(1-2+1-1,ny+2+0,indvars(4))-&
                    0.5_wp*(q(1-2+1-1,ny+2+0,indvars(2))*q(1-2+1-1,ny+2+0,indvars(2))+&
                    q(1-2+1-1,ny+2+0,indvars(3))*q(1-2+1-1,ny+2+0,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(1-2+1-1,ny+2+0,indvars(4))-&
                    0.5_wp*(q(1-2+1-1,ny+2+0,indvars(2))*q(1-2+1-1,ny+2+0,indvars(2))+&
                    q(1-2+1-1,ny+2+0,indvars(3))*q(1-2+1-1,ny+2+0,indvars(3)))))/param_float(4 + 5)**1.5*q(1-2+1-1,ny+2+0,indvars(1)))**3/((q(1-2+1-1,ny+2+0,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(1-2+1-1,ny+2+0,indvars(4))-&
                    0.5_wp*(q(1-2+1-1,ny+2+0,indvars(2))*q(1-2+1-1,ny+2+0,indvars(2))+&
                    q(1-2+1-1,ny+2+0,indvars(3))*q(1-2+1-1,ny+2+0,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(1-2+1-1,ny+2+0,indvars(4))-&
                    0.5_wp*(q(1-2+1-1,ny+2+0,indvars(2))*q(1-2+1-1,ny+2+0,indvars(2))+&
                    q(1-2+1-1,ny+2+0,indvars(3))*q(1-2+1-1,ny+2+0,indvars(3)))))/param_float(4 + 5)**1.5*q(1-2+1-1,ny+2+0,indvars(1)))**3+&
                    param_float(13 + 5)**3))*(q(1-2+1-1,ny+2+0,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(1-2+1-1,ny+2+0,indvars(4))-&
                    0.5_wp*(q(1-2+1-1,ny+2+0,indvars(2))*q(1-2+1-1,ny+2+0,indvars(2))+&
                    q(1-2+1-1,ny+2+0,indvars(3))*q(1-2+1-1,ny+2+0,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(1-2+1-1,ny+2+0,indvars(4))-&
                    0.5_wp*(q(1-2+1-1,ny+2+0,indvars(2))*q(1-2+1-1,ny+2+0,indvars(2))+&
                    q(1-2+1-1,ny+2+0,indvars(3))*q(1-2+1-1,ny+2+0,indvars(3)))))/param_float(4 + 5)**1.5*q(1-2+1-1,ny+2+0,indvars(1))))*param_float(1 + 5)*(2.0_wp*qst(1-2+1-1,ny+2+0,indvarsst(10))*(d2_rhs_rhou_dxdx_0_0_1m2p1m1nyp2p0k)-&
                    2.0_wp/3.0_wp*(qst(1-2+1-1,ny+2+0,indvarsst(10))*(d2_rhs_rhou_dxdx_0_0_1m2p1m1nyp2p0k)+&
                    qst(1-2+1-1,ny+2+0,indvarsst(11))*(d2_rhs_rhou_dxdy_0_0_1m2p1m1nyp2p0k))))-&
                    q(1-2+1-1,ny+2+0,indvars(3))*(((1+&
                    param_float(21 + 5))/(((q(1-2+1-1,ny+2+0,indvars(4))-&
                    0.5_wp*(q(1-2+1-1,ny+2+0,indvars(2))*q(1-2+1-1,ny+2+0,indvars(2))+&
                    q(1-2+1-1,ny+2+0,indvars(3))*q(1-2+1-1,ny+2+0,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(1-2+1-1,ny+2+0,indvars(4))-&
                    0.5_wp*(q(1-2+1-1,ny+2+0,indvars(2))*q(1-2+1-1,ny+2+0,indvars(2))+&
                    q(1-2+1-1,ny+2+0,indvars(3))*q(1-2+1-1,ny+2+0,indvars(3)))))/param_float(4 + 5)**1.5)*(1+&
                    ((q(1-2+1-1,ny+2+0,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(1-2+1-1,ny+2+0,indvars(4))-&
                    0.5_wp*(q(1-2+1-1,ny+2+0,indvars(2))*q(1-2+1-1,ny+2+0,indvars(2))+&
                    q(1-2+1-1,ny+2+0,indvars(3))*q(1-2+1-1,ny+2+0,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(1-2+1-1,ny+2+0,indvars(4))-&
                    0.5_wp*(q(1-2+1-1,ny+2+0,indvars(2))*q(1-2+1-1,ny+2+0,indvars(2))+&
                    q(1-2+1-1,ny+2+0,indvars(3))*q(1-2+1-1,ny+2+0,indvars(3)))))/param_float(4 + 5)**1.5*q(1-2+1-1,ny+2+0,indvars(1)))**3/((q(1-2+1-1,ny+2+0,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(1-2+1-1,ny+2+0,indvars(4))-&
                    0.5_wp*(q(1-2+1-1,ny+2+0,indvars(2))*q(1-2+1-1,ny+2+0,indvars(2))+&
                    q(1-2+1-1,ny+2+0,indvars(3))*q(1-2+1-1,ny+2+0,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(1-2+1-1,ny+2+0,indvars(4))-&
                    0.5_wp*(q(1-2+1-1,ny+2+0,indvars(2))*q(1-2+1-1,ny+2+0,indvars(2))+&
                    q(1-2+1-1,ny+2+0,indvars(3))*q(1-2+1-1,ny+2+0,indvars(3)))))/param_float(4 + 5)**1.5*q(1-2+1-1,ny+2+0,indvars(1)))**3+&
                    param_float(13 + 5)**3))*(q(1-2+1-1,ny+2+0,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(1-2+1-1,ny+2+0,indvars(4))-&
                    0.5_wp*(q(1-2+1-1,ny+2+0,indvars(2))*q(1-2+1-1,ny+2+0,indvars(2))+&
                    q(1-2+1-1,ny+2+0,indvars(3))*q(1-2+1-1,ny+2+0,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(1-2+1-1,ny+2+0,indvars(4))-&
                    0.5_wp*(q(1-2+1-1,ny+2+0,indvars(2))*q(1-2+1-1,ny+2+0,indvars(2))+&
                    q(1-2+1-1,ny+2+0,indvars(3))*q(1-2+1-1,ny+2+0,indvars(3)))))/param_float(4 + 5)**1.5*q(1-2+1-1,ny+2+0,indvars(1))))*param_float(1 + 5)*(qst(1-2+1-1,ny+2+0,indvarsst(11))*(d2_rhs_rhov_dxdy_0_0_1m2p1m1nyp2p0k)+&
                    qst(1-2+1-1,ny+2+0,indvarsst(10))*(d2_rhs_rhov_dxdx_0_0_1m2p1m1nyp2p0k)))

d1_rhs_et_dx_0_1m2p1p1nyp2p0k = -param_float(2 + 5)*qst(1-2+1+1,ny+2+0,indvarsst(10))*(d2_rhs_et_dxdx_0_0_1m2p1p1nyp2p0k)-&
                    q(1-2+1+1,ny+2+0,indvars(2))*(((1+&
                    param_float(21 + 5))/(((q(1-2+1+1,ny+2+0,indvars(4))-&
                    0.5_wp*(q(1-2+1+1,ny+2+0,indvars(2))*q(1-2+1+1,ny+2+0,indvars(2))+&
                    q(1-2+1+1,ny+2+0,indvars(3))*q(1-2+1+1,ny+2+0,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(1-2+1+1,ny+2+0,indvars(4))-&
                    0.5_wp*(q(1-2+1+1,ny+2+0,indvars(2))*q(1-2+1+1,ny+2+0,indvars(2))+&
                    q(1-2+1+1,ny+2+0,indvars(3))*q(1-2+1+1,ny+2+0,indvars(3)))))/param_float(4 + 5)**1.5)*(1+&
                    ((q(1-2+1+1,ny+2+0,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(1-2+1+1,ny+2+0,indvars(4))-&
                    0.5_wp*(q(1-2+1+1,ny+2+0,indvars(2))*q(1-2+1+1,ny+2+0,indvars(2))+&
                    q(1-2+1+1,ny+2+0,indvars(3))*q(1-2+1+1,ny+2+0,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(1-2+1+1,ny+2+0,indvars(4))-&
                    0.5_wp*(q(1-2+1+1,ny+2+0,indvars(2))*q(1-2+1+1,ny+2+0,indvars(2))+&
                    q(1-2+1+1,ny+2+0,indvars(3))*q(1-2+1+1,ny+2+0,indvars(3)))))/param_float(4 + 5)**1.5*q(1-2+1+1,ny+2+0,indvars(1)))**3/((q(1-2+1+1,ny+2+0,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(1-2+1+1,ny+2+0,indvars(4))-&
                    0.5_wp*(q(1-2+1+1,ny+2+0,indvars(2))*q(1-2+1+1,ny+2+0,indvars(2))+&
                    q(1-2+1+1,ny+2+0,indvars(3))*q(1-2+1+1,ny+2+0,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(1-2+1+1,ny+2+0,indvars(4))-&
                    0.5_wp*(q(1-2+1+1,ny+2+0,indvars(2))*q(1-2+1+1,ny+2+0,indvars(2))+&
                    q(1-2+1+1,ny+2+0,indvars(3))*q(1-2+1+1,ny+2+0,indvars(3)))))/param_float(4 + 5)**1.5*q(1-2+1+1,ny+2+0,indvars(1)))**3+&
                    param_float(13 + 5)**3))*(q(1-2+1+1,ny+2+0,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(1-2+1+1,ny+2+0,indvars(4))-&
                    0.5_wp*(q(1-2+1+1,ny+2+0,indvars(2))*q(1-2+1+1,ny+2+0,indvars(2))+&
                    q(1-2+1+1,ny+2+0,indvars(3))*q(1-2+1+1,ny+2+0,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(1-2+1+1,ny+2+0,indvars(4))-&
                    0.5_wp*(q(1-2+1+1,ny+2+0,indvars(2))*q(1-2+1+1,ny+2+0,indvars(2))+&
                    q(1-2+1+1,ny+2+0,indvars(3))*q(1-2+1+1,ny+2+0,indvars(3)))))/param_float(4 + 5)**1.5*q(1-2+1+1,ny+2+0,indvars(1))))*param_float(1 + 5)*(2.0_wp*qst(1-2+1+1,ny+2+0,indvarsst(10))*(d2_rhs_rhou_dxdx_0_0_1m2p1p1nyp2p0k)-&
                    2.0_wp/3.0_wp*(qst(1-2+1+1,ny+2+0,indvarsst(10))*(d2_rhs_rhou_dxdx_0_0_1m2p1p1nyp2p0k)+&
                    qst(1-2+1+1,ny+2+0,indvarsst(11))*(d2_rhs_rhou_dxdy_0_0_1m2p1p1nyp2p0k))))-&
                    q(1-2+1+1,ny+2+0,indvars(3))*(((1+&
                    param_float(21 + 5))/(((q(1-2+1+1,ny+2+0,indvars(4))-&
                    0.5_wp*(q(1-2+1+1,ny+2+0,indvars(2))*q(1-2+1+1,ny+2+0,indvars(2))+&
                    q(1-2+1+1,ny+2+0,indvars(3))*q(1-2+1+1,ny+2+0,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(1-2+1+1,ny+2+0,indvars(4))-&
                    0.5_wp*(q(1-2+1+1,ny+2+0,indvars(2))*q(1-2+1+1,ny+2+0,indvars(2))+&
                    q(1-2+1+1,ny+2+0,indvars(3))*q(1-2+1+1,ny+2+0,indvars(3)))))/param_float(4 + 5)**1.5)*(1+&
                    ((q(1-2+1+1,ny+2+0,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(1-2+1+1,ny+2+0,indvars(4))-&
                    0.5_wp*(q(1-2+1+1,ny+2+0,indvars(2))*q(1-2+1+1,ny+2+0,indvars(2))+&
                    q(1-2+1+1,ny+2+0,indvars(3))*q(1-2+1+1,ny+2+0,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(1-2+1+1,ny+2+0,indvars(4))-&
                    0.5_wp*(q(1-2+1+1,ny+2+0,indvars(2))*q(1-2+1+1,ny+2+0,indvars(2))+&
                    q(1-2+1+1,ny+2+0,indvars(3))*q(1-2+1+1,ny+2+0,indvars(3)))))/param_float(4 + 5)**1.5*q(1-2+1+1,ny+2+0,indvars(1)))**3/((q(1-2+1+1,ny+2+0,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(1-2+1+1,ny+2+0,indvars(4))-&
                    0.5_wp*(q(1-2+1+1,ny+2+0,indvars(2))*q(1-2+1+1,ny+2+0,indvars(2))+&
                    q(1-2+1+1,ny+2+0,indvars(3))*q(1-2+1+1,ny+2+0,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(1-2+1+1,ny+2+0,indvars(4))-&
                    0.5_wp*(q(1-2+1+1,ny+2+0,indvars(2))*q(1-2+1+1,ny+2+0,indvars(2))+&
                    q(1-2+1+1,ny+2+0,indvars(3))*q(1-2+1+1,ny+2+0,indvars(3)))))/param_float(4 + 5)**1.5*q(1-2+1+1,ny+2+0,indvars(1)))**3+&
                    param_float(13 + 5)**3))*(q(1-2+1+1,ny+2+0,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(1-2+1+1,ny+2+0,indvars(4))-&
                    0.5_wp*(q(1-2+1+1,ny+2+0,indvars(2))*q(1-2+1+1,ny+2+0,indvars(2))+&
                    q(1-2+1+1,ny+2+0,indvars(3))*q(1-2+1+1,ny+2+0,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(1-2+1+1,ny+2+0,indvars(4))-&
                    0.5_wp*(q(1-2+1+1,ny+2+0,indvars(2))*q(1-2+1+1,ny+2+0,indvars(2))+&
                    q(1-2+1+1,ny+2+0,indvars(3))*q(1-2+1+1,ny+2+0,indvars(3)))))/param_float(4 + 5)**1.5*q(1-2+1+1,ny+2+0,indvars(1))))*param_float(1 + 5)*(qst(1-2+1+1,ny+2+0,indvarsst(11))*(d2_rhs_rhov_dxdy_0_0_1m2p1p1nyp2p0k)+&
                    qst(1-2+1+1,ny+2+0,indvarsst(10))*(d2_rhs_rhov_dxdx_0_0_1m2p1p1nyp2p0k)))

d1_rhs_et_dx_0_1m2p1nyp2p0k = -&
          0.5_wp*d1_rhs_et_dx_0_1m2p1m1nyp2p0k+&
          0.5_wp*d1_rhs_et_dx_0_1m2p1p1nyp2p0k

d1_rhs_et_dx_0_1m2p1nyp2p0k = d1_rhs_et_dx_0_1m2p1nyp2p0k*param_float(1)

d1_rhs_et_dy_9_1m2p1nyp2p0p0k = (q(1-2+1,ny+2+0+0,indvars(1))*q(1-2+1,ny+2+0+0,indvars(4))+&
                    (param_float(3 + 5))*q(1-2+1,ny+2+0+0,indvars(1))*((q(1-2+1,ny+2+0+0,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2+0+0,indvars(2))*q(1-2+1,ny+2+0+0,indvars(2))+&
                    q(1-2+1,ny+2+0+0,indvars(3))*q(1-2+1,ny+2+0+0,indvars(3))))))*q(1-2+1,ny+2+0+0,indvars(3))

d1_rhs_et_dy_9_1m2p1nyp2p0m1k = (q(1-2+1,ny+2+0-1,indvars(1))*q(1-2+1,ny+2+0-1,indvars(4))+&
                    (param_float(3 + 5))*q(1-2+1,ny+2+0-1,indvars(1))*((q(1-2+1,ny+2+0-1,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2+0-1,indvars(2))*q(1-2+1,ny+2+0-1,indvars(2))+&
                    q(1-2+1,ny+2+0-1,indvars(3))*q(1-2+1,ny+2+0-1,indvars(3))))))*q(1-2+1,ny+2+0-1,indvars(3))

d1_rhs_et_dy_9_1m2p1nyp2p0m2k = (q(1-2+1,ny+2+0-2,indvars(1))*q(1-2+1,ny+2+0-2,indvars(4))+&
                    (param_float(3 + 5))*q(1-2+1,ny+2+0-2,indvars(1))*((q(1-2+1,ny+2+0-2,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2+0-2,indvars(2))*q(1-2+1,ny+2+0-2,indvars(2))+&
                    q(1-2+1,ny+2+0-2,indvars(3))*q(1-2+1,ny+2+0-2,indvars(3))))))*q(1-2+1,ny+2+0-2,indvars(3))

d1_rhs_et_dy_9_1m2p1nyp2p0k = 1.5_wp*d1_rhs_et_dy_9_1m2p1nyp2p0p0k-&
          2.0_wp*d1_rhs_et_dy_9_1m2p1nyp2p0m1k+&
          0.5_wp*d1_rhs_et_dy_9_1m2p1nyp2p0m2k

d1_rhs_et_dy_9_1m2p1nyp2p0k = d1_rhs_et_dy_9_1m2p1nyp2p0k*param_float(2)

d2_rhs_et_dydy_10_0_1m2p1nyp2p0p0k_1m2p1nyp2p0p0p0k = ((q(1-2+1,ny+2+0+0+0,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2+0+0+0,indvars(2))*q(1-2+1,ny+2+0+0+0,indvars(2))+&
                    q(1-2+1,ny+2+0+0+0,indvars(3))*q(1-2+1,ny+2+0+0+0,indvars(3)))))/param_float(4 + 5)

d2_rhs_et_dydy_10_0_1m2p1nyp2p0p0k_1m2p1nyp2p0p0m1k = ((q(1-2+1,ny+2+0+0-1,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2+0+0-1,indvars(2))*q(1-2+1,ny+2+0+0-1,indvars(2))+&
                    q(1-2+1,ny+2+0+0-1,indvars(3))*q(1-2+1,ny+2+0+0-1,indvars(3)))))/param_float(4 + 5)

d2_rhs_et_dydy_10_0_1m2p1nyp2p0p0k_1m2p1nyp2p0p0m2k = ((q(1-2+1,ny+2+0+0-2,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2+0+0-2,indvars(2))*q(1-2+1,ny+2+0+0-2,indvars(2))+&
                    q(1-2+1,ny+2+0+0-2,indvars(3))*q(1-2+1,ny+2+0+0-2,indvars(3)))))/param_float(4 + 5)

d2_rhs_et_dydy_10_0_1m2p1nyp2p0p0k = 1.5_wp*d2_rhs_et_dydy_10_0_1m2p1nyp2p0p0k_1m2p1nyp2p0p0p0k-&
          2.0_wp*d2_rhs_et_dydy_10_0_1m2p1nyp2p0p0k_1m2p1nyp2p0p0m1k+&
          0.5_wp*d2_rhs_et_dydy_10_0_1m2p1nyp2p0p0k_1m2p1nyp2p0p0m2k

d2_rhs_et_dydy_10_0_1m2p1nyp2p0p0k = d2_rhs_et_dydy_10_0_1m2p1nyp2p0p0k*param_float(2)

d2_rhs_et_dydy_10_0_1m2p1nyp2p0m1k_1m2p1nyp2p0m1m1k = ((q(1-2+1,ny+2+0-1-1,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2+0-1-1,indvars(2))*q(1-2+1,ny+2+0-1-1,indvars(2))+&
                    q(1-2+1,ny+2+0-1-1,indvars(3))*q(1-2+1,ny+2+0-1-1,indvars(3)))))/param_float(4 + 5)

d2_rhs_et_dydy_10_0_1m2p1nyp2p0m1k_1m2p1nyp2p0m1p1k = ((q(1-2+1,ny+2+0-1+1,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2+0-1+1,indvars(2))*q(1-2+1,ny+2+0-1+1,indvars(2))+&
                    q(1-2+1,ny+2+0-1+1,indvars(3))*q(1-2+1,ny+2+0-1+1,indvars(3)))))/param_float(4 + 5)

d2_rhs_et_dydy_10_0_1m2p1nyp2p0m1k = -&
          0.5_wp*d2_rhs_et_dydy_10_0_1m2p1nyp2p0m1k_1m2p1nyp2p0m1m1k+&
          0.5_wp*d2_rhs_et_dydy_10_0_1m2p1nyp2p0m1k_1m2p1nyp2p0m1p1k

d2_rhs_et_dydy_10_0_1m2p1nyp2p0m1k = d2_rhs_et_dydy_10_0_1m2p1nyp2p0m1k*param_float(2)

d2_rhs_et_dydy_10_0_1m2p1nyp2p0m2k_1m2p1nyp2p0m2m1k = ((q(1-2+1,ny+2+0-2-1,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2+0-2-1,indvars(2))*q(1-2+1,ny+2+0-2-1,indvars(2))+&
                    q(1-2+1,ny+2+0-2-1,indvars(3))*q(1-2+1,ny+2+0-2-1,indvars(3)))))/param_float(4 + 5)

d2_rhs_et_dydy_10_0_1m2p1nyp2p0m2k_1m2p1nyp2p0m2p1k = ((q(1-2+1,ny+2+0-2+1,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2+0-2+1,indvars(2))*q(1-2+1,ny+2+0-2+1,indvars(2))+&
                    q(1-2+1,ny+2+0-2+1,indvars(3))*q(1-2+1,ny+2+0-2+1,indvars(3)))))/param_float(4 + 5)

d2_rhs_et_dydy_10_0_1m2p1nyp2p0m2k = -&
          0.5_wp*d2_rhs_et_dydy_10_0_1m2p1nyp2p0m2k_1m2p1nyp2p0m2m1k+&
          0.5_wp*d2_rhs_et_dydy_10_0_1m2p1nyp2p0m2k_1m2p1nyp2p0m2p1k

d2_rhs_et_dydy_10_0_1m2p1nyp2p0m2k = d2_rhs_et_dydy_10_0_1m2p1nyp2p0m2k*param_float(2)

d1_rhs_et_dy_10_1m2p1nyp2p0p0k = -param_float(2 + 5)*qst(1-2+1,ny+2+0+0,indvarsst(11))*(d2_rhs_et_dydy_10_0_1m2p1nyp2p0p0k)-&
                    q(1-2+1,ny+2+0+0,indvars(2))*(((1+&
                    param_float(21 + 5))/(((q(1-2+1,ny+2+0+0,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2+0+0,indvars(2))*q(1-2+1,ny+2+0+0,indvars(2))+&
                    q(1-2+1,ny+2+0+0,indvars(3))*q(1-2+1,ny+2+0+0,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(1-2+1,ny+2+0+0,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2+0+0,indvars(2))*q(1-2+1,ny+2+0+0,indvars(2))+&
                    q(1-2+1,ny+2+0+0,indvars(3))*q(1-2+1,ny+2+0+0,indvars(3)))))/param_float(4 + 5)**1.5)*(1+&
                    ((q(1-2+1,ny+2+0+0,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(1-2+1,ny+2+0+0,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2+0+0,indvars(2))*q(1-2+1,ny+2+0+0,indvars(2))+&
                    q(1-2+1,ny+2+0+0,indvars(3))*q(1-2+1,ny+2+0+0,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(1-2+1,ny+2+0+0,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2+0+0,indvars(2))*q(1-2+1,ny+2+0+0,indvars(2))+&
                    q(1-2+1,ny+2+0+0,indvars(3))*q(1-2+1,ny+2+0+0,indvars(3)))))/param_float(4 + 5)**1.5*q(1-2+1,ny+2+0+0,indvars(1)))**3/((q(1-2+1,ny+2+0+0,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(1-2+1,ny+2+0+0,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2+0+0,indvars(2))*q(1-2+1,ny+2+0+0,indvars(2))+&
                    q(1-2+1,ny+2+0+0,indvars(3))*q(1-2+1,ny+2+0+0,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(1-2+1,ny+2+0+0,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2+0+0,indvars(2))*q(1-2+1,ny+2+0+0,indvars(2))+&
                    q(1-2+1,ny+2+0+0,indvars(3))*q(1-2+1,ny+2+0+0,indvars(3)))))/param_float(4 + 5)**1.5*q(1-2+1,ny+2+0+0,indvars(1)))**3+&
                    param_float(13 + 5)**3))*(q(1-2+1,ny+2+0+0,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(1-2+1,ny+2+0+0,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2+0+0,indvars(2))*q(1-2+1,ny+2+0+0,indvars(2))+&
                    q(1-2+1,ny+2+0+0,indvars(3))*q(1-2+1,ny+2+0+0,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(1-2+1,ny+2+0+0,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2+0+0,indvars(2))*q(1-2+1,ny+2+0+0,indvars(2))+&
                    q(1-2+1,ny+2+0+0,indvars(3))*q(1-2+1,ny+2+0+0,indvars(3)))))/param_float(4 + 5)**1.5*q(1-2+1,ny+2+0+0,indvars(1))))*param_float(1 + 5)*(qst(1-2+1,ny+2+0+0,indvarsst(11))*(d2_rhs_rhou_dydy_7_0_1m2p1nyp2p0p0k)+&
                    qst(1-2+1,ny+2+0+0,indvarsst(10))*(d2_rhs_rhou_dydx_7_0_1m2p1nyp2p0p0k)))-&
                    q(1-2+1,ny+2+0+0,indvars(3))*(((1+&
                    param_float(21 + 5))/(((q(1-2+1,ny+2+0+0,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2+0+0,indvars(2))*q(1-2+1,ny+2+0+0,indvars(2))+&
                    q(1-2+1,ny+2+0+0,indvars(3))*q(1-2+1,ny+2+0+0,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(1-2+1,ny+2+0+0,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2+0+0,indvars(2))*q(1-2+1,ny+2+0+0,indvars(2))+&
                    q(1-2+1,ny+2+0+0,indvars(3))*q(1-2+1,ny+2+0+0,indvars(3)))))/param_float(4 + 5)**1.5)*(1+&
                    ((q(1-2+1,ny+2+0+0,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(1-2+1,ny+2+0+0,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2+0+0,indvars(2))*q(1-2+1,ny+2+0+0,indvars(2))+&
                    q(1-2+1,ny+2+0+0,indvars(3))*q(1-2+1,ny+2+0+0,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(1-2+1,ny+2+0+0,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2+0+0,indvars(2))*q(1-2+1,ny+2+0+0,indvars(2))+&
                    q(1-2+1,ny+2+0+0,indvars(3))*q(1-2+1,ny+2+0+0,indvars(3)))))/param_float(4 + 5)**1.5*q(1-2+1,ny+2+0+0,indvars(1)))**3/((q(1-2+1,ny+2+0+0,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(1-2+1,ny+2+0+0,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2+0+0,indvars(2))*q(1-2+1,ny+2+0+0,indvars(2))+&
                    q(1-2+1,ny+2+0+0,indvars(3))*q(1-2+1,ny+2+0+0,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(1-2+1,ny+2+0+0,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2+0+0,indvars(2))*q(1-2+1,ny+2+0+0,indvars(2))+&
                    q(1-2+1,ny+2+0+0,indvars(3))*q(1-2+1,ny+2+0+0,indvars(3)))))/param_float(4 + 5)**1.5*q(1-2+1,ny+2+0+0,indvars(1)))**3+&
                    param_float(13 + 5)**3))*(q(1-2+1,ny+2+0+0,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(1-2+1,ny+2+0+0,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2+0+0,indvars(2))*q(1-2+1,ny+2+0+0,indvars(2))+&
                    q(1-2+1,ny+2+0+0,indvars(3))*q(1-2+1,ny+2+0+0,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(1-2+1,ny+2+0+0,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2+0+0,indvars(2))*q(1-2+1,ny+2+0+0,indvars(2))+&
                    q(1-2+1,ny+2+0+0,indvars(3))*q(1-2+1,ny+2+0+0,indvars(3)))))/param_float(4 + 5)**1.5*q(1-2+1,ny+2+0+0,indvars(1))))*param_float(1 + 5)*(2.0_wp*qst(1-2+1,ny+2+0+0,indvarsst(11))*(d2_rhs_rhov_dydy_6_0_1m2p1nyp2p0p0k)-&
                    2.0_wp/3.0_wp*(qst(1-2+1,ny+2+0+0,indvarsst(10))*(d2_rhs_rhov_dydx_6_0_1m2p1nyp2p0p0k)+&
                    qst(1-2+1,ny+2+0+0,indvarsst(11))*(d2_rhs_rhov_dydy_6_0_1m2p1nyp2p0p0k))))

d1_rhs_et_dy_10_1m2p1nyp2p0m1k = -param_float(2 + 5)*qst(1-2+1,ny+2+0-1,indvarsst(11))*(d2_rhs_et_dydy_10_0_1m2p1nyp2p0m1k)-&
                    q(1-2+1,ny+2+0-1,indvars(2))*(((1+&
                    param_float(21 + 5))/(((q(1-2+1,ny+2+0-1,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2+0-1,indvars(2))*q(1-2+1,ny+2+0-1,indvars(2))+&
                    q(1-2+1,ny+2+0-1,indvars(3))*q(1-2+1,ny+2+0-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(1-2+1,ny+2+0-1,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2+0-1,indvars(2))*q(1-2+1,ny+2+0-1,indvars(2))+&
                    q(1-2+1,ny+2+0-1,indvars(3))*q(1-2+1,ny+2+0-1,indvars(3)))))/param_float(4 + 5)**1.5)*(1+&
                    ((q(1-2+1,ny+2+0-1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(1-2+1,ny+2+0-1,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2+0-1,indvars(2))*q(1-2+1,ny+2+0-1,indvars(2))+&
                    q(1-2+1,ny+2+0-1,indvars(3))*q(1-2+1,ny+2+0-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(1-2+1,ny+2+0-1,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2+0-1,indvars(2))*q(1-2+1,ny+2+0-1,indvars(2))+&
                    q(1-2+1,ny+2+0-1,indvars(3))*q(1-2+1,ny+2+0-1,indvars(3)))))/param_float(4 + 5)**1.5*q(1-2+1,ny+2+0-1,indvars(1)))**3/((q(1-2+1,ny+2+0-1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(1-2+1,ny+2+0-1,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2+0-1,indvars(2))*q(1-2+1,ny+2+0-1,indvars(2))+&
                    q(1-2+1,ny+2+0-1,indvars(3))*q(1-2+1,ny+2+0-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(1-2+1,ny+2+0-1,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2+0-1,indvars(2))*q(1-2+1,ny+2+0-1,indvars(2))+&
                    q(1-2+1,ny+2+0-1,indvars(3))*q(1-2+1,ny+2+0-1,indvars(3)))))/param_float(4 + 5)**1.5*q(1-2+1,ny+2+0-1,indvars(1)))**3+&
                    param_float(13 + 5)**3))*(q(1-2+1,ny+2+0-1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(1-2+1,ny+2+0-1,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2+0-1,indvars(2))*q(1-2+1,ny+2+0-1,indvars(2))+&
                    q(1-2+1,ny+2+0-1,indvars(3))*q(1-2+1,ny+2+0-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(1-2+1,ny+2+0-1,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2+0-1,indvars(2))*q(1-2+1,ny+2+0-1,indvars(2))+&
                    q(1-2+1,ny+2+0-1,indvars(3))*q(1-2+1,ny+2+0-1,indvars(3)))))/param_float(4 + 5)**1.5*q(1-2+1,ny+2+0-1,indvars(1))))*param_float(1 + 5)*(qst(1-2+1,ny+2+0-1,indvarsst(11))*(d2_rhs_rhou_dydy_7_0_1m2p1nyp2p0m1k)+&
                    qst(1-2+1,ny+2+0-1,indvarsst(10))*(d2_rhs_rhou_dydx_7_0_1m2p1nyp2p0m1k)))-&
                    q(1-2+1,ny+2+0-1,indvars(3))*(((1+&
                    param_float(21 + 5))/(((q(1-2+1,ny+2+0-1,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2+0-1,indvars(2))*q(1-2+1,ny+2+0-1,indvars(2))+&
                    q(1-2+1,ny+2+0-1,indvars(3))*q(1-2+1,ny+2+0-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(1-2+1,ny+2+0-1,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2+0-1,indvars(2))*q(1-2+1,ny+2+0-1,indvars(2))+&
                    q(1-2+1,ny+2+0-1,indvars(3))*q(1-2+1,ny+2+0-1,indvars(3)))))/param_float(4 + 5)**1.5)*(1+&
                    ((q(1-2+1,ny+2+0-1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(1-2+1,ny+2+0-1,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2+0-1,indvars(2))*q(1-2+1,ny+2+0-1,indvars(2))+&
                    q(1-2+1,ny+2+0-1,indvars(3))*q(1-2+1,ny+2+0-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(1-2+1,ny+2+0-1,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2+0-1,indvars(2))*q(1-2+1,ny+2+0-1,indvars(2))+&
                    q(1-2+1,ny+2+0-1,indvars(3))*q(1-2+1,ny+2+0-1,indvars(3)))))/param_float(4 + 5)**1.5*q(1-2+1,ny+2+0-1,indvars(1)))**3/((q(1-2+1,ny+2+0-1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(1-2+1,ny+2+0-1,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2+0-1,indvars(2))*q(1-2+1,ny+2+0-1,indvars(2))+&
                    q(1-2+1,ny+2+0-1,indvars(3))*q(1-2+1,ny+2+0-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(1-2+1,ny+2+0-1,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2+0-1,indvars(2))*q(1-2+1,ny+2+0-1,indvars(2))+&
                    q(1-2+1,ny+2+0-1,indvars(3))*q(1-2+1,ny+2+0-1,indvars(3)))))/param_float(4 + 5)**1.5*q(1-2+1,ny+2+0-1,indvars(1)))**3+&
                    param_float(13 + 5)**3))*(q(1-2+1,ny+2+0-1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(1-2+1,ny+2+0-1,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2+0-1,indvars(2))*q(1-2+1,ny+2+0-1,indvars(2))+&
                    q(1-2+1,ny+2+0-1,indvars(3))*q(1-2+1,ny+2+0-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(1-2+1,ny+2+0-1,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2+0-1,indvars(2))*q(1-2+1,ny+2+0-1,indvars(2))+&
                    q(1-2+1,ny+2+0-1,indvars(3))*q(1-2+1,ny+2+0-1,indvars(3)))))/param_float(4 + 5)**1.5*q(1-2+1,ny+2+0-1,indvars(1))))*param_float(1 + 5)*(2.0_wp*qst(1-2+1,ny+2+0-1,indvarsst(11))*(d2_rhs_rhov_dydy_6_0_1m2p1nyp2p0m1k)-&
                    2.0_wp/3.0_wp*(qst(1-2+1,ny+2+0-1,indvarsst(10))*(d2_rhs_rhov_dydx_6_0_1m2p1nyp2p0m1k)+&
                    qst(1-2+1,ny+2+0-1,indvarsst(11))*(d2_rhs_rhov_dydy_6_0_1m2p1nyp2p0m1k))))

d1_rhs_et_dy_10_1m2p1nyp2p0m2k = -param_float(2 + 5)*qst(1-2+1,ny+2+0-2,indvarsst(11))*(d2_rhs_et_dydy_10_0_1m2p1nyp2p0m2k)-&
                    q(1-2+1,ny+2+0-2,indvars(2))*(((1+&
                    param_float(21 + 5))/(((q(1-2+1,ny+2+0-2,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2+0-2,indvars(2))*q(1-2+1,ny+2+0-2,indvars(2))+&
                    q(1-2+1,ny+2+0-2,indvars(3))*q(1-2+1,ny+2+0-2,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(1-2+1,ny+2+0-2,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2+0-2,indvars(2))*q(1-2+1,ny+2+0-2,indvars(2))+&
                    q(1-2+1,ny+2+0-2,indvars(3))*q(1-2+1,ny+2+0-2,indvars(3)))))/param_float(4 + 5)**1.5)*(1+&
                    ((q(1-2+1,ny+2+0-2,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(1-2+1,ny+2+0-2,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2+0-2,indvars(2))*q(1-2+1,ny+2+0-2,indvars(2))+&
                    q(1-2+1,ny+2+0-2,indvars(3))*q(1-2+1,ny+2+0-2,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(1-2+1,ny+2+0-2,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2+0-2,indvars(2))*q(1-2+1,ny+2+0-2,indvars(2))+&
                    q(1-2+1,ny+2+0-2,indvars(3))*q(1-2+1,ny+2+0-2,indvars(3)))))/param_float(4 + 5)**1.5*q(1-2+1,ny+2+0-2,indvars(1)))**3/((q(1-2+1,ny+2+0-2,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(1-2+1,ny+2+0-2,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2+0-2,indvars(2))*q(1-2+1,ny+2+0-2,indvars(2))+&
                    q(1-2+1,ny+2+0-2,indvars(3))*q(1-2+1,ny+2+0-2,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(1-2+1,ny+2+0-2,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2+0-2,indvars(2))*q(1-2+1,ny+2+0-2,indvars(2))+&
                    q(1-2+1,ny+2+0-2,indvars(3))*q(1-2+1,ny+2+0-2,indvars(3)))))/param_float(4 + 5)**1.5*q(1-2+1,ny+2+0-2,indvars(1)))**3+&
                    param_float(13 + 5)**3))*(q(1-2+1,ny+2+0-2,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(1-2+1,ny+2+0-2,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2+0-2,indvars(2))*q(1-2+1,ny+2+0-2,indvars(2))+&
                    q(1-2+1,ny+2+0-2,indvars(3))*q(1-2+1,ny+2+0-2,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(1-2+1,ny+2+0-2,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2+0-2,indvars(2))*q(1-2+1,ny+2+0-2,indvars(2))+&
                    q(1-2+1,ny+2+0-2,indvars(3))*q(1-2+1,ny+2+0-2,indvars(3)))))/param_float(4 + 5)**1.5*q(1-2+1,ny+2+0-2,indvars(1))))*param_float(1 + 5)*(qst(1-2+1,ny+2+0-2,indvarsst(11))*(d2_rhs_rhou_dydy_7_0_1m2p1nyp2p0m2k)+&
                    qst(1-2+1,ny+2+0-2,indvarsst(10))*(d2_rhs_rhou_dydx_7_0_1m2p1nyp2p0m2k)))-&
                    q(1-2+1,ny+2+0-2,indvars(3))*(((1+&
                    param_float(21 + 5))/(((q(1-2+1,ny+2+0-2,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2+0-2,indvars(2))*q(1-2+1,ny+2+0-2,indvars(2))+&
                    q(1-2+1,ny+2+0-2,indvars(3))*q(1-2+1,ny+2+0-2,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(1-2+1,ny+2+0-2,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2+0-2,indvars(2))*q(1-2+1,ny+2+0-2,indvars(2))+&
                    q(1-2+1,ny+2+0-2,indvars(3))*q(1-2+1,ny+2+0-2,indvars(3)))))/param_float(4 + 5)**1.5)*(1+&
                    ((q(1-2+1,ny+2+0-2,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(1-2+1,ny+2+0-2,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2+0-2,indvars(2))*q(1-2+1,ny+2+0-2,indvars(2))+&
                    q(1-2+1,ny+2+0-2,indvars(3))*q(1-2+1,ny+2+0-2,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(1-2+1,ny+2+0-2,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2+0-2,indvars(2))*q(1-2+1,ny+2+0-2,indvars(2))+&
                    q(1-2+1,ny+2+0-2,indvars(3))*q(1-2+1,ny+2+0-2,indvars(3)))))/param_float(4 + 5)**1.5*q(1-2+1,ny+2+0-2,indvars(1)))**3/((q(1-2+1,ny+2+0-2,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(1-2+1,ny+2+0-2,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2+0-2,indvars(2))*q(1-2+1,ny+2+0-2,indvars(2))+&
                    q(1-2+1,ny+2+0-2,indvars(3))*q(1-2+1,ny+2+0-2,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(1-2+1,ny+2+0-2,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2+0-2,indvars(2))*q(1-2+1,ny+2+0-2,indvars(2))+&
                    q(1-2+1,ny+2+0-2,indvars(3))*q(1-2+1,ny+2+0-2,indvars(3)))))/param_float(4 + 5)**1.5*q(1-2+1,ny+2+0-2,indvars(1)))**3+&
                    param_float(13 + 5)**3))*(q(1-2+1,ny+2+0-2,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(1-2+1,ny+2+0-2,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2+0-2,indvars(2))*q(1-2+1,ny+2+0-2,indvars(2))+&
                    q(1-2+1,ny+2+0-2,indvars(3))*q(1-2+1,ny+2+0-2,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(1-2+1,ny+2+0-2,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2+0-2,indvars(2))*q(1-2+1,ny+2+0-2,indvars(2))+&
                    q(1-2+1,ny+2+0-2,indvars(3))*q(1-2+1,ny+2+0-2,indvars(3)))))/param_float(4 + 5)**1.5*q(1-2+1,ny+2+0-2,indvars(1))))*param_float(1 + 5)*(2.0_wp*qst(1-2+1,ny+2+0-2,indvarsst(11))*(d2_rhs_rhov_dydy_6_0_1m2p1nyp2p0m2k)-&
                    2.0_wp/3.0_wp*(qst(1-2+1,ny+2+0-2,indvarsst(10))*(d2_rhs_rhov_dydx_6_0_1m2p1nyp2p0m2k)+&
                    qst(1-2+1,ny+2+0-2,indvarsst(11))*(d2_rhs_rhov_dydy_6_0_1m2p1nyp2p0m2k))))

d1_rhs_et_dy_10_1m2p1nyp2p0k = 1.5_wp*d1_rhs_et_dy_10_1m2p1nyp2p0p0k-&
          2.0_wp*d1_rhs_et_dy_10_1m2p1nyp2p0m1k+&
          0.5_wp*d1_rhs_et_dy_10_1m2p1nyp2p0m2k

d1_rhs_et_dy_10_1m2p1nyp2p0k = d1_rhs_et_dy_10_1m2p1nyp2p0k*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer 1 0 None d(rho et)/dt **********
!                                                           
!***********************************************************


rhs(1-2+1,ny+2+0,indvars(4)) =   -  ( 0.5_wp*(q(1-2+1,ny+2+0,indvars(2))**2+&
                    q(1-2+1,ny+2+0,indvars(3))**2)*(1.0_wp/(param_float(23 + 5)*(param_float(3 + 5))*q(1-2+1,ny+2+0,indvars(1))*((q(1-2+1,ny+2+0,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2+0,indvars(2))*q(1-2+1,ny+2+0,indvars(2))+&
                    q(1-2+1,ny+2+0,indvars(3))*q(1-2+1,ny+2+0,indvars(3)))))/q(1-2+1,ny+2+0,indvars(1)))**2*((q(1-2+1,ny+2+0,indvars(3))*(-&
                    (param_float(23 + 5)*(param_float(3 + 5))*q(1-2+1,ny+2+0,indvars(1))*((q(1-2+1,ny+2+0,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2+0,indvars(2))*q(1-2+1,ny+2+0,indvars(2))+&
                    q(1-2+1,ny+2+0,indvars(3))*q(1-2+1,ny+2+0,indvars(3)))))/q(1-2+1,ny+2+0,indvars(1)))**2*d1_rhs_rho_dy_1_1m2p1nyp2p0k+&
                    1.0_wp*d1_rhs_rho_dy_2_1m2p1nyp2p0k))*qst(1-2+1,ny+2+0,indvarsst(11))+&
                    0.5_wp*((((param_float(23 + 5)*(param_float(3 + 5))*q(1-2+1,ny+2+0,indvars(1))*((q(1-2+1,ny+2+0,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2+0,indvars(2))*q(1-2+1,ny+2+0,indvars(2))+&
                    q(1-2+1,ny+2+0,indvars(3))*q(1-2+1,ny+2+0,indvars(3)))))/q(1-2+1,ny+2+0,indvars(1)))+&
                    q(1-2+1,ny+2+0,indvars(3)))*((param_float(23 + 5)*(param_float(3 + 5))*q(1-2+1,ny+2+0,indvars(1))*((q(1-2+1,ny+2+0,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2+0,indvars(2))*q(1-2+1,ny+2+0,indvars(2))+&
                    q(1-2+1,ny+2+0,indvars(3))*q(1-2+1,ny+2+0,indvars(3)))))/q(1-2+1,ny+2+0,indvars(1)))*q(1-2+1,ny+2+0,indvars(1))*d1_rhs_rho_dy_3_1m2p1nyp2p0k+&
                    1.0_wp*d1_rhs_rho_dy_2_1m2p1nyp2p0k))*qst(1-2+1,ny+2+0,indvarsst(11))+&
                    param_float(19 + 5)*(param_float(23 + 5)*(param_float(3 + 5))*q(1-2+1,ny+2+0,indvars(1))*((q(1-2+1,ny+2+0,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2+0,indvars(2))*q(1-2+1,ny+2+0,indvars(2))+&
                    q(1-2+1,ny+2+0,indvars(3))*q(1-2+1,ny+2+0,indvars(3)))))/q(1-2+1,ny+2+0,indvars(1)))*(1-&
                    qface_j(1-2+1,1)*qface_j(1-2+1,1))/param_float(20 + 5)*((param_float(3 + 5))*q(1-2+1,ny+2+0,indvars(1))*((q(1-2+1,ny+2+0,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2+0,indvars(2))*q(1-2+1,ny+2+0,indvars(2))+&
                    q(1-2+1,ny+2+0,indvars(3))*q(1-2+1,ny+2+0,indvars(3)))))-&
                    param_float(26 + 5)))))+&
                    1.0_wp/param_float(3 + 5)*(0.5_wp*((((param_float(23 + 5)*(param_float(3 + 5))*q(1-2+1,ny+2+0,indvars(1))*((q(1-2+1,ny+2+0,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2+0,indvars(2))*q(1-2+1,ny+2+0,indvars(2))+&
                    q(1-2+1,ny+2+0,indvars(3))*q(1-2+1,ny+2+0,indvars(3)))))/q(1-2+1,ny+2+0,indvars(1)))+&
                    q(1-2+1,ny+2+0,indvars(3)))*((param_float(23 + 5)*(param_float(3 + 5))*q(1-2+1,ny+2+0,indvars(1))*((q(1-2+1,ny+2+0,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2+0,indvars(2))*q(1-2+1,ny+2+0,indvars(2))+&
                    q(1-2+1,ny+2+0,indvars(3))*q(1-2+1,ny+2+0,indvars(3)))))/q(1-2+1,ny+2+0,indvars(1)))*q(1-2+1,ny+2+0,indvars(1))*d1_rhs_rho_dy_3_1m2p1nyp2p0k+&
                    1.0_wp*d1_rhs_rho_dy_2_1m2p1nyp2p0k))*qst(1-2+1,ny+2+0,indvarsst(11))+&
                    param_float(19 + 5)*(param_float(23 + 5)*(param_float(3 + 5))*q(1-2+1,ny+2+0,indvars(1))*((q(1-2+1,ny+2+0,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2+0,indvars(2))*q(1-2+1,ny+2+0,indvars(2))+&
                    q(1-2+1,ny+2+0,indvars(3))*q(1-2+1,ny+2+0,indvars(3)))))/q(1-2+1,ny+2+0,indvars(1)))*(1-&
                    qface_j(1-2+1,1)*qface_j(1-2+1,1))/param_float(20 + 5)*((param_float(3 + 5))*q(1-2+1,ny+2+0,indvars(1))*((q(1-2+1,ny+2+0,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2+0,indvars(2))*q(1-2+1,ny+2+0,indvars(2))+&
                    q(1-2+1,ny+2+0,indvars(3))*q(1-2+1,ny+2+0,indvars(3)))))-&
                    param_float(26 + 5))))+&
                    q(1-2+1,ny+2+0,indvars(1))*q(1-2+1,ny+2+0,indvars(2))*(1.0_wp/(2*q(1-2+1,ny+2+0,indvars(1))*(param_float(23 + 5)*(param_float(3 + 5))*q(1-2+1,ny+2+0,indvars(1))*((q(1-2+1,ny+2+0,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2+0,indvars(2))*q(1-2+1,ny+2+0,indvars(2))+&
                    q(1-2+1,ny+2+0,indvars(3))*q(1-2+1,ny+2+0,indvars(3)))))/q(1-2+1,ny+2+0,indvars(1))))*((((param_float(23 + 5)*(param_float(3 + 5))*q(1-2+1,ny+2+0,indvars(1))*((q(1-2+1,ny+2+0,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2+0,indvars(2))*q(1-2+1,ny+2+0,indvars(2))+&
                    q(1-2+1,ny+2+0,indvars(3))*q(1-2+1,ny+2+0,indvars(3)))))/q(1-2+1,ny+2+0,indvars(1)))+&
                    q(1-2+1,ny+2+0,indvars(3)))*((param_float(23 + 5)*(param_float(3 + 5))*q(1-2+1,ny+2+0,indvars(1))*((q(1-2+1,ny+2+0,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2+0,indvars(2))*q(1-2+1,ny+2+0,indvars(2))+&
                    q(1-2+1,ny+2+0,indvars(3))*q(1-2+1,ny+2+0,indvars(3)))))/q(1-2+1,ny+2+0,indvars(1)))*q(1-2+1,ny+2+0,indvars(1))*d1_rhs_rho_dy_3_1m2p1nyp2p0k+&
                    1.0_wp*d1_rhs_rho_dy_2_1m2p1nyp2p0k))*qst(1-2+1,ny+2+0,indvarsst(11))-&
                    param_float(19 + 5)*(param_float(23 + 5)*(param_float(3 + 5))*q(1-2+1,ny+2+0,indvars(1))*((q(1-2+1,ny+2+0,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2+0,indvars(2))*q(1-2+1,ny+2+0,indvars(2))+&
                    q(1-2+1,ny+2+0,indvars(3))*q(1-2+1,ny+2+0,indvars(3)))))/q(1-2+1,ny+2+0,indvars(1)))*(1-&
                    qface_j(1-2+1,1)*qface_j(1-2+1,1))/param_float(20 + 5)*((param_float(3 + 5))*q(1-2+1,ny+2+0,indvars(1))*((q(1-2+1,ny+2+0,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2+0,indvars(2))*q(1-2+1,ny+2+0,indvars(2))+&
                    q(1-2+1,ny+2+0,indvars(3))*q(1-2+1,ny+2+0,indvars(3)))))-&
                    param_float(26 + 5))))+&
                    q(1-2+1,ny+2+0,indvars(1))*q(1-2+1,ny+2+0,indvars(3))*((1.0_wp*q(1-2+1,ny+2+0,indvars(3))*d1_rhs_rhov_dy_4_1m2p1nyp2p0k)*qst(1-2+1,ny+2+0,indvarsst(11)))+&
                    d1_rhs_et_dy_9_1m2p1nyp2p0k+&
                    qst(1-2+1,ny+2+0,indvarsst(10))*(d1_rhs_et_dx_0_1m2p1nyp2p0k)+&
                    qst(1-2+1,ny+2+0,indvarsst(11))*(d1_rhs_et_dy_10_1m2p1nyp2p0k) ) 

