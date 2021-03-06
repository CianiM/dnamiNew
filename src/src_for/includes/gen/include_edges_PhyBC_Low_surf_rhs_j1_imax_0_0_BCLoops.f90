

!***********************************************************
!                                                           
! Start building layers for BC : imax j1 None **************
!                                                           
!***********************************************************




!***********************************************************
!                                                           
! BC layer: 0 0 None ***************************************
!                                                           
!***********************************************************




!***********************************************************
!                                                           
! building source terms in RHS for layer 0 0 None d(rho)/dt 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (deltaxI*([rho*u]_1x))*symm+(rho*([v]_1y)*deltayI)*wall
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_rhs_rho_dx_0_nxp2p0p01m2p0k = q(nx+2+0+0,1-2+0,indvars(1))*q(nx+2+0+0,1-2+0,indvars(2))

d1_rhs_rho_dx_0_nxp2p0m11m2p0k = q(nx+2+0-1,1-2+0,indvars(1))*q(nx+2+0-1,1-2+0,indvars(2))

d1_rhs_rho_dx_0_nxp2p0m21m2p0k = q(nx+2+0-2,1-2+0,indvars(1))*q(nx+2+0-2,1-2+0,indvars(2))

d1_rhs_rho_dx_0_nxp2p01m2p0k = 1.5_wp*d1_rhs_rho_dx_0_nxp2p0p01m2p0k-&
          2.0_wp*d1_rhs_rho_dx_0_nxp2p0m11m2p0k+&
          0.5_wp*d1_rhs_rho_dx_0_nxp2p0m21m2p0k

d1_rhs_rho_dx_0_nxp2p01m2p0k = d1_rhs_rho_dx_0_nxp2p01m2p0k*param_float(1)

d1_rhs_rho_dy_0_nxp2p01m2p0p0k = q(nx+2+0,1-2+0+0,indvars(3))

d1_rhs_rho_dy_0_nxp2p01m2p0p1k = q(nx+2+0,1-2+0+1,indvars(3))

d1_rhs_rho_dy_0_nxp2p01m2p0p2k = q(nx+2+0,1-2+0+2,indvars(3))

d1_rhs_rho_dy_0_nxp2p01m2p0k = -&
          1.5_wp*d1_rhs_rho_dy_0_nxp2p01m2p0p0k+&
          2.0_wp*d1_rhs_rho_dy_0_nxp2p01m2p0p1k-&
          0.5_wp*d1_rhs_rho_dy_0_nxp2p01m2p0p2k

d1_rhs_rho_dy_0_nxp2p01m2p0k = d1_rhs_rho_dy_0_nxp2p01m2p0k*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer 0 0 None d(rho)/dt *************
!                                                           
!***********************************************************


rhs(nx+2+0,1-2+0,indvars(1)) =   -  ( (qst(nx+2+0,1-2+0,indvarsst(10))*(d1_rhs_rho_dx_0_nxp2p01m2p0k))*qst(nx+2+0,1-2+0,indvarsst(5))+&
                    (q(nx+2+0,1-2+0,indvars(1))*(d1_rhs_rho_dy_0_nxp2p01m2p0k)*qst(nx+2+0,1-2+0,indvarsst(11)))*qst(nx+2+0,1-2+0,indvarsst(19)) ) 



!***********************************************************
!                                                           
! building source terms in RHS for layer 0 0 None d(rho u)/dt 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (deltaxI*([-4.0_wp/3.0_wp*visc_t*({u}_1x)*deltaxI]_1x)+deltaxI*([rho*u*u+p]_1x))*symm
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d2_rhs_u_dxdx_0_0_nxp2p0p01m2p0k_nxp2p0p0p01m2p0k = q(nx+2+0+0+0,1-2+0,indvars(2))

d2_rhs_u_dxdx_0_0_nxp2p0p01m2p0k_nxp2p0p0m11m2p0k = q(nx+2+0+0-1,1-2+0,indvars(2))

d2_rhs_u_dxdx_0_0_nxp2p0p01m2p0k_nxp2p0p0m21m2p0k = q(nx+2+0+0-2,1-2+0,indvars(2))

d2_rhs_u_dxdx_0_0_nxp2p0p01m2p0k = 1.5_wp*d2_rhs_u_dxdx_0_0_nxp2p0p01m2p0k_nxp2p0p0p01m2p0k-&
          2.0_wp*d2_rhs_u_dxdx_0_0_nxp2p0p01m2p0k_nxp2p0p0m11m2p0k+&
          0.5_wp*d2_rhs_u_dxdx_0_0_nxp2p0p01m2p0k_nxp2p0p0m21m2p0k

d2_rhs_u_dxdx_0_0_nxp2p0p01m2p0k = d2_rhs_u_dxdx_0_0_nxp2p0p01m2p0k*param_float(1)

d2_rhs_u_dxdx_0_0_nxp2p0m11m2p0k_nxp2p0m1m11m2p0k = q(nx+2+0-1-1,1-2+0,indvars(2))

d2_rhs_u_dxdx_0_0_nxp2p0m11m2p0k_nxp2p0m1p11m2p0k = q(nx+2+0-1+1,1-2+0,indvars(2))

d2_rhs_u_dxdx_0_0_nxp2p0m11m2p0k = -&
          0.5_wp*d2_rhs_u_dxdx_0_0_nxp2p0m11m2p0k_nxp2p0m1m11m2p0k+&
          0.5_wp*d2_rhs_u_dxdx_0_0_nxp2p0m11m2p0k_nxp2p0m1p11m2p0k

d2_rhs_u_dxdx_0_0_nxp2p0m11m2p0k = d2_rhs_u_dxdx_0_0_nxp2p0m11m2p0k*param_float(1)

d2_rhs_u_dxdx_0_0_nxp2p0m21m2p0k_nxp2p0m2m11m2p0k = q(nx+2+0-2-1,1-2+0,indvars(2))

d2_rhs_u_dxdx_0_0_nxp2p0m21m2p0k_nxp2p0m2p11m2p0k = q(nx+2+0-2+1,1-2+0,indvars(2))

d2_rhs_u_dxdx_0_0_nxp2p0m21m2p0k = -&
          0.5_wp*d2_rhs_u_dxdx_0_0_nxp2p0m21m2p0k_nxp2p0m2m11m2p0k+&
          0.5_wp*d2_rhs_u_dxdx_0_0_nxp2p0m21m2p0k_nxp2p0m2p11m2p0k

d2_rhs_u_dxdx_0_0_nxp2p0m21m2p0k = d2_rhs_u_dxdx_0_0_nxp2p0m21m2p0k*param_float(1)

d1_rhs_u_dx_1_nxp2p0p01m2p0k = q(nx+2+0+0,1-2+0,indvars(1))*q(nx+2+0+0,1-2+0,indvars(2))*q(nx+2+0+0,1-2+0,indvars(2))+(param_float(3 + 5))*q(nx+2+0+0,1-2+0,indvars(1))*((q(nx+2+0+0,1-2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0+0,1-2+0,indvars(2))*q(nx+2+0+0,1-2+0,indvars(2))+&
                    q(nx+2+0+0,1-2+0,indvars(3))*q(nx+2+0+0,1-2+0,indvars(3)))))

d1_rhs_u_dx_1_nxp2p0m11m2p0k = q(nx+2+0-1,1-2+0,indvars(1))*q(nx+2+0-1,1-2+0,indvars(2))*q(nx+2+0-1,1-2+0,indvars(2))+(param_float(3 + 5))*q(nx+2+0-1,1-2+0,indvars(1))*((q(nx+2+0-1,1-2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0-1,1-2+0,indvars(2))*q(nx+2+0-1,1-2+0,indvars(2))+&
                    q(nx+2+0-1,1-2+0,indvars(3))*q(nx+2+0-1,1-2+0,indvars(3)))))

d1_rhs_u_dx_1_nxp2p0m21m2p0k = q(nx+2+0-2,1-2+0,indvars(1))*q(nx+2+0-2,1-2+0,indvars(2))*q(nx+2+0-2,1-2+0,indvars(2))+(param_float(3 + 5))*q(nx+2+0-2,1-2+0,indvars(1))*((q(nx+2+0-2,1-2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0-2,1-2+0,indvars(2))*q(nx+2+0-2,1-2+0,indvars(2))+&
                    q(nx+2+0-2,1-2+0,indvars(3))*q(nx+2+0-2,1-2+0,indvars(3)))))

d1_rhs_u_dx_1_nxp2p01m2p0k = 1.5_wp*d1_rhs_u_dx_1_nxp2p0p01m2p0k-&
          2.0_wp*d1_rhs_u_dx_1_nxp2p0m11m2p0k+&
          0.5_wp*d1_rhs_u_dx_1_nxp2p0m21m2p0k

d1_rhs_u_dx_1_nxp2p01m2p0k = d1_rhs_u_dx_1_nxp2p01m2p0k*param_float(1)

d1_rhs_u_dx_0_nxp2p0p01m2p0k = -4.0_wp/3.0_wp*((1+&
                    param_float(21 + 5))/(((q(nx+2+0+0,1-2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0+0,1-2+0,indvars(2))*q(nx+2+0+0,1-2+0,indvars(2))+&
                    q(nx+2+0+0,1-2+0,indvars(3))*q(nx+2+0+0,1-2+0,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2+0+0,1-2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0+0,1-2+0,indvars(2))*q(nx+2+0+0,1-2+0,indvars(2))+&
                    q(nx+2+0+0,1-2+0,indvars(3))*q(nx+2+0+0,1-2+0,indvars(3)))))/param_float(4 + 5)**1.5)*(1+&
                    ((q(nx+2+0+0,1-2+0,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+2+0+0,1-2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0+0,1-2+0,indvars(2))*q(nx+2+0+0,1-2+0,indvars(2))+&
                    q(nx+2+0+0,1-2+0,indvars(3))*q(nx+2+0+0,1-2+0,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2+0+0,1-2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0+0,1-2+0,indvars(2))*q(nx+2+0+0,1-2+0,indvars(2))+&
                    q(nx+2+0+0,1-2+0,indvars(3))*q(nx+2+0+0,1-2+0,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+2+0+0,1-2+0,indvars(1)))**3/((q(nx+2+0+0,1-2+0,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+2+0+0,1-2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0+0,1-2+0,indvars(2))*q(nx+2+0+0,1-2+0,indvars(2))+&
                    q(nx+2+0+0,1-2+0,indvars(3))*q(nx+2+0+0,1-2+0,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2+0+0,1-2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0+0,1-2+0,indvars(2))*q(nx+2+0+0,1-2+0,indvars(2))+&
                    q(nx+2+0+0,1-2+0,indvars(3))*q(nx+2+0+0,1-2+0,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+2+0+0,1-2+0,indvars(1)))**3+&
                    param_float(13 + 5)**3))*(q(nx+2+0+0,1-2+0,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+2+0+0,1-2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0+0,1-2+0,indvars(2))*q(nx+2+0+0,1-2+0,indvars(2))+&
                    q(nx+2+0+0,1-2+0,indvars(3))*q(nx+2+0+0,1-2+0,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2+0+0,1-2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0+0,1-2+0,indvars(2))*q(nx+2+0+0,1-2+0,indvars(2))+&
                    q(nx+2+0+0,1-2+0,indvars(3))*q(nx+2+0+0,1-2+0,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+2+0+0,1-2+0,indvars(1))))*param_float(1 + 5)*(d2_rhs_u_dxdx_0_0_nxp2p0p01m2p0k)*qst(nx+2+0+0,1-2+0,indvarsst(10))

d1_rhs_u_dx_0_nxp2p0m11m2p0k = -4.0_wp/3.0_wp*((1+&
                    param_float(21 + 5))/(((q(nx+2+0-1,1-2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0-1,1-2+0,indvars(2))*q(nx+2+0-1,1-2+0,indvars(2))+&
                    q(nx+2+0-1,1-2+0,indvars(3))*q(nx+2+0-1,1-2+0,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2+0-1,1-2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0-1,1-2+0,indvars(2))*q(nx+2+0-1,1-2+0,indvars(2))+&
                    q(nx+2+0-1,1-2+0,indvars(3))*q(nx+2+0-1,1-2+0,indvars(3)))))/param_float(4 + 5)**1.5)*(1+&
                    ((q(nx+2+0-1,1-2+0,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+2+0-1,1-2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0-1,1-2+0,indvars(2))*q(nx+2+0-1,1-2+0,indvars(2))+&
                    q(nx+2+0-1,1-2+0,indvars(3))*q(nx+2+0-1,1-2+0,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2+0-1,1-2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0-1,1-2+0,indvars(2))*q(nx+2+0-1,1-2+0,indvars(2))+&
                    q(nx+2+0-1,1-2+0,indvars(3))*q(nx+2+0-1,1-2+0,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+2+0-1,1-2+0,indvars(1)))**3/((q(nx+2+0-1,1-2+0,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+2+0-1,1-2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0-1,1-2+0,indvars(2))*q(nx+2+0-1,1-2+0,indvars(2))+&
                    q(nx+2+0-1,1-2+0,indvars(3))*q(nx+2+0-1,1-2+0,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2+0-1,1-2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0-1,1-2+0,indvars(2))*q(nx+2+0-1,1-2+0,indvars(2))+&
                    q(nx+2+0-1,1-2+0,indvars(3))*q(nx+2+0-1,1-2+0,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+2+0-1,1-2+0,indvars(1)))**3+&
                    param_float(13 + 5)**3))*(q(nx+2+0-1,1-2+0,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+2+0-1,1-2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0-1,1-2+0,indvars(2))*q(nx+2+0-1,1-2+0,indvars(2))+&
                    q(nx+2+0-1,1-2+0,indvars(3))*q(nx+2+0-1,1-2+0,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2+0-1,1-2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0-1,1-2+0,indvars(2))*q(nx+2+0-1,1-2+0,indvars(2))+&
                    q(nx+2+0-1,1-2+0,indvars(3))*q(nx+2+0-1,1-2+0,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+2+0-1,1-2+0,indvars(1))))*param_float(1 + 5)*(d2_rhs_u_dxdx_0_0_nxp2p0m11m2p0k)*qst(nx+2+0-1,1-2+0,indvarsst(10))

d1_rhs_u_dx_0_nxp2p0m21m2p0k = -4.0_wp/3.0_wp*((1+&
                    param_float(21 + 5))/(((q(nx+2+0-2,1-2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0-2,1-2+0,indvars(2))*q(nx+2+0-2,1-2+0,indvars(2))+&
                    q(nx+2+0-2,1-2+0,indvars(3))*q(nx+2+0-2,1-2+0,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2+0-2,1-2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0-2,1-2+0,indvars(2))*q(nx+2+0-2,1-2+0,indvars(2))+&
                    q(nx+2+0-2,1-2+0,indvars(3))*q(nx+2+0-2,1-2+0,indvars(3)))))/param_float(4 + 5)**1.5)*(1+&
                    ((q(nx+2+0-2,1-2+0,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+2+0-2,1-2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0-2,1-2+0,indvars(2))*q(nx+2+0-2,1-2+0,indvars(2))+&
                    q(nx+2+0-2,1-2+0,indvars(3))*q(nx+2+0-2,1-2+0,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2+0-2,1-2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0-2,1-2+0,indvars(2))*q(nx+2+0-2,1-2+0,indvars(2))+&
                    q(nx+2+0-2,1-2+0,indvars(3))*q(nx+2+0-2,1-2+0,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+2+0-2,1-2+0,indvars(1)))**3/((q(nx+2+0-2,1-2+0,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+2+0-2,1-2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0-2,1-2+0,indvars(2))*q(nx+2+0-2,1-2+0,indvars(2))+&
                    q(nx+2+0-2,1-2+0,indvars(3))*q(nx+2+0-2,1-2+0,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2+0-2,1-2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0-2,1-2+0,indvars(2))*q(nx+2+0-2,1-2+0,indvars(2))+&
                    q(nx+2+0-2,1-2+0,indvars(3))*q(nx+2+0-2,1-2+0,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+2+0-2,1-2+0,indvars(1)))**3+&
                    param_float(13 + 5)**3))*(q(nx+2+0-2,1-2+0,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+2+0-2,1-2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0-2,1-2+0,indvars(2))*q(nx+2+0-2,1-2+0,indvars(2))+&
                    q(nx+2+0-2,1-2+0,indvars(3))*q(nx+2+0-2,1-2+0,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2+0-2,1-2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0-2,1-2+0,indvars(2))*q(nx+2+0-2,1-2+0,indvars(2))+&
                    q(nx+2+0-2,1-2+0,indvars(3))*q(nx+2+0-2,1-2+0,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+2+0-2,1-2+0,indvars(1))))*param_float(1 + 5)*(d2_rhs_u_dxdx_0_0_nxp2p0m21m2p0k)*qst(nx+2+0-2,1-2+0,indvarsst(10))

d1_rhs_u_dx_0_nxp2p01m2p0k = 1.5_wp*d1_rhs_u_dx_0_nxp2p0p01m2p0k-&
          2.0_wp*d1_rhs_u_dx_0_nxp2p0m11m2p0k+&
          0.5_wp*d1_rhs_u_dx_0_nxp2p0m21m2p0k

d1_rhs_u_dx_0_nxp2p01m2p0k = d1_rhs_u_dx_0_nxp2p01m2p0k*param_float(1)



!***********************************************************
!                                                           
! Update BC terms for layer 0 0 None d(rho u)/dt ***********
!                                                           
!***********************************************************


rhs(nx+2+0,1-2+0,indvars(2)) =   -  ( (qst(nx+2+0,1-2+0,indvarsst(10))*(d1_rhs_u_dx_0_nxp2p01m2p0k)+&
                    qst(nx+2+0,1-2+0,indvarsst(10))*(d1_rhs_u_dx_1_nxp2p01m2p0k))*qst(nx+2+0,1-2+0,indvarsst(5)) ) 



!***********************************************************
!                                                           
! building source terms in RHS for layer 0 0 None d(rho v)/dt 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (deltaxI*([-visc_t*({v}_1x)*deltaxI]_1x)+deltaxI*([rho*u*v]_1x))*symm
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d2_rhs_v_dxdx_0_0_nxp2p0p01m2p0k_nxp2p0p0p01m2p0k = q(nx+2+0+0+0,1-2+0,indvars(3))

d2_rhs_v_dxdx_0_0_nxp2p0p01m2p0k_nxp2p0p0m11m2p0k = q(nx+2+0+0-1,1-2+0,indvars(3))

d2_rhs_v_dxdx_0_0_nxp2p0p01m2p0k_nxp2p0p0m21m2p0k = q(nx+2+0+0-2,1-2+0,indvars(3))

d2_rhs_v_dxdx_0_0_nxp2p0p01m2p0k = 1.5_wp*d2_rhs_v_dxdx_0_0_nxp2p0p01m2p0k_nxp2p0p0p01m2p0k-&
          2.0_wp*d2_rhs_v_dxdx_0_0_nxp2p0p01m2p0k_nxp2p0p0m11m2p0k+&
          0.5_wp*d2_rhs_v_dxdx_0_0_nxp2p0p01m2p0k_nxp2p0p0m21m2p0k

d2_rhs_v_dxdx_0_0_nxp2p0p01m2p0k = d2_rhs_v_dxdx_0_0_nxp2p0p01m2p0k*param_float(1)

d2_rhs_v_dxdx_0_0_nxp2p0m11m2p0k_nxp2p0m1m11m2p0k = q(nx+2+0-1-1,1-2+0,indvars(3))

d2_rhs_v_dxdx_0_0_nxp2p0m11m2p0k_nxp2p0m1p11m2p0k = q(nx+2+0-1+1,1-2+0,indvars(3))

d2_rhs_v_dxdx_0_0_nxp2p0m11m2p0k = -&
          0.5_wp*d2_rhs_v_dxdx_0_0_nxp2p0m11m2p0k_nxp2p0m1m11m2p0k+&
          0.5_wp*d2_rhs_v_dxdx_0_0_nxp2p0m11m2p0k_nxp2p0m1p11m2p0k

d2_rhs_v_dxdx_0_0_nxp2p0m11m2p0k = d2_rhs_v_dxdx_0_0_nxp2p0m11m2p0k*param_float(1)

d2_rhs_v_dxdx_0_0_nxp2p0m21m2p0k_nxp2p0m2m11m2p0k = q(nx+2+0-2-1,1-2+0,indvars(3))

d2_rhs_v_dxdx_0_0_nxp2p0m21m2p0k_nxp2p0m2p11m2p0k = q(nx+2+0-2+1,1-2+0,indvars(3))

d2_rhs_v_dxdx_0_0_nxp2p0m21m2p0k = -&
          0.5_wp*d2_rhs_v_dxdx_0_0_nxp2p0m21m2p0k_nxp2p0m2m11m2p0k+&
          0.5_wp*d2_rhs_v_dxdx_0_0_nxp2p0m21m2p0k_nxp2p0m2p11m2p0k

d2_rhs_v_dxdx_0_0_nxp2p0m21m2p0k = d2_rhs_v_dxdx_0_0_nxp2p0m21m2p0k*param_float(1)

d1_rhs_v_dx_1_nxp2p0p01m2p0k = q(nx+2+0+0,1-2+0,indvars(1))*q(nx+2+0+0,1-2+0,indvars(2))*q(nx+2+0+0,1-2+0,indvars(3))

d1_rhs_v_dx_1_nxp2p0m11m2p0k = q(nx+2+0-1,1-2+0,indvars(1))*q(nx+2+0-1,1-2+0,indvars(2))*q(nx+2+0-1,1-2+0,indvars(3))

d1_rhs_v_dx_1_nxp2p0m21m2p0k = q(nx+2+0-2,1-2+0,indvars(1))*q(nx+2+0-2,1-2+0,indvars(2))*q(nx+2+0-2,1-2+0,indvars(3))

d1_rhs_v_dx_1_nxp2p01m2p0k = 1.5_wp*d1_rhs_v_dx_1_nxp2p0p01m2p0k-&
          2.0_wp*d1_rhs_v_dx_1_nxp2p0m11m2p0k+&
          0.5_wp*d1_rhs_v_dx_1_nxp2p0m21m2p0k

d1_rhs_v_dx_1_nxp2p01m2p0k = d1_rhs_v_dx_1_nxp2p01m2p0k*param_float(1)

d1_rhs_v_dx_0_nxp2p0p01m2p0k = -((1+&
                    param_float(21 + 5))/(((q(nx+2+0+0,1-2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0+0,1-2+0,indvars(2))*q(nx+2+0+0,1-2+0,indvars(2))+&
                    q(nx+2+0+0,1-2+0,indvars(3))*q(nx+2+0+0,1-2+0,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2+0+0,1-2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0+0,1-2+0,indvars(2))*q(nx+2+0+0,1-2+0,indvars(2))+&
                    q(nx+2+0+0,1-2+0,indvars(3))*q(nx+2+0+0,1-2+0,indvars(3)))))/param_float(4 + 5)**1.5)*(1+&
                    ((q(nx+2+0+0,1-2+0,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+2+0+0,1-2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0+0,1-2+0,indvars(2))*q(nx+2+0+0,1-2+0,indvars(2))+&
                    q(nx+2+0+0,1-2+0,indvars(3))*q(nx+2+0+0,1-2+0,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2+0+0,1-2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0+0,1-2+0,indvars(2))*q(nx+2+0+0,1-2+0,indvars(2))+&
                    q(nx+2+0+0,1-2+0,indvars(3))*q(nx+2+0+0,1-2+0,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+2+0+0,1-2+0,indvars(1)))**3/((q(nx+2+0+0,1-2+0,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+2+0+0,1-2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0+0,1-2+0,indvars(2))*q(nx+2+0+0,1-2+0,indvars(2))+&
                    q(nx+2+0+0,1-2+0,indvars(3))*q(nx+2+0+0,1-2+0,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2+0+0,1-2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0+0,1-2+0,indvars(2))*q(nx+2+0+0,1-2+0,indvars(2))+&
                    q(nx+2+0+0,1-2+0,indvars(3))*q(nx+2+0+0,1-2+0,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+2+0+0,1-2+0,indvars(1)))**3+&
                    param_float(13 + 5)**3))*(q(nx+2+0+0,1-2+0,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+2+0+0,1-2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0+0,1-2+0,indvars(2))*q(nx+2+0+0,1-2+0,indvars(2))+&
                    q(nx+2+0+0,1-2+0,indvars(3))*q(nx+2+0+0,1-2+0,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2+0+0,1-2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0+0,1-2+0,indvars(2))*q(nx+2+0+0,1-2+0,indvars(2))+&
                    q(nx+2+0+0,1-2+0,indvars(3))*q(nx+2+0+0,1-2+0,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+2+0+0,1-2+0,indvars(1))))*param_float(1 + 5)*(d2_rhs_v_dxdx_0_0_nxp2p0p01m2p0k)*qst(nx+2+0+0,1-2+0,indvarsst(10))

d1_rhs_v_dx_0_nxp2p0m11m2p0k = -((1+&
                    param_float(21 + 5))/(((q(nx+2+0-1,1-2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0-1,1-2+0,indvars(2))*q(nx+2+0-1,1-2+0,indvars(2))+&
                    q(nx+2+0-1,1-2+0,indvars(3))*q(nx+2+0-1,1-2+0,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2+0-1,1-2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0-1,1-2+0,indvars(2))*q(nx+2+0-1,1-2+0,indvars(2))+&
                    q(nx+2+0-1,1-2+0,indvars(3))*q(nx+2+0-1,1-2+0,indvars(3)))))/param_float(4 + 5)**1.5)*(1+&
                    ((q(nx+2+0-1,1-2+0,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+2+0-1,1-2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0-1,1-2+0,indvars(2))*q(nx+2+0-1,1-2+0,indvars(2))+&
                    q(nx+2+0-1,1-2+0,indvars(3))*q(nx+2+0-1,1-2+0,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2+0-1,1-2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0-1,1-2+0,indvars(2))*q(nx+2+0-1,1-2+0,indvars(2))+&
                    q(nx+2+0-1,1-2+0,indvars(3))*q(nx+2+0-1,1-2+0,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+2+0-1,1-2+0,indvars(1)))**3/((q(nx+2+0-1,1-2+0,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+2+0-1,1-2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0-1,1-2+0,indvars(2))*q(nx+2+0-1,1-2+0,indvars(2))+&
                    q(nx+2+0-1,1-2+0,indvars(3))*q(nx+2+0-1,1-2+0,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2+0-1,1-2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0-1,1-2+0,indvars(2))*q(nx+2+0-1,1-2+0,indvars(2))+&
                    q(nx+2+0-1,1-2+0,indvars(3))*q(nx+2+0-1,1-2+0,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+2+0-1,1-2+0,indvars(1)))**3+&
                    param_float(13 + 5)**3))*(q(nx+2+0-1,1-2+0,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+2+0-1,1-2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0-1,1-2+0,indvars(2))*q(nx+2+0-1,1-2+0,indvars(2))+&
                    q(nx+2+0-1,1-2+0,indvars(3))*q(nx+2+0-1,1-2+0,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2+0-1,1-2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0-1,1-2+0,indvars(2))*q(nx+2+0-1,1-2+0,indvars(2))+&
                    q(nx+2+0-1,1-2+0,indvars(3))*q(nx+2+0-1,1-2+0,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+2+0-1,1-2+0,indvars(1))))*param_float(1 + 5)*(d2_rhs_v_dxdx_0_0_nxp2p0m11m2p0k)*qst(nx+2+0-1,1-2+0,indvarsst(10))

d1_rhs_v_dx_0_nxp2p0m21m2p0k = -((1+&
                    param_float(21 + 5))/(((q(nx+2+0-2,1-2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0-2,1-2+0,indvars(2))*q(nx+2+0-2,1-2+0,indvars(2))+&
                    q(nx+2+0-2,1-2+0,indvars(3))*q(nx+2+0-2,1-2+0,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2+0-2,1-2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0-2,1-2+0,indvars(2))*q(nx+2+0-2,1-2+0,indvars(2))+&
                    q(nx+2+0-2,1-2+0,indvars(3))*q(nx+2+0-2,1-2+0,indvars(3)))))/param_float(4 + 5)**1.5)*(1+&
                    ((q(nx+2+0-2,1-2+0,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+2+0-2,1-2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0-2,1-2+0,indvars(2))*q(nx+2+0-2,1-2+0,indvars(2))+&
                    q(nx+2+0-2,1-2+0,indvars(3))*q(nx+2+0-2,1-2+0,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2+0-2,1-2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0-2,1-2+0,indvars(2))*q(nx+2+0-2,1-2+0,indvars(2))+&
                    q(nx+2+0-2,1-2+0,indvars(3))*q(nx+2+0-2,1-2+0,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+2+0-2,1-2+0,indvars(1)))**3/((q(nx+2+0-2,1-2+0,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+2+0-2,1-2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0-2,1-2+0,indvars(2))*q(nx+2+0-2,1-2+0,indvars(2))+&
                    q(nx+2+0-2,1-2+0,indvars(3))*q(nx+2+0-2,1-2+0,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2+0-2,1-2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0-2,1-2+0,indvars(2))*q(nx+2+0-2,1-2+0,indvars(2))+&
                    q(nx+2+0-2,1-2+0,indvars(3))*q(nx+2+0-2,1-2+0,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+2+0-2,1-2+0,indvars(1)))**3+&
                    param_float(13 + 5)**3))*(q(nx+2+0-2,1-2+0,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+2+0-2,1-2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0-2,1-2+0,indvars(2))*q(nx+2+0-2,1-2+0,indvars(2))+&
                    q(nx+2+0-2,1-2+0,indvars(3))*q(nx+2+0-2,1-2+0,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2+0-2,1-2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0-2,1-2+0,indvars(2))*q(nx+2+0-2,1-2+0,indvars(2))+&
                    q(nx+2+0-2,1-2+0,indvars(3))*q(nx+2+0-2,1-2+0,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+2+0-2,1-2+0,indvars(1))))*param_float(1 + 5)*(d2_rhs_v_dxdx_0_0_nxp2p0m21m2p0k)*qst(nx+2+0-2,1-2+0,indvarsst(10))

d1_rhs_v_dx_0_nxp2p01m2p0k = 1.5_wp*d1_rhs_v_dx_0_nxp2p0p01m2p0k-&
          2.0_wp*d1_rhs_v_dx_0_nxp2p0m11m2p0k+&
          0.5_wp*d1_rhs_v_dx_0_nxp2p0m21m2p0k

d1_rhs_v_dx_0_nxp2p01m2p0k = d1_rhs_v_dx_0_nxp2p01m2p0k*param_float(1)



!***********************************************************
!                                                           
! Update BC terms for layer 0 0 None d(rho v)/dt ***********
!                                                           
!***********************************************************


rhs(nx+2+0,1-2+0,indvars(3)) =   -  ( (qst(nx+2+0,1-2+0,indvarsst(10))*(d1_rhs_v_dx_0_nxp2p01m2p0k)+&
                    qst(nx+2+0,1-2+0,indvarsst(10))*(d1_rhs_v_dx_1_nxp2p01m2p0k))*qst(nx+2+0,1-2+0,indvarsst(5)) ) 



!***********************************************************
!                                                           
! building source terms in RHS for layer 0 0 None d(rho et)/dt 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (deltaxI*([-4.0_wp/3.0_wp*visc_t*({u}_1x)*deltaxI*u-visc_t*({v}_1x)*deltaxI*v-kappa*({T}_1x)*deltaxI]_1x)+deltaxI*([(rho*et+p)*u]_1x))*symm+(-visc_t*([u]_1y)*([u]_1y)*deltayI**2-4.0_wp/3.0_wp*visc_t*([v]_1y)*([v]_1y)*deltayI**2-kappa*([({T}_1y)*deltayI]_1y)-kappa*([({T}_1x)*deltaxI]_1x)+(rho*et+p)*([v]_1y)*deltayI)*wall
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d2_rhs_et_dxdx_0_2_nxp2p0p01m2p0k_nxp2p0p0p01m2p0k = ((q(nx+2+0+0+0,1-2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0+0+0,1-2+0,indvars(2))*q(nx+2+0+0+0,1-2+0,indvars(2))+&
                    q(nx+2+0+0+0,1-2+0,indvars(3))*q(nx+2+0+0+0,1-2+0,indvars(3)))))/param_float(4 + 5)

d2_rhs_et_dxdx_0_2_nxp2p0p01m2p0k_nxp2p0p0m11m2p0k = ((q(nx+2+0+0-1,1-2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0+0-1,1-2+0,indvars(2))*q(nx+2+0+0-1,1-2+0,indvars(2))+&
                    q(nx+2+0+0-1,1-2+0,indvars(3))*q(nx+2+0+0-1,1-2+0,indvars(3)))))/param_float(4 + 5)

d2_rhs_et_dxdx_0_2_nxp2p0p01m2p0k_nxp2p0p0m21m2p0k = ((q(nx+2+0+0-2,1-2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0+0-2,1-2+0,indvars(2))*q(nx+2+0+0-2,1-2+0,indvars(2))+&
                    q(nx+2+0+0-2,1-2+0,indvars(3))*q(nx+2+0+0-2,1-2+0,indvars(3)))))/param_float(4 + 5)

d2_rhs_et_dxdx_0_2_nxp2p0p01m2p0k = 1.5_wp*d2_rhs_et_dxdx_0_2_nxp2p0p01m2p0k_nxp2p0p0p01m2p0k-&
          2.0_wp*d2_rhs_et_dxdx_0_2_nxp2p0p01m2p0k_nxp2p0p0m11m2p0k+&
          0.5_wp*d2_rhs_et_dxdx_0_2_nxp2p0p01m2p0k_nxp2p0p0m21m2p0k

d2_rhs_et_dxdx_0_2_nxp2p0p01m2p0k = d2_rhs_et_dxdx_0_2_nxp2p0p01m2p0k*param_float(1)

d2_rhs_et_dxdx_0_2_nxp2p0m11m2p0k_nxp2p0m1m11m2p0k = ((q(nx+2+0-1-1,1-2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0-1-1,1-2+0,indvars(2))*q(nx+2+0-1-1,1-2+0,indvars(2))+&
                    q(nx+2+0-1-1,1-2+0,indvars(3))*q(nx+2+0-1-1,1-2+0,indvars(3)))))/param_float(4 + 5)

d2_rhs_et_dxdx_0_2_nxp2p0m11m2p0k_nxp2p0m1p11m2p0k = ((q(nx+2+0-1+1,1-2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0-1+1,1-2+0,indvars(2))*q(nx+2+0-1+1,1-2+0,indvars(2))+&
                    q(nx+2+0-1+1,1-2+0,indvars(3))*q(nx+2+0-1+1,1-2+0,indvars(3)))))/param_float(4 + 5)

d2_rhs_et_dxdx_0_2_nxp2p0m11m2p0k = -&
          0.5_wp*d2_rhs_et_dxdx_0_2_nxp2p0m11m2p0k_nxp2p0m1m11m2p0k+&
          0.5_wp*d2_rhs_et_dxdx_0_2_nxp2p0m11m2p0k_nxp2p0m1p11m2p0k

d2_rhs_et_dxdx_0_2_nxp2p0m11m2p0k = d2_rhs_et_dxdx_0_2_nxp2p0m11m2p0k*param_float(1)

d2_rhs_et_dxdx_0_2_nxp2p0m21m2p0k_nxp2p0m2m11m2p0k = ((q(nx+2+0-2-1,1-2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0-2-1,1-2+0,indvars(2))*q(nx+2+0-2-1,1-2+0,indvars(2))+&
                    q(nx+2+0-2-1,1-2+0,indvars(3))*q(nx+2+0-2-1,1-2+0,indvars(3)))))/param_float(4 + 5)

d2_rhs_et_dxdx_0_2_nxp2p0m21m2p0k_nxp2p0m2p11m2p0k = ((q(nx+2+0-2+1,1-2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0-2+1,1-2+0,indvars(2))*q(nx+2+0-2+1,1-2+0,indvars(2))+&
                    q(nx+2+0-2+1,1-2+0,indvars(3))*q(nx+2+0-2+1,1-2+0,indvars(3)))))/param_float(4 + 5)

d2_rhs_et_dxdx_0_2_nxp2p0m21m2p0k = -&
          0.5_wp*d2_rhs_et_dxdx_0_2_nxp2p0m21m2p0k_nxp2p0m2m11m2p0k+&
          0.5_wp*d2_rhs_et_dxdx_0_2_nxp2p0m21m2p0k_nxp2p0m2p11m2p0k

d2_rhs_et_dxdx_0_2_nxp2p0m21m2p0k = d2_rhs_et_dxdx_0_2_nxp2p0m21m2p0k*param_float(1)

d1_rhs_et_dx_1_nxp2p0p01m2p0k = (q(nx+2+0+0,1-2+0,indvars(1))*q(nx+2+0+0,1-2+0,indvars(4))+&
                    (param_float(3 + 5))*q(nx+2+0+0,1-2+0,indvars(1))*((q(nx+2+0+0,1-2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0+0,1-2+0,indvars(2))*q(nx+2+0+0,1-2+0,indvars(2))+&
                    q(nx+2+0+0,1-2+0,indvars(3))*q(nx+2+0+0,1-2+0,indvars(3))))))*q(nx+2+0+0,1-2+0,indvars(2))

d1_rhs_et_dx_1_nxp2p0m11m2p0k = (q(nx+2+0-1,1-2+0,indvars(1))*q(nx+2+0-1,1-2+0,indvars(4))+&
                    (param_float(3 + 5))*q(nx+2+0-1,1-2+0,indvars(1))*((q(nx+2+0-1,1-2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0-1,1-2+0,indvars(2))*q(nx+2+0-1,1-2+0,indvars(2))+&
                    q(nx+2+0-1,1-2+0,indvars(3))*q(nx+2+0-1,1-2+0,indvars(3))))))*q(nx+2+0-1,1-2+0,indvars(2))

d1_rhs_et_dx_1_nxp2p0m21m2p0k = (q(nx+2+0-2,1-2+0,indvars(1))*q(nx+2+0-2,1-2+0,indvars(4))+&
                    (param_float(3 + 5))*q(nx+2+0-2,1-2+0,indvars(1))*((q(nx+2+0-2,1-2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0-2,1-2+0,indvars(2))*q(nx+2+0-2,1-2+0,indvars(2))+&
                    q(nx+2+0-2,1-2+0,indvars(3))*q(nx+2+0-2,1-2+0,indvars(3))))))*q(nx+2+0-2,1-2+0,indvars(2))

d1_rhs_et_dx_1_nxp2p01m2p0k = 1.5_wp*d1_rhs_et_dx_1_nxp2p0p01m2p0k-&
          2.0_wp*d1_rhs_et_dx_1_nxp2p0m11m2p0k+&
          0.5_wp*d1_rhs_et_dx_1_nxp2p0m21m2p0k

d1_rhs_et_dx_1_nxp2p01m2p0k = d1_rhs_et_dx_1_nxp2p01m2p0k*param_float(1)

d1_rhs_et_dx_0_nxp2p0p01m2p0k = -4.0_wp/3.0_wp*((1+&
                    param_float(21 + 5))/(((q(nx+2+0+0,1-2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0+0,1-2+0,indvars(2))*q(nx+2+0+0,1-2+0,indvars(2))+&
                    q(nx+2+0+0,1-2+0,indvars(3))*q(nx+2+0+0,1-2+0,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2+0+0,1-2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0+0,1-2+0,indvars(2))*q(nx+2+0+0,1-2+0,indvars(2))+&
                    q(nx+2+0+0,1-2+0,indvars(3))*q(nx+2+0+0,1-2+0,indvars(3)))))/param_float(4 + 5)**1.5)*(1+&
                    ((q(nx+2+0+0,1-2+0,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+2+0+0,1-2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0+0,1-2+0,indvars(2))*q(nx+2+0+0,1-2+0,indvars(2))+&
                    q(nx+2+0+0,1-2+0,indvars(3))*q(nx+2+0+0,1-2+0,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2+0+0,1-2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0+0,1-2+0,indvars(2))*q(nx+2+0+0,1-2+0,indvars(2))+&
                    q(nx+2+0+0,1-2+0,indvars(3))*q(nx+2+0+0,1-2+0,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+2+0+0,1-2+0,indvars(1)))**3/((q(nx+2+0+0,1-2+0,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+2+0+0,1-2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0+0,1-2+0,indvars(2))*q(nx+2+0+0,1-2+0,indvars(2))+&
                    q(nx+2+0+0,1-2+0,indvars(3))*q(nx+2+0+0,1-2+0,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2+0+0,1-2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0+0,1-2+0,indvars(2))*q(nx+2+0+0,1-2+0,indvars(2))+&
                    q(nx+2+0+0,1-2+0,indvars(3))*q(nx+2+0+0,1-2+0,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+2+0+0,1-2+0,indvars(1)))**3+&
                    param_float(13 + 5)**3))*(q(nx+2+0+0,1-2+0,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+2+0+0,1-2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0+0,1-2+0,indvars(2))*q(nx+2+0+0,1-2+0,indvars(2))+&
                    q(nx+2+0+0,1-2+0,indvars(3))*q(nx+2+0+0,1-2+0,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2+0+0,1-2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0+0,1-2+0,indvars(2))*q(nx+2+0+0,1-2+0,indvars(2))+&
                    q(nx+2+0+0,1-2+0,indvars(3))*q(nx+2+0+0,1-2+0,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+2+0+0,1-2+0,indvars(1))))*param_float(1 + 5)*(d2_rhs_u_dxdx_0_0_nxp2p0p01m2p0k)*qst(nx+2+0+0,1-2+0,indvarsst(10))*q(nx+2+0+0,1-2+0,indvars(2))-&
                    ((1+&
                    param_float(21 + 5))/(((q(nx+2+0+0,1-2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0+0,1-2+0,indvars(2))*q(nx+2+0+0,1-2+0,indvars(2))+&
                    q(nx+2+0+0,1-2+0,indvars(3))*q(nx+2+0+0,1-2+0,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2+0+0,1-2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0+0,1-2+0,indvars(2))*q(nx+2+0+0,1-2+0,indvars(2))+&
                    q(nx+2+0+0,1-2+0,indvars(3))*q(nx+2+0+0,1-2+0,indvars(3)))))/param_float(4 + 5)**1.5)*(1+&
                    ((q(nx+2+0+0,1-2+0,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+2+0+0,1-2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0+0,1-2+0,indvars(2))*q(nx+2+0+0,1-2+0,indvars(2))+&
                    q(nx+2+0+0,1-2+0,indvars(3))*q(nx+2+0+0,1-2+0,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2+0+0,1-2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0+0,1-2+0,indvars(2))*q(nx+2+0+0,1-2+0,indvars(2))+&
                    q(nx+2+0+0,1-2+0,indvars(3))*q(nx+2+0+0,1-2+0,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+2+0+0,1-2+0,indvars(1)))**3/((q(nx+2+0+0,1-2+0,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+2+0+0,1-2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0+0,1-2+0,indvars(2))*q(nx+2+0+0,1-2+0,indvars(2))+&
                    q(nx+2+0+0,1-2+0,indvars(3))*q(nx+2+0+0,1-2+0,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2+0+0,1-2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0+0,1-2+0,indvars(2))*q(nx+2+0+0,1-2+0,indvars(2))+&
                    q(nx+2+0+0,1-2+0,indvars(3))*q(nx+2+0+0,1-2+0,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+2+0+0,1-2+0,indvars(1)))**3+&
                    param_float(13 + 5)**3))*(q(nx+2+0+0,1-2+0,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+2+0+0,1-2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0+0,1-2+0,indvars(2))*q(nx+2+0+0,1-2+0,indvars(2))+&
                    q(nx+2+0+0,1-2+0,indvars(3))*q(nx+2+0+0,1-2+0,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2+0+0,1-2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0+0,1-2+0,indvars(2))*q(nx+2+0+0,1-2+0,indvars(2))+&
                    q(nx+2+0+0,1-2+0,indvars(3))*q(nx+2+0+0,1-2+0,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+2+0+0,1-2+0,indvars(1))))*param_float(1 + 5)*(d2_rhs_v_dxdx_0_0_nxp2p0p01m2p0k)*qst(nx+2+0+0,1-2+0,indvarsst(10))*q(nx+2+0+0,1-2+0,indvars(3))-&
                    param_float(2 + 5)*(d2_rhs_et_dxdx_0_2_nxp2p0p01m2p0k)*qst(nx+2+0+0,1-2+0,indvarsst(10))

d1_rhs_et_dx_0_nxp2p0m11m2p0k = -4.0_wp/3.0_wp*((1+&
                    param_float(21 + 5))/(((q(nx+2+0-1,1-2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0-1,1-2+0,indvars(2))*q(nx+2+0-1,1-2+0,indvars(2))+&
                    q(nx+2+0-1,1-2+0,indvars(3))*q(nx+2+0-1,1-2+0,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2+0-1,1-2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0-1,1-2+0,indvars(2))*q(nx+2+0-1,1-2+0,indvars(2))+&
                    q(nx+2+0-1,1-2+0,indvars(3))*q(nx+2+0-1,1-2+0,indvars(3)))))/param_float(4 + 5)**1.5)*(1+&
                    ((q(nx+2+0-1,1-2+0,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+2+0-1,1-2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0-1,1-2+0,indvars(2))*q(nx+2+0-1,1-2+0,indvars(2))+&
                    q(nx+2+0-1,1-2+0,indvars(3))*q(nx+2+0-1,1-2+0,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2+0-1,1-2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0-1,1-2+0,indvars(2))*q(nx+2+0-1,1-2+0,indvars(2))+&
                    q(nx+2+0-1,1-2+0,indvars(3))*q(nx+2+0-1,1-2+0,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+2+0-1,1-2+0,indvars(1)))**3/((q(nx+2+0-1,1-2+0,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+2+0-1,1-2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0-1,1-2+0,indvars(2))*q(nx+2+0-1,1-2+0,indvars(2))+&
                    q(nx+2+0-1,1-2+0,indvars(3))*q(nx+2+0-1,1-2+0,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2+0-1,1-2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0-1,1-2+0,indvars(2))*q(nx+2+0-1,1-2+0,indvars(2))+&
                    q(nx+2+0-1,1-2+0,indvars(3))*q(nx+2+0-1,1-2+0,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+2+0-1,1-2+0,indvars(1)))**3+&
                    param_float(13 + 5)**3))*(q(nx+2+0-1,1-2+0,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+2+0-1,1-2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0-1,1-2+0,indvars(2))*q(nx+2+0-1,1-2+0,indvars(2))+&
                    q(nx+2+0-1,1-2+0,indvars(3))*q(nx+2+0-1,1-2+0,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2+0-1,1-2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0-1,1-2+0,indvars(2))*q(nx+2+0-1,1-2+0,indvars(2))+&
                    q(nx+2+0-1,1-2+0,indvars(3))*q(nx+2+0-1,1-2+0,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+2+0-1,1-2+0,indvars(1))))*param_float(1 + 5)*(d2_rhs_u_dxdx_0_0_nxp2p0m11m2p0k)*qst(nx+2+0-1,1-2+0,indvarsst(10))*q(nx+2+0-1,1-2+0,indvars(2))-&
                    ((1+&
                    param_float(21 + 5))/(((q(nx+2+0-1,1-2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0-1,1-2+0,indvars(2))*q(nx+2+0-1,1-2+0,indvars(2))+&
                    q(nx+2+0-1,1-2+0,indvars(3))*q(nx+2+0-1,1-2+0,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2+0-1,1-2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0-1,1-2+0,indvars(2))*q(nx+2+0-1,1-2+0,indvars(2))+&
                    q(nx+2+0-1,1-2+0,indvars(3))*q(nx+2+0-1,1-2+0,indvars(3)))))/param_float(4 + 5)**1.5)*(1+&
                    ((q(nx+2+0-1,1-2+0,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+2+0-1,1-2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0-1,1-2+0,indvars(2))*q(nx+2+0-1,1-2+0,indvars(2))+&
                    q(nx+2+0-1,1-2+0,indvars(3))*q(nx+2+0-1,1-2+0,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2+0-1,1-2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0-1,1-2+0,indvars(2))*q(nx+2+0-1,1-2+0,indvars(2))+&
                    q(nx+2+0-1,1-2+0,indvars(3))*q(nx+2+0-1,1-2+0,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+2+0-1,1-2+0,indvars(1)))**3/((q(nx+2+0-1,1-2+0,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+2+0-1,1-2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0-1,1-2+0,indvars(2))*q(nx+2+0-1,1-2+0,indvars(2))+&
                    q(nx+2+0-1,1-2+0,indvars(3))*q(nx+2+0-1,1-2+0,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2+0-1,1-2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0-1,1-2+0,indvars(2))*q(nx+2+0-1,1-2+0,indvars(2))+&
                    q(nx+2+0-1,1-2+0,indvars(3))*q(nx+2+0-1,1-2+0,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+2+0-1,1-2+0,indvars(1)))**3+&
                    param_float(13 + 5)**3))*(q(nx+2+0-1,1-2+0,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+2+0-1,1-2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0-1,1-2+0,indvars(2))*q(nx+2+0-1,1-2+0,indvars(2))+&
                    q(nx+2+0-1,1-2+0,indvars(3))*q(nx+2+0-1,1-2+0,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2+0-1,1-2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0-1,1-2+0,indvars(2))*q(nx+2+0-1,1-2+0,indvars(2))+&
                    q(nx+2+0-1,1-2+0,indvars(3))*q(nx+2+0-1,1-2+0,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+2+0-1,1-2+0,indvars(1))))*param_float(1 + 5)*(d2_rhs_v_dxdx_0_0_nxp2p0m11m2p0k)*qst(nx+2+0-1,1-2+0,indvarsst(10))*q(nx+2+0-1,1-2+0,indvars(3))-&
                    param_float(2 + 5)*(d2_rhs_et_dxdx_0_2_nxp2p0m11m2p0k)*qst(nx+2+0-1,1-2+0,indvarsst(10))

d1_rhs_et_dx_0_nxp2p0m21m2p0k = -4.0_wp/3.0_wp*((1+&
                    param_float(21 + 5))/(((q(nx+2+0-2,1-2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0-2,1-2+0,indvars(2))*q(nx+2+0-2,1-2+0,indvars(2))+&
                    q(nx+2+0-2,1-2+0,indvars(3))*q(nx+2+0-2,1-2+0,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2+0-2,1-2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0-2,1-2+0,indvars(2))*q(nx+2+0-2,1-2+0,indvars(2))+&
                    q(nx+2+0-2,1-2+0,indvars(3))*q(nx+2+0-2,1-2+0,indvars(3)))))/param_float(4 + 5)**1.5)*(1+&
                    ((q(nx+2+0-2,1-2+0,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+2+0-2,1-2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0-2,1-2+0,indvars(2))*q(nx+2+0-2,1-2+0,indvars(2))+&
                    q(nx+2+0-2,1-2+0,indvars(3))*q(nx+2+0-2,1-2+0,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2+0-2,1-2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0-2,1-2+0,indvars(2))*q(nx+2+0-2,1-2+0,indvars(2))+&
                    q(nx+2+0-2,1-2+0,indvars(3))*q(nx+2+0-2,1-2+0,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+2+0-2,1-2+0,indvars(1)))**3/((q(nx+2+0-2,1-2+0,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+2+0-2,1-2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0-2,1-2+0,indvars(2))*q(nx+2+0-2,1-2+0,indvars(2))+&
                    q(nx+2+0-2,1-2+0,indvars(3))*q(nx+2+0-2,1-2+0,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2+0-2,1-2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0-2,1-2+0,indvars(2))*q(nx+2+0-2,1-2+0,indvars(2))+&
                    q(nx+2+0-2,1-2+0,indvars(3))*q(nx+2+0-2,1-2+0,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+2+0-2,1-2+0,indvars(1)))**3+&
                    param_float(13 + 5)**3))*(q(nx+2+0-2,1-2+0,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+2+0-2,1-2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0-2,1-2+0,indvars(2))*q(nx+2+0-2,1-2+0,indvars(2))+&
                    q(nx+2+0-2,1-2+0,indvars(3))*q(nx+2+0-2,1-2+0,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2+0-2,1-2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0-2,1-2+0,indvars(2))*q(nx+2+0-2,1-2+0,indvars(2))+&
                    q(nx+2+0-2,1-2+0,indvars(3))*q(nx+2+0-2,1-2+0,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+2+0-2,1-2+0,indvars(1))))*param_float(1 + 5)*(d2_rhs_u_dxdx_0_0_nxp2p0m21m2p0k)*qst(nx+2+0-2,1-2+0,indvarsst(10))*q(nx+2+0-2,1-2+0,indvars(2))-&
                    ((1+&
                    param_float(21 + 5))/(((q(nx+2+0-2,1-2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0-2,1-2+0,indvars(2))*q(nx+2+0-2,1-2+0,indvars(2))+&
                    q(nx+2+0-2,1-2+0,indvars(3))*q(nx+2+0-2,1-2+0,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2+0-2,1-2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0-2,1-2+0,indvars(2))*q(nx+2+0-2,1-2+0,indvars(2))+&
                    q(nx+2+0-2,1-2+0,indvars(3))*q(nx+2+0-2,1-2+0,indvars(3)))))/param_float(4 + 5)**1.5)*(1+&
                    ((q(nx+2+0-2,1-2+0,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+2+0-2,1-2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0-2,1-2+0,indvars(2))*q(nx+2+0-2,1-2+0,indvars(2))+&
                    q(nx+2+0-2,1-2+0,indvars(3))*q(nx+2+0-2,1-2+0,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2+0-2,1-2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0-2,1-2+0,indvars(2))*q(nx+2+0-2,1-2+0,indvars(2))+&
                    q(nx+2+0-2,1-2+0,indvars(3))*q(nx+2+0-2,1-2+0,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+2+0-2,1-2+0,indvars(1)))**3/((q(nx+2+0-2,1-2+0,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+2+0-2,1-2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0-2,1-2+0,indvars(2))*q(nx+2+0-2,1-2+0,indvars(2))+&
                    q(nx+2+0-2,1-2+0,indvars(3))*q(nx+2+0-2,1-2+0,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2+0-2,1-2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0-2,1-2+0,indvars(2))*q(nx+2+0-2,1-2+0,indvars(2))+&
                    q(nx+2+0-2,1-2+0,indvars(3))*q(nx+2+0-2,1-2+0,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+2+0-2,1-2+0,indvars(1)))**3+&
                    param_float(13 + 5)**3))*(q(nx+2+0-2,1-2+0,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+2+0-2,1-2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0-2,1-2+0,indvars(2))*q(nx+2+0-2,1-2+0,indvars(2))+&
                    q(nx+2+0-2,1-2+0,indvars(3))*q(nx+2+0-2,1-2+0,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2+0-2,1-2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0-2,1-2+0,indvars(2))*q(nx+2+0-2,1-2+0,indvars(2))+&
                    q(nx+2+0-2,1-2+0,indvars(3))*q(nx+2+0-2,1-2+0,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+2+0-2,1-2+0,indvars(1))))*param_float(1 + 5)*(d2_rhs_v_dxdx_0_0_nxp2p0m21m2p0k)*qst(nx+2+0-2,1-2+0,indvarsst(10))*q(nx+2+0-2,1-2+0,indvars(3))-&
                    param_float(2 + 5)*(d2_rhs_et_dxdx_0_2_nxp2p0m21m2p0k)*qst(nx+2+0-2,1-2+0,indvarsst(10))

d1_rhs_et_dx_0_nxp2p01m2p0k = 1.5_wp*d1_rhs_et_dx_0_nxp2p0p01m2p0k-&
          2.0_wp*d1_rhs_et_dx_0_nxp2p0m11m2p0k+&
          0.5_wp*d1_rhs_et_dx_0_nxp2p0m21m2p0k

d1_rhs_et_dx_0_nxp2p01m2p0k = d1_rhs_et_dx_0_nxp2p01m2p0k*param_float(1)

d1_rhs_et_dx_2_nxp2p0p01m2p0k = (d2_rhs_et_dxdx_0_2_nxp2p0p01m2p0k)*qst(nx+2+0+0,1-2+0,indvarsst(10))

d1_rhs_et_dx_2_nxp2p0m11m2p0k = (d2_rhs_et_dxdx_0_2_nxp2p0m11m2p0k)*qst(nx+2+0-1,1-2+0,indvarsst(10))

d1_rhs_et_dx_2_nxp2p0m21m2p0k = (d2_rhs_et_dxdx_0_2_nxp2p0m21m2p0k)*qst(nx+2+0-2,1-2+0,indvarsst(10))

d1_rhs_et_dx_2_nxp2p01m2p0k = 1.5_wp*d1_rhs_et_dx_2_nxp2p0p01m2p0k-&
          2.0_wp*d1_rhs_et_dx_2_nxp2p0m11m2p0k+&
          0.5_wp*d1_rhs_et_dx_2_nxp2p0m21m2p0k

d1_rhs_et_dx_2_nxp2p01m2p0k = d1_rhs_et_dx_2_nxp2p01m2p0k*param_float(1)

d1_rhs_et_dy_0_nxp2p01m2p0p0k = q(nx+2+0,1-2+0+0,indvars(2))

d1_rhs_et_dy_0_nxp2p01m2p0p1k = q(nx+2+0,1-2+0+1,indvars(2))

d1_rhs_et_dy_0_nxp2p01m2p0p2k = q(nx+2+0,1-2+0+2,indvars(2))

d1_rhs_et_dy_0_nxp2p01m2p0k = -&
          1.5_wp*d1_rhs_et_dy_0_nxp2p01m2p0p0k+&
          2.0_wp*d1_rhs_et_dy_0_nxp2p01m2p0p1k-&
          0.5_wp*d1_rhs_et_dy_0_nxp2p01m2p0p2k

d1_rhs_et_dy_0_nxp2p01m2p0k = d1_rhs_et_dy_0_nxp2p01m2p0k*param_float(2)

d2_rhs_et_dydy_4_0_nxp2p01m2p0p0k_nxp2p01m2p0p0p0k = ((q(nx+2+0,1-2+0+0+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0,1-2+0+0+0,indvars(2))*q(nx+2+0,1-2+0+0+0,indvars(2))+&
                    q(nx+2+0,1-2+0+0+0,indvars(3))*q(nx+2+0,1-2+0+0+0,indvars(3)))))/param_float(4 + 5)

d2_rhs_et_dydy_4_0_nxp2p01m2p0p0k_nxp2p01m2p0p0p1k = ((q(nx+2+0,1-2+0+0+1,indvars(4))-&
                    0.5_wp*(q(nx+2+0,1-2+0+0+1,indvars(2))*q(nx+2+0,1-2+0+0+1,indvars(2))+&
                    q(nx+2+0,1-2+0+0+1,indvars(3))*q(nx+2+0,1-2+0+0+1,indvars(3)))))/param_float(4 + 5)

d2_rhs_et_dydy_4_0_nxp2p01m2p0p0k_nxp2p01m2p0p0p2k = ((q(nx+2+0,1-2+0+0+2,indvars(4))-&
                    0.5_wp*(q(nx+2+0,1-2+0+0+2,indvars(2))*q(nx+2+0,1-2+0+0+2,indvars(2))+&
                    q(nx+2+0,1-2+0+0+2,indvars(3))*q(nx+2+0,1-2+0+0+2,indvars(3)))))/param_float(4 + 5)

d2_rhs_et_dydy_4_0_nxp2p01m2p0p0k = -&
          1.5_wp*d2_rhs_et_dydy_4_0_nxp2p01m2p0p0k_nxp2p01m2p0p0p0k+&
          2.0_wp*d2_rhs_et_dydy_4_0_nxp2p01m2p0p0k_nxp2p01m2p0p0p1k-&
          0.5_wp*d2_rhs_et_dydy_4_0_nxp2p01m2p0p0k_nxp2p01m2p0p0p2k

d2_rhs_et_dydy_4_0_nxp2p01m2p0p0k = d2_rhs_et_dydy_4_0_nxp2p01m2p0p0k*param_float(2)

d2_rhs_et_dydy_4_0_nxp2p01m2p0p1k_nxp2p01m2p0p1m1k = ((q(nx+2+0,1-2+0+1-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0,1-2+0+1-1,indvars(2))*q(nx+2+0,1-2+0+1-1,indvars(2))+&
                    q(nx+2+0,1-2+0+1-1,indvars(3))*q(nx+2+0,1-2+0+1-1,indvars(3)))))/param_float(4 + 5)

d2_rhs_et_dydy_4_0_nxp2p01m2p0p1k_nxp2p01m2p0p1p1k = ((q(nx+2+0,1-2+0+1+1,indvars(4))-&
                    0.5_wp*(q(nx+2+0,1-2+0+1+1,indvars(2))*q(nx+2+0,1-2+0+1+1,indvars(2))+&
                    q(nx+2+0,1-2+0+1+1,indvars(3))*q(nx+2+0,1-2+0+1+1,indvars(3)))))/param_float(4 + 5)

d2_rhs_et_dydy_4_0_nxp2p01m2p0p1k = -&
          0.5_wp*d2_rhs_et_dydy_4_0_nxp2p01m2p0p1k_nxp2p01m2p0p1m1k+&
          0.5_wp*d2_rhs_et_dydy_4_0_nxp2p01m2p0p1k_nxp2p01m2p0p1p1k

d2_rhs_et_dydy_4_0_nxp2p01m2p0p1k = d2_rhs_et_dydy_4_0_nxp2p01m2p0p1k*param_float(2)

d2_rhs_et_dydy_4_0_nxp2p01m2p0p2k_nxp2p01m2p0p2m1k = ((q(nx+2+0,1-2+0+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0,1-2+0+2-1,indvars(2))*q(nx+2+0,1-2+0+2-1,indvars(2))+&
                    q(nx+2+0,1-2+0+2-1,indvars(3))*q(nx+2+0,1-2+0+2-1,indvars(3)))))/param_float(4 + 5)

d2_rhs_et_dydy_4_0_nxp2p01m2p0p2k_nxp2p01m2p0p2p1k = ((q(nx+2+0,1-2+0+2+1,indvars(4))-&
                    0.5_wp*(q(nx+2+0,1-2+0+2+1,indvars(2))*q(nx+2+0,1-2+0+2+1,indvars(2))+&
                    q(nx+2+0,1-2+0+2+1,indvars(3))*q(nx+2+0,1-2+0+2+1,indvars(3)))))/param_float(4 + 5)

d2_rhs_et_dydy_4_0_nxp2p01m2p0p2k = -&
          0.5_wp*d2_rhs_et_dydy_4_0_nxp2p01m2p0p2k_nxp2p01m2p0p2m1k+&
          0.5_wp*d2_rhs_et_dydy_4_0_nxp2p01m2p0p2k_nxp2p01m2p0p2p1k

d2_rhs_et_dydy_4_0_nxp2p01m2p0p2k = d2_rhs_et_dydy_4_0_nxp2p01m2p0p2k*param_float(2)

d1_rhs_et_dy_4_nxp2p01m2p0p0k = (d2_rhs_et_dydy_4_0_nxp2p01m2p0p0k)*qst(nx+2+0,1-2+0+0,indvarsst(11))

d1_rhs_et_dy_4_nxp2p01m2p0p1k = (d2_rhs_et_dydy_4_0_nxp2p01m2p0p1k)*qst(nx+2+0,1-2+0+1,indvarsst(11))

d1_rhs_et_dy_4_nxp2p01m2p0p2k = (d2_rhs_et_dydy_4_0_nxp2p01m2p0p2k)*qst(nx+2+0,1-2+0+2,indvarsst(11))

d1_rhs_et_dy_4_nxp2p01m2p0k = -&
          1.5_wp*d1_rhs_et_dy_4_nxp2p01m2p0p0k+&
          2.0_wp*d1_rhs_et_dy_4_nxp2p01m2p0p1k-&
          0.5_wp*d1_rhs_et_dy_4_nxp2p01m2p0p2k

d1_rhs_et_dy_4_nxp2p01m2p0k = d1_rhs_et_dy_4_nxp2p01m2p0k*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer 0 0 None d(rho et)/dt **********
!                                                           
!***********************************************************


rhs(nx+2+0,1-2+0,indvars(4)) =   -  ( (qst(nx+2+0,1-2+0,indvarsst(10))*(d1_rhs_et_dx_0_nxp2p01m2p0k)+&
                    qst(nx+2+0,1-2+0,indvarsst(10))*(d1_rhs_et_dx_1_nxp2p01m2p0k))*qst(nx+2+0,1-2+0,indvarsst(5))+&
                    (-&
                    ((1+&
                    param_float(21 + 5))/(((q(nx+2+0,1-2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0,1-2+0,indvars(2))*q(nx+2+0,1-2+0,indvars(2))+&
                    q(nx+2+0,1-2+0,indvars(3))*q(nx+2+0,1-2+0,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2+0,1-2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0,1-2+0,indvars(2))*q(nx+2+0,1-2+0,indvars(2))+&
                    q(nx+2+0,1-2+0,indvars(3))*q(nx+2+0,1-2+0,indvars(3)))))/param_float(4 + 5)**1.5)*(1+&
                    ((q(nx+2+0,1-2+0,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+2+0,1-2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0,1-2+0,indvars(2))*q(nx+2+0,1-2+0,indvars(2))+&
                    q(nx+2+0,1-2+0,indvars(3))*q(nx+2+0,1-2+0,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2+0,1-2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0,1-2+0,indvars(2))*q(nx+2+0,1-2+0,indvars(2))+&
                    q(nx+2+0,1-2+0,indvars(3))*q(nx+2+0,1-2+0,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+2+0,1-2+0,indvars(1)))**3/((q(nx+2+0,1-2+0,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+2+0,1-2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0,1-2+0,indvars(2))*q(nx+2+0,1-2+0,indvars(2))+&
                    q(nx+2+0,1-2+0,indvars(3))*q(nx+2+0,1-2+0,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2+0,1-2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0,1-2+0,indvars(2))*q(nx+2+0,1-2+0,indvars(2))+&
                    q(nx+2+0,1-2+0,indvars(3))*q(nx+2+0,1-2+0,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+2+0,1-2+0,indvars(1)))**3+&
                    param_float(13 + 5)**3))*(q(nx+2+0,1-2+0,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+2+0,1-2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0,1-2+0,indvars(2))*q(nx+2+0,1-2+0,indvars(2))+&
                    q(nx+2+0,1-2+0,indvars(3))*q(nx+2+0,1-2+0,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2+0,1-2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0,1-2+0,indvars(2))*q(nx+2+0,1-2+0,indvars(2))+&
                    q(nx+2+0,1-2+0,indvars(3))*q(nx+2+0,1-2+0,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+2+0,1-2+0,indvars(1))))*param_float(1 + 5)*(d1_rhs_et_dy_0_nxp2p01m2p0k)*(d1_rhs_et_dy_0_nxp2p01m2p0k)*qst(nx+2+0,1-2+0,indvarsst(11))**2-&
                    4.0_wp/3.0_wp*((1+&
                    param_float(21 + 5))/(((q(nx+2+0,1-2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0,1-2+0,indvars(2))*q(nx+2+0,1-2+0,indvars(2))+&
                    q(nx+2+0,1-2+0,indvars(3))*q(nx+2+0,1-2+0,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2+0,1-2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0,1-2+0,indvars(2))*q(nx+2+0,1-2+0,indvars(2))+&
                    q(nx+2+0,1-2+0,indvars(3))*q(nx+2+0,1-2+0,indvars(3)))))/param_float(4 + 5)**1.5)*(1+&
                    ((q(nx+2+0,1-2+0,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+2+0,1-2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0,1-2+0,indvars(2))*q(nx+2+0,1-2+0,indvars(2))+&
                    q(nx+2+0,1-2+0,indvars(3))*q(nx+2+0,1-2+0,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2+0,1-2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0,1-2+0,indvars(2))*q(nx+2+0,1-2+0,indvars(2))+&
                    q(nx+2+0,1-2+0,indvars(3))*q(nx+2+0,1-2+0,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+2+0,1-2+0,indvars(1)))**3/((q(nx+2+0,1-2+0,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+2+0,1-2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0,1-2+0,indvars(2))*q(nx+2+0,1-2+0,indvars(2))+&
                    q(nx+2+0,1-2+0,indvars(3))*q(nx+2+0,1-2+0,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2+0,1-2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0,1-2+0,indvars(2))*q(nx+2+0,1-2+0,indvars(2))+&
                    q(nx+2+0,1-2+0,indvars(3))*q(nx+2+0,1-2+0,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+2+0,1-2+0,indvars(1)))**3+&
                    param_float(13 + 5)**3))*(q(nx+2+0,1-2+0,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+2+0,1-2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0,1-2+0,indvars(2))*q(nx+2+0,1-2+0,indvars(2))+&
                    q(nx+2+0,1-2+0,indvars(3))*q(nx+2+0,1-2+0,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2+0,1-2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0,1-2+0,indvars(2))*q(nx+2+0,1-2+0,indvars(2))+&
                    q(nx+2+0,1-2+0,indvars(3))*q(nx+2+0,1-2+0,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+2+0,1-2+0,indvars(1))))*param_float(1 + 5)*(d1_rhs_rho_dy_0_nxp2p01m2p0k)*(d1_rhs_rho_dy_0_nxp2p01m2p0k)*qst(nx+2+0,1-2+0,indvarsst(11))**2-&
                    param_float(2 + 5)*(d1_rhs_et_dy_4_nxp2p01m2p0k)-&
                    param_float(2 + 5)*(d1_rhs_et_dx_2_nxp2p01m2p0k)+&
                    (q(nx+2+0,1-2+0,indvars(1))*q(nx+2+0,1-2+0,indvars(4))+&
                    (param_float(3 + 5))*q(nx+2+0,1-2+0,indvars(1))*((q(nx+2+0,1-2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0,1-2+0,indvars(2))*q(nx+2+0,1-2+0,indvars(2))+&
                    q(nx+2+0,1-2+0,indvars(3))*q(nx+2+0,1-2+0,indvars(3))))))*(d1_rhs_rho_dy_0_nxp2p01m2p0k)*qst(nx+2+0,1-2+0,indvarsst(11)))*qst(nx+2+0,1-2+0,indvarsst(19)) ) 



!***********************************************************
!                                                           
! building source terms in RHS for layer 0 0 None d(rho nut)/dt 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (-ReI*Cb2*sigmaI*((deltaxI)**2*([rho*nut]_1x)*([nut]_1x))-Cb1*(1-ft2)*SS*rho*nut+ReI*(Cw1*fw-Cb1/k**2*ft2)*rho*nut**2/eta**2+deltaxI*([rho*u*nut-sigmaI*(visc+rho*nut)*({nut}_1x)*deltaxI]_1x))*symm
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_rhs_nut_dx_0_nxp2p0p01m2p0k = q(nx+2+0+0,1-2+0,indvars(1))*q(nx+2+0+0,1-2+0,indvars(5))

d1_rhs_nut_dx_0_nxp2p0m11m2p0k = q(nx+2+0-1,1-2+0,indvars(1))*q(nx+2+0-1,1-2+0,indvars(5))

d1_rhs_nut_dx_0_nxp2p0m21m2p0k = q(nx+2+0-2,1-2+0,indvars(1))*q(nx+2+0-2,1-2+0,indvars(5))

d1_rhs_nut_dx_0_nxp2p01m2p0k = 1.5_wp*d1_rhs_nut_dx_0_nxp2p0p01m2p0k-&
          2.0_wp*d1_rhs_nut_dx_0_nxp2p0m11m2p0k+&
          0.5_wp*d1_rhs_nut_dx_0_nxp2p0m21m2p0k

d1_rhs_nut_dx_0_nxp2p01m2p0k = d1_rhs_nut_dx_0_nxp2p01m2p0k*param_float(1)

d1_rhs_nut_dx_1_nxp2p0p01m2p0k = q(nx+2+0+0,1-2+0,indvars(5))

d1_rhs_nut_dx_1_nxp2p0m11m2p0k = q(nx+2+0-1,1-2+0,indvars(5))

d1_rhs_nut_dx_1_nxp2p0m21m2p0k = q(nx+2+0-2,1-2+0,indvars(5))

d1_rhs_nut_dx_1_nxp2p01m2p0k = 1.5_wp*d1_rhs_nut_dx_1_nxp2p0p01m2p0k-&
          2.0_wp*d1_rhs_nut_dx_1_nxp2p0m11m2p0k+&
          0.5_wp*d1_rhs_nut_dx_1_nxp2p0m21m2p0k

d1_rhs_nut_dx_1_nxp2p01m2p0k = d1_rhs_nut_dx_1_nxp2p01m2p0k*param_float(1)

d2_rhs_nut_dxdx_2_0_nxp2p0p01m2p0k_nxp2p0p0p01m2p0k = q(nx+2+0+0+0,1-2+0,indvars(5))

d2_rhs_nut_dxdx_2_0_nxp2p0p01m2p0k_nxp2p0p0m11m2p0k = q(nx+2+0+0-1,1-2+0,indvars(5))

d2_rhs_nut_dxdx_2_0_nxp2p0p01m2p0k_nxp2p0p0m21m2p0k = q(nx+2+0+0-2,1-2+0,indvars(5))

d2_rhs_nut_dxdx_2_0_nxp2p0p01m2p0k = 1.5_wp*d2_rhs_nut_dxdx_2_0_nxp2p0p01m2p0k_nxp2p0p0p01m2p0k-&
          2.0_wp*d2_rhs_nut_dxdx_2_0_nxp2p0p01m2p0k_nxp2p0p0m11m2p0k+&
          0.5_wp*d2_rhs_nut_dxdx_2_0_nxp2p0p01m2p0k_nxp2p0p0m21m2p0k

d2_rhs_nut_dxdx_2_0_nxp2p0p01m2p0k = d2_rhs_nut_dxdx_2_0_nxp2p0p01m2p0k*param_float(1)

d2_rhs_nut_dxdx_2_0_nxp2p0m11m2p0k_nxp2p0m1m11m2p0k = q(nx+2+0-1-1,1-2+0,indvars(5))

d2_rhs_nut_dxdx_2_0_nxp2p0m11m2p0k_nxp2p0m1p11m2p0k = q(nx+2+0-1+1,1-2+0,indvars(5))

d2_rhs_nut_dxdx_2_0_nxp2p0m11m2p0k = -&
          0.5_wp*d2_rhs_nut_dxdx_2_0_nxp2p0m11m2p0k_nxp2p0m1m11m2p0k+&
          0.5_wp*d2_rhs_nut_dxdx_2_0_nxp2p0m11m2p0k_nxp2p0m1p11m2p0k

d2_rhs_nut_dxdx_2_0_nxp2p0m11m2p0k = d2_rhs_nut_dxdx_2_0_nxp2p0m11m2p0k*param_float(1)

d2_rhs_nut_dxdx_2_0_nxp2p0m21m2p0k_nxp2p0m2m11m2p0k = q(nx+2+0-2-1,1-2+0,indvars(5))

d2_rhs_nut_dxdx_2_0_nxp2p0m21m2p0k_nxp2p0m2p11m2p0k = q(nx+2+0-2+1,1-2+0,indvars(5))

d2_rhs_nut_dxdx_2_0_nxp2p0m21m2p0k = -&
          0.5_wp*d2_rhs_nut_dxdx_2_0_nxp2p0m21m2p0k_nxp2p0m2m11m2p0k+&
          0.5_wp*d2_rhs_nut_dxdx_2_0_nxp2p0m21m2p0k_nxp2p0m2p11m2p0k

d2_rhs_nut_dxdx_2_0_nxp2p0m21m2p0k = d2_rhs_nut_dxdx_2_0_nxp2p0m21m2p0k*param_float(1)

d1_rhs_nut_dx_2_nxp2p0p01m2p0k = q(nx+2+0+0,1-2+0,indvars(1))*q(nx+2+0+0,1-2+0,indvars(2))*q(nx+2+0+0,1-2+0,indvars(5))-&
                    param_float(18 + 5)*((1+&
                    param_float(21 + 5))/(((q(nx+2+0+0,1-2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0+0,1-2+0,indvars(2))*q(nx+2+0+0,1-2+0,indvars(2))+&
                    q(nx+2+0+0,1-2+0,indvars(3))*q(nx+2+0+0,1-2+0,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2+0+0,1-2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0+0,1-2+0,indvars(2))*q(nx+2+0+0,1-2+0,indvars(2))+&
                    q(nx+2+0+0,1-2+0,indvars(3))*q(nx+2+0+0,1-2+0,indvars(3)))))/param_float(4 + 5)**1.5+&
                    q(nx+2+0+0,1-2+0,indvars(1))*q(nx+2+0+0,1-2+0,indvars(5)))*(d2_rhs_nut_dxdx_2_0_nxp2p0p01m2p0k)*qst(nx+2+0+0,1-2+0,indvarsst(10))

d1_rhs_nut_dx_2_nxp2p0m11m2p0k = q(nx+2+0-1,1-2+0,indvars(1))*q(nx+2+0-1,1-2+0,indvars(2))*q(nx+2+0-1,1-2+0,indvars(5))-&
                    param_float(18 + 5)*((1+&
                    param_float(21 + 5))/(((q(nx+2+0-1,1-2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0-1,1-2+0,indvars(2))*q(nx+2+0-1,1-2+0,indvars(2))+&
                    q(nx+2+0-1,1-2+0,indvars(3))*q(nx+2+0-1,1-2+0,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2+0-1,1-2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0-1,1-2+0,indvars(2))*q(nx+2+0-1,1-2+0,indvars(2))+&
                    q(nx+2+0-1,1-2+0,indvars(3))*q(nx+2+0-1,1-2+0,indvars(3)))))/param_float(4 + 5)**1.5+&
                    q(nx+2+0-1,1-2+0,indvars(1))*q(nx+2+0-1,1-2+0,indvars(5)))*(d2_rhs_nut_dxdx_2_0_nxp2p0m11m2p0k)*qst(nx+2+0-1,1-2+0,indvarsst(10))

d1_rhs_nut_dx_2_nxp2p0m21m2p0k = q(nx+2+0-2,1-2+0,indvars(1))*q(nx+2+0-2,1-2+0,indvars(2))*q(nx+2+0-2,1-2+0,indvars(5))-&
                    param_float(18 + 5)*((1+&
                    param_float(21 + 5))/(((q(nx+2+0-2,1-2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0-2,1-2+0,indvars(2))*q(nx+2+0-2,1-2+0,indvars(2))+&
                    q(nx+2+0-2,1-2+0,indvars(3))*q(nx+2+0-2,1-2+0,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2+0-2,1-2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0-2,1-2+0,indvars(2))*q(nx+2+0-2,1-2+0,indvars(2))+&
                    q(nx+2+0-2,1-2+0,indvars(3))*q(nx+2+0-2,1-2+0,indvars(3)))))/param_float(4 + 5)**1.5+&
                    q(nx+2+0-2,1-2+0,indvars(1))*q(nx+2+0-2,1-2+0,indvars(5)))*(d2_rhs_nut_dxdx_2_0_nxp2p0m21m2p0k)*qst(nx+2+0-2,1-2+0,indvarsst(10))

d1_rhs_nut_dx_2_nxp2p01m2p0k = 1.5_wp*d1_rhs_nut_dx_2_nxp2p0p01m2p0k-&
          2.0_wp*d1_rhs_nut_dx_2_nxp2p0m11m2p0k+&
          0.5_wp*d1_rhs_nut_dx_2_nxp2p0m21m2p0k

d1_rhs_nut_dx_2_nxp2p01m2p0k = d1_rhs_nut_dx_2_nxp2p01m2p0k*param_float(1)



!***********************************************************
!                                                           
! Update BC terms for layer 0 0 None d(rho nut)/dt *********
!                                                           
!***********************************************************


rhs(nx+2+0,1-2+0,indvars(5)) =   -  ( (-&
                    param_float(1 + 5)*param_float(7 + 5)*param_float(18 + 5)*((qst(nx+2+0,1-2+0,indvarsst(10)))**2*(d1_rhs_nut_dx_0_nxp2p01m2p0k)*(d1_rhs_nut_dx_1_nxp2p01m2p0k))-&
                    param_float(6 + 5)*(1-&
                    qst(nx+2+0,1-2+0,indvarsst(16)))*qst(nx+2+0,1-2+0,indvarsst(20))*q(nx+2+0,1-2+0,indvars(1))*q(nx+2+0,1-2+0,indvars(5))+&
                    param_float(1 + 5)*(param_float(10 + 5)*qst(nx+2+0,1-2+0,indvarsst(13))-&
                    param_float(6 + 5)/param_float(9 + 5)**2*qst(nx+2+0,1-2+0,indvarsst(16)))*q(nx+2+0,1-2+0,indvars(1))*q(nx+2+0,1-2+0,indvars(5))**2/qst(nx+2+0,1-2+0,indvarsst(2))**2+&
                    qst(nx+2+0,1-2+0,indvarsst(10))*(d1_rhs_nut_dx_2_nxp2p01m2p0k))*qst(nx+2+0,1-2+0,indvarsst(5)) ) 

