

!***********************************************************
!                                                           
! Start building layers for BC : None jmax None ************
!                                                           
!***********************************************************




!***********************************************************
!                                                           
! BC layer: None 1 None ************************************
!                                                           
!***********************************************************


 
      do i=idloop(1),idloop(2) 


!***********************************************************
!                                                           
! building source terms in RHS for layer None 1 None d(rho)/dt 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! deltaxI*([rho*u]_1x)+deltayI*([rho*v]_1y)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_conv_rho_dx_0_im1nyp2m1k = q(i-1,ny+2-1,indvars(1))*q(i-1,ny+2-1,indvars(2))

d1_conv_rho_dx_0_ip1nyp2m1k = q(i+1,ny+2-1,indvars(1))*q(i+1,ny+2-1,indvars(2))

d1_conv_rho_dx_0_inyp2m1k = -&
          0.5_wp*d1_conv_rho_dx_0_im1nyp2m1k+&
          0.5_wp*d1_conv_rho_dx_0_ip1nyp2m1k

d1_conv_rho_dx_0_inyp2m1k = d1_conv_rho_dx_0_inyp2m1k*param_float(1)

d1_conv_rho_dy_0_inyp2m1m1k = q(i,ny+2-1-1,indvars(1))*q(i,ny+2-1-1,indvars(3))

d1_conv_rho_dy_0_inyp2m1p1k = q(i,ny+2-1+1,indvars(1))*q(i,ny+2-1+1,indvars(3))

d1_conv_rho_dy_0_inyp2m1k = -&
          0.5_wp*d1_conv_rho_dy_0_inyp2m1m1k+&
          0.5_wp*d1_conv_rho_dy_0_inyp2m1p1k

d1_conv_rho_dy_0_inyp2m1k = d1_conv_rho_dy_0_inyp2m1k*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer None 1 None d(rho)/dt **********
!                                                           
!***********************************************************


rhs(i,ny+2-1,indvars(1)) =   -  ( qst(i,ny+2-1,indvarsst(10))*(d1_conv_rho_dx_0_inyp2m1k)+&
                    qst(i,ny+2-1,indvarsst(11))*(d1_conv_rho_dy_0_inyp2m1k) ) 



!***********************************************************
!                                                           
! building source terms in RHS for layer None 1 None d(rho u)/dt 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! deltaxI*([rho*u*u+p]_1x)+deltayI*([rho*v*u]_1y)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_conv_rhou_dx_0_im1nyp2m1k = q(i-1,ny+2-1,indvars(1))*q(i-1,ny+2-1,indvars(2))*q(i-1,ny+2-1,indvars(2))+(param_float(3 + 5))*q(i-1,ny+2-1,indvars(1))*((q(i-1,ny+2-1,indvars(4))-&
                    0.5_wp*(q(i-1,ny+2-1,indvars(2))*q(i-1,ny+2-1,indvars(2))+&
                    q(i-1,ny+2-1,indvars(3))*q(i-1,ny+2-1,indvars(3)))))

d1_conv_rhou_dx_0_ip1nyp2m1k = q(i+1,ny+2-1,indvars(1))*q(i+1,ny+2-1,indvars(2))*q(i+1,ny+2-1,indvars(2))+(param_float(3 + 5))*q(i+1,ny+2-1,indvars(1))*((q(i+1,ny+2-1,indvars(4))-&
                    0.5_wp*(q(i+1,ny+2-1,indvars(2))*q(i+1,ny+2-1,indvars(2))+&
                    q(i+1,ny+2-1,indvars(3))*q(i+1,ny+2-1,indvars(3)))))

d1_conv_rhou_dx_0_inyp2m1k = -&
          0.5_wp*d1_conv_rhou_dx_0_im1nyp2m1k+&
          0.5_wp*d1_conv_rhou_dx_0_ip1nyp2m1k

d1_conv_rhou_dx_0_inyp2m1k = d1_conv_rhou_dx_0_inyp2m1k*param_float(1)

d1_conv_rhou_dy_0_inyp2m1m1k = q(i,ny+2-1-1,indvars(1))*q(i,ny+2-1-1,indvars(3))*q(i,ny+2-1-1,indvars(2))

d1_conv_rhou_dy_0_inyp2m1p1k = q(i,ny+2-1+1,indvars(1))*q(i,ny+2-1+1,indvars(3))*q(i,ny+2-1+1,indvars(2))

d1_conv_rhou_dy_0_inyp2m1k = -&
          0.5_wp*d1_conv_rhou_dy_0_inyp2m1m1k+&
          0.5_wp*d1_conv_rhou_dy_0_inyp2m1p1k

d1_conv_rhou_dy_0_inyp2m1k = d1_conv_rhou_dy_0_inyp2m1k*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer None 1 None d(rho u)/dt ********
!                                                           
!***********************************************************


rhs(i,ny+2-1,indvars(2)) =   -  ( qst(i,ny+2-1,indvarsst(10))*(d1_conv_rhou_dx_0_inyp2m1k)+&
                    qst(i,ny+2-1,indvarsst(11))*(d1_conv_rhou_dy_0_inyp2m1k) ) 



!***********************************************************
!                                                           
! building source terms in RHS for layer None 1 None d(rho v)/dt 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! deltaxI*([rho*u*v]_1x)+deltayI*([rho*v*v+p]_1y)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_conv_rhov_dx_0_im1nyp2m1k = q(i-1,ny+2-1,indvars(1))*q(i-1,ny+2-1,indvars(2))*q(i-1,ny+2-1,indvars(3))

d1_conv_rhov_dx_0_ip1nyp2m1k = q(i+1,ny+2-1,indvars(1))*q(i+1,ny+2-1,indvars(2))*q(i+1,ny+2-1,indvars(3))

d1_conv_rhov_dx_0_inyp2m1k = -&
          0.5_wp*d1_conv_rhov_dx_0_im1nyp2m1k+&
          0.5_wp*d1_conv_rhov_dx_0_ip1nyp2m1k

d1_conv_rhov_dx_0_inyp2m1k = d1_conv_rhov_dx_0_inyp2m1k*param_float(1)

d1_conv_rhov_dy_0_inyp2m1m1k = q(i,ny+2-1-1,indvars(1))*q(i,ny+2-1-1,indvars(3))*q(i,ny+2-1-1,indvars(3))+(param_float(3 + 5))*q(i,ny+2-1-1,indvars(1))*((q(i,ny+2-1-1,indvars(4))-&
                    0.5_wp*(q(i,ny+2-1-1,indvars(2))*q(i,ny+2-1-1,indvars(2))+&
                    q(i,ny+2-1-1,indvars(3))*q(i,ny+2-1-1,indvars(3)))))

d1_conv_rhov_dy_0_inyp2m1p1k = q(i,ny+2-1+1,indvars(1))*q(i,ny+2-1+1,indvars(3))*q(i,ny+2-1+1,indvars(3))+(param_float(3 + 5))*q(i,ny+2-1+1,indvars(1))*((q(i,ny+2-1+1,indvars(4))-&
                    0.5_wp*(q(i,ny+2-1+1,indvars(2))*q(i,ny+2-1+1,indvars(2))+&
                    q(i,ny+2-1+1,indvars(3))*q(i,ny+2-1+1,indvars(3)))))

d1_conv_rhov_dy_0_inyp2m1k = -&
          0.5_wp*d1_conv_rhov_dy_0_inyp2m1m1k+&
          0.5_wp*d1_conv_rhov_dy_0_inyp2m1p1k

d1_conv_rhov_dy_0_inyp2m1k = d1_conv_rhov_dy_0_inyp2m1k*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer None 1 None d(rho v)/dt ********
!                                                           
!***********************************************************


rhs(i,ny+2-1,indvars(3)) =   -  ( qst(i,ny+2-1,indvarsst(10))*(d1_conv_rhov_dx_0_inyp2m1k)+&
                    qst(i,ny+2-1,indvarsst(11))*(d1_conv_rhov_dy_0_inyp2m1k) ) 



!***********************************************************
!                                                           
! building source terms in RHS for layer None 1 None d(rho et)/dt 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! deltaxI*([(rho*et+p)*u]_1x)+deltayI*([(rho*et+p)*v]_1y)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_conv_et_dx_0_im1nyp2m1k = (q(i-1,ny+2-1,indvars(1))*q(i-1,ny+2-1,indvars(4))+&
                    (param_float(3 + 5))*q(i-1,ny+2-1,indvars(1))*((q(i-1,ny+2-1,indvars(4))-&
                    0.5_wp*(q(i-1,ny+2-1,indvars(2))*q(i-1,ny+2-1,indvars(2))+&
                    q(i-1,ny+2-1,indvars(3))*q(i-1,ny+2-1,indvars(3))))))*q(i-1,ny+2-1,indvars(2))

d1_conv_et_dx_0_ip1nyp2m1k = (q(i+1,ny+2-1,indvars(1))*q(i+1,ny+2-1,indvars(4))+&
                    (param_float(3 + 5))*q(i+1,ny+2-1,indvars(1))*((q(i+1,ny+2-1,indvars(4))-&
                    0.5_wp*(q(i+1,ny+2-1,indvars(2))*q(i+1,ny+2-1,indvars(2))+&
                    q(i+1,ny+2-1,indvars(3))*q(i+1,ny+2-1,indvars(3))))))*q(i+1,ny+2-1,indvars(2))

d1_conv_et_dx_0_inyp2m1k = -&
          0.5_wp*d1_conv_et_dx_0_im1nyp2m1k+&
          0.5_wp*d1_conv_et_dx_0_ip1nyp2m1k

d1_conv_et_dx_0_inyp2m1k = d1_conv_et_dx_0_inyp2m1k*param_float(1)

d1_conv_et_dy_0_inyp2m1m1k = (q(i,ny+2-1-1,indvars(1))*q(i,ny+2-1-1,indvars(4))+&
                    (param_float(3 + 5))*q(i,ny+2-1-1,indvars(1))*((q(i,ny+2-1-1,indvars(4))-&
                    0.5_wp*(q(i,ny+2-1-1,indvars(2))*q(i,ny+2-1-1,indvars(2))+&
                    q(i,ny+2-1-1,indvars(3))*q(i,ny+2-1-1,indvars(3))))))*q(i,ny+2-1-1,indvars(3))

d1_conv_et_dy_0_inyp2m1p1k = (q(i,ny+2-1+1,indvars(1))*q(i,ny+2-1+1,indvars(4))+&
                    (param_float(3 + 5))*q(i,ny+2-1+1,indvars(1))*((q(i,ny+2-1+1,indvars(4))-&
                    0.5_wp*(q(i,ny+2-1+1,indvars(2))*q(i,ny+2-1+1,indvars(2))+&
                    q(i,ny+2-1+1,indvars(3))*q(i,ny+2-1+1,indvars(3))))))*q(i,ny+2-1+1,indvars(3))

d1_conv_et_dy_0_inyp2m1k = -&
          0.5_wp*d1_conv_et_dy_0_inyp2m1m1k+&
          0.5_wp*d1_conv_et_dy_0_inyp2m1p1k

d1_conv_et_dy_0_inyp2m1k = d1_conv_et_dy_0_inyp2m1k*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer None 1 None d(rho et)/dt *******
!                                                           
!***********************************************************


rhs(i,ny+2-1,indvars(4)) =   -  ( qst(i,ny+2-1,indvarsst(10))*(d1_conv_et_dx_0_inyp2m1k)+&
                    qst(i,ny+2-1,indvarsst(11))*(d1_conv_et_dy_0_inyp2m1k) ) 



!***********************************************************
!                                                           
! building source terms in RHS for layer None 1 None d(rho nut)/dt 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! deltaxI*([rho*u*nut-visc_t*sigmaI*deltaxI*({nut}_1x)]_1x)+deltayI*([rho*v*nut-visc_t*sigmaI*deltayI*({nut}_1y)]_1y)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d2_conv_nut_dxdx_0_0_im1nyp2m1k_im1m1nyp2m1k = q(i-1-1,ny+2-1,indvars(5))

d2_conv_nut_dxdx_0_0_im1nyp2m1k_im1p1nyp2m1k = q(i-1+1,ny+2-1,indvars(5))

d2_conv_nut_dxdx_0_0_im1nyp2m1k = -&
          0.5_wp*d2_conv_nut_dxdx_0_0_im1nyp2m1k_im1m1nyp2m1k+&
          0.5_wp*d2_conv_nut_dxdx_0_0_im1nyp2m1k_im1p1nyp2m1k

d2_conv_nut_dxdx_0_0_im1nyp2m1k = d2_conv_nut_dxdx_0_0_im1nyp2m1k*param_float(1)

d2_conv_nut_dxdx_0_0_ip1nyp2m1k_ip1m1nyp2m1k = q(i+1-1,ny+2-1,indvars(5))

d2_conv_nut_dxdx_0_0_ip1nyp2m1k_ip1p1nyp2m1k = q(i+1+1,ny+2-1,indvars(5))

d2_conv_nut_dxdx_0_0_ip1nyp2m1k = -&
          0.5_wp*d2_conv_nut_dxdx_0_0_ip1nyp2m1k_ip1m1nyp2m1k+&
          0.5_wp*d2_conv_nut_dxdx_0_0_ip1nyp2m1k_ip1p1nyp2m1k

d2_conv_nut_dxdx_0_0_ip1nyp2m1k = d2_conv_nut_dxdx_0_0_ip1nyp2m1k*param_float(1)

d1_conv_nut_dx_0_im1nyp2m1k = q(i-1,ny+2-1,indvars(1))*q(i-1,ny+2-1,indvars(2))*q(i-1,ny+2-1,indvars(5))-&
                    ((1+&
                    param_float(21 + 5))/(((q(i-1,ny+2-1,indvars(4))-&
                    0.5_wp*(q(i-1,ny+2-1,indvars(2))*q(i-1,ny+2-1,indvars(2))+&
                    q(i-1,ny+2-1,indvars(3))*q(i-1,ny+2-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(i-1,ny+2-1,indvars(4))-&
                    0.5_wp*(q(i-1,ny+2-1,indvars(2))*q(i-1,ny+2-1,indvars(2))+&
                    q(i-1,ny+2-1,indvars(3))*q(i-1,ny+2-1,indvars(3)))))/param_float(4 + 5)**1.5)*(1+&
                    ((q(i-1,ny+2-1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(i-1,ny+2-1,indvars(4))-&
                    0.5_wp*(q(i-1,ny+2-1,indvars(2))*q(i-1,ny+2-1,indvars(2))+&
                    q(i-1,ny+2-1,indvars(3))*q(i-1,ny+2-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(i-1,ny+2-1,indvars(4))-&
                    0.5_wp*(q(i-1,ny+2-1,indvars(2))*q(i-1,ny+2-1,indvars(2))+&
                    q(i-1,ny+2-1,indvars(3))*q(i-1,ny+2-1,indvars(3)))))/param_float(4 + 5)**1.5*q(i-1,ny+2-1,indvars(1)))**3/((q(i-1,ny+2-1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(i-1,ny+2-1,indvars(4))-&
                    0.5_wp*(q(i-1,ny+2-1,indvars(2))*q(i-1,ny+2-1,indvars(2))+&
                    q(i-1,ny+2-1,indvars(3))*q(i-1,ny+2-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(i-1,ny+2-1,indvars(4))-&
                    0.5_wp*(q(i-1,ny+2-1,indvars(2))*q(i-1,ny+2-1,indvars(2))+&
                    q(i-1,ny+2-1,indvars(3))*q(i-1,ny+2-1,indvars(3)))))/param_float(4 + 5)**1.5*q(i-1,ny+2-1,indvars(1)))**3+&
                    param_float(13 + 5)**3))*(q(i-1,ny+2-1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(i-1,ny+2-1,indvars(4))-&
                    0.5_wp*(q(i-1,ny+2-1,indvars(2))*q(i-1,ny+2-1,indvars(2))+&
                    q(i-1,ny+2-1,indvars(3))*q(i-1,ny+2-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(i-1,ny+2-1,indvars(4))-&
                    0.5_wp*(q(i-1,ny+2-1,indvars(2))*q(i-1,ny+2-1,indvars(2))+&
                    q(i-1,ny+2-1,indvars(3))*q(i-1,ny+2-1,indvars(3)))))/param_float(4 + 5)**1.5*q(i-1,ny+2-1,indvars(1))))*param_float(1 + 5)*param_float(18 + 5)*qst(i-1,ny+2-1,indvarsst(10))*(d2_conv_nut_dxdx_0_0_im1nyp2m1k)

d1_conv_nut_dx_0_ip1nyp2m1k = q(i+1,ny+2-1,indvars(1))*q(i+1,ny+2-1,indvars(2))*q(i+1,ny+2-1,indvars(5))-&
                    ((1+&
                    param_float(21 + 5))/(((q(i+1,ny+2-1,indvars(4))-&
                    0.5_wp*(q(i+1,ny+2-1,indvars(2))*q(i+1,ny+2-1,indvars(2))+&
                    q(i+1,ny+2-1,indvars(3))*q(i+1,ny+2-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(i+1,ny+2-1,indvars(4))-&
                    0.5_wp*(q(i+1,ny+2-1,indvars(2))*q(i+1,ny+2-1,indvars(2))+&
                    q(i+1,ny+2-1,indvars(3))*q(i+1,ny+2-1,indvars(3)))))/param_float(4 + 5)**1.5)*(1+&
                    ((q(i+1,ny+2-1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(i+1,ny+2-1,indvars(4))-&
                    0.5_wp*(q(i+1,ny+2-1,indvars(2))*q(i+1,ny+2-1,indvars(2))+&
                    q(i+1,ny+2-1,indvars(3))*q(i+1,ny+2-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(i+1,ny+2-1,indvars(4))-&
                    0.5_wp*(q(i+1,ny+2-1,indvars(2))*q(i+1,ny+2-1,indvars(2))+&
                    q(i+1,ny+2-1,indvars(3))*q(i+1,ny+2-1,indvars(3)))))/param_float(4 + 5)**1.5*q(i+1,ny+2-1,indvars(1)))**3/((q(i+1,ny+2-1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(i+1,ny+2-1,indvars(4))-&
                    0.5_wp*(q(i+1,ny+2-1,indvars(2))*q(i+1,ny+2-1,indvars(2))+&
                    q(i+1,ny+2-1,indvars(3))*q(i+1,ny+2-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(i+1,ny+2-1,indvars(4))-&
                    0.5_wp*(q(i+1,ny+2-1,indvars(2))*q(i+1,ny+2-1,indvars(2))+&
                    q(i+1,ny+2-1,indvars(3))*q(i+1,ny+2-1,indvars(3)))))/param_float(4 + 5)**1.5*q(i+1,ny+2-1,indvars(1)))**3+&
                    param_float(13 + 5)**3))*(q(i+1,ny+2-1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(i+1,ny+2-1,indvars(4))-&
                    0.5_wp*(q(i+1,ny+2-1,indvars(2))*q(i+1,ny+2-1,indvars(2))+&
                    q(i+1,ny+2-1,indvars(3))*q(i+1,ny+2-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(i+1,ny+2-1,indvars(4))-&
                    0.5_wp*(q(i+1,ny+2-1,indvars(2))*q(i+1,ny+2-1,indvars(2))+&
                    q(i+1,ny+2-1,indvars(3))*q(i+1,ny+2-1,indvars(3)))))/param_float(4 + 5)**1.5*q(i+1,ny+2-1,indvars(1))))*param_float(1 + 5)*param_float(18 + 5)*qst(i+1,ny+2-1,indvarsst(10))*(d2_conv_nut_dxdx_0_0_ip1nyp2m1k)

d1_conv_nut_dx_0_inyp2m1k = -&
          0.5_wp*d1_conv_nut_dx_0_im1nyp2m1k+&
          0.5_wp*d1_conv_nut_dx_0_ip1nyp2m1k

d1_conv_nut_dx_0_inyp2m1k = d1_conv_nut_dx_0_inyp2m1k*param_float(1)

d2_conv_nut_dydy_0_0_inyp2m1m1k_inyp2m1m1m1k = q(i,ny+2-1-1-1,indvars(5))

d2_conv_nut_dydy_0_0_inyp2m1m1k_inyp2m1m1p1k = q(i,ny+2-1-1+1,indvars(5))

d2_conv_nut_dydy_0_0_inyp2m1m1k = -&
          0.5_wp*d2_conv_nut_dydy_0_0_inyp2m1m1k_inyp2m1m1m1k+&
          0.5_wp*d2_conv_nut_dydy_0_0_inyp2m1m1k_inyp2m1m1p1k

d2_conv_nut_dydy_0_0_inyp2m1m1k = d2_conv_nut_dydy_0_0_inyp2m1m1k*param_float(2)

d2_conv_nut_dydy_0_0_inyp2m1p1k_inyp2m1p1p0k = q(i,ny+2-1+1+0,indvars(5))

d2_conv_nut_dydy_0_0_inyp2m1p1k_inyp2m1p1m1k = q(i,ny+2-1+1-1,indvars(5))

d2_conv_nut_dydy_0_0_inyp2m1p1k_inyp2m1p1m2k = q(i,ny+2-1+1-2,indvars(5))

d2_conv_nut_dydy_0_0_inyp2m1p1k = 1.5_wp*d2_conv_nut_dydy_0_0_inyp2m1p1k_inyp2m1p1p0k-&
          2.0_wp*d2_conv_nut_dydy_0_0_inyp2m1p1k_inyp2m1p1m1k+&
          0.5_wp*d2_conv_nut_dydy_0_0_inyp2m1p1k_inyp2m1p1m2k

d2_conv_nut_dydy_0_0_inyp2m1p1k = d2_conv_nut_dydy_0_0_inyp2m1p1k*param_float(2)

d1_conv_nut_dy_0_inyp2m1m1k = q(i,ny+2-1-1,indvars(1))*q(i,ny+2-1-1,indvars(3))*q(i,ny+2-1-1,indvars(5))-&
                    ((1+&
                    param_float(21 + 5))/(((q(i,ny+2-1-1,indvars(4))-&
                    0.5_wp*(q(i,ny+2-1-1,indvars(2))*q(i,ny+2-1-1,indvars(2))+&
                    q(i,ny+2-1-1,indvars(3))*q(i,ny+2-1-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(i,ny+2-1-1,indvars(4))-&
                    0.5_wp*(q(i,ny+2-1-1,indvars(2))*q(i,ny+2-1-1,indvars(2))+&
                    q(i,ny+2-1-1,indvars(3))*q(i,ny+2-1-1,indvars(3)))))/param_float(4 + 5)**1.5)*(1+&
                    ((q(i,ny+2-1-1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(i,ny+2-1-1,indvars(4))-&
                    0.5_wp*(q(i,ny+2-1-1,indvars(2))*q(i,ny+2-1-1,indvars(2))+&
                    q(i,ny+2-1-1,indvars(3))*q(i,ny+2-1-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(i,ny+2-1-1,indvars(4))-&
                    0.5_wp*(q(i,ny+2-1-1,indvars(2))*q(i,ny+2-1-1,indvars(2))+&
                    q(i,ny+2-1-1,indvars(3))*q(i,ny+2-1-1,indvars(3)))))/param_float(4 + 5)**1.5*q(i,ny+2-1-1,indvars(1)))**3/((q(i,ny+2-1-1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(i,ny+2-1-1,indvars(4))-&
                    0.5_wp*(q(i,ny+2-1-1,indvars(2))*q(i,ny+2-1-1,indvars(2))+&
                    q(i,ny+2-1-1,indvars(3))*q(i,ny+2-1-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(i,ny+2-1-1,indvars(4))-&
                    0.5_wp*(q(i,ny+2-1-1,indvars(2))*q(i,ny+2-1-1,indvars(2))+&
                    q(i,ny+2-1-1,indvars(3))*q(i,ny+2-1-1,indvars(3)))))/param_float(4 + 5)**1.5*q(i,ny+2-1-1,indvars(1)))**3+&
                    param_float(13 + 5)**3))*(q(i,ny+2-1-1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(i,ny+2-1-1,indvars(4))-&
                    0.5_wp*(q(i,ny+2-1-1,indvars(2))*q(i,ny+2-1-1,indvars(2))+&
                    q(i,ny+2-1-1,indvars(3))*q(i,ny+2-1-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(i,ny+2-1-1,indvars(4))-&
                    0.5_wp*(q(i,ny+2-1-1,indvars(2))*q(i,ny+2-1-1,indvars(2))+&
                    q(i,ny+2-1-1,indvars(3))*q(i,ny+2-1-1,indvars(3)))))/param_float(4 + 5)**1.5*q(i,ny+2-1-1,indvars(1))))*param_float(1 + 5)*param_float(18 + 5)*qst(i,ny+2-1-1,indvarsst(11))*(d2_conv_nut_dydy_0_0_inyp2m1m1k)

d1_conv_nut_dy_0_inyp2m1p1k = q(i,ny+2-1+1,indvars(1))*q(i,ny+2-1+1,indvars(3))*q(i,ny+2-1+1,indvars(5))-&
                    ((1+&
                    param_float(21 + 5))/(((q(i,ny+2-1+1,indvars(4))-&
                    0.5_wp*(q(i,ny+2-1+1,indvars(2))*q(i,ny+2-1+1,indvars(2))+&
                    q(i,ny+2-1+1,indvars(3))*q(i,ny+2-1+1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(i,ny+2-1+1,indvars(4))-&
                    0.5_wp*(q(i,ny+2-1+1,indvars(2))*q(i,ny+2-1+1,indvars(2))+&
                    q(i,ny+2-1+1,indvars(3))*q(i,ny+2-1+1,indvars(3)))))/param_float(4 + 5)**1.5)*(1+&
                    ((q(i,ny+2-1+1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(i,ny+2-1+1,indvars(4))-&
                    0.5_wp*(q(i,ny+2-1+1,indvars(2))*q(i,ny+2-1+1,indvars(2))+&
                    q(i,ny+2-1+1,indvars(3))*q(i,ny+2-1+1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(i,ny+2-1+1,indvars(4))-&
                    0.5_wp*(q(i,ny+2-1+1,indvars(2))*q(i,ny+2-1+1,indvars(2))+&
                    q(i,ny+2-1+1,indvars(3))*q(i,ny+2-1+1,indvars(3)))))/param_float(4 + 5)**1.5*q(i,ny+2-1+1,indvars(1)))**3/((q(i,ny+2-1+1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(i,ny+2-1+1,indvars(4))-&
                    0.5_wp*(q(i,ny+2-1+1,indvars(2))*q(i,ny+2-1+1,indvars(2))+&
                    q(i,ny+2-1+1,indvars(3))*q(i,ny+2-1+1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(i,ny+2-1+1,indvars(4))-&
                    0.5_wp*(q(i,ny+2-1+1,indvars(2))*q(i,ny+2-1+1,indvars(2))+&
                    q(i,ny+2-1+1,indvars(3))*q(i,ny+2-1+1,indvars(3)))))/param_float(4 + 5)**1.5*q(i,ny+2-1+1,indvars(1)))**3+&
                    param_float(13 + 5)**3))*(q(i,ny+2-1+1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(i,ny+2-1+1,indvars(4))-&
                    0.5_wp*(q(i,ny+2-1+1,indvars(2))*q(i,ny+2-1+1,indvars(2))+&
                    q(i,ny+2-1+1,indvars(3))*q(i,ny+2-1+1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(i,ny+2-1+1,indvars(4))-&
                    0.5_wp*(q(i,ny+2-1+1,indvars(2))*q(i,ny+2-1+1,indvars(2))+&
                    q(i,ny+2-1+1,indvars(3))*q(i,ny+2-1+1,indvars(3)))))/param_float(4 + 5)**1.5*q(i,ny+2-1+1,indvars(1))))*param_float(1 + 5)*param_float(18 + 5)*qst(i,ny+2-1+1,indvarsst(11))*(d2_conv_nut_dydy_0_0_inyp2m1p1k)

d1_conv_nut_dy_0_inyp2m1k = -&
          0.5_wp*d1_conv_nut_dy_0_inyp2m1m1k+&
          0.5_wp*d1_conv_nut_dy_0_inyp2m1p1k

d1_conv_nut_dy_0_inyp2m1k = d1_conv_nut_dy_0_inyp2m1k*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer None 1 None d(rho nut)/dt ******
!                                                           
!***********************************************************


rhs(i,ny+2-1,indvars(5)) =   -  ( qst(i,ny+2-1,indvarsst(10))*(d1_conv_nut_dx_0_inyp2m1k)+&
                    qst(i,ny+2-1,indvarsst(11))*(d1_conv_nut_dy_0_inyp2m1k) ) 

   enddo


!***********************************************************
!                                                           
! Start building layers for BC : None jmax None ************
!                                                           
!***********************************************************




!***********************************************************
!                                                           
! BC layer: None 1 None ************************************
!                                                           
!***********************************************************


 
      do i=idloop(1),idloop(2) 


!***********************************************************
!                                                           
! building source terms in RHS for layer None 1 None d(rho u)/dt 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! deltaxI*([-visc_t*(2.0_wp*deltaxI*({u}_1x)-2.0_wp/3.0_wp*(deltaxI*({u}_1x)+deltayI*({v}_1y)))]_1x)+deltayI*([-visc_t*(deltayI*({u}_1y)+deltaxI*({v}_1x))]_1y)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d2_dif_rhou_dxdx_0_0_im1nyp2m1k_im1m1nyp2m1k = q(i-1-1,ny+2-1,indvars(2))

d2_dif_rhou_dxdx_0_0_im1nyp2m1k_im1p1nyp2m1k = q(i-1+1,ny+2-1,indvars(2))

d2_dif_rhou_dxdx_0_0_im1nyp2m1k = -&
          0.5_wp*d2_dif_rhou_dxdx_0_0_im1nyp2m1k_im1m1nyp2m1k+&
          0.5_wp*d2_dif_rhou_dxdx_0_0_im1nyp2m1k_im1p1nyp2m1k

d2_dif_rhou_dxdx_0_0_im1nyp2m1k = d2_dif_rhou_dxdx_0_0_im1nyp2m1k*param_float(1)

d2_dif_rhou_dxdx_0_0_ip1nyp2m1k_ip1m1nyp2m1k = q(i+1-1,ny+2-1,indvars(2))

d2_dif_rhou_dxdx_0_0_ip1nyp2m1k_ip1p1nyp2m1k = q(i+1+1,ny+2-1,indvars(2))

d2_dif_rhou_dxdx_0_0_ip1nyp2m1k = -&
          0.5_wp*d2_dif_rhou_dxdx_0_0_ip1nyp2m1k_ip1m1nyp2m1k+&
          0.5_wp*d2_dif_rhou_dxdx_0_0_ip1nyp2m1k_ip1p1nyp2m1k

d2_dif_rhou_dxdx_0_0_ip1nyp2m1k = d2_dif_rhou_dxdx_0_0_ip1nyp2m1k*param_float(1)

d2_dif_rhou_dxdy_0_0_im1nyp2m1k_im1nyp2m1m1k = q(i-1,ny+2-1-1,indvars(3))

d2_dif_rhou_dxdy_0_0_im1nyp2m1k_im1nyp2m1p1k = q(i-1,ny+2-1+1,indvars(3))

d2_dif_rhou_dxdy_0_0_im1nyp2m1k = -&
          0.5_wp*d2_dif_rhou_dxdy_0_0_im1nyp2m1k_im1nyp2m1m1k+&
          0.5_wp*d2_dif_rhou_dxdy_0_0_im1nyp2m1k_im1nyp2m1p1k

d2_dif_rhou_dxdy_0_0_im1nyp2m1k = d2_dif_rhou_dxdy_0_0_im1nyp2m1k*param_float(2)

d2_dif_rhou_dxdy_0_0_ip1nyp2m1k_ip1nyp2m1m1k = q(i+1,ny+2-1-1,indvars(3))

d2_dif_rhou_dxdy_0_0_ip1nyp2m1k_ip1nyp2m1p1k = q(i+1,ny+2-1+1,indvars(3))

d2_dif_rhou_dxdy_0_0_ip1nyp2m1k = -&
          0.5_wp*d2_dif_rhou_dxdy_0_0_ip1nyp2m1k_ip1nyp2m1m1k+&
          0.5_wp*d2_dif_rhou_dxdy_0_0_ip1nyp2m1k_ip1nyp2m1p1k

d2_dif_rhou_dxdy_0_0_ip1nyp2m1k = d2_dif_rhou_dxdy_0_0_ip1nyp2m1k*param_float(2)

d1_dif_rhou_dx_0_im1nyp2m1k = -((1+&
                    param_float(21 + 5))/(((q(i-1,ny+2-1,indvars(4))-&
                    0.5_wp*(q(i-1,ny+2-1,indvars(2))*q(i-1,ny+2-1,indvars(2))+&
                    q(i-1,ny+2-1,indvars(3))*q(i-1,ny+2-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(i-1,ny+2-1,indvars(4))-&
                    0.5_wp*(q(i-1,ny+2-1,indvars(2))*q(i-1,ny+2-1,indvars(2))+&
                    q(i-1,ny+2-1,indvars(3))*q(i-1,ny+2-1,indvars(3)))))/param_float(4 + 5)**1.5)*(1+&
                    ((q(i-1,ny+2-1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(i-1,ny+2-1,indvars(4))-&
                    0.5_wp*(q(i-1,ny+2-1,indvars(2))*q(i-1,ny+2-1,indvars(2))+&
                    q(i-1,ny+2-1,indvars(3))*q(i-1,ny+2-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(i-1,ny+2-1,indvars(4))-&
                    0.5_wp*(q(i-1,ny+2-1,indvars(2))*q(i-1,ny+2-1,indvars(2))+&
                    q(i-1,ny+2-1,indvars(3))*q(i-1,ny+2-1,indvars(3)))))/param_float(4 + 5)**1.5*q(i-1,ny+2-1,indvars(1)))**3/((q(i-1,ny+2-1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(i-1,ny+2-1,indvars(4))-&
                    0.5_wp*(q(i-1,ny+2-1,indvars(2))*q(i-1,ny+2-1,indvars(2))+&
                    q(i-1,ny+2-1,indvars(3))*q(i-1,ny+2-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(i-1,ny+2-1,indvars(4))-&
                    0.5_wp*(q(i-1,ny+2-1,indvars(2))*q(i-1,ny+2-1,indvars(2))+&
                    q(i-1,ny+2-1,indvars(3))*q(i-1,ny+2-1,indvars(3)))))/param_float(4 + 5)**1.5*q(i-1,ny+2-1,indvars(1)))**3+&
                    param_float(13 + 5)**3))*(q(i-1,ny+2-1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(i-1,ny+2-1,indvars(4))-&
                    0.5_wp*(q(i-1,ny+2-1,indvars(2))*q(i-1,ny+2-1,indvars(2))+&
                    q(i-1,ny+2-1,indvars(3))*q(i-1,ny+2-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(i-1,ny+2-1,indvars(4))-&
                    0.5_wp*(q(i-1,ny+2-1,indvars(2))*q(i-1,ny+2-1,indvars(2))+&
                    q(i-1,ny+2-1,indvars(3))*q(i-1,ny+2-1,indvars(3)))))/param_float(4 + 5)**1.5*q(i-1,ny+2-1,indvars(1))))*param_float(1 + 5)*(2.0_wp*qst(i-1,ny+2-1,indvarsst(10))*(d2_dif_rhou_dxdx_0_0_im1nyp2m1k)-&
                    2.0_wp/3.0_wp*(qst(i-1,ny+2-1,indvarsst(10))*(d2_dif_rhou_dxdx_0_0_im1nyp2m1k)+&
                    qst(i-1,ny+2-1,indvarsst(11))*(d2_dif_rhou_dxdy_0_0_im1nyp2m1k)))

d1_dif_rhou_dx_0_ip1nyp2m1k = -((1+&
                    param_float(21 + 5))/(((q(i+1,ny+2-1,indvars(4))-&
                    0.5_wp*(q(i+1,ny+2-1,indvars(2))*q(i+1,ny+2-1,indvars(2))+&
                    q(i+1,ny+2-1,indvars(3))*q(i+1,ny+2-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(i+1,ny+2-1,indvars(4))-&
                    0.5_wp*(q(i+1,ny+2-1,indvars(2))*q(i+1,ny+2-1,indvars(2))+&
                    q(i+1,ny+2-1,indvars(3))*q(i+1,ny+2-1,indvars(3)))))/param_float(4 + 5)**1.5)*(1+&
                    ((q(i+1,ny+2-1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(i+1,ny+2-1,indvars(4))-&
                    0.5_wp*(q(i+1,ny+2-1,indvars(2))*q(i+1,ny+2-1,indvars(2))+&
                    q(i+1,ny+2-1,indvars(3))*q(i+1,ny+2-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(i+1,ny+2-1,indvars(4))-&
                    0.5_wp*(q(i+1,ny+2-1,indvars(2))*q(i+1,ny+2-1,indvars(2))+&
                    q(i+1,ny+2-1,indvars(3))*q(i+1,ny+2-1,indvars(3)))))/param_float(4 + 5)**1.5*q(i+1,ny+2-1,indvars(1)))**3/((q(i+1,ny+2-1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(i+1,ny+2-1,indvars(4))-&
                    0.5_wp*(q(i+1,ny+2-1,indvars(2))*q(i+1,ny+2-1,indvars(2))+&
                    q(i+1,ny+2-1,indvars(3))*q(i+1,ny+2-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(i+1,ny+2-1,indvars(4))-&
                    0.5_wp*(q(i+1,ny+2-1,indvars(2))*q(i+1,ny+2-1,indvars(2))+&
                    q(i+1,ny+2-1,indvars(3))*q(i+1,ny+2-1,indvars(3)))))/param_float(4 + 5)**1.5*q(i+1,ny+2-1,indvars(1)))**3+&
                    param_float(13 + 5)**3))*(q(i+1,ny+2-1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(i+1,ny+2-1,indvars(4))-&
                    0.5_wp*(q(i+1,ny+2-1,indvars(2))*q(i+1,ny+2-1,indvars(2))+&
                    q(i+1,ny+2-1,indvars(3))*q(i+1,ny+2-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(i+1,ny+2-1,indvars(4))-&
                    0.5_wp*(q(i+1,ny+2-1,indvars(2))*q(i+1,ny+2-1,indvars(2))+&
                    q(i+1,ny+2-1,indvars(3))*q(i+1,ny+2-1,indvars(3)))))/param_float(4 + 5)**1.5*q(i+1,ny+2-1,indvars(1))))*param_float(1 + 5)*(2.0_wp*qst(i+1,ny+2-1,indvarsst(10))*(d2_dif_rhou_dxdx_0_0_ip1nyp2m1k)-&
                    2.0_wp/3.0_wp*(qst(i+1,ny+2-1,indvarsst(10))*(d2_dif_rhou_dxdx_0_0_ip1nyp2m1k)+&
                    qst(i+1,ny+2-1,indvarsst(11))*(d2_dif_rhou_dxdy_0_0_ip1nyp2m1k)))

d1_dif_rhou_dx_0_inyp2m1k = -&
          0.5_wp*d1_dif_rhou_dx_0_im1nyp2m1k+&
          0.5_wp*d1_dif_rhou_dx_0_ip1nyp2m1k

d1_dif_rhou_dx_0_inyp2m1k = d1_dif_rhou_dx_0_inyp2m1k*param_float(1)

d2_dif_rhou_dydx_0_0_inyp2m1m1k_im1nyp2m1m1k = q(i-1,ny+2-1-1,indvars(3))

d2_dif_rhou_dydx_0_0_inyp2m1m1k_ip1nyp2m1m1k = q(i+1,ny+2-1-1,indvars(3))

d2_dif_rhou_dydx_0_0_inyp2m1m1k = -&
          0.5_wp*d2_dif_rhou_dydx_0_0_inyp2m1m1k_im1nyp2m1m1k+&
          0.5_wp*d2_dif_rhou_dydx_0_0_inyp2m1m1k_ip1nyp2m1m1k

d2_dif_rhou_dydx_0_0_inyp2m1m1k = d2_dif_rhou_dydx_0_0_inyp2m1m1k*param_float(1)

d2_dif_rhou_dydx_0_0_inyp2m1p1k_im1nyp2m1p1k = q(i-1,ny+2-1+1,indvars(3))

d2_dif_rhou_dydx_0_0_inyp2m1p1k_ip1nyp2m1p1k = q(i+1,ny+2-1+1,indvars(3))

d2_dif_rhou_dydx_0_0_inyp2m1p1k = -&
          0.5_wp*d2_dif_rhou_dydx_0_0_inyp2m1p1k_im1nyp2m1p1k+&
          0.5_wp*d2_dif_rhou_dydx_0_0_inyp2m1p1k_ip1nyp2m1p1k

d2_dif_rhou_dydx_0_0_inyp2m1p1k = d2_dif_rhou_dydx_0_0_inyp2m1p1k*param_float(1)

d2_dif_rhou_dydy_0_0_inyp2m1m1k_inyp2m1m1m1k = q(i,ny+2-1-1-1,indvars(2))

d2_dif_rhou_dydy_0_0_inyp2m1m1k_inyp2m1m1p1k = q(i,ny+2-1-1+1,indvars(2))

d2_dif_rhou_dydy_0_0_inyp2m1m1k = -&
          0.5_wp*d2_dif_rhou_dydy_0_0_inyp2m1m1k_inyp2m1m1m1k+&
          0.5_wp*d2_dif_rhou_dydy_0_0_inyp2m1m1k_inyp2m1m1p1k

d2_dif_rhou_dydy_0_0_inyp2m1m1k = d2_dif_rhou_dydy_0_0_inyp2m1m1k*param_float(2)

d2_dif_rhou_dydy_0_0_inyp2m1p1k_inyp2m1p1p0k = q(i,ny+2-1+1+0,indvars(2))

d2_dif_rhou_dydy_0_0_inyp2m1p1k_inyp2m1p1m1k = q(i,ny+2-1+1-1,indvars(2))

d2_dif_rhou_dydy_0_0_inyp2m1p1k_inyp2m1p1m2k = q(i,ny+2-1+1-2,indvars(2))

d2_dif_rhou_dydy_0_0_inyp2m1p1k = 1.5_wp*d2_dif_rhou_dydy_0_0_inyp2m1p1k_inyp2m1p1p0k-&
          2.0_wp*d2_dif_rhou_dydy_0_0_inyp2m1p1k_inyp2m1p1m1k+&
          0.5_wp*d2_dif_rhou_dydy_0_0_inyp2m1p1k_inyp2m1p1m2k

d2_dif_rhou_dydy_0_0_inyp2m1p1k = d2_dif_rhou_dydy_0_0_inyp2m1p1k*param_float(2)

d1_dif_rhou_dy_0_inyp2m1m1k = -((1+&
                    param_float(21 + 5))/(((q(i,ny+2-1-1,indvars(4))-&
                    0.5_wp*(q(i,ny+2-1-1,indvars(2))*q(i,ny+2-1-1,indvars(2))+&
                    q(i,ny+2-1-1,indvars(3))*q(i,ny+2-1-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(i,ny+2-1-1,indvars(4))-&
                    0.5_wp*(q(i,ny+2-1-1,indvars(2))*q(i,ny+2-1-1,indvars(2))+&
                    q(i,ny+2-1-1,indvars(3))*q(i,ny+2-1-1,indvars(3)))))/param_float(4 + 5)**1.5)*(1+&
                    ((q(i,ny+2-1-1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(i,ny+2-1-1,indvars(4))-&
                    0.5_wp*(q(i,ny+2-1-1,indvars(2))*q(i,ny+2-1-1,indvars(2))+&
                    q(i,ny+2-1-1,indvars(3))*q(i,ny+2-1-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(i,ny+2-1-1,indvars(4))-&
                    0.5_wp*(q(i,ny+2-1-1,indvars(2))*q(i,ny+2-1-1,indvars(2))+&
                    q(i,ny+2-1-1,indvars(3))*q(i,ny+2-1-1,indvars(3)))))/param_float(4 + 5)**1.5*q(i,ny+2-1-1,indvars(1)))**3/((q(i,ny+2-1-1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(i,ny+2-1-1,indvars(4))-&
                    0.5_wp*(q(i,ny+2-1-1,indvars(2))*q(i,ny+2-1-1,indvars(2))+&
                    q(i,ny+2-1-1,indvars(3))*q(i,ny+2-1-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(i,ny+2-1-1,indvars(4))-&
                    0.5_wp*(q(i,ny+2-1-1,indvars(2))*q(i,ny+2-1-1,indvars(2))+&
                    q(i,ny+2-1-1,indvars(3))*q(i,ny+2-1-1,indvars(3)))))/param_float(4 + 5)**1.5*q(i,ny+2-1-1,indvars(1)))**3+&
                    param_float(13 + 5)**3))*(q(i,ny+2-1-1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(i,ny+2-1-1,indvars(4))-&
                    0.5_wp*(q(i,ny+2-1-1,indvars(2))*q(i,ny+2-1-1,indvars(2))+&
                    q(i,ny+2-1-1,indvars(3))*q(i,ny+2-1-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(i,ny+2-1-1,indvars(4))-&
                    0.5_wp*(q(i,ny+2-1-1,indvars(2))*q(i,ny+2-1-1,indvars(2))+&
                    q(i,ny+2-1-1,indvars(3))*q(i,ny+2-1-1,indvars(3)))))/param_float(4 + 5)**1.5*q(i,ny+2-1-1,indvars(1))))*param_float(1 + 5)*(qst(i,ny+2-1-1,indvarsst(11))*(d2_dif_rhou_dydy_0_0_inyp2m1m1k)+&
                    qst(i,ny+2-1-1,indvarsst(10))*(d2_dif_rhou_dydx_0_0_inyp2m1m1k))

d1_dif_rhou_dy_0_inyp2m1p1k = -((1+&
                    param_float(21 + 5))/(((q(i,ny+2-1+1,indvars(4))-&
                    0.5_wp*(q(i,ny+2-1+1,indvars(2))*q(i,ny+2-1+1,indvars(2))+&
                    q(i,ny+2-1+1,indvars(3))*q(i,ny+2-1+1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(i,ny+2-1+1,indvars(4))-&
                    0.5_wp*(q(i,ny+2-1+1,indvars(2))*q(i,ny+2-1+1,indvars(2))+&
                    q(i,ny+2-1+1,indvars(3))*q(i,ny+2-1+1,indvars(3)))))/param_float(4 + 5)**1.5)*(1+&
                    ((q(i,ny+2-1+1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(i,ny+2-1+1,indvars(4))-&
                    0.5_wp*(q(i,ny+2-1+1,indvars(2))*q(i,ny+2-1+1,indvars(2))+&
                    q(i,ny+2-1+1,indvars(3))*q(i,ny+2-1+1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(i,ny+2-1+1,indvars(4))-&
                    0.5_wp*(q(i,ny+2-1+1,indvars(2))*q(i,ny+2-1+1,indvars(2))+&
                    q(i,ny+2-1+1,indvars(3))*q(i,ny+2-1+1,indvars(3)))))/param_float(4 + 5)**1.5*q(i,ny+2-1+1,indvars(1)))**3/((q(i,ny+2-1+1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(i,ny+2-1+1,indvars(4))-&
                    0.5_wp*(q(i,ny+2-1+1,indvars(2))*q(i,ny+2-1+1,indvars(2))+&
                    q(i,ny+2-1+1,indvars(3))*q(i,ny+2-1+1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(i,ny+2-1+1,indvars(4))-&
                    0.5_wp*(q(i,ny+2-1+1,indvars(2))*q(i,ny+2-1+1,indvars(2))+&
                    q(i,ny+2-1+1,indvars(3))*q(i,ny+2-1+1,indvars(3)))))/param_float(4 + 5)**1.5*q(i,ny+2-1+1,indvars(1)))**3+&
                    param_float(13 + 5)**3))*(q(i,ny+2-1+1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(i,ny+2-1+1,indvars(4))-&
                    0.5_wp*(q(i,ny+2-1+1,indvars(2))*q(i,ny+2-1+1,indvars(2))+&
                    q(i,ny+2-1+1,indvars(3))*q(i,ny+2-1+1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(i,ny+2-1+1,indvars(4))-&
                    0.5_wp*(q(i,ny+2-1+1,indvars(2))*q(i,ny+2-1+1,indvars(2))+&
                    q(i,ny+2-1+1,indvars(3))*q(i,ny+2-1+1,indvars(3)))))/param_float(4 + 5)**1.5*q(i,ny+2-1+1,indvars(1))))*param_float(1 + 5)*(qst(i,ny+2-1+1,indvarsst(11))*(d2_dif_rhou_dydy_0_0_inyp2m1p1k)+&
                    qst(i,ny+2-1+1,indvarsst(10))*(d2_dif_rhou_dydx_0_0_inyp2m1p1k))

d1_dif_rhou_dy_0_inyp2m1k = -&
          0.5_wp*d1_dif_rhou_dy_0_inyp2m1m1k+&
          0.5_wp*d1_dif_rhou_dy_0_inyp2m1p1k

d1_dif_rhou_dy_0_inyp2m1k = d1_dif_rhou_dy_0_inyp2m1k*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer None 1 None d(rho u)/dt ********
!                                                           
!***********************************************************


rhs(i,ny+2-1,indvars(2)) = rhs(i,ny+2-1,indvars(2))  -  ( qst(i,ny+2-1,indvarsst(10))*(d1_dif_rhou_dx_0_inyp2m1k)+&
                    qst(i,ny+2-1,indvarsst(11))*(d1_dif_rhou_dy_0_inyp2m1k) ) 



!***********************************************************
!                                                           
! building source terms in RHS for layer None 1 None d(rho v)/dt 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! deltaxI*([-visc_t*(deltayI*({u}_1y)+deltaxI*({v}_1x))]_1x)+deltayI*([-visc_t*(2.0_wp*deltayI*({v}_1y)-2.0_wp/3.0_wp*(deltaxI*({u}_1x)+deltayI*({v}_1y)))]_1y)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d2_dif_rhov_dxdx_0_0_im1nyp2m1k_im1m1nyp2m1k = q(i-1-1,ny+2-1,indvars(3))

d2_dif_rhov_dxdx_0_0_im1nyp2m1k_im1p1nyp2m1k = q(i-1+1,ny+2-1,indvars(3))

d2_dif_rhov_dxdx_0_0_im1nyp2m1k = -&
          0.5_wp*d2_dif_rhov_dxdx_0_0_im1nyp2m1k_im1m1nyp2m1k+&
          0.5_wp*d2_dif_rhov_dxdx_0_0_im1nyp2m1k_im1p1nyp2m1k

d2_dif_rhov_dxdx_0_0_im1nyp2m1k = d2_dif_rhov_dxdx_0_0_im1nyp2m1k*param_float(1)

d2_dif_rhov_dxdx_0_0_ip1nyp2m1k_ip1m1nyp2m1k = q(i+1-1,ny+2-1,indvars(3))

d2_dif_rhov_dxdx_0_0_ip1nyp2m1k_ip1p1nyp2m1k = q(i+1+1,ny+2-1,indvars(3))

d2_dif_rhov_dxdx_0_0_ip1nyp2m1k = -&
          0.5_wp*d2_dif_rhov_dxdx_0_0_ip1nyp2m1k_ip1m1nyp2m1k+&
          0.5_wp*d2_dif_rhov_dxdx_0_0_ip1nyp2m1k_ip1p1nyp2m1k

d2_dif_rhov_dxdx_0_0_ip1nyp2m1k = d2_dif_rhov_dxdx_0_0_ip1nyp2m1k*param_float(1)

d2_dif_rhov_dxdy_0_0_im1nyp2m1k_im1nyp2m1m1k = q(i-1,ny+2-1-1,indvars(2))

d2_dif_rhov_dxdy_0_0_im1nyp2m1k_im1nyp2m1p1k = q(i-1,ny+2-1+1,indvars(2))

d2_dif_rhov_dxdy_0_0_im1nyp2m1k = -&
          0.5_wp*d2_dif_rhov_dxdy_0_0_im1nyp2m1k_im1nyp2m1m1k+&
          0.5_wp*d2_dif_rhov_dxdy_0_0_im1nyp2m1k_im1nyp2m1p1k

d2_dif_rhov_dxdy_0_0_im1nyp2m1k = d2_dif_rhov_dxdy_0_0_im1nyp2m1k*param_float(2)

d2_dif_rhov_dxdy_0_0_ip1nyp2m1k_ip1nyp2m1m1k = q(i+1,ny+2-1-1,indvars(2))

d2_dif_rhov_dxdy_0_0_ip1nyp2m1k_ip1nyp2m1p1k = q(i+1,ny+2-1+1,indvars(2))

d2_dif_rhov_dxdy_0_0_ip1nyp2m1k = -&
          0.5_wp*d2_dif_rhov_dxdy_0_0_ip1nyp2m1k_ip1nyp2m1m1k+&
          0.5_wp*d2_dif_rhov_dxdy_0_0_ip1nyp2m1k_ip1nyp2m1p1k

d2_dif_rhov_dxdy_0_0_ip1nyp2m1k = d2_dif_rhov_dxdy_0_0_ip1nyp2m1k*param_float(2)

d1_dif_rhov_dx_0_im1nyp2m1k = -((1+&
                    param_float(21 + 5))/(((q(i-1,ny+2-1,indvars(4))-&
                    0.5_wp*(q(i-1,ny+2-1,indvars(2))*q(i-1,ny+2-1,indvars(2))+&
                    q(i-1,ny+2-1,indvars(3))*q(i-1,ny+2-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(i-1,ny+2-1,indvars(4))-&
                    0.5_wp*(q(i-1,ny+2-1,indvars(2))*q(i-1,ny+2-1,indvars(2))+&
                    q(i-1,ny+2-1,indvars(3))*q(i-1,ny+2-1,indvars(3)))))/param_float(4 + 5)**1.5)*(1+&
                    ((q(i-1,ny+2-1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(i-1,ny+2-1,indvars(4))-&
                    0.5_wp*(q(i-1,ny+2-1,indvars(2))*q(i-1,ny+2-1,indvars(2))+&
                    q(i-1,ny+2-1,indvars(3))*q(i-1,ny+2-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(i-1,ny+2-1,indvars(4))-&
                    0.5_wp*(q(i-1,ny+2-1,indvars(2))*q(i-1,ny+2-1,indvars(2))+&
                    q(i-1,ny+2-1,indvars(3))*q(i-1,ny+2-1,indvars(3)))))/param_float(4 + 5)**1.5*q(i-1,ny+2-1,indvars(1)))**3/((q(i-1,ny+2-1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(i-1,ny+2-1,indvars(4))-&
                    0.5_wp*(q(i-1,ny+2-1,indvars(2))*q(i-1,ny+2-1,indvars(2))+&
                    q(i-1,ny+2-1,indvars(3))*q(i-1,ny+2-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(i-1,ny+2-1,indvars(4))-&
                    0.5_wp*(q(i-1,ny+2-1,indvars(2))*q(i-1,ny+2-1,indvars(2))+&
                    q(i-1,ny+2-1,indvars(3))*q(i-1,ny+2-1,indvars(3)))))/param_float(4 + 5)**1.5*q(i-1,ny+2-1,indvars(1)))**3+&
                    param_float(13 + 5)**3))*(q(i-1,ny+2-1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(i-1,ny+2-1,indvars(4))-&
                    0.5_wp*(q(i-1,ny+2-1,indvars(2))*q(i-1,ny+2-1,indvars(2))+&
                    q(i-1,ny+2-1,indvars(3))*q(i-1,ny+2-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(i-1,ny+2-1,indvars(4))-&
                    0.5_wp*(q(i-1,ny+2-1,indvars(2))*q(i-1,ny+2-1,indvars(2))+&
                    q(i-1,ny+2-1,indvars(3))*q(i-1,ny+2-1,indvars(3)))))/param_float(4 + 5)**1.5*q(i-1,ny+2-1,indvars(1))))*param_float(1 + 5)*(qst(i-1,ny+2-1,indvarsst(11))*(d2_dif_rhov_dxdy_0_0_im1nyp2m1k)+&
                    qst(i-1,ny+2-1,indvarsst(10))*(d2_dif_rhov_dxdx_0_0_im1nyp2m1k))

d1_dif_rhov_dx_0_ip1nyp2m1k = -((1+&
                    param_float(21 + 5))/(((q(i+1,ny+2-1,indvars(4))-&
                    0.5_wp*(q(i+1,ny+2-1,indvars(2))*q(i+1,ny+2-1,indvars(2))+&
                    q(i+1,ny+2-1,indvars(3))*q(i+1,ny+2-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(i+1,ny+2-1,indvars(4))-&
                    0.5_wp*(q(i+1,ny+2-1,indvars(2))*q(i+1,ny+2-1,indvars(2))+&
                    q(i+1,ny+2-1,indvars(3))*q(i+1,ny+2-1,indvars(3)))))/param_float(4 + 5)**1.5)*(1+&
                    ((q(i+1,ny+2-1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(i+1,ny+2-1,indvars(4))-&
                    0.5_wp*(q(i+1,ny+2-1,indvars(2))*q(i+1,ny+2-1,indvars(2))+&
                    q(i+1,ny+2-1,indvars(3))*q(i+1,ny+2-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(i+1,ny+2-1,indvars(4))-&
                    0.5_wp*(q(i+1,ny+2-1,indvars(2))*q(i+1,ny+2-1,indvars(2))+&
                    q(i+1,ny+2-1,indvars(3))*q(i+1,ny+2-1,indvars(3)))))/param_float(4 + 5)**1.5*q(i+1,ny+2-1,indvars(1)))**3/((q(i+1,ny+2-1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(i+1,ny+2-1,indvars(4))-&
                    0.5_wp*(q(i+1,ny+2-1,indvars(2))*q(i+1,ny+2-1,indvars(2))+&
                    q(i+1,ny+2-1,indvars(3))*q(i+1,ny+2-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(i+1,ny+2-1,indvars(4))-&
                    0.5_wp*(q(i+1,ny+2-1,indvars(2))*q(i+1,ny+2-1,indvars(2))+&
                    q(i+1,ny+2-1,indvars(3))*q(i+1,ny+2-1,indvars(3)))))/param_float(4 + 5)**1.5*q(i+1,ny+2-1,indvars(1)))**3+&
                    param_float(13 + 5)**3))*(q(i+1,ny+2-1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(i+1,ny+2-1,indvars(4))-&
                    0.5_wp*(q(i+1,ny+2-1,indvars(2))*q(i+1,ny+2-1,indvars(2))+&
                    q(i+1,ny+2-1,indvars(3))*q(i+1,ny+2-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(i+1,ny+2-1,indvars(4))-&
                    0.5_wp*(q(i+1,ny+2-1,indvars(2))*q(i+1,ny+2-1,indvars(2))+&
                    q(i+1,ny+2-1,indvars(3))*q(i+1,ny+2-1,indvars(3)))))/param_float(4 + 5)**1.5*q(i+1,ny+2-1,indvars(1))))*param_float(1 + 5)*(qst(i+1,ny+2-1,indvarsst(11))*(d2_dif_rhov_dxdy_0_0_ip1nyp2m1k)+&
                    qst(i+1,ny+2-1,indvarsst(10))*(d2_dif_rhov_dxdx_0_0_ip1nyp2m1k))

d1_dif_rhov_dx_0_inyp2m1k = -&
          0.5_wp*d1_dif_rhov_dx_0_im1nyp2m1k+&
          0.5_wp*d1_dif_rhov_dx_0_ip1nyp2m1k

d1_dif_rhov_dx_0_inyp2m1k = d1_dif_rhov_dx_0_inyp2m1k*param_float(1)

d2_dif_rhov_dydx_0_0_inyp2m1m1k_im1nyp2m1m1k = q(i-1,ny+2-1-1,indvars(2))

d2_dif_rhov_dydx_0_0_inyp2m1m1k_ip1nyp2m1m1k = q(i+1,ny+2-1-1,indvars(2))

d2_dif_rhov_dydx_0_0_inyp2m1m1k = -&
          0.5_wp*d2_dif_rhov_dydx_0_0_inyp2m1m1k_im1nyp2m1m1k+&
          0.5_wp*d2_dif_rhov_dydx_0_0_inyp2m1m1k_ip1nyp2m1m1k

d2_dif_rhov_dydx_0_0_inyp2m1m1k = d2_dif_rhov_dydx_0_0_inyp2m1m1k*param_float(1)

d2_dif_rhov_dydx_0_0_inyp2m1p1k_im1nyp2m1p1k = q(i-1,ny+2-1+1,indvars(2))

d2_dif_rhov_dydx_0_0_inyp2m1p1k_ip1nyp2m1p1k = q(i+1,ny+2-1+1,indvars(2))

d2_dif_rhov_dydx_0_0_inyp2m1p1k = -&
          0.5_wp*d2_dif_rhov_dydx_0_0_inyp2m1p1k_im1nyp2m1p1k+&
          0.5_wp*d2_dif_rhov_dydx_0_0_inyp2m1p1k_ip1nyp2m1p1k

d2_dif_rhov_dydx_0_0_inyp2m1p1k = d2_dif_rhov_dydx_0_0_inyp2m1p1k*param_float(1)

d2_dif_rhov_dydy_0_0_inyp2m1m1k_inyp2m1m1m1k = q(i,ny+2-1-1-1,indvars(3))

d2_dif_rhov_dydy_0_0_inyp2m1m1k_inyp2m1m1p1k = q(i,ny+2-1-1+1,indvars(3))

d2_dif_rhov_dydy_0_0_inyp2m1m1k = -&
          0.5_wp*d2_dif_rhov_dydy_0_0_inyp2m1m1k_inyp2m1m1m1k+&
          0.5_wp*d2_dif_rhov_dydy_0_0_inyp2m1m1k_inyp2m1m1p1k

d2_dif_rhov_dydy_0_0_inyp2m1m1k = d2_dif_rhov_dydy_0_0_inyp2m1m1k*param_float(2)

d2_dif_rhov_dydy_0_0_inyp2m1p1k_inyp2m1p1p0k = q(i,ny+2-1+1+0,indvars(3))

d2_dif_rhov_dydy_0_0_inyp2m1p1k_inyp2m1p1m1k = q(i,ny+2-1+1-1,indvars(3))

d2_dif_rhov_dydy_0_0_inyp2m1p1k_inyp2m1p1m2k = q(i,ny+2-1+1-2,indvars(3))

d2_dif_rhov_dydy_0_0_inyp2m1p1k = 1.5_wp*d2_dif_rhov_dydy_0_0_inyp2m1p1k_inyp2m1p1p0k-&
          2.0_wp*d2_dif_rhov_dydy_0_0_inyp2m1p1k_inyp2m1p1m1k+&
          0.5_wp*d2_dif_rhov_dydy_0_0_inyp2m1p1k_inyp2m1p1m2k

d2_dif_rhov_dydy_0_0_inyp2m1p1k = d2_dif_rhov_dydy_0_0_inyp2m1p1k*param_float(2)

d1_dif_rhov_dy_0_inyp2m1m1k = -((1+&
                    param_float(21 + 5))/(((q(i,ny+2-1-1,indvars(4))-&
                    0.5_wp*(q(i,ny+2-1-1,indvars(2))*q(i,ny+2-1-1,indvars(2))+&
                    q(i,ny+2-1-1,indvars(3))*q(i,ny+2-1-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(i,ny+2-1-1,indvars(4))-&
                    0.5_wp*(q(i,ny+2-1-1,indvars(2))*q(i,ny+2-1-1,indvars(2))+&
                    q(i,ny+2-1-1,indvars(3))*q(i,ny+2-1-1,indvars(3)))))/param_float(4 + 5)**1.5)*(1+&
                    ((q(i,ny+2-1-1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(i,ny+2-1-1,indvars(4))-&
                    0.5_wp*(q(i,ny+2-1-1,indvars(2))*q(i,ny+2-1-1,indvars(2))+&
                    q(i,ny+2-1-1,indvars(3))*q(i,ny+2-1-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(i,ny+2-1-1,indvars(4))-&
                    0.5_wp*(q(i,ny+2-1-1,indvars(2))*q(i,ny+2-1-1,indvars(2))+&
                    q(i,ny+2-1-1,indvars(3))*q(i,ny+2-1-1,indvars(3)))))/param_float(4 + 5)**1.5*q(i,ny+2-1-1,indvars(1)))**3/((q(i,ny+2-1-1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(i,ny+2-1-1,indvars(4))-&
                    0.5_wp*(q(i,ny+2-1-1,indvars(2))*q(i,ny+2-1-1,indvars(2))+&
                    q(i,ny+2-1-1,indvars(3))*q(i,ny+2-1-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(i,ny+2-1-1,indvars(4))-&
                    0.5_wp*(q(i,ny+2-1-1,indvars(2))*q(i,ny+2-1-1,indvars(2))+&
                    q(i,ny+2-1-1,indvars(3))*q(i,ny+2-1-1,indvars(3)))))/param_float(4 + 5)**1.5*q(i,ny+2-1-1,indvars(1)))**3+&
                    param_float(13 + 5)**3))*(q(i,ny+2-1-1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(i,ny+2-1-1,indvars(4))-&
                    0.5_wp*(q(i,ny+2-1-1,indvars(2))*q(i,ny+2-1-1,indvars(2))+&
                    q(i,ny+2-1-1,indvars(3))*q(i,ny+2-1-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(i,ny+2-1-1,indvars(4))-&
                    0.5_wp*(q(i,ny+2-1-1,indvars(2))*q(i,ny+2-1-1,indvars(2))+&
                    q(i,ny+2-1-1,indvars(3))*q(i,ny+2-1-1,indvars(3)))))/param_float(4 + 5)**1.5*q(i,ny+2-1-1,indvars(1))))*param_float(1 + 5)*(2.0_wp*qst(i,ny+2-1-1,indvarsst(11))*(d2_dif_rhov_dydy_0_0_inyp2m1m1k)-&
                    2.0_wp/3.0_wp*(qst(i,ny+2-1-1,indvarsst(10))*(d2_dif_rhov_dydx_0_0_inyp2m1m1k)+&
                    qst(i,ny+2-1-1,indvarsst(11))*(d2_dif_rhov_dydy_0_0_inyp2m1m1k)))

d1_dif_rhov_dy_0_inyp2m1p1k = -((1+&
                    param_float(21 + 5))/(((q(i,ny+2-1+1,indvars(4))-&
                    0.5_wp*(q(i,ny+2-1+1,indvars(2))*q(i,ny+2-1+1,indvars(2))+&
                    q(i,ny+2-1+1,indvars(3))*q(i,ny+2-1+1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(i,ny+2-1+1,indvars(4))-&
                    0.5_wp*(q(i,ny+2-1+1,indvars(2))*q(i,ny+2-1+1,indvars(2))+&
                    q(i,ny+2-1+1,indvars(3))*q(i,ny+2-1+1,indvars(3)))))/param_float(4 + 5)**1.5)*(1+&
                    ((q(i,ny+2-1+1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(i,ny+2-1+1,indvars(4))-&
                    0.5_wp*(q(i,ny+2-1+1,indvars(2))*q(i,ny+2-1+1,indvars(2))+&
                    q(i,ny+2-1+1,indvars(3))*q(i,ny+2-1+1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(i,ny+2-1+1,indvars(4))-&
                    0.5_wp*(q(i,ny+2-1+1,indvars(2))*q(i,ny+2-1+1,indvars(2))+&
                    q(i,ny+2-1+1,indvars(3))*q(i,ny+2-1+1,indvars(3)))))/param_float(4 + 5)**1.5*q(i,ny+2-1+1,indvars(1)))**3/((q(i,ny+2-1+1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(i,ny+2-1+1,indvars(4))-&
                    0.5_wp*(q(i,ny+2-1+1,indvars(2))*q(i,ny+2-1+1,indvars(2))+&
                    q(i,ny+2-1+1,indvars(3))*q(i,ny+2-1+1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(i,ny+2-1+1,indvars(4))-&
                    0.5_wp*(q(i,ny+2-1+1,indvars(2))*q(i,ny+2-1+1,indvars(2))+&
                    q(i,ny+2-1+1,indvars(3))*q(i,ny+2-1+1,indvars(3)))))/param_float(4 + 5)**1.5*q(i,ny+2-1+1,indvars(1)))**3+&
                    param_float(13 + 5)**3))*(q(i,ny+2-1+1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(i,ny+2-1+1,indvars(4))-&
                    0.5_wp*(q(i,ny+2-1+1,indvars(2))*q(i,ny+2-1+1,indvars(2))+&
                    q(i,ny+2-1+1,indvars(3))*q(i,ny+2-1+1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(i,ny+2-1+1,indvars(4))-&
                    0.5_wp*(q(i,ny+2-1+1,indvars(2))*q(i,ny+2-1+1,indvars(2))+&
                    q(i,ny+2-1+1,indvars(3))*q(i,ny+2-1+1,indvars(3)))))/param_float(4 + 5)**1.5*q(i,ny+2-1+1,indvars(1))))*param_float(1 + 5)*(2.0_wp*qst(i,ny+2-1+1,indvarsst(11))*(d2_dif_rhov_dydy_0_0_inyp2m1p1k)-&
                    2.0_wp/3.0_wp*(qst(i,ny+2-1+1,indvarsst(10))*(d2_dif_rhov_dydx_0_0_inyp2m1p1k)+&
                    qst(i,ny+2-1+1,indvarsst(11))*(d2_dif_rhov_dydy_0_0_inyp2m1p1k)))

d1_dif_rhov_dy_0_inyp2m1k = -&
          0.5_wp*d1_dif_rhov_dy_0_inyp2m1m1k+&
          0.5_wp*d1_dif_rhov_dy_0_inyp2m1p1k

d1_dif_rhov_dy_0_inyp2m1k = d1_dif_rhov_dy_0_inyp2m1k*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer None 1 None d(rho v)/dt ********
!                                                           
!***********************************************************


rhs(i,ny+2-1,indvars(3)) = rhs(i,ny+2-1,indvars(3))  -  ( qst(i,ny+2-1,indvarsst(10))*(d1_dif_rhov_dx_0_inyp2m1k)+&
                    qst(i,ny+2-1,indvarsst(11))*(d1_dif_rhov_dy_0_inyp2m1k) ) 



!***********************************************************
!                                                           
! building source terms in RHS for layer None 1 None d(rho et)/dt 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! deltaxI*([-kappa*deltaxI*({T}_1x)-u*(visc_t*(2.0_wp*deltaxI*({u}_1x)-2.0_wp/3.0_wp*(deltaxI*({u}_1x)+deltayI*({v}_1y))))-v*(visc_t*(deltayI*({u}_1y)+deltaxI*({v}_1x)))]_1x)+deltayI*([-kappa*deltayI*({T}_1y)-u*(visc_t*(deltayI*({u}_1y)+deltaxI*({v}_1x)))-v*(visc_t*(2.0_wp*deltayI*({v}_1y)-2.0_wp/3.0_wp*(deltaxI*({u}_1x)+deltayI*({v}_1y))))]_1y)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d2_dif_et_dxdx_0_0_im1nyp2m1k_im1m1nyp2m1k = ((q(i-1-1,ny+2-1,indvars(4))-&
                    0.5_wp*(q(i-1-1,ny+2-1,indvars(2))*q(i-1-1,ny+2-1,indvars(2))+&
                    q(i-1-1,ny+2-1,indvars(3))*q(i-1-1,ny+2-1,indvars(3)))))/param_float(4 + 5)

d2_dif_et_dxdx_0_0_im1nyp2m1k_im1p1nyp2m1k = ((q(i-1+1,ny+2-1,indvars(4))-&
                    0.5_wp*(q(i-1+1,ny+2-1,indvars(2))*q(i-1+1,ny+2-1,indvars(2))+&
                    q(i-1+1,ny+2-1,indvars(3))*q(i-1+1,ny+2-1,indvars(3)))))/param_float(4 + 5)

d2_dif_et_dxdx_0_0_im1nyp2m1k = -&
          0.5_wp*d2_dif_et_dxdx_0_0_im1nyp2m1k_im1m1nyp2m1k+&
          0.5_wp*d2_dif_et_dxdx_0_0_im1nyp2m1k_im1p1nyp2m1k

d2_dif_et_dxdx_0_0_im1nyp2m1k = d2_dif_et_dxdx_0_0_im1nyp2m1k*param_float(1)

d2_dif_et_dxdx_0_0_ip1nyp2m1k_ip1m1nyp2m1k = ((q(i+1-1,ny+2-1,indvars(4))-&
                    0.5_wp*(q(i+1-1,ny+2-1,indvars(2))*q(i+1-1,ny+2-1,indvars(2))+&
                    q(i+1-1,ny+2-1,indvars(3))*q(i+1-1,ny+2-1,indvars(3)))))/param_float(4 + 5)

d2_dif_et_dxdx_0_0_ip1nyp2m1k_ip1p1nyp2m1k = ((q(i+1+1,ny+2-1,indvars(4))-&
                    0.5_wp*(q(i+1+1,ny+2-1,indvars(2))*q(i+1+1,ny+2-1,indvars(2))+&
                    q(i+1+1,ny+2-1,indvars(3))*q(i+1+1,ny+2-1,indvars(3)))))/param_float(4 + 5)

d2_dif_et_dxdx_0_0_ip1nyp2m1k = -&
          0.5_wp*d2_dif_et_dxdx_0_0_ip1nyp2m1k_ip1m1nyp2m1k+&
          0.5_wp*d2_dif_et_dxdx_0_0_ip1nyp2m1k_ip1p1nyp2m1k

d2_dif_et_dxdx_0_0_ip1nyp2m1k = d2_dif_et_dxdx_0_0_ip1nyp2m1k*param_float(1)

d1_dif_et_dx_0_im1nyp2m1k = -param_float(2 + 5)*qst(i-1,ny+2-1,indvarsst(10))*(d2_dif_et_dxdx_0_0_im1nyp2m1k)-&
                    q(i-1,ny+2-1,indvars(2))*(((1+&
                    param_float(21 + 5))/(((q(i-1,ny+2-1,indvars(4))-&
                    0.5_wp*(q(i-1,ny+2-1,indvars(2))*q(i-1,ny+2-1,indvars(2))+&
                    q(i-1,ny+2-1,indvars(3))*q(i-1,ny+2-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(i-1,ny+2-1,indvars(4))-&
                    0.5_wp*(q(i-1,ny+2-1,indvars(2))*q(i-1,ny+2-1,indvars(2))+&
                    q(i-1,ny+2-1,indvars(3))*q(i-1,ny+2-1,indvars(3)))))/param_float(4 + 5)**1.5)*(1+&
                    ((q(i-1,ny+2-1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(i-1,ny+2-1,indvars(4))-&
                    0.5_wp*(q(i-1,ny+2-1,indvars(2))*q(i-1,ny+2-1,indvars(2))+&
                    q(i-1,ny+2-1,indvars(3))*q(i-1,ny+2-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(i-1,ny+2-1,indvars(4))-&
                    0.5_wp*(q(i-1,ny+2-1,indvars(2))*q(i-1,ny+2-1,indvars(2))+&
                    q(i-1,ny+2-1,indvars(3))*q(i-1,ny+2-1,indvars(3)))))/param_float(4 + 5)**1.5*q(i-1,ny+2-1,indvars(1)))**3/((q(i-1,ny+2-1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(i-1,ny+2-1,indvars(4))-&
                    0.5_wp*(q(i-1,ny+2-1,indvars(2))*q(i-1,ny+2-1,indvars(2))+&
                    q(i-1,ny+2-1,indvars(3))*q(i-1,ny+2-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(i-1,ny+2-1,indvars(4))-&
                    0.5_wp*(q(i-1,ny+2-1,indvars(2))*q(i-1,ny+2-1,indvars(2))+&
                    q(i-1,ny+2-1,indvars(3))*q(i-1,ny+2-1,indvars(3)))))/param_float(4 + 5)**1.5*q(i-1,ny+2-1,indvars(1)))**3+&
                    param_float(13 + 5)**3))*(q(i-1,ny+2-1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(i-1,ny+2-1,indvars(4))-&
                    0.5_wp*(q(i-1,ny+2-1,indvars(2))*q(i-1,ny+2-1,indvars(2))+&
                    q(i-1,ny+2-1,indvars(3))*q(i-1,ny+2-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(i-1,ny+2-1,indvars(4))-&
                    0.5_wp*(q(i-1,ny+2-1,indvars(2))*q(i-1,ny+2-1,indvars(2))+&
                    q(i-1,ny+2-1,indvars(3))*q(i-1,ny+2-1,indvars(3)))))/param_float(4 + 5)**1.5*q(i-1,ny+2-1,indvars(1))))*param_float(1 + 5)*(2.0_wp*qst(i-1,ny+2-1,indvarsst(10))*(d2_dif_rhou_dxdx_0_0_im1nyp2m1k)-&
                    2.0_wp/3.0_wp*(qst(i-1,ny+2-1,indvarsst(10))*(d2_dif_rhou_dxdx_0_0_im1nyp2m1k)+&
                    qst(i-1,ny+2-1,indvarsst(11))*(d2_dif_rhou_dxdy_0_0_im1nyp2m1k))))-&
                    q(i-1,ny+2-1,indvars(3))*(((1+&
                    param_float(21 + 5))/(((q(i-1,ny+2-1,indvars(4))-&
                    0.5_wp*(q(i-1,ny+2-1,indvars(2))*q(i-1,ny+2-1,indvars(2))+&
                    q(i-1,ny+2-1,indvars(3))*q(i-1,ny+2-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(i-1,ny+2-1,indvars(4))-&
                    0.5_wp*(q(i-1,ny+2-1,indvars(2))*q(i-1,ny+2-1,indvars(2))+&
                    q(i-1,ny+2-1,indvars(3))*q(i-1,ny+2-1,indvars(3)))))/param_float(4 + 5)**1.5)*(1+&
                    ((q(i-1,ny+2-1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(i-1,ny+2-1,indvars(4))-&
                    0.5_wp*(q(i-1,ny+2-1,indvars(2))*q(i-1,ny+2-1,indvars(2))+&
                    q(i-1,ny+2-1,indvars(3))*q(i-1,ny+2-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(i-1,ny+2-1,indvars(4))-&
                    0.5_wp*(q(i-1,ny+2-1,indvars(2))*q(i-1,ny+2-1,indvars(2))+&
                    q(i-1,ny+2-1,indvars(3))*q(i-1,ny+2-1,indvars(3)))))/param_float(4 + 5)**1.5*q(i-1,ny+2-1,indvars(1)))**3/((q(i-1,ny+2-1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(i-1,ny+2-1,indvars(4))-&
                    0.5_wp*(q(i-1,ny+2-1,indvars(2))*q(i-1,ny+2-1,indvars(2))+&
                    q(i-1,ny+2-1,indvars(3))*q(i-1,ny+2-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(i-1,ny+2-1,indvars(4))-&
                    0.5_wp*(q(i-1,ny+2-1,indvars(2))*q(i-1,ny+2-1,indvars(2))+&
                    q(i-1,ny+2-1,indvars(3))*q(i-1,ny+2-1,indvars(3)))))/param_float(4 + 5)**1.5*q(i-1,ny+2-1,indvars(1)))**3+&
                    param_float(13 + 5)**3))*(q(i-1,ny+2-1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(i-1,ny+2-1,indvars(4))-&
                    0.5_wp*(q(i-1,ny+2-1,indvars(2))*q(i-1,ny+2-1,indvars(2))+&
                    q(i-1,ny+2-1,indvars(3))*q(i-1,ny+2-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(i-1,ny+2-1,indvars(4))-&
                    0.5_wp*(q(i-1,ny+2-1,indvars(2))*q(i-1,ny+2-1,indvars(2))+&
                    q(i-1,ny+2-1,indvars(3))*q(i-1,ny+2-1,indvars(3)))))/param_float(4 + 5)**1.5*q(i-1,ny+2-1,indvars(1))))*param_float(1 + 5)*(qst(i-1,ny+2-1,indvarsst(11))*(d2_dif_rhov_dxdy_0_0_im1nyp2m1k)+&
                    qst(i-1,ny+2-1,indvarsst(10))*(d2_dif_rhov_dxdx_0_0_im1nyp2m1k)))

d1_dif_et_dx_0_ip1nyp2m1k = -param_float(2 + 5)*qst(i+1,ny+2-1,indvarsst(10))*(d2_dif_et_dxdx_0_0_ip1nyp2m1k)-&
                    q(i+1,ny+2-1,indvars(2))*(((1+&
                    param_float(21 + 5))/(((q(i+1,ny+2-1,indvars(4))-&
                    0.5_wp*(q(i+1,ny+2-1,indvars(2))*q(i+1,ny+2-1,indvars(2))+&
                    q(i+1,ny+2-1,indvars(3))*q(i+1,ny+2-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(i+1,ny+2-1,indvars(4))-&
                    0.5_wp*(q(i+1,ny+2-1,indvars(2))*q(i+1,ny+2-1,indvars(2))+&
                    q(i+1,ny+2-1,indvars(3))*q(i+1,ny+2-1,indvars(3)))))/param_float(4 + 5)**1.5)*(1+&
                    ((q(i+1,ny+2-1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(i+1,ny+2-1,indvars(4))-&
                    0.5_wp*(q(i+1,ny+2-1,indvars(2))*q(i+1,ny+2-1,indvars(2))+&
                    q(i+1,ny+2-1,indvars(3))*q(i+1,ny+2-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(i+1,ny+2-1,indvars(4))-&
                    0.5_wp*(q(i+1,ny+2-1,indvars(2))*q(i+1,ny+2-1,indvars(2))+&
                    q(i+1,ny+2-1,indvars(3))*q(i+1,ny+2-1,indvars(3)))))/param_float(4 + 5)**1.5*q(i+1,ny+2-1,indvars(1)))**3/((q(i+1,ny+2-1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(i+1,ny+2-1,indvars(4))-&
                    0.5_wp*(q(i+1,ny+2-1,indvars(2))*q(i+1,ny+2-1,indvars(2))+&
                    q(i+1,ny+2-1,indvars(3))*q(i+1,ny+2-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(i+1,ny+2-1,indvars(4))-&
                    0.5_wp*(q(i+1,ny+2-1,indvars(2))*q(i+1,ny+2-1,indvars(2))+&
                    q(i+1,ny+2-1,indvars(3))*q(i+1,ny+2-1,indvars(3)))))/param_float(4 + 5)**1.5*q(i+1,ny+2-1,indvars(1)))**3+&
                    param_float(13 + 5)**3))*(q(i+1,ny+2-1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(i+1,ny+2-1,indvars(4))-&
                    0.5_wp*(q(i+1,ny+2-1,indvars(2))*q(i+1,ny+2-1,indvars(2))+&
                    q(i+1,ny+2-1,indvars(3))*q(i+1,ny+2-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(i+1,ny+2-1,indvars(4))-&
                    0.5_wp*(q(i+1,ny+2-1,indvars(2))*q(i+1,ny+2-1,indvars(2))+&
                    q(i+1,ny+2-1,indvars(3))*q(i+1,ny+2-1,indvars(3)))))/param_float(4 + 5)**1.5*q(i+1,ny+2-1,indvars(1))))*param_float(1 + 5)*(2.0_wp*qst(i+1,ny+2-1,indvarsst(10))*(d2_dif_rhou_dxdx_0_0_ip1nyp2m1k)-&
                    2.0_wp/3.0_wp*(qst(i+1,ny+2-1,indvarsst(10))*(d2_dif_rhou_dxdx_0_0_ip1nyp2m1k)+&
                    qst(i+1,ny+2-1,indvarsst(11))*(d2_dif_rhou_dxdy_0_0_ip1nyp2m1k))))-&
                    q(i+1,ny+2-1,indvars(3))*(((1+&
                    param_float(21 + 5))/(((q(i+1,ny+2-1,indvars(4))-&
                    0.5_wp*(q(i+1,ny+2-1,indvars(2))*q(i+1,ny+2-1,indvars(2))+&
                    q(i+1,ny+2-1,indvars(3))*q(i+1,ny+2-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(i+1,ny+2-1,indvars(4))-&
                    0.5_wp*(q(i+1,ny+2-1,indvars(2))*q(i+1,ny+2-1,indvars(2))+&
                    q(i+1,ny+2-1,indvars(3))*q(i+1,ny+2-1,indvars(3)))))/param_float(4 + 5)**1.5)*(1+&
                    ((q(i+1,ny+2-1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(i+1,ny+2-1,indvars(4))-&
                    0.5_wp*(q(i+1,ny+2-1,indvars(2))*q(i+1,ny+2-1,indvars(2))+&
                    q(i+1,ny+2-1,indvars(3))*q(i+1,ny+2-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(i+1,ny+2-1,indvars(4))-&
                    0.5_wp*(q(i+1,ny+2-1,indvars(2))*q(i+1,ny+2-1,indvars(2))+&
                    q(i+1,ny+2-1,indvars(3))*q(i+1,ny+2-1,indvars(3)))))/param_float(4 + 5)**1.5*q(i+1,ny+2-1,indvars(1)))**3/((q(i+1,ny+2-1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(i+1,ny+2-1,indvars(4))-&
                    0.5_wp*(q(i+1,ny+2-1,indvars(2))*q(i+1,ny+2-1,indvars(2))+&
                    q(i+1,ny+2-1,indvars(3))*q(i+1,ny+2-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(i+1,ny+2-1,indvars(4))-&
                    0.5_wp*(q(i+1,ny+2-1,indvars(2))*q(i+1,ny+2-1,indvars(2))+&
                    q(i+1,ny+2-1,indvars(3))*q(i+1,ny+2-1,indvars(3)))))/param_float(4 + 5)**1.5*q(i+1,ny+2-1,indvars(1)))**3+&
                    param_float(13 + 5)**3))*(q(i+1,ny+2-1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(i+1,ny+2-1,indvars(4))-&
                    0.5_wp*(q(i+1,ny+2-1,indvars(2))*q(i+1,ny+2-1,indvars(2))+&
                    q(i+1,ny+2-1,indvars(3))*q(i+1,ny+2-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(i+1,ny+2-1,indvars(4))-&
                    0.5_wp*(q(i+1,ny+2-1,indvars(2))*q(i+1,ny+2-1,indvars(2))+&
                    q(i+1,ny+2-1,indvars(3))*q(i+1,ny+2-1,indvars(3)))))/param_float(4 + 5)**1.5*q(i+1,ny+2-1,indvars(1))))*param_float(1 + 5)*(qst(i+1,ny+2-1,indvarsst(11))*(d2_dif_rhov_dxdy_0_0_ip1nyp2m1k)+&
                    qst(i+1,ny+2-1,indvarsst(10))*(d2_dif_rhov_dxdx_0_0_ip1nyp2m1k)))

d1_dif_et_dx_0_inyp2m1k = -&
          0.5_wp*d1_dif_et_dx_0_im1nyp2m1k+&
          0.5_wp*d1_dif_et_dx_0_ip1nyp2m1k

d1_dif_et_dx_0_inyp2m1k = d1_dif_et_dx_0_inyp2m1k*param_float(1)

d2_dif_et_dydy_0_0_inyp2m1m1k_inyp2m1m1m1k = ((q(i,ny+2-1-1-1,indvars(4))-&
                    0.5_wp*(q(i,ny+2-1-1-1,indvars(2))*q(i,ny+2-1-1-1,indvars(2))+&
                    q(i,ny+2-1-1-1,indvars(3))*q(i,ny+2-1-1-1,indvars(3)))))/param_float(4 + 5)

d2_dif_et_dydy_0_0_inyp2m1m1k_inyp2m1m1p1k = ((q(i,ny+2-1-1+1,indvars(4))-&
                    0.5_wp*(q(i,ny+2-1-1+1,indvars(2))*q(i,ny+2-1-1+1,indvars(2))+&
                    q(i,ny+2-1-1+1,indvars(3))*q(i,ny+2-1-1+1,indvars(3)))))/param_float(4 + 5)

d2_dif_et_dydy_0_0_inyp2m1m1k = -&
          0.5_wp*d2_dif_et_dydy_0_0_inyp2m1m1k_inyp2m1m1m1k+&
          0.5_wp*d2_dif_et_dydy_0_0_inyp2m1m1k_inyp2m1m1p1k

d2_dif_et_dydy_0_0_inyp2m1m1k = d2_dif_et_dydy_0_0_inyp2m1m1k*param_float(2)

d2_dif_et_dydy_0_0_inyp2m1p1k_inyp2m1p1p0k = ((q(i,ny+2-1+1+0,indvars(4))-&
                    0.5_wp*(q(i,ny+2-1+1+0,indvars(2))*q(i,ny+2-1+1+0,indvars(2))+&
                    q(i,ny+2-1+1+0,indvars(3))*q(i,ny+2-1+1+0,indvars(3)))))/param_float(4 + 5)

d2_dif_et_dydy_0_0_inyp2m1p1k_inyp2m1p1m1k = ((q(i,ny+2-1+1-1,indvars(4))-&
                    0.5_wp*(q(i,ny+2-1+1-1,indvars(2))*q(i,ny+2-1+1-1,indvars(2))+&
                    q(i,ny+2-1+1-1,indvars(3))*q(i,ny+2-1+1-1,indvars(3)))))/param_float(4 + 5)

d2_dif_et_dydy_0_0_inyp2m1p1k_inyp2m1p1m2k = ((q(i,ny+2-1+1-2,indvars(4))-&
                    0.5_wp*(q(i,ny+2-1+1-2,indvars(2))*q(i,ny+2-1+1-2,indvars(2))+&
                    q(i,ny+2-1+1-2,indvars(3))*q(i,ny+2-1+1-2,indvars(3)))))/param_float(4 + 5)

d2_dif_et_dydy_0_0_inyp2m1p1k = 1.5_wp*d2_dif_et_dydy_0_0_inyp2m1p1k_inyp2m1p1p0k-&
          2.0_wp*d2_dif_et_dydy_0_0_inyp2m1p1k_inyp2m1p1m1k+&
          0.5_wp*d2_dif_et_dydy_0_0_inyp2m1p1k_inyp2m1p1m2k

d2_dif_et_dydy_0_0_inyp2m1p1k = d2_dif_et_dydy_0_0_inyp2m1p1k*param_float(2)

d1_dif_et_dy_0_inyp2m1m1k = -param_float(2 + 5)*qst(i,ny+2-1-1,indvarsst(11))*(d2_dif_et_dydy_0_0_inyp2m1m1k)-&
                    q(i,ny+2-1-1,indvars(2))*(((1+&
                    param_float(21 + 5))/(((q(i,ny+2-1-1,indvars(4))-&
                    0.5_wp*(q(i,ny+2-1-1,indvars(2))*q(i,ny+2-1-1,indvars(2))+&
                    q(i,ny+2-1-1,indvars(3))*q(i,ny+2-1-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(i,ny+2-1-1,indvars(4))-&
                    0.5_wp*(q(i,ny+2-1-1,indvars(2))*q(i,ny+2-1-1,indvars(2))+&
                    q(i,ny+2-1-1,indvars(3))*q(i,ny+2-1-1,indvars(3)))))/param_float(4 + 5)**1.5)*(1+&
                    ((q(i,ny+2-1-1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(i,ny+2-1-1,indvars(4))-&
                    0.5_wp*(q(i,ny+2-1-1,indvars(2))*q(i,ny+2-1-1,indvars(2))+&
                    q(i,ny+2-1-1,indvars(3))*q(i,ny+2-1-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(i,ny+2-1-1,indvars(4))-&
                    0.5_wp*(q(i,ny+2-1-1,indvars(2))*q(i,ny+2-1-1,indvars(2))+&
                    q(i,ny+2-1-1,indvars(3))*q(i,ny+2-1-1,indvars(3)))))/param_float(4 + 5)**1.5*q(i,ny+2-1-1,indvars(1)))**3/((q(i,ny+2-1-1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(i,ny+2-1-1,indvars(4))-&
                    0.5_wp*(q(i,ny+2-1-1,indvars(2))*q(i,ny+2-1-1,indvars(2))+&
                    q(i,ny+2-1-1,indvars(3))*q(i,ny+2-1-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(i,ny+2-1-1,indvars(4))-&
                    0.5_wp*(q(i,ny+2-1-1,indvars(2))*q(i,ny+2-1-1,indvars(2))+&
                    q(i,ny+2-1-1,indvars(3))*q(i,ny+2-1-1,indvars(3)))))/param_float(4 + 5)**1.5*q(i,ny+2-1-1,indvars(1)))**3+&
                    param_float(13 + 5)**3))*(q(i,ny+2-1-1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(i,ny+2-1-1,indvars(4))-&
                    0.5_wp*(q(i,ny+2-1-1,indvars(2))*q(i,ny+2-1-1,indvars(2))+&
                    q(i,ny+2-1-1,indvars(3))*q(i,ny+2-1-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(i,ny+2-1-1,indvars(4))-&
                    0.5_wp*(q(i,ny+2-1-1,indvars(2))*q(i,ny+2-1-1,indvars(2))+&
                    q(i,ny+2-1-1,indvars(3))*q(i,ny+2-1-1,indvars(3)))))/param_float(4 + 5)**1.5*q(i,ny+2-1-1,indvars(1))))*param_float(1 + 5)*(qst(i,ny+2-1-1,indvarsst(11))*(d2_dif_rhou_dydy_0_0_inyp2m1m1k)+&
                    qst(i,ny+2-1-1,indvarsst(10))*(d2_dif_rhou_dydx_0_0_inyp2m1m1k)))-&
                    q(i,ny+2-1-1,indvars(3))*(((1+&
                    param_float(21 + 5))/(((q(i,ny+2-1-1,indvars(4))-&
                    0.5_wp*(q(i,ny+2-1-1,indvars(2))*q(i,ny+2-1-1,indvars(2))+&
                    q(i,ny+2-1-1,indvars(3))*q(i,ny+2-1-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(i,ny+2-1-1,indvars(4))-&
                    0.5_wp*(q(i,ny+2-1-1,indvars(2))*q(i,ny+2-1-1,indvars(2))+&
                    q(i,ny+2-1-1,indvars(3))*q(i,ny+2-1-1,indvars(3)))))/param_float(4 + 5)**1.5)*(1+&
                    ((q(i,ny+2-1-1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(i,ny+2-1-1,indvars(4))-&
                    0.5_wp*(q(i,ny+2-1-1,indvars(2))*q(i,ny+2-1-1,indvars(2))+&
                    q(i,ny+2-1-1,indvars(3))*q(i,ny+2-1-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(i,ny+2-1-1,indvars(4))-&
                    0.5_wp*(q(i,ny+2-1-1,indvars(2))*q(i,ny+2-1-1,indvars(2))+&
                    q(i,ny+2-1-1,indvars(3))*q(i,ny+2-1-1,indvars(3)))))/param_float(4 + 5)**1.5*q(i,ny+2-1-1,indvars(1)))**3/((q(i,ny+2-1-1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(i,ny+2-1-1,indvars(4))-&
                    0.5_wp*(q(i,ny+2-1-1,indvars(2))*q(i,ny+2-1-1,indvars(2))+&
                    q(i,ny+2-1-1,indvars(3))*q(i,ny+2-1-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(i,ny+2-1-1,indvars(4))-&
                    0.5_wp*(q(i,ny+2-1-1,indvars(2))*q(i,ny+2-1-1,indvars(2))+&
                    q(i,ny+2-1-1,indvars(3))*q(i,ny+2-1-1,indvars(3)))))/param_float(4 + 5)**1.5*q(i,ny+2-1-1,indvars(1)))**3+&
                    param_float(13 + 5)**3))*(q(i,ny+2-1-1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(i,ny+2-1-1,indvars(4))-&
                    0.5_wp*(q(i,ny+2-1-1,indvars(2))*q(i,ny+2-1-1,indvars(2))+&
                    q(i,ny+2-1-1,indvars(3))*q(i,ny+2-1-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(i,ny+2-1-1,indvars(4))-&
                    0.5_wp*(q(i,ny+2-1-1,indvars(2))*q(i,ny+2-1-1,indvars(2))+&
                    q(i,ny+2-1-1,indvars(3))*q(i,ny+2-1-1,indvars(3)))))/param_float(4 + 5)**1.5*q(i,ny+2-1-1,indvars(1))))*param_float(1 + 5)*(2.0_wp*qst(i,ny+2-1-1,indvarsst(11))*(d2_dif_rhov_dydy_0_0_inyp2m1m1k)-&
                    2.0_wp/3.0_wp*(qst(i,ny+2-1-1,indvarsst(10))*(d2_dif_rhov_dydx_0_0_inyp2m1m1k)+&
                    qst(i,ny+2-1-1,indvarsst(11))*(d2_dif_rhov_dydy_0_0_inyp2m1m1k))))

d1_dif_et_dy_0_inyp2m1p1k = -param_float(2 + 5)*qst(i,ny+2-1+1,indvarsst(11))*(d2_dif_et_dydy_0_0_inyp2m1p1k)-&
                    q(i,ny+2-1+1,indvars(2))*(((1+&
                    param_float(21 + 5))/(((q(i,ny+2-1+1,indvars(4))-&
                    0.5_wp*(q(i,ny+2-1+1,indvars(2))*q(i,ny+2-1+1,indvars(2))+&
                    q(i,ny+2-1+1,indvars(3))*q(i,ny+2-1+1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(i,ny+2-1+1,indvars(4))-&
                    0.5_wp*(q(i,ny+2-1+1,indvars(2))*q(i,ny+2-1+1,indvars(2))+&
                    q(i,ny+2-1+1,indvars(3))*q(i,ny+2-1+1,indvars(3)))))/param_float(4 + 5)**1.5)*(1+&
                    ((q(i,ny+2-1+1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(i,ny+2-1+1,indvars(4))-&
                    0.5_wp*(q(i,ny+2-1+1,indvars(2))*q(i,ny+2-1+1,indvars(2))+&
                    q(i,ny+2-1+1,indvars(3))*q(i,ny+2-1+1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(i,ny+2-1+1,indvars(4))-&
                    0.5_wp*(q(i,ny+2-1+1,indvars(2))*q(i,ny+2-1+1,indvars(2))+&
                    q(i,ny+2-1+1,indvars(3))*q(i,ny+2-1+1,indvars(3)))))/param_float(4 + 5)**1.5*q(i,ny+2-1+1,indvars(1)))**3/((q(i,ny+2-1+1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(i,ny+2-1+1,indvars(4))-&
                    0.5_wp*(q(i,ny+2-1+1,indvars(2))*q(i,ny+2-1+1,indvars(2))+&
                    q(i,ny+2-1+1,indvars(3))*q(i,ny+2-1+1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(i,ny+2-1+1,indvars(4))-&
                    0.5_wp*(q(i,ny+2-1+1,indvars(2))*q(i,ny+2-1+1,indvars(2))+&
                    q(i,ny+2-1+1,indvars(3))*q(i,ny+2-1+1,indvars(3)))))/param_float(4 + 5)**1.5*q(i,ny+2-1+1,indvars(1)))**3+&
                    param_float(13 + 5)**3))*(q(i,ny+2-1+1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(i,ny+2-1+1,indvars(4))-&
                    0.5_wp*(q(i,ny+2-1+1,indvars(2))*q(i,ny+2-1+1,indvars(2))+&
                    q(i,ny+2-1+1,indvars(3))*q(i,ny+2-1+1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(i,ny+2-1+1,indvars(4))-&
                    0.5_wp*(q(i,ny+2-1+1,indvars(2))*q(i,ny+2-1+1,indvars(2))+&
                    q(i,ny+2-1+1,indvars(3))*q(i,ny+2-1+1,indvars(3)))))/param_float(4 + 5)**1.5*q(i,ny+2-1+1,indvars(1))))*param_float(1 + 5)*(qst(i,ny+2-1+1,indvarsst(11))*(d2_dif_rhou_dydy_0_0_inyp2m1p1k)+&
                    qst(i,ny+2-1+1,indvarsst(10))*(d2_dif_rhou_dydx_0_0_inyp2m1p1k)))-&
                    q(i,ny+2-1+1,indvars(3))*(((1+&
                    param_float(21 + 5))/(((q(i,ny+2-1+1,indvars(4))-&
                    0.5_wp*(q(i,ny+2-1+1,indvars(2))*q(i,ny+2-1+1,indvars(2))+&
                    q(i,ny+2-1+1,indvars(3))*q(i,ny+2-1+1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(i,ny+2-1+1,indvars(4))-&
                    0.5_wp*(q(i,ny+2-1+1,indvars(2))*q(i,ny+2-1+1,indvars(2))+&
                    q(i,ny+2-1+1,indvars(3))*q(i,ny+2-1+1,indvars(3)))))/param_float(4 + 5)**1.5)*(1+&
                    ((q(i,ny+2-1+1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(i,ny+2-1+1,indvars(4))-&
                    0.5_wp*(q(i,ny+2-1+1,indvars(2))*q(i,ny+2-1+1,indvars(2))+&
                    q(i,ny+2-1+1,indvars(3))*q(i,ny+2-1+1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(i,ny+2-1+1,indvars(4))-&
                    0.5_wp*(q(i,ny+2-1+1,indvars(2))*q(i,ny+2-1+1,indvars(2))+&
                    q(i,ny+2-1+1,indvars(3))*q(i,ny+2-1+1,indvars(3)))))/param_float(4 + 5)**1.5*q(i,ny+2-1+1,indvars(1)))**3/((q(i,ny+2-1+1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(i,ny+2-1+1,indvars(4))-&
                    0.5_wp*(q(i,ny+2-1+1,indvars(2))*q(i,ny+2-1+1,indvars(2))+&
                    q(i,ny+2-1+1,indvars(3))*q(i,ny+2-1+1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(i,ny+2-1+1,indvars(4))-&
                    0.5_wp*(q(i,ny+2-1+1,indvars(2))*q(i,ny+2-1+1,indvars(2))+&
                    q(i,ny+2-1+1,indvars(3))*q(i,ny+2-1+1,indvars(3)))))/param_float(4 + 5)**1.5*q(i,ny+2-1+1,indvars(1)))**3+&
                    param_float(13 + 5)**3))*(q(i,ny+2-1+1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(i,ny+2-1+1,indvars(4))-&
                    0.5_wp*(q(i,ny+2-1+1,indvars(2))*q(i,ny+2-1+1,indvars(2))+&
                    q(i,ny+2-1+1,indvars(3))*q(i,ny+2-1+1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(i,ny+2-1+1,indvars(4))-&
                    0.5_wp*(q(i,ny+2-1+1,indvars(2))*q(i,ny+2-1+1,indvars(2))+&
                    q(i,ny+2-1+1,indvars(3))*q(i,ny+2-1+1,indvars(3)))))/param_float(4 + 5)**1.5*q(i,ny+2-1+1,indvars(1))))*param_float(1 + 5)*(2.0_wp*qst(i,ny+2-1+1,indvarsst(11))*(d2_dif_rhov_dydy_0_0_inyp2m1p1k)-&
                    2.0_wp/3.0_wp*(qst(i,ny+2-1+1,indvarsst(10))*(d2_dif_rhov_dydx_0_0_inyp2m1p1k)+&
                    qst(i,ny+2-1+1,indvarsst(11))*(d2_dif_rhov_dydy_0_0_inyp2m1p1k))))

d1_dif_et_dy_0_inyp2m1k = -&
          0.5_wp*d1_dif_et_dy_0_inyp2m1m1k+&
          0.5_wp*d1_dif_et_dy_0_inyp2m1p1k

d1_dif_et_dy_0_inyp2m1k = d1_dif_et_dy_0_inyp2m1k*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer None 1 None d(rho et)/dt *******
!                                                           
!***********************************************************


rhs(i,ny+2-1,indvars(4)) = rhs(i,ny+2-1,indvars(4))  -  ( qst(i,ny+2-1,indvarsst(10))*(d1_dif_et_dx_0_inyp2m1k)+&
                    qst(i,ny+2-1,indvarsst(11))*(d1_dif_et_dy_0_inyp2m1k) ) 



!***********************************************************
!                                                           
! building source terms in RHS for layer None 1 None d(rho nut)/dt 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! -ReI*Cb2*sigmaI*((deltaxI)**2*([rho*nut]_1x)*([nut]_1x)+(deltayI)**2*([rho*nut]_1y)*([nut]_1y))-Cb1*(1-ft2)*SS*rho*nut+ReI*(Cw1*fw-Cb1/k**2*ft2)*rho*nut**2/eta**2
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_dif_nut_dx_0_im1nyp2m1k = q(i-1,ny+2-1,indvars(1))*q(i-1,ny+2-1,indvars(5))

d1_dif_nut_dx_0_ip1nyp2m1k = q(i+1,ny+2-1,indvars(1))*q(i+1,ny+2-1,indvars(5))

d1_dif_nut_dx_0_inyp2m1k = -&
          0.5_wp*d1_dif_nut_dx_0_im1nyp2m1k+&
          0.5_wp*d1_dif_nut_dx_0_ip1nyp2m1k

d1_dif_nut_dx_0_inyp2m1k = d1_dif_nut_dx_0_inyp2m1k*param_float(1)

d1_dif_nut_dx_1_im1nyp2m1k = q(i-1,ny+2-1,indvars(5))

d1_dif_nut_dx_1_ip1nyp2m1k = q(i+1,ny+2-1,indvars(5))

d1_dif_nut_dx_1_inyp2m1k = -&
          0.5_wp*d1_dif_nut_dx_1_im1nyp2m1k+&
          0.5_wp*d1_dif_nut_dx_1_ip1nyp2m1k

d1_dif_nut_dx_1_inyp2m1k = d1_dif_nut_dx_1_inyp2m1k*param_float(1)

d1_dif_nut_dy_0_inyp2m1m1k = q(i,ny+2-1-1,indvars(1))*q(i,ny+2-1-1,indvars(5))

d1_dif_nut_dy_0_inyp2m1p1k = q(i,ny+2-1+1,indvars(1))*q(i,ny+2-1+1,indvars(5))

d1_dif_nut_dy_0_inyp2m1k = -&
          0.5_wp*d1_dif_nut_dy_0_inyp2m1m1k+&
          0.5_wp*d1_dif_nut_dy_0_inyp2m1p1k

d1_dif_nut_dy_0_inyp2m1k = d1_dif_nut_dy_0_inyp2m1k*param_float(2)

d1_dif_nut_dy_1_inyp2m1m1k = q(i,ny+2-1-1,indvars(5))

d1_dif_nut_dy_1_inyp2m1p1k = q(i,ny+2-1+1,indvars(5))

d1_dif_nut_dy_1_inyp2m1k = -&
          0.5_wp*d1_dif_nut_dy_1_inyp2m1m1k+&
          0.5_wp*d1_dif_nut_dy_1_inyp2m1p1k

d1_dif_nut_dy_1_inyp2m1k = d1_dif_nut_dy_1_inyp2m1k*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer None 1 None d(rho nut)/dt ******
!                                                           
!***********************************************************


rhs(i,ny+2-1,indvars(5)) = rhs(i,ny+2-1,indvars(5))  -  ( -param_float(1 + 5)*param_float(7 + 5)*param_float(18 + 5)*((qst(i,ny+2-1,indvarsst(10)))**2*(d1_dif_nut_dx_0_inyp2m1k)*(d1_dif_nut_dx_1_inyp2m1k)+&
                    (qst(i,ny+2-1,indvarsst(11)))**2*(d1_dif_nut_dy_0_inyp2m1k)*(d1_dif_nut_dy_1_inyp2m1k))-&
                    param_float(6 + 5)*(1-&
                    qst(i,ny+2-1,indvarsst(16)))*qst(i,ny+2-1,indvarsst(20))*q(i,ny+2-1,indvars(1))*q(i,ny+2-1,indvars(5))+&
                    param_float(1 + 5)*(param_float(10 + 5)*qst(i,ny+2-1,indvarsst(13))-&
                    param_float(6 + 5)/param_float(9 + 5)**2*qst(i,ny+2-1,indvarsst(16)))*q(i,ny+2-1,indvars(1))*q(i,ny+2-1,indvars(5))**2/qst(i,ny+2-1,indvarsst(2))**2 ) 

   enddo
