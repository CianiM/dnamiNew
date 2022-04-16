

!***********************************************************
!                                                           
! Start building layers for BC : imax jmax None ************
!                                                           
!***********************************************************




!***********************************************************
!                                                           
! BC layer: 0 0 None ***************************************
!                                                           
!***********************************************************




!***********************************************************
!                                                           
! building source terms in RHS for layer 0 0 None stemp ****
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (((0.5_wp*(deltayI*([u]_1y)-deltaxI*([v]_1x)))**2+(0.5_wp*(deltaxI*([v]_1x)-deltayI*([u]_1y)))**2)*2)**0.5
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_stemp_dx_0_nxp2p0p0nyp2p0k = q(nx+2+0+0,ny+2+0,indvars(3))

d1_stemp_dx_0_nxp2p0m1nyp2p0k = q(nx+2+0-1,ny+2+0,indvars(3))

d1_stemp_dx_0_nxp2p0m2nyp2p0k = q(nx+2+0-2,ny+2+0,indvars(3))

d1_stemp_dx_0_nxp2p0nyp2p0k = 1.5_wp*d1_stemp_dx_0_nxp2p0p0nyp2p0k-&
          2.0_wp*d1_stemp_dx_0_nxp2p0m1nyp2p0k+&
          0.5_wp*d1_stemp_dx_0_nxp2p0m2nyp2p0k

d1_stemp_dx_0_nxp2p0nyp2p0k = d1_stemp_dx_0_nxp2p0nyp2p0k*param_float(1)

d1_stemp_dy_0_nxp2p0nyp2p0p0k = q(nx+2+0,ny+2+0+0,indvars(2))

d1_stemp_dy_0_nxp2p0nyp2p0m1k = q(nx+2+0,ny+2+0-1,indvars(2))

d1_stemp_dy_0_nxp2p0nyp2p0m2k = q(nx+2+0,ny+2+0-2,indvars(2))

d1_stemp_dy_0_nxp2p0nyp2p0k = 1.5_wp*d1_stemp_dy_0_nxp2p0nyp2p0p0k-&
          2.0_wp*d1_stemp_dy_0_nxp2p0nyp2p0m1k+&
          0.5_wp*d1_stemp_dy_0_nxp2p0nyp2p0m2k

d1_stemp_dy_0_nxp2p0nyp2p0k = d1_stemp_dy_0_nxp2p0nyp2p0k*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer 0 0 None stemp *****************
!                                                           
!***********************************************************


qst(nx+2+0,ny+2+0,indvarsst(4)) =  (((0.5_wp*(qst(nx+2+0,ny+2+0,indvarsst(11))*(d1_stemp_dy_0_nxp2p0nyp2p0k)-&
                    qst(nx+2+0,ny+2+0,indvarsst(10))*(d1_stemp_dx_0_nxp2p0nyp2p0k)))**2+&
                    (0.5_wp*(qst(nx+2+0,ny+2+0,indvarsst(10))*(d1_stemp_dx_0_nxp2p0nyp2p0k)-&
                    qst(nx+2+0,ny+2+0,indvarsst(11))*(d1_stemp_dy_0_nxp2p0nyp2p0k)))**2)*2)**0.5



!***********************************************************
!                                                           
! building source terms in RHS for layer 0 0 None viscosità 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (1+sut)/(T+sut)*T**1.5
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 0 0 None viscosità *************
!                                                           
!***********************************************************


qst(nx+2+0,ny+2+0,indvarsst(12)) =  (1+&
                    param_float(21 + 5))/(((q(nx+2+0,ny+2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2+0,indvars(2))*q(nx+2+0,ny+2+0,indvars(2))+&
                    q(nx+2+0,ny+2+0,indvars(3))*q(nx+2+0,ny+2+0,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2+0,ny+2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2+0,indvars(2))*q(nx+2+0,ny+2+0,indvars(2))+&
                    q(nx+2+0,ny+2+0,indvars(3))*q(nx+2+0,ny+2+0,indvars(3)))))/param_float(4 + 5)**1.5



!***********************************************************
!                                                           
! building source terms in RHS for layer 0 0 None fw *******
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! gg*(1+Cw3**6)/(gg**6+Cw3**6)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 0 0 None fw ********************
!                                                           
!***********************************************************


qst(nx+2+0,ny+2+0,indvarsst(13)) =  qst(nx+2+0,ny+2+0,indvarsst(14))*(1+&
                    param_float(12 + 5)**6)/(qst(nx+2+0,ny+2+0,indvarsst(14))**6+&
                    param_float(12 + 5)**6)



!***********************************************************
!                                                           
! building source terms in RHS for layer 0 0 None gg *******
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (nut/(SS*k**2*eta**2))
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 0 0 None gg ********************
!                                                           
!***********************************************************


qst(nx+2+0,ny+2+0,indvarsst(14)) =  (q(nx+2+0,ny+2+0,indvars(5))/(qst(nx+2+0,ny+2+0,indvarsst(20))*param_float(9 + 5)**2*qst(nx+2+0,ny+2+0,indvarsst(2))**2))



!***********************************************************
!                                                           
! building source terms in RHS for layer 0 0 None fv2 ******
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (1-chi/(1+chi*fv1))
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 0 0 None fv2 *******************
!                                                           
!***********************************************************


qst(nx+2+0,ny+2+0,indvarsst(15)) =  (1-&
                    (q(nx+2+0,ny+2+0,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+2+0,ny+2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2+0,indvars(2))*q(nx+2+0,ny+2+0,indvars(2))+&
                    q(nx+2+0,ny+2+0,indvars(3))*q(nx+2+0,ny+2+0,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2+0,ny+2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2+0,indvars(2))*q(nx+2+0,ny+2+0,indvars(2))+&
                    q(nx+2+0,ny+2+0,indvars(3))*q(nx+2+0,ny+2+0,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+2+0,ny+2+0,indvars(1)))/(1+&
                    (q(nx+2+0,ny+2+0,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+2+0,ny+2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2+0,indvars(2))*q(nx+2+0,ny+2+0,indvars(2))+&
                    q(nx+2+0,ny+2+0,indvars(3))*q(nx+2+0,ny+2+0,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2+0,ny+2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2+0,indvars(2))*q(nx+2+0,ny+2+0,indvars(2))+&
                    q(nx+2+0,ny+2+0,indvars(3))*q(nx+2+0,ny+2+0,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+2+0,ny+2+0,indvars(1)))*((q(nx+2+0,ny+2+0,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+2+0,ny+2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2+0,indvars(2))*q(nx+2+0,ny+2+0,indvars(2))+&
                    q(nx+2+0,ny+2+0,indvars(3))*q(nx+2+0,ny+2+0,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2+0,ny+2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2+0,indvars(2))*q(nx+2+0,ny+2+0,indvars(2))+&
                    q(nx+2+0,ny+2+0,indvars(3))*q(nx+2+0,ny+2+0,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+2+0,ny+2+0,indvars(1)))**3/((q(nx+2+0,ny+2+0,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+2+0,ny+2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2+0,indvars(2))*q(nx+2+0,ny+2+0,indvars(2))+&
                    q(nx+2+0,ny+2+0,indvars(3))*q(nx+2+0,ny+2+0,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2+0,ny+2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2+0,indvars(2))*q(nx+2+0,ny+2+0,indvars(2))+&
                    q(nx+2+0,ny+2+0,indvars(3))*q(nx+2+0,ny+2+0,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+2+0,ny+2+0,indvars(1)))**3+&
                    param_float(13 + 5)**3))))



!***********************************************************
!                                                           
! building source terms in RHS for layer 0 0 None ft2 ******
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (Ct3*exp(-Ct4*chi**2))
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 0 0 None ft2 *******************
!                                                           
!***********************************************************


qst(nx+2+0,ny+2+0,indvarsst(16)) =  (param_float(16 + 5)*exp(-&
                    param_float(17 + 5)*(q(nx+2+0,ny+2+0,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+2+0,ny+2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2+0,indvars(2))*q(nx+2+0,ny+2+0,indvars(2))+&
                    q(nx+2+0,ny+2+0,indvars(3))*q(nx+2+0,ny+2+0,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2+0,ny+2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2+0,indvars(2))*q(nx+2+0,ny+2+0,indvars(2))+&
                    q(nx+2+0,ny+2+0,indvars(3))*q(nx+2+0,ny+2+0,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+2+0,ny+2+0,indvars(1)))**2))



!***********************************************************
!                                                           
! building source terms in RHS for layer 0 0 None fw2 ******
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (gg*(1+Cw3**6)/(gg**6+Cw3**6))
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 0 0 None fw2 *******************
!                                                           
!***********************************************************


qst(nx+2+0,ny+2+0,indvarsst(17)) =  (qst(nx+2+0,ny+2+0,indvarsst(14))*(1+&
                    param_float(12 + 5)**6)/(qst(nx+2+0,ny+2+0,indvarsst(14))**6+&
                    param_float(12 + 5)**6))



!***********************************************************
!                                                           
! building source terms in RHS for layer 0 0 None rr *******
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (nut/(SS*k**2*eta**2))
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 0 0 None rr ********************
!                                                           
!***********************************************************


qst(nx+2+0,ny+2+0,indvarsst(18)) =  (q(nx+2+0,ny+2+0,indvars(5))/(qst(nx+2+0,ny+2+0,indvarsst(20))*param_float(9 + 5)**2*qst(nx+2+0,ny+2+0,indvarsst(2))**2))



!***********************************************************
!                                                           
! building source terms in RHS for layer 0 0 None SS *******
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (stemp+ReI*nut/(k**2*eta**2))
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 0 0 None SS ********************
!                                                           
!***********************************************************


qst(nx+2+0,ny+2+0,indvarsst(20)) =  (qst(nx+2+0,ny+2+0,indvarsst(4))+&
                    param_float(1 + 5)*q(nx+2+0,ny+2+0,indvars(5))/(param_float(9 + 5)**2*qst(nx+2+0,ny+2+0,indvarsst(2))**2))

