

!***********************************************************
!                                                           
! Start building layers for BC : imax j1 None **************
!                                                           
!***********************************************************




!***********************************************************
!                                                           
! BC layer: 1 1 None ***************************************
!                                                           
!***********************************************************




!***********************************************************
!                                                           
! building source terms in RHS for layer 1 1 None stemp ****
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (((0.5_wp*(deltayI*([u]_1y)-deltaxI*([v]_1x)))**2+(0.5_wp*(deltaxI*([v]_1x)-deltayI*([u]_1y)))**2)*2)**0.5
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_stemp_dx_0_nxp2m1m11m2p1k = q(nx+2-1-1,1-2+1,indvars(3))

d1_stemp_dx_0_nxp2m1p11m2p1k = q(nx+2-1+1,1-2+1,indvars(3))

d1_stemp_dx_0_nxp2m11m2p1k = -&
          0.5_wp*d1_stemp_dx_0_nxp2m1m11m2p1k+&
          0.5_wp*d1_stemp_dx_0_nxp2m1p11m2p1k

d1_stemp_dx_0_nxp2m11m2p1k = d1_stemp_dx_0_nxp2m11m2p1k*param_float(1)

d1_stemp_dy_0_nxp2m11m2p1m1k = q(nx+2-1,1-2+1-1,indvars(2))

d1_stemp_dy_0_nxp2m11m2p1p1k = q(nx+2-1,1-2+1+1,indvars(2))

d1_stemp_dy_0_nxp2m11m2p1k = -&
          0.5_wp*d1_stemp_dy_0_nxp2m11m2p1m1k+&
          0.5_wp*d1_stemp_dy_0_nxp2m11m2p1p1k

d1_stemp_dy_0_nxp2m11m2p1k = d1_stemp_dy_0_nxp2m11m2p1k*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer 1 1 None stemp *****************
!                                                           
!***********************************************************


qst(nx+2-1,1-2+1,indvarsst(4)) =  (((0.5_wp*(qst(nx+2-1,1-2+1,indvarsst(11))*(d1_stemp_dy_0_nxp2m11m2p1k)-&
                    qst(nx+2-1,1-2+1,indvarsst(10))*(d1_stemp_dx_0_nxp2m11m2p1k)))**2+&
                    (0.5_wp*(qst(nx+2-1,1-2+1,indvarsst(10))*(d1_stemp_dx_0_nxp2m11m2p1k)-&
                    qst(nx+2-1,1-2+1,indvarsst(11))*(d1_stemp_dy_0_nxp2m11m2p1k)))**2)*2)**0.5



!***********************************************************
!                                                           
! building source terms in RHS for layer 1 1 None viscosità 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (1+sut)/(T+sut)*T**1.5
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 1 1 None viscosità *************
!                                                           
!***********************************************************


qst(nx+2-1,1-2+1,indvarsst(12)) =  (1+&
                    param_float(21 + 5))/(((q(nx+2-1,1-2+1,indvars(4))-&
                    0.5_wp*(q(nx+2-1,1-2+1,indvars(2))*q(nx+2-1,1-2+1,indvars(2))+&
                    q(nx+2-1,1-2+1,indvars(3))*q(nx+2-1,1-2+1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2-1,1-2+1,indvars(4))-&
                    0.5_wp*(q(nx+2-1,1-2+1,indvars(2))*q(nx+2-1,1-2+1,indvars(2))+&
                    q(nx+2-1,1-2+1,indvars(3))*q(nx+2-1,1-2+1,indvars(3)))))/param_float(4 + 5)**1.5



!***********************************************************
!                                                           
! building source terms in RHS for layer 1 1 None fw *******
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! gg*(1+Cw3**6)/(gg**6+Cw3**6)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 1 1 None fw ********************
!                                                           
!***********************************************************


qst(nx+2-1,1-2+1,indvarsst(13)) =  qst(nx+2-1,1-2+1,indvarsst(14))*(1+&
                    param_float(12 + 5)**6)/(qst(nx+2-1,1-2+1,indvarsst(14))**6+&
                    param_float(12 + 5)**6)



!***********************************************************
!                                                           
! building source terms in RHS for layer 1 1 None gg *******
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (nut/(SS*k**2*eta**2))
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 1 1 None gg ********************
!                                                           
!***********************************************************


qst(nx+2-1,1-2+1,indvarsst(14)) =  (q(nx+2-1,1-2+1,indvars(5))/(qst(nx+2-1,1-2+1,indvarsst(20))*param_float(9 + 5)**2*qst(nx+2-1,1-2+1,indvarsst(2))**2))



!***********************************************************
!                                                           
! building source terms in RHS for layer 1 1 None fv2 ******
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (1-chi/(1+chi*fv1))
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 1 1 None fv2 *******************
!                                                           
!***********************************************************


qst(nx+2-1,1-2+1,indvarsst(15)) =  (1-&
                    (q(nx+2-1,1-2+1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+2-1,1-2+1,indvars(4))-&
                    0.5_wp*(q(nx+2-1,1-2+1,indvars(2))*q(nx+2-1,1-2+1,indvars(2))+&
                    q(nx+2-1,1-2+1,indvars(3))*q(nx+2-1,1-2+1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2-1,1-2+1,indvars(4))-&
                    0.5_wp*(q(nx+2-1,1-2+1,indvars(2))*q(nx+2-1,1-2+1,indvars(2))+&
                    q(nx+2-1,1-2+1,indvars(3))*q(nx+2-1,1-2+1,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+2-1,1-2+1,indvars(1)))/(1+&
                    (q(nx+2-1,1-2+1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+2-1,1-2+1,indvars(4))-&
                    0.5_wp*(q(nx+2-1,1-2+1,indvars(2))*q(nx+2-1,1-2+1,indvars(2))+&
                    q(nx+2-1,1-2+1,indvars(3))*q(nx+2-1,1-2+1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2-1,1-2+1,indvars(4))-&
                    0.5_wp*(q(nx+2-1,1-2+1,indvars(2))*q(nx+2-1,1-2+1,indvars(2))+&
                    q(nx+2-1,1-2+1,indvars(3))*q(nx+2-1,1-2+1,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+2-1,1-2+1,indvars(1)))*((q(nx+2-1,1-2+1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+2-1,1-2+1,indvars(4))-&
                    0.5_wp*(q(nx+2-1,1-2+1,indvars(2))*q(nx+2-1,1-2+1,indvars(2))+&
                    q(nx+2-1,1-2+1,indvars(3))*q(nx+2-1,1-2+1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2-1,1-2+1,indvars(4))-&
                    0.5_wp*(q(nx+2-1,1-2+1,indvars(2))*q(nx+2-1,1-2+1,indvars(2))+&
                    q(nx+2-1,1-2+1,indvars(3))*q(nx+2-1,1-2+1,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+2-1,1-2+1,indvars(1)))**3/((q(nx+2-1,1-2+1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+2-1,1-2+1,indvars(4))-&
                    0.5_wp*(q(nx+2-1,1-2+1,indvars(2))*q(nx+2-1,1-2+1,indvars(2))+&
                    q(nx+2-1,1-2+1,indvars(3))*q(nx+2-1,1-2+1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2-1,1-2+1,indvars(4))-&
                    0.5_wp*(q(nx+2-1,1-2+1,indvars(2))*q(nx+2-1,1-2+1,indvars(2))+&
                    q(nx+2-1,1-2+1,indvars(3))*q(nx+2-1,1-2+1,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+2-1,1-2+1,indvars(1)))**3+&
                    param_float(13 + 5)**3))))



!***********************************************************
!                                                           
! building source terms in RHS for layer 1 1 None ft2 ******
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (Ct3*exp(-Ct4*chi**2))
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 1 1 None ft2 *******************
!                                                           
!***********************************************************


qst(nx+2-1,1-2+1,indvarsst(16)) =  (param_float(16 + 5)*exp(-&
                    param_float(17 + 5)*(q(nx+2-1,1-2+1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+2-1,1-2+1,indvars(4))-&
                    0.5_wp*(q(nx+2-1,1-2+1,indvars(2))*q(nx+2-1,1-2+1,indvars(2))+&
                    q(nx+2-1,1-2+1,indvars(3))*q(nx+2-1,1-2+1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2-1,1-2+1,indvars(4))-&
                    0.5_wp*(q(nx+2-1,1-2+1,indvars(2))*q(nx+2-1,1-2+1,indvars(2))+&
                    q(nx+2-1,1-2+1,indvars(3))*q(nx+2-1,1-2+1,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+2-1,1-2+1,indvars(1)))**2))



!***********************************************************
!                                                           
! building source terms in RHS for layer 1 1 None fw2 ******
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (gg*(1+Cw3**6)/(gg**6+Cw3**6))
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 1 1 None fw2 *******************
!                                                           
!***********************************************************


qst(nx+2-1,1-2+1,indvarsst(17)) =  (qst(nx+2-1,1-2+1,indvarsst(14))*(1+&
                    param_float(12 + 5)**6)/(qst(nx+2-1,1-2+1,indvarsst(14))**6+&
                    param_float(12 + 5)**6))



!***********************************************************
!                                                           
! building source terms in RHS for layer 1 1 None rr *******
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (nut/(SS*k**2*eta**2))
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 1 1 None rr ********************
!                                                           
!***********************************************************


qst(nx+2-1,1-2+1,indvarsst(18)) =  (q(nx+2-1,1-2+1,indvars(5))/(qst(nx+2-1,1-2+1,indvarsst(20))*param_float(9 + 5)**2*qst(nx+2-1,1-2+1,indvarsst(2))**2))



!***********************************************************
!                                                           
! building source terms in RHS for layer 1 1 None SS *******
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (stemp+ReI*nut/(k**2*eta**2))
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 1 1 None SS ********************
!                                                           
!***********************************************************


qst(nx+2-1,1-2+1,indvarsst(20)) =  (qst(nx+2-1,1-2+1,indvarsst(4))+&
                    param_float(1 + 5)*q(nx+2-1,1-2+1,indvars(5))/(param_float(9 + 5)**2*qst(nx+2-1,1-2+1,indvarsst(2))**2))

