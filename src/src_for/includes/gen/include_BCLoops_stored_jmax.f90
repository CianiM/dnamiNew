

!***********************************************************
!                                                           
! Start building layers for BC : None jmax None ************
!                                                           
!***********************************************************




!***********************************************************
!                                                           
! BC layer: None 0 None ************************************
!                                                           
!***********************************************************


 
      do i=idloop(1),idloop(2) 


!***********************************************************
!                                                           
! building source terms in RHS for layer None 0 None stemp *
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (((0.5_wp*(deltayI*([u]_1y)-deltaxI*([v]_1x)))**2+(0.5_wp*(deltaxI*([v]_1x)-deltayI*([u]_1y)))**2)*2)**0.5
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_stemp_dx_0_im1nyp2p0k = q(i-1,ny+2+0,indvars(3))

d1_stemp_dx_0_ip1nyp2p0k = q(i+1,ny+2+0,indvars(3))

d1_stemp_dx_0_inyp2p0k = -&
          0.5_wp*d1_stemp_dx_0_im1nyp2p0k+&
          0.5_wp*d1_stemp_dx_0_ip1nyp2p0k

d1_stemp_dx_0_inyp2p0k = d1_stemp_dx_0_inyp2p0k*param_float(1)

d1_stemp_dy_0_inyp2p0p0k = q(i,ny+2+0+0,indvars(2))

d1_stemp_dy_0_inyp2p0m1k = q(i,ny+2+0-1,indvars(2))

d1_stemp_dy_0_inyp2p0m2k = q(i,ny+2+0-2,indvars(2))

d1_stemp_dy_0_inyp2p0k = 1.5_wp*d1_stemp_dy_0_inyp2p0p0k-&
          2.0_wp*d1_stemp_dy_0_inyp2p0m1k+&
          0.5_wp*d1_stemp_dy_0_inyp2p0m2k

d1_stemp_dy_0_inyp2p0k = d1_stemp_dy_0_inyp2p0k*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer None 0 None stemp **************
!                                                           
!***********************************************************


qst(i,ny+2+0,indvarsst(4)) =  (((0.5_wp*(qst(i,ny+2+0,indvarsst(11))*(d1_stemp_dy_0_inyp2p0k)-&
                    qst(i,ny+2+0,indvarsst(10))*(d1_stemp_dx_0_inyp2p0k)))**2+&
                    (0.5_wp*(qst(i,ny+2+0,indvarsst(10))*(d1_stemp_dx_0_inyp2p0k)-&
                    qst(i,ny+2+0,indvarsst(11))*(d1_stemp_dy_0_inyp2p0k)))**2)*2)**0.5



!***********************************************************
!                                                           
! building source terms in RHS for layer None 0 None viscosità 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (1+sut)/(T+sut)*T**1.5
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer None 0 None viscosità **********
!                                                           
!***********************************************************


qst(i,ny+2+0,indvarsst(12)) =  (1+&
                    param_float(21 + 5))/(((q(i,ny+2+0,indvars(4))-&
                    0.5_wp*(q(i,ny+2+0,indvars(2))*q(i,ny+2+0,indvars(2))+&
                    q(i,ny+2+0,indvars(3))*q(i,ny+2+0,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(i,ny+2+0,indvars(4))-&
                    0.5_wp*(q(i,ny+2+0,indvars(2))*q(i,ny+2+0,indvars(2))+&
                    q(i,ny+2+0,indvars(3))*q(i,ny+2+0,indvars(3)))))/param_float(4 + 5)**1.5



!***********************************************************
!                                                           
! building source terms in RHS for layer None 0 None fw ****
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! gg*(1+Cw3**6)/(gg**6+Cw3**6)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer None 0 None fw *****************
!                                                           
!***********************************************************


qst(i,ny+2+0,indvarsst(13)) =  qst(i,ny+2+0,indvarsst(14))*(1+&
                    param_float(12 + 5)**6)/(qst(i,ny+2+0,indvarsst(14))**6+&
                    param_float(12 + 5)**6)



!***********************************************************
!                                                           
! building source terms in RHS for layer None 0 None gg ****
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (nut/(SS*k**2*eta**2))
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer None 0 None gg *****************
!                                                           
!***********************************************************


qst(i,ny+2+0,indvarsst(14)) =  (q(i,ny+2+0,indvars(5))/(qst(i,ny+2+0,indvarsst(20))*param_float(9 + 5)**2*qst(i,ny+2+0,indvarsst(2))**2))



!***********************************************************
!                                                           
! building source terms in RHS for layer None 0 None fv2 ***
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (1-chi/(1+chi*fv1))
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer None 0 None fv2 ****************
!                                                           
!***********************************************************


qst(i,ny+2+0,indvarsst(15)) =  (1-&
                    (q(i,ny+2+0,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(i,ny+2+0,indvars(4))-&
                    0.5_wp*(q(i,ny+2+0,indvars(2))*q(i,ny+2+0,indvars(2))+&
                    q(i,ny+2+0,indvars(3))*q(i,ny+2+0,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(i,ny+2+0,indvars(4))-&
                    0.5_wp*(q(i,ny+2+0,indvars(2))*q(i,ny+2+0,indvars(2))+&
                    q(i,ny+2+0,indvars(3))*q(i,ny+2+0,indvars(3)))))/param_float(4 + 5)**1.5*q(i,ny+2+0,indvars(1)))/(1+&
                    (q(i,ny+2+0,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(i,ny+2+0,indvars(4))-&
                    0.5_wp*(q(i,ny+2+0,indvars(2))*q(i,ny+2+0,indvars(2))+&
                    q(i,ny+2+0,indvars(3))*q(i,ny+2+0,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(i,ny+2+0,indvars(4))-&
                    0.5_wp*(q(i,ny+2+0,indvars(2))*q(i,ny+2+0,indvars(2))+&
                    q(i,ny+2+0,indvars(3))*q(i,ny+2+0,indvars(3)))))/param_float(4 + 5)**1.5*q(i,ny+2+0,indvars(1)))*((q(i,ny+2+0,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(i,ny+2+0,indvars(4))-&
                    0.5_wp*(q(i,ny+2+0,indvars(2))*q(i,ny+2+0,indvars(2))+&
                    q(i,ny+2+0,indvars(3))*q(i,ny+2+0,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(i,ny+2+0,indvars(4))-&
                    0.5_wp*(q(i,ny+2+0,indvars(2))*q(i,ny+2+0,indvars(2))+&
                    q(i,ny+2+0,indvars(3))*q(i,ny+2+0,indvars(3)))))/param_float(4 + 5)**1.5*q(i,ny+2+0,indvars(1)))**3/((q(i,ny+2+0,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(i,ny+2+0,indvars(4))-&
                    0.5_wp*(q(i,ny+2+0,indvars(2))*q(i,ny+2+0,indvars(2))+&
                    q(i,ny+2+0,indvars(3))*q(i,ny+2+0,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(i,ny+2+0,indvars(4))-&
                    0.5_wp*(q(i,ny+2+0,indvars(2))*q(i,ny+2+0,indvars(2))+&
                    q(i,ny+2+0,indvars(3))*q(i,ny+2+0,indvars(3)))))/param_float(4 + 5)**1.5*q(i,ny+2+0,indvars(1)))**3+&
                    param_float(13 + 5)**3))))



!***********************************************************
!                                                           
! building source terms in RHS for layer None 0 None ft2 ***
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (Ct3*exp(-Ct4*chi**2))
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer None 0 None ft2 ****************
!                                                           
!***********************************************************


qst(i,ny+2+0,indvarsst(16)) =  (param_float(16 + 5)*exp(-&
                    param_float(17 + 5)*(q(i,ny+2+0,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(i,ny+2+0,indvars(4))-&
                    0.5_wp*(q(i,ny+2+0,indvars(2))*q(i,ny+2+0,indvars(2))+&
                    q(i,ny+2+0,indvars(3))*q(i,ny+2+0,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(i,ny+2+0,indvars(4))-&
                    0.5_wp*(q(i,ny+2+0,indvars(2))*q(i,ny+2+0,indvars(2))+&
                    q(i,ny+2+0,indvars(3))*q(i,ny+2+0,indvars(3)))))/param_float(4 + 5)**1.5*q(i,ny+2+0,indvars(1)))**2))



!***********************************************************
!                                                           
! building source terms in RHS for layer None 0 None fw2 ***
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (gg*(1+Cw3**6)/(gg**6+Cw3**6))
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer None 0 None fw2 ****************
!                                                           
!***********************************************************


qst(i,ny+2+0,indvarsst(17)) =  (qst(i,ny+2+0,indvarsst(14))*(1+&
                    param_float(12 + 5)**6)/(qst(i,ny+2+0,indvarsst(14))**6+&
                    param_float(12 + 5)**6))



!***********************************************************
!                                                           
! building source terms in RHS for layer None 0 None rr ****
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (nut/(SS*k**2*eta**2))
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer None 0 None rr *****************
!                                                           
!***********************************************************


qst(i,ny+2+0,indvarsst(18)) =  (q(i,ny+2+0,indvars(5))/(qst(i,ny+2+0,indvarsst(20))*param_float(9 + 5)**2*qst(i,ny+2+0,indvarsst(2))**2))



!***********************************************************
!                                                           
! building source terms in RHS for layer None 0 None SS ****
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (stemp+ReI*nut/(k**2*eta**2))
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer None 0 None SS *****************
!                                                           
!***********************************************************


qst(i,ny+2+0,indvarsst(20)) =  (qst(i,ny+2+0,indvarsst(4))+&
                    param_float(1 + 5)*q(i,ny+2+0,indvars(5))/(param_float(9 + 5)**2*qst(i,ny+2+0,indvarsst(2))**2))

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
! building source terms in RHS for layer None 1 None stemp *
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (((0.5_wp*(deltayI*([u]_1y)-deltaxI*([v]_1x)))**2+(0.5_wp*(deltaxI*([v]_1x)-deltayI*([u]_1y)))**2)*2)**0.5
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_stemp_dx_0_im1nyp2m1k = q(i-1,ny+2-1,indvars(3))

d1_stemp_dx_0_ip1nyp2m1k = q(i+1,ny+2-1,indvars(3))

d1_stemp_dx_0_inyp2m1k = -&
          0.5_wp*d1_stemp_dx_0_im1nyp2m1k+&
          0.5_wp*d1_stemp_dx_0_ip1nyp2m1k

d1_stemp_dx_0_inyp2m1k = d1_stemp_dx_0_inyp2m1k*param_float(1)

d1_stemp_dy_0_inyp2m1m1k = q(i,ny+2-1-1,indvars(2))

d1_stemp_dy_0_inyp2m1p1k = q(i,ny+2-1+1,indvars(2))

d1_stemp_dy_0_inyp2m1k = -&
          0.5_wp*d1_stemp_dy_0_inyp2m1m1k+&
          0.5_wp*d1_stemp_dy_0_inyp2m1p1k

d1_stemp_dy_0_inyp2m1k = d1_stemp_dy_0_inyp2m1k*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer None 1 None stemp **************
!                                                           
!***********************************************************


qst(i,ny+2-1,indvarsst(4)) =  (((0.5_wp*(qst(i,ny+2-1,indvarsst(11))*(d1_stemp_dy_0_inyp2m1k)-&
                    qst(i,ny+2-1,indvarsst(10))*(d1_stemp_dx_0_inyp2m1k)))**2+&
                    (0.5_wp*(qst(i,ny+2-1,indvarsst(10))*(d1_stemp_dx_0_inyp2m1k)-&
                    qst(i,ny+2-1,indvarsst(11))*(d1_stemp_dy_0_inyp2m1k)))**2)*2)**0.5



!***********************************************************
!                                                           
! building source terms in RHS for layer None 1 None viscosità 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (1+sut)/(T+sut)*T**1.5
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer None 1 None viscosità **********
!                                                           
!***********************************************************


qst(i,ny+2-1,indvarsst(12)) =  (1+&
                    param_float(21 + 5))/(((q(i,ny+2-1,indvars(4))-&
                    0.5_wp*(q(i,ny+2-1,indvars(2))*q(i,ny+2-1,indvars(2))+&
                    q(i,ny+2-1,indvars(3))*q(i,ny+2-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(i,ny+2-1,indvars(4))-&
                    0.5_wp*(q(i,ny+2-1,indvars(2))*q(i,ny+2-1,indvars(2))+&
                    q(i,ny+2-1,indvars(3))*q(i,ny+2-1,indvars(3)))))/param_float(4 + 5)**1.5



!***********************************************************
!                                                           
! building source terms in RHS for layer None 1 None fw ****
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! gg*(1+Cw3**6)/(gg**6+Cw3**6)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer None 1 None fw *****************
!                                                           
!***********************************************************


qst(i,ny+2-1,indvarsst(13)) =  qst(i,ny+2-1,indvarsst(14))*(1+&
                    param_float(12 + 5)**6)/(qst(i,ny+2-1,indvarsst(14))**6+&
                    param_float(12 + 5)**6)



!***********************************************************
!                                                           
! building source terms in RHS for layer None 1 None gg ****
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (nut/(SS*k**2*eta**2))
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer None 1 None gg *****************
!                                                           
!***********************************************************


qst(i,ny+2-1,indvarsst(14)) =  (q(i,ny+2-1,indvars(5))/(qst(i,ny+2-1,indvarsst(20))*param_float(9 + 5)**2*qst(i,ny+2-1,indvarsst(2))**2))



!***********************************************************
!                                                           
! building source terms in RHS for layer None 1 None fv2 ***
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (1-chi/(1+chi*fv1))
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer None 1 None fv2 ****************
!                                                           
!***********************************************************


qst(i,ny+2-1,indvarsst(15)) =  (1-&
                    (q(i,ny+2-1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(i,ny+2-1,indvars(4))-&
                    0.5_wp*(q(i,ny+2-1,indvars(2))*q(i,ny+2-1,indvars(2))+&
                    q(i,ny+2-1,indvars(3))*q(i,ny+2-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(i,ny+2-1,indvars(4))-&
                    0.5_wp*(q(i,ny+2-1,indvars(2))*q(i,ny+2-1,indvars(2))+&
                    q(i,ny+2-1,indvars(3))*q(i,ny+2-1,indvars(3)))))/param_float(4 + 5)**1.5*q(i,ny+2-1,indvars(1)))/(1+&
                    (q(i,ny+2-1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(i,ny+2-1,indvars(4))-&
                    0.5_wp*(q(i,ny+2-1,indvars(2))*q(i,ny+2-1,indvars(2))+&
                    q(i,ny+2-1,indvars(3))*q(i,ny+2-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(i,ny+2-1,indvars(4))-&
                    0.5_wp*(q(i,ny+2-1,indvars(2))*q(i,ny+2-1,indvars(2))+&
                    q(i,ny+2-1,indvars(3))*q(i,ny+2-1,indvars(3)))))/param_float(4 + 5)**1.5*q(i,ny+2-1,indvars(1)))*((q(i,ny+2-1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(i,ny+2-1,indvars(4))-&
                    0.5_wp*(q(i,ny+2-1,indvars(2))*q(i,ny+2-1,indvars(2))+&
                    q(i,ny+2-1,indvars(3))*q(i,ny+2-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(i,ny+2-1,indvars(4))-&
                    0.5_wp*(q(i,ny+2-1,indvars(2))*q(i,ny+2-1,indvars(2))+&
                    q(i,ny+2-1,indvars(3))*q(i,ny+2-1,indvars(3)))))/param_float(4 + 5)**1.5*q(i,ny+2-1,indvars(1)))**3/((q(i,ny+2-1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(i,ny+2-1,indvars(4))-&
                    0.5_wp*(q(i,ny+2-1,indvars(2))*q(i,ny+2-1,indvars(2))+&
                    q(i,ny+2-1,indvars(3))*q(i,ny+2-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(i,ny+2-1,indvars(4))-&
                    0.5_wp*(q(i,ny+2-1,indvars(2))*q(i,ny+2-1,indvars(2))+&
                    q(i,ny+2-1,indvars(3))*q(i,ny+2-1,indvars(3)))))/param_float(4 + 5)**1.5*q(i,ny+2-1,indvars(1)))**3+&
                    param_float(13 + 5)**3))))



!***********************************************************
!                                                           
! building source terms in RHS for layer None 1 None ft2 ***
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (Ct3*exp(-Ct4*chi**2))
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer None 1 None ft2 ****************
!                                                           
!***********************************************************


qst(i,ny+2-1,indvarsst(16)) =  (param_float(16 + 5)*exp(-&
                    param_float(17 + 5)*(q(i,ny+2-1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(i,ny+2-1,indvars(4))-&
                    0.5_wp*(q(i,ny+2-1,indvars(2))*q(i,ny+2-1,indvars(2))+&
                    q(i,ny+2-1,indvars(3))*q(i,ny+2-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(i,ny+2-1,indvars(4))-&
                    0.5_wp*(q(i,ny+2-1,indvars(2))*q(i,ny+2-1,indvars(2))+&
                    q(i,ny+2-1,indvars(3))*q(i,ny+2-1,indvars(3)))))/param_float(4 + 5)**1.5*q(i,ny+2-1,indvars(1)))**2))



!***********************************************************
!                                                           
! building source terms in RHS for layer None 1 None fw2 ***
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (gg*(1+Cw3**6)/(gg**6+Cw3**6))
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer None 1 None fw2 ****************
!                                                           
!***********************************************************


qst(i,ny+2-1,indvarsst(17)) =  (qst(i,ny+2-1,indvarsst(14))*(1+&
                    param_float(12 + 5)**6)/(qst(i,ny+2-1,indvarsst(14))**6+&
                    param_float(12 + 5)**6))



!***********************************************************
!                                                           
! building source terms in RHS for layer None 1 None rr ****
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (nut/(SS*k**2*eta**2))
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer None 1 None rr *****************
!                                                           
!***********************************************************


qst(i,ny+2-1,indvarsst(18)) =  (q(i,ny+2-1,indvars(5))/(qst(i,ny+2-1,indvarsst(20))*param_float(9 + 5)**2*qst(i,ny+2-1,indvarsst(2))**2))



!***********************************************************
!                                                           
! building source terms in RHS for layer None 1 None SS ****
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (stemp+ReI*nut/(k**2*eta**2))
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer None 1 None SS *****************
!                                                           
!***********************************************************


qst(i,ny+2-1,indvarsst(20)) =  (qst(i,ny+2-1,indvarsst(4))+&
                    param_float(1 + 5)*q(i,ny+2-1,indvars(5))/(param_float(9 + 5)**2*qst(i,ny+2-1,indvarsst(2))**2))

   enddo
