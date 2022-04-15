

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
! building source terms in RHS for layer 0 1 None d ********
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! d
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 0 1 None d *********************
!                                                           
!***********************************************************


qst(nx+2+0,ny+2-1,indvarsst(1)) =  qst(nx+2+0,ny+2-1,indvarsst(1))



!***********************************************************
!                                                           
! building source terms in RHS for layer 0 1 None eta ******
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! eta
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 0 1 None eta *******************
!                                                           
!***********************************************************


qst(nx+2+0,ny+2-1,indvarsst(2)) =  qst(nx+2+0,ny+2-1,indvarsst(2))



!***********************************************************
!                                                           
! building source terms in RHS for layer 0 1 None ksi ******
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! ksi
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 0 1 None ksi *******************
!                                                           
!***********************************************************


qst(nx+2+0,ny+2-1,indvarsst(3)) =  qst(nx+2+0,ny+2-1,indvarsst(3))



!***********************************************************
!                                                           
! building source terms in RHS for layer 0 1 None stemp ****
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (((0.5_wp*(deltayI*([u]_1y)-deltaxI*([v]_1x)))**2+(0.5_wp*(deltaxI*([v]_1x)-deltayI*([u]_1y)))**2)*2)**0.5
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_stemp_dx_0_nxp2p0p0nyp2m1k = q(nx+2+0+0,ny+2-1,indvars(3))

d1_stemp_dx_0_nxp2p0m1nyp2m1k = q(nx+2+0-1,ny+2-1,indvars(3))

d1_stemp_dx_0_nxp2p0m2nyp2m1k = q(nx+2+0-2,ny+2-1,indvars(3))

d1_stemp_dx_0_nxp2p0nyp2m1k = 1.5_wp*d1_stemp_dx_0_nxp2p0p0nyp2m1k-&
          2.0_wp*d1_stemp_dx_0_nxp2p0m1nyp2m1k+&
          0.5_wp*d1_stemp_dx_0_nxp2p0m2nyp2m1k

d1_stemp_dx_0_nxp2p0nyp2m1k = d1_stemp_dx_0_nxp2p0nyp2m1k*param_float(1)

d1_stemp_dy_0_nxp2p0nyp2m1m1k = q(nx+2+0,ny+2-1-1,indvars(2))

d1_stemp_dy_0_nxp2p0nyp2m1p1k = q(nx+2+0,ny+2-1+1,indvars(2))

d1_stemp_dy_0_nxp2p0nyp2m1k = -&
          0.5_wp*d1_stemp_dy_0_nxp2p0nyp2m1m1k+&
          0.5_wp*d1_stemp_dy_0_nxp2p0nyp2m1p1k

d1_stemp_dy_0_nxp2p0nyp2m1k = d1_stemp_dy_0_nxp2p0nyp2m1k*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer 0 1 None stemp *****************
!                                                           
!***********************************************************


qst(nx+2+0,ny+2-1,indvarsst(4)) =  (((0.5_wp*(qst(nx+2+0,ny+2-1,indvarsst(11))*(d1_stemp_dy_0_nxp2p0nyp2m1k)-&
                    qst(nx+2+0,ny+2-1,indvarsst(10))*(d1_stemp_dx_0_nxp2p0nyp2m1k)))**2+&
                    (0.5_wp*(qst(nx+2+0,ny+2-1,indvarsst(10))*(d1_stemp_dx_0_nxp2p0nyp2m1k)-&
                    qst(nx+2+0,ny+2-1,indvarsst(11))*(d1_stemp_dy_0_nxp2p0nyp2m1k)))**2)*2)**0.5



!***********************************************************
!                                                           
! building source terms in RHS for layer 0 1 None symm *****
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! ((sign(1.0_wp,ksi)-1.0_wp)/(-2.0_wp))
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 0 1 None symm ******************
!                                                           
!***********************************************************


qst(nx+2+0,ny+2-1,indvarsst(5)) =  ((sign(1.0_wp,qst(nx+2+0,ny+2-1,indvarsst(3)))-&
                    1.0_wp)/(-&
                    2.0_wp))



!***********************************************************
!                                                           
! building source terms in RHS for layer 0 1 None detady ***
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! [eta]_1y
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_detady_dy_0_nxp2p0nyp2m1m1k = qst(nx+2+0,ny+2-1-1,indvarsst(2))

d1_detady_dy_0_nxp2p0nyp2m1p1k = qst(nx+2+0,ny+2-1+1,indvarsst(2))

d1_detady_dy_0_nxp2p0nyp2m1k = -&
          0.5_wp*d1_detady_dy_0_nxp2p0nyp2m1m1k+&
          0.5_wp*d1_detady_dy_0_nxp2p0nyp2m1p1k

d1_detady_dy_0_nxp2p0nyp2m1k = d1_detady_dy_0_nxp2p0nyp2m1k*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer 0 1 None detady ****************
!                                                           
!***********************************************************


qst(nx+2+0,ny+2-1,indvarsst(6)) =  d1_detady_dy_0_nxp2p0nyp2m1k



!***********************************************************
!                                                           
! building source terms in RHS for layer 0 1 None dksidy ***
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! [ksi]_1y
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_dksidy_dy_0_nxp2p0nyp2m1m1k = qst(nx+2+0,ny+2-1-1,indvarsst(3))

d1_dksidy_dy_0_nxp2p0nyp2m1p1k = qst(nx+2+0,ny+2-1+1,indvarsst(3))

d1_dksidy_dy_0_nxp2p0nyp2m1k = -&
          0.5_wp*d1_dksidy_dy_0_nxp2p0nyp2m1m1k+&
          0.5_wp*d1_dksidy_dy_0_nxp2p0nyp2m1p1k

d1_dksidy_dy_0_nxp2p0nyp2m1k = d1_dksidy_dy_0_nxp2p0nyp2m1k*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer 0 1 None dksidy ****************
!                                                           
!***********************************************************


qst(nx+2+0,ny+2-1,indvarsst(7)) =  d1_dksidy_dy_0_nxp2p0nyp2m1k



!***********************************************************
!                                                           
! building source terms in RHS for layer 0 1 None detadx ***
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! [eta]_1x
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_detadx_dx_0_nxp2p0p0nyp2m1k = qst(nx+2+0+0,ny+2-1,indvarsst(2))

d1_detadx_dx_0_nxp2p0m1nyp2m1k = qst(nx+2+0-1,ny+2-1,indvarsst(2))

d1_detadx_dx_0_nxp2p0m2nyp2m1k = qst(nx+2+0-2,ny+2-1,indvarsst(2))

d1_detadx_dx_0_nxp2p0nyp2m1k = 1.5_wp*d1_detadx_dx_0_nxp2p0p0nyp2m1k-&
          2.0_wp*d1_detadx_dx_0_nxp2p0m1nyp2m1k+&
          0.5_wp*d1_detadx_dx_0_nxp2p0m2nyp2m1k

d1_detadx_dx_0_nxp2p0nyp2m1k = d1_detadx_dx_0_nxp2p0nyp2m1k*param_float(1)



!***********************************************************
!                                                           
! Update BC terms for layer 0 1 None detadx ****************
!                                                           
!***********************************************************


qst(nx+2+0,ny+2-1,indvarsst(8)) =  d1_detadx_dx_0_nxp2p0nyp2m1k



!***********************************************************
!                                                           
! building source terms in RHS for layer 0 1 None dksidx ***
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! ([ksi]_1x)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_dksidx_dx_0_nxp2p0p0nyp2m1k = qst(nx+2+0+0,ny+2-1,indvarsst(3))

d1_dksidx_dx_0_nxp2p0m1nyp2m1k = qst(nx+2+0-1,ny+2-1,indvarsst(3))

d1_dksidx_dx_0_nxp2p0m2nyp2m1k = qst(nx+2+0-2,ny+2-1,indvarsst(3))

d1_dksidx_dx_0_nxp2p0nyp2m1k = 1.5_wp*d1_dksidx_dx_0_nxp2p0p0nyp2m1k-&
          2.0_wp*d1_dksidx_dx_0_nxp2p0m1nyp2m1k+&
          0.5_wp*d1_dksidx_dx_0_nxp2p0m2nyp2m1k

d1_dksidx_dx_0_nxp2p0nyp2m1k = d1_dksidx_dx_0_nxp2p0nyp2m1k*param_float(1)



!***********************************************************
!                                                           
! Update BC terms for layer 0 1 None dksidx ****************
!                                                           
!***********************************************************


qst(nx+2+0,ny+2-1,indvarsst(9)) =  (d1_dksidx_dx_0_nxp2p0nyp2m1k)



!***********************************************************
!                                                           
! building source terms in RHS for layer 0 1 None deltaxI **
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! 1.0_wp/(dksidx)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 0 1 None deltaxI ***************
!                                                           
!***********************************************************


qst(nx+2+0,ny+2-1,indvarsst(10)) =  1.0_wp/(qst(nx+2+0,ny+2-1,indvarsst(9)))



!***********************************************************
!                                                           
! building source terms in RHS for layer 0 1 None deltayI **
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! 1.0_wp/(detady)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 0 1 None deltayI ***************
!                                                           
!***********************************************************


qst(nx+2+0,ny+2-1,indvarsst(11)) =  1.0_wp/(qst(nx+2+0,ny+2-1,indvarsst(6)))



!***********************************************************
!                                                           
! building source terms in RHS for layer 0 1 None viscosità 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (1+sut)/(T+sut)*T**1.5
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 0 1 None viscosità *************
!                                                           
!***********************************************************


qst(nx+2+0,ny+2-1,indvarsst(12)) =  (1+&
                    param_float(21 + 5))/(((q(nx+2+0,ny+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2-1,indvars(2))*q(nx+2+0,ny+2-1,indvars(2))+&
                    q(nx+2+0,ny+2-1,indvars(3))*q(nx+2+0,ny+2-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2+0,ny+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2-1,indvars(2))*q(nx+2+0,ny+2-1,indvars(2))+&
                    q(nx+2+0,ny+2-1,indvars(3))*q(nx+2+0,ny+2-1,indvars(3)))))/param_float(4 + 5)**1.5



!***********************************************************
!                                                           
! building source terms in RHS for layer 0 1 None fw *******
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! gg*(1+Cw3**6)/(gg**6+Cw3**6)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 0 1 None fw ********************
!                                                           
!***********************************************************


qst(nx+2+0,ny+2-1,indvarsst(13)) =  ((q(nx+2+0,ny+2-1,indvars(5))/((qst(nx+2+0,ny+2-1,indvarsst(4))+&
                    param_float(1 + 5)*q(nx+2+0,ny+2-1,indvars(5))/(param_float(9 + 5)**2*qst(nx+2+0,ny+2-1,indvarsst(2))**2))*param_float(9 + 5)**2*qst(nx+2+0,ny+2-1,indvarsst(2))**2))+&
                    param_float(11 + 5)*((q(nx+2+0,ny+2-1,indvars(5))/((qst(nx+2+0,ny+2-1,indvarsst(4))+&
                    param_float(1 + 5)*q(nx+2+0,ny+2-1,indvars(5))/(param_float(9 + 5)**2*qst(nx+2+0,ny+2-1,indvarsst(2))**2))*param_float(9 + 5)**2*qst(nx+2+0,ny+2-1,indvarsst(2))**2))**6-&
                    (q(nx+2+0,ny+2-1,indvars(5))/((qst(nx+2+0,ny+2-1,indvarsst(4))+&
                    param_float(1 + 5)*q(nx+2+0,ny+2-1,indvars(5))/(param_float(9 + 5)**2*qst(nx+2+0,ny+2-1,indvarsst(2))**2))*param_float(9 + 5)**2*qst(nx+2+0,ny+2-1,indvarsst(2))**2))))*(1+&
                    param_float(12 + 5)**6)/(((q(nx+2+0,ny+2-1,indvars(5))/((qst(nx+2+0,ny+2-1,indvarsst(4))+&
                    param_float(1 + 5)*q(nx+2+0,ny+2-1,indvars(5))/(param_float(9 + 5)**2*qst(nx+2+0,ny+2-1,indvarsst(2))**2))*param_float(9 + 5)**2*qst(nx+2+0,ny+2-1,indvarsst(2))**2))+&
                    param_float(11 + 5)*((q(nx+2+0,ny+2-1,indvars(5))/((qst(nx+2+0,ny+2-1,indvarsst(4))+&
                    param_float(1 + 5)*q(nx+2+0,ny+2-1,indvars(5))/(param_float(9 + 5)**2*qst(nx+2+0,ny+2-1,indvarsst(2))**2))*param_float(9 + 5)**2*qst(nx+2+0,ny+2-1,indvarsst(2))**2))**6-&
                    (q(nx+2+0,ny+2-1,indvars(5))/((qst(nx+2+0,ny+2-1,indvarsst(4))+&
                    param_float(1 + 5)*q(nx+2+0,ny+2-1,indvars(5))/(param_float(9 + 5)**2*qst(nx+2+0,ny+2-1,indvarsst(2))**2))*param_float(9 + 5)**2*qst(nx+2+0,ny+2-1,indvarsst(2))**2))))**6+&
                    param_float(12 + 5)**6)

