

!***********************************************************
!                                                           
! Start building layers for BC : imax None None ************
!                                                           
!***********************************************************




!***********************************************************
!                                                           
! BC layer: 0 None None ************************************
!                                                           
!***********************************************************


     do j=idloop(3),idloop(4) 


!***********************************************************
!                                                           
! building source terms in RHS for layer 0 None None d *****
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! d
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 0 None None d ******************
!                                                           
!***********************************************************


qst(nx+2+0,j,indvarsst(1)) =  qst(nx+2+0,j,indvarsst(1))



!***********************************************************
!                                                           
! building source terms in RHS for layer 0 None None eta ***
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! eta
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 0 None None eta ****************
!                                                           
!***********************************************************


qst(nx+2+0,j,indvarsst(2)) =  qst(nx+2+0,j,indvarsst(2))



!***********************************************************
!                                                           
! building source terms in RHS for layer 0 None None ksi ***
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! ksi
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 0 None None ksi ****************
!                                                           
!***********************************************************


qst(nx+2+0,j,indvarsst(3)) =  qst(nx+2+0,j,indvarsst(3))



!***********************************************************
!                                                           
! building source terms in RHS for layer 0 None None stemp *
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (((0.5_wp*(deltayI*([u]_1y)-deltaxI*([v]_1x)))**2+(0.5_wp*(deltaxI*([v]_1x)-deltayI*([u]_1y)))**2)*2)**0.5
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_stemp_dx_0_nxp2p0p0jk = q(nx+2+0+0,j,indvars(3))

d1_stemp_dx_0_nxp2p0m1jk = q(nx+2+0-1,j,indvars(3))

d1_stemp_dx_0_nxp2p0m2jk = q(nx+2+0-2,j,indvars(3))

d1_stemp_dx_0_nxp2p0jk = 1.5_wp*d1_stemp_dx_0_nxp2p0p0jk-&
          2.0_wp*d1_stemp_dx_0_nxp2p0m1jk+&
          0.5_wp*d1_stemp_dx_0_nxp2p0m2jk

d1_stemp_dx_0_nxp2p0jk = d1_stemp_dx_0_nxp2p0jk*param_float(1)

d1_stemp_dy_0_nxp2p0jm1k = q(nx+2+0,j-1,indvars(2))

d1_stemp_dy_0_nxp2p0jp1k = q(nx+2+0,j+1,indvars(2))

d1_stemp_dy_0_nxp2p0jk = -&
          0.5_wp*d1_stemp_dy_0_nxp2p0jm1k+&
          0.5_wp*d1_stemp_dy_0_nxp2p0jp1k

d1_stemp_dy_0_nxp2p0jk = d1_stemp_dy_0_nxp2p0jk*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer 0 None None stemp **************
!                                                           
!***********************************************************


qst(nx+2+0,j,indvarsst(4)) =  (((0.5_wp*(qst(nx+2+0,j,indvarsst(11))*(d1_stemp_dy_0_nxp2p0jk)-&
                    qst(nx+2+0,j,indvarsst(10))*(d1_stemp_dx_0_nxp2p0jk)))**2+&
                    (0.5_wp*(qst(nx+2+0,j,indvarsst(10))*(d1_stemp_dx_0_nxp2p0jk)-&
                    qst(nx+2+0,j,indvarsst(11))*(d1_stemp_dy_0_nxp2p0jk)))**2)*2)**0.5



!***********************************************************
!                                                           
! building source terms in RHS for layer 0 None None symm **
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! ((sign(1.0_wp,ksi)-1.0_wp)/(-2.0_wp))
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 0 None None symm ***************
!                                                           
!***********************************************************


qst(nx+2+0,j,indvarsst(5)) =  ((sign(1.0_wp,qst(nx+2+0,j,indvarsst(3)))-&
                    1.0_wp)/(-&
                    2.0_wp))



!***********************************************************
!                                                           
! building source terms in RHS for layer 0 None None detady 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! [eta]_1y
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_detady_dy_0_nxp2p0jm1k = qst(nx+2+0,j-1,indvarsst(2))

d1_detady_dy_0_nxp2p0jp1k = qst(nx+2+0,j+1,indvarsst(2))

d1_detady_dy_0_nxp2p0jk = -&
          0.5_wp*d1_detady_dy_0_nxp2p0jm1k+&
          0.5_wp*d1_detady_dy_0_nxp2p0jp1k

d1_detady_dy_0_nxp2p0jk = d1_detady_dy_0_nxp2p0jk*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer 0 None None detady *************
!                                                           
!***********************************************************


qst(nx+2+0,j,indvarsst(6)) =  d1_detady_dy_0_nxp2p0jk



!***********************************************************
!                                                           
! building source terms in RHS for layer 0 None None dksidy 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! [ksi]_1y
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_dksidy_dy_0_nxp2p0jm1k = qst(nx+2+0,j-1,indvarsst(3))

d1_dksidy_dy_0_nxp2p0jp1k = qst(nx+2+0,j+1,indvarsst(3))

d1_dksidy_dy_0_nxp2p0jk = -&
          0.5_wp*d1_dksidy_dy_0_nxp2p0jm1k+&
          0.5_wp*d1_dksidy_dy_0_nxp2p0jp1k

d1_dksidy_dy_0_nxp2p0jk = d1_dksidy_dy_0_nxp2p0jk*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer 0 None None dksidy *************
!                                                           
!***********************************************************


qst(nx+2+0,j,indvarsst(7)) =  d1_dksidy_dy_0_nxp2p0jk



!***********************************************************
!                                                           
! building source terms in RHS for layer 0 None None detadx 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! [eta]_1x
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_detadx_dx_0_nxp2p0p0jk = qst(nx+2+0+0,j,indvarsst(2))

d1_detadx_dx_0_nxp2p0m1jk = qst(nx+2+0-1,j,indvarsst(2))

d1_detadx_dx_0_nxp2p0m2jk = qst(nx+2+0-2,j,indvarsst(2))

d1_detadx_dx_0_nxp2p0jk = 1.5_wp*d1_detadx_dx_0_nxp2p0p0jk-&
          2.0_wp*d1_detadx_dx_0_nxp2p0m1jk+&
          0.5_wp*d1_detadx_dx_0_nxp2p0m2jk

d1_detadx_dx_0_nxp2p0jk = d1_detadx_dx_0_nxp2p0jk*param_float(1)



!***********************************************************
!                                                           
! Update BC terms for layer 0 None None detadx *************
!                                                           
!***********************************************************


qst(nx+2+0,j,indvarsst(8)) =  d1_detadx_dx_0_nxp2p0jk



!***********************************************************
!                                                           
! building source terms in RHS for layer 0 None None dksidx 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! ([ksi]_1x)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_dksidx_dx_0_nxp2p0p0jk = qst(nx+2+0+0,j,indvarsst(3))

d1_dksidx_dx_0_nxp2p0m1jk = qst(nx+2+0-1,j,indvarsst(3))

d1_dksidx_dx_0_nxp2p0m2jk = qst(nx+2+0-2,j,indvarsst(3))

d1_dksidx_dx_0_nxp2p0jk = 1.5_wp*d1_dksidx_dx_0_nxp2p0p0jk-&
          2.0_wp*d1_dksidx_dx_0_nxp2p0m1jk+&
          0.5_wp*d1_dksidx_dx_0_nxp2p0m2jk

d1_dksidx_dx_0_nxp2p0jk = d1_dksidx_dx_0_nxp2p0jk*param_float(1)



!***********************************************************
!                                                           
! Update BC terms for layer 0 None None dksidx *************
!                                                           
!***********************************************************


qst(nx+2+0,j,indvarsst(9)) =  (d1_dksidx_dx_0_nxp2p0jk)



!***********************************************************
!                                                           
! building source terms in RHS for layer 0 None None deltaxI 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! 1.0_wp/(dksidx)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 0 None None deltaxI ************
!                                                           
!***********************************************************


qst(nx+2+0,j,indvarsst(10)) =  1.0_wp/(qst(nx+2+0,j,indvarsst(9)))



!***********************************************************
!                                                           
! building source terms in RHS for layer 0 None None deltayI 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! 1.0_wp/(detady)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 0 None None deltayI ************
!                                                           
!***********************************************************


qst(nx+2+0,j,indvarsst(11)) =  1.0_wp/(qst(nx+2+0,j,indvarsst(6)))



!***********************************************************
!                                                           
! building source terms in RHS for layer 0 None None viscosità 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (1+sut)/(T+sut)*T**1.5
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 0 None None viscosità **********
!                                                           
!***********************************************************


qst(nx+2+0,j,indvarsst(12)) =  (1+&
                    param_float(21 + 5))/(((q(nx+2+0,j,indvars(4))-&
                    0.5_wp*(q(nx+2+0,j,indvars(2))*q(nx+2+0,j,indvars(2))+&
                    q(nx+2+0,j,indvars(3))*q(nx+2+0,j,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2+0,j,indvars(4))-&
                    0.5_wp*(q(nx+2+0,j,indvars(2))*q(nx+2+0,j,indvars(2))+&
                    q(nx+2+0,j,indvars(3))*q(nx+2+0,j,indvars(3)))))/param_float(4 + 5)**1.5



!***********************************************************
!                                                           
! building source terms in RHS for layer 0 None None fw ****
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! gg*(1+Cw3**6)/(gg**6+Cw3**6)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 0 None None fw *****************
!                                                           
!***********************************************************


qst(nx+2+0,j,indvarsst(13)) =  ((q(nx+2+0,j,indvars(5))/((qst(nx+2+0,j,indvarsst(4))+&
                    param_float(1 + 5)*q(nx+2+0,j,indvars(5))/(param_float(9 + 5)**2*qst(nx+2+0,j,indvarsst(2))**2))*param_float(9 + 5)**2*qst(nx+2+0,j,indvarsst(2))**2))+&
                    param_float(11 + 5)*((q(nx+2+0,j,indvars(5))/((qst(nx+2+0,j,indvarsst(4))+&
                    param_float(1 + 5)*q(nx+2+0,j,indvars(5))/(param_float(9 + 5)**2*qst(nx+2+0,j,indvarsst(2))**2))*param_float(9 + 5)**2*qst(nx+2+0,j,indvarsst(2))**2))**6-&
                    (q(nx+2+0,j,indvars(5))/((qst(nx+2+0,j,indvarsst(4))+&
                    param_float(1 + 5)*q(nx+2+0,j,indvars(5))/(param_float(9 + 5)**2*qst(nx+2+0,j,indvarsst(2))**2))*param_float(9 + 5)**2*qst(nx+2+0,j,indvarsst(2))**2))))*(1+&
                    param_float(12 + 5)**6)/(((q(nx+2+0,j,indvars(5))/((qst(nx+2+0,j,indvarsst(4))+&
                    param_float(1 + 5)*q(nx+2+0,j,indvars(5))/(param_float(9 + 5)**2*qst(nx+2+0,j,indvarsst(2))**2))*param_float(9 + 5)**2*qst(nx+2+0,j,indvarsst(2))**2))+&
                    param_float(11 + 5)*((q(nx+2+0,j,indvars(5))/((qst(nx+2+0,j,indvarsst(4))+&
                    param_float(1 + 5)*q(nx+2+0,j,indvars(5))/(param_float(9 + 5)**2*qst(nx+2+0,j,indvarsst(2))**2))*param_float(9 + 5)**2*qst(nx+2+0,j,indvarsst(2))**2))**6-&
                    (q(nx+2+0,j,indvars(5))/((qst(nx+2+0,j,indvarsst(4))+&
                    param_float(1 + 5)*q(nx+2+0,j,indvars(5))/(param_float(9 + 5)**2*qst(nx+2+0,j,indvarsst(2))**2))*param_float(9 + 5)**2*qst(nx+2+0,j,indvarsst(2))**2))))**6+&
                    param_float(12 + 5)**6)

     enddo


!***********************************************************
!                                                           
! Start building layers for BC : imax None None ************
!                                                           
!***********************************************************




!***********************************************************
!                                                           
! BC layer: 1 None None ************************************
!                                                           
!***********************************************************


     do j=idloop(3),idloop(4) 


!***********************************************************
!                                                           
! building source terms in RHS for layer 1 None None d *****
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! d
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 1 None None d ******************
!                                                           
!***********************************************************


qst(nx+2-1,j,indvarsst(1)) =  qst(nx+2-1,j,indvarsst(1))



!***********************************************************
!                                                           
! building source terms in RHS for layer 1 None None eta ***
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! eta
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 1 None None eta ****************
!                                                           
!***********************************************************


qst(nx+2-1,j,indvarsst(2)) =  qst(nx+2-1,j,indvarsst(2))



!***********************************************************
!                                                           
! building source terms in RHS for layer 1 None None ksi ***
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! ksi
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 1 None None ksi ****************
!                                                           
!***********************************************************


qst(nx+2-1,j,indvarsst(3)) =  qst(nx+2-1,j,indvarsst(3))



!***********************************************************
!                                                           
! building source terms in RHS for layer 1 None None stemp *
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (((0.5_wp*(deltayI*([u]_1y)-deltaxI*([v]_1x)))**2+(0.5_wp*(deltaxI*([v]_1x)-deltayI*([u]_1y)))**2)*2)**0.5
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_stemp_dx_0_nxp2m1m1jk = q(nx+2-1-1,j,indvars(3))

d1_stemp_dx_0_nxp2m1p1jk = q(nx+2-1+1,j,indvars(3))

d1_stemp_dx_0_nxp2m1jk = -&
          0.5_wp*d1_stemp_dx_0_nxp2m1m1jk+&
          0.5_wp*d1_stemp_dx_0_nxp2m1p1jk

d1_stemp_dx_0_nxp2m1jk = d1_stemp_dx_0_nxp2m1jk*param_float(1)

d1_stemp_dy_0_nxp2m1jm1k = q(nx+2-1,j-1,indvars(2))

d1_stemp_dy_0_nxp2m1jp1k = q(nx+2-1,j+1,indvars(2))

d1_stemp_dy_0_nxp2m1jk = -&
          0.5_wp*d1_stemp_dy_0_nxp2m1jm1k+&
          0.5_wp*d1_stemp_dy_0_nxp2m1jp1k

d1_stemp_dy_0_nxp2m1jk = d1_stemp_dy_0_nxp2m1jk*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer 1 None None stemp **************
!                                                           
!***********************************************************


qst(nx+2-1,j,indvarsst(4)) =  (((0.5_wp*(qst(nx+2-1,j,indvarsst(11))*(d1_stemp_dy_0_nxp2m1jk)-&
                    qst(nx+2-1,j,indvarsst(10))*(d1_stemp_dx_0_nxp2m1jk)))**2+&
                    (0.5_wp*(qst(nx+2-1,j,indvarsst(10))*(d1_stemp_dx_0_nxp2m1jk)-&
                    qst(nx+2-1,j,indvarsst(11))*(d1_stemp_dy_0_nxp2m1jk)))**2)*2)**0.5



!***********************************************************
!                                                           
! building source terms in RHS for layer 1 None None symm **
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! ((sign(1.0_wp,ksi)-1.0_wp)/(-2.0_wp))
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 1 None None symm ***************
!                                                           
!***********************************************************


qst(nx+2-1,j,indvarsst(5)) =  ((sign(1.0_wp,qst(nx+2-1,j,indvarsst(3)))-&
                    1.0_wp)/(-&
                    2.0_wp))



!***********************************************************
!                                                           
! building source terms in RHS for layer 1 None None detady 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! [eta]_1y
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_detady_dy_0_nxp2m1jm1k = qst(nx+2-1,j-1,indvarsst(2))

d1_detady_dy_0_nxp2m1jp1k = qst(nx+2-1,j+1,indvarsst(2))

d1_detady_dy_0_nxp2m1jk = -&
          0.5_wp*d1_detady_dy_0_nxp2m1jm1k+&
          0.5_wp*d1_detady_dy_0_nxp2m1jp1k

d1_detady_dy_0_nxp2m1jk = d1_detady_dy_0_nxp2m1jk*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer 1 None None detady *************
!                                                           
!***********************************************************


qst(nx+2-1,j,indvarsst(6)) =  d1_detady_dy_0_nxp2m1jk



!***********************************************************
!                                                           
! building source terms in RHS for layer 1 None None dksidy 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! [ksi]_1y
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_dksidy_dy_0_nxp2m1jm1k = qst(nx+2-1,j-1,indvarsst(3))

d1_dksidy_dy_0_nxp2m1jp1k = qst(nx+2-1,j+1,indvarsst(3))

d1_dksidy_dy_0_nxp2m1jk = -&
          0.5_wp*d1_dksidy_dy_0_nxp2m1jm1k+&
          0.5_wp*d1_dksidy_dy_0_nxp2m1jp1k

d1_dksidy_dy_0_nxp2m1jk = d1_dksidy_dy_0_nxp2m1jk*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer 1 None None dksidy *************
!                                                           
!***********************************************************


qst(nx+2-1,j,indvarsst(7)) =  d1_dksidy_dy_0_nxp2m1jk



!***********************************************************
!                                                           
! building source terms in RHS for layer 1 None None detadx 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! [eta]_1x
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_detadx_dx_0_nxp2m1m1jk = qst(nx+2-1-1,j,indvarsst(2))

d1_detadx_dx_0_nxp2m1p1jk = qst(nx+2-1+1,j,indvarsst(2))

d1_detadx_dx_0_nxp2m1jk = -&
          0.5_wp*d1_detadx_dx_0_nxp2m1m1jk+&
          0.5_wp*d1_detadx_dx_0_nxp2m1p1jk

d1_detadx_dx_0_nxp2m1jk = d1_detadx_dx_0_nxp2m1jk*param_float(1)



!***********************************************************
!                                                           
! Update BC terms for layer 1 None None detadx *************
!                                                           
!***********************************************************


qst(nx+2-1,j,indvarsst(8)) =  d1_detadx_dx_0_nxp2m1jk



!***********************************************************
!                                                           
! building source terms in RHS for layer 1 None None dksidx 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! ([ksi]_1x)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_dksidx_dx_0_nxp2m1m1jk = qst(nx+2-1-1,j,indvarsst(3))

d1_dksidx_dx_0_nxp2m1p1jk = qst(nx+2-1+1,j,indvarsst(3))

d1_dksidx_dx_0_nxp2m1jk = -&
          0.5_wp*d1_dksidx_dx_0_nxp2m1m1jk+&
          0.5_wp*d1_dksidx_dx_0_nxp2m1p1jk

d1_dksidx_dx_0_nxp2m1jk = d1_dksidx_dx_0_nxp2m1jk*param_float(1)



!***********************************************************
!                                                           
! Update BC terms for layer 1 None None dksidx *************
!                                                           
!***********************************************************


qst(nx+2-1,j,indvarsst(9)) =  (d1_dksidx_dx_0_nxp2m1jk)



!***********************************************************
!                                                           
! building source terms in RHS for layer 1 None None deltaxI 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! 1.0_wp/(dksidx)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 1 None None deltaxI ************
!                                                           
!***********************************************************


qst(nx+2-1,j,indvarsst(10)) =  1.0_wp/(qst(nx+2-1,j,indvarsst(9)))



!***********************************************************
!                                                           
! building source terms in RHS for layer 1 None None deltayI 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! 1.0_wp/(detady)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 1 None None deltayI ************
!                                                           
!***********************************************************


qst(nx+2-1,j,indvarsst(11)) =  1.0_wp/(qst(nx+2-1,j,indvarsst(6)))



!***********************************************************
!                                                           
! building source terms in RHS for layer 1 None None viscosità 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (1+sut)/(T+sut)*T**1.5
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 1 None None viscosità **********
!                                                           
!***********************************************************


qst(nx+2-1,j,indvarsst(12)) =  (1+&
                    param_float(21 + 5))/(((q(nx+2-1,j,indvars(4))-&
                    0.5_wp*(q(nx+2-1,j,indvars(2))*q(nx+2-1,j,indvars(2))+&
                    q(nx+2-1,j,indvars(3))*q(nx+2-1,j,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+2-1,j,indvars(4))-&
                    0.5_wp*(q(nx+2-1,j,indvars(2))*q(nx+2-1,j,indvars(2))+&
                    q(nx+2-1,j,indvars(3))*q(nx+2-1,j,indvars(3)))))/param_float(4 + 5)**1.5



!***********************************************************
!                                                           
! building source terms in RHS for layer 1 None None fw ****
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! gg*(1+Cw3**6)/(gg**6+Cw3**6)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 1 None None fw *****************
!                                                           
!***********************************************************


qst(nx+2-1,j,indvarsst(13)) =  ((q(nx+2-1,j,indvars(5))/((qst(nx+2-1,j,indvarsst(4))+&
                    param_float(1 + 5)*q(nx+2-1,j,indvars(5))/(param_float(9 + 5)**2*qst(nx+2-1,j,indvarsst(2))**2))*param_float(9 + 5)**2*qst(nx+2-1,j,indvarsst(2))**2))+&
                    param_float(11 + 5)*((q(nx+2-1,j,indvars(5))/((qst(nx+2-1,j,indvarsst(4))+&
                    param_float(1 + 5)*q(nx+2-1,j,indvars(5))/(param_float(9 + 5)**2*qst(nx+2-1,j,indvarsst(2))**2))*param_float(9 + 5)**2*qst(nx+2-1,j,indvarsst(2))**2))**6-&
                    (q(nx+2-1,j,indvars(5))/((qst(nx+2-1,j,indvarsst(4))+&
                    param_float(1 + 5)*q(nx+2-1,j,indvars(5))/(param_float(9 + 5)**2*qst(nx+2-1,j,indvarsst(2))**2))*param_float(9 + 5)**2*qst(nx+2-1,j,indvarsst(2))**2))))*(1+&
                    param_float(12 + 5)**6)/(((q(nx+2-1,j,indvars(5))/((qst(nx+2-1,j,indvarsst(4))+&
                    param_float(1 + 5)*q(nx+2-1,j,indvars(5))/(param_float(9 + 5)**2*qst(nx+2-1,j,indvarsst(2))**2))*param_float(9 + 5)**2*qst(nx+2-1,j,indvarsst(2))**2))+&
                    param_float(11 + 5)*((q(nx+2-1,j,indvars(5))/((qst(nx+2-1,j,indvarsst(4))+&
                    param_float(1 + 5)*q(nx+2-1,j,indvars(5))/(param_float(9 + 5)**2*qst(nx+2-1,j,indvarsst(2))**2))*param_float(9 + 5)**2*qst(nx+2-1,j,indvarsst(2))**2))**6-&
                    (q(nx+2-1,j,indvars(5))/((qst(nx+2-1,j,indvarsst(4))+&
                    param_float(1 + 5)*q(nx+2-1,j,indvars(5))/(param_float(9 + 5)**2*qst(nx+2-1,j,indvarsst(2))**2))*param_float(9 + 5)**2*qst(nx+2-1,j,indvarsst(2))**2))))**6+&
                    param_float(12 + 5)**6)

     enddo
