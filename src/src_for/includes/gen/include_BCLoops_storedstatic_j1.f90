

!***********************************************************
!                                                           
! Start building layers for BC : None j1 None **************
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
! building source terms in RHS for layer None 0 None d *****
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! d
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer None 0 None d ******************
!                                                           
!***********************************************************


qst(i,1-2+0,indvarsst(1)) =  qst(i,1-2+0,indvarsst(1))



!***********************************************************
!                                                           
! building source terms in RHS for layer None 0 None eta ***
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! eta
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer None 0 None eta ****************
!                                                           
!***********************************************************


qst(i,1-2+0,indvarsst(2)) =  qst(i,1-2+0,indvarsst(2))



!***********************************************************
!                                                           
! building source terms in RHS for layer None 0 None ksi ***
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! ksi
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer None 0 None ksi ****************
!                                                           
!***********************************************************


qst(i,1-2+0,indvarsst(3)) =  qst(i,1-2+0,indvarsst(3))



!***********************************************************
!                                                           
! building source terms in RHS for layer None 0 None symm **
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! ((sign(1.0_wp,ksi)-1.0_wp)/(-2.0_wp))
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer None 0 None symm ***************
!                                                           
!***********************************************************


qst(i,1-2+0,indvarsst(5)) =  ((sign(1.0_wp,qst(i,1-2+0,indvarsst(3)))-&
                    1.0_wp)/(-&
                    2.0_wp))



!***********************************************************
!                                                           
! building source terms in RHS for layer None 0 None detady 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! [eta]_1y
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_detady_dy_0_i1m2p0p0k = qst(i,1-2+0+0,indvarsst(2))

d1_detady_dy_0_i1m2p0p1k = qst(i,1-2+0+1,indvarsst(2))

d1_detady_dy_0_i1m2p0p2k = qst(i,1-2+0+2,indvarsst(2))

d1_detady_dy_0_i1m2p0k = -&
          1.5_wp*d1_detady_dy_0_i1m2p0p0k+&
          2.0_wp*d1_detady_dy_0_i1m2p0p1k-&
          0.5_wp*d1_detady_dy_0_i1m2p0p2k

d1_detady_dy_0_i1m2p0k = d1_detady_dy_0_i1m2p0k*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer None 0 None detady *************
!                                                           
!***********************************************************


qst(i,1-2+0,indvarsst(6)) =  d1_detady_dy_0_i1m2p0k



!***********************************************************
!                                                           
! building source terms in RHS for layer None 0 None dksidy 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! [ksi]_1y
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_dksidy_dy_0_i1m2p0p0k = qst(i,1-2+0+0,indvarsst(3))

d1_dksidy_dy_0_i1m2p0p1k = qst(i,1-2+0+1,indvarsst(3))

d1_dksidy_dy_0_i1m2p0p2k = qst(i,1-2+0+2,indvarsst(3))

d1_dksidy_dy_0_i1m2p0k = -&
          1.5_wp*d1_dksidy_dy_0_i1m2p0p0k+&
          2.0_wp*d1_dksidy_dy_0_i1m2p0p1k-&
          0.5_wp*d1_dksidy_dy_0_i1m2p0p2k

d1_dksidy_dy_0_i1m2p0k = d1_dksidy_dy_0_i1m2p0k*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer None 0 None dksidy *************
!                                                           
!***********************************************************


qst(i,1-2+0,indvarsst(7)) =  d1_dksidy_dy_0_i1m2p0k



!***********************************************************
!                                                           
! building source terms in RHS for layer None 0 None detadx 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! [eta]_1x
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_detadx_dx_0_im11m2p0k = qst(i-1,1-2+0,indvarsst(2))

d1_detadx_dx_0_ip11m2p0k = qst(i+1,1-2+0,indvarsst(2))

d1_detadx_dx_0_i1m2p0k = -&
          0.5_wp*d1_detadx_dx_0_im11m2p0k+&
          0.5_wp*d1_detadx_dx_0_ip11m2p0k

d1_detadx_dx_0_i1m2p0k = d1_detadx_dx_0_i1m2p0k*param_float(1)



!***********************************************************
!                                                           
! Update BC terms for layer None 0 None detadx *************
!                                                           
!***********************************************************


qst(i,1-2+0,indvarsst(8)) =  d1_detadx_dx_0_i1m2p0k



!***********************************************************
!                                                           
! building source terms in RHS for layer None 0 None dksidx 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! ([ksi]_1x)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_dksidx_dx_0_im11m2p0k = qst(i-1,1-2+0,indvarsst(3))

d1_dksidx_dx_0_ip11m2p0k = qst(i+1,1-2+0,indvarsst(3))

d1_dksidx_dx_0_i1m2p0k = -&
          0.5_wp*d1_dksidx_dx_0_im11m2p0k+&
          0.5_wp*d1_dksidx_dx_0_ip11m2p0k

d1_dksidx_dx_0_i1m2p0k = d1_dksidx_dx_0_i1m2p0k*param_float(1)



!***********************************************************
!                                                           
! Update BC terms for layer None 0 None dksidx *************
!                                                           
!***********************************************************


qst(i,1-2+0,indvarsst(9)) =  (d1_dksidx_dx_0_i1m2p0k)



!***********************************************************
!                                                           
! building source terms in RHS for layer None 0 None deltaxI 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! 1.0_wp/(dksidx)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer None 0 None deltaxI ************
!                                                           
!***********************************************************


qst(i,1-2+0,indvarsst(10)) =  1.0_wp/(qst(i,1-2+0,indvarsst(9)))



!***********************************************************
!                                                           
! building source terms in RHS for layer None 0 None deltayI 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! 1.0_wp/(detady)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer None 0 None deltayI ************
!                                                           
!***********************************************************


qst(i,1-2+0,indvarsst(11)) =  1.0_wp/(qst(i,1-2+0,indvarsst(6)))



!***********************************************************
!                                                           
! building source terms in RHS for layer None 0 None wall **
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! dabs(1-symm)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer None 0 None wall ***************
!                                                           
!***********************************************************


qst(i,1-2+0,indvarsst(19)) =  dabs(1-&
                    qst(i,1-2+0,indvarsst(5)))

   enddo


!***********************************************************
!                                                           
! Start building layers for BC : None j1 None **************
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
! building source terms in RHS for layer None 1 None d *****
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! d
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer None 1 None d ******************
!                                                           
!***********************************************************


qst(i,1-2+1,indvarsst(1)) =  qst(i,1-2+1,indvarsst(1))



!***********************************************************
!                                                           
! building source terms in RHS for layer None 1 None eta ***
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! eta
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer None 1 None eta ****************
!                                                           
!***********************************************************


qst(i,1-2+1,indvarsst(2)) =  qst(i,1-2+1,indvarsst(2))



!***********************************************************
!                                                           
! building source terms in RHS for layer None 1 None ksi ***
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! ksi
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer None 1 None ksi ****************
!                                                           
!***********************************************************


qst(i,1-2+1,indvarsst(3)) =  qst(i,1-2+1,indvarsst(3))



!***********************************************************
!                                                           
! building source terms in RHS for layer None 1 None symm **
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! ((sign(1.0_wp,ksi)-1.0_wp)/(-2.0_wp))
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer None 1 None symm ***************
!                                                           
!***********************************************************


qst(i,1-2+1,indvarsst(5)) =  ((sign(1.0_wp,qst(i,1-2+1,indvarsst(3)))-&
                    1.0_wp)/(-&
                    2.0_wp))



!***********************************************************
!                                                           
! building source terms in RHS for layer None 1 None detady 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! [eta]_1y
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_detady_dy_0_i1m2p1m1k = qst(i,1-2+1-1,indvarsst(2))

d1_detady_dy_0_i1m2p1p1k = qst(i,1-2+1+1,indvarsst(2))

d1_detady_dy_0_i1m2p1k = -&
          0.5_wp*d1_detady_dy_0_i1m2p1m1k+&
          0.5_wp*d1_detady_dy_0_i1m2p1p1k

d1_detady_dy_0_i1m2p1k = d1_detady_dy_0_i1m2p1k*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer None 1 None detady *************
!                                                           
!***********************************************************


qst(i,1-2+1,indvarsst(6)) =  d1_detady_dy_0_i1m2p1k



!***********************************************************
!                                                           
! building source terms in RHS for layer None 1 None dksidy 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! [ksi]_1y
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_dksidy_dy_0_i1m2p1m1k = qst(i,1-2+1-1,indvarsst(3))

d1_dksidy_dy_0_i1m2p1p1k = qst(i,1-2+1+1,indvarsst(3))

d1_dksidy_dy_0_i1m2p1k = -&
          0.5_wp*d1_dksidy_dy_0_i1m2p1m1k+&
          0.5_wp*d1_dksidy_dy_0_i1m2p1p1k

d1_dksidy_dy_0_i1m2p1k = d1_dksidy_dy_0_i1m2p1k*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer None 1 None dksidy *************
!                                                           
!***********************************************************


qst(i,1-2+1,indvarsst(7)) =  d1_dksidy_dy_0_i1m2p1k



!***********************************************************
!                                                           
! building source terms in RHS for layer None 1 None detadx 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! [eta]_1x
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_detadx_dx_0_im11m2p1k = qst(i-1,1-2+1,indvarsst(2))

d1_detadx_dx_0_ip11m2p1k = qst(i+1,1-2+1,indvarsst(2))

d1_detadx_dx_0_i1m2p1k = -&
          0.5_wp*d1_detadx_dx_0_im11m2p1k+&
          0.5_wp*d1_detadx_dx_0_ip11m2p1k

d1_detadx_dx_0_i1m2p1k = d1_detadx_dx_0_i1m2p1k*param_float(1)



!***********************************************************
!                                                           
! Update BC terms for layer None 1 None detadx *************
!                                                           
!***********************************************************


qst(i,1-2+1,indvarsst(8)) =  d1_detadx_dx_0_i1m2p1k



!***********************************************************
!                                                           
! building source terms in RHS for layer None 1 None dksidx 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! ([ksi]_1x)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_dksidx_dx_0_im11m2p1k = qst(i-1,1-2+1,indvarsst(3))

d1_dksidx_dx_0_ip11m2p1k = qst(i+1,1-2+1,indvarsst(3))

d1_dksidx_dx_0_i1m2p1k = -&
          0.5_wp*d1_dksidx_dx_0_im11m2p1k+&
          0.5_wp*d1_dksidx_dx_0_ip11m2p1k

d1_dksidx_dx_0_i1m2p1k = d1_dksidx_dx_0_i1m2p1k*param_float(1)



!***********************************************************
!                                                           
! Update BC terms for layer None 1 None dksidx *************
!                                                           
!***********************************************************


qst(i,1-2+1,indvarsst(9)) =  (d1_dksidx_dx_0_i1m2p1k)



!***********************************************************
!                                                           
! building source terms in RHS for layer None 1 None deltaxI 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! 1.0_wp/(dksidx)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer None 1 None deltaxI ************
!                                                           
!***********************************************************


qst(i,1-2+1,indvarsst(10)) =  1.0_wp/(qst(i,1-2+1,indvarsst(9)))



!***********************************************************
!                                                           
! building source terms in RHS for layer None 1 None deltayI 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! 1.0_wp/(detady)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer None 1 None deltayI ************
!                                                           
!***********************************************************


qst(i,1-2+1,indvarsst(11)) =  1.0_wp/(qst(i,1-2+1,indvarsst(6)))



!***********************************************************
!                                                           
! building source terms in RHS for layer None 1 None wall **
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! dabs(1-symm)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer None 1 None wall ***************
!                                                           
!***********************************************************


qst(i,1-2+1,indvarsst(19)) =  dabs(1-&
                    qst(i,1-2+1,indvarsst(5)))

   enddo
