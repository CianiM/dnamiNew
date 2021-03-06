

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
! building source terms in RHS for layer None 0 None M_jmax 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (u/c)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer None 0 None M_jmax *************
!                                                           
!***********************************************************


qface_j(i,1) =  (q(i,ny+2+0,indvars(2))/(param_float(23 + 5)*(param_float(3 + 5))*q(i,ny+2+0,indvars(1))*((q(i,ny+2+0,indvars(4))-&
                    0.5_wp*(q(i,ny+2+0,indvars(2))*q(i,ny+2+0,indvars(2))+&
                    q(i,ny+2+0,indvars(3))*q(i,ny+2+0,indvars(3)))))/q(i,ny+2+0,indvars(1))))

   enddo
