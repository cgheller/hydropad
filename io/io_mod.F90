MODULE io_mod
!
	IMPLICIT NONE
	SAVE
!
        integer :: stepbox,initboxyn
        integer :: initcheckyn

        integer :: initproject,nprj
        integer :: nchk,npchk
	real*8  :: frome,toe
	
	integer, parameter :: nn=10	
	integer :: nlb(nn),nxb(nn),nyb(nn),nzb(nn)
	integer :: resident(nn,200)
	integer :: kbox(nn,200)
	integer :: nbox
	real*8  :: lb(nn),xb(nn),yb(nn),zb(nn)
	real*8  :: tbox
!
	character*60 :: cosmological_model,owner
	integer :: date_of_production,s_id
	
!
END MODULE io_mod
