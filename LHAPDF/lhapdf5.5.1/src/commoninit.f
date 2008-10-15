! -*- F90 -*-


subroutine commoninit()
  include 'common.inc'
  integer i

  if (commoninitflag .ne. 'commonblockinitdone') then
     !print *, "Initialising LHAPDF steering data"
     commoninitflag = 'commonblockinitdone'

     ! LHAPDF common block
     lhaname = ' '
     lhaset = 0
     lhamemb = 0
     
     ! LHASETS common block
     do i = 1, nmxset
        lhanames(i) = ' '
        lhanumbers(i) = 0
        lhamembers(i) = 0
     end do
     nsets = 0
     
     ! LHAPDFC common block
     lhapath = 'pdfsets'
     
     ! LHACONTROL common block

     do i = 1, 20
        lhaparm(i) = ' '
        lhavalue(i) = 0.0d0
     end do

     ! LHAGLSTA common block
     xminnum = 0.0d0
     xmaxnum = 0.0d0
     q2minnum = 0.0d0
     q2maxnum = 0.0d0
     totnum = 0.0d0
     xminnup = 0.0d0
     xmaxnup = 0.0d0
     q2minnup = 0.0d0
     q2maxnup = 0.0d0
     totnup = 0.0d0
  end if
  
end subroutine commoninit
