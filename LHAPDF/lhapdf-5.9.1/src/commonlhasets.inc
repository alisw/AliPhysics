! LHASETS common block
!implicit none
include 'parmsetup.inc'
character*1024 lhanames(nmxset)
integer lhanumbers(nmxset), lhamembers(nmxset)
integer nsets
common /lhasets/ lhanames, lhanumbers, lhamembers, nsets
save /lhasets/
