! LHAPDF common block
!implicit none
character*1024 lhaname
integer lhaset, lhamemb
common /lhapdf/ lhaname, lhaset, lhamemb
save /lhapdf/
