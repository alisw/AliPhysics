! LHACONTROL common block
!implicit none
character*20 lhaparm(20)
double precision lhavalue(20)
character*20 commoninitflag
common /lhacontrol/ lhaparm, lhavalue, commoninitflag
save /lhacontrol/
