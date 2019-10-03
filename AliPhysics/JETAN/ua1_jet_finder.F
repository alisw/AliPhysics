c
c => Original copy from /star/emc/shester/cl/pams/emc/jet/jet_finer_ua1.F
c
c:>------------------------------------------------------------------
C:ROUTINE:      subroutine jet_finder_ua1
C:DESCRIPTION:  UA1 jet algorithm from LUND JETSET
C:RETURN VALUE: ierror=1 on error
c:>------------------------------------------------------------------
      subroutine ua1_jet_finder(
     +        ncell, ncell_tot, etc, etac, phic,
     +        min_move, max_move, mode, prec_bg, ierror)
      implicit none ! 5-oct-2001

      real C_PI
      real C_2PI 
      real etaCellSize 
      real phiCellSize 
      real arg
      INTEGER NMAX,JMAX
      parameter(NMAX=60000,JMAX=100) ! 10-oct-2201
      integer ncell, ierror, mode, ncell_tot
      real etc(NMAX),etac(NMAX),phic(NMAX)
      real cone_rad, et_seed, ej_min, et_min
      real min_move, max_move, prec_bg
      integer flag(NMAX),index(NMAX)
      integer n_iter,i,j,k,l,nc
      real et_sum,et_ave,et_sum_old
      real et,eta,phi,eta0,phi0,etab,phib,ets,etas,phis
      real deta,dphi,r_int
!
      real occupationAll, occupationInJet  ! 3-oct-2001 by PAI 
      integer maxTowerInJet                ! 5-oct-2001 by PAI
      save    maxTowerInJet                ! This is obligatory
      integer idPerfomance
      data    idPerfomance /119/
      save    idPerfomance
!  print parameter - 27-sep-2001 by PAI
      integer kpri
      data    kpri /0/
      integer njet, ncellj
      real    etj, etaj, phij, etavg
*    Results
      COMMON /UA1JETS/ NJET, ETJ(100), ETAJ(100,2), PHIJ(100,2), 
     +     NCELLJ(100), ETAVG
*    Cell Geometry
      COMMON /UA1CELL/ etaCellSize, phiCellSize
*    Parameters
      COMMON /UA1PARA/ cone_rad, et_seed, ej_min, et_min

      C_PI  = 3.1415926
      C_2PI = 2.*C_PI

      if(kpri.eq.0) then
        kpri = 1
! for add. correction of jet energy if ocupation in jet not 100% !!
! may differ from real number because grid 
        maxTowerInJet=
     +       int((C_PI*cone_rad*cone_rad)/(etaCellSize*phiCellSize)+0.5)
        print 1, ncell_tot
     +       ,cone_rad, et_seed, ej_min, et_min
     +       ,min_move, max_move, mode, prec_bg
     +       ,maxTowerInJet
 1    format(/
     +    '    == jet_finder_UA1 ==  '/
     +    ' ncell_tot                   ', i5/
     +    ' cone rad                    ', f5.2/
     +    ' et_seed                     ', f5.2,' GeV/C'/
     +    ' ej_min                      ', f5.2,' GeV/C'/
     +    ' et_min(tower after bg sub.) ', f5.2,' GeV/C'/
     +    ' min_cone_move               ', f5.3/
     +    ' max_cone_move               ', f5.3/
     +    ' Mode for BG subtraction     ', i5/
     +    ' prec_bg                     ', f5.4/
     +    ' -- PAI"s addition --'/
     +    ' maxTowerInJet               ', i5/
     +    ' -- '/
     +    )
        if(NMAX .lt. ncell_tot) then
          print 2, NMAX, ncell_tot
 2    format('<E> you must increase memory -> NMAX ',i6
     +      ,' ncell_tot ', i6)
        endif 
      endif
      occupationAll = float(ncell) / float(ncell_tot)
! print parameter - 27-sep-2001 by PAI
      ierror =0
      n_iter =0
c*-sort cells by Et decending order, store the order in index
      call sortzv(etc,index,ncell,-1,1,0)   
c*-sum up total energy of all cells
      if(mode .eq. 1) then
         n_iter=1
         et_sum=0.0
         do i=1, ncell
            et_sum=et_sum+etc(i)
         enddo
         et_sum_old=et_sum
         et_ave=et_sum/float(ncell_tot)
      else
         et_ave=0.0
      endif
      print *,'Iter ', n_iter, ' et_ave ', et_ave, ' #cells ', ncell 
c*-Watch out!!! For mode=1, it can jump back to here for next iteration!!!
 999  continue
c*-kill cells (flag=2) with Et below ET_MIN after background subtraction
cfca      call vzero(flag,ncell)
      do i=1, ncell
         flag(i)=0
         if(etc(i)-et_ave .le. et_min) flag(i)=2
      enddo
      njet = 0
c*-Initiator cell is the one with highest Et of not yet used ones
      i=1
      j=index(i)
      if(i.eq.1. and. etc(j).lt.et_seed) then
        if(mode.eq.0 .or. (mode.eq.1 .and. n_iter.eq.1)) then
          print *,' no cell with Et higher than et_seed ->', etc(j)       
          return
        endif
      endif 
      do while(etc(j) .ge. et_seed)
         if(flag(j) .eq. 0) then
C         write(6,*) j, etc(j), et_seed, etac(j), phic(j), flag(j) 
         
            et =etc(j)-et_ave
            eta=etac(j)
            phi=phic(j)
            eta0=eta
            phi0=phi
            etab=eta
            phib=phi
            ets =0.0
            etas=0.0
            phis=0.0
c*-weighted eta and phi. 
            do k = 1, ncell
               l = index(k)
               if(flag(l).eq.0) then
                  deta = etac(l)-eta
c*-Is the cell is in the cone?
                  if(abs(deta).le.cone_rad)then
                     dphi=phic(l)-phi
                     do while(dphi .gt. C_PI)
                        dphi=dphi-C_2PI
                     enddo
                     do while(dphi .le. -C_PI)
                        dphi=dphi+C_2PI
                     enddo
                     if(abs(dphi).le.cone_rad) then
                        r_int=sqrt(deta**2+dphi**2)
                        if(r_int.le.cone_rad)then
c*-calculate offset from initiate cell
                           deta=etac(l)-eta0
                           dphi=phic(l)-phi0
                           do while(dphi .gt. C_PI)
                              dphi=dphi-C_2PI
                           enddo
                           do while(dphi .lt. -C_PI)
                              dphi=dphi+C_2PI
                           enddo
                           et=etc(l)-et_ave
                           etas=etas+abs(et)*deta
                           phis=phis+abs(et)*dphi
                           ets =ets +et
c*-New weighted eta and phi including this cell
                           eta=eta0+etas/ets
                           phi=phi0+phis/ets                          
c*-If cone does not move much from previous cone, just go next step
                           r_int=sqrt((eta-etab)**2+(phi-phib)**2)
                           if(r_int .le. min_move) then
                              goto 159
                           endif
c*-Cone should not move more than MAX_CONE_MOVE from initiator cell
                           r_int=sqrt((etas/ets)**2+(phis/ets)**2)              
                           if(r_int .ge. max_move) then
                              eta=etab
                              phi=phib
                              goto 159
                           endif
c*-Store this loop information
                           etab=eta
                           phib=phi
                        endif
                     endif
                  endif
               endif
            enddo
 159        continue 
            
c*-sum up unused cells within required distance of given eta/phi
            nc=0
            ets=0.0
            etas=0.0
            phis=0.0
            do k=1,ncell
               l=index(k)
               if(flag(l) .eq. 0) then
                  deta=etac(l)-eta
                  if(abs(deta).le.cone_rad)then
                     dphi=phic(l)-phi
                     do while(dphi .gt. C_PI)
                        dphi=dphi-C_2PI
                     enddo
                     do while(dphi .le. -C_PI)
                        dphi=dphi+C_2PI
                     enddo
                     if(abs(dphi).le.cone_rad) then
                        r_int=sqrt(deta**2+dphi**2)
                        if(r_int.le.cone_rad)then
                           flag(l)=-1
                           et  =etc(l)-et_ave
                           ets =ets +et
                           etas=etas+et*deta
                           phis=phis+et*dphi
                           nc  = nc + 1
                        endif
                     endif
                  endif
               endif
            enddo  ! do k=1,ncell
!  5-oct-2001 by PAI - remove 20-feb-2002 by PAI
! 20-feb-2002 - it is work if you apply cut on eT before jet finder !!!
!            if(maxTowerInJet .gt. nc) then
!              ets = ets - et_ave*(maxTowerInJet - nc) 
!            endif
! 5-oct-2001 by PAI
            
c*-reject cluster below minimum Ej_min
c* protection (am)

c            arg = 0.
c            if (ets .ne. 0.) then
c               if (abs(etas/ets) .lt. 23.719) then
c                  arg = ets * cosh(etas/ets)
c               else
c                  arg = 1.e10
c               endif
c            endif
            
            if(ets .lt. ej_min) then
               do k=1,ncell
                  if(flag(k).le.0) flag(k)=0
               enddo
            else
c*-eles, store flags and jet variables
               do k=1,ncell
                  if(flag(k).eq.-1) flag(k)=1
               enddo
               etas=eta+etas/ets
               phi=phi+phis/ets
               do while(phi .ge. C_2PI)
                  phi=phi-C_2PI
               enddo
               do while(phi .lt. 0.0)
                  phi=phi+C_2PI
               enddo
               njet=njet+1
               etj(njet) =ets
               etaj(njet,1)=eta0
               phij(njet,1)=phi0
               etaj(njet,2)=etas
               phij(njet,2)=phi
               ncellj(njet)=nc
               etavg = et_ave
            endif 
         endif
         i=i+1
         j=index(i)        
      enddo
      
c*-recalculate energy sum excluding used cells.
      if(mode.eq.1)then
         et_sum=0.0
         nc=0       ! #cells in jets
         do i=1,ncell
            if(flag(i).ne.1) then  ! 1 if cell in jet
               et_sum=et_sum+etc(i)
            else
               nc=nc+1
            endif
         enddo
c*-if background level changes more than prec_bg, go next iteration!!!
c*-after 10 iteration, stop working and finish
         if( (et_sum .gt. 0.) 
     +        .and. (abs(et_sum-et_sum_old)/et_sum.gt.prec_bg 
     +        .and. n_iter.le.10)
     +        .or. n_iter.eq.1) then ! minimum 2 iteration - 10-oct-2001 by pai
            et_ave=et_sum/float(ncell_tot-nc)
            n_iter=n_iter+1
            et_sum_old = et_sum
            print *,'End of iteration : et_ave ', et_ave, ' nc ', nc
     + , ' Jet energy ', ets     
            goto 999
c*-Watch out!!! Here is a big jump!!! 
         endif
         occupationInJet = float(ncell) / float(ncell_tot)
         write(*,*) njet,' jet found in ',n_iter,
     +     ' iteration(s) : EtSum, EtAve =',
     +     et_sum,et_sum/float(ncell_tot-nc)
     + , ' ncell_tot ', ncell_tot, ' #cell in jet ', nc
     + , ' occupationAll ', occupationAll
     + , ' occupationInJet ', occupationInJet
      endif
      
      if(njet.gt.100)then
         write(*,*)'UA1:Problem:More than 100 jets found!'
         ierror = 1
      elseif(njet.eq.0)then
         write(*,*)'UA1:Done:No jet found!'
      else
         write(*,*)
     +'UA1:Done:Found ',njet,' jet(s)'         
      end if
      return
      end



