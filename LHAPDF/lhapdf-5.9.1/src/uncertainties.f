! -*- F90 -*-

!--   31/05/2012 by Graeme Watt <Graeme.Watt(at)cern.ch>.
!--   Calculation of PDF uncertainties and PDF correlations.
!--
!--   Use formulae for PDF uncertainties and correlations in:
!--    G. Watt, JHEP 1109 (2011) 069 [arXiv:1106.5788 [hep-ph]].
!--   Code should distinguish between NNPDF (Monte Carlo approach),
!--   Alekhin02/ABKM09/ABM11 (symmetric Hessian approach).  Other
!--   PDF sets are assumed to use the asymmetric Hessian approach.
!--   List of subroutines in this file:
!--    GetMaxNumSets(MaxNumSets)
!--    GetPDFUncType(lMonteCarlo,lSymmetric)
!--    GetPDFUncTypeM(nset,lMonteCarlo,lSymmetric)
!--    GetPDFuncertainty(values,central,errplus,errminus,errsym)
!--    GetPDFuncertaintyM(nset,values,central,errplus,errminus,errsym)
!--    GetPDFcorrelation(valuesA,valuesB,correlation)
!--    GetPDFcorrelationM(nset,valuesA,valuesB,correlation)
!--   Input variables above are nset, values, valuesA, valuesB.
!--   Other arguments above are all output variables.


!-- Get flags indicating if Monte Carlo PDF set (NNPDF) and
!-- if should compute symmetric errors (NNPDF, Alekhin).

subroutine GetPDFUncType(lMonteCarlo,lSymmetric)
  implicit none
  logical lMonteCarlo,lSymmetric
  integer nset
  nset = 1
  call GetPDFUncTypeM(nset,lMonteCarlo,lSymmetric)
  return
end subroutine GetPDFUncType

subroutine GetPDFUncTypeM(nset,lMonteCarlo,lSymmetric)
  implicit none
  include 'parmsetup.inc'
  logical lMonteCarlo,lSymmetric
  character*16 name(nmxset)
  integer nset,nmem(nmxset),ndef(nmxset),mem
  common/NAME/name,nmem,ndef,mem
  if ((name(nset).eq.'NNPDF').or.(name(nset).eq.'NNPDFint').or. &
       (name(nset).eq.'NNPDF20int').or.(name(nset).eq.'NNPDF20intqed')) then ! NNPDF Monte Carlo PDF sets
     lMonteCarlo = .true.
     lSymmetric = .true.
  else if ((name(nset).eq.'A02M').or.(name(nset).eq.'ABKM09').or. &
       (name(nset).eq.'ABM11')) then ! symmetric eigenvector PDF sets
     lMonteCarlo = .false.
     lSymmetric = .true.
  else ! default: assume asymmetric Hessian eigenvector PDF sets
     lMonteCarlo = .false.
     lSymmetric = .false.
  end if
end subroutine GetPDFUncTypeM


!-- Calculate the PDF uncertainty using the appropriate formula for
!-- either the Hessian or Monte Carlo approach given an array
!-- "values(0:nmem)".  In the Monte Carlo approach, the uncertainty is
!-- given by the standard deviation, and the central (average) value
!-- is not necessarily "values(0)" for quantities with a non-linear
!-- dependence on PDFs.  In the Hessian approach, the central value is
!-- the best-fit "values(0)" and the uncertainty is given by either
!-- the symmetric or asymmetric formula using eigenvector PDF sets.

subroutine GetPDFuncertainty(values, &
     &     central,errplus,errminus,errsym)
  implicit none
  integer nset
  double precision values(0:*),central,errplus,errminus,errsym
  nset = 1
  call GetPDFuncertaintyM(nset,values, &
       &     central,errplus,errminus,errsym)
  return
end subroutine GetPDFuncertainty


subroutine GetPDFuncertaintyM(nset,values, &
     &     central,errplus,errminus,errsym)
  implicit none
  integer nset,nmem,imem
  double precision values(0:*),central,errplus,errminus,errsym
  logical lMonteCarlo,lSymmetric

  call numberPDFM(nset,nmem)
  call GetPDFUncTypeM(nset,lMonteCarlo,lSymmetric)

  central = 0.D0 ! central value
  errplus = 0.D0 ! positive uncertainty
  errminus = 0.D0 ! negative uncertainty
  errsym = 0.D0 ! symmetrised uncertainty

  if (lMonteCarlo) then ! calculate average and standard deviation

     do imem = 1, nmem
        central = central + values(imem)
        errsym = errsym + values(imem)**2
     end do
     central = central/nmem ! mean of values
     errsym = errsym/nmem ! mean of squared values
     errsym = nmem/(nmem-1.D0)*(errsym-central**2)
     if (errsym.gt.0.D0) then
        errsym = sqrt(errsym)
     else
        errsym = 0.D0
     end if
     errplus = errsym
     errminus = errsym

  else if (lSymmetric) then ! symmetric Hessian eigenvector PDF sets

     do imem = 1, nmem
        errsym = errsym + (values(imem)-values(0))**2
     end do
     errsym = sqrt(errsym)
     errplus = errsym
     errminus = errsym
     central = values(0)

  else ! default: assume asymmetric Hessian eigenvector PDF sets

     !-- check that nmem is non-zero and even
     if (nmem.ne.0.and.(mod(nmem,2).eq.0)) then
        do imem = 1, nmem/2 ! sum over eigenvectors
           errplus = errplus + max(0.D0, &
                &        values(2*imem-1)-values(0), &
                &        values(2*imem)-values(0))**2
           errminus = errminus + max(0.D0, &
                &        values(0)-values(2*imem-1), &
                &        values(0)-values(2*imem))**2
           errsym = errsym + (values(2*imem-1)-values(2*imem))**2
        end do
        errplus = sqrt(errplus)
        errminus = sqrt(errminus)
        errsym = 0.5D0*sqrt(errsym)
     end if
     central = values(0)

  end if

  return
end subroutine GetPDFuncertaintyM


!-- Calculate the PDF correlation using the appropriate formula for
!-- either the Hessian or Monte Carlo approach given two arrays
!-- "valuesA(0:nmem)" and "valuesB(0:nmem)".  The correlation can vary
!-- between -1 and +1 where values close to {-1,0,+1} mean that the two
!-- quantities A and B are {anticorrelated,uncorrelated,correlated}.

subroutine GetPDFcorrelation(valuesA,valuesB,correlation)
  implicit none
  integer nset
  double precision valuesA(0:*),valuesB(0:*),correlation
  nset = 1
  call GetPDFcorrelationM(nset,valuesA,valuesB,correlation)
  return
end subroutine GetPDFcorrelation

subroutine GetPDFcorrelationM(nset,valuesA,valuesB,correlation)
  implicit none
  integer nset,nmem,imem
  double precision valuesA(0:*),valuesB(0:*),correlation
  double precision A0,Ap,Am,As,B0,Bp,Bm,Bs
  logical lMonteCarlo,lSymmetric

  call numberPDFM(nset,nmem)
  call GetPDFUncTypeM(nset,lMonteCarlo,lSymmetric)

  call GetPDFuncertaintyM(nset,valuesA,A0,Ap,Am,As)
  call GetPDFuncertaintyM(nset,valuesB,B0,Bp,Bm,Bs)

  correlation = 0.D0

  if (lMonteCarlo) then ! calculate average and standard deviation

     do imem = 1, nmem
        correlation = correlation + valuesA(imem)*valuesB(imem)
     end do
     correlation = (correlation/nmem - A0*B0)/(As*Bs)*nmem/(nmem-1.D0)

  else if (lSymmetric) then ! symmetric Hessian eigenvector PDF sets

     do imem = 1, nmem
        correlation = correlation + &
             & (valuesA(imem)-A0)*(valuesB(imem)-B0)
     end do
     correlation = correlation/(As*Bs)

  else ! default: assume asymmetric Hessian eigenvector PDF sets

     !-- check that nmem is non-zero and even
     if (nmem.ne.0.and.(mod(nmem,2).eq.0)) then
        do imem = 1, nmem/2 ! sum over eigenvectors
           correlation = correlation + &
                & (valuesA(2*imem-1)-valuesA(2*imem)) * &
                & (valuesB(2*imem-1)-valuesB(2*imem))
        end do
        correlation = correlation/(4.D0*As*Bs)
     end if

  end if

  return
end subroutine GetPDFcorrelationM
