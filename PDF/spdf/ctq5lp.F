#include "pdf/pilot.h"
C   CTEQ5M1 and CTEQ5L Parton Distribution Functions in Parametrized Form
C                             
C               September 15, 1999
C
C   Ref: "GLOBAL QCD ANALYSIS OF PARTON STRUCTURE OF THE NUCLEON:
C         CTEQ5 PPARTON DISTRIBUTIONS"
C   hep-ph/9903282
C
C   The CTEQ5M1 set given here is an updated version of the original CTEQ5M
C     set posted, in the table version, on the Web page of CTEQ.
C     The differences between CTEQ5M and CTEQ5M1 are insignificant for almost
C     all applications. 
C   The improvement is in the QCD evolution which is now more accurate, and
C   which agrees completely with the benchmark work of the HERA 96/97 Workshop.
C
C   The differences between the parametrized and the corresponding table ver-
C sions (on which it is based) are of similar order as between the two version.
C    
C!! Because accurate parametrizations over a wide range of (x,Q) is hard to
C   obtain, only the most widely used sets CTEQ5M and CTEQ5L are available 
C   in parametrized form for now. 
C
C   These parametrizations were obtained by Jon Pumplin.
C
C                    ******************************
C  Iset   PDF        Description                 Alpha_s(Mz)  Lam4  Lam5
C ---------------------------------------------------------------------------
C   1    CTEQ5M1  Standard NLO MSbar scheme         0.118     326   226
C   3    CTEQ5L   Leading Order                     0.127     192   146
C ---------------------------------------------------------------------------
C   Note the Qcd-lambda values given for CTEQ5L is for the leading order
C     form of Alpha_s!!  Alpha_s(Mz) gives the absolute calibration.
C
C  The two Iset value are adopted to agree with the standard table versions.
C
C   The following user-callable routines are provided:
C 
C     FUNCTION Ctq5Pd (Iset, Iprtn, X, Q, Irt) 
C         returns the PROBABILITY density for a GIVEN flavor;
C
C     FUNCTION Ctq5df (Iset, Iprtn, X, Q, Irt)
C         returns the MOMENTUM density of a GIVEN valence or sea distribution.
C
C     SUBROUTINE Ctq5Pds(Iset, Pdf, X, Q, Irt)
C         returns an array of MOMENTUM densities for ALL flavors;
C
C   The arguments of these routines are as follows: 
C
C   Iset is the set number:  1 for CTEQ5M1 or 3 for CTEQ5L  
C
C   Iprtn  is the parton label (6, 5, 4, 3, 2, 1, 0, -1, ......, -6)
C                          for (t, b, c, s, d, u, g, u_bar, ..., t_bar)
C  *** WARNING: We use the parton label 2 as D-quark and 1 as U-quark, 
C               which might be different from your labels.
C
C   X, Q are the usual x, Q; 
C
C   Irt is an error code: 0 if there was no error; 1 or more if (x,q) was 
C   outside the range of validity of the parametrization.
C       
C  Range of validity:
C  
C     The range of (x, Q) covered by this parametrization of the QCD evolved
C     parton distributions is 1E-6 < x < 1 ; 1.1 GeV < Q < 10 TeV.  Of course,
C     the PDF's are constrained by data only in a subset of that region; and 
C     the assumed DGLAP evolution is unlikely to be valid for all of it either.
C
C     The range of (x, Q) used in the CTEQ5 round of global analysis is 
C     approximately 0.01 < x < 0.75 ; and 4 GeV^2 < Q^2 < 400 GeV^2 for 
C     fixed target experiments; 0.0001 < x < 0.3 from HERA data; and   
C     Q^2 up to 40,000 GeV^2 from Tevatron inclusive Jet data.
C
C   DOUBLE PRECISION is used throughout in these routines, but conversion to 
C   SINGLE PRECISION is possible by 
C                    removing the Implicit Double Precision statements. 
C
C **************************************************************************
C
C ********************************************************
c --------------------------------------------------------------------------
	double precision function ctq5Lp(ifl,x,q)
c Parametrization of cteq5L parton distribution functions (J. Pumplin 9/99).
c ifl: 1=u,2=d,3=s,4=c,5=b;0=gluon;-1=ubar,-2=dbar,-3=sbar,-4=cbar,-5=bbar.
c --------------------------------------------------------------------------
#include "pdf/impdp.inc"
	integer ifl

	ii = ifl
	if(ii .gt. 2) then
	   ii = -ii
	endif

	if(ii .eq. -1) then
	   sum = faux5L(-1,x,q)
	   ratio = faux5L(-2,x,q)
	   ctq5Lp = sum/(1.d0 + ratio)

	elseif(ii .eq. -2) then
	   sum = faux5L(-1,x,q)
	   ratio = faux5L(-2,x,q)
	   ctq5Lp = sum*ratio/(1.d0 + ratio)

	elseif(ii .ge. -5) then
	   ctq5Lp = faux5L(ii,x,q)

	else
	   ctq5Lp = 0.d0 

	endif

	return
	end
