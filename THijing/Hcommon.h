
#ifndef ROOT_HCommon
#define ROOT_HCommon

#ifndef __CFORTRAN_LOADED
//*KEEP,cfortran.
#include "cfortran.h"
//*KEND.
#endif

extern "C" {

/*=========================================================*/
/* COMMON/HIPARNT/HIPR1(100),IHPR2(50),HINT1(100),IHNT2(50)*/
/*---------------------------------------------------------*/
typedef struct {
   Float_t    hipr1[100];
   Int_t      ihpr2[50];
   Float_t    hint1[100];
   Int_t      ihnt2[50];
} HiparntCommon;

#define HIPARNT COMMON_BLOCK(HIPARNT,hiparnt)
COMMON_BLOCK_DEF(HiparntCommon,HIPARNT);

/**************************************************************************/
/*           D E S C R I P T I O N :                                      */
/*------------------------------------------------------------------------*/
/*COMMON/HIPARNT/HIPR1(100), IHPR2(50), HINT1(100), IHNT2(50)             */
/*Purpose: contains input parameters (HIPR1, IHPR2) for event options	  */
/*	   and some extra information (HINT1, IHNT2) of current event.	  */
/*HIPR1(1): (D=1.5 GeV/$c^2$) minimum value for the invariant mass of     */
/*	   the excited string system in a hadron-hadron interaction.	  */
/*HIPR1(2): (D=0.35 GeV) width of the Gaussian $P_T$ distribution of	  */
/*	   produced hadron in Lund string fragmentation 		  */
/*	   (PARJ(21) in JETSET 7.2).					  */
/*HIPR1(3), HIPR1(4): (D=0.5, 0.9 GeV$^{-2}$) give the $a$ and $b$	  */
/*	   parameters of the symmetric Lund fragmentation function	  */
/*	   (PARJ(41), PARJ(42) in JETSET 7.2).				  */
/*HIPR1(5): (D=2.0 GeV/$c^2$) invariant mass cut-off for the dipole	  */
/*	   radiation of a string system below which soft gluon		  */
/*	   radiations are terminated.					  */
/*HIPR1(6): (D=0.1) the depth of shadowing of structure functions	  */
/*	   at $x=0$:$\alpha_A=\mbox{HIPR1(6)}\times(A^{1/3}-1)$.	  */
/*HIPR1(7): not used							  */
/*HIPR1(8): (D=2.0 GeV/$c$) minimum $P_T$ transfer in hard or		  */
/*	   semihard scatterings.					  */
/*HIPR1(9): (D=$-1.0$ GeV/$c$) maximum $P_T$ transfer in hard or	  */
/*	   semihard scatterings. If negative, the limit is set		  */
/*	   by the colliding energy.					  */
/*HIPR1(10): (D=$-2.25$ GeV/$c$) specifies the value of $P_T$ for	  */
/*	   each triggered hard scattering generated per event		  */
/*	   (see Section \ref{sec:jet2}). If HIPR1(10) is negative,	  */
/*	   its absolute value gives the low limit of the		  */
/*	   $P_T$ of the triggered jets. 				  */
/*HIPR1(11): (D=2.0 GeV/$c$) minimum $P_T$ of a jet which will interact   */
/*	   with excited nuclear matter. When the $P_T$ of a jet 	  */
/*	   is smaller than HIPR1(11) it will stop interacting further.	  */
/*HIPR1(12): (D=1.0 fm) transverse distance between a traversing jet	  */
/*	   and an excited nucleon (string system) below which they	  */
/*	   will interact and the jet will lose energy and momentum	  */
/*	   to that string system.					  */
/*HIPR1(13): (D=1.0 fm) the mean free path of a jet when it goes	  */
/*	   through the excited nuclear matter.				  */
/*HIPR1(14): (D=2.0 GeV/fm) the energy loss $dE/dz$ of a gluon		  */
/*	   jet inside the excited nuclear matter. The energy loss	  */
/*	   for a quark jet is half of the energy loss of a gluon.	  */
/*HIPR1(15): (D=0.2 GeV/$c$) the scale $\Lambda$ in the 		  */
/*	   calculation of $\alpha_s$.					  */
/*HIPR1(16): (D=2.0 GeV/$c$) the initial scale $Q_0$ for the		  */
/*	   evolution of the structure functions.			  */
/*HIPR1(17): (D=2.0) $K$ factor for the differential jet cross		  */
/*	   sections in the lowest order pQCD calculation.		  */
/*HIPR1(18): not used							  */
/*HIPR1(19), HIPR1(20): (D=0.1, 1.4 GeV/$c$) parameters in the		  */
/*	   distribution for the $P_T$ kick from soft interactions ,	  */
/*	   $1/[(\mbox{HIPR1(19)}^2+P_T^2)(\mbox{HIPR1(20)}^2+P_T^2)]$.	  */
/*HIPR1(21): (D=1.6 GeV/$c$) the maximum $P_T$ for soft interactions,	  */
/*	   beyond which a Gaussian distribution as specified by 	  */
/*	   HIPR1(2) will be used.					  */
/*HIPR1(22): (D=2.0 GeV/$c$) the scale in the form factor to suppress	  */
/*	   the $P_T$ transfer to diquarks in hard scatterings,		  */
/*HIPR1(23)--HIPR1(28): not used.					  */
/*HIPR1(29): (D=0.4 fm) the minimum distance between two nucleons	  */
/*	   inside a nucleus when the coordinates of the nucleons	  */
/*	   inside a nucleus are initialized.				  */
/*HIPR1(30): (D=2$\times$HIPR1(31)=57.0 mb) the inclusive cross 	  */
/*	   section $\sigma_{soft}$ for soft interactions. The default	  */
/*	   value $\sigma_{soft}=2\sigma_0$ is used to ensure the	  */
/*	   geometrical scaling of $pp$ interaction cross sections	  */
/*	   at low energies.						  */
/*HIPR1(31): (D=28.5 mb) the cross section $\sigma_0$ which		  */
/*	   characterizes the geometrical size of a nucleon		  */
/*	   ($\pi b_0^2=\sigma_0$, see Eq.~\ref{eq:over2}).		  */
/*	   The default value is only for high-energy			  */
/*	   limit ($\sqrt{s}>200$ GeV). At lower energies, a slight	  */
/*	   decrease which depends on energy is parametrized in the	  */
/*	   program. The default values of the two parameters		  */
/*	   HIPR1(30), HIPR1(31) are only for $NN$ type interactions.	  */
/*	   For other kinds of projectile or target hadrons, users	  */
/*	   should change these values so that correct inelastic 	  */
/*	   and total cross sections (HINT1(12), HINT1(13)) are		  */
/*	   obtained by the program.					  */
/*HIPR1(32): (D=3.90) parameter $\mu_0$ in Eq.~\ref{eq:over2} for	  */
/*	   the scaled eikonal function. 				  */
/*HIPR1(33): fractional cross section of single-diffractive		  */
/*	   interaction as parametrized in Ref.~\cite{goulianos}.	  */
/*HIPR1(34): maximum radial coordinate for projectile nucleons		  */
/*	   to be given by the initialization program HIJSET.		  */
/*HIPR1(35): maximum radial coordinate for target nucleons		  */
/*	   to be given by the initialization program HIJSET.		  */
/*HIPR1(36)-HIPR1(39): not used.					  */
/*HIPR1(40): (D=3.141592654) value of $\pi$.				  */
/*HIPR1(41)--HIPR1(42): not used.					  */
/*HIPR1(43): (D=0.01) fractional energy error relative to the		  */
/*	   colliding energy permitted per nucleon-nucleon collision.	  */
/*HIPR1(44), HIPR1(45), HIPR1(46): (D=1.5, 0.1 GeV, 0.25) parameters	  */
/*	   $\alpha$, $c$ and $\beta$ in the valence quark		  */
/*	   distributions for soft string excitation,			  */
/*	   $(1-x)^{\alpha}/(x^2+c^2/s)^{\beta}$ for baryons,		  */
/*	   $1/{(x^2+c^2/s)[(1-x)^2+c^2/s)]}^{\beta}$ for mesons.	  */
/*HIPR1(47), HIPR1(48): (D=0.0, 0.5) parameters $\alpha$ and $\beta$	  */
/*	   in valence quark distribution,				  */
/*	   $(1-x)^{\alpha}/(x^2+c^2/s)^{\beta}$, for the		  */
/*	   disassociated excitation in a single diffractive collision.	  */
/*HIPR1(49)--HIPR1(100): not used.					  */
/*IHPR2(1): (D=1) switch for dipole-approximated QCD radiation		  */
/*	   of the string system in soft interactions.			  */
/*IHPR2(2): (D=3) option for initial and final state radiation in	  */
/*	   the hard scattering. 					  */
/*	   =0: both initial and final radiation are off. 	          */
/*	   =1: initial radiation on and final radiation off.	          */
/*	   =2: initial radiation off and final radiation on.	          */
/*	   =3: both initial and final radiation are on.		          */
/*IHPR2(3): (D=0) switch for triggered hard scattering with specified	  */
/*	   $P_T\geq$HIPR1(10).						  */
/*	   =0: no triggered jet production.			          */
/*	   =1: ordinary hard processes.				          */
/*	   =2: only direct photon production.			          */
/*IHPR2(4): (D=1) switch for jet quenching in the excited		  */
/*	   nuclear matter.						  */
/*IHPR2(5): (D=1) switch for the $P_T$ kick due to soft interactions.	  */
/*IHPR2(6): (D=1) switch for the nuclear effect on the parton		  */
/*	   distribution function such as shadowing.			  */
/*IHPR2(7): (D=1) selection of Duke-Owens set (1 or 2) of parametrization */
/*	   of nucleon structure functions.				  */
/*IHPR2(8): (D=10) maximum number of hard scatterings per		  */
/*	   nucleon-nucleon interaction. When IHPR2(8)=0, jet		  */
/*	   production will be turned off. When IHPR2(8)$<0$, the	  */
/*	   number of jet production will be fixed at its absolute	  */
/*	   value for each NN collision. 				  */
/*IHPR2(9): (D=0) switch to guarantee at least one pair of minijets	  */
/*	   production per event ($pp$, $pA$ or $AB$).			  */
/*IHPR2(10): (D=0) option to print warning messages about errors that	  */
/*	   might happen. When a fatal error happens the current event	  */
/*	   will be abandoned and a new one is generated.		  */
/*IHPR2(11): (D=1) choice of baryon production model.			  */
/*	   =0: no baryon-antibaryon pair production, initial	          */
/*		   diquark treated as a unit.				  */
/*	   =1: diquark-antidiquark pair production allowed,	          */
/*		   initial diquark treated as a unit.			  */
/*	   =2: diquark-antidiquark pair production allowed,	          */
/*		   with the possibility for diquark to split		  */
/*		   according to the ``popcorn'' scheme (see the 	  */
/*		   documentation of JETSET 7.2).			  */
/*IHPR2(12): (D=1) option to turn off the automatic decay of the	  */
/*	    following particles:					  */
/*	   $\pi^0$, $K^0_S$, $D^{\pm}$, $\Lambda$, $\Sigma^{\pm}$.	  */
/*IHPR2(13): (D=1) option to turn on single diffractive reactions.	  */
/*IHPR2(14): (D=1) option to turn on elastic scattering.		  */
/*IHPR2(15)--IHPR2(18): not used.					  */
/*IHPR2(19): (D=1) option to turn on initial state soft interaction.	  */
/*IHPR2(20): (D=1) switch for the final fragmentation.			  */
/*IHPR2(21): (D=0) option to keep the information of all particles	  */
/*	     including those which have decayed and the decay history	  */
/*	     in the common block HIMAIN2. The line number of the parent   */
/*	     particle is KATT(I,3). The status of a partcile,		  */
/*	     whether it is a finally produced particle (KATT(I,4)=1)	  */
/*	     or a decayed particle (KATT(I,4)=11) is also kept. 	  */
/*IHPR2(22)-IHPR2(50): not used.					  */
/*HINT1(1): (GeV) colliding energy in the c.m. frame of nucleon-nucleon   */
/*	   collisions.							  */
/*HINT1(2): Lorentz transformation variable $\beta$ from laboratory	  */
/*	   to c.m.  frame of nucleon nucleon collisions.		  */
/*HINT1(3): rapidity $y_{cm}$ of the c.m. frame 			  */
/*	   $\beta=\tanh y_{cm}$.					  */
/*HINT1(4): rapidity of projectile nucleons (hadron) $y_{proj}$.	  */
/*HINT1(5): rapidity of target nucleons (hadron) $y_{targ}$.		  */
/*HINT1(6): (GeV) energy of the projectile nucleons (hadron) in the	  */
/*	   given frame. 						  */
/*HINT1(7): (GeV) energy of the target nucleons (hadron) in the 	  */
/*	   given frame. 						  */
/*HINT1(8): (GeV) the rest mass of projectile particles.		  */
/*HINT1(9): (GeV) the rest mass of target particles.			  */
/*HINT1(10): (mb) the averaged cross section for jet production 	  */
/*	   per nucleon-nucleon collisions,				  */
/*	   $\int d^2b\{1-\exp[-\sigma_{jet}T_N(b)]\}$.			  */
/*HINT1(11): (mb) the averaged inclusive cross section $\sigma_{jet}$	  */
/*	   for jet production per nucleon-nucleon collisions.		  */
/*HINT1(12): (mb) the averaged inelastic cross section of		  */
/*	   nucleon-nucleon collisions.					  */
/*HINT1(13): (mb) the averaged total cross section of nucleon-nucleon	  */
/*	   collisions.							  */
/*HINT1(14): (mb) the jet production cross section without nuclear	  */
/*	   shadowing effect $\sigma_{jet}^0$ (see Eq.~\ref{eq:sjetab}).   */
/*HINT1(15): (mb) the cross section $\sigma_{jet}^A$ to account for	  */
/*	   the projectile shadowing correction term in the jet cross	  */
/*	   section (see Eq.~\ref{eq:sjetab}).				  */
/*HINT1(16): (mb) the cross section $\sigma_{jet}^B$ to account for	  */
/*	   the target shadowing correction term in the jet cross	  */
/*	   section (see Eq.~\ref{eq:sjetab}).				  */
/*HINT1(17): (mb) the cross section $\sigma_{jet}^{AB}$ to account	  */
/*	   for the cross term of shadowing correction in the jet	  */
/*	   cross section.						  */
/*HINT1(18): (mb) the effective cross section				  */
/*	   $\sigma_{jet}^{eff}(r_A,r_B)$ for jet production		  */
/*	   of the latest nucleon-nucleon collision which depends	  */
/*	   on the transverse coordinates of the colliding		  */
/*	   nucleons.							  */
/*HINT1(19): (fm) the (absolute value of) impact parameter of the	  */
/*	   latest event.						  */
/*HINT1(20): (radians) the azimuthal angle $\phi$ of the impact 	  */
/*	   parameter vector in the transverse plane of the latest	  */
/*	   event.							  */
/*HINT1(21)--HINT1(25): the four momentum and mass ($p_x,p_y,p_z,E,M$)	  */
/*	   (GeV/$c$, GeV, GeV/$c^2$) of the first scattered parton	  */
/*	   in the triggered hard scattering. This is before the final	  */
/*	   state radiation but after the initial state radiation.	  */
/*HINT1(26)--HINT1(30): not used.					  */
/*HINT1(31)--HINT1(35): the four momentum and mass ($p_x,p_y,p_z,E,M$)	  */
/*	   (GeV/$c$, GeV, GeV/$c^2$) of the second scattered parton	  */
/*	   in the triggered hard scattering. This is before the final	  */
/*	   state radiation but after the initial state radiation.	  */
/*HINT1(46)--HINT1(40): not used.					  */
/*HINT1(41)--HINT1(45): the four momentum and mass ($p_x,p_y,p_z,E,M$)	  */
/*	   (GeV/$c$, GeV, GeV/$c^2$) of the first scattered parton	  */
/*	   in the latest hard scattering of the latest event.		  */
/*HINT1(46): $P_T$ (GeV/$c$) of the first scattered parton in the	  */
/*	   latest hard scattering of the latest event.			  */
/*HINT1(47)--HINT1(50): not used.					  */
/*HINT1(51)--HINT1(55): the four momentum and mass ($p_x,p_y,p_z,E,M$)	  */
/*	   (GeV/$c$, GeV, GeV/$c^2$) of the second scattered parton	  */
/*	   in the latest hard scattering of the latest event.		  */
/*HINT1(56): $P_T$ (GeV/$c$) of the second scattered parton in the	  */
/*	   latest hard scattering of the latest event.			  */
/*HINT1(57)--HINT1(58): not used.					  */
/*HINT1(59): (mb) the averaged cross section of the			  */
/*	   triggered jet production (with $P_T$ specified by HIPR1(10)	  */
/*	   and with switch by IHPR2(3)) per nucleon-nucleon		  */
/*	   collision,							  */
/*	   $\int d^2b\{1-\exp[-\sigma_{jet}^{trig}T_N(b)]\}$		  */
/*HINT1(60): (mb) the averaged inclusive cross section of the		  */
/*	   triggered jet production $\sigma_{jet}^{trig}$		  */
/*	   (with $P_T$ specified by					  */
/*	   HIPR1(10) and with switch by IHPR2(3)) per			  */
/*	   nucleon-nucleon collision.					  */
/*HINT1(61): (mb) the triggered jet production cross section without	  */
/*	   nuclear shadowing effect (similar to HINT1(14)).		  */
/*HINT1(62): (mb) the cross section to account for the projectile	  */
/*	   shadowing correction term in the triggered jet cross 	  */
/*	   section (similar to HINT1(15)).				  */
/*HINT1(63): (mb) the cross section to account for the target		  */
/*	   shadowing correction term in the triggered jet cross 	  */
/*	   section (similar to HINT1(16)).				  */
/*HINT1(64): (mb) the cross section to account for the cross		  */
/*	   term of shadowing correction in the triggered jet		  */
/*	   cross section (similar to HINT1(17). 			  */
/*HINT1(65): (mb) the inclusive cross section for latest triggered	  */
/*	   jet production which depends on the transverse coordinates	  */
/*	   of the colliding nucleons (similar to HINT1(18)).		  */
/*HINT1(67)--HINT1(71): not used.					  */
/*HINT1(72)--HINT1(75): three parameters for the Wood-Saxon		  */
/*	   projectile nuclear distribution and the normalization	  */
/*	   read from a table inside the program,			  */
/*	   $\rho(r)=C[1+W(r/R_A)^2]/\{1+\exp[(r-R_A)/D]\}$,		  */
/*	   $R_A$=HINT1(72), $D$=HINT1(73), $W$=HINT1(74), $C$=HINT1(75).  */
/*HINT1(76)--HINT1(79): three parameters for the Wood-Saxon		  */
/*	   projectile nuclear distribution and the normalization	  */
/*	   read from a table inside the program,			  */
/*	   $\rho(r)=C[1+W(r/R_A)^2]/\{1+\exp[(r-R_A)/D]\}$,		  */
/*	   $R_A$=HINT1(76), $D$=HINT1(77), $W$=HINT1(78), $C$=HINT1(79).  */
/*HINT1(80)--HINT1(100): the probability of $j=0-20$ number of hard	  */
/*	   scatterings per nucleon-nucleon collisions.			  */
/*IHNT2(1): the mass number of the projectile nucleus (1 for a hadron).   */
/*IHNT2(2): the charge number of the projectile nucleus. If the 	  */
/*	   projectile is a hadron, it gives the charge of the hadron.	  */
/*IHNT2(3): the mass number of the target nucleus (1 for a hadron).	  */
/*IHNT2(4): the charge number of the target nucleus. If the target	  */
/*	   is a hadron, it gives the charge of the hadron.		  */
/*IHNT2(5): the flavor code of the projectile hadron (0 for nucleus).	  */
/*IHNT2(6): the flavor code of the target hadron (0 for nucleus).	  */
/*IHNT2(7)--IHNT2(8): not used. 					  */
/*IHNT2(9): the flavor code of the first scattered parton in the	  */
/*	   triggered hard scattering.					  */
/*IHNT2(10): the flavor code of the second scattered parton in the	  */
/*	   triggered hard scattering.					  */
/*IHNT2(11): the sequence number of the projectile nucleon in the	  */
/*	   latest nucleon-nucleon interaction of the latest event.	  */
/*IHNT2(12): the sequence number of the target nucleon in the latest	  */
/*	   nucleon-nucleon interaction of the latest event.		  */
/*IHNT2(13): status of the latest soft string excitation.		  */
/*	   =1: double diffractive.				          */
/*	   =2: single diffractive.				          */
/*	   =3: non-single diffractive.				          */
/*IHNT2(14): the flavor code of the first scattered parton in the	  */
/*	   latest hard scattering of the latest event.			  */
/*IHNT2(15): the flavor code of the second scattered parton in the	  */
/*	   latest hard scattering of the latest event.			  */
/*IHNT2(16)--IHNT2(50): not used.					  */
/*========================================================================*/

/*========================================================================*/
/* COMMON/HIMAIN1/ NATT,EATT,JATT,NT,NP,N0,N01,N10,N11,BB	          */
/*------------------------------------------------------------------------*/
typedef struct {
   Int_t      natt;
   Float_t    eatt;
   Int_t      jatt;
   Int_t      nt;
   Int_t      np;
   Int_t      n0;
   Int_t      n01;
   Int_t      n10;
   Int_t      n11;
   Float_t    bb;   
   } Himain1Common;

#define HIMAIN1 COMMON_BLOCK(HIMAIN1,himain1)
COMMON_BLOCK_DEF(Himain1Common,HIMAIN1);
/*************************************************************************/
/*           D E S C R I P T I O N :                                     */
/*-----------------------------------------------------------------------*/
/*COMMON/HIMAIN1/NATT, EATT, JATT, NT, NP, N0, N01, N10, N11             */
/*Purpose: to give the overall information of the generated event.	 */
/*NATT: total number of produced stable and undecayed particles of	 */
/*	   the current event.						 */
/*EATT: the total energy of the produced particles in c.m. frame	 */
/*	   of the collision to check energy conservation.		 */
/*JATT: the total number of hard scatterings in the current event.	 */
/*NP, NT: the number of participant projectile and target nucleons	 */
/*	   in the current event.					 */
/*N0, N01, N10, N11: number of $N$-$N$, $N$-$N_{wounded}$,		 */
/*	   $N_{wounded}$-$N$, and					 */
/*	   $N_{wounded}$-$N_{wounded}$ collisions in			 */
/*	   the current event ($N$, $N_{wounded}$ stand			 */
/*	   for nucleon and wounded nucleon respectively).		 */
/*=======================================================================*/


/*========================================================*/
/* COMMON/HIMAIN2/KATT(130000,4),PATT(130000,4)           */
/*--------------------------------------------------------*/
typedef struct {
   Int_t    katt[4][200000];
   Float_t  patt[4][200000];
} Himain2Common;

#define HIMAIN2 COMMON_BLOCK(HIMAIN2,himain2)
COMMON_BLOCK_DEF(Himain2Common,HIMAIN2);
/*************************************************************************/
/*           D E S C R I P T I O N :                                     */
/*-----------------------------------------------------------------------*/
/*Purpose: to give information of produced stable and undecayed          */
/*	   particles. Parent particles which decayed are not included	 */
/*	   here.							 */
/*KATT(I, 1): (I=1,$\cdots$,NATT) flavor codes (see appendix) of	 */
/*	   the produced particles.					 */
/*KATT(I, 2): (I=1,$\cdots$,NATT) status codes to identify the		 */
/*	   sources from which the particles come.			 */
/*	   =0: projectile nucleon (or hadron) which has 		 */
/*	    not interacted at all.					 */
/*	   =1: projectile nucleon (or hadron) which			 */
/*	    only suffers an elastic collision.				 */
/*	   =2: from a diffractive projectile nucleon (or hadron)	 */
/*	    in a single diffractive interaction.			 */
/*	   =3: from the fragmentation of a projectile string		 */
/*	    system (including gluon jets).				 */
/*	   =10 target nucleon (or hadron) which has not 		 */
/*	    interacted at all.						 */
/*	   =11: target nucleon (or hadron) which only			 */
/*	    suffers an elastic collision.				 */
/*	   =12: from a diffractive target nucleon (or hadron)		 */
/*	    in a single diffractive interaction.			 */
/*	   =13: from the fragmentation of a target string		 */
/*	    system (including gluon jets).				 */
/*	   =20: from scattered partons which form string		 */
/*	    systems themselves. 					 */
/*	   =40: from direct production in the hard processes		 */
/*		   ( currently, only direct photons are included).	 */
/*KATT(I,3): (I=1,$\cdots$,NATT) line number of the parent particle.	 */
/*		   For finally produced or directly produced (not from	 */
/*		   the decay of another particle) particles, it is set	 */
/*		   to 0 (The option to keep the information of all	 */
/*		   particles including the decayed ones is IHPR2(21)=1). */
/*KATT(I,4): (I=1,$\cdots$,NATT) status number of the particle. 	 */
/*	  =1: finally or directly produced particles.			 */
/*	  =11: particles which has already decayed.			 */
/*PATT(I, 1-4): (I=1,$\cdots$,NATT) four-momentum ($p_x,p_y,p_z,E$)	 */
/*	   (GeV/$c$, GeV) of the produced particles.			 */
/*									 */
/*=======================================================================*/

/*=======================================================================*/
/* COMMON/HIJJET1/NPJ(300),KFPJ(300,500),PJPX(300,500),PJPY(300,500)     */
/*     &         ,PJPZ(300,500),PJPE(300,500),PJPM(300,500)              */
/*     &	 ,NTJ(300),KFTJ(300,500),PJTX(300,500),PJTY(300,500)     */
/*     &	 ,PJTZ(300,500),PJTE(300,500),PJTM(300,500)              */
/*-----------------------------------------------------------------------*/
typedef struct {
   Int_t    npj[300];
   Int_t    kfpj[500][300];
   Float_t  pjpx[500][300];
   Float_t  pjpy[500][300];
   Float_t  pjpz[500][300];
   Float_t  pjpe[500][300];
   Float_t  pjpm[500][300];
   Int_t    ntj[300];
   Int_t    kftj[500][300];
   Float_t  pjtx[500][300];
   Float_t  pjty[500][300];
   Float_t  pjtz[500][300];
   Float_t  pjte[500][300];
   Float_t  pjtm[500][300];
} Hijjet1Common;

#define HIJJET1 COMMON_BLOCK(HIJJET1,hijjet1)
COMMON_BLOCK_DEF(Hijjet1Common,HIJJET1);
/*************************************************************************/
/*	     D E S C R I P T I O N :			                 */
/*-----------------------------------------------------------------------*/
/*Purpose: contains information about produced partons which are         */
/*	   connected with the valence quarks and diquarks of		 */
/*	   projectile or target nucleons (or hadron) to form		 */
/*	   string systems for fragmentation. The momentum and		 */
/*	   energy of all produced partons are calculated in		 */
/*	   the c.m. frame of the collision. IAP, IAT are the		 */
/*	   numbers of nucleons in projectile and target nucleus 	 */
/*	   respectively (IAP, IAT=1 for hadron projectile or target).	 */
/*NPJ(I): (I=1,$\cdots$,IAP) number of partons associated with projectile*/
/*	   nucleon I.							 */
/*KFPJ(I, J): (I=1,$\cdots$,IAP, J=1,$\cdots$,NPJ(I)) parton		 */
/*	   flavor code of the						 */
/*	   parton J associated with projectile nucleon I.		 */
/*PJPX(I, J), PJPY(I, J), PJPZ(I, J), PJPE(I, J), PJPM(I, J): the four	 */
/*	   momentum and mass ($p_x,p_y,p_z,E,M$)			 */
/*	   (GeV/$c$, GeV, GeV/$c^2$) of parton J associated with	 */
/*	   the projectile nucleon I.					 */
/*NTJ(I): (I=1,$\cdots$,IAT) number of partons associated with		 */
/*	   target nucleon I.						 */
/*KFTJ(I, J): (I=1,$\cdots$,IAT, J=1,$\cdots$,NTJ(I)): parton		 */
/*	   flavor code of the  parton J associated with 		 */
/*	   target nucleon I.						 */
/*PJTX(I, J), PJTY(I, J), PJTZ(I, J), PJTE(I, J), PJTM(I, J): the four	 */
/*	   momentum and mass ($p_x,p_y,p_z,E,M$)			 */
/*	   (GeV/$c$, GeV, GeV/$c^2$) of parton J associated with	 */
/*	   target nucleon I.						 */
/*									 */
/*=======================================================================*/

/*=======================================================================*/
/* COMMON/HIJJET2/NSG,NJSG(900),IASG(900,3),K1SG(900,100),               */
/*     &	  K2SG(900,100),PXSG(900,100),PYSG(900,100),             */
/*     &	  PZSG(900,100),PESG(900,100),PMSG(900,100)              */ 
/*-----------------------------------------------------------------------*/
typedef struct {
   Int_t      nsg;
   Int_t      njsg[900];
   Int_t      iasg[3][900];
   Int_t      k1sg[100][900];
   Int_t      k2sg[100][900];
   Float_t    pxsg[100][900];
   Float_t    pysg[100][900];
   Float_t    pzsg[100][900];
   Float_t    pesg[100][900];
   Float_t    pmsg[100][900];
} Hijjet2Common;

#define HIJJET2 COMMON_BLOCK(HIJJET2,hijjet2)
COMMON_BLOCK_DEF(Hijjet2Common,HIJJET2);
/*************************************************************************/
/*           D E S C R I P T I O N :                                     */
/*-----------------------------------------------------------------------*/
/*Purpose: contains information about the produced partons which         */
/*	   will form string systems themselves without being		 */
/*	   connected to valence quarks and diquarks.			 */
/*NSG: the total number of such string systems. 			 */
/*NJSG(I): (I=1,$\cdots$,NSG) number of partons in the string system I.  */
/*IASG(I, 1), IASG(I, 2): to specify which projectile and target	 */
/*	   nucleons produce string system I.				 */
/*IASG(I, 3): to indicate whether the jets will be quenched (0) 	 */
/*	   or will not be quenched (1). 				 */
/*K1SG(I, J): (J=1,$\cdots$,NJSG(I)) color flow information of parton J  */
/*	   in string system I (see JETSET 7.2 for detailed		 */
/*	   explanation).						 */
/*K2SG(I, J): (J=1,$\cdots$,NJSG(I)) flavor code of parton J in string	 */
/*	   system I.							 */
/*PXSG(I, J), PYSG(I, J), PZSG(I, J), PESG(I, J), PMSG(I, J): four	 */
/*	   momentum and mass ($p_x,p_y,p_z,E,M$)			 */
/*	   ( GeV/$c$, GeV, GeV/$c^2$) of parton J in string system I.	 */
/*=======================================================================*/

/*=======================================================================*/
/* COMMON/HISTRNG/NFP(300,15),PP(300,15),NFT(300,15),PT(300,15)          */
/*-----------------------------------------------------------------------*/
typedef struct {
   Int_t    nfp[15][300];
   Float_t  pp[15][300];
   Int_t    nft[15][300];
   Float_t  pt[15][300];
} HistrngCommon;

#define HISTRNG COMMON_BLOCK(HISTRNG,histrng)
COMMON_BLOCK_DEF(HistrngCommon,HISTRNG);
/*************************************************************************/
/*           D E S C R I P T I O N :                                     */
/*-----------------------------------------------------------------------*/
/*Purpose: contains information about the projectile and                 */
/* target nucleons (hadron) and the corresponding constituent		 */
/* quarks, diquarks. IAP, IAT are the numbers of nucleons in		 */
/* projectile and target nucleus respectively (IAP, IAT=1		 */
/* for hadron projectile or target).					 */
/*NFP(I, 1): (I=1,$\cdots$,IAP) flavor code of the valence quark in	 */
/*	   projectile nucleon (hadron) I.				 */
/*NFP(I, 2): flavor code of diquark in projectile nucleon (anti-quark	 */
/*	   in projectile meson) I.					 */
/*NFP(I, 3): present flavor code of the projectile nucleon (hadron) I	 */
/*	   ( a nucleon or meson can be excited to its vector resonance). */
/*NFP(I, 4): original flavor code of projectile nucleon (hadron) I.	 */
/*NFP(I, 5): collision status of projectile nucleon (hadron) I. 	 */
/*	   =0: suffered no collision.					 */
/*	   =1: suffered an elastic collision.				 */
/*	   =2: being the diffractive one in a single-diffractive	 */
/*	    collision.							 */
/*	   =3: became an excited string after an inelastic		 */
/*	    collision.							 */
/*NFP(I, 6): the total number of hard scatterings associated with	 */
/*	   projectile nucleon (hadron) I. If NFP(I,6)$<0$, it can not	 */
/*		   produce jets any more due to energy conservation.	 */
/*NFP(I, 10): to indicate whether the valence quarks or diquarks	 */
/*	   (anti-quarks) in projectile nucleon (hadron) I		 */
/*	   suffered a hard scattering,					 */
/*	   =0: has not  suffered a hard scattering.			 */
/*	   =1: suffered one or more hard scatterings in 		 */
/*	    current binary nucleon-nucleon collision.			 */
/*	   =-1: suffered one or more hard scatterings in		 */
/*	    previous binary nucleon-nucleon collisions. 		 */
/*NFP(I, 11): total number of interactions projectile nucleon (hadron)	 */
/*	   I  has suffered so far.					 */
/*PP(I, 1), PP(I, 2), PP(I, 3), PP(I, 4), PP(I, 5): four momentum and	 */
/*	   the invariant mass ($p_x,p_y,p_z,E,M$)			 */
/*	   (GeV/$c$, GeV, GeV/$c^2$) of projectile nucleon (hadron) I.	 */
/*PP(I, 6), PP(I, 7): transverse momentum ($p_x,p_y$) (GeV/$c$) of the	 */
/*	   valence quark in projectile nucleon (hadron) I.		 */
/*PP(I, 8), PP(I, 9): transverse momentum ($p_x,p_y$) (GeV/$c$) of the	 */
/*	   diquark (anti-quark) in projectile nucleon (hadron) I.	 */
/*PP(I, 10), PP(I, 11), PP(I, 12): three momentum ($p_x,p_y,p_z$)	 */
/*	   (GeV/$c$) transferred to the quark or diquark (anti-quark)	 */
/*	   in projectile nucleon (hadron) I from the last hard		 */
/*	   scattering.							 */
/*PP(I, 14): mass (GeV/$c^2$) of the quark in projectile nucleon	 */
/*	   (hadron) I.							 */
/*PP(I, 15): mass of the diquark (anti-quark) in projectile		 */
/*	   nucleon (hadron) I.						 */
/*NFT(I, 1--15), PT(I,1--15): give the same				 */
/*	   information for the target nucleons (hadron) and the 	 */
/*	   corresponding quarks and diquarks (anti-quarks) as for	 */
/*	   the projectile nucleons.					 */
/*									 */
/*=======================================================================*/
}

#endif
