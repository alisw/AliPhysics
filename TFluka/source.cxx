#define METHODDEBUG

// Fortran 
#include "TCallf77.h"

// Fluka commons
#include "Fdblprc.h"  //(DBLPRC) fluka common
#include "Fdimpar.h"  //(DIMPAR) fluka parameters
#include "Fepisor.h"  //(EPISOR) fluka common
#include "Fstack.h"   //(STACK)  fluka common
#include "Fstars.h"   //(STARS)  fluka common
#include "Fbeam.h"    //(BEAM)   fluka common
#include "Fpaprop.h"  //(PAPROP) fluka common
#include "Fltclcm.h"  //(LTCLCM) fluka common
//#include "Fcaslim.h"  //(CASLIM) fluka common

//Virutal MC
#include "TFluka.h"
#include "TVirtualMCStack.h"
#include "TVirtualMCApplication.h"
#include "TParticle.h"
#include "TVector3.h"

//Other
#include <Riostream.h>

#ifndef WIN32
# define source source_
# define geocrs geocrs_
# define georeg georeg_
# define geohsm geohsm_
# define soevsv soevsv_
#else
# define source SOURCE
# define geocrs GEOCRS
# define georeg GEOREG
# define geohsm GEOHSM
# define soevsv SOEVSV
#endif

extern "C" {
  //
  // Prototypes for FLUKA functions
  //
  void type_of_call geocrs(Double_t &, Double_t &, Double_t &);
  void type_of_call georeg(Double_t &, Double_t &, Double_t &, 
			   Int_t &, Int_t &);
  void type_of_call geohsm(Int_t &, Int_t &, Int_t &, Int_t &);
  void type_of_call soevsv();
 /*
   *----------------------------------------------------------------------*
   *                                                                      *
   *     Created on 07 january 1990   by    Alfredo Ferrari & Paola Sala  *
   *                                                   Infn - Milan       *
   *                                                                      *
   *     Last change on 21-jun-98     by    Alfredo Ferrari               *
   *                                                                      *
   *     C++ version on 27-sep-02     by    Isidro Gonzalez               *
   *                                                                      *
   *  This is just an example of a possible user written source routine.  *
   *  note that the beam card still has some meaning - in the scoring the *
   *  maximum momentum used in deciding the binning is taken from the     *
   *  beam momentum.  Other beam card parameters are obsolete.            *
   *                                                                      *
   *----------------------------------------------------------------------*/

  void source(Int_t& nomore) {
#ifdef METHODDEBUG
    cout << "==> source(" << nomore << ")" << endl;
#endif

    cout << "\t* EPISOR.lsouit = " << (EPISOR.lsouit?'T':'F') << endl;

    static Bool_t lfirst = true;
    /*======================================================================*
     *                                                                      *
     *                 BASIC VERSION                                        *
     *                                                                      *
     *======================================================================*/
    nomore = 0;
    /*  +-------------------------------------------------------------------*
     *  |  First call initializations:*/
    if (lfirst) {

      /*|  *** The following 3 cards are mandatory ***/
      
      EPISOR.tkesum = zerzer;
      lfirst = false;
      EPISOR.lussrc = true;
      /*|  *** User initialization ***/
    }
    /*  |
     *  +-------------------------------------------------------------------*
     *  Push one source particle to the stack. Note that you could as well
     *  push many but this way we reserve a maximum amount of space in the
     *  stack for the secondaries to be generated
     */

    // Get the pointer to the VMC
    TVirtualMC* fluka = TFluka::GetMC();
    // Get the stack produced from the generator
    TVirtualMCStack* cppstack = fluka->GetStack();
    //Get next particle
    if (STACK.lstack != 1) {
	TVirtualMCApplication::Instance()->PostTrack();
	TVirtualMCApplication::Instance()->FinishPrimary();
    }
    Int_t itrack = -1;
    TParticle* particle = cppstack->GetNextTrack(itrack);

    //Exit if itrack is negative (-1). Set lsouit to false to mark last track for
    //this event
    if (itrack<0) {
      nomore = 1;
      EPISOR.lsouit = false;
      cout << "\t* EPISOR.lsouit = " << (EPISOR.lsouit?'T':'F') << endl;
      cout << "\t* No more particles. Exiting..." << endl;
#ifdef METHODDEBUG
      cout << "<== source(" << nomore << ")" << endl;
#endif
      return;
    }

    //Get some info about the particle and print it
    TVector3 polarisation;
    particle->GetPolarisation(polarisation);
    cout << "\t* Particle " << itrack << " retrieved..." << endl;
    cout << "\t\t+ Name = " << particle->GetName() << endl;
    cout << "\t\t+ PDG/Fluka code = " << particle->GetPdgCode() 
	 << " / " << fluka->IdFromPDG(particle->GetPdgCode()) << endl;
    cout << "\t\t+ E = " << particle->Energy() << " GeV" << endl;
    cout << "\t\t+ P = (" 
	 << particle->Px() << " , "
	 << particle->Py() << " , "
	 << particle->Pz() << " ) --> "
	 << particle->P() << " GeV" << endl;
    cout << "\t\t+ M = " << particle->GetMass() << " GeV" << endl;
    cout << "\t\t+ Initial point = ( " 
	 << particle->Vx() << " , "
	 << particle->Vy() << " , "
	 << particle->Vz() << " )"
	 << endl;    
    cout << "\t\t+ Polarisation = ( " 
	 << polarisation.Px() << " , "
	 << polarisation.Py() << " , "
	 << polarisation.Pz() << " )"
	 << endl;    
    /* Lstack is the stack counter: of course any time source is called it
     * must be =0
     */
    
    STACK.lstack++;
    cout << "\t* Storing particle parameters in the stack, lstack = " 
	 << STACK.lstack << endl;
    /* Wt is the weight of the particle*/
    STACK.wt[STACK.lstack] = oneone;
    STARS.weipri += STACK.wt[STACK.lstack];
    /* Particle type (1=proton.....). Ijbeam is the type set by the BEAM
     * card
       */
    //STACK.ilo[STACK.lstack] = BEAM.ijbeam;
    STACK.ilo[STACK.lstack] = fluka-> IdFromPDG(particle->GetPdgCode());
    /* From this point .....
     * Particle generation (1 for primaries)
       */
    STACK.lo[STACK.lstack] = 1;
    /* User dependent flag:*/
    STACK.louse[STACK.lstack] = 0;
    /* User dependent spare variables:*/
    Int_t ispr = 0;
    for (ispr = 0; ispr < mkbmx1; ispr++)
      STACK.sparek[STACK.lstack][ispr] = zerzer;
    /* User dependent spare flags:*/
    for (ispr = 0; ispr < mkbmx2; ispr++)
	STACK.ispark[STACK.lstack][ispr] = 0;
    /* Save the track number of the stack particle:*/
    STACK.ispark[STACK.lstack][mkbmx2-1] = STACK.lstack;
    STACK.nparma++;
    STACK.numpar[STACK.lstack] = STACK.nparma;
    STACK.nevent[STACK.lstack] = 0;
    STACK.dfnear[STACK.lstack] = +zerzer;
      /* ... to this point: don't change anything
       * Particle age (s)
       */
    STACK.agestk[STACK.lstack] = +zerzer;
    STACK.aknshr[STACK.lstack] = -twotwo;
    /* Group number for "low" energy neutrons, set to 0 anyway*/
    STACK.igroup[STACK.lstack] = 0;
    /* Kinetic energy of the particle (GeV)*/
    //STACK.tke[STACK.lstack] = 
    //sqrt( BEAM.pbeam*BEAM.pbeam + 
    // PAPROP.am[BEAM.ijbeam+6]*PAPROP.am[BEAM.ijbeam+6] ) 
    //- PAPROP.am[BEAM.ijbeam+6];
    STACK.tke[STACK.lstack] = particle->Energy() - particle->GetMass();
    
    /* Particle momentum*/
    //STACK.pmom [STACK.lstack] = BEAM.pbeam;
    STACK.pmom [STACK.lstack] = particle->P();
    
    /*     PMOM (lstack) = SQRT ( TKE (stack) * ( TKE (lstack) + TWOTWO
     *    &                     * AM (ILO(lstack)) ) )
     * Cosines (tx,ty,tz)
     */
    //STACK.tx [STACK.lstack] = BEAM.tinx;
    //STACK.ty [STACK.lstack] = BEAM.tiny;
    //STACK.tz [STACK.lstack] = BEAM.tinz;
    Double_t cosx = particle->Px()/particle->P();
    Double_t cosy = particle->Py()/particle->P();
    Double_t cosz = sqrt(oneone - cosx*cosx - cosy*cosy);
    STACK.tx [STACK.lstack] = cosx;
    STACK.ty [STACK.lstack] = cosy;
    STACK.tz [STACK.lstack] = cosz;
    
    /* Polarization cosines:
     */
    //STACK.txpol [STACK.lstack] = -twotwo;
    //STACK.typol [STACK.lstack] = +zerzer;
    //STACK.tzpol [STACK.lstack] = +zerzer;
    if (polarisation.Mag()) {
      Double_t cospolx = polarisation.Px()/polarisation.Mag();
      Double_t cospoly = polarisation.Py()/polarisation.Mag();
      Double_t cospolz = sqrt(oneone - cospolx*cospolx - cospoly*cospoly);
      STACK.tx [STACK.lstack] = cospolx;
      STACK.ty [STACK.lstack] = cospoly;
      STACK.tz [STACK.lstack] = cospolz;
    }
    else {
      STACK.txpol [STACK.lstack] = -twotwo;
      STACK.typol [STACK.lstack] = +zerzer;
      STACK.tzpol [STACK.lstack] = +zerzer;
    }
    
    /* Particle coordinates*/
    //STACK.xa [STACK.lstack] = BEAM.xina;
    //STACK.ya [STACK.lstack] = BEAM.yina;
    //STACK.za [STACK.lstack] = BEAM.zina
      //Vertext coordinates;
    STACK.xa [STACK.lstack] = particle->Vx();
    STACK.ya [STACK.lstack] = particle->Vy();
    STACK.za [STACK.lstack] = particle->Vz();
    
    // Some printout
    cout << "\t* Particle information transfered to stack..." << endl;
    
    /*  Calculate the total kinetic energy of the primaries: don't change*/
    Int_t st_ilo =  STACK.ilo[STACK.lstack];
    if ( st_ilo != 0 )
      EPISOR.tkesum += 
	((STACK.tke[STACK.lstack] + PAPROP.amdisc[st_ilo+6])
	 * STACK.wt[STACK.lstack]);
    else
      EPISOR.tkesum += (STACK.tke[STACK.lstack] * STACK.wt[STACK.lstack]);
    
    /*  Here we ask for the region number of the hitting point.
     *     NREG (LSTACK) = ...
     *  The following line makes the starting region search much more
     *  robust if particles are starting very close to a boundary:
     */
    geocrs( STACK.tx[STACK.lstack], 
	    STACK.ty[STACK.lstack], 
	    STACK.tz[STACK.lstack] );
    Int_t idisc;
    georeg ( STACK.xa[STACK.lstack], 
	     STACK.ya[STACK.lstack], 
	     STACK.za[STACK.lstack],
	     STACK.nreg[STACK.lstack], 
	     idisc);//<-- dummy return variable not used
    
    /*  Do not change these cards:*/
    Int_t igeohsm1 = 1;
    Int_t igeohsm2 = -11;
    geohsm ( STACK.nhspnt[STACK.lstack], igeohsm1, igeohsm2, LTCLCM.mlattc );
    STACK.nlattc[STACK.lstack] = LTCLCM.mlattc;
    soevsv();
    
    cout << "\t* EPISOR.lsouit = " << (EPISOR.lsouit?'T':'F') << endl;
    cout << "\t* " << STACK.lstack << " particles in the event" << endl;
    TVirtualMCApplication::Instance()->PreTrack();
#ifdef METHODDEBUG
    cout << "<== source(" << nomore << ")" << endl;
#endif
  }
}
