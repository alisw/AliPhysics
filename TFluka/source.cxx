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
#include "Fopphst.h"  //(OPPHST) fluka common
//#include "Fcaslim.h"  //(CASLIM) fluka common

//Virutal MC
#include "TFluka.h"

#include "TVirtualMCStack.h"
//#include "TVirtualMCApplication.h"

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
// Get the pointer to TFluka
    TFluka* fluka = (TFluka*)gMC;

    Int_t verbosityLevel = fluka->GetVerbosityLevel();
    Bool_t debug = (verbosityLevel>=3)?kTRUE:kFALSE;
    if (debug) {
      cout << "==> source(" << nomore << ")" << endl;
      cout << "\t* EPISOR.lsouit = " << (EPISOR.lsouit?'T':'F') << endl;
    }  

    static Bool_t lfirst = true;
    static Bool_t particleIsPrimary = true;
    static Bool_t lastParticleWasPrimary = true;

    nomore = 0;
    
//  Get the stack 
    TVirtualMCStack* cppstack = fluka->GetStack();

    TParticle* particle;
    Int_t itrack = -1;
    Int_t  nprim  = cppstack->GetNprimary();
//  Get the next particle from the stack
    particle  = cppstack->PopNextTrack(itrack);
    fluka->SetTrackIsNew(kTRUE);

//  Is this a secondary not handled by Fluka, i.e. a particle added by user action ?
    lastParticleWasPrimary = particleIsPrimary;
    
    if (itrack >= nprim) {
	particleIsPrimary = kFALSE;
    } else {
	particleIsPrimary = kTRUE;
    }

    if (lfirst) {
	EPISOR.tkesum = zerzer;
	lfirst = false;
	EPISOR.lussrc = true;
    } else {
//
// Post-track actions for primary track
//
	if (particleIsPrimary) {
	    TVirtualMCApplication::Instance()->PostTrack();
	    TVirtualMCApplication::Instance()->FinishPrimary();
	    if ((itrack%10)==0) printf("=== TRACKING PRIMARY %d ===\n", itrack);
	}
    }

    // Exit if itrack is negative (-1). Set lsouit to false to mark last track for this event

    if (itrack<0) {
      nomore = 1;
      EPISOR.lsouit = false;
      if (debug) {
         cout << "\t* EPISOR.lsouit = " << (EPISOR.lsouit?'T':'F') << endl;
         cout << "\t* No more particles. Exiting..." << endl;
         cout << "<== source(" << nomore << ")" << endl;
      }   
      return;
    }
    
    //
    // Handle user event abortion
    if (fluka->EventIsStopped()) {
	printf("Event has been stopped by user !");
	fluka->SetStopEvent(kFALSE);
	nomore = 1;
	EPISOR.lsouit = false;
	return;
    }

    //Get some info about the particle and print it
    //
    //pdg code
    Int_t pdg = particle->GetPdgCode();
    TVector3 polarisation;
    particle->GetPolarisation(polarisation);
    if (debug) {
       cout << "\t* Particle " << itrack << " retrieved..." << endl;
       cout << "\t\t+ Name = " << particle->GetName() << endl;
       cout << "\t\t+ PDG/Fluka code = " << pdg 
	    << " / " << fluka->IdFromPDG(pdg) << endl;
       cout << "\t\t+ P = (" 
	    << particle->Px() << " , "
	    << particle->Py() << " , "
	    << particle->Pz() << " ) --> "
	    << particle->P() << " GeV " 
	    << particle->Energy() << " GeV "  
	    << particle->GetMass() << " GeV " << endl;
    }   
    /* Lstack is the stack counter: of course any time source is called it
     * must be =0
     */
    /* Cosines (tx,ty,tz)*/
    Double_t cosx = particle->Px()/particle->P();
    Double_t cosy = particle->Py()/particle->P();
    Double_t cosz = TMath::Sqrt(oneone - cosx*cosx - cosy*cosy);
    if (particle->Pz() < 0.) cosz = -cosz;    

    //STACK.ilo[STACK.lstack] = BEAM.ijbeam;
    if (pdg != 50000050 &&  pdg !=  50000051) {
	STACK.lstack++;
	Int_t ifl =  fluka-> IdFromPDG(pdg);
	STACK.ilo[STACK.lstack] = ifl;
	/* Wt is the weight of the particle*/
	STACK.wt[STACK.lstack] = oneone;
	STARS.weipri += STACK.wt[STACK.lstack];
	
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
	STACK.ispark[STACK.lstack][mkbmx2-1] = itrack;
	STACK.nparma++;
	STACK.numpar[STACK.lstack] = STACK.nparma;
	STACK.nevent[STACK.lstack] = 0;
	STACK.dfnear[STACK.lstack] = +zerzer;
	
	/* Particle age (s)*/
	STACK.agestk[STACK.lstack] = +zerzer;
	STACK.cmpath[STACK.lstack] = +zerzer;
	STACK.aknshr[STACK.lstack] = -twotwo;

	/* Group number for "low" energy neutrons, set to 0 anyway*/
	STACK.igroup[STACK.lstack] = 0;
	
	/* Kinetic energy */
	Double_t p    = particle->P();
	Double_t mass = PAPROP.am[ifl + 6];
	STACK.tke[STACK.lstack] =  TMath::Sqrt( p * p + mass * mass) - mass;
	/* Particle momentum*/
	STACK.pmom [STACK.lstack] = p;

	STACK.tx [STACK.lstack] = cosx;
	STACK.ty [STACK.lstack] = cosy;
	STACK.tz [STACK.lstack] = cosz;
    
	/* Polarization cosines:*/
	if (polarisation.Mag()) {
	    Double_t cospolx = polarisation.Px() / polarisation.Mag();
	    Double_t cospoly = polarisation.Py() / polarisation.Mag();
	    Double_t cospolz = sqrt(oneone - cospolx * cospolx - cospoly * cospoly);
	    STACK.txpol [STACK.lstack] = cospolx;
	    STACK.typol [STACK.lstack] = cospoly;
	    STACK.tzpol [STACK.lstack] = cospolz;
	}
	else {
	    STACK.txpol [STACK.lstack] = -twotwo;
	    STACK.typol [STACK.lstack] = +zerzer;
	    STACK.tzpol [STACK.lstack] = +zerzer;
	}
	
	/* Particle coordinates*/
	// Vertext coordinates;
	STACK.xa [STACK.lstack] = particle->Vx();
	STACK.ya [STACK.lstack] = particle->Vy();
	STACK.za [STACK.lstack] = particle->Vz();
    
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
    } else {
        //
	// Next particle is optical photon
	//
	OPPHST.lstopp++;
	OPPHST.donear [OPPHST.lstopp - 1] = 0.;

	OPPHST.xoptph [OPPHST.lstopp - 1] = particle->Vx();
	OPPHST.yoptph [OPPHST.lstopp - 1] = particle->Vy();
	OPPHST.zoptph [OPPHST.lstopp - 1] = particle->Vz();

	OPPHST.txopph [OPPHST.lstopp - 1] = cosx;
	OPPHST.tyopph [OPPHST.lstopp - 1] = cosy;
	OPPHST.tzopph [OPPHST.lstopp - 1] = cosz;


	if (polarisation.Mag()) {
	    Double_t cospolx = polarisation.Px() / polarisation.Mag();
	    Double_t cospoly = polarisation.Py() / polarisation.Mag();
	    Double_t cospolz = sqrt(oneone - cospolx * cospolx - cospoly * cospoly);
	    OPPHST.txpopp [OPPHST.lstopp - 1] = cospolx;
	    OPPHST.typopp [OPPHST.lstopp - 1] = cospoly;
	    OPPHST.tzpopp [OPPHST.lstopp - 1] = cospolz;
	}
	else {
	    OPPHST.txpopp [OPPHST.lstopp - 1] = -twotwo;
	    OPPHST.typopp [OPPHST.lstopp - 1] = +zerzer;
	    OPPHST.tzpopp [OPPHST.lstopp - 1] = +zerzer;
	}
	
	geocrs( OPPHST.txopph[OPPHST.lstopp - 1], 
		OPPHST.tyopph[OPPHST.lstopp - 1], 
		OPPHST.tzopph[OPPHST.lstopp - 1] );
    
	Int_t idisc;

	georeg ( OPPHST.xoptph[OPPHST.lstopp - 1], 
		 OPPHST.yoptph[OPPHST.lstopp - 1], 
		 OPPHST.zoptph[OPPHST.lstopp - 1],
		 OPPHST.nregop[OPPHST.lstopp - 1], 
		 idisc);//<-- dummy return variable not used
	
	OPPHST.wtopph [OPPHST.lstopp - 1] = particle->GetWeight();
	OPPHST.poptph [OPPHST.lstopp - 1] = particle->P();
	OPPHST.agopph [OPPHST.lstopp - 1] = particle->T();	
	OPPHST.cmpopp [OPPHST.lstopp - 1] = +zerzer;
	OPPHST.loopph [OPPHST.lstopp - 1] = 0;
	OPPHST.louopp [OPPHST.lstopp - 1] = itrack;
    }
    
//
//  Pre-track actions at for primary tracks
//
    if (particleIsPrimary) {
	TVirtualMCApplication::Instance()->BeginPrimary();
	TVirtualMCApplication::Instance()->PreTrack();
    }
    
//
    if (debug) cout << "<== source(" << nomore << ")" << endl;
  }
}
