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
#include "Fpart.h"  
//#include "Fcaslim.h"  //(CASLIM) fluka common

//Virutal MC
#include "AliGenerator.h"
#include "AliStack.h"
#include "../EVGEN/AliGenHaloProtvino.h"

#include "TVirtualMCStack.h"
#include "TParticle.h"
#include "TVector3.h"
#include "TRandom.h"

//Other
#include <Riostream.h>

#ifndef WIN32
# define source source_
# define geocrs geocrs_
# define georeg georeg_
# define geohsm geohsm_
# define soevsv soevsv_
# define mcihad mcihad_
# define source_halo source_halo__
#else
# define source SOURCE
# define geocrs GEOCRS
# define georeg GEOREG
# define geohsm GEOHSM
# define soevsv SOEVSV
# define mcihad MCIHAD
# define source_halo SOURCE_HALO
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
  int  type_of_call mcihad(const int&);
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

  void source_halo(Int_t& nomore) {

      static Bool_t lfirst               = true;
      static AliGenHaloProtvino* gener   = 0;
      static AliStack* stack             = 0;	

      nomore = 0;
      TParticle* particle;
      Int_t itrack = -1;
      if (lfirst) {
	  EPISOR.tkesum = zerzer;
	  lfirst = false;
	  EPISOR.lussrc = true;
//
// The generator
//
	  gener  = new AliGenHaloProtvino(10);
	  gener->SetFileName("/home/morsch/AliRoot/newio/data/ir.2.ring.one.b.10.m");
	  gener->SetSide(2);
	  gener->SetRunPeriod();
	  
	  
// The particle stack
	  stack = new AliStack(1000);
	  gener->SetStack(stack);
	  gener->Init();

      } else {
	  //
	  // Generate event
	  stack->Reset();
	  gener->Generate();
	  Int_t npart = stack->GetNprimary();
	  for (Int_t part = 1; part < npart; part += 2) {
	      particle = stack->Particle(part);
	      STACK.lstack++;
	      Int_t pdg = particle->GetPdgCode();
	      Int_t intfluka = mcihad(pdg);
	      Int_t ifl = GetFlukaKPTOIP(intfluka);
	      Float_t wgt = particle->GetWeight();
	      
	      printf("Particle %5d %10s %10.3f %10.3f %10.3f %10.3f\n", 
		     STACK.lstack, 
		     particle->GetName(), 
		     particle->Px(), particle->Py(), particle->Pz(), wgt);
	      

	      
	      /* Wt is the weight of the particle*/
	      STACK.wt[STACK.lstack] = wgt;
	      STARS.weipri += STACK.wt[STACK.lstack];
	      
	      STACK.ilo[STACK.lstack] = ifl;
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
	      STACK.ispark[STACK.lstack][mkbmx2-1] = itrack;
	      STACK.nparma++;
	      STACK.numpar[STACK.lstack] = STACK.nparma;
	      STACK.nevent[STACK.lstack] = 0;
	      STACK.dfnear[STACK.lstack] = +zerzer;
	      
	      /* Particle age (s)*/
	      STACK.agestk[STACK.lstack] = +zerzer;
	      STACK.aknshr[STACK.lstack] = -twotwo;
	      
	      /* Group number for "low" energy neutrons, set to 0 anyway*/
	      STACK.igroup[STACK.lstack] = 0;
	      
	      /* Kinetic energy */
	      STACK.tke[STACK.lstack] = particle->Energy() - particle->GetMass();
	      
	      
	      /* Particle momentum*/
	      STACK.pmom [STACK.lstack] = particle->P();
	      
	      /* Cosines (tx,ty,tz)*/
	      Double_t cosx = particle->Px()/particle->P();
	      Double_t cosy = particle->Py()/particle->P();
	      Double_t cosz = TMath::Sqrt(oneone - cosx*cosx - cosy*cosy);
	      if (particle->Pz() < 0.)   cosz = -cosz;
	      
	      STACK.tx [STACK.lstack] =   cosx;
	      STACK.ty [STACK.lstack] =   cosy;
	      STACK.tz [STACK.lstack] =   cosz;
	      
	      /* Polarization cosines:*/
	      
	      STACK.txpol [STACK.lstack] = -twotwo;
	      STACK.typol [STACK.lstack] = +zerzer;
	      STACK.tzpol [STACK.lstack] = +zerzer;
	      
	      /* Particle coordinates*/
	      STACK.xa [STACK.lstack] =   particle->Vx();
	      STACK.ya [STACK.lstack] =   particle->Vy();
	      STACK.za [STACK.lstack] =   particle->Vz();
	      
	      
	      
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
	      printf("Particle Vertex %10.3f %10.3f %10.3f %5d\n",  
		     STACK.xa [STACK.lstack],  STACK.ya [STACK.lstack],  STACK.za [STACK.lstack], 
		     STACK.nreg[STACK.lstack]);
	      

	      /*  Do not change these cards:*/
	      Int_t igeohsm1 = 1;
	      Int_t igeohsm2 = -11;
	      geohsm ( STACK.nhspnt[STACK.lstack], igeohsm1, igeohsm2, LTCLCM.mlattc );
	      STACK.nlattc[STACK.lstack] = LTCLCM.mlattc;
	      soevsv();
	  }
      }
  }
}

    
