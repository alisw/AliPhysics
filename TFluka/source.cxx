// Fortran 
#include "TCallf77.h"

// Fluka commons
#include "Fdblprc.h"  //(DBLPRC) fluka common
#include "Fdimpar.h"  //(DIMPAR) fluka parameters
#include "Fsourcm.h"  //(EPISOR) fluka common
#include "Fflkstk.h"  //(FLKSTK)  fluka common
#include "Fsumcou.h"  //(SUMCOU)  fluka common
#include "Fpaprop.h"  //(PAPROP) fluka common
#include "Fltclcm.h"  //(LTCLCM) fluka common
#include "Fopphst.h"  //(OPPHST) fluka common

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
      cout << "\t* SOURCM.lsouit = " << (SOURCM.lsouit?'T':'F') << endl;
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
    if (itrack == (nprim - 1)) lfirst = true;
//  Is this a secondary not handled by Fluka, i.e. a particle added by user action ?
    lastParticleWasPrimary = particleIsPrimary;
    
    if (itrack >= nprim) {
        particleIsPrimary = kFALSE;
    } else {
        particleIsPrimary = kTRUE;
    }

    if (lfirst) {
        SOURCM.tkesum = zerzer;
        lfirst = false;
        SOURCM.lussrc = true;
    } else {
//
// Post-track actions for primary track
//
        if (particleIsPrimary) {
            TVirtualMCApplication::Instance()->PostTrack();
            TVirtualMCApplication::Instance()->FinishPrimary();
            if ((itrack%10)==0)
                cout << "=== TRACKING PRIMARY "<< itrack <<" ===" << endl;
            //printf("=== TRACKING PRIMARY %d ===\n", itrack);
        }
    }

    // Exit if itrack is negative (-1). Set lsouit to false to mark last track for this event

    if (itrack<0) {
      nomore = 1;
      SOURCM.lsouit = false;
      if (debug) {
         cout << "\t* SOURCM.lsouit = " << (SOURCM.lsouit?'T':'F') << endl;
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
        SOURCM.lsouit = false;
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
    /* Npflka is the stack counter: of course any time source is called it
     * must be =0
     */
    /* Cosines (tx,ty,tz)*/
    Double_t cosx = particle->Px()/particle->P();
    Double_t cosy = particle->Py()/particle->P();
    Double_t cosz = TMath::Sqrt(oneone - cosx*cosx - cosy*cosy);
    if (particle->Pz() < 0.) cosz = -cosz;    

    if (pdg != 50000050 &&  pdg !=  50000051) {
        FLKSTK.npflka++;
        Int_t ifl =  fluka-> IdFromPDG(pdg);
        FLKSTK.iloflk[FLKSTK.npflka] = ifl;
        /* Wtflk is the weight of the particle*/
        FLKSTK.wtflk[FLKSTK.npflka] = oneone;
        SUMCOU.weipri += FLKSTK.wtflk[FLKSTK.npflka];
        
        FLKSTK.loflk[FLKSTK.npflka] = 1;
        
        /* User dependent flag:*/
        FLKSTK.louse[FLKSTK.npflka] = 0;

        /* User dependent spare variables:*/
        Int_t ispr = 0;
        for (ispr = 0; ispr < mkbmx1; ispr++)
            FLKSTK.sparek[FLKSTK.npflka][ispr] = zerzer;

        /* User dependent spare flags:*/
        for (ispr = 0; ispr < mkbmx2; ispr++)
            FLKSTK.ispark[FLKSTK.npflka][ispr] = 0;

        /* Save the track number of the stack particle:*/
        FLKSTK.ispark[FLKSTK.npflka][mkbmx2-1] = itrack;
        FLKSTK.nparma++;
        FLKSTK.numpar[FLKSTK.npflka] = FLKSTK.nparma;
        FLKSTK.nevent[FLKSTK.npflka] = 0;
        FLKSTK.dfnear[FLKSTK.npflka] = +zerzer;
        
        /* Particle age (s)*/
        FLKSTK.agestk[FLKSTK.npflka] = +zerzer;
        FLKSTK.cmpath[FLKSTK.npflka] = +zerzer;
        FLKSTK.aknshr[FLKSTK.npflka] = -twotwo;

        /* Group number for "low" energy neutrons, set to 0 anyway*/
        FLKSTK.igroup[FLKSTK.npflka] = 0;
        
        /* Kinetic energy */
        Double_t p    = particle->P();
        Double_t mass = PAPROP.am[ifl + 6];
        FLKSTK.tkeflk[FLKSTK.npflka] =  TMath::Sqrt( p * p + mass * mass) - mass;
        /* Particle momentum*/
        FLKSTK.pmoflk [FLKSTK.npflka] = p;

        FLKSTK.txflk [FLKSTK.npflka] = cosx;
        FLKSTK.tyflk [FLKSTK.npflka] = cosy;
        FLKSTK.tzflk [FLKSTK.npflka] = cosz;
    
        /* Polarization cosines:*/
        if (polarisation.Mag()) {
            Double_t cospolx = polarisation.Px() / polarisation.Mag();
            Double_t cospoly = polarisation.Py() / polarisation.Mag();
            Double_t cospolz = sqrt(oneone - cospolx * cospolx - cospoly * cospoly);
            FLKSTK.txpol [FLKSTK.npflka] = cospolx;
            FLKSTK.typol [FLKSTK.npflka] = cospoly;
            FLKSTK.tzpol [FLKSTK.npflka] = cospolz;
        }
        else {
            FLKSTK.txpol [FLKSTK.npflka] = -twotwo;
            FLKSTK.typol [FLKSTK.npflka] = +zerzer;
            FLKSTK.tzpol [FLKSTK.npflka] = +zerzer;
        }

        /* Particle coordinates*/
        // Vertext coordinates;
        FLKSTK.xflk [FLKSTK.npflka] = particle->Vx();
        FLKSTK.yflk [FLKSTK.npflka] = particle->Vy();
        FLKSTK.zflk [FLKSTK.npflka] = particle->Vz();

        /*  Calculate the total kinetic energy of the primaries: don't change*/
        Int_t st_ilo =  FLKSTK.iloflk[FLKSTK.npflka];
        if ( st_ilo != 0 )
            SOURCM.tkesum +=
                ((FLKSTK.tkeflk[FLKSTK.npflka] + PAPROP.amdisc[st_ilo+6])
                 * FLKSTK.wtflk[FLKSTK.npflka]);
        else
            SOURCM.tkesum += (FLKSTK.tkeflk[FLKSTK.npflka] * FLKSTK.wtflk[FLKSTK.npflka]);

        /*  Here we ask for the region number of the hitting point.
         *     NRGFLK (LFLKSTK) = ...
         *  The following line makes the starting region search much more
         *  robust if particles are starting very close to a boundary:
         */
        geocrs( FLKSTK.txflk[FLKSTK.npflka],
                FLKSTK.tyflk[FLKSTK.npflka],
                FLKSTK.tzflk[FLKSTK.npflka] );
    
        Int_t idisc;

        georeg ( FLKSTK.xflk[FLKSTK.npflka],
                 FLKSTK.yflk[FLKSTK.npflka],
                 FLKSTK.zflk[FLKSTK.npflka],
                 FLKSTK.nrgflk[FLKSTK.npflka],
                 idisc);//<-- dummy return variable not used
        /*  Do not change these cards:*/
        Int_t igeohsm1 = 1;
        Int_t igeohsm2 = -11;
        geohsm ( FLKSTK.nhspnt[FLKSTK.npflka], igeohsm1, igeohsm2, LTCLCM.mlattc );
        FLKSTK.nlattc[FLKSTK.npflka] = LTCLCM.mlattc;
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
        OPPHST.nlatop [OPPHST.lstopp - 1] = LTCLCM.mlattc;
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
