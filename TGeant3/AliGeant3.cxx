/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/*
$Log$
Revision 1.4  2000/04/26 10:17:32  fca
Changes in Lego for G4 compatibility

Revision 1.3  2000/04/07 11:12:35  fca
G4 compatibility changes

Revision 1.2  2000/02/29 19:11:17  fca
Move gucode into AliGeant3.cxx

Revision 1.1  2000/02/23 16:25:25  fca
AliVMC and AliGeant3 classes introduced
ReadEuclid moved from AliRun to AliModule

*/

#include <TParticle.h>

#include "AliGeant3.h"
#include "AliRun.h"
#include "TGeant3.h"
#include "AliCallf77.h" 

#ifndef WIN32

# define rxgtrak rxgtrak_
# define rxstrak rxstrak_
# define rxkeep  rxkeep_ 
# define rxouth  rxouth_
#else

# define rxgtrak RXGTRAK 
# define rxstrak RXSTRAK 
# define rxkeep  RXKEEP  
# define rxouth  RXOUTH
#endif

ClassImp(AliGeant3)

AliGeant3::AliGeant3(const char *title) : 
  TGeant3(title) {}

void AliGeant3::FinishGeometry()
{
  TGeant3::FinishGeometry();
  //Create the color table
  SetColors();
}

//____________________________________________________________________________
void AliGeant3::Init()
{
  //
  //=================Create Materials and geometry
  TObjArray *modules = gAlice->Modules();
  TIter next(modules);
  AliModule *detector;
  while((detector = (AliModule*)next())) {
    // Initialise detector materials and geometry
    detector->CreateMaterials();
    detector->CreateGeometry();
    detector->BuildGeometry();
    detector->Init();
  }

  //Terminate building of geometry
  FinishGeometry();
}

//____________________________________________________________________________
void AliGeant3::ProcessRun(Int_t nevent)
{
  Int_t todo = TMath::Abs(nevent);
  for (Int_t i=0; i<todo; i++) {
  // Process one run (one run = one event)
     gAlice->BeginEvent();
     ProcessEvent();
     gAlice->FinishEvent();
  }
}

void AliGeant3::ProcessEvent()
{
  Gtrigi();
  Gtrigc();
  Gtrig();
}

//_____________________________________________________________________________
void AliGeant3::SetColors()
{
  //
  // Set the colors for all the volumes
  // this is done sequentially for all volumes
  // based on the number of their medium
  //
  Int_t kv, icol;
  Int_t jvolum=fGclink->jvolum;
  //Int_t jtmed=fGclink->jtmed;
  //Int_t jmate=fGclink->jmate;
  Int_t nvolum=fGcnum->nvolum;
  char name[5];
  //
  //    Now for all the volumes
  for(kv=1;kv<=nvolum;kv++) {
    //     Get the tracking medium
    Int_t itm=Int_t (fZq[fZlq[jvolum-kv]+4]);
    //     Get the material
    //Int_t ima=Int_t (fZq[fZlq[jtmed-itm]+6]);
    //     Get z
    //Float_t z=fZq[fZlq[jmate-ima]+7];
    //     Find color number
    //icol = Int_t(z)%6+2;
    //icol = 17+Int_t(z*150./92.);
    //icol = kv%6+2;
    icol = itm%6+2;
    strncpy(name,(char*)&fZiq[jvolum+kv],4);
    name[4]='\0';
    Gsatt(name,"COLO",icol);
  }
}

//_____________________________________________________________________________
//
//                 Interfaces to Fortran
//
//_____________________________________________________________________________

extern "C" void type_of_call  rxgtrak (Int_t &mtrack, Int_t &ipart, Float_t *pmom, 
				       Float_t &e, Float_t *vpos, Float_t *polar,
				       Float_t &tof)
{
  //
  //     Fetches next track from the ROOT stack for transport. Called by the
  //     modified version of GTREVE.
  //
  //              Track number in the ROOT stack. If MTRACK=0 no
  //      mtrack  more tracks are left in the stack to be
  //              transported.
  //      ipart   Particle code in the GEANT conventions.
  //      pmom[3] Particle momentum in GeV/c
  //      e       Particle energy in GeV
  //      vpos[3] Particle position
  //      tof     Particle time of flight in seconds
  //
  Int_t pdg;
  gAlice->GetNextTrack(mtrack, pdg, pmom, e, vpos, polar, tof);
  ipart = gMC->IdFromPDG(pdg);
  mtrack++;
}

//_____________________________________________________________________________
extern "C" void type_of_call 
#ifndef WIN32
rxstrak (Int_t &keep, Int_t &parent, Int_t &ipart, Float_t *pmom, 
	       Float_t *vpos, Float_t &tof, const char* cmech, Int_t &ntr, const int cmlen)
#else
rxstrak (Int_t &keep, Int_t &parent, Int_t &ipart, Float_t *pmom,
	 Float_t *vpos, Float_t &tof, const char* cmech, const int cmlen,
	 Int_t &ntr)
#endif
{
  //
  //     Fetches next track from the ROOT stack for transport. Called by GUKINE
  //     and GUSTEP.
  //
  //              Status of the track. If keep=0 the track is put
  //      keep    on the ROOT stack but it is not fetched for
  //              transport.
  //      parent  Parent track. If parent=0 the track is a primary.
  //              In GUSTEP the routine is normally called to store
  //              secondaries generated by the current track whose
  //              ROOT stack number is MTRACK (common SCKINE.
  //      ipart   Particle code in the GEANT conventions.
  //      pmom[3] Particle momentum in GeV/c
  //      vpos[3] Particle position
  //      tof     Particle time of flight in seconds
  //
  //      cmech   (CHARACTER*10) Particle origin. This field is user
  //              defined and it is not used inside the GALICE code.
  //      ntr     Number assigned to the particle in the ROOT stack.
  //
  char mecha[11];
  Float_t polar[3]={0.,0.,0.};
  for(int i=0; i<10 && i<cmlen; i++) mecha[i]=cmech[i];
  mecha[10]=0;
  Int_t pdg=gMC->PDGFromId(ipart);
  gAlice->SetTrack(keep, parent-1, pdg, pmom, vpos, polar, tof, mecha, ntr);
  ntr++;
}

//_____________________________________________________________________________
extern "C" void type_of_call  rxkeep(const Int_t &n)
{
  if( NULL==gAlice ) exit(1);
  
  if( n<=0 || n>gAlice->Particles()->GetEntries() )
    {
      printf("  Bad index n=%d must be 0<n<=%d\n",
	     n,gAlice->Particles()->GetEntries());
      exit(1);
    }
  
  ((TParticle*)(gAlice->Particles()->UncheckedAt(n-1)))->SetBit(Keep_Bit);
}

//_____________________________________________________________________________
extern "C" void type_of_call  rxouth ()
{
  //
  // Called by Gtreve at the end of each primary track
  //
  gAlice->FinishPrimary();
}


#ifndef WIN32
#  define gudigi gudigi_
#  define guhadr guhadr_
#  define guout  guout_
#  define guphad guphad_
#  define gudcay gudcay_
#  define guiget guiget_
#  define guinme guinme_
#  define guinti guinti_
#  define gunear gunear_
#  define guskip guskip_
#  define guview guview_
#  define gupara gupara_
#  define gudtim gudtim_
#  define guplsh guplsh_
#  define gutrev gutrev_
#  define gutrak gutrak_
#  define guswim guswim_
#  define gufld  gufld_
#  define gustep gustep_
#  define gukine gukine_
#  define uglast uglast_

#  define gheish gheish_
#  define flufin flufin_
#  define gfmfin gfmfin_
#  define gpghei gpghei_
#  define fldist fldist_
#  define gfmdis gfmdis_
#  define ghelx3 ghelx3_
#  define ghelix ghelix_
#  define grkuta grkuta_
#  define gtrack gtrack_
#  define gtreve_root gtreve_root_
#  define glast  glast_

#else
#  define gudigi GUDIGI
#  define guhadr GUHADR
#  define guout  GUOUT
#  define guphad GUPHAD
#  define gudcay GUDCAY
#  define guiget GUIGET
#  define guinme GUINME
#  define guinti GUINTI
#  define gunear GUNEAR
#  define guskip GUSKIP
#  define guview GUVIEW
#  define gupara GUPARA
#  define gudtim GUDTIM
#  define guplsh GUPLSH
#  define gutrev GUTREV
#  define gutrak GUTRAK
#  define guswim GUSWIM
#  define gufld  GUFLD
#  define gustep GUSTEP
#  define gukine GUKINE
#  define uglast UGLAST

#  define gheish GHEISH
#  define flufin FLUFIN
#  define gfmfin GFMFIN
#  define gpghei GPGHEI
#  define fldist FLDIST
#  define gfmdis GFMDIS
#  define ghelx3 GHELX3
#  define ghelix GHELIX
#  define grkuta GRKUTA
#  define gtrack GTRACK
#  define gtreve_root GTREVE_ROOT
#  define glast  GLAST

#endif

extern "C" type_of_call void gheish();
extern "C" type_of_call void flufin();
extern "C" type_of_call void gfmfin();
extern "C" type_of_call void gpghei();
extern "C" type_of_call void fldist();
extern "C" type_of_call void gfmdis();
extern "C" type_of_call void ghelx3(Float_t&, Float_t&, Float_t*, Float_t*);
extern "C" type_of_call void ghelix(Float_t&, Float_t&, Float_t*, Float_t*);
extern "C" type_of_call void grkuta(Float_t&, Float_t&, Float_t*, Float_t*);
extern "C" type_of_call void gtrack();
extern "C" type_of_call void gtreve_root();
extern "C" type_of_call void glast();

extern "C" type_of_call {

//______________________________________________________________________
void gudigi() 
{
//
//    ******************************************************************
//    *                                                                *
//    *       User routine to digitize one event                       *
//    *                                                                *
//    *    ==>Called by : GTRIG                                        *
//    *                                                                *
//    ******************************************************************

}


//______________________________________________________________________
void guhadr()
{
//
//    ******************************************************************
//    *                                                                *
//    *       User routine to generate one hadronic interaction        *
//    *                                                                *
//    *    ==>Called by : GTHADR,GTNEUT                                *
//    *                                                                *
//    ******************************************************************
//
//
//    ------------------------------------------------------------------
//
      TGeant3* geant3 = (TGeant3*) gMC;
      Int_t ihadr=geant3->Gcphys()->ihadr;
      if (ihadr<4)       gheish();
      else if (ihadr==4) flufin();
      else               gfmfin();
}

//______________________________________________________________________
void guout()
{
//
//    ******************************************************************
//    *                                                                *
//    *       User routine called at the end of each event             *
//    *                                                                *
//    *    ==>Called by : GTRIG                                        *
//    *                                                                *
//    ******************************************************************
//
//
//    ------------------------------------------------------------------
//
}

//______________________________________________________________________
void guphad()
{
//
//    ******************************************************************
//    *                                                                *
//    *       User routine to compute Hadron. inter. probabilities     *
//    *                                                                *
//    *    ==>Called by : GTHADR,GTNEUT                                *
//    *                                                                *
//    ******************************************************************
//
//
//    ------------------------------------------------------------------
//
      TGeant3* geant3 = (TGeant3*) gMC;
      Int_t ihadr=geant3->Gcphys()->ihadr;
      if (ihadr<4)       gpghei();
      else if (ihadr==4) fldist();
      else               gfmdis();
}

//______________________________________________________________________
void gudcay()
{
//
//    ******************************************************************
//    *                                                                *
//    *       User routine to decay particles                          *
//    *                                                                *
//    *    ==>Called by : GDECAY                                       *
//    *                                                                *
//    ******************************************************************
//
//
//    ------------------------------------------------------------------
//
}

//______________________________________________________________________
void guiget(Int_t&, Int_t&, Int_t&)
{
//
//    ******************************************************************
//    *                                                                *
//    *       User routine for interactive control of GEANT            *
//    *                                                                *
//    *    ==>Called by : <GXINT>, GINCOM                              *
//    *                                                                *
//    ******************************************************************
//
//
//    ------------------------------------------------------------------
//
}

//______________________________________________________________________
void guinme(Float_t*, Int_t&, Float_t*, Int_t& IYES)
{
//
//    **********************************************
//    *                                            *
//    *    USER ROUTINE TO PROVIDE GINME FUNCTION  *
//    *    FOR ALL USER SHAPES IDENTIFIED BY THE   *
//    *    SHAPE NUMBER SH. POINT IS GIVEN IN X    *
//    *    THE PARAMETERS ARE GIVEN IN P. IYES IS  *
//    *    RETURNED 1 IF POINT IS IN, 0 IF POINT   *
//    *    IS OUT AND LESS THAN ZERO IF SHAPE      *
//    *    NUMBER IS NOT SUPPORTED.                *
//    *                                            *
//    *    ==>Called by : GINME                    *
//    *                                            *
//    **********************************************
//
      IYES=-1;
}

//______________________________________________________________________
void guinti()
{
//
//    ******************************************************************
//    *                                                                *
//    *       User routine for interactive version                     *
//    *                                                                *
//    *    ==>Called by : <GXINT>,  GINTRI                             *
//    *                                                                *
//    ******************************************************************
//
//
//    ------------------------------------------------------------------
//
}

//______________________________________________________________________
void gunear(Int_t&, Int_t&, Float_t*, Int_t&)
{
//
//    ******************************************************************
//    *                                                                *
//    *    User search                                                 *
//    *       ISEARC to identify the given volume                      *
//    *       ICALL  to identify the calling routine                   *
//    *              1 GMEDIA like                                     *
//    *              2 GNEXT like                                      *
//    *       X      coordinates (+direction for ICALL=2)              *
//    *       JNEAR  address of default list of neighbours             *
//    *              (list to be overwriten by user)                   *
//    *                                                                *
//    *    Called by : GFTRAC, GINVOL, GTMEDI, GTNEXT, GNEXT, GMEDIA   *
//    *                                                                *
//    ******************************************************************
//
//
//    ------------------------------------------------------------------
//
}

//______________________________________________________________________
void guskip(Int_t& ISKIP)
{
//
//    ******************************************************************
//    *                                                                *
//    *   User routine to skip unwanted tracks                         *
//    *                                                                *
//    *   Called by : GSSTAK                                           *
//    *   Author    : F.Bruyant                                        *
//    *                                                                *
//    ******************************************************************
//
//
//    ------------------------------------------------------------------
//
      ISKIP = 0;
}

//______________________________________________________________________
void guswim(Float_t& CHARGE, Float_t& STEP, Float_t* VECT, Float_t* VOUT)
{
//
//    ******************************************************************
//    *                                                                *
//    *       User routine to control tracking of one track            *
//    *       in a magnetic field                                      *
//    *                                                                *
//    *    ==>Called by : GTELEC,GTHADR,GTMUON                         *
//    *                                                                *
//    ******************************************************************
//
//
//    ------------------------------------------------------------------
//
  TGeant3* geant3 = (TGeant3*) gMC;
  Int_t ifield=geant3->Gctmed()->ifield;
  Float_t fieldm=geant3->Gctmed()->fieldm;

  if (ifield==3) {
    Float_t fldcharge = fieldm*CHARGE;
    ghelx3(fldcharge,STEP,VECT,VOUT);
  }
  else if (ifield==2) ghelix(CHARGE,STEP,VECT,VOUT);
  else                grkuta(CHARGE,STEP,VECT,VOUT);
}

//______________________________________________________________________
void guview(Int_t&, Int_t&, DEFCHARD, Int_t& DEFCHARL)
{
//
//    ******************************************************************
//    *                                                                *
//    *       User routine for interactive version                     *
//    *                                                                *
//    *    ==>Called by : <GXINT>, GINC1                               *
//    *                                                                *
//    ******************************************************************
//
//
//    ------------------------------------------------------------------
//
}

//______________________________________________________________________
void gupara()
{
//
//    ******************************************************************
//    *                                                                *
//    *       User routine called every time a particle falls below    *
//    *       parametrization threshold. This routine should create    *
//    *       the parametrization stack, and, when this is full,       *
//    *       parametrize the shower and track the geantinos.          *
//    *                                                                *
//    *    ==>Called by : GTRACK                                       *
//    *                                                                *
//    ******************************************************************
//
//
//    ------------------------------------------------------------------
//
}

//______________________________________________________________________
Float_t gudtim(Float_t&, Float_t&, Int_t&, Int_t&)
{
//
//    ******************************************************************
//    *                                                                *
//    *       User function called by GCDRIF to return drift time      *
//    *                                                                *
//    *    ==>Called by : GCDRIF                                       *
//    *                                                                *
//    ******************************************************************
//
//
//    ------------------------------------------------------------------
//
      return 0;
}


//______________________________________________________________________
Float_t guplsh(Int_t&, Int_t&)
{
//
//    ******************************************************************
//    *                                                                *
//    *                                                                *
//    *    ==>Called by : GLISUR                                       *
//    *                                                                *
//    ******************************************************************
//
//
//    ------------------------------------------------------------------
//
//
//*** By default this defines perfect smoothness
      return 1;
}

//______________________________________________________________________
void gutrak()
{
//
//    ******************************************************************
//    *                                                                *
//    *       User routine to control tracking of one track            *
//    *                                                                *
//    *    ==>Called by : GTREVE                                       *
//    *                                                                *
//    ******************************************************************
//
//
//    ------------------------------------------------------------------
//
     Int_t ndet = gAlice->Modules()->GetLast();
     TObjArray &dets = *gAlice->Modules();
     AliModule *module;
     Int_t i;

     for(i=0; i<=ndet; i++)
       if((module = (AliModule*)dets[i]))
	 module->PreTrack();

     gtrack();

     for(i=0; i<=ndet; i++)
       if((module = (AliModule*)dets[i]))
	 module->PostTrack();
}

//______________________________________________________________________
void gutrev()
{
//
//    ******************************************************************
//    *                                                                *
//    *       User routine to control tracking of one event            *
//    *                                                                *
//    *    ==>Called by : GTRIG                                        *
//    *                                                                *
//    ******************************************************************
//
//
//    ------------------------------------------------------------------
//
  gtreve_root();
}


//______________________________________________________________________
void gufld(Float_t *x, Float_t *b)
{
      if(gAlice->Field()) {
         gAlice->Field()->Field(x,b);
      } else {
         printf("No mag field defined!\n");
         b[0]=b[1]=b[2]=0;
      }
}

//______________________________________________________________________
void gustep()
{
//
//    ******************************************************************
//    *                                                                *
//    *       User routine called at the end of each tracking step     *
//    *       INWVOL is different from 0 when the track has reached    *
//    *              a volume boundary                                 *
//    *       ISTOP is different from 0 if the track has stopped       *
//    *                                                                *
//    *    ==>Called by : GTRACK                                       *
//    *                                                                *
//    ******************************************************************
//


  TLorentzVector x;
  Float_t r;
  Int_t ipp, jk, id, nt;
  Float_t polar[3]={0,0,0};
  Float_t mom[3];
  const char *chproc;
  
  TGeant3* geant3 = (TGeant3*) gMC;

  //     Stop particle if outside user defined tracking region 
  gMC->TrackPosition(x);
  r=TMath::Sqrt(x[0]*x[0]+x[1]*x[1]);
  if (r > gAlice->TrackingRmax() || TMath::Abs(x[2]) > gAlice->TrackingZmax()) {
	gMC->StopTrack();
  }

  // --- Add new created particles 
  if (gMC->NSecondaries() > 0) {
    chproc=gMC->ProdProcess();
    for (jk = 0; jk < geant3->Gcking()->ngkine; ++jk) {
      ipp = Int_t (geant3->Gcking()->gkin[jk][4]+0.5);
      // --- Skip neutrinos! 
      if (ipp != 4) {
	gAlice->SetTrack(1,gAlice->CurrentTrack(),gMC->PDGFromId(ipp), geant3->Gcking()->gkin[jk], 
			 geant3->Gckin3()->gpos[jk], polar,geant3->Gctrak()->tofg, chproc, nt);
      }
    }
  }
  // Cherenkov photons here
  if ( geant3->Gckin2()->ngphot ) {
    for (jk = 0; jk < geant3->Gckin2()->ngphot; ++jk) {
      mom[0]=geant3->Gckin2()->xphot[jk][3]*geant3->Gckin2()->xphot[jk][6];
      mom[1]=geant3->Gckin2()->xphot[jk][4]*geant3->Gckin2()->xphot[jk][6];
      mom[2]=geant3->Gckin2()->xphot[jk][5]*geant3->Gckin2()->xphot[jk][6];
      gAlice->SetTrack(1, gAlice->CurrentTrack(), gMC->PDGFromId(50),
		       mom,                             //momentum
		       geant3->Gckin2()->xphot[jk],     //position
		       &geant3->Gckin2()->xphot[jk][7], //polarisation
		       geant3->Gckin2()->xphot[jk][10], //time of flight
		       "Cherenkov", nt);
      }
  }
  // --- Particle leaving the setup ?
  if (!gMC->IsTrackOut()) 
    if ((id=gAlice->DetFromMate(geant3->Gctmed()->numed)) >= 0) gAlice->StepManager(id);

  // --- Standard GEANT debug routine 
  if(geant3->Gcflag()->idebug) geant3->Gdebug();
}

//______________________________________________________________________
void gukine ()
{
//
//    ******************************************************************
//    *                                                                *
//    *       Read or Generates Kinematics for primary tracks          *
//    *                                                                *
//    *    ==>Called by : GTRIG                                        *
//    *                                                                *
//    ******************************************************************
//
//
//    ------------------------------------------------------------------
//
  gAlice->Generator()->Generate();
}


//______________________________________________________________________
void uglast()
{
//
//    ******************************************************************
//    *                                                                *
//    *       User routine called at the end of the run                *
//    *                                                                *
//    *    ==>Called by : GRUN                                         *
//    *                                                                *
//    ******************************************************************
//
//
}
}

