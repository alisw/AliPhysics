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
Revision 1.8  2001/01/26 22:05:41  morsch
Unresolved conflicts resolved.

Revision 1.7  2001/01/26 21:54:46  morsch
Use access functions to AliMUONHit member data.


Revision 1.6  2001/01/26 20:00:53  hristov
Major upgrade of AliRoot code

Revision 1.4  2000/12/21 22:14:38  morsch
Clean-up of coding rule violations.

Revision 1.3  2000/10/02 21:28:09  fca
Removal of useless dependecies via forward declarations

Revision 1.2  2000/06/15 07:58:49  morsch
Code from MUON-dev joined

Revision 1.1.2.7  2000/06/09 22:06:29  morsch
Some coding rule violations corrected. Will soon be obsolete.

Revision 1.1.2.6  2000/05/02 07:15:29  morsch
Put back TH1.h and TH2.h includes.

Revision 1.1.2.5  2000/02/17 18:12:43  morsch
Corrections in trackf_read_spoint causing segmentation violations in previous version (I. Chevrot)
New histos (I. Chevrot)

Revision 1.1.2.4  2000/02/15 18:01:08  morsch
Log messages

Revision 1.1.2.3  2000/02/15 17:59:01  morsch
Log message added

Revision 1.1.2.2  2000/02/15 18:54:56  morsch
Reference between track contributing to reconstructed hit and particle corrected  
*/

#include "AliCallf77.h" 
#include "AliMUONTrackReconstructor.h" 
#include "AliRun.h"
#include "AliMUON.h"
#include "AliMC.h"

#include "AliMUONHit.h"
#include "AliMUONPadHit.h"
#include "AliMUONDigit.h"
#include "AliMUONRawCluster.h"
#include "AliMUONReconstHit.h"

#include "AliPDG.h"

#include <TRandom.h> 
#include <TFile.h> 
#include <TH1.h>
#include <TH2.h>  
#include <TTree.h> 
#include <TParticle.h> 
#include <TMinuit.h>
#include <iostream.h>

#ifndef WIN32 
# define reco_init       reco_init_
# define cutpxz          cutpxz_
# define sigmacut        sigmacut_
# define xpreci          xpreci_
# define ypreci          ypreci_
# define reconstmuon     reconstmuon_
# define reconstmuon2    reconstmuon2_
# define trackf_read_geant     trackf_read_geant_
# define trackf_read_fit     trackf_read_fit_
# define trackf_read_spoint     trackf_read_spoint_
# define chfill          chfill_
# define chfill2         chfill2_
# define chf1            chf1_
# define chfnt           chfnt_
# define hist_create     hist_create_
# define hist_closed     hist_closed_
# define rndm            rndm_
# define fcn             fcn_
# define trackf_fit      trackf_fit_
# define prec_fit        prec_fit_
# define fcnfit          fcnfit_
# define reco_term       reco_term_
#else 
# define reco_init       RECO_INIT
# define cutpxz          CUTPXZ
# define sigmacut        SIGMACUT
# define xpreci          XPRECI
# define ypreci          YPRECI
# define reconstmuon     RECONSTMUON
# define reconstmuon2    RECONSTMUON2
# define trackf_read_geant     TRACKF_READ_GEANT
# define trackf_read_fit     TRACKF_READ_FIT
# define trackf_read_spoint     TRACKF_READ_SPOINT
# define chfill          CHFILL
# define chfill2         CHFILL2
# define chf1            CHF1
# define chfnt           CHFNT
# define hist_create     HIST_CREATE
# define hist_closed     HIST_CLOSED
# define rndm            RNDM
# define fcn             FCN
# define trackf_fit      TRACKF_FIT
# define prec_fit        PREC_FIT
# define fcnfit          FCNFIT
# define reco_term       RECO_TERM
#endif 

extern "C"
{
void type_of_call reco_init(Double_t &, Double_t &, Double_t &);
void type_of_call reco_term();
void type_of_call cutpxz(Double_t &);
void type_of_call sigmacut(Double_t &);
void type_of_call xpreci(Double_t &);
void type_of_call ypreci(Double_t &);
void type_of_call reconstmuon(Int_t &, Int_t &, Int_t &, Int_t &, Int_t &);
void type_of_call reconstmuon2(Int_t &, Int_t &, Int_t &);
void type_of_call trackf_read_fit(Int_t &, Int_t &, Int_t &, Int_t *, Double_t *, Double_t *);
void type_of_call trackf_read_geant(Int_t *, Double_t *, Double_t *, Double_t *, Int_t *, Int_t *, Double_t *, Double_t *, Double_t *, Double_t *,Int_t &, Double_t *, Double_t *, Double_t *, Int_t &, Int_t &, Double_t *, Double_t *, Double_t *, Double_t *); 
void type_of_call trackf_read_spoint(Int_t *, Double_t *, Double_t *, Double_t *, Int_t *, Int_t *, Double_t *, Double_t *, Double_t *, Double_t *,Int_t &, Double_t *, Double_t *, Double_t *, Int_t &, Int_t &, Double_t *, Double_t *, Double_t *, Double_t *); 
void type_of_call chfill(Int_t &, Float_t &, Float_t &, Float_t &);
void type_of_call chfill2(Int_t &, Float_t &, Float_t &, Float_t &);
void type_of_call chf1(Int_t &, Float_t &, Float_t &);
void type_of_call chfnt(Int_t &, Int_t &, Int_t *, Int_t *, Float_t *, Float_t *, Float_t *, Float_t *, Float_t *, Float_t *, Float_t *, Float_t *);
void type_of_call hist_create();
void type_of_call hist_closed();
void type_of_call fcnf(Int_t &, Double_t *, Double_t &, Double_t *, Int_t);
void type_of_call fcn(Int_t &, Double_t *, Double_t &, Double_t *, Int_t &, Int_t &);
void type_of_call trackf_fit(Int_t &, Double_t *, Double_t *, Double_t &, Double_t &, Double_t &, Double_t &, Double_t &);
void type_of_call prec_fit(Double_t &, Double_t &, Double_t &, Double_t &, Double_t&, Double_t &, Double_t &, Double_t &, Double_t &, Double_t &, Double_t &, Double_t &, Double_t &, Double_t &, Double_t &);
void type_of_call fcnfitf(Int_t &, Double_t *, Double_t &, Double_t *, Int_t);
void type_of_call fcnfit(Int_t &, Double_t *, Double_t &, Double_t *, Int_t &, Int_t &);
Float_t type_of_call rndm() {return gRandom->Rndm();}
void type_of_call fit_trace(Float_t &, Float_t &, Float_t &, Float_t &, Float_t &, Float_t &, Int_t &, Int_t &, Int_t &, Int_t &);
}

static TTree *gAliNtupleGlobal;
static TFile *gAliFileGlobal;
static TTree *gAliTreeK1;
static TTree *gAliTrH1;
static TClonesArray *gAliHits2;        //List of hits for one track only
static TClonesArray *gAliParticles2;   //List of particles in the Kine tree


// variables of the tracking ntuple 
struct { 
  Int_t ievr;           // number of event 
  Int_t ntrackr;        // number of tracks per event
  Int_t istatr[500];    // 1 = good muon, 2 = ghost, 0 = something else
  Int_t isignr[500];    // sign of the track
  Float_t pxr[500];     // x momentum of the reconstructed track
  Float_t pyr[500];     // y momentum of the reconstructed track
  Float_t pzr[500];     // z momentum of the reconstructed track
  Float_t zvr[500];     // z vertex 
  Float_t chi2r[500];   // chi2 of the fit of the track with the field map
  Float_t pxv[500];     // x momentum at vertex
  Float_t pyv[500];     // y momentum at vertex
  Float_t pzv[500];     // z momentum at vertex
} NtupleSt;

ClassImp(AliMUONTrackReconstructor)

//___________________________________________________
AliMUONTrackReconstructor::AliMUONTrackReconstructor()
{
// Constructor
//
   fSPxzCut   = 3.0;
   fSSigmaCut = 4.0;
   fSXPrec    = 0.01; 
   fSYPrec    = 0.144;
}

//_____________________________________________________________________________
void AliMUONTrackReconstructor::Reconst2(Int_t &ifit, Int_t &idebug, Int_t &nev)
{
//
//
    reconstmuon2(ifit,idebug,nev);
}

//_____________________________________________________________________________
void AliMUONTrackReconstructor::Reconst(Int_t &ifit, Int_t &idebug, Int_t bgdEvent, Int_t &nev, Int_t &idres, Int_t &ireadgeant, Option_t *option,Text_t *filename)
{
  //
  // open kine and hits tree of background file for reconstruction of geant hits 
  // call tracking fortran program
  static Bool_t first=kTRUE;
  static TFile *pFile;
  const char *addBackground = strstr(option,"Add");
  
  if (addBackground ) { // only in case of background with geant hits 
    if(first) {
      fFileName=filename;
      cout<<"filename  "<<fFileName<<endl;
      pFile=new TFile(fFileName);
      cout<<"I have opened "<<fFileName<<" file "<<endl;
      gAliHits2     = new TClonesArray("AliMUONHit",1000);
      gAliParticles2 = new TClonesArray("TParticle",1000);
      first=kFALSE;
    }
    pFile->cd();
    if(gAliHits2) gAliHits2->Clear();
    if(gAliParticles2) gAliParticles2->Clear();
    if(gAliTrH1) delete gAliTrH1;
    gAliTrH1=0;
    if(gAliTreeK1) delete gAliTreeK1;
    gAliTreeK1=0;
    // Get Hits Tree header from file
    char treeName[20];
    sprintf(treeName,"TreeH%d",bgdEvent);
    gAliTrH1 = (TTree*)gDirectory->Get(treeName);
    //     printf("gAliTrH1 %p of treename %s for event %d \n",gAliTrH1,treeName,bgdEvent);
    if (!gAliTrH1) {
      printf("ERROR: cannot find Hits Tree for event:%d\n",bgdEvent);
    }
    // set branch addresses
    TBranch *branch;
    char branchname[30];
    AliMUON *pMUON  = (AliMUON*) gAlice->GetModule("MUON");  
    sprintf(branchname,"%s",pMUON->GetName());       
    if (gAliTrH1 && gAliHits2) {
      branch = gAliTrH1->GetBranch(branchname);
      if (branch) branch->SetAddress(&gAliHits2);
    }
    gAliTrH1->GetEntries();
    // get the Kine tree
    sprintf(treeName,"TreeK%d",bgdEvent);
    gAliTreeK1 = (TTree*)gDirectory->Get(treeName);
    if (!gAliTreeK1) {
      printf("ERROR: cannot find Kine Tree for event:%d\n",bgdEvent);
    }
    // set branch addresses
    if (gAliTreeK1) 
      gAliTreeK1->SetBranchAddress("Particles", &gAliParticles2);
    gAliTreeK1->GetEvent(0);
    
    // get back to the first file
    TTree *treeK = gAlice->TreeK();
    TFile *file1 = 0;
    if (treeK) file1 = treeK->GetCurrentFile();
    file1->cd();
    
  } // end if addBackground
  
  // call tracking fortran program
  reconstmuon(ifit,idebug,nev,idres,ireadgeant);
}

//________________________________________________________________________________
void AliMUONTrackReconstructor::Init(Double_t &seff, Double_t &sb0, Double_t &sbl3)
{
  //
  // introduce in fortran program somme parameters and cuts for tracking 
  // create output file "reconst.root" (histos + ntuple)
  cutpxz(fSPxzCut);          // Pxz cut (GeV/c) to begin the track finding
  sigmacut(fSSigmaCut);      // Number of sigmas delimiting the searching areas
  xpreci(fSXPrec);           // Chamber precision in X (cm) 
  ypreci(fSYPrec);           // Chamber precision in Y (cm)
  reco_init(seff,sb0,sbl3);
}

//__________________________________________
void AliMUONTrackReconstructor::FinishEvent()
{
   // Finish
   // TTree *treeK = gAlice->TreeK();
   // TFile *file1 = 0;
   // if (treeK) file1 = treeK->GetCurrentFile();
   // file1->cd();
}

//_____________________________________
void AliMUONTrackReconstructor::Close()
{
  //
  // write histos and ntuple to "reconst.root" file
    reco_term();
}

//________________________________________________________
void chfill(Int_t &id, Float_t &x, Float_t &y, Float_t &w)
{
  //
  // fill histo like hfill in fortran
    char name[5];
    sprintf(name,"h%d",id);
    TH1F *h1 = (TH1F*) gDirectory->Get(name);
    h1->Fill(x);
}

//_________________________________________________________
void chfill2(Int_t &id, Float_t &x, Float_t &y, Float_t &w)
{
  //
  // fill histo like hfill2 in fortran
    char name[5];
    sprintf(name,"h%d",id);
    TH2F *h2 = (TH2F*) gDirectory->Get(name);
    h2->Fill(x,y,w);
}

//__________________________________________
void chf1(Int_t &id, Float_t &x, Float_t &w)
{
  //
  // fill histo like hf1 in fortran
    char name[5];
    sprintf(name,"h%d",id);
    TH1F *h1 = (TH1F*) gDirectory->Get(name);
    h1->Fill(x,w);
}

//_________________
void hist_create()
{
  //
  // Create an output file ("reconst.root")
  // Create some histograms and an ntuple

    gAliFileGlobal = new TFile("reconst.root","RECREATE","Ntuple - reconstruction");

   gAliNtupleGlobal = new TTree("ntuple","Reconst ntuple");
   gAliNtupleGlobal->Branch("ievr",&NtupleSt.ievr,"ievr/I");
   gAliNtupleGlobal->Branch("ntrackr",&NtupleSt.ntrackr,"ntrackr/I");
   gAliNtupleGlobal->Branch("istatr",&NtupleSt.istatr[0],"istatr[500]/I");
   gAliNtupleGlobal->Branch("isignr",&NtupleSt.isignr[0],"isignr[500]/I");
   gAliNtupleGlobal->Branch("pxr",&NtupleSt.pxr[0],"pxr[500]/F");
   gAliNtupleGlobal->Branch("pyr",&NtupleSt.pyr[0],"pyr[500]/F");
   gAliNtupleGlobal->Branch("pzr",&NtupleSt.pzr[0],"pzr[500]/F");
   gAliNtupleGlobal->Branch("zvr",&NtupleSt.zvr[0],"zvr[500]/F");
   gAliNtupleGlobal->Branch("chi2r",&NtupleSt.chi2r[0],"chi2r[500]/F");
   gAliNtupleGlobal->Branch("pxv",&NtupleSt.pxv[0],"pxv[500]/F");
   gAliNtupleGlobal->Branch("pyv",&NtupleSt.pyv[0],"pyv[500]/F");
   gAliNtupleGlobal->Branch("pzv",&NtupleSt.pzv[0],"pzv[500]/F");

   // test aliroot

  new TH1F("h100","particule id du hit geant",20,0.,20.);
  new TH1F("h101","position en x du hit geant",100,-200.,200.);
  new TH1F("h102","position en y du hit geant",100,-200.,200.);
  new TH1F("h103","chambre de tracking concernee",15,0.,14.);
  new TH1F("h104","moment ptot du hit geant",50,0.,100.);
  new TH1F("h105","px au vertex",50,0.,20.);
  new TH1F("h106","py au vertex",50,0.,20.);
  new TH1F("h107","pz au vertex",50,0.,20.);
  new TH1F("h108","position zv",50,-15.,15.);
  new TH1F("h109","position en x du hit reconstruit",100,-300.,300.);
  new TH1F("h110","position en y du hit reconstruit",100,-300.,300.);
  new TH1F("h111","delta x station 1",100,-0.3,0.3);
  new TH1F("h112","delta x station 2",100,-0.3,0.3);
  new TH1F("h113","delta x station 3",100,-0.3,0.3);
  new TH1F("h114","delta x station 4",100,-0.5,0.5);
  new TH1F("h115","delta x station 5",100,-0.5,0.5);
  new TH1F("h116","delta x station 1",100,-2,2);
  new TH1F("h117","delta x station 2",100,-2,2);
  new TH1F("h121","delta y station 1",100,-0.04,0.04);
  new TH1F("h122","delta y station 2",100,-0.04,0.04);
  new TH1F("h123","delta y station 3",100,-0.04,0.04);
  new TH1F("h124","delta y station 4",100,-0.04,0.04);
  new TH1F("h125","delta y station 5",100,-0.04,0.04);

  /*  char hname[30];
      char hname1[30];
      for (int i=0;i<10;i++) {
      sprintf(hname,"deltax%d",i);
      sprintf(hname1,"h12%d",i);
      new TH1F(hname1,hname ,100,-0.4,0.4);
      sprintf(hname,"deltay%d",i);
      sprintf(hname1,"h13%d",i);
      new TH1F(hname1,hname ,100,-0.4,0.4);
      }
  */
  new TH2F("h2000","VAR X st. 5",30,3.0,183.0,100,0.,25.);
  new TH2F("h2001","VAR Y st. 5",30,3.0,183.0,100,0.,25.);

  new TH2F("h2500","P vs X HHIT",30,3.0,183.0,200,0.,200.);
  new TH2F("h2501","P vs X HHIT**2",30,3.0,183.0,200,0.,5000.);
  new TH2F("h2502","P vs X EPH2 st. 5",30,3.0,183.0,100,0.,0.000005);
  new TH2F("h2503","P vs X EAL2 st. 5",30,3.0,183.0,100,0.,0.01);
  //new TH2F("h2504","P vs X EXM2 st. 5",30,3.0,183.0,100,0.,1.5);
  new TH2F("h2504","P vs X EXM2 st. 5",30,3.0,183.0,100,0.,0.1);
  new TH2F("h2505","P vs X EYM2 st. 5",30,3.0,183.0,100,0.,30.);

  new TH2F("h2507","P vs X EPH st. 5",30,3.0,183.0,100,0.,0.003);
  new TH2F("h2508","P vs X EAL st. 5",30,3.0,183.0,100,0.,0.3);
  //new TH2F("h2509","P vs X EXM st. 5",30,3.0,183.0,100,0.,1.5);
  new TH2F("h2509","P vs X EXM st. 5",30,3.0,183.0,100,0.,0.4);
  new TH2F("h2510","P vs X EYM st. 5",30,3.0,183.0,100,0.,30.);

  new TH2F("h2511","P vs X EPH cut st. 5",30,3.0,183.0,100,0.,0.01);
  new TH2F("h2512","P vs X EAL cut st. 5",30,3.0,183.0,100,0.,0.3);
  //new TH2F("h2513","P vs X EXM cut st. 5",30,3.0,183.0,100,0.,1.5);
  new TH2F("h2513","P vs X EXM cut st. 5",30,3.0,183.0,100,0.,0.4);
  new TH2F("h2514","P vs X EYM cut st. 5",30,3.0,183.0,100,0.,30.);
  // 4
  new TH2F("h2400","P vs X HHIT",30,3.0,183.0,200,0.,200.);
  new TH2F("h2401","P vs X HHIT**2",30,3.0,183.0,200,0.,5000.);
  new TH2F("h2402","P vs X EPH2 st. 4",30,3.0,183.0,100,0.,0.000005);
  new TH2F("h2403","P vs X EAL2 st. 4",30,3.0,183.0,100,0.,0.05);
  //new TH2F("h2404","P vs X EXM2 st. 4",30,3.0,183.0,100,0.,1.5);
  new TH2F("h2404","P vs X EXM2 st. 4",30,3.0,183.0,100,0.,0.1);
  new TH2F("h2405","P vs X EYM2 st. 4",30,3.0,183.0,100,0.,30.);

  new TH2F("h2407","P vs X EPH st. 4",30,3.0,183.0,100,0.,0.003);
  new TH2F("h2408","P vs X EAL st. 4",30,3.0,183.0,100,0.,0.3);
  //new TH2F("h2409","P vs X EXM st. 4",30,3.0,183.0,100,0.,1.5);
  new TH2F("h2409","P vs X EXM st. 4",30,3.0,183.0,100,0.,0.1);
  new TH2F("h2410","P vs X EYM st. 4",30,3.0,183.0,100,0.,30.);

  new TH2F("h2411","P vs X EPH cut st. 4",30,3.0,183.0,100,0.,0.01);
  new TH2F("h2412","P vs X EAL cut st. 4",30,3.0,183.0,100,0.,0.3);
  //new TH2F("h2413","P vs X EXM cut st. 4",30,3.0,183.0,100,0.,1.5);
  new TH2F("h2413","P vs X EXM cut st. 4",30,3.0,183.0,100,0.,0.1);
  new TH2F("h2414","P vs X EYM cut st. 4",30,3.0,183.0,100,0.,30.);
  // 3
  new TH1F("h2301","P2",30,3.0,183.0);
  new TH2F("h2302","P2 vs X EPH2 st. 3",30,3.0,183.0,100,0.,0.0006);
  new TH2F("h2303","P2 vs X EAL2 st. 3",30,3.0,183.0,100,0.,0.0005);
  //new TH2F("h2304","P2 vs X EXM2 st. 3",30,3.0,183.0,100,0.,1.5);
  new TH2F("h2304","P2 vs X EXM2 st. 3",30,3.0,183.0,100,0.,2.);
  new TH2F("h2305","P2 vs X EYM2 st. 3",30,3.0,183.0,100,0.,3.);

  new TH2F("h2307","P vs X EPH2 st. 3",30,3.0,183.0,100,0.,0.0006);
  new TH2F("h2308","P vs X EAL2 st. 3",30,3.0,183.0,100,0.,0.005);
  //new TH2F("h2309","P vs X EXM2 st. 3",30,3.0,183.0,100,0.,1.5);
  new TH2F("h2309","P vs X EXM2 st. 3",30,3.0,183.0,100,0.,2.);
  new TH2F("h2310","P vs X EYM2 st. 3",30,3.0,183.0,100,0.,3.);

  new TH2F("h2311","P vs X EPH cut st. 3",30,3.0,183.0,100,0.,0.06);
  new TH2F("h2312","P vs X EAL cut st. 3",30,3.0,183.0,100,0.,0.05);
  //new TH2F("h2313","P vs X EXM cut st. 3",30,3.0,183.0,100,0.,1.5);
  new TH2F("h2313","P vs X EXM cut st. 3",30,3.0,183.0,100,0.,6.);
  new TH2F("h2314","P vs X EYM cut st. 3",30,3.0,183.0,100,0.,7.);

  new TH2F("h2315","P2 vs X EPH cut st. 3",30,3.0,183.0,100,0.,0.06);
  new TH2F("h2316","P2 vs X EAL cut st. 3",30,3.0,183.0,100,0.,0.05);
  //new TH2F("h2317","P2 vs X EXM cut st. 3",30,3.0,183.0,100,0.,1.5);
  new TH2F("h2317","P2 vs X EXM cut st. 3",30,3.0,183.0,100,0.,6.);
  new TH2F("h2318","P2 vs X EYM cut st. 3",30,3.0,183.0,100,0.,7.);
  
  // 2
  new TH1F("h2201","P2",30,3.0,183.0);
  new TH2F("h2202","P2 vs X EPH2 st. 2",30,3.0,183.0,100,0.,0.0006);
  new TH2F("h2203","P2 vs X EAL2 st. 2",30,3.0,183.0,100,0.,0.005);
  //new TH2F("h2204","P2 vs X EXM2 st. 2",30,3.0,183.0,100,0.,1.5);
  new TH2F("h2204","P2 vs X EXM2 st. 2",30,3.0,183.0,100,0.,7.);
  new TH2F("h2205","P2 vs X EYM2 st. 2",30,3.0,183.0,100,0.,5.);

  new TH2F("h2207","P vs X EPH2 st. 2",30,3.0,183.0,100,0.,0.0006);
  new TH2F("h2208","P vs X EAL2 st. 2",30,3.0,183.0,100,0.,0.005);
  //new TH2F("h2209","P vs X EXM2 st. 2",30,3.0,183.0,100,0.,1.5);
  new TH2F("h2209","P vs X EXM2 st. 2",30,3.0,183.0,100,0.,7.);
  new TH2F("h2210","P vs X EYM2 st. 2",30,3.0,183.0,100,0.,5.);

  new TH2F("h2211","P vs X EPH cut st. 2",30,3.0,183.0,100,0.,0.05);
  new TH2F("h2212","P vs X EAL cut st. 2",30,3.0,183.0,100,0.,0.2);
  //new TH2F("h2213","P vs X EXM cut st. 2",30,3.0,183.0,100,0.,1.5);
  new TH2F("h2213","P vs X EXM cut st. 2",30,3.0,183.0,100,0.,11.);
  new TH2F("h2214","P vs X EYM cut st. 2",30,3.0,183.0,100,0.,10.);

  new TH2F("h2215","P2 vs X EPH cut st. 2",30,3.0,183.0,100,0.,0.05);
  new TH2F("h2216","P2 vs X EAL cut st. 2",30,3.0,183.0,100,0.,0.2);
  //new TH2F("h2217","P2 vs X EXM cut st. 2",30,3.0,183.0,100,0.,1.5);
  new TH2F("h2217","P2 vs X EXM cut st. 2",30,3.0,183.0,100,0.,11.);
  new TH2F("h2218","P2 vs X EYM cut st. 2",30,3.0,183.0,100,0.,10.);

  // 1
  new TH2F("h2102","P2 vs X EPH2 st. 2",30,3.0,183.0,100,0.,0.0006);
  new TH2F("h2103","P2 vs X EAL2 st. 2",30,3.0,183.0,100,0.,0.005);
  //new TH2F("h2104","P2 vs X EXM2 st. 2",30,3.0,183.0,100,0.,1.5);
  new TH2F("h2104","P2 vs X EXM2 st. 2",30,3.0,183.0,100,0.,7.);
  new TH2F("h2105","P2 vs X EYM2 st. 2",30,3.0,183.0,100,0.,7.);

  new TH2F("h2107","P vs X EPH2 st. 2",30,3.0,183.0,100,0.,0.0006);
  new TH2F("h2108","P vs X EAL2 st. 2",30,3.0,183.0,100,0.,0.005);
  //new TH2F("h2109","P vs X EXM2 st. 2",30,3.0,183.0,100,0.,1.5);
  new TH2F("h2109","P vs X EXM2 st. 2",30,3.0,183.0,100,0.,7.);
  new TH2F("h2110","P vs X EYM2 st. 2",30,3.0,183.0,100,0.,7.);

  new TH2F("h2111","P vs X EPH cut st. 2",30,3.0,183.0,100,0.,0.1);
  new TH2F("h2112","P vs X EAL cut st. 2",30,3.0,183.0,100,0.,0.2);
  //new TH2F("h2113","P vs X EXM cut st. 2",30,3.0,183.0,100,0.,1.5);
  new TH2F("h2113","P vs X EXM cut st. 2",30,3.0,183.0,100,0.,11.);
  new TH2F("h2114","P vs X EYM cut st. 2",30,3.0,183.0,100,0.,11.);

  new TH2F("h2115","P2 vs X EPH cut st. 2",30,3.0,183.0,100,0.,0.1);
  new TH2F("h2116","P2 vs X EAL cut st. 2",30,3.0,183.0,100,0.,0.2);
  //new TH2F("h2117","P2 vs X EXM cut st. 2",30,3.0,183.0,100,0.,1.5);
  new TH2F("h2117","P2 vs X EXM cut st. 2",30,3.0,183.0,100,0.,11.);
  new TH2F("h2118","P2 vs X EYM cut st. 2",30,3.0,183.0,100,0.,11.);

  // 2,3,4,5
  new TH1F("h2701","P2 fit 2",30,3.0,183.0);
  new TH2F("h2702","P2 vs X EPH2 st. 1 fit 2",30,3.0,183.0,100,0.,0.0006);
  new TH2F("h2703","P2 vs X EAL2 st. 1 fit 2",30,3.0,183.0,100,0.,0.005);
  // new TH2F("h2704","P2 vs X EXM2 st. 1 fit 2",30,3.0,183.0,100,0.,1.5);
  new TH2F("h2704","P2 vs X EXM2 st. 1 fit 2",30,3.0,183.0,100,0.,2.);
  new TH2F("h2705","P2 vs X EYM2 st. 1 fit 2",30,3.0,183.0,100,0.,3.);

  new TH2F("h2707","P vs X EPH2 st. 1 fit 2",30,3.0,183.0,100,0.,0.0006);
  new TH2F("h2708","P vs X EAL2 st. 1 fit 2",30,3.0,183.0,100,0.,0.005);
  //new TH2F("h2709","P vs X EXM2 st. 1 fit 2",30,3.0,183.0,100,0.,1.5);
  new TH2F("h2709","P vs X EXM2 st. 1 fit 2",30,3.0,183.0,100,0.,2.);
  new TH2F("h2710","P vs X EYM2 st. 1 fit 2",30,3.0,183.0,100,0.,3.);

  new TH2F("h2711","P vs X EPH cut st. 1 fit 2",30,3.0,183.0,100,0.,0.07);
  new TH2F("h2712","P vs X EAL cut st. 1 fit 2",30,3.0,183.0,100,0.,0.2);
  //new TH2F("h2713","P vs X EXM cut st. 1 fit 2",30,3.0,183.0,100,0.,1.5);
  new TH2F("h2713","P vs X EXM cut st. 1 fit 2",30,3.0,183.0,100,0.,6.);
  new TH2F("h2714","P vs X EYM cut st. 1 fit 2",30,3.0,183.0,100,0.,7.);

  new TH2F("h2715","P2 vs X EPH cut st. 1 fit 2",30,3.0,183.0,100,0.,0.07);
  new TH2F("h2716","P2 vs X EAL cut st. 1 fit 2",30,3.0,183.0,100,0.,0.2);
  //new TH2F("h2717","P2 vs X EXM cut st. 1 fit 2",30,3.0,183.0,100,0.,1.5);
  new TH2F("h2717","P2 vs X EXM cut st. 1 fit 2",30,3.0,183.0,100,0.,6.);
  new TH2F("h2718","P2 vs X EYM cut st. 1 fit 2",30,3.0,183.0,100,0.,7.);

  // 1,3,4,5
  new TH1F("h2801","P2 fit 1",30,3.0,183.0);
  new TH2F("h2802","P2 vs X EPH2 st. 2 fit 1",30,3.0,183.0,100,0.,0.0006);
  new TH2F("h2803","P2 vs X EAL2 st. 2 fit 1",30,3.0,183.0,100,0.,0.005);
  //new TH2F("h2804","P2 vs X EXM2 st. 2 fit 1",30,3.0,183.0,100,0.,1.5);
  new TH2F("h2804","P2 vs X EXM2 st. 2 fit 1",30,3.0,183.0,100,0.,2.);
  new TH2F("h2805","P2 vs X EYM2 st. 2 fit 1",30,3.0,183.0,100,0.,3.);

  new TH2F("h2807","P vs X EPH2 st. 2 fit 1",30,3.0,183.0,100,0.,0.0006);
  new TH2F("h2808","P vs X EAL2 st. 2 fit 1",30,3.0,183.0,100,0.,0.005);
  //new TH2F("h2809","P vs X EXM2 st. 2 fit 1",30,3.0,183.0,100,0.,1.5);
  new TH2F("h2809","P vs X EXM2 st. 2 fit 1",30,3.0,183.0,100,0.,2.);
  new TH2F("h2810","P vs X EYM2 st. 2 fit 1",30,3.0,183.0,100,0.,3.);

  new TH2F("h2811","P vs X EPH cut st. 2 fit 1",30,3.0,183.0,100,0.,0.05);
  new TH2F("h2812","P vs X EAL cut st. 2 fit 1",30,3.0,183.0,100,0.,0.2);
  //new TH2F("h2813","P vs X EXM cut st. 2 fit 1",30,3.0,183.0,100,0.,1.5);
  new TH2F("h2813","P vs X EXM cut st. 2 fit 1",30,3.0,183.0,100,0.,5.);
  new TH2F("h2814","P vs X EYM cut st. 2 fit 1",30,3.0,183.0,100,0.,7.);

  new TH2F("h2815","P2 vs X EPH cut st. 2 fit 1",30,3.0,183.0,100,0.,0.05);
  new TH2F("h2816","P2 vs X EAL cut st. 2 fit 1",30,3.0,183.0,100,0.,0.2);
  //new TH2F("h2817","P2 vs X EXM cut st. 2 fit 1",30,3.0,183.0,100,0.,1.5);
  new TH2F("h2817","P2 vs X EXM cut st. 2 fit 1",30,3.0,183.0,100,0.,5.);
  new TH2F("h2818","P2 vs X EYM cut st. 2 fit 1",30,3.0,183.0,100,0.,7.);


  new TH2F("h1111","dx vs x station 1",30,-250.,250.,30,-0.5,0.5);
  new TH2F("h1112","dx vs x station 2",30,-250.,250.,30,-0.5,0.5);
  new TH2F("h1113","dx vs x station 3",30,-250.,250.,30,-0.5,0.5);
  new TH2F("h1114","dx vs x station 4",30,-250.,250.,30,-0.5,0.5);
  new TH2F("h1115","dx vs x station 5",30,-250.,250.,30,-0.5,0.5);
  new TH2F("h1121","dy vs y station 1",30,-250.,250.,30,-0.04,0.04);
  new TH2F("h1122","dy vs y station 2",30,-250.,250.,30,-0.04,0.04);
  new TH2F("h1123","dy vs y station 3",30,-250.,250.,30,-0.04,0.04);
  new TH2F("h1124","dy vs y station 4",30,-250.,250.,30,-0.04,0.04);
  new TH2F("h1125","dy vs y station 5",30,-250.,250.,30,-0.04,0.04);

  // fin de test

  new TH1F("h500","Acceptance en H st. 4",500,0.,500.);
  new TH1F("h600","Acceptance en H st. 5",500,0.,500.);
  new TH1F("h700","X vertex track found",200,-10.,10.);
  new TH1F("h701","Y vertex track found",200,-10.,10.);
  new TH1F("h800","Rap. muon gen.",100,0.,5.);
  new TH1F("h801","Rap. muon gen. recons.",100,0.,5.);
  new TH1F("h802","Rap. muon gen. ghost ",100,0.,5.);
  new TH1F("h900","Pt muon gen.",100,0.,20.);
  new TH1F("h901","Pt muon gen. recons.",100,0.,20.);
  new TH1F("h902","Pt muon gen. ghost",100,0.,20.);
  new TH1F("h910","phi muon gen.",100,-10.,10.);
  new TH1F("h911","phi muon gen. recons.",100,-10.,10.);
  new TH1F("h912","phi muon gen. ghost",100,-10.,10.);
  new TH2F("h1001","Y VS X hit st. 1",300,-300.,300.,300,-300.,300.);
  new TH2F("h1002","Y VS X hit st. 2",300,-300.,300.,300,-300.,300.);
  new TH2F("h1003","Y VS X hit st. 3",300,-300.,300.,300,-300.,300.);
  new TH2F("h1004","Y VS X hit st. 4",300,-300.,300.,300,-300.,300.);
  new TH2F("h1005","Y VS X hit st. 5",300,-300.,300.,300,-300.,300.);
  //  Histos variance dans 4      
  new TH2F("h11","VAR X st. 4",30,3.0,183.0,100,0.,2.);
  new TH2F("h12","VAR Y st. 4",30,3.0,183.0,100,0.,600.);
  new TH2F("h13","VAR PHI st. 4",30,3.0,183.0,100,0.,0.0001);
  new TH2F("h14","VAR ALM st. 4",30,3.0,183.0,100,0.,0.05);
  new TH1F("h15","P",30,3.0,183.0);
  new TH1F("h411","VAR X st. 4",100,-1.42,1.42);
  new TH1F("h412","VAR Y st. 4",100,-25.,25.);
  new TH1F("h413","VAR PHI st. 4",100,-0.01,0.01);
  new TH1F("h414","VAR ALM st. 4",100,-0.23,0.23);
  // histo2
  new TH2F("h211","histo2-VAR X st. 4",30,3.0,183.0,100,0.,2.);
  new TH2F("h212","histo2-VAR Y st. 4",30,3.0,183.0,100,0.,600.);
  new TH1F("h213","histo2-VAR X st. 4",100,-1.42,1.42);
  new TH1F("h214","histo2-VAR Y st. 4",100,-25.,25.);
  new TH1F("h215","histo2-P",30,3.0,183.0);

  //  Histos variance dans 2      
  new TH2F("h21","VAR X st. 2",30,3.0,183.0,100,0.,3.);
  new TH2F("h22","VAR Y st. 2",30,3.0,183.0,100,0.,7.);
  new TH2F("h23","VAR PHI st. 2",30,3.0,183.0,100,0.,0.006);
  new TH2F("h24","VAR ALM st. 2",30,3.0,183.0,100,0.,0.005);
  new TH1F("h25","P",30,3.0,183.0);
  new TH1F("h421","VAR X st. 2",100,-1.72,1.72);
  new TH1F("h422","VAR Y st. 2",100,-2.7,2.7);
  new TH1F("h423","VAR PHI st. 2",100,-0.08,0.08);
  new TH1F("h424","VAR ALM st. 2",100,-0.072,0.072);
  // histo2
  new TH2F("h221","histo2-VAR X st. 2",30,3.0,183.0,100,0.,3.);
  new TH2F("h222","histo2-VAR Y st. 2",30,3.0,183.0,100,0.,7.);
  new TH1F("h223","histo2-VAR X st. 2",100,-1.72,1.72);
  new TH1F("h224","histo2-VAR Y st. 2",100,-2.7,2.7);
  new TH1F("h225","histo2-P",30,3.0,183.0);

  //  Histos variance dans 1      
  new TH2F("h31","VAR X st. 1",30,3.0,183.0,100,0.,2.);
  new TH2F("h32","VAR Y st. 1",30,3.0,183.0,100,0.,0.5);
  new TH2F("h33","VAR PHI st. 1",30,3.0,183.0,100,0.,0.006);
  new TH2F("h34","VAR ALM st. 1",30,3.0,183.0,100,0.,0.005);
  new TH1F("h35","P",30,3.0,183.0);
  new TH1F("h431","VAR X st. 1",100,-1.42,1.42);
  new TH1F("h432","VAR Y st. 1",100,-0.72,0.72);
  new TH1F("h433","VAR PHI st. 1",100,-0.08,0.08);
  new TH1F("h434","VAR ALM st. 1",100,-0.072,0.072);
  //  Histos variance dans 1      
  new TH2F("h41","VAR X st. 1 fit 5,4,3,2,V",30,3.0,183.0,100,0.,4.);
  new TH2F("h42","VAR Y st. 1 fit 5,4,3,2,V",30,3.0,183.0,100,0.,20.);
  new TH2F("h43","VAR PHI st. 1 fit 5,4,3,2,V",30,3.0,183.0,100,0.,0.005);
  new TH2F("h44","VAR ALM st. 1 fit 5,4,3,2,V",30,3.0,183.0,100,0.,0.005);
  new TH1F("h45","P",30,3.0,183.0);
  new TH1F("h441","VAR X st. 1 fit 5,4,3,2,V",100,-2.,2.);
  new TH1F("h442","VAR Y st. 1 fit 5,4,3,2,V",100,-4.5,4.5);
  new TH1F("h443","VAR PHI st. 1 fit 5,4,3,2,V",100,-0.072,0.072);
  new TH1F("h444","VAR ALM st. 1 fit 5,4,3,2,V",100,-0.072,0.072);
  // histo2
  new TH2F("h241","histo2-VAR X st. 1 fit 5,4,3,2,V",30,3.0,183.0,100,0.,4.);
  new TH2F("h242","histo2-VAR Y st. 1 fit 5,4,3,2,V",30,3.0,183.0,100,0.,20.);
  new TH1F("h243","histo2-VAR X st. 1 fit 5,4,3,2,V",100,-2.,2.);
  new TH1F("h244","histo2-VAR Y st. 1 fit 5,4,3,2,V",100,-4.5,4.5);
  new TH1F("h245","histo2-P",30,3.0,183.0);

  //  Histos variance dans 2      
  new TH2F("h51","VAR X st. 2 fit 5,4,3,1,V",30,3.0,183.0,100,0.,0.5);
  new TH2F("h52","VAR Y st. 2 fit 5,4,3,1,V",30,3.0,183.0,100,0.,2.);
  new TH2F("h53","VAR PHI st. 2 fit 5,4,3,1,V",30,3.0,183.0,100,0.,0.005);
  new TH2F("h54","VAR ALM st. 2 fit 5,4,3,1,V",30,3.0,183.0,100,0.,0.01);
  new TH1F("h55","P",30,3.0,183.0);
  new TH1F("h451","VAR X st. 2 fit 5,4,3,1,V",100,-0.72,0.72);
  new TH1F("h452","VAR Y st. 2 fit 5,4,3,1,V",100,-1.42,1.42);
  new TH1F("h453","VAR PHI st. 2 fit 5,4,3,1,V",100,-0.072,0.072);
  new TH1F("h454","VAR ALM st. 2 fit 5,4,3,1,V",100,-0.1,0.1);
  new TH1F("h999","PTOT",30,3.0,183.0);
  // histo2
  new TH2F("h251","histo2-VAR X st. 2 fit 5,4,3,1,V",30,3.0,183.0,100,0.,0.5);
  new TH2F("h252","histo2-VAR Y st. 2 fit 5,4,3,1,V",30,3.0,183.0,100,0.,2.);
  new TH1F("h253","histo2-VAR X st. 2 fit 5,4,3,1,V",100,-0.72,0.72);
  new TH1F("h254","histo2-VAR Y st. 2 fit 5,4,3,1,V",100,-1.42,1.42);
  new TH1F("h255","histo2-P",30,3.0,183.0);
  //  Histos variance dans 3      
  new TH2F("h61","VAR X st. 3 fit 4,5,V",30,3.0,183.0,100,0.,5.);
  new TH2F("h62","VAR Y st. 3 fit 4,5,V",30,3.0,183.0,100,0.,2.);
  new TH2F("h63","VAR PHI st. 3 fit 4,5,V",30,3.0,183.0,100,0.,0.0006);
  new TH2F("h64","VAR ALM st. 3 fit 4,5,V",30,3.0,183.0,100,0.,0.0006);
  new TH1F("h65","P",30,3.0,183.0);
  new TH1F("h461","VAR X st. 3 fit 4,5,V",100,-2.25,2.25);
  new TH1F("h462","VAR Y st. 3 fit 4,5,V",100,-1.42,1.42);
  new TH1F("h463","VAR PHI st. 3 fit 4,5,V",100,-0.024,0.024);
  new TH1F("h464","VAR ALM st. 3 fit 4,5,V",100,-0.024,0.024);
  // histo2
  new TH2F("h261","histo2-VAR X st. 3 fit 4,5,V",30,3.0,183.0,100,0.,5.);
  new TH2F("h262","histo2-VAR Y st. 3 fit 4,5,V",30,3.0,183.0,100,0.,2.);
  new TH1F("h263","histo2-VAR X st. 3 fit 4,5,V",100,-2.25,2.25);
  new TH1F("h264","histo2-VAR Y st. 3 fit 4,5,V",100,-1.42,1.42);
  new TH1F("h265","Phisto2-",30,3.0,183.0);
  // Histos dx,dy distribution between chambers inside stations
  new TH1F("h71","DX in st. ID-70",100,-5.,5.);
  new TH1F("h81","DY in st. ID-80",100,-5.,5.);
  new TH1F("h72","DX in st. ID-70",100,-5.,5.);
  new TH1F("h82","DY in st. ID-80",100,-5.,5.);
  new TH1F("h73","DX in st. ID-70",100,-5.,5.);
  new TH1F("h83","DY in st. ID-80",100,-5.,5.);
  new TH1F("h74","DX in st. ID-70",100,-5.,5.);
  new TH1F("h84","DY in st. ID-80",100,-5.,5.);
  new TH1F("h75","DX in st. ID-70",100,-5.,5.);
  new TH1F("h85","DY in st. ID-80",100,-5.,5.);
}

//_____________________________________________________________________________
void chfnt(Int_t &ievr, Int_t &ntrackr, Int_t *istatr, Int_t *isignr, Float_t *pxr, Float_t *pyr, Float_t *pzr, Float_t *zvr, Float_t *chi2r,  Float_t *pxv, Float_t *pyv, Float_t *pzv)
{
  //
  // fill the ntuple 
    NtupleSt.ievr = ievr;
    NtupleSt.ntrackr = ntrackr;
    for (Int_t i=0; i<500; i++) {
	NtupleSt.istatr[i] = istatr[i];
	NtupleSt.isignr[i] = isignr[i]; 
	NtupleSt.pxr[i]    = pxr[i]; 
	NtupleSt.pyr[i]    = pyr[i];
	NtupleSt.pzr[i]    = pzr[i];
	NtupleSt.zvr[i]    = zvr[i];
	NtupleSt.chi2r[i]  = chi2r[i];
	NtupleSt.pxv[i]    = pxv[i]; 
	NtupleSt.pyv[i]    = pyv[i];
	NtupleSt.pzv[i]    = pzv[i];
    }
    gAliNtupleGlobal->Fill();   
}

//________________
void hist_closed()
{
  //
  // write histos and ntuple to "reconst.root" file
  gAliFileGlobal->Write();
}

//________________________________________________________________________
void trackf_read_fit(Int_t &ievr, Int_t &nev, Int_t &nhittot1, Int_t *izch, Double_t *xgeant, Double_t *ygeant) 
{

  // introduce aliroot variables in fortran common 
  // tracking study from geant hits 
  //

  AliMUON *pMUON  = (AliMUON*) gAlice->GetModule("MUON");
  TTree *treeH = gAlice->TreeH();
  Int_t ntracks = (Int_t)treeH->GetEntries();
  cout<<"ntrack="<<ntracks<<endl;

  nhittot1 = 0;

//  Loop over tracks
  for (Int_t track=0; track<ntracks;track++) {
      gAlice->ResetHits();
      treeH->GetEvent(track);
      
      if (pMUON)  {

//  Loop over hits
	  for(AliMUONHit* mHit=(AliMUONHit*)pMUON->FirstHit(-1); 
	      mHit;
	      mHit=(AliMUONHit*)pMUON->NextHit()) 
	  {
	      if (mHit->Chamber() > 10) continue;
	      Int_t ftrack = mHit->Track();
	      Int_t id = gAlice->Particle(ftrack)->GetPdgCode();
	      
	      if (id==kMuonPlus||id==kMuonMinus) {
		  xgeant[nhittot1]   = mHit->Y();
		  ygeant[nhittot1]   = mHit->X();
		  izch[nhittot1]     = mHit->Chamber();
//		  printf("id %d ch %d x %f y %f\n",id,izch[nhittot1],xgeant[nhittot1],ygeant[nhittot1]);  
		  nhittot1++;
	      }
	  } // hit loop
      } // if pMUON
  } // track loop

  ievr=nev;
  gAliFileGlobal->cd();    
}

//______________________________________________________________________________
void trackf_read_geant(Int_t *itypg, Double_t *xtrg, Double_t *ytrg, Double_t *ptotg, Int_t *idg, Int_t *izch, Double_t *pvert1g, Double_t *pvert2g, Double_t *pvert3g, Double_t *zvertg, Int_t &nhittot1, Double_t *cx, Double_t *cy, Double_t *cz, Int_t &ievr,Int_t &nev,Double_t *xgeant, Double_t *ygeant, Double_t *clsize1, Double_t *clsize2) 
{
  //
  // introduce aliroot variables in fortran common 
  // tracking study from geant hits 
  //

  AliMUON *pMUON  = (AliMUON*) gAlice->GetModule("MUON");
  
  //  TTree *treeK = gAlice->TreeK();
  TTree *treeH = gAlice->TreeH();
  Int_t ntracks = (Int_t)treeH->GetEntries();
  cout<<"ntrack="<<ntracks<<endl;

  Int_t maxidg = 0;
  Int_t nres=0;
  
//
//  Loop over tracks
//

  for (Int_t track=0; track<ntracks;track++) {
      gAlice->ResetHits();
      treeH->GetEvent(track);
      
      if (pMUON)  {
//
//  Loop over hits
//
	  for(AliMUONHit* mHit=(AliMUONHit*)pMUON->FirstHit(-1); 
	      mHit;
	      mHit=(AliMUONHit*)pMUON->NextHit()) 
	  {
	      if (maxidg<=20000) {
		
		if (mHit->Chamber() > 10) continue;
		TParticle *particle;
		Int_t ftrack = mHit->Track();
		Int_t id = gAlice->Particle(ftrack)->GetPdgCode();

//		if (id==kMuonPlus||id==kMuonMinus) {
		    
		    // inversion de x et y car le champ est inverse dans le programme tracking
		    xtrg[maxidg]   = 0;       
		    ytrg[maxidg]   = 0;       
		    xgeant[maxidg]   = mHit->Y();             // x-pos of hit
		    ygeant[maxidg]   = mHit->X();             // y-pos of hit
		    clsize1[maxidg]   = 0;     // cluster size on 1-st cathode
		    clsize2[maxidg]   = 0;     // cluster size on 2-nd cathode
		    cx[maxidg]     = mHit->Cy();              // Px/P of hit
		    cy[maxidg]     = mHit->Cx();              // Py/P of hit
		    cz[maxidg]     = mHit->Cz();              // Pz/P of hit
		    izch[maxidg]   = mHit->Chamber();         
		    /*      
		    Int_t pdgtype  = Int_t(mHit->fParticle); // particle number
		    itypg[maxidg]  = gMC->IdFromPDG(pdgtype);
		    */
		    itypg[maxidg] = 0;
                    if (id==kMuonPlus) itypg[maxidg]  = 5;
		    if (id==kMuonMinus) itypg[maxidg]  = 6;

                    //printf("ich, itypg[maxidg] %d %d\n",izch[maxidg],itypg[maxidg]);

		    ptotg[maxidg]  = mHit->Momentum();          // P of hit 
		    
		    particle = gAlice->Particle(ftrack);
		    Float_t thet = particle->Theta();
		    thet = thet*180./3.1416;
		    
		    //cout<<"chambre "<<izch[maxidg]<<"  ptot="<<ptotg[maxidg]<<"   theta="<<thet<<"   phi="<<mHit->fPhi<<" z="<<zz<<endl;	    
		    
		    Int_t iparent = particle->GetFirstMother();
		    if (iparent >= 0) {
			Int_t ip;
			while(1) {
			    ip=gAlice->Particle(iparent)->GetFirstMother();
			    if (ip < 0) {
				break;
			    } else {
				iparent = ip;
			    }
			}
		    }
		    //printf("iparent - %d\n",iparent);
		    Int_t id1  = ftrack; // numero de la particule generee au vertex
		    Int_t idum = track+1;
		    Int_t id2 = gAlice->Particle(iparent)->GetPdgCode();

		    if (id2==443) id2=114;
		    else id2=116;
                   
                    if (id2==116) {
		      nres++;
		    }
		    //printf("id2 %d\n",id2);
		    idg[maxidg] = 30000*id1+10000*idum+id2;
		    
		    pvert1g[maxidg] = particle->Py();      // Px vertex
		    pvert2g[maxidg] = particle->Px();      // Py vertex  
		    pvert3g[maxidg] = particle->Pz();      // Pz vertex
		    zvertg[maxidg]  = particle->Vz();      // z vertex 
	    
		    //	    cout<<"x="<<xgeant[maxidg]<<endl;
		    //cout<<"y="<<ygeant[maxidg]<<endl;
		    //cout<<"typ="<<itypg[maxidg]<<endl;

		    maxidg ++;

		}
	      }
	  } // hit loop
//      } // if pMUON
  } // track loop first file

  if (gAliTrH1 && gAliHits2 ) { // if background file
    ntracks =(Int_t)gAliTrH1->GetEntries();
    printf("Trackf_read - 2-nd file - ntracks %d\n",ntracks);

    //  Loop over tracks
    for (Int_t track=0; track<ntracks; track++) {
      
      if (gAliHits2) gAliHits2->Clear();
      gAliTrH1->GetEvent(track);

      //  Loop over hits
      AliMUONHit *mHit;
      for (int i=0;i<gAliHits2->GetEntriesFast();i++) 
	{
	  mHit=(AliMUONHit*) (*gAliHits2)[i];
	  if (mHit->Chamber() > 10) continue;
	  if (maxidg<=20000) {
	    
	    // inversion de x et y car le champ est inverse dans le programme tracking !!!!
	    xtrg[maxidg]   = 0;                    // only for reconstructed point
	    ytrg[maxidg]   = 0;                    // only for reconstructed point
	    xgeant[maxidg]   = mHit->Y();           // x-pos of hit
	    ygeant[maxidg]   = mHit->X();           // y-pos of hit
	    clsize1[maxidg]   = 0;           // cluster size on 1-st cathode
	    clsize2[maxidg]   = 0;           // cluster size on 2-nd cathode
	    cx[maxidg]     = mHit->Cy();            // Px/P of hit
	    cy[maxidg]     = mHit->Cx();            // Py/P of hit
	    cz[maxidg]     = mHit->Cz();            // Pz/P of hit
	    izch[maxidg]   = mHit->Chamber();       // chamber number
	    ptotg[maxidg]  = mHit->Momentum();      // P of hit 
	    
	    Int_t ftrack = mHit->Track();
	    Int_t id1  = ftrack;                   // track number 
	    Int_t idum = track+1;
	    
	    TClonesArray *fPartArray = gAliParticles2;
//	    TParticle *Part=NULL;
	    Int_t id = ((TParticle*) fPartArray->UncheckedAt(ftrack))->GetPdgCode();
	    if (id==kMuonPlus||id==kMuonMinus) {
	        if (id==kMuonPlus) itypg[maxidg]  = 5;
	        else  itypg[maxidg]  = 6;
	    } else itypg[maxidg]=0;
	    
	    Int_t id2=0; // set parent to 0 for background !!
	    idg[maxidg] = 30000*id1+10000*idum+id2;
	    

	    pvert1g[maxidg] = 0;      // Px vertex
	    pvert2g[maxidg] = 0;      // Py vertex  
	    pvert3g[maxidg] = 0;      // Pz vertex
	    zvertg[maxidg]  = 0;      // z vertex 
	    maxidg ++;

	  } // check limits (maxidg)
	} // hit loop 
    } // track loop
  } // if gAliTrH1

  ievr = nev;
  nhittot1 = maxidg ;
  cout<<"nhittot1="<<nhittot1<<endl;

  static Int_t nbres=0;
  if (nres>=19) nbres++;
  printf("nres, nbres %d %d \n",nres,nbres);

  gAliFileGlobal->cd();      

}

//________________________________________________________________________
void trackf_read_spoint(Int_t *itypg, Double_t *xtrg, Double_t *ytrg, Double_t *ptotg, Int_t *idg, Int_t *izch, Double_t *pvert1g, Double_t *pvert2g, Double_t *pvert3g, Double_t *zvertg, Int_t &nhittot1, Double_t *cx, Double_t *cy, Double_t *cz, Int_t &ievr,Int_t &nev,Double_t *xgeant, Double_t *ygeant,Double_t *clsize1, Double_t *clsize2) 

{
    //
    // introduce aliroot variables in fortran common 
    // tracking study from reconstructed points 
    //
    AliMUON *pMUON  = (AliMUON*) gAlice->GetModule("MUON");
    
    cout<<"numero de l'evenement "<<nev<<endl;
    
    TTree *treeR = gAlice->TreeR();
    Int_t nent=(Int_t)treeR->GetEntries();
    if (nev < 10) 
	printf("Found %d entries in the tree (must be one per cathode per event! + 1empty)\n",
	       nent);
//
    
    Int_t mult1, mult2;
    
    if (pMUON) {
	Int_t mpoi=0;
	for (Int_t ich=0;ich<10;ich++) {
	    printf("chambre %d\n",ich+1);
	    TClonesArray *reconstPoints  = pMUON->RawClustAddress(ich);

	    pMUON->ResetRawClusters();
	    treeR->GetEvent(nent-1);
	    Int_t npoints = (Int_t) reconstPoints->GetEntries();
	    if (!npoints) continue;
	    printf("\n ch %d npoints = %d\n",ich+1,npoints);
	    //  Loop over reconstruted points
	    for (Int_t ipoi=0; ipoi<npoints; ipoi++) {
		printf("  point %d\n",ipoi);
		AliMUONRawCluster* point = 
		    (AliMUONRawCluster*) reconstPoints->UncheckedAt(ipoi);
		
		mult1=point->fMultiplicity[0];
		mult2=point->fMultiplicity[1];	      
		xtrg[mpoi]=(Double_t) point->fY[0];
		ytrg[mpoi]=(Double_t) point->fX[0];
		izch[mpoi]=ich+1;
		Int_t itrack  = point->fTracks[1];
		Int_t ihit    = point->fTracks[0];
		xgeant[mpoi] = 0;
		ygeant[mpoi] = 0;
		clsize1[mpoi] = mult1;
		clsize2[mpoi] = mult2;
		Int_t id1, id2, idum;
		id1=id2=idum=-1;
		itypg[mpoi]=0;
		ihit = ihit-1;
		if (ihit >=0 && itrack >=0) {
		    gAlice->ResetHits();
		    gAlice->TreeH()->GetEvent(itrack);
		    TClonesArray *pMUONhits  = pMUON->Hits();
		    AliMUONHit* mHit;
		    mHit=(AliMUONHit*) (pMUONhits->UncheckedAt(ihit));
		    Int_t id = (Int_t) mHit->Particle();
		    xgeant[mpoi] = mHit->Y();          
		    ygeant[mpoi] = mHit->X(); 
		    if (id == kMuonPlus)  itypg[mpoi]=5;
		    if (id == kMuonMinus) itypg[mpoi]=6;
		    TParticle *particle;
		    particle = gAlice->Particle(mHit->Track());
		    TParticle* particleM=gAlice->Particle(particle->GetFirstMother());
		    Int_t iparent=particleM->GetPdgCode();
		    printf("\n Particle Id:%d %d \n", id, iparent);
		    if (iparent == 443) id2=114;
		    if (iparent == 553) id2=116;
		}
		id1=itrack;
		idum=itrack+1;
		idg[mpoi] = 30000*id1+10000*idum+id2;
		mpoi++;
	    } // loop over points
	} // loop over chamber
	ievr = nev;
	cout<<"evenement "<<ievr<<endl;
	nhittot1 = mpoi;
	cout<<"nhittot1="<<nhittot1<<endl;
	
	treeR->Reset();

	gAliFileGlobal->cd();

    } // if pMUON
}

//____________________________________________________________________________
void trackf_fit(Int_t &ivertex, Double_t *pest, Double_t *pstep, Double_t &pxzinv, Double_t &tphi, Double_t &talam, Double_t &xvert, Double_t &yvert)
{
  //
  //  Fit a track candidate with the following input parameters: 
  //  INPUT :  IVERTEX  : vertex flag, if IVERTEX=1 (XVERT,YVERT) are free paramaters
  //                                   if IVERTEX=1 (XVERT,YVERT)=(0.,0.) 
  //           PEST(5)  : starting value of parameters (minuit)
  //           PSTEP(5) : step size for parameters (minuit)
  //  OUTPUT : PXZINV,TPHI,TALAM,XVERT,YVERT : fitted value of the parameters

  static Double_t arglist[10];
  static Double_t c[5] = {0.4, 0.45, 0.45, 90., 90.};
  static Double_t b1, b2, epxz, efi, exs, exvert, eyvert;
  TString chname;
  Int_t ierflg = 0;
  
  TMinuit *gMinuit = new TMinuit(5);
  gMinuit->mninit(5,10,7);
  gMinuit->SetFCN(fcnf);  // constant m.f.

  arglist[0] = -1;
  
  gMinuit->mnexcm("SET PRINT", arglist, 1, ierflg);
  //      gMinuit->mnseti('track fitting');
  
  gMinuit->mnparm(0, "invmom",  pest[0], pstep[0], -c[0], c[0], ierflg);
  gMinuit->mnparm(1, "azimuth", pest[1], pstep[1], -c[1], c[1], ierflg);
  gMinuit->mnparm(2, "deep",    pest[2], pstep[2], -c[2], c[2], ierflg);
  if (ivertex==1) {
    gMinuit->mnparm(3, "x ", pest[3], pstep[3], -c[3], c[3], ierflg);
    gMinuit->mnparm(4, "y ", pest[4], pstep[4], -c[4], c[4], ierflg);  
  }   
  
  gMinuit->mnexcm("SET NOGR", arglist, 0, ierflg);
  gMinuit->mnexcm("MINIMIZE", arglist, 0, ierflg);
  gMinuit->mnexcm("EXIT" , arglist, 0, ierflg);
  
  gMinuit->mnpout(0, chname, pxzinv, epxz, b1, b2, ierflg);
  gMinuit->mnpout(1, chname, tphi, efi, b1, b2, ierflg);
  gMinuit->mnpout(2, chname, talam, exs, b1, b2, ierflg);
  if (ivertex==1) {
    gMinuit->mnpout(3, chname, xvert, exvert, b1, b2, ierflg);
    gMinuit->mnpout(4, chname, yvert, eyvert, b1, b2, ierflg);
  }   
  
  delete gMinuit;
}
	   
//________________________________________________________________________________
void fcnf(Int_t &npar, Double_t *grad, Double_t &fval, Double_t *pest, Int_t iflag)
{
  //
  // function called by trackf_fit
      Int_t futil = 0;
      fcn(npar,grad,fval,pest,iflag,futil);
}

//____________________________________________________________________________
void prec_fit(Double_t &pxzinv, Double_t &fis, Double_t &alams, Double_t &xvert, Double_t &yvert, Double_t &pxzinvf, Double_t &fif, Double_t &alf, Double_t &xvertf, Double_t &yvertf, Double_t &epxzinv, Double_t &efi, Double_t &exs, Double_t &exvert, Double_t &eyvert)
{
  // 
  // minuit fits for tracking finding 
                                                                            
      static Double_t arglist[10];
      static Double_t c1[5] = {0.001, 0.001, 0.001, 1., 1.};
      static Double_t c2[5] = {0.5, 0.5, 0.5, 120., 120.};
      static Double_t emat[9];
      static Double_t b1, b2;
      Double_t fmin, fedm, errdef; 
      Int_t npari, nparx, istat;

      TString chname;
      Int_t ierflg = 0;
      
      TMinuit *gMinuit = new TMinuit(5);
      gMinuit->mninit(5,10,7);
      gMinuit->SetFCN(fcnfitf);

      arglist[0] = -1.;
      gMinuit->mnexcm("SET PRINT", arglist, 1, ierflg);
      
      //      gMinuit->mnseti('track fitting');

      gMinuit->mnparm(0,"invmom",   pxzinv, c1[0], -c2[0], c2[0], ierflg); // 0.003, 0.5
      gMinuit->mnparm(1,"azimuth ", fis,    c1[1], -c2[1], c2[1], ierflg);
      gMinuit->mnparm(2,"deep    ", alams,  c1[2], -c2[2], c2[2], ierflg);
      gMinuit->mnparm(3,"xvert",    xvert,  c1[3], -c2[3], c2[3], ierflg);
      gMinuit->mnparm(4,"yvert",    yvert,  c1[4], -c2[4], c2[4], ierflg);

      gMinuit->mnexcm("SET NOGR", arglist, 0, ierflg);
      arglist[0] = 2.;
      gMinuit->mnexcm("MINIMIZE", arglist, 0, ierflg);
      gMinuit->mnexcm("EXIT", arglist, 0, ierflg);
 
      gMinuit->mnpout(0, chname, pxzinvf, epxzinv, b1, b2, ierflg);
      gMinuit->mnpout(1, chname, fif, efi, b1, b2, ierflg);
      gMinuit->mnpout(2, chname, alf, exs, b1, b2, ierflg);
      gMinuit->mnpout(3, chname, xvertf, exvert, b1, b2, ierflg);
      gMinuit->mnpout(4, chname, yvertf, eyvert, b1, b2, ierflg);
  
      gMinuit->mnemat(emat, 3);
      gMinuit->mnstat(fmin, fedm, errdef, npari, nparx, istat);

     delete gMinuit;
}

//____________________________________________________________________
void fcnfitf(Int_t &npar, Double_t *grad, Double_t &fval, Double_t *xval, Int_t iflag)
{
  //
  // function called by prec_fit 
      Int_t futil = 0;
      fcnfit(npar,grad,fval,xval,iflag,futil);
}
