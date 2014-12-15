/**************************************************************************
 * Copyright(c) 1998-2010, ALICE Experiment at CERN, All rights reserved. *
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

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Wsu Cosmic Ray SetUp                                                     //
//  This class contains the description of the  Wsu Cosmic Ray SetUp         //
//  external volume                                                          //
//                                                                           //
//Begin_Html
/*
<img src="picts/AliEMCALWsuCosmicRaySetUpClass.gif">
</pre>
<br clear=left>
<font size=+2 color=red>
<p>The responsible person for this module is
<a href="mailto:pavlinov@physics.wayne.edu">Alexei Pavlino, WSU</a>.
</font>
<pre>
*/
//End_Html
//                                                                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <TVirtualMC.h>
#include <TROOT.h>
#include <TDatabasePDG.h>
#include <TBrowser.h>
#include <TLorentzVector.h>
#include <TParticle.h>
#include <TList.h>
#include <TH1.h>
#include <TH2.h>
#include <TAxis.h>

#include "AliEMCALWsuCosmicRaySetUp.h"
//#include "AliMagF.h"
#include "AliStack.h"
#include "AliRun.h"
#include "AliMC.h"

#include "AliEMCALHistoUtilities.h"

using namespace std;

typedef AliEMCALHistoUtilities hist;

TDatabasePDG *pdg = 0; 

Int_t gid=0;

ClassImp(AliEMCALWsuCosmicRaySetUp)

//  TList *fLHists=0, *ll=0; 

//_____________________________________________________________________________
  AliEMCALWsuCosmicRaySetUp::AliEMCALWsuCosmicRaySetUp(): AliModule()
						      //,fMasterVolume()
						    ,fLHists(0),fMasterVolume()
{
  //
  // Default constructor
  //
}
 
//_____________________________________________________________________________
AliEMCALWsuCosmicRaySetUp::AliEMCALWsuCosmicRaySetUp(const char *name, const char *title)
       : AliModule(name,title)
						    //,fMasterVolume()
						    ,fLHists(0),fMasterVolume()
{
  //
  // Standard constructor of the  Wsu Cosmic Ray SetUp external volume
  //
  //PH  SetMarkerColor(7);
  //PH  SetMarkerStyle(2);
  //PH  SetMarkerSize(0.4);
}
 
//_____________________________________________________________________________
void AliEMCALWsuCosmicRaySetUp::CreateGeometry()
{
  //
  // Create the geometry of the Alice external body
  //
  //Begin_Html
  /*
    <img src="picts/AliEMCALWsuCosmicRaySetUpTree.gif">
  */
  //End_Html
  pdg = TDatabasePDG::Instance();

  // Master Volume
  fMasterVolume[0] = fMasterVolume[1] = 35.0;
  fMasterVolume[2] = 450.;

  Int_t *idtmed = fIdtmed->GetArray()+1;
  int idAir = idtmed[0];
  TVirtualMC::GetMC()->Gsvolu(GetName(),"BOX",idAir, fMasterVolume,3); // Master volume
  //
  // Sc counters
  //
  Float_t sc[3]; // tube
  sc[0] = 0.0;
  sc[1] = 10.0;
  sc[2] = 1.0; // thicness of Sc is 2 cm
  Float_t zsc[3] = {10.,330.6, 810.1}; 
  int idSC = idtmed[1];
  TVirtualMC::GetMC()->Gsvolu("SCOU","TUBE",idSC, sc,3); // Master volume
  printf(" idtmed[0] %i idtmed[1] %i \n", idtmed[0] , idtmed[1]); 
  Int_t idRot=0; // no rotation
  for(Int_t i=0; i<3; i++) {
    Float_t zpos = zsc[i] - fMasterVolume[2];
    TVirtualMC::GetMC()->Gspos("SCOU", i+1, "WSUC", 0.0, 0.0, zpos, idRot, "ONLY"); 
  }
  //
  // Dead end : Dec 2,2010
  //
  Float_t zbox[3]={30., 30.0, 0.1};
  TVirtualMC::GetMC()->Gsvolu("SEND","BOX",idAir, zbox,3); // Master volume
  TVirtualMC::GetMC()->Gspos("SEND", 1, "WSUC", 0.0, 0.0, 448.0, idRot, "ONLY"); 
  // Hists
  fLHists = new TList;
  fLHists->SetName("hists");
  //
  //AliMC *ALIMC  = dynamic_cast<AliMC *>(TVirtualMC::GetMC());
  //AliGenBox* gB = dynamic_cast<AliGenBox *>(ALIMC->Generator());
  //Double_t p = gB->
  Double_t pmom=1.5; 
  fLHists->Add(BookKineHists(pmom,"primeKineHists"));
  fLHists->Add(BookKineHists(pmom,"endKineHists"));
  fLHists->Add(BookKineHists(pmom,"secondaryKineHists"));
  //ll = BookKineHists(1.,"kineHists");
  //gROOT->GetListOfBrowsables()->Add(ll);
}
 
//_____________________________________________________________________________
void AliEMCALWsuCosmicRaySetUp::CreateMaterials()
{
// Create materials and media
  Int_t   isxfld = 0;
  Float_t sxmgmx = 0., deemax = 0.1;  
  // AIR
  Float_t aAir[4]={12.0107,14.0067,15.9994,39.948};
  Float_t zAir[4]={6.,7.,8.,18.};
  Float_t wAir[4]={0.000124,0.755267,0.231781,0.012827};
  Float_t dAir = 1.20479E-3;
  AliMixture(1,"Air     $",aAir,zAir,dAir,4,wAir);

  // --- The polysterene scintillator (CH) ---
  Float_t aP[2] = {12.011, 1.00794} ;
  Float_t zP[2] = {6.0, 1.0} ;
  Float_t wP[2] = {1.0, 1.0} ;
  Float_t dP = 1.032 ;
  AliMixture(2, "Polystyrene$", aP, zP, dP, -2, wP) ;

  //
  AliMedium(1,"Air     $",1,0,isxfld,sxmgmx,10,-1,-0.1,0.1 ,-10);
  AliMedium(2, "Scintillator$", 2, 1,
            isxfld, sxmgmx, 10.0, 0.001, deemax, 0.001, 0.001, 0, 0) ;
  // cuts
  DefineCuts(1);
  DefineCuts(2);
}
 
void AliEMCALWsuCosmicRaySetUp::DefineCuts(const Int_t idtmed)
{
  // Dec 2,2010 : it works
  Float_t cutgam=10.e-5; // 100 kev;
  Float_t cutele=10.e-5; // 100 kev;
  TVirtualMC::GetMC()->Gstpar(idtmed,"CUTGAM", cutgam);
  TVirtualMC::GetMC()->Gstpar(idtmed,"CUTELE", cutele); // 1MEV -> 0.1MEV; 15-aug-05
  TVirtualMC::GetMC()->Gstpar(idtmed,"BCUTE",  cutgam);  // BCUTE and BCUTM start from GUTGUM
  TVirtualMC::GetMC()->Gstpar(idtmed,"BCUTM",  cutgam);  // BCUTE and BCUTM start from GUTGUM
  // --- Generate explicitly delta rays in Lead ---
  TVirtualMC::GetMC()->Gstpar(idtmed, "LOSS", 3) ;
  TVirtualMC::GetMC()->Gstpar(idtmed, "DRAY", 1) ;
  TVirtualMC::GetMC()->Gstpar(idtmed, "DCUTE", cutele) ;
  TVirtualMC::GetMC()->Gstpar(idtmed, "DCUTM", cutele) ;
}
 
void AliEMCALWsuCosmicRaySetUp::StepManager(void)
{
  // Dec 1,2010
  static Int_t pri=1;
  static TString curVolName="";
  static TLorentzVector pos;  // Lorentz vector of the track current position.
  static TLorentzVector mom;  // Lorentz vector of the track current momentum.

  if(pri>=2) printf("<I> AliEMCALWsuCosmicRaySetUp::StepManager %s \n", TVirtualMC::GetMC()->CurrentVolName());
  Int_t tracknumber =  gAlice->GetMCApp()->GetCurrentTrackNumber();
  //  Int_t parent=0;
  TParticle* part=0;
  curVolName = TVirtualMC::GetMC()->CurrentVolName();
  if(curVolName.Contains("SEND")) {
    TVirtualMC::GetMC()->TrackMomentum(mom);
    TVirtualMC::GetMC()->TrackPosition(pos);
    if(pri>=2) printf(" %s tracknumber %i p %f \n", curVolName.Data(), tracknumber, mom.P());
    if(pri>=2) printf(" x %f y %f z %f \n", pos[0], pos[1], pos[2]);
    if(TVirtualMC::GetMC()->IsTrackEntering()) { // primary only TList *l = GetLhists(1);
      TList *l = 0;
      if(tracknumber==0){
        l = GetLhists(1);
        part=gAlice->GetMCApp()->Particle(0);
        gid = pdg->ConvertPdgToGeant3(part->GetPdgCode()); 
        hist::FillH1(l, 1, double(gid));
      } else {
        l = GetLhists(2);
        part=gAlice->GetMCApp()->Particle(tracknumber);
        gid = pdg->ConvertPdgToGeant3(part->GetPdgCode()); 
        hist::FillH1(l, 1, double(gid));
      }
      Int_t ic = 2;
      hist::FillH1(l, ic++, mom.P());
      hist::FillH1(l, ic++, mom.Eta());
      hist::FillH1(l, ic++, TVector2::Phi_0_2pi(mom.Phi())*TMath::RadToDeg() );
      hist::FillH1(l, ic++, mom.Theta()*TMath::RadToDeg());
    }
  }
}

void AliEMCALWsuCosmicRaySetUp::FinishEvent()
{
  // Dec 2,2010
  int ic=0;

  TList *l = GetLhists(0);
  //  TList *l = ll;
  AliStack* st =  AliRunLoader::Instance()->Stack();
  TParticle *p = st->Particle(0);
  gid = pdg->ConvertPdgToGeant3(p->GetPdgCode());

  ic=1;
  hist::FillH1(l, ic++, double(gid));
  hist::FillH1(l, ic++, p->P());
  hist::FillH1(l, ic++, p->Eta());
  hist::FillH1(l, ic++, TVector2::Phi_0_2pi(p->Phi())*TMath::RadToDeg() );
  hist::FillH1(l, ic++, p->Theta()*TMath::RadToDeg());
}

TList* AliEMCALWsuCosmicRaySetUp::BookKineHists(const Double_t p , const Char_t *opt)
{
  // Dec 2,2010
  gROOT->cd();
  TH1::AddDirectory(1);

  TH1 * hgid=0;
  Int_t nphi=180, nmax=1100;
  Double_t phimin=0.0, phimax=360.;
  Double_t pmax=110.;
  if(p>0.1) pmax = 1.1*p;
  new TH1F("00_hNPrim"," number of primary particles ", 10, 0.5, 10.5);
  hgid = new TH1F("01_hGidprim","Geant Id of primary particles", 16, 0.5, 16.5);
  new TH1F("02_hPmomPrim","p of primary particles", nmax, 0.0, pmax);
  new TH1F("03_hEtaPrim","#eta primary particles", 80, 0.0, 8.0);
  new TH1F("04_hPhiPrim","#phi primary particles", nphi,phimin,phimax);
  new TH1F("05_hThetaPrim","#theta primary particles", 90, 0.0, 90.);
  //
  TAxis *xax=hgid->GetXaxis();
  xax->SetBinLabel(1,"#gamma");
  xax->SetBinLabel(2,"e^{+}");
  xax->SetBinLabel(3,"e^{-}");
  xax->SetBinLabel(4,"#nu");
  xax->SetBinLabel(5,"#mu^{+}");
  xax->SetBinLabel(6,"#mu^{-}");
  xax->SetBinLabel(7,"#pi^{0}");
  xax->SetBinLabel(8,"#pi^{+}");
  xax->SetBinLabel(9,"#pi^{-}");
  xax->SetBinLabel(10,"K_{L}");
  xax->SetBinLabel(11,"K^{+}");
  xax->SetBinLabel(12,"K^{-}");
  xax->SetBinLabel(13,"n");
  xax->SetBinLabel(14,"p");
  xax->SetBinLabel(15,"#bar{p}");
  xax->SetBinLabel(16,"K_{s}");
  //  hgid->SetBinLabel(,"");

  TList *l = 0;
  l = hist::MoveHistsToList(opt, 1);
  if(strlen(opt)) hist::AddToNameAndTitleToList(l, opt, opt);

  TH1::AddDirectory(0);

  return l;
}

void AliEMCALWsuCosmicRaySetUp::Browse(TBrowser* b)
{
  if(fLHists) b->Add(fLHists);
  TObject::Browse(b);
}
/*
*/
