/**************************************************************************
 * Copyright(c) 1998-2004, ALICE Experiment at CERN, All rights reserved. *
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

//_________________________________________________________________________
//*-- Implementation version v2 of EMCAL Manager class; SHASHLYK version
//*-- An object of this class does not produce digits
//*-- It is the one to use if you do want to produce outputs in TREEH 
//*--                  
//*-- Author : Aleksei Pavlinov (WSU)

// This Class not stores information on all particles prior to EMCAL entry - in order to facilitate analysis.
// This is done by setting fIShunt =2, and flagging all parents of particles entering the EMCAL.

#include <cassert>
// --- ROOT system ---
#include <TBrowser.h>
#include <TClonesArray.h>
#include <TH1.h>
#include <TH2.h>
#include <TParticle.h>
#include <TROOT.h>
#include <TVirtualMC.h>

// --- Standard library ---

// --- AliRoot header files ---
#include "AliEMCALv2.h"
#include "AliEMCALHit.h"
#include "AliEMCALGeometry.h"
#include "AliRun.h"
#include "AliHeader.h"
#include "AliMC.h"
#include "AliPoints.h"
// for TRD1 case only; May 31,2006

ClassImp(AliEMCALv2)

//______________________________________________________________________
AliEMCALv2::AliEMCALv2()
  : AliEMCALv1(), 
    fGeometry(0),
    fHDe(0),
    fHNhits(0)
{
  // ctor
}

//______________________________________________________________________
AliEMCALv2::AliEMCALv2(const char *name, const char *title)
  : AliEMCALv1(name,title),
    fGeometry(0),
    fHDe(0),
    fHNhits(0)
{
    // Standard Creator.

    fHits= new TClonesArray("AliEMCALHit",1000);
    gAlice->GetMCApp()->AddHitList(fHits);

    fNhits    = 0;
    fIshunt   = 2; // All hits are associated with particles entering the calorimeter
    fTimeCut  = 30e-09;

    fGeometry = GetGeometry(); 
    fHDe = fHNhits = 0;
    //    if (gDebug>0){
    if (1){
      TH1::AddDirectory(0);
      fHDe    = new TH1F("fHDe","De in EMCAL", 1000, 0., 10.);
      fHNhits = new TH1F("fHNhits","#hits in EMCAL", 2001, -0.5, 2000.5);
      fHistograms->Add(fHDe);
      fHistograms->Add(fHNhits);
      TH1::AddDirectory(1);
    }
}

//______________________________________________________________________
AliEMCALv2::~AliEMCALv2(){
    // dtor

    if ( fHits) {
	fHits->Delete();
	delete fHits;
	fHits = 0;
    }
}

//______________________________________________________________________
void AliEMCALv2::AddHit(Int_t shunt, Int_t primary, Int_t tracknumber, Int_t iparent, Float_t ienergy, 
			Int_t id, Float_t * hits,Float_t * p){
    // Add a hit to the hit list.
    // An EMCAL hit is the sum of all hits in a tower section
    //   originating from the same entering particle 
    static Int_t hitCounter;
    static AliEMCALHit *newHit, *curHit;
    static Bool_t deja;

    deja = kFALSE;

    newHit = new AliEMCALHit(shunt, primary, tracknumber, iparent, ienergy, id, hits, p);
    for ( hitCounter = fNhits-1; hitCounter >= 0 && !deja; hitCounter-- ) {
	curHit = (AliEMCALHit*) (*fHits)[hitCounter];
	// We add hits with the same tracknumber, while GEANT treats
	// primaries succesively
	if(curHit->GetPrimary() != primary) 
	  break;
	if( *curHit == *newHit ) {
	    *curHit = *curHit + *newHit;
	    deja = kTRUE;
	    //            break; // 30-aug-04 by PAI 
	} // end if
    } // end for hitCounter
    
    if ( !deja ) {
	new((*fHits)[fNhits]) AliEMCALHit(*newHit);
	fNhits++;
    }
    //    printf(" fNhits %i \n", fNhits); 
    delete newHit;
}
//______________________________________________________________________
void AliEMCALv2::StepManager(void){
  // Accumulates hits as long as the track stays in a tower

  // position wrt MRS and energy deposited
  static Float_t        xyzte[5]={0.,0.,0.,0.,0.};// position wrt MRS, time and energy deposited
  static Float_t        pmom[4]={0.,0.,0.,0.};
  static TLorentzVector pos;  // Lorentz vector of the track current position.
  static TLorentzVector mom;  // Lorentz vector of the track current momentum.
  static Float_t ienergy = 0; // part->Energy();
  static TString curVolName;
  static int supModuleNumber, moduleNumber, yNumber, xNumber, absid;
  static int keyGeom=1;
  static char *vn = "SCMX"; // Apr 13, 2006 - only TRD1 case now
  static int nSMOP[7]={1,3,5,7,9,11}; // 30-mar-05
  static int nSMON[7]={2,4,6,8,10,12};
  static Float_t depositedEnergy=0.0; 

  if(keyGeom == 0) {
    keyGeom = 2;
    if(gMC->VolId("PBMO")==0 || gMC->VolId("WSUC")==1) {
      vn      = "SCMX";   // old TRD2(TRD1) or WSUC
      keyGeom = 1;
    }    
    printf("AliEMCALv2::StepManager():  keyGeom %i : Sensetive volume %s \n", 
    keyGeom, vn); 
    if(gMC->VolId("WSUC")==1) printf(" WSUC - cosmic ray stand geometry \n");
  }
  Int_t tracknumber =  gAlice->GetMCApp()->GetCurrentTrackNumber();

  curVolName = gMC->CurrentVolName();
  if(curVolName.Contains(vn) || curVolName.Contains("SCX")) { // We are in a scintillator layer; SCX for 3X3
    
    if( ((depositedEnergy = gMC->Edep()) > 0.)  && (gMC->TrackTime() < fTimeCut)){// Track is inside a scintillator and deposits some energy
      //       Info("StepManager "," entry %i DE %f",++ientry, depositedEnergy); // for testing
       if (fCurPrimary==-1) 
	fCurPrimary=gAlice->GetMCApp()->GetPrimary(tracknumber);

      if (fCurParent==-1 || tracknumber != fCurTrack) {
	// Check parentage
	Int_t parent=tracknumber;
	if (fCurParent != -1) {
	  while (parent != fCurParent && parent != -1) {
	    TParticle *part=gAlice->GetMCApp()->Particle(parent);
	    parent=part->GetFirstMother();
	  }
	}
	if (fCurParent==-1 || parent==-1) {
	  Int_t parent=tracknumber;
	  TParticle *part=gAlice->GetMCApp()->Particle(parent);
	  while (parent != -1 && fGeometry->IsInEMCAL(part->Vx(),part->Vy(),part->Vz())) {
	    parent=part->GetFirstMother();
	    if (parent!=-1) 
	      part=gAlice->GetMCApp()->Particle(parent);
	  } 
	  fCurParent=parent;
	  if (fCurParent==-1)
	    Error("StepManager","Cannot find parent");
	  else {
	    TParticle *part=gAlice->GetMCApp()->Particle(fCurParent);
	    ienergy = part->Energy(); 
	  }
	  while (parent != -1) {
	    part=gAlice->GetMCApp()->Particle(parent);
	    part->SetBit(kKeepBit);
	    parent=part->GetFirstMother();
	  }
	}
	fCurTrack=tracknumber;
      }    
      gMC->TrackPosition(pos);
      xyzte[0] = pos[0];
      xyzte[1] = pos[1];
      xyzte[2] = pos[2];
      xyzte[3] = gMC->TrackTime() ;       
      
      gMC->TrackMomentum(mom);
      pmom[0] = mom[0];
      pmom[1] = mom[1];
      pmom[2] = mom[2];
      pmom[3] = mom[3];
      
      //      if(ientry%200 > 0) return; // testing
      supModuleNumber = moduleNumber = yNumber = xNumber = absid = 0;
      if(keyGeom >= 1) { // TRD1 case now
        gMC->CurrentVolOffID(4, supModuleNumber);
        gMC->CurrentVolOffID(3, moduleNumber);
        gMC->CurrentVolOffID(1, yNumber);
        gMC->CurrentVolOffID(0, xNumber); // really x number now
        if(strcmp(gMC->CurrentVolOffName(4),"SM10")==0) supModuleNumber += 10; // 13-oct-05
	// Nov 10,2006
        if(strcmp(gMC->CurrentVolOffName(0),vn) != 0) { // 3X3 case
          if     (strcmp(gMC->CurrentVolOffName(0),"SCX1")==0) xNumber=1;
          else if(strcmp(gMC->CurrentVolOffName(0),"SCX2")==0) xNumber=2;
          else if(strcmp(gMC->CurrentVolOffName(0),"SCX3")==0) xNumber=3;
          else Fatal("StepManager()", "Wrong name of sensetive volume in 3X3 case : %s ", gMC->CurrentVolOffName(0));
	}
      } else {
        gMC->CurrentVolOffID(5, supModuleNumber);
        gMC->CurrentVolOffID(4, moduleNumber);
        gMC->CurrentVolOffID(1, yNumber);
        gMC->CurrentVolOffID(0, xNumber);
        if     (strcmp(gMC->CurrentVolOffName(5),"SMOP")==0) supModuleNumber = nSMOP[supModuleNumber-1];
        else if(strcmp(gMC->CurrentVolOffName(5),"SMON")==0) supModuleNumber = nSMON[supModuleNumber-1];
        else   assert(0); // something wrong
      }
      absid = fGeometry->GetAbsCellId(supModuleNumber-1, moduleNumber-1, yNumber-1, xNumber-1);
    
      if (absid < 0) {
        printf(" supModuleNumber %i : moduleNumber %i : yNumber %i : xNumber %i \n",
        supModuleNumber, moduleNumber, yNumber, xNumber); 
	Fatal("StepManager()", "Wrong id : %i ", absid) ; 
      }

      Float_t lightYield =  depositedEnergy ;
      // Apply Birk's law (copied from G3BIRK)

      if (gMC->TrackCharge()!=0) { // Check
	  Float_t birkC1Mod = 0;
	if (fBirkC0==1){ // Apply correction for higher charge states
	  if (TMath::Abs(gMC->TrackCharge())>=2) birkC1Mod = fBirkC1*7.2/12.6;
	  else                                   birkC1Mod = fBirkC1;
	}

	Float_t dedxcm;
	if (gMC->TrackStep()>0)  dedxcm=1000.*gMC->Edep()/gMC->TrackStep();
	else                     dedxcm=0;
	lightYield=lightYield/(1.+birkC1Mod*dedxcm+fBirkC2*dedxcm*dedxcm);
      } 

      // use sampling fraction to get original energy --HG
      //      xyzte[4] = lightYield * fGeometry->GetSampling();
      xyzte[4] = lightYield; // 15-dec-04
        
      if (gDebug == -2) 
      printf("#sm %2i #m %3i #x %1i #z %1i -> absid %i : xyzte[4] = %f\n",
      supModuleNumber,moduleNumber,yNumber,xNumber,absid, xyzte[4]);

      AddHit(fIshunt, fCurPrimary,tracknumber, fCurParent, ienergy, absid,  xyzte, pmom);
    } // there is deposited energy
  }
}

void AliEMCALv2::FinishEvent()
{ 
  // Calculate deposit energy and fill control histogram; 26-may-05
  static double de=0.;
  fHNhits->Fill(double(fHits->GetEntries()));
  de = GetDepositEnergy(0);
  if(fHDe) fHDe->Fill(de);
}

Double_t AliEMCALv2::GetDepositEnergy(int print)
{ 
  // 23-mar-05 - for testing
  if(fHits == 0) return 0.;
  AliEMCALHit  *hit=0;
  Double_t de=0.;
  for(int ih=0; ih<fHits->GetEntries(); ih++) {
    hit = (AliEMCALHit*)fHits->UncheckedAt(ih);
    de += hit->GetEnergy();
  }
  if(print>0) {
    cout<<"AliEMCALv2::GetDepositEnergy() : fHits "<<fHits<<endl; 
    printf(" #hits %i de %f \n", fHits->GetEntries(), de);
    if(print>1) {
      printf(" #primary particles %i\n", gAlice->GetHeader()->GetNprimary()); 
    }
  }
  return de;
}

void AliEMCALv2::Browse(TBrowser* b)
{
  TObject::Browse(b);
}

void AliEMCALv2::DrawCalorimeterCut(const char *name, int axis, double dcut)
{ 
  // Size of tower is 5.6x5.6x24.8 (25.0); cut on Z axiz
  TString g(fGeometry->GetName());
  g.ToUpper();
  gMC->Gsatt("*", "seen", 0);

  int fill = 1;

  if(axis!=1) SetVolumeAttributes("STPL", 1, 1, fill);

  int colo=4;
  TString sn(name);
  if(sn.Contains("SCM")) colo=5;
  SetVolumeAttributes(name, 1, colo, fill);
  if(g.Contains("110DEG") && sn=="SMOD") SetVolumeAttributes("SM10", 1, colo, fill);

  TString st(GetTitle());
  st += ", zcut, ";
  st += name;

  char *optShad = "on", *optHide="on";
  double cxy=0.02;
  if     (axis==1) {
    dcut = 0.;
    //    optHide = optShad = "off";
  } else if(axis==2){
    if(dcut==0.) SetVolumeAttributes("SCM0", 1, 5, fill); // yellow    
  }

  gMC->Gdopt("hide", optHide);
  gMC->Gdopt("shad", optShad);

  gROOT->ProcessLine("TGeant3 *g3 = (TGeant3*)gMC");
  char cmd[200];
  sprintf(cmd,"g3->Gdrawc(\"XEN1\", %i, %5.1f, 12., 8., %3.2f, %3.2f)\n", axis, dcut, cxy,cxy);
  printf("\n%s\n",cmd); gROOT->ProcessLine(cmd);

  sprintf(cmd,"gMC->Gdhead(1111, \"%s\")\n", st.Data());
  printf("%s\n",cmd); gROOT->ProcessLine(cmd);
}

void AliEMCALv2::DrawSuperModuleCut(const char *name, int axis, double dcut, int fill)
{ 
 // Size of tower is 5.6x5.6x24.8 (25.0); cut on Z axiz
  TString sn(GetGeometry()->GetName());
  sn.ToUpper();
  char *tit[3]={"xcut", "ycut", "zcut"};
  if(axis<1) axis=1; if(axis>3) axis=3;

  gMC->Gsatt("*", "seen", 0);
  //  int fill = 6;

  // SetVolumeAttributes("SMOD", 1, 1, fill);  // black
  SetVolumeAttributes(name, 1, 5, fill);    // yellow 

  double cxy=0.055, x0=10., y0=10.;
  char *optShad = "on", *optHide="on";
  if(sn.Contains("TRD1")) {
    SetVolumeAttributes("STPL", 1, 3, fill);  // green 
    if     (axis==1) {
      gMC->Gsatt("STPL", "seen", 0);
      dcut = 0.;
      optHide = "off";
      optShad = "off";
    } else if(axis==3) cxy = 0.1;
  } else if(sn.Contains("TRD2")) {
    y0 = -10.;
    if (axis==2) cxy=0.06;
  }
  gMC->Gdopt("hide", optHide);
  gMC->Gdopt("shad", optShad);

  TString st("Shish-Kebab, Compact, SMOD, ");
  if(sn.Contains("TWIST")) st = "Shish-Kebab, Twist, SMOD, ";
  st += tit[axis-1];;

  gROOT->ProcessLine("TGeant3 *g3 = (TGeant3*)gMC");
  char cmd[200];
  sprintf(cmd,"g3->Gdrawc(\"SMOD\", %i, %6.3f, %5.1f, %5.1f, %5.3f, %5.3f)\n",axis,dcut,x0,y0,cxy,cxy);
  printf("\n%s\n",cmd); gROOT->ProcessLine(cmd);

  sprintf(cmd,"gMC->Gdhead(1111, \"%s\")\n", st.Data());
  printf("%s\n",cmd); gROOT->ProcessLine(cmd);
  // hint for testing
  if(sn.Contains("TRD1")){
    printf("Begin of super module\n");
    printf("g3->Gdrawc(\"SMOD\", 2,  0.300, 89., 10., 0.5, 0.5)\n");
    printf("Center of super modules\n");
    printf("g3->Gdrawc(\"SMOD\", 2,  0.300, 0., 10., 0.5, 0.5)\n");
    printf("end of super modules\n");
    printf("g3->Gdrawc(\"SMOD\", 2,  0.300, -70., 10., 0.5, 0.5)\n");
  } else if(sn.Contains("TRD2")){
    printf("Begin of super module  ** TRD2 ** \n");
    printf("g3->Gdrawc(\"SMOD\", 2,  0.00,  40., -80, 0.2, 0.2)\n");
    printf("end of super modules\n");
    printf("g3->Gdrawc(\"SMOD\", 2,  0.00, -20., -80, 0.2, 0.2)\n");

    printf(" ***  Z cut (Y|X plane)\n scale 0.4\n");
    printf("g3->Gdrawc(\"SMOD\", 3,  -170.,-165.,10.,0.4,0.4)\n");
    printf(" scale 0.2\n");
    printf("g3->Gdrawc(\"SMOD\", 3,  -170.,-80.,10.,0.2,0.2)\n");
    printf(" scale 0.12\n");
    printf("g3->Gdrawc(\"SMOD\", 3,   -170.,-45.,10.,0.12,0.12)\n");
  }
}

void AliEMCALv2::DrawTowerCut(const char *name, int axis, double dcut, int fill, char *optShad)
{ 
  // Size of tower is 5.6x5.6x24.8 (25.0); cut on Z axiz
  if(axis<1) axis=1; if(axis>3) axis=3;
  TString mn(name); mn.ToUpper();
  TString sn(GetGeometry()->GetName());

  gMC->Gsatt("*", "seen", 0);
  gMC->Gsatt("*", "fill", fill);

  // int mainColo=5; // yellow
  if(mn == "EMOD") {
    SetVolumeAttributes(mn.Data(), 1, 1, fill);
    SetVolumeAttributes("SCM0", 1, 5, fill); // yellow
    SetVolumeAttributes("SCPA", 1, 3, fill); // green - 13-may-05
  } else if(mn == "SCM0") { // first division 
    SetVolumeAttributes(mn.Data(), 1, 1, fill);
    SetVolumeAttributes("SCMY", 1, 5, fill); // yellow
  } else if(mn == "SCMY") { // first division 
    SetVolumeAttributes(mn.Data(), 1, 1, fill); 
    if(sn.Contains("TEST") && sn.Contains("3X3")) {
      SetVolumeAttributes("SCX1", 1, 5, fill); // yellow
      SetVolumeAttributes("SCX2", 1, 2, fill); // red
      SetVolumeAttributes("SCX3", 1, 3, fill); // green
    } else {
      SetVolumeAttributes("SCMX", 1, 5, fill); // yellow
    }
  } else if(mn == "SCMX" || mn.Contains("SCX")) {
    SetVolumeAttributes(mn.Data(), 1, 5, fill);// yellow
    SetVolumeAttributes("PBTI", 1, 4, fill);
  } else {
    printf("<W> for volume |%s| did not defined volume attributes\n", mn.Data());
  }

  //  TString st("Shish-Kebab, 2x2 mm sampling, 62 layers, ");
  TString st("Shashlyk, 2x2 mm sampling, 62 layers, ");
  if    (sn.Contains("25"))   st = "Shish-Kebab, 5x5 mm sampling, 25 layers, ";
  else if(sn.Contains("MAY05"))   st = "Shish-Kebab, 5x5 mm sampling, 77 layers, ";
  if(sn.Contains("TRD1")) st += " TRD1, ";
  if(sn.Contains("3X3"))  st += " 3x3, ";
  st += name;

  gROOT->ProcessLine("TGeant3 *g3 = (TGeant3*)gMC");
  double cx=0.78, cy=2.;
  if     (axis==1 && sn.Contains("TRD")==0) {
   cx = cy = 1.;
   gMC->Gsatt("PBTI", "seen", 0);
  } else if (sn.Contains("TEST") && sn.Contains("3X3")) {
    cx = cy = 0.7;
    if (axis==3) cx = cy = 1.;
  } else if (sn.Contains("TRD2")) {
    if (axis==3) cx = cy = 2.;
  } else if (mn.Contains("EMOD")) {
    cx = cy = 0.5;
  }
  char cmd[200];

  gMC->Gdopt("hide", optShad);
  gMC->Gdopt("shad", optShad);

  sprintf(cmd,"g3->Gdrawc(\"%s\", %i, %6.2f, 10., 10., %5.3f, %5.3f)\n",name, axis,dcut,cx,cy);
  printf("\n%s\n",cmd); gROOT->ProcessLine(cmd);

  sprintf(cmd,"gMC->Gdhead(1111, \"%s\")\n", st.Data());
  printf("%s\n",cmd); gROOT->ProcessLine(cmd);
}
  
void AliEMCALv2::DrawAlicWithHits(int mode)
{ 
 // 20-sep-04; does not work now
  static TH2F *h2;
  if(h2==0) h2 = new TH2F("h2","test fo hits", 60,0.5,60.5, 28,0.5,28.5);
  else      h2->Reset();
  if(mode==0) {
    gROOT->ProcessLine("TGeant3 *g3 = (TGeant3*)gMC");
    gMC->Gsatt("*","seen",0);
    gMC->Gsatt("scm0","seen",5);

    gROOT->ProcessLine("g3->Gdrawc(\"ALIC\", 1,   2.0, 12., -130, 0.3, 0.3)");
    // g3->Gdrawc("ALIC", 1,   2.0, 10., -14, 0.05, 0.05)
  }

  TClonesArray *hits = Hits();
  Int_t nhits = hits->GetEntries(), absId, nSupMod, nModule, nIphi, nIeta, iphi, ieta;
  AliEMCALHit *hit = 0;
  Double_t de, des=0.;
  for(Int_t i=0; i<nhits; i++) {
    hit   = (AliEMCALHit*)hits->UncheckedAt(i);
    absId = hit->GetId();
    de    = hit->GetEnergy();
    des += de;
    if(fGeometry->GetCellIndex(absId, nSupMod, nModule, nIphi, nIeta)){
      //      printf(" de %f abs id %i smod %i tower %i | cell iphi %i : ieta %i\n",
      // de, absId, nSupMod, nModule, nIphi, nIeta); 
      if(nSupMod==3) {
        fGeometry->GetCellPhiEtaIndexInSModule(nSupMod,nModule,nIphi,nIeta, iphi,ieta);
        // printf(" iphi %i : ieta %i\n", iphi,ieta);
        h2->Fill(double(ieta),double(iphi), de);
      }
    }
    //    hit->Dump();
  }
  printf(" #hits %i : sum de %f -> %f \n", nhits, des, h2->Integral());
  if(mode>0 && h2->Integral()>0.) h2->Draw();
}

void AliEMCALv2::SetVolumeAttributes(const char *name, int seen, int color, int fill)
{
 /* seen=-2:volume is visible but none of its descendants;
         -1:volume is not visible together with all its descendants;
          0:volume is not visible;
	  1:volume is     visible. */
  gMC->Gsatt(name, "seen", seen);
  gMC->Gsatt(name, "colo", color);
  gMC->Gsatt(name, "fill", fill); 
  printf(" %s : seen %i color %i fill %i \n", name, seen, color, fill);
} 

void AliEMCALv2::TestIndexTransition(int pri, int idmax)
{ 
 // Test for EMCAL_SHISH geometry
  TString sn(fGeometry->GetName());
  if(!sn.Contains("SHISH")) {
    printf("Wrong geometry |%s| ! Bye \n", sn.Data());
    return; 
  }

  Int_t nSupMod, nModule, nIphi, nIeta, idNew, nGood=0;
  if(idmax==0) idmax = fGeometry->GetNCells();
  for(Int_t id=1; id<=idmax; id++) {
    if(!fGeometry->GetCellIndex(id, nSupMod, nModule, nIphi, nIeta)){
      printf(" Wrong abs ID %i : #cells %i\n", id, fGeometry->GetNCells());
      break;  
    }
    idNew = fGeometry->GetAbsCellId(nSupMod, nModule, nIphi, nIeta);
    if(id != idNew || pri>0) {
      printf(" ID %i : %i <- new id\n", id, idNew);
      printf(" nSupMod   %i ",  nSupMod);
      printf(" nModule    %i ", nModule);
      printf(" nIphi     %i ", nIphi);
      printf(" nIeta     %i \n", nIeta);
     
    } else nGood++;
  }
  printf(" Good decoding %i : %i <- #cells \n", nGood, fGeometry->GetNCells());
}
