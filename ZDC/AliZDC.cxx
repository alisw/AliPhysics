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

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  			Zero Degree Calorimeter			             //
//  	     This class contains the basic functions for the ZDCs;           //
//            functions specific to one particular geometry are              //
//                      contained in the derived classes                     //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <stdlib.h>
#include <Riostream.h>

// --- ROOT system
#include <TBRIK.h>
#include <TDirectory.h>
#include <TF1.h>
#include <TFile.h>
#include <TGeometry.h>
#include <TNode.h>
#include <TTree.h>
#include <TVirtualMC.h>

// --- AliRoot header files
#include "AliDetector.h"
#include "AliZDC.h"
#include "AliZDCDigit.h"
#include "AliZDCHit.h"
#include "AliZDCMergedHit.h"
#include "AliZDCMerger.h"
#include "AliZDCReco.h"

#include "AliConst.h"

#include "AliHeader.h"
#include "AliLoader.h"
#include "AliRun.h"
#include "AliMC.h"

 
ClassImp(AliZDC)
 
//_____________________________________________________________________________
AliZDC::AliZDC()
{
  //
  // Default constructor for the Zero Degree Calorimeter base class
  //
  
  fIshunt     = 1;
  fNoShower   = 0;

  fMerger     = 0;
  fHits       = 0;
  fNhits      = 0;

  fDigits     = 0;
  fNdigits    = 0;

  fMergedHits = 0;

  fNRecPoints = 0;
  fRecPoints  = 0;
  
}
 
//_____________________________________________________________________________
AliZDC::AliZDC(const char *name, const char *title)
  : AliDetector(name,title)
{
  //
  // Standard constructor for the Zero Degree Calorimeter base class
  //

  fIshunt   = 1;
  fNoShower = 0;
  fMerger   = 0;

  // Allocate the hits array  
  fHits   = new TClonesArray("AliZDCHit",1000);
  gAlice->GetMCApp()->AddHitList(fHits);
  // Allocate the merged hits array  
  fMergedHits = new TClonesArray("AliZDCMergedHit",1000);

  // Allocate the digits array  
  fDigits = new TClonesArray("AliZDCDigit",1000);
  
  fNRecPoints = 0;
  fRecPoints = 0;  

}
//____________________________________________________________________________ 
AliZDC::~AliZDC()
{
  //
  // ZDC destructor
  //

  fIshunt   = 0;
  
  if(fMerger) delete fMerger;

}
//_____________________________________________________________________________
void AliZDC::AddHit(Int_t track, Int_t *vol, Float_t *hits)
{
  //
  // 		Add a ZDC hit to the hit list.
  // -> We make use of 2 array of hits:
  // [1]  fHits (the usual one) that contains hits for each PRIMARY
  // [2]  fStHits that contains hits for each EVENT and is used to
  //	  obtain digits at the end of each event
  //
  
  static Float_t primKinEn, xImpact, yImpact, sFlag;

  AliZDCHit *newquad, *curprimquad;
  newquad = new AliZDCHit(fIshunt, track, vol, hits);
  TClonesArray &lhits = *fHits;
  
  if(fNhits==0){
      // First hit -> setting flag for primary or secondary particle
      Int_t primary = gAlice->GetMCApp()->GetPrimary(track);     
      if(track != primary){
        newquad->SetSFlag(1);  // SECONDARY particle entering the ZDC
      }
      else if(track == primary){
        newquad->SetSFlag(0);  // PRIMARY particle entering the ZDC
      }  
      sFlag 	= newquad->GetSFlag();
      primKinEn = newquad->GetPrimKinEn();
      xImpact 	= newquad->GetXImpact();
      yImpact 	= newquad->GetYImpact();
   }
   else{       
      newquad->SetPrimKinEn(primKinEn);
      newquad->SetXImpact(xImpact);
      newquad->SetYImpact(yImpact);
      newquad->SetSFlag(sFlag);
   }
 
  Int_t j;
  for(j=0; j<fNhits; j++){
    // If hits are equal (same track, same volume), sum them.
     curprimquad = (AliZDCHit*) lhits[j];
     if(*curprimquad == *newquad){
        *curprimquad = *curprimquad+*newquad;
	delete newquad;
	return;
     } 
  }

    //Otherwise create a new hit
    new(lhits[fNhits]) AliZDCHit(newquad);
    fNhits++;
    
    delete newquad;
}

//_____________________________________________________________________________
void  AliZDC::AddDigit(Int_t *sect, Int_t digit)
{
//
  AliZDCDigit *newdigit;
  newdigit = new AliZDCDigit(sect, digit);
  
//  printf("\n	AddDigit -> sector[0] = %d, sector[1] = %d, digit = %d",
//         sect[0], sect[1], digit);
  new((*fDigits)[fNdigits]) AliZDCDigit(*newdigit);
  fNdigits++;
  delete newdigit;
}
      
//_____________________________________________________________________________
void AliZDC::BuildGeometry()
{
  //
  // Build the ROOT TNode geometry for event display 
  // in the Zero Degree Calorimeter
  // This routine is dummy for the moment
  //

  TNode *node, *top;
  TBRIK *brik;
  const int kColorZDC  = kBlue;
  
  //
  top=gAlice->GetGeometry()->GetNode("alice");
  
  // ZDC
    brik = new TBRIK("S_ZDC","ZDC box","void",300,300,5);
    top->cd();
    node = new TNode("ZDC","ZDC","S_ZDC",0,0,600,"");
    node->SetLineColor(kColorZDC);
    fNodes->Add(node);
}

//_____________________________________________________________________________
Int_t AliZDC::DistancetoPrimitive(Int_t , Int_t )
{
  //
  // Distance from the mouse to the Zero Degree Calorimeter
  // Dummy routine
  //
  return 9999;
}

//____________________________________________________________________________
Float_t AliZDC::ZMin(void) const
{
  // Minimum dimension of the ZDC module in z
  return 11600.;
}

//____________________________________________________________________________
Float_t AliZDC::ZMax(void) const
{
  // Maximum dimension of the ZDC module in z
  return  11750.;
}
  

//_____________________________________________________________________________
 void AliZDC::MakeBranch(Option_t *opt, const char *file)
{
  //
  // Create Tree branches for the ZDC
  //

  char branchname[10];
  sprintf(branchname,"%s",GetName());

  const char *cH = strstr(opt,"H");
  
  if (cH && fLoader->TreeH())
   fHits   = new TClonesArray("AliZDCHit",1000); 
  
  AliDetector::MakeBranch(opt);

  const char *cS = strstr(opt,"S");

  if (fLoader->TreeS() && cS) {
    if(fMergedHits!=0) fMergedHits->Clear();
    else fMergedHits = new TClonesArray ("AliZDCMergedHit",1000);
    MakeBranchInTree(fLoader->TreeS(), 
                     branchname, &fMergedHits, fBufferSize, file) ;
    if (GetDebug()) printf("* AliZDC::MakeBranch    * Making Branch %s for SDigits\n\n",branchname);
  }

    
  const char *cD = strstr(opt,"D");

  if (fLoader->TreeD() && cD) {
    if(fDigits!=0) fDigits->Clear();
    else fDigits = new TClonesArray ("AliZDCDigit",1000);
    MakeBranchInTree(fLoader->TreeD(), 
                     branchname, &fDigits, fBufferSize, file) ;
    if (GetDebug()) printf("* AliZDC::MakeBranch    * Making Branch %s for Digits\n\n",branchname);
  }

  
  const char *cR = strstr(opt,"R");

  if (fLoader->TreeR() && cR) {
    if(fRecPoints==0) fRecPoints = new TClonesArray("AliZDCReco",1000);
    MakeBranchInTree(fLoader->TreeR(), 
		     branchname, &fRecPoints, fBufferSize, file) ;
    if (GetDebug()) printf("* AliZDC::MakeBranch    * Making Branch %s for RecPoints\n\n",branchname);   }
          
}

//_____________________________________________________________________________
 void AliZDC::MakeBranchInTreeS(TTree *treeS, const char *file)
{
  // MakeBranchInTree
  const Int_t kBufferSize = 4000;
  char  branchname[20];
  sprintf(branchname,"%s",GetName());
  if (fMergedHits==0x0) fMergedHits = new TClonesArray("AliZDCMergedHit",1000);
  MakeBranchInTree(treeS, branchname, &fMergedHits, kBufferSize, file) ;
  if (GetDebug()) printf("* AliZDC::MakeBranch    * Making Branch %s for SDigits\n\n",branchname);

}
//_____________________________________________________________________________
 void AliZDC::MakeBranchInTreeD(TTree *treeD, const char *file)
{
  // MakeBranchInTree
  const Int_t kBufferSize = 4000;
  char  branchname[20];
  sprintf(branchname,"%s",GetName());
  if (fDigits == 0x0) fDigits = new TClonesArray("AliZDCDigit",1000);
  MakeBranchInTree(treeD, branchname, &fDigits, kBufferSize, file) ;
  if (GetDebug()) printf("* AliZDC::MakeBranch    * Making Branch %s for Digits\n\n",branchname);

}
//_____________________________________________________________________________
 void AliZDC::MakeBranchInTreeR(TTree *treeR, const char *file)
{
  // MakeBranchInTree
  const Int_t kBufferSize = 4000;
  char  branchname[20];
  sprintf(branchname,"%s",GetName());
  MakeBranchInTree(treeR, branchname, &fRecPoints, kBufferSize, file) ;
  if (GetDebug()) printf("* AliZDC::MakeBranch    * Making Branch %s for RecPoints\n\n",branchname);

}
//_____________________________________________________________________________
void AliZDC::Hits2SDigits()
{
  if (GetDebug()) printf("\n	Entering AliZDC::SDigits2Digits() ");
  
  fLoader->LoadHits("read");
  fLoader->LoadSDigits("recreate");
  AliRunLoader* runLoader = fLoader->GetRunLoader(); 

  for (Int_t iEvent = 0; iEvent < runLoader->GetNumberOfEvents(); iEvent++) {
    runLoader->GetEvent(iEvent);
    if (!fLoader->TreeS()) fLoader->MakeTree("S");
    MakeBranch("S");

  //----------------------------------------------------------------
  if(!fMerger){ 
    if (GetDebug()) printf("	ZDC digitization (without merging)\n");

    AliZDCMergedHit *mHit;
    Int_t j, sector[2];
    Float_t mHits[7];
    fNMergedhits = 0;

    TTree *treeH = TreeH();
    Int_t ntracks = (Int_t) treeH->GetEntries();
    gAlice->ResetHits();
  
    // Tracks loop
    for(Int_t itrack=0; itrack<ntracks; itrack++){
       treeH->GetEvent(itrack);
       for(AliZDCHit* zdcHit=(AliZDCHit*)this->FirstHit(-1); zdcHit;
                      zdcHit = (AliZDCHit*)this->NextHit()){ 
		      
	   for(j=0; j<2; j++) sector[j] = zdcHit->GetVolume(j);
	   mHits[0] = zdcHit->GetPrimKinEn();
	   mHits[1] = zdcHit->GetXImpact();
	   mHits[2] = zdcHit->GetYImpact();
	   mHits[3] = zdcHit->GetSFlag();
	   mHits[4] = zdcHit->GetLightPMQ();
	   mHits[5] = zdcHit->GetLightPMC();
	   mHits[6] = zdcHit->GetEnergy();
       }//Hits loop
       
	  mHit = new AliZDCMergedHit(sector, mHits);
	  new((*fMergedHits)[fNMergedhits]) AliZDCMergedHit(*mHit);	  
	  TClonesArray &sdigits = *fMergedHits;
	  new (sdigits[fNMergedhits]) AliZDCMergedHit(*mHit);
	  fNMergedhits++;
	  delete mHit;
    }
      fLoader->TreeS()->Fill();
      fLoader->TreeS()->AutoSave(); 
      fLoader->TreeS()->Reset();  
  }
  //----------------------------------------------------------------
  else if(fMerger){
    if (GetDebug()) printf("	ZDC merging and digitization\n");
    // ### Initialise merging
    fMerger -> InitMerging();

    // SDigits tree



    TTree *treeS = fLoader->TreeS();
    if (treeS == 0x0)
     {
      Int_t retval = fLoader->LoadSDigits();
      if (retval)
       {
         Error("Hits2SDigits","Error while loading S. Digits");
         return;
       }
      treeS = fLoader->TreeS();
     }
    
    if(!treeS){
      if (GetDebug()) printf("\n ERROR -> Can't find TreeS%d in background file\n",fMerger->EvNum());
    }	 

    // ### Get TCA of MergedHits from AliZDCMerger
    fMergedHits  = fMerger->MergedHits();
    fNMergedhits = fMerger->GetNMhits();

    // Branch address
    char branchSDname[20];
    sprintf(branchSDname,"%s",GetName());
    if(treeS && fMergedHits){
      TBranch *branchSD = treeS->GetBranch(branchSDname);
      if(branchSD) branchSD->SetAddress(&fMergedHits);
      else if(!branchSD) MakeBranchInTreeS(treeS);
    }
    AliZDCMergedHit *mHit;
    TClonesArray &sdigits = *fMergedHits;
    Int_t imhit;
    //Merged Hits loop
    for(imhit=0; imhit<fNMergedhits; imhit++){
       mHit = (AliZDCMergedHit*) fMergedHits->UncheckedAt(imhit);
       new (sdigits[imhit]) AliZDCMergedHit(*mHit);
    }
    treeS->Fill();
    treeS->AutoSave();
  }
  
  }

  fLoader->UnloadHits();
  fLoader->UnloadSDigits();
}

//_____________________________________________________________________________
void AliZDC::SDigits2Digits()
{
  if(!fMerger){ // Only digitization
    if (GetDebug()) printf("	ZDC digitization (no merging) \n");
    fMerger = new AliZDCMerger();    
    fMerger->Digitize(fNMergedhits, fMergedHits);

    char hname[30];
    AliRunLoader * rl = fLoader->GetRunLoader();
    sprintf(hname,"TreeD%d",rl->GetHeader()->GetEvent());
    fLoader->TreeD()->Fill();
    fLoader->TreeD()->AutoSave();
    fLoader->TreeD()->Reset();  
  }
  else if(fMerger){	// Merging and digitization
    if (GetDebug()) printf("	ZDC merging and digitization\n");
    fMerger->Digitize(fNMergedhits, fMergedHits);

    // Digits tree

    TTree *treeD = fLoader->TreeD();
    if (treeD == 0x0)
     {
      Int_t retval = fLoader->LoadDigits();
      if (retval)
       {
         Error("SDigits2Digits","Error while loading Digits");
         return;
       }
      treeD = fLoader->TreeD();
     }



    if(!treeD){
      if (GetDebug()) printf("\n ERROR -> Can't find TreeD%d in background file\n",fMerger->EvNum());
    }	 
    // Branch address
    char branchDname[20];
    sprintf(branchDname,"%s",GetName());
    if(treeD && fDigits){
      TBranch *branchD = treeD->GetBranch(branchDname);
      if(branchD) branchD->SetAddress(&fDigits);
      else if(!branchD) MakeBranchInTreeD(treeD);
    }
    treeD->Fill();
    treeD->AutoSave();
  }
  
  
}
//_____________________________________________________________________________
void AliZDC::Hits2Digits()
{
    gAlice->Hits2SDigits();
    gAlice->SDigits2Digits();
}

//_____________________________________________________________________________
void AliZDC::Digits2Reco()
{
    if (GetDebug()) printf("	Entering AliZDC::Digits2Reco\n");
    AliDetector *zdcd  = gAlice->GetDetector("ZDC");
    TClonesArray *zdcdigits = zdcd->Digits();
    
    TTree *td = fLoader->TreeD();
    if (td == 0x0)
     {
      Int_t retval = fLoader->LoadDigits();
      if (retval)
       {
         Error("Digits2Reco","Error while loading Digits");
         return;
       }
      td = fLoader->TreeD();
     }


    if(td){
      char brname[20];
      sprintf(brname,"%s",zdcd->GetName());
      TBranch *br = td->GetBranch(brname);
      if(br) br->SetAddress(&zdcdigits);
    }
    else if(!td) printf("	ERROR -> TreeD NOT found in gAlice object\n");
    
    Int_t nt = (Int_t) (td->GetEntries());
    gAlice->ResetDigits();    
    
    AliZDCDigit *dig;
    Int_t j, idig, ndigits, znraw=0, zpraw=0, zemraw=0;
    //	---	 Summing raw ADCs for each detector to obtain total light
    for(j=0; j<nt; j++){
      td->GetEvent(j);
      ndigits = zdcdigits->GetEntries();
      znraw=0;
      zpraw=0; 
      zemraw=0;
      //  ---  Loop over event digits
      for(idig=0; idig<ndigits; idig++){
         dig = (AliZDCDigit*) zdcdigits->UncheckedAt(idig);
         if(dig->GetSector(0) == 1)	 znraw  += dig->GetADCValue();
         else if(dig->GetSector(0) == 2) zpraw  += dig->GetADCValue();
         else if(dig->GetSector(0) == 3) zemraw += dig->GetADCValue();
      } // Digits loop
    } //  TreeD entries loop
    if (GetDebug()) printf("\n	---	znraw = %d, zpraw = %d, zemraw = %d\n",znraw, zpraw, zemraw);
    
  //  ---      Pedestal subtraction
  Int_t zncorr, zpcorr, zemcorr, meanPed=50;
  zncorr  = znraw  - 5*meanPed;
  zpcorr  = zpraw  - 5*meanPed;
  zemcorr = zemraw - 2*meanPed;
  if(zncorr<0)  zncorr=0;
  if(zpcorr<0)  zpcorr=0;
  if(zemcorr<0) zemcorr=0;
 if (GetDebug()) printf("\n    zncorr = %d, zpcorr = %d, zemcorr = %d\n",zncorr,zpcorr,zemcorr);
  
  //  ---      ADCchannel -> photoelectrons
  // NB-> PM gain = 10^(5), ADC resolution = 6.4*10^(-7)
  Float_t znphe, zpphe, zemphe, convFactor = 0.064;
  znphe  = zncorr/convFactor;
  zpphe  = zpcorr/convFactor;
  zemphe = zemcorr/convFactor;
  if (GetDebug()) printf("\n    znphe = %f, zpphe = %f, zemphe = %f\n",znphe, zpphe, zemphe);
  
  //  ---      Energy calibration
  // Conversion factors for hadronic ZDCs goes from phe yield to TRUE incident 
  //  energy (conversion from GeV to TeV is included); while for EM calos 
  // conversion is from light yield to detected energy calculated by GEANT
  // NB -> ZN and ZP conversion factors are constant since incident spectators
  // have all the same energy, ZEM energy is obtained through a fit over the whole
  // range of incident particle energies (obtained with full HIJING simulations) 
  Float_t znenergy, zpenergy, zemenergy, zdcenergy;
  Float_t znphexTeV=329., zpphexTeV=369.;
  znenergy  = znphe/znphexTeV;
  zpenergy  = zpphe/zpphexTeV;
  zdcenergy = znenergy+zpenergy;
  zemenergy = -4.81+0.3238*zemphe;
  if(zemenergy<0) zemenergy=0;
  if (GetDebug()) printf("    znenergy = %f TeV, zpenergy = %f TeV, zdcenergy = %f GeV, "
         "\n		zemenergy = %f TeV\n", znenergy, zpenergy, 
	 zdcenergy, zemenergy);
  
  if(zdcenergy==0)
    if (GetDebug()) printf("\n\n	###	ATTENZIONE!!! -> ev# %d: znenergy = %f TeV, zpenergy = %f TeV, zdcenergy = %f GeV, "
         " zemenergy = %f TeV\n\n", fMerger->EvNum(), znenergy, zpenergy, zdcenergy, zemenergy); 
  
  //  ---      Number of incident spectator nucleons
  Int_t nDetSpecN, nDetSpecP;
  nDetSpecN = (Int_t) (znenergy/2.760);
  nDetSpecP = (Int_t) (zpenergy/2.760);
  if (GetDebug()) printf("\n    nDetSpecN = %d, nDetSpecP = %d\n",nDetSpecN, nDetSpecP);
  
  //  ---      Number of generated spectator nucleons and impact parameter
  // --------------------------------------------------------------------------------------------------
  // [1] ### Results in Chiara's PhD thesis -> 0<b<15 fm (Dec 2001)
  /*// Fit results for neutrons (Nspectator n true vs. EZN)
  TF1 *fZNCen = new TF1("fZNCen",
      "(-2.116909+sqrt(2.116909*2.116909-4*(-0.00651)*(14.556798-x)))/(2*(-0.00651))",0.,158.5);
  TF1 *fZNPer = new TF1("fZNPer",
      "(-34.695134-sqrt(34.695134*34.695134-4*(-0.174780)*(-1562.283443-x)))/(2*(-0.174780))",0.,158.5);
  // Fit results for protons (Nspectator p true vs. EZP)
  TF1 *fZPCen = new TF1("fZPCen",
      "(-1.3217+sqrt(1.3217*1.3217-4*(-0.007934)*(4.742873-x)))/(2*(-0.007934))",0.,58.91);
  TF1 *fZPPer = new TF1("fZPPer",
      "(-15.788267-sqrt(15.788267*15.788267-4*(-0.133359)*(-383.800673-x)))/(2*(-0.133359))",0.,58.91);
  // Fit results for total number of spectators (Nspectators true vs. EZDC)
  TF1 *fZDCCen = new TF1("fZDCCen",
      "(-1.867335+sqrt(1.867335*1.867335-4*(-0.004119)*(19.100289-x)))/(2*(-0.004119))",0.,220.4);
  TF1 *fZDCPer = new TF1("fZDCPer",
      "(-22.429097-sqrt(22.429097*22.429097-4*(-0.072435)*(-1482.034526-x)))/(2*(-0.072435))",0.,220.4);*/
  // --------------------------------------------------------------------------------------------------
  // [1] ### Results from a new production  -> 0<b<18 fm (Apr 2002)
  // Fit results for neutrons (Nspectator n true vs. EZN)
  TF1 *fZNCen = new TF1("fZNCen",
      "(-2.287920+sqrt(2.287920*2.287920-4*(-0.007629)*(11.921710-x)))/(2*(-0.007629))",0.,164.);
  TF1 *fZNPer = new TF1("fZNPer",
      "(-37.812280-sqrt(37.812280*37.812280-4*(-0.190932)*(-1709.249672-x)))/(2*(-0.190932))",0.,164.);
  // Fit results for protons (Nspectator p true vs. EZP)
  TF1 *fZPCen = new TF1("fZPCen",
       "(-1.321353+sqrt(1.321353*1.321353-4*(-0.007283)*(3.550697-x)))/(2*(-0.007283))",0.,60.);
  TF1 *fZPPer = new TF1("fZPPer",
      "(-42.643308-sqrt(42.643308*42.643308-4*(-0.310786)*(-1402.945615-x)))/(2*(-0.310786))",0.,60.);
  // Fit results for total number of spectators (Nspectators true vs. EZDC)
  TF1 *fZDCCen = new TF1("fZDCCen",
      "(-1.934991+sqrt(1.934991*1.934991-4*(-0.004080)*(15.111124-x)))/(2*(-0.004080))",0.,225.);
  TF1 *fZDCPer = new TF1("fZDCPer",
      "(-34.380639-sqrt(34.380639*34.380639-4*(-0.104251)*(-2612.189017-x)))/(2*(-0.104251))",0.,225.);
  // --------------------------------------------------------------------------------------------------
  // [1] ### Results in Chiara's PhD thesis -> 0<b<15 fm (Dec 2001)
  /*// Fit results for b (b vs. EZDC)
  //TF1 *fbCen = new TF1("fbCen","0.611543+0.052231*x-0.000112*x*x+0.000000374*x*x*x",0.,222.);
  //TF1 *fbPer = new TF1("fbPer","16.552010-0.023866*x-0.00001*x*x",0.,222.);
  TF1 *fbCen = new TF1("fbCen","0.612769+0.051929*x-0.0001074*x*x+0.0000003724*x*x*x",0.,225.);
  TF1 *fbPer = new TF1("fbPer","16.6131016-0.026053*x+0.000006893*x*x",0.,225.);*/
  // --------------------------------------------------------------------------------------------------
  // [2] ### Results from a new production  -> 0<b<18 fm (Apr 2002)
  TF1 *fbCen = new TF1("fbCen","-0.056923+0.079703*x-0.0004301*x*x+0.000001366*x*x*x",0.,220.);
  TF1 *fbPer = new TF1("fbPer","17.943998-0.046846*x+0.000074*x*x",0.,220.);
  // --------------------------------------------------------------------------------------------------
  // Evaluating Nspectators and b from ZEM energy
  // [1] ### Results in Chiara's PhD thesis -> 0<b<15 fm (Dec 2001)
  /*TF1 *fZEMn  = new TF1("fZEMn","124.2-0.0566*x+0.000006014*x*x",0.,3500.);
  TF1 *fZEMp  = new TF1("fZEMp","81.3-0.03834*x+0.000004359*x*x",0.,3500.);
  TF1 *fZEMsp = new TF1("fZEMsp","205.6-0.09567*x+0.00001056*x*x",0.,3500.);
  TF1 *fZEMb  = new TF1("fZEMb","15.8-0.02084*x+2.802e-5*x*x-2.007e-8*x*x*x+6.586e-12*x*x*x*x-8.042e-16*x*x*x*x*x",0.,3500.);*/
  // --------------------------------------------------------------------------------------------------
  // [2] ### Results from a new production  -> 0<b<18 fm (Apr 2002)
  TF1 *fZEMn  = new TF1("fZEMn","126.2-0.05399*x+0.000005679*x*x",0.,4000.);
  TF1 *fZEMp  = new TF1("fZEMp","82.49-0.03611*x+0.00000385*x*x",0.,4000.);
  TF1 *fZEMsp = new TF1("fZEMsp","208.7-0.09006*x+0.000009526*x*x",0.,4000.);
  TF1 *fZEMb  = new TF1("fZEMb","16.06-0.01633*x+1.44e-5*x*x-6.778e-9*x*x*x+1.438e-12*x*x*x*x-1.112e-16*x*x*x*x*x",0.,4000.);
  
  Int_t nGenSpecN=0, nGenSpecP=0, nGenSpec=0;
  Double_t impPar=0;
  // Cut value for Ezem (GeV)
  // [1] ### Results in Chiara's PhD thesis -> 0<b<15 fm (Dec 2001)
  //Float_t eZEMCut = 360.; 
  // [2] ### Results from a new production  -> 0<b<18 fm (Apr 2002)
  Float_t eZEMCut = 420.;
  Float_t deltaEZEMSup = 690.; 
  Float_t deltaEZEMInf = 270.; 
  if(zemenergy > (eZEMCut+deltaEZEMSup)){
    nGenSpecN = (Int_t) (fZNCen->Eval(znenergy));
    nGenSpecP = (Int_t) (fZPCen->Eval(zpenergy));
    nGenSpec  = (Int_t) (fZDCCen->Eval(zdcenergy));
    impPar    = fbCen->Eval(zdcenergy);
    //printf("    fZNCen = %f, fZPCen = %f, fZDCCen = %f\n",fZNCen->Eval(znenergy),
    //            fZPCen->Eval(zpenergy),fZDCCen->Eval(zdcenergy));
  }
  else if(zemenergy < (eZEMCut-deltaEZEMInf)){
    nGenSpecN = (Int_t) (fZNPer->Eval(znenergy)); 
    nGenSpecP = (Int_t) (fZPPer->Eval(zpenergy));
    nGenSpec  = (Int_t) (fZDCPer->Eval(zdcenergy));
    impPar    = fbPer->Eval(zdcenergy);
    //printf("    fZNPer = %f, fZPPer = %f, fZDCPer = %f\n",fZNPer->Eval(znenergy),
    //            fZPPer->Eval(zpenergy),fZDCPer->Eval(zdcenergy));
  }
  else if(zemenergy >= (eZEMCut-deltaEZEMInf) && zemenergy <= (eZEMCut+deltaEZEMSup)){
    nGenSpecN = (Int_t) (fZEMn->Eval(zemenergy));
    nGenSpecP = (Int_t) (fZEMp->Eval(zemenergy));
    nGenSpec  = (Int_t)(fZEMsp->Eval(zemenergy));
    impPar    =  fZEMb->Eval(zemenergy);
    //printf("    Nspec ZEM = %f, Nspec ZDC = %f\n",fZEMsp->Eval(znenergy),fZDCPer->Eval(zdcenergy));
  }
  // [1] ### Results in Chiara's PhD thesis -> 0<b<15 fm (Dec 2001)
  /*if(znenergy>158.5)  nGenSpecN = (Int_t) (fZEMn->Eval(zemenergy));
  if(zpenergy>58.91)  nGenSpecP = (Int_t) (fZEMp->Eval(zemenergy));
  if(zdcenergy>220.4) nGenSpec  = (Int_t)(fZEMsp->Eval(zemenergy));
  if(zdcenergy>225.)  impPar    =          fZEMb->Eval(zemenergy);*/
  // [2] ### Results from a new production  -> 0<b<18 fm (Apr 2002)
  if(znenergy>162.)  nGenSpecN = (Int_t) (fZEMn->Eval(zemenergy));
  if(zpenergy>59.75)  nGenSpecP = (Int_t) (fZEMp->Eval(zemenergy));
  if(zdcenergy>221.5) nGenSpec  = (Int_t)(fZEMsp->Eval(zemenergy));
  if(zdcenergy>220.)  impPar    =  fZEMb->Eval(zemenergy);
  
  if(nGenSpecN>125)    nGenSpecN=125;
  else if(nGenSpecN<0) nGenSpecN=0;
  if(nGenSpecP>82)     nGenSpecP=82;
  else if(nGenSpecP<0) nGenSpecP=0;
  if(nGenSpec>207)     nGenSpec=207;
  else if(nGenSpec<0)  nGenSpec=0;
  //printf("    NRecSpecN = %d, NRecSpecP = %d, NRecSpec = %d\n",nGenSpecN,nGenSpecP,nGenSpec);
  
  //  ---      Number of participants
  Int_t nPart, nPartTot;
  nPart = 207-nGenSpecN-nGenSpecP;
  nPartTot = 207-nGenSpec;
  //printf("	###	nPart(ZP+ZN) = %d, nPart(ZDC) = %d, b = %f fm\n",nPart,nPartTot,impPar);
  if (GetDebug()) printf("	###	nPart = %d, b = %f fm\n",nPartTot,impPar);
  
  //  ---     Writing RecPoints TCA
  // Allocate the RecPoints TCA 
  fRecPoints = new TClonesArray("AliZDCReco",1000);
  AliZDCReco *reco = new AliZDCReco(znenergy,zpenergy,zdcenergy,zemenergy,
  	      nDetSpecN,nDetSpecP,nGenSpecN,nGenSpecP,nGenSpec,nPartTot,impPar);
  new((*fRecPoints)[fNRecPoints]) AliZDCReco(*reco);
  //fNRecPoints++;
  //fRecPoints->Dump();
  delete reco;
  
  // TreeR
  TTree *treeR = fLoader->TreeR();
  if(!treeR) printf("\n ERROR -> Can't find TreeR%d in background file\n",fMerger->EvNum());
  // Branch address
  char branchRname[20];
  sprintf(branchRname,"%s",GetName());
  if(fRecPoints){
    TBranch *branchR = treeR->GetBranch(branchRname);
    if(branchR) branchR->SetAddress(&fRecPoints);
    else if(!branchR) MakeBranchInTreeR(treeR);
  }
  treeR->Fill();
  treeR->AutoSave();
  treeR->Reset();
}

//______________________________________________________________________
void AliZDC::SetTreeAddress(){
  // Set branch address for the Trees.
  // Inputs:
  //      none.
  // Outputs:
  //      none.
  // Return:
  //      none.
  if (fLoader->TreeH() && (fHits == 0x0))
    fHits   = new TClonesArray("AliZDCHit",1000);
      
  if (fLoader->TreeD() && (fDigits == 0x0))
    fDigits = new TClonesArray("AliZDCDigit",1000);
      
  AliDetector::SetTreeAddress();
}
 
