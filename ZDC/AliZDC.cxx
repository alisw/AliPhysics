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
Revision 1.28.6.2  2002/07/24 10:10:13  alibrary
Updating VirtualMC

Revision 1.28.6.1  2002/06/10 15:29:36  hristov
Merged with v3-08-02

Revision 1.30  2002/06/07 10:19:23  coppedis
TreeS, TreeD and TreeR for ZDC can be written in a separate file

Revision 1.29  2002/06/04 08:17:04  coppedis
Reconstruction method improved

Revision 1.28  2002/02/04 09:18:08  coppedis
Merging and reconstruction code review

Revision 1.26.2.2  2001/11/12 18:41:44  hristov
All the changes from the head are merged to the release

Revision 1.27  2001/10/21 18:27:45  hristov
Several pointers were set to zero in the default constructors to avoid memory management problems

Revision 1.26  2001/10/04 14:30:28  coppedis
Event merging for ZDC

Revision 1.25  2001/10/04 14:24:15  coppedis
Event merging for ZDC

Revision 1.24  2001/09/26 16:03:41  coppedis
Merging implemented

Revision 1.23  2001/05/15 13:44:57  coppedis
Changes in AddHit method

Revision 1.22  2001/05/14 09:53:32  coppedis
Adding functions ZMin and ZMax

Revision 1.21  2001/04/20 10:05:02  coppedis
Minor changes

Revision 1.20  2001/03/26 13:39:20  coppedis
Comment prints

Revision 1.19  2001/03/26 09:10:23  coppedis
Corrected bug in constructor (fIshunt has to be =1)

Revision 1.18  2001/03/20 08:21:55  coppedis
ZDC needs PIPE, ABSO, DIPO and SHIL

Revision 1.17  2001/03/16 16:18:03  coppedis
Correction for superposition of ZDC volumes with MUON arm one

Revision 1.16  2001/03/15 16:01:11  coppedis
Code review

Revision 1.15  2001/01/26 19:56:27  hristov
Major upgrade of AliRoot code

Revision 1.14  2000/12/12 13:17:01  coppedis
Minor corrections suggested by P. Hristov

Revision 1.12  2000/12/01 08:19:01  coppedis
Adding a message error if ZDC is constructed without DIPO

Revision 1.11  2000/11/30 17:21:03  coppedis
Introduce hit array fStHits reset only at the end of the event (for digitization)

Revision 1.10  2000/11/22 11:32:58  coppedis
Major code revision

Revision 1.9  2000/10/02 21:28:20  fca
Removal of useless dependecies via forward declarations

Revision 1.8  2000/07/10 13:58:01  fca
New version of ZDC from E.Scomparin & C.Oppedisano

Revision 1.7  2000/01/19 17:17:40  fca

Revision 1.6  1999/09/29 09:24:35  fca
Introduction of the Copyright and cvs Log

*/

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  			Zero Degree Calorimeter			             //
//  	     This class contains the basic functions for the ZDCs;           //
//            functions specific to one particular geometry are              //
//                      contained in the derived classes                     //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <stdlib.h>
#include <iostream.h>

// --- ROOT system
#include <TBRIK.h>
#include <TNode.h>
#include <TGeometry.h>
#include <TFile.h>
#include <TTree.h>
#include <TDirectory.h>
#include <TF1.h>

// --- AliRoot header files
#include "AliZDC.h"
#include "AliZDCHit.h"
#include "AliZDCMergedHit.h"
#include "AliZDCMerger.h"
#include "AliZDCDigit.h"
#include "AliZDCReco.h"
#include "AliDetector.h"
//#include "AliCallf77.h"
#include "AliConst.h"
#include "AliMC.h"
#include "AliRun.h"
#include "AliHeader.h"

 
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
  gAlice->AddHitList(fHits);
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
      Int_t primary = gAlice->GetPrimary(track);     
      if(track != primary){
        newquad->fSFlag = 1;  // SECONDARY particle entering the ZDC
      }
      else if(track == primary){
        newquad->fSFlag = 0;  // PRIMARY particle entering the ZDC
      }  
      sFlag 	= newquad->fSFlag;
      primKinEn = newquad->fPrimKinEn;
      xImpact 	= newquad->fXImpact;
      yImpact 	= newquad->fYImpact;
   }
   else{       
      newquad->fPrimKinEn = primKinEn;
      newquad->fXImpact	= xImpact;
      newquad->fYImpact = yImpact;
      newquad->fSFlag 	= sFlag;
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
  
  AliDetector::MakeBranch(opt);

  const char *cS = strstr(opt,"S");

  if (gAlice->TreeS() && cS) {
    if(fMergedHits!=0) fMergedHits->Clear();
    else fMergedHits = new TClonesArray ("AliZDCMergedHit",1000);
    MakeBranchInTree(gAlice->TreeS(), 
                     branchname, &fMergedHits, fBufferSize, file) ;
    printf("* AliZDC::MakeBranch    * Making Branch %s for SDigits\n\n",branchname);
  }

    
  const char *cD = strstr(opt,"D");

  if (gAlice->TreeD() && cD) {
    if(fDigits!=0) fDigits->Clear();
    else fDigits = new TClonesArray ("AliZDCDigit",1000);
    MakeBranchInTree(gAlice->TreeD(), 
                     branchname, &fDigits, fBufferSize, file) ;
    printf("* AliZDC::MakeBranch    * Making Branch %s for Digits\n\n",branchname);
  }

  
  const char *cR = strstr(opt,"R");

  if (gAlice->TreeR() && cR) {
    if(fRecPoints==0) fRecPoints = new TClonesArray("AliZDCReco",1000);
    MakeBranchInTree(gAlice->TreeR(), 
		     branchname, &fRecPoints, fBufferSize, file) ;
    printf("* AliZDC::MakeBranch    * Making Branch %s for RecPoints\n\n",branchname);   }
          
}

//_____________________________________________________________________________
 void AliZDC::MakeBranchInTreeS(TTree *treeS, const char *file)
{
  // MakeBranchInTree
  const Int_t kBufferSize = 4000;
  char  branchname[20];
  sprintf(branchname,"%s",GetName());
  MakeBranchInTree(treeS, branchname, &fMergedHits, kBufferSize, file) ;
  printf("* AliZDC::MakeBranch    * Making Branch %s for SDigits\n\n",branchname);

}
//_____________________________________________________________________________
 void AliZDC::MakeBranchInTreeD(TTree *treeD, const char *file)
{
  // MakeBranchInTree
  const Int_t kBufferSize = 4000;
  char  branchname[20];
  sprintf(branchname,"%s",GetName());
  MakeBranchInTree(treeD, branchname, &fDigits, kBufferSize, file) ;
  printf("* AliZDC::MakeBranch    * Making Branch %s for Digits\n\n",branchname);

}
//_____________________________________________________________________________
 void AliZDC::MakeBranchInTreeR(TTree *treeR, const char *file)
{
  // MakeBranchInTree
  const Int_t kBufferSize = 4000;
  char  branchname[20];
  sprintf(branchname,"%s",GetName());
  MakeBranchInTree(treeR, branchname, &fRecPoints, kBufferSize, file) ;
  printf("* AliZDC::MakeBranch    * Making Branch %s for RecPoints\n\n",branchname);

}
//_____________________________________________________________________________
void AliZDC::Hits2SDigits()
{
  printf("\n	Entering AliZDC::SDigits2Digits() ");
  
  //----------------------------------------------------------------
  if(!fMerger){ 
    printf("	ZDC digitization (without merging)\n");

    AliZDCMergedHit *MHit;
    Int_t j, sector[2];
    Float_t MHits[7];
    fNMergedhits = 0;

    TTree *treeH = gAlice->TreeH();
    Int_t ntracks = (Int_t) treeH->GetEntries();
    gAlice->ResetHits();
  
    // Tracks loop
    for(Int_t itrack=0; itrack<ntracks; itrack++){
       treeH->GetEvent(itrack);
       for(AliZDCHit* zdcHit=(AliZDCHit*)this->FirstHit(-1); zdcHit;
                      zdcHit = (AliZDCHit*)this->NextHit()){ 
		      
	   for(j=0; j<2; j++) sector[j] = zdcHit->GetVolume(j);
	   MHits[0] = zdcHit->GetPrimKinEn();
	   MHits[1] = zdcHit->GetXImpact();
	   MHits[2] = zdcHit->GetYImpact();
	   MHits[3] = zdcHit->GetSFlag();
	   MHits[4] = zdcHit->GetLightPMQ();
	   MHits[5] = zdcHit->GetLightPMC();
	   MHits[6] = zdcHit->GetEnergy();
       }//Hits loop
       
	  MHit = new AliZDCMergedHit(sector, MHits);
	  new((*fMergedHits)[fNMergedhits]) AliZDCMergedHit(*MHit);	  
	  TClonesArray &sdigits = *fMergedHits;
	  new (sdigits[fNMergedhits]) AliZDCMergedHit(*MHit);
	  fNMergedhits++;
	  delete MHit;
    }
    gAlice->TreeS()->Fill();
    gAlice->TreeS()->AutoSave(); 
    gAlice->TreeS()->Reset();  
  }
  //----------------------------------------------------------------
  else if(fMerger){
    printf("	ZDC merging and digitization\n");
    // ### Initialise merging
    fMerger -> InitMerging();

    // SDigits tree
    TTree *treeS = 0;
    char treeSName[20];
    sprintf(treeSName,"TreeS%d",fMerger->EvNum());
    if(gAlice->GetTreeSFile()){
      gAlice->GetTreeSFile()->cd();
      treeS = (TTree*)gAlice->GetTreeSFile()->Get(treeSName);
    }
    else {
      treeS = gAlice->TreeS();
    }
    if(!treeS){
      printf("\n ERROR -> Can't find TreeS%d in background file\n",fMerger->EvNum());
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
    AliZDCMergedHit *MHit;
    TClonesArray &sdigits = *fMergedHits;
    Int_t imhit;
    //Merged Hits loop
    for(imhit=0; imhit<fNMergedhits; imhit++){
       MHit = (AliZDCMergedHit*) fMergedHits->UncheckedAt(imhit);
       new (sdigits[imhit]) AliZDCMergedHit(*MHit);
    }
    treeS->Fill();
    treeS->AutoSave();
  }
  
}

//_____________________________________________________________________________
void AliZDC::SDigits2Digits()
{
  if(!fMerger){ // Only digitization
    printf("	ZDC digitization (no merging) \n");
    fMerger = new AliZDCMerger();    
    fMerger->Digitize(fNMergedhits, fMergedHits);

    char hname[30];
    sprintf(hname,"TreeD%d",gAlice->GetHeader()->GetEvent());
    gAlice->TreeD()->Fill();
    gAlice->TreeD()->AutoSave();
    gAlice->TreeD()->Reset();  
  }
  else if(fMerger){	// Merging and digitization
    printf("	ZDC merging and digitization\n");
    fMerger->Digitize(fNMergedhits, fMergedHits);

    // Digits tree
    TTree *treeD = 0;
    char treeDName[20];
    sprintf(treeDName,"TreeD%d",fMerger->EvNum());  
    if(gAlice->GetTreeDFile()){
      treeD = (TTree*)gAlice->GetTreeDFile()->Get(treeDName);
    }
    else {
      treeD = gAlice->TreeD();
    }
    if(!treeD){
      printf("\n ERROR -> Can't find TreeD%d in background file\n",fMerger->EvNum());
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
    printf("	Entering AliZDC::Digits2Reco\n");
    AliDetector *ZDC  = gAlice->GetDetector("ZDC");
    TClonesArray *ZDCdigits = ZDC->Digits();
    
    TTree *TD = 0;
    char treeDame[20];
    sprintf(treeDame,"TreeD%d",fMerger->EvNum());
    if(gAlice->GetTreeDFile()){
      TD = (TTree*)gAlice->GetTreeDFile()->Get(treeDame);
    }
    else {
      TD = gAlice->TreeD();
    }
    if(TD){
      char brname[20];
      sprintf(brname,"%s",ZDC->GetName());
      TBranch *br = TD->GetBranch(brname);
      if(br) br->SetAddress(&ZDCdigits);
    }
    else if(!TD) printf("	ERROR -> TreeD NOT found in gAlice object\n");
    
    Int_t nt = (Int_t) (TD->GetEntries());
    gAlice->ResetDigits();    
    
    AliZDCDigit *dig;
    Int_t j, idig, ndigits, ZNraw=0, ZPraw=0, ZEMraw=0;
    //	---	 Summing raw ADCs for each detector to obtain total light
    for(j=0; j<nt; j++){
      TD->GetEvent(j);
      ndigits = ZDCdigits->GetEntries();
      ZNraw=0;
      ZPraw=0; 
      ZEMraw=0;
      //  ---  Loop over event digits
      for(idig=0; idig<ndigits; idig++){
         dig = (AliZDCDigit*) ZDCdigits->UncheckedAt(idig);
         if(dig->GetSector(0) == 1)	 ZNraw  += dig->GetADCValue();
         else if(dig->GetSector(0) == 2) ZPraw  += dig->GetADCValue();
         else if(dig->GetSector(0) == 3) ZEMraw += dig->GetADCValue();
      } // Digits loop
    } //  TreeD entries loop
    printf("\n	---	ZNraw = %d, ZPraw = %d, ZEMraw = %d\n",ZNraw, ZPraw, ZEMraw);
    
  //  ---      Pedestal subtraction
  Int_t ZNcorr, ZPcorr, ZEMcorr, MeanPed=50;
  ZNcorr  = ZNraw  - 5*MeanPed;
  ZPcorr  = ZPraw  - 5*MeanPed;
  ZEMcorr = ZEMraw - 2*MeanPed;
  if(ZNcorr<0)  ZNcorr=0;
  if(ZPcorr<0)  ZPcorr=0;
  if(ZEMcorr<0) ZEMcorr=0;
 printf("\n    ZNcorr = %d, ZPcorr = %d, ZEMcorr = %d\n",ZNcorr,ZPcorr,ZEMcorr);
  
  //  ---      ADCchannel -> photoelectrons
  // NB-> PM gain = 10^(5), ADC resolution = 6.4*10^(-7)
  Float_t ZNphe, ZPphe, ZEMphe, ConvFactor = 0.064;
  ZNphe  = ZNcorr/ConvFactor;
  ZPphe  = ZPcorr/ConvFactor;
  ZEMphe = ZEMcorr/ConvFactor;
  printf("\n    ZNphe = %f, ZPphe = %f, ZEMphe = %f\n",ZNphe, ZPphe, ZEMphe);
  
  //  ---      Energy calibration
  // Conversion factors for hadronic ZDCs goes from phe yield to TRUE incident 
  //  energy (conversion from GeV to TeV is included); while for EM calos 
  // conversion is from light yield to detected energy calculated by GEANT
  // NB -> ZN and ZP conversion factors are constant since incident spectators
  // have all the same energy, ZEM energy is obtained through a fit over the whole
  // range of incident particle energies (obtained with full HIJING simulations) 
  Float_t ZNenergy, ZPenergy, ZEMenergy, ZDCenergy;
  Float_t ZNphexTeV=329., ZPphexTeV=369.;
  ZNenergy  = ZNphe/ZNphexTeV;
  ZPenergy  = ZPphe/ZPphexTeV;
  ZDCenergy = ZNenergy+ZPenergy;
  ZEMenergy = -4.81+0.3238*ZEMphe;
  if(ZEMenergy<0) ZEMenergy=0;
  printf("    ZNenergy = %f TeV, ZPenergy = %f TeV, ZDCenergy = %f GeV, "
         "\n		ZEMenergy = %f TeV\n", ZNenergy, ZPenergy, 
	 ZDCenergy, ZEMenergy);
  
  if(ZDCenergy==0)
    printf("\n\n	###	ATTENZIONE!!! -> ev# %d: ZNenergy = %f TeV, ZPenergy = %f TeV, ZDCenergy = %f GeV, "
         " ZEMenergy = %f TeV\n\n", fMerger->EvNum(), ZNenergy, ZPenergy, ZDCenergy, ZEMenergy); 
  
  //  ---      Number of incident spectator nucleons
  Int_t NDetSpecN, NDetSpecP;
  NDetSpecN = (Int_t) (ZNenergy/2.760);
  NDetSpecP = (Int_t) (ZPenergy/2.760);
  printf("\n    NDetSpecN = %d, NDetSpecP = %d\n",NDetSpecN, NDetSpecP);
  
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
  
  Int_t NGenSpecN=0, NGenSpecP=0, NGenSpec=0;
  Double_t ImpPar=0;
  // Cut value for Ezem (GeV)
  // [1] ### Results in Chiara's PhD thesis -> 0<b<15 fm (Dec 2001)
  //Float_t EZEMCut = 360.; 
  // [2] ### Results from a new production  -> 0<b<18 fm (Apr 2002)
  Float_t EZEMCut = 420.;
  Float_t DeltaEZEMSup = 690.; 
  Float_t DeltaEZEMInf = 270.; 
  if(ZEMenergy > (EZEMCut+DeltaEZEMSup)){
    NGenSpecN = (Int_t) (fZNCen->Eval(ZNenergy));
    NGenSpecP = (Int_t) (fZPCen->Eval(ZPenergy));
    NGenSpec  = (Int_t) (fZDCCen->Eval(ZDCenergy));
    ImpPar    = fbCen->Eval(ZDCenergy);
    //printf("    fZNCen = %f, fZPCen = %f, fZDCCen = %f\n",fZNCen->Eval(ZNenergy),
    //            fZPCen->Eval(ZPenergy),fZDCCen->Eval(ZDCenergy));
  }
  else if(ZEMenergy < (EZEMCut-DeltaEZEMInf)){
    NGenSpecN = (Int_t) (fZNPer->Eval(ZNenergy)); 
    NGenSpecP = (Int_t) (fZPPer->Eval(ZPenergy));
    NGenSpec  = (Int_t) (fZDCPer->Eval(ZDCenergy));
    ImpPar    = fbPer->Eval(ZDCenergy);
    //printf("    fZNPer = %f, fZPPer = %f, fZDCPer = %f\n",fZNPer->Eval(ZNenergy),
    //            fZPPer->Eval(ZPenergy),fZDCPer->Eval(ZDCenergy));
  }
  else if(ZEMenergy >= (EZEMCut-DeltaEZEMInf) && ZEMenergy <= (EZEMCut+DeltaEZEMSup)){
    NGenSpecN = (Int_t) (fZEMn->Eval(ZEMenergy));
    NGenSpecP = (Int_t) (fZEMp->Eval(ZEMenergy));
    NGenSpec  = (Int_t)(fZEMsp->Eval(ZEMenergy));
    ImpPar    =  fZEMb->Eval(ZEMenergy);
    //printf("    Nspec ZEM = %f, Nspec ZDC = %f\n",fZEMsp->Eval(ZNenergy),fZDCPer->Eval(ZDCenergy));
  }
  // [1] ### Results in Chiara's PhD thesis -> 0<b<15 fm (Dec 2001)
  /*if(ZNenergy>158.5)  NGenSpecN = (Int_t) (fZEMn->Eval(ZEMenergy));
  if(ZPenergy>58.91)  NGenSpecP = (Int_t) (fZEMp->Eval(ZEMenergy));
  if(ZDCenergy>220.4) NGenSpec  = (Int_t)(fZEMsp->Eval(ZEMenergy));
  if(ZDCenergy>225.)  ImpPar    =          fZEMb->Eval(ZEMenergy);*/
  // [2] ### Results from a new production  -> 0<b<18 fm (Apr 2002)
  if(ZNenergy>162.)  NGenSpecN = (Int_t) (fZEMn->Eval(ZEMenergy));
  if(ZPenergy>59.75)  NGenSpecP = (Int_t) (fZEMp->Eval(ZEMenergy));
  if(ZDCenergy>221.5) NGenSpec  = (Int_t)(fZEMsp->Eval(ZEMenergy));
  if(ZDCenergy>220.)  ImpPar    =  fZEMb->Eval(ZEMenergy);
  
  if(NGenSpecN>125)    NGenSpecN=125;
  else if(NGenSpecN<0) NGenSpecN=0;
  if(NGenSpecP>82)     NGenSpecP=82;
  else if(NGenSpecP<0) NGenSpecP=0;
  if(NGenSpec>207)     NGenSpec=207;
  else if(NGenSpec<0)  NGenSpec=0;
  //printf("    NRecSpecN = %d, NRecSpecP = %d, NRecSpec = %d\n",NGenSpecN,NGenSpecP,NGenSpec);
  
  //  ---      Number of participants
  Int_t NPart, NPartTot;
  NPart = 207-NGenSpecN-NGenSpecP;
  NPartTot = 207-NGenSpec;
  //printf("	###	NPart(ZP+ZN) = %d, NPart(ZDC) = %d, b = %f fm\n",NPart,NPartTot,ImpPar);
  printf("	###	NPart = %d, b = %f fm\n",NPartTot,ImpPar);
  
  //  ---     Writing RecPoints TCA
  // Allocate the RecPoints TCA 
  fRecPoints = new TClonesArray("AliZDCReco",1000);
  AliZDCReco *reco = new AliZDCReco(ZNenergy,ZPenergy,ZDCenergy,ZEMenergy,
  	      NDetSpecN,NDetSpecP,NGenSpecN,NGenSpecP,NGenSpec,NPartTot,ImpPar);
  new((*fRecPoints)[fNRecPoints]) AliZDCReco(*reco);
  //fNRecPoints++;
  //fRecPoints->Dump();
  delete reco;
  
  // TreeR
  TTree *treeR = gAlice->TreeR();
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

 
