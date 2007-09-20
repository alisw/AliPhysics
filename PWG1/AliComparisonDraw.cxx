// The class definition in esdClus.h has been generated automatically
// by the ROOT utility TTree::MakeSelector(). This class is derived
// from the ROOT class TSelector. For more information on the TSelector
// framework see $ROOTSYS/README/README.SELECTOR or the ROOT User Manual.

// The following methods are defined in this file:
//    Begin():        called everytime a loop on the tree starts,
//                    a convenient place to create your histograms.
//    SlaveBegin():   called after Begin(), when on PROOF called only on the
//                    slave servers.
//    Process():      called for each event, in this function you decide what
//                    to read and fill your histograms.
//    SlaveTerminate: called at the end of the loop on the tree, when on PROOF
//                    called only on the slave servers.
//    Terminate():    called at the end of the loop on the tree,
//                    a convenient place to draw/fit your histograms.
//

//
// Comaprison draw
// Comapre the MC information with the reconstructed 
//



#include "TSystem.h"
#include <TPDGCode.h>
#include <TStyle.h>
#include "TCint.h"
#include "TH1I.h"
#include "TTimeStamp.h"
#include "TProof.h"
#include "TTree.h"
#include "TH2F.h"
#include "TH1F.h"

//
#include "AliTracker.h"
#include "AliMagF.h"
// 
#include "AliESDEvent.h"   // new container
#include "AliESD.h"
#include "AliESDtrack.h"
#include "AliESDfriend.h"
#include "AliESDfriendTrack.h"
#include "AliTPCseed.h"
#include "AliTPCclusterMI.h"
//
#include "AliComparisonDraw.h" 


ClassImp(AliComparisonDraw)



AliComparisonDraw::AliComparisonDraw():
  TObject(),
  fPtResolLPT(0),
  fPtResolHPT(0)
{
  InitHisto();
}

void AliComparisonDraw::InitHisto(){
  //
  //
  //
  fPtResolLPT = new TH2F("Pt resol","pt resol",10, 0.1,3,200,-0.2,0.2);
  fPtResolHPT = new TH2F("Pt resol","pt resol",10, 2,100,200,-0.3,0.3);  
  //
  fPtPoolLPT = new TH2F("Pt pool","pt pool",10, 0.1,3,200,-6,6);
  fPtPoolHPT = new TH2F("Pt pool","pt pool",10, 2,100,200,-6,6);  
}

void AliComparisonDraw::Process(AliMCInfo* infoMC, AliESDRecInfo *infoRC){
  //
  // 
  //
  Float_t mcpt = infoMC->GetParticle().Pt();

  //
  //
  if (infoRC->GetStatus(1)==0) return;
  if (!infoRC->GetESDtrack()) return;  //buggy line


  if (infoRC->GetESDtrack()->GetTPCNcls()<10) return;

  //  printf("Pt\t%f\t%f\n",mcpt, infoRC->GetESDtrack()->Pt());
  
  Float_t deltaPt= (mcpt-infoRC->GetESDtrack()->Pt())/mcpt;  
  Float_t poolPt= (1/mcpt-infoRC->GetESDtrack()->OneOverPt())/
    TMath::Sqrt(infoRC->GetESDtrack()->GetSigma1Pt2());  
  fPtResolLPT->Fill(mcpt,deltaPt);
  fPtResolHPT->Fill(mcpt,deltaPt);
  fPtPoolLPT->Fill(mcpt,poolPt);
  fPtPoolHPT->Fill(mcpt,poolPt);
  
}

