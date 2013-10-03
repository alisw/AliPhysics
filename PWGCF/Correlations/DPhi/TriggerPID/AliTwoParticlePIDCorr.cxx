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

#include "AliTwoParticlePIDCorr.h"
#include "AliUEHistograms.h"
#include "AliLog.h"
#include "AliVParticle.h"
#include "AliCFContainer.h"
#include "TFormula.h"

#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TCanvas.h"
#include "TList.h"
#include "TKey.h"
#include "TFile.h"

#include "AliCentrality.h"
#include "Riostream.h"

#include <TSpline.h>
#include <AliPID.h>
#include "AliESDpid.h"
#include "AliAODpidUtil.h"
#include <AliPIDResponse.h>
#include <AliITSPIDResponse.h>
#include <AliTPCPIDResponse.h>
#include <AliTRDPIDResponse.h>
#include <AliTOFPIDResponse.h>

#include <AliAnalysisManager.h>
#include <AliInputEventHandler.h>
#include "AliAODInputHandler.h"
#include "AliESDInputHandler.h"

#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliESDtrackCuts.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliCentrality.h"
#include "AliPIDResponse.h"

#include "AliVEvent.h"
#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "AliVTrack.h"

#include "THnSparse.h"

#include "AliAODMCHeader.h"
#include "AliAODMCParticle.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliMCParticle.h"
#include "AliStack.h"
#include "TParticle.h"
#include <TDatabasePDG.h>
#include <TParticlePDG.h>

#include "AliGenCocktailEventHeader.h"
#include "AliGenEventHeader.h"

#include "AliESDPmdTrack.h"
#include "AliEventPoolManager.h"
//#include "AliAnalysisUtils.h"
using namespace AliPIDNameSpace;
using namespace std;

ClassImp(AliTwoParticlePIDCorr)
ClassImp(LRCParticlePID)

//*********Debojit Sarkar******************
//________________________________________________________________________
AliTwoParticlePIDCorr::AliTwoParticlePIDCorr() // All data members should be initialised here
:AliAnalysisTaskSE(),
  fOutput(0),
  fhistcentrality(0),
  fEventCounter(0),
  fEtaSpectrasso(0),
  fphiSpectraasso(0),
  fEtaSpectraTrigall(0),
  fCentralityCorrelation(0x0),
  fHistoTPCdEdx(0x0),
  fHistoTOFbeta(0x0),
// fHistocentNSigmaTPC(0x0),
//fHistocentNSigmaTOF(0x0),
  fsame(0x0),
  fmix(0x0),
  fdeletasame(0),
  fdelphisame(0),
  fdeletamixed(0),
  fdelphimixed(0),
 fdeletamixedproton(0),
  fdelphimixedproton(0),
  fdeletamixedkaonpion(0),
  fdelphimixedkaonpion(0),
  fPoolMgr(0x0),
  fArrayMC(0),
  fAnalysisType("AOD"), 
  twoTrackEfficiencyCutValue(0.02),
//fControlConvResoncances(0),
  fPID(NULL),
 eventno(0),
  fPtTOFPID(.6),
  fRequestTOFPID(kTRUE),
  fPIDType(NSigmaTPCTOF),
  fNSigmaPID(3),
  fNSigmaPIDtrig1(4),
  fUseExclusiveNSigma(kFALSE),
  fRemoveTracksT0Fill(kFALSE),
fSelectCharge(0),
fTriggerSelectCharge(0),
fAssociatedSelectCharge(0),
fTriggerRestrictEta(-1),
fEtaOrdering(kFALSE),
fCutConversions(kFALSE),
fCutResonances(kFALSE),
fRejectResonanceDaughters(-1),
  fOnlyOneEtaSide(0),
  fPtOrder(kTRUE),
fInjectedSignals(kFALSE),
  fRemoveWeakDecays(kFALSE),
fRemoveDuplicates(kFALSE),
applyefficiency(kFALSE),
  fDCAXYCut(0)     

{
  
for(Int_t jj=0;jj<4;jj++)
    {
for(Int_t kk=0;kk<10;kk++)
   {
     //event no counting
     fEventno[jj][kk]=0;
     fEventnobaryon[jj][kk]=0;
     fEventnomeson[jj][kk]=0;
   }
    }

 
 for ( Int_t i = 0; i < 16; i++) { 
    fHistQA[i] = NULL;
  }

for(Int_t jj=0;jj<9;jj++)
    {
      fTPCTOFpion2d[jj]=0;
  }

for(Int_t jj=0;jj<5;jj++)
    {
      fhistoassopioncont[jj]=0;
      fhistoassokaoncont[jj]=0;
      fhistoassoprotoncont[jj]=0;
    }

for(Int_t jj=0;jj<4;jj++)
    {
      fhistotrigbaryoncont[jj]=0;
      fhistotrigmesoncont[jj]=0;
    }


for(Int_t kk=0;kk<4;kk++)
    {
for(Int_t jj=0;jj<9;jj++)
    {
//binning acording to pttrig only,no zvtx  dependance, but include centrality dependence
      //data
       fHistoNSigmaTPCpion[kk][jj]=0; 
       fHistoNSigmaTOFpion[kk][jj]=0;
       fHistoNSigmaTPCTOFpion[kk][jj]=0; 
       fhistopionnsigmaTPCMC[kk][jj]=0; 
       fhistopionnsigmaTOFMC[kk][jj]=0;
       fhistopionnsigmaTPCTOFMC[kk][jj]=0;
       fhistokaonnsigmaTPCMC[kk][jj]=0; 
       fhistokaonnsigmaTOFMC[kk][jj]=0;
       fhistokaonnsigmaTPCTOFMC[kk][jj]=0; 
       fhistoprotonnsigmaTPCMC[kk][jj]=0; 
       fhistoprotonnsigmaTOFMC[kk][jj]=0;
       fhistoprotonnsigmaTPCTOFMC[kk][jj]=0; 
       fhistoelectronnsigmaTPCMC[kk][jj]=0; 
       fhistoelectronnsigmaTOFMC[kk][jj]=0;
       fhistoelectronnsigmaTPCTOFMC[kk][jj]=0;  
       }
    }

for(Int_t jj=0;jj<4;jj++)
    {
for(Int_t kk=0;kk<10;kk++)
    {
      //for(Int_t ii=0;ii<4;ii++)
      // { //trigger particle counting with centrality,zvtx,pttrig binning
     fEtaSpectraTrig[jj][kk]=0;
     fEtaSpectraTrigbaryon[jj][kk]=0;
     fEtaSpectraTrigmeson[jj][kk]=0;
     // }
     }
     }

for(Int_t jj=0;jj<4;jj++)
    {
for(Int_t ii=0;ii<2;ii++) //associated pt binning
   {
for(Int_t kk=0;kk<10;kk++)
   {
     //for(Int_t pp=0;pp<4;pp++) //trigger pt binning ,do weighted average
       // {
    falltrigallasso[jj][ii][kk]=0x0;
    falltrigpionasso[jj][ii][kk]=0x0;
    falltrigkaonasso[jj][ii][kk]=0x0; 
    falltrigprotonasso[jj][ii][kk]=0x0;   
    fbaryontrigallasso[jj][ii][kk]=0x0; 
    fbaryontrigpionasso[jj][ii][kk]=0x0;
    fbaryontrigkaonasso[jj][ii][kk]=0x0;
    fbaryontrigprotonasso[jj][ii][kk]=0x0;
    fmesontrigallasso[jj][ii][kk]=0x0;
    fmesontrigpionasso[jj][ii][kk]=0x0;   
    fmesontrigkaonasso[jj][ii][kk]=0x0; 
    fmesontrigprotonasso[jj][ii][kk]=0x0;

    falltrigallassomix[jj][ii][kk]=0x0;
    falltrigpionassomix[jj][ii][kk]=0x0;
    falltrigkaonassomix[jj][ii][kk]=0x0; 
    falltrigprotonassomix[jj][ii][kk]=0x0;    
    fbaryontrigallassomix[jj][ii][kk]=0x0; 
    fbaryontrigpionassomix[jj][ii][kk]=0x0;
    fbaryontrigkaonassomix[jj][ii][kk]=0x0;
    fbaryontrigprotonassomix[jj][ii][kk]=0x0;
    fmesontrigallassomix[jj][ii][kk]=0x0;
    fmesontrigpionassomix[jj][ii][kk]=0x0;   
    fmesontrigkaonassomix[jj][ii][kk]=0x0; 
    fmesontrigprotonassomix[jj][ii][kk]=0x0;
    
    //}
    }
  }
}


 for ( Int_t i = 0; i < 6; i++ ){
    fTHnrecoallPid[i] = NULL;
    fTHngenprimPidTruth[i] = NULL;
    effcorection[i]=NULL;
    //effmap[i]=NULL;

  }



for(Int_t jj=0;jj<4;jj++)
  {//centrality binning
     recoallpt[jj]=0;
     recoalleta[jj]=0;
     alltrigeta[jj]=0;
     allassoeta[jj]=0;
     baryontrigeta[jj]=0;
     mesontrigeta[jj]=0;
     pionassoeta[jj]=0;
     kaonassoeta[jj]=0;
     protonassoeta[jj]=0;
     recoallphi[jj]=0;
     MCrecomatchedprimpt[jj]=0;
     MCrecomatchedprimeta[jj]=0;
     MCrecomatchedprimphi[jj]=0;
     MCtruthpt[jj]=0;
     MCtrutheta[jj]=0;
     MCtruthphi[jj]=0;

MCrecomatchedprimpionpt[jj]=0;
     MCrecomatchedprimpioneta[jj]=0;
     MCrecomatchedprimpionphi[jj]=0;

MCrecomatchedprimkaonpt[jj]=0;
     MCrecomatchedprimkaoneta[jj]=0;
     MCrecomatchedprimkaonphi[jj]=0;

MCrecomatchedprimprotonpt[jj]=0;
     MCrecomatchedprimprotoneta[jj]=0;
     MCrecomatchedprimprotonphi[jj]=0;

     MCtruthpionpt[jj]=0;
     MCtruthpioneta[jj]=0;
     MCtruthpionphi[jj]=0;

     MCtruthkaonpt[jj]=0;
     MCtruthkaoneta[jj]=0;
     MCtruthkaonphi[jj]=0;

     MCtruthprotonpt[jj]=0;
     MCtruthprotoneta[jj]=0;
     MCtruthprotonphi[jj]=0;
   }


  }
//________________________________________________________________________
AliTwoParticlePIDCorr::AliTwoParticlePIDCorr(const char *name) // All data members should be initialised here
  :AliAnalysisTaskSE(name),
 fOutput(0),
  fhistcentrality(0),
  fEventCounter(0),
  fEtaSpectrasso(0),
  fphiSpectraasso(0),
  fEtaSpectraTrigall(0),
  fCentralityCorrelation(0x0),
  fHistoTPCdEdx(0x0),
  fHistoTOFbeta(0x0),
// fHistocentNSigmaTPC(0x0),
//fHistocentNSigmaTOF(0x0),
  fsame(0x0),
  fmix(0x0),
  fdeletasame(0),
  fdelphisame(0),
  fdeletamixed(0),
  fdelphimixed(0),
 fdeletamixedproton(0),
  fdelphimixedproton(0),
  fdeletamixedkaonpion(0),
  fdelphimixedkaonpion(0),
   fPoolMgr(0x0),
  fArrayMC(0),
  fAnalysisType("AOD"), 
  twoTrackEfficiencyCutValue(0.02),
//fControlConvResoncances(0),
  fPID(NULL),
 eventno(0),
  fPtTOFPID(.6),
  fRequestTOFPID(kTRUE),
  fPIDType(NSigmaTPCTOF),
  fNSigmaPID(3),
  fNSigmaPIDtrig1(4),
  fUseExclusiveNSigma(kFALSE),
  fRemoveTracksT0Fill(kFALSE),
fSelectCharge(0),
fTriggerSelectCharge(0),
fAssociatedSelectCharge(0),
fTriggerRestrictEta(-1),
fEtaOrdering(kFALSE),
fCutConversions(kFALSE),
fCutResonances(kFALSE),
fRejectResonanceDaughters(-1),
  fOnlyOneEtaSide(0),
  fPtOrder(kTRUE),
fInjectedSignals(kFALSE),
  fRemoveWeakDecays(kFALSE),
fRemoveDuplicates(kFALSE),
applyefficiency(kFALSE),
  fDCAXYCut(0)     

    
{

for(Int_t jj=0;jj<4;jj++)
    {
for(Int_t kk=0;kk<10;kk++)
   {
     //event no counting
     fEventno[jj][kk]=0;
     fEventnobaryon[jj][kk]=0;
     fEventnomeson[jj][kk]=0;
   }
    }
  
   for ( Int_t i = 0; i < 16; i++) { 
    fHistQA[i] = NULL;
  }


for(Int_t jj=0;jj<9;jj++)
    {
      fTPCTOFpion2d[jj]=0;
  }

for(Int_t jj=0;jj<5;jj++)
    {
      fhistoassopioncont[jj]=0;
      fhistoassokaoncont[jj]=0;
      fhistoassoprotoncont[jj]=0;
    }

for(Int_t jj=0;jj<4;jj++)
    {
      fhistotrigbaryoncont[jj]=0;
      fhistotrigmesoncont[jj]=0;
    }

for(Int_t kk=0;kk<4;kk++)
    {
for(Int_t jj=0;jj<9;jj++)
    {
//binning acording to pttrig only,no zvtx  dependance, but include centrality dependence
      //data
       fHistoNSigmaTPCpion[kk][jj]=0; 
       fHistoNSigmaTOFpion[kk][jj]=0;
       fHistoNSigmaTPCTOFpion[kk][jj]=0; 
       fhistopionnsigmaTPCMC[kk][jj]=0; 
       fhistopionnsigmaTOFMC[kk][jj]=0;
       fhistopionnsigmaTPCTOFMC[kk][jj]=0;
       fhistokaonnsigmaTPCMC[kk][jj]=0; 
       fhistokaonnsigmaTOFMC[kk][jj]=0;
       fhistokaonnsigmaTPCTOFMC[kk][jj]=0; 
       fhistoprotonnsigmaTPCMC[kk][jj]=0; 
       fhistoprotonnsigmaTOFMC[kk][jj]=0;
       fhistoprotonnsigmaTPCTOFMC[kk][jj]=0; 
       fhistoelectronnsigmaTPCMC[kk][jj]=0; 
       fhistoelectronnsigmaTOFMC[kk][jj]=0;
       fhistoelectronnsigmaTPCTOFMC[kk][jj]=0;  
       }
    }

for(Int_t jj=0;jj<4;jj++)
    {
for(Int_t kk=0;kk<10;kk++)
   {
     //for(Int_t ii=0;ii<4;ii++)
     //{
     //trigger particle counting with centrality,zvtx,pttrig binning
     fEtaSpectraTrig[jj][kk]=0;
     fEtaSpectraTrigbaryon[jj][kk]=0;
     fEtaSpectraTrigmeson[jj][kk]=0;
     //}
   }
 }


for(Int_t jj=0;jj<4;jj++)
   {
for(Int_t ii=0;ii<2;ii++)//asso particle pt binning
   {
for(Int_t kk=0;kk<10;kk++)
   {
     //for(Int_t pp=0;pp<4;pp++) //trigger pt binning ,do weighted average
       //{
    falltrigallasso[jj][ii][kk]=0x0;
    falltrigpionasso[jj][ii][kk]=0x0;
    falltrigkaonasso[jj][ii][kk]=0x0; 
    falltrigprotonasso[jj][ii][kk]=0x0;     
    fbaryontrigallasso[jj][ii][kk]=0x0; 
    fbaryontrigpionasso[jj][ii][kk]=0x0;
    fbaryontrigkaonasso[jj][ii][kk]=0x0;
    fbaryontrigprotonasso[jj][ii][kk]=0x0;
    fmesontrigallasso[jj][ii][kk]=0x0;
    fmesontrigpionasso[jj][ii][kk]=0x0;   
    fmesontrigkaonasso[jj][ii][kk]=0x0; 
    fmesontrigprotonasso[jj][ii][kk]=0x0;

    falltrigallassomix[jj][ii][kk]=0x0;
    falltrigpionassomix[jj][ii][kk]=0x0;
    falltrigkaonassomix[jj][ii][kk]=0x0; 
    falltrigprotonassomix[jj][ii][kk]=0x0;   
    fbaryontrigallassomix[jj][ii][kk]=0x0; 
    fbaryontrigpionassomix[jj][ii][kk]=0x0;
    fbaryontrigkaonassomix[jj][ii][kk]=0x0;
    fbaryontrigprotonassomix[jj][ii][kk]=0x0;
    fmesontrigallassomix[jj][ii][kk]=0x0;
    fmesontrigpionassomix[jj][ii][kk]=0x0;   
    fmesontrigkaonassomix[jj][ii][kk]=0x0; 
    fmesontrigprotonassomix[jj][ii][kk]=0x0;
     //}
    }
  }
}


for ( Int_t i = 0; i < 6; i++ ){
    fTHnrecoallPid[i] = NULL;
    fTHngenprimPidTruth[i] = NULL;
    effcorection[i]=NULL;
    //effmap[i]=NULL;

  }

for(Int_t jj=0;jj<4;jj++)
  {//centrality binning
     recoallpt[jj]=0;
     recoalleta[jj]=0;
     alltrigeta[jj]=0;
     allassoeta[jj]=0;
     baryontrigeta[jj]=0;
     mesontrigeta[jj]=0;
     pionassoeta[jj]=0;
     kaonassoeta[jj]=0;
     protonassoeta[jj]=0;
     recoallphi[jj]=0;
     MCrecomatchedprimpt[jj]=0;
     MCrecomatchedprimeta[jj]=0;
     MCrecomatchedprimphi[jj]=0;
     MCtruthpt[jj]=0;
     MCtrutheta[jj]=0;
     MCtruthphi[jj]=0;

MCrecomatchedprimpionpt[jj]=0;
     MCrecomatchedprimpioneta[jj]=0;
     MCrecomatchedprimpionphi[jj]=0;

MCrecomatchedprimkaonpt[jj]=0;
     MCrecomatchedprimkaoneta[jj]=0;
     MCrecomatchedprimkaonphi[jj]=0;

MCrecomatchedprimprotonpt[jj]=0;
     MCrecomatchedprimprotoneta[jj]=0;
     MCrecomatchedprimprotonphi[jj]=0;

     MCtruthpionpt[jj]=0;
     MCtruthpioneta[jj]=0;
     MCtruthpionphi[jj]=0;

     MCtruthkaonpt[jj]=0;
     MCtruthkaoneta[jj]=0;
     MCtruthkaonphi[jj]=0;

     MCtruthprotonpt[jj]=0;
     MCtruthprotoneta[jj]=0;
     MCtruthprotonphi[jj]=0;

   }
  // The last in the above list should not have a comma after it
  
  // Constructor
  // Define input and output slots here (never in the dummy constructor)
  // Input slot #0 works with a TChain - it is connected to the default input container
  // Output slot #1 writes into a TH1 container
 
  DefineOutput(1, TList::Class());                                        // for output list

}

//________________________________________________________________________
AliTwoParticlePIDCorr::~AliTwoParticlePIDCorr()
{
  // Destructor. Clean-up the output list, but not the histograms that are put inside
  // (the list is owner and will clean-up these histograms). Protect in PROOF case.
  if (fOutput && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()) {
    delete fOutput;

  }
  if (fPID) delete fPID;
   
  }
//________________________________________________________________________
Float_t AliTwoParticlePIDCorr::PhiRange(Float_t DPhi)

{
	//
	// Puts the argument in the range [-pi/2,3 pi/2].
	//
	
	if (DPhi < -TMath::Pi()/2) DPhi += 2*TMath::Pi();
	if (DPhi > 3*TMath::Pi()/2) DPhi -= 2*TMath::Pi();	

	return DPhi;
	
}
//________________________________________________________________________
void AliTwoParticlePIDCorr::UserCreateOutputObjects()
{
  // Create histograms
  // Called once (on the worker node)
  AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
  fPID = inputHandler->GetPIDResponse();

  //AliAnalysisUtils *fUtils = new AliAnalysisUtils();

//get the efficiency correction map


  fOutput = new TList();
  fOutput->SetOwner();  // IMPORTANT!  

fhistcentrality=new TH1F("fhistcentrality","fhistcentrality",100,0.,100.);
fOutput->Add(fhistcentrality);

  fEventCounter = new TH1F("fEventCounter","EventCounter", 10, 0.5,10.5);
  fEventCounter->GetXaxis()->SetBinLabel(1,"Event Accesed");
  fEventCounter->GetXaxis()->SetBinLabel(2,"Within 0-100% centrality");
  fEventCounter->GetXaxis()->SetBinLabel(5,"Have a vertex");
  fEventCounter->GetXaxis()->SetBinLabel(6,"After vertex Cut");
  fEventCounter->GetXaxis()->SetBinLabel(7,"Event Analyzed");
  //fEventCounter->GetXaxis()->SetBinLabel(8,"Event Analysis finished");
  fOutput->Add(fEventCounter);
  
fEtaSpectrasso=new TH2F("fEtaSpectraasso","fEtaSpectraasso",180,-0.9,0.9,10,0,5 );
fOutput->Add(fEtaSpectrasso);

fphiSpectraasso=new TH2F("fphiSpectraasso","fphiSpectraasso",72,0,2*TMath::Pi(),10,0,5);
fOutput->Add(fphiSpectraasso);

fEtaSpectraTrigall=new TH1F("fEtaSpectraTrigall","fEtaSpectraTrigall",180,-0.9,0.9);
fOutput->Add(fEtaSpectraTrigall);

 fCentralityCorrelation = new TH2F("fCentralityCorrelation", ";centrality;multiplicity", 101, 0, 101, 200, 0, 4000);
      fOutput->Add(fCentralityCorrelation);


fHistoTPCdEdx = new TH2F("hHistoTPCdEdx", ";p_{T} (GeV/c);dE/dx (au.)",70, 0., 7., 500, 0., 500.);
fOutput->Add(fHistoTPCdEdx);
fHistoTOFbeta = new TH2F(Form("hHistoTOFbeta"), ";p_{T} (GeV/c);v/c",70, 0., 7., 500, 0.1, 1.1);
  fOutput->Add(fHistoTOFbeta);
  
  fsame=new TH2F ("fsame","#delta#eta-#delta#phi corr",72,-TMath::Pi()/2,3*TMath::Pi()/2,36,-1.8,1.8);
fOutput->Add(fsame);

 fmix=new TH2F ("fmix","#delta#eta-#delta#phi corr",72,-TMath::Pi()/2,3*TMath::Pi()/2,36,-1.8,1.8);
fOutput->Add(fmix);

fdeletasame=new TH1F("fdeletasame","#delta#eta mixed event distribution",36,-1.8,1.8);
fOutput->Add(fdeletasame);

fdelphisame=new TH1F("fdelphisame","#delta#phi mixed event distribution",72,-TMath::Pi()/2,3*TMath::Pi()/2);
fOutput->Add(fdelphisame);

fdeletamixed=new TH1F("fdeletamixed","#delta#eta mixed event distribution",36,-1.8,1.8);
fOutput->Add(fdeletamixed);

fdelphimixed=new TH1F("fdelphimixed","#delta#phi mixed event distribution",72,-TMath::Pi()/2,3*TMath::Pi()/2);
fOutput->Add(fdelphimixed);

fdeletamixedproton=new TH1F("fdeletamixedproton","#delta#eta mixed event distribution",36,-1.8,1.8);
fOutput->Add(fdeletamixedproton);

fdelphimixedproton=new TH1F("fdelphimixedproton","#delta#phi mixed event distribution",72,-TMath::Pi()/2,3*TMath::Pi()/2);
fOutput->Add(fdelphimixedproton);

fdeletamixedkaonpion=new TH1F("fdeletamixedkaonpion","#delta#eta mixed event distribution",36,-1.8,1.8);
fOutput->Add(fdeletamixedkaonpion);

fdelphimixedkaonpion=new TH1F("fdelphimixedkaonpion","#delta#phi mixed event distribution",72,-TMath::Pi()/2,3*TMath::Pi()/2);
fOutput->Add(fdelphimixedkaonpion);

TString Histpname;
for(Int_t jj=0;jj<4;jj++)
 {
for(Int_t kk=0;kk<10;kk++)
   {
  Histpname="fEventno";Histpname+=jj;Histpname+=kk;
  fEventno[jj][kk]=new TH1F (Histpname.Data(),"eventno",50,0.5,50.5);
  fOutput->Add(fEventno[jj][kk]);

 Histpname="fEventnobaryon";Histpname+=jj;Histpname+=kk;
  fEventnobaryon[jj][kk]=new TH1F (Histpname.Data(),"eventnobaryon",50,0.5,50.5);
  fOutput->Add(fEventnobaryon[jj][kk]);

 Histpname="fEventnomeson";Histpname+=jj;Histpname+=kk;
  fEventnomeson[jj][kk]=new TH1F (Histpname.Data(),"eventnomeson",50,0.5,50.5);
  fOutput->Add(fEventnomeson[jj][kk]);
   }
 }
 
  fHistQA[0] = new TH1F("fHistQAvx", "Histo Vx All ", 50, -5., 5.);
  fHistQA[1] = new TH1F("fHistQAvy", "Histo Vy All", 50, -5., 5.);
  fHistQA[2] = new TH1F("fHistQAvz", "Histo Vz All", 50, -25., 25.);  
  fHistQA[3] = new TH1F("fHistQAvxA", "Histo Vx  After Cut ", 50, -5., 5.);
  fHistQA[4] = new TH1F("fHistQAvyA", "Histo Vy After Cut", 50, -5., 5.);
  fHistQA[5] = new TH1F("fHistQAvzA", "Histo Vz After Cut", 50, -25., 25.);
  fHistQA[6] = new TH1F("fHistQADcaXyC", "Histo DCAxy after cut", 50, -5., 5.);
  fHistQA[7] = new TH1F("fHistQADcaZC", "Histo DCAz after cut", 50, -5., 5.);   
  fHistQA[8] = new TH1F("fHistQAPt","p_{T} distribution",900,0.,9.);
  fHistQA[9] = new TH1F("fHistQAEta","#eta distribution",360,-1.8,1.8);
  fHistQA[10] = new TH1F("fHistQAPhi","#phi distribution",340,0,6.8);
  fHistQA[11] = new TH1F("fHistQANCls","Number of TPC cluster",200,0,200);
  fHistQA[13] = new TH1F("fHistQAChi2","Chi2 per NDF",100,0,10);
 fHistQA[12] = new TH1F("fHistQANCls1","Number of TPC cluster1",200,0,200);
 fHistQA[14] = new TH1F("nCrossedRowsTPC","Number of TPC ccrossed rows",200,0,200);
 fHistQA[15] = new TH1F("ratioCrossedRowsOverFindableClustersTPC","Number of TPC ccrossed rows find clusters",200,0,2);



for(Int_t i = 0; i < 16; i++)
    {
      fOutput->Add(fHistQA[i]);
    }

TString Histpvrkname;
for(Int_t jj=0;jj<9;jj++)
  {//pt binning(asso)0,1,2,3,4,5,6,7,8
Histpvrkname="fTPCTOFpion2d";Histpvrkname+=jj;//Histtname+=ii;//Histtname+=kk;
 fTPCTOFpion2d[jj]=new TH2F (Histpvrkname.Data(),"fTPCTOFpion2d",1600,-80.,80.,1600,-80.,80);
   fOutput->Add(fTPCTOFpion2d[jj]);
  }

TString Histpvrname;
for(Int_t jj=0;jj<5;jj++)
  {//pt binning(asso) 4,5,6,7,8
    
Histpvrname="fhistoassoprotoncont";Histpvrname+=jj;//Histtname+=ii;//Histtname+=kk;
   fhistoassoprotoncont[jj]=new TH1F (Histpvrname.Data(),"fhistoassoprotoncont",50,0.,50.);
   fOutput->Add(fhistoassoprotoncont[jj]);
    
   Histpvrname="fhistoassopioncont";Histpvrname+=jj;//Histtname+=ii;//Histtname+=kk;
   fhistoassopioncont[jj]=new TH1F (Histpvrname.Data(),"fhistoassopioncont",50,0.,50.);
   fOutput->Add(fhistoassopioncont[jj]);

   Histpvrname="fhistoassokaoncont";Histpvrname+=jj;//Histtname+=ii;//Histtname+=kk;
   fhistoassokaoncont[jj]=new TH1F (Histpvrname.Data(),"fhistoassokaoncont",50,0.,50.);
   fOutput->Add(fhistoassokaoncont[jj]);

  }

TString Histpvname;
for(Int_t jj=0;jj<4;jj++)
  {//pt binning(trigger)0,1,2,3
   
   Histpvname="fhistotrigbaryoncont";Histpvname+=jj;//Histtname+=ii;//Histtname+=kk;
   fhistotrigbaryoncont[jj]=new TH1F (Histpvname.Data(),"fhistotrigbaryoncont",50,0.,50.);
   fOutput->Add(fhistotrigbaryoncont[jj]);

   Histpvname="fhistotrigmesoncont";Histpvname+=jj;//Histtname+=ii;//Histtname+=kk;
   fhistotrigmesoncont[jj]=new TH1F (Histpvname.Data(),"fhistotrigmesoncont",50,0.,50.);
   fOutput->Add(fhistotrigmesoncont[jj]);

 }



TString Histtname;
for(Int_t kk=0;kk<4;kk++)
    {
for(Int_t jj=0;jj<9;jj++)
    {
//data
   Histtname="fHistoNSigmaTPCpion";Histtname+=kk;Histtname+=jj;//Histtname+=ii;//Histtname+=kk;
   fHistoNSigmaTPCpion[kk][jj]=new TH1F (Histtname.Data(),"fHistoNSigmaTPCpion",1200, -60., 60.);
   fOutput->Add(fHistoNSigmaTPCpion[kk][jj]);

   Histtname="fHistoNSigmaTOFpion";Histtname+=kk;Histtname+=jj;//Histtname+=ii;//Histtname+=kk;
   fHistoNSigmaTOFpion[kk][jj]=new TH1F (Histtname.Data(),"fHistoNSigmaTOFpion",1200, -60., 60.);
   fOutput->Add(fHistoNSigmaTOFpion[kk][jj]);

 Histtname="fHistoNSigmaTPCTOFpion";Histtname+=kk;Histtname+=jj;//Histtname+=ii;//Histtname+=kk;
   fHistoNSigmaTPCTOFpion[kk][jj]=new TH1F (Histtname.Data(),"fHistoNSigmaTPCTOFpion",1200, -60., 60.);
   fOutput->Add(fHistoNSigmaTPCTOFpion[kk][jj]);

 Histtname="fhistopionnsigmaTPCMC";Histtname+=kk;Histtname+=jj;//Histtname+=ii;//Histtname+=kk;
   fhistopionnsigmaTPCMC[kk][jj]=new TH1F (Histtname.Data(),"fhistopionnsigmaTPCMC",1200, -60., 60.);
   fOutput->Add(fhistopionnsigmaTPCMC[kk][jj]);

 Histtname="fhistopionnsigmaTOFMC";Histtname+=kk;Histtname+=jj;//Histtname+=ii;//Histtname+=kk;
   fhistopionnsigmaTOFMC[kk][jj]=new TH1F (Histtname.Data(),"fhistopionnsigmaTOFMC",1200, -60., 60.);
   fOutput->Add(fhistopionnsigmaTOFMC[kk][jj]);

 Histtname="fhistopionnsigmaTPCTOFMC";Histtname+=kk;Histtname+=jj;//Histtname+=ii;//Histtname+=kk;
 fhistopionnsigmaTPCTOFMC[kk][jj]=new TH1F (Histtname.Data(),"fhistopionnsigmaTPCTOFMC",1200, -60., 60.);
 fOutput->Add(fhistopionnsigmaTPCTOFMC[kk][jj]);
  
 Histtname="fhistokaonnsigmaTPCMC";Histtname+=kk;Histtname+=jj;//Histtname+=ii;//Histtname+=kk;
   fhistokaonnsigmaTPCMC[kk][jj]=new TH1F (Histtname.Data(),"fhistokaonnsigmaTPCMC",1200, -60., 60.);
   fOutput->Add(fhistokaonnsigmaTPCMC[kk][jj]);

 Histtname="fhistokaonnsigmaTOFMC";Histtname+=kk;Histtname+=jj;//Histtname+=ii;//Histtname+=kk;
   fhistokaonnsigmaTOFMC[kk][jj]=new TH1F (Histtname.Data(),"fhistokaonnsigmaTOFMC",1200, -60., 60.);
   fOutput->Add(fhistokaonnsigmaTOFMC[kk][jj]);

 Histtname="fhistokaonnsigmaTPCTOFMC";Histtname+=kk;Histtname+=jj;//Histtname+=ii;//Histtname+=kk;
 fhistokaonnsigmaTPCTOFMC[kk][jj]=new TH1F (Histtname.Data(),"fhistokaonnsigmaTPCTOFMC",1200, -60., 60.);
 fOutput->Add(fhistokaonnsigmaTPCTOFMC[kk][jj]);
  
 Histtname="fhistoprotonnsigmaTPCMC";Histtname+=kk;Histtname+=jj;//Histtname+=ii;//Histtname+=kk;
   fhistoprotonnsigmaTPCMC[kk][jj]=new TH1F (Histtname.Data(),"fhistoprotonnsigmaTPCMC",1200, -60., 60.);
   fOutput->Add(fhistoprotonnsigmaTPCMC[kk][jj]);

 Histtname="fhistoprotonnsigmaTOFMC";Histtname+=kk;Histtname+=jj;//Histtname+=ii;//Histtname+=kk;
   fhistoprotonnsigmaTOFMC[kk][jj]=new TH1F (Histtname.Data(),"fhistoprotonnsigmaTOFMC",1200, -60., 60.);
   fOutput->Add(fhistoprotonnsigmaTOFMC[kk][jj]);

 Histtname="fhistoprotonnsigmaTPCTOFMC";Histtname+=kk;Histtname+=jj;//Histtname+=ii;//Histtname+=kk;
 fhistoprotonnsigmaTPCTOFMC[kk][jj]=new TH1F (Histtname.Data(),"fhistoprotonnsigmaTPCTOFMC",1200, -60., 60.);
 fOutput->Add(fhistoprotonnsigmaTPCTOFMC[kk][jj]);
  
 Histtname="fhistoelectronnsigmaTPCMC";Histtname+=kk;Histtname+=jj;//Histtname+=ii;//Histtname+=kk;
   fhistoelectronnsigmaTPCMC[kk][jj]=new TH1F (Histtname.Data(),"fhistoelectronnsigmaTPCMC",1200, -60., 60.);
   fOutput->Add(fhistoelectronnsigmaTPCMC[kk][jj]);

 Histtname="fhistoelectronnsigmaTOFMC";Histtname+=kk;Histtname+=jj;//Histtname+=ii;//Histtname+=kk;
   fhistoelectronnsigmaTOFMC[kk][jj]=new TH1F (Histtname.Data(),"fhistoelectronnsigmaTOFMC",1200, -60., 60.);
   fOutput->Add(fhistoelectronnsigmaTOFMC[kk][jj]);

 Histtname="fhistoelectronnsigmaTPCTOFMC";Histtname+=kk;Histtname+=jj;//Histtname+=ii;//Histtname+=kk;
 fhistoelectronnsigmaTPCTOFMC[kk][jj]=new TH1F (Histtname.Data(),"fhistoelectronnsigmaTPCTOFMC",1200, -60., 60.);
 fOutput->Add(fhistoelectronnsigmaTPCTOFMC[kk][jj]);
   
    }

    }

TString Histoname;
for(Int_t jj=0;jj<4;jj++)
    {
for(Int_t kk=0;kk<10;kk++)
   {
     //for(Int_t ii=0;ii<4;ii++)
     //{
     //trigger particle counting with centrality,zvtx,pttrig binning(use for weighted average)
     Histoname="fEtaSpectraTrig";Histoname+=jj;Histoname+=kk;//Histoname+=ii;
  fEtaSpectraTrig[jj][kk]=new TH1F (Histoname.Data(),"#eta distribution",180,-0.9,0.9);
  fOutput->Add(fEtaSpectraTrig[jj][kk]);

  Histoname="fEtaSpectraTrigbaryon";Histoname+=jj;Histoname+=kk;//Histoname+=ii;
 fEtaSpectraTrigbaryon[jj][kk]=new TH1F (Histoname.Data(),"#eta distribution",180,-0.9,0.9);
  fOutput->Add(fEtaSpectraTrigbaryon[jj][kk]);

  Histoname="fEtaSpectraTrigmeson";Histoname+=jj;Histoname+=kk;//Histoname+=ii;
 fEtaSpectraTrigmeson[jj][kk]=new TH1F (Histoname.Data(),"#eta distribution",180,-0.9,0.9);
  fOutput->Add(fEtaSpectraTrigmeson[jj][kk]);
  //}
  }
 }

TString Histiname;
for(Int_t jj=0;jj<4;jj++)//centrality binning
    {
for(Int_t ii=0;ii<2;ii++)//asso pt binning
   {
for(Int_t kk=0;kk<10;kk++)//zvtx binning
   {
     //for(Int_t pp=0;pp<4;pp++) //trigger pt binning ,do weighted average
     //{
     Histiname="falltrigallasso";Histiname+=jj;Histiname+=ii;Histiname+=kk;//Histiname+=pp;
falltrigallasso[jj][ii][kk]=new TH2F (Histiname.Data(),"#delta#eta-#delta#phi corr",72,-TMath::Pi()/2,3*TMath::Pi()/2,36,-1.8,1.8);
 falltrigallasso[jj][ii][kk]->Sumw2();
fOutput->Add(falltrigallasso[jj][ii][kk]);

Histiname="falltrigpionasso";Histiname+=jj;Histiname+=ii;Histiname+=kk;//Histiname+=pp;
falltrigpionasso[jj][ii][kk]=new TH2F (Histiname.Data(),"#delta#eta-#delta#phi corr",72,-TMath::Pi()/2,3*TMath::Pi()/2,36,-1.8,1.8);
 falltrigpionasso[jj][ii][kk]->Sumw2();
fOutput->Add(falltrigpionasso[jj][ii][kk]);

Histiname="falltrigkaonasso";Histiname+=jj;Histiname+=ii;Histiname+=kk;//Histiname+=pp;
falltrigkaonasso[jj][ii][kk]=new TH2F (Histiname.Data(),"#delta#eta-#delta#phi corr",72,-TMath::Pi()/2,3*TMath::Pi()/2,36,-1.8,1.8);
 falltrigkaonasso[jj][ii][kk]->Sumw2();
fOutput->Add(falltrigkaonasso[jj][ii][kk]);

Histiname="falltrigprotonasso";Histiname+=jj;Histiname+=ii;Histiname+=kk;//Histiname+=pp;
falltrigprotonasso[jj][ii][kk]=new TH2F (Histiname.Data(),"#delta#eta-#delta#phi corr",72,-TMath::Pi()/2,3*TMath::Pi()/2,36,-1.8,1.8);
 falltrigpionasso[jj][ii][kk]->Sumw2();
fOutput->Add(falltrigprotonasso[jj][ii][kk]);

Histiname="fbaryontrigallasso";Histiname+=jj;Histiname+=ii;Histiname+=kk;//Histiname+=pp;
fbaryontrigallasso[jj][ii][kk]=new TH2F (Histiname.Data(),"#delta#eta-#delta#phi corr",72,-TMath::Pi()/2,3*TMath::Pi()/2,36,-1.8,1.8);
 fbaryontrigallasso[jj][ii][kk]->Sumw2();
   fOutput->Add(fbaryontrigallasso[jj][ii][kk]);

 Histiname="fbaryontrigpionasso";Histiname+=jj;Histiname+=ii;Histiname+=kk;//Histiname+=pp;
fbaryontrigpionasso[jj][ii][kk]=new TH2F (Histiname.Data(),"#delta#eta-#delta#phi corr",72,-TMath::Pi()/2,3*TMath::Pi()/2,36,-1.8,1.8);
 fbaryontrigpionasso[jj][ii][kk]->Sumw2();
   fOutput->Add(fbaryontrigpionasso[jj][ii][kk]);

   Histiname="fbaryontrigkaonasso";Histiname+=jj;Histiname+=ii;Histiname+=kk;//Histiname+=pp;
fbaryontrigkaonasso[jj][ii][kk]=new TH2F (Histiname.Data(),"#delta#eta-#delta#phi corr",72,-TMath::Pi()/2,3*TMath::Pi()/2,36,-1.8,1.8);
 fbaryontrigkaonasso[jj][ii][kk]->Sumw2();
fOutput->Add(fbaryontrigkaonasso[jj][ii][kk]);

 Histiname="fbaryontrigprotonasso";Histiname+=jj;Histiname+=ii;Histiname+=kk;//Histiname+=pp;
fbaryontrigprotonasso[jj][ii][kk]=new TH2F (Histiname.Data(),"#delta#eta-#delta#phi corr",72,-TMath::Pi()/2,3*TMath::Pi()/2,36,-1.8,1.8);
 fbaryontrigprotonasso[jj][ii][kk]->Sumw2();
fOutput->Add(fbaryontrigprotonasso[jj][ii][kk]);

Histiname="fmesontrigallasso";Histiname+=jj;Histiname+=ii;Histiname+=kk;//Histiname+=pp;
fmesontrigallasso[jj][ii][kk]=new TH2F (Histiname.Data(),"#delta#eta-#delta#phi corr",72,-TMath::Pi()/2,3*TMath::Pi()/2,36,-1.8,1.8);
 fmesontrigallasso[jj][ii][kk]->Sumw2();
fOutput->Add(fmesontrigallasso[jj][ii][kk]);

 Histiname="fmesontrigpionasso";Histiname+=jj;Histiname+=ii;Histiname+=kk;//Histiname+=pp;
fmesontrigpionasso[jj][ii][kk]=new TH2F (Histiname.Data(),"#delta#eta-#delta#phi corr",72,-TMath::Pi()/2,3*TMath::Pi()/2,36,-1.8,1.8);
 fmesontrigpionasso[jj][ii][kk]->Sumw2();
fOutput->Add(fmesontrigpionasso[jj][ii][kk]);
 
 Histiname="fmesontrigkaonasso";Histiname+=jj;Histiname+=ii;Histiname+=kk;//Histiname+=pp;
fmesontrigkaonasso[jj][ii][kk]=new TH2F (Histiname.Data(),"#delta#eta-#delta#phi corr",72,-TMath::Pi()/2,3*TMath::Pi()/2,36,-1.8,1.8);
 fmesontrigkaonasso[jj][ii][kk]->Sumw2();
fOutput->Add(fmesontrigkaonasso[jj][ii][kk]);

 Histiname="fmesontrigprotonasso";Histiname+=jj;Histiname+=ii;Histiname+=kk;//Histiname+=pp;
fmesontrigprotonasso[jj][ii][kk]=new TH2F (Histiname.Data(),"#delta#eta-#delta#phi corr",72,-TMath::Pi()/2,3*TMath::Pi()/2,36,-1.8,1.8);
 fmesontrigprotonasso[jj][ii][kk]->Sumw2();
fOutput->Add(fmesontrigprotonasso[jj][ii][kk]);

Histiname="falltrigallassomix";Histiname+=jj;Histiname+=ii;Histiname+=kk;//Histiname+=pp;
falltrigallassomix[jj][ii][kk]=new TH2F (Histiname.Data(),"#delta#eta-#delta#phi corr",72,-TMath::Pi()/2,3*TMath::Pi()/2,36,-1.8,1.8);
 falltrigallassomix[jj][ii][kk]->Sumw2();
fOutput->Add(falltrigallassomix[jj][ii][kk]);

Histiname="falltrigpionassomix";Histiname+=jj;Histiname+=ii;Histiname+=kk;//Histiname+=pp;
falltrigpionassomix[jj][ii][kk]=new TH2F (Histiname.Data(),"#delta#eta-#delta#phi corr",72,-TMath::Pi()/2,3*TMath::Pi()/2,36,-1.8,1.8);
 falltrigpionassomix[jj][ii][kk]->Sumw2();
fOutput->Add(falltrigpionassomix[jj][ii][kk]);

Histiname="falltrigkaonassomix";Histiname+=jj;Histiname+=ii;Histiname+=kk;//Histiname+=pp;
falltrigkaonassomix[jj][ii][kk]=new TH2F (Histiname.Data(),"#delta#eta-#delta#phi corr",72,-TMath::Pi()/2,3*TMath::Pi()/2,36,-1.8,1.8);
 falltrigkaonassomix[jj][ii][kk]->Sumw2();
fOutput->Add(falltrigkaonassomix[jj][ii][kk]);

Histiname="falltrigprotonassomix";Histiname+=jj;Histiname+=ii;Histiname+=kk;//Histiname+=pp;
falltrigprotonassomix[jj][ii][kk]=new TH2F (Histiname.Data(),"#delta#eta-#delta#phi corr",72,-TMath::Pi()/2,3*TMath::Pi()/2,36,-1.8,1.8);
 falltrigprotonassomix[jj][ii][kk]->Sumw2();
fOutput->Add(falltrigprotonassomix[jj][ii][kk]);


Histiname="fbaryontrigallassomix";Histiname+=jj;Histiname+=ii;Histiname+=kk;//Histiname+=pp;
fbaryontrigallassomix[jj][ii][kk]=new TH2F (Histiname.Data(),"#delta#eta-#delta#phi corr",72,-TMath::Pi()/2,3*TMath::Pi()/2,36,-1.8,1.8);
 fbaryontrigallassomix[jj][ii][kk]->Sumw2();
   fOutput->Add(fbaryontrigallassomix[jj][ii][kk]);

Histiname="fbaryontrigpionassomix";Histiname+=jj;Histiname+=ii;Histiname+=kk;//Histiname+=pp;
fbaryontrigpionassomix[jj][ii][kk]=new TH2F (Histiname.Data(),"#delta#eta-#delta#phi corr",72,-TMath::Pi()/2,3*TMath::Pi()/2,36,-1.8,1.8);
 fbaryontrigpionassomix[jj][ii][kk]->Sumw2();
   fOutput->Add(fbaryontrigpionassomix[jj][ii][kk]);

   Histiname="fbaryontrigkaonassomix";Histiname+=jj;Histiname+=ii;Histiname+=kk;//Histiname+=pp;
fbaryontrigkaonassomix[jj][ii][kk]=new TH2F (Histiname.Data(),"#delta#eta-#delta#phi corr",72,-TMath::Pi()/2,3*TMath::Pi()/2,36,-1.8,1.8);
 fbaryontrigkaonassomix[jj][ii][kk]->Sumw2();
fOutput->Add(fbaryontrigkaonassomix[jj][ii][kk]);

 Histiname="fbaryontrigprotonassomix";Histiname+=jj;Histiname+=ii;Histiname+=kk;//Histiname+=pp;
fbaryontrigprotonassomix[jj][ii][kk]=new TH2F (Histiname.Data(),"#delta#eta-#delta#phi corr",72,-TMath::Pi()/2,3*TMath::Pi()/2,36,-1.8,1.8);
 fbaryontrigprotonassomix[jj][ii][kk]->Sumw2();
fOutput->Add(fbaryontrigprotonassomix[jj][ii][kk]);

Histiname="fmesontrigallassomix";Histiname+=jj;Histiname+=ii;Histiname+=kk;//Histiname+=pp;
fmesontrigallassomix[jj][ii][kk]=new TH2F (Histiname.Data(),"#delta#eta-#delta#phi corr",72,-TMath::Pi()/2,3*TMath::Pi()/2,36,-1.8,1.8);
 fmesontrigallassomix[jj][ii][kk]->Sumw2();
fOutput->Add(fmesontrigallassomix[jj][ii][kk]);

 Histiname="fmesontrigpionassomix";Histiname+=jj;Histiname+=ii;Histiname+=kk;//Histiname+=pp;
fmesontrigpionassomix[jj][ii][kk]=new TH2F (Histiname.Data(),"#delta#eta-#delta#phi corr",72,-TMath::Pi()/2,3*TMath::Pi()/2,36,-1.8,1.8);
 fmesontrigpionassomix[jj][ii][kk]->Sumw2();
fOutput->Add(fmesontrigpionassomix[jj][ii][kk]);
 
 Histiname="fmesontrigkaonassomix";Histiname+=jj;Histiname+=ii;Histiname+=kk;//Histiname+=pp;
fmesontrigkaonassomix[jj][ii][kk]=new TH2F (Histiname.Data(),"#delta#eta-#delta#phi corr",72,-TMath::Pi()/2,3*TMath::Pi()/2,36,-1.8,1.8);
 fmesontrigkaonassomix[jj][ii][kk]->Sumw2();
fOutput->Add(fmesontrigkaonassomix[jj][ii][kk]);

 Histiname="fmesontrigprotonassomix";Histiname+=jj;Histiname+=ii;Histiname+=kk;//Histiname+=pp;
fmesontrigprotonassomix[jj][ii][kk]=new TH2F (Histiname.Data(),"#delta#eta-#delta#phi corr",72,-TMath::Pi()/2,3*TMath::Pi()/2,36,-1.8,1.8);
 fmesontrigprotonassomix[jj][ii][kk]->Sumw2();
fOutput->Add(fmesontrigprotonassomix[jj][ii][kk]);
//} 
    }
   }
 }

if(fAnalysisType == "MCAOD") {
     const Int_t nDim = 4;//       cent zvtx  pt   eta
     Int_t fBinsCh[nDim] = {20, 10, 90 , 16};//******************************************change it
     Double_t fMinCh[nDim] = { 0.0, -10.0 , 0.0,-0.8 };
     Double_t fMaxCh[nDim] = { 100.0, 10.0, 9.0,0.8};

TString Histrename;
for(Int_t jj=0;jj<6;jj++)//centrality binning
    {
   Histrename="fTHnrecoallPid";Histrename+=jj;
  fTHnrecoallPid[jj] = new THnF(Histrename.Data(),"cent:zvtx::Pt:eta", nDim, fBinsCh, fMinCh, fMaxCh); 
 fTHnrecoallPid[jj]->Sumw2();
  fTHnrecoallPid[jj]->GetAxis(0)->SetTitle("Centrality");
  fTHnrecoallPid[jj]->GetAxis(1)->SetTitle("zvtx");
  fTHnrecoallPid[jj]->GetAxis(2)->SetTitle("Pt");
  fTHnrecoallPid[jj]->GetAxis(3)->SetTitle("eta");
  fOutput->Add(fTHnrecoallPid[jj]);

Histrename="fTHngenprimPidTruth";Histrename+=jj;
  fTHngenprimPidTruth[jj] = new THnF(Histrename.Data(),"cent:zvtx::Pt:eta", nDim, fBinsCh, fMinCh, fMaxCh);
  fTHngenprimPidTruth[jj]->Sumw2(); 
  fTHngenprimPidTruth[jj]->GetAxis(0)->SetTitle("Centrality");
  fTHngenprimPidTruth[jj]->GetAxis(1)->SetTitle("zvtx");
  fTHngenprimPidTruth[jj]->GetAxis(2)->SetTitle("Pt");
  fTHngenprimPidTruth[jj]->GetAxis(3)->SetTitle("eta");
  fOutput->Add(fTHngenprimPidTruth[jj]);
    }
 }

     const Int_t nDim = 4;//       cent zvtx  pt   eta
     Int_t fBinsCh[nDim] = {20, 10, 90 , 16};//******************************************change it
     Double_t fMinCh[nDim] = { 0.0, -10.0 , 0.0,-0.8 };
     Double_t fMaxCh[nDim] = { 100.0, 10.0, 9.0,0.8};

  TString Histrexname;
for(Int_t jj=0;jj<6;jj++)//centrality binning
    {
  Histrexname="effcorection";Histrexname+=jj;
  effcorection[jj] = new THnF(Histrexname.Data(),"cent:zvtx::Pt:eta", nDim, fBinsCh, fMinCh, fMaxCh);
  effcorection[jj]->Sumw2(); 
  effcorection[jj]->GetAxis(0)->SetTitle("Centrality");
  effcorection[jj]->GetAxis(1)->SetTitle("zvtx");
  effcorection[jj]->GetAxis(2)->SetTitle("Pt");
  effcorection[jj]->GetAxis(3)->SetTitle("eta");
  fOutput->Add(effcorection[jj]);
  /*
Histrexname="effmap";Histrexname+=jj;
  effmap[jj] = new THnF(Histrexname.Data(),"cent:zvtx::Pt:eta", nDim, fBinsCh, fMinCh, fMaxCh);
  effmap[jj]->Sumw2(); 
  effmap[jj]->GetAxis(0)->SetTitle("Centrality");
  effmap[jj]->GetAxis(1)->SetTitle("zvtx");
  effmap[jj]->GetAxis(2)->SetTitle("Pt");
  effmap[jj]->GetAxis(3)->SetTitle("eta");
  //fOutput->Add(effcorection[jj]);
  */
    }


TString Histmcname;
for(Int_t jj=0;jj<4;jj++)
  {//centrality binning 

Histmcname="recoallpt";Histmcname+=jj;
  recoallpt[jj]=new TH1F (Histmcname.Data(),"ptdistributionrecoall",900,0.,9.);
  fOutput->Add(recoallpt[jj]);
  
Histmcname="MCrecomatchedprimpt";Histmcname+=jj;
  MCrecomatchedprimpt[jj]=new TH1F (Histmcname.Data(),"ptdistributionrecomatchedprim",900,0.,9.);
  fOutput->Add(MCrecomatchedprimpt[jj]);

Histmcname="recoalleta";Histmcname+=jj;
  recoalleta[jj]=new TH1F (Histmcname.Data(),"etadistributionrecoall",360,-1.8,1.8);
  fOutput->Add(recoalleta[jj]);

Histmcname="alltrigeta";Histmcname+=jj;
  alltrigeta[jj]=new TH1F (Histmcname.Data(),"alltrigeta",360,-1.8,1.8);
  fOutput->Add(alltrigeta[jj]);

Histmcname="allassoeta";Histmcname+=jj;
  allassoeta[jj]=new TH1F (Histmcname.Data(),"allassoeta",360,-1.8,1.8);
  fOutput->Add(allassoeta[jj]);

Histmcname="baryontrigeta";Histmcname+=jj;
  baryontrigeta[jj]=new TH1F (Histmcname.Data(),"baryontrigeta",360,-1.8,1.8);
  fOutput->Add(baryontrigeta[jj]);


Histmcname="mesontrigeta";Histmcname+=jj;
  mesontrigeta[jj]=new TH1F (Histmcname.Data(),"mesontrigeta",360,-1.8,1.8);
  fOutput->Add(mesontrigeta[jj]);


Histmcname="pionassoeta";Histmcname+=jj;
  pionassoeta[jj]=new TH1F (Histmcname.Data(),"pionassoeta",360,-1.8,1.8);
  fOutput->Add(pionassoeta[jj]);


Histmcname="kaonassoeta";Histmcname+=jj;
  kaonassoeta[jj]=new TH1F (Histmcname.Data(),"kaonassoeta",360,-1.8,1.8);
  fOutput->Add(kaonassoeta[jj]);

Histmcname="protonassoeta";Histmcname+=jj;
  protonassoeta[jj]=new TH1F (Histmcname.Data(),"protonassoeta",360,-1.8,1.8);
  fOutput->Add(protonassoeta[jj]);



Histmcname="MCrecomatchedprimeta";Histmcname+=jj;
  MCrecomatchedprimeta[jj]=new TH1F (Histmcname.Data(),"etadistributionrecomatchedprim",360,-1.8,1.8);
  fOutput->Add(MCrecomatchedprimeta[jj]);

Histmcname="recoallphi";Histmcname+=jj;
  recoallphi[jj]=new TH1F (Histmcname.Data(),"phidistrecoall",340,0,6.8);
  fOutput->Add(recoallphi[jj]);
  Histmcname="MCrecomatchedprimphi";Histmcname+=jj;
  MCrecomatchedprimphi[jj]=new TH1F (Histmcname.Data(),"phidistrecomatchedprim",340,0,6.8);
  fOutput->Add(MCrecomatchedprimphi[jj]);

Histmcname="MCtruthpt";Histmcname+=jj;
  MCtruthpt[jj]=new TH1F (Histmcname.Data(),"ptdistributiontruthprim",900,0.,9.);
  fOutput->Add(MCtruthpt[jj]);

Histmcname="MCtrutheta";Histmcname+=jj;
  MCtrutheta[jj]=new TH1F (Histmcname.Data(),"etadistributiontruthprim",360,-1.8,1.8);
  fOutput->Add(MCtrutheta[jj]);

Histmcname="MCtruthphi";Histmcname+=jj;
  MCtruthphi[jj]=new TH1F (Histmcname.Data(),"phidisttruthprim",340,0,6.8);
  fOutput->Add(MCtruthphi[jj]);


Histmcname="MCtruthpionpt";Histmcname+=jj;
  MCtruthpionpt[jj]=new TH1F (Histmcname.Data(),"MCtruthpionpt",900,0.,9.);
  fOutput->Add(MCtruthpionpt[jj]);

Histmcname="MCtruthpioneta";Histmcname+=jj;
  MCtruthpioneta[jj]=new TH1F (Histmcname.Data(),"MCtruthpioneta",360,-1.8,1.8);
  fOutput->Add(MCtruthpioneta[jj]);


Histmcname="MCtruthpionphi";Histmcname+=jj;
  MCtruthpionphi[jj]=new TH1F (Histmcname.Data(),"MCtruthpionphi",340,0,6.8);
  fOutput->Add(MCtruthpionphi[jj]);


Histmcname="MCtruthkaonpt";Histmcname+=jj;
  MCtruthkaonpt[jj]=new TH1F (Histmcname.Data(),"MCtruthkaonpt",900,0.,9.);
  fOutput->Add(MCtruthkaonpt[jj]);

Histmcname="MCtruthkaoneta";Histmcname+=jj;
  MCtruthkaoneta[jj]=new TH1F (Histmcname.Data(),"MCtruthkaoneta",360,-1.8,1.8);
  fOutput->Add(MCtruthkaoneta[jj]);

Histmcname="MCtruthkaonphi";Histmcname+=jj;
  MCtruthkaonphi[jj]=new TH1F (Histmcname.Data(),"MCtruthkaonphi",340,0,6.8);
  fOutput->Add(MCtruthkaonphi[jj]);


Histmcname="MCtruthprotonpt";Histmcname+=jj;
  MCtruthprotonpt[jj]=new TH1F (Histmcname.Data(),"MCtruthprotonpt",900,0.,9.);
  fOutput->Add(MCtruthprotonpt[jj]);

Histmcname="MCtruthprotoneta";Histmcname+=jj;
  MCtruthprotoneta[jj]=new TH1F (Histmcname.Data(),"MCtruthprotoneta",360,-1.8,1.8);
  fOutput->Add(MCtruthprotoneta[jj]);

Histmcname="MCtruthprotonphi";Histmcname+=jj;
  MCtruthprotonphi[jj]=new TH1F (Histmcname.Data(),"MCtruthprotonphi",340,0,6.8);
  fOutput->Add(MCtruthprotonphi[jj]);

Histmcname="MCrecomatchedprimpionpt";Histmcname+=jj;
  MCrecomatchedprimpionpt[jj]=new TH1F (Histmcname.Data(),"MCrecomatchedprimpionpt",900,0.,9.);
  fOutput->Add(MCrecomatchedprimpionpt[jj]);

Histmcname="MCrecomatchedprimpioneta";Histmcname+=jj;
  MCrecomatchedprimpioneta[jj]=new TH1F (Histmcname.Data(),"MCrecomatchedprimpioneta",360,-1.8,1.8);
  fOutput->Add(MCrecomatchedprimpioneta[jj]);

Histmcname="MCrecomatchedprimpionphi";Histmcname+=jj;
  MCrecomatchedprimpionphi[jj]=new TH1F (Histmcname.Data(),"MCrecomatchedprimpionphi",340,0,6.8);
  fOutput->Add(MCrecomatchedprimpionphi[jj]);

Histmcname="MCrecomatchedprimkaonpt";Histmcname+=jj;
  MCrecomatchedprimkaonpt[jj]=new TH1F (Histmcname.Data(),"MCrecomatchedprimkaonpt",900,0.,9.);
  fOutput->Add(MCrecomatchedprimkaonpt[jj]);

Histmcname="MCrecomatchedprimkaoneta";Histmcname+=jj;
  MCrecomatchedprimkaoneta[jj]=new TH1F (Histmcname.Data(),"MCrecomatchedprimkaoneta",360,-1.8,1.8);
  fOutput->Add(MCrecomatchedprimkaoneta[jj]);

Histmcname="MCrecomatchedprimkaonphi";Histmcname+=jj;
  MCrecomatchedprimkaonphi[jj]=new TH1F (Histmcname.Data(),"MCrecomatchedprimkaonphi",340,0,6.8);
  fOutput->Add(MCrecomatchedprimkaonphi[jj]);


Histmcname="MCrecomatchedprimprotonpt";Histmcname+=jj;
  MCrecomatchedprimprotonpt[jj]=new TH1F (Histmcname.Data(),"MCrecomatchedprimprotonpt",900,0.,9.);
  fOutput->Add(MCrecomatchedprimprotonpt[jj]);

Histmcname="MCrecomatchedprimprotoneta";Histmcname+=jj;
  MCrecomatchedprimprotoneta[jj]=new TH1F (Histmcname.Data(),"MCrecomatchedprimprotoneta",360,-1.8,1.8);
  fOutput->Add(MCrecomatchedprimprotoneta[jj]);

Histmcname="MCrecomatchedprimprotonphi";Histmcname+=jj;
  MCrecomatchedprimprotonphi[jj]=new TH1F (Histmcname.Data(),"MCrecomatchedprimprotonphi",340,0,6.8);
  fOutput->Add(MCrecomatchedprimprotonphi[jj]);


 }

//Mixing
DefineEventPool();

  if(applyefficiency)
   {
     TFile *fsifile = new TFile("map32.root","READ");
 TString Nameg;
for(Int_t jj=0;jj<6;jj++)//type binning
    {
Nameg="effmap";Nameg+=jj;
effcorection[jj] = (THnF*)fsifile->Get(Nameg.Data());
//effcorection[jj]->SetDirectory(0);//****************************not present in case oh THnF
    }
fsifile->Close();
   }
    
//fControlConvResoncances = new TH2F("fControlConvResoncances", ";id;delta mass", 3, -0.5, 2.5, 100, -0.1, 0.1);
// fOutput->Add(fControlConvResoncances);

 
  PostData(1, fOutput);              // Post data for ALL output slots >0 here, to get at least an empty histogram
}
//-------------------------------------------------------------------------------
void AliTwoParticlePIDCorr::UserExec( Option_t * ){

 
  if(fAnalysisType == "AOD") {

    doAODevent();
    
  }//AOD--analysis-----

  else if(fAnalysisType == "MCAOD") {
  
    doMCAODevent();
    
  }
  
  else return;
  
}
//-------------------------------------------------------------------------
void AliTwoParticlePIDCorr::doMCAODevent() 
{
  AliVEvent *event = InputEvent();
  if (!event) { Printf("ERROR: Could not retrieve event"); return; }
  AliAODEvent* aod = dynamic_cast<AliAODEvent*>(event);
  if (!aod) {
    AliError("Cannot get the AOD event");
    return;
  }
 
// count all events   
  fEventCounter->Fill(1);

 
// get centrality object and check quality
  Float_t cent_v0m=-999;
  AliCentrality *centrality=0;
  if(aod) 
 centrality = aod->GetHeader()->GetCentralityP();
  
if(centrality)
  { 
    // if (centrality->GetQuality() != 0) return ;
  cent_v0m = centrality->GetCentralityPercentile("V0A");
  //AliInfo(Form("Centrality is %f", cent_v0m));
  }
 else
    {
      Printf("WARNING: Centrality object is 0");
      cent_v0m = -1;
     }

if (cent_v0m < 0)   return;

//do centrality binning(for 0-100% centrality events,reject events outside this range)
 Int_t centbin=Getcentbin(cent_v0m);
 if(centbin==-999) return;

 //check the PIDResponse handler
  if (!fPID) return;

// get mag. field required for twotrack efficiency cut
 if(!aod) return; //for safety
 Float_t bSign = 0;
 bSign = (aod->GetMagneticField() > 0) ? 1 : -1;

 //check for TClonesArray(truth track MC information)
fArrayMC = dynamic_cast<TClonesArray*>(aod->FindListObject(AliAODMCParticle::StdBranchName()));
  if (!fArrayMC) {
    AliFatal("Error: MC particles branch not found!\n");
    return;
  }
  
  //check for AliAODMCHeader(truth event MC information)
  AliAODMCHeader *header=NULL;
  header=(AliAODMCHeader*)aod->GetList()->FindObject(AliAODMCHeader::StdBranchName());  
  if(!header) {
    printf("MC header branch not found!\n");
    return;
  }
 


 //count events having centrality betn 0-100%
  fEventCounter->Fill(2);


//Only consider MC events within the vtx-z region used also as cut on the reconstructed vertex
Float_t zVtxmc =header->GetVtxZ();
 if(TMath::Abs(zVtxmc)>10.) return;

 // For productions with injected signals, figure out above which label to skip particles/tracks
  Int_t skipParticlesAbove = 0;

 if (fInjectedSignals)
  {
    AliGenEventHeader* eventHeader = 0;
    Int_t headers = 0;

// AOD
      if (!header)
      AliFatal("fInjectedSignals set but no MC header found");
      
      headers = header->GetNCocktailHeaders();
      eventHeader = header->GetCocktailHeader(0);

 if (!eventHeader)
    {
      // We avoid AliFatal here, because the AOD productions sometimes have events where the MC header is missing 
      // (due to unreadable Kinematics) and we don't want to loose the whole job because of a few events
      AliError("First event header not found. Skipping this event.");
      //fHistos->FillEvent(centrality, AliUEHist::kCFStepAnaTopology);
      return;
    }
skipParticlesAbove = eventHeader->NProduced();
    AliInfo(Form("Injected signals in this event (%d headers). Keeping events of %s. Will skip particles/tracks above %d.", headers, eventHeader->ClassName(), skipParticlesAbove));
  }

 //determine the two particle correlation function with generator level (cleaned up) primary particles

  // Trigger selection ************************************************



// Vertex selection *************************************************
   AliAODVertex* trkVtx = aod->GetPrimaryVertex();
  if (!trkVtx || trkVtx->GetNContributors()<=0) return;
  TString vtxTtl = trkVtx->GetTitle();
  if (!vtxTtl.Contains("VertexerTracks")) return;
  Float_t zvtx = trkVtx->GetZ();
  const AliAODVertex* spdVtx = aod->GetPrimaryVertexSPD();
  if (spdVtx->GetNContributors()<=0) return;
  TString vtxTyp = spdVtx->GetTitle();
  Double_t cov[6]={0};
  spdVtx->GetCovarianceMatrix(cov);
  Double_t zRes = TMath::Sqrt(cov[5]);
  if (vtxTyp.Contains("vertexer:Z") && (zRes>0.25)) return;
  if (TMath::Abs(spdVtx->GetZ() - trkVtx->GetZ())>0.5) return;

  fHistQA[0]->Fill((trkVtx->GetX()));fHistQA[1]->Fill((trkVtx->GetY()));fHistQA[2]->Fill((trkVtx->GetZ()));   //for trkVtx only before vertex cut |zvtx|<10 cm

  //count events having a proper vertex
   fEventCounter->Fill(5);

 if (TMath::Abs(zvtx) > 10) return;

fHistQA[3]->Fill((trkVtx->GetX()));fHistQA[4]->Fill((trkVtx->GetY()));fHistQA[5]->Fill((trkVtx->GetZ()));//after vertex cut for trkVtx only

//do zvtx binning
  Int_t vtx=Getzbin(zvtx);
  if(vtx==-999) return;//all the events outside the defined vtx range will be ignored

  //now we have events passed physics trigger, centrality,zvtx cut 

  //count events after vertex cut
  fEventCounter->Fill(6);
 
  //centrality dist. of accepted events
 fhistcentrality->Fill(cent_v0m);

  eventno++;
    
if(!aod) return;  //for safety

   TObjArray* trackstrig = new TObjArray;
   TObjArray* tracksasso = new TObjArray;
   trackstrig->SetOwner(kTRUE);  // IMPORTANT!
   tracksasso->SetOwner(kTRUE);  // IMPORTANT!

   Float_t bSign1=aod->GetHeader()->GetMagneticField() ;//used for reconstructed track dca cut
   Int_t nooftracks=0;

// loop over reconstructed tracks 
  for (Int_t itrk = 0; itrk < aod->GetNumberOfTracks(); itrk++) 
{ // reconstructed track loop starts
  AliAODTrack* track = dynamic_cast<AliAODTrack*>(aod->GetTrack(itrk));
  if (!track) continue;
 //get the corresponding MC track at the truth level 
  AliAODMCParticle* recomatched = static_cast<AliAODMCParticle*>(fArrayMC->At(TMath::Abs(track->GetLabel())));
  if(!recomatched) continue;//if a reco track doesn't have corresponding truth track at generated level is a fake track, ignore it

//remove injected signals 
 if(fInjectedSignals)
   {
    AliAODMCParticle* mother = recomatched;

      while (!mother->IsPhysicalPrimary())
      {// find the primary mother;the first stable mother is searched and checked if it is <= <maxLabel>
	if (mother->GetMother() < 0)
	{
	  mother = 0;
	  break;
	}
	  
   mother =(AliAODMCParticle*) fArrayMC->At(((AliAODMCParticle*)mother)->GetMother());
	if (!mother)
	  break;
      }
 if (!mother)
    {
      Printf("WARNING: No mother found for particle %d:", recomatched->GetLabel());
      continue;
    }
 if (mother->GetLabel() >= skipParticlesAbove) continue;//remove injected signals(primaries above <maxLabel>)
   }//remove injected signals

 if (fRemoveWeakDecays && ((AliAODMCParticle*) recomatched)->IsSecondaryFromWeakDecay()) continue;//remove weak decays
	
  Bool_t isduplicate2=kFALSE;
if (fRemoveDuplicates)
   {
  for (Int_t j =itrk+1; j < aod->GetNumberOfTracks(); j++) 
    {//2nd loop starts
 AliAODTrack* track2 = dynamic_cast<AliAODTrack*>(aod->GetTrack(j));
 if (!track2) continue;
 AliAODMCParticle* recomatched2 = static_cast<AliAODMCParticle*>(fArrayMC->At(TMath::Abs(track2->GetLabel())));
if(!recomatched2) continue;

if (track->GetLabel() == track2->GetLabel())
   {
isduplicate2=kTRUE;
 break;  
   }
    }//2nd loop ends
   }
 if(fRemoveDuplicates && isduplicate2) continue;//remove duplicates
     
  fHistQA[11]->Fill(track->GetTPCNcls());
   Int_t tracktype=ClassifyTrack(track,centbin,trkVtx,bSign1,kFALSE);//dcacut=kFALSE

 if(tracktype==0) continue; 
 if(tracktype==1)//tracks "not" passed AliAODTrack::kPrimary at reconstructed level & have proper TPC PID response(?)
{
 //accepted all(primaries+secondary) reconstructed tracks(pt 0.2 to 9.0,,eta -0.8 to 0.8)
  nooftracks++;

 MCrecomatchedprimpt[centbin]->Fill(recomatched->Pt());
 MCrecomatchedprimeta[centbin]->Fill(recomatched->Eta());
 MCrecomatchedprimphi[centbin]->Fill(recomatched->Phi());


Float_t dEdx = track->GetTPCsignal();
 fHistoTPCdEdx->Fill(track->Pt(), dEdx);

 if(HasTOFPID(track))
{
Float_t beta = GetBeta(track);
fHistoTOFbeta->Fill(track->Pt(), beta);
 }

 Float_t effmatrix=1.;

 //get the pdg code of the corresponding truth particle
 Int_t pdgCode = ((AliAODMCParticle*)recomatched)->GetPdgCode();
 
if (TMath::Abs(pdgCode)==211)
  {
 MCrecomatchedprimpionpt[centbin]->Fill(recomatched->Pt());
 MCrecomatchedprimpioneta[centbin]->Fill(recomatched->Eta());
 MCrecomatchedprimpionphi[centbin]->Fill(recomatched->Phi());
  }

if(TMath::Abs(pdgCode)==321)
  {
 MCrecomatchedprimkaonpt[centbin]->Fill(recomatched->Pt());
 MCrecomatchedprimkaoneta[centbin]->Fill(recomatched->Eta());
 MCrecomatchedprimkaonphi[centbin]->Fill(recomatched->Phi());
  }

if(TMath::Abs(pdgCode)==2212)
  {
 MCrecomatchedprimprotonpt[centbin]->Fill(recomatched->Pt());
 MCrecomatchedprimprotoneta[centbin]->Fill(recomatched->Eta());
 MCrecomatchedprimprotonphi[centbin]->Fill(recomatched->Phi());
  }
  


//now we have only those reconstructed particles each of them have a corresponding  particle(primary+secondary) at the truth level


// -- Fill THnSparse 
 Double_t allrecomatchedpid[4] = {cent_v0m, zVtxmc,recomatched->Pt(), recomatched->Eta()};
 fTHnrecoallPid[5]->Fill(allrecomatchedpid);//for all

 //do track identification(nsigma method)
 Int_t particletypeMC=GetParticle(track,centbin);//******************************problem is here

 //fill tracking efficiency
 if(particletypeMC==SpPion || particletypeMC==SpKaon)
   {
if(TMath::Abs(pdgCode)==211 ||  TMath::Abs(pdgCode)==321)
  {   
  fTHnrecoallPid[4]->Fill(allrecomatchedpid);//for mesons
  if(TMath::Abs(pdgCode)==211)  fTHnrecoallPid[0]->Fill(allrecomatchedpid);//for pions
  if(TMath::Abs(pdgCode)==321)  fTHnrecoallPid[1]->Fill(allrecomatchedpid);//for kaons
  }
 }

 if(particletypeMC==SpProton && TMath::Abs(pdgCode)==2212 )
   fTHnrecoallPid[2]->Fill(allrecomatchedpid);//for protons
  
 if(particletypeMC==SpUndefined && TMath::Abs(pdgCode)!=211 && TMath::Abs(pdgCode)!=321 && TMath::Abs(pdgCode)!=2212) fTHnrecoallPid[3]->Fill(allrecomatchedpid);//for others
 
  
Int_t ptmc1=Getptbin(track->Pt());//trig--0,1,2,3; asso--4,5,6,7,8
if(ptmc1==-999) continue;//remove particles with pt<1.0 Gev/c, now only particleswithin 1.0<=pt<=4.0Gev/c are present

//Fill the nsigma histograms
   Float_t nsigmaTPCPionmc    =  fnsigmas[SpPion][NSigmaTPC];
   Float_t nsigmaTOFPionmc    =  fnsigmas[SpPion][NSigmaTOF];
   Float_t nsigmaTPCTOFPionmc    = TMath::Abs(fnsigmas[SpPion][NSigmaTPCTOF]);//it may be negative if TOF pid is not available(i.e. this is basically nsigmaTPCPionmc)

if (TMath::Abs(pdgCode)==211)//pions
   {
     fhistopionnsigmaTPCMC[centbin][ptmc1]->Fill(nsigmaTPCPionmc); 
     if(nsigmaTOFPionmc!=999) fhistopionnsigmaTOFMC[centbin][ptmc1]->Fill(nsigmaTOFPionmc);
     fhistopionnsigmaTPCTOFMC[centbin][ptmc1]->Fill(nsigmaTPCTOFPionmc); 
   }
 if(TMath::Abs(pdgCode)==321)//kaons 
    {
     fhistokaonnsigmaTPCMC[centbin][ptmc1]->Fill(nsigmaTPCPionmc); 
     if(nsigmaTOFPionmc!=999) fhistokaonnsigmaTOFMC[centbin][ptmc1]->Fill(nsigmaTOFPionmc);
     fhistokaonnsigmaTPCTOFMC[centbin][ptmc1]->Fill(nsigmaTPCTOFPionmc); 
    }
if(TMath::Abs(pdgCode)==2212)//protons 
    {
     fhistoprotonnsigmaTPCMC[centbin][ptmc1]->Fill(nsigmaTPCPionmc); 
     if(nsigmaTOFPionmc!=999) fhistoprotonnsigmaTOFMC[centbin][ptmc1]->Fill(nsigmaTOFPionmc);
     fhistoprotonnsigmaTPCTOFMC[centbin][ptmc1]->Fill(nsigmaTPCTOFPionmc); 
    }
if(TMath::Abs(pdgCode)==11)//electrons 
    {
     fhistoelectronnsigmaTPCMC[centbin][ptmc1]->Fill(nsigmaTPCPionmc); 
     if(nsigmaTOFPionmc!=999) fhistoelectronnsigmaTOFMC[centbin][ptmc1]->Fill(nsigmaTOFPionmc);
     fhistoelectronnsigmaTPCTOFMC[centbin][ptmc1]->Fill(nsigmaTPCTOFPionmc); 
    }


//2-d TPCTOF map(for each Pt interval)
 if(HasTOFPID(track)) fTPCTOFpion2d[ptmc1]->Fill(nsigmaTPCPionmc,nsigmaTOFPionmc);
 
//now depending on the pt of the current particle it will either go into the asso or trigger particle loop

if ((track->Pt()>=1.0) && (track->Pt()<2.0)) 
  {//asso filling loop starts

Int_t ptmc2;
if(ptmc1==4) ptmc2=0;
if(ptmc1==5) ptmc2=1;
if(ptmc1==6) ptmc2=2;
if(ptmc1==7) ptmc2=3;
if(ptmc1==8) ptmc2=4; 

//for purity check
if(particletypeMC==SpPion)//should be pions
   {
  if (TMath::Abs(pdgCode)==211) fhistoassopioncont[ptmc2]->Fill(1);//pions
   
else if(TMath::Abs(pdgCode)==321) fhistoassopioncont[ptmc2]->Fill(3);//kaons 
    
else if(TMath::Abs(pdgCode)==2212)  fhistoassopioncont[ptmc2]->Fill(5);//protons 
   
else if(TMath::Abs(pdgCode)==11) fhistoassopioncont[ptmc2]->Fill(7);//electrons 
  
else  fhistoassopioncont[ptmc2]->Fill(9);//anything else(contamination)

   }//if(particletypeMC==10) condition ends


if(particletypeMC==SpKaon)//should be kaons
   {
if (TMath::Abs(pdgCode)==321) fhistoassokaoncont[ptmc2]->Fill(1);//kaons
   
else if(TMath::Abs(pdgCode)==211) fhistoassokaoncont[ptmc2]->Fill(3);//pions 
   
else if(TMath::Abs(pdgCode)==2212) fhistoassokaoncont[ptmc2]->Fill(5);//protons 
    
else if(TMath::Abs(pdgCode)==11) fhistoassokaoncont[ptmc2]->Fill(7);//electrons 
    
else  fhistoassokaoncont[ptmc2]->Fill(9);//anything else(contamination)

   }//if(particletypeMC==20) condition ends

if(particletypeMC==SpProton)//should be protons
   {
if(TMath::Abs(pdgCode)==2212) fhistoassoprotoncont[ptmc2]->Fill(1);//protons
   
else if(TMath::Abs(pdgCode)==321) fhistoassoprotoncont[ptmc2]->Fill(3);//kaons 
   
else if(TMath::Abs(pdgCode)==211) fhistoassoprotoncont[ptmc2]->Fill(5);//pions 
    
else if(TMath::Abs(pdgCode)==11)  fhistoassoprotoncont[ptmc2]->Fill(7);//electrons 
    
else  fhistoassoprotoncont[ptmc2]->Fill(9);//anything else(contamination)

   }//if(particletypeMC==30) condition ends
fEtaSpectrasso->Fill(track->Eta(),track->Pt());
fphiSpectraasso->Fill(track->Phi(),track->Pt());
if (applyefficiency)
  effmatrix=GetTrackbyTrackeffvalue(track,cent_v0m,zvtx,particletypeMC);
 LRCParticlePID* copy = new LRCParticlePID(particletypeMC,track->Charge(),track->Pt(),track->Eta(), track->Phi(),centbin,vtx,effmatrix);
 copy->SetUniqueID(track->GetUniqueID());
 tracksasso->Add(copy);//fill it with either asso pions,kaons or protons in case of identified associated particles,now all particles are present even e,muon etc

if(particletypeMC==SpPion) pionassoeta[centbin]->Fill(track->Eta());
if(particletypeMC==SpKaon) kaonassoeta[centbin]->Fill(track->Eta());
if(particletypeMC==SpProton) protonassoeta[centbin]->Fill(track->Eta());

  }//asso filling loop ends

//now deal with trigger particles only
 
if ((track->Pt()>=2.0) && (track->Pt()<=4.0))  
  {//trigger filling  loop starts
 if(particletypeMC==SpProton)//only for 2 to 4Gev,should be protons
   {
if (TMath::Abs(pdgCode)==2212)  fhistotrigbaryoncont[ptmc1]->Fill(1);//protons

else if(TMath::Abs(pdgCode)==321) fhistotrigbaryoncont[ptmc1]->Fill(3);//kaons 
    
else if(TMath::Abs(pdgCode)==211)  fhistotrigbaryoncont[ptmc1]->Fill(5);//pions 
   
else if(TMath::Abs(pdgCode)==11) fhistotrigbaryoncont[ptmc1]->Fill(7);//electrons 
    
else  fhistotrigbaryoncont[ptmc1]->Fill(9);//anything else(contamination)

   }//if(particletypeMC==1) condition ends
 
 if(particletypeMC==SpPion || particletypeMC==SpKaon) //only for 2 to 4 Gev,should be either kaons or pions 
  {
if(TMath::Abs(pdgCode)==211) fhistotrigmesoncont[ptmc1]->Fill(1);//pions 
   
else if (TMath::Abs(pdgCode)==321) fhistotrigmesoncont[ptmc1]->Fill(1);//kaons
    
else if(TMath::Abs(pdgCode)==2212) fhistotrigmesoncont[ptmc1]->Fill(3);//protons
    
else if(TMath::Abs(pdgCode)==11) fhistotrigmesoncont[ptmc1]->Fill(5);//electrons 

else fhistotrigmesoncont[ptmc1]->Fill(7);//anything else(contamination)
  }//if(particletypeMC==2) condition ends

//fill up the trigger particle container

     fEtaSpectraTrigall->Fill(track->Eta());   //to know the number of trigger particls
if (applyefficiency)
  effmatrix=GetTrackbyTrackeffvalue(track,cent_v0m,zvtx,particletypeMC);
 LRCParticlePID* copy1 = new LRCParticlePID(particletypeMC,track->Charge(),track->Pt(),track->Eta(), track->Phi(),centbin,vtx,effmatrix);
    copy1->SetUniqueID(track->GetUniqueID());
    trackstrig->Add(copy1);
if(particletypeMC==SpProton) baryontrigeta[centbin]->Fill(track->Eta());
if(particletypeMC==SpPion || particletypeMC==SpKaon) mesontrigeta[centbin]->Fill(track->Eta());
     }//trigger filling loop ends

  }// if(tracktype==1) condition structure ands

}//reco track loop ends

//still in main event loop
 fCentralityCorrelation->Fill(cent_v0m, nooftracks);

Bool_t isbaryontrig=kFALSE;
Bool_t ismesontrig=kFALSE;

if(trackstrig && tracksasso && trackstrig->GetEntriesFast()>0)
  {//same event calculation starts
//calculate no of events in each centrality for different zvtx bins
    fEventno[centbin][vtx]->Fill(5);//only those events which have at least one trigger particle
    Fillcorrelation(trackstrig,tracksasso,centbin,vtx,bSign,kTRUE,kFALSE);//mixcase=kFALSE for same event case

for(Int_t i=0;i<trackstrig->GetEntriesFast();i++)
    {  //trigger loop starts
      LRCParticlePID *trig=(LRCParticlePID*)(trackstrig->UncheckedAt(i));
      if(!trig) continue;
      Int_t particlepidtrig=trig->getparticle(); //either 1 or 2
      if(particlepidtrig==SpProton) isbaryontrig=kTRUE;
      if(particlepidtrig==SpPion || particlepidtrig==SpKaon) ismesontrig=kTRUE;
    }//trig loop ends
 if (isbaryontrig) fEventnobaryon[centbin][vtx]->Fill(5); 
 if (ismesontrig) fEventnomeson[centbin][vtx]->Fill(5);
  }//same event calculation ends

//start mixing
AliEventPool* pool = fPoolMgr->GetEventPool(cent_v0m, zvtx);
//if (!pool)
//AliFatal(Form("No pool found for centrality = %f, zVtx = %f", cent_v0m, zvtx));
if (pool && pool->IsReady())
  {//start mixing only when pool->IsReady
if(trackstrig && trackstrig->GetEntriesFast()>0)
  {//proceed only when no. of trigger particles >0 in current event
    //TObjArray* bgTracks =new TObjArray();
for (Int_t jMix=0; jMix<pool->GetCurrentNEvents(); jMix++) 
  { //pool event loop start
 TObjArray* bgTracks = pool->GetEvent(jMix);
  if(!bgTracks) continue;
  Fillcorrelation(trackstrig,bgTracks,centbin,vtx,bSign,kTRUE,kTRUE);//mixcase=kTRUE for mixing case
  
   }// pool event loop ends mixing case
 }//if(trackstrig && trackstrig->GetEntriesFast()>0) condition ends mixing case
} //if pool->IsReady() condition ends mixing case


 //still in main event loop
if(pool)
 {
if(tracksasso)
pool->UpdatePool(tracksasso);//ownership of tracksasso is with pool now, don't delete it
 }
if(trackstrig)
delete trackstrig;

 //still in main event loop

//now process the truth particles

Int_t nMCTrack = fArrayMC->GetEntriesFast();
  
for (Int_t iMC = 0; iMC < nMCTrack; iMC++) 
{      //MC truth track loop starts
    
AliAODMCParticle *partMC = (AliAODMCParticle*) fArrayMC->At(iMC);
    
    if(!partMC){
      AliError(Form("ERROR: Could not retrieve AODMCtrack %d",iMC));
      continue;
    }

//consider only charged particles
    if(partMC->Charge() == 0) continue;

//consider only primary particles; neglect all secondary particles including from weak decays
    if(!partMC->IsPhysicalPrimary()) continue;


//remove injected signals(primaries above <maxLabel>)
 if (fInjectedSignals && partMC->GetLabel() >= skipParticlesAbove) continue;


  Bool_t isduplicate=kFALSE;
 if (fRemoveDuplicates)
   { 
 for (Int_t j=iMC+1; j<nMCTrack; ++j) 
   {//2nd trutuh loop starts
AliAODMCParticle *partMC2 = (AliAODMCParticle*) fArrayMC->At(j);
   if(!partMC2){
      AliError(Form("ERROR: Could not retrieve AODMCtrack %d",j));
      continue;
    }    
 if (partMC->GetLabel() == partMC2->GetLabel())
   {
isduplicate=kTRUE;
 break;  
   }    
   }//2nd truth loop ends
   }
 if(fRemoveDuplicates && isduplicate) continue;//remove duplicates

 //kinematic cuts    
 if (partMC->Eta() < -0.8 || partMC->Eta() > 0.8) continue;
 if (partMC->Pt() < 0.2 ||  partMC->Pt() > 9.0) continue;

//only physical primary(all)  
 MCtruthpt[centbin]->Fill(partMC->Pt());
 MCtrutheta[centbin]->Fill(partMC->Eta());
 MCtruthphi[centbin]->Fill(partMC->Phi());

 //get particle ID
Int_t pdgtruth=((AliAODMCParticle*)partMC)->GetPdgCode();

 if (TMath::Abs(pdgtruth)==211)
   {
MCtruthpionpt[centbin]->Fill(partMC->Pt());
 MCtruthpioneta[centbin]->Fill(partMC->Eta());
 MCtruthpionphi[centbin]->Fill(partMC->Phi());
      }

 if (TMath::Abs(pdgtruth)==321)
   {
 MCtruthkaonpt[centbin]->Fill(partMC->Pt());
 MCtruthkaoneta[centbin]->Fill(partMC->Eta());
 MCtruthkaonphi[centbin]->Fill(partMC->Phi());
  }
if(TMath::Abs(pdgtruth)==2212)
  {
 MCtruthprotonpt[centbin]->Fill(partMC->Pt());
 MCtruthprotoneta[centbin]->Fill(partMC->Eta());
 MCtruthprotonphi[centbin]->Fill(partMC->Phi());
  }

 // -- Fill THnSparse
 Double_t primmctruth[4] = {cent_v0m, zVtxmc,partMC->Pt(), partMC->Eta()};

  fTHngenprimPidTruth[5]->Fill(primmctruth);//for all

if (TMath::Abs(pdgtruth)==211 || TMath::Abs(pdgtruth)==321)
  {
    fTHngenprimPidTruth[4]->Fill(primmctruth);//for mesons
  if (TMath::Abs(pdgtruth)==211)  fTHngenprimPidTruth[0]->Fill(primmctruth);//for pions
  if (TMath::Abs(pdgtruth)==321)  fTHngenprimPidTruth[1]->Fill(primmctruth);//for kaons
  }
 else if(TMath::Abs(pdgtruth)==2212)  fTHngenprimPidTruth[2]->Fill(primmctruth);//for protons
 else fTHngenprimPidTruth[3]->Fill(primmctruth);//for others
 
  }//MC truth track loop ends

fEventCounter->Fill(7);


PostData(1, fOutput);


}
//________________________________________________________________________
void AliTwoParticlePIDCorr::doAODevent() 
{

  //get AOD
  AliVEvent *event = InputEvent();
  if (!event) { Printf("ERROR: Could not retrieve event"); return; }
  AliAODEvent* aod = dynamic_cast<AliAODEvent*>(event);
  if (!aod) {
    AliError("Cannot get the AOD event");
    return;
  }

// count all events   
  fEventCounter->Fill(1);

  // get centrality object and check quality
  Float_t cent_v0m=0;
  AliCentrality *centrality=0;
  if(aod) 
 centrality = aod->GetHeader()->GetCentralityP();
  //if(!centrality) return;
  // if (centrality->GetQuality() != 0) return ;

if(centrality)
  { 
  cent_v0m = centrality->GetCentralityPercentile("V0A");
  }
 else
    {
      cent_v0m = -1;
     }

  if (!fPID) return;

if (cent_v0m < 0)  return;
//do centrality binning
 Int_t centbin=Getcentbin(cent_v0m);
 if(centbin==-999) return;

 fhistcentrality->Fill(cent_v0m);

//count events having centrality betn 0-100%
  fEventCounter->Fill(2);

  // Pileup selection ************************************************

  //if(fUtils->IsPileUpEvent(aod)) return;  //applicable only for TPC only tracks,not for hybrid tracks

  
  //vertex selection
  AliAODVertex* trkVtx = aod->GetPrimaryVertex();
  if (!trkVtx || trkVtx->GetNContributors()<=0) return;
  TString vtxTtl = trkVtx->GetTitle();
  if (!vtxTtl.Contains("VertexerTracks")) return;
  Float_t zvtx = trkVtx->GetZ();
  const AliAODVertex* spdVtx = aod->GetPrimaryVertexSPD();
  if (!spdVtx || spdVtx->GetNContributors()<=0) return;
  TString vtxTyp = spdVtx->GetTitle();
  Double_t cov[6]={0};
  spdVtx->GetCovarianceMatrix(cov);
  Double_t zRes = TMath::Sqrt(cov[5]);
  if (vtxTyp.Contains("vertexer:Z") && (zRes>0.25)) return;
   if (TMath::Abs(spdVtx->GetZ() - trkVtx->GetZ())>0.5) return;

  fHistQA[0]->Fill((trkVtx->GetX()));fHistQA[1]->Fill((trkVtx->GetY()));fHistQA[2]->Fill((trkVtx->GetZ()));   //for trkVtx only before vertex cut |zvtx|<10 cm

//count events having a proper vertex
   fEventCounter->Fill(5);

 if (TMath::Abs(zvtx) > 10) return;

//count events after vertex cut
  fEventCounter->Fill(6);


  //if(!fUtils->IsVertexSelected2013pA(aod)) return;
  
 fHistQA[3]->Fill((trkVtx->GetX()));fHistQA[4]->Fill((trkVtx->GetY()));fHistQA[5]->Fill((trkVtx->GetZ()));//after vertex cut,for trkVtx only

//do zvtx binning
  Int_t vtx=Getzbin(zvtx);
  if(vtx==-999) return;//all the events outside the defined vtx range will be ignored
    
  eventno++;

 if(!aod) return; //for safety
  
 Float_t bSign = (aod->GetMagneticField() > 0) ? 1 : -1;//for two track efficiency cut
 Float_t bSign1=aod->GetHeader()->GetMagneticField() ;//for dca cut


   TObjArray* trackstrig = new TObjArray;
   TObjArray* tracksasso = new TObjArray;
   trackstrig->SetOwner(kTRUE);  // IMPORTANT!
   tracksasso->SetOwner(kTRUE);  // IMPORTANT!
 
   Int_t nooftracks=0;
// loop over tracks 
  for (Int_t itrk = 0; itrk < aod->GetNumberOfTracks(); itrk++) 
{ //track loop starts
  AliAODTrack* track = dynamic_cast<AliAODTrack*>(aod->GetTrack(itrk));
  if (!track) continue;
  fHistQA[11]->Fill(track->GetTPCNcls());
  Int_t tracktype=ClassifyTrack(track,centbin,trkVtx,bSign1,kFALSE);//dcacut=kFALSE
  if(tracktype==0) continue; 
  if(tracktype==1)//tracks not passed AliAODTrack::kPrimary at reconstructed level & have proper TPC PID response
{
  //accepted tracks having pt (0.2 to 9.0 Gev) & eta (-0.8 to 0.8)

  nooftracks++;
Float_t pt=track->Pt();
 
Float_t dEdx = track->GetTPCsignal();
fHistoTPCdEdx->Fill(pt, dEdx);

 if(HasTOFPID(track))
{
Float_t beta = GetBeta(track);
fHistoTOFbeta->Fill(pt, beta);
 }

Int_t pt1=Getptbin(track->Pt());//trig--0,1,2,3; asso--4,5,6,7,8
if(pt1==-999) continue;//remove particles with pt<1.0 Gev/c, now only particleswithin 1.0<=pt<=4.0Gev/c are present

//track identification(using nsigma method)
Int_t particletype=GetParticle(track,centbin);

//Fill the nsigma histograms
   Float_t nsigmaTPCPion    =  fnsigmas[SpPion][NSigmaTPC];
   Float_t nsigmaTOFPion    =  fnsigmas[SpPion][NSigmaTOF];
   Float_t nsigmaTPCTOFPion    = TMath::Abs(fnsigmas[SpPion][NSigmaTPCTOF]);//it may be negative if TOF pid is not available(i.e. this is basically nsigmaTPCPionmc)

if (particletype==SpPion)//pions like
   {
     fhistopionnsigmaTPCMC[centbin][pt1]->Fill(nsigmaTPCPion); 
     fhistopionnsigmaTOFMC[centbin][pt1]->Fill(nsigmaTOFPion);
     fhistopionnsigmaTPCTOFMC[centbin][pt1]->Fill(nsigmaTPCTOFPion); 
   }
 if(particletype==SpKaon)//kaons like 
    {
     fhistokaonnsigmaTPCMC[centbin][pt1]->Fill(nsigmaTPCPion); 
     fhistokaonnsigmaTOFMC[centbin][pt1]->Fill(nsigmaTOFPion);
     fhistokaonnsigmaTPCTOFMC[centbin][pt1]->Fill(nsigmaTPCTOFPion); 
    }
if(particletype==SpProton)//protons like
    {
     fhistoprotonnsigmaTPCMC[centbin][pt1]->Fill(nsigmaTPCPion); 
     fhistoprotonnsigmaTOFMC[centbin][pt1]->Fill(nsigmaTOFPion);
     fhistoprotonnsigmaTPCTOFMC[centbin][pt1]->Fill(nsigmaTPCTOFPion); 
    }

//2-d TPCTOF map(for each specified Pt intervals)
 if (HasTOFPID(track)) fTPCTOFpion2d[pt1]->Fill(nsigmaTPCPion,nsigmaTOFPion);
   

 Float_t effmatrix=1.;

 if ((pt>=1.0) && (pt<2.0)) 
	  {
        //if(particletype==kSpUndefined) continue;      
        fEtaSpectrasso->Fill(track->Eta(),track->Pt());
	fphiSpectraasso->Fill(track->Phi(),track->Pt());
if (applyefficiency)
  effmatrix=GetTrackbyTrackeffvalue(track,cent_v0m,zvtx,particletype);
 LRCParticlePID* copy = new LRCParticlePID(particletype,track->Charge(),pt,track->Eta(), track->Phi(),centbin,vtx,effmatrix);
       copy->SetUniqueID(track->GetUniqueID());
       tracksasso->Add(copy);
    //kSpUndefined particles are also added

if(particletype==SpPion) pionassoeta[centbin]->Fill(track->Eta());
if(particletype==SpKaon) kaonassoeta[centbin]->Fill(track->Eta());
if(particletype==SpProton) protonassoeta[centbin]->Fill(track->Eta());	
	  }

if ((pt>=2.0)&&(pt<=4.0))  
          {
   //if(particletype==kSpUndefined) continue; 
   fEtaSpectraTrigall->Fill(track->Eta());//to know the number of trigger particls
if (applyefficiency)
 effmatrix=GetTrackbyTrackeffvalue(track,cent_v0m,zvtx,particletype);
//cout<<effmatrix<<endl;
 LRCParticlePID* copy1 = new LRCParticlePID(particletype,track->Charge(),pt,track->Eta(), track->Phi(),centbin,vtx,effmatrix);
    copy1->SetUniqueID(track->GetUniqueID());
    trackstrig->Add(copy1);
if(particletype==SpProton) baryontrigeta[centbin]->Fill(track->Eta());
if(particletype==SpPion || particletype==SpKaon) mesontrigeta[centbin]->Fill(track->Eta());
	  }
   }// selected particle condition structure ends

} //track loop ends but still in event loop

 fCentralityCorrelation->Fill(cent_v0m, nooftracks);
  
Bool_t isbaryontrig=kFALSE;
Bool_t ismesontrig=kFALSE;

//same event delta-eta-deltaphi plot

if(trackstrig && tracksasso && trackstrig->GetEntriesFast()>0)
  {//same event calculation starts
//calculate no of events in each centrality for different zvtx bins
    fEventno[centbin][vtx]->Fill(5);//only those events which have at least one trigger particle
    Fillcorrelation(trackstrig,tracksasso,centbin,vtx,bSign,kTRUE,kFALSE);//mixcase=kFALSE 

for(Int_t i=0;i<trackstrig->GetEntriesFast();i++)
    {  //trigger loop starts
      LRCParticlePID *trig=(LRCParticlePID*)(trackstrig->UncheckedAt(i));
      if(!trig) continue;
      Int_t particlepidtrig=trig->getparticle(); //either 1 or 2
      if(particlepidtrig==SpProton) isbaryontrig=kTRUE;
      if(particlepidtrig==SpPion || particlepidtrig==SpKaon) ismesontrig=kTRUE;
    }//trig loop ends
 if (isbaryontrig) fEventnobaryon[centbin][vtx]->Fill(5); 
 if (ismesontrig) fEventnomeson[centbin][vtx]->Fill(5);
  }//same event calculation ends

//still in  main event loop
//start mixing
AliEventPool* pool = fPoolMgr->GetEventPool(cent_v0m, zvtx);
if (!pool)
AliFatal(Form("No pool found for centrality = %f, zVtx = %f", cent_v0m, zvtx));
if (pool->IsReady())
  {//start mixing only when pool->IsReady
if(trackstrig && trackstrig->GetEntriesFast()>0)
  {//proceed only when no. of trigger particles >0 in current event
    //TObjArray* bgTracks =new TObjArray();
for (Int_t jMix=0; jMix<pool->GetCurrentNEvents(); jMix++) 
  { //pool event loop start
 TObjArray* bgTracks = pool->GetEvent(jMix);
  if(!bgTracks) continue;
  Fillcorrelation(trackstrig,bgTracks,centbin,vtx,bSign,kTRUE,kTRUE);//mixcase=kTRUE
  
   }// pool event loop ends mixing case
 }//if(trackstrig && trackstrig->GetEntriesFast()>0) condition ends mixing case
} //if pool->IsReady() condition ends mixing case

//still in main event loop
if(pool)
 {
if(tracksasso)
pool->UpdatePool(tracksasso);//ownership of tracksasso is with pool now, don't delete it
 }
if(trackstrig)
delete trackstrig;

fEventCounter->Fill(7);


PostData(1, fOutput);

} // *************************event loop ends******************************************//_______________________________________________________________________

//--------------------------------------------------------------------------------
void AliTwoParticlePIDCorr::Fillcorrelation(TObjArray *trackstrig,TObjArray *tracksasso,Int_t centbin,Int_t vtx,Float_t bSign,Bool_t twoTrackEfficiencyCut,Bool_t mixcase)
{

  //before calling this function check that both trackstrig & tracksasso are available 

 // Eta() is extremely time consuming, therefore cache it for the inner loop here:
  TObjArray* input = (tracksasso) ? tracksasso : trackstrig;
  TArrayF eta(input->GetEntriesFast());
  for (Int_t i=0; i<input->GetEntriesFast(); i++)
    eta[i] = ((LRCParticlePID*) input->UncheckedAt(i))->Eta();
  

// identify K, Lambda candidates and flag those particles
    // a TObject bit is used for this
const UInt_t kResonanceDaughterFlag = 1 << 14;
    if (fRejectResonanceDaughters > 0)
    {
      Double_t resonanceMass = -1;
      Double_t massDaughter1 = -1;
      Double_t massDaughter2 = -1;
      const Double_t interval = 0.02;
 switch (fRejectResonanceDaughters)
      {
	case 1: resonanceMass = 0.9; massDaughter1 = 0.1396; massDaughter2 = 0.9383; break; // method test
	case 2: resonanceMass = 0.4976; massDaughter1 = 0.1396; massDaughter2 = massDaughter1; break; // k0
	case 3: resonanceMass = 1.115; massDaughter1 = 0.1396; massDaughter2 = 0.9383; break; // lambda
	default: AliFatal(Form("Invalid setting %d", fRejectResonanceDaughters));
      }      

for (Int_t i=0; i<trackstrig->GetEntriesFast(); i++)
	trackstrig->UncheckedAt(i)->ResetBit(kResonanceDaughterFlag);
 for (Int_t i=0; tracksasso->GetEntriesFast(); i++)
	  tracksasso->UncheckedAt(i)->ResetBit(kResonanceDaughterFlag);

 for (Int_t i=0; i<trackstrig->GetEntriesFast(); i++)
      {
      LRCParticlePID *trig=(LRCParticlePID*)(trackstrig->UncheckedAt(i));
for (Int_t j=0; tracksasso->GetEntriesFast(); j++)
	{
        LRCParticlePID *asso=(LRCParticlePID*)(tracksasso->UncheckedAt(j));

 // check if both particles point to the same element (does not occur for mixed events, but if subsets are mixed within the same event)
if (trig->IsEqual(asso)) continue;

if (trig->Charge() * asso->Charge() > 0) continue;

 Float_t mass = GetInvMassSquaredCheap(trig->Pt(), trig->Eta(), trig->Phi(), asso->Pt(), asso->Eta(), asso->Phi(), massDaughter1, massDaughter2);
     
if (TMath::Abs(mass - resonanceMass*resonanceMass) < interval*5)
	  {
	    mass = GetInvMassSquared(trig->Pt(), trig->Eta(), trig->Phi(), asso->Pt(), asso->Eta(), asso->Phi(), massDaughter1, massDaughter2);

	    if (mass > (resonanceMass-interval)*(resonanceMass-interval) && mass < (resonanceMass+interval)*(resonanceMass+interval))
	    {
	      trig->SetBit(kResonanceDaughterFlag);
	      asso->SetBit(kResonanceDaughterFlag);
	      
// 	      Printf("Flagged %d %d %f", i, j, TMath::Sqrt(mass));
	    }
	  }
	}
      }
    }

    //two particle correlation filling

for(Int_t i=0;i<trackstrig->GetEntriesFast();i++)
    {  //trigger loop starts
      LRCParticlePID *trig=(LRCParticlePID*)(trackstrig->UncheckedAt(i));
      if(!trig) continue;
      Float_t trigpt=trig->Pt();
      if(trigpt<2.0 || trigpt>4.0) continue;//for safety
      Float_t trigeta=trig->Eta();

      // some optimization
 if (fTriggerRestrictEta > 0 && TMath::Abs(trigeta) > fTriggerRestrictEta)
	continue;

if (fOnlyOneEtaSide != 0)
      {
	if (fOnlyOneEtaSide * trigeta < 0)
	  continue;
      }
  if (fTriggerSelectCharge != 0)
	if (trig->Charge() * fTriggerSelectCharge < 0)
	  continue;
	
      if (fRejectResonanceDaughters > 0)
	if (trig->TestBit(kResonanceDaughterFlag)) continue;

      Float_t trigphi=trig->Phi();
      Int_t particlepidtrig=trig->getparticle(); //either 1 or 2
      Float_t trackefftrig=trig->geteffcorrectionval();

      Int_t pttrig2=Getptbin(trigpt);//pt binning of trigger particles:0,1,2,3
      if(pttrig2==-999) continue;//neglect particles outside the defined range of trigger particles(2 to 4 Gev)
      //filling only for same event case(mixcase=kFALSE)
if(mixcase==kFALSE)  
  {
    fEtaSpectraTrig[centbin][vtx]->Fill(trigeta,1./trackefftrig);
    if(particlepidtrig==SpProton)//proton
      fEtaSpectraTrigbaryon[centbin][vtx]->Fill(trigeta,1./trackefftrig);
      if(particlepidtrig==SpPion || particlepidtrig==SpKaon)//kaon+pion
	fEtaSpectraTrigmeson[centbin][vtx]->Fill(trigeta,1./trackefftrig);
  }
    //asso loop starts within trigger loop
   for(Int_t j=0;j<tracksasso->GetEntriesFast();j++)
             {
        LRCParticlePID *asso=(LRCParticlePID*)(tracksasso->UncheckedAt(j));
	if(!asso) continue;
Float_t assoPt=asso->Pt();
if(assoPt<1.0 || assoPt>=2.0) continue;//for safety

if (trig->IsEqual(asso)) continue;
 
if (fPtOrder)
 if (asso->Pt() >= trig->Pt()) continue;
	    

if (fAssociatedSelectCharge != 0)
if (asso->Charge() * fAssociatedSelectCharge < 0) continue;
	    
 if (fSelectCharge > 0)
        {
          // skip like sign
          if (fSelectCharge == 1 && asso->Charge() * trig->Charge() > 0)
            continue;
            
          // skip unlike sign
          if (fSelectCharge == 2 && asso->Charge() * trig->Charge() < 0)
            continue;
        }

if (fEtaOrdering)
	{
	  if (trigeta < 0 && asso->Eta() < trigeta)
	    continue;
	  if (trigeta > 0 && asso->Eta() > trigeta)
	    continue;
	}

if (fRejectResonanceDaughters > 0)
	  if (asso->TestBit(kResonanceDaughterFlag))
	  {
// 	    Printf("Skipped j=%d", j);
	    continue;
	  }

	// conversions
	if (fCutConversions && asso->Charge() * trig->Charge() < 0)
	{
	  Float_t mass = GetInvMassSquaredCheap(trig->Pt(), trigeta, trig->Phi(), asso->Pt(),eta[j], asso->Phi(), 0.510e-3, 0.510e-3);
	  
	  if (mass < 0.1)
	  {
	    mass = GetInvMassSquared(trig->Pt(), trigeta, trig->Phi(), asso->Pt(), eta[j], asso->Phi(), 0.510e-3, 0.510e-3);
	    
	    //fControlConvResoncances->Fill(0.0, mass);

	    if (mass < 0.04*0.04) 
	      continue;
	  }
	}

	// K0s
	if (fCutResonances && asso->Charge() * trig->Charge() < 0)
	{
	  Float_t mass = GetInvMassSquaredCheap(trig->Pt(), trigeta, trig->Phi(), asso->Pt(), eta[j], asso->Phi(), 0.1396, 0.1396);
	  
	  const Float_t kK0smass = 0.4976;
	  
	  if (TMath::Abs(mass - kK0smass*kK0smass) < 0.1)
	  {
	    mass = GetInvMassSquared(trig->Pt(), trigeta, trig->Phi(), asso->Pt(),eta[j], asso->Phi(), 0.1396, 0.1396);
	    
	    //fControlConvResoncances->Fill(1, mass - kK0smass*kK0smass);

	    if (mass > (kK0smass-0.02)*(kK0smass-0.02) && mass < (kK0smass+0.02)*(kK0smass+0.02))
	      continue;
	  }
	}
	
	// Lambda
	if (fCutResonances && asso->Charge() * trig->Charge() < 0)
	{
	  Float_t mass1 = GetInvMassSquaredCheap(trig->Pt(), trigeta, trig->Phi(), asso->Pt(), eta[j], asso->Phi(), 0.1396, 0.9383);
	  Float_t mass2 = GetInvMassSquaredCheap(trig->Pt(), trigeta, trig->Phi(), asso->Pt(),eta[j] , asso->Phi(), 0.9383, 0.1396);
	  
	  const Float_t kLambdaMass = 1.115;

	  if (TMath::Abs(mass1 - kLambdaMass*kLambdaMass) < 0.1)
	  {
	    mass1 = GetInvMassSquared(trig->Pt(), trigeta, trig->Phi(), asso->Pt(),eta[j], asso->Phi(), 0.1396, 0.9383);

	    //fControlConvResoncances->Fill(2, mass1 - kLambdaMass*kLambdaMass);
	    
	    if (mass1 > (kLambdaMass-0.02)*(kLambdaMass-0.02) && mass1 < (kLambdaMass+0.02)*(kLambdaMass+0.02))
	      continue;
	  }
	  if (TMath::Abs(mass2 - kLambdaMass*kLambdaMass) < 0.1)
	  {
	    mass2 = GetInvMassSquared(trig->Pt(), trigeta, trig->Phi(), asso->Pt(),eta[j] , asso->Phi(), 0.9383, 0.1396);

	    //fControlConvResoncances->Fill(2, mass2 - kLambdaMass*kLambdaMass);

	    if (mass2 > (kLambdaMass-0.02)*(kLambdaMass-0.02) && mass2 < (kLambdaMass+0.02)*(kLambdaMass+0.02))
	      continue;
	  }
	}

	if (twoTrackEfficiencyCut)
	{
	  // the variables & cuthave been developed by the HBT group 
	  // see e.g. https://indico.cern.ch/materialDisplay.py?contribId=36&sessionId=6&materialId=slides&confId=142700
	  Float_t phi1 = trig->Phi();
	  Float_t pt1 = trig->Pt();
	  Float_t charge1 = trig->Charge();
	  Float_t phi2 = asso->Phi();
	  Float_t pt2 = asso->Pt();
	  Float_t charge2 = asso->Charge();

	  Float_t deta= trigeta - eta[j]; 
    
 // optimization
	  if (TMath::Abs(deta) < twoTrackEfficiencyCutValue * 2.5 * 3)
	  {

  // check first boundaries to see if is worth to loop and find the minimum
	    Float_t dphistar1 = GetDPhiStar(phi1, pt1, charge1, phi2, pt2, charge2, 0.8, bSign);
	    Float_t dphistar2 = GetDPhiStar(phi1, pt1, charge1, phi2, pt2, charge2, 2.5, bSign);

 const Float_t kLimit = twoTrackEfficiencyCutValue * 3;

	    Float_t dphistarminabs = 1e5;
	    Float_t dphistarmin = 1e5;

 if (TMath::Abs(dphistar1) < kLimit || TMath::Abs(dphistar2) < kLimit || dphistar1 * dphistar2 < 0)
	    {
	      for (Double_t rad=0.8; rad<2.51; rad+=0.01) 
	      {
		Float_t dphistar = GetDPhiStar(phi1, pt1, charge1, phi2, pt2, charge2, rad, bSign);

		Float_t dphistarabs = TMath::Abs(dphistar);

	if (dphistarabs < dphistarminabs)
		{
		  dphistarmin = dphistar;
		  dphistarminabs = dphistarabs;
		}
	      }

if (dphistarminabs < twoTrackEfficiencyCutValue && TMath::Abs(deta) < twoTrackEfficiencyCutValue)
	      {
// 		Printf("Removed track pair %d %d with %f %f %f %f %f %f %f %f %f", i, j, deta, dphistarminabs, phi1, pt1, charge1, phi2, pt2, charge2, bSign);
		continue;
	      }
//fTwoTrackDistancePt[1]->Fill(deta, dphistarmin, TMath::Abs(pt1 - pt2));

	    }
	  }
	}

        Int_t particlepidasso=asso->getparticle(); //either 10 or 20 or 30
       Float_t trackeffasso=asso->geteffcorrectionval();

        Float_t deleta=trigeta-eta[j];
	Float_t delphi=PhiRange(trigphi-asso->Phi()); 
        Int_t ptbin=(Int_t)asso->Pt(); //asso pt binning either 0 or 1


	//here get the two particle efficiency correction factor
	Float_t effweight=trackefftrig*trackeffasso;
	//	cout<<effweight<<endl;
//same event calculation(mixcase=kFALSE)
if(mixcase==kFALSE)
  {	
fdeletasame->Fill(deleta,1./effweight);
fdelphisame->Fill(delphi,1./effweight);
 fsame->Fill(delphi,deleta,1./effweight);
// correlation histograms same event
 falltrigallasso[centbin][ptbin][vtx]->Fill(delphi,deleta,1./effweight);
if (particlepidasso==SpPion ) falltrigpionasso[centbin][ptbin][vtx]->Fill(delphi,deleta,1./effweight);
if (particlepidasso==SpKaon ) falltrigkaonasso[centbin][ptbin][vtx]->Fill(delphi,deleta,1./effweight);
if (particlepidasso==SpProton ) falltrigprotonasso[centbin][ptbin][vtx]->Fill(delphi,deleta,1./effweight);
 if(particlepidtrig==SpProton)//proton trig
   { 
fbaryontrigallasso[centbin][ptbin][vtx]->Fill(delphi,deleta,1./effweight);
if (particlepidasso==SpPion ) 
fbaryontrigpionasso[centbin][ptbin][vtx]->Fill(delphi,deleta,1./effweight);
if (particlepidasso==SpKaon ) 
fbaryontrigkaonasso[centbin][ptbin][vtx]->Fill(delphi,deleta,1./effweight);
if (particlepidasso==SpProton ) 
fbaryontrigprotonasso[centbin][ptbin][vtx]->Fill(delphi,deleta,1./effweight);
   }
 if(particlepidtrig==SpPion || particlepidtrig==SpKaon)//meson trig
   {
fmesontrigallasso[centbin][ptbin][vtx]->Fill(delphi,deleta,1./effweight); 
if (particlepidasso==SpPion) 
fmesontrigpionasso[centbin][ptbin][vtx]->Fill(delphi,deleta,1./effweight);
if (particlepidasso==SpKaon) 
fmesontrigkaonasso[centbin][ptbin][vtx]->Fill(delphi,deleta,1./effweight);
if (particlepidasso==SpProton) 
fmesontrigprotonasso[centbin][ptbin][vtx]->Fill(delphi,deleta,1./effweight);
   }
   }
 if(mixcase==kTRUE)//for mixing case
   {
fmix->Fill(delphi,deleta,1./effweight);
fdeletamixed->Fill(deleta,1./effweight);
fdelphimixed->Fill(delphi,1./effweight);
falltrigallassomix[centbin][ptbin][vtx]->Fill(delphi,deleta,1./effweight);
if (particlepidasso==SpPion ) falltrigpionassomix[centbin][ptbin][vtx]->Fill(delphi,deleta,1./effweight);
if (particlepidasso==SpKaon ) falltrigkaonassomix[centbin][ptbin][vtx]->Fill(delphi,deleta,1./effweight);
if (particlepidasso==SpProton ) falltrigprotonassomix[centbin][ptbin][vtx]->Fill(delphi,deleta,1./effweight);
if(particlepidtrig==SpProton)
  {
 fdeletamixedproton->Fill(deleta,1./effweight);
 fdelphimixedproton->Fill(delphi,1./effweight); 
fbaryontrigallassomix[centbin][ptbin][vtx]->Fill(delphi,deleta,1./effweight);
if (particlepidasso==SpPion ) 
fbaryontrigpionassomix[centbin][ptbin][vtx]->Fill(delphi,deleta,1./effweight);
if (particlepidasso==SpKaon ) 
fbaryontrigkaonassomix[centbin][ptbin][vtx]->Fill(delphi,deleta,1./effweight);
if (particlepidasso==SpProton ) 
fbaryontrigprotonassomix[centbin][ptbin][vtx]->Fill(delphi,deleta,1./effweight);
  }
if(particlepidtrig==SpPion || particlepidtrig==SpKaon)
  { 
fdeletamixedkaonpion->Fill(deleta,1./effweight);
fdelphimixedkaonpion->Fill(delphi,1./effweight);
fmesontrigallassomix[centbin][ptbin][vtx]->Fill(delphi,deleta,1./effweight);
if (particlepidasso==SpPion) 
fmesontrigpionassomix[centbin][ptbin][vtx]->Fill(delphi,deleta,1./effweight);
if (particlepidasso==SpKaon) 
fmesontrigkaonassomix[centbin][ptbin][vtx]->Fill(delphi,deleta,1./effweight);
if (particlepidasso==SpProton) 
fmesontrigprotonassomix[centbin][ptbin][vtx]->Fill(delphi,deleta,1./effweight);
 }

   }

   }//asso loop ends 
 }//trigger loop ends 

}

//--------------------------------------------------------------------------------
Float_t AliTwoParticlePIDCorr::GetTrackbyTrackeffvalue(AliAODTrack* track,Float_t cent,Float_t evzvtx, Int_t parpid)
{
  //This function is called only when applyefficiency=kTRUE
 Int_t effVars[4];
 Float_t effvalue=1.; 

  if(parpid==SpUndefined)
            {
	    effVars[0] = effcorection[5]->GetAxis(0)->FindBin(cent);
	    effVars[1] = effcorection[5]->GetAxis(1)->FindBin(evzvtx); 
	    effVars[2] = effcorection[5]->GetAxis(2)->FindBin(track->Pt()); 
	    effVars[3] = effcorection[5]->GetAxis(3)->FindBin(track->Eta()); 
            effvalue=effcorection[5]->GetBinContent(effVars);
	    }
if(parpid==SpPion || parpid==SpKaon)
            {
	      if(track->Pt()>=2.0)
		{
	    effVars[0] = effcorection[4]->GetAxis(0)->FindBin(cent);
	    effVars[1] = effcorection[4]->GetAxis(1)->FindBin(evzvtx); 
	    effVars[2] = effcorection[4]->GetAxis(2)->FindBin(track->Pt()); 
	    effVars[3] = effcorection[4]->GetAxis(3)->FindBin(track->Eta());
            effvalue=effcorection[4]->GetBinContent(effVars);
		}
if(track->Pt()<2.0){
 if(parpid==SpPion)
            {
	    effVars[0] = effcorection[0]->GetAxis(0)->FindBin(cent);
	    effVars[1] = effcorection[0]->GetAxis(1)->FindBin(evzvtx); 
	    effVars[2] = effcorection[0]->GetAxis(2)->FindBin(track->Pt()); 
	    effVars[3] = effcorection[0]->GetAxis(3)->FindBin(track->Eta()); 
            effvalue=effcorection[0]->GetBinContent(effVars);
	    }
	    
 if(parpid==SpKaon)
            {
	    effVars[0] = effcorection[1]->GetAxis(0)->FindBin(cent);
	    effVars[1] =  effcorection[1]->GetAxis(1)->FindBin(evzvtx); 
	    effVars[2] =  effcorection[1]->GetAxis(2)->FindBin(track->Pt()); 
	    effVars[3] =  effcorection[1]->GetAxis(3)->FindBin(track->Eta()); 
            effvalue=effcorection[1]->GetBinContent(effVars);
	    }
		}
	    }  
 if(parpid==SpProton)
            {
	    effVars[0] =  effcorection[2]->GetAxis(0)->FindBin(cent);
	    effVars[1] = effcorection[2]->GetAxis(1)->FindBin(evzvtx); 
	    effVars[2] = effcorection[2]->GetAxis(2)->FindBin(track->Pt()); 
	    effVars[3] = effcorection[2]->GetAxis(3)->FindBin(track->Eta()); 
            effvalue=effcorection[2]->GetBinContent(effVars);
	    }	    
	    // 	  Printf("%d %d %d %d %f", effVars[0], effVars[1], effVars[2], effVars[3], fEfficiencyCorrectionAssociated->GetBinContent(effVars));
     if(effvalue==0.) effvalue=1.;

     return effvalue; 

}
//-----------------------------------------------------------------------

Int_t AliTwoParticlePIDCorr::ClassifyTrack(AliAODTrack* track,Int_t centbin,AliAODVertex* vertex,Float_t magfield,Bool_t dcacut)
{  
 
  if(!track) return 0;
  Bool_t trackOK = track->TestFilterBit(768);
  if(!trackOK) return 0;
  //select only primary traks 
  //if(track->GetType()!=AliAODTrack::kPrimary) return 0;
  if(track->Charge()==0) return 0;
  fHistQA[12]->Fill(track->GetTPCNcls());  
  Float_t dxy, dz;		  
  dxy = track->DCA();
  dz = track->ZAtDCA();
  fHistQA[6]->Fill(dxy);
  fHistQA[7]->Fill(dz);
  Float_t chi2ndf = track->Chi2perNDF();
  fHistQA[13]->Fill(chi2ndf);  
  Float_t nCrossedRowsTPC = track->GetTPCClusterInfo(2,1);
  fHistQA[14]->Fill(nCrossedRowsTPC); 
  //Float_t  ratioCrossedRowsOverFindableClustersTPC = 1.0;
  if (track->GetTPCNclsF()>0) {
   Float_t  ratioCrossedRowsOverFindableClustersTPC = nCrossedRowsTPC/track->GetTPCNclsF();
   fHistQA[15]->Fill(ratioCrossedRowsOverFindableClustersTPC);
    }
//accepted tracks  
     Float_t pt=track->Pt();
     if(pt<0.2 || pt>9.0) return 0;
     if(TMath::Abs(track->Eta())>0.8) return 0;
     if(track->Phi()<0. || track->Phi()>2*TMath::Pi()) return 0;
     //if (!HasTPCPID(track)) return 0;//trigger & associated particles must have TPC PID
// DCA XY
	if (dcacut && fDCAXYCut)
	{
	  if (!vertex)
	    return 0;
	  
	  Double_t pos[2];
	  Double_t covar[2];
	  AliAODTrack* clone =(AliAODTrack*) track->Clone();
	  Bool_t success = clone->PropagateToDCA(vertex, magfield, 3, pos, covar);
	  delete clone;
	  if (!success)
	    return 0;

// 	  Printf("%f", ((AliAODTrack*)part)->DCA());
// 	  Printf("%f", pos[0]);
	  if (TMath::Abs(pos[0]) > fDCAXYCut->Eval(track->Pt()))
	    return 0;
	}

     fHistQA[8]->Fill(pt);
     recoallpt[centbin]->Fill(track->Pt());
     fHistQA[9]->Fill(track->Eta());
     recoalleta[centbin]->Fill(track->Eta());
     fHistQA[10]->Fill(track->Phi());
     recoallphi[centbin]->Fill(track->Phi());
     if(pt>=1.0 && pt<2.0) allassoeta[centbin]->Fill(track->Eta());
     if(pt>=2.0 && pt<=4.0) alltrigeta[centbin]->Fill(track->Eta());
     return 1;
  }
  //________________________________________________________________________________
void AliTwoParticlePIDCorr::CalculateNSigmas(AliAODTrack *track,Int_t centbin) 
{
//This function is called within the func GetParticle() for accepted tracks only i.e.after call of Classifytrack() & for those tracks which have proper TPC PID response & only for particles having pt within 1.0 to 4.0 Gev/c i.e. while filling the the TObjArray for trig & asso 
Float_t pt=track->Pt();

//it is assumed that every track that passed the filterbit have proper TPC response(!!)
Float_t nsigmaTPCkPion =fPID->NumberOfSigmasTPC(track, AliPID::kPion);
Float_t nsigmaTPCkKaon =fPID->NumberOfSigmasTPC(track, AliPID::kKaon);
Float_t nsigmaTPCkProton =fPID->NumberOfSigmasTPC(track, AliPID::kProton);

Float_t nsigmaTOFkProton=999.,nsigmaTOFkKaon=999.,nsigmaTOFkPion=999.;
Float_t nsigmaTPCTOFkProton=999.,nsigmaTPCTOFkKaon=999.,nsigmaTPCTOFkPion=999.;

 if(HasTOFPID(track) && pt>fPtTOFPID)
   {

nsigmaTOFkPion =fPID->NumberOfSigmasTOF(track, AliPID::kPion);
nsigmaTOFkKaon =fPID->NumberOfSigmasTOF(track, AliPID::kKaon);
nsigmaTOFkProton =fPID->NumberOfSigmasTOF(track, AliPID::kProton);
//---combined
nsigmaTPCTOFkPion   = TMath::Sqrt(nsigmaTPCkPion*nsigmaTPCkPion+nsigmaTOFkPion*nsigmaTOFkPion);
nsigmaTPCTOFkKaon   = TMath::Sqrt(nsigmaTPCkKaon*nsigmaTPCkKaon+nsigmaTOFkKaon*nsigmaTOFkKaon);
nsigmaTPCTOFkProton = TMath::Sqrt(nsigmaTPCkProton*nsigmaTPCkProton+nsigmaTOFkProton*nsigmaTOFkProton);

//fill the nsigma pion histograms
 if(pt>=1.0 && pt<=4.0)
   {
Int_t ptbinval=Getptbin(track->Pt());//trig->0,1,2,3---ass0->4,5,6,7,8       
fHistoNSigmaTPCpion[centbin][ptbinval]->Fill(nsigmaTPCkPion);
fHistoNSigmaTOFpion[centbin][ptbinval]->Fill(nsigmaTOFkPion);
fHistoNSigmaTPCTOFpion[centbin][ptbinval]->Fill(nsigmaTPCTOFkPion);
//fHistocentNSigmaTPC->Fill(centrality, pt, nsigmaTPCkPion);
//fHistocentNSigmaTOF->Fill(centrality, pt, nsigmaTOFkPion);
   }
   }
else{
    // --- combined
    // if TOF is missing and below fPtTOFPID only the TPC information is used
    nsigmaTPCTOFkProton = TMath::Abs(nsigmaTPCkProton);
    nsigmaTPCTOFkKaon   = TMath::Abs(nsigmaTPCkKaon);
    nsigmaTPCTOFkPion   = TMath::Abs(nsigmaTPCkPion);
if(pt>=1.0 && pt<=4.0)
   {
     Int_t ptbinval=Getptbin(track->Pt());//trig->0,1,2,3---ass0->4,5,6,7,8       
     fHistoNSigmaTPCpion[centbin][ptbinval]->Fill(nsigmaTPCkPion);
     fHistoNSigmaTPCTOFpion[centbin][ptbinval]->Fill(nsigmaTPCTOFkPion);
   //fHistocentNSigmaTPC->Fill(centrality, pt, nsigmaTPCkPion);
   }
  }

//set data member fnsigmas
  fnsigmas[SpPion][NSigmaTPC]=nsigmaTPCkPion;
  fnsigmas[SpKaon][NSigmaTPC]=nsigmaTPCkKaon;
  fnsigmas[SpProton][NSigmaTPC]=nsigmaTPCkProton;
  fnsigmas[SpPion][NSigmaTOF]=nsigmaTOFkPion;
  fnsigmas[SpKaon][NSigmaTOF]=nsigmaTOFkKaon;
  fnsigmas[SpProton][NSigmaTOF]=nsigmaTOFkProton;
  fnsigmas[SpPion][NSigmaTPCTOF]=nsigmaTPCTOFkPion;
  fnsigmas[SpKaon][NSigmaTPCTOF]=nsigmaTPCTOFkKaon;
  fnsigmas[SpProton][NSigmaTPCTOF]=nsigmaTPCTOFkProton;


}
//----------------------------------------------------------------------------
Int_t AliTwoParticlePIDCorr::FindMinNSigma(AliAODTrack *track) 
{
  //this function is always called after calling the function CalculateNSigmas(AliAODTrack *track)
if(fRequestTOFPID && (!HasTOFPID(track)) && track->Pt()>fPtTOFPID)return SpUndefined;
//get the identity of the particle with the minimum Nsigma
  Float_t nsigmaPion=999., nsigmaKaon=999., nsigmaProton=999.;
  switch (fPIDType){
  case NSigmaTPC:
    nsigmaProton  =  TMath::Abs(fnsigmas[SpProton][NSigmaTPC]);
    nsigmaKaon	  =  TMath::Abs(fnsigmas[SpKaon][NSigmaTPC])  ;
    nsigmaPion    =  TMath::Abs(fnsigmas[SpPion][NSigmaTPC])  ;
    break;
  case NSigmaTOF:
    nsigmaProton  =  TMath::Abs(fnsigmas[SpProton][NSigmaTOF]);
    nsigmaKaon	  =  TMath::Abs(fnsigmas[SpKaon][NSigmaTOF])  ;
    nsigmaPion    =  TMath::Abs(fnsigmas[SpPion][NSigmaTOF])  ;
    break;
  case NSigmaTPCTOF://In case of no TOF matching the combined nsigma is the TPC one
    nsigmaProton  =  TMath::Abs(fnsigmas[SpProton][NSigmaTPCTOF]);
    nsigmaKaon	  =  TMath::Abs(fnsigmas[SpKaon][NSigmaTPCTOF])  ;
    nsigmaPion    =  TMath::Abs(fnsigmas[SpPion][NSigmaTPCTOF])  ;
    break;
  }
  
 // guess the particle based on the smaller nsigma (within fNSigmaPID)
  if( ( nsigmaKaon==nsigmaPion ) && ( nsigmaKaon==nsigmaProton )) return SpUndefined;//it is the default value for the three
if( ( nsigmaKaon   < nsigmaPion ) && ( nsigmaKaon < nsigmaProton ) && (nsigmaKaon < fNSigmaPID)) return SpKaon;
if( ( nsigmaPion   < nsigmaKaon ) && ( nsigmaPion < nsigmaProton ) && (nsigmaPion < fNSigmaPID)) return SpPion;
if( ( nsigmaProton < nsigmaKaon ) && ( nsigmaProton < nsigmaPion ) && (nsigmaProton < fNSigmaPID)) return SpProton;

// else, return undefined
  return SpUndefined;

}

//------------------------------------------------------------------------------------------
Bool_t* AliTwoParticlePIDCorr::GetDoubleCounting(AliAODTrack * trk){ 
  //this function is always called after calling the function CalculateNSigmas(AliAODTrack *track)

  //if a particle has double counting set fHasDoubleCounting[ipart]=kTRUE
  //fill DC histos
  for(Int_t ipart=0;ipart<NSpecies;ipart++)fHasDoubleCounting[ipart]=kFALSE;//array with kTRUE for second (or third) identity of the track
  
  Int_t MinNSigma=FindMinNSigma(trk);//not filling the NSigmaRec histos
  
  
  if(MinNSigma==SpUndefined)return fHasDoubleCounting;//in case of undefined no Double counting
  
  Float_t nsigmaPion=999., nsigmaKaon=999., nsigmaProton=999.;
  switch (fPIDType) {
  case NSigmaTPC:
    nsigmaProton  =  TMath::Abs(fnsigmas[SpProton][NSigmaTPC]);
    nsigmaKaon	  =  TMath::Abs(fnsigmas[SpKaon][NSigmaTPC])  ;
    nsigmaPion    =  TMath::Abs(fnsigmas[SpPion][NSigmaTPC])  ;
    break;
  case NSigmaTOF:
    nsigmaProton  =  TMath::Abs(fnsigmas[SpProton][NSigmaTOF]);
    nsigmaKaon	  =  TMath::Abs(fnsigmas[SpKaon][NSigmaTOF])  ;
    nsigmaPion    =  TMath::Abs(fnsigmas[SpPion][NSigmaTOF])  ;
    break;
  case NSigmaTPCTOF://In case of no TOF matching the combined nsigma is the TPC one
    nsigmaProton  =  TMath::Abs(fnsigmas[SpProton][NSigmaTPCTOF]);
    nsigmaKaon	  =  TMath::Abs(fnsigmas[SpKaon][NSigmaTPCTOF])  ;
    nsigmaPion    =  TMath::Abs(fnsigmas[SpPion][NSigmaTPCTOF])  ;
    break;
  }

if(trk->Pt()<2.0)//only associated particles
    {
  if(nsigmaPion<fNSigmaPID && MinNSigma!=SpPion)fHasDoubleCounting[SpPion]=kTRUE;
  if(nsigmaKaon<fNSigmaPID && MinNSigma!=SpKaon)fHasDoubleCounting[SpKaon]=kTRUE;
  if(nsigmaProton<fNSigmaPID && MinNSigma!=SpProton)fHasDoubleCounting[SpProton]=kTRUE;
    } 

 if (trk->Pt()>=2.0 )// just consider overlapping between meson & baryon(proton)
    {
  if(nsigmaPion<fNSigmaPID && MinNSigma==SpProton)fHasDoubleCounting[SpPion]=kTRUE;
  if(nsigmaKaon<fNSigmaPID && MinNSigma==SpProton)fHasDoubleCounting[SpKaon]=kTRUE;
  if(nsigmaProton<fNSigmaPID && MinNSigma!=SpProton)fHasDoubleCounting[SpProton]=kTRUE;
    } 

 
  return fHasDoubleCounting;
}

//////////////////////////////////////////////////////////////////////////////////////////////////
Int_t AliTwoParticlePIDCorr::GetParticle(AliAODTrack * trk,Int_t centbin){ 
  //return the specie according to the minimum nsigma value
  //no double counting, this has to be evaluated using CheckDoubleCounting()
  //Printf("fPtTOFPID %.1f, fRequestTOFPID %d, fNSigmaPID %.1f, fPIDType %d",fPtTOFPID,fRequestTOFPID,fNSigmaPID,fPIDType);
  
  CalculateNSigmas(trk,centbin);//fill the data member fnsigmas with the nsigmas value [ipart][iPID]
  
  if(fUseExclusiveNSigma){
    Bool_t *HasDC;
    HasDC=GetDoubleCounting(trk);
    for(Int_t ipart=0;ipart<NSpecies;ipart++){
      if(HasDC[ipart]==kTRUE)  return SpUndefined;
    }
    return FindMinNSigma(trk);//NSigmaRec distr filled here
  }
  else return FindMinNSigma(trk);//NSigmaRec distr filled here
  
}
//////////////////////////////////////////////////////////////////////////////////////////////////

  Int_t AliTwoParticlePIDCorr::Getzbin(Float_t z)
{
  Int_t z1=-999;
  if(z>=-10. && z<-8.) z1=0;
  if(z>=-8. && z<-6.) z1=1;
  if(z>=-6. && z<-4.) z1=2;
  if(z>=-4. && z<-2.) z1=3;
  if(z>=-2. && z<0.) z1=4;
  if(z>=0. && z<2.) z1=5;
  if(z>=2. && z<4.) z1=6;
  if(z>=4. && z<6.) z1=7;
  if(z>=6. && z<=8.) z1=8;
  if(z>8. && z<=10.) z1=9;
  return z1;
 }
//-------------------------------------------------------------------------------
Int_t AliTwoParticlePIDCorr::Getptbin(Float_t pt)
  {
 Int_t ptval=-999;
 if ((pt>=2.0) && (pt<=2.5))  ptval=0;
 if ((pt>2.5) && (pt<=3.0))  ptval=1;
 if ((pt>3.0) && (pt<=3.5))  ptval=2;
 if ((pt>3.5) && (pt<=4.0))  ptval=3; 
 if ((pt>=1.0) && (pt<=1.2)) ptval=4;
 if ((pt>1.2) && (pt<=1.4))  ptval=5;
 if ((pt>1.4) && (pt<=1.6))  ptval=6;
 if ((pt>1.6) && (pt<=1.8))  ptval=7;
 if ((pt>1.8) && (pt<2.0))  ptval=8;
 return ptval;//will return -999 if the pt value passed as argument is not within the above mentioned range
  }


//------------------------------------------------------------------------------
Int_t AliTwoParticlePIDCorr::Getcentbin(Float_t cent_v0m)
{
Int_t centbinn=-999;
 if(cent_v0m>0. && cent_v0m<=20.) centbinn=0;
 if(cent_v0m>20. && cent_v0m<=40.) centbinn=1;
 if(cent_v0m>40. && cent_v0m<=60.) centbinn=2;
 if(cent_v0m>60. && cent_v0m<=100.) centbinn=3;

 return centbinn;
}

//-------------------------------------------------------------------------------------
Bool_t
AliTwoParticlePIDCorr::HasTPCPID(AliAODTrack *track) const
{  
  // check PID signal 
   AliPIDResponse::EDetPidStatus statustpc = fPID->CheckPIDStatus(AliPIDResponse::kTPC,track);
   if(statustpc!=AliPIDResponse::kDetPidOk) return kFALSE;
   //ULong_t status=track->GetStatus();
   //if  (!( (status & AliAODTrack::kTPCpid  ) == AliAODTrack::kTPCpid  )) return kFALSE;//remove light nuclei
   //if (track->GetTPCsignal() <= 0.) return kFALSE;
   // if(track->GetTPCsignalN() < 60) return kFALSE;//tracks with TPCsignalN< 60 have questionable dEdx,cutting on TPCsignalN > 70 or > 60 shouldn't make too much difference in statistics,also  it is IMO safe to use TPC also for MIPs.
   
  return kTRUE;  
}
//___________________________________________________________

Bool_t
AliTwoParticlePIDCorr::HasTOFPID(AliAODTrack *track) const
{
  // check TOF matched track 
  //ULong_t status=track->GetStatus();
  //if  (!( (status & AliAODTrack::kITSin  ) == AliAODTrack::kITSin  )) return kFALSE;
 AliPIDResponse::EDetPidStatus statustof = fPID->CheckPIDStatus(AliPIDResponse::kTOF,track);
 if(statustof!= AliPIDResponse::kDetPidOk) return kFALSE;
 //if(!((status & AliAODTrack::kTOFpid  ) == AliAODTrack::kTOFpid  )) return kFALSE;
 //Float_t probMis = fPIDresponse->GetTOFMismatchProbability(track);
 // if (probMis > 0.01) return kFALSE;
if(fRemoveTracksT0Fill)
    {
Int_t startTimeMask = fPID->GetTOFResponse().GetStartTimeMask(track->P());
      if (startTimeMask < 0)return kFALSE; 
    }
  return kTRUE;
}

//________________________________________________________________________
Float_t AliTwoParticlePIDCorr :: GetBeta(AliAODTrack *track)
{
  //it is called only when TOF PID is available
  Double_t p = track->P();
  Double_t time=track->GetTOFsignal()-fPID->GetTOFResponse().GetStartTime(p);
  Double_t timei[5];
  track->GetIntegratedTimes(timei);
  return timei[0]/time;
}
//------------------------------------------------------------------------------------------------------

Float_t AliTwoParticlePIDCorr::GetInvMassSquared(Float_t pt1, Float_t eta1, Float_t phi1, Float_t pt2, Float_t eta2, Float_t phi2, Float_t m0_1, Float_t m0_2)
{
  // calculate inv mass squared
  // same can be achieved, but with more computing time with
  /*TLorentzVector photon, p1, p2;
  p1.SetPtEtaPhiM(triggerParticle->Pt(), triggerEta, triggerParticle->Phi(), 0.510e-3);
  p2.SetPtEtaPhiM(particle->Pt(), eta[j], particle->Phi(), 0.510e-3);
  photon = p1+p2;
  photon.M()*/
  
  Float_t tantheta1 = 1e10;
  
  if (eta1 < -1e-10 || eta1 > 1e-10)
    tantheta1 = 2 * TMath::Exp(-eta1) / ( 1 - TMath::Exp(-2*eta1));
  
  Float_t tantheta2 = 1e10;
  if (eta2 < -1e-10 || eta2 > 1e-10)
    tantheta2 = 2 * TMath::Exp(-eta2) / ( 1 - TMath::Exp(-2*eta2));
  
  Float_t e1squ = m0_1 * m0_1 + pt1 * pt1 * (1.0 + 1.0 / tantheta1 / tantheta1);
  Float_t e2squ = m0_2 * m0_2 + pt2 * pt2 * (1.0 + 1.0 / tantheta2 / tantheta2);
  
  Float_t mass2 = m0_1 * m0_1 + m0_2 * m0_2 + 2 * ( TMath::Sqrt(e1squ * e2squ) - ( pt1 * pt2 * ( TMath::Cos(phi1 - phi2) + 1.0 / tantheta1 / tantheta2 ) ) );
  
  return mass2;
}
//---------------------------------------------------------------------------------

Float_t AliTwoParticlePIDCorr::GetInvMassSquaredCheap(Float_t pt1, Float_t eta1, Float_t phi1, Float_t pt2, Float_t eta2, Float_t phi2, Float_t m0_1, Float_t m0_2)
{
  // calculate inv mass squared approximately
  
  Float_t tantheta1 = 1e10;
  
  if (eta1 < -1e-10 || eta1 > 1e-10)
  {
    Float_t expTmp = 1.0-eta1+eta1*eta1/2-eta1*eta1*eta1/6+eta1*eta1*eta1*eta1/24;
    tantheta1 = 2.0 * expTmp / ( 1.0 - expTmp*expTmp);
  }
  
  Float_t tantheta2 = 1e10;
  if (eta2 < -1e-10 || eta2 > 1e-10)
  {
    Float_t expTmp = 1.0-eta2+eta2*eta2/2-eta2*eta2*eta2/6+eta2*eta2*eta2*eta2/24;
    tantheta2 = 2.0 * expTmp / ( 1.0 - expTmp*expTmp);
  }
  
  Float_t e1squ = m0_1 * m0_1 + pt1 * pt1 * (1.0 + 1.0 / tantheta1 / tantheta1);
  Float_t e2squ = m0_2 * m0_2 + pt2 * pt2 * (1.0 + 1.0 / tantheta2 / tantheta2);
  
  // fold onto 0...pi
  Float_t deltaPhi = TMath::Abs(phi1 - phi2);
  while (deltaPhi > TMath::TwoPi())
    deltaPhi -= TMath::TwoPi();
  if (deltaPhi > TMath::Pi())
    deltaPhi = TMath::TwoPi() - deltaPhi;
  
  Float_t cosDeltaPhi = 0;
  if (deltaPhi < TMath::Pi()/3)
    cosDeltaPhi = 1.0 - deltaPhi*deltaPhi/2 + deltaPhi*deltaPhi*deltaPhi*deltaPhi/24;
  else if (deltaPhi < 2*TMath::Pi()/3)
    cosDeltaPhi = -(deltaPhi - TMath::Pi()/2) + 1.0/6 * TMath::Power((deltaPhi - TMath::Pi()/2), 3);
  else
    cosDeltaPhi = -1.0 + 1.0/2.0*(deltaPhi - TMath::Pi())*(deltaPhi - TMath::Pi()) - 1.0/24.0 * TMath::Power(deltaPhi - TMath::Pi(), 4);
  
  Float_t mass2 = m0_1 * m0_1 + m0_2 * m0_2 + 2 * ( TMath::Sqrt(e1squ * e2squ) - ( pt1 * pt2 * ( cosDeltaPhi + 1.0 / tantheta1 / tantheta2 ) ) );
  
//   Printf(Form("%f %f %f %f %f %f %f %f %f", pt1, eta1, phi1, pt2, eta2, phi2, m0_1, m0_2, mass2));
  
  return mass2;
}
//--------------------------------------------------------------------------------
Float_t  AliTwoParticlePIDCorr::GetDPhiStar(Float_t phi1, Float_t pt1, Float_t charge1, Float_t phi2, Float_t pt2, Float_t charge2, Float_t radius, Float_t bSign)
{ 
  //
  // calculates dphistar
  //
  
  Float_t dphistar = phi1 - phi2 - charge1 * bSign * TMath::ASin(0.075 * radius / pt1) + charge2 * bSign * TMath::ASin(0.075 * radius / pt2);
  
  static const Double_t kPi = TMath::Pi();
  
  // circularity
//   if (dphistar > 2 * kPi)
//     dphistar -= 2 * kPi;
//   if (dphistar < -2 * kPi)
//     dphistar += 2 * kPi;
  
  if (dphistar > kPi)
    dphistar = kPi * 2 - dphistar;
  if (dphistar < -kPi)
    dphistar = -kPi * 2 - dphistar;
  if (dphistar > kPi) // might look funny but is needed
    dphistar = kPi * 2 - dphistar;
  
  return dphistar;
}
//_________________________________________________________________________
void AliTwoParticlePIDCorr ::DefineEventPool()
{
const Int_t MaxNofEvents=1000;
const Int_t MaxNofTracks=50000;
const Int_t NofCentBins=11;
//Double_t CentralityBins[]={0, 20, 40, 60,80,100.1};
Double_t CentralityBins[NofCentBins+1]={ 0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100.1};
const Int_t NofVrtxBins=10;
Double_t ZvrtxBins[]={-10,-8,-6,-4,-2,0,2,4,6,8,10};


fPoolMgr = new AliEventPoolManager(MaxNofEvents,MaxNofTracks,NofCentBins,CentralityBins,NofVrtxBins,ZvrtxBins);
fPoolMgr->SetTargetValues(MaxNofTracks, 0.1, 5);

//if(!fPoolMgr) return kFALSE;
//return kTRUE;

}
//------------------------------------------------------------------------

//________________________________________________________________________
void AliTwoParticlePIDCorr::Terminate(Option_t *) 
{
  // Draw result to screen, or perform fitting, normalizations
  // Called once at the end of the query
  fOutput = dynamic_cast<TList*> (GetOutputData(1));
  if(!fOutput) { Printf("ERROR: could not retrieve TList fOutput"); return; }
  
  
}
//------------------------------------------------------------------ 
 
