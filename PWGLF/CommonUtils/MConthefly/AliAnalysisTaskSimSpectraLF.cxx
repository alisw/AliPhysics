/**************************************************************************
 * Copyright(c) 1998-2017, ALICE Experiment at CERN, All rights reserved. *
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
 *  
 *  AliAnalysisTaskSimSpectraLF.cxx
 *  
 *  Minimal analysis task meant to be used for on-the-fly simulations in LEGO trains (for generating light falvor particle pt spectra)
 *
 *  Gyula BENCEDI  <Gyula.Bencedi@cern.ch>, WIGNER RCP
 * 
 */

//_____ ROOT headers
#include <TList.h>
#include <TTree.h>
#include <TMath.h>
#include <TFile.h>
#include <TTreeStream.h>
#include "TChain.h"
#include <THnSparse.h>
#include <TDatabasePDG.h>
#include "TObjArray.h"
#include <TClonesArray.h>


//_____ ALIROOT headers
#include "AliAnalysisTask.h"
#include "AliAnalysisFilter.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliESDInputHandler.h"
#include "AliESDEvent.h"
#include "AliLog.h"
#include "AliVParticle.h"
#include "AliStack.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliHeader.h"
#include "AliGenPythiaEventHeader.h"
#include "AliAODInputHandler.h"
#include "AliAODHandler.h" 
#include "AliAODMCHeader.h"
#include "AliAODEvent.h"
#include "AliAODMCParticle.h"

//_____ Additional includes
#include "AliVEvent.h"
#include "AliGenEventHeader.h"
// #include "AliAnalysisUtils.h"

//_____ AnalysisTask headers
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisTaskSimSpectraLF.h"

//_____ STL includes
#include <iostream>
using namespace std;

ClassImp( AliAnalysisTaskSimSpectraLF )

//_____________________________________________________________________________

AliAnalysisTaskSimSpectraLF::AliAnalysisTaskSimSpectraLF():
  AliAnalysisTaskSE(),
  fMcEvent(0x0),
  fMcHandler(0x0),
  fStack(0),
  fY(0.5),
  fHistEvt(0x0),
  fHistPart(0x0),
  fListOfObjects(0)
{
  // Default constructor (should not be used)
}

//______________________________________________________________________________

AliAnalysisTaskSimSpectraLF::AliAnalysisTaskSimSpectraLF(const char *name):
  AliAnalysisTaskSE(name),
  fMcEvent(0x0),
  fMcHandler(0x0),
  fStack(0),
  fY(0.5),
  fHistEvt(0x0),
  fHistPart(0x0),
  fListOfObjects(0)
{
  DefineInput( 0, TChain::Class());
  DefineOutput(1, TList::Class() ); // Basic output slot 
}


//_____________________________________________________________________________

AliAnalysisTaskSimSpectraLF::~AliAnalysisTaskSimSpectraLF(){
  // Destructor
  // histograms are in the output list and deleted when the output
  // list is deleted by the TSelector dtor

  if (fHistEvt) { delete fHistEvt; fHistEvt=0x0; }
  if (fHistPart) { delete fHistPart; fHistPart=0x0; }
  if (fListOfObjects) { delete fListOfObjects; fListOfObjects=0x0; }
}

//______________________________________________________________________________

void AliAnalysisTaskSimSpectraLF::UserCreateOutputObjects(){
    
  // ### Analysis output
  fListOfObjects = new TList();
  fListOfObjects->SetOwner(kTRUE);

  TString pidNames[11] = { "Pion", "Kaon", "Proton", "K0Short", "Lambda", "Xi", "Omega", "Phi", "KStar", "KStarPM", "SigmaZero" };

  // ### Create histograms
  fHistEvt = new TH1I("fHistEvt","fHistEvt",2,0,2) ;
  fHistEvt->GetXaxis()->SetBinLabel(1,"All events");
  fHistEvt->GetXaxis()->SetBinLabel(2,"All particles");
  fHistEvt->Sumw2();
  fListOfObjects->Add(fHistEvt);

  InitHisto<TH1F>("fHistMultPrimary","Multiplicity Primary", 100, 0., 2000., "N_{prim.}", "Entries");
  InitHisto<TH1F>("fHistEta","Eta Distr.", 200, -1., 1., "#eta", "N_{part}");
  InitHisto<TH1F>("fHistY", "Y Distr.", 200, -1., 1., "#it{y}", "N_{part}");
  
  for(Int_t i=0; i<11; i++)
    InitHisto<TH1F>(Form("fHistPt_%s",pidNames[i].Data()), "Generated #it{p}_{T} distribution",2000,0.,20., "#it{p}_{T} (GeV/#it{c})", "Entries");

  // ### List of outputs
  PostData(1, fListOfObjects);

}

//______________________________________________________________________________

inline void AliAnalysisTaskSimSpectraLF::FillHisto(const char* objkey, Double_t x)
{
  TH1* hTmp = static_cast<TH1*>(fListOfObjects->FindObject(objkey));
  if(!hTmp){
    AliError(Form("Cannot find histogram: %s",objkey)) ;
    return;
  }
  hTmp->Fill(x);
}

//______________________________________________________________________________

template <class T> T* AliAnalysisTaskSimSpectraLF::InitHisto(const char* hname, const char* htitle, Int_t nxbins, Double_t xmin, Double_t xmax, const char* xtitle, const char* ytitle)
{
  T* hTmp = new T(hname, htitle, nxbins, xmin, xmax);
  hTmp->GetXaxis()->SetTitle(xtitle);
  hTmp->GetYaxis()->SetTitle(ytitle);
  hTmp->SetMarkerStyle(kFullCircle);
  hTmp->Sumw2();
  fListOfObjects->Add(hTmp);

  return hTmp;
}

//______________________________________________________________________________

void AliAnalysisTaskSimSpectraLF::Init(){
  //
  fMcHandler = dynamic_cast<AliInputEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
}

//______________________________________________________________________________

void AliAnalysisTaskSimSpectraLF::UserExec(Option_t *){

  // ### Initialize
  Init();
  
  // ### MC handler
  if(fMcHandler)
    fMcEvent = fMcHandler->MCEvent();
  else { if(fDebug > 1) printf("AliAnalysisTaskSimSpectraLF::Handler() fMcHandler = NULL\n"); return; }

  // ### MC event
  if( !fMcEvent ) { if(fDebug > 1) printf("AliAnalysisTaskSimSpectraLF::UserExec() fMcEvent = NULL \n"); return; }

  fStack = ((AliMCEvent*)fMcEvent)->Stack();
  if (!fStack) {
    Printf("ERROR: Could not retrieve MC stack \n");
    cout << "Name of the file with pb :" <<  fInputHandler->GetTree()->GetCurrentFile()->GetName() << endl;
    return;
  }
    
  // ### MC event selection
  Bool_t isEventMCSelected = IsMCEventSelected(fMcEvent);
  if( !isEventMCSelected ) return;
    
  // ### Event and particle selection
  EventSel( fMcEvent );
  ParticleSel( fMcEvent );
    
  // ### Post data for all output slots
  PostData(1, fListOfObjects);

  return;
}

//______________________________________________________________________________

Bool_t AliAnalysisTaskSimSpectraLF::IsMCEventSelected(TObject* obj){

  Bool_t isSelected = kTRUE;
  
  AliMCEvent *event = 0x0;
  event = dynamic_cast<AliMCEvent*>(obj);
  if( !event ) 
    isSelected = kFALSE;

  if( isSelected ) 
    FillHisto("fHistEvt",0.5);

  return isSelected;
}

//______________________________________________________________________________

void AliAnalysisTaskSimSpectraLF::EventSel(TObject* obj){

  if( !obj ) return;
  AliMCEvent *event = dynamic_cast<AliMCEvent*>(obj);
  if( !event ) return;

  Int_t multPrimary = fStack->GetNprimary();
  FillHisto("fHistMultPrimary",multPrimary);

}

//______________________________________________________________________________

void AliAnalysisTaskSimSpectraLF::ParticleSel(TObject* obj){

  TString pidNames[11] = { "Pion", "Kaon", "Proton", "K0Short", "Lambda", "Xi", "Omega", "Phi", "KStar", "KStarPM", "SigmaZero" };

  if ( !obj ) return;

  AliMCEvent *event = dynamic_cast<AliMCEvent*>(obj);
  if ( !event ) return;
    
  Bool_t isPrimary[11] = { kTRUE, kTRUE, kTRUE, kTRUE, kTRUE, kTRUE, kTRUE, kFALSE, kFALSE, kFALSE, kFALSE };

  Int_t pidCodeMC = 0;
  Double_t ipt = 0.;

  Bool_t isPhysPrim = kFALSE;
  
  // ### particle loop
  for (Int_t ipart = 0; ipart < event->GetNumberOfTracks(); ipart++) {

    TParticle* mcPart         = 0x0;
    mcPart                    = (TParticle *)event->Particle(ipart);
    if (!mcPart) continue;
    
    FillHisto("fHistEvt",1.5);

    Int_t pPDG = TMath::Abs(mcPart->GetPdgCode());
    pidCodeMC = GetPidCode(pPDG);

    Bool_t isSelectedPart = kTRUE;
    for(Int_t i=0; i<11; i++) 
      if( pidCodeMC == i ) 
	isSelectedPart = kFALSE;
    if ( isSelectedPart ) continue;

    FillHisto("fHistEta",mcPart->Eta());
    
    if (!(TMath::Abs(mcPart->Energy()-mcPart->Pz())>0.)) continue;
    Double_t myY = (mcPart->Energy()+mcPart->Pz())/(mcPart->Energy()-mcPart->Pz());
    if( myY <= 0 ) continue;

    Double_t y = 0.5*TMath::Log(myY);
    
    ipt = mcPart->Pt();

    isPhysPrim = event->IsPhysicalPrimary(ipart);

    for(Int_t i=0; i<11; i++)
    {
      if( pidCodeMC == i && TMath::Abs(y) < fY)
      {
	if( isPrimary[i] == kTRUE && isPhysPrim == kFALSE ) 
	  continue;

	if(!i) FillHisto("fHistY",y);
	FillHisto(Form("fHistPt_%s",pidNames[i].Data()),ipt);
	
      }
    }
    
  } // particle loop
}

//______________________________________________________________________________

void AliAnalysisTaskSimSpectraLF::Terminate(Option_t*){

  fListOfObjects = dynamic_cast<TList*> (GetOutputData(1));
  if (!fListOfObjects) { Printf("ERROR: Output list not available"); return; }

  return;
}

//_____________________________________________________________________________

Int_t AliAnalysisTaskSimSpectraLF::GetPidCode(Int_t pdgCode) const  {

  Int_t pidCode = 999;

  switch (TMath::Abs(pdgCode)) {
    case 211:
      pidCode = 0; // pion
    break;
    case 321:
      pidCode = 1; // kaon
    break;
    case 2212:
      pidCode = 2; // proton
    break;
    case 310:
      pidCode = 3; // K0s
    break;
    case 3122:
      pidCode = 4; // Lambda
    break;
    case 3312:
      pidCode = 5; // Xi-
    break;
    case 3334:
      pidCode = 6; // Omega-
    break;
    case 333:
      pidCode = 7; // phi(1020)
    break;
    case 313:
      pidCode = 8; // K*(892)0
    break;
    case 323:
      pidCode = 9; // K*(892) +-
    break;
    case 3212:
      pidCode = 10; // Sigma 0
    break;    
    default:
    break;
  };
  
  return pidCode;
}


