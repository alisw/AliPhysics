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
  fEta(12),
  fY(12),
  fPtMin(0.),
  fPtMax(20.),
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
  fEta(12),
  fY(12),
  fPtMin(0.),
  fPtMax(20.),
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

  TString pidNames[9] = { "Pion", "Kaon", "Proton", "K0Short", "Lambda", "Xi", "Omega", "Phi", "KStar" };

  // ### Create histograms
  fHistEvt = new TH1I("fHistEvt","fHistEvt",5,0,5) ;
  fHistEvt->GetXaxis()->SetBinLabel(1,"All events");
  fHistEvt->GetXaxis()->SetBinLabel(2,"Event with Zvtx");
  fHistEvt->GetXaxis()->SetBinLabel(3,"All Particle");
  fHistEvt->GetXaxis()->SetBinLabel(4,"Only Physical Particle");
  fHistEvt->GetXaxis()->SetBinLabel(5,"After Kinematic Cuts");
  fHistEvt->Sumw2();
  fListOfObjects->Add(fHistEvt);

  fHistPart = new TH1I("fHistPart","Particle Counts",4,0,4);
  fHistPart->GetXaxis()->SetBinLabel(1,"All Particle");
  fHistPart->GetXaxis()->SetBinLabel(2,"Pos Charge");
  fHistPart->GetXaxis()->SetBinLabel(3,"Neg Charge");
  fHistPart->GetXaxis()->SetBinLabel(4,"Neutral Charge");
  fHistPart->Sumw2();
  fListOfObjects->Add(fHistPart);


  InitHisto<TH1F>("fHistZvtx","ZVtx distribution", 30, -30., 30., "ZVtx", "Entries");
  InitHisto<TH1F>("fHistMultPrimary","Multiplicity Primary", 100, 0., 2000., "N_{prim.}", "Entries");
  InitHisto<TH1F>("fHistMultPhysPrimary","Multiplicity Physical Primary", 200, 0, 400, "N_{phys.prim.}", "Entries");
  InitHisto<TH1F>("fHistMultNotPhysPrimary","Multiplicity Not Physical Primary", 200, 0, 400., "N_{not phys.prim.}", "Entries");
  InitHisto<TH1F>("fHistMultPhysPrimaryCharged","Multiplicity Charged Physical Primary", 150, 0, 150, "N_{charged phys.prim.}", "Entries");
  InitHisto<TH1F>("fHistEta","Eta Distr.", 200, -1., 1., "#eta", "N_{part}");
  InitHisto<TH1F>("fHistY", "Y Distr.", 200, -1., 1., "#it{y}", "N_{part}");
  
  for(Int_t i=0; i<9; i++)
    InitHisto<TH1F>(Form("fHistPt_%s",pidNames[i].Data()), "Generated #it{p}_{T} distribution",400,0.,20., "#it{p}_{T} (GeV/#it{c})", "Entries");


  Int_t    fBins[2] = {  9  , 400  };
  Double_t  fMin[2] = {  0. ,   0. };
  Double_t  fMax[2] = {  9. ,  20. };
  
  THnSparseD *fHistQAPart   = new THnSparseD("fHistNParticle",  "LF particles",  2, fBins, fMin, fMax); 
  fHistQAPart->Sumw2();
  fHistQAPart->GetAxis(0)->SetBinLabel( 1, "Pions"	);
  fHistQAPart->GetAxis(0)->SetBinLabel( 2, "Kaons"	);
  fHistQAPart->GetAxis(0)->SetBinLabel( 3, "Protons"	);
  fHistQAPart->GetAxis(0)->SetBinLabel( 4, "K0Short"	);
  fHistQAPart->GetAxis(0)->SetBinLabel( 5, "Lambda"	);
  fHistQAPart->GetAxis(0)->SetBinLabel( 6, "Xi"		);
  fHistQAPart->GetAxis(0)->SetBinLabel( 7, "Omega"	);
  fHistQAPart->GetAxis(0)->SetBinLabel( 8 ,"Phi"	);
  fHistQAPart->GetAxis(0)->SetBinLabel( 9, "Kstar"	);
  fHistQAPart->GetAxis(1)->SetTitle("pt");
  fListOfObjects->Add(fHistQAPart);  
  

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
  
  if( isSelected )
    FillHisto("fHistEvt",0.5);

  AliMCEvent *event = 0x0;
  event = dynamic_cast<AliMCEvent*>(obj);
  if( !event ) 
    return kFALSE;

  const AliVVertex *vtxMC = event->GetPrimaryVertex();
    
  if( TMath::Abs(vtxMC->GetZ()) > 10. ) 
    return kFALSE;
  
  if( isSelected ) 
    FillHisto("fHistEvt",1.5);

  return isSelected;
}

//______________________________________________________________________________

Bool_t AliAnalysisTaskSimSpectraLF::IsMCParticleGenerated(TObject* obj, AliStack *fstack){

  if ( !obj ) return kFALSE;
  
  if ( !obj->InheritsFrom("AliVParticle") ) { printf("Object must derived from AliVParticle ! \n"); return kFALSE; }

  AliVParticle* vpart = dynamic_cast<AliVParticle *>(obj);
  if ( !vpart ) return kFALSE;

  Bool_t isSelected = kTRUE;

  // ### not needed for resonances, but it is kept for later reference
  /*  
  if ( !fstack->IsPhysicalPrimary(vpart->GetLabel()) )
  {
    if( vpart->PdgCode() != 11 )
      isSelected = kFALSE;
  }
  */
  return isSelected;
}

//______________________________________________________________________________

Bool_t AliAnalysisTaskSimSpectraLF::IsMCParticleInKinematicRange(TObject *obj){

  if ( !obj ) return  kFALSE;
  
  if ( !obj->InheritsFrom("AliVParticle") ) { printf("Object must derived from AliVParticle ! \n"); return kFALSE; }

  AliVParticle* vpart = dynamic_cast<AliVParticle *>(obj);
  if ( !vpart ) return kFALSE;

  Bool_t isSelected = kTRUE;

  if ( vpart->Pt() < fPtMin || vpart->Pt() > fPtMax ) // set pt cut
    isSelected = kFALSE;
  if ( TMath::Abs(vpart->Eta()) > fEta ) // set pseudorapidity cut
    isSelected = kFALSE;

  return isSelected;
}

//______________________________________________________________________________

void AliAnalysisTaskSimSpectraLF::EventSel(TObject* obj){

  if( !obj ) return;
  AliMCEvent *event = dynamic_cast<AliMCEvent*>(obj);
  if( !event ) return;

  const AliVVertex *vtxMC = event->GetPrimaryVertex();
  Float_t zVtx = vtxMC->GetZ();
  FillHisto("fHistZvtx",zVtx);

  Int_t multPrimary = fStack->GetNprimary();
  FillHisto("fHistMultPrimary",multPrimary);

  Int_t multPrimaryYes=0, multPrimaryCharged=0, multPrimaryNo=0;
  
  // ### particle loop
  for(Int_t i=0; i < event->GetNumberOfTracks(); i++) {

    AliVParticle *mcPart = (AliVParticle*)event->GetTrack(i);

    Int_t pdgcode = TMath::Abs(mcPart->PdgCode());
    
    if ( TMath::Abs(mcPart->Eta()) > fEta ) continue;

    if ( event->IsPhysicalPrimary(i) )
    {
      FillHisto("fHistPart",0.5);

      if ( mcPart->Charge() > 0 ) 
	FillHisto("fHistPart",1.5);
      else if ( mcPart->Charge() < 0 )
	FillHisto("fHistPart",2.5);
      else if ( mcPart->Charge() == 0 )
	FillHisto("fHistPart",3.5);

      multPrimaryYes++;

      if(pdgcode==11||pdgcode==13||pdgcode==211||pdgcode==321||pdgcode==2212)
	multPrimaryCharged++;
    }
    else 
      multPrimaryNo++;
  }

  if ( multPrimaryYes != 0 )
     FillHisto("fHistMultPhysPrimary",multPrimaryYes);
  if ( multPrimaryCharged != 0 )
    FillHisto("fHistMultPhysPrimaryCharged",multPrimaryCharged);
  if ( multPrimaryNo != 0 )
    FillHisto("fHistMultNotPhysPrimary",multPrimaryNo);
    
}

//______________________________________________________________________________

void AliAnalysisTaskSimSpectraLF::ParticleSel(TObject* obj){

  TString pidNames[9] = { "Pion", "Kaon", "Proton", "K0Short", "Lambda", "Xi", "Omega", "Phi", "KStar" };

  if ( !obj ) return;

  AliMCEvent *event = dynamic_cast<AliMCEvent*>(obj);
  if ( !event ) return;
    
  Bool_t isPrimary[9] = { kTRUE, kTRUE, kTRUE, kTRUE, kTRUE, kTRUE, kTRUE, kFALSE, kFALSE };

  Short_t pidCodeMC = 0;
  Double_t myY = 0.;
  Double_t ipt = 0.;

  Bool_t isPhysPrim = kFALSE;
  
  // ### particle loop
  for (Int_t ipart = 0; ipart < event->GetNumberOfTracks(); ipart++) {

    AliVParticle *mcPart = (AliVParticle*)event->GetTrack(ipart);
    if ( !mcPart ) continue;

    FillHisto("fHistEvt",2.5);

    Bool_t isParticleMCSelected = IsMCParticleGenerated(mcPart, fStack);
    if ( !isParticleMCSelected ) continue;
        
    FillHisto("fHistEvt",3.5);

    Bool_t isPartKineAccept = IsMCParticleInKinematicRange(mcPart);
    if ( !isPartKineAccept ) continue;

    FillHisto("fHistEvt",4.5);

    
    Int_t pPDG = TMath::Abs(mcPart->PdgCode());
    Double_t MCPartSel[2] = { static_cast<Double_t>(pPDG), mcPart->Pt() };
    
    isPhysPrim = event->IsPhysicalPrimary(ipart);

    if ( pPDG==211 || pPDG==321 || pPDG==2212 || pPDG==310 || pPDG==3122 || pPDG==3312 || pPDG==3334 || pPDG==333 || pPDG==313 )
    {
      if( isPhysPrim )
      {
	if (pPDG==211)	MCPartSel[0]	= 0.5;	// 	pi
	if (pPDG==321)	MCPartSel[0]	= 1.5;	//	K
	if (pPDG==2212)	MCPartSel[0]	= 2.5;	//	p
	if (pPDG==310)	MCPartSel[0]	= 3.5;	// 	K0Short
	if (pPDG==3122)	MCPartSel[0]	= 4.5;	//	Lambda
	if (pPDG==3312) MCPartSel[0]	= 5.5;	//	Xi -
	if (pPDG==3334) MCPartSel[0]	= 6.5;	//	Omega -
      }
      else{
	if (pPDG==333)	MCPartSel[0]	= 7.5;	//	Phi(1020)
	if (pPDG==313)	MCPartSel[0]	= 8.5;	//	K*(892)
      }
    }

    ((THnSparseD*)fListOfObjects->FindObject(Form("fHistNParticle")))->Fill(MCPartSel);
    FillHisto("fHistEta",mcPart->Eta());

    
    pidCodeMC = GetPidCode(pPDG);
    
    Bool_t isSelectedPart = kTRUE;
    for(Int_t i=0; i<9; i++) 
      if( pidCodeMC == i ) 
	isSelectedPart = kFALSE;
    if ( isSelectedPart ) continue;

    myY = myRap(mcPart->E(),mcPart->Pz());
    ipt = mcPart->Pt();
    
    for(Int_t i=0; i<9; i++)
    {
      if( pidCodeMC == i && TMath::Abs(myY) < fY)
      {
	if( isPrimary[i] == kTRUE && isPhysPrim == kFALSE ) 
	  continue;

	FillHisto(Form("fHistPt_%s",pidNames[i].Data()),ipt);
	
	if(!i)
	  FillHisto("fHistY",myY);
	
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

Double_t AliAnalysisTaskSimSpectraLF::myRap(Double_t rE, Double_t rPz) const {

  Double_t rRap = -999.;
  if( (rE-rPz+1.e-13) != 0 && (rE+rPz) != 0 )
    rRap =  0.5*TMath::Log((rE+rPz)/(rE-rPz+1.e-13));

  return rRap;
}

//_____________________________________________________________________________

Short_t AliAnalysisTaskSimSpectraLF::GetPidCode(Int_t pdgCode) const  {

  Short_t pidCode = 999;

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
    default:
      pidCode = 999;  // something else
  };
  
  return pidCode;
}


