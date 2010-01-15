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
//
//
// First implementation of a class
// to reject tagged electron coming from conversion, pi0 and eta
// by calculating the e+e- invariant mass 
// of the tagged electron with other tracks
// after looser cuts for the partner.
// PostProcess should extract the background yield
// If running with MC, it can be compared to the expected background yield 
//
// Authors:
//   Raphaelle Bailhache <rbailhache@ikf.uni-frankfurt.de > <R.Bailhache@gsi.de >
//
//

#include <THnSparse.h>
#include <TParticle.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TString.h>
#include <TLegend.h>

#include <AliESDEvent.h>
#include <AliAODEvent.h>
#include <AliESDtrack.h>
#include <AliAODTrack.h>
#include "AliHFEelecbackground.h"
#include <AliMCEvent.h>
#include <AliStack.h>


ClassImp(AliHFEelecbackground)

Bool_t AliHFEelecbackground::fgUseMCPID = kFALSE;
const Double_t    AliHFEelecbackground::fgkMe = 0.00051099892;

//___________________________________________________________________________________________
AliHFEelecbackground::AliHFEelecbackground():
  fESD1(0x0)
  ,fAOD1(0x0)
  ,fMCEvent(0x0)
  ,fBz(0)
  ,fkVertex(0x0)
  ,fPtESD(0.0)
  ,fIndexTrack(0)
  ,fPdg(0)
  ,fLabMother(-1)
  ,fIsFrom(-1)
  ,fMotherGamma(-1)
  ,fMotherPi0(-1)
  ,fMotherC(-1)
  ,fMotherB(-1)
  ,fMotherEta(-1)
  ,fIsPartner(kFALSE)
  ,fIsSplittedTrack(kFALSE)
  ,fList(0x0)
  ,fListPostProcess(0x0)
{ 
  //
  // Default constructor
  //
  
}

//_______________________________________________________________________________________________
AliHFEelecbackground::AliHFEelecbackground(const AliHFEelecbackground &p):
  TObject(p)
  ,fESD1(0x0)
  ,fAOD1(0x0)
  ,fMCEvent(0x0)
  ,fBz(p.fBz)
  ,fkVertex(p.fkVertex)  
  ,fPtESD(p.fPtESD)
  ,fIndexTrack(0)
  ,fPdg(0)
  ,fLabMother(-1)
  ,fIsFrom(-1)
  ,fMotherGamma(-1)
  ,fMotherPi0(-1)
  ,fMotherC(-1)
  ,fMotherB(-1)
  ,fMotherEta(-1)
  ,fIsPartner(kFALSE)
  ,fIsSplittedTrack(kFALSE)
  ,fList(0x0)  
  ,fListPostProcess(0x0)
{ 
  //
  // Copy constructor
  //
}

//_______________________________________________________________________________________________
AliHFEelecbackground&
AliHFEelecbackground::operator=(const AliHFEelecbackground &)
{
  //
  // Assignment operator
  //

  AliInfo("Not yet implemented.");
  return *this;
}

//_______________________________________________________________________________________________
AliHFEelecbackground::~AliHFEelecbackground()
{
  //
  // Destructor
  //

  if(fList){
    fList->Clear();
    delete fList;
  }

  if(fListPostProcess){
    fListPostProcess->Clear();
    delete fListPostProcess;
  }
}
//___________________________________________________________________________________________
Bool_t AliHFEelecbackground::Load(const Char_t *filename)
{
  //
  // Generic container loader
  //

  if(!TFile::Open(filename)){
    return kFALSE;
  }
  TList *o = 0x0;
  if(!(o = (TList*)gFile->Get("Results"))){
    return kFALSE;
  }
  TList *oe = 0x0;
  if(!(oe = (TList*)dynamic_cast<TList *>(o->FindObject("HFEelecbackground")))){
    return kFALSE;
  }
  fList = (TList*)oe->Clone("HFEelecbackground");
  gFile->Close();
  return kTRUE;
}
//_______________________________________________________________________________________________
void AliHFEelecbackground::Reset()
{
  //
  // Reset variables
  //
  fPtESD = 0.0;
  fIndexTrack = -1;
  fPdg = -1;
  fLabMother = -1;
  fIsFrom = -1;
  fMotherGamma = -1;
  fMotherPi0 = -1;
  fMotherC = -1;
  fMotherB = -1;
  fMotherEta = -1;
  fIsPartner = kFALSE;
  fIsSplittedTrack = kFALSE;
 
}
//_______________________________________________________________________________________________
void AliHFEelecbackground::CreateHistograms(TList* const qaList)
{ 
  //
  // create histograms
  //
  if(!qaList) return;
  
  fList = qaList;
  fList->SetName("HFEelecbackground");
  
  // bins

  Int_t nBinsPt = 25;
  Double_t minPt = 0.001;
  Double_t maxPt = 10.0;
  
  Int_t nBinsInv = 50;
  Double_t minInv = 0.0;
  Double_t maxInv = 0.2;
  
  Int_t nBinsOp = 50;
  Double_t minOp = 0.0;
  Double_t maxOp = 2;

  Double_t *binLimLogPt = new Double_t[nBinsPt+1];
  Double_t *binLimPt    = new Double_t[nBinsPt+1];
  for(Int_t i=0; i<=nBinsPt; i++) binLimLogPt[i]=(Double_t)TMath::Log10(minPt) + (TMath::Log10(maxPt)-TMath::Log10(minPt))/nBinsPt*(Double_t)i ;
  for(Int_t i=0; i<=nBinsPt; i++) binLimPt[i]=(Double_t)TMath::Power(10,binLimLogPt[i]);

  Double_t *binLimInv    = new Double_t[nBinsInv+1];
  for(Int_t i=0; i<=nBinsInv; i++) binLimInv[i]=(Double_t)minInv  + (maxInv-minInv)  /nBinsInv*(Double_t)i ;
  
  Double_t *binLimOp    = new Double_t[nBinsOp+1];
  for(Int_t i=0; i<=nBinsOp; i++) binLimOp[i]=(Double_t)minOp  + (maxOp-minOp) /nBinsOp*(Double_t)i ;
  
  
  const Int_t nvarO = 3; // ptrectaggede, ptrecmother, openingangle or invmass

  Int_t iBinOInv[nvarO];
  iBinOInv[0]=nBinsPt;
  iBinOInv[1]=nBinsPt;
  iBinOInv[2]=nBinsInv;
  
  Int_t iBinOOp[nvarO];
  iBinOOp[0]=nBinsPt;
  iBinOOp[1]=nBinsPt;
  iBinOOp[2]=nBinsOp;
  
  //
  //
  //
  
  THnSparseF *openinganglepp = new THnSparseF("openinganglepp","",nvarO,iBinOOp);
  openinganglepp->SetBinEdges(0,binLimPt);
  openinganglepp->SetBinEdges(1,binLimPt);
  openinganglepp->SetBinEdges(2,binLimOp);
  openinganglepp->Sumw2();
  
  THnSparseF *openinganglenn = new THnSparseF("openinganglenn","",nvarO,iBinOOp);
  openinganglenn->SetBinEdges(0,binLimPt);
  openinganglenn->SetBinEdges(1,binLimPt);
  openinganglenn->SetBinEdges(2,binLimOp);
  openinganglenn->Sumw2();
  
  THnSparseF *openingangless = new THnSparseF("openingangless","",nvarO,iBinOOp);
  openingangless->SetBinEdges(0,binLimPt);
  openingangless->SetBinEdges(1,binLimPt);
  openingangless->SetBinEdges(2,binLimOp);
  openingangless->Sumw2();
  
  THnSparseF *openingangler = new THnSparseF("openingangler","",nvarO,iBinOOp);
  openingangler->SetBinEdges(0,binLimPt);
  openingangler->SetBinEdges(1,binLimPt);
  openingangler->SetBinEdges(2,binLimOp);
  openingangler->Sumw2();
  
  THnSparseF *openingangleos = new THnSparseF("openingangleos","",nvarO,iBinOOp);
  openingangleos->SetBinEdges(0,binLimPt);
  openingangleos->SetBinEdges(1,binLimPt);
  openingangleos->SetBinEdges(2,binLimOp);
  openingangleos->Sumw2();

  THnSparseF *openinganglepi0=0x0;
  THnSparseF *openingangleeta=0x0;
  THnSparseF *openinganglegamma=0x0;
  THnSparseF *openingangleC=0x0;
  THnSparseF *openingangleB=0x0;
  THnSparseF *openingangleSplittedTrackss=0x0;
  THnSparseF *openingangleSplittedTrackos=0x0;

  if(HasMCData()) {

    openinganglepi0 = new THnSparseF("openinganglepi0","",nvarO,iBinOOp);
    openinganglepi0->SetBinEdges(0,binLimPt);
    openinganglepi0->SetBinEdges(1,binLimPt);
    openinganglepi0->SetBinEdges(2,binLimOp);
    openinganglepi0->Sumw2();
    
    openingangleeta = new THnSparseF("openingangleeta","",nvarO,iBinOOp);
    openingangleeta->SetBinEdges(0,binLimPt);
    openingangleeta->SetBinEdges(1,binLimPt);
    openingangleeta->SetBinEdges(2,binLimOp);
    openingangleeta->Sumw2();    

    openinganglegamma = new THnSparseF("openinganglegamma","",nvarO,iBinOOp);
    openinganglegamma->SetBinEdges(0,binLimPt);
    openinganglegamma->SetBinEdges(1,binLimPt);
    openinganglegamma->SetBinEdges(2,binLimOp);
    openinganglegamma->Sumw2();

    openingangleC = new THnSparseF("openingangleC","",nvarO,iBinOOp);
    openingangleC->SetBinEdges(0,binLimPt);
    openingangleC->SetBinEdges(1,binLimPt);
    openingangleC->SetBinEdges(2,binLimOp);
    openingangleC->Sumw2();

    openingangleB = new THnSparseF("openingangleB","",nvarO,iBinOOp);
    openingangleB->SetBinEdges(0,binLimPt);
    openingangleB->SetBinEdges(1,binLimPt);
    openingangleB->SetBinEdges(2,binLimOp);
    openingangleB->Sumw2();

    openingangleSplittedTrackss = new THnSparseF("openingangleSplittedTrackss","",nvarO,iBinOOp);
    openingangleSplittedTrackss->SetBinEdges(0,binLimPt);
    openingangleSplittedTrackss->SetBinEdges(1,binLimPt);
    openingangleSplittedTrackss->SetBinEdges(2,binLimOp);
    openingangleSplittedTrackss->Sumw2();

    openingangleSplittedTrackos = new THnSparseF("openingangleSplittedTrackos","",nvarO,iBinOOp);
    openingangleSplittedTrackos->SetBinEdges(0,binLimPt);
    openingangleSplittedTrackos->SetBinEdges(1,binLimPt);
    openingangleSplittedTrackos->SetBinEdges(2,binLimOp);
    openingangleSplittedTrackos->Sumw2();

  }
  
  //
  
  THnSparseF *invmasspp = new THnSparseF("invmasspp","",nvarO,iBinOInv);
  invmasspp->SetBinEdges(0,binLimPt);
  invmasspp->SetBinEdges(1,binLimPt);
  invmasspp->SetBinEdges(2,binLimInv);
  invmasspp->Sumw2();
  
  THnSparseF *invmassnn = new THnSparseF("invmassnn","",nvarO,iBinOInv);
  invmassnn->SetBinEdges(0,binLimPt);
  invmassnn->SetBinEdges(1,binLimPt);
  invmassnn->SetBinEdges(2,binLimInv);
  invmassnn->Sumw2();
  
  THnSparseF *invmassss = new THnSparseF("invmassss","",nvarO,iBinOInv);
  invmassss->SetBinEdges(0,binLimPt);
  invmassss->SetBinEdges(1,binLimPt);
  invmassss->SetBinEdges(2,binLimInv);
  invmassss->Sumw2();
  
  THnSparseF *invmassr = new THnSparseF("invmassr","",nvarO,iBinOInv);
  invmassr->SetBinEdges(0,binLimPt);
  invmassr->SetBinEdges(1,binLimPt);
  invmassr->SetBinEdges(2,binLimInv);
  invmassr->Sumw2();
  
  THnSparseF *invmassos = new THnSparseF("invmassos","",nvarO,iBinOInv);
  invmassos->SetBinEdges(0,binLimPt);
  invmassos->SetBinEdges(1,binLimPt);
  invmassos->SetBinEdges(2,binLimInv);
  invmassos->Sumw2();

  THnSparseF *invmasspi0=0x0;
  THnSparseF *invmasseta=0x0;
  THnSparseF *invmassgamma=0x0;
  THnSparseF *invmassC=0x0;
  THnSparseF *invmassB=0x0;
  THnSparseF *invmassSplittedTrackss=0x0;
  THnSparseF *invmassSplittedTrackos=0x0;
  
  if(HasMCData()) {
    
    invmasspi0 = new THnSparseF("invmasspi0","",nvarO,iBinOInv);
    invmasspi0->SetBinEdges(0,binLimPt);
    invmasspi0->SetBinEdges(1,binLimPt);
    invmasspi0->SetBinEdges(2,binLimInv);
    invmasspi0->Sumw2();
    
    invmasseta = new THnSparseF("invmasseta","",nvarO,iBinOInv);
    invmasseta->SetBinEdges(0,binLimPt);
    invmasseta->SetBinEdges(1,binLimPt);
    invmasseta->SetBinEdges(2,binLimInv);
    invmasseta->Sumw2();
    
    invmassgamma = new THnSparseF("invmassgamma","",nvarO,iBinOInv);
    invmassgamma->SetBinEdges(0,binLimPt);
    invmassgamma->SetBinEdges(1,binLimPt);
    invmassgamma->SetBinEdges(2,binLimInv);
    invmassgamma->Sumw2();

    invmassC = new THnSparseF("invmassC","",nvarO,iBinOInv);
    invmassC->SetBinEdges(0,binLimPt);
    invmassC->SetBinEdges(1,binLimPt);
    invmassC->SetBinEdges(2,binLimInv);
    invmassC->Sumw2();

    invmassB = new THnSparseF("invmassB","",nvarO,iBinOInv);
    invmassB->SetBinEdges(0,binLimPt);
    invmassB->SetBinEdges(1,binLimPt);
    invmassB->SetBinEdges(2,binLimInv);
    invmassB->Sumw2();

    invmassSplittedTrackss = new THnSparseF("invmassSplittedTrackss","",nvarO,iBinOInv);
    invmassSplittedTrackss->SetBinEdges(0,binLimPt);
    invmassSplittedTrackss->SetBinEdges(1,binLimPt);
    invmassSplittedTrackss->SetBinEdges(2,binLimInv);
    invmassSplittedTrackss->Sumw2();

    invmassSplittedTrackos = new THnSparseF("invmassSplittedTrackos","",nvarO,iBinOInv);
    invmassSplittedTrackos->SetBinEdges(0,binLimPt);
    invmassSplittedTrackos->SetBinEdges(1,binLimPt);
    invmassSplittedTrackos->SetBinEdges(2,binLimInv);
    invmassSplittedTrackos->Sumw2();
  
  }
  
  //
  //
  //

  fList->AddAt(openinganglepp,kPp);
  fList->AddAt(openinganglenn,kNn);
  fList->AddAt(openingangless,kSs);
  fList->AddAt(openingangler,kR);
  fList->AddAt(openingangleos,kOs);
  
  fList->AddAt(invmasspp,kNSignComb+kPp);
  fList->AddAt(invmassnn,kNSignComb+kNn);
  fList->AddAt(invmassss,kNSignComb+kSs);
  fList->AddAt(invmassr,kNSignComb+kR);
  fList->AddAt(invmassos,kNSignComb+kOs);

  if(HasMCData()) {
    fList->AddAt(openinganglegamma,2*kNSignComb+kElectronFromGamma);
    fList->AddAt(openinganglepi0,2*kNSignComb+kElectronFromPi0);
    fList->AddAt(openingangleC,2*kNSignComb+kElectronFromC);
    fList->AddAt(openingangleB,2*kNSignComb+kElectronFromB);
    fList->AddAt(openingangleeta,2*kNSignComb+kElectronFromEta);
    fList->AddAt(openingangleSplittedTrackss,2*kNSignComb+kSplittedTrackss);
    fList->AddAt(openingangleSplittedTrackos,2*kNSignComb+kSplittedTrackos);

    fList->AddAt(invmassgamma,2*kNSignComb+kNMCInfo+kElectronFromGamma);
    fList->AddAt(invmasspi0,2*kNSignComb+kNMCInfo+kElectronFromPi0);
    fList->AddAt(invmassC,2*kNSignComb+kNMCInfo+kElectronFromC);
    fList->AddAt(invmassB,2*kNSignComb+kNMCInfo+kElectronFromB);
    fList->AddAt(invmasseta,2*kNSignComb+kNMCInfo+kElectronFromEta);
    fList->AddAt(invmassSplittedTrackss,2*kNSignComb+kNMCInfo+kSplittedTrackss);
    fList->AddAt(invmassSplittedTrackos,2*kNSignComb+kNMCInfo+kSplittedTrackos);
    
  }

}
//_______________________________________________________________________________________________
void AliHFEelecbackground::PairAnalysis(AliESDtrack* const track, AliESDtrack* const trackPart)
{
  //
  // calculate (tagged e-partner) dca, opening angle, invariant mass 
  //
  
  ////////////////////////
  // Partner track cut
  ////////////////////////
  if(!SingleTrackCut(trackPart)) return;

  ////////////////////////
  // Take label
  ////////////////////////
  Int_t indexTrack = TMath::Abs(track->GetLabel());
  Int_t indexTrackPart = TMath::Abs(trackPart->GetLabel());

  /////////////////////////
  // If MC data
  ////////////////////////
  
  if(HasMCData()) {
    
    // Take info track if not already done 
    if(indexTrack!= fIndexTrack) {
      
      fIsFrom = -1;
          
      fPdg = GetPdg(indexTrack); 
      fLabMother = GetLabMother(indexTrack);
      
      fMotherGamma = IsMotherGamma(indexTrack);
      if((fMotherGamma != -1) && ((TMath::Abs(fPdg)) == 11)) fIsFrom = kElectronFromGamma;
      fMotherPi0 = IsMotherPi0(indexTrack);
      if((fMotherPi0 != -1) && ((TMath::Abs(fPdg)) == 11)) fIsFrom = kElectronFromPi0;
      fMotherC = IsMotherC(indexTrack);
      if((fMotherC != -1) && ((TMath::Abs(fPdg)) == 11)) fIsFrom = kElectronFromC;
      fMotherB = IsMotherB(indexTrack);
      if((fMotherB != -1) && ((TMath::Abs(fPdg)) == 11)) fIsFrom = kElectronFromB;
      fMotherEta = IsMotherEta(indexTrack);
      if((fMotherEta != -1) && ((TMath::Abs(fPdg)) == 11)) fIsFrom = kElectronFromEta;
      
      fIndexTrack = indexTrack;
      
    }

    // MC PID for tagged
    if(fgUseMCPID) {
      if(TMath::Abs(fPdg) != 11) return;
    }
    
    // Look at trackPart
    fIsPartner = kFALSE;
    Int_t pdgPart = GetPdg(indexTrackPart);
    if(TMath::Abs(pdgPart) == 11) {
      Int_t labMotherPart = GetLabMother(indexTrackPart);
      if((labMotherPart == fLabMother) && (indexTrack != indexTrackPart) && (TMath::Abs(fPdg) == 11) && (fPdg*pdgPart < 0) && (fLabMother >=0) && (fLabMother < (((AliStack *)fMCEvent->Stack())->GetNtrack()))) fIsPartner = kTRUE;
      // special case of c and b
      Int_t motherCPart = IsMotherC(indexTrackPart);
      if((motherCPart != -1) && (fIsFrom == kElectronFromC) && (fPdg*pdgPart < 0)){
	fIsPartner = kTRUE;	
      }
      Int_t motherBPart = IsMotherB(indexTrackPart);
      if((motherBPart != -1) && (fIsFrom == kElectronFromB) && (fPdg*pdgPart < 0)){
	fIsPartner = kTRUE;	
      }
    }

    // Look at splitted tracks
    fIsSplittedTrack = kFALSE;
    if(indexTrackPart == fIndexTrack) fIsSplittedTrack = kTRUE;
    
  }
 
  //////////////////////
  // Sign
  /////////////////////
  Int_t sign = -1;
  if((track->Charge() > 0.0) && (trackPart->Charge() > 0.0)) sign = kPp; 
  if((track->Charge() < 0.0) && (trackPart->Charge() < 0.0)) sign = kNn; 
  if(((track->Charge() > 0.0) && (trackPart->Charge() < 0.0)) || ((track->Charge() < 0.0) && (trackPart->Charge() > 0.0))) sign = kOs; 
       
  //////////////////////
  // DCA
  /////////////////////
  
  Double_t xthis,xp;
  Double_t dca = track->GetDCA(trackPart,fBz,xthis,xp);
  if(TMath::Abs(dca) > 3.0) return;
  
  /////////////////////////////
  // Propagate
  ////////////////////////////
      
  Double_t norradius = TMath::Sqrt(fkVertex->GetX()*fkVertex->GetX()+fkVertex->GetY()*fkVertex->GetY());
  
  AliESDtrack *trackCopy = new AliESDtrack(*track);
  AliESDtrack *trackPartCopy = new AliESDtrack(*trackPart);
  Bool_t propagateok = kTRUE;
  if((!(trackPartCopy->PropagateTo(norradius,fBz))) || (!(trackCopy->PropagateTo(norradius,fBz)))) propagateok = kFALSE;
  if(!propagateok) {
    if(trackCopy) delete trackCopy;
    if(trackPartCopy) delete trackPartCopy;
    return;
  }  
  
  ///////////////////////////////////
  // Calcul mother variables
  ///////////////////////////////////
  Double_t results[5];
  Double_t resultsr[5];
  
  CalculateMotherVariable(trackCopy,trackPartCopy,&results[0]);
  CalculateMotherVariableR(trackCopy,trackPartCopy,&resultsr[0]);
  
  /////////////////////////////////////
  // Fill
  /////////////////////////////////////
   
  FillOutput(results, resultsr, sign);
  
  if(trackCopy) delete trackCopy;
  if(trackPartCopy) delete trackPartCopy;
  
}
//_____________________________________________________________________________________
void AliHFEelecbackground::CalculateMotherVariable(AliESDtrack* const track, AliESDtrack* const trackpart, Double_t *results)
{
  //
  // variables mother and take the pt of the track
  //
  // results contain: ptmother, invmass, etamother, phimother, openingangle
  //
  
  TVector3 v3Dtagged;
  TVector3 v3Dpart;
  
  Double_t *pxyz = new Double_t[3];
  track->PxPyPz(pxyz);
  v3Dtagged.SetXYZ(pxyz[0],pxyz[1],pxyz[2]);
  fPtESD = TMath::Sqrt(pxyz[0]*pxyz[0]+pxyz[1]*pxyz[1]); 

  Double_t *pxyzpart = new Double_t[3];
  trackpart->PxPyPz(pxyzpart);
  v3Dpart.SetXYZ(pxyzpart[0],pxyzpart[1],pxyzpart[2]);
 

  
  TVector3 motherrec = v3Dtagged + v3Dpart;
  
  Double_t etaESDmother = motherrec.Eta();
  Double_t ptESDmother  = motherrec.Pt();
  Double_t phiESDmother = motherrec.Phi();
  if(phiESDmother > TMath::Pi()) phiESDmother = phiESDmother - (2*TMath::Pi());
  
  
  // openinganglepropagated
  Double_t openingangle = v3Dtagged.Angle(v3Dpart);
  
  // invmass
  Double_t pESD      = TMath::Sqrt(pxyz[0]*pxyz[0]+pxyz[1]*pxyz[1]+pxyz[2]*pxyz[2]);
  Double_t pESDpart  = TMath::Sqrt(pxyzpart[0]*pxyzpart[0]+pxyzpart[1]*pxyzpart[1]+pxyzpart[2]*pxyzpart[2]);
  
  // e propagate
  Double_t eESD     = TMath::Sqrt(pESD*pESD+fgkMe*fgkMe);
  Double_t eESDpart = TMath::Sqrt(pESDpart*pESDpart+fgkMe*fgkMe);
	    
  Double_t invmass = TMath::Sqrt((eESD+eESDpart)*(eESD+eESDpart)-(motherrec.Px()*motherrec.Px()+motherrec.Py()*motherrec.Py()+motherrec.Pz()*motherrec.Pz()));

  if(!results) {
    results = new Double_t[5];
  }
  results[0] = ptESDmother;
  results[1] = etaESDmother;
  results[2] = phiESDmother;
  results[3] = invmass;
  results[4] = openingangle;
  
}
//_____________________________________________________________________________________
void AliHFEelecbackground::CalculateMotherVariableR(AliESDtrack* const track, AliESDtrack* const trackpart, Double_t *results)
{
  //
  // variables mother
  //
  // results contain: ptmother, invmass, etamother, phimother, openingangle
  //
  
  TVector3 v3Dtagged;
  TVector3 v3Dpart;
  
  Double_t *pxyz = new Double_t[3];
  track->PxPyPz(pxyz);
  v3Dtagged.SetXYZ(pxyz[0],pxyz[1],pxyz[2]);
  Double_t *pxyzpart = new Double_t[3];
  trackpart->PxPyPz(pxyzpart);
  v3Dpart.SetXYZ(pxyzpart[0],pxyzpart[1],pxyzpart[2]);

  // rotate the partner
  v3Dpart.RotateZ(TMath::Pi());
  v3Dpart.GetXYZ(pxyzpart);

  
  TVector3 motherrec = v3Dtagged + v3Dpart;
  
  Double_t etaESDmother = motherrec.Eta();
  Double_t ptESDmother  = motherrec.Pt();
  Double_t phiESDmother = motherrec.Phi();
  if(phiESDmother > TMath::Pi()) phiESDmother = phiESDmother - (2*TMath::Pi());
  
  
  // openinganglepropagated
  Double_t openingangle = v3Dtagged.Angle(v3Dpart);
  //printf("Openingangle %f\n",openingangle);
  
  // invmass
  Double_t pESD      = TMath::Sqrt(pxyz[0]*pxyz[0]+pxyz[1]*pxyz[1]+pxyz[2]*pxyz[2]);
  Double_t pESDpart  = TMath::Sqrt(pxyzpart[0]*pxyzpart[0]+pxyzpart[1]*pxyzpart[1]+pxyzpart[2]*pxyzpart[2]);
  // e propagate
  Double_t eESD     = TMath::Sqrt(pESD*pESD+fgkMe*fgkMe);
  Double_t eESDpart = TMath::Sqrt(pESDpart*pESDpart+fgkMe*fgkMe);
  
  Double_t invmass = TMath::Sqrt((eESD+eESDpart)*(eESD+eESDpart)-(motherrec.Px()*motherrec.Px()+motherrec.Py()*motherrec.Py()+motherrec.Pz()*motherrec.Pz()));

  if(!results) {
    results = new Double_t[5];
  }
  results[0] = ptESDmother;
  results[1] = etaESDmother;
  results[2] = phiESDmother;
  results[3] = invmass;
  results[4] = openingangle;

  
}
//_________________________________________________________________________________
void AliHFEelecbackground::FillOutput(Double_t *results, Double_t *resultsr, Int_t sign) 
{
  //
  // Fill the invariant mass and opening angle distributions 
  //
  
  Double_t co[3];
  co[0] = fPtESD;
   
  switch(sign){
	
      case kPp:
	co[1] = results[0];
	co[2] = TMath::Abs(results[4]);
	(dynamic_cast<THnSparseF *>(fList->At(kPp)))->Fill(co);
	(dynamic_cast<THnSparseF *>(fList->At(kSs)))->Fill(co);
	if(HasMCData()){
	  if(fIsSplittedTrack) (dynamic_cast<THnSparseF *>(fList->At(2*kNSignComb+kSplittedTrackss)))->Fill(co);
	}


	co[2] = results[3];	
	if(TMath::Abs(results[4]) < 0.8){
	  (dynamic_cast<THnSparseF *>(fList->At(kPp+kNSignComb)))->Fill(co);
	  (dynamic_cast<THnSparseF *>(fList->At(kSs+kNSignComb)))->Fill(co);
	  if(HasMCData()){
	    if(fIsSplittedTrack) (dynamic_cast<THnSparseF *>(fList->At(2*kNSignComb+kSplittedTrackss+kNMCInfo)))->Fill(co);
	  }
	}

	break;
	
      case kNn:
	co[1] = results[0];
	co[2] = TMath::Abs(results[4]);	
	(dynamic_cast<THnSparseF *>(fList->At(kNn)))->Fill(co);
	(dynamic_cast<THnSparseF *>(fList->At(kSs)))->Fill(co);
	if(HasMCData()){
	  if(fIsSplittedTrack) (dynamic_cast<THnSparseF *>(fList->At(2*kNSignComb+kSplittedTrackss)))->Fill(co);
	}
	
       	co[2] = results[3];	
	if(TMath::Abs(results[4]) < 0.8){
	  (dynamic_cast<THnSparseF *>(fList->At(kNn+kNSignComb)))->Fill(co);
	  (dynamic_cast<THnSparseF *>(fList->At(kSs+kNSignComb)))->Fill(co);
	  if(HasMCData()){
	    if(fIsSplittedTrack) (dynamic_cast<THnSparseF *>(fList->At(2*kNSignComb+kSplittedTrackss+kNMCInfo)))->Fill(co);
	  }
	}
	break;
	
      case kOs:
	co[1] = results[0];
	co[2] = TMath::Abs(results[4]);	
	(dynamic_cast<THnSparseF *>(fList->At(kOs)))->Fill(co);
	if(HasMCData()) {
	  if((fIsFrom == kElectronFromPi0) && (fIsPartner)) (dynamic_cast<THnSparseF *>(fList->At(2*kNSignComb+kElectronFromPi0)))->Fill(co);
	  if((fIsFrom == kElectronFromEta) && (fIsPartner)) (dynamic_cast<THnSparseF *>(fList->At(2*kNSignComb+kElectronFromEta)))->Fill(co);
	  if((fIsFrom == kElectronFromGamma) && (fIsPartner)) (dynamic_cast<THnSparseF *>(fList->At(2*kNSignComb+kElectronFromGamma)))->Fill(co);
	  if((fIsFrom == kElectronFromC) && (fIsPartner)) (dynamic_cast<THnSparseF *>(fList->At(2*kNSignComb+kElectronFromC)))->Fill(co);
	  if((fIsFrom == kElectronFromB) && (fIsPartner)) (dynamic_cast<THnSparseF *>(fList->At(2*kNSignComb+kElectronFromB)))->Fill(co);
	  if(fIsSplittedTrack) (dynamic_cast<THnSparseF *>(fList->At(2*kNSignComb+kSplittedTrackos)))->Fill(co);
	}	
	
	co[2] = results[3];	
	if(TMath::Abs(results[4]) < 0.8){
	  (dynamic_cast<THnSparseF *>(fList->At(kOs+kNSignComb)))->Fill(co);
	  if(HasMCData()) {
	    if((fIsFrom == kElectronFromPi0) && (fIsPartner)) (dynamic_cast<THnSparseF *>(fList->At(2*kNSignComb+kElectronFromPi0+kNMCInfo)))->Fill(co);
	    if((fIsFrom == kElectronFromEta) && (fIsPartner)) (dynamic_cast<THnSparseF *>(fList->At(2*kNSignComb+kElectronFromEta+kNMCInfo)))->Fill(co);
	    if((fIsFrom == kElectronFromGamma) && (fIsPartner)) (dynamic_cast<THnSparseF *>(fList->At(2*kNSignComb+kElectronFromGamma+kNMCInfo)))->Fill(co);
	    if((fIsFrom == kElectronFromC) && (fIsPartner)) (dynamic_cast<THnSparseF *>(fList->At(2*kNSignComb+kElectronFromC+kNMCInfo)))->Fill(co);
	    if((fIsFrom == kElectronFromB) && (fIsPartner)) (dynamic_cast<THnSparseF *>(fList->At(2*kNSignComb+kElectronFromB+kNMCInfo)))->Fill(co);
	    if(fIsSplittedTrack) (dynamic_cast<THnSparseF *>(fList->At(2*kNSignComb+kSplittedTrackos+kNMCInfo)))->Fill(co);
	  }	
	}

	// rotated
	co[1] = resultsr[0];
	co[2] = TMath::Abs(resultsr[4]);

	(dynamic_cast<THnSparseF *>(fList->At(kR)))->Fill(co);

	co[2] = resultsr[3];	
	if(TMath::Abs(resultsr[4]) < 0.8){
	  (dynamic_cast<THnSparseF *>(fList->At(kR+kNSignComb)))->Fill(co);
	}
	break;
	
  default:
    
    break;
    
  }
}      
//_______________________________________________________________________________________________
Bool_t AliHFEelecbackground::SingleTrackCut(AliESDtrack* const track) const
{
  //
  // Return minimum quality for the partner
  //
  
  //if(track->GetKinkIndex(0)>0) return kFALSE;

  UInt_t status = track->GetStatus();
  
  if(((status & AliESDtrack::kTPCin)==0) && (status & AliESDtrack::kITSin)) {
    
    Int_t nbcl = track->GetITSclusters(0);
    if(nbcl > 1) return kTRUE;
    else return kFALSE;
  
  }
  
  if(status & AliESDtrack::kTPCin) {
  
    if(status & AliESDtrack::kTPCrefit)  return kTRUE;
    else return kFALSE;
  
  }
  
  return kFALSE;
  
}
//______________________________________________________________________________________________
void AliHFEelecbackground::SetEvent(AliESDEvent* const ESD)
{
  //
  // Set the AliESD Event, the magnetic field and the primary vertex
  //
  
  fESD1=ESD;
  fBz=fESD1->GetMagneticField();
  fkVertex = fESD1->GetPrimaryVertex();

}
//____________________________________________________________________________________________________________
Int_t AliHFEelecbackground::IsMotherGamma(Int_t tr) {

  //
  // Return the lab of gamma mother or -1 if not gamma
  //

  AliStack* stack = fMCEvent->Stack();
  if((tr < 0) || (tr >= stack->GetNtrack())) return -1;  

  // Take mother
  TParticle * particle = stack->Particle(tr);
  if(!particle) return -1;
  Int_t imother   = particle->GetFirstMother(); 
  if((imother < 0) || (imother >= stack->GetNtrack())) return -1;  
  TParticle * mother = stack->Particle(imother);
  if(!mother) return -1;

  // Check gamma    
  Int_t pdg = mother->GetPdgCode();
  if(TMath::Abs(pdg) == 22) return imother;
  if(TMath::Abs(pdg) == 11) {
    return IsMotherGamma(imother);
  }
  return -1;
 
}
//
//____________________________________________________________________________________________________________
Int_t AliHFEelecbackground::IsMotherPi0(Int_t tr) {

  //
  // Return the lab of pi0 mother or -1 if not pi0
  //

  AliStack* stack = fMCEvent->Stack();
  if((tr < 0) || (tr >= stack->GetNtrack())) return -1;  

  // Take mother
  TParticle * particle = stack->Particle(tr);
  if(!particle) return -1;
  Int_t imother   = particle->GetFirstMother(); 
  if((imother < 0) || (imother >= stack->GetNtrack())) return -1;  
  TParticle * mother = stack->Particle(imother);
  if(!mother) return -1;

  // Check gamma    
  Int_t pdg = mother->GetPdgCode();
  if(TMath::Abs(pdg) == 111) return imother;
  if(TMath::Abs(pdg) == 11) {
    return IsMotherPi0(imother);
  }
  return -1;
 
}
//____________________________________________________________________________________________________________
Int_t AliHFEelecbackground::IsMotherEta(Int_t tr) {

  //
  // Return the lab of pi0 mother or -1 if not pi0
  //

  AliStack* stack = fMCEvent->Stack();
  if((tr < 0) || (tr >= stack->GetNtrack())) return -1;  

  // Take mother
  TParticle * particle = stack->Particle(tr);
  if(!particle) return -1;
  Int_t imother   = particle->GetFirstMother(); 
  if((imother < 0) || (imother >= stack->GetNtrack())) return -1;  
  TParticle * mother = stack->Particle(imother);
  if(!mother) return -1;

  // Check gamma    
  Int_t pdg = mother->GetPdgCode();
  if(TMath::Abs(pdg) == 221) return imother;
  if(TMath::Abs(pdg) == 11) {
    return IsMotherEta(imother);
  }
  return -1;
 
}
//____________________________________________________________________________________________________________
Int_t AliHFEelecbackground::IsMotherC(Int_t tr) {

  //
  // Return the lab of signal mother or -1 if not signal
  //

  AliStack* stack = fMCEvent->Stack();
  if((tr < 0) || (tr >= stack->GetNtrack())) return -1;  

  // Take mother
  TParticle * particle = stack->Particle(tr);
  if(!particle) return -1;
  Int_t imother   = particle->GetFirstMother(); 
  if((imother < 0) || (imother >= stack->GetNtrack())) return -1;  
  TParticle * mother = stack->Particle(imother);
  if(!mother) return -1;

  // Check gamma    
  Int_t pdg = mother->GetPdgCode();
  if((TMath::Abs(pdg)==411) || (TMath::Abs(pdg)==421) || (TMath::Abs(pdg)==431) || (TMath::Abs(pdg)==4122) || (TMath::Abs(pdg)==4132) || (TMath::Abs(pdg)==4232) || (TMath::Abs(pdg)==43320)) return imother;
  if(TMath::Abs(pdg) == 11) {
    return IsMotherC(imother);
  }
  return -1;
 
}
//____________________________________________________________________________________________________________
Int_t AliHFEelecbackground::IsMotherB(Int_t tr) {

  //
  // Return the lab of signal mother or -1 if not signal
  //

  AliStack* stack = fMCEvent->Stack();
  if((tr < 0) || (tr >= stack->GetNtrack())) return -1;  

  // Take mother
  TParticle * particle = stack->Particle(tr);
  if(!particle) return -1;
  Int_t imother   = particle->GetFirstMother(); 
  if((imother < 0) || (imother >= stack->GetNtrack())) return -1;  
  TParticle * mother = stack->Particle(imother);
  if(!mother) return -1;

  // Check gamma    
  Int_t pdg = mother->GetPdgCode();
  if((TMath::Abs(pdg)==511) || (TMath::Abs(pdg)==521) || (TMath::Abs(pdg)==531) || (TMath::Abs(pdg)==5122) || (TMath::Abs(pdg)==5132) || (TMath::Abs(pdg)==5232) || (TMath::Abs(pdg)==53320)) return imother;
  if(TMath::Abs(pdg) == 11) {
    return IsMotherB(imother);
  }
  return -1;
 
}
//____________________________________________________________________________________________________________
Int_t AliHFEelecbackground::GetPdg(Int_t tr) {

  //
  // Simply pdg code
  //

  AliStack* stack = fMCEvent->Stack();
  if((tr < 0) || (tr >= stack->GetNtrack())) return -1;  

  // MC Information
  TParticle * particle = stack->Particle(tr);
  if(!particle) return -1;
  Int_t pdg = particle->GetPdgCode();

  return pdg;
 
}
//____________________________________________________________________________________________________________
Int_t AliHFEelecbackground::GetLabMother(Int_t tr) {

  //
  // Simply lab mother
  //

  AliStack* stack = fMCEvent->Stack();
  if((tr < 0) || (tr >= stack->GetNtrack())) return -1;  

  // MC Information
  TParticle * particle = stack->Particle(tr);
  if(!particle) return -1;
  Int_t imother = particle->GetFirstMother(); 

  return imother;
 
}
//_______________________________________________________________________________________________
void AliHFEelecbackground::PostProcess()
{
  //
  // Post process the histos and extract the background pt spectra
  //

  if(!fList) return;

  // invariant mass input spectra
  THnSparseF *invmassss = dynamic_cast<THnSparseF *>(fList->FindObject("invmassss"));
  THnSparseF *invmassr  = dynamic_cast<THnSparseF *>(fList->FindObject("invmassr"));
  THnSparseF *invmassos = dynamic_cast<THnSparseF *>(fList->FindObject("invmassos"));
  THnSparseF *invmassgamma = dynamic_cast<THnSparseF *>(fList->FindObject("invmassgamma"));
  THnSparseF *invmasspi0 = dynamic_cast<THnSparseF *>(fList->FindObject("invmasspi0"));
  THnSparseF *invmasseta = dynamic_cast<THnSparseF *>(fList->FindObject("invmasseta"));
  THnSparseF *invmassC = dynamic_cast<THnSparseF *>(fList->FindObject("invmassC"));
  THnSparseF *invmassB = dynamic_cast<THnSparseF *>(fList->FindObject("invmassB"));
  
  TAxis *ptaxisinvmass = invmassss->GetAxis(0);
  Int_t  nbinsptinvmass = ptaxisinvmass->GetNbins();  
  
  // outputs
  TH1D **invmassosptproj = new TH1D*[nbinsptinvmass];
  TH1D **invmassssptproj = new TH1D*[nbinsptinvmass];
  TH1D **invmassrptproj = new TH1D*[nbinsptinvmass];
  TH1D **invmassdiffptproj = new TH1D*[nbinsptinvmass];
  TH1D **invmassgammaptproj = new TH1D*[nbinsptinvmass];
  TH1D **invmasspi0ptproj = new TH1D*[nbinsptinvmass];
  TH1D **invmassetaptproj = new TH1D*[nbinsptinvmass];
  TH1D **invmassCptproj = new TH1D*[nbinsptinvmass];
  TH1D **invmassBptproj = new TH1D*[nbinsptinvmass];

  TH1D *yieldPtFound = (TH1D *) invmassss->Projection(0);
  yieldPtFound->SetName("Found yield");
  yieldPtFound->Reset();

  TH1D *yieldPtSourcesMC = 0x0;
  if(invmasspi0 && invmasseta && invmassgamma) {
    yieldPtSourcesMC = (TH1D *) invmassss->Projection(0);
    yieldPtSourcesMC->SetName("Found yield");
    yieldPtSourcesMC->Reset();
  }

  TH1D *yieldPtSignalCutMC = 0x0;
  if(invmassC && invmassB) {
    yieldPtSignalCutMC = (TH1D *) invmassss->Projection(0);
    yieldPtSignalCutMC->SetName("Found yield");
    yieldPtSignalCutMC->Reset();
  }

  // canvas
  Int_t nbrow = (Int_t) (nbinsptinvmass/5);
  TString namecanvas("Invmassnamecanvas");
  TCanvas * canvas =new TCanvas(namecanvas,namecanvas,800,800);
  canvas->Divide(5,nbrow+1);


  // loop on pt bins
  for(Int_t k=1; k <= nbinsptinvmass; k++){

    Double_t lowedge = ptaxisinvmass->GetBinLowEdge(k);
    Double_t upedge  = ptaxisinvmass->GetBinUpEdge(k);
    
    ((TAxis *)invmassss->GetAxis(0))->SetRange(k,k);
    ((TAxis *)invmassr->GetAxis(0))->SetRange(k,k);
    ((TAxis *)invmassos->GetAxis(0))->SetRange(k,k);
     
    invmassosptproj[k-1] = invmassos->Projection(2);
    invmassssptproj[k-1] = invmassss->Projection(2);
    invmassrptproj[k-1]  = invmassr->Projection(2);
    invmassgammaptproj[k-1] = 0x0;
    invmasspi0ptproj[k-1] = 0x0;
    invmassetaptproj[k-1] = 0x0;
    invmassCptproj[k-1] = 0x0;
    invmassBptproj[k-1] = 0x0;
    if(invmassgamma) invmassgammaptproj[k-1] = invmassgamma->Projection(2);
    if(invmasspi0) invmasspi0ptproj[k-1] = invmasspi0->Projection(2);
    if(invmasseta) invmassetaptproj[k-1] = invmasseta->Projection(2);
    if(invmassC) invmassCptproj[k-1] = invmassC->Projection(2);
    if(invmassB) invmassBptproj[k-1] = invmassB->Projection(2);
    
    invmassdiffptproj[k-1] = (TH1D *) invmassosptproj[k-1]->Clone();
    TString name("Invmassdiffptbin");
    name += k;
    invmassdiffptproj[k-1]->SetName(name);
    invmassdiffptproj[k-1]->Add(invmassssptproj[k-1],-1.0);

    TString namee("p_{T} tagged from ");
    namee += lowedge;
    namee += " GeV/c to ";
    namee += upedge;
    namee += " GeV/c";

    invmassosptproj[k-1]->SetTitle((const char*)namee);
    invmassssptproj[k-1]->SetTitle((const char*)namee);
    invmassrptproj[k-1]->SetTitle((const char*)namee);
    invmassdiffptproj[k-1]->SetTitle((const char*)namee);
    if(invmassgammaptproj[k-1]) invmassgammaptproj[k-1]->SetTitle((const char*)namee);
    if(invmasspi0ptproj[k-1]) invmasspi0ptproj[k-1]->SetTitle((const char*)namee);
    if(invmassetaptproj[k-1]) invmassetaptproj[k-1]->SetTitle((const char*)namee);
    if(invmassCptproj[k-1]) invmassCptproj[k-1]->SetTitle((const char*)namee);
    if(invmassBptproj[k-1]) invmassBptproj[k-1]->SetTitle((const char*)namee);
        


    invmassosptproj[k-1]->SetStats(0);
    invmassssptproj[k-1]->SetStats(0);
    invmassrptproj[k-1]->SetStats(0);
    invmassdiffptproj[k-1]->SetStats(0);
    if(invmassgammaptproj[k-1]) invmassgammaptproj[k-1]->SetStats(0);
    if(invmasspi0ptproj[k-1]) invmasspi0ptproj[k-1]->SetStats(0);
    if(invmassetaptproj[k-1]) invmassetaptproj[k-1]->SetStats(0);
    if(invmassCptproj[k-1]) invmassCptproj[k-1]->SetStats(0);
    if(invmassBptproj[k-1]) invmassBptproj[k-1]->SetStats(0);
        
    Double_t yieldf = invmassdiffptproj[k-1]->Integral();
    if(invmassetaptproj[k-1] && invmasspi0ptproj[k-1] && invmassgammaptproj[k-1] && invmassCptproj[k-1] && invmassBptproj[k-1]) {
      Double_t yieldg = invmassetaptproj[k-1]->Integral() + invmasspi0ptproj[k-1]->Integral() + invmassgammaptproj[k-1]->Integral();
      yieldPtSourcesMC->SetBinContent(k,yieldg);
      
      Double_t yieldsignal = invmassCptproj[k-1]->Integral() + invmassBptproj[k-1]->Integral();
      yieldPtSignalCutMC->SetBinContent(k,yieldsignal);
    }

    yieldPtFound->SetBinContent(k,yieldf);
    
    canvas->cd(k);
    invmassosptproj[k-1]->Draw();
    invmassssptproj[k-1]->Draw("same");
    invmassdiffptproj[k-1]->Draw("same");
    invmassrptproj[k-1]->Draw("same");
    TLegend *legiv = new TLegend(0.4,0.6,0.89,0.89);
    legiv->AddEntry(invmassosptproj[k-1],"Opposite signs","p"); 
    legiv->AddEntry(invmassssptproj[k-1],"Same signs","p"); 
    legiv->AddEntry(invmassdiffptproj[k-1],"(Opposite - Same) signs","p"); 
    legiv->AddEntry(invmassrptproj[k-1],"rotated","p"); 
    if(invmassgammaptproj[k-1]) legiv->AddEntry(invmassgammaptproj[k-1],"e^{+}e^{-} from #gamma","p"); 
    if(invmasspi0ptproj[k-1]) legiv->AddEntry(invmasspi0ptproj[k-1],"e^{+}e^{-} from #pi^{0}","p"); 
    if(invmassetaptproj[k-1]) legiv->AddEntry(invmassetaptproj[k-1],"e^{+}e^{-} from #eta","p"); 
    legiv->Draw("same");
  
  }

  yieldPtFound->SetStats(0);
  if(yieldPtSourcesMC) yieldPtSourcesMC->SetStats(0); 
  if(yieldPtSignalCutMC) yieldPtSignalCutMC->SetStats(0);  

  TCanvas * canvasfin =new TCanvas("results","results",800,800);
  canvasfin->cd(1);
  yieldPtFound->Draw();
  if(yieldPtSourcesMC && yieldPtSignalCutMC) {
    yieldPtSourcesMC->Draw("same");
    yieldPtSignalCutMC->Draw("same");
    TLegend *lega = new TLegend(0.4,0.6,0.89,0.89);
    lega->AddEntry(yieldPtFound,"Contributions found","p"); 
    lega->AddEntry(yieldPtSourcesMC,"Contributions of e^{+}e^{-} from #gamma, #pi^{0} and #eta","p"); 
    lega->AddEntry(yieldPtSignalCutMC,"Contributions of e^{+}e^{-} from C and B","p"); 
    lega->Draw("same");
  }
  

  
  if(!fListPostProcess) fListPostProcess = new TList();
  fListPostProcess->SetName("ListPostProcess");
  
  for(Int_t k=0; k < nbinsptinvmass; k++){
    fListPostProcess->AddAt(invmassosptproj[k],kOos+kNOutput*k);
    fListPostProcess->AddAt(invmassssptproj[k],kOss+kNOutput*k);
    fListPostProcess->AddAt(invmassrptproj[k],kOr+kNOutput*k);
    fListPostProcess->AddAt(invmassdiffptproj[k],kOdiff+kNOutput*k);
    if(invmassgammaptproj[k]) fListPostProcess->AddAt(invmassgammaptproj[k],kNOutput*nbinsptinvmass+kNMCInfo*k+kElectronFromGamma);
    if(invmasspi0ptproj[k]) fListPostProcess->AddAt(invmasspi0ptproj[k],kNOutput*nbinsptinvmass+kNMCInfo*k+kElectronFromPi0);
    if(invmassetaptproj[k]) fListPostProcess->AddAt(invmassetaptproj[k],kNOutput*nbinsptinvmass+kNMCInfo*k+kElectronFromEta);
    if(invmassCptproj[k]) fListPostProcess->AddAt(invmassCptproj[k],kNOutput*nbinsptinvmass+kNMCInfo*k+kElectronFromC);
    if(invmassBptproj[k]) fListPostProcess->AddAt(invmassBptproj[k],kNOutput*nbinsptinvmass+kNMCInfo*k+kElectronFromB);
  }

  fListPostProcess->AddAt(yieldPtFound,kNOutput*nbinsptinvmass+kNMCInfo*nbinsptinvmass);
  if(yieldPtSourcesMC) fListPostProcess->AddAt(yieldPtSourcesMC,kNOutput*nbinsptinvmass+kNMCInfo*nbinsptinvmass+1);
  if(yieldPtSignalCutMC) fListPostProcess->AddAt(yieldPtSignalCutMC,kNOutput*nbinsptinvmass+kNMCInfo*nbinsptinvmass+2);
  
  // delete dynamic array
  delete[] invmassosptproj;
  delete[] invmassssptproj;
  delete[] invmassrptproj;
  delete[] invmassdiffptproj;
  delete[] invmassgammaptproj;
  delete[] invmasspi0ptproj;
  delete[] invmassetaptproj;
  delete[] invmassCptproj;
  delete[] invmassBptproj;

}
//_______________________________________________________________________________________________
void AliHFEelecbackground::Plot() const
{
  //
  // Plot the output
  //
  
  if(!fList) return;
  
  // opening angle 
  THnSparseF *openinganglepp = dynamic_cast<THnSparseF *>(fList->FindObject("openinganglepp"));
  THnSparseF *openinganglenn = dynamic_cast<THnSparseF *>(fList->FindObject("openinganglenn"));
  THnSparseF *openingangless = dynamic_cast<THnSparseF *>(fList->FindObject("openingangless"));
  THnSparseF *openingangler  = dynamic_cast<THnSparseF *>(fList->FindObject("openingangler"));
  THnSparseF *openingangleos = dynamic_cast<THnSparseF *>(fList->FindObject("openingangleos"));
  
  THnSparseF *openinganglegamma = dynamic_cast<THnSparseF *>(fList->FindObject("openinganglegamma"));
  THnSparseF *openinganglepi0 = dynamic_cast<THnSparseF *>(fList->FindObject("openinganglepi0"));
  THnSparseF *openingangleC = dynamic_cast<THnSparseF *>(fList->FindObject("openingangleC"));
  THnSparseF *openingangleB = dynamic_cast<THnSparseF *>(fList->FindObject("openingangleB"));
  THnSparseF *openingangleeta = dynamic_cast<THnSparseF *>(fList->FindObject("openingangleeta"));  

  THnSparseF *openingangleSplittedTrackss = dynamic_cast<THnSparseF *>(fList->FindObject("openingangleSplittedTrackss"));  
  THnSparseF *openingangleSplittedTrackos = dynamic_cast<THnSparseF *>(fList->FindObject("openingangleSplittedTrackos"));  


  // invariant mass
  THnSparseF *invmasspp = dynamic_cast<THnSparseF *>(fList->FindObject("invmasspp"));
  THnSparseF *invmassnn = dynamic_cast<THnSparseF *>(fList->FindObject("invmassnn"));
  THnSparseF *invmassss = dynamic_cast<THnSparseF *>(fList->FindObject("invmassss"));
  THnSparseF *invmassr  = dynamic_cast<THnSparseF *>(fList->FindObject("invmassr"));
  THnSparseF *invmassos = dynamic_cast<THnSparseF *>(fList->FindObject("invmassos"));
  
  THnSparseF *invmassgamma = dynamic_cast<THnSparseF *>(fList->FindObject("invmassgamma"));
  THnSparseF *invmasspi0 = dynamic_cast<THnSparseF *>(fList->FindObject("invmasspi0"));
  THnSparseF *invmassC = dynamic_cast<THnSparseF *>(fList->FindObject("invmassC"));
  THnSparseF *invmassB = dynamic_cast<THnSparseF *>(fList->FindObject("invmassB"));
  THnSparseF *invmasseta = dynamic_cast<THnSparseF *>(fList->FindObject("invmasseta"));

  THnSparseF *invmassSplittedTrackss = dynamic_cast<THnSparseF *>(fList->FindObject("invmassSplittedTrackss"));
  THnSparseF *invmassSplittedTrackos = dynamic_cast<THnSparseF *>(fList->FindObject("invmassSplittedTrackos"));

  // Projection over all pt
  TH1D *openingangleppproj = openinganglepp->Projection(2);
  TH1D *openinganglennproj = openinganglenn->Projection(2);
  TH1D *openinganglessproj = openingangless->Projection(2);
  TH1D *openinganglerproj  = openingangler->Projection(2);
  TH1D *openingangleosproj = openingangleos->Projection(2);

  TH1D *openinganglegammaproj = 0x0;
  TH1D *openinganglepi0proj = 0x0;
  TH1D *openingangleCproj = 0x0;
  TH1D *openingangleBproj = 0x0;
  TH1D *openingangleetaproj = 0x0;
  TH1D *openingangleSplittedTrackssproj = 0x0;
  TH1D *openingangleSplittedTrackosproj = 0x0;
  if(openinganglegamma) openinganglegammaproj = openinganglegamma->Projection(2);
  if(openinganglepi0) openinganglepi0proj = openinganglepi0->Projection(2);
  if(openingangleC) openingangleCproj = openingangleC->Projection(2);
  if(openingangleB) openingangleBproj = openingangleB->Projection(2);
  if(openingangleeta) openingangleetaproj = openingangleeta->Projection(2);
  if(openingangleSplittedTrackss) openingangleSplittedTrackssproj = openingangleSplittedTrackss->Projection(2);
  if(openingangleSplittedTrackos) openingangleSplittedTrackosproj = openingangleSplittedTrackos->Projection(2);


  TH1D *invmassppproj = invmasspp->Projection(2);
  TH1D *invmassnnproj = invmassnn->Projection(2);
  TH1D *invmassssproj = invmassss->Projection(2);
  TH1D *invmassrproj  = invmassr->Projection(2);
  TH1D *invmassosproj = invmassos->Projection(2);

  TH1D *invmassgammaproj = 0x0;
  TH1D *invmasspi0proj = 0x0;
  TH1D *invmassCproj = 0x0;
  TH1D *invmassBproj = 0x0;
  TH1D *invmassetaproj = 0x0;
  TH1D *invmassSplittedTrackssproj = 0x0;
  TH1D *invmassSplittedTrackosproj = 0x0;
  if(invmassgamma) invmassgammaproj = invmassgamma->Projection(2);
  if(invmasspi0) invmasspi0proj = invmasspi0->Projection(2);
  if(invmassC) invmassCproj = invmassC->Projection(2);
  if(invmassB) invmassBproj = invmassB->Projection(2);
  if(invmasseta) invmassetaproj = invmasseta->Projection(2);
  if(invmassSplittedTrackss) invmassSplittedTrackssproj = invmassSplittedTrackss->Projection(2);
  if(invmassSplittedTrackos) invmassSplittedTrackosproj = invmassSplittedTrackos->Projection(2);
  
  openingangleppproj->SetStats(0);
  openinganglennproj->SetStats(0);
  openinganglessproj->SetStats(0);
  openinganglerproj->SetStats(0);
  openingangleosproj->SetStats(0);
  if(openinganglegammaproj) openinganglegammaproj->SetStats(0);
  if(openinganglepi0proj) openinganglepi0proj->SetStats(0);
  if(openingangleCproj) openingangleCproj->SetStats(0);
  if(openingangleBproj) openingangleBproj->SetStats(0);
  if(openingangleetaproj) openingangleetaproj->SetStats(0);
  if(openingangleSplittedTrackssproj) openingangleSplittedTrackssproj->SetStats(0);
  if(openingangleSplittedTrackosproj) openingangleSplittedTrackosproj->SetStats(0);
  
  invmassppproj->SetStats(0);
  invmassnnproj->SetStats(0);
  invmassssproj->SetStats(0);
  invmassrproj->SetStats(0);
  invmassosproj->SetStats(0);
  if(invmassgammaproj) invmassgammaproj->SetStats(0);
  if(invmasspi0proj) invmasspi0proj->SetStats(0);
  if(invmassCproj) invmassCproj->SetStats(0);
  if(invmassBproj) invmassBproj->SetStats(0);
  if(invmassetaproj) invmassetaproj->SetStats(0);
  if(invmassSplittedTrackssproj) invmassSplittedTrackssproj->SetStats(0);
  if(invmassSplittedTrackosproj) invmassSplittedTrackosproj->SetStats(0);

  openingangleppproj->SetTitle("");
  openinganglennproj->SetTitle("");
  openinganglessproj->SetTitle("");
  openinganglerproj->SetTitle("");
  openingangleosproj->SetTitle("");
  if(openinganglegammaproj) openinganglegammaproj->SetTitle("");
  if(openinganglepi0proj) openinganglepi0proj->SetTitle("");
  if(openingangleCproj) openingangleCproj->SetTitle("");
  if(openingangleBproj) openingangleBproj->SetTitle("");
  if(openingangleetaproj) openingangleetaproj->SetTitle("");
  if(openingangleSplittedTrackssproj) openingangleSplittedTrackssproj->SetTitle("");
  if(openingangleSplittedTrackosproj) openingangleSplittedTrackosproj->SetTitle("");

  invmassppproj->SetTitle("");
  invmassnnproj->SetTitle("");
  invmassssproj->SetTitle("");
  invmassrproj->SetTitle("");
  invmassosproj->SetTitle("");
  if(invmassgammaproj) invmassgammaproj->SetTitle("");
  if(invmasspi0proj) invmasspi0proj->SetTitle("");
  if(invmassCproj) invmassCproj->SetTitle("");
  if(invmassBproj) invmassBproj->SetTitle("");
  if(invmassetaproj) invmassetaproj->SetTitle("");
  if(invmassSplittedTrackssproj) invmassSplittedTrackssproj->SetTitle("");
  if(invmassSplittedTrackosproj) invmassSplittedTrackosproj->SetTitle("");

  // Projection pttagged variable
  TH2D *openingangleppproj2D = openinganglepp->Projection(0,2);
  TH2D *openinganglennproj2D = openinganglenn->Projection(0,2);
  TH2D *openinganglessproj2D = openingangless->Projection(0,2);
  TH2D *openinganglerproj2D  = openingangler->Projection(0,2);
  TH2D *openingangleosproj2D = openingangleos->Projection(0,2);

  TH2D *openinganglegammaproj2D = 0x0;
  TH2D *openinganglepi0proj2D = 0x0;
  TH2D *openingangleCproj2D = 0x0;
  TH2D *openingangleBproj2D = 0x0;
  TH2D *openingangleetaproj2D = 0x0;
  TH2D *openingangleSplittedTrackssproj2D = 0x0;
  TH2D *openingangleSplittedTrackosproj2D = 0x0;
  if(openinganglegamma) openinganglegammaproj2D = openinganglegamma->Projection(0,2);
  if(openinganglepi0) openinganglepi0proj2D = openinganglepi0->Projection(0,2);
  if(openingangleC) openingangleCproj2D = openingangleC->Projection(0,2);
  if(openingangleB) openingangleBproj2D = openingangleB->Projection(0,2);
  if(openingangleeta) openingangleetaproj2D = openingangleeta->Projection(0,2);
  if(openingangleSplittedTrackss) openingangleSplittedTrackssproj2D = openingangleSplittedTrackss->Projection(0,2);
  if(openingangleSplittedTrackos) openingangleSplittedTrackosproj2D = openingangleSplittedTrackos->Projection(0,2);


  TH2D *invmassppproj2D = invmasspp->Projection(0,2);
  TH2D *invmassnnproj2D = invmassnn->Projection(0,2);
  TH2D *invmassssproj2D = invmassss->Projection(0,2);
  TH2D *invmassrproj2D  = invmassr->Projection(0,2);
  TH2D *invmassosproj2D = invmassos->Projection(0,2);

  TH2D *invmassgammaproj2D = 0x0;
  TH2D *invmasspi0proj2D = 0x0;
  TH2D *invmassCproj2D = 0x0;
  TH2D *invmassBproj2D = 0x0;
  TH2D *invmassetaproj2D = 0x0;
  TH2D *invmassSplittedTrackssproj2D = 0x0;
  TH2D *invmassSplittedTrackosproj2D = 0x0;
  if(invmassgamma) invmassgammaproj2D = invmassgamma->Projection(0,2);
  if(invmasspi0) invmasspi0proj2D = invmasspi0->Projection(0,2);
  if(invmassC) invmassCproj2D = invmassC->Projection(0,2);
  if(invmassB) invmassBproj2D = invmassB->Projection(0,2);
  if(invmasseta) invmassetaproj2D = invmasseta->Projection(0,2);
  if(invmassSplittedTrackss) invmassSplittedTrackssproj2D = invmassSplittedTrackss->Projection(0,2);
  if(invmassSplittedTrackos) invmassSplittedTrackosproj2D = invmassSplittedTrackos->Projection(0,2);
  
  openingangleppproj2D->SetStats(0);
  openinganglennproj2D->SetStats(0);
  openinganglessproj2D->SetStats(0);
  openinganglerproj2D->SetStats(0);
  openingangleosproj2D->SetStats(0);
  if(openinganglegammaproj2D) openinganglegammaproj2D->SetStats(0);
  if(openinganglepi0proj2D) openinganglepi0proj2D->SetStats(0);
  if(openingangleCproj2D) openingangleCproj2D->SetStats(0);
  if(openingangleBproj2D) openingangleBproj2D->SetStats(0);
  if(openingangleetaproj2D) openingangleetaproj2D->SetStats(0);
  if(openingangleSplittedTrackssproj2D) openingangleSplittedTrackssproj2D->SetStats(0);
  if(openingangleSplittedTrackosproj2D) openingangleSplittedTrackosproj2D->SetStats(0);
  
  invmassppproj2D->SetStats(0);
  invmassnnproj2D->SetStats(0);
  invmassssproj2D->SetStats(0);
  invmassrproj2D->SetStats(0);
  invmassosproj2D->SetStats(0);
  if(invmassgammaproj2D) invmassgammaproj2D->SetStats(0);
  if(invmasspi0proj2D) invmasspi0proj2D->SetStats(0);
  if(invmassCproj2D) invmassCproj2D->SetStats(0);
  if(invmassBproj2D) invmassBproj2D->SetStats(0);
  if(invmassetaproj2D) invmassetaproj2D->SetStats(0);
  if(invmassSplittedTrackssproj2D) invmassSplittedTrackssproj2D->SetStats(0);
  if(invmassSplittedTrackosproj2D) invmassSplittedTrackosproj2D->SetStats(0);

  openingangleppproj2D->SetTitle("openingangleppproj2D");
  openinganglennproj2D->SetTitle("openinganglennproj2D");
  openinganglessproj2D->SetTitle("openinganglessproj2D");
  openinganglerproj2D->SetTitle("openinganglerproj2D");
  openingangleosproj2D->SetTitle("openingangleosproj2D");
  if(openinganglegammaproj2D) openinganglegammaproj2D->SetTitle("openinganglegammaproj2D");
  if(openinganglepi0proj2D) openinganglepi0proj2D->SetTitle("openinganglepi0proj2D");
  if(openingangleCproj2D) openingangleCproj2D->SetTitle("openingangleCproj2D");
  if(openingangleBproj2D) openingangleBproj2D->SetTitle("openingangleBproj2D");
  if(openingangleetaproj2D) openingangleetaproj2D->SetTitle("openingangleetaproj2D");
  if(openingangleSplittedTrackssproj2D) openingangleSplittedTrackssproj2D->SetTitle("openingangleSplittedTrackssproj2D");
  if(openingangleSplittedTrackosproj2D) openingangleSplittedTrackosproj2D->SetTitle("openingangleSplittedTrackosproj2D");

  invmassppproj2D->SetTitle("invmassppproj2D");
  invmassnnproj2D->SetTitle("invmassnnproj2D");
  invmassssproj2D->SetTitle("invmassssproj2D");
  invmassrproj2D->SetTitle("invmassrproj2D");
  invmassosproj2D->SetTitle("invmassosproj2D");
  if(invmassgammaproj2D) invmassgammaproj2D->SetTitle("invmassgammaproj2D");
  if(invmasspi0proj2D) invmasspi0proj2D->SetTitle("invmasspi0proj2D");
  if(invmassCproj2D) invmassCproj2D->SetTitle("invmassCproj2D");
  if(invmassBproj2D) invmassBproj2D->SetTitle("invmassBproj2D");
  if(invmassetaproj2D) invmassetaproj2D->SetTitle("invmassetaproj2D");
  if(invmassSplittedTrackssproj2D) invmassSplittedTrackssproj2D->SetTitle("invmassSplittedTrackssproj2D");
  if(invmassSplittedTrackosproj2D) invmassSplittedTrackosproj2D->SetTitle("invmassSplittedTrackosproj2D");


  // Draw histograms for opening angle
  TCanvas * copeningangle =new TCanvas("openingangle","Openingangle",800,800);
  copeningangle->cd();
  openingangleppproj->Draw();
  openinganglennproj->Draw("same");
  openinganglessproj->Draw("same");
  openinganglerproj->Draw("same");
  openingangleosproj->Draw("same");
  if(openinganglegammaproj) openinganglegammaproj->Draw("same");
  if(openinganglepi0proj) openinganglepi0proj->Draw("same");
  if(openingangleCproj) openingangleCproj->Draw("same");
  if(openingangleBproj) openingangleBproj->Draw("same");
  if(openingangleetaproj) openingangleetaproj->Draw("same");
  if(openingangleSplittedTrackssproj) openingangleSplittedTrackssproj->Draw("same");
  if(openingangleSplittedTrackosproj) openingangleSplittedTrackosproj->Draw("same");
  TLegend *lego = new TLegend(0.4,0.6,0.89,0.89);
  lego->AddEntry(openingangleppproj,"positive-positive","p");
  lego->AddEntry(openinganglennproj,"negative-negative","p");
  lego->AddEntry(openinganglessproj,"same-sign","p");
  lego->AddEntry(openinganglerproj,"rotated","p");
  lego->AddEntry(openingangleosproj,"positive-negative","p");
  if(openinganglegammaproj) lego->AddEntry(openinganglegammaproj,"e^{+}e^{-} from #gamma","p");
  if(openinganglepi0proj) lego->AddEntry(openinganglepi0proj,"e^{+}e^{-} from #pi^{0}","p");
  if(openingangleCproj) lego->AddEntry(openingangleCproj,"e^{+}e^{-} from c","p");
  if(openingangleBproj) lego->AddEntry(openingangleBproj,"e^{+}e^{-} from b","p");
  if(openingangleetaproj) lego->AddEntry(openingangleetaproj,"e^{+}e^{-} from #eta","p");
  if(openingangleSplittedTrackssproj) lego->AddEntry(openingangleSplittedTrackssproj,"Splitted tracks same sign","p");
  if(openingangleSplittedTrackosproj) lego->AddEntry(openingangleSplittedTrackosproj,"Splitted tracks opposite sign","p");
  lego->Draw("same");

  // Draw histograms for invariant mass
  TCanvas * cinvmass =new TCanvas("invmass","Invmass",800,800);
  cinvmass->cd();
  invmassppproj->Draw();
  invmassnnproj->Draw("same");
  invmassssproj->Draw("same");
  invmassrproj->Draw("same");
  invmassosproj->Draw("same");
  if(invmassgammaproj) invmassgammaproj->Draw("same");
  if(invmasspi0proj) invmasspi0proj->Draw("same");
  if(invmassCproj) invmassCproj->Draw("same");
  if(invmassBproj) invmassBproj->Draw("same");
  if(invmassetaproj) invmassetaproj->Draw("same");
  if(invmassSplittedTrackssproj) invmassSplittedTrackssproj->Draw("same");
  if(invmassSplittedTrackosproj) invmassSplittedTrackosproj->Draw("same");
  TLegend *legi = new TLegend(0.4,0.6,0.89,0.89);
  legi->AddEntry(invmassppproj,"positive-positive","p");
  legi->AddEntry(invmassnnproj,"negative-negative","p");
  legi->AddEntry(invmassssproj,"same-sign","p");
  legi->AddEntry(invmassrproj,"rotated","p");
  legi->AddEntry(invmassosproj,"positive-negative","p");
  if(invmassgammaproj) legi->AddEntry(invmassgammaproj,"e^{+}e^{-} from #gamma","p");
  if(invmasspi0proj) legi->AddEntry(invmasspi0proj,"e^{+}e^{-} from #pi^{0}","p");
  if(invmassCproj) legi->AddEntry(invmassCproj,"e^{+}e^{-} from c","p");
  if(invmassBproj) legi->AddEntry(invmassBproj,"e^{+}e^{-} from b","p");
  if(invmassetaproj) legi->AddEntry(invmassetaproj,"e^{+}e^{-} from #eta","p");
  if(invmassSplittedTrackssproj) legi->AddEntry(invmassSplittedTrackssproj,"Splitted tracks same sign","p");
  if(invmassSplittedTrackosproj) legi->AddEntry(invmassSplittedTrackosproj,"Splitted tracks opposite sign","p");
  legi->Draw("same");

  // Draw histograms for opening angle 2D
  TCanvas * copeningangle2D =new TCanvas("openingangle2D","Openingangle2D",800,800);
  copeningangle2D->Divide(6,2);
  copeningangle2D->cd(1);
  openingangleppproj2D->Draw("lego");
  copeningangle2D->cd(2);
  openinganglennproj2D->Draw("lego");
  copeningangle2D->cd(3);
  openinganglessproj2D->Draw("lego");
  copeningangle2D->cd(4);
  openinganglerproj2D->Draw("lego");
  copeningangle2D->cd(5);
  openingangleosproj2D->Draw("lego");
  copeningangle2D->cd(6);
  if(openinganglegammaproj2D) openinganglegammaproj2D->Draw("lego");
  copeningangle2D->cd(7);
  if(openinganglepi0proj2D) openinganglepi0proj2D->Draw("lego");
  copeningangle2D->cd(8);
  if(openingangleCproj2D) openingangleCproj2D->Draw("lego");
  copeningangle2D->cd(9);
  if(openingangleBproj2D) openingangleBproj2D->Draw("lego");
  copeningangle2D->cd(10);
  if(openingangleetaproj2D) openingangleetaproj2D->Draw("lego");
  copeningangle2D->cd(11);
  if(openingangleSplittedTrackssproj2D) openingangleSplittedTrackssproj2D->Draw("lego");
  copeningangle2D->cd(12);
  if(openingangleSplittedTrackosproj2D) openingangleSplittedTrackosproj2D->Draw("lego");
  
  // Draw histograms for invariant mass 2D
  TCanvas * cinvmass2D =new TCanvas("invmass2D","Invmass2D",800,800);
  cinvmass2D->Divide(6,2);
  cinvmass2D->cd(1);
  invmassppproj2D->Draw("lego");
  cinvmass2D->cd(2);
  invmassnnproj2D->Draw("lego");
  cinvmass2D->cd(3);
  invmassssproj2D->Draw("lego");
  cinvmass2D->cd(4);
  invmassrproj2D->Draw("lego");
  cinvmass2D->cd(5);
  invmassosproj2D->Draw("lego");
  cinvmass2D->cd(6);
  if(invmassgammaproj2D) invmassgammaproj2D->Draw("lego");
  cinvmass2D->cd(7);
  if(invmasspi0proj2D) invmasspi0proj2D->Draw("lego");
  cinvmass2D->cd(8);
  if(invmassCproj2D) invmassCproj2D->Draw("lego");
  cinvmass2D->cd(9);
  if(invmassBproj2D) invmassBproj2D->Draw("lego");
  cinvmass2D->cd(10);
  if(invmassetaproj2D) invmassetaproj2D->Draw("lego");
  cinvmass2D->cd(11);
  if(invmassSplittedTrackssproj2D) invmassSplittedTrackssproj2D->Draw("lego");
  cinvmass2D->cd(12);
  if(invmassSplittedTrackosproj2D) invmassSplittedTrackosproj2D->Draw("lego");
 
}
