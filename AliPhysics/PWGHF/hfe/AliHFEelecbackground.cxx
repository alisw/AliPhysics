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
#include <TH1F.h>
#include <TH2F.h>
#include <TF1.h>
#include <TString.h>
#include <TLegend.h>
#include <TStyle.h>

#include <AliESDEvent.h>
#include <AliAODEvent.h>
#include <AliESDtrack.h>
#include <AliAODTrack.h>
#include "AliHFEelecbackground.h"
#include <AliMCEvent.h>
#include <AliStack.h>
#include <AliKFParticle.h>
#include "AliCFContainer.h"
#include "AliHFEpid.h"
#include "AliESDpid.h"
#include "AliLog.h"
#include "AliITSPIDResponse.h"
#include "AliTPCPIDResponse.h"

ClassImp(AliHFEelecbackground)

Bool_t AliHFEelecbackground::fgUseMCPID = kFALSE;
const Double_t AliHFEelecbackground::fgkMe= 0.0005109989;


//___________________________________________________________________________________________
AliHFEelecbackground::AliHFEelecbackground():
  fhtmp(0x0)
  ,fhtmpf(0x0)
  ,fhtmpp(0x0)
  ,fESD1(0x0)
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
  ,fOpeningAngleCut(0.35)
  ,fInvMassCut(140.0)
  ,fChi2NdfCut(999999999.0)
  ,fUseAliKFCode(kTRUE)
  ,fSharedClusterCut(kFALSE)
  ,fRequireITSStandalone(0)
  ,fMinNbCls(2)
  ,fMinITSChi2(10.0)
  ,fMinNbClsSDDSPD(2)
  ,fPIDPartner(kFALSE)
  ,fPIDMethodPartner(0x0)
  ,fPIDMethodPartnerITS(0x0)
  ,fDebugLevel(0)
  ,fList(0x0)
  ,fListPostProcess(0x0)
{ 
  //
  // Default constructor
  //
  for(Int_t k =0; k < 10; k++) {
    fCuts[k] = kFALSE;
  }
  
}

//_______________________________________________________________________________________________
AliHFEelecbackground::AliHFEelecbackground(const AliHFEelecbackground &p):
  TObject(p)
  ,fhtmp(0x0)
  ,fhtmpf(0x0)
  ,fhtmpp(0x0)
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
  ,fOpeningAngleCut(0.35)
  ,fInvMassCut(140.0)
  ,fChi2NdfCut(999999999.0)
  ,fUseAliKFCode(kTRUE)
  ,fSharedClusterCut(kFALSE)
  ,fRequireITSStandalone(0)
  ,fMinNbCls(2)
  ,fMinITSChi2(10.0)
  ,fMinNbClsSDDSPD(2)
  ,fPIDPartner(kFALSE)
  ,fPIDMethodPartner(0x0)
  ,fPIDMethodPartnerITS(0x0)
  ,fDebugLevel(0)
  ,fList(0x0)  
  ,fListPostProcess(0x0)
{ 
  //
  // Copy constructor
  //
  for(Int_t k =0; k < 10; k++) {
    fCuts[k] = kFALSE;
  }
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
  if(fPIDMethodPartner) delete fPIDMethodPartner;
  if(fPIDMethodPartnerITS) delete fPIDMethodPartnerITS;

  if(fListPostProcess){
    fListPostProcess->SetOwner(kTRUE);
    delete fListPostProcess;
  }

/*
  if(fhtmp) delete fhtmp;
  if(fhtmpf) delete fhtmpf;
  if(fhtmpp) delete fhtmpp;
*/

}
//___________________________________________________________________________________________
Bool_t AliHFEelecbackground::Load(const Char_t * filename)
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
//___________________________________________________________________________________________
Bool_t AliHFEelecbackground::Load(TList * const outputlist)
{
  //
  // Generic container loader
  //
  if(!outputlist) return kFALSE;
  else   fList = (TList*)outputlist->Clone("HFEelecbackground");
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
  for(Int_t id = 0; id < 10; id++) {
    fCuts[id] = kFALSE;
  }
 
}
//_______________________________________________________________________________________________
void AliHFEelecbackground::CreateHistograms(TList * const qaList)
{ 
  //
  // create histograms
  //
  if(!qaList) return;

  fList = qaList;
  fList->SetName("HFEelecbackground");  

  //////////
  // bins
  /////////

  const Int_t nBinsPt = 25;
  Double_t minPt = 0.01;
  Double_t maxPt = 10.0;
  
  const Int_t nBinsPtMore = 100;
  Double_t minPtMore = 0.01;
  Double_t maxPtMore = 10.0;
  
  const Int_t nBinsInv = 50;
  Double_t minInv = 0.0;
  Double_t maxInv = 0.2;
  
  const Int_t nBinsOp = 50;
  Double_t minOp = 0.0;
  Double_t maxOp = 2;

  const Int_t nBinsCh = 4;
  Double_t minCh = 0.0;
  Double_t maxCh = 4.0;
  
  Double_t binLimLogPt[nBinsPt+1];
  Double_t binLimPt[nBinsPt+1];
  for(Int_t i=0; i<=nBinsPt; i++) binLimLogPt[i]=(Double_t)TMath::Log10(minPt) + (TMath::Log10(maxPt)-TMath::Log10(minPt))/nBinsPt*(Double_t)i ;
  for(Int_t i=0; i<=nBinsPt; i++) binLimPt[i]=(Double_t)TMath::Power(10,binLimLogPt[i]);

  Double_t binLimLogPtMore[nBinsPtMore+1];
  Double_t binLimPtMore[nBinsPtMore+1];
  for(Int_t i=0; i<=nBinsPtMore; i++) binLimLogPtMore[i]=(Double_t)TMath::Log10(minPtMore) + (TMath::Log10(maxPtMore)-TMath::Log10(minPtMore))/nBinsPtMore*(Double_t)i ;
  for(Int_t i=0; i<=nBinsPtMore; i++) binLimPtMore[i]=(Double_t)TMath::Power(10,binLimLogPtMore[i]);

  Double_t binLimInv[nBinsInv+1];
  for(Int_t i=0; i<=nBinsInv; i++) binLimInv[i]=(Double_t)minInv  + (maxInv-minInv)  /nBinsInv*(Double_t)i ;
  
  Double_t binLimOp[nBinsOp+1];
  for(Int_t i=0; i<=nBinsOp; i++) binLimOp[i]=(Double_t)minOp  + (maxOp-minOp) /nBinsOp*(Double_t)i ;
  
  Double_t binLimCh[nBinsCh+1];
  for(Int_t i=0; i<=nBinsCh; i++) binLimCh[i]=(Double_t)minCh  + (maxCh-minCh) /nBinsCh*(Double_t)i ;
  
  const Int_t nvarData = 5;
  // pt reconstructed tagged e
  // pt reconstructed mother
  // opening angle
  // invariant mass 
  // Data: charge (0-opposite sign, 1-like sign ++, 2-like sign --, 3-rotated tracks)

  Int_t iBinData[nvarData];
  iBinData[0]=nBinsPt;
  iBinData[1]=nBinsPt;
  iBinData[2]=nBinsOp;
  iBinData[3]=nBinsInv;
  iBinData[4]=nBinsCh;
  
  //
  // Opening angle and invariant mass
  //
  
  THnSparseF *hsSparseData = new THnSparseF("OpeningangleinvmassData","",nvarData,iBinData);
  hsSparseData->SetBinEdges(0,&binLimPt[0]);
  hsSparseData->SetBinEdges(1,&binLimPt[0]);
  hsSparseData->SetBinEdges(2,&binLimOp[0]);
  hsSparseData->SetBinEdges(3,&binLimInv[0]);
  hsSparseData->SetBinEdges(4,&binLimCh[0]);
  hsSparseData->Sumw2();

  fList->AddAt(hsSparseData,kDatai);

  //
  // Radius, DCA and Chi2Ndf
  //

  TH1F *dataRadiusHisto = new TH1F("DataRadius","", 200, 0.0, 200.0); // recontructed radius from the AliKF of the e+e- pair
  fList->AddAt(dataRadiusHisto,kDatar);

  TH1F *dataDcaHisto = new TH1F("DataDCA","", 100, 0.0, 6.0); // dca distribution
  fList->AddAt(dataDcaHisto,kDatadca); 
  
  TH1F *dataChi2NdfHisto = new TH1F("DataChi2Ndf","", 100, 0.0, 5.0); // chi2Ndf distribution    
  fList->AddAt(dataChi2NdfHisto,kDatachi2Ndf); 
  

  if(HasMCData()) {

    //
    // Opening angle and invariant mass with MC infos
    //

    const Int_t nvarMCo = 6;
    // pt reconstructed tagged e
    // pt reconstructed mother
    // opening angle
    // invariant mass 
    // MC: 0-FromBackground, 1-FromGamma, 2-FromPi0, 3-FromEta, 4-FromC, 5-FromB
    // MCSplitted: 0-not, 1-splittedOs, 2-ksplittedSs
    

    const Int_t nBinsMCOrigin = 6;
    Double_t minMCOrigin = 0.0;
    Double_t maxMCOrigin = 6.0;
    
    Double_t binLimMCOrigin[nBinsMCOrigin+1];
    for(Int_t i=0; i<=nBinsMCOrigin; i++) binLimMCOrigin[i]=(Double_t)minMCOrigin  + (maxMCOrigin-minMCOrigin) /nBinsMCOrigin*(Double_t)i ;

    const Int_t nBinsMCSplitted = 3;
    Double_t minMCSplitted = 0.0;
    Double_t maxMCSplitted = 3.0;
    
    Double_t binLimMCSplitted[nBinsMCSplitted+1];
    for(Int_t i=0; i<=nBinsMCSplitted; i++) binLimMCSplitted[i]=(Double_t)minMCSplitted  + (maxMCSplitted-minMCSplitted) /nBinsMCSplitted*(Double_t)i ;
    
    Int_t iBinMCo[nvarMCo];
    iBinMCo[0]=nBinsPt;
    iBinMCo[1]=nBinsPt;
    iBinMCo[2]=nBinsOp;
    iBinMCo[3]=nBinsInv;
    iBinMCo[4]=nBinsMCOrigin;
    iBinMCo[5]=nBinsMCSplitted;
        
    THnSparseF *hsSparseMCo = new THnSparseF("OpeningangleinvmassMC","",nvarMCo,iBinMCo);
    hsSparseMCo->SetBinEdges(0,&binLimPt[0]);
    hsSparseMCo->SetBinEdges(1,&binLimPt[0]);
    hsSparseMCo->SetBinEdges(2,&binLimOp[0]);
    hsSparseMCo->SetBinEdges(3,&binLimInv[0]);
    hsSparseMCo->SetBinEdges(4,&binLimMCOrigin[0]);
    hsSparseMCo->SetBinEdges(5,&binLimMCSplitted[0]);
    hsSparseMCo->Sumw2();

    fList->AddAt(hsSparseMCo,kMCo);

    //
    // Radius, DCA and Chi2Ndf with MC info
    //

    TH2F *mcRadiusHisto = new TH2F("MCRadius","", 200, 0.0, 200.0,6,-0.5,5.5); // recontructed radius from the AliKF of the e+e- pair
    fList->AddAt(mcRadiusHisto,kMCr);
    
    TH2F *mcDcaHisto = new TH2F("MCDCA","", 100, 0.0, 6.0,6,-0.5,5.5); // dca distribution
    fList->AddAt(mcDcaHisto,kMCdca); 
    
    TH2F *mcChi2NdfHisto = new TH2F("MCChi2Ndf","", 100, 0.0, 5.0,6,-0.5,5.5); // chi2Ndf distribution    
    fList->AddAt(mcChi2NdfHisto,kMCchi2Ndf); 

    //////////////////////////////////////////////////////////
    // if fDebugLevel 1: Rejection efficiencies of the cuts
    //////////////////////////////////////////////////////////

    if(fDebugLevel > 0) {

      if(HasMCData()) {
	
	const Int_t nvarMCe = 3;
	// pt reconstructed tagged e
	// cut passed: 0-all, 1-Partner tracked, 2-Opposite-sign, 3-SingleTrackCutPart, 4-ShareCluster, 5-PID, 6-DCA, 7-chi2Ndf AliKF, 8-Openingangle, 9-Invmass
	// MC: 0-FromBackground, 1-FromGamma, 2-FromPi0, 3-FromEta, 4-FromC, 5-FromB
      
      const Int_t nBinsMCCutPassed = 10;
      Double_t minMCCutPassed = -0.5;
      Double_t maxMCCutPassed = 9.5;
      
      Double_t binLimMCCutPassed[nBinsMCCutPassed+1];
      for(Int_t i=0; i<=nBinsMCCutPassed; i++) binLimMCCutPassed[i]=(Double_t)minMCCutPassed  + (maxMCCutPassed-minMCCutPassed) /nBinsMCCutPassed*(Double_t)i ;
      
      Int_t iBinMCe[nvarMCe];
      iBinMCe[0]=nBinsPt;
      iBinMCe[1]=nBinsMCCutPassed;
      iBinMCe[2]=nBinsMCOrigin;
      
      THnSparseF *hsSparseMCe = new THnSparseF("CutPassedMC","",nvarMCe,iBinMCe);
      hsSparseMCe->SetBinEdges(0,&binLimPt[0]);
      hsSparseMCe->SetBinEdges(1,&binLimMCCutPassed[0]);
      hsSparseMCe->SetBinEdges(2,&binLimMCOrigin[0]);
      hsSparseMCe->Sumw2();
      
      fList->AddAt(hsSparseMCe,kMCe); 
      
      }
    }

    /////////////////////////////////////////////////////////////////
    // if fDebugLevel 1: PIDPartCut and ShareClusters
    /////////////////////////////////////////////////////////////////

    if(fDebugLevel > 1) {

      if(!fRequireITSStandalone){
	
	//
	// TPC 
	//

	TH2F *tpcPartner0 = new TH2F("TPCPartner0","", nBinsPtMore, binLimPtMore, 200, 0.0, 700.0); 
	fList->AddAt(tpcPartner0,kMCcutPart0); 
	TH2F *tpcPartner1 = new TH2F("TPCPartner1","", nBinsPtMore, binLimPtMore, 200, 0.0, 700.0); 
	fList->AddAt(tpcPartner1,kMCcutPart1); 
      
      }
      else {

	//
	// ITS
	//

	TH2F *itsPartner0 = new TH2F("ITSPartner0","", nBinsPtMore, binLimPtMore, 200, 0.0, 700.0); 
	fList->AddAt(itsPartner0,kMCcutPart0); 
	TH2F *itsPartner1 = new TH2F("ITSPartner1","", nBinsPtMore, binLimPtMore, 200, 0.0, 700.0); 
	fList->AddAt(itsPartner1,kMCcutPart1); 

	/////////////////////////////////////////////////////
       	// dEdx of the four layers for the track partner
	/////////////////////////////////////////////////////
	const Int_t nvarITSsignal = 5;
	
	const Int_t nBinsITSsignal  = 100;
	Double_t minITSsignal = 0.0;
	Double_t maxITSsignal = 350.0;
	
	Double_t binLimITSsignal[nBinsITSsignal+1];
	for(Int_t i=0; i<=nBinsITSsignal; i++) binLimITSsignal[i]=(Double_t)minITSsignal  + (maxITSsignal-minITSsignal) /nBinsITSsignal*(Double_t)i ;
	
	Int_t iBinITSsignal[nvarITSsignal];
	iBinITSsignal[0]=nBinsPt;
	iBinITSsignal[1]=nBinsITSsignal;
	iBinITSsignal[2]=nBinsITSsignal;
	iBinITSsignal[3]=nBinsITSsignal;
	iBinITSsignal[4]=nBinsITSsignal;
	
	THnSparseF *hsSparseITSpid = new THnSparseF("SparseITSsignal","",nvarITSsignal,iBinITSsignal);
	hsSparseITSpid->SetBinEdges(0,&binLimPt[0]);
	hsSparseITSpid->SetBinEdges(1,&binLimITSsignal[0]);
	hsSparseITSpid->SetBinEdges(2,&binLimITSsignal[0]);
	hsSparseITSpid->SetBinEdges(3,&binLimITSsignal[0]);
	hsSparseITSpid->SetBinEdges(4,&binLimITSsignal[0]);
	hsSparseITSpid->Sumw2();
	
	fList->AddAt(hsSparseITSpid,kMCcutPart2); 

	///////////////////////////////////////////////////////////////////////////////////////
 	// dEdx of the four layers for the track partner and track to reject splitted track
	///////////////////////////////////////////////////////////////////////////////////////
	const Int_t nvarITSsignalSplit = 5;

	const Int_t nBinsITSSplit  = 2;
	Double_t minITSSplit = 0.0;
	Double_t maxITSSplit = 2.0;
	
	Double_t binLimITSSplit[nBinsITSSplit+1];
	for(Int_t i=0; i<=nBinsITSSplit; i++) binLimITSSplit[i]=(Double_t)minITSSplit  + (maxITSSplit-minITSSplit) /nBinsITSSplit*(Double_t)i ;


	const Int_t nBinsITSsignalSplit  = 50;
	Double_t minITSsignalSplit = -25.0;
	Double_t maxITSsignalSplit = 25.0;
	
	Double_t binLimITSsignalSplit[nBinsITSsignalSplit+1];
	for(Int_t i=0; i<=nBinsITSsignalSplit; i++) binLimITSsignalSplit[i]=(Double_t)minITSsignalSplit  + (maxITSsignalSplit-minITSsignalSplit) /nBinsITSsignalSplit*(Double_t)i ;
	
	Int_t iBinITSsignalSplit[nvarITSsignalSplit];
	iBinITSsignalSplit[0]=nBinsITSSplit;
	for(Int_t k = 1; k < 5; k++){
	  iBinITSsignalSplit[k]=nBinsITSsignalSplit;
	}
	
	THnSparseF *hsSparseITSpidSplit = new THnSparseF("SparseITSsignalSplit","",nvarITSsignalSplit,iBinITSsignalSplit);
	hsSparseITSpidSplit->SetBinEdges(0,&binLimITSSplit[0]);
	for(Int_t k = 1; k < 5; k++) {
	  hsSparseITSpidSplit->SetBinEdges(k,&binLimITSsignalSplit[0]);
	}
	hsSparseITSpidSplit->Sumw2();
	
	fList->AddAt(hsSparseITSpidSplit,kMCcutPart3); 
	
      }

    }
    
  }

  //qaList->Add(fList);

}
//_______________________________________________________________________________________________
void AliHFEelecbackground::PairAnalysis(AliESDtrack* const track, AliESDtrack* const trackPart)
{
  //
  // calculate (tagged e-partner) dca, opening angle, invariant mass 
  //

  /////////////////////
  // pt tagged
  //////////////////////
  TVector3 v3Dtagged;
  Double_t pxyz[3];
  track->PxPyPz(&pxyz[0]);
  v3Dtagged.SetXYZ(pxyz[0],pxyz[1],pxyz[2]);
  fPtESD = TMath::Sqrt(pxyz[0]*pxyz[0]+pxyz[1]*pxyz[1]); 
  
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

      for(Int_t id = 0; id < 10; id++) {
	fCuts[id] = kFALSE;
      }
      
      fIsFrom = kElectronFromBackground;
          
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
  
  /////////////////////////
  // Cut effects
  ////////////////////////   
  Double_t cuteffect[3];

  if(fDebugLevel > 0) {  
    if(HasMCData()) {
      cuteffect[0] = fPtESD;
      cuteffect[1] = 0.0;
      cuteffect[2] = fIsFrom;
      if(!fCuts[0]){
	if((fhtmp = dynamic_cast<THnSparseF *>(fList->At(kMCe)))) fhtmp->Fill(cuteffect);
	//if(fList->At(kMCe)) (dynamic_cast<THnSparseF *>(fList->At(kMCe)))->Fill(cuteffect);
	fCuts[0] = kTRUE;
      }
    }
  }


  ///////////////////////////////
  // Cut effect: Partner track 
  ///////////////////////////////   
  
  if(fDebugLevel > 0) {  
    if(HasMCData() && fIsPartner) {
      cuteffect[1] = 1.0;
      if(!fCuts[1]) {
	if((fhtmp = dynamic_cast<THnSparseF *>(fList->At(kMCe)))) fhtmp->Fill(cuteffect);
	//if(fList->At(kMCe)) (dynamic_cast<THnSparseF *>(fList->At(kMCe)))->Fill(cuteffect);
	fCuts[1] = kTRUE;
      }
    }
  }

  ///////////////////////////////
  // Cut effect: Opposite sign 
  ///////////////////////////////   
  
  if(fDebugLevel > 0) {  
    if(HasMCData() && fIsPartner && (sign == kOs)) {
      cuteffect[1] = 2.0;
      if(!fCuts[2]) {
	if((fhtmp = dynamic_cast<THnSparseF *>(fList->At(kMCe)))) fhtmp->Fill(cuteffect);
	//if(fList->At(kMCe)) (dynamic_cast<THnSparseF *>(fList->At(kMCe)))->Fill(cuteffect);
	fCuts[2] = kTRUE;
      }
    }
  }

  ////////////////////////
  // Partner track cut
  ////////////////////////
  if(!SingleTrackCut(trackPart)) return;

  if(fDebugLevel > 0) {  
    if(HasMCData() && fIsPartner && (sign==kOs)) {
      cuteffect[1] = 3.0;
      if(!fCuts[3]) {
	if((fhtmp = dynamic_cast<THnSparseF *>(fList->At(kMCe)))) fhtmp->Fill(cuteffect);
	//if(fList->At(kMCe)) (dynamic_cast<THnSparseF *>(fList->At(kMCe)))->Fill(cuteffect);
	fCuts[3] = kTRUE;
      }
    }
  }
  
  /////////////////////////
  // shared clusters cut
  /////////////////////////
  if(fSharedClusterCut && ShareCluster(track,trackPart)) return;

  if(fDebugLevel > 0) {  
    if(HasMCData() && fIsPartner && (sign==kOs)) {
      cuteffect[1] = 4.0;
      if(!fCuts[4]) {
	if((fhtmp = dynamic_cast<THnSparseF *>(fList->At(kMCe)))) fhtmp->Fill(cuteffect);
	//if(fList->At(kMCe)) (dynamic_cast<THnSparseF *>(fList->At(kMCe)))->Fill(cuteffect);
	fCuts[4] = kTRUE;
      }
    } 
  }

  ////////////////////////
  // PID Partner track 
  ////////////////////////
  if(!PIDTrackCut(trackPart)) return;

  if(fDebugLevel > 0) {  
    if(HasMCData() && fIsPartner && (sign==kOs)) {
      cuteffect[1] = 5.0;
      if(!fCuts[5]) {
	if((fhtmp = dynamic_cast<THnSparseF *>(fList->At(kMCe)))) fhtmp->Fill(cuteffect);
	//if(fList->At(kMCe)) (dynamic_cast<THnSparseF *>(fList->At(kMCe)))->Fill(cuteffect);
	fCuts[5] = kTRUE;
      }
    }
  }

  
  //////////////////////
  // DCA
  /////////////////////
  
  Double_t xthis,xp;
  Double_t dca = track->GetDCA(trackPart,fBz,xthis,xp);
  if((fhtmpp = dynamic_cast<TH1F *>(fList->At(kDatadca)))) fhtmpp->Fill(dca);
  if(HasMCData()) {
    //printf("has MC data for DCA\n");
    //printf("IsPartner %d and isfrom %d\n",fIsPartner,fIsFrom);
    if(fIsFrom==kElectronFromBackground) {
      if((fhtmpf = dynamic_cast<TH2F *>(fList->At(kMCdca)))) fhtmpf->Fill(dca,fIsFrom);
    }
    else {
      if(fIsPartner){
	if((fhtmpf = dynamic_cast<TH2F *>(fList->At(kMCdca)))) fhtmpf->Fill(dca,fIsFrom);
      }
    }
  }
   
  if(TMath::Abs(dca) > 3.0) return;
  
  if(fDebugLevel > 0) {  
    if(HasMCData() && fIsPartner && (sign==kOs)) {
      cuteffect[1] = 6.0;
      if(!fCuts[6]) {
	if((fhtmp = dynamic_cast<THnSparseF *>(fList->At(kMCe)))) fhtmp->Fill(cuteffect);
	//if(fList->At(kMCe)) (dynamic_cast<THnSparseF *>(fList->At(kMCe)))->Fill(cuteffect);
	fCuts[6] = kTRUE;
      }
    }
  }

  ///////////////////////////////////
  // Calcul mother variables
  ///////////////////////////////////
  Double_t results[5];
  Double_t resultsr[5];

  
  if(!fUseAliKFCode) {
    
    /////////////////////////////
    // Propagate
    ////////////////////////////
    
    Double_t norradius = TMath::Sqrt(fkVertex->GetX()*fkVertex->GetX()+fkVertex->GetY()*fkVertex->GetY());
    
    AliESDtrack trackCopy = AliESDtrack(*track);
    AliESDtrack trackPartCopy = AliESDtrack(*trackPart);
    Bool_t propagateok = kTRUE;
    if((!(trackPartCopy.PropagateTo(norradius,fBz))) || (!(trackCopy.PropagateTo(norradius,fBz)))) propagateok = kFALSE;
    if(!propagateok) {
      //if(trackCopy) delete trackCopy;
      //if(trackPartCopy) delete trackPartCopy;
      return;
    }  
  
    CalculateMotherVariable(&trackCopy,&trackPartCopy,&results[0]);
    CalculateMotherVariableR(&trackCopy,&trackPartCopy,&resultsr[0]);
    
    //if(trackCopy) delete trackCopy;
    //if(trackPartCopy) delete trackPartCopy;
    
  }
  else {
    
    if(!CalculateMotherVariable(track,trackPart,&results[0])) return;
    if(fDebugLevel > 0) {     
      if(HasMCData() && fIsPartner && (sign==kOs)) {
	cuteffect[1] = 7.0;
	if(!fCuts[7]) {
	  if((fhtmp = dynamic_cast<THnSparseF *>(fList->At(kMCe)))) fhtmp->Fill(cuteffect);
	  //if(fList->At(kMCe)) (dynamic_cast<THnSparseF *>(fList->At(kMCe)))->Fill(cuteffect);
	  fCuts[7] = kTRUE;
	}
      }
    }
  
  }
  
  /////////////////////////////////////
  // Fill
  /////////////////////////////////////
   
  FillOutput(results, resultsr, sign);

  if(fDebugLevel > 0) {  
    if(HasMCData() && fIsPartner && (sign==kOs)) {
      
      if(TMath::Abs(results[4]) < fOpeningAngleCut) {
	
	cuteffect[1] = 8.0;
	if(!fCuts[8]) {
	  if((fhtmp = dynamic_cast<THnSparseF *>(fList->At(kMCe)))) fhtmp->Fill(cuteffect);
	  //if(fList->At(kMCe)) (dynamic_cast<THnSparseF *>(fList->At(kMCe)))->Fill(cuteffect);
	  fCuts[8] = kTRUE;
	}
	if(TMath::Abs(results[1]) < fInvMassCut) {
	  cuteffect[1] = 9.0;
	  if(!fCuts[9]) {
	    if((fhtmp = dynamic_cast<THnSparseF *>(fList->At(kMCe)))) fhtmp->Fill(cuteffect);
	    //if(fList->At(kMCe)) (dynamic_cast<THnSparseF *>(fList->At(kMCe)))->Fill(cuteffect);
	    fCuts[9] = kTRUE;
	  }
	}
      }
    }    
  }
  
 
}
//_____________________________________________________________________________________
Bool_t AliHFEelecbackground::CalculateMotherVariable(AliESDtrack* const track, AliESDtrack* const trackpart, Double_t *results)
{
  //
  // variables mother and take the pt of the track
  //
  // results contain: ptmother, invmass, etamother, phimother, openingangle
  //
  // with a chi2Ndf cut for AliKF code
  //
  
  if(!fUseAliKFCode) {
    
    TVector3 v3Dtagged;
    TVector3 v3Dpart;
    
    Double_t pxyz[3];
    track->PxPyPz(&pxyz[0]);
    v3Dtagged.SetXYZ(pxyz[0],pxyz[1],pxyz[2]);
    
    Double_t pxyzpart[3];
    trackpart->PxPyPz(&pxyzpart[0]);
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
    
    if(!results) return kFALSE;

    results[0] = ptESDmother;
    results[1] = etaESDmother;
    results[2] = phiESDmother;
    results[3] = invmass;
    results[4] = openingangle;
    
    return kTRUE;

  }
  else {
    
    AliKFParticle pair;
    pair.Initialize();
    
    // pid
    Int_t pid1 = -11;
    if(track->Charge() > 0.0) pid1 = 11;
    Int_t pid2 = -11;
    if(trackpart->Charge() > 0.0) pid2 = 11;
    
    
    // daughters
    AliKFParticle kf(*track,pid1);
    AliKFParticle kfpart(*trackpart,pid2);
    
    pair.AddDaughter(kf);
    pair.AddDaughter(kfpart);
    
    // variables
    Double_t openingangle = kf.GetAngle(kfpart);
    Double_t chi2ndf = pair.GetChi2()/pair.GetNDF();
    //Double_t decayLength = pair.GetDecayLength();
    Double_t radius = pair.GetR();
    //Double_t masserror = pair.GetErrMass()>0?pair.GetErrMass()/pair.GetMass():1000000;
    Double_t ptpair = pair.GetPt();
    Double_t etapair = pair.GetEta();
    Double_t phipair = pair.GetPhi();
    Double_t masspair = pair.GetMass();
    
    // Put them
    if(!results) return kFALSE;

    results[0] = ptpair;
    results[1] = etapair;
    results[2] = phipair;
    results[3] = masspair;
    results[4] = openingangle;

    // chi2Ndf cut
    if((fhtmpp = dynamic_cast<TH1F *>(fList->At(kDatachi2Ndf)))) fhtmpp->Fill(chi2ndf);
    //if(fList->At(kDatachi2Ndf)) (dynamic_cast<TH1F *>(fList->At(kDatachi2Ndf)))->Fill(chi2ndf);
    if(HasMCData()){
      if(fIsFrom==kElectronFromBackground) {
	if((fhtmpf = dynamic_cast<TH2F *>(fList->At(kMCchi2Ndf)))) fhtmpf->Fill(chi2ndf,fIsFrom); 
      }
      else {
	if(fIsPartner){
	  if((fhtmpf = dynamic_cast<TH2F *>(fList->At(kMCchi2Ndf)))) fhtmpf->Fill(chi2ndf,fIsFrom); 
	}
      }
    }
    if(chi2ndf > fChi2NdfCut) return kFALSE;
    else {
      if((fhtmpp = dynamic_cast<TH1F *>(fList->At(kDatar)))) fhtmpp->Fill(radius); 
      //if(fList->At(kDatar)) (dynamic_cast<TH1F *>(fList->At(kDatar)))->Fill(radius);
      if(HasMCData()) {
	if(fIsFrom==kElectronFromBackground) {
	  if((fhtmpf = dynamic_cast<TH2F *>(fList->At(kMCr)))) fhtmpf->Fill(radius,fIsFrom); 
	  //if(fList->At(kMCr))) (dynamic_cast<TH2F *>(fList->At(kMCr)))->Fill(radius,fIsFrom);
	}
	else {
	  if(fIsPartner) {
	    if((fhtmpf = dynamic_cast<TH2F *>(fList->At(kMCr)))) fhtmpf->Fill(radius,fIsFrom); 
	  }
	}
      }
      return kTRUE;
    }
    
  }
  
}
//_____________________________________________________________________________________
void AliHFEelecbackground::CalculateMotherVariableR(AliESDtrack* const track, AliESDtrack* const trackpart, Double_t *results)
{
  //
  // variables mother
  //
  // results contain: ptmother, invmass, etamother, phimother, openingangle
  // Implemented only for no AliKF
  //

  if(!fUseAliKFCode) {
    
    TVector3 v3Dtagged;
    TVector3 v3Dpart;
    
    Double_t pxyz[3];
    track->PxPyPz(&pxyz[0]);
    v3Dtagged.SetXYZ(pxyz[0],pxyz[1],pxyz[2]);
    Double_t pxyzpart[3];
    trackpart->PxPyPz(&pxyzpart[0]);
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
    
    if(!results) return;
    
    results[0] = ptESDmother;
    results[1] = etaESDmother;
    results[2] = phiESDmother;
    results[3] = invmass;
    results[4] = openingangle;
    
  }
  
}
//_________________________________________________________________________________
void AliHFEelecbackground::FillOutput(const Double_t *results, const Double_t *resultsr,Int_t sign) 
{
  //
  // Fill the Data and MC THnSparseF 
  //

  if((!results) || (!resultsr)) return;
  
  Double_t co[6];
  co[0] = fPtESD;
  co[1] = results[0];
  co[2] = TMath::Abs(results[4]);
  co[3] = results[3];
  co[4] = sign;
  co[5] = 0.0;

  if((fhtmp = dynamic_cast<THnSparseF *>(fList->At(kDatai)))) fhtmp->Fill(co);
  //if(fList->At(kDatai))(dynamic_cast<THnSparseF *>(fList->At(kDatai)))->Fill(co);

  if((sign==kOs) && (!fUseAliKFCode)) {
    
    co[1] = resultsr[0];
    co[2] = TMath::Abs(resultsr[4]);
    co[3] = resultsr[3];
    co[4] = kR;
    co[5] = 0.0;

    if((fhtmp = dynamic_cast<THnSparseF *>(fList->At(kDatai)))) fhtmp->Fill(co);
    //if(fList->At(kDatai)) (dynamic_cast<THnSparseF *>(fList->At(kDatai)))->Fill(co);
    
  }
  
  if(HasMCData()){

    // Reset
    co[1] = results[0];
    co[2] = TMath::Abs(results[4]);
    co[3] = results[3];

    // Origin
    co[4] = kElectronFromBackground;
    if((sign==kOs) && fIsPartner) co[4] = fIsFrom;
    
    // Splitted tracks
    co[5] = kNotSplitted;
    if(fIsSplittedTrack) {
      if(sign==kOs){
	co[5] = kSplittedOs;
      }
      else {
	co[5] = kSplittedSs;
      }
    }

    if((fhtmp = dynamic_cast<THnSparseF *>(fList->At(kMCo)))) fhtmp->Fill(co);
    //if(fList->At(kMCo)) (dynamic_cast<THnSparseF *>(fList->At(kMCo)))->Fill(co);

  }
 
}    
//_______________________________________________________________________________________________
Bool_t AliHFEelecbackground::SingleTrackCut(const AliESDtrack* const trackPart) const
{
  //
  // Return minimum quality for the partner
  //
  
  if(trackPart->GetKinkIndex(0)>0) return kFALSE;


  UInt_t status = trackPart->GetStatus();

  if(fRequireITSStandalone > 0) {

    /////////////////////
    // ITS Standalone
    ////////////////////
    
    if(fRequireITSStandalone==1) {
      if(((status & AliESDtrack::kITSin) == 0 || (trackPart->IsPureITSStandalone()) || ((status&AliESDtrack::kITSrefit)==0))) return kFALSE;
    }
    
    if(fRequireITSStandalone==2) {
      if(!trackPart->IsPureITSStandalone() || ((status&AliESDtrack::kITSrefit)==0)) return kFALSE;
    }
    
    // Chi2
    Double_t chi2 = trackPart->GetITSchi2();
    if(chi2 > fMinITSChi2) return kFALSE;

    // Min Nb of clusters
    Int_t nbcl = trackPart->GetITSclusters(0);
    if(nbcl < fMinNbCls)  return kFALSE;  

    // Min Nb of points in SDD and SPD
    Int_t nbSDDSPD = 0;
    for(Int_t layer = 0; layer < 4; layer++){
      if(trackPart->HasPointOnITSLayer(layer)) nbSDDSPD++;
    }
    if(nbSDDSPD < fMinNbClsSDDSPD) return kFALSE;
    
    
  }
  else {

    /////////
    // TPC
    /////////
    
    if((status&AliESDtrack::kTPCrefit)==0) return kFALSE;

    // Min Nb of clusters
    Int_t nbcl = trackPart->GetTPCclusters(0);
    if(nbcl < fMinNbCls) return kFALSE;   
    
  }

  return kTRUE;

}  
//_______________________________________________________________________________________________
Bool_t AliHFEelecbackground::PIDTrackCut(AliESDtrack* const trackPart)
{
  //
  // PID for the partner using TPC or ITS
  //
  
  if(fRequireITSStandalone > 0) {

    /////////////////////
    // ITS Standalone
    ////////////////////
    
    // signal
    Double_t itsSignal = trackPart->GetITSsignal();
    Double_t p = trackPart->P();
    
    if(fDebugLevel > 1) {    
      if((fhtmpf = dynamic_cast<TH2F *>(fList->At(kMCcutPart0)))) fhtmpf->Fill(p,itsSignal);     
      //if(fList->At(kMCcutPart0)) (dynamic_cast<TH2F *>(fList->At(kMCcutPart0)))->Fill(p,itsSignal);
    }

    ///////////
    // PID
    //////////
    if(fPIDPartner) {
      
      // Take signal trackPart
      Double_t dEdxSamplesPart[4];
      trackPart->GetITSdEdxSamples(dEdxSamplesPart);

      // Cut at 2 sigma
      if(!fPIDMethodPartnerITS) fPIDMethodPartnerITS = new AliESDpid;
      Float_t nsigma = fPIDMethodPartnerITS->NumberOfSigmasITS(trackPart, AliPID::kElectron);
      if(TMath::Abs(nsigma) > 2.0) return kFALSE;
      
      // fill signal
      if(fDebugLevel > 1) {        
	
	Double_t entries[5];
	entries[0] = p;
	entries[1] = dEdxSamplesPart[0];
	entries[2] = dEdxSamplesPart[1];
	entries[3] = dEdxSamplesPart[2];
	entries[4] = dEdxSamplesPart[3];

	//if(fList->At(kMCcutPart1)) (dynamic_cast<TH2F *>(fList->At(kMCcutPart1)))->Fill(p,itsSignal);
	if((fhtmpf = dynamic_cast<TH2F *>(fList->At(kMCcutPart1)))) fhtmpf->Fill(p,itsSignal); 
	if((fhtmp = dynamic_cast<THnSparseF *>(fList->At(kMCcutPart2)))) fhtmp->Fill(entries);
	//if(fList->At(kMCcutPart2)) (dynamic_cast<THnSparseF *>(fList->At(kMCcutPart2)))->Fill(entries);
	
      }
      
    }
    
  }
  else {

    /////////
    // TPC
    /////////
    
    Double_t tpcSignal = trackPart->GetTPCsignal();
    Double_t p = trackPart->GetInnerParam() ? trackPart->GetInnerParam()->P() : trackPart->P();

    if(fDebugLevel > 1) {        
      //printf("tpcSignal %f\n",tpcSignal);
      //if(fList->At(kMCcutPart0)) (dynamic_cast<TH2F *>(fList->At(kMCcutPart0)))->Fill(p,tpcSignal);
      if((fhtmpf = dynamic_cast<TH2F *>(fList->At(kMCcutPart0)))) fhtmpf->Fill(p,tpcSignal); 
    }

    // PID
    if(fPIDPartner) {
      if(!fPIDMethodPartner) return kFALSE;
      AliHFEpidObject hfetrack;
      hfetrack.SetAnalysisType(AliHFEpidObject::kESDanalysis);
      hfetrack.SetRecTrack(trackPart);
      //if(HasMCData()) hfetrack.fMCtrack = mctrack;
      if(!fPIDMethodPartner->IsSelected(&hfetrack)) return kFALSE;
      
      if(fDebugLevel > 1) {  
	if((fhtmpf = dynamic_cast<TH2F *>(fList->At(kMCcutPart1)))) fhtmpf->Fill(p,tpcSignal); 
	//if(fList->At(kMCcutPart1)) (dynamic_cast<TH2F *>(fList->At(kMCcutPart1)))->Fill(p,tpcSignal);
      }
    }   
    
  }

  return kTRUE;

}
//__________________________________________________________________________________________
Bool_t AliHFEelecbackground::ShareCluster(AliESDtrack * const track1,AliESDtrack * const track2) 
{
  //
  // Look if the two tracks shared clusters in the TPC or in the ITS depending on the method
  //
  // For TPC:
  // hsfval: number of shared clusters 
  // hsmval: quality of the tracks
  //
  // For ITS:
  // compare the dEdx in the ITS
  //

  if(!fRequireITSStandalone) {

    //////////
    // TPC
    //////////

    
    Int_t nh = 0;
    Int_t an = 0;
    Int_t ns = 0;
    Float_t hsmval = 0.0;
    Float_t hsfval = 0.0;
   
    for(unsigned int imap=0;imap<track1->GetTPCClusterMap().GetNbits(); imap++) {
      if(track1->GetTPCClusterMap().TestBitNumber(imap) &&
         track2->GetTPCClusterMap().TestBitNumber(imap)) {
        // Do they share it ?
        if (track1->GetTPCSharedMap().TestBitNumber(imap) &&
            track2->GetTPCSharedMap().TestBitNumber(imap))
          {
            an++;
            nh+=2;
            ns+=2;
          }
	else {
          an--;
          nh+=2;
        }
      }
      else if (track1->GetTPCClusterMap().TestBitNumber(imap) ||
               track2->GetTPCClusterMap().TestBitNumber(imap)) {
	an++;
        nh++;
      }
    }

    if (nh >0) {
      hsmval = an*1.0/nh;
      hsfval = ns*1.0/nh;
    }
    
    
    if((hsfval > 0.15) || (hsmval > -0.15)) return kTRUE; //they share cluster 
    else return kFALSE;   
 
  
  }
  else {

    
    //////////
    // ITS
    /////////

    // Take signals 
    Double_t dEdxSamples1[4];
    track1->GetITSdEdxSamples(dEdxSamples1);
    Double_t dEdxSamples2[4];
    track2->GetITSdEdxSamples(dEdxSamples2);
    
    // If there are matching
    Int_t nbClusters = 0;
    Bool_t match[4] = {kTRUE,kTRUE,kTRUE,kTRUE};
    Double_t limit[4] = {1.5,1.5,1.5,1.5};
    for(Int_t layer = 0; layer < 4; layer++) {
      if(track1->HasPointOnITSLayer(layer+2) && track2->HasPointOnITSLayer(layer+2)) {
	if(TMath::Abs(dEdxSamples1[layer]-dEdxSamples2[layer])>limit[layer]) match[layer] = kFALSE;
	nbClusters++;
      }
    }
    //printf("nbClusters %d\n",nbClusters);
    
    // fill signal   
    if(fDebugLevel > 1) {        
      Double_t entriesSplit[5];
      entriesSplit[0] = 0.0;
      if(fIsSplittedTrack) entriesSplit[0] = 1.0; 
      
      for(Int_t layer = 0; layer < 4; layer++) {
	if(track1->HasPointOnITSLayer(layer+2) && track2->HasPointOnITSLayer(layer+2)) {
	  entriesSplit[layer+1] = dEdxSamples1[layer]-dEdxSamples2[layer];
	}
	else entriesSplit[layer+1] = -100.0;
      }
      if((fhtmp = dynamic_cast<THnSparseF *>(fList->At(kMCcutPart3)))) fhtmp->Fill(entriesSplit);
      //if(fList->At(kMCcutPart3)) (dynamic_cast<THnSparseF *>(fList->At(kMCcutPart3)))->Fill(entriesSplit);
    }

    // Return
    Int_t nbClustersNotClose = 0;
    for(Int_t layer = 0; layer < 4; layer++) {
      if(!match[layer]) nbClustersNotClose++;
    }
    if((nbClusters > 1) && (nbClustersNotClose > 0.75*nbClusters)) return kFALSE;
    else return kTRUE;
    
  }   
  
}
//____________________________________________________________________________________________________________
void AliHFEelecbackground::SetPIDPartner() {

  //
  // Init the stuff for PID on the partner track
  //

  fPIDPartner = kTRUE;

  if(fRequireITSStandalone == 0) {
    
    if(!fPIDMethodPartner) {
      fPIDMethodPartner = new AliHFEpid();
      fPIDMethodPartner->AddDetector("TPC", 0);
      fPIDMethodPartner->InitializePID();     // 3 sigma cut in TPC
    }

  }
  else {
    
    if(!fPIDMethodPartnerITS) fPIDMethodPartnerITS = new AliESDpid;
    
  }
  
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

  gStyle->SetPalette(1);
  gStyle->SetOptStat(1111);
  gStyle->SetPadBorderMode(0);
  gStyle->SetCanvasColor(10);
  gStyle->SetPadLeftMargin(0.13);
  gStyle->SetPadRightMargin(0.13);

  /////////////////////////
  // Take the THnSparseF
  /////////////////////////
  THnSparseF *hsSparseData = dynamic_cast<THnSparseF *>(fList->FindObject("OpeningangleinvmassData")); 
  THnSparseF *hsSparseMC = dynamic_cast<THnSparseF *>(fList->FindObject("OpeningangleinvmassMC")); 
  THnSparseF *hsSparseCutPassedMC = dynamic_cast<THnSparseF *>(fList->FindObject("CutPassedMC")); 

  /////////////////////////////////
  // Cuts on the opening angle
  ////////////////////////////////
  if(!hsSparseData) return;
  TAxis *axisOpeningAngleData = hsSparseData->GetAxis(2);
  Int_t binCutData = axisOpeningAngleData->FindBin(fOpeningAngleCut);
  hsSparseData->GetAxis(2)->SetRange(1,binCutData);

  if(hsSparseMC) {
    TAxis *axisOpeningAngleMC = hsSparseMC->GetAxis(2);
    Int_t binCutMC = axisOpeningAngleMC->FindBin(fOpeningAngleCut);
    hsSparseMC->GetAxis(2)->SetRange(1,binCutMC);
  }
  
  /////////////////////////
  // Prepare the histos
  ////////////////////////  

  TAxis *ptaxisinvmass = hsSparseData->GetAxis(3);
  Int_t  nbinsptinvmass = ptaxisinvmass->GetNbins();  
  
  TH1D **invmassosptproj = new TH1D*[nbinsptinvmass];
  TH1D **invmassssptproj = new TH1D*[nbinsptinvmass];
  TH1D **invmassrptproj = new TH1D*[nbinsptinvmass];
  TH1D **invmassdiffptproj = new TH1D*[nbinsptinvmass];
  TH1D **invmassgammaptproj = new TH1D*[nbinsptinvmass];
  TH1D **invmasspi0ptproj = new TH1D*[nbinsptinvmass];
  TH1D **invmassetaptproj = new TH1D*[nbinsptinvmass];
  TH1D **invmassCptproj = new TH1D*[nbinsptinvmass];
  TH1D **invmassBptproj = new TH1D*[nbinsptinvmass];

  TH1D *yieldPtFound = (TH1D *) hsSparseData->Projection(0);
  yieldPtFound->SetName("Found yield");
  yieldPtFound->Reset();

  TH1D *yieldPtSourcesMC = 0x0;
  TH1D *yieldPtSignalCutMC = 0x0;
  if(hsSparseMC) {
    yieldPtSourcesMC = (TH1D *) hsSparseMC->Projection(0);
    yieldPtSourcesMC->SetName("Found yield");
    yieldPtSourcesMC->Reset();
    
    yieldPtSignalCutMC = (TH1D *) hsSparseMC->Projection(0);
    yieldPtSignalCutMC->SetName("Found yield");
    yieldPtSignalCutMC->Reset();
  }
  
  ////////////
  // canvas
  ///////////
  Int_t nbrow = (Int_t) (nbinsptinvmass/5);
  TString namecanvas("InvMassSpectra");
  TCanvas * canvas =new TCanvas(namecanvas,namecanvas,800,800);
  canvas->Divide(5,nbrow+1);

  /////////////////////////////
  // Loop on pt bins
  /////////////////////////////

  for(Int_t k=1; k <= nbinsptinvmass; k++){

    Double_t lowedge = ptaxisinvmass->GetBinLowEdge(k);
    Double_t upedge  = ptaxisinvmass->GetBinUpEdge(k);

    // Pt bin
    hsSparseData->GetAxis(0)->SetRange(k,k);
    if(hsSparseMC) hsSparseMC->GetAxis(0)->SetRange(k,k);
    
    // 
    hsSparseData->GetAxis(4)->SetRange(kOs+1,kOs+1);
    invmassosptproj[k-1] = hsSparseData->Projection(3);
    hsSparseData->GetAxis(4)->SetRange(kPp+1,kNn+1);
    invmassssptproj[k-1] = hsSparseData->Projection(3);
    hsSparseData->GetAxis(4)->SetRange(kR+1,kR+1);
    invmassrptproj[k-1]  = hsSparseData->Projection(3);
    hsSparseData->GetAxis(4)->SetRange(1,hsSparseData->GetAxis(4)->GetNbins()); 
    invmassgammaptproj[k-1] = 0x0;
    invmasspi0ptproj[k-1] = 0x0;
    invmassetaptproj[k-1] = 0x0;
    invmassCptproj[k-1] = 0x0;
    invmassBptproj[k-1] = 0x0;
    if(hsSparseMC) {
      hsSparseMC->GetAxis(4)->SetRange(kElectronFromGamma+1,kElectronFromGamma+1);
      invmassgammaptproj[k-1] = hsSparseMC->Projection(3);
      hsSparseMC->GetAxis(4)->SetRange(kElectronFromPi0+1,kElectronFromPi0+1);
      invmasspi0ptproj[k-1] = hsSparseMC->Projection(3);
      hsSparseMC->GetAxis(4)->SetRange(kElectronFromEta+1,kElectronFromEta+1);
      invmassetaptproj[k-1] = hsSparseMC->Projection(3);
      hsSparseMC->GetAxis(4)->SetRange(kElectronFromC+1,kElectronFromC+1);
      invmassCptproj[k-1] = hsSparseMC->Projection(3);
      hsSparseMC->GetAxis(4)->SetRange(kElectronFromB+1,kElectronFromB+1);
      invmassBptproj[k-1] = hsSparseMC->Projection(3);
      hsSparseMC->GetAxis(4)->SetRange(1,hsSparseMC->GetAxis(4)->GetNbins()); 
    }      

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
      if(yieldPtSourcesMC) yieldPtSourcesMC->SetBinContent(k,yieldg);
      
      Double_t yieldsignal = invmassCptproj[k-1]->Integral() + invmassBptproj[k-1]->Integral();
      if(yieldPtSignalCutMC) yieldPtSignalCutMC->SetBinContent(k,yieldsignal);
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
  
    hsSparseData->GetAxis(0)->SetRange(1,hsSparseData->GetAxis(0)->GetNbins()); 
    if(hsSparseMC) hsSparseMC->GetAxis(0)->SetRange(1,hsSparseMC->GetAxis(0)->GetNbins()); 

  }

  ////////////////////////////////////////////////////
  // End of plotting: do subtraction of background
  ///////////////////////////////////////////////////

  yieldPtFound->SetStats(0);
  if(yieldPtSourcesMC) yieldPtSourcesMC->SetStats(0); 
  if(yieldPtSignalCutMC) yieldPtSignalCutMC->SetStats(0);  

  TCanvas * canvasfin =new TCanvas("ResultsElecBackGround","ResultsElecBackGround",800,800);
  canvasfin->cd(1);
  yieldPtFound->Draw();
  if(yieldPtSourcesMC && yieldPtSignalCutMC) {
    yieldPtSourcesMC->Draw("same");
    yieldPtSignalCutMC->Draw("same");
    TLegend *lega = new TLegend(0.4,0.6,0.89,0.89);
    lega->AddEntry(yieldPtFound,"Contributions found","l"); 
    lega->AddEntry(yieldPtSourcesMC,"Contributions of e^{+}e^{-} from #gamma, #pi^{0} and #eta","l"); 
    lega->AddEntry(yieldPtSignalCutMC,"Contributions of e^{+}e^{-} from C and B","l"); 
    lega->Draw("same");
  }
  
  if(hsSparseCutPassedMC){
    hsSparseCutPassedMC->GetAxis(1)->SetRange(1,1); 
    hsSparseCutPassedMC->GetAxis(2)->SetRange(1,4); 
    TH1D *hsSparseCutPassedMCproj = hsSparseCutPassedMC->Projection(0);

    TH1D *cYieldPtFound = (TH1D*)yieldPtFound->Clone("RatioEfficiency");
    if(hsSparseCutPassedMCproj->Integral() > 0.0) cYieldPtFound->Divide(hsSparseCutPassedMCproj);

    TCanvas * canvasfratio =new TCanvas("RatioEfficiency","RatioEfficiency",800,800);
    canvasfratio->cd(1);
    cYieldPtFound->Draw();
  }

  //////////////////////////
  // fListPostProcess
  /////////////////////////
  
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

  gStyle->SetPalette(1);
  gStyle->SetOptStat(1111);
  gStyle->SetPadBorderMode(0);
  gStyle->SetCanvasColor(10);
  gStyle->SetPadLeftMargin(0.13);
  gStyle->SetPadRightMargin(0.13);


  /////////////////////////
  // Take the THnSparseF
  /////////////////////////
  THnSparseF *hsSparseData = dynamic_cast<THnSparseF *>(fList->FindObject("OpeningangleinvmassData")); 
  THnSparseF *hsSparseMC = dynamic_cast<THnSparseF *>(fList->FindObject("OpeningangleinvmassMC")); 
  if(!hsSparseData) return;
  
  ////////////////////
  // Opening angle
  ////////////////////

  // Opening angle one direction
  hsSparseData->GetAxis(4)->SetRange(kPp+1,kPp+1);  
  TH1D *openingangleppproj = hsSparseData->Projection(2);
  hsSparseData->GetAxis(4)->SetRange(kNn+1,kNn+1);  
  TH1D *openinganglennproj = hsSparseData->Projection(2);
  hsSparseData->GetAxis(4)->SetRange(kPp+1,kNn+1);  
  TH1D *openinganglessproj = hsSparseData->Projection(2);
  hsSparseData->GetAxis(4)->SetRange(kR+1,kR+1);  
  TH1D *openinganglerproj  = hsSparseData->Projection(2);
  hsSparseData->GetAxis(4)->SetRange(kOs+1,kOs+1);  
  TH1D *openingangleosproj = hsSparseData->Projection(2);
  hsSparseData->GetAxis(4)->SetRange(1,hsSparseData->GetAxis(4)->GetNbins()); 

  TH1D *openinganglegammaproj = 0x0;
  TH1D *openinganglepi0proj = 0x0;
  TH1D *openingangleCproj = 0x0;
  TH1D *openingangleBproj = 0x0;
  TH1D *openingangleetaproj = 0x0;
  TH1D *openingangleSplittedTrackssproj = 0x0;
  TH1D *openingangleSplittedTrackosproj = 0x0;
  if(hsSparseMC) {
    hsSparseMC->GetAxis(4)->SetRange(kElectronFromGamma+1,kElectronFromGamma+1);  
    openinganglegammaproj = hsSparseMC->Projection(2);
    hsSparseMC->GetAxis(4)->SetRange(kElectronFromPi0+1,kElectronFromPi0+1);  
    openinganglepi0proj = hsSparseMC->Projection(2);
    hsSparseMC->GetAxis(4)->SetRange(kElectronFromEta+1,kElectronFromEta+1);  
    openingangleetaproj = hsSparseMC->Projection(2);
    hsSparseMC->GetAxis(4)->SetRange(kElectronFromC+1,kElectronFromC+1);  
    openingangleCproj = hsSparseMC->Projection(2);
    hsSparseMC->GetAxis(4)->SetRange(kElectronFromB+1,kElectronFromB+1);  
    openingangleBproj = hsSparseMC->Projection(2);
    hsSparseMC->GetAxis(4)->SetRange(1,hsSparseMC->GetAxis(4)->GetNbins());  
    hsSparseMC->GetAxis(5)->SetRange(kSplittedSs+1,kSplittedSs+1);  
    openingangleSplittedTrackssproj = hsSparseMC->Projection(2);
    hsSparseMC->GetAxis(5)->SetRange(kSplittedOs+1,kSplittedOs+1);  
    openingangleSplittedTrackosproj = hsSparseMC->Projection(2);
    hsSparseMC->GetAxis(5)->SetRange(1,hsSparseMC->GetAxis(5)->GetNbins()); 
  }

  // Projection pt-opening angle
  hsSparseData->GetAxis(4)->SetRange(kPp+1,kPp+1); 
  TH2D *openingangleppproj2D = hsSparseData->Projection(0,2);
  hsSparseData->GetAxis(4)->SetRange(kNn+1,kNn+1); 
  TH2D *openinganglennproj2D = hsSparseData->Projection(0,2);
  hsSparseData->GetAxis(4)->SetRange(kPp+1,kNn+1); 
  TH2D *openinganglessproj2D = hsSparseData->Projection(0,2);
  hsSparseData->GetAxis(4)->SetRange(kR+1,kR+1); 
  TH2D *openinganglerproj2D  = hsSparseData->Projection(0,2);
  hsSparseData->GetAxis(4)->SetRange(kOs+1,kOs+1); 
  TH2D *openingangleosproj2D = hsSparseData->Projection(0,2);
  hsSparseData->GetAxis(4)->SetRange(1,hsSparseData->GetAxis(4)->GetNbins()); 

  TH2D *openinganglegammaproj2D = 0x0;
  TH2D *openinganglepi0proj2D = 0x0;
  TH2D *openingangleCproj2D = 0x0;
  TH2D *openingangleBproj2D = 0x0;
  TH2D *openingangleetaproj2D = 0x0;
  TH2D *openingangleSplittedTrackssproj2D = 0x0;
  TH2D *openingangleSplittedTrackosproj2D = 0x0;
  if(hsSparseMC) {
    hsSparseMC->GetAxis(4)->SetRange(kElectronFromGamma+1,kElectronFromGamma+1); 
    openinganglegammaproj2D = hsSparseMC->Projection(0,2);
    hsSparseMC->GetAxis(4)->SetRange(kElectronFromPi0+1,kElectronFromPi0+1); 
    openinganglepi0proj2D = hsSparseMC->Projection(0,2);
    hsSparseMC->GetAxis(4)->SetRange(kElectronFromEta+1,kElectronFromEta+1); 
    openingangleetaproj2D = hsSparseMC->Projection(0,2);
    hsSparseMC->GetAxis(4)->SetRange(kElectronFromC+1,kElectronFromC+1); 
    openingangleCproj2D = hsSparseMC->Projection(0,2);
    hsSparseMC->GetAxis(4)->SetRange(kElectronFromB+1,kElectronFromB+1); 
    openingangleBproj2D = hsSparseMC->Projection(0,2);
    hsSparseMC->GetAxis(4)->SetRange(1, hsSparseMC->GetAxis(4)->GetNbins());  
    hsSparseMC->GetAxis(5)->SetRange(kSplittedSs+1,kSplittedSs+1); 
    openingangleSplittedTrackssproj2D = hsSparseMC->Projection(0,2);
    hsSparseMC->GetAxis(5)->SetRange(kSplittedOs+1,kSplittedOs+1); 
    openingangleSplittedTrackosproj2D = hsSparseMC->Projection(0,2);
    hsSparseMC->GetAxis(5)->SetRange(1,hsSparseMC->GetAxis(5)->GetNbins());
  }

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

  ////////////////////////////
  // Invariant mass
  ///////////////////////////

  // Cuts on the opening angle
  TAxis *axisOpeningAngleData = hsSparseData->GetAxis(2);
  Int_t binCutData = axisOpeningAngleData->FindBin(fOpeningAngleCut);
  hsSparseData->GetAxis(2)->SetRange(1,binCutData);
  
  // Debug
  //printf("Get Bin low edge %f, Get Bin Up edge %f for hsSparseData\n",axisOpeningAngleData->GetBinLowEdge(binCutData),axisOpeningAngleData->GetBinUpEdge(binCutData));

  // Invariant mass
  hsSparseData->GetAxis(4)->SetRange(kPp+1,kPp+1);  
  TH1D *invmassppproj = hsSparseData->Projection(3);
  hsSparseData->GetAxis(4)->SetRange(kNn+1,kNn+1);  
  TH1D *invmassnnproj = hsSparseData->Projection(3);
  hsSparseData->GetAxis(4)->SetRange(kPp+1,kNn+1);  
  TH1D *invmassssproj = hsSparseData->Projection(3);
  hsSparseData->GetAxis(4)->SetRange(kR+1,kR+1);  
  TH1D *invmassrproj  = hsSparseData->Projection(3);
  hsSparseData->GetAxis(4)->SetRange(kOs+1,kOs+1);  
  TH1D *invmassosproj = hsSparseData->Projection(3);
  hsSparseData->GetAxis(4)->SetRange(1,hsSparseData->GetAxis(4)->GetNbins()); 

  TH1D *invmassgammaproj = 0x0;
  TH1D *invmasspi0proj = 0x0;
  TH1D *invmassCproj = 0x0;
  TH1D *invmassBproj = 0x0;
  TH1D *invmassetaproj = 0x0;
  TH1D *invmassSplittedTrackssproj = 0x0;
  TH1D *invmassSplittedTrackosproj = 0x0;
  if(hsSparseMC) {
    TAxis *axisOpeningAngleMC = hsSparseMC->GetAxis(2);
    Int_t binCutMC = axisOpeningAngleMC->FindBin(fOpeningAngleCut);
    hsSparseMC->GetAxis(2)->SetRange(1,binCutMC);
    
    // Debug
    //printf("Get Bin low edge %f, Get Bin Up edge %f for hsSparseMC\n",axisOpeningAngleMC->GetBinLowEdge(binCutMC),axisOpeningAngleMC->GetBinUpEdge(binCutMC));
    
    hsSparseMC->GetAxis(4)->SetRange(kElectronFromGamma+1,kElectronFromGamma+1);  
    invmassgammaproj = hsSparseMC->Projection(3);
    hsSparseMC->GetAxis(4)->SetRange(kElectronFromPi0+1,kElectronFromPi0+1);  
    invmasspi0proj = hsSparseMC->Projection(3);
    hsSparseMC->GetAxis(4)->SetRange(kElectronFromEta+1,kElectronFromEta+1);  
    invmassetaproj = hsSparseMC->Projection(3);
    hsSparseMC->GetAxis(4)->SetRange(kElectronFromC+1,kElectronFromC+1);  
    invmassCproj = hsSparseMC->Projection(3);
    hsSparseMC->GetAxis(4)->SetRange(kElectronFromB+1,kElectronFromB+1);  
    invmassBproj = hsSparseMC->Projection(3);
    hsSparseMC->GetAxis(4)->SetRange(1,hsSparseMC->GetAxis(4)->GetNbins()); 
    hsSparseMC->GetAxis(5)->SetRange(kSplittedSs+1,kSplittedSs+1);  
    invmassSplittedTrackssproj = hsSparseMC->Projection(3);
    hsSparseMC->GetAxis(5)->SetRange(kSplittedOs+1,kSplittedOs+1);  
    invmassSplittedTrackosproj = hsSparseMC->Projection(3);
    hsSparseMC->GetAxis(5)->SetRange(1,hsSparseMC->GetAxis(5)->GetNbins());
  }
  
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

  // Projection pt-invariant mass angle
  hsSparseData->GetAxis(4)->SetRange(kPp+1,kPp+1); 
  TH2D *invmassppproj2D = hsSparseData->Projection(0,3);
  hsSparseData->GetAxis(4)->SetRange(kNn+1,kNn+1); 
  TH2D *invmassnnproj2D = hsSparseData->Projection(0,3);
  hsSparseData->GetAxis(4)->SetRange(kPp+1,kNn+1); 
  TH2D *invmassssproj2D = hsSparseData->Projection(0,3);
  hsSparseData->GetAxis(4)->SetRange(kR+1,kR+1); 
  TH2D *invmassrproj2D  = hsSparseData->Projection(0,3);
  hsSparseData->GetAxis(4)->SetRange(kOs+1,kOs+1); 
  TH2D *invmassosproj2D = hsSparseData->Projection(0,3);
  hsSparseData->GetAxis(4)->SetRange(1,hsSparseData->GetAxis(4)->GetNbins()); 

  TH2D *invmassgammaproj2D = 0x0;
  TH2D *invmasspi0proj2D = 0x0;
  TH2D *invmassCproj2D = 0x0;
  TH2D *invmassBproj2D = 0x0;
  TH2D *invmassetaproj2D = 0x0;
  TH2D *invmassSplittedTrackssproj2D = 0x0;
  TH2D *invmassSplittedTrackosproj2D = 0x0;
  if(hsSparseMC) {
    hsSparseMC->GetAxis(4)->SetRange(kElectronFromGamma+1,kElectronFromGamma+1); 
    invmassgammaproj2D = hsSparseMC->Projection(0,3);
    hsSparseMC->GetAxis(4)->SetRange(kElectronFromPi0+1,kElectronFromPi0+1); 
    invmasspi0proj2D = hsSparseMC->Projection(0,3);
    hsSparseMC->GetAxis(4)->SetRange(kElectronFromEta+1,kElectronFromEta+1); 
    invmassetaproj2D = hsSparseMC->Projection(0,3);
    hsSparseMC->GetAxis(4)->SetRange(kElectronFromC+1,kElectronFromC+1); 
    invmassCproj2D = hsSparseMC->Projection(0,3);
    hsSparseMC->GetAxis(4)->SetRange(kElectronFromB+1,kElectronFromB+1); 
    invmassBproj2D = hsSparseMC->Projection(0,3);
    hsSparseMC->GetAxis(4)->SetRange(1,hsSparseMC->GetAxis(4)->GetNbins()); 
    hsSparseMC->GetAxis(5)->SetRange(kSplittedSs+1,kSplittedSs+1); 
    invmassSplittedTrackssproj2D = hsSparseMC->Projection(0,3);
    hsSparseMC->GetAxis(5)->SetRange(kSplittedOs+1,kSplittedOs+1); 
    invmassSplittedTrackosproj2D = hsSparseMC->Projection(0,3);
    hsSparseMC->GetAxis(5)->SetRange(1,hsSparseMC->GetAxis(5)->GetNbins()); 
  }

  
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


  /////////////
  // Plot 
  ////////////

  // Draw histograms for opening angle
  TCanvas * copeningangle =new TCanvas("openingangle","Openingangle",800,800);
  copeningangle->cd();
  //openingangleppproj->Draw();
  //openinganglennproj->Draw("same");
  openinganglessproj->Draw();
  //openinganglerproj->Draw("same");
  openingangleosproj->Draw("same");
  if(openinganglegammaproj) openinganglegammaproj->Draw("same");
  if(openinganglepi0proj) openinganglepi0proj->Draw("same");
  //if(openingangleCproj) openingangleCproj->Draw("same");
  //if(openingangleBproj) openingangleBproj->Draw("same");
  if(openingangleetaproj) openingangleetaproj->Draw("same");
  if(openingangleSplittedTrackssproj) openingangleSplittedTrackssproj->Draw("same");
  if(openingangleSplittedTrackosproj) openingangleSplittedTrackosproj->Draw("same");
  TLegend *lego = new TLegend(0.4,0.6,0.89,0.89);
  //lego->AddEntry(openingangleppproj,"positive-positive","p");
  //lego->AddEntry(openinganglennproj,"negative-negative","p");
  lego->AddEntry(openinganglessproj,"same-sign","p");
  //lego->AddEntry(openinganglerproj,"rotated","p");
  lego->AddEntry(openingangleosproj,"positive-negative","p");
  if(openinganglegammaproj) lego->AddEntry(openinganglegammaproj,"e^{+}e^{-} from #gamma","p");
  if(openinganglepi0proj) lego->AddEntry(openinganglepi0proj,"e^{+}e^{-} from #pi^{0}","p");
  //if(openingangleCproj) lego->AddEntry(openingangleCproj,"e^{+}e^{-} from c","p");
  //if(openingangleBproj) lego->AddEntry(openingangleBproj,"e^{+}e^{-} from b","p");
  if(openingangleetaproj) lego->AddEntry(openingangleetaproj,"e^{+}e^{-} from #eta","p");
  if(openingangleSplittedTrackssproj) lego->AddEntry(openingangleSplittedTrackssproj,"Splitted tracks same sign","p");
  if(openingangleSplittedTrackosproj) lego->AddEntry(openingangleSplittedTrackosproj,"Splitted tracks opposite sign","p");
  lego->Draw("same");
  
  // Draw histograms for invariant mass
  TCanvas * cinvmass =new TCanvas("invmass","Invmass",800,800);
  cinvmass->cd();
  //invmassppproj->Draw();
  //invmassnnproj->Draw("same");
  invmassssproj->Draw();
  //invmassrproj->Draw("same");
  invmassosproj->Draw("same");
  if(invmassgammaproj) invmassgammaproj->Draw("same");
  if(invmasspi0proj) invmasspi0proj->Draw("same");
  //if(invmassCproj) invmassCproj->Draw("same");
  //if(invmassBproj) invmassBproj->Draw("same");
  if(invmassetaproj) invmassetaproj->Draw("same");
  if(invmassSplittedTrackssproj) invmassSplittedTrackssproj->Draw("same");
  if(invmassSplittedTrackosproj) invmassSplittedTrackosproj->Draw("same");
  TLegend *legi = new TLegend(0.4,0.6,0.89,0.89);
  //legi->AddEntry(invmassppproj,"positive-positive","p");
  //legi->AddEntry(invmassnnproj,"negative-negative","p");
  legi->AddEntry(invmassssproj,"same-sign","p");
  //legi->AddEntry(invmassrproj,"rotated","p");
  legi->AddEntry(invmassosproj,"positive-negative","p");
  if(invmassgammaproj) legi->AddEntry(invmassgammaproj,"e^{+}e^{-} from #gamma","p");
  if(invmasspi0proj) legi->AddEntry(invmasspi0proj,"e^{+}e^{-} from #pi^{0}","p");
  //if(invmassCproj) legi->AddEntry(invmassCproj,"e^{+}e^{-} from c","p");
  //if(invmassBproj) legi->AddEntry(invmassBproj,"e^{+}e^{-} from b","p");
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

  /////////////////////////////////////
  // Data Radius and chi2Ndf if AliKF
  ////////////////////////////////////
  
  TH1F *hDataRadius = dynamic_cast<TH1F *>(fList->FindObject("DataRadius")); 
  TH1F *hDataChi2Ndf = dynamic_cast<TH1F *>(fList->FindObject("DataChi2Ndf")); 

  if(hDataRadius || hDataChi2Ndf) {
    TCanvas * cDataRadiusChi2Ndf =new TCanvas("CanvasDataRadiusChi2Ndf","CanvasDataRadiusChi2Ndf",800,800);
    cDataRadiusChi2Ndf->Divide(2,1);
    cDataRadiusChi2Ndf->cd(1);
    if(hDataRadius) hDataRadius->Draw();
    cDataRadiusChi2Ndf->cd(2);
    if(hDataChi2Ndf) hDataChi2Ndf->Draw();
  }

  ///////////////////////
  // Data DCA
  //////////////////////
 
  TH1F *hDataDCA = dynamic_cast<TH1F *>(fList->FindObject("DataDCA")); 
  
  if(hDataDCA) {
    TCanvas * cDataDCA =new TCanvas("CanvasDataDCA","CanvasDataDCA",800,800);
    cDataDCA->cd(1);
    hDataDCA->Draw();
  }

  /////////////////////////////////////
  // MC Radius and chi2Ndf if AliKF
  ////////////////////////////////////
  
  TH2F *hMCRadius = dynamic_cast<TH2F *>(fList->FindObject("MCRadius")); 
  TH2F *hMCChi2Ndf = dynamic_cast<TH2F *>(fList->FindObject("MCChi2Ndf")); 

  if(hMCRadius || hMCChi2Ndf) {
    TCanvas * cMCRadiusChi2Ndf =new TCanvas("CanvasMCRadiusChi2Ndf","CanvasMCRadiusChi2Ndf",800,800);
    cMCRadiusChi2Ndf->Divide(2,1);
    cMCRadiusChi2Ndf->cd(1);
    //TH1D *hMCRadiusBackground = hMCRadius->ProjectionX("MCRadiusBackGround",1,1,"e");
    TH1D *hMCRadiusGamma = hMCRadius->ProjectionX("MCRadiusGamma",2,2,"e");
    TH1D *hMCRadiusPi0 = hMCRadius->ProjectionX("MCRadiusPi0",3,3,"e");
    TH1D *hMCRadiusEta = hMCRadius->ProjectionX("MCRadiusEta",4,4,"e");
    TH1D *hMCRadiusC = hMCRadius->ProjectionX("MCRadiusC",5,5,"e");
    TH1D *hMCRadiusB = hMCRadius->ProjectionX("MCRadiusB",6,6,"e");
    //hMCRadiusBackground->Draw();
    hMCRadiusGamma->Draw();
    hMCRadiusPi0->Draw("same");
    hMCRadiusEta->Draw("same");
    hMCRadiusC->Draw("same");
    hMCRadiusB->Draw("same");
    TLegend *legRadius = new TLegend(0.4,0.6,0.89,0.89);
    //legRadius->AddEntry(hMCRadiusBackground,"Background","p");
    legRadius->AddEntry(hMCRadiusGamma,"#gamma","p");
    legRadius->AddEntry(hMCRadiusPi0,"#pi^{0}","p");
    legRadius->AddEntry(hMCRadiusEta,"#eta","p");
    legRadius->AddEntry(hMCRadiusC,"c","p");
    legRadius->AddEntry(hMCRadiusB,"b","p");
    legRadius->Draw("same");
    cMCRadiusChi2Ndf->cd(2);
    //TH1D *hMCChi2NdfBackground = hMCChi2Ndf->ProjectionX("MCChi2NdfBackGround",1,1,"e");
    TH1D *hMCChi2NdfGamma = hMCChi2Ndf->ProjectionX("MCChi2NdfGamma",2,2,"e");
    TH1D *hMCChi2NdfPi0 = hMCChi2Ndf->ProjectionX("MCChi2NdfPi0",3,3,"e");
    TH1D *hMCChi2NdfEta = hMCChi2Ndf->ProjectionX("MCChi2NdfEta",4,4,"e");
    TH1D *hMCChi2NdfC = hMCChi2Ndf->ProjectionX("MCChi2NdfC",5,5,"e");
    TH1D *hMCChi2NdfB = hMCChi2Ndf->ProjectionX("MCChi2NdfB",6,6,"e");
    //hMCChi2NdfBackground->Draw();
    hMCChi2NdfGamma->Draw();
    hMCChi2NdfPi0->Draw("same");
    hMCChi2NdfEta->Draw("same");
    hMCChi2NdfC->Draw("same");
    hMCChi2NdfB->Draw("same");
    TLegend *legChi2Ndf = new TLegend(0.4,0.6,0.89,0.89);
    //legChi2Ndf->AddEntry(hMCChi2NdfBackground,"Background","p");
    legChi2Ndf->AddEntry(hMCChi2NdfGamma,"#gamma","p");
    legChi2Ndf->AddEntry(hMCChi2NdfPi0,"#pi^{0}","p");
    legChi2Ndf->AddEntry(hMCChi2NdfEta,"#eta","p");
    legChi2Ndf->AddEntry(hMCChi2NdfC,"c","p");
    legChi2Ndf->AddEntry(hMCChi2NdfB,"b","p");
    legChi2Ndf->Draw("same");
  }

  ///////////////////////
  // MC DCA
  //////////////////////
  
  TH2F *hMCDCA = dynamic_cast<TH2F *>(fList->FindObject("MCDCA")); 
  
  if(hMCDCA) {
    TCanvas * cMCDCA =new TCanvas("CanvasMCDCA","CanvasMCDCA",800,800);
    cMCDCA->cd(1);
    //TH1D *hMCDCABackground = hMCDCA->ProjectionX("MCDCABackGround",1,1,"e");
    TH1D *hMCDCAGamma = hMCDCA->ProjectionX("MCDCAGamma",2,2,"e");
    TH1D *hMCDCAPi0 = hMCDCA->ProjectionX("MCDCAPi0",3,3,"e");
    TH1D *hMCDCAEta = hMCDCA->ProjectionX("MCDCAEta",4,4,"e");
    TH1D *hMCDCAC = hMCDCA->ProjectionX("MCDCAC",5,5,"e");
    TH1D *hMCDCAB = hMCDCA->ProjectionX("MCDCAB",6,6,"e");
    //hMCDCABackground->Draw();
    hMCDCAGamma->Draw();
    hMCDCAPi0->Draw("same");
    hMCDCAEta->Draw("same");
    hMCDCAC->Draw("same");
    hMCDCAB->Draw("same");
    TLegend *legDCA = new TLegend(0.4,0.6,0.89,0.89);
    //legDCA->AddEntry(hMCDCABackground,"Background","p");
    legDCA->AddEntry(hMCDCAGamma,"#gamma","p");
    legDCA->AddEntry(hMCDCAPi0,"#pi^{0}","p");
    legDCA->AddEntry(hMCDCAEta,"#eta","p");
    legDCA->AddEntry(hMCDCAC,"c","p");
    legDCA->AddEntry(hMCDCAB,"b","p");
    legDCA->Draw("same");
  }


  /////////////////////////
  // PID Partner Signal
  /////////////////////////
  TF1 *betheBlochElectron = 0x0;
  TF1 *betheBlochMuon = 0x0;
  TF1 *betheBlochPion = 0x0;
  TF1 *betheBlochKaon = 0x0;
  TF1 *betheBlochProton = 0x0;

  TH2F *hsignalPidPartner0 = dynamic_cast<TH2F *>(fList->FindObject("TPCPartner0")); 
  TH2F *hsignalPidPartner1 = dynamic_cast<TH2F *>(fList->FindObject("TPCPartner1"));
  if((!hsignalPidPartner0) && (!hsignalPidPartner1)) {
    hsignalPidPartner0 = dynamic_cast<TH2F *>(fList->FindObject("ITSPartner0")); 
    hsignalPidPartner1 = dynamic_cast<TH2F *>(fList->FindObject("ITSPartner1"));
    
    betheBlochElectron = new TF1("betheBlochElectron",BetheBlochElectronITS,0.1,10.0,0);
    betheBlochMuon = new TF1("betheBlochMuon",BetheBlochMuonITS,0.1,10.0,0);
    betheBlochPion = new TF1("betheBlochPion",BetheBlochPionITS,0.1,10.0,0);
    betheBlochKaon = new TF1("betheBlochKaon",BetheBlochKaonITS,0.1,10.0,0);
    betheBlochProton = new TF1("betheBlochProton",BetheBlochProtonITS,0.1,10.0,0);

  }
  else {

    betheBlochElectron = new TF1("betheBlochElectron",BetheBlochElectronTPC,0.1,10.0,0);
    betheBlochMuon = new TF1("betheBlochMuon",BetheBlochMuonTPC,0.1,10.0,0);
    betheBlochPion = new TF1("betheBlochPion",BetheBlochPionTPC,0.1,10.0,0);
    betheBlochKaon = new TF1("betheBlochKaon",BetheBlochKaonTPC,0.1,10.0,0);
    betheBlochProton = new TF1("betheBlochProton",BetheBlochProtonTPC,0.1,10.0,0);

  }


  if((hsignalPidPartner0) || (hsignalPidPartner1)) {
    TCanvas * cPidSignal =new TCanvas("cPidSignal","cPidSignal",800,800);
    cPidSignal->Divide(2,1);
    cPidSignal->cd(1);
    if(hsignalPidPartner0) hsignalPidPartner0->Draw("colz");
    if(betheBlochElectron) betheBlochElectron->Draw("same");
    if(betheBlochMuon) betheBlochMuon->Draw("same");
    if(betheBlochPion) betheBlochPion->Draw("same");
    if(betheBlochKaon) betheBlochKaon->Draw("same");
    if(betheBlochProton) betheBlochProton->Draw("same");
    cPidSignal->cd(2);
    if(hsignalPidPartner1) hsignalPidPartner1->Draw("colz");
    if(betheBlochElectron) betheBlochElectron->Draw("same");
    if(betheBlochMuon) betheBlochMuon->Draw("same");
    if(betheBlochPion) betheBlochPion->Draw("same");
    if(betheBlochKaon) betheBlochKaon->Draw("same");
    if(betheBlochProton) betheBlochProton->Draw("same");
  }
  
  THnSparseF *hsSparseITSsignal = dynamic_cast<THnSparseF *>(fList->FindObject("SparseITSsignal"));
  if(hsSparseITSsignal) {

   
    TH2D *sddsdd = hsSparseITSsignal->Projection(1,2);
    TH2D *ssdssd = hsSparseITSsignal->Projection(3,4);
    TH2D *sddssda = hsSparseITSsignal->Projection(1,3);    
    TH2D *sddssdb = hsSparseITSsignal->Projection(2,4);    
    TH2D *sddssdc = hsSparseITSsignal->Projection(1,4);   
    TH2D *sddssdd = hsSparseITSsignal->Projection(2,3);    

    TCanvas * cITSSignal =new TCanvas("cITSSignal","cITSSignal",800,800);
    cITSSignal->Divide(2,3);
    cITSSignal->cd(1);
    sddsdd->Draw("colz");
    cITSSignal->cd(2);
    ssdssd->Draw("colz");
    cITSSignal->cd(3);
    sddssda->Draw("colz");
    cITSSignal->cd(4);
    sddssdb->Draw("colz");
    cITSSignal->cd(5);
    sddssdc->Draw("colz");
    cITSSignal->cd(6);
    sddssdd->Draw("colz");
  
  } 

  THnSparseF *hsSparseITSsignalSplit = dynamic_cast<THnSparseF *>(fList->FindObject("SparseITSsignalSplit"));
  if(hsSparseITSsignalSplit) {

    // no splitted
    hsSparseITSsignalSplit->GetAxis(0)->SetRange(1,1);
    
    TH1D *layerITS2 = hsSparseITSsignalSplit->Projection(1);
    TH1D *layerITS3 = hsSparseITSsignalSplit->Projection(2);
    TH1D *layerITS4 = hsSparseITSsignalSplit->Projection(3);
    TH1D *layerITS5 = hsSparseITSsignalSplit->Projection(4);

    // splitted
    hsSparseITSsignalSplit->GetAxis(0)->SetRange(2,2);
    
    TH1D *layerITS2s = hsSparseITSsignalSplit->Projection(1);
    TH1D *layerITS3s = hsSparseITSsignalSplit->Projection(2);
    TH1D *layerITS4s = hsSparseITSsignalSplit->Projection(3);
    TH1D *layerITS5s = hsSparseITSsignalSplit->Projection(4);
    
    TCanvas * cITSSignalSplit =new TCanvas("cITSSignalSplit","cITSSignalSplit",800,800);
    cITSSignalSplit->Divide(2,2);
    cITSSignalSplit->cd(1);
    layerITS2->Draw();
    layerITS2s->Draw("same");
    TLegend *legITS2 = new TLegend(0.4,0.6,0.89,0.89);
    legITS2->AddEntry(layerITS2,"No splitted","p");
    legITS2->AddEntry(layerITS2s,"Splitted","p");
    legITS2->Draw("same");
    cITSSignalSplit->cd(2);
    layerITS3->Draw();
    layerITS3s->Draw("same");
    TLegend *legITS3 = new TLegend(0.4,0.6,0.89,0.89);
    legITS3->AddEntry(layerITS3,"No splitted","p");
    legITS3->AddEntry(layerITS3s,"Splitted","p");
    legITS3->Draw("same");
    cITSSignalSplit->cd(3);
    layerITS4->Draw();
    layerITS4s->Draw("same");
    TLegend *legITS4 = new TLegend(0.4,0.6,0.89,0.89);
    legITS4->AddEntry(layerITS4,"No splitted","p");
    legITS4->AddEntry(layerITS4s,"Splitted","p");
    legITS4->Draw("same");
    cITSSignalSplit->cd(4);
    layerITS5->Draw();
    layerITS5s->Draw("same");
    TLegend *legITS5 = new TLegend(0.4,0.6,0.89,0.89);
    legITS5->AddEntry(layerITS5,"No splitted","p");
    legITS5->AddEntry(layerITS5s,"Splitted","p");
    legITS5->Draw("same");

  
  }

  ////////////////////////
  // Cut efficiencies
  ////////////////////////
  
  THnSparseF *hsSparseMCe = dynamic_cast<THnSparseF *>(fList->FindObject("CutPassedMC"));
  if(!hsSparseMCe) return;

  // init histos
  TAxis *axissources = hsSparseMCe->GetAxis(2);
  Int_t  nbsources = axissources->GetNbins(); 
  TAxis *axiscuts = hsSparseMCe->GetAxis(1);
  Int_t  nbcuts = axiscuts->GetNbins(); 
  TH1D **histopassedcuts = new TH1D*[nbsources*nbcuts];  
  Double_t *nbEntriesCuts = new Double_t[nbsources*nbcuts]; 
  for(Int_t k =0; k < nbsources*nbcuts; k++){
    nbEntriesCuts[k] = 0.0;
    histopassedcuts[k] = 0x0;
  }
  
  //printf("Number of cuts %d\n",nbcuts);
  
  // canvas
  TCanvas * chsSparseMCeeff =new TCanvas("hsSparseMCeeffDebug","hsSparseMCeeffDebug",800,800);
  chsSparseMCeeff->Divide(3,1);
  
  // histos
  for(Int_t sourceid = 0; sourceid < nbsources; sourceid++) {
    hsSparseMCe->GetAxis(2)->SetRange(sourceid+1,sourceid+1);  
    for(Int_t cut = 0; cut < nbcuts; cut++){
      hsSparseMCe->GetAxis(1)->SetRange(cut+1,cut+1); 
      histopassedcuts[sourceid*nbcuts+cut] = hsSparseMCe->Projection(0);
      hsSparseMCe->GetAxis(1)->SetRange(1,hsSparseMCe->GetAxis(1)->GetNbins()); 
    }
    hsSparseMCe->GetAxis(2)->SetRange(1,hsSparseMCe->GetAxis(2)->GetNbins());  
  }
  
  // calcul efficiencies
  
  // histos
  for(Int_t sourceid = 0; sourceid < nbsources; sourceid++) {
    // Next is compared to the partner tracked
    for(Int_t cut = 2; cut < nbcuts; cut++){
      nbEntriesCuts[sourceid*nbcuts+cut] = histopassedcuts[sourceid*nbcuts+cut]->GetEntries();
      if(histopassedcuts[sourceid*nbcuts+1]->GetEntries() > 0.0) histopassedcuts[sourceid*nbcuts+cut]->Divide(histopassedcuts[sourceid*nbcuts+1]);
    }
    // First one is if the partner is tracked.
    nbEntriesCuts[sourceid*nbcuts+1] = histopassedcuts[sourceid*nbcuts+1]->GetEntries();
    if(histopassedcuts[sourceid*nbcuts]->GetEntries() > 0.0) histopassedcuts[sourceid*nbcuts+1]->Divide(histopassedcuts[sourceid*nbcuts]);
    // First one is input
    nbEntriesCuts[sourceid*nbcuts] = histopassedcuts[sourceid*nbcuts]->GetEntries();
  }
  
  // ratios
  for(Int_t sourceid = 0; sourceid < nbsources; sourceid++) {
    for(Int_t cut = 1; cut < nbcuts; cut++){
      if(nbEntriesCuts[sourceid*nbcuts] > 0.0) nbEntriesCuts[sourceid*nbcuts+cut] = nbEntriesCuts[sourceid*nbcuts+cut]/nbEntriesCuts[sourceid*nbcuts]; 
    }
  }
  TH1F *ratioHistoEntriesGamma = new TH1F("ratioHistoEntriesGamma","", nbcuts-1, 0.0, nbcuts-1.0);
  TH1F *ratioHistoEntriesPi0 = new TH1F("ratioHistoEntriesPi0","", nbcuts-1, 0.0, nbcuts-1.0);
  TH1F *ratioHistoEntriesC = new TH1F("ratioHistoEntriesC","", nbcuts-1, 0.0, nbcuts-1.0);
  for(Int_t k = 1; k < nbcuts; k++){
    if((nbcuts+k) < (nbsources*nbcuts)) ratioHistoEntriesGamma->SetBinContent(k,nbEntriesCuts[nbcuts+k]);
    if((2*nbcuts+k) < (nbsources*nbcuts)) ratioHistoEntriesPi0->SetBinContent(k,nbEntriesCuts[2*nbcuts+k]);
    if((4*nbcuts+k) < (nbsources*nbcuts)) ratioHistoEntriesC->SetBinContent(k,nbEntriesCuts[4*nbcuts+k]);
  }     
  //
  TAxis *xAxisGamma = ratioHistoEntriesGamma->GetXaxis();
  xAxisGamma->SetBinLabel(1,"Partner tracked");
  xAxisGamma->SetBinLabel(2,"Opposite sign");
  xAxisGamma->SetBinLabel(3,"Single Track Cut");
  xAxisGamma->SetBinLabel(4,"Shared Clusters");
  xAxisGamma->SetBinLabel(5,"PID");
  xAxisGamma->SetBinLabel(6,"DCA");
  xAxisGamma->SetBinLabel(7,"Chi^{2}/Ndf");
  xAxisGamma->SetBinLabel(8,"Opening angle");
  xAxisGamma->SetBinLabel(9,"Invariant mass");
  //
  TAxis *xAxisPi0 = ratioHistoEntriesPi0->GetXaxis();
  xAxisPi0->SetBinLabel(1,"Partner tracked");
  xAxisPi0->SetBinLabel(2,"Opposite sign");
  xAxisPi0->SetBinLabel(3,"Single Track Cut");
  xAxisPi0->SetBinLabel(4,"Shared Clusters");
  xAxisPi0->SetBinLabel(5,"PID");
  xAxisPi0->SetBinLabel(6,"DCA");
  xAxisPi0->SetBinLabel(7,"Chi^{2}/Ndf");
  xAxisPi0->SetBinLabel(8,"Opening angle");
  xAxisPi0->SetBinLabel(9,"Invariant mass");
  //
  TAxis *xAxisC = ratioHistoEntriesC->GetXaxis();
  xAxisC->SetBinLabel(1,"Partner tracked");
  xAxisC->SetBinLabel(2,"Opposite sign");
  xAxisC->SetBinLabel(3,"Single Track Cut");
  xAxisC->SetBinLabel(4,"Shared Clusters");
  xAxisC->SetBinLabel(5,"PID");
  xAxisC->SetBinLabel(6,"DCA");
  xAxisC->SetBinLabel(7,"Chi^{2}/Ndf");
  xAxisC->SetBinLabel(8,"Opening angle");
  xAxisC->SetBinLabel(9,"Invariant mass");
  //
  TCanvas * cRatioHistoEntries =new TCanvas("cRatioHistoEntries","cRatioHistoEntries",800,800);
  cRatioHistoEntries->cd(1);
  ratioHistoEntriesGamma->SetStats(0);
  ratioHistoEntriesGamma->Draw();
  ratioHistoEntriesPi0->SetStats(0);
  ratioHistoEntriesPi0->Draw("same");
  ratioHistoEntriesC->SetStats(0);
  //ratioHistoEntriesC->Draw("same");
  TLegend *legEntries = new TLegend(0.4,0.6,0.89,0.89);
  legEntries->AddEntry(ratioHistoEntriesGamma,"#gamma","l");
  legEntries->AddEntry(ratioHistoEntriesPi0,"#pi^{0}","l");
  //legEntries->AddEntry(ratioHistoEntriesC,"c","p");
  legEntries->Draw("same"); 
  
  // plot Debug
  Int_t source = 1;
  chsSparseMCeeff->cd(1);
  if(((source*nbcuts+0)> (nbsources*nbcuts-1)) || ((source*nbcuts+1)> (nbsources*nbcuts-1)) || ((source*nbcuts+2)> (nbsources*nbcuts-1)) || ((source*nbcuts+3)> (nbsources*nbcuts-1)) || ((source*nbcuts+4)> (nbsources*nbcuts-1)) || ((source*nbcuts+5)> (nbsources*nbcuts-1)) || ((source*nbcuts+6)> (nbsources*nbcuts-1)) || ((source*nbcuts+7)> (nbsources*nbcuts-1)) || ((source*nbcuts+8)> (nbsources*nbcuts-1)) || ((source*nbcuts+9)> (nbsources*nbcuts-1))) {
    delete [] histopassedcuts;
    delete [] nbEntriesCuts;
    return;
  }
  if((!histopassedcuts[source*nbcuts+0]) || (!histopassedcuts[source*nbcuts+1]) || (!histopassedcuts[source*nbcuts+2]) || (!histopassedcuts[source*nbcuts+3]) || (!histopassedcuts[source*nbcuts+4]) || (!histopassedcuts[source*nbcuts+5]) || (!histopassedcuts[source*nbcuts+6]) || (!histopassedcuts[source*nbcuts+7]) || (!histopassedcuts[source*nbcuts+8]) || (!histopassedcuts[source*nbcuts+9])) {
    delete [] histopassedcuts;
    delete [] nbEntriesCuts;
    return;
  }
  histopassedcuts[source*nbcuts+0]->SetTitle("#gamma");
  histopassedcuts[source*nbcuts+1]->SetTitle("#gamma");
  histopassedcuts[source*nbcuts+2]->SetTitle("#gamma");
  histopassedcuts[source*nbcuts+3]->SetTitle("#gamma");
  histopassedcuts[source*nbcuts+4]->SetTitle("#gamma");
  histopassedcuts[source*nbcuts+5]->SetTitle("#gamma");
  histopassedcuts[source*nbcuts+6]->SetTitle("#gamma");
  histopassedcuts[source*nbcuts+7]->SetTitle("#gamma");
  histopassedcuts[source*nbcuts+8]->SetTitle("#gamma");
  histopassedcuts[source*nbcuts+9]->SetTitle("#gamma");
  //histopassedcuts[source*nbcuts+0]->SetStats(0);
  histopassedcuts[source*nbcuts+1]->SetStats(0);
  histopassedcuts[source*nbcuts+2]->SetStats(0);
  histopassedcuts[source*nbcuts+3]->SetStats(0);
  histopassedcuts[source*nbcuts+4]->SetStats(0);
  histopassedcuts[source*nbcuts+5]->SetStats(0);
  histopassedcuts[source*nbcuts+6]->SetStats(0);
  histopassedcuts[source*nbcuts+7]->SetStats(0);
  histopassedcuts[source*nbcuts+8]->SetStats(0);
  histopassedcuts[source*nbcuts+9]->SetStats(0);
  //histopassedcuts[source*nbcuts+0]->Draw();
  //histopassedcuts[source*nbcuts+1]->Draw("");
  histopassedcuts[source*nbcuts+2]->Draw();
  histopassedcuts[source*nbcuts+3]->Draw("same");
  //histopassedcuts[source*nbcuts+4]->Draw("same");
  histopassedcuts[source*nbcuts+5]->Draw("same");
  histopassedcuts[source*nbcuts+6]->Draw("same");
  //histopassedcuts[source*nbcuts+7]->Draw("same");
  histopassedcuts[source*nbcuts+8]->Draw("same");
  histopassedcuts[source*nbcuts+9]->Draw("same");
  TLegend *legb = new TLegend(0.4,0.6,0.89,0.89);
  //legb->AddEntry(histopassedcuts[source*nbcuts+0],"all","p");
  //legb->AddEntry(histopassedcuts[source*nbcuts+1],"Partner tracked","p");
  legb->AddEntry(histopassedcuts[source*nbcuts+2],"Opposite sign","p");
  legb->AddEntry(histopassedcuts[source*nbcuts+3],"SingleTrackPart","p");
  //legb->AddEntry(histopassedcuts[source*nbcuts+4],"SharedCluster","p");
  legb->AddEntry(histopassedcuts[source*nbcuts+5],"PID","p");
  legb->AddEntry(histopassedcuts[source*nbcuts+6],"DCA","p");
  //legb->AddEntry(histopassedcuts[source*nbcuts+7],"Chi2Ndf","p");
  legb->AddEntry(histopassedcuts[source*nbcuts+8],"OpeningAngle","p");
  legb->AddEntry(histopassedcuts[source*nbcuts+9],"InvMass","p");
  legb->Draw("same");
  
  source = 2;
  if(((source*nbcuts+0)> (nbsources*nbcuts-1)) || ((source*nbcuts+1)> (nbsources*nbcuts-1)) || ((source*nbcuts+2)> (nbsources*nbcuts-1)) || ((source*nbcuts+3)> (nbsources*nbcuts-1)) || ((source*nbcuts+4)> (nbsources*nbcuts-1)) || ((source*nbcuts+5)> (nbsources*nbcuts-1)) || ((source*nbcuts+6)> (nbsources*nbcuts-1)) || ((source*nbcuts+7)> (nbsources*nbcuts-1)) || ((source*nbcuts+8)> (nbsources*nbcuts-1)) || ((source*nbcuts+9)> (nbsources*nbcuts-1))) {
    delete [] histopassedcuts;
    delete [] nbEntriesCuts;
    return;
  }
  if((!histopassedcuts[source*nbcuts+0]) || (!histopassedcuts[source*nbcuts+1]) || (!histopassedcuts[source*nbcuts+2]) || (!histopassedcuts[source*nbcuts+3]) || (!histopassedcuts[source*nbcuts+4]) || (!histopassedcuts[source*nbcuts+5]) || (!histopassedcuts[source*nbcuts+6]) || (!histopassedcuts[source*nbcuts+7]) || (!histopassedcuts[source*nbcuts+8]) || (!histopassedcuts[source*nbcuts+9])) {
    delete [] histopassedcuts;
    delete [] nbEntriesCuts;
    return;
  }
  chsSparseMCeeff->cd(2);
  histopassedcuts[source*nbcuts+0]->SetTitle("#pi^{0}");
  histopassedcuts[source*nbcuts+1]->SetTitle("#pi^{0}");
  histopassedcuts[source*nbcuts+2]->SetTitle("#pi^{0}");
  histopassedcuts[source*nbcuts+3]->SetTitle("#pi^{0}");
  histopassedcuts[source*nbcuts+4]->SetTitle("#pi^{0}");
  histopassedcuts[source*nbcuts+5]->SetTitle("#pi^{0}");
  histopassedcuts[source*nbcuts+6]->SetTitle("#pi^{0}");
  histopassedcuts[source*nbcuts+7]->SetTitle("#pi^{0}");
  histopassedcuts[source*nbcuts+8]->SetTitle("#pi^{0}");
  histopassedcuts[source*nbcuts+9]->SetTitle("#pi^{0}");
  //histopassedcuts[source*nbcuts+0]->SetStats(0);
  histopassedcuts[source*nbcuts+1]->SetStats(0);
  histopassedcuts[source*nbcuts+2]->SetStats(0);
  histopassedcuts[source*nbcuts+3]->SetStats(0);
  histopassedcuts[source*nbcuts+4]->SetStats(0);
  histopassedcuts[source*nbcuts+5]->SetStats(0);
  histopassedcuts[source*nbcuts+6]->SetStats(0);
  histopassedcuts[source*nbcuts+7]->SetStats(0);
  histopassedcuts[source*nbcuts+8]->SetStats(0);
  histopassedcuts[source*nbcuts+9]->SetStats(0);
  //histopassedcuts[source*nbcuts+0]->Draw();
  //histopassedcuts[source*nbcuts+1]->Draw();
  histopassedcuts[source*nbcuts+2]->Draw();
  histopassedcuts[source*nbcuts+3]->Draw("same");
  //histopassedcuts[source*nbcuts+4]->Draw("same");
  histopassedcuts[source*nbcuts+5]->Draw("same");
  histopassedcuts[source*nbcuts+6]->Draw("same");
  //histopassedcuts[source*nbcuts+7]->Draw("same");
  histopassedcuts[source*nbcuts+8]->Draw("same");
  histopassedcuts[source*nbcuts+9]->Draw("same");
  TLegend *legc = new TLegend(0.4,0.6,0.89,0.89);
  //legc->AddEntry(histopassedcuts[source*nbcuts+0],"all","p");
  //legc->AddEntry(histopassedcuts[source*nbcuts+1],"Partner tracked","p");
  legc->AddEntry(histopassedcuts[source*nbcuts+2],"Opposite sign","p");
  legc->AddEntry(histopassedcuts[source*nbcuts+3],"SingleTrackPart","p");
  //legc->AddEntry(histopassedcuts[source*nbcuts+4],"SharedCluster","p");
  legc->AddEntry(histopassedcuts[source*nbcuts+5],"PID","p");
  legc->AddEntry(histopassedcuts[source*nbcuts+6],"DCA","p");
  //legc->AddEntry(histopassedcuts[source*nbcuts+7],"Chi2Ndf","p");
  legc->AddEntry(histopassedcuts[source*nbcuts+8],"OpeningAngle","p");
  legc->AddEntry(histopassedcuts[source*nbcuts+9],"InvMass","p");
  legc->Draw("same");
  
  source = 4;
  if(((source*nbcuts+0)> (nbsources*nbcuts-1)) || ((source*nbcuts+1)> (nbsources*nbcuts-1)) || ((source*nbcuts+2)> (nbsources*nbcuts-1)) || ((source*nbcuts+3)> (nbsources*nbcuts-1)) || ((source*nbcuts+4)> (nbsources*nbcuts-1)) || ((source*nbcuts+5)> (nbsources*nbcuts-1)) || ((source*nbcuts+6)> (nbsources*nbcuts-1)) || ((source*nbcuts+7)> (nbsources*nbcuts-1)) || ((source*nbcuts+8)> (nbsources*nbcuts-1)) || ((source*nbcuts+9)> (nbsources*nbcuts-1))) {
    delete [] histopassedcuts;
    delete [] nbEntriesCuts;
    return;
  }
  if((!histopassedcuts[source*nbcuts+0]) || (!histopassedcuts[source*nbcuts+1]) || (!histopassedcuts[source*nbcuts+2]) || (!histopassedcuts[source*nbcuts+3]) || (!histopassedcuts[source*nbcuts+4]) || (!histopassedcuts[source*nbcuts+5]) || (!histopassedcuts[source*nbcuts+6]) || (!histopassedcuts[source*nbcuts+7]) || (!histopassedcuts[source*nbcuts+8]) || (!histopassedcuts[source*nbcuts+9])) {
    delete [] histopassedcuts;
    delete [] nbEntriesCuts;
    return;
  }
  chsSparseMCeeff->cd(3);
  histopassedcuts[source*nbcuts+0]->SetTitle("C");
  histopassedcuts[source*nbcuts+1]->SetTitle("C");
  histopassedcuts[source*nbcuts+2]->SetTitle("C");
  histopassedcuts[source*nbcuts+3]->SetTitle("C");
  histopassedcuts[source*nbcuts+4]->SetTitle("C");
  histopassedcuts[source*nbcuts+5]->SetTitle("C");
  histopassedcuts[source*nbcuts+6]->SetTitle("C");
  histopassedcuts[source*nbcuts+7]->SetTitle("C");
  histopassedcuts[source*nbcuts+8]->SetTitle("C");
  histopassedcuts[source*nbcuts+9]->SetTitle("C");
  //histopassedcuts[source*nbcuts+0]->SetStats(0);
  histopassedcuts[source*nbcuts+1]->SetStats(0);
  histopassedcuts[source*nbcuts+2]->SetStats(0);
  histopassedcuts[source*nbcuts+3]->SetStats(0);
  histopassedcuts[source*nbcuts+4]->SetStats(0);
  histopassedcuts[source*nbcuts+5]->SetStats(0);
  histopassedcuts[source*nbcuts+6]->SetStats(0);
  histopassedcuts[source*nbcuts+7]->SetStats(0);
  histopassedcuts[source*nbcuts+8]->SetStats(0);
  histopassedcuts[source*nbcuts+9]->SetStats(0);
  //histopassedcuts[source*nbcuts+0]->Draw();
  //histopassedcuts[source*nbcuts+1]->Draw();
  histopassedcuts[source*nbcuts+2]->Draw();
  histopassedcuts[source*nbcuts+3]->Draw("same");
  //histopassedcuts[source*nbcuts+4]->Draw("same");
  histopassedcuts[source*nbcuts+5]->Draw("same");
  histopassedcuts[source*nbcuts+6]->Draw("same");
  //histopassedcuts[source*nbcuts+7]->Draw("same");
  histopassedcuts[source*nbcuts+8]->Draw("same");
  histopassedcuts[source*nbcuts+9]->Draw("same");
  TLegend *lege = new TLegend(0.4,0.6,0.89,0.89);
  //lege->AddEntry(histopassedcuts[source*nbcuts+0],"all","p");
  //lege->AddEntry(histopassedcuts[source*nbcuts+1],"Partner tracked","p");
  lege->AddEntry(histopassedcuts[source*nbcuts+2],"Opposite sign","p");
  lege->AddEntry(histopassedcuts[source*nbcuts+3],"SingleTrackPart","p");
  //lege->AddEntry(histopassedcuts[source*nbcuts+4],"SharedCluster","p");
  lege->AddEntry(histopassedcuts[source*nbcuts+5],"PID","p");
  lege->AddEntry(histopassedcuts[source*nbcuts+6],"DCA","p");
  //lege->AddEntry(histopassedcuts[source*nbcuts+7],"Chi2Ndf","p");
  lege->AddEntry(histopassedcuts[source*nbcuts+8],"OpeningAngle","p");
  lege->AddEntry(histopassedcuts[source*nbcuts+9],"InvMass","p");
  lege->Draw("same");
  
  //////////////////////
  // Input
  //////////////////////
  
  TCanvas * chsSparseMCein =new TCanvas("hsSparseMCeinput","hsSparseMCeinput",800,800);
  chsSparseMCein->cd(1);
  Double_t nbGamma = 0.0;
  source = 1;
  if(((source*nbcuts+0)> (nbsources*nbcuts-1)) || (!histopassedcuts[source*nbcuts+0])) {
    delete [] histopassedcuts;
    delete [] nbEntriesCuts;
    return;
  }
  nbGamma = histopassedcuts[source*nbcuts+0]->GetEntries();
  histopassedcuts[source*nbcuts+0]->SetStats(0);
  histopassedcuts[source*nbcuts+0]->Draw();
  TLegend *leginput = new TLegend(0.4,0.6,0.89,0.89);
  leginput->AddEntry(histopassedcuts[source*nbcuts+0],"#gamma","p");
  Double_t nbPi0 = 0.0;
  source = 2;
  if(((source*nbcuts+0)> (nbsources*nbcuts-1)) || (!histopassedcuts[source*nbcuts+0])) {
    delete [] histopassedcuts;
    delete [] nbEntriesCuts;
    return;
  }
  nbPi0 = histopassedcuts[source*nbcuts+0]->GetEntries();
  histopassedcuts[source*nbcuts+0]->SetStats(0);
  histopassedcuts[source*nbcuts+0]->Draw("same");
  leginput->AddEntry(histopassedcuts[source*nbcuts+0],"#pi^{0}","p");
  Double_t nbEta = 0.0;
  source = 3;
  if(((source*nbcuts+0)> (nbsources*nbcuts-1)) || (!histopassedcuts[source*nbcuts+0])) {
    delete [] histopassedcuts;
    delete [] nbEntriesCuts;
    return;
  }
  nbEta = histopassedcuts[source*nbcuts+0]->GetEntries();
  histopassedcuts[source*nbcuts+0]->SetStats(0);
  histopassedcuts[source*nbcuts+0]->Draw("same");
  leginput->AddEntry(histopassedcuts[source*nbcuts+0],"#eta","p");
  Double_t nbC = 0.0;
  source = 4;
  if(((source*nbcuts+0)> (nbsources*nbcuts-1)) || (!histopassedcuts[source*nbcuts+0])) {
    delete [] histopassedcuts;
    delete [] nbEntriesCuts;
    return;
  }
  nbC = histopassedcuts[source*nbcuts+0]->GetEntries();
  histopassedcuts[source*nbcuts+0]->SetStats(0);
  histopassedcuts[source*nbcuts+0]->Draw("same");
  leginput->AddEntry(histopassedcuts[source*nbcuts+0],"c","p");
  leginput->Draw("same");
  
  AliInfo(Form("Gamma %f, pi^{0} %f and #eta %f, c %f",nbGamma,nbPi0,nbEta,nbC));
  
  //////////////////////
  // Tracked
  //////////////////////
  
  TCanvas * cTracked = new TCanvas("cTracked","cTracked",800,800);
  cTracked->cd(1);
  source = 1;
  if(((source*nbcuts+1)> (nbsources*nbcuts-1)) || (!histopassedcuts[source*nbcuts+1])) {
    delete [] histopassedcuts;
    delete [] nbEntriesCuts;
    return;
  }
  histopassedcuts[source*nbcuts+1]->Draw();
  TLegend *legTracked = new TLegend(0.4,0.6,0.89,0.89);
  legTracked->AddEntry(histopassedcuts[source*nbcuts+1],"#gamma","p");
  source = 2;
  if(((source*nbcuts+1)> (nbsources*nbcuts-1)) || (!histopassedcuts[source*nbcuts+1])) {
    delete [] histopassedcuts;
    delete [] nbEntriesCuts;
    return;
  }
  histopassedcuts[source*nbcuts+1]->Draw("same");
  legTracked->AddEntry(histopassedcuts[source*nbcuts+1],"#pi^{0}","p");
  legTracked->Draw("same");
  
  delete [] histopassedcuts;
  delete [] nbEntriesCuts;
  
}
//_____________________________________________________________________________
Double_t AliHFEelecbackground::BetheBlochElectronITS(const Double_t *x, const Double_t * /*par*/) 
{
  //
  // Bethe Bloch for ITS
  //
  static AliITSPIDResponse itsPidResponse;
  return itsPidResponse.Bethe(x[0],AliPID::ParticleMass(0));
}
//_____________________________________________________________________________
Double_t AliHFEelecbackground::BetheBlochMuonITS(const Double_t *x, const Double_t * /*par*/) 
{
  //
  // Bethe Bloch for ITS
  //
  static AliITSPIDResponse itsPidResponse;
  return itsPidResponse.Bethe(x[0],AliPID::ParticleMass(1));
}
//_____________________________________________________________________________
Double_t AliHFEelecbackground::BetheBlochPionITS(const Double_t *x, const Double_t * /*par*/) 
{
  //
  // Bethe Bloch for ITS
  //
  static AliITSPIDResponse itsPidResponse;
  return itsPidResponse.Bethe(x[0],AliPID::ParticleMass(2));
}
//_____________________________________________________________________________
Double_t AliHFEelecbackground::BetheBlochKaonITS(const Double_t *x, const Double_t * /*par*/) 
{
  //
  // Bethe Bloch for ITS
  //
  static AliITSPIDResponse itsPidResponse;
  return itsPidResponse.Bethe(x[0],AliPID::ParticleMass(3));
}
//_____________________________________________________________________________
Double_t AliHFEelecbackground::BetheBlochProtonITS(const Double_t *x, const Double_t * /*par*/) 
{
  //
  // Bethe Bloch for ITS
  //
  static AliITSPIDResponse itsPidResponse;
  return itsPidResponse.Bethe(x[0],AliPID::ParticleMass(4));
}
//_____________________________________________________________________________
Double_t AliHFEelecbackground::BetheBlochElectronTPC(const Double_t *x, const Double_t * /*par*/) 
{
  //
  // Bethe Bloch for TPC
  //
  static AliTPCPIDResponse tpcPidResponse;
  return tpcPidResponse.GetExpectedSignal(x[0],AliPID::kElectron);
}
//_____________________________________________________________________________
Double_t AliHFEelecbackground::BetheBlochMuonTPC(const Double_t *x, const Double_t * /*par*/) 
{
  //
  // Bethe Bloch for TPC
  //
  static AliTPCPIDResponse tpcPidResponse;
  return tpcPidResponse.GetExpectedSignal(x[0],AliPID::kMuon);
}
//_____________________________________________________________________________
Double_t AliHFEelecbackground::BetheBlochPionTPC(const Double_t *x, const Double_t * /*par*/) 
{
  //
  // Bethe Bloch for TPC
  //
  static AliTPCPIDResponse tpcPidResponse;
  return tpcPidResponse.GetExpectedSignal(x[0],AliPID::kPion);
}
//_____________________________________________________________________________
Double_t AliHFEelecbackground::BetheBlochKaonTPC(const Double_t *x, const Double_t * /*par*/) 
{
  //
  // Bethe Bloch for TPC
  //
  static AliTPCPIDResponse tpcPidResponse;
  return tpcPidResponse.GetExpectedSignal(x[0],AliPID::kKaon);
}
//_____________________________________________________________________________
Double_t AliHFEelecbackground::BetheBlochProtonTPC(const Double_t *x, const Double_t * /*par*/) 
{
  //
  // Bethe Bloch for TPC
  //
  static AliTPCPIDResponse tpcPidResponse;
  return tpcPidResponse.GetExpectedSignal(x[0],AliPID::kProton);
}



