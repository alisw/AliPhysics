/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
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

/* $Id: AliAnalysisTaskSELambdac.cxx 63644 2013-07-23 09:04:11Z zconesa $ */

/////////////////////////////////////////////////////////////
//
// AliAnalysisTaskSE for the extraction of signal(e.g Lambdac) of heavy flavor
// decay candidates with the MC truth.
// Authors: r.romita@gsi.de
/////////////////////////////////////////////////////////////

#include <TClonesArray.h>
#include <TNtuple.h>
#include <TCanvas.h>
#include <TList.h>
#include <TString.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TDatabasePDG.h>

#include <AliAnalysisDataSlot.h>
#include <AliAnalysisDataContainer.h>
#include "AliAnalysisManager.h"
#include "AliAODHandler.h"
#include "AliAODEvent.h"
#include "AliAODVertex.h"
#include "AliAODTrack.h"
#include "AliAODMCHeader.h"
#include "AliAODMCParticle.h"
#include "AliAODRecoDecayHF3Prong.h"
#include "AliAnalysisVertexingHF.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisTaskSELambdac.h"
#include "AliKFParticle.h"
#include "AliAODPidHF.h"
#include "AliRDHFCutsLctopKpi.h"
#include "AliRDHFCuts.h"
#include "AliKFVertex.h"
#include "AliESDVertex.h"
//#include "AliAODpidUtil.h"
#include "AliAODPid.h"
#include "AliInputEventHandler.h"
#include "AliPID.h"
#include "AliNormalizationCounter.h"
#include "AliVertexingHFUtils.h"


/// \cond CLASSIMP
ClassImp(AliAnalysisTaskSELambdac);
/// \endcond

//________________________________________________________________________
AliAnalysisTaskSELambdac::AliAnalysisTaskSELambdac():
AliAnalysisTaskSE(),
  fOutput(0), 
  fHistNEvents(0),
  fhChi2(0),
  fhMassPtGreater3(0),
  fhMassPtGreater3TC(0),
  fhMassPtGreater3Kp(0),
  fhMassPtGreater3KpTC(0),
  fhMassPtGreater3Lpi(0),
  fhMassPtGreater3LpiTC(0),
  fhMassPtGreater3Dk(0),
  fhMassPtGreater3DkTC(0),
  fhMassPtGreater33Pr(0),
  fhMassPtGreater33PrTC(0),
  fhMassPtGreater2(0),
  fhMassPtGreater2TC(0),
  fhMassPtGreater2Kp(0),
  fhMassPtGreater2KpTC(0),
  fhMassPtGreater2Lpi(0),
  fhMassPtGreater2LpiTC(0),
  fhMassPtGreater2Dk(0),
  fhMassPtGreater2DkTC(0),
  fhMassPtGreater23Pr(0),
  fhMassPtGreater23PrTC(0),
  fhMassLcPt(0),
  fhMassLcplusPt(0),
  fhMassLcminusPt(0),
  fhEta3Prong(0),
  fhEta3ProngAcc(0),
  fhEta3ProngProd(0),
  fhEta3ProngAn(0),
  fhRap3Prong(0),
  fhRap3ProngAcc(0),
  fhRap3ProngProd(0),
  fhRap3ProngAn(0),
  fhSelectBit(0),
  fhProtonPtProngLcPt(0),
  fhBProtonPtProngLcPt(0),
  fhProtond0ProngLcPt(0),
  fhBProtond0ProngLcPt(0),
  fhKaonPtProngLcPt(0),
  fhBKaonPtProngLcPt(0),
  fhKaond0ProngLcPt(0),
  fhBKaond0ProngLcPt(0),
  fhPionPtProngLcPt(0),
  fhBPionPtProngLcPt(0),
  fhPiond0ProngLcPt(0),
  fhBPiond0ProngLcPt(0),
  fhDist12PrimLcPt(0),
  fhBDist12PrimLcPt(0),
  fhSigmaVertLcPt(0),
  fhBSigmaVertLcPt(0),
  fhdcasLcPt(0),
  fhBdcasLcPt(0),
  fhCosPointingAngleLcPt(0),
  fhBCosPointingAngleLcPt(0),
  fhDecayLengthLcPt(0),
  fhBDecayLengthLcPt(0),
  fhSum2LcPt(0),
  fhBSum2LcPt(0),
  fhPtMaxLcPt(0),
  fhBPtMaxLcPt(0),
  fNtupleLambdac(0),
  fUpmasslimit(2.486),
  fLowmasslimit(2.086),
  fNPtBins(0),
  fRDCutsAnalysis(0),
  fRDCutsProduction(0),
  fListCuts(0),
  fFillNtuple(kFALSE),
  fReadMC(kFALSE),
  fMCPid(kFALSE),
  fRealPid(kFALSE),
  fResPid(kTRUE),
  fUseKF(kFALSE),
  fAnalysis(kFALSE),
  fVHF(0),
  fFillVarHists(kFALSE),
  fMultiplicityHists(kFALSE),
  fPriorsHists(kFALSE),
  fLcCut(kFALSE),
  fLcPIDCut(kFALSE),    
  fNentries(0),
  fOutputMC(0),
  fAPriori(0),
  fMultiplicity(0),
//fUtilPid(0),
  fPIDResponse(0),
  fCounter(0)
{
  // Default constructor
  Float_t ptlims[7]={0.,2.,4.,6.,8.,12.,24.};
  SetPtBinLimit(7,ptlims);
   for(Int_t icut=0; icut<2; icut++) fCutsKF[icut]=0.;
   for(Int_t j=0; j<3*kMaxPtBins; j++){
    fMassHist[j]=0x0;
    fMassHistTC[j]=0x0;
    fMassHistLpi[j]=0x0;
    fMassHistLpiTC[j]=0x0;
    fMassHistKp[j]=0x0;
    fMassHistKpTC[j]=0x0;
    fMassHistDk[j]=0x0;
    fMassHistDkTC[j]=0x0;
    fMassHist3Pr[j]=0x0;
    fMassHist3PrTC[j]=0x0;
   }
}

//________________________________________________________________________
AliAnalysisTaskSELambdac::AliAnalysisTaskSELambdac(const char *name,Bool_t fillNtuple,AliRDHFCutsLctopKpi *lccutsana,AliRDHFCutsLctopKpi *lccutsprod):
  AliAnalysisTaskSE(name),
  fOutput(0),
  fHistNEvents(0),
  fhChi2(0),
  fhMassPtGreater3(0),
  fhMassPtGreater3TC(0),
  fhMassPtGreater3Kp(0),
  fhMassPtGreater3KpTC(0),
  fhMassPtGreater3Lpi(0),
  fhMassPtGreater3LpiTC(0),
  fhMassPtGreater3Dk(0),
  fhMassPtGreater3DkTC(0),
  fhMassPtGreater33Pr(0),
  fhMassPtGreater33PrTC(0),
  fhMassPtGreater2(0),
  fhMassPtGreater2TC(0),
  fhMassPtGreater2Kp(0),
  fhMassPtGreater2KpTC(0),
  fhMassPtGreater2Lpi(0),
  fhMassPtGreater2LpiTC(0),
  fhMassPtGreater2Dk(0),
  fhMassPtGreater2DkTC(0),
  fhMassPtGreater23Pr(0),
  fhMassPtGreater23PrTC(0),
  fhMassLcPt(0),
  fhMassLcplusPt(0),
  fhMassLcminusPt(0),
  fhEta3Prong(0),
  fhEta3ProngAcc(0),
  fhEta3ProngProd(0),
  fhEta3ProngAn(0),
  fhRap3Prong(0),
  fhRap3ProngAcc(0),
  fhRap3ProngProd(0),
  fhRap3ProngAn(0),
  fhSelectBit(0),
  fhProtonPtProngLcPt(0),
  fhBProtonPtProngLcPt(0),
  fhProtond0ProngLcPt(0),
  fhBProtond0ProngLcPt(0),
  fhKaonPtProngLcPt(0),
  fhBKaonPtProngLcPt(0),
  fhKaond0ProngLcPt(0),
  fhBKaond0ProngLcPt(0),
  fhPionPtProngLcPt(0),
  fhBPionPtProngLcPt(0),
  fhPiond0ProngLcPt(0),
  fhBPiond0ProngLcPt(0),
  fhDist12PrimLcPt(0),
  fhBDist12PrimLcPt(0),
  fhSigmaVertLcPt(0),
  fhBSigmaVertLcPt(0),
  fhdcasLcPt(0),
  fhBdcasLcPt(0),
  fhCosPointingAngleLcPt(0),
  fhBCosPointingAngleLcPt(0),
  fhDecayLengthLcPt(0),
  fhBDecayLengthLcPt(0),
  fhSum2LcPt(0),
  fhBSum2LcPt(0),
  fhPtMaxLcPt(0),
  fhBPtMaxLcPt(0),
  fNtupleLambdac(0),
  fUpmasslimit(2.486),
  fLowmasslimit(2.086),
  fNPtBins(0),
  fRDCutsAnalysis(lccutsana),
  fRDCutsProduction(lccutsprod),
  fListCuts(0),
  fFillNtuple(fillNtuple),
  fReadMC(kFALSE),
  fMCPid(kFALSE),
  fRealPid(kTRUE),
  fResPid(kFALSE),
  fUseKF(kFALSE),
  fAnalysis(kFALSE),
  fVHF(0),
  fFillVarHists(kFALSE),
  fMultiplicityHists(kFALSE),
  fPriorsHists(kFALSE),
  fLcCut(kFALSE),
  fLcPIDCut(kFALSE),
  fNentries(0),
  fOutputMC(0),
  fAPriori(0),
  fMultiplicity(0),
  //fUtilPid(0),
  fPIDResponse(0),
  fCounter(0)
{
  SetPtBinLimit(fRDCutsAnalysis->GetNPtBins()+1,
		fRDCutsAnalysis->GetPtBinLimits());
   for(Int_t icut=0; icut<2; icut++) fCutsKF[icut]=0.;
   for(Int_t j=0; j<3*kMaxPtBins; j++){
    fMassHist[j]=0x0;
    fMassHistTC[j]=0x0;
    fMassHistLpi[j]=0x0;
    fMassHistLpiTC[j]=0x0;
    fMassHistKp[j]=0x0;
    fMassHistKpTC[j]=0x0;
    fMassHistDk[j]=0x0;
    fMassHistDkTC[j]=0x0;
    fMassHist3Pr[j]=0x0;
    fMassHist3PrTC[j]=0x0;
   }
  // Default constructor
  // Output slot #1 writes into a TList container
  DefineOutput(1,TList::Class());  //My private output
  DefineOutput(2,TList::Class());
  DefineOutput(3,TList::Class());
  DefineOutput(4,TH1F::Class());
  DefineOutput(5,TList::Class());
  DefineOutput(6,TList::Class());
  DefineOutput(7,AliNormalizationCounter::Class());
  if (fFillNtuple) {
    // Output slot #2 writes into a TNtuple container
    DefineOutput(8,TNtuple::Class());  //My private output
  }
}

//________________________________________________________________________
AliAnalysisTaskSELambdac::~AliAnalysisTaskSELambdac()
{
  // Destructor
  if (fOutput) {
    delete fOutput;
    fOutput = 0;
  }
  if (fOutputMC) {
    delete fOutputMC;
    fOutputMC = 0;
  }
  if (fAPriori) {
    delete fAPriori;
    fAPriori = 0;
  }
  if (fMultiplicity) {
    delete fMultiplicity;
    fMultiplicity = 0;
  }
  
  if (fVHF) {
    delete fVHF;
    fVHF = 0;
  }
  
  if(fRDCutsAnalysis){
    delete fRDCutsAnalysis;
    fRDCutsAnalysis = 0;
  }
  if(fRDCutsProduction){
    delete fRDCutsProduction;
    fRDCutsProduction = 0;
  }

  if (fListCuts) {
    delete fListCuts;
    fListCuts = 0;
  }
  if (fNentries){
    delete fNentries;
    fNentries = 0;
  }
  /*
    if (fUtilPid){
    delete fUtilPid;
    fUtilPid = 0;
    }
  */
  if (fPIDResponse) {
    delete  fPIDResponse;
  }
  if(fCounter){
   delete fCounter;
   fCounter = 0;
  }

}  
//_________________________________________________________________
void  AliAnalysisTaskSELambdac::SetMassLimits(Float_t range){
  fUpmasslimit = 2.286+range;
  fLowmasslimit = 2.286-range;
}
//_________________________________________________________________
void  AliAnalysisTaskSELambdac::SetMassLimits(Float_t lowlimit, Float_t uplimit){
  if(uplimit>lowlimit)
    {
      fUpmasslimit = lowlimit;
      fLowmasslimit = uplimit;
    }
}


//________________________________________________________________________
void AliAnalysisTaskSELambdac::SetPtBinLimit(Int_t n, Float_t* lim){
  // define pt bins for analysis
  if(n>kMaxPtBins){
    printf("Max. number of Pt bins = %d\n",kMaxPtBins);
    fNPtBins=kMaxPtBins;
    fArrayBinLimits[0]=0.;
    fArrayBinLimits[1]=2.;
    fArrayBinLimits[2]=3.;
    fArrayBinLimits[3]=4.;
    for(Int_t i=4; i<kMaxPtBins+1; i++) fArrayBinLimits[i]=99999999.;
  }else{
    fNPtBins=n-1;
    fArrayBinLimits[0]=lim[0];
    for(Int_t i=1; i<fNPtBins+1; i++) 
      if(lim[i]>fArrayBinLimits[i-1]){
	fArrayBinLimits[i]=lim[i];
      }
      else {
	fArrayBinLimits[i]=fArrayBinLimits[i-1];
      }
    for(Int_t i=fNPtBins; i<kMaxPtBins+1; i++) fArrayBinLimits[i]=99999999.;
  }
  if(fDebug > 1){
    printf("Number of Pt bins = %d\n",fNPtBins);
    for(Int_t i=0; i<fNPtBins; i++) printf(" Bin%d = %8.2f-%8.2f\n",i,fArrayBinLimits[i],fArrayBinLimits[i+1]);    
  }
}
//_________________________________________________________________
Double_t  AliAnalysisTaskSELambdac::GetPtBinLimit(Int_t ibin) const{
  if(ibin>fNPtBins)return -1;
  return fArrayBinLimits[ibin];
} 

//_________________________________________________________________
void AliAnalysisTaskSELambdac::Init()
{
  // Initialization

  if (fDebug > 1) printf("AnalysisTaskSELambdac::Init() \n");

  fListCuts=new TList();
  fListCuts->SetOwner();

  fListCuts->Add(new AliRDHFCutsLctopKpi(*fRDCutsAnalysis));
  fListCuts->Add(new AliRDHFCutsLctopKpi(*fRDCutsProduction));
  PostData(2,fListCuts);
  return;
}

//________________________________________________________________________
void AliAnalysisTaskSELambdac::UserCreateOutputObjects()
{
  // Create the output container
  //
  if (fDebug > 1) printf("AnalysisTaskSELambdac::UserCreateOutputObjects() \n");

  // Several histograms are more conveniently managed in a TList
  fOutput = new TList();
  fOutput->SetOwner();
  fOutput->SetName("OutputHistos");

  TString hisname;
  Int_t index=0;
  //  Int_t indexLS=0;
  for(Int_t i=0;i<fNPtBins;i++){

    index=GetHistoIndex(i);
    //    indexLS=GetLSHistoIndex(i);

    hisname.Form("hMassPt%d",i);
    fMassHist[index]=new TH1F(hisname.Data(),hisname.Data(),200,fLowmasslimit,fUpmasslimit);
    fMassHist[index]->Sumw2();
    hisname.Form("hMassPt%dTC",i);
    fMassHistTC[index]=new TH1F(hisname.Data(),hisname.Data(),200,fLowmasslimit,fUpmasslimit);
    fMassHistTC[index]->Sumw2();

    hisname.Form("hMassPtLpi%d",i);
    fMassHistLpi[index]=new TH1F(hisname.Data(),hisname.Data(),200,fLowmasslimit,fUpmasslimit);
    fMassHistLpi[index]->Sumw2();
    hisname.Form("hMassPtLpi%dTC",i);
    fMassHistLpiTC[index]=new TH1F(hisname.Data(),hisname.Data(),200,fLowmasslimit,fUpmasslimit);
    fMassHistLpiTC[index]->Sumw2();

    hisname.Form("hMassPtKp%d",i);
    fMassHistKp[index]=new TH1F(hisname.Data(),hisname.Data(),200,fLowmasslimit,fUpmasslimit);
    fMassHistKp[index]->Sumw2();
    hisname.Form("hMassPtKp%dTC",i);
    fMassHistKpTC[index]=new TH1F(hisname.Data(),hisname.Data(),200,fLowmasslimit,fUpmasslimit);
    fMassHistKpTC[index]->Sumw2();
    hisname.Form("hMassPtDk%d",i);
    fMassHistDk[index]=new TH1F(hisname.Data(),hisname.Data(),200,fLowmasslimit,fUpmasslimit);
    fMassHistDk[index]->Sumw2();
    hisname.Form("hMassPtDk%dTC",i);
    fMassHistDkTC[index]=new TH1F(hisname.Data(),hisname.Data(),200,fLowmasslimit,fUpmasslimit);
    fMassHistDkTC[index]->Sumw2();

    hisname.Form("hMassPt3Pr%d",i);
    fMassHist3Pr[index]=new TH1F(hisname.Data(),hisname.Data(),200,fLowmasslimit,fUpmasslimit);
    fMassHist3Pr[index]->Sumw2();
    hisname.Form("hMassPt3Pr%dTC",i);
    fMassHist3PrTC[index]=new TH1F(hisname.Data(),hisname.Data(),200,fLowmasslimit,fUpmasslimit);
    fMassHist3PrTC[index]->Sumw2();
   //signal
    index=GetSignalHistoIndex(i);    
    hisname.Form("hSigPt%d",i);
    fMassHist[index]=new TH1F(hisname.Data(),hisname.Data(),200,fLowmasslimit,fUpmasslimit);
    fMassHist[index]->Sumw2();
    hisname.Form("hSigPt%dTC",i);
    fMassHistTC[index]=new TH1F(hisname.Data(),hisname.Data(),200,fLowmasslimit,fUpmasslimit);
    fMassHistTC[index]->Sumw2();
    hisname.Form("hSigPtLpi%d",i);
    fMassHistLpi[index]=new TH1F(hisname.Data(),hisname.Data(),200,fLowmasslimit,fUpmasslimit);
    fMassHistLpi[index]->Sumw2();
    hisname.Form("hSigPtLpi%dTC",i);
    fMassHistLpiTC[index]=new TH1F(hisname.Data(),hisname.Data(),200,fLowmasslimit,fUpmasslimit);
    fMassHistLpiTC[index]->Sumw2();

    hisname.Form("hSigPtKp%d",i);
    fMassHistKp[index]=new TH1F(hisname.Data(),hisname.Data(),200,fLowmasslimit,fUpmasslimit);
    fMassHistKp[index]->Sumw2();
    hisname.Form("hSigPtKp%dTC",i);
    fMassHistKpTC[index]=new TH1F(hisname.Data(),hisname.Data(),200,fLowmasslimit,fUpmasslimit);
    fMassHistKpTC[index]->Sumw2();

    hisname.Form("hSigPtDk%d",i);
    fMassHistDk[index]=new TH1F(hisname.Data(),hisname.Data(),200,fLowmasslimit,fUpmasslimit);
    fMassHistDk[index]->Sumw2();
    hisname.Form("hSigPtDk%dTC",i);
    fMassHistDkTC[index]=new TH1F(hisname.Data(),hisname.Data(),200,fLowmasslimit,fUpmasslimit);
    fMassHistDkTC[index]->Sumw2();

    hisname.Form("hSigPt3Pr%d",i);
    fMassHist3Pr[index]=new TH1F(hisname.Data(),hisname.Data(),200,fLowmasslimit,fUpmasslimit);
    fMassHist3Pr[index]->Sumw2();
    hisname.Form("hSigPt3Pr%dTC",i);
    fMassHist3PrTC[index]=new TH1F(hisname.Data(),hisname.Data(),200,fLowmasslimit,fUpmasslimit);
    fMassHist3PrTC[index]->Sumw2();

    index=GetLbHistoIndex(i); 
    hisname.Form("hSigLbPt%d",i);
    fMassHist[index]=new TH1F(hisname.Data(),hisname.Data(),200,fLowmasslimit,fUpmasslimit);
    fMassHist[index]->Sumw2();
    hisname.Form("hSigLbPt%dTC",i);
    fMassHistTC[index]=new TH1F(hisname.Data(),hisname.Data(),200,fLowmasslimit,fUpmasslimit);
    fMassHistTC[index]->Sumw2();

    index=GetcOnlyHistoIndex(i); 
    hisname.Form("hSigcOnlyPt%d",i);
    fMassHist[index]=new TH1F(hisname.Data(),hisname.Data(),200,fLowmasslimit,fUpmasslimit);
    fMassHist[index]->Sumw2();
    hisname.Form("hSigcOnlyPt%dTC",i);
    fMassHistTC[index]=new TH1F(hisname.Data(),hisname.Data(),200,fLowmasslimit,fUpmasslimit);
    fMassHistTC[index]->Sumw2();
   
    index=GetNoQuarkHistoIndex(i); 
    hisname.Form("hSigNoQuarkPt%d",i);
    fMassHist[index]=new TH1F(hisname.Data(),hisname.Data(),200,fLowmasslimit,fUpmasslimit);
    fMassHist[index]->Sumw2();
    hisname.Form("hSigNoQuarkPt%dTC",i);
    fMassHistTC[index]=new TH1F(hisname.Data(),hisname.Data(),200,fLowmasslimit,fUpmasslimit);
    fMassHistTC[index]->Sumw2();

    index=GetBackgroundHistoIndex(i); 
    hisname.Form("hBkgPt%d",i);
    fMassHist[index]=new TH1F(hisname.Data(),hisname.Data(),200,fLowmasslimit,fUpmasslimit);
    fMassHist[index]->Sumw2();
    hisname.Form("hBkgPt%dTC",i);
    fMassHistTC[index]=new TH1F(hisname.Data(),hisname.Data(),200,fLowmasslimit,fUpmasslimit);
    fMassHistTC[index]->Sumw2();
    hisname.Form("hBkgPtLpi%d",i);
    fMassHistLpi[index]=new TH1F(hisname.Data(),hisname.Data(),200,fLowmasslimit,fUpmasslimit);
    fMassHistLpi[index]->Sumw2();
    hisname.Form("hBkgPtLpi%dTC",i);
    fMassHistLpiTC[index]=new TH1F(hisname.Data(),hisname.Data(),200,fLowmasslimit,fUpmasslimit);
    fMassHistLpiTC[index]->Sumw2();

    hisname.Form("hBkgPtKp%d",i);
    fMassHistKp[index]=new TH1F(hisname.Data(),hisname.Data(),200,fLowmasslimit,fUpmasslimit);
    fMassHistKp[index]->Sumw2();
    hisname.Form("hBkgPtKp%dTC",i);
    fMassHistKpTC[index]=new TH1F(hisname.Data(),hisname.Data(),200,fLowmasslimit,fUpmasslimit);
    fMassHistKpTC[index]->Sumw2();

    hisname.Form("hBkgPtDk%d",i);
    fMassHistDk[index]=new TH1F(hisname.Data(),hisname.Data(),200,fLowmasslimit,fUpmasslimit);
    fMassHistDk[index]->Sumw2();
    hisname.Form("hBkgPtDk%dTC",i);
    fMassHistDkTC[index]=new TH1F(hisname.Data(),hisname.Data(),200,fLowmasslimit,fUpmasslimit);
    fMassHistDkTC[index]->Sumw2();

    hisname.Form("hBkgPt3Pr%d",i);
    fMassHist3Pr[index]=new TH1F(hisname.Data(),hisname.Data(),200,fLowmasslimit,fUpmasslimit);
    fMassHist3Pr[index]->Sumw2();
    hisname.Form("hBkgPt3Pr%dTC",i);
    fMassHist3PrTC[index]=new TH1F(hisname.Data(),hisname.Data(),200,fLowmasslimit,fUpmasslimit);
    fMassHist3PrTC[index]->Sumw2();
  }
  
   for(Int_t ii=0; ii<6*fNPtBins; ii++){
    fOutput->Add(fMassHist[ii]);
    fOutput->Add(fMassHistTC[ii]);
  }
  for(Int_t i=0; i<3*fNPtBins; i++){
    fOutput->Add(fMassHistLpi[i]);
    fOutput->Add(fMassHistLpiTC[i]);
    fOutput->Add(fMassHistKp[i]);
    fOutput->Add(fMassHistKpTC[i]);
    fOutput->Add(fMassHistDk[i]);
    fOutput->Add(fMassHistDkTC[i]);
    fOutput->Add(fMassHist3Pr[i]);
    fOutput->Add(fMassHist3PrTC[i]);
  }
  
  fhMassLcPt        = new TH2F("hMassLcPt","hMassLcPt;3-Prong p_{T} GeV/c;3-Prong Mass GeV/c^2",150,0.,15.,200,fLowmasslimit,fUpmasslimit);
  fhMassLcplusPt    = new TH2F("hMassLcplusPt","hMassLcplusPt;3-Prong p_{T} GeV/c;3-Prong Mass GeV/c^2",150,0.,15.,200,fLowmasslimit,fUpmasslimit);
  fhMassLcminusPt   = new TH2F("hMassLcminusPt","hMassLcminusPt;3-Prong p_{T} GeV/c;3-Prong Mass GeV/c^2",150,0.,15.,200,fLowmasslimit,fUpmasslimit);
  fOutput->Add(fhMassLcPt);  
  fOutput->Add(fhMassLcplusPt);  
  fOutput->Add(fhMassLcminusPt);  
  
  fhEta3Prong		= new TH2F("hEta3Prong","hEta3Prong;3-Prong p_{T} GeV/c;3-Prong #eta",75,0.,15.,50,-2.0,2.0);
  fhEta3ProngAcc	= new TH2F("hEta3ProngAcc","hEta3ProngAcc;3-Prong p_{T} GeV/c;3-Prong #eta",75,0.,15.,50,-2.0,2.0);
  fhEta3ProngProd	= new TH2F("hEta3ProngProd","hEta3ProngProd;3-Prong p_{T} GeV/c;3-Prong #eta",75,0.,15.,50,-2.0,2.0);
  fhEta3ProngAn	    = new TH2F("hEta3ProngAn","hEta3ProngAn;3-Prong p_{T} GeV/c;3-Prong #eta",75,0.,15.,50,-2.0,2.0);
  fhRap3Prong		= new TH2F("hRap3Prong","hRap3Prong;3-Prong p_{T} GeV/c;3-Prong y",75,0.,15.,50,-2.0,2.0);
  fhRap3ProngAcc	= new TH2F("hRap3ProngAcc","hRap3ProngAcc;3-Prong p_{T} GeV/c;3-Prong y",75,0.,15.,50,-2.0,2.0);
  fhRap3ProngProd	= new TH2F("hRap3ProngProd","hRap3ProngProd;3-Prong p_{T} GeV/c;3-Prong y",75,0.,15.,50,-2.0,2.0);
  fhRap3ProngAn	    = new TH2F("hRap3ProngAn","hRap3ProngAn;3-Prong p_{T} GeV/c;3-Prong y",75,0.,15.,50,-2.0,2.0);
    
  fhSelectBit			= new TH1F("hSelectBit","hSelectBit",5,-0.5,5.5);
  fhSelectBit->GetXaxis()->SetBinLabel(2,"All");
  fhSelectBit->GetXaxis()->SetBinLabel(3,"SelectionMap");
  fhSelectBit->GetXaxis()->SetBinLabel(4,"!LcCut");
  fhSelectBit->GetXaxis()->SetBinLabel(5,"!LcPID");

  fOutput->Add(fhEta3Prong);
  fOutput->Add(fhEta3ProngAcc);
  fOutput->Add(fhEta3ProngProd);
  fOutput->Add(fhEta3ProngAn);
  fOutput->Add(fhRap3Prong);
  fOutput->Add(fhRap3ProngAcc);
  fOutput->Add(fhRap3ProngProd);
  fOutput->Add(fhRap3ProngAn);
  fOutput->Add(fhSelectBit);
  
  fHistNEvents = new TH1F("fHistNEvents", "Number of processed events; ; Events",3,-1.5,1.5);
  fHistNEvents->Sumw2();
  fHistNEvents->SetMinimum(0);
  fOutput->Add(fHistNEvents);

  fhChi2 = new TH1F("fhChi2", "Chi2",100,0.,10.);
  fhChi2->Sumw2();
  fOutput->Add(fhChi2);

  fhMassPtGreater3=new TH1F("fhMassPtGreater3","Pt > 3 GeV/c",200,fLowmasslimit,fUpmasslimit);
  fhMassPtGreater3->Sumw2();
  fOutput->Add(fhMassPtGreater3);
  fhMassPtGreater3TC=new TH1F("fhMassPtGreater3TC","Pt > 3 GeV/c",200,fLowmasslimit,fUpmasslimit);
  fhMassPtGreater3TC->Sumw2();
  fOutput->Add(fhMassPtGreater3TC);
  fhMassPtGreater3Kp=new TH1F("fhMassPtGreater3Kp","Pt > 3 GeV/c",200,fLowmasslimit,fUpmasslimit);
  fhMassPtGreater3Kp->Sumw2();
  fOutput->Add(fhMassPtGreater3Kp);
  fhMassPtGreater3KpTC=new TH1F("fhMassPtGreater3KpTC","Pt > 3 GeV/c",200,fLowmasslimit,fUpmasslimit);
  fhMassPtGreater3KpTC->Sumw2();
  fOutput->Add(fhMassPtGreater3KpTC);
  fhMassPtGreater3Lpi=new TH1F("fhMassPtGreater3Lpi","Pt > 3 GeV/c",200,fLowmasslimit,fUpmasslimit);
  fhMassPtGreater3Lpi->Sumw2();
  fOutput->Add(fhMassPtGreater3Lpi);
  fhMassPtGreater3LpiTC=new TH1F("fhMassPtGreater3LpiTC","Pt > 3 GeV/c",200,fLowmasslimit,fUpmasslimit);
  fhMassPtGreater3LpiTC->Sumw2();
  fOutput->Add(fhMassPtGreater3LpiTC);
  fhMassPtGreater3Dk=new TH1F("fhMassPtGreater3Dk","Pt > 3 GeV/c",200,fLowmasslimit,fUpmasslimit);
  fhMassPtGreater3Dk->Sumw2();
  fOutput->Add(fhMassPtGreater3Dk);
  fhMassPtGreater3DkTC=new TH1F("fhMassPtGreater3DkTC","Pt > 3 GeV/c",200,fLowmasslimit,fUpmasslimit);
  fhMassPtGreater3DkTC->Sumw2();
  fOutput->Add(fhMassPtGreater3DkTC);

  fhMassPtGreater33Pr=new TH1F("fhMassPtGreater33Pr","Pt > 3 GeV/c",200,fLowmasslimit,fUpmasslimit);
  fhMassPtGreater33Pr->Sumw2();
  fOutput->Add(fhMassPtGreater33Pr);
  fhMassPtGreater33PrTC=new TH1F("fhMassPtGreater33PrTC","Pt > 3 GeV/c",200,fLowmasslimit,fUpmasslimit);
  fhMassPtGreater33PrTC->Sumw2();
  fOutput->Add(fhMassPtGreater33PrTC);


  fhMassPtGreater2=new TH1F("fhMassPtGreater2","Pt > 2 GeV/c",200,fLowmasslimit,fUpmasslimit);
  fhMassPtGreater2->Sumw2();
  fOutput->Add(fhMassPtGreater2);
  fhMassPtGreater2TC=new TH1F("fhMassPtGreater2TC","Pt > 2 GeV/c",200,fLowmasslimit,fUpmasslimit);
  fhMassPtGreater2TC->Sumw2();
  fOutput->Add(fhMassPtGreater2TC);
  fhMassPtGreater2Kp=new TH1F("fhMassPtGreater2Kp","Pt > 2 GeV/c",200,fLowmasslimit,fUpmasslimit);
  fhMassPtGreater2Kp->Sumw2();
  fOutput->Add(fhMassPtGreater2Kp);
  fhMassPtGreater2KpTC=new TH1F("fhMassPtGreater2KpTC","Pt > 2 GeV/c",200,fLowmasslimit,fUpmasslimit);
  fhMassPtGreater2KpTC->Sumw2();
  fOutput->Add(fhMassPtGreater2KpTC);
  fhMassPtGreater2Lpi=new TH1F("fhMassPtGreater2Lpi","Pt > 2 GeV/c",200,fLowmasslimit,fUpmasslimit);
  fhMassPtGreater2Lpi->Sumw2();
  fOutput->Add(fhMassPtGreater2Lpi);
  fhMassPtGreater2LpiTC=new TH1F("fhMassPtGreater2LpiTC","Pt > 2 GeV/c",200,fLowmasslimit,fUpmasslimit);
  fhMassPtGreater2LpiTC->Sumw2();
  fOutput->Add(fhMassPtGreater2LpiTC);
  fhMassPtGreater2Dk=new TH1F("fhMassPtGreater2Dk","Pt > 2 GeV/c",200,fLowmasslimit,fUpmasslimit);
  fhMassPtGreater2Dk->Sumw2();
  fOutput->Add(fhMassPtGreater2Dk);
  fhMassPtGreater2DkTC=new TH1F("fhMassPtGreater2DkTC","Pt > 2 GeV/c",200,fLowmasslimit,fUpmasslimit);
  fhMassPtGreater2DkTC->Sumw2();
  fOutput->Add(fhMassPtGreater2DkTC);
  fhMassPtGreater23Pr=new TH1F("fhMassPtGreater23Pr","Pt > 2 GeV/c",200,fLowmasslimit,fUpmasslimit);
  fhMassPtGreater23Pr->Sumw2();
  fOutput->Add(fhMassPtGreater23Pr);
  fhMassPtGreater23PrTC=new TH1F("fhMassPtGreater23PrTC","Pt > 2 GeV/c",200,fLowmasslimit,fUpmasslimit);
  fhMassPtGreater23PrTC->Sumw2();
  fOutput->Add(fhMassPtGreater23PrTC);
  
  fOutputMC = new TList();
  fOutputMC->SetOwner();
  fOutputMC->SetName("QAMCHistos");

  //  const char* nameoutput=GetOutputSlot(4)->GetContainer()->GetName();

  fNentries=new TH1F("fNentries", "Integral(1,2) = number of AODs *** Integral(2,3) = number of candidates selected with cuts *** Integral(3,4) = number of Lc selected with cuts *** Integral(4,5) = events with good vertex ***  Integral(5,6) = pt out of bounds", 12,-0.5,14.5);

  //ROS: qui il bin assignment e' modellato su D0 ma sicuramente iv arie
  fNentries->GetXaxis()->SetBinLabel(1,"nEventsAnal");
  fNentries->GetXaxis()->SetBinLabel(2,"nCandSel(Cuts)");
  fNentries->GetXaxis()->SetBinLabel(3,"nLcSelected");
  fNentries->GetXaxis()->SetBinLabel(4,"nEventsGoodVtxS");
  fNentries->GetXaxis()->SetBinLabel(5,"ptbin = -1");
  fNentries->GetXaxis()->SetBinLabel(6,"no daughter");
  fNentries->GetXaxis()->SetBinLabel(7,"nCandSel(Tr)");
  fNentries->GetXaxis()->SetBinLabel(8,"PID=0");
  fNentries->GetXaxis()->SetBinLabel(9,"PID=1");
  fNentries->GetXaxis()->SetBinLabel(10,"PID=2");
  fNentries->GetXaxis()->SetBinLabel(11,"PID=3");
  fNentries->GetXaxis()->SetBinLabel(15,"no. of not on-the-fly rec Lc");
  fNentries->GetXaxis()->SetNdivisions(1,kFALSE);

  hisname.Form("hMass");
  TH1F *hMassInv=new TH1F(hisname.Data(),hisname.Data(),200,fLowmasslimit,fUpmasslimit);
  fOutputMC->Add(hMassInv);
  hisname.Form("hbMass");
  TH1F *hBMassInv=new TH1F(hisname.Data(),hisname.Data(),200,fLowmasslimit,fUpmasslimit);
  fOutputMC->Add(hBMassInv);

  // proton specific
  hisname.Form("hpTOFSignal");
  TH1F *hProtonTOFSignal=new TH1F(hisname.Data(),hisname.Data(),100,12000.,50000.0);
  fOutputMC->Add(hProtonTOFSignal);
  hisname.Form("hbpTOFSignal");
  TH1F *hBProtonTOFSignal=new TH1F(hisname.Data(),hisname.Data(),100,12000.,50000.0);
  fOutputMC->Add(hBProtonTOFSignal);

  hisname.Form("hpTPCSignal");
  TH1F *hProtonTPCSignal=new TH1F(hisname.Data(),hisname.Data(),150,0.,150.0);
  fOutputMC->Add(hProtonTPCSignal);
  hisname.Form("hbpTPCSignal");
  TH1F *hBProtonTPCSignal=new TH1F(hisname.Data(),hisname.Data(),150,0.,150.0);
  fOutputMC->Add(hBProtonTPCSignal);
  
  hisname.Form("hpptProng");
  TH1F *hProtonPtProng=new TH1F(hisname.Data(),hisname.Data(),100,0.,5.0);
  fOutputMC->Add(hProtonPtProng);
  hisname.Form("hbpptProng");
  TH1F *hBProtonPtProng=new TH1F(hisname.Data(),hisname.Data(),100,0.,5.0);
  fOutputMC->Add(hBProtonPtProng);

  hisname.Form("hpRealTot");
  TH1F *hProtonRealTot=new TH1F(hisname.Data(),hisname.Data(),100,0.,10.0);
  fOutputMC->Add(hProtonRealTot);
  hisname.Form("hbpRealTot");
  TH1F *hBProtonRealTot=new TH1F(hisname.Data(),hisname.Data(),100,0.,10.0);
  fOutputMC->Add(hBProtonRealTot);
  hisname.Form("hpIDTot");
  TH1F *hProtonIDTot=new TH1F(hisname.Data(),hisname.Data(),100,0.,10.0);
  fOutputMC->Add(hProtonIDTot);
  hisname.Form("hpIDGood");
  TH1F *hProtonIDGood=new TH1F(hisname.Data(),hisname.Data(),100,0.,10.0);
  fOutputMC->Add(hProtonIDGood);
  hisname.Form("hbpIDGood");
  TH1F *hBProtonIDGood=new TH1F(hisname.Data(),hisname.Data(),100,0.,10.0);
  fOutputMC->Add(hBProtonIDGood);
  hisname.Form("hbpIDTot");
  TH1F *hBProtonIDTot=new TH1F(hisname.Data(),hisname.Data(),100,0.,10.0);
  fOutputMC->Add(hBProtonIDTot);
  hisname.Form("hnopIDp");
  TH1F *hnoProtonIDpTot=new TH1F(hisname.Data(),hisname.Data(),100,0.,10.0);
  fOutputMC->Add(hnoProtonIDpTot);
  hisname.Form("hbnopIDp");
  TH1F *hBnoProtonIDpTot=new TH1F(hisname.Data(),hisname.Data(),100,0.,10.0);
  fOutputMC->Add(hBnoProtonIDpTot);

  hisname.Form("hpd0Prong");
  TH1F *hProtond0Prong=new TH1F(hisname.Data(),hisname.Data(),100,-0.1,0.1);
  fOutputMC->Add(hProtond0Prong);
  hisname.Form("hbpd0Prong");
  TH1F *hBProtond0Prong=new TH1F(hisname.Data(),hisname.Data(),100,-0.1,0.1);
  fOutputMC->Add(hBProtond0Prong);
  hisname.Form("hbpSignalVspTOF");
  TH2F *hBpSignalVspTOF=new TH2F(hisname.Data(),hisname.Data(),100,0.,5.0,100,12000.,50000.0);
  fOutputMC->Add(hBpSignalVspTOF);
  hisname.Form("hbpSignalVspTPC");
  TH2F *hBpSignalVspTPC=new TH2F(hisname.Data(),hisname.Data(),100,0.,5.0,150,0.,150.0);
  fOutputMC->Add(hBpSignalVspTPC);
  hisname.Form("hpSignalVspTOF");
  TH2F *hpSignalVspTOF=new TH2F(hisname.Data(),hisname.Data(),100,0.,5.0,100,12000.,50000.0);
  fOutputMC->Add(hpSignalVspTOF);
  hisname.Form("hpSignalVspTPC");
  TH2F *hpSignalVspTPC=new TH2F(hisname.Data(),hisname.Data(),100,0.,5.0,150,0.,150.0);
  fOutputMC->Add(hpSignalVspTPC);

  hisname.Form("hpSigmaVspTOF");
  TH2F *hpSigmaVspTOF=new TH2F(hisname.Data(),hisname.Data(),100,0.,5.0,100,-10.0,10.0);
  fOutputMC->Add(hpSigmaVspTOF);
  hisname.Form("hpSigmaVspTPC");
  TH2F *hpSigmaVspTPC=new TH2F(hisname.Data(),hisname.Data(),100,0.,5.0,100,-10.0,10.0);
  fOutputMC->Add(hpSigmaVspTPC);
  hisname.Form("hbpSigmaVspTOF");
  TH2F *hBpSigmaVspTOF=new TH2F(hisname.Data(),hisname.Data(),100,0.,5.0,100,-10.0,10.0);
  fOutputMC->Add(hBpSigmaVspTOF);
  hisname.Form("hbpSigmaVspTPC");
  TH2F *hBpSigmaVspTPC=new TH2F(hisname.Data(),hisname.Data(),100,0.,5.0,100,-10.0,10.0);
  fOutputMC->Add(hBpSigmaVspTPC);


  //kaon specific
  hisname.Form("hKTOFSignal");
  TH1F *hKaonTOFSignal=new TH1F(hisname.Data(),hisname.Data(),100,12000.,50000.0);
  fOutputMC->Add(hKaonTOFSignal);
  hisname.Form("hbKTOFSignal");
  TH1F *hBKaonTOFSignal=new TH1F(hisname.Data(),hisname.Data(),100,12000.,50000.0);
  fOutputMC->Add(hBKaonTOFSignal);
  hisname.Form("hKTPCSignal");
  TH1F *hKaonTPCSignal=new TH1F(hisname.Data(),hisname.Data(),150,0.,150.0);
  fOutputMC->Add(hKaonTPCSignal);
  hisname.Form("hbKTPCSignal");
  TH1F *hBKaonTPCSignal=new TH1F(hisname.Data(),hisname.Data(),150,0.,150.0);
  fOutputMC->Add(hBKaonTPCSignal);

  hisname.Form("hKptProng");
  TH1F *hKaonPtProng=new TH1F(hisname.Data(),hisname.Data(),100,0.,5.0);
  fOutputMC->Add(hKaonPtProng);
  hisname.Form("hbKptProng");
  TH1F *hBKaonPtProng=new TH1F(hisname.Data(),hisname.Data(),100,0.,5.0);
  fOutputMC->Add(hBKaonPtProng);
  hisname.Form("hKRealTot");
  TH1F *hKaonRealTot=new TH1F(hisname.Data(),hisname.Data(),100,0.,10.0);
  fOutputMC->Add(hKaonRealTot);
  hisname.Form("hbKRealTot");
  TH1F *hBKaonRealTot=new TH1F(hisname.Data(),hisname.Data(),100,0.,10.0);
  fOutputMC->Add(hBKaonRealTot);
  hisname.Form("hKIDGood");
  TH1F *hKaonIDGood=new TH1F(hisname.Data(),hisname.Data(),100,0.,10.0);
  fOutputMC->Add(hKaonIDGood);
  hisname.Form("hKIDTot");
  TH1F *hKaonIDTot=new TH1F(hisname.Data(),hisname.Data(),100,0.,10.0);
  fOutputMC->Add(hKaonIDTot);
  hisname.Form("hbKIDGood");
  TH1F *hBKaonIDGood=new TH1F(hisname.Data(),hisname.Data(),100,0.,10.0);
  fOutputMC->Add(hBKaonIDGood);
  hisname.Form("hbKIDTot");
  TH1F *hBKaonIDTot=new TH1F(hisname.Data(),hisname.Data(),100,0.,10.0);
  fOutputMC->Add(hBKaonIDTot);
  hisname.Form("hnokIDk");
  TH1F *hnoKaonIDkTot=new TH1F(hisname.Data(),hisname.Data(),100,0.,10.0);
  fOutputMC->Add(hnoKaonIDkTot);
  hisname.Form("hbnokIDk");
  TH1F *hBnoKaonIDkTot=new TH1F(hisname.Data(),hisname.Data(),100,0.,10.0);
  fOutputMC->Add(hBnoKaonIDkTot);


  hisname.Form("hKd0Prong");
  TH1F *hKaond0Prong=new TH1F(hisname.Data(),hisname.Data(),100,-0.1,0.1);
  fOutputMC->Add(hKaond0Prong);
  hisname.Form("hbKd0Prong");
  TH1F *hBKaond0Prong=new TH1F(hisname.Data(),hisname.Data(),100,-0.1,0.1);
  fOutputMC->Add(hBKaond0Prong);
  hisname.Form("hbKSignalVspTOF");
  TH2F *hbKSignalVspTOF=new TH2F(hisname.Data(),hisname.Data(),100,0.,5.0,100,12000.,50000.0);
  fOutputMC->Add(hbKSignalVspTOF);
  hisname.Form("hbKSignalVspTPC");
  TH2F *hbKSignalVspTPC=new TH2F(hisname.Data(),hisname.Data(),100,0.,5.0,150,0.,150.0);
  fOutputMC->Add(hbKSignalVspTPC);
  hisname.Form("hKSignalVspTOF");
  TH2F *hKSignalVspTOF=new TH2F(hisname.Data(),hisname.Data(),100,0.,5.0,100,12000.,50000.0);
  fOutputMC->Add(hKSignalVspTOF);
  hisname.Form("hKSignalVspTPC");
  TH2F *hKSignalVspTPC=new TH2F(hisname.Data(),hisname.Data(),100,0.,5.0,150,0.,150.0);
  fOutputMC->Add(hKSignalVspTPC);
  hisname.Form("hKSigmaVspTOF");
  TH2F *hKSigmaVspTOF=new TH2F(hisname.Data(),hisname.Data(),100,0.,5.0,100,-10.0,10.0);
  fOutputMC->Add(hKSigmaVspTOF);
  hisname.Form("hKSigmaVspTPC");
  TH2F *hKSigmaVspTPC=new TH2F(hisname.Data(),hisname.Data(),100,0.,5.0,100,-10.0,10.0);
  fOutputMC->Add(hKSigmaVspTPC);
  hisname.Form("hbKSigmaVspTOF");
  TH2F *hBKSigmaVspTOF=new TH2F(hisname.Data(),hisname.Data(),100,0.,5.0,100,-10.0,10.0);
  fOutputMC->Add(hBKSigmaVspTOF);
  hisname.Form("hbKSigmaVspTPC");
  TH2F *hBKSigmaVspTPC=new TH2F(hisname.Data(),hisname.Data(),100,0.,5.0,100,-10.0,10.0);
  fOutputMC->Add(hBKSigmaVspTPC);


  // pion specific
  hisname.Form("hpiTOFSignal");
  TH1F *hPionTOFSignal=new TH1F(hisname.Data(),hisname.Data(),100,12000.,50000.0);
  fOutputMC->Add(hPionTOFSignal);
  hisname.Form("hbpiTOFSignal");
  TH1F *hBPionTOFSignal=new TH1F(hisname.Data(),hisname.Data(),100,12000.,50000.0);
  fOutputMC->Add(hBPionTOFSignal);
  hisname.Form("hpiTPCSignal");
  TH1F *hPionTPCSignal=new TH1F(hisname.Data(),hisname.Data(),100,30.,100.0);
  fOutputMC->Add(hPionTPCSignal);
  hisname.Form("hbpiTPCSignal");
  TH1F *hBPionTPCSignal=new TH1F(hisname.Data(),hisname.Data(),100,30.,100.0);
  fOutputMC->Add(hBPionTPCSignal);
  
  hisname.Form("hpiptProng");
  TH1F *hPionPtProng=new TH1F(hisname.Data(),hisname.Data(),100,0.,5.0);
  fOutputMC->Add(hPionPtProng);
  hisname.Form("hbpiptProng");
  TH1F *hBPionPtProng=new TH1F(hisname.Data(),hisname.Data(),100,0.,5.0);
  fOutputMC->Add(hBPionPtProng);

  hisname.Form("hpiRealTot");
  TH1F *hPionRealTot=new TH1F(hisname.Data(),hisname.Data(),100,0.,10.0);
  fOutputMC->Add(hPionRealTot);
  hisname.Form("hbpiRealTot");
  TH1F *hBPionRealTot=new TH1F(hisname.Data(),hisname.Data(),100,0.,10.0);
  fOutputMC->Add(hBPionRealTot);
  hisname.Form("hpiIDGood");
  TH1F *hPionIDGood=new TH1F(hisname.Data(),hisname.Data(),100,0.,10.0);
  fOutputMC->Add(hPionIDGood);
  hisname.Form("hpiIDTot");
  TH1F *hPionIDTot=new TH1F(hisname.Data(),hisname.Data(),100,0.,10.0);
  fOutputMC->Add(hPionIDTot);
  hisname.Form("hbpiIDTot");
  TH1F *hBPionIDTot=new TH1F(hisname.Data(),hisname.Data(),100,0.,10.0);
  fOutputMC->Add(hBPionIDTot);
  hisname.Form("hbpiIDGood");
  TH1F *hBPionIDGood=new TH1F(hisname.Data(),hisname.Data(),100,0.,10.0);
  fOutputMC->Add(hBPionIDGood);
  hisname.Form("hnopiIDpi");
  TH1F *hnoPionIDpiTot=new TH1F(hisname.Data(),hisname.Data(),100,0.,10.0);
  fOutputMC->Add(hnoPionIDpiTot);
  hisname.Form("hbnopiIDpi");
  TH1F *hBnoPionIDpiTot=new TH1F(hisname.Data(),hisname.Data(),100,0.,10.0);
  fOutputMC->Add(hBnoPionIDpiTot);

  hisname.Form("hpid0Prong");
  TH1F *hPiond0Prong=new TH1F(hisname.Data(),hisname.Data(),100,-0.1,-0.1);
  fOutputMC->Add(hPiond0Prong);
  hisname.Form("hbpid0Prong");
  TH1F *hBPiond0Prong=new TH1F(hisname.Data(),hisname.Data(),100,-0.1,0.1);
  fOutputMC->Add(hBPiond0Prong);

  hisname.Form("hpiSignalVspTOF");
  TH2F *hpiSignalVspTOF=new TH2F(hisname.Data(),hisname.Data(),100,0.,5.0,100,12000.,50000.0);
  fOutputMC->Add(hpiSignalVspTOF);
  hisname.Form("hpiSignalVspTPC");
  TH2F *hpiSignalVspTPC=new TH2F(hisname.Data(),hisname.Data(),100,0.,5.0,150,0.,150.0);
  fOutputMC->Add(hpiSignalVspTPC);
  hisname.Form("hbpiSignalVspTOF");
  TH2F *hbpiSignalVspTOF=new TH2F(hisname.Data(),hisname.Data(),100,0.,5.0,100,12000.,50000.0);
  fOutputMC->Add(hbpiSignalVspTOF);
  hisname.Form("hbpiSignalVspTPC");
  TH2F *hbpiSignalVspTPC=new TH2F(hisname.Data(),hisname.Data(),100,0.,5.0,150,0.,150.0);
  fOutputMC->Add(hbpiSignalVspTPC);
  hisname.Form("hpiSigmaVspTOF");
  TH2F *hpiSigmaVspTOF=new TH2F(hisname.Data(),hisname.Data(),100,0.,5.0,100,-10.0,10.0);
  fOutputMC->Add(hpiSigmaVspTOF);
  hisname.Form("hpiSigmaVspTPC");
  TH2F *hpiSigmaVspTPC=new TH2F(hisname.Data(),hisname.Data(),100,0.,5.0,100,-10.0,10.0);
  fOutputMC->Add(hpiSigmaVspTPC);
  hisname.Form("hbpiSigmaVspTOF");
  TH2F *hBpiSigmaVspTOF=new TH2F(hisname.Data(),hisname.Data(),100,0.,5.0,100,-10.0,10.0);
  fOutputMC->Add(hBpiSigmaVspTOF);
  hisname.Form("hbpiSigmaVspTPC");
  TH2F *hBpiSigmaVspTPC=new TH2F(hisname.Data(),hisname.Data(),100,0.,5.0,100,-10.0,10.0);
  fOutputMC->Add(hBpiSigmaVspTPC);

  //Jaime Lc specific
  hisname.Form("hLcRealTot");
  TH1F *hLambdaRealTot=new TH1F(hisname.Data(),hisname.Data(),100,0.,15.0);
  fOutputMC->Add(hLambdaRealTot);
  hisname.Form("hbLcRealTot");
  TH1F *hbLambdaRealTot=new TH1F(hisname.Data(),hisname.Data(),100,0.,15.0);
  fOutputMC->Add(hbLambdaRealTot);
  hisname.Form("hLcIDTot");
  TH1F *hLambdaIDTot=new TH1F(hisname.Data(),hisname.Data(),100,0.,15.0);
  fOutputMC->Add(hLambdaIDTot);
  hisname.Form("hbLcIDTot");
  TH1F *hbLambdaIDTot=new TH1F(hisname.Data(),hisname.Data(),100,0.,15.0);
  fOutputMC->Add(hbLambdaIDTot);
  hisname.Form("hLcIDGood");
  TH1F *hLambdaIDGood=new TH1F(hisname.Data(),hisname.Data(),100,0.,15.0);
  fOutputMC->Add(hLambdaIDGood);
  hisname.Form("hbLcIDGood");
  TH1F *hbLambdaIDGood=new TH1F(hisname.Data(),hisname.Data(),100,0.,15.0);
  fOutputMC->Add(hbLambdaIDGood);
  hisname.Form("hLcnoID");
  TH1F *hLambdanoID=new TH1F(hisname.Data(),hisname.Data(),100,0.,15.0);
  fOutputMC->Add(hLambdanoID);
  hisname.Form("hbLcnoID");
  TH1F *hbLambdanoID=new TH1F(hisname.Data(),hisname.Data(),100,0.,15.0);
  fOutputMC->Add(hbLambdanoID);
  
  // other generic 
  hisname.Form("hLcpt");
  TH1F *hLcPt=new TH1F(hisname.Data(),hisname.Data(),100,0.,5.0);
  fOutputMC->Add(hLcPt);
  hisname.Form("hbLcpt");
  TH1F *hBLcPt=new TH1F(hisname.Data(),hisname.Data(),100,0.,5.0);
  fOutputMC->Add(hBLcPt);
  hisname.Form("hDist12toPrim");
  TH1F *hDist12Prim=new TH1F(hisname.Data(),hisname.Data(),100,0.,1.0);
  fOutputMC->Add(hDist12Prim);
  hisname.Form("hbDist12toPrim");
  TH1F *hBDist12Prim=new TH1F(hisname.Data(),hisname.Data(),100,0.,1.0);
  fOutputMC->Add(hBDist12Prim);

  hisname.Form("hSigmaVert");
  TH1F *hSigmaVert=new TH1F(hisname.Data(),hisname.Data(),60,0.,0.06);
  fOutputMC->Add(hSigmaVert);
  hisname.Form("hbSigmaVert");
  TH1F *hBSigmaVert=new TH1F(hisname.Data(),hisname.Data(),60,0.,0.06);
  fOutputMC->Add(hBSigmaVert);

  hisname.Form("hDCAs");
  TH1F *hdcas=new TH1F(hisname.Data(),hisname.Data(),200,0.,0.1);
  fOutputMC->Add(hdcas);
  hisname.Form("hbDCAs");
  TH1F *hBdcas=new TH1F(hisname.Data(),hisname.Data(),200,0.,0.1);
  fOutputMC->Add(hBdcas);

  hisname.Form("hCosPointingAngle");
  TH1F *hCosPointingAngle=new TH1F(hisname.Data(),hisname.Data(),40,0.,1.);
  fOutputMC->Add(hCosPointingAngle);
  hisname.Form("hbCosPointingAngle");
  TH1F *hBCosPointingAngle=new TH1F(hisname.Data(),hisname.Data(),40,0.,1.);
  fOutputMC->Add(hBCosPointingAngle);

  hisname.Form("hDecayLength");
  TH1F *hDecayLength=new TH1F(hisname.Data(),hisname.Data(),100,0.,0.1);
  fOutputMC->Add(hDecayLength);
  hisname.Form("hbDecayLength");
  TH1F *hBDecayLength=new TH1F(hisname.Data(),hisname.Data(),100,0.,0.1);
  fOutputMC->Add(hBDecayLength);

  hisname.Form("hSum2");
  TH1F *hSum2=new TH1F(hisname.Data(),hisname.Data(),100,0.,0.1);
  fOutputMC->Add(hSum2);
  hisname.Form("hbSum2");
  TH1F *hBSum2=new TH1F(hisname.Data(),hisname.Data(),100,0.,0.1);
  fOutputMC->Add(hBSum2);

  hisname.Form("hptmax");
  TH1F *hPtMax=new TH1F(hisname.Data(),hisname.Data(),100,0.,5.);
  fOutputMC->Add(hPtMax);
  hisname.Form("hbptmax");
  TH1F *hBPtMax=new TH1F(hisname.Data(),hisname.Data(),100,0.,5.);
  fOutputMC->Add(hBPtMax);

//2D Var x Lc Pt
  hisname.Form("hpptProngLcPt");
  fhProtonPtProngLcPt=new TH2F(hisname.Data(),hisname.Data(),75,0.,15.,100,0.,5.0);
  fOutputMC->Add(fhProtonPtProngLcPt);
  hisname.Form("hbpptProngLcPt");
  fhBProtonPtProngLcPt=new TH2F(hisname.Data(),hisname.Data(),75,0.,15.,100,0.,5.0);
  fOutputMC->Add(fhBProtonPtProngLcPt);
  
  hisname.Form("hpd0ProngLcPt");
  fhProtond0ProngLcPt=new TH2F(hisname.Data(),hisname.Data(),75,0.,15.,100,-0.1,0.1);
  fOutputMC->Add(fhProtond0ProngLcPt);
  hisname.Form("hbpd0ProngLcPt");
  fhBProtond0ProngLcPt=new TH2F(hisname.Data(),hisname.Data(),75,0.,15.,100,-0.1,0.1);
  fOutputMC->Add(fhBProtond0ProngLcPt);
 
  hisname.Form("hKptProngLcPt");
  fhKaonPtProngLcPt=new TH2F(hisname.Data(),hisname.Data(),75,0.,15.,100,0.,5.0);
  fOutputMC->Add(fhKaonPtProngLcPt);
  hisname.Form("hbKptProngLcPt");
  fhBKaonPtProngLcPt=new TH2F(hisname.Data(),hisname.Data(),75,0.,15.,100,0.,5.0);
  fOutputMC->Add(fhBKaonPtProngLcPt);
   
  hisname.Form("hKd0ProngLcPt");
  fhKaond0ProngLcPt=new TH2F(hisname.Data(),hisname.Data(),75,0.,15.,100,-0.1,0.1);
  fOutputMC->Add(fhKaond0ProngLcPt);
  hisname.Form("hbKd0ProngLcPt");
  fhBKaond0ProngLcPt=new TH2F(hisname.Data(),hisname.Data(),75,0.,15.,100,-0.1,0.1);
  fOutputMC->Add(fhBKaond0ProngLcPt);

  hisname.Form("hpiptProngLcPt");
  fhPionPtProngLcPt=new TH2F(hisname.Data(),hisname.Data(),75,0.,15.,100,0.,5.0);
  fOutputMC->Add(fhPionPtProngLcPt);
  hisname.Form("hbpiptProngLcPt");
  fhBPionPtProngLcPt=new TH2F(hisname.Data(),hisname.Data(),75,0.,15.,100,0.,5.0);
  fOutputMC->Add(fhBPionPtProngLcPt);

  hisname.Form("hpid0ProngLcPt");
  fhPiond0ProngLcPt=new TH2F(hisname.Data(),hisname.Data(),75,0.,15.,100,-0.1,-0.1);
  fOutputMC->Add(fhPiond0ProngLcPt);
  hisname.Form("hbpid0ProngLcPt");
  fhBPiond0ProngLcPt=new TH2F(hisname.Data(),hisname.Data(),75,0.,15.,100,-0.1,0.1);
  fOutputMC->Add(fhBPiond0ProngLcPt);
 
  hisname.Form("hDist12toPrimLcPt");
  fhDist12PrimLcPt = new TH2F(hisname.Data(),hisname.Data(),75,0.,15.,100,0.,1.0);
  fOutputMC->Add(fhDist12PrimLcPt);
  hisname.Form("hbDist12toPrimLcPt");
  fhBDist12PrimLcPt = new TH2F(hisname.Data(),hisname.Data(),75,0.,15.,100,0.,1.0);
  fOutputMC->Add(fhBDist12PrimLcPt);

  hisname.Form("hSigmaVertLcPt");
  fhSigmaVertLcPt=new TH2F(hisname.Data(),hisname.Data(),75,0.,15.,60,0.,0.06);
  fOutputMC->Add(fhSigmaVertLcPt);
  hisname.Form("hbSigmaVertLcPt");
  fhBSigmaVertLcPt=new TH2F(hisname.Data(),hisname.Data(),75,0.,15.,60,0.,0.06);
  fOutputMC->Add(fhBSigmaVertLcPt);

  hisname.Form("hDCAsLcPt");
  fhdcasLcPt=new TH2F(hisname.Data(),hisname.Data(),75,0.,15.,200,0.,0.1);
  fOutputMC->Add(fhdcasLcPt);
  hisname.Form("hbDCAsLcPt");
  fhBdcasLcPt=new TH2F(hisname.Data(),hisname.Data(),75,0.,15.,200,0.,0.1);
  fOutputMC->Add(fhBdcasLcPt);
  
  hisname.Form("hCosPointingAngleLcPt");
  fhCosPointingAngleLcPt = new TH2F(hisname.Data(),hisname.Data(),75,0.,15.,40,0.,1.);
  fOutputMC->Add(fhCosPointingAngleLcPt);
  hisname.Form("hbCosPointingAngleLcPt");
  fhBCosPointingAngleLcPt = new TH2F(hisname.Data(),hisname.Data(),75,0.,15.,40,0.,1.);
  fOutputMC->Add(fhBCosPointingAngleLcPt);

  hisname.Form("hDecayLengthLcPt");
  fhDecayLengthLcPt	= new TH2F(hisname.Data(),hisname.Data(),75,0.,15.,100,0.,0.1);
  fOutputMC->Add(fhDecayLengthLcPt);
  hisname.Form("hbDecayLengthLcPt");
  fhBDecayLengthLcPt=new TH2F(hisname.Data(),hisname.Data(),75,0.,15.,100,0.,0.1);
  fOutputMC->Add(fhBDecayLengthLcPt);

  hisname.Form("hSum2LcPt");
  fhSum2LcPt=new TH2F(hisname.Data(),hisname.Data(),75,0.,15.,100,0.,0.1);
  fOutputMC->Add(fhSum2LcPt);
  hisname.Form("hbSum2LcPt");
  fhBSum2LcPt=new TH2F(hisname.Data(),hisname.Data(),75,0.,15.,100,0.,0.1);
  fOutputMC->Add(fhBSum2LcPt);

  hisname.Form("hptmaxLcPt");
  fhPtMaxLcPt=new TH2F(hisname.Data(),hisname.Data(),75,0.,15.,100,0.,5.);
  fOutputMC->Add(fhPtMaxLcPt);
  hisname.Form("hbptmaxLcPt");
  fhBPtMaxLcPt=new TH2F(hisname.Data(),hisname.Data(),75,0.,15.,100,0.,5.);
  fOutputMC->Add(fhBPtMaxLcPt); 

  fAPriori = new TList(); // AdC
  fAPriori->SetOwner(); // AdC
  fAPriori->SetName("APrioriMCHistos"); // AdC

  hisname.Form("hElIn3Prong");
  TH1F *hElIn3Prong=new TH1F(hisname.Data(),hisname.Data(),100,0.,5.);
  fAPriori->Add(hElIn3Prong);
  hisname.Form("hMuIn3Prong");
  TH1F *hMuIn3Prong=new TH1F(hisname.Data(),hisname.Data(),100,0.,5.);
  fAPriori->Add(hMuIn3Prong);
  hisname.Form("hPiIn3Prong");
  TH1F *hPiIn3Prong=new TH1F(hisname.Data(),hisname.Data(),100,0.,5.);
  fAPriori->Add(hPiIn3Prong);
  hisname.Form("hKaIn3Prong");
  TH1F *hKaIn3Prong=new TH1F(hisname.Data(),hisname.Data(),100,0.,5.);
  fAPriori->Add(hKaIn3Prong);
  hisname.Form("hPrIn3Prong");
  TH1F *hPrIn3Prong=new TH1F(hisname.Data(),hisname.Data(),100,0.,5.);
  fAPriori->Add(hPrIn3Prong);

  hisname.Form("hElIn3Prong1");
  TH1F *hElIn3Prong1=new TH1F(hisname.Data(),hisname.Data(),100,0.,5.);
  fAPriori->Add(hElIn3Prong1);
  hisname.Form("hMuIn3Prong1");
  TH1F *hMuIn3Prong1=new TH1F(hisname.Data(),hisname.Data(),100,0.,5.);
  fAPriori->Add(hMuIn3Prong1);
  hisname.Form("hPiIn3Prong1");
  TH1F *hPiIn3Prong1=new TH1F(hisname.Data(),hisname.Data(),100,0.,5.);
  fAPriori->Add(hPiIn3Prong1);
  hisname.Form("hKaIn3Prong1");
  TH1F *hKaIn3Prong1=new TH1F(hisname.Data(),hisname.Data(),100,0.,5.);
  fAPriori->Add(hKaIn3Prong1);
  hisname.Form("hPrIn3Prong1");
  TH1F *hPrIn3Prong1=new TH1F(hisname.Data(),hisname.Data(),100,0.,5.);
  fAPriori->Add(hPrIn3Prong1);

  hisname.Form("hElIn3Prong2");
  TH1F *hElIn3Prong2=new TH1F(hisname.Data(),hisname.Data(),100,0.,5.);
  fAPriori->Add(hElIn3Prong2);
  hisname.Form("hMuIn3Prong2");
  TH1F *hMuIn3Prong2=new TH1F(hisname.Data(),hisname.Data(),100,0.,5.);
  fAPriori->Add(hMuIn3Prong2);
  hisname.Form("hPiIn3Prong2");
  TH1F *hPiIn3Prong2=new TH1F(hisname.Data(),hisname.Data(),100,0.,5.);
  fAPriori->Add(hPiIn3Prong2);
  hisname.Form("hKaIn3Prong2");
  TH1F *hKaIn3Prong2=new TH1F(hisname.Data(),hisname.Data(),100,0.,5.);
  fAPriori->Add(hKaIn3Prong2);
  hisname.Form("hPrIn3Prong2");
  TH1F *hPrIn3Prong2=new TH1F(hisname.Data(),hisname.Data(),100,0.,5.);
  fAPriori->Add(hPrIn3Prong2);

  hisname.Form("hElIn3Prong3");
  TH1F *hElIn3Prong3=new TH1F(hisname.Data(),hisname.Data(),100,0.,5.);
  fAPriori->Add(hElIn3Prong3);
  hisname.Form("hMuIn3Prong3");
  TH1F *hMuIn3Prong3=new TH1F(hisname.Data(),hisname.Data(),100,0.,5.);
  fAPriori->Add(hMuIn3Prong3);
  hisname.Form("hPiIn3Prong3");
  TH1F *hPiIn3Prong3=new TH1F(hisname.Data(),hisname.Data(),100,0.,5.);
  fAPriori->Add(hPiIn3Prong3);
  hisname.Form("hKaIn3Prong3");
  TH1F *hKaIn3Prong3=new TH1F(hisname.Data(),hisname.Data(),100,0.,5.);
  fAPriori->Add(hKaIn3Prong3);
  hisname.Form("hPrIn3Prong3");
  TH1F *hPrIn3Prong3=new TH1F(hisname.Data(),hisname.Data(),100,0.,5.);
  fAPriori->Add(hPrIn3Prong3);

  hisname.Form("hElIn3Prong4");
  TH1F *hElIn3Prong4=new TH1F(hisname.Data(),hisname.Data(),100,0.,5.);
  fAPriori->Add(hElIn3Prong4);
  hisname.Form("hMuIn3Prong4");
  TH1F *hMuIn3Prong4=new TH1F(hisname.Data(),hisname.Data(),100,0.,5.);
  fAPriori->Add(hMuIn3Prong4);
  hisname.Form("hPiIn3Prong4");
  TH1F *hPiIn3Prong4=new TH1F(hisname.Data(),hisname.Data(),100,0.,5.);
  fAPriori->Add(hPiIn3Prong4);
  hisname.Form("hKaIn3Prong4");
  TH1F *hKaIn3Prong4=new TH1F(hisname.Data(),hisname.Data(),100,0.,5.);
  fAPriori->Add(hKaIn3Prong4);
  hisname.Form("hPrIn3Prong4");
  TH1F *hPrIn3Prong4=new TH1F(hisname.Data(),hisname.Data(),100,0.,5.);
  fAPriori->Add(hPrIn3Prong4);

  hisname.Form("hElIn3Prong5");
  TH1F *hElIn3Prong5=new TH1F(hisname.Data(),hisname.Data(),100,0.,5.);
  fAPriori->Add(hElIn3Prong5);
  hisname.Form("hMuIn3Prong5");
  TH1F *hMuIn3Prong5=new TH1F(hisname.Data(),hisname.Data(),100,0.,5.);
  fAPriori->Add(hMuIn3Prong5);
  hisname.Form("hPiIn3Prong5");
  TH1F *hPiIn3Prong5=new TH1F(hisname.Data(),hisname.Data(),100,0.,5.);
  fAPriori->Add(hPiIn3Prong5);
  hisname.Form("hKaIn3Prong5");
  TH1F *hKaIn3Prong5=new TH1F(hisname.Data(),hisname.Data(),100,0.,5.);
  fAPriori->Add(hKaIn3Prong5);
  hisname.Form("hPrIn3Prong5");
  TH1F *hPrIn3Prong5=new TH1F(hisname.Data(),hisname.Data(),100,0.,5.);
  fAPriori->Add(hPrIn3Prong5);

  hisname.Form("hElIn3Prong6");
  TH1F *hElIn3Prong6=new TH1F(hisname.Data(),hisname.Data(),100,0.,5.);
  fAPriori->Add(hElIn3Prong6);
  hisname.Form("hMuIn3Prong6");
  TH1F *hMuIn3Prong6=new TH1F(hisname.Data(),hisname.Data(),100,0.,5.);
  fAPriori->Add(hMuIn3Prong6);
  hisname.Form("hPiIn3Prong6");
  TH1F *hPiIn3Prong6=new TH1F(hisname.Data(),hisname.Data(),100,0.,5.);
  fAPriori->Add(hPiIn3Prong6);
  hisname.Form("hKaIn3Prong6");
  TH1F *hKaIn3Prong6=new TH1F(hisname.Data(),hisname.Data(),100,0.,5.);
  fAPriori->Add(hKaIn3Prong6);
  hisname.Form("hPrIn3Prong6");
  TH1F *hPrIn3Prong6=new TH1F(hisname.Data(),hisname.Data(),100,0.,5.);
  fAPriori->Add(hPrIn3Prong6);

  
  fMultiplicity = new TList(); // AdC
  fMultiplicity->SetOwner(); // AdC
  fMultiplicity->SetName("MultiplicityMCHistos"); // AdC

  hisname.Form("hLcinEvent");
  TH1I*hLcinEvent=new TH1I(hisname.Data(),hisname.Data(),3,-1,2);
  fMultiplicity->Add(hLcinEvent);

  hisname.Form("hPrimariesvsAOD");
  TH2I *hPrimariesvsAOD=new TH2I(hisname.Data(),hisname.Data(),300,0,300,300,0,300);
  fMultiplicity->Add(hPrimariesvsAOD);

  hisname.Form("hMultiplicityInLcEvent");
  TH1I *hMultiplicityInLcEvent=new TH1I(hisname.Data(),hisname.Data(),300,0,300);
  fMultiplicity->Add(hMultiplicityInLcEvent);
  hisname.Form("h2MultiplicityInLcEvent");
  TH1I *h2MultiplicityInLcEvent=new TH1I(hisname.Data(),hisname.Data(),300,0,300);
  fMultiplicity->Add(h2MultiplicityInLcEvent);
  hisname.Form("hAll2MultiplicityInEvent");
  TH1I *hAll2MultiplicityInEvent=new TH1I(hisname.Data(),hisname.Data(),300,0,300);
  fMultiplicity->Add(hAll2MultiplicityInEvent);
  hisname.Form("hAllMultiplicityInEvent");
  TH1I *hAllMultiplicityInEvent=new TH1I(hisname.Data(),hisname.Data(),300,0,300);
  fMultiplicity->Add(hAllMultiplicityInEvent);
  hisname.Form("hAllMultiplicityPrimaryInEvent");
  TH1I *hAllMultiplicityPrimaryInEvent=new TH1I(hisname.Data(),hisname.Data(),300,0,300);
  fMultiplicity->Add(hAllMultiplicityPrimaryInEvent);
  hisname.Form("hAll2MultiplicityPrimaryInEvent");
  TH1I *hAll2MultiplicityPrimaryInEvent=new TH1I(hisname.Data(),hisname.Data(),300,0,300);
  fMultiplicity->Add(hAll2MultiplicityPrimaryInEvent);
  hisname.Form("hMultiplicityInEvent");
  TH1I *hMultiplicityInEvent=new TH1I(hisname.Data(),hisname.Data(),300,0.,300);
  fMultiplicity->Add(hMultiplicityInEvent);
  hisname.Form("hMultiplicityIn3ProngLC");
  TH1I *hMultiplicityIn3ProngLC=new TH1I(hisname.Data(),hisname.Data(),300,0,300);
  fMultiplicity->Add(hMultiplicityIn3ProngLC);
  hisname.Form("hMultiplicityInLCpid");
  TH1I *hMultiplicityInLCpid=new TH1I(hisname.Data(),hisname.Data(),300,0,300);
  fMultiplicity->Add(hMultiplicityInLCpid);
  hisname.Form("hMultiplicityInLCmc");
  TH1I *hMultiplicityInLCmc=new TH1I(hisname.Data(),hisname.Data(),300,0,300);
  fMultiplicity->Add(hMultiplicityInLCmc);
  hisname.Form("hMultiplicityInLCNomc");
  TH1I *hMultiplicityInLCNomc=new TH1I(hisname.Data(),hisname.Data(),300,0,300);
  fMultiplicity->Add(hMultiplicityInLCNomc);
  hisname.Form("hMultiplicityYesC");
  TH1I *hMultiplicityYesC=new TH1I(hisname.Data(),hisname.Data(),300,0,300);
  fMultiplicity->Add(hMultiplicityYesC);
  hisname.Form("hMultiplicityYesB");
  TH1I *hMultiplicityYesB=new TH1I(hisname.Data(),hisname.Data(),300,0,300);
  fMultiplicity->Add(hMultiplicityYesB);
  hisname.Form("hMultiplicityJPsi");
  TH1I *hMultiplicityJPsi=new TH1I(hisname.Data(),hisname.Data(),300,0,300);
  fMultiplicity->Add(hMultiplicityJPsi);

  
  AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
  fPIDResponse = inputHandler->GetPIDResponse();

  if(fRDCutsProduction->GetIsUsePID()){
    fRDCutsProduction->GetPidHF()->SetPidResponse(fPIDResponse);
    fRDCutsProduction->GetPidpion()->SetPidResponse(fPIDResponse);
    fRDCutsProduction->GetPidprot()->SetPidResponse(fPIDResponse);
    //fUtilPid=new AliAODpidUtil(fPIDResponse);
    fRDCutsProduction->GetPidHF()->SetOldPid(kFALSE);
    fRDCutsProduction->GetPidpion()->SetOldPid(kFALSE);
    fRDCutsProduction->GetPidprot()->SetOldPid(kFALSE);
  }
  if(fRDCutsAnalysis->GetIsUsePID()){
    fRDCutsAnalysis->GetPidHF()->SetPidResponse(fPIDResponse);
    fRDCutsAnalysis->GetPidpion()->SetPidResponse(fPIDResponse);
    fRDCutsAnalysis->GetPidprot()->SetPidResponse(fPIDResponse);
    fRDCutsAnalysis->GetPidHF()->SetOldPid(kFALSE);
    fRDCutsAnalysis->GetPidpion()->SetOldPid(kFALSE);
    fRDCutsAnalysis->GetPidprot()->SetOldPid(kFALSE);
  }

  PostData(1,fOutput);
  if (fFillVarHists) PostData(3,fOutputMC);
  PostData(4,fNentries);
  if (fPriorsHists) PostData(5,fAPriori);
  if (fMultiplicityHists) PostData(6,fMultiplicity);
   TString normName="NormalizationCounter";
   AliAnalysisDataContainer *cont = GetOutputSlot(7)->GetContainer();
  if(cont)normName=(TString)cont->GetName(); 
  fCounter = new AliNormalizationCounter(normName.Data());
   fCounter->Init();
    PostData(7,fCounter);
  if (fFillNtuple) {
    //OpenFile(3); // 2 is the slot number of the ntuple
    fNtupleLambdac = new TNtuple("fNtupleLambdac","D +",
				 "pdg:Px:Py:Pz:PtTrue:VxTrue:VyTrue:VzTrue:Ptpi:PtK:Ptpi2:PtRec:PointingAngle:DecLeng:VxRec:VyRec:VzRec:InvMass:sigvert:d0Pi:d0K:d0Pi2:dca:d0square");  
    PostData(8,fNtupleLambdac);
  }

  return;
}

//________________________________________________________________________
void AliAnalysisTaskSELambdac::UserExec(Option_t */*option*/)
{
  // Execute analysis for current event:
  // heavy flavor candidates association to MC truth

  AliAODEvent *aod = dynamic_cast<AliAODEvent*> (InputEvent());
  //tmp
  fHistNEvents->Fill(0); // count event
  // Post the data already here

  TClonesArray *array3Prong = 0;
  TClonesArray *arrayLikeSign =0;
  if(!aod && AODEvent() && IsStandardAOD()) {
    // In case there is an AOD handler writing a standard AOD, use the AOD 
    // event in memory rather than the input (ESD) event.    
    aod = dynamic_cast<AliAODEvent*> (AODEvent());
    // in this case the braches in the deltaAOD (AliAOD.VertexingHF.root)
    // have to taken from the AOD event hold by the AliAODExtension
    AliAODHandler* aodHandler = (AliAODHandler*) 
      ((AliAnalysisManager::GetAnalysisManager())->GetOutputEventHandler());
    if(aodHandler->GetExtensions()) {
      AliAODExtension *ext = (AliAODExtension*)aodHandler->GetExtensions()->FindObject("AliAOD.VertexingHF.root");
      AliAODEvent *aodFromExt = ext->GetAOD();
      array3Prong=(TClonesArray*)aodFromExt->GetList()->FindObject("Charm3Prong");
      arrayLikeSign=(TClonesArray*)aodFromExt->GetList()->FindObject("LikeSign3Prong");
    }
  } else if(aod) {
    array3Prong=(TClonesArray*)aod->GetList()->FindObject("Charm3Prong");
    arrayLikeSign=(TClonesArray*)aod->GetList()->FindObject("LikeSign3Prong");
  }

  if(!aod) return;

  TClonesArray *arrayMC=0;
  AliAODMCHeader *mcHeader=0;

  // load MC particles
  if(fReadMC){
    
    arrayMC =  (TClonesArray*)aod->GetList()->FindObject(AliAODMCParticle::StdBranchName());
    if(!arrayMC) {
      AliError("AliAnalysisTaskSELambdac::UserExec: MC particles branch not found!\n");
      return;
    }
 

    // load MC header
    mcHeader =  (AliAODMCHeader*)aod->GetList()->FindObject(AliAODMCHeader::StdBranchName());
    if(!mcHeader) {
      AliError("AliAnalysisTaskSELambdac::UserExec: MC header branch not found!\n");
      return;
    }
  }

  TString fillthis="";
  Int_t numberOfPrimaries= NumberPrimaries(aod);

  if (fMultiplicityHists && fReadMC) {
    fillthis="hPrimariesvsAOD";
    ((TH1I*)fMultiplicity->FindObject(fillthis))->Fill(aod->GetNumberOfTracks(),numberOfPrimaries);

    fillthis="hAll2MultiplicityInEvent";
    ((TH1I*)fMultiplicity->FindObject(fillthis))->Fill(aod->GetNumberOfTracks());

    fillthis="hAll2MultiplicityPrimaryInEvent";
    ((TH1I*)fMultiplicity->FindObject(fillthis))->Fill(numberOfPrimaries);
   
    if (IsThereAGeneratedLc(arrayMC)) {
      fillthis="h2MultiplicityInLcEvent";
      ((TH1I*)fMultiplicity->FindObject(fillthis))->Fill(numberOfPrimaries);
    }

  }

  if(!array3Prong || !aod) {
    AliError("AliAnalysisTaskSELambdac::UserExec: Charm3Prong branch not found!\n");
    return;
  }
  if(!arrayLikeSign) {
    AliDebug(2,"AliAnalysisTaskSELambdac::UserExec: LikeSign3Prong branch not found!\n");
    //  return;
  }

  // fix for temporary bug in ESDfilter
  // the AODs with null vertex pointer didn't pass the PhysSel
  AliAODVertex *vtx1 = (AliAODVertex*)aod->GetPrimaryVertex();
  if(!vtx1 || TMath::Abs(aod->GetMagneticField())<0.001) return;

  fNentries->Fill(0);
  fCounter->StoreEvent(aod,fRDCutsProduction,fReadMC);
  TString trigclass=aod->GetFiredTriggerClasses();
  if(trigclass.Contains("C0SMH-B-NOPF-ALLNOTRD") || trigclass.Contains("C0SMH-B-NOPF-ALL")) fNentries->Fill(14);
  Bool_t isEvSelProdCuts=fRDCutsProduction->IsEventSelected(aod);
  Bool_t isEvSelAnCuts=fRDCutsAnalysis->IsEventSelected(aod);
  if(!isEvSelProdCuts || !isEvSelAnCuts){
    if(fRDCutsProduction->GetWhyRejection()==1) // rejected for pileup
      fNentries->Fill(13);
    return;
  }

  Bool_t isThereA3prongWithGoodTracks = kFALSE;
  Bool_t isThereA3ProngLcKine = kFALSE;
  Bool_t isThereA3ProngLcKineANDpid = kFALSE;
  Bool_t isThereA3ProngLcMC = kFALSE;
  Bool_t isThereA3ProngCyes = kFALSE;
  Bool_t isThereA3ProngByes = kFALSE;
  Bool_t isThereA3ProngJPsi = kFALSE;

  Int_t n3Prong = array3Prong->GetEntriesFast();
  Int_t nSelectedloose[1]={0};
  Int_t nSelectedtight[1]={0};

  // vHF object is needed to call the method that refills the missing info of the candidates
  // if they have been deleted in dAOD reconstruction phase
  // in order to reduce the size of the file
  AliAnalysisVertexingHF *vHF=new AliAnalysisVertexingHF();
  for (Int_t i3Prong = 0; i3Prong < n3Prong; i3Prong++) {
    AliAODRecoDecayHF3Prong *d = (AliAODRecoDecayHF3Prong*)array3Prong->UncheckedAt(i3Prong);

    //Filter bit selection and QA:
    fhSelectBit->Fill(1);
    if(d->GetSelectionMap()){
      fhSelectBit->Fill(2);
      if(!d->HasSelectionBit(AliRDHFCuts::kLcCuts))  fhSelectBit->Fill(3);
      if(!d->HasSelectionBit(AliRDHFCuts::kLcPID))   fhSelectBit->Fill(4);
      if(fLcCut&&!d->HasSelectionBit(AliRDHFCuts::kLcCuts))		continue;
      if(fLcPIDCut&&!d->HasSelectionBit(AliRDHFCuts::kLcPID))		continue;
    }

    if(!(vHF->FillRecoCand(aod,d))) {////Fill the data members of the candidate only if they are empty.
      fNentries->Fill(15); //monitor how often this fails 
      continue;
    }
     Bool_t unsetvtx=kFALSE;
    if(!d->GetOwnPrimaryVtx()){
      d->SetOwnPrimaryVtx(vtx1);
      unsetvtx=kTRUE;
    }


    Int_t isSelectedTracks = fRDCutsProduction->IsSelected(d,AliRDHFCuts::kTracks,aod);
    if(!isSelectedTracks) continue;
    fhEta3Prong->Fill(d->Pt(),d->Eta());
    fhRap3Prong->Fill(d->Pt(),d->Y(4122));
    
    isThereA3prongWithGoodTracks=kTRUE;

    if (fRDCutsProduction->IsInFiducialAcceptance(d->Pt(),d->Y(4122))) {fNentries->Fill(6);}else{continue;}
    fhEta3ProngAcc->Fill(d->Pt(),d->Eta());
    fhRap3ProngAcc->Fill(d->Pt(),d->Y(4122));
    
    Int_t ptbin=fRDCutsProduction->PtBin(d->Pt());
    if(ptbin==-1) {fNentries->Fill(4); continue;} //out of bounds

    FillMassHists(aod,d,arrayMC,fRDCutsProduction,nSelectedloose,nSelectedtight);
    if(fFillVarHists) FillVarHists(d,arrayMC,fRDCutsProduction,/*fOutputMC,*/aod);

    if (fPriorsHists && fReadMC)
      FillAPrioriConcentrations(d, fRDCutsProduction, aod, arrayMC);
    if (fMultiplicityHists && fReadMC)
      MultiplicityStudies(d, fRDCutsProduction, aod, arrayMC,
			  isThereA3ProngLcKine,isThereA3ProngLcKineANDpid,isThereA3ProngLcMC,
			  isThereA3ProngCyes,isThereA3ProngByes,isThereA3ProngJPsi);

    /*
    //start OS analysis
    if(labDp<0)fHistOSbkg->Fill(d->InvMassDplus());
    fHistOS->Fill(d->InvMassDplus());
    */

    if(unsetvtx) d->UnsetOwnPrimaryVtx();
  }
  fCounter->StoreCandidates(aod,nSelectedloose[0],kTRUE);
  fCounter->StoreCandidates(aod,nSelectedtight[0],kFALSE);

  if (fMultiplicityHists && fReadMC) {
    fillthis="hAllMultiplicityInEvent";
    ((TH1I*)fMultiplicity->FindObject(fillthis))->Fill(aod->GetNumberOfTracks());

    fillthis="hAllMultiplicityPrimaryInEvent";
    ((TH1I*)fMultiplicity->FindObject(fillthis))->Fill(numberOfPrimaries);
   
    if (IsThereAGeneratedLc(arrayMC)) {
      fillthis="hMultiplicityInLcEvent";
      ((TH1I*)fMultiplicity->FindObject(fillthis))->Fill(numberOfPrimaries);
    }

    if (isThereA3prongWithGoodTracks) {
      fillthis="hMultiplicityInEvent";
      ((TH1I*)fMultiplicity->FindObject(fillthis))->Fill(numberOfPrimaries);
    }
    if (isThereA3ProngLcKine) {
      fillthis="hMultiplicityIn3ProngLC";
      ((TH1I*)fMultiplicity->FindObject(fillthis))->Fill(numberOfPrimaries);
    }
    if (isThereA3ProngLcKineANDpid) {
      fillthis="hMultiplicityInLCpid";
      ((TH1I*)fMultiplicity->FindObject(fillthis))->Fill(numberOfPrimaries);
    }
    if (isThereA3ProngLcMC) {
      fillthis="hMultiplicityInLCmc";
      ((TH1I*)fMultiplicity->FindObject(fillthis))->Fill(numberOfPrimaries);
    }
    if (isThereA3ProngLcKine && !isThereA3ProngLcMC) {
      fillthis="hMultiplicityInLCNomc";
      ((TH1I*)fMultiplicity->FindObject(fillthis))->Fill(numberOfPrimaries);
    }

    if (isThereA3ProngCyes) {
      fillthis="hMultiplicityYesC";
      ((TH1I*)fMultiplicity->FindObject(fillthis))->Fill(numberOfPrimaries);
    }
    if (isThereA3ProngByes) {
      fillthis="hMultiplicityYesB";
      ((TH1I*)fMultiplicity->FindObject(fillthis))->Fill(numberOfPrimaries);
    }
    if (isThereA3ProngJPsi) {
      fillthis="hMultiplicityJPsi";
      ((TH1I*)fMultiplicity->FindObject(fillthis))->Fill(numberOfPrimaries);
    }

  }

  delete vHF;
   
  PostData(1,fOutput); 
  if (fFillVarHists) PostData(3,fOutputMC);
  if (fPriorsHists) PostData(5,fAPriori);
  if (fMultiplicityHists) PostData(6,fMultiplicity);

  PostData(4,fNentries);
  PostData(7,fCounter);
      
  return;
}



//________________________________________________________________________
void AliAnalysisTaskSELambdac::Terminate(Option_t */*option*/)
{
  // Terminate analysis
  //
  if (fDebug > 1) printf("AnalysisTaskSELambdac: Terminate() \n");

  fOutput = dynamic_cast<TList*> (GetOutputData(1));
  if (!fOutput) {     
    AliError("ERROR: fOutput not available\n");
    return;
  }
  //fHistNEvents = dynamic_cast<TH1F*>(fOutput->FindObject("fHistNEvents"));

  if(fFillNtuple){
    fNtupleLambdac = dynamic_cast<TNtuple*>(GetOutputData(3));
  }

 
  return;
}

//________________________________________________________________________
Int_t AliAnalysisTaskSELambdac::MatchToMCLambdac(AliAODRecoDecayHF3Prong *d,TClonesArray *arrayMC) const{
  // check if the candidate is a Lambdac decaying in pKpi or in the resonant channels
  Int_t lambdacLab[3]={0,0,0};
  //  Int_t pdgs[3]={0,0,0};
  for(Int_t i=0;i<3;i++){
    AliAODTrack *daugh=(AliAODTrack*)d->GetDaughter(i);
    Int_t lab=daugh->GetLabel();
    if(lab<0) return 0;
    AliAODMCParticle *part= (AliAODMCParticle*)arrayMC->At(lab);
    if(!part) continue;
    //    pdgs[i]=part->GetPdgCode();
    Int_t partPdgcode = TMath::Abs(part->GetPdgCode());
    if(partPdgcode==211 || partPdgcode==321 || partPdgcode==2212){
      Int_t motherLabel=part->GetMother();
      if(motherLabel<0) return 0;
      AliAODMCParticle *motherPart = (AliAODMCParticle*)arrayMC->At(motherLabel);
      if(!motherPart) continue;
      Int_t motherPdg = TMath::Abs(motherPart->GetPdgCode());
      if(motherPdg==4122) {
        if(GetLambdacDaugh(motherPart,arrayMC)){lambdacLab[i]=motherLabel;continue;}
      }
      if(motherPdg==313 || motherPdg==2224 || motherPdg==3124){
	Int_t granMotherLabel=motherPart->GetMother();
	if(granMotherLabel<0) return 0;
	AliAODMCParticle *granMotherPart = (AliAODMCParticle*)arrayMC->At(granMotherLabel);
	if(!granMotherPart) continue;
	Int_t granMotherPdg  = TMath::Abs(granMotherPart->GetPdgCode());
	if(granMotherPdg ==4122) {
	  if(GetLambdacDaugh(granMotherPart,arrayMC)) {lambdacLab[i]=granMotherLabel;continue;}
	}
      }
    }
  }

  if(lambdacLab[0]==lambdacLab[1] && lambdacLab[1]==lambdacLab[2]) {return lambdacLab[0];}
  return 0;

}
//------------------------
Bool_t AliAnalysisTaskSELambdac::GetLambdacDaugh(AliAODMCParticle *part,TClonesArray *arrayMC) const{
  // check if the particle is a lambdac and if its decay mode is the correct one 
  Int_t numberOfLambdac=0;
  if(TMath::Abs(part->GetPdgCode())!=4122) return kFALSE;
  // Int_t daughTmp[2];
  // daughTmp[0]=part->GetDaughterLabel(0);
  // daughTmp[1]=part->GetDaughterLabel(1);
  Int_t nDaugh = (Int_t)part->GetNDaughters();
  if(nDaugh<2) return kFALSE;
  if(nDaugh>3) return kFALSE;
  AliAODMCParticle* pdaugh1 = (AliAODMCParticle*)arrayMC->At(part->GetDaughterLabel(0));
  if(!pdaugh1) {return kFALSE;}
  Int_t number1 = TMath::Abs(pdaugh1->GetPdgCode());
  AliAODMCParticle* pdaugh2 = (AliAODMCParticle*)arrayMC->At(part->GetDaughterLabel(1));
  if(!pdaugh2) {return kFALSE;}
  Int_t number2 = TMath::Abs(pdaugh2->GetPdgCode());

  if(nDaugh==3){
    Int_t thirdDaugh=part->GetDaughterLabel(1)-1;
    AliAODMCParticle* pdaugh3 = (AliAODMCParticle*)arrayMC->At(thirdDaugh);
    Int_t number3 = TMath::Abs(pdaugh3->GetPdgCode());
    if((number1==321 && number2==211 && number3==2212) ||
       (number1==211 && number2==321 && number3==2212) ||
       (number1==211 && number2==2212 && number3==321) ||
       (number1==321 && number2==2212 && number3==211) ||
       (number1==2212 && number2==321 && number3==211) ||
       (number1==2212 && number2==211 && number3==321)) {
      numberOfLambdac++;
    } 
  }

  if(nDaugh==2){

    //Lambda resonant
  
    //Lambda -> p K*0
    //
    Int_t nfiglieK=0;

    if((number1==2212 && number2==313)){
      nfiglieK=pdaugh2->GetNDaughters();
      if(nfiglieK!=2) return kFALSE;
      AliAODMCParticle* pdaughK1 = (AliAODMCParticle*)arrayMC->At(pdaugh2->GetDaughterLabel(0));
      AliAODMCParticle* pdaughK2 = (AliAODMCParticle*)arrayMC->At(pdaugh2->GetDaughterLabel(1));
      if(!pdaughK1) return kFALSE;
      if(!pdaughK2) return kFALSE;
      if((TMath::Abs(pdaughK1->GetPdgCode())==211 && TMath::Abs(pdaughK2->GetPdgCode())==321) || (TMath::Abs(pdaughK1->GetPdgCode())==321 && TMath::Abs(pdaughK2->GetPdgCode())==211)) numberOfLambdac++;
    }

    if((number1==313 && number2==2212)){
      nfiglieK=pdaugh1->GetNDaughters();
      if(nfiglieK!=2) return kFALSE;
      AliAODMCParticle* pdaughK1 = (AliAODMCParticle*)arrayMC->At(pdaugh1->GetDaughterLabel(0));
      AliAODMCParticle* pdaughK2 = (AliAODMCParticle*)arrayMC->At(pdaugh1->GetDaughterLabel(1));
      if(!pdaughK1) return kFALSE;
      if(!pdaughK2) return kFALSE;
      if((TMath::Abs(pdaughK1->GetPdgCode())==211 && TMath::Abs(pdaughK2->GetPdgCode())==321) || (TMath::Abs(pdaughK1->GetPdgCode())==321 && TMath::Abs(pdaughK2->GetPdgCode())==211)) numberOfLambdac++;
    }

    //Lambda -> Delta++ k
    Int_t nfiglieDelta=0;
    if(number1==321 && number2==2224){
      nfiglieDelta=pdaugh2->GetNDaughters();
      if(nfiglieDelta!=2) return kFALSE;
      AliAODMCParticle *pdaughD1=(AliAODMCParticle*)arrayMC->At(pdaugh2->GetDaughterLabel(0));
      AliAODMCParticle *pdaughD2=(AliAODMCParticle*)arrayMC->At(pdaugh2->GetDaughterLabel(1));
      if(!pdaughD1) return kFALSE;
      if(!pdaughD2) return kFALSE;
      if((TMath::Abs(pdaughD1->GetPdgCode())==211 && TMath::Abs(pdaughD2->GetPdgCode())==2212) || (TMath::Abs(pdaughD1->GetPdgCode())==2212 && TMath::Abs(pdaughD2->GetPdgCode())==211)) numberOfLambdac++;
    }
    if(number1==2224 && number2==321){
      nfiglieDelta=pdaugh1->GetNDaughters();
      if(nfiglieDelta!=2) return kFALSE;
      AliAODMCParticle* pdaughD1 = (AliAODMCParticle*)arrayMC->At(pdaugh1->GetDaughterLabel(0));
      AliAODMCParticle* pdaughD2 = (AliAODMCParticle*)arrayMC->At(pdaugh1->GetDaughterLabel(1));
      if(!pdaughD1) return kFALSE;
      if(!pdaughD2) return kFALSE;
      if((TMath::Abs(pdaughD1->GetPdgCode())==211 && TMath::Abs(pdaughD2->GetPdgCode())==2212) || (TMath::Abs(pdaughD1->GetPdgCode())==2212 && TMath::Abs(pdaughD2->GetPdgCode())==211)) numberOfLambdac++;
    }
    

    //Lambdac -> Lambda(1520) pi
    Int_t nfiglieLa=0;
    if(number1==3124 && number2==211){
      nfiglieLa=pdaugh1->GetNDaughters();
      if(nfiglieLa!=2) return kFALSE;
      AliAODMCParticle *pdaughL1=(AliAODMCParticle*)arrayMC->At(pdaugh1->GetDaughterLabel(0));
      AliAODMCParticle *pdaughL2=(AliAODMCParticle*)arrayMC->At(pdaugh1->GetDaughterLabel(1));
      if(!pdaughL1) return kFALSE;
      if(!pdaughL2) return kFALSE;
      if((TMath::Abs(pdaughL1->GetPdgCode())==321 && TMath::Abs(pdaughL2->GetPdgCode())==2212) || (TMath::Abs(pdaughL1->GetPdgCode())==2212 && TMath::Abs(pdaughL2->GetPdgCode())==321)) numberOfLambdac++;
    }
    if(number1==211 && number2==3124){
      nfiglieLa=pdaugh2->GetNDaughters();
      if(nfiglieLa!=2) return kFALSE;
      AliAODMCParticle *pdaughL1=(AliAODMCParticle*)arrayMC->At(pdaugh2->GetDaughterLabel(0));
      AliAODMCParticle *pdaughL2=(AliAODMCParticle*)arrayMC->At(pdaugh2->GetDaughterLabel(1));
      if(!pdaughL1) return kFALSE;
      if(!pdaughL2) return kFALSE;
      if((TMath::Abs(pdaughL1->GetPdgCode())==321 && TMath::Abs(pdaughL2->GetPdgCode())==2212) || (TMath::Abs(pdaughL1->GetPdgCode())==2212 && TMath::Abs(pdaughL2->GetPdgCode())==321)) numberOfLambdac++;
   
    }
  }

  if(numberOfLambdac>0) {return kTRUE;}
  return kFALSE;
}
//-----------------------------
Bool_t AliAnalysisTaskSELambdac::IspKpiMC(AliAODRecoDecayHF3Prong *d,TClonesArray *arrayMC) const{
  // Apply MC PID
  Int_t lab[3]={0,0,0},pdgs[3]={0,0,0};
  for(Int_t i=0;i<3;i++){
    AliAODTrack *daugh=(AliAODTrack*)d->GetDaughter(i);
    lab[i]=daugh->GetLabel();
    if(lab[i]<0) return kFALSE;
    AliAODMCParticle *part= (AliAODMCParticle*)arrayMC->At(lab[i]);
    if(!part) return kFALSE;
    pdgs[i]=TMath::Abs(part->GetPdgCode());
  }

  if(pdgs[0]==2212 && pdgs[1]==321 && pdgs[2]==211) return kTRUE;

  return kFALSE;
}
//-----------------------------
Bool_t AliAnalysisTaskSELambdac::IspiKpMC(AliAODRecoDecayHF3Prong *d,TClonesArray *arrayMC) const{

  // Apply MC PID
  Int_t lab[3]={0,0,0},pdgs[3]={0,0,0};
  for(Int_t i=0;i<3;i++){
    AliAODTrack *daugh=(AliAODTrack*)d->GetDaughter(i);
    lab[i]=daugh->GetLabel();
    if(lab[i]<0) return kFALSE;
    AliAODMCParticle *part= (AliAODMCParticle*)arrayMC->At(lab[i]);
    if(!part) return kFALSE;
    pdgs[i]=TMath::Abs(part->GetPdgCode());
  }

  if(pdgs[2]==2212 && pdgs[1]==321 && pdgs[0]==211) {return kTRUE;}

  return kFALSE;
}
//--------------------------------------
Bool_t AliAnalysisTaskSELambdac::VertexingKF(AliAODRecoDecayHF3Prong *d,Int_t *pdgs,Double_t field) const{
  // apply vertexing KF 
  Int_t iprongs[3]={0,1,2};
  Double_t mass[2]={0.,0.};
  Bool_t constraint=kFALSE;
  if(fCutsKF[0]>0.)constraint=kTRUE;
  //topological constr
  AliKFParticle *lambdac=d->ApplyVertexingKF(iprongs,3,pdgs,constraint,field,mass);
  if(!lambdac) return kFALSE;
  if(lambdac->GetChi2()/lambdac->GetNDF()>fCutsKF[1]) return kFALSE;
  return kTRUE;
}
//-------------------------------------
void AliAnalysisTaskSELambdac::IspiKpResonant(AliAODRecoDecayHF3Prong *d,Double_t field,Int_t *resNumber) const{
  
  // apply PID using the resonant channels 
  //if lambda* -> pk
  Double_t mass[2]={1.520,0.005};
  Int_t ipRes[2]={1,2};
  Int_t pdgres[2]={321,2212};
  AliKFParticle *lambda1520=d->ApplyVertexingKF(ipRes,2,pdgres,kFALSE,field,mass);
  Double_t probLa=TMath::Prob(lambda1520->GetChi2(),lambda1520->GetNDF());
  if(probLa>0.1) resNumber[0]=1;
  //K* -> kpi
  mass[0]=0.8961;mass[1]=0.03;
  ipRes[0]=0;ipRes[1]=1;
  pdgres[0]=211;pdgres[1]=321;
  AliKFParticle *kstar=d->ApplyVertexingKF(ipRes,2,pdgres,kFALSE,field,mass);
  Double_t probKa=TMath::Prob(kstar->GetChi2(),kstar->GetNDF());
  if(probKa>0.1) resNumber[1]=1;

 // Delta++
  mass[0]=1.232;mass[1]=0.15;
  ipRes[0]=0;ipRes[1]=2;
  pdgres[0]=211;pdgres[1]=2122;
  AliKFParticle *delta=d->ApplyVertexingKF(ipRes,2,pdgres,kFALSE,field,mass);
  Double_t probDe=TMath::Prob(delta->GetChi2(),delta->GetNDF());
  if(probDe>0.1) resNumber[2]=1;

  return ;

}
//-------------------------------------
void AliAnalysisTaskSELambdac::IspKpiResonant(AliAODRecoDecayHF3Prong *d,Double_t field,Int_t *resNumber) const{
   
  // apply PID using the resonant channels 
  //if lambda* -> pk
  Double_t mass[2]={1.520,0.005};
  Int_t ipRes[2]={0,1};
  Int_t pdgres[2]={2212,321};
  AliKFParticle *lambda1520=d->ApplyVertexingKF(ipRes,2,pdgres,kFALSE,field,mass);
  Double_t probLa=TMath::Prob(lambda1520->GetChi2(),lambda1520->GetNDF());
  if(probLa>0.1) resNumber[0]=1;
  //K* -> kpi
  mass[0]=0.8961;mass[1]=0.03;
  ipRes[0]=1;ipRes[1]=2;
  pdgres[1]=211;pdgres[0]=321;
  AliKFParticle *kstar=d->ApplyVertexingKF(ipRes,2,pdgres,kFALSE,field,mass);
  Double_t probKa=TMath::Prob(kstar->GetChi2(),kstar->GetNDF());
  if(probKa>0.1) resNumber[1]=1;

  //	Delta++
  mass[0]=1.232;mass[1]=0.15;
  ipRes[0]=0;ipRes[1]=2;
  pdgres[0]=2122;pdgres[1]=211;
  AliKFParticle *delta=d->ApplyVertexingKF(ipRes,2,pdgres,kFALSE,field,mass);
  Double_t probDe=TMath::Prob(delta->GetChi2(),delta->GetNDF());
  if(probDe>0.1) resNumber[2]=1;

  return ;

}
//---------------------------
void AliAnalysisTaskSELambdac::FillMassHists(AliAODEvent *aod,AliAODRecoDecayHF3Prong *part,
						TClonesArray *arrayMC, AliRDHFCutsLctopKpi *cuts,Int_t *nSelectedloose,Int_t *nSelectedtight)
{

  //if MC PID or no PID, unset PID
  //if(!fRealPid) cuts->SetUsePID(kFALSE);
  Int_t selection=cuts->IsSelected(part,AliRDHFCuts::kCandidate,aod);
  if(selection>0){
    nSelectedloose[0]=nSelectedloose[0]+1;
    Int_t iPtBin = -1;
    Double_t ptCand = part->Pt();
    Int_t index=0;

    for(Int_t ibin=0;ibin<fNPtBins&&iPtBin<0&&ptCand>fArrayBinLimits[0]&&ptCand<fArrayBinLimits[fNPtBins];ibin++){
      if(ptCand<fArrayBinLimits[ibin+1])iPtBin=ibin;
    }

    Float_t pdgCode=-2;
    Float_t pdgCode1=-1;
    Float_t pdgCode2=-1;
    Int_t labDp=-1;
    Float_t deltaPx=0.;
    Float_t deltaPy=0.;
    Float_t deltaPz=0.;
    Float_t truePt=0.;
    Float_t xDecay=0.;
    Float_t yDecay=0.;
    Float_t zDecay=0.;
    //    Bool_t IsInjected   = -1; 
    Bool_t IsLc		= 0;
    Bool_t IsLcfromLb	= 0;
    Bool_t IsLcfromc   = 0;
    Bool_t IsLcFromq = 0;

    if(fReadMC){
      //      AliAODMCHeader *mcHeader2 = (AliAODMCHeader*)aod->GetList()->FindObject(AliAODMCHeader::StdBranchName());
      AliVertexingHFUtils *util = new AliVertexingHFUtils();
      //      IsInjected = util->IsCandidateInjected(part,mcHeader2,arrayMC)?1:0;
   
			Int_t pdgCand =4122;
			Int_t pdgDaughter[3]={-1,-1,-1};
			pdgDaughter[0]=2212;
			pdgDaughter[1]=321;
			pdgDaughter[2]=211;   
			//labDp = MatchToMCLambdac(part,arrayMC);
      labDp = part->MatchToMC(pdgCand,arrayMC,3,pdgDaughter);
      if(labDp>=0){
	IsLc = 1;
	AliAODMCParticle *partDp = (AliAODMCParticle*)arrayMC->At(labDp);
	AliAODMCParticle *dg0 = (AliAODMCParticle*)arrayMC->At(partDp->GetDaughterLabel(0));
	AliAODMCParticle *dg1 = (AliAODMCParticle*)arrayMC->At(partDp->GetDaughterLabel(1));
	deltaPx=partDp->Px()-part->Px();
	deltaPy=partDp->Py()-part->Py();
	deltaPz=partDp->Pz()-part->Pz();
	truePt=partDp->Pt();
	xDecay=dg0->Xv();
	yDecay=dg0->Yv();
	zDecay=dg0->Zv();
	pdgCode=TMath::Abs(partDp->GetPdgCode());
	pdgCode1=TMath::Abs(dg0->GetPdgCode());
	pdgCode2=TMath::Abs(dg1->GetPdgCode());
        Int_t pdgMom=util->CheckOrigin(arrayMC,partDp,kFALSE);
          if(pdgMom == 5) IsLcfromLb =1;
          if(pdgMom == 4) IsLcfromc =1; 
        Int_t isThereaQuark=util->CheckOrigin(arrayMC,partDp,kTRUE);
	if(isThereaQuark>0) IsLcFromq = 1;
          
	
      }else{
	pdgCode=-1;
      }
      delete util;
    }

    Double_t invMasspKpi=-1.;
    Double_t invMasspiKp=-1.;
    Double_t invMasspKpiLpi=-1.;
    Double_t invMasspiKpLpi=-1.;
    Double_t invMasspKpiKp=-1.;
    Double_t invMasspiKpKp=-1.;
    Double_t invMasspiKpDk=-1.;
    Double_t invMasspKpiDk=-1.;
    Double_t invMasspiKp3Pr=-1.;
    Double_t invMasspKpi3Pr=-1.;
    Double_t field=aod->GetMagneticField();
    //apply MC PID
    if(fReadMC && fMCPid){
      if(IspKpiMC(part,arrayMC)) invMasspKpi=part->InvMassLcpKpi();
      if(IspiKpMC(part,arrayMC)) invMasspiKp=part->InvMassLcpiKp();
    }

    // apply realistic PID
    if(fRealPid){
      if(selection==1 || selection==3) invMasspKpi=part->InvMassLcpKpi();
      if(selection>=2) invMasspiKp=part->InvMassLcpiKp();
    }
    //apply PID using resonances
    if (fResPid && fRealPid) {
      Int_t ispKpi[3]={0,0,0};
      IspKpiResonant(part,field,ispKpi);
      Int_t ispiKp[3]={0,0,0};
      IspiKpResonant(part,field,ispiKp);
      if (selection==3 || selection==1) {
       if(ispKpi[0]==1){
        invMasspKpiLpi=part->InvMassLcpKpi();
       }
       if(ispKpi[1]==1){
        invMasspKpiKp=part->InvMassLcpKpi();
       }
       if(ispKpi[2]==1){
        invMasspKpiDk=part->InvMassLcpKpi();
       }
       if(ispKpi[2]==0 && ispKpi[1]==0 && ispKpi[0]==0){
        invMasspKpi3Pr=part->InvMassLcpKpi();
       }
      }
      if(selection>=2) {
       if(ispiKp[0]==1){
        invMasspiKpLpi=part->InvMassLcpiKp();
       }
       if(ispiKp[1]==1){
        invMasspiKpKp=part->InvMassLcpiKp();
       }
       if(ispiKp[2]==1){
        invMasspiKpDk=part->InvMassLcpiKp();
       }
       if(ispiKp[2]==0 && ispiKp[1]==0 && ispiKp[0]==0){
        invMasspiKp3Pr=part->InvMassLcpiKp();
       }
      }

     }
  

    if(invMasspiKp<0. && invMasspKpi<0.) return;
   
    Int_t passTightCuts=0;
    if(fAnalysis) {
      passTightCuts=fRDCutsAnalysis->IsSelected(part,AliRDHFCuts::kCandidate,aod);
      if(fUseKF){
	Int_t pdgs[3]={0,0,0};
	if(invMasspKpi>0.){
	  pdgs[0]=2212;pdgs[1]=321;pdgs[2]=211;
	  if(!VertexingKF(part,pdgs,field)) {
	    invMasspKpi=-1.;
	    invMasspKpi3Pr=-1.;
	    invMasspKpiDk=-1.;
	    invMasspKpiLpi=-1.;
	    invMasspKpiKp=-1.;
	  }
	}
	if(invMasspiKp>0.){
	  pdgs[0]=211;pdgs[1]=321;pdgs[2]=2212;
	  if(!VertexingKF(part,pdgs,field)) {
	    invMasspiKp=-1.;
	    invMasspiKp3Pr=-1.;
	    invMasspiKpDk=-1.;
	    invMasspiKpLpi=-1.;
	    invMasspiKpKp=-1.;
	  }
	}
      }
      if(passTightCuts>0) nSelectedtight[0]=nSelectedtight[0]+1;
    }
    // Eta and y plots:
    fhEta3ProngProd->Fill(part->Pt(),part->Eta());
    fhRap3ProngProd->Fill(part->Pt(),part->Y(4122));
    if(passTightCuts>0){
      fhEta3ProngAn->Fill(part->Pt(),part->Eta());
      fhRap3ProngAn->Fill(part->Pt(),part->Y(4122));
      if(invMasspiKp>0. && invMasspKpi>0.){
        fhMassLcPt->Fill(part->Pt(),invMasspKpi,0.5);
        fhMassLcPt->Fill(part->Pt(),invMasspiKp,0.5);
		if(part->Charge()==1){
		  fhMassLcplusPt->Fill(part->Pt(),invMasspKpi,0.5);
          fhMassLcplusPt->Fill(part->Pt(),invMasspiKp,0.5);
		}
		else if(part->Charge()==-1){
		  fhMassLcminusPt->Fill(part->Pt(),invMasspKpi,0.5);
          fhMassLcminusPt->Fill(part->Pt(),invMasspiKp,0.5);
		}
      }
      else if(invMasspiKp>0.){
        fhMassLcPt->Fill(part->Pt(),invMasspiKp);
        if(part->Charge()==1)       fhMassLcplusPt->Fill(part->Pt(),invMasspiKp);
	    else if(part->Charge()==-1) fhMassLcminusPt->Fill(part->Pt(),invMasspiKp);
	  }
	  else if(invMasspKpi>0.){
	    fhMassLcPt->Fill(part->Pt(),invMasspKpi);
        if(part->Charge()==1)       fhMassLcplusPt->Fill(part->Pt(),invMasspKpi);
	    else if(part->Charge()==-1) fhMassLcminusPt->Fill(part->Pt(),invMasspKpi);
	  } 
      
    }
    Float_t tmp[24];
    if (fFillNtuple) {
      tmp[0]=pdgCode;
      tmp[1]=deltaPx;
      tmp[2]=deltaPy;
      tmp[3]=deltaPz;
      tmp[4]=truePt;
      tmp[5]=xDecay;
      tmp[6]=yDecay;
      tmp[7]=zDecay;
      if(pdgCode1==2212) {tmp[8]=part->PtProng(0);}else{tmp[8]=0.;}
      if(pdgCode1==211) {tmp[10]=part->PtProng(0);}else{tmp[10]=0.;}
      tmp[9]=part->PtProng(1);
      if(pdgCode2==211) {tmp[10]=part->PtProng(2);}else{tmp[10]=0.;}
      tmp[11]=part->Pt();
      tmp[12]=part->CosPointingAngle();
      tmp[13]=part->DecayLength();
      tmp[14]=part->Xv();
      tmp[15]=part->Yv();
      tmp[16]=part->Zv();
      if(invMasspiKp>0.) tmp[17]=invMasspiKp;
      if(invMasspKpi>0.) tmp[17]=invMasspKpi;
      tmp[18]=part->GetSigmaVert();
      tmp[19]=part->Getd0Prong(0);
      tmp[20]=part->Getd0Prong(1);
      tmp[21]=part->Getd0Prong(2);
      tmp[22]=part->GetDCA();
      tmp[23]=part->Prodd0d0();
      fNtupleLambdac->Fill(tmp);
      PostData(7,fNtupleLambdac);
    }
    
    
    
    if(part->Pt()>3.&& part->Pt()<=6.){
      if(invMasspiKp>0. && invMasspKpi>0.){
	if(invMasspiKp>0.) fhMassPtGreater3->Fill(invMasspiKp,0.5);
	if(invMasspKpi>0.) fhMassPtGreater3->Fill(invMasspKpi,0.5);
      }else{
	if(invMasspiKp>0.) fhMassPtGreater3->Fill(invMasspiKp);
	if(invMasspKpi>0.) fhMassPtGreater3->Fill(invMasspKpi);
      }
      if(invMasspiKpLpi>0. && invMasspKpiLpi>0.){
	if(invMasspiKpLpi>0.) fhMassPtGreater3Lpi->Fill(invMasspiKpLpi,0.5);
	if(invMasspKpiLpi>0.) fhMassPtGreater3Lpi->Fill(invMasspKpiLpi,0.5);
      }else{
	if(invMasspiKpLpi>0.) fhMassPtGreater3Lpi->Fill(invMasspiKpLpi);
	if(invMasspKpiLpi>0.) fhMassPtGreater3Lpi->Fill(invMasspKpiLpi);
      }
      if(invMasspiKpKp>0. && invMasspKpiKp>0.){
	if(invMasspiKpKp>0.) fhMassPtGreater3Kp->Fill(invMasspiKpKp,0.5);
	if(invMasspKpiKp>0.) fhMassPtGreater3Kp->Fill(invMasspKpiKp,0.5);
      }else{
	if(invMasspiKpKp>0.) fhMassPtGreater3Kp->Fill(invMasspiKpKp);
	if(invMasspKpiKp>0.) fhMassPtGreater3Kp->Fill(invMasspKpiKp);
      }
      if(invMasspiKpDk>0. && invMasspKpiDk>0.){
       if(invMasspiKpDk>0.) fhMassPtGreater3Dk->Fill(invMasspiKpDk,0.5);
       if(invMasspKpiDk>0.) fhMassPtGreater3Dk->Fill(invMasspKpiDk,0.5);
      }else{
       if(invMasspiKpDk>0.) fhMassPtGreater3Dk->Fill(invMasspiKpDk);
       if(invMasspKpiDk>0.) fhMassPtGreater3Dk->Fill(invMasspKpiDk);
      }
      if(invMasspiKp3Pr>0. && invMasspKpi3Pr>0.){
        if(invMasspiKp3Pr>0.) fhMassPtGreater33Pr->Fill(invMasspiKp3Pr,0.5);
        if(invMasspKpi3Pr>0.) fhMassPtGreater33Pr->Fill(invMasspKpi3Pr,0.5);
      }else{
        if(invMasspiKp3Pr>0.) fhMassPtGreater33Pr->Fill(invMasspiKp3Pr);
        if(invMasspKpi3Pr>0.) fhMassPtGreater33Pr->Fill(invMasspKpi3Pr);
      }

      if(passTightCuts>0){
	if(invMasspiKp>0. && invMasspKpi>0.){
	  if(invMasspiKp>0.) fhMassPtGreater3TC->Fill(invMasspiKp,0.5);
	  if(invMasspKpi>0.) fhMassPtGreater3TC->Fill(invMasspKpi,0.5);
	}else{
	  if(invMasspiKp>0.) fhMassPtGreater3TC->Fill(invMasspiKp);
	  if(invMasspKpi>0.) fhMassPtGreater3TC->Fill(invMasspKpi);
	}
	if(invMasspiKpLpi>0. && invMasspKpiLpi>0.){
	  if(invMasspiKpLpi>0.) fhMassPtGreater3LpiTC->Fill(invMasspiKpLpi,0.5);
	  if(invMasspKpiLpi>0.) fhMassPtGreater3LpiTC->Fill(invMasspKpiLpi,0.5);
	}else{
	  if(invMasspiKpLpi>0.) fhMassPtGreater3LpiTC->Fill(invMasspiKpLpi);
	  if(invMasspKpiLpi>0.) fhMassPtGreater3LpiTC->Fill(invMasspKpiLpi);
	}
	if(invMasspiKpKp>0. && invMasspKpiKp>0.){
	  if(invMasspiKpKp>0.) fhMassPtGreater3KpTC->Fill(invMasspiKpKp,0.5);
	  if(invMasspKpiKp>0.) fhMassPtGreater3KpTC->Fill(invMasspKpiKp,0.5);
	}else{
	  if(invMasspiKpKp>0.) fhMassPtGreater3KpTC->Fill(invMasspiKpKp);
	  if(invMasspKpiKp>0.) fhMassPtGreater3KpTC->Fill(invMasspKpiKp);
	}
       
       if(invMasspiKpDk>0. && invMasspKpiDk>0.){
        if(invMasspiKpDk>0.) fhMassPtGreater3DkTC->Fill(invMasspiKpDk,0.5);
        if(invMasspKpiDk>0.) fhMassPtGreater3DkTC->Fill(invMasspKpiDk,0.5);
       }else{
        if(invMasspiKpDk>0.) fhMassPtGreater3DkTC->Fill(invMasspiKpDk);
        if(invMasspKpiDk>0.) fhMassPtGreater3DkTC->Fill(invMasspKpiDk);
       }
      if(invMasspiKp3Pr>0. && invMasspKpi3Pr>0.){
       if(invMasspiKp3Pr>0.) fhMassPtGreater33PrTC->Fill(invMasspiKp3Pr,0.5);
       if(invMasspKpi3Pr>0.) fhMassPtGreater33PrTC->Fill(invMasspKpi3Pr,0.5);
      }else{
       if(invMasspiKp3Pr>0.) fhMassPtGreater33PrTC->Fill(invMasspiKp3Pr);
       if(invMasspKpi3Pr>0.) fhMassPtGreater33PrTC->Fill(invMasspKpi3Pr);
      }

     }
    }
    if(part->Pt()>2.&& part->Pt()<=6.){
      if(invMasspiKp>0. && invMasspKpi>0.){
	if(invMasspiKp>0.) fhMassPtGreater2->Fill(invMasspiKp,0.5);
	if(invMasspKpi>0.) fhMassPtGreater2->Fill(invMasspKpi,0.5);
      }else{
	if(invMasspiKp>0.) fhMassPtGreater2->Fill(invMasspiKp);
	if(invMasspKpi>0.) fhMassPtGreater2->Fill(invMasspKpi);
      }
      if(invMasspiKpLpi>0. && invMasspKpiLpi>0.){
	if(invMasspiKpLpi>0.) fhMassPtGreater2Lpi->Fill(invMasspiKpLpi,0.5);
	if(invMasspKpiLpi>0.) fhMassPtGreater2Lpi->Fill(invMasspKpiLpi,0.5);
      }else{
	if(invMasspiKpLpi>0.) fhMassPtGreater2Lpi->Fill(invMasspiKpLpi);
	if(invMasspKpiLpi>0.) fhMassPtGreater2Lpi->Fill(invMasspKpiLpi);
      }
      if(invMasspiKpKp>0. && invMasspKpiKp>0.){
	if(invMasspiKpKp>0.) fhMassPtGreater2Kp->Fill(invMasspiKpKp,0.5);
	if(invMasspKpiKp>0.) fhMassPtGreater2Kp->Fill(invMasspKpiKp,0.5);
      }else{
	if(invMasspiKpKp>0.) fhMassPtGreater2Kp->Fill(invMasspiKpKp);
	if(invMasspKpiKp>0.) fhMassPtGreater2Kp->Fill(invMasspKpiKp);
      }
      if(invMasspiKpDk>0. && invMasspKpiDk>0.){
       if(invMasspiKpDk>0.) fhMassPtGreater2Dk->Fill(invMasspiKpDk,0.5);
       if(invMasspKpiDk>0.) fhMassPtGreater2Dk->Fill(invMasspKpiDk,0.5);
      }else{
       if(invMasspiKpDk>0.) fhMassPtGreater2Dk->Fill(invMasspiKpDk);
       if(invMasspKpiDk>0.) fhMassPtGreater2Dk->Fill(invMasspKpiDk);
      }
     if(invMasspiKp3Pr>0. && invMasspKpi3Pr>0.){
      if(invMasspiKp3Pr>0.) fhMassPtGreater23Pr->Fill(invMasspiKp3Pr,0.5);
      if(invMasspKpi3Pr>0.) fhMassPtGreater23Pr->Fill(invMasspKpi3Pr,0.5);
     }else{
      if(invMasspiKp3Pr>0.) fhMassPtGreater23Pr->Fill(invMasspiKp3Pr);
      if(invMasspKpi3Pr>0.) fhMassPtGreater23Pr->Fill(invMasspKpi3Pr);
    }

      if(passTightCuts>0){
	if(invMasspiKp>0. && invMasspKpi>0.){
	  if(invMasspiKp>0.) fhMassPtGreater2TC->Fill(invMasspiKp,0.5);
	  if(invMasspKpi>0.) fhMassPtGreater2TC->Fill(invMasspKpi,0.5);
	}else{
	  if(invMasspiKp>0.) fhMassPtGreater2TC->Fill(invMasspiKp);
	  if(invMasspKpi>0.) fhMassPtGreater2TC->Fill(invMasspKpi);
	}
	if(invMasspiKpLpi>0. && invMasspKpiLpi>0.){
	  if(invMasspiKpLpi>0.) fhMassPtGreater2LpiTC->Fill(invMasspiKpLpi,0.5);
	  if(invMasspKpiLpi>0.) fhMassPtGreater2LpiTC->Fill(invMasspKpiLpi,0.5);
	}else{
	  if(invMasspiKpLpi>0.) fhMassPtGreater2LpiTC->Fill(invMasspiKpLpi);
	  if(invMasspKpiLpi>0.) fhMassPtGreater2LpiTC->Fill(invMasspKpiLpi);
	}
	if(invMasspiKpKp>0. && invMasspKpiKp>0.){
	  if(invMasspiKpKp>0.) fhMassPtGreater2KpTC->Fill(invMasspiKpKp,0.5);
	  if(invMasspKpiKp>0.) fhMassPtGreater2KpTC->Fill(invMasspKpiKp,0.5);
	}else{
	  if(invMasspiKpKp>0.) fhMassPtGreater2KpTC->Fill(invMasspiKpKp);
	  if(invMasspKpiKp>0.) fhMassPtGreater2KpTC->Fill(invMasspKpiKp);
	}
       if(invMasspiKpDk>0. && invMasspKpiDk>0.){
        if(invMasspiKpDk>0.) fhMassPtGreater2DkTC->Fill(invMasspiKpDk,0.5);
        if(invMasspKpiDk>0.) fhMassPtGreater2DkTC->Fill(invMasspKpiDk,0.5);
       }else{
        if(invMasspiKpDk>0.) fhMassPtGreater2DkTC->Fill(invMasspiKpDk);
        if(invMasspKpiDk>0.) fhMassPtGreater2DkTC->Fill(invMasspKpiDk);
       }
       if(invMasspiKp3Pr>0. && invMasspKpi3Pr>0.){
        if(invMasspiKp3Pr>0.) fhMassPtGreater23PrTC->Fill(invMasspiKp3Pr,0.5);
        if(invMasspKpi3Pr>0.) fhMassPtGreater23PrTC->Fill(invMasspKpi3Pr,0.5);
      }else{
       if(invMasspiKp3Pr>0.) fhMassPtGreater23PrTC->Fill(invMasspiKp3Pr);
       if(invMasspKpi3Pr>0.) fhMassPtGreater23PrTC->Fill(invMasspKpi3Pr);
      }

    }
   }

    if(iPtBin>=0){
      index=GetHistoIndex(iPtBin);
      if(invMasspiKp>0. && invMasspKpi>0.){
	if(invMasspiKp>0.) fMassHist[index]->Fill(invMasspiKp,0.5);
	if(invMasspKpi>0.) fMassHist[index]->Fill(invMasspKpi,0.5);
      }else{
	if(invMasspiKp>0.) fMassHist[index]->Fill(invMasspiKp);
	if(invMasspKpi>0.) fMassHist[index]->Fill(invMasspKpi);
      }
      if(invMasspiKpLpi>0. && invMasspKpiLpi>0.){
	if(invMasspiKpLpi>0.) fMassHistLpi[index]->Fill(invMasspiKpLpi,0.5);
	if(invMasspKpiLpi>0.) fMassHistLpi[index]->Fill(invMasspKpiLpi,0.5);
      }else{
	if(invMasspiKpLpi>0.) fMassHistLpi[index]->Fill(invMasspiKpLpi);
	if(invMasspKpiLpi>0.) fMassHistLpi[index]->Fill(invMasspKpiLpi);
      }
      if(invMasspiKpKp>0. && invMasspKpiKp>0.){
	if(invMasspiKpKp>0.) fMassHistKp[index]->Fill(invMasspiKpKp,0.5);
	if(invMasspKpiKp>0.) fMassHistKp[index]->Fill(invMasspKpiKp,0.5);
      }else{
	if(invMasspiKpKp>0.) fMassHistKp[index]->Fill(invMasspiKpKp);
	if(invMasspKpiKp>0.) fMassHistKp[index]->Fill(invMasspKpiKp);
      }
      if(invMasspiKpDk>0. && invMasspKpiDk>0.){
       if(invMasspiKpDk>0.) fMassHistDk[index]->Fill(invMasspiKpDk,0.5);
       if(invMasspKpiDk>0.) fMassHistDk[index]->Fill(invMasspKpiDk,0.5);
      }else{
       if(invMasspiKpDk>0.) fMassHistDk[index]->Fill(invMasspiKpDk);
       if(invMasspKpiDk>0.) fMassHistDk[index]->Fill(invMasspKpiDk);
      }
      if(invMasspiKp3Pr>0. && invMasspKpi3Pr>0.){
       if(invMasspiKp3Pr>0.) fMassHist3Pr[index]->Fill(invMasspiKp3Pr,0.5);
       if(invMasspKpi3Pr>0.) fMassHist3Pr[index]->Fill(invMasspKpi3Pr,0.5);
      }else{
       if(invMasspiKp3Pr>0.) fMassHist3Pr[index]->Fill(invMasspiKp3Pr);
       if(invMasspKpi3Pr>0.) fMassHist3Pr[index]->Fill(invMasspKpi3Pr);
     }

      if(passTightCuts>0){
	if(invMasspiKp>0. && invMasspKpi>0. && passTightCuts==3){
	  if(invMasspiKp>0.) fMassHistTC[index]->Fill(invMasspiKp,0.5);
	  if(invMasspKpi>0.) fMassHistTC[index]->Fill(invMasspKpi,0.5);
	}else{
	  if(invMasspiKp>0. && passTightCuts==2) fMassHistTC[index]->Fill(invMasspiKp);
	  if(invMasspKpi>0. && passTightCuts==1) fMassHistTC[index]->Fill(invMasspKpi);
	}
	if(invMasspiKpLpi>0. && invMasspKpiLpi>0. && passTightCuts==3){
	  if(invMasspiKpLpi>0.) fMassHistLpiTC[index]->Fill(invMasspiKpLpi,0.5);
	  if(invMasspKpiLpi>0.) fMassHistLpiTC[index]->Fill(invMasspKpiLpi,0.5);
	}else{
	  if(invMasspiKpLpi>0. && passTightCuts==2) fMassHistLpiTC[index]->Fill(invMasspiKpLpi);
	  if(invMasspKpiLpi>0.&& passTightCuts==1) fMassHistLpiTC[index]->Fill(invMasspKpiLpi);
	}
	if(invMasspiKpKp>0. && invMasspKpiKp>0. && passTightCuts==3){
	  if(invMasspiKpKp>0.) fMassHistKpTC[index]->Fill(invMasspiKpKp,0.5);
	  if(invMasspKpiKp>0.) fMassHistKpTC[index]->Fill(invMasspKpiKp,0.5);
	}else{
	  if(invMasspiKpKp>0. && passTightCuts==2) fMassHistKpTC[index]->Fill(invMasspiKpKp);
	  if(invMasspKpiKp>0.&& passTightCuts==1) fMassHistKpTC[index]->Fill(invMasspKpiKp);
	}
       if(invMasspiKpDk>0. && invMasspKpiDk>0. && passTightCuts==3){
        if(invMasspiKpDk>0.) fMassHistDkTC[index]->Fill(invMasspiKpDk,0.5);
        if(invMasspKpiDk>0.) fMassHistDkTC[index]->Fill(invMasspKpiDk,0.5);
       }else{
        if(invMasspiKpDk>0. && passTightCuts==2) fMassHistDkTC[index]->Fill(invMasspiKpDk);
        if(invMasspKpiDk>0.&& passTightCuts==1) fMassHistDkTC[index]->Fill(invMasspKpiDk);
       }
       if(invMasspiKp3Pr>0. && invMasspKpi3Pr>0. && passTightCuts==3){
        if(invMasspiKp3Pr>0.) fMassHist3PrTC[index]->Fill(invMasspiKp3Pr,0.5);
        if(invMasspKpi3Pr>0.) fMassHist3PrTC[index]->Fill(invMasspKpi3Pr,0.5);
       }else{
        if(invMasspiKp3Pr>0. && passTightCuts==2) fMassHist3PrTC[index]->Fill(invMasspiKp3Pr);
        if(invMasspKpi3Pr>0.&& passTightCuts==1) fMassHist3PrTC[index]->Fill(invMasspKpi3Pr);
      }

      }

      if(fReadMC){
	if(labDp>=0) {
	  index=GetSignalHistoIndex(iPtBin);
	  if(invMasspiKp>0. && invMasspKpi>0.){
	    if(invMasspiKp>0.) fMassHist[index]->Fill(invMasspiKp,0.5);
	    if(invMasspKpi>0.) fMassHist[index]->Fill(invMasspKpi,0.5);
	  }else{
	    if(invMasspiKp>0.) fMassHist[index]->Fill(invMasspiKp);
	    if(invMasspKpi>0.) fMassHist[index]->Fill(invMasspKpi);
	  }
	  if(invMasspiKpLpi>0. && invMasspKpiLpi>0.){
	    if(invMasspiKpLpi>0.) fMassHistLpi[index]->Fill(invMasspiKpLpi,0.5);
	    if(invMasspKpiLpi>0.) fMassHistLpi[index]->Fill(invMasspKpiLpi,0.5);
	  }else{
	    if(invMasspiKpLpi>0.) fMassHistLpi[index]->Fill(invMasspiKpLpi);
	    if(invMasspKpiLpi>0.) fMassHistLpi[index]->Fill(invMasspKpiLpi);
	  }
	  if(invMasspiKpKp>0. && invMasspKpiKp>0.){
	    if(invMasspiKpKp>0.) fMassHistKp[index]->Fill(invMasspiKpKp,0.5);
	    if(invMasspKpiKp>0.) fMassHistKp[index]->Fill(invMasspKpiKp,0.5);
	  }else{
	    if(invMasspiKpKp>0.) fMassHistKp[index]->Fill(invMasspiKpKp);
	    if(invMasspKpiKp>0.) fMassHistKp[index]->Fill(invMasspKpiKp);
	  }
	  if(invMasspiKpDk>0. && invMasspKpiDk>0.){
	    if(invMasspiKpDk>0.) fMassHistDk[index]->Fill(invMasspiKpDk,0.5);
	    if(invMasspKpiDk>0.) fMassHistDk[index]->Fill(invMasspKpiDk,0.5);
	  }else{
	    if(invMasspiKpDk>0.) fMassHistDk[index]->Fill(invMasspiKpDk);
	    if(invMasspKpiDk>0.) fMassHistDk[index]->Fill(invMasspKpiDk);
	  }
	  if(invMasspiKp3Pr>0. && invMasspKpi3Pr>0.){
	    if(invMasspiKp3Pr>0.) fMassHist3Pr[index]->Fill(invMasspiKp3Pr,0.5);
	    if(invMasspKpi3Pr>0.) fMassHist3Pr[index]->Fill(invMasspKpi3Pr,0.5);
	  }else{
	    if(invMasspiKp3Pr>0.) fMassHist3Pr[index]->Fill(invMasspiKp3Pr);
	    if(invMasspKpi3Pr>0.) fMassHistDk[index]->Fill(invMasspKpi3Pr);
	  }

	  if(passTightCuts>0){
	    if(invMasspiKp>0. && invMasspKpi>0. && passTightCuts==3){
	      if(invMasspiKp>0.) fMassHistTC[index]->Fill(invMasspiKp,0.5);
	      if(invMasspKpi>0.) fMassHistTC[index]->Fill(invMasspKpi,0.5);
	    }else{
	      if(invMasspiKp>0. && passTightCuts==2) fMassHistTC[index]->Fill(invMasspiKp);
	      if(invMasspKpi>0.&& passTightCuts==1) fMassHistTC[index]->Fill(invMasspKpi);
	    }
	    if(invMasspiKpLpi>0. && invMasspKpiLpi>0. && passTightCuts==3){
	      if(invMasspiKpLpi>0.) fMassHistLpiTC[index]->Fill(invMasspiKpLpi,0.5);
	      if(invMasspKpiLpi>0.) fMassHistLpiTC[index]->Fill(invMasspKpiLpi,0.5);
	    }else{
	      if(invMasspiKpLpi>0. && passTightCuts==2) fMassHistLpiTC[index]->Fill(invMasspiKpLpi);
	      if(invMasspKpiLpi>0.&& passTightCuts==1) fMassHistLpiTC[index]->Fill(invMasspKpiLpi);
	    } 
	    if(invMasspiKpKp>0. && invMasspKpiKp>0. && passTightCuts==3){
	      if(invMasspiKpKp>0.) fMassHistKpTC[index]->Fill(invMasspiKpKp,0.5);
	      if(invMasspKpiKp>0.) fMassHistKpTC[index]->Fill(invMasspKpiKp,0.5);
	    }else{
	      if(invMasspiKpKp>0. && passTightCuts==2) fMassHistKpTC[index]->Fill(invMasspiKpKp);
	      if(invMasspKpiKp>0.&& passTightCuts==1) fMassHistKpTC[index]->Fill(invMasspKpiKp);
	    }
	    if(invMasspiKpDk>0. && invMasspKpiDk>0. && passTightCuts==3){
	      if(invMasspiKpDk>0.) fMassHistDkTC[index]->Fill(invMasspiKpDk,0.5);
	      if(invMasspKpiDk>0.) fMassHistDkTC[index]->Fill(invMasspKpiDk,0.5);
	    }else{
	      if(invMasspiKpDk>0. && passTightCuts==2) fMassHistDkTC[index]->Fill(invMasspiKpDk);
	      if(invMasspKpiDk>0.&& passTightCuts==1) fMassHistDkTC[index]->Fill(invMasspKpiDk);
	    }
	    if(invMasspiKp3Pr>0. && invMasspKpi3Pr>0. && passTightCuts==3){
	      if(invMasspiKp3Pr>0.) fMassHist3PrTC[index]->Fill(invMasspiKp3Pr,0.5);
	      if(invMasspKpi3Pr>0.) fMassHist3PrTC[index]->Fill(invMasspKpi3Pr,0.5);
	    }else{
	      if(invMasspiKp3Pr>0. && passTightCuts==2) fMassHist3PrTC[index]->Fill(invMasspiKp3Pr);
	      if(invMasspKpi3Pr>0.&& passTightCuts==1) fMassHist3PrTC[index]->Fill(invMasspKpi3Pr);
	    }
	  }      
    
         if(IsLc && IsLcfromLb){
          index=GetLbHistoIndex(iPtBin);
          if(invMasspiKp>0. && invMasspKpi>0.){
	    	if(invMasspiKp>0.) fMassHist[index]->Fill(invMasspiKp,0.5);
	        if(invMasspKpi>0.) fMassHist[index]->Fill(invMasspKpi,0.5);
	  }else{
	        if(invMasspiKp>0.) fMassHist[index]->Fill(invMasspiKp);
	        if(invMasspKpi>0.) fMassHist[index]->Fill(invMasspKpi);
	  } 
          if(passTightCuts>0){ 	
                if(invMasspiKp>0. && invMasspKpi>0. && passTightCuts==3){
        		if(invMasspiKp>0.) fMassHistTC[index]->Fill(invMasspiKp,0.5);
       			if(invMasspKpi>0.) fMassHistTC[index]->Fill(invMasspKpi,0.5);
	    	}else{
	       	       	if(invMasspiKp>0. && passTightCuts==2) fMassHistTC[index]->Fill(invMasspiKp); 
                        if(invMasspKpi>0.&& passTightCuts==1) fMassHistTC[index]->Fill(invMasspKpi);
                }
	  }
	 }
        if(IsLc && !IsLcfromLb && IsLcfromc)  {
           index=GetcOnlyHistoIndex(iPtBin);
           if(invMasspiKp>0. && invMasspKpi>0.){
	         	if(invMasspiKp>0.) fMassHist[index]->Fill(invMasspiKp,0.5);
	   	        if(invMasspKpi>0.) fMassHist[index]->Fill(invMasspKpi,0.5);
	   }else{
	         	    if(invMasspiKp>0.) fMassHist[index]->Fill(invMasspiKp);
			    if(invMasspKpi>0.) fMassHist[index]->Fill(invMasspKpi);
	   } 
           if(passTightCuts>0){
	       if(invMasspiKp>0. && invMasspKpi>0. && passTightCuts==3){
	       if(invMasspiKp>0.) fMassHistTC[index]->Fill(invMasspiKp,0.5);
	       if(invMasspKpi>0.) fMassHistTC[index]->Fill(invMasspKpi,0.5);
	    }else{
	       	if(invMasspiKp>0. && passTightCuts==2) fMassHistTC[index]->Fill(invMasspiKp);
	    	if(invMasspKpi>0.&& passTightCuts==1) fMassHistTC[index]->Fill(invMasspKpi);
	    }
	   }

	}

         if(IsLc && !IsLcFromq)  {
           index=GetNoQuarkHistoIndex(iPtBin);
           if(invMasspiKp>0. && invMasspKpi>0.){
	         	if(invMasspiKp>0.) fMassHist[index]->Fill(invMasspiKp,0.5);
	   	        if(invMasspKpi>0.) fMassHist[index]->Fill(invMasspKpi,0.5);
	   }else{
	         	    if(invMasspiKp>0.) fMassHist[index]->Fill(invMasspiKp);
			    if(invMasspKpi>0.) fMassHist[index]->Fill(invMasspKpi);
	   } 
           if(passTightCuts>0){
	       if(invMasspiKp>0. && invMasspKpi>0. && passTightCuts==3){
	       if(invMasspiKp>0.) fMassHistTC[index]->Fill(invMasspiKp,0.5);
	       if(invMasspKpi>0.) fMassHistTC[index]->Fill(invMasspKpi,0.5);
	    }else{
	       	if(invMasspiKp>0. && passTightCuts==2) fMassHistTC[index]->Fill(invMasspiKp);
	    	if(invMasspKpi>0.&& passTightCuts==1) fMassHistTC[index]->Fill(invMasspKpi);
	    }
	   }

	 }
     
	}else{
	  index=GetBackgroundHistoIndex(iPtBin);
	  if(invMasspiKp>0. && invMasspKpi>0.){
	    fMassHist[index]->Fill(invMasspiKp,0.5);
	    fMassHist[index]->Fill(invMasspKpi,0.5);
	  }else{
	    if(invMasspiKp>0.) fMassHist[index]->Fill(invMasspiKp);
	    if(invMasspKpi>0.) fMassHist[index]->Fill(invMasspKpi);
	  }
	  if(invMasspiKpLpi>0. && invMasspKpiLpi>0.){
	    if(invMasspiKpLpi>0.) fMassHistLpi[index]->Fill(invMasspiKpLpi,0.5);
	    if(invMasspKpiLpi>0.) fMassHistLpi[index]->Fill(invMasspKpiLpi,0.5);
	  }else{
	    if(invMasspiKpLpi>0.) fMassHistLpi[index]->Fill(invMasspiKpLpi);
	    if(invMasspKpiLpi>0.) fMassHistLpi[index]->Fill(invMasspKpiLpi);
	  }
	  if(invMasspiKpKp>0. && invMasspKpiKp>0.){
	    if(invMasspiKpKp>0.) fMassHistKp[index]->Fill(invMasspiKpKp,0.5);
	    if(invMasspKpiKp>0.) fMassHistKp[index]->Fill(invMasspKpiKp,0.5);
	  }else{
	    if(invMasspiKpKp>0.) fMassHistKp[index]->Fill(invMasspiKpKp);
	    if(invMasspKpiKp>0.) fMassHistKp[index]->Fill(invMasspKpiKp);
	  }
	  if(invMasspiKpDk>0. && invMasspKpiDk>0.){
	    if(invMasspiKpDk>0.) fMassHistDk[index]->Fill(invMasspiKpDk,0.5);
	    if(invMasspKpiDk>0.) fMassHistDk[index]->Fill(invMasspKpiDk,0.5);
	  }else{
	    if(invMasspiKpDk>0.) fMassHistDk[index]->Fill(invMasspiKpDk);
	    if(invMasspKpiDk>0.) fMassHistDk[index]->Fill(invMasspKpiDk);
	  }
	  if(invMasspiKp3Pr>0. && invMasspKpi3Pr>0.){
	    if(invMasspiKp3Pr>0.) fMassHist3Pr[index]->Fill(invMasspiKp3Pr,0.5);
	    if(invMasspKpi3Pr>0.) fMassHist3Pr[index]->Fill(invMasspKpi3Pr,0.5);
	  }else{
	    if(invMasspiKp3Pr>0.) fMassHist3Pr[index]->Fill(invMasspiKp3Pr);
	    if(invMasspKpi3Pr>0.) fMassHistDk[index]->Fill(invMasspKpi3Pr);
	  }
	  if(invMasspiKp>0. && invMasspKpi>0. && passTightCuts==3){
	    fMassHistTC[index]->Fill(invMasspiKp,0.5);
	    fMassHistTC[index]->Fill(invMasspKpi,0.5);
	  }else{
	    if(invMasspiKp>0. && passTightCuts==2) fMassHistTC[index]->Fill(invMasspiKp);
	    if(invMasspKpi>0. && passTightCuts==1) fMassHistTC[index]->Fill(invMasspKpi);
	  }
	  if(invMasspiKpLpi>0. && invMasspKpiLpi>0. && passTightCuts==3){
	    if(invMasspiKpLpi>0.) fMassHistLpiTC[index]->Fill(invMasspiKpLpi,0.5);
	    if(invMasspKpiLpi>0.) fMassHistLpiTC[index]->Fill(invMasspKpiLpi,0.5);
	  }else{
	    if(invMasspiKpLpi>0. && passTightCuts==2) fMassHistLpiTC[index]->Fill(invMasspiKpLpi);
	    if(invMasspKpiLpi>0.&& passTightCuts==1) fMassHistLpiTC[index]->Fill(invMasspKpiLpi);
	  } 
	  if(invMasspiKpKp>0. && invMasspKpiKp>0. && passTightCuts==3){
	    if(invMasspiKpKp>0.) fMassHistKpTC[index]->Fill(invMasspiKpKp,0.5);
	    if(invMasspKpiKp>0.) fMassHistKpTC[index]->Fill(invMasspKpiKp,0.5);
	  }else{
	    if(invMasspiKpKp>0. && passTightCuts==2) fMassHistKpTC[index]->Fill(invMasspiKpKp);
	    if(invMasspKpiKp>0.&& passTightCuts==1) fMassHistKpTC[index]->Fill(invMasspKpiKp);
	  }
	  if(invMasspiKpDk>0. && invMasspKpiDk>0. && passTightCuts==3){
	    if(invMasspiKpDk>0.) fMassHistDkTC[index]->Fill(invMasspiKpDk,0.5);
	    if(invMasspKpiDk>0.) fMassHistDkTC[index]->Fill(invMasspKpiDk,0.5);
	  }else{
	    if(invMasspiKpDk>0. && passTightCuts==2) fMassHistDkTC[index]->Fill(invMasspiKpDk);
	    if(invMasspKpiDk>0.&& passTightCuts==1) fMassHistDkTC[index]->Fill(invMasspKpiDk);
	   }
	   if(invMasspiKp3Pr>0. && invMasspKpi3Pr>0. && passTightCuts==3){
	    if(invMasspiKp3Pr>0.) fMassHist3PrTC[index]->Fill(invMasspiKp3Pr,0.5);
	    if(invMasspKpi3Pr>0.) fMassHist3PrTC[index]->Fill(invMasspKpi3Pr,0.5);
	   }else{
	    if(invMasspiKp3Pr>0. && passTightCuts==2) fMassHist3PrTC[index]->Fill(invMasspiKp3Pr);
	    if(invMasspKpi3Pr>0.&& passTightCuts==1) fMassHist3PrTC[index]->Fill(invMasspKpi3Pr);
	  }
	}

      }
    }
  }
  return;
}
//-----------------------
void AliAnalysisTaskSELambdac::FillVarHists(AliAODRecoDecayHF3Prong *part,
					       TClonesArray *arrMC,
					       AliRDHFCutsLctopKpi *cuts,
					       /* TList *listout,*/
					       AliAODEvent* aod) {
  //
  // function used in UserExec to fill variable histograms analysing MC
  //

  TString fillthis="";

  Double_t mPDG=TDatabasePDG::Instance()->GetParticle(4122)->Mass();
  Double_t invmasscut=0.05;

  Int_t pdgDgLctopKpi[3]={2212,321,211};
  Int_t lab=-9999;
  //  Bool_t IsLcfromLb=0;
  Bool_t IsLcfromc=0;
  if(fReadMC){
    lab=part->MatchToMC(4122,arrMC,3,pdgDgLctopKpi); //return MC particle label if the array corresponds to a Lc, -1 if not (cf. AliAODRecoDecay.cxx)
    if(lab>=0){
      AliAODMCParticle *partDp = (AliAODMCParticle*)arrMC->At(lab);
      AliVertexingHFUtils *util = new AliVertexingHFUtils();
      Int_t pdgMom=util->CheckOrigin(arrMC,partDp,kFALSE);
      //      if(pdgMom == 5) IsLcfromLb =1;
      if(pdgMom == 4) IsLcfromc =1; 
      delete util;
    }
  }
 // Int_t isSelectedPID=cuts->IsSelectedPID(part); // 0 rejected, 1 Lc -> p K- pi+ (K at center because different sign is there),
                                                 // 2 Lc -> pi+ K- p (K at center because different sign is there),
                                                 // 3 both (it should never happen...)
  Int_t isSelectedPID=cuts->IsSelected(part,AliRDHFCuts::kPID,aod); // 0 rejected, 1 Lc -> p K- pi+ (K at center because different sign is there),
                                                 // 2 Lc -> pi+ K- p (K at center because different sign is there),
                                                 // 3 both (it should never happen...)

  if (isSelectedPID==0)fNentries->Fill(7);
  if (isSelectedPID==1)fNentries->Fill(8);
  if (isSelectedPID==2)fNentries->Fill(9);
  if (isSelectedPID==3)fNentries->Fill(10);

  AliAODTrack *prong0=(AliAODTrack*)part->GetDaughter(0);
  AliAODTrack *prong1=(AliAODTrack*)part->GetDaughter(1);
  AliAODTrack *prong2=(AliAODTrack*)part->GetDaughter(2);
  if (!prong0 || !prong1 || !prong2) {
    fNentries->Fill(5);
    return;
  }

  Double_t minvLcpKpi = part->InvMassLcpKpi();
  Double_t minvLcpiKp = part->InvMassLcpiKp();


  //check pdg of the prongs
  Int_t labprong[3]={-1,-1,-1};
  if(fReadMC){
    labprong[0]=prong0->GetLabel();
    labprong[1]=prong1->GetLabel();
    labprong[2]=prong2->GetLabel();
  }

  AliAODMCParticle *mcprong=0;
  Int_t pdgProngMC[3]={-1,-1,-1};
  if(fReadMC) {
    for (Int_t iprong=0;iprong<3;iprong++){
      if(labprong[iprong]<0) continue;
      mcprong = (AliAODMCParticle*)arrMC->At(labprong[iprong]);
      pdgProngMC[iprong]=TMath::Abs(mcprong->GetPdgCode());
    }
  }

  Int_t pdgProngPID[3]={-1,-1,-1};
  if(isSelectedPID>0){
    pdgProngPID[1]=321;
    if(isSelectedPID==1) {pdgProngPID[0]=2212;pdgProngPID[2]=211;}
    if(isSelectedPID==2) {pdgProngPID[0]=211;pdgProngPID[2]=2212;}
  }

  //fill hRealTot ---> PID efficiency denominators
  //cuts->SetUsePID(kFALSE); //PAOLA
  Int_t selectionTrack = cuts->IsSelected(part,AliRDHFCuts::kTracks,aod);
  
  //IsInjected check	 
  Bool_t IsInjected = 0;
  if(fReadMC){
    AliAODMCHeader *mcHeader2 = (AliAODMCHeader*)aod->GetList()->FindObject(AliAODMCHeader::StdBranchName());
    AliVertexingHFUtils *util = new AliVertexingHFUtils();
    IsInjected = util->IsCandidateInjected(part,mcHeader2,arrMC)?1:0;
    delete util;
  }
  if(fReadMC && selectionTrack>0) { // 3prongs candidate x Lc (only track selection) Jaime
    Int_t isReal=0;
    if(lab>=0 &&  IsLcfromc){ // Signal
      for (Int_t iprong=0; iprong<3; iprong++) {
	switch (pdgProngMC[iprong]) {
	case 2212:
	  fillthis="hpRealTot";
	  ((TH1F*)fOutputMC->FindObject(fillthis))->Fill(part->PtProng(iprong));
	  isReal++;
	  break;
	case 321:
	  fillthis="hKRealTot";
	  ((TH1F*)fOutputMC->FindObject(fillthis))->Fill(part->PtProng(iprong));
	  isReal++;
	  break;
	case 211:
	  fillthis="hpiRealTot";
	  ((TH1F*)fOutputMC->FindObject(fillthis))->Fill(part->PtProng(iprong));
	  isReal++;
	  break;
	default:
	  break;
	}
      }
	fillthis="hLcRealTot";
	if(isReal==3) ((TH1F*)fOutputMC->FindObject(fillthis))->Fill(part->Pt());
    //Marcel: Should we zero isReal
    } else if(!IsInjected) { // bkg  // bkg
      for (Int_t iprong=0; iprong<3; iprong++) {
	switch (pdgProngMC[iprong]) {
	case 2212:
	  fillthis="hbpRealTot";
	  ((TH1F*)fOutputMC->FindObject(fillthis))->Fill(part->PtProng(iprong));
	  isReal++;
	  break;
	case 321:
	  fillthis="hbKRealTot";
	  ((TH1F*)fOutputMC->FindObject(fillthis))->Fill(part->PtProng(iprong));
	  isReal++;
	  break;
	case 211:
	  fillthis="hbpiRealTot";
	  ((TH1F*)fOutputMC->FindObject(fillthis))->Fill(part->PtProng(iprong));
	  isReal++;
	  break;
	default:
	  break;
	}
      }
    fillthis="hbLcRealTot";
	if(isReal==3) ((TH1F*)fOutputMC->FindObject(fillthis))->Fill(part->Pt()); 
    }
  }


  //cuts->SetUsePID(kTRUE); //PAOLA
  Int_t selection=cuts->IsSelected(part,AliRDHFCuts::kCandidate,aod);

  if ( (lab>=0 && fReadMC && IsLcfromc) ||             // signal X MC
       (isSelectedPID>0 && !fReadMC) ) { // signal+bkg X real data

    fillthis="hMass";
    if ( (fReadMC && ((AliAODMCParticle*)arrMC->At(lab))->GetPdgCode() == 4122) ||
	 (!fReadMC && (isSelectedPID>0 && part->Charge()>0)) ) { // Lc
      if ( (pdgProngPID[0]==2212) && (pdgProngPID[1]==321) && (pdgProngPID[2]==211) )
	((TH1F*)fOutputMC->FindObject(fillthis))->Fill(minvLcpKpi); // signal+bkg X MC and real data
      else if ( (pdgProngPID[0]==211) && (pdgProngPID[1]==321) && (pdgProngPID[2]==2212) )
	((TH1F*)fOutputMC->FindObject(fillthis))->Fill(minvLcpiKp); // signal+bkg X MC and real data
    }
    else if ( (fReadMC && ((AliAODMCParticle*)arrMC->At(lab))->GetPdgCode() == -4122) ||
	      (!fReadMC && (isSelectedPID>0 && part->Charge()<0)) ) { // anti-Lc
      if ( (pdgProngPID[0]==2212) && (pdgProngPID[1]==321) && (pdgProngPID[2]==211) )
	((TH1F*)fOutputMC->FindObject(fillthis))->Fill(minvLcpKpi); // signal+bkg X MC and real data
      else if ( (pdgProngPID[0]==211) && (pdgProngPID[1]==321) && (pdgProngPID[2]==2212) )
	((TH1F*)fOutputMC->FindObject(fillthis))->Fill(minvLcpiKp); // signal+bkg X MC and real data
    }

    if (selection>0) { // 3prongs candidate x Lc (yes PID)

      Double_t ptmax=0.;
      Int_t isID=0;
      Int_t isCorrect=0;
      for (Int_t iprong=0; iprong<3; iprong++) {
	if (part->PtProng(iprong)>ptmax) ptmax=part->PtProng(iprong);

	AliAODTrack *prong = (AliAODTrack*)part->GetDaughter(iprong);
	AliAODPid *pidObjtrk = (AliAODPid*)prong->GetDetPid();
	AliAODPidHF *pidObj = (AliAODPidHF*)cuts->GetPidHF();
	Bool_t hasTOF=pidObj->CheckStatus(prong,"TOF");
	Bool_t hasTPC=pidObj->CheckStatus(prong,"TPC");
	Double_t tofSignal=0.;
	Double_t dedxTPC=0.;
	Double_t momTOF=0.;
	Double_t momTPC=0.;
	if(hasTOF) {
	  momTOF = prong->P();
	  tofSignal=pidObjtrk->GetTOFsignal();
	}
	if(hasTPC) {
	  momTPC = pidObjtrk->GetTPCmomentum();
	  dedxTPC=pidObjtrk->GetTPCsignal();
	}
	switch (pdgProngPID[iprong]) {
	case 2212:
	  fillthis="hpTOFSignal";
	  ((TH1F*)fOutputMC->FindObject(fillthis))->Fill(tofSignal);
	  fillthis="hpTPCSignal";
	  ((TH1F*)fOutputMC->FindObject(fillthis))->Fill(dedxTPC);
	  fillthis="hpptProng";
	  ((TH1F*)fOutputMC->FindObject(fillthis))->Fill(part->PtProng(iprong));
          fillthis="hpptProngLcPt";
	  ((TH1F*)fOutputMC->FindObject(fillthis))->Fill(part->Pt(),part->PtProng(iprong));
	  fillthis="hpd0Prong";
	  ((TH1F*)fOutputMC->FindObject(fillthis))->Fill(part->Getd0Prong(iprong));
          fillthis="hpd0ProngLcPt";
	  ((TH1F*)fOutputMC->FindObject(fillthis))->Fill(part->Pt(),part->Getd0Prong(iprong));
	  fillthis="hpSignalVspTPC";
	  ((TH2F*)fOutputMC->FindObject(fillthis))->Fill(momTPC,dedxTPC);
	  fillthis="hpSignalVspTOF";
	  ((TH2F*)fOutputMC->FindObject(fillthis))->Fill(momTOF,tofSignal);
	  AliPID::EParticleType typep;
	  typep=AliPID::EParticleType(4);
	  if(hasTPC) {
	    //Double_t nsigmap = fUtilPid->NumberOfSigmasTPC(prong,typep);
	    Double_t nsigmap = fPIDResponse->NumberOfSigmasTPC(prong,typep);
	    fillthis="hpSigmaVspTPC";
	    ((TH2F*)fOutputMC->FindObject(fillthis))->Fill(momTPC,nsigmap);
	  }
	  if(hasTOF){
	    //Double_t nsigma=fUtilPid->NumberOfSigmasTOF(prong,typep);
	    Double_t nsigma=fPIDResponse->NumberOfSigmasTOF(prong,typep);
	    fillthis="hpSigmaVspTOF";
	    ((TH2F*)fOutputMC->FindObject(fillthis))->Fill(momTOF,nsigma);
	  }

	  if (fReadMC) { // fill hpIDTot ---> PID contamination denominator
	                 //      hIDGood, hnoID ---> PID numerators
	    fillthis="hpIDTot";
	    ((TH1F*)fOutputMC->FindObject(fillthis))->Fill(part->PtProng(iprong));
	    isID++;
	    if(pdgProngMC[iprong]==2212) { // id protons
	      fillthis="hpIDGood";
	      ((TH1F*)fOutputMC->FindObject(fillthis))->Fill(part->PtProng(iprong));
	      isCorrect++;  
	    }
	    else { // misidentified electrons, muons, pions and kaons
	      fillthis="hnopIDp";
	      ((TH1F*)fOutputMC->FindObject(fillthis))->Fill(part->PtProng(iprong));
	    }
	  }
	  break;

	case 321:
	  fillthis="hKTOFSignal";
	  ((TH1F*)fOutputMC->FindObject(fillthis))->Fill(tofSignal);
	  fillthis="hKTPCSignal";
	  ((TH1F*)fOutputMC->FindObject(fillthis))->Fill(dedxTPC);
	  fillthis="hKptProng";
	  ((TH1F*)fOutputMC->FindObject(fillthis))->Fill(part->PtProng(iprong));
          fillthis="hKptProngLcPt";
	  ((TH1F*)fOutputMC->FindObject(fillthis))->Fill(part->Pt(),part->PtProng(iprong));
	  fillthis="hKd0Prong";
	  ((TH1F*)fOutputMC->FindObject(fillthis))->Fill(part->Getd0Prong(iprong));
          fillthis="hKd0ProngLcPt";
	  ((TH1F*)fOutputMC->FindObject(fillthis))->Fill(part->Pt(),part->Getd0Prong(iprong));
	  fillthis="hKSignalVspTPC";
	  ((TH2F*)fOutputMC->FindObject(fillthis))->Fill(momTPC,dedxTPC);
	  fillthis="hKSignalVspTOF";
	  ((TH2F*)fOutputMC->FindObject(fillthis))->Fill(momTOF,tofSignal);
	  AliPID::EParticleType typek;
	  typek=AliPID::EParticleType(3);
	  if(hasTPC) {
	    //Double_t nsigmap = fUtilPid->NumberOfSigmasTPC(prong,typek);
	    Double_t nsigmap = fPIDResponse->NumberOfSigmasTPC(prong,typek);
	    fillthis="hKSigmaVspTPC";
	    ((TH2F*)fOutputMC->FindObject(fillthis))->Fill(momTPC,nsigmap);
	  }
	  if(hasTOF){
	    //Double_t nsigma = fUtilPid->NumberOfSigmasTOF(prong,typek);
	    Double_t nsigma=fPIDResponse->NumberOfSigmasTOF(prong,typek); // protection against 'old' AODs...
	    fillthis="hKSigmaVspTOF";
	    ((TH2F*)fOutputMC->FindObject(fillthis))->Fill(momTOF,nsigma);
	  }

	  if (fReadMC) { // fill hKIDTot ---> PID contamination denominator
	                 //      hIDGood, hnoID ---> PID numerators
	    fillthis="hKIDTot";
	    ((TH1F*)fOutputMC->FindObject(fillthis))->Fill(part->PtProng(iprong));
	    isID++;
	    if(pdgProngMC[iprong]==321) { // id kaons
	      fillthis="hKIDGood";
	      ((TH1F*)fOutputMC->FindObject(fillthis))->Fill(part->PtProng(iprong));
	      isCorrect++;
	    }
	    else { // misidentified electrons, muons, pions and protons
	      fillthis="hnokIDk";
	      ((TH1F*)fOutputMC->FindObject(fillthis))->Fill(part->PtProng(iprong));
	    }
	  }
	  break;

	case 211:
	  fillthis="hpiTOFSignal";
	  ((TH1F*)fOutputMC->FindObject(fillthis))->Fill(tofSignal);
	  fillthis="hpiTPCSignal";
	  ((TH1F*)fOutputMC->FindObject(fillthis))->Fill(dedxTPC);
	  fillthis="hpiptProng";
	  ((TH1F*)fOutputMC->FindObject(fillthis))->Fill(part->PtProng(iprong));
          fillthis="hpiptProngLcPt";
	  ((TH1F*)fOutputMC->FindObject(fillthis))->Fill(part->Pt(),part->PtProng(iprong));
	  fillthis="hpid0Prong";
	  ((TH1F*)fOutputMC->FindObject(fillthis))->Fill(part->Getd0Prong(iprong));
          fillthis="hpid0ProngLcPt";
	  ((TH1F*)fOutputMC->FindObject(fillthis))->Fill(part->Pt(),part->Getd0Prong(iprong));
          
	  fillthis="hpiSignalVspTPC";
	  ((TH2F*)fOutputMC->FindObject(fillthis))->Fill(momTPC,dedxTPC);
	  fillthis="hpiSignalVspTOF";
	  ((TH2F*)fOutputMC->FindObject(fillthis))->Fill(momTOF,tofSignal);
	  AliPID::EParticleType typepi;
	  typepi=AliPID::EParticleType(2);
	  if(hasTPC) {
	    //Double_t nsigmap = fUtilPid->NumberOfSigmasTPC(prong,typepi);
	    Double_t nsigmap = fPIDResponse->NumberOfSigmasTPC(prong,typepi);
	    fillthis="hpiSigmaVspTPC";
	    ((TH2F*)fOutputMC->FindObject(fillthis))->Fill(momTPC,nsigmap);
	  }
	  if(hasTOF){
	    //Double_t nsigma = fUtilPid->NumberOfSigmasTOF(prong,typepi);
	    Double_t nsigma=fPIDResponse->NumberOfSigmasTOF(prong,typepi); // protection against 'old' AODs...
	    fillthis="hpiSigmaVspTOF";
	    ((TH2F*)fOutputMC->FindObject(fillthis))->Fill(momTOF,nsigma);
	  }

	  if (fReadMC) { // fill hpiIDTot ---> PID contamination denominator
	                 //      hIDGood, hnoID ---> PID numerators
	    fillthis="hpiIDTot";
	    ((TH1F*)fOutputMC->FindObject(fillthis))->Fill(part->PtProng(iprong));
	    isID++;
	    if(pdgProngMC[iprong]==211) { // id pions
	      fillthis="hpiIDGood";
	      ((TH1F*)fOutputMC->FindObject(fillthis))->Fill(part->PtProng(iprong));
	      isCorrect++;
	    }
	    else { // misidentified electrons, muons, kaons and protons
	      fillthis="hnopiIDpi";
	      ((TH1F*)fOutputMC->FindObject(fillthis))->Fill(part->PtProng(iprong));
	    }
	  }
	  break;

	default:
	  break;
	}

      }//end loop on prongs

      //Jaime Lc checks
      fillthis="hLcIDTot";
      if(isID==3) ((TH1F*)fOutputMC->FindObject(fillthis))->Fill(part->Pt());
      fillthis="hLcIDGood";
      if(isCorrect==3) ((TH1F*)fOutputMC->FindObject(fillthis))->Fill(part->Pt());
      fillthis="hLcnoID";
      if(isCorrect<3) ((TH1F*)fOutputMC->FindObject(fillthis))->Fill(part->Pt());

      //now histograms where we don't need to check identity
      fillthis="hLcpt";
      ((TH1F*)fOutputMC->FindObject(fillthis))->Fill(part->Pt());
      fillthis = "hDist12toPrim";
      ((TH1F*)fOutputMC->FindObject(fillthis))->Fill(part->GetDist12toPrim());
      ((TH1F*)fOutputMC->FindObject(fillthis))->Fill(part->GetDist23toPrim());
      fillthis = "hDist12toPrimLcPt";
      ((TH1F*)fOutputMC->FindObject(fillthis))->Fill(part->Pt(),part->GetDist12toPrim());
      ((TH1F*)fOutputMC->FindObject(fillthis))->Fill(part->Pt(),part->GetDist23toPrim());
      fillthis = "hSigmaVert";
      ((TH1F*)fOutputMC->FindObject(fillthis))->Fill(part->GetSigmaVert());
      fillthis = "hSigmaVertLcPt";
      ((TH1F*)fOutputMC->FindObject(fillthis))->Fill(part->Pt(),part->GetSigmaVert());
      fillthis = "hDCAs";
      Double_t dcas[3];
      part->GetDCAs(dcas);
      for (Int_t idca=0;idca<3;idca++) ((TH1F*)fOutputMC->FindObject(fillthis))->Fill(dcas[idca]);
      fillthis = "hDCAsLcPt";
      for (Int_t idca=0;idca<3;idca++) ((TH1F*)fOutputMC->FindObject(fillthis))->Fill(part->Pt(),dcas[idca]);
      fillthis = "hCosPointingAngle";
      ((TH1F*)fOutputMC->FindObject(fillthis))->Fill(part->CosPointingAngle());
      fillthis = "hCosPointingAngleLcPt";
      ((TH1F*)fOutputMC->FindObject(fillthis))->Fill(part->Pt(),part->CosPointingAngle());
      fillthis = "hDecayLength";
      ((TH1F*)fOutputMC->FindObject(fillthis))->Fill(part->DecayLength());
      fillthis = "hDecayLengthLcPt";
      ((TH1F*)fOutputMC->FindObject(fillthis))->Fill(part->Pt(),part->DecayLength());
      Double_t sum2=part->Getd0Prong(0)*part->Getd0Prong(0)+
	part->Getd0Prong(1)*part->Getd0Prong(1)+
	part->Getd0Prong(2)*part->Getd0Prong(2);
      fillthis = "hSum2";
      ((TH1F*)fOutputMC->FindObject(fillthis))->Fill(sum2);
      fillthis = "hSum2LcPt";
      ((TH1F*)fOutputMC->FindObject(fillthis))->Fill(part->Pt(),sum2);
      fillthis = "hptmax";
      ((TH1F*)fOutputMC->FindObject(fillthis))->Fill(ptmax);
      fillthis = "hptmaxLcPt";
      ((TH1F*)fOutputMC->FindObject(fillthis))->Fill(part->Pt(),ptmax);


    } // end if (selection)

  }else if( lab<0 && fReadMC && !IsInjected ) { // bkg x MC

    fillthis="hbMass";

    if (part->Charge()>0) {      //Lc background
      if ( (pdgProngPID[0]==2212) && (pdgProngPID[1]==321) && (pdgProngPID[2]==211) )
	((TH1F*)fOutputMC->FindObject(fillthis))->Fill(minvLcpKpi); // bkg X MC
      else if ( (pdgProngPID[0]==211) && (pdgProngPID[1]==321) && (pdgProngPID[2]==2212) )
	((TH1F*)fOutputMC->FindObject(fillthis))->Fill(minvLcpiKp); // bkg X MC
    }
    else if (part->Charge()<0){ //anti-Lc background
      if ( (pdgProngPID[0]==2212) && (pdgProngPID[1]==321) && (pdgProngPID[2]==211) )
	((TH1F*)fOutputMC->FindObject(fillthis))->Fill(minvLcpKpi); // bkg X MC
      else if ( (pdgProngPID[0]==211) && (pdgProngPID[1]==321) && (pdgProngPID[2]==2212) )
	((TH1F*)fOutputMC->FindObject(fillthis))->Fill(minvLcpiKp); // bkg X MC
    }


    //apply cut on invariant mass on the pair
    if (TMath::Abs(minvLcpKpi-mPDG)<invmasscut || TMath::Abs(minvLcpiKp-mPDG)<invmasscut) {

      if (selection>0) { // 3prongs candidate x Lc (yes PID)

	Double_t ptmax=0.;
	Int_t isID=0;
	Int_t isCorrect=0;
	for (Int_t iprong=0; iprong<3; iprong++) {
	  if(part->PtProng(iprong)>ptmax)ptmax=part->PtProng(iprong);

	  AliAODTrack *prong = (AliAODTrack*)part->GetDaughter(iprong);
	  AliAODPid *pidObjtrk = (AliAODPid*)prong->GetDetPid();
	  AliAODPidHF *pidObj = (AliAODPidHF*)cuts->GetPidHF();
	  Bool_t hasTOF=pidObj->CheckStatus(prong,"TOF");
	  Bool_t hasTPC=pidObj->CheckStatus(prong,"TPC");
	  Double_t tofSignal=0.;
	  Double_t dedxTPC=0.;
	  Double_t momTOF=0.;
	  Double_t momTPC=0.;
	  if(hasTOF) {
	    momTOF = prong->P();
	    tofSignal=pidObjtrk->GetTOFsignal();
	  }
	  if(hasTPC) {
	    momTPC = pidObjtrk->GetTPCmomentum();
	    dedxTPC=pidObjtrk->GetTPCsignal();
	  }

	  switch (pdgProngPID[iprong]) {
	  case 2212:
	    fillthis="hbpTOFSignal";
	    ((TH1F*)fOutputMC->FindObject(fillthis))->Fill(tofSignal);
	    fillthis="hbpTPCSignal";
	    ((TH1F*)fOutputMC->FindObject(fillthis))->Fill(dedxTPC);
	    fillthis="hbpptProng";
	    ((TH1F*)fOutputMC->FindObject(fillthis))->Fill(part->PtProng(iprong));
            fillthis="hbpptProngLcPt";
	    ((TH1F*)fOutputMC->FindObject(fillthis))->Fill(part->Pt(),part->PtProng(iprong));
	    fillthis="hbpd0Prong";
	    ((TH1F*)fOutputMC->FindObject(fillthis))->Fill(part->Getd0Prong(iprong));
            fillthis="hbpd0ProngLcPt";
	    ((TH1F*)fOutputMC->FindObject(fillthis))->Fill(part->Pt(),part->Getd0Prong(iprong));
	    fillthis="hbpSignalVspTPC";
	    ((TH2F*)fOutputMC->FindObject(fillthis))->Fill(momTPC,dedxTPC);
	    fillthis="hbpSignalVspTOF";
	    ((TH2F*)fOutputMC->FindObject(fillthis))->Fill(momTOF,tofSignal);
	    AliPID::EParticleType typep;
	    typep=AliPID::EParticleType(4);
	    if(hasTPC) {
	      //Double_t nsigmap = fUtilPid->NumberOfSigmasTPC(prong,typep);
	      Double_t nsigmap = fPIDResponse->NumberOfSigmasTPC(prong,typep);
	      fillthis="hbpSigmaVspTPC";
	      ((TH2F*)fOutputMC->FindObject(fillthis))->Fill(momTPC,nsigmap);
	    }
	    if(hasTOF){
	      //Double_t nsigma = fUtilPid->NumberOfSigmasTOF(prong,typep);
	      Double_t nsigma=fPIDResponse->NumberOfSigmasTOF(prong,typep); // protection against 'old' AODs...
	      fillthis="hbpSigmaVspTOF";
	      ((TH2F*)fOutputMC->FindObject(fillthis))->Fill(momTOF,nsigma);
	    }
	    if(fReadMC){ // fill hbpIDTot ---> PID contamination denominator
	                 //      hbIDGood, hbnoID ---> PID numerators
	      fillthis="hbpIDTot";
	      ((TH1F*)fOutputMC->FindObject(fillthis))->Fill(part->PtProng(iprong));
	      isID++;
	      if(pdgProngMC[iprong]==2212) { // id protons
		fillthis="hbpIDGood";
		((TH1F*)fOutputMC->FindObject(fillthis))->Fill(part->PtProng(iprong));
		isCorrect++;
	      }
	      else { // misidentified electrons, muons, pions and kaons
		fillthis="hbnopIDp";
		((TH1F*)fOutputMC->FindObject(fillthis))->Fill(part->PtProng(iprong));
	      }
	    }
	    break;

	  case 321:
	    fillthis="hbKTOFSignal";
	    ((TH1F*)fOutputMC->FindObject(fillthis))->Fill(tofSignal);
	    fillthis="hbKTPCSignal";
	    ((TH1F*)fOutputMC->FindObject(fillthis))->Fill(dedxTPC);
	    fillthis="hbKptProng";
	    ((TH1F*)fOutputMC->FindObject(fillthis))->Fill(part->PtProng(iprong));
            fillthis="hbKptProngLcPt";
	    ((TH1F*)fOutputMC->FindObject(fillthis))->Fill(part->Pt(),part->PtProng(iprong));
	    fillthis="hbKd0Prong";
	    ((TH1F*)fOutputMC->FindObject(fillthis))->Fill(part->Getd0Prong(iprong));
            fillthis="hbKd0ProngLcPt";
	    ((TH1F*)fOutputMC->FindObject(fillthis))->Fill(part->Pt(),part->Getd0Prong(iprong));
	    fillthis="hbKSignalVspTPC";
	    ((TH2F*)fOutputMC->FindObject(fillthis))->Fill(momTPC,dedxTPC);
	    fillthis="hbKSignalVspTOF";
	    ((TH2F*)fOutputMC->FindObject(fillthis))->Fill(momTOF,tofSignal);
	    AliPID::EParticleType typek;
	    typek=AliPID::EParticleType(3);
	    if(hasTPC) {
	      //Double_t nsigmap = fUtilPid->NumberOfSigmasTPC(prong,typek);
	      Double_t nsigmap = fPIDResponse->NumberOfSigmasTPC(prong,typek);
	      fillthis="hbKSigmaVspTPC";
	      ((TH2F*)fOutputMC->FindObject(fillthis))->Fill(momTPC,nsigmap);
	    }
	    if(hasTOF){
	      //Double_t nsigma = fUtilPid->NumberOfSigmasTOF(prong,typek);
	      Double_t nsigma=fPIDResponse->NumberOfSigmasTOF(prong,typek); // protection against 'old' AODs...
	      fillthis="hbKSigmaVspTOF";
	      ((TH2F*)fOutputMC->FindObject(fillthis))->Fill(momTOF,nsigma);
	    }
	    if (fReadMC) { // fill hbKIDTot ---> PID contamination denominator
	                   //      hIDGood, hnoID ---> PID numerators
	      fillthis="hbKIDTot";
	      ((TH1F*)fOutputMC->FindObject(fillthis))->Fill(part->PtProng(iprong));
	      isID++;
	      if(pdgProngMC[iprong]==321) { // id kaons
		fillthis="hbKIDGood";
		((TH1F*)fOutputMC->FindObject(fillthis))->Fill(part->PtProng(iprong));
		isCorrect++;
	      }
	      else { // misidentified electrons, muons, pions and protons
		fillthis="hbnokIDk";
		((TH1F*)fOutputMC->FindObject(fillthis))->Fill(part->PtProng(iprong));
	      }
	    }
	    break;

	  case 211:
	    fillthis="hbpiTOFSignal";
	    ((TH1F*)fOutputMC->FindObject(fillthis))->Fill(tofSignal);
	    fillthis="hbpiTPCSignal";
	    ((TH1F*)fOutputMC->FindObject(fillthis))->Fill(dedxTPC);
	    fillthis="hbpiptProng";
	    ((TH1F*)fOutputMC->FindObject(fillthis))->Fill(part->PtProng(iprong));
	    fillthis="hbpiptProngLcPt";
	    ((TH1F*)fOutputMC->FindObject(fillthis))->Fill(part->PtProng(iprong));
	    fillthis="hbpid0Prong";
	    ((TH1F*)fOutputMC->FindObject(fillthis))->Fill(part->Getd0Prong(iprong));
            fillthis="hbpid0ProngLcPt";
	    ((TH1F*)fOutputMC->FindObject(fillthis))->Fill(part->Pt(),part->Getd0Prong(iprong));
	    fillthis="hbpiSignalVspTPC";
	    ((TH2F*)fOutputMC->FindObject(fillthis))->Fill(momTPC,dedxTPC);
	    fillthis="hbpiSignalVspTOF";
	    ((TH2F*)fOutputMC->FindObject(fillthis))->Fill(momTOF,tofSignal);
	    AliPID::EParticleType typepi;
	    typepi=AliPID::EParticleType(2);
	    if(hasTPC) {
	      //Double_t nsigmap = fUtilPid->NumberOfSigmasTPC(prong,typepi);
	      Double_t nsigmap = fPIDResponse->NumberOfSigmasTPC(prong,typepi);
	      fillthis="hbpiSigmaVspTPC";
	      ((TH2F*)fOutputMC->FindObject(fillthis))->Fill(momTPC,nsigmap);
	    }
	    if(hasTOF){
	      //Double_t nsigma = fUtilPid->NumberOfSigmasTOF(prong,typepi);
	      Double_t nsigma=fPIDResponse->NumberOfSigmasTOF(prong,typepi); // protection against 'old' AODs...
	      fillthis="hbpiSigmaVspTOF";
	      ((TH2F*)fOutputMC->FindObject(fillthis))->Fill(momTOF,nsigma);
	    }
	    if (fReadMC) { // fill hbpiIDTot ---> PID contamination denominator
	                   //      hIDGood, hnonIDTot ---> PID numerators
	      fillthis="hbpiIDTot";
	      ((TH1F*)fOutputMC->FindObject(fillthis))->Fill(part->PtProng(iprong));
	      isID++;
	      if(pdgProngMC[iprong]==211) { // id pions
		fillthis="hbpiIDGood";
		((TH1F*)fOutputMC->FindObject(fillthis))->Fill(part->PtProng(iprong));
		isCorrect++;
	      }
	      else { // misidentified electrons, muons, kaons and protons
		fillthis="hbnopiIDpi";
		((TH1F*)fOutputMC->FindObject(fillthis))->Fill(part->PtProng(iprong));
	      }
	    }
	    break;

	  default:
	    break;
	  }

	} //end loop on prongs
	
	fillthis="hbLcIDTot";
	if(isID==3) ((TH1F*)fOutputMC->FindObject(fillthis))->Fill(part->Pt());	
	fillthis="hbLcIDGood";
	if(isCorrect==3) ((TH1F*)fOutputMC->FindObject(fillthis))->Fill(part->Pt());
	fillthis="hbLcnoID";
	if(isCorrect<3) ((TH1F*)fOutputMC->FindObject(fillthis))->Fill(part->Pt());

	//now histograms where we don't need to check identity
	fillthis="hbLcpt";
	((TH1F*)fOutputMC->FindObject(fillthis))->Fill(part->Pt());
	
	fillthis = "hbDist12toPrim";
	((TH1F*)fOutputMC->FindObject(fillthis))->Fill(part->GetDist12toPrim());
	((TH1F*)fOutputMC->FindObject(fillthis))->Fill(part->GetDist23toPrim());
	fillthis = "hbDist12toPrimLcPt";
	((TH1F*)fOutputMC->FindObject(fillthis))->Fill(part->Pt(),part->GetDist12toPrim());
	((TH1F*)fOutputMC->FindObject(fillthis))->Fill(part->Pt(),part->GetDist23toPrim());
	fillthis = "hbSigmaVert";
	((TH1F*)fOutputMC->FindObject(fillthis))->Fill(part->GetSigmaVert());
	fillthis = "hbSigmaVertLcPt";
	((TH1F*)fOutputMC->FindObject(fillthis))->Fill(part->Pt(),part->GetSigmaVert());
	fillthis = "hbDCAs";
	Double_t dcas[3];
	part->GetDCAs(dcas);
	for (Int_t idca=0;idca<3;idca++) ((TH1F*)fOutputMC->FindObject(fillthis))->Fill(dcas[idca]);
	fillthis = "hbDCAsLcPt";
	for (Int_t idca=0;idca<3;idca++) ((TH1F*)fOutputMC->FindObject(fillthis))->Fill(part->Pt(),dcas[idca]);
	fillthis = "hbCosPointingAngle";
	((TH1F*)fOutputMC->FindObject(fillthis))->Fill(part->CosPointingAngle());
	fillthis = "hbCosPointingAngleLcPt";
	((TH1F*)fOutputMC->FindObject(fillthis))->Fill(part->Pt(),part->CosPointingAngle());
	fillthis = "hbDecayLength";
	((TH1F*)fOutputMC->FindObject(fillthis))->Fill(part->DecayLength());
	fillthis = "hbDecayLengthLcPt";
	((TH1F*)fOutputMC->FindObject(fillthis))->Fill(part->Pt(),part->DecayLength());
	Double_t sum2=part->Getd0Prong(0)*part->Getd0Prong(0)+
	part->Getd0Prong(1)*part->Getd0Prong(1)+
	part->Getd0Prong(2)*part->Getd0Prong(2);
	fillthis = "hbSum2";
	((TH1F*)fOutputMC->FindObject(fillthis))->Fill(sum2);
	fillthis = "hbSum2LcPt";
	((TH1F*)fOutputMC->FindObject(fillthis))->Fill(part->Pt(),sum2);
	fillthis = "hbptmax";
	((TH1F*)fOutputMC->FindObject(fillthis))->Fill(ptmax);
	fillthis = "hbptmaxLcPt";
	((TH1F*)fOutputMC->FindObject(fillthis))->Fill(part->Pt(),ptmax);

      } // end if (selection)

    } // end mass cut

  } // end background case

  return;

}

//-------------------------
void AliAnalysisTaskSELambdac::FillAPrioriConcentrations(AliAODRecoDecayHF3Prong *part,
							    AliRDHFCutsLctopKpi *cuts,
							    AliAODEvent* aod,
							    TClonesArray *arrMC)
{

  // FillAPrioriConcentrations
  //cuts->SetUsePID(kFALSE); //Annalisa
  Int_t isSelected3ProngByLc=cuts->IsSelected(part,AliRDHFCuts::kCandidate,aod);

  TString fillthis="";

  if(isSelected3ProngByLc>0 && fReadMC) {

    for (Int_t ii=0; ii<3; ii++) {
      AliAODTrack *prongTest=(AliAODTrack*)part->GetDaughter(ii);
      if (!prongTest) continue;
      Int_t labprongTest = prongTest->GetLabel();
      if(labprongTest<0) continue;
      AliAODMCParticle *mcprongTest = (AliAODMCParticle*)arrMC->At(labprongTest);

      switch (TMath::Abs(mcprongTest->GetPdgCode())) {
      case 11:
	fillthis="hElIn3Prong";
	((TH1F*)fAPriori->FindObject(fillthis))->Fill(part->PtProng(ii));
	break;
      case 13:
	fillthis="hMuIn3Prong";
	((TH1F*)fAPriori->FindObject(fillthis))->Fill(part->PtProng(ii));
	break;
      case 211:
	fillthis="hPiIn3Prong";
	((TH1F*)fAPriori->FindObject(fillthis))->Fill(part->PtProng(ii));
	break;
      case 321:
	fillthis="hKaIn3Prong";
	((TH1F*)fAPriori->FindObject(fillthis))->Fill(part->PtProng(ii));
	break;
      case 2212:
	fillthis="hPrIn3Prong";
	((TH1F*)fAPriori->FindObject(fillthis))->Fill(part->PtProng(ii));
	break;
      default:
	break;
      }
    }
  }

  //cuts->SetUsePID(kTRUE); //Annalisa

}

//-----------------------
Bool_t AliAnalysisTaskSELambdac::Is3ProngFromPDG(AliAODRecoDecayHF3Prong *d,
						    TClonesArray *arrayMC,
						    Int_t pdgToBeCompared)
{
  //
  // Returs kTRUE if at least one of tracks in the current 3prong
  // cames from particle with pdgCode=pdgToBeCompared
  //

  Bool_t localFlag = kFALSE;

  for (Int_t ii=0;ii<3;ii++) {
    AliAODTrack *daugh = (AliAODTrack*)d->GetDaughter(ii);

    Int_t lab = daugh->GetLabel();
    if(lab<0) continue;
    AliAODMCParticle *part= (AliAODMCParticle*)arrayMC->At(lab);
    //Int_t partPdg = part->GetPdgCode();
    //printf(" -------------------------------------- partLab=%d partPdg=%d ---\n", lab, partPdg);

    Int_t motherLabel = part->GetMother();
    if(motherLabel<0) continue;

    AliAODMCParticle *motherPart = 0;
    Int_t motherPdg = 0;
    while (!localFlag && motherLabel>=0) {
      motherPart = (AliAODMCParticle*)arrayMC->At(motherLabel);
      motherPdg = motherPart->GetPdgCode();
      //printf(" ------------- motherLab=%d motherPdg=%d ---\n", motherLabel, motherPdg);
      if (TMath::Abs(motherPdg)==pdgToBeCompared) {
	//printf("-------------------------- trovato quark c/cbar\n");
	localFlag = kTRUE;
      }
      else {
	motherLabel = motherPart->GetMother();
      }	
    }

  }

  return localFlag;

}

//-----------------------
Bool_t AliAnalysisTaskSELambdac::IsTrackFromPDG(const AliAODTrack *daugh,
						   TClonesArray *arrayMC,
						   Int_t pdgToBeCompared)
{
  //
  // Returs kTRUE if at the tracks comes
  // from particle with pdgCode=pdgToBeCompared
  //

  Bool_t localFlag = kFALSE;

  Int_t lab = daugh->GetLabel();
  if(lab<0) return localFlag;
  AliAODMCParticle *part= (AliAODMCParticle*)arrayMC->At(lab);
  //Int_t partPdg = part->GetPdgCode();
  //printf(" -------------------------------------- partLab=%d partPdg=%d ---\n", lab, partPdg);

  Int_t motherLabel = part->GetMother();
  if(motherLabel<0) return localFlag;

  AliAODMCParticle *motherPart = 0;
  Int_t motherPdg = 0;
  while (!localFlag && motherLabel>=0) {
    motherPart = (AliAODMCParticle*)arrayMC->At(motherLabel);
    motherPdg = motherPart->GetPdgCode();
    //printf(" ------------- motherLab=%d motherPdg=%d ---\n", motherLabel, motherPdg);
    if (TMath::Abs(motherPdg)==pdgToBeCompared) {
      //printf("-------------------------- trovato quark c/cbar\n");
      localFlag = kTRUE;
    }
    else {
      motherLabel = motherPart->GetMother();
    }	
  }

  return localFlag;

}




//-----------------------
Bool_t AliAnalysisTaskSELambdac::IsThereAGeneratedLc(TClonesArray *arrayMC)
{
  //
  // Returs kTRUE if there is a lepton related to the LambdaC
  //

  Bool_t localFlag = kFALSE;
 
  AliAODMCParticle *searchLc;
  
  for (Int_t iii=0; iii<arrayMC->GetEntries(); iii++) {
    searchLc = (AliAODMCParticle*)arrayMC->At(iii);
    Int_t searchLcpdg =  searchLc->GetPdgCode();
    if (TMath::Abs(searchLcpdg) == 4122){
      localFlag = kTRUE;
      break;  
    }
  }
  

  return localFlag;

}
//_________________________________
Int_t AliAnalysisTaskSELambdac::NumberPrimaries(const AliAODEvent *aods)
{
///////////////estimate primaries
  Int_t counter=0;

  
  TClonesArray *aodtracks=(TClonesArray *)aods->GetTracks();
 
  // for(Int_t ji=0;ji<aods->GetNumberOfTracks();ji++)
  for(Int_t ji=0;ji<aodtracks->GetEntriesFast();ji++)
    {
      AliAODTrack*aodTrack=(AliAODTrack*)aodtracks->UncheckedAt(ji);
      if(aodTrack->IsPrimaryCandidate()) counter++;
    }
  return counter;
  
}  
//_________________________________

void AliAnalysisTaskSELambdac::MultiplicityStudies(AliAODRecoDecayHF3Prong *part,
						      AliRDHFCutsLctopKpi *cuts,
						      AliAODEvent* aod,
						      TClonesArray *arrMC,
						      Bool_t &flag1,Bool_t &flag2,Bool_t &flag3,
						      Bool_t &flag4, Bool_t &flag5, Bool_t &flag6)
{


// Multiplicity studies

  TString fillthis="";

  Int_t pdgDgLctopKpi[3]={2212,321,211};
  Int_t lab=-9999;
  if(fReadMC)
    lab=part->MatchToMC(4122,arrMC,3,pdgDgLctopKpi); //return MC particle label if the array corresponds to a Lc, -1 if not (cf. AliAODRecoDecay.cxx)

  //cuts->SetUsePID(kFALSE); //Annalisa
  Int_t isSelected3ProngByLc=cuts->IsSelected(part,AliRDHFCuts::kCandidate,aod);

  if(isSelected3ProngByLc>0 && fReadMC) {
    flag1 = kTRUE;
    if (lab>=0) {
      flag3 = kTRUE;
    }

    Bool_t is3ProngFromJPsi = Is3ProngFromPDG(part,arrMC,443);
    if (is3ProngFromJPsi) flag6=is3ProngFromJPsi;

    Bool_t is3ProngFromC = Is3ProngFromPDG(part,arrMC,4);
    if (is3ProngFromC) flag4=is3ProngFromC;

    Bool_t is3ProngFromB = Is3ProngFromPDG(part,arrMC,5);
    if (is3ProngFromB) flag5=is3ProngFromB;

    for (Int_t ii=0; ii<3; ii++) {
      AliAODTrack *prongTest=(AliAODTrack*)part->GetDaughter(ii);
      if (!prongTest) continue;
      Int_t labprongTest = prongTest->GetLabel();
      if(labprongTest<0) continue;
      AliAODMCParticle *mcprongTest = (AliAODMCParticle*)arrMC->At(labprongTest);

      switch (TMath::Abs(mcprongTest->GetPdgCode())) {
      case 11:
	fillthis="hElIn3Prong";
	((TH1F*)fAPriori->FindObject(fillthis))->Fill(part->PtProng(ii));
	if (IsTrackFromPDG(prongTest,arrMC,443)) {
	  fillthis="hElIn3Prong6";
	  ((TH1F*)fAPriori->FindObject(fillthis))->Fill(part->PtProng(ii));
	}
	break;
      case 13:
	fillthis="hMuIn3Prong";
	((TH1F*)fAPriori->FindObject(fillthis))->Fill(part->PtProng(ii));
	if (IsTrackFromPDG(prongTest,arrMC,443)) {
	  fillthis="hMuIn3Prong6";
	  ((TH1F*)fAPriori->FindObject(fillthis))->Fill(part->PtProng(ii));
	}
	break;
      case 211:
	fillthis="hPiIn3Prong";
	((TH1F*)fAPriori->FindObject(fillthis))->Fill(part->PtProng(ii));
	if (IsTrackFromPDG(prongTest,arrMC,443)) {
	  fillthis="hPiIn3Prong6";
	  ((TH1F*)fAPriori->FindObject(fillthis))->Fill(part->PtProng(ii));
	}
	break;
      case 321:
	fillthis="hKaIn3Prong";
	((TH1F*)fAPriori->FindObject(fillthis))->Fill(part->PtProng(ii));
	if (IsTrackFromPDG(prongTest,arrMC,443)) {
	  fillthis="hKaIn3Prong6";
	  ((TH1F*)fAPriori->FindObject(fillthis))->Fill(part->PtProng(ii));
	}
	break;
      case 2212:
	fillthis="hPrIn3Prong";
	((TH1F*)fAPriori->FindObject(fillthis))->Fill(part->PtProng(ii));
	if (IsTrackFromPDG(prongTest,arrMC,443)) {
	  fillthis="hPrIn3Prong6";
	  ((TH1F*)fAPriori->FindObject(fillthis))->Fill(part->PtProng(ii));
	}
	break;
      default:
	break;
      }

      if (lab>=0) {
	switch (TMath::Abs(mcprongTest->GetPdgCode())) {
	case 11:
	  fillthis="hElIn3Prong1";
	  ((TH1F*)fAPriori->FindObject(fillthis))->Fill(part->PtProng(ii));
	  break;
	case 13:
	  fillthis="hMuIn3Prong1";
	  ((TH1F*)fAPriori->FindObject(fillthis))->Fill(part->PtProng(ii));
	  break;
	case 211:
	  fillthis="hPiIn3Prong1";
	  ((TH1F*)fAPriori->FindObject(fillthis))->Fill(part->PtProng(ii));
	  break;
	case 321:
	  fillthis="hKaIn3Prong1";
	  ((TH1F*)fAPriori->FindObject(fillthis))->Fill(part->PtProng(ii));
	  break;
	case 2212:
	  fillthis="hPrIn3Prong1";
	  ((TH1F*)fAPriori->FindObject(fillthis))->Fill(part->PtProng(ii));
	  break;
	default:
	  break;
	}
      } else {
	switch (TMath::Abs(mcprongTest->GetPdgCode())) {
	case 11:
	  fillthis="hElIn3Prong2";
	  ((TH1F*)fAPriori->FindObject(fillthis))->Fill(part->PtProng(ii));
	  if (IsTrackFromPDG(prongTest,arrMC,4)) {
	    fillthis="hElIn3Prong3";
	    ((TH1F*)fAPriori->FindObject(fillthis))->Fill(part->PtProng(ii));
	  }
	  if (IsTrackFromPDG(prongTest,arrMC,5)) {
	    fillthis="hElIn3Prong4";
	    ((TH1F*)fAPriori->FindObject(fillthis))->Fill(part->PtProng(ii));
	  }
	  if (!IsTrackFromPDG(prongTest,arrMC,4) && !IsTrackFromPDG(prongTest,arrMC,5)) {
	    fillthis="hElIn3Prong5";
	    ((TH1F*)fAPriori->FindObject(fillthis))->Fill(part->PtProng(ii));
	  }
	  break;
	case 13:
	  fillthis="hMuIn3Prong2";
	  ((TH1F*)fAPriori->FindObject(fillthis))->Fill(part->PtProng(ii));
	  if (IsTrackFromPDG(prongTest,arrMC,4)) {
	    fillthis="hMuIn3Prong3";
	    ((TH1F*)fAPriori->FindObject(fillthis))->Fill(part->PtProng(ii));
	  }
	  if (IsTrackFromPDG(prongTest,arrMC,5)) {
	    fillthis="hMuIn3Prong4";
	    ((TH1F*)fAPriori->FindObject(fillthis))->Fill(part->PtProng(ii));
	  }
	  if (!IsTrackFromPDG(prongTest,arrMC,4) && !IsTrackFromPDG(prongTest,arrMC,5)) {
	    fillthis="hMuIn3Prong5";
	    ((TH1F*)fAPriori->FindObject(fillthis))->Fill(part->PtProng(ii));
	  }
	  break;
	case 211:
	  fillthis="hPiIn3Prong2";
	  ((TH1F*)fAPriori->FindObject(fillthis))->Fill(part->PtProng(ii));
	  if (IsTrackFromPDG(prongTest,arrMC,4)) {
	    fillthis="hPiIn3Prong3";
	    ((TH1F*)fAPriori->FindObject(fillthis))->Fill(part->PtProng(ii));
	  }
	  if (IsTrackFromPDG(prongTest,arrMC,5)) {
	    fillthis="hPiIn3Prong4";
	    ((TH1F*)fAPriori->FindObject(fillthis))->Fill(part->PtProng(ii));
	  }
	  if (!IsTrackFromPDG(prongTest,arrMC,4) && !IsTrackFromPDG(prongTest,arrMC,5)) {
	    fillthis="hPiIn3Prong5";
	    ((TH1F*)fAPriori->FindObject(fillthis))->Fill(part->PtProng(ii));
	  }
	  break;
	case 321:
	  fillthis="hKaIn3Prong2";
	  ((TH1F*)fAPriori->FindObject(fillthis))->Fill(part->PtProng(ii));
	  if (IsTrackFromPDG(prongTest,arrMC,4)) {
	    fillthis="hKaIn3Prong3";
	    ((TH1F*)fAPriori->FindObject(fillthis))->Fill(part->PtProng(ii));
	  }
	  if (IsTrackFromPDG(prongTest,arrMC,5)) {
	    fillthis="hKaIn3Prong4";
	    ((TH1F*)fAPriori->FindObject(fillthis))->Fill(part->PtProng(ii));
	  }
	  if (!IsTrackFromPDG(prongTest,arrMC,4) && !IsTrackFromPDG(prongTest,arrMC,5)) {
	    fillthis="hKaIn3Prong5";
	    ((TH1F*)fAPriori->FindObject(fillthis))->Fill(part->PtProng(ii));
	  }
	  break;
	case 2212:
	  fillthis="hPrIn3Prong2";
	  ((TH1F*)fAPriori->FindObject(fillthis))->Fill(part->PtProng(ii));
	  if (IsTrackFromPDG(prongTest,arrMC,4)) {
	    fillthis="hPrIn3Prong3";
	    ((TH1F*)fAPriori->FindObject(fillthis))->Fill(part->PtProng(ii));
	  }
	  if (IsTrackFromPDG(prongTest,arrMC,5)) {
	    fillthis="hPrIn3Prong4";
	    ((TH1F*)fAPriori->FindObject(fillthis))->Fill(part->PtProng(ii));
	  }
	  if (!IsTrackFromPDG(prongTest,arrMC,4) && !IsTrackFromPDG(prongTest,arrMC,5)) {
	    fillthis="hPrIn3Prong5";
	    ((TH1F*)fAPriori->FindObject(fillthis))->Fill(part->PtProng(ii));
	  }
	  break;
	default:
	  break;
	}
      }
      /*
	if (is3ProngFromC) {
	switch (TMath::Abs(mcprongTest->GetPdgCode())) {
	case 11:
	fillthis="hElIn3Prong3";
	((TH1F*)fAPriori->FindObject(fillthis))->Fill(part->PtProng(ii));
	break;
	case 13:
	fillthis="hMuIn3Prong3";
	((TH1F*)fAPriori->FindObject(fillthis))->Fill(part->PtProng(ii));
	break;
	case 211:
	fillthis="hPiIn3Prong3";
	((TH1F*)fAPriori->FindObject(fillthis))->Fill(part->PtProng(ii));
	break;
	case 321:
	fillthis="hKaIn3Prong3";
	((TH1F*)fAPriori->FindObject(fillthis))->Fill(part->PtProng(ii));
	break;
	case 2212:
	fillthis="hPrIn3Prong3";
	((TH1F*)fAPriori->FindObject(fillthis))->Fill(part->PtProng(ii));
	break;
	default:
	break;
	}
	} else { // !is3ProngFromC
	if (is3ProngFromB) {
	switch (TMath::Abs(mcprongTest->GetPdgCode())) {
	case 11:
	fillthis="hElIn3Prong4";
	((TH1F*)fAPriori->FindObject(fillthis))->Fill(part->PtProng(ii));
	break;
	case 13:
	fillthis="hMuIn3Prong4";
	((TH1F*)fAPriori->FindObject(fillthis))->Fill(part->PtProng(ii));
	break;
	case 211:
	fillthis="hPiIn3Prong4";
	((TH1F*)fAPriori->FindObject(fillthis))->Fill(part->PtProng(ii));
	break;
	case 321:
	fillthis="hKaIn3Prong4";
	((TH1F*)fAPriori->FindObject(fillthis))->Fill(part->PtProng(ii));
	break;
	case 2212:
	fillthis="hPrIn3Prong4";
	((TH1F*)fAPriori->FindObject(fillthis))->Fill(part->PtProng(ii));
	break;
	default:
	break;
	}
	} else {//!is3ProngFromB
	switch (TMath::Abs(mcprongTest->GetPdgCode())) {
	case 11:
	fillthis="hElIn3Prong5";
	((TH1F*)fAPriori->FindObject(fillthis))->Fill(part->PtProng(ii));
	break;
	case 13:
	fillthis="hMuIn3Prong5";
	((TH1F*)fAPriori->FindObject(fillthis))->Fill(part->PtProng(ii));
	break;
	case 211:
	fillthis="hPiIn3Prong5";
	((TH1F*)fAPriori->FindObject(fillthis))->Fill(part->PtProng(ii));
	break;
	case 321:
	fillthis="hKaIn3Prong5";
	((TH1F*)fAPriori->FindObject(fillthis))->Fill(part->PtProng(ii));
	break;
	case 2212:
	fillthis="hPrIn3Prong5";
	((TH1F*)fAPriori->FindObject(fillthis))->Fill(part->PtProng(ii));
	break;
	default:
	break;
	}
	}
	}
      */
    }
  }

  Double_t mPDG=TDatabasePDG::Instance()->GetParticle(4122)->Mass();
  Double_t invmasscut=0.05;

  Double_t minvLcpKpi = part->InvMassLcpKpi();
  Double_t minvLcpiKp = part->InvMassLcpiKp();

  //cuts->SetUsePID(kTRUE); //Annalisa
  Int_t isSelected3ProngByLcPID=cuts->IsSelected(part,AliRDHFCuts::kCandidate,aod);
  if (isSelected3ProngByLcPID>0) {
    if (TMath::Abs(minvLcpKpi-mPDG)<invmasscut || TMath::Abs(minvLcpiKp-mPDG)<invmasscut) {
      flag2 = kTRUE;
    }
  }


}
