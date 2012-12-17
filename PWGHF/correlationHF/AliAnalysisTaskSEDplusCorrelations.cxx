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

/* $Id$ */

//*************************************************************************
// Class AliAnalysisTaskSEDplusCorrelations
// AliAnalysisTaskSE for the D+ candidates Invariant Mass Histogram and 
//comparison of heavy-flavour decay candidates
// to MC truth (kinematics stored in the AOD)
// Authors: Sadhana Dash (correlation)
/////////////////////////////////////////////////////////////

#include <TClonesArray.h>
#include <TNtuple.h>
#include <TCanvas.h>
#include <TList.h>
#include <TString.h>
#include <TDatabasePDG.h>

#include <AliAnalysisDataSlot.h>
#include <AliAnalysisDataContainer.h>
#include "AliAnalysisManager.h"
#include "AliRDHFCutsDplustoKpipi.h"
#include "AliAODHandler.h"
#include "AliAODEvent.h"
#include "AliAODPidHF.h"
#include "AliAODVertex.h"
#include "AliAODTrack.h"
#include "AliAODRecoDecayHF3Prong.h"
#include "AliAnalysisVertexingHF.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisTaskSEDplusCorrelations.h"
#include "AliNormalizationCounter.h"
#include "AliVParticle.h"
#include "AliHFAssociatedTrackCuts.h"
#include "AliReducedParticle.h"
#include "AliHFCorrelator.h"


using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskSEDplusCorrelations)


//________________________________________________________________________
AliAnalysisTaskSEDplusCorrelations::AliAnalysisTaskSEDplusCorrelations():
AliAnalysisTaskSE(),
  fOutput(0),
  fCorrelator(0),
  fSelect(0),
  fDisplacement(0),
  fHistNEvents(0),
  fMCSources(0),
  fK0Origin(0),
  fKaonOrigin(0),
  fInvMassK0S(0),
  fEventMix(0),
  fPtVsMass(0),
  fPtVsMassTC(0),
  fYVsPt(0),
  fYVsPtTC(0),
  fYVsPtSig(0),
  fYVsPtSigTC(0),
  fUpmasslimit(1.965),
  fLowmasslimit(1.765),
  fNPtBins(0),
  fBinWidth(0.002),
  fListCuts(0),
  fListCutsAsso(0), 
  fRDCutsAnalysis(0),
  fCuts(0),
  fCounter(0),
  fReadMC(kFALSE),
  fUseBit(kTRUE),
  fMixing(kFALSE),
  fSystem(kFALSE)
{
  // Default constructor
   
  for(Int_t i=0;i<3*kMaxPtBins;i++){
    fMassVsdPhiHistHad[i]=0;
    fMassVsdEtaHistHad[i]=0;
    fMassVsdPhiHistKaon[i]=0;
    fMassVsdEtaHistKaon[i]=0;
    fMassVsdPhiHistKshort[i]=0;
    fMassVsdEtaHistKshort[i]=0;
    fMassVsdPhiHistLeadHad[i]=0;
    fMassVsdEtaHistLeadHad[i]=0;
    fMassHistK0S[i]=0;
    fLeadPt[i]=0;
    fMassHist[i]=0;
    fMassHistTC[i]=0;
    fMassHistTCPlus[i]=0;
    fMassHistTCMinus[i]=0;
  }

  for(Int_t i=0;i<kMaxPtBins+1;i++){
    fArrayBinLimits[i]=0;
  }

}

//________________________________________________________________________
AliAnalysisTaskSEDplusCorrelations::AliAnalysisTaskSEDplusCorrelations(const char *name,AliRDHFCutsDplustoKpipi *dpluscutsana, AliHFAssociatedTrackCuts *asstrkcuts):
  AliAnalysisTaskSE(name),
  fOutput(0),
  fCorrelator(0),
  fSelect(0),
  fDisplacement(0),
  fHistNEvents(0),
  fMCSources(0),
  fK0Origin(0),
  fKaonOrigin(0),
  fInvMassK0S(0),
  fEventMix(0),
  fPtVsMass(0),
  fPtVsMassTC(0),
  fYVsPt(0),
  fYVsPtTC(0),
  fYVsPtSig(0),
  fYVsPtSigTC(0),
  fUpmasslimit(1.965),
  fLowmasslimit(1.765),
  fNPtBins(0),
  fBinWidth(0.002),
  fListCuts(0),
  fListCutsAsso(0),
  fRDCutsAnalysis(dpluscutsana),
  fCuts(asstrkcuts),
  fCounter(0),
  fReadMC(kFALSE),
  fUseBit(kTRUE),
  fMixing(kFALSE),
  fSystem(kFALSE)
{
  // 
  // Standrd constructor
  //
  fNPtBins=fRDCutsAnalysis->GetNPtBins();
    
  for(Int_t i=0;i<3*kMaxPtBins;i++){
    fMassVsdPhiHistHad[i]=0;
    fMassVsdEtaHistHad[i]=0;
    fMassVsdPhiHistKaon[i]=0;
    fMassVsdEtaHistKaon[i]=0;
    fMassVsdPhiHistKshort[i]=0;
    fMassVsdEtaHistKshort[i]=0;
    fMassVsdPhiHistLeadHad[i]=0;
    fMassVsdEtaHistLeadHad[i]=0;
    fMassHistK0S[i]=0;
    fLeadPt[i]=0;
    fMassHist[i]=0;
    fMassHistTC[i]=0;
    fMassHistTCPlus[i]=0;
    fMassHistTCMinus[i]=0;
  }
   
  
  for(Int_t i=0;i<kMaxPtBins+1;i++){
    fArrayBinLimits[i]=0;
  }
  
  
  // Default constructor
  // Output slot #1 writes into a TList container
  DefineOutput(1,TList::Class());  //My private output
  // Output slot #2 writes cut to private output
  //  DefineOutput(2,AliRDHFCutsDplusCorrelationstoKpipi::Class());
  DefineOutput(2,TList::Class());
  // Output slot #3 writes cut to private output
  DefineOutput(3,AliNormalizationCounter::Class());
  DefineOutput(4,AliHFAssociatedTrackCuts::Class());
  
 
}

//________________________________________________________________________
AliAnalysisTaskSEDplusCorrelations::~AliAnalysisTaskSEDplusCorrelations()
{
  //
  // Destructor
  //
  delete fOutput;
  delete fListCuts;
  delete fRDCutsAnalysis;
  delete fCuts;
  delete fCounter;
  delete fCorrelator; 


}  
//_________________________________________________________________
void  AliAnalysisTaskSEDplusCorrelations::SetMassLimits(Float_t range){
  // set invariant mass limits
  Float_t bw=GetBinWidth();
  fUpmasslimit = 1.865+range;
  fLowmasslimit = 1.865-range;
  SetBinWidth(bw);
}
//_________________________________________________________________
void  AliAnalysisTaskSEDplusCorrelations::SetMassLimits(Float_t lowlimit, Float_t uplimit){
  // set invariant mass limits
  if(uplimit>lowlimit)
    {
      Float_t bw=GetBinWidth();
      fUpmasslimit = lowlimit;
      fLowmasslimit = uplimit;
      SetBinWidth(bw);
    }
}
//________________________________________________________________
void AliAnalysisTaskSEDplusCorrelations::SetBinWidth(Float_t w){
  // Define width of mass bins
  Float_t width=w;
  Int_t nbins=(Int_t)((fUpmasslimit-fLowmasslimit)/width+0.5);
  Int_t missingbins=4-nbins%4;
  nbins=nbins+missingbins;
  width=(fUpmasslimit-fLowmasslimit)/nbins;
  if(missingbins!=0){
    printf("AliAnalysisTaskSEDplusCorrelations::SetBinWidth: W-bin width of %f will produce histograms not rebinnable by 4. New width set to %f\n",w,width);
  }
  else{
    if(fDebug>1) printf("AliAnalysisTaskSEDplusCorrelations::SetBinWidth: width set to %f\n",width);
  }
  fBinWidth=width;
}
//_________________________________________________________________
Int_t AliAnalysisTaskSEDplusCorrelations::GetNBinsHistos(){
  // Compute number of mass bins
  return (Int_t)((fUpmasslimit-fLowmasslimit)/fBinWidth+0.5);
}


//__________________________________________
void AliAnalysisTaskSEDplusCorrelations::Init(){
  //
  // Initialization
  //
  if(fDebug > 1) printf("AnalysisTaskSEDplusCorrelations::Init() \n");
  
  //PostData(2,fRDCutsloose);//we should then put those cuts in a tlist if we have more than 1
  fListCuts=new TList();
  // fListCutsAsso=new TList();
  
  AliRDHFCutsDplustoKpipi *analysis = new AliRDHFCutsDplustoKpipi(*fRDCutsAnalysis);
  analysis->SetName("AnalysisCuts");

  // AliHFAssociatedTrackCuts *trkcuts = new AliHFAssociatedTrackCuts(*fCuts);
  //trkcuts->SetName("Assotrkcuts");
  
  fListCuts->Add(analysis);
  //fListCuts->Add(trkcuts);

   

  PostData(2,fListCuts);
  PostData(4,fCuts);
  
  return;
}

//________________________________________________________________________
void AliAnalysisTaskSEDplusCorrelations::UserCreateOutputObjects()
{
  // Create the output container
  //
  if(fDebug > 1) printf("AnalysisTaskSEDplusCorrelations::UserCreateOutputObjects() \n");
  // correlator creation and definition

  Double_t Pi = TMath::Pi();
  fCorrelator = new AliHFCorrelator("Correlator",fCuts,fSystem); // fCuts is the hadron cut object, fSystem to switch between pp or PbPb
  fCorrelator->SetDeltaPhiInterval((-0.5-1./32)*Pi,(1.5-1./32)*Pi); // set correct phi interval
  //fCorrelator->SetDeltaPhiInterval(-1.57,4.71);
  fCorrelator->SetEventMixing(fMixing); //set kFALSE/kTRUE for mixing Off/On
  fCorrelator->SetAssociatedParticleType(fSelect); // set 1/2/3 for hadron/kaons/kzeros
  fCorrelator->SetApplyDisplacementCut(fDisplacement); //set kFALSE/kTRUE for using the displacement cut
  fCorrelator->SetUseMC(fReadMC);
  fCorrelator->SetPIDmode(2);

  Bool_t pooldef = fCorrelator->DefineEventPool();
  
  if(!pooldef) AliInfo("Warning:: Event pool not defined properly");


  // Several histograms are more conveniently managed in a TList
  fOutput = new TList();
  fOutput->SetOwner();
  fOutput->SetName("OutputHistos");

  TString hisname;
  Int_t index=0;
  Int_t nbins=GetNBinsHistos();


   Int_t nbinsphi = 32;
   Double_t philow = -0.5*Pi - Pi/32; // shift the bin by half the width so that at 0 is it the bin center
   Double_t phiup = 1.5*Pi - Pi/32;	

   Int_t nbinseta = 16;
   Double_t etalow = -1.6; // shift the bin by half the width so that at 0 is it the bin center
   Double_t etaup = +1.6;	
     

  
  for(Int_t i=0;i<fNPtBins;i++){

    index=GetHistoIndex(i);


    hisname.Form("hMassVsdPhiHad%d",i);
    fMassVsdPhiHistHad[index]=new TH2F(hisname.Data(),hisname.Data(),nbins,fLowmasslimit,fUpmasslimit,nbinsphi,philow,phiup);
    fMassVsdPhiHistHad[index]->Sumw2();

    hisname.Form("hMassVsdEtaHad%d",i);
    fMassVsdEtaHistHad[index]=new TH3D(hisname.Data(),hisname.Data(),nbins,fLowmasslimit,fUpmasslimit,nbinsphi,philow,phiup,nbinseta,etalow,etaup);
    fMassVsdEtaHistHad[index]->Sumw2();

    hisname.Form("hMassVsdPhiKaon%d",i);
    fMassVsdPhiHistKaon[index]=new TH2F(hisname.Data(),hisname.Data(),nbins,fLowmasslimit,fUpmasslimit,nbinsphi,philow,phiup);
    fMassVsdPhiHistKaon[index]->Sumw2();

    hisname.Form("hMassVsdEtaKaon%d",i);
    fMassVsdEtaHistKaon[index]=new TH3D(hisname.Data(),hisname.Data(),nbins,fLowmasslimit,fUpmasslimit,nbinsphi,philow,phiup,nbinseta,etalow,etaup);
    fMassVsdEtaHistKaon[index]->Sumw2();

    hisname.Form("hMassVsdPhiK0%d",i);
    fMassVsdPhiHistKshort[index]=new TH2F(hisname.Data(),hisname.Data(),nbins,fLowmasslimit,fUpmasslimit,nbinsphi,philow,phiup);
    fMassVsdPhiHistKshort[index]->Sumw2();

    hisname.Form("hMassVsdEtaK0%d",i);
    fMassVsdEtaHistKshort[index]=new TH3D(hisname.Data(),hisname.Data(),nbins,fLowmasslimit,fUpmasslimit,nbinsphi,philow,phiup,nbinseta,etalow,etaup);
    fMassVsdEtaHistKshort[index]->Sumw2();
    
    hisname.Form("hMassVsdPhiLeadHad%d",i);
    fMassVsdPhiHistLeadHad[index]=new TH2F(hisname.Data(),hisname.Data(),nbins,fLowmasslimit,fUpmasslimit,nbinsphi,philow,phiup);
    fMassVsdPhiHistLeadHad[index]->Sumw2();

    hisname.Form("hMassVsdEtaLeadHad%d",i);
    fMassVsdEtaHistLeadHad[index]=new TH3D(hisname.Data(),hisname.Data(),nbins,fLowmasslimit,fUpmasslimit,nbinsphi,philow,phiup,nbinseta,etalow,etaup);
    fMassVsdEtaHistLeadHad[index]->Sumw2();

    hisname.Form("hMassPtK0S%d",i);
    fMassHistK0S[index]=new TH1F(hisname.Data(),hisname.Data(),500,0.3,0.8);
    fMassHistK0S[index]->Sumw2();

    hisname.Form("hLeadPt%d",i);
    fLeadPt[index]=new TH1F(hisname.Data(),hisname.Data(),500,0.0,50);
    fLeadPt[index]->Sumw2();

       
    hisname.Form("hMassPt%d",i);
    fMassHist[index]=new TH1F(hisname.Data(),hisname.Data(),nbins,fLowmasslimit,fUpmasslimit);
    fMassHist[index]->Sumw2();
    
   
    hisname.Form("hMassPt%dTC",i);
    fMassHistTC[index]=new TH1F(hisname.Data(),hisname.Data(),nbins,fLowmasslimit,fUpmasslimit);
    fMassHistTC[index]->Sumw2();

    hisname.Form("hMassPt%dTCPlus",i);
    fMassHistTCPlus[index]=new TH1F(hisname.Data(),hisname.Data(),nbins,fLowmasslimit,fUpmasslimit);
    fMassHistTCPlus[index]->Sumw2();

    hisname.Form("hMassPt%dTCMinus",i);
    fMassHistTCMinus[index]=new TH1F(hisname.Data(),hisname.Data(),nbins,fLowmasslimit,fUpmasslimit);
    fMassHistTCMinus[index]->Sumw2();


    
    index=GetSignalHistoIndex(i); 

    hisname.Form("hMassVsdPhiHadSig%d",i);
    fMassVsdPhiHistHad[index]=new TH2F(hisname.Data(),hisname.Data(),nbins,fLowmasslimit,fUpmasslimit,nbinsphi,philow,phiup);
    fMassVsdPhiHistHad[index]->Sumw2();

    hisname.Form("hMassVsdEtaHadSig%d",i);
    fMassVsdEtaHistHad[index]=new TH3D(hisname.Data(),hisname.Data(),nbins,fLowmasslimit,fUpmasslimit,nbinsphi,philow,phiup,nbinseta,etalow,etaup);
    fMassVsdEtaHistHad[index]->Sumw2();

    hisname.Form("hMassVsdPhiKaonSig%d",i);
    fMassVsdPhiHistKaon[index]=new TH2F(hisname.Data(),hisname.Data(),nbins,fLowmasslimit,fUpmasslimit,nbinsphi,philow,phiup);
    fMassVsdPhiHistKaon[index]->Sumw2();

    hisname.Form("hMassVsdEtaKaonSig%d",i);
    fMassVsdEtaHistKaon[index]=new TH3D(hisname.Data(),hisname.Data(),nbins,fLowmasslimit,fUpmasslimit,nbinsphi,philow,phiup,nbinseta,etalow,etaup);
    fMassVsdEtaHistKaon[index]->Sumw2();

    hisname.Form("hMassVsdPhiK0Sig%d",i);
    fMassVsdPhiHistKshort[index]=new TH2F(hisname.Data(),hisname.Data(),nbins,fLowmasslimit,fUpmasslimit,nbinsphi,philow,phiup);
    fMassVsdPhiHistKshort[index]->Sumw2();

    hisname.Form("hMassVsdEtaK0Sig%d",i);
    fMassVsdEtaHistKshort[index]=new TH3D(hisname.Data(),hisname.Data(),nbins,fLowmasslimit,fUpmasslimit,nbinsphi,philow,phiup,nbinseta,etalow,etaup);
    fMassVsdEtaHistKshort[index]->Sumw2();
    
    hisname.Form("hMassVsdPhiLeadHadSig%d",i);
    fMassVsdPhiHistLeadHad[index]=new TH2F(hisname.Data(),hisname.Data(),nbins,fLowmasslimit,fUpmasslimit,nbinsphi,philow,phiup);
    fMassVsdPhiHistLeadHad[index]->Sumw2();

    hisname.Form("hMassVsdEtaLeadHadSig%d",i);
    fMassVsdEtaHistLeadHad[index]=new TH3D(hisname.Data(),hisname.Data(),nbins,fLowmasslimit,fUpmasslimit,nbinsphi,philow,phiup,nbinseta,etalow,etaup);
    fMassVsdEtaHistLeadHad[index]->Sumw2();

 
    hisname.Form("hSigPtK0S%d",i);
    fMassHistK0S[index]=new TH1F(hisname.Data(),hisname.Data(),500,0.3,0.8);
    fMassHistK0S[index]->Sumw2();

    hisname.Form("hSigLeadPt%d",i);
    fLeadPt[index]=new TH1F(hisname.Data(),hisname.Data(),500,0.0,50);
    fLeadPt[index]->Sumw2();
    
    hisname.Form("hSigPt%d",i);
    fMassHist[index]=new TH1F(hisname.Data(),hisname.Data(),nbins,fLowmasslimit,fUpmasslimit);
    fMassHist[index]->Sumw2();

    hisname.Form("hSigPt%dTC",i);
    fMassHistTC[index]=new TH1F(hisname.Data(),hisname.Data(),nbins,fLowmasslimit,fUpmasslimit);
    fMassHistTC[index]->Sumw2();
    hisname.Form("hSigPt%dTCPlus",i);
    fMassHistTCPlus[index]=new TH1F(hisname.Data(),hisname.Data(),nbins,fLowmasslimit,fUpmasslimit);
    fMassHistTCPlus[index]->Sumw2();
    hisname.Form("hSigPt%dTCMinus",i);
    fMassHistTCMinus[index]=new TH1F(hisname.Data(),hisname.Data(),nbins,fLowmasslimit,fUpmasslimit);
    fMassHistTCMinus[index]->Sumw2();


   
    index=GetBackgroundHistoIndex(i);

    hisname.Form("hMassVsdPhiBkgHad%d",i);
    fMassVsdPhiHistHad[index]=new TH2F(hisname.Data(),hisname.Data(),nbins,fLowmasslimit,fUpmasslimit,nbinsphi,philow,phiup);
    fMassVsdPhiHistHad[index]->Sumw2();

    hisname.Form("hMassVsdEtaBkgHad%d",i);
    fMassVsdEtaHistHad[index]=new TH3D(hisname.Data(),hisname.Data(),nbins,fLowmasslimit,fUpmasslimit,nbinsphi,philow,phiup,nbinseta,etalow,etaup);
    fMassVsdEtaHistHad[index]->Sumw2();

    hisname.Form("hMassVsdPhiBkgKaon%d",i);
    fMassVsdPhiHistKaon[index]=new TH2F(hisname.Data(),hisname.Data(),nbins,fLowmasslimit,fUpmasslimit,nbinsphi,philow,phiup);
    fMassVsdPhiHistKaon[index]->Sumw2();

    hisname.Form("hMassVsdEtaBkgKaon%d",i);
    fMassVsdEtaHistKaon[index]=new TH3D(hisname.Data(),hisname.Data(),nbins,fLowmasslimit,fUpmasslimit,nbinsphi,philow,phiup,nbinseta,etalow,etaup);
    fMassVsdEtaHistKaon[index]->Sumw2();

    hisname.Form("hMassVsdPhiBkgKshort%d",i);
    fMassVsdPhiHistKshort[index]=new TH2F(hisname.Data(),hisname.Data(),nbins,fLowmasslimit,fUpmasslimit,nbinsphi,philow,phiup);
    fMassVsdPhiHistKshort[index]->Sumw2();


    hisname.Form("hMassVsdPhiBkgKshort%d",i);
    fMassVsdEtaHistKshort[index]=new TH3D(hisname.Data(),hisname.Data(),nbins,fLowmasslimit,fUpmasslimit,nbinsphi,philow,phiup,nbinseta,etalow,etaup);
    fMassVsdEtaHistKshort[index]->Sumw2();

    hisname.Form("hMassVsdPhiBkgLeadHad%d",i);
    fMassVsdPhiHistLeadHad[index]=new TH2F(hisname.Data(),hisname.Data(),nbins,fLowmasslimit,fUpmasslimit,nbinsphi,philow,phiup);
    fMassVsdPhiHistLeadHad[index]->Sumw2();

    hisname.Form("hMassVsdPhiBkgKshort%d",i);
    fMassVsdEtaHistLeadHad[index]=new TH3D(hisname.Data(),hisname.Data(),nbins,fLowmasslimit,fUpmasslimit,nbinsphi,philow,phiup,nbinseta,etalow,etaup);
    fMassVsdEtaHistLeadHad[index]->Sumw2();
    
    hisname.Form("hBkgPtK0S%d",i);
    fMassHistK0S[index]=new TH1F(hisname.Data(),hisname.Data(),500,0.3,0.8);
    fMassHistK0S[index]->Sumw2();

    hisname.Form("hLeadBkgPt%d",i);
    fLeadPt[index]=new TH1F(hisname.Data(),hisname.Data(),500,0.0,50);
    fLeadPt[index]->Sumw2();

    hisname.Form("hBkgPt%dTC",i);
    fMassHistTC[index]=new TH1F(hisname.Data(),hisname.Data(),nbins,fLowmasslimit,fUpmasslimit);
    fMassHistTC[index]->Sumw2();

    hisname.Form("hBkgPt%dTCPlus",i);
    fMassHistTCPlus[index]=new TH1F(hisname.Data(),hisname.Data(),nbins,fLowmasslimit,fUpmasslimit);
    fMassHistTCPlus[index]->Sumw2();

    hisname.Form("hBkgPt%dTCMinus",i);
    fMassHistTCMinus[index]=new TH1F(hisname.Data(),hisname.Data(),nbins,fLowmasslimit,fUpmasslimit);
    fMassHistTCMinus[index]->Sumw2();
  }
    

  for(Int_t i=0; i<3*fNPtBins; i++){
    fOutput->Add(fMassVsdPhiHistHad[i]);
    fOutput->Add(fMassVsdEtaHistHad[i]);
    fOutput->Add(fMassVsdPhiHistKaon[i]);
    fOutput->Add(fMassVsdEtaHistKaon[i]);
    fOutput->Add(fMassVsdPhiHistKshort[i]);
    fOutput->Add(fMassVsdEtaHistKshort[i]);
    fOutput->Add(fMassVsdPhiHistLeadHad[i]);
    fOutput->Add(fMassVsdEtaHistLeadHad[i]);
    fOutput->Add(fMassHistK0S[i]);
    fOutput->Add(fLeadPt[i]);
    fOutput->Add(fMassHist[i]);
    fOutput->Add(fMassHistTC[i]);
    fOutput->Add(fMassHistTCPlus[i]);
    fOutput->Add(fMassHistTCMinus[i]);
    


  }
  fInvMassK0S = new TH2F("K0S","K0S", 500,0.3,0.8,500,0,50);
  fInvMassK0S->GetXaxis()->SetTitle("Invariant Mass (#pi #pi) (GeV/c^{2})");
  fInvMassK0S->GetYaxis()->SetTitle("K0S pt (GeV/c)");
  fOutput->Add(fInvMassK0S);

  fMCSources = new TH1F("MCSources","Origin of associated particles in MC", 10, -0.5, 9.5);
  fMCSources->GetXaxis()->SetBinLabel(1,"All ");
  fMCSources->GetXaxis()->SetBinLabel(2," from Heavy flavour");
  fMCSources->GetXaxis()->SetBinLabel(3," from c->D");
  fMCSources->GetXaxis()->SetBinLabel(4," from b->D");
  fMCSources->GetXaxis()->SetBinLabel(5," from b->B");
  if(fReadMC) fOutput->Add(fMCSources);

  fK0Origin = new TH1F("K0Origin","Origin of K0 ", 10, -0.5, 5.5);
  fK0Origin->GetXaxis()->SetBinLabel(1,"All K0s");
  fK0Origin->GetXaxis()->SetBinLabel(2,"K0s from Heavy flavour");
  fK0Origin->GetXaxis()->SetBinLabel(3,"K0s from Charm");
  fK0Origin->GetXaxis()->SetBinLabel(4,"K0s from Beauty");
  if(fReadMC) fOutput->Add(fK0Origin);

  fKaonOrigin = new TH1F("K0Origin","Origin of Kaon ", 10, -0.5, 5.5);
  fKaonOrigin->GetXaxis()->SetBinLabel(1,"All Kaons");
  fKaonOrigin->GetXaxis()->SetBinLabel(2,"Kaons from Heavy flavour");
  fKaonOrigin->GetXaxis()->SetBinLabel(3,"Kaons from Charm");
  fKaonOrigin->GetXaxis()->SetBinLabel(4,"Kaons from Beauty");
  if(fReadMC) fOutput->Add(fKaonOrigin);
  
	
  
  fEventMix = new TH2F("EventMixingCheck","EventMixingCheck",5,-0.5,4.5,7,-0.5,6.5);
  if(fMixing)fOutput->Add(fEventMix);
	
  
  fHistNEvents = new TH1F("fHistNEvents", "number of events ",11,-0.5,10.5);
  fHistNEvents->GetXaxis()->SetBinLabel(1,"nEventsAnal");
  fHistNEvents->GetXaxis()->SetBinLabel(2,"nEvents accepted");
  fHistNEvents->GetXaxis()->SetBinLabel(3,"Rejected due to trigger");
  fHistNEvents->GetXaxis()->SetBinLabel(4,"Rejected pileup events");
  fHistNEvents->GetXaxis()->SetBinLabel(5,"Rejected due to centrality"); 
  fHistNEvents->GetXaxis()->SetBinLabel(6,"Rejected due to vtxz");
  fHistNEvents->GetXaxis()->SetBinLabel(7,"Rejected due to Physics Sel");
  fHistNEvents->GetXaxis()->SetBinLabel(8,"Total no. of candidate");
  fHistNEvents->GetXaxis()->SetBinLabel(9,"no. of cand wo bitmask");
  fHistNEvents->GetXaxis()->SetBinLabel(10,"D+ after loose cuts");
  fHistNEvents->GetXaxis()->SetBinLabel(11,"D+ after tight cuts");
 
  fHistNEvents->GetXaxis()->SetNdivisions(1,kFALSE);  
  fHistNEvents->Sumw2();
  fHistNEvents->SetMinimum(0);
  fOutput->Add(fHistNEvents);

  fPtVsMass=new TH2F("hPtVsMass","PtVsMass (prod. cuts)",nbins,fLowmasslimit,fUpmasslimit,200,0.,20.);
  fPtVsMassTC=new TH2F("hPtVsMassTC","PtVsMass (analysis cuts)",nbins,fLowmasslimit,fUpmasslimit,200,0.,20.);  
  fYVsPt=new TH2F("hYVsPt","YvsPt (prod. cuts)",40,0.,20.,80,-2.,2.);
  fYVsPtTC=new TH2F("hYVsPtTC","YvsPt (analysis cuts)",40,0.,20.,80,-2.,2.);
  fYVsPtSig=new TH2F("hYVsPtSig","YvsPt (MC, only sig., prod. cuts)",40,0.,20.,80,-2.,2.);
  fYVsPtSigTC=new TH2F("hYVsPtSigTC","YvsPt (MC, only Sig, analysis cuts)",40,0.,20.,80,-2.,2.);

  fOutput->Add(fPtVsMass);
  fOutput->Add(fPtVsMassTC);
  fOutput->Add(fYVsPt);
  fOutput->Add(fYVsPtTC);
  fOutput->Add(fYVsPtSig);
  fOutput->Add(fYVsPtSigTC);


  // Counter for Normalization
  TString normName="NormalizationCounter";
  AliAnalysisDataContainer *cont = GetOutputSlot(3)->GetContainer();
  if(cont)normName=(TString)cont->GetName();
  fCounter = new AliNormalizationCounter(normName.Data());
  fCounter->Init();

  

  PostData(1,fOutput);      
  PostData(3,fCounter);      
  return;
}

//________________________________________________________________________
void AliAnalysisTaskSEDplusCorrelations::UserExec(Option_t */*option*/)
{
  // Do the analysis
  AliAODEvent *aod = dynamic_cast<AliAODEvent*> (InputEvent());
  TClonesArray *array3Prong = 0;

  if(!fMixing){
    if(fSelect==1) cout << "TASK::Correlation with hadrons on SE "<< endl;
    if(fSelect==2) cout << "TASK::Correlation with kaons on SE "<< endl;
    if(fSelect==3) cout << "TASK::Correlation with kzeros on SE "<< endl;
  }
  if(fMixing){
    if(fSelect==1) cout << "TASK::Correlation with hadrons on ME "<< endl;
    if(fSelect==2) cout << "TASK::Correlation with kaons on ME "<< endl;
    if(fSelect==3) cout << "TASK::Correlation with kzeros on ME "<< endl;
  }
  
  
  if(!aod && AODEvent() && IsStandardAOD()) {

    aod = dynamic_cast<AliAODEvent*> (AODEvent());
    AliAODHandler* aodHandler = (AliAODHandler*) 
      ((AliAnalysisManager::GetAnalysisManager())->GetOutputEventHandler());
    if(aodHandler->GetExtensions()) {
      AliAODExtension *ext = (AliAODExtension*)aodHandler->GetExtensions()->FindObject("AliAOD.VertexingHF.root");
      AliAODEvent *aodFromExt = ext->GetAOD();
      array3Prong=(TClonesArray*)aodFromExt->GetList()->FindObject("Charm3Prong");
    }
  } else if(aod) {
    array3Prong=(TClonesArray*)aod->GetList()->FindObject("Charm3Prong");
  }

  if(!aod || !array3Prong) {
    printf("AliAnalysisTaskSEDplusCorrelations::UserExec: Charm3Prong branch not found!\n");
    return;
  }

  
  // the AODs with null vertex pointer didn't pass the PhysSel
  if(!aod->GetPrimaryVertex()||TMath::Abs(aod->GetMagneticField())<0.001) return;
  fCounter->StoreEvent(aod,fRDCutsAnalysis,fReadMC);
  fHistNEvents->Fill(0); // count event
  

  Bool_t isEvSel=fRDCutsAnalysis->IsEventSelected(aod);  
   

  
  PostData(1,fOutput);
  if(!isEvSel)return;
  fHistNEvents->Fill(1);

  // set PIDResponse for associated tracks
  fCorrelator->SetPidAssociated(); 

  TClonesArray *arrayMC=0;
  AliAODMCHeader *mcHeader=0;
  // AOD primary vertex
  AliAODVertex *vtx1 = (AliAODVertex*)aod->GetPrimaryVertex();
  //    vtx1->Print();
  TString primTitle = vtx1->GetTitle();
  //if(primTitle.Contains("VertexerTracks") && vtx1->GetNContributors()>0)fHistNEvents->Fill(2);
   
  // load MC particles
  if(fReadMC){
    
    arrayMC =  (TClonesArray*)aod->GetList()->FindObject(AliAODMCParticle::StdBranchName());
    if(!arrayMC) {
      printf("AliAnalysisTaskSEDplusCorrelations::UserExec: MC particles branch not found!\n");
      return;
    }
    
    // load MC header
    mcHeader =  (AliAODMCHeader*)aod->GetList()->FindObject(AliAODMCHeader::StdBranchName());
    if(!mcHeader) {
      printf("AliAnalysisTaskSEDplusCorrelations::UserExec: MC header branch not found!\n");
      return;
    }
  }

  //HFCorrelators initialization (for this event)

  fCorrelator->SetAODEvent(aod); // set the event to be processedfCorrelator->
  Bool_t correlatorON = fCorrelator->Initialize(); //define the pool for mixing
  if(!correlatorON) {
    AliInfo("AliHFCorrelator didn't initialize the pool correctly or processed a bad event");
    return;
  }
  if(fReadMC) fCorrelator->SetMCArray(arrayMC);


 
  Int_t n3Prong = array3Prong->GetEntriesFast();
   printf("Number of D+->Kpipi: %d and of tracks: %d\n",n3Prong,aod->GetNumberOfTracks());  
  Int_t nOS=0;
  Int_t index;
  Int_t pdgDgDplustoKpipi[3]={321,211,211};
  Int_t nSelectedloose=0,nSelectedtight=0;

     
   
  for (Int_t i3Prong = 0; i3Prong < n3Prong; i3Prong++) {
    AliAODRecoDecayHF3Prong *d = (AliAODRecoDecayHF3Prong*)array3Prong->UncheckedAt(i3Prong);
    fHistNEvents->Fill(7);
    
    if(fUseBit && !d->HasSelectionBit(AliRDHFCuts::kDplusCuts)){
      fHistNEvents->Fill(8);
      continue;
    }

    Int_t passTightCuts=fRDCutsAnalysis->IsSelected(d,AliRDHFCuts::kAll,aod);
     
     if(!fRDCutsAnalysis->GetIsSelectedCuts()) continue;
      
    Bool_t unsetvtx=kFALSE;
    if(!d->GetOwnPrimaryVtx()){
      d->SetOwnPrimaryVtx(vtx1);
      unsetvtx=kTRUE;
    }

      
    Double_t ptCand = d->Pt();
    Int_t iPtBin = fRDCutsAnalysis->PtBin(ptCand);
    Bool_t recVtx=kFALSE;
    AliAODVertex *origownvtx=0x0; 
    if(fRDCutsAnalysis->GetIsPrimaryWithoutDaughters()){
      if(d->GetOwnPrimaryVtx()) origownvtx=new AliAODVertex(*d->GetOwnPrimaryVtx());	
      if(fRDCutsAnalysis->RecalcOwnPrimaryVtx(d,aod))recVtx=kTRUE;
      else fRDCutsAnalysis->CleanOwnPrimaryVtx(d,aod,origownvtx);
    }
      
     
      
                
    Int_t labDp=-1;
    Bool_t isDplus = kFALSE;
      
    if(fReadMC){
      labDp = d->MatchToMC(411,arrayMC,3,pdgDgDplustoKpipi);
      if(labDp>=0){
	isDplus = kTRUE;
      }
    }

          
    Double_t invMass=d->InvMassDplus();
    Double_t rapid=d->YDplus();
    fYVsPt->Fill(ptCand,rapid);
    if(passTightCuts>0.) {fYVsPtTC->Fill(ptCand,rapid);nOS++;}
    //printf("****************: %d and of tracks: %d\n",n3Prong,nOS);
    Bool_t isFidAcc=fRDCutsAnalysis->IsInFiducialAcceptance(ptCand,rapid);
    if(isFidAcc){
      fPtVsMass->Fill(invMass,ptCand);
      if(passTightCuts>0) fPtVsMassTC->Fill(invMass,ptCand);
    }
            
           
    if(iPtBin>=0){
      index=GetHistoIndex(iPtBin);
      if(isFidAcc){
	fHistNEvents->Fill(9);
	nSelectedloose++;
	fMassHist[index]->Fill(invMass);
	
	// loop for tight cuts
       	if(passTightCuts>0){   
	  fHistNEvents->Fill(10);
	  nSelectedtight++;
	  fMassHistTC[index]->Fill(invMass);

	  //Dplus info

	  Double_t phiDplus = fCorrelator->SetCorrectPhiRange(d->Phi());
	  fCorrelator->SetTriggerParticleProperties(d->Pt(),phiDplus,d->Eta());
	    
	  Int_t nIDs[3] = {-9999999};
	  nIDs[0] = ((AliAODTrack*)d->GetDaughter(0))->GetID();
	  nIDs[1] = ((AliAODTrack*)d->GetDaughter(1))->GetID();
	  nIDs[2] = ((AliAODTrack*)d->GetDaughter(2))->GetID();

	  Double_t ptlead = 0;
	  Double_t philead = 0;
	  Double_t etalead = 0;
	  Double_t refpt = 0;
	    


	  Bool_t execPool = fCorrelator->ProcessEventPool();

	  //     printf("*************: %d\n",execPool);
	  if(fMixing && !execPool) {
	    AliInfo("Mixed event analysis: pool is not ready");
	    continue;
	  }
	  Int_t NofEventsinPool = 1;
	  if(fMixing) {
	    NofEventsinPool = fCorrelator->GetNofEventsInPool();
	  }
		
		
		
	  for (Int_t jMix =0; jMix < NofEventsinPool; jMix++){// loop on events in the pool; if it is SE analysis, stops at one

	    Bool_t analyzetracks = fCorrelator->ProcessAssociatedTracks(jMix);
	    if(!analyzetracks) {
	      AliInfo("AliHFCorrelator::Cannot process the track array");
	      continue;
	    }
		 
	    //start the track loop
		  
	    // Int_t NofTracks = fCorrelator->GetNofTracks();

	    //cout<<"*******"<<NofTracks<<endl;

	    for (Int_t iTrack = 0;iTrack<fCorrelator->GetNofTracks();iTrack++) {	               
	      Bool_t runcorrelation = fCorrelator->Correlate(iTrack);
		       
	      if(!runcorrelation) continue;
	      Double_t DeltaPhi = fCorrelator->GetDeltaPhi();
	      Double_t DeltaEta = fCorrelator->GetDeltaEta();
		       
	      AliReducedParticle* redpart = fCorrelator->GetAssociatedParticle();
	      Double_t phiHad=redpart->Phi();
	      Double_t ptHad=redpart->Pt();
	      Double_t etaHad=redpart->Eta();
	      Int_t label = redpart->GetLabel();
	      Int_t trackid = redpart->GetID();
	      phiHad = fCorrelator->SetCorrectPhiRange(phiHad);

		       
	      //  discard the dplus daughters
	      if (!fMixing){
		if( trackid == nIDs[0] || trackid == nIDs[1] || trackid == nIDs[2]) continue;
	      }
	      // discard the negative id tracks 
	      if(trackid < 0) continue;

		    
	      FillCorrelations(d,DeltaPhi,DeltaEta,index,fSelect);
		    
	      // For leading particle
		    
	      if (ptHad > refpt) {
		refpt = ptHad; ptlead = ptHad;
		philead = d->Phi() - phiHad;
		etalead = d->Eta() - etaHad;
		if (philead < (-1)*TMath::Pi()/2) philead += 2*TMath::Pi();
		if (philead > 3*TMath::Pi()/2) philead -= 2*TMath::Pi();
		      
	      }

	      // montecarlo
	       
	      if(fReadMC && isDplus) {

		index=GetSignalHistoIndex(iPtBin);

		Bool_t* partSource = fCuts->IsMCpartFromHF(label,arrayMC); // check source of associated particle (hadron/kaon/K0)
		FillMCCorrelations(d,DeltaPhi,DeltaEta,index,partSource,fSelect);   
		delete partSource;
		
			 
	      } // readMC 

	    }//count good tracks

	    // For leading particle		  		  
	    fMassVsdPhiHistLeadHad[index]->Fill(invMass,philead);	
	    fMassVsdEtaHistLeadHad[index]->Fill(invMass,philead,etalead);

	    fLeadPt[index]->Fill(ptlead);
	     
	    if(fReadMC && isDplus) {
	      index=GetSignalHistoIndex(iPtBin);
	      fMassVsdPhiHistLeadHad[index]->Fill(invMass,philead);
	      fMassVsdEtaHistLeadHad[index]->Fill(invMass,philead,etalead);	
	      fLeadPt[index]->Fill(ptlead);
		     
	    }
	      
	  }//jmix
		
	  }// tc 
      }//fid
		
		
		
		
		
    }
    if(recVtx)fRDCutsAnalysis->CleanOwnPrimaryVtx(d,aod,origownvtx);
      
    if(unsetvtx) d->UnsetOwnPrimaryVtx();
    
  }

  if(fMixing){
    Bool_t updated = fCorrelator->PoolUpdate();
	
    if(!updated) AliInfo("Pool was not updated");
  }
  fCounter->StoreCandidates(aod,nSelectedloose,kTRUE);
  fCounter->StoreCandidates(aod,nSelectedtight,kFALSE);
  	
      
  PostData(1,fOutput); 
  PostData(2,fListCuts);
  PostData(3,fCounter);    
  return;
}

//________________________________________________________________________
void AliAnalysisTaskSEDplusCorrelations::Terminate(Option_t */*option*/)
{
  // Terminate analysis
  //
  if(fDebug > 1) printf("AnalysisTaskSEDplusCorrelations: Terminate() \n");

  fOutput = dynamic_cast<TList*> (GetOutputData(1));
  if (!fOutput) {     
    printf("ERROR: fOutput not available\n");
    return;
  }
  fHistNEvents = dynamic_cast<TH1F*>(fOutput->FindObject("fHistNEvents"));

  TString hisname;
  Int_t index=0;

  for(Int_t i=0;i<fNPtBins;i++){
    index=GetHistoIndex(i);
    
    hisname.Form("hMassPt%dTC",i);
    fMassHistTC[index]=dynamic_cast<TH1F*>(fOutput->FindObject(hisname.Data()));
  } 
    
  TCanvas *c1=new TCanvas("c1","D+ invariant mass distribution",500,500);
  c1->cd();
  TH1F *hMassPt=(TH1F*)fOutput->FindObject("hMassPt3TC");
  hMassPt->SetLineColor(kBlue);
  hMassPt->SetXTitle("M[GeV/c^{2}]"); 
  hMassPt->Draw();
 
  return;
}


//________________________________________________________________________
void AliAnalysisTaskSEDplusCorrelations::FillCorrelations(AliAODRecoDecayHF3Prong* d, Double_t deltaPhi, Double_t deltaEta, Int_t ind, Int_t sel) const{
  //Filling histogams
  
  Double_t invMass=d->InvMassDplus();
	

  if(sel==1){	  	
    fMassVsdPhiHistHad[ind]->Fill(invMass,deltaPhi);
    fMassVsdEtaHistHad[ind]->Fill(invMass,deltaPhi,deltaEta);
  }
  if(sel==2){
    fMassVsdPhiHistKaon[ind]->Fill(invMass,deltaPhi);
    fMassVsdEtaHistKaon[ind]->Fill(invMass,deltaPhi,deltaEta);
  }
  if(sel==3){
    fMassVsdPhiHistKshort[ind]->Fill(invMass,deltaPhi);
    fMassVsdEtaHistKshort[ind]->Fill(invMass,deltaPhi,deltaEta);
  }
	 
  return;
}



//________________________________________________________________________
void AliAnalysisTaskSEDplusCorrelations::FillMCCorrelations(AliAODRecoDecayHF3Prong* d, Double_t deltaPhi, Double_t deltaEta, Int_t ind,Bool_t* mcSource, Int_t sel) const{
  // Filling histos with MC information

  Double_t invMass=d->InvMassDplus();
  

  if(sel==1){
    fMassVsdPhiHistHad[ind]->Fill(invMass,deltaPhi);
    fMassVsdEtaHistHad[ind]->Fill(invMass,deltaPhi,deltaEta);

  	
    fMCSources->Fill(0);
	
    if(mcSource[2]&&mcSource[0]){ // is from charm ->D
      fMCSources->Fill(1);
      fMCSources->Fill(2);
    }
	if(mcSource[2]&&mcSource[1]){ // is from beauty -> D
      fMCSources->Fill(1);
      fMCSources->Fill(3);
    }	
    if(mcSource[3]&&mcSource[1]){ // is from beauty ->B
      fMCSources->Fill(1);
      fMCSources->Fill(4);
    }
  }

  if(sel==2){
    fMassVsdPhiHistKaon[ind]->Fill(invMass,deltaPhi);
    fMassVsdEtaHistKaon[ind]->Fill(invMass,deltaPhi,deltaEta);
    fKaonOrigin->Fill(0);
    if(mcSource[2]&&mcSource[0]){ // is from charm ->D
      fKaonOrigin->Fill(1);
      fKaonOrigin->Fill(2);
    }	
    if(mcSource[2]&&mcSource[1]){ // is from beauty -> D
      fKaonOrigin->Fill(1);
      fKaonOrigin->Fill(3);
    }	
    if(mcSource[3]&&mcSource[1]){ // is from beauty ->B
      fKaonOrigin->Fill(1);
      fKaonOrigin->Fill(4);
    }
  }
  if(sel==3){
    fMassVsdPhiHistKshort[ind]->Fill(invMass,deltaPhi);
    fMassVsdEtaHistKshort[ind]->Fill(invMass,deltaPhi,deltaEta);
    fK0Origin->Fill(0);
    if(mcSource[2]&&mcSource[0]){ // is from charm ->D
      fK0Origin->Fill(1);
      fK0Origin->Fill(2);
    }	
    if(mcSource[2]&&mcSource[1]){ // is from beauty -> D
      fK0Origin->Fill(1);
      fK0Origin->Fill(3);
    }
    if(mcSource[3]&&mcSource[1]){ // is from beauty ->B
      fK0Origin->Fill(1);
      fK0Origin->Fill(4);
    }
	  
  }

  return;
}










