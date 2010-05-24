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

//*************************************************************************
// Class AliAnalysisTaskSEDplus
// AliAnalysisTaskSE for the D+ candidates Invariant Mass Histogram and 
//comparison of heavy-flavour decay candidates
// to MC truth (kinematics stored in the AOD)
// Authors: Renu Bala, bala@to.infn.it
// F. Prino, prino@to.infn.it
// G. Ortona, ortona@to.infn.it
/////////////////////////////////////////////////////////////

#include <TClonesArray.h>
#include <TNtuple.h>
#include <TCanvas.h>
#include <TList.h>
#include <TString.h>
#include <TH1F.h>
#include <TDatabasePDG.h>

#include "AliAnalysisManager.h"
#include "AliRDHFCutsDplustoKpipi.h"
#include "AliAODHandler.h"
#include "AliAODEvent.h"
#include "AliAODVertex.h"
#include "AliAODTrack.h"
#include "AliAODMCHeader.h"
#include "AliAODMCParticle.h"
#include "AliAODRecoDecayHF3Prong.h"
#include "AliAnalysisVertexingHF.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisTaskSEDplus.h"

ClassImp(AliAnalysisTaskSEDplus)


//________________________________________________________________________
AliAnalysisTaskSEDplus::AliAnalysisTaskSEDplus():
AliAnalysisTaskSE(),
fOutput(0), 
fHistNEvents(0),
fNtupleDplus(0),
fUpmasslimit(1.965),
fLowmasslimit(1.765),
fNPtBins(0),
fBinWidth(0.002),
fListCuts(0),
fRDCutsProduction(0),
fRDCutsAnalysis(0),
fFillNtuple(kFALSE),
fReadMC(kFALSE),
fDoLS(kFALSE)
{
   // Default constructor
}

//________________________________________________________________________
AliAnalysisTaskSEDplus::AliAnalysisTaskSEDplus(const char *name,AliRDHFCutsDplustoKpipi *dpluscutsana,AliRDHFCutsDplustoKpipi *dpluscutsprod,Bool_t fillNtuple):
AliAnalysisTaskSE(name),
fOutput(0),
fHistNEvents(0),
fNtupleDplus(0),
fUpmasslimit(1.965),
fLowmasslimit(1.765),
fNPtBins(0),
fBinWidth(0.002),
fListCuts(0),
fRDCutsProduction(dpluscutsprod),
fRDCutsAnalysis(dpluscutsana),
fFillNtuple(fillNtuple),
fReadMC(kFALSE),
fDoLS(kFALSE)
{
  // 
  // Standrd constructor
  //
  //Double_t ptlim[5]={0.,2.,3.,5,9999999.};
   //SetPtBinLimit(5, ptlim);
  SetPtBinLimit(fRDCutsAnalysis->GetNPtBins()+1,fRDCutsAnalysis->GetPtBinLimits());
  // Default constructor
   // Output slot #1 writes into a TList container
  DefineOutput(1,TList::Class());  //My private output
 // Output slot #2 writes cut to private output
  //  DefineOutput(2,AliRDHFCutsDplustoKpipi::Class());
  DefineOutput(2,TList::Class());
  if(fFillNtuple){
    // Output slot #3 writes into a TNtuple container
    DefineOutput(3,TNtuple::Class());  //My private output
  }
}

//________________________________________________________________________
AliAnalysisTaskSEDplus::~AliAnalysisTaskSEDplus()
{
  //
  // Destructor
  //
  if (fOutput) {
    delete fOutput;
    fOutput = 0;
  }

  if (fListCuts) {
    delete fListCuts;
    fListCuts = 0;
  }

  if(fRDCutsProduction){
    delete fRDCutsProduction;
    fRDCutsProduction = 0;
  }

   if(fRDCutsAnalysis){
    delete fRDCutsAnalysis;
    fRDCutsAnalysis = 0;
  }

}  
//_________________________________________________________________
void  AliAnalysisTaskSEDplus::SetMassLimits(Float_t range){
  // set invariant mass limits
  Float_t bw=GetBinWidth();
  fUpmasslimit = 1.865+range;
  fLowmasslimit = 1.865-range;
  SetBinWidth(bw);
}
//_________________________________________________________________
void  AliAnalysisTaskSEDplus::SetMassLimits(Float_t lowlimit, Float_t uplimit){
  // set invariant mass limits
  if(uplimit>lowlimit)
    {
      Float_t bw=GetBinWidth();
      fUpmasslimit = lowlimit;
      fLowmasslimit = uplimit;
      SetBinWidth(bw);
    }
}
//________________________________________________________________________
void AliAnalysisTaskSEDplus::SetPtBinLimit(Int_t n, Float_t* lim){
  // define pt bins for analysis
  if(n>kMaxPtBins){
    printf("Max. number of Pt bins = %d\n",kMaxPtBins);
    fNPtBins=kMaxPtBins;
    fArrayBinLimits[0]=0.;
    fArrayBinLimits[1]=2.;
    fArrayBinLimits[2]=3.;
    fArrayBinLimits[3]=5.;
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
    for(Int_t i=0; i<fNPtBins+1; i++) printf(" Bin%d = %8.2f-%8.2f\n",i,fArrayBinLimits[i],fArrayBinLimits[i+1]);    
  }
}
//________________________________________________________________
void AliAnalysisTaskSEDplus::SetBinWidth(Float_t w){
  Float_t width=w;
  Int_t nbins=(Int_t)((fUpmasslimit-fLowmasslimit)/width+0.5);
  Int_t missingbins=4-nbins%4;
  nbins=nbins+missingbins;
  width=(fUpmasslimit-fLowmasslimit)/nbins;
  if(missingbins!=0){
    printf("AliAnalysisTaskSEDplus::SetBinWidth: W-bin width of %f will produce histograms not rebinnable by 4. New width set to %f\n",w,width);
  }
  else{
    if(fDebug>1) printf("AliAnalysisTaskSEDplus::SetBinWidth: width set to %f\n",width);
  }
  fBinWidth=width;
}
//_________________________________________________________________
Double_t  AliAnalysisTaskSEDplus::GetPtBinLimit(Int_t ibin){
  // get pt bin limit
  if(ibin>fNPtBins)return -1;
  return fArrayBinLimits[ibin];
} 
//_________________________________________________________________
Int_t AliAnalysisTaskSEDplus::GetNBinsHistos(){
  return (Int_t)((fUpmasslimit-fLowmasslimit)/fBinWidth+0.5);
}
//_________________________________________________________________
void AliAnalysisTaskSEDplus::LSAnalysis(TClonesArray *arrayOppositeSign,TClonesArray *arrayLikeSign,AliAODEvent *aod,AliAODVertex *vtx1, Int_t nDplusOS){
  //
  //
  // Fill the Like Sign histograms
  //

  //count pos/neg tracks
  Int_t nPosTrks=0,nNegTrks=0;
 //counter for particles passing single particle cuts
  Int_t nspcplus=0;
  Int_t nspcminus=0;

  for(Int_t it=0;it<aod->GetNumberOfTracks();it++) {
    AliAODTrack *track = aod->GetTrack(it);
    if(track->Charge()>0){
      nPosTrks++;
      if(track->Pt()>=0.4){
	nspcplus++;
      }
    }
    if(track->Charge()<0)
      {
	nNegTrks++;
	if(track->Pt()>=0.4){
	  nspcminus++;
	}
      }
  }

  Int_t nOStriplets = arrayOppositeSign->GetEntriesFast();

  Int_t nDplusLS=0;
  Int_t nLikeSign = arrayLikeSign->GetEntriesFast();
  Int_t index; 

  for(Int_t iLikeSign = 0; iLikeSign < nLikeSign; iLikeSign++) {
    AliAODRecoDecayHF3Prong *d = (AliAODRecoDecayHF3Prong*)arrayLikeSign->UncheckedAt(iLikeSign);
    Bool_t unsetvtx=kFALSE;
    if(!d->GetOwnPrimaryVtx()) {
      d->SetOwnPrimaryVtx(vtx1); // needed to compute all variables
      unsetvtx=kTRUE;
    }
    if(fRDCutsProduction->IsSelected(d,AliRDHFCuts::kCandidate))nDplusLS++;
    if(unsetvtx) d->UnsetOwnPrimaryVtx();
  }

 Float_t wei2=0;
 if(nLikeSign!=0)wei2 = (Float_t)nOStriplets/(Float_t)nLikeSign;
 Float_t wei3=0;
 if(nDplusLS!=0)wei3 = (Float_t)nDplusOS/(Float_t)nDplusLS;

 // loop over like sign candidates
  for(Int_t iLikeSign = 0; iLikeSign < nLikeSign; iLikeSign++) {
    AliAODRecoDecayHF3Prong *d = (AliAODRecoDecayHF3Prong*)arrayLikeSign->UncheckedAt(iLikeSign);
    Bool_t unsetvtx=kFALSE;
    if(!d->GetOwnPrimaryVtx()) {
      d->SetOwnPrimaryVtx(vtx1); // needed to compute all variables
      unsetvtx=kTRUE;
    }
 
    if(fRDCutsProduction->IsSelected(d,AliRDHFCuts::kCandidate)){

      //set tight cuts values
      Int_t iPtBin=-1;
      Double_t ptCand = d->Pt();
      for(Int_t ibin=0;ibin<fNPtBins&&iPtBin<0&&ptCand>fArrayBinLimits[0]&&ptCand<fArrayBinLimits[fNPtBins];ibin++){
      	if(ptCand<fArrayBinLimits[ibin+1])iPtBin=ibin;
      }
      
      if(iPtBin<0){
	return;
      }

      Int_t passTightCuts=fRDCutsAnalysis->IsSelected(d,AliRDHFCuts::kCandidate);

      Int_t sign= d->GetCharge();
      Float_t wei=1;
      Float_t wei4=1;
      if(sign>0&&nPosTrks>2&&nspcplus>2) {  //wei* should be automatically protected, since to get a triplet there must be at least 3 good tracks in the event
	
      	wei=3.*(Float_t)nNegTrks/((Float_t)nPosTrks-2.);
	wei4=3.*(Float_t)nspcminus/((Float_t)nspcplus-2.);
      }
      
      if(sign<0&&nNegTrks>2&&nspcminus>2){     
 	wei=3.*(Float_t)nPosTrks/((Float_t)nNegTrks-2.);
	wei4=3.*(Float_t)nspcplus/((Float_t)nspcminus-2.);

      }

      Float_t invMass = d->InvMassDplus();
      Double_t dlen=d->DecayLength();
      Double_t cosp=d->CosPointingAngle();
      Double_t sumD02=d->Getd0Prong(0)*d->Getd0Prong(0)+d->Getd0Prong(1)*d->Getd0Prong(1)+d->Getd0Prong(2)*d->Getd0Prong(2);
     Double_t dca=d->GetDCA();   
     Double_t ptmax=0;
     for(Int_t i=0;i<3;i++){
       if(d->PtProng(i)>ptmax)ptmax=d->PtProng(i);
     }
     
      index=GetLSHistoIndex(iPtBin);
      fMassHistLS[index]->Fill(invMass,wei);
      fMassHistLS[index+1]->Fill(invMass);
      fMassHistLS[index+2]->Fill(invMass,wei2);
      fMassHistLS[index+3]->Fill(invMass,wei3);
      fMassHistLS[index+4]->Fill(invMass,wei4);
      
      Int_t indexcut=GetHistoIndex(iPtBin);
      fCosPHistLS[indexcut]->Fill(cosp);
      fDLenHistLS[indexcut]->Fill(dlen);
      fSumd02HistLS[indexcut]->Fill(sumD02);
      fPtMaxHistLS[indexcut]->Fill(ptmax);
      fDCAHistLS[indexcut]->Fill(dca);
      
      if(passTightCuts==1){
	fMassHistLSTC[index]->Fill(invMass,wei);
	fMassHistLSTC[index+1]->Fill(invMass);
	fMassHistLSTC[index+2]->Fill(invMass,wei2);
	fMassHistLSTC[index+3]->Fill(invMass,wei3);
	fMassHistLSTC[index+4]->Fill(invMass,wei4);
      }
    }
    if(unsetvtx) d->UnsetOwnPrimaryVtx();
  }
  
  //printf("------------ N. of positive tracks in Event ----- %d \n", nPosTrks);
  //printf("------------ N. of negative tracks in Event ----- %d \n", nNegTrks);

  //  printf("LS analysis...done\n");

}


//__________________________________________
void AliAnalysisTaskSEDplus::Init(){
  //
  // Initialization
  //
  if(fDebug > 1) printf("AnalysisTaskSEDplus::Init() \n");
  
  //PostData(2,fRDCutsloose);//we should then put those cuts in a tlist if we have more than 1
  fListCuts=new TList();
  AliRDHFCutsDplustoKpipi *production = new AliRDHFCutsDplustoKpipi();
  production=fRDCutsProduction;
  AliRDHFCutsDplustoKpipi *analysis = new AliRDHFCutsDplustoKpipi();
  analysis=fRDCutsAnalysis;
  
  fListCuts->Add(production);
  fListCuts->Add(analysis);
  PostData(2,fListCuts);
  
  return;
}

//________________________________________________________________________
void AliAnalysisTaskSEDplus::UserCreateOutputObjects()
{
  // Create the output container
  //
  if(fDebug > 1) printf("AnalysisTaskSEDplus::UserCreateOutputObjects() \n");

  // Several histograms are more conveniently managed in a TList
  fOutput = new TList();
  fOutput->SetOwner();
  fOutput->SetName("OutputHistos");

  TString hisname;
  Int_t index=0;
  Int_t indexLS=0;
  Int_t nbins=GetNBinsHistos();
  for(Int_t i=0;i<fNPtBins;i++){

    index=GetHistoIndex(i);
    indexLS=GetLSHistoIndex(i);

    hisname.Form("hMassPt%d",i);
    fMassHist[index]=new TH1F(hisname.Data(),hisname.Data(),nbins,fLowmasslimit,fUpmasslimit);
    fMassHist[index]->Sumw2();
    hisname.Form("hCosPAllPt%d",i);
    fCosPHist[index]=new TH1F(hisname.Data(),hisname.Data(),nbins,0.5,1.);
    fCosPHist[index]->Sumw2();
    hisname.Form("hDLenAllPt%d",i);
    fDLenHist[index]=new TH1F(hisname.Data(),hisname.Data(),nbins,0.,0.5);
    fDLenHist[index]->Sumw2();
    hisname.Form("hSumd02AllPt%d",i);
    fSumd02Hist[index]=new TH1F(hisname.Data(),hisname.Data(),nbins,0.,1.);
    fSumd02Hist[index]->Sumw2();
    hisname.Form("hSigVertAllPt%d",i);
    fSigVertHist[index]=new TH1F(hisname.Data(),hisname.Data(),nbins,0.,0.1);
    fSigVertHist[index]->Sumw2();
    hisname.Form("hPtMaxAllPt%d",i);
    fPtMaxHist[index]=new TH1F(hisname.Data(),hisname.Data(),nbins,0.5,5.);
    fPtMaxHist[index]->Sumw2();

    hisname.Form("hDCAAllPt%d",i);
    fDCAHist[index]=new TH1F(hisname.Data(),hisname.Data(),nbins,0.,0.1);
    fDCAHist[index]->Sumw2();



    hisname.Form("hMassPt%dTC",i);
    fMassHistTC[index]=new TH1F(hisname.Data(),hisname.Data(),nbins,fLowmasslimit,fUpmasslimit);
    fMassHistTC[index]->Sumw2();




    
    hisname.Form("hCosPAllPt%dLS",i);
    fCosPHistLS[index]=new TH1F(hisname.Data(),hisname.Data(),nbins,0.5,1.);
    fCosPHistLS[index]->Sumw2();
    hisname.Form("hDLenAllPt%dLS",i);
    fDLenHistLS[index]=new TH1F(hisname.Data(),hisname.Data(),nbins,0.,0.5);
    fDLenHistLS[index]->Sumw2();
    hisname.Form("hSumd02AllPt%dLS",i);
    fSumd02HistLS[index]=new TH1F(hisname.Data(),hisname.Data(),nbins,0.,1.);
    fSumd02HistLS[index]->Sumw2();
    hisname.Form("hSigVertAllPt%dLS",i);
    fSigVertHistLS[index]=new TH1F(hisname.Data(),hisname.Data(),nbins,0.,0.1);
    fSigVertHistLS[index]->Sumw2();
    hisname.Form("hPtMaxAllPt%dLS",i);
    fPtMaxHistLS[index]=new TH1F(hisname.Data(),hisname.Data(),nbins,0.5,5.);
    fPtMaxHistLS[index]->Sumw2();
    
    hisname.Form("hDCAAllPt%dLS",i);
    fDCAHistLS[index]=new TH1F(hisname.Data(),hisname.Data(),nbins,0.,0.1);
    fDCAHistLS[index]->Sumw2();
    
    hisname.Form("hLSPt%dLC",i);
    fMassHistLS[indexLS] = new TH1F(hisname.Data(),hisname.Data(),nbins,fLowmasslimit,fUpmasslimit);
    fMassHistLS[indexLS]->Sumw2();
    
    hisname.Form("hLSPt%dTC",i);
    fMassHistLSTC[indexLS] = new TH1F(hisname.Data(),hisname.Data(),nbins,fLowmasslimit,fUpmasslimit);
    fMassHistLSTC[indexLS]->Sumw2();


    
    index=GetSignalHistoIndex(i);    
    indexLS++;
    hisname.Form("hSigPt%d",i);
    fMassHist[index]=new TH1F(hisname.Data(),hisname.Data(),nbins,fLowmasslimit,fUpmasslimit);
    fMassHist[index]->Sumw2();
    hisname.Form("hCosPSigPt%d",i);
    fCosPHist[index]=new TH1F(hisname.Data(),hisname.Data(),nbins,0.5,1.);
    fCosPHist[index]->Sumw2();
    hisname.Form("hDLenSigPt%d",i);
    fDLenHist[index]=new TH1F(hisname.Data(),hisname.Data(),nbins,0.,0.5);
    fDLenHist[index]->Sumw2();
    hisname.Form("hSumd02SigPt%d",i);
    fSumd02Hist[index]=new TH1F(hisname.Data(),hisname.Data(),nbins,0.,1.);
    fSumd02Hist[index]->Sumw2();
    hisname.Form("hSigVertSigPt%d",i);
    fSigVertHist[index]=new TH1F(hisname.Data(),hisname.Data(),nbins,0.,0.1);
    fSigVertHist[index]->Sumw2();
    hisname.Form("hPtMaxSigPt%d",i);
    fPtMaxHist[index]=new TH1F(hisname.Data(),hisname.Data(),nbins,0.5,5.);
    fPtMaxHist[index]->Sumw2();    

    hisname.Form("hDCASigPt%d",i);
    fDCAHist[index]=new TH1F(hisname.Data(),hisname.Data(),nbins,0.,0.1);
    fDCAHist[index]->Sumw2();    


    hisname.Form("hSigPt%dTC",i);
    fMassHistTC[index]=new TH1F(hisname.Data(),hisname.Data(),nbins,fLowmasslimit,fUpmasslimit);
    fMassHistTC[index]->Sumw2();

    hisname.Form("hLSPt%dLCnw",i);
    fMassHistLS[indexLS]=new TH1F(hisname.Data(),hisname.Data(),nbins,fLowmasslimit,fUpmasslimit);
    fMassHistLS[indexLS]->Sumw2();
    hisname.Form("hLSPt%dTCnw",i);
    fMassHistLSTC[indexLS]=new TH1F(hisname.Data(),hisname.Data(),nbins,fLowmasslimit,fUpmasslimit);
    fMassHistLSTC[indexLS]->Sumw2();


    
    hisname.Form("hCosPSigPt%dLS",i);
    fCosPHistLS[index]=new TH1F(hisname.Data(),hisname.Data(),nbins,0.5,1.);
    fCosPHistLS[index]->Sumw2();
    hisname.Form("hDLenSigPt%dLS",i);
    fDLenHistLS[index]=new TH1F(hisname.Data(),hisname.Data(),nbins,0.,0.5);
    fDLenHistLS[index]->Sumw2();
    hisname.Form("hSumd02SigPt%dLS",i);
    fSumd02HistLS[index]=new TH1F(hisname.Data(),hisname.Data(),nbins,0.,1.);
    fSumd02HistLS[index]->Sumw2();
    hisname.Form("hSigVertSigPt%dLS",i);
    fSigVertHistLS[index]=new TH1F(hisname.Data(),hisname.Data(),nbins,0.,0.1);
    fSigVertHistLS[index]->Sumw2();
    hisname.Form("hPtMaxSigPt%dLS",i);
    fPtMaxHistLS[index]=new TH1F(hisname.Data(),hisname.Data(),nbins,0.5,5.);
    fPtMaxHistLS[index]->Sumw2();

    hisname.Form("hDCASigPt%dLS",i);
    fDCAHistLS[index]=new TH1F(hisname.Data(),hisname.Data(),nbins,0.,0.1);
    fDCAHistLS[index]->Sumw2();
    


    index=GetBackgroundHistoIndex(i); 
    indexLS++;
    hisname.Form("hBkgPt%d",i);
    fMassHist[index]=new TH1F(hisname.Data(),hisname.Data(),nbins,fLowmasslimit,fUpmasslimit);
    fMassHist[index]->Sumw2();
    hisname.Form("hCosPBkgPt%d",i);
    fCosPHist[index]=new TH1F(hisname.Data(),hisname.Data(),nbins,0.5,1.);
    fCosPHist[index]->Sumw2();
    hisname.Form("hDLenBkgPt%d",i);
    fDLenHist[index]=new TH1F(hisname.Data(),hisname.Data(),nbins,0.,0.5);
    fDLenHist[index]->Sumw2();
    hisname.Form("hSumd02BkgPt%d",i);
    fSumd02Hist[index]=new TH1F(hisname.Data(),hisname.Data(),nbins,0.,1.);
    fSumd02Hist[index]->Sumw2();
    hisname.Form("hSigVertBkgPt%d",i);
    fSigVertHist[index]=new TH1F(hisname.Data(),hisname.Data(),nbins,0.,0.1);
    fSigVertHist[index]->Sumw2();
    hisname.Form("hPtMaxBkgPt%d",i);
    fPtMaxHist[index]=new TH1F(hisname.Data(),hisname.Data(),nbins,0.5,5.);
    fPtMaxHist[index]->Sumw2();

    hisname.Form("hDCABkgPt%d",i);
    fDCAHist[index]=new TH1F(hisname.Data(),hisname.Data(),nbins,0.,0.1);
    fDCAHist[index]->Sumw2();


    hisname.Form("hBkgPt%dTC",i);
    fMassHistTC[index]=new TH1F(hisname.Data(),hisname.Data(),nbins,fLowmasslimit,fUpmasslimit);
    fMassHistTC[index]->Sumw2();

    hisname.Form("hLSPt%dLCntrip",i);
    fMassHistLS[indexLS]=new TH1F(hisname.Data(),hisname.Data(),nbins,fLowmasslimit,fUpmasslimit);
    fMassHistLS[indexLS]->Sumw2();
    hisname.Form("hLSPt%dTCntrip",i);
    fMassHistLSTC[indexLS]=new TH1F(hisname.Data(),hisname.Data(),nbins,fLowmasslimit,fUpmasslimit);
    fMassHistLSTC[indexLS]->Sumw2();

    
    hisname.Form("hCosPBkgPt%dLS",i);
    fCosPHistLS[index]=new TH1F(hisname.Data(),hisname.Data(),nbins,0.5,1.);
    fCosPHistLS[index]->Sumw2();
    hisname.Form("hDLenBkgPt%dLS",i);
    fDLenHistLS[index]=new TH1F(hisname.Data(),hisname.Data(),nbins,0.,0.5);
    fDLenHistLS[index]->Sumw2();
    hisname.Form("hSumd02BkgPt%dLS",i);
    fSumd02HistLS[index]=new TH1F(hisname.Data(),hisname.Data(),nbins,0.,1.);
    fSumd02HistLS[index]->Sumw2();
    hisname.Form("hSigVertBkgPt%dLS",i);
    fSigVertHistLS[index]=new TH1F(hisname.Data(),hisname.Data(),nbins,0.,0.1);
    fSigVertHistLS[index]->Sumw2();
    hisname.Form("hPtMaxBkgPt%dLS",i);
    fPtMaxHistLS[index]=new TH1F(hisname.Data(),hisname.Data(),nbins,0.5,5.);
    fPtMaxHistLS[index]->Sumw2();
    hisname.Form("hDCABkgPt%dLS",i);
    fDCAHistLS[index]=new TH1F(hisname.Data(),hisname.Data(),nbins,0.,0.1);
    fDCAHistLS[index]->Sumw2();
    

    indexLS++;
    hisname.Form("hLSPt%dLCntripsinglecut",i);
    fMassHistLS[indexLS]=new TH1F(hisname.Data(),hisname.Data(),nbins,fLowmasslimit,fUpmasslimit);
    fMassHistLS[indexLS]->Sumw2();
    hisname.Form("hLSPt%dTCntripsinglecut",i);
    fMassHistLSTC[indexLS]=new TH1F(hisname.Data(),hisname.Data(),nbins,fLowmasslimit,fUpmasslimit);
    fMassHistLSTC[indexLS]->Sumw2();

    indexLS++;
    hisname.Form("hLSPt%dLCspc",i);
    fMassHistLS[indexLS]=new TH1F(hisname.Data(),hisname.Data(),nbins,fLowmasslimit,fUpmasslimit);
    fMassHistLS[indexLS]->Sumw2();
    hisname.Form("hLSPt%dTCspc",i);
    fMassHistLSTC[indexLS]=new TH1F(hisname.Data(),hisname.Data(),nbins,fLowmasslimit,fUpmasslimit);
    fMassHistLSTC[indexLS]->Sumw2();
  }
  
  for(Int_t i=0; i<3*fNPtBins; i++){
    fOutput->Add(fMassHist[i]);
     fOutput->Add(fCosPHist[i]);
    fOutput->Add(fDLenHist[i]);
    fOutput->Add(fSumd02Hist[i]);
    fOutput->Add(fSigVertHist[i]);
    fOutput->Add(fPtMaxHist[i]);
    fOutput->Add(fDCAHist[i]);
    fOutput->Add(fMassHistTC[i]);
  }
  for(Int_t i=0; i<3*fNPtBins&&fDoLS; i++){
    fOutput->Add(fCosPHistLS[i]);
    fOutput->Add(fDLenHistLS[i]);
    fOutput->Add(fSumd02HistLS[i]);
    fOutput->Add(fSigVertHistLS[i]);
    fOutput->Add(fPtMaxHistLS[i]);  
    fOutput->Add(fDCAHistLS[i]);  
  }
  for(Int_t i=0; i<5*fNPtBins&&fDoLS; i++){
    fOutput->Add(fMassHistLS[i]);
    fOutput->Add(fMassHistLSTC[i]);
  }
  

  fHistNEvents = new TH1F("fHistNEvents", "Number of processed events; ; Events",3,-1.5,1.5);
  fHistNEvents->Sumw2();
  fHistNEvents->SetMinimum(0);
  fOutput->Add(fHistNEvents);
  
  

  if(fFillNtuple){
    OpenFile(2); // 2 is the slot number of the ntuple
   
    fNtupleDplus = new TNtuple("fNtupleDplus","D +","pdg:Px:Py:Pz:PtTrue:VxTrue:VyTrue:VzTrue:Ptpi:PtK:Ptpi2:PtRec:PointingAngle:DecLeng:VxRec:VyRec:VzRec:InvMass:sigvert:d0Pi:d0K:d0Pi2:dca:d0square");  
    
  }
  
  return;
}

//________________________________________________________________________
void AliAnalysisTaskSEDplus::UserExec(Option_t */*option*/)
{
  // Execute analysis for current event:
  // heavy flavor candidates association to MC truth

  AliAODEvent *aod = dynamic_cast<AliAODEvent*> (InputEvent());
  fHistNEvents->Fill(0); // count event
  // Post the data already here
  PostData(1,fOutput);
  
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
  } else {
    array3Prong=(TClonesArray*)aod->GetList()->FindObject("Charm3Prong");
    arrayLikeSign=(TClonesArray*)aod->GetList()->FindObject("LikeSign3Prong");
  }

  if(!array3Prong) {
    printf("AliAnalysisTaskSEDplus::UserExec: Charm3Prong branch not found!\n");
    return;
  }
  if(!arrayLikeSign) {
    printf("AliAnalysisTaskSEDplus::UserExec: LikeSign3Prong branch not found!\n");
    return;
  }

 
  TClonesArray *arrayMC=0;
  AliAODMCHeader *mcHeader=0;

  // AOD primary vertex
  AliAODVertex *vtx1 = (AliAODVertex*)aod->GetPrimaryVertex();
  //    vtx1->Print();
  
  // load MC particles
  if(fReadMC){
    
    arrayMC =  (TClonesArray*)aod->GetList()->FindObject(AliAODMCParticle::StdBranchName());
    if(!arrayMC) {
      printf("AliAnalysisTaskSEDplus::UserExec: MC particles branch not found!\n");
      //    return;
    }
    
  // load MC header
    mcHeader =  (AliAODMCHeader*)aod->GetList()->FindObject(AliAODMCHeader::StdBranchName());
    if(!mcHeader) {
    printf("AliAnalysisTaskSEDplus::UserExec: MC header branch not found!\n");
    return;
    }
  }
  
  Int_t n3Prong = array3Prong->GetEntriesFast();
  if(fDebug>1) printf("Number of D+->Kpipi: %d\n",n3Prong);
  
  
  Int_t nOS=0;
  Int_t index;
  Int_t pdgDgDplustoKpipi[3]={321,211,211};
  // Double_t cutsDplus[12]={0.2,0.4,0.4,0.,0.,0.01,0.06,0.02,0.,0.85,0.,10000000000.};//TO REMOVE
  //Double_t *cutsDplus = new (Double_t*)fRDCuts->GetCuts();
  for (Int_t i3Prong = 0; i3Prong < n3Prong; i3Prong++) {
    AliAODRecoDecayHF3Prong *d = (AliAODRecoDecayHF3Prong*)array3Prong->UncheckedAt(i3Prong);
    
    
    Bool_t unsetvtx=kFALSE;
    if(!d->GetOwnPrimaryVtx()){
      d->SetOwnPrimaryVtx(vtx1);
      unsetvtx=kTRUE;
    }

    if(fRDCutsProduction->IsSelected(d,AliRDHFCuts::kCandidate)) {
      Int_t iPtBin = -1;
      Double_t ptCand = d->Pt();

      for(Int_t ibin=0;ibin<fNPtBins&&iPtBin<0&&ptCand>fArrayBinLimits[0]&&ptCand<fArrayBinLimits[fNPtBins];ibin++){
	if(ptCand<fArrayBinLimits[ibin+1])iPtBin=ibin;
      }
      
      Int_t passTightCuts=fRDCutsAnalysis->IsSelected(d,AliRDHFCuts::kCandidate);

      Int_t labDp=-1;
      Float_t deltaPx=0.;
      Float_t deltaPy=0.;
      Float_t deltaPz=0.;
      Float_t truePt=0.;
      Float_t xDecay=0.;
      Float_t yDecay=0.;
      Float_t zDecay=0.;
      Float_t pdgCode=-2;
      if(fReadMC){
	labDp = d->MatchToMC(411,arrayMC,3,pdgDgDplustoKpipi);
	if(labDp>=0){
	  AliAODMCParticle *partDp = (AliAODMCParticle*)arrayMC->At(labDp);
	  AliAODMCParticle *dg0 = (AliAODMCParticle*)arrayMC->At(partDp->GetDaughter(0));
	  deltaPx=partDp->Px()-d->Px();
	  deltaPy=partDp->Py()-d->Py();
	  deltaPz=partDp->Pz()-d->Pz();
	  truePt=partDp->Pt();
	  xDecay=dg0->Xv();	  
	  yDecay=dg0->Yv();	  
	  zDecay=dg0->Zv();
	  pdgCode=TMath::Abs(partDp->GetPdgCode());
	}else{
	  pdgCode=-1;
	}
      }
      Double_t invMass=d->InvMassDplus();

      Float_t tmp[24];
      if(fFillNtuple){  	  
	tmp[0]=pdgCode;
	tmp[1]=deltaPx;
	tmp[2]=deltaPy;
	tmp[3]=deltaPz;
	tmp[4]=truePt;
	tmp[5]=xDecay;	  
	tmp[6]=yDecay;	  
	tmp[7]=zDecay;	  
	tmp[8]=d->PtProng(0);
	tmp[9]=d->PtProng(1);
	tmp[10]=d->PtProng(2);
	tmp[11]=d->Pt();
	tmp[12]=d->CosPointingAngle();
	tmp[13]=d->DecayLength();
	tmp[14]=d->Xv();
	tmp[15]=d->Yv();
	tmp[16]=d->Zv();
	tmp[17]=d->InvMassDplus();
	tmp[18]=d->GetSigmaVert();
	tmp[19]=d->Getd0Prong(0);
	tmp[20]=d->Getd0Prong(1);
	tmp[21]=d->Getd0Prong(2);
	tmp[22]=d->GetDCA();
	tmp[23]=d->Prodd0d0(); 
	fNtupleDplus->Fill(tmp);
	PostData(3,fNtupleDplus);
      }
      Double_t dlen=d->DecayLength();
      Double_t cosp=d->CosPointingAngle();
      Double_t sumD02=d->Getd0Prong(0)*d->Getd0Prong(0)+d->Getd0Prong(1)*d->Getd0Prong(1)+d->Getd0Prong(2)*d->Getd0Prong(2);
      Double_t dca=d->GetDCA();      
Double_t ptmax=0;
      for(Int_t i=0;i<3;i++){
	if(d->PtProng(i)>ptmax)ptmax=d->PtProng(i);
      }
      if(iPtBin>=0){
      
	index=GetHistoIndex(iPtBin);
	fMassHist[index]->Fill(invMass);
	fCosPHist[index]->Fill(cosp);
	fDLenHist[index]->Fill(dlen);
	fSumd02Hist[index]->Fill(sumD02);
	fPtMaxHist[index]->Fill(ptmax);
	fDCAHist[index]->Fill(dca);
	
	if(passTightCuts==1){
	  fMassHistTC[index]->Fill(invMass);
	}
	
	if(fReadMC){
	  if(labDp>=0) {
	    index=GetSignalHistoIndex(iPtBin);
	    fMassHist[index]->Fill(invMass);
	    fCosPHist[index]->Fill(cosp);
	    fDLenHist[index]->Fill(dlen);
	    fSumd02Hist[index]->Fill(sumD02);
	    fPtMaxHist[index]->Fill(ptmax);
	    fDCAHist[index]->Fill(dca);
	    if(passTightCuts==1){
	      fMassHistTC[index]->Fill(invMass);

	    }
	    
	  }else{
	    index=GetBackgroundHistoIndex(iPtBin);
	    fMassHist[index]->Fill(invMass);
	    fCosPHist[index]->Fill(cosp);
	    fDLenHist[index]->Fill(dlen);
	    fSumd02Hist[index]->Fill(sumD02);
	    fPtMaxHist[index]->Fill(ptmax);
	    fDCAHist[index]->Fill(dca);
	    if(passTightCuts==1){
	      fMassHistTC[index]->Fill(invMass);

	    }	
	  }
	}
      }
    }
    if(unsetvtx) d->UnsetOwnPrimaryVtx();
  }
  
  //start LS analysis
  if(fDoLS && arrayLikeSign) LSAnalysis(array3Prong,arrayLikeSign,aod,vtx1,nOS);
  
  PostData(1,fOutput);    
  return;
}



//________________________________________________________________________
void AliAnalysisTaskSEDplus::Terminate(Option_t */*option*/)
{
  // Terminate analysis
  //
  if(fDebug > 1) printf("AnalysisTaskSEDplus: Terminate() \n");

  fOutput = dynamic_cast<TList*> (GetOutputData(1));
  if (!fOutput) {     
    printf("ERROR: fOutput not available\n");
    return;
  }
 fHistNEvents = dynamic_cast<TH1F*>(fOutput->FindObject("fHistNEvents"));

 TString hisname;
 Int_t index=0;
 

 Int_t indexLS=0;
 for(Int_t i=0;i<fNPtBins;i++){
    index=GetHistoIndex(i);
    if(fDoLS)indexLS=GetLSHistoIndex(i);
    hisname.Form("hMassPt%d",i);
    fMassHist[index]=dynamic_cast<TH1F*>(fOutput->FindObject(hisname.Data()));
     hisname.Form("hCosPAllPt%d",i);
    fCosPHist[index]=dynamic_cast<TH1F*>(fOutput->FindObject(hisname.Data()));
     hisname.Form("hDLenAllPt%d",i);
    fDLenHist[index]=dynamic_cast<TH1F*>(fOutput->FindObject(hisname.Data()));
     hisname.Form("hSumd02AllPt%d",i);
     fSumd02Hist[index]=dynamic_cast<TH1F*>(fOutput->FindObject(hisname.Data()));
     hisname.Form("hSigVertAllPt%d",i);
     fSigVertHist[index]=dynamic_cast<TH1F*>(fOutput->FindObject(hisname.Data()));
     hisname.Form("hPtMaxAllPt%d",i);
     fPtMaxHist[index]=dynamic_cast<TH1F*>(fOutput->FindObject(hisname.Data()));
     hisname.Form("hDCAAllPt%d",i);
     fDCAHist[index]=dynamic_cast<TH1F*>(fOutput->FindObject(hisname.Data()));
     hisname.Form("hMassPt%dTC",i);
    fMassHistTC[index]=dynamic_cast<TH1F*>(fOutput->FindObject(hisname.Data()));
    if(fDoLS){
      hisname.Form("hLSPt%dLC",i);
      fMassHistLS[indexLS]=dynamic_cast<TH1F*>(fOutput->FindObject(hisname.Data()));
      hisname.Form("hCosPAllPt%dLS",i);
      fCosPHistLS[index]=dynamic_cast<TH1F*>(fOutput->FindObject(hisname.Data()));
      hisname.Form("hDLenAllPt%dLS",i);
      fDLenHistLS[index]=dynamic_cast<TH1F*>(fOutput->FindObject(hisname.Data()));
      hisname.Form("hSumd02AllPt%dLS",i);
      fSumd02HistLS[index]=dynamic_cast<TH1F*>(fOutput->FindObject(hisname.Data()));
      hisname.Form("hSigVertAllPt%dLS",i);
      fSigVertHistLS[index]=dynamic_cast<TH1F*>(fOutput->FindObject(hisname.Data()));
      hisname.Form("hPtMaxAllPt%dLS",i);
      fPtMaxHistLS[index]=dynamic_cast<TH1F*>(fOutput->FindObject(hisname.Data()));
      hisname.Form("hDCAAllPt%dLS",i);
      fDCAHistLS[index]=dynamic_cast<TH1F*>(fOutput->FindObject(hisname.Data()));
      
      hisname.Form("hLSPt%dTC",i);
      fMassHistLSTC[indexLS]=dynamic_cast<TH1F*>(fOutput->FindObject(hisname.Data()));
      
    } 
    
    index=GetSignalHistoIndex(i);    
    if(fDoLS)indexLS++;
    hisname.Form("hSigPt%d",i);
    fMassHist[index]=dynamic_cast<TH1F*>(fOutput->FindObject(hisname.Data()));
    hisname.Form("hCosPSigPt%d",i);
    fCosPHist[index]=dynamic_cast<TH1F*>(fOutput->FindObject(hisname.Data()));
    hisname.Form("hDLenSigPt%d",i);
    fDLenHist[index]=dynamic_cast<TH1F*>(fOutput->FindObject(hisname.Data()));
    hisname.Form("hSumd02SigPt%d",i);
    fSumd02Hist[index]=dynamic_cast<TH1F*>(fOutput->FindObject(hisname.Data()));
    hisname.Form("hSigVertSigPt%d",i);
    fSigVertHist[index]=dynamic_cast<TH1F*>(fOutput->FindObject(hisname.Data()));
    hisname.Form("hPtMaxSigPt%d",i);
    fPtMaxHist[index]=dynamic_cast<TH1F*>(fOutput->FindObject(hisname.Data()));
    hisname.Form("hDCASigPt%d",i);
    fDCAHist[index]=dynamic_cast<TH1F*>(fOutput->FindObject(hisname.Data()));
    
    hisname.Form("hSigPt%dTC",i);
    fMassHistTC[index]=dynamic_cast<TH1F*>(fOutput->FindObject(hisname.Data()));
    if(fDoLS){
      hisname.Form("hLSPt%dLCnw",i);
      fMassHistLS[indexLS]=dynamic_cast<TH1F*>(fOutput->FindObject(hisname.Data()));
      hisname.Form("hCosPSigPt%dLS",i);
      fCosPHistLS[index]=dynamic_cast<TH1F*>(fOutput->FindObject(hisname.Data()));
      hisname.Form("hDLenSigPt%dLS",i);
      fDLenHistLS[index]=dynamic_cast<TH1F*>(fOutput->FindObject(hisname.Data()));
      hisname.Form("hSumd02SigPt%dLS",i);
      fSumd02HistLS[index]=dynamic_cast<TH1F*>(fOutput->FindObject(hisname.Data()));
      hisname.Form("hSigVertSigPt%dLS",i);
      fSigVertHistLS[index]=dynamic_cast<TH1F*>(fOutput->FindObject(hisname.Data()));
      hisname.Form("hPtMaxSigPt%dLS",i);
      fPtMaxHistLS[index]=dynamic_cast<TH1F*>(fOutput->FindObject(hisname.Data()));
      hisname.Form("hDCASigPt%dLS",i);
      fDCAHistLS[index]=dynamic_cast<TH1F*>(fOutput->FindObject(hisname.Data()));

      hisname.Form("hLSPt%dTCnw",i);
      fMassHistLSTC[indexLS]=dynamic_cast<TH1F*>(fOutput->FindObject(hisname.Data()));
         }
    
    index=GetBackgroundHistoIndex(i); 
    if(fDoLS)indexLS++;
    hisname.Form("hBkgPt%d",i);
    fMassHist[index]=dynamic_cast<TH1F*>(fOutput->FindObject(hisname.Data()));
    hisname.Form("hCosPBkgPt%d",i);
    fCosPHist[index]=dynamic_cast<TH1F*>(fOutput->FindObject(hisname.Data()));
    hisname.Form("hDLenBkgPt%d",i);
    fDLenHist[index]=dynamic_cast<TH1F*>(fOutput->FindObject(hisname.Data()));
    hisname.Form("hSumd02BkgPt%d",i);
    fSumd02Hist[index]=dynamic_cast<TH1F*>(fOutput->FindObject(hisname.Data()));
    hisname.Form("hSigVertBkgPt%d",i);
    fSigVertHist[index]=dynamic_cast<TH1F*>(fOutput->FindObject(hisname.Data()));
    hisname.Form("hPtMaxBkgPt%d",i);
    fPtMaxHist[index]=dynamic_cast<TH1F*>(fOutput->FindObject(hisname.Data()));
    hisname.Form("hDCABkgPt%d",i);
    fDCAHist[index]=dynamic_cast<TH1F*>(fOutput->FindObject(hisname.Data()));
    hisname.Form("hBkgPt%dTC",i);
    fMassHistTC[index]=dynamic_cast<TH1F*>(fOutput->FindObject(hisname.Data()));
    if(fDoLS){
      hisname.Form("hLSPt%dLCntrip",i);
      fMassHistLS[indexLS]=dynamic_cast<TH1F*>(fOutput->FindObject(hisname.Data()));
 

      hisname.Form("hCosPBkgPt%dLS",i);
      fCosPHistLS[index]=dynamic_cast<TH1F*>(fOutput->FindObject(hisname.Data()));
      hisname.Form("hDLenBkgPt%dLS",i);
      fDLenHistLS[index]=dynamic_cast<TH1F*>(fOutput->FindObject(hisname.Data()));
      hisname.Form("hSumd02BkgPt%dLS",i);
      fSumd02HistLS[index]=dynamic_cast<TH1F*>(fOutput->FindObject(hisname.Data()));
      hisname.Form("hSigVertBkgPt%dLS",i);
      fSigVertHistLS[index]=dynamic_cast<TH1F*>(fOutput->FindObject(hisname.Data()));
      hisname.Form("hPtMaxBkgPt%dLS",i);
      fPtMaxHistLS[index]=dynamic_cast<TH1F*>(fOutput->FindObject(hisname.Data()));
      hisname.Form("hDCABkgPt%dLS",i);
      fDCAHistLS[index]=dynamic_cast<TH1F*>(fOutput->FindObject(hisname.Data()));
      
      hisname.Form("hLSPt%dTCntrip",i);
      fMassHistLSTC[indexLS]=dynamic_cast<TH1F*>(fOutput->FindObject(hisname.Data()));
      
      
      indexLS++;
      hisname.Form("hLSPt%dLCntripsinglecut",i);
      fMassHistLS[indexLS]=dynamic_cast<TH1F*>(fOutput->FindObject(hisname.Data()));
      
      hisname.Form("hLSPt%dTCntripsinglecut",i);
      fMassHistLSTC[indexLS]=dynamic_cast<TH1F*>(fOutput->FindObject(hisname.Data()));
      
      
      indexLS++;
      hisname.Form("hLSPt%dLCspc",i);
      fMassHistLS[indexLS]=dynamic_cast<TH1F*>(fOutput->FindObject(hisname.Data()));
      
      hisname.Form("hLSPt%dTCspc",i);
      fMassHistLSTC[indexLS]=dynamic_cast<TH1F*>(fOutput->FindObject(hisname.Data()));
    }
 
 }

  if(fFillNtuple){
    fNtupleDplus = dynamic_cast<TNtuple*>(GetOutputData(3));
  }

  TCanvas *c1=new TCanvas("c1","D+ invariant mass distribution",500,500);
  c1->cd();
  TH1F *hMassPt=(TH1F*)fOutput->FindObject("hMassPt3TC");
  hMassPt->SetLineColor(kBlue);
  hMassPt->SetXTitle("M[GeV/c^{2}]"); 
  hMassPt->Draw();
 
 return;
}
