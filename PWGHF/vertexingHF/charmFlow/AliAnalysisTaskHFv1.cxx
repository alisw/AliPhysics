/**************************************************************************
 * Copyright(c) 1998-2010, ALICE Experiment at CERN, All rights reserved. *
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

/////////////////////////////////////////////////////////////
//
// AliAnalysisTaskHFv1 directed flow of D mesons with scalar
// product method (modified from AliAnalysisTaskHFv1)
// Authors: Andrea Dubla, Jacopo Margutti
//
// AliAnalysisTaskHFv1 gives the needed tools for the D
// mesons vn analysis with event plane method
// Authors: Chiara Bianchin, Robert Grajcarek, Giacomo Ortona, 
//          Carlos Perez Lara, Francesco Prino, Anastasia Barbano,
//          Fabrizio Grosa, Andrea Festanti
// 
/////////////////////////////////////////////////////////////

/* $Id$ */

#include <Riostream.h>
#include <TClonesArray.h>
#include <TCanvas.h>
#include <TList.h>
#include <TFile.h>
#include <TString.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TGraphErrors.h>
#include <TGraph.h>
#include <TDatabasePDG.h>
#include <TRandom3.h>
#include <TVector2.h>
#include <TArrayF.h>
#include <TAxis.h>

#include <AliLog.h>
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
#include "AliAODRecoDecayHF.h"
#include "AliAODRecoDecayHF2Prong.h"
#include "AliAODRecoDecayHF4Prong.h"
#include "AliAODRecoCascadeHF.h"

#include "AliAnalysisTaskSE.h"
#include "AliRDHFCutsDplustoKpipi.h"
#include "AliRDHFCutsD0toKpipipi.h"
#include "AliRDHFCutsDstoKKpi.h"
#include "AliRDHFCutsDStartoKpipi.h"
#include "AliRDHFCutsD0toKpi.h"
#include "AliRDHFCutsLctopKpi.h"

#include "AliEventplane.h"
#include "AliFlowTrack.h"
#include "AliFlowVector.h"
#include "AliFlowTrackCuts.h"
#include "AliFlowEvent.h"

#include "AliAnalysisVertexingHF.h"
#include "AliEventplane.h"
#include "AliFlowTrack.h"
#include "AliFlowVector.h"
#include "AliFlowTrackCuts.h"
#include "AliFlowEvent.h"
#include "AliMultSelection.h"
#include "AliQnCorrectionsManager.h"
#include "AliAnalysisTaskFlowVectorCorrections.h"
#include "AliAnalysisTaskZDCEP.h"

#include "AliAnalysisTaskHFv1.h"

ClassImp(AliAnalysisTaskHFv1)


//________________________________________________________________________
AliAnalysisTaskHFv1::AliAnalysisTaskHFv1():
AliAnalysisTaskSE(),
  fhEventsInfo(0),
  fOutput(0),
  fRDCuts(0),
  fLowmasslimit(1.669),
  fUpmasslimit(2.069),
  fLowEtaLimit(-0.8),
  fUpEtaLimit(0.8),
  fNPtBins(1),
  fNEtaBins(4),
  fNMassBins(200),
  fReadMC(kFALSE),
  fUseAfterBurner(kFALSE),
  fDecChannel(0),
  fAfterBurner(0),
  fEventPlanesComp(10),
  fHarmonic(2),
  fMinCentr(20),
  fMaxCentr(80),
  fEtaGap(kFALSE),
  fEvPlaneDet(kFullV0),
  fSubEvDetA(kPosTPC),
  fSubEvDetB(kNegTPC),
  fCentBinSizePerMil(25),
  fAODProtection(1),
  fUseNewQnCorrFw(kTRUE),
  fCentrBinName(""),
  fFlowMethod(kEP),
  fNormMethod("QoverQlength"),
  fHistMassPtPhiq2Centr(0x0),
  fq2Meth(kq2TPC),
  fSeparateD0D0bar(kFALSE),
  fOnTheFlyTPCEP(kFALSE),
  fEtaGapInTPCHalves(-1.),
  fUsePtWeights(kFALSE),
  fOnTheFlyTPCq2(kFALSE),
  fFractionOfTracksForTPCq2(1.1),
  fq2SmearingHisto(0x0),
  fq2Smearing(kFALSE),
  fq2SmearingAxis(1),
  fScalProdLimit(0.3)
{
  // Default constructor
  for(int i = 0; i < 3; i++) {
    fHistCentrality[i]     = 0x0;
    fHistEvPlaneQncorrTPC[i]   = 0x0;
    fHistEvPlaneQncorrVZERO[i] = 0x0;
  }
  
}

//________________________________________________________________________
AliAnalysisTaskHFv1::AliAnalysisTaskHFv1(const char *name,AliRDHFCuts *rdCuts,Int_t decaychannel):
  AliAnalysisTaskSE(name),
  fhEventsInfo(0),
  fOutput(0),
  fRDCuts(rdCuts),
  fLowmasslimit(0),
  fUpmasslimit(0),
  fLowEtaLimit(-0.8),
  fUpEtaLimit(0.8),
  fNPtBins(1),
  fNEtaBins(4),
  fNMassBins(200),
  fReadMC(kFALSE),
  fUseAfterBurner(kFALSE),
  fDecChannel(decaychannel),
  fAfterBurner(0),
  fEventPlanesComp(10),
  fHarmonic(2),
  fMinCentr(20),
  fMaxCentr(80),
  fEtaGap(kFALSE),
  fEvPlaneDet(kFullV0),
  fSubEvDetA(kPosTPC),
  fSubEvDetB(kNegTPC),
  fCentBinSizePerMil(25),
  fAODProtection(1),
  fUseNewQnCorrFw(kTRUE),
  fCentrBinName(""),
  fFlowMethod(kEP),
  fNormMethod("QoverQlength"),
  fHistMassPtPhiq2Centr(0x0),
  fq2Meth(kq2TPC),
  fSeparateD0D0bar(kFALSE),
  fOnTheFlyTPCEP(kFALSE),
  fEtaGapInTPCHalves(-1.),
  fUsePtWeights(kFALSE),
  fOnTheFlyTPCq2(kFALSE),
  fFractionOfTracksForTPCq2(1.1),
  fq2SmearingHisto(0x0),
  fq2Smearing(kFALSE),
  fq2SmearingAxis(1),
  fScalProdLimit(0.3)
{
  // standard constructor
  for(int i = 0; i < 3; i++) {
    fHistCentrality[i]     = 0x0;
    fHistEvPlaneQncorrTPC[i]   = 0x0;
    fHistEvPlaneQncorrVZERO[i] = 0x0;
  }

  Int_t pdg=421;
  switch(fDecChannel){
  case 0:
    pdg=411;
    break;
  case 1:
    pdg=421;
    break;
  case 2:
    pdg=413;
    break;
  case 3:
    pdg=431;
    break;
  case 4:
    pdg=421;
    break;
  case 5:
    pdg=4122;
    break;
  }
  fAfterBurner = new AliHFAfterBurner(fDecChannel);
  if(pdg==413) SetMassLimits((Float_t)0.1,(Float_t)0.2);
  else SetMassLimits((Float_t)0.2,pdg); //check range
  fNPtBins=fRDCuts->GetNPtBins();
  
  DefineInput(1,AliFlowEventSimple::Class());

  if(fDebug>1)fRDCuts->PrintAll();
  // Output slot #1 writes into a TH1F container
  DefineOutput(1,TH1F::Class());   //Info on the number of events etc.
  // Output slot #2 writes into a TList container
  DefineOutput(2,TList::Class());  //Main output
  // Output slot #3 writes into a AliRDHFCuts container (cuts)
  switch(fDecChannel){
  case 0:
    DefineOutput(3,AliRDHFCutsDplustoKpipi::Class());  //Cut object for Dplus
    break;
  case 1:
    DefineOutput(3,AliRDHFCutsD0toKpi::Class());  //Cut object for D0
    break;
  case 2:
    DefineOutput(3,AliRDHFCutsDStartoKpipi::Class());  //Cut object for D*
    break;
  case 3:
    DefineOutput(3,AliRDHFCutsDstoKKpi::Class());  //Cut object for Ds
    break;
  }

  fDetTPCConfName[0] = "TPC";
  fDetTPCConfName[1] = "TPCNegEta";
  fDetTPCConfName[2] = "TPCPosEta";
    
  fDetV0ConfName[0]  = "VZERO";
  fDetV0ConfName[1]  = "VZEROA";
  fDetV0ConfName[2]  = "VZEROC";
  
}

//________________________________________________________________________
AliAnalysisTaskHFv1::~AliAnalysisTaskHFv1()
{
  // Destructor
  if(fOutput && !fOutput->IsOwner()){
    for(Int_t i=0;i<3;i++){
      delete fHistCentrality[i];
      delete fHistEvPlaneQncorrTPC[i];
      delete fHistEvPlaneQncorrVZERO[i];
    }
  }
  delete fOutput;
  delete fhEventsInfo;
  delete fRDCuts;
  delete fAfterBurner;
  if(fq2SmearingHisto) {delete fq2SmearingHisto;}
}
//_________________________________________________________________
void  AliAnalysisTaskHFv1::SetMassLimits(Float_t range, Int_t pdg){
  // Set limits for mass spectra plots
  Float_t mass=0;
  Int_t abspdg=TMath::Abs(pdg);
  mass=TDatabasePDG::Instance()->GetParticle(abspdg)->Mass();
  fUpmasslimit = mass+range;
  fLowmasslimit = mass-range;
}
//_________________________________________________________________
void  AliAnalysisTaskHFv1::SetMassLimits(Float_t lowlimit, Float_t uplimit){
  // Set limits for mass spectra plots
  if(uplimit>lowlimit)
    {
      fUpmasslimit = uplimit;
      fLowmasslimit = lowlimit;
    }
}
//________________________________________________________________________
void AliAnalysisTaskHFv1::LocalInit()
{
  // Initialization

  if(fDebug > 1) printf("AnalysisTaskSEHFvn::Init() \n");

  fRDCuts->SetMinCentrality(fMinCentr);
  fRDCuts->SetMaxCentrality(fMaxCentr);

  switch(fDecChannel){
  case 0:
    {
      AliRDHFCutsDplustoKpipi* copycut=new AliRDHFCutsDplustoKpipi(*(static_cast<AliRDHFCutsDplustoKpipi*>(fRDCuts)));
      // Post the data
      PostData(3,copycut);
    }
    break;
  case 1:
    {
      AliRDHFCutsD0toKpi* copycut=new AliRDHFCutsD0toKpi(*(static_cast<AliRDHFCutsD0toKpi*>(fRDCuts)));
      // Post the data
      PostData(3,copycut);
    }
    break;
  case 2:
    {
      AliRDHFCutsDStartoKpipi* copycut=new AliRDHFCutsDStartoKpipi(*(static_cast<AliRDHFCutsDStartoKpipi*>(fRDCuts)));
      // Post the data
      PostData(3,copycut);
    }
    break;
  case 3:
    {
      AliRDHFCutsDstoKKpi* copycut=new AliRDHFCutsDstoKKpi(*(static_cast<AliRDHFCutsDstoKKpi*>(fRDCuts)));
      // Post the data
      PostData(3,copycut);
    }
    break;
  default:
    return;
  }
  return;
}
//________________________________________________________________________
void AliAnalysisTaskHFv1::UserCreateOutputObjects()
{
  // Create the output container
 
  if(fDebug > 1) printf("AnalysisTaskSEHFvn::UserCreateOutputObjects() \n");

  fhEventsInfo = new TH1F(GetOutputSlot(1)->GetContainer()->GetName(), "Number of AODs scanned",15,-0.5,14.5);
  fhEventsInfo->GetXaxis()->SetBinLabel(1,"nEventsRead");
  fhEventsInfo->GetXaxis()->SetBinLabel(2,"nEvents Matched dAOD");
  fhEventsInfo->GetXaxis()->SetBinLabel(3,"nEvents Mismatched dAOD");
  fhEventsInfo->GetXaxis()->SetBinLabel(4,"nEventsAnal");
  fhEventsInfo->GetXaxis()->SetBinLabel(5,"n. passing IsEvSelected");
  fhEventsInfo->GetXaxis()->SetBinLabel(6,"n. rejected due to trigger");
  fhEventsInfo->GetXaxis()->SetBinLabel(7,"n. rejected due to not reco vertex");
  fhEventsInfo->GetXaxis()->SetBinLabel(8,"n. rejected for contr vertex");
  fhEventsInfo->GetXaxis()->SetBinLabel(9,"n. rejected for vertex out of accept");
  fhEventsInfo->GetXaxis()->SetBinLabel(10,"n. rejected for pileup events");
  fhEventsInfo->GetXaxis()->SetBinLabel(11,Form("no. of out %.0f-%.0f%s centrality events",fRDCuts->GetMinCentrality(),fRDCuts->GetMaxCentrality(),"%"));
  fhEventsInfo->GetXaxis()->SetBinLabel(12,"non valid TPC EP");
  fhEventsInfo->GetXaxis()->SetBinLabel(13,"bad event plane");
  fhEventsInfo->GetXaxis()->SetBinLabel(14,"no. of sel. candidates");
  fhEventsInfo->GetXaxis()->SetBinLabel(15,"no. cand. out of pt bounds");
    
  fhEventsInfo->GetXaxis()->SetNdivisions(1,kFALSE);

  // Several histograms are more conveniently managed in a TList
  fOutput = new TList();
  fOutput->SetOwner();
  fOutput->SetName("MainOutput");
 
  fHistCentrality[0]=new TH1F("hCentr","centrality",10000,0.,100.);
  fHistCentrality[1]=new TH1F("hCentr(selectedCent)","centrality(selectedCent)",10000,0.,100.);
  fHistCentrality[2]=new TH1F("hCentr(OutofCent)","centrality(OutofCent)",10000,0.,100.);
  for(Int_t i=0;i<3;i++){
    fHistCentrality[i]->Sumw2();
    fOutput->Add(fHistCentrality[i]);
  }
   
  int index=0;
  for(int iDet = 0; iDet < 3; iDet++) {
    fHistEvPlaneQncorrTPC[iDet]   = new TH1F(Form("hEvPlaneQncorr%s%s",fDetTPCConfName[iDet].Data(),fNormMethod.Data()),Form("hEvPlaneQncorr%s%s;#phi Ev Plane;Entries",fDetTPCConfName[iDet].Data(),fNormMethod.Data()),200,0.,TMath::Pi());
    fHistEvPlaneQncorrVZERO[iDet] = new TH1F(Form("hEvPlaneQncorr%s%s",fDetV0ConfName[iDet].Data(),fNormMethod.Data()),Form("hEvPlaneQncorr%s%s;#phi Ev Plane;Entries",fDetV0ConfName[iDet].Data(),fNormMethod.Data()),200,0.,TMath::Pi());
    fOutput->Add(fHistEvPlaneQncorrTPC[iDet]);
    fOutput->Add(fHistEvPlaneQncorrVZERO[iDet]);
  }
  if(fFlowMethod==kSP) {
    TH1F* hNormQA = new TH1F("hNormQA","hNormQA;QAx;QAy",100,0.,1);
    fOutput->Add(hNormQA);
    TH1F* hNormQB = new TH1F("hNormQB","hNormQB;QBx;QBy",100,0.,1);
    fOutput->Add(hNormQB);
  }

  for(Int_t icentr=fMinCentr*10+fCentBinSizePerMil;icentr<=fMaxCentr*10;icentr=icentr+fCentBinSizePerMil){
    TString centrname;centrname.Form("centr%d_%d",icentr-fCentBinSizePerMil,icentr);

    if(fFlowMethod!=kEvShape) {
      TH2F* hMPtCand=new TH2F(Form("hMPtCand%s",centrname.Data()),Form("Mass vs pt %s;p_{t} (GeV/c);M (GeV/c^{2})",centrname.Data()),500,0.,50.,fNMassBins,fLowmasslimit,fUpmasslimit);
      fOutput->Add(hMPtCand);//For <pt> calculation
    }
    if(fFlowMethod==kEP) {
      //Candidate distributions
      for(Int_t i=0;i<fNPtBins;i++){
        // Delta Phi histograms
        Double_t maxphi = TMath::Pi();
        if (fDecChannel == 2) maxphi = TMath::PiOver2();
        TH2F* hMdeltaphi=new TH2F(Form("hMdeltaphi_pt%d%s",i,centrname.Data()),Form("Mass vs #Delta#phi %s;#Delta#phi;M (GeV/c^{2})",centrname.Data()),96,0,maxphi,fNMassBins,fLowmasslimit,fUpmasslimit);
        fOutput->Add(hMdeltaphi);//for phi bins analysis
        TH2F* hMc2deltaphi=new TH2F(Form("hMc2deltaphi_pt%d%s",i,centrname.Data()),Form("Mass vs cos2#Delta#phi (p_{t} bin %d %s);cos2#Delta#phi;M (GeV/c^{2})",i,centrname.Data()),100,-1.,1.,fNMassBins,fLowmasslimit,fUpmasslimit);
        fOutput->Add(hMc2deltaphi);
        TH2F* hMs2deltaphi=new TH2F(Form("hMs2deltaphi_pt%d%s",i,centrname.Data()),Form("Mass vs sin2#Delta#phi (p_{t} bin %d %s);sin2#Delta#phi;M (GeV/c^{2})",i,centrname.Data()),100,-1.,1.,fNMassBins,fLowmasslimit,fUpmasslimit);
        fOutput->Add(hMs2deltaphi);
      
        // phi histograms (for systematics)
        TH2F* hCos2PhiMass=new TH2F(Form("hCos2phiMass_pt%d%s",i,centrname.Data()),Form("Mass vs cos(2#phi) %s;cos(2#phi);M (GeV/c^{2})",centrname.Data()),100,-1.,1.,fNMassBins,fLowmasslimit,fUpmasslimit);
        fOutput->Add(hCos2PhiMass);
        TH2F* hSin2PhiMass=new TH2F(Form("hSin2phiMass_pt%d%s",i,centrname.Data()),Form("Mass vs sin(2#phi) %s;sin(2#phi);M (GeV/c^{2})",centrname.Data()),100,-1.,1.,fNMassBins,fLowmasslimit,fUpmasslimit);
        fOutput->Add(hSin2PhiMass);

        // Histos using MC truth
        if (fReadMC){
          TH2F* hMc2deltaphiS=new TH2F(Form("hMc2deltaphiS_pt%d%s",i,centrname.Data()),Form("Mass vs cos2#Delta#phi (p_{t} bin %d %s);cos2#Delta#phi;M (GeV/c^{2})",i,centrname.Data()),100,-1.,1.,fNMassBins,fLowmasslimit,fUpmasslimit);
          fOutput->Add(hMc2deltaphiS);
          TH2F * hMdeltaphiS=new TH2F(Form("hMdeltaphiS_pt%d%s",i,centrname.Data()),Form("Mass vs #Delta#phi (p_{t} bin %d %s);#Delta#phi;M (GeV/c^{2})",i,centrname.Data()),96,0,2*TMath::Pi(),fNMassBins,fLowmasslimit,fUpmasslimit);
          fOutput->Add(hMdeltaphiS);
          TH2F* hMc2deltaphiB=new TH2F(Form("hMc2deltaphiB_pt%d%s",i,centrname.Data()),Form("Mass vs cos2#Delta#phi (p_{t} bin %d %s);cos2#Delta#phi;M (GeV/c^{2})",i,centrname.Data()),100,-1.,1.,fNMassBins,fLowmasslimit,fUpmasslimit);
          fOutput->Add(hMc2deltaphiB);
          TH2F * hMdeltaphiB=new TH2F(Form("hMdeltaphiB_pt%d%s",i,centrname.Data()),Form("Mass vs #Delta#phi (p_{t} bin %d %s);#Delta#phi;M (GeV/c^{2})",i,centrname.Data()),96,0,2*TMath::Pi(),fNMassBins,fLowmasslimit,fUpmasslimit);
          fOutput->Add(hMdeltaphiB);
          if((fDecChannel != AliAnalysisTaskHFv1::kDplustoKpipi) &&(fDecChannel != AliAnalysisTaskHFv1::kDstartoKpipi)){
            TH2F* hMc2deltaphiR=new TH2F(Form("hMc2deltaphiR_pt%d%s",i,centrname.Data()),Form("Mass vs cos2#Delta#phi (p_{t} bin %d %s);cos2#Delta#phi;M (GeV/c^{2})",i,centrname.Data()),100,-1.,1.,fNMassBins,fLowmasslimit,fUpmasslimit);
            fOutput->Add(hMc2deltaphiR);
            TH2F* hMdeltaphiR=new TH2F(Form("hMdeltaphiR_pt%d%s",i,centrname.Data()),Form("Mass vs #Delta#phi (p_{t} bin %d %s);#Delta#phi;M (GeV/c^{2})",i,centrname.Data()),96,0,2*TMath::Pi(),fNMassBins,fLowmasslimit,fUpmasslimit);
            fOutput->Add(hMdeltaphiR);
          }
        }
      }

      // Event Plane
      TH2F* hEvPlane=new TH2F(Form("hEvPlane%s",centrname.Data()),Form("VZERO/TPC Event plane angle %s;#phi Ev Plane (TPC);#phi Ev Plane (VZERO);Entries",centrname.Data()),200,0.,TMath::Pi(),200,0.,TMath::Pi());
      fOutput->Add(hEvPlane);
      TH2F* hEvPlaneRecomp=new TH2F(Form("hEvPlaneRecomp%s",centrname.Data()),Form("TPC Event plane angle before after recomputation %s;#phi Ev Plane (TPC);#phi Ev Plane (TPC recomp);Entries",centrname.Data()),200,0.,TMath::Pi(),200,0.,TMath::Pi());
      fOutput->Add(hEvPlaneRecomp);
      TH1F* hEvPlaneA=new TH1F(Form("hEvPlaneA%s",centrname.Data()),Form("Event plane angle %s;#phi Ev Plane;Entries",centrname.Data()),200,0.,TMath::Pi());
      fOutput->Add(hEvPlaneA);
      TH1F* hEvPlaneB=new TH1F(Form("hEvPlaneB%s",centrname.Data()),Form("Event plane angle %s;#phi Ev Plane;Entries",centrname.Data()),200,0.,TMath::Pi());
      fOutput->Add(hEvPlaneB);
      TH1F* hEvPlaneCand=new TH1F(Form("hEvPlaneCand%s",centrname.Data()),Form("Event plane angle - Event plane angle per candidate %s;#phi(Ev Plane Candidate);Entries",centrname.Data()),200,-TMath::Pi(),TMath::Pi());
      fOutput->Add(hEvPlaneCand);
      
      // histos for EP resolution
      TH1F* hEvPlaneReso=new TH1F(Form("hEvPlaneReso1%s",centrname.Data()),Form("Event plane angle Resolution %s;cos2(#psi_{A}-#psi_{B});Entries",centrname.Data()),220,-1.1,1.1);
      fOutput->Add(hEvPlaneReso);
      TH1F* hEvPlaneReso2=new TH1F(Form("hEvPlaneReso2%s",centrname.Data()),Form("Event plane angle Resolution %s;cos2(#psi_{A}-#psi_{C});Entries",centrname.Data()),220,-1.1,1.1);
      fOutput->Add(hEvPlaneReso2);
      TH1F* hEvPlaneReso3=new TH1F(Form("hEvPlaneReso3%s",centrname.Data()),Form("Event plane angle Resolution %s;cos2(#psi_{B}-#psi_{C});Entries",centrname.Data()),220,-1.1,1.1);
      fOutput->Add(hEvPlaneReso3);

      // histos for EPsystematics
      TH1F *hCos2EP=new TH1F(Form("hCos2EP%s",centrname.Data()),Form("cos(2PsiEP) %s;cos2(#psi_{EP});Entries",centrname.Data()),100,-1.,1.);
      fOutput->Add(hCos2EP);
      TH1F *hSin2EP=new TH1F(Form("hSin2EP%s",centrname.Data()),Form("sin(2PsiEP) %s;sin2(#psi_{EP});Entries",centrname.Data()),100,-1.,1.);
      fOutput->Add(hSin2EP);
    }
    //Scalar Product Histograms
    else if(fFlowMethod==kSP) {
      TString SPsuffix1="";
      TString SPsuffix2="";
      TString SPsuffix3="QAQC";
      if(fEvPlaneDet==kFullTPC || fEvPlaneDet==kPosTPC || fEvPlaneDet==kNegTPC) { //two sub-event formula
        SPsuffix1="uAQB";
        SPsuffix2="uBQA";
      }
      else { // three sub-event formula A,B -> sub-events for Q-vetors C -> sub-event for D-mesons
        SPsuffix1="uQA";
        SPsuffix2="uQC";
      }
      Int_t nBins[5] = {100,fNMassBins,50,4,3};
      Double_t xmin[5] = {-fScalProdLimit,fLowmasslimit,0.,-0.8,0.};
      Double_t xmax[5] = {fScalProdLimit,fUpmasslimit,50.,0.8,3.};
      THnSparseD* hMassScalProd1 = new THnSparseD(Form("hMassScalProd%s_%s",SPsuffix1.Data(),centrname.Data()),Form("Mass vs hScalProd%s ;%s;M (GeV/c^{2})",SPsuffix1.Data(),centrname.Data()),5,nBins,xmin,xmax);
      hMassScalProd1->Sumw2();
      fOutput->Add(hMassScalProd1);
      THnSparseD* hMassScalProd2 = new THnSparseD(Form("hMassScalProd%s_%s",SPsuffix2.Data(),centrname.Data()),Form("Mass vs hScalProd%s ;%s;M (GeV/c^{2})",SPsuffix2.Data(),centrname.Data()),5,nBins,xmin,xmax);
      hMassScalProd2->Sumw2();
      fOutput->Add(hMassScalProd2);
      TProfile* hScalProdDenom = new TProfile(Form("hScalProd%s",SPsuffix3.Data()),Form("hScalProd%s",SPsuffix3.Data()),(Int_t)1000./fCentBinSizePerMil,0.,100.);
      fOutput->Add(hScalProdDenom);
    }
    else if (fFlowMethod==kEvShape){
      TString q2axisname="q_{2}^{TPC}";
      if(fq2Meth==kq2PosTPC) {q2axisname="q_{2}^{pos TPC}";}
      else if(fq2Meth==kq2NegTPC) {q2axisname="q_{2}^{neg TPC}";}
      else if(fq2Meth==kq2VZERO) {q2axisname="q_{2}^{V0}";}
      else if(fq2Meth==kq2VZEROA) {q2axisname="q_{2}^{V0A}";}
      else if(fq2Meth==kq2VZEROC) {q2axisname="q_{2}^{V0C}";}
      
      //q2TPC correlations
      TH2F* hq2TPCPosEtaVsNegEta = new TH2F("hq2TPCPosEtaVsNegEta","q_{2}^{pos TPC} vs. q_{2}^{neg TPC};q_{2}^{neg TPC};q_{2}^{pos TPC}",500,0.,10.,500,0.,10.);
      TH2F* hq2TPCFullEtaVsNegEta = new TH2F("hq2TPCFullEtaVsNegEta","q_{2}^{full TPC} vs. q_{2}^{neg TPC};q_{2}^{neg TPC};q_{2}^{full TPC}",500,0.,10.,500,0.,10.);
      TH2F* hq2TPCFullEtaVsPosEta = new TH2F("hq2TPCFullEtaVsPosEta","q_{2}^{full TPC} vs. q_{2}^{pos TPC};q_{2}^{pos TPC};q_{2}^{full TPC}",500,0.,10.,500,0.,10.);
      fOutput->Add(hq2TPCPosEtaVsNegEta);
      fOutput->Add(hq2TPCFullEtaVsNegEta);
      fOutput->Add(hq2TPCFullEtaVsPosEta);
      
      // histos for EP resolution vs q2
      TH2F* hEvPlaneResoVsq2=new TH2F(Form("hEvPlaneReso1Vsq2%s",centrname.Data()),Form("Event plane angle Resolution vs. %s %s;cos2(#psi_{A}-#psi_{B});%s;Entries",q2axisname.Data(),centrname.Data(),q2axisname.Data()),220,-1.1,1.1,500,0,10.);
      fOutput->Add(hEvPlaneResoVsq2);
      TH2F* hEvPlaneReso2Vsq2=new TH2F(Form("hEvPlaneReso2Vsq2%s",centrname.Data()),Form("Event plane angle Resolution vs. %s %s;cos2(#psi_{A}-#psi_{C});%s;Entries",q2axisname.Data(),centrname.Data(),q2axisname.Data()),220,-1.1,1.1,500,0,10.);
      fOutput->Add(hEvPlaneReso2Vsq2);
      TH2F* hEvPlaneReso3Vsq2=new TH2F(Form("hEvPlaneReso3Vsq2%s",centrname.Data()),Form("Event plane angle Resolution vs. %s %s;cos2(#psi_{B}-#psi_{C});%s;Entries",q2axisname.Data(),centrname.Data(),q2axisname.Data()),220,-1.1,1.1,500,0.,10.);
      fOutput->Add(hEvPlaneReso3Vsq2);
    }
  }
  
  if(fFlowMethod==kEvShape) {CreateSparseForEvShapeAnalysis();}
  
  PostData(1,fhEventsInfo);
  PostData(2,fOutput);

  return;
}

//________________________________________________________________________
void AliAnalysisTaskHFv1::UserExec(Option_t */*option*/)
{
  // Execute analysis for current event:
  // heavy flavor candidates association to MC truth
  AliAODEvent *aod = dynamic_cast<AliAODEvent*> (InputEvent());
  
  fhEventsInfo->Fill(0);

  //   Protection against the mismatch of candidate TRefs:
  //   Check if AOD and corresponding deltaAOD files contain the same number of events.
  //   In case of discrepancy the event is rejected.
  if(fAODProtection>=0){
    //   Protection against different number of events in the AOD and deltaAOD
    //   In case of discrepancy the event is rejected.
    Int_t matchingAODdeltaAODlevel = AliRDHFCuts::CheckMatchingAODdeltaAODevents();
    if (matchingAODdeltaAODlevel<0 || (matchingAODdeltaAODlevel==0 && fAODProtection==1)) {
      // AOD/deltaAOD trees have different number of entries || TProcessID do not match while it was required
      fhEventsInfo->Fill(2);
      return;
    }
    fhEventsInfo->Fill(1);
  }
    
  if(fDebug>2) printf("Analysing decay %d\n",fDecChannel);
  // Post the data already here
  PostData(1,fhEventsInfo);
  PostData(2,fOutput);

  TClonesArray *arrayProng =0;
  Int_t absPdgMom=0;
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
   
      
      if(fDecChannel==0){
	absPdgMom=411;
	arrayProng=(TClonesArray*)aodFromExt->GetList()->FindObject("Charm3Prong");
      }
      if(fDecChannel==1){
	absPdgMom=421;
	arrayProng=(TClonesArray*)aodFromExt->GetList()->FindObject("D0toKpi");
      }
      if(fDecChannel==2){
	absPdgMom=413;
	arrayProng=(TClonesArray*)aodFromExt->GetList()->FindObject("Dstar");
      }
      if(fDecChannel==3){
	absPdgMom=431;
	arrayProng=(TClonesArray*)aodFromExt->GetList()->FindObject("Charm3Prong");
      }
    }
  } else if(aod){
    if(fDecChannel==0){
      absPdgMom=411;
      arrayProng=(TClonesArray*)aod->GetList()->FindObject("Charm3Prong");
    }
    if(fDecChannel==1){
      absPdgMom=421;
      arrayProng=(TClonesArray*)aod->GetList()->FindObject("D0toKpi");
    }
    if(fDecChannel==2){
      absPdgMom=413;
      arrayProng=(TClonesArray*)aod->GetList()->FindObject("Dstar");
    }
    if(fDecChannel==3){
      absPdgMom=431;
      arrayProng=(TClonesArray*)aod->GetList()->FindObject("Charm3Prong");
    }
  }

  if(!aod || !arrayProng) {
    AliError("AliAnalysisTaskHFv1::UserExec:Branch not found!\n");
    return;
  }
  
  // fix for temporary bug in ESDfilter 
  // the AODs with null vertex pointer didn't pass the PhysSel
  if(!aod->GetPrimaryVertex() || TMath::Abs(aod->GetMagneticField())<0.001) return;

  fhEventsInfo->Fill(3);

  TClonesArray *arrayMC=0;
  AliAODMCHeader *mcHeader=0;

  // load MC particles
  if(fReadMC){
    
    arrayMC =  (TClonesArray*)aod->GetList()->FindObject(AliAODMCParticle::StdBranchName());
    if(!arrayMC) {
      AliWarning("AliAnalysisTaskHFv1::UserExec:MC particles branch not found!\n");
      return;
    }
    
    // load MC header
    mcHeader =  (AliAODMCHeader*)aod->GetList()->FindObject(AliAODMCHeader::StdBranchName());
    if(!mcHeader) {
      AliError("AliAnalysisTaskHFv1::UserExec:MC header branch not found!\n");
      return;
    }
  }

  AliAODRecoDecayHF *d=0;
  
  Int_t nCand = arrayProng->GetEntriesFast();
  if(fDebug>2) printf("Number of D2H: %d\n",nCand);
  
  Bool_t isEvSel=fRDCuts->IsEventSelected(aod);
  Float_t evCentr=fRDCuts->GetCentrality(aod);
  fHistCentrality[0]->Fill(evCentr);
  if(!isEvSel){
    if(fRDCuts->IsEventRejectedDueToTrigger())fhEventsInfo->Fill(5);
    if(fRDCuts->IsEventRejectedDueToNotRecoVertex())fhEventsInfo->Fill(6);
    if(fRDCuts->IsEventRejectedDueToVertexContributors())fhEventsInfo->Fill(7);
    if(fRDCuts->IsEventRejectedDueToZVertexOutsideFiducialRegion())fhEventsInfo->Fill(8);
    if(fRDCuts->IsEventRejectedDueToPileup())fhEventsInfo->Fill(9);
    if(fRDCuts->IsEventRejectedDueToCentrality()){
      fhEventsInfo->Fill(10);
      fHistCentrality[2]->Fill(evCentr);
    }
    return;
  }
  
  fhEventsInfo->Fill(4);
  fHistCentrality[1]->Fill(evCentr);
  
  Double_t eventplaneqncorrTPC[3];
  Double_t eventplaneqncorrVZERO[3];
  TList *qnlist = 0x0;
  Double_t QA[2]={-2.-2};
  Double_t QB[2]={-2.-2};
 
  if(fUseNewQnCorrFw) {
    /////////////////////////////////////////////////////////////
    //////////////          GET Qn vectors         //////////////
    /////////////////////////////////////////////////////////////
    
    AliQnCorrectionsManager *flowQnVectorMgr;
    AliAnalysisTaskFlowVectorCorrections *flowQnVectorTask = dynamic_cast<AliAnalysisTaskFlowVectorCorrections*>(AliAnalysisManager::GetAnalysisManager()->GetTask("FlowQnVectorCorrections"));
    if (flowQnVectorTask != NULL) {
      flowQnVectorMgr = flowQnVectorTask->GetAliQnCorrectionsManager();
    }
    else {
      AliWarning("This task needs the Flow Qn vector corrections framework and it is not present. Aborting!!!\n");
      return;
    }
    qnlist = flowQnVectorMgr->GetQnVectorList();
    if (!qnlist) {
      return;
    }
    const AliQnCorrectionsQnVector* qnVectTPC[3];
    const AliQnCorrectionsQnVector* qnVectV0[3];
      
    for(int iDet = 0; iDet < 3; iDet++) {
      qnVectTPC[iDet]  = GetQnVectorFromList(qnlist, Form("%s%s",fDetTPCConfName[iDet].Data(),fNormMethod.Data()), "latest", "plain");
      qnVectV0[iDet]   = GetQnVectorFromList(qnlist, Form("%s%s",fDetV0ConfName[iDet].Data(),fNormMethod.Data()), "latest", "raw");
      if(!qnVectTPC[iDet] || !qnVectV0[iDet]) return;
      eventplaneqncorrTPC[iDet]   = qnVectTPC[iDet]->EventPlane(fHarmonic);
      eventplaneqncorrVZERO[iDet] = qnVectV0[iDet]->EventPlane(fHarmonic);
      if(eventplaneqncorrTPC[iDet]<0.)   eventplaneqncorrTPC[iDet]    += 2.*(TMath::Pi())/fHarmonic;
      if(eventplaneqncorrVZERO[iDet]<0.) eventplaneqncorrVZERO[iDet]  += 2.*(TMath::Pi())/fHarmonic;
      fHistEvPlaneQncorrTPC[iDet]->Fill(eventplaneqncorrTPC[iDet]);
      fHistEvPlaneQncorrVZERO[iDet]->Fill(eventplaneqncorrVZERO[iDet]);
    }
    
    //get Qx and Qy for SP analysis
    if(fEvPlaneDet==kFullTPC || fEvPlaneDet==kPosTPC || fEvPlaneDet==kNegTPC) {
      QA[0]=qnVectTPC[1]->Qx(fHarmonic);//NegTPC
      QA[1]=qnVectTPC[1]->Qy(fHarmonic);//NegTPC
      QB[0]=qnVectTPC[2]->Qx(fHarmonic);//PosTPC
      QB[1]=qnVectTPC[2]->Qy(fHarmonic);//PosTPC
    }
    else {
      QA[0]=qnVectV0[1]->Qx(fHarmonic);//V0A
      QA[1]=qnVectV0[1]->Qy(fHarmonic);//V0A
      QB[0]=qnVectV0[2]->Qx(fHarmonic);//V0C
      QB[1]=qnVectV0[2]->Qy(fHarmonic);//V0C
    }
    if(fFlowMethod==kSP) {
      ((TH1F*)fOutput->FindObject("hNormQA"))->Fill(TMath::Sqrt(QA[0]*QA[0]+QA[1]*QA[1]));
      ((TH1F*)fOutput->FindObject("hNormQB"))->Fill(TMath::Sqrt(QB[0]*QB[0]+QB[1]*QB[1]));
    }
  }
  else {
    //get Qx and Qy for SP analysis
    if(fEvPlaneDet==kZDCA || fEvPlaneDet==kZDCC) {
      AliAnalysisTaskZDCEP *fZDCEPTask = dynamic_cast<AliAnalysisTaskZDCEP*>(AliAnalysisManager::GetAnalysisManager()->GetTask("AnalysisTaskZDCEP"));
      if (fZDCEPTask != NULL) {
        // get ZDC Q-vectors
        AliFlowEvent* anEvent = dynamic_cast<AliFlowEvent*>(GetInputData(1));
        if(anEvent) {
          // Get Q vectors for the subevents
          AliFlowVector vQarray[2];
          anEvent->GetZDC2Qsub(vQarray);
          QA[0] = vQarray[0].X();
          QA[1] = vQarray[0].Y();
          QB[0] = vQarray[1].X();
          QB[1] = vQarray[1].Y();
        } else {
          printf("flowevent not found!!\n");
        }
      }
      else {
        AliWarning("This task needs AliAnalysisTaskZDCEP and it is not present. Aborting!!!\n");
        return;
      }
    }
    if(fFlowMethod==kSP) {
      ((TH1F*)fOutput->FindObject("hNormQA"))->Fill(TMath::Sqrt(QA[0]*QA[0]+QA[1]*QA[1]));
      ((TH1F*)fOutput->FindObject("hNormQB"))->Fill(TMath::Sqrt(QB[0]*QB[0]+QB[1]*QB[1]));
    }
  }
  
  
  //get q2 for event shape anaylsis
  Double_t q2=-1;
  Double_t q2PosTPC=-1;
  Double_t q2NegTPC=-1;
  Double_t q2FullTPC=-1;
  if(fFlowMethod==kEvShape) {
    if(fOnTheFlyTPCq2){
      q2=ComputeTPCq2(aod,q2FullTPC,q2PosTPC,q2NegTPC);
    }else{
      q2 = Getq2(qnlist,fq2Meth);
      q2PosTPC = Getq2(qnlist,kq2PosTPC);
      q2NegTPC = Getq2(qnlist,kq2NegTPC);
      q2FullTPC = Getq2(qnlist,kq2TPC);
    }
    if(q2<0 || q2PosTPC<0 || q2NegTPC<0 || q2FullTPC<0) return;
    ((TH1F*)fOutput->FindObject("hq2TPCPosEtaVsNegEta"))->Fill(q2NegTPC,q2PosTPC);
    ((TH1F*)fOutput->FindObject("hq2TPCFullEtaVsNegEta"))->Fill(q2NegTPC,q2FullTPC);
    ((TH1F*)fOutput->FindObject("hq2TPCFullEtaVsPosEta"))->Fill(q2PosTPC,q2FullTPC);
  }
  if(fq2Smearing && fq2SmearingHisto) {
    TAxis* ax=0x0;
    if(fq2SmearingAxis==1) {ax=(TAxis*)fq2SmearingHisto->GetYaxis();}
    else {ax=(TAxis*)fq2SmearingHisto->GetXaxis();}
    Int_t bin = ax->FindBin(q2);
    TH1F* hq2Slice = 0x0;
    if(fq2SmearingAxis==1) {hq2Slice = (TH1F*)fq2SmearingHisto->ProjectionX("hq2Slice",bin,bin);}
    else {hq2Slice = (TH1F*)fq2SmearingHisto->ProjectionY("hq2Slice",bin,bin);}
    if(hq2Slice->GetEntries()>10) {q2 = hq2Slice->GetRandom();}
    delete hq2Slice;
  }
  
  AliEventplane *pl=aod->GetEventplane();
  if(!pl){
    Printf("AliAnalysisTaskHFv1::UserExec:no eventplane! v2 analysis without eventplane not possible!\n");
    fhEventsInfo->Fill(11);
    return;
  }

  //determine centrality bin
  Float_t centr=fRDCuts->GetCentrality(aod);
  Float_t centrPerMil=centr*10.;
  Int_t icentr=-1;
  for(Int_t ic=fMinCentr*10+fCentBinSizePerMil;ic<=fMaxCentr*10;ic=ic+fCentBinSizePerMil){
    if(ic>centrPerMil){
      icentr=ic;
      break;
    }
  }
  if(icentr==-1) return;

  fCentrBinName=Form("centr%d_%d",icentr-fCentBinSizePerMil,icentr);
  Double_t eventplane=0;

  if(fReadMC){
    TRandom3 *g = new TRandom3(0);
    eventplane=g->Rndm()*TMath::Pi();
    delete g;g=0x0;
    if(fUseAfterBurner)fAfterBurner->SetEventPlane((Double_t)eventplane);
  }else{
    if(fFlowMethod!=kSP) {
      eventplane=GetEventPlane(aod,pl,eventplaneqncorrTPC,eventplaneqncorrVZERO,q2);
      if(eventplane<-999){
        Printf("Bad event plane calculation\n");
        fhEventsInfo->Fill(12);
        return;
      }
    }
    else {
      //fill histograms for denominator in SP formula
      Double_t scalprod=QA[0]*QB[0]+QA[1]*QB[1];
      ((TProfile*)fOutput->FindObject("hScalProdQAQC"))->Fill(centr,scalprod);
    }
  }
  if(fDebug>2)printf("Loop on D candidates\n");
  
  AliAnalysisVertexingHF *vHF=new AliAnalysisVertexingHF();

  //Loop on D candidates
  for (Int_t iCand = 0; iCand < nCand; iCand++) {
      
    d=(AliAODRecoDecayHF*)arrayProng->UncheckedAt(iCand);
    Bool_t isSelBit=kTRUE;
    if(fDecChannel==0) isSelBit=d->HasSelectionBit(AliRDHFCuts::kDplusCuts);
    if(fDecChannel==1) isSelBit=d->HasSelectionBit(AliRDHFCuts::kD0toKpiCuts);
    if(fDecChannel==2) isSelBit=kTRUE;
    if(fDecChannel==3) isSelBit=d->HasSelectionBit(AliRDHFCuts::kDsCuts);
    if(!isSelBit)continue;
      
    if(fDecChannel==0 || fDecChannel==3) {
      if(!vHF->FillRecoCand(aod,(AliAODRecoDecayHF3Prong*)d))continue;
    }
    else if(fDecChannel == 1) {
      if(!vHF->FillRecoCand(aod,(AliAODRecoDecayHF2Prong*)d))continue;
    }
    else if(fDecChannel == 2) {
      if(!vHF->FillRecoCasc(aod,((AliAODRecoCascadeHF*)d),kTRUE))continue; 
    }
    
    Int_t ptbin=fRDCuts->PtBin(d->Pt());
    if(ptbin<0) {
      fhEventsInfo->Fill(14);
      continue;
    }
    // hard-coded: force ptbin=0
    ptbin=0;
    
    Bool_t isFidAcc = fRDCuts->IsInFiducialAcceptance(d->Pt(),d->Y(absPdgMom));
    if(!isFidAcc)continue;    
    Int_t isSelected= fRDCuts->IsSelected(d,AliRDHFCuts::kAll,aod);
    if(!isSelected)continue;
    if(fDecChannel==3) {
      Int_t isDsPhiKKpi = isSelected&4;
      Int_t isDsPhipiKK = isSelected&8;
      if(!isDsPhiKKpi & !isDsPhipiKK) continue;
    }

    fhEventsInfo->Fill(13); // candidate selected
    if(fDebug>3) printf("+++++++Is Selected\n");
      
    Float_t* invMass=0x0;
    Int_t nmasses;
    CalculateInvMasses(d,invMass,nmasses);
    
    if(fEvPlaneDet==kFullTPC || fEvPlaneDet==kPosTPC || fEvPlaneDet==kNegTPC){
      Float_t eventplaneOld=eventplane;
      if(fUseNewQnCorrFw) {
        Float_t eventplaneOld=eventplane;
        eventplane = GetEventPlaneForCandidateNewQnFw(d,qnlist); // remove autocorrelations for new Qn fw
        if(eventplane==-999) eventplane = eventplaneOld;
      }
      else eventplane = GetEventPlaneForCandidate(d,pl); // remove autocorrelations
      if(fFlowMethod==kEP) ((TH1F*)fOutput->FindObject(Form("hEvPlaneCand%s",fCentrBinName.Data())))->Fill(eventplaneOld-eventplane);
    }
    
    Double_t phi=d->Phi();
    Double_t eta=d->Eta();
    Double_t pt=d->Pt();
    if(fReadMC&&fUseAfterBurner)phi=fAfterBurner->GetNewAngle(d,arrayMC);
    Float_t deltaphi  = GetPhiInRange(phi-eventplane);
    
    //fill the histograms with the appropriate method
    if(fFlowMethod!=kEvShape) {
      if(fDecChannel==0)FillDplus(d,arrayMC,ptbin,deltaphi,invMass,isSelected,icentr,phi,eta,pt,QA,QB);
      else if(fDecChannel==1)FillD02p(d,arrayMC,ptbin,deltaphi,invMass,isSelected,icentr,phi,eta,pt,QA,QB);
      else if(fDecChannel==2)FillDstar(d,arrayMC,ptbin,deltaphi,invMass,isSelected,icentr,phi,eta,pt,QA,QB);
      else if(fDecChannel==3)FillDs(d,arrayMC,ptbin,deltaphi,invMass,isSelected,icentr,phi,eta,pt,QA,QB);
    }
    else {      
      if(!fUseNewQnCorrFw){
        printf("ERROR: event shape analysis implemented only with new Qn correction framework\n");
        return;
      }
      
      //fill the THnSparseF for event-shape engineering
      if((fDecChannel==0 || fDecChannel==2) && isSelected) {
        Double_t sparsearray[5] = {invMass[0],d->Pt(),deltaphi,q2,centr};
        fHistMassPtPhiq2Centr->Fill(sparsearray);
      }
      else if(fDecChannel==1) {
        if(fSeparateD0D0bar) {
          if(isSelected==1 || isSelected==3) {
            Double_t sparsearray1[6] = {invMass[0],d->Pt(),deltaphi,q2,centr,(Double_t)isSelected};
            fHistMassPtPhiq2Centr->Fill(sparsearray1);
          }
          if(isSelected>=2) {
            Double_t sparsearray2[6] = {invMass[1],d->Pt(),deltaphi,q2,centr,(Double_t)isSelected};
            fHistMassPtPhiq2Centr->Fill(sparsearray2);
          }
        }
        else {
          if(isSelected==1 || isSelected==3) {
            Double_t sparsearray1[5] = {invMass[0],d->Pt(),deltaphi,q2,centr};
            fHistMassPtPhiq2Centr->Fill(sparsearray1);
          }
          if(isSelected>=2) {
            Double_t sparsearray2[5] = {invMass[1],d->Pt(),deltaphi,q2,centr};
            fHistMassPtPhiq2Centr->Fill(sparsearray2);
          }
        }
      }
      else if(fDecChannel==3) {
        if(isSelected==1 || isSelected==3) {
          Double_t sparsearray1[5] = {invMass[0],d->Pt(),deltaphi,q2,centr};
          fHistMassPtPhiq2Centr->Fill(sparsearray1);
        }
        if(isSelected>=2) {
          Double_t sparsearray2[5] = {invMass[1],d->Pt(),deltaphi,q2,centr};
          fHistMassPtPhiq2Centr->Fill(sparsearray2);
        }
      }
    }
    delete [] invMass;
  }
 
  delete vHF;
  PostData(1,fhEventsInfo);
  PostData(2,fOutput);

  return;
}

void AliAnalysisTaskHFv1::CreateSparseForEvShapeAnalysis() {
  /// Sparse for v2 analysis as a function of q2
 
  Int_t nptbins=200;
  Double_t ptmin=0.;
  Double_t ptmax=50.;
  
  Int_t nphibins=96;
  Double_t minphi=0.;
  Double_t maxphi = TMath::Pi();
  
  Int_t nq2bins=500;
  Double_t minq2=0.;
  Double_t maxq2=10.;

  Int_t ncentbins=(fMaxCentr-fMinCentr)/(fCentBinSizePerMil/10);
  
  TString massaxisname;
  if(fDecChannel==0) massaxisname = "M_{K#pi#pi} (GeV/c^{2})";
  else if(fDecChannel==1) massaxisname = "M_{K#pi} (GeV/c^{2})";
  else if(fDecChannel==2) massaxisname = "M_{K#pi#pi}-M_{K#pi} (GeV/c^{2})";
  else if(fDecChannel==3) massaxisname = "M_{KK#pi} (GeV/c^{2})";
  
  Int_t naxes=5;
  
  if(fSeparateD0D0bar && fDecChannel==1) {
    Int_t npartantipartbins=3;
    Double_t minpartantipart=0.5;
    Double_t maxpartantipart=3.5;
    
    Int_t nbins[6]={fNMassBins,nptbins,nphibins,nq2bins,ncentbins,npartantipartbins};
    Double_t xmin[6]={fLowmasslimit,ptmin,minphi,minq2,(Double_t)fMinCentr,(Double_t)minpartantipart};
    Double_t xmax[6]={fUpmasslimit,ptmax,maxphi,maxq2,(Double_t)fMaxCentr,(Double_t)maxpartantipart};
    
    naxes=6;
    
    fHistMassPtPhiq2Centr=new THnSparseF("hMassPtPhiq2Centr","Mass vs. pt vs. #Delta#phi vs. q_{2} vs. centr vs. D0-D0bar",naxes,nbins,xmin,xmax);
  }
  else {
    Int_t nbins[5]={fNMassBins,nptbins,nphibins,nq2bins,ncentbins};
    Double_t xmin[5]={fLowmasslimit,ptmin,minphi,minq2,(Double_t)fMinCentr};
    Double_t xmax[5]={fUpmasslimit,ptmax,maxphi,maxq2,(Double_t)fMaxCentr};
    
    fHistMassPtPhiq2Centr=new THnSparseF("hMassPtPhiq2Centr","Mass vs. pt vs. #Delta#phi vs. q_{2} vs. centr",naxes,nbins,xmin,xmax);
  }
  
  TString q2axisname="q_{2}^{TPC}";
  if(fq2Meth==kq2PosTPC) {q2axisname="q_{2}^{pos TPC}";}
  else if(fq2Meth==kq2NegTPC) {q2axisname="q_{2}^{neg TPC}";}
  else if(fq2Meth==kq2VZERO) {q2axisname="q_{2}^{V0}";}
  else if(fq2Meth==kq2VZEROA) {q2axisname="q_{2}^{V0A}";}
  else if(fq2Meth==kq2VZEROC) {q2axisname="q_{2}^{V0C}";}
  
  TString axTit[6]={massaxisname,"p_{T} (GeV/c)","#Delta#phi",q2axisname,"Centrality (%)","part-antipart"};
  
  for(Int_t iax=0; iax<naxes; iax++)
    fHistMassPtPhiq2Centr->GetAxis(iax)->SetTitle(axTit[iax].Data());
  
  fOutput->Add(fHistMassPtPhiq2Centr);
  
}

//***************************************************************************

// Methods used in the UserExec

void AliAnalysisTaskHFv1::CalculateInvMasses(AliAODRecoDecayHF* d,Float_t*& masses,Int_t& nmasses){
  //Calculates all the possible invariant masses for each candidate
  //NB: the implementation for each candidate is responsibility of the corresponding developer

  if(fDecChannel==0){
    //D+ -- Giacomo
    nmasses=1;
    masses=new Float_t[nmasses];
    Int_t pdgdaughters[3] = {211,321,211};
    masses[0]=d->InvMass(3,(UInt_t*)pdgdaughters);
  }
  if(fDecChannel==1){
    //D0 (Kpi)  -- Chiara
    const Int_t ndght=2;
    nmasses=2;
    masses=new Float_t[nmasses];
    Int_t pdgdaughtersD0[ndght]={211,321};//pi,K 
    masses[0]=d->InvMass(ndght,(UInt_t*)pdgdaughtersD0); //D0
    Int_t pdgdaughtersD0bar[ndght]={321,211};//K,pi 
    masses[1]=d->InvMass(ndght,(UInt_t*)pdgdaughtersD0bar); //D0bar
  }
  if(fDecChannel==2){
    //D* -- Robert,Yifei, Alessandro
    nmasses=1;
    masses=new Float_t[nmasses];
    masses[0]=((AliAODRecoCascadeHF*)d)->DeltaInvMass();
  } 
  if(fDecChannel==3){
    //Ds -- Anastasia
    nmasses=2;
    masses=new Float_t[nmasses];
    UInt_t pdgdaughtersKKpi[3] = {321,321,211};
    UInt_t pdgdaughterspiKK[3] = {211,321,321};
    masses[0]=d->InvMass(3,pdgdaughtersKKpi);
    masses[1]=d->InvMass(3,pdgdaughterspiKK);
  }
}

//******************************************************************************

//Methods to fill the histograms, one for each channel
//NB: the implementation for each candidate is responsibility of the corresponding developer

//******************************************************************************
void AliAnalysisTaskHFv1::FillDplus(AliAODRecoDecayHF* d,TClonesArray *arrayMC,Int_t ptbin,Float_t deltaphi, const Float_t* masses,Int_t isSel,Int_t icentr, Double_t phiD, Double_t etaD, Double_t ptD, Double_t QA[2], Double_t QB[2]){
  //D+ channel
  if(!isSel){
    if(fDebug>3)AliWarning("Candidate not selected\n");
    return;
  }
  if(!masses){
    if(fDebug>3)AliWarning("Masses not calculated\n");
    return;
  }
  Int_t icentrmin=icentr-fCentBinSizePerMil;
  if(fFlowMethod==kSP) {
    Double_t scalprod[2]={-2.,-2.};
    scalprod[0] = TMath::Cos(fHarmonic*phiD)*QB[0]+TMath::Sin(fHarmonic*phiD)*QB[1];
    Double_t sparsearray[5] = {scalprod[0],masses[0],ptD,etaD,0.5};
    ((THnSparseD*)fOutput->FindObject(Form("hMassScalProduQA_centr%d_%d",icentrmin,icentr)))->Fill(sparsearray);
    scalprod[1] = TMath::Cos(fHarmonic*phiD)*QA[0]+TMath::Sin(fHarmonic*phiD)*QA[1];
    sparsearray[0] = scalprod[1];
    ((THnSparseD*)fOutput->FindObject(Form("hMassScalProduQC_centr%d_%d",icentrmin,icentr)))->Fill(sparsearray);
  }
  else {
    ((TH2F*)fOutput->FindObject(Form("hMc2deltaphi_pt%dcentr%d_%d",ptbin,icentrmin,icentr)))->Fill(TMath::Cos(fHarmonic*deltaphi),masses[0]);
    ((TH2F*)fOutput->FindObject(Form("hMs2deltaphi_pt%dcentr%d_%d",ptbin,icentrmin,icentr)))->Fill(TMath::Sin(fHarmonic*deltaphi),masses[0]);
    ((TH2F*)fOutput->FindObject(Form("hMdeltaphi_pt%dcentr%d_%d",ptbin,icentrmin,icentr)))->Fill(deltaphi,masses[0]);
    ((TH2F*)fOutput->FindObject(Form("hCos2phiMass_pt%dcentr%d_%d",ptbin,icentrmin,icentr)))->Fill(TMath::Cos(fHarmonic*phiD),masses[0]);
    ((TH2F*)fOutput->FindObject(Form("hSin2phiMass_pt%dcentr%d_%d",ptbin,icentrmin,icentr)))->Fill(TMath::Sin(fHarmonic*phiD),masses[0]);
  }
  ((TH2F*)fOutput->FindObject(Form("hMPtCandcentr%d_%d",icentrmin,icentr)))->Fill(d->Pt(),masses[0]);

  Int_t pdgdaughters[3] = {211,321,211};

  if(fReadMC && fFlowMethod==kEP){
    Int_t lab=-1;
    if(fUseAfterBurner){
      Bool_t isSignal=fAfterBurner->GetIsSignal();
      if(isSignal)lab=10;
    }else {
      lab = d->MatchToMC(411,arrayMC,3,pdgdaughters);
    }
    if(lab>=0){ //signal
      ((TH2F*)fOutput->FindObject(Form("hMc2deltaphiS_pt%dcentr%d_%d",ptbin,icentrmin,icentr)))->Fill(TMath::Cos(fHarmonic*deltaphi),masses[0]);
      ((TH2F*)fOutput->FindObject(Form("hMdeltaphiS_pt%dcentr%d_%d",ptbin,icentrmin,icentr)))->Fill(deltaphi,masses[0]);
    } else{ //background
      ((TH2F*)fOutput->FindObject(Form("hMc2deltaphiB_pt%dcentr%d_%d",ptbin,icentrmin,icentr)))->Fill(TMath::Cos(fHarmonic*deltaphi),masses[0]);
      ((TH2F*)fOutput->FindObject(Form("hMdeltaphiB_pt%dcentr%d_%d",ptbin,icentrmin,icentr)))->Fill(deltaphi,masses[0]);
    } 
  }   
}

//******************************************************************************
void AliAnalysisTaskHFv1::FillD02p(AliAODRecoDecayHF* d,TClonesArray *arrayMC,Int_t ptbin,Float_t deltaphi, const Float_t* masses,Int_t isSel,Int_t icentr,Double_t phiD, Double_t etaD, Double_t ptD, Double_t QA[2], Double_t QB[2]){

  //D0->Kpi channel

  //mass histograms
  if(!masses){
    if(fDebug>3)AliWarning("Masses not calculated\n");
    return;
  }
  Int_t icentrmin=icentr-fCentBinSizePerMil;
  Double_t scalprod[2]={-2.,-2.};
  Double_t FillisSel = (Double_t)(isSel-0.5);
  if(isSel==1 || isSel==3) {
    if(fFlowMethod==kSP) {
      scalprod[0] = TMath::Cos(fHarmonic*phiD)*QB[0]+TMath::Sin(fHarmonic*phiD)*QB[1];
      Double_t sparsearray[5] = {scalprod[0],masses[0],ptD,etaD,FillisSel};
      ((THnSparseD*)fOutput->FindObject(Form("hMassScalProduQA_centr%d_%d",icentrmin,icentr)))->Fill(sparsearray);
      scalprod[1] = TMath::Cos(fHarmonic*phiD)*QA[0]+TMath::Sin(fHarmonic*phiD)*QA[1];
      sparsearray[0] = scalprod[1];
      ((THnSparseD*)fOutput->FindObject(Form("hMassScalProduQC_centr%d_%d",icentrmin,icentr)))->Fill(sparsearray);
    }
    else {
      ((TH2F*)fOutput->FindObject(Form("hMc2deltaphi_pt%dcentr%d_%d",ptbin,icentrmin,icentr)))->Fill(TMath::Cos(fHarmonic*deltaphi),masses[0]);
      ((TH2F*)fOutput->FindObject(Form("hMs2deltaphi_pt%dcentr%d_%d",ptbin,icentrmin,icentr)))->Fill(TMath::Sin(fHarmonic*deltaphi),masses[0]);
      ((TH2F*)fOutput->FindObject(Form("hMdeltaphi_pt%dcentr%d_%d",ptbin,icentrmin,icentr)))->Fill(deltaphi,masses[0]);
      ((TH2F*)fOutput->FindObject(Form("hCos2phiMass_pt%dcentr%d_%d",ptbin,icentrmin,icentr)))->Fill(TMath::Cos(fHarmonic*phiD),masses[0]);
      ((TH2F*)fOutput->FindObject(Form("hSin2phiMass_pt%dcentr%d_%d",ptbin,icentrmin,icentr)))->Fill(TMath::Sin(fHarmonic*phiD),masses[0]);
    }
    ((TH2F*)fOutput->FindObject(Form("hMPtCandcentr%d_%d",icentrmin,icentr)))->Fill(d->Pt(),masses[0]);
  }
  if(isSel>=2) {
    if(fFlowMethod==kSP) {
      scalprod[0] = TMath::Cos(fHarmonic*phiD)*QB[0]+TMath::Sin(fHarmonic*phiD)*QB[1];
      Double_t sparsearray[5] = {scalprod[0],masses[1],ptD,etaD,FillisSel};
      ((THnSparseD*)fOutput->FindObject(Form("hMassScalProduQA_centr%d_%d",icentrmin,icentr)))->Fill(sparsearray);
      scalprod[1] = TMath::Cos(fHarmonic*phiD)*QA[0]+TMath::Sin(fHarmonic*phiD)*QA[1];
      sparsearray[0] = scalprod[1];
      ((THnSparseD*)fOutput->FindObject(Form("hMassScalProduQC_centr%d_%d",icentrmin,icentr)))->Fill(sparsearray);
    }
    else {
      ((TH2F*)fOutput->FindObject(Form("hMc2deltaphi_pt%dcentr%d_%d",ptbin,icentrmin,icentr)))->Fill(TMath::Cos(fHarmonic*deltaphi),masses[1]);
      ((TH2F*)fOutput->FindObject(Form("hMs2deltaphi_pt%dcentr%d_%d",ptbin,icentrmin,icentr)))->Fill(TMath::Sin(fHarmonic*deltaphi),masses[1]);
      ((TH2F*)fOutput->FindObject(Form("hMdeltaphi_pt%dcentr%d_%d",ptbin,icentrmin,icentr)))->Fill(deltaphi,masses[1]);
      ((TH2F*)fOutput->FindObject(Form("hCos2phiMass_pt%dcentr%d_%d",ptbin,icentrmin,icentr)))->Fill(TMath::Cos(fHarmonic*phiD),masses[1]);
      ((TH2F*)fOutput->FindObject(Form("hSin2phiMass_pt%dcentr%d_%d",ptbin,icentrmin,icentr)))->Fill(TMath::Sin(fHarmonic*phiD),masses[1]);
    }
    ((TH2F*)fOutput->FindObject(Form("hMPtCandcentr%d_%d",icentrmin,icentr)))->Fill(d->Pt(),masses[1]);
  }
  
  //MC histograms
  if(fReadMC && fFlowMethod==kEP){

    Int_t matchtoMC=-1;

    //D0
    Int_t pdgdaughters[2];
    pdgdaughters[0]=211;//pi 
    pdgdaughters[1]=321;//K
    Int_t nprongs=2;
    Int_t absPdgMom=421;

    matchtoMC = d->MatchToMC(absPdgMom,arrayMC,nprongs,pdgdaughters);

    Int_t prongPdgPlus=421,prongPdgMinus=(-1)*421;
    if((isSel==1 || isSel==3)){ //D0
      if(matchtoMC>=0){
        AliAODMCParticle *dMC = (AliAODMCParticle*)arrayMC->At(matchtoMC);
        Int_t pdgMC = dMC->GetPdgCode();
	
        if(pdgMC==prongPdgPlus) {
          ((TH2F*)fOutput->FindObject(Form("hMc2deltaphiS_pt%dcentr%d_%d",ptbin,icentrmin,icentr)))->Fill(TMath::Cos(fHarmonic*deltaphi),masses[0]);
          ((TH2F*)fOutput->FindObject(Form("hMdeltaphiS_pt%dcentr%d_%d",ptbin,icentrmin,icentr)))->Fill(deltaphi,masses[0]);
        }
        else {
          ((TH2F*)fOutput->FindObject(Form("hMc2deltaphiR_pt%dcentr%d_%d",ptbin,icentrmin,icentr)))->Fill(TMath::Cos(fHarmonic*deltaphi),masses[0]);
          ((TH2F*)fOutput->FindObject(Form("hMdeltaphiR_pt%dcentr%d_%d",ptbin,icentrmin,icentr)))->Fill(deltaphi,masses[0]);
        }
      } else {
        ((TH2F*)fOutput->FindObject(Form("hMc2deltaphiB_pt%dcentr%d_%d",ptbin,icentrmin,icentr)))->Fill(TMath::Cos(fHarmonic*deltaphi),masses[0]);
        ((TH2F*)fOutput->FindObject(Form("hMdeltaphi_pt%dcentr%d_%d",ptbin,icentrmin,icentr)))->Fill(deltaphi,masses[0]);
      }
    }
    if(isSel>=2){ //D0bar
      if(matchtoMC>=0){
        AliAODMCParticle *dMC = (AliAODMCParticle*)arrayMC->At(matchtoMC);
        Int_t pdgMC = dMC->GetPdgCode();
	
        if(pdgMC==prongPdgMinus) {
          ((TH2F*)fOutput->FindObject(Form("hMc2deltaphiS_pt%dcentr%d_%d",ptbin,icentrmin,icentr)))->Fill(TMath::Cos(fHarmonic*deltaphi),masses[1]);
          ((TH2F*)fOutput->FindObject(Form("hMdeltaphiS_pt%dcentr%d_%d",ptbin,icentrmin,icentr)))->Fill(deltaphi,masses[1]);
        }
        else {
          ((TH2F*)fOutput->FindObject(Form("hMc2deltaphiR_pt%dcentr%d_%d",ptbin,icentrmin,icentr)))->Fill(TMath::Cos(fHarmonic*deltaphi),masses[1]);
          ((TH2F*)fOutput->FindObject(Form("hMdeltaphi_pt%dcentr%d_%d",ptbin,icentrmin,icentr)))->Fill(deltaphi,masses[1]);
        }
      } else {
        ((TH2F*)fOutput->FindObject(Form("hMc2deltaphiB_pt%dcentr%d_%d",ptbin,icentrmin,icentr)))->Fill(TMath::Cos(fHarmonic*deltaphi),masses[1]);
        ((TH2F*)fOutput->FindObject(Form("hMdeltaphiB_pt%dcentr%d_%d",ptbin,icentrmin,icentr)))->Fill(deltaphi,masses[1]);
      }
    }
  }
}

//******************************************************************************
void AliAnalysisTaskHFv1::FillDstar(AliAODRecoDecayHF* d,TClonesArray *arrayMC,Int_t ptbin,Float_t deltaphi, const Float_t* masses,Int_t isSel,Int_t icentr, Double_t phiD, Double_t etaD, Double_t ptD, Double_t QA[2], Double_t QB[2]){
  //D* channel
  if(!isSel){
    if(fDebug>3)AliWarning("Candidate not selected\n");
    return;
  }
  if(!masses){
    if(fDebug>3)AliWarning("Masses not calculated\n");
    return;
  }

  Float_t deltaphi1 = deltaphi;
  if(deltaphi1 > TMath::PiOver2()) deltaphi1 = TMath::Pi() - deltaphi1; 
  Int_t icentrmin=icentr-fCentBinSizePerMil;
  if(fFlowMethod==kSP) {
    Double_t scalprod[2]={-2.,-2.};
    scalprod[0] = TMath::Cos(fHarmonic*phiD)*QB[0]+TMath::Sin(fHarmonic*phiD)*QB[1];
    Double_t sparsearray[5] = {scalprod[0],masses[0],ptD,etaD,0.5};
    ((THnSparseD*)fOutput->FindObject(Form("hMassScalProduQA_centr%d_%d",icentrmin,icentr)))->Fill(sparsearray);
    scalprod[1] = TMath::Cos(fHarmonic*phiD)*QA[0]+TMath::Sin(fHarmonic*phiD)*QA[1];
    sparsearray[0] = scalprod[1];
    ((THnSparseD*)fOutput->FindObject(Form("hMassScalProduQC_centr%d_%d",icentrmin,icentr)))->Fill(sparsearray);
  }
  else {
    ((TH2F*)fOutput->FindObject(Form("hMc2deltaphi_pt%dcentr%d_%d",ptbin,icentrmin,icentr)))->Fill(TMath::Cos(fHarmonic*deltaphi),masses[0]);
    ((TH2F*)fOutput->FindObject(Form("hMs2deltaphi_pt%dcentr%d_%d",ptbin,icentrmin,icentr)))->Fill(TMath::Sin(fHarmonic*deltaphi),masses[0]);
    ((TH2F*)fOutput->FindObject(Form("hMdeltaphi_pt%dcentr%d_%d",ptbin,icentrmin,icentr)))->Fill(deltaphi1,masses[0]);
    ((TH2F*)fOutput->FindObject(Form("hCos2phiMass_pt%dcentr%d_%d",ptbin,icentrmin,icentr)))->Fill(TMath::Cos(fHarmonic*phiD),masses[0]);
    ((TH2F*)fOutput->FindObject(Form("hSin2phiMass_pt%dcentr%d_%d",ptbin,icentrmin,icentr)))->Fill(TMath::Sin(fHarmonic*phiD),masses[0]);
  }
  ((TH2F*)fOutput->FindObject(Form("hMPtCandcentr%d_%d",icentrmin,icentr)))->Fill(d->Pt(),masses[0]);

  Int_t pdgDgDStartoD0pi[2]={421,211};
  Int_t pdgDgD0toKpi[2]={321,211};
  
  if(fReadMC && fFlowMethod==kEP){
    Int_t lab=-1;
    lab = ((AliAODRecoCascadeHF*)d)->MatchToMC(413,421,pdgDgDStartoD0pi,pdgDgD0toKpi,arrayMC);
    if(lab>=0){ //signal
      ((TH2F*)fOutput->FindObject(Form("hMc2deltaphiS_pt%dcentr%d_%d",ptbin,icentrmin,icentr)))->Fill(TMath::Cos(fHarmonic*deltaphi),masses[0]);
      ((TH2F*)fOutput->FindObject(Form("hMdeltaphiS_pt%dcentr%d_%d",ptbin,icentrmin,icentr)))->Fill(deltaphi,masses[0]);
    } else{ //background
      ((TH2F*)fOutput->FindObject(Form("hMc2deltaphiB_pt%dcentr%d_%d",ptbin,icentrmin,icentr)))->Fill(TMath::Cos(fHarmonic*deltaphi),masses[0]);
      ((TH2F*)fOutput->FindObject(Form("hMdeltaphiB_pt%dcentr%d_%d",ptbin,icentrmin,icentr)))->Fill(deltaphi,masses[0]);
    } 
  }
}

//******************************************************************************
void AliAnalysisTaskHFv1::FillDs(AliAODRecoDecayHF* d,TClonesArray *arrayMC,Int_t ptbin,Float_t deltaphi, const Float_t* masses,Int_t isSel,Int_t icentr, Double_t phiD, Double_t etaD, Double_t ptD, Double_t QA[2], Double_t QB[2]){

  //Ds channel
  if(!isSel){
    if(fDebug>3)AliWarning("Candidate not selected\n");
    return;
  }
  if(!masses){
    if(fDebug>3)AliWarning("Masses not calculated\n");
    return;
  }
  Int_t icentrmin=icentr-fCentBinSizePerMil;
  Int_t isDsPhiKKpi = isSel&4;
  Int_t isDsPhipiKK = isSel&8;
  Double_t scalprod[2]={-2.,-2.};

  if(isDsPhiKKpi){
    if(fFlowMethod==kSP) {
      scalprod[0] = TMath::Cos(fHarmonic*phiD)*QB[0]+TMath::Sin(fHarmonic*phiD)*QB[1];
      Double_t sparsearray[5] = {scalprod[0],masses[0],ptD,etaD,0.5};
      ((THnSparseD*)fOutput->FindObject(Form("hMassScalProduQA_centr%d_%d",icentrmin,icentr)))->Fill(sparsearray);
      scalprod[1] = TMath::Cos(fHarmonic*phiD)*QA[0]+TMath::Sin(fHarmonic*phiD)*QA[1];
      sparsearray[0] = scalprod[1];
      ((THnSparseD*)fOutput->FindObject(Form("hMassScalProduQC_centr%d_%d",icentrmin,icentr)))->Fill(sparsearray);
    }
    else {
      ((TH2F*)fOutput->FindObject(Form("hMc2deltaphi_pt%dcentr%d_%d",ptbin,icentrmin,icentr)))->Fill(TMath::Cos(fHarmonic*deltaphi),masses[0]);
      ((TH2F*)fOutput->FindObject(Form("hMs2deltaphi_pt%dcentr%d_%d",ptbin,icentrmin,icentr)))->Fill(TMath::Sin(fHarmonic*deltaphi),masses[0]);
      ((TH2F*)fOutput->FindObject(Form("hMdeltaphi_pt%dcentr%d_%d",ptbin,icentrmin,icentr)))->Fill(deltaphi,masses[0]);
      ((TH2F*)fOutput->FindObject(Form("hCos2phiMass_pt%dcentr%d_%d",ptbin,icentrmin,icentr)))->Fill(TMath::Cos(fHarmonic*phiD),masses[0]);
      ((TH2F*)fOutput->FindObject(Form("hSin2phiMass_pt%dcentr%d_%d",ptbin,icentrmin,icentr)))->Fill(TMath::Sin(fHarmonic*phiD),masses[0]);
    }
    ((TH2F*)fOutput->FindObject(Form("hMPtCandcentr%d_%d",icentrmin,icentr)))->Fill(d->Pt(),masses[0]);
  }
  if(isDsPhipiKK){
    if(fFlowMethod==kSP) {
      scalprod[0] = TMath::Cos(fHarmonic*phiD)*QB[0]+TMath::Sin(fHarmonic*phiD)*QB[1];
      Double_t sparsearray[5] = {scalprod[0],masses[1],ptD,etaD,0.5};
      ((THnSparseD*)fOutput->FindObject(Form("hMassScalProduQA_centr%d_%d",icentrmin,icentr)))->Fill(sparsearray);
      scalprod[1] = TMath::Cos(fHarmonic*phiD)*QA[0]+TMath::Sin(fHarmonic*phiD)*QA[1];
      sparsearray[0] = scalprod[1];
      ((THnSparseD*)fOutput->FindObject(Form("hMassScalProduQC_centr%d_%d",icentrmin,icentr)))->Fill(sparsearray);
    }
    else {
      ((TH2F*)fOutput->FindObject(Form("hMc2deltaphi_pt%dcentr%d_%d",ptbin,icentrmin,icentr)))->Fill(TMath::Cos(fHarmonic*deltaphi),masses[1]);
      ((TH2F*)fOutput->FindObject(Form("hMs2deltaphi_pt%dcentr%d_%d",ptbin,icentrmin,icentr)))->Fill(TMath::Sin(fHarmonic*deltaphi),masses[1]);
      ((TH2F*)fOutput->FindObject(Form("hMdeltaphi_pt%dcentr%d_%d",ptbin,icentrmin,icentr)))->Fill(deltaphi,masses[1]);
      ((TH2F*)fOutput->FindObject(Form("hCos2phiMass_pt%dcentr%d_%d",ptbin,icentrmin,icentr)))->Fill(TMath::Cos(fHarmonic*phiD),masses[1]);
      ((TH2F*)fOutput->FindObject(Form("hSin2phiMass_pt%dcentr%d_%d",ptbin,icentrmin,icentr)))->Fill(TMath::Sin(fHarmonic*phiD),masses[1]);
    }
    ((TH2F*)fOutput->FindObject(Form("hMPtCandcentr%d_%d",icentrmin,icentr)))->Fill(d->Pt(),masses[1]);
  }

  if(fReadMC && fFlowMethod==kEP){

    Int_t pdgdaughters[3] = {321,321,211};

    Int_t lab=-1;
    Int_t pdgCode0=-1;

    lab = d->MatchToMC(431,arrayMC,3,pdgdaughters);
    if(lab>=0){ //signal
      AliAODMCParticle *dMC = (AliAODMCParticle*)arrayMC->At(lab);
      Int_t labDau0=((AliAODTrack*)d->GetDaughter(0))->GetLabel();
      AliAODMCParticle* p=(AliAODMCParticle*)arrayMC->UncheckedAt(TMath::Abs(labDau0));
      pdgCode0=TMath::Abs(p->GetPdgCode());
    }
    if(isDsPhiKKpi){
      if(lab>=0){
        if(pdgCode0==321){
          ((TH2F*)fOutput->FindObject(Form("hMc2deltaphiS_pt%dcentr%d_%d",ptbin,icentrmin,icentr)))->Fill(TMath::Cos(fHarmonic*deltaphi),masses[0]);
          ((TH2F*)fOutput->FindObject(Form("hMdeltaphiS_pt%dcentr%d_%d",ptbin,icentrmin,icentr)))->Fill(deltaphi,masses[0]);
        }else{
          ((TH2F*)fOutput->FindObject(Form("hMc2deltaphiR_pt%dcentr%d_%d",ptbin,icentrmin,icentr)))->Fill(TMath::Cos(fHarmonic*deltaphi),masses[0]);
          ((TH2F*)fOutput->FindObject(Form("hMdeltaphiR_pt%dcentr%d_%d",ptbin,icentrmin,icentr)))->Fill(deltaphi,masses[0]);
        }
      }else{
        ((TH2F*)fOutput->FindObject(Form("hMc2deltaphiB_pt%dcentr%d_%d",ptbin,icentrmin,icentr)))->Fill(TMath::Cos(fHarmonic*deltaphi),masses[0]);
        ((TH2F*)fOutput->FindObject(Form("hMdeltaphiB_pt%dcentr%d_%d",ptbin,icentrmin,icentr)))->Fill(deltaphi,masses[0]);
      }
    }
    if(isDsPhipiKK){
      if(lab>=0){
        if(pdgCode0==211){
          ((TH2F*)fOutput->FindObject(Form("hMc2deltaphiS_pt%dcentr%d_%d",ptbin,icentrmin,icentr)))->Fill(TMath::Cos(fHarmonic*deltaphi),masses[1]);
          ((TH2F*)fOutput->FindObject(Form("hMdeltaphiS_pt%dcentr%d_%d",ptbin,icentrmin,icentr)))->Fill(deltaphi,masses[1]);
        }else{
          ((TH2F*)fOutput->FindObject(Form("hMc2deltaphiR_pt%dcentr%d_%d",ptbin,icentrmin,icentr)))->Fill(TMath::Cos(fHarmonic*deltaphi),masses[1]);
          ((TH2F*)fOutput->FindObject(Form("hMdeltaphiR_pt%dcentr%d_%d",ptbin,icentrmin,icentr)))->Fill(deltaphi,masses[1]);
        }
      }else{
        ((TH2F*)fOutput->FindObject(Form("hMc2deltaphiB_pt%dcentr%d_%d",ptbin,icentrmin,icentr)))->Fill(TMath::Cos(fHarmonic*deltaphi),masses[1]);
        ((TH2F*)fOutput->FindObject(Form("hMdeltaphiB_pt%dcentr%d_%d",ptbin,icentrmin,icentr)))->Fill(deltaphi,masses[1]);
      }
    }
  }
}

//________________________________________________________________________
void AliAnalysisTaskHFv1::SetEventPlaneMethod(Int_t method){
  // Interface method to set some specific cases

  switch (method){
  case kTPC:
    {
      fEvPlaneDet=kPosTPC;
      fSubEvDetA=kPosTPC;
      fSubEvDetB=kNegTPC;
      break;
    }
  case kTPCVZERO:
    {
      fEvPlaneDet=kFullTPC;
      fSubEvDetA=kV0A;
      fSubEvDetB=kV0C;
      break;
    }
  case kVZERO:
    {
      fEvPlaneDet=kFullV0;
      fSubEvDetA=kPosTPC;
      fSubEvDetB=kNegTPC;
      break;
    }
  case kVZEROA:
    {
      fEvPlaneDet=kV0A;
      fSubEvDetA=kV0C;
      fSubEvDetB=kFullTPC;
      break;
    }
  case kVZEROC:
    {
      fEvPlaneDet=kV0C;
      fSubEvDetA=kV0A;
      fSubEvDetB=kFullTPC;
      break;
    }
  case kPosTPCVZERO:
    {
      fEvPlaneDet=kPosTPC;
      fSubEvDetA=kV0A;
      fSubEvDetB=kV0C;
      break;
    }
  case kNegTPCVZERO:
    {
      fEvPlaneDet=kNegTPC;
      fSubEvDetA=kV0A;
      fSubEvDetB=kV0C;
      break;
    }
  case kZDC:
    {
      fEvPlaneDet=kZDCA;
      fSubEvDetA=kZDCC;
      fSubEvDetB=kFullTPC;
      break;
    }
  default:
    {
      AliWarning("Case not defined -> use full V0");
      fEvPlaneDet=kFullV0;
      fSubEvDetA=kPosTPC;
      fSubEvDetB=kNegTPC;
      break;
    }
  }
}
//________________________________________________________________________
Float_t AliAnalysisTaskHFv1::GetPhiInRange(Float_t phi){
  // Sets the phi angle in the range 0-pi
  Float_t result=phi;
  while(result<0){
    result=result+2.*TMath::Pi()/fHarmonic;
  }
  while(result>2.*TMath::Pi()/fHarmonic){
    result=result-2.*TMath::Pi()/fHarmonic;
  }
  return result;
}

//________________________________________________________________________
Float_t AliAnalysisTaskHFv1::GetEventPlane(AliAODEvent* aod, AliEventplane *pl, Double_t eventplaneqncorrTPC[3], Double_t eventplaneqncorrVZERO[3], Double_t q2){
  //Event plane

  Double_t rpangleTPC=0;
  Double_t rpangleTPCpos=0;
  Double_t rpangleTPCneg=0;
  Double_t rpangleVZERO=0;
  Double_t rpangleVZEROA=0;
  Double_t rpangleVZEROC=0;
  Double_t eventplane=-9999.;
  if(fUseNewQnCorrFw) {
    rpangleTPC = eventplaneqncorrTPC[0];
    rpangleTPCpos = eventplaneqncorrTPC[1];
    rpangleTPCneg = eventplaneqncorrTPC[2];
    rpangleVZERO=eventplaneqncorrVZERO[0];
    rpangleVZEROA=eventplaneqncorrVZERO[1];
    rpangleVZEROC=eventplaneqncorrVZERO[2];
  }else{ 
    rpangleTPC = pl->GetEventplane("Q");
    rpangleVZERO=GetPhiInRange(pl->GetEventplane("V0",aod,fHarmonic));
    rpangleVZEROA= GetPhiInRange(pl->GetEventplane("V0A",aod,fHarmonic));
    rpangleVZEROC= GetPhiInRange(pl->GetEventplane("V0C",aod,fHarmonic));
    TVector2* qsub1 = pl->GetQsub1();
    TVector2* qsub2 = pl->GetQsub2();
    if(!qsub1 || !qsub2){
      fhEventsInfo->Fill(11);
      return-9999.;
    }
    rpangleTPCpos = qsub1->Phi()/2.;
    rpangleTPCneg = qsub2->Phi()/2.;
  }
  if(fOnTheFlyTPCEP){
    Double_t rpangleTPCold=rpangleTPC;
    ComputeTPCEventPlane(aod,rpangleTPC,rpangleTPCpos,rpangleTPCneg);
    ((TH2F*)fOutput->FindObject(Form("hEvPlaneRecomp%s",fCentrBinName.Data())))->Fill(rpangleTPCold,rpangleTPC);
  }
  if(rpangleTPC<0){
    fhEventsInfo->Fill(11);
    return -9999.;
  }

  if(TMath::Abs(rpangleTPC-rpangleVZERO)>fEventPlanesComp) return -9999.;

  Bool_t set=kFALSE;
  Double_t rpangleeventC=0;
  Double_t rpangleeventB=0;
  Double_t rpangleeventA=0;
  Int_t nSubEvents=2;
  if(fEvPlaneDet==kFullTPC){
    eventplane=rpangleTPC;
    if((fSubEvDetA==kV0A && fSubEvDetB==kV0C)||
       (fSubEvDetB==kV0A && fSubEvDetA==kV0C)){
      nSubEvents=3;
      rpangleeventA=rpangleVZEROA;
      rpangleeventB=rpangleVZEROC;
      rpangleeventC=rpangleTPC;
      set=kTRUE;
    }else if((fSubEvDetA==kPosTPC && fSubEvDetB==kNegTPC)||
       (fSubEvDetB==kPosTPC && fSubEvDetA==kNegTPC)){
      nSubEvents=2;
      eventplane=rpangleTPC;
      rpangleeventA=rpangleTPCpos;
      rpangleeventB=rpangleTPCneg;
      set=kTRUE;
    }
  }else if(fEvPlaneDet==kPosTPC){
    eventplane=rpangleTPCpos;
    if(!fUseNewQnCorrFw) eventplane=rpangleTPC;
    if((fSubEvDetA==kV0A && fSubEvDetB==kV0C)||
       (fSubEvDetB==kV0A && fSubEvDetA==kV0C)){
      nSubEvents=3;
      rpangleeventA=rpangleVZEROA;
      rpangleeventB=rpangleVZEROC;
      rpangleeventC=rpangleTPCpos;
      set=kTRUE;
    }else if((fSubEvDetA==kPosTPC && fSubEvDetB==kNegTPC)||
       (fSubEvDetB==kPosTPC && fSubEvDetA==kNegTPC)){
      nSubEvents=2;
      rpangleeventA=rpangleTPCneg;
      rpangleeventB=rpangleTPCpos;
      set=kTRUE;
    }
  }else if(fEvPlaneDet==kNegTPC){
    eventplane=rpangleTPCneg;
    if(!fUseNewQnCorrFw) eventplane=rpangleTPC;
    if((fSubEvDetA==kV0A && fSubEvDetB==kV0C)||
       (fSubEvDetB==kV0A && fSubEvDetA==kV0C)){
      nSubEvents=3;
      rpangleeventA=rpangleVZEROA;
      rpangleeventB=rpangleVZEROC;
      rpangleeventC=rpangleTPCneg;
      set=kTRUE;
    }else if((fSubEvDetA==kPosTPC && fSubEvDetB==kNegTPC)||
       (fSubEvDetB==kPosTPC && fSubEvDetA==kNegTPC)){
      nSubEvents=2;
      rpangleeventA=rpangleTPCpos;
      rpangleeventB=rpangleTPCneg;
      set=kTRUE;
    }
  }else if(fEvPlaneDet==kFullV0){
    eventplane=rpangleVZERO;
    if((fSubEvDetA==kPosTPC && fSubEvDetB==kNegTPC)||
       (fSubEvDetB==kPosTPC && fSubEvDetA==kNegTPC)){
      nSubEvents=3;
      rpangleeventA=rpangleTPCpos;
      rpangleeventB=rpangleTPCneg;
      rpangleeventC=rpangleVZERO;
      set=kTRUE;
    }
  }else if(fEvPlaneDet==kV0A){
    eventplane=rpangleVZEROA;
    if((fSubEvDetA==kPosTPC && fSubEvDetB==kNegTPC)||
       (fSubEvDetB==kPosTPC && fSubEvDetA==kNegTPC)){
      nSubEvents=3;
      rpangleeventA=rpangleTPCpos;
      rpangleeventB=rpangleTPCneg;
      rpangleeventC=rpangleVZEROA;
      set=kTRUE;
    }
    else if((fSubEvDetA==kFullTPC && fSubEvDetB==kV0C) ||
       (fSubEvDetB==kFullTPC && fSubEvDetA==kV0C)){
      nSubEvents=3;
      rpangleeventA=rpangleVZEROC;
      rpangleeventB=rpangleTPC;
      rpangleeventC=rpangleVZEROA;
      set=kTRUE;
    }
  }else if(fEvPlaneDet==kV0C){
    eventplane=rpangleVZEROC;
    if((fSubEvDetA==kPosTPC && fSubEvDetB==kNegTPC)||
       (fSubEvDetB==kPosTPC && fSubEvDetA==kNegTPC)){
      nSubEvents=3;
      rpangleeventA=rpangleTPCpos;
      rpangleeventB=rpangleTPCneg;
      rpangleeventC=rpangleVZEROC;
      set=kTRUE;
    }
    else if((fSubEvDetA==kFullTPC && fSubEvDetB==kV0A) ||
       (fSubEvDetB==kFullTPC && fSubEvDetA==kV0A)){
      nSubEvents=3;
      rpangleeventA=rpangleVZEROA;
      rpangleeventB=rpangleTPC;
      rpangleeventC=rpangleVZEROC;
      set=kTRUE;
    }
  }
  if(!set){
    AliError("The requested sub event configuration is not implemented!\n");
    return -999.;
  }

  Double_t deltaPsi =rpangleeventA-rpangleeventB;  
  if(TMath::Abs(deltaPsi)>TMath::Pi()/fHarmonic){
    if(deltaPsi>0.) deltaPsi-=2.*TMath::Pi()/fHarmonic;
    else deltaPsi+=2.*TMath::Pi()/fHarmonic;
  } // difference of subevents reaction plane angle cannot be bigger than pi/2 (pi/3)
  Double_t planereso = TMath::Cos(fHarmonic*deltaPsi); // reaction plane resolution

  if(fDebug>2)printf("Filling EP-related histograms\n");
  //Filling EP-related histograms
  if(fFlowMethod!=kEvShape) {
    ((TH2F*)fOutput->FindObject(Form("hEvPlane%s",fCentrBinName.Data())))->Fill(rpangleTPC,rpangleVZERO); // reaction plane angle without autocorrelations removal
    ((TH1F*)fOutput->FindObject(Form("hEvPlaneReso1%s",fCentrBinName.Data())))->Fill(planereso); //RP resolution
    ((TH1F*)fOutput->FindObject(Form("hCos2EP%s",fCentrBinName.Data())))->Fill(TMath::Cos(fHarmonic*eventplane)); // syst check
    ((TH1F*)fOutput->FindObject(Form("hSin2EP%s",fCentrBinName.Data())))->Fill(TMath::Sin(fHarmonic*eventplane)); // syst check    
    ((TH1F*)fOutput->FindObject(Form("hEvPlaneA%s",fCentrBinName.Data())))->Fill(rpangleeventA); //Angle of first subevent
    ((TH1F*)fOutput->FindObject(Form("hEvPlaneB%s",fCentrBinName.Data())))->Fill(rpangleeventB); //Angle of second subevent
  }
  else {
    ((TH1F*)fOutput->FindObject(Form("hEvPlaneReso1Vsq2%s",fCentrBinName.Data())))->Fill(planereso,q2); //RP resolution vs q2
  }
  
  if(nSubEvents==3){
    Double_t deltaSubAC=rpangleeventA-rpangleeventC;
    if(TMath::Abs(deltaSubAC)>TMath::Pi()/fHarmonic){
      // difference of subevents reaction plane angle cannot be bigger than phi/2
      if(deltaSubAC>0.) deltaSubAC-=2.*TMath::Pi()/fHarmonic;
      else deltaSubAC +=2.*TMath::Pi()/fHarmonic;
    } 
    Double_t deltaSubBC=rpangleeventB-rpangleeventC;
    if(TMath::Abs(deltaSubBC)>TMath::Pi()){
      // difference of subevents reaction plane angle cannot be bigger than phi/2
      if(deltaSubBC>0.) deltaSubBC-=2.*TMath::Pi()/fHarmonic;
      else deltaSubBC +=2.*TMath::Pi()/fHarmonic;
    } 
    if(fFlowMethod!=kEvShape) {
      TH1F* htmp1=(TH1F*)fOutput->FindObject(Form("hEvPlaneReso2%s",fCentrBinName.Data()));
      if(htmp1) htmp1->Fill(TMath::Cos(fHarmonic*deltaSubAC)); //RP resolution
      TH1F* htmp2=(TH1F*)fOutput->FindObject(Form("hEvPlaneReso3%s",fCentrBinName.Data()));
      if(htmp2) htmp2->Fill(TMath::Cos(fHarmonic*deltaSubBC)); //RP resolution
    }
    else {
      TH2F* htmp1=(TH2F*)fOutput->FindObject(Form("hEvPlaneReso2Vsq2%s",fCentrBinName.Data()));
      if(htmp1) htmp1->Fill(TMath::Cos(fHarmonic*deltaSubAC),q2); //RP resolution vs q2
      TH2F* htmp2=(TH2F*)fOutput->FindObject(Form("hEvPlaneReso3Vsq2%s",fCentrBinName.Data()));
      if(htmp2) htmp2->Fill(TMath::Cos(fHarmonic*deltaSubBC),q2); //RP resolution vs q2
    }
  }
  return eventplane;
}
//________________________________________________________________________
void AliAnalysisTaskHFv1::ComputeTPCEventPlane(AliAODEvent* aod, Double_t &rpangleTPC, Double_t &rpangleTPCpos,Double_t &rpangleTPCneg) const {
  Int_t nTracks=aod->GetNumberOfTracks();
  Double_t qVec[2]={0.,0.};
  Double_t qVecPosEta[2]={0.,0.};
  Double_t qVecNegEta[2]={0.,0.};
  for(Int_t it=0; it<nTracks; it++){
    AliAODTrack* track=(AliAODTrack*)aod->GetTrack(it);
    if(!track) continue;
    if(track->TestFilterBit(BIT(8))||track->TestFilterBit(BIT(9))) {
      Double_t eta=track->Eta();
      if(fEtaGapInTPCHalves>0. && TMath::Abs(eta)<fEtaGapInTPCHalves) continue;
      Double_t phi=track->Phi();
      Double_t wi=1.;
      if(fUsePtWeights){
	Double_t pt=track->Pt();
	wi=pt;
	if(pt>2) wi=2.;
      }
      Double_t qx=wi*TMath::Cos(fHarmonic*phi);
      Double_t qy=wi*TMath::Sin(fHarmonic*phi);
      qVec[0]+=qx;
      qVec[1]+=qy;
      if(eta>0){
	qVecPosEta[0]+=qx;
	qVecPosEta[1]+=qy;
      }else{
	qVecNegEta[0]+=qx;
	qVecNegEta[1]+=qy;
      }
    }
  }
  rpangleTPC = (TMath::Pi()+TMath::ATan2(-qVec[1],-qVec[0]))/2;
  rpangleTPCpos = (TMath::Pi()+TMath::ATan2(-qVecPosEta[1],-qVecPosEta[0]))/2;
  rpangleTPCneg = (TMath::Pi()+TMath::ATan2(-qVecNegEta[1],-qVecNegEta[0]))/2;
  return;
}
//________________________________________________________________________
Float_t AliAnalysisTaskHFv1::GetEventPlaneForCandidate(AliAODRecoDecayHF* d, AliEventplane *pl){
  // remove autocorrelations 
 
  TArrayF* qx = 0x0;
  TArrayF* qy = 0x0;
  TVector2* q = pl->GetQVector();
  TVector2* qsub1=pl->GetQsub1();
  TVector2* qsub2=pl->GetQsub2();
  TVector2 qcopy;   

  if(!fEtaGap){
    qx = pl->GetQContributionXArray();
    qy = pl->GetQContributionYArray();
    qcopy = *q;
  }else {
    if(d->Eta()<0.){
      qx = pl->GetQContributionXArraysub1();
      qy = pl->GetQContributionYArraysub1();
      qcopy = *qsub1;
    }else{
      qx = pl->GetQContributionXArraysub2();
      qy = pl->GetQContributionYArraysub2();
      qcopy = *qsub2;
    }
  }
  
  if(fDecChannel==2){
    //D* -- Yifei, Alessandro,Robert
    AliAODRecoDecayHF2Prong* theD0particle = ((AliAODRecoCascadeHF*)d)->Get2Prong();
    AliAODTrack *track0 = (AliAODTrack*)theD0particle->GetDaughter(0);
    AliAODTrack *track1 = (AliAODTrack*)theD0particle->GetDaughter(1);  
    AliAODTrack *track2 = ((AliAODRecoCascadeHF*)d)->GetBachelor(); 
    // reduce global q vector

    TVector2 q0;
    if((track0->GetID()) < qx->fN){
      q0.Set(qx->At(track0->GetID()),qy->At(track0->GetID()));}
	
    TVector2 q1;
    if((track1->GetID()) < qx->fN){
      q1.Set(qx->At(track1->GetID()),qy->At(track1->GetID()));}
	
    TVector2 q2;
    if((track2->GetID()) < qx->fN){
      q2.Set(qx->At(track2->GetID()),qy->At(track2->GetID()));}
      
    qcopy = qcopy -(q0+q1+q2);
  
  }
  
  // reduce Q vector for D+ and D0
  
  if(fDecChannel==1){    
    //D0 -- Chiara
    AliAODTrack *track0 = (AliAODTrack*)d->GetDaughter(0);
    AliAODTrack *track1 = (AliAODTrack*)d->GetDaughter(1);  

    TVector2 q0;
    if((track0->GetID()) < qx->fN){
      q0.Set(qx->At(track0->GetID()),qy->At(track0->GetID()));}
	
    TVector2 q1;
    if((track1->GetID()) < qx->fN){
      q1.Set(qx->At(track1->GetID()),qy->At(track1->GetID()));}
	
    qcopy = qcopy -(q0+q1);
  }

  if(fDecChannel==0){
    //D+ -- Giacomo
    AliAODTrack *track0 = (AliAODTrack*)d->GetDaughter(0);
    AliAODTrack *track1 = (AliAODTrack*)d->GetDaughter(1);  
    AliAODTrack *track2 = (AliAODTrack*)d->GetDaughter(2);  
    
    TVector2 q0;
    if((track0->GetID()) < qx->fN){
      q0.Set(qx->At(track0->GetID()),qy->At(track0->GetID()));}
	
    TVector2 q1;
    if((track1->GetID()) < qx->fN){
      q1.Set(qx->At(track1->GetID()),qy->At(track1->GetID()));}
	
    TVector2 q2;
    if((track2->GetID()) < qx->fN){
      q2.Set(qx->At(track2->GetID()),qy->At(track2->GetID()));}
      
    qcopy = qcopy -(q0+q1+q2);
	
  }

  if(fDecChannel==3){
    //Ds -- Anastasia
    AliAODTrack *track0 = (AliAODTrack*)d->GetDaughter(0);
    AliAODTrack *track1 = (AliAODTrack*)d->GetDaughter(1);
    AliAODTrack *track2 = (AliAODTrack*)d->GetDaughter(2);

    TVector2 q0;
    if((track0->GetID()) < qx->fN){
      q0.Set(qx->At(track0->GetID()),qy->At(track0->GetID()));}

    TVector2 q1;
    if((track1->GetID()) < qx->fN){
      q1.Set(qx->At(track1->GetID()),qy->At(track1->GetID()));}

    TVector2 q2;
    if((track2->GetID()) < qx->fN){
      q2.Set(qx->At(track2->GetID()),qy->At(track2->GetID()));}

    qcopy = qcopy -(q0+q1+q2);

  }

  return qcopy.Phi()/2.;

}
// //________________________________________________________________________
// Float_t AliAnalysisTaskHFv1::GetEventPlaneFromV0(AliAODEvent *aodEvent){
//   // Compute event plane for VZERO - Obsolete, used for 2010 data

//   Int_t centr=fRDCuts->GetCentrality(aodEvent);
//   centr=centr-centr%10;
//   //temporary fix
//   if(centr<20)centr=20;
//   if(centr>70)centr=70;
//   //end temporary fix
//   Int_t binx=0;
//   Int_t iParHist=(centr-20)/10;

//   TString name;name.Form("parhist%d_%d",centr,centr+10);

//   if(fDebug>15)printf("EPfromV0 centr %d, iparhist %d (%p-%p)\n",centr,iParHist,fParHist->FindObject(name.Data()),fParHist->At(iParHist));

//   Int_t runnumber=aodEvent->GetRunNumber();
//   if(fParHist->At(iParHist)){
//     for(Int_t i=1;i<=((TH2D*)fParHist->At(iParHist))->GetNbinsX()&&binx<=0;i++){
//       Int_t run=atoi(((TH2D*)fParHist->At(iParHist))->GetXaxis()->GetBinLabel(i));
//       if(run>=runnumber)binx=i;
//     }
//   }else{
//     fhEventsInfo->Fill(7);
//   }

//   AliFlowTrackCuts* cutsRP = AliFlowTrackCuts::GetStandardVZEROOnlyTrackCuts();
//   cutsRP->SetEvent(aodEvent, MCEvent());//, 0x0);
//   cutsRP->SetName( Form("rp_cuts") );
//   AliFlowTrackCuts* dummy = new AliFlowTrackCuts("null_cuts");
//   dummy->SetParamType(AliFlowTrackCuts::kGlobal);
//   dummy->SetPtRange(+1,-1); // select nothing QUICK
//   dummy->SetEtaRange(+1,-1); // select nothing VZERO
//   dummy->SetEvent(aodEvent,MCEvent());

//   //////////////// construct the flow event container ////////////
//   AliFlowEvent flowEvent(cutsRP,dummy);
//   flowEvent.SetReferenceMultiplicity( 64 );
//   for(Int_t i=0;i<64&&binx>0;i++){
//     AliFlowTrack *flowTrack=flowEvent.GetTrack(i);
//     Double_t inte=((TH2D*)fParHist->At(iParHist))->Integral(binx,binx,i+1,i+1);
//     if(inte>0)flowTrack->SetWeight(flowTrack->Weight()/inte);
//   }
//   if(fDebug>15)printf("EPfromV0 flow tracks weights done\n");

//   AliFlowVector qvec=flowEvent.GetQ(fHarmonic);
//   Double_t angleEP=(1./(Double_t)fHarmonic)*qvec.Phi();
//   if(fDebug>15)printf("EPfromV0 phi %f\n",angleEP);
//   return angleEP;
// }


//________________________________________________________________________
Float_t AliAnalysisTaskHFv1::GetEventPlaneForCandidateNewQnFw(AliAODRecoDecayHF* d, const TList *qnlist) {
    
  AliQnCorrectionsQnVector *theQnVectorCorr   = NULL;
  AliQnCorrectionsQnVector *theQnVectorUncorr = NULL;
  Float_t corrX, corrY;
  Float_t eventplane = -999.;
  Float_t ux = 0., uy = 0.; // to store qx, qy of the candidate's prongs
  Float_t qxRec = -999., qyRec = -999.;     // to store Qx, Qy at the rec step
  
  TList *pQvecList = dynamic_cast<TList*> (qnlist->FindObject("TPC"));
  if (pQvecList != NULL) {
    /* the detector is present */
    theQnVectorUncorr = (AliQnCorrectionsQnVector*) pQvecList->FindObject("plain"); //raw step for TPC
    if (theQnVectorUncorr == NULL || !(theQnVectorUncorr->IsGoodQuality()) || theQnVectorUncorr->GetN() == 0) return eventplane;
    Float_t qxPlain = theQnVectorUncorr->Qx(fHarmonic);
    Float_t qyPlain = theQnVectorUncorr->Qy(fHarmonic);
    
    //auto-correlation removal from rec step (if present)
    theQnVectorCorr   = (AliQnCorrectionsQnVector*) pQvecList->FindObject("rec"); //rec step for TPC
    if (theQnVectorCorr == NULL || !(theQnVectorCorr->IsGoodQuality()) || theQnVectorUncorr->GetN() == 0) {
      corrX = 0.;
      corrY = 0.;
      
    }
    else {
      qxRec = theQnVectorCorr->Qx(2);
      qyRec = theQnVectorCorr->Qy(2);
      corrX = qxPlain - qxRec;
      corrY = qyPlain - qyRec;
    }
    Int_t M = theQnVectorUncorr->GetN(); //# of tracks
    
    Int_t nProngs = 0;
    if(fDecChannel==0) nProngs = 3;      //D+
    else if(fDecChannel==1) nProngs = 2; //D0
    else if(fDecChannel==2) nProngs = 3; //D*
    else if(fDecChannel==3) nProngs = 3; //Ds
    
    AliAODTrack *track;
    for(int i = 0; i < nProngs; i++) {
      track = (AliAODTrack*)d->GetDaughter(i);
      if(track->TestFilterBit(BIT(8))||track->TestFilterBit(BIT(9))) {
	ux += TMath::Cos(fHarmonic*track->Phi());
	uy += TMath::Sin(fHarmonic*track->Phi());
      }
    }
    Float_t qxRecSub = (qxPlain*M - ux)/(M-nProngs) - corrX;
    Float_t qyRecSub = (qyPlain*M - uy)/(M-nProngs) - corrY;
        
    if(fHarmonic==2) eventplane = (TMath::Pi()+TMath::ATan2(-qyRecSub,-qxRecSub))/2;
    if(fHarmonic==3){
      eventplane = TMath::ATan2(qyRecSub,qxRecSub)/fHarmonic;
      if(eventplane<0.) eventplane+=2.*(TMath::Pi())/3.;
    }
    //auto-correlation removal from twist step (if present)
    theQnVectorCorr   = (AliQnCorrectionsQnVector*) pQvecList->FindObject("twist"); //twist step for TPC
    if (theQnVectorCorr == NULL || !(theQnVectorCorr->IsGoodQuality()) || theQnVectorUncorr->GetN() == 0) {
      return eventplane;
    }
    else {
      Float_t qxTwist = theQnVectorCorr->Qx(fHarmonic);
      Float_t qyTwist = theQnVectorCorr->Qy(fHarmonic);
      Float_t Lbplus  = (qyRec - qyTwist)/qxTwist;
      Float_t Lbminus = (qxRec - qxTwist)/qyTwist;
      Float_t qxTwistSub = (qxRecSub - Lbminus*qyRecSub)/(1-Lbminus*Lbplus);
      Float_t qyTwistSub = (qyRecSub - Lbplus*qxRecSub)/(1-Lbminus*Lbplus);
      
      if(fHarmonic==2) eventplane = (TMath::Pi()+TMath::ATan2(-qyTwistSub,-qxTwistSub))/2;
      if(fHarmonic==3){
	eventplane = TMath::ATan2(qyTwistSub,qxTwistSub)/fHarmonic;
	if(eventplane<0.) eventplane+=2.*(TMath::Pi())/3.;
      }
    }
  }
  return eventplane;
}

//________________________________________________________________________
const AliQnCorrectionsQnVector *AliAnalysisTaskHFv1::GetQnVectorFromList(const TList *list,
									   const char *subdetector,
									   const char *expectedstep,
									   const char *altstep)
{    
  AliQnCorrectionsQnVector *theQnVector = NULL;
  
  TList *pQvecList = dynamic_cast<TList*> (list->FindObject(subdetector));
  if (pQvecList != NULL) {
    /* the detector is present */
    if (TString(expectedstep).EqualTo("latest"))
      theQnVector = (AliQnCorrectionsQnVector*) pQvecList->First();
    else
      theQnVector = (AliQnCorrectionsQnVector*) pQvecList->FindObject(expectedstep);
    
    if (theQnVector == NULL || !(theQnVector->IsGoodQuality()) || !(theQnVector->GetN() != 0)) {
      /* the Qn vector for the expected step was not there */
      if (TString(altstep).EqualTo("latest"))
	theQnVector = (AliQnCorrectionsQnVector*) pQvecList->First();
      else
	theQnVector = (AliQnCorrectionsQnVector*) pQvecList->FindObject(altstep);
    }
  }
  if (theQnVector != NULL) {
    /* check the Qn vector quality */
    if (!(theQnVector->IsGoodQuality()) || !(theQnVector->GetN() != 0))
      /* not good quality, discarded */
      theQnVector = NULL;
  }
  return theQnVector;
}

//________________________________________________________________________
Double_t AliAnalysisTaskHFv1::Getq2(TList* qnlist, Int_t q2meth)
{
  if(!qnlist) {return -1;}
  
  const AliQnCorrectionsQnVector* qnVect = 0x0;
  
  if(q2meth==kq2TPC)
    qnVect = GetQnVectorFromList(qnlist, Form("%sQoverSqrtM",fDetTPCConfName[0].Data()), "latest", "plain");
  else if(q2meth==kq2NegTPC)
    qnVect = GetQnVectorFromList(qnlist, Form("%sQoverSqrtM",fDetTPCConfName[1].Data()), "latest", "plain");
  else if(q2meth==kq2PosTPC)
    qnVect = GetQnVectorFromList(qnlist, Form("%sQoverSqrtM",fDetTPCConfName[2].Data()), "latest", "plain");
  else if(q2meth==kq2VZERO)
    qnVect = GetQnVectorFromList(qnlist, Form("%sQoverSqrtM",fDetV0ConfName[0].Data()), "latest", "raw");
  else if(q2meth==kq2VZEROA)
    qnVect = GetQnVectorFromList(qnlist, Form("%sQoverSqrtM",fDetV0ConfName[1].Data()), "latest", "raw");
  else if(q2meth==kq2VZEROC)
    qnVect = GetQnVectorFromList(qnlist, Form("%sQoverSqrtM",fDetV0ConfName[2].Data()), "latest", "raw");
  else {return -1;}
  
  if(!qnVect) {return -1;}
  
  Double_t q2 = TMath::Sqrt(qnVect->Qx(2)*qnVect->Qx(2)+qnVect->Qy(2)*qnVect->Qy(2)); //qnVect->Length();
  return q2;
}

//________________________________________________________________________
void AliAnalysisTaskHFv1::Setq2Smearing(TString smearingfilepath, TString histoname, Int_t smearingaxis) {
  fq2Smearing=kTRUE;
  fq2SmearingAxis=smearingaxis;
  
  TFile* smearingfile = TFile::Open(smearingfilepath.Data(),"READ");
  fq2SmearingHisto=(TH2F*)smearingfile->Get(histoname.Data());
  if(fq2SmearingHisto) {fq2SmearingHisto->SetDirectory(0);}
  smearingfile->Close();
}

//________________________________________________________________________
Double_t AliAnalysisTaskHFv1::ComputeTPCq2(AliAODEvent* aod, Double_t &q2TPCfull, Double_t &q2TPCpos,Double_t &q2TPCneg) const {
  /// Compute the q2 for ESE starting from TPC tracks
  /// Option to reject a fraction of tracks to emulate resolution effects

  if(fq2Meth==kq2VZERO || fq2Meth==kq2VZEROA || fq2Meth==kq2VZEROC){
    AliWarning("The recalculation of q2 is implemented only for TPC\n");
    return 0.;
  }

  Int_t nTracks=aod->GetNumberOfTracks();
  Double_t qVec[2]={0.,0.};
  Double_t qVecPosEta[2]={0.,0.};
  Double_t qVecNegEta[2]={0.,0.};
  Double_t nHarmonic=2.;
  Int_t multQvec=0;
  Int_t multQvecPosEta=0;
  Int_t multQvecNegEta=0;
  for(Int_t it=0; it<nTracks; it++){
    AliAODTrack* track=(AliAODTrack*)aod->GetTrack(it);
    if(!track) continue;
    if(track->TestFilterBit(BIT(8))||track->TestFilterBit(BIT(9))) {
      Double_t pt=track->Pt();
      Double_t pseudoRand=pt*1000.-(Long_t)(pt*1000);
      Double_t eta=track->Eta();
      Double_t phi=track->Phi();
      Double_t qx=TMath::Cos(nHarmonic*phi);
      Double_t qy=TMath::Sin(nHarmonic*phi);
      if(pseudoRand<fFractionOfTracksForTPCq2){
	qVec[0]+=qx;
	qVec[1]+=qy;
	multQvec++;
      }
      if(eta>0){
	qVecPosEta[0]+=qx;
	qVecPosEta[1]+=qy;
	multQvecPosEta++;
      }else{
	qVecNegEta[0]+=qx;
	qVecNegEta[1]+=qy;
	multQvecNegEta++;
      }
    }
  }

  q2TPCfull = 0.;
  if(multQvec>0) q2TPCfull = TMath::Sqrt(qVec[0]*qVec[0]+qVec[1]*qVec[1])/TMath::Sqrt(multQvec);
  q2TPCpos = 0.;
  if(multQvecPosEta>0) q2TPCpos = TMath::Sqrt(qVecPosEta[0]*qVecPosEta[0]+qVecPosEta[1]*qVecPosEta[1])/TMath::Sqrt(multQvecPosEta);
  q2TPCneg = 0.;
  if(multQvecNegEta>0) q2TPCneg = TMath::Sqrt(qVecNegEta[0]*qVecNegEta[0]+qVecNegEta[1]*qVecNegEta[1])/TMath::Sqrt(multQvecNegEta);

  if(fq2Meth==kq2TPC) return q2TPCfull;
  else if(fq2Meth==kq2PosTPC) return q2TPCpos;
  else if(fq2Meth==kq2NegTPC) return q2TPCneg;
  else return 0.;
}


//________________________________________________________________________
void AliAnalysisTaskHFv1::Terminate(Option_t */*option*/)
{
  // Terminate analysis
  //
  if(fDebug > 1) printf("AnalysisTaskSEHFvn: Terminate() \n");
 
  return;
}

