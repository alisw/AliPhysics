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
// AliAnalysisTaskSEHFvn gives the needed tools for the D
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
#include <TH3F.h>
#include <TGraphErrors.h>
#include <TGraph.h>
#include <TDatabasePDG.h>
#include <TRandom3.h>
#include <TVector2.h>
#include <TArrayF.h>
#include <TAxis.h>
#include <TSpline.h>
#include <vector>

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
#include "AliVertexingHFUtils.h"

#include "AliAnalysisTaskSEHFvn.h"

ClassImp(AliAnalysisTaskSEHFvn)


//________________________________________________________________________
AliAnalysisTaskSEHFvn::AliAnalysisTaskSEHFvn():
AliAnalysisTaskSE(),
  fhEventsInfo(0),
  fHistCandVsCent(0x0),
  fHistCandMassRangeVsCent(0x0),
  fOutput(0),
  fRDCuts(0),
  fLowmasslimit(1.669),
  fUpmasslimit(2.069),
  fNPtBins(1),
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
  fScalProdLimit(0.3),
  fRemoveDauFromq2(0),
  fTPCEtaMin(-0.8),
  fTPCEtaMax(0.8),
  fRemoveNdauRandomTracks(kFALSE),
  fUseQnFrameworkCorrq2(kTRUE),
  fRequireMassForDauRemFromq2(kFALSE),
  fDeltaEtaDmesonq2(0.),
  fEPVsq2VsCent(kFALSE),
  fEnableNtrklHistos(kFALSE),
  fRemoverSoftPionFromq2(kFALSE),
  fPercentileq2(kFALSE),
  fEnableCentralityCorrCuts(kFALSE),
  fEnableCentralityMultiplicityCorrStrongCuts(kFALSE),
  fOptD0FromDstar(0),
  fUseFiltBit4SoftPion(kFALSE),
  fCutsSoftPion(0x0),
  fFlowMethod(kEP)
{
  // Default constructor
  for(int i = 0; i < 3; i++) {
    fHistCentrality[i]     = 0x0;
    fHistEvPlaneQncorrTPC[i]   = 0x0;
    fHistEvPlaneQncorrVZERO[i] = 0x0;
  }
  
  for(int i = 0; i < 6; i++) {
    fq2SplinesList[i] = 0x0;
  }
  
  if(fFlowMethod==kEvShape) {
    fCentBinSizePerMil=10.;
  }
}

//________________________________________________________________________
AliAnalysisTaskSEHFvn::AliAnalysisTaskSEHFvn(const char *name,AliRDHFCuts *rdCuts,Int_t decaychannel):
  AliAnalysisTaskSE(name),
  fhEventsInfo(0),
  fHistCandVsCent(0x0),
  fHistCandMassRangeVsCent(0x0),
  fOutput(0),
  fRDCuts(rdCuts),
  fLowmasslimit(0),
  fUpmasslimit(0),
  fNPtBins(1),
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
  fScalProdLimit(0.3),
  fRemoveDauFromq2(0),
  fTPCEtaMin(-0.8),
  fTPCEtaMax(0.8),
  fRemoveNdauRandomTracks(kFALSE),
  fUseQnFrameworkCorrq2(kTRUE),
  fRequireMassForDauRemFromq2(kFALSE),
  fDeltaEtaDmesonq2(0.),
  fEPVsq2VsCent(kFALSE),
  fEnableNtrklHistos(kFALSE),
  fRemoverSoftPionFromq2(kFALSE),
  fPercentileq2(kFALSE),
  fEnableCentralityCorrCuts(kFALSE),
  fEnableCentralityMultiplicityCorrStrongCuts(kFALSE),
  fOptD0FromDstar(0),
  fUseFiltBit4SoftPion(kFALSE),
  fCutsSoftPion(0x0),
  fFlowMethod(kEP)
{
  // standard constructor
  for(int i = 0; i < 3; i++) {
    fHistCentrality[i]     = 0x0;
    fHistEvPlaneQncorrTPC[i]   = 0x0;
    fHistEvPlaneQncorrVZERO[i] = 0x0;
  }

  for(int i = 0; i < 6; i++) {
    fq2SplinesList[i] = 0x0;
  }

  if(fFlowMethod==kEvShape) {
    fCentBinSizePerMil=10.;
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
  case 4:
    DefineOutput(3,AliRDHFCutsD0toKpi::Class());  //Cut object for D0
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
AliAnalysisTaskSEHFvn::~AliAnalysisTaskSEHFvn()
{
  // Destructor
  if(fOutput && !fOutput->IsOwner()){
    for(Int_t i=0;i<3;i++){
      delete fHistCentrality[i];
      delete fHistEvPlaneQncorrTPC[i];
      delete fHistEvPlaneQncorrVZERO[i];
    }
    delete fHistCandVsCent;
    delete fHistCandMassRangeVsCent;
  }
  delete fOutput;
  delete fhEventsInfo;
  delete fRDCuts;
  delete fAfterBurner;
  if(fq2SmearingHisto) {delete fq2SmearingHisto;}
  for(Int_t i=0; i<6; i++) {
    if(fq2SplinesList[i]) {delete fq2SplinesList[i];}
  }
  if(fDecChannel==kD0toKpiFromDstar && fOptD0FromDstar==1 && !fUseFiltBit4SoftPion) {
    delete fCutsSoftPion;
  }
}
//_________________________________________________________________
void  AliAnalysisTaskSEHFvn::SetMassLimits(Float_t range, Int_t pdg){
  // Set limits for mass spectra plots
  Float_t mass=0;
  Int_t abspdg=TMath::Abs(pdg);
  mass=TDatabasePDG::Instance()->GetParticle(abspdg)->Mass();
  fUpmasslimit = mass+range;
  fLowmasslimit = mass-range;
}
//_________________________________________________________________
void  AliAnalysisTaskSEHFvn::SetMassLimits(Float_t lowlimit, Float_t uplimit){
  // Set limits for mass spectra plots
  if(uplimit>lowlimit)
    {
      fUpmasslimit = uplimit;
      fLowmasslimit = lowlimit;
    }
}
//________________________________________________________________________
void AliAnalysisTaskSEHFvn::LocalInit()
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
    case 4:
    {
      AliRDHFCutsD0toKpi* copycut=new AliRDHFCutsD0toKpi(*(static_cast<AliRDHFCutsD0toKpi*>(fRDCuts)));
      // Post the data
      PostData(3,copycut);
    }
      break;
  default:
    return;
  }
  
  if(fDecChannel==kD0toKpiFromDstar && fOptD0FromDstar==1 && !fUseFiltBit4SoftPion) {
    fCutsSoftPion = new AliESDtrackCuts();
    fCutsSoftPion->SetRequireSigmaToVertex(kFALSE);
    //default
    fCutsSoftPion->SetRequireTPCRefit(kFALSE);
    fCutsSoftPion->SetRequireITSRefit(kTRUE);
    fCutsSoftPion->SetMinNClustersITS(3);
    fCutsSoftPion->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
    fCutsSoftPion->SetPtRange(0.1,1.e10);
    fCutsSoftPion->SetEtaRange(-0.8,0.8);
  }

  return;
}
//________________________________________________________________________
void AliAnalysisTaskSEHFvn::UserCreateOutputObjects()
{
  // Create the output container

  if(fDebug > 1) printf("AnalysisTaskSEHFvn::UserCreateOutputObjects() \n");

  fhEventsInfo = new TH1F(GetOutputSlot(1)->GetContainer()->GetName(), "Number of AODs scanned",16,-0.5,15.5);
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
  fhEventsInfo->GetXaxis()->SetBinLabel(12,"n. rejected for bad cent corr");
  fhEventsInfo->GetXaxis()->SetBinLabel(13,"non valid TPC EP");
  fhEventsInfo->GetXaxis()->SetBinLabel(14,"bad event plane");
  fhEventsInfo->GetXaxis()->SetBinLabel(15,"no. of sel. candidates");
  fhEventsInfo->GetXaxis()->SetBinLabel(16,"no. cand. out of pt bounds");

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

  if(fEnableCentralityCorrCuts) {
    fEventCuts.AddQAplotsToList(fOutput,true);
  }
  
  const Int_t ncentbins = (fMaxCentr-fMinCentr)*10/fCentBinSizePerMil;
  
  fHistCandVsCent=new TH2F("hCandVsCent","number of selected candidates vs. centrality;centrality(%);number of candidates",ncentbins,fMinCentr,fMaxCentr,101,-0.5,100.5);
  fOutput->Add(fHistCandVsCent);
  fHistCandMassRangeVsCent=new TH2F("hCandMassRangeVsCent","number of selected candidates vs. centrality;centrality(%);number of candidates",ncentbins,fMinCentr,fMaxCentr,101,-0.5,100.5);
  fOutput->Add(fHistCandMassRangeVsCent);

  for(int iDet = 0; iDet < 3; iDet++) {
    fHistEvPlaneQncorrTPC[iDet]   = new TH1F(Form("hEvPlaneQncorr%s%s",fDetTPCConfName[iDet].Data(),fNormMethod.Data()),Form("hEvPlaneQncorr%s%s;#phi Ev Plane;Entries",fDetTPCConfName[iDet].Data(),fNormMethod.Data()),200,0.,TMath::Pi());
    fHistEvPlaneQncorrVZERO[iDet] = new TH1F(Form("hEvPlaneQncorr%s%s",fDetV0ConfName[iDet].Data(),fNormMethod.Data()),Form("hEvPlaneQncorr%s%s;#phi Ev Plane;Entries",fDetV0ConfName[iDet].Data(),fNormMethod.Data()),200,0.,TMath::Pi());
    fOutput->Add(fHistEvPlaneQncorrTPC[iDet]);
    fOutput->Add(fHistEvPlaneQncorrVZERO[iDet]);
  }
  if(fFlowMethod==kSP) {
    TH1F* hNormQ1 = new TH1F("hNormQ1","hNormQ1;|Q_{1}|;Entries",100,0.,1);
    fOutput->Add(hNormQ1);
    TH1F* hNormQ2 = new TH1F("hNormQ2","hNormQ2;|Q_{2}|;Entries",100,0.,1);
    fOutput->Add(hNormQ2);
    TH1F* hNormQ3 = new TH1F("hNormQ3","hNormQ3;|Q_{3}|;Entries",100,0.,1);
    fOutput->Add(hNormQ3);
  }

  TString q2axisname="q_{2}^{TPC}";
  if(fq2Meth==kq2PosTPC) {q2axisname="q_{2}^{pos TPC}";}
  else if(fq2Meth==kq2NegTPC) {q2axisname="q_{2}^{neg TPC}";}
  else if(fq2Meth==kq2VZERO) {q2axisname="q_{2}^{V0}";}
  else if(fq2Meth==kq2VZEROA) {q2axisname="q_{2}^{V0A}";}
  else if(fq2Meth==kq2VZEROC) {q2axisname="q_{2}^{V0C}";}
  
  TString q2percaxisname = q2axisname + " (%)";
  TString q2axisnamefill = q2axisname;
  if(fPercentileq2) {q2axisnamefill=q2percaxisname;}
  
  Int_t nq2bins=500;
  Double_t q2min = 0.;
  Double_t q2max = 10.;
  if(fPercentileq2) {
    nq2bins=100;
    q2min = 0.;
    q2max = 100.;
  }

  if(fFlowMethod==kEvShape) {
    //EP angle vs q2
    TH3F* hEvPlaneQncorrTPCVsq2VsCent[3];
    TH3F* hEvPlaneQncorrVZEROVsq2VsCent[3];
    if(fEPVsq2VsCent) {
      for(Int_t iDet=0; iDet<3; iDet++) {
        hEvPlaneQncorrTPCVsq2VsCent[iDet] = new TH3F(Form("hEvPlaneQncorr%s%sVsq2VsCent",fDetTPCConfName[iDet].Data(),fNormMethod.Data()),Form("hEvPlaneQncorr%s%sVsq2VsCent;centrality(%%);%s;#phi Ev Plane",fDetTPCConfName[iDet].Data(),fNormMethod.Data(),q2axisnamefill.Data()),ncentbins,fMinCentr,fMaxCentr,nq2bins,q2min,q2max,100,0.,TMath::Pi());
        hEvPlaneQncorrVZEROVsq2VsCent[iDet] = new TH3F(Form("hEvPlaneQncorr%s%sVsq2VsCent",fDetV0ConfName[iDet].Data(),fNormMethod.Data()),Form("hEvPlaneQncorr%s%sVsq2VsCent;centrality(%%);%s;#phi Ev Plane",fDetV0ConfName[iDet].Data(),fNormMethod.Data(),q2axisnamefill.Data()),ncentbins,fMinCentr,fMaxCentr,nq2bins,q2min,q2max,100,0.,TMath::Pi());
        fOutput->Add(hEvPlaneQncorrTPCVsq2VsCent[iDet]);
        fOutput->Add(hEvPlaneQncorrVZEROVsq2VsCent[iDet]);
      }
    }

    // histos for q2 vs. centrality with fine binning (for q2 percentiles calibration)
    for(Int_t iDet=0; iDet<3; iDet++) {
      TH2F* hq2vsCentrTPC=new TH2F(Form("hq2vsCentr%s",fDetTPCConfName[iDet].Data()),Form("q_{2}^{%s} vs. centrality;centrality(%%);q_{2}^{%s}",fDetTPCConfName[iDet].Data(),fDetTPCConfName[iDet].Data()),ncentbins,fMinCentr,fMaxCentr,10000,0.,15.);
      fOutput->Add(hq2vsCentrTPC);
      TH2F* hq2vsCentrV0=new TH2F(Form("hq2vsCentr%s",fDetV0ConfName[iDet].Data()),Form("q_{2}^{%s} vs. centrality;centrality(%%);q_{2}^{%s}",fDetV0ConfName[iDet].Data(),fDetV0ConfName[iDet].Data()),ncentbins,fMinCentr,fMaxCentr,10000,0.,15.);
      fOutput->Add(hq2vsCentrV0);
    }
    
    TH3F* hPercq2vsq2vsCentr = new TH3F("hPercq2vsq2vsCentr",Form("%s vs. %s vs. centrality;centrality (%%);%s;%s",q2percaxisname.Data(),q2axisname.Data(),q2axisname.Data(),q2percaxisname.Data()),ncentbins,fMinCentr,fMaxCentr,600,0.,12.,100,0.,100.);
    fOutput->Add(hPercq2vsq2vsCentr);
    
    //multiplicity used for q2 vs. centrality (TPC)
    TH2F* hMultVsCentFullTPC = new TH2F("hMultVsCentFullTPC","Multiplicity for q_{2} vs. centrality (full TPC);centrality(%);M",ncentbins,fMinCentr,fMaxCentr,100,-0.5,4999.5);
    TH2F* hMultVsCentPosTPC = new TH2F("hMultVsCentPosTPC","Multiplicity for q_{2} vs. centrality (pos TPC);centrality(%);M",ncentbins,fMinCentr,fMaxCentr,100,-0.5,4999.5);
    TH2F* hMultVsCentNegTPC = new TH2F("hMultVsCentNegTPC","Multiplicity for q_{2} vs. centrality (neg TPC);centrality(%);M",ncentbins,fMinCentr,fMaxCentr,100,-0.5,4999.5);
    fOutput->Add(hMultVsCentFullTPC);
    fOutput->Add(hMultVsCentPosTPC);
    fOutput->Add(hMultVsCentNegTPC);

    //q2TPC correlations
    TH2F* hq2TPCPosEtaVsNegEta = new TH2F("hq2TPCPosEtaVsNegEta","q_{2}^{pos TPC} vs. q_{2}^{neg TPC};q_{2}^{neg TPC};q_{2}^{pos TPC}",500,0.,10.,500,0.,10.);
    TH2F* hq2TPCFullEtaVsNegEta = new TH2F("hq2TPCFullEtaVsNegEta","q_{2}^{full TPC} vs. q_{2}^{neg TPC};q_{2}^{neg TPC};q_{2}^{full TPC}",500,0.,10.,500,0.,10.);
    TH2F* hq2TPCFullEtaVsPosEta = new TH2F("hq2TPCFullEtaVsPosEta","q_{2}^{full TPC} vs. q_{2}^{pos TPC};q_{2}^{pos TPC};q_{2}^{full TPC}",500,0.,10.,500,0.,10.);
    fOutput->Add(hq2TPCPosEtaVsNegEta);
    fOutput->Add(hq2TPCFullEtaVsNegEta);
    fOutput->Add(hq2TPCFullEtaVsPosEta);
    
    //correlation between q2TPC and q2VZEROC/A (all events, events with candidates, events with candidates in mass range)
    TH2F* hq2TPCVsq2VZEROC = new TH2F("hq2TPCVsq2VZEROC","q_{2}^{V0C} vs. q_{2}^{TPC};q_{2}^{TPC};q_{2}^{V0C}",nq2bins,q2min,q2max,nq2bins,q2min,q2max);
    TH2F* hq2TPCVsq2VZEROCCand = new TH2F("hq2TPCVsq2VZEROCCand","q_{2}^{V0C} vs. q_{2}^{TPC};q_{2}^{TPC};q_{2}^{V0C}",nq2bins,q2min,q2max,nq2bins,q2min,q2max);
    TH2F* hq2TPCVsq2VZEROCCandInMass = new TH2F("hq2TPCVsq2VZEROCCandInMass","q_{2}^{V0C} vs. q_{2}^{TPC};q_{2}^{TPC};q_{2}^{V0C}",nq2bins,q2min,q2max,nq2bins,q2min,q2max);
    TH2F* hq2TPCVsq2VZEROA = new TH2F("hq2TPCVsq2VZEROA","q_{2}^{V0A} vs. q_{2}^{TPC};q_{2}^{TPC};q_{2}^{V0A}",nq2bins,q2min,q2max,nq2bins,q2min,q2max);
    TH2F* hq2TPCVsq2VZEROACand = new TH2F("hq2TPCVsq2VZEROACand","q_{2}^{V0A} vs. q_{2}^{TPC};q_{2}^{TPC};q_{2}^{V0A}",nq2bins,q2min,q2max,nq2bins,q2min,q2max);
    TH2F* hq2TPCVsq2VZEROACandInMass = new TH2F("hq2TPCVsq2VZEROACandInMass","q_{2}^{V0A} vs. q_{2}^{TPC};q_{2}^{TPC};q_{2}^{V0A}",nq2bins,q2min,q2max,nq2bins,q2min,q2max);

    fOutput->Add(hq2TPCVsq2VZEROC);
    fOutput->Add(hq2TPCVsq2VZEROCCand);
    fOutput->Add(hq2TPCVsq2VZEROCCandInMass);
    fOutput->Add(hq2TPCVsq2VZEROA);
    fOutput->Add(hq2TPCVsq2VZEROACand);
    fOutput->Add(hq2TPCVsq2VZEROACandInMass);

    //Ntracklets vs. q2 vs. centrality histos
    if(fEnableNtrklHistos) {
      TH3F* hNtrklVsq2VsCent = new TH3F("hNtrklVsq2VsCent",Form("N_{tracklets} vs. %s vs. centrality;centrality(%%);%s;N_{tracklets}",q2axisnamefill.Data(),q2axisnamefill.Data()),ncentbins,fMinCentr,fMaxCentr,nq2bins,q2min,q2max,500,-0.5,4999.5);
      TH3F* hNtrklVsq2VsCentCand = new TH3F("hNtrklVsq2VsCentCand",Form("N_{tracklets} vs. %s vs. centrality (cand);centrality(%%);%s;N_{tracklets}",q2axisnamefill.Data(),q2axisnamefill.Data()),ncentbins,fMinCentr,fMaxCentr,nq2bins,q2min,q2max,500,-0.5,4999.5);
      TH3F* hNtrklVsq2VsCentCandInMass = new TH3F("hNtrklVsq2VsCentCandInMass",Form("N_{tracklets} vs. %s vs. centrality (cand in mass);centrality(%%);%s;N_{tracklets}",q2axisnamefill.Data(),q2axisnamefill.Data()),ncentbins,fMinCentr,fMaxCentr,nq2bins,q2min,q2max,500,-0.5,4999.5);
      fOutput->Add(hNtrklVsq2VsCent);
      fOutput->Add(hNtrklVsq2VsCentCand);
      fOutput->Add(hNtrklVsq2VsCentCandInMass);
    }
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
          if((fDecChannel != AliAnalysisTaskSEHFvn::kDplustoKpipi) &&(fDecChannel != AliAnalysisTaskSEHFvn::kDstartoKpipi)){
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
      TString SPsuffix3="";
      TString SPsuffix4="";
      TString SPsuffix5="";
      if(fEvPlaneDet==kFullTPC || fEvPlaneDet==kPosTPC || fEvPlaneDet==kNegTPC) { //two sub-event formula with two TPC halves
        SPsuffix1="u1Q2"; //1->neg TPC, 2->pos TPC
        SPsuffix2="u2Q1";
        SPsuffix3="Q1Q2";
      }
      else if(fEvPlaneDet==kFullV0){ // three sub-event formula (two terms at numerator with V0A and V0C)
        SPsuffix1="u3Q1"; //3->full TPC, 1->V0A, 2->V0C
        SPsuffix2="u3Q2";
        SPsuffix3="Q1Q2";
      }
      else if(fEvPlaneDet==kV0A || fEvPlaneDet==kV0C){ //three sub-event formula (one term at numerator)
        SPsuffix1="u3Q1"; //3->full TPC, 1->main V0 (V0A or V0C), 2->V0 for resolution term (V0C or V0A)
        SPsuffix3="Q1Q2";
        SPsuffix4="Q1Q3";
        SPsuffix5="Q2Q3";
      }
      for(Int_t i=0;i<fNPtBins;i++){
        TH2F* hMassScalProd1=new TH2F(Form("hMassScalProd%s_pt%d%s",SPsuffix1.Data(),i,centrname.Data()),Form("Mass vs hScalProd%s (p_{t} bin %d %s);%s;M (GeV/c^{2})",SPsuffix1.Data(),i,centrname.Data(),SPsuffix1.Data()),100,-fScalProdLimit,fScalProdLimit,fNMassBins,fLowmasslimit,fUpmasslimit);
        fOutput->Add(hMassScalProd1);
        if(fEvPlaneDet==kFullTPC || fEvPlaneDet==kPosTPC || fEvPlaneDet==kNegTPC || fEvPlaneDet==kFullV0) {
          TH2F* hMassScalProd2=new TH2F(Form("hMassScalProd%s_pt%d%s",SPsuffix2.Data(),i,centrname.Data()),Form("Mass vs hScalProd%s (p_{t} bin %d %s);%s;M (GeV/c^{2})",SPsuffix2.Data(),i,centrname.Data(),SPsuffix2.Data()),100,-fScalProdLimit,fScalProdLimit,fNMassBins,fLowmasslimit,fUpmasslimit);
          fOutput->Add(hMassScalProd2);
        }
      }
      TH1F* hScalProdDenom1=new TH1F(Form("hScalProd%s_%s",SPsuffix3.Data(),centrname.Data()),Form("hScalProd%s (%s);%s",SPsuffix3.Data(),centrname.Data(),SPsuffix3.Data()),100,-fScalProdLimit*fScalProdLimit,fScalProdLimit*fScalProdLimit);
      fOutput->Add(hScalProdDenom1);
      if(fEvPlaneDet==kV0A || fEvPlaneDet==kV0C) {
        TH1F* hScalProdDenom2=new TH1F(Form("hScalProd%s_%s",SPsuffix4.Data(),centrname.Data()),Form("hScalProd%s (%s);%s",SPsuffix4.Data(),centrname.Data(),SPsuffix4.Data()),100,-fScalProdLimit*fScalProdLimit,fScalProdLimit*fScalProdLimit);
        fOutput->Add(hScalProdDenom2);
        TH1F* hScalProdDenom3=new TH1F(Form("hScalProd%s_%s",SPsuffix5.Data(),centrname.Data()),Form("hScalProd%s (%s);%s",SPsuffix5.Data(),centrname.Data(),SPsuffix5.Data()),100,-fScalProdLimit*fScalProdLimit,fScalProdLimit*fScalProdLimit);
        fOutput->Add(hScalProdDenom3);
      }
    }
    else if (fFlowMethod==kEvShape){
      // histos for EP resolution vs q2
      TH2F* hEvPlaneResoVsq2=new TH2F(Form("hEvPlaneReso1Vsq2%s",centrname.Data()),Form("Event plane angle Resolution vs. %s %s;cos2(#psi_{A}-#psi_{B});%s;Entries",q2axisnamefill.Data(),centrname.Data(),q2axisnamefill.Data()),220,-1.1,1.1,nq2bins,q2min,q2max);
      fOutput->Add(hEvPlaneResoVsq2);
      TH2F* hEvPlaneReso2Vsq2=new TH2F(Form("hEvPlaneReso2Vsq2%s",centrname.Data()),Form("Event plane angle Resolution vs. %s %s;cos2(#psi_{A}-#psi_{C});%s;Entries",q2axisnamefill.Data(),centrname.Data(),q2axisnamefill.Data()),220,-1.1,1.1,nq2bins,q2min,q2max);
      fOutput->Add(hEvPlaneReso2Vsq2);
      TH2F* hEvPlaneReso3Vsq2=new TH2F(Form("hEvPlaneReso3Vsq2%s",centrname.Data()),Form("Event plane angle Resolution vs. %s %s;cos2(#psi_{B}-#psi_{C});%s;Entries",q2axisnamefill.Data(),centrname.Data(),q2axisnamefill.Data()),220,-1.1,1.1,nq2bins,q2min,q2max);
      fOutput->Add(hEvPlaneReso3Vsq2);
    }
  }

  if(fFlowMethod==kEvShape) {
    //q2 candidate vs. q2 global event
    TH2F* hq2CandVsq2Event=new TH2F("hq2CandVsq2Event",Form("%s candidate vs. %s global event;%s global event;%s candidate;Entries",q2axisname.Data(),q2axisname.Data(),q2axisname.Data(),q2axisname.Data()),nq2bins,q2min,q2max,nq2bins,q2min,q2max);
    fOutput->Add(hq2CandVsq2Event);

    //histo of number of Dstar built using from the same Dzero (only in case of Dzero<-Dstar analysis)
    if(fDecChannel==kD0toKpiFromDstar && fOptD0FromDstar==0) {
      TH1F* hNumDstarCandFromSameDzeroCand = new TH1F("hNumDstarCandFromSameDzeroCand","Number of Dstar candidates built using same Dzero candidate",101,-0.5,100.5);
      fOutput->Add(hNumDstarCandFromSameDzeroCand);
    }
    
    CreateSparseForEvShapeAnalysis();
  }

  PostData(1,fhEventsInfo);
  PostData(2,fOutput);

  return;
}

//________________________________________________________________________
void AliAnalysisTaskSEHFvn::UserExec(Option_t */*option*/)
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
      if(fDecChannel==4){
        absPdgMom=421;
        if(fOptD0FromDstar==0) {
          arrayProng=(TClonesArray*)aodFromExt->GetList()->FindObject("Dstar");
        }
        else {
          arrayProng=(TClonesArray*)aodFromExt->GetList()->FindObject("D0toKpi");
        }
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
    if(fDecChannel==4){
      absPdgMom=421;
      if(fOptD0FromDstar==0) {
        arrayProng=(TClonesArray*)aod->GetList()->FindObject("Dstar");
      }
      else {
        arrayProng=(TClonesArray*)aod->GetList()->FindObject("D0toKpi");
      }
    }
  }

  if(!aod || !arrayProng) {
    AliError("AliAnalysisTaskSEHFvn::UserExec:Branch not found!\n");
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
      AliWarning("AliAnalysisTaskSEHFvn::UserExec:MC particles branch not found!\n");
      return;
    }

    // load MC header
    mcHeader =  (AliAODMCHeader*)aod->GetList()->FindObject(AliAODMCHeader::StdBranchName());
    if(!mcHeader) {
      AliError("AliAnalysisTaskSEHFvn::UserExec:MC header branch not found!\n");
      return;
    }
  }

  AliAODRecoDecayHF *d=0;

  Int_t nCand = arrayProng->GetEntriesFast();
  if(fDebug>2) printf("Number of D2H: %d\n",nCand);

  Double_t evCentr = fRDCuts->GetCentrality(aod);
  Bool_t isEvSel=fRDCuts->IsEventSelected(aod);
  
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

  fEventCuts.AcceptEvent(aod);
  if(!fEventCuts.PassedCut(AliEventCuts::kCorrelations) || !fEventCuts.PassedCut(AliEventCuts::kMultiplicity)) {
    fhEventsInfo->Fill(11);
    return;
  }

  fhEventsInfo->Fill(4);
  fHistCentrality[1]->Fill(evCentr);

  Double_t eventplaneqncorrTPC[3];
  Double_t eventplaneqncorrVZERO[3];
  TList *qnlist = 0x0;
  Double_t Q1[2]={-999.,-999.};
  Double_t Q2[2]={-999.,-999.};
  Double_t Q3[2]={-999.,-999.};

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
      Q1[0]=qnVectTPC[1]->Qx(fHarmonic);//NegTPC
      Q1[1]=qnVectTPC[1]->Qy(fHarmonic);//NegTPC
      Q2[0]=qnVectTPC[2]->Qx(fHarmonic);//PosTPC
      Q2[1]=qnVectTPC[2]->Qy(fHarmonic);//PosTPC
    }
    else if(fEvPlaneDet==kFullV0) {
      Q1[0]=qnVectV0[1]->Qx(fHarmonic);//V0A
      Q1[1]=qnVectV0[1]->Qy(fHarmonic);//V0A
      Q2[0]=qnVectV0[2]->Qx(fHarmonic);//V0C
      Q2[1]=qnVectV0[2]->Qy(fHarmonic);//V0C
    }
    else if(fEvPlaneDet==kV0A) {
      Q1[0]=qnVectV0[1]->Qx(fHarmonic);//V0A
      Q1[1]=qnVectV0[1]->Qy(fHarmonic);//V0A
      Q2[0]=qnVectV0[2]->Qx(fHarmonic);//V0C
      Q2[1]=qnVectV0[2]->Qy(fHarmonic);//V0C
      Q3[0]=qnVectTPC[0]->Qx(fHarmonic);//FullTPC
      Q3[1]=qnVectTPC[0]->Qy(fHarmonic);//FullTPC
    }
    else if(fEvPlaneDet==kV0C) {
      Q1[0]=qnVectV0[2]->Qx(fHarmonic);//V0C
      Q1[1]=qnVectV0[2]->Qy(fHarmonic);//V0C
      Q2[0]=qnVectV0[1]->Qx(fHarmonic);//V0A
      Q2[1]=qnVectV0[1]->Qy(fHarmonic);//V0A
      Q3[0]=qnVectTPC[0]->Qx(fHarmonic);//FullTPC
      Q3[1]=qnVectTPC[0]->Qy(fHarmonic);//FullTPC
    }
    if(fFlowMethod==kSP) {
      ((TH1F*)fOutput->FindObject("hNormQ1"))->Fill(TMath::Sqrt(Q1[0]*Q1[0]+Q1[1]*Q1[1]));
      ((TH1F*)fOutput->FindObject("hNormQ2"))->Fill(TMath::Sqrt(Q2[0]*Q2[0]+Q2[1]*Q2[1]));
      if(fEvPlaneDet==kV0A || fEvPlaneDet==kV0C) {((TH1F*)fOutput->FindObject("hNormQ3"))->Fill(TMath::Sqrt(Q3[0]*Q3[0]+Q3[1]*Q3[1]));}
    }
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

  //get q2 for event shape anaylsis
  Double_t q2=-1;
  Double_t q2PosTPC=-1;
  Double_t q2NegTPC=-1;
  Double_t q2FullTPC=-1;
  Double_t q2VZERO=-1;
  Double_t q2VZEROA=-1;
  Double_t q2VZEROC=-1;
  //keep q vector and multiplicity for daughter removals (if activated)
  Double_t qVecFullTPC[2] = {0.,0.};
  Double_t qVecPosTPC[2] = {0.,0.};
  Double_t qVecNegTPC[2] = {0.,0.};
  Double_t multQvecTPC[3] = {0.,0.,0.}; //{full TPC, pos TPC, neg TPC}
  Double_t qVecDefault[2]={0.,0.};
  Double_t multQvecDefault=0;
  std::vector<Int_t> labrejtracks;
  if(fFlowMethod==kEvShape) {
    if(!fUseNewQnCorrFw && !fOnTheFlyTPCq2){
      AliWarning("AliAnalysisTaskSEHFvn: if you do not want to use the new Qn-correction framework for kEvShape, set q2 on-the-fly!\n");
      return;
    }
    if(fOnTheFlyTPCq2){
      q2=ComputeTPCq2(aod,q2FullTPC,q2PosTPC,q2NegTPC,qVecFullTPC,qVecPosTPC,qVecNegTPC,multQvecTPC,labrejtracks);
    }else{
      q2 = Getq2(qnlist,fq2Meth,multQvecDefault);
      q2PosTPC = Getq2(qnlist,kq2PosTPC,multQvecTPC[1]);
      q2NegTPC = Getq2(qnlist,kq2NegTPC,multQvecTPC[2]);
      q2FullTPC = Getq2(qnlist,kq2TPC,multQvecTPC[0]);
    }
    Double_t multV0=-1;
    q2VZERO = Getq2(qnlist,kq2VZERO,multV0);
    q2VZEROA = Getq2(qnlist,kq2VZEROA,multV0);
    q2VZEROC = Getq2(qnlist,kq2VZEROC,multV0);

    if(q2<0 || q2PosTPC<0 || q2NegTPC<0 || q2FullTPC<0 || q2VZERO<0 || q2VZEROA<0 || q2VZEROC<0) return;
  }

  Float_t planereso=0.;
  Float_t deltaSubAC=0.;
  Float_t deltaSubBC=0.;
  Int_t nSubEvents=2;

  AliEventplane *pl=aod->GetEventplane();
  if(!pl){
    Printf("AliAnalysisTaskSEHFvn::UserExec:no eventplane! v2 analysis without eventplane not possible!\n");
    fhEventsInfo->Fill(12);
    return;
  }

  if(fReadMC){
    TRandom3 *g = new TRandom3(0);
    eventplane=g->Rndm()*TMath::Pi();
    delete g;g=0x0;
    if(fUseAfterBurner)fAfterBurner->SetEventPlane((Double_t)eventplane);
  }else{
    if(fFlowMethod!=kSP) {
      eventplane=GetEventPlane(aod,pl,eventplaneqncorrTPC,eventplaneqncorrVZERO,planereso,deltaSubAC,deltaSubBC,nSubEvents);
      if(eventplane<-999){
        Printf("Bad event plane calculation\n");
        fhEventsInfo->Fill(13);
        return;
      }
    }
    else {
      //fill histograms for denominator in SP formula
      Double_t scalprod=-999.;
      if(fEvPlaneDet==kFullTPC || fEvPlaneDet==kPosTPC || fEvPlaneDet==kNegTPC || fEvPlaneDet==kFullV0) {
        scalprod=Q1[0]*Q2[0]+Q1[1]*Q2[1];
        ((TH1F*)fOutput->FindObject(Form("hScalProdQ1Q2_%s",fCentrBinName.Data())))->Fill(scalprod);
      }
      else if(fEvPlaneDet==kV0A || fEvPlaneDet==kV0C) {
        scalprod=Q1[0]*Q2[0]+Q1[1]*Q2[1];
        ((TH1F*)fOutput->FindObject(Form("hScalProdQ1Q2_%s",fCentrBinName.Data())))->Fill(scalprod);
        scalprod=Q1[0]*Q3[0]+Q1[1]*Q3[1];
        ((TH1F*)fOutput->FindObject(Form("hScalProdQ1Q3_%s",fCentrBinName.Data())))->Fill(scalprod);
        scalprod=Q2[0]*Q3[0]+Q2[1]*Q3[1];
        ((TH1F*)fOutput->FindObject(Form("hScalProdQ2Q3_%s",fCentrBinName.Data())))->Fill(scalprod);
      }
    }
  }
  if(fDebug>2)printf("Loop on D candidates\n");

  //quantities used to remove daughter tracks or tracks close in eta from q2 on-the-fly
  Double_t qVecRemDauFullTPC[2]={qVecFullTPC[0],qVecFullTPC[1]};
  Double_t qVecRemDauPosTPC[2]={qVecPosTPC[0],qVecPosTPC[1]};
  Double_t qVecRemDauNegTPC[2]={qVecNegTPC[0],qVecNegTPC[1]};
  Double_t multQvecRemDauTPC[3]={multQvecTPC[0],multQvecTPC[1],multQvecTPC[2]};

  //quantities used to remove daughter tracks from q2 with Qn-framework
  TList *pQvecList[3] = {0x0,0x0,0x0};
  AliQnCorrectionsQnVector *theQnVectorCorr[3]   = {0x0,0x0,0x0};
  AliQnCorrectionsQnVector *theQnVectorUncorr[3] = {0x0,0x0,0x0};
  Double_t dauqFullTPC[2] = {0.,0.};
  Double_t dauqPosTPC[2] = {0.,0.};
  Double_t dauqNegTPC[2] = {0.,0.};
  Double_t corrFullTPC[2] = {0.,0.};
  Double_t corrPosTPC[2] = {0.,0.};
  Double_t corrNegTPC[2] = {0.,0.};
  Double_t LbFullTPC[2] = {0.,0.};
  Double_t LbPosTPC[2] = {0.,0.};
  Double_t LbNegTPC[2] = {0.,0.};
  Double_t qRecFullTPC[2] = {0.,0.};
  Double_t qRecPosTPC[2] = {0.,0.};
  Double_t qRecNegTPC[2] = {0.,0.};
  Double_t qTwistFullTPC[2] = {0.,0.};
  Double_t qTwistPosTPC[2] = {0.,0.};
  Double_t qTwistNegTPC[2] = {0.,0.};
  Bool_t isTwistApplied[3] = {kFALSE,kFALSE,kFALSE};
  Double_t nDauRemoved[3]={0.,0.,0.};

  if((fq2Meth==kq2TPC || fq2Meth==kq2PosTPC || fq2Meth==kq2NegTPC) && fRemoveDauFromq2>0 && !fOnTheFlyTPCq2) {
    pQvecList[0] = dynamic_cast<TList*> (qnlist->FindObject(Form("%sQoverSqrtM",fDetTPCConfName[0].Data()))); //full TPC
    pQvecList[1] = dynamic_cast<TList*> (qnlist->FindObject(Form("%sQoverSqrtM",fDetTPCConfName[2].Data()))); //pos TPC
    pQvecList[2] = dynamic_cast<TList*> (qnlist->FindObject(Form("%sQoverSqrtM",fDetTPCConfName[1].Data()))); //neg TPC

    for(Int_t iDet=0; iDet<3; iDet++) {
      /* the detector is present */
      theQnVectorUncorr[iDet] = (AliQnCorrectionsQnVector*) pQvecList[iDet]->FindObject("plain"); //raw step for TPC
      if (theQnVectorUncorr[iDet] && theQnVectorUncorr[iDet]->IsGoodQuality() && theQnVectorUncorr[iDet]->GetN() != 0) {
        if(iDet==0) {
          qVecFullTPC[0] = theQnVectorUncorr[iDet]->Qx(2);
          qVecFullTPC[1] = theQnVectorUncorr[iDet]->Qy(2);
        }
        else if(iDet==1) {
          qVecPosTPC[0] = theQnVectorUncorr[iDet]->Qx(2);
          qVecPosTPC[1] = theQnVectorUncorr[iDet]->Qy(2);
        }
        else if(iDet==2) {
          qVecNegTPC[0] = theQnVectorUncorr[iDet]->Qx(2);
          qVecNegTPC[1] = theQnVectorUncorr[iDet]->Qy(2);
        }

        if(fUseQnFrameworkCorrq2) {
          theQnVectorCorr[iDet] = (AliQnCorrectionsQnVector*) pQvecList[iDet]->FindObject("rec"); //rec step for TPC
          if (theQnVectorCorr[iDet] && theQnVectorCorr[iDet]->IsGoodQuality() && theQnVectorCorr[iDet]->GetN() != 0) {
            if(iDet==0) {
              qRecFullTPC[0] = theQnVectorCorr[iDet]->Qx(2);
              qRecFullTPC[1] = theQnVectorCorr[iDet]->Qy(2);
              corrFullTPC[0] = qVecFullTPC[0] - qRecFullTPC[0];
              corrFullTPC[1] = qVecFullTPC[1] - qRecFullTPC[1];
            }
            else if(iDet==1) {
              qRecPosTPC[0] = theQnVectorCorr[iDet]->Qx(2);
              qRecPosTPC[1] = theQnVectorCorr[iDet]->Qy(2);
              corrPosTPC[0] = qVecPosTPC[0] - qRecPosTPC[0];
              corrPosTPC[1] = qVecPosTPC[1] - qRecPosTPC[1];
            }
            else if(iDet==2) {
              qRecNegTPC[0] = theQnVectorCorr[iDet]->Qx(2);
              qRecNegTPC[1] = theQnVectorCorr[iDet]->Qy(2);
              corrNegTPC[0] = qVecNegTPC[0] - qRecNegTPC[0];
              corrNegTPC[1] = qVecNegTPC[1] - qRecNegTPC[1];
            }
          }
          theQnVectorCorr[iDet]=0x0;
          theQnVectorCorr[iDet] = (AliQnCorrectionsQnVector*) pQvecList[iDet]->FindObject("twist"); //twist step for TPC
          if (theQnVectorCorr[iDet] && theQnVectorCorr[iDet]->IsGoodQuality() && theQnVectorUncorr[iDet]->GetN() != 0) {
            isTwistApplied[iDet]=kTRUE;
            if(iDet==0) {
              qTwistFullTPC[0] = theQnVectorCorr[iDet]->Qx(2);
              qTwistFullTPC[1] = theQnVectorCorr[iDet]->Qy(2);
              LbFullTPC[0] = (qRecFullTPC[0]-qTwistFullTPC[0])/qTwistFullTPC[1];
              LbFullTPC[1] = (qRecFullTPC[1]-qTwistFullTPC[1])/qTwistFullTPC[0];
            }
            else if(iDet==1) {
              qTwistPosTPC[0] = theQnVectorCorr[iDet]->Qx(2);
              qTwistPosTPC[1] = theQnVectorCorr[iDet]->Qy(2);
              LbPosTPC[0] = (qRecPosTPC[0]-qTwistPosTPC[0])/qTwistPosTPC[1];
              LbPosTPC[1] = (qRecPosTPC[1]-qTwistPosTPC[1])/qTwistPosTPC[0];
            }
            else if(iDet==2) {
              qTwistNegTPC[0] = theQnVectorCorr[iDet]->Qx(2);
              qTwistNegTPC[1] = theQnVectorCorr[iDet]->Qy(2);
              LbNegTPC[0] = (qRecNegTPC[0]-qTwistNegTPC[0])/qTwistNegTPC[1];
              LbNegTPC[1] = (qRecNegTPC[1]-qTwistNegTPC[1])/qTwistNegTPC[0];
            }
          }
        }
      }
    }
  }

  AliAnalysisVertexingHF *vHF=new AliAnalysisVertexingHF();
  std::vector<Double_t> invMassCand;
  std::vector<Double_t> invMassCand2;
  std::vector<Double_t> ptCand;
  std::vector<Double_t> deltaphiCand;
  std::vector<Double_t> q2Cand;
  std::vector<Double_t> cosnphiCand;
  std::vector<Double_t> sinnphiCand;
  std::vector<Double_t> phiCand;
  std::vector<Int_t> isSelectedCand;

  Int_t nSelCandInMassRange=0;
  std::vector<Int_t> dzerodaulab1; //vector of labels of D0 daugheter 1 necessary if kD0toKpiFromDstar is used with fOptD0FromDstar==0
  std::vector<Int_t> dzerodaulab2; //vector of labels of D0 daugheter 2 necessary if kD0toKpiFromDstar is used with fOptD0FromDstar==0
  Int_t nD0usedfordiffDstar=0;
  //Loop on D candidates
  for (Int_t iCand = 0; iCand < nCand; iCand++) {
    Bool_t isD0new = kTRUE;
    if(!(fDecChannel==4 && fOptD0FromDstar==0)) {
      d=(AliAODRecoDecayHF*)arrayProng->UncheckedAt(iCand);
      Bool_t isSelBit=kTRUE;
      if(fDecChannel==0) isSelBit=d->HasSelectionBit(AliRDHFCuts::kDplusCuts);
      if(fDecChannel==1 || fDecChannel==4) isSelBit=d->HasSelectionBit(AliRDHFCuts::kD0toKpiCuts);
      if(fDecChannel==2) isSelBit=kTRUE;
      if(fDecChannel==3) isSelBit=d->HasSelectionBit(AliRDHFCuts::kDsCuts);
      if(!isSelBit)continue;
      
      if(fDecChannel==0 || fDecChannel==3) {
        if(!vHF->FillRecoCand(aod,(AliAODRecoDecayHF3Prong*)d))continue;
      }
      else if(fDecChannel == 1 || fDecChannel==4) {
        if(!vHF->FillRecoCand(aod,(AliAODRecoDecayHF2Prong*)d))continue;
      }
      else if(fDecChannel == 2) {
        if(!vHF->FillRecoCasc(aod,((AliAODRecoCascadeHF*)d),kTRUE))continue;
      }
    }
    else if (fDecChannel == 4 && fOptD0FromDstar==0) {
      //Get the D* candidate, see if it is in the mass range and if it is the case assign the D0 daughter to d
      AliAODRecoDecayHF* dstar=(AliAODRecoDecayHF*)arrayProng->UncheckedAt(iCand);
      if(!vHF->FillRecoCasc(aod,((AliAODRecoCascadeHF*)dstar),kTRUE))continue;
      Double_t invMassDstar = ((AliAODRecoCascadeHF*)dstar)->DeltaInvMass();
      Double_t deltamassPDG=(TDatabasePDG::Instance()->GetParticle(413)->Mass())-(TDatabasePDG::Instance()->GetParticle(421)->Mass());
      Double_t sigma = 0.0008;
      if(invMassDstar<deltamassPDG-3*sigma || invMassDstar>deltamassPDG+3*sigma) continue;
      d = ((AliAODRecoCascadeHF*)dstar)->Get2Prong();
      AliAODTrack* dautrack1 = (AliAODTrack*)d->GetDaughter(0);
      AliAODTrack* dautrack2 = (AliAODTrack*)d->GetDaughter(1);
      Int_t lab1 = dautrack1->GetLabel();
      Int_t lab2 = dautrack2->GetLabel();
      for(UInt_t iDzero=0; iDzero<dzerodaulab1.size(); iDzero++) {
        if(lab1 == dzerodaulab1[iDzero] && lab2 == dzerodaulab2[iDzero]) {
          isD0new=kFALSE;
          break;
        }
      }
      if(isD0new) {
        dzerodaulab1.push_back(lab1);
        dzerodaulab2.push_back(lab2);
      }
    }
    
    Int_t ptbin=fRDCuts->PtBin(d->Pt());
    if(ptbin<0) {
      fhEventsInfo->Fill(15);
      continue;
    }
    Bool_t isFidAcc = fRDCuts->IsInFiducialAcceptance(d->Pt(),d->Y(absPdgMom));
    if(!isFidAcc)continue;
    Int_t isSelected = fRDCuts->IsSelected(d,AliRDHFCuts::kAll,aod);
    if(!isSelected)continue;
    if(fDecChannel==3) {
      Int_t isDsPhiKKpi = isSelected&4;
      Int_t isDsPhipiKK = isSelected&8;
      if(!isDsPhiKKpi & !isDsPhipiKK) continue;
    }
    if(fDecChannel==4 && fOptD0FromDstar==0 && !isD0new) {
      nD0usedfordiffDstar++;
      continue; //if D0 candidate already used for previous Dstar candidates, continue (the same D0 candidate can build many D* candidates)
    }
    if(fDecChannel==4 && fOptD0FromDstar==1) {
      
      Double_t deltamassPDG=(TDatabasePDG::Instance()->GetParticle(413)->Mass())-(TDatabasePDG::Instance()->GetParticle(421)->Mass());
      Int_t isD0fromDstar=kFALSE;
      Double_t invmassD0=-1;

      if(isSelected==1) {
        Int_t pdgdaughtersD0[2]={211,321};//pi+,K-
        invmassD0 = d->InvMass(2,(UInt_t*)pdgdaughtersD0); //D0
      }
      else if(isSelected==2) {
        Int_t pdgdaughtersD0bar[2]={321,211};//K+,pi-
        invmassD0 = d->InvMass(2,(UInt_t*)pdgdaughtersD0bar); //D0bar
      }
      
      AliAODTrack* dautrack1 = (AliAODTrack*)d->GetDaughter(0);
      AliAODTrack* dautrack2 = (AliAODTrack*)d->GetDaughter(1);
      Double_t mpion=TDatabasePDG::Instance()->GetParticle(211)->Mass();
      Double_t mkaon=TDatabasePDG::Instance()->GetParticle(321)->Mass();

      for(Int_t iTrk=0; iTrk<aod->GetNumberOfTracks(); iTrk++) {
        //check whether D0/D0bar comes from a Dstar, combining it with soft pions
        AliAODTrack* track=(AliAODTrack*)aod->GetTrack(iTrk);
        if(!track) continue;
        Short_t charge = track->Charge();
        //wrong charge sign --> continue
        if(isSelected==1 && charge<0) continue;
        if(isSelected==2 && charge>0) continue;
        
        //track does not pass soft pion selections --> continue
        if(!IsSoftPionSelected(track)) continue;
        
        if(isSelected==3 && charge>0) {
          Int_t pdgdaughtersD0[2]={211,321};//pi+,K-
          invmassD0 = d->InvMass(2,(UInt_t*)pdgdaughtersD0); //D0
        }
        else if(isSelected==3 && charge<0) {
          Int_t pdgdaughtersD0bar[2]={321,211};//K+,pi-
          invmassD0 = d->InvMass(2,(UInt_t*)pdgdaughtersD0bar); //D0bar
        }
        
        //compute invariant mass of D0 + soft pion
        Double_t energysum = 0;
        if(charge>0) {
          energysum = track->E(mpion)+dautrack1->E(mpion)+dautrack2->E(mkaon);
        }
        else if(charge<0) {
          energysum = track->E(mpion)+dautrack1->E(mkaon)+dautrack2->E(mpion);
        }
        Double_t psum2 = (track->Px()+dautrack1->Px()+dautrack2->Px())*(track->Px()+dautrack1->Px()+dautrack2->Px())+(track->Py()+dautrack1->Py()+dautrack2->Py())*(track->Py()+dautrack1->Py()+dautrack2->Py())+(track->Pz()+dautrack1->Pz()+dautrack2->Pz())*(track->Pz()+dautrack1->Pz()+dautrack2->Pz());
        Double_t invmassDstar = TMath::Sqrt(energysum*energysum-psum2);
        
        Double_t deltainvmassKpipi = invmassDstar-invmassD0;
        
        Double_t sigma = 0.0008;
        if(deltainvmassKpipi<deltamassPDG+3*sigma && deltainvmassKpipi>deltamassPDG-3*sigma) {
          isD0fromDstar=kTRUE;
          break;
        }
      }
      if(!isD0fromDstar) continue;
    }
    isSelectedCand.push_back(isSelected);
    
    fhEventsInfo->Fill(14); // candidate selected
    if(fDebug>3) printf("+++++++Is Selected\n");

    Float_t* invMass=0x0;
    Int_t nmasses;
    CalculateInvMasses(d,invMass,nmasses);

    invMassCand.push_back(invMass[0]);
    if(nmasses>1) {invMassCand2.push_back(invMass[1]);}
    ptCand.push_back(d->Pt());

    Bool_t ismassrange = kFALSE;
    if(isInMassRange(invMass[0],d->Pt()) || (nmasses>1 && isInMassRange(invMass[1],d->Pt()))) {
      ismassrange = kTRUE;
      nSelCandInMassRange++;
    }
    if(!fRequireMassForDauRemFromq2) {ismassrange=kTRUE;}

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
    if(fReadMC&&fUseAfterBurner)phi=fAfterBurner->GetNewAngle(d,arrayMC);
    Double_t deltaphi=GetPhiInRange(phi-eventplane);
    deltaphiCand.push_back(deltaphi);
    cosnphiCand.push_back(TMath::Cos(fHarmonic*phi));
    sinnphiCand.push_back(TMath::Sin(fHarmonic*phi));
    phiCand.push_back(phi);

    //fill the histograms with the appropriate method
    if(fFlowMethod!=kEvShape) {
      if(fDecChannel==0)FillDplus(d,arrayMC,ptbin,deltaphi,invMass,isSelected,icentr,phi,eta,Q1,Q2);
      else if(fDecChannel==1 || fDecChannel==4)FillD02p(d,arrayMC,ptbin,deltaphi,invMass,isSelected,icentr,phi,eta,Q1,Q2);
      else if(fDecChannel==2)FillDstar(d,arrayMC,ptbin,deltaphi,invMass,isSelected,icentr,phi,eta,Q1,Q2);
      else if(fDecChannel==3)FillDs(d,arrayMC,ptbin,deltaphi,invMass,isSelected,icentr,phi,eta,Q1,Q2);
    }
    else {
      std::vector<Int_t> daulab;
      if((fRemoveDauFromq2==0 || (fRemoveDauFromq2==1 && !ismassrange))) {
        if(fOnTheFlyTPCq2 && fDeltaEtaDmesonq2>0.) { //if deltaeta>0, remove tracks
          Double_t etaLims[2] = {fTPCEtaMin,fTPCEtaMax};
          if(fq2Meth==kq2TPC) {
            RemoveTracksInDeltaEtaFromOnTheFlyTPCq2(aod,eta,etaLims,qVecRemDauFullTPC,multQvecRemDauTPC[0],daulab);
            if(multQvecRemDauTPC[0]>0) {q2Cand.push_back(TMath::Sqrt(qVecRemDauFullTPC[0]*qVecRemDauFullTPC[0]+qVecRemDauFullTPC[1]*qVecRemDauFullTPC[1])/TMath::Sqrt(multQvecRemDauTPC[0]));}
            else {q2Cand.push_back(0);}
            //reset with all tracks for next candidate
            qVecRemDauFullTPC[0]=qVecFullTPC[0];
            qVecRemDauFullTPC[1]=qVecFullTPC[1];
            multQvecRemDauTPC[0]=multQvecTPC[0];
          }
          else if(fq2Meth==kq2PosTPC) {
            etaLims[0]=0.;
            RemoveTracksInDeltaEtaFromOnTheFlyTPCq2(aod,eta,etaLims,qVecRemDauPosTPC,multQvecRemDauTPC[1],daulab);
            if(multQvecRemDauTPC[1]>0) {q2Cand.push_back(TMath::Sqrt(qVecRemDauPosTPC[0]*qVecRemDauPosTPC[0]+qVecRemDauPosTPC[1]*qVecRemDauPosTPC[1])/TMath::Sqrt(multQvecRemDauTPC[1]));}
            else {q2Cand.push_back(0);}
            //reset with all tracks for next candidate
            qVecRemDauPosTPC[0]=qVecPosTPC[0];
            qVecRemDauPosTPC[1]=qVecPosTPC[1];
            multQvecRemDauTPC[1]=multQvecTPC[1];
          }
          else if(fq2Meth==kq2NegTPC) {
            etaLims[1]=0.;
            RemoveTracksInDeltaEtaFromOnTheFlyTPCq2(aod,eta,etaLims,qVecRemDauNegTPC,multQvecRemDauTPC[2],daulab);
            if(multQvecRemDauTPC[2]>0) {q2Cand.push_back(TMath::Sqrt(qVecRemDauNegTPC[0]*qVecRemDauNegTPC[0]+qVecRemDauNegTPC[1]*qVecRemDauNegTPC[1])/TMath::Sqrt(multQvecRemDauTPC[2]));}
            else {q2Cand.push_back(0);}
            //reset with all tracks for next candidate
            qVecRemDauNegTPC[0]=qVecNegTPC[0];
            qVecRemDauNegTPC[1]=qVecNegTPC[1];
            multQvecRemDauTPC[2]=multQvecTPC[2];
          }
        }
        else{ //if removal of daughter tracks disabled or enabled for single candidate and candidate not in mass range, assign global q2
          q2Cand.push_back(q2);
        }
      }
      else if((fq2Meth==kq2TPC || fq2Meth==kq2PosTPC || fq2Meth==kq2NegTPC) && fRemoveDauFromq2>0 && ismassrange) { //remove daughter tracks from q2
        Int_t nDau=3;
        if(fDecChannel==1 || fDecChannel==4 || (fDecChannel==2 && !fRemoverSoftPionFromq2)) {nDau=2;}
        AliAODTrack *dautrack[3] = {0x0,0x0,0x0};
        Double_t daueta = 0.;
        Double_t dauphi = 0.;
        Double_t daupt = 0.;
        Double_t dauqx = 0.;
        Double_t dauqy = 0.;

        if(fDecChannel!=2) { //D0, Dplus, Ds
          for(Int_t iDauTrk=0; iDauTrk<nDau; iDauTrk++) {
            dautrack[iDauTrk] = (AliAODTrack*)d->GetDaughter(iDauTrk);
          }
        }
        else {//D*
          AliAODRecoDecayHF2Prong* theD0particle = ((AliAODRecoCascadeHF*)d)->Get2Prong();
          dautrack[0] = (AliAODTrack*)theD0particle->GetDaughter(0);
          dautrack[1] = (AliAODTrack*)theD0particle->GetDaughter(1);
          if(fRemoverSoftPionFromq2) dautrack[2] = ((AliAODRecoCascadeHF*)d)->GetBachelor();
        }
        if(fOnTheFlyTPCq2) {
          for(Int_t iDauTrk=0; iDauTrk<nDau; iDauTrk++) { //remove daughter tracks from q2
            if(dautrack[iDauTrk]->TestFilterBit(BIT(8)) || dautrack[iDauTrk]->TestFilterBit(BIT(9))) { //if passes track cuts is used for q2 -> has to be removed
              Int_t lab = dautrack[iDauTrk]->GetLabel();
              daulab.push_back(lab);
              daueta = dautrack[iDauTrk]->Eta();
              dauphi = dautrack[iDauTrk]->Phi();
              daupt=dautrack[iDauTrk]->Pt();
              dauqx=TMath::Cos(2*dauphi);
              dauqy=TMath::Sin(2*dauphi);
              std::vector<Int_t>::iterator it = find(labrejtracks.begin(),labrejtracks.end(),lab);
              if((daupt>0.2 && daupt<5) && (daueta>fTPCEtaMin && daueta<fTPCEtaMax) && it==labrejtracks.end()) {//if is in right eta and pt region w.r.t. q2, and it was not removed by the random downsamping, remove
                qVecRemDauFullTPC[0] -= dauqx;
                qVecRemDauFullTPC[1] -= dauqy;
                multQvecRemDauTPC[0]--;
                if(daueta>0) {
                  qVecRemDauPosTPC[0] -= dauqx;
                  qVecRemDauPosTPC[1] -= dauqy;
                  multQvecRemDauTPC[1]--;
                }
                else {
                  qVecRemDauNegTPC[0] -= dauqx;
                  qVecRemDauNegTPC[1] -= dauqy;
                  multQvecRemDauTPC[2]--;
                }
              }
            }
          }
          if(fRemoveDauFromq2==1) { //remove only for the single candidate
            if(fq2Meth==kq2TPC) {
              if(fDeltaEtaDmesonq2>0) { //if deltaeta>0, remove also tracks in etagap
                Double_t etaLims[2] = {fTPCEtaMin,fTPCEtaMax};
                RemoveTracksInDeltaEtaFromOnTheFlyTPCq2(aod,eta,etaLims,qVecRemDauFullTPC,multQvecRemDauTPC[0],daulab);
              }
              qVecDefault[0]=qVecRemDauFullTPC[0];
              qVecDefault[1]=qVecRemDauFullTPC[1];
              multQvecDefault=multQvecRemDauTPC[0];
            }
            else if(fq2Meth==kq2PosTPC) {
              if(fDeltaEtaDmesonq2>0) { //if deltaeta>0, remove also tracks in etagap
                Double_t etaLims[2] = {0.,fTPCEtaMax};
                RemoveTracksInDeltaEtaFromOnTheFlyTPCq2(aod,eta,etaLims,qVecRemDauPosTPC,multQvecRemDauTPC[1],daulab);
              }
              qVecDefault[0]=qVecRemDauPosTPC[0];
              qVecDefault[1]=qVecRemDauPosTPC[1];
              multQvecDefault=multQvecRemDauTPC[1];
            }
            else if(fq2Meth==kq2NegTPC) {
              if(fDeltaEtaDmesonq2>0) { //if deltaeta>0, remove also tracks in etagap
                Double_t etaLims[2] = {fTPCEtaMin,0.};
                RemoveTracksInDeltaEtaFromOnTheFlyTPCq2(aod,eta,etaLims,qVecRemDauNegTPC,multQvecRemDauTPC[2],daulab);
              }
              qVecDefault[0]=qVecRemDauNegTPC[0];
              qVecDefault[1]=qVecRemDauNegTPC[1];
              multQvecDefault=multQvecRemDauTPC[2];
            }
            if(multQvecDefault>0) {q2Cand.push_back(TMath::Sqrt(qVecDefault[0]*qVecDefault[0]+qVecDefault[1]*qVecDefault[1])/TMath::Sqrt(multQvecDefault));}
            else {q2Cand.push_back(0);}

            //reset Q-vector with all tracks
            qVecRemDauFullTPC[0]=qVecFullTPC[0];
            qVecRemDauFullTPC[1]=qVecFullTPC[1];
            qVecRemDauPosTPC[0]=qVecPosTPC[0];
            qVecRemDauPosTPC[1]=qVecPosTPC[1];
            qVecRemDauNegTPC[0]=qVecNegTPC[0];
            qVecRemDauNegTPC[1]=qVecNegTPC[1];
            multQvecRemDauTPC[0]=multQvecTPC[0];
            multQvecRemDauTPC[1]=multQvecTPC[1];
            multQvecRemDauTPC[2]=multQvecTPC[2];
          }
        }
        else {
          for(Int_t iDauTrk = 0; iDauTrk < nDau; iDauTrk++) {
            if(dautrack[iDauTrk]->TestFilterBit(BIT(8))||dautrack[iDauTrk]->TestFilterBit(BIT(9))) {
              Int_t lab = dautrack[iDauTrk]->GetLabel();
              daueta = dautrack[iDauTrk]->Eta();
              dauphi = dautrack[iDauTrk]->Phi();
              daupt=dautrack[iDauTrk]->Pt();
              std::vector<Int_t>::iterator it = find(labrejtracks.begin(),labrejtracks.end(),lab);
              if((daupt>0.2 && daupt<5) && (daueta>fTPCEtaMin && daueta<fTPCEtaMax) && it==labrejtracks.end()) {
                dauqFullTPC[0] += TMath::Cos(2*dauphi);
                dauqFullTPC[1] += TMath::Sin(2*dauphi);
                nDauRemoved[0]++;
                if(daueta>0) {
                  dauqPosTPC[0] += TMath::Cos(2*dauphi);
                  dauqPosTPC[1] += TMath::Sin(2*dauphi);
                  nDauRemoved[1]++;
                }
                else {
                  dauqNegTPC[0] += TMath::Cos(2*dauphi);
                  dauqNegTPC[1] += TMath::Sin(2*dauphi);
                  nDauRemoved[2]++;
                }
              }
            }
          }
          if(fRemoveDauFromq2==1) { //remove only for the single candidate
            if(fq2Meth==kq2TPC) {
              q2Cand.push_back(GetTPCq2DauSubQnFramework(qVecFullTPC,multQvecTPC[0],nDauRemoved[0],dauqFullTPC,corrFullTPC,LbFullTPC,isTwistApplied[0]));
            }
            else if(fq2Meth==kq2PosTPC) {
              q2Cand.push_back(GetTPCq2DauSubQnFramework(qVecPosTPC,multQvecTPC[1],nDauRemoved[1],dauqPosTPC,corrPosTPC,LbPosTPC,isTwistApplied[1]));
            }
            else if(fq2Meth==kq2NegTPC) {
              q2Cand.push_back(GetTPCq2DauSubQnFramework(qVecNegTPC,multQvecTPC[2],nDauRemoved[2],dauqNegTPC,corrNegTPC,LbNegTPC,isTwistApplied[2]));
            }
            //reset Q-vectors of daughter tracks
            for(Int_t iDet=0; iDet<3; iDet++) {nDauRemoved[iDet]=0;}
            dauqFullTPC[0]=0.;
            dauqFullTPC[1]=0.;
            dauqPosTPC[0]=0.;
            dauqPosTPC[1]=0.;
            dauqNegTPC[0]=0.;
            dauqNegTPC[1]=0.;
          }
        }
      }
      daulab.clear(); //reset vector of daughter labels
    }
    delete [] invMass;
  }

  UInt_t nSelCand = invMassCand.size();
  std::vector<Double_t> q2CandFill;
  
  if(fFlowMethod==kEvShape) {
    //if removed the daughter tracks for all the candidates, recompute q2 from q-vector w/o daughter tracks
    if((fq2Meth==kq2TPC || fq2Meth==kq2PosTPC || fq2Meth==kq2NegTPC) && fRemoveDauFromq2==2) {
      if(fOnTheFlyTPCq2) {
        if(multQvecRemDauTPC[0]>0) {q2FullTPC = TMath::Sqrt(qVecRemDauFullTPC[0]*qVecRemDauFullTPC[0]+qVecRemDauFullTPC[1]*qVecRemDauFullTPC[1])/TMath::Sqrt(multQvecRemDauTPC[0]);}
        else {q2FullTPC=0.;}
        if(multQvecRemDauTPC[1]>0) {q2PosTPC = TMath::Sqrt(qVecRemDauPosTPC[0]*qVecRemDauPosTPC[0]+qVecRemDauPosTPC[1]*qVecRemDauPosTPC[1])/TMath::Sqrt(multQvecRemDauTPC[1]);}
        else {q2PosTPC=0.;}
        if(multQvecRemDauTPC[2]>0) {q2NegTPC = TMath::Sqrt(qVecRemDauNegTPC[0]*qVecRemDauNegTPC[0]+qVecRemDauNegTPC[1]*qVecRemDauNegTPC[1])/TMath::Sqrt(multQvecRemDauTPC[2]);}
        else {q2NegTPC=0.;}
      }
      else {
        q2FullTPC=GetTPCq2DauSubQnFramework(qVecFullTPC,multQvecTPC[0],nDauRemoved[0],dauqFullTPC,corrFullTPC,LbFullTPC,isTwistApplied[0]);
        q2PosTPC=GetTPCq2DauSubQnFramework(qVecPosTPC,multQvecTPC[1],nDauRemoved[1],dauqPosTPC,corrPosTPC,LbPosTPC,isTwistApplied[1]);
        q2NegTPC=GetTPCq2DauSubQnFramework(qVecNegTPC,multQvecTPC[2],nDauRemoved[2],dauqNegTPC,corrNegTPC,LbNegTPC,isTwistApplied[2]);

        multQvecRemDauTPC[0]=multQvecTPC[0]-nDauRemoved[0];
        multQvecRemDauTPC[1]=multQvecTPC[1]-nDauRemoved[1];
        multQvecRemDauTPC[2]=multQvecTPC[2]-nDauRemoved[2];
      }
      if(fq2Meth==kq2TPC) {q2=q2FullTPC;}
      else if(fq2Meth==kq2PosTPC) {q2=q2PosTPC;}
      else if(fq2Meth==kq2NegTPC) {q2=q2NegTPC;}
    }
    else if((fq2Meth==kq2TPC || fq2Meth==kq2PosTPC || fq2Meth==kq2NegTPC) && fOnTheFlyTPCq2 && fDeltaEtaDmesonq2>0) { //remove a random delta eta from global q2 to have consistent multiplicity
      std::vector<Int_t> daulabdummy; //dummy -> not used
      Double_t etaLims[2] = {fTPCEtaMin,fTPCEtaMax};
      if(fq2Meth==kq2TPC) {
        Double_t etarndm = fTPCEtaMin+gRandom->Rndm()*(fTPCEtaMax-fTPCEtaMin);
        RemoveTracksInDeltaEtaFromOnTheFlyTPCq2(aod,etarndm,etaLims,qVecFullTPC,multQvecTPC[0],daulabdummy);
        if(multQvecTPC[0]>0) {q2=TMath::Sqrt(qVecFullTPC[0]*qVecFullTPC[0]+qVecFullTPC[1]*qVecFullTPC[1])/TMath::Sqrt(multQvecTPC[0]);}
        else {q2=0;}
      }
      else if(fq2Meth==kq2PosTPC) {
        Double_t etarndm = gRandom->Rndm()*fTPCEtaMax;
        etaLims[0]=0.;
        RemoveTracksInDeltaEtaFromOnTheFlyTPCq2(aod,etarndm,etaLims,qVecPosTPC,multQvecTPC[1],daulabdummy);
        if(multQvecTPC[1]>0) {q2=TMath::Sqrt(qVecPosTPC[0]*qVecPosTPC[0]+qVecPosTPC[1]*qVecPosTPC[1])/TMath::Sqrt(multQvecTPC[1]);}
        else {q2=0;}
      }
      else if(fq2Meth==kq2NegTPC) {
        Double_t etarndm = fTPCEtaMin+gRandom->Rndm()*(0-fTPCEtaMin);
        etaLims[1]=0.;
        RemoveTracksInDeltaEtaFromOnTheFlyTPCq2(aod,etarndm,etaLims,qVecNegTPC,multQvecTPC[2],daulabdummy);
        if(multQvecTPC[2]>0) {q2=TMath::Sqrt(qVecNegTPC[0]*qVecNegTPC[0]+qVecNegTPC[1]*qVecNegTPC[1])/TMath::Sqrt(multQvecTPC[2]);}
        else {q2=0;}
      }
    }

    Double_t q2fill = q2;
    Double_t q2percentile = -1.;
    Double_t q2percentileFullTPC = -1.;
    Double_t q2percentileVZEROC = -1.;
    Double_t q2percentileVZEROA = -1.;
    const Int_t ncentbins = (fMaxCentr-fMinCentr)*10/fCentBinSizePerMil;
    TSpline3* q2spline=0x0;
    TSpline3* q2splineFullTPC=0x0;
    TSpline3* q2splineVZEROC=0x0;
    TSpline3* q2splineVZEROA=0x0;
    Int_t splineindex = 0;
    if(fq2Meth==kq2NegTPC) {splineindex=1;}
    else if(fq2Meth==kq2PosTPC) {splineindex=2;}
    else if(fq2Meth==kq2VZERO) {splineindex=3;}
    else if(fq2Meth==kq2VZEROA) {splineindex=4;}
    else if(fq2Meth==kq2VZEROC) {splineindex=5;}

    if(fPercentileq2) {
      for(Int_t iCentr=0; iCentr<ncentbins; iCentr++) {
        if(centr>fMinCentr+iCentr*fCentBinSizePerMil/10 && centr<fMinCentr+(iCentr+1)*fCentBinSizePerMil/10) {
          q2spline = (TSpline3*)fq2SplinesList[splineindex]->FindObject(Form("sq2Int_centr_%d_%d",fMinCentr+iCentr*fCentBinSizePerMil/10,fMinCentr+(iCentr+1)*fCentBinSizePerMil/10));
          q2splineFullTPC = (TSpline3*)fq2SplinesList[0]->FindObject(Form("sq2Int_centr_%d_%d",fMinCentr+iCentr*fCentBinSizePerMil/10,fMinCentr+(iCentr+1)*fCentBinSizePerMil/10));
          q2splineVZEROC = (TSpline3*)fq2SplinesList[5]->FindObject(Form("sq2Int_centr_%d_%d",fMinCentr+iCentr*fCentBinSizePerMil/10,fMinCentr+(iCentr+1)*fCentBinSizePerMil/10));
          q2splineVZEROA = (TSpline3*)fq2SplinesList[4]->FindObject(Form("sq2Int_centr_%d_%d",fMinCentr+iCentr*fCentBinSizePerMil/10,fMinCentr+(iCentr+1)*fCentBinSizePerMil/10));
          break;
        }
      }
      if(!q2spline || !q2splineFullTPC || !q2splineVZEROC) {AliFatal("Centrality binning and centrality intervals of q2 splines do not match!");}
      q2percentile=q2spline->Eval(q2);
      q2percentileFullTPC=q2splineFullTPC->Eval(q2FullTPC);
      q2percentileVZEROC=q2splineVZEROC->Eval(q2VZEROC);
      q2percentileVZEROA=q2splineVZEROA->Eval(q2VZEROA);
      q2fill=q2percentile;
      if(fRemoveDauFromq2==1 || fDeltaEtaDmesonq2>0) {
        for(UInt_t iCand=0; iCand<q2Cand.size(); iCand++) {
          q2CandFill.push_back(q2spline->Eval(q2Cand[iCand]));
        }
      }
      else {
        for(UInt_t iCand=0; iCand<nSelCand; iCand++) {
          q2CandFill.push_back(q2fill);
        }
      }
    }
    else {
      if(fRemoveDauFromq2==1 || fDeltaEtaDmesonq2>0) {
        for(UInt_t iCand=0; iCand<q2Cand.size(); iCand++) {
          q2CandFill.push_back(q2Cand[iCand]);
        }
      }
      else {
        for(UInt_t iCand=0; iCand<nSelCand; iCand++) {
          q2CandFill.push_back(q2fill);
        }
      }
    }
    
    //fill percentile q2 vs. q2 vs. centrality histogram
    ((TH3F*)fOutput->FindObject("hPercq2vsq2vsCentr"))->Fill(centr,q2,q2percentile);

    //fill q2 vs. centrality histograms (q2 percentile calibration)
    ((TH2F*)fOutput->FindObject(Form("hq2vsCentr%s",fDetTPCConfName[0].Data())))->Fill(centr,q2FullTPC);
    ((TH2F*)fOutput->FindObject(Form("hq2vsCentr%s",fDetTPCConfName[1].Data())))->Fill(centr,q2NegTPC);
    ((TH2F*)fOutput->FindObject(Form("hq2vsCentr%s",fDetTPCConfName[2].Data())))->Fill(centr,q2PosTPC);
    ((TH2F*)fOutput->FindObject(Form("hq2vsCentr%s",fDetV0ConfName[0].Data())))->Fill(centr,q2VZERO);
    ((TH2F*)fOutput->FindObject(Form("hq2vsCentr%s",fDetV0ConfName[1].Data())))->Fill(centr,q2VZEROA);
    ((TH2F*)fOutput->FindObject(Form("hq2vsCentr%s",fDetV0ConfName[2].Data())))->Fill(centr,q2VZEROC);

    //If enabled, fill EP angle vs. q2 vs. centrality histograms
    if(fEPVsq2VsCent) {
      for(Int_t iDet=0; iDet<3; iDet++) {
        ((TH3F*)fOutput->FindObject(Form("hEvPlaneQncorr%s%sVsq2VsCent",fDetTPCConfName[iDet].Data(),fNormMethod.Data())))->Fill(centr,q2fill,eventplaneqncorrTPC[iDet]);
        ((TH3F*)fOutput->FindObject(Form("hEvPlaneQncorr%s%sVsq2VsCent",fDetV0ConfName[iDet].Data(),fNormMethod.Data())))->Fill(centr,q2fill,eventplaneqncorrVZERO[iDet]);
      }
    }

    //fill mult vs. centrality histo (EvShape)
    ((TH2F*)fOutput->FindObject("hMultVsCentFullTPC"))->Fill(centr,multQvecRemDauTPC[0]);
    ((TH2F*)fOutput->FindObject("hMultVsCentPosTPC"))->Fill(centr,multQvecRemDauTPC[1]);
    ((TH2F*)fOutput->FindObject("hMultVsCentNegTPC"))->Fill(centr,multQvecRemDauTPC[2]);
    
    //fill q2 correlation histograms
    ((TH2F*)fOutput->FindObject("hq2TPCPosEtaVsNegEta"))->Fill(q2NegTPC,q2PosTPC);
    ((TH2F*)fOutput->FindObject("hq2TPCFullEtaVsNegEta"))->Fill(q2NegTPC,q2FullTPC);
    ((TH2F*)fOutput->FindObject("hq2TPCFullEtaVsPosEta"))->Fill(q2PosTPC,q2FullTPC);

    if(!fPercentileq2) {
      ((TH2F*)fOutput->FindObject("hq2TPCVsq2VZEROC"))->Fill(q2FullTPC,q2VZEROC);
      ((TH2F*)fOutput->FindObject("hq2TPCVsq2VZEROA"))->Fill(q2FullTPC,q2VZEROA);
      if(nSelCand>0) {
        ((TH2F*)fOutput->FindObject("hq2TPCVsq2VZEROCCand"))->Fill(q2FullTPC,q2VZEROC);
        ((TH2F*)fOutput->FindObject("hq2TPCVsq2VZEROACand"))->Fill(q2FullTPC,q2VZEROA);
      }
      if(nSelCandInMassRange>0) {
        ((TH2F*)fOutput->FindObject("hq2TPCVsq2VZEROCCandInMass"))->Fill(q2FullTPC,q2VZEROC);
        ((TH2F*)fOutput->FindObject("hq2TPCVsq2VZEROACandInMass"))->Fill(q2FullTPC,q2VZEROA);
      }
    }
    else {
      ((TH2F*)fOutput->FindObject("hq2TPCVsq2VZEROC"))->Fill(q2percentileFullTPC,q2percentileVZEROC);
      ((TH2F*)fOutput->FindObject("hq2TPCVsq2VZEROA"))->Fill(q2percentileFullTPC,q2percentileVZEROA);
      if(nSelCand>0) {
        ((TH2F*)fOutput->FindObject("hq2TPCVsq2VZEROCCand"))->Fill(q2percentileFullTPC,q2percentileVZEROC);
        ((TH2F*)fOutput->FindObject("hq2TPCVsq2VZEROACand"))->Fill(q2percentileFullTPC,q2percentileVZEROA);
      }
      if(nSelCandInMassRange>0) {
        ((TH2F*)fOutput->FindObject("hq2TPCVsq2VZEROCCandInMass"))->Fill(q2percentileFullTPC,q2percentileVZEROC);
        ((TH2F*)fOutput->FindObject("hq2TPCVsq2VZEROACandInMass"))->Fill(q2percentileFullTPC,q2percentileVZEROA);
      }
    }
    
    Int_t tracklets=AliVertexingHFUtils::GetNumberOfTrackletsInEtaRange(aod,-1.,1.);
    if(fEnableNtrklHistos) {
      ((TH3F*)fOutput->FindObject("hNtrklVsq2VsCent"))->Fill(centr,q2fill,(Double_t)tracklets);
      if(nSelCand>0) {
        ((TH3F*)fOutput->FindObject("hNtrklVsq2VsCentCand"))->Fill(centr,q2fill,(Double_t)tracklets);
        if(nSelCandInMassRange>0) {
          ((TH3F*)fOutput->FindObject("hNtrklVsq2VsCentCandInMass"))->Fill(centr,q2fill,(Double_t)tracklets);
        }
      }
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

    //fill resolution histograms
    ((TH2F*)fOutput->FindObject(Form("hEvPlaneReso1Vsq2%s",fCentrBinName.Data())))->Fill(planereso,q2fill); //RP resolution vs q2
    if(nSubEvents==3){
      ((TH2F*)fOutput->FindObject(Form("hEvPlaneReso2Vsq2%s",fCentrBinName.Data())))->Fill(TMath::Cos(fHarmonic*deltaSubAC),q2fill);
      ((TH2F*)fOutput->FindObject(Form("hEvPlaneReso3Vsq2%s",fCentrBinName.Data())))->Fill(TMath::Cos(fHarmonic*deltaSubBC),q2fill);
    }

    //fill histo for number of Dstar built using the same Dzero (only in case of Dzero from Dstar analysis)
    if(fDecChannel==kD0toKpiFromDstar && fOptD0FromDstar==0) {
      ((TH1F*)fOutput->FindObject("hNumDstarCandFromSameDzeroCand"))->Fill(nD0usedfordiffDstar);
    }
    
    //fill q2Cand vs. q2Event histo and THnSparseF for event-shape engineering
    for(UInt_t iSelCand=0; iSelCand<nSelCand; iSelCand++) {
      ((TH2F*)fOutput->FindObject("hq2CandVsq2Event"))->Fill(q2fill,q2CandFill[iSelCand]);

      if((fDecChannel==0 || fDecChannel==2) && isSelectedCand[iSelCand]) {
        Double_t sparsearray[9] = {invMassCand[iSelCand],ptCand[iSelCand],deltaphiCand[iSelCand],q2CandFill[iSelCand],centr,cosnphiCand[iSelCand],sinnphiCand[iSelCand],phiCand[iSelCand],(Double_t)tracklets};
        fHistMassPtPhiq2Centr->Fill(sparsearray);
      }
      else if(fDecChannel==1 || fDecChannel==4) {
        if(fSeparateD0D0bar) {
          if(isSelectedCand[iSelCand]==1 || isSelectedCand[iSelCand]==3) {
            Double_t sparsearray1[10] = {invMassCand[iSelCand],ptCand[iSelCand],deltaphiCand[iSelCand],q2CandFill[iSelCand],centr,cosnphiCand[iSelCand],sinnphiCand[iSelCand],phiCand[iSelCand],(Double_t)tracklets,(Double_t)isSelectedCand[iSelCand]};
            fHistMassPtPhiq2Centr->Fill(sparsearray1);
          }
          if(isSelectedCand[iSelCand]>=2) {
            Double_t sparsearray2[10] = {invMassCand2[iSelCand],ptCand[iSelCand],deltaphiCand[iSelCand],q2CandFill[iSelCand],centr,cosnphiCand[iSelCand],sinnphiCand[iSelCand],phiCand[iSelCand],(Double_t)tracklets,(Double_t)isSelectedCand[iSelCand]};
            fHistMassPtPhiq2Centr->Fill(sparsearray2);
          }
        }
        else {
          if(isSelectedCand[iSelCand]==1 || isSelectedCand[iSelCand]==3) {
            Double_t sparsearray1[9] = {invMassCand[iSelCand],ptCand[iSelCand],deltaphiCand[iSelCand],q2CandFill[iSelCand],centr,cosnphiCand[iSelCand],sinnphiCand[iSelCand],phiCand[iSelCand],(Double_t)tracklets};
            fHistMassPtPhiq2Centr->Fill(sparsearray1);
          }
          if(isSelectedCand[iSelCand]>=2) {
            Double_t sparsearray2[9] = {invMassCand2[iSelCand],ptCand[iSelCand],deltaphiCand[iSelCand],q2CandFill[iSelCand],centr,cosnphiCand[iSelCand],sinnphiCand[iSelCand],phiCand[iSelCand],(Double_t)tracklets};
            fHistMassPtPhiq2Centr->Fill(sparsearray2);
          }
        }
      }
      else if(fDecChannel==3) {
        if(isSelectedCand[iSelCand]==1 || isSelectedCand[iSelCand]==3) {
          Double_t sparsearray1[9] = {invMassCand[iSelCand],ptCand[iSelCand],deltaphiCand[iSelCand],q2CandFill[iSelCand],centr,cosnphiCand[iSelCand],sinnphiCand[iSelCand],phiCand[iSelCand],(Double_t)tracklets};
          fHistMassPtPhiq2Centr->Fill(sparsearray1);
        }
        if(isSelectedCand[iSelCand]>=2) {
          Double_t sparsearray2[9] = {invMassCand2[iSelCand],ptCand[iSelCand],deltaphiCand[iSelCand],q2CandFill[iSelCand],centr,cosnphiCand[iSelCand],sinnphiCand[iSelCand],phiCand[iSelCand],(Double_t)tracklets};
          fHistMassPtPhiq2Centr->Fill(sparsearray2);
        }
      }
    }
  }

  fHistCandVsCent->Fill(evCentr,nSelCand);
  fHistCandMassRangeVsCent->Fill(evCentr,nSelCandInMassRange);

  invMassCand.clear();
  invMassCand2.clear();
  ptCand.clear();
  deltaphiCand.clear();
  q2Cand.clear();
  q2CandFill.clear();
  cosnphiCand.clear();
  sinnphiCand.clear();
  isSelectedCand.clear();
  phiCand.clear();
  labrejtracks.clear();
  dzerodaulab1.clear();
  dzerodaulab2.clear();

  delete vHF;
  PostData(1,fhEventsInfo);
  PostData(2,fOutput);

  return;
}

void AliAnalysisTaskSEHFvn::CreateSparseForEvShapeAnalysis() {
  /// Sparse for v2 analysis as a function of q2

  Int_t nptbins=200;
  Double_t ptmin=0.;
  Double_t ptmax=50.;

  Int_t ndeltaphibins=96;
  Double_t mindeltaphi=0.;
  Double_t maxdeltaphi = TMath::Pi();

  Int_t nq2bins=500;
  Double_t minq2=0.;
  Double_t maxq2=10.;
  if(fPercentileq2) {
    nq2bins=100;
    minq2=0.;
    maxq2=100.;
  }
  
  Int_t ncentbins=(fMaxCentr-fMinCentr)*10/fCentBinSizePerMil;

  Int_t nphibins=18;
  Double_t phimin=0.;
  Double_t phimax=2*TMath::Pi();

  Int_t nNtrkBins=100;
  Double_t Ntrkmin=0.;
  Double_t Ntrkmax=5000.;

  TString massaxisname;
  if(fDecChannel==0) massaxisname = "M_{K#pi#pi} (GeV/c^{2})";
  else if(fDecChannel==1 || fDecChannel==4) massaxisname = "M_{K#pi} (GeV/c^{2})";
  else if(fDecChannel==2) massaxisname = "M_{K#pi#pi}-M_{K#pi} (GeV/c^{2})";
  else if(fDecChannel==3) massaxisname = "M_{KK#pi} (GeV/c^{2})";

  Int_t naxes=9;

  if(fSeparateD0D0bar && (fDecChannel==1 || fDecChannel==4)) {
    Int_t npartantipartbins=3;
    Double_t minpartantipart=0.5;
    Double_t maxpartantipart=3.5;

    Int_t nbins[10]={fNMassBins,nptbins,ndeltaphibins,nq2bins,ncentbins,100,100,nphibins,nNtrkBins,npartantipartbins};
    Double_t xmin[10]={fLowmasslimit,ptmin,mindeltaphi,minq2,(Double_t)fMinCentr,-1.,-1.,phimin,Ntrkmin,(Double_t)minpartantipart};
    Double_t xmax[10]={fUpmasslimit,ptmax,maxdeltaphi,maxq2,(Double_t)fMaxCentr,1.,1.,phimax,Ntrkmax,(Double_t)maxpartantipart};

    naxes=10;

    fHistMassPtPhiq2Centr=new THnSparseF("hMassPtPhiq2Centr","Mass vs. pt vs. #Delta#phi vs. q_{2} vs. centr vs. D0-D0bar",naxes,nbins,xmin,xmax);
  }
  else {
    Int_t nbins[9]={fNMassBins,nptbins,ndeltaphibins,nq2bins,ncentbins,100,100,nphibins,nNtrkBins};
    Double_t xmin[9]={fLowmasslimit,ptmin,mindeltaphi,minq2,(Double_t)fMinCentr,-1.,-1.,phimin,Ntrkmin};
    Double_t xmax[9]={fUpmasslimit,ptmax,maxdeltaphi,maxq2,(Double_t)fMaxCentr,1.,1.,phimax,Ntrkmax};

    fHistMassPtPhiq2Centr=new THnSparseF("hMassPtPhiq2Centr","Mass vs. pt vs. #Delta#phi vs. q_{2} vs. centr",naxes,nbins,xmin,xmax);
  }

  TString q2axisname="q_{2}^{TPC}";
  if(fq2Meth==kq2PosTPC) {q2axisname="q_{2}^{pos TPC}";}
  else if(fq2Meth==kq2NegTPC) {q2axisname="q_{2}^{neg TPC}";}
  else if(fq2Meth==kq2VZERO) {q2axisname="q_{2}^{V0}";}
  else if(fq2Meth==kq2VZEROA) {q2axisname="q_{2}^{V0A}";}
  else if(fq2Meth==kq2VZEROC) {q2axisname="q_{2}^{V0C}";}

  if(fPercentileq2) {q2axisname += " (%)";}
  
  TString axTit[10]={massaxisname,"p_{T} (GeV/c)","#Delta#varphi",q2axisname,"Centrality (%)",Form("Cos(%d#varphi_{D})",fHarmonic),Form("Sin(%d#varphi_{D})",fHarmonic),"#varphi_{D}","N_{tracklets}","part-antipart"};

  for(Int_t iax=0; iax<naxes; iax++)
    fHistMassPtPhiq2Centr->GetAxis(iax)->SetTitle(axTit[iax].Data());

  fOutput->Add(fHistMassPtPhiq2Centr);
}

//***************************************************************************

// Methods used in the UserExec

void AliAnalysisTaskSEHFvn::CalculateInvMasses(AliAODRecoDecayHF* d,Float_t*& masses,Int_t& nmasses){
  //Calculates all the possible invariant masses for each candidate
  //NB: the implementation for each candidate is responsibility of the corresponding developer

  if(fDecChannel==0){
    //D+ -- Giacomo
    nmasses=1;
    masses=new Float_t[nmasses];
    Int_t pdgdaughters[3] = {211,321,211};
    masses[0]=d->InvMass(3,(UInt_t*)pdgdaughters);
  }
  if(fDecChannel==1 || fDecChannel==4){
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
void AliAnalysisTaskSEHFvn::FillDplus(AliAODRecoDecayHF* d,TClonesArray *arrayMC,Int_t ptbin,Float_t deltaphi, const Float_t* masses,Int_t isSel,Int_t icentr, Double_t phiD, Double_t etaD, Double_t Q1[2], Double_t Q2[2]){
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
    if(fEvPlaneDet==kFullTPC || fEvPlaneDet==kPosTPC || fEvPlaneDet==kNegTPC) {
      if(etaD>0) {
        scalprod[0] = TMath::Cos(fHarmonic*phiD)*Q2[0]+TMath::Sin(fHarmonic*phiD)*Q2[1];
        ((TH2F*)fOutput->FindObject(Form("hMassScalProdu1Q2_pt%dcentr%d_%d",ptbin,icentrmin,icentr)))->Fill(scalprod[0],masses[0]);
      }
      else {
        scalprod[1] = TMath::Cos(fHarmonic*phiD)*Q1[0]+TMath::Sin(fHarmonic*phiD)*Q1[1];
        ((TH2F*)fOutput->FindObject(Form("hMassScalProdu2Q1_pt%dcentr%d_%d",ptbin,icentrmin,icentr)))->Fill(scalprod[1],masses[0]);
      }
    }
    else if(fEvPlaneDet==kFullV0){
      scalprod[0] = TMath::Cos(fHarmonic*phiD)*Q2[0]+TMath::Sin(fHarmonic*phiD)*Q2[1];
      scalprod[1] = TMath::Cos(fHarmonic*phiD)*Q1[0]+TMath::Sin(fHarmonic*phiD)*Q1[1];
      ((TH2F*)fOutput->FindObject(Form("hMassScalProdu3Q1_pt%dcentr%d_%d",ptbin,icentrmin,icentr)))->Fill(scalprod[0],masses[0]);
      ((TH2F*)fOutput->FindObject(Form("hMassScalProdu3Q2_pt%dcentr%d_%d",ptbin,icentrmin,icentr)))->Fill(scalprod[1],masses[0]);
    }
    else if(fEvPlaneDet==kV0A || fEvPlaneDet==kV0C){
      scalprod[0] = TMath::Cos(fHarmonic*phiD)*Q1[0]+TMath::Sin(fHarmonic*phiD)*Q1[1];
      ((TH2F*)fOutput->FindObject(Form("hMassScalProdu3Q1_pt%dcentr%d_%d",ptbin,icentrmin,icentr)))->Fill(scalprod[0],masses[0]);
    }
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
void AliAnalysisTaskSEHFvn::FillD02p(AliAODRecoDecayHF* d,TClonesArray *arrayMC,Int_t ptbin,Float_t deltaphi, const Float_t* masses,Int_t isSel,Int_t icentr,Double_t phiD, Double_t etaD, Double_t Q1[2], Double_t Q2[2]){

  //D0->Kpi channel

  //mass histograms
  if(!masses){
    if(fDebug>3)AliWarning("Masses not calculated\n");
    return;
  }
  Int_t icentrmin=icentr-fCentBinSizePerMil;
  Double_t scalprod[2]={-2.,-2.};
  if(isSel==1 || isSel==3) {
    if(fFlowMethod==kSP) {
      if(fEvPlaneDet==kFullTPC || fEvPlaneDet==kPosTPC || fEvPlaneDet==kNegTPC) {
        if(etaD>0) {
          scalprod[0] = TMath::Cos(fHarmonic*phiD)*Q2[0]+TMath::Sin(fHarmonic*phiD)*Q2[1];
          ((TH2F*)fOutput->FindObject(Form("hMassScalProdu1Q2_pt%dcentr%d_%d",ptbin,icentrmin,icentr)))->Fill(scalprod[0],masses[0]);
        }
        else {
          scalprod[1] = TMath::Cos(fHarmonic*phiD)*Q1[0]+TMath::Sin(fHarmonic*phiD)*Q1[1];
          ((TH2F*)fOutput->FindObject(Form("hMassScalProdu2Q1_pt%dcentr%d_%d",ptbin,icentrmin,icentr)))->Fill(scalprod[1],masses[0]);
        }
      }
      else if(fEvPlaneDet==kFullV0){
        scalprod[0] = TMath::Cos(fHarmonic*phiD)*Q2[0]+TMath::Sin(fHarmonic*phiD)*Q2[1];
        scalprod[1] = TMath::Cos(fHarmonic*phiD)*Q1[0]+TMath::Sin(fHarmonic*phiD)*Q1[1];
        ((TH2F*)fOutput->FindObject(Form("hMassScalProdu3Q1_pt%dcentr%d_%d",ptbin,icentrmin,icentr)))->Fill(scalprod[0],masses[0]);
        ((TH2F*)fOutput->FindObject(Form("hMassScalProdu3Q2_pt%dcentr%d_%d",ptbin,icentrmin,icentr)))->Fill(scalprod[1],masses[0]);
      }
      else if(fEvPlaneDet==kV0A || fEvPlaneDet==kV0C){
        scalprod[0] = TMath::Cos(fHarmonic*phiD)*Q1[0]+TMath::Sin(fHarmonic*phiD)*Q1[1];
        ((TH2F*)fOutput->FindObject(Form("hMassScalProdu3Q1_pt%dcentr%d_%d",ptbin,icentrmin,icentr)))->Fill(scalprod[0],masses[0]);
      }
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
      if(fEvPlaneDet==kFullTPC || fEvPlaneDet==kPosTPC || fEvPlaneDet==kNegTPC) {
        if(etaD>0) {
          scalprod[0] = TMath::Cos(fHarmonic*phiD)*Q2[0]+TMath::Sin(fHarmonic*phiD)*Q2[1];
          ((TH2F*)fOutput->FindObject(Form("hMassScalProdu1Q2_pt%dcentr%d_%d",ptbin,icentrmin,icentr)))->Fill(scalprod[0],masses[1]);
        }
        else {
          scalprod[1] = TMath::Cos(fHarmonic*phiD)*Q1[0]+TMath::Sin(fHarmonic*phiD)*Q1[1];
          ((TH2F*)fOutput->FindObject(Form("hMassScalProdu2Q1_pt%dcentr%d_%d",ptbin,icentrmin,icentr)))->Fill(scalprod[1],masses[1]);
        }
      }
      else if(fEvPlaneDet==kFullV0){
        scalprod[0] = TMath::Cos(fHarmonic*phiD)*Q2[0]+TMath::Sin(fHarmonic*phiD)*Q2[1];
        scalprod[1] = TMath::Cos(fHarmonic*phiD)*Q1[0]+TMath::Sin(fHarmonic*phiD)*Q1[1];
        ((TH2F*)fOutput->FindObject(Form("hMassScalProdu3Q1_pt%dcentr%d_%d",ptbin,icentrmin,icentr)))->Fill(scalprod[0],masses[1]);
        ((TH2F*)fOutput->FindObject(Form("hMassScalProdu3Q2_pt%dcentr%d_%d",ptbin,icentrmin,icentr)))->Fill(scalprod[1],masses[1]);
      }
      else if(fEvPlaneDet==kV0A || fEvPlaneDet==kV0C){
        scalprod[0] = TMath::Cos(fHarmonic*phiD)*Q1[0]+TMath::Sin(fHarmonic*phiD)*Q1[1];
        ((TH2F*)fOutput->FindObject(Form("hMassScalProdu3Q1_pt%dcentr%d_%d",ptbin,icentrmin,icentr)))->Fill(scalprod[0],masses[1]);
      }
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
void AliAnalysisTaskSEHFvn::FillDstar(AliAODRecoDecayHF* d,TClonesArray *arrayMC,Int_t ptbin,Float_t deltaphi, const Float_t* masses,Int_t isSel,Int_t icentr, Double_t phiD, Double_t etaD, Double_t Q1[2], Double_t Q2[2]){
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
    if(fEvPlaneDet==kFullTPC || fEvPlaneDet==kPosTPC || fEvPlaneDet==kNegTPC) {
      if(etaD>0) {
        scalprod[0] = TMath::Cos(fHarmonic*phiD)*Q2[0]+TMath::Sin(fHarmonic*phiD)*Q2[1];
        ((TH2F*)fOutput->FindObject(Form("hMassScalProdu1Q2_pt%dcentr%d_%d",ptbin,icentrmin,icentr)))->Fill(scalprod[0],masses[0]);
      }
      else {
        scalprod[1] = TMath::Cos(fHarmonic*phiD)*Q1[0]+TMath::Sin(fHarmonic*phiD)*Q1[1];
        ((TH2F*)fOutput->FindObject(Form("hMassScalProdu2Q1_pt%dcentr%d_%d",ptbin,icentrmin,icentr)))->Fill(scalprod[1],masses[0]);
      }
    }
    else if(fEvPlaneDet==kFullV0){
      scalprod[0] = TMath::Cos(fHarmonic*phiD)*Q2[0]+TMath::Sin(fHarmonic*phiD)*Q2[1];
      scalprod[1] = TMath::Cos(fHarmonic*phiD)*Q1[0]+TMath::Sin(fHarmonic*phiD)*Q1[1];
      ((TH2F*)fOutput->FindObject(Form("hMassScalProdu3Q1_pt%dcentr%d_%d",ptbin,icentrmin,icentr)))->Fill(scalprod[0],masses[0]);
      ((TH2F*)fOutput->FindObject(Form("hMassScalProdu3Q2_pt%dcentr%d_%d",ptbin,icentrmin,icentr)))->Fill(scalprod[1],masses[0]);
    }
    else if(fEvPlaneDet==kV0A || fEvPlaneDet==kV0C){
      scalprod[0] = TMath::Cos(fHarmonic*phiD)*Q1[0]+TMath::Sin(fHarmonic*phiD)*Q1[1];
      ((TH2F*)fOutput->FindObject(Form("hMassScalProdu3Q1_pt%dcentr%d_%d",ptbin,icentrmin,icentr)))->Fill(scalprod[0],masses[0]);
    }
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
void AliAnalysisTaskSEHFvn::FillDs(AliAODRecoDecayHF* d,TClonesArray *arrayMC,Int_t ptbin,Float_t deltaphi, const Float_t* masses,Int_t isSel,Int_t icentr, Double_t phiD, Double_t etaD, Double_t Q1[2], Double_t Q2[2]){

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
      if(fEvPlaneDet==kFullTPC || fEvPlaneDet==kPosTPC || fEvPlaneDet==kNegTPC) {
        if(etaD>0) {
          scalprod[0] = TMath::Cos(fHarmonic*phiD)*Q2[0]+TMath::Sin(fHarmonic*phiD)*Q2[1];
          ((TH2F*)fOutput->FindObject(Form("hMassScalProdu1Q2_pt%dcentr%d_%d",ptbin,icentrmin,icentr)))->Fill(scalprod[0],masses[0]);
        }
        else {
          scalprod[1] = TMath::Cos(fHarmonic*phiD)*Q1[0]+TMath::Sin(fHarmonic*phiD)*Q1[1];
          ((TH2F*)fOutput->FindObject(Form("hMassScalProdu2Q1_pt%dcentr%d_%d",ptbin,icentrmin,icentr)))->Fill(scalprod[1],masses[0]);
        }
      }
      else if(fEvPlaneDet==kFullV0){
        scalprod[0] = TMath::Cos(fHarmonic*phiD)*Q2[0]+TMath::Sin(fHarmonic*phiD)*Q2[1];
        scalprod[1] = TMath::Cos(fHarmonic*phiD)*Q1[0]+TMath::Sin(fHarmonic*phiD)*Q1[1];
        ((TH2F*)fOutput->FindObject(Form("hMassScalProdu3Q1_pt%dcentr%d_%d",ptbin,icentrmin,icentr)))->Fill(scalprod[0],masses[0]);
        ((TH2F*)fOutput->FindObject(Form("hMassScalProdu3Q2_pt%dcentr%d_%d",ptbin,icentrmin,icentr)))->Fill(scalprod[1],masses[0]);
      }
      else if(fEvPlaneDet==kV0A || fEvPlaneDet==kV0C){
        scalprod[0] = TMath::Cos(fHarmonic*phiD)*Q1[0]+TMath::Sin(fHarmonic*phiD)*Q1[1];
        ((TH2F*)fOutput->FindObject(Form("hMassScalProdu3Q1_pt%dcentr%d_%d",ptbin,icentrmin,icentr)))->Fill(scalprod[0],masses[0]);
      }
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
      if(fEvPlaneDet==kFullTPC || fEvPlaneDet==kPosTPC || fEvPlaneDet==kNegTPC) {
        if(etaD>0) {
          scalprod[0] = TMath::Cos(fHarmonic*phiD)*Q2[0]+TMath::Sin(fHarmonic*phiD)*Q2[1];
          ((TH2F*)fOutput->FindObject(Form("hMassScalProdu1Q2_pt%dcentr%d_%d",ptbin,icentrmin,icentr)))->Fill(scalprod[0],masses[1]);
        }
        else {
          scalprod[1] = TMath::Cos(fHarmonic*phiD)*Q1[0]+TMath::Sin(fHarmonic*phiD)*Q1[1];
          ((TH2F*)fOutput->FindObject(Form("hMassScalProdu2Q1_pt%dcentr%d_%d",ptbin,icentrmin,icentr)))->Fill(scalprod[1],masses[1]);
        }
      }
      else if(fEvPlaneDet==kFullV0){
        scalprod[0] = TMath::Cos(fHarmonic*phiD)*Q2[0]+TMath::Sin(fHarmonic*phiD)*Q2[1];
        scalprod[1] = TMath::Cos(fHarmonic*phiD)*Q1[0]+TMath::Sin(fHarmonic*phiD)*Q1[1];
        ((TH2F*)fOutput->FindObject(Form("hMassScalProdu3Q1_pt%dcentr%d_%d",ptbin,icentrmin,icentr)))->Fill(scalprod[0],masses[1]);
        ((TH2F*)fOutput->FindObject(Form("hMassScalProdu3Q2_pt%dcentr%d_%d",ptbin,icentrmin,icentr)))->Fill(scalprod[1],masses[1]);
      }
      else if(fEvPlaneDet==kV0A || fEvPlaneDet==kV0C){
        scalprod[0] = TMath::Cos(fHarmonic*phiD)*Q1[0]+TMath::Sin(fHarmonic*phiD)*Q1[1];
        ((TH2F*)fOutput->FindObject(Form("hMassScalProdu3Q1_pt%dcentr%d_%d",ptbin,icentrmin,icentr)))->Fill(scalprod[0],masses[1]);
      }
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
void AliAnalysisTaskSEHFvn::SetEventPlaneMethod(Int_t method){
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
Float_t AliAnalysisTaskSEHFvn::GetPhiInRange(Float_t phi){
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
Float_t AliAnalysisTaskSEHFvn::GetEventPlane(AliAODEvent* aod, AliEventplane *pl, Double_t eventplaneqncorrTPC[3], Double_t eventplaneqncorrVZERO[3], Float_t &planereso, Float_t &deltaSubAC, Float_t &deltaSubBC, Int_t &nSubEvents) {
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
      fhEventsInfo->Fill(12);
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
    fhEventsInfo->Fill(12);
    return -9999.;
  }

  if(TMath::Abs(rpangleTPC-rpangleVZERO)>fEventPlanesComp) return -9999.;

  Bool_t set=kFALSE;
  Double_t rpangleeventC=0;
  Double_t rpangleeventB=0;
  Double_t rpangleeventA=0;
  nSubEvents=2;
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
  planereso = TMath::Cos(fHarmonic*deltaPsi); // reaction plane resolution

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
  //if kEvShape histograms filled at the end of UserExec (because of removal of daughters from q2)

  if(nSubEvents==3){
    deltaSubAC=rpangleeventA-rpangleeventC;
    if(TMath::Abs(deltaSubAC)>TMath::Pi()/fHarmonic){
      // difference of subevents reaction plane angle cannot be bigger than phi/2
      if(deltaSubAC>0.) deltaSubAC-=2.*TMath::Pi()/fHarmonic;
      else deltaSubAC +=2.*TMath::Pi()/fHarmonic;
    }
    deltaSubBC=rpangleeventB-rpangleeventC;
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
    //if kEvShape histograms filled at the end of UserExec (because of removal of daughters from q2)
  }
  return eventplane;
}
//________________________________________________________________________
void AliAnalysisTaskSEHFvn::ComputeTPCEventPlane(AliAODEvent* aod, Double_t &rpangleTPC, Double_t &rpangleTPCpos,Double_t &rpangleTPCneg) const {
  Int_t nTracks=aod->GetNumberOfTracks();
  Double_t qVec[2]={0.,0.};
  Double_t qVecPosEta[2]={0.,0.};
  Double_t qVecNegEta[2]={0.,0.};
  for(Int_t it=0; it<nTracks; it++){
    AliAODTrack* track=(AliAODTrack*)aod->GetTrack(it);
    if(!track) continue;
    if(track->TestFilterBit(BIT(8))||track->TestFilterBit(BIT(9))) {
      Double_t eta=track->Eta();
      Double_t pt=track->Pt();
      if((fEtaGapInTPCHalves>0. && TMath::Abs(eta)<fEtaGapInTPCHalves) || eta<fTPCEtaMin || eta>fTPCEtaMax) continue;
      if(pt<0.2 || pt>5) {continue;}
      Double_t phi=track->Phi();
      Double_t wi=1.;
      if(fUsePtWeights){
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
Float_t AliAnalysisTaskSEHFvn::GetEventPlaneForCandidate(AliAODRecoDecayHF* d, AliEventplane *pl){
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

  if(fDecChannel==1 || fDecChannel==4){
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
// Float_t AliAnalysisTaskSEHFvn::GetEventPlaneFromV0(AliAODEvent *aodEvent){
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
Float_t AliAnalysisTaskSEHFvn::GetEventPlaneForCandidateNewQnFw(AliAODRecoDecayHF* d, const TList *qnlist) {

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
    else if(fDecChannel==1 || fDecChannel==4) nProngs = 2; //D0
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
const AliQnCorrectionsQnVector *AliAnalysisTaskSEHFvn::GetQnVectorFromList(const TList *list,
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
Double_t AliAnalysisTaskSEHFvn::Getq2(TList* qnlist, Int_t q2meth, Double_t &mult)
{
  if(!qnlist) {return -1;}

  TString expectedstepTPC="latest";
  TString altstepTPC="plain";
  TString expectedstepV0="latest";
  TString altstepV0="raw";
  if(!fUseQnFrameworkCorrq2) {
    expectedstepTPC="plain";
    altstepTPC="plain";
    expectedstepV0="raw";
    altstepV0="raw";
  }

  const AliQnCorrectionsQnVector* qnVect = 0x0;

  if(q2meth==kq2TPC)
    qnVect = GetQnVectorFromList(qnlist, Form("%sQoverSqrtM",fDetTPCConfName[0].Data()), expectedstepTPC.Data(), altstepTPC.Data());
  else if(q2meth==kq2NegTPC)
    qnVect = GetQnVectorFromList(qnlist, Form("%sQoverSqrtM",fDetTPCConfName[1].Data()), expectedstepTPC.Data(), altstepTPC.Data());
  else if(q2meth==kq2PosTPC)
    qnVect = GetQnVectorFromList(qnlist, Form("%sQoverSqrtM",fDetTPCConfName[2].Data()), expectedstepTPC.Data(), altstepTPC.Data());
  else if(q2meth==kq2VZERO)
    qnVect = GetQnVectorFromList(qnlist, Form("%sQoverSqrtM",fDetV0ConfName[0].Data()), expectedstepV0.Data(), altstepV0.Data());
  else if(q2meth==kq2VZEROA)
    qnVect = GetQnVectorFromList(qnlist, Form("%sQoverSqrtM",fDetV0ConfName[1].Data()), expectedstepV0.Data(), altstepV0.Data());
  else if(q2meth==kq2VZEROC)
    qnVect = GetQnVectorFromList(qnlist, Form("%sQoverSqrtM",fDetV0ConfName[2].Data()), expectedstepV0.Data(), altstepV0.Data());
  else {return -1;}

  if(!qnVect) {return -1;}

  mult = qnVect->GetSumOfWeights();
  Double_t q2 = TMath::Sqrt(qnVect->Qx(2)*qnVect->Qx(2)+qnVect->Qy(2)*qnVect->Qy(2)); //qnVect->Length();
  return q2;
}

//________________________________________________________________________
void AliAnalysisTaskSEHFvn::Setq2Smearing(TString smearingfilepath, TString histoname, Int_t smearingaxis) {
  fq2Smearing=kTRUE;
  fq2SmearingAxis=smearingaxis;

  TFile* smearingfile = TFile::Open(smearingfilepath.Data(),"READ");
  fq2SmearingHisto=(TH2F*)smearingfile->Get(histoname.Data());
  if(fq2SmearingHisto) {fq2SmearingHisto->SetDirectory(0);}
  smearingfile->Close();
}

//________________________________________________________________________
void AliAnalysisTaskSEHFvn::Setq2PercentileSelection(TString splinesfilepath) {
  fPercentileq2=kTRUE;
  
  TString listname[6];
  for(Int_t iDet=0; iDet<3; iDet++) {
    listname[iDet] = "SplineListq2"+fDetTPCConfName[iDet];
    listname[iDet+3] = "SplineListq2"+fDetV0ConfName[iDet];
  }
  
  TFile* splinesfile = TFile::Open(splinesfilepath.Data(),"READ");
  if(!splinesfile) {AliFatal("File with splines for q2 percentiles not found!");}
  for(Int_t iDet=0; iDet<6; iDet++) {
    fq2SplinesList[iDet] = (TList*)splinesfile->Get(listname[iDet].Data());
    if(!fq2SplinesList[iDet]) {AliFatal("TList with splines for q2 percentiles not found in the spline file!");}
    fq2SplinesList[iDet]->SetOwner(0);
  }
  splinesfile->Close();
}

//________________________________________________________________________
Double_t AliAnalysisTaskSEHFvn::ComputeTPCq2(AliAODEvent* aod, Double_t &q2TPCfull, Double_t &q2TPCpos,Double_t &q2TPCneg, Double_t q2VecFullTPC[2], Double_t q2VecPosTPC[2], Double_t q2VecNegTPC[2], Double_t multQvecTPC[3], std::vector<Int_t>& labrejtracks) const {
  /// Compute the q2 for ESE starting from TPC tracks
  /// Option to reject a fraction of tracks to emulate resolution effects

  if(fq2Meth==kq2VZERO || fq2Meth==kq2VZEROA || fq2Meth==kq2VZEROC){
    AliWarning("The recalculation of q2 is implemented only for TPC\n");
    return 0.;
  }

  Int_t nTracks=aod->GetNumberOfTracks();
  Double_t nHarmonic=2.;
  multQvecTPC[0]=0; //full TPC
  multQvecTPC[1]=0; //pos TPC
  multQvecTPC[2]=0; //neg TPC
  Int_t nDau=3;
  if(fDecChannel==1 || fDecChannel==4) {nDau=2;}
  Int_t RandTracks[3]={-1,-1,-1};
  if(fRemoveNdauRandomTracks) {
    if(nTracks>=nDau) {
      for(Int_t iRand=0; iRand<nDau; iRand++) {
        RandTracks[iRand]=gRandom->Integer((UInt_t)(nTracks+1));
      }
    }
    else return 0.; //remove all tracks if they are <= number of daughters
  }
  for(Int_t it=0; it<nTracks; it++){
    if(fRemoveNdauRandomTracks && (it==RandTracks[0] || it==RandTracks[1] || it==RandTracks[2])) {continue;}
    AliAODTrack* track=(AliAODTrack*)aod->GetTrack(it);
    if(!track) continue;
    if(track->TestFilterBit(BIT(8))||track->TestFilterBit(BIT(9))) {
      Int_t lab=track->GetLabel();
      Double_t pt=track->Pt();
      Double_t pseudoRand=pt*1000.-(Long_t)(pt*1000);
      Double_t eta=track->Eta();
      Double_t phi=track->Phi();
      Double_t qx=TMath::Cos(nHarmonic*phi);
      Double_t qy=TMath::Sin(nHarmonic*phi);
      if(eta<fTPCEtaMax && eta>fTPCEtaMin && pt>0.2 && pt<5) {
        if(pseudoRand<fFractionOfTracksForTPCq2){
          q2VecFullTPC[0]+=qx;
          q2VecFullTPC[1]+=qy;
          multQvecTPC[0]++;
        }
        else {
          labrejtracks.push_back(lab);
        }
        if(eta>0){
          q2VecPosTPC[0]+=qx;
          q2VecPosTPC[1]+=qy;
          multQvecTPC[1]++;
        }else{
          q2VecNegTPC[0]+=qx;
          q2VecNegTPC[1]+=qy;
          multQvecTPC[2]++;
        }
      }
    }
  }

  q2TPCfull = 0.;
  if(multQvecTPC[0]>0) q2TPCfull = TMath::Sqrt(q2VecFullTPC[0]*q2VecFullTPC[0]+q2VecFullTPC[1]*q2VecFullTPC[1])/TMath::Sqrt(multQvecTPC[0]);
  q2TPCpos = 0.;
  if(multQvecTPC[1]>0) q2TPCpos = TMath::Sqrt(q2VecPosTPC[0]*q2VecPosTPC[0]+q2VecPosTPC[1]*q2VecPosTPC[1])/TMath::Sqrt(multQvecTPC[1]);
  q2TPCneg = 0.;
  if(multQvecTPC[2]>0) q2TPCneg = TMath::Sqrt(q2VecNegTPC[0]*q2VecNegTPC[0]+q2VecNegTPC[1]*q2VecNegTPC[1])/TMath::Sqrt(multQvecTPC[2]);

  if(fq2Meth==kq2TPC) {return q2TPCfull;}
  else if(fq2Meth==kq2PosTPC) {return q2TPCpos;}
  else if(fq2Meth==kq2NegTPC) {return q2TPCneg;}
  else return 0.;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskSEHFvn::isInMassRange(Double_t massCand, Double_t pt) {
  if(fDecChannel==0) {
    Double_t mass=TDatabasePDG::Instance()->GetParticle(411)->Mass();
    Double_t sigma = 0.01+0.0005*pt; //GeV
    if(massCand>mass-3*sigma && massCand<mass+3*sigma) {return kTRUE;}
  }
  else if(fDecChannel==1 || fDecChannel==4) {
    Double_t mass=TDatabasePDG::Instance()->GetParticle(421)->Mass();
    Double_t sigma = 0.01+0.0005*pt; //GeV
    if(massCand>mass-3*sigma && massCand<mass+3*sigma) {return kTRUE;}
  }
  else if(fDecChannel==2) {
    Double_t deltamass=(TDatabasePDG::Instance()->GetParticle(413)->Mass())-(TDatabasePDG::Instance()->GetParticle(421)->Mass());
    Double_t sigma = 0.0008; //GeV
    if(massCand>deltamass-3*sigma && massCand<deltamass+3*sigma) {return kTRUE;}
  }
  else if(fDecChannel==3) {
    Double_t mass=TDatabasePDG::Instance()->GetParticle(431)->Mass();
    Double_t sigma = 0.01+0.0005*pt; //GeV
    if(massCand>mass-3*sigma && massCand<mass+3*sigma) {return kTRUE;}
  }

  return kFALSE;
}

//________________________________________________________________________
Double_t AliAnalysisTaskSEHFvn::GetTPCq2DauSubQnFramework(Double_t qVectWOcorr[2], Double_t multQvec, Int_t nDauRemoved, Double_t qVecDau[2], Double_t corrRec[2], Double_t LbTwist[2], Bool_t isTwistApplied) {

  if(multQvec<=0) {return 0.;}

  Double_t qVecRemDau[2];
  qVecRemDau[0] = (qVectWOcorr[0]*multQvec - qVecDau[0])/(multQvec-nDauRemoved) - corrRec[0];
  qVecRemDau[1] = (qVectWOcorr[1]*multQvec - qVecDau[1])/(multQvec-nDauRemoved) - corrRec[1];
  if(isTwistApplied) {
    Double_t qRemDauRec[2] = {qVecRemDau[0],qVecRemDau[1]};
    qVecRemDau[0] = (qRemDauRec[0]-LbTwist[0]*qRemDauRec[1])/(1-LbTwist[0]*LbTwist[1]);
    qVecRemDau[1] = (qRemDauRec[1]-LbTwist[1]*qRemDauRec[0])/(1-LbTwist[0]*LbTwist[1]);
  }

  return TMath::Sqrt(qVecRemDau[0]*qVecRemDau[0]+qVecRemDau[1]*qVecRemDau[1]); //already normalised to sqrtM
}

//________________________________________________________________________
void AliAnalysisTaskSEHFvn::RemoveTracksInDeltaEtaFromOnTheFlyTPCq2(AliAODEvent* aod, Double_t etaD, Double_t etaLims[2], Double_t qVec[2], Double_t &M, std::vector<Int_t> daulab) {

  if(2*fDeltaEtaDmesonq2>=etaLims[1]-etaLims[0]) {return;}
  if(etaD<etaLims[0] || etaD>etaLims[1]) {return;}

  Double_t etagaplow=fDeltaEtaDmesonq2;
  Double_t etagaphigh=fDeltaEtaDmesonq2;
  Double_t etagaplowlimit = etaD-fDeltaEtaDmesonq2;
  Double_t etagaphighlimit = etaD+fDeltaEtaDmesonq2;
  //Asymmetric eta gap if meson close to acceptance limit
  if(etagaplowlimit < etaLims[0]) {
    etagaplowlimit = etaLims[0];
    etagaplow = etaD-etagaplowlimit;
    etagaphighlimit += (fDeltaEtaDmesonq2-etagaplow);
  }
  if(etagaphighlimit > etaLims[1]) {
    etagaphighlimit = etaLims[1];
    etagaphigh = etagaphighlimit-etaD;
    etagaplowlimit -= (fDeltaEtaDmesonq2-etagaphigh);
  }

  //loop on tracks and remove those inside eta gap between D meson and q2
  Int_t nTracks=aod->GetNumberOfTracks();
  for(Int_t iTrack=0; iTrack<nTracks; iTrack++) {
    AliAODTrack* track = (AliAODTrack*)aod->GetTrack(iTrack);
    if(!track || (!track->TestFilterBit(BIT(8)) && !track->TestFilterBit(BIT(9)))) {continue;}
    Int_t lab = track->GetLabel();
    Bool_t isdau=kFALSE;
    for(UInt_t iDau=0; iDau<daulab.size(); iDau++) {
      if(lab==daulab[iDau]) {isdau=kTRUE; break;}
    }
    if(isdau) {continue;} //do not remove daughters twice!
    Double_t eta = track->Eta();
    Double_t phi = track->Phi();
    Double_t pt = track->Pt();
    Double_t qx=TMath::Cos(2*phi);
    Double_t qy=TMath::Sin(2*phi);
    if((pt>0.2 && pt<5) && (eta>etagaplowlimit && eta<etagaphighlimit)) {//if is in right eta and pt regioxn, remove
      qVec[0] -= qx;
      qVec[1] -= qy;
      M--;
    }
  }
}

//________________________________________________________________________
Bool_t AliAnalysisTaskSEHFvn::IsSoftPionSelected(AliAODTrack* track)
{
  //track cuts
  if(fUseFiltBit4SoftPion) {
    if(!track->TestFilterBit(4)) return kFALSE;
  }
  else {
    //applying ESDtrackCut
    if(!fCutsSoftPion->IsSelected(track)) return kFALSE;
  }
  
  return kTRUE;
}

//________________________________________________________________________
void AliAnalysisTaskSEHFvn::Terminate(Option_t */*option*/)
{
  // Terminate analysis
  //
  if(fDebug > 1) printf("AnalysisTaskSEHFvn: Terminate() \n");

  return;
}
