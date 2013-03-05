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
// Class AliAnalysisTaskSEDvsMultiplicity
// AliAnalysisTaskSE for the D meson vs. multiplcity analysis
// Authors: Renu Bala, Zaida Conesa del Valle, Francesco Prino
/////////////////////////////////////////////////////////////

#include <TClonesArray.h>
#include <TCanvas.h>
#include <TList.h>
#include <TString.h>
#include <TDatabasePDG.h>
#include <TH1F.h>
#include <TH2F.h>     
#include <TH3F.h>
#include <THnSparse.h>
#include <TProfile.h>
#include "AliAnalysisManager.h"
#include "AliRDHFCuts.h"
#include "AliRDHFCutsDplustoKpipi.h"
#include "AliRDHFCutsDStartoKpipi.h"
#include "AliRDHFCutsD0toKpi.h"
#include "AliAODHandler.h"
#include "AliAODEvent.h"
#include "AliAODVertex.h"
#include "AliAODTrack.h"
#include "AliAODRecoDecayHF.h"
#include "AliAODRecoCascadeHF.h"
#include "AliAnalysisVertexingHF.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisTaskSEDvsMultiplicity.h"
#include "AliNormalizationCounter.h"
#include "AliVertexingHFUtils.h"
ClassImp(AliAnalysisTaskSEDvsMultiplicity)


//________________________________________________________________________
AliAnalysisTaskSEDvsMultiplicity::AliAnalysisTaskSEDvsMultiplicity():
AliAnalysisTaskSE(),
  fOutput(0),
  fListCuts(0),
  fOutputCounters(0),
  fListProfiles(0),
  fHistNEvents(0),
  fHistNtrEta16vsNtrEta1(0),
  fHistNtrCorrEta1vsNtrRawEta1(0),
  fHistNtrVsZvtx(0),
  fHistNtrCorrVsZvtx(0),
  fHistNtrVsNchMC(0),
  fHistNtrCorrVsNchMC(0),
  fHistNtrVsNchMCPrimary(0),
  fHistNtrCorrVsNchMCPrimary(0),
  fHistNtrVsNchMCPhysicalPrimary(0),
  fHistNtrCorrVsNchMCPhysicalPrimary(0),
  fHistGenPrimaryParticlesInelGt0(0),
  fHistNchMCVsNchMCPrimaryVsNchMCPhysicalPrimary(0),
  fHistNtrUnCorrEvSel(0),
  fHistNtrCorrEvSel(0),
  fHistNtrCorrEvWithCand(0),
  fHistNtrCorrEvWithD(0),
  fPtVsMassVsMult(0),
  fPtVsMassVsMultNoPid(0),
  fPtVsMassVsMultUncorr(0),
  fPtVsMassVsMultPart(0),
  fPtVsMassVsMultAntiPart(0),
  fUpmasslimit(1.965),
  fLowmasslimit(1.765),
  fNMassBins(200),
  fRDCutsAnalysis(0),
  fCounter(0),
  fCounterU(0),
  fDoImpPar(kFALSE),
  fNImpParBins(400),
  fLowerImpPar(-2000.),
  fHigherImpPar(2000.),
  fReadMC(kFALSE),
  fMCOption(0),
  fUseBit(kTRUE),
  fSubtractTrackletsFromDau(kFALSE),
  fRefMult(9.26),
  fPdgMeson(411)
{
   // Default constructor
  for(Int_t i=0; i<5; i++) fHistMassPtImpPar[i]=0;
  for(Int_t i=0; i<4; i++) fMultEstimatorAvg[i]=0;
}

//________________________________________________________________________
AliAnalysisTaskSEDvsMultiplicity::AliAnalysisTaskSEDvsMultiplicity(const char *name, Int_t pdgMeson,AliRDHFCuts *cuts):
  AliAnalysisTaskSE(name),
  fOutput(0),
  fListCuts(0),
  fOutputCounters(0),
  fListProfiles(0),
  fHistNEvents(0),
  fHistNtrEta16vsNtrEta1(0),
  fHistNtrCorrEta1vsNtrRawEta1(0),
  fHistNtrVsZvtx(0),
  fHistNtrCorrVsZvtx(0),
  fHistNtrVsNchMC(0),
  fHistNtrCorrVsNchMC(0),
  fHistNtrVsNchMCPrimary(0),
  fHistNtrCorrVsNchMCPrimary(0),
  fHistNtrVsNchMCPhysicalPrimary(0),
  fHistNtrCorrVsNchMCPhysicalPrimary(0),
  fHistGenPrimaryParticlesInelGt0(0),
  fHistNchMCVsNchMCPrimaryVsNchMCPhysicalPrimary(0),
  fHistNtrUnCorrEvSel(0),
  fHistNtrCorrEvSel(0),
  fHistNtrCorrEvWithCand(0),
  fHistNtrCorrEvWithD(0),
  fPtVsMassVsMult(0),
  fPtVsMassVsMultNoPid(0),
  fPtVsMassVsMultUncorr(0),
  fPtVsMassVsMultPart(0),
  fPtVsMassVsMultAntiPart(0),
  fUpmasslimit(1.965),
  fLowmasslimit(1.765),
  fNMassBins(200),
  fRDCutsAnalysis(cuts),
  fCounter(0),
  fCounterU(0),
  fDoImpPar(kFALSE),
  fNImpParBins(400),
  fLowerImpPar(-2000.),
  fHigherImpPar(2000.),
  fReadMC(kFALSE),
  fMCOption(0),
  fUseBit(kTRUE),
  fSubtractTrackletsFromDau(kFALSE),
  fRefMult(9.26),
  fPdgMeson(pdgMeson)
{
  // 
  // Standard constructor
  //
  for(Int_t i=0; i<5; i++) fHistMassPtImpPar[i]=0;
  for(Int_t i=0; i<4; i++) fMultEstimatorAvg[i]=0;
  if(fPdgMeson==413){
    fNMassBins=200; // FIXME
    SetMassLimits(0.,0.2); // FIXME
  }else{ 
    fNMassBins=200; 
    SetMassLimits(fPdgMeson,0.1);
  }
  // Default constructor
   // Output slot #1 writes into a TList container
  DefineOutput(1,TList::Class());  //My private output
  // Output slot #2 writes cut to private output
  DefineOutput(2,TList::Class());
  // Output slot #3 writes cut to private output
  DefineOutput(3,TList::Class()); 
  // Output slot #4 writes cut to private output
  DefineOutput(4,TList::Class()); 
}
//________________________________________________________________________
AliAnalysisTaskSEDvsMultiplicity::~AliAnalysisTaskSEDvsMultiplicity()
{
  //
  // Destructor
  //
  delete fOutput;
  delete fHistNEvents;
  delete fListCuts;
  delete fListProfiles;
  delete fRDCutsAnalysis;
  delete fCounter;
  delete fCounterU;
  for(Int_t i=0; i<4; i++) delete fMultEstimatorAvg[i];
  for(Int_t i=0; i<5; i++){
    delete fHistMassPtImpPar[i];
  }
}  

//_________________________________________________________________
void  AliAnalysisTaskSEDvsMultiplicity::SetMassLimits(Double_t lowlimit, Double_t uplimit){
  // set invariant mass limits
  if(uplimit>lowlimit){
    fLowmasslimit = lowlimit;
    fUpmasslimit = uplimit;
  }else{
    AliError("Wrong mass limits: upper value should be larger than lower one");
  }
}
//_________________________________________________________________
void  AliAnalysisTaskSEDvsMultiplicity::SetMassLimits(Int_t pdg, Double_t range){
  // set invariant mass limits
  Double_t mass=TDatabasePDG::Instance()->GetParticle(TMath::Abs(pdg))->Mass();
  SetMassLimits(mass-range,mass+range);
}
//________________________________________________________________________
void AliAnalysisTaskSEDvsMultiplicity::Init(){
  //
  // Initialization
  //
  printf("AnalysisTaskSEDvsMultiplicity::Init() \n");
  
  fListCuts=new TList();

  if(fPdgMeson==411){
     AliRDHFCutsDplustoKpipi* copycut=new AliRDHFCutsDplustoKpipi(*(static_cast<AliRDHFCutsDplustoKpipi*>(fRDCutsAnalysis)));
    copycut->SetName("AnalysisCutsDplus");
    fListCuts->Add(copycut);
  }else if(fPdgMeson==421){
    AliRDHFCutsD0toKpi* copycut=new AliRDHFCutsD0toKpi(*(static_cast<AliRDHFCutsD0toKpi*>(fRDCutsAnalysis)));
    copycut->SetName("AnalysisCutsDzero");
    fListCuts->Add(copycut);
  }else if(fPdgMeson==413){
     AliRDHFCutsDStartoKpipi* copycut=new AliRDHFCutsDStartoKpipi(*(static_cast<AliRDHFCutsDStartoKpipi*>(fRDCutsAnalysis)));
    copycut->SetName("AnalysisCutsDStar");
    fListCuts->Add(copycut);
  }
  PostData(2,fListCuts);
  
  fListProfiles = new TList();
  fListProfiles->SetOwner();
  TString period[4]={"LHC10b","LHC10c","LHC10d","LHC10e"};
  for(Int_t i=0; i<4; i++){
    if(fMultEstimatorAvg[i]){
      TProfile* hprof=new TProfile(*fMultEstimatorAvg[i]);
      hprof->SetName(Form("ProfileTrkVsZvtx%s\n",period[i].Data()));
      fListProfiles->Add(hprof);
    }
  }
  PostData(4,fListProfiles);

  return;
}

//________________________________________________________________________
void AliAnalysisTaskSEDvsMultiplicity::UserCreateOutputObjects()
{
  // Create the output container
  //
  if(fDebug > 1) printf("AnalysisTaskSEDvsMultiplicity::UserCreateOutputObjects() \n");

  // Several histograms are more conveniently managed in a TList
  fOutput = new TList();
  fOutput->SetOwner();
  fOutput->SetName("OutputHistos");

  fHistNtrUnCorrEvSel = new TH1F("hNtrUnCorrEvSel","Uncorrected tracklets multiplicity for selected events; Tracklets ; Entries",200,-0.5,199.5);
  fHistNtrCorrEvSel = new TH1F("hNtrCorrEvSel","Corrected tracklets multiplicity for selected events; Tracklets ; Entries",200,-0.5,199.5);
  fHistNtrCorrEvWithCand = new TH1F("hNtrCorrEvWithCand", "Tracklets multiplicity for events with D candidates; Tracklets ; Entries",200,-0.5,199.5);// Total multiplicity
  fHistNtrCorrEvWithD = new TH1F("hNtrCorrEvWithD", "Tracklets multiplicity for events with D in mass region ; Tracklets ; Entries",200,-0.5,199.5); // 
  fHistNtrEta16vsNtrEta1 = new TH2F("hNtrEta16vsNtrEta1","Uncorrected Eta1.6 vs Eta1.0; Ntracklets #eta<1.0; Ntracklets #eta<1.6",200,-0.5,199.5,200,-0.5,199.5); //eta 1.6 vs eta 1.0 histogram 
  fHistNtrCorrEta1vsNtrRawEta1 = new TH2F("hNtrCorrEta1vsNtrRawEta1","Corrected Eta1 vs Eta1.0; Ntracklets #eta<1.0 corrected; Ntracklets #eta<1",200,-0.5,199.5,200,-0.5,199.5); //eta 1.6 vs eta 1.0 histogram 
  fHistNtrVsZvtx = new TH2F("hNtrVsZvtx","Ntracklet vs VtxZ; VtxZ;N_{tracklet};",300,-15,15,200,-0.5,199.5); // 
  fHistNtrCorrVsZvtx = new TH2F("hNtrCorrVsZvtx","Ntracklet vs VtxZ; VtxZ;N_{tracklet};",300,-15,15,200,-0.5,199.5); // 

  fHistNtrVsNchMC = new TH2F("hNtrVsNchMC","Ntracklet vs NchMC; Nch;N_{tracklet};",200,-0.5,199.5,200,-0.5,199.5); // 
  fHistNtrCorrVsNchMC = new TH2F("hNtrCorrVsNchMC","Ntracklet vs Nch; Nch;N_{tracklet};",200,-0.5,199.5,200,-0.5,199.5); // 
  
  fHistNtrVsNchMCPrimary = new TH2F("hNtrVsNchMCPrimary","Ntracklet vs Nch (Primary); Nch (Primary);N_{tracklet};",200,-0.5,199.5,200,-0.5,199.5); // 
  fHistNtrCorrVsNchMCPrimary = new TH2F("hNtrCorrVsNchMCPrimary","Ntracklet vs Nch (Primary); Nch(Primary) ;N_{tracklet};",200,-0.5,199.5,200,-0.5,199.5); // 

  fHistNtrVsNchMCPhysicalPrimary = new TH2F("hNtrVsNchMCPhysicalPrimary","Ntracklet vs Nch (Physical Primary); Nch (Physical Primary);N_{tracklet};",200,-0.5,199.5,200,-0.5,199.5); // 
  fHistNtrCorrVsNchMCPhysicalPrimary = new TH2F("hNtrCorrVsMCPhysicalPrimary","Ntracklet vs Nch (Physical Primary); Nch (Physical Primary);N_{tracklet};",200,-0.5,199.5,200,-0.5,199.5); // 
  
  fHistGenPrimaryParticlesInelGt0 = new TH1F("hGenPrimaryParticlesInelGt0","Multiplcity of generated charged particles ; Nparticles ; Entries",200,-0.5,199.5);

  fHistNchMCVsNchMCPrimaryVsNchMCPhysicalPrimary = new TH3F("fHistNchMCVsNchMCPrimaryVsNchMCPhysicalPrimary", "MC: Nch (Physical Primary) vs Nch (Primary) vs Nch (Generated); Nch (Generated); Nch (Primary); Nch (Physical Primary)",200,-0.5,199.5,200,-0.5,199.5,200,-0.5,199.5);

  fHistNtrUnCorrEvSel->Sumw2();
  fHistNtrCorrEvSel->Sumw2();
  fHistNtrCorrEvWithCand->Sumw2();
  fHistNtrCorrEvWithD->Sumw2();
  fHistGenPrimaryParticlesInelGt0->Sumw2();
  fOutput->Add(fHistNtrUnCorrEvSel);
  fOutput->Add(fHistNtrCorrEvSel);
  fOutput->Add(fHistNtrCorrEvWithCand);
  fOutput->Add(fHistNtrCorrEvWithD);
  fOutput->Add(fHistNtrEta16vsNtrEta1);
  fOutput->Add(fHistNtrCorrEta1vsNtrRawEta1);
  fOutput->Add(fHistNtrVsZvtx);
  fOutput->Add(fHistNtrCorrVsZvtx);

  fOutput->Add(fHistNtrVsNchMC);
  fOutput->Add(fHistNtrCorrVsNchMC);
  fOutput->Add(fHistNtrVsNchMCPrimary);
  fOutput->Add(fHistNtrCorrVsNchMCPrimary);
  fOutput->Add(fHistNtrVsNchMCPhysicalPrimary);
  fOutput->Add(fHistNtrCorrVsNchMCPhysicalPrimary);
  fOutput->Add(fHistGenPrimaryParticlesInelGt0);
  fOutput->Add(fHistNchMCVsNchMCPrimaryVsNchMCPhysicalPrimary);

  
  fHistNEvents = new TH1F("fHistNEvents", "number of events ",11,-0.5,10.5);
  fHistNEvents->GetXaxis()->SetBinLabel(1,"nEvents total");
  fHistNEvents->GetXaxis()->SetBinLabel(2,"nEvents with Z vertex");
  fHistNEvents->GetXaxis()->SetBinLabel(3,"nEvents selected");
  fHistNEvents->GetXaxis()->SetBinLabel(4,"Rejected due to trigger");
  fHistNEvents->GetXaxis()->SetBinLabel(5,"Rejected due to phys sel");
  fHistNEvents->GetXaxis()->SetBinLabel(6,"Rejected due to vertex cuts");
  fHistNEvents->GetXaxis()->SetBinLabel(7,"Rejected due to pileup");
  fHistNEvents->GetXaxis()->SetBinLabel(8,"Total no. of candidate");
  fHistNEvents->GetXaxis()->SetBinLabel(9,"no. of cand wo bitmask");
  fHistNEvents->GetXaxis()->SetBinLabel(10,"D after cuts (No PID)");
  fHistNEvents->GetXaxis()->SetBinLabel(11,"D after cuts + PID)"); 
  fHistNEvents->GetXaxis()->SetNdivisions(1,kFALSE);  
  fHistNEvents->Sumw2();
  fHistNEvents->SetMinimum(0);
  fOutput->Add(fHistNEvents);

  fPtVsMassVsMult=new TH3F("hPtVsMassvsMult", "D candidates: p_{t} vs mass vs tracklets multiplicity; Tracklets; Mass M [GeV/c^{2}]; p_{t} [GeV/c]",200,-0.5,199.5,fNMassBins,fLowmasslimit,fUpmasslimit,48,0.,24.);
 
  fPtVsMassVsMultNoPid=new TH3F("hPtVsMassvsMultNoPid", "D candidates: p_{t} vs mass vs tracklets multiplicity; Tracklets; Mass M [GeV/c^{2}]; p_{t} [GeV/c]",200,-0.5,199.5,fNMassBins,fLowmasslimit,fUpmasslimit,48,0.,24.); 

  fPtVsMassVsMultUncorr=new TH3F("hPtVsMassvsMultUncorr", "D candidates: p_{t} vs mass vs tracklets multiplicity; Tracklets; Mass M [GeV/c^{2}]; p_{t} [GeV/c]",200,-0.5,199.5,fNMassBins,fLowmasslimit,fUpmasslimit,48,0.,24.);

  fPtVsMassVsMultPart=new TH3F("hPtVsMassvsMultPart", "D candidates: p_{t} vs mass vs tracklets multiplicity; Tracklets; Mass M [GeV/c^{2}]; p_{t} [GeV/c]",200,-0.5,199.5,fNMassBins,fLowmasslimit,fUpmasslimit,48,0.,24.);

  fPtVsMassVsMultAntiPart=new TH3F("hPtVsMassvsMultAntiPart", "D candidates: p_{t} vs mass vs tracklets multiplicity; Tracklets; Mass M [GeV/c^{2}]; p_{t} [GeV/c]",200,-0.5,199.5,fNMassBins,fLowmasslimit,fUpmasslimit,48,0.,24.);

  fOutput->Add(fPtVsMassVsMult);
  fOutput->Add(fPtVsMassVsMultUncorr);
  fOutput->Add(fPtVsMassVsMultNoPid);
  fOutput->Add(fPtVsMassVsMultPart);
  fOutput->Add(fPtVsMassVsMultAntiPart);

  if(fDoImpPar) CreateImpactParameterHistos();

  fCounter = new AliNormalizationCounter("NormCounterCorrMult");
  fCounter->SetStudyMultiplicity(kTRUE,1.);
  fCounter->Init(); 

  fCounterU = new AliNormalizationCounter("NormCounterUnCorrMult");
  fCounterU->SetStudyMultiplicity(kTRUE,1.);
  fCounterU->Init(); 

  fOutputCounters = new TList();
  fOutputCounters->SetOwner();
  fOutputCounters->SetName("OutputCounters");
  fOutputCounters->Add(fCounter);
  fOutputCounters->Add(fCounterU);
  
  PostData(1,fOutput); 
  PostData(2,fListCuts);
  PostData(3,fOutputCounters);
  PostData(4,fListProfiles);

  

  return;
}



//________________________________________________________________________
void AliAnalysisTaskSEDvsMultiplicity::UserExec(Option_t */*option*/)
{
  // Execute analysis for current event:
  // heavy flavor candidates association to MC truth

  AliAODEvent *aod = dynamic_cast<AliAODEvent*> (InputEvent());
  
  //  AliAODTracklets* tracklets = aod->GetTracklets();
  //Int_t ntracklets = tracklets->GetNumberOfTracklets();
 
  
  TClonesArray *arrayCand = 0;
  TString arrayName="";
  UInt_t pdgDau[3];
  Int_t nDau=0;
  Int_t selbit=0;
  if(fPdgMeson==411){
    arrayName="Charm3Prong";
    pdgDau[0]=211; pdgDau[1]=321; pdgDau[2]=211; 
    nDau=3;
    selbit=AliRDHFCuts::kDplusCuts;
  }else if(fPdgMeson==421){
    arrayName="D0toKpi";
    pdgDau[0]=211; pdgDau[1]=321; pdgDau[2]=0;
    nDau=2;
    selbit=AliRDHFCuts::kD0toKpiCuts;
  }else if(fPdgMeson==413){
    arrayName="Dstar";
    pdgDau[0]=321; pdgDau[1]=211; pdgDau[2]=211; 
    nDau=3;
    selbit=AliRDHFCuts::kDstarCuts;
  }

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
      arrayCand=(TClonesArray*)aodFromExt->GetList()->FindObject(arrayName.Data());
    }
  } else if(aod) {
    arrayCand=(TClonesArray*)aod->GetList()->FindObject(arrayName.Data());
  }

  if(!aod || !arrayCand) {
    printf("AliAnalysisTaskSEDvsMultiplicity::UserExec: Charm3Prong branch not found!\n");
    return;
  }



  // fix for temporary bug in ESDfilter 
  // the AODs with null vertex pointer didn't pass the PhysSel
  if(!aod->GetPrimaryVertex()||TMath::Abs(aod->GetMagneticField())<0.001) return;

  Int_t countTreta1=AliVertexingHFUtils::GetNumberOfTrackletsInEtaRange(aod,-1.,1.);
  Int_t countTr=AliVertexingHFUtils::GetNumberOfTrackletsInEtaRange(aod,-1.6,1.6);

  fCounterU->StoreEvent(aod,fRDCutsAnalysis,fReadMC,countTreta1);
  fHistNEvents->Fill(0); // count event

  Double_t countTreta1corr=countTreta1;
  AliAODVertex *vtx1 = (AliAODVertex*)aod->GetPrimaryVertex();
  if(vtx1){
    if(vtx1->GetNContributors()>0){    
      fHistNEvents->Fill(1); 
      TProfile* estimatorAvg = GetEstimatorHistogram(aod);
      if(estimatorAvg){
	countTreta1corr=AliVertexingHFUtils::GetCorrectedNtracklets(estimatorAvg,countTreta1,vtx1->GetZ(),fRefMult); 
      }
    }
  }
   
  fCounter->StoreEvent(aod,fRDCutsAnalysis,fReadMC,(Int_t)countTreta1corr);

  Bool_t isEvSel=fRDCutsAnalysis->IsEventSelected(aod);

  if(fRDCutsAnalysis->GetWhyRejection()==5) fHistNEvents->Fill(3);
  if(fRDCutsAnalysis->GetWhyRejection()==7) fHistNEvents->Fill(4); 
  if(fRDCutsAnalysis->GetWhyRejection()==6) fHistNEvents->Fill(5);
  if(fRDCutsAnalysis->GetWhyRejection()==1) fHistNEvents->Fill(6);
  
  
  if(!isEvSel)return;
  fHistNtrEta16vsNtrEta1->Fill(countTreta1,countTr);
  fHistNtrCorrEta1vsNtrRawEta1->Fill(countTreta1,countTreta1corr);
  if(vtx1){
    fHistNtrVsZvtx->Fill(vtx1->GetZ(),countTreta1);
    fHistNtrCorrVsZvtx->Fill(vtx1->GetZ(),countTreta1corr);
  }

  TClonesArray *arrayMC=0;
  AliAODMCHeader *mcHeader=0;

  // load MC particles
  if(fReadMC){
     
    arrayMC =  (TClonesArray*)aod->GetList()->FindObject(AliAODMCParticle::StdBranchName());
    if(!arrayMC) {
      printf("AliAnalysisTaskSEDvsMultiplicity::UserExec: MC particles branch not found!\n");
      return;
    }  
    // load MC header
    mcHeader =  (AliAODMCHeader*)aod->GetList()->FindObject(AliAODMCHeader::StdBranchName());
    if(!mcHeader) {
      printf("AliAnalysisTaskSEDvsMultiplicity::UserExec: MC header branch not found!\n");
      return;
     }
  

    Int_t nChargedMC=AliVertexingHFUtils::GetGeneratedMultiplicityInEtaRange(arrayMC,-1.0,1.0);
    Int_t nChargedMCPrimary=AliVertexingHFUtils::GetGeneratedPrimariesInEtaRange(arrayMC,-1.0,1.0);
    Int_t nChargedMCPhysicalPrimary=AliVertexingHFUtils::GetGeneratedPhysicalPrimariesInEtaRange(arrayMC,-1.0,1.0);
    if(nChargedMCPhysicalPrimary>0){ // INEL>0 for |eta|<1
      fHistGenPrimaryParticlesInelGt0->Fill(nChargedMCPhysicalPrimary);
    }
    fHistNtrVsNchMC->Fill(nChargedMC,countTreta1);
    fHistNtrCorrVsNchMC->Fill(nChargedMC,countTreta1corr);

    fHistNtrVsNchMCPrimary->Fill(nChargedMCPrimary,countTreta1);
    fHistNtrCorrVsNchMCPrimary->Fill(nChargedMCPrimary,countTreta1corr);

    fHistNtrVsNchMCPhysicalPrimary->Fill(nChargedMCPhysicalPrimary,countTreta1);
    fHistNtrCorrVsNchMCPhysicalPrimary->Fill(nChargedMCPhysicalPrimary,countTreta1corr);

    fHistNchMCVsNchMCPrimaryVsNchMCPhysicalPrimary->Fill(nChargedMC,nChargedMCPrimary,nChargedMCPhysicalPrimary);
  }
  
  Int_t nCand = arrayCand->GetEntriesFast(); 
  Int_t nSelectedNoPID=0,nSelectedPID=0,nSelectedInMassPeak=0;
  Double_t mD0PDG = TDatabasePDG::Instance()->GetParticle(421)->Mass();
  Double_t mDplusPDG = TDatabasePDG::Instance()->GetParticle(411)->Mass();
  Double_t mDstarPDG = TDatabasePDG::Instance()->GetParticle(413)->Mass();

  for (Int_t iCand = 0; iCand < nCand; iCand++) {
    AliAODRecoDecayHF *d = (AliAODRecoDecayHF*)arrayCand->UncheckedAt(iCand);
    fHistNEvents->Fill(7);
    if(fUseBit && !d->HasSelectionBit(selbit)){
      fHistNEvents->Fill(8);
      continue;
    }
    
    Double_t ptCand = d->Pt();
    Double_t rapid=d->Y(fPdgMeson);
    Bool_t isFidAcc=fRDCutsAnalysis->IsInFiducialAcceptance(ptCand,rapid);
    if(!isFidAcc) continue;
    Int_t passAllCuts=fRDCutsAnalysis->IsSelected(d,AliRDHFCuts::kAll,aod);
    Int_t passTopolCuts=fRDCutsAnalysis->GetIsSelectedCuts();
    if(passTopolCuts==0) continue;
    nSelectedNoPID++;
    fHistNEvents->Fill(9);
    if(passAllCuts){
      nSelectedPID++;
      fHistNEvents->Fill(10);
    }
    Double_t multForCand=countTreta1corr;
    if(fSubtractTrackletsFromDau){
      for(Int_t iDau=0; iDau<nDau; iDau++){
	AliAODTrack *t = (AliAODTrack*)d->GetDaughter(iDau);
	if(!t) continue;
	if(t->HasPointOnITSLayer(0) && t->HasPointOnITSLayer(1)){
	  if(multForCand>0) multForCand-=1;
	}
      }
    }
    Bool_t isPrimary=kTRUE;
    Int_t labD=-1;
    Double_t trueImpParXY=9999.;
    Double_t impparXY=d->ImpParXY()*10000.;
    Double_t dlen=0.1; //FIXME
    Double_t mass[2]={-1.,-1.};
    if(fPdgMeson==411){
      mass[0]=d->InvMass(nDau,pdgDau);
      mass[1]=-1.;
      if(TMath::Abs(mass[0]-mDplusPDG)<0.02) nSelectedInMassPeak++; //20 MeV for now... FIXME
    }else if(fPdgMeson==421){
      UInt_t pdgdaughtersD0[2]={211,321};//pi,K 
      UInt_t pdgdaughtersD0bar[2]={321,211};//K,pi 
      mass[0]=d->InvMass(2,pdgdaughtersD0);
      mass[1]=d->InvMass(2,pdgdaughtersD0bar);
      if(TMath::Abs(mass[0]-mD0PDG)<0.02 || TMath::Abs(mass[1]-mD0PDG)<0.02 ) nSelectedInMassPeak++; //20 MeV for now... FIXME
    }else if(fPdgMeson==413){
      // FIXME
      AliAODRecoCascadeHF* temp = (AliAODRecoCascadeHF*)d;
      mass[0]=temp->DeltaInvMass();
      mass[1]=-1.;
      if(TMath::Abs(mass[0]-(mDstarPDG-mD0PDG))<0.001) nSelectedInMassPeak++; //1 MeV for now... FIXME
    }
    for(Int_t iHyp=0; iHyp<2; iHyp++){
      if(mass[iHyp]<0.) continue; // for D+ and D* we have 1 mass hypothesis
      Double_t invMass=mass[iHyp];
      Double_t arrayForSparse[5]={invMass,ptCand,impparXY,dlen,multForCand};

      if(fReadMC){
	
	labD = d->MatchToMC(fPdgMeson,arrayMC,nDau,(Int_t*)pdgDau);
	Bool_t fillHisto=fDoImpPar;
	if(labD>=0){
	  AliAODMCParticle *partD = (AliAODMCParticle*)arrayMC->At(labD);
	  Int_t code=partD->GetPdgCode();
	  if(CheckOrigin(arrayMC,partD)==5) isPrimary=kFALSE;
	  if(code<0 && iHyp==0) fillHisto=kFALSE;
	  if(code>0 && iHyp==1) fillHisto=kFALSE;
	  if(!isPrimary){
	    if(fPdgMeson==411){
	      trueImpParXY=AliVertexingHFUtils::GetTrueImpactParameterDplus(mcHeader,arrayMC,partD)*10000.;
	    }else if(fPdgMeson==421){
	      trueImpParXY=AliVertexingHFUtils::GetTrueImpactParameterDzero(mcHeader,arrayMC,partD)*10000.;
	    }else if(fPdgMeson==413){
	      trueImpParXY=0.; /// FIXME
	    }
	    Double_t arrayForSparseTrue[5]={invMass,ptCand,trueImpParXY,dlen,multForCand};
	    if(fillHisto && passAllCuts){
	      fHistMassPtImpPar[2]->Fill(arrayForSparse);
	      fHistMassPtImpPar[3]->Fill(arrayForSparseTrue);
	    }
	  }else{
	    if(fillHisto && passAllCuts) fHistMassPtImpPar[1]->Fill(arrayForSparse);
	  }
	}else{
	  if(fillHisto && passAllCuts)fHistMassPtImpPar[4]->Fill(arrayForSparse);
	}
	if(fPdgMeson==421){
	  if(TMath::Abs(labD)==fPdgMeson && fMCOption==2) continue;
	  if(TMath::Abs(labD)!=fPdgMeson && fMCOption==1) continue;      
	}
      }
      
      if(fPdgMeson==421){
	if(iHyp==0 && !(passTopolCuts&1)) continue; // candidate not passing as D0
	if(iHyp==1 && !(passTopolCuts&2)) continue; // candidate not passing as D0bar
      }

      fPtVsMassVsMultNoPid->Fill(multForCand,invMass,ptCand);

      if(fPdgMeson==421){
	if(iHyp==0 && !(passAllCuts&1)) continue; // candidate not passing as D0
	if(iHyp==1 && !(passAllCuts&2)) continue; // candidate not passing as D0bar
      }
      if(passAllCuts){
	fPtVsMassVsMult->Fill(multForCand,invMass,ptCand);	   	      
	fPtVsMassVsMultUncorr->Fill(countTreta1,invMass,ptCand);
	// Add separation between part antipart
	if(fPdgMeson==411){
	  if(d->GetCharge()>0) fPtVsMassVsMultPart->Fill(multForCand,invMass,ptCand);
	  else fPtVsMassVsMultAntiPart->Fill(multForCand,invMass,ptCand);
	}else if(fPdgMeson==421){
	  if(passAllCuts&1) fPtVsMassVsMultPart->Fill(multForCand,invMass,ptCand);
	  if(passAllCuts&2) fPtVsMassVsMultAntiPart->Fill(multForCand,invMass,ptCand);
	}else if(fPdgMeson==413){
	  // FIXME ADD Dstar!!!!!!!!
	}
      	
	if(fDoImpPar){
	  fHistMassPtImpPar[0]->Fill(arrayForSparse);
	}
	
      }

    }
  }
  fCounter->StoreCandidates(aod,nSelectedNoPID,kTRUE);
  fCounter->StoreCandidates(aod,nSelectedPID,kFALSE);
  fHistNtrUnCorrEvSel->Fill(countTreta1);
  fHistNtrCorrEvSel->Fill(countTreta1corr);
  if(nSelectedPID>0) fHistNtrCorrEvWithCand->Fill(countTreta1corr);
  if(nSelectedInMassPeak>0) fHistNtrCorrEvWithD->Fill(countTreta1corr);

  PostData(1,fOutput); 
  PostData(2,fListCuts);
  PostData(3,fOutputCounters);
    
  return;
}
//________________________________________________________________________
void AliAnalysisTaskSEDvsMultiplicity::CreateImpactParameterHistos(){
  // Histos for impact paramter study
  // mass . pt , impact parameter , decay length , multiplicity

  Int_t nbins[5]={fNMassBins,200,fNImpParBins,50,100};
  Double_t xmin[5]={fLowmasslimit,0.,fLowerImpPar,0.,0.};
  Double_t xmax[5]={fUpmasslimit,20.,fHigherImpPar,1.,100.};

  fHistMassPtImpPar[0]=new THnSparseF("hMassPtImpParAll",
					"Mass vs. pt vs.imppar - All",
					5,nbins,xmin,xmax);
  fHistMassPtImpPar[1]=new THnSparseF("hMassPtImpParPrompt",
					"Mass vs. pt vs.imppar - promptD",
					5,nbins,xmin,xmax);
  fHistMassPtImpPar[2]=new THnSparseF("hMassPtImpParBfeed",
					"Mass vs. pt vs.imppar - DfromB",
					5,nbins,xmin,xmax);
  fHistMassPtImpPar[3]=new THnSparseF("hMassPtImpParTrueBfeed",
					"Mass vs. pt vs.true imppar -DfromB",
					5,nbins,xmin,xmax);
  fHistMassPtImpPar[4]=new THnSparseF("hMassPtImpParBkg",
				        "Mass vs. pt vs.imppar - backgr.",
					5,nbins,xmin,xmax);
  for(Int_t i=0; i<5;i++){
    fOutput->Add(fHistMassPtImpPar[i]);
  }
}

//________________________________________________________________________
void AliAnalysisTaskSEDvsMultiplicity::Terminate(Option_t */*option*/)
{
  // Terminate analysis
  //
  if(fDebug > 1) printf("AnalysisTaskSEDvsMultiplicity: Terminate() \n");

  fOutput = dynamic_cast<TList*> (GetOutputData(1));
  if (!fOutput) {     
    printf("ERROR: fOutput not available\n");
    return;
  }

  fHistNEvents = dynamic_cast<TH1F*>(fOutput->FindObject("fHistNEvents"));
  if(!fHistNEvents){
    printf("ERROR: fHistNEvents not available\n");
    return;    
  }
  printf("Number of Analyzed Events = %d\n",(Int_t)fHistNEvents->GetBinContent(3));
 
  return;
}
//_________________________________________________________________________________________________
Int_t AliAnalysisTaskSEDvsMultiplicity::CheckOrigin(TClonesArray* arrayMC, AliAODMCParticle *mcPartCandidate) const {		
  //
  // checking whether the mother of the particles come from a charm or a bottom quark
  //
	
  Int_t pdgGranma = 0;
  Int_t mother = 0;
  mother = mcPartCandidate->GetMother();
  Int_t istep = 0;
  Int_t abspdgGranma =0;
  Bool_t isFromB=kFALSE;
  Bool_t isQuarkFound=kFALSE;
  while (mother >0 ){
    istep++;
    AliAODMCParticle* mcGranma = dynamic_cast<AliAODMCParticle*>(arrayMC->At(mother));
    if (mcGranma){
      pdgGranma = mcGranma->GetPdgCode();
      abspdgGranma = TMath::Abs(pdgGranma);
      if ((abspdgGranma > 500 && abspdgGranma < 600) || (abspdgGranma > 5000 && abspdgGranma < 6000)){
	isFromB=kTRUE;
      }
      if(abspdgGranma==4 || abspdgGranma==5) isQuarkFound=kTRUE;
      mother = mcGranma->GetMother();
    }else{
      AliError("Failed casting the mother particle!");
      break;
    }
  }
  
  if(isFromB) return 5;
  else return 4;
}



//____________________________________________________________________________
TProfile* AliAnalysisTaskSEDvsMultiplicity::GetEstimatorHistogram(const AliVEvent* event){
  // Get Estimator Histogram from period event->GetRunNumber();
  //
  // If you select SPD tracklets in |eta|<1 you should use type == 1
  //

  Int_t runNo  = event->GetRunNumber();
  Int_t period = -1;   // 0-LHC10b, 1-LHC10c, 2-LHC10d, 3-LHC10e
  if(runNo>114930 && runNo<117223) period = 0;
  if(runNo>119158 && runNo<120830) period = 1;
  if(runNo>122373 && runNo<126438) period = 2;
  if(runNo>127711 && runNo<130841) period = 3;
  if(period<0 || period>3) return 0;

  return fMultEstimatorAvg[period];
}
