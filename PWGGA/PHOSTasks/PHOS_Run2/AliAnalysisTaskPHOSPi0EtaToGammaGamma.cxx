#include "TChain.h"
#include "TTree.h"
#include "TObjArray.h"
#include "TClonesArray.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TProfile.h"
#include "TVector2.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TParticle.h"
#include "THnSparse.h"
#include "THashList.h"
#include "TMath.h"

#include "AliAnalysisManager.h"
#include "AliAnalysisTaskSE.h"
#include "AliLog.h"

#include "AliCaloPhoton.h"
#include "AliPHOSGeometry.h"

#include "AliESDtrackCuts.h"
#include "AliESDHeader.h"
#include "AliESDEvent.h"
#include "AliESDCaloCells.h"
#include "AliESDCaloCluster.h"
#include "AliESDVertex.h"
#include "AliESDtrack.h"

#include "AliVEvent.h"
#include "AliVHeader.h"
#include "AliVTrack.h"
#include "AliVCluster.h"
#include "AliVCaloCells.h"

#include "AliMultSelection.h"
#include "AliEventplane.h"
#include "AliQnCorrectionsManager.h"
#include "AliAnalysisTaskFlowVectorCorrections.h"
#include "AliQnCorrectionsQnVector.h"

#include "AliOADBContainer.h"

#include "AliPID.h"
#include "AliPIDResponse.h"

#include "AliAODMCHeader.h"
#include "AliAODHeader.h"
#include "AliAODEvent.h"
#include "AliAODCaloCells.h"
#include "AliAODCaloCluster.h"
#include "AliAODVertex.h"
#include "AliAODTrack.h"
#include "AliAODMCParticle.h"
#include "AliAODInputHandler.h"

#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliGenEventHeader.h"
#include "AliGenPythiaEventHeader.h"
#include "AliGenHijingEventHeader.h"
#include "AliGenCocktailEventHeader.h"

#include "AliPHOSTenderTask.h"
#include "AliPHOSTenderSupply.h"
#include "AliMagF.h"
#include "TGeoGlobalMagField.h"

#include "AliPHOSEventCuts.h"
#include "AliPHOSClusterCuts.h"
#include "AliPHOSJetJetMC.h"
#include "AliPHOSTriggerHelper.h"
#include "AliAnalysisTaskPHOSPi0EtaToGammaGamma.h"

// Author: Daiki Sekihata (Hiroshima University)

ClassImp(AliAnalysisTaskPHOSPi0EtaToGammaGamma)

//________________________________________________________________________
AliAnalysisTaskPHOSPi0EtaToGammaGamma::AliAnalysisTaskPHOSPi0EtaToGammaGamma(const char *name):
  AliAnalysisTaskSE(name),
	fIsMC(kFALSE),
  fIsJJMC(kFALSE),
  fPtHardBin(-1),
  fMCType("MBMC"),
	fUsePHOSTender(kFALSE),
  fUseCoreEnergy(kFALSE),
  fPHOSEventCuts(0x0),
  fPHOSClusterCuts(0x0),
  fBunchSpace(25.),
  fCollisionSystem(-1),
  fTOFEfficiency(0x0),
  fTriggerEfficiency(0x0),
  fESDtrackCutsGlobal(0x0),
  fESDtrackCutsGlobalConstrained(0x0),
  fCentArrayPi0(0x0),
  fCentArrayEta(0x0),
  fCentArrayGamma(0x0),
  fCentArrayK0S(0x0),
  fCentArrayL0(0x0),

  fOutputContainer(0x0),
  fEvent(0x0),
  fESDEvent(0x0),
  fAODEvent(0x0),
  fMCArrayESD(0x0),
  fMCArrayAOD(0x0),
  fJJMCHandler(0x0),
  fPtHardAndJetPtFactor(-1),
  fPtHardAndSinglePtFactor(-1),
  fRunNumber(0),
  fPHOSGeo(0x0),
  fPHOSClusterArray(NULL),
  fEstimator("V0M"),
  fMultSelection(0x0),
  fCentralityMain(0),
  fCentralityMin(0.),
  fCentralityMax(90.),
  fNMixed(10),
  fZvtx(-1),
  fEPBin(-1),
  fIsFlowTask(kFALSE),
  fHarmonics(-1),
  fQNormalization("QoverM"),
  fFM(-1),
  fQnDetectorMain(AliAnalysisTaskPHOSPi0EtaToGammaGamma::kFullV0),
  fQnDetectorSub1(AliAnalysisTaskPHOSPi0EtaToGammaGamma::kTPCPosEta),
  fQnDetectorSub2(AliAnalysisTaskPHOSPi0EtaToGammaGamma::kTPCNegEta),
  fFlowQnVectorMgr(0x0),
  fEventPlane(-999.),
  fQVector1(0,0),
  fQnStep(5),
  fNHybridTrack(0),
  fIsPHOSTriggerAnalysis(kFALSE),
  fEnergyThreshold(0.),
  fPHOSTriggerHelper(0x0),
  fPHOSTriggerHelperL0(0x0),
  fPHOSTriggerHelperL1H(0x0),
  fPHOSTriggerHelperL1M(0x0),
  fPHOSTriggerHelperL1L(0x0),
  fForceActiveTRU(kFALSE),
  fTRFM(-1),
  fPIDResponse(0x0),
  fIsNonLinStudy(kFALSE),
  fGlobalEScale(1.0),
  fEmin(0.2),
  fIsOAStudy(kFALSE),
  fNMixTrack(5),
  fMatchingR(2.),
  fAnaOmega3Pi(kFALSE),
  fMinPtPi0(0),
  fMinPtChPi(0),
  fMaxR(999.),
  fPIDStudy(kFALSE)
{
  // Constructor

  for(Int_t i=0;i<10;i++){
    for(Int_t j=0;j<12;j++){
      fPHOSEvents[i][j] = 0x0;
    }
  }

  for(Int_t i=0;i<3;i++){
    fVertex[i] = 0;
  }

  for(Int_t i=0;i<7;i++){
    for(Int_t j=0;j<7;j++){
      fNonLin[i][j] = 0x0;
    }
  }

  for(Int_t i=0;i<11;i++){
    fAdditionalPi0PtWeight[i]   = 0x0;
    fAdditionalEtaPtWeight[i]   = 0x0;
    fAdditionalGammaPtWeight[i] = 0x0;
    fAdditionalK0SPtWeight[i]   = 0x0;
    fAdditionalL0PtWeight[i]    = 0x0;
  }

  fTPCEPName[0] = "TPC";
  fTPCEPName[1] = "TPCNegEta";
  fTPCEPName[2] = "TPCPosEta";

  fV0EPName[0] = "VZERO";
  fV0EPName[1] = "VZEROA";
  fV0EPName[2] = "VZEROC";

  for(Int_t i=0;i<11;i++){
    fAdditionalPi0PtWeight[i]   = new TF1(Form("fAdditionalPi0PtWeight_%d",i)  ,"1.",0,100);
    fAdditionalEtaPtWeight[i]   = new TF1(Form("fAdditionalEtaPtWeight_%d",i)  ,"1.",0,100);
    fAdditionalGammaPtWeight[i] = new TF1(Form("fAdditionalGammaPtWeight_%d",i),"1.",0,100);
    fAdditionalK0SPtWeight[i]   = new TF1(Form("fAdditionalK0SPtWeight_%d",i)  ,"1.",0,100);
    fAdditionalL0PtWeight[i]    = new TF1(Form("fAdditionalL0PtWeight_%d",i)   ,"1.",0,100);
  }

  fTOFEfficiency = new TF1("fTOFEfficiency","1.",0,100);
  fTriggerEfficiency = new TF1("fTriggerEfficiency","1.",0,100);

  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #0 id reserved by the base class for AOD
  // Output slot #1 writes into a TH1 container
  DefineOutput(1, THashList::Class());

}
//________________________________________________________________________
AliAnalysisTaskPHOSPi0EtaToGammaGamma::~AliAnalysisTaskPHOSPi0EtaToGammaGamma()
{
  for(Int_t i=0;i<10;i++){
    for(Int_t j=0;j<12;j++){

      if(fPHOSEvents[i][j]){
        delete fPHOSEvents[i][j];
        fPHOSEvents[i][j] = 0x0;
      }

    }
  }

  if(fPHOSTriggerHelper){
    delete fPHOSTriggerHelper;
    fPHOSTriggerHelper  = 0x0;
  }

  if(fPHOSTriggerHelperL0 ){
    delete fPHOSTriggerHelperL0;
    fPHOSTriggerHelperL0  = 0x0;
  }
  if(fPHOSTriggerHelperL1H){
    delete fPHOSTriggerHelperL1H;
    fPHOSTriggerHelperL1H = 0x0;
  }

  if(fPHOSTriggerHelperL1M){
    delete fPHOSTriggerHelperL1M;
    fPHOSTriggerHelperL1M = 0x0;
  }

  if(fPHOSTriggerHelperL1L){
    delete fPHOSTriggerHelperL1L;
    fPHOSTriggerHelperL1L = 0x0;
  }

  for(Int_t i=0;i<11;i++){
    if(fAdditionalPi0PtWeight[i]){
      delete fAdditionalPi0PtWeight[i];
      fAdditionalPi0PtWeight[i] = 0x0;
    }

    if(fAdditionalK0SPtWeight[i]){
      delete fAdditionalK0SPtWeight[i];
      fAdditionalK0SPtWeight[i] = 0x0;
    }

    if(fAdditionalEtaPtWeight[i]){
      delete fAdditionalEtaPtWeight[i];
      fAdditionalEtaPtWeight[i] = 0x0;
    }

    if(fAdditionalGammaPtWeight[i]){
      delete fAdditionalGammaPtWeight[i];
      fAdditionalGammaPtWeight[i] = 0x0;
    }

  }

  if(fCentArrayPi0){
    delete fCentArrayPi0;
    fCentArrayPi0 = 0x0;
  }
  if(fCentArrayEta){
    delete fCentArrayEta;
    fCentArrayEta = 0x0;
  }
  if(fCentArrayGamma){
    delete fCentArrayGamma;
    fCentArrayGamma = 0x0;
  }
  if(fCentArrayK0S){
    delete fCentArrayK0S;
    fCentArrayK0S = 0x0;
  }
  if(fCentArrayL0){
    delete fCentArrayL0;
    fCentArrayL0 = 0x0;
  }

  delete fTOFEfficiency;
  fTOFEfficiency = 0x0;

  delete fTriggerEfficiency;
  fTriggerEfficiency = 0x0;

  if(fPHOSEventCuts){
    delete fPHOSEventCuts;
    fPHOSEventCuts = 0x0;
  }

  if(fPHOSClusterCuts){
    delete fPHOSClusterCuts;
    fPHOSClusterCuts = 0x0;
  }

  for(Int_t i=0;i<7;i++){
    for(Int_t j=0;j<7;j++){
      if(fNonLin[i][j]){
        delete fNonLin[i][j];
        fNonLin[i][j] = 0x0;
      }
    }
  }

}
//________________________________________________________________________
void AliAnalysisTaskPHOSPi0EtaToGammaGamma::UserCreateOutputObjects()
{
  // Create histograms
  // Called once

  fOutputContainer = new THashList();
  fOutputContainer->SetOwner(kTRUE);

  TH1F *hEventSummary = new TH1F("hEventSummary","Event Summary",10,0.5,10.5);
  hEventSummary->GetXaxis()->SetBinLabel(1 ,"all");
  hEventSummary->GetXaxis()->SetBinLabel(2 ,"selected");
  hEventSummary->GetXaxis()->SetBinLabel(3 ,"0PH0 fired");
  hEventSummary->GetXaxis()->SetBinLabel(4 ,"1PHL fired");
  hEventSummary->GetXaxis()->SetBinLabel(5 ,"1PHM fired");
  hEventSummary->GetXaxis()->SetBinLabel(6 ,"1PHH fired");
  hEventSummary->GetXaxis()->SetBinLabel(7 ,"0PH0 fired & matched");
  hEventSummary->GetXaxis()->SetBinLabel(8 ,"1PHL fired & matched");
  hEventSummary->GetXaxis()->SetBinLabel(9 ,"1PHM fired & matched");
  hEventSummary->GetXaxis()->SetBinLabel(10,"1PHH fired & matched");
  fOutputContainer->Add(hEventSummary);

  TH1F *hPHI7Summary = new TH1F("hPHI7Summary","PHI7 Summary",4,-1.5,2.5);
  hPHI7Summary->GetXaxis()->SetBinLabel(1,"L0");//-1
  hPHI7Summary->GetXaxis()->SetBinLabel(2,"1PHH");//0
  hPHI7Summary->GetXaxis()->SetBinLabel(3,"1PHM");//1
  hPHI7Summary->GetXaxis()->SetBinLabel(4,"1PHL");//2
  fOutputContainer->Add(hPHI7Summary);

  //event character histogram
  fOutputContainer->Add(new TH1F("hVertexZ","VertexZ;Zvtx (cm)",100,-50.,50.));
  fOutputContainer->Add(new TH1F("hVertexZSelectEvent","VertexZ SelectEvent;Zvtx (cm)",100,-50.,50.));

  fOutputContainer->Add(new TH2F(Form("hCentrality%svsNContributor",fEstimator.Data()),Form("Centrality %s vs. Ncontributor;centrality (%%);N_{contributor}",fEstimator.Data()),100,0.,100,2001,-0.5,2000.5));

  fOutputContainer->Add(new TH2F("hCentralityV0MvsCL0","Centrality V0M vs. CL0;V0M;CL0",100,0.,100,100,0.,100.));
  fOutputContainer->Add(new TH2F("hCentralityV0MvsCL1","Centrality V0M vs. CL1;V0M;CL1",100,0.,100,100,0.,100.));
  fOutputContainer->Add(new TH2F("hCentralityCL0vsCL1","Centrality CL0 vs. CL1;CL0;CL1",100,0.,100,100,0.,100.));
  fOutputContainer->Add(new TH2F("hCentralityV0AvsV0C","Centrality V0A vs. V0C;V0A;V0C",100,0.,100,100,0.,100.));
  fOutputContainer->Add(new TH2F("hCentralityZNAvsZNC","Centrality ZNA vs. ZNC;ZNA;ZNC",100,0.,100,100,0.,100.));

  const Double_t Pi = TMath::Pi();
  const Double_t TwoPi = TMath::TwoPi();

  for(Int_t i=0;i<3;i++){
    fOutputContainer->Add(new TH2F(Form("hCentrality%svsEventPlane%s%s",fEstimator.Data(),fTPCEPName[i].Data(),fQNormalization.Data()),Form("Centrality %s vs. EP %s %s;centrality (%%);#Psi_{EP}",fEstimator.Data(),fTPCEPName[i].Data(),fQNormalization.Data()),100,0,100,30,0,Pi));
  }
  for(Int_t i=0;i<3;i++){
    fOutputContainer->Add(new TH2F(Form("hCentrality%svsEventPlane%s%s",fEstimator.Data(),fV0EPName[i].Data(),fQNormalization.Data()),Form("Centrality %s vs. EP %s %s;centrality (%%);#Psi_{EP}",fEstimator.Data(),fV0EPName[i].Data(),fQNormalization.Data()),100,0,100,30,0,Pi));
  }
  Double_t SPlimit = 0.2;
  Int_t NbinSP     = 200;
  if(fQnStep == 5){//In Qn correction task, rescaling of Q vector is applied for V0 after step4
    SPlimit = 15.;
    NbinSP  = 600;
  }
  fOutputContainer->Add(new TH2F(Form("hCentrality%svsSPQ1Q2",fEstimator.Data()),Form("Centrality %s vs. SP #vec{Q_{1}} #upoint #vec{Q_{2}};centrality (%%);SP #vec{Q_{1}} #upoint #vec{Q_{2}}",fEstimator.Data()),100,0,100,NbinSP,-SPlimit,SPlimit));
  fOutputContainer->Add(new TH2F(Form("hCentrality%svsSPQ2Q3",fEstimator.Data()),Form("Centrality %s vs. SP #vec{Q_{2}} #upoint #vec{Q_{3}};centrality (%%);SP #vec{Q_{2}} #upoint #vec{Q_{3}}",fEstimator.Data()),100,0,100,NbinSP,-SPlimit,SPlimit));
  fOutputContainer->Add(new TH2F(Form("hCentrality%svsSPQ3Q1",fEstimator.Data()),Form("Centrality %s vs. SP #vec{Q_{3}} #upoint #vec{Q_{1}};centrality (%%);SP #vec{Q_{3}} #upoint #vec{Q_{1}}",fEstimator.Data()),100,0,100,NbinSP,-SPlimit,SPlimit));

  Double_t Qxylimit = 1;
  Int_t NbinQxy     = 200;
  if(fQnStep == 5){//In Qn correction task, rescaling of Q vector is applied for V0 after step4
    Qxylimit = 10.;
    NbinQxy  = 200;
  }

  fOutputContainer->Add(new TH2F(Form("hCentrality%svsNormQ1",fEstimator.Data()),Form("Centrality %s vs. |Q_{1}|;centrality (%%);|Q_{1}|",fEstimator.Data()),100,0,100,NbinQxy,0,2. * Qxylimit));
  fOutputContainer->Add(new TH2F(Form("hCentrality%svsNormQ2",fEstimator.Data()),Form("Centrality %s vs. |Q_{2}|;centrality (%%);|Q_{2}|",fEstimator.Data()),100,0,100,NbinQxy,0,2. * Qxylimit));
  fOutputContainer->Add(new TH2F(Form("hCentrality%svsNormQ3",fEstimator.Data()),Form("Centrality %s vs. |Q_{3}|;centrality (%%);|Q_{3}|",fEstimator.Data()),100,0,100,NbinQxy,0,2. * Qxylimit));

  fOutputContainer->Add(new TH2F(Form("hCentrality%svsQ1x",fEstimator.Data()),Form("Centrality %s vs. Q_{1x};centrality (%%);Q_{1x}",fEstimator.Data()),100,0,100,NbinQxy,-Qxylimit,Qxylimit));
  fOutputContainer->Add(new TH2F(Form("hCentrality%svsQ1y",fEstimator.Data()),Form("Centrality %s vs. Q_{1y};centrality (%%);Q_{1y}",fEstimator.Data()),100,0,100,NbinQxy,-Qxylimit,Qxylimit));
  fOutputContainer->Add(new TH2F(Form("hCentrality%svsQ2x",fEstimator.Data()),Form("Centrality %s vs. Q_{2x};centrality (%%);Q_{2x}",fEstimator.Data()),100,0,100,NbinQxy,-Qxylimit,Qxylimit));
  fOutputContainer->Add(new TH2F(Form("hCentrality%svsQ2y",fEstimator.Data()),Form("Centrality %s vs. Q_{2y};centrality (%%);Q_{2y}",fEstimator.Data()),100,0,100,NbinQxy,-Qxylimit,Qxylimit));
  fOutputContainer->Add(new TH2F(Form("hCentrality%svsQ3x",fEstimator.Data()),Form("Centrality %s vs. Q_{3x};centrality (%%);Q_{3x}",fEstimator.Data()),100,0,100,NbinQxy,-Qxylimit,Qxylimit));
  fOutputContainer->Add(new TH2F(Form("hCentrality%svsQ3y",fEstimator.Data()),Form("Centrality %s vs. Q_{3y};centrality (%%);Q_{3y}",fEstimator.Data()),100,0,100,NbinQxy,-Qxylimit,Qxylimit));

  fOutputContainer->Add(new TH2F(Form("hCentrality%svsCosDeltaEventPlane12",fEstimator.Data()),Form("Centrality %s vs. cos(%d #Delta#Psi_{EP});centrality (%%);cos(%d #Delta#Psi_{EP})",fEstimator.Data(),fHarmonics,fHarmonics),100,0,100,20,-1,1));
  fOutputContainer->Add(new TH2F(Form("hCentrality%svsCosDeltaEventPlane23",fEstimator.Data()),Form("Centrality %s vs. cos(%d #Delta#Psi_{EP});centrality (%%);cos(%d #Delta#Psi_{EP})",fEstimator.Data(),fHarmonics,fHarmonics),100,0,100,20,-1,1));
  fOutputContainer->Add(new TH2F(Form("hCentrality%svsCosDeltaEventPlane31",fEstimator.Data()),Form("Centrality %s vs. cos(%d #Delta#Psi_{EP});centrality (%%);cos(%d #Delta#Psi_{EP})",fEstimator.Data(),fHarmonics,fHarmonics),100,0,100,20,-1,1));

  fOutputContainer->Add(new TH2F(Form("hCentrality%svsSinDeltaEventPlane12",fEstimator.Data()),Form("Centrality %s vs. sin(%d #Delta#Psi_{EP});centrality (%%);sin(%d #Delta#Psi_{EP})",fEstimator.Data(),fHarmonics,fHarmonics),100,0,100,20,-1,1));
  fOutputContainer->Add(new TH2F(Form("hCentrality%svsSinDeltaEventPlane23",fEstimator.Data()),Form("Centrality %s vs. sin(%d #Delta#Psi_{EP});centrality (%%);sin(%d #Delta#Psi_{EP})",fEstimator.Data(),fHarmonics,fHarmonics),100,0,100,20,-1,1));
  fOutputContainer->Add(new TH2F(Form("hCentrality%svsSinDeltaEventPlane31",fEstimator.Data()),Form("Centrality %s vs. sin(%d #Delta#Psi_{EP});centrality (%%);sin(%d #Delta#Psi_{EP})",fEstimator.Data(),fHarmonics,fHarmonics),100,0,100,20,-1,1));

  fOutputContainer->Add(new TH2F(Form("hCentrality%svsPHOSClusterMultiplicityMC" ,fEstimator.Data()),Form("Centrality %s vs. Cluster Multiplicity;centrality (%%);Ncluster"             ,fEstimator.Data()),100,0,100,201,-0.5,200.5));
  fOutputContainer->Add(new TH2F(Form("hCentrality%svsPHOSClusterMultiplicity"   ,fEstimator.Data()),Form("Centrality %s vs. Cluster Multiplicity;centrality (%%);Ncluster"             ,fEstimator.Data()),100,0,100,201,-0.5,200.5));
  fOutputContainer->Add(new TH2F(Form("hCentrality%svsPHOSClusterMultiplicityTOF",fEstimator.Data()),Form("Centrality %s vs. Cluster Multiplicity with TOF cut;centrality (%%);Ncluster",fEstimator.Data()),100,0,100,201,-0.5,200.5));

  fOutputContainer->Add(new TH2F(Form("hCentrality%svsSPDTracklet",fEstimator.Data()),Form("Centrality %s vs. SPD tracklet;centrality (%%);SPD tracklets",fEstimator.Data()),100,0,100,600,0,6000));
  fOutputContainer->Add(new TH2F(Form("hCentrality%svsTrackMultiplicity",fEstimator.Data()),Form("Centrality %s vs. track Multiplicity;centrality (%%);track multiplicity",fEstimator.Data()),100,0,100,400,0,4000));
  fOutputContainer->Add(new TH2F("hPHOSClusterMultiplicityvsTrackMultiplicity"   ,"cluster multiplicity vs. track multiplicity;track multiplicity;Ncluster"         ,400,0,4000,201,-0.5,200.5));
  fOutputContainer->Add(new TH2F("hPHOSClusterMultiplicityTOFvsTrackMultiplicity","cluster multiplicity with TOF vs. track multiplicity;track multiplicity;Ncluster",400,0,4000,201,-0.5,200.5));

  fOutputContainer->Add(new TH2F("hMultiplicityV0AvsTrackMultiplicity","V0A multiplicity vs. Track Multiplicity;V0A multiplicity;Ntrack"                   ,80,0,40000,400,0,4000));
  fOutputContainer->Add(new TH2F("hMultiplicityV0CvsTrackMultiplicity","V0C multiplicity vs. Track Multiplicity;V0C multiplicity;Ntrack"                   ,80,0,40000,400,0,4000));
  fOutputContainer->Add(new TH2F("hMultiplicityV0vsTrackMultiplicity" ,"V0  multiplicity vs. Track Multiplicity;V0A + V0C multiplicity;Ntrack"             ,80,0,40000,400,0,4000));

  fOutputContainer->Add(new TH2F("hMultiplicityV0AvsPHOSClusterMultiplicity"   ,"V0C multiplicity vs. Cluster Multiplicity;V0A multiplicity;Ncluster"                   ,80,0,40000,201,-0.5,200.5));
  fOutputContainer->Add(new TH2F("hMultiplicityV0CvsPHOSClusterMultiplicity"   ,"V0C multiplicity vs. Cluster Multiplicity;V0C multiplicity;Ncluster"                   ,80,0,40000,201,-0.5,200.5));
  fOutputContainer->Add(new TH2F("hMultiplicityV0vsPHOSClusterMultiplicity"    ,"V0  multiplicity vs. Cluster Multiplicity;V0A + V0C multiplicity;Ncluster"             ,80,0,40000,201,-0.5,200.5));
  fOutputContainer->Add(new TH2F("hMultiplicityV0AvsPHOSClusterMultiplicityTOF","V0C multiplicity vs. Cluster Multiplicity with TOF cut;V0A multiplicity;Ncluster"      ,80,0,40000,201,-0.5,200.5));
  fOutputContainer->Add(new TH2F("hMultiplicityV0CvsPHOSClusterMultiplicityTOF","V0C multiplicity vs. Cluster Multiplicity with TOF cut;V0C multiplicity;Ncluster"      ,80,0,40000,201,-0.5,200.5));
  fOutputContainer->Add(new TH2F("hMultiplicityV0vsPHOSClusterMultiplicityTOF" ,"V0  multiplicity vs. Cluster Multiplicity with TOF cut;V0A + V0C multiplicity;Ncluster",80,0,40000,201,-0.5,200.5));

  //track QA histograms
  const Int_t Ntype=2;
  const TString tracktype[Ntype] = {"Global","GlobalNoDCA"};

  for(Int_t itype=0;itype<Ntype;itype++){
    fOutputContainer->Add(new TH1F(Form("h%sTrackMult"   ,tracktype[itype].Data()),Form("Number of %s track;track multiplicity"                         ,tracktype[itype].Data()),400,0,4000));
    fOutputContainer->Add(new TH1F(Form("h%sTrackPt"     ,tracktype[itype].Data()),Form("%s track p_{T};p_{T} (GeV/c)"                                  ,tracktype[itype].Data()),100,0,100));
    fOutputContainer->Add(new TH2F(Form("h%sTrackEtaPhi" ,tracktype[itype].Data()),Form("%s track #eta vs. #phi;#phi;#eta"                              ,tracktype[itype].Data()),100,0,TwoPi,40,-1,1));
    fOutputContainer->Add(new TH2F(Form("h%sTrackDCA"    ,tracktype[itype].Data()),Form("%s track DCA;DCA_{xy} (cm);DCA_{z} (cm)"                       ,tracktype[itype].Data()),100,-5,5,100,-5,5));
    fOutputContainer->Add(new TH2F(Form("h%sTrackTPCdEdx",tracktype[itype].Data()),Form("%s TPC dE/dx vs. track momentum;p^{track} (GeV/c);dE/dx (a.u.)",tracktype[itype].Data()),200,0,20,200,0,200));
  }

  //cell QA histograms
  const Int_t Nmod=5;
  const Int_t Ntru=8;
  for(Int_t imod=1;imod<Nmod;imod++) fOutputContainer->Add(new TH2F(Form("hCellNXZM%d",imod),Form("Cell N(X,Z) M%d",imod),64,0.5,64.5,56,0.5,56.5));
  for(Int_t imod=1;imod<Nmod;imod++) fOutputContainer->Add(new TH2F(Form("hCellEXZM%d",imod),Form("Cell E(X,Z) M%d",imod),64,0.5,64.5,56,0.5,56.5));
  for(Int_t imod=1;imod<Nmod;imod++) fOutputContainer->Add(new TH2F(Form("hCellAmpTimeM%d_LG",imod),Form("Cell Amplitude vs. Time LG M%d",imod),200,0,20,1000,-500,500));
  for(Int_t imod=1;imod<Nmod;imod++) fOutputContainer->Add(new TH2F(Form("hCellAmpTimeM%d_HG",imod),Form("Cell Amplitude vs. Time HG M%d",imod),200,0,20,1000,-500,500));
  for(Int_t imod=1;imod<Nmod;imod++) fOutputContainer->Add(new TH1F(Form("hCellMultEventM%d",imod),Form("PHOS cell multiplicity per event M%d",imod),1001,-0.5,1000.5));

  //cluster QA histograms
  for(Int_t imod=1;imod<Nmod;imod++) fOutputContainer->Add(new TH1F(Form("hPHOSClusterMultM%d",imod)   ,Form("PHOS cluster multiplicity M%d",imod)             ,101,-0.5,100.5));
  for(Int_t imod=1;imod<Nmod;imod++) fOutputContainer->Add(new TH1F(Form("hPHOSClusterMultTOFM%d",imod),Form("PHOS cluster multiplicity M%d with TOF cut",imod),101,-0.5,100.5));
  for(Int_t imod=1;imod<Nmod;imod++) fOutputContainer->Add(new TH1F(Form("hPHOSClusterEnergyM%d",imod)   ,Form("PHOS cluster energy M%d",imod)             ,500,0.,100));
  for(Int_t imod=1;imod<Nmod;imod++) fOutputContainer->Add(new TH1F(Form("hPHOSClusterEnergyTOFM%d",imod),Form("PHOS cluster energy M%d with TOF cut",imod),500,0.,100));
  for(Int_t imod=1;imod<Nmod;imod++) fOutputContainer->Add(new TH2F(Form("hCluNXZM%d",imod)    ,Form("Cluster N(X,Z) M%d",imod)     ,64,0.5,64.5, 56,0.5,56.5));
  for(Int_t imod=1;imod<Nmod;imod++) fOutputContainer->Add(new TH2F(Form("hCluEXZM%d",imod)    ,Form("Cluster E(X,Z) M%d",imod)     ,64,0.5,64.5, 56,0.5,56.5));
  for(Int_t imod=1;imod<Nmod;imod++) fOutputContainer->Add(new TH2F(Form("hCluLowNXZM%d",imod) ,Form("Cluster Low N(X,Z) M%d",imod) ,64,0.5,64.5, 56,0.5,56.5));
  for(Int_t imod=1;imod<Nmod;imod++) fOutputContainer->Add(new TH2F(Form("hCluLowEXZM%d",imod) ,Form("Cluster Low E(X,Z) M%d",imod) ,64,0.5,64.5, 56,0.5,56.5));
  for(Int_t imod=1;imod<Nmod;imod++) fOutputContainer->Add(new TH2F(Form("hCluHighNXZM%d",imod),Form("Cluster High N(X,Z) M%d",imod),64,0.5,64.5, 56,0.5,56.5));
  for(Int_t imod=1;imod<Nmod;imod++) fOutputContainer->Add(new TH2F(Form("hCluHighEXZM%d",imod),Form("Cluster High E(X,Z) M%d",imod),64,0.5,64.5, 56,0.5,56.5));
  for(Int_t imod=1;imod<Nmod;imod++) fOutputContainer->Add(new TH2F(Form("hCluNXZTOFM%d",imod) ,Form("Cluster N(X,Z) M%d TOF",imod) ,64,0.5,64.5, 56,0.5,56.5));
  for(Int_t imod=1;imod<Nmod;imod++) fOutputContainer->Add(new TH2F(Form("hCluEXZTOFM%d",imod) ,Form("Cluster E(X,Z) M%d TOF",imod) ,64,0.5,64.5, 56,0.5,56.5));

  for(Int_t imod=1;imod<Nmod;imod++) fOutputContainer->Add(new TH2F(Form("hClusterEvsTM%d",imod),Form("Cluster E vs TOF M%d;E (GeV);TOF (ns)",imod)     ,500,0,50, 1000,-500,500));

  fOutputContainer->Add(new TH2F("hClusterEvsN","Cluster E vs N_{cell};E (GeV);N_{cell}",500,0,50,100,0.5,100.5));
  fOutputContainer->Add(new TH2F("hClusterEvsM02","Cluster E vs M02;E (GeV);M02 (cm)",500,0,50,50,0,5));
  fOutputContainer->Add(new TH2F("hClusterEvsM20","Cluster E vs M20;E (GeV);M20 (cm)",500,0,50,50,0,5));
  fOutputContainer->Add(new TH2F("hFullDispvsFullE","full dispersion vs full E;E (GeV);dispersion (#sigma)",100,0,50,100,0,10));
  fOutputContainer->Add(new TH2F("hCoreDispvsCoreE","core dispersion vs core E;E (GeV);dispersion (#sigma)",100,0,50,100,0,10));
  fOutputContainer->Add(new TH2F("hFullDispvsCoreE","full dispersion vs full E;E (GeV);dispersion (#sigma)",100,0,50,100,0,10));
  fOutputContainer->Add(new TH2F("hCoreDispvsFullE","core dispersion vs core E;E (GeV);dispersion (#sigma)",100,0,50,100,0,10));

  for(Int_t imod=1;imod<Nmod;imod++) fOutputContainer->Add(new TH3F(Form("hdZvsZvsTrackPt_M%d",imod)        ,"dZ vs. Z;Z (cm);dZ (cm);p_{T}^{track} (GeV/c)"           ,160,-80,80,80,-20,20,40,0,20));
  for(Int_t imod=1;imod<Nmod;imod++) fOutputContainer->Add(new TH3F(Form("hdXvsXvsTrackPt_plus_M%d",imod)   ,"dX vs. X positive;X (cm);dX (cm);p_{T}^{track +} (GeV/c)",160,-80,80,80,-20,20,40,0,20));
  for(Int_t imod=1;imod<Nmod;imod++) fOutputContainer->Add(new TH3F(Form("hdXvsXvsTrackPt_minus_M%d",imod)  ,"dX vs. X negative;X (cm);dX (cm);p_{T}^{track -} (GeV/c)",160,-80,80,80,-20,20,40,0,20));
  for(Int_t imod=1;imod<Nmod;imod++) fOutputContainer->Add(new TH3F(Form("hdZvsZvsTrackPtElectron_M%d",imod),"dZ vs. Z of e^{#pm};Z (cm);dZ (cm);p_{T}^{track} (GeV/c)",160,-80,80,80,-20,20,40,0,20));//for radial displacement

  const Int_t NpTggModule = 71;
  Double_t pTggModule[NpTggModule]={};
  for(Int_t i=0;i<50;i++)     pTggModule[i] = 0.1 * i;            //every 0.1 GeV/c, up to 5 GeV/c
  for(Int_t i=50;i<60;i++)    pTggModule[i] = 0.5 * (i-50) + 5.0; //every 0.5 GeV/c, up to 10 GeV/c
  for(Int_t i=60;i<NpTggModule;i++) pTggModule[i] = 1.0 * (i-60) + 10.0;//every 1.0 GeV/c, up to 20 GeV/c

  for(Int_t imod=1;imod<Nmod;imod++){
    fOutputContainer->Add(new TH2F(Form("hEpRatiovsEnergy_M%d_Electron",imod) ,Form("E/p ratio vs. E_{cluster} M%d;E/p;E_{cluster} (GeV)",imod)      ,50,0,2,NpTggModule-1,pTggModule));
    fOutputContainer->Add(new TH2F(Form("hEpRatiovsEnergy_M%d_Others",imod)   ,Form("E/p ratio vs. E_{cluster} M%d;E/p;E_{cluster} (GeV)",imod)      ,50,0,2,NpTggModule-1,pTggModule));
    fOutputContainer->Add(new TH2F(Form("hEpRatiovsTrackPt_M%d_Electron",imod),Form("E/p ratio vs. p_{T}^{track} M%d;E/p;p_{T}^{track} (GeV/c)",imod),50,0,2,NpTggModule-1,pTggModule));
    fOutputContainer->Add(new TH2F(Form("hEpRatiovsTrackPt_M%d_Others",imod)  ,Form("E/p ratio vs. p_{T}^{track} M%d;E/p;p_{T}^{track} (GeV/c)",imod),50,0,2,NpTggModule-1,pTggModule));
  }

  fOutputContainer->Add(new TH2F("hEpRatiovsNsigmaElectronTPC","E/p ratio vs. N_{#sigma}^{e};E/p;n#sigma^{e}",50,0,2,20,-5,5));
  fOutputContainer->Add(new TH2F("hTPCdEdx_Electron","TPC dE/dx vs. electron momentum;p^{track} (GeV/c);dE/dx (a.u.)"    ,200,0,20,200,0,200));
  fOutputContainer->Add(new TH2F("hTPCdEdx_Others"  ,"TPC dE/dx vs. non-electron momentum;p^{track} (GeV/c);dE/dx (a.u.)",200,0,20,200,0,200));

  fOutputContainer->Add(new TH2F("hClusterEtaPhi","Cluster eta vs. phi;#phi;#eta",100,0,TwoPi,200,-1,1));
  fOutputContainer->Add(new TH2F("hEnergyvsDistanceToBadChannel","distance to closest bad channel;E (GeV);distance in cell",100,0,50,10,0,5));

  //<- histograms for QA
  //histograms for physics analysis ->

  const Int_t NpTgg = 101;
  Double_t pTgg[NpTgg]={};

  for(Int_t i=0;i<50;i++)     pTgg[i] = 0.1 * i;            //every 0.1 GeV/c, up to 5 GeV/c
  for(Int_t i=50;i<60;i++)    pTgg[i] = 0.5 * (i-50) + 5.0; //every 0.5 GeV/c, up to 10 GeV/c
  for(Int_t i=60;i<NpTgg;i++) pTgg[i] = 1.0 * (i-60) + 10.0;//every 1.0 GeV/c, up to 50 GeV/c

  Int_t NbinQ   =  1;
  Double_t Qmin = -1;
  Double_t Qmax = +1;
  TString axistitle = Form("cos(%d #Delta#phi)",fHarmonics);
  if(fFM == AliAnalysisTaskPHOSPi0EtaToGammaGamma::kEP){
    NbinQ = 24;
    Qmin  = -1;
    Qmax  = +1; 
    axistitle = Form("cos(%d #Delta#phi)",fHarmonics);
  }
  else if(fFM == AliAnalysisTaskPHOSPi0EtaToGammaGamma::kSP){
    NbinQ = 50;
    Qmin  = -0.5;
    Qmax  = +0.5;
    if(fQnStep == 5){//In Qn correction task, rescaling of Q vector is applied for V0 after step4
      //adjusted for centrality 10-60, where event plane is well defined.
      NbinQ = 40;
      Qmin  = -10;
      Qmax  = +10;
    }
    axistitle = "#vec{u} #upoint #vec{Q_{1}}";
  }

  const Int_t Ndimg = 2;
  const Int_t Nbing[Ndimg]    = {500, NbinQ};
  const Double_t xming[Ndimg] = {  0,  Qmin};
  const Double_t xmaxg[Ndimg] = { 50,  Qmax};
  //+/-0.3 is optimized for centrality 10-30, 30-50.

  //photon pT
  THnSparseF *hs_Photon = new THnSparseF("hSparsePhoton",Form("Photon;p_{T} (GeV/c);%s;",axistitle.Data()),Ndimg,Nbing,xming,xmaxg);
  hs_Photon->Sumw2();
  fOutputContainer->Add(hs_Photon);

  THnSparseF *hs_Photon_TOF = new THnSparseF("hSparsePhoton_TOF",Form("Photon with TOF;p_{T} (GeV/c);%s;",axistitle.Data()),Ndimg,Nbing,xming,xmaxg);
  hs_Photon_TOF->Sumw2();
  fOutputContainer->Add(hs_Photon_TOF);

  const Int_t Ndim = 4;
  const Int_t Nbin[Ndim]    = { 240, 500, 10, NbinQ};
  const Double_t xmin[Ndim] = {   0,   0,  0,  Qmin};
  const Double_t xmax[Ndim] = {0.96,  50,  1,  Qmax};

  //same event for pi0 region
  THnSparseF *hs_Mgg = new THnSparseF("hSparseMgg",Form("M_{#gamma#gamma};M_{#gamma#gamma} (GeV/c^{2});p_{T} (GeV/c);asymmetry;%s;",axistitle.Data()),Ndim,Nbin,xmin,xmax);
  hs_Mgg->Sumw2();
  fOutputContainer->Add(hs_Mgg);

  THnSparseF *hs_Mgg_TOF = new THnSparseF("hSparseMgg_TOF",Form("M_{#gamma#gamma} with TOF;M_{#gamma#gamma} (GeV/c^{2});p_{T} (GeV/c);asymmetry;%s;",axistitle.Data()),Ndim,Nbin,xmin,xmax);
  hs_Mgg_TOF->Sumw2();
  fOutputContainer->Add(hs_Mgg_TOF);

  //mixed event for pi0 region
  THnSparseF *hs_MixMgg = new THnSparseF("hSparseMixMgg",Form("M_{#gamma#gamma}^{mix};M_{#gamma#gamma} (GeV/c^{2});p_{T} (GeV/c);asymmetry;%s;",axistitle.Data()),Ndim,Nbin,xmin,xmax);
  hs_MixMgg->Sumw2();
  fOutputContainer->Add(hs_MixMgg);

  THnSparseF *hs_MixMgg_TOF = new THnSparseF("hSparseMixMgg_TOF",Form("M_{#gamma#gamma}^{mix} with TOF;M_{#gamma#gamma} (GeV/c^{2});p_{T} (GeV/c);asymmetry;%s;",axistitle.Data()),Ndim,Nbin,xmin,xmax);
  hs_MixMgg_TOF->Sumw2();
  fOutputContainer->Add(hs_MixMgg_TOF);

  if(fAnaOmega3Pi){

    const Int_t Ndim_omega = 4;
    const Int_t Nbin_omega[Ndim_omega]    = { 240, 500, 10, NbinQ};
    const Double_t xmin_omega[Ndim_omega] = { 0.4,   0,  0,  Qmin};
    const Double_t xmax_omega[Ndim_omega] = {1.36,  50,  1,  Qmax};

    //same event
    THnSparseF *hs_M3pi = new THnSparseF("hSparseM3pi",Form("M_{#pi^{0}#pi^{+}#pi^{-}};M_{#pi^{0}#pi^{+}#pi^{-}} (GeV/c^{2});p_{T} (GeV/c);asymmetry;%s;",axistitle.Data()),Ndim_omega,Nbin_omega,xmin_omega,xmax_omega);
    hs_M3pi->Sumw2();
    fOutputContainer->Add(hs_M3pi);

    THnSparseF *hs_M3pi_TOF = new THnSparseF("hSparseM3pi_TOF",Form("M_{#pi^{0}#pi^{+}#pi^{-}} with TOF;M_{#pi^{0}#pi^{+}#pi^{-}} (GeV/c^{2});p_{T} (GeV/c);asymmetry;%s;",axistitle.Data()),Ndim_omega,Nbin_omega,xmin_omega,xmax_omega);
    hs_M3pi_TOF->Sumw2();
    fOutputContainer->Add(hs_M3pi_TOF);

  }

  //<- histograms for physics analysis

  for(Int_t imod1=1;imod1<Nmod;imod1++){
    for(Int_t imod2=imod1;imod2<Nmod;imod2++){
      if(imod2 - imod1 > 1) continue;

      TH2F *h2 = new TH2F(Form("hMgg_M%d%d",imod1,imod2),"M_{#gamma#gamma} vs p_{T}",180,0,0.72,NpTggModule-1,pTggModule);
      h2->Sumw2();
      fOutputContainer->Add(h2);
    }
  }
  for(Int_t imod1=1;imod1<Nmod;imod1++){
    for(Int_t imod2=imod1;imod2<Nmod;imod2++){
      if(imod2 - imod1 > 1) continue;

      TH2F *h2 = new TH2F(Form("hMixMgg_M%d%d",imod1,imod2),"M_{#gamma#gamma} vs p_{T}",180,0,0.72,NpTggModule-1,pTggModule);
      h2->Sumw2();
      fOutputContainer->Add(h2);
    }
  }

  for(Int_t imod1=1;imod1<Nmod;imod1++){
    for(Int_t imod2=imod1;imod2<Nmod;imod2++){
      if(imod2 - imod1 > 1) continue;

      TH2F *h2 = new TH2F(Form("hMgg_M%d%d_TOF",imod1,imod2),"M_{#gamma#gamma} vs p_{T}",180,0,0.72,NpTggModule-1,pTggModule);
      h2->Sumw2();
      fOutputContainer->Add(h2);
    }
  }

  for(Int_t imod1=1;imod1<Nmod;imod1++){
    for(Int_t imod2=imod1;imod2<Nmod;imod2++){
      if(imod2 - imod1 > 1) continue;

      TH2F *h2 = new TH2F(Form("hMixMgg_M%d%d_TOF",imod1,imod2),"M_{#gamma#gamma} vs p_{T}",180,0,0.72,NpTggModule-1,pTggModule);
      h2->Sumw2();
      fOutputContainer->Add(h2);
    }
  }
  //for PID cut study

  if(fPIDStudy){
    const Int_t Ndim_PID = 5;
    const Int_t Nbin_PID[Ndim_PID]    = { 180, 500,   50, 50, 50};//Mgg vs. pT vs. Ncell vs. M02 vs. M20
    const Double_t xmin_PID[Ndim_PID] = {   0,   0,  0.5,  0,  0};//Mgg vs. pT vs. Ncell vs. M02 vs. M20
    const Double_t xmax_PID[Ndim_PID] = {0.72,  50, 50.5,  5,  5};//Mgg vs. pT vs. Ncell vs. M02 vs. M20

    //same event for PID study at low pT
    THnSparseF *hs_Mgg_PID = new THnSparseF("hSparseMgg_PID","M_{#gamma#gamma} for PID;M_{#gamma#gamma} (GeV/c^{2});E_{#gamma} (GeV);N_{cell};M20 (cm);M02 (cm);",Ndim_PID,Nbin_PID,xmin_PID,xmax_PID);
    hs_Mgg_PID->Sumw2();
    fOutputContainer->Add(hs_Mgg_PID);

    //mixed event for PID study at low pT
    THnSparseF *hs_MixMgg_PID = new THnSparseF("hSparseMixMgg_PID","M_{#gamma#gamma}^{mix} for PID;M_{#gamma#gamma} (GeV/c^{2});E_{#gamma} (GeV);N_{cell};M020 (cm);M02 (cm);",Ndim_PID,Nbin_PID,xmin_PID,xmax_PID);
    hs_MixMgg_PID->Sumw2();
    fOutputContainer->Add(hs_MixMgg_PID);
  }

  //for PID cut efficiency
  const TString PIDtype[4] = {"noPID","CPV","Disp","PID"};//str PID is for main PID cut efficiency used in this analysis

  TH2F *h2_p_PID     = new TH2F(Form("hMgg_Probe_%s"          ,PIDtype[3].Data()),Form("Probe #gamma %s;M_{#gamma#gamma} (GeV/c^{2});p_{T}^{#gamma} (GeV/c)"            ,PIDtype[3].Data()),180,0,0.72,NpTgg-1,pTgg);
  h2_p_PID->Sumw2();
  fOutputContainer->Add(h2_p_PID);
  TH2F *h2mix_p_PID  = new TH2F(Form("hMixMgg_Probe_%s"       ,PIDtype[3].Data()),Form("Mix Probe #gamma %s;M_{#gamma#gamma} (GeV/c^{2});p_{T}^{#gamma} (GeV/c)"        ,PIDtype[3].Data()),180,0,0.72,NpTgg-1,pTgg);
  h2mix_p_PID->Sumw2();
  fOutputContainer->Add(h2mix_p_PID);
  for(Int_t ip=1;ip<4;ip++){

    TH2F *h2_pp_PID    = new TH2F(Form("hMgg_PassingProbe_%s"   ,PIDtype[ip].Data()),Form("Passing Probe #gamma %s;M_{#gamma#gamma} (GeV/c^{2});p_{T}^{#gamma} (GeV/c)"    ,PIDtype[ip].Data()),180,0,0.72,NpTgg-1,pTgg);
    TH2F *h2mix_pp_PID = new TH2F(Form("hMixMgg_PassingProbe_%s",PIDtype[ip].Data()),Form("Mix Passing Probe #gamma %s;M_{#gamma#gamma} (GeV/c^{2});p_{T}^{#gamma} (GeV/c)",PIDtype[ip].Data()),180,0,0.72,NpTgg-1,pTgg);

    h2_pp_PID->Sumw2();
    fOutputContainer->Add(h2_pp_PID);

    h2mix_pp_PID->Sumw2();
    fOutputContainer->Add(h2mix_pp_PID);

  }//end of PID loop

  //for TOF cut efficiency
  TH2F *h2_p_TOF = new TH2F("hMgg_Probe_TOF"              ,"Probe #gamma TOF;M_{#gamma#gamma} (GeV/c^{2});E_{#gamma} (GeV)"            ,60,0,0.24,NpTgg-1,pTgg);
  h2_p_TOF->Sumw2();                                                                                                                    
  fOutputContainer->Add(h2_p_TOF);                                                                                                      
                                                                                                                                        
  TH2F *h2_pp_TOF = new TH2F("hMgg_PassingProbe_TOF"      ,"Passing Probe #gamma TOF;M_{#gamma#gamma} (GeV/c^{2});E_{#gamma} (GeV)"    ,60,0,0.24,NpTgg-1,pTgg);
  h2_pp_TOF->Sumw2();                                                                                                                   
  fOutputContainer->Add(h2_pp_TOF);                                                                                                     
                                                                                                                                        
  TH2F *h2mix_p_TOF = new TH2F("hMixMgg_Probe_TOF"        ,"Mix Probe #gamma TOF;M_{#gamma#gamma} (GeV/c^{2});E_{#gamma} (GeV)"        ,60,0,0.24,NpTgg-1,pTgg);
  h2mix_p_TOF->Sumw2();                                                                                                                 
  fOutputContainer->Add(h2mix_p_TOF);                                                                                                   
                                                                                                                                        
  TH2F *h2mix_pp_TOF = new TH2F("hMixMgg_PassingProbe_TOF","Mix Passing Probe #gamma TOF;M_{#gamma#gamma} (GeV/c^{2});E_{#gamma} (GeV)",60,0,0.24,NpTgg-1,pTgg);
  h2mix_pp_TOF->Sumw2();
  fOutputContainer->Add(h2mix_pp_TOF);

  //for photon purity by DDA

  for(Int_t ip=0;ip<4;ip++){
    TH1F *h1PhotonPt = new TH1F(Form("hPhotonPt_%s",PIDtype[ip].Data()),Form("#gamma p_{T} %s;p_{T}^{#gamma} (GeV/c)",PIDtype[ip].Data()),NpTgg-1,pTgg);
    h1PhotonPt->Sumw2();
    fOutputContainer->Add(h1PhotonPt);
  }//end of PID loop

  const TString PIDName[] = {"Electron","Pion","Kaon","Proton","AntiProton","K0L","Neutron","AntiNeutron","Gamma","Others"};//AntiProton is needed for antineutron study.
  const Int_t Npid = sizeof(PIDName)/sizeof(PIDName[0]);

  for(Int_t ip=0;ip<5;ip++){
    fOutputContainer->Add(new TH2F(Form("hRvsTrackPt_%s"     ,PIDName[ip].Data()),Form("r vs track p_{T} %s;p_{T}^{track} (GeV/c);cpv (#sigma)"    ,PIDName[ip].Data()),500,0,50,200,0,20));
    fOutputContainer->Add(new TH2F(Form("hRvsClusterPt_%s"   ,PIDName[ip].Data()),Form("r vs cluster p_{T} %s;p_{T}^{cluster} (GeV/c);cpv (#sigma)",PIDName[ip].Data()),500,0,50,200,0,20));
    fOutputContainer->Add(new TH2F(Form("hMixRvsTrackPt_%s"  ,PIDName[ip].Data()),Form("r vs track p_{T} %s;p_{T}^{track} (GeV/c);cpv (#sigma)"    ,PIDName[ip].Data()),500,0,50,200,0,20));
    fOutputContainer->Add(new TH2F(Form("hMixRvsClusterPt_%s",PIDName[ip].Data()),Form("r vs cluster p_{T} %s;p_{T}^{cluster} (GeV/c);cpv (#sigma)",PIDName[ip].Data()),500,0,50,200,0,20));
  }

  for(Int_t ip=0;ip<5;ip++){
    TH1F *h1noPID = new TH1F(Form("hMatched%s",PIDName[ip].Data()),Form("p_{T} of %s in clusters for purity no PID;p_{T} (GeV/c)",PIDName[ip].Data()),NpTgg-1,pTgg);
    h1noPID->Sumw2();
    fOutputContainer->Add(h1noPID);
    
    TH1F *h1PID = new TH1F(Form("hMatched%s_Disp",PIDName[ip].Data()),Form("p_{T} of %s in clusters for purity Disp;p_{T} (GeV/c)",PIDName[ip].Data()),NpTgg-1,pTgg);
    h1PID->Sumw2();
    fOutputContainer->Add(h1PID);
  };

  const TString str[4] = {"L0","L1L","L1M","L1H"};
  TH1F *h1all = new TH1F("hAllClusterEnergy","all cluster energy;E_{cluster} (GeV)",500,0,50);//for MB
  h1all->Sumw2();
  fOutputContainer->Add(h1all);

  TH1F *h1all_GT = new TH1F("hAllClusterEnergyGoodTRU","allcluster energy on good TRU;E_{cluster} (GeV)",500,0,50);//for MB
  h1all_GT->Sumw2();
  fOutputContainer->Add(h1all_GT);

  for(Int_t it=0;it<4;it++){
    TH1F *h1matched = new TH1F(Form("hMatchedClusterEnergy%s",str[it].Data()),"matched cluster energy;E_{cluster} (GeV)",500,0,50);//for MB
    h1matched->Sumw2();
    fOutputContainer->Add(h1matched);
  }//end of trigger type loop

  for(Int_t imod=1;imod<Nmod;imod++){
    for(Int_t itru=1;itru<=Ntru;itru++){
      TH1F *h1all_mt     = new TH1F(Form("hAllClusterEnergyM%dTRU%d",imod,itru)    ,Form("all cluster energy M%d TRU%d;E_{cluster} (GeV)",imod,itru)    ,500,0,50);//for MB
      h1all_mt->Sumw2();
      fOutputContainer->Add(h1all_mt);
    }//end of TRU loop
  }//end of module loop

  for(Int_t it=0;it<4;it++){
    for(Int_t imod=1;imod<Nmod;imod++){
      for(Int_t itru=1;itru<=Ntru;itru++){

        TH1F *h1matched_mt = new TH1F(Form("hMatchedClusterEnergy%sM%dTRU%d",str[it].Data(),imod,itru),Form("matched cluster energy %s M%d TRU%d;E_{cluster} (GeV)",str[it].Data(),imod,itru),500,0,50);//for MB
        h1matched_mt->Sumw2();
        fOutputContainer->Add(h1matched_mt);

      }//end of TRU loop
    }//end of module loop
  }//end of trigger type loop

  if(fIsPHOSTriggerAnalysis){//for your trigger configuration

    //for Trigger QA
    fOutputContainer->Add(new TH1F("hNFiredTRUChannel","Fired TRU channel per event",101,-0.5,100.5));
    for(Int_t imod=1;imod<Nmod;imod++) fOutputContainer->Add(new TH2F(Form("hFiredTRUChannelM%d",imod),Form("Fired TRU channel M%d",imod),64,0.5,64.5, 56,0.5,56.5));
    for(Int_t imod=1;imod<Nmod;imod++) fOutputContainer->Add(new TH3F(Form("hMatchedFiredTRUChannelM%d",imod),Form("Fired TRU channel M%d",imod),64,0.5,64.5, 56,0.5,56.5,100,0,50));

    for(Int_t imod=1;imod<Nmod;imod++){
      fOutputContainer->Add(new TH3F(Form("hClusterHitM%d",imod)       ,Form("cluster hit M%d",imod)         ,64,0.5,64.5,56,0.5,56.5,100,0,50));
      fOutputContainer->Add(new TH3F(Form("hMatchedClusterHitM%d",imod),Form("matched cluster hit M%d",imod) ,64,0.5,64.5,56,0.5,56.5,100,0,50));

      for(Int_t itru=1;itru<=Ntru;itru++){
        fOutputContainer->Add(new TH1F(Form("hClusterEnergyM%dTRU%d",imod,itru)       ,Form("cluster energy M%d",imod),500,0,50));
        fOutputContainer->Add(new TH1F(Form("hMatchedClusterEnergyM%dTRU%d",imod,itru),Form("cluster energy M%d",imod),500,0,50));
        fOutputContainer->Add(new TH3F(Form("hDistanceXZM%dTRU%d",imod,itru)        ,Form("distance TRUch-clu M%d TRU%d",imod,itru),21,-10.5,+10.5,21,-10.5,+10.5,100,0,50));
        fOutputContainer->Add(new TH3F(Form("hMatchedDistanceXZM%dTRU%d",imod,itru) ,Form("distance TRUch-clu M%d TRU%d",imod,itru),21,-10.5,+10.5,21,-10.5,+10.5,100,0,50));
        fOutputContainer->Add(new TH2F(Form("hDeltaRM%dTRU%d",imod,itru)        ,Form("#DeltaR TRUch-clu M%d TRU%d",imod,itru),100,0,0.1,100,0,50));
        fOutputContainer->Add(new TH2F(Form("hMatchedDeltaRM%dTRU%d",imod,itru) ,Form("#DeltaR TRUch-clu M%d TRU%d",imod,itru),100,0,0.1,100,0,50));
      }
    }

    for(Int_t imod=1;imod<Nmod;imod++){
      for(Int_t itru=1;itru<=Ntru;itru++){

        TH1F *h1all     = new TH1F(Form("hClusterE_M%d_TRU%d"         ,imod,itru)     ,Form("cluster energy M%d TRU%d"          ,imod,itru),NpTgg-1,pTgg);
        TH1F *h1trg     = new TH1F(Form("hTriggeredClusterE_M%d_TRU%d",imod,itru)     ,Form("triggered cluster energy M%d TRU%d",imod,itru),NpTgg-1,pTgg);
        TH1F *h1all_tof = new TH1F(Form("hClusterE_TOF_M%d_TRU%d"         ,imod,itru) ,Form("cluster energy M%d TRU%d"          ,imod,itru),NpTgg-1,pTgg);
        TH1F *h1trg_tof = new TH1F(Form("hTriggeredClusterE_TOF_M%d_TRU%d",imod,itru) ,Form("triggered cluster energy M%d TRU%d",imod,itru),NpTgg-1,pTgg);

        h1all->Sumw2();
        h1trg->Sumw2();
        h1all_tof->Sumw2();
        h1trg_tof->Sumw2();

        fOutputContainer->Add(h1all);
        fOutputContainer->Add(h1trg);
        fOutputContainer->Add(h1all_tof);
        fOutputContainer->Add(h1trg_tof);

      }
    }

    //for trigger efficiency

    TH2F *h2_p_Trg = new TH2F("hMgg_Probe_Trg"              ,"Probe #gamma Trg;M_{#gamma#gamma} (GeV/c^{2});E_{#gamma} (GeV)"            ,180,0,0.72,NpTgg-1,pTgg);
    h2_p_Trg->Sumw2();                                                                                                                    
    fOutputContainer->Add(h2_p_Trg);                                                                                                      

    TH2F *h2_pp_Trg = new TH2F("hMgg_PassingProbe_Trg"      ,"Passing Probe #gamma Trg;M_{#gamma#gamma} (GeV/c^{2});E_{#gamma} (GeV)"    ,180,0,0.72,NpTgg-1,pTgg);
    h2_pp_Trg->Sumw2();                                                                                                                   
    fOutputContainer->Add(h2_pp_Trg);                                                                                                     

    TH2F *h2mix_p_Trg = new TH2F("hMixMgg_Probe_Trg"       ,"Mix Probe #gamma Trg;M_{#gamma#gamma} (GeV/c^{2});E_{#gamma} (GeV)"         ,180,0,0.72,NpTgg-1,pTgg);
    h2mix_p_Trg->Sumw2();                                                                                                                 
    fOutputContainer->Add(h2mix_p_Trg);                                                                                                   

    TH2F *h2mix_pp_Trg = new TH2F("hMixMgg_PassingProbe_Trg","Mix Passing Probe #gamma Trg;M_{#gamma#gamma} (GeV/c^{2});E_{#gamma} (GeV)",180,0,0.72,NpTgg-1,pTgg);
    h2mix_pp_Trg->Sumw2();
    fOutputContainer->Add(h2mix_pp_Trg);

    for(Int_t imod=1;imod<Nmod;imod++){
      for(Int_t itru=1;itru<=Ntru;itru++){
        TH2F *h2 = new TH2F(Form("hMgg_Probe_Trg_M%d_TRU%d"          ,imod,itru),"Probe #gamma Trg;M_{#gamma#gamma} (GeV/c^{2});E_{#gamma} (GeV)"            ,180,0,0.72,NpTgg-1,pTgg);
        h2->Sumw2();
        fOutputContainer->Add(h2);
      }
    }

    for(Int_t imod=1;imod<Nmod;imod++){
      for(Int_t itru=1;itru<=Ntru;itru++){
        TH2F *h2 = new TH2F(Form("hMgg_PassingProbe_Trg_M%d_TRU%d"   ,imod,itru),"Passing Probe #gamma Trg;M_{#gamma#gamma} (GeV/c^{2});E_{#gamma} (GeV)"    ,180,0,0.72,NpTgg-1,pTgg);
        h2->Sumw2();
        fOutputContainer->Add(h2);
      }
    }

    for(Int_t imod=1;imod<Nmod;imod++){
      for(Int_t itru=1;itru<=Ntru;itru++){
        TH2F *h2 = new TH2F(Form("hMixMgg_Probe_Trg_M%d_TRU%d"       ,imod,itru),"Mix Probe #gamma Trg;M_{#gamma#gamma} (GeV/c^{2});E_{#gamma} (GeV)"        ,180,0,0.72,NpTgg-1,pTgg);
        h2->Sumw2();
        fOutputContainer->Add(h2);
      }
    }

    for(Int_t imod=1;imod<Nmod;imod++){
      for(Int_t itru=1;itru<=Ntru;itru++){
        TH2F *h2 = new TH2F(Form("hMixMgg_PassingProbe_Trg_M%d_TRU%d",imod,itru),"Mix Passing Probe #gamma Trg;M_{#gamma#gamma} (GeV/c^{2});E_{#gamma} (GeV)",180,0,0.72,NpTgg-1,pTgg);
        h2->Sumw2();
        fOutputContainer->Add(h2);
      }
    }

  }//fIsPHOSTriggerAnalysis ends

  if(fIsOAStudy){
    //opening angle cut study
    TH3F *h3MggvsPt_OA = new TH3F("hMgg_OA","M_{#gamma#gamma} for opening angle;M_{#gamma#gamma} (GeV/c^{2});p_{T} (GeV/c);opening angle (mrad.)" ,180,0,0.72,100,0,100,20,0,20);
    h3MggvsPt_OA->Sumw2();
    fOutputContainer->Add(h3MggvsPt_OA);

    TH3F *h3MixMggvsPt_OA = new TH3F("hMixMgg_OA","M_{#gamma#gamma} for opening angle;M_{#gamma#gamma} (GeV/c^{2});p_{T} (GeV/c);opening angle (mrad.)" ,180,0,0.72,100,0,100,20,0,20);
    h3MixMggvsPt_OA->Sumw2();
    fOutputContainer->Add(h3MixMggvsPt_OA);
  }

  if(fIsMC){

    //fOutputContainer->Add(new TH1F("hPDGPhysicalPrimary","PDG code of Phyical Primary",8001,-4000-0.5,4000+0.5));
    //fOutputContainer->Add(new TH1F("hPDGPhysicalPrimaryStable","PDG code of Phyical Primary",8001,-4000-0.5,4000+0.5));

    const TString parname[] = {"Pi0","Eta","Gamma","Omega","ChargedPion","ChargedKaon","K0S","K0L","Lambda0","Sigma0","Proton","AntiProton","Neutron","AntiNeutron"};
    const Int_t Npar = sizeof(parname)/sizeof(parname[0]);

    for(Int_t ipar=0;ipar<Npar;ipar++){
      TH1F *h1Pt = new TH1F(Form("hGen%sPt",parname[ipar].Data()),Form("generated %s pT;p_{T} (GeV/c)",parname[ipar].Data()),NpTgg-1,pTgg);
      h1Pt->Sumw2();
      fOutputContainer->Add(h1Pt);

      TH2F *h2EtaPhi = new TH2F(Form("hGen%sEtaPhi",parname[ipar].Data()),Form("generated %s y vs phi;#phi (rad);rapidity",parname[ipar].Data()),60,0,TwoPi,200,-1,1);
      h2EtaPhi->Sumw2();
      fOutputContainer->Add(h2EtaPhi);

      TH2F *h2EtaPt = new TH2F(Form("hGen%sEtaPt",parname[ipar].Data()),Form("generated %s y vs pT;rapidity;p_{T} (GeV/c)",parname[ipar].Data()),200,-1,1,NpTgg-1,pTgg);
      h2EtaPt->Sumw2();
      fOutputContainer->Add(h2EtaPt);

      if(parname[ipar]=="Pi0" || parname[ipar]=="Eta"){
        //check detector acceptance
        //these histograms are filled only when 2 gammas are in PHOS acceptance.

        TH1F *h1PtACC = new TH1F(Form("hGen%sPtACC",parname[ipar].Data()),Form("generated %s pT;p_{T} (GeV/c)",parname[ipar].Data()),NpTgg-1,pTgg);
        h1PtACC->Sumw2();
        fOutputContainer->Add(h1PtACC);

        TH2F *h2EtaPhiACC = new TH2F(Form("hGen%sEtaPhiACC",parname[ipar].Data()),Form("generated %s y vs phi;#phi (rad);rapidity",parname[ipar].Data()),60,0,TwoPi,200,-1,1);
        h2EtaPhiACC->Sumw2();
        fOutputContainer->Add(h2EtaPhiACC);

        TH2F *h2EtaPtACC = new TH2F(Form("hGen%sEtaPtACC",parname[ipar].Data()),Form("generated %s y vs pT;rapidity;p_{T} (GeV/c)",parname[ipar].Data()),200,-1,1,NpTgg-1,pTgg);
        h2EtaPtACC->Sumw2();
        fOutputContainer->Add(h2EtaPtACC);
      }

    }//end of particle loop


    //const TString PIDtype[4] = {"noPID","CPV","Disp","PID"};//str PID is for main PID cut efficiency used in this analysis
    //for purity based on MC truth
    for(Int_t ip=0;ip<Npid;ip++){
      for(Int_t jp=0;jp<4;jp++){

        TH1F *h1purity = new TH1F(Form("hPurity%s_%s",PIDName[ip].Data(),PIDtype[jp].Data()),Form("p_{T} of true %s in clusters for purity %s;p_{T} (GeV/c)",PIDName[ip].Data(),PIDtype[jp].Data()),NpTgg-1,pTgg);
        h1purity->Sumw2();
        fOutputContainer->Add(h1purity);

      }
    }

    //converted electrons.
    for(Int_t jp=0;jp<4;jp++){
      TH1F *h1purity_lce = new TH1F(Form("hPurityLCE_%s",PIDtype[jp].Data()),Form("p_{T} of true late conversion Electron in clusters for purity %s;p_{T} (GeV/c)",PIDtype[jp].Data()),NpTgg-1,pTgg);
      h1purity_lce->Sumw2();
      fOutputContainer->Add(h1purity_lce);
    }

    for(Int_t jp=0;jp<4;jp++){
      TH2F *h2e = new TH2F(Form("hElectronRxy_%s",PIDtype[jp].Data()),Form("p_{T} of true Electron in clusters for purity %s vs. V_{prod};p_{T} (GeV/c);R_{xy} (cm)",PIDtype[jp].Data()),NpTgg-1,pTgg,500,0,500);
      h2e->Sumw2();
      fOutputContainer->Add(h2e);

      TH2F *h2ce = new TH2F(Form("hConvertedElectronRxy_%s",PIDtype[jp].Data()),Form("p_{T} of true converted Electron in clusters for purity %s vs. V_{prod};p_{T} (GeV/c);R_{xy} (cm)",PIDtype[jp].Data()),NpTgg-1,pTgg,500,0,500);
      h2ce->Sumw2();
      fOutputContainer->Add(h2ce);
    }
   

    //for feed down correction
    TH1F *h1gamma_K0S = new TH1F("hGammaFromK0S","#gamma from K^{0}_{S};p_{T} (GeV/c)",NpTgg-1,pTgg);
    h1gamma_K0S->Sumw2();
    fOutputContainer->Add(h1gamma_K0S);
    TH1F *h1gamma_L0 = new TH1F("hGammaFromL0","#gamma from #Lambda;p_{T} (GeV/c)",NpTgg-1,pTgg);
    h1gamma_L0->Sumw2();
    fOutputContainer->Add(h1gamma_L0);
    TH1F *h1gamma_S0 = new TH1F("hGammaFromS0","#gamma from #Sigma^{0};p_{T} (GeV/c)",NpTgg-1,pTgg);
    h1gamma_S0->Sumw2();
    fOutputContainer->Add(h1gamma_S0);

    const TString Asym[] = {"","_asym08"};
    const Int_t Nasym = sizeof(Asym)/sizeof(Asym[0]);

    for(Int_t iasym=0;iasym<Nasym;iasym++){
      TH2F *h2 = new TH2F(Form("hMggFromPi0%s",Asym[iasym].Data()),"M_{#gamma#gamma} from #pi^{0};M_{#gamma#gamma} (GeV/c^{2});p_{T} (GeV/c)",60,0,0.24,NpTgg-1,pTgg);
      h2->Sumw2();
      fOutputContainer->Add(h2);
    }//end of asym loop

    for(Int_t iasym=0;iasym<Nasym;iasym++){
      TH2F *h2 = new TH2F(Form("hMggFromK0S%s",Asym[iasym].Data()),"M_{#gamma#gamma} from K_{S}^{0};M_{#gamma#gamma} (GeV/c^{2});p_{T} (GeV/c)",60,0,0.24,NpTgg-1,pTgg);
      h2->Sumw2();
      fOutputContainer->Add(h2);
    }//end of asym loop

    for(Int_t iasym=0;iasym<Nasym;iasym++){
      TH2F *h2 = new TH2F(Form("hMggFromLambda0%s",Asym[iasym].Data()),"M_{#gamma#gamma} from #Lambda;M_{#gamma#gamma} (GeV/c^{2});p_{T} (GeV/c);",60,0,0.24,NpTgg-1,pTgg);
      h2->Sumw2();
      fOutputContainer->Add(h2);

    }//end of asym loop

    if(fIsNonLinStudy){
      //for systematic non-linearity study
      const Int_t Na = 7;
      const Int_t Nb = 7;

      for(Int_t ia=0;ia<Na;ia++){
        for(Int_t ib=0;ib<Nb;ib++){
          Double_t a = -0.09 + 0.01*ia;
          Double_t b =  0.4  + 0.1 *ib;

          fNonLin[ia][ib] = new TF1(Form("fNonLin_a%d_b%d",ia,ib),"[2]*(1. + [0] / (1. + TMath::Power(x/[1],2)))",0,100);
          fNonLin[ia][ib]->SetParNames("a","b (GeV)");
          fNonLin[ia][ib]->FixParameter(0,a);
          fNonLin[ia][ib]->FixParameter(1,b);
          fNonLin[ia][ib]->FixParameter(2,fGlobalEScale);
 
          TH2F *h2same_NL = new TH2F(Form("hMgg_a%d_b%d",ia,ib)   ,Form("M_{#gamma#gamma} vs. p_{T} a = %3.2f , b = %3.2f;M_{#gamma#gamma} (GeV/c^{2});p_{T} (GeV/c)",a,b)      ,240,0,0.96,NpTgg-1,pTgg);
          TH2F *h2mix_NL  = new TH2F(Form("hMixMgg_a%d_b%d",ia,ib),Form("M_{#gamma#gamma}^{mix} vs. p_{T} a = %3.2f , b = %3.2f;M_{#gamma#gamma} (GeV/c^{2});p_{T} (GeV/c)",a,b),240,0,0.96,NpTgg-1,pTgg);
          h2same_NL->Sumw2();
          h2mix_NL ->Sumw2();
 
          fOutputContainer->Add(h2same_NL);
          fOutputContainer->Add(h2mix_NL);
        }//end of ib
      }//end of ia
    }

    //for JJMC
    if(fIsJJMC){
      fOutputContainer->Add(new TH1F("hPtHard","pT hard in GeV/c;p_{T} hard (GeV/c)",1000,0,1000));
      fOutputContainer->Add(new TH1F("hNTrial","nTrial;p_{T} hard bin",20,0.5,20.5));
      fOutputContainer->Add(new TProfile("hProfCrossSection","inelastic cross section;p_{T} hard bin;<#sigma^{INEL}>",20,0.5,20.5));

      TH1F *hNMerged = new TH1F("hNMerged","N merged",20,0.5,20.5);
      hNMerged->SetYTitle("number of merged files");
      fOutputContainer->Add(hNMerged);
    }//end of IsJJMC

  }//end of IsMC

  if(fIsFlowTask){
    AliAnalysisTaskFlowVectorCorrections *flowQnVectorTask = dynamic_cast<AliAnalysisTaskFlowVectorCorrections*>(AliAnalysisManager::GetAnalysisManager()->GetTask("FlowQnVectorCorrections"));
    if(flowQnVectorTask != NULL){
      fFlowQnVectorMgr = flowQnVectorTask->GetAliQnCorrectionsManager();
    }
    else {
      AliFatal("Flow Qn vector corrections framework needed but it is not present. ABORTING!!!");
    }
  }

  PostData(1,fOutputContainer);

}
//________________________________________________________________________
void AliAnalysisTaskPHOSPi0EtaToGammaGamma::UserExec(Option_t *option) 
{
  // Main loop
  // Called for each event

  if(!fPHOSEventCuts){
    AliError("fPHOSEventCuts is not set! return");
    return;
  }
  if(!fPHOSClusterCuts){
    AliError("fPHOSClusterCuts is not set! return");
    return;
  }

  fEvent = dynamic_cast<AliVEvent*>(InputEvent());
  if(!fEvent){
    AliError("event is not available.");
    return;
  }

  fESDEvent = dynamic_cast<AliESDEvent*>(fEvent);
  fAODEvent = dynamic_cast<AliAODEvent*>(fEvent);

  //select event which PHOS was readout from trigger cluster point of view.
  //for example, PHOS is not in MUFAST cluster.
  TString trigClasses = fEvent->GetFiredTriggerClasses();
  //printf("%s\n",trigClasses.Data());
  if(!fIsMC 
      && !trigClasses.Contains("-CENT") //accept CENT, CENTNO[TRD|PMD]
      && !trigClasses.Contains("-FAST") //accept FAST, not MUFAST
      && !trigClasses.Contains("-CALO") //accept CALO, CALOFAST
    ){
    AliWarning(Form("Skip event with triggers %s",trigClasses.Data()));
    return;
  }

  fPIDResponse = fInputHandler->GetPIDResponse();
  if(!fPIDResponse){
    AliFatal("fPIDResponse does not exist!");
    return;
  }

  if(fRunNumber != fEvent->GetRunNumber()){ // Check run number
    fRunNumber = fEvent->GetRunNumber();
    fPHOSGeo = GetPHOSGeometry();
  }

  UInt_t fSelectMask = fInputHandler->IsEventSelected();
  Bool_t isINT7selected = fSelectMask & AliVEvent::kINT7;

  //reject pile up by physics selection
  if(!fIsPHOSTriggerAnalysis && !isINT7selected){
    AliInfo("INT7 Event is rejected by physics selection.");
    return;
  }

  //Bool_t isPHI7selected = fSelectMask & AliVEvent::kPHI7;
  Bool_t isPHI7selected = fSelectMask & (AliVEvent::kPHI7|AliVEvent::kCaloOnly);
  if(!fIsMC && fIsPHOSTriggerAnalysis && !isPHI7selected){
    AliInfo("PHI7 Event is rejected by physics selection.");
    return;
  }

  Int_t L0input  = 17;//for LHC17
  Int_t L1Hinput = 7;
  Int_t L1Minput = 6;
  Int_t L1Linput = 5;
  if(fRunNumber <= 246994)  L0input = 9;//for LHC15

  Bool_t Is0PH0fired = fEvent->GetHeader()->GetL0TriggerInputs() & 1 << (L0input  - 1);//trigger input -1
  Bool_t Is1PHHfired = fEvent->GetHeader()->GetL1TriggerInputs() & 1 << (L1Hinput - 1);//trigger input -1
  Bool_t Is1PHMfired = fEvent->GetHeader()->GetL1TriggerInputs() & 1 << (L1Minput - 1);//trigger input -1
  Bool_t Is1PHLfired = fEvent->GetHeader()->GetL1TriggerInputs() & 1 << (L1Linput - 1);//trigger input -1

  //As of 20180617, PHI7 in [CALO/CALOFAST] is moved AliVEvent::kCaloOnly (or AliVEvent::kMuonCalo).
  //EMC triggers and PHOS triggers are merged to 1 bit.
  //First, select PHOS triggered event by bit operation
  if(!fIsMC && fIsPHOSTriggerAnalysis && 
      (!Is0PH0fired && !Is1PHHfired && !Is1PHMfired && !Is1PHLfired)
    ){
    AliInfo("PHI7 Event is rejected by bit operation.");
    return;
  }

  const AliVVertex *vVertex = fEvent->GetPrimaryVertex();
  fVertex[0] = vVertex->GetX();
  fVertex[1] = vVertex->GetY();
  fVertex[2] = vVertex->GetZ();
  fZvtx = (Int_t)((fVertex[2]+10.)/2.);//it should be 0-9.
  if(fZvtx < 0) fZvtx = 0;//protection to avoid fZvtx = -1.
  if(fZvtx > 9) fZvtx = 9;//protection to avoid fZvtx = 10.

  //for JJMC
  if(fIsMC && fIsJJMC){

    //get pT hard bin from file path

    AliAODInputHandler *aodH = dynamic_cast<AliAODInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
    TTree* tr = aodH->GetTree();
    TDirectory* dr = tr->GetDirectory();
    TString DataPath = dr->GetPath();

    TObjArray *tx = DataPath.Tokenize("/");
    //tx->Print();
    //for (Int_t i = 0; i < tx->GetEntries(); i++) std::cout << "i = " << i << " , " << ((TObjString *)(tx->At(i)))->String() << std::endl;

    TString binstr = ((TObjString *)(tx->At(5)))->String();
    Int_t pThardbin = binstr.Atoi();
    fPtHardBin = pThardbin;

    fJJMCHandler = new AliPHOSJetJetMC(pThardbin);
    fJJMCHandler->ConfigureJetJetMC(fEvent);//this should be called in each event to get pythia event header

    if(fPtHardAndJetPtFactor > 0)    fJJMCHandler->SetJetPtFactor(fPtHardAndJetPtFactor);
    if(fPtHardAndSinglePtFactor > 0) fJJMCHandler->SetSingleParticlePtFactor(fPtHardAndSinglePtFactor);

    AliMCEvent *mcevent = MCEvent();
    AliGenPythiaEventHeader* pythiaGenHeader = 0x0;
    if(fESDEvent) pythiaGenHeader      = fJJMCHandler->GetPythiaEventHeader(mcevent);
    else if(fAODEvent) pythiaGenHeader = fJJMCHandler->GetPythiaEventHeader(fAODEvent);

    Float_t pThard   = pythiaGenHeader->GetPtHard();
    Float_t xsection = pythiaGenHeader->GetXsection();
    Int_t nTrial     = pythiaGenHeader->Trials();

    AliInfo(Form("pT-hard-bin = %d , pThard = %f GeV/c , xsection = %f mb , nTrial = %d.",fPtHardBin, pThard, xsection, nTrial));
    FillHistogramTH1(fOutputContainer,"hNTrial",fPtHardBin,nTrial);//this should be filled in all events (no rejected/selected events).
    FillProfile(fOutputContainer,"hProfCrossSection",fPtHardBin,xsection);//this should be filled in all events (no rejected/selected events).

    if(!(fJJMCHandler->ComparePtHardBin(fEvent))) return;
    if(!(fJJMCHandler->ComparePtHardWithJet(fEvent))) return;
    if(!(fJJMCHandler->ComparePtHardWithSingleParticle(fEvent))) return;
  }

  Int_t Ncontributor  = vVertex->GetNContributors();

  //centrality estimation

  Float_t fCentralityV0M = -1.;
  Float_t fCentralityCL0 = -1.;
  Float_t fCentralityCL1 = -1.;
  Float_t fCentralityV0A = -1.;
  Float_t fCentralityV0C = -1.;
  Float_t fCentralityZNA = -1.;
  Float_t fCentralityZNC = -1.;

  if(fEstimator.Contains("V0") || fEstimator.Contains("ZN") || fEstimator.Contains("CL")){

    //Get Centrality
    fMultSelection = (AliMultSelection*)fEvent->FindListObject("MultSelection");
    if(!fMultSelection){
      //If you get this warning (and fCentralityV0M 300) please check that the AliMultSelectionTask actually ran (before your task)
      AliWarning("AliMultSelection object not found!");
      return;
    }
    else{
      fCentralityV0M  = fMultSelection->GetMultiplicityPercentile("V0M");
      fCentralityCL0  = fMultSelection->GetMultiplicityPercentile("CL0");
      fCentralityCL1  = fMultSelection->GetMultiplicityPercentile("CL1");
      fCentralityV0A  = fMultSelection->GetMultiplicityPercentile("V0A");
      fCentralityV0C  = fMultSelection->GetMultiplicityPercentile("V0C");
      fCentralityZNA  = fMultSelection->GetMultiplicityPercentile("ZNA");
      fCentralityZNC  = fMultSelection->GetMultiplicityPercentile("ZNC");
      fCentralityMain = fMultSelection->GetMultiplicityPercentile(fEstimator);
    }

  }
  else if(fEstimator.Contains("HybridTrack")){
    //hybrid track multiplicity 
    //done manually in this task.
    Int_t NHybrid = 0; 
    const Int_t trackMult = fEvent->GetNumberOfTracks();
    if(fESDEvent){
      for(Int_t itrack=0;itrack<trackMult;itrack++){
        AliESDtrack *esdtrack = (AliESDtrack*)fEvent->GetTrack(itrack);
        if(TMath::Abs(esdtrack->Eta()) > 0.8) continue;

        if(fESDtrackCutsGlobal           ->AcceptTrack(esdtrack)//select global track
        || fESDtrackCutsGlobalConstrained->AcceptTrack(esdtrack)){//select complementary track
          NHybrid++;
        }

      }//end of track loop
    }//end of ESD
    else if(fAODEvent){
      for(Int_t itrack=0;itrack<trackMult;itrack++){
        AliAODTrack *aodtrack = (AliAODTrack*)fEvent->GetTrack(itrack);
        if(TMath::Abs(aodtrack->Eta()) > 0.8) continue;

        if(aodtrack->IsHybridGlobalConstrainedGlobal()){//hybrid track
          NHybrid++;
        }

      }//end of track loop
    }//end of AOD
    fCentralityMain = (Float_t)NHybrid;
  }
  else if(fEstimator.Contains("SPDTracklet")){
    //hybrid track multiplicity 
    fCentralityMain = (Float_t)(fEvent->GetMultiplicity()->GetNumberOfTracklets());
  }
  else{
    AliInfo(Form("%s is not supported. return",fEstimator.Data()));
    return;
  }

  if(fCentralityMain < fCentralityMin || fCentralityMax < fCentralityMain){
    AliInfo(Form("Reject this event because centrality %s %f %% is out of the configuration of this task.", fEstimator.Data(),fCentralityMain));
    return;
  }

  FillHistogramTH1(fOutputContainer,"hEventSummary",1);//all
  FillHistogramTH1(fOutputContainer,"hVertexZ" ,fVertex[2]);
  FillHistogramTH2(fOutputContainer,Form("hCentrality%svsNContributor",fEstimator.Data()),fCentralityMain,Ncontributor);

  fPHOSClusterArray = (TClonesArray*)fEvent->FindListObject("PHOSClusterArray");
  if(!fPHOSClusterArray){
    AliWarning("fPHOSClusterArray object not found!");
    return;
  }

  if(!fIsMC && !fIsPHOSTriggerAnalysis && (244917 <= fRunNumber && fRunNumber <= 246994)){
    //only for PbPb MB analysis
    const Int_t multClust = fPHOSClusterArray->GetEntriesFast();
    Int_t NclusterTOF = 0;
    for(Int_t i1=0;i1<multClust;i1++){
      AliCaloPhoton *ph = (AliCaloPhoton*)fPHOSClusterArray->At(i1);
      if(!CheckMinimumEnergy(ph)) continue;

      if(ph->IsTOFOK()) NclusterTOF++;

    }//end of cluster loop

    AliVVZERO *V0info = (AliVVZERO*)fEvent->GetVZEROData();
    Float_t mtotV0A = V0info->GetMTotV0A();
    Float_t mtotV0C = V0info->GetMTotV0C();
    Float_t mtotV0 = mtotV0A + mtotV0C;
    if(NclusterTOF < 1.2e-3 * mtotV0 - 12) return;

  }

  if(fIsFlowTask){
    Bool_t QnOK = ExtractQnVector();
    if(!QnOK){
      AliInfo("Event is rejected by Qn vector quality.");
      return;
    }
  }
  else fEPBin = 0;

  //event selection
  if(!(fPHOSEventCuts->AcceptEvent(fEvent))){
    AliInfo("event is rejected.");
    return;
  }

  //fill fired trigger statistics right after basic event selection, but before PHI7 event selection.
  if(Is0PH0fired) FillHistogramTH1(fOutputContainer,"hEventSummary",3);//0PH0
  if(Is1PHLfired) FillHistogramTH1(fOutputContainer,"hEventSummary",4);//1PHL
  if(Is1PHMfired) FillHistogramTH1(fOutputContainer,"hEventSummary",5);//1PHM
  if(Is1PHHfired) FillHistogramTH1(fOutputContainer,"hEventSummary",6);//1PHH

  if(fIsPHOSTriggerAnalysis && !(fPHOSTriggerHelper->IsPHI7(fEvent,fPHOSClusterCuts,fEmin,fEnergyThreshold,fUseCoreEnergy))){
    AliInfo("event is rejected. IsPHI7 = kFALSE.");
    return;
  }

  fPHOSTriggerHelper->IsPHI7(fEvent,fPHOSClusterCuts,fEmin,fEnergyThreshold,fUseCoreEnergy);//only to set event in fPHOSTriggerHelper. //do nothing. just return kTRUE/kFALSE

  //<- end of event selection
  //-> start physics analysis

  FillHistogramTH2(fOutputContainer,"hCentralityV0MvsCL0",fCentralityV0M,fCentralityCL0);
  FillHistogramTH2(fOutputContainer,"hCentralityV0MvsCL1",fCentralityV0M,fCentralityCL1);
  FillHistogramTH2(fOutputContainer,"hCentralityCL0vsCL1",fCentralityCL0,fCentralityCL1);
  FillHistogramTH2(fOutputContainer,"hCentralityV0AvsV0C",fCentralityV0A,fCentralityV0C);
  FillHistogramTH2(fOutputContainer,"hCentralityZNAvsZNC",fCentralityZNA,fCentralityZNC);

  FillHistogramTH1(fOutputContainer,"hVertexZSelectEvent" ,fVertex[2]);
  FillHistogramTH1(fOutputContainer,"hEventSummary",2);//selected event
 
  AliInfo(Form("Collision system = %d | fCentralityMain estimated by %s = %f %% | Zvtx = %f cm , fZvtx = %d | Harmonics = %d , fEventPlane = %f (rad.) , fEPBin = %d |",fCollisionSystem,fEstimator.Data(),fCentralityMain,fVertex[2],fZvtx,fHarmonics,fEventPlane,fEPBin));
 
  if(!fIsMC) FillRejectionFactorMB();

  //track QA
  TrackQA();

  if(fIsMC) ProcessMC();

  if(fIsMC){//fill cluster occupancy in M.C. to obtain total(i.e. HIJNG/HIJING+JJ) Ncluster.
    Int_t multPHOSClustAll = 0;
    for(Int_t i1=0;i1<fPHOSClusterArray->GetEntriesFast();i1++){
      AliCaloPhoton *ph = (AliCaloPhoton*)fPHOSClusterArray->At(i1);
      if(!fPHOSClusterCuts->AcceptPhoton(ph)) continue;
      if(!CheckMinimumEnergy(ph)) continue;

      multPHOSClustAll++;
    }
    FillHistogramTH2(fOutputContainer,Form("hCentrality%svsPHOSClusterMultiplicityMC",fEstimator.Data()),fCentralityMain,multPHOSClustAll);
  }

  if(fIsMC && fIsJJMC){
    //This should be called after ProcessMC(), because fMCArrayAOD is used.
    //it is possible to analyze HIJING event for efficiency calculation even in HIJING + JJ.

    Int_t firstJetindex = fJJMCHandler->GetFirstJetIndex();
    Int_t lastJetindex  = fJJMCHandler->GetLastJetIndex();
    Int_t genIDJet      = fJJMCHandler->GetGeneratorJetIndex();
    Int_t firstUEindex  = fJJMCHandler->GetFirstUEIndex();
    Int_t lastUEindex   = fJJMCHandler->GetLastUEIndex();
    Int_t genIDUE       = fJJMCHandler->GetGeneratorUEIndex();

    if(fESDEvent){
      const Int_t multClust = fPHOSClusterArray->GetEntriesFast();
      for(Int_t i=0;i<multClust;i++){
        AliCaloPhoton *ph = (AliCaloPhoton*)fPHOSClusterArray->At(i);
        Int_t primary = FindPrimaryMotherESD(ph->GetPrimary());

        if(fMCType.Contains("JJMC") && (primary < firstJetindex || lastJetindex < primary)) fPHOSClusterArray->Remove(ph);
        if(fMCType.Contains("MBMC") && (primary < firstUEindex  || lastUEindex  < primary)) fPHOSClusterArray->Remove(ph);

      }
      fPHOSClusterArray->Compress();

    }
    else if(fAODEvent){
      genIDJet = fJJMCHandler->GetGeneratorJetIndex();
      genIDUE  = fJJMCHandler->GetGeneratorUEIndex();

      const Int_t multClust = fPHOSClusterArray->GetEntriesFast();
      for(Int_t i=0;i<multClust;i++){
        AliCaloPhoton *ph = (AliCaloPhoton*)fPHOSClusterArray->At(i);
        Int_t primary = ph->GetPrimary();
        AliAODMCParticle *p = (AliAODMCParticle*)fMCArrayAOD->At(primary);
        Int_t genID = p->GetGeneratorIndex();

        if(fMCType.Contains("JJMC") && genID != genIDJet) fPHOSClusterArray->Remove(ph);
        if(fMCType.Contains("MBMC") && genID != genIDUE)  fPHOSClusterArray->Remove(ph);
      }
      fPHOSClusterArray->Compress();
    }
  }
  if(!fPHOSEvents[fZvtx][fEPBin]) fPHOSEvents[fZvtx][fEPBin] = new TList();
  TList *prevPHOS = fPHOSEvents[fZvtx][fEPBin];

  if(!fIsMC && fIsPHOSTriggerAnalysis){
    AliInfo(Form("PHOS trigger analysis is ON! RF method = %d",fTRFM));
    TriggerQA();
    SelectTriggeredCluster();
    EstimateTriggerEfficiency();
  }

  //cell QA
  CellQA();
  //cluster QA
  ClusterQA();

  if(fIsMC) SetMCWeight();

  FillPhoton();
  FillEpRatio();
  FillMgg();
  FillMixMgg();
  if(fAnaOmega3Pi) FillM3pi();

  EstimatePIDCutEfficiency();
  EstimateTOFCutEfficiency();
  DDAPhotonPurity();
  FillTrackMatching();
  FillMixTrackMatching();

  if(fIsNonLinStudy) DoNonLinearityStudy();

  //Now we either add current events to stack or remove
  //If no photons in current event - no need to add it to mixed
  if(fPHOSClusterArray->GetEntriesFast() > 0){
    //don't call fPHOSClucster=0; this will affect original array provided from PHOSbjectCreator.
    //prevPHOS->AddFirst(fPHOSClusterArray);
    //fPHOSClusterArray=0;

    TClonesArray *clone = new TClonesArray(*fPHOSClusterArray);
    prevPHOS->AddFirst(clone);
    //delete clone;
    clone = 0;

    if(prevPHOS->GetSize() > fNMixed){//Remove redundant events
      TClonesArray * tmp = static_cast<TClonesArray*>(prevPHOS->Last());
      prevPHOS->RemoveLast();
      delete tmp;
      tmp = NULL;
    }
  }

  if(fJJMCHandler){
    delete fJJMCHandler;
    fJJMCHandler = 0x0;
  }

  PostData(1, fOutputContainer);
}
//________________________________________________________________________
void AliAnalysisTaskPHOSPi0EtaToGammaGamma::Terminate(Option_t *option) 
{
  //Called once at the end of the query
  //In principle, this function is not needed...

  AliInfo(Form("%s is done.",GetName()));

}
//________________________________________________________________________
void AliAnalysisTaskPHOSPi0EtaToGammaGamma::TrackQA() 
{
  const Int_t trackMult = fEvent->GetNumberOfTracks();
  const Int_t Ntracklet = fEvent->GetMultiplicity()->GetNumberOfTracklets();

  Double_t pT=0, eta=0, phi=0, dEdx=0, p=0;
  Int_t NGlobal=0;
  Int_t NGlobalNoDCA=0;
  Int_t NHybrid=0;
  Int_t NComplementary=0;
  Float_t dca_xy = 999, dca_z = 999;

  if(fESDEvent){
    for(Int_t itrack=0;itrack<trackMult;itrack++){
      AliESDtrack *esdtrack = (AliESDtrack*)fEvent->GetTrack(itrack);
      if(TMath::Abs(esdtrack->Eta()) > 0.8) continue;

      if(!fESDtrackCutsGlobal           ->AcceptTrack(esdtrack)//select global track
      && !fESDtrackCutsGlobalConstrained->AcceptTrack(esdtrack)){//select complementary track
        //cout << "itrack = " << itrack << " , rejected" << endl;
        continue;
      }
      NHybrid++;
      //cout << "itrack = " << itrack << " , accepted" << endl;
      pT  = esdtrack->Pt();
      eta = esdtrack->Eta();
      phi = esdtrack->Phi();
      if(phi<0) phi += TMath::TwoPi();
      p = esdtrack->P();
      dEdx = esdtrack->GetTPCsignal();
      FillHistogramTH2(fOutputContainer,"hTrackTPCdEdx",p,dEdx);

      FillHistogramTH1(fOutputContainer,"hHybridTrackPt",pT);
      FillHistogramTH2(fOutputContainer,"hHybridTrackEtaPhi",phi,eta);

      if(fESDtrackCutsGlobal->AcceptTrack(esdtrack)){//global track
        NGlobal++;
        FillHistogramTH1(fOutputContainer,"hGlobalTrackPt",pT);
        FillHistogramTH2(fOutputContainer,"hGlobalTrackEtaPhi",phi,eta);
      }
      if(fESDtrackCutsGlobalConstrained->AcceptTrack(esdtrack)){//global-constrained = complementary track
        NComplementary++;
        FillHistogramTH1(fOutputContainer,"hComplementaryTrackPt",pT);
        FillHistogramTH2(fOutputContainer,"hComplementaryTrackEtaPhi",phi,eta);
      }
    }//end of track loop
  }//end of ESD
  else if(fAODEvent){
    for(Int_t itrack=0;itrack<trackMult;itrack++){
      AliAODTrack *aodtrack = (AliAODTrack*)fEvent->GetTrack(itrack);
      if(TMath::Abs(aodtrack->Eta()) > 0.8) continue;
      if(aodtrack->Pt() < 0.15) continue;

      pT = aodtrack->Pt();
      eta = aodtrack->Eta();
      phi = aodtrack->Phi();
      if(phi<0) phi += TMath::TwoPi();
      p = aodtrack->P();
      dEdx = aodtrack->GetTPCsignal();

      dca_xy = 999; dca_z = 999;
      aodtrack->GetImpactParameters(dca_xy,dca_z);

      if(aodtrack->TestFilterMask(AliAODTrack::kTrkGlobalNoDCA)){//standard cuts with very loose DCA cut
        NGlobalNoDCA++;
        FillHistogramTH1(fOutputContainer,"hGlobalNoDCATrackPt",pT);
        FillHistogramTH2(fOutputContainer,"hGlobalNoDCATrackEtaPhi",phi,eta);
        FillHistogramTH2(fOutputContainer,"hGlobalNoDCATrackDCA",dca_xy,dca_z);
        FillHistogramTH2(fOutputContainer,"hGlobalNoDCATrackTPCdEdx",p,dEdx);
      }//end of loose DCA
      if(aodtrack->TestFilterMask(AliAODTrack::kTrkGlobal)){//standard cuts with tight DCA cut (AliESDtrackCuts::GetStandardITSTPCTrackCuts2011())
        NGlobal++;
        FillHistogramTH1(fOutputContainer,"hGlobalTrackPt",pT);
        FillHistogramTH2(fOutputContainer,"hGlobalTrackEtaPhi",phi,eta);
        FillHistogramTH2(fOutputContainer,"hGlobalTrackDCA",dca_xy,dca_z);
        FillHistogramTH2(fOutputContainer,"hGlobalTrackTPCdEdx",p,dEdx);
      }//end of tight DCA

    }//end of track loop

  }//end of AOD

  FillHistogramTH1(fOutputContainer,"hGlobalTrackMult",NGlobal);
  FillHistogramTH1(fOutputContainer,"hGlobalNoDCATrackMult",NGlobalNoDCA);
  FillHistogramTH2(fOutputContainer,Form("hCentrality%svsTrackMultiplicity",fEstimator.Data()),fCentralityMain,NGlobalNoDCA);
  FillHistogramTH2(fOutputContainer,Form("hCentrality%svsSPDTracklet",fEstimator.Data()),fCentralityMain,Ntracklet);
  fNHybridTrack = NGlobalNoDCA;

}
//________________________________________________________________________
void AliAnalysisTaskPHOSPi0EtaToGammaGamma::CellQA() 
{
  AliVCaloCells *cells = dynamic_cast<AliVCaloCells*>(fEvent->GetPHOSCells());
  const Int_t multCells = cells->GetNumberOfCells();

  Int_t relId[4]={};
  Int_t nCellModule[5] = {};
  Double_t cellamp=0;
  Double_t celltime=0;
  Int_t module=0,cellx=0,cellz=0;
  Int_t cellAbsId=0;
  Bool_t isHG=kFALSE;

  for(Int_t iCell=0; iCell<multCells; iCell++){
    cellAbsId = cells->GetCellNumber(iCell);

    if(cellAbsId<0) continue;//reject CPV
    isHG = cells->GetCellHighGain(cellAbsId);
    cellamp = cells->GetCellAmplitude(cellAbsId);
    celltime = cells->GetCellTime(cellAbsId);

    fPHOSGeo->AbsToRelNumbering(cellAbsId,relId);

    module = relId[0];
    cellx  = relId[2];
    cellz  = relId[3];

    if(module < 1 || module > 4){
      AliError(Form("Wrong module number %d",module));
      return;
    }

    if(isHG) FillHistogramTH2(fOutputContainer,Form("hCellAmpTimeM%d_HG",module),cellamp,celltime*1e+9);
    else     FillHistogramTH2(fOutputContainer,Form("hCellAmpTimeM%d_LG",module),cellamp,celltime*1e+9);

    nCellModule[module]++;
    FillHistogramTH2(fOutputContainer,Form("hCellNXZM%d",module),cellx,cellz);
    FillHistogramTH2(fOutputContainer,Form("hCellEXZM%d",module),cellx,cellz,cellamp);
  }

  FillHistogramTH1(fOutputContainer,"hCellMultEventM1",nCellModule[1]);
  FillHistogramTH1(fOutputContainer,"hCellMultEventM2",nCellModule[2]);
  FillHistogramTH1(fOutputContainer,"hCellMultEventM3",nCellModule[3]);
  FillHistogramTH1(fOutputContainer,"hCellMultEventM4",nCellModule[4]);

}
//________________________________________________________________________
void AliAnalysisTaskPHOSPi0EtaToGammaGamma::TriggerQA() 
{
  const Int_t multClust = fPHOSClusterArray->GetEntriesFast();
  //AliVCaloCells *cells = dynamic_cast<AliVCaloCells*>(fEvent->GetPHOSCells());
  AliVCaloTrigger* trg = fEvent->GetCaloTrigger("PHOS");
  trg->Reset();

  const Int_t Nfired = trg->GetEntries();
  FillHistogramTH1(fOutputContainer,"hNFiredTRUChannel",Nfired);

  Int_t trgrelId[4]={};
  Int_t tmod=0; //"Online" module number
  Int_t trgabsId=0; //bottom-left 4x4 edge cell absId [1-64], [1-56]
  Int_t trgmodule=0,trgcellx=0,trgcellz=0;
  Float_t x=0, y=0; // cell coordinates in local reference system of a module

  Int_t relId[4]={};
  Int_t module=0,cellx=0,cellz=0;
  Double_t energy=0;
  Float_t position[3] = {};
  Double_t eta_cluster = 0;
  Double_t phi_cluster = 0;
  Double_t eta_trigger = 0;
  Double_t phi_trigger = 0;
  Double_t DeltaR = -1;

  Bool_t IsFilledOnce = kFALSE;//flag to avoid double counting

  const Int_t L1input = fPHOSTriggerHelper->GetL1TriggerInput();
  const Int_t L0input = fPHOSTriggerHelper->GetL0TriggerInput();

  Int_t L1=-999;
  if(L1input > 0){
    if     (L1input == 7) L1 = 0;//L1 high
    else if(L1input == 6) L1 = 1;//L1 medium
    else if(L1input == 5) L1 = 2;//L1 low
  }
  else if(L0input > 0){
    //if(L0input == 9) L1 = -1;//L0
    L1 = -1;//L0
  }

  Int_t kUsedCluster[] = {multClust*0};

  while(trg->Next()){

    FillHistogramTH1(fOutputContainer,"hPHI7Summary",trg->GetL1TimeSum());

    if(trg->GetL1TimeSum() != L1) continue;
    trg->GetPosition(tmod,trgabsId);//tmod is online module numbering. i.e., tmod=1 means half module.

    fPHOSGeo->AbsToRelNumbering(trgabsId,trgrelId);
    //for offline numbering, relId should be used.

    trgmodule = trgrelId[0];
    trgcellx  = trgrelId[2];
    trgcellz  = trgrelId[3];

    fPHOSGeo->RelPosInModule(trgrelId,x,y); 
    TVector3 trgglobal;
    fPHOSGeo->Local2Global(trgrelId[0],x,y,trgglobal);
    eta_trigger = trgglobal.Eta();
    phi_trigger = trgglobal.Phi();
    if(phi_trigger < 0) phi_trigger += TMath::TwoPi();

    //cout << "trgmodule = "<< trgmodule << endl;

    if(trgmodule < 1 || trgmodule > 4){
      AliError(Form("Wrong module number %d",trgmodule));
      return;
    }

    FillHistogramTH1(fOutputContainer,Form("hFiredTRUChannelM%d",trgmodule),trgcellx,trgcellz);

    for(Int_t i=0;i<multClust;i++){
      AliCaloPhoton *ph = (AliCaloPhoton*)fPHOSClusterArray->At(i);
      if(!fPHOSClusterCuts->AcceptPhoton(ph)) continue;
      //AliVCluster *clu1 = (AliVCluster*)ph->GetCluster();

      energy = ph->Energy();

      position[0] = ph->EMCx();
      position[1] = ph->EMCy();
      position[2] = ph->EMCz();

      relId[0] = 0; relId[1] = 0; relId[2] = 0; relId[3] = 0;
      TVector3 global1(position);
      eta_cluster = global1.Eta();
      phi_cluster = global1.Phi();
      if(phi_cluster < 0) phi_cluster += TMath::TwoPi();

      fPHOSGeo->GlobalPos2RelId(global1,relId);

      //Int_t maxAbsId = fPHOSTriggerHelper->FindHighestAmplitudeCellAbsId(clu1, cells);
      //fPHOSGeo->AbsToRelNumbering(maxAbsId,relId);

      module = relId[0];
      cellx  = relId[2];
      cellz  = relId[3];

      if(module < 1 || module > 4){
        AliError(Form("Wrong module number %d",module));
        return;
      }

      Int_t diffx = trgcellx - cellx;
      Int_t diffz = trgcellz - cellz;
      DeltaR = TMath::Sqrt( TMath::Power(eta_trigger - eta_cluster,2) +  TMath::Power(phi_trigger - phi_cluster,2) );
      //DeltaR = fPHOSTriggerHelper->GetDistanceToClosestTRUChannel(ph);

      FillHistogramTH3(fOutputContainer,Form("hDistanceXZM%dTRU%d",trgmodule,fPHOSTriggerHelper->WhichTRU(trgcellx,trgcellz)),diffx,diffz,energy);
      FillHistogramTH2(fOutputContainer,Form("hDeltaRM%dTRU%d",trgmodule,fPHOSTriggerHelper->WhichTRU(trgcellx,trgcellz)),DeltaR,energy);

      if(!IsFilledOnce){//to avoid double counting of energy
        FillHistogramTH1(fOutputContainer,Form("hClusterEnergyM%dTRU%d",module,fPHOSTriggerHelper->WhichTRU(cellx,cellz)),energy);
        FillHistogramTH3(fOutputContainer,Form("hClusterHitM%d",module),cellx,cellz,energy);
      }

      if(fPHOSTriggerHelper->IsDeltaRUsed()){
        if(fPHOSTriggerHelper->IsMatchedDeltaR(trgrelId,global1)){
          FillHistogramTH3(fOutputContainer,Form("hMatchedFiredTRUChannelM%d",trgmodule),trgcellx,trgcellz,energy);
          FillHistogramTH3(fOutputContainer,Form("hMatchedDistanceXZM%dTRU%d",trgmodule,fPHOSTriggerHelper->WhichTRU(trgcellx,trgcellz)),diffx,diffz,energy);
          FillHistogramTH2(fOutputContainer,Form("hMatchedDeltaRM%dTRU%d",trgmodule,fPHOSTriggerHelper->WhichTRU(trgcellx,trgcellz)),DeltaR,energy);

          if(kUsedCluster[i] == 0){
            FillHistogramTH1(fOutputContainer,Form("hMatchedClusterEnergyM%dTRU%d",trgmodule,fPHOSTriggerHelper->WhichTRU(trgcellx,trgcellz)),energy);
            FillHistogramTH3(fOutputContainer,Form("hMatchedClusterHitM%d",module),cellx,cellz,energy);
          }
          kUsedCluster[i] = 1;
        }
      }
      else{
        if(fPHOSTriggerHelper->IsMatched(trgrelId,relId)){
          FillHistogramTH3(fOutputContainer,Form("hMatchedFiredTRUChannelM%d",trgmodule),trgcellx,trgcellz,energy);
          FillHistogramTH3(fOutputContainer,Form("hMatchedDistanceXZM%dTRU%d",trgmodule,fPHOSTriggerHelper->WhichTRU(trgcellx,trgcellz)),diffx,diffz,energy);
          FillHistogramTH2(fOutputContainer,Form("hMatchedDeltaRM%dTRU%d",trgmodule,fPHOSTriggerHelper->WhichTRU(trgcellx,trgcellz)),DeltaR,energy);

          if(kUsedCluster[i] == 0){
            FillHistogramTH1(fOutputContainer,Form("hMatchedClusterEnergyM%dTRU%d",trgmodule,fPHOSTriggerHelper->WhichTRU(trgcellx,trgcellz)),energy);
            FillHistogramTH3(fOutputContainer,Form("hMatchedClusterHitM%d",module),cellx,cellz,energy);
          }
          kUsedCluster[i] = 1;
        }
      }

    }//end of cluster loop

    IsFilledOnce = kTRUE;

  }//end of while trigger loop

}
//________________________________________________________________________
void AliAnalysisTaskPHOSPi0EtaToGammaGamma::ClusterQA() 
{
  const Int_t multClust = fPHOSClusterArray->GetEntriesFast();

  Int_t multPHOSClust[5]={};
  Int_t multPHOSClustTOF[5]={};

  Int_t module=0,cellx=0,cellz=0;
  Int_t relId[4]={};
  Double_t position[3] = {};
  Int_t digMult=0;
  Double_t energy=0,tof=0,eta=0,phi=0;
  Double_t M02=0, M20=0;
  Double_t R = 0, coreR=0;
  Double_t coreE = 0;

  for(Int_t i1=0;i1<multClust;i1++){
    AliCaloPhoton *ph = (AliCaloPhoton*)fPHOSClusterArray->At(i1);
    if(!fPHOSClusterCuts->AcceptPhoton(ph)) continue;
    if(!CheckMinimumEnergy(ph)) continue;

    if(fIsPHOSTriggerAnalysis){
      if( fIsMC && !fPHOSTriggerHelper->IsOnActiveTRUChannel(ph)) continue;//only for MC
      if(!fIsMC && !ph->IsTrig()) continue;//it is meaningless to focus on photon without fired trigger in PHOS triggered data.
    }
    if(fForceActiveTRU && !fPHOSTriggerHelper->IsOnActiveTRUChannel(ph)) continue;

    energy  = ph->Energy();
    if(fUseCoreEnergy){
      energy = (ph->GetMomV2())->Energy();
    }

    digMult = ph->GetNCells();
    tof     = ph->GetTime();//unit is second.
    M20 = ph->GetLambda1();
    M02 = ph->GetLambda2();

    position[0] = ph->EMCx();
    position[1] = ph->EMCy();
    position[2] = ph->EMCz();

    relId[0] = 0; relId[1] = 0; relId[2] = 0; relId[3] = 0;
    TVector3 global1(position);
    fPHOSGeo->GlobalPos2RelId(global1,relId);

    module = relId[0];
    cellx  = relId[2];
    cellz  = relId[3];

    if(module < 1 || module > 4){
      AliError(Form("Wrong module number %d",module));
      return;
    }

    if(ph->IsTOFOK()){
      Int_t tru = fPHOSTriggerHelper->WhichTRU(cellx,cellz);
      FillHistogramTH1(fOutputContainer,"hAllClusterEnergy",energy);//all cluster
      if(fPHOSTriggerHelper->IsOnActiveTRUChannel(ph)){
        FillHistogramTH1(fOutputContainer,"hAllClusterEnergyGoodTRU",energy);//all cluster on good TRU acceptance
        FillHistogramTH1(fOutputContainer,Form("hAllClusterEnergyM%dTRU%d",module,tru),energy);//all cluster
      }
    }
    FillHistogramTH2(fOutputContainer,"hEnergyvsDistanceToBadChannel",energy,ph->DistToBadfp()); //in unit of cell with floating point

    eta = global1.Eta();
    phi = global1.Phi();
    if(phi<0) phi += TMath::TwoPi();
    FillHistogramTH2(fOutputContainer,"hClusterEtaPhi",phi,eta);
    FillHistogramTH2(fOutputContainer,Form("hClusterEvsTM%d",module),energy,tof*1e+9);
    FillHistogramTH2(fOutputContainer,Form("hCluNXZM%d",module),cellx,cellz);
    FillHistogramTH2(fOutputContainer,Form("hCluEXZM%d",module),cellx,cellz,energy);

    R     = ph->GetNsigmaFullDisp();
    coreR = ph->GetNsigmaCoreDisp();
    coreE = (ph->GetMomV2())->Energy();
 
    FillHistogramTH2(fOutputContainer,"hClusterEvsN",energy,digMult);
    FillHistogramTH2(fOutputContainer,"hClusterEvsM02",energy,M02);
    FillHistogramTH2(fOutputContainer,"hClusterEvsM20",energy,M20);
    FillHistogramTH2(fOutputContainer,"hFullDispvsFullE",energy,R);
    FillHistogramTH2(fOutputContainer,"hCoreDispvsCoreE",coreE,coreR);
    FillHistogramTH2(fOutputContainer,"hFullDispvsCoreE",coreE,R);
    FillHistogramTH2(fOutputContainer,"hCoreDispvsFullE",energy,coreR);

    multPHOSClust[0]++;
    multPHOSClust[module]++;

    FillHistogramTH1(fOutputContainer,Form("hPHOSClusterEnergyM%d",module),energy);
    if(ph->IsTOFOK()){// unit in ns.
      multPHOSClustTOF[0]++;
      multPHOSClustTOF[module]++;
      FillHistogramTH2(fOutputContainer,Form("hCluNXZTOFM%d",module),cellx,cellz);
      FillHistogramTH2(fOutputContainer,Form("hCluEXZTOFM%d",module),cellx,cellz,energy);
      FillHistogramTH1(fOutputContainer,Form("hPHOSClusterEnergyTOFM%d",module),energy);
    }

    if(energy > 0.5){
      if(energy < 2.){// 0.5 < energy < 2GeV
        FillHistogramTH2(fOutputContainer,Form("hCluLowNXZM%d",module),cellx,cellz);
        FillHistogramTH2(fOutputContainer,Form("hCluLowEXZM%d",module),cellx,cellz,energy);
      }
      else{//energy>=2GeV
        FillHistogramTH2(fOutputContainer,Form("hCluHighNXZM%d",module),cellx,cellz);
        FillHistogramTH2(fOutputContainer,Form("hCluHighEXZM%d",module),cellx,cellz,energy);
      }
    }

  }//end of cluster loop

  FillHistogramTH1(fOutputContainer,"hPHOSClusterMultM1",multPHOSClust[1]);
  FillHistogramTH1(fOutputContainer,"hPHOSClusterMultM2",multPHOSClust[2]);
  FillHistogramTH1(fOutputContainer,"hPHOSClusterMultM3",multPHOSClust[3]);
  FillHistogramTH1(fOutputContainer,"hPHOSClusterMultM4",multPHOSClust[4]);

  //for TOF cut
  FillHistogramTH1(fOutputContainer,"hPHOSClusterMultTOFM1",multPHOSClustTOF[1]);
  FillHistogramTH1(fOutputContainer,"hPHOSClusterMultTOFM2",multPHOSClustTOF[2]);
  FillHistogramTH1(fOutputContainer,"hPHOSClusterMultTOFM3",multPHOSClustTOF[3]);
  FillHistogramTH1(fOutputContainer,"hPHOSClusterMultTOFM4",multPHOSClustTOF[4]);
  //for tof cut end

  FillHistogramTH2(fOutputContainer,Form("hCentrality%svsPHOSClusterMultiplicity"   ,fEstimator.Data()),fCentralityMain,multPHOSClust[0]);
  FillHistogramTH2(fOutputContainer,Form("hCentrality%svsPHOSClusterMultiplicityTOF",fEstimator.Data()),fCentralityMain,multPHOSClustTOF[0]);

  FillHistogramTH2(fOutputContainer,"hPHOSClusterMultiplicityvsTrackMultiplicity",fNHybridTrack,multPHOSClust[0]);
  FillHistogramTH2(fOutputContainer,"hPHOSClusterMultiplicityTOFvsTrackMultiplicity",fNHybridTrack,multPHOSClustTOF[0]);

  AliVVZERO *V0info = (AliVVZERO*)fEvent->GetVZEROData();
  Float_t mtotV0A = V0info->GetMTotV0A();
  Float_t mtotV0C = V0info->GetMTotV0C();
  Float_t mtotV0 = mtotV0A + mtotV0C;

  FillHistogramTH2(fOutputContainer,"hMultiplicityV0AvsPHOSClusterMultiplicity",mtotV0A,multPHOSClust[0]);
  FillHistogramTH2(fOutputContainer,"hMultiplicityV0CvsPHOSClusterMultiplicity",mtotV0C,multPHOSClust[0]);
  FillHistogramTH2(fOutputContainer,"hMultiplicityV0vsPHOSClusterMultiplicity" ,mtotV0 ,multPHOSClust[0]);
  FillHistogramTH2(fOutputContainer,"hMultiplicityV0AvsPHOSClusterMultiplicityTOF",mtotV0A,multPHOSClustTOF[0]);
  FillHistogramTH2(fOutputContainer,"hMultiplicityV0CvsPHOSClusterMultiplicityTOF",mtotV0C,multPHOSClustTOF[0]);
  FillHistogramTH2(fOutputContainer,"hMultiplicityV0vsPHOSClusterMultiplicityTOF" ,mtotV0 ,multPHOSClustTOF[0]);

  FillHistogramTH2(fOutputContainer,"hMultiplicityV0AvsTrackMultiplicity",mtotV0A,fNHybridTrack);
  FillHistogramTH2(fOutputContainer,"hMultiplicityV0CvsTrackMultiplicity",mtotV0C,fNHybridTrack);
  FillHistogramTH2(fOutputContainer,"hMultiplicityV0vsTrackMultiplicity" ,mtotV0 ,fNHybridTrack);

}
//________________________________________________________________________
void AliAnalysisTaskPHOSPi0EtaToGammaGamma::FillPhoton() 
{
  const Int_t multClust = fPHOSClusterArray->GetEntriesFast();

  Double_t pT=0,energy=0;
  Double_t phi = -999, dphi = -999.;
  Double_t eff=1;
  Double_t trgeff=1;
  TF1 *f1tof = GetTOFCutEfficiencyFunction();
  TF1 * f1trg = GetTriggerEfficiencyFunction();
  Double_t value[2] = {};
  Double_t sp1 = -999;
  Int_t primary = -1;
  Double_t weight = 1.;
  Double_t TrueK0SPt = 0;
  Double_t TrueL0Pt = 0;
  Double_t TrueS0Pt = 0;

  for(Int_t iph=0;iph<multClust;iph++){
    AliCaloPhoton *ph = (AliCaloPhoton*)fPHOSClusterArray->At(iph);
    if(!CheckMinimumEnergy(ph)) continue;

    if(fIsPHOSTriggerAnalysis){
      if( fIsMC && fTRFM == AliAnalysisTaskPHOSPi0EtaToGammaGamma::kRFE && !fPHOSTriggerHelper->IsOnActiveTRUChannel(ph)) continue;//keep same TRU acceptance only in kRFE.
      if(!fIsMC && !ph->IsTrig()) continue;//it is meaningless to focus on photon without fired trigger in PHOS triggered data.
    }

    if(fForceActiveTRU && !fPHOSTriggerHelper->IsOnActiveTRUChannel(ph)) continue;//criterion fTRFM == kRFE is not needed.

    weight = 1.;
    if(fIsMC){
      primary = ph->GetPrimary();
      weight = ph->GetWeight();
    }

    pT = ph->Pt();
    energy = ph->Energy();
    phi = ph->Phi();

    if(fUseCoreEnergy){
      pT = (ph->GetMomV2())->Pt();
      energy = (ph->GetMomV2())->Energy();
      phi = (ph->GetMomV2())->Phi();
    }

    eff = f1tof->Eval(energy);
    if(!fIsMC && fIsPHOSTriggerAnalysis && fTRFM == AliAnalysisTaskPHOSPi0EtaToGammaGamma::kTAP){
      trgeff = f1trg->Eval(energy);
      AliInfo(Form("E = %4.3f GeV : trigger efficiency = %4.3f",energy,trgeff));
    }

    if(fIsMC || (!fIsMC && ph->IsTOFOK())){
                                             FillHistogramTH1(fOutputContainer,"hPhotonPt_noPID",pT,1/eff * weight * 1/trgeff);
      if(fPHOSClusterCuts->IsNeutral(ph))    FillHistogramTH1(fOutputContainer,"hPhotonPt_CPV"  ,pT,1/eff * weight * 1/trgeff);
      if(fPHOSClusterCuts->AcceptDisp(ph))   FillHistogramTH1(fOutputContainer,"hPhotonPt_Disp" ,pT,1/eff * weight * 1/trgeff);
      if(fPHOSClusterCuts->AcceptPhoton(ph)) FillHistogramTH1(fOutputContainer,"hPhotonPt_PID"  ,pT,1/eff * weight * 1/trgeff);
    }

    if(!fPHOSClusterCuts->AcceptPhoton(ph)) continue;

    //0 < photon phi < 2pi
    if(phi < 0) phi += TMath::TwoPi();
    TVector2 vg(TMath::Cos(fHarmonics * phi),TMath::Sin(fHarmonics * phi));

    if(fIsFlowTask){
      dphi = DeltaPhiIn0Pi(phi - fEventPlane);
      sp1 = vg * fQVector1;
      if(fFM == AliAnalysisTaskPHOSPi0EtaToGammaGamma::kEP)      value[1] = TMath::Cos(fHarmonics * dphi);
      else if(fFM == AliAnalysisTaskPHOSPi0EtaToGammaGamma::kSP) value[1] = sp1;
      else                                                       value[1] = 0;
    }
    else{
      //dphi = phi;
      //sp1 = 0;
      value[1] = 0;
    }

    if(fIsMC){
      if(IsFrom(primary,TrueK0SPt,310)  && IsPhoton(primary)) FillHistogramTH1(fOutputContainer,"hGammaFromK0S",pT,weight); 

      if(IsFrom(primary,TrueS0Pt ,3212) && IsPhoton(primary)) FillHistogramTH1(fOutputContainer,"hGammaFromS0" ,pT,weight);//Sigma0->L0 + gamma (Br = 100%)
      else if(IsFrom(primary,TrueL0Pt ,3122) && IsPhoton(primary)) FillHistogramTH1(fOutputContainer,"hGammaFromL0" ,pT,weight);//L0->n + pi0 (Br = 35.8%)
    }

    value[0] = pT;

    FillSparse(fOutputContainer,"hSparsePhoton",value,weight * 1/trgeff);

    if(ph->IsTOFOK()){
      FillSparse(fOutputContainer,"hSparsePhoton_TOF",value,1/eff * weight * 1/trgeff);
    }

  }//end of cluster loop

}
//________________________________________________________________________
void AliAnalysisTaskPHOSPi0EtaToGammaGamma::FillMgg() 
{
  const Int_t multClust = fPHOSClusterArray->GetEntriesFast();
  TLorentzVector p12, p12core;

  Double_t m12=0,pt12=0,asym=0;
  Double_t e1=0,e2=0;
  Double_t phi = -999, dphi = -999.;
  Double_t value[4] = {};
  Double_t sp1 = -999;

  Double_t eff1=1, eff2=1, eff12=1;
  TF1 *f1tof = GetTOFCutEfficiencyFunction();
  Double_t weight = 1., w1 = 1., w2 = 1.;

  Double_t trgeff1=1;
  Double_t trgeff2=1;
  Double_t trgeff12=1;
  TF1 *f1trg = GetTriggerEfficiencyFunction();

  Int_t primary1 = -1;
  Int_t primary2 = -1;
  Double_t TruePi0Pt = 0;
  Double_t TrueK0SPt = 0;
  Double_t TrueL0Pt = 0;
  Int_t commonID = -1;

  for(Int_t i1=0;i1<multClust-1;i1++){
    AliCaloPhoton *ph1 = (AliCaloPhoton*)fPHOSClusterArray->At(i1);
    if(!fPHOSClusterCuts->AcceptPhoton(ph1)) continue;
    if(!CheckMinimumEnergy(ph1)) continue;

    for(Int_t i2=i1+1;i2<multClust;i2++){
      AliCaloPhoton *ph2 = (AliCaloPhoton*)fPHOSClusterArray->At(i2);
      if(!fPHOSClusterCuts->AcceptPhoton(ph2)) continue;
      if(!CheckMinimumEnergy(ph2)) continue;

      if(fIsPHOSTriggerAnalysis){
        if(!fIsMC && (!ph1->IsTrig() && !ph2->IsTrig())) continue;//it is meaningless to reconstruct invariant mass with FALSE-FALSE combination in PHOS triggered data.
        if(fTRFM == AliAnalysisTaskPHOSPi0EtaToGammaGamma::kRFE && (!fPHOSTriggerHelper->IsOnActiveTRUChannel(ph1) || !fPHOSTriggerHelper->IsOnActiveTRUChannel(ph2))) continue;//use cluster pairs only on active TRU both in data and M.C.
      }

      if(fForceActiveTRU 
          && (!fPHOSTriggerHelper->IsOnActiveTRUChannel(ph1) || !fPHOSTriggerHelper->IsOnActiveTRUChannel(ph2))
        ) continue;//only for kINT7

      e1 = ph1->Energy();
      e2 = ph2->Energy();

      p12  = *ph1 + *ph2;
      m12  = p12.M();
      pt12 = p12.Pt();
      phi  = p12.Phi();
      asym = TMath::Abs((ph1->Energy()-ph2->Energy())/(ph1->Energy()+ph2->Energy()));

      if(fUseCoreEnergy){
        p12core = *(ph1->GetMomV2()) + *(ph2->GetMomV2());
        m12     = p12core.M();
        pt12    = p12core.Pt();
        phi     = p12core.Phi();

        e1 = (ph1->GetMomV2())->Energy();
        e2 = (ph2->GetMomV2())->Energy();
        asym = TMath::Abs(e1 - e2) / (e1 + e2);
      }

      eff1 = f1tof->Eval(e1);
      eff2 = f1tof->Eval(e2);
      eff12 = eff1 * eff2;

      if(!fIsMC && fIsPHOSTriggerAnalysis && fTRFM == AliAnalysisTaskPHOSPi0EtaToGammaGamma::kTAP){
        if((!ph1->IsTrig() || e1 < fEnergyThreshold)
        && (!ph2->IsTrig() || e2 < fEnergyThreshold)
          ) continue;

        trgeff1  = f1trg->Eval(e1);
        trgeff2  = f1trg->Eval(e2);
        trgeff12 = trgeff1 + trgeff2 - (trgeff1 * trgeff2);//logical OR//this is true only when occupancy is uniformed.
      }

      weight = 1.;
      if(fIsMC){
        w1= ph1->GetWeight();
        primary1 = ph1->GetPrimary();

        w2 = ph2->GetWeight();
        primary2 = ph2->GetPrimary();

        commonID = FindCommonParent(primary1,primary2);
        if(commonID > -1) AliInfo(Form("commonID is found. primary1 = %d, primary2 = %d, commonID = %d, w1 = %e , w2 = %e",primary1,primary2,commonID,w1,w2));

        if(commonID > -1) weight = w1;
        else weight = w1*w2;

        if(IsFrom(primary1,TruePi0Pt,111) && IsFrom(primary2,TruePi0Pt,111) && commonID > -1){
          AliAODMCParticle *p1 = (AliAODMCParticle*)fMCArrayAOD->At(primary1);
          AliAODMCParticle *p2 = (AliAODMCParticle*)fMCArrayAOD->At(primary2);
          Int_t pdg1 = p1->PdgCode(); 
          Int_t pdg2 = p2->PdgCode(); 

          AliAODMCParticle *p = (AliAODMCParticle*)fMCArrayAOD->At(commonID);
          Int_t pdg = p->PdgCode(); 
          if(pdg1 == 22 && pdg2 == 22 && pdg == 111){
            FillHistogramTH2(fOutputContainer,"hMggFromPi0",m12,pt12,weight);
            if(asym < 0.8) FillHistogramTH2(fOutputContainer,"hMggFromPi0_asym08",m12,pt12,weight);
          }
        }
        if(IsFrom(primary1,TrueK0SPt,310) && IsFrom(primary2,TrueK0SPt,310) && commonID > -1){
          //for feed down correction from K0S->pi0 + pi0
          AliAODMCParticle *p1 = (AliAODMCParticle*)fMCArrayAOD->At(primary1);
          AliAODMCParticle *p2 = (AliAODMCParticle*)fMCArrayAOD->At(primary2);
          Int_t pdg1 = p1->PdgCode(); 
          Int_t pdg2 = p2->PdgCode(); 
          AliAODMCParticle *p = (AliAODMCParticle*)fMCArrayAOD->At(commonID);
          Int_t pdg = p->PdgCode(); 
          if(pdg1 == 22 && pdg2 == 22 && pdg == 111){
            FillHistogramTH2(fOutputContainer,"hMggFromK0S",m12,pt12,weight);
            if(asym < 0.8) FillHistogramTH2(fOutputContainer,"hMggFromK0S_asym08",m12,pt12,weight);
          }
        }
        if(IsFrom(primary1,TrueL0Pt,3122) && IsFrom(primary2,TrueL0Pt,3122) && commonID > -1){
          //for feed down correction from L0->pi0 + neutron
          AliAODMCParticle *p1 = (AliAODMCParticle*)fMCArrayAOD->At(primary1);
          AliAODMCParticle *p2 = (AliAODMCParticle*)fMCArrayAOD->At(primary2);
          Int_t pdg1 = p1->PdgCode(); 
          Int_t pdg2 = p2->PdgCode(); 
          AliAODMCParticle *p = (AliAODMCParticle*)fMCArrayAOD->At(commonID);
          Int_t pdg = p->PdgCode(); 
          if(pdg1 == 22 && pdg2 == 22 && pdg == 111){
            FillHistogramTH2(fOutputContainer,"hMggFromLambda0",m12,pt12,weight);
            if(asym < 0.8) FillHistogramTH2(fOutputContainer,"hMggFromLambda0_asym08",m12,pt12,weight);
          }
        }

      }//end of if fIsMC

      if(phi < 0) phi += TMath::TwoPi();

      TVector2 vgg(TMath::Cos(fHarmonics * phi),TMath::Sin(fHarmonics * phi));
      if(fIsFlowTask){
        dphi = DeltaPhiIn0Pi(phi - fEventPlane);
        sp1 = vgg * fQVector1;
        if(fFM == AliAnalysisTaskPHOSPi0EtaToGammaGamma::kEP)      value[3] = TMath::Cos(fHarmonics * dphi);
        else if(fFM == AliAnalysisTaskPHOSPi0EtaToGammaGamma::kSP) value[3] = sp1;
        else                                                       value[3] = 0;
      }
      else{
        dphi = phi;
        sp1 = 0;
        value[3] = 0;
      }

      value[0] = m12;
      value[1] = pt12;
      value[2] = asym;

      if(fIsOAStudy){
        Double_t oa = TMath::Abs(ph1->Angle(ph2->Vect())) * 1e+3;//rad->mrad
        FillHistogramTH3(fOutputContainer,"hMgg_OA",m12,pt12,oa,weight);
      }
      if(m12 > 0.96) continue;//reduce entry in THnSparse

      if(TMath::Abs(ph1->Module()-ph2->Module()) < 2) FillHistogramTH2(fOutputContainer,Form("hMgg_M%d%d",TMath::Min(ph1->Module(),ph2->Module()), TMath::Max(ph1->Module(),ph2->Module())),m12,pt12,weight * 1/trgeff12);
      FillSparse(fOutputContainer,"hSparseMgg",value,weight * 1/trgeff12);

      if(ph1->IsTOFOK() && ph2->IsTOFOK()){
        FillSparse(fOutputContainer,"hSparseMgg_TOF",value,1/eff12 * weight * 1/trgeff12);

        if(TMath::Abs(ph1->Module()-ph2->Module()) < 2) FillHistogramTH2(fOutputContainer,Form("hMgg_M%d%d_TOF",TMath::Min(ph1->Module(),ph2->Module()), TMath::Max(ph1->Module(),ph2->Module())),m12,pt12,1/eff12 * weight * 1/trgeff12);

      }//end of TOF cut


    }//end of ph2

  }//end of ph1

}
//________________________________________________________________________
void AliAnalysisTaskPHOSPi0EtaToGammaGamma::FillMixMgg() 
{
  TList *prevPHOS = fPHOSEvents[fZvtx][fEPBin];

  const Int_t multClust = fPHOSClusterArray->GetEntriesFast();

  TLorentzVector p12, p12core;
  Double_t m12=0,pt12=0,asym=0;
  Double_t phi = -999, dphi = -999.;
  Double_t weight = 1., w1 = 1., w2 = 1.;
  Double_t value[4] = {};
  Double_t sp1 = -999;

  Double_t eff1=1, eff2=1, eff12=1;
  Double_t e1=0,e2=0;
  TF1 *f1tof = GetTOFCutEfficiencyFunction();

  Double_t trgeff1=1;
  Double_t trgeff2=1;
  Double_t trgeff12=1;
  TF1 *f1trg = GetTriggerEfficiencyFunction();

  for(Int_t i1=0;i1<multClust;i1++){
    AliCaloPhoton *ph1 = (AliCaloPhoton*)fPHOSClusterArray->At(i1);
    if(!fPHOSClusterCuts->AcceptPhoton(ph1)) continue;
    if(!CheckMinimumEnergy(ph1)) continue;

    for(Int_t ev=0;ev<prevPHOS->GetSize();ev++){
      TClonesArray *mixPHOS = static_cast<TClonesArray*>(prevPHOS->At(ev));

      for(Int_t i2=0;i2<mixPHOS->GetEntriesFast();i2++){
        AliCaloPhoton *ph2 = (AliCaloPhoton*)mixPHOS->At(i2);
        if(!fPHOSClusterCuts->AcceptPhoton(ph2)) continue;
        if(!CheckMinimumEnergy(ph2)) continue;

        if(fIsPHOSTriggerAnalysis){
          if(!fIsMC && (!ph1->IsTrig() && !ph2->IsTrig())) continue;//it is meaningless to reconstruct invariant mass with FALSE-FALSE combination in PHOS triggered data.
          if(fTRFM == AliAnalysisTaskPHOSPi0EtaToGammaGamma::kRFE && (!fPHOSTriggerHelper->IsOnActiveTRUChannel(ph1) || !fPHOSTriggerHelper->IsOnActiveTRUChannel(ph2))) continue;//use cluster pairs only on active TRU both in data and M.C.
        }

        if(fForceActiveTRU 
            && (!fPHOSTriggerHelper->IsOnActiveTRUChannel(ph1) || !fPHOSTriggerHelper->IsOnActiveTRUChannel(ph2))
          ) continue;//only for kINT7

        e1 = ph1->Energy();
        e2 = ph2->Energy();

        p12  = *ph1 + *ph2;
        m12  = p12.M();
        pt12 = p12.Pt();
        phi  = p12.Phi();
        asym = TMath::Abs((ph1->Energy()-ph2->Energy())/(ph1->Energy()+ph2->Energy()));

        if(fUseCoreEnergy){
          p12core = *(ph1->GetMomV2()) + *(ph2->GetMomV2());
          m12     = p12core.M();
          pt12    = p12core.Pt();
          phi     = p12core.Phi();

          e1 = (ph1->GetMomV2())->Energy();
          e2 = (ph2->GetMomV2())->Energy();
          asym = TMath::Abs(e1 - e2) / (e1 + e2);

        }

        eff1 = f1tof->Eval(e1);
        eff2 = f1tof->Eval(e2);
        eff12 = eff1 * eff2;
        weight = 1.;

        if(!fIsMC && fIsPHOSTriggerAnalysis && fTRFM == AliAnalysisTaskPHOSPi0EtaToGammaGamma::kTAP){
          if((!ph1->IsTrig() || e1 < fEnergyThreshold)
          && (!ph2->IsTrig() || e2 < fEnergyThreshold)
            ) continue;

          trgeff1  = f1trg->Eval(e1);
          trgeff2  = f1trg->Eval(e2);
          trgeff12 = trgeff1 + trgeff2 - (trgeff1 * trgeff2);//logical OR

        }

        if(fIsMC){
          w1= ph1->GetWeight();
          w2 = ph2->GetWeight();

          weight = w1*w2;

        }//end of if fIsMC

        if(phi < 0) phi += TMath::TwoPi();

        TVector2 vgg(TMath::Cos(fHarmonics * phi),TMath::Sin(fHarmonics * phi));
        if(fIsFlowTask){
          dphi = DeltaPhiIn0Pi(phi - fEventPlane);
          sp1 = vgg * fQVector1;
          if(fFM == AliAnalysisTaskPHOSPi0EtaToGammaGamma::kEP)      value[3] = TMath::Cos(fHarmonics * dphi);
          else if(fFM == AliAnalysisTaskPHOSPi0EtaToGammaGamma::kSP) value[3] = sp1;
          else                                                       value[3] = 0;
        }
        else{
          dphi = phi;
          sp1 = 0;
          value[3] = 0;
        }
        value[0] = m12;
        value[1] = pt12;
        value[2] = asym;

        if(fIsOAStudy){
          Double_t oa = TMath::Abs(ph1->Angle(ph2->Vect())) * 1e+3;//rad->mrad
          FillHistogramTH3(fOutputContainer,"hMixMgg_OA",m12,pt12,oa,weight);
        }

        if(m12 > 0.96) continue;//reduce entry in THnSparse

        FillSparse(fOutputContainer,"hSparseMixMgg",value,weight * 1/trgeff12);

        if(TMath::Abs(ph1->Module()-ph2->Module()) < 2) FillHistogramTH2(fOutputContainer,Form("hMixMgg_M%d%d",TMath::Min(ph1->Module(),ph2->Module()), TMath::Max(ph1->Module(),ph2->Module())),m12,pt12, weight * 1/trgeff12);

        if(ph1->IsTOFOK() && ph2->IsTOFOK()){
          FillSparse(fOutputContainer,"hSparseMixMgg_TOF",value,1/eff12 * weight * 1/trgeff12);

          if(TMath::Abs(ph1->Module()-ph2->Module()) < 2) FillHistogramTH2(fOutputContainer,Form("hMixMgg_M%d%d_TOF",TMath::Min(ph1->Module(),ph2->Module()), TMath::Max(ph1->Module(),ph2->Module())),m12,pt12,1/eff12 * weight * 1/trgeff12);

        }//end of TOF cut

      }//end of ph2

    }//end of mix

  }//end of ph1

}
//________________________________________________________________________
void AliAnalysisTaskPHOSPi0EtaToGammaGamma::FillM3pi() 
{
  const Int_t trackMult = fEvent->GetNumberOfTracks();

  TObjArray *array_pip = new TObjArray(trackMult);
  TObjArray *array_pim = new TObjArray(trackMult);
  Double_t trackPx = 0;
  Double_t trackPy = 0;
  Double_t trackPz = 0;
  Double_t trackE  = 0;
  const Double_t Mpi  = 0.13957018;//in GeV/c2
  const Double_t Mpi0 = 0.1349766;//in GeV/c2
  //Double_t nsigmaPion = 999;

  for(Int_t itrack=0;itrack<trackMult;itrack++){
    AliAODTrack *aodtrack = (AliAODTrack*)fEvent->GetTrack(itrack);

    if(TMath::Abs(aodtrack->Eta()) > 0.8) continue;
    if(aodtrack->Pt() < fMinPtChPi) continue;

    //if(!aodtrack->IsHybridGlobalConstrainedGlobal()) continue; //select only hybrid track//loose DCA cut.
    if(!aodtrack->TestFilterMask(AliAODTrack::kTrkGlobal)) continue;//standard cuts with tight DCA cut (AliESDtrackCuts::GetStandardITSTPCTrackCuts2011())
    //if(!aodtrack->TestFilterMask(AliAODTrack::kTrkGlobalNoDCA)) continue;//standard cuts with very loose DCA cut
   
    //nsigmaPion = fPIDResponse->NumberOfSigmasTPC(dynamic_cast<AliVTrack*>(aodtrack),AliPID::kPion);
    //if(nsigmaPion < -3 || 3 < nsigmaPion) continue;

    if(TMath::Abs(aodtrack->Charge()) != 1) continue;//select pi+ / pi-//reject deuteron, triton
 
    trackPx = aodtrack->Px();
    trackPy = aodtrack->Py();
    trackPz = aodtrack->Pz();
    trackE  = TMath::Sqrt(TMath::Power(trackPx,2) + TMath::Power(trackPy,2) + TMath::Power(trackPz,2) + TMath::Power(Mpi,2));
    TLorentzVector *p1track = new TLorentzVector(trackPx,trackPy,trackPz,trackE);

    if(aodtrack->Charge() > 0){//pi+
      array_pip->Add(p1track);
    }
    else{//pi-
      array_pim->Add(p1track);
    }

  }//end of track loop

  const Int_t Npip = array_pip->GetEntriesFast();
  const Int_t Npim = array_pim->GetEntriesFast();
  AliInfo(Form("Npip = %d , Npim = %d",Npip,Npim));

  const Int_t multClust = fPHOSClusterArray->GetEntriesFast();
  TLorentzVector p12, p12core;

  Double_t m12=0,pt12=0,asym=0;
  Double_t pi0Px = 0;
  Double_t pi0Py = 0;
  Double_t pi0Pz = 0;
  Double_t e1=0,e2=0;
  Double_t phi = -999, dphi = -999.;
  Double_t value[4] = {};
  Double_t sp1 = -999;

  Double_t eff1=1, eff2=1, eff12=1;
  TF1 *f1tof = GetTOFCutEfficiencyFunction();
  Double_t weight = 1., w1 = 1., w2 = 1.;

  Double_t trgeff1=1;
  Double_t trgeff2=1;
  Double_t trgeff12=1;
  TF1 *f1trg = GetTriggerEfficiencyFunction();

  Int_t primary1 = -1;
  Int_t primary2 = -1;
  Int_t commonID = -1;

  Double_t m3pi=0,pt3pi=0;
  Double_t phipi0 = 0, etapi0 = 0;
  Double_t phipip = 0, etapip = 0;
  Double_t phipim = 0, etapim = 0;
  Double_t dphi_pp = 999, deta_pp = 999, dR = 999;

  for(Int_t i1=0;i1<multClust-1;i1++){
    AliCaloPhoton *ph1 = (AliCaloPhoton*)fPHOSClusterArray->At(i1);
    if(!fPHOSClusterCuts->AcceptPhoton(ph1)) continue;
    if(!CheckMinimumEnergy(ph1)) continue;

    for(Int_t i2=i1+1;i2<multClust;i2++){
      AliCaloPhoton *ph2 = (AliCaloPhoton*)fPHOSClusterArray->At(i2);
      if(!fPHOSClusterCuts->AcceptPhoton(ph2)) continue;
      if(!CheckMinimumEnergy(ph2)) continue;

      if(fIsPHOSTriggerAnalysis){
        if(!fIsMC && (!ph1->IsTrig() && !ph2->IsTrig())) continue;//it is meaningless to reconstruct invariant mass with FALSE-FALSE combination in PHOS triggered data.
        if(fTRFM == AliAnalysisTaskPHOSPi0EtaToGammaGamma::kRFE && (!fPHOSTriggerHelper->IsOnActiveTRUChannel(ph1) || !fPHOSTriggerHelper->IsOnActiveTRUChannel(ph2))) continue;//use cluster pairs only on active TRU both in data and M.C.
      }

      if(fForceActiveTRU 
          && (!fPHOSTriggerHelper->IsOnActiveTRUChannel(ph1) || !fPHOSTriggerHelper->IsOnActiveTRUChannel(ph2))
        ) continue;//only for kINT7

      p12  = *ph1 + *ph2;
      m12  = p12.M();
      pt12 = p12.Pt();
      pi0Px = p12.Px();
      pi0Py = p12.Py();
      pi0Pz = p12.Pz();

      e1 = ph1->Energy();
      e2 = ph2->Energy();
      asym = TMath::Abs((ph1->Energy()-ph2->Energy())/(ph1->Energy()+ph2->Energy()));

      if(fUseCoreEnergy){
        p12core = *(ph1->GetMomV2()) + *(ph2->GetMomV2());
        m12     = p12core.M();
        pt12    = p12core.Pt();

        pi0Px = p12core.Px();
        pi0Py = p12core.Py();
        pi0Pz = p12core.Pz();

        e1 = (ph1->GetMomV2())->Energy();
        e2 = (ph2->GetMomV2())->Energy();
        asym = TMath::Abs(e1 - e2) / (e1 + e2);
      }

      eff1 = f1tof->Eval(e1);
      eff2 = f1tof->Eval(e2);
      eff12 = eff1 * eff2;

      if(!fIsMC && fIsPHOSTriggerAnalysis && fTRFM == AliAnalysisTaskPHOSPi0EtaToGammaGamma::kTAP){
        trgeff1  = f1trg->Eval(e1);
        trgeff2  = f1trg->Eval(e2);
        trgeff12 = trgeff1 + trgeff2 - (trgeff1 * trgeff2);//logical OR//this is true only when occupancy is uniformed.
      }

      if(m12 < 0.12 || 0.15 < m12) continue;
      if(pt12 < fMinPtPi0) continue;

      Double_t pi0E = TMath::Sqrt(TMath::Power(pi0Px,2) + TMath::Power(pi0Py,2) + TMath::Power(pi0Pz,2) + TMath::Power(Mpi0,2));
      TLorentzVector *p1pi0 = new TLorentzVector(pi0Px,pi0Py,pi0Pz,pi0E);
      phipi0 = p1pi0->Phi();
      if(phipi0 < 0) phipi0 += TMath::TwoPi();
      etapi0 = p1pi0->Eta();

      for(Int_t itrack1=0;itrack1<Npip;itrack1++){
        TLorentzVector *p1pip = (TLorentzVector*)array_pip->At(itrack1);
        phipip = p1pip->Phi();
        if(phipip < 0) phipip += TMath::TwoPi();
        etapip = p1pip->Eta();

        dphi_pp = phipip - phipi0;
        if(dphi_pp >  TMath::Pi()) dphi_pp -= TMath::TwoPi();
        if(dphi_pp < -TMath::Pi()) dphi_pp += TMath::TwoPi();
        deta_pp = etapip - etapi0;
        dR = TMath::Sqrt(TMath::Power(dphi_pp,2) + TMath::Power(deta_pp,2));
        if(dR > fMaxR) continue;

        //printf("pip | px = %e , py = %e , pz = %e\n",p1pip->Px(),p1pip->Py(),p1pip->Pz());

        for(Int_t itrack2=0;itrack2<Npim;itrack2++){
          TLorentzVector *p1pim = (TLorentzVector*)array_pim->At(itrack2);
          phipim = p1pim->Phi();
          if(phipim < 0) phipim += TMath::TwoPi();
          etapim = p1pim->Eta();

          dphi_pp = phipim - phipi0;
          if(dphi_pp >  TMath::Pi()) dphi_pp -= TMath::TwoPi();
          if(dphi_pp < -TMath::Pi()) dphi_pp += TMath::TwoPi();
          deta_pp = etapim - etapi0;
          dR = TMath::Sqrt(TMath::Power(dphi_pp,2) + TMath::Power(deta_pp,2));
          if(dR > fMaxR) continue;

          //printf("pim | px = %e , py = %e , pz = %e\n",p1pim->Px(),p1pim->Py(),p1pim->Pz());
          //printf("Mpi0 = %e , Mpip = %e , Mpim = %e\n",p1pi0->M(),p1pip->M(),p1pim->M());

          TLorentzVector p3pi = *p1pi0 + *p1pip + *p1pim;

          m3pi  = p3pi.M();
          pt3pi = p3pi.Pt();

          if(TMath::Abs(p3pi.Y()) > 0.5) continue;//measure omega at |rapidity| < 0.5

          weight = 1.;
          if(fIsMC){
            w1= ph1->GetWeight();
            primary1 = ph1->GetPrimary();

            w2 = ph2->GetWeight();
            primary2 = ph2->GetPrimary();

            commonID = FindCommonParent(primary1,primary2);
            if(commonID > -1) weight = w1;
            else weight = w1*w2;

          }//end of if fIsMC

          phi = p3pi.Phi();
          if(phi < 0) phi += TMath::TwoPi();

          TVector2 vgg(TMath::Cos(fHarmonics * phi),TMath::Sin(fHarmonics * phi));
          if(fIsFlowTask){
            dphi = DeltaPhiIn0Pi(phi - fEventPlane);
            sp1 = vgg * fQVector1;
            if(fFM == AliAnalysisTaskPHOSPi0EtaToGammaGamma::kEP)      value[3] = TMath::Cos(fHarmonics * dphi);
            else if(fFM == AliAnalysisTaskPHOSPi0EtaToGammaGamma::kSP) value[3] = sp1;
            else                                                       value[3] = 0;
          }
          else{
            dphi = phi;
            sp1 = 0;
            value[3] = 0;
          }

          value[0] = m3pi;
          value[1] = pt3pi;
          value[2] = asym;

          FillSparse(fOutputContainer,"hSparseM3pi",value,weight * 1/trgeff12);

          if(ph1->IsTOFOK() && ph2->IsTOFOK()){
            FillSparse(fOutputContainer,"hSparseM3pi_TOF",value,1/eff12 * weight * 1/trgeff12);

          }//end of TOF cut

        }//end of track loop2

      }//end of track loop1

      delete p1pi0;
      p1pi0 = 0x0;

    }//end of ph2

  }//end of ph1

  array_pip->Clear();
  array_pim->Clear();

  delete array_pip;
  array_pip = 0x0;
  delete array_pim;
  array_pim = 0x0;

}
//________________________________________________________________________
void AliAnalysisTaskPHOSPi0EtaToGammaGamma::FillEpRatio() 
{
  //since this E/p is for global energy scale and mis-alignment study, neither TOF cut nor additional weight in M.C. are applied.

  const Int_t multClust = fPHOSClusterArray->GetEntriesFast();

  const Double_t NsigmaCPV = 2.0;
  const Double_t NsigmaDisp = fPHOSClusterCuts->GetDispParameter();

  AliPHOSClusterCuts *cuts = new AliPHOSClusterCuts("CutsForCharged");//for charged tracks
  cuts->SetUseCoreDispersion(fPHOSClusterCuts->IsCoreDisp());
  cuts->SetNsigmaCPV(NsigmaCPV);
  cuts->SetNsigmaDisp(NsigmaDisp);

  //Double_t eff = 1;
  //TF1 *f1tof = GetTOFCutEfficiencyFunction();

  TVector3 localPos;
  Int_t module = 0;
  Double_t trackP=0, trackPt=0;
  Double_t dEdx=0;
  Double_t energy=0;
  Int_t trackC = 0;
  Double_t trackDx=0, trackDz=0;
  Float_t position[3] = {};

  AliVCluster *cluster=0x0;
  AliVTrack *track = 0x0;
  Double_t nsigmaElectron = 999;
  Bool_t PID_ele = kFALSE;
  Bool_t isGlobal = kFALSE;
  //Double_t weight = 1.;

  for(Int_t iph=0;iph<multClust;iph++){
    AliCaloPhoton *ph = (AliCaloPhoton*)fPHOSClusterArray->At(iph);
    if(!cuts->AcceptChargedParticle(ph)) continue;
    if(!CheckMinimumEnergy(ph)) continue;

    energy = ph->Energy();
    if(fUseCoreEnergy) energy = (ph->GetMomV2())->Energy();

    //weight = 1.;
    //if(fIsMC) weight = ph->GetWeight();

    //eff = 1.;
    //if(!fIsMC) eff = f1tof->Eval(energy);

    cluster = (AliVCluster*)ph->GetCluster();
    PID_ele = kFALSE;
    isGlobal = kFALSE;
    track = 0x0;

    if(fESDEvent){
      Int_t trackindex = cluster->GetTrackMatchedIndex();
      if(trackindex > 0){
        track = (AliVTrack*)(fEvent->GetTrack(trackindex));
        isGlobal = fESDtrackCutsGlobal->AcceptTrack(dynamic_cast<AliESDtrack*>(track));
      }
    }//end of ESD
    else if(fAODEvent){
      if(cluster->GetNTracksMatched() > 0){
        track = dynamic_cast<AliVTrack*>(cluster->GetTrackMatched(0));
        isGlobal = dynamic_cast<AliAODTrack*>(track)->TestFilterMask(AliAODTrack::kTrkGlobalNoDCA);//standard cuts with very loose DCA cut

      }//end of track matching
    }//end of AOD

    if(track && isGlobal){
      trackP  = track->P();
      trackPt = track->Pt();
      dEdx = track->GetTPCsignal();
      trackC  = track->Charge();

      trackDx = cluster->GetTrackDx();
      trackDz = cluster->GetTrackDz();
      module = ph->Module();

      position[0] = ph->EMCx();
      position[1] = ph->EMCy();
      position[2] = ph->EMCz();

      TVector3 global1(position);
      fPHOSGeo->Global2Local(localPos,global1,module);

      FillHistogramTH3(fOutputContainer,Form("hdZvsZvsTrackPt_M%d",module),localPos.Z(),trackDz,trackPt);//tracks bended to only phi direction

      if(trackC > 0) FillHistogramTH3(fOutputContainer,Form("hdXvsXvsTrackPt_plus_M%d",module) ,localPos.X(),trackDx,trackPt);//positive
      else           FillHistogramTH3(fOutputContainer,Form("hdXvsXvsTrackPt_minus_M%d",module),localPos.X(),trackDx,trackPt);//negative

      if(!cuts->AcceptElectron(ph)) continue;

      nsigmaElectron = fPIDResponse->NumberOfSigmasTPC(track,AliPID::kElectron);
      FillHistogramTH2(fOutputContainer,"hEpRatiovsNsigmaElectronTPC",energy/trackP,nsigmaElectron);
      if(-2 < nsigmaElectron && nsigmaElectron < 3) PID_ele = kTRUE;

      if(PID_ele){
        FillHistogramTH2(fOutputContainer,Form("hEpRatiovsEnergy_M%d_Electron" ,module),energy/trackP,energy);
        FillHistogramTH2(fOutputContainer,Form("hEpRatiovsTrackPt_M%d_Electron",module),energy/trackP,trackPt);
        FillHistogramTH2(fOutputContainer,"hTPCdEdx_Electron",trackP,dEdx);
        if(0.8 < energy/trackP && energy/trackP < 1.2) FillHistogramTH3(fOutputContainer,Form("hdZvsZvsTrackPtElectron_M%d",module),localPos.Z(),trackDz,trackPt);//for alignment study
      }
      if(nsigmaElectron < -3 || 5 < nsigmaElectron){//non-electron and far from border
        FillHistogramTH2(fOutputContainer,Form("hEpRatiovsEnergy_M%d_Others" ,module),energy/trackP,energy);
        FillHistogramTH2(fOutputContainer,Form("hEpRatiovsTrackPt_M%d_Others",module),energy/trackP,trackPt);
        FillHistogramTH2(fOutputContainer,"hTPCdEdx_Others",trackP,dEdx);
      }

    }//end of track matching

  }//end of cluster loop

  delete cuts;
  cuts = 0x0;

}
//________________________________________________________________________
void AliAnalysisTaskPHOSPi0EtaToGammaGamma::EstimatePIDCutEfficiency()
{
  //tag and probe method is used.

  const Int_t multClust = fPHOSClusterArray->GetEntriesFast();

  TLorentzVector p12, p12core;
  Double_t m12=0;
  Double_t pT=0;
  Double_t energy=0;
  Double_t weight = 1., w1 = 1., w2 = 1.;
  Int_t primary1 = -1;
  Int_t primary2 = -1;
  Int_t commonID = -1;

  Double_t value[5] = {0,0,0,0,0};

  for(Int_t i1=0;i1<multClust;i1++){
    AliCaloPhoton *ph1 = (AliCaloPhoton*)fPHOSClusterArray->At(i1);
    if(!fIsMC && fIsPHOSTriggerAnalysis && !ph1->IsTrig()) continue;//take trigger bias into account.

    if(!CheckMinimumEnergy(ph1)) continue;
    //if(!fPHOSClusterCuts->AcceptPhoton(ph1)) continue;

    //apply tight cut to photon1
    if(ph1->Energy() < 0.5 || ph1->GetNsigmaCPV() < 4 || ph1->GetNsigmaCoreDisp() > 2.5) continue;

    for(Int_t i2=0;i2<multClust;i2++){
      AliCaloPhoton *ph2 = (AliCaloPhoton*)fPHOSClusterArray->At(i2);
      if(!CheckMinimumEnergy(ph2)) continue;

      if(i2==i1) continue;//reject same cluster combination

      p12 = *ph1 + *ph2;
      m12 = p12.M();
      pT = ph2->Pt();
      energy = ph2->Energy();

      if(fUseCoreEnergy){
        p12core = *(ph1->GetMomV2()) + *(ph2->GetMomV2());
        m12 = p12core.M();
        pT = (ph2->GetMomV2())->Pt();
        energy = (ph2->GetMomV2())->Energy();
      }

      weight = 1.;
      if(fIsMC){
        w1= ph1->GetWeight();
        primary1 = ph1->GetPrimary();

        w2 = ph2->GetWeight();
        primary2 = ph2->GetPrimary();

        commonID = FindCommonParent(primary1,primary2);
        if(commonID > -1) weight = w1;
        else weight = w1*w2;

      }//end of if fIsMC


      value[0] = m12;
      value[1] = energy;
      value[2] = ph2->GetNCells();
      value[3] = ph2->GetLambda1();
      value[4] = ph2->GetLambda2();
      if(fPIDStudy) FillSparse(fOutputContainer,"hSparseMgg_PID",value,weight);

      FillHistogramTH2(fOutputContainer,"hMgg_Probe_PID",m12,pT,weight);

      if(fPHOSClusterCuts->IsNeutral(ph2))    FillHistogramTH2(fOutputContainer,"hMgg_PassingProbe_CPV" ,m12,pT,weight);
      if(fPHOSClusterCuts->AcceptDisp(ph2))   FillHistogramTH2(fOutputContainer,"hMgg_PassingProbe_Disp",m12,pT,weight);
      if(fPHOSClusterCuts->AcceptPhoton(ph2)) FillHistogramTH2(fOutputContainer,"hMgg_PassingProbe_PID" ,m12,pT,weight);

    }//end of ph2

  }//end of ph1

  //next mixed event
  TList *prevPHOS = fPHOSEvents[fZvtx][fEPBin];

  for(Int_t i1=0;i1<multClust;i1++){
    AliCaloPhoton *ph1 = (AliCaloPhoton*)fPHOSClusterArray->At(i1);
    if(!fIsMC && fIsPHOSTriggerAnalysis && !ph1->IsTrig()) continue;//take trigger bias into account.

    if(!CheckMinimumEnergy(ph1)) continue;
    //if(!fPHOSClusterCuts->AcceptPhoton(ph1)) continue;

    //apply tight cut to photon1
    if(ph1->Energy() < 0.5 || ph1->GetNsigmaCPV() < 4 || ph1->GetNsigmaCoreDisp() > 2.5) continue;

    for(Int_t ev=0;ev<prevPHOS->GetSize();ev++){
      TClonesArray *mixPHOS = static_cast<TClonesArray*>(prevPHOS->At(ev));

      for(Int_t i2=0;i2<mixPHOS->GetEntriesFast();i2++){
        AliCaloPhoton *ph2 = (AliCaloPhoton*)mixPHOS->At(i2);
        if(!CheckMinimumEnergy(ph2)) continue;

        p12 = *ph1 + *ph2;
        m12 = p12.M();
        pT = ph2->Pt();
        energy = ph2->Energy();

        if(fUseCoreEnergy){
          p12core = *(ph1->GetMomV2()) + *(ph2->GetMomV2());
          m12 = p12core.M();
          pT = (ph2->GetMomV2())->Pt();
          energy = (ph2->GetMomV2())->Energy();
        }

        weight = 1.;
        if(fIsMC){
          w1= ph1->GetWeight();
          w2 = ph2->GetWeight();
          weight = w1*w2;
        }//end of if fIsMC

        value[0] = m12;
        value[1] = energy;
        value[2] = ph2->GetNCells();
        value[3] = ph2->GetLambda1();
        value[4] = ph2->GetLambda2();
        if(fPIDStudy) FillSparse(fOutputContainer,"hSparseMixMgg_PID",value,weight);

        FillHistogramTH2(fOutputContainer,"hMixMgg_Probe_PID",m12,pT,weight);
        if(fPHOSClusterCuts->IsNeutral(ph2))    FillHistogramTH2(fOutputContainer,"hMixMgg_PassingProbe_CPV" ,m12,pT,weight);
        if(fPHOSClusterCuts->AcceptDisp(ph2))   FillHistogramTH2(fOutputContainer,"hMixMgg_PassingProbe_Disp",m12,pT,weight);
        if(fPHOSClusterCuts->AcceptPhoton(ph2)) FillHistogramTH2(fOutputContainer,"hMixMgg_PassingProbe_PID" ,m12,pT,weight);

      }//end of mix

    }//end of ph2

  }//end of ph1

}
//________________________________________________________________________
void AliAnalysisTaskPHOSPi0EtaToGammaGamma::EstimateTOFCutEfficiency()
{
  //tag and probe method is used.

  const Int_t multClust = fPHOSClusterArray->GetEntriesFast();

  TLorentzVector p12, p12core;
  Double_t m12=0;
  Double_t energy=0;

  for(Int_t i1=0;i1<multClust;i1++){
    AliCaloPhoton *ph1 = (AliCaloPhoton*)fPHOSClusterArray->At(i1);
    if(!fPHOSClusterCuts->AcceptPhoton(ph1)) continue;
    if(!ph1->IsTOFOK()) continue;
    if(!CheckMinimumEnergy(ph1)) continue;

    if(!fIsMC && fIsPHOSTriggerAnalysis && !ph1->IsTrig()) continue;//take trigger bias into account.
    //triggered clusters have very likely their hit timing within 12.5 ns in kPHI7 events.
 
    for(Int_t i2=0;i2<multClust;i2++){
      AliCaloPhoton *ph2 = (AliCaloPhoton*)fPHOSClusterArray->At(i2);
      if(!fPHOSClusterCuts->AcceptPhoton(ph2)) continue;
      if(!CheckMinimumEnergy(ph2)) continue;

      if(i2==i1) continue;//reject same cluster combination

      p12 = *ph1 + *ph2;
      m12 = p12.M();
      energy = ph2->Energy();

      if(fUseCoreEnergy){
        p12core = *(ph1->GetMomV2()) + *(ph2->GetMomV2());
        m12 = p12core.M();
        energy = (ph2->GetMomV2())->Energy();
      }

      FillHistogramTH2(fOutputContainer,"hMgg_Probe_TOF",m12,energy);
      if(ph2->IsTOFOK()) FillHistogramTH2(fOutputContainer,"hMgg_PassingProbe_TOF",m12,energy);

    }//end of ph2

  }//end of ph1

  //next mixed event
  TList *prevPHOS = fPHOSEvents[fZvtx][fEPBin];

  for(Int_t i1=0;i1<multClust;i1++){
    AliCaloPhoton *ph1 = (AliCaloPhoton*)fPHOSClusterArray->At(i1);
    if(!fPHOSClusterCuts->AcceptPhoton(ph1)) continue;
    if(!ph1->IsTOFOK()) continue;
    if(!CheckMinimumEnergy(ph1)) continue;

    if(!fIsMC && fIsPHOSTriggerAnalysis && !ph1->IsTrig()) continue;//take trigger bias into account.

    for(Int_t ev=0;ev<prevPHOS->GetSize();ev++){
      TClonesArray *mixPHOS = static_cast<TClonesArray*>(prevPHOS->At(ev));

      for(Int_t i2=0;i2<mixPHOS->GetEntriesFast();i2++){
        AliCaloPhoton *ph2 = (AliCaloPhoton*)mixPHOS->At(i2);
        if(!fPHOSClusterCuts->AcceptPhoton(ph2)) continue;
        if(!CheckMinimumEnergy(ph2)) continue;

        p12 = *ph1 + *ph2;
        m12 = p12.M();
        energy = ph2->Energy();

        if(fUseCoreEnergy){
          p12core = *(ph1->GetMomV2()) + *(ph2->GetMomV2());
          m12 = p12core.M();
          energy = (ph2->GetMomV2())->Energy();
        }

        FillHistogramTH2(fOutputContainer,"hMixMgg_Probe_TOF",m12,energy);
        if(ph2->IsTOFOK()) FillHistogramTH2(fOutputContainer,"hMixMgg_PassingProbe_TOF",m12,energy);

      }//end of mix

    }//end of ph2

  }//end of ph1

}
//________________________________________________________________________
void AliAnalysisTaskPHOSPi0EtaToGammaGamma::DDAPhotonPurity()
{
  //Data Driven Approach to measure photon purity.
  //identified charged particle pi/K/p are used.
  //shower shape is same in pi/K/p, but different from EM (e/gamma) shower shape.-> hadronic contamination after shower shape cut is measured.
  //photon should be neutral hit, so, clusters matching with track is not photon.-> charged particle contamination after CPV cut is measured.
  //an electrons is EM particle, but charged.-> electrons are rejected by CPV.

  const Int_t multClust = fPHOSClusterArray->GetEntriesFast();

  Double_t pT=0, cluE = 0;
  Double_t weight = 1.;
  Int_t primary = -1;

  AliVCluster *cluster = 0x0;
  AliVTrack *track = 0x0;
  Int_t charge = 0;
  Double_t trackP = 0;
  Double_t nsigmaElectron = 999;
  Double_t nsigmaPion = 999;
  Double_t nsigmaKaon = 999;
  Double_t nsigmaProton = 999;
  Bool_t pidPion = kFALSE;
  Bool_t pidKaon = kFALSE;
  Bool_t pidProton = kFALSE;
  Bool_t pidElectron = kFALSE;
  const Double_t NsigmaCut = 3;
  Bool_t isGlobal = kFALSE;
  const Double_t Rcut_CE = 240.;

  const Double_t NsigmaDisp = fPHOSClusterCuts->GetDispParameter();
  AliInfo(Form("NsigmaDisp = %2.1f sigma",NsigmaDisp));

  for(Int_t i1=0;i1<multClust;i1++){
    AliCaloPhoton *ph = (AliCaloPhoton*)fPHOSClusterArray->At(i1);
    if(!fIsMC && fIsPHOSTriggerAnalysis && !ph->IsTrig()) continue;//it is meaningless to look at non-triggered cluster in PHOS trigger analysis.
    if(!CheckMinimumEnergy(ph)) continue;

    //if(!fIsMC && !ph->IsTOFOK()) continue;

    pT = ph->Pt();
    cluE = ph->Energy();

    if(fUseCoreEnergy){
      pT = (ph->GetMomV2())->Pt();
      cluE = (ph->GetMomV2())->Energy();
    }

    weight = 1.;
    if(fIsMC){
      primary = ph->GetPrimary();
      weight = ph->GetWeight();
    }

    cluster = 0x0;
    track = 0x0;

    cluster = (AliVCluster*)ph->GetCluster();
    isGlobal = kFALSE;

    if(fESDEvent){
      Int_t trackindex = cluster->GetTrackMatchedIndex();
      if(trackindex > 0){
        track = (AliVTrack*)(fEvent->GetTrack(trackindex));
        isGlobal = fESDtrackCutsGlobal->AcceptTrack(dynamic_cast<AliESDtrack*>(track));
      }
    }//end of ESD
    else if(fAODEvent){
      if(cluster->GetNTracksMatched() > 0){
        track = dynamic_cast<AliVTrack*>(cluster->GetTrackMatched(0));
        isGlobal = dynamic_cast<AliAODTrack*>(track)->TestFilterMask(AliAODTrack::kTrkGlobalNoDCA);//standard cuts with very loose DCA cut
      }//end of track matching
    }//end of AOD

    pidPion = kFALSE;
    pidKaon = kFALSE;
    pidProton = kFALSE;
    pidElectron = kFALSE;

    if(track){
      //Note that histograms are filled as a function of cluster pT.
      charge = track->Charge();
      trackP = track->P();
      nsigmaElectron = fPIDResponse->NumberOfSigmasTPC(track,AliPID::kElectron);//(-2,3) is for real electron analysis.
      nsigmaPion     = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track,AliPID::kPion));
      nsigmaKaon     = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track,AliPID::kKaon));
      nsigmaProton   = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track,AliPID::kProton));

      if((TMath::Abs(nsigmaElectron) < nsigmaPion)     && (TMath::Abs(nsigmaElectron) < nsigmaKaon) && (TMath::Abs(nsigmaElectron) < nsigmaProton) && (-2 < nsigmaElectron && nsigmaElectron < 3)) pidElectron = kTRUE;
      if((nsigmaPion     < TMath::Abs(nsigmaElectron)) && (nsigmaPion     < nsigmaKaon) && (nsigmaPion     < nsigmaProton) && (nsigmaPion     < NsigmaCut)) pidPion = kTRUE;
      if((nsigmaKaon     < TMath::Abs(nsigmaElectron)) && (nsigmaKaon     < nsigmaPion) && (nsigmaKaon     < nsigmaProton) && (nsigmaKaon     < NsigmaCut)) pidKaon = kTRUE;
      if((nsigmaProton   < TMath::Abs(nsigmaElectron)) && (nsigmaProton   < nsigmaPion) && (nsigmaProton   < nsigmaKaon)   && (nsigmaProton   < NsigmaCut)) pidProton = kTRUE;

      if(ph->GetNsigmaCPV() < fMatchingR){//matching cut to evaluate dispersion cut efficiency for charged particle.
        if(pidElectron && (0.8 < cluE/trackP && cluE/trackP < 1.2))  FillHistogramTH1(fOutputContainer,"hMatchedElectron",pT,weight);
        else if(pidPion) FillHistogramTH1(fOutputContainer,"hMatchedPion",pT,weight);
        else if(pidKaon) FillHistogramTH1(fOutputContainer,"hMatchedKaon",pT,weight);
        else if(pidProton){
          if(charge > 0) FillHistogramTH1(fOutputContainer,"hMatchedProton",pT,weight);
          else           FillHistogramTH1(fOutputContainer,"hMatchedAntiProton",pT,weight);
        }

        if(fPHOSClusterCuts->AcceptDisp(ph)){
          if(pidElectron && (0.8 < cluE/trackP && cluE/trackP < 1.2))  FillHistogramTH1(fOutputContainer,"hMatchedElectron_Disp",pT,weight);
          else if(pidPion) FillHistogramTH1(fOutputContainer,"hMatchedPion_Disp",pT,weight);
          else if(pidKaon) FillHistogramTH1(fOutputContainer,"hMatchedKaon_Disp",pT,weight);
          else if(pidProton){
            if(charge > 0) FillHistogramTH1(fOutputContainer,"hMatchedProton_Disp",pT,weight);
            else           FillHistogramTH1(fOutputContainer,"hMatchedAntiProton_Disp",pT,weight);
          }
        }//end of disp

      }//end of track matching

    }//end of track

    if(fIsMC){
      AliAODMCParticle *p = (AliAODMCParticle*)fMCArrayAOD->At(primary);
      Int_t pdg = p->PdgCode(); 
      Double_t x = p->Xv();//absolute coordinate in ALICE
      Double_t y = p->Yv();//absolute coordinate in ALICE
      Double_t Rxy = TMath::Sqrt(x*x + y*y);

      Int_t motherid = p->GetMother();
      Int_t pdg_mother = 0;//0 is not assiend to any particle

      if(motherid > -1){
        AliAODMCParticle *mp = (AliAODMCParticle*)fMCArrayAOD->At(motherid);
        pdg_mother = mp->PdgCode();
      }

      //border of Rxy is 250 cm from (0,0,0) where TPC outer frame is.
      //only for safety mergin, 240 cm is used.

      if(pdg == 22) FillHistogramTH1(fOutputContainer,"hPurityGamma_noPID",pT,weight);

      else if(TMath::Abs(pdg) == 11){
        FillHistogramTH2(fOutputContainer,"hElectronRxy_noPID",pT,Rxy,weight);
        if(motherid > -1 && pdg_mother == 22){//conversion gamma->ee
          FillHistogramTH2(fOutputContainer,"hConvertedElectronRxy_noPID",pT,Rxy,weight);
          if(Rxy < Rcut_CE) FillHistogramTH1(fOutputContainer,"hPurityElectron_noPID",pT,weight);
          else              FillHistogramTH1(fOutputContainer,"hPurityLCE_noPID",pT,weight);
        }
        else FillHistogramTH1(fOutputContainer,"hPurityElectron_noPID",pT,weight);
      }
      else if(TMath::Abs(pdg) == 211) FillHistogramTH1(fOutputContainer,"hPurityPion_noPID",pT,weight);
      else if(TMath::Abs(pdg) == 321) FillHistogramTH1(fOutputContainer,"hPurityKaon_noPID",pT,weight);
      else if(TMath::Abs(pdg) == 130) FillHistogramTH1(fOutputContainer,"hPurityK0L_noPID",pT,weight);
      else if(pdg ==  2212)           FillHistogramTH1(fOutputContainer,"hPurityProton_noPID",pT,weight);
      else if(pdg == -2212)           FillHistogramTH1(fOutputContainer,"hPurityAntiProton_noPID",pT,weight);
      else if(pdg ==  2112)           FillHistogramTH1(fOutputContainer,"hPurityNeutron_noPID",pT,weight);
      else if(pdg == -2112)           FillHistogramTH1(fOutputContainer,"hPurityAntiNeutron_noPID",pT,weight);
      else{
        if(pdg == 111){//hadronic interaction
          //printf("mother pdg of %d is %d and production vertex = %4.3f\n",pdg,pdg_mother,Rxy);
          if(pdg_mother ==  2212)      FillHistogramTH1(fOutputContainer,"hPurityProton_noPID",pT,weight);
          else if(pdg_mother == -2212) FillHistogramTH1(fOutputContainer,"hPurityAntiProton_noPID",pT,weight);
          else if(pdg_mother ==  2112) FillHistogramTH1(fOutputContainer,"hPurityNeutron_noPID",pT,weight);
          else if(pdg_mother == -2112) FillHistogramTH1(fOutputContainer,"hPurityAntiNeutron_noPID",pT,weight);
          else                         FillHistogramTH1(fOutputContainer,"hPurityOthers_noPID",pT,weight);
        }
        else                           FillHistogramTH1(fOutputContainer,"hPurityOthers_noPID",pT,weight);
      }

      if(fPHOSClusterCuts->IsNeutral(ph)){
        if(pdg == 22) FillHistogramTH1(fOutputContainer,"hPurityGamma_CPV",pT,weight);
        else if(TMath::Abs(pdg) == 11){
          FillHistogramTH2(fOutputContainer,"hElectronRxy_CPV",pT,Rxy,weight);
          if(motherid > -1 && pdg_mother == 22){//conversion gamma->ee
            FillHistogramTH2(fOutputContainer,"hConvertedElectronRxy_CPV",pT,Rxy,weight);
            if(Rxy < Rcut_CE) FillHistogramTH1(fOutputContainer,"hPurityElectron_CPV",pT,weight);
            else              FillHistogramTH1(fOutputContainer,"hPurityLCE_CPV",pT,weight);
          }
          else FillHistogramTH1(fOutputContainer,"hPurityElectron_CPV",pT,weight);
        }
        else if(TMath::Abs(pdg) == 211) FillHistogramTH1(fOutputContainer,"hPurityPion_CPV",pT,weight);
        else if(TMath::Abs(pdg) == 321) FillHistogramTH1(fOutputContainer,"hPurityKaon_CPV",pT,weight);
        else if(TMath::Abs(pdg) == 130) FillHistogramTH1(fOutputContainer,"hPurityK0L_CPV",pT,weight);
        else if(pdg ==  2212)           FillHistogramTH1(fOutputContainer,"hPurityProton_CPV",pT,weight);
        else if(pdg == -2212)           FillHistogramTH1(fOutputContainer,"hPurityAntiProton_CPV",pT,weight);
        else if(pdg ==  2112)           FillHistogramTH1(fOutputContainer,"hPurityNeutron_CPV",pT,weight);
        else if(pdg == -2112)           FillHistogramTH1(fOutputContainer,"hPurityAntiNeutron_CPV",pT,weight);
        else{
          if(pdg == 111){//hadronic interaction
            if(pdg_mother ==  2212)      FillHistogramTH1(fOutputContainer,"hPurityProton_CPV",pT,weight);
            else if(pdg_mother == -2212) FillHistogramTH1(fOutputContainer,"hPurityAntiProton_CPV",pT,weight);
            else if(pdg_mother ==  2112) FillHistogramTH1(fOutputContainer,"hPurityNeutron_CPV",pT,weight);
            else if(pdg_mother == -2112) FillHistogramTH1(fOutputContainer,"hPurityAntiNeutron_CPV",pT,weight);
            else                         FillHistogramTH1(fOutputContainer,"hPurityOthers_CPV",pT,weight);
          }
          else                           FillHistogramTH1(fOutputContainer,"hPurityOthers_CPV",pT,weight);
        }

      }//end of CPV

      if(fPHOSClusterCuts->AcceptDisp(ph)){
        if(pdg == 22)                   FillHistogramTH1(fOutputContainer,"hPurityGamma_Disp",pT,weight);
        else if(TMath::Abs(pdg) == 11){
          FillHistogramTH2(fOutputContainer,"hElectronRxy_Disp",pT,Rxy,weight);
          if(motherid > -1 && pdg_mother == 22){//conversion gamma->ee
            FillHistogramTH2(fOutputContainer,"hConvertedElectronRxy_Disp",pT,Rxy,weight);
            if(Rxy < Rcut_CE) FillHistogramTH1(fOutputContainer,"hPurityElectron_Disp",pT,weight);
            else              FillHistogramTH1(fOutputContainer,"hPurityLCE_Disp",pT,weight);
          }
          else FillHistogramTH1(fOutputContainer,"hPurityElectron_Disp",pT,weight);
        }
        else if(TMath::Abs(pdg) == 211) FillHistogramTH1(fOutputContainer,"hPurityPion_Disp",pT,weight);
        else if(TMath::Abs(pdg) == 321) FillHistogramTH1(fOutputContainer,"hPurityKaon_Disp",pT,weight);
        else if(TMath::Abs(pdg) == 130) FillHistogramTH1(fOutputContainer,"hPurityK0L_Disp",pT,weight);
        else if(pdg ==  2212)           FillHistogramTH1(fOutputContainer,"hPurityProton_Disp",pT,weight);
        else if(pdg == -2212)           FillHistogramTH1(fOutputContainer,"hPurityAntiProton_Disp",pT,weight);
        else if(pdg ==  2112)           FillHistogramTH1(fOutputContainer,"hPurityNeutron_Disp",pT,weight);
        else if(pdg == -2112)           FillHistogramTH1(fOutputContainer,"hPurityAntiNeutron_Disp",pT,weight);
        else{
          if(pdg == 111){//hadronic interaction
            if(pdg_mother ==  2212)      FillHistogramTH1(fOutputContainer,"hPurityProton_Disp",pT,weight);
            else if(pdg_mother == -2212) FillHistogramTH1(fOutputContainer,"hPurityAntiProton_Disp",pT,weight);
            else if(pdg_mother ==  2112) FillHistogramTH1(fOutputContainer,"hPurityNeutron_Disp",pT,weight);
            else if(pdg_mother == -2112) FillHistogramTH1(fOutputContainer,"hPurityAntiNeutron_Disp",pT,weight);
            else                         FillHistogramTH1(fOutputContainer,"hPurityOthers_Disp",pT,weight);
          }
          else                           FillHistogramTH1(fOutputContainer,"hPurityOthers_Disp",pT,weight);
        }

      }//end of Disp

      if(fPHOSClusterCuts->AcceptPhoton(ph)){
        if(pdg == 22)                   FillHistogramTH1(fOutputContainer,"hPurityGamma_PID",pT,weight);
        else if(TMath::Abs(pdg) == 11){
          FillHistogramTH2(fOutputContainer,"hElectronRxy_PID",pT,Rxy,weight);
          if(motherid > -1 && pdg_mother == 22){//conversion gamma->ee
            FillHistogramTH2(fOutputContainer,"hConvertedElectronRxy_PID",pT,Rxy,weight);
            if(Rxy < Rcut_CE) FillHistogramTH1(fOutputContainer,"hPurityElectron_PID",pT,weight);
            else              FillHistogramTH1(fOutputContainer,"hPurityLCE_PID",pT,weight);
          }
          else FillHistogramTH1(fOutputContainer,"hPurityElectron_PID",pT,weight);
        }
        else if(TMath::Abs(pdg) == 211) FillHistogramTH1(fOutputContainer,"hPurityPion_PID",pT,weight);
        else if(TMath::Abs(pdg) == 321) FillHistogramTH1(fOutputContainer,"hPurityKaon_PID",pT,weight);
        else if(TMath::Abs(pdg) == 130) FillHistogramTH1(fOutputContainer,"hPurityK0L_PID",pT,weight);
        else if(pdg ==  2212)           FillHistogramTH1(fOutputContainer,"hPurityProton_PID",pT,weight);
        else if(pdg == -2212)           FillHistogramTH1(fOutputContainer,"hPurityAntiProton_PID",pT,weight);
        else if(pdg ==  2112)           FillHistogramTH1(fOutputContainer,"hPurityNeutron_PID",pT,weight);
        else if(pdg == -2112)           FillHistogramTH1(fOutputContainer,"hPurityAntiNeutron_PID",pT,weight);
        else{
          if(pdg == 111){//hadronic interaction
            if(pdg_mother ==  2212)      FillHistogramTH1(fOutputContainer,"hPurityProton_PID",pT,weight);
            else if(pdg_mother == -2212) FillHistogramTH1(fOutputContainer,"hPurityAntiProton_PID",pT,weight);
            else if(pdg_mother ==  2112) FillHistogramTH1(fOutputContainer,"hPurityNeutron_PID",pT,weight);
            else if(pdg_mother == -2112) FillHistogramTH1(fOutputContainer,"hPurityAntiNeutron_PID",pT,weight);
            else                         FillHistogramTH1(fOutputContainer,"hPurityOthers_PID",pT,weight);
          }
          else                           FillHistogramTH1(fOutputContainer,"hPurityOthers_PID",pT,weight);
        }

      }//end of PID

    }//end of M.C.

  }//end of ph1

}
//________________________________________________________________________
void AliAnalysisTaskPHOSPi0EtaToGammaGamma::SelectTriggeredCluster()
{

  Int_t trgrelId[4]={};
  Int_t tmod=0; //"Online" module number
  Int_t trgabsId=0; //bottom-left 4x4 edge cell absId [1-64], [1-56]
  Int_t trgmodule=0,trgcellx=0,trgcellz=0;

  Int_t relId[4]={};
  Int_t module=0,cellx=0,cellz=0,tru=0;
  Double_t energy=0;
  Float_t position[3] = {}; 
 
  //AliVCaloCells *cells = dynamic_cast<AliVCaloCells*>(fEvent->GetPHOSCells());
  const Int_t L1input = fPHOSTriggerHelper->GetL1TriggerInput();
  const Int_t L0input = fPHOSTriggerHelper->GetL0TriggerInput();

  Int_t L1=-999;
  if(L1input > 0){
    if     (L1input == 7) L1 = 0;//L1 high
    else if(L1input == 6) L1 = 1;//L1 medium
    else if(L1input == 5) L1 = 2;//L1 low
  }
  else if(L0input > 0){
    //if(L0input == 9) L1 = -1;//L0
    L1 = -1;//L0
  }

  const Int_t multClust = fPHOSClusterArray->GetEntriesFast();
  //set label to all clusters. PID/TOF cut is evaluated EstimateTOFCutEfficiency, EstimatePIDCutEfficiency.
  for(Int_t i=0;i<multClust;i++){
    AliCaloPhoton *ph = (AliCaloPhoton*)fPHOSClusterArray->At(i);
    //AliVCluster *clu1 = (AliVCluster*)ph->GetCluster();

    position[0] = ph->EMCx();
    position[1] = ph->EMCy();
    position[2] = ph->EMCz();
    TVector3 global1(position);
    fPHOSGeo->GlobalPos2RelId(global1,relId);

    //Int_t maxAbsId = fPHOSTriggerHelper->FindHighestAmplitudeCellAbsId(clu1, cells);
    //fPHOSGeo->AbsToRelNumbering(maxAbsId,relId);

    module = relId[0];
    cellx  = relId[2];
    cellz  = relId[3];

    if(module < 1 || module > 4){
      AliError(Form("Wrong module number %d",module));
      return;
    }
    AliInfo(Form("cluster position M:%d, X:%d , Z:%d",module,cellx,cellz));

    AliVCaloTrigger* trg = fEvent->GetCaloTrigger("PHOS");
    trg->Reset();

    while(trg->Next()){

      if(trg->GetL1TimeSum() != L1) continue;

      trg->GetPosition(tmod,trgabsId);//tmod is online module numbering. i.e., tmod=1 means half module.

      fPHOSGeo->AbsToRelNumbering(trgabsId,trgrelId);
      //for offline numbering, relId should be used.

      trgmodule = trgrelId[0];
      trgcellx  = trgrelId[2];
      trgcellz  = trgrelId[3];

      AliInfo(Form("fired TRU channel M:%d, X:%d , Z:%d",trgmodule,trgcellx,trgcellz));

      if(trgmodule < 1 || trgmodule > 4){
        AliError(Form("Wrong module number %d",trgmodule));
        return;
      }

      if(fPHOSTriggerHelper->IsDeltaRUsed()){//decision is taken by DeltaR = sqrt(deta^2 + dphi^2)
        if(fPHOSTriggerHelper->IsMatchedDeltaR(trgrelId,global1)){
          ph->SetTrig(kTRUE);
          break;
        }
      }
      else{//decision is taken by dx and dz
        if(fPHOSTriggerHelper->IsMatched(trgrelId,relId)){
          ph->SetTrig(kTRUE);
          break;
        }
      }
    }//end of while trigger patch loop

  }//end of cluster loop

  relId[0] = 0; relId[1] = 0; relId[2] = 0; relId[3] = 0;
  module = 0; cellx = 0; cellz = 0; tru = 0;

  for(Int_t i=0;i<multClust;i++){
    AliCaloPhoton *ph = (AliCaloPhoton*)fPHOSClusterArray->At(i);
    if(!fPHOSClusterCuts->AcceptPhoton(ph)) continue;
    if(!CheckMinimumEnergy(ph)) continue;
    if(!fPHOSTriggerHelper->IsOnActiveTRUChannel(ph)) continue;

    energy = ph->Energy();
    if(fUseCoreEnergy) energy = (ph->GetMomV2())->Energy();

    position[0] = ph->EMCx();
    position[1] = ph->EMCy();
    position[2] = ph->EMCz();
    TVector3 global1(position);
    fPHOSGeo->GlobalPos2RelId(global1,relId);

    //AliVCluster *clu1 = (AliVCluster*)ph->GetCluster();
    //Int_t maxAbsId = fPHOSTriggerHelper->FindHighestAmplitudeCellAbsId(clu1, cells);
    //fPHOSGeo->AbsToRelNumbering(maxAbsId,relId);

    module = relId[0];
    cellx  = relId[2];
    cellz  = relId[3];
    tru = fPHOSTriggerHelper->WhichTRU(cellx,cellz);

    FillHistogramTH1(fOutputContainer,Form("hClusterE_M%d_TRU%d",module,tru),energy);
    if(ph->IsTrig()) FillHistogramTH1(fOutputContainer,Form("hTriggeredClusterE_M%d_TRU%d",module,tru),energy);
    //for the QA purpose, TOF cut efficiency is not corrected.
    if(ph->IsTOFOK()){
      FillHistogramTH1(fOutputContainer,Form("hClusterE_TOF_M%d_TRU%d",module,tru),energy);
      if(ph->IsTrig()) FillHistogramTH1(fOutputContainer,Form("hTriggeredClusterE_TOF_M%d_TRU%d",module,tru),energy);
    }
  }//end of cluster loop


}
//________________________________________________________________________
void AliAnalysisTaskPHOSPi0EtaToGammaGamma::EstimateTriggerEfficiency()
{
  //tag and probe method is used.

  //AliVCaloCells *cells = dynamic_cast<AliVCaloCells*>(fEvent->GetPHOSCells());
  const Int_t multClust = fPHOSClusterArray->GetEntriesFast();

  Double_t position[3] = {};
  Double_t energy=0;

  Int_t module_tag=0,cellx_tag=0,cellz_tag=0,tru_tag=0;
  Int_t module=0,cellx=0,cellz=0,tru=0;
  Int_t relId[4]={};
  Int_t truch_tag = -1, chX_tag=-1, chZ_tag=-1;
  Int_t truch = -1, chX=-1, chZ=-1;

  //Double_t eta_tag = -999;
  //Double_t phi_tag = -999;
  //Double_t DeltaR_tag = -999;
  //Double_t eta_probe = -999;
  //Double_t phi_probe = -999;
  //Double_t DeltaR_probe = -999;
  //Double_t DeltaRgg = -999;
  //Double_t Rlimit = 0;

  TLorentzVector p12, p12core;
  Double_t m12=0;
  //Int_t maxAbsId = -1;

  for(Int_t i1=0;i1<multClust;i1++){
    AliCaloPhoton *ph1 = (AliCaloPhoton*)fPHOSClusterArray->At(i1);
    if(!fPHOSClusterCuts->AcceptPhoton(ph1)) continue;
    if(!CheckMinimumEnergy(ph1)) continue;

    if(!ph1->IsTrig()) continue;
    if(ph1->Energy() < fEnergyThreshold) continue;

    relId[0] = 0; relId[1] = 0; relId[2] = 0; relId[3] = 0;

    position[0] = ph1->EMCx();
    position[1] = ph1->EMCy();
    position[2] = ph1->EMCz();
    TVector3 global_tag(position);
    fPHOSGeo->GlobalPos2RelId(global_tag,relId);

    //eta_tag = global_tag.Eta();
    //phi_tag = global_tag.Phi();
    //if(phi_tag < 0) phi_tag += TMath::TwoPi();

    //maxAbsId = fPHOSTriggerHelper->FindHighestAmplitudeCellAbsId(ph1->GetCluster(), cells);
    //fPHOSGeo->AbsToRelNumbering(maxAbsId,relId);

    module_tag = relId[0];
    cellx_tag  = relId[2];
    cellz_tag  = relId[3];
    tru_tag = fPHOSTriggerHelper->WhichTRU(cellx_tag,cellz_tag);
    truch_tag = fPHOSTriggerHelper->WhichTRUChannel(cellx_tag,cellz_tag,chX_tag,chZ_tag);

    //DeltaR_tag = fPHOSTriggerHelper->GetDistanceToClosestTRUChannel(ph1);

    for(Int_t i2=0;i2<multClust;i2++){
      AliCaloPhoton *ph2 = (AliCaloPhoton*)fPHOSClusterArray->At(i2);
      if(!fPHOSClusterCuts->AcceptPhoton(ph2)) continue;
      if(!CheckMinimumEnergy(ph2)) continue;

      if(!fPHOSTriggerHelper->IsOnActiveTRUChannel(ph2)) continue;

      if(i2==i1) continue;//reject same cluster combination

      p12 = *ph1 + *ph2;
      m12 = p12.M();
      energy = ph2->Energy();

      if(fUseCoreEnergy){
        p12core = *(ph1->GetMomV2()) + *(ph2->GetMomV2());
        m12 = p12core.M();
        energy = (ph2->GetMomV2())->Energy();
      }

      relId[0] = 0; relId[1] = 0; relId[2] = 0; relId[3] = 0;

      position[0] = ph2->EMCx();
      position[1] = ph2->EMCy();
      position[2] = ph2->EMCz();
      TVector3 global_probe(position);
      fPHOSGeo->GlobalPos2RelId(global_probe,relId);
      //eta_probe = global_probe.Eta();
      //phi_probe = global_probe.Phi();
      //if(phi_probe < 0) phi_probe += TMath::TwoPi();

      //maxAbsId = fPHOSTriggerHelper->FindHighestAmplitudeCellAbsId(ph2->GetCluster(), cells);
      //fPHOSGeo->AbsToRelNumbering(maxAbsId,relId);

      module = relId[0];
      cellx  = relId[2];
      cellz  = relId[3];
      tru = fPHOSTriggerHelper->WhichTRU(cellx,cellz);
      truch = fPHOSTriggerHelper->WhichTRUChannel(cellx,cellz,chX,chZ);

      Int_t dm = module_tag - module;
      Int_t dt = tru_tag - tru;
      Int_t dx = chX_tag - chX;
      Int_t dz = chZ_tag - chZ;

      if( (dm == 0)
       && (dt == 0)
       && (TMath::Abs(dx) < 4 && TMath::Abs(dz) < 4)
       ) continue;//reject cluster pair where they belong to same 4x4 region.//also considuer leakage of EM shower at high energy

      AliInfo(Form("dm = %d , dt = %d , dx = %d , dz = %d, truch = %d , truch_tag = %d.",dm,dt,dx,dz,truch,truch_tag));
      //DeltaR_probe = fPHOSTriggerHelper->GetDistanceToClosestTRUChannel(ph2);

      //Rlimit = DeltaR_tag + DeltaR_probe;
      //DeltaRgg = TMath::Sqrt(TMath::Power(eta_tag - eta_probe,2) + TMath::Power(phi_tag - phi_probe,2));
      //AliInfo(Form("DeltaR between 2 gammas = %e, DeltaR_tag = %e , DeltaR_probe = = %e.",DeltaRgg,DeltaR_tag,DeltaR_probe));
      //if(DeltaRgg < Rlimit) continue;//efficiency can be measured at low pT by pi0 and at high pT by eta meson.

      FillHistogramTH2(fOutputContainer,"hMgg_Probe_Trg",m12,energy);
      FillHistogramTH2(fOutputContainer,Form("hMgg_Probe_Trg_M%d_TRU%d",module,tru),m12,energy);

      if(ph2->IsTrig()){
        FillHistogramTH2(fOutputContainer,"hMgg_PassingProbe_Trg",m12,energy);
        FillHistogramTH2(fOutputContainer,Form("hMgg_PassingProbe_Trg_M%d_TRU%d",module,tru),m12,energy);
      }

    }//end of ph2

  }//end of ph1

  //next mixed event
  TList *prevPHOS = fPHOSEvents[fZvtx][fEPBin];

  for(Int_t i1=0;i1<multClust;i1++){
    AliCaloPhoton *ph1 = (AliCaloPhoton*)fPHOSClusterArray->At(i1);
    if(!fPHOSClusterCuts->AcceptPhoton(ph1)) continue;
    if(!CheckMinimumEnergy(ph1)) continue;

    if(!ph1->IsTrig()) continue;
    if(ph1->Energy() < fEnergyThreshold) continue;

    //position[0] = ph1->EMCx();
    //position[1] = ph1->EMCy();
    //position[2] = ph1->EMCz();
    //TVector3 global_tag(position);
    //fPHOSGeo->GlobalPos2RelId(global_tag,relId);

    //eta_tag = global_tag.Eta();
    //phi_tag = global_tag.Phi();
    //if(phi_tag < 0) phi_tag += TMath::TwoPi();

    for(Int_t ev=0;ev<prevPHOS->GetSize();ev++){
      TClonesArray *mixPHOS = static_cast<TClonesArray*>(prevPHOS->At(ev));

      for(Int_t i2=0;i2<mixPHOS->GetEntriesFast();i2++){
        AliCaloPhoton *ph2 = (AliCaloPhoton*)mixPHOS->At(i2);
        if(!fPHOSClusterCuts->AcceptPhoton(ph2)) continue;
        if(!CheckMinimumEnergy(ph2)) continue;

        if(!fPHOSTriggerHelper->IsOnActiveTRUChannel(ph2)) continue;

        p12 = *ph1 + *ph2;
        m12 = p12.M();
        energy = ph2->Energy();

        if(fUseCoreEnergy){
          p12core = *(ph1->GetMomV2()) + *(ph2->GetMomV2());
          m12 = p12core.M();
          energy = (ph2->GetMomV2())->Energy();
        }

        position[0] = ph2->EMCx();
        position[1] = ph2->EMCy();
        position[2] = ph2->EMCz();

        relId[0] = 0; relId[1] = 0; relId[2] = 0; relId[3] = 0;
        TVector3 global_probe(position);
        fPHOSGeo->GlobalPos2RelId(global_probe,relId);
        //eta_probe = global_probe.Eta();
        //phi_probe = global_probe.Phi();
        //if(phi_probe < 0) phi_probe += TMath::TwoPi();

        module = relId[0];
        cellx  = relId[2];
        cellz  = relId[3];
        tru = fPHOSTriggerHelper->WhichTRU(cellx,cellz);

        //DeltaRgg = TMath::Sqrt(TMath::Power(eta_tag - eta_probe,2) + TMath::Power(phi_tag - phi_probe,2));
        //if(DeltaRgg < Rlimit) continue;

        FillHistogramTH2(fOutputContainer,"hMixMgg_Probe_Trg",m12,energy);
        FillHistogramTH2(fOutputContainer,Form("hMixMgg_Probe_Trg_M%d_TRU%d",module,tru),m12,energy);
        if(ph2->IsTrig()){
          FillHistogramTH2(fOutputContainer,"hMixMgg_PassingProbe_Trg",m12,energy);
          FillHistogramTH2(fOutputContainer,Form("hMixMgg_PassingProbe_Trg_M%d_TRU%d",module,tru),m12,energy);
        }
      }//end of mix

    }//end of ph2

  }//end of ph1

}
//________________________________________________________________________
void AliAnalysisTaskPHOSPi0EtaToGammaGamma::ProcessMC()
{
  //This is for analyzing general purpose MC such as pure PYTHIA, HIJING, DPMJET, PHOJET and so on.
  //get MC information

  Int_t firstJetindex = -1;
  Int_t lastJetindex  = -1;
  Int_t genIDJet      = -1;
  Int_t firstUEindex  = -1;
  Int_t lastUEindex   = -1;
  Int_t genIDUE       = -1;

  if(fIsJJMC){
    AliMCEvent *mcevent = MCEvent();
    AliGenPythiaEventHeader* pythiaGenHeader = 0x0;
    if(fESDEvent) pythiaGenHeader      = fJJMCHandler->GetPythiaEventHeader(mcevent);
    else if(fAODEvent) pythiaGenHeader = fJJMCHandler->GetPythiaEventHeader(fAODEvent);

    Float_t pThard   = pythiaGenHeader->GetPtHard();
    FillHistogramTH1(fOutputContainer,"hPtHard",pThard);

    firstJetindex = fJJMCHandler->GetFirstJetIndex();
    lastJetindex  = fJJMCHandler->GetLastJetIndex();
    genIDJet      = fJJMCHandler->GetGeneratorJetIndex();

    firstUEindex  = fJJMCHandler->GetFirstUEIndex();
    lastUEindex   = fJJMCHandler->GetLastUEIndex();
    genIDUE       = fJJMCHandler->GetGeneratorUEIndex();

    AliInfo(Form("genIDUE = %d , genIDJet = %d , firstindexJet = %d , lastindexJet = %d.",genIDUE,genIDJet,firstJetindex,lastJetindex));
  }

  TF1 *f1Pi0Weight   = (TF1*)GetAdditionalPi0PtWeightFunction(fCentralityMain);
  TF1 *f1EtaWeight   = (TF1*)GetAdditionalEtaPtWeightFunction(fCentralityMain);
  //TF1 *f1GammaWeight = (TF1*)GetAdditionalGammaPtWeightFunction(fCentralityMain);
  TF1 *f1K0SWeight   = (TF1*)GetAdditionalK0SPtWeightFunction(fCentralityMain);
  TF1 *f1L0Weight    = (TF1*)GetAdditionalL0PtWeightFunction(fCentralityMain);

  Double_t TruePi0Pt = 0.;
  Double_t TrueK0SPt = 0.;
  Double_t TrueL0Pt  = 0.;
  Double_t TrueEtaPt = 0.;

  Int_t genID = -1;
  Double_t pT=0, rapidity=0, phi=0;
  Double_t weight = 1;
  Int_t pdg = 0;
  TString parname = "";
  TString genname = "";
  //Int_t motherid = -1;
  //AliAODMCParticle *mp = 0x0;//mother particle
  //Double_t motherpT = 0;

  if(fESDEvent){//for ESD
    fMCArrayESD = (AliStack*)GetMCInfoESD();
    if(!fMCArrayESD){
      AliError("Could not get MC Stack!");
      return;
    }

    //const Int_t Ntrack = fMCArrayESD->GetNtrack();//this is the number of all particles (event geneartor + GEANT).
    const Int_t Ntrack = fMCArrayESD->GetNprimary();//this is the number of generated particles by event generator.
    for(Int_t i=0;i<Ntrack;i++){
      TParticle *p = (TParticle*)fMCArrayESD->Particle(i);

      if(fIsJJMC){
        Int_t primary = FindPrimaryMotherESD(i);
        if(fMCType.Contains("JJMC") && (primary < firstJetindex || lastJetindex < primary)) continue;
        if(fMCType.Contains("MBMC") && (primary < firstUEindex  || lastUEindex  < primary)) continue;
      }

      pT = p->Pt();
      rapidity = p->Y();
      phi = p->Phi();
      pdg = p->GetPdgCode();

//      if(fMCArrayESD->IsPhysicalPrimary(i))
//      FillHistogramTH1(fOutputContainer,"hPDGPhysicalPrimary",pdg);

      //rapidity is Y(), but, pseudo-rapidity is Eta();

      weight = 1.;
      if(pT < 1e-3) continue;//reject below 1 MeV
      if(TMath::Abs(rapidity) > 0.5) continue;

      Double32_t x = p->Vx() - fVertex[0];
      Double32_t y = p->Vy() - fVertex[1];
      //Double32_t z = p->Vz() - fVertex[2];
      //Double32_t Rho = sqrt(x*x + y*y + z*z);
      Double32_t R = sqrt(x*x + y*y);

      //if(Rho > 1.0) continue;//3D
      if(R > 1.0) continue;//2D
      weight = 1.;

      if(pdg==111){//pi0
        parname = "Pi0";
        weight = f1Pi0Weight->Eval(pT);
        if(IsFrom(i,TrueK0SPt,310)) weight = f1K0SWeight->Eval(TrueK0SPt) * f1Pi0Weight->Eval(TrueK0SPt);
        if(IsFrom(i,TrueEtaPt,221)) weight = f1EtaWeight->Eval(TrueEtaPt) * f1Pi0Weight->Eval(TrueEtaPt);
        if(IsFrom(i,TrueL0Pt,3122)) weight = f1L0Weight->Eval(TrueL0Pt)   * f1Pi0Weight->Eval(TrueL0Pt);

      }
      else if(pdg==221){//eta
        parname = "Eta";
        weight = f1EtaWeight->Eval(pT) * f1Pi0Weight->Eval(pT);

      }
      else if(pdg==22){//gamma
        parname = "Gamma";
        if(IsFrom(i,TruePi0Pt,111)) weight = f1Pi0Weight->Eval(TruePi0Pt);
        if(IsFrom(i,TrueK0SPt,310)) weight = f1K0SWeight->Eval(TrueK0SPt) * f1Pi0Weight->Eval(TrueK0SPt);
        if(IsFrom(i,TrueEtaPt,221)) weight = f1EtaWeight->Eval(TrueEtaPt) * f1Pi0Weight->Eval(TrueEtaPt);
        if(IsFrom(i,TrueL0Pt,3122)) weight = f1L0Weight->Eval(TrueL0Pt)   * f1Pi0Weight->Eval(TrueL0Pt);
      }
  
      else if(pdg==223){//omega 782 meson
        parname = "Omega";
        weight = 1.;
      }

      else if(pdg==211 || pdg==-211){//pi+ or pi-
        //c x tau = 7.8m
        parname = "ChargedPion";
        weight = f1Pi0Weight->Eval(pT);
        if(IsFrom(i,TrueK0SPt,310)) weight = f1K0SWeight->Eval(TrueK0SPt) * f1Pi0Weight->Eval(TrueK0SPt);
        if(IsFrom(i,TrueEtaPt,221)) weight = f1EtaWeight->Eval(TrueEtaPt) * f1Pi0Weight->Eval(TrueEtaPt);
        if(IsFrom(i,TrueL0Pt,3122)) weight = f1L0Weight->Eval(TrueL0Pt)   * f1Pi0Weight->Eval(TrueL0Pt);
      }
      else if(pdg==321 || pdg==-321){//K+ or K-
        //c x tau = 3.7m
        parname = "ChargedKaon";
        weight = f1K0SWeight->Eval(pT) * f1Pi0Weight->Eval(pT);
      }
      else if(pdg==310){//K0S
        //c x tau = 2.7cm
        parname = "K0S";
        weight = f1K0SWeight->Eval(pT) * f1Pi0Weight->Eval(pT);
      }
      else if(pdg==130){//K0L
        //c x tau = 15.34m
        parname = "K0L";
        weight = f1K0SWeight->Eval(pT) * f1Pi0Weight->Eval(pT);
      }
      else if(TMath::Abs(pdg) == 3122){//Lmabda0
        //c x tau = 7.89cm
        parname = "Lambda0";
        weight = f1L0Weight->Eval(pT) * f1Pi0Weight->Eval(pT);
      }
      else if(TMath::Abs(pdg) == 3212){//Sigma0
        parname = "Sigma0";
        weight = 1.;
      }
      else if(pdg == 2212){//proton
        parname = "Proton";
        weight = 1.;
      }
      else if(pdg == -2212){//anti-proton
        parname = "AntiProton";
        weight = 1.;
      }
      else if(pdg == 2112){//neutron
        parname = "Neutron";
        weight = 1.;
      }
      else if(pdg == -2112){//anti-neutron
        parname = "AntiNeutron";
        weight = 1.;
      }
      else{
        continue;
      }

      FillHistogramTH1(fOutputContainer,Form("hGen%sPt"    ,parname.Data()),pT          ,weight);
      FillHistogramTH2(fOutputContainer,Form("hGen%sEtaPhi",parname.Data()),phi,rapidity,weight);
      FillHistogramTH2(fOutputContainer,Form("hGen%sEtaPt" ,parname.Data()),rapidity,pT ,weight);

    }//end of generated particle loop

  }//end of ESD
  else if(fAODEvent){//for AOD
    fMCArrayAOD = (TClonesArray*)GetMCInfoAOD();
    if(!fMCArrayAOD){
      AliError("Could not retrieve AOD event!");
      return;
    }

    const Int_t Ntrack = fMCArrayAOD->GetEntriesFast();
    for(Int_t i=0;i<Ntrack;i++){
      AliAODMCParticle *p = (AliAODMCParticle*)fMCArrayAOD->At(i);
      genID = p->GetGeneratorIndex();

      if(fIsJJMC){
        if(fMCType.Contains("JJMC") && genID != genIDJet) continue;
        if(fMCType.Contains("MBMC") && genID != genIDUE ) continue;
      }

      //if(fIsJJMC && (i < firstJetindex || lastJetindex < i) ) continue;

      pT = p->Pt();
      rapidity = p->Y();
      phi = p->Phi();
      pdg = p->PdgCode();

      //if(p->IsPhysicalPrimary())
      //  FillHistogramTH1(fOutputContainer,"hPDGPhysicalPrimary",pdg);

      //rapidity is Y(), but, pseudo-rapidity is Eta();

      weight = 1.;
      if(pT < 1e-3) continue;//reject below 1 MeV
      if(TMath::Abs(rapidity) > 0.5) continue;

      //if(Rho(p) > 1.0) continue;//3D
      if(R(p) > 1.0) continue;//2D
      weight = 1.;

      if(pdg==111){//pi0
        parname = "Pi0";
        weight = f1Pi0Weight->Eval(pT);
        if(IsFrom(i,TrueK0SPt,310)) weight = f1K0SWeight->Eval(TrueK0SPt) * f1Pi0Weight->Eval(TrueK0SPt);
        if(IsFrom(i,TrueEtaPt,221)) weight = f1EtaWeight->Eval(TrueEtaPt) * f1Pi0Weight->Eval(TrueEtaPt);
        if(IsFrom(i,TrueL0Pt,3122)) weight = f1L0Weight->Eval(TrueL0Pt)   * f1Pi0Weight->Eval(TrueL0Pt);

        if(Are2GammasInPHOSAcceptance(i)){
          FillHistogramTH1(fOutputContainer,Form("hGen%sPtACC"    ,parname.Data()),pT          ,weight);
          FillHistogramTH2(fOutputContainer,Form("hGen%sEtaPhiACC",parname.Data()),phi,rapidity,weight);
          FillHistogramTH2(fOutputContainer,Form("hGen%sEtaPtACC" ,parname.Data()),rapidity,pT ,weight);
        }

      }
      else if(pdg==221){//eta
        parname = "Eta";
        weight = f1EtaWeight->Eval(pT) * f1Pi0Weight->Eval(pT);
        if(Are2GammasInPHOSAcceptance(i)){
          FillHistogramTH1(fOutputContainer,Form("hGen%sPtACC"    ,parname.Data()),pT          ,weight);
          FillHistogramTH2(fOutputContainer,Form("hGen%sEtaPhiACC",parname.Data()),phi,rapidity,weight);
          FillHistogramTH2(fOutputContainer,Form("hGen%sEtaPtACC" ,parname.Data()),rapidity,pT ,weight);
        }
      }
      else if(pdg==22){//gamma
        parname = "Gamma";
        if(IsFrom(i,TruePi0Pt,111)) weight = f1Pi0Weight->Eval(TruePi0Pt);
        if(IsFrom(i,TrueK0SPt,310)) weight = f1K0SWeight->Eval(TrueK0SPt) * f1Pi0Weight->Eval(TrueK0SPt);
        if(IsFrom(i,TrueEtaPt,221)) weight = f1EtaWeight->Eval(TrueEtaPt) * f1Pi0Weight->Eval(TrueEtaPt);
        if(IsFrom(i,TrueL0Pt,3122)) weight = f1L0Weight->Eval(TrueL0Pt)   * f1Pi0Weight->Eval(TrueL0Pt);
      }
      else if(pdg==223){//omega 782 meson
        parname = "Omega";
        weight = 1.;
      }
      else if(pdg==211 || pdg==-211){//pi+ or pi-
        parname = "ChargedPion";
        weight = f1Pi0Weight->Eval(pT);
        if(IsFrom(i,TrueK0SPt,310)) weight = f1K0SWeight->Eval(TrueK0SPt) * f1Pi0Weight->Eval(TrueK0SPt);
        if(IsFrom(i,TrueEtaPt,221)) weight = f1EtaWeight->Eval(TrueEtaPt) * f1Pi0Weight->Eval(TrueEtaPt);
        if(IsFrom(i,TrueL0Pt,3122)) weight = f1L0Weight->Eval(TrueL0Pt)   * f1Pi0Weight->Eval(TrueL0Pt);
      }
      else if(pdg==321 || pdg==-321){//K+ or K-
        parname = "ChargedKaon";
        weight = f1K0SWeight->Eval(pT) * f1Pi0Weight->Eval(pT);
      }
      else if(pdg==310){//K0S
        parname = "K0S";
        weight = f1K0SWeight->Eval(pT) * f1Pi0Weight->Eval(pT);
      }
      else if(pdg==130){//K0L
        parname = "K0L";
        weight = f1K0SWeight->Eval(pT) * f1Pi0Weight->Eval(pT);
      }
      else if(TMath::Abs(pdg) == 3122){//Lmabda0
        parname = "Lambda0";
        weight = f1L0Weight->Eval(pT) * f1Pi0Weight->Eval(pT);
      }
      else if(TMath::Abs(pdg) == 3212){//Sigma0
        parname = "Sigma0";
      }
      else if(pdg == 2212){//proton
        parname = "Proton";
        weight = 1.;
      }
      else if(pdg == -2212){//anti-proton
        parname = "AntiProton";
        weight = 1.;
      }
      else if(pdg == 2112){//neutron
        parname = "Neutron";
        weight = 1.;
      }
      else if(pdg == -2112){//anti-neutron
        parname = "AntiNeutron";
        weight = 1.;
      }
      else{
        continue;
      }

      FillHistogramTH1(fOutputContainer,Form("hGen%sPt"    ,parname.Data()),pT          ,weight);
      FillHistogramTH2(fOutputContainer,Form("hGen%sEtaPhi",parname.Data()),phi,rapidity,weight);
      FillHistogramTH2(fOutputContainer,Form("hGen%sEtaPt" ,parname.Data()),rapidity,pT ,weight);

    }//end of generated particle loop
   
    //AliGenHijingEventHeader* hijingGenHeader = 0x0;
    //if(fAODEvent) {
    //  AliAODMCHeader* mcHeader = (AliAODMCHeader*) fAODEvent->GetList()->FindObject(AliAODMCHeader::StdBranchName());
    //  TList* headerList = mcHeader->GetCocktailHeaders();
    //  for (Int_t i=0; i<headerList->GetEntries(); i++) {
    //    hijingGenHeader = dynamic_cast<AliGenHijingEventHeader*>(headerList->At(i));
    //    if(hijingGenHeader) break;
    //  }
    //}
    //if(hijingGenHeader){
    //  const Int_t Nprimary = hijingGenHeader->NProduced();
    //  for(Int_t i=0;i<Nprimary;i++){
    //  AliAODMCParticle *p = (AliAODMCParticle*)fMCArrayAOD->At(i);
    //  genID = p->GetGeneratorIndex();

    //  pdg = p->PdgCode();
    //  FillHistogramTH1(fOutputContainer,"hPDGPhysicalPrimary",pdg);
    //  if(p->MCStatusCode() == 1) FillHistogramTH1(fOutputContainer,"hPDGPhysicalPrimaryStable",pdg);

    //  }

    //}
  }


}
//________________________________________________________________________
Double_t AliAnalysisTaskPHOSPi0EtaToGammaGamma::R(AliAODMCParticle* p)
{
  //Radius of vertex in cylindrical system.
  //relative position to vertex is NOT mandatory, if you want to use cylindrical system.
  //just for consistency, relative position is used here.

  Double32_t x = p->Xv() - fVertex[0];
  Double32_t y = p->Yv() - fVertex[1];
  return sqrt(x*x + y*y);
}
//________________________________________________________________________
Double_t AliAnalysisTaskPHOSPi0EtaToGammaGamma::Rho(AliAODMCParticle* p)
{
  //Radius of vertex in spherical system.
  //relative position to vertex is necessary, if you want to use spherical system.

  Double32_t x = p->Xv() - fVertex[0];
  Double32_t y = p->Yv() - fVertex[1];
  Double32_t z = p->Zv() - fVertex[2];
  return sqrt(x*x + y*y + z*z);
}
//________________________________________________________________________
Double_t AliAnalysisTaskPHOSPi0EtaToGammaGamma::RAbs(AliAODMCParticle* p)
{
  //Radius of vertex in cylindrical system.

  Double32_t x = p->Xv();
  Double32_t y = p->Yv();
  return sqrt(x*x + y*y);
}
//________________________________________________________________________
Double_t AliAnalysisTaskPHOSPi0EtaToGammaGamma::RhoAbs(AliAODMCParticle* p)
{
  //Radius of vertex in spherical system.

  Double32_t x = p->Xv();
  Double32_t y = p->Yv();
  Double32_t z = p->Zv();
  return sqrt(x*x + y*y + z*z);
}
//________________________________________________________________________
Double_t AliAnalysisTaskPHOSPi0EtaToGammaGamma::DeltaPhiIn0Pi(Double_t dphi)
{
  //this returns dphi in 0-pi range in unit of radian.
  Double_t tmp = dphi;
  while(tmp < 0)                                     tmp += 2./(Double_t)fHarmonics * TMath::Pi();
  while(tmp > 2./(Double_t)fHarmonics * TMath::Pi()) tmp -= 2./(Double_t)fHarmonics * TMath::Pi();

  return tmp;
}
//________________________________________________________________________
Bool_t AliAnalysisTaskPHOSPi0EtaToGammaGamma::IsFrom(Int_t label, Double_t &TruePt, const Int_t target_pdg)
{
  AliAODEvent *aod = dynamic_cast<AliAODEvent*>(fEvent);
  AliESDEvent *esd = dynamic_cast<AliESDEvent*>(fEvent);

  Int_t motherid = -1;
  Int_t pdg=0;
  Double_t pT = 0;

  if(esd){
    const Int_t Nprimary = fMCArrayESD->GetNprimary();//this number contains only generated particles by event generator.
    TParticle *p = (TParticle*)fMCArrayESD->Particle(label);
    TParticle *mp = 0x0;
    motherid = p->GetFirstMother();

    Double32_t x = 999;
    Double32_t y = 999;
    //Double32_t z = 999;
    //Double32_t Rho = 999;
    Double32_t R = 999;

    //printf("Nprimary = %d\n",Nprimary);

    while(motherid >= Nprimary){
      mp = (TParticle*)fMCArrayESD->Particle(motherid);
      pT = mp->Pt();
      pdg = mp->GetPdgCode(); 

      //printf("IsFrom::motherid = %d\n",motherid);

      x = mp->Vx() - fVertex[0];
      y = mp->Vy() - fVertex[1];
      //z = mp->Vz() - fVertex[2];
      //Rho = sqrt(x*x + y*y + z*z);
      R = sqrt(x*x + y*y);

      //if(TMath::Abs(pdg) == target_pdg && Rho < 1.0){//pi0 from primary vertex
      if(TMath::Abs(pdg) == target_pdg && R < 1.0){//pi0 from primary vertex
        TruePt = pT;
        return kTRUE;
      }

      motherid = mp->GetFirstMother();
    }

  }
  else if(aod){
    AliAODMCParticle *p1 = (AliAODMCParticle*)fMCArrayAOD->At(label);
    AliAODMCParticle *mp = 0x0;
    motherid = p1->GetMother();

    while(motherid > -1){
      mp = (AliAODMCParticle*)fMCArrayAOD->At(motherid);
      pT = mp->Pt();
      pdg = mp->PdgCode(); 

      //if(TMath::Abs(pdg) == target_pdg && Rho(mp) < 1.0){//pi0 from primary vertex
      if(TMath::Abs(pdg) == target_pdg && R(mp) < 1.0){//pi0 from primary vertex
        TruePt = pT;
        return kTRUE;
      }
      motherid = mp->GetMother();
    }
  }
  return kFALSE;

}
//________________________________________________________________________
Bool_t AliAnalysisTaskPHOSPi0EtaToGammaGamma::IsPhoton(Int_t label)
{
  Int_t pdg = 0;
  if(fESDEvent){
    TParticle *p = (TParticle*)fMCArrayESD->Particle(label);
    pdg = p->GetPdgCode(); 
  }
  else if(fAODEvent){
    AliAODMCParticle *p = (AliAODMCParticle*)fMCArrayAOD->At(label);
    pdg = p->PdgCode();
  }

  if(pdg == 22) return kTRUE;
  else          return kFALSE; 

}
//________________________________________________________________________
Int_t AliAnalysisTaskPHOSPi0EtaToGammaGamma::FindCommonParent(Int_t iPart, Int_t jPart)
{
  //check if there is a common parent for particles i and j
  // -1: no common parent or wrong iPart/jPart

  Int_t ntrack = fMCArrayAOD->GetEntriesFast();
  if(iPart==-1 || iPart>=ntrack || jPart==-1 || jPart>=ntrack) return -1;

  Int_t iprim1 = iPart;

  while(iprim1>-1){
    AliAODMCParticle *p1 = dynamic_cast<AliAODMCParticle*>(fMCArrayAOD->At(iprim1));
    Int_t pdg1 = p1->GetPdgCode();
    if(TMath::Abs(pdg1) <= 6 || TMath::Abs(pdg1) == 21){//reject quarks and gluons as parents.
      iprim1 = p1->GetMother();
      continue;
    }
    if(iprim1 == 0 || iprim1 == 1){//colliding protons
      iprim1 = p1->GetMother();
      continue;
    }

    Int_t iprim2 = jPart;
    while(iprim2>-1){
      AliAODMCParticle *p2 = dynamic_cast<AliAODMCParticle*>(fMCArrayAOD->At(iprim2));
      Int_t pdg2 = p2->GetPdgCode();

      if(TMath::Abs(pdg2) <= 6 || TMath::Abs(pdg2) == 21){//reject quarks and gluons as parents.
        iprim2 = p2->GetMother();
        continue;
      }
      if(iprim2 == 0 || iprim2 == 1){//colliding protons
        iprim2 = p2->GetMother();
        continue;
      }

      if(iprim1==iprim2) return iprim1;
      iprim2 = p2->GetMother();

    }

    iprim1 = p1->GetMother();
  }

  return -1;
}
//________________________________________________________________________
void AliAnalysisTaskPHOSPi0EtaToGammaGamma::FillHistogramTH1(TList *list, const Char_t *hname, Double_t x, Double_t w, Option_t *opt) const
{
  //FillHistogram
  TH1 * hist = dynamic_cast<TH1*>(list->FindObject(hname));
  if(!hist){
    AliError(Form("can not find histogram (of instance TH1) <%s> ",hname));
    return;
  }
  else{
    TString optstring(opt);
    Double_t myweight = optstring.Contains("w") ? 1. : w;
    if(optstring.Contains("wx")){
      // use bin width as weight
      Int_t bin = hist->GetXaxis()->FindBin(x);
      // check if not overflow or underflow bin
      if(bin != 0 && bin != hist->GetXaxis()->GetNbins()){
        Double_t binwidth = hist->GetXaxis()->GetBinWidth(bin);
        myweight = w/binwidth;
      }
    }
    hist->Fill(x, myweight);
    return;
  }
}
//_____________________________________________________________________________
void AliAnalysisTaskPHOSPi0EtaToGammaGamma::FillHistogramTH2(TList *list, const Char_t *name, Double_t x, Double_t y, Double_t w, Option_t *opt) const
{
  TH2 * hist = dynamic_cast<TH2*>(list->FindObject(name));
  if(!hist){
    AliError(Form("can not find histogram (of instance TH2) <%s> ",name));
    return;
  }
  else{
    TString optstring(opt);
    Double_t myweight = optstring.Contains("w") ? 1. : w;
    if(optstring.Contains("wx")){
      Int_t binx = hist->GetXaxis()->FindBin(x);
      if(binx != 0 && binx != hist->GetXaxis()->GetNbins()) myweight *= 1./hist->GetXaxis()->GetBinWidth(binx);
    }
    if(optstring.Contains("wy")){
      Int_t biny = hist->GetYaxis()->FindBin(y);
      if(biny != 0 && biny != hist->GetYaxis()->GetNbins()) myweight *= 1./hist->GetYaxis()->GetBinWidth(biny);
    }
    hist->Fill(x, y, myweight);
    return;
  }
}
//_____________________________________________________________________________
void AliAnalysisTaskPHOSPi0EtaToGammaGamma::FillHistogramTH3(TList *list, const Char_t *name, Double_t x, Double_t y, Double_t z, Double_t w, Option_t *opt) const
{
  TH3 * hist = dynamic_cast<TH3*>(list->FindObject(name));
  if(!hist){
    AliError(Form("can not find histogram (of instance TH3) <%s> ",name));
    return;
  }
  else{
    TString optstring(opt);
    Double_t myweight = optstring.Contains("w") ? 1. : w;
    if(optstring.Contains("wx")){
      Int_t binx = hist->GetXaxis()->FindBin(x);
      if(binx != 0 && binx != hist->GetXaxis()->GetNbins()) myweight *= 1./hist->GetXaxis()->GetBinWidth(binx);
    }
    if(optstring.Contains("wy")){
      Int_t biny = hist->GetYaxis()->FindBin(y);
      if(biny != 0 && biny != hist->GetYaxis()->GetNbins()) myweight *= 1./hist->GetYaxis()->GetBinWidth(biny);
    }
    if(optstring.Contains("wz")){
      Int_t binz = hist->GetZaxis()->FindBin(z);
      if(binz != 0 && binz != hist->GetZaxis()->GetNbins()) myweight *= 1./hist->GetZaxis()->GetBinWidth(binz);
    }
    hist->Fill(x, y, z, myweight);
    return;
  }
}
//_____________________________________________________________________________
void AliAnalysisTaskPHOSPi0EtaToGammaGamma::FillProfile(TList *list, const Char_t *name, Double_t x, Double_t y) const
{
  TProfile * hist = dynamic_cast<TProfile*>(list->FindObject(name));
  if(!hist){
    AliError(Form("can not find histogram (of instance TH3) <%s> ",name));
    return;
  }
  else{
    hist->Fill(x,y);
    return;
  }

}
//_____________________________________________________________________________
void AliAnalysisTaskPHOSPi0EtaToGammaGamma::FillSparse(TList *list, const Char_t *name, Double_t *x, Double_t w) const
{
  THnSparse * hist = dynamic_cast<THnSparse*>(list->FindObject(name));
  if(!hist){
    AliError(Form("can not find histogram (of instance THnSparse) <%s> ",name));
    return;
  }
  else{
    hist->Fill(x,w);
    return;
  }

}
//_____________________________________________________________________________
AliPHOSGeometry *AliAnalysisTaskPHOSPi0EtaToGammaGamma::GetPHOSGeometry()
{
  AliAODEvent *aod = dynamic_cast<AliAODEvent*>(fEvent);
  AliESDEvent *esd = dynamic_cast<AliESDEvent*>(fEvent);
  Int_t RunNumber = fEvent->GetRunNumber();

  //Initialize the PHOS geometry

  fPHOSGeo = 0x0;

  if(fUsePHOSTender){
    if(RunNumber < 209122)//Run1
      fPHOSGeo = AliPHOSGeometry::GetInstance("IHEP") ;
    else//Run2
      fPHOSGeo = AliPHOSGeometry::GetInstance("Run2") ;
  }
  else{//PHOSTender is not applied.
    AliOADBContainer geomContainer("phosGeo");

    if(RunNumber < 209122)//Run1
      fPHOSGeo = AliPHOSGeometry::GetInstance("IHEP") ;
    else//Run2
      fPHOSGeo = AliPHOSGeometry::GetInstance("Run2") ;

    if(fIsMC){ //use excatly the same geometry as in simulation, stored in esd
      if(esd){
        for(Int_t mod=0; mod<6; mod++) {
          const TGeoHMatrix * m = esd->GetPHOSMatrix(mod);
          if(m){
            fPHOSGeo->SetMisalMatrix(m,mod) ;
            printf(".........Adding Matrix(%d), geo=%p\n",mod,fPHOSGeo);
            m->Print();
          }
        }
      }
      else if(aod){ //To be fixed
        geomContainer.InitFromFile("$ALICE_PHYSICS/OADB/PHOS/PHOSMCGeometry.root","PHOSMCRotationMatrixes");
        TObjArray *matrixes = (TObjArray*)geomContainer.GetObject(RunNumber,"PHOSRotationMatrixes");
        for(Int_t mod=0; mod<6; mod++){
          if(!matrixes->At(mod)) continue;
          fPHOSGeo->SetMisalMatrix(((TGeoHMatrix*)matrixes->At(mod)),mod);
          printf(".........Adding Matrix(%d), geo=%p\n",mod,fPHOSGeo);
          ((TGeoHMatrix*)matrixes->At(mod))->Print();
        }
      }
    }//end of MC
    else{ //Use best approaximation to real geometry
      geomContainer.InitFromFile("$ALICE_PHYSICS/OADB/PHOS/PHOSGeometry.root","PHOSRotationMatrixes");
      TObjArray *matrixes = (TObjArray*)geomContainer.GetObject(RunNumber,"PHOSRotationMatrixes");
      for(Int_t mod=0; mod<6; mod++){
        if(!matrixes->At(mod)) continue;
        fPHOSGeo->SetMisalMatrix(((TGeoHMatrix*)matrixes->At(mod)),mod);
        printf(".........Adding Matrix(%d), geo=%p\n",mod,fPHOSGeo);
        ((TGeoHMatrix*)matrixes->At(mod))->Print();
      }
    }//end of real data
  }

  return fPHOSGeo;

}
//_______________________________________________________________________________
AliStack *AliAnalysisTaskPHOSPi0EtaToGammaGamma::GetMCInfoESD() 
{
  AliStack *fStack = 0x0;
  AliVEventHandler* eventHandler = AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler();
  if(eventHandler){
    AliMCEventHandler* mcEventHandler = dynamic_cast<AliMCEventHandler*> (eventHandler);
    if(mcEventHandler) fStack = static_cast<AliMCEventHandler*>(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler())->MCEvent()->Stack();
  }

  if(!fStack) AliError("Could not get MC Stack!");

  return fStack;

}
//_______________________________________________________________________________
TClonesArray *AliAnalysisTaskPHOSPi0EtaToGammaGamma::GetMCInfoAOD() 
{
  TClonesArray *fMCArray = 0x0;
  AliAODInputHandler* aodHandler=dynamic_cast<AliAODInputHandler*>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());

  if(aodHandler){
    AliAODEvent *aod=aodHandler->GetEvent();
    if(aod){
      fMCArray = dynamic_cast<TClonesArray*>(aod->FindListObject(AliAODMCParticle::StdBranchName()));
      if (!fMCArray) AliError("Could not retrieve MC array!");
    }
    else AliError("Could not retrieve AOD event!");
  }

  return fMCArray;

}
//_______________________________________________________________________________
Int_t AliAnalysisTaskPHOSPi0EtaToGammaGamma::FindPrimaryMotherESD(Int_t label)
{
  const Int_t Nprimary = fMCArrayESD->GetNprimary();//this number contains only generated particles by event generator.
  //const Int_t Ntrack   = fMCArrayESD->GetNtrack();//this number contains generated particles by event generator + GEANT.
  Int_t tmp = label;
  while(tmp >= Nprimary){
    TParticle *p = (TParticle*)fMCArrayESD->Particle(tmp);
    tmp = p->GetMother(0);//first mother.
  }

  return tmp;
}
//_______________________________________________________________________________
void AliAnalysisTaskPHOSPi0EtaToGammaGamma::SetMCWeight()
{
  const Int_t multClust = fPHOSClusterArray->GetEntriesFast();

  Double_t weight = 1.;
  TF1 *f1Pi0Weight   = (TF1*)GetAdditionalPi0PtWeightFunction(fCentralityMain);
  TF1 *f1K0SWeight   = (TF1*)GetAdditionalK0SPtWeightFunction(fCentralityMain);
  TF1 *f1L0Weight    = (TF1*)GetAdditionalL0PtWeightFunction(fCentralityMain);
  TF1 *f1EtaWeight   = (TF1*)GetAdditionalEtaPtWeightFunction(fCentralityMain);
  //TF1 *f1GammaWeight = (TF1*)GetAdditionalGammaPtWeightFunction(fCentralityMain);

  Int_t primary = -1;
  Double_t TruePi0Pt = 0.;
  Double_t TrueEtaPt = 0.;
  //Double_t TrueGammaPt = 0.;
  Double_t TrueK0SPt = 0;
  Double_t TrueL0Pt = 0;

  for(Int_t iph=0;iph<multClust;iph++){
    AliCaloPhoton *ph = (AliCaloPhoton*)fPHOSClusterArray->At(iph);
    primary = ph->GetPrimary();
    weight = 1.;

    //weight is always defined as relative X/pi0 ratio, and pi0 weight is absolute number.
    if(IsFrom(primary,TruePi0Pt,111)) weight = f1Pi0Weight->Eval(TruePi0Pt);//pi0
    if(IsFrom(primary,TruePi0Pt,211)) weight = f1Pi0Weight->Eval(TruePi0Pt);//pi+/-

    if(IsFrom(primary,TrueK0SPt,310)) weight = f1K0SWeight->Eval(TrueK0SPt) * f1Pi0Weight->Eval(TrueK0SPt);//K0S
    if(IsFrom(primary,TrueK0SPt,130)) weight = f1K0SWeight->Eval(TrueK0SPt) * f1Pi0Weight->Eval(TrueK0SPt);//K0L
    if(IsFrom(primary,TrueK0SPt,321)) weight = f1K0SWeight->Eval(TrueK0SPt) * f1Pi0Weight->Eval(TrueK0SPt);//K+/-

    if(IsFrom(primary,TrueEtaPt,221)) weight = f1EtaWeight->Eval(TrueEtaPt) * f1Pi0Weight->Eval(TrueEtaPt);//eta

    if(IsFrom(primary,TrueL0Pt,3122)) weight = f1L0Weight->Eval(TrueL0Pt)   * f1Pi0Weight->Eval(TrueL0Pt);//Lambda

    ph->SetWeight(weight);

  }//end of cluster loop

}
//_______________________________________________________________________________
void AliAnalysisTaskPHOSPi0EtaToGammaGamma::FillRejectionFactorMB()
{
  //search for 0PH0 event in MB for the rejection factor.
  //the best senario of how to estimate trigger rejection factor is reading ALICE electric log book.
  //but, fake trigger might be there.

  Int_t L1Hinput = 7;
  Int_t L1Minput = 6;
  Int_t L1Linput = 5;
  Int_t L0input  = 17;//LHC17

  if(fRunNumber <= 246994)  L0input = 9;
  const Int_t Nmod=5;//number of mpdules
  const Bool_t IsPrivateTRUBadMap = fPHOSTriggerHelper->IsPrivateTRUBadMap();

  if(!fPHOSTriggerHelperL0 ){
    fPHOSTriggerHelperL0  = new AliPHOSTriggerHelper(-1      ,L0input,kFALSE);
    if(IsPrivateTRUBadMap){
      for(Int_t imod=1;imod<Nmod;imod++) fPHOSTriggerHelperL0->SetPHOSTRUBadMap(imod,fPHOSTriggerHelper->GetPHOSTRUBadMap(imod));
    }
  }
  if(!fPHOSTriggerHelperL1H){
    fPHOSTriggerHelperL1H = new AliPHOSTriggerHelper(L1Hinput,-1     ,kFALSE);
    if(IsPrivateTRUBadMap){
      for(Int_t imod=1;imod<Nmod;imod++) fPHOSTriggerHelperL1H->SetPHOSTRUBadMap(imod,fPHOSTriggerHelper->GetPHOSTRUBadMap(imod));
    }
  }
  if(!fPHOSTriggerHelperL1M){
    fPHOSTriggerHelperL1M = new AliPHOSTriggerHelper(L1Minput,-1     ,kFALSE);
    if(IsPrivateTRUBadMap){
      for(Int_t imod=1;imod<Nmod;imod++) fPHOSTriggerHelperL1M->SetPHOSTRUBadMap(imod,fPHOSTriggerHelper->GetPHOSTRUBadMap(imod));
    }
  }
  if(!fPHOSTriggerHelperL1L){
    fPHOSTriggerHelperL1L = new AliPHOSTriggerHelper(L1Linput,-1     ,kFALSE);
    if(IsPrivateTRUBadMap){
      for(Int_t imod=1;imod<Nmod;imod++) fPHOSTriggerHelperL1L->SetPHOSTRUBadMap(imod,fPHOSTriggerHelper->GetPHOSTRUBadMap(imod));
    }
  }

  Bool_t Is0PH0fired = fEvent->GetHeader()->GetL0TriggerInputs() & 1 << (L0input  - 1);//trigger input -1
  Bool_t Is1PHHfired = fEvent->GetHeader()->GetL1TriggerInputs() & 1 << (L1Hinput - 1);//trigger input -1
  Bool_t Is1PHMfired = fEvent->GetHeader()->GetL1TriggerInputs() & 1 << (L1Minput - 1);//trigger input -1
  Bool_t Is1PHLfired = fEvent->GetHeader()->GetL1TriggerInputs() & 1 << (L1Linput - 1);//trigger input -1

  AliInfo(Form("Energy threshold of PHOS trigger is %3.2f GeV",fEnergyThreshold));

  Bool_t Is0PH0matched = fPHOSTriggerHelperL0 ->IsPHI7(fEvent,fPHOSClusterCuts,fEmin,fEnergyThreshold,fUseCoreEnergy);
  Bool_t Is1PHHmatched = fPHOSTriggerHelperL1H->IsPHI7(fEvent,fPHOSClusterCuts,fEmin,fEnergyThreshold,fUseCoreEnergy);
  Bool_t Is1PHMmatched = fPHOSTriggerHelperL1M->IsPHI7(fEvent,fPHOSClusterCuts,fEmin,fEnergyThreshold,fUseCoreEnergy);
  Bool_t Is1PHLmatched = fPHOSTriggerHelperL1L->IsPHI7(fEvent,fPHOSClusterCuts,fEmin,fEnergyThreshold,fUseCoreEnergy);

  if(Is0PH0fired && Is0PH0matched) FillHistogramTH1(fOutputContainer,"hEventSummary", 7);//0PH0
  if(Is1PHLfired && Is1PHHmatched) FillHistogramTH1(fOutputContainer,"hEventSummary", 8);//1PHL
  if(Is1PHMfired && Is1PHMmatched) FillHistogramTH1(fOutputContainer,"hEventSummary", 9);//1PHM
  if(Is1PHHfired && Is1PHLmatched) FillHistogramTH1(fOutputContainer,"hEventSummary",10);//1PHH

  if(Is0PH0fired) FillTriggerInfoMB(fPHOSTriggerHelperL0,-1);
  if(Is1PHLfired) FillTriggerInfoMB(fPHOSTriggerHelperL1L,2);
  if(Is1PHMfired) FillTriggerInfoMB(fPHOSTriggerHelperL1M,1);
  if(Is1PHHfired) FillTriggerInfoMB(fPHOSTriggerHelperL1H,0);

}
//_______________________________________________________________________________
void AliAnalysisTaskPHOSPi0EtaToGammaGamma::FillTriggerInfoMB(AliPHOSTriggerHelper *helper, const Int_t L1)
{
  Double_t energy = 0;
  //Double_t pT = 0;
  Float_t position[3] = {};
  Int_t trgrelId[4]={};
  Int_t tmod=0; //"Online" module number
  Int_t trgabsId=0; //bottom-left 4x4 edge cell absId [1-64], [1-56]
  Int_t trgmodule=0,trgcellx=0,trgcellz=0;


  Int_t relId[4]={};
  Int_t module=0,cellx=0,cellz=0,tru=0;

  TString type = "";
  if(L1== -1) type = "L0";
  else if(L1== 0) type = "L1H";
  else if(L1== 1) type = "L1M";
  else if(L1== 2) type = "L1L";
  else{
    AliInfo(Form("Your choice of trigger type %d is wrong. return.",L1));
    return;
  }

  const Int_t multClust = fPHOSClusterArray->GetEntriesFast();

  for(Int_t i=0;i<multClust;i++){
    AliCaloPhoton *ph = (AliCaloPhoton*)fPHOSClusterArray->At(i);
    if(!fPHOSClusterCuts->AcceptPhoton(ph)) continue;
    if(!CheckMinimumEnergy(ph)) continue;

    if(!ph->IsTOFOK()) continue;//MB data requires timing cut.

    //pT = ph->Pt();
    energy = ph->Energy();

    if(fUseCoreEnergy){
      //pT = (ph->GetMomV2())->Pt();
      energy = (ph->GetMomV2())->Energy();
    }

    position[0] = ph->EMCx();
    position[1] = ph->EMCy();
    position[2] = ph->EMCz();
    TVector3 global1(position);
    fPHOSGeo->GlobalPos2RelId(global1,relId);

    module = relId[0];
    cellx  = relId[2];
    cellz  = relId[3];
    tru = helper->WhichTRU(cellx,cellz);

    if(module < 1 || module > 4){
      AliError(Form("Wrong module number %d",module));
      return;
    }
    AliInfo(Form("cluster position M:%d, X:%d , Z:%d , TRU:%d",module,cellx,cellz,tru));

    AliVCaloTrigger* trg = fEvent->GetCaloTrigger("PHOS");
    trg->Reset();

    while(trg->Next()){

      // L1 threshold: -1-L0, 0-PHH, 1-PHM, 2-PHL
      if(trg->GetL1TimeSum() != L1) continue;

      trg->GetPosition(tmod,trgabsId);//tmod is online module numbering. i.e., tmod=1 means half module.

      fPHOSGeo->AbsToRelNumbering(trgabsId,trgrelId);
      //for offline numbering, relId should be used.

      trgmodule = trgrelId[0];
      trgcellx  = trgrelId[2];
      trgcellz  = trgrelId[3];

      AliInfo(Form("fired TRU channel M:%d, X:%d , Z:%d",trgmodule,trgcellx,trgcellz));

      if(trgmodule < 1 || trgmodule > 4){
        AliError(Form("Wrong module number %d",trgmodule));
        return;
      }

      if(helper->IsMatched(trgrelId,relId)){
        FillHistogramTH1(fOutputContainer,Form("hMatchedClusterEnergy%s",type.Data()),energy);//matched cluster
        FillHistogramTH1(fOutputContainer,Form("hMatchedClusterEnergy%sM%dTRU%d",type.Data(),module,tru),energy);//matched cluster
        break;//exit from trigger patch loop
      }
    }//end of while trigger patch loop

  }//end of cluster loop

}
//_______________________________________________________________________________
const AliQnCorrectionsQnVector *AliAnalysisTaskPHOSPi0EtaToGammaGamma::GetQnVectorFromList(const TList *list, const char* subdetector, const char *expcorr, const char *altcorr)
{
  AliQnCorrectionsQnVector *theQnVector = NULL;

  TList *pQvecList = dynamic_cast<TList*> (list->FindObject(subdetector));
  if(pQvecList != NULL){//sub detector is found
    if(TString(expcorr).EqualTo("latest")) theQnVector = (AliQnCorrectionsQnVector*) pQvecList->First();
    else                                   theQnVector = (AliQnCorrectionsQnVector*) pQvecList->FindObject(expcorr);
    //theQnVector = (AliQnCorrectionsQnVector*) pQvecList->FindObject(expcorr);

    if(theQnVector == NULL || !(theQnVector->IsGoodQuality()) || !(theQnVector->GetN() != 0)){ //the Qn vector for the expected correction is not found
      AliInfo(Form("expected correction (%s) is not found. use %s as an alternative step in %s.",expcorr,altcorr,subdetector));
      if(TString(altcorr).EqualTo("latest")) theQnVector = (AliQnCorrectionsQnVector*) pQvecList->First();
      else                                   theQnVector = (AliQnCorrectionsQnVector*) pQvecList->FindObject(altcorr);
        //theQnVector = (AliQnCorrectionsQnVector*) pQvecList->FindObject(altcorr);
    }

    //check the Qn vector quality
    if(!(theQnVector->IsGoodQuality()) || !(theQnVector->GetN() != 0)) theQnVector = NULL; //bad quality, discarded

  }
  return theQnVector;

}
//_______________________________________________________________________________
void AliAnalysisTaskPHOSPi0EtaToGammaGamma::DoNonLinearityStudy()
{
  const Int_t multClust = fPHOSClusterArray->GetEntriesFast();
  TLorentzVector p12, p12core;
  Double_t m12=0,pt12=0;
  Double_t e1=0,e2=0;
  Double_t weight = 1., w1 = 1., w2 = 1.;
  Int_t primary1 = -1;
  Int_t primary2 = -1;
  Int_t commonID = -1;

  for(Int_t ia=0;ia<7;ia++){
    for(Int_t ib=0;ib<7;ib++){

      for(Int_t i1=0;i1<multClust-1;i1++){
        AliCaloPhoton *ph1 = (AliCaloPhoton*)fPHOSClusterArray->At(i1);
        if(!fPHOSClusterCuts->AcceptPhoton(ph1)) continue;
        if(!CheckMinimumEnergy(ph1)) continue;

        for(Int_t i2=i1+1;i2<multClust;i2++){
          AliCaloPhoton *ph2 = (AliCaloPhoton*)fPHOSClusterArray->At(i2);
          if(!fPHOSClusterCuts->AcceptPhoton(ph2)) continue;
          if(!CheckMinimumEnergy(ph2)) continue;

          e1 = ph1->Energy();
          e2 = ph2->Energy();

          p12 = *ph1*fNonLin[ia][ib]->Eval(e1) + *ph2*fNonLin[ia][ib]->Eval(e2);
          m12 = p12.M();
          pt12 = p12.Pt();

          if(fUseCoreEnergy){
            e1 = (ph1->GetMomV2())->Energy();
            e2 = (ph2->GetMomV2())->Energy();

            p12core = *(ph1->GetMomV2())*fNonLin[ia][ib]->Eval(e1) + *(ph2->GetMomV2())*fNonLin[ia][ib]->Eval(e2);
            m12 = p12core.M();
            pt12 = p12core.Pt();
          }

          weight = 1.;
          if(fIsMC){
            w1= ph1->GetWeight();
            primary1 = ph1->GetPrimary();

            w2 = ph2->GetWeight();
            primary2 = ph2->GetPrimary();

            commonID = FindCommonParent(primary1,primary2);
            if(commonID > -1) weight = w1;
            else weight = w1*w2;

          }//end of if fIsMC

          FillHistogramTH2(fOutputContainer,Form("hMgg_a%d_b%d",ia,ib),m12,pt12,weight);

        }//end of ph2
      }//end of ph1
    }//end of ib
  }//end of ia

  TList *prevPHOS = fPHOSEvents[fZvtx][fEPBin];

  for(Int_t ia=0;ia<7;ia++){
    for(Int_t ib=0;ib<7;ib++){

      for(Int_t i1=0;i1<multClust;i1++){
        AliCaloPhoton *ph1 = (AliCaloPhoton*)fPHOSClusterArray->At(i1);
        if(!fPHOSClusterCuts->AcceptPhoton(ph1)) continue;
        if(!CheckMinimumEnergy(ph1)) continue;

        for(Int_t ev=0;ev<prevPHOS->GetSize();ev++){
          TClonesArray *mixPHOS = static_cast<TClonesArray*>(prevPHOS->At(ev));

          for(Int_t i2=0;i2<mixPHOS->GetEntriesFast();i2++){
            AliCaloPhoton *ph2 = (AliCaloPhoton*)mixPHOS->At(i2);
            if(!fPHOSClusterCuts->AcceptPhoton(ph2)) continue;
            if(!CheckMinimumEnergy(ph2)) continue;

            e1 = ph1->Energy();
            e2 = ph2->Energy();

            p12 = *ph1*fNonLin[ia][ib]->Eval(e1) + *ph2*fNonLin[ia][ib]->Eval(e2);
            m12 = p12.M();
            pt12 = p12.Pt();

            if(fUseCoreEnergy){
              e1 = (ph1->GetMomV2())->Energy();
              e2 = (ph2->GetMomV2())->Energy();
              p12core = *(ph1->GetMomV2())*fNonLin[ia][ib]->Eval(e1) + *(ph2->GetMomV2())*fNonLin[ia][ib]->Eval(e2);
              m12 = p12core.M();
              pt12 = p12core.Pt();
            }

            FillHistogramTH2(fOutputContainer,Form("hMixMgg_a%d_b%d",ia,ib),m12,pt12);

          }//end of mix
        }//end of ph2
      }//end of ph1
    }//end of ib
  }//end of ia

}
//_______________________________________________________________________________
Bool_t AliAnalysisTaskPHOSPi0EtaToGammaGamma::ExtractQnVector()
{
  if(fHarmonics < 0){
    AliError(Form("Qn Flow vector correction flag is ON, but fHarmonics is not set. (it is %d now).",fHarmonics));
    return kFALSE;
  }

  TList* qnlist = fFlowQnVectorMgr->GetQnVectorList();

  const AliQnCorrectionsQnVector *QnVectorTPCDet[3];
  Double_t TPCEP[3] = {};
  for(Int_t i=0;i<3;i++){
    QnVectorTPCDet[i] = GetQnVectorFromList(qnlist,Form("%s%s",fTPCEPName[i].Data(),fQNormalization.Data()),"latest","plain");
    if(!QnVectorTPCDet[i]){
      AliInfo("Event is rejected because event plane is not found or bad event plane quality in TPC.");
      return kFALSE;//Qn vector correction does not exist or bad quality.
    }
    TPCEP[i] = QnVectorTPCDet[i]->EventPlane(fHarmonics);
    if(TPCEP[i] < 0) TPCEP[i] += 2./(Double_t) fHarmonics * TMath::Pi();
    FillHistogramTH2(fOutputContainer,Form("hCentrality%svsEventPlane%s%s",fEstimator.Data(),fTPCEPName[i].Data(),fQNormalization.Data()),fCentralityMain,TPCEP[i]);
    AliInfo(Form("harmonics %d | TPC sub detector name %s%s : event plane = %f (rad).",fHarmonics,fTPCEPName[i].Data(),fQNormalization.Data(),TPCEP[i]));
  }

  const AliQnCorrectionsQnVector *QnVectorV0Det[3];
  Double_t V0EP[3]  = {};
  for(Int_t i=0;i<3;i++){
    QnVectorV0Det[i]  = GetQnVectorFromList(qnlist,Form("%s%s",fV0EPName[i].Data(),fQNormalization.Data()),"latest","raw");
    if(!QnVectorV0Det[i]){
      AliInfo("Event is rejected because event plane is not found or bad event plane quality in VZERO.");
      return kFALSE;//Qn vector correction does not exist or bad quality.
    }
    V0EP[i] = QnVectorV0Det[i]->EventPlane(fHarmonics);
    if(V0EP[i] < 0)  V0EP[i]  += 2./(Double_t) fHarmonics * TMath::Pi();
    FillHistogramTH2(fOutputContainer,Form("hCentrality%svsEventPlane%s%s",fEstimator.Data(),fV0EPName[i].Data(),fQNormalization.Data()),fCentralityMain,V0EP[i]);
    AliInfo(Form("harmonics %d | V0  sub detector name %s%s : event plane = %f (rad).",fHarmonics,fV0EPName[i].Data(),fQNormalization.Data(),V0EP[i]));
  }

  //0 < event plane < 2*pi/fHarmonics.
  Double_t EP2 = -999; 
  Double_t EP3 = -999; 

  Double_t Q1[2] = {};//for Main
  Double_t Q2[2] = {};//for Sub1
  Double_t Q3[2] = {};//for Sub2

  if(fQnDetectorMain == AliAnalysisTaskPHOSPi0EtaToGammaGamma::kFullTPC){
    Q1[0] = QnVectorTPCDet[0]->Qx(fHarmonics);//FullTPC
    Q1[1] = QnVectorTPCDet[0]->Qy(fHarmonics);//FullTPC
    Q2[0] = QnVectorV0Det[1]->Qx(fHarmonics);//V0A
    Q2[1] = QnVectorV0Det[1]->Qy(fHarmonics);//V0A
    Q3[0] = QnVectorV0Det[2]->Qx(fHarmonics);//V0C
    Q3[1] = QnVectorV0Det[2]->Qy(fHarmonics);//V0C
    fEventPlane = TPCEP[0];
    EP2 = V0EP[1];
    EP3 = V0EP[2];
  }
  else if(fQnDetectorMain == AliAnalysisTaskPHOSPi0EtaToGammaGamma::kTPCNegEta){
    Q1[0] = QnVectorTPCDet[1]->Qx(fHarmonics);//TPCNegEta
    Q1[1] = QnVectorTPCDet[1]->Qy(fHarmonics);//TPCNegEta
    Q2[0] = QnVectorV0Det[1]->Qx(fHarmonics);//V0A
    Q2[1] = QnVectorV0Det[1]->Qy(fHarmonics);//V0A
    Q3[0] = QnVectorV0Det[2]->Qx(fHarmonics);//V0C
    Q3[1] = QnVectorV0Det[2]->Qy(fHarmonics);//V0C
    fEventPlane = TPCEP[2];
    EP2 = V0EP[1];
    EP3 = V0EP[2];
  }
  else if(fQnDetectorMain == AliAnalysisTaskPHOSPi0EtaToGammaGamma::kTPCPosEta){
    Q1[0] = QnVectorTPCDet[2]->Qx(fHarmonics);//TPCPosEta
    Q1[1] = QnVectorTPCDet[2]->Qy(fHarmonics);//TPCPosEta
    Q2[0] = QnVectorV0Det[1]->Qx(fHarmonics);//V0A
    Q2[1] = QnVectorV0Det[1]->Qy(fHarmonics);//V0A
    Q3[0] = QnVectorV0Det[2]->Qx(fHarmonics);//V0C
    Q3[1] = QnVectorV0Det[2]->Qy(fHarmonics);//V0C
    fEventPlane = TPCEP[1];
    EP2 = V0EP[1];
    EP3 = V0EP[2];
  }
  else if(fQnDetectorMain == AliAnalysisTaskPHOSPi0EtaToGammaGamma::kFullV0){
    Q1[0] = QnVectorV0Det[0]->Qx(fHarmonics);//FullV0
    Q1[1] = QnVectorV0Det[0]->Qy(fHarmonics);//FullV0
    Q2[0] = QnVectorTPCDet[1]->Qx(fHarmonics);//TPCNegEta
    Q2[1] = QnVectorTPCDet[1]->Qy(fHarmonics);//TPCNegEta
    Q3[0] = QnVectorTPCDet[2]->Qx(fHarmonics);//TPCPosEta
    Q3[1] = QnVectorTPCDet[2]->Qy(fHarmonics);//TPCPosEta
    fEventPlane = V0EP[0];
    EP2 = TPCEP[1];
    EP3 = TPCEP[2];
  }
  else if(fQnDetectorMain == AliAnalysisTaskPHOSPi0EtaToGammaGamma::kV0A){
    Q1[0] = QnVectorV0Det[1]->Qx(fHarmonics);//V0A
    Q1[1] = QnVectorV0Det[1]->Qy(fHarmonics);//V0A
    Q2[0] = QnVectorV0Det[2]->Qx(fHarmonics);//V0C
    Q2[1] = QnVectorV0Det[2]->Qy(fHarmonics);//V0C
    Q3[0] = QnVectorTPCDet[0]->Qx(fHarmonics);//full acceptance of TPC
    Q3[1] = QnVectorTPCDet[0]->Qy(fHarmonics);//full acceptance of TPC
    fEventPlane = V0EP[1];
    EP2 = V0EP[2];
    EP3 = TPCEP[0];
  }
  else if(fQnDetectorMain == AliAnalysisTaskPHOSPi0EtaToGammaGamma::kV0C){
    Q1[0] = QnVectorV0Det[2]->Qx(fHarmonics);//V0C
    Q1[1] = QnVectorV0Det[2]->Qy(fHarmonics);//V0C
    Q2[0] = QnVectorV0Det[1]->Qx(fHarmonics);//V0A
    Q2[1] = QnVectorV0Det[1]->Qy(fHarmonics);//V0A
    Q3[0] = QnVectorTPCDet[0]->Qx(fHarmonics);//full acceptance of TPC
    Q3[1] = QnVectorTPCDet[0]->Qy(fHarmonics);//full acceptance of TPC
    fEventPlane = V0EP[2];
    EP2 = V0EP[1];
    EP3 = TPCEP[0];
  }

  fQVector1.Set(Q1[0],Q1[1]);
  TVector2 QVector2(Q2[0],Q2[1]);
  TVector2 QVector3(Q3[0],Q3[1]);


  FillHistogramTH2(fOutputContainer,Form("hCentrality%svsNormQ1",fEstimator.Data()),fCentralityMain,fQVector1.Mod());
  FillHistogramTH2(fOutputContainer,Form("hCentrality%svsNormQ2",fEstimator.Data()),fCentralityMain, QVector2.Mod());
  FillHistogramTH2(fOutputContainer,Form("hCentrality%svsNormQ3",fEstimator.Data()),fCentralityMain, QVector3.Mod());


  FillHistogramTH2(fOutputContainer,Form("hCentrality%svsQ1x",fEstimator.Data()),fCentralityMain,Q1[0]);
  FillHistogramTH2(fOutputContainer,Form("hCentrality%svsQ1y",fEstimator.Data()),fCentralityMain,Q1[1]);
  FillHistogramTH2(fOutputContainer,Form("hCentrality%svsQ2x",fEstimator.Data()),fCentralityMain,Q2[0]);
  FillHistogramTH2(fOutputContainer,Form("hCentrality%svsQ2y",fEstimator.Data()),fCentralityMain,Q2[1]);
  FillHistogramTH2(fOutputContainer,Form("hCentrality%svsQ3x",fEstimator.Data()),fCentralityMain,Q3[0]);
  FillHistogramTH2(fOutputContainer,Form("hCentrality%svsQ3y",fEstimator.Data()),fCentralityMain,Q3[1]);

  Double_t sp12 = fQVector1 *  QVector2;//scalar product between Q1 vector and Q2 vector
  Double_t sp23 =  QVector2 *  QVector3;//scalar product between Q2 vector and Q3 vector
  Double_t sp31 =  QVector3 * fQVector1;//scalar product between Q3 vector and Q1 vector

  FillHistogramTH2(fOutputContainer,Form("hCentrality%svsSPQ1Q2",fEstimator.Data()),fCentralityMain,sp12);//mean value of sp is denominator of SP method
  FillHistogramTH2(fOutputContainer,Form("hCentrality%svsSPQ2Q3",fEstimator.Data()),fCentralityMain,sp23);//mean value of sp is denominator of SP method
  FillHistogramTH2(fOutputContainer,Form("hCentrality%svsSPQ3Q1",fEstimator.Data()),fCentralityMain,sp31);//mean value of sp is denominator of SP method
  AliInfo(Form("Q1x = %e , Q1y = %e , Q2x = %e , Q2y = %e , Q3x = %e , Q3y = %e ,  SP12 = %e ,  SP23 = %e ,  SP31 = %e",Q1[0],Q1[1],Q2[0],Q2[1],Q3[0],Q3[1],sp12,sp23,sp31));
  
  const Double_t delta = 2. * TMath::Pi() / Double_t(fHarmonics) / 12.;
  fEPBin = (Int_t)((fEventPlane) / delta);//it should be 0-11.
  if(fEPBin < 0)  fEPBin =  0;//protection to avoid fEPBin = -1.
  if(fEPBin > 11) fEPBin = 11;//protection to avoid fEPBin = 12.
  
  //for event plane resolution
  FillHistogramTH2(fOutputContainer,Form("hCentrality%svsCosDeltaEventPlane12",fEstimator.Data()),fCentralityMain,TMath::Cos(fHarmonics * (fEventPlane - EP2)));
  FillHistogramTH2(fOutputContainer,Form("hCentrality%svsCosDeltaEventPlane23",fEstimator.Data()),fCentralityMain,TMath::Cos(fHarmonics * (EP2         - EP3)));
  FillHistogramTH2(fOutputContainer,Form("hCentrality%svsCosDeltaEventPlane31",fEstimator.Data()),fCentralityMain,TMath::Cos(fHarmonics * (EP3 - fEventPlane)));
  FillHistogramTH2(fOutputContainer,Form("hCentrality%svsSinDeltaEventPlane12",fEstimator.Data()),fCentralityMain,TMath::Sin(fHarmonics * (fEventPlane - EP2)));
  FillHistogramTH2(fOutputContainer,Form("hCentrality%svsSinDeltaEventPlane23",fEstimator.Data()),fCentralityMain,TMath::Sin(fHarmonics * (EP2         - EP3)));
  FillHistogramTH2(fOutputContainer,Form("hCentrality%svsSinDeltaEventPlane31",fEstimator.Data()),fCentralityMain,TMath::Sin(fHarmonics * (EP3 - fEventPlane)));

  return kTRUE;
}
//_______________________________________________________________________________
Bool_t AliAnalysisTaskPHOSPi0EtaToGammaGamma::Are2GammasInPHOSAcceptance(Int_t id)
{
  if(!fMCArrayAOD) return kFALSE;
  AliAODMCParticle *p = (AliAODMCParticle*)fMCArrayAOD->At(id);
 
  Int_t id0 = p->GetDaughterLabel(0);
  Int_t id1 = p->GetDaughterLabel(1);

  if(id0 < 0 || id1 < 0) return kFALSE;
 
  AliAODMCParticle *d0 = (AliAODMCParticle*)fMCArrayAOD->At(id0);
  //printf("d0 | pdg = %d , eta = %e , phi = %e\n",d0->GetPdgCode(),d0->Eta(),d0->Phi());

  AliAODMCParticle *d1 = (AliAODMCParticle*)fMCArrayAOD->At(id1);
  //printf("d1 | pdg = %d , eta = %e , phi = %e\n",d1->GetPdgCode(),d1->Eta(),d1->Phi());

  if(d0->GetPdgCode() != 22 || d1->GetPdgCode() != 22) return kFALSE;
  if(TMath::Abs(d0->Eta()) > 0.13 || TMath::Abs(d1->Eta()) > 0.13) return kFALSE;

  if(d0->Phi() < 260.*TMath::DegToRad() || 320.*TMath::DegToRad() < d0->Phi()) return kFALSE;
  if(d1->Phi() < 260.*TMath::DegToRad() || 320.*TMath::DegToRad() < d1->Phi()) return kFALSE;

  return kTRUE;
}
//_______________________________________________________________________________
void AliAnalysisTaskPHOSPi0EtaToGammaGamma::FillTrackMatching()
{

  const Int_t multClust = fPHOSClusterArray->GetEntriesFast();

  Double_t r=999, trackPt=0;
  Double_t pT = 0;
  Int_t charge = 0;
  Double_t nsigmaElectron = 999;
  Double_t nsigmaPion = 999;
  Double_t nsigmaKaon = 999;
  Double_t nsigmaProton = 999;
  Bool_t pidPion = kFALSE;
  Bool_t pidKaon = kFALSE;
  Bool_t pidProton = kFALSE;
  Bool_t pidElectron = kFALSE;
  const Double_t NsigmaCut = 3;
  Bool_t isGlobal = kFALSE;

  for(Int_t i1=0;i1<multClust;i1++){
    AliCaloPhoton *ph = (AliCaloPhoton*)fPHOSClusterArray->At(i1);
    //if(!fPHOSClusterCuts->AcceptPhoton(ph)) continue;
    if(!CheckMinimumEnergy(ph)) continue;
    //if(!fIsMC && !ph->IsTOFOK()) continue;

    if(fIsPHOSTriggerAnalysis){
      if( fIsMC && !fPHOSTriggerHelper->IsOnActiveTRUChannel(ph)) continue;//only for MC
      if(!fIsMC && !ph->IsTrig()) continue;//it is meaningless to focus on photon without fired trigger in PHOS triggered data.
    }
    if(fForceActiveTRU && !fPHOSTriggerHelper->IsOnActiveTRUChannel(ph)) continue;

    pT = ph->Pt();
    if(fUseCoreEnergy) pT = (ph->GetMomV2())->Pt();

    AliVCluster *clu = (AliVCluster*)ph->GetCluster();
    AliVTrack *track = 0x0;
    r = 999;
    isGlobal = kFALSE;

    if(fESDEvent){
      Int_t trackindex = clu->GetTrackMatchedIndex();
      if(trackindex > 0){
        track = (AliVTrack*)(fEvent->GetTrack(trackindex));
        isGlobal = fESDtrackCutsGlobal->AcceptTrack(dynamic_cast<AliESDtrack*>(track));
      }//end of track matching
    }//end of ESD
    else if(fAODEvent){
      if(clu->GetNTracksMatched() > 0){
        track = dynamic_cast<AliVTrack*>(clu->GetTrackMatched(0));
        isGlobal = dynamic_cast<AliAODTrack*>(track)->TestFilterMask(AliAODTrack::kTrkGlobalNoDCA);//standard cuts with very loose DCA cut
      }//end of track matching
    }//end of AOD

    pidPion = kFALSE;
    pidKaon = kFALSE;
    pidProton = kFALSE;
    pidElectron = kFALSE;
    r = ph->GetNsigmaCPV();

    if(track){
      trackPt = track->Pt();
      charge = track->Charge();

      nsigmaElectron = fPIDResponse->NumberOfSigmasTPC(track,AliPID::kElectron);//(-2,3) is for real electron analysis.
      nsigmaPion     = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track,AliPID::kPion));
      nsigmaKaon     = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track,AliPID::kKaon));
      nsigmaProton   = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track,AliPID::kProton));

      if((TMath::Abs(nsigmaElectron) < nsigmaPion)     && (TMath::Abs(nsigmaElectron) < nsigmaKaon) && (TMath::Abs(nsigmaElectron) < nsigmaProton) && (-2 < nsigmaElectron && nsigmaElectron < 3)) pidElectron = kTRUE;
      if((nsigmaPion     < TMath::Abs(nsigmaElectron)) && (nsigmaPion     < nsigmaKaon) && (nsigmaPion     < nsigmaProton) && (nsigmaPion     < NsigmaCut)) pidPion = kTRUE;
      if((nsigmaKaon     < TMath::Abs(nsigmaElectron)) && (nsigmaKaon     < nsigmaPion) && (nsigmaKaon     < nsigmaProton) && (nsigmaKaon     < NsigmaCut)) pidKaon = kTRUE;
      if((nsigmaProton   < TMath::Abs(nsigmaElectron)) && (nsigmaProton   < nsigmaPion) && (nsigmaProton   < nsigmaKaon)   && (nsigmaProton   < NsigmaCut)) pidProton = kTRUE;

      if(pidElectron){
        FillHistogramTH2(fOutputContainer,"hRvsTrackPt_Electron",trackPt,r);
        FillHistogramTH2(fOutputContainer,"hRvsClusterPt_Electron",pT,r);
      }
      else if(pidPion){
        FillHistogramTH2(fOutputContainer,"hRvsTrackPt_Pion",trackPt,r);
        FillHistogramTH2(fOutputContainer,"hRvsClusterPt_Pion",pT,r);
      }
      else if(pidKaon){
        FillHistogramTH2(fOutputContainer,"hRvsTrackPt_Kaon",trackPt,r);
        FillHistogramTH2(fOutputContainer,"hRvsClusterPt_Kaon",pT,r);
      }
      else if(pidProton){
        if(charge > 0){
          FillHistogramTH2(fOutputContainer,"hRvsTrackPt_Proton",trackPt,r);
          FillHistogramTH2(fOutputContainer,"hRvsClusterPt_Proton",pT,r);
        }
        else{
          FillHistogramTH2(fOutputContainer,"hRvsTrackPt_AntiProton",trackPt,r);
          FillHistogramTH2(fOutputContainer,"hRvsClusterPt_AntiProton",pT,r);
        }
      }
    }//end of track matching

    //FillHistogramTH2(fOutputContainer,"hRvsClusterPt",pT,r);
 
  }//end of cluster loop

}
//_______________________________________________________________________________
void AliAnalysisTaskPHOSPi0EtaToGammaGamma::FillMixTrackMatching()
{
  TList *prevPHOS = fPHOSEvents[fZvtx][fEPBin];
  Float_t position[3] = {0,0,0};
  Int_t relId[4]={0,0,0,0};
  Double_t pT = 0;
  Double_t pttrack=0.;
  Int_t charge=0;
  Double_t cpv = 999;
  Double_t dx=999.,dz=999.;

  Double_t nsigmaElectron = 999;
  Double_t nsigmaPion = 999;
  Double_t nsigmaKaon = 999;
  Double_t nsigmaProton = 999;
  Bool_t pidPion = kFALSE;
  Bool_t pidKaon = kFALSE;
  Bool_t pidProton = kFALSE;
  Bool_t pidElectron = kFALSE;
  const Double_t NsigmaCut = 3;
  Bool_t isGlobal = kFALSE;

  AliPHOSTenderTask *PHOSTenderTask = dynamic_cast<AliPHOSTenderTask*>(AliAnalysisManager::GetAnalysisManager()->GetTask("PHOSTenderTask"));
  if(!PHOSTenderTask){
    AliInfo("PHOSTenderTask does not exist!");
    return;
  }
  else{
    AliPHOSTenderSupply *supply = (AliPHOSTenderSupply*)PHOSTenderTask->GetPHOSTenderSupply();
    //track matching by AliPHOSTenderSupply

    Int_t Nev = prevPHOS->GetSize();
    if(Nev > fNMixTrack) Nev = fNMixTrack;//set max N events.

    for(Int_t iev=0;iev<Nev;iev++){
      TClonesArray *mixPHOS = static_cast<TClonesArray*>(prevPHOS->At(iev));

      for(Int_t iph=0;iph<mixPHOS->GetEntriesFast();iph++){
        AliCaloPhoton *ph = (AliCaloPhoton*)mixPHOS->At(iph);
        //if(!fPHOSClusterCuts->AcceptPhoton(ph)) continue;
        if(!CheckMinimumEnergy(ph)) continue;
        //if(!fIsMC && !ph->IsTOFOK()) continue;

        if(fIsPHOSTriggerAnalysis){
          if( fIsMC && fTRFM == AliAnalysisTaskPHOSPi0EtaToGammaGamma::kRFE && !fPHOSTriggerHelper->IsOnActiveTRUChannel(ph)) continue;//keep same TRU acceptance only in kRFE.
          if(!fIsMC && !ph->IsTrig()) continue;//it is meaningless to focus on photon without fired trigger in PHOS triggered data.
        }

        if(fForceActiveTRU && !fPHOSTriggerHelper->IsOnActiveTRUChannel(ph)) continue;//criterion fTRFM == kRFE is not needed.

        pT = ph->Pt();
        if(fUseCoreEnergy) pT = (ph->GetMomV2())->Pt();

        position[0] = ph->EMCx();
        position[1] = ph->EMCy();
        position[2] = ph->EMCz();

        TVector3 global(position);
        relId[0] = 0; relId[1] = 0; relId[2] = 0; relId[3] = 0;
        fPHOSGeo->GlobalPos2RelId(global,relId);
        Int_t module = relId[0];
        TVector3 locPos;
        fPHOSGeo->Global2Local(locPos,global,module);
        cpv = 999.; dx  = 999.; dz  = 999.;
        Int_t itr = supply->FindTrackMatching(module,&locPos,dx,dz,pttrack,charge);

        if(itr > 0){
          cpv = supply->TestCPV(dx,dz,pttrack,charge);
          //AliInfo(Form("dx = %3.2f cm , dz = %3.2f cm , distance in sigma = %3.2f",dx,dz,cpv));

          //FillHistogramTH2(fOutputContainer,"hMixRvsTrackPt",pttrack,cpv);
          //FillHistogramTH2(fOutputContainer,"hMixRvsClusterPt",pTcluster,cpv);

          pidPion = kFALSE;
          pidKaon = kFALSE;
          pidProton = kFALSE;
          pidElectron = kFALSE;
          isGlobal = kFALSE;

          AliVTrack *track = (AliVTrack*)fEvent->GetTrack(itr);
          if(fESDEvent){
            isGlobal = fESDtrackCutsGlobal->AcceptTrack(dynamic_cast<AliESDtrack*>(track));
          }//end of ESD
          else if(fAODEvent){
            isGlobal = dynamic_cast<AliAODTrack*>(track)->TestFilterMask(AliAODTrack::kTrkGlobalNoDCA);//standard cuts with very loose DCA cut
          }//end of AOD

          nsigmaElectron = fPIDResponse->NumberOfSigmasTPC(track,AliPID::kElectron);//(-2,3) is for real electron analysis.
          nsigmaPion     = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track,AliPID::kPion));
          nsigmaKaon     = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track,AliPID::kKaon));
          nsigmaProton   = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track,AliPID::kProton));

          if((TMath::Abs(nsigmaElectron) < nsigmaPion)     && (TMath::Abs(nsigmaElectron) < nsigmaKaon) && (TMath::Abs(nsigmaElectron) < nsigmaProton) && (-2 < nsigmaElectron && nsigmaElectron < 3)) pidElectron = kTRUE;
          if((nsigmaPion     < TMath::Abs(nsigmaElectron)) && (nsigmaPion     < nsigmaKaon) && (nsigmaPion     < nsigmaProton) && (nsigmaPion     < NsigmaCut)) pidPion = kTRUE;
          if((nsigmaKaon     < TMath::Abs(nsigmaElectron)) && (nsigmaKaon     < nsigmaPion) && (nsigmaKaon     < nsigmaProton) && (nsigmaKaon     < NsigmaCut)) pidKaon = kTRUE;
          if((nsigmaProton   < TMath::Abs(nsigmaElectron)) && (nsigmaProton   < nsigmaPion) && (nsigmaProton   < nsigmaKaon)   && (nsigmaProton   < NsigmaCut)) pidProton = kTRUE;

          if(pidElectron){
            FillHistogramTH2(fOutputContainer,"hMixRvsTrackPt_Electron",pttrack,cpv);
            FillHistogramTH2(fOutputContainer,"hMixRvsClusterPt_Electron",pT,cpv);
          }
          else if(pidPion){
            FillHistogramTH2(fOutputContainer,"hMixRvsTrackPt_Pion",pttrack,cpv);
            FillHistogramTH2(fOutputContainer,"hMixRvsClusterPt_Pion",pT,cpv);
          }
          else if(pidKaon){
            FillHistogramTH2(fOutputContainer,"hMixRvsTrackPt_Kaon",pttrack,cpv);
            FillHistogramTH2(fOutputContainer,"hMixRvsClusterPt_Kaon",pT,cpv);
          }
          else if(pidProton){
            if(charge > 0){
              FillHistogramTH2(fOutputContainer,"hMixRvsTrackPt_Proton",pttrack,cpv);
              FillHistogramTH2(fOutputContainer,"hMixRvsClusterPt_Proton",pT,cpv);
            }
            else{
              FillHistogramTH2(fOutputContainer,"hMixRvsTrackPt_AntiProton",pttrack,cpv);
              FillHistogramTH2(fOutputContainer,"hMixRvsClusterPt_AntiProton",pT,cpv);
            }
          }
        }//end of track matching

      }//end of cluster loop

    }//end of mixed event loop

  }//end of accessing PHOSTender

}
//_______________________________________________________________________________
//_______________________________________________________________________________

