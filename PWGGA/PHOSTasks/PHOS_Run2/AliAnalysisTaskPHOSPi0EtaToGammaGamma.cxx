#include "TChain.h"
#include "TTree.h"
#include "TObjArray.h"
#include "TClonesArray.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TProfile.h"
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
#include "AliGenCocktailEventHeader.h"


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
  fMinDistBCM(-1),
  fCollisionSystem(-1),
  fTOFEfficiency(0x0),
  fESDtrackCutsGlobal(0x0),
  fESDtrackCutsGlobalConstrained(0x0),
  fCentArrayPi0(0x0),
  fCentArrayK0S(0x0),

  fOutputContainer(0x0),
  fEvent(0x0),
  fESDEvent(0x0),
  fAODEvent(0x0),
  fMCArrayESD(0x0),
  fMCArrayAOD(0x0),
  fJJMCHandler(0x0),
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
  fIsFlowTask(kFALSE),
  fQnEstimator("QoverM"),
  fFlowQnVectorMgr(0x0),
  fEventPlane(-999.),
  fNHybridTrack(0),
  fIsPHOSTriggerAnalysis(kFALSE),
  fPHOSTriggerHelper(0x0),
  fPHOSTriggerHelperL0(0x0),
  fPHOSTriggerHelperL1H(0x0),
  fPHOSTriggerHelperL1M(0x0),
  fPHOSTriggerHelperL1L(0x0),
  fPIDResponse(0x0)
{
  // Constructor

  for(Int_t i=0;i<10;i++){
      fPHOSEvents[i] = 0;
  }

  for(Int_t i=0;i<3;i++){
    fVertex[i] = 0;
  }


  fTPCEPName[0] = "TPC";
  fTPCEPName[1] = "TPCNegEta";
  fTPCEPName[2] = "TPCPosEta";

  fV0EPName[0] = "VZERO";
  fV0EPName[1] = "VZEROA";
  fV0EPName[2] = "VZEROC";

  fPHOSTriggerHelperL0  = new AliPHOSTriggerHelper("L0");
  fPHOSTriggerHelperL1H = new AliPHOSTriggerHelper("L1H");
  fPHOSTriggerHelperL1M = new AliPHOSTriggerHelper("L1M");
  fPHOSTriggerHelperL1L = new AliPHOSTriggerHelper("L1L");

  for(Int_t i=0;i<10;i++){
    fAdditionalPi0PtWeight[i] = new TF1("fAdditionalPi0PtWeight","1.",0,100);
    fAdditionalK0SPtWeight[i] = new TF1("fAdditionalK0SPtWeight","1.",0,100);
  }

  fTOFEfficiency = new TF1("fTOFEfficiency","1.",0,100);

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
    if(fPHOSEvents[i]){
      delete fPHOSEvents[i];
      fPHOSEvents[i] = 0;
    }
  }


  if(fPHOSTriggerHelper){
    delete fPHOSTriggerHelper;
    fPHOSTriggerHelper  = 0x0;
  }

  delete fPHOSTriggerHelperL0;
  delete fPHOSTriggerHelperL1H;
  delete fPHOSTriggerHelperL1M;
  delete fPHOSTriggerHelperL1L;

  fPHOSTriggerHelperL0  = 0x0;
  fPHOSTriggerHelperL1H = 0x0;
  fPHOSTriggerHelperL1M = 0x0;
  fPHOSTriggerHelperL1L = 0x0;

  for(Int_t i=0;i<10;i++){
    if(fAdditionalPi0PtWeight[i]){
      delete fAdditionalPi0PtWeight[i];
      fAdditionalPi0PtWeight[i] = 0x0;
    }

    if(fAdditionalK0SPtWeight[i]){
      delete fAdditionalK0SPtWeight[i];
      fAdditionalK0SPtWeight[i] = 0x0;
    }
  }


  delete fTOFEfficiency;
  fTOFEfficiency = 0x0;

  if(fPHOSEventCuts){
    delete fPHOSEventCuts;
    fPHOSEventCuts = 0x0;
  }

  if(fPHOSClusterCuts){
    delete fPHOSClusterCuts;
    fPHOSClusterCuts = 0x0;
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
  fOutputContainer->Add(new TH1F("hVertexZ","VertexZ",1000,-50.,50.));
  fOutputContainer->Add(new TH1F("hVertexZSelectEvent","VertexZ SelectEvent",1000,-50.,50.));

  fOutputContainer->Add(new TH2F("hVertexXY","VertexXY",100,-0.5,0.5,100,-0.5,0.5));
  fOutputContainer->Add(new TH2F("hVertexXYSelectEvent","VertexXY SelectEvent",100,-0.5,0.5,100,-0.5,0.5));

  fOutputContainer->Add(new TH2F("hCentralityV0MvsNContibutor","Centrality V0M vs. Ncontributor",100,0.,100,101,-0.5,100.5));
  fOutputContainer->Add(new TH2F("hCentralityCL0vsNContibutor","Centrality CL0 vs. Ncontributor",100,0.,100,101,-0.5,100.5));
  fOutputContainer->Add(new TH2F("hCentralityCL1vsNContibutor","Centrality CL1 vs. Ncontributor",100,0.,100,101,-0.5,100.5));
  fOutputContainer->Add(new TH2F("hCentralityV0AvsNContibutor","Centrality V0A vs. Ncontributor",100,0.,100,101,-0.5,100.5));
  fOutputContainer->Add(new TH2F("hCentralityV0CvsNContibutor","Centrality V0C vs. Ncontributor",100,0.,100,101,-0.5,100.5));
  fOutputContainer->Add(new TH2F("hCentralityZNAvsNContibutor","Centrality ZNA vs. Ncontributor",100,0.,100,101,-0.5,100.5));
  fOutputContainer->Add(new TH2F("hCentralityZNCvsNContibutor","Centrality ZNC vs. Ncontributor",100,0.,100,101,-0.5,100.5));

  fOutputContainer->Add(new TH2F("hCentralityV0MvsCL0","Centrality V0M vs. CL0",100,0.,100,100,0.,100.));
  fOutputContainer->Add(new TH2F("hCentralityV0MvsCL1","Centrality V0M vs. CL1",100,0.,100,100,0.,100.));
  fOutputContainer->Add(new TH2F("hCentralityCL0vsCL1","Centrality CL0 vs. CL1",100,0.,100,100,0.,100.));
  fOutputContainer->Add(new TH2F("hCentralityV0AvsV0C","Centrality V0A vs. V0C",100,0.,100,100,0.,100.));
  fOutputContainer->Add(new TH2F("hCentralityZNAvsZNC","Centrality ZNA vs. ZNC",100,0.,100,100,0.,100.));

  fOutputContainer->Add(new TH2F(Form("hCentrality%svsDeltaEventPlane%s%s",fEstimator.Data(),fV0EPName[0].Data() ,fTPCEPName[1].Data()),Form("Centrality %s vs. cos(#Delta#Psi_{EP});centrality (%%);cos(#Delta#Psi_{EP})",fEstimator.Data()),100,0,100,200,-1,1));
  fOutputContainer->Add(new TH2F(Form("hCentrality%svsDeltaEventPlane%s%s",fEstimator.Data(),fV0EPName[0].Data() ,fTPCEPName[2].Data()),Form("Centrality %s vs. cos(#Delta#Psi_{EP});centrality (%%);cos(#Delta#Psi_{EP})",fEstimator.Data()),100,0,100,200,-1,1));
  fOutputContainer->Add(new TH2F(Form("hCentrality%svsDeltaEventPlane%s%s",fEstimator.Data(),fTPCEPName[1].Data(),fTPCEPName[2].Data()),Form("Centrality %s vs. cos(#Delta#Psi_{EP});centrality (%%);cos(#Delta#Psi_{EP})",fEstimator.Data()),100,0,100,200,-1,1));

  fOutputContainer->Add(new TH2F(Form("hCentrality%svsDeltaEventPlane%s%s",fEstimator.Data(),fTPCEPName[0].Data(),fV0EPName[1].Data()),Form("Centrality %s vs. cos(#Delta#Psi_{EP});centrality (%%);cos(#Delta#Psi_{EP})",fEstimator.Data()),100,0,100,200,-1,1));
  fOutputContainer->Add(new TH2F(Form("hCentrality%svsDeltaEventPlane%s%s",fEstimator.Data(),fTPCEPName[0].Data(),fV0EPName[2].Data()),Form("Centrality %s vs. cos(#Delta#Psi_{EP});centrality (%%);cos(#Delta#Psi_{EP})",fEstimator.Data()),100,0,100,200,-1,1));
  fOutputContainer->Add(new TH2F(Form("hCentrality%svsDeltaEventPlane%s%s",fEstimator.Data(),fV0EPName[1].Data() ,fV0EPName[2].Data()),Form("Centrality %s vs. cos(#Delta#Psi_{EP});centrality (%%);cos(#Delta#Psi_{EP})",fEstimator.Data()),100,0,100,200,-1,1));

  for(Int_t i=0;i<3;i++){
    fOutputContainer->Add(new TH2F(Form("hCentrality%svsEventPlane%s%s",fEstimator.Data(),fTPCEPName[i].Data(),fQnEstimator.Data()),Form("Centrality %s vs. EP %s %s;centrality (%%);#Psi_{EP}",fEstimator.Data(),fTPCEPName[i].Data(),fQnEstimator.Data()),100,0,100,320,0,3.2));
  }
  for(Int_t i=0;i<3;i++){
    fOutputContainer->Add(new TH2F(Form("hCentrality%svsEventPlane%s%s",fEstimator.Data(),fV0EPName[i].Data(),fQnEstimator.Data()),Form("Centrality %s vs. EP %s %s;centrality (%%);#Psi_{EP}",fEstimator.Data(),fV0EPName[i].Data(),fQnEstimator.Data()),100,0,100,320,0,3.2));
  }


  fOutputContainer->Add(new TH2F(Form("hCentrality%svsPHOSClusterMultiplicityMC" ,fEstimator.Data()),Form("Centrality %s vs. Cluster Multiplicity;centrality (%%);Ncluster"                  ,fEstimator.Data()),100,0,100,201,-0.5,200.5));
  fOutputContainer->Add(new TH2F(Form("hCentrality%svsPHOSClusterMultiplicity"   ,fEstimator.Data()),Form("Centrality %s vs. Cluster Multiplicity;centrality (%%);Ncluster"                  ,fEstimator.Data()),100,0,100,201,-0.5,200.5));
  fOutputContainer->Add(new TH2F(Form("hCentrality%svsPHOSClusterMultiplicityTOF",fEstimator.Data()),Form("Centrality %s vs. Cluster Multiplicity with TOF cut;centrality (%%);Ncluster"     ,fEstimator.Data()),100,0,100,201,-0.5,200.5));
  fOutputContainer->Add(new TH2F(Form("hCentrality%svsPHOSClusterMultiplicityNoPID",fEstimator.Data()),Form("Centrality %s vs. Cluster Multiplicity;centrality (%%);Ncluster"                ,fEstimator.Data()),100,0,100,201,-0.5,200.5));
  fOutputContainer->Add(new TH2F(Form("hCentrality%svsPHOSClusterMultiplicityNoPIDTOF",fEstimator.Data()),Form("Centrality %s vs. Cluster Multiplicity with TOF cut;centrality (%%);Ncluster",fEstimator.Data()),100,0,100,201,-0.5,200.5));

  fOutputContainer->Add(new TH2F(Form("hCentrality%svsTrackMultiplicity",fEstimator.Data()),Form("Centrality %s vs. track Multiplicity",fEstimator.Data()),100,0,100,5000,0,5000));
  fOutputContainer->Add(new TH2F("hPHOSClusterMultiplicityvsTrackMultiplicity"   ,"cluster multiplicity vs. track multiplicity"         ,500,0,5000,201,-0.5,200.5));
  fOutputContainer->Add(new TH2F("hPHOSClusterMultiplicityTOFvsTrackMultiplicity","cluster multiplicity with TOF vs. track multiplicity",500,0,5000,201,-0.5,200.5));

  //track QA histograms
  const Int_t Ntype=3;
  const TString tracktype[Ntype] = {"Hybrid","Global","Complementary"};

  for(Int_t itype=0;itype<Ntype;itype++){
    fOutputContainer->Add(new TH1F(Form("h%sTrackMult",tracktype[itype].Data())  ,Form("Number of %s track",tracktype[itype].Data()),500,0,5000));
    fOutputContainer->Add(new TH1F(Form("h%sTrackPt",tracktype[itype].Data())    ,Form("%s track p_{T}",tracktype[itype].Data()),100,0,100));
    fOutputContainer->Add(new TH2F(Form("h%sTrackEtaPhi",tracktype[itype].Data()),Form("%s track #eta vs. #phi;#phi;#eta",tracktype[itype].Data()),63,0,6.3,20,-1,1));
  }
  fOutputContainer->Add(new TH2F("hTrackTPCdEdx","TPC dE/dx vs. track momentum;p^{track} (GeV/c);dE/dx (a.u.)",40,0,20,200,0,200));

  //cell QA histograms
  const Int_t Nmod=5;
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

  for(Int_t imod=1;imod<Nmod;imod++) fOutputContainer->Add(new TH2F(Form("hClusterEvsNM%d",imod),Form("Cluster E vs N_{cell} M%d;E (GeV);N_{cell}",imod),500,0,50,100,0.5,100.5));
  for(Int_t imod=1;imod<Nmod;imod++) fOutputContainer->Add(new TH2F(Form("hClusterEvsTM%d",imod),Form("Cluster E vs TOF M%d;E (GeV);TOF (ns)",imod)     ,500,0,50, 1000,-500,500));
  for(Int_t imod=1;imod<Nmod;imod++) fOutputContainer->Add(new TH2F(Form("hClusterEvsM02M%d",imod),Form("Cluster E vs M02 M%d;E (GeV);M02 (cm)",imod),500,0,50,100,0,10));
  for(Int_t imod=1;imod<Nmod;imod++) fOutputContainer->Add(new TH2F(Form("hClusterNvsM02M%d",imod),Form("Cluster N vs M02 M%d;N_{cell};M02 (cm)",imod),100,0.5,100.5,100,0,10));
  for(Int_t imod=1;imod<Nmod;imod++) fOutputContainer->Add(new TH2F(Form("hFullDispvsFullEM%d",imod),Form("full dispersion vs full E M%d;E (GeV);dipersion (#sigma)",imod),100,0,50,100,0,10));
  for(Int_t imod=1;imod<Nmod;imod++) fOutputContainer->Add(new TH2F(Form("hCoreDispvsCoreEM%d",imod),Form("core dispersion vs core E M%d;E (GeV);dipersion (#sigma)",imod),100,0,50,100,0,10));
  for(Int_t imod=1;imod<Nmod;imod++) fOutputContainer->Add(new TH2F(Form("hFullDispvsCoreEM%d",imod),Form("full dispersion vs full E M%d;E (GeV);dipersion (#sigma)",imod),100,0,50,100,0,10));
  for(Int_t imod=1;imod<Nmod;imod++) fOutputContainer->Add(new TH2F(Form("hCoreDispvsFullEM%d",imod),Form("core dispersion vs core E M%d;E (GeV);dipersion (#sigma)",imod),100,0,50,100,0,10));
  for(Int_t imod=1;imod<Nmod;imod++) fOutputContainer->Add(new TH2F(Form("hRvsTrackPtM%d",imod),Form("r vs track pT M%d;p_{T}^{track} (GeV/c);cpv (#sigma)",imod)         ,100,0,50,100,0,10));

  for(Int_t imod=1;imod<Nmod;imod++) fOutputContainer->Add(new TH3F(Form("hdZvsZvsTrackPt_M%d",imod)        ,"dZ vs. Z;Z (cm);dZ (cm);p_{T}^{track} (GeV/c)"           ,160,-80,80,80,-20,20,40,0,20));
  for(Int_t imod=1;imod<Nmod;imod++) fOutputContainer->Add(new TH3F(Form("hdXvsXvsTrackPt_plus_M%d",imod)   ,"dX vs. X positive;X (cm);dX (cm);p_{T}^{track +} (GeV/c)",160,-80,80,80,-20,20,40,0,20));
  for(Int_t imod=1;imod<Nmod;imod++) fOutputContainer->Add(new TH3F(Form("hdXvsXvsTrackPt_minus_M%d",imod)  ,"dX vs. X negative;X (cm);dX (cm);p_{T}^{track -} (GeV/c)",160,-80,80,80,-20,20,40,0,20));
  for(Int_t imod=1;imod<Nmod;imod++) fOutputContainer->Add(new TH3F(Form("hdZvsZvsTrackPtElectron_M%d",imod),"dZ vs. Z of e^{#pm};Z (cm);dZ (cm);p_{T}^{track} (GeV/c)",160,-80,80,80,-20,20,40,0,20));//for radial displacement

  for(Int_t imod=1;imod<Nmod;imod++){
    fOutputContainer->Add(new TH2F(Form("hEpRatiovsEnergy_M%d_Electron",imod) ,Form("E/p ratio vs. E_{cluster} M%d;E/p;E_{cluster} (GeV)",imod)    ,100,0,2,100,0,50));
    fOutputContainer->Add(new TH2F(Form("hEpRatiovsEnergy_M%d_Others",imod)   ,Form("E/p ratio vs. E_{cluster} M%d;E/p;E_{cluster} (GeV)",imod)    ,100,0,2,100,0,50));
    fOutputContainer->Add(new TH2F(Form("hEpRatiovsTrackPt_M%d_Electron",imod),Form("E/p ratio vs. E_{cluster} M%d;E/p;p_{T}^{track} (GeV/c)",imod),100,0,2,100,0,50));
    fOutputContainer->Add(new TH2F(Form("hEpRatiovsTrackPt_M%d_Others",imod)  ,Form("E/p ratio vs. E_{cluster} M%d;E/p;p_{T}^{track} (GeV/c)",imod),100,0,2,100,0,50));
  }

  fOutputContainer->Add(new TH2F("hEpRatiovsNsigmaElectronTPC","E/p ratio vs. N_{#sigma}^{e};E/p;n#sigma^{e}",100,0,2,100,-5,5));
  fOutputContainer->Add(new TH2F("hTPCdEdx_Electron","TPC dEdx vs. electron momentum;p^{track} (GeV/c);dE/dx (a.u.)"    ,40,0,20,200,0,200));
  fOutputContainer->Add(new TH2F("hTPCdEdx_Others"  ,"TPC dEdx vs. non-electron momentum;p^{track} (GeV/c);dE/dx (a.u.)",40,0,20,200,0,200));

  fOutputContainer->Add(new TH2F("hClusterEtaPhi","Cluster eta vs. phi;#phi;#eta",63,0,6.3,200,-1,1));
  fOutputContainer->Add(new TH2F("hEnergyvsDistanceToBadChannel","distance to closest bad channel;E (GeV);distance in cell",100,0,50,20,0,10));
  fOutputContainer->Add(new TH2F("hAsymvsMgg","asymmetry vs. M_{#gamma#gamma};asymmetry;M_{#gamma#gamma} (GeV/c^{2})",10,0,1,60,0,0.24));
  fOutputContainer->Add(new TH2F("hAsymvsPt" ,"asymmetry vs. p_{T}^{#gamma#gamma};asymmetry;p_{T} (GeV/c)",10,0,1,100,0,50));

  //<- histograms for QA
  //histograms for physics analysis ->

  const Int_t NpTgg = 101;
  Double_t pTgg[NpTgg]={};

  for(Int_t i=0;i<50;i++)     pTgg[i] = 0.1 * i;            //every 0.1 GeV/c, up to 5 GeV/c
  for(Int_t i=50;i<60;i++)    pTgg[i] = 0.5 * (i-50) + 5.0; //every 0.5 GeV/c, up to 10 GeV/c
  for(Int_t i=60;i<NpTgg;i++) pTgg[i] = 1.0 * (i-60) + 10.0;//every 1.0 GeV/c, up to 50 GeV/c

  const Int_t Nm = 241;
  Double_t mgg[Nm] = {};
  for(Int_t i=0;i<Nm;i++) mgg[i] = 0.004*i;

  const Double_t EventPlane[] = {0,15,30,45,60,75,90,105,120,135,150,165,180};
  const Int_t Nep = sizeof(EventPlane)/sizeof(EventPlane[0]);

  //photon pT histograms
  TH2F *h2PhotonPt = new TH2F("hPhotonPt","Photon Pt;p_{T} (GeV/c);#Delta#phi (deg)",NpTgg-1,pTgg,Nep-1,EventPlane);
  h2PhotonPt->Sumw2();
  fOutputContainer->Add(h2PhotonPt);

  TH2F *h2PhotonPt_TOF = new TH2F("hPhotonPt_TOF","Photon Pt with TOF cut;p_{T} (GeV/c);#Delta#phi (deg)",NpTgg-1,pTgg,Nep-1,EventPlane);
  h2PhotonPt_TOF->Sumw2();
  fOutputContainer->Add(h2PhotonPt_TOF);

  const Int_t Ndimg = 4;
  const Int_t    Nbing[Ndimg] = {10,10,10,500};//cpv1,disp1,distBC,pt
  const Double_t xming[Ndimg] = { 0, 0, 0,  0};//cpv1,disp1,distBC,pt
  const Double_t xmaxg[Ndimg] = { 5, 5, 5, 50};//cpv1,disp1,distBC,pt

  THnSparse *hsPhoton = new THnSparseF("hSparsePhoton","Photon;CPV (#sigma);Disp (#sigma);distance to BC (cell);p_{T} (GeV/c);",Ndimg,Nbing,xming,xmaxg);
  hsPhoton->Sumw2();
  fOutputContainer->Add(hsPhoton);

  const TString Asym[] = {"","_asym08"};
  const Int_t Nasym = sizeof(Asym)/sizeof(Asym[0]);

  //Mgg vs. pT histogram
  for(Int_t iasym=0;iasym<Nasym;iasym++){
    TH3F *h3 = new TH3F(Form("hMgg%s",Asym[iasym].Data()),"Invariant Mass with 2#gamma;M_{#gamma#gamma} (GeV/c^{2});p_{T} (GeV/c);#Delta#phi (deg.)",Nm-1,mgg,NpTgg-1,pTgg,Nep-1,EventPlane);
    h3->Sumw2();
    fOutputContainer->Add(h3);

    TH3F *h3_TOF = new TH3F(Form("hMgg_TOF%s",Asym[iasym].Data()),"Invariant Mass of 2#gamma with TOF cut;M_{#gamma#gamma} (GeV/c^{2});p_{T} (GeV/c);#Delta#phi (deg.)",Nm-1,mgg,NpTgg-1,pTgg,Nep-1,EventPlane);
    h3_TOF->Sumw2();
    fOutputContainer->Add(h3_TOF);


    TH3F *h3Mix = new TH3F(Form("hMixMgg%s",Asym[iasym].Data()),"Mix Invariant Mass with 2#gamma;M_{#gamma#gamma} (GeV/c^{2});p_{T} (GeV/c);#Delta#phi (deg.)",Nm-1,mgg,NpTgg-1,pTgg,Nep-1,EventPlane);
    h3Mix->Sumw2();
    fOutputContainer->Add(h3Mix);

    TH3F *h3Mix_TOF = new TH3F(Form("hMixMgg_TOF%s",Asym[iasym].Data()),"Mix Invariant Mass of 2#gamma with TOF cut;M_{#gamma#gamma} (GeV/c^{2});p_{T} (GeV/c);#Delta#phi (deg.)",Nm-1,mgg,NpTgg-1,pTgg,Nep-1,EventPlane);
    h3Mix_TOF->Sumw2();
    fOutputContainer->Add(h3Mix_TOF);
  }//end of asymmetry loop


  //<- histograms for physics analysis


  const Int_t NpTggModule = 71;
  Double_t pTggModule[NpTggModule]={};
  for(Int_t i=0;i<50;i++)     pTggModule[i] = 0.1 * i;            //every 0.1 GeV/c, up to 5 GeV/c
  for(Int_t i=50;i<60;i++)    pTggModule[i] = 0.5 * (i-50) + 5.0; //every 0.5 GeV/c, up to 10 GeV/c
  for(Int_t i=60;i<NpTggModule;i++) pTggModule[i] = 1.0 * (i-60) + 10.0;//every 1.0 GeV/c, up to 20 GeV/c

  for(Int_t imod1=1;imod1<Nmod;imod1++){
    for(Int_t imod2=imod1;imod2<Nmod;imod2++){
      if(imod2 - imod1 > 1) continue;

      TH2F *h2 = new TH2F(Form("hMgg_M%d%d",imod1,imod2),"M_{#gamma#gamma} vs p_{T}",60,0,0.24,NpTggModule-1,pTggModule);
      h2->Sumw2();
      fOutputContainer->Add(h2);
    }
  }
  for(Int_t imod1=1;imod1<Nmod;imod1++){
    for(Int_t imod2=imod1;imod2<Nmod;imod2++){
      if(imod2 - imod1 > 1) continue;

      TH2F *h2 = new TH2F(Form("hMixMgg_M%d%d",imod1,imod2),"M_{#gamma#gamma} vs p_{T}",60,0,0.24,NpTggModule-1,pTggModule);
      h2->Sumw2();
      fOutputContainer->Add(h2);
    }
  }

  for(Int_t imod1=1;imod1<Nmod;imod1++){
    for(Int_t imod2=imod1;imod2<Nmod;imod2++){
      if(imod2 - imod1 > 1) continue;

      TH2F *h2 = new TH2F(Form("hMgg_M%d%d_TOF",imod1,imod2),"M_{#gamma#gamma} vs p_{T}",60,0,0.24,NpTggModule-1,pTggModule);
      h2->Sumw2();
      fOutputContainer->Add(h2);
    }
  }

  for(Int_t imod1=1;imod1<Nmod;imod1++){
    for(Int_t imod2=imod1;imod2<Nmod;imod2++){
      if(imod2 - imod1 > 1) continue;

      TH2F *h2 = new TH2F(Form("hMixMgg_M%d%d_TOF",imod1,imod2),"M_{#gamma#gamma} vs p_{T}",60,0,0.24,NpTggModule-1,pTggModule);
      h2->Sumw2();
      fOutputContainer->Add(h2);
    }
  }

  //for PID cut efficiency
  TH2F *h2_p_PID = new TH2F("hMgg_Probe_PID"              ,"Probe #gamma PID;M_{#gamma#gamma} (GeV/c^{2});E_{#gamma} (GeV)"            ,60,0,0.24,NpTgg-1,pTgg);
  h2_p_PID->Sumw2();
  fOutputContainer->Add(h2_p_PID);

  TH2F *h2_pp_PID = new TH2F("hMgg_PassingProbe_PID"      ,"Passing Probe #gamma PID;M_{#gamma#gamma} (GeV/c^{2});E_{#gamma} (GeV)"    ,60,0,0.24,NpTgg-1,pTgg);
  h2_pp_PID->Sumw2();
  fOutputContainer->Add(h2_pp_PID);

  TH2F *h2mix_p_PID = new TH2F("hMixMgg_Probe_PID"        ,"Mix Probe #gamma PID;M_{#gamma#gamma} (GeV/c^{2});E_{#gamma} (GeV)"        ,60,0,0.24,NpTgg-1,pTgg);
  h2mix_p_PID->Sumw2();
  fOutputContainer->Add(h2mix_p_PID);

  TH2F *h2mix_pp_PID = new TH2F("hMixMgg_PassingProbe_PID","Mix Passing Probe #gamma PID;M_{#gamma#gamma} (GeV/c^{2});E_{#gamma} (GeV)",60,0,0.24,NpTgg-1,pTgg);
  h2mix_pp_PID->Sumw2();
  fOutputContainer->Add(h2mix_pp_PID);

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

  //for nonlinearity
  TH2F *h2MggvsE = new TH2F("hMggvsE_asym01","Invariant Mass with 2#gamma;M_{#gamma#gamma} (GeV/c^{2});E_{cluster} (GeV)",60,0,0.24,NpTgg-1,pTgg);
  h2MggvsE->Sumw2();
  fOutputContainer->Add(h2MggvsE);

  TH2F *h2MixMggvsE = new TH2F("hMixMggvsE_asym01","Mix Invariant Mass with 2#gamma;M_{#gamma#gamma} (GeV/c^{2});E_{cluster} (GeV)",60,0,0.24,NpTgg-1,pTgg);
  h2MixMggvsE->Sumw2();
  fOutputContainer->Add(h2MixMggvsE);

  if(fIsPHOSTriggerAnalysis){

    //for Trigger QA
    fOutputContainer->Add(new TH1F("hNFiredTRUChannel","Fired TRU channel per event",101,-0.5,100.5));
    for(Int_t imod=1;imod<Nmod;imod++) fOutputContainer->Add(new TH2F(Form("hFiredTRUChannelM%d",imod),Form("Fired TRU channel M%d",imod),64,0.5,64.5, 56,0.5,56.5));
    for(Int_t imod=1;imod<Nmod;imod++) fOutputContainer->Add(new TH3F(Form("hMatchedFiredTRUChannelM%d",imod),Form("Fired TRU channel M%d",imod),64,0.5,64.5, 56,0.5,56.5,100,0,50));

    const Int_t Ntru=8;
    for(Int_t imod=1;imod<Nmod;imod++){
      fOutputContainer->Add(new TH1F(Form("hClusterEnergyM%d",imod)       ,Form("cluster energy M%d",imod),500,0,50));
      fOutputContainer->Add(new TH1F(Form("hMatchedClusterEnergyM%d",imod),Form("cluster energy M%d",imod),500,0,50));

      fOutputContainer->Add(new TH3F(Form("hClusterHitM%d",imod)       ,Form("cluster hit M%d",imod)         ,64,0.5,64.5,56,0.5,56.5,100,0,50));
      fOutputContainer->Add(new TH3F(Form("hMatchedClusterHitM%d",imod),Form("matched cluster hit M%d",imod) ,64,0.5,64.5,56,0.5,56.5,100,0,50));

      for(Int_t itru=1;itru<=Ntru;itru++){
        fOutputContainer->Add(new TH1F(Form("hClusterEnergyM%dTRU%d",imod,itru)       ,Form("cluster energy M%d",imod),500,0,50));
        fOutputContainer->Add(new TH1F(Form("hMatchedClusterEnergyM%dTRU%d",imod,itru),Form("cluster energy M%d",imod),500,0,50));
        fOutputContainer->Add(new TH3F(Form("hDistanceXZM%dTRU%d",imod,itru)        ,Form("distance TRUch-clu M%d TRU%d",imod,itru),21,-10.5,+10.5,21,-10.5,+10.5,100,0,50));
        fOutputContainer->Add(new TH3F(Form("hMatchedDistanceXZM%dTRU%d",imod,itru) ,Form("distance TRUch-clu M%d TRU%d",imod,itru),21,-10.5,+10.5,21,-10.5,+10.5,100,0,50));
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

    TH2F *h2_p_Trg = new TH2F("hMgg_Probe_Trg"              ,"Probe #gamma Trg;M_{#gamma#gamma} (GeV/c^{2});E_{#gamma} (GeV)"            ,60,0,0.24,NpTgg-1,pTgg);
    h2_p_Trg->Sumw2();                                                                                                                    
    fOutputContainer->Add(h2_p_Trg);                                                                                                      

    TH2F *h2_pp_Trg = new TH2F("hMgg_PassingProbe_Trg"      ,"Passing Probe #gamma Trg;M_{#gamma#gamma} (GeV/c^{2});E_{#gamma} (GeV)"    ,60,0,0.24,NpTgg-1,pTgg);
    h2_pp_Trg->Sumw2();                                                                                                                   
    fOutputContainer->Add(h2_pp_Trg);                                                                                                     

    TH2F *h2mix_p_Trg = new TH2F("hMixMgg_Probe_Trg"       ,"Mix Probe #gamma Trg;M_{#gamma#gamma} (GeV/c^{2});E_{#gamma} (GeV)"        ,60,0,0.24,NpTgg-1,pTgg);
    h2mix_p_Trg->Sumw2();                                                                                                                 
    fOutputContainer->Add(h2mix_p_Trg);                                                                                                   

    TH2F *h2mix_pp_Trg = new TH2F("hMixMgg_PassingProbe_Trg","Mix Passing Probe #gamma Trg;M_{#gamma#gamma} (GeV/c^{2});E_{#gamma} (GeV)",60,0,0.24,NpTgg-1,pTgg);
    h2mix_pp_Trg->Sumw2();
    fOutputContainer->Add(h2mix_pp_Trg);

    for(Int_t imod=1;imod<Nmod;imod++){
      for(Int_t itru=1;itru<=Ntru;itru++){
        TH2F *h2 = new TH2F(Form("hMgg_Probe_Trg_M%d_TRU%d"          ,imod,itru),"Probe #gamma Trg;M_{#gamma#gamma} (GeV/c^{2});E_{#gamma} (GeV)"            ,60,0,0.24,NpTgg-1,pTgg);
        h2->Sumw2();
        fOutputContainer->Add(h2);
      }
    }

    for(Int_t imod=1;imod<Nmod;imod++){
      for(Int_t itru=1;itru<=Ntru;itru++){
        TH2F *h2 = new TH2F(Form("hMgg_PassingProbe_Trg_M%d_TRU%d"   ,imod,itru),"Passing Probe #gamma Trg;M_{#gamma#gamma} (GeV/c^{2});E_{#gamma} (GeV)"    ,60,0,0.24,NpTgg-1,pTgg);
        h2->Sumw2();
        fOutputContainer->Add(h2);
      }
    }

    for(Int_t imod=1;imod<Nmod;imod++){
      for(Int_t itru=1;itru<=Ntru;itru++){
        TH2F *h2 = new TH2F(Form("hMixMgg_Probe_Trg_M%d_TRU%d"       ,imod,itru),"Mix Probe #gamma Trg;M_{#gamma#gamma} (GeV/c^{2});E_{#gamma} (GeV)"        ,60,0,0.24,NpTgg-1,pTgg);
        h2->Sumw2();
        fOutputContainer->Add(h2);
      }
    }

    for(Int_t imod=1;imod<Nmod;imod++){
      for(Int_t itru=1;itru<=Ntru;itru++){
        TH2F *h2 = new TH2F(Form("hMixMgg_PassingProbe_Trg_M%d_TRU%d",imod,itru),"Mix Passing Probe #gamma Trg;M_{#gamma#gamma} (GeV/c^{2});E_{#gamma} (GeV)",60,0,0.24,NpTgg-1,pTgg);
        h2->Sumw2();
        fOutputContainer->Add(h2);
      }
    }

  }//fIsPHOSTriggerAnalysis ends

  if(fIsMC){
    const TString parname[] = {"Pi0","Eta","Gamma","ChargedPion","ChargedKaon","K0S","K0L","Lambda0","Sigma0"};
    const Int_t Npar = sizeof(parname)/sizeof(parname[0]);

    for(Int_t ipar=0;ipar<Npar;ipar++){
      TH1F *h1Pt = new TH1F(Form("hGen%sPt",parname[ipar].Data()),Form("generated %s pT;p_{T} (GeV/c)",parname[ipar].Data()),NpTgg-1,pTgg);
      h1Pt->Sumw2();
      fOutputContainer->Add(h1Pt);

      TH2F *h2EtaPhi = new TH2F(Form("hGen%sEtaPhi",parname[ipar].Data()),Form("generated %s eta vs phi;#eta;#phi (rad)",parname[ipar].Data()),63,0,6.3,200,-1,1);
      h2EtaPhi->Sumw2();
      fOutputContainer->Add(h2EtaPhi);

      TH2F *h2EtaPt = new TH2F(Form("hGen%sEtaPt",parname[ipar].Data()),Form("generated %s eta vs pT;#eta;p_{T} (GeV/c)",parname[ipar].Data()),200,-1,1,NpTgg-1,pTgg);
      h2EtaPt->Sumw2();
      fOutputContainer->Add(h2EtaPt);

    }//end of particle loop

    TH1F *h1TrueGamma = new TH1F("hPurityGamma","pT of true photon in clusters for purity;p_{T} (GeV/c)",NpTgg-1,pTgg);
    h1TrueGamma->Sumw2();
    fOutputContainer->Add(h1TrueGamma);

    for(Int_t iasym=0;iasym<Nasym;iasym++){
      TH2F *h2 = new TH2F(Form("hMggFromK0S%s",Asym[iasym].Data()),"M_{#gamma#gamma} from K_{S}^{0};M_{#gamma#gamma} (GeV/c^{2});p_{T} (GeV/c)",240,0,0.96,NpTgg-1,pTgg);
      h2->Sumw2();
      fOutputContainer->Add(h2);

    }//end of asym loop

    for(Int_t iasym=0;iasym<Nasym;iasym++){
      TH2F *h2 = new TH2F(Form("hMggFromLambda0%s",Asym[iasym].Data()),"M_{#gamma#gamma} from K_{S}^{0};M_{#gamma#gamma} (GeV/c^{2});p_{T} (GeV/c);",240,0,0.96,NpTgg-1,pTgg);
      h2->Sumw2();
      fOutputContainer->Add(h2);

    }//end of asym loop

    //for JJMC
    if(fIsJJMC){
      fOutputContainer->Add(new TH1F("hPtHard","pT hard in GeV/c",1000,0,1000));
      fOutputContainer->Add(new TH1F("hNTrial","nTrial",20,0.5,20.5));
      fOutputContainer->Add(new TProfile("hProfCrossSection","inelastic cross section",20,0.5,20.5));

      TH1F *hNMerged = new TH1F("hNMerged","N merged",20,0.5,20.5);
      //hNMerged->Fill(fPtHardBin);
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

  FillHistogramTH1(fOutputContainer,"hEventSummary",1);//all

  const AliVVertex *vVertex = fEvent->GetPrimaryVertex();
  fVertex[0] = vVertex->GetX();
  fVertex[1] = vVertex->GetY();
  fVertex[2] = vVertex->GetZ();
  FillHistogramTH1(fOutputContainer,"hVertexZ" ,fVertex[2]);
  FillHistogramTH2(fOutputContainer,"hVertexXY",fVertex[0],fVertex[1]);
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
  Float_t fCentralityV0M = -1.;
  Float_t fCentralityCL0 = -1.;
  Float_t fCentralityCL1 = -1.;
  Float_t fCentralityV0A = -1.;
  Float_t fCentralityV0C = -1.;
  Float_t fCentralityZNA = -1.;
  Float_t fCentralityZNC = -1.;

  //const Int_t system = GetCollisionSystem();
  if(fCollisionSystem==0){
    //not centrality, but track multiplicity
    //track QA to extrack number of hybrid track for centrality binning

  }
  else if(fCollisionSystem==1){//for PbPb
    //Get Centrality
    fMultSelection = (AliMultSelection*)fEvent->FindListObject("MultSelection");
    if(!fMultSelection){
      //If you get this warning (and fCentralityV0M 300) please check that the AliMultSelectionTask actually ran (before your task)
      AliWarning("AliMultSelection object not found!");
      return;
    }
    else{
      fCentralityV0M = fMultSelection->GetMultiplicityPercentile("V0M");
      fCentralityCL0 = fMultSelection->GetMultiplicityPercentile("CL0");
      fCentralityCL1 = fMultSelection->GetMultiplicityPercentile("CL1");
      fCentralityV0A = fMultSelection->GetMultiplicityPercentile("V0A");
      fCentralityV0C = fMultSelection->GetMultiplicityPercentile("V0C");
      fCentralityZNA = fMultSelection->GetMultiplicityPercentile("ZNA");
      fCentralityZNC = fMultSelection->GetMultiplicityPercentile("ZNC");
      fCentralityMain = fMultSelection->GetMultiplicityPercentile(fEstimator);
    }


    if(fCentralityMain < fCentralityMin || fCentralityMax < fCentralityMain){
      AliInfo(Form("Reject this event because centrality %s %f %% is out of the configuration of this task.", fEstimator.Data(),fCentralityMain));
      return;
    }

    AliInfo(Form("fCentralityMain estimated by %s = %f %%,Zvtx = %f cm, fZvtx = %d.",fEstimator.Data(),fCentralityMain,fVertex[2],fZvtx));

  }
  else if(fCollisionSystem==2){//for pPb,Pbp to be implemented
    //Get Centrality
    fMultSelection = (AliMultSelection*)fEvent->FindListObject("MultSelection");
    if(!fMultSelection){
      //If you get this warning (and fCentralityV0M 300) please check that the AliMultSelectionTask actually ran (before your task)
      AliWarning("AliMultSelection object not found!");
      return;
    }
    else{
      fCentralityV0M = fMultSelection->GetMultiplicityPercentile("V0M");
      fCentralityCL0 = fMultSelection->GetMultiplicityPercentile("CL0");
      fCentralityCL1 = fMultSelection->GetMultiplicityPercentile("CL1");
      fCentralityV0A = fMultSelection->GetMultiplicityPercentile("V0A");
      fCentralityV0C = fMultSelection->GetMultiplicityPercentile("V0C");
      fCentralityZNA = fMultSelection->GetMultiplicityPercentile("ZNA");
      fCentralityZNC = fMultSelection->GetMultiplicityPercentile("ZNC");
      fCentralityMain = fMultSelection->GetMultiplicityPercentile(fEstimator);
    }

    if(fCentralityMain < fCentralityMin || fCentralityMax < fCentralityMain){
      AliInfo(Form("Reject this event because centrality %s %f %% is out of the configuration of this task.", fEstimator.Data(),fCentralityMain));
      return;
    }

    AliInfo(Form("fCentralityMain estimated by %s = %f %%, Zvtx = %f cm, fZvtx = %d.",fEstimator.Data(),fCentralityMain,fVertex[2],fZvtx));

  }
  else{
    AliInfo("collision system shuold be 0-pp , 1-PbPb, 2-pPb/Pbp.");
    return;
  }

  FillHistogramTH2(fOutputContainer,"hCentralityV0MvsNContibutor",fCentralityV0M,Ncontributor);
  FillHistogramTH2(fOutputContainer,"hCentralityCL0vsNContibutor",fCentralityCL0,Ncontributor);
  FillHistogramTH2(fOutputContainer,"hCentralityCL1vsNContibutor",fCentralityCL1,Ncontributor);
  FillHistogramTH2(fOutputContainer,"hCentralityV0AvsNContibutor",fCentralityV0A,Ncontributor);
  FillHistogramTH2(fOutputContainer,"hCentralityV0CvsNContibutor",fCentralityV0C,Ncontributor);
  FillHistogramTH2(fOutputContainer,"hCentralityZNAvsNContibutor",fCentralityZNA,Ncontributor);
  FillHistogramTH2(fOutputContainer,"hCentralityZNCvsNContibutor",fCentralityZNC,Ncontributor);

  UInt_t fSelectMask = fInputHandler->IsEventSelected();
  Bool_t isINT7selected = fSelectMask & AliVEvent::kINT7;

  if(!fIsPHOSTriggerAnalysis && !isINT7selected){
    AliInfo("INT7 Event is rejected by IsEventSelected()");
    return;
  }

  Bool_t isPHI7selected = fSelectMask & AliVEvent::kPHI7;
  if(fIsPHOSTriggerAnalysis && !isPHI7selected){
    AliInfo("PHI7 Event is rejected by IsEventSelected()");
    return;
  }

  //event selection
  if(!(fPHOSEventCuts->AcceptEvent(fEvent))){
    AliInfo("event is rejected.");
    return;
  }

  //fill fired trigger statistics right after basic event selection, but before PHI7 event selection.
  Bool_t Is0PH0fired = fEvent->GetHeader()->GetL0TriggerInputs() & 1 << (9 -1);//trigger input -1
  Bool_t Is1PHHfired = fEvent->GetHeader()->GetL1TriggerInputs() & 1 << (7 -1);//trigger input -1
  Bool_t Is1PHMfired = fEvent->GetHeader()->GetL1TriggerInputs() & 1 << (6 -1);//trigger input -1
  Bool_t Is1PHLfired = fEvent->GetHeader()->GetL1TriggerInputs() & 1 << (5 -1);//trigger input -1

  if(Is0PH0fired) FillHistogramTH1(fOutputContainer,"hEventSummary",3);//0PH0
  if(Is1PHLfired) FillHistogramTH1(fOutputContainer,"hEventSummary",4);//1PHL
  if(Is1PHMfired) FillHistogramTH1(fOutputContainer,"hEventSummary",5);//1PHM
  if(Is1PHHfired) FillHistogramTH1(fOutputContainer,"hEventSummary",6);//1PHH

  if(fIsPHOSTriggerAnalysis && !(fPHOSTriggerHelper->IsPHI7(fEvent,fPHOSClusterCuts))){
    AliInfo("event is rejected. IsPHI7 = kFALSE.");
    return;
  }

  fPHOSClusterArray = (TClonesArray*)fEvent->FindListObject("PHOSClusterArray");
  if(!fPHOSClusterArray){
    AliWarning("fPHOSClusterArray object not found!");
    return;
  }

  //<- end of event selection
  //-> start physics analysis

  FillHistogramTH2(fOutputContainer,"hCentralityV0MvsCL0",fCentralityV0M,fCentralityCL0);
  FillHistogramTH2(fOutputContainer,"hCentralityV0MvsCL1",fCentralityV0M,fCentralityCL1);
  FillHistogramTH2(fOutputContainer,"hCentralityCL0vsCL1",fCentralityCL0,fCentralityCL1);
  FillHistogramTH2(fOutputContainer,"hCentralityV0AvsV0C",fCentralityV0A,fCentralityV0C);
  FillHistogramTH2(fOutputContainer,"hCentralityZNAvsZNC",fCentralityZNA,fCentralityZNC);

  FillHistogramTH1(fOutputContainer,"hVertexZSelectEvent" ,fVertex[2]);
  FillHistogramTH2(fOutputContainer,"hVertexXYSelectEvent",fVertex[0],fVertex[1]);

  FillHistogramTH1(fOutputContainer,"hEventSummary",2);//selected event

  if(fIsFlowTask){
    //fFlowQnVectorMgr->GetQnVectorList()->Print("",-1);
    const Int_t myHarmonic = 2;
    TList* qnlist = fFlowQnVectorMgr->GetQnVectorList();

    const AliQnCorrectionsQnVector *QnVectorTPCDet[3];
    Double_t TPCEP[3] = {};
    for(Int_t i=0;i<3;i++){
      QnVectorTPCDet[i] = GetQnVectorFromList(qnlist,Form("%s%s",fTPCEPName[i].Data(),fQnEstimator.Data()));
      if(QnVectorTPCDet[i]) TPCEP[i] = QnVectorTPCDet[i]->EventPlane(myHarmonic);
      if(TPCEP[i] < 0) TPCEP[i] += 2./(Double_t) myHarmonic * TMath::Pi();
      FillHistogramTH2(fOutputContainer,Form("hCentrality%svsEventPlane%s%s",fEstimator.Data(),fTPCEPName[i].Data(),fQnEstimator.Data()),fCentralityMain,TPCEP[i]);
      AliInfo(Form("TPC sub detector name %s : event plane = %f (rad).",fTPCEPName[i].Data(),TPCEP[i]));
    }

    const AliQnCorrectionsQnVector *QnVectorV0Det[3];
    Double_t V0EP[3]  = {};
    for(Int_t i=0;i<3;i++){
      QnVectorV0Det[i]  = GetQnVectorFromList(qnlist,Form("%s%s",fV0EPName[i].Data(),fQnEstimator.Data()));
      if(QnVectorV0Det[i]) V0EP[i] = QnVectorV0Det[i]->EventPlane(myHarmonic);
      if(V0EP[i] < 0)  V0EP[i]  += 2./(Double_t) myHarmonic * TMath::Pi();
      FillHistogramTH2(fOutputContainer,Form("hCentrality%svsEventPlane%s%s",fEstimator.Data(),fV0EPName[i].Data(),fQnEstimator.Data()),fCentralityMain,V0EP[i]);
      AliInfo(Form("V0 sub detector name %s : event plane = %f (rad)." ,fV0EPName[i].Data() ,V0EP[i]));
    }

    //0 < event plane < pi.
    //fEventPlane = TPCEP[0];//full acceptance of TPC
    fEventPlane = V0EP[0];//full V0

    //for event plane resolution
    //V0(main)-TPCA-TPCC
    FillHistogramTH2(fOutputContainer,Form("hCentrality%svsDeltaEventPlane%s%s",fEstimator.Data(),fV0EPName[0].Data() ,fTPCEPName[1].Data()),fCentralityMain,TMath::Cos(myHarmonic * (V0EP[0]  - TPCEP[1])));
    FillHistogramTH2(fOutputContainer,Form("hCentrality%svsDeltaEventPlane%s%s",fEstimator.Data(),fV0EPName[0].Data() ,fTPCEPName[2].Data()),fCentralityMain,TMath::Cos(myHarmonic * (V0EP[0]  - TPCEP[2])));
    FillHistogramTH2(fOutputContainer,Form("hCentrality%svsDeltaEventPlane%s%s",fEstimator.Data(),fTPCEPName[1].Data(),fTPCEPName[2].Data()),fCentralityMain,TMath::Cos(myHarmonic * (TPCEP[1] - TPCEP[2])));

    //TPC(main)-V0A-V0C
    FillHistogramTH2(fOutputContainer,Form("hCentrality%svsDeltaEventPlane%s%s",fEstimator.Data(),fTPCEPName[0].Data(),fV0EPName[1].Data()),fCentralityMain,TMath::Cos(myHarmonic * (TPCEP[0] - V0EP[1])));
    FillHistogramTH2(fOutputContainer,Form("hCentrality%svsDeltaEventPlane%s%s",fEstimator.Data(),fTPCEPName[0].Data(),fV0EPName[2].Data()),fCentralityMain,TMath::Cos(myHarmonic * (TPCEP[0] - V0EP[2])));
    FillHistogramTH2(fOutputContainer,Form("hCentrality%svsDeltaEventPlane%s%s",fEstimator.Data(),fV0EPName[1].Data() ,fV0EPName[2].Data()),fCentralityMain,TMath::Cos(myHarmonic * (V0EP[1]  - V0EP[2])));

  }

  
  if(!fIsMC) FillRejectionFactorMB();

  if(fRunNumber != fEvent->GetRunNumber()){ // Check run number
    fRunNumber = fEvent->GetRunNumber();
    fPHOSGeo = GetPHOSGeometry();
  }

  fPIDResponse = fInputHandler->GetPIDResponse();
  if(!fPIDResponse){
    AliWarning("fPIDResponse does not exist! This is not crucial in photon analysis.");
  }

  //track QA
  TrackQA();
  if(fCollisionSystem==0){
    //not centrality, but track multiplicity
    //track QA to extrack number of hybrid track for centrality binning
    AliInfo(Form("This is pp analysis. Nybrid = %d , Zvtx = %f cm, fZvtx = %d.",fNHybridTrack,fVertex[2],fZvtx));
  }

  if(fIsMC) ProcessMC();

  if(fIsMC){//fill cluster occupancy in M.C. to obtain total(i.e. HIJNG/HIJING+JJ) Ncluster.
    Int_t multPHOSClustAll = 0;
    for(Int_t i1=0;i1<fPHOSClusterArray->GetEntriesFast();i1++){
      AliCaloPhoton *ph = (AliCaloPhoton*)fPHOSClusterArray->At(i1);
      if(!fPHOSClusterCuts->AcceptPhoton(ph)) continue;
      multPHOSClustAll++;
    }
    FillHistogramTH2(fOutputContainer,Form("hCentrality%svsPHOSClusterMultiplicityMC",fEstimator.Data()),fCentralityMain  ,multPHOSClustAll);
  }

  Int_t multPHOSClustnoPID = 0;
  Int_t multPHOSClustnoPIDTOF = 0;
  for(Int_t i1=0;i1<fPHOSClusterArray->GetEntriesFast();i1++){
    AliCaloPhoton *ph = (AliCaloPhoton*)fPHOSClusterArray->At(i1);
    multPHOSClustnoPID++;
    if(ph->IsTOFOK()) multPHOSClustnoPIDTOF++;
  }
  FillHistogramTH2(fOutputContainer,Form("hCentrality%svsPHOSClusterMultiplicityNoPID",fEstimator.Data()),fCentralityMain   ,multPHOSClustnoPID);
  FillHistogramTH2(fOutputContainer,Form("hCentrality%svsPHOSClusterMultiplicityNoPIDTOF",fEstimator.Data()),fCentralityMain,multPHOSClustnoPIDTOF);

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

        if(fIsJJMC){
          if(fMCType.Contains("JJMC") && (primary < firstJetindex || lastJetindex < primary)) fPHOSClusterArray->Remove(ph);
          if(fMCType.Contains("MBMC") && (primary < firstUEindex  || lastUEindex  < primary)) fPHOSClusterArray->Remove(ph);
        }

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

        if(fMCType.Contains("JJMC") && p->GetGeneratorIndex() != genIDJet) fPHOSClusterArray->Remove(ph);
        if(fMCType.Contains("MBMC") && p->GetGeneratorIndex() != genIDUE)  fPHOSClusterArray->Remove(ph);
      }
      fPHOSClusterArray->Compress();
    }
  }
  if(!fPHOSEvents[fZvtx]) fPHOSEvents[fZvtx] = new TList();
  TList *prevPHOS = fPHOSEvents[fZvtx];

  if(fIsPHOSTriggerAnalysis){
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

  EstimatePIDCutEfficiency();
  EstimateTOFCutEfficiency();

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

  if(fJJMCHandler) delete fJJMCHandler;

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

  Double_t pT=0, eta=0, phi=0, dEdx=0, p=0;
  Int_t NHybrid=0;
  Int_t NGlobal=0;
  Int_t NComplementary=0;

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

      if(aodtrack->IsHybridGlobalConstrainedGlobal()){//hybrid track
        NHybrid++;
        pT = aodtrack->Pt();
        eta = aodtrack->Eta();
        phi = aodtrack->Phi();
        if(phi<0) phi += TMath::TwoPi();
        p = aodtrack->P();
        dEdx = aodtrack->GetTPCsignal();

        FillHistogramTH2(fOutputContainer,"hTrackTPCdEdx",p,dEdx);

        FillHistogramTH1(fOutputContainer,"hHybridTrackPt",pT);
        FillHistogramTH2(fOutputContainer,"hHybridTrackEtaPhi",phi,eta);

        if(aodtrack->IsGlobalConstrained()){//constrained to primary vertex, instead of SPD hits.//complementary track
          NComplementary++;
          FillHistogramTH1(fOutputContainer,"hComplementaryTrackPt",pT);
          FillHistogramTH2(fOutputContainer,"hComplementaryTrackEtaPhi",phi,eta);
        }
        else{//global track
          NGlobal++;
          FillHistogramTH1(fOutputContainer,"hGlobalTrackPt",pT);
          FillHistogramTH2(fOutputContainer,"hGlobalTrackEtaPhi",phi,eta);
        }

      }
    }//end of track loop


  }//end of AOD

  FillHistogramTH1(fOutputContainer,"hHybridTrackMult",NHybrid);
  FillHistogramTH1(fOutputContainer,"hGlobalTrackMult",NGlobal);
  FillHistogramTH1(fOutputContainer,"hComplementaryTrackMult",NComplementary);

  FillHistogramTH2(fOutputContainer,Form("hCentrality%svsTrackMultiplicity",fEstimator.Data()),fCentralityMain,NHybrid);
  fNHybridTrack = NHybrid;

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
  AliVCaloCells *cells = dynamic_cast<AliVCaloCells*>(fEvent->GetPHOSCells());
  AliVCaloTrigger* trg = fEvent->GetCaloTrigger("PHOS");
  trg->Reset();

  const Int_t Nfired = trg->GetEntries();
  //if(fPHOSTriggerAnalysis) cout << "Number of fired region = " <<trg->GetEntries() << endl;
  FillHistogramTH1(fOutputContainer,"hNFiredTRUChannel",Nfired);

  Int_t trgrelId[4]={};
  Int_t tmod=0; //"Online" module number
  Int_t trgabsId=0; //bottom-left 4x4 edge cell absId [1-64], [1-56]
  Int_t trgmodule=0,trgcellx=0,trgcellz=0;

  Int_t relId[4]={};
  Int_t module=0,cellx=0,cellz=0;
  Double_t energy=0;

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
    if(L0input == 9) L1 = -1;//L0
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

    //cout << "trgmodule = "<< trgmodule << endl;

    if(trgmodule < 1 || trgmodule > 4){
      AliError(Form("Wrong module number %d",trgmodule));
      return;
    }

    FillHistogramTH1(fOutputContainer,Form("hFiredTRUChannelM%d",trgmodule),trgcellx,trgcellz);

    for(Int_t i=0;i<multClust;i++){
      AliCaloPhoton *ph = (AliCaloPhoton*)fPHOSClusterArray->At(i);
      AliVCluster *clu1 = (AliVCluster*)ph->GetCluster();
      energy = ph->Energy();

      Int_t maxAbsId = fPHOSTriggerHelper->FindHighestAmplitudeCellAbsId(clu1, cells);
      //cout << "maxAbsId = " << maxAbsId << endl;

      fPHOSGeo->AbsToRelNumbering(maxAbsId,relId);
      module = relId[0];
      cellx  = relId[2];
      cellz  = relId[3];

      if(module < 1 || module > 4){
        AliError(Form("Wrong module number %d",module));
        return;
      }

      //cout << "module = " << module << " , cellx = " << cellx << " , cellz = " << cellz << endl;
      Int_t diffx = trgcellx - cellx;
      Int_t diffz = trgcellz - cellz;
      FillHistogramTH3(fOutputContainer,Form("hDistanceXZM%dTRU%d",trgmodule,fPHOSTriggerHelper->WhichTRU(trgcellx,trgcellz)),diffx,diffz,energy);

      if(!IsFilledOnce){//to avoid double counting of energy
        FillHistogramTH1(fOutputContainer,Form("hClusterEnergyM%d",module),energy);
        FillHistogramTH1(fOutputContainer,Form("hClusterEnergyM%dTRU%d",module,fPHOSTriggerHelper->WhichTRU(cellx,cellz)),energy);
        FillHistogramTH3(fOutputContainer,Form("hClusterHitM%d",module),cellx,cellz,energy);
      }

      if(fPHOSTriggerHelper->IsMatched(trgrelId,relId)){
        FillHistogramTH3(fOutputContainer,Form("hMatchedFiredTRUChannelM%d",trgmodule),trgcellx,trgcellz,energy);
        FillHistogramTH3(fOutputContainer,Form("hMatchedDistanceXZM%dTRU%d",trgmodule,fPHOSTriggerHelper->WhichTRU(trgcellx,trgcellz)),diffx,diffz,energy);

        if(kUsedCluster[i] == 0){
          FillHistogramTH1(fOutputContainer,Form("hMatchedClusterEnergyM%d",module),energy);
          FillHistogramTH1(fOutputContainer,Form("hMatchedClusterEnergyM%dTRU%d",trgmodule,fPHOSTriggerHelper->WhichTRU(trgcellx,trgcellz)),energy);
          FillHistogramTH3(fOutputContainer,Form("hMatchedClusterHitM%d",module),cellx,cellz,energy);
        }
        kUsedCluster[i] = 1;
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
  Double_t DistToBadChannel = 0;
  Double_t M02=0;
  Double_t R = 0, coreR=0;
  Double_t coreE = 0;
  Double_t r=999, trackPt=0;

  for(Int_t i1=0;i1<multClust;i1++){
    AliCaloPhoton *ph = (AliCaloPhoton*)fPHOSClusterArray->At(i1);
    if(!fPHOSClusterCuts->AcceptPhoton(ph)) continue;

    energy  = ph->Energy();
    if(fUseCoreEnergy){
      energy = (ph->GetMomV2())->Energy();
   }

    digMult = ph->GetNCells();
    tof     = ph->GetTime();//unit is second.
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

    DistToBadChannel = ph->DistToBadfp();//in unit of cell with floating point
    FillHistogramTH2(fOutputContainer,"hEnergyvsDistanceToBadChannel",energy,DistToBadChannel);

    eta = global1.Eta();
    phi = global1.Phi();
    if(phi<0) phi += TMath::TwoPi();
    FillHistogramTH2(fOutputContainer,"hClusterEtaPhi",phi,eta);

    FillHistogramTH2(fOutputContainer,Form("hClusterEvsNM%d",module),energy,digMult);
    FillHistogramTH2(fOutputContainer,Form("hClusterEvsTM%d",module),energy,tof*1e+9);
    FillHistogramTH2(fOutputContainer,Form("hCluNXZM%d",module),cellx,cellz);
    FillHistogramTH2(fOutputContainer,Form("hCluEXZM%d",module),cellx,cellz,energy);
   
    R     = ph->GetNsigmaFullDisp();
    coreR = ph->GetNsigmaCoreDisp();
    coreE = (ph->GetMomV2())->Energy();
    AliVCluster *clu = (AliVCluster*)ph->GetCluster();
    AliVTrack *track = 0x0;

    if(fESDEvent){
      Int_t trackindex = clu->GetTrackMatchedIndex();
      if(trackindex > 0){
        track = (AliVTrack*)(fEvent->GetTrack(trackindex));
      }//end of track matching
    }//end of ESD
    else if(fAODEvent){
      if(clu->GetNTracksMatched() > 0){
        track = dynamic_cast<AliVTrack*>(clu->GetTrackMatched(0));
      }//end of track matching
    }//end of AOD

    if(track){
      trackPt = track->Pt();
      r = ph->GetNsigmaCPV();
      FillHistogramTH2(fOutputContainer,Form("hRvsTrackPtM%d",module),trackPt,r);
    } 
 
    FillHistogramTH2(fOutputContainer,Form("hFullDispvsFullEM%d",module),energy,R);
    FillHistogramTH2(fOutputContainer,Form("hCoreDispvsCoreEM%d",module),coreE,coreR);
    FillHistogramTH2(fOutputContainer,Form("hFullDispvsCoreEM%d",module),coreE,R);
    FillHistogramTH2(fOutputContainer,Form("hCoreDispvsFullEM%d",module),energy,coreR);


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

    FillHistogramTH2(fOutputContainer,Form("hClusterEvsM02M%d",module),energy,M02);
    FillHistogramTH2(fOutputContainer,Form("hClusterNvsM02M%d",module),digMult,M02);
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

}
//________________________________________________________________________
void AliAnalysisTaskPHOSPi0EtaToGammaGamma::FillPhoton() 
{
  const Int_t multClust = fPHOSClusterArray->GetEntriesFast();

  Double_t pT=0,energy=0;
  Double_t phi = -999, dphi = -999.;
  Double_t eff=1;
  TF1 *f1tof = GetTOFCutEfficiencyFunction();
  Double_t value[4] = {};

  Double_t weight = 1.;
  Int_t primary = -1;

  for(Int_t iph=0;iph<multClust;iph++){
    AliCaloPhoton *ph = (AliCaloPhoton*)fPHOSClusterArray->At(iph);
    weight = 1.;

    if(!fIsMC && fIsPHOSTriggerAnalysis && !ph->IsTrig()) continue;//it is meaningless to focus on photon without fired trigger in PHOS triggered data.

    pT = ph->Pt();
    energy = ph->Energy();
    phi = ph->Phi();

    if(fUseCoreEnergy){
      pT = (ph->GetMomV2())->Pt();
      energy = (ph->GetMomV2())->Energy();
      phi = (ph->GetMomV2())->Phi();
    }

    eff = f1tof->Eval(energy);

    value[0] = ph->GetNsigmaCPV();
    value[1] = ph->GetNsigmaCoreDisp();
    value[2] = ph->DistToBadfp();
    value[3] = pT;

    //for PID systematic study
    if(fIsMC) FillSparse(fOutputContainer,"hSparsePhoton",value,weight);
    else{
      if(ph->IsTOFOK()){
        FillSparse(fOutputContainer,"hSparsePhoton",value,1/eff * weight);
      }
    }

    if(!fPHOSClusterCuts->AcceptPhoton(ph)) continue;

    if(fIsMC){
      primary = ph->GetPrimary();
      weight = ph->GetWeight();
      if(IsPhoton(primary)) FillHistogramTH1(fOutputContainer,"hPurityGamma",pT,weight);
    }

    //0 < photon phi < 2pi
    if(phi < 0) phi += TMath::TwoPi();
    if(fIsFlowTask){
      dphi = DeltaPhiIn0Pi(phi - fEventPlane);
    }
    else dphi = DeltaPhiIn0Pi(phi);

    FillHistogramTH2(fOutputContainer,"hPhotonPt",pT,dphi*TMath::RadToDeg(),weight);

    if(ph->IsTOFOK()){
      FillHistogramTH2(fOutputContainer,"hPhotonPt_TOF",pT,dphi*TMath::RadToDeg(),1/eff * weight);
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

  Double_t eff1=1, eff2=1, eff12=1;
  TF1 *f1tof = GetTOFCutEfficiencyFunction();

  Double_t weight = 1., w1 = 1., w2 = 1.;

  Int_t primary1 = -1;
  Int_t primary2 = -1;
  Double_t TrueK0SPt = 0;
  Double_t TrueL0Pt = 0;

  Int_t commonID = -1;

  for(Int_t i1=0;i1<multClust-1;i1++){
    AliCaloPhoton *ph1 = (AliCaloPhoton*)fPHOSClusterArray->At(i1);
    if(!fPHOSClusterCuts->AcceptPhoton(ph1)) continue;

    for(Int_t i2=i1+1;i2<multClust;i2++){
      AliCaloPhoton *ph2 = (AliCaloPhoton*)fPHOSClusterArray->At(i2);
      if(!fPHOSClusterCuts->AcceptPhoton(ph2)) continue;

      if(!fIsMC && fIsPHOSTriggerAnalysis && (!ph1->IsTrig() && !ph2->IsTrig())) continue;//it is meaningless to reconstruct invariant mass with FALSE-FALSE combination in PHOS triggered data.

      e1 = ph1->Energy();
      e2 = ph2->Energy();

      p12  = *ph1 + *ph2;
      m12  = p12.M();
      pt12 = p12.Pt();
      phi  = p12.Phi();
      asym = TMath::Abs((ph1->Energy()-ph2->Energy())/(ph1->Energy()+ph2->Energy()));//always full energy

      if(fUseCoreEnergy){
        p12core = *(ph1->GetMomV2()) + *(ph2->GetMomV2());
        m12     = p12core.M();
        pt12    = p12core.Pt();
        phi     = p12core.Phi();

        e1 = (ph1->GetMomV2())->Energy();
        e2 = (ph2->GetMomV2())->Energy();
      }

      eff1 = f1tof->Eval(e1);
      eff2 = f1tof->Eval(e2);
      eff12 = eff1 * eff2;

      weight = 1.;
      if(fIsMC){
        w1= ph1->GetWeight();
        primary1 = ph1->GetPrimary();

        w2 = ph2->GetWeight();
        primary2 = ph2->GetPrimary();

        commonID = FindCommonParent(primary1,primary2);
        if(commonID > -1) weight = w1;
        else weight = w1*w2;

        if(IsFrom(primary1,TrueK0SPt,310) && IsFrom(primary2,TrueK0SPt,310) && commonID > -1){
          //for feed down correction from K0S->pi0 + pi0
          FillHistogramTH2(fOutputContainer,"hMggFromK0S",m12,pt12,weight);
          if(asym < 0.8) FillHistogramTH2(fOutputContainer,"hMggFromK0S_asym08",m12,pt12,weight);
        }
        if(IsFrom(primary1,TrueL0Pt,3122) && IsFrom(primary2,TrueL0Pt,3122) && commonID > -1){//this contibution is neglegible
          //for feed down correction from L0->pi0 + neutron
          FillHistogramTH2(fOutputContainer,"hMggFromLambda0",m12,pt12,weight);
          if(asym < 0.8) FillHistogramTH2(fOutputContainer,"hMggFromLambda0_asym08",m12,pt12,weight);
        }

      }//end of if fIsMC

      if(phi < 0) phi += TMath::TwoPi();
      if(fIsFlowTask){
        dphi = DeltaPhiIn0Pi(phi - fEventPlane);
      }
      else dphi = DeltaPhiIn0Pi(phi);

      if(TMath::Abs(ph1->Module()-ph2->Module()) < 2) FillHistogramTH2(fOutputContainer,Form("hMgg_M%d%d",TMath::Min(ph1->Module(),ph2->Module()), TMath::Max(ph1->Module(),ph2->Module())),m12,pt12,weight);
      FillHistogramTH3(fOutputContainer,"hMgg",m12,pt12,dphi*TMath::RadToDeg(),weight);

      FillHistogramTH2(fOutputContainer,"hAsymvsMgg",asym,m12,weight);
      if(0.12 < m12 && m12 < 0.15) FillHistogramTH2(fOutputContainer,"hAsymvsPt",asym,pt12,weight);

      if(asym < 0.8) FillHistogramTH3(fOutputContainer,"hMgg_asym08",m12,pt12,dphi*TMath::RadToDeg(),weight);

      if(ph1->IsTOFOK() && ph2->IsTOFOK()){
        FillHistogramTH3(fOutputContainer,"hMgg_TOF",m12,pt12,dphi*TMath::RadToDeg(),1/eff12 * weight);

        if(TMath::Abs(ph1->Module()-ph2->Module()) < 2) FillHistogramTH2(fOutputContainer,Form("hMgg_M%d%d_TOF",TMath::Min(ph1->Module(),ph2->Module()), TMath::Max(ph1->Module(),ph2->Module())),m12,pt12,1/eff12 * weight);

        if(asym < 0.8) FillHistogramTH3(fOutputContainer,"hMgg_TOF_asym08",m12,pt12,dphi*TMath::RadToDeg(),1/eff12 * weight);

      }//end of TOF cut

      if(asym < 0.1){
        FillHistogramTH2(fOutputContainer,"hMggvsE_asym01",m12,e1,weight);
        FillHistogramTH2(fOutputContainer,"hMggvsE_asym01",m12,e2,weight);
      }//end of asym01

    }//end of ph2

  }//end of ph1

}
//________________________________________________________________________
void AliAnalysisTaskPHOSPi0EtaToGammaGamma::FillMixMgg() 
{
  TList *prevPHOS = fPHOSEvents[fZvtx];

  const Int_t multClust = fPHOSClusterArray->GetEntriesFast();

  TLorentzVector p12, p12core;
  Double_t m12=0,pt12=0,asym=0;
  Double_t phi = -999, dphi = -999.;
  Double_t weight = 1., w1 = 1., w2 = 1.;

  Double_t eff1=1, eff2=1, eff12=1;
  Double_t e1=0,e2=0;
  TF1 *f1tof = GetTOFCutEfficiencyFunction();

  for(Int_t i1=0;i1<multClust;i1++){
    AliCaloPhoton *ph1 = (AliCaloPhoton*)fPHOSClusterArray->At(i1);
    if(!fPHOSClusterCuts->AcceptPhoton(ph1)) continue;

    for(Int_t ev=0;ev<prevPHOS->GetSize();ev++){
      TClonesArray *mixPHOS = static_cast<TClonesArray*>(prevPHOS->At(ev));

      for(Int_t i2=0;i2<mixPHOS->GetEntriesFast();i2++){
        AliCaloPhoton *ph2 = (AliCaloPhoton*)mixPHOS->At(i2);
        if(!fPHOSClusterCuts->AcceptPhoton(ph2)) continue;

        if(!fIsMC && fIsPHOSTriggerAnalysis && (!ph1->IsTrig() && !ph2->IsTrig())) continue;//it is meaningless to reconstruct invariant mass with FALSE-FALSE combination in PHOS triggered data.

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
        }

        eff1 = f1tof->Eval(e1);
        eff2 = f1tof->Eval(e2);
        eff12 = eff1 * eff2;
        weight = 1.;

        if(fIsMC){
          w1= ph1->GetWeight();
          w2 = ph2->GetWeight();

          weight = w1*w2;

        }//end of if fIsMC

        if(phi < 0) phi += TMath::TwoPi();

        if(fIsFlowTask){
          dphi = DeltaPhiIn0Pi(phi - fEventPlane);
        }
        else dphi = DeltaPhiIn0Pi(phi);

        FillHistogramTH3(fOutputContainer,"hMixMgg",m12,pt12,dphi*TMath::RadToDeg(),weight);
        if(TMath::Abs(ph1->Module()-ph2->Module()) < 2) FillHistogramTH2(fOutputContainer,Form("hMixMgg_M%d%d",TMath::Min(ph1->Module(),ph2->Module()), TMath::Max(ph1->Module(),ph2->Module())),m12,pt12);

        if(asym < 0.8) FillHistogramTH3(fOutputContainer,"hMixMgg_asym08",m12,pt12,dphi*TMath::RadToDeg(),weight);

        if(ph1->IsTOFOK() && ph2->IsTOFOK()){
          FillHistogramTH3(fOutputContainer,"hMixMgg_TOF",m12,pt12,dphi*TMath::RadToDeg(),1/eff12 * weight);
          if(TMath::Abs(ph1->Module()-ph2->Module()) < 2) FillHistogramTH2(fOutputContainer,Form("hMixMgg_M%d%d_TOF",TMath::Min(ph1->Module(),ph2->Module()), TMath::Max(ph1->Module(),ph2->Module())),m12,pt12,1/eff12);

          if(asym < 0.8) FillHistogramTH3(fOutputContainer,"hMixMgg_TOF_asym08",m12,pt12,dphi*TMath::RadToDeg(),1/eff12 * weight);

        }//end of TOF cut

        if(asym < 0.1){
          FillHistogramTH2(fOutputContainer,"hMixMggvsE_asym01",m12,e1,weight);
          FillHistogramTH2(fOutputContainer,"hMixMggvsE_asym01",m12,e2,weight);
        }//end of asym01

      }//end of ph2

    }//end of mix

  }//end of ph1

}
//________________________________________________________________________
void AliAnalysisTaskPHOSPi0EtaToGammaGamma::FillEpRatio() 
{
  //since this E/p is for global energy scale and mis-alignment study, neither TOF cut nor additional weight in M.C. are applied.

  const Int_t multClust = fPHOSClusterArray->GetEntriesFast();

  const Double_t NsigmaCPV = 2.0;
  const Double_t NsigmaDisp = fPHOSClusterCuts->GetDispParameter();

  AliPHOSClusterCuts *cuts = new AliPHOSClusterCuts("CutsForCharged");//for charged tracks
  cuts->SetUseCoreDispersion(fUseCoreEnergy);
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
  Bool_t isHybridTrack = kFALSE;

  //Double_t weight = 1.;

  for(Int_t iph=0;iph<multClust;iph++){
    AliCaloPhoton *ph = (AliCaloPhoton*)fPHOSClusterArray->At(iph);
    if(!cuts->AcceptChargedParticle(ph)) continue;

    //weight = 1.;
    //if(fIsMC) weight = ph->GetWeight();

    energy = ph->Energy();
    if(fUseCoreEnergy) energy = (ph->GetMomV2())->Energy();

    //eff = f1tof->Eval(energy);

    cluster = (AliVCluster*)ph->GetCluster();
    PID_ele = kFALSE;
    isHybridTrack = kFALSE;
    track = 0x0;

    if(fESDEvent){
      Int_t trackindex = cluster->GetTrackMatchedIndex();
      if(trackindex > 0){
        track = (AliVTrack*)(fEvent->GetTrack(trackindex));
        isHybridTrack = fESDtrackCutsGlobal->AcceptTrack(dynamic_cast<AliESDtrack*>(track)) || fESDtrackCutsGlobalConstrained->AcceptTrack(dynamic_cast<AliESDtrack*>(track));
      }
    }//end of ESD
    else if(fAODEvent){
      if(cluster->GetNTracksMatched() > 0){
        track = dynamic_cast<AliVTrack*>(cluster->GetTrackMatched(0));
        isHybridTrack = dynamic_cast<AliAODTrack*>(track)->IsHybridGlobalConstrainedGlobal();//hybrid track

      }//end of track matching
    }//end of AOD

    if(track && isHybridTrack){
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
      else if(nsigmaElectron < -3 || 5 < nsigmaElectron){//non-electron and far from border
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
  Double_t energy=0;

  for(Int_t i1=0;i1<multClust;i1++){
    AliCaloPhoton *ph1 = (AliCaloPhoton*)fPHOSClusterArray->At(i1);
    if(!fPHOSClusterCuts->AcceptPhoton(ph1)) continue;

    for(Int_t i2=0;i2<multClust;i2++){
      AliCaloPhoton *ph2 = (AliCaloPhoton*)fPHOSClusterArray->At(i2);

      if(i2==i1) continue;//reject same cluster combination

      p12 = *ph1 + *ph2;
      m12 = p12.M();
      energy = ph2->Energy();

      if(fUseCoreEnergy){
        p12core = *(ph1->GetMomV2()) + *(ph2->GetMomV2());
        m12 = p12core.M();
        energy = (ph2->GetMomV2())->Energy();
      }

      FillHistogramTH2(fOutputContainer,"hMgg_Probe_PID",m12,energy);
      if(fPHOSClusterCuts->AcceptPhoton(ph2))
        FillHistogramTH2(fOutputContainer,"hMgg_PassingProbe_PID",m12,energy);

    }//end of ph2

  }//end of ph1

  //next mixed event
  TList *prevPHOS = fPHOSEvents[fZvtx];

  for(Int_t i1=0;i1<multClust;i1++){
    AliCaloPhoton *ph1 = (AliCaloPhoton*)fPHOSClusterArray->At(i1);
    if(!fPHOSClusterCuts->AcceptPhoton(ph1)) continue;

    for(Int_t ev=0;ev<prevPHOS->GetSize();ev++){
      TClonesArray *mixPHOS = static_cast<TClonesArray*>(prevPHOS->At(ev));

      for(Int_t i2=0;i2<mixPHOS->GetEntriesFast();i2++){
        AliCaloPhoton *ph2 = (AliCaloPhoton*)mixPHOS->At(i2);

        p12 = *ph1 + *ph2;
        m12 = p12.M();
        energy = ph2->Energy();

        if(fUseCoreEnergy){
          p12core = *(ph1->GetMomV2()) + *(ph2->GetMomV2());
          m12 = p12core.M();
          energy = (ph2->GetMomV2())->Energy();
        }

        FillHistogramTH2(fOutputContainer,"hMixMgg_Probe_PID",m12,energy);
        if(fPHOSClusterCuts->AcceptPhoton(ph2))
          FillHistogramTH2(fOutputContainer,"hMixMgg_PassingProbe_PID",m12,energy);

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

    for(Int_t i2=0;i2<multClust;i2++){
      AliCaloPhoton *ph2 = (AliCaloPhoton*)fPHOSClusterArray->At(i2);

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
  TList *prevPHOS = fPHOSEvents[fZvtx];

  for(Int_t i1=0;i1<multClust;i1++){
    AliCaloPhoton *ph1 = (AliCaloPhoton*)fPHOSClusterArray->At(i1);
    if(!fPHOSClusterCuts->AcceptPhoton(ph1)) continue;
    if(!ph1->IsTOFOK()) continue;

    for(Int_t ev=0;ev<prevPHOS->GetSize();ev++){
      TClonesArray *mixPHOS = static_cast<TClonesArray*>(prevPHOS->At(ev));

      for(Int_t i2=0;i2<mixPHOS->GetEntriesFast();i2++){
        AliCaloPhoton *ph2 = (AliCaloPhoton*)mixPHOS->At(i2);

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
void AliAnalysisTaskPHOSPi0EtaToGammaGamma::SelectTriggeredCluster()
{

  Int_t trgrelId[4]={};
  Int_t tmod=0; //"Online" module number
  Int_t trgabsId=0; //bottom-left 4x4 edge cell absId [1-64], [1-56]
  Int_t trgmodule=0,trgcellx=0,trgcellz=0;

  Int_t relId[4]={};
  Int_t module=0,cellx=0,cellz=0,tru=0;
  Double_t energy=0;
  
  AliVCaloCells *cells = dynamic_cast<AliVCaloCells*>(fEvent->GetPHOSCells());
  const Int_t L1input = fPHOSTriggerHelper->GetL1TriggerInput();
  const Int_t L0input = fPHOSTriggerHelper->GetL0TriggerInput();

  Int_t L1=-999;
  if(L1input > 0){
    if     (L1input == 7) L1 = 0;//L1 high
    else if(L1input == 6) L1 = 1;//L1 medium
    else if(L1input == 5) L1 = 2;//L1 low
  }
  else if(L0input > 0){
    if(L0input == 9) L1 = -1;//L0
  }

  const Int_t multClust = fPHOSClusterArray->GetEntriesFast();
  //set label to all clusters. PID/TOF cut is evaluated EstimateTOFCutEfficiency, EstimatePIDCutEfficiency.
  for(Int_t i=0;i<multClust;i++){
    AliCaloPhoton *ph = (AliCaloPhoton*)fPHOSClusterArray->At(i);
    AliVCluster *clu1 = (AliVCluster*)ph->GetCluster();

    Int_t maxAbsId = fPHOSTriggerHelper->FindHighestAmplitudeCellAbsId(clu1, cells);

    fPHOSGeo->AbsToRelNumbering(maxAbsId,relId);
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

      if(fPHOSTriggerHelper->IsMatched(trgrelId,relId)){
        ph->SetTrig(kTRUE);
        break;
      }
    }//end of while trigger patch loop

  }//end of cluster loop

  relId[0] = 0; relId[1] = 0; relId[2] = 0; relId[3] = 0;
  module = 0; cellx = 0; cellz = 0; tru = 0;

  for(Int_t i=0;i<multClust;i++){
    AliCaloPhoton *ph = (AliCaloPhoton*)fPHOSClusterArray->At(i);
    if(!fPHOSClusterCuts->AcceptPhoton(ph)) continue;

    energy = ph->Energy();
    if(fUseCoreEnergy) energy = (ph->GetMomV2())->Energy();

    AliVCluster *clu1 = (AliVCluster*)ph->GetCluster();
    Int_t maxAbsId = fPHOSTriggerHelper->FindHighestAmplitudeCellAbsId(clu1, cells);
    fPHOSGeo->AbsToRelNumbering(maxAbsId,relId);
    module = relId[0];
    cellx  = relId[2];
    cellz  = relId[3];
    tru = fPHOSTriggerHelper->WhichTRU(cellx,cellz);

    FillHistogramTH1(fOutputContainer,Form("hClusterE_M%d_TRU%d",module,tru),energy);
    if(ph->IsTrig()) FillHistogramTH1(fOutputContainer,Form("hTriggeredClusterE_M%d_TRU%d",module,tru),energy);
    //for QA purpose TOF cut efficiency is not corrected.
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

  AliVCaloCells *cells = dynamic_cast<AliVCaloCells*>(fEvent->GetPHOSCells());
  const Int_t multClust = fPHOSClusterArray->GetEntriesFast();

  Double_t position[3] = {};
  Double_t energy=0;
  Int_t module_tag=0,cellx_tag=0,cellz_tag=0,tru_tag=0;
  Int_t module=0,cellx=0,cellz=0,tru=0;
  Int_t relId[4]={};

  Int_t truch_tag = -1, chX_tag=-1, chZ_tag=-1;
  Int_t truch = -1, chX=-1, chZ=-1;

  TLorentzVector p12, p12core;
  Double_t m12=0;
  Int_t maxAbsId = -1;

  for(Int_t i1=0;i1<multClust;i1++){
    AliCaloPhoton *ph1 = (AliCaloPhoton*)fPHOSClusterArray->At(i1);
    if(!fPHOSClusterCuts->AcceptPhoton(ph1)) continue;
    if(!ph1->IsTrig()) continue;
    //if(!ph1->IsTOFOK()) continue;

    relId[0] = 0; relId[1] = 0; relId[2] = 0; relId[3] = 0;

    //position[0] = ph1->EMCx();
    //position[1] = ph1->EMCy();
    //position[2] = ph1->EMCz();
    //TVector3 global_tag(position);
    //fPHOSGeo->GlobalPos2RelId(global_tag,relId);

    maxAbsId = fPHOSTriggerHelper->FindHighestAmplitudeCellAbsId(ph1->GetCluster(), cells);
    fPHOSGeo->AbsToRelNumbering(maxAbsId,relId);

    module_tag = relId[0];
    cellx_tag  = relId[2];
    cellz_tag  = relId[3];
    tru_tag = fPHOSTriggerHelper->WhichTRU(cellx_tag,cellz_tag);
    truch_tag = fPHOSTriggerHelper->WhichTRUChannel(cellx_tag,cellz_tag,chX_tag,chZ_tag);
    AliInfo(Form("module_tag = %d , tru_tag = %d , cellx_tag = %d , cellz_tag = %d , chX_tag = %d , chZ_tag = %d , truch_tag = %d.",module_tag, tru_tag, cellx_tag, cellz_tag, chX_tag, chZ_tag, truch_tag));

    for(Int_t i2=0;i2<multClust;i2++){
      AliCaloPhoton *ph2 = (AliCaloPhoton*)fPHOSClusterArray->At(i2);
      //if(!ph2->IsTOFOK()) continue;

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

      //position[0] = ph2->EMCx();
      //position[1] = ph2->EMCy();
      //position[2] = ph2->EMCz();
      //TVector3 global1(position);
      //fPHOSGeo->GlobalPos2RelId(global1,relId);

      maxAbsId = fPHOSTriggerHelper->FindHighestAmplitudeCellAbsId(ph2->GetCluster(), cells);
      fPHOSGeo->AbsToRelNumbering(maxAbsId,relId);

      module = relId[0];
      cellx  = relId[2];
      cellz  = relId[3];
      tru = fPHOSTriggerHelper->WhichTRU(cellx,cellz);
      truch = fPHOSTriggerHelper->WhichTRUChannel(cellx,cellz,chX,chZ);

      Int_t dm = module_tag - module;
      Int_t dt = tru_tag - tru;
      Int_t dx = chX_tag - chX;
      Int_t dz = chZ_tag - chZ;

      //if( (dm == 0)
      // && (dt == 0)
      // && (TMath::Abs(dx) < 2)
      // && (TMath::Abs(dz) < 2)
      // ) continue;//reject cluster pair where they belong to same 4x4 region.

      AliInfo(Form("dm = %d , dt = %d , dx = %d , dz = %d, truch = %d.",dm,dt,dx,dz,truch));

      FillHistogramTH2(fOutputContainer,"hMgg_Probe_Trg",m12,energy);
      FillHistogramTH2(fOutputContainer,Form("hMgg_Probe_Trg_M%d_TRU%d",module,tru),m12,energy);

      if(ph2->IsTrig()){
        FillHistogramTH2(fOutputContainer,"hMgg_PassingProbe_Trg",m12,energy);
        FillHistogramTH2(fOutputContainer,Form("hMgg_PassingProbe_Trg_M%d_TRU%d",module,tru),m12,energy);
      }

    }//end of ph2

  }//end of ph1

  //next mixed event
  TList *prevPHOS = fPHOSEvents[fZvtx];

  for(Int_t i1=0;i1<multClust;i1++){
    AliCaloPhoton *ph1 = (AliCaloPhoton*)fPHOSClusterArray->At(i1);
    if(!fPHOSClusterCuts->AcceptPhoton(ph1)) continue;
    if(!ph1->IsTrig()) continue;
    if(!ph1->IsTOFOK()) continue;

    for(Int_t ev=0;ev<prevPHOS->GetSize();ev++){
      TClonesArray *mixPHOS = static_cast<TClonesArray*>(prevPHOS->At(ev));

      for(Int_t i2=0;i2<mixPHOS->GetEntriesFast();i2++){
        AliCaloPhoton *ph2 = (AliCaloPhoton*)mixPHOS->At(i2);
        if(!ph2->IsTOFOK()) continue;

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
        TVector3 global1(position);
        fPHOSGeo->GlobalPos2RelId(global1,relId);

        module = relId[0];
        cellx  = relId[2];
        cellz  = relId[3];
        tru = fPHOSTriggerHelper->WhichTRU(cellx,cellz);

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

  TF1 *f1Pi0Weight = (TF1*)GetAdditionalPi0PtWeightFunction(fCentralityMain);
  TF1 *f1K0SWeight = (TF1*)GetAdditionalK0SPtWeightFunction(fCentralityMain);

  Double_t TruePi0Pt = 1.;
  Double_t TrueK0SPt = 0;

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
      Int_t primary = FindPrimaryMotherESD(i);

      if(fIsJJMC){
        if(fMCType.Contains("JJMC") && (primary < firstJetindex || lastJetindex < primary)) continue;
        if(fMCType.Contains("MBMC") && (primary < firstUEindex  || lastUEindex  < primary)) continue;
      }

      pT = p->Pt();
      rapidity = p->Y();
      phi = p->Phi();
      pdg = p->GetPdgCode();

      //rapidity is Y(), but, pseudo-rapidity is Eta();

      if(pT < 1e-3) continue;//reject below 1 MeV
      if(TMath::Abs(rapidity) > 0.5) continue;

      if(p->Rho() > 1.0) continue;
      weight = 1.;

      if(pdg==111){//pi0
        parname = "Pi0";
        if(IsFrom(i,TrueK0SPt,310)) weight = f1K0SWeight->Eval(TrueK0SPt);
        else                        weight = f1Pi0Weight->Eval(pT);
      }
      else if(pdg==221){//eta
        parname = "Eta";
      }
      else if(pdg==22){//gamma
        parname = "Gamma";
        if(IsFrom(i,TruePi0Pt,111)){
          if(IsFrom(i,TrueK0SPt,310)) weight = f1K0SWeight->Eval(TrueK0SPt);
          else                        weight = f1Pi0Weight->Eval(TruePi0Pt);
        }
      }
      else if(pdg==211 || pdg==-211){//pi+ or pi-
        //c x tau = 7.8m
        parname = "ChargedPion";
      }
      else if(pdg==321 || pdg==-321){//K+ or K-
        //c x tau = 3.7m
        parname = "ChargedKaon";
        weight = f1K0SWeight->Eval(pT);
      }
      else if(pdg==310){//K0S
        parname = "K0S";
        weight = f1K0SWeight->Eval(pT);
      }
      else if(pdg==130){//K0L
        parname = "K0L";
        weight = f1K0SWeight->Eval(pT);
      }
      else if(pdg==3122){//Lmabda0
        parname = "Lambda0";
      }
      else if(pdg==3212){//Sigma0
        parname = "Sigma0";
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

      //rapidity is Y(), but, pseudo-rapidity is Eta();

      if(pT < 1e-3) continue;//reject below 1 MeV
      if(TMath::Abs(rapidity) > 0.5) continue;

      if(Rho(p) > 1.0) continue;
      weight = 1.;

      if(pdg==111){//pi0
        parname = "Pi0";
        if(IsFrom(i,TrueK0SPt,310)) weight = f1K0SWeight->Eval(TrueK0SPt);
        else                        weight = f1Pi0Weight->Eval(pT);
      }
      else if(pdg==221){//eta
        parname = "Eta";
      }
      else if(pdg==22){//gamma
        parname = "Gamma";
        if(IsFrom(i,TruePi0Pt,111)){
          if(IsFrom(i,TrueK0SPt,310)) weight = f1K0SWeight->Eval(TrueK0SPt);
          else                        weight = f1Pi0Weight->Eval(TruePi0Pt);
        }

      }
      else if(pdg==211 || pdg==-211){//pi+ or pi-
        parname = "ChargedPion";
      }
      else if(pdg==321 || pdg==-321){//K+ or K-
        parname = "ChargedKaon";
        weight = f1K0SWeight->Eval(pT);
      }
      else if(pdg==310){//K0S
        parname = "K0S";
        weight = f1K0SWeight->Eval(pT);
      }
      else if(pdg==130){//K0L
        parname = "K0L";
        weight = f1K0SWeight->Eval(pT);
      }
      else if(pdg==3122){//Lmabda0
        parname = "Lambda0";
      }
      else if(pdg==3212){//Sigma0
        parname = "Sigma0";
      }
      else{
        continue;
      }

      FillHistogramTH1(fOutputContainer,Form("hGen%sPt"    ,parname.Data()),pT          ,weight);
      FillHistogramTH2(fOutputContainer,Form("hGen%sEtaPhi",parname.Data()),phi,rapidity,weight);
      FillHistogramTH2(fOutputContainer,Form("hGen%sEtaPt" ,parname.Data()),rapidity,pT ,weight);

    }//end of generated particle loop


  }


}
//________________________________________________________________________
Double_t AliAnalysisTaskPHOSPi0EtaToGammaGamma::R(AliAODMCParticle* p)
{
  //Radius of vertex in cylindrical system.

//  if(p->PdgCode() == 111){
//    cout << "Xv = " << p->Xv() << " , fVertex[0] = " << fVertex[0] << endl;
//    cout << "Yv = " << p->Yv() << " , fVertex[1] = " << fVertex[1] << endl;
//    cout << "Zv = " << p->Zv() << " , fVertex[2] = " << fVertex[2] << endl;
//  }

  Double32_t x = p->Xv() - fVertex[0];
  Double32_t y = p->Yv() - fVertex[1];
  return sqrt(x*x + y*y);

  //Double32_t x = p->Xv();
  //Double32_t y = p->Yv();
  //Double32_t z = p->Zv();
  //return sqrt(x*x + y*y + z*z);
}
//________________________________________________________________________
Double_t AliAnalysisTaskPHOSPi0EtaToGammaGamma::Rho(AliAODMCParticle* p)
{
  //Radius of vertex in spherical system.

//  if(p->PdgCode() == 111){
//    cout << "Xv = " << p->Xv() << " , fVertex[0] = " << fVertex[0] << endl;
//    cout << "Yv = " << p->Yv() << " , fVertex[1] = " << fVertex[1] << endl;
//    cout << "Zv = " << p->Zv() << " , fVertex[2] = " << fVertex[2] << endl;
//  }

  Double32_t x = p->Xv() - fVertex[0];
  Double32_t y = p->Yv() - fVertex[1];
  Double32_t z = p->Zv() - fVertex[2];
  return sqrt(x*x + y*y + z*z);

  //Double32_t x = p->Xv();
  //Double32_t y = p->Yv();
  //Double32_t z = p->Zv();
  //return sqrt(x*x + y*y + z*z);
}
//________________________________________________________________________
Double_t AliAnalysisTaskPHOSPi0EtaToGammaGamma::DeltaPhiIn0Pi(Double_t dphi)
{
  //this returns dphi in 0-pi range.
  const Int_t myHarmonic = 2;
  Double_t tmp = dphi;
  while(tmp < 0)           tmp += 2./(Double_t)myHarmonic * TMath::Pi();
  while(tmp > TMath::Pi()) tmp -= 2./(Double_t)myHarmonic * TMath::Pi();

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
    TParticle *p = (TParticle*)fMCArrayESD->Particle(label);
    TParticle *mp = 0x0;
    motherid = p->GetFirstMother();
    while(motherid > -1){
      mp = (TParticle*)fMCArrayESD->Particle(motherid);
      pT = mp->Pt();
      pdg = mp->GetPdgCode(); 

      if(pdg == target_pdg && mp->Rho() < 1.0){//pi0 from primary vertex
        TruePt = pT;
        return kTRUE;
      }

      motherid = p->GetFirstMother();
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

      if(pdg == target_pdg && Rho(mp) < 1.0){//pi0 from primary vertex
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

  AliAODMCParticle *p = (AliAODMCParticle*)fMCArrayAOD->At(label);
  Int_t pdg = p->PdgCode(); 

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
    Int_t iprim2=jPart;

    while(iprim2>-1){
      if(iprim1==iprim2) return iprim1;
      iprim2 = dynamic_cast<AliAODMCParticle*>(fMCArrayAOD->At(iprim2))->GetMother();
    }

    iprim1 = dynamic_cast<AliAODMCParticle*>(fMCArrayAOD->At(iprim1))->GetMother();
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
  TF1 *f1Pi0Weight = (TF1*)GetAdditionalPi0PtWeightFunction(fCentralityMain);
  TF1 *f1K0SWeight = (TF1*)GetAdditionalK0SPtWeightFunction(fCentralityMain);
  Int_t primary = -1;
  Double_t TruePi0Pt = 1.;
  Double_t TrueK0SPt = 0;
  Double_t TrueL0Pt = 0;

  for(Int_t iph=0;iph<multClust;iph++){
    AliCaloPhoton *ph = (AliCaloPhoton*)fPHOSClusterArray->At(iph);
    primary = ph->GetPrimary();
    weight = 1.;

    if(IsFrom(primary,TruePi0Pt,111)){//pi0
      //for feed down correction from K0S->pi0 + pi0
      if(IsFrom(primary,TrueK0SPt,310)) weight *= f1K0SWeight->Eval(TrueK0SPt);
      else                              weight *= f1Pi0Weight->Eval(TruePi0Pt);
    }

    if(IsFrom(primary,TrueL0Pt,3122)){//lambda0
      //for feed down correction from L0->pi0 + neutron
      //weight *= f1K0SWeight->Eval(TrueK0SPt);
      weight *= 1.;
    }

    ph->SetWeight(weight);

  }//end of cluster loop

}
//_______________________________________________________________________________
void AliAnalysisTaskPHOSPi0EtaToGammaGamma::FillRejectionFactorMB()
{
  //search for 0PH0 event in MB for the rejection factor.
  //the best senario of how to estimate trigger rejection factor is reading ALICE electric log book.
  //but, fake trigger might be there.

  Bool_t Is0PH0fired = fEvent->GetHeader()->GetL0TriggerInputs() & 1 << (9 -1);//trigger input -1
  Bool_t Is1PHHfired = fEvent->GetHeader()->GetL1TriggerInputs() & 1 << (7 -1);//trigger input -1
  Bool_t Is1PHMfired = fEvent->GetHeader()->GetL1TriggerInputs() & 1 << (6 -1);//trigger input -1
  Bool_t Is1PHLfired = fEvent->GetHeader()->GetL1TriggerInputs() & 1 << (5 -1);//trigger input -1

  Bool_t Is0PH0matched = fPHOSTriggerHelperL0 ->IsPHI7(fEvent,fPHOSClusterCuts);
  Bool_t Is1PHHmatched = fPHOSTriggerHelperL1H->IsPHI7(fEvent,fPHOSClusterCuts);
  Bool_t Is1PHMmatched = fPHOSTriggerHelperL1M->IsPHI7(fEvent,fPHOSClusterCuts);
  Bool_t Is1PHLmatched = fPHOSTriggerHelperL1L->IsPHI7(fEvent,fPHOSClusterCuts);

  if(Is0PH0fired && Is0PH0matched) FillHistogramTH1(fOutputContainer,"hEventSummary", 7);//0PH0
  if(Is1PHLfired && Is1PHHmatched) FillHistogramTH1(fOutputContainer,"hEventSummary", 8);//1PHL
  if(Is1PHMfired && Is1PHMmatched) FillHistogramTH1(fOutputContainer,"hEventSummary", 9);//1PHM
  if(Is1PHHfired && Is1PHLmatched) FillHistogramTH1(fOutputContainer,"hEventSummary",10);//1PHH

}
//_______________________________________________________________________________
const AliQnCorrectionsQnVector *AliAnalysisTaskPHOSPi0EtaToGammaGamma::GetQnVectorFromList(const TList *list, const char* subdetector)
{
  AliQnCorrectionsQnVector *theQnVector = NULL;

  TList *pQvecList = dynamic_cast<TList*> (list->FindObject(subdetector));
  if (pQvecList != NULL){
    theQnVector = (AliQnCorrectionsQnVector*) pQvecList->First();//always get latest
    return theQnVector;
  }
  else return 0x0;

//  AliQnCorrectionsQnVector *theQnVector = NULL;
//
//  TList *pQvecList = dynamic_cast<TList*> (list->FindObject(subdetector));
//  if (pQvecList != NULL) {
//    /* the detector is present */
//    if (TString(expectedstep).EqualTo("latest"))
//      theQnVector = (AliQnCorrectionsQnVector*) pQvecList->First();
//    else
//      theQnVector = (AliQnCorrectionsQnVector*) pQvecList->FindObject(expectedstep);
//
//    if (theQnVector == NULL || !(theQnVector->IsGoodQuality()) || !(theQnVector->GetN() != 0)) {
//      /* the Qn vector for the expected step was not there */
//      if (TString(altstep).EqualTo("latest"))
//        theQnVector = (AliQnCorrectionsQnVector*) pQvecList->First();
//      else
//        theQnVector = (AliQnCorrectionsQnVector*) pQvecList->FindObject(altstep);
//    }
//  }
//  if (theQnVector != NULL) {
//    /* check the Qn vector quality */
//    if (!(theQnVector->IsGoodQuality()) || !(theQnVector->GetN() != 0))
//      /* not good quality, discarded */
//      theQnVector = NULL;
//  }
//  return theQnVector;



//  return 0x0;
}
//_______________________________________________________________________________

