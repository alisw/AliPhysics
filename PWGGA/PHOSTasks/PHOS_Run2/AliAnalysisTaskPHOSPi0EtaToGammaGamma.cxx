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
#include "TList.h"
#include "THashList.h"
#include "TMath.h"

#include "AliInputEventHandler.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskSE.h"
#include "AliLog.h"

#include "AliCaloPhoton.h"
#include "AliPHOSGeometry.h"

#include "AliESDtrackCuts.h"
#include "AliESDHeader.h"
#include "AliESDEvent.h"
#include "AliESDVZERO.h"
#include "AliESDTZERO.h"
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

#include "AliPID.h"
#include "AliPIDResponse.h"

#include "AliT0digit.h"
#include "AliAODMCHeader.h"
#include "AliAODHeader.h"
#include "AliAODEvent.h"
#include "AliAODVZERO.h"
#include "AliAODTZERO.h"
#include "AliAODCaloCells.h"
#include "AliAODCaloCluster.h"
#include "AliAODVertex.h"
#include "AliAODTrack.h"
#include "AliAODMCParticle.h"
#include "AliAODInputHandler.h"

#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliStack.h"

#include "AliPHOSCalibData.h"
#include "AliCDBManager.h"
#include "AliCDBEntry.h"
#include "AliCDBMetaData.h"
#include "AliPHOSAodCluster.h"
#include "AliPHOSEsdCluster.h"
#include "AliOADBContainer.h"
//#include "AliEventCuts.h"

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
	fUsePHOSTender(kFALSE),
  fUseCoreEnergy(kFALSE),
  fPHOSEventCuts(0x0),
  fPHOSClusterCuts(0x0),
  fBunchSpace(25.),
  fCentEdges(0),
  fCentNMixed(0),
  fCollisionSystem(-1),
  fTOFEfficiency(0x0),
  fESDtrackCutsGlobal(0x0),
  fESDtrackCutsGlobalConstrained(0x0),
  fAdditionalPi0PtWeight(0x0),
  fAdditionalK0SPtWeight(0x0),

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
  fCentralityV0M(0),
  fCentralityCL0(0),
  fCentralityCL1(0),
  fCentralityV0A(0),
  fCentralityV0C(0),
  fCentralityZNA(0),
  fCentralityZNC(0),
  fCentralityMain(0),
  fZvtx(-1),
  fCenBin(-1),
  fNHybridTrack(0),
  fIsNonLinStudyNeeded(kFALSE)
 // fEventCuts(0x0)
{
  // Constructor

  for(Int_t i=0;i<10;i++){
    for(Int_t j=0;j<10;j++){
      fPHOSEvents[i][j] = 0;
    }
  }

  for(Int_t i=0;i<3;i++){
    fVertex[i] = 0;
  }

  for(Int_t i=0;i<7;i++){
    for(Int_t j=0;j<7;j++){
      for(Int_t k=0;k<7;k++){
        fNonLin[i][j][k] = 0x0;
      }
    }
  }
  //fEventCuts = new AliEventCuts(kFALSE); 
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
    for(Int_t j=0;j<10;j++){
      if(fPHOSEvents[i][j]){
        delete fPHOSEvents[i][j];
        fPHOSEvents[i][j] = 0;
      }
    }
  }

  if(fIsNonLinStudyNeeded){
    for(Int_t i=0;i<7;i++){
      for(Int_t j=0;j<7;j++){
        for(Int_t k=0;k<7;k++){
          delete fNonLin[i][j][k];
          fNonLin[i][j][k] = 0;
        }
      }
    }
  }

  //delete fEventCuts;
  if(fJJMCHandler) delete fJJMCHandler;

}
//________________________________________________________________________
void AliAnalysisTaskPHOSPi0EtaToGammaGamma::UserCreateOutputObjects()
{
  // Create histograms
  // Called once

  fOutputContainer = new THashList();
  fOutputContainer->SetOwner(kTRUE);

  //fEventCuts.AddQAplotsToList(fOutputContainer); /// fList is your output TList

  TH1F *hEventSummary = new TH1F("hEventSummary","Event Summary",10,0.5,10.5);
  hEventSummary->GetXaxis()->SetBinLabel(1,"all");
  hEventSummary->GetXaxis()->SetBinLabel(2,"selected");
  hEventSummary->GetXaxis()->SetBinLabel(3,"kINT7");
  hEventSummary->GetXaxis()->SetBinLabel(4,"kPHI7 0PH0");
  hEventSummary->GetXaxis()->SetBinLabel(5,"kPHI7 1PHL");
  hEventSummary->GetXaxis()->SetBinLabel(6,"kPHI7 1PHM");
  hEventSummary->GetXaxis()->SetBinLabel(7,"kPHI7 1PHH");
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

  fOutputContainer->Add(new TH2F("hCentralityvsPHOSClusterMultiplicity"      ,"Centrality vs. Cluster Multiplicity"                ,100,0,100,101,-0.5,100.5));
  fOutputContainer->Add(new TH2F("hCentralityvsPHOSClusterMultiplicityTOF"   ,"Centrality vs. Cluster Multiplicity with TOF cut"   ,100,0,100,101,-0.5,100.5));
  fOutputContainer->Add(new TH2F("hCentralityV0MvsPHOSClusterMultiplicity"   ,"CentralityV0M vs. Cluster Multiplicity"             ,100,0,100,101,-0.5,100.5));
  fOutputContainer->Add(new TH2F("hCentralityV0MvsPHOSClusterMultiplicityTOF","CentralityV0M vs. Cluster Multiplicity with TOF cut",100,0,100,101,-0.5,100.5));
  fOutputContainer->Add(new TH2F("hCentralityZNAvsPHOSClusterMultiplicity"   ,"CentralityZNA vs. Cluster Multiplicity"             ,100,0,100,101,-0.5,100.5));
  fOutputContainer->Add(new TH2F("hCentralityZNAvsPHOSClusterMultiplicityTOF","CentralityZNA vs. Cluster Multiplicity with TOF cut",100,0,100,101,-0.5,100.5));
  fOutputContainer->Add(new TH2F("hCentralityZNCvsPHOSClusterMultiplicity"   ,"CentralityZNC vs. Cluster Multiplicity"             ,100,0,100,101,-0.5,100.5));
  fOutputContainer->Add(new TH2F("hCentralityZNCvsPHOSClusterMultiplicityTOF","CentralityZNC vs. Cluster Multiplicity with TOF cut",100,0,100,101,-0.5,100.5));

  fOutputContainer->Add(new TH2F("hCentralityvsTrackMultiplicity"   ,"Centrality vs. track Multiplicity"   ,100,0,100,5000,0,5000));
  fOutputContainer->Add(new TH2F("hCentralityV0MvsTrackMultiplicity","CentralityV0M vs. track Multiplicity",100,0,100,5000,0,5000));
  fOutputContainer->Add(new TH2F("hCentralityZNAvsTrackMultiplicity","CentralityZNA vs. track Multiplicity",100,0,100,5000,0,5000));
  fOutputContainer->Add(new TH2F("hCentralityZNCvsTrackMultiplicity","CentralityZNC vs. track Multiplicity",100,0,100,5000,0,5000));
  fOutputContainer->Add(new TH2F("hPHOSClusterMultiplicityvsTrackMultiplicity"   ,"cluster multiplicity vs. track multiplicity"         ,500,0,5000,101,-0.5,100.5));
  fOutputContainer->Add(new TH2F("hPHOSClusterMultiplicityTOFvsTrackMultiplicity","cluster multiplicity with TOF vs. track multiplicity",500,0,5000,101,-0.5,100.5));

  //track QA histograms
  const Int_t Ntype=3;
  const TString tracktype[Ntype] = {"Hybrid","Global","Complementary"};

  for(Int_t itype=0;itype<Ntype;itype++){
    fOutputContainer->Add(new TH1F(Form("h%sTrackMult",tracktype[itype].Data())  ,Form("Number of %s track",tracktype[itype].Data()),500,0,5000));
    fOutputContainer->Add(new TH1F(Form("h%sTrackPt",tracktype[itype].Data())    ,Form("%s track p_{T}",tracktype[itype].Data()),100,0,100));
    fOutputContainer->Add(new TH2F(Form("h%sTrackDCA",tracktype[itype].Data())   ,Form("%s track DCA_{xy} vs. DCA_{z}",tracktype[itype].Data()),100,-5,5,100,-5,5));//in cm.
    fOutputContainer->Add(new TH2F(Form("h%sTrackEtaPhi",tracktype[itype].Data()),Form("%s track #eta vs. #phi",tracktype[itype].Data()),100,0,6.3,200,-1,1));
  }
  fOutputContainer->Add(new TH2F("hTrackTPCdEdx","TPC dE/dx vs. track momentum",200,0,20,200,0,200));

  //cell QA histograms
  const Int_t Nmod=5;
  for(Int_t imod=1;imod<Nmod;imod++) fOutputContainer->Add(new TH2F(Form("hCellNXZM%d",imod),Form("Cell N(X,Z) M%d",imod),64,0.5,64.5,56,0.5,56.5));
  for(Int_t imod=1;imod<Nmod;imod++) fOutputContainer->Add(new TH2F(Form("hCellEXZM%d",imod),Form("Cell E(X,Z) M%d",imod),64,0.5,64.5,56,0.5,56.5));
  for(Int_t imod=1;imod<Nmod;imod++) fOutputContainer->Add(new TH2F(Form("hCellAmpTimeM%d_LG",imod),Form("Cell Amplitude vs. Time LG M%d",imod),200,0,20,1000,-500,500));
  for(Int_t imod=1;imod<Nmod;imod++) fOutputContainer->Add(new TH2F(Form("hCellAmpTimeM%d_HG",imod),Form("Cell Amplitude vs. Time HG M%d",imod),200,0,20,1000,-500,500));
  for(Int_t imod=1;imod<Nmod;imod++) fOutputContainer->Add(new TH1F(Form("hCellMultEventM%d",imod),Form("PHOS cell multiplicity per event M%d",imod),1001,-0.5,1000.5));

  //cluster QA histograms
  for(Int_t imod=1;imod<Nmod;imod++) fOutputContainer->Add(new TH1F(Form("hPHOSClusterMultM%d",imod)   ,Form("PHOS cluster multiplicity M%d",imod)             ,201,-0.5,200.5));
  for(Int_t imod=1;imod<Nmod;imod++) fOutputContainer->Add(new TH1F(Form("hPHOSClusterMultTOFM%d",imod),Form("PHOS cluster multiplicity M%d with TOF cut",imod),201,-0.5,200.5));
  for(Int_t imod=1;imod<Nmod;imod++) fOutputContainer->Add(new TH1F(Form("hPHOSClusterEnergyM%d",imod)   ,Form("PHOS cluster energy M%d",imod)             ,1000,0.,100));
  for(Int_t imod=1;imod<Nmod;imod++) fOutputContainer->Add(new TH1F(Form("hPHOSClusterEnergyTOFM%d",imod),Form("PHOS cluster energy M%d with TOF cut",imod),1000,0.,100));
  for(Int_t imod=1;imod<Nmod;imod++) fOutputContainer->Add(new TH2F(Form("hCluNXZM%d",imod)    ,Form("Cluster N(X,Z) M%d",imod)     ,64,0.5,64.5, 56,0.5,56.5));
  for(Int_t imod=1;imod<Nmod;imod++) fOutputContainer->Add(new TH2F(Form("hCluEXZM%d",imod)    ,Form("Cluster E(X,Z) M%d",imod)     ,64,0.5,64.5, 56,0.5,56.5));
  for(Int_t imod=1;imod<Nmod;imod++) fOutputContainer->Add(new TH2F(Form("hCluNXZTOFM%d",imod) ,Form("Cluster N(X,Z) M%d TOF",imod) ,64,0.5,64.5, 56,0.5,56.5));
  for(Int_t imod=1;imod<Nmod;imod++) fOutputContainer->Add(new TH2F(Form("hCluEXZTOFM%d",imod) ,Form("Cluster E(X,Z) M%d TOF",imod) ,64,0.5,64.5, 56,0.5,56.5));
  for(Int_t imod=1;imod<Nmod;imod++) fOutputContainer->Add(new TH2F(Form("hCluLowNXZM%d",imod) ,Form("Cluster Low N(X,Z) M%d",imod) ,64,0.5,64.5, 56,0.5,56.5));
  for(Int_t imod=1;imod<Nmod;imod++) fOutputContainer->Add(new TH2F(Form("hCluLowEXZM%d",imod) ,Form("Cluster Low E(X,Z) M%d",imod) ,64,0.5,64.5, 56,0.5,56.5));
  for(Int_t imod=1;imod<Nmod;imod++) fOutputContainer->Add(new TH2F(Form("hCluHighNXZM%d",imod),Form("Cluster High N(X,Z) M%d",imod),64,0.5,64.5, 56,0.5,56.5));
  for(Int_t imod=1;imod<Nmod;imod++) fOutputContainer->Add(new TH2F(Form("hCluHighEXZM%d",imod),Form("Cluster High E(X,Z) M%d",imod),64,0.5,64.5, 56,0.5,56.5));

  for(Int_t imod=1;imod<Nmod;imod++) fOutputContainer->Add(new TH2F(Form("hClusterEvsNM%d",imod),Form("Cluster E vs N_{cell} M%d",imod),500,0,50,100,0.5,100.5));
  for(Int_t imod=1;imod<Nmod;imod++) fOutputContainer->Add(new TH2F(Form("hClusterEvsTM%d",imod),Form("Cluster E vs TOF M%d",imod)     ,500,0,50, 1000,-500,500));
  for(Int_t imod=1;imod<Nmod;imod++) fOutputContainer->Add(new TH2F(Form("hClusterEvsM02M%d",imod),Form("Cluster E vs M02 M%d",imod),500,0,50,100,0,10));
  for(Int_t imod=1;imod<Nmod;imod++) fOutputContainer->Add(new TH2F(Form("hClusterNvsM02M%d",imod),Form("Cluster N vs M02 M%d",imod),100,0.5,100.5,100,0,10));
  for(Int_t imod=1;imod<Nmod;imod++) fOutputContainer->Add(new TH2F(Form("hFullDispvsFullEM%d",imod),Form("full dispersion vs full E M%d",imod),500,0,50,100,0,10));
  for(Int_t imod=1;imod<Nmod;imod++) fOutputContainer->Add(new TH2F(Form("hCoreDispvsCoreEM%d",imod),Form("core dispersion vs core E M%d",imod),500,0,50,100,0,10));
  for(Int_t imod=1;imod<Nmod;imod++) fOutputContainer->Add(new TH2F(Form("hFullDispvsCoreEM%d",imod),Form("full dispersion vs full E M%d",imod),500,0,50,100,0,10));
  for(Int_t imod=1;imod<Nmod;imod++) fOutputContainer->Add(new TH2F(Form("hCoreDispvsFullEM%d",imod),Form("core dispersion vs core E M%d",imod),500,0,50,100,0,10));
  for(Int_t imod=1;imod<Nmod;imod++) fOutputContainer->Add(new TH2F(Form("hRvsTrackPtM%d",imod),Form("r vs track pT M%d",imod),100,0,10,100,0,50));

  fOutputContainer->Add(new TH2F("hClusterEtaPhi","Cluster eta vs.phi",200,-1,1,100,0,6.3));
  fOutputContainer->Add(new TH2F("hEnergyvsDistanceToBadChannel","distance to closest bad channel",100,0,50,200,0,20));
  fOutputContainer->Add(new TH3F("hEvsNvsM02","energy vs. N vs. M02",100,0,50,100,0.5,100.5,100,0,10));

  fOutputContainer->Add(new TH2F("hAsymvsMgg","asymmetry vs. M_{#gamma#gamma}",100,0,1,60,0,0.24));
  fOutputContainer->Add(new TH2F("hAsymvsPt" ,"asymmetry vs. p_{T}^{#gamma#gamma}",100,0,1,500,0,50));

  const Int_t Ncen = fCentEdges.GetSize();

  Int_t CenMax=0;
  if(fCollisionSystem==0)      CenMax = Ncen;//for pp
  else if(fCollisionSystem==1) CenMax = Ncen-1;//for PbPb
  else if(fCollisionSystem==2) CenMax = Ncen-1;//for pPb

  Int_t nM = 240;
  Double_t Mmin = 0.;
  Double_t Mmax = 0.96;

  const Int_t NpTgg = 101;
  Double_t pTgg[NpTgg]={};

  for(Int_t i=0;i<50;i++)     pTgg[i] = 0.1 * i;            //every 0.1 GeV/c, up to 5 GeV/c
  for(Int_t i=50;i<60;i++)    pTgg[i] = 0.5 * (i-50) + 5.0; //every 0.5 GeV/c, up to 10 GeV/c
  for(Int_t i=60;i<NpTgg;i++) pTgg[i] = 1.0 * (i-60) + 10.0;//every 1.0 GeV/c, up to 50 GeV/c

  //In order to support pp and PbPb in this same code, name for centrality binning is just 0,1,2,3...

  //photon pT histograms
  for(Int_t icen=0;icen<CenMax;icen++){
    TH1F *h1 = new TH1F(Form("hPhotonPt_Cen%d",icen),Form("Photon pT cenbin %d",icen),NpTgg-1,pTgg);
    h1->Sumw2();
    fOutputContainer->Add(h1);
  }

  for(Int_t icen=0;icen<CenMax;icen++){
    TH1F *h1 = new TH1F(Form("hPhotonPt_Cen%d_TOF",icen),Form("Photon pT with TOF cut cenbin %d",icen),NpTgg-1,pTgg);
    h1->Sumw2();
    fOutputContainer->Add(h1);
  }

  //const TString Asym[] = {"","_asym08","_asym07"};
  const TString Asym[] = {"","_asym08"};
  const Int_t Nasym = sizeof(Asym)/sizeof(Asym[0]);

  for(Int_t iasym=0;iasym<Nasym;iasym++){

    //Mgg vs. pT histogram
    for(Int_t icen=0;icen<CenMax;icen++){
      TH2F *h2 = new TH2F(Form("hMgg_Cen%d%s",icen,Asym[iasym].Data()),"M_{#gamma#gamma} vs p_{T}",240,0,0.96,NpTgg-1,pTgg);
      h2->Sumw2();
      fOutputContainer->Add(h2);
    }
    for(Int_t icen=0;icen<CenMax;icen++){
      TH2F *h2 = new TH2F(Form("hMgg_Cen%d_TOF%s",icen,Asym[iasym].Data()),"M_{#gamma#gamma} vs p_{T}",240,0,0.96,NpTgg-1,pTgg);
      h2->Sumw2();
      fOutputContainer->Add(h2);
    }//end of centrality bin loop

    //Mix Mgg vs. pT histogram
    for(Int_t icen=0;icen<CenMax;icen++){
      TH2F *h2 = new TH2F(Form("hMixMgg_Cen%d%s",icen,Asym[iasym].Data()),"M_{#gamma#gamma}^{mix} vs p_{T}",240,0,0.96,NpTgg-1,pTgg);
      h2->Sumw2();
      fOutputContainer->Add(h2);
    }
    for(Int_t icen=0;icen<CenMax;icen++){
      TH2F *h2 = new TH2F(Form("hMixMgg_Cen%d_TOF%s",icen,Asym[iasym].Data()),"M_{#gamma#gamma}^{mix} vs p_{T}",240,0,0.96,NpTgg-1,pTgg);
      h2->Sumw2();
      fOutputContainer->Add(h2);
    }//end of centrality bin loop

  }//end of asymmetry loop


  const Int_t NpTggModule = 71;
  Double_t pTggModule[NpTggModule]={};
  for(Int_t i=0;i<50;i++)     pTggModule[i] = 0.1 * i;            //every 0.1 GeV/c, up to 5 GeV/c
  for(Int_t i=50;i<60;i++)    pTggModule[i] = 0.5 * (i-50) + 5.0; //every 0.5 GeV/c, up to 10 GeV/c
  for(Int_t i=60;i<NpTggModule;i++) pTggModule[i] = 1.0 * (i-60) + 10.0;//every 1.0 GeV/c, up to 20 GeV/c

//  for(Int_t icen=0;icen<CenMax;icen++){
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

//  }//end of centrality loop

  //for PID cut efficiency
  for(Int_t icen=0;icen<CenMax;icen++){
    TH2F *h2 = new TH2F(Form("hMgg_Probe_PID_Cen%d"         ,icen),"M_{#gamma#gamma} vs E_{#gamma}^{probe}"             ,60,0,0.24,NpTgg-1,pTgg);
    h2->Sumw2();
    fOutputContainer->Add(h2);
  }
  for(Int_t icen=0;icen<CenMax;icen++){
    TH2F *h2 = new TH2F(Form("hMgg_PassingProbe_PID_Cen%d"   ,icen),"M_{#gamma#gamma} vs E_{#gamma}^{passing probe}"      ,60,0,0.24,NpTgg-1,pTgg);
    h2->Sumw2();
    fOutputContainer->Add(h2);
  }
  for(Int_t icen=0;icen<CenMax;icen++){
    TH2F *h2 = new TH2F(Form("hMixMgg_Probe_PID_Cen%d"      ,icen),"M_{#gamma#gamma}^{mix} vs E_{#gamma}^{probe}"       ,60,0,0.24,NpTgg-1,pTgg);
    h2->Sumw2();
    fOutputContainer->Add(h2);
  }
  for(Int_t icen=0;icen<CenMax;icen++){
    TH2F *h2 = new TH2F(Form("hMixMgg_PassingProbe_PID_Cen%d",icen),"M_{#gamma#gamma}^{mix} vs E_{#gamma}^{passing probe}",60,0,0.24,NpTgg-1,pTgg);
    h2->Sumw2();
    fOutputContainer->Add(h2);
  }

/*
  //for Ncell cut efficiency
  for(Int_t icen=0;icen<CenMax;icen++){
    TH2F *h2 = new TH2F(Form("hMgg_Probe_Ncell_Cen%d"         ,icen),"M_{#gamma#gamma} vs E_{#gamma}^{probe}"             ,60,0,0.24,NpTgg-1,pTgg);
    h2->Sumw2();
    fOutputContainer->Add(h2);
  }
  for(Int_t icen=0;icen<CenMax;icen++){
    TH2F *h2 = new TH2F(Form("hMgg_PassingProbe_Ncell_Cen%d"   ,icen),"M_{#gamma#gamma} vs E_{#gamma}^{passing probe}"      ,60,0,0.24,NpTgg-1,pTgg);
    h2->Sumw2();
    fOutputContainer->Add(h2);
  }
  for(Int_t icen=0;icen<CenMax;icen++){
    TH2F *h2 = new TH2F(Form("hMixMgg_Probe_Ncell_Cen%d"      ,icen),"M_{#gamma#gamma}^{mix} vs E_{#gamma}^{probe}"       ,60,0,0.24,NpTgg-1,pTgg);
    h2->Sumw2();
    fOutputContainer->Add(h2);
  }
  for(Int_t icen=0;icen<CenMax;icen++){
    TH2F *h2 = new TH2F(Form("hMixMgg_PassingProbe_Ncell_Cen%d",icen),"M_{#gamma#gamma}^{mix} vs E_{#gamma}^{passing probe}",60,0,0.24,NpTgg-1,pTgg);
    h2->Sumw2();
    fOutputContainer->Add(h2);
  }

  //for M20 cut efficiency
  for(Int_t icen=0;icen<CenMax;icen++){
    TH2F *h2 = new TH2F(Form("hMgg_Probe_M20_Cen%d"         ,icen),"M_{#gamma#gamma} vs E_{#gamma}^{probe}"             ,60,0,0.24,NpTgg-1,pTgg);
    h2->Sumw2();
    fOutputContainer->Add(h2);
  }
  for(Int_t icen=0;icen<CenMax;icen++){
    TH2F *h2 = new TH2F(Form("hMgg_PassingProbe_M20_Cen%d"   ,icen),"M_{#gamma#gamma} vs E_{#gamma}^{passing probe}"      ,60,0,0.24,NpTgg-1,pTgg);
    h2->Sumw2();
    fOutputContainer->Add(h2);
  }
  for(Int_t icen=0;icen<CenMax;icen++){
    TH2F *h2 = new TH2F(Form("hMixMgg_Probe_M20_Cen%d"      ,icen),"M_{#gamma#gamma}^{mix} vs E_{#gamma}^{probe}"       ,60,0,0.24,NpTgg-1,pTgg);
    h2->Sumw2();
    fOutputContainer->Add(h2);
  }
  for(Int_t icen=0;icen<CenMax;icen++){
    TH2F *h2 = new TH2F(Form("hMixMgg_PassingProbe_M20_Cen%d",icen),"M_{#gamma#gamma}^{mix} vs E_{#gamma}^{passing probe}",60,0,0.24,NpTgg-1,pTgg);
    h2->Sumw2();
    fOutputContainer->Add(h2);
  }

  //for M02 cut efficiency
  for(Int_t icen=0;icen<CenMax;icen++){
    TH2F *h2 = new TH2F(Form("hMgg_Probe_M02_Cen%d"         ,icen),"M_{#gamma#gamma} vs E_{#gamma}^{probe}"             ,60,0,0.24,NpTgg-1,pTgg);
    h2->Sumw2();
    fOutputContainer->Add(h2);
  }
  for(Int_t icen=0;icen<CenMax;icen++){
    TH2F *h2 = new TH2F(Form("hMgg_PassingProbe_M02_Cen%d"   ,icen),"M_{#gamma#gamma} vs E_{#gamma}^{passing probe}"      ,60,0,0.24,NpTgg-1,pTgg);
    h2->Sumw2();
    fOutputContainer->Add(h2);
  }
  for(Int_t icen=0;icen<CenMax;icen++){
    TH2F *h2 = new TH2F(Form("hMixMgg_Probe_M02_Cen%d"      ,icen),"M_{#gamma#gamma}^{mix} vs E_{#gamma}^{probe}"       ,60,0,0.24,NpTgg-1,pTgg);
    h2->Sumw2();
    fOutputContainer->Add(h2);
  }
  for(Int_t icen=0;icen<CenMax;icen++){
    TH2F *h2 = new TH2F(Form("hMixMgg_PassingProbe_M02_Cen%d",icen),"M_{#gamma#gamma}^{mix} vs E_{#gamma}^{passing probe}",60,0,0.24,NpTgg-1,pTgg);
    h2->Sumw2();
    fOutputContainer->Add(h2);
  }

  //for STD cut efficiency
  for(Int_t icen=0;icen<CenMax;icen++){
    TH2F *h2 = new TH2F(Form("hMgg_Probe_STDCut_Cen%d"         ,icen),"M_{#gamma#gamma} vs E_{#gamma}^{probe}"             ,60,0,0.24,NpTgg-1,pTgg);
    h2->Sumw2();
    fOutputContainer->Add(h2);
  }
  for(Int_t icen=0;icen<CenMax;icen++){
    TH2F *h2 = new TH2F(Form("hMgg_PassingProbe_STDCut_Cen%d"   ,icen),"M_{#gamma#gamma} vs E_{#gamma}^{passing probe}"      ,60,0,0.24,NpTgg-1,pTgg);
    h2->Sumw2();
    fOutputContainer->Add(h2);
  }
  for(Int_t icen=0;icen<CenMax;icen++){
    TH2F *h2 = new TH2F(Form("hMixMgg_Probe_STDCut_Cen%d"      ,icen),"M_{#gamma#gamma}^{mix} vs E_{#gamma}^{probe}"       ,60,0,0.24,NpTgg-1,pTgg);
    h2->Sumw2();
    fOutputContainer->Add(h2);
  }
  for(Int_t icen=0;icen<CenMax;icen++){
    TH2F *h2 = new TH2F(Form("hMixMgg_PassingProbe_STDCut_Cen%d",icen),"M_{#gamma#gamma}^{mix} vs E_{#gamma}^{passing probe}",60,0,0.24,NpTgg-1,pTgg);
    h2->Sumw2();
    fOutputContainer->Add(h2);
  }
*/

  //for TOF cut efficiency
  for(Int_t icen=0;icen<CenMax;icen++){
    TH2F *h2 = new TH2F(Form("hMgg_Probe_TOF_Cen%d"         ,icen),"M_{#gamma#gamma} vs E_{#gamma}^{probe}"             ,60,0,0.24,NpTgg-1,pTgg);
    h2->Sumw2();
    fOutputContainer->Add(h2);
  }
  for(Int_t icen=0;icen<CenMax;icen++){
    TH2F *h2 = new TH2F(Form("hMgg_PassingProbe_TOF_Cen%d"   ,icen),"M_{#gamma#gamma} vs E_{#gamma}^{passing probe}"      ,60,0,0.24,NpTgg-1,pTgg);
    h2->Sumw2();
    fOutputContainer->Add(h2);
  }
  for(Int_t icen=0;icen<CenMax;icen++){
    TH2F *h2 = new TH2F(Form("hMixMgg_Probe_TOF_Cen%d"      ,icen),"M_{#gamma#gamma}^{mix} vs E_{#gamma}^{probe}"       ,60,0,0.24,NpTgg-1,pTgg);
    h2->Sumw2();
    fOutputContainer->Add(h2);
  }
  for(Int_t icen=0;icen<CenMax;icen++){
    TH2F *h2 = new TH2F(Form("hMixMgg_PassingProbe_TOF_Cen%d",icen),"M_{#gamma#gamma}^{mix} vs E_{#gamma}^{passing probe}",60,0,0.24,NpTgg-1,pTgg);
    h2->Sumw2();
    fOutputContainer->Add(h2);
  }

  Bool_t IsPHOSTriggerAnalysis = fPHOSEventCuts->IsPHOSTriggerAnalysis();

  if(IsPHOSTriggerAnalysis){
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
        fOutputContainer->Add(new TH3F(Form("hDistanceXZM%dTRU%d",imod,itru)        ,Form("distance TRUch-clu M%d TRU%d",imod,itru),41,-20.5,+20.5,41,-20.5,+20.5,100,0,50));
        fOutputContainer->Add(new TH3F(Form("hMatchedDistanceXZM%dTRU%d",imod,itru) ,Form("distance TRUch-clu M%d TRU%d",imod,itru),41,-20.5,+20.5,41,-20.5,+20.5,100,0,50));
      }
    }

    for(Int_t icen=0;icen<CenMax;icen++){
      for(Int_t imod=1;imod<Nmod;imod++){
        for(Int_t itru=1;itru<=Ntru;itru++){
          fOutputContainer->Add(new TH1F(Form("hClusterE_Cen%d_M%d_TRU%d",icen,imod,itru)         ,Form("cluster energy M%d TRU%d",imod,itru)          ,NpTgg-1,pTgg));
          fOutputContainer->Add(new TH1F(Form("hTriggeredClusterE_Cen%d_M%d_TRU%d",icen,imod,itru),Form("triggered cluster energy M%d TRU%d",imod,itru),NpTgg-1,pTgg));
        }
      }
    }

    //for trigger efficiency
    for(Int_t icen=0;icen<CenMax;icen++){
      TH2F *h2 = new TH2F(Form("hMgg_Probe_Trg_Cen%d"         ,icen),"M_{#gamma#gamma} vs E_{#gamma}^{probe}"             ,60,0,0.24,NpTgg-1,pTgg);
      h2->Sumw2();
      fOutputContainer->Add(h2);
    }
    for(Int_t icen=0;icen<CenMax;icen++){
      TH2F *h2 = new TH2F(Form("hMgg_PassingProbe_Trg_Cen%d"   ,icen),"M_{#gamma#gamma} vs E_{#gamma}^{passing probe}"      ,60,0,0.24,NpTgg-1,pTgg);
      h2->Sumw2();
      fOutputContainer->Add(h2);
    }
    for(Int_t icen=0;icen<CenMax;icen++){
      TH2F *h2 = new TH2F(Form("hMixMgg_Probe_Trg_Cen%d"      ,icen),"M_{#gamma#gamma}^{mix} vs E_{#gamma}^{probe}"       ,60,0,0.24,NpTgg-1,pTgg);
      h2->Sumw2();
      fOutputContainer->Add(h2);
    }
    for(Int_t icen=0;icen<CenMax;icen++){
      TH2F *h2 = new TH2F(Form("hMixMgg_PassingProbe_Trg_Cen%d",icen),"M_{#gamma#gamma}^{mix} vs E_{#gamma}^{passing probe}",60,0,0.24,NpTgg-1,pTgg);
      h2->Sumw2();
      fOutputContainer->Add(h2);
    }

    for(Int_t icen=0;icen<CenMax;icen++){
      for(Int_t imod=1;imod<Nmod;imod++){
        for(Int_t itru=1;itru<=Ntru;itru++){
          TH2F *h2 = new TH2F(Form("hMgg_Probe_Trg_Cen%d_M%d_TRU%d"         ,icen,imod,itru),"M_{#gamma#gamma} vs E_{#gamma}^{probe}",60,0,0.24,NpTgg-1,pTgg);
          h2->Sumw2();
          fOutputContainer->Add(h2);
        }
      }
    }

    for(Int_t icen=0;icen<CenMax;icen++){
      for(Int_t imod=1;imod<Nmod;imod++){
        for(Int_t itru=1;itru<=Ntru;itru++){
          TH2F *h2 = new TH2F(Form("hMgg_PassingProbe_Trg_Cen%d_M%d_TRU%d"         ,icen,imod,itru),"M_{#gamma#gamma} vs E_{#gamma}^{passing probe}",60,0,0.24,NpTgg-1,pTgg);
          h2->Sumw2();
          fOutputContainer->Add(h2);
        }
      }
    }

    for(Int_t icen=0;icen<CenMax;icen++){
      for(Int_t imod=1;imod<Nmod;imod++){
        for(Int_t itru=1;itru<=Ntru;itru++){
          TH2F *h2 = new TH2F(Form("hMixMgg_Probe_Trg_Cen%d_M%d_TRU%d"         ,icen,imod,itru),"M_{#gamma#gamma}^{mix} vs E_{#gamma}^{probe}",60,0,0.24,NpTgg-1,pTgg);
          h2->Sumw2();
          fOutputContainer->Add(h2);
        }
      }
    }

    for(Int_t icen=0;icen<CenMax;icen++){
      for(Int_t imod=1;imod<Nmod;imod++){
        for(Int_t itru=1;itru<=Ntru;itru++){
          TH2F *h2 = new TH2F(Form("hMixMgg_PassingProbe_Trg_Cen%d_M%d_TRU%d"         ,icen,imod,itru),"M_{#gamma#gamma}^{mix} vs E_{#gamma}^{passing probe}",60,0,0.24,NpTgg-1,pTgg);
          h2->Sumw2();
          fOutputContainer->Add(h2);
        }
      }
    }

  }//IsPHOSTriggerAnalysis ends

  if(fIsMC){
    //const Int_t Npar = 3;
    const TString parname[] = {"Pi0","Eta","Gamma","ChargedPion","ChargedKaon","K0S","K0L","Lambda0","Sigma0"};
    const Int_t Npar = sizeof(parname)/sizeof(parname[0]);

    for(Int_t icen=0;icen<CenMax;icen++){
      for(Int_t ipar=0;ipar<Npar;ipar++){
        TH1F *h1Pt = new TH1F(Form("hGen%sPt_Cen%d",parname[ipar].Data()        ,icen),Form("generated %s pT cenbin %d",parname[ipar].Data()        ,icen),NpTgg-1,pTgg);
        h1Pt->Sumw2();
        fOutputContainer->Add(h1Pt);

        TH2F *h2EtaPhi = new TH2F(Form("hGen%sEtaPhi_Cen%d",parname[ipar].Data(),icen),Form("generated %s eta vs phi cenbin %d",parname[ipar].Data(),icen),200,-1,1,100,0,6.3);
        h2EtaPhi->Sumw2();
        fOutputContainer->Add(h2EtaPhi);

        TH2F *h2EtaPt = new TH2F(Form("hGen%sEtaPt_Cen%d",parname[ipar].Data()  ,icen),Form("generated %s eta vs pT cenbin %d",parname[ipar].Data() ,icen),200,-1,1,NpTgg-1,pTgg);
        h2EtaPt->Sumw2();
        fOutputContainer->Add(h2EtaPt);

      }//end of particle loop
    }//end of centrality loop

    for(Int_t icen=0;icen<CenMax;icen++){
      TH1F *h1Pt = new TH1F(Form("hGenGammaFromPi0Pt_Cen%d",icen),Form("generated #gamma from #pi^{0} pT cenbin %d",icen),NpTgg-1,pTgg);
      h1Pt->Sumw2();
      fOutputContainer->Add(h1Pt);

      TH2F *h2EtaPhi = new TH2F(Form("hGenGammaFromPi0EtaPhi_Cen%d",icen),Form("generated #gamma from #pi^{0} eta vs phi cenbin %d",icen),200,-1,1,100,0,6.3);
      h2EtaPhi->Sumw2();
      fOutputContainer->Add(h2EtaPhi);

      TH2F *h2EtaPt = new TH2F(Form("hGenGammaFromPi0EtaPt_Cen%d",icen),Form("generated #gamma from #pi^{0} eta vs pT cenbin %d",icen),200,-1,1,NpTgg-1,pTgg);
      h2EtaPt->Sumw2();
      fOutputContainer->Add(h2EtaPt);

    }//end of centrality loop

    for(Int_t icen=0;icen<CenMax;icen++){
      TH1F *h1Pt = new TH1F(Form("hPurityGamma_Cen%d",icen),Form("pT of true photon in clusters for purity cenbin %d",icen),NpTgg-1,pTgg);
      h1Pt->Sumw2();
      fOutputContainer->Add(h1Pt);
    }

    if(fIsNonLinStudyNeeded){
      //for non-linearity study
      const Int_t Na=7;
      const Int_t Nb=7;

      for(Int_t icen=0;icen<CenMax;icen++){
        for(Int_t ia=0;ia<Na;ia++){
          for(Int_t ib=0;ib<Nb;ib++){
            Double_t a = -0.07 + 0.01*ia;
            Double_t b =  0.4  + 0.1 *ib;

            //fNonLin[ia][ib] = new TF1(Form("fNonLin_a%d_b%d"),"(1+[0]/(1+pow([1]*x,2.)))*1.007",0,100);
            fNonLin[icen][ia][ib] = new TF1(Form("fNonLin_Cen%d_a%d_b%d",icen,ia,ib),"1.00*(1.+[0]*TMath::Exp(-x*x/2./[1]/[1]))",0,100);
            fNonLin[icen][ia][ib]->SetParNames("a","b");
            fNonLin[icen][ia][ib]->FixParameter(0,a);
            fNonLin[icen][ia][ib]->FixParameter(1,b);

            fOutputContainer->Add(new TH2F(Form("hMgg_Cen%d_a%d_b%d",icen,ia,ib),Form("a = %3.2f , b = %3.2f",a,b),60,0,0.24,NpTggModule-1,pTggModule));
            fOutputContainer->Add(new TH2F(Form("hMixMgg_Cen%d_a%d_b%d",icen,ia,ib),Form("a = %3.2f , b = %3.2f",a,b),60,0,0.24,NpTggModule-1,pTggModule));
          }
        }
      }
    }

    for(Int_t icen=0;icen<CenMax;icen++){
      TH2F *h2 = new TH2F(Form("hPi0FromK0S_Cen%d",icen),"#pi^{0} from K_{S}^{0}",NpTgg-1,pTgg,500,0,500);
      h2->Sumw2();
      fOutputContainer->Add(h2);
    }

    for(Int_t iasym=0;iasym<Nasym;iasym++){
      for(Int_t icen=0;icen<CenMax;icen++){
        TH2F *h2 = new TH2F(Form("hMggFromK0S_Cen%d%s",icen,Asym[iasym].Data()),"M_{#gamma#gamma} from K_{S}^{0}",240,0,0.96,NpTgg-1,pTgg);
        h2->Sumw2();
        fOutputContainer->Add(h2);
      }//end of centrality loop

    }//end of asym loop

    for(Int_t iasym=0;iasym<Nasym;iasym++){
      for(Int_t icen=0;icen<CenMax;icen++){
        TH2F *h2 = new TH2F(Form("hMggFromLambda0_Cen%d%s",icen,Asym[iasym].Data()),"M_{#gamma#gamma} from K_{S}^{0}",240,0,0.96,NpTgg-1,pTgg);
        h2->Sumw2();
        fOutputContainer->Add(h2);
      }//end of centrality loop

    }//end of asym loop


    //for JJMC
    if(fIsJJMC){
      fOutputContainer->Add(new TH1F("hPtHard","pT hard in GeV/c",1000,0,1000));
      fOutputContainer->Add(new TH1F("hNTrial","nTrial",20,0.5,20.5));
      fOutputContainer->Add(new TProfile("hProfCrossSection","inelastic cross section",20,0.5,20.5));

      TH1F *hNMerged = new TH1F("hNMerged","N merged",20,0.5,20.5);
      hNMerged->Fill(fPtHardBin);
      hNMerged->SetYTitle("number of merged files");
      fOutputContainer->Add(hNMerged);
    }//end of IsJJMC

  }//end of IsMC

//  TH1F *h1test = new TH1F("h1test","h1test",1000,0,100);
//  h1test->Sumw2();
//  TF1 *f1 = GetTOFCutEfficiencyFunction();
//
//  for(Int_t i=1; i <= 1000; i++){
//    Double_t e = h1test->GetBinCenter(i);
//    Double_t eff = f1->Eval(e);
//    h1test->SetBinContent(i,eff);
//  }
//  fOutputContainer->Add(h1test);




  PostData(1,fOutputContainer);

}
//________________________________________________________________________
void AliAnalysisTaskPHOSPi0EtaToGammaGamma::UserExec(Option_t *option) 
{
  // Main loop
  // Called for each event

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

  //for JJMC
  if(fIsMC && fIsJJMC){
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

    const Int_t Ncen = fCentEdges.GetSize();

    if(fCentralityMain < (Float_t)fCentEdges[0] || (Float_t)fCentEdges[Ncen-1] < fCentralityMain){
      AliInfo(Form("Reject this event because centrality %s %f %% is out of the configuration of this task.", fEstimator.Data(),fCentralityMain));
      return;
    }

    fCenBin = GetCentralityBinPbPb(fCentralityMain);
    AliInfo(Form("fCentralityMain estimated by %s = %f %, fCenbin = %d , Zvtx = %f cm, fZvtx = %d.",fEstimator.Data(),fCentralityMain,fCenBin,fVertex[2],fZvtx));

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

    const Int_t Ncen = fCentEdges.GetSize();

    if(fCentralityMain < (Float_t)fCentEdges[0] || (Float_t)fCentEdges[Ncen-1] < fCentralityMain){
      AliInfo(Form("Reject this event because centrality %s %f %% is out of the configuration of this task.", fEstimator.Data(),fCentralityMain));
      return;
    }

    fCenBin = GetCentralityBinPbPb(fCentralityMain);
    AliInfo(Form("fCentralityMain estimated by %s = %f %, fCenbin = %d , Zvtx = %f cm, fZvtx = %d.",fEstimator.Data(),fCentralityMain,fCenBin,fVertex[2],fZvtx));

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

  //event selection
  if(!(fPHOSEventCuts->AcceptEvent(fEvent))){
    AliInfo("event is rejected.");
    return;
  }

  Bool_t IsPHOSTriggerAnalysis = fPHOSEventCuts->IsPHOSTriggerAnalysis();
  UInt_t fSelectMask = fInputHandler->IsEventSelected();
  Bool_t isINT7selected = fSelectMask & AliVEvent::kINT7;

  if(!IsPHOSTriggerAnalysis && !isINT7selected){
    AliInfo("INT7 Event is rejected by IsEventSelected()");
    return;
  }

  fPHOSClusterArray = (TClonesArray*)fEvent->FindListObject("PHOSClusterArray");
  if(!fPHOSClusterArray){
    //If you get this warning (and fCentralityV0M 300) please check that the AliMultSelectionTask actually ran (before your task)
    AliWarning("fPHOSClusterArray object not found!");
    return;
  }


  //Bool_t IsPHOSTriggerAnalysis = fPHOSEventCuts->IsPHOSTriggerAnalysis();
  Bool_t isPHI7selected = fSelectMask & AliVEvent::kPHI7;
  if(IsPHOSTriggerAnalysis && !isPHI7selected){
    AliInfo("PHI7 Event is rejected by IsEventSelected()");
    return;
  }

  //cout << "fCentralityV0M = " << fCentralityV0M << endl;

  //<- end of event selection
  //-> start physics analysis

  FillHistogramTH2(fOutputContainer,"hCentralityV0MvsCL0",fCentralityV0M,fCentralityCL0);
  FillHistogramTH2(fOutputContainer,"hCentralityV0MvsCL1",fCentralityV0M,fCentralityCL1);
  FillHistogramTH2(fOutputContainer,"hCentralityCL0vsCL1",fCentralityCL0,fCentralityCL1);
  FillHistogramTH2(fOutputContainer,"hCentralityV0AvsV0C",fCentralityV0A,fCentralityV0C);
  FillHistogramTH2(fOutputContainer,"hCentralityZNAvsZNC",fCentralityZNA,fCentralityZNC);

  FillHistogramTH1(fOutputContainer,"hVertexZSelectEvent" ,fVertex[2]);
  FillHistogramTH2(fOutputContainer,"hVertexXYSelectEvent",fVertex[0],fVertex[1]);
  FillHistogramTH1(fOutputContainer,"hEventSummary",2);//select event

  Bool_t IsPH0fired = fEvent->GetHeader()->GetL0TriggerInputs() & 1 << 8;//trigger input -1
  Bool_t IsPHHfired = fEvent->GetHeader()->GetL1TriggerInputs() & 1 << 6;//trigger input -1
  Bool_t IsPHMfired = fEvent->GetHeader()->GetL1TriggerInputs() & 1 << 5;//trigger input -1
  Bool_t IsPHLfired = fEvent->GetHeader()->GetL1TriggerInputs() & 1 << 4;//trigger input -1

  if(IsPH0fired) FillHistogramTH1(fOutputContainer,"hEventSummary",4);//0PH0
  if(IsPHLfired) FillHistogramTH1(fOutputContainer,"hEventSummary",5);//1PHL
  if(IsPHMfired) FillHistogramTH1(fOutputContainer,"hEventSummary",6);//1PHM
  if(IsPHHfired) FillHistogramTH1(fOutputContainer,"hEventSummary",7);//1PHH

  if(fRunNumber != fEvent->GetRunNumber()){ // Check run number
    fRunNumber = fEvent->GetRunNumber();
    fPHOSGeo = GetPHOSGeometry();
  }

  //track QA
  TrackQA();
  if(fCollisionSystem==0){
    //not centrality, but track multiplicity
    //track QA to extrack number of hybrid track for centrality binning
    fCenBin = GetCentralityBinPP(fNHybridTrack);
    AliInfo(Form("This is pp analysis. Nybrid = %d , fCenbin = %d , Zvtx = %f cm, fZvtx = %d.",fNHybridTrack,fCenBin,fVertex[2],fZvtx));
  }

  if(fIsMC) ProcessMC();

  if(fIsMC && fIsJJMC){
    //This should be called after ProcessMC(), because fMCArrayAOD is used.
    Int_t firstJetindex = fJJMCHandler->GetFirstJetIndex();
    Int_t lastJetindex  = fJJMCHandler->GetLastJetIndex();
    Int_t genIDJet      = fJJMCHandler->GetGeneratorJetIndex();

    const Int_t multClust = fPHOSClusterArray->GetEntriesFast();
    for(Int_t i=0;i<multClust;i++){
      AliCaloPhoton *ph = (AliCaloPhoton*)fPHOSClusterArray->At(i);
      Int_t primary = ph->GetPrimary();
      AliAODMCParticle *p = (AliAODMCParticle*)fMCArrayAOD->At(primary);

      if(p->GetGeneratorIndex() != genIDJet) fPHOSClusterArray->Remove(ph);
    }
    fPHOSClusterArray->Compress();
  }

  if(!fPHOSEvents[fZvtx][fCenBin]) fPHOSEvents[fZvtx][fCenBin] = new TList();
  TList *prevPHOS = fPHOSEvents[fZvtx][fCenBin];

  if(IsPHOSTriggerAnalysis){
    TriggerQA();
    SelectTriggeredCluster();
    EstimateTriggerEfficiency();
  }

  //cell QA
  CellQA();
  //cluster QA
  ClusterQA();

  FillPhoton();
  FillMgg();
  FillMixMgg();

  EstimatePIDCutEfficiency();
  EstimateTOFCutEfficiency();

  //EstimateNcellCutEfficiency();
  //EstimateM20CutEfficiency();
  //EstimateM02CutEfficiency();
  //EstimateSTDCutEfficiency();

  if(fIsNonLinStudyNeeded) DoNonLinearityStudy();

  const Int_t Nmixed = fCentNMixed.GetSize();

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

    if(prevPHOS->GetSize() > fCentNMixed[fCenBin]){//Remove redundant events
      TClonesArray * tmp = static_cast<TClonesArray*>(prevPHOS->Last());
      prevPHOS->RemoveLast();
      delete tmp;
      tmp = NULL;
    }
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

  Double_t dcaXY=0, dcaZ=0, pT=0, eta=0, phi=0, dEdx=0, p=0;
  Int_t NHybrid=0;
  Int_t NGlobal=0;
  Int_t NComplementary=0;
  Int_t PID=-1, charge=0;

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
      FillHistogramTH2(fOutputContainer,"hHybridTrackDCA",dcaXY,dcaZ);//in cm.

      if(fESDtrackCutsGlobal->AcceptTrack(esdtrack)){//global track
        NGlobal++;
        FillHistogramTH1(fOutputContainer,"hGlobalTrackPt",pT);
        FillHistogramTH2(fOutputContainer,"hGlobalTrackEtaPhi",phi,eta);
        FillHistogramTH2(fOutputContainer,"hGlobalTrackDCA",dcaXY,dcaZ);//in cm.
      }
      if(fESDtrackCutsGlobalConstrained->AcceptTrack(esdtrack)){//global-constrained = complementary track
        NComplementary++;
        FillHistogramTH1(fOutputContainer,"hComplementaryTrackPt",pT);
        FillHistogramTH2(fOutputContainer,"hComplementaryTrackEtaPhi",phi,eta);
        FillHistogramTH2(fOutputContainer,"hComplementaryTrackDCA",dcaXY,dcaZ);//in cm.
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
        dcaXY = aodtrack->DCA();
        dcaZ  = aodtrack->ZAtDCA();
        eta = aodtrack->Eta();
        phi = aodtrack->Phi();
        if(phi<0) phi += TMath::TwoPi();
        p = aodtrack->P();
        dEdx = aodtrack->GetTPCsignal();

        FillHistogramTH2(fOutputContainer,"hTrackTPCdEdx",p,dEdx);

        FillHistogramTH1(fOutputContainer,"hHybridTrackPt",pT);
        FillHistogramTH2(fOutputContainer,"hHybridTrackEtaPhi",phi,eta);
        FillHistogramTH2(fOutputContainer,"hHybridTrackDCA",dcaXY,dcaZ);//in cm.

        if(aodtrack->IsGlobalConstrained()){//constrained to primary vertex, instead of SPD hits.//complementary track
          NComplementary++;
          FillHistogramTH1(fOutputContainer,"hComplementaryTrackPt",pT);
          FillHistogramTH2(fOutputContainer,"hComplementaryTrackEtaPhi",phi,eta);
          FillHistogramTH2(fOutputContainer,"hComplementaryTrackDCA",dcaXY,dcaZ);//in cm.
        }
        else{//global track
          NGlobal++;
          FillHistogramTH1(fOutputContainer,"hGlobalTrackPt",pT);
          FillHistogramTH2(fOutputContainer,"hGlobalTrackEtaPhi",phi,eta);
          FillHistogramTH2(fOutputContainer,"hGlobalTrackDCA",dcaXY,dcaZ);//in cm.
        }

      }
    }//end of track loop


  }//end of AOD

  FillHistogramTH1(fOutputContainer,"hHybridTrackMult",NHybrid);
  FillHistogramTH1(fOutputContainer,"hGlobalTrackMult",NGlobal);
  FillHistogramTH1(fOutputContainer,"hComplementaryTrackMult",NComplementary);

  FillHistogramTH2(fOutputContainer,"hCentralityvsTrackMultiplicity",fCentralityMain,NHybrid);
  FillHistogramTH2(fOutputContainer,"hCentralityV0MvsTrackMultiplicity",fCentralityV0M, NHybrid);
  FillHistogramTH2(fOutputContainer,"hCentralityZNAvsTrackMultiplicity",fCentralityZNA, NHybrid);
  FillHistogramTH2(fOutputContainer,"hCentralityZNCvsTrackMultiplicity",fCentralityZNC, NHybrid);
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

  AliAODEvent *aod = dynamic_cast<AliAODEvent*>(fEvent);

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
  Double_t energy=0, tof=0;

  Bool_t IsMatched = kFALSE;
  Bool_t IsFilledOnce = kFALSE;//flag to avoid double counting

  AliPHOSTriggerHelper *helper = (AliPHOSTriggerHelper*)(fPHOSEventCuts->GetPHOSTriggerHelper());
  const Int_t L1input = helper->GetL1TriggerInput();
  const Int_t L0input = helper->GetL0TriggerInput();

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

    //cout << "tmod = "<< tmod << " , trgabsId = " << trgabsId << endl;

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

      Int_t maxAbsId = helper->FindHighestAmplitudeCellAbsId(clu1, cells);
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
      FillHistogramTH3(fOutputContainer,Form("hDistanceXZM%dTRU%d",trgmodule,helper->WhichTRU(trgcellx,trgcellz)),diffx,diffz,energy);

      if(!IsFilledOnce){//to avoid double counting of energy
        FillHistogramTH1(fOutputContainer,Form("hClusterEnergyM%d",module),energy);
        FillHistogramTH1(fOutputContainer,Form("hClusterEnergyM%dTRU%d",module,helper->WhichTRU(cellx,cellz)),energy);
        FillHistogramTH3(fOutputContainer,Form("hClusterHitM%d",module),cellx,cellz,energy);
      }

      if(helper->IsMatched(trgrelId,relId)){
//        ph->SetTrig(kTRUE);
        FillHistogramTH3(fOutputContainer,Form("hMatchedFiredTRUChannelM%d",trgmodule),trgcellx,trgcellz,energy);
        FillHistogramTH3(fOutputContainer,Form("hMatchedDistanceXZM%dTRU%d",trgmodule,helper->WhichTRU(trgcellx,trgcellz)),diffx,diffz,energy);

        if(kUsedCluster[i] == 0){
          FillHistogramTH1(fOutputContainer,Form("hMatchedClusterEnergyM%d",module),energy);
          FillHistogramTH1(fOutputContainer,Form("hMatchedClusterEnergyM%dTRU%d",trgmodule,helper->WhichTRU(trgcellx,trgcellz)),energy);
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

  Double_t cellamp=0;
  Int_t module=0,cellx=0,cellz=0;
  Int_t relId[4]={};
  Double_t position[3] = {};
  Int_t digMult=0;
  Double_t energy=0,tof=0,eta=0,phi=0;
  Double_t DistToBadChannel = 0;
  Double_t M20=0, M02=0;
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

    eta = global1.Eta();
    phi = global1.Phi();
    if(phi<0) phi += TMath::TwoPi();
    FillHistogramTH2(fOutputContainer,"hClusterEtaPhi",eta,phi);

    FillHistogramTH2(fOutputContainer,Form("hClusterEvsNM%d",module),energy,digMult);
    FillHistogramTH2(fOutputContainer,Form("hClusterEvsTM%d",module),energy,tof*1e+9);
    FillHistogramTH2(fOutputContainer,Form("hCluNXZM%d",module),cellx,cellz);
    FillHistogramTH2(fOutputContainer,Form("hCluEXZM%d",module),cellx,cellz,energy);
    FillHistogramTH3(fOutputContainer,"hEvsNvsM02",energy,digMult,M02);
   
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
      FillHistogramTH2(fOutputContainer,Form("hRvsTrackPtM%d",module),r,trackPt);
    } 
 
    FillHistogramTH2(fOutputContainer,Form("hFullDispvsFullEM%d",module),energy,R);
    FillHistogramTH2(fOutputContainer,Form("hCoreDispvsCoreEM%d",module),coreE,coreR);
    FillHistogramTH2(fOutputContainer,Form("hFullDispvsCoreEM%d",module),coreE,R);
    FillHistogramTH2(fOutputContainer,Form("hCoreDispvsFullEM%d",module),energy,coreR);

    DistToBadChannel = ph->DistToBad()/100.;
    FillHistogramTH2(fOutputContainer,"hEnergyvsDistanceToBadChannel",energy,DistToBadChannel);

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

  FillHistogramTH2(fOutputContainer,"hCentralityvsPHOSClusterMultiplicity"   ,fCentralityMain  ,multPHOSClust[0]);
  FillHistogramTH2(fOutputContainer,"hCentralityvsPHOSClusterMultiplicityTOF",fCentralityMain  ,multPHOSClustTOF[0]);

  FillHistogramTH2(fOutputContainer,"hCentralityV0MvsPHOSClusterMultiplicity"   ,fCentralityV0M,multPHOSClust[0]);
  FillHistogramTH2(fOutputContainer,"hCentralityV0MvsPHOSClusterMultiplicityTOF",fCentralityV0M,multPHOSClustTOF[0]);

  FillHistogramTH2(fOutputContainer,"hCentralityZNAvsPHOSClusterMultiplicity"   ,fCentralityZNA,multPHOSClust[0]);
  FillHistogramTH2(fOutputContainer,"hCentralityZNAvsPHOSClusterMultiplicityTOF",fCentralityZNA,multPHOSClustTOF[0]);

  FillHistogramTH2(fOutputContainer,"hCentralityZNCvsPHOSClusterMultiplicity"   ,fCentralityZNA,multPHOSClust[0]);
  FillHistogramTH2(fOutputContainer,"hCentralityZNCvsPHOSClusterMultiplicityTOF",fCentralityZNA,multPHOSClustTOF[0]);

  FillHistogramTH2(fOutputContainer,"hPHOSClusterMultiplicityvsTrackMultiplicity",fNHybridTrack,multPHOSClust[0]);
  FillHistogramTH2(fOutputContainer,"hPHOSClusterMultiplicityTOFvsTrackMultiplicity",fNHybridTrack,multPHOSClustTOF[0]);

}
//________________________________________________________________________
void AliAnalysisTaskPHOSPi0EtaToGammaGamma::FillPhoton() 
{
  
  AliESDEvent *esd = dynamic_cast<AliESDEvent*>(fEvent);
  AliAODEvent *aod = dynamic_cast<AliAODEvent*>(fEvent);

  const Int_t multClust = fPHOSClusterArray->GetEntriesFast();
  AliAODMCParticle *p = 0x0;

  Double_t pT=0,energy=0;
  Double_t eff=1;
  TF1 *f1tof = GetTOFCutEfficiencyFunction();

  Int_t firstJetindex = -1;
  Int_t lastJetindex  = -1;
  Int_t genIDJet = -1;

  if(fIsJJMC){
    firstJetindex = fJJMCHandler->GetFirstJetIndex();
    lastJetindex  = fJJMCHandler->GetLastJetIndex();
    genIDJet = fJJMCHandler->GetGeneratorJetIndex();
  }
  AliInfo(Form("genIDJet = %d , firstindexJet = %d , lastindexJet = %d.",genIDJet,firstJetindex,lastJetindex));

  Double_t weight = 1.;
  TF1 *f1Pi0Weight = 0x0;
  if(fIsMC) f1Pi0Weight = (TF1*)GetAdditionalPi0PtWeightFunction();
  Int_t primary = -1;
  Double_t TruePi0Pt = 1.;

  for(Int_t iph=0;iph<multClust;iph++){
    AliCaloPhoton *ph = (AliCaloPhoton*)fPHOSClusterArray->At(iph);
    if(!fPHOSClusterCuts->AcceptPhoton(ph)) continue;
    pT = ph->Pt();
    energy = ph->Energy();

    if(fUseCoreEnergy){
      pT = (ph->GetMomV2())->Pt();
      energy = (ph->GetMomV2())->Energy();
    }

    eff = f1tof->Eval(energy);

    if(fIsMC){
      primary = ph->GetPrimary();
      Bool_t isFromPi0 = IsFromPi0(primary,TruePi0Pt);
      if(isFromPi0){
        weight = f1Pi0Weight->Eval(TruePi0Pt);
        //cout << "true pi0 pT = " << TruePi0Pt << endl;
      }
      if(fIsJJMC){
        if(esd){
          if(primary < firstJetindex || lastJetindex < primary) continue;
        }
        else if(aod){
          //if(primary < firstJetindex || lastJetindex < primary) continue;
          p = (AliAODMCParticle*)fMCArrayAOD->At(primary);
          if(p->GetGeneratorIndex() != genIDJet) continue;
        }
      }
    }

    FillHistogramTH1(fOutputContainer,Form("hPhotonPt_Cen%d",fCenBin),pT,weight);
    //FillHistogramTH1(fOutputContainer,Form("hPhotonPt_M%d_Cen%d",ph->Module(),fCenBin),pT,weight);

    if(ph->IsTOFOK()){
      FillHistogramTH1(fOutputContainer,Form("hPhotonPt_Cen%d_TOF",fCenBin),pT,1/eff * weight);
      //FillHistogramTH1(fOutputContainer,Form("hPhotonPt_M%d_Cen%d_TOF",ph->Module(),fCenBin),pT,1/eff * weight);
    }

    if(fIsMC){
      if(IsPhoton(primary)) FillHistogramTH1(fOutputContainer,Form("hPurityGamma_Cen%d",fCenBin),pT,weight);
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

  Double_t eff1=1, eff2=1, eff12=1;
  TF1 *f1tof = GetTOFCutEfficiencyFunction();

  Double_t weight = 1., w1 = 1., w2 = 1.;
  TF1 *f1Pi0Weight = 0x0;
  TF1 *f1K0SWeight = 0x0;

  if(fIsMC){
    f1Pi0Weight = (TF1*)GetAdditionalPi0PtWeightFunction();
    f1K0SWeight = (TF1*)GetAdditionalK0SPtWeightFunction();
  }

  Int_t primary1 = -1;
  Int_t primary2 = -1;
  Double_t TruePi0Pt1 = 1.;
  Double_t TruePi0Pt2 = 1.;
  Int_t commonID = -1;
  AliAODMCParticle *p = 0x0;

  Int_t genIDJet      = -1;
  Int_t firstJetindex = -1;
  Int_t lastJetindex  = -1;

  if(fIsJJMC){
    genIDJet      = fJJMCHandler->GetGeneratorJetIndex();
    firstJetindex = fJJMCHandler->GetFirstJetIndex();
    lastJetindex  = fJJMCHandler->GetLastJetIndex();
    AliInfo(Form("firstindexJet = %d , lastindexJet = %d.",firstJetindex,lastJetindex));
  }
  Double_t TrueK0SPt = 0;
  Double_t TrueL0Pt = 0;

  for(Int_t i1=0;i1<multClust-1;i1++){
    AliCaloPhoton *ph1 = (AliCaloPhoton*)fPHOSClusterArray->At(i1);
    if(!fPHOSClusterCuts->AcceptPhoton(ph1)) continue;

    if(fIsMC){
      w1=1.;
      primary1 = ph1->GetPrimary();
      Bool_t isFromPi0 = IsFromPi0(primary1,TruePi0Pt1);
      if(isFromPi0){
        w1 = f1Pi0Weight->Eval(TruePi0Pt1);
        //cout << "true pi0 pT1 = " << TruePi0Pt1 << endl;
      }

      if(fIsJJMC){
        p = (AliAODMCParticle*)fMCArrayAOD->At(primary1);
        if(p->GetGeneratorIndex() != genIDJet) continue;
        //if(primary1 < firstJetindex || lastJetindex < primary1) continue;
      }
    }

    for(Int_t i2=i1+1;i2<multClust;i2++){
      AliCaloPhoton *ph2 = (AliCaloPhoton*)fPHOSClusterArray->At(i2);
      if(!fPHOSClusterCuts->AcceptPhoton(ph2)) continue;

      if(fIsMC){
        w2=1.;
        primary2 = ph2->GetPrimary();
        Bool_t isFromPi0 = IsFromPi0(primary2,TruePi0Pt2);
        if(isFromPi0){
          w2 = f1Pi0Weight->Eval(TruePi0Pt2);
          //cout << "true pi0 pT2 = " << TruePi0Pt2 << endl;
        }
        if(fIsJJMC){
          p = (AliAODMCParticle*)fMCArrayAOD->At(primary2);
          if(p->GetGeneratorIndex() != genIDJet) continue;
          //if(primary2 < firstJetindex || lastJetindex < primary2) continue;
        }
      }

      e1 = ph1->Energy();
      e2 = ph2->Energy();

      p12 = *ph1 + *ph2;
      m12 = p12.M();
      pt12 = p12.Pt();
      asym  = TMath::Abs((ph1->Energy()-ph2->Energy())/(ph1->Energy()+ph2->Energy()));//always full energy

      if(fUseCoreEnergy){
        p12core = *(ph1->GetMomV2()) + *(ph2->GetMomV2());
        m12 = p12core.M();
        pt12 = p12core.Pt();

        e1 = (ph1->GetMomV2())->Energy();
        e2 = (ph2->GetMomV2())->Energy();

        //cout << "core e1 = " << e1 << " , core e2 = " << e2 << endl;
      }

      eff1 = f1tof->Eval(e1);
      eff2 = f1tof->Eval(e2);
      eff12 = eff1 * eff2;

      if(fIsMC){
        commonID = FindCommonParent(primary1,primary2);
        if(commonID > -1) weight = w1;
        else weight = w1*w2;

        //for feed down correction from K0S->pi0 + pi0
        if(IsFromK0S(primary1,primary2,TrueK0SPt)){
          weight *= f1K0SWeight->Eval(TrueK0SPt);
          FillHistogramTH2(fOutputContainer,Form("hMggFromK0S_Cen%d",fCenBin),m12,pt12,weight);
          if(asym < 0.8) FillHistogramTH2(fOutputContainer,Form("hMggFromK0S_Cen%d_asym08",fCenBin),m12,pt12,weight);
          //          if(asym < 0.7) FillHistogramTH2(fOutputContainer,Form("hMggFromK0S_Cen%d_asym07",fCenBin),m12,pt12,weight);
        }
        if(IsFromLambda0(primary1,primary2,TrueL0Pt)){//this contibution is neglegible
          weight *= 1.;
          FillHistogramTH2(fOutputContainer,Form("hMggFromLambda0_Cen%d",fCenBin),m12,pt12,weight);
          if(asym < 0.8) FillHistogramTH2(fOutputContainer,Form("hMggFromLambda0_Cen%d_asym08",fCenBin),m12,pt12,weight);
          //          if(asym < 0.7) FillHistogramTH2(fOutputContainer,Form("hMggFromLambda0_Cen%d_asym07",fCenBin),m12,pt12,weight);
        }

      }//end of if fIsMC

      if(TMath::Abs(ph1->Module()-ph2->Module()) < 2) FillHistogramTH2(fOutputContainer,Form("hMgg_M%d%d",TMath::Min(ph1->Module(),ph2->Module()), TMath::Max(ph1->Module(),ph2->Module())),m12,pt12,weight);
      FillHistogramTH2(fOutputContainer,Form("hMgg_Cen%d",fCenBin),m12,pt12,weight);

      FillHistogramTH2(fOutputContainer,"hAsymvsMgg",asym,m12,weight);
      if(0.12 < m12 && m12 < 0.15) FillHistogramTH2(fOutputContainer,"hAsymvsPt",asym,pt12,weight);

      if(asym < 0.8) FillHistogramTH2(fOutputContainer,Form("hMgg_Cen%d_asym08",fCenBin),m12,pt12,weight);
//      if(asym < 0.7) FillHistogramTH2(fOutputContainer,Form("hMgg_Cen%d_asym07",fCenBin),m12,pt12,weight);

      if(ph1->IsTOFOK() && ph2->IsTOFOK()){
        FillHistogramTH2(fOutputContainer,Form("hMgg_Cen%d_TOF",fCenBin),m12,pt12,1/eff12 * weight);

        if(TMath::Abs(ph1->Module()-ph2->Module()) < 2) FillHistogramTH2(fOutputContainer,Form("hMgg_M%d%d_TOF",TMath::Min(ph1->Module(),ph2->Module()), TMath::Max(ph1->Module(),ph2->Module())),m12,pt12,1/eff12 * weight);

        if(asym < 0.8) FillHistogramTH2(fOutputContainer,Form("hMgg_Cen%d_TOF_asym08",fCenBin),m12,pt12,1/eff12 * weight);
//        if(asym < 0.7) FillHistogramTH2(fOutputContainer,Form("hMgg_Cen%d_TOF_asym07",fCenBin),m12,pt12,1/eff12 * weight);

      }//end of TOF cut

    }//end of ph2

  }//end of ph1

}
//________________________________________________________________________
void AliAnalysisTaskPHOSPi0EtaToGammaGamma::FillMixMgg() 
{
  TList *prevPHOS = fPHOSEvents[fZvtx][fCenBin];

  const Int_t multClust = fPHOSClusterArray->GetEntriesFast();
  AliAODMCParticle *p = 0x0;

  TLorentzVector p12, p12core;
  Double_t m12=0,pt12=0,asym=0,m12core=0,pt12core=0;
  Double_t eff1=1, eff2=1, eff12=1;
  Double_t e1=0,e2=0;
  TF1 *f1tof = GetTOFCutEfficiencyFunction();

  Int_t primary1 = -1;
  Int_t primary2 = -1;

  Int_t genIDJet      = -1;
  Int_t firstJetindex = -1;
  Int_t lastJetindex  = -1;

  if(fIsJJMC){
    genIDJet      = fJJMCHandler->GetGeneratorJetIndex();
    firstJetindex = fJJMCHandler->GetFirstJetIndex();
    lastJetindex  = fJJMCHandler->GetLastJetIndex();
    AliInfo(Form("firstindexJet = %d , lastindexJet = %d.",firstJetindex,lastJetindex));
  }

  for(Int_t i1=0;i1<multClust;i1++){
    AliCaloPhoton *ph1 = (AliCaloPhoton*)fPHOSClusterArray->At(i1);
    if(!fPHOSClusterCuts->AcceptPhoton(ph1)) continue;

    if(fIsJJMC){
      primary1 = ph1->GetPrimary();
      p = (AliAODMCParticle*)fMCArrayAOD->At(primary1);
      if(p->GetGeneratorIndex() != genIDJet) continue;
      //if(primary1 < firstJetindex || lastJetindex < primary1) continue;
    }

    for(Int_t ev=0;ev<prevPHOS->GetSize();ev++){
      TClonesArray *mixPHOS = static_cast<TClonesArray*>(prevPHOS->At(ev));

      for(Int_t i2=0;i2<mixPHOS->GetEntriesFast();i2++){
        AliCaloPhoton *ph2 = (AliCaloPhoton*)mixPHOS->At(i2);
        if(!fPHOSClusterCuts->AcceptPhoton(ph2)) continue;

        //if(fIsJJMC){
        //  primary2 = ph2->GetPrimary();
        //  p = (AliAODMCParticle*)fMCArrayAOD->At(primary2);
        //  if(p->GetGeneratorIndex() != genIDJet) continue;
        //  //if(primary2 < firstJetindex || lastJetindex < primary2) continue;
        //}

        e1 = ph1->Energy();
        e2 = ph2->Energy();

        p12 = *ph1 + *ph2;
        m12 = p12.M();
        pt12 = p12.Pt();
        asym  = TMath::Abs((ph1->Energy()-ph2->Energy())/(ph1->Energy()+ph2->Energy()));

        if(fUseCoreEnergy){
          p12core = *(ph1->GetMomV2()) + *(ph2->GetMomV2());
          m12 = p12core.M();
          pt12 = p12core.Pt();

          e1 = (ph1->GetMomV2())->Energy();
          e2 = (ph2->GetMomV2())->Energy();
        }

        eff1 = f1tof->Eval(e1);
        eff2 = f1tof->Eval(e2);
        eff12 = eff1 * eff2;


        FillHistogramTH2(fOutputContainer,Form("hMixMgg_Cen%d",fCenBin),m12,pt12);
        if(TMath::Abs(ph1->Module()-ph2->Module())<2)FillHistogramTH2(fOutputContainer,Form("hMixMgg_M%d%d",TMath::Min(ph1->Module(),ph2->Module()), TMath::Max(ph1->Module(),ph2->Module())),m12,pt12);

        if(asym < 0.8) FillHistogramTH2(fOutputContainer,Form("hMixMgg_Cen%d_asym08",fCenBin),m12,pt12);
//        if(asym < 0.7) FillHistogramTH2(fOutputContainer,Form("hMixMgg_Cen%d_asym07",fCenBin),m12,pt12);

        if(ph1->IsTOFOK() && ph2->IsTOFOK()){
          FillHistogramTH2(fOutputContainer,Form("hMixMgg_Cen%d_TOF",fCenBin),m12,pt12,1/eff12);
          if(TMath::Abs(ph1->Module()-ph2->Module())<2)FillHistogramTH2(fOutputContainer,Form("hMixMgg_M%d%d_TOF",TMath::Min(ph1->Module(),ph2->Module()), TMath::Max(ph1->Module(),ph2->Module())),m12,pt12,1/eff12);

          if(asym < 0.8) FillHistogramTH2(fOutputContainer,Form("hMixMgg_Cen%d_TOF_asym08",fCenBin),m12,pt12,1/eff12);
//          if(asym < 0.7) FillHistogramTH2(fOutputContainer,Form("hMixMgg_Cen%d_TOF_asym07",fCenBin),m12,pt12,1/eff12);

        }//end of TOF cut

      }//end of mix

    }//end of ph2

  }//end of ph1

}
//________________________________________________________________________
void AliAnalysisTaskPHOSPi0EtaToGammaGamma::EstimatePIDCutEfficiency()
{
  //tag and probe method is used.

  const Int_t multClust = fPHOSClusterArray->GetEntriesFast();

  TLorentzVector p12, p12core;
  Double_t m12=0,pt12=0;
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

      FillHistogramTH2(fOutputContainer,Form("hMgg_Probe_PID_Cen%d",fCenBin),m12,energy);
      if(fPHOSClusterCuts->AcceptPhoton(ph2))
        FillHistogramTH2(fOutputContainer,Form("hMgg_PassingProbe_PID_Cen%d",fCenBin),m12,energy);

    }//end of ph2

  }//end of ph1

  //next mixed event
  TList *prevPHOS = fPHOSEvents[fZvtx][fCenBin];

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

        FillHistogramTH2(fOutputContainer,Form("hMixMgg_Probe_PID_Cen%d",fCenBin),m12,energy);
        if(fPHOSClusterCuts->AcceptPhoton(ph2))
          FillHistogramTH2(fOutputContainer,Form("hMixMgg_PassingProbe_PID_Cen%d",fCenBin),m12,energy);

      }//end of mix

    }//end of ph2

  }//end of ph1

}
//________________________________________________________________________
void AliAnalysisTaskPHOSPi0EtaToGammaGamma::EstimateNcellCutEfficiency()
{
  //tag and probe method is used.

  const Int_t multClust = fPHOSClusterArray->GetEntriesFast();

  TLorentzVector p12, p12core;
  Double_t m12=0,pt12=0;
  Double_t energy=0;

  for(Int_t i1=0;i1<multClust;i1++){
    AliCaloPhoton *ph1 = (AliCaloPhoton*)fPHOSClusterArray->At(i1);
    if(!fPHOSClusterCuts->AcceptPhoton(ph1)) continue;
    if(ph1->GetNCells() < 1) continue;

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

      FillHistogramTH2(fOutputContainer,Form("hMgg_Probe_Ncell_Cen%d",fCenBin),m12,energy);
      if(fPHOSClusterCuts->AcceptPhoton(ph2)){
        if(ph2->GetNCells() > 1)
          FillHistogramTH2(fOutputContainer,Form("hMgg_PassingProbe_Ncell_Cen%d",fCenBin),m12,energy);
      }
    }//end of ph2

  }//end of ph1

  //next mixed event
  TList *prevPHOS = fPHOSEvents[fZvtx][fCenBin];

  for(Int_t i1=0;i1<multClust;i1++){
    AliCaloPhoton *ph1 = (AliCaloPhoton*)fPHOSClusterArray->At(i1);
    if(!fPHOSClusterCuts->AcceptPhoton(ph1)) continue;
    if(ph1->GetNCells() < 2) continue;

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

        FillHistogramTH2(fOutputContainer,Form("hMixMgg_Probe_Ncell_Cen%d",fCenBin),m12,energy);
        if(fPHOSClusterCuts->AcceptPhoton(ph2)){
          if(ph2->GetNCells() > 1)
            FillHistogramTH2(fOutputContainer,Form("hMixMgg_PassingProbe_Ncell_Cen%d",fCenBin),m12,energy);
        }
      }//end of mix

    }//end of ph2

  }//end of ph1

}
//________________________________________________________________________
void AliAnalysisTaskPHOSPi0EtaToGammaGamma::EstimateM20CutEfficiency()
{
  //tag and probe method is used.

  const Int_t multClust = fPHOSClusterArray->GetEntriesFast();

  TLorentzVector p12, p12core;
  Double_t m12=0,pt12=0;
  Double_t energy=0;

  for(Int_t i1=0;i1<multClust;i1++){
    AliCaloPhoton *ph1 = (AliCaloPhoton*)fPHOSClusterArray->At(i1);
    if(!fPHOSClusterCuts->AcceptPhoton(ph1)) continue;
    if(ph1->GetLambda1() < 0.2) continue;

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

      FillHistogramTH2(fOutputContainer,Form("hMgg_Probe_M20_Cen%d",fCenBin),m12,energy);
      if(fPHOSClusterCuts->AcceptPhoton(ph2)){
        if(ph2->GetLambda1() > 0.2)
          FillHistogramTH2(fOutputContainer,Form("hMgg_PassingProbe_M20_Cen%d",fCenBin),m12,energy);
      }
    }//end of ph2

  }//end of ph1

  //next mixed event
  TList *prevPHOS = fPHOSEvents[fZvtx][fCenBin];

  for(Int_t i1=0;i1<multClust;i1++){
    AliCaloPhoton *ph1 = (AliCaloPhoton*)fPHOSClusterArray->At(i1);
    if(!fPHOSClusterCuts->AcceptPhoton(ph1)) continue;
    if(ph1->GetLambda1() < 0.2) continue;

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

        FillHistogramTH2(fOutputContainer,Form("hMixMgg_Probe_M20_Cen%d",fCenBin),m12,energy);
        if(fPHOSClusterCuts->AcceptPhoton(ph2)){
          if(ph2->GetLambda1() > 0.2)
            FillHistogramTH2(fOutputContainer,Form("hMixMgg_PassingProbe_M20_Cen%d",fCenBin),m12,energy);
        }
      }//end of mix

    }//end of ph2

  }//end of ph1

}
//________________________________________________________________________
void AliAnalysisTaskPHOSPi0EtaToGammaGamma::EstimateM02CutEfficiency()
{
  //tag and probe method is used.

  const Int_t multClust = fPHOSClusterArray->GetEntriesFast();

  TLorentzVector p12, p12core;
  Double_t m12=0,pt12=0;
  Double_t energy=0;

  for(Int_t i1=0;i1<multClust;i1++){
    AliCaloPhoton *ph1 = (AliCaloPhoton*)fPHOSClusterArray->At(i1);
    if(!fPHOSClusterCuts->AcceptPhoton(ph1)) continue;
    if(ph1->GetLambda2() < 0.2) continue;

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

      FillHistogramTH2(fOutputContainer,Form("hMgg_Probe_M02_Cen%d",fCenBin),m12,energy);
      if(fPHOSClusterCuts->AcceptPhoton(ph2)){
        if(ph2->GetLambda2() > 0.2)
          FillHistogramTH2(fOutputContainer,Form("hMgg_PassingProbe_M02_Cen%d",fCenBin),m12,energy);
      }
    }//end of ph2

  }//end of ph1

  //next mixed event
  TList *prevPHOS = fPHOSEvents[fZvtx][fCenBin];

  for(Int_t i1=0;i1<multClust;i1++){
    AliCaloPhoton *ph1 = (AliCaloPhoton*)fPHOSClusterArray->At(i1);
    if(!fPHOSClusterCuts->AcceptPhoton(ph1)) continue;
    if(ph1->GetLambda2() < 0.2) continue;

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

        FillHistogramTH2(fOutputContainer,Form("hMixMgg_Probe_M02_Cen%d",fCenBin),m12,energy);
        if(fPHOSClusterCuts->AcceptPhoton(ph2)){
          if(ph2->GetLambda2() > 0.2)
            FillHistogramTH2(fOutputContainer,Form("hMixMgg_PassingProbe_M02_Cen%d",fCenBin),m12,energy);
        }
      }//end of mix

    }//end of ph2

  }//end of ph1

}
//________________________________________________________________________
void AliAnalysisTaskPHOSPi0EtaToGammaGamma::EstimateSTDCutEfficiency()
{
  //tag and probe method is used.

  const Int_t multClust = fPHOSClusterArray->GetEntriesFast();

  TLorentzVector p12, p12core;
  Double_t m12=0,pt12=0;
  Double_t energy=0;

  for(Int_t i1=0;i1<multClust;i1++){
    AliCaloPhoton *ph1 = (AliCaloPhoton*)fPHOSClusterArray->At(i1);
    if(!fPHOSClusterCuts->AcceptPhoton(ph1)) continue;
    if( ph1->Energy() < 0.3 // MIP cut
        || ph1->GetNCells() < 3 //accidental noise
        || ph1->GetLambda2() < 0.2 //too small shower shape
      ) continue;

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

      FillHistogramTH2(fOutputContainer,Form("hMgg_Probe_STDCut_Cen%d",fCenBin),m12,energy);
      if(fPHOSClusterCuts->AcceptPhoton(ph2)){
        if(ph2->Energy() > 0.3 // MIP cut
        && ph2->GetNCells() > 2 //accidental noise
        && ph2->GetLambda2() > 0.2 //too small shower shape. Lambda1 = M02 (short axis), Lambda2 = M02 (long axis)
          )
          FillHistogramTH2(fOutputContainer,Form("hMgg_PassingProbe_STDCut_Cen%d",fCenBin),m12,energy);
      }
    }//end of ph2

  }//end of ph1

  //next mixed event
  TList *prevPHOS = fPHOSEvents[fZvtx][fCenBin];

  for(Int_t i1=0;i1<multClust;i1++){
    AliCaloPhoton *ph1 = (AliCaloPhoton*)fPHOSClusterArray->At(i1);
    if(!fPHOSClusterCuts->AcceptPhoton(ph1)) continue;
    if( ph1->Energy() < 0.3 // MIP cut
        || ph1->GetNCells() < 3 //accidental noise
        || ph1->GetLambda2() < 0.2 //too small shower shape
      ) continue;

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

        FillHistogramTH2(fOutputContainer,Form("hMixMgg_Probe_STDCut_Cen%d",fCenBin),m12,energy);
        if(fPHOSClusterCuts->AcceptPhoton(ph2)){
          if(ph2->Energy() > 0.3 // MIP cut
          && ph2->GetNCells() > 2 //accidental noise
          && ph2->GetLambda2() > 0.2 //too small shower shape
            )
            FillHistogramTH2(fOutputContainer,Form("hMixMgg_PassingProbe_STDCut_Cen%d",fCenBin),m12,energy);
        }
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
  Double_t m12=0,pt12=0;
  Double_t energy=0;

  for(Int_t i1=0;i1<multClust;i1++){
    AliCaloPhoton *ph1 = (AliCaloPhoton*)fPHOSClusterArray->At(i1);
    if(!fPHOSClusterCuts->AcceptPhoton(ph1)) continue;
    if(!ph1->IsTOFOK()) continue;

    for(Int_t i2=0;i2<multClust;i2++){
      AliCaloPhoton *ph2 = (AliCaloPhoton*)fPHOSClusterArray->At(i2);
      if(!fPHOSClusterCuts->AcceptPhoton(ph2)) continue;

      if(i2==i1) continue;//reject same cluster combination

      p12 = *ph1 + *ph2;
      m12 = p12.M();
      energy = ph2->Energy();

      if(fUseCoreEnergy){
        p12core = *(ph1->GetMomV2()) + *(ph2->GetMomV2());
        m12 = p12core.M();
        energy = (ph2->GetMomV2())->Energy();
      }

      FillHistogramTH2(fOutputContainer,Form("hMgg_Probe_TOF_Cen%d",fCenBin),m12,energy);
      if(ph2->IsTOFOK()) FillHistogramTH2(fOutputContainer,Form("hMgg_PassingProbe_TOF_Cen%d",fCenBin),m12,energy);

    }//end of ph2

  }//end of ph1

  //next mixed event
  TList *prevPHOS = fPHOSEvents[fZvtx][fCenBin];

  for(Int_t i1=0;i1<multClust;i1++){
    AliCaloPhoton *ph1 = (AliCaloPhoton*)fPHOSClusterArray->At(i1);
    if(!fPHOSClusterCuts->AcceptPhoton(ph1)) continue;
    if(!ph1->IsTOFOK()) continue;

    for(Int_t ev=0;ev<prevPHOS->GetSize();ev++){
      TClonesArray *mixPHOS = static_cast<TClonesArray*>(prevPHOS->At(ev));

      for(Int_t i2=0;i2<mixPHOS->GetEntriesFast();i2++){
        AliCaloPhoton *ph2 = (AliCaloPhoton*)mixPHOS->At(i2);
        if(!fPHOSClusterCuts->AcceptPhoton(ph2)) continue;

        p12 = *ph1 + *ph2;
        m12 = p12.M();
        energy = ph2->Energy();

        if(fUseCoreEnergy){
          p12core = *(ph1->GetMomV2()) + *(ph2->GetMomV2());
          m12 = p12core.M();
          energy = (ph2->GetMomV2())->Energy();
        }

        FillHistogramTH2(fOutputContainer,Form("hMixMgg_Probe_TOF_Cen%d",fCenBin),m12,energy);
        if(ph2->IsTOFOK()) FillHistogramTH2(fOutputContainer,Form("hMixMgg_PassingProbe_TOF_Cen%d",fCenBin),m12,energy);

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

  AliVCaloTrigger* trg = fEvent->GetCaloTrigger("PHOS");
  trg->Reset();
  
  AliVCaloCells *cells = dynamic_cast<AliVCaloCells*>(fEvent->GetPHOSCells());
  AliPHOSTriggerHelper *helper = (AliPHOSTriggerHelper*)(fPHOSEventCuts->GetPHOSTriggerHelper());
  const Int_t L1input = helper->GetL1TriggerInput();
  const Int_t L0input = helper->GetL0TriggerInput();

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

  for(Int_t i=0;i<multClust;i++){
    AliCaloPhoton *ph = (AliCaloPhoton*)fPHOSClusterArray->At(i);
    AliVCluster *clu1 = (AliVCluster*)ph->GetCluster();

    Int_t maxAbsId = helper->FindHighestAmplitudeCellAbsId(clu1, cells);
    //cout << "maxAbsId = " << maxAbsId << endl;

    fPHOSGeo->AbsToRelNumbering(maxAbsId,relId);
    module = relId[0];
    cellx  = relId[2];
    cellz  = relId[3];

    if(module < 1 || module > 4){
      AliError(Form("Wrong module number %d",module));
      return;
    }
    AliInfo(Form("cluster position M:%d, X:%d , Z:%d",module,cellx,cellz));

    while(trg->Next()){

      if(trg->GetL1TimeSum() != L1) continue;

      trg->GetPosition(tmod,trgabsId);//tmod is online module numbering. i.e., tmod=1 means half module.

      //cout << "tmod = "<< tmod << " , trgabsId = " << trgabsId << endl;

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
        ph->SetTrig(kTRUE);
        break;
      }
    }//end of while trigger patch loop

  }//end of cluster loop

  for(Int_t i=0;i<multClust;i++){
    AliCaloPhoton *ph = (AliCaloPhoton*)fPHOSClusterArray->At(i);
    if(!fPHOSClusterCuts->AcceptPhoton(ph)) continue;

    energy = ph->Energy();
    if(fUseCoreEnergy) energy = (ph->GetMomV2())->Energy();

    tru = helper->WhichTRU(cellx,cellz);
    FillHistogramTH1(fOutputContainer,Form("hClusterE_Cen%d_M%d_TRU%d",fCenBin,module,tru),energy);
    if(ph->IsTrig()) FillHistogramTH1(fOutputContainer,Form("hTriggeredClusterE_Cen%d_M%d_TRU%d",fCenBin,module,tru),energy);
  }//end of cluster loop


}
//________________________________________________________________________
void AliAnalysisTaskPHOSPi0EtaToGammaGamma::EstimateTriggerEfficiency()
{
  //tag and probe method is used.

  const Int_t multClust = fPHOSClusterArray->GetEntriesFast();

  Double_t position[3] = {};
  Double_t energy=0;
  Int_t module_tag=0,cellx_tag=0,cellz_tag=0,tru_tag=0;
  Int_t module=0,cellx=0,cellz=0,tru=0;
  Int_t relId[4]={};

  Int_t truch_tag = -1, chX_tag=-1, chZ_tag=-1;
  Int_t truch = -1, chX=-1, chZ=-1;

  TLorentzVector p12, p12core;
  Double_t m12=0,pt12=0;

  AliPHOSTriggerHelper *helper = (AliPHOSTriggerHelper*)(fPHOSEventCuts->GetPHOSTriggerHelper());

  for(Int_t i1=0;i1<multClust;i1++){
    AliCaloPhoton *ph1 = (AliCaloPhoton*)fPHOSClusterArray->At(i1);
    if(!fPHOSClusterCuts->AcceptPhoton(ph1)) continue;
    if(!ph1->IsTrig()) continue;

    position[0] = ph1->EMCx();
    position[1] = ph1->EMCy();
    position[2] = ph1->EMCz();
    relId[0] = 0; relId[1] = 0; relId[2] = 0; relId[3] = 0;
    TVector3 global_tag(position);
    fPHOSGeo->GlobalPos2RelId(global_tag,relId);

    module_tag = relId[0];
    cellx_tag  = relId[2];
    cellz_tag  = relId[3];
    tru_tag = helper->WhichTRU(cellx_tag,cellz_tag);
    truch_tag = helper->WhichTRUChannel(cellx_tag,cellz_tag,chX_tag,chZ_tag);
    //printf("cellx_tag = %d , cellz_tag = %d , chX_tag = %d , chZ_tag = %d , truch_tag = %d.\n",cellx_tag,cellz_tag,chX_tag,chZ_tag,truch_tag);

    for(Int_t i2=0;i2<multClust;i2++){
      AliCaloPhoton *ph2 = (AliCaloPhoton*)fPHOSClusterArray->At(i2);
      if(!fPHOSClusterCuts->AcceptPhoton(ph2)) continue;

      if(i2==i1) continue;//reject same cluster combination

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
      tru = helper->WhichTRU(cellx,cellz);
      truch = helper->WhichTRUChannel(cellx,cellz,chX,chZ);

      Int_t dm = module_tag - module;
      Int_t dt = tru_tag - tru;
      Int_t dx = chX_tag - chX;
      Int_t dz = chZ_tag - chZ;

      //if( (dm == 0)
      //    && (dt == 0)
      //    && (TMath::Abs(dx) < 2)
      //    && (TMath::Abs(dz) < 2)
      //  ) continue;//reject cluster pair where they belong to same 4x4 region.

      FillHistogramTH2(fOutputContainer,Form("hMgg_Probe_Trg_Cen%d",fCenBin),m12,energy);
      FillHistogramTH2(fOutputContainer,Form("hMgg_Probe_Trg_Cen%d_M%d_TRU%d",fCenBin,module,tru),m12,energy);

      if(ph2->IsTrig()){
        FillHistogramTH2(fOutputContainer,Form("hMgg_PassingProbe_Trg_Cen%d",fCenBin),m12,energy);
        FillHistogramTH2(fOutputContainer,Form("hMgg_PassingProbe_Trg_Cen%d_M%d_TRU%d",fCenBin,module,tru),m12,energy);
      }

    }//end of ph2

  }//end of ph1

  //next mixed event
  TList *prevPHOS = fPHOSEvents[fZvtx][fCenBin];

  for(Int_t i1=0;i1<multClust;i1++){
    AliCaloPhoton *ph1 = (AliCaloPhoton*)fPHOSClusterArray->At(i1);
    if(!fPHOSClusterCuts->AcceptPhoton(ph1)) continue;
    if(!ph1->IsTrig()) continue;

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

        position[0] = ph2->EMCx();
        position[1] = ph2->EMCy();
        position[2] = ph2->EMCz();

        relId[0] = 0; relId[1] = 0; relId[2] = 0; relId[3] = 0;
        TVector3 global1(position);
        fPHOSGeo->GlobalPos2RelId(global1,relId);

        module = relId[0];
        cellx  = relId[2];
        cellz  = relId[3];
        tru = helper->WhichTRU(cellx,cellz);

        FillHistogramTH2(fOutputContainer,Form("hMixMgg_Probe_Trg_Cen%d",fCenBin),m12,energy);
        FillHistogramTH2(fOutputContainer,Form("hMixMgg_Probe_Trg_Cen%d_M%d_TRU%d",fCenBin,module,tru),m12,energy);
        if(fPHOSClusterCuts->AcceptPhoton(ph2)){
          if(ph2->IsTrig()){
            FillHistogramTH2(fOutputContainer,Form("hMixMgg_PassingProbe_Trg_Cen%d",fCenBin),m12,energy);
            FillHistogramTH2(fOutputContainer,Form("hMixMgg_PassingProbe_Trg_Cen%d_M%d_TRU%d",fCenBin,module,tru),m12,energy);
          }
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
    AliInfo(Form("genIDJet = %d , firstindexJet = %d , lastindexJet = %d.",genIDJet,firstJetindex,lastJetindex));
  }

  TF1 *f1Pi0Weight = (TF1*)GetAdditionalPi0PtWeightFunction();
  TF1 *f1K0SWeight = (TF1*)GetAdditionalK0SPtWeightFunction();

  Int_t genID = -1;
  Double_t pT=0, rapidity=0, phi=0;
  Double_t weight = 1;
  Int_t pdg = 0;
  TString parname = "";
  TString genname = "";
  Int_t motherid = -1;
  AliAODMCParticle *mp = 0x0;//mother particle
  Double_t motherpT = 0;

  if(fESDEvent){//for ESD
    fMCArrayESD = (AliStack*)GetMCInfoESD();
    if(!fMCArrayESD){
      AliError("Could not get MC Stack!");
      return;
    }

    const Int_t Ntrack = fMCArrayESD->GetNtrack();
    for(Int_t i=0;i<Ntrack;i++){
      TParticle *p = (TParticle*)fMCArrayESD->Particle(i);


    }//end of generated particle loop

  }
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

      //if(fIsJJMC && (i < firstJetindex || lastJetindex < i) ) continue;
      if(fIsJJMC && genID != genIDJet) continue;

      pT = p->Pt();
      rapidity = p->Y();
      phi = p->Phi();
      pdg = p->PdgCode();

      //rapidity is Y(), but, pseudo-rapidity is Eta();

      if(pT < 1e-3) continue;//reject below 1 MeV
      if(TMath::Abs(rapidity) > 0.5) continue;

      motherid = p->GetMother();
      if(motherid > -1 && pdg == 111){//pi0 has mother
        mp = (AliAODMCParticle*)fMCArrayAOD->At(motherid);
        if(mp->PdgCode() == 310)//mother is K0S
          FillHistogramTH2(fOutputContainer,Form("hPi0FromK0S_Cen%d",fCenBin),pT,Rho(p),1.);
      }

      if(Rho(p) > 1.0) continue;
      weight = 1.;

      if(pdg==111){//pi0
        parname = "Pi0";
        weight = f1Pi0Weight->Eval(pT);
        //printf("This is pi0. pT = %e GeV/c, weight = %e.\n",pT, weight);

      }
      else if(pdg==221){//eta
        parname = "Eta";
      }
      else if(pdg==22){//gamma
        parname = "Gamma";

        motherid = p->GetMother();
        while(motherid > -1){//this particle has mother.
          mp = (AliAODMCParticle*)fMCArrayAOD->At(motherid);
          if(mp->PdgCode() == 111){ //mother is pi0.
            motherpT = mp->Pt();
            weight = f1Pi0Weight->Eval(motherpT);
            //printf("The mother of this gamma is pi0. motherpT = %e GeV/c, gamma pT = %e , weight = %e.\n",motherpT, pT, weight);
            FillHistogramTH1(fOutputContainer,Form("hGenGammaFromPi0Pt_Cen%d"    ,fCenBin),pT          ,weight);
            FillHistogramTH2(fOutputContainer,Form("hGenGammaFromPi0EtaPhi_Cen%d",fCenBin),rapidity,phi,weight);
            FillHistogramTH2(fOutputContainer,Form("hGenGammaFromPi0EtaPt_Cen%d" ,fCenBin),rapidity,pT ,weight);
            break;
          }
          motherid = mp->GetMother();
        }//end of while

      }
      else if(pdg==211 || pdg==-211){//pi+ or pi-
        parname = "ChargedPion";
      }
      else if(pdg==321 || pdg==-321){//K+ or K-
        parname = "ChargedKaon";
        weight = f1K0SWeight->Eval(pT);
        //printf("This is K+/-. pT = %e GeV/c, weight = %e.\n",pT, weight);
      }
      else if(pdg==310){//K0S
        parname = "K0S";
        weight = f1K0SWeight->Eval(pT);
        //printf("This is K0S. pT = %e GeV/c, weight = %e.\n",pT, weight);
      }
      else if(pdg==130){//K0L
        parname = "K0L";
        weight = f1K0SWeight->Eval(pT);
        //printf("This is K0L. pT = %e GeV/c, weight = %e.\n",pT, weight);
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

      FillHistogramTH1(fOutputContainer,Form("hGen%sPt_Cen%d"    ,parname.Data(),fCenBin),pT          ,weight);
      FillHistogramTH2(fOutputContainer,Form("hGen%sEtaPhi_Cen%d",parname.Data(),fCenBin),rapidity,phi,weight);
      FillHistogramTH2(fOutputContainer,Form("hGen%sEtaPt_Cen%d" ,parname.Data(),fCenBin),rapidity,pT ,weight);

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
Bool_t AliAnalysisTaskPHOSPi0EtaToGammaGamma::IsFromLambda0(Int_t label1, Int_t label2, Double_t &TrueL0Pt)
{
  AliAODMCParticle *p1 = (AliAODMCParticle*)fMCArrayAOD->At(label1);
  AliAODMCParticle *p2 = (AliAODMCParticle*)fMCArrayAOD->At(label2);

  AliAODMCParticle *mp = 0x0;
  Int_t motherid = -1;
  Int_t pdg=0;
  Double_t pT = 0;

  Bool_t IsK0S1 = kFALSE;
  Bool_t IsK0S2 = kFALSE;
  Int_t motherID1 = -1;
  Int_t motherID2 = -1;

  motherid = p1->GetMother();
  while(motherid > -1){

    mp = (AliAODMCParticle*)fMCArrayAOD->At(motherid);
    pdg = mp->PdgCode(); 
    pT = mp->Pt();

    if(pdg==3122){//Lambda0
      IsK0S1 = kTRUE;
      TrueL0Pt = pT;
      motherID1 = motherid;
      break;
    }

    motherid = mp->GetMother();

  }

  motherid = p2->GetMother();
  while(motherid > -1){
    mp = (AliAODMCParticle*)fMCArrayAOD->At(motherid);

    pdg = mp->PdgCode();
    if(pdg==3122){//Lambda0
      IsK0S2 = kTRUE;
      TrueL0Pt = pT;
      motherID2 = motherid;
      break;
    }

    motherid = mp->GetMother();

  }

  return (motherID1 == motherID2) & (IsK0S1 & IsK0S2);

}
//________________________________________________________________________
Bool_t AliAnalysisTaskPHOSPi0EtaToGammaGamma::IsFromK0S(Int_t label1, Int_t label2, Double_t &TrueK0SPt)
{
  AliAODMCParticle *p1 = (AliAODMCParticle*)fMCArrayAOD->At(label1);
  AliAODMCParticle *p2 = (AliAODMCParticle*)fMCArrayAOD->At(label2);

  AliAODMCParticle *mp = 0x0;
  Int_t motherid = -1;
  Int_t pdg=0;
  Double_t pT = 0;

  Bool_t IsK0S1 = kFALSE;
  Bool_t IsK0S2 = kFALSE;

  Int_t motherID1 = -1;
  Int_t motherID2 = -1;

  motherid = p1->GetMother();
  while(motherid > -1){

    mp = (AliAODMCParticle*)fMCArrayAOD->At(motherid);
    pdg = mp->PdgCode(); 
    pT = mp->Pt();

    if(pdg==310){//K0S
      IsK0S1 = kTRUE;
      TrueK0SPt = pT;
      motherID1 = motherid;
      break;
    }

    motherid = mp->GetMother();

  }

  motherid = p2->GetMother();
  while(motherid > -1){
    mp = (AliAODMCParticle*)fMCArrayAOD->At(motherid);
    pdg = mp->PdgCode();
    pT = mp->Pt();

    if(pdg==310){//K0S
      IsK0S2 = kTRUE;
      TrueK0SPt = pT;
      motherID2 = motherid;
      break;
    }

    motherid = mp->GetMother();

  }

  return (motherID1 == motherID2) & (IsK0S1 & IsK0S2);

}
//________________________________________________________________________
Bool_t AliAnalysisTaskPHOSPi0EtaToGammaGamma::IsFromPi0(Int_t label1, Double_t &TruePi0Pt)
{

  AliAODMCParticle *p1 = (AliAODMCParticle*)fMCArrayAOD->At(label1);
 
  AliAODMCParticle *mp = 0x0;
  Int_t motherid = -1;
  Int_t pdg=0;
  Double_t pT = 0;

  Bool_t IsPi0 = kFALSE;
  motherid = p1->GetMother();

  while(motherid > -1){
    mp = (AliAODMCParticle*)fMCArrayAOD->At(motherid);
    pT = mp->Pt();
    pdg = mp->PdgCode(); 

    if(pdg==111 && Rho(mp) < 1.0){//pi0 from primary vertex
      IsPi0 = kTRUE;
      TruePi0Pt = pT;
      break;
    }

    motherid = mp->GetMother();

  }

  return IsPi0;

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
void AliAnalysisTaskPHOSPi0EtaToGammaGamma::DoNonLinearityStudy()
{

  const Int_t multClust = fPHOSClusterArray->GetEntriesFast();
  TLorentzVector p12, p12core;

  Double_t m12=0,pt12=0,asym=0;
  Double_t e1=0,e2=0;
  for(Int_t ia=0;ia<7;ia++){
    for(Int_t ib=0;ib<7;ib++){

      for(Int_t i1=0;i1<multClust-1;i1++){
        AliCaloPhoton *ph1 = (AliCaloPhoton*)fPHOSClusterArray->At(i1);
        if(!fPHOSClusterCuts->AcceptPhoton(ph1)) continue;

        for(Int_t i2=i1+1;i2<multClust;i2++){
          AliCaloPhoton *ph2 = (AliCaloPhoton*)fPHOSClusterArray->At(i2);
          if(!fPHOSClusterCuts->AcceptPhoton(ph2)) continue;

          e1 = ph1->Energy();
          e2 = ph2->Energy();

          p12 = *ph1*fNonLin[fCenBin][ia][ib]->Eval(e1) + *ph2*fNonLin[fCenBin][ia][ib]->Eval(e2);
          m12 = p12.M();
          pt12 = p12.Pt();

          if(fUseCoreEnergy){
            e1 = (ph1->GetMomV2())->Energy();
            e2 = (ph2->GetMomV2())->Energy();

            p12core = *(ph1->GetMomV2())*fNonLin[fCenBin][ia][ib]->Eval(e1) + *(ph2->GetMomV2())*fNonLin[fCenBin][ia][ib]->Eval(e2);
            m12 = p12core.M();
            pt12 = p12core.Pt();
          }

          FillHistogramTH2(fOutputContainer,Form("hMgg_Cen%d_a%d_b%d",fCenBin,ia,ib),m12,pt12);

        }//end of ph2
      }//end of ph1
    }
  }

  TList *prevPHOS = fPHOSEvents[fZvtx][fCenBin];

  for(Int_t ia=0;ia<7;ia++){
    for(Int_t ib=0;ib<7;ib++){

      for(Int_t i1=0;i1<multClust;i1++){
        AliCaloPhoton *ph1 = (AliCaloPhoton*)fPHOSClusterArray->At(i1);
        if(!fPHOSClusterCuts->AcceptPhoton(ph1)) continue;

        for(Int_t ev=0;ev<prevPHOS->GetSize();ev++){
          TClonesArray *mixPHOS = static_cast<TClonesArray*>(prevPHOS->At(ev));

          for(Int_t i2=0;i2<mixPHOS->GetEntriesFast();i2++){
            AliCaloPhoton *ph2 = (AliCaloPhoton*)mixPHOS->At(i2);
            if(!fPHOSClusterCuts->AcceptPhoton(ph2)) continue;

            e1 = ph1->Energy();
            e2 = ph2->Energy();

            p12 = *ph1*fNonLin[fCenBin][ia][ib]->Eval(e1) + *ph2*fNonLin[fCenBin][ia][ib]->Eval(e2);
            m12 = p12.M();
            pt12 = p12.Pt();

            if(fUseCoreEnergy){
              e1 = (ph1->GetMomV2())->Energy();
              e2 = (ph2->GetMomV2())->Energy();
              p12core = *(ph1->GetMomV2())*fNonLin[fCenBin][ia][ib]->Eval(e1) + *(ph2->GetMomV2())*fNonLin[fCenBin][ia][ib]->Eval(e2);
              m12 = p12core.M();
              pt12 = p12core.Pt();
            }

            FillHistogramTH2(fOutputContainer,Form("hMixMgg_Cen%d_a%d_b%d",fCenBin,ia,ib),m12,pt12);

          }//end of mix

        }//end of ph2

      }//end of ph1

    }

  }

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
AliPHOSGeometry *AliAnalysisTaskPHOSPi0EtaToGammaGamma::GetPHOSGeometry()
{
  AliAODEvent *aod = dynamic_cast<AliAODEvent*>(fEvent);
  AliESDEvent *esd = dynamic_cast<AliESDEvent*>(fEvent);
  Int_t RunNumber = fEvent->GetRunNumber();

  //Initialize the PHOS geometry

  fPHOSGeo = 0x0;

  if(fUsePHOSTender){
    AliInfo("PHOSTender is used.");

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
