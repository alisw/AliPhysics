#include "AliAnalysisTaskPi0V2.h"

#include <Riostream.h>
#include <TCanvas.h>
#include <TChain.h>
#include <TClonesArray.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <THnSparse.h>
#include <TList.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TString.h>
#include <TTree.h>
#include <TRandom3.h>

#include "AliAODEvent.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskSE.h"
#include "AliEMCALGeometry.h"
#include "AliEPFlattener.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDtrackCuts.h"
#include "AliEventplane.h"
#include "AliMCEvent.h"
#include "AliOADBContainer.h"
#include "AliStack.h"
#include "AliVCluster.h"
#include "AliCalorimeterUtils.h"
#include "AliFiducialCut.h" // Needed for detector flag enum kEMCAL, kPHOS
#include "AliAODCaloCluster.h"

using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskPi0V2)

//-----------------------------------------------------------------
AliAnalysisTaskPi0V2::AliAnalysisTaskPi0V2(const char *name) : AliAnalysisTaskSE(name)
{
  fDebug=0;

  fUseV2Cluster=kFALSE;
  fUseV1Cluster=kFALSE;
  fUseTrk=kFALSE;
  fV2ClusterName="";
  fV1ClusterName="";
  fTrackName="";
  fV2Cluster=NULL;
  fV1Cluster=NULL;
  fTrack=NULL;

  fOutput=NULL;
  fAODEvent=NULL;
  fGeom=NULL;
  fGeomName="EMCAL_COMPLETEV1";
  fPhosEPCaliContainer=NULL;
  fPhosEPCaliFileName="$ALICE_PHYSICS/OADB/PHOS/PHOSflat.root";
  fUsePhosEPCali=kTRUE;
  fEventPlane=NULL;
  fCaloUtils=NULL;
  fRunNum=-999;
  fInterRunNum=-999;
  fVzCut=10.;
  fVzBin=-1;
  fCentMin=0;
  fCentMax=100;
  fCentDetector="V0M";
  fCentrality=99.;
  fCentBin=-1;
  fFlattenSemiCent=kFALSE;

  fNCellCut=2.;
  fECut=1.;
  fEtaCut=0.65;
  fV2M02Cut=0.5;
  fV1M02Cut=0.3;
  fDrCut=0.025;
  fPi0AsyCut=kFALSE;
  fNLMCutMin=1;
  fNLMCutMax=10;
  fApplySSCut=kTRUE;
  fSplitV1Cluster=kFALSE;
  for (int i = 0; i < 10; ++i) {
    for (int j = 0; j < 2; ++j) {
      fBufferIndex[i][j] = 0;

      for (int k = 0; k < 10; ++k) {
        fBufferSplitV1[i][j][k].SetPxPyPzE(0,0,0,0);
      }
    }
  }

  fEPTPC=-999.;
  fEPTPCReso=0.;
  fEPV0=-999.;
  fEPV0A=-999.;
  fEPV0C=-999.;
  fEPV0AR=-999.;
  fEPV0CR=-999.;
  fEPV0R=-999.;
  fEPV0AR4=-999.;
  fEPV0AR5=-999.;
  fEPV0AR6=-999.;
  fEPV0AR7=-999.;
  fEPV0CR0=-999.;
  fEPV0CR1=-999.;
  fEPV0CR2=-999.;
  fEPV0CR3=-999.;
  fEPTPCFlat=NULL;
  fEPV0AFlat=NULL;
  fEPV0CFlat=NULL;

  hEvtCount=NULL;
  hCentA=NULL;
  hCentB=NULL;

  hEPTPC=NULL;
  hEPTPCReso=NULL;
  hEPV0A=NULL;
  hEPV0C=NULL;
  hEPTPCFlat=NULL;
  hEPV0AFlat=NULL;
  hEPV0CFlat=NULL;
  hEPV0=NULL;
  hEPV0AR=NULL;
  hEPV0CR=NULL;
  hEPV0R=NULL;
  hEPV0AR4=NULL;
  hEPV0AR7=NULL;
  hEPV0CR0=NULL;
  hEPV0CR3=NULL;
  hEPDiffV0A_V0CR0=NULL;
  hEPDiffV0A_V0CR3=NULL;
  hEPDiffV0CR0_V0CR3=NULL;
  hEPDiffV0C_V0AR4=NULL;
  hEPDiffV0C_V0AR7=NULL;
  hEPDiffV0AR4_V0AR7=NULL;
  hEPRbrCosV0A=NULL;
  hEPRbrSinV0A=NULL;
  hEPRbrCosV0C=NULL;
  hEPRbrSinV0C=NULL;
  hEPRbrCosTPC=NULL;
  hEPRbrSinTPC=NULL;

  hV2ClusterDxDzA=NULL;
  hV2ClusterDxDzB=NULL;
  hV2ClusterDphiV0A=NULL;
  hV2ClusterDphiV0C=NULL;
  hV2ClusterCos2phiV0A=NULL;
  hV2ClusterCos2phiV0C=NULL;
  hV1ClusterDxDzA=NULL;
  hV1ClusterDxDzB=NULL;
  hV1ClusterM02EA=NULL;
  hV1ClusterM02EB=NULL;
  hV1ClusterNlmA=NULL;
  hV1ClusterNlmB=NULL;

  hTrkPhiEta=NULL;
  hTrkPt=NULL;
  hTrkDphiEmcV0A=NULL;
  hTrkDphiEmcV0C=NULL;
  hTrkCos2phiEmcV0A=NULL;
  hTrkCos2phiEmcV0C=NULL;
  hTrkDphiOutEmcV0A=NULL;
  hTrkDphiOutEmcV0C=NULL;
  hTrkCos2phiOutEmcV0A=NULL;
  hTrkCos2phiOutEmcV0C=NULL;

  fV2ClusterV0A=NULL;
  fV2ClusterV0C=NULL;
  fV2ClusterTPC=NULL;
  fV1ClusterV0A=NULL;
  fV1ClusterV0C=NULL;
  fV1ClusterTPC=NULL;

  // Dummy constructor ALWAYS needed for I/O.
  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
}

//-----------------------------------------------------------------
AliAnalysisTaskPi0V2::AliAnalysisTaskPi0V2() : AliAnalysisTaskSE("default_name")
{
  fDebug=0;

  fUseV2Cluster=kFALSE;
  fUseV1Cluster=kFALSE;
  fUseTrk=kFALSE;
  fV2ClusterName="";
  fV1ClusterName="";
  fTrackName="";
  fV2Cluster=NULL;
  fV1Cluster=NULL;
  fTrack=NULL;

  fOutput=NULL;
  fAODEvent=NULL;
  fGeom=NULL;
  fGeomName="EMCAL_COMPLETEV1";
  fPhosEPCaliContainer=NULL;
  fPhosEPCaliFileName="$ALICE_PHYSICS/OADB/PHOS/PHOSflat.root";
  fUsePhosEPCali=kTRUE;
  fEventPlane=NULL;
  fCaloUtils=NULL;
  fRunNum=-999;
  fInterRunNum=-999;
  fVzCut=10.;
  fVzBin=-1;
  fCentMin=0;
  fCentMax=100;
  fCentDetector="V0M";
  fCentrality=99.;
  fCentBin=-1;
  fFlattenSemiCent=kFALSE;

  fNCellCut=2.;
  fECut=1.;
  fEtaCut=0.65;
  fV2M02Cut=0.5;
  fV1M02Cut=0.3;
  fDrCut=0.025;
  fPi0AsyCut=kFALSE;
  fNLMCutMin=1;
  fNLMCutMax=10;
  fApplySSCut=kTRUE;
  fSplitV1Cluster=kFALSE;
  for (int i = 0; i < 10; ++i) {
    for (int j = 0; j < 2; ++j) {
      fBufferIndex[i][j] = 0;

      for (int k = 0; k < 10; ++k) {
        fBufferSplitV1[i][j][k].SetPxPyPzE(0,0,0,0);
      }
    }
  }

  fEPTPC=-999.;
  fEPTPCReso=0.;
  fEPV0=-999.;
  fEPV0A=-999.;
  fEPV0C=-999.;
  fEPV0AR=-999.;
  fEPV0CR=-999.;
  fEPV0R=-999.;
  fEPV0AR4=-999.;
  fEPV0AR5=-999.;
  fEPV0AR6=-999.;
  fEPV0AR7=-999.;
  fEPV0CR0=-999.;
  fEPV0CR1=-999.;
  fEPV0CR2=-999.;
  fEPV0CR3=-999.;
  fEPTPCFlat=NULL;
  fEPV0AFlat=NULL;
  fEPV0CFlat=NULL;

  hEvtCount=NULL;
  hCentA=NULL;
  hCentB=NULL;

  hEPTPC=NULL;
  hEPTPCReso=NULL;
  hEPV0A=NULL;
  hEPV0C=NULL;
  hEPTPCFlat=NULL;
  hEPV0AFlat=NULL;
  hEPV0CFlat=NULL;
  hEPV0=NULL;
  hEPV0AR=NULL;
  hEPV0CR=NULL;
  hEPV0R=NULL;
  hEPV0AR4=NULL;
  hEPV0AR7=NULL;
  hEPV0CR0=NULL;
  hEPV0CR3=NULL;
  hEPDiffV0A_V0CR0=NULL;
  hEPDiffV0A_V0CR3=NULL;
  hEPDiffV0CR0_V0CR3=NULL;
  hEPDiffV0C_V0AR4=NULL;
  hEPDiffV0C_V0AR7=NULL;
  hEPDiffV0AR4_V0AR7=NULL;
  hEPRbrCosV0A=NULL;
  hEPRbrSinV0A=NULL;
  hEPRbrCosV0C=NULL;
  hEPRbrSinV0C=NULL;
  hEPRbrCosTPC=NULL;
  hEPRbrSinTPC=NULL;

  hV2ClusterDxDzA=NULL;
  hV2ClusterDxDzB=NULL;
  hV2ClusterDphiV0A=NULL;
  hV2ClusterDphiV0C=NULL;
  hV2ClusterCos2phiV0A=NULL;
  hV2ClusterCos2phiV0C=NULL;
  hV1ClusterDxDzA=NULL;
  hV1ClusterDxDzB=NULL;
  hV1ClusterM02EA=NULL;
  hV1ClusterM02EB=NULL;
  hV1ClusterNlmA=NULL;
  hV1ClusterNlmB=NULL;

  hTrkPhiEta=NULL;
  hTrkPt=NULL;
  hTrkDphiEmcV0A=NULL;
  hTrkDphiEmcV0C=NULL;
  hTrkCos2phiEmcV0A=NULL;
  hTrkCos2phiEmcV0C=NULL;
  hTrkDphiOutEmcV0A=NULL;
  hTrkDphiOutEmcV0C=NULL;
  hTrkCos2phiOutEmcV0A=NULL;
  hTrkCos2phiOutEmcV0C=NULL;

  fV2ClusterV0A=NULL;
  fV2ClusterV0C=NULL;
  fV2ClusterTPC=NULL;
  fV1ClusterV0A=NULL;
  fV1ClusterV0C=NULL;
  fV1ClusterTPC=NULL;

  // Constructor
  // Define input and output slots here (never in the dummy constructor)
  // Input slot #0 works with a TChain - it is connected to the default input container
  // Output slot #1 writes into a TH1 container
  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
}

//-----------------------------------------------------------------
AliAnalysisTaskPi0V2::~AliAnalysisTaskPi0V2()
{
  // Destructor. Clean-up the output list, but not the histograms that are put inside
  // (the list is owner and will clean-up these histograms). Protect in PROOF case.
  if (fEPTPCFlat) delete fEPTPCFlat;
  fEPTPCFlat=0x0;
  if (fEPV0AFlat) delete fEPV0AFlat;
  fEPV0AFlat=0x0;
  if (fEPV0CFlat) delete fEPV0CFlat;
  fEPV0CFlat=0x0;
  delete fOutput;
}

//-----------------------------------------------------------------
void AliAnalysisTaskPi0V2::UserCreateOutputObjects()
{
  // Create histograms
  // Called once (on the worker node)
  fGeom = AliEMCALGeometry::GetInstance(fGeomName.Data());
  if (!fGeom) {
    AliError(Form("%s: Could not retrieve AliEMCALGeometry", GetName()));
    return;
  }

  fPhosEPCaliContainer = new AliOADBContainer("phosFlat");
  fPhosEPCaliContainer->InitFromFile(fPhosEPCaliFileName.Data(), "phosFlat");
  if (fUsePhosEPCali && !fPhosEPCaliContainer) {
    AliError(Form("%s: Could not retrieve PhosEPCaliFile", GetName()));
    return;
  }

  fOutput = new TList();
  fOutput->SetOwner();  // IMPORTANT!

  hEvtCount = new TH1F("hEvtCount", " Event Plane", 9, 0.5, 9.5);
  hEvtCount->GetXaxis()->SetBinLabel(1,"All");
  hEvtCount->GetXaxis()->SetBinLabel(2,"Evt.");
  hEvtCount->GetXaxis()->SetBinLabel(3,"Trg.");
  hEvtCount->GetXaxis()->SetBinLabel(4,"Vtx.");
  hEvtCount->GetXaxis()->SetBinLabel(5,"Cent.");
  hEvtCount->GetXaxis()->SetBinLabel(6,"EP.");
  hEvtCount->GetXaxis()->SetBinLabel(7,"Clust.");
  hEvtCount->GetXaxis()->SetBinLabel(8,"Trk.");
  fOutput->Add(hEvtCount);

  hCentA = new TH1F("hCentA", "Centrality dist. before App. flat cut", 100, 0., 100.);
  hCentB = new TH1F("hCentB", "Centrality dist. after App. flat cut", 100, 0., 100.);
  fOutput->Add(hCentA);
  fOutput->Add(hCentB);

  hEPTPC     = new TH2F("hEPTPC",     "EPTPC     vs Cent.", 100, 0., 100., 100, 0., TMath::Pi());
  hEPTPCReso = new TH2F("hEPTPCReso", "TPC Reso. vs Cent.", 100, 0., 100., 100, 0., 1.);
  hEPV0A     = new TH2F("hEPV0A",     "EPV0A     vs Cent.", 100, 0., 100., 100, 0., TMath::Pi());
  hEPV0C     = new TH2F("hEPV0C",     "EPV0C     vs Cent.", 100, 0., 100., 100, 0., TMath::Pi());
  hEPTPCFlat = new TH2F("hEPTPCFlat", "EPTPC vs Cent. after PHOS Flattening", 100, 0., 100., 100, 0., TMath::Pi());
  hEPV0AFlat = new TH2F("hEPV0AFlat", "EPV0A vs Cent. after PHOS Flattening", 100, 0., 100., 100, 0., TMath::Pi());
  hEPV0CFlat = new TH2F("hEPV0CFlat", "EPV0C vs Cent. after PHOS Flattening", 100, 0., 100., 100, 0., TMath::Pi());
  fOutput->Add(hEPTPC);
  fOutput->Add(hEPTPCReso);
  fOutput->Add(hEPV0A);
  fOutput->Add(hEPV0C);
  fOutput->Add(hEPTPCFlat);
  fOutput->Add(hEPV0AFlat);
  fOutput->Add(hEPV0CFlat);

  hEPV0    = new TH2F("hEPV0",    "EPV0    vs Cent.", 100, 0., 100., 100, 0., TMath::Pi());
  hEPV0AR  = new TH2F("hEPV0AR",  "EPV0AR  vs Cent.", 100, 0., 100., 100, 0., TMath::Pi());
  hEPV0CR  = new TH2F("hEPV0CR",  "EPV0CR  vs Cent.", 100, 0., 100., 100, 0., TMath::Pi());
  hEPV0R   = new TH2F("hEPV0R",   "EPV0R   vs Cent.", 100, 0., 100., 100, 0., TMath::Pi());
  hEPV0AR4 = new TH2F("hEPV0AR4", "EPV0AR4 vs Cent.", 100, 0., 100., 100, 0., TMath::Pi());
  hEPV0AR7 = new TH2F("hEPV0AR7", "EPV0AR7 vs Cent.", 100, 0., 100., 100, 0., TMath::Pi());
  hEPV0CR0 = new TH2F("hEPV0CR0", "EPV0CR0 vs Cent.", 100, 0., 100., 100, 0., TMath::Pi());
  hEPV0CR3 = new TH2F("hEPV0CR3", "EPV0CR3 vs Cent.", 100, 0., 100., 100, 0., TMath::Pi());
  fOutput->Add(hEPV0);
  fOutput->Add(hEPV0AR);
  fOutput->Add(hEPV0CR);
  fOutput->Add(hEPV0R);
  fOutput->Add(hEPV0AR4);
  fOutput->Add(hEPV0AR7);
  fOutput->Add(hEPV0CR0);
  fOutput->Add(hEPV0CR3);

  hEPDiffV0A_V0CR0    = new TH2F("hEPDiffV0A_V0CR0",    "EP A-R0 ",  100, 0., 100., 100, -1., 1.);    
  hEPDiffV0A_V0CR3    = new TH2F("hEPDiffV0A_V0CR3",    "EP A-R3 ",  100, 0., 100., 100, -1., 1.);    
  hEPDiffV0CR0_V0CR3  = new TH2F("hEPDiffV0CR0_V0CR3",  "EP R0-R3 ", 100, 0., 100., 100, -1., 1.);    
  hEPDiffV0C_V0AR4    = new TH2F("hEPDiffV0C_V0AR4",    "EP C-R4 ",  100, 0., 100., 100, -1., 1.);    
  hEPDiffV0C_V0AR7    = new TH2F("hEPDiffV0C_V0AR7",    "EP C-R7 ",  100, 0., 100., 100, -1., 1.);    
  hEPDiffV0AR4_V0AR7  = new TH2F("hEPDiffV0AR4_V0AR7",  "EP R4-R7 ", 100, 0., 100., 100, -1., 1.);    
  fOutput->Add(hEPDiffV0A_V0CR0);
  fOutput->Add(hEPDiffV0A_V0CR3);
  fOutput->Add(hEPDiffV0CR0_V0CR3);
  fOutput->Add(hEPDiffV0C_V0AR4);
  fOutput->Add(hEPDiffV0C_V0AR7);
  fOutput->Add(hEPDiffV0AR4_V0AR7);

  hEPRbrCosV0A = new TProfile2D("hEPRbrCosV0A", "", 100, 0, 100, 200, 0, 200, -1., 1.);
  hEPRbrSinV0A = new TProfile2D("hEPRbrSinV0A", "", 100, 0, 100, 200, 0, 200, -1., 1.);
  hEPRbrCosV0C = new TProfile2D("hEPRbrCosV0C", "", 100, 0, 100, 200, 0, 200, -1., 1.);
  hEPRbrSinV0C = new TProfile2D("hEPRbrSinV0C", "", 100, 0, 100, 200, 0, 200, -1., 1.);
  hEPRbrCosTPC = new TProfile2D("hEPRbrCosTPC", "", 100, 0, 100, 200, 0, 200, -1., 1.);
  hEPRbrSinTPC = new TProfile2D("hEPRbrSinTPC", "", 100, 0, 100, 200, 0, 200, -1., 1.);
  fOutput->Add(hEPRbrCosV0A);
  fOutput->Add(hEPRbrSinV0A);
  fOutput->Add(hEPRbrCosV0C);
  fOutput->Add(hEPRbrSinV0C);
  fOutput->Add(hEPRbrCosTPC);
  fOutput->Add(hEPRbrSinTPC);

  if (fUseV2Cluster) {
    hV2ClusterDxDzA = new TH2F("hV2ClusterDxDzA", "clust. Dx vs Dz before cut", 1000, -1., 1., 1000, -1., 1);  
    hV2ClusterDxDzB = new TH2F("hV2ClusterDxDzB", "clust. Dx vs Dz after cut",  1000, -1., 1., 1000, -1., 1);
    fOutput->Add(hV2ClusterDxDzA);
    fOutput->Add(hV2ClusterDxDzB);

    hV2ClusterDphiV0A = new TH3F("hV2ClusterDphiV0A", "V2 clust. dphi with EPV0A", 100, 0., 100., 100, 0., TMath::Pi(), 15, 0., 15.);
    hV2ClusterDphiV0C = new TH3F("hV2ClusterDphiV0C", "V2 clust. dphi with EPV0C", 100, 0., 100., 100, 0., TMath::Pi(), 15, 0., 15.);
    fOutput->Add(hV2ClusterDphiV0A);
    fOutput->Add(hV2ClusterDphiV0C);
    
    hV2ClusterCos2phiV0A = new TH3F("hV2ClusterCos2phiV0A", "V2 Clust. raw v2 with EPV0A", 100, 0, 100, 50, -1., 1., 15, 0., 15.);
    hV2ClusterCos2phiV0C = new TH3F("hV2ClusterCos2phiV0C", "V2 Clust. raw v2 with EPV0C", 100, 0, 100, 50, -1., 1., 15, 0., 15.);
    fOutput->Add(hV2ClusterCos2phiV0A);
    fOutput->Add(hV2ClusterCos2phiV0C);

    const Int_t n   = 5;
    Int_t   bins[n] = {500, 40,  100, 315,  200};
    Double_t min[n] = {0,   0.,  0,   0.,   -1.};
    Double_t max[n] = {0.5, 40., 100, 3.15, 1.};

    fV2ClusterV0A = new THnSparseF("fV2ClusterV0A", "Flow histogram EPV0A", n, bins, min, max);
    fV2ClusterV0C = new THnSparseF("fV2ClusterV0C", "Flow histogram EPV0C", n, bins, min, max);
    fV2ClusterTPC = new THnSparseF("fV2ClusterTPC", "Flow histogram EPTPC", n, bins, min, max);
    
    fV2ClusterV0A->GetAxis(0)->SetTitle("m_{#gamma#gamma}"); 
    fV2ClusterV0A->GetAxis(1)->SetTitle("p_{T}"); 
    fV2ClusterV0A->GetAxis(2)->SetTitle("Centrality");
    fV2ClusterV0A->GetAxis(3)->SetTitle("#Delta#phi");
    fV2ClusterV0A->GetAxis(4)->SetTitle("cos(2#Delta#phi)");

    fV2ClusterV0C->GetAxis(0)->SetTitle("m_{#gamma#gamma}"); 
    fV2ClusterV0C->GetAxis(1)->SetTitle("p_{T}"); 
    fV2ClusterV0C->GetAxis(2)->SetTitle("Centrality");
    fV2ClusterV0C->GetAxis(3)->SetTitle("#Delta#phi");
    fV2ClusterV0C->GetAxis(4)->SetTitle("cos(2#Delta#phi)");

    fV2ClusterTPC->GetAxis(0)->SetTitle("m_{#gamma#gamma}"); 
    fV2ClusterTPC->GetAxis(1)->SetTitle("p_{T}"); 
    fV2ClusterTPC->GetAxis(2)->SetTitle("Centrality");
    fV2ClusterTPC->GetAxis(3)->SetTitle("#Delta#phi");
    fV2ClusterTPC->GetAxis(4)->SetTitle("cos(2#Delta#phi)");

    fOutput->Add(fV2ClusterV0A);
    fOutput->Add(fV2ClusterV0C);
    fOutput->Add(fV2ClusterTPC);
  }

  if (fUseV1Cluster) {
    hV1ClusterDxDzA = new TH2F("hV1ClusterDxDzA", "clust. Dx vs Dz before cut", 1000, -1., 1., 1000, -1., 1);  
    hV1ClusterDxDzB = new TH2F("hV1ClusterDxDzB", "clust. Dx vs Dz after cut",  1000, -1., 1., 1000, -1., 1);
    fOutput->Add(hV1ClusterDxDzA);
    fOutput->Add(hV1ClusterDxDzB);

    hV1ClusterM02EA = new TH2F("hV1ClusterM02EA", "M02 vs E before SS cut", 5000, 0, 50, 500, 0, 5.);
    hV1ClusterM02EB = new TH2F("hV1ClusterM02EB", "M02 vs E after SS cut",  5000, 0, 50, 500, 0, 5.);
    fOutput->Add(hV1ClusterM02EA);
    fOutput->Add(hV1ClusterM02EB);

    hV1ClusterNlmA = new TH1F("hV1ClusterNlmA", "", 10, 0, 10);
    hV1ClusterNlmB = new TH1F("hV1ClusterNlmB", "", 10, 0, 10);
    fOutput->Add(hV1ClusterNlmA);
    fOutput->Add(hV1ClusterNlmB);

    if (!fSplitV1Cluster) {
      const Int_t n    = 6;
      Int_t   bins[n]  = {40,  40, 350, 100, 315,  200};
      Double_t min[n]  = {0.,  0., 0.,  0,   0.,  -1.};
      Double_t max[n]  = {40., 40, 3.5, 100, 3.15, 1.};
  
      fV1ClusterV0A = new THnSparseF("fV1ClusterV0A", "", n, bins, min, max);
      fV1ClusterV0C = new THnSparseF("fV1ClusterV0C", "", n, bins, min, max);
      fV1ClusterTPC = new THnSparseF("fV1ClusterTPC", "", n, bins, min, max);
  
      fV1ClusterV0A->GetAxis(0)->SetTitle("p_{T} [GeV]");
      fV1ClusterV0A->GetAxis(1)->SetTitle("E [GeV]"); 
      fV1ClusterV0A->GetAxis(2)->SetTitle("M02"); 
      fV1ClusterV0A->GetAxis(3)->SetTitle("Centrality");
      fV1ClusterV0A->GetAxis(4)->SetTitle("#Delta#phi"); 
      fV1ClusterV0A->GetAxis(5)->SetTitle("cos(2#Delta#phi)");
  
      fV1ClusterV0C->GetAxis(0)->SetTitle("p_{T} [GeV]"); 
      fV1ClusterV0C->GetAxis(1)->SetTitle("E [GeV]"); 
      fV1ClusterV0C->GetAxis(2)->SetTitle("M02"); 
      fV1ClusterV0C->GetAxis(3)->SetTitle("Centrality");
      fV1ClusterV0C->GetAxis(4)->SetTitle("#Delta#phi"); 
      fV1ClusterV0C->GetAxis(5)->SetTitle("cos(2#Delta#phi)");
  
      fV1ClusterTPC->GetAxis(0)->SetTitle("p_{T} [GeV]"); 
      fV1ClusterTPC->GetAxis(1)->SetTitle("E [GeV]"); 
      fV1ClusterTPC->GetAxis(2)->SetTitle("M02"); 
      fV1ClusterTPC->GetAxis(3)->SetTitle("Centrality");
      fV1ClusterTPC->GetAxis(4)->SetTitle("#Delta#phi"); 
      fV1ClusterTPC->GetAxis(5)->SetTitle("cos(2#Delta#phi)");
  
      fOutput->Add(fV1ClusterV0A);
      fOutput->Add(fV1ClusterV0C);
      fOutput->Add(fV1ClusterTPC);
    } else {      
      const Int_t n   = 5;
      Int_t   bins[n] = {500, 40,  100, 315,  200};
      Double_t min[n] = {0,   0.,  0,   0.,   -1.};
      Double_t max[n] = {0.5, 40., 100, 3.15, 1.};

      fV1ClusterV0A = new THnSparseF("fV1ClusterV0A", "", n, bins, min, max);
      fV1ClusterV0C = new THnSparseF("fV1ClusterV0C", "", n, bins, min, max);
      fV1ClusterTPC = new THnSparseF("fV1ClusterTPC", "", n, bins, min, max);
      
      fV1ClusterV0A->GetAxis(0)->SetTitle("m_{#gamma#gamma}"); 
      fV1ClusterV0A->GetAxis(1)->SetTitle("p_{T}"); 
      fV1ClusterV0A->GetAxis(2)->SetTitle("Centrality");
      fV1ClusterV0A->GetAxis(3)->SetTitle("#Delta#phi");
      fV1ClusterV0A->GetAxis(4)->SetTitle("cos(2#Delta#phi)");

      fV1ClusterV0C->GetAxis(0)->SetTitle("m_{#gamma#gamma}"); 
      fV1ClusterV0C->GetAxis(1)->SetTitle("p_{T}"); 
      fV1ClusterV0C->GetAxis(2)->SetTitle("Centrality");
      fV1ClusterV0C->GetAxis(3)->SetTitle("#Delta#phi");
      fV1ClusterV0C->GetAxis(4)->SetTitle("cos(2#Delta#phi)");

      fV1ClusterTPC->GetAxis(0)->SetTitle("m_{#gamma#gamma}"); 
      fV1ClusterTPC->GetAxis(1)->SetTitle("p_{T}"); 
      fV1ClusterTPC->GetAxis(2)->SetTitle("Centrality");
      fV1ClusterTPC->GetAxis(3)->SetTitle("#Delta#phi");
      fV1ClusterTPC->GetAxis(4)->SetTitle("cos(2#Delta#phi)");

      fOutput->Add(fV1ClusterV0A);
      fOutput->Add(fV1ClusterV0C);
      fOutput->Add(fV1ClusterTPC);
    }
  }

  if (fUseTrk) {
    hTrkPhiEta = new TH2F("hTrkPhiEta","ch. trk. #phi vs #eta", 100, -1.0, 1.0, 100, 0.0, 6.29);
    hTrkPt     = new TH1F("hTrkPt",    "ch. trk. p_{T}", 100, 0.0, 30.0);
    fOutput->Add(hTrkPhiEta);
    fOutput->Add(hTrkPt);

    hTrkDphiEmcV0A = new TH3F("hTrkDphiEmcV0A", "ch. trk. dphi in EMC with EPV0A", 100, 0., 100., 100, 0., TMath::Pi(), 15, 0., 15.);
    hTrkDphiEmcV0C = new TH3F("hTrkDphiEmcV0C", "ch. trk. dphi in EMC with EPV0C", 100, 0., 100., 100, 0., TMath::Pi(), 15, 0., 15.);
    fOutput->Add(hTrkDphiEmcV0A);
    fOutput->Add(hTrkDphiEmcV0C);
    
    hTrkCos2phiEmcV0A = new TH3F("hTrkCos2phiEmcV0A", "ch. trk. raw v2 in EMC with EPV0A", 100, 0, 100, 50, -1., 1., 15, 0., 15.);
    hTrkCos2phiEmcV0C = new TH3F("hTrkCos2phiEmcV0C", "ch. trk. raw v2 in EMC with EPV0C", 100, 0, 100, 50, -1., 1., 15, 0., 15.);
    fOutput->Add(hTrkCos2phiEmcV0A);
    fOutput->Add(hTrkCos2phiEmcV0C);

    hTrkDphiOutEmcV0A = new TH3F("hTrkDphiOutEmcV0A", "ch. trk. dphi NOT in EMC with EPV0A", 100, 0., 100., 100, 0., TMath::Pi(), 15, 0., 15.);
    hTrkDphiOutEmcV0C = new TH3F("hTrkDphiOutEmcV0C", "ch. trk. dphi NOT in EMC with EPV0C", 100, 0., 100., 100, 0., TMath::Pi(), 15, 0., 15.);
    fOutput->Add(hTrkDphiOutEmcV0A);
    fOutput->Add(hTrkDphiOutEmcV0C);

    hTrkCos2phiOutEmcV0A = new TH3F("hTrkCos2phiOutEmcV0A", "ch. trk. raw v2 NOT in EMC with EPV0A", 100, 0, 100, 50, -1., 1., 15, 0., 15.);
    hTrkCos2phiOutEmcV0C = new TH3F("hTrkCos2phiOutEmcV0C", "ch. trk. raw v2 NOT in EMC with EPV0A", 100, 0, 100, 50, -1., 1., 15, 0., 15.);
    fOutput->Add(hTrkCos2phiOutEmcV0A);
    fOutput->Add(hTrkCos2phiOutEmcV0C);
  }


  PostData(1, fOutput); // Post data for ALL output slots >0 here, to get at least an empty histogram
}

//-----------------------------------------------------------------
void AliAnalysisTaskPi0V2::UserExec(Option_t *) 
{
  // Main loop
  // Called for each event

  hEvtCount->Fill(1);

  AliVEvent* event = InputEvent();
  if (!event) {
    AliError(Form("%s: Could not retrieve event", GetName()));
    return;
  }

  TString type = AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()->GetDataType();
  if (type!="AOD") {
    AliError(Form("%s: Could not handle non-AOD event", GetName()));
    return;
  } else {
    fAODEvent = dynamic_cast<AliAODEvent*>(event);
    if (!fAODEvent) {
      AliError(Form("%s: Could not retrieve AOD event", GetName()));
      return;
    }
  }

  if (!fCaloUtils) {
    AliError(Form("%s: Could not retrieve CaloUtils", GetName()));
  }
  fCaloUtils->AccessGeometry(fAODEvent);
  // fCaloUtils->AccessOADB(fAODEvent);

  hEvtCount->Fill(2);

  //--------------------------------
  // event
  //--------------------------------

  // trig class
  hEvtCount->Fill(3);

  // run num.
  fRunNum = event->GetRunNumber();
  fInterRunNum = ConvertToInternalRunNum(fRunNum);

  // vertex
  const AliVVertex* fvertex;
  fvertex = event->GetPrimaryVertex();
  if (fabs(fvertex->GetZ()) > fVzCut) return;
  if (fabs(fvertex->GetZ() - event->GetPrimaryVertexSPD()->GetZ()) > 0.1) return;
  Double_t vertex[3] = {fvertex->GetX(), fvertex->GetY(), fvertex->GetZ()};

  fVzBin = -1;
  if (fvertex->GetZ() < 0) fVzBin = 0;
  else fVzBin = 1;

  hEvtCount->Fill(4);

  // centrality
  if      (fCentDetector=="V0M") fCentrality = event->GetCentrality()->GetCentralityPercentile("V0M");
  else if (fCentDetector=="CL1") fCentrality = event->GetCentrality()->GetCentralityPercentile("CL1");
  else if (fCentDetector=="TRK") fCentrality = event->GetCentrality()->GetCentralityPercentile("TRK");
  else {
   AliError(Form("%s: No such centrality definition", GetName()));
   return; 
  }
  hCentA->Fill(fCentrality);
  if (!IsCentAccepted()) return;
  hCentB->Fill(fCentrality);

  fCentBin = -1;
  for (int i = 0; i < 10; ++i) {
    if (fCentrality>i*10) {
      fCentBin++;
    } else break;
  }

  hEvtCount->Fill(5);

  // event plane
  TObjArray *maps = (TObjArray*)fPhosEPCaliContainer->GetObject(fRunNum, "phosFlat");
  fEPTPCFlat = (AliEPFlattener*)maps->At(0);
  fEPV0AFlat = (AliEPFlattener*)maps->At(1);
  fEPV0CFlat = (AliEPFlattener*)maps->At(2);
  fEventPlane = event->GetEventplane();
  if (fEventPlane) {
    VZEROEventPlane(fUsePhosEPCali);
    hEvtCount->Fill(6);
  } else {
    AliError(Form("%s: Could not retrieve the event plane", GetName()));
    return;
  }

  //--------------------------------
  // V2 Cluster
  //--------------------------------

  if (fUseV2Cluster) {
    fV2Cluster = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fV2ClusterName));
    if (!fV2Cluster) {
      AliError(Form("%s: Could not retrieve V2 cluster %s", GetName(), fV2ClusterName.Data()));
    } else {
      Int_t nCluster = fV2Cluster->GetEntries();
      if (fDebug) cout << "Num. of V2 Cluster = " << nCluster << endl;

      for (Int_t i=0; i<nCluster; ++i) {
        AliVCluster* c1 = static_cast<AliVCluster*>(fV2Cluster->At(i));
        if (!c1) continue;
        if (!c1->IsEMCAL()) continue;
  
        hV2ClusterDxDzA->Fill(c1->GetTrackDz(), c1->GetTrackDx());
        if (!IsGoodCluster(c1)) continue;
        if (c1->GetM02() > fV2M02Cut) continue;
        hV2ClusterDxDzB->Fill(c1->GetTrackDz(), c1->GetTrackDx());
  
        TLorentzVector p1;
        GetMom(p1, c1, vertex);
        Double_t cPhi = p1.Phi();
        Double_t cPt  = p1.Pt();
        Double_t dphiV0A = cPhi-fEPV0A;
        Double_t dphiV0C = cPhi-fEPV0C;
        if (dphiV0A<0.) dphiV0A+=TMath::Pi(); if(dphiV0A>TMath::Pi()) dphiV0A-=TMath::Pi();
        if (dphiV0C<0.) dphiV0C+=TMath::Pi(); if(dphiV0C>TMath::Pi()) dphiV0C-=TMath::Pi();
        hV2ClusterDphiV0A->Fill(fCentrality, dphiV0A, cPt);
        hV2ClusterDphiV0C->Fill(fCentrality, dphiV0A, cPt);
        hV2ClusterCos2phiV0A ->Fill(fCentrality, TMath::Cos(2.*dphiV0A), cPt);
        hV2ClusterCos2phiV0C ->Fill(fCentrality, TMath::Cos(2.*dphiV0C), cPt);
  
        for (Int_t j=i+1; j<nCluster; ++j) {
         AliVCluster* c2 = static_cast<AliVCluster*>(fV2Cluster->At(j));      
         if (!c2) continue;
         if (!c2->IsEMCAL()) continue;
         if (!IsGoodCluster(c2)) continue;
         if (c2->GetM02() > fV2M02Cut) continue;

         TLorentzVector p2;
         GetMom(p2, c2, vertex);
         FillPionFromV2(p1, p2);
        }
      }
    }
  }

  //--------------------------------
  // V1 cluster
  //--------------------------------

  if (fUseV1Cluster) {
    fV1Cluster = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fV1ClusterName));
    if (!fV1Cluster) {
      AliError(Form("%s: Could not retrieve V1 cluster %s", GetName(), fV1ClusterName.Data()));
      return;
    } 

    Int_t nCluster = fV1Cluster->GetEntries();
    if (fDebug) cout << "Num. of V1 Cluster = " << nCluster << endl;

    for (Int_t i=0; i<nCluster; ++i) {
      AliVCluster* cluster = dynamic_cast<AliVCluster*>(fV1Cluster->At(i));      
      if (!cluster) continue;
      if (!cluster->IsEMCAL()) continue;

      hV1ClusterDxDzA->Fill(cluster->GetTrackDz(), cluster->GetTrackDx());
      if (!IsGoodCluster(cluster)) continue;
      hV1ClusterDxDzB->Fill(cluster->GetTrackDz(), cluster->GetTrackDx());

      double M02 = cluster->GetM02();
      int    nlm = cluster->GetNExMax();
      double   E = cluster->E();
      int     nc = cluster->GetNCells();

      hV1ClusterNlmA->Fill(nlm);
      hV1ClusterM02EA->Fill(E, M02);
      if (M02<fV1M02Cut) continue;
      if (nlm<fNLMCutMin || nlm>fNLMCutMax) continue;
      if (fApplySSCut && !PassPi0SSCut(cluster)) continue;
      hV1ClusterNlmB->Fill(nlm);
      hV1ClusterM02EB->Fill(E, M02);

      //--------------------------------
      // fill pion from V1 cluster
      //--------------------------------

      if (!fSplitV1Cluster) {
        TLorentzVector p;
        GetMom(p, cluster, vertex);
        FillPionFromV1(p, cluster);
      } else {

        //--------------------------------
        // split cluster by CaloUtils
        //--------------------------------
        // copy from CaloPID

        AliVCaloCells* caloCells = (AliVCaloCells*)fAODEvent->GetEMCALCells();
        if (!caloCells) continue;

        int   absIdList[nc];
        float maxEList [nc];
  
        int nMax = fCaloUtils->GetNumberOfLocalMaxima(cluster, caloCells, absIdList, maxEList); 
        if (nMax!=2) continue;
  
        int absId1 = absIdList[0];
        int absId2 = absIdList[1];
  
        // order in energy
        Int_t  calorimeter = AliFiducialCut::kEMCAL;
  
        float en1 = caloCells->GetCellAmplitude(absId1);
        fCaloUtils->RecalibrateCellAmplitude(en1,calorimeter,absId1);
        float en2 = caloCells->GetCellAmplitude(absId2);
        fCaloUtils->RecalibrateCellAmplitude(en2,calorimeter,absId2);
        if (en1 < en2) {
          absId2 = absIdList[0];
          absId1 = absIdList[1];
        }
        if (absId2<0 || absId1<0) continue;
  
        AliAODCaloCluster cluster1(0, 0,NULL,0.,NULL,NULL,1,0);
        AliAODCaloCluster cluster2(1, 0,NULL,0.,NULL,NULL,1,0);
        
        fCaloUtils->SplitEnergy(absId1, absId2, cluster, caloCells, &cluster1, &cluster2, nMax);
        // caloutils->GetEMCALRecoUtils()->RecalculateClusterDistanceToBadChannel(caloutils->GetEMCALGeometry(),caloCells,&cluster1);
        // caloutils->GetEMCALRecoUtils()->RecalculateClusterDistanceToBadChannel(caloutils->GetEMCALGeometry(),caloCells,&cluster2);
  
        TLorentzVector fMomentum1;
        TLorentzVector fMomentum2;
        cluster1.GetMomentum(fMomentum1, vertex);
        cluster2.GetMomentum(fMomentum2, vertex);

        //--------------------------------
        // fill buffer for bg
        //--------------------------------

        if (fBufferIndex[fCentBin][fVzBin]%10==0) fBufferIndex[fCentBin][fVzBin] = 0;
        fBufferSplitV1[fCentBin][fVzBin][fBufferIndex[fCentBin][fVzBin]] = fMomentum2;
        fBufferIndex[fCentBin][fVzBin]++;

        //--------------------------------
        // fill pion from split V1 cluster
        //--------------------------------

        FillPionFromSplitV1(fMomentum1);
      }
    }
  }
  hEvtCount->Fill(7);

  //--------------------------------
  // Track
  //--------------------------------

  if (fUseTrk) {
    fTrack = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fTrackName));
    if (!fTrack) {
      AliError(Form("%s: Could not retrieve tracks %s", GetName(), fTrackName.Data())); 
    } else {
      Int_t nTracks = fTrack->GetEntries();
      for (Int_t i=0; i<nTracks; ++i){
        AliVTrack* track = static_cast<AliVTrack*>(fTrack->At(i));
        if (!track) continue;
        Double_t tPhi = track->Phi();
        Double_t tPt  = track->Pt();
        Double_t tEta  = track->Eta();
        hTrkPhiEta->Fill(tEta, tPhi);
        hTrkPt->Fill(tPt);
  
        Double_t dphiV0A = TVector2::Phi_0_2pi(tPhi-fEPV0A); if (dphiV0A>TMath::Pi()) dphiV0A-=TMath::Pi();
        Double_t dphiV0C = TVector2::Phi_0_2pi(tPhi-fEPV0C); if (dphiV0C>TMath::Pi()) dphiV0C-=TMath::Pi();
  
        if (fabs(tEta)>0.7) continue;
        if (tPhi*TMath::RadToDeg()>80. && tPhi*TMath::RadToDeg()<180.) {
          hTrkDphiEmcV0A->Fill(fCentrality, dphiV0A, tPt);
          hTrkDphiEmcV0C->Fill(fCentrality, dphiV0C, tPt);
          hTrkCos2phiEmcV0A->Fill(fCentrality, TMath::Cos(2.*dphiV0A), tPt);
          hTrkCos2phiEmcV0C->Fill(fCentrality, TMath::Cos(2.*dphiV0C), tPt);
        } else {
          hTrkDphiOutEmcV0A->Fill(fCentrality, dphiV0A, tPt);
          hTrkDphiOutEmcV0C->Fill(fCentrality, dphiV0C, tPt);
          hTrkCos2phiOutEmcV0A->Fill(fCentrality, TMath::Cos(2.*dphiV0A), tPt);
          hTrkCos2phiOutEmcV0C->Fill(fCentrality, TMath::Cos(2.*dphiV0C), tPt);
        }
      }
    }
  }
  hEvtCount->Fill(8);

  // NEW HISTO should be filled before this point, as PostData puts the
  // information for this iteration of the UserExec in the container
  PostData(1, fOutput);
}

//-----------------------------------------------------------------
Bool_t AliAnalysisTaskPi0V2::IsCentAccepted()
{
  if (fCentrality<fCentMin || fCentrality>fCentMax) return kFALSE;

  if (10<fCentrality && fCentrality<=50) { // 10-50%
    if (fFlattenSemiCent) {
      TString centfired = fAODEvent->GetFiredTriggerClasses();
      if (!centfired.Contains("CVLN_B2-B-NOPF-ALLNOTRD") && 
          !centfired.Contains("CVLN_R1-B-NOPF-ALLNOTRD") && 
          !centfired.Contains("CSEMI_R1-B-NOPF-ALLNOTRD")) return kFALSE;
      else return kTRUE;
    } else return kTRUE;
  } else return kTRUE; // other%

  return kFALSE;
}

//-----------------------------------------------------------------
void AliAnalysisTaskPi0V2::VZEROEventPlane(Bool_t flattenEP)
{ // Calculate the V0 Event Plane
  if (fEventPlane->GetQVector()) { 
    // fEPTPC = TVector2::Phi_0_2pi(fEventPlane->GetQVector()->Phi())/2.0; if(fEPTPC>TMath::Pi()) fEPTPC-=TMath::Pi();
    fEPTPC = fEventPlane->GetEventplane("Q");
  } else fEPTPC = -999;

  if (fEventPlane->GetQsub1() && fEventPlane->GetQsub2()) {
    // fEPTPCReso = TMath::Cos(2.0*(fEventPlane->GetQsub1()->Phi()/2.0-fEventPlane->GetQsub2()->Phi()/2.0));
    fEPTPCReso = TMath::Cos(2.*(fEventPlane->GetQsubRes()));
  } else fEPTPCReso = -1;

  fEPV0  = TVector2::Phi_0_2pi(fEventPlane->GetEventplane("V0",  fAODEvent)); if (fEPV0>TMath::Pi()) fEPV0-=TMath::Pi();
  fEPV0A = TVector2::Phi_0_2pi(fEventPlane->GetEventplane("V0A", fAODEvent)); if (fEPV0A>TMath::Pi()) fEPV0A-=TMath::Pi();
  fEPV0C = TVector2::Phi_0_2pi(fEventPlane->GetEventplane("V0C", fAODEvent)); if (fEPV0C>TMath::Pi()) fEPV0C-=TMath::Pi();

  if (fDebug==2) cout << "EPTPC/EPV0A/EPV0C = " << fEPTPC << "/" << fEPV0A << "/" << fEPV0C << endl;

  Double_t qx=0, qy=0, qxr=0, qyr=0;
  fEPV0AR = TVector2::Phi_0_2pi(fEventPlane->CalculateVZEROEventPlane(fAODEvent, 4, 5, 2, qxr, qyr)); if (fEPV0AR>TMath::Pi()) fEPV0AR-=TMath::Pi();
  fEPV0CR = TVector2::Phi_0_2pi(fEventPlane->CalculateVZEROEventPlane(fAODEvent, 2, 3, 2, qx,  qy)); if (fEPV0CR>TMath::Pi()) fEPV0CR-=TMath::Pi();
  qxr += qx; qyr += qy;
  fEPV0R   = TVector2::Phi_0_2pi(TMath::ATan2(qyr,qxr))/2.0; // equals to ring 2-5
  fEPV0AR4 = TVector2::Phi_0_2pi(fEventPlane->CalculateVZEROEventPlane(fAODEvent, 4, 2, qx, qy)); if (fEPV0AR4>TMath::Pi()) fEPV0AR4-=TMath::Pi();
  fEPV0AR5 = TVector2::Phi_0_2pi(fEventPlane->CalculateVZEROEventPlane(fAODEvent, 5, 2, qx, qy)); if (fEPV0AR5>TMath::Pi()) fEPV0AR5-=TMath::Pi();
  fEPV0AR6 = TVector2::Phi_0_2pi(fEventPlane->CalculateVZEROEventPlane(fAODEvent, 6, 2, qx, qy)); if (fEPV0AR6>TMath::Pi()) fEPV0AR6-=TMath::Pi();
  fEPV0AR7 = TVector2::Phi_0_2pi(fEventPlane->CalculateVZEROEventPlane(fAODEvent, 7, 2, qx, qy)); if (fEPV0AR7>TMath::Pi()) fEPV0AR7-=TMath::Pi();
  fEPV0CR0 = TVector2::Phi_0_2pi(fEventPlane->CalculateVZEROEventPlane(fAODEvent, 0, 2, qx, qy)); if (fEPV0CR0>TMath::Pi()) fEPV0CR0-=TMath::Pi();
  fEPV0CR1 = TVector2::Phi_0_2pi(fEventPlane->CalculateVZEROEventPlane(fAODEvent, 1, 2, qx, qy)); if (fEPV0CR1>TMath::Pi()) fEPV0CR1-=TMath::Pi();
  fEPV0CR2 = TVector2::Phi_0_2pi(fEventPlane->CalculateVZEROEventPlane(fAODEvent, 2, 2, qx, qy)); if (fEPV0CR2>TMath::Pi()) fEPV0CR2-=TMath::Pi();
  fEPV0CR3 = TVector2::Phi_0_2pi(fEventPlane->CalculateVZEROEventPlane(fAODEvent, 3, 2, qx, qy)); if (fEPV0CR3>TMath::Pi()) fEPV0CR3-=TMath::Pi();

  hEPTPC->Fill(fCentrality,  fEPTPC); 
  if (fEPTPCReso!=-1) hEPTPCReso->Fill(fCentrality, fEPTPCReso);
  hEPV0->Fill(fCentrality, fEPV0);
  hEPV0A->Fill(fCentrality, fEPV0A);
  hEPV0C->Fill(fCentrality, fEPV0C);
  hEPV0AR->Fill(fCentrality, fEPV0AR);
  hEPV0CR->Fill(fCentrality, fEPV0CR);
  hEPV0R->Fill(fCentrality, fEPV0R);
  hEPV0AR4->Fill(fCentrality, fEPV0AR4);
  hEPV0AR7->Fill(fCentrality, fEPV0AR7);
  hEPV0CR0->Fill(fCentrality, fEPV0CR0);
  hEPV0CR3->Fill(fCentrality, fEPV0CR3);

  if (flattenEP) {
    fEPV0A = FlattenV0A(fEPV0A, fCentrality);
    fEPV0C = FlattenV0C(fEPV0C, fCentrality);
    if (fEPTPC != -999.) fEPTPC = FlattenTPC(fEPTPC, fCentrality);
    fEPV0A = TVector2::Phi_0_2pi(fEPV0A); if (fEPV0A>TMath::Pi()) fEPV0A-=TMath::Pi();
    fEPV0C = TVector2::Phi_0_2pi(fEPV0C); if (fEPV0C>TMath::Pi()) fEPV0C-=TMath::Pi();
    fEPTPC = TVector2::Phi_0_2pi(fEPTPC); if (fEPTPC>TMath::Pi()) fEPTPC-=TMath::Pi();

    hEPTPCFlat->Fill(fCentrality, fEPTPC);
    hEPV0AFlat->Fill(fCentrality, fEPV0A);
    hEPV0CFlat->Fill(fCentrality, fEPV0C);
  }

  hEPDiffV0A_V0CR0->Fill(fCentrality, TMath::Cos(2.0*(fEPV0A - fEPV0CR0)));
  hEPDiffV0A_V0CR3->Fill(fCentrality, TMath::Cos(2.0*(fEPV0A - fEPV0CR3)));
  hEPDiffV0CR0_V0CR3->Fill(fCentrality, TMath::Cos(2.0*(fEPV0CR0 - fEPV0CR3)));
  hEPDiffV0C_V0AR4->Fill(fCentrality, TMath::Cos(2.0*(fEPV0C - fEPV0AR4)));
  hEPDiffV0C_V0AR7->Fill(fCentrality, TMath::Cos(2.0*(fEPV0C - fEPV0AR7)));
  hEPDiffV0AR4_V0AR7->Fill(fCentrality, TMath::Cos(2.0*(fEPV0AR4 - fEPV0AR7)));   

  // run-by-run QA
  hEPRbrCosV0A->Fill(fCentrality, fInterRunNum, TMath::Cos(2*fEPV0A));
  hEPRbrSinV0A->Fill(fCentrality, fInterRunNum, TMath::Sin(2*fEPV0A));
  hEPRbrCosV0C->Fill(fCentrality, fInterRunNum, TMath::Cos(2*fEPV0C));
  hEPRbrSinV0C->Fill(fCentrality, fInterRunNum, TMath::Sin(2*fEPV0C));
  hEPRbrCosTPC->Fill(fCentrality, fInterRunNum, TMath::Cos(2*fEPTPC));
  hEPRbrSinTPC->Fill(fCentrality, fInterRunNum, TMath::Sin(2*fEPTPC));
}

//-----------------------------------------------------------------
Double_t  AliAnalysisTaskPi0V2::FlattenV0A(Double_t phi, Double_t c)
{
  if (fEPV0AFlat) return fEPV0AFlat->MakeFlat(phi,c);
  return phi;
}

Double_t  AliAnalysisTaskPi0V2::FlattenV0C(Double_t phi, Double_t c)
{ 
  if (fEPV0CFlat) return fEPV0CFlat->MakeFlat(phi,c);
  return phi;
}

Double_t  AliAnalysisTaskPi0V2::FlattenTPC(Double_t phi, Double_t c)
{
  if (fEPTPCFlat) return fEPTPCFlat->MakeFlat(phi,c);
  return phi;
}

//-----------------------------------------------------------------
Bool_t AliAnalysisTaskPi0V2::IsGoodCluster(const AliVCluster* c) const
{
  if (!c) return kFALSE;
  if (c->GetNCells() < fNCellCut) return kFALSE;
  if (c->E() < fECut) return kFALSE;

  Short_t id = -1;
  Double_t maxE = GetMaxCellEnergy(c, id); 
  if ((1. - double(GetCrossEnergy(c,id))/maxE) > 0.97) return kFALSE;

  Float_t pos[3] = {0,0,0};
  c->GetPosition(pos);
  TVector3 cPos(pos);
  Double_t eta = cPos.Eta();

  if (fabs(eta) > fEtaCut) return kFALSE;  
  if (!IsWithinFiducialVolume(id)) return kFALSE;

  Double_t dr = TMath::Sqrt(c->GetTrackDx()*c->GetTrackDx() + c->GetTrackDz()*c->GetTrackDz());
  if (dr < fDrCut) return kFALSE;

  return kTRUE;
}

Double_t AliAnalysisTaskPi0V2::GetMaxCellEnergy(const AliVCluster *c, Short_t &id) const
{
  AliVCaloCells* caloCells = (AliVCaloCells*)fAODEvent->GetEMCALCells();
  if (!caloCells) return 0;

  // Get maximum energy of attached cell.
  id = -1;
  Double_t maxe = 0;
  const Int_t ncells = c->GetNCells();

  for (Int_t i=0; i<ncells; i++) {
    Double_t e = caloCells->GetCellAmplitude(fabs(c->GetCellAbsId(i)));
    if (e>maxe) {
      maxe = e;
      id = c->GetCellAbsId(i);
    }
  }
  return maxe;
}

Double_t AliAnalysisTaskPi0V2::GetCrossEnergy(const AliVCluster *c, Short_t &idmax) const
{ 
  AliVCaloCells* caloCells = (AliVCaloCells*)fAODEvent->GetEMCALCells();
  if (!caloCells) return 0;

  // Calculate the energy of cross cells around the leading cell.
  Int_t iSupMod = -1;
  Int_t iTower  = -1;
  Int_t iIphi   = -1;
  Int_t iIeta   = -1;
  Int_t iphi    = -1;
  Int_t ieta    = -1;
  Int_t iphis   = -1;
  Int_t ietas   = -1;

  Double_t crossEnergy = 0.;
  fGeom->GetCellIndex(idmax,iSupMod,iTower,iIphi,iIeta);
  fGeom->GetCellPhiEtaIndexInSModule(iSupMod,iTower,iIphi, iIeta,iphis,ietas);

  Int_t ncells = c->GetNCells();
  for (Int_t i=0; i<ncells; i++) {
    Int_t cellAbsId = c->GetCellAbsId(i);
    fGeom->GetCellIndex(cellAbsId,iSupMod,iTower,iIphi,iIeta);
    fGeom->GetCellPhiEtaIndexInSModule(iSupMod,iTower,iIphi, iIeta,iphi,ieta);
    Int_t aphidiff = fabs(iphi-iphis);
    if (aphidiff>1) continue;
    Int_t aetadiff = fabs(ieta-ietas);
    if (aetadiff>1) continue;
    if ((aphidiff==1 && aetadiff==0) ||
        (aphidiff==0 && aetadiff==1) ) {
      crossEnergy += caloCells->GetCellAmplitude(cellAbsId);
    }
  }

  return crossEnergy;
}

Bool_t AliAnalysisTaskPi0V2::IsWithinFiducialVolume(Short_t id) const
{ 
  // Check if cell is within given fiducial volume.
  Double_t fNFiducial = 1;

  Int_t iSupMod = -1;
  Int_t iTower  = -1;
  Int_t iIphi   = -1;
  Int_t iIeta   = -1;
  Int_t iphi    = -1;
  Int_t ieta    = -1;

  Bool_t okrow = kFALSE;
  Bool_t okcol = kFALSE;

  Int_t cellAbsId = id;
  fGeom->GetCellIndex(cellAbsId,iSupMod,iTower,iIphi,iIeta);
  fGeom->GetCellPhiEtaIndexInSModule(iSupMod,iTower,iIphi, iIeta,iphi,ieta);

  // Check rows/phi
  if (iSupMod < 10) {
    if (iphi >= fNFiducial && iphi < 24-fNFiducial) okrow = kTRUE;
  } else {
    if (iphi >= fNFiducial && iphi < 12-fNFiducial) okrow = kTRUE;
  }
  // Check columns/eta
  Bool_t noEMCALBorderAtEta0 = kTRUE;
  if (!noEMCALBorderAtEta0) {
    if (ieta > fNFiducial && ieta < 48-fNFiducial) okcol = kTRUE;
  } else {
    if (iSupMod%2==0) {
      if (ieta >= fNFiducial) okcol = kTRUE;
    } else {
      if (ieta < 48-fNFiducial) okcol = kTRUE;
    }
  }
  if (okrow && okcol) return kTRUE;

  return kFALSE;
}

//-----------------------------------------------------------------
Bool_t AliAnalysisTaskPi0V2::PassPi0SSCut(const AliVCluster *c)
{
  Double_t E = c->E();
  Double_t M02 = c->GetM02();
  Int_t nlm = c->GetNExMax();

  Double_t M02Min = exp(2.135-0.245*E);
  Double_t M02Max = 0.;
  if      (nlm == 1) M02Max = exp(0.0662-0.0201*E) - 0.0955 + 0.00186*E + 9.91/E;
  else if (nlm == 2) M02Max = exp(0.353-0.0264*E) - 0.524 + 0.00559*E + 21.9/E;
  if (nlm > 2) M02Max += 0.75;

  if (M02>M02Max || M02<M02Min) return kFALSE;
  if (fDebug) cout << "E/M02/NLM = " << E << "/" << M02 << "/" << nlm 
                   << " (" << M02Min << ", " << M02Max << ")" << endl;

  return kTRUE;
}

//-----------------------------------------------------------------
Bool_t AliAnalysisTaskPi0V2::IsInPi0SplitAsymmetryRange(Float_t energy, Float_t asy, Int_t nlm) const
{
  Float_t fAsyMinParam[2][4];
  fAsyMinParam[0][0] = 0.96 ;
  fAsyMinParam[0][1] = 0    ;
  fAsyMinParam[0][2] =-879  ;
  fAsyMinParam[0][3] = 0.96 ; // Absolute max

  //TF1 *lAsyNLM2 = new TF1("lAsyNLM2","0.95+0.0015*x-233/(x*x*x)",5,100);
  fAsyMinParam[1][0] = 0.95  ;
  fAsyMinParam[1][1] = 0.0015;
  fAsyMinParam[1][2] =-233   ;
  fAsyMinParam[1][3] = 1.0   ; // Absolute max

  float abasy = TMath::Abs(asy);

  int inlm = nlm-1;
  if (nlm > 2) inlm=1; // only 2 cases defined nlm=1 and nlm>=2
  
  // Get the parametrized min cut of asymmetry for NLM=2 up to 11 GeV
  float cut = fAsyMinParam[inlm][0] + fAsyMinParam[inlm][1]*energy + fAsyMinParam[inlm][2]/energy/energy/energy ;
  
  // In any case and beyond validity energy range of the function,
  // the parameter cannot be smaller than 1
  if (cut > fAsyMinParam[inlm][3]) cut = fAsyMinParam[inlm][3];
  
  //printf("energy %2.2f - nlm: %d (%d)- p0 %f, p1 %f, p2 %f, p3 %f ; cut: %2.2f\n",energy,nlm,inlm,
  //       fAsyMinParam[inlm][0],fAsyMinParam[inlm][1],fAsyMinParam[inlm][2],fAsyMinParam[inlm][3],cut);
  
  if (abasy < cut) return kTRUE;
  else             return kFALSE;
}

//-----------------------------------------------------------------
void AliAnalysisTaskPi0V2::FillPionFromV2(const TLorentzVector& p1, const TLorentzVector& p2)
{
  if (fPi0AsyCut) {
    Double_t asym = fabs(p1.E()-p2.E())/(p1.E()+p2.E());
    if (asym>0.7) return;
  }

  TLorentzVector pion;
  pion = p1 + p2;
  Double_t eta = pion.Eta();
  if (fabs(eta) > fEtaCut) return;

  Double_t mass = pion.M();
  Double_t pt   = pion.Pt();
  Double_t phi  = pion.Phi();

  Double_t dphiV0A = TVector2::Phi_0_2pi(phi-fEPV0A); if(dphiV0A>TMath::Pi()) {dphiV0A-=TMath::Pi();}
  Double_t dphiV0C = TVector2::Phi_0_2pi(phi-fEPV0C); if(dphiV0C>TMath::Pi()) {dphiV0C-=TMath::Pi();}
  Double_t dphiTPC = TVector2::Phi_0_2pi(phi-fEPTPC); if(dphiTPC>TMath::Pi()) {dphiTPC-=TMath::Pi();}

  Double_t DataV0A[5];
  DataV0A[0] = mass;
  DataV0A[1] = pt;
  DataV0A[2] = fCentrality;
  DataV0A[3] = dphiV0A;
  DataV0A[4] = TMath::Cos(2.*(dphiV0A));
  fV2ClusterV0A->Fill(DataV0A);

  Double_t DataV0C[5];
  DataV0C[0] = mass;
  DataV0C[1] = pt;
  DataV0C[2] = fCentrality;
  DataV0C[3] = dphiV0C;
  DataV0C[4] = TMath::Cos(2.*(dphiV0C));
  fV2ClusterV0C->Fill(DataV0C);

  Double_t DataTPC[5];
  DataTPC[0] = mass;
  DataTPC[1] = pt;
  DataTPC[2] = fCentrality;
  DataTPC[3] = dphiTPC;
  DataTPC[4] = TMath::Cos(2.*(dphiTPC));
  fV2ClusterTPC->Fill(DataTPC);
}

//-----------------------------------------------------------------
void AliAnalysisTaskPi0V2::FillPionFromV1(const TLorentzVector& p1, AliVCluster *c)
{
  Double_t Et   = p1.Et();
  Double_t phi  = p1.Phi();
  Double_t M02  = c->GetM02();
  Double_t E    = c->E();

  Double_t dphiV0A = TVector2::Phi_0_2pi(phi-fEPV0A); if (dphiV0A>TMath::Pi()) dphiV0A-=TMath::Pi();
  Double_t dphiV0C = TVector2::Phi_0_2pi(phi-fEPV0C); if (dphiV0C>TMath::Pi()) dphiV0C-=TMath::Pi();
  Double_t dphiTPC = TVector2::Phi_0_2pi(phi-fEPTPC); if (dphiTPC>TMath::Pi()) dphiTPC-=TMath::Pi();

  Double_t DataV0A[6];
  DataV0A[0] = Et;
  DataV0A[1] = E;
  DataV0A[2] = M02;
  DataV0A[3] = fCentrality;
  DataV0A[4] = dphiV0A;
  DataV0A[5] = TMath::Cos(2.0*(dphiV0A));
  fV1ClusterV0A->Fill(DataV0A);

  Double_t DataV0C[6];
  DataV0C[0] = Et;
  DataV0C[1] = E;
  DataV0C[2] = M02;
  DataV0C[3] = fCentrality;
  DataV0C[4] = dphiV0C;
  DataV0C[5] = TMath::Cos(2.0*(dphiV0C));
  fV1ClusterV0C->Fill(DataV0C);

  Double_t DataTPC[6];
  DataTPC[0] = Et;
  DataTPC[1] = E;
  DataTPC[2] = M02;
  DataTPC[3] = fCentrality;
  DataTPC[4] = dphiTPC;
  DataTPC[5] = TMath::Cos(2.0*(dphiTPC));
  fV1ClusterTPC->Fill(DataTPC);
}

void AliAnalysisTaskPi0V2::FillPionFromSplitV1(const TLorentzVector& p1)
{
  Double_t DataV0A[5];
  Double_t DataV0C[5];
  Double_t DataTPC[5];

  TLorentzVector pion;

  for (int i = 0; i < 10; ++i) {
    if (!fBufferSplitV1[fCentBin][fVzBin][i].E()  ||
        !fBufferSplitV1[fCentBin][fVzBin][i].Pt() ||
        !fBufferSplitV1[fCentBin][fVzBin][i].M()) continue;

    pion = p1 + fBufferSplitV1[fCentBin][fVzBin][i];

    Double_t eta  = pion.Eta();
    Double_t mass = pion.M();
    Double_t pt   = pion.Pt();
    Double_t phi  = pion.Phi();

    if (pt<5) continue;

    Double_t dphiV0A = TVector2::Phi_0_2pi(phi-fEPV0A); if(dphiV0A>TMath::Pi()) {dphiV0A-=TMath::Pi();}
    Double_t dphiV0C = TVector2::Phi_0_2pi(phi-fEPV0C); if(dphiV0C>TMath::Pi()) {dphiV0C-=TMath::Pi();}
    Double_t dphiTPC = TVector2::Phi_0_2pi(phi-fEPTPC); if(dphiTPC>TMath::Pi()) {dphiTPC-=TMath::Pi();}

    DataV0A[0] = mass;
    DataV0A[1] = pt;
    DataV0A[2] = fCentrality;
    DataV0A[3] = dphiV0A;
    DataV0A[4] = TMath::Cos(2.*(dphiV0A));
    fV1ClusterV0A->Fill(DataV0A);

    DataV0C[0] = mass;
    DataV0C[1] = pt;
    DataV0C[2] = fCentrality;
    DataV0C[3] = dphiV0C;
    DataV0C[4] = TMath::Cos(2.*(dphiV0C));
    fV1ClusterV0C->Fill(DataV0C);

    DataTPC[0] = mass;
    DataTPC[1] = pt;
    DataTPC[2] = fCentrality;
    DataTPC[3] = dphiTPC;
    DataTPC[4] = TMath::Cos(2.*(dphiTPC));
    fV1ClusterTPC->Fill(DataTPC);
  }
}

//-----------------------------------------------------------------
void AliAnalysisTaskPi0V2::GetMom(TLorentzVector& p, const AliVCluster* c, Double_t* vertex)
{
  // Calculate momentum.
  Float_t pos[3] = {0,0,0};
  c->GetPosition(pos);
  TVector3 cPos(pos);

  Double_t e   = c->E();
  Double_t r   = cPos.Perp();
  Double_t eta = cPos.Eta();
  Double_t phi = cPos.Phi();

  TVector3 posNew;
  posNew.SetPtEtaPhi(r,eta,phi);

  if (vertex) posNew -= vertex;

  Double_t rad = posNew.Mag();
  p.SetPxPyPzE(e*posNew.x()/rad, e*posNew.y()/rad, e*posNew.z()/rad, e);
}

//-----------------------------------------------------------------
Int_t AliAnalysisTaskPi0V2::ConvertToInternalRunNum(Int_t n)
{
  switch (n) {
    case  170593 : return 179;
    case  170572 : return 178;
    case  170556 : return 177;
    case  170552 : return 176;
    case  170546 : return 175;
    case  170390 : return 174;
    case  170389 : return 173;
    case  170388 : return 172;
    case  170387 : return 171;
    case  170315 : return 170;
    case  170313 : return 169;
    case  170312 : return 168;
    case  170311 : return 167;
    case  170309 : return 166;
    case  170308 : return 165;
    case  170306 : return 164;
    case  170270 : return 163;
    case  170269 : return 162;
    case  170268 : return 161;
    case  170267 : return 160;
    case  170264 : return 159;
    case  170230 : return 158;
    case  170228 : return 157;
    case  170208 : return 156;
    case  170207 : return 155;
    case  170205 : return 154;
    case  170204 : return 153;
    case  170203 : return 152;
    case  170195 : return 151;
    case  170193 : return 150;
    case  170163 : return 149;
    case  170162 : return 148;
    case  170159 : return 147;
    case  170155 : return 146;
    case  170152 : return 145;
    case  170091 : return 144;
    case  170089 : return 143;
    case  170088 : return 142;
    case  170085 : return 141;
    case  170084 : return 140;
    case  170083 : return 139;
    case  170081 : return 138;
    case  170040 : return 137;
    case  170038 : return 136;
    case  170036 : return 135;
    case  170027 : return 134;
    case  169981 : return 133;
    case  169975 : return 132;
    case  169969 : return 131;
    case  169965 : return 130;
    case  169961 : return 129;
    case  169956 : return 128;
    case  169926 : return 127;
    case  169924 : return 126;
    case  169923 : return 125;
    case  169922 : return 124;
    case  169919 : return 123;
    case  169918 : return 122;
    case  169914 : return 121;
    case  169859 : return 120;
    case  169858 : return 119;
    case  169855 : return 118;
    case  169846 : return 117;
    case  169838 : return 116;
    case  169837 : return 115;
    case  169835 : return 114;
    case  169683 : return 113;
    case  169628 : return 112;
    case  169591 : return 111;
    case  169590 : return 110;
    case  169588 : return 109;
    case  169587 : return 108;
    case  169586 : return 107;
    case  169584 : return 106;
    case  169557 : return 105;
    case  169555 : return 104;
    case  169554 : return 103;
    case  169553 : return 102;
    case  169550 : return 101;
    case  169515 : return 100;
    case  169512 : return 99;
    case  169506 : return 98;
    case  169504 : return 97;
    case  169498 : return 96;
    case  169475 : return 95;
    case  169420 : return 94;
    case  169419 : return 93;
    case  169418 : return 92;
    case  169417 : return 91;
    case  169415 : return 90;
    case  169411 : return 89;
    case  169238 : return 88;
    case  169236 : return 87;
    case  169167 : return 86;
    case  169160 : return 85;
    case  169156 : return 84;
    case  169148 : return 83;
    case  169145 : return 82;
    case  169144 : return 81;
    case  169143 : return 80;
    case  169138 : return 79;
    case  169099 : return 78;
    case  169094 : return 77;
    case  169091 : return 76;
    case  169045 : return 75;
    case  169044 : return 74;
    case  169040 : return 73;
    case  169035 : return 72;
    case  168992 : return 71;
    case  168988 : return 70;
    case  168984 : return 69;
    case  168826 : return 68;
    case  168777 : return 67;
    case  168514 : return 66;
    case  168512 : return 65;
    case  168511 : return 64;
    case  168467 : return 63;
    case  168464 : return 62;
    case  168461 : return 61;
    case  168460 : return 60;
    case  168458 : return 59;
    case  168362 : return 58;
    case  168361 : return 57;
    case  168356 : return 56;
    case  168342 : return 55;
    case  168341 : return 54;
    case  168325 : return 53;
    case  168322 : return 52;
    case  168318 : return 51;
    case  168311 : return 50;
    case  168310 : return 49;
    case  168213 : return 48;
    case  168212 : return 47;
    case  168208 : return 46;
    case  168207 : return 45;
    case  168206 : return 44;
    case  168205 : return 43;
    case  168204 : return 42;
    case  168203 : return 41;
    case  168181 : return 40;
    case  168177 : return 39;
    case  168175 : return 38;
    case  168173 : return 37;
    case  168172 : return 36;
    case  168171 : return 35;
    case  168115 : return 34;
    case  168108 : return 33;
    case  168107 : return 32;
    case  168105 : return 31;
    case  168104 : return 30;
    case  168103 : return 29;
    case  168076 : return 28;
    case  168069 : return 27;
    case  168068 : return 26;
    case  168066 : return 25;
    case  167988 : return 24;
    case  167987 : return 23;
    case  167986 : return 22;
    case  167985 : return 21;
    case  167921 : return 20;
    case  167920 : return 19;
    case  167915 : return 18;
    case  167909 : return 17;
    case  167903 : return 16;
    case  167902 : return 15;
    case  167818 : return 14;
    case  167814 : return 13;
    case  167813 : return 12;
    case  167808 : return 11;
    case  167807 : return 10;
    case  167806 : return 9;
    case  167713 : return 8;
    case  167712 : return 7;
    case  167711 : return 6;
    case  167706 : return 5;
    case  167693 : return 4;
    case  166532 : return 3;
    case  166530 : return 2;
    case  166529 : return 1;

    default : return 199;
  }
}

//-----------------------------------------------------------------
void AliAnalysisTaskPi0V2::Terminate(Option_t *) 
{
}
