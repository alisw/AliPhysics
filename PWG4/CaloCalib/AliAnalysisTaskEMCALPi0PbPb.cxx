// $Id$

#include "AliAnalysisTaskEMCALPi0PbPb.h"
#include <TChain.h>
#include <TClonesArray.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TList.h>
#include <TLorentzVector.h>
#include "AliAODEvent.h"
#include "AliAODVertex.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskEMCALClusterizeFast.h"
#include "AliCentrality.h"
#include "AliEMCALGeoUtils.h"
#include "AliESDEvent.h"
#include "AliESDVertex.h"

ClassImp(AliAnalysisTaskEMCALPi0PbPb)

//________________________________________________________________________
AliAnalysisTaskEMCALPi0PbPb::AliAnalysisTaskEMCALPi0PbPb() 
  : AliAnalysisTaskSE(),
    fCentVar(),
    fCentFrom(0),
    fCentTo(100),
    fVtxZMin(-7),
    fVtxZMax(+7),
    fClusName(),
    fOutput(0),
    fEsdEv(0),
    fAodEv(0),
    fEsdClusters(0),
    fEsdCells(0),
    fAodClusters(0),
    fAodCells(0),
    fHcuts(0),
    fHvertexZ(0),
    fHcent(0)
{
  // ROOT constructor.
}

//________________________________________________________________________
AliAnalysisTaskEMCALPi0PbPb::AliAnalysisTaskEMCALPi0PbPb(const char *name) 
  : AliAnalysisTaskSE(name),
    fCentVar("V0M"),
    fCentFrom(0),
    fCentTo(100),
    fVtxZMin(-7),
    fVtxZMax(+7),
    fClusName(),
    fOutput(0),
    fEsdEv(0),
    fAodEv(0),
    fEsdClusters(0),
    fEsdCells(0),
    fAodClusters(0),
    fAodCells(0),
    fHcuts(0),
    fHvertexZ(0),
    fHcent(0)
{
  // Constructor.

  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
  fBranchNames="ESD:AliESDRun.,AliESDHeader.,PrimaryVertex,EMCALCells.,CaloClusters "
               "AOD:header,vertices,emcalCells,caloClusters";
}

//________________________________________________________________________
AliAnalysisTaskEMCALPi0PbPb::~AliAnalysisTaskEMCALPi0PbPb()
{
  // Destructor.

  delete fOutput; fOutput = 0;
}

//________________________________________________________________________
void AliAnalysisTaskEMCALPi0PbPb::UserCreateOutputObjects()
{
  // Create user objects here.

  fOutput = new TList();
  fOutput->SetOwner();

  fHcuts = new TH1F("hCuts","",3,0.5,3.5);
  fHcuts->GetXaxis()->SetBinLabel(1,"All (PS)");
  fHcuts->GetXaxis()->SetBinLabel(2,Form("%s: %.0f-%.0f",fCentVar.Data(),fCentFrom,fCentTo));
  fHcuts->GetXaxis()->SetBinLabel(3,Form("zvtx: %.0f-%.0f",fVtxZMin,fVtxZMax));
  fOutput->Add(fHcuts);
  fHvertexZ = new TH1F("hVertexZ",";z [cm];",100,-25,25);
  fOutput->Add(fHvertexZ);
  fHcent = new TH1F("hCent",Form(";%s;",fCentVar.Data()),101,-1,100);
  fOutput->Add(fHcent);

  PostData(1, fOutput); 
}

//________________________________________________________________________
void AliAnalysisTaskEMCALPi0PbPb::UserExec(Option_t *) 
{
  // Called for each event.

  if (!InputEvent())
    return;

  Int_t cut = 1;
  fHcuts->Fill(cut++);

  const AliCentrality *centP = InputEvent()->GetCentrality();
  Double_t cent = centP->GetCentralityPercentileUnchecked(fCentVar);
  fHcent->Fill(cent);
  if (cent<fCentFrom||cent>fCentTo)
    return;
  
  // todo test quality flag (optional)

  fHcuts->Fill(cut++);

  const AliVVertex *vertex =  InputEvent()->GetPrimaryVertex();
  if (!vertex)
    return;

  fHvertexZ->Fill(vertex->GetZ());

  if(vertex->GetZ()<fVtxZMin||vertex->GetZ()>fVtxZMax)
    return;

  fHcuts->Fill(cut++);

  AliAnalysisManager *am = AliAnalysisManager::GetAnalysisManager();
  fEsdEv = dynamic_cast<AliESDEvent*>(InputEvent());
  fAodEv = dynamic_cast<AliAODEvent*>(InputEvent());

  fEsdClusters = 0;
  fEsdCells    = 0;
  fAodClusters = 0;
  fAodCells    = 0;

  // deal with AOD output first
  if (AODEvent() && !fClusName.IsNull()) {
    TList *l = AODEvent()->GetList();
    TClonesArray *clus = 0;
    if (l) {
      clus = dynamic_cast<TClonesArray*>(l->FindObject(fClusName));
      fAodClusters = clus;
    }
  }

  if (fEsdEv) { // ESD input mode
    if (!fEsdClusters) {
      AliAnalysisTaskEMCALClusterizeFast *cltask = 0;
      TObjArray *ts = am->GetTasks();
      for (Int_t i=0; i<ts->GetEntries(); ++i) {
        cltask = dynamic_cast<AliAnalysisTaskEMCALClusterizeFast*>(ts->At(i));
        if (cltask && cltask->GetClusters()) {
          fEsdClusters = cltask->GetClusters();
          break;
        }
      }
    }
    if (!fEsdClusters) {
      am->LoadBranch("CaloClusters");
      TList *l = fEsdEv->GetList();
      TClonesArray *clus = 0;
      if (l) {
        clus = dynamic_cast<TClonesArray*>(l->FindObject("CaloClusters"));
        fEsdClusters = clus;
      }
    }
    if (!fEsdCells) {
      am->LoadBranch("EMCALCells.");
      fEsdCells = fEsdEv->GetEMCALCells();
    }
  } else if (fAodEv) { // AOD input mode
    if (!fAodClusters) {
      am->LoadBranch("caloClusters");
      TList *l = fAodEv->GetList();
      TClonesArray *clus = 0;
      if (l) {
        clus = dynamic_cast<TClonesArray*>(l->FindObject("caloClusters"));
        fAodClusters = clus;
      }
    }
    if (!fAodCells) {
      am->LoadBranch("emcalCells");
      fAodCells = fAodEv->GetEMCALCells();
    }
  } else {
    AliFatal("Impossible to not have either pointer to ESD or AOD event");
  }

  FillCellHists();
  FillClusHists();
  FillPionHists();

  PostData(1, fOutput);
}      

//________________________________________________________________________
void AliAnalysisTaskEMCALPi0PbPb::Terminate(Option_t *) 
{
  // Terminate called at the end of analysis.
}

//________________________________________________________________________
void AliAnalysisTaskEMCALPi0PbPb::FillCellHists()
{
  // Fill histograms related to cell properties.
}

//________________________________________________________________________
void AliAnalysisTaskEMCALPi0PbPb::FillClusHists()
{
  // Fill histograms related to cluster properties.
}

//________________________________________________________________________
void AliAnalysisTaskEMCALPi0PbPb::FillPionHists()
{
  // Fill histograms related to pions.
}

