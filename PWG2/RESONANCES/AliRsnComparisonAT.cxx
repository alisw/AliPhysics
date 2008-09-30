//
// Class AliRsnComparisonAT
//
// TODO
//
// authors: Martin Vala (martin.vala@cern.ch)
//          Alberto Pulvirenti (alberto.pulvirenti@ct.infn.it)
//
#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TFolder.h"
#include "TCanvas.h"
#include "TClonesArray.h"

#include "AliLog.h"

#include "AliESDtrack.h"
#include "AliMCEvent.h"

#include "AliRsnMCInfo.h"
#include "AliRsnComparisonObj.h"
#include "AliRsnComparisonAT.h"

ClassImp(AliRsnComparisonAT)

//________________________________________________________________________
AliRsnComparisonAT::AliRsnComparisonAT(const char*name)
    : AliRsnAnalysisTaskBase(name),fOutList(0x0),fMyInputNum(1),fMyPIDInputNum(0)
{
//=========================================================
// Default constructor
//=========================================================

  InitIOVars();

  DefineOutput(0, TList::Class());
}

//________________________________________________________________________
void AliRsnComparisonAT::InitIOVars()
{
//=========================================================
// Sets default values for input and output
//=========================================================
  AliDebug(AliLog::kDebug, "<-");
  AliRsnAnalysisTaskBase::InitIOVars();
  fOutList = 0;
  AliDebug(AliLog::kDebug, "->");
}

//________________________________________________________________________
void AliRsnComparisonAT::LocalInit()
{
//=========================================================
// Local init
//=========================================================
}

//________________________________________________________________________
void AliRsnComparisonAT::CreateOutputObjects()
{
//=========================================================
// CreateOutputObjects()
//=========================================================
  AliDebug(AliLog::kDebug, "<-");

//   AliRsnDaughter::SetPIDMethod ( AliRsnDaughter::kRealistic );
//   fReader = dynamic_cast<AliRsnReader*> ( GetInputData ( 1 ) );
//   fPID = dynamic_cast<AliRsnPID*> ( GetInputData ( 2 ) );


  OpenFile(0);
  fOutList = new TList();
  AliRsnComparisonObj* obj;

  TList *listMainParticle,*listMainPID,*listPID,*listParticle,*listTmp;

  listMainParticle = new TList();
  listMainParticle->SetName("MCParticles");
  for (Int_t i=0;i<fMyInput[0].GetEntries();i++)
  {
    obj = (AliRsnComparisonObj*) fMyInput[0].At(i);
    if (!obj) continue;
    listParticle = obj->GenerateParticleInfoHistogramList("");
    if (!listParticle)
    {
      AliError(Form("List not crated for i=%d",i));
      continue;
    }

    listMainParticle->Add(listParticle);
  }
  fOutList->Add(listMainParticle);
  if (fMyPIDInputNum>0)
  {
    listMainPID = new TList();
    listMainPID->SetName("PID");
    for (Int_t j=0;j<fMyPIDInputNum;j++)
    {
      listPID = new TList();
      listPID->SetName(fMyPIDInput[j].GetName());
      for (Int_t i=0;i<fMyPIDInput[j].GetEntries();i++)
      {
        obj = (AliRsnComparisonObj*) fMyPIDInput[j].At(i);
        if (!obj) continue;
        listTmp = obj->GeneratePIDHistogramList(fMyPIDInput[j].GetName());
        if (!listTmp)
        {
          AliError(Form("List not crated for i=%d",i));
          continue;
        }
        listPID->Add(listTmp);

      }
      listMainPID->Add(listPID);
    }
    fOutList->Add(listMainPID);
//     listMainPID->Print();
  }
  AliDebug(AliLog::kDebug, "->");
}

//________________________________________________________________________
void AliRsnComparisonAT::Exec(Option_t *)
{
//=========================================================
// Exec()
//=========================================================

  if (fInputType[0]==kESDMC)
  {
    fRSN[0] = GetRsnEventFromInputType(0);
    LoopOverESDtracks();
    LoopOverMCtracks();
  }


  if (!fRSN[0]) return;
  LoopOverRSNDaughters();

  if (fInputType[0]==kESDMC)
  {
    delete fRSN[0];
    fRSN[0] = 0;
  }

//   fNumEventsProcessed++;
  PostData(0, fOutList);
}

//________________________________________________________________________
void AliRsnComparisonAT::Terminate(Option_t *)
{
//=========================================================
// Terminate()
//=========================================================
  AliDebug(AliLog::kDebug, "<-");
  fOutList = dynamic_cast<TList*>(GetOutputData(0));
  if (!fOutList)
  {
    AliError(" fOutList not available");
    return;
  }

  AliDebug(AliLog::kDebug, "->");
}

//________________________________________________________________________
void AliRsnComparisonAT::LoopOverESDtracks()
{
//=========================================================
// Loop over all ESD tracks
//=========================================================
  if (!fESD[0]) return;
//   AliInfo(Form("%d",fESD[0]->GetNumberOfTracks()));
  AliRsnComparisonObj* input=0;
  for (Int_t i=0;i< fESD[0]->GetNumberOfTracks(); i++)
  {
    AliESDtrack *esdtrack = fESD[0]->GetTrack(i);
    for (Int_t j=0;j<fMyInputNum;j++)
    {
      for (Int_t i=0;i<fMyInput[j].GetEntries();i++)
      {
        input = (AliRsnComparisonObj*) fMyInput[j].At(i);
//         input->FillPIDHistograms(esdtrack, fMC[0]);
      }
    }
    for (Int_t j=0;j<fMyPIDInputNum;j++)
    {
      for (Int_t i=0;i<fMyPIDInput[j].GetEntries();i++)
      {
        input = (AliRsnComparisonObj*) fMyPIDInput[j].At(i);
        input->FillPIDHistograms(esdtrack, fMC[0]);
      }
    }
  }
}

//________________________________________________________________________
void AliRsnComparisonAT::LoopOverMCtracks()
{
//=========================================================
// Loop over all MC tracks
//=========================================================

  if (!fMC[0]) return;

  AliRsnComparisonObj* input=0;
  AliMCParticle *mctrack;
  for (Int_t i=0;i<fMC[0]->GetNumberOfTracks(); i++)
  {
    mctrack = fMC[0]->GetTrack(i);
    if (!mctrack) {AliInfo("mctrack == null");continue;}
    for (Int_t j=0;j<fMyInputNum;j++)
    {
      for (Int_t i=0;i<fMyInput[j].GetEntries();i++)
      {
        input = (AliRsnComparisonObj*) fMyInput[j].At(i);
        input->FillHistograms(mctrack);
      }
    }
    for (Int_t j=0;j<fMyPIDInputNum;j++)
    {
      for (Int_t i=0;i<fMyPIDInput[j].GetEntries();i++)
      {
        input = (AliRsnComparisonObj*) fMyPIDInput[j].At(i);
        input->FillPIDHistograms(mctrack);
      }
    }
  }
}

//________________________________________________________________________
void AliRsnComparisonAT::LoopOverRSNDaughters()
{
//=========================================================
// Loop over all rsn daughters
//=========================================================

  TClonesArray *tracks = (TClonesArray*) fRSN[0]->GetTracks();
  AliRsnDaughter *daughter=0;
  AliRsnComparisonObj* input=0;
  for (Int_t i=0;i< tracks->GetEntriesFast(); i++)
  {
    daughter = (AliRsnDaughter*) tracks->At(i);
    for (Int_t j=0;j<fMyInputNum;j++)
    {
      for (Int_t i=0;i<fMyInput[j].GetEntries();i++)
      {
        input = (AliRsnComparisonObj*) fMyInput[j].At(i);
//         input->FillPIDHistograms(daughter);
      }
    }
    for (Int_t j=0;j<fMyPIDInputNum;j++)
    {
      for (Int_t i=0;i<fMyPIDInput[j].GetEntries();i++)
      {
        input = (AliRsnComparisonObj*) fMyPIDInput[j].At(i);
        input->FillPIDHistograms(daughter);
      }
    }
  }
}

//________________________________________________________________________
void AliRsnComparisonAT::PrintStat()
{
//=========================================================
// Print stat
//=========================================================
  AliDebug(AliLog::kDebug, "<-");

  AliDebug(AliLog::kDebug, "->");
}

void AliRsnComparisonAT::AddMyInput(AliRsnComparisonObj * obj, const Int_t & index)
{
//   if ( ( index<0 ) || ( index>kMyInputNum ) ) return;

  fMyInput[index].Add(obj);
//   AliInfo(Form("%d %s",fMyInput[index].GetEntries(),fMyInput[index].GetName()));
}

void AliRsnComparisonAT::AddMyPIDInput(AliRsnComparisonObj * obj, const Int_t & index)
{
//   if ( ( index<0 ) || ( index>kMyInputNum ) ) return;

  fMyPIDInput[index].Add(obj);
//   AliInfo(Form("%d %s",fMyInput[index].GetEntries(),fMyInput[index].GetName()));
}

void AliRsnComparisonAT::SetMyPIDInputName(TString name, const Int_t & index)
{
  fMyPIDInput[index].SetName(name.Data());
}

