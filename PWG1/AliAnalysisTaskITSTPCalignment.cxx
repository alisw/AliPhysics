////////////////////////////////////////////////////////////////////////////////
//  AliAnalysisTaskITSTPCalignment
//  Runs the relative ITS TPC alignment procedure and TPC vdrift calib
//  Origin: Mikolaj Krzewicki, mikolaj.krzewicki@cern.ch
////////////////////////////////////////////////////////////////////////////////

#include "TChain.h"
#include "TTree.h"
#include "TH2.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliRelAlignerKalman.h"
#include "AliRelAlignerKalmanArray.h"
#include "AliAnalysisTaskITSTPCalignment.h"

ClassImp(AliAnalysisTaskITSTPCalignment)

//________________________________________________________________________
AliAnalysisTaskITSTPCalignment::AliAnalysisTaskITSTPCalignment():
    AliAnalysisTask(),
    fESD(0),
    fArray(0),
    fYZResidualsHist(0),
    fPLResidualsHist(0),
    fListOfHistos(0),
    fSaveInterval(600),
    fTimeMatchingTolerance(20),
    fDoQA(kFALSE)
{
  //dummy ctor
}

//________________________________________________________________________
AliAnalysisTaskITSTPCalignment::AliAnalysisTaskITSTPCalignment(const char *name):
    AliAnalysisTask(name,name),
    fESD(0),
    fArray(0),
    fYZResidualsHist(0),
    fPLResidualsHist(0),
    fListOfHistos(0),
    fSaveInterval(600),
    fTimeMatchingTolerance(20),
    fDoQA(kFALSE)
{
  // Constructor

  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  DefineOutput(0, AliRelAlignerKalmanArray::Class());
  DefineOutput(1, TList::Class());
}

//________________________________________________________________________
void AliAnalysisTaskITSTPCalignment::ConnectInputData(Option_t *)
{
  // Called once
  TTree* tree = dynamic_cast<TTree*> (GetInputData(0));
  if (!tree)
  {
    Printf("ERROR: Could not read chain from input slot 0");
  }
  else
  {
    AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
    if (!esdH)
    {
      Printf("ERROR: Could not get ESDInputHandler");
    }
    else
      fESD = esdH->GetEvent();
  }
}

//________________________________________________________________________
void AliAnalysisTaskITSTPCalignment::CreateOutputObjects()
{
  // Create output objects
  // Called once

  fArray = new AliRelAlignerKalmanArray("ITSTPCalignmentArray");
  fArray->SetSaveInterval(fSaveInterval);
  fArray->SetTimeMatchingTolerance(fTimeMatchingTolerance);

  fYZResidualsHist = new TH2F("fYZResidualsHist", "YZ residuals", 50, -0.5, 0.5, 50, -2., 2. );
  fPLResidualsHist = new TH2F("fPLResidualsHist", "sin(phi) tan(lambda) residuals", 50, -.05, 0.05, 50, -0.05, 0.05 );

  fListOfHistos = new TList();
  fListOfHistos->Add(fYZResidualsHist);
  fListOfHistos->Add(fPLResidualsHist);
}

//________________________________________________________________________
void AliAnalysisTaskITSTPCalignment::Exec(Option_t *)
{
  // Main loop
  // Called for each event
  if (!fESD)
  {
    Printf("ERROR: fESD not available");
    return;
  }

  AliRelAlignerKalman* aligner = fArray->GetAligner();

  Int_t lastrunnumber = aligner->GetRunNumber();
  UInt_t currentTimeStamp = fESD->GetTimeStamp();

  //for a new run reset TPC errors
  if (lastrunnumber != fESD->GetRunNumber())
  {
    aligner->ResetTPCparamsCovariance();
    fArray->SetCurrentTimeBin(currentTimeStamp);
  }

  //if time jumps back reset all
  if (currentTimeStamp < aligner->GetTimeStamp())
  {
    aligner->Reset();
    fArray->SetCurrentTimeBin(currentTimeStamp);
  }

  //Update the parmeters
  if (fArray->AddCosmicEvent(fESD))
    if (fDoQA)
    {
      //fill the QA histograms
      TArrayI trackTArrITS(1);
      TArrayI trackTArrTPC(1);
      if (aligner->FindCosmicTrackletNumbersInEvent(
            trackTArrITS, trackTArrTPC, fESD ))
      {
        AliESDtrack* ptrack;
        const AliExternalTrackParam* pconstparams1;
        const AliExternalTrackParam* pconstparams2;
        AliExternalTrackParam params1;
        AliExternalTrackParam params2;

        ////////////////////////////////
        for (Int_t i=0;i<trackTArrITS.GetSize();i++)
        {
          //ITS track
          ptrack = fESD->GetTrack(trackTArrITS[i]);
          pconstparams1 = ptrack->GetOuterParam();
          if (!pconstparams1) continue;
          params1 = *pconstparams1; //make copy to be safe

          //TPC track
          ptrack = fESD->GetTrack(trackTArrTPC[i]);
          pconstparams2 = ptrack->GetInnerParam();
          if (!pconstparams2) continue;
          params2 = *pconstparams2; //make copy
          params2.Rotate(params1.GetAlpha());
          params2.PropagateTo( params1.GetX(), aligner->GetMagField() );

          Float_t resy = params2.GetY() - params1.GetY();
          Float_t resz = params2.GetZ() - params1.GetZ();
          Float_t ressnp = params2.GetSnp() - params1.GetSnp();
          Float_t restgl = params2.GetTgl() - params1.GetTgl();
          fYZResidualsHist->Fill(resy,resz);
          fPLResidualsHist->Fill(ressnp,restgl);
        }
      }//if DoQA

    }//if AddEvent

  // Post output data.
  PostData(0, fArray);
  PostData(1, fListOfHistos);
}

//________________________________________________________________________
void AliAnalysisTaskITSTPCalignment::Terminate(Option_t *)
{
  // Called once at the end of the query
  fArray->Merge(new TList()); //final cleanup
}

