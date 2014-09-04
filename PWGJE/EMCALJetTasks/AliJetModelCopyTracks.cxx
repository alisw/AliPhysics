// $Id$
//
// Jet model task to copy tracks while making small change
// - make particles massless
//
// Author: M. Verweij

#include "AliJetModelCopyTracks.h"

#include <TClonesArray.h>
#include <TFolder.h>
#include <TLorentzVector.h>
#include <TParticle.h>
#include <TParticlePDG.h>
#include <TRandom3.h>
#include <TProfile.h>
#include <TGrid.h>
#include <TFile.h>
#include <TF1.h>
#include "AliAnalysisManager.h"
#include "AliEMCALDigit.h"
#include "AliEMCALGeometry.h"
#include "AliEMCALRecPoint.h"
#include "AliGenerator.h"
#include "AliHeader.h"
#include "AliLog.h"
#include "AliPicoTrack.h"
#include "AliRun.h"
#include "AliRunLoader.h"
#include "AliStack.h"
#include "AliStack.h"
#include "AliVCluster.h"
#include "AliVEvent.h"

ClassImp(AliJetModelCopyTracks)

//________________________________________________________________________
AliJetModelCopyTracks::AliJetModelCopyTracks() : 
AliAnalysisTaskEmcal("AliJetModelCopyTracks",kTRUE),
  fTracksOutName(""),
  fTracksOut(0x0),
  fParticleMass(kMassive),
  fHistPtOut(0)
{
  // Default constructor.
  SetMakeGeneralHistograms(kTRUE);
}

//________________________________________________________________________
AliJetModelCopyTracks::AliJetModelCopyTracks(const char *name) : 
  AliAnalysisTaskEmcal(name,kTRUE),
  fTracksOutName(""),
  fTracksOut(0x0),
  fParticleMass(kMassive),
  fHistPtOut(0)
{
  // Standard constructor.
  SetMakeGeneralHistograms(kTRUE);
}

//________________________________________________________________________
AliJetModelCopyTracks::~AliJetModelCopyTracks()
{
  // Destructor

}

//________________________________________________________________________
void AliJetModelCopyTracks::ExecOnce() 
{
  // Exec only once.

  AliAnalysisTaskEmcal::ExecOnce();

  if (!fTracksOutName.IsNull()) {
    fTracksOut = new TClonesArray("AliPicoTrack");
    fTracksOut->SetName(fTracksOutName);
    if (InputEvent()->FindListObject(fTracksOutName)) {
      AliFatal(Form("%s: Collection %s is already present in the event!", GetName(), fTracksOutName.Data()));
      return;
    }
    else {
      InputEvent()->AddObject(fTracksOut);
    }
  }
}

//________________________________________________________________________
void AliJetModelCopyTracks::UserCreateOutputObjects() 
{
  AliAnalysisTaskEmcal::UserCreateOutputObjects();

  const Int_t nBinPt = 100;
  Double_t binLimitsPt[nBinPt+1];
  for(Int_t iPt = 0;iPt <= nBinPt;iPt++){
    if(iPt == 0){
      binLimitsPt[iPt] = 0.0;
    } else {// 1.0
      binLimitsPt[iPt] =  binLimitsPt[iPt-1] + 1.0;
    }
  }

  fHistPtOut = new TH1F("fHistPtOut","fHistPtOut;#it{p}_{T};N",nBinPt,binLimitsPt);
  fOutput->Add(fHistPtOut);

  PostData(1, fOutput); // Post data for ALL output slots > 0 here.
}


//________________________________________________________________________
Bool_t AliJetModelCopyTracks::Run() 
{
  CopyTracks();
  return kTRUE;
}

//________________________________________________________________________
void AliJetModelCopyTracks::CopyTracks()
{
  //Apply toy detector simulation to tracks
  fTracksOut->Delete();
  
  Int_t nt = 0;
  const Int_t nTracks = fTracks->GetEntriesFast();
   for (Int_t i = 0; i < nTracks; ++i) {
    AliPicoTrack *picotrack = static_cast<AliPicoTrack*>(fTracks->At(i));
    if (!picotrack)
      continue;

    Double_t mass = picotrack->M();
    if(fParticleMass==kMassless) mass = 0.;
    if(fParticleMass==kPionMass) mass = 0.13957;

    AliPicoTrack *track = new ((*fTracksOut)[nt]) AliPicoTrack(picotrack->Pt(),
							       picotrack->Eta(),
							       picotrack->Phi(),
							       picotrack->Charge(),
							       picotrack->GetLabel(),
							       AliPicoTrack::GetTrackType(picotrack),
							       picotrack->GetTrackEtaOnEMCal(),
							       picotrack->GetTrackPhiOnEMCal(),
							       picotrack->GetTrackPtOnEMCal(),
							       picotrack->IsEMCAL(),
							       mass); 
    track->SetBit(TObject::kBitMask,1);
    fHistPtOut->Fill(track->Pt());
    nt++;
   }
}





