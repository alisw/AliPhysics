#include "TPDGCode.h"
#include "TFile.h"
#include "TTree.h"
#include "TEventList.h"
#include "TObjArray.h"
#include "TH2.h"

#include "AliLog.h"
#include "AliESDtrack.h"
#include "AliTrackReference.h"

#include "AliTRDReconstructor.h"
#include "AliTRDtrackV1.h"
#include "AliTRDseedV1.h"
#include "AliTRDpidRefMaker.h"
#include "info/AliTRDv0Info.h"


// Defines and implements basic functionality for building reference data for TRD PID. 
// 
// Here is the list of functionality provided by this class
// 1.
// 2.
// 
// Authors:
//   Alex Bercuci <A.Bercuci@gsi.de>
//   Alex Wilk    <wilka@uni-muenster.de>
//   Markus Fasel <mfasel@gsi.de>
//   Markus Heide <mheide@uni-muenster.de>
// 

ClassImp(AliTRDpidRefMaker)

//________________________________________________________________________
AliTRDpidRefMaker::AliTRDpidRefMaker(const char *name, const char *title) 
  :AliTRDrecoTask(name, title)
  ,fReconstructor(0x0)
  ,fV0s(0x0)
  ,fData(0x0)
  ,fRefPID(kMC)
  ,fRefP(kMC)
  ,fTrainFreq(1.)
  ,fTestFreq(0.)
  ,fLy(-1)
  ,fP(-1.)
{
  //
  // Default constructor
  //

  fReconstructor = new AliTRDReconstructor();
  fReconstructor->SetRecoParam(AliTRDrecoParam::GetLowFluxParam());
  memset(fdEdx, 0, 10*sizeof(Float_t));
  memset(fPID, 0, AliPID::kSPECIES*sizeof(Float_t));

  DefineInput(1, TObjArray::Class());
  DefineOutput(1, TTree::Class());
}


//________________________________________________________________________
AliTRDpidRefMaker::~AliTRDpidRefMaker() 
{
  if(fReconstructor) delete fReconstructor;
}


//________________________________________________________________________
void AliTRDpidRefMaker::ConnectInputData(Option_t *opt)
{
  AliTRDrecoTask::ConnectInputData(opt);
  fV0s = dynamic_cast<TObjArray*>(GetInputData(1));
}

//________________________________________________________________________
void AliTRDpidRefMaker::CreateOutputObjects()
{
  // Create histograms
  // Called once

  OpenFile(0, "RECREATE");
  fContainer = new TObjArray();
  fContainer->SetName(Form("Moni%s", GetName()));

  TH2 *h2 = new TH2I("hPDG","Particle abundance", AliPID::kSPECIES, -0.5, 4.5, AliTRDCalPID::kNMom, -0.5, AliTRDCalPID::kNMom-0.5);
  TAxis *ax = h2->GetXaxis();
  ax->SetNdivisions(505);
  ax->SetTitle("Particle species");
  for(Int_t is=AliPID::kSPECIES; is--;) ax->SetBinLabel(is+1, AliPID::ParticleShortName(is));
  h2->GetYaxis()->SetTitle("P bins");
  h2->GetYaxis()->SetNdivisions(511);

  fContainer->AddAt(h2, 0);
}

//________________________________________________________________________
void AliTRDpidRefMaker::Exec(Option_t *) 
{
  // Main loop
  // Called for each event

  



//   DEBUG ?!
//   for(Int_t iv0=0; iv0<fV0s->GetEntriesFast(); iv0++){
//     v0 = dynamic_cast<AliTRDv0Info*>(fV0s->At(iv0));
//     v0->Print();
//   }
  
  Int_t labelsacc[10000]; // MC labels
  memset(labelsacc, 0, sizeof(Int_t) * 10000);

  AliTRDtrackInfo     *track = 0x0;
  //AliTRDv0Info           *v0 = 0x0;
  AliTRDtrackV1    *trackTRD = 0x0;
  AliTrackReference     *ref = 0x0;
  //AliExternalTrackParam *esd = 0x0;
  AliTRDseedV1 *trackletTRD = 0x0;
  for(Int_t itrk=0, nTRD=0; itrk<fTracks->GetEntriesFast(); itrk++){
    track = (AliTRDtrackInfo*)fTracks->UncheckedAt(itrk);
    if(!track->HasESDtrack()) continue;
    ULong_t status = track->GetStatus();
    if(!(status&AliESDtrack::kTPCout)) continue;

    if(!(trackTRD = track->GetTrack())) continue; 
    //&&(track->GetNumberOfClustersRefit()

    // TOO STRONG and might introduce a bias if short 
    // tracks are to be analysed (A.Bercuci 23.09.09) 
    // use only tracks that hit 6 chambers
    //if(!(TRDtrack->GetNumberOfTracklets() == AliTRDgeometry::kNlayer)) continue;
     

    if(HasMCdata()) labelsacc[nTRD++] = track->GetLabel();

    // fill the pid information
    memset(fPID, 0, AliPID::kSPECIES*sizeof(Float_t));
    switch(fRefPID){
    case kV0: SetRefPID(kV0, track, fPID); break;
    case kMC: SetRefPID(kMC, track, fPID); break;
    case kRec: SetRefPID(kRec, trackTRD, fPID); break;
    }

    // fill the momentum and dE/dx information
    trackTRD -> SetReconstructor(fReconstructor);
    for(Int_t ily = 0; ily < AliTRDgeometry::kNlayer; ily++){
      if(!(trackletTRD = trackTRD -> GetTracklet(ily))) continue;
      if(!CookdEdx(trackletTRD)) continue;
      switch(fRefP){
      case kMC:
        if(!HasMCdata()){
          AliError("Could not retrive reference momentum from MC");
          return;
        }
        if(!(ref = track->GetTrackRef(trackletTRD))) continue;
        fP = ref->P();
        break;
      case kRec:
        fP = trackletTRD->GetMomentum();
        break;
      default:
        AliWarning("Momentum reference source not implemented");
        return;
      }
      fLy = ily;
      Fill();
    }

    for(Int_t iPart = 0; iPart < AliPID::kSPECIES; iPart++){
      AliDebug(4, Form("PDG is %d %f", iPart, fPID[iPart]));
    }
  }

  PostData(0, fContainer);
  PostData(1, fData);
}


//________________________________________________________________________
void AliTRDpidRefMaker::Fill() 
{
// Fill data tree
  fData->Fill();
}


//________________________________________________________________________
void AliTRDpidRefMaker::SetRefPID(ETRDpidRefMakerSource select, void *source, Float_t *pid) 
{
// Fill the reference PID values "pid" from "source" object
// according to the option "select". Possible options are
// - kV0  - v0 based PID
// - kMC  - MC truth [default]
// - kRec - outside detectors

  switch(select){ 
  case kV0:
    {
      AliTRDtrackInfo *track = static_cast<AliTRDtrackInfo*>(source);
      if(!track){
        AliError("No trackInfo found");
        return;
      }

      //Get V0 PID decisions from the AliTRDv0Info for all particle species (implemented so far : electrons from conversions, pions from K0s and protons from Lambdas) :
      AliTRDv0Info v0;
      for(Int_t is=AliPID::kSPECIES; is--;) fPID[is] = v0.GetV0PID(is, track);
    }
    break;
  case kMC:
    if(!HasMCdata()){
      AliError("Could not retrive reference PID from MC");
      return;
    }
    {
      AliTRDtrackInfo *track = static_cast<AliTRDtrackInfo*>(source);
      switch(track->GetPDG()){
      case kElectron:
      case kPositron:
        fPID[AliPID::kElectron] = 1.;
        break;
      case kMuonPlus:
      case kMuonMinus:
        fPID[AliPID::kMuon] = 1.;
        break;
      case kPiPlus:
      case kPiMinus:
        fPID[AliPID::kPion] = 1.;
        break;
      case kKPlus:
      case kKMinus:
        fPID[AliPID::kKaon] = 1.;
        break;
      case kProton:
      case kProtonBar:
        fPID[AliPID::kProton] = 1.;
        break;
      }
    }
    break;
  case kRec:
    { 
      AliTRDtrackV1 *trackTRD = static_cast<AliTRDtrackV1*>(source);
      trackTRD -> SetReconstructor(fReconstructor);
      //fReconstructor -> SetOption("nn");
      trackTRD -> CookPID();
      for(Int_t iPart = 0; iPart < AliPID::kSPECIES; iPart++){
        pid[iPart] = trackTRD -> GetPID(iPart);
        AliDebug(4, Form("PDG is (in V0info) %d %f", iPart, pid[iPart]));
      }
    }
    break;
  default:
    AliWarning("PID reference source not implemented");
    return;
  }
}

//________________________________________________________________________
void AliTRDpidRefMaker::SetAbundance(Float_t train, Float_t test) 
{
// Split data sample between trainning and testing

  if(fTrainFreq<0. || fTrainFreq >1. ||
     fTestFreq<0.  || fTestFreq >1.){
    AliWarning("The input data should be in the interval [0, 1]");
    return;
  }
  if(TMath::Abs(fTrainFreq+fTestFreq - 1.) > 0.001){
    AliWarning("The sum of the 2 abundances should pe one.");
    return;
  }

  fTrainFreq = train; 
  fTestFreq = test;
}

