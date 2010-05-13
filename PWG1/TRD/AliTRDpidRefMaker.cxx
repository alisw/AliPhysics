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
#include "info/AliTRDpidInfo.h"


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
AliTRDpidRefMaker::AliTRDpidRefMaker() 
  :AliTRDrecoTask()
  ,fReconstructor(NULL)
  ,fV0s(NULL)
  ,fData(NULL)
  ,fInfo(NULL)
  ,fPIDdataArray(NULL)
  ,fRefPID(kMC)
  ,fRefP(kMC)
  ,fFreq(1.)
  ,fP(-1.)
  ,fPthreshold(0.)
{
  //
  // Default constructor
  //
  SetNameTitle("PIDrefMaker", "PID Reference Maker");
}

//________________________________________________________________________
AliTRDpidRefMaker::AliTRDpidRefMaker(const char *name, const char *title) 
  :AliTRDrecoTask(name, title)
  ,fReconstructor(NULL)
  ,fV0s(NULL)
  ,fData(NULL)
  ,fInfo(NULL)
  ,fPIDdataArray(NULL)
  ,fRefPID(kMC)
  ,fRefP(kMC)
  ,fFreq(1.)
  ,fP(-1.)
  ,fPthreshold(0.5)
{
  //
  // Default constructor
  //

  fReconstructor = new AliTRDReconstructor();
  fReconstructor->SetRecoParam(AliTRDrecoParam::GetLowFluxParam());
  memset(fdEdx, 0, 10*sizeof(Float_t));
  memset(fPID, 0, AliPID::kSPECIES*sizeof(Float_t));

  DefineInput(2, TObjArray::Class()); // v0 list
  DefineInput(3, TObjArray::Class()); // pid info list 
  DefineOutput(2, TTree::Class());
}


//________________________________________________________________________
AliTRDpidRefMaker::~AliTRDpidRefMaker() 
{
  if(fPIDdataArray) delete fPIDdataArray;
  if(fReconstructor) delete fReconstructor;
}

//________________________________________________________________________
Bool_t AliTRDpidRefMaker::CheckQuality(AliTRDseedV1* /*trklt*/)
{
// Place holder for checking tracklet quality for PID.
  return kTRUE;  
}


//________________________________________________________________________
Float_t* AliTRDpidRefMaker::CookdEdx(AliTRDseedV1 *trklt)
{
  trklt->CookdEdx(AliTRDpidUtil::kNNslices);
  memcpy(fdEdx, trklt->GetdEdx(), AliTRDpidUtil::kNNslices*sizeof(Float_t));
  return fdEdx;
}


//________________________________________________________________________
void AliTRDpidRefMaker::UserCreateOutputObjects()
{
  // Create histograms
  // Called once

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

  fData = new TTree("RefPID", "Reference data for PID");
  fData->Branch("data", &fPIDdataArray);
}

//________________________________________________________________________
void AliTRDpidRefMaker::UserExec(Option_t *) 
{
  // Main loop
  // Called for each event

  if(!(fTracks = dynamic_cast<TObjArray*>(GetInputData(1)))) return;
  if(!(fV0s    = dynamic_cast<TObjArray*>(GetInputData(2)))) return;
  if(!(fInfo   = dynamic_cast<TObjArray*>(GetInputData(3)))) return;

  AliDebug(1, Form("Entries: Tracks[%d] V0[%d] PID[%d]", fTracks->GetEntriesFast(), fV0s->GetEntriesFast(), fInfo->GetEntriesFast()));
  AliTRDtrackInfo     *track = NULL;
  AliTRDtrackV1    *trackTRD = NULL;
  AliTrackReference     *ref = NULL;
  const AliTRDtrackInfo::AliESDinfo *infoESD = NULL;
  for(Int_t itrk=0; itrk<fTracks->GetEntriesFast(); itrk++){
    track = (AliTRDtrackInfo*)fTracks->UncheckedAt(itrk);
    if(!track->HasESDtrack()) continue;
    trackTRD = track->GetTrack();
    infoESD  = track->GetESDinfo();
    Double32_t *infoPID = infoESD->GetSliceIter();
    Int_t n = infoESD->GetNSlices() - AliTRDgeometry::kNlayer;
    if(n==0){
      AliWarning(Form("dEdx info missing in ESD track %d", itrk));
      continue;
    }
    Double32_t *p = &infoPID[n];
    AliDebug(4, Form("n[%d] p[GeV/c]{%6.2f %6.2f %6.2f %6.2f %6.2f %6.2f}", n, p[0], p[1], p[2], p[3], p[4], p[5]));

    ULong_t status = track->GetStatus();
    if(!(status&AliESDtrack::kTRDpid)) continue;

    // fill the pid information
    SetRefPID(fRefPID, track, fPID);
    // get particle type
    Int_t idx(TMath::LocMax(AliPID::kSPECIES, fPID)); 
    if(idx == 0 && fPID[0]<1.e-5) continue;

    // prepare PID data array
    if(!fPIDdataArray){ 
      fPIDdataArray = new AliTRDpidInfo();
    } else fPIDdataArray->Reset();
    fPIDdataArray->SetPID(idx);

    // fill PID information
    for(Int_t ily = 0; ily < AliTRDgeometry::kNlayer; ily++){

      // fill P & dE/dx information
      if(HasFriends()){ // from TRD track
        if(!trackTRD) continue;
        AliTRDseedV1 *trackletTRD(NULL);
        trackTRD -> SetReconstructor(fReconstructor);
        if(!(trackletTRD = trackTRD -> GetTracklet(ily))) continue;
        if(!CheckQuality(trackletTRD)) continue;
        if(!CookdEdx(trackletTRD)) continue;

        // fill momentum information
        fP = 0.;
        switch(fRefP){
        case kMC:
          if(!(ref = track->GetTrackRef(trackletTRD))) continue;
          fP = ref->P();
          break;
        case kRec:
          fP = trackletTRD->GetMomentum();
          break;
        default: continue;
        }
      } else { // from ESD track
        // fill momentum information
        switch(fRefP){
        case kMC:
          if(!(ref = track->GetTrackRef(ily))) continue;
          fP = ref->P();
          break;
        case kRec:
          fP = p[ily];
          break;
        default: continue;
        } 
        Double32_t *it = &infoPID[ily*AliTRDCalPID::kNSlicesNN];
        for(Int_t is=AliTRDCalPID::kNSlicesNN; is--; it++) fdEdx[is] = (*it);
      }

      // momentum threshold
      if(fP < fPthreshold) continue;

      // store information
      fPIDdataArray->PushBack(ily, AliTRDpidUtil::GetMomentumBin(fP), fdEdx);
    }

    Fill();
  }

  PostData(1, fContainer);
  PostData(2, fData);
}


//________________________________________________________________________
void AliTRDpidRefMaker::Fill() 
{
// Fill data tree

  if(!fPIDdataArray->GetNtracklets()) return;
  // Fill data tree
  fData->Fill();

  
  // fill monitor
  for(Int_t itrklt=fPIDdataArray->GetNtracklets(); itrklt--;){
    Int_t pBin = fPIDdataArray->GetData(itrklt)->Momentum();
    ((TH2*)fContainer->At(0))->Fill(fPIDdataArray->GetPID(), pBin);
  }
}

//________________________________________________________________________
void AliTRDpidRefMaker::LinkPIDdata() 
{
// Link data tree to data members
  fData->SetBranchAddress("data", &fPIDdataArray);
}

//________________________________________________________________________
void AliTRDpidRefMaker::SetRefPID(ETRDpidRefMakerSource select, AliTRDtrackInfo *track, Float_t *pid) 
{
// Fill the reference PID values "pid" from "source" object
// according to the option "select". Possible options are
// - kV0  - v0 based PID
// - kMC  - MC truth [default]
// - kRec - outside detectors

  if(!track){
    AliError("No trackInfo found");
    return;
  }
  memset(fPID, 0, AliPID::kSPECIES*sizeof(Float_t));
  switch(select){ 
  case kV0:
    {
      //Get V0 PID decisions from the AliTRDv0Info for all particle species (implemented so far : electrons from conversions, pions from K0s and protons from Lambdas) :
      AliTRDv0Info *v0(NULL);
      for(Int_t iv(0); iv<fV0s->GetEntriesFast(); iv++){
        if(!(v0 = (AliTRDv0Info*)fV0s->At(iv))) continue;
        if(!v0->HasTrack(track)) continue;
        for(Int_t is=AliPID::kSPECIES; is--;) fPID[is] = v0->GetPID(is, track);
        break;
      }
    }
    break;
  case kMC:
    if(!HasMCdata()){
      AliError("Could not retrive reference PID from MC");
      return;
    }
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
    break;
  case kRec:
    { 
      AliTRDtrackV1 *trackTRD = track->GetTrack();
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
  AliDebug(4, Form("Ref PID [%] : %s[%5.2f] %s[%5.2f] %s[%5.2f] %s[%5.2f] %s[%5.2f]"
    ,AliPID::ParticleShortName(0), 1.e2*fPID[0]
    ,AliPID::ParticleShortName(1), 1.e2*fPID[1]
    ,AliPID::ParticleShortName(2), 1.e2*fPID[2]
    ,AliPID::ParticleShortName(3), 1.e2*fPID[3]
    ,AliPID::ParticleShortName(4), 1.e2*fPID[4]
  ));
}

//________________________________________________________________________
void AliTRDpidRefMaker::SetAbundance(Float_t train) 
{
// Split data sample between trainning and testing

  if(train<0. || train >1.){
    AliWarning("The input data should be in the interval [0, 1]");
    return;
  }

  fFreq = train; 
}

