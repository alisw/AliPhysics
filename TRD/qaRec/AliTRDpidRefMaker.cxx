#include "TPDGCode.h"
#include "TH1F.h"
#include "TTreeStream.h"

#include "AliPID.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliTrackReference.h"

#include "AliAnalysisTask.h"

#include "AliTRDtrackV1.h"
#include "AliTRDReconstructor.h"
#include "../Cal/AliTRDCalPID.h"
#include "../Cal/AliTRDCalPIDNN.h"

#include "AliTRDpidRefMaker.h"
#include "AliTRDtrackInfo/AliTRDtrackInfo.h"

// builds the reference tree for the training of neural networks


ClassImp(AliTRDpidRefMaker)

//________________________________________________________________________
AliTRDpidRefMaker::AliTRDpidRefMaker() 
  :AliTRDrecoTask("PIDR", "PID Reference Tree Maker")
  ,fReconstructor(0x0)
{
  //
  // Default constructor
  //

  fReconstructor = new AliTRDReconstructor();
  fReconstructor->SetRecoParam(AliTRDrecoParam::GetLowFluxParam());
}


//________________________________________________________________________
AliTRDpidRefMaker::~AliTRDpidRefMaker() 
{
  if(fReconstructor) delete fReconstructor;
}


//________________________________________________________________________
void AliTRDpidRefMaker::CreateOutputObjects()
{
  // Create histograms
  // Called once

  OpenFile(0, "RECREATE");
  fContainer = new TObjArray();

  fContainer->AddAt(new TH1F("hPDG","hPDG",AliPID::kSPECIES,-0.5,5.5),0);
}


//________________________________________________________________________
void AliTRDpidRefMaker::Exec(Option_t *) 
{
  // Main loop
  // Called for each event

  Int_t labelsacc[10000]; 
  memset(labelsacc, 0, sizeof(Int_t) * 10000);
  
  struct tChamb{
    Float_t Slide[AliTRDReconstructor::kNNslices];
  };
  tChamb Chamb[AliTRDCalPID::kNPlane];

  Float_t *fdEdx;

  Float_t mom;
  ULong_t status;
  Int_t nTRD = 0;
//   Float_t fdEdx[AliTRDCalPID::kNPlane][AliTRDReconstructor::kNNslices];
  Float_t v0pdg[AliPID::kSPECIES];

  AliTRDtrackInfo     *track = 0x0;
  AliTRDtrackV1    *TRDtrack = 0x0;
  AliTrackReference     *ref = 0x0;
  AliExternalTrackParam *esd = 0x0;

  AliTRDseedV1 *TRDtracklet[AliTRDCalPID::kNPlane];
  for(Int_t iChamb = 0; iChamb < AliTRDCalPID::kNPlane; iChamb++) TRDtracklet[iChamb] = 0x0;

  //AliTRDcluster *TRDcluster = 0x0;

  for(Int_t itrk=0; itrk<fTracks->GetEntriesFast(); itrk++){

    // reset the pid information
    for(Int_t iPart = 0; iPart < AliPID::kSPECIES; iPart++)
      v0pdg[iPart] = 0.;

    track = (AliTRDtrackInfo*)fTracks->UncheckedAt(itrk);
    if(!track->HasESDtrack()) continue;
    status = track->GetStatus();
    if(!(status&AliESDtrack::kTPCout)) continue;

    if(!(TRDtrack = track->GetTRDtrack())) continue; 
    //&&(track->GetNumberOfClustersRefit()

    // use only tracks that hit 6 chambers
    if(!(TRDtrack->GetNumberOfTracklets() == AliTRDCalPID::kNPlane)) continue;
     
    ref = track->GetTrackRef(0);
    esd = track->GetOuterParam();
    mom = ref ? ref->P(): esd->P();

    labelsacc[nTRD] = track->GetLabel();
    nTRD++;
      
    // if no monte carlo data available -> use V0 information
    if(!HasMCdata()){
      GetV0info(TRDtrack,v0pdg);
    }
    // else use the MC info
    else{
      switch(track -> GetPDG()){
      case kElectron:
      case kPositron:
        v0pdg[AliPID::kElectron] = 1.;
        break;
      case kMuonPlus:
      case kMuonMinus:
        v0pdg[AliPID::kMuon] = 1.;
        break;
      case kPiPlus:
      case kPiMinus:
        v0pdg[AliPID::kPion] = 1.;
        break;
      case kKPlus:
      case kKMinus:
        v0pdg[AliPID::kKaon] = 1.;
        break;
      case kProton:
      case kProtonBar:
        v0pdg[AliPID::kProton] = 1.;
        break;
      }
    }


    // set reconstructor
    TRDtrack -> SetReconstructor(fReconstructor);
    fReconstructor -> SetOption("nn");

    // fill the dE/dx information
    for(Int_t iChamb = 0; iChamb < AliTRDCalPID::kNPlane; iChamb++){
      TRDtracklet[iChamb] = TRDtrack -> GetTracklet(iChamb);
      fdEdx = TRDtracklet[iChamb] -> GetdEdx();
      for(Int_t iSlide = 0; iSlide < AliTRDReconstructor::kNNslices; iSlide++)
      Chamb[iChamb].Slide[iSlide] = fdEdx[iSlide]/AliTRDCalPIDNN::kMLPscale;
    }
    

    // fill the debug streams
    if(fDebugLevel>=2){
      (*fDebugStream) << "TreeInfo"
        << "isele="   << v0pdg[0]
        << "ismuo="   << v0pdg[1]
        << "ispio="   << v0pdg[2]
        << "iskao="   << v0pdg[3]
        << "ispro="   << v0pdg[4]
        << "\n";
      for(Int_t iChamb = 0; iChamb < AliTRDCalPID::kNPlane; iChamb++){
      (*fDebugStream) << Form("Chamb%d", iChamb)
          << "Slide0=" << Chamb[iChamb].Slide[0]
          << "Slide1=" << Chamb[iChamb].Slide[1]
          << "Slide2=" << Chamb[iChamb].Slide[2]
          << "Slide3=" << Chamb[iChamb].Slide[3]
          << "Slide4=" << Chamb[iChamb].Slide[4]
          << "Slide5=" << Chamb[iChamb].Slide[5]
          << "Slide6=" << Chamb[iChamb].Slide[6]
          << "Slide7=" << Chamb[iChamb].Slide[7]
          << "\n";
      }
    }

    for(Int_t iPart = 0; iPart < AliPID::kSPECIES; iPart++){
      if(fDebugLevel>=4) Printf("PDG is %d %f", iPart, v0pdg[iPart]);
    }
  }

  PostData(0, fContainer);
}


//________________________________________________________
void AliTRDpidRefMaker::GetRefFigure(Int_t /*ifig*/, Int_t &/*first*/, Int_t &/*last*/, Option_t */*opt*/)
{
  
}


//________________________________________________________________________
Bool_t AliTRDpidRefMaker::PostProcess()
{
  // Draw result to the screen
  // Called once at the end of the query

  return kTRUE; // testing protection
}


//________________________________________________________________________
void AliTRDpidRefMaker::Terminate(Option_t *) 
{
  // Draw result to the screen
  // Called once at the end of the query

  fContainer = dynamic_cast<TObjArray*>(GetOutputData(0));
  if (!fContainer) {
    Printf("ERROR: list not available");
    return;
  }
}


//________________________________________________________________________
void AliTRDpidRefMaker::GetV0info(AliTRDtrackV1 *TRDtrack, Float_t *v0pdg) 
{

  // !!!! PREMILMINARY FUNCTION !!!!
  //
  // this is the place for the V0 procedure
  // as long as there is no implemented, just the probabilities
  // of the TRDtrack is used!

  TRDtrack -> SetReconstructor(fReconstructor);
  fReconstructor -> SetOption("nn");
  TRDtrack -> CookPID();
  for(Int_t iPart = 0; iPart < AliPID::kSPECIES; iPart++){
    v0pdg[iPart] = TRDtrack -> GetPID(iPart);
    if(fDebugLevel>=4) Printf("PDG is (in V0info) %d %f", iPart, v0pdg[iPart]);
  }
}
