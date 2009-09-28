#include "TPDGCode.h"
#include "TFile.h"
#include "TTree.h"
#include "TEventList.h"
#include "TObjArray.h"
#include "TH1.h"

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
AliTRDpidRefMaker::AliTRDpidRefMaker(UChar_t n, const char *name, const char *title) 
  :AliTRDrecoTask(name, title)
  ,fReconstructor(0x0)
  ,fV0s(0x0)
  ,fData(0x0)
  ,fRefPID(kMC)
  ,fRefP(kMC)
  ,fNslices(n)
  ,fTrainFreq(1.)
  ,fTestFreq(0.)
  ,fScale(1.)
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

  memset(fTrain, 0, AliTRDCalPID::kNMom*AliTRDgeometry::kNlayer*sizeof(TEventList*));
  memset(fTest, 0, AliTRDCalPID::kNMom*AliTRDgeometry::kNlayer*sizeof(TEventList*));

  DefineInput(1, TObjArray::Class());
  DefineOutput(1, TTree::Class());
}


//________________________________________________________________________
AliTRDpidRefMaker::~AliTRDpidRefMaker() 
{
  if(fReconstructor) delete fReconstructor;
  for(Int_t ip=AliTRDCalPID::kNMom; ip--;){
    for(Int_t il=AliTRDgeometry::kNlayer; il--;){
      if(fTrain[ip][il]) delete fTrain[ip][il];
      if(fTest[ip][il]) delete fTest[ip][il];
    }
  }
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
  fContainer->AddAt(new TH1F("hPDG","hPDG",AliPID::kSPECIES,-0.5,5.5),0);

  // open reference TTree
  OpenFile(1, "RECREATE");
  fData = new TTree(GetName(), Form("Reference data for %s", GetName()));
  fData->Branch("l"  , &fLy  , "l/b");
  fData->Branch("p"  , &fP  , "p/F");
  fData->Branch("pid", fPID, Form("fPID[%d]/F", AliPID::kSPECIES));
  fData->Branch("dEdx" , fdEdx , Form("fdEdx[%d]/F", fNslices));
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
  AliTRDtrackV1    *TRDtrack = 0x0;
  AliTrackReference     *ref = 0x0;
  //AliExternalTrackParam *esd = 0x0;
  AliTRDseedV1 *TRDtracklet = 0x0;
  for(Int_t itrk=0, nTRD=0; itrk<fTracks->GetEntriesFast(); itrk++){
    // reset the pid information
    memset(fPID, 0, AliPID::kSPECIES*sizeof(Float_t));

    track = (AliTRDtrackInfo*)fTracks->UncheckedAt(itrk);
    if(!track->HasESDtrack()) continue;
    ULong_t status = track->GetStatus();
    if(!(status&AliESDtrack::kTPCout)) continue;

    if(!(TRDtrack = track->GetTrack())) continue; 
    //&&(track->GetNumberOfClustersRefit()

    // TOO STRONG and might introduce a bias if short 
    // tracks are to be analysed (A.Bercuci 23.09.09) 
    // use only tracks that hit 6 chambers
    //if(!(TRDtrack->GetNumberOfTracklets() == AliTRDgeometry::kNlayer)) continue;
     

    if(HasMCdata()) labelsacc[nTRD++] = track->GetLabel();
      

    switch(fRefPID){ 
    case kV0:
      SetRefPID(TRDtrack,fPID);
      break;
    case kMC:
      if(!HasMCdata()){
        AliError("Could not retrive reference PID from MC");
        return;
      }
      switch(track -> GetPDG()){
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
    default:
      AliWarning("PID reference source not implemented");
      return;
    }

    // set reconstructor
    Float_t *dedx;
    TRDtrack -> SetReconstructor(fReconstructor);

    // fill the momentum and dE/dx information
    for(Int_t ily = 0; ily < AliTRDgeometry::kNlayer; ily++){
      if(!(TRDtracklet = TRDtrack -> GetTracklet(ily))) continue;
      TRDtracklet->CookdEdx(fNslices);
      dedx = const_cast<Float_t*>(TRDtracklet->GetdEdx());
      for(Int_t is=fNslices; is--;) fdEdx[is] =  dedx[is]/fScale;

      switch(fRefP){
      case kMC:
        if(!HasMCdata()){
          AliError("Could not retrive reference momentum from MC");
          return;
        }
        if(!(ref = track->GetTrackRef(TRDtracklet))) continue;
        fP = ref->P();
        break;
      case kRec:
        fP = TRDtracklet->GetMomentum();
        break;
      default:
        AliWarning("Momentum reference source not implemented");
        return;
      }
      fLy = ily;
      fData->Fill();
    }

    for(Int_t iPart = 0; iPart < AliPID::kSPECIES; iPart++){
      AliDebug(4, Form("PDG is %d %f", iPart, fPID[iPart]));
    }
  }

  PostData(0, fContainer);
  PostData(1, fData);
}



//________________________________________________________________________
void AliTRDpidRefMaker::Terminate(Option_t *) 
{
  // Draw result to the screen
  // Called once at the end of the query

  fContainer = dynamic_cast<TObjArray*>(GetOutputData(0));
  if (!fContainer) {
    AliError("List not available");
    return;
  }
}



//________________________________________________________________________
void AliTRDpidRefMaker::MakeTrainingList() 
{
  //
  // build the training lists for the neural networks
  //

  if (!fData) LoadFile(Form("TRD.%s.root", GetName()));
  if (!fData) {
    AliError("Tree for training list not available");
    return;
  }

  Int_t nPart[AliTRDCalPID::kNMom][AliTRDgeometry::kNlayer][AliPID::kSPECIES];
  memset(nPart, 0, AliPID::kSPECIES*AliTRDCalPID::kNMom*AliTRDgeometry::kNlayer*sizeof(Int_t));

  // set needed branches
  fData -> SetBranchAddress("l", &fLy);
  fData -> SetBranchAddress("dEdx", &fdEdx);
  fData -> SetBranchAddress("pid", &fPID);
  fData -> SetBranchAddress("p", &fP);

  AliTRDpidUtil *util = new AliTRDpidUtil();

  // start first loop to check total number of each particle type
  Int_t n(0);
  for(Int_t iEv=0; iEv < fData -> GetEntries(); iEv++){
    fData -> GetEntry(iEv);
    
    // get particle type
    Int_t pidPos = TMath::LocMax(AliPID::kSPECIES, fPID);

    // get momentum bin
    Int_t ip = util->GetMomentumBin(fP);
    if(ip<0){
      AliWarning(Form("Wrong momentum %f.", fP));
      continue;
    } 
    nPart[ip][fLy][pidPos]++;
    n++;
  }

  AliDebug(2, Form("Particle multiplicities [%d]:", n));
  n=0;
  for(Int_t ip = 0; ip <AliTRDCalPID::kNMom; ip++){
    printf("  P[%2d] ", ip);
    for(Int_t is = 0; is <AliPID::kSPECIES; is++){  
      printf("%s[", AliPID::ParticleShortName(is));
      for(Int_t il = 0; il <AliTRDgeometry::kNlayer; il++){ 
        printf("%d ", nPart[ip][il][is]);
        n+=nPart[ip][il][is];
      }
      printf("] ");
    }
    printf("\n");
  }
  AliDebug(2, Form("Particle multiplicities check [%d]:", n));


  // set training/test sample size per momentum interval
  Int_t iTrain[AliTRDCalPID::kNMom], iTest[AliTRDCalPID::kNMom];
  for(Int_t ip = 0; ip < AliTRDCalPID::kNMom; ip++){
    Int_t min = 1000000;//TMath::MinElement(AliPID::kSPECIES, nPart[ip][5]);
    AliDebug(10, Form("Ref Stat in pBin[%d]=%d", ip, min));
    iTrain[ip] = Int_t(min * fTrainFreq);
    iTest[ip] = Int_t(min * fTestFreq);
    AliDebug(10, Form("P_bin[%2d]  Train[%d] Test[%d]", ip, iTrain[ip], iTest[ip]));
  }


  // start second loop to set the event lists
  // reset couters
  memset(nPart, 0, AliPID::kSPECIES*AliTRDCalPID::kNMom*AliTRDgeometry::kNlayer*sizeof(Int_t));
  for(Int_t iEv = 0; iEv < fData->GetEntries(); iEv++){
    fData->GetEntry(iEv);
    // get PID position
    Int_t pidPos = TMath::LocMax(AliPID::kSPECIES, fPID);


    Int_t iMomBin = util->GetMomentumBin(fP);
    if(nPart[iMomBin][fLy][pidPos] < iTrain[iMomBin]){
      fTrain[iMomBin][fLy]->Enter(iEv);
      nPart[iMomBin][fLy][pidPos]++;
    } else if(nPart[iMomBin][fLy][pidPos] < iTest[iMomBin]+iTrain[iMomBin]){
      fTest[iMomBin][fLy]->Enter(iEv);
      nPart[iMomBin][fLy][pidPos]++;
    } else continue;
  }
  
  AliDebug(2, "Particle multiplicities in both lists:");
  for(Int_t ip = 0; ip <AliTRDCalPID::kNMom; ip++){
    printf("  P[%2d] ", ip);
    for(Int_t is = 0; is <AliPID::kSPECIES; is++){  
      Int_t m(0);
      for(Int_t il = 0; il <AliTRDgeometry::kNlayer; il++) m+=nPart[ip][il][is];
      printf("%s[%4d] ", AliPID::ParticleShortName(is), m);
    }
    printf("\n");
  }

  util->Delete();
  //delete util;
}




//________________________________________________________________________
void AliTRDpidRefMaker::LoadFile(const Char_t *InFile) 
{
  //
  // Loads the files and sets the event list
  // for neural network training and 
  // building of the 2-dim reference histograms.
  // Usable for training outside of the makeResults.C macro
  //

  TFile::Open(InFile, "READ");
  fData = (TTree*)gFile->Get(GetName());

  for(Int_t iMom = 0; iMom < AliTRDCalPID::kNMom; iMom++){
    for(Int_t iLy = 0; iLy < AliTRDgeometry::kNlayer; iLy++){
      fTrain[iMom][iLy] = new TEventList(Form("fTrain_P%d_L%d", iMom, iLy), Form("Training list for momentum bin %d layer %d", iMom, iLy));
      fTest[iMom][iLy] = new TEventList(Form("fTest_P%d_L%d", iMom, iLy), Form("Test list for momentum bin %d layer %d", iMom, iLy));
    }
  }
}


//________________________________________________________________________
void AliTRDpidRefMaker::LoadContainer(const Char_t *InFileCont) 
{

  //
  // Loads the container if no container is there.
  // Useable for training outside of the makeResults.C macro
  //

  TFile::Open(InFileCont, "READ");
  fContainer = (TObjArray*)gFile->Get(GetName());
}


//________________________________________________________________________
void AliTRDpidRefMaker::SetRefPID(void *source, Float_t *pid) 
{
  // !!!! PREMILMINARY FUNCTION !!!!
  //
  // this is the place for the V0 procedure
  // as long as there is no one implemented, 
  // just the probabilities
  // of the TRDtrack are used!

  AliTRDtrackV1 *TRDtrack = static_cast<AliTRDtrackV1*>(source);
  TRDtrack -> SetReconstructor(fReconstructor);
  //fReconstructor -> SetOption("nn");
  TRDtrack -> CookPID();
  for(Int_t iPart = 0; iPart < AliPID::kSPECIES; iPart++){
    pid[iPart] = TRDtrack -> GetPID(iPart);
    AliDebug(4, Form("PDG is (in V0info) %d %f", iPart, pid[iPart]));
  }
}

//________________________________________________________________________
void AliTRDpidRefMaker::SetAbundance(Float_t train, Float_t test) 
{
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

