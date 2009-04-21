#include "TSystem.h"
#include "TDatime.h"
#include "TPDGCode.h"
#include "TH1F.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TTree.h"
#include "TTreeStream.h"
#include "TEventList.h"
#include "TMultiLayerPerceptron.h"

#include "AliPID.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliTrackReference.h"

#include "AliAnalysisTask.h"

#include "AliTRDtrackV1.h"
#include "AliTRDReconstructor.h"
#include "../Cal/AliTRDCalPID.h"
#include "../Cal/AliTRDCalPIDNN.h"

#include "AliTRDpidUtil.h"

#include "AliTRDpidRefMaker.h"
#include "info/AliTRDtrackInfo.h"

// builds the reference tree for the training of neural networks


ClassImp(AliTRDpidRefMaker)

//________________________________________________________________________
AliTRDpidRefMaker::AliTRDpidRefMaker() 
  :AliTRDrecoTask("PidRefMaker", "PID(NN) Reference Maker")
  ,fReconstructor(0x0)
  ,fNN(0x0)
  ,fLQ(0x0)
  ,fLayer(0xff)
  ,fTrainMomBin(kAll)
  ,fEpochs(1000)
  ,fMinTrain(100)
  ,fDate(0)
  ,fMom(0.)
  ,fDoTraining(0)
  ,fContinueTraining(0)
  ,fTrainPath(0x0)
{
  //
  // Default constructor
  //

  fReconstructor = new AliTRDReconstructor();
  fReconstructor->SetRecoParam(AliTRDrecoParam::GetLowFluxParam());
  memset(fv0pid, 0, AliPID::kSPECIES*sizeof(Float_t));
  memset(fdEdx, 0, 10*sizeof(Float_t));

  const Int_t nnSize = AliTRDCalPID::kNMom * AliTRDgeometry::kNlayer;
  memset(fTrain, 0, nnSize*sizeof(TEventList*));
  memset(fTest, 0, nnSize*sizeof(TEventList*));
  memset(fNet, 0, AliTRDgeometry::kNlayer*sizeof(TMultiLayerPerceptron*));

  TDatime datime;
  fDate = datime.GetDate();

  DefineOutput(1, TTree::Class());
  DefineOutput(2, TTree::Class());
}


//________________________________________________________________________
AliTRDpidRefMaker::~AliTRDpidRefMaker() 
{
  if(fReconstructor) delete fReconstructor;
  //if(fNN) delete fNN;
  //if(fLQ) delete fLQ;
}


//________________________________________________________________________
void AliTRDpidRefMaker::CreateOutputObjects()
{
  // Create histograms
  // Called once

  OpenFile(0, "RECREATE");
  fContainer = new TObjArray();
  fContainer->AddAt(new TH1F("hPDG","hPDG",AliPID::kSPECIES,-0.5,5.5),0);

  TGraphErrors *gEffisTrain = new TGraphErrors(kMoniTrain);
  gEffisTrain -> SetLineColor(4);
  gEffisTrain -> SetMarkerColor(4);
  gEffisTrain -> SetMarkerStyle(29);
  gEffisTrain -> SetMarkerSize(1);

  TGraphErrors *gEffisTest = new TGraphErrors(kMoniTrain);
  gEffisTest -> SetLineColor(2);
  gEffisTest -> SetMarkerColor(2);
  gEffisTest -> SetMarkerStyle(29);
  gEffisTest -> SetMarkerSize(1);

  fContainer -> AddAt(gEffisTrain,kGraphTrain);
  fContainer -> AddAt(gEffisTest,kGraphTest);

  // open reference TTree for NN
  OpenFile(1, "RECREATE");
  fNN = new TTree("NN", "Reference data for NN");
  fNN->Branch("fLayer", &fLayer, "fLayer/I");
  fNN->Branch("fMom", &fMom, "fMom/F");
  fNN->Branch("fv0pid", fv0pid, Form("fv0pid[%d]/F", AliPID::kSPECIES));
  fNN->Branch("fdEdx", fdEdx, Form("fdEdx[%d]/F", AliTRDpidUtil::kNNslices));

  // open reference TTree for LQ
  OpenFile(2, "RECREATE");
  fLQ = new TTree("LQ", "Reference data for LQ");
  fLQ->Branch("fLayer", &fLayer, "fLayer/I");
  fLQ->Branch("fMom", &fMom, "fMom/F");
  fLQ->Branch("fv0pid", fv0pid, Form("fv0pid[%d]/F", AliPID::kSPECIES));
  fLQ->Branch("fdEdx", fdEdx, Form("fdEdx[%d]/F", AliTRDpidUtil::kLQslices));
}


//________________________________________________________________________
void AliTRDpidRefMaker::Exec(Option_t *) 
{
  // Main loop
  // Called for each event

  Int_t labelsacc[10000]; 
  memset(labelsacc, 0, sizeof(Int_t) * 10000);
  
  Float_t mom;
  ULong_t status;
  Int_t nTRD = 0;

  AliTRDtrackInfo     *track = 0x0;
  AliTRDtrackV1    *TRDtrack = 0x0;
  AliTrackReference     *ref = 0x0;
  AliExternalTrackParam *esd = 0x0;

  AliTRDseedV1 *TRDtracklet = 0x0;

  //AliTRDcluster *TRDcluster = 0x0;

  for(Int_t itrk=0; itrk<fTracks->GetEntriesFast(); itrk++){

    // reset the pid information
    for(Int_t iPart = 0; iPart < AliPID::kSPECIES; iPart++)
      fv0pid[iPart] = 0.;

    track = (AliTRDtrackInfo*)fTracks->UncheckedAt(itrk);
    if(!track->HasESDtrack()) continue;
    status = track->GetStatus();
    if(!(status&AliESDtrack::kTPCout)) continue;

    if(!(TRDtrack = track->GetTrack())) continue; 
    //&&(track->GetNumberOfClustersRefit()

    // use only tracks that hit 6 chambers
    if(!(TRDtrack->GetNumberOfTracklets() == AliTRDgeometry::kNlayer)) continue;
     
    ref = track->GetTrackRef(0);
    esd = track->GetESDinfo()->GetOuterParam();
    mom = ref ? ref->P(): esd->P();
    fMom = mom;


    labelsacc[nTRD] = track->GetLabel();
    nTRD++;
      
    // if no monte carlo data available -> use V0 information
    if(!HasMCdata()){
      GetV0info(TRDtrack,fv0pid);
    }
    // else use the MC info
    else{
      switch(track -> GetPDG()){
      case kElectron:
      case kPositron:
        fv0pid[AliPID::kElectron] = 1.;
        break;
      case kMuonPlus:
      case kMuonMinus:
        fv0pid[AliPID::kMuon] = 1.;
        break;
      case kPiPlus:
      case kPiMinus:
        fv0pid[AliPID::kPion] = 1.;
        break;
      case kKPlus:
      case kKMinus:
        fv0pid[AliPID::kKaon] = 1.;
        break;
      case kProton:
      case kProtonBar:
        fv0pid[AliPID::kProton] = 1.;
        break;
      }
    }

    // set reconstructor
    Float_t *dedx;
    TRDtrack -> SetReconstructor(fReconstructor);

    // fill the dE/dx information for NN
    fReconstructor -> SetOption("nn");
    for(Int_t ily = 0; ily < AliTRDgeometry::kNlayer; ily++){
      if(!(TRDtracklet = TRDtrack -> GetTracklet(ily))) continue;
      TRDtracklet->CookdEdx(AliTRDpidUtil::kNNslices);
      dedx = TRDtracklet->GetdEdx();
      for(Int_t iSlice = 0; iSlice < AliTRDpidUtil::kNNslices; iSlice++)
	dedx[iSlice] = dedx[iSlice]/AliTRDCalPIDNN::kMLPscale;
      memcpy(fdEdx, dedx, AliTRDpidUtil::kNNslices*sizeof(Float_t));
      if(fDebugLevel>=2) Printf("LayerNN : %d", ily);
      fLayer = ily;
      fNN->Fill();
    }
    

    // fill the dE/dx information for LQ
    fReconstructor -> SetOption("!nn");
    for(Int_t ily = 0; ily < AliTRDgeometry::kNlayer; ily++){
      if(!(TRDtracklet = TRDtrack -> GetTracklet(ily))) continue;
      TRDtracklet->CookdEdx(AliTRDpidUtil::kLQslices);
      dedx = TRDtracklet->GetdEdx();
      memcpy(fdEdx, dedx, AliTRDpidUtil::kLQslices*sizeof(Float_t));
      if(fDebugLevel>=2) Printf("LayerLQ : %d", ily);
      fLayer = ily;
      fLQ->Fill();
    }
    

    for(Int_t iPart = 0; iPart < AliPID::kSPECIES; iPart++){
      if(fDebugLevel>=4) Printf("PDG is %d %f", iPart, fv0pid[iPart]);
    }
  }

  PostData(0, fContainer);
  PostData(1, fNN);
  PostData(2, fLQ);
}


//________________________________________________________________________
Bool_t AliTRDpidRefMaker::PostProcess()
{
  // Draw result to the screen
  // Called once at the end of the query

  // build the training andthe test list for the neural networks
  MakeTrainingLists();        
  if(!fDoTraining) return kTRUE;

  // train the neural networks and build the refrence histos for 2-dim LQ
  gSystem->Exec(Form("mkdir ./Networks_%d/",fDate));
  if(fDebugLevel>=2) Printf("TrainMomBin [%d] [%d]", fTrainMomBin, kAll);

  // train single network for a single momentum (recommended)
  if(!(fTrainMomBin == kAll)){
    if(fTrain[fTrainMomBin][0] -> GetN() < fMinTrain){
      if(fDebugLevel>=2) Printf("Warning in AliTRDpidRefMaker::PostProcess : Not enough events for training available! Please check Data sample!");
      return kFALSE;
    }
    TrainNetworks(fTrainMomBin);
    BuildLQRefs(fTrainMomBin);
    MonitorTraining(fTrainMomBin);
  }
  // train all momenta
  else{
    for(Int_t iMomBin = 0; iMomBin < AliTRDCalPID::kNMom; iMomBin++){
      if(fTrain[iMomBin][0] -> GetN() < fMinTrain){
	if(fDebugLevel>=2) Printf("Warning in AliTRDpidRefMaker::PostProcess : Not enough events for training available for momentum bin [%d]! Please check Data sample!", iMomBin);
	continue;
      }
      TrainNetworks(iMomBin);
      BuildLQRefs(fTrainMomBin);
      MonitorTraining(iMomBin);
    }
  }

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
void AliTRDpidRefMaker::GetV0info(AliTRDtrackV1 *TRDtrack, Float_t *v0pid) 
{
  // !!!! PREMILMINARY FUNCTION !!!!
  //
  // this is the place for the V0 procedure
  // as long as there is no one implemented, 
  // just the probabilities
  // of the TRDtrack are used!

  TRDtrack -> SetReconstructor(fReconstructor);
  fReconstructor -> SetOption("nn");
  TRDtrack -> CookPID();
  for(Int_t iPart = 0; iPart < AliPID::kSPECIES; iPart++){
    v0pid[iPart] = TRDtrack -> GetPID(iPart);
    if(fDebugLevel>=4) Printf("PDG is (in V0info) %d %f", iPart, v0pid[iPart]);
  }
}


//________________________________________________________________________
void AliTRDpidRefMaker::MakeTrainingLists() 
{
  //
  // build the training lists for the neural networks
  //

  if (!fNN) {
    LoadFiles("TRD.TaskPidRefMakerNN.root", "TRD.TaskPidRefMakerLQ.root");
  }

  if (!fNN) {
    Printf("ERROR tree for training list not available");
    return;
  }

  if(fDebugLevel>=2) Printf("\n Making training lists! \n");

  Int_t nPart[AliPID::kSPECIES][AliTRDCalPID::kNMom];
  memset(nPart, 0, AliPID::kSPECIES*AliTRDCalPID::kNMom*sizeof(Int_t));

  // set needed branches
  fNN -> SetBranchAddress("fv0pid", &fv0pid);
  fNN -> SetBranchAddress("fMom", &fMom);
  fNN -> SetBranchAddress("fLayer", &fLayer);

  AliTRDpidUtil *util = new AliTRDpidUtil();

  // start first loop to check total number of each particle type
  for(Int_t iEv=0; iEv < fNN -> GetEntries(); iEv++){
    fNN -> GetEntry(iEv);

    // use only events with goes through 6 layers TRD
    if(!fLayer == 0)
      continue;

    // set the 11 momentum bins
    Int_t iMomBin = -1;
    iMomBin = util -> GetMomentumBin(fMom);
    
    // check PID information and count particle types per momentum interval
    if(fv0pid[AliPID::kElectron] == 1)
      nPart[AliPID::kElectron][iMomBin]++;
    else if(fv0pid[AliPID::kMuon] == 1)
      nPart[AliPID::kMuon][iMomBin]++;
    else if(fv0pid[AliPID::kPion] == 1)
      nPart[AliPID::kPion][iMomBin]++;
    else if(fv0pid[AliPID::kKaon] == 1)
      nPart[AliPID::kKaon][iMomBin]++;
    else if(fv0pid[AliPID::kProton] == 1)
      nPart[AliPID::kProton][iMomBin]++;
  }

  if(fDebugLevel>=2){ 
    Printf("Particle multiplicities:");
    for(Int_t iMomBin = 0; iMomBin <AliTRDCalPID::kNMom; iMomBin++)
      Printf("Momentum[%d]  Elecs[%d] Muons[%d] Pions[%d] Kaons[%d] Protons[%d]", iMomBin, nPart[AliPID::kElectron][iMomBin], nPart[AliPID::kMuon][iMomBin], nPart[AliPID::kPion][iMomBin], nPart[AliPID::kKaon][iMomBin], nPart[AliPID::kProton][iMomBin]);
    Printf("\n");
  }

  // implement counter of training and test sample size
  Int_t iTrain[AliTRDCalPID::kNMom], iTest[AliTRDCalPID::kNMom];
  memset(iTrain, 0, AliTRDCalPID::kNMom*sizeof(Int_t));
  memset(iTest, 0, AliTRDCalPID::kNMom*sizeof(Int_t));

  // set training sample size per momentum interval to 2/3 
  // of smallest particle counter and test sample to 1/3
  for(Int_t iMomBin = 0; iMomBin < AliTRDCalPID::kNMom; iMomBin++){
    iTrain[iMomBin] = nPart[0][iMomBin];
    for(Int_t iPart = 1; iPart < AliPID::kSPECIES; iPart++){
      if(iTrain[iMomBin] > nPart[iPart][iMomBin])
	iTrain[iMomBin] = nPart[iPart][iMomBin];
    } 
    iTrain[iMomBin] = Int_t(iTrain[iMomBin] * .66);
    iTest[iMomBin] = Int_t( iTrain[iMomBin] * .5);
    if(fDebugLevel>=2) Printf("Momentum[%d]  Train[%d] Test[%d]", iMomBin, iTrain[iMomBin], iTest[iMomBin]);
  }
  if(fDebugLevel>=2) Printf("\n");


  // reset couters
  memset(nPart, 0, AliPID::kSPECIES*AliTRDCalPID::kNMom*sizeof(Int_t));

  // start second loop to set the event lists
  for(Int_t iEv = 0; iEv < fNN -> GetEntries(); iEv++){
    fNN -> GetEntry(iEv);

    // use only events with goes through 6 layers TRD
    if(!fLayer == 0)
      continue;

    // set the 11 momentum bins
    Int_t iMomBin = -1;
    iMomBin = util -> GetMomentumBin(fMom);
    
    // set electrons
    if(fv0pid[AliPID::kElectron] == 1){
      if(nPart[AliPID::kElectron][iMomBin] < iTrain[iMomBin]){
	for(Int_t ily = 0; ily < AliTRDgeometry::kNlayer; ily++)
	  fTrain[iMomBin][ily] -> Enter(iEv + ily);
	nPart[AliPID::kElectron][iMomBin]++;
      }
      else if(nPart[AliPID::kElectron][iMomBin] < iTest[iMomBin]+iTrain[iMomBin]){
	for(Int_t ily = 0; ily < AliTRDgeometry::kNlayer; ily++)
	  fTest[iMomBin][ily] -> Enter(iEv + ily);
	nPart[AliPID::kElectron][iMomBin]++;
      }
      else
	continue;
    }
    // set muons
    else if(fv0pid[AliPID::kMuon] == 1){
      if(nPart[AliPID::kMuon][iMomBin] < iTrain[iMomBin]){
	for(Int_t ily = 0; ily < AliTRDgeometry::kNlayer; ily++)
	  fTrain[iMomBin][ily] -> Enter(iEv + ily);
	nPart[AliPID::kMuon][iMomBin]++;
      }
      else if(nPart[AliPID::kMuon][iMomBin] < iTest[iMomBin]+iTrain[iMomBin]){
	for(Int_t ily = 0; ily < AliTRDgeometry::kNlayer; ily++)
	  fTest[iMomBin][ily] -> Enter(iEv + ily);
	nPart[AliPID::kMuon][iMomBin]++;
      }
      else
	continue;
    }
    // set pions
    else if(fv0pid[AliPID::kPion] == 1){
      if(nPart[AliPID::kPion][iMomBin] < iTrain[iMomBin]){
	for(Int_t ily = 0; ily < AliTRDgeometry::kNlayer; ily++)
	  fTrain[iMomBin][ily] -> Enter(iEv + ily);
	nPart[AliPID::kPion][iMomBin]++;
      }
      else if(nPart[AliPID::kPion][iMomBin] < iTest[iMomBin]+iTrain[iMomBin]){
	for(Int_t ily = 0; ily < AliTRDgeometry::kNlayer; ily++)
	  fTest[iMomBin][ily] -> Enter(iEv + ily);
	nPart[AliPID::kPion][iMomBin]++;
      }
      else
	continue;
    }
    // set kaons
    else if(fv0pid[AliPID::kKaon] == 1){
      if(nPart[AliPID::kKaon][iMomBin] < iTrain[iMomBin]){
	for(Int_t ily = 0; ily < AliTRDgeometry::kNlayer; ily++)
	  fTrain[iMomBin][ily] -> Enter(iEv + ily);
	nPart[AliPID::kKaon][iMomBin]++;
      }
      else if(nPart[AliPID::kKaon][iMomBin] < iTest[iMomBin]+iTrain[iMomBin]){
	for(Int_t ily = 0; ily < AliTRDgeometry::kNlayer; ily++)
	  fTest[iMomBin][ily] -> Enter(iEv + ily);
	nPart[AliPID::kKaon][iMomBin]++;
      }
      else
	continue;
    }
    // set protons
    else if(fv0pid[AliPID::kProton] == 1){
      if(nPart[AliPID::kProton][iMomBin] < iTrain[iMomBin]){
	for(Int_t ily = 0; ily < AliTRDgeometry::kNlayer; ily++)
	  fTrain[iMomBin][ily] -> Enter(iEv + ily);
	nPart[AliPID::kProton][iMomBin]++;
      }
      else if(nPart[AliPID::kProton][iMomBin] < iTest[iMomBin]+iTrain[iMomBin]){
	for(Int_t ily = 0; ily < AliTRDgeometry::kNlayer; ily++)
	  fTest[iMomBin][ily] -> Enter(iEv + ily);
	nPart[AliPID::kProton][iMomBin]++;
      }
      else
	continue;
    }
  }
  
  if(fDebugLevel>=2){ 
    Printf("Particle multiplicities in both lists:");
    for(Int_t iMomBin = 0; iMomBin <AliTRDCalPID::kNMom; iMomBin++)
      Printf("Momentum[%d]  Elecs[%d] Muons[%d] Pions[%d] Kaons[%d] Protons[%d]", iMomBin, nPart[AliPID::kElectron][iMomBin], nPart[AliPID::kMuon][iMomBin], nPart[AliPID::kPion][iMomBin], nPart[AliPID::kKaon][iMomBin], nPart[AliPID::kProton][iMomBin]);
    Printf("\n");
  }

  util -> Delete();
}


//________________________________________________________________________
void AliTRDpidRefMaker::TrainNetworks(Int_t mombin) 
{
  //
  // train the neural networks
  //
  
  
  if (!fNN) {
    LoadFiles("TRD.TaskPidRefMakerNN.root", "TRD.TaskPidRefMakerLQ.root");
  }

  if (!fNN) {
    Printf("ERROR tree for training list not available");
    return;
  }

  TDatime datime;
  fDate = datime.GetDate();

  if(fDebugLevel>=2) Printf("Training momentum bin %d", mombin);

  // set variable to monitor the training and to save the development of the networks
  Int_t nEpochs = fEpochs/kMoniTrain;       
  if(fDebugLevel>=2) Printf("Training %d times %d epochs", kMoniTrain, nEpochs);

  // make directories to save the networks 
  gSystem->Exec(Form("rm -r ./Networks_%d/MomBin_%d",fDate, mombin));
  gSystem->Exec(Form("mkdir ./Networks_%d/MomBin_%d",fDate, mombin));

  // variable to check if network can load weights from previous training
  Bool_t bFirstLoop[AliTRDgeometry::kNlayer];
  memset(bFirstLoop, kTRUE, AliTRDgeometry::kNlayer*sizeof(Bool_t));
 
  // train networks over several loops and save them after each loop
  for(Int_t iLoop = 0; iLoop < kMoniTrain; iLoop++){
    // loop over chambers
    for(Int_t iChamb = 0; iChamb < AliTRDgeometry::kNlayer; iChamb++){
      // set the event lists
      fNN -> SetEventList(fTrain[mombin][iChamb]);
      fNN -> SetEventList(fTest[mombin][iChamb]);
      
      if(fDebugLevel>=2) Printf("Trainingloop[%d] Chamber[%d]", iLoop, iChamb);
      
      // check if network is already implemented
      if(bFirstLoop[iChamb] == kTRUE){
	fNet[iChamb] = new TMultiLayerPerceptron("fdEdx[0],fdEdx[1],fdEdx[2],fdEdx[3],fdEdx[4],fdEdx[5],fdEdx[6],fdEdx[7]:15:7:fv0pid[0],fv0pid[1],fv0pid[2],fv0pid[3],fv0pid[4]!",fNN,fTrain[mombin][iChamb],fTest[mombin][iChamb]);
	fNet[iChamb] -> SetLearningMethod(TMultiLayerPerceptron::kStochastic);       // set learning method
	fNet[iChamb] -> TMultiLayerPerceptron::SetEta(0.001);                        // set learning speed
	if(!fContinueTraining){
	  if(fDebugLevel>=2) fNet[iChamb] -> Train(nEpochs,"text update=10, graph");
	  else fNet[iChamb] -> Train(nEpochs,"");
	}
	else{
	  fNet[iChamb] -> LoadWeights(Form("./Networks_%d/MomBin_%d/Net%d_%d",fTrainPath, mombin, iChamb, kMoniTrain - 1));
	  if(fDebugLevel>=2) fNet[iChamb] -> Train(nEpochs,"text update=10, graph+");      
	  else fNet[iChamb] -> Train(nEpochs,"+");                   
	}
	bFirstLoop[iChamb] = kFALSE;
      }
      else{    
	if(fDebugLevel>=2) fNet[iChamb] -> Train(nEpochs,"text update=10, graph+");      
	else fNet[iChamb] -> Train(nEpochs,"+");                   
      }
      
      // save weights for monitoring of the training
      fNet[iChamb] -> DumpWeights(Form("./Networks_%d/MomBin_%d/Net%d_%d",fDate, mombin, iChamb, iLoop));
    } // end chamber loop
  }   // end training loop
}


//________________________________________________________________________
void AliTRDpidRefMaker::BuildLQRefs(Int_t mombin) 
{
  //
  // build the 2-dim LQ reference histograms
  //

  if(fDebugLevel>=2) Printf("Building LQRefs for momentum bin %d", mombin);
}


//________________________________________________________________________
void AliTRDpidRefMaker::MonitorTraining(Int_t mombin) 
{
  //
  // train the neural networks
  //
  
  if(!fContainer){
    LoadContainer("TRD.TaskPidRefMaker.root");
  }
  if(!fContainer){
    Printf("ERROR container not available");
    return;
  }

  if (!fNN) {
    LoadFiles("TRD.TaskPidRefMakerNN.root", "TRD.TaskPidRefMakerLQ.root");
  }
  if (!fNN) {
    Printf("ERROR tree for training list not available");
    return;
  }

  // init networks and set event list
  for(Int_t iChamb = 0; iChamb < AliTRDgeometry::kNlayer; iChamb++){
	fNet[iChamb] = new TMultiLayerPerceptron("fdEdx[0],fdEdx[1],fdEdx[2],fdEdx[3],fdEdx[4],fdEdx[5],fdEdx[6],fdEdx[7]:15:7:fv0pid[0],fv0pid[1],fv0pid[2],fv0pid[3],fv0pid[4]!",fNN,fTrain[mombin][iChamb],fTest[mombin][iChamb]);   
	fNN -> SetEventList(fTrain[mombin][iChamb]);
	fNN -> SetEventList(fTest[mombin][iChamb]);
  }

  // implement variables for likelihoods
  Float_t Like[AliPID::kSPECIES][AliTRDgeometry::kNlayer];
  memset(Like, 0, AliPID::kSPECIES*AliTRDgeometry::kNlayer*sizeof(Float_t));
  Float_t LikeAll[AliPID::kSPECIES], TotProb;

  Double_t PionEffiTrain[kMoniTrain], PionEffiErrTrain[kMoniTrain];
  Double_t PionEffiTest[kMoniTrain], PionEffiErrTest[kMoniTrain];
  memset(PionEffiTrain, 0, kMoniTrain*sizeof(Double_t));
  memset(PionEffiErrTrain, 0, kMoniTrain*sizeof(Double_t));
  memset(PionEffiTest, 0, kMoniTrain*sizeof(Double_t));
  memset(PionEffiErrTest, 0, kMoniTrain*sizeof(Double_t));

  // init histos
  const Float_t epsilon = 1/(2*(AliTRDpidUtil::kBins-1));     // get nice histos with bin center at 0 and 1
  TH1F *hElecs, *hPions;
  hElecs = new TH1F("hElecs","Likelihood for electrons", AliTRDpidUtil::kBins, 0.-epsilon, 1.+epsilon);
  hPions = new TH1F("hPions","Likelihood for pions", AliTRDpidUtil::kBins, 0.-epsilon, 1.+epsilon);

  TGraphErrors *gEffisTrain=0x0, *gEffisTest=0x0;
  gEffisTrain = (TGraphErrors*)fContainer->At(kGraphTrain);
  gEffisTest = (TGraphErrors*)fContainer->At(kGraphTest);

  AliTRDpidUtil *util = new AliTRDpidUtil();
  
  // monitor training progress
  for(Int_t iLoop = 0; iLoop < kMoniTrain; iLoop++){

    // load weights
    for(Int_t iChamb = 0; iChamb < AliTRDgeometry::kNlayer; iChamb++){
      fNet[iChamb] -> LoadWeights(Form("./Networks_%d/MomBin_%d/Net%d_%d",fDate, mombin, iChamb, iLoop));
    }

    // event loop training list
    for(Int_t iEvent = 0; iEvent < fTrain[mombin][0] -> GetN(); iEvent++ ){

      // reset particle probabilities
      for(Int_t iPart = 0; iPart < AliPID::kSPECIES; iPart++){
	LikeAll[iPart] = 1./AliPID::kSPECIES;
      }
      TotProb = 0.;

      fNN -> GetEntry(fTrain[mombin][0] -> GetEntry(iEvent));
      // use event only if it is electron or pion
      if(!((fv0pid[AliPID::kElectron] == 1.0) || (fv0pid[AliPID::kPion] == 1.0))) continue;

      // get the probabilities for each particle type in each chamber
      for(Int_t iChamb = 0; iChamb < AliTRDgeometry::kNlayer; iChamb++){
	for(Int_t iPart = 0; iPart < AliPID::kSPECIES; iPart++){
	  Like[iPart][iChamb] = fNet[iChamb] -> Result(fTrain[mombin][iChamb] -> GetEntry(iEvent), iPart);
	  LikeAll[iPart] *=  Like[iPart][iChamb];
	}
      }

      // get total probability and normalize it
      for(Int_t iPart = 0; iPart < AliPID::kSPECIES; iPart++){
	TotProb += LikeAll[iPart];
      }
      for(Int_t iPart = 0; iPart < AliPID::kSPECIES; iPart++){
	LikeAll[iPart] /= TotProb;
      }

      // fill likelihood distributions
      if(fv0pid[AliPID::kElectron] == 1)      
	hElecs -> Fill(LikeAll[AliPID::kElectron]);
      if(fv0pid[AliPID::kPion] == 1)      
	hPions -> Fill(LikeAll[AliPID::kElectron]);
    } // end event loop


    // calculate the pion efficiency and fill the graph
    util -> CalculatePionEffi(hElecs, hPions);
    PionEffiTrain[iLoop] = util -> GetPionEfficiency();
    PionEffiErrTrain[iLoop] = util -> GetError();

    gEffisTrain -> SetPoint(iLoop, iLoop+1, PionEffiTrain[iLoop]);
    gEffisTrain -> SetPointError(iLoop, 0, PionEffiErrTrain[iLoop]);
    hElecs -> Reset();
    hPions -> Reset();
    if(fDebugLevel>=2) Printf("TrainingLoop[%d] PionEfficiency[%f +/- %f]", iLoop, PionEffiTrain[iLoop], PionEffiErrTrain[iLoop]);
    // end training loop
    


    // event loop test list
    for(Int_t iEvent = 0; iEvent < fTest[mombin][0] -> GetN(); iEvent++ ){

      // reset particle probabilities
      for(Int_t iPart = 0; iPart < AliPID::kSPECIES; iPart++){
	LikeAll[iPart] = 1./AliTRDgeometry::kNlayer;
      }
      TotProb = 0.;

      fNN -> GetEntry(fTest[mombin][0] -> GetEntry(iEvent));
      // use event only if it is electron or pion
      if(!((fv0pid[AliPID::kElectron] == 1.0) || (fv0pid[AliPID::kPion] == 1.0))) continue;

      // get the probabilities for each particle type in each chamber
      for(Int_t iChamb = 0; iChamb < AliTRDgeometry::kNlayer; iChamb++){
	for(Int_t iPart = 0; iPart < AliPID::kSPECIES; iPart++){
	  Like[iPart][iChamb] = fNet[iChamb] -> Result(fTest[mombin][iChamb] -> GetEntry(iEvent), iPart);
	  LikeAll[iPart] *=  Like[iPart][iChamb];
	}
      }

      // get total probability and normalize it
      for(Int_t iPart = 0; iPart < AliPID::kSPECIES; iPart++){
	TotProb += LikeAll[iPart];
      }
      for(Int_t iPart = 0; iPart < AliPID::kSPECIES; iPart++){
	LikeAll[iPart] /= TotProb;
      }

      // fill likelihood distributions
      if(fv0pid[AliPID::kElectron] == 1)      
	hElecs -> Fill(LikeAll[AliPID::kElectron]);
      if(fv0pid[AliPID::kPion] == 1)      
	hPions -> Fill(LikeAll[AliPID::kElectron]);
    } // end event loop

    // calculate the pion efficiency and fill the graph
    util -> CalculatePionEffi(hElecs, hPions);
    PionEffiTest[iLoop] = util -> GetPionEfficiency();
    PionEffiErrTest[iLoop] = util -> GetError();

    gEffisTest -> SetPoint(iLoop, iLoop+1, PionEffiTest[iLoop]);
    gEffisTest -> SetPointError(iLoop, 0, PionEffiErrTest[iLoop]);
    hElecs -> Reset();
    hPions -> Reset();
    if(fDebugLevel>=2) Printf("TestLoop[%d] PionEfficiency[%f +/- %f] \n", iLoop, PionEffiTest[iLoop], PionEffiErrTest[iLoop]);
    
  } //   end training loop
 
  util -> Delete();

  gEffisTest -> Draw("PAL");
  gEffisTrain -> Draw("PL");

}


//________________________________________________________________________
void AliTRDpidRefMaker::LoadFiles(const Char_t *InFileNN, const Char_t *InFileLQ) 
{
  //
  // Loads the files and sets the event list
  // for neural network training and 
  // building of the 2-dim reference histograms.
  // Useable for training outside of the makeResults.C macro
  //

  TFile *fInFileNN;
  fInFileNN = new TFile(InFileNN, "READ");
  fNN = (TTree*)fInFileNN -> Get("NN");

  TFile *fInFileLQ;
  fInFileLQ = new TFile(InFileLQ, "READ");
  fLQ = (TTree*)fInFileLQ -> Get("LQ");

  for(Int_t iMom = 0; iMom < AliTRDCalPID::kNMom; iMom++){
    for(Int_t ily = 0; ily < AliTRDgeometry::kNlayer; ily++){
      fTrain[iMom][ily] = new TEventList(Form("fTrainMom%d_%d", iMom, ily), Form("Training list for momentum intervall %d and plane %d", iMom, ily));
      fTest[iMom][ily] = new TEventList(Form("fTestMom%d_%d", iMom, ily), Form("Test list for momentum intervall %d and plane %d", iMom, ily));
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

  TFile *fInFileCont;
  fInFileCont = new TFile(InFileCont, "READ");
  fContainer = (TObjArray*)fInFileCont -> Get("PidRefMaker");

}


// //________________________________________________________________________
// void AliTRDpidRefMaker::CreateGraphs()
// {
//   // Create histograms
//   // Called once

//   OpenFile(0, "RECREATE");
//   fContainer = new TObjArray();
//   fContainer->AddAt(new TH1F("hPDG","hPDG",AliPID::kSPECIES,-0.5,5.5),0);

//   TGraphErrors *gEffisTrain = new TGraphErrors(kMoniTrain);
//   gEffisTrain -> SetLineColor(4);
//   gEffisTrain -> SetMarkerColor(4);
//   gEffisTrain -> SetMarkerStyle(29);
//   gEffisTrain -> SetMarkerSize(2);

//   TGraphErrors *gEffisTest = new TGraphErrors(kMoniTrain);
//   gEffisTest -> SetLineColor(2);
//   gEffisTest -> SetMarkerColor(2);
//   gEffisTest -> SetMarkerSize(2);

//   fContainer -> AddAt(gEffisTrain,kGraphTrain);
//   fContainer -> AddAt(gEffisTest,kGraphTest);
// }

