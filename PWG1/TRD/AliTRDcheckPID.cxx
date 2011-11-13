//////////////////////////////////////////////////////
//
// PID performance checker of the TRD
//
// Performs checks of ESD information for TRD-PID and recalculate PID based on OCDB information
// Also provides performance plots on detector based on the PID information - for the moment only 
// MC source is used but V0 is also possible. Here is a list of detector performance checks
//   - Integrated dE/dx per chamber
//   - <PH> as function of time bin and local radial position
//   - number of clusters/tracklet 
//   - number of tracklets/track 
//
// Author : Alex Wilk <wilka@uni-muenster.de>
//          Alex Bercuci <A.Bercuci@gsi.de>
//          Markus Fasel <M.Fasel@gsi.de>
//
///////////////////////////////////////////////////////

#include "TAxis.h"
#include "TROOT.h"
#include "TPDGCode.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TLegend.h"

#include <TClonesArray.h>
#include <TObjArray.h>
#include <TList.h>

#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliTrackReference.h"

#include "AliAnalysisTask.h"

#include "AliTRDtrackerV1.h"
#include "AliTRDtrackV1.h"
#include "AliTRDcluster.h"
#include "AliTRDReconstructor.h"
#include "AliCDBManager.h"
#include "AliTRDpidUtil.h"

#include "Cal/AliTRDCalPID.h"
#include "Cal/AliTRDCalPIDNN.h"
#include "AliTRDcheckPID.h"
#include "AliTRDinfoGen.h"
#include "AliAnalysisManager.h"
#include "info/AliTRDtrackInfo.h"
#include "info/AliTRDpidInfo.h"
#include "info/AliTRDv0Info.h"

Char_t const * AliTRDcheckPID::fgMethod[3] = {"LQ", "NN", "ESD"};

ClassImp(AliTRDcheckPID)

//________________________________________________________________________
AliTRDcheckPID::AliTRDcheckPID() 
  :AliTRDrecoTask()
  ,fUtil(NULL)
  ,fGraph(NULL)
  ,fPID(NULL)
  ,fV0s(NULL)
  ,fMomentumAxis(NULL)
  ,fMinNTracklets(AliTRDgeometry::kNlayer)
  ,fMaxNTracklets(AliTRDgeometry::kNlayer)
 {
  //
  // Default constructor
  //
  SetNameTitle("TRDcheckPID", "Check TRD PID");
  LocalInit();
}

//________________________________________________________________________
AliTRDcheckPID::AliTRDcheckPID(char* name ) 
  :AliTRDrecoTask(name, "Check TRD PID")
  ,fUtil(NULL)
  ,fGraph(NULL)
  ,fPID(NULL)
  ,fV0s(NULL)
  ,fMomentumAxis(NULL)
  ,fMinNTracklets(AliTRDgeometry::kNlayer)
  ,fMaxNTracklets(AliTRDgeometry::kNlayer)
 {
  //
  // Default constructor
  //

  LocalInit();
  InitFunctorList();

  DefineInput(3, TObjArray::Class());  // v0 list
  DefineOutput(2, TObjArray::Class()); // pid info list
}


//________________________________________________________________________
void AliTRDcheckPID::LocalInit() 
{
// Initialize working data

  // Initialize momentum axis with default values
  Double_t defaultMomenta[AliTRDCalPID::kNMom+1];
  for(Int_t imom = 0; imom < AliTRDCalPID::kNMom+1; imom++)
    defaultMomenta[imom] = AliTRDCalPID::GetMomentumBinning(imom);
  SetMomentumBinning(AliTRDCalPID::kNMom, defaultMomenta);

  memset(fEfficiency, 0, AliPID::kSPECIES*sizeof(TObjArray*));

  fUtil = new AliTRDpidUtil();
}

//________________________________________________________________________
AliTRDcheckPID::~AliTRDcheckPID() 
{
  if (AliAnalysisManager::GetAnalysisManager()->IsProofMode()) return;
  if(fPID){fPID->Delete(); delete fPID;}
  if(fGraph){fGraph->Delete(); delete fGraph;}
  if(fUtil) delete fUtil;
}


//________________________________________________________________________
void AliTRDcheckPID::UserCreateOutputObjects()
{
  // Create histograms
  // Called once
  
  AliTRDrecoTask::UserCreateOutputObjects();
  fPID = new TObjArray();
  fPID->SetOwner(kTRUE);
  PostData(2, fPID);
}

//________________________________________________________
void AliTRDcheckPID::UserExec(Option_t *opt)
{
  //
  // Execution part
  //

  fV0s = dynamic_cast<TObjArray *>(GetInputData(3));
  fPID->Delete();

  AliTRDrecoTask::UserExec(opt);
}


//_______________________________________________________
TObjArray * AliTRDcheckPID::Histos(){

  //
  // Create QA histograms
  //
  if(fContainer) return fContainer;

  Int_t xBins = AliPID::kSPECIES*fMomentumAxis->GetNbins(); 
  fContainer = new TObjArray(); fContainer->Expand(kNPlots); fContainer->SetOwner(kTRUE);

  const Float_t epsilon = 1./(2*(AliTRDpidUtil::kBins-1));     // get nice histos with bin center at 0 and 1
  TH1 *h = NULL;

  const Int_t kNmethodsPID=Int_t(sizeof(fgMethod)/sizeof(Char_t*));
  // histos with posterior probabilities for all particle species
  for(Int_t is=AliPID::kSPECIES; is--;){
    fEfficiency[is] = new TObjArray(kNmethodsPID);
    fEfficiency[is]->SetOwner();
    fEfficiency[is]->SetName(Form("Eff_%s", AliPID::ParticleShortName(is)));
    fContainer->AddAt(fEfficiency[is], is?kEfficiencyMu+is-1:kEfficiency);
    for(Int_t im=kNmethodsPID; im--;){
      if(!(h = (TH2F*)gROOT->FindObject(Form("PID_%s_%s", fgMethod[im], AliPID::ParticleShortName(is))))) {
        h = new TH2F(Form("PID_%s_%s", fgMethod[im], AliPID::ParticleShortName(is)), "", xBins, -0.5, xBins - 0.5,
          AliTRDpidUtil::kBins, 0.-epsilon, 1.+epsilon);
      } else h->Reset();
      fEfficiency[is]->AddAt(h, im);
    }
  }

  // histos of the dE/dx distribution for all 5 particle species and 11 momenta 
  if(!(h = (TH2F*)gROOT->FindObject("dEdx"))){
    h = new TH2F("dEdx", "", 
      xBins, -0.5, xBins - 0.5,
      200, 0, 15);
//       200, 0, 10000);
  } else h->Reset();
  fContainer->AddAt(h, kdEdx);

  // histos of the dE/dx slices for all 5 particle species and 11 momenta 
  if(!(h = (TH2F*)gROOT->FindObject("dEdxSlice"))){
    h = new TH2F("dEdxSlice", "", 
      xBins*AliTRDpidUtil::kLQslices, -0.5, xBins*AliTRDpidUtil::kLQslices - 0.5,
      200, 0, 5000);
  } else h->Reset();
  fContainer->AddAt(h, kdEdxSlice);

  // histos of the pulse height distribution for all 5 particle species and 11 momenta 
  TObjArray *fPH = new TObjArray(2);
  fPH->SetOwner(); fPH->SetName("PH");
  fContainer->AddAt(fPH, kPH);
  if(!(h = (TProfile2D*)gROOT->FindObject("PHT"))){
    h = new TProfile2D("PHT", "<PH>(tb);p*species;tb [100*ns];entries", 
      xBins, -0.5, xBins - 0.5,
      AliTRDseedV1::kNtb, -0.5, AliTRDseedV1::kNtb - 0.5);
  } else h->Reset();
  fPH->AddAt(h, 0);
  if(!(h = (TProfile2D*)gROOT->FindObject("PHX"))){
    h = new TProfile2D("PHX", "<PH>(x);p*species;x_{drift} [cm];entries", 
      xBins, -0.5, xBins - 0.5,
      40, 0., 4.5);
  } else h->Reset();
  fPH->AddAt(h, 1);

  // histos of the number of clusters distribution for all 5 particle species and 11 momenta 
  if(!(h = (TH2F*)gROOT->FindObject("NClus"))){
    h = new TH2F("NClus", "", 
      xBins, -0.5, xBins - 0.5,
      50, -0.5, 49.5);
  } else h->Reset();
  fContainer->AddAt(h, kNClus);


  // momentum distributions - absolute and in momentum bins
  if(!(h = (TH1F*)gROOT->FindObject("hMom"))){
    h = new TH1F("hMom", "momentum distribution", fMomentumAxis->GetNbins(), fMomentumAxis->GetXmin(), fMomentumAxis->GetXmax());
  } else h->Reset();
  fContainer->AddAt(h, kMomentum);
  
  if(!(h = (TH1F*)gROOT->FindObject("hMomBin"))){
    h = new TH1F("hMomBin", "momentum distribution in momentum bins", fMomentumAxis->GetNbins(), fMomentumAxis->GetXmin(), fMomentumAxis->GetXmax());
  } else h->Reset();
  fContainer->AddAt(h, kMomentumBin);

  // Number of tracklets per track
  if(!(h = (TH2F*)gROOT->FindObject("nTracklets"))){
    h = new TH2F("nTracklets", "", 
      xBins, -0.5, xBins - 0.5,
      AliTRDgeometry::kNlayer, 0.5, AliTRDgeometry::kNlayer+.5);
  } else h->Reset();
  fContainer->AddAt(h, kNTracklets);

  // V0 performance
  if(!(h = (TH1F*)gROOT->FindObject("nV0"))){
    h = new TH1F("nV0", "V0s/track;v0/track;entries", 
      6, -0.5, 5.5);
  } else h->Reset();
  fContainer->AddAt(h, kV0);

  // dQ/dl for 1D-Likelihood
  if(!(h = (TH1F *)gROOT->FindObject("dQdl"))){
    h = new TH2F("dQdl", "dQ/dl per layer;p*species;dQ/dl [a.u.]", xBins, -0.5, xBins - 0.5, 800, 0., 40000.);
  } else h->Reset();
  fContainer->AddAt(h, kdQdl);

  return fContainer;
}


//________________________________________________________________________
Bool_t AliTRDcheckPID::CheckTrackQuality(const AliTRDtrackV1* track) const
{
  //
  // Check if the track is ok for PID
  //
  
  Int_t ntracklets = track->GetNumberOfTracklets();
  if(ntracklets >= fMinNTracklets && ntracklets <= fMaxNTracklets) return 1;
//   if(!fESD)
//     return 0;

  return 0;
}

//________________________________________________________________________
Int_t AliTRDcheckPID::CalcPDG(AliTRDtrackV1* track) 
{
// Documentation to come

 /* track -> SetReconstructor(AliTRDinfoGen::Reconstructor());
  (const_cast<AliTRDrecoParam*>(AliTRDinfoGen::Reconstructor()->GetRecoParam()))->SetPIDNeuralNetwork();
  track -> CookPID();*/

  if(track -> GetPID(AliPID::kElectron) > track -> GetPID(AliPID::kMuon) + track -> GetPID(AliPID::kPion)  + track -> GetPID(AliPID::kKaon) + track -> GetPID(AliPID::kProton)){
    return kElectron;
  }
  else if(track -> GetPID(kProton) > track -> GetPID(AliPID::kPion)  && track -> GetPID(AliPID::kProton) > track -> GetPID(AliPID::kKaon)  && track -> GetPID(AliPID::kProton) > track -> GetPID(AliPID::kMuon)){
    return kProton;
  }
  else if(track -> GetPID(AliPID::kKaon) > track -> GetPID(AliPID::kMuon)  && track -> GetPID(AliPID::kKaon) > track -> GetPID(AliPID::kPion)){
    return kKPlus;
  }
  else if(track -> GetPID(AliPID::kMuon) > track -> GetPID(AliPID::kPion)){
    return kMuonPlus;
  }
  else{
    return kPiPlus;
  }
}


//_______________________________________________________
TH1 *AliTRDcheckPID::PlotLQ(const AliTRDtrackV1 *track)
{
  //
  // Plot the probabilities for electrons using 2-dim LQ
  //

  if(!fkESD){
    AliDebug(2, "No ESD info available.");
    return NULL;
  }

  if(track) fkTrack = track;
  if(!fkTrack){
    AliDebug(2, "No Track defined.");
    return NULL;
  }

  ULong_t status;
  status = fkESD -> GetStatus();
  if(!(status&AliESDtrack::kTRDin)) return NULL;

  if(!CheckTrackQuality(fkTrack)) return NULL;

  AliTRDtrackV1 cTrack(*fkTrack);
  cTrack.SetReconstructor(AliTRDinfoGen::Reconstructor());

  Int_t pdg = 0;
  Float_t momentum = 0.;

  if(fkMC){
    if(fkMC->GetTrackRef()) momentum = fkMC->GetTrackRef()->P();
    pdg = fkMC->GetPDG();
  } else{
    //AliWarning("No MC info available!");
    momentum = cTrack.GetMomentum(0);
    pdg = CalcPDG(&cTrack);
  }
  if(!IsInRange(momentum)) return NULL;

  (const_cast<AliTRDrecoParam*>(AliTRDinfoGen::Reconstructor()->GetRecoParam()))->SetPIDNeuralNetwork(kFALSE);
  cTrack.CookPID();
  if(cTrack.GetNumberOfTrackletsPID() < fMinNTracklets) return NULL;
  Int_t species = AliTRDpidUtil::Pdg2Pid(pdg);

  TH2F *hPIDLQ(NULL);
  TObjArray *eff(NULL);
  for(Int_t is=AliPID::kSPECIES; is--;){
    if(!(eff = dynamic_cast<TObjArray *>(fContainer->At(is?kEfficiencyMu+is-1:kEfficiency)))){
      AliWarning("No Histogram List defined.");
      return NULL;
    }
    if(!(hPIDLQ = dynamic_cast<TH2F *>(eff->At(AliTRDpidUtil::kLQ)))){
      AliWarning("No Histogram defined.");
      return NULL;
    }
  
    hPIDLQ -> Fill(FindBin(species, momentum), cTrack.GetPID(is));
  }
  return hPIDLQ;
}




//_______________________________________________________
TH1 *AliTRDcheckPID::PlotNN(const AliTRDtrackV1 *track)
{
  //
  // Plot the probabilities for electrons using 2-dim LQ
  //

  if(!fkESD){
    AliDebug(2, "No ESD info available.");
    return NULL;
  }

  if(track) fkTrack = track;
  if(!fkTrack){
    AliDebug(2, "No Track defined.");
    return NULL;
  }
  
  ULong_t status;
  status = fkESD -> GetStatus();
  if(!(status&AliESDtrack::kTRDin)) return NULL;

  if(!CheckTrackQuality(fkTrack)) return NULL;
  
  AliTRDtrackV1 cTrack(*fkTrack);
  cTrack.SetReconstructor(AliTRDinfoGen::Reconstructor());

  Int_t pdg = 0;
  Float_t momentum = 0.;
  if(fkMC){
    if(fkMC->GetTrackRef()) momentum = fkMC->GetTrackRef()->P();
    pdg = fkMC->GetPDG();
  } else {
    //AliWarning("No MC info available!");
    momentum = cTrack.GetMomentum(0);
    pdg = CalcPDG(&cTrack);
  }
  if(!IsInRange(momentum)) return NULL;

  (const_cast<AliTRDrecoParam*>(AliTRDinfoGen::Reconstructor()->GetRecoParam()))->SetPIDNeuralNetwork();
  cTrack.CookPID();
  if(cTrack.GetNumberOfTrackletsPID() < fMinNTracklets) return NULL;

  Int_t species = AliTRDpidUtil::Pdg2Pid(pdg);


  TH2F *hPIDNN(NULL);
  TObjArray *eff(NULL);
  for(Int_t is=AliPID::kSPECIES; is--;){
    if(!(eff = dynamic_cast<TObjArray *>(fContainer->At(is?kEfficiencyMu+is-1:kEfficiency)))){
      AliWarning("No Histogram List defined.");
      return NULL;
    }
    if(!(hPIDNN = dynamic_cast<TH2F *>(eff->At(AliTRDpidUtil::kNN)))){
      AliWarning("No Histogram defined.");
      return NULL;
    }
  
    hPIDNN->Fill(FindBin(species, momentum), cTrack.GetPID(is));
  }
  return hPIDNN;
}



//_______________________________________________________
TH1 *AliTRDcheckPID::PlotESD(const AliTRDtrackV1 *track)
{
  //
  // Plot the probabilities for electrons using 2-dim LQ
  //

  if(!fkESD){
    AliDebug(2, "No ESD info available.");
    return NULL;
  }

  if(track) fkTrack = track;
  if(!fkTrack){
    AliDebug(2, "No Track defined.");
    return NULL;
  }
  
  ULong_t status;
  status = fkESD -> GetStatus();
  if(!(status&AliESDtrack::kTRDin)) return NULL;

  if(!CheckTrackQuality(fkTrack)) return NULL;
  if(fkTrack->GetNumberOfTrackletsPID() < fMinNTracklets) return NULL;
  
  Int_t pdg = 0;
  Float_t momentum = 0.;
  if(fkMC){
    if(fkMC->GetTrackRef()) momentum = fkMC->GetTrackRef()->P();
    pdg = fkMC->GetPDG();
  } else {
    //AliWarning("No MC info available!");
    AliTRDtrackV1 cTrack(*fkTrack);
    cTrack.SetReconstructor(AliTRDinfoGen::Reconstructor());
    momentum = cTrack.GetMomentum(0);
    pdg = CalcPDG(&cTrack);
  }
  if(!IsInRange(momentum)) return NULL;

  const Double32_t *pidESD = fkESD->GetResponseIter();
  Int_t species = AliTRDpidUtil::Pdg2Pid(pdg);

  TH2F *hPID(NULL);
  TObjArray *eff(NULL);
  for(Int_t is=AliPID::kSPECIES; is--;){
    if(!(eff = dynamic_cast<TObjArray *>(fContainer->At(is?kEfficiencyMu+is-1:kEfficiency)))){
      AliWarning("No Histogram List defined.");
      return NULL;
    }
    if(!(hPID = dynamic_cast<TH2F *>(eff->At(AliTRDpidUtil::kESD)))){
      AliWarning("No Histogram defined.");
      return NULL;
    }

    hPID->Fill(FindBin(species, momentum), pidESD[is]);
  }
  return hPID;
}



//_______________________________________________________
TH1 *AliTRDcheckPID::PlotdQdl(const AliTRDtrackV1 *track){
  //
  // Plot the total charge for the 1D Likelihood method
  //
  if(track) fkTrack = track;
  if(!fkTrack){
    AliDebug(2, "No Track defined");
    return NULL;
  }
  TH2 *hdQdl = dynamic_cast<TH2F *>(fContainer->At(kdQdl));
  if(!hdQdl){
    AliWarning("No Histogram defined");
    return NULL;
  }

  if(!CheckTrackQuality(fkTrack)) return NULL;

  Int_t pdg = 0;
  Float_t momentum = 0.;
  AliTRDtrackV1 cTrack(*fkTrack);
  if(fkMC){
    if(fkMC->GetTrackRef()) momentum = fkMC->GetTrackRef()->P();
    pdg = fkMC->GetPDG();
  } else {
    //AliWarning("No MC info available!");
    momentum = cTrack.GetMomentum(0);
    pdg = CalcPDG(&cTrack);
  }
  if(!IsInRange(momentum)) return NULL;

  // Init exchange container
  Int_t s(AliTRDpidUtil::Pdg2Pid(pdg));
  Int_t ibin = FindBin(s, momentum);

  AliTRDseedV1 *tracklet = NULL;
  for(Int_t iseed = 0; iseed < 6; iseed++){
    if(!((tracklet = fkTrack->GetTracklet(iseed)) && tracklet->IsOK())) continue;
    hdQdl->Fill(ibin, tracklet->GetdQdl());
  }
  return hdQdl;
}

//_______________________________________________________
TH1 *AliTRDcheckPID::PlotdEdx(const AliTRDtrackV1 *track)
{
  //
  // Plot the probabilities for electrons using 2-dim LQ
  //

  if(track) fkTrack = track;
  if(!fkTrack){
    AliDebug(2, "No Track defined.");
    return NULL;
  }
  
  if(!CheckTrackQuality(fkTrack)) return NULL;
  
  TH2F *hdEdx;
  if(!(hdEdx = dynamic_cast<TH2F *>(fContainer->At(kdEdx)))){
    AliWarning("No Histogram defined.");
    return NULL;
  }

  AliTRDtrackV1 cTrack(*fkTrack);
  cTrack.SetReconstructor(AliTRDinfoGen::Reconstructor());
  Int_t pdg = 0;
  Float_t momentum = 0.;
  if(fkMC){
    if(fkMC->GetTrackRef()) momentum = fkMC->GetTrackRef()->P();
    pdg = fkMC->GetPDG();
  } else {
    //AliWarning("No MC info available!");
    momentum = cTrack.GetMomentum(0);
    pdg = CalcPDG(&cTrack);
  }
  if(!IsInRange(momentum)) return NULL;

  // Init exchange container
  Int_t s(AliTRDpidUtil::Pdg2Pid(pdg));
  AliTRDpidInfo *pid = new AliTRDpidInfo(s);

  (const_cast<AliTRDrecoParam*>(AliTRDinfoGen::Reconstructor()->GetRecoParam()))->SetPIDNeuralNetwork(kTRUE);

  Float_t sumdEdx(0.);
  Int_t iBin = FindBin(s, momentum);
  AliTRDseedV1 *tracklet = NULL;
  for(Int_t ily = 0; ily < AliTRDgeometry::kNlayer; ily++){
    tracklet = cTrack.GetTracklet(ily);
    if(!tracklet) continue;
    tracklet -> CookdEdx(AliTRDpidUtil::kNNslices);

    // fill exchange container
    pid->PushBack(tracklet->GetPlane(), 
                  AliTRDpidUtil::GetMomentumBin(tracklet->GetMomentum()), tracklet->GetdEdx());

    sumdEdx = 0.;
    for(Int_t i = AliTRDpidUtil::kNNslices; i--;) sumdEdx += tracklet->GetdEdx()[i];
    sumdEdx /= AliTRDCalPIDNN::kMLPscale;
    hdEdx -> Fill(iBin, sumdEdx);
  }
  fPID->Add(pid);

  return hdEdx;
}


//_______________________________________________________
TH1 *AliTRDcheckPID::PlotdEdxSlice(const AliTRDtrackV1 *track)
{
  //
  // Plot the probabilities for electrons using 2-dim LQ
  //

  if(track) fkTrack = track;
  if(!fkTrack){
    AliDebug(2, "No Track defined.");
    return NULL;
  }
  
  if(!CheckTrackQuality(fkTrack)) return NULL;
  
  TH2F *hdEdxSlice;
  if(!(hdEdxSlice = dynamic_cast<TH2F *>(fContainer->At(kdEdxSlice)))){
    AliWarning("No Histogram defined.");
    return NULL;
  }

  AliTRDtrackV1 cTrack(*fkTrack);
  cTrack.SetReconstructor(AliTRDinfoGen::Reconstructor());
  Int_t pdg = 0;
  Float_t momentum = 0.;
  if(fkMC){
    if(fkMC->GetTrackRef()) momentum = fkMC->GetTrackRef()->P();
    pdg = fkMC->GetPDG();
  } else {
    //AliWarning("No MC info available!");
    momentum = cTrack.GetMomentum(0);
    pdg = CalcPDG(&cTrack);
  }
  if(!IsInRange(momentum)) return NULL;

  (const_cast<AliTRDrecoParam*>(AliTRDinfoGen::Reconstructor()->GetRecoParam()))->SetPIDNeuralNetwork(kFALSE);
  Int_t iMomBin = fMomentumAxis->FindBin(momentum);
  Int_t species = AliTRDpidUtil::Pdg2Pid(pdg);
  Float_t *fdEdx;
  AliTRDseedV1 *tracklet = NULL;
  for(Int_t iChamb = 0; iChamb < AliTRDgeometry::kNlayer; iChamb++){
    tracklet = cTrack.GetTracklet(iChamb);
    if(!tracklet) continue;
    tracklet -> CookdEdx(AliTRDpidUtil::kLQslices);
    fdEdx = const_cast<Float_t *>(tracklet->GetdEdx());
    for(Int_t iSlice = 0; iSlice < AliTRDpidUtil::kLQslices; iSlice++){
      hdEdxSlice -> Fill(species * fMomentumAxis->GetNbins() * AliTRDpidUtil::kLQslices + (iMomBin-1) * AliTRDpidUtil::kLQslices + iSlice, fdEdx[iSlice]);
    }
  }  

  return hdEdxSlice;
}


//_______________________________________________________
TH1 *AliTRDcheckPID::PlotPH(const AliTRDtrackV1 *track)
{
  //
  // Plot the probabilities for electrons using 2-dim LQ
  //

  if(track) fkTrack = track;
  if(!fkTrack){
    AliDebug(2, "No Track defined.");
    return NULL;
  }
  
  if(!CheckTrackQuality(fkTrack)) return NULL;
  
  TObjArray *arr = NULL;
  TProfile2D *hPHX, *hPHT;
  if(!(arr = dynamic_cast<TObjArray *>(fContainer->At(kPH)))){
    AliWarning("No Histogram defined.");
    return NULL;
  }
  hPHT = (TProfile2D*)arr->At(0);
  hPHX = (TProfile2D*)arr->At(1);

  Int_t pdg = 0;
  Float_t momentum = 0.;
  if(fkMC){
    if(fkMC->GetTrackRef()) momentum = fkMC->GetTrackRef()->P();
    pdg = fkMC->GetPDG();
  } else {
    //AliWarning("No MC info available!");
    AliTRDtrackV1 cTrack(*fkTrack);
    cTrack.SetReconstructor(AliTRDinfoGen::Reconstructor());
    momentum = cTrack.GetMomentum(0);
    pdg = CalcPDG(&cTrack);
  }
  if(!IsInRange(momentum)) return NULL;;

  AliTRDseedV1 *tracklet = NULL;
  AliTRDcluster *cluster = NULL;
  Int_t species = AliTRDpidUtil::Pdg2Pid(pdg);
  Int_t iBin = FindBin(species, momentum);
  for(Int_t iChamb = 0; iChamb < AliTRDgeometry::kNlayer; iChamb++){
    tracklet = fkTrack->GetTracklet(iChamb);
    if(!tracklet) continue;
    Float_t x0 = tracklet->GetX0(); 
    for(Int_t ic = 0; ic < AliTRDseedV1::kNclusters; ic++){
      if(!(cluster = tracklet->GetClusters(ic))) continue;
      hPHT -> Fill(iBin, cluster->GetLocalTimeBin(), TMath::Abs(cluster->GetQ()));
      if(ic<AliTRDseedV1::kNtb) hPHX -> Fill(iBin, x0 - cluster->GetX(), tracklet->GetdQdl(ic));
    }
  }
  return hPHT;
}


//_______________________________________________________
TH1 *AliTRDcheckPID::PlotNClus(const AliTRDtrackV1 *track)
{
  //
  // Plot the probabilities for electrons using 2-dim LQ
  //

  if(track) fkTrack = track;
  if(!fkTrack){
    AliDebug(2, "No Track defined.");
    return NULL;
  }
  
  if(!CheckTrackQuality(fkTrack)) return NULL;
  
  TH2F *hNClus;
  if(!(hNClus = dynamic_cast<TH2F *>(fContainer->At(kNClus)))){
    AliWarning("No Histogram defined.");
    return NULL;
  }


  Int_t pdg = 0;
  Float_t momentum = 0.;
  if(fkMC){
    if(fkMC->GetTrackRef()) momentum = fkMC->GetTrackRef()->P();
    pdg = fkMC->GetPDG();
  } else {
    //AliWarning("No MC info available!");
    AliTRDtrackV1 cTrack(*fkTrack);
    cTrack.SetReconstructor(AliTRDinfoGen::Reconstructor());
    momentum = cTrack.GetMomentum(0);
    pdg = CalcPDG(&cTrack);
  }
  if(!IsInRange(momentum)) return NULL;

  Int_t species = AliTRDpidUtil::AliTRDpidUtil::Pdg2Pid(pdg);
  Int_t iBin = FindBin(species, momentum); 
  AliTRDseedV1 *tracklet = NULL;
  for(Int_t iChamb = 0; iChamb < AliTRDgeometry::kNlayer; iChamb++){
    tracklet = fkTrack->GetTracklet(iChamb);
    if(!tracklet) continue;
    hNClus -> Fill(iBin, tracklet->GetN());
  }

  return hNClus;
}

//_______________________________________________________
TH1 *AliTRDcheckPID::PlotNTracklets(const AliTRDtrackV1 *track)
{
  //
  // Plot the probabilities for electrons using 2-dim LQ
  //

  if(track) fkTrack = track;
  if(!fkTrack){
    AliDebug(2, "No Track defined.");
    return NULL;
  }
  
  TH2F *hTracklets;
  if(!(hTracklets = dynamic_cast<TH2F *>(fContainer->At(kNTracklets)))){
    AliWarning("No Histogram defined.");
    return NULL;
  }

  AliTRDtrackV1 cTrack(*fkTrack);
  cTrack.SetReconstructor(AliTRDinfoGen::Reconstructor());
  Int_t pdg = 0;
  Float_t momentum = 0.;
  if(fkMC){
    if(fkMC->GetTrackRef()) momentum = fkMC->GetTrackRef()->P();
    pdg = fkMC->GetPDG();
  } else {
    //AliWarning("No MC info available!");
    momentum = cTrack.GetMomentum();
    pdg = CalcPDG(&cTrack);
  }
  Int_t species = AliTRDpidUtil::Pdg2Pid(pdg);
  if(!IsInRange(momentum)) return NULL;

  Int_t iBin = FindBin(species, momentum);
  hTracklets -> Fill(iBin, cTrack.GetNumberOfTracklets());
  return hTracklets;
}

//_______________________________________________________
TH1 *AliTRDcheckPID::PlotMom(const AliTRDtrackV1 *track)
{
  //
  // Plot the probabilities for electrons using 2-dim LQ
  //

  if(track) fkTrack = track;
  if(!fkTrack){
    AliDebug(2, "No Track defined.");
    return NULL;
  }
  
  if(!CheckTrackQuality(fkTrack)) return NULL;
  
  TH1F *hMom;
  if(!(hMom = dynamic_cast<TH1F *>(fContainer->At(kMomentum)))){
    AliWarning("No Histogram defined.");
    return NULL;
  }


  Int_t pdg = 0;
  Float_t momentum = 0.;
  if(fkMC){
    if(fkMC->GetTrackRef()) momentum = fkMC->GetTrackRef()->P();
    pdg = fkMC->GetPDG();
  } else {
    //AliWarning("No MC info available!");
    AliTRDtrackV1 cTrack(*fkTrack);
    cTrack.SetReconstructor(AliTRDinfoGen::Reconstructor());
    momentum = cTrack.GetMomentum(0);
    pdg = CalcPDG(&cTrack);
  }
  if(IsInRange(momentum)) hMom -> Fill(momentum);
  return hMom;
}


//_______________________________________________________
TH1 *AliTRDcheckPID::PlotMomBin(const AliTRDtrackV1 *track)
{
  //
  // Plot the probabilities for electrons using 2-dim LQ
  //

  if(track) fkTrack = track;
  if(!fkTrack){
    AliDebug(2, "No Track defined.");
    return NULL;
  }
  
  if(!CheckTrackQuality(fkTrack)) return NULL;
  
  TH1F *hMomBin;
  if(!(hMomBin = dynamic_cast<TH1F *>(fContainer->At(kMomentumBin)))){
    AliWarning("No Histogram defined.");
    return NULL;
  }


  Int_t pdg = 0;
  Float_t momentum = 0.;

  if(fkMC){
    if(fkMC->GetTrackRef()) momentum = fkMC->GetTrackRef()->P();
    pdg = fkMC->GetPDG();
  } else {
    //AliWarning("No MC info available!");
    AliTRDtrackV1 cTrack(*fkTrack);
    cTrack.SetReconstructor(AliTRDinfoGen::Reconstructor());
    momentum = cTrack.GetMomentum(0);
  }
  if(IsInRange(momentum)) hMomBin -> Fill(fMomentumAxis->FindBin(momentum));
  return hMomBin;
}

//_______________________________________________________
TH1 *AliTRDcheckPID::PlotV0(const AliTRDtrackV1 *track)
{
  //
  // Plot the V0 performance against MC
  //

  if(track) fkTrack = track;
  if(!fkTrack){
    AliDebug(2, "No Track defined.");
    return NULL;
  }  
  if(!fkESD->HasV0()) return NULL;
  if(!HasMCdata()){ 
    AliDebug(1, "No MC defined.");
    return NULL;
  }
  if(!fContainer){
    AliWarning("No output container defined.");
    return NULL;
  }
  AliDebug(2, Form("TRACK[%d] species[%s][%d]\n", fkESD->GetId(), fkMC->GetPID()>=0?AliPID::ParticleShortName(fkMC->GetPID()):"none", fkMC->GetPDG()));

  TH1 *h(NULL);
  if(!(h = dynamic_cast<TH1F*>(fContainer->At(kV0)))) return NULL;
  Int_t sgn(0), n(0); AliTRDv0Info *v0(NULL);
  for(Int_t iv0(fV0s->GetEntriesFast()); iv0--;){
    if(!(v0=(AliTRDv0Info*)fV0s->At(iv0))) continue;
    if(!(sgn = v0->HasTrack(fkESD->GetId()))) continue;
    //for(Int_t is=AliPID::kSPECIES; is--;) v0->GetPID(is, track);
    //v0->Print();
    n++;
    //break;
  }
  h->Fill(n);
  return h;
}

//________________________________________________________
Bool_t AliTRDcheckPID::GetRefFigure(Int_t ifig)
{
// Steering function to retrieve performance plots

  Bool_t kFIRST = kTRUE;
  TGraphErrors *g = NULL;
  TAxis *ax = NULL;
  TObjArray *arr = NULL;
  TH1 *h1 = NULL, *h=NULL;
  TH2 *h2 = NULL;
  TList *content = NULL;
  switch(ifig){
  case kEfficiency:{
    gPad->Divide(2, 1, 1.e-5, 1.e-5);
    TList *l=gPad->GetListOfPrimitives();
    TVirtualPad *pad = ((TVirtualPad*)l->At(0));pad->cd();
    pad->SetMargin(0.1, 0.01, 0.1, 0.01);

    TLegend *legEff = new TLegend(.64, .84, .98, .98);
    legEff->SetBorderSize(1);legEff->SetTextSize(0.03255879);
    legEff->SetFillColor(0);
    h=new TH1S("hEff", "", 1, .5, 11.);
    h->SetLineColor(1);h->SetLineWidth(1);
    ax = h->GetXaxis();
    ax->SetTitle("p [GeV/c]");
    ax->SetRangeUser(.5, 11.);
    ax->SetMoreLogLabels();
    ax = h->GetYaxis();
    ax->SetTitle("#epsilon_{#pi} [%]");
    ax->CenterTitle();
    ax->SetRangeUser(1.e-2, 10.);
    h->Draw();
    content = (TList *)fGraph->FindObject(Form("Eff_%s", AliTRDCalPID::GetPartName(AliPID::kPion)));
    if(!(g = (TGraphErrors*)content->At(AliTRDpidUtil::kLQ))) break;
    if(!g->GetN()) break;
    legEff->SetHeader("PID Method [PION]");
    g->Draw("pc"); legEff->AddEntry(g, "LQ 2D", "pl");
    if(! (g = (TGraphErrors*)content->At(AliTRDpidUtil::kNN))) break;
    g->Draw("pc"); legEff->AddEntry(g, "NN", "pl");
    if(! (g = (TGraphErrors*)content->At(AliTRDpidUtil::kESD))) break;
    g->Draw("p"); legEff->AddEntry(g, "ESD", "pl");
    legEff->Draw();
    gPad->SetLogy();
    gPad->SetLogx();
    gPad->SetGridy();
    gPad->SetGridx();


    pad = ((TVirtualPad*)l->At(1));pad->cd();
    pad->SetMargin(0.1, 0.01, 0.1, 0.01);
    h=new TH1S("hThr", "", 1, .5, 11.);
    h->SetLineColor(1);h->SetLineWidth(1);
    ax = h->GetXaxis();
    ax->SetTitle("p [GeV/c]");
    ax->SetMoreLogLabels();
    ax = h->GetYaxis();
    ax->SetTitle("Threshold [%]");
    ax->SetRangeUser(5.e-2, 1.);
    h->Draw();
    content = (TList *)fGraph->FindObject("Thres");
    if(!(g = (TGraphErrors*)content->At(AliTRDpidUtil::kLQ))) break;
    if(!g->GetN()) break;
    g->Draw("pc");
    if(!(g = (TGraphErrors*)content->At(AliTRDpidUtil::kNN))) break;
    g->Draw("pc");
    if(!(g = (TGraphErrors*)content->At(AliTRDpidUtil::kESD))) break;
    g->Draw("p");
    gPad->SetLogx();
    gPad->SetGridy();
    gPad->SetGridx();
    return kTRUE;
  }
  case kEfficiencyKa:{
    gPad->Divide(2, 1, 1.e-5, 1.e-5);
    TList *l=gPad->GetListOfPrimitives();
    TVirtualPad *pad = ((TVirtualPad*)l->At(0));pad->cd();
    pad->SetMargin(0.1, 0.01, 0.1, 0.01);

    TLegend *legEff = new TLegend(.64, .84, .98, .98);
    legEff->SetBorderSize(1);legEff->SetTextSize(0.03255879);
    legEff->SetFillColor(0);
    h = (TH1S*)gROOT->FindObject("hEff");
    h=(TH1S*)h->Clone("hEff_K");
    h->SetYTitle("#epsilon_{K} [%]");
    h->GetYaxis()->SetRangeUser(1.e-2, 1.e2);
    h->Draw();
    content = (TList *)fGraph->FindObject(Form("Eff_%s", AliTRDCalPID::GetPartName(AliPID::kKaon)));
    if(!(g = (TGraphErrors*)content->At(AliTRDpidUtil::kLQ))) break;
    if(!g->GetN()) break;
    legEff->SetHeader("PID Method [KAON]");
    g->Draw("pc"); legEff->AddEntry(g, "LQ 2D", "pl");
    if(! (g = (TGraphErrors*)content->At(AliTRDpidUtil::kNN))) break;
    g->Draw("pc"); legEff->AddEntry(g, "NN", "pl");
    if(! (g = (TGraphErrors*)content->At(AliTRDpidUtil::kESD))) break;
    g->Draw("p"); legEff->AddEntry(g, "ESD", "pl");
    legEff->Draw();
    gPad->SetLogy();
    gPad->SetLogx();
    gPad->SetGridy();
    gPad->SetGridx();

    TLegend *legEff2 = new TLegend(.64, .84, .98, .98);
    legEff2->SetBorderSize(1);legEff2->SetTextSize(0.03255879);
    legEff2->SetFillColor(0);
    pad = ((TVirtualPad*)l->At(1));pad->cd();
    pad->SetMargin(0.1, 0.01, 0.1, 0.01);
    h=(TH1S*)h->Clone("hEff_p");
    h->SetYTitle("#epsilon_{p} [%]");
    h->GetYaxis()->SetRangeUser(1.e-2, 1.e2);
    h->Draw();
    content = (TList *)fGraph->FindObject(Form("Eff_%s", AliTRDCalPID::GetPartName(AliPID::kProton)));
    if(!(g = (TGraphErrors*)content->At(AliTRDpidUtil::kLQ))) break;
    if(!g->GetN()) break;
    legEff2->SetHeader("PID Method [PROTON]");
    g->Draw("pc"); legEff2->AddEntry(g, "LQ 2D", "pl");
    if(! (g = (TGraphErrors*)content->At(AliTRDpidUtil::kNN))) break;
    g->Draw("pc"); legEff2->AddEntry(g, "NN", "pl");
    if(! (g = (TGraphErrors*)content->At(AliTRDpidUtil::kESD))) break;
    g->Draw("p"); legEff2->AddEntry(g, "ESD", "pl");
    legEff2->Draw();
    gPad->SetLogy();
    gPad->SetLogx();
    gPad->SetGridy();
    gPad->SetGridx();
    return kTRUE;
  }
  case kdEdx:{
    // save 2.0 GeV projection as reference
    TLegend *legdEdx = new TLegend(.7, .7, .98, .98);
    legdEdx->SetBorderSize(1);
    kFIRST = kTRUE;
    if(!(h2 = (TH2F*)(fContainer->At(kdEdx)))) break;
    legdEdx->SetHeader("Particle Species");
    gPad->SetMargin(0.1, 0.01, 0.1, 0.01);
    for(Int_t is = AliPID::kSPECIES-1; is>=0; is--){
      Int_t bin = FindBin(is, 2.);
      h1 = h2->ProjectionY(Form("px%d", is), bin, bin);
      if(!h1->GetEntries()) continue;
      h1->Scale(1./h1->Integral());
      h1->SetLineColor(AliTRDCalPID::GetPartColor(is));
      if(kFIRST){
        h1->GetXaxis()->SetTitle("dE/dx (a.u.)");
        h1->GetYaxis()->SetTitle("<Entries>");
      }
      h = (TH1F*)h1->DrawClone(kFIRST ? "c" : "samec");
      legdEdx->AddEntry(h, Form("%s", AliTRDCalPID::GetPartName(is)), "l");
      kFIRST = kFALSE;
    }
    if(kFIRST) break;
    legdEdx->Draw();
    gPad->SetLogy();
    gPad->SetLogx(0);
    gPad->SetGridy();
    gPad->SetGridx();
    return kTRUE;
  }
  case kdEdxSlice:
    break;
  case kPH:{
    gPad->Divide(2, 1, 1.e-5, 1.e-5);
    TList *l=gPad->GetListOfPrimitives();

    // save 2.0 GeV projection as reference
    TLegend *legPH = new TLegend(.4, .7, .68, .98);
    legPH->SetBorderSize(1);legPH->SetFillColor(0);
    legPH->SetHeader("Particle Species");
    if(!(arr = (TObjArray*)(fContainer->At(kPH)))) break;
    if(!(h2 = (TProfile2D*)(arr->At(0)))) break;

    TVirtualPad *pad = ((TVirtualPad*)l->At(0));pad->cd();
    pad->SetMargin(0.1, 0.01, 0.1, 0.01);
    kFIRST = kTRUE;
    for(Int_t is=0; is<AliPID::kSPECIES; is++){
      Int_t bin = FindBin(is, 2.);
      h1 = h2->ProjectionY(Form("pyt%d", is), bin, bin);
      if(!h1->GetEntries()) continue;
      h1->SetMarkerStyle(24);
      h1->SetMarkerColor(AliTRDCalPID::GetPartColor(is));
      h1->SetLineColor(AliTRDCalPID::GetPartColor(is));
      if(kFIRST){
        h1->GetXaxis()->SetTitle("t_{drift} [100*ns]");
        h1->GetYaxis()->SetTitle("<dQ/dt> [a.u.]");
      }
      h = (TH1F*)h1->DrawClone(kFIRST ? "c" : "samec");
      legPH->AddEntry(h, Form("%s", AliTRDCalPID::GetPartName(is)), "pl");
      kFIRST = kFALSE;
    }

    pad = ((TVirtualPad*)l->At(1));pad->cd();
    pad->SetMargin(0.1, 0.01, 0.1, 0.01);
    if(!(h2 = (TProfile2D*)(arr->At(1)))) break;
    kFIRST = kTRUE;
    for(Int_t is=0; is<AliPID::kSPECIES; is++){
      Int_t bin = FindBin(is, 2.);
      h1 = h2->ProjectionY(Form("pyx%d", is), bin, bin);
      if(!h1->GetEntries()) continue;
      h1->SetMarkerStyle(24);
      h1->SetMarkerColor(AliTRDCalPID::GetPartColor(is));
      h1->SetLineColor(AliTRDCalPID::GetPartColor(is));
      if(kFIRST){
        h1->GetXaxis()->SetTitle("x_{drift} [cm]");
        h1->GetYaxis()->SetTitle("<dQ/dl> [a.u./cm]");
      }
      h1->DrawClone(kFIRST ? "c" : "samec");
      kFIRST = kFALSE;
    }

    if(kFIRST) break;
    legPH->Draw();
    gPad->SetLogy(0);
    gPad->SetLogx(0);
    gPad->SetGridy();
    gPad->SetGridx();
    return kTRUE;
  }
  case kNClus:{
    // save 2.0 GeV projection as reference
    TLegend *legNClus = new TLegend(.13, .7, .4, .98);
    legNClus->SetBorderSize(1);
    legNClus->SetFillColor(0);

    kFIRST = kTRUE;
    if(!(h2 = (TH2F*)(fContainer->At(kNClus)))) break;
    legNClus->SetHeader("Particle Species");
    for(Int_t is=0; is<AliPID::kSPECIES; is++){
      Int_t bin = FindBin(is, 2.);
      h1 = h2->ProjectionY(Form("pyNClus%d", is), bin, bin);
      if(!h1->GetEntries()) continue;
      h1->Scale(100./h1->Integral());
      //h1->SetMarkerStyle(24);
      //h1->SetMarkerColor(AliTRDCalPID::GetPartColor(is));
      h1->SetLineColor(AliTRDCalPID::GetPartColor(is));
      if(kFIRST){ 
        h1->GetXaxis()->SetTitle("N^{cl}/tracklet");
        h1->GetYaxis()->SetTitle("Prob. [%]");
        h = (TH1F*)h1->DrawClone("c");
        h->SetMaximum(20.);
        h->GetXaxis()->SetRangeUser(0., 35.);
        kFIRST = kFALSE;
      } else h = (TH1F*)h1->DrawClone("samec");

      legNClus->AddEntry(h, Form("%s", AliTRDCalPID::GetPartName(is)), "l");
    }
    if(kFIRST) break;
    legNClus->Draw();
    gPad->SetLogy(0);
    gPad->SetLogx(0);
    gPad->SetGridy();
    gPad->SetGridx();
    return kTRUE;
  }
  case kMomentum:
  case kMomentumBin:
    break; 
  case kNTracklets:{
    TLegend *legNClus = new TLegend(.4, .7, .68, .98);
    legNClus->SetBorderSize(1);
    kFIRST = kTRUE;
    if(!(h2 = (TH2F*)(fContainer->At(kNTracklets)))) break;
    legNClus->SetHeader("Particle Species");
    for(Int_t is=0; is<AliPID::kSPECIES; is++){
      Int_t bin = FindBin(is, 2.);
      h1 = h2->ProjectionY(Form("pyNTracklets%d", is), bin, bin);
      if(!h1->GetEntries()) continue;
      h1->Scale(100./h1->Integral());
      //h1->SetMarkerStyle(24);
      //h1->SetMarkerColor(AliTRDCalPID::GetPartColor(is));
      h1->SetLineColor(AliTRDCalPID::GetPartColor(is));
      if(kFIRST){ 
        h1->GetXaxis()->SetTitle("N^{trklt}/track");
        h1->GetXaxis()->SetRangeUser(1.,6.);
        h1->GetYaxis()->SetTitle("Prob. [%]");
      }
      h = (TH1F*)h1->DrawClone(kFIRST ? "c" : "samec");
      legNClus->AddEntry(h, Form("%s", AliTRDCalPID::GetPartName(is)), "l");
      kFIRST = kFALSE;
    }
    if(kFIRST) break;
    legNClus->Draw();
    gPad->SetLogy(0);
    gPad->SetLogx(0);
    gPad->SetGridy();
    gPad->SetGridx();
    return kTRUE;
  }
  }
  AliInfo(Form("Reference plot [%d] missing result", ifig));
  return kFALSE;
}

//________________________________________________________________________
Bool_t AliTRDcheckPID::PostProcess()
{
  // Draw result to the screen
  // Called once at the end of the query

  if (!fContainer) {
    Printf("ERROR: list not available");
    return kFALSE;
  }

  TObjArray *eff(NULL);
  if(!fGraph){ 
    fGraph = new TObjArray(2*AliPID::kSPECIES);
    fGraph->SetOwner();
    
    if(!(eff = dynamic_cast<TObjArray *>(fContainer->At(kEfficiency)))){
      AliError("Efficiency container for Electrons missing.");
      return kFALSE;
    }
    EvaluateEfficiency(eff, fGraph, AliPID::kPion, 0.9);
    EvaluateEfficiency(eff, fGraph, AliPID::kKaon, 0.9);
    EvaluateEfficiency(eff, fGraph, AliPID::kProton, 0.9);
  }
  fNRefFigures = 12;
  return kTRUE;
}

//________________________________________________________________________
void AliTRDcheckPID::EvaluateEfficiency(const TObjArray * const histoContainer, TObjArray *results, Int_t species, Float_t electronEfficiency){
// Process PID information for pion efficiency

  fUtil->SetElectronEfficiency(electronEfficiency);


  const Int_t kNmethodsPID=Int_t(sizeof(fgMethod)/sizeof(Char_t*));
  Color_t colors[kNmethodsPID] = {kBlue, kGreen+2, kRed};
  Int_t markerStyle[kNmethodsPID] = {7, 7, 24};
  // efficiency graphs
  TGraphErrors *g(NULL);
  TObjArray *eff = new TObjArray(kNmethodsPID); eff->SetOwner(); eff->SetName(Form("Eff_%s", AliTRDCalPID::GetPartName(species)));
  results->AddAt(eff, species);
  for(Int_t iMethod = 0; iMethod < kNmethodsPID; iMethod++){
    eff->AddAt(g = new TGraphErrors(), iMethod);
    g->SetName(Form("%s", fgMethod[iMethod]));
    g->SetLineColor(colors[iMethod]);
    g->SetMarkerColor(colors[iMethod]);
    g->SetMarkerStyle(markerStyle[iMethod]);
  }

  // Threshold graphs if not already
  TObjArray *thres(NULL);
  if(!(results->At(AliPID::kSPECIES))){
    thres = new TObjArray(kNmethodsPID); thres->SetOwner(); 
    thres->SetName("Thres");
    results->AddAt(thres, AliPID::kSPECIES);
    for(Int_t iMethod = 0; iMethod < kNmethodsPID; iMethod++){
      thres->AddAt(g = new TGraphErrors(), iMethod);
      g->SetName(Form("%s", fgMethod[iMethod]));
      g->SetLineColor(colors[iMethod]);
      g->SetMarkerColor(colors[iMethod]);
      g->SetMarkerStyle(markerStyle[iMethod]);
    }
  }

  TH2F *hPtr[kNmethodsPID]={
    (TH2F*)histoContainer->At(AliTRDpidUtil::kLQ),
    (TH2F*)histoContainer->At(AliTRDpidUtil::kNN),
    (TH2F*)histoContainer->At(AliTRDpidUtil::kESD)
  };
  for(Int_t iMom = 0; iMom < fMomentumAxis->GetNbins(); iMom++){
    Float_t mom(fMomentumAxis->GetBinCenter(iMom+1));

    Int_t binEl(fMomentumAxis->GetNbins() * AliPID::kElectron + iMom + 1), 
	        binXX(fMomentumAxis->GetNbins() * species + iMom + 1);

    for(Int_t iMethod = 0; iMethod < kNmethodsPID; iMethod++){
      // Calculate the Species Efficiency at electronEfficiency% electron efficiency for each Method
  
      TH1D *histo1 = hPtr[iMethod] -> ProjectionY(Form("%s_el", fgMethod[iMethod]), binEl, binEl);
      TH1D *histo2 = hPtr[iMethod] -> ProjectionY(Form("%s_%s", fgMethod[iMethod], AliTRDCalPID::GetPartName(species)), binXX, binXX);

      if(!fUtil->CalculatePionEffi(histo1, histo2)) continue;
     
      g=(TGraphErrors*)eff->At(iMethod);
      g->SetPoint(iMom, mom, 1.e2*fUtil->GetPionEfficiency());
      g->SetPointError(iMom, 0., 1.e2*fUtil->GetError());
      AliDebug(2, Form("%s Efficiency for %s is : %f +/- %f", AliTRDCalPID::GetPartName(species), fgMethod[iMethod], fUtil->GetPionEfficiency(), fUtil->GetError()));

      if(!thres) continue;
      g=(TGraphErrors*)thres->At(iMethod);
      g->SetPoint(iMom, mom, fUtil->GetThreshold());
      g->SetPointError(iMom, 0., 0.);
    }
  }
}
