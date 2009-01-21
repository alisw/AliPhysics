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

#include "AliTRDpidChecker.h"
#include "AliTRDtrackInfo/AliTRDtrackInfo.h"

// calculate pion efficiency at 90% electron efficiency for 11 momentum bins
// this task should be used with simulated data only

ClassImp(AliTRDpidChecker)

//________________________________________________________________________
AliTRDpidChecker::AliTRDpidChecker() 
  :AliTRDrecoTask("PID", "PID Checker")
  ,fReconstructor(0x0)
  ,fUtil(0x0)
  ,fGraph(0x0)
  ,fEfficiency(0x0)
{
  //
  // Default constructor
  //

  fReconstructor = new AliTRDReconstructor();
  fReconstructor->SetRecoParam(AliTRDrecoParam::GetLowFluxParam());

  fUtil = new AliTRDpidUtil();

  InitFunctorList();
}


//________________________________________________________________________
AliTRDpidChecker::~AliTRDpidChecker() 
{
 if(fGraph){fGraph->Delete(); delete fGraph;}
 if(fReconstructor) delete fReconstructor;
  if(fUtil) delete fUtil;
}


//________________________________________________________________________
void AliTRDpidChecker::CreateOutputObjects()
{
  // Create histograms
  // Called once

  OpenFile(0, "RECREATE");
  fContainer = Histos();
}


//_______________________________________________________
TObjArray * AliTRDpidChecker::Histos(){

  //
  // Create QA histograms
  //
  if(fContainer) return fContainer;

  Int_t xBins = AliPID::kSPECIES*AliTRDCalPID::kNMom; 
  fContainer = new TObjArray(); fContainer->Expand(7);

  const Float_t epsilon = 1./(2*(AliTRDpidUtil::kBins-1));     // get nice histos with bin center at 0 and 1

  // histos of the electron probability of all 5 particle species and 11 momenta for the 2-dim LQ method 
  fEfficiency = new TObjArray(3);
  fEfficiency->SetOwner(); fEfficiency->SetName("Efficiency");
  fContainer->AddAt(fEfficiency, kEfficiency);
  
  TH1 *h = 0x0;
  if(!(h = (TH2F*)gROOT->FindObject("PID_LQ"))){
    h = new TH2F("PID_LQ", "", xBins, -0.5, xBins - 0.5,
      AliTRDpidUtil::kBins, 0.-epsilon, 1.+epsilon);
  } else h->Reset();
  fEfficiency->AddAt(h, kLQ);

  // histos of the electron probability of all 5 particle species and 11 momenta for the neural network method
  if(!(h = (TH2F*)gROOT->FindObject("PID_NN"))){
    h = new TH2F("PID_NN", "", 
      xBins, -0.5, xBins - 0.5,
      AliTRDpidUtil::kBins, 0.-epsilon, 1.+epsilon);
  } else h->Reset();
  fEfficiency->AddAt(h, kNN);

  // histos of the electron probability of all 5 particle species and 11 momenta for the ESD output
  if(!(h = (TH2F*)gROOT->FindObject("PID_ESD"))){
    h = new TH2F("PID_ESD", "", 
      xBins, -0.5, xBins - 0.5,
      AliTRDpidUtil::kBins, 0.-epsilon, 1.+epsilon);
  } else h->Reset();
  fEfficiency->AddAt(h, kESD);

  // histos of the dE/dx distribution for all 5 particle species and 11 momenta 
  if(!(h = (TH2F*)gROOT->FindObject("dEdx"))){
    h = new TH2F("dEdx", "", 
      xBins, -0.5, xBins - 0.5,
      200, 0, 10000);
  } else h->Reset();
  fContainer->AddAt(h, kdEdx);

  // histos of the dE/dx slices for all 5 particle species and 11 momenta 
  if(!(h = (TH2F*)gROOT->FindObject("dEdxSlice"))){
    h = new TH2F("dEdxSlice", "", 
      xBins*AliTRDReconstructor::kLQslices, -0.5, xBins*AliTRDReconstructor::kLQslices - 0.5,
      200, 0, 5000);
  } else h->Reset();
  fContainer->AddAt(h, kdEdxSlice);

  // histos of the pulse height distribution for all 5 particle species and 11 momenta 
  if(!(h = (TH2F*)gROOT->FindObject("PH"))){
    h = new TProfile2D("PH", "", 
      xBins, -0.5, xBins - 0.5,
      AliTRDtrackerV1::GetNTimeBins(), -0.5, AliTRDtrackerV1::GetNTimeBins() - 0.5);
  } else h->Reset();
  fContainer->AddAt(h, kPH);

  // histos of the number of clusters distribution for all 5 particle species and 11 momenta 
  if(!(h = (TH2F*)gROOT->FindObject("NClus"))){
    h = new TH2F("NClus", "", 
      xBins, -0.5, xBins - 0.5,
      AliTRDtrackerV1::GetNTimeBins(), -0.5, AliTRDtrackerV1::GetNTimeBins() - 0.5);
  } else h->Reset();
  fContainer->AddAt(h, kNClus);


  // momentum distributions - absolute and in momentum bins
  if(!(h = (TH1F*)gROOT->FindObject("hMom"))){
    h = new TH1F("hMom", "momentum distribution", 100, 0., 12.);
  } else h->Reset();
  fContainer->AddAt(h, kMomentum);
  
  if(!(h = (TH1F*)gROOT->FindObject("hMomBin"))){
    h = new TH1F("hMomBin", "momentum distribution in momentum bins", AliTRDCalPID::kNMom, 0.5, 11.5);
  } else h->Reset();
  fContainer->AddAt(h, kMomentumBin);


  return fContainer;
}


//________________________________________________________________________
Bool_t AliTRDpidChecker::CheckTrackQuality(const AliTRDtrackV1* track) 
{
  //
  // Check if the track is ok for PID
  //
  
  if(track->GetNumberOfTracklets() == AliTRDgeometry::kNlayer) return 1;
//   if(!fESD)
//     return 0;

  return 0;
}


//________________________________________________________________________
Int_t AliTRDpidChecker::CalcPDG(AliTRDtrackV1* track) 
{

  track -> SetReconstructor(fReconstructor);

  fReconstructor -> SetOption("nn");
  track -> CookPID();

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
TH1 *AliTRDpidChecker::PlotLQ(const AliTRDtrackV1 *track)
{
  //
  // Plot the probabilities for electrons using 2-dim LQ
  //

  if(track) fTrack = track;
  if(!fTrack){
    AliWarning("No Track defined.");
    return 0x0;
  }

  if(!CheckTrackQuality(fTrack)) return 0x0;
  
  if(!(fEfficiency = dynamic_cast<TObjArray *>(fContainer->At(kEfficiency)))){
    AliWarning("No Histogram defined.");
    return 0x0;
  }
  TH2F *hPIDLQ = 0x0;
  if(!(hPIDLQ = dynamic_cast<TH2F *>(fEfficiency->At(kLQ)))){
    AliWarning("No Histogram defined.");
    return 0x0;
  }

  AliTRDtrackV1 cTrack(*fTrack);
  cTrack.SetReconstructor(fReconstructor);

  Int_t pdg = 0;
  Float_t momentum = 0.;

  if(fMC){
    if(fMC->GetTrackRef()) momentum = fMC->GetTrackRef()->P();
    pdg = fMC->GetPDG();
  } else{
    //AliWarning("No MC info available!");
    momentum = cTrack.GetMomentum(0);
    pdg = CalcPDG(&cTrack);
  }
  if(momentum < 0.4) return 0x0;;
  if(momentum > 12.) return 0x0;;

  fReconstructor -> SetOption("!nn");
  cTrack.CookPID();
  Int_t iMomBin = fUtil->GetMomentumBin(momentum);

  switch(pdg){
  case kElectron:
  case kPositron:
    hPIDLQ -> Fill(AliPID::kElectron * AliTRDCalPID::kNMom + iMomBin, cTrack.GetPID(AliPID::kElectron));
    break;
  case kMuonPlus:
  case kMuonMinus:
    hPIDLQ -> Fill(AliPID::kMuon * AliTRDCalPID::kNMom + iMomBin, cTrack .GetPID(AliPID::kElectron));
    break;
  case kPiPlus:
  case kPiMinus:
    hPIDLQ -> Fill(AliPID::kPion * AliTRDCalPID::kNMom + iMomBin, cTrack .GetPID(AliPID::kElectron));
    break;
  case kKPlus:
  case kKMinus:
    hPIDLQ -> Fill(AliPID::kKaon * AliTRDCalPID::kNMom + iMomBin, cTrack .GetPID(AliPID::kElectron));
    break;
  case kProton:
  case kProtonBar:
    hPIDLQ -> Fill(AliPID::kProton * AliTRDCalPID::kNMom + iMomBin, cTrack.GetPID(AliPID::kElectron));
    break;
  }

  return hPIDLQ;
}


//_______________________________________________________
TH1 *AliTRDpidChecker::PlotNN(const AliTRDtrackV1 *track)
{
  //
  // Plot the probabilities for electrons using 2-dim LQ
  //

  if(track) fTrack = track;
  if(!fTrack){
    AliWarning("No Track defined.");
    return 0x0;
  }
  
  if(!CheckTrackQuality(fTrack)) return 0x0;
  
  if(!(fEfficiency = dynamic_cast<TObjArray *>(fContainer->At(kEfficiency)))){
    AliWarning("No Histogram defined.");
    return 0x0;
  }
  TH2F *hPIDNN;
  if(!(hPIDNN = dynamic_cast<TH2F *>(fEfficiency->At(kNN)))){
    AliWarning("No Histogram defined.");
    return 0x0;
  }


  AliTRDtrackV1 cTrack(*fTrack);
  cTrack.SetReconstructor(fReconstructor);

  Int_t pdg = 0;
  Float_t momentum = 0.;
  if(fMC){
    if(fMC->GetTrackRef()) momentum = fMC->GetTrackRef()->P();
    pdg = fMC->GetPDG();
  } else {
    //AliWarning("No MC info available!");
    momentum = cTrack.GetMomentum(0);
    pdg = CalcPDG(&cTrack);
  }
  if(momentum < 0.4) return 0x0;;
  if(momentum > 12.) return 0x0;;

  fReconstructor -> SetOption("nn");
  cTrack.CookPID();
  Int_t iMomBin = fUtil -> GetMomentumBin(momentum);

  switch(pdg){
  case kElectron:
  case kPositron:
    hPIDNN -> Fill(AliPID::kElectron * AliTRDCalPID::kNMom + iMomBin, cTrack.GetPID(AliPID::kElectron));
    break;
  case kMuonPlus:
  case kMuonMinus:
    hPIDNN -> Fill(AliPID::kMuon * AliTRDCalPID::kNMom + iMomBin, cTrack.GetPID(AliPID::kElectron));
    break;
  case kPiPlus:
  case kPiMinus:
    hPIDNN -> Fill(AliPID::kPion * AliTRDCalPID::kNMom + iMomBin, cTrack.GetPID(AliPID::kElectron));
    break;
  case kKPlus:
  case kKMinus:
    hPIDNN -> Fill(AliPID::kKaon * AliTRDCalPID::kNMom + iMomBin, cTrack.GetPID(AliPID::kElectron));
    break;
  case kProton:
  case kProtonBar:
    hPIDNN -> Fill(AliPID::kProton * AliTRDCalPID::kNMom + iMomBin, cTrack.GetPID(AliPID::kElectron));
    break;
  }
  return hPIDNN;
}


//_______________________________________________________
TH1 *AliTRDpidChecker::PlotESD(const AliTRDtrackV1 *track)
{
  //
  // Plot the probabilities for electrons using 2-dim LQ
  //

  if(!fESD){
    AliWarning("No ESD info available.");
    return 0x0;
  }

  if(track) fTrack = track;
  if(!fTrack){
    AliWarning("No Track defined.");
    return 0x0;
  }
  
  if(!CheckTrackQuality(fTrack)) return 0x0;
  
  if(!(fEfficiency = dynamic_cast<TObjArray *>(fContainer->At(kEfficiency)))){
    AliWarning("No Histogram defined.");
    return 0x0;
  }
  TH2F *hPIDESD = 0x0;
  if(!(hPIDESD = dynamic_cast<TH2F *>(fEfficiency->At(kESD)))){
    AliWarning("No Histogram defined.");
    return 0x0;
  }


  Int_t pdg = 0;
  Float_t momentum = 0.;
  if(fMC){
    if(fMC->GetTrackRef()) momentum = fMC->GetTrackRef()->P();
    pdg = fMC->GetPDG();
  } else {
    //AliWarning("No MC info available!");
    AliTRDtrackV1 cTrack(*fTrack);
    cTrack.SetReconstructor(fReconstructor);
    momentum = cTrack.GetMomentum(0);
    pdg = CalcPDG(&cTrack);
  }
  if(momentum < 0.4) return 0x0;;
  if(momentum > 12.) return 0x0;;
  

  Int_t iMomBin = fUtil->GetMomentumBin(momentum);


//   Double32_t pidESD[AliPID::kSPECIES];
  const Double32_t *pidESD = fESD->GetResponseIter();

  switch(pdg){
  case kElectron:
  case kPositron:
    hPIDESD -> Fill(AliPID::kElectron * AliTRDCalPID::kNMom + iMomBin, pidESD[0]);
    break;
  case kMuonPlus:
  case kMuonMinus:
    hPIDESD -> Fill(AliPID::kMuon * AliTRDCalPID::kNMom + iMomBin, pidESD[0]);
    break;
  case kPiPlus:
  case kPiMinus:
    hPIDESD -> Fill(AliPID::kPion * AliTRDCalPID::kNMom + iMomBin, pidESD[0]);
    break;
  case kKPlus:
  case kKMinus:
    hPIDESD -> Fill(AliPID::kKaon * AliTRDCalPID::kNMom + iMomBin, pidESD[0]);
    break;
  case kProton:
  case kProtonBar:
    hPIDESD -> Fill(AliPID::kProton * AliTRDCalPID::kNMom + iMomBin, pidESD[0]);
    break;
  }
  return hPIDESD;
}


//_______________________________________________________
TH1 *AliTRDpidChecker::PlotdEdx(const AliTRDtrackV1 *track)
{
  //
  // Plot the probabilities for electrons using 2-dim LQ
  //

  if(track) fTrack = track;
  if(!fTrack){
    AliWarning("No Track defined.");
    return 0x0;
  }
  
  if(!CheckTrackQuality(fTrack)) return 0x0;
  
  TH2F *hdEdx;
  if(!(hdEdx = dynamic_cast<TH2F *>(fContainer->At(kdEdx)))){
    AliWarning("No Histogram defined.");
    return 0x0;
  }


  Int_t pdg = 0;
  Float_t momentum = 0.;
  if(fMC){
    if(fMC->GetTrackRef()) momentum = fMC->GetTrackRef()->P();
    pdg = fMC->GetPDG();
  } else {
    //AliWarning("No MC info available!");
    AliTRDtrackV1 cTrack(*fTrack);
    cTrack.SetReconstructor(fReconstructor);
    momentum = cTrack.GetMomentum(0);
    pdg = CalcPDG(&cTrack);
  }
  if(momentum < 0.4) return 0x0;;
  if(momentum > 12.) return 0x0;;

  Int_t iMomBin = fUtil -> GetMomentumBin(momentum);



  Float_t SumdEdx[AliTRDgeometry::kNlayer];
  AliTRDseedV1 *TRDtracklet[AliTRDgeometry::kNlayer];
  for(Int_t iChamb = 0; iChamb < AliTRDgeometry::kNlayer; iChamb++) TRDtracklet[iChamb] = 0x0;
  for(Int_t iChamb = 0; iChamb < AliTRDgeometry::kNlayer; iChamb++) TRDtracklet[iChamb] = fTrack -> GetTracklet(iChamb);

  Float_t *fdEdx;
  Float_t dEdxSlice[AliTRDgeometry::kNlayer][AliTRDReconstructor::kLQslices];

  for(Int_t iChamb = 0; iChamb < AliTRDgeometry::kNlayer; iChamb++){
    SumdEdx[iChamb] = 0.;
    fdEdx = TRDtracklet[iChamb] -> GetdEdx();
    SumdEdx[iChamb] += fdEdx[0] + fdEdx[1] + fdEdx[2];
    for(Int_t iSlice = 0; iSlice < AliTRDReconstructor::kLQslices; iSlice++){
      dEdxSlice[iChamb][iSlice] = fdEdx[iSlice];
    }
  }

  switch(pdg){
  case kElectron:
  case kPositron:
    for(Int_t iChamb = 0; iChamb < AliTRDgeometry::kNlayer; iChamb++)
      hdEdx -> Fill(AliPID::kElectron * AliTRDCalPID::kNMom + iMomBin, SumdEdx[iChamb]);
    break;
  case kMuonPlus:
  case kMuonMinus:
    for(Int_t iChamb = 0; iChamb < AliTRDgeometry::kNlayer; iChamb++)
      hdEdx -> Fill(AliPID::kMuon * AliTRDCalPID::kNMom + iMomBin, SumdEdx[iChamb]);
    break;
  case kPiPlus:
  case kPiMinus:
    for(Int_t iChamb = 0; iChamb < AliTRDgeometry::kNlayer; iChamb++)
      hdEdx -> Fill(AliPID::kPion * AliTRDCalPID::kNMom + iMomBin, SumdEdx[iChamb]);
    break;
  case kKPlus:
  case kKMinus:
    for(Int_t iChamb = 0; iChamb < AliTRDgeometry::kNlayer; iChamb++)
      hdEdx -> Fill(AliPID::kKaon * AliTRDCalPID::kNMom + iMomBin, SumdEdx[iChamb]);
    break;
  case kProton:
  case kProtonBar:
    for(Int_t iChamb = 0; iChamb < AliTRDgeometry::kNlayer; iChamb++)
      hdEdx -> Fill(AliPID::kProton * AliTRDCalPID::kNMom + iMomBin, SumdEdx[iChamb]);
    break;
  }

  return hdEdx;
}


//_______________________________________________________
TH1 *AliTRDpidChecker::PlotdEdxSlice(const AliTRDtrackV1 *track)
{
  //
  // Plot the probabilities for electrons using 2-dim LQ
  //

  if(track) fTrack = track;
  if(!fTrack){
    AliWarning("No Track defined.");
    return 0x0;
  }
  
  if(!CheckTrackQuality(fTrack)) return 0x0;
  
  TH2F *hdEdxSlice;
  if(!(hdEdxSlice = dynamic_cast<TH2F *>(fContainer->At(kdEdxSlice)))){
    AliWarning("No Histogram defined.");
    return 0x0;
  }

  Int_t pdg = 0;
  Float_t momentum = 0.;
  if(fMC){
    if(fMC->GetTrackRef()) momentum = fMC->GetTrackRef()->P();
    pdg = fMC->GetPDG();
  } else {
    //AliWarning("No MC info available!");
    AliTRDtrackV1 cTrack(*fTrack);
    cTrack.SetReconstructor(fReconstructor);
    momentum = cTrack.GetMomentum(0);
    pdg = CalcPDG(&cTrack);
  }
  if(momentum < 0.4) return 0x0;
  if(momentum > 12.) return 0x0;;

  Int_t iMomBin = fUtil -> GetMomentumBin(momentum);



  AliTRDseedV1 *TRDtracklet[AliTRDgeometry::kNlayer];
  for(Int_t iChamb = 0; iChamb < AliTRDgeometry::kNlayer; iChamb++) TRDtracklet[iChamb] = 0x0;
  for(Int_t iChamb = 0; iChamb < AliTRDgeometry::kNlayer; iChamb++) TRDtracklet[iChamb] = fTrack -> GetTracklet(iChamb);

  Float_t *fdEdx;
  Float_t dEdxSlice[AliTRDgeometry::kNlayer][AliTRDReconstructor::kLQslices];

  for(Int_t iChamb = 0; iChamb < AliTRDgeometry::kNlayer; iChamb++){
    fdEdx = TRDtracklet[iChamb] -> GetdEdx();
    for(Int_t iSlice = 0; iSlice < AliTRDReconstructor::kLQslices; iSlice++){
      dEdxSlice[iChamb][iSlice] = fdEdx[iSlice];
    }
  }

  switch(pdg){
  case kElectron:
  case kPositron:
    for(Int_t iChamb = 0; iChamb < AliTRDgeometry::kNlayer; iChamb++){
      for(Int_t iSlice = 0; iSlice < AliTRDReconstructor::kLQslices; iSlice++){
        hdEdxSlice -> Fill(AliPID::kElectron * AliTRDCalPID::kNMom * AliTRDReconstructor::kLQslices + iMomBin * AliTRDReconstructor::kLQslices + iSlice, dEdxSlice[iChamb][iSlice]);
      }
    }  
    break;
  case kMuonPlus:
  case kMuonMinus:
    for(Int_t iChamb = 0; iChamb < AliTRDgeometry::kNlayer; iChamb++){
      for(Int_t iSlice = 0; iSlice < AliTRDReconstructor::kLQslices; iSlice++){
        hdEdxSlice -> Fill(AliPID::kMuon * AliTRDCalPID::kNMom * AliTRDReconstructor::kLQslices + iMomBin * AliTRDReconstructor::kLQslices + iSlice,
        dEdxSlice[iChamb][iSlice]);
      }
    }
    break;
  case kPiPlus:
  case kPiMinus:
    for(Int_t iChamb = 0; iChamb < AliTRDgeometry::kNlayer; iChamb++)
      for(Int_t iSlice = 0; iSlice < AliTRDReconstructor::kLQslices; iSlice++)
	hdEdxSlice -> Fill(AliPID::kPion * AliTRDCalPID::kNMom * AliTRDReconstructor::kLQslices + iMomBin * AliTRDReconstructor::kLQslices + iSlice,
			   dEdxSlice[iChamb][iSlice]);
    break;
  case kKPlus:
  case kKMinus:
    for(Int_t iChamb = 0; iChamb < AliTRDgeometry::kNlayer; iChamb++)
      for(Int_t iSlice = 0; iSlice < AliTRDReconstructor::kLQslices; iSlice++)
	hdEdxSlice -> Fill(AliPID::kKaon * AliTRDCalPID::kNMom * AliTRDReconstructor::kLQslices + iMomBin * AliTRDReconstructor::kLQslices + iSlice,
			   dEdxSlice[iChamb][iSlice]);
    break;
  case kProton:
  case kProtonBar:
    for(Int_t iChamb = 0; iChamb < AliTRDgeometry::kNlayer; iChamb++)
      for(Int_t iSlice = 0; iSlice < AliTRDReconstructor::kLQslices; iSlice++)
	hdEdxSlice -> Fill(AliPID::kProton * AliTRDCalPID::kNMom * AliTRDReconstructor::kLQslices + iMomBin * AliTRDReconstructor::kLQslices + iSlice,
			   dEdxSlice[iChamb][iSlice]);
    break;
  }

  return hdEdxSlice;

}


//_______________________________________________________
TH1 *AliTRDpidChecker::PlotPH(const AliTRDtrackV1 *track)
{
  //
  // Plot the probabilities for electrons using 2-dim LQ
  //

  if(track) fTrack = track;
  if(!fTrack){
    AliWarning("No Track defined.");
    return 0x0;
  }
  
  if(!CheckTrackQuality(fTrack)) return 0x0;
  
  TProfile2D *hPH;
  if(!(hPH = dynamic_cast<TProfile2D *>(fContainer->At(kPH)))){
    AliWarning("No Histogram defined.");
    return 0x0;
  }

  Int_t pdg = 0;
  Float_t momentum = 0.;
  if(fMC){
    if(fMC->GetTrackRef()) momentum = fMC->GetTrackRef()->P();
    pdg = fMC->GetPDG();
  } else {
    //AliWarning("No MC info available!");
    AliTRDtrackV1 cTrack(*fTrack);
    cTrack.SetReconstructor(fReconstructor);
    momentum = cTrack.GetMomentum(0);
    pdg = CalcPDG(&cTrack);
  }
  if(momentum < 0.4) return 0x0;;
  if(momentum > 12.) return 0x0;;

  Int_t iMomBin = fUtil -> GetMomentumBin(momentum);

  AliTRDseedV1 *TRDtracklet[AliTRDgeometry::kNlayer];
  for(Int_t iChamb = 0; iChamb < AliTRDgeometry::kNlayer; iChamb++) TRDtracklet[iChamb] = 0x0;
  for(Int_t iChamb = 0; iChamb < AliTRDgeometry::kNlayer; iChamb++) TRDtracklet[iChamb] = fTrack -> GetTracklet(iChamb);

  AliTRDcluster *TRDcluster = 0x0;

  switch(pdg){
  case kElectron:
  case kPositron:
    for(Int_t iChamb = 0; iChamb < AliTRDgeometry::kNlayer; iChamb++){
      for(Int_t iClus = 0; iClus < AliTRDtrackerV1::GetNTimeBins(); iClus++){
	if(!(TRDcluster = (AliTRDcluster*)TRDtracklet[iChamb] -> GetClusters(iClus)))
	  continue;
	hPH -> Fill(AliPID::kElectron * AliTRDCalPID::kNMom + iMomBin, TRDcluster -> GetLocalTimeBin(), TRDtracklet[iChamb] -> GetdQdl(iClus));
      }
    }
    break;
  case kMuonPlus:
  case kMuonMinus:
    for(Int_t iChamb = 0; iChamb < AliTRDgeometry::kNlayer; iChamb++){
      for(Int_t iClus = 0; iClus < AliTRDtrackerV1::GetNTimeBins(); iClus++){
	if(!(TRDcluster = (AliTRDcluster*)TRDtracklet[iChamb] -> GetClusters(iClus)))
	  continue;
	hPH -> Fill(AliPID::kMuon * AliTRDCalPID::kNMom + iMomBin, TRDcluster -> GetLocalTimeBin(), TRDtracklet[iChamb] -> GetdQdl(iClus));
      }
    }
    break;
  case kPiPlus:
  case kPiMinus:
    for(Int_t iChamb = 0; iChamb < AliTRDgeometry::kNlayer; iChamb++){
      for(Int_t iClus = 0; iClus < AliTRDtrackerV1::GetNTimeBins(); iClus++){
	if(!(TRDcluster = (AliTRDcluster*)TRDtracklet[iChamb] -> GetClusters(iClus)))
	  continue;
	hPH -> Fill(AliPID::kPion * AliTRDCalPID::kNMom + iMomBin, TRDcluster -> GetLocalTimeBin(), TRDtracklet[iChamb] -> GetdQdl(iClus));
      }
    }
    break;
  case kKPlus:
  case kKMinus:
    for(Int_t iChamb = 0; iChamb < AliTRDgeometry::kNlayer; iChamb++){
      for(Int_t iClus = 0; iClus < AliTRDtrackerV1::GetNTimeBins(); iClus++){
	if(!(TRDcluster = (AliTRDcluster*)TRDtracklet[iChamb] -> GetClusters(iClus)))
	  continue;
	hPH -> Fill(AliPID::kKaon * AliTRDCalPID::kNMom + iMomBin, TRDcluster -> GetLocalTimeBin(), TRDtracklet[iChamb] -> GetdQdl(iClus));
      }
    }
    break;
  case kProton:
  case kProtonBar:
    for(Int_t iChamb = 0; iChamb < AliTRDgeometry::kNlayer; iChamb++){
      for(Int_t iClus = 0; iClus < AliTRDtrackerV1::GetNTimeBins(); iClus++){
	if(!(TRDcluster = (AliTRDcluster*)TRDtracklet[iChamb] -> GetClusters(iClus)))
	  continue;
	hPH -> Fill(AliPID::kProton * AliTRDCalPID::kNMom + iMomBin, TRDcluster -> GetLocalTimeBin(), TRDtracklet[iChamb] -> GetdQdl(iClus));
      }
    }
    break; 
  }
  
  return hPH;
}


//_______________________________________________________
TH1 *AliTRDpidChecker::PlotNClus(const AliTRDtrackV1 *track)
{
  //
  // Plot the probabilities for electrons using 2-dim LQ
  //

  if(track) fTrack = track;
  if(!fTrack){
    AliWarning("No Track defined.");
    return 0x0;
  }
  
  if(!CheckTrackQuality(fTrack)) return 0x0;
  
  TH2F *hNClus;
  if(!(hNClus = dynamic_cast<TH2F *>(fContainer->At(kNClus)))){
    AliWarning("No Histogram defined.");
    return 0x0;
  }


  Int_t pdg = 0;
  Float_t momentum = 0.;
  if(fMC){
    if(fMC->GetTrackRef()) momentum = fMC->GetTrackRef()->P();
    pdg = fMC->GetPDG();
  } else {
    //AliWarning("No MC info available!");
    AliTRDtrackV1 cTrack(*fTrack);
    cTrack.SetReconstructor(fReconstructor);
    momentum = cTrack.GetMomentum(0);
    pdg = CalcPDG(&cTrack);
  }
  if(momentum < 0.4) return 0x0;;
  if(momentum > 12.) return 0x0;;

  Int_t iMomBin = fUtil -> GetMomentumBin(momentum);


  Int_t iNClus[AliTRDgeometry::kNlayer]; 
  memset(iNClus, 0, sizeof(Int_t) * AliTRDgeometry::kNlayer);

  AliTRDseedV1 *TRDtracklet[AliTRDgeometry::kNlayer];
  for(Int_t iChamb = 0; iChamb < AliTRDgeometry::kNlayer; iChamb++) TRDtracklet[iChamb] = 0x0;
  for(Int_t iChamb = 0; iChamb < AliTRDgeometry::kNlayer; iChamb++){
    TRDtracklet[iChamb] = fTrack -> GetTracklet(iChamb);
    iNClus[iChamb] = TRDtracklet[iChamb] -> GetN();
  }

  switch(pdg){
  case kElectron:
  case kPositron:
    for(Int_t iChamb = 0; iChamb < AliTRDgeometry::kNlayer; iChamb++)
      hNClus -> Fill(AliPID::kElectron * AliTRDCalPID::kNMom + iMomBin, iNClus[iChamb]);
    break;
  case kMuonPlus:
  case kMuonMinus:
    for(Int_t iChamb = 0; iChamb < AliTRDgeometry::kNlayer; iChamb++)
      hNClus -> Fill(AliPID::kMuon * AliTRDCalPID::kNMom + iMomBin, iNClus[iChamb]);
    break;
  case kPiPlus:
  case kPiMinus:
    for(Int_t iChamb = 0; iChamb < AliTRDgeometry::kNlayer; iChamb++)
      hNClus -> Fill(AliPID::kPion * AliTRDCalPID::kNMom + iMomBin, iNClus[iChamb]);
    break;
  case kKPlus:
  case kKMinus:
    for(Int_t iChamb = 0; iChamb < AliTRDgeometry::kNlayer; iChamb++)
      hNClus -> Fill(AliPID::kKaon * AliTRDCalPID::kNMom + iMomBin, iNClus[iChamb]);
    break;
  case kProton:
  case kProtonBar:
    for(Int_t iChamb = 0; iChamb < AliTRDgeometry::kNlayer; iChamb++)
      hNClus -> Fill(AliPID::kProton * AliTRDCalPID::kNMom + iMomBin, iNClus[iChamb]);
    break;
  }

  return hNClus;
}


//_______________________________________________________
TH1 *AliTRDpidChecker::PlotMom(const AliTRDtrackV1 *track)
{
  //
  // Plot the probabilities for electrons using 2-dim LQ
  //

  if(track) fTrack = track;
  if(!fTrack){
    AliWarning("No Track defined.");
    return 0x0;
  }
  
  if(!CheckTrackQuality(fTrack)) return 0x0;
  
  TH1F *hMom;
  if(!(hMom = dynamic_cast<TH1F *>(fContainer->At(kMomentum)))){
    AliWarning("No Histogram defined.");
    return 0x0;
  }


  Int_t pdg = 0;
  Float_t momentum = 0.;
  if(fMC){
    if(fMC->GetTrackRef()) momentum = fMC->GetTrackRef()->P();
    pdg = fMC->GetPDG();
  } else {
    //AliWarning("No MC info available!");
    AliTRDtrackV1 cTrack(*fTrack);
    cTrack.SetReconstructor(fReconstructor);
    momentum = cTrack.GetMomentum(0);
    pdg = CalcPDG(&cTrack);
  }
  if(momentum < 0.4) return 0x0;
  if(momentum > 12.) return 0x0;;

  hMom -> Fill(momentum);
  return hMom;
}


//_______________________________________________________
TH1 *AliTRDpidChecker::PlotMomBin(const AliTRDtrackV1 *track)
{
  //
  // Plot the probabilities for electrons using 2-dim LQ
  //

  if(track) fTrack = track;
  if(!fTrack){
    AliWarning("No Track defined.");
    return 0x0;
  }
  
  if(!CheckTrackQuality(fTrack)) return 0x0;
  
  TH1F *hMomBin;
  if(!(hMomBin = dynamic_cast<TH1F *>(fContainer->At(kMomentumBin)))){
    AliWarning("No Histogram defined.");
    return 0x0;
  }


  Int_t pdg = 0;
  Float_t momentum = 0.;

  if(fMC){
    if(fMC->GetTrackRef()) momentum = fMC->GetTrackRef()->P();
    pdg = fMC->GetPDG();
  } else {
    //AliWarning("No MC info available!");
    AliTRDtrackV1 cTrack(*fTrack);
    cTrack.SetReconstructor(fReconstructor);
    momentum = cTrack.GetMomentum(0);
    pdg = CalcPDG(&cTrack);
  }
  if(momentum < 0.4) return 0x0;
  if(momentum > 12.) return 0x0;;

  Int_t iMomBin = fUtil -> GetMomentumBin(momentum);
  hMomBin -> Fill(iMomBin);
  return hMomBin;
}


//________________________________________________________
Bool_t AliTRDpidChecker::GetRefFigure(Int_t ifig)
{
  Bool_t FIRST = kTRUE;
  TLegend *leg = new TLegend(.7, .7, .98, .98);
  leg->SetBorderSize(1);
  TGraphErrors *g = 0x0;
  TAxis *ax = 0x0;
  TH1 *h1 = 0x0, *h=0x0;
  TH2 *h2 = 0x0;
  TObjArray *arr = 0x0;
  switch(ifig){
  case kEfficiency:
    arr = (TObjArray*)fGraph->At(ifig);
    if(!(g = (TGraphErrors*)arr->At(kLQ))) break;
    if(!g->GetN()) break;
    leg->SetHeader("PID Method");
    g->Draw("apl");
    ax = g->GetHistogram()->GetXaxis();
    ax->SetTitle("p [GeV/c]");
    ax->SetRangeUser(.6, 10.5);
    ax->SetMoreLogLabels();
    ax = g->GetHistogram()->GetYaxis();
    ax->SetTitle("#epsilon_{#pi} [%]");
    ax->SetRangeUser(1.e-3, 1.e-1);
    leg->AddEntry(g, "2D LQ", "pl");
    if(! (g = (TGraphErrors*)arr->At(kNN))) break;
    g->Draw("pl");
    leg->AddEntry(g, "NN", "pl");
    if(! (g = (TGraphErrors*)arr->At(kESD))) break;
    g->Draw("p");
    leg->AddEntry(g, "ESD", "pl");
    leg->Draw();
    gPad->SetLogy();
    gPad->SetLogx();
    gPad->SetGridy();
    gPad->SetGridx();
    return kTRUE;
  case kdEdx:
    // save 2.0 GeV projection as reference
    FIRST = kTRUE;
    if(!(h2 = (TH2F*)(fContainer->At(kdEdx)))) break;
    leg->SetHeader("Particle Species");
    for(Int_t is = AliPID::kSPECIES-1; is>=0; is--){
      Int_t bin = is*AliTRDCalPID::kNMom+4;
      h1 = h2->ProjectionY("px", bin, bin);
      if(!h1->GetEntries()) continue;
      h1->Scale(1./h1->Integral());
      h1->SetLineColor(AliTRDCalPID::GetPartColor(is));
      h = (TH1F*)h1->DrawClone(FIRST ? "c" : "samec");
      leg->AddEntry(h, Form("%s", AliTRDCalPID::GetPartName(is)), "l");
      FIRST = kFALSE;
    }
    if(FIRST) break;
    leg->Draw();
    gPad->SetLogy();
    gPad->SetLogx(0);
    gPad->SetGridy();
    gPad->SetGridx();
    return kTRUE;
  case kdEdxSlice:
    break;
  case kPH:
    // save 2.0 GeV projection as reference
    FIRST = kTRUE;
    if(!(h2 = (TH2F*)(fContainer->At(kPH)))) break;;
    leg->SetHeader("Particle Species");
    for(Int_t is=0; is<AliPID::kSPECIES; is++){
      Int_t bin = is*AliTRDCalPID::kNMom+4;
      h1 = h2->ProjectionY("py", bin, bin);
      if(!h1->GetEntries()) continue;
      h1->SetMarkerStyle(24);
      h1->SetMarkerColor(AliTRDCalPID::GetPartColor(is));
      h1->SetLineColor(AliTRDCalPID::GetPartColor(is));
      if(FIRST){
        h1->GetXaxis()->SetTitle("tb[1/100 ns^{-1}]");
        h1->GetYaxis()->SetTitle("<PH> [a.u.]");
      }
      h = (TH1F*)h1->DrawClone(FIRST ? "c" : "samec");
      leg->AddEntry(h, Form("%s", AliTRDCalPID::GetPartName(is)), "pl");
      FIRST = kFALSE;
    }
    if(FIRST) break;
    leg->Draw();
    gPad->SetLogy(0);
    gPad->SetLogx(0);
    gPad->SetGridy();
    gPad->SetGridx();
    return kTRUE;
  case kNClus:
    // save 2.0 GeV projection as reference
    FIRST = kTRUE;
    if(!(h2 = (TH2F*)(fContainer->At(kNClus)))) break;
    leg->SetHeader("Particle Species");
    for(Int_t is=0; is<AliPID::kSPECIES; is++){
      Int_t bin = is*AliTRDCalPID::kNMom+4;
      h1 = h2->ProjectionY("py", bin, bin);
      if(!h1->GetEntries()) continue;
      //h1->SetMarkerStyle(24);
      //h1->SetMarkerColor(AliTRDCalPID::GetPartColor(is));
      h1->SetLineColor(AliTRDCalPID::GetPartColor(is));
      if(FIRST) h1->GetXaxis()->SetTitle("N^{cl}/tracklet");
      h = (TH1F*)h1->DrawClone(FIRST ? "c" : "samec");
      leg->AddEntry(h, Form("%s", AliTRDCalPID::GetPartName(is)), "l");
      FIRST = kFALSE;
    }
    if(FIRST) break;
    leg->Draw();
    gPad->SetLogy();
    gPad->SetLogx(0);
    gPad->SetGridy();
    gPad->SetGridx();
    return kTRUE;
  case kMomentum:
  case kMomentumBin:
    break; 
  case kThresh:
    arr = (TObjArray*)fGraph->FindObject("Thresholds");
    if(!(g = (TGraphErrors*)arr->At(kLQ))) break;
    if(!g->GetN()) break;
    leg->SetHeader("PID Method");
    g->Draw("apl");
    ax = g->GetHistogram()->GetXaxis();
    ax->SetTitle("p [GeV/c]");
    ax->SetRangeUser(.6, 10.5);
    ax->SetMoreLogLabels();
    ax = g->GetHistogram()->GetYaxis();
    ax->SetTitle("threshold");
    ax->SetRangeUser(5.e-2, 1.);
    leg->AddEntry(g, "2D LQ", "pl");
    if(!(g = (TGraphErrors*)arr->At(kNN))) break;
    g->Draw("pl");
    leg->AddEntry(g, "NN", "pl");
    if(!(g = (TGraphErrors*)arr->At(kESD))) break;
    g->Draw("p");
    leg->AddEntry(g, "ESD", "pl");
    leg->Draw();
    gPad->SetLogx();
    gPad->SetGridy();
    gPad->SetGridx();
    return kTRUE;
  }
  AliInfo(Form("Reference plot [%d] missing result", ifig));
  return kFALSE;
}

//________________________________________________________________________
Bool_t AliTRDpidChecker::PostProcess()
{
  // Draw result to the screen
  // Called once at the end of the query

  if (!fContainer) {
    Printf("ERROR: list not available");
    return kFALSE;
  }
  if(!(fEfficiency = dynamic_cast<TObjArray *>(fContainer->At(kEfficiency)))){
    AliError("Efficiency container missing.");
    return 0x0;
  }
  if(!fGraph){ 
    fGraph = new TObjArray(2);
    fGraph->SetOwner();
  }
  EvaluatePionEfficiency(fEfficiency, fGraph, 0.9);
  fNRefFigures = 8;
  return kTRUE;
}

//________________________________________________________________________
void AliTRDpidChecker::EvaluatePionEfficiency(TObjArray *histoContainer, TObjArray *results, Float_t electron_efficiency)
{
  if(!histoContainer || !results) return;

  TGraphErrors *g = 0x0;
  fUtil->SetElectronEfficiency(electron_efficiency);

  // efficiency graphs
  TObjArray *eff = new TObjArray(3); eff->SetOwner(); eff->SetName("Efficiencies");
  results->AddAt(eff, 0);
  eff->AddAt(g = new TGraphErrors(), kLQ);
  g->SetLineColor(kBlue);
  g->SetMarkerColor(kBlue);
  g->SetMarkerStyle(7);
  eff->AddAt(g = new TGraphErrors(), kNN);
  g->SetLineColor(kGreen);
  g->SetMarkerColor(kGreen);
  g->SetMarkerStyle(7);
  eff -> AddAt(g = new TGraphErrors(), kESD);
  g->SetLineColor(kRed);
  g->SetMarkerColor(kRed);
  g->SetMarkerStyle(24);

  // Threshold graphs
  TObjArray *thres = new TObjArray(3); thres->SetOwner(); thres->SetName("Thresholds");
  results->AddAt(thres, 1);
  thres->AddAt(g = new TGraphErrors(), kLQ);
  g->SetLineColor(kBlue);
  g->SetMarkerColor(kBlue);
  g->SetMarkerStyle(7);
  thres->AddAt(g = new TGraphErrors(), kNN);
  g->SetLineColor(kGreen);
  g->SetMarkerColor(kGreen);
  g->SetMarkerStyle(7);
  thres -> AddAt(g = new TGraphErrors(), kESD);
  g->SetLineColor(kRed);
  g->SetMarkerColor(kRed);
  g->SetMarkerStyle(24);
  
  Float_t mom = 0.;
  TH1D *Histo1=0x0, *Histo2=0x0;

  // calculate the pion efficiencies and the errors for 90% electron efficiency (2-dim LQ)
  TH2F *hPIDLQ = (TH2F*)histoContainer->At(kLQ);
  for(Int_t iMom = 0; iMom < AliTRDCalPID::kNMom; iMom++){
    mom = AliTRDCalPID::GetMomentum(iMom);

    Histo1 = hPIDLQ -> ProjectionY("LQ_ele",AliTRDCalPID::kNMom*AliPID::kElectron+iMom+1,AliTRDCalPID::kNMom*AliPID::kElectron+iMom+1);
    Histo2 = hPIDLQ -> ProjectionY("LQ_pio",AliTRDCalPID::kNMom*AliPID::kPion+iMom+1,AliTRDCalPID::kNMom*AliPID::kPion+iMom+1);

    if(!fUtil->CalculatePionEffi(Histo1, Histo2)) continue;

    g = (TGraphErrors*)eff->At(kLQ);
    g->SetPoint(iMom, mom, fUtil->GetPionEfficiency());
    g->SetPointError(iMom, 0., fUtil->GetError());
    g = (TGraphErrors*)thres->At(kLQ);
    g->SetPoint(iMom, mom, fUtil->GetThreshold());
    g->SetPointError(iMom, 0., 0.);

    if(fDebugLevel>=2) Printf("Pion Efficiency for 2-dim LQ is : %f +/- %f\n\n", fUtil->GetPionEfficiency(), fUtil->GetError());
  }
  

  // calculate the pion efficiencies and the errors for 90% electron efficiency (NN)
  TH2F *hPIDNN = (TH2F*)histoContainer->At(kNN);
  for(Int_t iMom = 0; iMom < AliTRDCalPID::kNMom; iMom++){
    mom = AliTRDCalPID::GetMomentum(iMom);

    Histo1 = hPIDNN -> ProjectionY("NN_ele",AliTRDCalPID::kNMom*AliPID::kElectron+iMom+1,AliTRDCalPID::kNMom*AliPID::kElectron+iMom+1);
    Histo2 = hPIDNN -> ProjectionY("NN_pio",AliTRDCalPID::kNMom*AliPID::kPion+iMom+1,AliTRDCalPID::kNMom*AliPID::kPion+iMom+1);

    if(!fUtil -> CalculatePionEffi(Histo1, Histo2)) continue;

    g = (TGraphErrors*)eff->At(kNN);
    g->SetPoint(iMom, mom, fUtil->GetPionEfficiency());
    g->SetPointError(iMom, 0., fUtil->GetError());
    g = (TGraphErrors*)thres->At(3+kNN);
    g->SetPoint(iMom, mom, fUtil->GetThreshold());
    g->SetPointError(iMom, 0., 0.);

    if(fDebugLevel>=2) Printf("Pion Efficiency for NN is : %f +/- %f\n\n", fUtil->GetPionEfficiency(), fUtil->GetError());
  }


  // calculate the pion efficiencies and the errors for 90% electron efficiency (ESD)
  TH2F *hPIDESD = (TH2F*)histoContainer->At(kESD);
  for(Int_t iMom = 0; iMom < AliTRDCalPID::kNMom; iMom++){
    mom = AliTRDCalPID::GetMomentum(iMom);

    Histo1 = hPIDESD -> ProjectionY("NN_ele",AliTRDCalPID::kNMom*AliPID::kElectron+iMom+1,AliTRDCalPID::kNMom*AliPID::kElectron+iMom+1);
    Histo2 = hPIDESD -> ProjectionY("NN_pio",AliTRDCalPID::kNMom*AliPID::kPion+iMom+1,AliTRDCalPID::kNMom*AliPID::kPion+iMom+1);

    if(!fUtil->CalculatePionEffi(Histo1, Histo2)) continue;

    g = (TGraphErrors*)eff->At(kESD);
    g->SetPoint(iMom, mom, fUtil->GetPionEfficiency());
    g->SetPointError(iMom, 0., fUtil->GetError());
    g = (TGraphErrors*)thres->At(3+kESD);
    g->SetPoint(iMom, mom, fUtil->GetThreshold());
    g->SetPointError(iMom, 0., 0.);

    if(fDebugLevel>=2) Printf("Pion Efficiency for ESD is : %f +/- %f\n\n", fUtil->GetPionEfficiency(), fUtil->GetError());
  }

}

//________________________________________________________________________
void AliTRDpidChecker::Terminate(Option_t *) 
{
  // Draw result to the screen
  // Called once at the end of the query

  fContainer = dynamic_cast<TObjArray*>(GetOutputData(0));
  if (!fContainer) {
    Printf("ERROR: list not available");
    return;
  }
}


