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
  fContainer = new TObjArray();
  fContainer -> Expand(kGraphNN + 1);

  const Float_t epsilon = 1./(2*(AliTRDpidUtil::kBins-1));     // get nice histos with bin center at 0 and 1

  // histos of the electron probability of all 5 particle species and 11 momenta for the 2-dim LQ method 
  fContainer->AddAt(
    new TH2F("PID_LQ", "", 
      xBins, -0.5, xBins - 0.5,
      AliTRDpidUtil::kBins, 0.-epsilon, 1.+epsilon)
  ,kLQlikelihood);


  // histos of the electron probability of all 5 particle species and 11 momenta for the neural network method
  fContainer->AddAt(
    new TH2F("PID_NN", "", 
      xBins, -0.5, xBins - 0.5,
      AliTRDpidUtil::kBins, 0.-epsilon, 1.+epsilon)
  ,kNNlikelihood);

  // histos of the dE/dx distribution for all 5 particle species and 11 momenta 
  fContainer->AddAt(
    new TH2F("dEdx", "", 
      xBins, -0.5, xBins - 0.5,
      200, 0, 10000)
    ,kdEdx);

  // histos of the dE/dx distribution for all 5 particle species and 11 momenta 
  fContainer->AddAt(
    new TH2F("dEdxSlice", "", 
      xBins*AliTRDReconstructor::kLQslices, -0.5, xBins*AliTRDReconstructor::kLQslices - 0.5,
      200, 0, 5000)
    ,kdEdxSlice);

  // histos of the pulse height distribution for all 5 particle species and 11 momenta 
  fContainer->AddAt(
    new TProfile2D("PH", "", 
      xBins, -0.5, xBins - 0.5,
      AliTRDtrackerV1::GetNTimeBins(), -0.5, AliTRDtrackerV1::GetNTimeBins() - 0.5)
    ,kPH);

  // histos of the number of clusters distribution for all 5 particle species and 11 momenta 
  fContainer->AddAt(
    new TH2F("NClus", "", 
      xBins, -0.5, xBins - 0.5,
      AliTRDtrackerV1::GetNTimeBins(), -0.5, AliTRDtrackerV1::GetNTimeBins() - 0.5)
    ,kNClus);


  // momentum distributions - absolute and in momentum bins
  fContainer->AddAt(new TH1F("hMom", "momentum distribution", 100, 0., 12.),kMomentum);
  fContainer->AddAt(new TH1F("hMomBin", "momentum distribution in momentum bins", AliTRDCalPID::kNMom, 0.5, 11.5),kMomentumBin);


  // TGraph of the pion efficiencies

  TGraphErrors *gEffisLQ = 0x0;
  TGraphErrors *gEffisNN = 0x0;

  fContainer->AddAt(gEffisLQ = new TGraphErrors(), kGraphLQ);
  gEffisLQ->SetLineColor(kBlue);
  gEffisLQ->SetMarkerColor(kBlue);
  gEffisLQ->SetMarkerStyle(29);

  fContainer -> AddAt(gEffisNN = new TGraphErrors(),kGraphNN);
  gEffisNN->SetLineColor(kRed);
  gEffisNN->SetMarkerColor(kRed);
  gEffisNN->SetMarkerStyle(29);

  gEffisLQ -> SetNameTitle("gEffisLQErr", "Efficiencies and Errors of the 2-dim LQ method");
  gEffisNN -> SetNameTitle("gEffisNNErr", "Efficiencies and Errors of the NN method");

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
  
  TH2F *hPIDLQ;
  if(!(hPIDLQ = dynamic_cast<TH2F *>(fContainer->At(kLQlikelihood)))){
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


  fReconstructor -> SetOption("!nn");
  cTrack.CookPID();

  Int_t iMomBin = -1;
  iMomBin = fUtil->GetMomentumBin(momentum);

  if(momentum < 0.4) return 0x0;;

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
  
  TH2F *hPIDNN;
  if(!(hPIDNN = dynamic_cast<TH2F *>(fContainer->At(kNNlikelihood)))){
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


  fReconstructor -> SetOption("nn");
  cTrack.CookPID();

  Int_t iMomBin = -1;
  iMomBin = fUtil -> GetMomentumBin(momentum);

  if(momentum < 0.4) return 0x0;;

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

  Int_t iMomBin = -1;
  iMomBin = fUtil -> GetMomentumBin(momentum);

  if(momentum < 0.4) return 0x0;;

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

  Int_t iMomBin = -1;
  iMomBin = fUtil -> GetMomentumBin(momentum);

  if(momentum < 0.4) return 0x0;

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
    for(Int_t iChamb = 0; iChamb < AliTRDgeometry::kNlayer; iChamb++)
      for(Int_t iSlice = 0; iSlice < AliTRDReconstructor::kLQslices; iSlice++)
	hdEdxSlice -> Fill(AliPID::kElectron * AliTRDCalPID::kNMom * AliTRDReconstructor::kLQslices + iMomBin * AliTRDReconstructor::kLQslices + iSlice,
			   dEdxSlice[iChamb][iSlice]);
    break;
  case kMuonPlus:
  case kMuonMinus:
    for(Int_t iChamb = 0; iChamb < AliTRDgeometry::kNlayer; iChamb++)
      for(Int_t iSlice = 0; iSlice < AliTRDReconstructor::kLQslices; iSlice++)
	hdEdxSlice -> Fill(AliPID::kMuon * AliTRDCalPID::kNMom * AliTRDReconstructor::kLQslices + iMomBin * AliTRDReconstructor::kLQslices + iSlice,
			   dEdxSlice[iChamb][iSlice]);
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

  Int_t iMomBin = -1;
  iMomBin = fUtil -> GetMomentumBin(momentum);

  if(momentum < 0.4) return 0x0;;

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

  Int_t iMomBin = -1;
  iMomBin = fUtil -> GetMomentumBin(momentum);

  if(momentum < 0.4) return 0x0;;

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

  if(momentum < 0.4) return 0x0;

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

  Int_t iMomBin = -1;
  iMomBin = fUtil -> GetMomentumBin(momentum);

  if(momentum < 0.4) return 0x0;

  hMomBin -> Fill(iMomBin);
  return hMomBin;
}


//________________________________________________________
void AliTRDpidChecker::GetRefFigure(Int_t ifig)
{
  Bool_t FIRST = kTRUE;
  TGraphErrors *g = 0x0;
  TH1 *h1 = 0x0;
  TH2 *h2 = 0x0;
  switch(ifig){
  case 0:
    g = (TGraphErrors*)fContainer->At(kGraphStart);
    g->Draw("apl");
    g->GetHistogram()->GetXaxis()->SetTitle("p [GeV/c]");
    g->GetHistogram()->GetXaxis()->SetRangeUser(.6, 10.5);
    g->GetHistogram()->GetYaxis()->SetTitle("#epsilon_{#pi} [%]");
    ((TGraphErrors*)fContainer->At(kGraphStart+1))->Draw("pl");
    gPad->SetLogy();
    gPad->SetLogx();
    gPad->SetGridy();
    gPad->SetGridx();
    break;
  case 1:
    // save 2.0 GeV projection as reference
    FIRST = kTRUE;
    h2 = (TH2F*)(fContainer->At(kdEdx));
    for(Int_t is = AliPID::kSPECIES-1; is>=0; is--){
      Int_t bin = is*AliTRDCalPID::kNMom+4;
      h1 = h2->ProjectionY("px", bin, bin);
      if(!h1->GetEntries()) continue;
      h1->Scale(1./h1->Integral());
      h1->SetLineColor(AliTRDCalPID::GetPartColor(is));
      h1->DrawClone(FIRST ? "c" : "samec");
      FIRST = kFALSE;
    }
    gPad->SetLogy();
    gPad->SetLogx(0);
    gPad->SetGridy();
    gPad->SetGridx();
    break;
  case 2:
    // save 2.0 GeV projection as reference
    FIRST = kTRUE;
    h2 = (TH2F*)(fContainer->At(kPH));
    for(Int_t is=0; is<AliPID::kSPECIES; is++){
      Int_t bin = is*AliTRDCalPID::kNMom+4;
      h1 = h2->ProjectionY("py", bin, bin);
      if(!h1->GetEntries()) continue;
      h1->SetMarkerStyle(24);
      h1->SetMarkerColor(AliTRDCalPID::GetPartColor(is));
      h1->SetLineColor(AliTRDCalPID::GetPartColor(is));
      h1->DrawClone(FIRST ? "e2" : "same e2");
      FIRST = kFALSE;
    }
    gPad->SetLogy(0);
    gPad->SetLogx(0);
    gPad->SetGridy();
    gPad->SetGridx();
    break;
  case 3:
    // save 2.0 GeV projection as reference
    FIRST = kTRUE;
    h2 = (TH2F*)(fContainer->At(kNClus));
    for(Int_t is=0; is<AliPID::kSPECIES; is++){
      Int_t bin = is*AliTRDCalPID::kNMom+4;
      h1 = h2->ProjectionY("py", bin, bin);
      if(!h1->GetEntries()) continue;
      h1->SetMarkerStyle(24);
      h1->SetMarkerColor(AliTRDCalPID::GetPartColor(is));
      h1->SetLineColor(AliTRDCalPID::GetPartColor(is));
      h1->DrawClone(FIRST ? "e2" : "same e2");
      FIRST = kFALSE;
    }
    gPad->SetLogy(0);
    gPad->SetLogx(0);
    gPad->SetGridy();
    gPad->SetGridx();
    break;
  }
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
//   return kTRUE; // testing protection


  // container for the pion efficiencies and the errors
  Double_t  PionEffiLQ[AliTRDCalPID::kNMom], 
            PionEffiErrorLQ[AliTRDCalPID::kNMom], 
            EleEffiLQ[AliTRDCalPID::kNMom],
            ThresholdLQ[AliTRDCalPID::kNMom];

  Double_t  PionEffiNN[AliTRDCalPID::kNMom],
            PionEffiErrorNN[AliTRDCalPID::kNMom],
            EleEffiNN[AliTRDCalPID::kNMom],
            ThresholdNN[AliTRDCalPID::kNMom];

  Float_t mom = 0.;

  TH1D *Histo1=0x0, *Histo2=0x0;

  TH2F *hPIDLQ=0x0, *hPIDNN=0x0;
  hPIDLQ = (TH2F*)fContainer->At(kLQlikelihood);
  hPIDNN = (TH2F*)fContainer->At(kNNlikelihood);

  // calculate the pion efficiencies and the errors for 90% electron efficiency (2-dim LQ)
  for(Int_t iMom = 0; iMom < AliTRDCalPID::kNMom; iMom++){

    AliTRDpidUtil *util = new AliTRDpidUtil();
    mom = AliTRDCalPID::GetMomentum(iMom);

    Histo1 = hPIDLQ -> ProjectionY("LQ_ele",AliTRDCalPID::kNMom*AliPID::kElectron+iMom+1,AliTRDCalPID::kNMom*AliPID::kElectron+iMom+1);
    Histo2 = hPIDLQ -> ProjectionY("LQ_pio",AliTRDCalPID::kNMom*AliPID::kPion+iMom+1,AliTRDCalPID::kNMom*AliPID::kPion+iMom+1);

    util -> CalculatePionEffi(Histo1, Histo2);

    PionEffiLQ[iMom] = util -> GetPionEfficiency();
    PionEffiErrorLQ[iMom] = util -> GetError();
    EleEffiLQ[iMom] = util -> GetCalcElectronEfficiency();
    ThresholdLQ[iMom] = util -> GetThreshold();

    if(fDebugLevel>=1) Printf("Pion Efficiency for 2-dim LQ is : %f +/- %f\n\n", PionEffiLQ[iMom], PionEffiErrorLQ[iMom]);
    util -> Delete();
  }
  


  // calculate the pion efficiencies and the errors for 90% electron efficiency (NN)
  for(Int_t iMom = 0; iMom < AliTRDCalPID::kNMom; iMom++){

    AliTRDpidUtil *util = new AliTRDpidUtil();
    mom = AliTRDCalPID::GetMomentum(iMom);

    Histo1 = hPIDNN -> ProjectionY("NN_ele",AliTRDCalPID::kNMom*AliPID::kElectron+iMom+1,AliTRDCalPID::kNMom*AliPID::kElectron+iMom+1);
    Histo2 = hPIDNN -> ProjectionY("NN_pio",AliTRDCalPID::kNMom*AliPID::kPion+iMom+1,AliTRDCalPID::kNMom*AliPID::kPion+iMom+1);

    util -> CalculatePionEffi(Histo1, Histo2);

    PionEffiNN[iMom] = util -> GetPionEfficiency();
    PionEffiErrorNN[iMom] = util -> GetError();
    EleEffiNN[iMom] = util -> GetCalcElectronEfficiency();
    ThresholdNN[iMom] = util -> GetThreshold();

    if(fDebugLevel>=1) Printf("Pion Efficiency for NN is : %f +/- %f\n\n", PionEffiNN[iMom], PionEffiErrorNN[iMom]);

    util -> Delete();
  }
  

  // create TGraph to plot the pion efficiencies
  TGraphErrors *gEffisLQ=0x0, *gEffisNN=0x0;
  gEffisLQ = (TGraphErrors*)fContainer->At(kGraphLQ);
  gEffisNN = (TGraphErrors*)fContainer->At(kGraphNN);


  for(Int_t iBin = 0; iBin < AliTRDCalPID::kNMom; iBin++){

    Float_t momentum = AliTRDCalPID::GetMomentum(iBin);
    gEffisLQ->SetPoint(iBin, momentum, PionEffiLQ[iBin]);
    gEffisLQ->SetPointError(iBin, 0., PionEffiErrorLQ[iBin]);

    gEffisNN->SetPoint(iBin, momentum, PionEffiNN[iBin]);
    gEffisNN->SetPointError(iBin, 0., PionEffiErrorNN[iBin]);
  }


  fNRefFigures = 4/*1*/;
  return kTRUE; // testing protection
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


