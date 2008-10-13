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

// #include "AliPID.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliTrackReference.h"

#include "AliAnalysisTask.h"
// #include "AliAnalysisManager.h"

#include "AliTRDtrackerV1.h"
#include "AliTRDtrackV1.h"
#include "AliTRDcluster.h"
#include "AliTRDReconstructor.h"
#include "AliCDBManager.h"
// #include "../Cal/AliTRDCalPID.h"
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
{
  //
  // Default constructor
  //

  fReconstructor = new AliTRDReconstructor();
  fReconstructor->SetRecoParam(AliTRDrecoParam::GetLowFluxParam());
}


//________________________________________________________________________
AliTRDpidChecker::~AliTRDpidChecker() 
{
  if(fReconstructor) delete fReconstructor;
}


//________________________________________________________________________
void AliTRDpidChecker::CreateOutputObjects()
{
  // Create histograms
  // Called once

  OpenFile(0, "RECREATE");
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

  // histos of the pulse height distribution for all 5 particle species and 11 momenta 
  fContainer->AddAt(
    new TProfile2D("PH", "", 
      xBins, -0.5, xBins - 0.5,
      AliTRDtrackerV1::GetNTimeBins(), -0.5, AliTRDtrackerV1::GetNTimeBins() - 0.5)
    ,kPH);


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

}

//________________________________________________________________________
void AliTRDpidChecker::Exec(Option_t *) 
{
  // Main loop
  // Called for each event


//   if(!AliTracker::GetFieldMap()){
//     AliMagFMaps* field = new AliMagFMaps("Maps","Maps", 2, 1., 10., AliMagFMaps::k5kG);
//     AliTracker::SetFieldMap(field, kTRUE);
//   }

  TH2F *hPIDLQ;
  TH2F *hPIDNN;
  TH2F *hdEdx;
  TProfile2D *hPH;

  hPIDLQ = (TH2F*)fContainer->At(kLQlikelihood);
  hPIDNN = (TH2F*)fContainer->At(kNNlikelihood);
  hdEdx  = (TH2F*)fContainer->At(kdEdx);
  hPH    = (TProfile2D*)fContainer->At(kPH);
  
  TH1F *hMom    = (TH1F*)fContainer->At(kMomentum);	
  TH1F *hMomBin = (TH1F*)fContainer->At(kMomentumBin);	
  
  Int_t labelsacc[10000]; 
  memset(labelsacc, 0, sizeof(Int_t) * 10000);
  
  Float_t mom;
  ULong_t status;
  Int_t nTRD = 0;
  Float_t *fdEdx;       

  AliTRDtrackInfo     *track = 0x0;
  AliTRDtrackV1    *TRDtrack = 0x0;
  AliTrackReference     *ref = 0x0;
  AliExternalTrackParam *esd = 0x0;

  AliTRDseedV1 *TRDtracklet[AliTRDgeometry::kNlayer];
  for(Int_t iChamb = 0; iChamb < AliTRDgeometry::kNlayer; iChamb++)
    TRDtracklet[iChamb] = 0x0;

  AliTRDcluster *TRDcluster = 0x0;

  AliTRDpidUtil *util = new AliTRDpidUtil();
  for(Int_t itrk=0; itrk<fTracks->GetEntriesFast(); itrk++){
    track = (AliTRDtrackInfo*)fTracks->UncheckedAt(itrk);
    if(!track->HasESDtrack()) continue;
    status = track->GetStatus();
    if(!(status&AliESDtrack::kTPCout)) continue;

    if(!(TRDtrack = track->GetTRDtrack())) continue; 
    //&&(track->GetNumberOfClustersRefit()

    // use only tracks that hit 6 chambers
    if(!(TRDtrack->GetNumberOfTracklets() == AliTRDgeometry::kNlayer)) continue;
     
    ref = track->GetTrackRef(0);
    esd = track->GetOuterParam();
    mom = ref ? ref->P(): esd->P();

    labelsacc[nTRD] = track->GetLabel();
    nTRD++;
      

    // set the 11 momentum bins
    Int_t iMomBin = -1;
    iMomBin = util -> GetMomentumBin(mom);
    if(fDebugLevel>=4) Printf("MomBin[%d] MomTot[%f]", iMomBin, mom);


    // fill momentum histo to have the momentum distribution
    hMom -> Fill(mom);
    hMomBin -> Fill(iMomBin);


    // set the reconstructor
    TRDtrack -> SetReconstructor(fReconstructor);

    
    // if no monte carlo data available -> use TRDpid
    if(!HasMCdata()){
      fReconstructor -> SetOption("nn");
      TRDtrack -> CookPID();
      if(TRDtrack -> GetPID(0) > TRDtrack -> GetPID(1) + TRDtrack -> GetPID(2)  + TRDtrack -> GetPID(3) + TRDtrack -> GetPID(4)){
	track -> SetPDG(kElectron);
      }
      else if(TRDtrack -> GetPID(4) > TRDtrack -> GetPID(2)  && TRDtrack -> GetPID(4) > TRDtrack -> GetPID(3)  && TRDtrack -> GetPID(4) > TRDtrack -> GetPID(1)){
	track -> SetPDG(kProton);
      }
      else if(TRDtrack -> GetPID(3) > TRDtrack -> GetPID(1)  && TRDtrack -> GetPID(3) > TRDtrack -> GetPID(2)){
	track -> SetPDG(kKPlus);
      }
      else if(TRDtrack -> GetPID(1) > TRDtrack -> GetPID(2)){
	track -> SetPDG(kMuonPlus);
      }
      else{
	track -> SetPDG(kPiPlus);
      }
    }


    // calculate the probabilities for electron probability using 2-dim LQ, the deposited charge per chamber and the pulse height spectra and fill histograms
    fReconstructor -> SetOption("!nn");
    TRDtrack -> CookPID();

    if(fDebugLevel>=4) Printf("PIDmethod[%d] Slices[%d] PDG[%d] LQLike[%f]", fReconstructor->GetPIDMethod(), fReconstructor->GetNdEdxSlices(), track->GetPDG(), TRDtrack -> GetPID(0));


    Float_t SumdEdx[AliTRDgeometry::kNlayer];
    for(Int_t iChamb = 0; iChamb < AliTRDgeometry::kNlayer; iChamb++){
      TRDtracklet[iChamb] = TRDtrack -> GetTracklet(iChamb);
      SumdEdx[iChamb] = 0.;
      fdEdx = TRDtracklet[iChamb] -> GetdEdx();
      SumdEdx[iChamb] += fdEdx[0] + fdEdx[1] + fdEdx[2]; 
    }

    switch(track->GetPDG()){
    case kElectron:
    case kPositron:
      hPIDLQ -> Fill(AliPID::kElectron * AliTRDCalPID::kNMom + iMomBin, TRDtrack -> GetPID(0));
      for(Int_t iChamb = 0; iChamb < AliTRDgeometry::kNlayer; iChamb++){
        hdEdx -> Fill(AliPID::kElectron * AliTRDCalPID::kNMom + iMomBin, SumdEdx[iChamb]);
        for(Int_t iClus = 0; iClus < AliTRDtrackerV1::GetNTimeBins(); iClus++){
          if(!(TRDcluster = (AliTRDcluster*)TRDtracklet[iChamb] -> GetClusters(iClus)))
            continue;
          hPH -> Fill(AliPID::kElectron * AliTRDCalPID::kNMom + iMomBin, TRDcluster -> GetLocalTimeBin(), TRDtracklet[iChamb] -> GetdQdl(iClus));
        }
      }
      break;
    case kMuonPlus:
    case kMuonMinus:
      hPIDLQ -> Fill(AliPID::kMuon * AliTRDCalPID::kNMom + iMomBin, TRDtrack -> GetPID(0));
      for(Int_t iChamb = 0; iChamb < AliTRDgeometry::kNlayer; iChamb++){
        hdEdx -> Fill(AliPID::kMuon * AliTRDCalPID::kNMom + iMomBin, SumdEdx[iChamb]);
        for(Int_t iClus = 0; iClus < AliTRDtrackerV1::GetNTimeBins(); iClus++){
          if(!(TRDcluster = (AliTRDcluster*)TRDtracklet[iChamb] -> GetClusters(iClus)))
            continue;
          hPH -> Fill(AliPID::kMuon * AliTRDCalPID::kNMom + iMomBin, TRDcluster -> GetLocalTimeBin(), TRDtracklet[iChamb] -> GetdQdl(iClus));
        }
      }
      break;
    case kPiPlus:
    case kPiMinus:
      hPIDLQ -> Fill(AliPID::kPion * AliTRDCalPID::kNMom + iMomBin, TRDtrack -> GetPID(0));
      for(Int_t iChamb = 0; iChamb < AliTRDgeometry::kNlayer; iChamb++){
        hdEdx -> Fill(AliPID::kPion * AliTRDCalPID::kNMom + iMomBin, SumdEdx[iChamb]);
        for(Int_t iClus = 0; iClus < AliTRDtrackerV1::GetNTimeBins(); iClus++){
          if(!(TRDcluster = (AliTRDcluster*)TRDtracklet[iChamb] -> GetClusters(iClus)))
            continue;
          hPH -> Fill(AliPID::kPion * AliTRDCalPID::kNMom + iMomBin, TRDcluster -> GetLocalTimeBin(), TRDtracklet[iChamb] -> GetdQdl(iClus));
        }
      }
      break;
    case kKPlus:
    case kKMinus:
      hPIDLQ -> Fill(AliPID::kKaon * AliTRDCalPID::kNMom + iMomBin, TRDtrack -> GetPID(0));
      for(Int_t iChamb = 0; iChamb < AliTRDgeometry::kNlayer; iChamb++){
        hdEdx -> Fill(AliPID::kKaon * AliTRDCalPID::kNMom + iMomBin, SumdEdx[iChamb]);
        for(Int_t iClus = 0; iClus < AliTRDtrackerV1::GetNTimeBins(); iClus++){
          if(!(TRDcluster = (AliTRDcluster*)TRDtracklet[iChamb] -> GetClusters(iClus)))
            continue;
          hPH -> Fill(AliPID::kKaon * AliTRDCalPID::kNMom + iMomBin, TRDcluster -> GetLocalTimeBin(), TRDtracklet[iChamb] -> GetdQdl(iClus));
        }
      }
      break;
    case kProton:
    case kProtonBar:
      hPIDLQ -> Fill(AliPID::kProton * AliTRDCalPID::kNMom + iMomBin, TRDtrack -> GetPID(0));
      for(Int_t iChamb = 0; iChamb < AliTRDgeometry::kNlayer; iChamb++){
        hdEdx -> Fill(AliPID::kProton * AliTRDCalPID::kNMom + iMomBin, SumdEdx[iChamb]);
        for(Int_t iClus = 0; iClus < AliTRDtrackerV1::GetNTimeBins(); iClus++){
          if(!(TRDcluster = (AliTRDcluster*)TRDtracklet[iChamb] -> GetClusters(iClus)))
            continue;
          hPH -> Fill(AliPID::kProton * AliTRDCalPID::kNMom + iMomBin, TRDcluster -> GetLocalTimeBin(), TRDtracklet[iChamb] -> GetdQdl(iClus));
        }
      }
      break;
    }


    // calculate the probabilities and fill histograms for electrons using NN
    fReconstructor -> SetOption("nn");
    TRDtrack->CookPID();


    if(fDebugLevel>=4) Printf("PIDmethod[%d] Slices[%d] PDG[%d] NNLike[%f]", fReconstructor->GetPIDMethod(), fReconstructor->GetNdEdxSlices(), track->GetPDG(), TRDtrack -> GetPID(0));


    switch(track->GetPDG()){
    case kElectron:
    case kPositron:
      hPIDNN -> Fill(AliPID::kElectron * AliTRDCalPID::kNMom + iMomBin, TRDtrack -> GetPID(0));
      break;
    case kMuonPlus:
    case kMuonMinus:
      hPIDNN -> Fill(AliPID::kMuon * AliTRDCalPID::kNMom + iMomBin, TRDtrack -> GetPID(0));
      break;
    case kPiPlus:
    case kPiMinus:
      hPIDNN -> Fill(AliPID::kPion * AliTRDCalPID::kNMom + iMomBin, TRDtrack -> GetPID(0));
      break;
    case kKPlus:
    case kKMinus:
      hPIDNN -> Fill(AliPID::kKaon * AliTRDCalPID::kNMom + iMomBin, TRDtrack -> GetPID(0));
      break;
    case kProton:
    case kProtonBar:
      hPIDNN -> Fill(AliPID::kProton * AliTRDCalPID::kNMom + iMomBin, TRDtrack -> GetPID(0));
      break;
    }

  }

  util -> Delete();
  PostData(0, fContainer);
}


//________________________________________________________
void AliTRDpidChecker::GetRefFigure(Int_t ifig)
{
  Bool_t FIRST = kTRUE;
  TH1 *h1 = 0x0;
  TH2 *h2 = 0x0;
  switch(ifig){
  case 0:
    ((TGraphErrors*)fContainer->At(kGraphStart))->Draw("apl");
    ((TGraphErrors*)fContainer->At(kGraphStart+1))->Draw("pl");
    gPad->SetLogy();
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


  fNRefFigures = 3/*1*/;
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


