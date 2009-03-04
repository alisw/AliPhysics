
/**************************************************************************
* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
*                                                                        *
* Author: The ALICE Off-line Project.                                    *
* Contributors are mentioned in the code where appropriate.              *
*                                                                        *
* Permission to use, copy, modify and distribute this software and its   *
* documentation strictly for non-commercial purposes is hereby granted   *
* without fee, provided that the above copyright notice appears in all   *
* copies and that both the copyright notice and this permission notice   *
* appear in the supporting documentation. The authors make no claims     *
* about the suitability of this software for any purpose. It is          *
* provided "as is" without express or implied warranty.                  *
**************************************************************************/

// macro for very simple analysis



#include "Riostream.h"
#include "TChain.h"
#include "TTree.h"
#include "TProfile2D.h"
#include "TH2I.h"
#include "TH1F.h"
#include "TList.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TObject.h"
#include "TFile.h"
#include "TObjArray.h"

#include "AliTRDrecoTask.h"
#include "AliAnalysisManager.h"

#include "AliESDInputHandler.h"
#include "AliTRDtrackV1.h"
#include "AliTRDseedV1.h"
#include "AliTRDcluster.h"
#include "AliTRDtrackInfo/AliTRDtrackInfo.h"
#include "AliTRDcalibDB.h"

#include "AliTRDCalibraFillHisto.h"

#include "AliLog.h"

#include "AliTRDcalibration.h"


ClassImp(AliTRDcalibration)

//________________________________________________________________________
AliTRDcalibration::AliTRDcalibration() 
  :AliTRDrecoTask("Calibration", "Calibration on tracks")
  ,fTrackInfo(0)
  ,ftrdTrack(0)
  ,fcl(0)
  ,fTRDCalibraFillHisto(0)
  ,fNbTRDTrackUsed(0)
  ,fNbTimeBin(0x0)
  ,fNbClusters(0)
  ,fPHSum(0)
  ,fCHSum(0)
  ,flow(0)
  ,fhigh(30)
  ,fNbTimeBins(30)
  ,ffillZero(kFALSE)
{
  // Constructor
}  

//________________________________________________________________________
void AliTRDcalibration::CreateOutputObjects() 
{
  OpenFile(0, "RECREATE");

  // Number of time bins
  AliTRDcalibDB *cal = AliTRDcalibDB::Instance();
  fNbTimeBins = cal->GetNumberOfTimeBins();
  
  // instance calibration
  fTRDCalibraFillHisto = AliTRDCalibraFillHisto::Instance();
  fTRDCalibraFillHisto->SetHisto2d(); // choose to use histograms
  fTRDCalibraFillHisto->SetCH2dOn();  // choose to calibrate the gain
  fTRDCalibraFillHisto->SetPH2dOn();  // choose to calibrate the drift velocity
  fTRDCalibraFillHisto->SetPRF2dOn(); // choose to look at the PRF
  fTRDCalibraFillHisto->Init2Dhistos(); // initialise the histos
  fTRDCalibraFillHisto->SetNumberClusters(flow); // At least 11 clusters
  fTRDCalibraFillHisto->SetNumberClustersf(fhigh); // At least 11 clusters
  fTRDCalibraFillHisto->SetFillWithZero(ffillZero); // Fill zeros

  fTRDCalibraFillHisto->SetDebugLevel(fDebugLevel); //debug stuff

  fContainer = new TObjArray();
  fContainer->Add(fTRDCalibraFillHisto->GetCH2d()); //TH2I
  fContainer->Add(fTRDCalibraFillHisto->GetPH2d()); //TProfile2D
  fContainer->Add(fTRDCalibraFillHisto->GetPRF2d()); //TProfile2D

  //printf("create output objects 1\n");

  fNbTRDTrackUsed = new TH1F("TRDTrackUsed","TRDTrackUsed",50,0,50);
  fNbTRDTrackUsed->Sumw2();
  //
  fNbTimeBin = new TH1F("TimeBin","TimeBin",35,0,35);
  fNbTimeBin->Sumw2();
  //
  fNbClusters = new TH1F("NbClusters","",35,0,35);
  fNbClusters->Sumw2();
  //
  fPHSum = new TProfile2D("PH2dSum","Nz0Nrphi0"
      ,30,-0.05,(Double_t)fNbTimeBins/10.0-0.05
      ,540,0,540);
  fPHSum->SetYTitle("Det/pad groups");
  fPHSum->SetXTitle("time [#mus]");
  fPHSum->SetZTitle("<PH> [a.u.]");
  fPHSum->SetStats(0);
  //
  fCHSum = new TH2I("CH2dSum","Nz0Nrphi0",100,0,300,540,0,540);
  fCHSum->SetYTitle("Det/pad groups");
  fCHSum->SetXTitle("charge deposit [a.u]");
  fCHSum->SetZTitle("counts");
  fCHSum->SetStats(0);
  fCHSum->Sumw2();

  fContainer->Add(fNbTRDTrackUsed);
  fContainer->Add(fNbTimeBin);
  fContainer->Add(fNbClusters);
  fContainer->Add(fPHSum);
  fContainer->Add(fCHSum);
  //printf("create output objects 2\n");
  
}

//________________________________________________________________________
void AliTRDcalibration::Exec(Option_t *) 
{
  
  Int_t nbTrdTracksUsed = 0;
  
  for(Int_t itrk=0; itrk < fTracks->GetEntriesFast(); itrk++){
    
    fTrackInfo = (AliTRDtrackInfo*)fTracks->UncheckedAt(itrk);
    ftrdTrack = fTrackInfo->GetTrack();
    if(!ftrdTrack) continue;
    nbTrdTracksUsed++;
    fTRDCalibraFillHisto->UpdateHistogramsV1(ftrdTrack);
    
    const AliTRDseedV1 *tracklet = 0x0;
    for(Int_t itr = 0; itr < AliTRDgeometry::kNlayer; itr++){
      if(!(tracklet = ftrdTrack->GetTracklet(itr))) continue;
      if(!tracklet->IsOK()) continue;
      Int_t nbclusters = 0;
      // For PH
      Double_t phtb[AliTRDseedV1::kNtb];
      memset(phtb, 0, AliTRDseedV1::kNtb*sizeof(Double_t));

      // For CH
      Double_t sum = 0.0;
      // normalisation
      Float_t normalisation = 6.67;
      Int_t detector = 0;
      for(int ic=0; ic<AliTRDseedV1::kNclusters; ic++){
        if(!(fcl = tracklet->GetClusters(ic))) continue;
        nbclusters++;
        Int_t time = fcl->GetPadTime();
        Float_t ch =  tracklet->GetdQdl(ic);
        detector = fcl->GetDetector();	  
        if((time>-1) && (time<fNbTimeBins)) phtb[time]=ch/normalisation;
        sum += ch/normalisation;
        //printf("time %d\n",time);
        fNbTimeBin->Fill(time);
      }
      fNbClusters->Fill(nbclusters);
      if((nbclusters > flow) && (nbclusters < fhigh)){
        fCHSum->Fill(sum/20.0,0.0);
        for(int ic=0; ic<fNbTimeBins; ic++){
          if(ffillZero) fPHSum->Fill((Double_t)ic/10.0,0.0,(Double_t)phtb[ic]);
          else {
            if(phtb[ic] > 0.0) fPHSum->Fill((Double_t)ic/10.0,0.0,(Double_t)phtb[ic]);
          }
        }
      }
    }
  }
  
  //Fill Histos
  fNbTRDTrackUsed->Fill(nbTrdTracksUsed);
  
  
  // Post output data
  PostData(0, fContainer);
  
}      
//________________________________________________________________________
void AliTRDcalibration::Terminate(Option_t *) 
{
  // Draw result to the screen
  // Called once at the end of the query

  //printf("terminate\n");

  if(fTRDCalibraFillHisto) fTRDCalibraFillHisto->DestroyDebugStreamer();

  /*
  fContainer = dynamic_cast<TList*> (GetOutputData(0));
  if (!fContainer) {
    Printf("ERROR: fList not available");
    return;
  }
  */
}
