/* $Id$ */

// This class contains a number of histograms for diagnostics of a TPC
// read out chamber from the reconstructed clusters.
//
// TODO:
//  
//
//

#include "AliTPCClusterHistograms.h"

#include <TStyle.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TH2F.h>
#include <TProfile2D.h>
#include <TObjArray.h>
#include <TLatex.h>
#include <TTimeStamp.h>
#include <TRandom.h>

#include <AliTPCclusterMI.h>
#include <AliTPCseed.h>

#include <AliLog.h>


//____________________________________________________________________
ClassImp(AliTPCClusterHistograms)

//____________________________________________________________________
AliTPCClusterHistograms::AliTPCClusterHistograms() 
  : TNamed(),
  fhQmaxVsRow(0),          
  fhQtotVsRow(0),          
  fhQtotProfileVsRow(0),   
  fhQmaxProfileVsRow(0),
  fhNClustersYVsRow(0),  
  fhNClustersZVsRow(0),
  fhSigmaYVsRow(0),        
  fhSigmaZVsRow(0),          			
  fhQmaxProfileYVsRow(0), 
  fhQtotProfileYVsRow(0),
  fhSigmaYProfileYVsRow(0),
  fhSigmaZProfileYVsRow(0),
  fhQmaxProfileZVsRow(0), 
  fhQtotProfileZVsRow(0),
  fhSigmaYProfileZVsRow(0),
  fhSigmaZProfileZVsRow(0),
  fhMeanQtotVsTime(0),  
  fhQtotVsTime(0),
  fhMeanNClustersVsTime(0),    
  fhNClustersVsTime(0),        
  fhTrackQtotPerCluster(0),
  fhTrackQtotPerClusterVsPhi(0),
  fhTrackQtotPerClusterVsTheta(0),
  fhTrackMeanQtotPerClusterVsPhi(0),
  fhTrackMeanQtotPerClusterVsTheta(0),
  fhMeanNTracksVsTime(),
  fhNEventsVsTime(),
  fIsIROC(kFALSE),
  fEdgeSuppression(kFALSE)
{
  // default constructor
}

//____________________________________________________________________
AliTPCClusterHistograms::AliTPCClusterHistograms(Int_t detector, const Char_t* comment, Int_t timeStart, Int_t timeStop, Bool_t edgeSuppression)
  : TNamed(),
  fhQmaxVsRow(0),          
  fhQtotVsRow(0),          
  fhQtotProfileVsRow(0),   
  fhQmaxProfileVsRow(0),
  fhNClustersYVsRow(0),  
  fhNClustersZVsRow(0),
  fhSigmaYVsRow(0),        
  fhSigmaZVsRow(0),          			
  fhQmaxProfileYVsRow(0), 
  fhQtotProfileYVsRow(0),
  fhSigmaYProfileYVsRow(0),
  fhSigmaZProfileYVsRow(0),
  fhQmaxProfileZVsRow(0), 
  fhQtotProfileZVsRow(0),
  fhSigmaYProfileZVsRow(0),
  fhSigmaZProfileZVsRow(0),
  fhMeanQtotVsTime(0),  
  fhQtotVsTime(0),
  fhMeanNClustersVsTime(0),    
  fhNClustersVsTime(0),        
  fhTrackQtotPerCluster(0),
  fhTrackQtotPerClusterVsPhi(0),
  fhTrackQtotPerClusterVsTheta(0),
  fhTrackMeanQtotPerClusterVsPhi(0),
  fhTrackMeanQtotPerClusterVsTheta(0),
  fhMeanNTracksVsTime(0),
  fhNEventsVsTime(0),
  fIsIROC(kFALSE),
  fEdgeSuppression(edgeSuppression)
{
  // constructor 
  
  // make name and title
  if (detector < 0 || detector >= 72) {
    AliDebug(AliLog::kError, Form("Detector %d does not exist", detector));
    return;
  }
      
  TString name(FormDetectorName(detector, edgeSuppression, comment));

  fNClustersInEvent = 0;
  fQtotInEvent      = 0;
  fMaxQtotInEvent   = 0;

  fKeepEvent = kFALSE;
  fWhyKeepEvent = TString("hi");

  fDetector = detector;
  if (detector < 36)
    fIsIROC = kTRUE; 
  
  SetName(name);
  SetTitle(Form("%s (detector %d)",name.Data(), detector));

  // rounding down to the closest hour and starting 10 hours before
  fTimeStart = 3600*UInt_t(timeStart/3600) - 36000;
  // rounding up to the closest hour
  fTimeStop  = 3600*UInt_t((3600 + timeStop)/3600);
  // each time bin covers 5 min
  Int_t nTimeBins = (fTimeStop-fTimeStart)/300;
  
  //  printf(Form(" start time: %d,  stop time: %d \n",fTimeStart, fTimeStop));

  #define BINNING_Z 250, 0, 250
  
  Float_t yRange   = 45;
  Int_t nPadRows   = 96;
  
  if (fIsIROC)
  {
    yRange   = 25;
    nPadRows = 63;
  }
  
  // 1 bin for each 0.5 cm
  Int_t nBinsY = Int_t(4*yRange);

  // do not add this hists to the directory
  Bool_t oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);

  //defining histograms and profile plots
  fhQmaxVsRow  = new TH2F("QmaxVsPadRow", "Qmax vs. pad row;Pad row;Qmax", nPadRows+2, -1.5, nPadRows+0.5, 500,  0,  500);
  fhQtotVsRow  = new TH2F("QtotVsPadRow", "Qtot vs. pad row;Pad row;Qtot", nPadRows+2, -1.5, nPadRows+0.5, 400,  0,  4000);

  fhQmaxProfileVsRow = new TProfile("MeanQmaxVsPadRow","Mean Qmax vs. pad row;Pad row;Mean Qmax",nPadRows+2, -1.5, nPadRows+0.5);
  fhQtotProfileVsRow = new TProfile("MeanQtotVsPadRow","Mean Qtot vs. pad row;Pad row;Mean Qtot",nPadRows+2, -1.5, nPadRows+0.5);
  
  fhNClustersYVsRow = new TH2F("NClusters y vs pad row","N clusters y vs pad;Pad row;y",nPadRows+2, -1.5, nPadRows+0.5, nBinsY, -yRange, yRange);
  fhNClustersZVsRow = new TH2F("NClusters z vs pad row","N clusters z vs pad;Pad row;z",nPadRows+2, -1.5, nPadRows+0.5, BINNING_Z);

  fhSigmaYVsRow = new TH2F("SigmaYVsPadRow", "Sigma Y vs. pad row;Pad row;#sigma_{Y}", nPadRows+2, -1.5, nPadRows+0.5, 100,  0,  0.5);
  fhSigmaZVsRow = new TH2F("SigmaZVsPadRow", "Sigma Z vs. pad row;Pad row;#sigma_{Z}", nPadRows+2, -1.5, nPadRows+0.5, 100,  0,  0.5);
  
  fhQmaxProfileYVsRow = new TProfile2D("MeanQmaxYVsPadRow","Mean Qmax, y vs pad row;Pad row;y",nPadRows+2, -1.5, nPadRows+0.5, nBinsY, -yRange, yRange);
  fhQtotProfileYVsRow = new TProfile2D("MeanQtotYVsPadRow","Mean Qtot, y vs pad row;Pad row;y",nPadRows+2, -1.5, nPadRows+0.5, nBinsY, -yRange, yRange);
  fhSigmaYProfileYVsRow = new TProfile2D("MeanSigmaYYVsPadRow","Mean Sigma y, y vs pad row;Pad row;y",nPadRows+2, -1.5, nPadRows+0.5, nBinsY, -yRange, yRange);
  fhSigmaZProfileYVsRow = new TProfile2D("MeanSigmaZYVsPadRow","Mean Sigma z, y vs pad row;Pad row;y",nPadRows+2, -1.5, nPadRows+0.5, nBinsY, -yRange, yRange);

  fhQmaxProfileZVsRow = new TProfile2D("MeanQmaxZVsPadRow","Mean Qmax, z vs pad row;Pad row;z",nPadRows+2, -1.5, nPadRows+0.5, BINNING_Z);
  fhQtotProfileZVsRow = new TProfile2D("MeanQtotZVsPadRow","Mean Qtot, z vs pad row;Pad row;z",nPadRows+2, -1.5, nPadRows+0.5, BINNING_Z);
  fhSigmaYProfileZVsRow = new TProfile2D("MeanSigmaYZVsPadRow","Mean Sigma y, z vs pad row;Pad row;z",nPadRows+2, -1.5, nPadRows+0.5, BINNING_Z);
  fhSigmaZProfileZVsRow = new TProfile2D("MeanSigmaZZVsPadRow","Mean Sigma z, z vs pad row;Pad row;z",nPadRows+2, -1.5, nPadRows+0.5, BINNING_Z);


  TString start(TTimeStamp(fTimeStart).AsString());
  //  TString stop(TTimeStamp(fTimeStart).AsString());
  start.Remove(26);
  
  fhMeanQtotVsTime = new TProfile("MeanQtotVsTime",Form("Mean Qtot vs. time (start %s , 1 min bins); time; Qtot",start.Data()),5*nTimeBins, fTimeStart, fTimeStop);
  fhQtotVsTime     = new TH2F("QtotVsTime",Form("Qtot vs. time (start %s , 1 min bins); time; Qtot",start.Data()),5*nTimeBins, fTimeStart, fTimeStop,400,0,2000);

  fhMeanNClustersVsTime = new TProfile("MeanNClustersVsTime",Form("Mean N Cluster vs. time (start %s , 5 min bins); time; NClusters",start.Data()),nTimeBins, fTimeStart, fTimeStop);
  fhNClustersVsTime = new TH2F("NClustersVsTime",Form("N Clusters vs. time (start %s , 5 min bins); time; NClusters",start.Data()),nTimeBins, fTimeStart, fTimeStop,400,-0.5,3999.5);

  fhQmaxProfileVsRow->SetLineWidth(2);
  fhQtotProfileVsRow->SetLineWidth(2);

  fhMeanQtotVsTime->SetLineWidth(2);

  // histograms related to tracks

  fhTrackQtotPerCluster = new TH1F("QtotPerCluster","Qtot per cluster; (Sum Qtot)/clusters",200,0,2000);
  fhTrackQtotPerCluster->SetMarkerStyle(22);
  fhTrackQtotPerCluster->SetMarkerSize(1);

  fhTrackQtotPerClusterVsPhi = new TH2F("QtotPerClusterVsPhi","QtotPerCluster vs Phi; Phi; (Sum Qtot)/clusters",40,-2,2,200,0,2000);
  fhTrackQtotPerClusterVsTheta = new TH2F("QtotPerClusterVsTheta","QtotPerCluster vs Theta; Theta; (Sum Qtot)/clusters",40,-2,2,200,0,2000);

  fhTrackMeanQtotPerClusterVsPhi = new TProfile("MeanQtotPerClusterVsPhi", "QtotPerCluster vs Phi; Phi; Mean (Sum Qtot)/clusters",40,-2,2);
  fhTrackMeanQtotPerClusterVsTheta = new TProfile("MeanQtotPerClusterVsTheta", "QtotPerCluster vs Theta; Theta; Mean (Sum Qtot)/clusters",40,-2,2);

  fhTrackMeanQtotPerClusterVsPhi->SetLineWidth(2);
  fhTrackMeanQtotPerClusterVsTheta->SetLineWidth(2);

  fhMeanNTracksVsTime = new TProfile("MeanNTracksVsTime",Form("Mean n tracks vs. time (start %s , 5 min bins); time; N tracks",start.Data()),nTimeBins, fTimeStart, fTimeStop);

  fhNEventsVsTime = new TH1F("NEventsVsTime",Form("N events vs. time (start %s , 5 min bins); time; N events",start.Data()),nTimeBins, fTimeStart, fTimeStop);

  TH1::AddDirectory(oldStatus);
}

//____________________________________________________________________
AliTPCClusterHistograms::AliTPCClusterHistograms(const AliTPCClusterHistograms& c) : TNamed(c)
{
  // copy constructor
  ((AliTPCClusterHistograms &)c).Copy(*this);
}

//____________________________________________________________________
AliTPCClusterHistograms::~AliTPCClusterHistograms()
{
  //
  // destructor
  //

  if (fhQmaxVsRow) {
    delete fhQmaxVsRow;
    fhQmaxVsRow = 0;
  }
  if (fhQtotVsRow) {
    delete fhQtotVsRow;
    fhQtotVsRow = 0; 
  }
  if (fhQmaxProfileVsRow) {
    delete fhQmaxProfileVsRow;
    fhQmaxProfileVsRow = 0;
  }
  if (fhQtotProfileVsRow) {
    delete fhQtotProfileVsRow;
    fhQtotProfileVsRow = 0;
  }
  if (fhNClustersYVsRow) {
    delete fhNClustersYVsRow;
    fhNClustersYVsRow = 0;
  }
  if (fhNClustersZVsRow) {
    delete fhNClustersZVsRow;
    fhNClustersZVsRow = 0;
  }
  if (fhSigmaYVsRow) {
    delete fhSigmaYVsRow;
    fhSigmaYVsRow = 0;
  } 
  if (fhSigmaZVsRow) {
    delete fhSigmaZVsRow;
    fhSigmaZVsRow = 0; 
  }
  if (fhQmaxProfileYVsRow) {
    delete fhQmaxProfileYVsRow;
    fhQmaxProfileYVsRow = 0;
  }
  if (fhQtotProfileYVsRow) {
    delete fhQtotProfileYVsRow;
    fhQtotProfileYVsRow = 0;
  }
  if (fhSigmaYProfileYVsRow) {
    delete fhSigmaYProfileYVsRow;
    fhSigmaYProfileYVsRow = 0;
  }
  if (fhSigmaZProfileYVsRow) {
    delete fhSigmaZProfileYVsRow;
    fhSigmaZProfileYVsRow = 0;
  }
  if (fhQmaxProfileZVsRow) {
    delete fhQmaxProfileZVsRow;
    fhQmaxProfileZVsRow = 0;
  }
  if (fhQtotProfileZVsRow) {
    delete fhQtotProfileZVsRow;
    fhQtotProfileZVsRow = 0;
  }
  if (fhSigmaYProfileZVsRow) {
    delete fhSigmaYProfileZVsRow;
    fhSigmaYProfileZVsRow = 0;
  }
  if (fhSigmaZProfileZVsRow) {
    delete fhSigmaZProfileZVsRow;
    fhSigmaZProfileZVsRow = 0;
  }
  if (fhMeanQtotVsTime) {
    delete fhMeanQtotVsTime;
    fhMeanQtotVsTime = 0;
  }
  if (fhQtotVsTime) {
    delete fhQtotVsTime;
    fhQtotVsTime = 0;
  }
  if (fhMeanNClustersVsTime) {
    delete fhMeanNClustersVsTime;
    fhMeanNClustersVsTime = 0;
  }
  if (fhNClustersVsTime) {
    delete fhNClustersVsTime;
    fhNClustersVsTime = 0;
  }
  if (fhTrackQtotPerCluster) {
    delete fhTrackQtotPerCluster;
    fhTrackQtotPerCluster = 0;
  }
  if (fhTrackQtotPerClusterVsPhi) {
    delete fhTrackQtotPerClusterVsPhi;
    fhTrackQtotPerClusterVsPhi = 0;
  }
  if (fhTrackQtotPerClusterVsTheta) {
    delete fhTrackQtotPerClusterVsTheta;
    fhTrackQtotPerClusterVsTheta = 0;
  }
  if (fhTrackMeanQtotPerClusterVsPhi) {
    delete fhTrackMeanQtotPerClusterVsPhi;
    fhTrackMeanQtotPerClusterVsPhi = 0;
  }
  if (fhTrackMeanQtotPerClusterVsTheta) {
    delete fhTrackMeanQtotPerClusterVsTheta;
    fhTrackMeanQtotPerClusterVsTheta = 0;
  }
  if (fhMeanNTracksVsTime) {
    delete fhMeanNTracksVsTime;
    fhMeanNTracksVsTime = 0;
  }  
  if (fhNEventsVsTime) {
    delete fhNEventsVsTime;
    fhNEventsVsTime = 0;
  }  
}

//____________________________________________________________________
AliTPCClusterHistograms &AliTPCClusterHistograms::operator=(const AliTPCClusterHistograms &c)
{
  // assigment operator

  if (this != &c)
    ((AliTPCClusterHistograms &) c).Copy(*this);

  return *this;
}

//____________________________________________________________________
const char* AliTPCClusterHistograms::FormDetectorName(Int_t detector, Bool_t edgeSuppression, const char* comment)
{
  //
  // creates a readable name from the detector number
  //   
  
  Int_t sector = detector%18;
  TString side;
  TString inout;
  
  if (detector<18 || ( detector>=36 && detector<54))
    side.Form("A");
  else 
    side.Form("C");
  
  if (detector<36)
    inout.Form("IROC");
  else 
    inout.Form("OROC");

  TString name;
  name.Form("sector_%s%d_%s", side.Data(), sector, inout.Data());

  if (edgeSuppression)
    name += "_noedge";
  
  if (comment)
    name += comment;

  return name; 
}

//____________________________________________________________________
Long64_t AliTPCClusterHistograms::Merge(TCollection* list)
{
  // Merge a list of AliTPCClusterHistograms objects with this (needed for
  // PROOF). 
  // Returns the number of merged objects (including this).

  if (!list)
    return 0;
  
  if (list->IsEmpty())
    return 1;

  TIterator* iter = list->MakeIterator();
  TObject* obj;

  // collections of measured and generated histograms
  TList* collectionQmaxVsRow     = new TList;
  TList* collectionQtotVsRow	 = new TList;

  TList* collectionQmaxProfileVsRow = new TList;
  TList* collectionQtotProfileVsRow = new TList;

  TList* collectionNClustersYVsRow = new TList;
  TList* collectionNClustersZVsRow = new TList;

  TList* collectionSigmaYVsRow	 = new TList;
  TList* collectionSigmaZVsRow	 = new TList;
		   			
  TList* collectionQmaxProfileYVsRow    = new TList;
  TList* collectionQtotProfileYVsRow    = new TList;
  TList* collectionSigmaYProfileYVsRow  = new TList;
  TList* collectionSigmaZProfileYVsRow  = new TList;

  TList* collectionQmaxProfileZVsRow    = new TList;
  TList* collectionQtotProfileZVsRow    = new TList;
  TList* collectionSigmaYProfileZVsRow  = new TList;
  TList* collectionSigmaZProfileZVsRow  = new TList;

  TList* collectionMeanQtotVsTime  = new TList;
  TList* collectionQtotVsTime      = new TList;

  TList* collectionMeanNClustersVsTime = new TList; 
  TList* collectionNClustersVsTime     = new TList;

  TList* collectionTrackQtotPerCluster = new TList;

  TList* collectionTrackQtotPerClusterVsPhi = new TList;
  TList* collectionTrackQtotPerClusterVsTheta = new TList;

  TList* collectionTrackMeanQtotPerClusterVsPhi = new TList;
  TList* collectionTrackMeanQtotPerClusterVsTheta = new TList;

  TList* collectionMeanNTracksVsTime = new TList;
  TList* collectionNEventsVsTime = new TList;

  Int_t count = 0;
  while ((obj = iter->Next())) {
    
    AliTPCClusterHistograms* entry = dynamic_cast<AliTPCClusterHistograms*> (obj);
    if (entry == 0) 
      continue;
    
    collectionQmaxVsRow          ->Add(entry->fhQmaxVsRow	   );
    collectionQtotVsRow	  ->Add(entry->fhQtotVsRow	   );
    
    collectionQmaxProfileVsRow   ->Add(entry->fhQmaxProfileVsRow  );
    collectionQtotProfileVsRow	  ->Add(entry->fhQtotProfileVsRow  );
    
    collectionNClustersYVsRow    ->Add(entry->fhNClustersYVsRow);
    collectionNClustersZVsRow    ->Add(entry->fhNClustersZVsRow);
    
    collectionSigmaYVsRow	  ->Add(entry->fhSigmaYVsRow	   );
    collectionSigmaZVsRow	  ->Add(entry->fhSigmaZVsRow	   );
    
    collectionQmaxProfileYVsRow  ->Add(entry->fhQmaxProfileYVsRow );
    collectionQtotProfileYVsRow  ->Add(entry->fhQtotProfileYVsRow );
    collectionSigmaYProfileYVsRow->Add(entry->fhSigmaYProfileYVsRow);
    collectionSigmaZProfileYVsRow->Add(entry->fhSigmaZProfileYVsRow);
    
    collectionQmaxProfileZVsRow  ->Add(entry->fhQmaxProfileZVsRow );
    collectionQtotProfileZVsRow  ->Add(entry->fhQtotProfileZVsRow );
    collectionSigmaYProfileZVsRow->Add(entry->fhSigmaYProfileZVsRow);
    collectionSigmaZProfileZVsRow->Add(entry->fhSigmaZProfileZVsRow);
    
    collectionMeanQtotVsTime     ->Add(entry->fhMeanQtotVsTime);
    collectionQtotVsTime         ->Add(entry->fhQtotVsTime);

    collectionMeanNClustersVsTime->Add(entry->fhMeanNClustersVsTime);  
    collectionNClustersVsTime    ->Add(entry->fhNClustersVsTime);
    
    collectionTrackQtotPerCluster->Add(entry->fhTrackQtotPerCluster);
    
    collectionTrackQtotPerClusterVsPhi->Add(entry->fhTrackQtotPerClusterVsPhi);
    collectionTrackQtotPerClusterVsTheta->Add(entry->fhTrackQtotPerClusterVsTheta);
    
    collectionTrackMeanQtotPerClusterVsPhi->Add(entry->fhTrackMeanQtotPerClusterVsPhi);
    collectionTrackMeanQtotPerClusterVsTheta->Add(entry->fhTrackMeanQtotPerClusterVsTheta);
    
    collectionMeanNTracksVsTime->Add(entry->fhMeanNTracksVsTime);
    collectionNEventsVsTime->Add(entry->fhNEventsVsTime);
    
    count++;
  }
  
  fhQmaxVsRow          ->Merge(collectionQmaxVsRow       );	   
  fhQtotVsRow          ->Merge(collectionQtotVsRow	  );	   
  
  fhQmaxProfileVsRow   ->Merge(collectionQtotProfileVsRow);
  fhQtotProfileVsRow   ->Merge(collectionQtotProfileVsRow);
  
  fhNClustersYVsRow    ->Merge(collectionNClustersYVsRow);
  fhNClustersZVsRow    ->Merge(collectionNClustersZVsRow);
  
  fhSigmaYVsRow        ->Merge(collectionSigmaYVsRow	  );	   
  fhSigmaZVsRow        ->Merge(collectionSigmaZVsRow	  );	   
  
  fhQmaxProfileYVsRow  ->Merge(collectionQmaxProfileYVsRow  ); 
  fhQtotProfileYVsRow  ->Merge(collectionQtotProfileYVsRow  );
  fhSigmaYProfileYVsRow->Merge(collectionSigmaYProfileYVsRow);
  fhSigmaZProfileYVsRow->Merge(collectionSigmaZProfileYVsRow);
  
  fhQmaxProfileZVsRow  ->Merge(collectionQmaxProfileZVsRow  ); 
  fhQtotProfileZVsRow  ->Merge(collectionQtotProfileZVsRow  );
  fhSigmaYProfileZVsRow->Merge(collectionSigmaYProfileZVsRow);
  fhSigmaZProfileZVsRow->Merge(collectionSigmaZProfileZVsRow);
  
  fhMeanQtotVsTime     ->Merge(collectionMeanQtotVsTime);
  fhQtotVsTime         ->Merge(collectionQtotVsTime);

  fhMeanNClustersVsTime->Merge(collectionMeanNClustersVsTime );
  fhNClustersVsTime    ->Merge(collectionNClustersVsTime);
  
  fhTrackQtotPerCluster->Merge(collectionTrackQtotPerCluster);
  
  fhTrackQtotPerClusterVsPhi->Merge(collectionTrackQtotPerClusterVsPhi);
  fhTrackQtotPerClusterVsTheta->Merge(collectionTrackQtotPerClusterVsTheta);
  
  fhTrackMeanQtotPerClusterVsPhi->Merge(collectionTrackMeanQtotPerClusterVsPhi);
  fhTrackMeanQtotPerClusterVsTheta->Merge(collectionTrackMeanQtotPerClusterVsTheta);
  
  fhMeanNTracksVsTime->Merge(collectionMeanNTracksVsTime);
  fhNEventsVsTime->Merge(collectionNEventsVsTime);
  
  delete collectionQmaxVsRow;          
  delete collectionQtotVsRow;  
  
  delete collectionQmaxProfileVsRow;
  delete collectionQtotProfileVsRow;
  
  delete collectionNClustersYVsRow;
  delete collectionNClustersZVsRow;
  
  delete collectionSigmaYVsRow;	  
  delete collectionSigmaZVsRow;	  
  
  delete collectionQmaxProfileYVsRow;  
  delete collectionQtotProfileYVsRow;  
  delete collectionSigmaYProfileYVsRow;
  delete collectionSigmaZProfileYVsRow;
  
  delete collectionQmaxProfileZVsRow;  
  delete collectionQtotProfileZVsRow;  
  delete collectionSigmaYProfileZVsRow;
  delete collectionSigmaZProfileZVsRow;
  
  delete collectionMeanQtotVsTime;
  delete collectionQtotVsTime;

  delete collectionMeanNClustersVsTime; 
  delete collectionNClustersVsTime;     
  
  delete collectionTrackQtotPerCluster;
  
  delete collectionTrackQtotPerClusterVsPhi; 
  delete collectionTrackQtotPerClusterVsTheta;
  
  delete collectionTrackMeanQtotPerClusterVsPhi; 
  delete collectionTrackMeanQtotPerClusterVsTheta;
  
  delete collectionMeanNTracksVsTime;
  delete collectionNEventsVsTime;

  return count+1;
}


//____________________________________________________________________
void AliTPCClusterHistograms::FillCluster(AliTPCclusterMI* cluster, Int_t time) {
  //
  // Fills the different histograms with the information from the cluster.
  //

  Int_t padRow =   cluster->GetRow(); 
  Float_t qMax =   cluster->GetMax();
  Float_t qTot =   cluster->GetQ();
  Float_t sigmaY = cluster->GetSigmaY2();
  Float_t sigmaZ = cluster->GetSigmaZ2();
  Float_t y      = cluster->GetY();
  Float_t z      = cluster->GetZ();

  // check if this is ok!!!
  z = TMath::Abs(z);

  if (qMax<=0) {
    printf(Form("\n WARNING: Hi Marian! How can we have Qmax = %f ??? \n \n", qMax));
  }
  if (qTot<=0) {
    printf(Form("\n WARNING: Hi Marian! How can we have Qtot = %f ??? \n \n ", qTot));
  }   
  
  // check if the cluster is accepted
  if (fEdgeSuppression)
    if (IsClusterOnEdge(cluster))
      return;

  fNClustersInEvent++;
  fQtotInEvent = fQtotInEvent + qTot;
  if (qTot > fMaxQtotInEvent)
    fMaxQtotInEvent = qTot;

  fhQmaxVsRow           ->Fill(padRow, qMax);
  fhQtotVsRow           ->Fill(padRow, qTot);

  fhQmaxProfileVsRow    ->Fill(padRow, qMax);
  fhQtotProfileVsRow    ->Fill(padRow, qTot);

  fhNClustersYVsRow     ->Fill(padRow, y, 1);
  fhNClustersZVsRow     ->Fill(padRow, z, 1);
  			
  fhSigmaYVsRow         ->Fill(padRow, sigmaY);
  fhSigmaZVsRow         ->Fill(padRow, sigmaZ);
  			
  fhQmaxProfileYVsRow   ->Fill(padRow, y, qMax);
  fhQtotProfileYVsRow   ->Fill(padRow, y, qTot); 
  fhSigmaYProfileYVsRow ->Fill(padRow, y, sigmaY);
  fhSigmaZProfileYVsRow ->Fill(padRow, y, sigmaZ);

  fhQmaxProfileZVsRow   ->Fill(padRow, z, qMax);
  fhQtotProfileZVsRow   ->Fill(padRow, z, qTot); 
  fhSigmaYProfileZVsRow ->Fill(padRow, z, sigmaY);
  fhSigmaZProfileZVsRow ->Fill(padRow, z, sigmaZ);

  if (time>0 & fTimeStart>0 & fTimeStop>0 & time>fTimeStart) {
    //Float_t timeFraction = (time - fTimeStart)/(fTimeStop-fTimeStart); 

    fhMeanQtotVsTime->Fill(time,qTot);
    fhQtotVsTime->Fill(time,qTot);
  }
}

//____________________________________________________________________
void AliTPCClusterHistograms::FillTrack(const AliTPCseed* seed) {
  //
  // fill histograms related to tracks
  //

  Float_t totalQtot = 0;
  Int_t   nClusters = 0;
  for (Int_t clusterID = 0; clusterID < 160; clusterID++) {
    AliTPCclusterMI* cluster = seed->GetClusterPointer(clusterID);
    if (!cluster)
      continue;
    
    // only use clusters within this detector
    if (cluster->GetDetector()!=fDetector)
      continue;
    
    // check if the cluster is accepted
    if (fEdgeSuppression)
      if (IsClusterOnEdge(cluster))
	return;

    Int_t padRow =   cluster->GetRow(); 
    Float_t qMax =   cluster->GetMax();
    Float_t qTot =   cluster->GetQ();    

    nClusters++;
    totalQtot += qTot;
    
  }
  if (nClusters==0) 
    return;
  
  Float_t meanQtot = totalQtot/nClusters;
  
  // azimuthal angle 
  Float_t phi    =  TMath::ASin(seed->GetSnp() + seed->GetAlpha());
  // angle with respect to the central membrane
  Float_t theta  =  TMath::ATan(seed->GetTgl());

  fhTrackQtotPerCluster->Fill(meanQtot);

  fhTrackMeanQtotPerClusterVsPhi->Fill(phi,   meanQtot);
  fhTrackMeanQtotPerClusterVsTheta->Fill(theta, meanQtot);

  fhTrackQtotPerClusterVsPhi->Fill(phi, meanQtot);
  fhTrackQtotPerClusterVsTheta->Fill(theta, meanQtot);
}

//____________________________________________________________________
void AliTPCClusterHistograms::FillEvent(Int_t time, Int_t nTracks) {
  //
  // fill event
  //

  fhMeanNTracksVsTime->Fill(time, nTracks);

  //  fhNEventsVsTime->Fill(time);
}


//____________________________________________________________________
Bool_t AliTPCClusterHistograms::IsClusterOnEdge(AliTPCclusterMI* clusterMI) {
  //
  // check if the cluster is on the edge
  //

  Int_t padRow =   clusterMI->GetRow(); 
  Float_t y      = clusterMI->GetY();
  
  Float_t limit = 0;
  if (fIsIROC)
    {
      limit = 12 + padRow * (20.0 - 12.0) / 63; 
    }
  else
    limit = 16 + padRow * (36.0 - 16.0) / 96;
  
  if (TMath::Abs(y) > limit)
    return kTRUE;
  
  return kFALSE;
}

//____________________________________________________________________
Float_t AliTPCClusterHistograms::DistanceToEdge(AliTPCclusterMI* clusterMI) {
  //
  // get the y-distance to closest edge 
  //
  
  Int_t detector = clusterMI->GetDetector();
  Int_t padRow   = clusterMI->GetRow(); 
  Float_t y      = clusterMI->GetY();
  
  Float_t yEdge = -9999;
  Float_t d     = 0;
  
  // IROC
  if (detector < 36) {    
    yEdge = 14 + padRow * 0.1333;
    
  }    
  else { // OROC
    if (padRow<64) // small pads
      yEdge = 22.5 + padRow * 0.1746;
    else           // large pads
      yEdge = 34.0 + (padRow-64) * 0.2581;      
  }
  if (y<=0) yEdge = -yEdge;
  
  d = yEdge - y;
  
  return d;
}


//____________________________________________________________________
Bool_t AliTPCClusterHistograms::KeepThisEvent(TString& why) {
  //
  // is this event interesting? 
  //
  // the criteria are ...
  // 
  
  if (fKeepEvent) {
    why = TString(fWhyKeepEvent);
    return kTRUE;
  }

  if (fNClustersInEvent>20000) {
    why.Append("_moreThan20000clusters");
    fWhyKeepEvent = TString(why);
    fKeepEvent = kTRUE;
    return kTRUE;
  }
  
  if (fMaxQtotInEvent>10000) {
    why.Append("_clusterWithQtot20000plus");
    fWhyKeepEvent = TString(why);
    fKeepEvent = kTRUE;
    return kTRUE;
  }

  if (gRandom->Uniform()<0.001) {
    why.Append("_random");
    fWhyKeepEvent = TString(why);
    fKeepEvent = kTRUE;
    return kTRUE;
  }

  return kFALSE;
}

//____________________________________________________________________
void AliTPCClusterHistograms::StartEvent() {
  //
  // reset counters 
  // 
 
  fNClustersInEvent  = 0; 
  fQtotInEvent       = 0; 
  fMaxQtotInEvent    = 0; 
  fKeepEvent         = kFALSE; 
  fWhyKeepEvent      = TString("");
 
}


//____________________________________________________________________
void AliTPCClusterHistograms::FinishEvent(Int_t timeStamp) {
  //
  // fill histograms related to the event
  //
 
  fhMeanNClustersVsTime->Fill(timeStamp, fNClustersInEvent); 
  fhNClustersVsTime    ->Fill(timeStamp, fNClustersInEvent); 

  fhNEventsVsTime->Fill(timeStamp);
 
}


//____________________________________________________________________
void AliTPCClusterHistograms::SaveHistograms()
{
  //
  // saves the histograms
  //

  gDirectory->mkdir(fName.Data());
  gDirectory->cd(fName.Data());

  fhQmaxVsRow           ->Write();
  fhQtotVsRow           ->Write();

  fhQmaxProfileVsRow    ->Write();
  fhQtotProfileVsRow    ->Write();

  fhNClustersYVsRow     ->Write();
  fhNClustersZVsRow     ->Write();
  			
  fhSigmaYVsRow         ->Write();
  fhSigmaZVsRow         ->Write();
  			
  fhQmaxProfileYVsRow   ->Write();
  fhQtotProfileYVsRow   ->Write();
  fhSigmaYProfileYVsRow ->Write();
  fhSigmaZProfileYVsRow ->Write();

  fhQmaxProfileZVsRow   ->Write();
  fhQtotProfileZVsRow   ->Write();
  fhSigmaYProfileZVsRow ->Write();
  fhSigmaZProfileZVsRow ->Write();

  TNamed* comment = new TNamed("comment", fCommentToHistograms.Data());
  comment->Write();

  if (fhMeanQtotVsTime->GetEntries()>0)
    fhMeanQtotVsTime->Write();

  if (fhQtotVsTime->GetEntries()>0)
    fhQtotVsTime->Write();

  if (fhMeanNClustersVsTime->GetEntries()>0)
    fhMeanNClustersVsTime->Write();

  if (fhNClustersVsTime->GetEntries()>0)
    fhNClustersVsTime->Write();

  if (fhNEventsVsTime->GetEntries()>0)
    fhNEventsVsTime->Write();

  gDirectory->mkdir("track_hists");
  gDirectory->cd("track_hists");

  fhTrackQtotPerCluster->Write();

  fhTrackQtotPerClusterVsPhi->Write();
  fhTrackQtotPerClusterVsTheta->Write();

  fhTrackMeanQtotPerClusterVsPhi->Write();
  fhTrackMeanQtotPerClusterVsTheta->Write();

  fhMeanNTracksVsTime->Write();

  gDirectory->cd("../");

  gDirectory->cd("../");

}


//____________________________________________________________________
TCanvas* AliTPCClusterHistograms::DrawHistograms(const Char_t* /*opt*/) {
  //
  // Draws some histograms and save the canvas as eps and gif file.
  //  

  TCanvas* c = new TCanvas(fName.Data(), fName.Data(), 1200, 1000);

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  gStyle->SetPadLeftMargin(0.1);

  c->Divide(3,3);

  c->Draw();  

  c->cd(1);
  
  // this is not really a nice way to do it...
  c->GetPad(1)->Delete();
  
  TLatex* tName = new TLatex(0.05,0.9,fName.Data());
  tName->SetTextSize(0.02);
  tName->DrawClone();
  
  TLatex* tEdge;
  if (fEdgeSuppression) 
    tEdge = new TLatex(0.05,0.85,"(edges cut)");
  else 
    tEdge = new TLatex(0.05,0.85,"(no edge cut)");
  
  tEdge->SetTextSize(0.015);
  tEdge->DrawClone();
  


  c->cd(2);
  fhNClustersYVsRow->Draw("colz");

  c->cd(3);
  fhNClustersZVsRow->Draw("colz");

  c->cd(4);
  fhQmaxVsRow->Draw("colz");
  fhQmaxProfileVsRow->Draw("same");

  c->cd(5);
  fhQtotVsRow->Draw("colz"); 
  fhQtotProfileVsRow->Draw("same");       
  			
  c->cd(6);
  fhQtotProfileYVsRow   ->Draw("colz");

  c->cd(7);
  fhQtotProfileZVsRow   ->Draw("colz");

  c->cd(8);
  fhQmaxProfileYVsRow   ->Draw("colz");

  c->cd(9);
  fhQmaxProfileZVsRow   ->Draw("colz");    
  			
  return c;
}
