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

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// TRD cluster finder                                                        //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <TF1.h>
#include <TTree.h>
#include <TH1.h>
#include <TFile.h>
#include <TClonesArray.h>
#include <TObjArray.h>

#include "AliRunLoader.h"
#include "AliLoader.h"
#include "AliRawReader.h"
#include "AliLog.h"
#include "AliAlignObj.h"

#include "AliTRDclusterizer.h"
#include "AliTRDcluster.h"
#include "AliTRDReconstructor.h"
#include "AliTRDgeometry.h"
#include "AliTRDarraySignal.h"
#include "AliTRDarrayDictionary.h"
#include "AliTRDarrayADC.h"
#include "AliTRDdigitsManager.h"
#include "AliTRDrawData.h"
#include "AliTRDcalibDB.h"
#include "AliTRDrecoParam.h"
#include "AliTRDCommonParam.h"
#include "AliTRDtransform.h"
#include "AliTRDSignalIndex.h"
#include "AliTRDrawStreamBase.h"
#include "AliTRDfeeParam.h"

#include "Cal/AliTRDCalROC.h"
#include "Cal/AliTRDCalDet.h"
#include "Cal/AliTRDCalSingleChamberStatus.h"

ClassImp(AliTRDclusterizer)

//_____________________________________________________________________________
AliTRDclusterizer::AliTRDclusterizer(AliTRDReconstructor *rec)
  :TNamed()
  ,fReconstructor(rec)  
  ,fRunLoader(NULL)
  ,fClusterTree(NULL)
  ,fRecPoints(NULL)
  ,fTrackletTree(NULL)
  ,fDigitsManager(NULL)
  ,fTrackletContainer(NULL)
  ,fAddLabels(kTRUE)
  ,fRawVersion(2)
  ,fIndexesOut(NULL)
  ,fIndexesMaxima(NULL)
  ,fTransform(new AliTRDtransform(0))
  ,fLUTbin(0)
  ,fLUT(NULL)
{
  //
  // AliTRDclusterizer default constructor
  //

  AliTRDcalibDB *trd = 0x0;
  if (!(trd = AliTRDcalibDB::Instance())) {
    AliFatal("Could not get calibration object");
  }

  fRawVersion = AliTRDfeeParam::Instance()->GetRAWversion();

  // Initialize debug stream
  if(fReconstructor){
    if(fReconstructor->GetStreamLevel(AliTRDReconstructor::kClusterizer) > 1){
      TDirectory *savedir = gDirectory; 
      //fgDebugStreamer    = new TTreeSRedirector("TRD.ClusterizerDebug.root");
      savedir->cd();
    }
  }

}

//_____________________________________________________________________________
AliTRDclusterizer::AliTRDclusterizer(const Text_t *name, const Text_t *title, AliTRDReconstructor *rec)
  :TNamed(name,title)
  ,fReconstructor(rec)
  ,fRunLoader(NULL)
  ,fClusterTree(NULL)
  ,fRecPoints(NULL)
  ,fTrackletTree(NULL)
  ,fDigitsManager(new AliTRDdigitsManager())
  ,fTrackletContainer(NULL)
  ,fAddLabels(kTRUE)
  ,fRawVersion(2)
  ,fIndexesOut(NULL)
  ,fIndexesMaxima(NULL)
  ,fTransform(new AliTRDtransform(0))
  ,fLUTbin(0)
  ,fLUT(NULL)
{
  //
  // AliTRDclusterizer constructor
  //

  AliTRDcalibDB *trd = 0x0;
  if (!(trd = AliTRDcalibDB::Instance())) {
    AliFatal("Could not get calibration object");
  }

  // Initialize debug stream
  if (fReconstructor){
    if (fReconstructor->GetStreamLevel(AliTRDReconstructor::kClusterizer) > 1){
      TDirectory *savedir = gDirectory; 
      //fgDebugStreamer    = new TTreeSRedirector("TRD.ClusterizerDebug.root");
      savedir->cd();
    }
  }

  fDigitsManager->CreateArrays();

  fRawVersion = AliTRDfeeParam::Instance()->GetRAWversion();

  FillLUT();

}

//_____________________________________________________________________________
AliTRDclusterizer::AliTRDclusterizer(const AliTRDclusterizer &c)
  :TNamed(c)
  ,fReconstructor(c.fReconstructor)
  ,fRunLoader(NULL)
  ,fClusterTree(NULL)
  ,fRecPoints(NULL)
  ,fTrackletTree(NULL)
  ,fDigitsManager(NULL)
  ,fTrackletContainer(NULL)
  ,fAddLabels(kTRUE)
  ,fRawVersion(2)
  ,fIndexesOut(NULL)
  ,fIndexesMaxima(NULL)
  ,fTransform(NULL)
  ,fLUTbin(0)
  ,fLUT(0)
{
  //
  // AliTRDclusterizer copy constructor
  //

  FillLUT();

}

//_____________________________________________________________________________
AliTRDclusterizer::~AliTRDclusterizer()
{
  //
  // AliTRDclusterizer destructor
  //

  if (fRecPoints/* && IsClustersOwner()*/){
    fRecPoints->Delete();
    delete fRecPoints;
  }

  if (fDigitsManager) {
    delete fDigitsManager;
    fDigitsManager = NULL;
  }

  if (fTrackletContainer){
    delete fTrackletContainer;
    fTrackletContainer = NULL;
  }

  if (fIndexesOut){
    delete fIndexesOut;
    fIndexesOut    = NULL;
  }

  if (fIndexesMaxima){
    delete fIndexesMaxima;
    fIndexesMaxima = NULL;
  }

  if (fTransform){
    delete fTransform;
    fTransform     = NULL;
  }

  if (fLUT) {
    delete [] fLUT;
    fLUT           = NULL;
  }

}

//_____________________________________________________________________________
AliTRDclusterizer &AliTRDclusterizer::operator=(const AliTRDclusterizer &c)
{
  //
  // Assignment operator
  //

  if (this != &c) 
    {
      ((AliTRDclusterizer &) c).Copy(*this);
    }

  return *this;

}

//_____________________________________________________________________________
void AliTRDclusterizer::Copy(TObject &c) const
{
  //
  // Copy function
  //

  ((AliTRDclusterizer &) c).fClusterTree   = NULL;
  ((AliTRDclusterizer &) c).fRecPoints     = NULL;  
  ((AliTRDclusterizer &) c).fTrackletTree  = NULL;
  ((AliTRDclusterizer &) c).fDigitsManager = NULL;
  ((AliTRDclusterizer &) c).fTrackletContainer = NULL;
  ((AliTRDclusterizer &) c).fAddLabels     = fAddLabels;
  ((AliTRDclusterizer &) c).fRawVersion    = fRawVersion;
  ((AliTRDclusterizer &) c).fIndexesOut    = NULL;
  ((AliTRDclusterizer &) c).fIndexesMaxima = NULL;
  ((AliTRDclusterizer &) c).fTransform     = NULL;
  ((AliTRDclusterizer &) c).fLUTbin        = 0;
  ((AliTRDclusterizer &) c).fLUT           = NULL;

}

//_____________________________________________________________________________
Bool_t AliTRDclusterizer::Open(const Char_t *name, Int_t nEvent)
{
  //
  // Opens the AliROOT file. Output and input are in the same file
  //

  TString evfoldname = AliConfig::GetDefaultEventFolderName();
  fRunLoader         = AliRunLoader::GetRunLoader(evfoldname);

  if (!fRunLoader) {
    fRunLoader = AliRunLoader::Open(name);
  }

  if (!fRunLoader) {
    AliError(Form("Can not open session for file %s.",name));
    return kFALSE;
  }

  OpenInput(nEvent);
  OpenOutput();

  return kTRUE;

}

//_____________________________________________________________________________
Bool_t AliTRDclusterizer::OpenOutput()
{
  //
  // Open the output file
  //

  if (!fReconstructor->IsWritingClusters()) return kTRUE;

  TObjArray *ioArray = 0x0; 

  AliLoader* loader = fRunLoader->GetLoader("TRDLoader");
  loader->MakeTree("R");

  fClusterTree = loader->TreeR();
  fClusterTree->Branch("TRDcluster", "TObjArray", &ioArray, 32000, 0);

  return kTRUE;

}

//_____________________________________________________________________________
Bool_t AliTRDclusterizer::OpenOutput(TTree *clusterTree)
{
  //
  // Connect the output tree
  //

  // clusters writing
  if (fReconstructor->IsWritingClusters()){
    TObjArray *ioArray = 0x0;
    fClusterTree = clusterTree;
    fClusterTree->Branch("TRDcluster", "TObjArray", &ioArray, 32000, 0);
  }

  // tracklet writing
  if (fReconstructor->IsWritingTracklets()){
    TString evfoldname = AliConfig::GetDefaultEventFolderName();
    fRunLoader         = AliRunLoader::GetRunLoader(evfoldname);

    if (!fRunLoader) {
      fRunLoader = AliRunLoader::Open("galice.root");
    }
    if (!fRunLoader) {
      AliError(Form("Can not open session for file galice.root."));
      return kFALSE;
    }

    UInt_t **leaves = new UInt_t *[2];
    AliDataLoader *dl = fRunLoader->GetLoader("TRDLoader")->GetDataLoader("tracklets");
    if (!dl) {
      AliError("Could not get the tracklets data loader!");
      dl = new AliDataLoader("TRD.Tracklets.root","tracklets", "tracklets");
      fRunLoader->GetLoader("TRDLoader")->AddDataLoader(dl);
    }
    else {
      fTrackletTree = dl->Tree();
      if (!fTrackletTree)
        {
        dl->MakeTree();
        fTrackletTree = dl->Tree();
        }
      TBranch *trkbranch = fTrackletTree->GetBranch("trkbranch");
      if (!trkbranch)
        fTrackletTree->Branch("trkbranch",leaves[0],"det/i:side/i:tracklets[256]/i");
    }
  }

  return kTRUE;

}

//_____________________________________________________________________________
Bool_t AliTRDclusterizer::OpenInput(Int_t nEvent)
{
  //
  // Opens a ROOT-file with TRD-hits and reads in the digits-tree
  //

  // Import the Trees for the event nEvent in the file
  fRunLoader->GetEvent(nEvent);
  
  return kTRUE;

}

//_____________________________________________________________________________
Bool_t AliTRDclusterizer::WriteClusters(Int_t det)
{
  //
  // Fills TRDcluster branch in the tree with the clusters 
  // found in detector = det. For det=-1 writes the tree. 
  //

  if ((det <                      -1) || 
      (det >= AliTRDgeometry::Ndet())) {
    AliError(Form("Unexpected detector index %d.\n",det));
    return kFALSE;
  }

  TObjArray *ioArray = new TObjArray(400);
  TBranch *branch = fClusterTree->GetBranch("TRDcluster");
  if (!branch) {
    branch = fClusterTree->Branch("TRDcluster","TObjArray",&ioArray,32000,0);
  } else branch->SetAddress(&ioArray);
  
  Int_t nRecPoints = RecPoints()->GetEntriesFast();
  if(det >= 0){
    for (Int_t i = 0; i < nRecPoints; i++) {
      AliTRDcluster *c = (AliTRDcluster *) RecPoints()->UncheckedAt(i);
      if(det != c->GetDetector()) continue;
      ioArray->AddLast(c);
    }
    fClusterTree->Fill();
  } else {
    
    Int_t detOld = -1;
    for (Int_t i = 0; i < nRecPoints; i++) {
      AliTRDcluster *c = (AliTRDcluster *) RecPoints()->UncheckedAt(i);
      if(c->GetDetector() != detOld){
        fClusterTree->Fill();
        ioArray->Clear();
        detOld = c->GetDetector();
      } 
      ioArray->AddLast(c);
    }
  }
  delete ioArray;

  return kTRUE;  

}

//_____________________________________________________________________________
Bool_t AliTRDclusterizer::WriteTracklets(Int_t det)
{
  //
  // Write the raw data tracklets into seperate file
  //

  UInt_t **leaves = new UInt_t *[2];
  for (Int_t i=0; i<2 ;i++){
    leaves[i] = new UInt_t[258];
    leaves[i][0] = det; // det
    leaves[i][1] = i;   // side
    memcpy(leaves[i]+2, fTrackletContainer[i], sizeof(UInt_t) * 256);
  }

  if (!fTrackletTree){
    AliDataLoader *dl = fRunLoader->GetLoader("TRDLoader")->GetDataLoader("tracklets");
    dl->MakeTree();
    fTrackletTree = dl->Tree();
  }

  TBranch *trkbranch = fTrackletTree->GetBranch("trkbranch");
  if (!trkbranch) {
    trkbranch = fTrackletTree->Branch("trkbranch",leaves[0],"det/i:side/i:tracklets[256]/i");
  }

  for (Int_t i=0; i<2; i++){
    if (leaves[i][2]>0) {
      trkbranch->SetAddress(leaves[i]);
      fTrackletTree->Fill();
    }
  }

  AliDataLoader *dl = fRunLoader->GetLoader("TRDLoader")->GetDataLoader("tracklets");
  dl->WriteData("OVERWRITE");
  //dl->Unload();
  delete [] leaves;

  return kTRUE;

}

//_____________________________________________________________________________
void AliTRDclusterizer::ResetHelperIndexes(AliTRDSignalIndex *indexesIn)
{
  // 
  // Reset the helper indexes
  //

  if (fIndexesOut)
    {
      // carefull here - we assume that only row number may change - most probable
      if (indexesIn->GetNrow() <= fIndexesOut->GetNrow())
  fIndexesOut->ResetContent();
      else
  fIndexesOut->ResetContentConditional(indexesIn->GetNrow()
                                          , indexesIn->GetNcol()
                                          , indexesIn->GetNtime());
    }
  else
    {
      fIndexesOut = new AliTRDSignalIndex(indexesIn->GetNrow()
                                        , indexesIn->GetNcol() 
                                        , indexesIn->GetNtime());
    }
  
  if (fIndexesMaxima)
    {
      // carefull here - we assume that only row number may change - most probable
      if (indexesIn->GetNrow() <= fIndexesMaxima->GetNrow()) 
        {
    fIndexesMaxima->ResetContent();
  }
      else
  {
    fIndexesMaxima->ResetContentConditional(indexesIn->GetNrow()
                                                , indexesIn->GetNcol()
                                                , indexesIn->GetNtime());
  }
    }
  else
    {
      fIndexesMaxima = new AliTRDSignalIndex(indexesIn->GetNrow()
                                          , indexesIn->GetNcol()
                                          , indexesIn->GetNtime());
    }

}

//_____________________________________________________________________________
Bool_t AliTRDclusterizer::ReadDigits()
{
  //
  // Reads the digits arrays from the input aliroot file
  //

  if (!fRunLoader) {
    AliError("No run loader available");
    return kFALSE;
  }

  AliLoader* loader = fRunLoader->GetLoader("TRDLoader");
  if (!loader->TreeD()) {
    loader->LoadDigits();
  }

  // Read in the digit arrays
  return (fDigitsManager->ReadDigits(loader->TreeD()));

}

//_____________________________________________________________________________
Bool_t AliTRDclusterizer::ReadDigits(TTree *digitsTree)
{
  //
  // Reads the digits arrays from the input tree
  //

  // Read in the digit arrays
  return (fDigitsManager->ReadDigits(digitsTree));

}

//_____________________________________________________________________________
Bool_t AliTRDclusterizer::ReadDigits(AliRawReader *rawReader)
{
  //
  // Reads the digits arrays from the ddl file
  //

  AliTRDrawData raw;
  fDigitsManager = raw.Raw2Digits(rawReader);

  return kTRUE;

}

//_____________________________________________________________________________
Bool_t AliTRDclusterizer::MakeClusters()
{
  //
  // Creates clusters from digits
  //

  // Propagate info from the digits manager
  if (fAddLabels == kTRUE){
    fAddLabels = fDigitsManager->UsesDictionaries();
  }
  
  Bool_t fReturn = kTRUE;
  for (Int_t i = 0; i < AliTRDgeometry::kNdet; i++){
  
    AliTRDarrayADC *digitsIn = (AliTRDarrayADC*) fDigitsManager->GetDigits(i); //mod     
    // This is to take care of switched off super modules
    if (!digitsIn->HasData()) continue;
    digitsIn->Expand();
    digitsIn->DeleteNegatives();  // Restore digits array to >=0 values
    AliTRDSignalIndex* indexes = fDigitsManager->GetIndexes(i);
    if (indexes->IsAllocated() == kFALSE){
      fDigitsManager->BuildIndexes(i);
    }
  
    Bool_t fR = kFALSE;
    if (indexes->HasEntry()){
      if (fAddLabels){
        for (Int_t iDict = 0; iDict < AliTRDdigitsManager::kNDict; iDict++){
          AliTRDarrayDictionary *tracksIn = 0; //mod
          tracksIn = (AliTRDarrayDictionary *) fDigitsManager->GetDictionary(i,iDict);  //mod
          tracksIn->Expand();
        }
      }
      fR = MakeClusters(i);
      fReturn = fR && fReturn;
    }
  
    //if (fR == kFALSE){
    //  if(IsWritingClusters()) WriteClusters(i);
    //  ResetRecPoints();
    //}
        
    // No compress just remove
    fDigitsManager->RemoveDigits(i);
    fDigitsManager->RemoveDictionaries(i);      
    fDigitsManager->ClearIndexes(i);  
  }
  
  if(fReconstructor->IsWritingClusters()) WriteClusters(-1);

  AliInfo(Form("Number of found clusters : %d", RecPoints()->GetEntriesFast())); 

  return fReturn;

}

//_____________________________________________________________________________
Bool_t AliTRDclusterizer::Raw2Clusters(AliRawReader *rawReader)
{
  //
  // Creates clusters from raw data
  //

  return Raw2ClustersChamber(rawReader);

}

//_____________________________________________________________________________
Bool_t AliTRDclusterizer::Raw2ClustersChamber(AliRawReader *rawReader)
{
  //
  // Creates clusters from raw data
  //

  // Create the digits manager
  if (!fDigitsManager){
    fDigitsManager = new AliTRDdigitsManager();
    fDigitsManager->CreateArrays();
  }

  fDigitsManager->SetUseDictionaries(fAddLabels);

  // tracklet container for raw tracklet writing
  if (!fTrackletContainer && fReconstructor->IsWritingTracklets()) {
    // maximum tracklets for one HC
    const Int_t kTrackletChmb=256;
    fTrackletContainer = new UInt_t *[2];
    fTrackletContainer[0] = new UInt_t[kTrackletChmb]; 
    fTrackletContainer[1] = new UInt_t[kTrackletChmb]; 
  }


  AliTRDrawStreamBase *pinput = AliTRDrawStreamBase::GetRawStream(rawReader);
  AliTRDrawStreamBase &input = *pinput;

  AliInfo(Form("Stream version: %s", input.IsA()->GetName()));
  
  Int_t det    = 0;
  while ((det = input.NextChamber(fDigitsManager,fTrackletContainer)) >= 0){
    Bool_t iclusterBranch = kFALSE;
    if (fDigitsManager->GetIndexes(det)->HasEntry()){
    iclusterBranch = MakeClusters(det);
    }

    fDigitsManager->RemoveDigits(det);
    fDigitsManager->RemoveDictionaries(det);      
    fDigitsManager->ClearIndexes(det);

    if (!fReconstructor->IsWritingTracklets()) continue;
    if (*(fTrackletContainer[0]) > 0 || *(fTrackletContainer[1]) > 0) WriteTracklets(det);
  }
  if (fReconstructor->IsWritingTracklets()){
    delete [] fTrackletContainer[0];
    delete [] fTrackletContainer[1];
    delete [] fTrackletContainer;
    fTrackletContainer = NULL;
  }

  if(fReconstructor->IsWritingClusters()) WriteClusters(-1);

  delete fDigitsManager;
  fDigitsManager = NULL;

  delete pinput;
  pinput = NULL;

  AliInfo(Form("Number of found clusters : %d", RecPoints()->GetEntriesFast())); 
  return kTRUE;

}

//_____________________________________________________________________________
UChar_t AliTRDclusterizer::GetStatus(Short_t &signal)
{
  //
  // Check if a pad is masked
  //

  UChar_t status = 0;

  if(signal>0 && TESTBIT(signal, 10)){
    CLRBIT(signal, 10);
    for(int ibit=0; ibit<4; ibit++){
      if(TESTBIT(signal, 11+ibit)){
        SETBIT(status, ibit);
        CLRBIT(signal, 11+ibit);
      } 
    }
  }
  return status;
}

//_____________________________________________________________________________
void AliTRDclusterizer::SetPadStatus(UChar_t status, UChar_t &out){
  //
  // Set the pad status into out
  // First three bits are needed for the position encoding
  //
  status = status << 3;
  out |= status;
}

//_____________________________________________________________________________
UChar_t AliTRDclusterizer::GetPadStatus(UChar_t encoding){
  //
  // return the staus encoding of the corrupted pad
  //
  return static_cast<UChar_t>(encoding >> 3);
}

//_____________________________________________________________________________
Int_t AliTRDclusterizer::GetCorruption(UChar_t encoding){
  //
  // Return the position of the corruption
  //
  return encoding & 7;
}

//_____________________________________________________________________________
Bool_t AliTRDclusterizer::MakeClusters(Int_t det)
{
  //
  // Generates the cluster.
  //

  // Get the digits
  //   digits should be expanded beforehand!
  //   digitsIn->Expand();
  AliTRDarrayADC *digitsIn = (AliTRDarrayADC *) fDigitsManager->GetDigits(det); //mod     
  
  // This is to take care of switched off super modules
  if (!digitsIn->HasData()) 
    {
      return kFALSE;
    }

  AliTRDSignalIndex *indexesIn = fDigitsManager->GetIndexes(det);
  if (indexesIn->IsAllocated() == kFALSE)
    {
      AliError("Indexes do not exist!");
      return kFALSE;      
    }
    
  AliTRDcalibDB  *calibration = AliTRDcalibDB::Instance();
  if (!calibration) 
    {
      AliFatal("No AliTRDcalibDB instance available\n");
      return kFALSE;  
    }

  // ADC thresholds
  // There is no ADC threshold anymore, and simParam should not be used in clusterizer. KO
  Float_t adcThreshold   = 0; 

  if (!fReconstructor){
    AliError("Reconstructor not set\n");
    return kFALSE;
  }

  // Threshold value for the maximum
  Float_t maxThresh            = fReconstructor->GetRecoParam()->GetClusMaxThresh();
  // Threshold value for the digit signal
  Float_t sigThresh            = fReconstructor->GetRecoParam()->GetClusSigThresh();

  // Threshold value for the maximum ( cut noise)
  Float_t minMaxCutSigma       = fReconstructor->GetRecoParam()->GetMinMaxCutSigma();
  // Threshold value for the sum pad ( cut noise)
  Float_t minLeftRightCutSigma = fReconstructor->GetRecoParam()->GetMinLeftRightCutSigma();

  // Iteration limit for unfolding procedure
  const Float_t kEpsilon = 0.01;             
  const Int_t   kNclus   = 3;  
  const Int_t   kNsig    = 5;

  Int_t    iUnfold       = 0;  
  Double_t ratioLeft     = 1.0;
  Double_t ratioRight    = 1.0;

  Double_t padSignal[kNsig];   
  Double_t clusterSignal[kNclus];

  Int_t istack  = indexesIn->GetStack();
  Int_t ilayer  = indexesIn->GetLayer();
  Int_t isector = indexesIn->GetSM();

  // Start clustering in the chamber

  Int_t idet  = AliTRDgeometry::GetDetector(ilayer,istack,isector);
  if (idet != det) {
    AliError("Strange Detector number Missmatch!");
    return kFALSE;
  }

  // TRD space point transformation
  fTransform->SetDetector(det);

  Int_t    iGeoLayer  = AliGeomManager::kTRD1 + ilayer;
  Int_t    iGeoModule = istack + AliTRDgeometry::Nstack() * isector;
  UShort_t volid      = AliGeomManager::LayerToVolUID(iGeoLayer,iGeoModule); 

  Int_t nColMax    = digitsIn->GetNcol();
  Int_t nRowMax    = digitsIn->GetNrow();
  Int_t nTimeTotal = digitsIn->GetNtime();

  // Detector wise calibration object for the gain factors
  const AliTRDCalDet           *calGainFactorDet      = calibration->GetGainFactorDet();
  // Calibration object with pad wise values for the gain factors
  AliTRDCalROC                 *calGainFactorROC      = calibration->GetGainFactorROC(idet);
  // Calibration value for chamber wise gain factor
  Float_t                       calGainFactorDetValue = calGainFactorDet->GetValue(idet);

  // Detector wise calibration object for the noise
  const AliTRDCalDet           *calNoiseDet           = calibration->GetNoiseDet();
  // Calibration object with pad wise values for the noise
  AliTRDCalROC                 *calNoiseROC           = calibration->GetNoiseROC(idet);
  // Calibration value for chamber wise noise
  Float_t                       calNoiseDetValue      = calNoiseDet->GetValue(idet);

  Int_t nClusters = 0;

  AliTRDarraySignal *digitsOut = new AliTRDarraySignal(nRowMax, nColMax, nTimeTotal);
  AliTRDarrayADC     padStatus(nRowMax, nColMax, nTimeTotal);

  ResetHelperIndexes(indexesIn);

  // Apply the gain and the tail cancelation via digital filter
  TailCancelation(digitsIn
                 ,digitsOut  
                 ,indexesIn
                 ,fIndexesOut
                 ,nTimeTotal
                 ,adcThreshold
                 ,calGainFactorROC
                 ,calGainFactorDetValue);	
  
  Int_t row  = 0;
  Int_t col  = 0;
  Int_t time = 0;
  Int_t iPad = 0;
    
  UChar_t status[3]={0, 0, 0}, ipos = 0;
  fIndexesOut->ResetCounters();
  while (fIndexesOut->NextRCTbinIndex(row, col, time)) {
    // reset pad status
    ipos = 0; for(Int_t is=3; is--;) status[is] = 0;

    Float_t signalM = TMath::Abs(digitsOut->GetData(row,col,time));
    status[1] = digitsIn->GetPadStatus(row,col,time);
    if(status[1]) SETBIT(ipos, AliTRDcluster::kMaskedCenter);

    if(signalM < maxThresh) continue; 

    Float_t  noiseMiddleThresh = minMaxCutSigma*calNoiseDetValue*calNoiseROC->GetValue(col,row);
    if (signalM < noiseMiddleThresh) continue;

    if (col + 1 >= nColMax || col-1 < 0) continue;
    
    Float_t signalL = TMath::Abs(digitsOut->GetData(row,col+1,time));
    status[0] = digitsIn->GetPadStatus(row,col+1,time);
    if(status[0]) SETBIT(ipos, AliTRDcluster::kMaskedLeft);
    
    Float_t signalR = TMath::Abs(digitsOut->GetData(row,col-1,time));
    status[2] = digitsIn->GetPadStatus(row,col-1,time);
    if(status[2]) SETBIT(ipos, AliTRDcluster::kMaskedRight);
    
    // reject candidates with more than 1 problematic pad
    if(ipos == 3 || ipos > 4) continue;
    
    if (!status[1]) { // good central pad
      if (!ipos) { // all pads are OK
        if ((signalL <= signalM) && (signalR <  signalM)) {
          if ((signalL >= sigThresh) || (signalR >= sigThresh)) {
             Float_t  noiseSumThresh = minLeftRightCutSigma
                                     * calNoiseDetValue
                                     * calNoiseROC->GetValue(col,row);
            if ((signalL+signalR+signalM) >= noiseSumThresh) {
              // Maximum found, mark the position by a negative signal
              digitsOut->SetData(row,col,time,-signalM);
              fIndexesMaxima->AddIndexTBin(row,col,time);
              padStatus.SetData(row, col, time, ipos); // No corruption
            }
          }
        }
      } else { // one of the neighbouring pads are bad
        if (status[0] && signalR < signalM && signalR >= sigThresh) {
          digitsOut->SetData(row,col,time,-signalM);
          digitsOut->SetData(row, col, time+1, 0.);
          fIndexesMaxima->AddIndexTBin(row,col,time);
          SetPadStatus(status[0], ipos);
          padStatus.SetData(row, col, time, ipos);
        } 
        else if (status[2] && signalL <= signalM && signalL >= sigThresh) {
          digitsOut->SetData(row,col,time,-signalM);
          digitsOut->SetData(row, col, time-1, 0.);
          fIndexesMaxima->AddIndexTBin(row,col,time);
          SetPadStatus(status[2], ipos);
          padStatus.SetData(row, col, time, ipos);
        }
      }
    } 
    else { // wrong maximum pad
      if ((signalL >= sigThresh) || (signalR >= sigThresh)) {
        // Maximum found, mark the position by a negative signal
        digitsOut->SetData(row,col,time,-maxThresh);
        fIndexesMaxima->AddIndexTBin(row,col,time);
        SetPadStatus(status[1], ipos);
        padStatus.SetData(row, col, time, ipos);
      }
    }
  }

  // The index to the first cluster of a given ROC
  Int_t firstClusterROC = -1;
  // The number of cluster in a given ROC
  Int_t nClusterROC     =  0;

  // Now check the maxima and calculate the cluster position
  fIndexesMaxima->ResetCounters();
  while (fIndexesMaxima->NextRCTbinIndex(row, col, time)) {

    // Maximum found ?             
    if (digitsOut->GetData(row,col,time) < 0.0) {

      for (iPad = 0; iPad < kNclus; iPad++) {
        Int_t iPadCol = col - 1 + iPad;
        clusterSignal[iPad] = TMath::Abs(digitsOut->GetData(row,iPadCol,time));
      }

      // Count the number of pads in the cluster
      Int_t nPadCount = 0;
      Int_t ii;
      // Look to the right
      ii = 0;
      while (TMath::Abs(digitsOut->GetData(row,col-ii  ,time)) >= sigThresh) {
        nPadCount++;
        ii++;
        if (col-ii   <        0) break;
      }
      // Look to the left
      ii = 0;
      while (TMath::Abs(digitsOut->GetData(row,col+ii+1,time)) >= sigThresh) {
        nPadCount++;
        ii++;
        if (col+ii+1 >= nColMax) break;
      }
      nClusters++;

      // Look for 5 pad cluster with minimum in the middle
      Bool_t fivePadCluster = kFALSE;
      if (col < (nColMax - 3)){
        if (digitsOut->GetData(row,col+2,time) < 0) {
          fivePadCluster = kTRUE;
        }
        if ((fivePadCluster) && (col < (nColMax - 5))) {
          if (digitsOut->GetData(row,col+4,time) >= sigThresh) {
            fivePadCluster = kFALSE;
          }
        }
        if ((fivePadCluster) && (col >             1)) {
          if (digitsOut->GetData(row,col-2,time) >= sigThresh) {
            fivePadCluster = kFALSE;
          }
        }
      }

      // 5 pad cluster
      // Modify the signal of the overlapping pad for the left part 
      // of the cluster which remains from a previous unfolding
      if (iUnfold) {
        clusterSignal[0] *= ratioLeft;
        iUnfold = 0;
      }

      // Unfold the 5 pad cluster
      if (fivePadCluster) {
        for (iPad = 0; iPad < kNsig; iPad++) {
          padSignal[iPad] = TMath::Abs(digitsOut->GetData(row
                                                         ,col-1+iPad
                                                         ,time));
        }
        // Unfold the two maxima and set the signal on 
        // the overlapping pad to the ratio
        ratioRight        = Unfold(kEpsilon,ilayer,padSignal);
        ratioLeft         = 1.0 - ratioRight; 
        clusterSignal[2] *= ratioRight;
        iUnfold = 1;
      }

      // The position of the cluster in COL direction relative to the center pad (pad units)
      Double_t clusterPosCol = 0.0;
      if (fReconstructor->GetRecoParam()->IsLUT()) {
        // Calculate the position of the cluster by using the
        // lookup table method
        clusterPosCol = LUTposition(ilayer,clusterSignal[0]
                                          ,clusterSignal[1]
                                          ,clusterSignal[2]);
      } 
      else {
        // Calculate the position of the cluster by using the
        // center of gravity method
        for (Int_t i = 0; i < kNsig; i++) {
          padSignal[i] = 0.0;
        }
        padSignal[2] = TMath::Abs(digitsOut->GetData(row,col  ,time)); // Central pad
        padSignal[3] = TMath::Abs(digitsOut->GetData(row,col+1,time)); // Left    pad
        padSignal[1] = TMath::Abs(digitsOut->GetData(row,col-1,time)); // Right   pad
        if ((col >           2) && 
            (TMath::Abs(digitsOut->GetData(row,col-2,time)) < padSignal[1])) {
	  padSignal[4] = TMath::Abs(digitsOut->GetData(row,col-2,time));
        }
        if ((col < nColMax - 3) &&
            (TMath::Abs(digitsOut->GetData(row,col+2,time)) < padSignal[3])) {
          padSignal[0] = TMath::Abs(digitsOut->GetData(row,col+2,time));
        }
        clusterPosCol = GetCOG(padSignal);
      }

      // Store the amplitudes of the pads in the cluster for later analysis
      // and check whether one of these pads is masked in the database
      Short_t signals[7] = { 0, 0, 0, 0, 0, 0, 0 };
      for (Int_t jPad = col-3; jPad <= col+3; jPad++) {
        if ((jPad <          0) || 
            (jPad >= nColMax-1)) {
          continue;
        }
        signals[jPad-col+3] = TMath::Nint(TMath::Abs(digitsOut->GetData(row,jPad,time)));
      }

      // Transform the local cluster coordinates into calibrated 
      // space point positions defined in the local tracking system.
      // Here the calibration for T0, Vdrift and ExB is applied as well.
      Double_t clusterXYZ[6];
      clusterXYZ[0] = clusterPosCol;
      clusterXYZ[1] = clusterSignal[2];
      clusterXYZ[2] = clusterSignal[1];
      clusterXYZ[3] = clusterSignal[0];
      clusterXYZ[4] = 0.0;
      clusterXYZ[5] = 0.0;
      Int_t    clusterRCT[3];
      clusterRCT[0] = row;
      clusterRCT[1] = col;
      clusterRCT[2] = 0;

      Bool_t out = kTRUE;
      if (fTransform->Transform(clusterXYZ,clusterRCT,((UInt_t) time),out,0)) {

        // Add the cluster to the output array
        // The track indices will be stored later 
        Float_t clusterPos[3];
        clusterPos[0] = clusterXYZ[0];
        clusterPos[1] = clusterXYZ[1];
        clusterPos[2] = clusterXYZ[2];
        Float_t clusterSig[2];
        clusterSig[0] = clusterXYZ[4];
        clusterSig[1] = clusterXYZ[5];
        Double_t clusterCharge  = clusterXYZ[3];
        Char_t   clusterTimeBin = ((Char_t) clusterRCT[2]);

        Int_t n = RecPoints()->GetEntriesFast();
        AliTRDcluster *cluster = new((*RecPoints())[n]) AliTRDcluster(
          idet,
          clusterCharge, clusterPos, clusterSig,
          0x0,
          ((Char_t) nPadCount),
          signals,
          ((UChar_t) col), ((UChar_t) row), ((UChar_t) time),
          clusterTimeBin, clusterPosCol,
          volid);
        cluster->SetInChamber(!out);

        UChar_t maskPosition = GetCorruption(padStatus.GetData(row, col, time));
        UChar_t padstatus = GetPadStatus(padStatus.GetData(row, col, time));
        if (maskPosition) { 
          cluster->SetPadMaskedPosition(maskPosition);
          cluster->SetPadMaskedStatus(padstatus);
        }

        // Temporarily store the row, column and time bin of the center pad
        // Used to later on assign the track indices
        cluster->SetLabel( row,0);
        cluster->SetLabel( col,1);
        cluster->SetLabel(time,2);

        // Store the index of the first cluster in the current ROC
        if (firstClusterROC < 0) {
          firstClusterROC = RecPoints()->GetEntriesFast() - 1;
        }

        // Count the number of cluster in the current ROC
        nClusterROC++;

      } // if: Transform ok ?

    } // if: Maximum found ?

  }

  delete digitsOut;

  if (fAddLabels) {
    AddLabels(idet,firstClusterROC,nClusterROC);
  }

  return kTRUE;

}

//_____________________________________________________________________________
Bool_t AliTRDclusterizer::AddLabels(Int_t idet, Int_t firstClusterROC, Int_t nClusterROC)
{
  //
  // Add the track indices to the found clusters
  //
  
  const Int_t   kNclus  = 3;  
  const Int_t   kNdict  = AliTRDdigitsManager::kNDict;
  const Int_t   kNtrack = kNdict * kNclus;

  Int_t iClusterROC = 0;

  Int_t row  = 0;
  Int_t col  = 0;
  Int_t time = 0;
  Int_t iPad = 0;

  // Temporary array to collect the track indices
  Int_t *idxTracks = new Int_t[kNtrack*nClusterROC];

  // Loop through the dictionary arrays one-by-one
  // to keep memory consumption low
  AliTRDarrayDictionary *tracksIn = 0;  //mod
  for (Int_t iDict = 0; iDict < kNdict; iDict++) {

    // tracksIn should be expanded beforehand!
    tracksIn = (AliTRDarrayDictionary *) fDigitsManager->GetDictionary(idet,iDict);

    // Loop though the clusters found in this ROC
    for (iClusterROC = 0; iClusterROC < nClusterROC; iClusterROC++) {

      AliTRDcluster *cluster = (AliTRDcluster *)
                               RecPoints()->UncheckedAt(firstClusterROC+iClusterROC);
      row  = cluster->GetLabel(0);
      col  = cluster->GetLabel(1);
      time = cluster->GetLabel(2);

      for (iPad = 0; iPad < kNclus; iPad++) {
        Int_t iPadCol = col - 1 + iPad;
        Int_t index   = tracksIn->GetData(row,iPadCol,time);  //Modification of -1 in Track
        idxTracks[3*iPad+iDict + iClusterROC*kNtrack] = index;     
      }

    }

  }

  // Copy the track indices into the cluster
  // Loop though the clusters found in this ROC
  for (iClusterROC = 0; iClusterROC < nClusterROC; iClusterROC++) {

    AliTRDcluster *cluster = (AliTRDcluster *)
      RecPoints()->UncheckedAt(firstClusterROC+iClusterROC);
    cluster->SetLabel(-9999,0);
    cluster->SetLabel(-9999,1);
    cluster->SetLabel(-9999,2);
  
    cluster->AddTrackIndex(&idxTracks[iClusterROC*kNtrack]);

  }

  delete [] idxTracks;

  return kTRUE;

}

//_____________________________________________________________________________
Double_t AliTRDclusterizer::GetCOG(Double_t signal[5]) const
{
  //
  // Get COG position
  // Used for clusters with more than 3 pads - where LUT not applicable
  //

  Double_t sum = signal[0]
              + signal[1]
              + signal[2] 
              + signal[3]
              + signal[4];

  // ???????????? CBL
  // Go to 3 pad COG ????
  // ???????????? CBL
  Double_t res = (0.0 * (-signal[0] + signal[4])
                      + (-signal[1] + signal[3])) / sum;

  return res;		  

}

//_____________________________________________________________________________
Double_t AliTRDclusterizer::Unfold(Double_t eps, Int_t layer, Double_t *padSignal)
{
  //
  // Method to unfold neighbouring maxima.
  // The charge ratio on the overlapping pad is calculated
  // until there is no more change within the range given by eps.
  // The resulting ratio is then returned to the calling method.
  //

  AliTRDcalibDB *calibration = AliTRDcalibDB::Instance();
  if (!calibration) {
    AliError("No AliTRDcalibDB instance available\n");
    return kFALSE;  
  }
  
  Int_t   irc                = 0;
  Int_t   itStep             = 0;                 // Count iteration steps

  Double_t ratio             = 0.5;               // Start value for ratio
  Double_t prevRatio         = 0.0;               // Store previous ratio

  Double_t newLeftSignal[3]  = { 0.0, 0.0, 0.0 }; // Array to store left cluster signal
  Double_t newRightSignal[3] = { 0.0, 0.0, 0.0 }; // Array to store right cluster signal
  Double_t newSignal[3]      = { 0.0, 0.0, 0.0 };

  // Start the iteration
  while ((TMath::Abs(prevRatio - ratio) > eps) && (itStep < 10)) {

    itStep++;
    prevRatio = ratio;

    // Cluster position according to charge ratio
    Double_t maxLeft  = (ratio*padSignal[2] - padSignal[0]) 
                      / (padSignal[0] + padSignal[1] + ratio*padSignal[2]);
    Double_t maxRight = (padSignal[4] - (1-ratio)*padSignal[2]) 
                      / ((1.0 - ratio)*padSignal[2] + padSignal[3] + padSignal[4]);

    // Set cluster charge ratio
    irc = calibration->PadResponse(1.0,maxLeft ,layer,newSignal);
    Double_t ampLeft  = padSignal[1] / newSignal[1];
    irc = calibration->PadResponse(1.0,maxRight,layer,newSignal);
    Double_t ampRight = padSignal[3] / newSignal[1];

    // Apply pad response to parameters
    irc = calibration->PadResponse(ampLeft ,maxLeft ,layer,newLeftSignal );
    irc = calibration->PadResponse(ampRight,maxRight,layer,newRightSignal);

    // Calculate new overlapping ratio
    ratio = TMath::Min((Double_t) 1.0
                      ,newLeftSignal[2] / (newLeftSignal[2] + newRightSignal[0]));

  }

  return ratio;

}

//_____________________________________________________________________________
void AliTRDclusterizer::TailCancelation(AliTRDarrayADC *digitsIn
		       		      , AliTRDarraySignal *digitsOut   
				      , AliTRDSignalIndex *indexesIn
				      , AliTRDSignalIndex *indexesOut
				      , Int_t nTimeTotal
				      , Float_t adcThreshold
				      , AliTRDCalROC *calGainFactorROC
				      , Float_t calGainFactorDetValue)
{
  //
  // Applies the tail cancelation and gain factors: 
  // Transform digitsIn to digitsOut
  //

  Int_t iRow  = 0;
  Int_t iCol  = 0;
  Int_t iTime = 0;

  Double_t *inADC  = new Double_t[nTimeTotal];  // ADC data before tail cancellation
  Double_t *outADC = new Double_t[nTimeTotal];  // ADC data after tail cancellation
  indexesIn->ResetCounters();
  while (indexesIn->NextRCIndex(iRow, iCol))
    {
      Float_t  calGainFactorROCValue = calGainFactorROC->GetValue(iCol,iRow);
      Double_t gain                  = calGainFactorDetValue 
                                    * calGainFactorROCValue;

      Bool_t corrupted = kFALSE;
      for (iTime = 0; iTime < nTimeTotal; iTime++) 
        {	  
          // Apply gain gain factor
          inADC[iTime]   = digitsIn->GetDataB(iRow,iCol,iTime);
          if (digitsIn->GetPadStatus(iRow, iCol, iTime)) corrupted = kTRUE;
          inADC[iTime]  /= gain;
          outADC[iTime]  = inADC[iTime];
        }
      if (!corrupted)
        {
          // Apply the tail cancelation via the digital filter
          // (only for non-coorupted pads)
          if (fReconstructor->GetRecoParam() ->IsTailCancelation()) 
            {
              DeConvExp(inADC,outADC,nTimeTotal,fReconstructor->GetRecoParam() ->GetTCnexp());
            } 
        }

      indexesIn->ResetTbinCounter();

      while (indexesIn->NextTbinIndex(iTime))
        {
          // Store the amplitude of the digit if above threshold
          if (outADC[iTime] > adcThreshold) 
            {
              digitsOut->SetData(iRow,iCol,iTime,outADC[iTime]);
              indexesOut->AddIndexTBin(iRow,iCol,iTime);
            }	  
        } // while itime

    } // while irow icol
  
  delete [] inADC;
  delete [] outADC;

  return;

}

//_____________________________________________________________________________
void AliTRDclusterizer::DeConvExp(Double_t *source, Double_t *target
                                , Int_t n, Int_t nexp) 
{
  //
  // Tail cancellation by deconvolution for PASA v4 TRF
  //

  Double_t rates[2];
  Double_t coefficients[2];

  // Initialization (coefficient = alpha, rates = lambda)
  Double_t r1 = 1.0;
  Double_t r2 = 1.0;
  Double_t c1 = 0.5;
  Double_t c2 = 0.5;

  if (nexp == 1) {   // 1 Exponentials
    r1 = 1.156;
    r2 = 0.130;
    c1 = 0.066;
    c2 = 0.000;
  }
  if (nexp == 2) {   // 2 Exponentials
    Double_t par[4];
    fReconstructor->GetTCParams(par);
    r1 = par[0];//1.156;
    r2 = par[1];//0.130;
    c1 = par[2];//0.114;
    c2 = par[3];//0.624;
  }

  coefficients[0] = c1;
  coefficients[1] = c2;

  Double_t dt = 0.1;

  rates[0] = TMath::Exp(-dt/(r1));
  rates[1] = TMath::Exp(-dt/(r2));
  
  Int_t i = 0;
  Int_t k = 0;

  Double_t reminder[2];
  Double_t correction = 0.0;
  Double_t result     = 0.0;

  // Attention: computation order is important
  for (k = 0; k < nexp; k++) {
    reminder[k] = 0.0;
  }

  for (i = 0; i < n; i++) {

    result    = (source[i] - correction);    // No rescaling
    target[i] = result;

    for (k = 0; k < nexp; k++) {
      reminder[k] = rates[k] * (reminder[k] + coefficients[k] * result);
    }

    correction = 0.0;
    for (k = 0; k < nexp; k++) {
      correction += reminder[k];
    }

  }

}

//_____________________________________________________________________________
void AliTRDclusterizer::ResetRecPoints() 
{
  //
  // Resets the list of rec points
  //

  if (fRecPoints) {
    fRecPoints->Delete();
    delete fRecPoints;
  }
}

//_____________________________________________________________________________
TClonesArray *AliTRDclusterizer::RecPoints() 
{
  //
  // Returns the list of rec points
  //

  if (!fRecPoints) {
    if(!(fRecPoints = AliTRDReconstructor::GetClusters())){
      // determine number of clusters which has to be allocated
      Float_t nclusters = fReconstructor->GetRecoParam()->GetNClusters();
      if(fReconstructor->IsHLT()) nclusters /= AliTRDgeometry::kNsector;

      fRecPoints = new TClonesArray("AliTRDcluster", Int_t(nclusters));
    }
    //SetClustersOwner(kTRUE);
    AliTRDReconstructor::SetClusters(0x0);
  }
  return fRecPoints;

}

//_____________________________________________________________________________
void AliTRDclusterizer::FillLUT()
{
  //
  // Create the LUT
  //

  const Int_t kNlut = 128;

  fLUTbin = AliTRDgeometry::kNlayer * kNlut;

  // The lookup table from Bogdan
  Float_t lut[AliTRDgeometry::kNlayer][kNlut] = {  
    {
      0.0070, 0.0150, 0.0224, 0.0298, 0.0374, 0.0454, 0.0533, 0.0611, 
      0.0684, 0.0755, 0.0827, 0.0900, 0.0975, 0.1049, 0.1120, 0.1187, 
      0.1253, 0.1318, 0.1385, 0.1453, 0.1519, 0.1584, 0.1646, 0.1704, 
      0.1762, 0.1821, 0.1879, 0.1938, 0.1996, 0.2053, 0.2108, 0.2160, 
      0.2210, 0.2260, 0.2310, 0.2361, 0.2411, 0.2461, 0.2509, 0.2557, 
      0.2602, 0.2646, 0.2689, 0.2732, 0.2774, 0.2816, 0.2859, 0.2901, 
      0.2942, 0.2983, 0.3022, 0.3061, 0.3099, 0.3136, 0.3172, 0.3207, 
      0.3242, 0.3278, 0.3312, 0.3347, 0.3382, 0.3416, 0.3450, 0.3483, 
      0.3515, 0.3547, 0.3579, 0.3609, 0.3639, 0.3669, 0.3698, 0.3727, 
      0.3756, 0.3785, 0.3813, 0.3842, 0.3870, 0.3898, 0.3926, 0.3952, 
      0.3979, 0.4005, 0.4032, 0.4057, 0.4082, 0.4108, 0.4132, 0.4157, 
      0.4181, 0.4205, 0.4228, 0.4252, 0.4275, 0.4299, 0.4322, 0.4345, 
      0.4367, 0.4390, 0.4412, 0.4434, 0.4456, 0.4478, 0.4499, 0.4520, 
      0.4541, 0.4562, 0.4583, 0.4603, 0.4623, 0.4643, 0.4663, 0.4683, 
      0.4702, 0.4722, 0.4741, 0.4758, 0.4774, 0.4790, 0.4805, 0.4824, 
      0.4844, 0.4863, 0.4883, 0.4902, 0.4921, 0.4940, 0.4959, 0.4978 
    },
    {
      0.0072, 0.0156, 0.0235, 0.0313, 0.0394, 0.0478, 0.0561, 0.0642, 
      0.0718, 0.0792, 0.0868, 0.0947, 0.1025, 0.1101, 0.1172, 0.1241, 
      0.1309, 0.1378, 0.1449, 0.1518, 0.1586, 0.1650, 0.1710, 0.1770, 
      0.1830, 0.1891, 0.1952, 0.2011, 0.2070, 0.2125, 0.2177, 0.2229, 
      0.2280, 0.2332, 0.2383, 0.2435, 0.2484, 0.2533, 0.2581, 0.2627, 
      0.2670, 0.2714, 0.2757, 0.2799, 0.2842, 0.2884, 0.2927, 0.2968, 
      0.3008, 0.3048, 0.3086, 0.3123, 0.3159, 0.3195, 0.3231, 0.3266, 
      0.3301, 0.3335, 0.3370, 0.3404, 0.3438, 0.3471, 0.3504, 0.3536, 
      0.3567, 0.3598, 0.3628, 0.3657, 0.3686, 0.3715, 0.3744, 0.3772, 
      0.3800, 0.3828, 0.3856, 0.3884, 0.3911, 0.3938, 0.3965, 0.3991, 
      0.4016, 0.4042, 0.4067, 0.4092, 0.4116, 0.4140, 0.4164, 0.4187, 
      0.4211, 0.4234, 0.4257, 0.4280, 0.4302, 0.4325, 0.4347, 0.4369, 
      0.4391, 0.4413, 0.4434, 0.4456, 0.4477, 0.4497, 0.4518, 0.4538, 
      0.4558, 0.4578, 0.4598, 0.4618, 0.4637, 0.4656, 0.4675, 0.4694, 
      0.4713, 0.4732, 0.4750, 0.4766, 0.4781, 0.4797, 0.4813, 0.4832, 
      0.4851, 0.4870, 0.4888, 0.4906, 0.4925, 0.4942, 0.4960, 0.4978
    },
    {
      0.0075, 0.0163, 0.0246, 0.0328, 0.0415, 0.0504, 0.0592, 0.0674, 
      0.0753, 0.0832, 0.0914, 0.0996, 0.1077, 0.1154, 0.1225, 0.1296, 
      0.1369, 0.1442, 0.1515, 0.1585, 0.1652, 0.1714, 0.1776, 0.1839, 
      0.1902, 0.1965, 0.2025, 0.2085, 0.2141, 0.2194, 0.2247, 0.2299, 
      0.2352, 0.2405, 0.2457, 0.2507, 0.2557, 0.2604, 0.2649, 0.2693, 
      0.2737, 0.2780, 0.2823, 0.2867, 0.2909, 0.2951, 0.2992, 0.3033, 
      0.3072, 0.3110, 0.3146, 0.3182, 0.3218, 0.3253, 0.3288, 0.3323, 
      0.3357, 0.3392, 0.3426, 0.3459, 0.3492, 0.3524, 0.3555, 0.3586, 
      0.3616, 0.3645, 0.3674, 0.3703, 0.3731, 0.3759, 0.3787, 0.3815, 
      0.3843, 0.3870, 0.3897, 0.3925, 0.3950, 0.3976, 0.4002, 0.4027, 
      0.4052, 0.4076, 0.4101, 0.4124, 0.4148, 0.4171, 0.4194, 0.4217, 
      0.4239, 0.4262, 0.4284, 0.4306, 0.4328, 0.4350, 0.4371, 0.4393, 
      0.4414, 0.4435, 0.4455, 0.4476, 0.4496, 0.4516, 0.4536, 0.4555, 
      0.4575, 0.4594, 0.4613, 0.4632, 0.4650, 0.4669, 0.4687, 0.4705, 
      0.4723, 0.4741, 0.4758, 0.4773, 0.4789, 0.4804, 0.4821, 0.4839, 
      0.4857, 0.4875, 0.4893, 0.4910, 0.4928, 0.4945, 0.4961, 0.4978
    },
    {
      0.0078, 0.0171, 0.0258, 0.0345, 0.0438, 0.0532, 0.0624, 0.0708, 
      0.0791, 0.0875, 0.0962, 0.1048, 0.1130, 0.1206, 0.1281, 0.1356, 
      0.1432, 0.1508, 0.1582, 0.1651, 0.1716, 0.1780, 0.1845, 0.1910, 
      0.1975, 0.2038, 0.2099, 0.2155, 0.2210, 0.2263, 0.2317, 0.2371, 
      0.2425, 0.2477, 0.2528, 0.2578, 0.2626, 0.2671, 0.2715, 0.2759, 
      0.2803, 0.2846, 0.2890, 0.2933, 0.2975, 0.3016, 0.3056, 0.3095, 
      0.3132, 0.3168, 0.3204, 0.3239, 0.3274, 0.3309, 0.3344, 0.3378, 
      0.3412, 0.3446, 0.3479, 0.3511, 0.3543, 0.3574, 0.3603, 0.3633, 
      0.3662, 0.3690, 0.3718, 0.3747, 0.3774, 0.3802, 0.3829, 0.3857, 
      0.3883, 0.3910, 0.3936, 0.3962, 0.3987, 0.4012, 0.4037, 0.4061, 
      0.4085, 0.4109, 0.4132, 0.4155, 0.4177, 0.4200, 0.4222, 0.4244, 
      0.4266, 0.4288, 0.4309, 0.4331, 0.4352, 0.4373, 0.4394, 0.4414, 
      0.4435, 0.4455, 0.4475, 0.4494, 0.4514, 0.4533, 0.4552, 0.4571, 
      0.4590, 0.4608, 0.4626, 0.4645, 0.4662, 0.4680, 0.4698, 0.4715, 
      0.4733, 0.4750, 0.4766, 0.4781, 0.4796, 0.4812, 0.4829, 0.4846, 
      0.4863, 0.4880, 0.4897, 0.4914, 0.4930, 0.4946, 0.4963, 0.4979
    },
    {
      0.0081, 0.0178, 0.0270, 0.0364, 0.0463, 0.0562, 0.0656, 0.0744, 
      0.0831, 0.0921, 0.1013, 0.1102, 0.1183, 0.1261, 0.1339, 0.1419, 
      0.1499, 0.1576, 0.1648, 0.1715, 0.1782, 0.1849, 0.1917, 0.1984, 
      0.2048, 0.2110, 0.2167, 0.2223, 0.2278, 0.2333, 0.2389, 0.2444, 
      0.2497, 0.2548, 0.2598, 0.2645, 0.2691, 0.2735, 0.2780, 0.2824, 
      0.2868, 0.2912, 0.2955, 0.2997, 0.3038, 0.3078, 0.3116, 0.3152, 
      0.3188, 0.3224, 0.3259, 0.3294, 0.3329, 0.3364, 0.3398, 0.3432, 
      0.3465, 0.3497, 0.3529, 0.3561, 0.3591, 0.3620, 0.3649, 0.3677, 
      0.3705, 0.3733, 0.3761, 0.3788, 0.3816, 0.3843, 0.3869, 0.3896, 
      0.3922, 0.3948, 0.3973, 0.3998, 0.4022, 0.4047, 0.4070, 0.4094, 
      0.4117, 0.4139, 0.4162, 0.4184, 0.4206, 0.4227, 0.4249, 0.4270, 
      0.4291, 0.4313, 0.4334, 0.4354, 0.4375, 0.4395, 0.4415, 0.4435, 
      0.4455, 0.4474, 0.4493, 0.4512, 0.4531, 0.4550, 0.4568, 0.4586, 
      0.4604, 0.4622, 0.4639, 0.4657, 0.4674, 0.4691, 0.4708, 0.4725, 
      0.4742, 0.4758, 0.4773, 0.4788, 0.4803, 0.4819, 0.4836, 0.4852, 
      0.4869, 0.4885, 0.4901, 0.4917, 0.4933, 0.4948, 0.4964, 0.4979
    },
    {
      0.0085, 0.0189, 0.0288, 0.0389, 0.0497, 0.0603, 0.0699, 0.0792, 
      0.0887, 0.0985, 0.1082, 0.1170, 0.1253, 0.1336, 0.1421, 0.1505, 
      0.1587, 0.1662, 0.1733, 0.1803, 0.1874, 0.1945, 0.2014, 0.2081, 
      0.2143, 0.2201, 0.2259, 0.2316, 0.2374, 0.2431, 0.2487, 0.2541, 
      0.2593, 0.2642, 0.2689, 0.2735, 0.2781, 0.2826, 0.2872, 0.2917, 
      0.2961, 0.3003, 0.3045, 0.3086, 0.3125, 0.3162, 0.3198, 0.3235, 
      0.3270, 0.3306, 0.3342, 0.3377, 0.3411, 0.3446, 0.3479, 0.3511, 
      0.3543, 0.3575, 0.3605, 0.3634, 0.3663, 0.3691, 0.3720, 0.3748, 
      0.3775, 0.3803, 0.3830, 0.3857, 0.3884, 0.3911, 0.3937, 0.3962, 
      0.3987, 0.4012, 0.4036, 0.4060, 0.4084, 0.4107, 0.4129, 0.4152, 
      0.4174, 0.4196, 0.4218, 0.4239, 0.4261, 0.4282, 0.4303, 0.4324, 
      0.4344, 0.4365, 0.4385, 0.4405, 0.4425, 0.4445, 0.4464, 0.4483, 
      0.4502, 0.4521, 0.4539, 0.4558, 0.4576, 0.4593, 0.4611, 0.4629, 
      0.4646, 0.4663, 0.4680, 0.4697, 0.4714, 0.4730, 0.4747, 0.4759, 
      0.4769, 0.4780, 0.4790, 0.4800, 0.4811, 0.4827, 0.4843, 0.4859, 
      0.4874, 0.4889, 0.4905, 0.4920, 0.4935, 0.4950, 0.4965, 0.4979
    }
  }; 

  if (fLUT) {
    delete [] fLUT;
  }
  fLUT = new Double_t[fLUTbin];

  for (Int_t ilayer = 0; ilayer < AliTRDgeometry::kNlayer; ilayer++) {
    for (Int_t ilut  = 0; ilut  <  kNlut; ilut++  ) {
      fLUT[ilayer*kNlut+ilut] = lut[ilayer][ilut];
    }
  }

}

//_____________________________________________________________________________
Double_t AliTRDclusterizer::LUTposition(Int_t ilayer
                                      , Double_t ampL
                                      , Double_t ampC
                                      , Double_t ampR) const
{
  //
  // Calculates the cluster position using the lookup table.
  // Method provided by Bogdan Vulpescu.
  //

  const Int_t kNlut = 128;

  Double_t pos;
  Double_t x    = 0.0;
  Double_t xmin;
  Double_t xmax;
  Double_t xwid;

  Int_t    side = 0;
  Int_t    ix;

  Double_t xMin[AliTRDgeometry::kNlayer] = { 0.006492, 0.006377, 0.006258
                                          , 0.006144, 0.006030, 0.005980 };
  Double_t xMax[AliTRDgeometry::kNlayer] = { 0.960351, 0.965870, 0.970445
                                          , 0.974352, 0.977667, 0.996101 };

  if      (ampL > ampR) {
    x    = (ampL - ampR) / ampC;
    side = -1;
  } 
  else if (ampL < ampR) {
    x    = (ampR - ampL) / ampC;
    side = +1;
  }

  if (ampL != ampR) {

    xmin = xMin[ilayer] + 0.000005;
    xmax = xMax[ilayer] - 0.000005;
    xwid = (xmax - xmin) / 127.0;

    if      (x < xmin) {
      pos = 0.0000;
    } 
    else if (x > xmax) {
      pos = side * 0.5000;
    } 
    else {
      ix  = (Int_t) ((x - xmin) / xwid);
      pos = side * fLUT[ilayer*kNlut+ix];
    }
      
  } 
  else {

    pos = 0.0;

  }

  return pos;

}
