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

#include "TTreeStream.h"

#include "Cal/AliTRDCalROC.h"
#include "Cal/AliTRDCalDet.h"
#include "Cal/AliTRDCalSingleChamberStatus.h"

ClassImp(AliTRDclusterizer)

//_____________________________________________________________________________
AliTRDclusterizer::AliTRDclusterizer(const AliTRDReconstructor *const rec)
  :TNamed()
  ,fReconstructor(rec)  
  ,fRunLoader(NULL)
  ,fClusterTree(NULL)
  ,fRecPoints(NULL)
  ,fTrackletTree(NULL)
  ,fDigitsManager(new AliTRDdigitsManager())
  ,fTrackletContainer(NULL)
  ,fRawVersion(2)
  ,fTransform(new AliTRDtransform(0))
  ,fDigits(NULL)
  ,fIndexes(NULL)
  ,fADCthresh(0)
  ,fMaxThresh(0)
  ,fSigThresh(0)
  ,fMinMaxCutSigma(0)
  ,fMinLeftRightCutSigma(0)
  ,fLayer(0)
  ,fDet(0)
  ,fVolid(0)
  ,fColMax(0)
  ,fTimeTotal(0)
  ,fCalGainFactorROC(NULL)
  ,fCalGainFactorDetValue(0)
  ,fCalNoiseROC(NULL)
  ,fCalNoiseDetValue(0)
  ,fCalPadStatusROC(NULL)
  ,fClusterROC(0)
  ,firstClusterROC(0)
  ,fNoOfClusters(0)
{
  //
  // AliTRDclusterizer default constructor
  //

  SetBit(kLabels, kTRUE);

  AliTRDcalibDB *trd = 0x0;
  if (!(trd = AliTRDcalibDB::Instance())) {
    AliFatal("Could not get calibration object");
  }

  fRawVersion = AliTRDfeeParam::Instance()->GetRAWversion();

  // Initialize debug stream
  if(fReconstructor){
    if(fReconstructor->GetStreamLevel(AliTRDReconstructor::kClusterizer) > 1){
      TDirectory *savedir = gDirectory; 
      //fgGetDebugStream    = new TTreeSRedirector("TRD.ClusterizerDebug.root");
      savedir->cd();
    }
  }

}

//_____________________________________________________________________________
AliTRDclusterizer::AliTRDclusterizer(const Text_t *name, const Text_t *title, const AliTRDReconstructor *const rec)
  :TNamed(name,title)
  ,fReconstructor(rec)
  ,fRunLoader(NULL)
  ,fClusterTree(NULL)
  ,fRecPoints(NULL)
  ,fTrackletTree(NULL)
  ,fDigitsManager(new AliTRDdigitsManager())
  ,fTrackletContainer(NULL)
  ,fRawVersion(2)
  ,fTransform(new AliTRDtransform(0))
  ,fDigits(NULL)
  ,fIndexes(NULL)
  ,fADCthresh(0)
  ,fMaxThresh(0)
  ,fSigThresh(0)
  ,fMinMaxCutSigma(0)
  ,fMinLeftRightCutSigma(0)
  ,fLayer(0)
  ,fDet(0)
  ,fVolid(0)
  ,fColMax(0)
  ,fTimeTotal(0)
  ,fCalGainFactorROC(NULL)
  ,fCalGainFactorDetValue(0)
  ,fCalNoiseROC(NULL)
  ,fCalNoiseDetValue(0)
  ,fCalPadStatusROC(NULL)
  ,fClusterROC(0)
  ,firstClusterROC(0)
  ,fNoOfClusters(0)
{
  //
  // AliTRDclusterizer constructor
  //

  SetBit(kLabels, kTRUE);

  AliTRDcalibDB *trd = 0x0;
  if (!(trd = AliTRDcalibDB::Instance())) {
    AliFatal("Could not get calibration object");
  }

  fDigitsManager->CreateArrays();

  fRawVersion = AliTRDfeeParam::Instance()->GetRAWversion();

  //FillLUT();

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
  ,fRawVersion(2)
  ,fTransform(NULL)
  ,fDigits(NULL)
  ,fIndexes(NULL)
  ,fADCthresh(0)
  ,fMaxThresh(0)
  ,fSigThresh(0)
  ,fMinMaxCutSigma(0)
  ,fMinLeftRightCutSigma(0)
  ,fLayer(0)
  ,fDet(0)
  ,fVolid(0)
  ,fColMax(0)
  ,fTimeTotal(0)
  ,fCalGainFactorROC(NULL)
  ,fCalGainFactorDetValue(0)
  ,fCalNoiseROC(NULL)
  ,fCalNoiseDetValue(0)
  ,fCalPadStatusROC(NULL)
  ,fClusterROC(0)
  ,firstClusterROC(0)
  ,fNoOfClusters(0)
{
  //
  // AliTRDclusterizer copy constructor
  //

  SetBit(kLabels, kTRUE);

  //FillLUT();

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

  if (fTransform){
    delete fTransform;
    fTransform     = NULL;
  }

//   if (fLUT) {
//     delete [] fLUT;
//     fLUT           = NULL;
//   }

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
  ((AliTRDclusterizer &) c).fRawVersion    = fRawVersion;
  ((AliTRDclusterizer &) c).fTransform     = NULL;
  ((AliTRDclusterizer &) c).fDigits      = NULL;
  ((AliTRDclusterizer &) c).fIndexes       = NULL;
  ((AliTRDclusterizer &) c).fADCthresh     = 0;
  ((AliTRDclusterizer &) c).fMaxThresh     = 0;
  ((AliTRDclusterizer &) c).fSigThresh     = 0;
  ((AliTRDclusterizer &) c).fMinMaxCutSigma= 0;
  ((AliTRDclusterizer &) c).fMinLeftRightCutSigma = 0;
  ((AliTRDclusterizer &) c).fLayer         = 0;
  ((AliTRDclusterizer &) c).fDet           = 0;
  ((AliTRDclusterizer &) c).fVolid         = 0;
  ((AliTRDclusterizer &) c).fColMax        = 0;
  ((AliTRDclusterizer &) c).fTimeTotal     = 0;
  ((AliTRDclusterizer &) c).fCalGainFactorROC = NULL;
  ((AliTRDclusterizer &) c).fCalGainFactorDetValue = 0;
  ((AliTRDclusterizer &) c).fCalNoiseROC   = NULL;
  ((AliTRDclusterizer &) c).fCalNoiseDetValue = 0;
  ((AliTRDclusterizer &) c).fCalPadStatusROC = NULL;
  ((AliTRDclusterizer &) c).fClusterROC    = 0;
  ((AliTRDclusterizer &) c).firstClusterROC= 0;
  ((AliTRDclusterizer &) c).fNoOfClusters  = 0;
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
  if (TestBit(kLabels)){
    SetBit(kLabels, fDigitsManager->UsesDictionaries());
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
      if (TestBit(kLabels)){
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
    fDigitsManager = new AliTRDdigitsManager(kTRUE);
    fDigitsManager->CreateArrays();
  }

  fDigitsManager->SetUseDictionaries(TestBit(kLabels));

  // tracklet container for raw tracklet writing
  if (!fTrackletContainer && fReconstructor->IsWritingTracklets()) {
    // maximum tracklets for one HC
    const Int_t kTrackletChmb=256;
    fTrackletContainer = new UInt_t *[2];
    fTrackletContainer[0] = new UInt_t[kTrackletChmb]; 
    fTrackletContainer[1] = new UInt_t[kTrackletChmb]; 
  }

  AliTRDrawStreamBase *input = AliTRDrawStreamBase::GetRawStream(rawReader);

  AliInfo(Form("Stream version: %s", input->IsA()->GetName()));
  
  Int_t det    = 0;
  while ((det = input->NextChamber(fDigitsManager,fTrackletContainer)) >= 0){
    Bool_t iclusterBranch = kFALSE;
    if (fDigitsManager->GetIndexes(det)->HasEntry()){
      iclusterBranch = MakeClusters(det);
    }

    fDigitsManager->ResetArrays(det);
    
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

  delete input;
  input = NULL;

  AliInfo(Form("Number of found clusters : %d", fNoOfClusters)); 
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
void AliTRDclusterizer::SetPadStatus(const UChar_t status, UChar_t &out){
  //
  // Set the pad status into out
  // First three bits are needed for the position encoding
  //
  out |= status << 3;
}

//_____________________________________________________________________________
UChar_t AliTRDclusterizer::GetPadStatus(UChar_t encoding) const {
  //
  // return the staus encoding of the corrupted pad
  //
  return static_cast<UChar_t>(encoding >> 3);
}

//_____________________________________________________________________________
Int_t AliTRDclusterizer::GetCorruption(UChar_t encoding) const {
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
  fDigits = (AliTRDarrayADC *) fDigitsManager->GetDigits(det); //mod     
  
  // This is to take care of switched off super modules
  if (!fDigits->HasData()) return kFALSE;

  fIndexes = fDigitsManager->GetIndexes(det);
  if (fIndexes->IsAllocated() == kFALSE) {
    AliError("Indexes do not exist!");
    return kFALSE;      
  }

  AliTRDcalibDB  *calibration = AliTRDcalibDB::Instance();
  if (!calibration) {
    AliFatal("No AliTRDcalibDB instance available\n");
    return kFALSE;  
  }

  fADCthresh = 0; 

  if (!fReconstructor){
    AliError("Reconstructor not set\n");
    return kFALSE;
  }

  TTreeSRedirector *fDebugStream = fReconstructor->GetDebugStream(AliTRDReconstructor::kClusterizer);

  fMaxThresh            = fReconstructor->GetRecoParam()->GetClusMaxThresh();
  fSigThresh            = fReconstructor->GetRecoParam()->GetClusSigThresh();
  fMinMaxCutSigma       = fReconstructor->GetRecoParam()->GetMinMaxCutSigma();
  fMinLeftRightCutSigma = fReconstructor->GetRecoParam()->GetMinLeftRightCutSigma();

  Int_t istack  = fIndexes->GetStack();
  fLayer  = fIndexes->GetLayer();
  Int_t isector = fIndexes->GetSM();

  // Start clustering in the chamber

  fDet  = AliTRDgeometry::GetDetector(fLayer,istack,isector);
  if (fDet != det) {
    AliError("Strange Detector number Missmatch!");
    return kFALSE;
  }

  // TRD space point transformation
  fTransform->SetDetector(det);

  Int_t    iGeoLayer  = AliGeomManager::kTRD1 + fLayer;
  Int_t    iGeoModule = istack + AliTRDgeometry::Nstack() * isector;
  fVolid      = AliGeomManager::LayerToVolUID(iGeoLayer,iGeoModule); 

  fColMax    = fDigits->GetNcol();
  //Int_t nRowMax    = fDigits->GetNrow();
  fTimeTotal = fDigits->GetNtime();

  // Detector wise calibration object for the gain factors
  const AliTRDCalDet *calGainFactorDet = calibration->GetGainFactorDet();
  // Calibration object with pad wise values for the gain factors
  fCalGainFactorROC      = calibration->GetGainFactorROC(fDet);
  // Calibration value for chamber wise gain factor
  fCalGainFactorDetValue = calGainFactorDet->GetValue(fDet);

  // Detector wise calibration object for the noise
  const AliTRDCalDet *calNoiseDet = calibration->GetNoiseDet();
  // Calibration object with pad wise values for the noise
  fCalNoiseROC           = calibration->GetNoiseROC(fDet);
  // Calibration value for chamber wise noise
  fCalNoiseDetValue      = calNoiseDet->GetValue(fDet);
  
  // Calibration object with the pad status
  fCalPadStatusROC       = calibration->GetPadStatusROC(fDet);
  
  SetBit(kLUT, fReconstructor->UseLUT());
  SetBit(kGAUS, fReconstructor->UseGAUS());
  SetBit(kHLT, fReconstructor->IsHLT());

  firstClusterROC = -1;
  fClusterROC     =  0;

  // Apply the gain and the tail cancelation via digital filter
  if(fReconstructor->UseTailCancelation()) TailCancelation();

  MaxStruct curr, last;
  Int_t nMaximas = 0, nCorrupted = 0;

  // Here the clusterfining is happening
  
  for(curr.Time = 0; curr.Time < fTimeTotal; curr.Time++){
    while(fIndexes->NextRCIndex(curr.Row, curr.Col)){
      //printf("\nCHECK r[%2d] c[%3d] t[%d]\n", curr.Row, curr.Col, curr.Time);
      if(IsMaximum(curr, curr.padStatus, &curr.Signals[0])){
        //printf("\tMAX s[%d %d %d]\n", curr.Signals[0], curr.Signals[1], curr.Signals[2]);
        if(last.Row>-1){
          if(curr.Time==last.Time && curr.Row==last.Row && curr.Col==last.Col+2) FivePadCluster(last, curr);
          CreateCluster(last);
        }
        last=curr; curr.FivePad=kFALSE;
      }
      //printf("\t--- s[%d %d %d]\n", curr.Signals[0], curr.Signals[1], curr.Signals[2]);
    }
  }
  if(last.Row>-1) CreateCluster(last);

  if(fReconstructor->GetStreamLevel(AliTRDReconstructor::kClusterizer) > 2){
    (*fDebugStream) << "MakeClusters"
      << "Detector="   << det
      << "NMaxima="    << nMaximas
      << "NClusters="  << fClusterROC
      << "NCorrupted=" << nCorrupted
      << "\n";
  }
  if (TestBit(kLabels)) AddLabels();

  return kTRUE;

}

//_____________________________________________________________________________
Bool_t AliTRDclusterizer::IsMaximum(const MaxStruct &Max, UChar_t &padStatus, Short_t *const Signals) 
{
  //
  // Returns true if this row,col,time combination is a maximum. 
  // Gives back the padStatus and the signals of the center pad and the two neighbouring pads.
  //

  Signals[1] = fDigits->GetData(Max.Row, Max.Col, Max.Time);
  if(Signals[1] < fMaxThresh) return kFALSE;

  Float_t  noiseMiddleThresh = fMinMaxCutSigma*fCalNoiseDetValue*fCalNoiseROC->GetValue(Max.Col, Max.Row);
  if (Signals[1] < noiseMiddleThresh) return kFALSE;

  if (Max.Col + 1 >= fColMax || Max.Col < 1) return kFALSE;

  UChar_t status[3]={
    fCalPadStatusROC->GetStatus(Max.Col-1, Max.Row)
   ,fCalPadStatusROC->GetStatus(Max.Col,   Max.Row)
   ,fCalPadStatusROC->GetStatus(Max.Col+1, Max.Row)
  };

  Signals[0] = fDigits->GetData(Max.Row, Max.Col-1, Max.Time);
  Signals[2] = fDigits->GetData(Max.Row, Max.Col+1, Max.Time);  

  if(!(status[0] | status[1] | status[2])) {//all pads are good
    if ((Signals[2] <= Signals[1]) && (Signals[0] <  Signals[1])) {
      if ((Signals[2] >= fSigThresh) || (Signals[0] >= fSigThresh)) {
        Float_t  noiseSumThresh = fMinLeftRightCutSigma
          * fCalNoiseDetValue
          * fCalNoiseROC->GetValue(Max.Col, Max.Row);
        if ((Signals[2]+Signals[0]+Signals[1]) < noiseSumThresh) return kFALSE;
        padStatus = 0;
        return kTRUE;
      }
    }
  } else { // at least one of the pads is bad, and reject candidates with more than 1 problematic pad
    if (status[2] && (!(status[0] || status[1])) && Signals[1] > Signals[0] && Signals[0] >= fSigThresh) { 
      Signals[2]=0;
      SetPadStatus(status[2], padStatus);
      return kTRUE;
    } 
    else if (status[0] && (!(status[1] || status[2])) && Signals[1] >= Signals[2] && Signals[2] >= fSigThresh) {
      Signals[0]=0;
      SetPadStatus(status[0], padStatus);
      return kTRUE;
    }
    else if (status[1] && (!(status[0] || status[2])) && ((Signals[2] >= fSigThresh) || (Signals[0] >= fSigThresh))) {
      Signals[1]=TMath::Nint(fMaxThresh);
      SetPadStatus(status[1], padStatus);
      return kTRUE;
    }
  }
  return kFALSE;
}

//_____________________________________________________________________________
Bool_t AliTRDclusterizer::FivePadCluster(MaxStruct &ThisMax, MaxStruct &NeighbourMax)
{
  //
  // Look for 5 pad cluster with minimum in the middle
  // Gives back the ratio
  //
  if (ThisMax.Col >= fColMax - 3) return kFALSE;
  if (ThisMax.Col < fColMax - 5){
    if (fDigits->GetData(ThisMax.Row, ThisMax.Col+4, ThisMax.Time) >= fSigThresh)
      return kFALSE;
  }
  if (ThisMax.Col > 1) {
    if (fDigits->GetData(ThisMax.Row, ThisMax.Col-2, ThisMax.Time) >= fSigThresh)
      return kFALSE;
  }
  
  const Float_t kEpsilon = 0.01;
  Double_t padSignal[5] = {ThisMax.Signals[0], ThisMax.Signals[1], ThisMax.Signals[2],
      NeighbourMax.Signals[1], NeighbourMax.Signals[2]};
  
  // Unfold the two maxima and set the signal on 
  // the overlapping pad to the ratio
  Float_t ratio = Unfold(kEpsilon,fLayer,padSignal);
  ThisMax.Signals[2] = TMath::Nint(ThisMax.Signals[2]*ratio);
  NeighbourMax.Signals[0] = TMath::Nint(NeighbourMax.Signals[0]*(1-ratio));
  ThisMax.FivePad=kTRUE;
  NeighbourMax.FivePad=kTRUE;
  return kTRUE;

}

//_____________________________________________________________________________
void AliTRDclusterizer::CreateCluster(const MaxStruct &Max)
{
  //
  // Creates a cluster at the given position and saves it in fRecPoints
  //

  Int_t nPadCount = 1;
  Short_t signals[7] = { 0, 0, Max.Signals[0], Max.Signals[1], Max.Signals[2], 0, 0 };
  if(!TestBit(kHLT)) CalcAdditionalInfo(Max, signals, nPadCount);

  AliTRDcluster cluster(fDet, ((UChar_t) Max.Col), ((UChar_t) Max.Row), ((UChar_t) Max.Time), signals, fVolid);
  cluster.SetNPads(nPadCount);
  if(TestBit(kLUT)) cluster.SetRPhiMethod(AliTRDcluster::kLUT);
  else if(TestBit(kGAUS)) cluster.SetRPhiMethod(AliTRDcluster::kGAUS);
  else cluster.SetRPhiMethod(AliTRDcluster::kCOG);

  cluster.SetFivePad(Max.FivePad);
  // set pads status for the cluster
  UChar_t maskPosition = GetCorruption(Max.padStatus);
  if (maskPosition) { 
    cluster.SetPadMaskedPosition(maskPosition);
    cluster.SetPadMaskedStatus(GetPadStatus(Max.padStatus));
  }

  // Transform the local cluster coordinates into calibrated 
  // space point positions defined in the local tracking system.
  // Here the calibration for T0, Vdrift and ExB is applied as well.
  if(!fTransform->Transform(&cluster)) return;
  // Temporarily store the Max.Row, column and time bin of the center pad
  // Used to later on assign the track indices
  cluster.SetLabel(Max.Row, 0);
  cluster.SetLabel(Max.Col, 1);
  cluster.SetLabel(Max.Time,2);

  //needed for HLT reconstruction
  AddClusterToArray(&cluster); 

  // Store the index of the first cluster in the current ROC
  if (firstClusterROC < 0) firstClusterROC = fNoOfClusters;
  
  fNoOfClusters++;
  fClusterROC++;
}

//_____________________________________________________________________________
void AliTRDclusterizer::CalcAdditionalInfo(const MaxStruct &Max, Short_t *const signals, Int_t &nPadCount)
{
  // Look to the right
  Int_t ii = 1;
  while (fDigits->GetData(Max.Row, Max.Col-ii, Max.Time) >= fSigThresh) {
    nPadCount++;
    ii++;
    if (Max.Col < ii) break;
  }
  // Look to the left
  ii = 1;
  while (fDigits->GetData(Max.Row, Max.Col+ii, Max.Time) >= fSigThresh) {
    nPadCount++;
    ii++;
    if (Max.Col+ii >= fColMax) break;
  }

  // Store the amplitudes of the pads in the cluster for later analysis
  // and check whether one of these pads is masked in the database
  signals[2]=Max.Signals[0];
  signals[3]=Max.Signals[1];
  signals[4]=Max.Signals[2];
  for(Int_t i = 0; i<2; i++)
    {
      if(Max.Col+i >= 3)
	signals[i] = fDigits->GetData(Max.Row, Max.Col-3+i, Max.Time);
      if(Max.Col+3-i < fColMax)
	signals[6-i] = fDigits->GetData(Max.Row, Max.Col+3-i, Max.Time);
    }
  /*for (Int_t jPad = Max.Col-3; jPad <= Max.Col+3; jPad++) {
    if ((jPad >= 0) && (jPad < fColMax))
      signals[jPad-Max.Col+3] = TMath::Nint(fDigits->GetData(Max.Row,jPad,Max.Time));
      }*/
}

//_____________________________________________________________________________
void AliTRDclusterizer::AddClusterToArray(AliTRDcluster *cluster)
{
  //
  // Add a cluster to the array
  //

  Int_t n = RecPoints()->GetEntriesFast();
  if(n!=fNoOfClusters)AliError(Form("fNoOfClusters != RecPoints()->GetEntriesFast %i != %i \n", fNoOfClusters, n));
  new((*RecPoints())[n]) AliTRDcluster(*cluster);
}

//_____________________________________________________________________________
Bool_t AliTRDclusterizer::AddLabels()
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
  Int_t *idxTracks = new Int_t[kNtrack*fClusterROC];

  // Loop through the dictionary arrays one-by-one
  // to keep memory consumption low
  AliTRDarrayDictionary *tracksIn = 0;  //mod
  for (Int_t iDict = 0; iDict < kNdict; iDict++) {

    // tracksIn should be expanded beforehand!
    tracksIn = (AliTRDarrayDictionary *) fDigitsManager->GetDictionary(fDet,iDict);

    // Loop though the clusters found in this ROC
    for (iClusterROC = 0; iClusterROC < fClusterROC; iClusterROC++) {

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
  for (iClusterROC = 0; iClusterROC < fClusterROC; iClusterROC++) {

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
Float_t AliTRDclusterizer::Unfold(Double_t eps, Int_t layer, Double_t *padSignal) const
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
                      / (padSignal[0] + padSignal[1] + ratio * padSignal[2]);
    Double_t maxRight = (padSignal[4] - (1-ratio)*padSignal[2]) 
                      / ((1.0 - ratio)*padSignal[2] + padSignal[3] + padSignal[4]);

    // Set cluster charge ratio
    irc = calibration->PadResponse(1.0, maxLeft, layer, newSignal);
    Double_t ampLeft  = padSignal[1] / newSignal[1];
    irc = calibration->PadResponse(1.0, maxRight, layer, newSignal);
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
void AliTRDclusterizer::TailCancelation()
{
  //
  // Applies the tail cancelation and gain factors: 
  // Transform fDigits to fDigits
  //

  Int_t iRow  = 0;
  Int_t iCol  = 0;
  Int_t iTime = 0;

  Double_t *inADC = new Double_t[fTimeTotal];  // ADC data before tail cancellation
  Double_t *outADC = new Double_t[fTimeTotal];  // ADC data after tail cancellation

  fIndexes->ResetCounters();
  TTreeSRedirector *fDebugStream = fReconstructor->GetDebugStream(AliTRDReconstructor::kClusterizer);
  while(fIndexes->NextRCIndex(iRow, iCol))
    {
      Float_t  fCalGainFactorROCValue = fCalGainFactorROC->GetValue(iCol,iRow);
      Double_t gain                  = fCalGainFactorDetValue 
                                    * fCalGainFactorROCValue;

      Bool_t corrupted = kFALSE;
      for (iTime = 0; iTime < fTimeTotal; iTime++) 
        {	  
          // Apply gain gain factor
          inADC[iTime]   = fDigits->GetData(iRow,iCol,iTime);
          if (fCalPadStatusROC->GetStatus(iCol, iRow)) corrupted = kTRUE;
          inADC[iTime]  /= gain;
          outADC[iTime]  = inADC[iTime];
      	  if(fReconstructor->GetStreamLevel(AliTRDReconstructor::kClusterizer) > 7){
      	    (*fDebugStream) << "TailCancellation"
  			      << "col="  << iCol
  			      << "row="  << iRow
  			      << "time=" << iTime
  			      << "inADC=" << inADC[iTime]
  			      << "gain=" << gain
  			      << "outADC=" << outADC[iTime]
  			      << "corrupted=" << corrupted
  			      << "\n";
      	  }
        }
      if (!corrupted)
        {
          // Apply the tail cancelation via the digital filter
          // (only for non-coorupted pads)
	  DeConvExp(&inADC[0],&outADC[0],fTimeTotal,fReconstructor->GetRecoParam() ->GetTCnexp());
        }

      for(iTime = 0; iTime < fTimeTotal; iTime++)//while (fIndexes->NextTbinIndex(iTime))
        {
          // Store the amplitude of the digit if above threshold
          if (outADC[iTime] > fADCthresh)
	    fDigits->SetData(iRow,iCol,iTime,TMath::Nint(outADC[iTime]));
	  else
	    fDigits->SetData(iRow,iCol,iTime,0);
        } // while itime

    } // while irow icol

  delete [] inADC;
  delete [] outADC;

  return;

}

//_____________________________________________________________________________
void AliTRDclusterizer::DeConvExp(const Double_t *const source, Double_t *const target
				  ,const Int_t n, const Int_t nexp) 
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

      fRecPoints = new TClonesArray("AliTRDcluster", Int_t(nclusters));
    }
    //SetClustersOwner(kTRUE);
    AliTRDReconstructor::SetClusters(0x0);
  }
  return fRecPoints;

}

