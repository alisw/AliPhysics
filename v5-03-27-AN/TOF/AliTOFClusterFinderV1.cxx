/***************************************************************************
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

/*
$Log: AliTOFClusterFinderV1.cxx,v $
Revision 2009/04/20 A. De Caro
    - added two new global variables, called fTOFGeometry and fTOFdigits;
    - added a new method, called FindClustersWithoutTOT,
      to transform TOF digits with fTOT=0 in one pad clusters;
    - update of the covariance matrix elements for the TOF clusters

Revision 0.01  2008/05/10 A. De Caro
 */

/////////////////////////////////////////
//                                     //
//  Class for TOF cluster finder (V1)  //
//                                     //
//  Input data: Raw Data or Digits;    //
//  Output data: Digits or Rec Points  //
//                                     //
/////////////////////////////////////////

#include "Riostream.h"

#include "TClonesArray.h"
#include "TStopwatch.h"
#include "TTree.h"
#include "TGeoMatrix.h"
#include "TString.h"

#include "AliDAQ.h"
#include "AliLog.h"
#include "AliRawReader.h"
#include "AliLoader.h"
#include "AliRunLoader.h"
#include "AliGeomManager.h"

#include "AliTOFcalib.h"
#include "AliTOFChannelOnlineArray.h"
#include "AliTOFChannelOnlineStatusArray.h"
#include "AliTOFChannelOffline.h"
#include "AliTOFClusterFinderV1.h"
#include "AliTOFcluster.h"
#include "AliTOFdigit.h"
#include "AliTOFDigitMap.h"
#include "AliTOFrawData.h"
#include "AliTOFReconstructor.h"
#include "AliTOFRecoParam.h"

ClassImp(AliTOFClusterFinderV1)

//_____________________________________________________________________________
AliTOFClusterFinderV1::AliTOFClusterFinderV1(AliTOFcalib *calib):
  TNamed("AliTOFClusterFinderV1",""),
  fRunLoader(0),
  fDigits(new TClonesArray("AliTOFdigit", 4000)),
  fRecPoints(new TClonesArray("AliTOFcluster", 4000)),
  fNumberOfTofClusters(0),
  fNumberOfTofDigits(0),
  fkRecoParam(0),//AliTOFReconstructor::GetRecoParam()),
  fMaxDeltaTime(0),//fkRecoParam->GetMaxDeltaTime()),
  fVerbose(0),
  fDecoderVersion(0),
  fTOFcalib(calib),
  fTOFdigitMap(new AliTOFDigitMap()),
  fTOFGeometry(new AliTOFGeometry()),
  fTOFdigits(new TTree()),
  fTOFRawStream(AliTOFRawStream()),
  fCalibrateTOFtimes(1)
{
//
// Constructor
//

  for (Int_t ii=0; ii<kTofMaxCluster; ii++) fTofClusters[ii]=0x0;

  if (AliTOFReconstructor::GetRecoParam()) {
    fkRecoParam = AliTOFReconstructor::GetRecoParam();
    fMaxDeltaTime = fkRecoParam->GetMaxDeltaTime();
  }
  else
    fMaxDeltaTime = 2;

  TString validity = (TString)fTOFcalib->GetOfflineValidity();
  if (validity.CompareTo("valid")==0) {
    AliInfo(Form(" validity = %s - Using offline calibration parameters", validity.Data()));
  } else {
    AliInfo(Form(" validity = %s - Using online calibration parameters", validity.Data()));
  }

}

//_____________________________________________________________________________
AliTOFClusterFinderV1::AliTOFClusterFinderV1(AliRunLoader* runLoader, AliTOFcalib *calib):
  TNamed("AliTOFClusterFinderV1",""),
  fRunLoader(runLoader),
  fDigits(new TClonesArray("AliTOFdigit", 4000)),
  fRecPoints(new TClonesArray("AliTOFcluster", 4000)),
  fNumberOfTofClusters(0),
  fNumberOfTofDigits(0),
  fkRecoParam(0),//AliTOFReconstructor::GetRecoParam()),
  fMaxDeltaTime(0),//fkRecoParam->GetMaxDeltaTime()),
  fVerbose(0),
  fDecoderVersion(0),
  fTOFcalib(calib),
  fTOFdigitMap(new AliTOFDigitMap()),
  fTOFGeometry(new AliTOFGeometry()),
  fTOFdigits(new TTree()),
  fTOFRawStream(AliTOFRawStream()),
  fCalibrateTOFtimes(1)
{
//
// Constructor
//

  for (Int_t ii=0; ii<kTofMaxCluster; ii++) fTofClusters[ii]=0x0;

  if (AliTOFReconstructor::GetRecoParam()) {
    fkRecoParam = AliTOFReconstructor::GetRecoParam();
    fMaxDeltaTime = fkRecoParam->GetMaxDeltaTime();
  }
  else
    fMaxDeltaTime = 2;

  TString validity = (TString)fTOFcalib->GetOfflineValidity();
  if (validity.CompareTo("valid")==0) {
    AliInfo(Form(" validity = %s - Using offline calibration parameters", validity.Data()));
  } else {
    AliInfo(Form(" validity = %s - Using online calibration parameters", validity.Data()));
  }

}
//_____________________________________________________________________________

AliTOFClusterFinderV1::AliTOFClusterFinderV1(const AliTOFClusterFinderV1 &source)
  :TNamed(source),
   fRunLoader(0),
   fDigits(source.fDigits),
   fRecPoints(source.fRecPoints),
   fNumberOfTofClusters(0),
   fNumberOfTofDigits(0),
   fkRecoParam(0),//AliTOFReconstructor::GetRecoParam()),
   fMaxDeltaTime(0),//fkRecoParam->GetMaxDeltaTime()),
   fVerbose(0),
   fDecoderVersion(source.fDecoderVersion),
   fTOFcalib(source.fTOFcalib),
   fTOFdigitMap(new AliTOFDigitMap()),
   fTOFGeometry(new AliTOFGeometry()),
   fTOFdigits(source.fTOFdigits),
   fTOFRawStream(source.fTOFRawStream),
   fCalibrateTOFtimes(1)
{
  // copy constructor

  for (Int_t ii=0; ii<kTofMaxCluster; ii++) fTofClusters[ii]=source.fTofClusters[ii];

  if (AliTOFReconstructor::GetRecoParam()) {
    fkRecoParam = AliTOFReconstructor::GetRecoParam();
    fMaxDeltaTime = fkRecoParam->GetMaxDeltaTime();
  }
  else
    fMaxDeltaTime = 2;

}
//_____________________________________________________________________________

AliTOFClusterFinderV1& AliTOFClusterFinderV1::operator=(const AliTOFClusterFinderV1 &source)
{
  // ass. op.

  if (this == &source)
    return *this;

  TObject::operator=(source);
  for (Int_t ii=0; ii<kTofMaxCluster; ii++) fTofClusters[ii]=source.fTofClusters[ii];
  fDigits=source.fDigits;
  fRecPoints=source.fRecPoints;
  fVerbose=source.fVerbose;
  fDecoderVersion=source.fDecoderVersion;
  fTOFcalib=source.fTOFcalib;
  fTOFdigitMap=source.fTOFdigitMap;
  fTOFGeometry=source.fTOFGeometry;
  fTOFdigits=source.fTOFdigits;
  fTOFRawStream=source.fTOFRawStream;
  fCalibrateTOFtimes=source.fCalibrateTOFtimes;
  return *this;

}
//_____________________________________________________________________________

AliTOFClusterFinderV1::~AliTOFClusterFinderV1()
{

  //
  // Destructor
  //

  if (fDigits)
    {
      fDigits->Delete();
      delete fDigits;
      fDigits=0;
    }
  if (fRecPoints)
    {
      fRecPoints->Delete();
      delete fRecPoints;
      fRecPoints=0;
    }

  delete fTOFdigitMap;

  delete fTOFGeometry;

  delete fTOFdigits;

  //if (fTofClusters || fNumberOfTofClusters) {
  if (fNumberOfTofClusters) {
    for (Int_t ii=0; ii<fNumberOfTofClusters; ii++)
      if (fTofClusters[ii]) fTofClusters[ii]->Delete();
    fNumberOfTofClusters = 0;
  }

}
//_____________________________________________________________________________

void AliTOFClusterFinderV1::Digits2RecPoints(TTree* digitsTree, TTree* clusterTree)
{
  //
  // Converts digits to recPoints for TOF
  //

  TStopwatch stopwatch;
  stopwatch.Start();

  Int_t inholes = 0;

  fDigits->Clear();
  TClonesArray &aDigits = *fDigits;

  if (digitsTree == 0x0) {
    AliFatal("Can not get TreeD for TOF");
    return;
  }

  TBranch *branch = digitsTree->GetBranch("TOF");
  if (!branch) {
    AliError("Can not get branch with the TOF digits !");
    return;
  }

  TClonesArray staticDigits("AliTOFdigit",10000);
  staticDigits.Clear();
  TClonesArray *digits = &staticDigits;
  branch->SetAddress(&digits);
  digitsTree->GetEvent(0);
  AliDebug(1,Form("Number of TOF digits: %d", digits->GetEntriesFast()));

  AliTOFdigit *tofDigit;

  Int_t jj = 0;
  Int_t detectorIndex[5];
  for (jj=0; jj<5; jj++) detectorIndex[jj] = -1;
  Int_t info[4];
  for (jj=0; jj<4; jj++) info[jj] = -1;
  Int_t *tracks;
  Int_t tdcCorr;
  Int_t dummy = -1;
  Int_t last = -1;

  Bool_t status = kTRUE;

  AliDebug(1," Calibrating TOF Digits");
  /*
  TString validity = (TString)fTOFcalib->GetOfflineValidity();
  if (validity.CompareTo("valid")==0) {
    AliInfo(Form(" validity = %s - Using offline calibration parameters", validity.Data()));
  } else
    AliInfo(Form(" validity = %s - Using online calibration parameters", validity.Data()));
  */

  Int_t ii = 0;
  for (ii=0; ii<digits->GetEntriesFast(); ii++) {
    tofDigit = (AliTOFdigit*)digits->UncheckedAt(ii);
    detectorIndex[0] = tofDigit->GetSector();
    detectorIndex[1] = tofDigit->GetPlate();
    detectorIndex[2] = tofDigit->GetStrip();
    detectorIndex[3] = tofDigit->GetPadz();
    detectorIndex[4] = tofDigit->GetPadx();

    if (detectorIndex[0]==13 || detectorIndex[0]==14 || detectorIndex[0]==15 ) { // sectors with holes
      if (detectorIndex[1]==2) { // plate with holes
	inholes++;
	continue;
      }
    }

    tdcCorr = tofDigit->GetTdc();
    status = MakeSlewingCorrection(detectorIndex, tofDigit->GetToT(), tofDigit->GetTdc(), tdcCorr);

    for (jj=0; jj<4; jj++) info[jj] = -1;
    info[0] = tdcCorr;//tofDigit->GetTdc();
    info[1] = tofDigit->GetAdc();
    info[2] = tofDigit->GetToT();
    info[3] = tofDigit->GetTdcND();//tofDigit->GetTdc();//
    tracks  = tofDigit->GetTracks();

    dummy = detectorIndex[3];
    detectorIndex[3] = detectorIndex[4];//padx
    detectorIndex[4] = dummy;//padz
    last = fDigits->GetEntriesFast();
    new (aDigits[last]) AliTOFdigit(tracks, detectorIndex, info);
    if (status) fTOFdigitMap->AddDigit(detectorIndex, last);

    AliDebug(2, Form(" Digits reading %2d -> %2d %1d %2d %1d %2d (%d, %d, %d)",
		     last,
		     detectorIndex[0], detectorIndex[1], detectorIndex[2], detectorIndex[3], detectorIndex[4],
		     info[0], info[1], info[3]));

  }

  fNumberOfTofDigits = fDigits->GetEntriesFast();

  ResetRecpoint();

  Int_t bufsize = 32000;
  clusterTree->Branch("TOF", &fRecPoints, bufsize);

  FillRecPoint();
  clusterTree->Fill();

  AliDebug(1,Form("Number of found clusters: %d", fNumberOfTofClusters));

  ResetRecpoint();

  fTOFdigitMap->Clear();

  ResetDigits();

  AliDebug(1,Form("Execution time to read TOF digits and to write TOF clusters : R:%.4fs C:%.4fs",
		  stopwatch.RealTime(),stopwatch.CpuTime()));

  if (inholes) AliWarning(Form("Clusters in the TOF holes: %d",inholes));

}
//_____________________________________________________________________________

void AliTOFClusterFinderV1::Digits2RecPoints(AliRawReader *rawReader, TTree *clustersTree)
{
  //
  // Converts raw data to recPoints for TOF
  //

  TStopwatch stopwatch;
  stopwatch.Start();


  AliDebug(2, "TreeD re-creation");
  //TTree *digitsTree = new TTree();
  //Raw2Digits(rawReader, digitsTree);

  Raw2Digits(rawReader, fTOFdigits);

  AliDebug(1,Form("Number of TOF digits: %d", fNumberOfTofDigits));
  ResetRecpoint();

  Int_t bufsize = 32000;
  clustersTree->Branch("TOF", &fRecPoints, bufsize);
  FillRecPoint();

  clustersTree->Fill();

  AliDebug(1,Form("Number of found clusters: %d", fNumberOfTofClusters));

  ResetRecpoint();

  fTOFdigitMap->Clear();

  ResetDigits();

  AliDebug(1,Form("Execution time to read TOF raw data and to write TOF clusters : R:%.4fs C:%.4fs",
		  stopwatch.RealTime(),stopwatch.CpuTime()));

}

//_____________________________________________________________________________

void AliTOFClusterFinderV1::Raw2Digits(AliRawReader *rawReader, TTree* digitsTree)
{
  //
  // Converts raw data to digits for TOF
  //

  TStopwatch stopwatch;
  stopwatch.Start();

  Int_t inholes = 0;

  const Int_t kMaxNumberOfTracksPerDigit = 3;
  const Int_t kDDL = AliDAQ::NumberOfDdls("TOF");

  digitsTree->Branch("TOF", &fDigits);
  TClonesArray &aDigits = *fDigits;

  fTOFRawStream.Clear();
  fTOFRawStream.SetRawReader(rawReader);

  ofstream ftxt;
  if (fVerbose==2) ftxt.open("TOFdigitsRead.txt",ios::app);

  TClonesArray staticRawData("AliTOFrawData",10000);
  staticRawData.Clear();
  TClonesArray * clonesRawData = &staticRawData;

  Int_t dummy = -1;
  Int_t detectorIndex[5] = {-1, -1, -1, -1, -1};
  Int_t digit[4];
  Int_t tracks[kMaxNumberOfTracksPerDigit];
  for (Int_t ii=0; ii<kMaxNumberOfTracksPerDigit; ii++)
    tracks[ii] = -1;
  Int_t last = -1;
  Int_t tdcCorr = 0;

  Bool_t status = kTRUE;

  /*
  TString validity = (TString)fTOFcalib->GetOfflineValidity();
  if (validity.CompareTo("valid")==0) {
    AliInfo(Form(" validity = %s - Using offline calibration parameters", validity.Data()));
  } else
    AliInfo(Form(" validity = %s - Using online calibration parameters", validity.Data()));
  */

  if (fDecoderVersion)
    AliInfo("Using New Decoder");

  Int_t indexDDL = 0;
  Int_t iRawData = 0;
  for (indexDDL=0; indexDDL<kDDL; indexDDL++) {

    rawReader->Reset();
    if (fDecoderVersion)
      fTOFRawStream.LoadRawDataBuffers(indexDDL,fVerbose);
    else fTOFRawStream.LoadRawData(indexDDL);

    clonesRawData = (TClonesArray*)fTOFRawStream.GetRawData();
    if (clonesRawData->GetEntriesFast()!=0) AliDebug(2,Form(" TOF raw data number = %3d", clonesRawData->GetEntriesFast()));
    for (iRawData = 0; iRawData<clonesRawData->GetEntriesFast(); iRawData++) {

      AliTOFrawData *tofRawDatum = (AliTOFrawData*)clonesRawData->UncheckedAt(iRawData);

      //if (tofRawDatum->GetTOT()==-1 || tofRawDatum->GetTOF()==-1) continue;
      if (tofRawDatum->GetTOF()==-1) continue;

      fTOFRawStream.EquipmentId2VolumeId(indexDDL, tofRawDatum->GetTRM(), tofRawDatum->GetTRMchain(),
					 tofRawDatum->GetTDC(), tofRawDatum->GetTDCchannel(), detectorIndex);

      tdcCorr = 0;
      dummy = detectorIndex[3];
      detectorIndex[3] = detectorIndex[4];//padz
      detectorIndex[4] = dummy;//padx

      tdcCorr = tofRawDatum->GetTOF();
      status = MakeSlewingCorrection(detectorIndex, tofRawDatum->GetTOT(), tofRawDatum->GetTOF(), tdcCorr);

      digit[0] = tdcCorr;
      digit[1] = tofRawDatum->GetTOT();
      digit[2] = tofRawDatum->GetTOT();
      digit[3] = -1;//tofRawDatum->GetTOF(); //tofND

      dummy = detectorIndex[3];
      detectorIndex[3] = detectorIndex[4];//padx
      detectorIndex[4] = dummy;//padz

      /* check valid index */
      if (detectorIndex[0]==-1||detectorIndex[1]==-1||detectorIndex[2]==-1||detectorIndex[3]==-1||detectorIndex[4]==-1) continue;

      // Do not reconstruct anything in the holes
      if (detectorIndex[0]==13 || detectorIndex[0]==14 || detectorIndex[0]==15 ) { // sectors with holes
	if (detectorIndex[1]==2) { // plate with holes
	  inholes++;
	  continue;
	}
      }

      last = fDigits->GetEntriesFast();
      new (aDigits[last]) AliTOFdigit(tracks, detectorIndex, digit);
      if (status) fTOFdigitMap->AddDigit(detectorIndex, last);

      if (fVerbose==2) {
	if (indexDDL<10) ftxt << "  " << indexDDL;
	else         ftxt << " " << indexDDL;
	if (tofRawDatum->GetTRM()<10) ftxt << "  " << tofRawDatum->GetTRM();
	else         ftxt << " " << tofRawDatum->GetTRM();
	ftxt << "  " << tofRawDatum->GetTRMchain();
	if (tofRawDatum->GetTDC()<10) ftxt << "  " << tofRawDatum->GetTDC();
	else         ftxt << " " << tofRawDatum->GetTDC();
	ftxt << "  " << tofRawDatum->GetTDCchannel();

	if (detectorIndex[0]<10) ftxt  << "  ->  " << detectorIndex[0];
	else              ftxt  << "  -> " << detectorIndex[0];
	ftxt << "  " << detectorIndex[1];
	if (detectorIndex[2]<10) ftxt << "  " << detectorIndex[2];
	else              ftxt << " " << detectorIndex[2];
	ftxt << "  " << detectorIndex[4];
	if (detectorIndex[4]<10) ftxt << "  " << detectorIndex[3];
	else              ftxt << " " << detectorIndex[3];

	if (digit[1]<10)ftxt << "        " << digit[1];
	else if (digit[1]>=10 && digit[1]<100) ftxt << "      " << digit[1];
	else ftxt << "      " << digit[1];
	if (digit[0]<10) ftxt << "      " << digit[0] << endl;
	else if (digit[0]>=10 && digit[0]<100)   ftxt << "    " << digit[0] << endl;
	else if (digit[0]>=100 && digit[0]<1000) ftxt << "    " << digit[0] << endl;
	else ftxt << "   " << digit[3] << endl;
      }

      AliDebug(2, Form(" Raw data reading %2d -> %2d %1d %2d %1d %2d (%d, %d, %d)",
		       last,
		       detectorIndex[0], detectorIndex[1], detectorIndex[2], detectorIndex[4], detectorIndex[3],
		       digit[0], digit[1], digit[3]));

    } // while loop

    clonesRawData->Clear();

  } // DDL Loop

  if (fVerbose==2) ftxt.close();

  digitsTree->Fill();

  fNumberOfTofDigits = fDigits->GetEntries();

  AliDebug(1, Form("Got %d TOF digits", fNumberOfTofDigits));
  AliDebug(1, Form("Execution time to read TOF raw data and fill TOF digit tree : R:%.2fs C:%.2fs",
		   stopwatch.RealTime(),stopwatch.CpuTime()));

  if (inholes) AliWarning(Form("Clusters in the TOF holes: %d",inholes));

}
//_____________________________________________________________________________

void AliTOFClusterFinderV1::FillRecPoint()
{
  //
  // Fill the global TClonesArray of AliTOFcluster objects,
  // i.e. fRecPoints
  //

  Int_t dummy4 = -1;
  Int_t dummy3 = -1;
  Int_t dummy2 = -1;
  Int_t dummy  = -1;

  for(Int_t iPlate=AliTOFGeometry::NPlates()-1; iPlate>=0; iPlate--) {
    for(Int_t iStrip=AliTOFGeometry::NStrip(iPlate)-1; iStrip>=0; iStrip--) {
      //for (Int_t iSector=AliTOFGeometry::NSectors()-1; iSector>=0; iSector--) {
      for (Int_t iSector=0; iSector<AliTOFGeometry::NSectors(); iSector++) {


	if (fTOFdigitMap->StripDigitCheck(iSector,iPlate,iStrip))
	  AliDebug(1,Form(" Number of TOF digits in (%2d,%1d,%2d) -> %d",
			  iSector,iPlate,iStrip,fTOFdigitMap->FilledCellsInStrip(iSector,iPlate,iStrip)));

	if (!(fTOFdigitMap->StripDigitCheck(iSector,iPlate,iStrip))) continue;
	FindClustersWithoutTOT(iSector, iPlate, iStrip); // clusters coming from digits without TOT measurement

	if (fMaxDeltaTime>0) {

	  if (!(fTOFdigitMap->StripDigitCheck(iSector,iPlate,iStrip))) continue;
	  //if (fTOFdigitMap->FilledCellsInStrip(iSector,iPlate,iStrip)>=4)
	  FindClustersPerStrip(iSector, iPlate, iStrip, 4); // 4 pads clusters
	  if (!(fTOFdigitMap->StripDigitCheck(iSector,iPlate,iStrip))) continue;

	  dummy4 = fNumberOfTofClusters;
	  FindClustersPerStrip(iSector, iPlate, iStrip, 4); // 4 pads clusters
	  if (fNumberOfTofClusters!=dummy4)
	    AliDebug(2, Form(" (4): n1= %5d, n2 = %5d", dummy4, fNumberOfTofClusters));


	  if (!(fTOFdigitMap->StripDigitCheck(iSector,iPlate,iStrip))) continue;
	  FindClustersPerStrip(iSector, iPlate, iStrip, 3); // 3 pads clusters
	  if (!(fTOFdigitMap->StripDigitCheck(iSector,iPlate,iStrip))) continue;

	  dummy3 = fNumberOfTofClusters;
	  FindClustersPerStrip(iSector, iPlate, iStrip, 3); // 3 pads clusters
	  if (fNumberOfTofClusters!=dummy3)
	    AliDebug(2, Form(" (3): n1= %5d, n2 = %5d", dummy3, fNumberOfTofClusters));


	  if (!(fTOFdigitMap->StripDigitCheck(iSector,iPlate,iStrip))) continue;
	  FindClustersPerStrip(iSector, iPlate, iStrip, 2); // 2 pads clusters
	  if (!(fTOFdigitMap->StripDigitCheck(iSector,iPlate,iStrip))) continue;

	  dummy2 = fNumberOfTofClusters;
	  FindClustersPerStrip(iSector, iPlate, iStrip, 2); // 2 pads clusters
	  if (fNumberOfTofClusters!=dummy2)
	    AliDebug(2, Form(" (2): n1= %5d, n2 =%5d", dummy2, fNumberOfTofClusters));


	  if (!(fTOFdigitMap->StripDigitCheck(iSector,iPlate,iStrip))) continue;
	  dummy = fNumberOfTofClusters;
	  FindClusters34(iSector, iPlate, iStrip); // 3 pads clusters between 4 hit pads
	  if (fNumberOfTofClusters!=dummy)
	    AliDebug(2, Form(" (3 between 4): n1 = %5d, n2 = %5d", fNumberOfTofClusters, dummy));


	  if (!(fTOFdigitMap->StripDigitCheck(iSector,iPlate,iStrip))) continue;
	  dummy = fNumberOfTofClusters;
	  FindClusters23(iSector, iPlate, iStrip); // 2 pads clusters between 3 hit pads
	  if (fNumberOfTofClusters!=dummy)
	    AliDebug(2, Form(" (2 between 3): n1 = %5d, n2 = %5d", fNumberOfTofClusters, dummy));

	  if (!(fTOFdigitMap->StripDigitCheck(iSector,iPlate,iStrip))) continue;
	  dummy = fNumberOfTofClusters;
	  FindClusters24(iSector, iPlate, iStrip); // 2 pads clusters between 4 hit pads
	  if (fNumberOfTofClusters!=dummy)
	    AliDebug(2, Form(" (2 between 4): n1 = %5d, n2 = %5d", fNumberOfTofClusters, dummy));


	  if (!(fTOFdigitMap->StripDigitCheck(iSector,iPlate,iStrip))) continue;
	  dummy = fNumberOfTofClusters;
	  FindOnePadClusterPerStrip(iSector, iPlate, iStrip); // 1 pad clusters
	  if (fNumberOfTofClusters!=dummy)
	    AliDebug(2,Form(" (1): n1 = %5d, n2 = %5d", fNumberOfTofClusters, dummy));

	  if (fTOFdigitMap->DigitInStrip(iSector,iPlate,iStrip)>0)
	    AliDebug(2, Form(" (1): number of clusters = %5d (remaining digit %2d), -%2d %1d %2d-",
			     fNumberOfTofClusters, fTOFdigitMap->DigitInStrip(iSector,iPlate,iStrip),
			     iSector, iPlate, iStrip));

	}
	else {
	  if (!(fTOFdigitMap->StripDigitCheck(iSector,iPlate,iStrip))) continue;
	  dummy = fNumberOfTofClusters;
	  FindOnePadClusterPerStrip(iSector, iPlate, iStrip); // 1 pad clusters
	  if (fNumberOfTofClusters!=dummy)
	    AliDebug(2,Form(" (1): n1 = %5d, n2 = %5d", fNumberOfTofClusters, dummy));

	  if (fTOFdigitMap->DigitInStrip(iSector,iPlate,iStrip)>0)
	    AliDebug(2, Form(" (1): number of clusters = %5d (remaining digit %2d), -%2d %1d %2d-",
			     fNumberOfTofClusters, fTOFdigitMap->DigitInStrip(iSector,iPlate,iStrip),
			     iSector, iPlate, iStrip));

	}


      }
    }
  }


  TClonesArray &lRecPoints = *fRecPoints;
  
  Int_t ii, jj;

  Int_t detectorIndex[5];
  for (jj=0; jj<5; jj++) detectorIndex[jj] = -1;
  Int_t parTOF[7];
  for (jj=0; jj<7; jj++) parTOF[jj] = -1;
  Int_t trackLabels[3];
  for (jj=0; jj<3; jj++) trackLabels[jj] = -1;
  Int_t digitIndex = -1;
  Bool_t status = kTRUE;
  Float_t posClus[3];
  for (ii=0; ii<3; ii++) posClus[ii] = 0.;
  //Float_t covClus[6];
  //for (ii=0; ii<6; ii++) covClus[ii] = 0.;
  UShort_t volIdClus;

  for (ii=0; ii<fNumberOfTofClusters; ii++) {

    digitIndex = fTofClusters[ii]->GetIndex();
    for(jj=0; jj<5; jj++) detectorIndex[jj] = fTofClusters[ii]->GetDetInd(jj);
    volIdClus = fTOFGeometry->GetAliSensVolIndex(detectorIndex[0],detectorIndex[1],detectorIndex[2]);
    //volIdClus = GetClusterVolIndex(detectorIndex);
    for(jj=0; jj<3; jj++) trackLabels[jj] = fTofClusters[ii]->GetLabel(jj);
    parTOF[0] = fTofClusters[ii]->GetTDC(); // TDC
    parTOF[1] = fTofClusters[ii]->GetToT(); // TOT
    parTOF[2] = fTofClusters[ii]->GetADC(); // ADC=TOT
    parTOF[3] = fTofClusters[ii]->GetTDCND(); // TDCND
    parTOF[4] = fTofClusters[ii]->GetTDCRAW();//RAW
    parTOF[5] = 0;
    parTOF[6] = 0;
    status = fTofClusters[ii]->GetStatus();

    posClus[0] = fTofClusters[ii]->GetX();
    posClus[1] = fTofClusters[ii]->GetY();
    posClus[2] = fTofClusters[ii]->GetZ();

    //for (jj=0; jj<6; jj++) covClus[jj] = 0.;
    //((AliCluster*)fTofClusters[ii])->GetGlobalCov(covClus);

    new(lRecPoints[ii]) AliTOFcluster(volIdClus, (Double_t)posClus[0], (Double_t)posClus[1], (Double_t)posClus[2],
				      (Double_t)(fTofClusters[ii]->GetSigmaX2()),
				      (Double_t)(fTofClusters[ii]->GetSigmaXY()),
				      (Double_t)(fTofClusters[ii]->GetSigmaXZ()),
				      (Double_t)(fTofClusters[ii]->GetSigmaY2()),
				      (Double_t)(fTofClusters[ii]->GetSigmaYZ()),
				      (Double_t)(fTofClusters[ii]->GetSigmaZ2()),
				      //(Double_t)covClus[0], (Double_t)covClus[1], (Double_t)covClus[2],
				      //(Double_t)covClus[3], (Double_t)covClus[4], (Double_t)covClus[5],
				      trackLabels, detectorIndex, parTOF, status, digitIndex);

    AliDebug(2, Form(" %4d  %4d  %f %f %f  %f %f %f %f %f %f  %3d %3d %3d  %2d %1d %2d %1d %2d  %4d %3d %3d %4d %4d  %1d  %4d", 
		     ii, volIdClus, posClus[0], posClus[1], posClus[2],
		     fTofClusters[ii]->GetSigmaX2(),
		     fTofClusters[ii]->GetSigmaXY(),
		     fTofClusters[ii]->GetSigmaXZ(),
		     fTofClusters[ii]->GetSigmaY2(),
		     fTofClusters[ii]->GetSigmaYZ(),
		     fTofClusters[ii]->GetSigmaZ2(),
		     trackLabels[0], trackLabels[1], trackLabels[2],
		     detectorIndex[0], detectorIndex[1], detectorIndex[2], detectorIndex[3], detectorIndex[4],
		     parTOF[0], parTOF[1], parTOF[2], parTOF[3], parTOF[4],
		     status, digitIndex));


  }

  for (Int_t iSector=0; iSector<AliTOFGeometry::NSectors(); iSector++)
    for(Int_t iPlate=0; iPlate<AliTOFGeometry::NPlates(); iPlate++) {
      for(Int_t iStrip=0; iStrip<AliTOFGeometry::NStrip(iPlate); iStrip++) {
	if (!(fTOFdigitMap->StripDigitCheck(iSector,iPlate,iStrip))) continue;
	AliDebug(2, Form(" END %2d %1d %2d   %5d",
			 iSector, iPlate, iStrip, fTOFdigitMap->DigitInStrip(iSector,iPlate,iStrip)));
      }
    }

}
//_____________________________________________________________________________

void AliTOFClusterFinderV1::FindOnePadClusterPerStrip(Int_t nSector,
						      Int_t nPlate,
						      Int_t nStrip)
{
  //
  // This function searches the isolated digits (stored in the fDigits object),
  // to perform clusters (stored in the fTofClusters array).
  // This research has been made by checking the fTOFdigitMap object,
  // filled at digits/raw-data reading time.
  //

  const Int_t kMaxNumberOfTracksPerDigit = 3;
  const Int_t kMaxNumberOfDigitsPerVolume = 10;

  Int_t jj = 0;

  Int_t det[5] = {nSector,nPlate,nStrip,-1,-1};//sector,plate,strip,padZ,padX
  Int_t vol[5] = {nSector,nPlate,nStrip,-1,-1};//sector,plate,strip,padX,padZ
  UShort_t volIdClus = 0;

  Float_t pos[3];
  for (jj=0; jj<3; jj++) pos[jj] = 0.;
  Double_t posClus[3];
  for (jj=0; jj<3; jj++) posClus[jj] = 0.;

  Double_t covClus[6];
  for (jj=0; jj<6; jj++) covClus[jj] = 0.;

  Int_t parTOF[7];
  for (jj=0; jj<7; jj++) parTOF[jj] = 0;

  Bool_t status = kTRUE; //assume all sim channels ok in the beginning...

  Int_t tracks[kMaxNumberOfTracksPerDigit];
  for (jj=0; jj<kMaxNumberOfTracksPerDigit; jj++) tracks[jj] = -1;

  Int_t dummyCounter=-1;

  AliTOFdigit *digitInteresting;

  Int_t iPadX = -1;
  Int_t iPadZ = -1;
  for (iPadX=0; iPadX<AliTOFGeometry::NpadX(); iPadX++) {
    for (iPadZ=0; iPadZ<AliTOFGeometry::NpadZ(); iPadZ++) {
      vol[4] = iPadZ  , vol[3]  = iPadX;

      AliDebug(3, Form(" %1d %2d\n", iPadZ, iPadX));

      if (fTOFdigitMap->GetNumberOfDigits(vol)==0) continue;

      for(Int_t digIndex=0; digIndex<kMaxNumberOfDigitsPerVolume; digIndex++) {
	if (fTOFdigitMap->GetDigitIndex(vol,digIndex)<0) continue;
	digitInteresting = (AliTOFdigit*)fDigits->UncheckedAt(fTOFdigitMap->GetDigitIndex(vol,digIndex));

	AliDebug(2, Form(" %3d  %5d    %2d %1d %2d %1d %2d  %d %d %d  %5d  %5d %5d %5d",
			 fTOFdigitMap->GetNumberOfDigits(vol), digIndex,
			 vol[0], vol[1], vol[2] ,vol[4], vol[3],
			 digitInteresting->GetTdc(), digitInteresting->GetAdc(),
			 digitInteresting->GetToT(),
			 fTOFdigitMap->GetDigitIndex(vol,digIndex),
			 digitInteresting->GetTrack(0), digitInteresting->GetTrack(1), digitInteresting->GetTrack(2)));
	
	det[3] = vol[4]; // padz
	det[4] = vol[3]; // padx
	fTOFGeometry->GetPosPar(det,pos);
	AliDebug(1,Form(" %f %f %f", pos[0], pos[1], pos[2]));

	//insert cluster
	for (jj=0; jj<3; jj++) posClus[jj] = pos[jj];

	parTOF[0] = Int_t(digitInteresting->GetTdc());
	parTOF[1] = Int_t(digitInteresting->GetToT());
	parTOF[2] = Int_t(digitInteresting->GetAdc());
	parTOF[3] = Int_t(digitInteresting->GetTdcND());
	parTOF[4] = Int_t(digitInteresting->GetTdc());
	parTOF[5] = 0;
	parTOF[6] = 0;

	volIdClus = fTOFGeometry->GetAliSensVolIndex(det[0],det[1],det[2]);
	//volIdClus = GetClusterVolIndex(det);

	for (jj=0; jj<6; jj++) covClus[jj] = 0.;
	GetClusterPars(det, posClus, covClus);

	// To fill the track index array
	dummyCounter=-1;
	for (jj=0; jj<kMaxNumberOfTracksPerDigit; jj++) tracks[jj] = -1;
	for (jj=0; jj<kMaxNumberOfTracksPerDigit; jj++) { // three is the max number of tracks associated to one digit
	  if (digitInteresting->GetTrack(jj)==-1) continue;
	  else {
	    dummyCounter++;
	    tracks[dummyCounter] = digitInteresting->GetTrack(jj);
	  }
	}

	AliTOFcluster *tofCluster =
	  new AliTOFcluster(volIdClus, posClus[0], posClus[1], posClus[2],
			    covClus[0], covClus[1], covClus[2], covClus[3], covClus[4], covClus[5],
			    tracks, det, parTOF, status, fTOFdigitMap->GetDigitIndex(vol,digIndex));
	InsertCluster(tofCluster);

	AliDebug(2, Form("       %4d  %f %f %f  %f %f %f %f %f %f  %3d %3d %3d  %2d %1d %2d %1d %2d  %4d %3d %3d %4d %4d  %1d  %4d", 
			 volIdClus, posClus[0], posClus[1], posClus[2],
			 covClus[0], covClus[1], covClus[2], covClus[3], covClus[4], covClus[5],
			 tracks[0], tracks[1], tracks[2],
			 det[0], det[1], det[2], det[3], det[4],
			 parTOF[0], parTOF[1], parTOF[2], parTOF[3], parTOF[4],
			 status, fTOFdigitMap->GetDigitIndex(vol,digIndex)));

	AliDebug(2, Form("        %f %f %f", pos[0], pos[1], pos[2]));
	AliDebug(2, Form("           %d %d", parTOF[0], parTOF[2]));

	fTOFdigitMap->ResetDigit(vol, digIndex);

      }

    }
  }

}
//_____________________________________________________________________________

void AliTOFClusterFinderV1::FindClustersWithoutTOT(Int_t nSector,
						   Int_t nPlate,
						   Int_t nStrip)
{
  //
  // This function searches the isolated digits without TOT
  // measurement (stored in the fDigits object), to perform clusters
  // (stored in the fTofClusters array). This research has been made
  // by checking the fTOFdigitMap object, filled at digits/raw-data
  // reading time.
  //

  const Int_t kMaxNumberOfTracksPerDigit = 3;
  const Int_t kMaxNumberOfDigitsPerVolume = 10;

  Int_t jj = 0;

  Int_t det[5] = {nSector,nPlate,nStrip,-1,-1};//sector,plate,strip,padZ,padX
  Int_t vol[5] = {nSector,nPlate,nStrip,-1,-1};//sector,plate,strip,padX,padZ
  UShort_t volIdClus = 0;

  Float_t pos[3];
  for (jj=0; jj<3; jj++) pos[jj] = 0.;
  Double_t posClus[3];
  for (jj=0; jj<3; jj++) posClus[jj] = 0.;

  Double_t covClus[6];
  for (jj=0; jj<6; jj++) covClus[jj] = 0.;

  Int_t parTOF[7];
  for (jj=0; jj<7; jj++) parTOF[jj] = 0;

  Bool_t status = kTRUE; //assume all sim channels ok in the beginning...
  Int_t tracks[kMaxNumberOfTracksPerDigit];
  for (jj=0; jj<kMaxNumberOfTracksPerDigit; jj++) tracks[jj] = -1;

  Int_t dummyCounter=-1;

  AliTOFdigit *digitInteresting;

  Int_t iPadX = -1;
  Int_t iPadZ = -1;
  for (iPadX=0; iPadX<AliTOFGeometry::NpadX(); iPadX++) {
    for (iPadZ=0; iPadZ<AliTOFGeometry::NpadZ(); iPadZ++) {
      vol[4] = iPadZ  , vol[3]  = iPadX;

      AliDebug(3, Form(" %1d %2d\n", iPadZ, iPadX));

      if (fTOFdigitMap->GetNumberOfDigits(vol)==0) continue;

      for(Int_t digIndex=0; digIndex<kMaxNumberOfDigitsPerVolume; digIndex++) {
	if (fTOFdigitMap->GetDigitIndex(vol,digIndex)<0) continue;
	digitInteresting = (AliTOFdigit*)fDigits->UncheckedAt(fTOFdigitMap->GetDigitIndex(vol,digIndex));
	if (digitInteresting->GetToT()>0) continue; // AdC

	AliDebug(2, Form(" %3d  %5d    %2d %1d %2d %1d %2d  %d %d %d  %5d  %5d %5d %5d",
			 fTOFdigitMap->GetNumberOfDigits(vol), digIndex,
			 vol[0], vol[1], vol[2] ,vol[4], vol[3],
			 digitInteresting->GetTdc(), digitInteresting->GetAdc(),
			 digitInteresting->GetToT(),
			 fTOFdigitMap->GetDigitIndex(vol,digIndex),
			 digitInteresting->GetTrack(0), digitInteresting->GetTrack(1), digitInteresting->GetTrack(2)));
	
	det[3] = vol[4]; // padz
	det[4] = vol[3]; // padx
	fTOFGeometry->GetPosPar(det,pos);
	AliDebug(1,Form(" %f %f %f", pos[0], pos[1], pos[2]));

	//insert cluster
	for (jj=0; jj<3; jj++) posClus[jj] = pos[jj];

	parTOF[0] = Int_t(digitInteresting->GetTdc());
	parTOF[1] = Int_t(digitInteresting->GetToT());
	parTOF[2] = Int_t(digitInteresting->GetAdc());
	parTOF[3] = Int_t(digitInteresting->GetTdcND());
	parTOF[4] = Int_t(digitInteresting->GetTdc());
	parTOF[5] = 0;
	parTOF[6] = 0;

	volIdClus = fTOFGeometry->GetAliSensVolIndex(det[0],det[1],det[2]);
	//volIdClus = GetClusterVolIndex(det);

	for (jj=0; jj<6; jj++) covClus[jj] = 0.;
	GetClusterPars(det, posClus, covClus);

	// To fill the track index array
	dummyCounter=-1;
	for (jj=0; jj<kMaxNumberOfTracksPerDigit; jj++) tracks[jj] = -1;
	for (jj=0; jj<kMaxNumberOfTracksPerDigit; jj++) { // three is the max number of tracks associated to one digit
	  if (digitInteresting->GetTrack(jj)==-1) continue;
	  else {
	    dummyCounter++;
	    tracks[dummyCounter] = digitInteresting->GetTrack(jj);
	  }
	}

	AliTOFcluster *tofCluster =
	  new AliTOFcluster(volIdClus, posClus[0], posClus[1], posClus[2],
			    covClus[0], covClus[1], covClus[2], covClus[3], covClus[4], covClus[5],
			    tracks, det, parTOF, status, fTOFdigitMap->GetDigitIndex(vol,digIndex));
	InsertCluster(tofCluster);

	AliDebug(2, Form("       %4d  %f %f %f  %f %f %f %f %f %f  %3d %3d %3d  %2d %1d %2d %1d %2d  %4d %3d %3d %4d %4d  %1d  %4d", 
			 volIdClus, posClus[0], posClus[1], posClus[2],
			 covClus[0], covClus[1], covClus[2], covClus[3], covClus[4], covClus[5],
			 tracks[0], tracks[1], tracks[2],
			 det[0], det[1], det[2], det[3], det[4],
			 parTOF[0], parTOF[1], parTOF[2], parTOF[3], parTOF[4],
			 status, fTOFdigitMap->GetDigitIndex(vol,digIndex)));

	AliDebug(2, Form("        %f %f %f", pos[0], pos[1], pos[2]));
	AliDebug(2, Form("           %d %d", parTOF[0], parTOF[2]));

	fTOFdigitMap->ResetDigit(vol, digIndex);

      }

    }
  }

}
//_____________________________________________________________________________

void AliTOFClusterFinderV1::FindClusters34(Int_t nSector,
					   Int_t nPlate,
					   Int_t nStrip)
{
  //
  // This function searches the neighbouring digits (stored in the fDigits object),
  // to perform clusters (stored in the fTofClusters array).
  //
  // This research has been made by checking the fTOFdigitMap object,
  // filled at digits/raw-data reading time.
  //

  const Int_t kMaxNumberOfInterestingPads = 4;
  const Int_t kMaxNumberOfTracksPerDigit = 3;
  const Int_t kMaxNumberOfDigitsPerVolume = 10;

  Int_t ii = 0;

  Int_t digitsInVolumeIndices[kMaxNumberOfDigitsPerVolume];
  for (ii=0; ii<kMaxNumberOfDigitsPerVolume; ii++)
    digitsInVolumeIndices[ii] = -1;

  Int_t vol[5] = {nSector,nPlate,nStrip,-1,-1};

  Float_t pos[3] = {0.,0.,0.};

  Int_t jj = 0;
  Int_t interestingPadX[kMaxNumberOfInterestingPads];
  for (jj=0; jj<kMaxNumberOfInterestingPads; jj++) interestingPadX[jj] = -1;
  Int_t interestingPadZ[kMaxNumberOfInterestingPads];
  for (jj=0; jj<kMaxNumberOfInterestingPads; jj++) interestingPadZ[jj] = -1;
  Double_t interestingTOT[kMaxNumberOfInterestingPads];
  for (jj=0; jj<kMaxNumberOfInterestingPads; jj++) interestingTOT[jj] = 0;
  Double_t interestingADC[kMaxNumberOfInterestingPads];
  for (jj=0; jj<kMaxNumberOfInterestingPads; jj++) interestingADC[jj] = 0;
  Double_t interestingTOF[kMaxNumberOfInterestingPads];
  for (jj=0; jj<kMaxNumberOfInterestingPads; jj++) interestingTOF[jj] = 0;
  Double_t interestingWeight[kMaxNumberOfInterestingPads];
  for (jj=0; jj<kMaxNumberOfInterestingPads; jj++) interestingWeight[jj] = 0;

  Float_t interestingX[kMaxNumberOfInterestingPads];
  for (jj=0; jj<kMaxNumberOfInterestingPads; jj++) interestingX[jj] = 0;
  Float_t interestingY[kMaxNumberOfInterestingPads];
  for (jj=0; jj<kMaxNumberOfInterestingPads; jj++) interestingY[jj] = 0;
  Float_t interestingZ[kMaxNumberOfInterestingPads];
  for (jj=0; jj<kMaxNumberOfInterestingPads; jj++) interestingZ[jj] = 0;

  Float_t interDigit[kMaxNumberOfInterestingPads];
  for (jj=0; jj<kMaxNumberOfInterestingPads; jj++) interDigit[jj] = 0;

  Int_t padsCluster[11];
  padsCluster[0] = nSector;
  padsCluster[1] = nPlate;
  padsCluster[2] = nStrip;
  for (jj=3; jj<11; jj++) padsCluster[jj] = -1;

  Int_t interestingCounter=-1;
  Int_t  digitIndexLocal=-1; // AdC
  Int_t iPad  = -1;
  Int_t iPadX = -1;
  Int_t iPadZ = -1;

  Int_t parTOF[7];
  for (jj=0; jj<7; jj++) parTOF[jj] = 0;
  Double_t posClus[3];
  for (jj=0; jj<3; jj++) posClus[jj] = 0.;
  Int_t det[5];
  for (jj=0; jj<5; jj++) det[jj] = -1;
  Float_t posF[3];
  for (jj=0; jj<3; jj++) posF[jj] = 0.;
  UShort_t volIdClus = 0;
  Bool_t check = kFALSE;
  Bool_t status = kTRUE; //assume all sim channels ok in the beginning...
  Double_t covClus[6];
  for (jj=0; jj<6; jj++) covClus[jj] = 0.;
  Int_t tracks[kMaxNumberOfTracksPerDigit];
  for (jj=0; jj<kMaxNumberOfTracksPerDigit; jj++) tracks[jj] = -1;
  Int_t dummyCounter=-1;
  Bool_t alreadyStored = kFALSE;

  AliTOFselectedDigit ***selectedDigit = new AliTOFselectedDigit**[kMaxNumberOfInterestingPads];
  for (ii=0; ii<kMaxNumberOfInterestingPads; ii++)
    selectedDigit[ii] = new AliTOFselectedDigit*[kMaxNumberOfDigitsPerVolume];

  for (ii=0; ii<kMaxNumberOfInterestingPads; ii++)
    for (jj=0; jj<kMaxNumberOfDigitsPerVolume; jj++) selectedDigit[ii][jj] = 0x0;

  AliTOFdigit *digitInteresting;

  for (iPad=0; iPad<AliTOFGeometry::NpadZ()*AliTOFGeometry::NpadX()-3; iPad+=2) {

    iPadZ = iPad%AliTOFGeometry::NpadZ(); //iPad%2;
    iPadX = iPad/AliTOFGeometry::NpadZ(); //iPad/2;

    AliDebug(3, Form("%2d %1d %2d\n", iPad, iPadZ, iPadX));






    interestingCounter=-1;

    vol[4] = iPadZ  , vol[3]  = iPadX;
    if (fTOFdigitMap->GetNumberOfDigits(vol)>0)
      interestingCounter++;

    vol[4] = iPadZ, vol[3] = iPadX+1;
    if (fTOFdigitMap->GetNumberOfDigits(vol)>0)
      interestingCounter++;

    vol[4] = AliTOFGeometry::NpadZ()-(iPadZ+1), vol[3] = iPadX;
    if (fTOFdigitMap->GetNumberOfDigits(vol)>0)
      interestingCounter++;

    vol[4] = AliTOFGeometry::NpadZ()-(iPadZ+1), vol[3] = iPadX+1;
    if (fTOFdigitMap->GetNumberOfDigits(vol)>0)
      interestingCounter++;

    if (interestingCounter+1!=4) continue; // the hit pads have to be 4
    else interestingCounter=-1;


    for (ii=0; ii<kMaxNumberOfInterestingPads; ii++)
      for (jj=0; jj<kMaxNumberOfDigitsPerVolume; jj++)
	selectedDigit[ii][jj] = 0x0;


    vol[4] = iPadZ, vol[3] = iPadX;

    if (fTOFdigitMap->GetNumberOfDigits(vol)>0) {
      interestingCounter++;
      for (ii=0; ii<kMaxNumberOfDigitsPerVolume; ii++)
	digitsInVolumeIndices[ii] = -1;
      fTOFdigitMap->GetDigitIndex(vol, digitsInVolumeIndices);
      digitIndexLocal=-1; // AdC
      for(Int_t digIndex=0; digIndex<kMaxNumberOfDigitsPerVolume; digIndex++) {
	if (digitsInVolumeIndices[digIndex]<0) continue;
	digitInteresting = (AliTOFdigit*)fDigits->UncheckedAt(digitsInVolumeIndices[digIndex]);
	if (digitInteresting->GetToT()<=0) continue; // AdC
	digitIndexLocal++; // AdC

	AliDebug(1,Form(" %2d %1d %2d %1d %2d  %d %d %d %d  %5d  %5d %5d %5d",
			vol[0], vol[1], vol[2] ,vol[4], vol[3],
			digitInteresting->GetTdc(), digitInteresting->GetAdc(),
			digitInteresting->GetToT(), digitInteresting->GetToT()*digitInteresting->GetToT(),
			digitsInVolumeIndices[digIndex],
			digitInteresting->GetTrack(0), digitInteresting->GetTrack(1), digitInteresting->GetTrack(2)));


	selectedDigit[interestingCounter][digitIndexLocal] = new // AdC
	  AliTOFselectedDigit(vol, digitInteresting->GetTdc(),
			      digitInteresting->GetAdc(), digitInteresting->GetToT(),
			      digitInteresting->GetToT()*digitInteresting->GetToT(),
			      digitsInVolumeIndices[digIndex],
			      digitInteresting->GetTracks());
      }
      if (digitIndexLocal==-1) interestingCounter--; // AdC
    }


    vol[4] = iPadZ, vol[3] = iPadX+1;

    if (fTOFdigitMap->GetNumberOfDigits(vol)>0) {
      interestingCounter++;
      for (ii=0; ii<kMaxNumberOfDigitsPerVolume; ii++)
	digitsInVolumeIndices[ii] = -1;
      fTOFdigitMap->GetDigitIndex(vol, digitsInVolumeIndices);
      digitIndexLocal=-1; // AdC
      for(Int_t digIndex=0; digIndex<kMaxNumberOfDigitsPerVolume; digIndex++) {
	if (digitsInVolumeIndices[digIndex]<0) continue;
	digitInteresting = (AliTOFdigit*)fDigits->UncheckedAt(digitsInVolumeIndices[digIndex]);
	if (digitInteresting->GetToT()<=0) continue; // AdC
	digitIndexLocal++; // AdC

	AliDebug(1,Form(" %2d %1d %2d %1d %2d  %d %d %d %d  %5d  %5d %5d %5d",
			vol[0], vol[1], vol[2] ,vol[4], vol[3],
			digitInteresting->GetTdc(), digitInteresting->GetAdc(),
			digitInteresting->GetToT(), digitInteresting->GetToT()*digitInteresting->GetToT(),
			digitsInVolumeIndices[digIndex],
			digitInteresting->GetTrack(0), digitInteresting->GetTrack(1), digitInteresting->GetTrack(2)));


	selectedDigit[interestingCounter][digitIndexLocal] = new // AdC
	  AliTOFselectedDigit(vol, digitInteresting->GetTdc(),
			      digitInteresting->GetAdc(), digitInteresting->GetToT(),
			      digitInteresting->GetToT()*digitInteresting->GetToT(),
			      digitsInVolumeIndices[digIndex],
			      digitInteresting->GetTracks());
      }
      if (digitIndexLocal==-1) interestingCounter--; // AdC
    }


    vol[4] = AliTOFGeometry::NpadZ()-(iPadZ+1), vol[3] = iPadX;

    if (fTOFdigitMap->GetNumberOfDigits(vol)>0) {
      interestingCounter++;
      for (ii=0; ii<kMaxNumberOfDigitsPerVolume; ii++)
	digitsInVolumeIndices[ii] = -1;
      fTOFdigitMap->GetDigitIndex(vol, digitsInVolumeIndices);
      digitIndexLocal=-1; // AdC
      for(Int_t digIndex=0; digIndex<kMaxNumberOfDigitsPerVolume; digIndex++) {
	if (digitsInVolumeIndices[digIndex]<0) continue;
	digitInteresting = (AliTOFdigit*)fDigits->UncheckedAt(digitsInVolumeIndices[digIndex]);
	if (digitInteresting->GetToT()<=0) continue; // AdC
	digitIndexLocal++; // AdC

	AliDebug(1,Form(" %2d %1d %2d %1d %2d  %d %d %d %d  %5d  %5d %5d %5d",
			vol[0], vol[1], vol[2] ,vol[4], vol[3],
			digitInteresting->GetTdc(), digitInteresting->GetAdc(),
			digitInteresting->GetToT(), digitInteresting->GetToT()*digitInteresting->GetToT(),
			digitsInVolumeIndices[digIndex],
			digitInteresting->GetTrack(0), digitInteresting->GetTrack(1), digitInteresting->GetTrack(2)));


	selectedDigit[interestingCounter][digitIndexLocal] = new // AdC
	  AliTOFselectedDigit(vol, digitInteresting->GetTdc(),
			      digitInteresting->GetAdc(), digitInteresting->GetToT(),
			      digitInteresting->GetToT()*digitInteresting->GetToT(),
			      digitsInVolumeIndices[digIndex],
			      digitInteresting->GetTracks());
      }
      if (digitIndexLocal==-1) interestingCounter--; // AdC
    }


    vol[4] = AliTOFGeometry::NpadZ()-(iPadZ+1), vol[3] = iPadX+1;

    if (fTOFdigitMap->GetNumberOfDigits(vol)>0) {
      interestingCounter++;
      for (ii=0; ii<kMaxNumberOfDigitsPerVolume; ii++)
	digitsInVolumeIndices[ii] = -1;
      fTOFdigitMap->GetDigitIndex(vol, digitsInVolumeIndices);
      digitIndexLocal=-1;
      for(Int_t digIndex=0; digIndex<kMaxNumberOfDigitsPerVolume; digIndex++) {
	if (digitsInVolumeIndices[digIndex]<0) continue;
	digitInteresting = (AliTOFdigit*)fDigits->UncheckedAt(digitsInVolumeIndices[digIndex]);
	if (digitInteresting->GetToT()<=0) continue; // AdC
	digitIndexLocal++; // AdC

	AliDebug(1,Form(" %2d %1d %2d %1d %2d  %d %d %d %d  %5d  %5d %5d %5d",
			vol[0], vol[1], vol[2] ,vol[4], vol[3],
			digitInteresting->GetTdc(), digitInteresting->GetAdc(),
			digitInteresting->GetToT(), digitInteresting->GetToT()*digitInteresting->GetToT(),
			digitsInVolumeIndices[digIndex],
			digitInteresting->GetTrack(0), digitInteresting->GetTrack(1), digitInteresting->GetTrack(2)));

	selectedDigit[interestingCounter][digitIndexLocal] = new // AdC
	  AliTOFselectedDigit(vol, digitInteresting->GetTdc(),
			      digitInteresting->GetAdc(), digitInteresting->GetToT(),
			      digitInteresting->GetToT()*digitInteresting->GetToT(),
			      digitsInVolumeIndices[digIndex],
			      digitInteresting->GetTracks());
      }
      if (digitIndexLocal==-1) interestingCounter--; // AdC
    }

    AliDebug(1,Form(" e adesso %1d", interestingCounter+1));

    for (Int_t adesso1=0; adesso1<interestingCounter+1; adesso1++) {
      for (Int_t firstIndex=0; firstIndex<kMaxNumberOfDigitsPerVolume; firstIndex++) {
	if (selectedDigit[adesso1][firstIndex]==0x0) continue;

	for (Int_t adesso2=adesso1+1; adesso2<interestingCounter+1; adesso2++) {
	  for (Int_t secondIndex=0; secondIndex<kMaxNumberOfDigitsPerVolume; secondIndex++) {
	    if (selectedDigit[adesso2][secondIndex]==0x0) continue;

	    for (Int_t adesso3=adesso2+1; adesso3<interestingCounter+1; adesso3++) {
	      for (Int_t thirdIndex=0; thirdIndex<kMaxNumberOfDigitsPerVolume; thirdIndex++) {
		if (selectedDigit[adesso3][thirdIndex]==0x0) continue;


		if (TMath::Abs(selectedDigit[adesso1][firstIndex]->GetTDC()-selectedDigit[adesso2][secondIndex]->GetTDC())>fMaxDeltaTime
		    ||
		    TMath::Abs(selectedDigit[adesso1][firstIndex]->GetTDC()-selectedDigit[adesso3][thirdIndex]->GetTDC())>fMaxDeltaTime
		    ||
		    TMath::Abs(selectedDigit[adesso2][secondIndex]->GetTDC()-selectedDigit[adesso3][thirdIndex]->GetTDC())>fMaxDeltaTime) continue;

		interestingTOF[0] = selectedDigit[adesso1][firstIndex]->GetTDC();
		interestingTOT[0] = selectedDigit[adesso1][firstIndex]->GetTOT();
		interestingADC[0] = selectedDigit[adesso1][firstIndex]->GetADC();
		interestingWeight[0] = selectedDigit[adesso1][firstIndex]->GetWeight();
		Int_t vol1[5]; for(jj=0; jj<5; jj++) vol1[jj] = selectedDigit[adesso1][firstIndex]->GetDetectorIndex(jj);
		AliDebug(1,Form(" %2d %1d %2d %1d %2d", vol1[0], vol1[1], vol1[2], vol1[4], vol1[3]));
		Int_t volDum = vol1[3];
		vol1[3] = vol1[4];
		vol1[4] = volDum;
		fTOFGeometry->GetPosPar(vol1,pos);
		interestingX[0] = pos[0];
		interestingY[0] = pos[1];
		interestingZ[0] = pos[2];

		interestingTOF[1] = selectedDigit[adesso2][secondIndex]->GetTDC();
		interestingTOT[1] = selectedDigit[adesso2][secondIndex]->GetTOT();
		interestingADC[1] = selectedDigit[adesso2][secondIndex]->GetADC();
		interestingWeight[1] = selectedDigit[adesso2][secondIndex]->GetWeight();
		Int_t vol2[5]; for(jj=0; jj<5; jj++) vol2[jj] = selectedDigit[adesso2][secondIndex]->GetDetectorIndex(jj);
		AliDebug(1,Form(" %2d %1d %2d %1d %2d", vol2[0], vol2[1], vol2[2], vol2[4], vol2[3]));
		volDum = vol2[3];
		vol2[3] = vol2[4];
		vol2[4] = volDum;
		fTOFGeometry->GetPosPar(vol2,pos);
		interestingX[1] = pos[0];
		interestingY[1] = pos[1];
		interestingZ[1] = pos[2];

		interestingTOF[2] = selectedDigit[adesso3][thirdIndex]->GetTDC();
		interestingTOT[2] = selectedDigit[adesso3][thirdIndex]->GetTOT();
		interestingADC[2] = selectedDigit[adesso3][thirdIndex]->GetADC();
		interestingWeight[2] = selectedDigit[adesso3][thirdIndex]->GetWeight();
		Int_t vol3[5]; for(jj=0; jj<5; jj++) vol3[jj] = selectedDigit[adesso3][thirdIndex]->GetDetectorIndex(jj);
		AliDebug(1,Form(" %2d %1d %2d %1d %2d", vol3[0], vol3[1], vol3[2], vol3[4], vol3[3]));
		volDum = vol3[3];
		vol3[3] = vol3[4];
		vol3[4] = volDum;
		fTOFGeometry->GetPosPar(vol3,pos);
		interestingX[2] = pos[0];
		interestingY[2] = pos[1];
		interestingZ[2] = pos[2];


		AverageCalculations(3, interestingX, interestingY, interestingZ,
				    interestingTOF, interestingTOT, interestingADC,
				    interestingWeight,
				    parTOF, posClus, check);


		for (jj=0; jj<5; jj++) det[jj] = -1;
		for (jj=0; jj<3; jj++) posF[jj] = posClus[jj];
		fTOFGeometry->GetDetID(posF, det);

		volIdClus = fTOFGeometry->GetAliSensVolIndex(det[0],det[1],det[2]);
		//volIdClus = GetClusterVolIndex(det);

		for (jj=3; jj<11; jj++) padsCluster[jj] = -1;
		padsCluster[3] = selectedDigit[adesso1][firstIndex]->GetDetectorIndex(4);
		padsCluster[4] = selectedDigit[adesso1][firstIndex]->GetDetectorIndex(3);
		padsCluster[5] = selectedDigit[adesso2][secondIndex]->GetDetectorIndex(4);
		padsCluster[6] = selectedDigit[adesso2][secondIndex]->GetDetectorIndex(3);
		padsCluster[7] = selectedDigit[adesso3][thirdIndex]->GetDetectorIndex(4);
		padsCluster[8] = selectedDigit[adesso3][thirdIndex]->GetDetectorIndex(3);

		for (jj=0; jj<6; jj++) covClus[jj] = 0.;
		Int_t ** indDet = new Int_t*[3];
		for (jj=0; jj<3; jj++) indDet[jj] = new Int_t [5];
		for (jj=0; jj<3; jj++) indDet[jj][0] = nSector;
		for (jj=0; jj<3; jj++) indDet[jj][1] = nPlate;
		for (jj=0; jj<3; jj++) indDet[jj][2] = nStrip;
		for (jj=0; jj<3; jj++) indDet[jj][3] = padsCluster[2*jj+3];
		for (jj=0; jj<3; jj++) indDet[jj][4] = padsCluster[2*jj+1+3];
		GetClusterPars(/*check,*/ 3, indDet, interestingWeight, posClus, covClus);
		for (jj=0; jj<3; jj++) delete [] indDet[jj];
		delete [] indDet;

		// To fill the track index array
		dummyCounter=-1;
		for (jj=0; jj<kMaxNumberOfTracksPerDigit; jj++) tracks[jj] = -1;
		for (Int_t kk=0; kk<kMaxNumberOfTracksPerDigit; kk++) { // three is the max number of tracks associated to one digit
		  if (selectedDigit[adesso1][firstIndex]->GetTrackLabel(kk)==-1) continue;
		  else {
		    dummyCounter++;
		    tracks[dummyCounter] = selectedDigit[adesso1][firstIndex]->GetTrackLabel(kk);
		  }
		}
		for (Int_t kk=0; kk<kMaxNumberOfTracksPerDigit; kk++) { // three is the max number of tracks associated to one digit
		  if (selectedDigit[adesso2][secondIndex]->GetTrackLabel(kk)==-1) continue;
		  else {

		    alreadyStored = kFALSE;
		    for (jj=0; jj<dummyCounter+1; jj++)
		      alreadyStored = alreadyStored || (tracks[jj]!=-1 && tracks[jj]==selectedDigit[adesso2][secondIndex]->GetTrackLabel(kk));

		    if (alreadyStored) continue;
		    if (dummyCounter==2) { // three is the max number of tracks associated to one cluster
		      AliWarning("  Siamo al limite!");
		      continue;
		    }

		    dummyCounter++;
		    tracks[dummyCounter] = selectedDigit[adesso2][secondIndex]->GetTrackLabel(kk);

		  }

		}
		for (Int_t kk=0; kk<kMaxNumberOfTracksPerDigit; kk++) { // three is the max number of tracks associated to one digit
		  if (selectedDigit[adesso3][thirdIndex]->GetTrackLabel(kk)==-1) continue;
		  else {

		    alreadyStored = kFALSE;
		    for (jj=0; jj<dummyCounter+1; jj++)
		      alreadyStored = alreadyStored || (tracks[jj]!=-1 && tracks[jj]==selectedDigit[adesso3][thirdIndex]->GetTrackLabel(kk));

		    if (alreadyStored) continue;
		    if (dummyCounter==2) { // three is the max number of tracks associated to one cluster
		      AliWarning("  Siamo al limite!");
		      continue;
		    }

		    dummyCounter++;
		    tracks[dummyCounter] = selectedDigit[adesso3][thirdIndex]->GetTrackLabel(kk);

		  }

		}


		AliTOFcluster *tofCluster =
		  new AliTOFcluster(volIdClus, posClus[0], posClus[1], posClus[2],
				    covClus[0], covClus[1], covClus[2], covClus[3], covClus[4], covClus[5],
				    tracks, det, parTOF, status, selectedDigit[adesso1][firstIndex]->GetIndex()); // to be updated
		InsertCluster(tofCluster);

		AliDebug(2, Form("       %4d  %f %f %f  %f %f %f %f %f %f  %3d %3d %3d  %2d %1d %2d %1d %2d  %4d %3d %3d %4d %4d  %1d  %4d", 
				 volIdClus, posClus[0], posClus[1], posClus[2],
				 covClus[0], covClus[1], covClus[2], covClus[3], covClus[4], covClus[5],
				 tracks[0], tracks[1], tracks[2],
				 det[0], det[1], det[2], det[3], det[4],
				 parTOF[0], parTOF[1], parTOF[2], parTOF[3], parTOF[4],
				 status, selectedDigit[adesso1][firstIndex]->GetIndex()));


		volDum = vol1[3];
		vol1[3] = vol1[4];
		vol1[4] = volDum;
		fTOFdigitMap->ResetDigitNumber(vol1,selectedDigit[adesso1][firstIndex]->GetIndex());
		volDum = vol2[3];
		vol2[3] = vol2[4];
		vol2[4] = volDum;
		fTOFdigitMap->ResetDigitNumber(vol2,selectedDigit[adesso2][secondIndex]->GetIndex());
		volDum = vol3[3];
		vol3[3] = vol3[4];
		vol3[4] = volDum;
		fTOFdigitMap->ResetDigitNumber(vol3,selectedDigit[adesso3][thirdIndex]->GetIndex());


	      } // close loop on third digit
	    } // close loop on adesso3

	  } // close loop on second digit
	} // close loop on adesso2

      } // close loop on first digit
    } // close loop on adesso1

    for (ii=0; ii<kMaxNumberOfInterestingPads; ii++)
      for (jj=0; jj<kMaxNumberOfDigitsPerVolume; jj++)
	selectedDigit[ii][jj] = 0x0;

  } // loop on iPad

  for (ii=0; ii<kMaxNumberOfInterestingPads; ii++) {
    for (jj=0; jj<kMaxNumberOfDigitsPerVolume; jj++) {
      delete [] selectedDigit[ii][jj];
      selectedDigit[ii][jj] = 0x0;
    }
    delete [] selectedDigit[ii];
    selectedDigit[ii] = 0x0;
  }
  delete [] selectedDigit;
  selectedDigit = 0x0;

}
//_____________________________________________________________________________

void AliTOFClusterFinderV1::FindClusters23(Int_t nSector,
					   Int_t nPlate,
					   Int_t nStrip)
{
  //
  // This function searches the neighbouring digits (stored in the fDigits object),
  // to perform clusters (stored in the fTofClusters array).
  //
  // This research has been made by checking the fTOFdigitMap object,
  // filled at digits/raw-data reading time.
  //

  const Int_t kMaxNumberOfInterestingPads = 4;
  const Int_t kMaxNumberOfTracksPerDigit = 3;
  const Int_t kMaxNumberOfDigitsPerVolume = 10;

  Int_t ii = 0;

  Int_t digitsInVolumeIndices[kMaxNumberOfDigitsPerVolume];
  for (ii=0; ii<kMaxNumberOfDigitsPerVolume; ii++)
    digitsInVolumeIndices[ii] = -1;

  Int_t vol[5] = {nSector,nPlate,nStrip,-1,-1};

  Float_t pos[3] = {0.,0.,0.};

  Int_t jj = 0;
  Int_t interestingPadX[kMaxNumberOfInterestingPads];
  for (jj=0; jj<kMaxNumberOfInterestingPads; jj++) interestingPadX[jj] = -1;
  Int_t interestingPadZ[kMaxNumberOfInterestingPads];
  for (jj=0; jj<kMaxNumberOfInterestingPads; jj++) interestingPadZ[jj] = -1;
  Double_t interestingTOT[kMaxNumberOfInterestingPads];
  for (jj=0; jj<kMaxNumberOfInterestingPads; jj++) interestingTOT[jj] = 0;
  Double_t interestingADC[kMaxNumberOfInterestingPads];
  for (jj=0; jj<kMaxNumberOfInterestingPads; jj++) interestingADC[jj] = 0;
  Double_t interestingTOF[kMaxNumberOfInterestingPads];
  for (jj=0; jj<kMaxNumberOfInterestingPads; jj++) interestingTOF[jj] = 0;
  Double_t interestingWeight[kMaxNumberOfInterestingPads];
  for (jj=0; jj<kMaxNumberOfInterestingPads; jj++) interestingWeight[jj] = 0;

  Float_t interestingX[kMaxNumberOfInterestingPads];
  for (jj=0; jj<kMaxNumberOfInterestingPads; jj++) interestingX[jj] = 0;
  Float_t interestingY[kMaxNumberOfInterestingPads];
  for (jj=0; jj<kMaxNumberOfInterestingPads; jj++) interestingY[jj] = 0;
  Float_t interestingZ[kMaxNumberOfInterestingPads];
  for (jj=0; jj<kMaxNumberOfInterestingPads; jj++) interestingZ[jj] = 0;

  Float_t interDigit[kMaxNumberOfInterestingPads];
  for (jj=0; jj<kMaxNumberOfInterestingPads; jj++) interDigit[jj] = 0;

  Int_t padsCluster[11];
  padsCluster[0] = nSector;
  padsCluster[1] = nPlate;
  padsCluster[2] = nStrip;
  for (jj=3; jj<11; jj++) padsCluster[jj] = -1;

  Int_t interestingCounter=-1;
  Int_t digitIndexLocal = -1;
  Int_t iPad  = -1;
  Int_t iPadX = -1;
  Int_t iPadZ = -1;

  Bool_t check = kFALSE;
  Int_t parTOF[7];
  for (jj=0; jj<7; jj++) parTOF[jj] = 0;
  Double_t posClus[3];
  for (jj=0; jj<3; jj++) posClus[jj] = 0.;
  Int_t det[5];
  for (jj=0; jj<5; jj++) det[jj] = -1;
  Float_t posF[3];
  for (jj=0; jj<3; jj++) posF[jj] = 0.;
  UShort_t volIdClus = 0;
  Bool_t status = kTRUE; //assume all sim channels ok in the beginning...
  Double_t covClus[6];
  for (jj=0; jj<6; jj++) covClus[jj] = 0.;
  Int_t tracks[kMaxNumberOfTracksPerDigit];
  for (jj=0; jj<kMaxNumberOfTracksPerDigit; jj++) tracks[jj] = -1;
  Int_t dummyCounter=-1;
  Bool_t alreadyStored = kFALSE;

  AliTOFselectedDigit ***selectedDigit = new AliTOFselectedDigit**[kMaxNumberOfInterestingPads];
  for (ii=0; ii<kMaxNumberOfInterestingPads; ii++)
    selectedDigit[ii] = new AliTOFselectedDigit*[kMaxNumberOfDigitsPerVolume];

  for (ii=0; ii<kMaxNumberOfInterestingPads; ii++)
    for (jj=0; jj<kMaxNumberOfDigitsPerVolume; jj++) selectedDigit[ii][jj] = 0x0;

  AliTOFdigit *digitInteresting;

  for (iPad=0; iPad<AliTOFGeometry::NpadZ()*AliTOFGeometry::NpadX()-3; iPad+=2) {

    iPadZ = iPad%AliTOFGeometry::NpadZ(); //iPad%2;
    iPadX = iPad/AliTOFGeometry::NpadZ(); //iPad/2;

    AliDebug(3, Form("%2d %1d %2d\n", iPad, iPadZ, iPadX));






    interestingCounter=-1;

    vol[4] = iPadZ  , vol[3]  = iPadX;
    if (fTOFdigitMap->GetNumberOfDigits(vol)>0)
      interestingCounter++;

    vol[4] = iPadZ, vol[3] = iPadX+1;
    if (fTOFdigitMap->GetNumberOfDigits(vol)>0)
      interestingCounter++;

    vol[4] = AliTOFGeometry::NpadZ()-(iPadZ+1), vol[3] = iPadX;
    if (fTOFdigitMap->GetNumberOfDigits(vol)>0)
      interestingCounter++;

    vol[4] = AliTOFGeometry::NpadZ()-(iPadZ+1), vol[3] = iPadX+1;
    if (fTOFdigitMap->GetNumberOfDigits(vol)>0)
      interestingCounter++;

    if (interestingCounter+1!=3) continue; // the hit pads have to be 3
    else interestingCounter=-1;


    for (ii=0; ii<kMaxNumberOfInterestingPads; ii++)
      for (jj=0; jj<kMaxNumberOfDigitsPerVolume; jj++)
	selectedDigit[ii][jj] = 0x0;


    vol[4] = iPadZ, vol[3] = iPadX;

    if (fTOFdigitMap->GetNumberOfDigits(vol)>0) {
      interestingCounter++;
      for (ii=0; ii<kMaxNumberOfDigitsPerVolume; ii++)
	digitsInVolumeIndices[ii] = -1;
      fTOFdigitMap->GetDigitIndex(vol, digitsInVolumeIndices);
      digitIndexLocal=-1;
      for(Int_t digIndex=0; digIndex<kMaxNumberOfDigitsPerVolume; digIndex++) {
	if (digitsInVolumeIndices[digIndex]<0) continue;
	digitInteresting = (AliTOFdigit*)fDigits->UncheckedAt(digitsInVolumeIndices[digIndex]);
	if (digitInteresting->GetToT()<=0) continue; // AdC
	digitIndexLocal++; // AdC

	AliDebug(1,Form(" %2d %1d %2d %1d %2d  %d %d %d %d  %5d  %5d %5d %5d",
			vol[0], vol[1], vol[2] ,vol[4], vol[3],
			digitInteresting->GetTdc(), digitInteresting->GetAdc(),
			digitInteresting->GetToT(), digitInteresting->GetToT()*digitInteresting->GetToT(),
			digitsInVolumeIndices[digIndex],
			digitInteresting->GetTrack(0), digitInteresting->GetTrack(1), digitInteresting->GetTrack(2)));

	selectedDigit[interestingCounter][digitIndexLocal] = new // AdC
	  AliTOFselectedDigit(vol, digitInteresting->GetTdc(),
			      digitInteresting->GetAdc(), digitInteresting->GetToT(),
			      digitInteresting->GetToT()*digitInteresting->GetToT(),
			      digitsInVolumeIndices[digIndex],
			      digitInteresting->GetTracks());
      }
      if (digitIndexLocal==-1) interestingCounter--; // AdC
    }


    vol[4] = iPadZ, vol[3] = iPadX+1;

    if (fTOFdigitMap->GetNumberOfDigits(vol)>0) {
      interestingCounter++;
      for (ii=0; ii<kMaxNumberOfDigitsPerVolume; ii++)
	digitsInVolumeIndices[ii] = -1;
      fTOFdigitMap->GetDigitIndex(vol, digitsInVolumeIndices);
      digitIndexLocal=-1; // AdC
      for(Int_t digIndex=0; digIndex<kMaxNumberOfDigitsPerVolume; digIndex++) {
	if (digitsInVolumeIndices[digIndex]<0) continue;
	digitInteresting = (AliTOFdigit*)fDigits->UncheckedAt(digitsInVolumeIndices[digIndex]);
	if (digitInteresting->GetToT()<=0) continue; // AdC
	digitIndexLocal++; // AdC

	AliDebug(1,Form(" %2d %1d %2d %1d %2d  %d %d %d %d  %5d  %5d %5d %5d",
			vol[0], vol[1], vol[2] ,vol[4], vol[3],
			digitInteresting->GetTdc(), digitInteresting->GetAdc(),
			digitInteresting->GetToT(), digitInteresting->GetToT()*digitInteresting->GetToT(),
			digitsInVolumeIndices[digIndex],
			digitInteresting->GetTrack(0), digitInteresting->GetTrack(1), digitInteresting->GetTrack(2)));

	selectedDigit[interestingCounter][digitIndexLocal] = new // AdC
	  AliTOFselectedDigit(vol, digitInteresting->GetTdc(),
			      digitInteresting->GetAdc(), digitInteresting->GetToT(),
			      digitInteresting->GetToT()*digitInteresting->GetToT(),
			      digitsInVolumeIndices[digIndex],
			      digitInteresting->GetTracks());
      }
      if (digitIndexLocal==-1) interestingCounter--; // AdC
    }


    vol[4] = AliTOFGeometry::NpadZ()-(iPadZ+1), vol[3] = iPadX;

    if (fTOFdigitMap->GetNumberOfDigits(vol)>0) {
      interestingCounter++;
      for (ii=0; ii<kMaxNumberOfDigitsPerVolume; ii++)
	digitsInVolumeIndices[ii] = -1;
      fTOFdigitMap->GetDigitIndex(vol, digitsInVolumeIndices);
      digitIndexLocal=-1; // AdC
      for(Int_t digIndex=0; digIndex<kMaxNumberOfDigitsPerVolume; digIndex++) {
	if (digitsInVolumeIndices[digIndex]<0) continue;
	digitInteresting = (AliTOFdigit*)fDigits->UncheckedAt(digitsInVolumeIndices[digIndex]);
	if (digitInteresting->GetToT()<=0) continue; // AdC
	digitIndexLocal++; // AdC

	AliDebug(1,Form(" %2d %1d %2d %1d %2d  %d %d %d %d  %5d  %5d %5d %5d",
			vol[0], vol[1], vol[2] ,vol[4], vol[3],
			digitInteresting->GetTdc(), digitInteresting->GetAdc(),
			digitInteresting->GetToT(), digitInteresting->GetToT()*digitInteresting->GetToT(),
			digitsInVolumeIndices[digIndex],
			digitInteresting->GetTrack(0), digitInteresting->GetTrack(1), digitInteresting->GetTrack(2)));

	selectedDigit[interestingCounter][digitIndexLocal] = new // AdC
	  AliTOFselectedDigit(vol, digitInteresting->GetTdc(),
			      digitInteresting->GetAdc(), digitInteresting->GetToT(),
			      digitInteresting->GetToT()*digitInteresting->GetToT(),
			      digitsInVolumeIndices[digIndex],
			      digitInteresting->GetTracks());
      }
      if (digitIndexLocal==-1) interestingCounter--; // AdC
    }


    vol[4] = AliTOFGeometry::NpadZ()-(iPadZ+1), vol[3] = iPadX+1;

    if (fTOFdigitMap->GetNumberOfDigits(vol)>0) {
      interestingCounter++;
      for (ii=0; ii<kMaxNumberOfDigitsPerVolume; ii++)
	digitsInVolumeIndices[ii] = -1;
      fTOFdigitMap->GetDigitIndex(vol, digitsInVolumeIndices);
      digitIndexLocal=-1; // AdC
      for(Int_t digIndex=0; digIndex<kMaxNumberOfDigitsPerVolume; digIndex++) {
	if (digitsInVolumeIndices[digIndex]<0) continue;
	digitInteresting = (AliTOFdigit*)fDigits->UncheckedAt(digitsInVolumeIndices[digIndex]);
	if (digitInteresting->GetToT()<=0) continue; // AdC
	digitIndexLocal++; // AdC

	AliDebug(1,Form(" %2d %1d %2d %1d %2d  %d %d %d %d  %5d  %5d %5d %5d",
			vol[0], vol[1], vol[2] ,vol[4], vol[3],
			digitInteresting->GetTdc(), digitInteresting->GetAdc(),
			digitInteresting->GetToT(), digitInteresting->GetToT()*digitInteresting->GetToT(),
			digitsInVolumeIndices[digIndex],
			digitInteresting->GetTrack(0), digitInteresting->GetTrack(1), digitInteresting->GetTrack(2)));

	selectedDigit[interestingCounter][digitIndexLocal] = new // AdC
	  AliTOFselectedDigit(vol, digitInteresting->GetTdc(),
			      digitInteresting->GetAdc(), digitInteresting->GetToT(),
			      digitInteresting->GetToT()*digitInteresting->GetToT(),
			      digitsInVolumeIndices[digIndex],
			      digitInteresting->GetTracks());
      }
      if (digitIndexLocal==-1) interestingCounter--; // AdC
    }

    AliDebug(1,Form(" e adesso %1d", interestingCounter+1));

    for (Int_t adesso1=0; adesso1<interestingCounter+1; adesso1++) {
      for (Int_t firstIndex=0; firstIndex<kMaxNumberOfDigitsPerVolume; firstIndex++) {
	if (selectedDigit[adesso1][firstIndex]==0x0) continue;

	for (Int_t adesso2=adesso1+1; adesso2<interestingCounter+1; adesso2++) {
	  for (Int_t secondIndex=0; secondIndex<kMaxNumberOfDigitsPerVolume; secondIndex++) {
	    if (selectedDigit[adesso2][secondIndex]==0x0) continue;

	    if (TMath::Abs(selectedDigit[adesso1][firstIndex]->GetTDC()-selectedDigit[adesso2][secondIndex]->GetTDC())>fMaxDeltaTime) continue;

	    interestingTOF[0] = selectedDigit[adesso1][firstIndex]->GetTDC();
	    interestingTOT[0] = selectedDigit[adesso1][firstIndex]->GetTOT();
	    interestingADC[0] = selectedDigit[adesso1][firstIndex]->GetADC();
	    interestingWeight[0] = selectedDigit[adesso1][firstIndex]->GetWeight();
	    Int_t vol1[5]; for(jj=0; jj<5; jj++) vol1[jj] = selectedDigit[adesso1][firstIndex]->GetDetectorIndex(jj);
	    AliDebug(1,Form(" %2d %1d %2d %1d %2d", vol1[0], vol1[1], vol1[2], vol1[4], vol1[3]));
	    Int_t volDum = vol1[3];
	    vol1[3] = vol1[4];
	    vol1[4] = volDum;
	    fTOFGeometry->GetPosPar(vol1,pos);
	    interestingX[0] = pos[0];
	    interestingY[0] = pos[1];
	    interestingZ[0] = pos[2];

	    interestingTOF[1] = selectedDigit[adesso2][secondIndex]->GetTDC();
	    interestingTOT[1] = selectedDigit[adesso2][secondIndex]->GetTOT();
	    interestingADC[1] = selectedDigit[adesso2][secondIndex]->GetADC();
	    interestingWeight[1] = selectedDigit[adesso2][secondIndex]->GetWeight();
	    Int_t vol2[5]; for(jj=0; jj<5; jj++) vol2[jj] = selectedDigit[adesso2][secondIndex]->GetDetectorIndex(jj);
	    AliDebug(1,Form(" %2d %1d %2d %1d %2d", vol2[0], vol2[1], vol2[2], vol2[4], vol2[3]));
	    volDum = vol2[3];
	    vol2[3] = vol2[4];
	    vol2[4] = volDum;
	    fTOFGeometry->GetPosPar(vol2,pos);
	    interestingX[1] = pos[0];
	    interestingY[1] = pos[1];
	    interestingZ[1] = pos[2];

	    AverageCalculations(2, interestingX, interestingY, interestingZ,
				interestingTOF, interestingTOT, interestingADC,
				interestingWeight,
				parTOF, posClus, check);

	    for (jj=0; jj<5; jj++) det[jj] = -1;
	    for (jj=0; jj<3; jj++) posF[jj] = posClus[jj];
	    fTOFGeometry->GetDetID(posF, det);

	    volIdClus = fTOFGeometry->GetAliSensVolIndex(det[0],det[1],det[2]);
	    //volIdClus = GetClusterVolIndex(det);

	    for (jj=3; jj<11; jj++) padsCluster[jj] = -1;
	    padsCluster[3] = selectedDigit[adesso1][firstIndex]->GetDetectorIndex(4);
	    padsCluster[4] = selectedDigit[adesso1][firstIndex]->GetDetectorIndex(3);
	    padsCluster[5] = selectedDigit[adesso2][secondIndex]->GetDetectorIndex(4);
	    padsCluster[6] = selectedDigit[adesso2][secondIndex]->GetDetectorIndex(3);

	    for (jj=0; jj<6; jj++) covClus[jj] = 0.;
	    Int_t ** indDet = new Int_t*[2];
	    for (jj=0; jj<2; jj++) indDet[jj] = new Int_t [5];
	    for (jj=0; jj<2; jj++) indDet[jj][0] = nSector;
	    for (jj=0; jj<2; jj++) indDet[jj][1] = nPlate;
	    for (jj=0; jj<2; jj++) indDet[jj][2] = nStrip;
	    for (jj=0; jj<2; jj++) indDet[jj][3] = padsCluster[2*jj+3];
	    for (jj=0; jj<2; jj++) indDet[jj][4] = padsCluster[2*jj+1+3];
	    GetClusterPars(/*check,*/ 2, indDet, interestingWeight, posClus, covClus);
	    for (jj=0; jj<2; jj++) delete [] indDet[jj];
	    delete [] indDet;

	    // To fill the track index array
	    dummyCounter=-1;
	    for (jj=0; jj<kMaxNumberOfTracksPerDigit; jj++) tracks[jj] = -1;
	    for (Int_t kk=0; kk<kMaxNumberOfTracksPerDigit; kk++) { // three is the max number of tracks associated to one digit
	      if (selectedDigit[adesso1][firstIndex]->GetTrackLabel(kk)==-1) continue;
	      else {
		dummyCounter++;
		tracks[dummyCounter] = selectedDigit[adesso1][firstIndex]->GetTrackLabel(kk);
	      }
	    }
	    for (Int_t kk=0; kk<kMaxNumberOfTracksPerDigit; kk++) { // three is the max number of tracks associated to one digit
	      if (selectedDigit[adesso2][secondIndex]->GetTrackLabel(kk)==-1) continue;
	      else {

		alreadyStored = kFALSE;
		for (jj=0; jj<dummyCounter+1; jj++)
		  alreadyStored = alreadyStored || (tracks[jj]!=-1 && tracks[jj]==selectedDigit[adesso2][secondIndex]->GetTrackLabel(kk));

		if (alreadyStored) continue;
		if (dummyCounter==2) { // three is the max number of tracks associated to one cluster
		  AliWarning("  Siamo al limite!");
		  continue;
		}

		dummyCounter++;
		tracks[dummyCounter] = selectedDigit[adesso2][secondIndex]->GetTrackLabel(kk);

	      }

	    }


	    AliTOFcluster *tofCluster =
	      new AliTOFcluster(volIdClus, posClus[0], posClus[1], posClus[2],
				covClus[0], covClus[1], covClus[2], covClus[3], covClus[4], covClus[5],
				tracks, det, parTOF, status, selectedDigit[adesso1][firstIndex]->GetIndex()); // to be updated
	    InsertCluster(tofCluster);

	    AliDebug(2, Form("       %4d  %f %f %f  %f %f %f %f %f %f  %3d %3d %3d  %2d %1d %2d %1d %2d  %4d %3d %3d %4d %4d  %1d  %4d", 
			     volIdClus, posClus[0], posClus[1], posClus[2],
			     covClus[0], covClus[1], covClus[2], covClus[3], covClus[4], covClus[5],
			     tracks[0], tracks[1], tracks[2],
			     det[0], det[1], det[2], det[3], det[4],
			     parTOF[0], parTOF[1], parTOF[2], parTOF[3], parTOF[4],
			     status, selectedDigit[adesso1][firstIndex]->GetIndex()));

	    volDum = vol1[3];
	    vol1[3] = vol1[4];
	    vol1[4] = volDum;
	    fTOFdigitMap->ResetDigitNumber(vol1,selectedDigit[adesso1][firstIndex]->GetIndex());
	    volDum = vol2[3];
	    vol2[3] = vol2[4];
	    vol2[4] = volDum;
	    fTOFdigitMap->ResetDigitNumber(vol2,selectedDigit[adesso2][secondIndex]->GetIndex());


	  } // close loop on second digit
	} // close loop on adesso2

      } // close loop on first digit
    } // close loop on adesso1

    for (ii=0; ii<kMaxNumberOfInterestingPads; ii++)
      for (jj=0; jj<kMaxNumberOfDigitsPerVolume; jj++)
	selectedDigit[ii][jj] = 0x0;

  } // loop on iPad

  for (ii=0; ii<kMaxNumberOfInterestingPads; ii++) {
    for (jj=0; jj<kMaxNumberOfDigitsPerVolume; jj++) {
      delete [] selectedDigit[ii][jj];
      selectedDigit[ii][jj] = 0x0;
    }
    delete [] selectedDigit[ii];
    selectedDigit[ii] = 0x0;
  }
  delete [] selectedDigit;
  selectedDigit = 0x0;

}
//_____________________________________________________________________________

void AliTOFClusterFinderV1::FindClusters24(Int_t nSector,
					   Int_t nPlate,
					   Int_t nStrip)
{
  //
  // This function searches the neighbouring digits (stored in the fDigits object),
  // to perform clusters (stored in the fTofClusters array).
  //
  // This research has been made by checking the fTOFdigitMap object,
  // filled at digits/raw-data reading time.
  //

  const Int_t kMaxNumberOfInterestingPads = 4;
  const Int_t kMaxNumberOfTracksPerDigit = 3;
  const Int_t kMaxNumberOfDigitsPerVolume = 10;

  Int_t ii = 0;

  Int_t digitsInVolumeIndices[kMaxNumberOfDigitsPerVolume];
  for (ii=0; ii<kMaxNumberOfDigitsPerVolume; ii++)
    digitsInVolumeIndices[ii] = -1;

  Int_t vol[5] = {nSector,nPlate,nStrip,-1,-1};

  Float_t pos[3] = {0.,0.,0.};

  Int_t jj = 0;
  Int_t interestingPadX[kMaxNumberOfInterestingPads];
  for (jj=0; jj<kMaxNumberOfInterestingPads; jj++) interestingPadX[jj] = -1;
  Int_t interestingPadZ[kMaxNumberOfInterestingPads];
  for (jj=0; jj<kMaxNumberOfInterestingPads; jj++) interestingPadZ[jj] = -1;
  Double_t interestingTOT[kMaxNumberOfInterestingPads];
  for (jj=0; jj<kMaxNumberOfInterestingPads; jj++) interestingTOT[jj] = 0;
  Double_t interestingADC[kMaxNumberOfInterestingPads];
  for (jj=0; jj<kMaxNumberOfInterestingPads; jj++) interestingADC[jj] = 0;
  Double_t interestingTOF[kMaxNumberOfInterestingPads];
  for (jj=0; jj<kMaxNumberOfInterestingPads; jj++) interestingTOF[jj] = 0;
  Double_t interestingWeight[kMaxNumberOfInterestingPads];
  for (jj=0; jj<kMaxNumberOfInterestingPads; jj++) interestingWeight[jj] = 0;

  Float_t interestingX[kMaxNumberOfInterestingPads];
  for (jj=0; jj<kMaxNumberOfInterestingPads; jj++) interestingX[jj] = 0;
  Float_t interestingY[kMaxNumberOfInterestingPads];
  for (jj=0; jj<kMaxNumberOfInterestingPads; jj++) interestingY[jj] = 0;
  Float_t interestingZ[kMaxNumberOfInterestingPads];
  for (jj=0; jj<kMaxNumberOfInterestingPads; jj++) interestingZ[jj] = 0;

  Float_t interDigit[kMaxNumberOfInterestingPads];
  for (jj=0; jj<kMaxNumberOfInterestingPads; jj++) interDigit[jj] = 0;

  Int_t padsCluster[11];
  padsCluster[0] = nSector;
  padsCluster[1] = nPlate;
  padsCluster[2] = nStrip;
  for (jj=3; jj<11; jj++) padsCluster[jj] = -1;

  Int_t interestingCounter=-1;
  Int_t digitIndexLocal = -1;
  Int_t iPad  = -1;
  Int_t iPadX = -1;
  Int_t iPadZ = -1;

  Bool_t check = kFALSE;
  Int_t parTOF[7];
  for (jj=0; jj<7; jj++) parTOF[jj] = 0;
  Double_t posClus[3];
  for (jj=0; jj<3; jj++) posClus[jj] = 0.;
  Int_t det[5];
  for (jj=0; jj<5; jj++) det[jj] = -1;
  Float_t posF[3];
  for (jj=0; jj<3; jj++) posF[jj] = 0.;
  UShort_t volIdClus = 0;
  Bool_t status = kTRUE; //assume all sim channels ok in the beginning...
  Double_t covClus[6];
  for (jj=0; jj<6; jj++) covClus[jj] = 0.;
  Int_t tracks[kMaxNumberOfTracksPerDigit];
  for (jj=0; jj<kMaxNumberOfTracksPerDigit; jj++) tracks[jj] = -1;
  Int_t dummyCounter=-1;
  Bool_t alreadyStored = kFALSE;

  AliTOFselectedDigit ***selectedDigit = new AliTOFselectedDigit**[kMaxNumberOfInterestingPads];
  for (ii=0; ii<kMaxNumberOfInterestingPads; ii++)
    selectedDigit[ii] = new AliTOFselectedDigit*[kMaxNumberOfDigitsPerVolume];

  for (ii=0; ii<kMaxNumberOfInterestingPads; ii++)
    for (jj=0; jj<kMaxNumberOfDigitsPerVolume; jj++) selectedDigit[ii][jj] = 0x0;

  AliTOFdigit *digitInteresting;

  for (iPad=0; iPad<AliTOFGeometry::NpadZ()*AliTOFGeometry::NpadX()-3; iPad+=2) {

    iPadZ = iPad%AliTOFGeometry::NpadZ(); //iPad%2;
    iPadX = iPad/AliTOFGeometry::NpadZ(); //iPad/2;

    AliDebug(3, Form("%2d %1d %2d\n", iPad, iPadZ, iPadX));






    interestingCounter=-1;

    vol[4] = iPadZ  , vol[3]  = iPadX;
    if (fTOFdigitMap->GetNumberOfDigits(vol)>0)
      interestingCounter++;

    vol[4] = iPadZ, vol[3] = iPadX+1;
    if (fTOFdigitMap->GetNumberOfDigits(vol)>0)
      interestingCounter++;

    vol[4] = AliTOFGeometry::NpadZ()-(iPadZ+1), vol[3] = iPadX;
    if (fTOFdigitMap->GetNumberOfDigits(vol)>0)
      interestingCounter++;

    vol[4] = AliTOFGeometry::NpadZ()-(iPadZ+1), vol[3] = iPadX+1;
    if (fTOFdigitMap->GetNumberOfDigits(vol)>0)
      interestingCounter++;

    if (interestingCounter+1!=4) continue; // the hit pads have to be 4
    else interestingCounter=-1;


    for (ii=0; ii<kMaxNumberOfInterestingPads; ii++)
      for (jj=0; jj<kMaxNumberOfDigitsPerVolume; jj++)
	selectedDigit[ii][jj] = 0x0;


    vol[4] = iPadZ, vol[3] = iPadX;

    if (fTOFdigitMap->GetNumberOfDigits(vol)>0) {
      interestingCounter++;
      for (ii=0; ii<kMaxNumberOfDigitsPerVolume; ii++)
	digitsInVolumeIndices[ii] = -1;
      fTOFdigitMap->GetDigitIndex(vol, digitsInVolumeIndices);
      digitIndexLocal=-1; // AdC
      for(Int_t digIndex=0; digIndex<kMaxNumberOfDigitsPerVolume; digIndex++) {
	if (digitsInVolumeIndices[digIndex]<0) continue;
	digitInteresting = (AliTOFdigit*)fDigits->UncheckedAt(digitsInVolumeIndices[digIndex]);
	if (digitInteresting->GetToT()<=0) continue; // AdC
	digitIndexLocal++; // AdC

	AliDebug(1,Form(" %2d %1d %2d %1d %2d  %d %d %d %d  %5d  %5d %5d %5d",
			vol[0], vol[1], vol[2] ,vol[4], vol[3],
			digitInteresting->GetTdc(), digitInteresting->GetAdc(),
			digitInteresting->GetToT(), digitInteresting->GetToT()*digitInteresting->GetToT(),
			digitsInVolumeIndices[digIndex],
			digitInteresting->GetTrack(0), digitInteresting->GetTrack(1), digitInteresting->GetTrack(2)));

	selectedDigit[interestingCounter][digitIndexLocal] = new // AdC
	  AliTOFselectedDigit(vol, digitInteresting->GetTdc(),
			      digitInteresting->GetAdc(), digitInteresting->GetToT(),
			      digitInteresting->GetToT()*digitInteresting->GetToT(),
			      digitsInVolumeIndices[digIndex],
			      digitInteresting->GetTracks());
      }
      if (digitIndexLocal==-1) interestingCounter--; // AdC
    }


    vol[4] = iPadZ, vol[3] = iPadX+1;

    if (fTOFdigitMap->GetNumberOfDigits(vol)>0) {
      interestingCounter++;
      for (ii=0; ii<kMaxNumberOfDigitsPerVolume; ii++)
	digitsInVolumeIndices[ii] = -1;
      fTOFdigitMap->GetDigitIndex(vol, digitsInVolumeIndices);
      digitIndexLocal=-1; // AdC
      for(Int_t digIndex=0; digIndex<kMaxNumberOfDigitsPerVolume; digIndex++) {
	if (digitsInVolumeIndices[digIndex]<0) continue;
	digitInteresting = (AliTOFdigit*)fDigits->UncheckedAt(digitsInVolumeIndices[digIndex]);
	if (digitInteresting->GetToT()<=0) continue; // AdC
	digitIndexLocal++; // AdC

	AliDebug(1,Form(" %2d %1d %2d %1d %2d  %d %d %d %d  %5d  %5d %5d %5d",
			vol[0], vol[1], vol[2] ,vol[4], vol[3],
			digitInteresting->GetTdc(), digitInteresting->GetAdc(),
			digitInteresting->GetToT(), digitInteresting->GetToT()*digitInteresting->GetToT(),
			digitsInVolumeIndices[digIndex],
			digitInteresting->GetTrack(0), digitInteresting->GetTrack(1), digitInteresting->GetTrack(2)));

	selectedDigit[interestingCounter][digitIndexLocal] = new // AdC
	  AliTOFselectedDigit(vol, digitInteresting->GetTdc(),
			      digitInteresting->GetAdc(), digitInteresting->GetToT(),
			      digitInteresting->GetToT()*digitInteresting->GetToT(),
			      digitsInVolumeIndices[digIndex],
			      digitInteresting->GetTracks());
      }
      if (digitIndexLocal==-1) interestingCounter--; // AdC
    }


    vol[4] = AliTOFGeometry::NpadZ()-(iPadZ+1), vol[3] = iPadX;

    if (fTOFdigitMap->GetNumberOfDigits(vol)>0) {
      interestingCounter++;
      for (ii=0; ii<kMaxNumberOfDigitsPerVolume; ii++)
	digitsInVolumeIndices[ii] = -1;
      fTOFdigitMap->GetDigitIndex(vol, digitsInVolumeIndices);
      digitIndexLocal=-1; // AdC
      for(Int_t digIndex=0; digIndex<kMaxNumberOfDigitsPerVolume; digIndex++) {
	if (digitsInVolumeIndices[digIndex]<0) continue;
	digitInteresting = (AliTOFdigit*)fDigits->UncheckedAt(digitsInVolumeIndices[digIndex]);
	if (digitInteresting->GetToT()<=0) continue; // AdC
	digitIndexLocal++; // AdC

	AliDebug(1,Form(" %2d %1d %2d %1d %2d  %d %d %d %d  %5d  %5d %5d %5d",
			vol[0], vol[1], vol[2] ,vol[4], vol[3],
			digitInteresting->GetTdc(), digitInteresting->GetAdc(),
			digitInteresting->GetToT(), digitInteresting->GetToT()*digitInteresting->GetToT(),
			digitsInVolumeIndices[digIndex],
			digitInteresting->GetTrack(0), digitInteresting->GetTrack(1), digitInteresting->GetTrack(2)));

	selectedDigit[interestingCounter][digitIndexLocal] = new // AdC
	  AliTOFselectedDigit(vol, digitInteresting->GetTdc(),
			      digitInteresting->GetAdc(), digitInteresting->GetToT(),
			      digitInteresting->GetToT()*digitInteresting->GetToT(),
			      digitsInVolumeIndices[digIndex],
			      digitInteresting->GetTracks());
      }
      if (digitIndexLocal==-1) interestingCounter--; // AdC
    }


    vol[4] = AliTOFGeometry::NpadZ()-(iPadZ+1), vol[3] = iPadX+1;

    if (fTOFdigitMap->GetNumberOfDigits(vol)>0) {
      interestingCounter++;
      for (ii=0; ii<kMaxNumberOfDigitsPerVolume; ii++)
	digitsInVolumeIndices[ii] = -1;
      fTOFdigitMap->GetDigitIndex(vol, digitsInVolumeIndices);
      digitIndexLocal=-1; // AdC
      for(Int_t digIndex=0; digIndex<kMaxNumberOfDigitsPerVolume; digIndex++) {
	if (digitsInVolumeIndices[digIndex]<0) continue;
	digitInteresting = (AliTOFdigit*)fDigits->UncheckedAt(digitsInVolumeIndices[digIndex]);
	if (digitInteresting->GetToT()<=0) continue; // AdC
	digitIndexLocal++; // AdC

	AliDebug(1,Form(" %2d %1d %2d %1d %2d  %d %d %d %d  %5d  %5d %5d %5d",
			vol[0], vol[1], vol[2] ,vol[4], vol[3],
			digitInteresting->GetTdc(), digitInteresting->GetAdc(),
			digitInteresting->GetToT(), digitInteresting->GetToT()*digitInteresting->GetToT(),
			digitsInVolumeIndices[digIndex],
			digitInteresting->GetTrack(0), digitInteresting->GetTrack(1), digitInteresting->GetTrack(2)));

	selectedDigit[interestingCounter][digitIndexLocal] = new // AdC
	  AliTOFselectedDigit(vol, digitInteresting->GetTdc(),
			      digitInteresting->GetAdc(), digitInteresting->GetToT(),
			      digitInteresting->GetToT()*digitInteresting->GetToT(),
			      digitsInVolumeIndices[digIndex],
			      digitInteresting->GetTracks());
      }
      if (digitIndexLocal==-1) interestingCounter--; // AdC
    }

    AliDebug(1,Form(" e adesso %1d", interestingCounter+1));

    for (Int_t adesso1=0; adesso1<interestingCounter+1; adesso1++) {
      for (Int_t firstIndex=0; firstIndex<kMaxNumberOfDigitsPerVolume; firstIndex++) {
	if (selectedDigit[adesso1][firstIndex]==0x0) continue;

	for (Int_t adesso2=adesso1+1; adesso2<interestingCounter+1; adesso2++) {
	  for (Int_t secondIndex=0; secondIndex<kMaxNumberOfDigitsPerVolume; secondIndex++) {
	    if (selectedDigit[adesso2][secondIndex]==0x0) continue;

	    if (TMath::Abs(selectedDigit[adesso1][firstIndex]->GetTDC()-selectedDigit[adesso2][secondIndex]->GetTDC())>fMaxDeltaTime) continue;

	    interestingTOF[0] = selectedDigit[adesso1][firstIndex]->GetTDC();
	    interestingTOT[0] = selectedDigit[adesso1][firstIndex]->GetTOT();
	    interestingADC[0] = selectedDigit[adesso1][firstIndex]->GetADC();
	    interestingWeight[0] = selectedDigit[adesso1][firstIndex]->GetWeight();
	    Int_t vol1[5]; for(jj=0; jj<5; jj++) vol1[jj] = selectedDigit[adesso1][firstIndex]->GetDetectorIndex(jj);
	    AliDebug(1,Form(" %2d %1d %2d %1d %2d", vol1[0], vol1[1], vol1[2], vol1[4], vol1[3]));
	    Int_t volDum = vol1[3];
	    vol1[3] = vol1[4];
	    vol1[4] = volDum;
	    fTOFGeometry->GetPosPar(vol1,pos);
	    interestingX[0] = pos[0];
	    interestingY[0] = pos[1];
	    interestingZ[0] = pos[2];

	    interestingTOF[1] = selectedDigit[adesso2][secondIndex]->GetTDC();
	    interestingTOT[1] = selectedDigit[adesso2][secondIndex]->GetTOT();
	    interestingADC[1] = selectedDigit[adesso2][secondIndex]->GetADC();
	    interestingWeight[1] = selectedDigit[adesso2][secondIndex]->GetWeight();
	    Int_t vol2[5]; for(jj=0; jj<5; jj++) vol2[jj] = selectedDigit[adesso2][secondIndex]->GetDetectorIndex(jj);
	    AliDebug(1,Form(" %2d %1d %2d %1d %2d", vol2[0], vol2[1], vol2[2], vol2[4], vol2[3]));
	    volDum = vol2[3];
	    vol2[3] = vol2[4];
	    vol2[4] = volDum;
	    fTOFGeometry->GetPosPar(vol2,pos);
	    interestingX[1] = pos[0];
	    interestingY[1] = pos[1];
	    interestingZ[1] = pos[2];


	    AverageCalculations(2, interestingX, interestingY, interestingZ,
				interestingTOF, interestingTOT, interestingADC,
				interestingWeight,
				parTOF, posClus, check);

	    for (jj=0; jj<5; jj++) det[jj] = -1;
	    for (jj=0; jj<3; jj++) posF[jj] = posClus[jj];
	    fTOFGeometry->GetDetID(posF, det);

	    volIdClus = fTOFGeometry->GetAliSensVolIndex(det[0],det[1],det[2]);
	    //volIdClus = GetClusterVolIndex(det);

	    for (jj=3; jj<11; jj++) padsCluster[jj] = -1;
	    padsCluster[3] = selectedDigit[adesso1][firstIndex]->GetDetectorIndex(4);
	    padsCluster[4] = selectedDigit[adesso1][firstIndex]->GetDetectorIndex(3);
	    padsCluster[5] = selectedDigit[adesso2][secondIndex]->GetDetectorIndex(4);
	    padsCluster[6] = selectedDigit[adesso2][secondIndex]->GetDetectorIndex(3);

	    for (jj=0; jj<6; jj++) covClus[jj] = 0.;
	    Int_t ** indDet = new Int_t*[2];
	    for (jj=0; jj<2; jj++) indDet[jj] = new Int_t [5];
	    for (jj=0; jj<2; jj++) indDet[jj][0] = nSector;
	    for (jj=0; jj<2; jj++) indDet[jj][1] = nPlate;
	    for (jj=0; jj<2; jj++) indDet[jj][2] = nStrip;
	    for (jj=0; jj<2; jj++) indDet[jj][3] = padsCluster[2*jj+3];
	    for (jj=0; jj<2; jj++) indDet[jj][4] = padsCluster[2*jj+1+3];
	    GetClusterPars(/*check,*/ 2, indDet, interestingWeight, posClus, covClus);
	    for (jj=0; jj<2; jj++) delete [] indDet[jj];
	    delete [] indDet;

	    // To fill the track index array
	    dummyCounter=-1;
	    for (jj=0; jj<kMaxNumberOfTracksPerDigit; jj++) tracks[jj] = -1;
	    for (Int_t kk=0; kk<kMaxNumberOfTracksPerDigit; kk++) { // three is the max number of tracks associated to one digit
	      if (selectedDigit[adesso1][firstIndex]->GetTrackLabel(kk)==-1) continue;
	      else {
		dummyCounter++;
		tracks[dummyCounter] = selectedDigit[adesso1][firstIndex]->GetTrackLabel(kk);
	      }
	    }
	    for (Int_t kk=0; kk<kMaxNumberOfTracksPerDigit; kk++) { // three is the max number of tracks associated to one digit
	      if (selectedDigit[adesso2][secondIndex]->GetTrackLabel(kk)==-1) continue;
	      else {

		alreadyStored = kFALSE;
		for (jj=0; jj<dummyCounter+1; jj++)
		  alreadyStored = alreadyStored || (tracks[jj]!=-1 && tracks[jj]==selectedDigit[adesso2][secondIndex]->GetTrackLabel(kk));

		if (alreadyStored) continue;
		if (dummyCounter==2) { // three is the max number of tracks associated to one cluster
		  AliWarning("  Siamo al limite!");
		  continue;
		}

		dummyCounter++;
		tracks[dummyCounter] = selectedDigit[adesso2][secondIndex]->GetTrackLabel(kk);

	      }

	    }


	    AliTOFcluster *tofCluster =
	      new AliTOFcluster(volIdClus, posClus[0], posClus[1], posClus[2],
				covClus[0], covClus[1], covClus[2], covClus[3], covClus[4], covClus[5],
				tracks, det, parTOF, status, selectedDigit[adesso1][firstIndex]->GetIndex()); // to be updated
	    InsertCluster(tofCluster);

	    AliDebug(2, Form("       %4d  %f %f %f  %f %f %f %f %f %f  %3d %3d %3d  %2d %1d %2d %1d %2d  %4d %3d %3d %4d %4d  %1d  %4d", 
			     volIdClus, posClus[0], posClus[1], posClus[2],
			     covClus[0], covClus[1], covClus[2], covClus[3], covClus[4], covClus[5],
			     tracks[0], tracks[1], tracks[2],
			     det[0], det[1], det[2], det[3], det[4],
			     parTOF[0], parTOF[1], parTOF[2], parTOF[3], parTOF[4],
			     status, selectedDigit[adesso1][firstIndex]->GetIndex()));

	    volDum = vol1[3];
	    vol1[3] = vol1[4];
	    vol1[4] = volDum;
	    fTOFdigitMap->ResetDigitNumber(vol1,selectedDigit[adesso1][firstIndex]->GetIndex());
	    volDum = vol2[3];
	    vol2[3] = vol2[4];
	    vol2[4] = volDum;
	    fTOFdigitMap->ResetDigitNumber(vol2,selectedDigit[adesso2][secondIndex]->GetIndex());


	  } // close loop on second digit
	} // close loop on adesso2

      } // close loop on first digit
    } // close loop on adesso1

    for (ii=0; ii<kMaxNumberOfInterestingPads; ii++)
      for (jj=0; jj<kMaxNumberOfDigitsPerVolume; jj++)
	selectedDigit[ii][jj] = 0x0;

  } // loop on iPad

  for (ii=0; ii<kMaxNumberOfInterestingPads; ii++) {
    for (jj=0; jj<kMaxNumberOfDigitsPerVolume; jj++) {
      delete [] selectedDigit[ii][jj];
      selectedDigit[ii][jj] = 0x0;
    }
    delete [] selectedDigit[ii];
    selectedDigit[ii] = 0x0;
  }
  delete [] selectedDigit;
  selectedDigit = 0x0;

}
//_____________________________________________________________________________

void AliTOFClusterFinderV1::FindClustersPerStrip(Int_t nSector,
						 Int_t nPlate,
						 Int_t nStrip,
						 Int_t group)
{
  //
  // This function searches the neighbouring digits (stored in the fDigits object),
  // to perform clusters (stored in the fTofClusters array).
  //
  // Each strip is read four times:
  //  - 1st time: it searches possible clusters formed by four
  //              neighbouring digits;
  //  - 2nd time: it searches possible clusters formed by three
  //              neighbouring digits;
  //  - 3rd time: it searches possible clusters formed by two
  //              neighbouring digits;
  //  - 4th time: the remaining isolated digits have been transformed
  //              in clusters.
  // This research has been made by checking the fTOFdigitMap object,
  // filled at digits/raw-data reading time.
  //

  const Int_t kMaxNumberOfInterestingPads = 4;
  const Int_t kMaxNumberOfTracksPerDigit = 3;
  const Int_t kMaxNumberOfDigitsPerVolume = 10;

  Int_t ii = 0;

  Int_t digitsInVolumeIndices[kMaxNumberOfDigitsPerVolume];
  for (ii=0; ii<kMaxNumberOfDigitsPerVolume; ii++)
    digitsInVolumeIndices[ii] = -1;

  Int_t vol[5] = {nSector,nPlate,nStrip,-1,-1};

  Float_t pos[3] = {0.,0.,0.};

  Int_t jj = 0;
  Int_t interestingPadX[kMaxNumberOfInterestingPads];
  for (jj=0; jj<kMaxNumberOfInterestingPads; jj++) interestingPadX[jj] = -1;
  Int_t interestingPadZ[kMaxNumberOfInterestingPads];
  for (jj=0; jj<kMaxNumberOfInterestingPads; jj++) interestingPadZ[jj] = -1;
  Double_t interestingTOT[kMaxNumberOfInterestingPads];
  for (jj=0; jj<kMaxNumberOfInterestingPads; jj++) interestingTOT[jj] = 0;
  Double_t interestingADC[kMaxNumberOfInterestingPads];
  for (jj=0; jj<kMaxNumberOfInterestingPads; jj++) interestingADC[jj] = 0;
  Double_t interestingTOF[kMaxNumberOfInterestingPads];
  for (jj=0; jj<kMaxNumberOfInterestingPads; jj++) interestingTOF[jj] = 0;
  Double_t interestingWeight[kMaxNumberOfInterestingPads];
  for (jj=0; jj<kMaxNumberOfInterestingPads; jj++) interestingWeight[jj] = 0;

  Float_t interestingX[kMaxNumberOfInterestingPads];
  for (jj=0; jj<kMaxNumberOfInterestingPads; jj++) interestingX[jj] = 0;
  Float_t interestingY[kMaxNumberOfInterestingPads];
  for (jj=0; jj<kMaxNumberOfInterestingPads; jj++) interestingY[jj] = 0;
  Float_t interestingZ[kMaxNumberOfInterestingPads];
  for (jj=0; jj<kMaxNumberOfInterestingPads; jj++) interestingZ[jj] = 0;

  Float_t interDigit[kMaxNumberOfInterestingPads];
  for (jj=0; jj<kMaxNumberOfInterestingPads; jj++) interDigit[jj] = 0;

  Int_t padsCluster[11];
  padsCluster[0] = nSector;
  padsCluster[1] = nPlate;
  padsCluster[2] = nStrip;
  for (jj=3; jj<11; jj++) padsCluster[jj] = -1;

  Int_t interestingCounter=-1;
  Int_t digitIndexLocal = -1;
  Int_t iPad  = -1;
  Int_t iPadX = -1;
  Int_t iPadZ = -1;

  Bool_t check = kFALSE;
  Int_t parTOF[7];
  for (jj=0; jj<7; jj++) parTOF[jj] = 0;
  Double_t posClus[3];
  for (jj=0; jj<3; jj++) posClus[jj] = 0.;
  Int_t det[5];
  for (jj=0; jj<5; jj++) det[jj] = -1;
  Float_t posF[3];
  for (jj=0; jj<3; jj++) posF[jj] = 0.;
  UShort_t volIdClus = 0;
  Bool_t status = kTRUE; //assume all sim channels ok in the beginning...
  Double_t covClus[6];
  for (jj=0; jj<6; jj++) covClus[jj] = 0.;
  Int_t tracks[kMaxNumberOfTracksPerDigit];
  for (jj=0; jj<kMaxNumberOfTracksPerDigit; jj++) tracks[jj] = -1;
  Int_t dummyCounter=-1;
  Bool_t alreadyStored = kFALSE;

  AliTOFselectedDigit ***selectedDigit = new AliTOFselectedDigit**[kMaxNumberOfInterestingPads];
  for (ii=0; ii<kMaxNumberOfInterestingPads; ii++)
    selectedDigit[ii] = new AliTOFselectedDigit*[kMaxNumberOfDigitsPerVolume];

  for (ii=0; ii<kMaxNumberOfInterestingPads; ii++)
    for (jj=0; jj<kMaxNumberOfDigitsPerVolume; jj++) selectedDigit[ii][jj] = 0x0;

  AliTOFdigit *digitInteresting;

  group = group-1;

  for (iPad=0; iPad<AliTOFGeometry::NpadZ()*AliTOFGeometry::NpadX()-3; iPad+=2) {

    iPadZ = iPad%AliTOFGeometry::NpadZ(); //iPad%2;
    iPadX = iPad/AliTOFGeometry::NpadZ(); //iPad/2;

    AliDebug(3, Form("%2d %1d %2d\n", iPad, iPadZ, iPadX));






    interestingCounter=-1;

    vol[4] = iPadZ  , vol[3]  = iPadX;
    if (fTOFdigitMap->GetNumberOfDigits(vol)>0)
      interestingCounter++;

    vol[4] = iPadZ, vol[3] = iPadX+1;
    if (fTOFdigitMap->GetNumberOfDigits(vol)>0)
      interestingCounter++;

    vol[4] = AliTOFGeometry::NpadZ()-(iPadZ+1), vol[3] = iPadX;
    if (fTOFdigitMap->GetNumberOfDigits(vol)>0)
      interestingCounter++;

    vol[4] = AliTOFGeometry::NpadZ()-(iPadZ+1), vol[3] = iPadX+1;
    if (fTOFdigitMap->GetNumberOfDigits(vol)>0)
      interestingCounter++;

    if (interestingCounter!=group) continue;
    else interestingCounter=-1;


    for (ii=0; ii<kMaxNumberOfInterestingPads; ii++)
      for (jj=0; jj<kMaxNumberOfDigitsPerVolume; jj++)
	selectedDigit[ii][jj] = 0x0;


    vol[4] = iPadZ, vol[3] = iPadX;

    if (fTOFdigitMap->GetNumberOfDigits(vol)>0) {
      interestingCounter++;
      for (ii=0; ii<kMaxNumberOfDigitsPerVolume; ii++)
	digitsInVolumeIndices[ii] = -1;
      fTOFdigitMap->GetDigitIndex(vol, digitsInVolumeIndices);
      digitIndexLocal=-1; // AdC
      for(Int_t digIndex=0; digIndex<kMaxNumberOfDigitsPerVolume; digIndex++) {
	if (digitsInVolumeIndices[digIndex]<0) continue;
	digitInteresting = (AliTOFdigit*)fDigits->UncheckedAt(digitsInVolumeIndices[digIndex]);
	if (digitInteresting->GetToT()<=0) continue; // AdC
	digitIndexLocal++; // AdC

	AliDebug(1,Form(" %2d %1d %2d %1d %2d  %d %d %d %d  %5d  %5d %5d %5d",
			vol[0], vol[1], vol[2] ,vol[4], vol[3],
			digitInteresting->GetTdc(), digitInteresting->GetAdc(),
			digitInteresting->GetToT(), digitInteresting->GetToT()*digitInteresting->GetToT(),
			digitsInVolumeIndices[digIndex],
			digitInteresting->GetTrack(0), digitInteresting->GetTrack(1), digitInteresting->GetTrack(2)));
	AliDebug(1,Form("   to try   %d %d ",digIndex, interestingCounter));

	selectedDigit[interestingCounter][digitIndexLocal] = new // AdC
	  AliTOFselectedDigit(vol, digitInteresting->GetTdc(),
			      digitInteresting->GetAdc(), digitInteresting->GetToT(),
			      digitInteresting->GetToT()*digitInteresting->GetToT(),
			      digitsInVolumeIndices[digIndex],
			      digitInteresting->GetTracks());
      }
      if (digitIndexLocal==-1) interestingCounter--; // AdC
    }


    vol[4] = iPadZ, vol[3] = iPadX+1;

    if (fTOFdigitMap->GetNumberOfDigits(vol)>0) {
      interestingCounter++;
      for (ii=0; ii<kMaxNumberOfDigitsPerVolume; ii++)
	digitsInVolumeIndices[ii] = -1;
      fTOFdigitMap->GetDigitIndex(vol, digitsInVolumeIndices);
      digitIndexLocal=-1; // AdC
      for(Int_t digIndex=0; digIndex<kMaxNumberOfDigitsPerVolume; digIndex++) {
	if (digitsInVolumeIndices[digIndex]<0) continue;
	digitInteresting = (AliTOFdigit*)fDigits->UncheckedAt(digitsInVolumeIndices[digIndex]);
	if (digitInteresting->GetToT()<=0) continue; // AdC
	digitIndexLocal++; // AdC

	AliDebug(1,Form(" %2d %1d %2d %1d %2d  %d %d %d %d  %5d  %5d %5d %5d",
			vol[0], vol[1], vol[2] ,vol[4], vol[3],
			digitInteresting->GetTdc(), digitInteresting->GetAdc(),
			digitInteresting->GetToT(), digitInteresting->GetToT()*digitInteresting->GetToT(),
			digitsInVolumeIndices[digIndex],
			digitInteresting->GetTrack(0), digitInteresting->GetTrack(1), digitInteresting->GetTrack(2)));
	AliDebug(1,Form("   to try   %d %d ",digIndex, interestingCounter));

	selectedDigit[interestingCounter][digitIndexLocal] = new // AdC
	  AliTOFselectedDigit(vol, digitInteresting->GetTdc(),
			      digitInteresting->GetAdc(), digitInteresting->GetToT(),
			      digitInteresting->GetToT()*digitInteresting->GetToT(),
			      digitsInVolumeIndices[digIndex],
			      digitInteresting->GetTracks());
      }
      if (digitIndexLocal==-1) interestingCounter--; // AdC
    }


    vol[4] = AliTOFGeometry::NpadZ()-(iPadZ+1), vol[3] = iPadX;

    if (fTOFdigitMap->GetNumberOfDigits(vol)>0) {
      interestingCounter++;
      for (ii=0; ii<kMaxNumberOfDigitsPerVolume; ii++)
	digitsInVolumeIndices[ii] = -1;
      fTOFdigitMap->GetDigitIndex(vol, digitsInVolumeIndices);
      digitIndexLocal=-1; // AdC
      for(Int_t digIndex=0; digIndex<kMaxNumberOfDigitsPerVolume; digIndex++) {
	if (digitsInVolumeIndices[digIndex]<0) continue;
	digitInteresting = (AliTOFdigit*)fDigits->UncheckedAt(digitsInVolumeIndices[digIndex]);
	if (digitInteresting->GetToT()<=0) continue; // AdC
	digitIndexLocal++; // AdC

	AliDebug(1,Form(" %2d %1d %2d %1d %2d  %d %d %d %d  %5d  %5d %5d %5d",
			vol[0], vol[1], vol[2] ,vol[4], vol[3],
			digitInteresting->GetTdc(), digitInteresting->GetAdc(),
			digitInteresting->GetToT(), digitInteresting->GetToT()*digitInteresting->GetToT(),
			digitsInVolumeIndices[digIndex],
			digitInteresting->GetTrack(0), digitInteresting->GetTrack(1), digitInteresting->GetTrack(2)));
	AliDebug(1,Form("   to try   %d %d ",digIndex, interestingCounter));

	selectedDigit[interestingCounter][digitIndexLocal] = new // AdC
	  AliTOFselectedDigit(vol, digitInteresting->GetTdc(),
			      digitInteresting->GetAdc(), digitInteresting->GetToT(),
			      digitInteresting->GetToT()*digitInteresting->GetToT(),
			      digitsInVolumeIndices[digIndex],
			      digitInteresting->GetTracks());
      }
      if (digitIndexLocal==-1) interestingCounter--; // AdC
    }


    vol[4] = AliTOFGeometry::NpadZ()-(iPadZ+1), vol[3] = iPadX+1;

    if (fTOFdigitMap->GetNumberOfDigits(vol)>0) {
      interestingCounter++;
      for (ii=0; ii<kMaxNumberOfDigitsPerVolume; ii++)
	digitsInVolumeIndices[ii] = -1;
      fTOFdigitMap->GetDigitIndex(vol, digitsInVolumeIndices);
      digitIndexLocal=-1; // AdC
      for(Int_t digIndex=0; digIndex<kMaxNumberOfDigitsPerVolume; digIndex++) {
	if (digitsInVolumeIndices[digIndex]<0) continue;
	digitInteresting = (AliTOFdigit*)fDigits->UncheckedAt(digitsInVolumeIndices[digIndex]);
	if (digitInteresting->GetToT()<=0) continue; // AdC
	digitIndexLocal++; // AdC

	AliDebug(1,Form(" %2d %1d %2d %1d %2d  %d %d %d %d  %5d  %5d %5d %5d",
			vol[0], vol[1], vol[2] ,vol[4], vol[3],
			digitInteresting->GetTdc(), digitInteresting->GetAdc(),
			digitInteresting->GetToT(), digitInteresting->GetToT()*digitInteresting->GetToT(),
			digitsInVolumeIndices[digIndex],
			digitInteresting->GetTrack(0), digitInteresting->GetTrack(1), digitInteresting->GetTrack(2)));
	AliDebug(1,Form("   to try   %d %d ",digIndex, interestingCounter));

	selectedDigit[interestingCounter][digitIndexLocal] = new // AdC
	  AliTOFselectedDigit(vol, digitInteresting->GetTdc(),
			      digitInteresting->GetAdc(), digitInteresting->GetToT(),
			      digitInteresting->GetToT()*digitInteresting->GetToT(),
			      digitsInVolumeIndices[digIndex],
			      digitInteresting->GetTracks());
      }
      if (digitIndexLocal==-1) interestingCounter--; // AdC
    }

    AliDebug(1,Form(" e adesso %1d", interestingCounter+1));

    Int_t adesso1 = -1;
    Int_t adesso2 = -1;
    Int_t adesso3 = -1;
    Int_t adesso4 = -1;

    switch(interestingCounter+1) {

    case 2:

      //for (adesso1=0; adesso1<interestingCounter+1; adesso1++) {
      adesso1 = 0;
	for (Int_t firstIndex=0; firstIndex<kMaxNumberOfDigitsPerVolume; firstIndex++) {
	  if (selectedDigit[adesso1][firstIndex]==0x0) continue;

	  //for (adesso2=adesso1+1; adesso2<interestingCounter+1; adesso2++) {
	  adesso2 = 1;
	    for (Int_t secondIndex=0; secondIndex<kMaxNumberOfDigitsPerVolume; secondIndex++) {
	      if (selectedDigit[adesso2][secondIndex]==0x0) continue;

	      if (TMath::Abs(selectedDigit[adesso1][firstIndex]->GetTDC()-selectedDigit[adesso2][secondIndex]->GetTDC())>fMaxDeltaTime) {
		AliDebug(1,Form(" selD1[%d][%d]->GetTDC()=%d selD2[%d][%d]->GetTDC()=%d -- %d ",
				adesso1,firstIndex,(Int_t)selectedDigit[adesso1][firstIndex]->GetTDC(),
				adesso2,secondIndex,(Int_t)selectedDigit[adesso2][secondIndex]->GetTDC(),
				fMaxDeltaTime));
		continue;
	      }

	      AliDebug(1, Form(" %1d %1d (0x%p) %1d %1d (0x%p)", adesso1, firstIndex,selectedDigit[adesso1][firstIndex],
			       adesso2, secondIndex,selectedDigit[adesso2][secondIndex]));

	      interestingTOF[0] = selectedDigit[adesso1][firstIndex]->GetTDC();
	      interestingTOT[0] = selectedDigit[adesso1][firstIndex]->GetTOT();
	      interestingADC[0] = selectedDigit[adesso1][firstIndex]->GetADC();
	      interestingWeight[0] = selectedDigit[adesso1][firstIndex]->GetWeight();
	      Int_t vol1[5]; for(jj=0; jj<5; jj++) vol1[jj] = selectedDigit[adesso1][firstIndex]->GetDetectorIndex(jj);
	      AliDebug(1,Form(" %2d %1d %2d %1d %2d", vol1[0], vol1[1], vol1[2], vol1[4], vol1[3]));
	      Int_t volDum = vol1[3];
	      vol1[3] = vol1[4];
	      vol1[4] = volDum;
	      fTOFGeometry->GetPosPar(vol1,pos);
	      AliDebug(1,Form(" %f %f %f", pos[0], pos[1], pos[2]));
	      interestingX[0] = pos[0];
	      interestingY[0] = pos[1];
	      interestingZ[0] = pos[2];

	      interestingTOF[1] = selectedDigit[adesso2][secondIndex]->GetTDC();
	      interestingTOT[1] = selectedDigit[adesso2][secondIndex]->GetTOT();
	      interestingADC[1] = selectedDigit[adesso2][secondIndex]->GetADC();
	      interestingWeight[1] = selectedDigit[adesso2][secondIndex]->GetWeight();
	      Int_t vol2[5]; for(jj=0; jj<5; jj++) vol2[jj] = selectedDigit[adesso2][secondIndex]->GetDetectorIndex(jj);
	      AliDebug(1,Form(" %2d %1d %2d %1d %2d", vol2[0], vol2[1], vol2[2], vol2[4], vol2[3]));
	      volDum = vol2[3];
	      vol2[3] = vol2[4];
	      vol2[4] = volDum;
	      fTOFGeometry->GetPosPar(vol2,pos);
	      AliDebug(1,Form(" %f %f %f", pos[0], pos[1], pos[2]));
	      interestingX[1] = pos[0];
	      interestingY[1] = pos[1];
	      interestingZ[1] = pos[2];


	      AverageCalculations(interestingCounter+1,
				  interestingX, interestingY, interestingZ,
				  interestingTOF, interestingTOT, interestingADC,
				  interestingWeight,
				  parTOF, posClus, check);

	      for (jj=0; jj<5; jj++) det[jj] = -1;
	      for (jj=0; jj<3; jj++) posF[jj] = posClus[jj];

	      AliDebug(1,Form(" %f %f %f", posF[0], posF[1], posF[2]));
	      fTOFGeometry->GetDetID(posF, det);
	      AliDebug(1,Form(" %2d %1d %2d %1d %2d", det[0], det[1], det[2], det[3], det[4]));

	      volIdClus = fTOFGeometry->GetAliSensVolIndex(det[0],det[1],det[2]);
	      //volIdClus = GetClusterVolIndex(det);

	      for (jj=3; jj<11; jj++) padsCluster[jj] = -1;
	      padsCluster[3] = selectedDigit[adesso1][firstIndex]->GetDetectorIndex(4);
	      padsCluster[4] = selectedDigit[adesso1][firstIndex]->GetDetectorIndex(3);
	      padsCluster[5] = selectedDigit[adesso2][secondIndex]->GetDetectorIndex(4);
	      padsCluster[6] = selectedDigit[adesso2][secondIndex]->GetDetectorIndex(3);

	      for (jj=0; jj<6; jj++) covClus[jj] = 0.;
	      Int_t ** indDet = new Int_t*[interestingCounter+1];
	      for (jj=0; jj<interestingCounter+1; jj++) indDet[jj] = new Int_t [5];
	      for (jj=0; jj<interestingCounter+1; jj++) indDet[jj][0] = nSector;
	      for (jj=0; jj<interestingCounter+1; jj++) indDet[jj][1] = nPlate;
	      for (jj=0; jj<interestingCounter+1; jj++) indDet[jj][2] = nStrip;
	      for (jj=0; jj<interestingCounter+1; jj++) indDet[jj][3] = padsCluster[2*jj+3];
	      for (jj=0; jj<interestingCounter+1; jj++) indDet[jj][4] = padsCluster[2*jj+1+3];
	      GetClusterPars(/*check,*/ interestingCounter+1, indDet, interestingWeight, posClus, covClus);
	      for (jj=0; jj<interestingCounter+1; jj++) delete [] indDet[jj];
	      delete [] indDet;

	      // To fill the track index array
	      dummyCounter=-1;
	      for (jj=0; jj<kMaxNumberOfTracksPerDigit; jj++) tracks[jj] = -1;
	      for (Int_t kk=0; kk<kMaxNumberOfTracksPerDigit; kk++) { // three is the max number of tracks associated to one digit
		if (selectedDigit[adesso1][firstIndex]->GetTrackLabel(kk)==-1) continue;
		else {
		  dummyCounter++;
		  tracks[dummyCounter] = selectedDigit[adesso1][firstIndex]->GetTrackLabel(kk);
		}
	      }
	      for (Int_t kk=0; kk<kMaxNumberOfTracksPerDigit; kk++) { // three is the max number of tracks associated to one digit
		if (selectedDigit[adesso2][secondIndex]->GetTrackLabel(kk)==-1) continue;
		else {

		  alreadyStored = kFALSE;
		  for (jj=0; jj<dummyCounter+1; jj++)
		    alreadyStored = alreadyStored || (tracks[jj]!=-1 && tracks[jj]==selectedDigit[adesso2][secondIndex]->GetTrackLabel(kk));

		  if (alreadyStored) continue;
		  if (dummyCounter==2) { // three is the max number of tracks associated to one cluster
		    AliWarning("  Siamo al limite!");
		    continue;
		  }

		  dummyCounter++;
		  tracks[dummyCounter] = selectedDigit[adesso2][secondIndex]->GetTrackLabel(kk);

		}

	      }


	      AliTOFcluster *tofCluster =
		new AliTOFcluster(volIdClus, posClus[0], posClus[1], posClus[2],
				  covClus[0], covClus[1], covClus[2], covClus[3], covClus[4], covClus[5],
				  tracks, det, parTOF, status, selectedDigit[adesso1][firstIndex]->GetIndex()); //to updated
	      InsertCluster(tofCluster);

	      AliDebug(2, Form("       %4d  %f %f %f  %f %f %f %f %f %f  %3d %3d %3d  %2d %1d %2d %1d %2d  %4d %3d %3d %4d %4d  %1d  %4d", 
			       volIdClus, posClus[0], posClus[1], posClus[2],
			       covClus[0], covClus[1], covClus[2], covClus[3], covClus[4], covClus[5],
			       tracks[0], tracks[1], tracks[2],
			       det[0], det[1], det[2], det[3], det[4],
			       parTOF[0], parTOF[1], parTOF[2], parTOF[3], parTOF[4],
			       status, selectedDigit[adesso1][firstIndex]->GetIndex()));

	      volDum = vol1[3];
	      vol1[3] = vol1[4];
	      vol1[4] = volDum;
	      fTOFdigitMap->ResetDigitNumber(vol1,selectedDigit[adesso1][firstIndex]->GetIndex());
	      volDum = vol2[3];
	      vol2[3] = vol2[4];
	      vol2[4] = volDum;
	      fTOFdigitMap->ResetDigitNumber(vol2,selectedDigit[adesso2][secondIndex]->GetIndex());


	    } // close loop on second digit
	    //} // close loop on adesso2

	} // close loop on first digit
	  //} // close loop on adesso1


      break;

    case 3:

      //for (adesso1=0; adesso1<interestingCounter+1; adesso1++) {
      adesso1 = 0;
	for (Int_t firstIndex=0; firstIndex<kMaxNumberOfDigitsPerVolume; firstIndex++) {
	  if (selectedDigit[adesso1][firstIndex]==0x0) continue;

	  //for (adesso2=adesso1+1; adesso2<interestingCounter+1; adesso2++) {
	  adesso2 = 1;
	    for (Int_t secondIndex=0; secondIndex<kMaxNumberOfDigitsPerVolume; secondIndex++) {
	      if (selectedDigit[adesso2][secondIndex]==0x0) continue;

	      //for (adesso3=adesso2+1; adesso3<interestingCounter+1; adesso3++) {
	      adesso3 = 2;
		for (Int_t thirdIndex=0; thirdIndex<kMaxNumberOfDigitsPerVolume; thirdIndex++) {
		  if (selectedDigit[adesso3][thirdIndex]==0x0) continue;


		  if (TMath::Abs(selectedDigit[adesso1][firstIndex]->GetTDC()-selectedDigit[adesso2][secondIndex]->GetTDC())>fMaxDeltaTime
		      ||
		      TMath::Abs(selectedDigit[adesso1][firstIndex]->GetTDC()-selectedDigit[adesso3][thirdIndex]->GetTDC())>fMaxDeltaTime
		      ||
		      TMath::Abs(selectedDigit[adesso2][secondIndex]->GetTDC()-selectedDigit[adesso3][thirdIndex]->GetTDC())>fMaxDeltaTime) continue;

		  interestingTOF[0] = selectedDigit[adesso1][firstIndex]->GetTDC();
		  interestingTOT[0] = selectedDigit[adesso1][firstIndex]->GetTOT();
		  interestingADC[0] = selectedDigit[adesso1][firstIndex]->GetADC();
		  interestingWeight[0] = selectedDigit[adesso1][firstIndex]->GetWeight();
		  Int_t vol1[5]; for(jj=0; jj<5; jj++) vol1[jj] = selectedDigit[adesso1][firstIndex]->GetDetectorIndex(jj);
		  AliDebug(1,Form(" %2d %1d %2d %1d %2d", vol1[0], vol1[1], vol1[2], vol1[4], vol1[3]));
		  Int_t volDum = vol1[3];
		  vol1[3] = vol1[4];
		  vol1[4] = volDum;
		  fTOFGeometry->GetPosPar(vol1,pos);
		  interestingX[0] = pos[0];
		  interestingY[0] = pos[1];
		  interestingZ[0] = pos[2];

		  interestingTOF[1] = selectedDigit[adesso2][secondIndex]->GetTDC();
		  interestingTOT[1] = selectedDigit[adesso2][secondIndex]->GetTOT();
		  interestingADC[1] = selectedDigit[adesso2][secondIndex]->GetADC();
		  interestingWeight[1] = selectedDigit[adesso2][secondIndex]->GetWeight();
		  Int_t vol2[5]; for(jj=0; jj<5; jj++) vol2[jj] = selectedDigit[adesso2][secondIndex]->GetDetectorIndex(jj);
		  AliDebug(1,Form(" %2d %1d %2d %1d %2d", vol2[0], vol2[1], vol2[2], vol2[4], vol2[3]));
		  volDum = vol2[3];
		  vol2[3] = vol2[4];
		  vol2[4] = volDum;
		  fTOFGeometry->GetPosPar(vol2,pos);
		  interestingX[1] = pos[0];
		  interestingY[1] = pos[1];
		  interestingZ[1] = pos[2];

		  interestingTOF[2] = selectedDigit[adesso3][thirdIndex]->GetTDC();
		  interestingTOT[2] = selectedDigit[adesso3][thirdIndex]->GetTOT();
		  interestingADC[2] = selectedDigit[adesso3][thirdIndex]->GetADC();
		  interestingWeight[2] = selectedDigit[adesso3][thirdIndex]->GetWeight();
		  Int_t vol3[5]; for(jj=0; jj<5; jj++) vol3[jj] = selectedDigit[adesso3][thirdIndex]->GetDetectorIndex(jj);
		  AliDebug(1,Form(" %2d %1d %2d %1d %2d", vol3[0], vol3[1], vol3[2], vol3[4], vol3[3]));
		  volDum = vol3[3];
		  vol3[3] = vol3[4];
		  vol3[4] = volDum;
		  fTOFGeometry->GetPosPar(vol3,pos);
		  interestingX[2] = pos[0];
		  interestingY[2] = pos[1];
		  interestingZ[2] = pos[2];


		  AverageCalculations(interestingCounter+1,
				      interestingX, interestingY, interestingZ,
				      interestingTOF, interestingTOT, interestingADC,
				      interestingWeight,
				      parTOF, posClus, check);

		  for (jj=0; jj<5; jj++) det[jj] = -1;
		  for (jj=0; jj<3; jj++) posF[jj] = posClus[jj];
		  fTOFGeometry->GetDetID(posF, det);

		  volIdClus = fTOFGeometry->GetAliSensVolIndex(det[0],det[1],det[2]);
		  //volIdClus = GetClusterVolIndex(det);

		  for (jj=3; jj<11; jj++) padsCluster[jj] = -1;
		  padsCluster[3] = selectedDigit[adesso1][firstIndex]->GetDetectorIndex(4);
		  padsCluster[4] = selectedDigit[adesso1][firstIndex]->GetDetectorIndex(3);
		  padsCluster[5] = selectedDigit[adesso2][secondIndex]->GetDetectorIndex(4);
		  padsCluster[6] = selectedDigit[adesso2][secondIndex]->GetDetectorIndex(3);
		  padsCluster[7] = selectedDigit[adesso3][thirdIndex]->GetDetectorIndex(4);
		  padsCluster[8] = selectedDigit[adesso3][thirdIndex]->GetDetectorIndex(3);

		  for (jj=0; jj<6; jj++) covClus[jj] = 0.;
		  Int_t ** indDet = new Int_t*[interestingCounter+1];
		  for (jj=0; jj<interestingCounter+1; jj++) indDet[jj] = new Int_t [5];
		  for (jj=0; jj<interestingCounter+1; jj++) indDet[jj][0] = nSector;
		  for (jj=0; jj<interestingCounter+1; jj++) indDet[jj][1] = nPlate;
		  for (jj=0; jj<interestingCounter+1; jj++) indDet[jj][2] = nStrip;
		  for (jj=0; jj<interestingCounter+1; jj++) indDet[jj][3] = padsCluster[2*jj+3];
		  for (jj=0; jj<interestingCounter+1; jj++) indDet[jj][4] = padsCluster[2*jj+1+3];
		  GetClusterPars(/*check,*/ interestingCounter+1, indDet, interestingWeight, posClus, covClus);
		  for (jj=0; jj<interestingCounter+1; jj++) delete [] indDet[jj];
		  delete [] indDet;


		  // To fill the track index array
		  dummyCounter=-1;
		  for (jj=0; jj<kMaxNumberOfTracksPerDigit; jj++) tracks[jj] = -1;
		  for (Int_t kk=0; kk<kMaxNumberOfTracksPerDigit; kk++) { // three is the max number of tracks associated to one digit
		    if (selectedDigit[adesso1][firstIndex]->GetTrackLabel(kk)==-1) continue;
		    else {
		      dummyCounter++;
		      tracks[dummyCounter] = selectedDigit[adesso1][firstIndex]->GetTrackLabel(kk);
		    }
		  }
		  for (Int_t kk=0; kk<kMaxNumberOfTracksPerDigit; kk++) { // three is the max number of tracks associated to one digit
		    if (selectedDigit[adesso2][secondIndex]->GetTrackLabel(kk)==-1) continue;
		    else {

		      alreadyStored = kFALSE;
		      for (jj=0; jj<dummyCounter+1; jj++)
			alreadyStored = alreadyStored || (tracks[jj]!=-1 && tracks[jj]==selectedDigit[adesso2][secondIndex]->GetTrackLabel(kk));

		      if (alreadyStored) continue;
		      if (dummyCounter==2) { // three is the max number of tracks associated to one cluster
			AliWarning("  Siamo al limite!");
			continue;
		      }

		      dummyCounter++;
		      tracks[dummyCounter] = selectedDigit[adesso2][secondIndex]->GetTrackLabel(kk);

		    }

		  }
		  for (Int_t kk=0; kk<kMaxNumberOfTracksPerDigit; kk++) { // three is the max number of tracks associated to one digit
		    if (selectedDigit[adesso3][thirdIndex]->GetTrackLabel(kk)==-1) continue;
		    else {

		      alreadyStored = kFALSE;
		      for (jj=0; jj<dummyCounter+1; jj++)
			alreadyStored = alreadyStored || (tracks[jj]!=-1 && tracks[jj]==selectedDigit[adesso3][thirdIndex]->GetTrackLabel(kk));

		      if (alreadyStored) continue;
		      if (dummyCounter==2) { // three is the max number of tracks associated to one cluster
			AliWarning("  Siamo al limite!");
			continue;
		      }

		      dummyCounter++;
		      tracks[dummyCounter] = selectedDigit[adesso3][thirdIndex]->GetTrackLabel(kk);

		    }

		  }


		  AliTOFcluster *tofCluster =
		    new AliTOFcluster(volIdClus, posClus[0], posClus[1], posClus[2],
				      covClus[0], covClus[1], covClus[2], covClus[3], covClus[4], covClus[5],
				      tracks, det, parTOF, status, selectedDigit[adesso1][firstIndex]->GetIndex()); // to be updated
		  InsertCluster(tofCluster);

		  AliDebug(2, Form("       %4d  %f %f %f  %f %f %f %f %f %f  %3d %3d %3d  %2d %1d %2d %1d %2d  %4d %3d %3d %4d %4d  %1d  %4d", 
				   volIdClus, posClus[0], posClus[1], posClus[2],
				   covClus[0], covClus[1], covClus[2], covClus[3], covClus[4], covClus[5],
				   tracks[0], tracks[1], tracks[2],
				   det[0], det[1], det[2], det[3], det[4],
				   parTOF[0], parTOF[1], parTOF[2], parTOF[3], parTOF[4],
				   status, selectedDigit[adesso1][firstIndex]->GetIndex()));

		  volDum = vol1[3];
		  vol1[3] = vol1[4];
		  vol1[4] = volDum;
		  fTOFdigitMap->ResetDigitNumber(vol1,selectedDigit[adesso1][firstIndex]->GetIndex());
		  volDum = vol2[3];
		  vol2[3] = vol2[4];
		  vol2[4] = volDum;
		  fTOFdigitMap->ResetDigitNumber(vol2,selectedDigit[adesso2][secondIndex]->GetIndex());
		  volDum = vol3[3];
		  vol3[3] = vol3[4];
		  vol3[4] = volDum;
		  fTOFdigitMap->ResetDigitNumber(vol3,selectedDigit[adesso3][thirdIndex]->GetIndex());


		} // close loop on third digit
		//} // close loop on adesso3

	    } // close loop on second digit
	    //} // close loop on adesso2

	} // close loop on first digit
	//} // close loop on adesso1


      break;

    case 4:

      adesso1 = 0;
      for (Int_t firstIndex=0; firstIndex<kMaxNumberOfDigitsPerVolume; firstIndex++) {
	if (selectedDigit[adesso1][firstIndex]==0x0) continue;

	adesso2 = 1;
	for (Int_t secondIndex=0; secondIndex<kMaxNumberOfDigitsPerVolume; secondIndex++) {
	  if (selectedDigit[adesso2][secondIndex]==0x0) continue;

	  adesso3 = 2;
	  for (Int_t thirdIndex=0; thirdIndex<kMaxNumberOfDigitsPerVolume; thirdIndex++) {
	    if (selectedDigit[adesso3][thirdIndex]==0x0) continue;

	    adesso4 = 3;
	    for (Int_t fourthIndex=0; fourthIndex<kMaxNumberOfDigitsPerVolume; fourthIndex++) {
	      if (selectedDigit[adesso4][fourthIndex]==0x0) continue;


	      if (TMath::Abs(selectedDigit[adesso1][firstIndex]->GetTDC()-selectedDigit[adesso2][secondIndex]->GetTDC())>fMaxDeltaTime
		  ||
		  TMath::Abs(selectedDigit[adesso1][firstIndex]->GetTDC()-selectedDigit[adesso3][thirdIndex]->GetTDC())>fMaxDeltaTime
		  ||
		  TMath::Abs(selectedDigit[adesso1][firstIndex]->GetTDC()-selectedDigit[adesso4][fourthIndex]->GetTDC())>fMaxDeltaTime
		  ||
		  TMath::Abs(selectedDigit[adesso2][secondIndex]->GetTDC()-selectedDigit[adesso3][thirdIndex]->GetTDC())>fMaxDeltaTime
		  ||
		  TMath::Abs(selectedDigit[adesso2][secondIndex]->GetTDC()-selectedDigit[adesso4][fourthIndex]->GetTDC())>fMaxDeltaTime
		  ||
		  TMath::Abs(selectedDigit[adesso3][thirdIndex]->GetTDC()-selectedDigit[adesso4][fourthIndex]->GetTDC())>fMaxDeltaTime) continue;

	      interestingTOF[0] = selectedDigit[adesso1][firstIndex]->GetTDC();
	      interestingTOT[0] = selectedDigit[adesso1][firstIndex]->GetTOT();
	      interestingADC[0] = selectedDigit[adesso1][firstIndex]->GetADC();
	      interestingWeight[0] = selectedDigit[adesso1][firstIndex]->GetWeight();
	      Int_t vol1[5]; for(jj=0; jj<5; jj++) vol1[jj] = selectedDigit[adesso1][firstIndex]->GetDetectorIndex(jj);
	      AliDebug(1,Form(" %2d %1d %2d %1d %2d", vol1[0], vol1[1], vol1[2], vol1[4], vol1[3]));
	      Int_t volDum = vol1[3];
	      vol1[3] = vol1[4];
	      vol1[4] = volDum;
	      fTOFGeometry->GetPosPar(vol1,pos);
	      interestingX[0] = pos[0];
	      interestingY[0] = pos[1];
	      interestingZ[0] = pos[2];

	      interestingTOF[1] = selectedDigit[adesso2][secondIndex]->GetTDC();
	      interestingTOT[1] = selectedDigit[adesso2][secondIndex]->GetTOT();
	      interestingADC[1] = selectedDigit[adesso2][secondIndex]->GetADC();
	      interestingWeight[1] = selectedDigit[adesso2][secondIndex]->GetWeight();
	      Int_t vol2[5]; for(jj=0; jj<5; jj++) vol2[jj] = selectedDigit[adesso2][secondIndex]->GetDetectorIndex(jj);
	      AliDebug(1,Form(" %2d %1d %2d %1d %2d", vol2[0], vol2[1], vol2[2], vol2[4], vol2[3]));
	      volDum = vol2[3];
	      vol2[3] = vol2[4];
	      vol2[4] = volDum;
	      fTOFGeometry->GetPosPar(vol2,pos);
	      interestingX[1] = pos[0];
	      interestingY[1] = pos[1];
	      interestingZ[1] = pos[2];

	      interestingTOF[2] = selectedDigit[adesso3][thirdIndex]->GetTDC();
	      interestingTOT[2] = selectedDigit[adesso3][thirdIndex]->GetTOT();
	      interestingADC[2] = selectedDigit[adesso3][thirdIndex]->GetADC();
	      interestingWeight[2] = selectedDigit[adesso3][thirdIndex]->GetWeight();
	      Int_t vol3[5]; for(jj=0; jj<5; jj++) vol3[jj] = selectedDigit[adesso3][thirdIndex]->GetDetectorIndex(jj);
	      AliDebug(1,Form(" %2d %1d %2d %1d %2d", vol3[0], vol3[1], vol3[2], vol3[4], vol3[3]));
	      volDum = vol3[3];
	      vol3[3] = vol3[4];
	      vol3[4] = volDum;
	      fTOFGeometry->GetPosPar(vol3,pos);
	      interestingX[2] = pos[0];
	      interestingY[2] = pos[1];
	      interestingZ[2] = pos[2];

	      interestingTOF[3] = selectedDigit[adesso4][fourthIndex]->GetTDC();
	      interestingTOT[3] = selectedDigit[adesso4][fourthIndex]->GetTOT();
	      interestingADC[3] = selectedDigit[adesso4][fourthIndex]->GetADC();
	      interestingWeight[3] = selectedDigit[adesso4][fourthIndex]->GetWeight();
	      Int_t vol4[5]; for(jj=0; jj<5; jj++) vol4[jj] = selectedDigit[adesso4][fourthIndex]->GetDetectorIndex(jj);
	      AliDebug(1,Form(" %2d %1d %2d %1d %2d", vol4[0], vol4[1], vol4[2], vol4[4], vol4[3]));
	      volDum = vol4[3];
	      vol4[3] = vol4[4];
	      vol4[4] = volDum;
	      fTOFGeometry->GetPosPar(vol4,pos);
	      interestingX[3] = pos[0];
	      interestingY[3] = pos[1];
	      interestingZ[3] = pos[2];


	      AverageCalculations(interestingCounter+1,
				  interestingX, interestingY, interestingZ,
				  interestingTOF, interestingTOT, interestingADC,
				  interestingWeight,
				  parTOF, posClus, check);

	      for (jj=0; jj<5; jj++) det[jj] = -1;
	      for (jj=0; jj<3; jj++) posF[jj] = posClus[jj];
	      fTOFGeometry->GetDetID(posF, det);

	      volIdClus = fTOFGeometry->GetAliSensVolIndex(det[0],det[1],det[2]);
	      //volIdClus = GetClusterVolIndex(det);

	      for (jj=3; jj<11; jj++) padsCluster[jj] = -1;
	      padsCluster[3] = selectedDigit[adesso1][firstIndex]->GetDetectorIndex(4);
	      padsCluster[4] = selectedDigit[adesso1][firstIndex]->GetDetectorIndex(3);
	      padsCluster[5] = selectedDigit[adesso2][secondIndex]->GetDetectorIndex(4);
	      padsCluster[6] = selectedDigit[adesso2][secondIndex]->GetDetectorIndex(3);
	      padsCluster[7] = selectedDigit[adesso3][thirdIndex]->GetDetectorIndex(4);
	      padsCluster[8] = selectedDigit[adesso3][thirdIndex]->GetDetectorIndex(3);
	      padsCluster[9] = selectedDigit[adesso4][fourthIndex]->GetDetectorIndex(4);
	      padsCluster[10] = selectedDigit[adesso4][fourthIndex]->GetDetectorIndex(3);

	      for (jj=0; jj<6; jj++) covClus[jj] = 0.;
	      Int_t ** indDet = new Int_t*[interestingCounter+1];
	      for (jj=0; jj<interestingCounter+1; jj++) indDet[jj] = new Int_t [5];
	      for (jj=0; jj<interestingCounter+1; jj++) indDet[jj][0] = nSector;
	      for (jj=0; jj<interestingCounter+1; jj++) indDet[jj][1] = nPlate;
	      for (jj=0; jj<interestingCounter+1; jj++) indDet[jj][2] = nStrip;
	      for (jj=0; jj<interestingCounter+1; jj++) indDet[jj][3] = padsCluster[2*jj+3];
	      for (jj=0; jj<interestingCounter+1; jj++) indDet[jj][4] = padsCluster[2*jj+1+3];
	      GetClusterPars(/*check,*/ interestingCounter+1, indDet, interestingWeight, posClus, covClus);
	      for (jj=0; jj<interestingCounter+1; jj++) delete [] indDet[jj];
	      delete [] indDet;

	      // To fill the track index array
	      dummyCounter=-1;
	      for (jj=0; jj<kMaxNumberOfTracksPerDigit; jj++) tracks[jj] = -1;
	      for (Int_t kk=0; kk<kMaxNumberOfTracksPerDigit; kk++) { // three is the max number of tracks associated to one digit
		if (selectedDigit[adesso1][firstIndex]->GetTrackLabel(kk)==-1) continue;
		else {
		  dummyCounter++;
		  tracks[dummyCounter] = selectedDigit[adesso1][firstIndex]->GetTrackLabel(kk);
		}
	      }
	      for (Int_t kk=0; kk<kMaxNumberOfTracksPerDigit; kk++) { // three is the max number of tracks associated to one digit
		if (selectedDigit[adesso2][secondIndex]->GetTrackLabel(kk)==-1) continue;
		else {

		  alreadyStored = kFALSE;
		  for (jj=0; jj<dummyCounter+1; jj++)
		    alreadyStored = alreadyStored || (tracks[jj]!=-1 && tracks[jj]==selectedDigit[adesso2][secondIndex]->GetTrackLabel(kk));

		  if (alreadyStored) continue;
		  if (dummyCounter==2) { // three is the max number of tracks associated to one cluster
		    AliWarning("  Siamo al limite!");
		    continue;
		  }

		  dummyCounter++;
		  tracks[dummyCounter] = selectedDigit[adesso2][secondIndex]->GetTrackLabel(kk);

		}

	      }
	      for (Int_t kk=0; kk<kMaxNumberOfTracksPerDigit; kk++) { // three is the max number of tracks associated to one digit
		if (selectedDigit[adesso3][thirdIndex]->GetTrackLabel(kk)==-1) continue;
		else {

		  alreadyStored = kFALSE;
		  for (jj=0; jj<dummyCounter+1; jj++)
		    alreadyStored = alreadyStored || (tracks[jj]!=-1 && tracks[jj]==selectedDigit[adesso3][thirdIndex]->GetTrackLabel(kk));

		  if (alreadyStored) continue;
		  if (dummyCounter==2) { // three is the max number of tracks associated to one cluster
		    AliWarning("  Siamo al limite!");
		    continue;
		  }

		  dummyCounter++;
		  tracks[dummyCounter] = selectedDigit[adesso3][thirdIndex]->GetTrackLabel(kk);

		}

	      }
	      for (Int_t kk=0; kk<kMaxNumberOfTracksPerDigit; kk++) { // three is the max number of tracks associated to one digit
		if (selectedDigit[adesso4][fourthIndex]->GetTrackLabel(kk)==-1) continue;
		else {

		  alreadyStored = kFALSE;
		  for (jj=0; jj<dummyCounter+1; jj++)
		    alreadyStored = alreadyStored || (tracks[jj]!=-1 && tracks[jj]==selectedDigit[adesso4][fourthIndex]->GetTrackLabel(kk));

		  if (alreadyStored) continue;
		  if (dummyCounter==2) { // three is the max number of tracks associated to one cluster
		    AliWarning("  Siamo al limite!");
		    continue;
		  }

		  dummyCounter++;
		  tracks[dummyCounter] = selectedDigit[adesso4][fourthIndex]->GetTrackLabel(kk);

		}

	      }


	      AliTOFcluster *tofCluster =
		new AliTOFcluster(volIdClus, posClus[0], posClus[1], posClus[2],
				  covClus[0], covClus[1], covClus[2], covClus[3], covClus[4], covClus[5],
				  tracks, det, parTOF, status, selectedDigit[adesso1][firstIndex]->GetIndex()); // to be updated
	      InsertCluster(tofCluster);

	      AliDebug(2, Form("       %4d  %f %f %f  %f %f %f %f %f %f  %3d %3d %3d  %2d %1d %2d %1d %2d  %4d %3d %3d %4d %4d  %1d  %4d", 
			       volIdClus, posClus[0], posClus[1], posClus[2],
			       covClus[0], covClus[1], covClus[2], covClus[3], covClus[4], covClus[5],
			       tracks[0], tracks[1], tracks[2],
			       det[0], det[1], det[2], det[3], det[4],
			       parTOF[0], parTOF[1], parTOF[2], parTOF[3], parTOF[4],
			       status, selectedDigit[adesso1][firstIndex]->GetIndex()));

	      volDum = vol1[3];
	      vol1[3] = vol1[4];
	      vol1[4] = volDum;
	      fTOFdigitMap->ResetDigitNumber(vol1,selectedDigit[adesso1][firstIndex]->GetIndex());
	      volDum = vol2[3];
	      vol2[3] = vol2[4];
	      vol2[4] = volDum;
	      fTOFdigitMap->ResetDigitNumber(vol2,selectedDigit[adesso2][secondIndex]->GetIndex());
	      volDum = vol3[3];
	      vol3[3] = vol3[4];
	      vol3[4] = volDum;
	      fTOFdigitMap->ResetDigitNumber(vol3,selectedDigit[adesso3][thirdIndex]->GetIndex());
	      volDum = vol4[3];
	      vol4[3] = vol4[4];
	      vol4[4] = volDum;
	      fTOFdigitMap->ResetDigitNumber(vol4,selectedDigit[adesso4][fourthIndex]->GetIndex());


	    } // close loop on fourth digit

	  } // close loop on third digit

	} // close loop on second digit

      } // close loop on first digit

      break;

    }

    for (ii=0; ii<kMaxNumberOfInterestingPads; ii++)
      for (jj=0; jj<kMaxNumberOfDigitsPerVolume; jj++)
	selectedDigit[ii][jj] = 0x0;

  } // loop on iPad

  for (ii=0; ii<kMaxNumberOfInterestingPads; ii++) {
    for (jj=0; jj<kMaxNumberOfDigitsPerVolume; jj++) {
      delete [] selectedDigit[ii][jj];
      selectedDigit[ii][jj] = 0x0;
    }
    delete [] selectedDigit[ii];
    selectedDigit[ii] = 0x0;
  }
  delete [] selectedDigit;
  selectedDigit = 0x0;

}
//_____________________________________________________________________________

Int_t AliTOFClusterFinderV1::InsertCluster(AliTOFcluster *tofCluster)
{
  //
  // This function adds a TOF cluster to the array of TOF clusters
  // sorted in Z, i.e. fTofClusters
  //

  if (fNumberOfTofClusters==kTofMaxCluster) {
    AliError("Too many clusters !");
    return 1;
  }

  if (fNumberOfTofClusters==0) {
    fTofClusters[fNumberOfTofClusters++] = tofCluster;
    return 0;
  }

  Int_t ii = FindClusterIndex(tofCluster->GetZ());
  memmove(fTofClusters+ii+1 ,fTofClusters+ii,(fNumberOfTofClusters-ii)*sizeof(AliTOFcluster*));
  fTofClusters[ii] = tofCluster;
  fNumberOfTofClusters++;
  
  return 0;

}
//_____________________________________________________________________________

Int_t AliTOFClusterFinderV1::FindClusterIndex(Double_t z) const
{
  //
  // This function returns the index of the nearest cluster in z
  //

  if (fNumberOfTofClusters==0) return 0;
  if (z <= fTofClusters[0]->GetZ()) return 0;
  if (z > fTofClusters[fNumberOfTofClusters-1]->GetZ()) return fNumberOfTofClusters;
  Int_t b = 0, e = fNumberOfTofClusters-1, m = (b+e)/2;
  for (; b<e; m=(b+e)/2) {
    if (z > fTofClusters[m]->GetZ()) b=m+1;
    else e=m;
  }

  return m;

}
//_____________________________________________________________________________

void AliTOFClusterFinderV1::ResetRecpoint()
{
  //
  // Clear the list of reconstructed points
  //

  fNumberOfTofClusters = 0;
  if (fRecPoints) fRecPoints->Clear();

}
//_____________________________________________________________________________

void AliTOFClusterFinderV1::ResetDigits()
{
  //
  // Clear the list of digits
  //

  fNumberOfTofDigits = 0;
  if (fDigits) fDigits->Clear();

}
//_____________________________________________________________________________
//UShort_t AliTOFClusterFinderV1::GetClusterVolIndex(Int_t *ind) const
//{
  //
  // Get the volume ID to retrieve the l2t transformation
  //

  // Detector numbering scheme
/*
  Int_t nSector = AliTOFGeometry::NSectors();
  Int_t nPlate  = AliTOFGeometry::NPlates();
  Int_t nStripA = AliTOFGeometry::NStripA();
  Int_t nStripB = AliTOFGeometry::NStripB();
  Int_t nStripC = AliTOFGeometry::NStripC();

  Int_t isector =ind[0];
  if (isector >= nSector)
    AliError(Form("Wrong sector number in TOF (%d) !", isector));
  Int_t iplate = ind[1];
  if (iplate >= nPlate)
    AliError(Form("Wrong plate number in TOF (%d) !", iplate));
  Int_t istrip = ind[2];

  Int_t stripOffset = 0;
  switch (iplate) {
  case 0:
    stripOffset = 0;
    break;
  case 1:
    stripOffset = nStripC;
    break;
  case 2:
    stripOffset = nStripC+nStripB;
    break;
  case 3:
    stripOffset = nStripC+nStripB+nStripA;
    break;
  case 4:
    stripOffset = nStripC+nStripB+nStripA+nStripB;
    break;
  default:
    AliError(Form("Wrong plate number in TOF (%d) !", iplate));
    break;
  };

  Int_t index= (2*(nStripC+nStripB)+nStripA)*isector +
               stripOffset +
               istrip;

  UShort_t volIndex = AliGeomManager::LayerToVolUID(AliGeomManager::kTOF, index);
  return volIndex;

}
*/
//_____________________________________________________________________________

void AliTOFClusterFinderV1::GetClusterPars(Int_t *ind, Double_t* pos, Double_t* cov) const
{
  //
  // Starting from the volume indices (ind[5]), for a cluster coming from
  // a isolated digits, this function returns:
  //   the cluster position (pos),
  //   the cluster covariance matrix elements (cov)
  // 

  //
  //we now go in the system of the strip: determine the local coordinates
  //
  //
  // 47---------------------------------------------------0  ^ z
  // | | | | | | | | | | | | | | | | | | | | | | | | | | | 1 |
  // -----------------------------------------------------   | y going outwards
  // | | | | | | | | | | | | | | | | | | | | | | | | | | | 0 |  par[0]=0;

  // -----------------------------------------------------   |
  // x <-----------------------------------------------------

  //move to the tracking ref system
  Double_t lpos[3] = { (ind[4]-23.5)*AliTOFGeometry::XPad(),
		       0.,
		       (ind[3]- 0.5)*AliTOFGeometry::ZPad() };   
  AliDebug(1, Form(" %f %f %f", lpos[0], lpos[1], lpos[2]));

  // Volume ID
  UShort_t volIndex = fTOFGeometry->GetAliSensVolIndex(ind[0],ind[1],ind[2]);
  //UShort_t volIndex = GetClusterVolIndex(ind);
  const TGeoHMatrix *l2t = AliGeomManager::GetTracking2LocalMatrix(volIndex);

  // Get the position in the track ref system
  Double_t tpos[3];
  l2t->MasterToLocal(lpos,tpos);
  pos[0] = tpos[0];
  pos[1] = tpos[1];
  pos[2] = tpos[2];

  //Get the cluster covariance in the track ref system
  Double_t lcov[9];
  for (Int_t ii=0; ii<9; ii++) lcov[ii] = 0.;

  //cluster covariance in the local system:
  // sx2   0   0
  // 0     0   0
  // 0     0 sz2
  /*
  lcov[4] = 0.42*0.42/3.;
                       // = ( 5*0.025 (gas gaps thikness)
                       //   + 4*0.040 (internal glasses thickness)
                       //   + 0.5*0.160 (internl PCB)
                       //   + 1*0.055 (external red glass))
  */

  lcov[0] = 0.499678;//AliTOFGeometry::XPad()*AliTOFGeometry::XPad()/12.;
  lcov[8] = 0.992429;//AliTOFGeometry::ZPad()*AliTOFGeometry::ZPad()/12.;

  //cluster covariance in the tracking system:
  TGeoHMatrix m;
  m.SetRotation(lcov);
  m.Multiply(l2t);
  m.MultiplyLeft(&l2t->Inverse());
  Double_t *tcov = m.GetRotationMatrix();
  cov[0] = tcov[0]; cov[1] = tcov[1]; cov[2] = tcov[2];
  cov[3] = tcov[4]; cov[4] = tcov[5]; cov[5] = tcov[8];

  return;

}
//_____________________________________________________________________________

void AliTOFClusterFinderV1::GetClusterPars(/*Bool_t check,*/ Int_t counter,
					   Int_t **ind, Double_t *weight,
					   Double_t *pos, Double_t *cov) const
{
  //
  // Starting from:
  //               the volumes indices (ind[counter][5]), for a
  //                  cluster coming from a collection of 'counter'
  //                  digits,
  //               the volumes weights (weight[counter]), -controlled
  //                  by the 'check' variable control-, for a cluster
  //                  coming from a collection of 'counter' digits,
  //               the cluster position (pos),
  // this function returns:
  //   the covariance matrix elements (cov) for the found cluster
  //

  //
  // we now go in the system of the strip: determine the local coordinates
  //
  // 47---------------------------------------------------0  ^ z
  // | | | | | | | | | | | | | | | | | | | | | | | | | | | 1 |
  // -----------------------------------------------------   | y going outwards
  // | | | | | | | | | | | | | | | | | | | | | | | | | | | 0 |  par[0]=0;

  // -----------------------------------------------------   |
  // x <-----------------------------------------------------

  for (Int_t ii=0; ii<counter; ii++)
    AliDebug(1, Form(" %2d  %2d %1d %2d %1d %2d ",
		     ii, ind[ii][0], ind[ii][1], ind[ii][2], ind[ii][3], ind[ii][4]));

  Float_t posF[3]; for (Int_t ii=0; ii<3; ii++) posF[ii] = (Float_t)pos[ii];
  AliDebug(1, Form(" %f %f %f", pos[0], pos[1], pos[2]));

  Int_t detClus[5] = {-1, -1, -1, -1, -1};
  fTOFGeometry->GetDetID(posF, detClus);

  // Volume ID
  UShort_t volIndex = fTOFGeometry->GetAliSensVolIndex(detClus[0],detClus[1],detClus[2]);
  //UShort_t volIndex = GetClusterVolIndex(detClus);
  AliDebug(1, Form(" %2d %1d %2d %1d %2d  %7i",
		   detClus[0], detClus[1], detClus[2], detClus[3], detClus[4], volIndex));

  // Get the position in the TOF strip ref system
  const TGeoHMatrix *alice2strip = AliGeomManager::GetOrigGlobalMatrix(volIndex);
  Double_t ppos[3] = {-1, -1, -1};
  alice2strip->MasterToLocal(pos,ppos);
  AliDebug(1, Form(" %f %f %f", ppos[0], ppos[1], ppos[2]));


  // Get the position in the tracking ref system
  const TGeoHMatrix *g2l = AliGeomManager::GetTracking2LocalMatrix(volIndex);
  Double_t lpos[3] = {-1, -1, -1};
  g2l->MasterToLocal(ppos,lpos);
  AliDebug(1, Form(" %f %f %f", lpos[0], lpos[1], lpos[2]));
  for (Int_t ii=0; ii<3; ii++) pos[ii] = lpos[ii];

  //Get the cluster covariance in the track ref system
  Double_t lcov[9];
  for (Int_t ii=0; ii<9; ii++) lcov[ii] = 0.;

  //cluster covariance in the local system:
  // sx2   0   0
  // 0     0   0
  // 0     0 sz2

  // Evaluation of the ovariance matrix elements
  TOFclusterError(/*check,*/ counter, ind, weight, ppos, lcov);

  AliDebug(1, Form("lcov[0] = %f, lcov[8] = %f", lcov[0], lcov[8]));

  //cluster covariance in the tracking system:
  TGeoHMatrix m;
  m.SetRotation(lcov);
  m.Multiply(g2l);
  m.MultiplyLeft(&g2l->Inverse());
  Double_t *tcov = m.GetRotationMatrix();
  cov[0] = tcov[0]; cov[1] = tcov[1]; cov[2] = tcov[2];
  cov[3] = tcov[4]; cov[4] = tcov[5]; cov[5] = tcov[8];

  return;

}
//_____________________________________________________________________________

void AliTOFClusterFinderV1::TOFclusterError(/*Bool_t check,*/ Int_t counter,
					    Int_t **ind, Double_t *weight,
					    Double_t ppos[], Double_t lcov[]) const
{
  //
  //
  //

  //lcov[4] = 0.42*0.42/3.; // cm2
                       // = ( 5*0.025 (gas gaps thikness)
                       //   + 4*0.040 (internal glasses thickness)
                       //   + 0.5*0.160 (internl PCB)
                       //   + 1*0.055 (external red glass))


  Float_t *delta2X = new Float_t[counter];
  for (Int_t ii=0; ii<counter; ii++)
    delta2X[ii] =
      ((ind[ii][4]-23.5)*AliTOFGeometry::XPad() - ppos[0])*((ind[ii][4]-23.5)*AliTOFGeometry::XPad() - ppos[0]);

  Float_t *delta2Z = new Float_t[counter];
  for (Int_t ii=0; ii<counter; ii++)
    delta2Z[ii] =
      ((ind[ii][3]- 0.5)*AliTOFGeometry::ZPad() - ppos[2])*((ind[ii][3]- 0.5)*AliTOFGeometry::ZPad() - ppos[2]);

  for (Int_t ii=0; ii<counter; ii++)
    AliDebug(1, Form("x[%d] = %f, z[%d] = %f, weight[%d] = %f",
		     ii, (ind[ii][4]-23.5)*AliTOFGeometry::XPad(),
		     ii, (ind[ii][3]- 0.5)*AliTOFGeometry::ZPad(),
		     ii, weight[ii]
		     ));
  AliDebug(1, Form("xMean = %f, zMean = %f",ppos[0], ppos[2]));


  switch (counter)
    {

    case 2:

      if (ind[0][3]==ind[1][3] && TMath::Abs(ind[0][4]-ind[1][4])==1) { //

	lcov[8] = 1.02039; // cm2
	lcov[0] = 0.0379409; // cm2
	/*
	if (check)
	  lcov[0] = 0.5*0.5; // cm2
	else {
	  if (weight[0]==weight[1])
	    lcov[0] = 0.0379409; // cm2
	  else
	    lcov[0] = TMath::Mean(counter, delta2X, weight); // cm2
	}
	*/

      }

      else if (ind[0][4]==ind[1][4] && TMath::Abs(ind[0][3]-ind[1][3])==1) {//

	lcov[0] = 0.505499; // cm2
	lcov[8] = 0.0422046; // cm2
	/*
	if (check)
	  lcov[8] = 0.5*0.5; // cm2
	else {
	  if (weight[0]==weight[1])
	    lcov[8] = 0.0422046; // cm2
	  else
	    lcov[8] = TMath::Mean(counter, delta2Z, weight); // cm2
	}
	*/

      }

      break;

    case 3:
      lcov[0] = 0.0290677; // cm2
      lcov[8] = 0.0569726; // cm2
      /*
      if (check) {
	lcov[0] = 0.5*0.5; // cm2
	lcov[8] = 0.5*0.5; // cm2
      }
      else {
      if (weight[0]==weight[1] && weight[0]==weight[2]) {
	  lcov[0] = 0.0290677; // cm2
	  lcov[8] = 0.0569726; // cm2
	  }
	else {
	  lcov[0] = TMath::Mean(counter, delta2X, weight); // cm2
	  lcov[8] = TMath::Mean(counter, delta2Z, weight); // cm2
	  }

	}
      */

      break;

    case 4:
      lcov[0] = 0.0223807; // cm2
      lcov[8] = 0.0438662; // cm2
      /*
      if (check) {
	lcov[0] = 0.5*0.5; // cm2
	lcov[8] = 0.5*0.5; // cm2
      }
      else {
      if (weight[0]==weight[1] && weight[0]==weight[2] && weight[0]==weight[3]) {
	  lcov[0] = 0.0223807; // cm2
	  lcov[8] = 0.0438662; // cm2
	  }
	else {
	  lcov[0] = TMath::Mean(counter, delta2X, weight); // cm2
	  lcov[8] = TMath::Mean(counter, delta2Z, weight); // cm2
	  }

	}
      */

      break;

    }

  delete [] delta2Z;
  delete [] delta2X;

}
//_____________________________________________________________________________

Bool_t AliTOFClusterFinderV1::MakeSlewingCorrection(Int_t *detectorIndex,
						    Int_t tofDigitToT,
						    Int_t tofDigitTdc,
						    Int_t &tdcCorr)
{
  //
  // This funtion makes the following:
  //
  //      - if at least one of the three status (Pulser/Noise/HW) is
  //        bad, is sets the status of electronic channel, corresponding to the
  //        volume identified by detectorIndex, as kFALSE;
  //      - if offline calibration is in the valid status, it performs the
  //        slewing correction. In particular, by taking into account:
  //          * the measured tot and tof values (tofDigitToT and tofDigitTdc,
  //            respectively);
  //          * the six parameters of 5th order polynomial used
  //            to fit the tofVStot scatter plot,
  //         it returns the corrected tof value, i.e. tdcCorr value.
  //

  Bool_t output = kTRUE;

  Double_t timeCorr;
  Int_t jj;

  //AliInfo(" Calibrating TOF Digits: ");
  
  AliTOFChannelOnlineArray *calDelay = fTOFcalib->GetTOFOnlineDelay();
  AliTOFChannelOnlineStatusArray *calStatus = fTOFcalib->GetTOFOnlineStatus();

  TObjArray *calTOFArrayOffline = fTOFcalib->GetTOFCalArrayOffline();

  Int_t index = AliTOFGeometry::GetIndex(detectorIndex);

  UChar_t statusPulser = calStatus->GetPulserStatus(index);
  UChar_t statusNoise  = calStatus->GetNoiseStatus(index);
  UChar_t statusHW     = calStatus->GetHWStatus(index);
  UChar_t status       = calStatus->GetStatus(index);

  //check the status, also unknown is fine!!!!!!!

  AliDebug(2, Form(" Status for channel %d = %d",index, (Int_t)status));
  if((statusPulser & AliTOFChannelOnlineStatusArray::kTOFPulserBad)==(AliTOFChannelOnlineStatusArray::kTOFPulserBad)||(statusNoise & AliTOFChannelOnlineStatusArray::kTOFNoiseBad)==(AliTOFChannelOnlineStatusArray::kTOFNoiseBad)||(statusHW & AliTOFChannelOnlineStatusArray::kTOFHWBad)==(AliTOFChannelOnlineStatusArray::kTOFHWBad)){
    AliDebug(2, Form(" Bad Status for channel %d",index));
    //fTofClusters[ii]->SetStatus(kFALSE); //odd convention, to avoid conflict with calibration objects currently in the db (temporary solution).
    output = kFALSE;
  }
  else
    AliDebug(2, Form(" Good Status for channel %d",index));


  if (fCalibrateTOFtimes) { // AdC

  // Get Rough channel online equalization 
  Double_t roughDelay = (Double_t)calDelay->GetDelay(index);  // in ns
  AliDebug(2,Form(" channel delay (ns) = %f", roughDelay));
  // Get Refined channel offline calibration parameters
  TString validity = (TString)fTOFcalib->GetOfflineValidity();
  if (validity.CompareTo("valid")==0) {
    AliTOFChannelOffline * calChannelOffline = (AliTOFChannelOffline*)calTOFArrayOffline->At(index);
    Double_t par[6];
    for (jj = 0; jj<6; jj++)
      par[jj] = (Double_t)calChannelOffline->GetSlewPar(jj);

    AliDebug(2,Form(" Calib Pars = %f, %f, %f, %f, %f, %f ",par[0],par[1],par[2],par[3],par[4],par[5]));
    AliDebug(2,Form(" The ToT and Time, uncorr (counts) = %d , %d", tofDigitToT, tofDigitTdc));
    Double_t tToT = (Double_t)(tofDigitToT*AliTOFGeometry::ToTBinWidth());    
    tToT*=1.E-3; //ToT in ns
    AliDebug(2,Form(" The ToT and Time, uncorr (ns)= %e, %e",tofDigitTdc*AliTOFGeometry::TdcBinWidth()*1.E-3,tToT));
    timeCorr = par[0] + tToT*(par[1] + tToT*(par[2] + tToT*(par[3] + tToT*(par[4] + tToT*par[5])))); // the time correction (ns)
  }
  else
    timeCorr = roughDelay; // correction in ns

  AliDebug(2,Form(" The ToT and Time, uncorr (ns)= %e, %e",tofDigitTdc*AliTOFGeometry::TdcBinWidth()*1.E-3,tofDigitToT*AliTOFGeometry::ToTBinWidth()));
  AliDebug(2,Form(" The time correction (ns) = %f", timeCorr));
  timeCorr = (Double_t)(tofDigitTdc)*AliTOFGeometry::TdcBinWidth()*1.E-3-timeCorr;//redefine the time
  timeCorr *= 1.E3;
  AliDebug(2,Form(" The channel time, corr (ps)= %e",timeCorr ));
  //tdcCorr = (Int_t)(timeCorr/AliTOFGeometry::TdcBinWidth()); //the corrected time (tdc counts)
  tdcCorr = TMath::Nint(timeCorr/AliTOFGeometry::TdcBinWidth()); //the corrected time (tdc counts)
  
  } // AdC

  return output;

}
//______________________________________________________________________________

void AliTOFClusterFinderV1::Digits2RecPoints(Int_t iEvent)
{
  //
  // Converts digits to recpoints for TOF
  //

  TStopwatch stopwatch;
  stopwatch.Start();

  fRunLoader->GetEvent(iEvent);

  AliLoader *localTOFLoader = (AliLoader*)fRunLoader->GetLoader("TOFLoader");

  TTree * localTreeD = (TTree*)localTOFLoader->TreeD();
  if (localTreeD == 0x0) {
    AliFatal("Can not get TreeD");
    return;
  }

  TBranch *branch = localTreeD->GetBranch("TOF");
  if (!branch) {
    AliError("Can't get the branch with the TOF digits !");
    return;
  }

  TTree *localTreeR = (TTree*)localTOFLoader->TreeR();
  if (localTreeR == 0x0)
    {
      localTOFLoader->MakeTree("R");
      localTreeR = localTOFLoader->TreeR();
    }

  Digits2RecPoints(localTreeD, localTreeR);

  //localTOFLoader = fRunLoader->GetLoader("TOFLoader");  
  localTOFLoader->WriteRecPoints("OVERWRITE");

  AliDebug(1, Form("Execution time to read TOF digits and to write TOF clusters for the event number %d: R:%.4fs C:%.4fs",
		   iEvent, stopwatch.RealTime(),stopwatch.CpuTime()));

}
//______________________________________________________________________________

void AliTOFClusterFinderV1::Digits2RecPoints(Int_t iEvent, AliRawReader *rawReader)
{
  //
  // Converts RAW data to recpoints for TOF
  //

  TStopwatch stopwatch;
  stopwatch.Start();

  fRunLoader->GetEvent(iEvent);

  AliDebug(2,Form(" Event number %2d ", iEvent));

  AliLoader *localTOFLoader = (AliLoader*)fRunLoader->GetLoader("TOFLoader");

  TTree *localTreeR = localTOFLoader->TreeR();

  if (localTreeR == 0x0){
    localTOFLoader->MakeTree("R");
    localTreeR = localTOFLoader->TreeR();
  }

  Digits2RecPoints(rawReader, localTreeR);

  AliDebug(1, Form("Execution time to read TOF raw data and to write TOF clusters for the event number %d: R:%.4fs C:%.4fs",
		   iEvent, stopwatch.RealTime(),stopwatch.CpuTime()));

}
//______________________________________________________________________________

void AliTOFClusterFinderV1::Raw2Digits(Int_t iEvent, AliRawReader *rawReader)
{
  //
  // Converts RAW data to MC digits for TOF
  //


  TStopwatch stopwatch;
  stopwatch.Start();

  fRunLoader->GetEvent(iEvent);

  AliDebug(2,Form(" Event number %2d ", iEvent));

  AliLoader *localTOFLoader = (AliLoader*)fRunLoader->GetLoader("TOFLoader");

  TTree *localTreeD = localTOFLoader->TreeD();

  if (localTreeD == 0x0){
    localTOFLoader->MakeTree("D");
    localTreeD = localTOFLoader->TreeD();
  }

  Raw2Digits(rawReader, localTreeD);

  AliDebug(1, Form("Execution time to read TOF raw data and to write TOF clusters for the event number %d: R:%.4fs C:%.4fs",
		   iEvent, stopwatch.RealTime(),stopwatch.CpuTime()));

}
//______________________________________________________________________________

void AliTOFClusterFinderV1::AverageCalculations(Int_t number, Float_t *interestingX,
						Float_t *interestingY, Float_t *interestingZ,
						Double_t *interestingTOF, Double_t *interestingTOT,
						Double_t *interestingADC, Double_t *interestingWeight,
						Int_t *parTOF, Double_t *posClus, Bool_t &check)
{
  //
  // Calculates the mean values for cluster position (x,y,z),
  //  TOF charge and time
  //

  Double_t tofAverage = 0.;
  Double_t totAverage = 0.;
  Double_t adcAverage = 0.;

  check = kFALSE;
  Int_t ii=-1;
  for (ii=number-1; ii>=0; ii--) check=check||(interestingWeight[ii]==0 || interestingWeight[ii]==-1);

  if (check) {
		  
    posClus[0] = TMath::Mean(number, interestingX);
    posClus[1] = TMath::Mean(number, interestingY);
    posClus[2] = TMath::Mean(number, interestingZ);
    tofAverage = TMath::Mean(number, interestingTOF);
    totAverage = TMath::Mean(number, interestingTOT);
    adcAverage = TMath::Mean(number, interestingADC);

  }
  else {

    posClus[0] = TMath::Mean(number, interestingX, interestingWeight);
    posClus[1] = TMath::Mean(number, interestingY, interestingWeight);
    posClus[2] = TMath::Mean(number, interestingZ, interestingWeight);
    tofAverage = TMath::Mean(number, interestingTOF, interestingWeight);
    totAverage = TMath::Mean(number, interestingTOT, interestingWeight);
    adcAverage = TMath::Mean(number, interestingADC, interestingWeight);

  }

  parTOF[0] = Int_t(tofAverage);
  parTOF[1] = Int_t(totAverage);
  parTOF[2] = Int_t(adcAverage);
  parTOF[3] = Int_t(tofAverage);//tofND
  parTOF[4] = Int_t(tofAverage);//tofRAW
  parTOF[5] = 0;
  parTOF[6] = 0;

}
