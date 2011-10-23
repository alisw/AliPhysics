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
$Log: AliTOFClusterFinder.cxx,v $
Revision 1.31  2007/11/24 14:53:19  zampolli
Status flag implemented as UChar_t

Revision 1.30  2007/10/04 13:08:52  arcelli
updates to comply with AliTOFGeometryV5 becoming AliTOFGeometry

Revision 1.29  2007/10/03 10:42:33  arcelli
updates to handle new AliTOFcluster, inheriting form AliCluster3D

Revision 1.28  2007/05/31 16:06:05  arcelli
move instance of AliRawStream outside loop on DDL

Revision 1.27  2007/05/02 16:31:49  arcelli
Add methods to handle single event reconstruction.  retrieval of Calib
info moved to AliTOFReconstructor ctor, and passed via a pointer to
AliTOFcalib

Revision 1.26  2007/04/30 19:02:24  arcelli
hopefully the last refinements for correct type conversion in calibration

Revision 1.25  2007/04/30 15:22:17  arcelli
Change TOF digit Time, Tot etc to int type

Revision 1.24  2007/04/27 11:19:31  arcelli
updates for the new decoder

Revision 1.23  2007/04/23 16:51:39  decaro
Digits-to-raw_data conversion: correction for a more real description
(A.De Caro, R.Preghenella)

Revision 1.22  2007/04/19 17:26:32  arcelli
Fix a bug (add some debug printout

Revision 1.21  2007/04/18 17:28:12  arcelli
Set the ToT bin width to the one actually used...

Revision 1.20  2007/03/09 09:57:23  arcelli
 Remove a forgotten include of Riostrem

Revision 1.19  2007/03/08 15:41:20  arcelli
set uncorrected times when filling RecPoints

Revision 1.18  2007/03/06 16:31:20  arcelli
Add Uncorrected TOF Time signal

Revision 1.17  2007/02/28 18:09:11  arcelli
Add protection against failed retrieval of the CDB cal object,
now Reconstruction exits with AliFatal

Revision 1.16  2007/02/20 15:57:00  decaro
Raw data update: to read the TOF raw data defined in UNPACKED mode


Revision 0.03  2005/07/28 A. De Caro
         Implement public method
	 Raw2Digits(Int_t, AliRawReader *)
	 to convert digits from raw data in MC digits
	 (temporary solution)

Revision 0.02  2005/07/27 A. De Caro
         Implement public method
	 Digits2RecPoint(Int_t)
	 to convert digits in clusters

Revision 0.02  2005/07/26 A. De Caro
         Implement private methods
	 InsertCluster(AliTOFcluster *)
	 FindClusterIndex(Double_t)
	 originally implemented in AliTOFtracker
	 by S. Arcelli and C. Zampolli

Revision 0.01  2005/07/25 A. De Caro
         Implement public methods
	 Digits2RecPoint(AliRawReader *, TTree *)
	 Digits2RecPoint(Int_t, AliRawReader *)
	 to convert raw data in clusters
 */

////////////////////////////////////////////////////////////////
//                                                            //
//         Class for TOF cluster finder                       //
//                                                            //
// Starting from Raw Data, create rec points,                 //
//                         fill TreeR for TOF,                //
//                         write TOF.RecPoints.root file      //
//                                                            //
////////////////////////////////////////////////////////////////

#include "Riostream.h"

#include "TClonesArray.h"
#include "TStopwatch.h"
#include "TTree.h"
//#include <TGeoManager.h>
#include <TGeoMatrix.h>
//#include <TGeoPhysicalNode.h>

#include "AliDAQ.h"
#include "AliLoader.h"
#include "AliLog.h"
#include "AliRawReader.h"
#include "AliRunLoader.h"
//#include "AliAlignObj.h"
#include <AliGeomManager.h>

#include "AliTOFcalib.h"
#include "AliTOFChannelOnlineArray.h"
#include "AliTOFChannelOnlineStatusArray.h"
#include "AliTOFChannelOffline.h"
#include "AliTOFClusterFinder.h"
#include "AliTOFcluster.h"
#include "AliTOFdigit.h"
#include "AliTOFGeometry.h"
#include "AliTOFrawData.h"

#include "AliTOFDeltaBCOffset.h"
#include "AliTOFCTPLatency.h"
#include "AliTOFRunParams.h"

//extern TFile *gFile;

ClassImp(AliTOFClusterFinder)

AliTOFClusterFinder::AliTOFClusterFinder(AliTOFcalib *calib):
  TNamed("AliTOFClusterFinder",""),
  fRunLoader(0),
  fTOFLoader(0),
  fTreeD(0),
  fTreeR(0),
  fDigits(new TClonesArray("AliTOFdigit", 4000)),
  fRecPoints(new TClonesArray("AliTOFcluster", 4000)),
  fNumberOfTofClusters(0),
  fVerbose(0),
  fDecoderVersion(0),
  fTOFcalib(calib),
  fTOFRawStream(AliTOFRawStream())
{
//
// Constructor
//

  for (Int_t ii=0; ii<kTofMaxCluster; ii++) fTofClusters[ii]=0x0;

  TString validity = (TString)fTOFcalib->GetOfflineValidity();
  if (validity.CompareTo("valid")==0) {
    AliInfo(Form(" validity = %s - Using offline calibration parameters",validity.Data()));
  } else {
    AliInfo(Form(" validity = %s - Using online calibration parameters",validity.Data()));
  }

}

//______________________________________________________________________________

AliTOFClusterFinder::AliTOFClusterFinder(AliRunLoader* runLoader, AliTOFcalib *calib):
  TNamed("AliTOFClusterFinder",""),
  fRunLoader(runLoader),
  fTOFLoader(runLoader->GetLoader("TOFLoader")),
  fTreeD(0),
  fTreeR(0),
  fDigits(new TClonesArray("AliTOFdigit", 4000)),
  fRecPoints(new TClonesArray("AliTOFcluster", 4000)),
  fNumberOfTofClusters(0),
  fVerbose(0),
  fDecoderVersion(0),
  fTOFcalib(calib),
  fTOFRawStream(AliTOFRawStream())
{
//
// Constructor
//

  for (Int_t ii=0; ii<kTofMaxCluster; ii++) fTofClusters[ii]=0x0;

  TString validity = (TString)fTOFcalib->GetOfflineValidity();
  if (validity.CompareTo("valid")==0) {
    AliInfo(Form(" validity = %s - Using offline calibration parameters",validity.Data()));
  } else {
    AliInfo(Form(" validity = %s - Using online calibration parameters",validity.Data()));
  }

}

//------------------------------------------------------------------------
AliTOFClusterFinder::AliTOFClusterFinder(const AliTOFClusterFinder &source) :
  TNamed(source),
  fRunLoader(0),
  fTOFLoader(0),
  fTreeD(0),
  fTreeR(0),
  fDigits(source.fDigits),
  fRecPoints(source.fRecPoints),
  fNumberOfTofClusters(0),
  fVerbose(0),
  fDecoderVersion(source.fDecoderVersion),
  fTOFcalib(source.fTOFcalib),
  fTOFRawStream(source.fTOFRawStream)
{
  // copy constructor

  for (Int_t ii=0; ii<kTofMaxCluster; ii++) fTofClusters[ii]=source.fTofClusters[ii];

}

//------------------------------------------------------------------------
AliTOFClusterFinder& AliTOFClusterFinder::operator=(const AliTOFClusterFinder &source)
{
  // ass. op.

  if (this == &source)
    return *this;

  TNamed::operator=(source);  
  fDigits=source.fDigits;
  fRecPoints=source.fRecPoints;
  fVerbose=source.fVerbose;
  fDecoderVersion=source.fDecoderVersion;
  fTOFcalib=source.fTOFcalib;
  fTOFRawStream=source.fTOFRawStream;
  for (Int_t ii=0; ii<kTofMaxCluster; ii++) fTofClusters[ii]=source.fTofClusters[ii];

  return *this;

}
//______________________________________________________________________________

AliTOFClusterFinder::~AliTOFClusterFinder()
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

  //if (fTofClusters || fNumberOfTofClusters) {
  if (fNumberOfTofClusters) {
    for (Int_t ii=0; ii<kTofMaxCluster; ii++) {
      if (fTofClusters[ii]) fTofClusters[ii]->Delete();
      delete fTofClusters[ii];
    }
    fNumberOfTofClusters=0;
   }

}
//______________________________________________________________________________

void AliTOFClusterFinder::Digits2RecPoints(Int_t iEvent)
{
  //
  // Converts digits to recpoints for TOF
  //

  TStopwatch stopwatch;
  stopwatch.Start();

  Int_t inholes = 0;

  fRunLoader->GetEvent(iEvent);

  fTreeD = fTOFLoader->TreeD();
  if (fTreeD == 0x0) {
    AliFatal("AliTOFClusterFinder: Can not get TreeD");
    return;
  }

  fDigits->Clear();
  fTreeD->GetBranch("TOF")->SetAutoDelete(kFALSE);
  fTreeD->SetBranchAddress("TOF",&fDigits);

  ResetRecpoint();

  fTreeR = fTOFLoader->TreeR();
  if (fTreeR == 0x0)
    {
      fTOFLoader->MakeTree("R");
      fTreeR = fTOFLoader->TreeR();
    }

  Int_t bufsize = 32000;
  fTreeR->Branch("TOF", &fRecPoints, bufsize);

  fTreeD->GetEvent(0);
  Int_t nDigits = fDigits->GetEntriesFast();
  AliDebug(2,Form("Number of TOF digits: %d",nDigits));

  Int_t ii;
  Int_t dig[5]={-1,-1,-1,-1,-1}; //cluster detector indeces
  Int_t parTOF[7]={0,0,0,0,0,0,0}; //The TOF signal parameters
  Bool_t status=kTRUE; // assume all sim channels ok in the beginning...
  for (ii=0; ii<nDigits; ii++) {
    AliTOFdigit *d = (AliTOFdigit*)fDigits->UncheckedAt(ii);
    dig[0]=d->GetSector();
    dig[1]=d->GetPlate();
    dig[2]=d->GetStrip();
    dig[3]=d->GetPadz();
    dig[4]=d->GetPadx();

    /* check valid index */
    if (dig[0]==-1||dig[1]==-1||dig[2]==-1||dig[3]==-1||dig[4]==-1) continue;

    // Do not reconstruct anything in the holes
    if (dig[0]==13 || dig[0]==14 || dig[0]==15 ) { // sectors with holes
      if (dig[1]==2) { // plate with holes
	inholes++;
	continue;
      }
    }

    AliDebug(2,Form(" %2d  %1d  %2d  %1d  %2d ",dig[0],dig[1],dig[2],dig[3],dig[4]));

    parTOF[0] = d->GetTdc(); //the TDC signal
    parTOF[1] = d->GetToT(); //the ToT signal
    parTOF[2] = d->GetAdc(); // the adc charge
    parTOF[3] = d->GetTdcND(); // non decalibrated sim time
    parTOF[4] = d->GetTdc(); // raw time, == Tdc time for the moment
    parTOF[5] = 0; // deltaBC
    parTOF[6] = 0; // L0-L1 latency
    Double_t posClus[3];
    Double_t covClus[6];
    UShort_t volIdClus=GetClusterVolIndex(dig);
    GetClusterPars(dig, posClus,covClus);
    AliTOFcluster *tofCluster = new AliTOFcluster(volIdClus,posClus[0],posClus[1],posClus[2],covClus[0],covClus[1],covClus[2],covClus[3],covClus[4],covClus[5],d->GetTracks(),dig,parTOF,status,ii);
    InsertCluster(tofCluster);

  }

  AliDebug(1,Form("Number of found clusters: %d for event: %d", fNumberOfTofClusters, iEvent));

  CalibrateRecPoint();
  FillRecPoint();

  fTreeR->Fill();
  ResetRecpoint();

  fTOFLoader = fRunLoader->GetLoader("TOFLoader");  
  fTOFLoader->WriteRecPoints("OVERWRITE");

  AliDebug(1,Form("Execution time to read TOF digits and to write TOF clusters : R:%.4fs C:%.4fs",
		  stopwatch.RealTime(),stopwatch.CpuTime()));
  if (inholes) AliWarning(Form("Clusters in the TOF holes: %d",inholes));

}

//______________________________________________________________________________

void AliTOFClusterFinder::Digits2RecPoints(TTree* digitsTree, TTree* clusterTree)
{
  //
  // Converts digits to recpoints for TOF
  //

  TStopwatch stopwatch;
  stopwatch.Start();

  Int_t inholes = 0;

  if (digitsTree == 0x0) {
    AliFatal("AliTOFClusterFinder: Can not get TreeD");
    return;
  }

  fDigits->Clear();
  digitsTree->GetBranch("TOF")->SetAutoDelete(kFALSE);
  digitsTree->SetBranchAddress("TOF",&fDigits);

  ResetRecpoint();
  Int_t bufsize = 32000;
  clusterTree->Branch("TOF", &fRecPoints, bufsize);

  digitsTree->GetEvent(0);
  Int_t nDigits = fDigits->GetEntriesFast();
  AliDebug(2,Form("Number of TOF digits: %d",nDigits));

  Int_t ii;
  Int_t dig[5]={-1,-1,-1,-1,-1}; //cluster detector indeces
  Int_t parTOF[7]={0,0,0,0,0,0,0}; //The TOF signal parameters
  Bool_t status=kTRUE; // assume all sim channels ok in the beginning...
  for (ii=0; ii<nDigits; ii++) {
    AliTOFdigit *d = (AliTOFdigit*)fDigits->UncheckedAt(ii);
    dig[0]=d->GetSector();
    dig[1]=d->GetPlate();
    dig[2]=d->GetStrip();
    dig[3]=d->GetPadz();
    dig[4]=d->GetPadx();

    /* check valid index */
    if (dig[0]==-1||dig[1]==-1||dig[2]==-1||dig[3]==-1||dig[4]==-1) continue;

    // Do not reconstruct anything in the holes
    if (dig[0]==13 || dig[0]==14 || dig[0]==15 ) { // sectors with holes
      if (dig[1]==2) { // plate with holes
	inholes++;
	continue;
      }
    }

    //    AliDebug(2,Form(" %2d  %1d  %2d  %1d  %2d ",dig[0],dig[1],dig[2],dig[3],dig[4]));

    parTOF[0] = d->GetTdc(); //the TDC signal
    parTOF[1] = d->GetToT(); //the ToT signal
    parTOF[2] = d->GetAdc(); // the adc charge
    parTOF[3] = d->GetTdcND(); // non decalibrated sim time
    parTOF[4] = d->GetTdc(); // raw time, == Tdc time for the moment
    parTOF[5] = 0; // deltaBC
    parTOF[6] = 0; // L0-L1 latency
    
    Double_t posClus[3];
    Double_t covClus[6];
    UShort_t volIdClus=GetClusterVolIndex(dig);
    GetClusterPars(dig,posClus,covClus);
    AliTOFcluster *tofCluster = new AliTOFcluster(volIdClus,posClus[0],posClus[1],posClus[2],covClus[0],covClus[1],covClus[2],covClus[3],covClus[4],covClus[5],d->GetTracks(),dig,parTOF,status,ii);
    InsertCluster(tofCluster);

  }

  AliDebug(1,Form("Number of found clusters: %d", fNumberOfTofClusters));

  CalibrateRecPoint();
  FillRecPoint();

  clusterTree->Fill();
  ResetRecpoint();

  AliDebug(1,Form("Execution time to read TOF digits and to write TOF clusters : R:%.4fs C:%.4fs",
		  stopwatch.RealTime(),stopwatch.CpuTime()));
  if (inholes) AliWarning(Form("Clusters in the TOF holes: %d",inholes));

}
//______________________________________________________________________________

void AliTOFClusterFinder::Digits2RecPoints(AliRawReader *rawReader,
					   TTree *clustersTree)
{
  //
  // Converts RAW data to recpoints for TOF
  //

  TStopwatch stopwatch;
  stopwatch.Start();

  Int_t inholes = 0;

  const Int_t kDDL = AliDAQ::NumberOfDdls("TOF");

  ResetRecpoint();

  Int_t bufsize = 32000;
  clustersTree->Branch("TOF", &fRecPoints, bufsize);

  TClonesArray * clonesRawData;
  Int_t dummy = -1;

  Int_t detectorIndex[5];
  Int_t parTOF[7];

  ofstream ftxt;
  if (fVerbose==2) ftxt.open("TOFdigitsRead.txt",ios::app);

  fTOFRawStream.Clear();
  fTOFRawStream.SetRawReader(rawReader);

  if (fDecoderVersion == 1) {
    AliInfo("Using New Decoder");
  }
  else if (fDecoderVersion == 2) {
    AliInfo("Using Enhanced Decoder");
  }
  else {
    AliInfo("Using Old Decoder");
  }

  Int_t indexDDL = 0;
  for (indexDDL = 0; indexDDL < kDDL; indexDDL++) {

    rawReader->Reset();
    if (fDecoderVersion == 1) {
      fTOFRawStream.LoadRawDataBuffers(indexDDL,fVerbose);
    }
    else if (fDecoderVersion == 2) {
      fTOFRawStream.LoadRawDataBuffersV2(indexDDL,fVerbose);
    }
    else  {
      fTOFRawStream.LoadRawData(indexDDL);
    }
    
    clonesRawData = (TClonesArray*)fTOFRawStream.GetRawData();

    for (Int_t iRawData = 0; iRawData<clonesRawData->GetEntriesFast(); iRawData++) {

      AliTOFrawData *tofRawDatum = (AliTOFrawData*)clonesRawData->UncheckedAt(iRawData);

      //if (tofRawDatum->GetTOT()==-1 || tofRawDatum->GetTOF()==-1) continue;
      if (tofRawDatum->GetTOF()==-1) continue;

      if (fVerbose==2) {
	if (indexDDL<10) ftxt << "  " << indexDDL;
	else         ftxt << " " << indexDDL;
	if (tofRawDatum->GetTRM()<10) ftxt << "  " << tofRawDatum->GetTRM();
	else         ftxt << " " << tofRawDatum->GetTRM();
	ftxt << "  " << tofRawDatum->GetTRMchain();
	if (tofRawDatum->GetTDC()<10) ftxt << "  " << tofRawDatum->GetTDC();
	else         ftxt << " " << tofRawDatum->GetTDC();
	ftxt << "  " << tofRawDatum->GetTDCchannel();
      }

      fTOFRawStream.EquipmentId2VolumeId(indexDDL, tofRawDatum->GetTRM(), tofRawDatum->GetTRMchain(),
					 tofRawDatum->GetTDC(), tofRawDatum->GetTDCchannel(), detectorIndex);
      dummy = detectorIndex[3];
      detectorIndex[3] = detectorIndex[4];
      detectorIndex[4] = dummy;

      if (fVerbose==2) {
	if (detectorIndex[0]<10) ftxt  << "  ->  " << detectorIndex[0];
	else              ftxt  << "  -> " << detectorIndex[0];
	ftxt << "  " << detectorIndex[1];
	if (detectorIndex[2]<10) ftxt << "  " << detectorIndex[2];
	else              ftxt << " " << detectorIndex[2];
	ftxt << "  " << detectorIndex[3];
	if (detectorIndex[4]<10) ftxt << "  " << detectorIndex[4];
	else              ftxt << " " << detectorIndex[4];
      }

    /* check valid index */
    if (detectorIndex[0]==-1||detectorIndex[1]==-1||detectorIndex[2]==-1||detectorIndex[3]==-1||detectorIndex[4]==-1) continue;

      // Do not reconstruct anything in the holes
      if (detectorIndex[0]==13 || detectorIndex[0]==14 || detectorIndex[0]==15 ) { // sectors with holes
	if (detectorIndex[1]==2) { // plate with holes
	  inholes++;
	  continue;
	}
      }

      parTOF[0] = tofRawDatum->GetTOF(); //TDC
      parTOF[1] = tofRawDatum->GetTOT(); // TOT
      parTOF[2] = tofRawDatum->GetTOT(); //ADC==TOF
      parTOF[3] = 0;//raw data: no track of undecalib sim time
      parTOF[4] = tofRawDatum->GetTOF(); // RAW time
      parTOF[5] = tofRawDatum->GetDeltaBC(); // deltaBC
      parTOF[6] = tofRawDatum->GetL0L1Latency(); // L0-L1 latency
      Double_t posClus[3];
      Double_t covClus[6];
      UShort_t volIdClus=GetClusterVolIndex(detectorIndex);
      Int_t lab[3]={-1,-1,-1};
      Bool_t status=kTRUE;
      GetClusterPars(detectorIndex,posClus,covClus);
      AliTOFcluster *tofCluster = new AliTOFcluster(volIdClus,posClus[0],posClus[1],posClus[2],covClus[0],covClus[1],covClus[2],covClus[3],covClus[4],covClus[5],lab,detectorIndex,parTOF,status,-1);
      InsertCluster(tofCluster);

      if (fVerbose==2) {
	if (parTOF[1]<10)ftxt << "        " << parTOF[1];
	else if (parTOF[1]>=10 && parTOF[1]<100) ftxt << "      " << parTOF[1];
	else ftxt << "      " << parTOF[1];
	if (parTOF[0]<10) ftxt << "      " << parTOF[0] << endl;
	else if (parTOF[0]>=10 && parTOF[0]<100)   ftxt << "    " << parTOF[0] << endl;
	else if (parTOF[0]>=100 && parTOF[0]<1000) ftxt << "    " << parTOF[0] << endl;
	else ftxt << "   " << parTOF[0] << endl;
      }

    } // closed loop on TOF raw data per current DDL file

    clonesRawData->Clear("C");

  } // closed loop on DDL index

  if (fVerbose==2) ftxt.close();

  AliDebug(1,Form("Number of found clusters: %d", fNumberOfTofClusters));

  CalibrateRecPoint(rawReader->GetTimestamp());
  FillRecPoint();

  clustersTree->Fill();

  ResetRecpoint();

  AliDebug(1, Form("Execution time to read TOF raw data and to write TOF clusters : R:%.4fs C:%.4fs",
		   stopwatch.RealTime(),stopwatch.CpuTime()));
  if (inholes) AliWarning(Form("Clusters in the TOF holes: %d",inholes));

}
//______________________________________________________________________________

void AliTOFClusterFinder::Digits2RecPoints(Int_t iEvent, AliRawReader *rawReader)
{
  //
  // Converts RAW data to recpoints for TOF
  //

  TStopwatch stopwatch;
  stopwatch.Start();

  Int_t inholes = 0;

  const Int_t kDDL = AliDAQ::NumberOfDdls("TOF");

  fRunLoader->GetEvent(iEvent);

  AliDebug(2,Form(" Event number %2d ", iEvent));

  fTreeR = fTOFLoader->TreeR();

  if (fTreeR == 0x0){
    fTOFLoader->MakeTree("R");
    fTreeR = fTOFLoader->TreeR();
  }

  Int_t bufsize = 32000;
  fTreeR->Branch("TOF", &fRecPoints, bufsize);

  TClonesArray * clonesRawData;

  Int_t dummy = -1;

  Int_t detectorIndex[5] = {-1, -1, -1, -1, -1};
  Int_t parTOF[7];
  ofstream ftxt;
  if (fVerbose==2) ftxt.open("TOFdigitsRead.txt",ios::app);

  fTOFRawStream.Clear();
  fTOFRawStream.SetRawReader(rawReader);

  if (fDecoderVersion == 1) {
    AliInfo("Using New Decoder");
  }
  else if (fDecoderVersion == 2) {
    AliInfo("Using Enhanced Decoder");
  }
  else {
    AliInfo("Using Old Decoder");
  }

  Int_t indexDDL = 0;
  for (indexDDL = 0; indexDDL < kDDL; indexDDL++) {

    rawReader->Reset();
    if (fDecoderVersion == 1) {
      fTOFRawStream.LoadRawDataBuffers(indexDDL,fVerbose);
    }
    else if (fDecoderVersion == 2) {
      fTOFRawStream.LoadRawDataBuffersV2(indexDDL,fVerbose);
    }
    else {
      fTOFRawStream.LoadRawData(indexDDL);
    }

    clonesRawData = (TClonesArray*)fTOFRawStream.GetRawData();

    for (Int_t iRawData = 0; iRawData<clonesRawData->GetEntriesFast(); iRawData++) {

      AliTOFrawData *tofRawDatum = (AliTOFrawData*)clonesRawData->UncheckedAt(iRawData);

      //if (tofRawDatum->GetTOT()==-1 || tofRawDatum->GetTOF()==-1) continue;
      if (tofRawDatum->GetTOF()==-1) continue;

      if (fVerbose==2) {
	if (indexDDL<10) ftxt << "  " << indexDDL;
	else         ftxt << " " << indexDDL;
	if (tofRawDatum->GetTRM()<10) ftxt << "  " << tofRawDatum->GetTRM();
	else         ftxt << " " << tofRawDatum->GetTRM();
	ftxt << "  " << tofRawDatum->GetTRMchain();
	if (tofRawDatum->GetTDC()<10) ftxt << "  " << tofRawDatum->GetTDC();
	else         ftxt << " " << tofRawDatum->GetTDC();
	ftxt << "  " << tofRawDatum->GetTDCchannel();
      }

      fTOFRawStream.EquipmentId2VolumeId(indexDDL, tofRawDatum->GetTRM(), tofRawDatum->GetTRMchain(),
					 tofRawDatum->GetTDC(), tofRawDatum->GetTDCchannel(), detectorIndex);
      dummy = detectorIndex[3];
      detectorIndex[3] = detectorIndex[4];
      detectorIndex[4] = dummy;

      if (fVerbose==2) {
	if (detectorIndex[0]<10) ftxt  << "  ->  " << detectorIndex[0];
	else              ftxt  << "  -> " << detectorIndex[0];
	ftxt << "  " << detectorIndex[1];
	if (detectorIndex[2]<10) ftxt << "  " << detectorIndex[2];
	else              ftxt << " " << detectorIndex[2];
	ftxt << "  " << detectorIndex[3];
	if (detectorIndex[4]<10) ftxt << "  " << detectorIndex[4];
	else              ftxt << " " << detectorIndex[4];
      }

    /* check valid index */
    if (detectorIndex[0]==-1||detectorIndex[1]==-1||detectorIndex[2]==-1||detectorIndex[3]==-1||detectorIndex[4]==-1) continue;

      // Do not reconstruct anything in the holes
      if (detectorIndex[0]==13 || detectorIndex[0]==14 || detectorIndex[0]==15 ) { // sectors with holes
	if (detectorIndex[1]==2) { // plate with holes
	  inholes++;
	  continue;
	}
      }

      parTOF[0] = tofRawDatum->GetTOF(); // TDC
      parTOF[1] = tofRawDatum->GetTOT(); // TOT
      parTOF[2] = tofRawDatum->GetTOT(); // raw data have ADC=TOT
      parTOF[3] = 0; //raw data: no track of the undecalib sim time
      parTOF[4] = tofRawDatum->GetTOF(); // Raw time == TDC
      parTOF[5] = tofRawDatum->GetDeltaBC(); // deltaBC
      parTOF[6] = tofRawDatum->GetL0L1Latency(); // L0-L1 latency
      Double_t posClus[3];
      Double_t covClus[6];
      UShort_t volIdClus=GetClusterVolIndex(detectorIndex);
      Int_t lab[3]={-1,-1,-1};
      Bool_t status=kTRUE;
      GetClusterPars(detectorIndex,posClus,covClus);
      AliTOFcluster *tofCluster = new AliTOFcluster(volIdClus,posClus[0],posClus[1],posClus[2],covClus[0],covClus[1],covClus[2],covClus[3],covClus[4],covClus[5],lab,detectorIndex,parTOF,status,-1);
      InsertCluster(tofCluster);

      if (fVerbose==2) {
	if (parTOF[1]<10)ftxt << "        " << parTOF[1];
	else if (parTOF[1]>=10 && parTOF[1]<100) ftxt << "      " << parTOF[1];
	else ftxt << "      " << parTOF[1];
	if (parTOF[0]<10) ftxt << "      " << parTOF[0] << endl;
	else if (parTOF[0]>=10 && parTOF[0]<100)   ftxt << "    " << parTOF[0] << endl;
	else if (parTOF[0]>=100 && parTOF[0]<1000) ftxt << "    " << parTOF[0] << endl;
	else ftxt << "   " << parTOF[0] << endl;
      }

    } // closed loop on TOF raw data per current DDL file

    clonesRawData->Clear("C");

  } // closed loop on DDL index

  if (fVerbose==2) ftxt.close();

  AliDebug(1,Form("Number of found clusters: %d for event: %d", fNumberOfTofClusters, iEvent));

  CalibrateRecPoint(rawReader->GetTimestamp());
  FillRecPoint();

  fTreeR->Fill();
  ResetRecpoint();

  fTOFLoader = fRunLoader->GetLoader("TOFLoader");
  fTOFLoader->WriteRecPoints("OVERWRITE");
  
  AliDebug(1, Form("Execution time to read TOF raw data and to write TOF clusters : R:%.4fs C:%.4fs",
	       stopwatch.RealTime(),stopwatch.CpuTime()));
  if (inholes) AliWarning(Form("Clusters in the TOF holes: %d",inholes));

}
//______________________________________________________________________________

void AliTOFClusterFinder::Raw2Digits(Int_t iEvent, AliRawReader *rawReader)
{
  //
  // Converts RAW data to MC digits for TOF
  //
  //             (temporary solution)
  //

  TStopwatch stopwatch;
  stopwatch.Start();

  const Int_t kDDL = AliTOFGeometry::NDDL()*AliTOFGeometry::NSectors();

  fRunLoader->GetEvent(iEvent);

  fTreeD = fTOFLoader->TreeD();
  if (fTreeD)
    {
    AliInfo("TreeD re-creation");
    fTreeD = 0x0;
    fTOFLoader->MakeTree("D");
    fTreeD = fTOFLoader->TreeD();
    }
  else {
    AliFatal("Can not get TreeD");
    return;
  }

  Int_t bufsize = 32000;
  fDigits->Clear();
  fTreeD->Branch("TOF", &fDigits, bufsize);

  fRunLoader->GetEvent(iEvent);

  AliDebug(2,Form(" Event number %2d ", iEvent));

  TClonesArray * clonesRawData;

  Int_t dummy = -1;

  Int_t detectorIndex[5];
  Int_t digit[4];

  fTOFRawStream.Clear();
  fTOFRawStream.SetRawReader(rawReader);

  if (fDecoderVersion == 1) {
    AliInfo("Using New Decoder");
  }
  else if (fDecoderVersion == 2) {
    AliInfo("Using Enhanced Decoder");
  }
  else {
    AliInfo("Using Old Decoder");
  }

  Int_t indexDDL = 0;
  for (indexDDL = 0; indexDDL < kDDL; indexDDL++) {

    rawReader->Reset();
    if (fDecoderVersion == 1) {
      fTOFRawStream.LoadRawDataBuffers(indexDDL,fVerbose);
    }
    else if (fDecoderVersion == 2) {
      fTOFRawStream.LoadRawDataBuffersV2(indexDDL,fVerbose);
    }
    else {
      fTOFRawStream.LoadRawData(indexDDL);
    }

    clonesRawData = (TClonesArray*)fTOFRawStream.GetRawData();

    for (Int_t iRawData = 0; iRawData<clonesRawData->GetEntriesFast(); iRawData++) {

      AliTOFrawData *tofRawDatum = (AliTOFrawData*)clonesRawData->UncheckedAt(iRawData);

      //if (!tofRawDatum->GetTOT() || !tofRawDatum->GetTOF()) continue;
      if (tofRawDatum->GetTOF()==-1) continue;

      fTOFRawStream.EquipmentId2VolumeId(indexDDL, tofRawDatum->GetTRM(), tofRawDatum->GetTRMchain(),
					 tofRawDatum->GetTDC(), tofRawDatum->GetTDCchannel(), detectorIndex);
      dummy = detectorIndex[3];
      detectorIndex[3] = detectorIndex[4];
      detectorIndex[4] = dummy;

      digit[0] = fTOFRawStream.GetTofBin();
      digit[1] = fTOFRawStream.GetToTbin();
      digit[2] = fTOFRawStream.GetToTbin();
      digit[3] = 0;

      Int_t tracknum[3]={-1,-1,-1};

      TClonesArray &aDigits = *fDigits;
      Int_t last=fDigits->GetEntriesFast();
      new (aDigits[last]) AliTOFdigit(tracknum, detectorIndex, digit);

    } // while loop

    clonesRawData->Clear("C");

  } // DDL Loop

  fTreeD->Fill();

  fTOFLoader = fRunLoader->GetLoader("TOFLoader");
  fTOFLoader->WriteDigits("OVERWRITE");

  AliDebug(1, Form("Execution time to read TOF raw data and to write TOF clusters : R:%.2fs C:%.2fs",
		   stopwatch.RealTime(),stopwatch.CpuTime()));

}

//______________________________________________________________________________

void AliTOFClusterFinder::Raw2Digits(AliRawReader *rawReader, TTree* digitsTree)
{
  //
  // Converts RAW data to MC digits for TOF for the current event
  //
  //

  TStopwatch stopwatch;
  stopwatch.Start();

  const Int_t kDDL = AliTOFGeometry::NDDL()*AliTOFGeometry::NSectors();

  if (!digitsTree)
    {
    AliError("No input digits Tree");
    return;
    }

  Int_t bufsize = 32000;
  digitsTree->Branch("TOF", &fDigits, bufsize);

  TClonesArray * clonesRawData;
  Int_t dummy = -1;

  Int_t detectorIndex[5];
  Int_t digit[4];

  fTOFRawStream.Clear();
  fTOFRawStream.SetRawReader(rawReader);

  if (fDecoderVersion == 1) {
    AliInfo("Using New Decoder");
  }
  else if (fDecoderVersion == 2) {
    AliInfo("Using Enhanced Decoder");
  }
  else {
    AliInfo("Using Old Decoder");
  }

  Int_t indexDDL = 0;
  for (indexDDL = 0; indexDDL < kDDL; indexDDL++) {

    rawReader->Reset();
    if (fDecoderVersion == 1) {
      fTOFRawStream.LoadRawDataBuffers(indexDDL,fVerbose);
    }
    else if (fDecoderVersion == 2) {
      fTOFRawStream.LoadRawDataBuffersV2(indexDDL,fVerbose);
    }
    else {
      fTOFRawStream.LoadRawData(indexDDL);
    }

    clonesRawData = (TClonesArray*)fTOFRawStream.GetRawData();

    for (Int_t iRawData = 0; iRawData<clonesRawData->GetEntriesFast(); iRawData++) {

      AliTOFrawData *tofRawDatum = (AliTOFrawData*)clonesRawData->UncheckedAt(iRawData);

      //if (!tofRawDatum->GetTOT() || !tofRawDatum->GetTOF()) continue;
      if (tofRawDatum->GetTOF()==-1) continue;

      fTOFRawStream.EquipmentId2VolumeId(indexDDL, tofRawDatum->GetTRM(), tofRawDatum->GetTRMchain(),
					 tofRawDatum->GetTDC(), tofRawDatum->GetTDCchannel(), detectorIndex);
      dummy = detectorIndex[3];
      detectorIndex[3] = detectorIndex[4];
      detectorIndex[4] = dummy;

      digit[0] = fTOFRawStream.GetTofBin();
      digit[1] = fTOFRawStream.GetToTbin();
      digit[2] = fTOFRawStream.GetToTbin();
      digit[3] = 0;

      Int_t tracknum[3]={-1,-1,-1};

      TClonesArray &aDigits = *fDigits;
      Int_t last=fDigits->GetEntriesFast();
      new (aDigits[last]) AliTOFdigit(tracknum, detectorIndex, digit);

    } // while loop

    clonesRawData->Clear("C");

  } // DDL Loop

  digitsTree->Fill();

  AliDebug(1, Form("Got %d digits: ", fDigits->GetEntries()));
  AliDebug(1, Form("Execution time to read TOF raw data and fill TOF digit tree : R:%.2fs C:%.2fs",
		   stopwatch.RealTime(),stopwatch.CpuTime()));

}
//______________________________________________________________________________

Int_t AliTOFClusterFinder::InsertCluster(AliTOFcluster *tofCluster) {
  //---------------------------------------------------------------------------//
  // This function adds a TOF cluster to the array of TOF clusters sorted in Z //
  //---------------------------------------------------------------------------//
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
//_________________________________________________________________________

Int_t AliTOFClusterFinder::FindClusterIndex(Double_t z) const {
  //--------------------------------------------------------------------
  // This function returns the index of the nearest cluster in z
  //--------------------------------------------------------------------
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
//_________________________________________________________________________

void AliTOFClusterFinder::FillRecPoint()
{
  //
  // Copy the global array of AliTOFcluster, i.e. fTofClusters (sorted
  // in Z) in the global TClonesArray of AliTOFcluster,
  // i.e. fRecPoints.
  //

  Int_t ii, jj;

  Int_t detectorIndex[5];
  Int_t parTOF[7];
  Int_t trackLabels[3];
  Int_t digitIndex = -1;
  Bool_t status=kTRUE;

  TClonesArray &lRecPoints = *fRecPoints;
  
  for (ii=0; ii<fNumberOfTofClusters; ii++) {

    digitIndex = fTofClusters[ii]->GetIndex();
    for(jj=0; jj<5; jj++) detectorIndex[jj] = fTofClusters[ii]->GetDetInd(jj);
    for(jj=0; jj<3; jj++) trackLabels[jj] = fTofClusters[ii]->GetLabel(jj);
    parTOF[0] = fTofClusters[ii]->GetTDC(); // TDC
    parTOF[1] = fTofClusters[ii]->GetToT(); // TOT
    parTOF[2] = fTofClusters[ii]->GetADC(); // ADC=TOT
    parTOF[3] = fTofClusters[ii]->GetTDCND(); // TDCND
    parTOF[4] = fTofClusters[ii]->GetTDCRAW();//RAW
    parTOF[5] = fTofClusters[ii]->GetDeltaBC();//deltaBC
    parTOF[6] = fTofClusters[ii]->GetL0L1Latency();//L0-L1 latency
    status=fTofClusters[ii]->GetStatus();
    Double_t posClus[3];
    Double_t covClus[6];
    UShort_t volIdClus=GetClusterVolIndex(detectorIndex);
    GetClusterPars(detectorIndex,posClus,covClus);
    new(lRecPoints[ii]) AliTOFcluster(volIdClus,posClus[0],posClus[1],posClus[2],covClus[0],covClus[1],covClus[2],covClus[3],covClus[4],covClus[5],trackLabels,detectorIndex, parTOF,status,digitIndex);

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

  } // loop on clusters

}

//_________________________________________________________________________

/*
 * OLD CALIBRATE REC POINTS FUNCTION
 */

#if 0
void AliTOFClusterFinder::CalibrateRecPoint(UInt_t timestamp)
{
  //
  // Copy the global array of AliTOFcluster, i.e. fTofClusters (sorted
  // in Z) in the global TClonesArray of AliTOFcluster,
  // i.e. fRecPoints.
  //

  Int_t ii, jj;

  Int_t detectorIndex[5];
  Int_t digitIndex = -1;
  Double_t tToT;
  Double_t timeCorr;
  Int_t   tdcCorr;
  Float_t tdcLatencyWindow;
  AliDebug(1," Calibrating TOF Clusters");

  AliTOFChannelOnlineArray *calDelay = fTOFcalib->GetTOFOnlineDelay();  
  AliTOFChannelOnlineStatusArray *calStatus = fTOFcalib->GetTOFOnlineStatus();  
  TObjArray *calTOFArrayOffline = fTOFcalib->GetTOFCalArrayOffline();
  
  AliTOFDeltaBCOffset *deltaBCOffsetObj = fTOFcalib->GetDeltaBCOffset();
  Int_t deltaBCOffset = deltaBCOffsetObj->GetDeltaBCOffset();
  AliTOFCTPLatency *ctpLatencyObj = fTOFcalib->GetCTPLatency();
  Float_t ctpLatency = ctpLatencyObj->GetCTPLatency();
  AliTOFRunParams *runParamsObj = fTOFcalib->GetRunParams();
  Float_t t0 = runParamsObj->EvalT0(timestamp);

  TString validity = (TString)fTOFcalib->GetOfflineValidity();
  Int_t calibration = -1;
  if (validity.CompareTo("valid")==0) {
    //AliInfo(Form(" validity = %s - Using offline calibration parameters",validity.Data()));
    calibration = 1;
  } else {
    //AliInfo(Form(" validity = %s - Using online calibration parameters",validity.Data()));
    calibration = 0 ;
  }

  for (ii=0; ii<fNumberOfTofClusters; ii++) {
    digitIndex = fTofClusters[ii]->GetIndex();
    for(jj=0; jj<5; jj++) detectorIndex[jj] = fTofClusters[ii]->GetDetInd(jj);

    Int_t index = AliTOFGeometry::GetIndex(detectorIndex);
     
    UChar_t statusPulser=calStatus->GetPulserStatus(index);
    UChar_t statusNoise=calStatus->GetNoiseStatus(index);
    UChar_t statusHW=calStatus->GetHWStatus(index);
    UChar_t status=calStatus->GetStatus(index);
    tdcLatencyWindow = calStatus->GetLatencyWindow(index) * 1.e3; /* ns -> ps */
    
    //check the status, also unknown is fine!!!!!!!

    AliDebug(2, Form(" Status for channel %d = %d",index, (Int_t)status));
    if((statusPulser & AliTOFChannelOnlineStatusArray::kTOFPulserBad)==(AliTOFChannelOnlineStatusArray::kTOFPulserBad)||(statusNoise & AliTOFChannelOnlineStatusArray::kTOFNoiseBad)==(AliTOFChannelOnlineStatusArray::kTOFNoiseBad)||(statusHW & AliTOFChannelOnlineStatusArray::kTOFHWBad)==(AliTOFChannelOnlineStatusArray::kTOFHWBad)){
      AliDebug(2, Form(" Bad Status for channel %d",index));
      fTofClusters[ii]->SetStatus(kFALSE); //odd convention, to avoid conflict with calibration objects currently in the db (temporary solution).
    }
    else {
      AliDebug(2, Form(" Good Status for channel %d",index));
    }
    // Get Rough channel online equalization 
    Double_t roughDelay=(Double_t)calDelay->GetDelay(index);  // in ns
    AliDebug(2,Form(" channel delay (ns) = %f", roughDelay));
    // Get Refined channel offline calibration parameters
    if (calibration ==1){
      AliTOFChannelOffline * calChannelOffline = (AliTOFChannelOffline*)calTOFArrayOffline->At(index);
      Double_t par[6];
      for (Int_t j = 0; j<6; j++){
	par[j]=(Double_t)calChannelOffline->GetSlewPar(j);
     } 
      AliDebug(2,Form(" Calib Pars = %f, %f, %f, %f, %f, %f ",par[0],par[1],par[2],par[3],par[4],par[5]));
      AliDebug(2,Form(" The ToT and Time, uncorr (counts) = %d , %d", fTofClusters[ii]->GetToT(),fTofClusters[ii]->GetTDC()));
      tToT = (Double_t)(fTofClusters[ii]->GetToT()*AliTOFGeometry::ToTBinWidth());
      tToT*=1.E-3; //ToT in ns

      /* check TOT limits and set new TOT in case */
      if (tToT < AliTOFGeometry::SlewTOTMin()) tToT = AliTOFGeometry::SlewTOTMin();
      if (tToT > AliTOFGeometry::SlewTOTMax()) tToT = AliTOFGeometry::SlewTOTMax();

      AliDebug(2,Form(" The ToT and Time, uncorr (ns)= %e, %e",fTofClusters[ii]->GetTDC()*AliTOFGeometry::TdcBinWidth()*1.E-3,tToT));
      timeCorr=par[0]+par[1]*tToT+par[2]*tToT*tToT+par[3]*tToT*tToT*tToT+par[4]*tToT*tToT*tToT*tToT+par[5]*tToT*tToT*tToT*tToT*tToT; // the time correction (ns)
    }
    else {
      timeCorr = roughDelay; // correction in ns
    }
    AliDebug(2,Form(" The ToT and Time, uncorr (ns)= %e, %e",fTofClusters[ii]->GetTDC()*AliTOFGeometry::TdcBinWidth()*1.E-3,fTofClusters[ii]->GetToT()*AliTOFGeometry::ToTBinWidth()));
    AliDebug(2,Form(" The time correction (ns) = %f", timeCorr));
    timeCorr=(Double_t)(fTofClusters[ii]->GetTDC())*AliTOFGeometry::TdcBinWidth()*1.E-3-timeCorr;//redefine the time
    timeCorr*=1.E3;
    AliDebug(2,Form(" The channel time, corr (ps)= %e",timeCorr ));

    /* here timeCorr should be already corrected for calibration. 
     * we now go into further corrections keeping in mind that timeCorr
     * is in ps.
     *
     * the following corrections are performed in this way:
     *
     *    time = timeRaw - deltaBC + L0L1Latency + CTPLatency - TDCLatencyWindow - T0
     *
     */

    AliDebug(2, Form("applying further corrections (DeltaBC): DeltaBC=%d (BC bins), DeltaBCoffset=%d (BC bins)", fTofClusters[ii]->GetDeltaBC(), deltaBCOffset));
    AliDebug(2, Form("applying further corrections (L0L1Latency): L0L1Latency=%d (BC bins)", fTofClusters[ii]->GetL0L1Latency()));
    AliDebug(2, Form("applying further corrections (CTPLatency): CTPLatency=%f (ps)", ctpLatency));
    AliDebug(2, Form("applying further corrections (TDCLatencyWindow): TDCLatencyWindow=%f (ps)", tdcLatencyWindow));
    AliDebug(2, Form("applying further corrections (T0): T0=%f (ps)", t0));

    /* deltaBC correction (inhibited for the time being) */
    //    timeCorr -= (fTofClusters[ii]->GetDeltaBC() - deltaBCOffset) * AliTOFGeometry::BunchCrossingBinWidth();
    /* L0L1-latency correction */
    timeCorr += fTofClusters[ii]->GetL0L1Latency() * AliTOFGeometry::BunchCrossingBinWidth();
    /* CTP-latency correction (from OCDB) */
    timeCorr += ctpLatency;
    /* TDC latency-window correction (from OCDB) */
    timeCorr -= tdcLatencyWindow;
    /* T0 correction (from OCDB) */
    timeCorr -= t0;

    /*
     * end of furhter corrections
     */

    tdcCorr=(Int_t)(timeCorr/AliTOFGeometry::TdcBinWidth()); //the corrected time (tdc counts)
    fTofClusters[ii]->SetTDC(tdcCorr);
  } // loop on clusters

}
#endif

//_________________________________________________________________________

void AliTOFClusterFinder::CalibrateRecPoint(UInt_t timestamp)
{
  //
  // Copy the global array of AliTOFcluster, i.e. fTofClusters (sorted
  // in Z) in the global TClonesArray of AliTOFcluster,
  // i.e. fRecPoints.
  //

  Int_t detectorIndex[5];
  Double_t time, tot, corr;
  Int_t deltaBC, l0l1, tdcBin;
  for (Int_t ii = 0; ii < fNumberOfTofClusters; ii++) {
    for(Int_t jj = 0; jj < 5; jj++) detectorIndex[jj] = fTofClusters[ii]->GetDetInd(jj);

    Int_t index = AliTOFGeometry::GetIndex(detectorIndex);

    /* check channel enabled */
    if (!fTOFcalib->IsChannelEnabled(index)) fTofClusters[ii]->SetStatus(kFALSE);
    
    /* get cluster info */
    time = fTofClusters[ii]->GetTDC() * AliTOFGeometry::TdcBinWidth(); /* ps */
    tot = fTofClusters[ii]->GetToT() * AliTOFGeometry::ToTBinWidth() * 1.e-3; /* ns */
    deltaBC = fTofClusters[ii]->GetDeltaBC();
    l0l1 = fTofClusters[ii]->GetL0L1Latency();

    /* get correction */
    corr = fTOFcalib->GetTimeCorrection(index, tot, deltaBC, l0l1, timestamp); /* ps */
    AliDebug(2, Form("calibrate index %d: time=%f (ps) tot=%f (ns) deltaBC=%d l0l1=%d timestamp=%d corr=%f (ps)", index, time, tot, deltaBC, l0l1, timestamp, corr));

    /* apply time correction */
    time -= corr;

    /* convert in TDC bins and set cluster */
    //tdcBin = (Int_t)(time / AliTOFGeometry::TdcBinWidth()); //the corrected time (tdc counts)
    tdcBin = TMath::Nint(time / AliTOFGeometry::TdcBinWidth()); //the corrected time (tdc counts)
    fTofClusters[ii]->SetTDC(tdcBin);

  } // loop on clusters

}

//______________________________________________________________________________

void AliTOFClusterFinder::ResetRecpoint()
{
  //
  // Clear the list of reconstructed points
  //

  fNumberOfTofClusters = 0;
  if (fRecPoints) fRecPoints->Clear();

}
//______________________________________________________________________________

void AliTOFClusterFinder::Load()
{
  //
  // Load TOF.Digits.root and TOF.RecPoints.root files
  //

  fTOFLoader->LoadDigits("READ");
  fTOFLoader->LoadRecPoints("recreate");

}
//______________________________________________________________________________

void AliTOFClusterFinder::LoadClusters()
{
  //
  // Load TOF.RecPoints.root file
  //

  fTOFLoader->LoadRecPoints("recreate");

}
//______________________________________________________________________________

void AliTOFClusterFinder::UnLoad()
{
  //
  // Unload TOF.Digits.root and TOF.RecPoints.root files
  //

  fTOFLoader->UnloadDigits();
  fTOFLoader->UnloadRecPoints();

}
//______________________________________________________________________________

void AliTOFClusterFinder::UnLoadClusters()
{
  //
  // Unload TOF.RecPoints.root file
  //

  fTOFLoader->UnloadRecPoints();

}
//-------------------------------------------------------------------------
UShort_t AliTOFClusterFinder::GetClusterVolIndex(const Int_t * const ind) const {

  //First of all get the volume ID to retrieve the l2t transformation...
  //
  // Detector numbering scheme
  Int_t nSector = 18;
  Int_t nPlate  = 5;
  Int_t nStripA = 15;
  Int_t nStripB = 19;
  Int_t nStripC = 19;

  Int_t isector =ind[0];
  if (isector >= nSector)
    AliError(Form("Wrong sector number in TOF (%d) !",isector));
  Int_t iplate = ind[1];
  if (iplate >= nPlate)
    AliError(Form("Wrong plate number in TOF (%d) !",iplate));
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
    AliError(Form("Wrong plate number in TOF (%d) !",iplate));
    break;
  };

  Int_t index= (2*(nStripC+nStripB)+nStripA)*isector +
               stripOffset +
               istrip;

  UShort_t volIndex = AliGeomManager::LayerToVolUID(AliGeomManager::kTOF,index);
  return volIndex;
}
//
//-------------------------------------------------------------------------
void AliTOFClusterFinder::GetClusterPars(Int_t *ind, Double_t* pos,Double_t* cov) const {

  //First of all get the volume ID to retrieve the l2t transformation...
  //
  UShort_t volIndex = GetClusterVolIndex(ind);
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

  /*
  Float_t localX = (ind[4]-23.5)*AliTOFGeometry::XPad();
  Float_t localY = 0;
  Float_t localZ = (ind[3]- 0.5)*AliTOFGeometry::ZPad();
  */
  Float_t localX = (ind[4]-AliTOFGeometry::NpadX()/2)*AliTOFGeometry::XPad()+AliTOFGeometry::XPad()/2.;
  Float_t localY = 0;
  Float_t localZ = (ind[3]-AliTOFGeometry::NpadZ()/2)*AliTOFGeometry::ZPad()+AliTOFGeometry::ZPad()/2.;
  //move to the tracking ref system

  Double_t lpos[3];
  lpos[0] = localX;
  lpos[1] = localY;
  lpos[2] = localZ; 

  const TGeoHMatrix *l2t= AliGeomManager::GetTracking2LocalMatrix(volIndex);
  // Get The position in the track ref system
  Double_t tpos[3];
  l2t->MasterToLocal(lpos,tpos);
  pos[0] = tpos[0];
  pos[1] = tpos[1];
  pos[2] = tpos[2];

  //Get The cluster covariance in the track ref system
  Double_t lcov[9];

  //cluster covariance in the local system:
  // sx2   0   0
  // 0     0   0
  // 0     0   sz2

  lcov[0] = AliTOFGeometry::XPad()*AliTOFGeometry::XPad()/12.;
  lcov[1] = 0;
  lcov[2] = 0;
  lcov[3] = 0;
  lcov[4] = 0;
  lcov[5] = 0;
  lcov[6] = 0;
  lcov[7] = 0;
  lcov[8] = AliTOFGeometry::ZPad()*AliTOFGeometry::ZPad()/12.;

  //cluster covariance in the tracking system:
  TGeoHMatrix m;
  m.SetRotation(lcov);
  m.Multiply(l2t);
  m.MultiplyLeft(&l2t->Inverse());
  Double_t *tcov = m.GetRotationMatrix();
  cov[0] = tcov[0]; cov[1] = tcov[1]; cov[2] = tcov[2];
  cov[3] = tcov[4]; cov[4] = tcov[5];
  cov[5] = tcov[8];

  return;

}
