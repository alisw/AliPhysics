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

/*
$Log: AliTOFcalib.cxx,v $
Revision 1.21  2007/11/02 15:41:49  hristov
Provide return value if the function is not void

Revision 1.20  2007/10/26 15:13:50  zampolli
Using a TChain instead of a TTree

Revision 1.19  2007/10/23 15:27:38  zampolli
Rearrangement of Calibration objects for simulation

Revision 1.16  2007/10/08 10:13:26  zampolli
First Run and Last Run members added, infinite validity of calib obj implemented.

Revision 1.15  2007/10/04 13:23:28  zampolli
Updates to handle functionalities in TOF online/offline calibration according to the latest schema

Revision 1.14  2007/06/06 16:26:30  arcelli
remove fall-back call to local CDB storage

Revision 1.13  2007/04/20 13:59:40  arcelli
make protections agains failed retrieval of the CDB object in a proper way

Revision 1.12  2007/03/23 11:31:16  arcelli
CDB Entry for TOF Reconstruction Parameters

Revision 1.11  2007/02/28 18:08:26  arcelli
Add protection against failed retrieval of the CDB cal object

Revision 1.10  2006/08/22 13:30:49  arcelli
removal of effective c++ warnings (C.Zampolli)

Revision 1.9  2006/04/20 22:30:50  hristov
Coding conventions (Annalisa)

Revision 1.8  2006/04/16 22:29:05  hristov
Coding conventions (Annalisa)

Revision 1.7  2006/04/16 20:12:46  hristov
Removing memory leak in case of cached CDB entries

Revision 1.6  2006/04/11 15:28:32  hristov
Checks on cache status before deleting calibration objects (A.Colla)

Revision 1.5  2006/04/05 08:35:38  hristov
Coding conventions (S.Arcelli, C.Zampolli)

Revision 1.4  2006/03/31 11:26:46  arcelli
 changing CDB Ids according to standard convention

Revision 1.3  2006/03/28 14:57:02  arcelli
updates to handle new V5 geometry & some re-arrangements

Revision 1.2  2006/02/13 17:22:26  arcelli
just Fixing Log info

Revision 1.1  2006/02/13 16:10:48  arcelli
Add classes for TOF Calibration (C.Zampolli)

author: Chiara Zampolli, zampolli@bo.infn.it
*/  

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// class for TOF calibration                                                 //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "TF1.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH1C.h"
#include "TH2F.h"
//#include "TList.h"
//#include "TROOT.h"
//#include "TStyle.h"
#include "TTree.h"
#include "TChain.h"
#include "TProfile.h"
#include "TGrid.h"
#include "TMath.h"
#include "TMap.h"

#include "AliCDBEntry.h"
#include "AliCDBRunRange.h"
#include "AliCDBId.h"
#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#include "AliCDBMetaData.h"
//#include "AliESDtrack.h"
//#include "AliESD.h"
#include "AliLog.h"

#include "AliTOFcalib.h"
#include "AliTOFChannelOnlineArray.h"
#include "AliTOFChannelOnline.h"
#include "AliTOFChannelOnlineStatus.h"
#include "AliTOFChannelOnlineStatusArray.h"
#include "AliTOFChannelOffline.h"
#include "AliTOFGeometry.h"
#include "AliTOFRecoParam.h"
#include "AliTOFDeltaBCOffset.h"
#include "AliTOFCTPLatency.h"
#include "AliTOFT0Fill.h"
#include "AliTOFRunParams.h"
#include "AliLHCClockPhase.h"
#include "AliTOFResponseParams.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "TRandom.h"

class TROOT;
class TStyle;

extern TROOT *gROOT;
extern TStyle *gStyle;

ClassImp(AliTOFcalib)

//_______________________________________________________________________
AliTOFcalib::AliTOFcalib():
  TTask("AliTOFcalib",""),
  fNChannels(-1),
  fTOFCalOnline(0x0),
  fTOFCalOnlinePulser(0x0),
  fTOFCalOnlineNoise(0x0),
  fTOFCalOnlineHW(0x0),
  fTOFCalOffline(0x0),
  fCal(0x0),
  fStatus(0x0),
  fTOFSimToT(0x0),
  fkValidity(0x0),
  fTree(0x0),
  fChain(0x0),
  fNruns(0),
  fFirstRun(0),
  fLastRun(AliCDBRunRange::Infinity()),
  fConfigMap(new TMap),
  fDeltaBCOffset(NULL),
  fCTPLatency(NULL),
  fT0Fill(NULL),
  fRunParams(NULL),
  fLHCClockPhase(NULL),
  fResponseParams(NULL),
  fReadoutEfficiency(NULL),
  fProblematic(NULL),
  fInitFlag(kFALSE),
  fRemoveMeanT0(kTRUE),
  fUseLHCClockPhase(kFALSE),
  fCalibrateTOFsignal(kTRUE),
  fCorrectTExp(kFALSE)
{ 
  //TOF Calibration Class ctor
  fNChannels = AliTOFGeometry::NSectors()*(2*(AliTOFGeometry::NStripC()+AliTOFGeometry::NStripB())+AliTOFGeometry::NStripA())*AliTOFGeometry::NpadZ()*AliTOFGeometry::NpadX();

  gRandom->SetSeed(123456789);
}
//____________________________________________________________________________ 

AliTOFcalib::AliTOFcalib(const AliTOFcalib & calib):
  TTask(calib),
  fNChannels(calib.fNChannels),
  fTOFCalOnline(0x0),
  fTOFCalOnlinePulser(0x0),
  fTOFCalOnlineNoise(0x0),
  fTOFCalOnlineHW(0x0),
  fTOFCalOffline(0x0),
  fCal(calib.fCal),
  fStatus(calib.fStatus),
  fTOFSimToT(calib.fTOFSimToT),
  fkValidity(calib.fkValidity),
  fTree(calib.fTree),
  fChain(calib.fChain),
  fNruns(calib.fNruns),
  fFirstRun(calib.fFirstRun),
  fLastRun(calib.fLastRun),
  fConfigMap(calib.fConfigMap),
  fDeltaBCOffset(NULL),
  fCTPLatency(NULL),
  fT0Fill(NULL),
  fRunParams(NULL),
  fLHCClockPhase(NULL),
  fResponseParams(NULL),
  fReadoutEfficiency(NULL),
  fProblematic(NULL),
  fInitFlag(calib.fInitFlag),
  fRemoveMeanT0(calib.fRemoveMeanT0),
  fUseLHCClockPhase(calib.fUseLHCClockPhase),
  fCalibrateTOFsignal(calib.fCalibrateTOFsignal),
  fCorrectTExp(calib.fCorrectTExp)
{

  fTOFCalOnline = new TObjArray(fNChannels);
  fTOFCalOnlinePulser = new TObjArray(fNChannels);
  fTOFCalOnlineNoise = new TObjArray(fNChannels);
  fTOFCalOnlineHW = new TObjArray(fNChannels);
  fTOFCalOffline = new TObjArray(fNChannels);
  fTOFCalOnline->SetOwner();
  fTOFCalOnlinePulser->SetOwner();
  fTOFCalOnlineNoise->SetOwner();
  fTOFCalOnlineHW->SetOwner();
  fTOFCalOffline->SetOwner();

  //TOF Calibration Class copy ctor
  for (Int_t iarray = 0; iarray<fNChannels; iarray++){
    AliTOFChannelOnline * calChOnline = (AliTOFChannelOnline*)calib.fTOFCalOnline->At(iarray);
    AliTOFChannelOnlineStatus * calChOnlineStPulser = (AliTOFChannelOnlineStatus*)calib.fTOFCalOnlinePulser->At(iarray);
    AliTOFChannelOnlineStatus * calChOnlineStNoise = (AliTOFChannelOnlineStatus*)calib.fTOFCalOnlineNoise->At(iarray);
    AliTOFChannelOnlineStatus * calChOnlineStHW = (AliTOFChannelOnlineStatus*)calib.fTOFCalOnlineHW->At(iarray);
    AliTOFChannelOffline * calChOffline = (AliTOFChannelOffline*)calib.fTOFCalOffline->At(iarray);
    fTOFCalOnline->AddAt(calChOnline,iarray);
    fTOFCalOnlinePulser->AddAt(calChOnlineStPulser,iarray);
    fTOFCalOnlineNoise->AddAt(calChOnlineStNoise,iarray);
    fTOFCalOnlineHW->AddAt(calChOnlineStHW,iarray);
    fTOFCalOffline->AddAt(calChOffline,iarray);
  }

  if (calib.fDeltaBCOffset) fDeltaBCOffset = new AliTOFDeltaBCOffset(*calib.fDeltaBCOffset);
  if (calib.fCTPLatency) fCTPLatency = new AliTOFCTPLatency(*calib.fCTPLatency);
  if (calib.fT0Fill) fT0Fill = new AliTOFT0Fill(*calib.fT0Fill);
  if (calib.fRunParams) fRunParams = new AliTOFRunParams(*calib.fRunParams);
  if (calib.fResponseParams) fResponseParams = new AliTOFResponseParams(*calib.fResponseParams);
  if (calib.fReadoutEfficiency) fReadoutEfficiency = new TH1F(*calib.fReadoutEfficiency);
  if (calib.fProblematic) fProblematic = new TH1C(*calib.fProblematic);

  gRandom->SetSeed(123456789);
}

//____________________________________________________________________________ 

AliTOFcalib& AliTOFcalib::operator=(const AliTOFcalib &calib)
{
  //TOF Calibration Class assignment operator

  if (this == &calib)
    return *this;
  
  TTask::operator=(calib);
  fNChannels = calib.fNChannels;
  fCal = calib.fCal;
  fStatus = calib.fStatus;
  fTOFSimToT = calib.fTOFSimToT;
  fkValidity = calib.fkValidity;
  fTree = calib.fTree;
  fChain = calib.fChain;
  fNruns = calib.fNruns;
  fFirstRun = calib.fFirstRun;
  fLastRun = calib.fLastRun;
  for (Int_t iarray = 0; iarray<fNChannels; iarray++){
    AliTOFChannelOnline * calChOnline = (AliTOFChannelOnline*)calib.fTOFCalOnline->At(iarray);
    AliTOFChannelOffline * calChOffline = (AliTOFChannelOffline*)calib.fTOFCalOffline->At(iarray);
    AliTOFChannelOnlineStatus * calChOnlineStPulser = (AliTOFChannelOnlineStatus*)calib.fTOFCalOnlinePulser->At(iarray);
    AliTOFChannelOnlineStatus * calChOnlineStNoise = (AliTOFChannelOnlineStatus*)calib.fTOFCalOnlineNoise->At(iarray);
    AliTOFChannelOnlineStatus * calChOnlineStHW = (AliTOFChannelOnlineStatus*)calib.fTOFCalOnlineHW->At(iarray);
    fTOFCalOnline->AddAt(calChOnline,iarray);
    fTOFCalOnlinePulser->AddAt(calChOnlineStPulser,iarray);
    fTOFCalOnlineNoise->AddAt(calChOnlineStNoise,iarray);
    fTOFCalOnlineHW->AddAt(calChOnlineStHW,iarray);
    fTOFCalOffline->AddAt(calChOffline,iarray);
  }

  if (calib.fDeltaBCOffset) {
    if (fDeltaBCOffset) *fDeltaBCOffset = *calib.fDeltaBCOffset;
    else fDeltaBCOffset = new AliTOFDeltaBCOffset(*calib.fDeltaBCOffset);
  }

  if (calib.fCTPLatency) {
    if (fCTPLatency) *fCTPLatency = *calib.fCTPLatency;
    else fCTPLatency = new AliTOFCTPLatency(*calib.fCTPLatency);
  }

  if (calib.fT0Fill) {
    if (fT0Fill) *fT0Fill = *calib.fT0Fill;
    else fT0Fill = new AliTOFT0Fill(*calib.fT0Fill);
  }
  if (calib.fRunParams) {
    if (fRunParams) *fRunParams = *calib.fRunParams;
    else fRunParams = new AliTOFRunParams(*calib.fRunParams);
  }
  if (calib.fResponseParams) {
    if (fResponseParams) *fResponseParams = *calib.fResponseParams;
    else fResponseParams = new AliTOFResponseParams(*calib.fResponseParams);
  }
  if (calib.fReadoutEfficiency) {
    if (fReadoutEfficiency) *fReadoutEfficiency = *calib.fReadoutEfficiency;
    else fReadoutEfficiency = new TH1F(*calib.fReadoutEfficiency);
  }
  if (calib.fProblematic) {
    if (fProblematic) *fProblematic = *calib.fProblematic;
    else fProblematic = new TH1C(*calib.fProblematic);
  }
  fInitFlag = calib.fInitFlag;
  fRemoveMeanT0 = calib.fRemoveMeanT0;
  fUseLHCClockPhase = calib.fUseLHCClockPhase;
  fCalibrateTOFsignal = calib.fCalibrateTOFsignal;
  fCorrectTExp = calib.fCorrectTExp;

  return *this;
}

//____________________________________________________________________________ 

AliTOFcalib::~AliTOFcalib()
{
  //TOF Calibration Class dtor
  if(!(AliCDBManager::Instance()->GetCacheFlag())){ // CDB objects must NOT be deleted if cache is active!
    if (fTOFCalOnline){
      delete fTOFCalOnline;
    }
    if (fTOFCalOnlinePulser){
      delete fTOFCalOnlinePulser;
    }
    if (fTOFCalOnlineNoise){
      delete fTOFCalOnlineNoise;
    }
    if (fTOFCalOnlineHW){
      delete fTOFCalOnlineHW;
    }
    if (fTOFCalOffline){
      delete fTOFCalOffline;
    }
    if (fCal){
      delete fCal;
    }
    if (fStatus){
      delete fStatus;
    }
    if (fConfigMap){
      delete fConfigMap;
    }
    if (fDeltaBCOffset) delete fDeltaBCOffset;
    if (fCTPLatency) delete fCTPLatency;
    if (fT0Fill) delete fT0Fill;
    if (fRunParams) delete fRunParams;
    if (fResponseParams) delete fResponseParams;
    if (fReadoutEfficiency) delete fReadoutEfficiency;
    if (fProblematic) delete fProblematic;
  }
  if (fTree!=0x0) delete fTree;
  if (fChain!=0x0) delete fChain;

}
//_____________________________________________________________________________
void AliTOFcalib::CreateCalArrays(){

  // creating arrays for online/offline calibration objs

  fTOFCalOnline = new TObjArray(fNChannels);
  fTOFCalOnlinePulser = new TObjArray(fNChannels);
  fTOFCalOnlineNoise = new TObjArray(fNChannels);
  fTOFCalOnlineHW = new TObjArray(fNChannels);
  fTOFCalOffline = new TObjArray(fNChannels);
  fTOFCalOnline->SetOwner();
  fTOFCalOnlinePulser->SetOwner();
  fTOFCalOnlineNoise->SetOwner();
  fTOFCalOnlineHW->SetOwner();
  fTOFCalOffline->SetOwner();
  for (Int_t iarray = 0; iarray<fNChannels; iarray++){
    AliTOFChannelOnline * calChOnline = new AliTOFChannelOnline();
    AliTOFChannelOnlineStatus * calChOnlineStPulser = new AliTOFChannelOnlineStatus();
    AliTOFChannelOnlineStatus * calChOnlineStNoise = new AliTOFChannelOnlineStatus();
    AliTOFChannelOnlineStatus * calChOnlineStHW = new AliTOFChannelOnlineStatus();
    AliTOFChannelOffline * calChOffline = new AliTOFChannelOffline();
    fTOFCalOnline->AddAt(calChOnline,iarray);
    fTOFCalOnlinePulser->AddAt(calChOnlineStPulser,iarray);
    fTOFCalOnlineNoise->AddAt(calChOnlineStNoise,iarray);
    fTOFCalOnlineHW->AddAt(calChOnlineStHW,iarray);
    fTOFCalOffline->AddAt(calChOffline,iarray);
  }
  fCal = new AliTOFChannelOnlineArray(fNChannels);
  fStatus = new AliTOFChannelOnlineStatusArray(fNChannels);
}
//_____________________________________________________________________________
void AliTOFcalib::CreateCalObjects(){

  // creating arrays for online/offline calibration objs

  fTOFCalOffline = new TObjArray(fNChannels);
  fTOFCalOffline->SetOwner();
  for (Int_t iarray = 0; iarray<fNChannels; iarray++){
    AliTOFChannelOffline * calChOffline = new AliTOFChannelOffline();
    fTOFCalOffline->AddAt(calChOffline,iarray);
  }
  fCal = new AliTOFChannelOnlineArray(fNChannels);
  fStatus = new AliTOFChannelOnlineStatusArray(fNChannels);
}
//_____________________________________________________________________________
void AliTOFcalib::WriteConfigMapOnCDB(const Char_t *sel, Int_t minrun, Int_t maxrun)
{
  //Write calibration parameters to the CDB
  SetFirstRun(minrun);
  SetLastRun(maxrun);
  AliCDBManager *man = AliCDBManager::Instance();
  const Char_t *sel1 = "Config" ;  // to be consistent with TOFPreprocessor
  TString out(Form("%s/%s",sel,sel1));
  AliDebug(2,Form("Writing TOF configuration map for online calib on CDB with run range [%i, %i] ",fFirstRun,fLastRun));
  AliCDBId id(out,fFirstRun,fLastRun);
  AliCDBMetaData *md = new AliCDBMetaData();
  md->SetResponsible("Chiara Zampolli");
  if (!fConfigMap) {
    // deve uscire!!
  }
  man->Put(fConfigMap,id,md);
  delete md;
}
//_____________________________________________________________________________

void AliTOFcalib::WriteConfigMapOnCDB(const Char_t *sel)
{
  //Write calibration parameters to the CDB with infinite validity
  AliCDBManager *man = AliCDBManager::Instance();
  const Char_t *sel1 = "Config" ;  // to be consistent with TOFPreprocessor
  TString out(Form("%s/%s",sel,sel1));
  AliCDBRunRange runrange(fFirstRun,fLastRun);
  AliDebug(2,Form("Writing TOF config map for online calib on CDB with run range [%i, %i] ",runrange.GetFirstRun(),runrange.GetLastRun()));
  AliCDBId id(out,runrange);
  AliCDBMetaData *md = new AliCDBMetaData();
  md->SetResponsible("Chiara Zampolli");
  if (!fConfigMap) {
    // deve uscire!!
  }
  man->Put(fConfigMap,id,md);
  delete md;
}
//_____________________________________________________________________________
void AliTOFcalib::WriteParOnlineDelayOnCDB(const Char_t *sel, Int_t minrun, Int_t maxrun)
{
  //Write calibration parameters to the CDB -------> new calib objs!!!!!
  SetFirstRun(minrun);
  SetLastRun(maxrun);
  AliCDBManager *man = AliCDBManager::Instance();
  const Char_t *sel1 = "ParOnlineDelay" ;  // to be consistent with TOFPreprocessor
  TString out(Form("%s/%s",sel,sel1));
  AliDebug(2,Form("Writing TOF online calib obj on CDB with run range [%i, %i] ",fFirstRun,fLastRun));
  AliCDBId id(out,fFirstRun,fLastRun);
  AliCDBMetaData *md = new AliCDBMetaData();
  md->SetResponsible("Chiara Zampolli");
  if (!fCal) {
    // deve uscire!!
  }
  man->Put(fCal,id,md);
  delete md;
}
//_____________________________________________________________________________
void AliTOFcalib::WriteParOnlineStatusOnCDB(const Char_t *sel, Int_t minrun, Int_t maxrun)
{
  //Write calibration parameters to the CDB -------> new calib objs!!!!!
  SetFirstRun(minrun);
  SetLastRun(maxrun);
  AliCDBManager *man = AliCDBManager::Instance();
  const Char_t *sel1 = "Status" ;  // to be consistent with TOFPreprocessor
  TString out(Form("%s/%s",sel,sel1));
  AliDebug(2,Form("Writing TOF online status calib obj on CDB with run range [%i, %i] ",fFirstRun,fLastRun));
  AliCDBId id(out,fFirstRun,fLastRun);
  AliCDBMetaData *md = new AliCDBMetaData();
  md->SetResponsible("Chiara Zampolli");
  if (!fStatus) {
    // deve uscire!!
  }
  man->Put(fStatus,id,md);
  delete md;
}
//_____________________________________________________________________________

void AliTOFcalib::WriteParOnlineDelayOnCDB(const Char_t *sel)
{
  //Write calibration parameters to the CDB with infinite validity -------> new calib objs!!!!!
  AliCDBManager *man = AliCDBManager::Instance();
  const Char_t *sel1 = "ParOnlineDelay" ;  // to be consistent with TOFPreprocessor
  TString out(Form("%s/%s",sel,sel1));
  AliCDBRunRange runrange(fFirstRun,fLastRun);
  AliDebug(2,Form("Writing TOF online calib obj on CDB with run range [%i, %i] ",runrange.GetFirstRun(),runrange.GetLastRun()));
  AliCDBId id(out,runrange);
  AliCDBMetaData *md = new AliCDBMetaData();
  md->SetResponsible("Chiara Zampolli");
  if (!fCal) {
    // deve uscire!!
  }
  man->Put(fCal,id,md);
  delete md;
}
//_____________________________________________________________________________

void AliTOFcalib::WriteParOnlineStatusOnCDB(const Char_t *sel)
{
  //Write calibration parameters to the CDB with infinite validity -------> new calib objs!!!!!
  AliCDBManager *man = AliCDBManager::Instance();
  const Char_t *sel1 = "Status" ;  // to be consistent with TOFPreprocessor
  TString out(Form("%s/%s",sel,sel1));
  AliCDBRunRange runrange(fFirstRun,fLastRun);
  AliDebug(2,Form("Writing TOF online status calib obj on CDB with run range [%i, %i] ",runrange.GetFirstRun(),runrange.GetLastRun()));
  AliCDBId id(out,runrange);
  AliCDBMetaData *md = new AliCDBMetaData();
  md->SetResponsible("Chiara Zampolli");
  if (!fStatus) {
    // deve uscire!!
  }
  man->Put(fStatus,id,md);
  delete md;
}
//_____________________________________________________________________________
void AliTOFcalib::WriteParOnlineOnCDB(const Char_t *sel, Int_t minrun, Int_t maxrun)
{
  //Write calibration parameters to the CDB
  SetFirstRun(minrun);
  SetLastRun(maxrun);
  AliCDBManager *man = AliCDBManager::Instance();
  const Char_t *sel1 = "ParOnline" ;  // to be consistent with TOFPreprocessor
  TString out(Form("%s/%s",sel,sel1));
  AliDebug(2,Form("Writing TOF online calib obj on CDB with run range [%i, %i] ",fFirstRun,fLastRun));
  AliCDBId id(out,fFirstRun,fLastRun);
  AliCDBMetaData *md = new AliCDBMetaData();
  md->SetResponsible("Chiara Zampolli");
  if (!fTOFCalOnline) {
    // deve uscire!!
  }
  man->Put(fTOFCalOnline,id,md);
  delete md;
}
//_____________________________________________________________________________
void AliTOFcalib::WriteParOnlinePulserOnCDB(const Char_t *sel, Int_t minrun, Int_t maxrun)
{
  //Write calibration parameters from pulser to the CDB
  SetFirstRun(minrun);
  SetLastRun(maxrun);
  AliCDBManager *man = AliCDBManager::Instance();
  const Char_t *sel1 = "Pulser" ;  // to be consistent with TOFPreprocessor
  TString out(Form("%s/%s",sel,sel1));
  AliDebug(2,Form("Writing TOF online calib obj from pulser on CDB with run range [%i, %i] ",fFirstRun,fLastRun));
  AliCDBId id(out,fFirstRun,fLastRun);
  AliCDBMetaData *md = new AliCDBMetaData();
  md->SetResponsible("Chiara Zampolli");
  if (!fTOFCalOnlinePulser) {
    // deve uscire!!
  }
  man->Put(fTOFCalOnlinePulser,id,md);
  delete md;
}
//_____________________________________________________________________________
void AliTOFcalib::WriteParOnlineNoiseOnCDB(const Char_t *sel, Int_t minrun, Int_t maxrun)
{
  //Write calibration parameters from noise to the CDB
  SetFirstRun(minrun);
  SetLastRun(maxrun);
  AliCDBManager *man = AliCDBManager::Instance();
  const Char_t *sel1 = "Noise" ;  // to be consistent with TOFPreprocessor
  TString out(Form("%s/%s",sel,sel1));
  AliDebug(2,Form("Writing TOF online calib obj from noise on CDB with run range [%i, %i] ",fFirstRun,fLastRun));
  AliCDBId id(out,fFirstRun,fLastRun);
  AliCDBMetaData *md = new AliCDBMetaData();
  md->SetResponsible("Chiara Zampolli");
  if (!fTOFCalOnlineNoise) {
    // deve uscire!!
  }
  man->Put(fTOFCalOnlineNoise,id,md);
  delete md;
}
//_____________________________________________________________________________
void AliTOFcalib::WriteParOnlineHWOnCDB(const Char_t *sel, Int_t minrun, Int_t maxrun)
{
  //Write calibration parameters from hardware to the CDB
  SetFirstRun(minrun);
  SetLastRun(maxrun);
  AliCDBManager *man = AliCDBManager::Instance();
  const Char_t *sel1 = "HW" ;  // to be consistent with TOFPreprocessor
  TString out(Form("%s/%s",sel,sel1));
  AliDebug(2,Form("Writing TOF online calib obj from hardware on CDB with run range [%i, %i] ",fFirstRun,fLastRun));
  AliCDBId id(out,fFirstRun,fLastRun);
  AliCDBMetaData *md = new AliCDBMetaData();
  md->SetResponsible("Chiara Zampolli");
  if (!fTOFCalOnlineHW) {
    // deve uscire!!
  }
  man->Put(fTOFCalOnlineHW,id,md);
  delete md;
}
//_____________________________________________________________________________

void AliTOFcalib::WriteParOnlineOnCDB(const Char_t *sel)
{
  //Write calibration parameters to the CDB with infinite validity
  AliCDBManager *man = AliCDBManager::Instance();
  const Char_t *sel1 = "ParOnline" ;  // to be consistent with TOFPreprocessor
  TString out(Form("%s/%s",sel,sel1));
  AliCDBRunRange runrange(fFirstRun,fLastRun);
  AliDebug(2,Form("Writing TOF online calib obj on CDB with run range [%i, %i] ",runrange.GetFirstRun(),runrange.GetLastRun()));
  AliCDBId id(out,runrange);
  AliCDBMetaData *md = new AliCDBMetaData();
  md->SetResponsible("Chiara Zampolli");
  if (!fTOFCalOnline) {
    // deve uscire!!
  }
  man->Put(fTOFCalOnline,id,md);
  delete md;
}
//_____________________________________________________________________________

void AliTOFcalib::WriteParOnlinePulserOnCDB(const Char_t *sel)
{
  //Write calibration parameters from pulser to the CDB with infinite validity
  AliCDBManager *man = AliCDBManager::Instance();
  const Char_t *sel1 = "Pulser" ;  // to be consistent with TOFPreprocessor
  TString out(Form("%s/%s",sel,sel1));
  AliCDBRunRange runrange(fFirstRun,fLastRun);
  AliDebug(2,Form("Writing TOF online calib obj from pulser on CDB with run range [%i, %i] ",runrange.GetFirstRun(),runrange.GetLastRun()));
  AliCDBId id(out,runrange);
  AliCDBMetaData *md = new AliCDBMetaData();
  md->SetResponsible("Chiara Zampolli");
  if (!fTOFCalOnlinePulser) {
    // deve uscire!!
  }
  man->Put(fTOFCalOnlinePulser,id,md);
  delete md;
}
//_____________________________________________________________________________

void AliTOFcalib::WriteParOnlineNoiseOnCDB(const Char_t *sel)
{
  //Write calibration parameters from noise to the CDB with infinite validity
  AliCDBManager *man = AliCDBManager::Instance();
  const Char_t *sel1 = "Noise" ;  // to be consistent with TOFPreprocessor
  TString out(Form("%s/%s",sel,sel1));
  AliCDBRunRange runrange(fFirstRun,fLastRun);
  AliDebug(2,Form("Writing TOF online calib obj from noise on CDB with run range [%i, %i] ",runrange.GetFirstRun(),runrange.GetLastRun()));
  AliCDBId id(out,runrange);
  AliCDBMetaData *md = new AliCDBMetaData();
  md->SetResponsible("Chiara Zampolli");
  if (!fTOFCalOnlineNoise) {
    // deve uscire!!
  }
  man->Put(fTOFCalOnlineNoise,id,md);
  delete md;
}
//_____________________________________________________________________________

void AliTOFcalib::WriteParOnlineHWOnCDB(const Char_t *sel)
{
  //Write calibration parameters from hardware to the CDB with infinite validity
  AliCDBManager *man = AliCDBManager::Instance();
  const Char_t *sel1 = "HW" ;  // to be consistent with TOFPreprocessor
  TString out(Form("%s/%s",sel,sel1));
  AliCDBRunRange runrange(fFirstRun,fLastRun);
  AliDebug(2,Form("Writing TOF online calib obj from harware on CDB with run range [%i, %i] ",runrange.GetFirstRun(),runrange.GetLastRun()));
  AliCDBId id(out,runrange);
  AliCDBMetaData *md = new AliCDBMetaData();
  md->SetResponsible("Chiara Zampolli");
  if (!fTOFCalOnlineHW) {
    // deve uscire!!
  }
  man->Put(fTOFCalOnlineHW,id,md);
  delete md;
}
//_____________________________________________________________________________

void AliTOFcalib::WriteParOfflineOnCDB(const Char_t *sel, const Char_t *validity, Int_t minrun, Int_t maxrun)
{
  //Write calibration parameters to the CDB
  SetFirstRun(minrun);
  SetLastRun(maxrun);
  AliCDBManager *man = AliCDBManager::Instance();
  const Char_t *sel1 = "ParOffline" ;
  TString out(Form("%s/%s",sel,sel1));
  AliDebug(2,Form("Writing TOF offline calib obj on CDB with run range [%i, %i] ",fFirstRun,fLastRun));
  AliCDBId id(out,fFirstRun,fLastRun);
  AliCDBMetaData *md = new AliCDBMetaData();
  md->SetResponsible("Chiara Zampolli");
  md->SetComment(validity);
  man->Put(fTOFCalOffline,id,md);
  delete md;
}
//_____________________________________________________________________________

void AliTOFcalib::WriteParOfflineOnCDB(const Char_t *sel, const Char_t *validity)
{
  //Write calibration parameters to the CDB with infinite validity
  AliCDBManager *man = AliCDBManager::Instance();
  const Char_t *sel1 = "ParOffline" ;
  TString out(Form("%s/%s",sel,sel1));
  AliCDBRunRange runrange(fFirstRun,fLastRun);
  AliDebug(2,Form("Writing TOF offline calib obj on CDB with run range [%i, %i] ",runrange.GetFirstRun(),runrange.GetLastRun()));
  AliCDBId id(out,runrange);
  AliCDBMetaData *md = new AliCDBMetaData();
  md->SetResponsible("Chiara Zampolli");
  md->SetComment(validity);
  man->Put(fTOFCalOffline,id,md);
  delete md;
}
//_____________________________________________________________________________

Bool_t AliTOFcalib::ReadConfigMapFromCDB(const Char_t *sel, Int_t nrun)
{
  //Read calibration parameters from the CDB
  AliCDBManager *man = AliCDBManager::Instance();
  const Char_t *sel1 = "Config" ;
  TString out(Form("%s/%s",sel,sel1));
  AliCDBEntry *entry = man->Get(out,nrun);
  if (!entry) { 
    AliFatal("Exiting, no CDB object (ConfigMap) found!!!");
    exit(0);  
  }
  if(!entry->GetObject()){
    AliFatal("Exiting, no CDB object (ConfigMap) found!!!");
    exit(0);  
  }  
  
  fConfigMap =(TMap*)entry->GetObject();

  return kTRUE; 
   
}
//_____________________________________________________________________________

Bool_t AliTOFcalib::ReadParOnlineDelayFromCDB(const Char_t *sel, Int_t nrun)
{
  //Read calibration parameters from the CDB -------> new calib objs!!!!!
  AliCDBManager *man = AliCDBManager::Instance();
  const Char_t *sel1 = "ParOnlineDelay" ;
  TString out(Form("%s/%s",sel,sel1));
  AliCDBEntry *entry = man->Get(out,nrun);
  if (!entry) { 
    AliFatal("Exiting, no CDB object (ParOnlineDelay) found!!!");
    exit(0);  
  }
  if(!entry->GetObject()){
    AliFatal("Exiting, no CDB object (ParOnlineDelay) found!!!");
    exit(0);  
  }  
  
  fCal =(AliTOFChannelOnlineArray*)entry->GetObject();

  return kTRUE; 
   
}
//_____________________________________________________________________________

Bool_t AliTOFcalib::ReadParOnlineStatusFromCDB(const Char_t *sel, Int_t nrun)
{
  //Read calibration parameters from the CDB -------> new calib objs!!!!!
  AliCDBManager *man = AliCDBManager::Instance();
  const Char_t *sel1 = "Status" ;
  TString out(Form("%s/%s",sel,sel1));
  AliCDBEntry *entry = man->Get(out,nrun);
  if (!entry) { 
    AliFatal("Exiting, no CDB object (Status) found!!!");
    exit(0);  
  }
  if(!entry->GetObject()){
    AliFatal("Exiting, no CDB object (Status) found!!!");
    exit(0);  
  }  
  
  fStatus =(AliTOFChannelOnlineStatusArray*)entry->GetObject();

  return kTRUE; 
   
}
//_____________________________________________________________________________

Bool_t AliTOFcalib::ReadParOnlineFromCDB(const Char_t *sel, Int_t nrun)
{
  //Read calibration parameters from the CDB
  AliCDBManager *man = AliCDBManager::Instance();
  const Char_t *sel1 = "ParOnline" ;
  TString out(Form("%s/%s",sel,sel1));
  AliCDBEntry *entry = man->Get(out,nrun);
  if (!entry) { 
    AliFatal("Exiting, no CDB object (ParOnline) found!!!");
    exit(0);  
  }
  if(!entry->GetObject()){
    AliFatal("Exiting, no CDB object (ParOnline) found!!!");
    exit(0);  
  }  
  
  fTOFCalOnline =(TObjArray*)entry->GetObject();

  return kTRUE; 
   
}
//_____________________________________________________________________________

Bool_t AliTOFcalib::ReadParOnlinePulserFromCDB(const Char_t *sel, Int_t nrun)
{
  //Read calibration parameters from pulser from the CDB
  AliCDBManager *man = AliCDBManager::Instance();
  const Char_t *sel1 = "Pulser" ;
  TString out(Form("%s/%s",sel,sel1));
  AliCDBEntry *entry = man->Get(out,nrun);
  if (!entry) { 
    AliFatal("Exiting, no CDB object (Pulser) found!!!");
    exit(0);  
  }
  if(!entry->GetObject()){
    AliFatal("Exiting, no CDB object (Pulser) found!!!");
    exit(0);  
  }  
  
  fTOFCalOnlinePulser =(TObjArray*)entry->GetObject();

  return kTRUE; 
   
}
//_____________________________________________________________________________

Bool_t AliTOFcalib::ReadParOnlineNoiseFromCDB(const Char_t *sel, Int_t nrun)
{
  //Read calibration parameters from noise from the CDB
  AliCDBManager *man = AliCDBManager::Instance();
  const Char_t *sel1 = "Noise" ;
  TString out(Form("%s/%s",sel,sel1));
  AliCDBEntry *entry = man->Get(out,nrun);
  if (!entry) { 
    AliFatal("Exiting, no CDB object (Noise) found!!!");
    exit(0);  
  }
  if(!entry->GetObject()){
    AliFatal("Exiting, no CDB object (Noise) found!!!");
    exit(0);  
  }  
  
  fTOFCalOnlineNoise =(TObjArray*)entry->GetObject();

  return kTRUE; 
   
}
//_____________________________________________________________________________

Bool_t AliTOFcalib::ReadParOnlineHWFromCDB(const Char_t *sel, Int_t nrun)
{
  //Read calibration parameters from hardware from the CDB
  AliCDBManager *man = AliCDBManager::Instance();
  const Char_t *sel1 = "HW" ;
  TString out(Form("%s/%s",sel,sel1));
  AliCDBEntry *entry = man->Get(out,nrun);
  if (!entry) { 
    AliFatal("Exiting, no CDB object (HW map) found!!!");
    exit(0);  
  }
  if(!entry->GetObject()){
    AliFatal("Exiting, no CDB object (HW map) found!!!");
    exit(0);  
  }  
  
  fTOFCalOnlineHW =(TObjArray*)entry->GetObject();

  return kTRUE; 
   
}
//_____________________________________________________________________________

Bool_t AliTOFcalib::ReadParOfflineFromCDB(const Char_t *sel, Int_t nrun)
{
  //Read calibration parameters from the CDB
  AliCDBManager *man = AliCDBManager::Instance();
  const Char_t *sel1 = "ParOffline" ;
  TString out(Form("%s/%s",sel,sel1));
  AliCDBEntry *entry = man->Get(out,nrun);
  if (!entry) { 
    AliFatal("Exiting, no CDB object (ParOffline) found!!!");
    exit(0);  
  }
  if(!entry->GetObject()){
    AliFatal("Exiting, no CDB object (ParOffline) found!!!");
    exit(0);  
  }  
  AliCDBMetaData * md = entry->GetMetaData();
  fkValidity = md->GetComment();  
  fTOFCalOffline =(TObjArray*)entry->GetObject();

  return kTRUE; 
   
}
//_____________________________________________________________________________
void AliTOFcalib::WriteSimHistoOnCDB(const Char_t *sel, Int_t minrun, Int_t maxrun, TH1F *histo){
  //Write Sim miscalibration parameters to the CDB

  fTOFSimToT=histo;
  AliCDBManager *man = AliCDBManager::Instance();
  const Char_t *sel1 = "SimHisto" ;
  TString out(Form("%s/%s",sel,sel1));
  AliCDBMetaData *mdhisto = new AliCDBMetaData();
  mdhisto->SetResponsible("Chiara Zampolli");
  AliCDBId id(out,minrun,maxrun);
  man->Put(fTOFSimToT,id,mdhisto);
  delete mdhisto;
}
//_____________________________________________________________________________
Bool_t AliTOFcalib::ReadSimHistoFromCDB(const Char_t *sel, Int_t nrun)
{
  //Read miscalibration parameters from the CDB
  AliCDBManager *man = AliCDBManager::Instance();

  // The Tot Histo

  const Char_t *sel1 = "SimHisto" ;
  TString out(Form("%s/%s",sel,sel1));
  AliCDBEntry *entry = man->Get(out,nrun);
  if (!entry) { 
    AliFatal("Exiting, no CDB object (SimHisto) found!!!");
    exit(0);  
  }
  if(!entry->GetObject()){
    AliFatal("Exiting, no CDB object (SimHisto) found!!!");
    exit(0);  
  }  
  TH1F *histo =(TH1F*)entry->GetObject();
  fTOFSimToT=histo;
  return kTRUE;
}
//_____________________________________________________________________________
void AliTOFcalib::WriteRecParOnCDB(const Char_t *sel, Int_t minrun, Int_t maxrun, AliTOFRecoParam *param){
  //Write reconstruction parameters to the CDB

  AliCDBManager *man = AliCDBManager::Instance();
  AliCDBMetaData *md = new AliCDBMetaData();
  md->SetResponsible("Silvia Arcelli");
  const Char_t *sel1 = "RecoParam" ;
  TString out(Form("%s/%s",sel,sel1));
  AliCDBId id(out,minrun,maxrun);

  TObjArray *arr=new TObjArray(1);
  arr->AddLast(param);
  man->Put(arr,id,md);
  //man->Put(param,id,md);
  delete md;
}
//_____________________________________________________________________________
void AliTOFcalib::WriteRecParOnCDB(const Char_t *sel, Int_t minrun, Int_t maxrun, TObjArray *arr){
  //Write reconstruction parameters to the CDB

  AliCDBManager *man = AliCDBManager::Instance();
  AliCDBMetaData *md = new AliCDBMetaData();
  md->SetResponsible("Silvia Arcelli");
  const Char_t *sel1 = "RecoParam" ;
  TString out(Form("%s/%s",sel,sel1));
  AliCDBId id(out,minrun,maxrun);
  man->Put(arr,id,md);
  delete md;
}
//_____________________________________________________________________________
AliTOFRecoParam * AliTOFcalib::ReadRecParFromCDB(const Char_t *sel, Int_t nrun, Int_t eventType)
{
  //Read reconstruction parameters from the CDB
  AliCDBManager *man = AliCDBManager::Instance();
  const Char_t *sel1 = "RecoParam" ;
  TString out(Form("%s/%s",sel,sel1));
  AliCDBEntry *entry = man->Get(out,nrun);
  if (!entry) { 
    AliFatal("Exiting, no CDB object (RecoParam) found!!!");
    exit(0);  
  }  
  if(!entry->GetObject()){
    AliFatal("Exiting, no CDB object (RecoParam) found!!!");
    exit(0);  
  }  

  TObjArray *array = (TObjArray*)entry->GetObject();
  AliTOFRecoParam *param=0x0;
  if (eventType>=0 || eventType<array->GetEntries())
    param=(AliTOFRecoParam*)array->At(eventType);
  return param;

}
//-----------------------------------------------------------------------------
// Calibration methods
//-----------------------------------------------------------------------------
void AliTOFcalib::CreateTreeFromCDB(Int_t minrun, Int_t maxrun){

  // creating the chain with the trees for calibration
  // collecting them from reference data 
  // from minrun to maxrun

  Float_t p[CHENTRIESSMALL];
  Int_t nentries;
  fTree = new TTree("TOFCalib","Tree for TOF Calibration");
  fTree->Branch("nentries",&nentries,"nentries/I");
  fTree->Branch("TOFentries",p,"TOFentries[nentries]/F");
  AliCDBManager *man = AliCDBManager::Instance();
  AliCDBStorage *aStorage = man->GetStorage("local://$ALICE_ROOT/OCDB");
  for (Int_t irun = minrun;irun<=maxrun;irun++){
    AliCDBEntry *entry = aStorage->Get("TOF/RefData/TreeForCalib",irun);
    if (!entry){
      AliInfo(Form("No entry found for run %i",irun));
    }
    else{
      TTree *tree = (TTree*)entry->GetObject();
      tree->SetBranchAddress("nentries",&nentries);
      tree->SetBranchAddress("TOFentries",p);      
      fTree->CopyEntries(tree);
      fNruns++;
    }
  }
  AliInfo(Form("Number of runs being analyzed %i",fNruns));
}
//-----------------------------------------------------------------------------
void AliTOFcalib::CreateTreeFromGrid(Int_t minrun, Int_t maxrun){

  // creating the chain with the trees for calibration
  // collecting them from the Grid 
  // from minrun to maxrun

  Float_t p[CHENTRIESSMALL];
  Int_t nentries;
  fTree = new TTree("TOFCalib","Tree for TOF Calibration");
  fTree->SetDirectory(0);
  fTree->Branch("nentries",&nentries,"nentries/I");
  fTree->Branch("TOFentries",p,"TOFentries[nentries]/F");
  AliInfo("connected to alien");
  TGrid::Connect("alien://");
  
  TString filename;
  for (Int_t irun = minrun;irun<=maxrun;irun++){
    filename = Form("alien:///alice/cern.ch/user/c/czampolli/TOFCalibReference_%i.root",irun);
    TFile *filegrid = TFile::Open(filename.Data(),"READ");
    TTree *tree = (TTree*)filegrid->Get("T");
    tree->SetBranchAddress("nentries",&nentries);
    tree->SetBranchAddress("TOFentries",p);      
    fTree->CopyEntries(tree);
    delete tree;
    fNruns++;    
  }
  
  AliInfo(Form("Number of runs being analyzed %i",fNruns));
}
//-----------------------------------------------------------------------------
void AliTOFcalib::CreateTreeFromFile(Int_t minrun, Int_t maxrun){

  // creating the tree with the trees for calibration
  // collecting them from reference data (from file)
  // from minrun to maxrun

  Float_t p[CHENTRIESSMALL];
  Int_t nentries;
  fTree = new TTree("TOFCalib","Tree for TOF Calibration");
  fTree->SetDirectory(0);
  fTree->Branch("nentries",&nentries,"nentries/I");
  fTree->Branch("TOFentries",p,"TOFentries[nentries]/F");
  TString filename;
  for (Int_t irun = minrun;irun<=maxrun;irun++){
    filename = Form("$ALICE_ROOT/TOF/RefData/TreeForCalib/fileout_%i.root",irun);
    TFile *file = new TFile(filename.Data(),"READ");
    TTree *tree = (TTree*)file->Get("T");
    tree->SetBranchAddress("nentries",&nentries);
    tree->SetBranchAddress("TOFentries",p);      
    fTree->CopyEntries(tree);
    delete tree;
    delete file;
    file = 0x0;
    fNruns++;
  }

  AliInfo(Form("Number of runs being analyzed %i",fNruns));
}
//-----------------------------------------------------------------------------
void AliTOFcalib::CreateChainFromGrid(Int_t minrun, Int_t maxrun){

  // creating the chain with the trees for calibration
  // collecting them from the Grid 
  // from minrun to maxrun

  fChain = new TChain("T");
  AliInfo("connected to alien");
  TGrid::Connect("alien://");
  
  TString filename;
  for (Int_t irun = minrun;irun<=maxrun;irun++){
    filename = Form("alien:///alice/cern.ch/user/c/czampolli/TOFCalibReference_%i.root",irun);
    fChain->Add(filename.Data());
    fNruns++;    
  }
  
  AliInfo(Form("Number of runs being analyzed %i",fNruns));
}
//-----------------------------------------------------------------------------
Int_t AliTOFcalib::Calibrate(Int_t ichmin, Int_t ichmax, Option_t *optionSave, Option_t *optionFit){

  // calibrating summing more than one channels
  // computing calibration parameters
  // Returning codes:
  // 0 -> everything was ok
  // 1 -> no tree for calibration found
  // 2 -> not enough statistics to perform calibration
  // 3 -> problems with arrays
  
  TH1::AddDirectory(0);

  AliInfo(Form("*** Calibrating Histograms %s, summing more channels, from channel %i, to channel %i, storing Calib Pars in channel %i", GetName(),ichmin,ichmax,ichmin)) ; 
  AliInfo(Form("Option for Saving histos = %s",optionSave )) ; 
  AliInfo(Form("Option for Fitting Profile histos = %s",optionFit )) ; 

  Float_t p[CHENTRIESSMALL];
  Int_t nentries;
  //fTree->SetBranchAddress("nentries",&nentries);
  //fTree->SetBranchAddress("TOFentries",p);
  fChain->SetBranchAddress("nentries",&nentries);
  fChain->SetBranchAddress("TOFentries",p);

  Float_t ntracksTotalmean =0;
  for (Int_t i=ichmin; i<ichmax; i++){
    Int_t ientry = -1;
    for (Int_t irun=0;irun<fNruns;irun++){
      ientry = i+irun*fNChannels;
      //fTree->GetEntry(ientry);
      fChain->GetEntry(ientry);
      Int_t ntracksRun=nentries/3;
      ntracksTotalmean+=ntracksRun;
    }
  }
  
  if (ntracksTotalmean < MEANENTRIES) {
    AliInfo(Form(" Too small mean number of entires per channel (mean number = %f) not calibrating and exiting.....",ntracksTotalmean));
    return 2;
  }

  //filling ToT and Time arrays

  Int_t nbinToT = 100;  // ToT bin width in Profile = 48.8 ps 
  Float_t minToT = 0;   // ns
  Float_t maxToT = 4.88;  // ns

  TH1F *hToT = new TH1F("htot","htot",nbinToT, minToT, maxToT);
  TH1F *hdeltaTime = new TH1F("hdeltaTime","hdeltaTime",200,2,4);
  Int_t ntracksTotal = 0;
  Int_t ntracksRun = 0;
  Double_t binsProfile[101]; // sized larger than necessary, the correct 
                             // dim being set in the booking of the profile
  Int_t nusefulbins=0;
  Float_t meantime=0;
  for (Int_t i = ichmin;i<ichmax;i++){
    Int_t ientry = -1;
    for (Int_t irun=0;irun<fNruns;irun++){
      ientry = i+irun*fNChannels;
      //fTree->GetEntry(ientry);
      fChain->GetEntry(ientry);
      ntracksTotal+=nentries/3;
      ntracksRun=nentries/3;
      AliDebug(2,Form("run %i, channel %i, nentries = %i, ntracks = %i",irun,i,nentries, ntracksRun));
      for (Int_t j=0;j<ntracksRun;j++){
	Int_t idxexToT = (j* NIDXSMALL)+DELTAIDXTOT; 
	Int_t idxexTime = (j* NIDXSMALL)+DELTAIDXTIME; 
	Int_t idxexExTime = (j* NIDXSMALL)+DELTAIDXPID; 
	Float_t tot = p[idxexToT];
	hdeltaTime->Fill(p[idxexTime]-p[idxexExTime]);
	meantime+=p[idxexTime]-p[idxexExTime];
	hToT->Fill(tot);
      }
    }
  }
  nusefulbins = FindBins(hToT,&binsProfile[0]);
  meantime/=ntracksTotal;
  AliDebug(2, Form("meantime = %f",meantime));
  
  for (Int_t j=1;j<=nusefulbins;j++) {
    AliDebug(2,Form(" summing channels from %i to %i, nusefulbins = %i, bin %i = %f",ichmin,ichmax,nusefulbins,j,binsProfile[j])); 
  }

  TProfile* hSlewingProf = new TProfile("hSlewingProf", "hSlewingProf",nusefulbins, binsProfile, "G");  // CHECK THE BUILD OPTION, PLEASE!!!!!!
  TH2F * htimetot = new TH2F("htimetot","htimetot",nbinToT, minToT, maxToT,600,-5,10);

  for (Int_t irun=0;irun<fNruns;irun++){
    Int_t ientry = -1;
    for (Int_t i=ichmin; i<ichmax; i++){
      ientry = i+irun*fNChannels;
      //fTree->GetEntry(ientry);
      fChain->GetEntry(ientry);
      ntracksRun=nentries/3;
      for (Int_t j=0;j<ntracksRun;j++){
	Int_t idxexToT = (j* NIDXSMALL)+DELTAIDXTOT; 
	Int_t idxexTime = (j* NIDXSMALL)+DELTAIDXTIME; 
	Int_t idxexExTime = (j* NIDXSMALL)+DELTAIDXPID; 
	Float_t tot = p[idxexToT];
	Float_t time = p[idxexTime]-p[idxexExTime];
	AliDebug (2, Form("track = %i, time = %f, tot = %f, time-meantime = %f",j,time, tot, time-meantime));
	hSlewingProf->Fill(tot,time);  // if meantime is not used, the fill may be moved in the loop above
	htimetot->Fill(tot,time-meantime);  // if meantime is not used, the fill may be moved in the loop above
      }
    }
  }

  hSlewingProf->Fit("pol5",optionFit,"",0,4);
  TF1 * calibfunc = (TF1*)hSlewingProf->GetFunction("pol5");
  Float_t par[6];    
  for(Int_t kk=0;kk<6;kk++){
    par[kk]=calibfunc->GetParameter(kk);
    AliDebug(2,Form("parameter %i = %f",kk,par[kk]));
  }

  if(strstr(optionSave,"save")){
    TFile * fileProf = new TFile("TOFCalibSave.root","recreate");
    fileProf->cd(); 
    TString profName=Form("Profile%06i_%06i",ichmin,ichmax);
    TString timeTotName=Form("TimeTot%06i_%06i",ichmin,ichmax);
    TString totName=Form("Tot%06i_%06i",ichmin,ichmax);
    TString deltaName=Form("Delta%06i_%06i",ichmin,ichmax);
    hSlewingProf->Write(profName);
    htimetot->Write(timeTotName);
    hToT->Write(totName);
    hdeltaTime->Write(deltaName);
    fileProf->Close();
    delete fileProf;
    fileProf=0x0;
  }

  delete hToT;
  hToT=0x0;
  delete hSlewingProf;
  hSlewingProf=0x0;
  delete htimetot;
  htimetot=0x0;
  delete hdeltaTime;
  hdeltaTime=0x0;

  AliTOFChannelOffline * calChannel = (AliTOFChannelOffline*)fTOFCalOffline->At(ichmin);
  calChannel->SetSlewPar(par);
  WriteParOfflineOnCDB("TOF/Calib","valid");
  return 0;
}
//----------------------------------------------------------------------------
Int_t AliTOFcalib::Calibrate(Int_t i, Option_t *optionSave, Option_t *optionFit){

  // computing calibration parameters for channel i
  // Returning codes:
  // 0 -> everything was ok
  // 1 -> no tree for calibration found
  // 2 -> not enough statistics to perform calibration
  // 3 -> problems with arrays

  TH1::AddDirectory(0);
  
  AliInfo(Form("*** Calibrating Histograms (one channel) %s", GetName())) ; 
  AliInfo(Form("Option for Saving histos = %s",optionSave )) ; 
  AliInfo(Form("Option for Fitting Profile histos = %s",optionFit )) ; 

  Float_t p[MAXCHENTRIESSMALL];
  Int_t nentries;
  //fTree->SetBranchAddress("nentries",&nentries);
  //fTree->SetBranchAddress("TOFentries",p);
  fChain->SetBranchAddress("nentries",&nentries);
  fChain->SetBranchAddress("TOFentries",p);

  Float_t ntracksTotal =0;
  for (Int_t irun=0;irun<fNruns;irun++){
    Int_t ientry = -1;
    ientry = i+irun*fNChannels;
    //fTree->GetEntry(ientry);
    fChain->GetEntry(ientry);
    ntracksTotal+=nentries/3;    
  }
  
  if (ntracksTotal < MEANENTRIES) {  
    AliInfo(Form(" Too small mean number of entires per channel (mean number = %f) not calibrating and exiting.....",ntracksTotal));
    return 2;
  }

  //filling ToT and Time arrays

  Int_t nbinToT = 100;  // ToT bin width in Profile = 48.8 ps 
  Float_t minToT = 0;   // ns
  Float_t maxToT = 4.88;  // ns

  TH1F *hToT = new TH1F("htot","htot",nbinToT, minToT, maxToT);
  TH1F *hdeltaTime = new TH1F("hdeltaTime","hdeltaTime",200,2,4);
  Int_t ntracksRun = 0;
  Double_t binsProfile[101]; // sized larger than necessary, the correct 
                             // dim being set in the booking of the profile
  Int_t nusefulbins=0;
  Float_t meantime=0;
  for (Int_t irun=0;irun<fNruns;irun++){
    Int_t ientry = -1;
    ientry = i+irun*fNChannels;
    //fTree->GetEntry(ientry);
    fChain->GetEntry(ientry);
    ntracksRun=nentries/3;
    AliDebug(2,Form("run %i, channel %i, nentries = %i, ntracksRun = %i",irun, i ,nentries, ntracksRun));
    for (Int_t j=0;j<ntracksRun;j++){
      Int_t idxexToT = (j* NIDXSMALL)+DELTAIDXTOT; 
      Int_t idxexTime = (j* NIDXSMALL)+DELTAIDXTIME; 
      Int_t idxexExTime = (j* NIDXSMALL)+DELTAIDXPID; 
      Float_t tot = p[idxexToT];
      meantime+=p[idxexTime]-p[idxexExTime];
      hdeltaTime->Fill(p[idxexTime]-p[idxexExTime]);
      hToT->Fill(tot);
    }
  }

  nusefulbins = FindBins(hToT,&binsProfile[0]);
  meantime/=ntracksTotal;
  AliDebug(2,Form("meantime = %f",meantime));
  
  for (Int_t j=1;j<=nusefulbins;j++) {
    AliDebug(2,Form(" channel %i, nusefulbins = %i, bin %i = %f",i,nusefulbins,j,binsProfile[j])); 
  }

  TProfile* hSlewingProf = new TProfile("hSlewingProf", "hSlewingProf",nusefulbins, binsProfile, "G");  // CHECK THE BUILD OPTION, PLEASE!!!!!!
  TH2F * htimetot = new TH2F("htimetot","htimetot",nbinToT, minToT, maxToT,600,-5,10);
  for (Int_t irun=0;irun<fNruns;irun++){
    Int_t ientry = -1;
    ientry = i+irun*fNChannels;
    //fTree->GetEntry(ientry);
    fChain->GetEntry(ientry);
    ntracksRun=nentries/3;
    for (Int_t j=0;j<ntracksRun;j++){
      Int_t idxexToT = (j* NIDXSMALL)+DELTAIDXTOT; 
      Int_t idxexTime = (j* NIDXSMALL)+DELTAIDXTIME; 
      Int_t idxexExTime = (j* NIDXSMALL)+DELTAIDXPID; 
      Float_t tot = p[idxexToT];
      Float_t time = p[idxexTime]-p[idxexExTime];
      AliDebug (2,Form("track = %i, time = %f, tot = %f, time-meantime = %f",j,time, tot, time-meantime));
      hSlewingProf->Fill(tot,time);  // if meantime is not used, the fill may be moved in the loop above
      htimetot->Fill(tot,time-meantime);  // if meantime is not used, the fill may be moved in the loop above
    }
  }

  hSlewingProf->Fit("pol5",optionFit,"",0,4);
  TF1 * calibfunc = (TF1*)hSlewingProf->GetFunction("pol5");
  Float_t par[6];    
  for(Int_t kk=0;kk<6;kk++){
    par[kk]=calibfunc->GetParameter(kk);
    AliDebug(2,Form("parameter %i = %f",kk,par[kk]));
  }


  if(strstr(optionSave,"save")){
    TFile * fileProf = new TFile("TOFCalibSave.root","recreate");
    fileProf->cd();   
    TString profName=Form("Profile%06i",i);
    TString timeTotName=Form("TimeTot%06i",i);
    TString totName=Form("Tot%06i",i);
    TString deltaName=Form("Delta%06i",i);
    hSlewingProf->Write(profName);
    htimetot->Write(timeTotName);
    hToT->Write(totName);
    hdeltaTime->Write(deltaName);
    fileProf->Close();
    delete fileProf;
    fileProf=0x0;
  }

  delete hToT;
  hToT=0x0; 
  delete hSlewingProf;
  hSlewingProf=0x0;
  delete htimetot;
  htimetot=0x0;
  delete hdeltaTime;
  hdeltaTime=0x0;

  AliTOFChannelOffline * calChannel = (AliTOFChannelOffline*)fTOFCalOffline->At(i);
  calChannel->SetSlewPar(par);
  WriteParOfflineOnCDB("TOF/Calib","valid");
  return 0;
}
//----------------------------------------------------------------------------
Int_t AliTOFcalib::Calibrate(Int_t nch, Int_t *ch, Option_t *optionSave, Option_t *optionFit){

  // calibrating an array of channels
  // computing calibration parameters
  // Returning codes:
  // 0 -> everything was ok
  // 1 -> no tree for calibration found
  // 2 -> not enough statistics to perform calibration
  // 3 -> problems with arrays
  
  TH1::AddDirectory(0);

  AliInfo(Form("*** Calibrating Histograms %s, number of channels = %i", GetName(),nch)) ; 
  AliInfo(Form("Option for Saving histos = %s",optionSave )) ; 
  AliInfo(Form("Option for Fitting Profile histos = %s",optionFit )) ; 
  for (Int_t ich=0; ich<nch; ich++){
    Int_t i = ch[ich];
    AliInfo(Form("Calibrating channel = %i",i )) ; 
  }
  Float_t p[MAXCHENTRIESSMALL];
  Int_t nentries;
  //fTree->SetBranchAddress("nentries",&nentries);
  //fTree->SetBranchAddress("TOFentries",p);
  fChain->SetBranchAddress("nentries",&nentries);
  fChain->SetBranchAddress("TOFentries",p);

  Float_t ntracksTotalmean =0;
  for (Int_t ich=0; ich<nch; ich++){
    Int_t ientry = -1;
      Int_t i = ch[ich];
      for (Int_t irun=0;irun<fNruns;irun++){
      ientry = i+irun*fNChannels;
      //fTree->GetEntry(ientry);
      fChain->GetEntry(ientry);
      ntracksTotalmean+=nentries/3;
    }
  }

  ntracksTotalmean/=nch;
  if (ntracksTotalmean < MEANENTRIES) { 
    AliInfo(Form(" Too small mean number of entires per channel (mean number = %f) not calibrating and exiting.....",ntracksTotalmean));
    return 2;
  }

  //filling ToT and Time arrays

  Int_t nbinToT = 100;  // ToT bin width in Profile = 48.8 ps 
  Float_t minToT = 0;   // ns
  Float_t maxToT = 4.88;  // ns
  TFile * fileProf=0x0;
  if(strstr(optionSave,"save")){
    fileProf = new TFile("TOFCalibSave.root","recreate");
  }
  for (Int_t ich=0; ich<nch; ich++) {
    TH1F *hToT = new TH1F("htot","htot",nbinToT, minToT, maxToT);
    TH1F *hdeltaTime = new TH1F("hdeltaTime","hdeltaTime",200,2,4);
    Double_t binsProfile[101]; // sized larger than necessary, the correct 
    // dim being set in the booking of the profile
    TH2F * htimetot = new TH2F("htimetot","htimetot",nbinToT, minToT, maxToT,600,-5,10);
    Int_t ntracksTotal = 0;
    Int_t ntracksRun = 0;
    Int_t nusefulbins=0;
    Float_t meantime=0;
    Int_t i=-1;
    for (Int_t irun=0;irun<fNruns;irun++){
      i = ch[ich]+irun*fNChannels;
      AliDebug(2,Form("Calibrating channel %i",i));
      //fTree->GetEntry(i);
      fChain->GetEntry(i);
      ntracksTotal+=nentries/3;
    }
    if (ntracksTotal < MEANENTRIES) {
      AliInfo(Form(" Too small mean number of entires in channel %i (number of tracks = %d), not calibrating channel and continuing.....",i,ntracksTotal));
      continue;
    }
  
    for (Int_t irun=0;irun<fNruns;irun++){
      i = ch[ich]+irun*fNChannels;
      //fTree->GetEntry(i);
      fChain->GetEntry(i);
      ntracksRun=nentries/3;
      for (Int_t j=0;j<ntracksRun;j++){
	Int_t idxexToT = (j* NIDXSMALL)+DELTAIDXTOT; 
	Int_t idxexTime = (j* NIDXSMALL)+DELTAIDXTIME; 
	Int_t idxexExTime = (j* NIDXSMALL)+DELTAIDXPID; 
	Float_t tot = p[idxexToT];
	hdeltaTime->Fill(p[idxexTime]-p[idxexExTime]);
	meantime+=p[idxexTime]-p[idxexExTime];
	hToT->Fill(tot);
      }
    }

    nusefulbins = FindBins(hToT,&binsProfile[0]);
    meantime/=ntracksTotal;
    for (Int_t j=1;j<=nusefulbins;j++) {
      AliDebug(2,Form(" channel %i, nusefulbins = %i, bin %i = %f",i,nusefulbins,j,binsProfile[j])); 
    }

    TProfile* hSlewingProf = new TProfile("hSlewingProf", "hSlewingProf",nusefulbins, binsProfile, "G");  // CHECK THE BUILD OPTION, PLEASE!!!!!!
    for (Int_t irun=0;irun<fNruns;irun++){
      i = ch[ich]+irun*fNChannels;
      //fTree->GetEntry(i);
      fChain->GetEntry(i);
      ntracksRun=nentries/3;
      for (Int_t j=0;j<ntracksRun;j++){
	Int_t idxexToT = (j* NIDXSMALL)+DELTAIDXTOT; 
	Int_t idxexTime = (j* NIDXSMALL)+DELTAIDXTIME; 
	Int_t idxexExTime = (j* NIDXSMALL)+DELTAIDXPID; 
	Float_t tot = p[idxexToT];
	Float_t time = p[idxexTime]-p[idxexExTime];
	AliDebug(2,Form("track = %i, time = %f, tot = %f, time-meantime = %f",j,time, tot, time-meantime));
	hSlewingProf->Fill(tot,time);  // if meantime is not used, the fill may be moved in the loop above
	htimetot->Fill(tot,time-meantime);  // if meantime is not used, the fill may be moved in the loop above
      }
    }
    
    hSlewingProf->Fit("pol5",optionFit,"",1,4);
    TF1 * calibfunc = (TF1*)hSlewingProf->GetFunction("pol5");
    Float_t par[6];    
    for(Int_t kk=0;kk<6;kk++){
      par[kk]=calibfunc->GetParameter(kk);
      AliDebug(2,Form("parameter %i = %f",kk,par[kk]));
    }
    
    if(strstr(optionSave,"save") && fileProf){
      TString profName=Form("Profile%06i",i);
      TString timeTotName=Form("TimeTot%06i",i);
      TString totName=Form("Tot%06i",i);
      TString deltaName=Form("Delta%06i",i);
      fileProf->cd();
      hSlewingProf->Write(profName);
      htimetot->Write(timeTotName);
      hToT->Write(totName);
      hdeltaTime->Write(deltaName);
    }

    AliTOFChannelOffline * calChannel = (AliTOFChannelOffline*)fTOFCalOffline->At(i);
    calChannel->SetSlewPar(par);
    delete hToT;
    hToT=0x0;
    delete hSlewingProf;
    hSlewingProf=0x0;
    delete htimetot;
    htimetot=0x0;
    delete hdeltaTime;
    hdeltaTime=0x0;
  }

  if(strstr(optionSave,"save") && fileProf){
    fileProf->Close();
    delete fileProf;
    fileProf=0x0;
  }
  WriteParOfflineOnCDB("TOF/Calib","valid");

  return 0;
}
//----------------------------------------------------------------------------
Int_t AliTOFcalib::CalibrateFromProfile(Int_t ich, Option_t *optionSave, Option_t *optionFit){

  // computing calibration parameters using the old profiling algo
  // Returning codes:
  // 0 -> everything was ok
  // 1 -> no tree for calibration found
  // 2 -> not enough statistics to perform calibration
  // 3 -> problems with arrays

  TH1::AddDirectory(0);

  AliInfo(Form("*** Calibrating Histograms From Profile %s", GetName())) ; 
  AliInfo(Form("Option for Saving histos = %s",optionSave )) ; 
  AliInfo(Form("Option for Fitting Profile histos = %s",optionFit )) ; 
  Float_t p[MAXCHENTRIESSMALL];
  Int_t nentries;
  Int_t ntracksTotal=0;
  //fTree->SetBranchAddress("nentries",&nentries);
  //fTree->SetBranchAddress("TOFentries",p);
  fChain->SetBranchAddress("nentries",&nentries);
  fChain->SetBranchAddress("TOFentries",p);

  for (Int_t irun=0;irun<fNruns;irun++){
    Int_t i = ich+irun*fNChannels;
    //fTree->GetEntry(i);
    fChain->GetEntry(i);
    ntracksTotal+=nentries/3;
  }

  if (ntracksTotal < MEANENTRIES) {  
    AliInfo(Form(" Too small mean number of entires per channel (mean number = %d) not calibrating and exiting.....",ntracksTotal));
    return 2;
  }

  TH1F * hProf = Profile(ich);
  hProf->Fit("pol5",optionFit,"",0,4);
  TF1 * calibfunc = (TF1*)hProf->GetFunction("pol5");
  Float_t par[6];    
  for(Int_t kk=0;kk<6;kk++){
    par[kk]=calibfunc->GetParameter(kk);
    AliDebug(2,Form("parameter %i = %f",kk,par[kk]));
  }

  if(strstr(optionSave,"save")){
    TFile * fileProf = new TFile("TOFCalibSave.root","recreate");
    fileProf->cd(); 
    TString profName=Form("Profile%06i",ich);
    hProf->Write(profName);
    fileProf->Close();
    delete fileProf;
    fileProf=0x0;
  }

  delete hProf;
  hProf=0x0;
  AliTOFChannelOffline * calChannel = (AliTOFChannelOffline*)fTOFCalOffline->At(ich);
  calChannel->SetSlewPar(par);
  WriteParOfflineOnCDB("TOF/Calib","valid");
  return 0;
}
//----------------------------------------------------------------------------
Int_t AliTOFcalib::Calibrate(Option_t *optionSave, Option_t *optionFit){

  // calibrating the whole TOF
  // computing calibration parameters
  // Returning codes:
  // 0 -> everything was ok
  // 1 -> no tree for calibration found
  // 2 -> not enough statistics to perform calibration
  // 3 -> problems with arrays

  TH1::AddDirectory(0);

  AliInfo(Form("*** Calibrating Histograms %s, all channels", GetName())) ; 
  AliInfo(Form("Option for Saving histos = %s",optionSave )) ; 
  AliInfo(Form("Option for Fitting Profile histos = %s",optionFit )) ; 

  TFile * fileProf=0x0;
  if(strstr(optionSave,"save")){
    fileProf = new TFile("TOFCalibSave.root","recreate");
  }

  Float_t p[MAXCHENTRIESSMALL];
  Int_t nentries;
  //fTree->SetBranchAddress("nentries",&nentries);
  //fTree->SetBranchAddress("TOFentries",p);
  fChain->SetBranchAddress("nentries",&nentries);
  fChain->SetBranchAddress("TOFentries",p);

  Float_t ntracksTotalmean =0;
  for (Int_t ii=0; ii<fNChannels; ii++){
    for (Int_t irun=0;irun<fNruns;irun++){
      Int_t i = ii+irun*fNChannels;
      //fTree->GetEntry(i);
      fChain->GetEntry(i);
      ntracksTotalmean+=nentries/3;
    }
  }

  ntracksTotalmean/=fNChannels;
  if (ntracksTotalmean < MEANENTRIES) {
    AliInfo(Form(" Too small mean number of entires per channel (mean number = %f) not calibrating and exiting.....",ntracksTotalmean));
    return 2;
  }

  //filling ToT and Time arrays

  Int_t nbinToT = 100;  // ToT bin width in Profile = 50.0 ps 
  Float_t minToT = 0;   // ns
  Float_t maxToT = 4.88;// ns
  for (Int_t ii=0; ii<fNChannels; ii++) {
    TH1F *hToT = new TH1F("htot","htot",nbinToT, minToT, maxToT);
    TH1F *hdeltaTime = new TH1F("hdeltaTime","hdeltaTime",200,2,4);
    TH2F * htimetot = new TH2F("htimetot","htimetot",nbinToT, minToT, maxToT,600,-5,10);
    if (ii%1000 == 0) AliDebug(1,Form("Calibrating channel %i ",ii));
    //Int_t i = 3;
    Int_t nusefulbins=0;
    Double_t binsProfile[101]; // sized larger than necessary, the correct 
                              // dim being set in the booking of the profile
    Int_t ntracksRun = 0;
    Int_t ntracksTotal = 0;
    for (Int_t irun=0;irun<fNruns;irun++){
      Int_t i = ii+irun*fNChannels;
      //fTree->GetEntry(i);
      fChain->GetEntry(i);
      ntracksTotal+=nentries/3;
    }
    if (ntracksTotal < MEANENTRIES) {
      AliInfo(Form(" Too small mean number of entires in channel %i (number of tracks = %d), not calibrating channel and continuing.....",ii,ntracksTotal));
      continue;
    }
    Float_t meantime=0;
    for (Int_t irun=0;irun<fNruns;irun++){
      Int_t i = ii+irun*fNChannels;
      //fTree->GetEntry(i);
      fChain->GetEntry(i);
      ntracksRun=nentries/3;
      for (Int_t j=0;j<ntracksRun;j++){
	Int_t idxexToT = (j* NIDXSMALL)+DELTAIDXTOT; 
	Int_t idxexTime = (j* NIDXSMALL)+DELTAIDXTIME; 
	Int_t idxexExTime = (j* NIDXSMALL)+DELTAIDXPID; 
	Float_t tot = p[idxexToT];
	hdeltaTime->Fill(p[idxexTime]-p[idxexExTime]);
	meantime+=p[idxexTime]-p[idxexExTime];
	hToT->Fill(tot);
      }
    }
    nusefulbins = FindBins(hToT,&binsProfile[0]);
    meantime/=ntracksTotal;
    for (Int_t j=0;j<nusefulbins;j++) {
      AliDebug(2,Form(" channel %i, usefulbin = %i, low edge = %f",ii,j,binsProfile[j])); 
    }
    TProfile* hSlewingProf = new TProfile("hSlewingProf", "hSlewingProf",nusefulbins, binsProfile, "G");  // CHECK THE BUILD OPTION, PLEASE!!!!!!
    for (Int_t irun=0;irun<fNruns;irun++){
      Int_t i = ii+irun*fNChannels;
      //fTree->GetEntry(i);
      fChain->GetEntry(i);
      ntracksRun=nentries/3;
      for (Int_t j=0;j<ntracksRun;j++){
	Int_t idxexToT = (j* NIDXSMALL)+DELTAIDXTOT; 
	Int_t idxexTime = (j* NIDXSMALL)+DELTAIDXTIME; 
	Int_t idxexExTime = (j* NIDXSMALL)+DELTAIDXPID; 
	Float_t tot = p[idxexToT];
	Float_t time = p[idxexTime]-p[idxexExTime];
	AliDebug (2,Form("track = %i, time = %f, tot = %f, time-meantime = %f",j,time, tot, time-meantime));
	hSlewingProf->Fill(tot,time);  // if meantime is not used, the fill may be moved in the loop above
	htimetot->Fill(tot,time-meantime);  // if meantime is not used, the fill may be moved in the loop above
      }
    }
    hSlewingProf->Fit("pol5",optionFit,"",1,4);
    TF1 * calibfunc = (TF1*)hSlewingProf->GetFunction("pol5");
    Float_t par[6];    
    for(Int_t kk=0;kk<6;kk++){
      par[kk]=calibfunc->GetParameter(kk);
      AliDebug(2,Form("parameter %i = %f",kk,par[kk]));
    }

    if(strstr(optionSave,"save") && fileProf){
      TString profName=Form("Profile%06i",ii);
      TString timeTotName=Form("TimeTot%06i",ii);
      TString totName=Form("Tot%06i",ii);
      TString deltaName=Form("Delta%06i",ii);
      fileProf->cd();
      hSlewingProf->Write(profName);
      htimetot->Write(timeTotName);
      hToT->Write(totName);
      hdeltaTime->Write(deltaName);
    }
    AliTOFChannelOffline * calChannel = (AliTOFChannelOffline*)fTOFCalOffline->At(ii);
    calChannel->SetSlewPar(par);

    delete hToT;
    hToT=0x0;
    delete hSlewingProf;
    hSlewingProf=0x0;
    delete htimetot;
    htimetot=0x0;
    delete hdeltaTime;
    hdeltaTime=0x0;
  }

  if(strstr(optionSave,"save")){
    fileProf->Close();
    delete fileProf;
    fileProf=0x0;
  }
  WriteParOfflineOnCDB("TOF/Calib","valid");
  return 0;
}

//-----------------------------------------------------------------------
TH1F* AliTOFcalib::Profile(Int_t ich)
{
  // profiling algo

  Float_t p[MAXCHENTRIESSMALL];
  Int_t nentries;
  //fTree->SetBranchAddress("nentries",&nentries);
  //fTree->SetBranchAddress("TOFentries",p);
  fChain->SetBranchAddress("nentries",&nentries);
  fChain->SetBranchAddress("TOFentries",p);

  //Prepare histograms for Slewing Correction
  const Int_t knbinToT = 100;
  Int_t nbinTime = 200;
  Float_t minTime = -5.5; //ns
  Float_t maxTime = 5.5; //ns
  Float_t minToT = 0; //ns
  Float_t maxToT = 5.; //ns
  Float_t deltaToT = (maxToT-minToT)/knbinToT;
  Double_t mTime[knbinToT+1],mToT[knbinToT+1],meanTime[knbinToT+1], meanTime2[knbinToT+1],vToT[knbinToT+1], vToT2[knbinToT+1],meanToT[knbinToT+1],meanToT2[knbinToT+1],vTime[knbinToT+1],vTime2[knbinToT+1],xlow[knbinToT+1],sigmaTime[knbinToT+1];
  Int_t n[knbinToT+1], nentrx[knbinToT+1];
  Double_t sigmaToT[knbinToT+1];
  for (Int_t i = 0; i < knbinToT+1 ; i++){
    mTime[i]=0;
    mToT[i]=0;
    n[i]=0;
    meanTime[i]=0;
    meanTime2[i]=0;
    vToT[i]=0;
    vToT2[i]=0;
    meanToT[i]=0;
    meanToT2[i]=0;
    vTime[i]=0;
    vTime2[i]=0;
    xlow[i]=0;
    sigmaTime[i]=0;
    sigmaToT[i]=0;
    n[i]=0;
    nentrx[i]=0;
  }
  TH2F* hSlewing = new TH2F("hSlewing", "hSlewing", knbinToT, minToT, maxToT, nbinTime, minTime, maxTime);
  Int_t ntracksRun = 0;
  TH1F *histo = new TH1F("histo", "1D Time vs ToT", knbinToT, minToT, maxToT);
  for (Int_t irun=0;irun<fNruns;irun++){
    Int_t i = ich+irun*fNChannels;
    //fTree->GetEntry(i);
    fChain->GetEntry(i);
    ntracksRun=nentries/3;
    for (Int_t j=0;j<ntracksRun;j++){
      Int_t idxexToT = (j* NIDXSMALL)+DELTAIDXTOT; 
      Int_t idxexTime = (j* NIDXSMALL)+DELTAIDXTIME; 
      Int_t idxexExTime = (j* NIDXSMALL)+DELTAIDXPID; 
      Float_t tot = p[idxexToT];
      Float_t time = p[idxexTime]-p[idxexExTime];
      Int_t nx = (Int_t)((tot-minToT)/deltaToT)+1;
      if ((tot != 0) && ( time!= 0)){
	vTime[nx]+=time;
	vTime2[nx]+=time*time;
	vToT[nx]+=tot;
	vToT2[nx]+=tot*tot;
	nentrx[nx]++;
	hSlewing->Fill(tot,time);
      }
    }
  }
  Int_t nbinsToT=hSlewing->GetNbinsX();
  if (nbinsToT != knbinToT) {
    AliError("Profile :: incompatible numbers of bins");
    return 0x0;
  }
  
  Int_t usefulBins=0;
  for (Int_t i=1;i<=nbinsToT;i++){
    if (nentrx[i]!=0){
      n[usefulBins]+=nentrx[i];
      if (n[usefulBins]==0 && i == nbinsToT) {
	break;
      }
      meanTime[usefulBins]+=vTime[i];
      meanTime2[usefulBins]+=vTime2[i];
      meanToT[usefulBins]+=vToT[i];
      meanToT2[usefulBins]+=vToT2[i];
      if (n[usefulBins]<10 && i!=nbinsToT) continue; 
      mTime[usefulBins]=meanTime[usefulBins]/n[usefulBins];
      mToT[usefulBins]=meanToT[usefulBins]/n[usefulBins];
      sigmaTime[usefulBins]=TMath::Sqrt(1./n[usefulBins]/n[usefulBins]
			    *(meanTime2[usefulBins]-meanTime[usefulBins]
			    *meanTime[usefulBins]/n[usefulBins]));
      if ((1./n[usefulBins]/n[usefulBins]
	   *(meanToT2[usefulBins]-meanToT[usefulBins]
	     *meanToT[usefulBins]/n[usefulBins]))< 0) {
	AliError(" too small radical" );
	sigmaToT[usefulBins]=0;
      }
      else{       
	sigmaToT[usefulBins]=TMath::Sqrt(1./n[usefulBins]/n[usefulBins]
			     *(meanToT2[usefulBins]-meanToT[usefulBins]
			     *meanToT[usefulBins]/n[usefulBins]));
      }
      usefulBins++;
    }
  }
  for (Int_t i=0;i<usefulBins;i++){
    Int_t binN = (Int_t)((mToT[i]-minToT)/deltaToT)+1;
    histo->Fill(mToT[i],mTime[i]);
    histo->SetBinError(binN,sigmaTime[i]);
  } 
  delete hSlewing;
  hSlewing=0x0;

  return histo;
}
//----------------------------------------------------------------------------
Int_t AliTOFcalib::FindBins(TH1F* h, Double_t *binsProfile) const{

  // to determine the bins for ToT histo

  Int_t cont = 0;
  Int_t startBin = 1;
  Int_t nbin = h->GetNbinsX();
  Int_t nentries = (Int_t)h->GetEntries();
  Float_t max = h->GetBinLowEdge(nbin);
  Int_t nusefulbins=0;
  Int_t maxcont=0;
  // setting maxvalue of entries per bin
  if (nentries <= 60) maxcont = 2;
  else  if (nentries <= 100) maxcont = 5;
  else  if (nentries <= 500) maxcont = 10;
  else  maxcont = 20;
  for (Int_t j=1;j<=nbin;j++) {
    cont += (Int_t)h->GetBinContent(j);
    if (j<nbin){
      if (cont>=maxcont){
	nusefulbins++;
	binsProfile[nusefulbins-1]=h->GetBinLowEdge(startBin);
	cont=0;
	startBin=j+1;
	continue;
      }
    }
    else{
      if (cont>=maxcont){
	nusefulbins++;
	binsProfile[nusefulbins-1]=h->GetBinLowEdge(startBin);
	binsProfile[nusefulbins]=max;
      }
      else {
	binsProfile[nusefulbins]=h->GetBinLowEdge(startBin);
      }
    }
  }
  return nusefulbins;
}


//----------------------------------------------------------------------------

void
AliTOFcalib::CreateDeltaBCOffset()
{
  /*
   * create deltaBC offset
   */

  if (fDeltaBCOffset) {
    AliWarning("DeltaBCOffset object already defined, cannot create a new one");
    return;
  }
  fDeltaBCOffset = new AliTOFDeltaBCOffset();
}
  
//----------------------------------------------------------------------------

void
AliTOFcalib::CreateCTPLatency()
{
  /*
   * create CTP latency
   */

  if (fCTPLatency) {
    AliWarning("CTPLatency object already defined, cannot create a new one");
    return;
  }
  fCTPLatency = new AliTOFCTPLatency();
}
  
//----------------------------------------------------------------------------

void
AliTOFcalib::CreateT0Fill()
{
  /*
   * create event-time
   */

  if (fT0Fill) {
    AliWarning("T0Fill object already defined, cannot create a new one");
    return;
  }
  fT0Fill = new AliTOFT0Fill();
}
  
//----------------------------------------------------------------------------

void
AliTOFcalib::CreateRunParams()
{
  /*
   * create run params
   */

  if (fRunParams) {
    AliWarning("RunParams object already defined, cannot create a new one");
    return;
  }
  fRunParams = new AliTOFRunParams();
}
  
//----------------------------------------------------------------------------

void
AliTOFcalib::WriteDeltaBCOffsetOnCDB(const Char_t *sel , Int_t minrun, Int_t maxrun)
{
  /*
   * deltaBC offset on CDB 
   */
  
  if (!fDeltaBCOffset) return;
  AliCDBId id(Form("%s/DeltaBCOffset", sel), minrun, maxrun);
  AliCDBMetaData *md = new AliCDBMetaData();
  md->SetResponsible("Roberto Preghenella");
  AliCDBManager *man = AliCDBManager::Instance();
  man->Put(fDeltaBCOffset, id, md);
  AliDebug(2,Form("DeltaBCOffset written on CDB with run range [%i, %i] ",minrun ,maxrun));
  delete md;
}

//----------------------------------------------------------------------------

void
AliTOFcalib::WriteCTPLatencyOnCDB(const Char_t *sel , Int_t minrun, Int_t maxrun)
{
  /*
   * write CTP latency on CDB 
   */
  
  if (!fCTPLatency) return;
  AliCDBId id(Form("%s/CTPLatency", sel), minrun, maxrun);
  AliCDBMetaData *md = new AliCDBMetaData();
  md->SetResponsible("Roberto Preghenella");
  AliCDBManager *man = AliCDBManager::Instance();
  man->Put(fCTPLatency, id, md);
  AliDebug(2,Form("CTPLatency written on CDB with run range [%i, %i] ",minrun ,maxrun));
  delete md;
}

//----------------------------------------------------------------------------

void
AliTOFcalib::WriteT0FillOnCDB(const Char_t *sel , Int_t minrun, Int_t maxrun)
{
  /*
   * write event-time on CDB 
   */
  
  if (!fT0Fill) return;
  AliCDBId id(Form("%s/T0Fill", sel), minrun, maxrun);
  AliCDBMetaData *md = new AliCDBMetaData();
  md->SetResponsible("Roberto Preghenella");
  AliCDBManager *man = AliCDBManager::Instance();
  man->Put(fT0Fill, id, md);
  AliDebug(2,Form("T0Fill written on CDB with run range [%i, %i] ",minrun ,maxrun));
  delete md;
}

//----------------------------------------------------------------------------

void
AliTOFcalib::WriteRunParamsOnCDB(const Char_t *sel , Int_t minrun, Int_t maxrun)
{
  /*
   * write run params on CDB 
   */
  
  if (!fRunParams) return;
  AliCDBId id(Form("%s/RunParams", sel), minrun, maxrun);
  AliCDBMetaData *md = new AliCDBMetaData();
  md->SetResponsible("Roberto Preghenella");
  AliCDBManager *man = AliCDBManager::Instance();
  man->Put(fRunParams, id, md);
  AliDebug(2,Form("RunParams written on CDB with run range [%i, %i] ",minrun ,maxrun));
  delete md;
}

//----------------------------------------------------------------------------

void
AliTOFcalib::WriteReadoutEfficiencyOnCDB(const Char_t *sel , Int_t minrun, Int_t maxrun)
{
  /*
   * write readout efficiency on CDB 
   */
  
  if (!fReadoutEfficiency) return;
  AliCDBId id(Form("%s/ReadoutEfficiency", sel), minrun, maxrun);
  AliCDBMetaData *md = new AliCDBMetaData();
  md->SetResponsible("Roberto Preghenella");
  AliCDBManager *man = AliCDBManager::Instance();
  man->Put(fReadoutEfficiency, id, md);
  AliDebug(2,Form("ReadoutEfficiency written on CDB with run range [%i, %i] ",minrun ,maxrun));
  delete md;
}

//----------------------------------------------------------------------------

void
AliTOFcalib::WriteProblematicOnCDB(const Char_t *sel , Int_t minrun, Int_t maxrun)
{
  /*
   * write problematic on CDB 
   */
  
  if (!fProblematic) return;
  AliCDBId id(Form("%s/Problematic", sel), minrun, maxrun);
  AliCDBMetaData *md = new AliCDBMetaData();
  md->SetResponsible("Roberto Preghenella");
  AliCDBManager *man = AliCDBManager::Instance();
  man->Put(fProblematic, id, md);
  AliDebug(2,Form("Problematic written on CDB with run range [%i, %i] ",minrun ,maxrun));
  delete md;
}

//----------------------------------------------------------------------------

Bool_t
AliTOFcalib::ReadDeltaBCOffsetFromCDB(const Char_t *sel , Int_t nrun)
{
  /*
   * read deltaBC offset from CDB
   */
  
  AliCDBManager *man = AliCDBManager::Instance();
  AliCDBEntry *entry = man->Get(Form("%s/DeltaBCOffset", sel),nrun);
  if (!entry) { 
    AliFatal("No DeltaBCOffset entry found in CDB");
    exit(0);  
  }
  fDeltaBCOffset =(AliTOFDeltaBCOffset *)entry->GetObject();
  if(!fDeltaBCOffset){
    AliFatal("No DeltaBCOffset object found in CDB entry");
    exit(0);  
  }  
  return kTRUE; 
}

//----------------------------------------------------------------------------

Bool_t
AliTOFcalib::ReadCTPLatencyFromCDB(const Char_t *sel , Int_t nrun)
{
  /*
   * read CTP latency from CDB
   */
  
  AliCDBManager *man = AliCDBManager::Instance();
  AliCDBEntry *entry = man->Get(Form("%s/CTPLatency", sel),nrun);
  if (!entry) { 
    AliFatal("No CTPLatency entry found in CDB");
    exit(0);  
  }
  fCTPLatency =(AliTOFCTPLatency *)entry->GetObject();
  if(!fCTPLatency){
    AliFatal("No CTPLatency object found in CDB entry");
    exit(0);  
  }  
  return kTRUE; 
}

//----------------------------------------------------------------------------

Bool_t
AliTOFcalib::ReadT0FillFromCDB(const Char_t *sel , Int_t nrun)
{
  /*
   * read event-time from CDB
   */
  
  AliCDBManager *man = AliCDBManager::Instance();
  AliCDBEntry *entry = man->Get(Form("%s/T0Fill", sel),nrun);
  if (!entry) { 
    AliFatal("No T0Fill entry found in CDB");
    exit(0);  
  }
  fT0Fill =(AliTOFT0Fill *)entry->GetObject();
  if(!fT0Fill){
    AliFatal("No T0Fill object found in CDB entry");
    exit(0);  
  }  
  return kTRUE; 
}

//----------------------------------------------------------------------------

Bool_t
AliTOFcalib::ReadRunParamsFromCDB(const Char_t *sel , Int_t nrun)
{
  /*
   * read run params from CDB
   */
  
  AliCDBManager *man = AliCDBManager::Instance();
  AliCDBEntry *entry = man->Get(Form("%s/RunParams", sel),nrun);
  if (!entry) { 
    AliFatal("No RunParams entry found in CDB");
    exit(0);  
  }
  fRunParams =(AliTOFRunParams *)entry->GetObject();
  if(!fRunParams){
    AliFatal("No RunParams object found in CDB entry");
    exit(0);  
  }  
  return kTRUE; 
}

//----------------------------------------------------------------------------

Bool_t
AliTOFcalib::ReadLHCClockPhaseFromCDB(const Char_t *sel , Int_t nrun)
{
  /*
   * read LHC clock-phase from CDB
   */
  
  AliCDBManager *man = AliCDBManager::Instance();
  AliCDBEntry *entry = man->Get(Form("%s/LHCClockPhase", sel),nrun);
  if (!entry) { 
    AliFatal("No LHCClockPhase entry found in CDB");
    exit(0);  
  }
  fLHCClockPhase =(AliLHCClockPhase *)entry->GetObject();
  if(!fRunParams){
    AliFatal("No LHCClockPhase object found in CDB entry");
    exit(0);  
  }  
  return kTRUE; 
}

//----------------------------------------------------------------------------

Bool_t
AliTOFcalib::ReadReadoutEfficiencyFromCDB(const Char_t *sel , Int_t nrun)
{
  /*
   * read readout efficiency from CDB
   */
  
  AliCDBManager *man = AliCDBManager::Instance();
  AliCDBEntry *entry = man->Get(Form("%s/ReadoutEfficiency", sel),nrun);
  if (!entry) { 
    AliFatal("No ReadoutEfficiency entry found in CDB");
    exit(0);  
  }
  fReadoutEfficiency = (TH1F *)entry->GetObject();
  if(!fReadoutEfficiency){
    AliFatal("No ReadoutEfficiency object found in CDB entry");
    exit(0);  
  }  
  return kTRUE; 
}

//----------------------------------------------------------------------------

Bool_t
AliTOFcalib::ReadProblematicFromCDB(const Char_t *sel , Int_t nrun)
{
  /*
   * read problematic from CDB
   */
  
  AliCDBManager *man = AliCDBManager::Instance();
  AliCDBEntry *entry = man->Get(Form("%s/Problematic", sel),nrun);
  if (!entry) { 
    AliFatal("No Problematic entry found in CDB");
    exit(0);  
  }
  fProblematic = (TH1C *)entry->GetObject();
  if(!fProblematic){
    AliFatal("No Problematic object found in CDB entry");
    exit(0);  
  }  
  return kTRUE; 
}

//----------------------------------------------------------------------------

Bool_t 
AliTOFcalib::Init(Int_t run)
{
  /*
   * init
   */

  if (fInitFlag) {
    AliWarning("the class was already initialized, re-initialize it");
    fInitFlag = kFALSE;
  }
  
  /* read channel status array */
  if (!ReadParOnlineStatusFromCDB("TOF/Calib", run)) {
    AliError("cannot get \"Status\" object from OCDB");
    return kFALSE;
  }
  /* get par offline array */
  if (!ReadParOfflineFromCDB("TOF/Calib", run)) {
    AliError("cannot get \"ParOffline\" object from OCDB");
    return kFALSE;
  }
  /* get deltaBC offset obj */
  if (!ReadDeltaBCOffsetFromCDB("TOF/Calib", run)) {
    AliError("cannot get \"DeltaBCOffset\" object from OCDB");
    return kFALSE;
  }
  /* get CTP latency obj */
  if (!ReadCTPLatencyFromCDB("TOF/Calib", run)) {
    AliError("cannot get \"CTPLatency\" object from OCDB");
    return kFALSE;
  }
  /* get run params obj */
  if (!ReadRunParamsFromCDB("TOF/Calib", run)) {
    AliError("cannot get \"RunParams\" object from OCDB");
    return kFALSE;
  }
  /* get LHC clock-phase obj */
  if (!ReadLHCClockPhaseFromCDB("GRP/Calib", run)) {
    AliError("cannot get \"LHCClockPhase\" object from OCDB");
    return kFALSE;
  }
  /* get readout efficiency obj */
  if (!ReadReadoutEfficiencyFromCDB("TOF/Calib", run)) {
    AliError("cannot get \"ReadoutEfficiency\" object from OCDB");
    return kFALSE;
  }
  /* get readout efficiency obj */
  if (!ReadProblematicFromCDB("TOF/Calib", run)) {
    AliError("cannot get \"Problematic\" object from OCDB");
    return kFALSE;
  }
  /* get response params */
  TFile *responseFile = TFile::Open("$ALICE_ROOT/TOF/data/AliTOFresponsePar.root");
  if (!responseFile || !responseFile->IsOpen()) {
    AliError("cannot open \"ResponseParams\" local file");
    return kFALSE;
  }
  fResponseParams = (AliTOFResponseParams *)responseFile->Get("ResponseParams");
  if (!fResponseParams) {
    AliError("cannot get \"ResponseParams\" object from local file");
    return kFALSE;
  }
  responseFile->Close();

  /* check whether to use the clock phase */
  if (fRunParams->GetUseLHCClockPhase())
    fUseLHCClockPhase = kTRUE;
 
  if (fUseLHCClockPhase)
    AliInfo("calibration using BPTX LHC clock-phase");

  /* all done */
  fInitFlag = kTRUE;
  return kTRUE;

}

//----------------------------------------------------------------------------

Double_t
AliTOFcalib::GetTimeCorrection(Int_t index, Double_t tot, Int_t deltaBC, Int_t l0l1, UInt_t timestamp)
{
  /*
   * get time correction
   */

  if (!fInitFlag) {
    AliError("class not yet initialized. Initialize it before.");
    return 0.;
  }

  /* deal with L0-L1 orbit crossing (negative values) */
  if (l0l1 < 0) l0l1 += 3564;

  /* get calibration params */
  AliTOFChannelOffline *parOffline = (AliTOFChannelOffline *)fTOFCalOffline->At(index);
  Int_t deltaBCOffset = fDeltaBCOffset->GetDeltaBCOffset();
  Float_t ctpLatency = fCTPLatency->GetCTPLatency();
  Float_t tdcLatencyWindow = fStatus->GetLatencyWindow(index) * 1.e3;
  Float_t timezero = fRunParams->EvalT0(timestamp);
  Float_t clockphase = fLHCClockPhase->GetPhase(timestamp);
  /* check whether to remove mean T0.
   * useful when one wants to compute mean T0 */
  if (!fRemoveMeanT0) timezero = 0.;
  /* check whether to use the clock phase */
  if (fUseLHCClockPhase) timezero -= 1.e3 * clockphase;

  /* compute correction */
  Double_t corr = 0.;
  /* deltaBC correction */
  deltaBC = deltaBCOffset; /* inhibit deltaBC correction for the time being */
  corr += (deltaBC - deltaBCOffset) * AliTOFGeometry::BunchCrossingBinWidth();
  /* L0-L1 latency correction */
  corr -= l0l1 * AliTOFGeometry::BunchCrossingBinWidth();
  /* CTP latency correction */
  corr -= ctpLatency;
  /* TDC latency window correction */
  corr += tdcLatencyWindow;
  /* time-zero correction */
  corr += timezero;
  /* time calibration correction */
  if (tot < AliTOFGeometry::SlewTOTMin()) 
    tot = AliTOFGeometry::SlewTOTMin();
  if (tot > AliTOFGeometry::SlewTOTMax()) 
    tot = AliTOFGeometry::SlewTOTMax();
  for (Int_t islew = 0; islew < 6; islew++)
    corr += parOffline->GetSlewPar(islew) * TMath::Power(tot, islew) * 1.e3;

  /* return correction */
  return corr;
}

//----------------------------------------------------------------------------

void
AliTOFcalib::CalibrateESD(AliESDEvent *event)
{
  /*
   * calibrate ESD
   */

  if (!fInitFlag) {
    AliError("class not yet initialized. Initialize it before.");
    return;
  }

  /* loop over tracks */
  AliESDtrack *track = NULL;
  Int_t index, l0l1, deltaBC;
  Double_t time, tot, corr, texp[AliPID::kSPECIES];
  UInt_t timestamp = event->GetTimeStamp();
  for (Int_t itrk = 0; itrk < event->GetNumberOfTracks(); itrk++) {

    /* get track */
    track = event->GetTrack(itrk);
    if (!track || !(track->GetStatus() & AliESDtrack::kTOFout)) continue;

    /* calibrate TOF signal */
    if (fCalibrateTOFsignal) {
      /* get info */
      index = track->GetTOFCalChannel();
      time = track->GetTOFsignalRaw();
      tot = track->GetTOFsignalToT();
      l0l1 = track->GetTOFL0L1();
      deltaBC = track->GetTOFDeltaBC();
      /* get correction */
      corr = GetTimeCorrection(index, tot, deltaBC, l0l1, timestamp);
      /* apply correction */
      time -= corr;
      /* set new TOF signal */
      track->SetTOFsignal(time);
    }

    /* correct expected time */
    if (fCorrectTExp) {
      /* get integrated times */
      track->GetIntegratedTimes(texp);
      /* loop over particle types and correct expected time */
      for (Int_t ipart = 0; ipart < AliPID::kSPECIES; ipart++)
	texp[ipart] += fResponseParams->EvalTExpCorr(ipart, track->P());
      /* set integrated times */
      track->SetIntegratedTimes(texp);
    }

  }

}

//----------------------------------------------------------------------------

Bool_t
AliTOFcalib::IsChannelEnabled(Int_t index, Bool_t checkEfficiency, Bool_t checkProblematic)
{
  /*
   * is channel enabled
   */

  if (!fInitFlag) {
    AliError("class not yet initialized. Initialize it before.");
    return kTRUE;
  }

  /* check bad status */
  if (fStatus->GetPulserStatus(index) == AliTOFChannelOnlineStatusArray::kTOFPulserBad) return kFALSE;
  if (fStatus->GetNoiseStatus(index) == AliTOFChannelOnlineStatusArray::kTOFNoiseBad) return kFALSE;
  if (fStatus->GetHWStatus(index) == AliTOFChannelOnlineStatusArray::kTOFHWBad) return kFALSE;
  if (checkEfficiency && !IsChannelEfficient(index)) return kFALSE;
  if (checkProblematic && IsChannelProblematic(index)) return kFALSE;
  
  /* good status */
  return kTRUE;

}

//----------------------------------------------------------------------------

Bool_t
AliTOFcalib::IsChannelEfficient(Int_t index)
{
  /*
   * is channel efficient
   */

  if (!fInitFlag) {
    AliError("class not yet initialized. Initialize it before.");
    return kTRUE;
  }

  /* check efficiency */
  if (fReadoutEfficiency->GetBinContent(index + 1) < 0.95) return kFALSE;
  return kTRUE;

}

//----------------------------------------------------------------------------

Bool_t
AliTOFcalib::IsChannelProblematic(Int_t index)
{
  /*
   * is channel problematic
   */

  if (!fInitFlag) {
    AliError("class not yet initialized. Initialize it before.");
    return kTRUE;
  }

  /* check problematic */
  if (fProblematic->GetBinContent(index + 1) != 0) return kTRUE;
  return kFALSE;

}

//----------------------------------------------------------------------------

void
AliTOFcalib::CalibrateTExp(AliESDEvent *event) const
{
  /*
   * calibrate TExp
   */

  if (!fInitFlag) {
    AliError("class not yet initialized. Initialize it before.");
    return;
  }

  /* loop over tracks */
  AliESDtrack *track = NULL;
  Double_t texp[AliPID::kSPECIES];
  for (Int_t itrk = 0; itrk < event->GetNumberOfTracks(); itrk++) {

    /* get track */
    track = event->GetTrack(itrk);
    if (!track || !(track->GetStatus() & AliESDtrack::kTOFout)) continue;

    /* get integrated times */
    track->GetIntegratedTimes(texp);
    /* loop over particle types and correct expected time */
    for (Int_t ipart = 0; ipart < AliPID::kSPECIES; ipart++)
      texp[ipart] += fResponseParams->EvalTExpCorr(ipart, track->P());
    /* set integrated times */
    track->SetIntegratedTimes(texp);

  }

}

//----------------------------------------------------------------------------

Double_t
AliTOFcalib::TuneForMC(AliESDEvent *event, Double_t resolution)
{
  /*
   * tune for MC
   */

  /* get vertex spread and define T0-spread */
  Double_t diamond2 = TMath::Abs(event->GetSigma2DiamondZ());
  Double_t t0spread = TMath::Sqrt(diamond2) / 2.99792457999999984e-02;
  /* generate random startTime */
  Double_t startTime = gRandom->Gaus(0., t0spread);
  /* define extra smearing for resolution */
  Double_t defaultResolution = 80.;
  Double_t extraSmearing = 0.;
  if (resolution > defaultResolution)
    extraSmearing = TMath::Sqrt(resolution * resolution - defaultResolution * defaultResolution);

  /* loop over tracks */
  AliESDtrack *track = NULL;
  Double_t time;
  for (Int_t itrk = 0; itrk < event->GetNumberOfTracks(); itrk++) {
    /* get track */
    track = event->GetTrack(itrk);
    if (!track) continue;
    /* check TOF match */
    if (!track->IsOn(AliESDtrack::kTOFout)) continue;
    /* check if channel is enabled */
    if (!IsChannelEnabled(track->GetTOFCalChannel())) {
      /* reset TOF status */
      track->ResetStatus(AliESDtrack::kTOFin);
      track->ResetStatus(AliESDtrack::kTOFout);
      track->ResetStatus(AliESDtrack::kTOFmismatch);
      track->ResetStatus(AliESDtrack::kTOFpid);
    }
    /* get original time and manipulate it */
    time = track->GetTOFsignal();
    time += startTime; /* add start time */
    time += gRandom->Gaus(0., extraSmearing); /* extra smearing */
    time -= 25.; /* remove 25 ps to center the signal */
    track->SetTOFsignal(time);
  }

  return startTime;
}
