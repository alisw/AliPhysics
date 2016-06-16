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
  Responsible: marian.ivanov@cern.ch 
  Code to analyze the TPC calibration and to produce OCDB entries  


   .x ~/rootlogon.C
   gSystem->Load("libANALYSIS");
   gSystem->Load("libTPCcalib");

   AliTPCPreprocessorOffline proces;
   TString ocdbPath="local:////"
   ocdbPath+=gSystem->GetFromPipe("pwd");

   proces.CalibTimeGain("CalibObjects.root",run0,run1,ocdbPath);
   proces.CalibTimeVdrift("CalibObjects.root",run0,run1,ocdbPath);
  // take the raw calibration data from the file CalibObjects.root 
  // and make a OCDB entry with run  validity run0-run1
  // results are stored at the ocdbPath - local or alien ...
  // default storage ""- data stored at current working directory 
 
  e.g.
  gSystem->Load("libANALYSIS");
  gSystem->Load("libTPCcalib");
  AliTPCPreprocessorOffline proces;
  proces.CalibTimeGain("TPCMultObjects.root",114000,140040,0);
  TFile oo("OCDB/TPC/Calib/TimeGain/Run114000_121040_v0_s0.root")
 TObjArray * arr = AliCDBEntry->GetObject()
  arr->At(4)->Draw("alp")

*/
#include "Riostream.h"
#include <fstream>
#include "TMap.h"
#include "TGraphErrors.h"
#include "AliExternalTrackParam.h"
#include "TROOT.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TCanvas.h"
#include "THnSparse.h"
#include "TLegend.h"
#include "TPad.h"
#include "TH2D.h"
#include "TH3D.h"
#include "AliTPCROC.h"
#include "AliTPCCalROC.h"
//#include "AliESDfriend.h"
#include "AliTPCcalibTime.h"
#include "AliSplineFit.h"
#include "AliCDBMetaData.h"
#include "AliCDBId.h"
#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#include "AliTPCcalibBase.h"
#include "AliTPCcalibDB.h"
#include "AliTPCcalibDButil.h"
#include "AliRelAlignerKalman.h"
#include "AliTPCParamSR.h"
#include "AliTPCcalibTimeGain.h"
#include "AliTPCcalibGainMult.h"
#include "AliTPCcalibAlign.h"
#include "AliSplineFit.h"
#include "AliTPCComposedCorrection.h"
#include "AliTPCExBTwist.h"
#include "AliTPCCalibGlobalMisalignment.h"
#include "TStatToolkit.h"
#include "TChain.h"
#include "TCut.h"
#include "AliTrackerBase.h"
#include "AliTracker.h"
#include "AliTPCPreprocessorOffline.h"
#include "AliTPCCorrectionFit.h"
#include "AliCDBEntry.h"

#include "AliTPCClusterParam.h"
#include "AliTPCRecoParam.h"

using std::endl;
using std::cout;

ClassImp(AliTPCPreprocessorOffline)
//_____________________________________________________________________________
AliTPCPreprocessorOffline::AliTPCPreprocessorOffline():
  TNamed("TPCPreprocessorOffline","TPCPreprocessorOffline"),
  fNormaliseQA(kTRUE),
  fGainCalibrationType(kFullGainCalib),  // gain calibration type
  fMinEntries(500),                      // minimal number of entries for fit
  fStartRun(0),                         // start Run - used to make fast selection in THnSparse
  fEndRun(0),                           // end   Run - used to make fast selection in THnSparse
  fStartTime(0),                        // fStartTime - used to make fast selection in THnSparse
  fEndTime(0),                          // fEndTime   - used to make fast selection in THnSparse
  fOCDBstorage(0),                       // OCDB storage
  fVdriftArray(new TObjArray),
  fTimeDrift(0),
  fGraphMIP(0),                // graph time dependence of MIP
  fGraphCosmic(0),             // graph time dependence at Plateu
  fGraphAttachmentMIP(0),
  fFitMIP(0),                  // fit of dependence - MIP
  fFitCosmic(0),               // fit of dependence - Plateu
  fGainArray(new TObjArray),               // array to be stored in the OCDB
  fGainArrayCombined(0x0),
  fArrQAhist(0x0),
  fGainMIP(0),          // calibration component for MIP
  fGainCosmic(0),       // calibration component for cosmic
  fGainMult(0),
  fAlignTree(0),        // alignment tree
  fSwitchOnValidation(kFALSE), // flag to switch on validation of OCDB parameters
  fMinGain(1.5),
  fMaxGain(4.5),
  fMaxVdriftCorr(0.03),
  fNtracksVdrift(0),
  fMinTracksVdrift(0),
  fNeventsVdrift(0),
  fMinEventsVdrift(0),
  fCalibrationStatus(0),
  fDriftCDBentry(NULL)
{
  //
  // default constructor
  //
}

//_____________________________________________________________________________
AliTPCPreprocessorOffline::~AliTPCPreprocessorOffline() {
  //
  // Destructor
  //
  delete fDriftCDBentry;
}

//_____________________________________________________________________________
void AliTPCPreprocessorOffline::GetRunRange(AliTPCcalibTime * const  timeDrift){
  //
  // find the fist and last run
  //
  TObjArray *hisArray =timeDrift->GetHistoDrift();
  {for (Int_t i=0; i<hisArray->GetEntriesFast(); i++){
    THnSparse* addHist=(THnSparse*)hisArray->UncheckedAt(i);
    if (!addHist) continue;
    if (addHist->GetEntries()<fMinEntries) continue;
    TH1D* histo    =addHist->Projection(3);
    TH1D* histoTime=addHist->Projection(0);
    AliInfo(Form("%s\t%f\t%d\t%d",histo->GetName(), histo->GetEntries(),histo->FindFirstBinAbove(0),histo->FindLastBinAbove(0)));

    if (fStartRun<=0){ 
      fStartRun=histo->FindFirstBinAbove(0);
      fEndRun  =histo->FindLastBinAbove(0);
    }else{
      fStartRun=TMath::Min(histo->FindFirstBinAbove(0),fStartRun);
      fEndRun  =TMath::Max(histo->FindLastBinAbove(0),fEndRun);
    }
    if (fStartTime==0){ 
      fStartTime=histoTime->FindFirstBinAbove(0);
      fEndTime  =histoTime->FindLastBinAbove(0);
    }else{
      fStartTime=TMath::Min(histoTime->FindFirstBinAbove(0),fStartTime);
      fEndTime  =TMath::Max(histoTime->FindLastBinAbove(0),fEndTime);
    }
    delete histo;
    delete histoTime;
  }}
  if (fStartRun<0) fStartRun=0;
  if (fEndRun<0) fEndRun=100000000;
  AliInfo(Form("Run range  :\t%d-%d", fStartRun, fEndRun));
  AliInfo(Form("Time range :\t%d-%d", fStartTime, fEndTime));

}

//_____________________________________________________________________________
Int_t AliTPCPreprocessorOffline::CalibTimeVdrift(AliTPCcalibTime* timeDrift, Int_t ustartRun, Int_t uendRun)
{
  // make calibration of the drift velocity
  // Input parameters:
  //      timeDrift              - the calibration object
  //      ustartRun, uendRun     - run validity period 
  // return 0 on success, 1 on failure

  fTimeDrift=timeDrift;
  fStartRun=ustartRun;
  fEndRun=ustartRun; 
  GetRunRange(fTimeDrift);
  
  //extract statistics
  fNtracksVdrift = TMath::Nint(fTimeDrift->GetResHistoTPCITS(0)->GetEntries());
  //if we have 0 ITS TPC matches it means we have no ITS tracks and we try to use TPC-TOF matching for calibration
  if (fNtracksVdrift==0) fNtracksVdrift=TMath::Nint(fTimeDrift->GetResHistoTPCTOF(0)->GetEntries());
  fNeventsVdrift = TMath::Nint(fTimeDrift->GetTPCVertexHisto(0)->GetEntries());

  TObjArray *hisArray =fTimeDrift->GetHistoDrift();  
  for (Int_t i=0; i<hisArray->GetEntriesFast(); i++){
    THnSparse* addHist=(THnSparse*)hisArray->At(i);
    if (!addHist) continue;
    if (fStartTime<fEndTime) addHist->GetAxis(0)->SetRange(fStartTime-1,fEndTime+1);
    if (fStartRun<fEndRun) addHist->GetAxis(3)->SetRange(fStartRun-1,fEndRun+1);
  }
  //
  //
  // 2. extraction of the information
  //
  if (fVdriftArray) 
  {
    fVdriftArray->Delete();
    delete fVdriftArray;
  }
  fVdriftArray = new TObjArray();
  AddAlignmentGraphs(fVdriftArray,fTimeDrift);
  AddHistoGraphs(fVdriftArray,fTimeDrift,fMinEntries);
  AddLaserGraphs(fVdriftArray,fTimeDrift);
  
  //
  // 3. Append QA plots
  //
  MakeDefaultPlots(fVdriftArray,fVdriftArray);

  //
  // 4. validate OCDB entries
  //
  if(fSwitchOnValidation==kTRUE && ValidateTimeDrift()==kFALSE) { 
    AliWarning("TPC time drift OCDB parameters out of range!");
    return 1;
  }
  //
  //4.b make alignment
  //
  MakeFitTime();
  TFile * ftime= TFile::Open("fitITSVertex.root");
  if (ftime){
    TObject * alignmentTime=ftime->Get("FitCorrectionTime");
    if (alignmentTime) fVdriftArray->AddLast(alignmentTime);
  }
  //
  // 5.) Add the RecoParam and ClusterParam - for compatibility checks -different sets of parameters can invalidate calibration 
  //
  AliTPCClusterParam *clParam =   AliTPCcalibDB::Instance()->GetClusterParam();
  TObjArray *recoParams = new TObjArray(4) ;
  for (Int_t i=0;i<4;i++) recoParams->AddAt(AliTPCcalibDB::Instance()->GetRecoParam(i),i);
  fVdriftArray->AddLast(clParam);
  fVdriftArray->AddLast(recoParams);
  return 0;
}

//_____________________________________________________________________________
void AliTPCPreprocessorOffline::CalibTimeVdrift(const Char_t* file, Int_t ustartRun, Int_t uendRun, AliCDBStorage* pocdbStorage){
  //
  // make calibration of the drift velocity
  // Input parameters:
  //      file                   - the location of input file
  //      ustartRun, uendRun     - run validity period 
  //      pocdbStorage           - path to hte OCDB storage
  //                             - if empty - local storage 'pwd' uesed
  if (pocdbStorage) fOCDBstorage=pocdbStorage;
  else {
    TString localStorage = "local://"+gSystem->GetFromPipe("pwd")+"/OCDB"; 
    fOCDBstorage=AliCDBManager::Instance()->GetStorage(localStorage.Data());
  }

  //
  // 1. Extract the calibration object form file, may have a number of layouts
  TFile fcalib(file);
  TObject* obj = dynamic_cast<TObject*>(fcalib.Get("TPCCalib"));
  TObjArray* array = dynamic_cast<TObjArray*>(obj);
  TDirectory* dir = dynamic_cast<TDirectory*>(obj);
  AliTPCcalibTime* timeDrift = NULL;
  if (dir) {
    timeDrift = dynamic_cast<AliTPCcalibTime*>(dir->Get("calibTime"));
  }
  else if (array){
    timeDrift = (AliTPCcalibTime *)array->FindObject("calibTime");
  } else {
    timeDrift = (AliTPCcalibTime*)fcalib.Get("calibTime");
  }
  if(!timeDrift) return;

  //calculate the calibration
  if (CalibTimeVdrift(timeDrift,ustartRun,uendRun)!=0) return;

  //
  //
  // 6. update of OCDB
  //
  //
  UpdateOCDBDrift(ustartRun,uendRun,fOCDBstorage);
}

//_____________________________________________________________________________
AliCDBEntry* AliTPCPreprocessorOffline::CreateDriftCDBentryObject(Int_t ustartRun, Int_t uendRun)
{
  //create the CDB entry (and cache internally)
  
  //reset if we dont have anything
  if (!fVdriftArray) 
  {
    delete fDriftCDBentry; fDriftCDBentry=NULL;
    return NULL;
  }

  //will be owned by the cdb entry once attached
  AliCDBMetaData* metaData = new AliCDBMetaData;
  metaData->SetObjectClassName("TObjArray");
  metaData->SetResponsible("Marian Ivanov");
  metaData->SetBeamPeriod(1);
  metaData->SetAliRootVersion("05-25-01"); //root version
  metaData->SetComment("Calibration of the time dependence of the drift velocity");
  
  AliCDBId id1("TPC/Calib/TimeDrift", ustartRun, uendRun);

  //now the entry owns the metadata, but NOT the data
  delete fDriftCDBentry;
  fDriftCDBentry=new AliCDBEntry(fVdriftArray,id1,metaData,kFALSE);
  
  return fDriftCDBentry;
}

//_____________________________________________________________________________
void AliTPCPreprocessorOffline::UpdateOCDBDrift( Int_t ustartRun, Int_t uendRun,  AliCDBStorage* storage ){
  //
  // Update OCDB 
  //
  Bool_t status=kFALSE;
  if (CreateDriftCDBentryObject(ustartRun,uendRun))
  {
    status=storage->Put(fDriftCDBentry); 
  }
  if (status==kFALSE) fCalibrationStatus|=kCalibFailedExport ;
}

void AliTPCPreprocessorOffline::TakeOwnershipDriftCDBEntry(){
	fVdriftArray = NULL;
	fDriftCDBentry = NULL;
}

//_____________________________________________________________________________
Bool_t AliTPCPreprocessorOffline::ValidateTimeGain()
{
  //
  // Validate time gain corrections 
  //
  AliInfo("ValidateTimeGain..." );
  Float_t minGain = fMinGain;
  Float_t maxGain = fMaxGain;

  TGraphErrors *gr = (TGraphErrors*)fGainArray->FindObject("TGRAPHERRORS_MEAN_GAIN_BEAM_ALL");

  // ===| treat the case if the combined calibration should be used |==========
  if (fGainCalibrationType==kCombinedGainCalib) {
    gr = (TGraphErrors*)fGainArrayCombined->FindObject("TGRAPHERRORS_MEAN_GAIN_BEAM_ALL");
  }

  if (!gr) {
    gr = (TGraphErrors*)fGainArray->FindObject("TGRAPHERRORS_MEAN_GAIN_COSMIC_ALL");
    if (fGainCalibrationType==kCombinedGainCalib) {
      gr = (TGraphErrors*)fGainArrayCombined->FindObject("TGRAPHERRORS_MEAN_GAIN_COSMIC_ALL");
    }
    if (!gr) 
    { 
      fCalibrationStatus |= kCalibFailedTimeGain;
      return kFALSE;
    }
    AliInfo("Assuming given run is a cosmic run. Using gain calibration from Fermi-plateau muons.");
  }
  if(gr->GetN()<1) 
  { 
    fCalibrationStatus |= kCalibFailedTimeGain;
    return kFALSE;
  }

  // check whether gain in the range
  for(Int_t iPoint=0; iPoint<gr->GetN(); iPoint++) 
  {
    if(gr->GetY()[iPoint] < minGain || gr->GetY()[iPoint] > maxGain)  
    { 
      fCalibrationStatus |= kCalibFailedTimeGain;
      AliError(TString::Format("Gain point outside of range %f != (%f,%f)",gr->GetY()[iPoint], minGain, maxGain ).Data());
      return kFALSE;
    }
  }

  AliInfo("Validation successful");
  return kTRUE;
}

//_____________________________________________________________________________
AliTPCPreprocessorOffline::EGainCalibType AliTPCPreprocessorOffline::GetGainCalibrationTypeFromString(const TString& type)
{
  //
  // return the gain calibration type analysing the string 'type'
  // if an error occurs, return kNGainCalibTypes
  //
  if (type.IsNull()) {
    ::Warning("AliTPCPreprocessorOffline::GetGainCalibrationTypeFromString","Empty gain calibration type string");
    return kNGainCalibTypes;
  }

  if (!type.IsDigit()) {
    ::Warning("AliTPCPreprocessorOffline::GetGainCalibrationTypeFromString","Gain calibration string '%s' is not a digit", type.Data());
    return kNGainCalibTypes;
  }

  const Int_t gainType = type.Atoi();
  if (gainType<0 || gainType>=Int_t(kNGainCalibTypes)) {
    ::Warning("AliTPCPreprocessorOffline::GetGainCalibrationTypeFromString","Gain calibration string '%s' is not a valid gain calibration type (0-%d)", type.Data(), kNGainCalibTypes-1);
    return kNGainCalibTypes;
  }

  return EGainCalibType(gainType);
}

//_____________________________________________________________________________
Bool_t AliTPCPreprocessorOffline::SetGainCalibrationType(const TString& type)
{
  //
  // set the gain calibration type analysing the string 'type'
  //
  EGainCalibType gainType=GetGainCalibrationTypeFromString(type);
  if (type==kNGainCalibTypes) return kFALSE;

  fGainCalibrationType=gainType;

  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliTPCPreprocessorOffline::ProduceCombinedGainCalibration()
{
  //
  // This function will produce the combined calibratin of the objects presently stored in the OCDB
  // and the calibration extracted in this calibration
  //

  // ===| write some info |=====================================================
  AliCDBManager *cdbMan = AliCDBManager::Instance();
  TString calibTimeGainStorage=cdbMan->GetDefaultStorage()->GetURI();
  const AliCDBEntry *e=cdbMan->Get("TPC/Calib/TimeGain");
  const TString timeGainID=e->GetId().ToString();
  if (cdbMan->GetSpecificStorage("TPC/Calib/TimeGain")) {
    calibTimeGainStorage=cdbMan->GetSpecificStorage("TPC/Calib/TimeGain")->GetURI();
  }
  AliInfoF("Using TimeGain '%s','%s' for combined gain calibration", calibTimeGainStorage.Data(), timeGainID.Data());

  // ===| get latest gain calibration from OCDB |===============================
  TObjArray *gainOCDB = AliTPCcalibDB::Instance()->GetTimeGainSplines();
  if (!gainOCDB) {
    AliError("Could not retrieve gain calibration from OCDB. Cannot perform combined calibration.");
    return kFALSE;
  }

  const Int_t nOCDB=gainOCDB->GetEntriesFast();
  const Int_t nThis=fGainArray->GetEntriesFast();

//   // ===| check consistency of entries |=======================================
//   TString type("this");
//   if ( nOCDB != nThis ) {
//     AliError(TString::Format("Entries in present OCDB calibration and this pass differ: %d != %d", nOCDB, nThis));
//     TObjArray *a1=gainOCDB->GetEntriesFast();
//     TObjArray *a2=fGainArray->GetEntriesFast();
//     if (nThis>nOCDB) {
//       a2=gainOCDB->GetEntriesFast();
//       a1=fGainArray->GetEntriesFast();
//       type="OCDB";
//     }
//
//     for (Int_t icalib=0; icalib<a1->GetEntriesFast(); ++icalib) {
//       const TObject *o=a1->At(icalib);
//       if (!o) continue;
//       if (!a2->FindObject(o->GetName())) {
//         AliError(TString::Format("Could not find '%s' in %s calibration", o->GetName(), type.Data()));
//       }
//     }
//     return kFALSE;
//   }

  // ===| create combined gain array if needed and reset |=====================
  if (!fGainArrayCombined) fGainArrayCombined=new TObjArray(nThis);
  //fGainArrayCombined->SetOwner(); // either not owner, or all steering objects must be cloned
  fGainArrayCombined->Clear();

  // ===| explicitly treat entries |===========================================
  TGraphErrors *grOCDB     = 0x0;
  TGraphErrors *grThis     = 0x0;
  TGraphErrors *grCombined = 0x0;
  AliSplineFit *splineFit  = 0x0;

  Bool_t error=kFALSE;

  // ---| gain vs time |--------------------------------------------------------
  GetGraphs("TGRAPHERRORS_MEAN_GAIN_BEAM_ALL", grOCDB, grThis);
  AliInfo(Form("Graphs: %p, %p: %s, %s", grOCDB, grThis, grOCDB?grOCDB->GetName():"", grThis?grThis->GetName():""));
  if (!grOCDB || !grThis) {
    AliError("Gain vs. time for beam cannot be processed");
    if (fGainMIP) {
      error=kTRUE;
    }
  } else {
    grCombined = CombineGraphs(grOCDB, grThis, 1 );
    AliInfo(Form("Combined: %p: %s", grCombined, grCombined?grCombined->GetName():""));
    if (grCombined) {
      splineFit = AliTPCcalibTimeGain::MakeSplineFit(grCombined);
      fGainArrayCombined->AddAt(splineFit ,0);
      fGainArrayCombined->AddAt(grCombined,2);
    } else {
      AliError("Gain vs. time for beam cannot be processed");
      error=kTRUE;
    }
  }

  GetGraphs("TGRAPHERRORS_MEAN_GAIN_COSMIC_ALL", grOCDB, grThis);
  AliInfo(Form("Graphs: %p, %p: %s, %s", grOCDB, grThis, grOCDB?grOCDB->GetName():"", grThis?grThis->GetName():""));
  if (!grOCDB || !grThis) {
    AliError("Gain vs. time from cosmics cannot be processed");
    if (fGainCosmic) {
      error=kTRUE;
    }
  } else {
    grCombined = CombineGraphs(grOCDB, grThis, 1);
    AliInfo(Form("Combined: %p: %s", grCombined, grCombined?grCombined->GetName():""));
    if (grCombined) {
      splineFit = AliTPCcalibTimeGain::MakeSplineFit(grCombined);
      fGainArrayCombined->AddAt(splineFit ,1);
      fGainArrayCombined->AddAt(grCombined,3);
    } else {
      AliError("Gain vs. time from cosmics cannot be processed");
      error=kTRUE;
    }

  }

  // ===| steering objects |====================================================
  TObjArray steeringObjectNames;
  steeringObjectNames.Add(new TNamed("GainSlopesHV","1"));
  steeringObjectNames.Add(new TNamed("GainSlopesPT","1"));
  steeringObjectNames.Add(new TNamed("AliTPCClusterParam","1"));
  steeringObjectNames.Add(new TNamed("TObjArray","1"));

  for (Int_t isteer=0; isteer<steeringObjectNames.GetEntriesFast(); ++isteer) {
    TObject *objName  = steeringObjectNames.At(isteer);
    TObject *steerObj = fGainArray->FindObject(objName->GetName());
    if (!steerObj) {
      AliErrorF("%s cannot be processed", objName->GetName());
      // ---| check if missing object should produce an error |---
      if (TString(objName->GetTitle()).Atoi()) {
        error=kTRUE;
      }
      continue;
    }
    fGainArrayCombined->AddLast(steerObj);
  }

  // ===| attachement is not applied, simply copy |=============================
  fGainArrayCombined->AddLast(fGainArray->FindObject("TGRAPHERRORS_MEAN_ATTACHMENT_BEAM_ALL"));

  // ===| multiplicative corrections |==========================================
  TObjArray graphsToCombine;
  graphsToCombine.Add(new TNamed("TGRAPHERRORS_MEANQTOT_PADREGIONGAIN_BEAM_ALL","1"));
  graphsToCombine.Add(new TNamed("TGRAPHERRORS_MEANQMAX_PADREGIONGAIN_BEAM_ALL","1"));
  graphsToCombine.Add(new TNamed("TGRAPHERRORS_MEANQMAX_MULTIPLICITYDEPENDENCE_BEAM_ALL","0"));
  graphsToCombine.Add(new TNamed("TGRAPHERRORS_MEANQTOT_MULTIPLICITYDEPENDENCE_BEAM_ALL","0"));
  graphsToCombine.Add(new TNamed("TGRAPHERRORS_MEAN_CHAMBERGAIN_SHORT_BEAM_ALL","1"));
  graphsToCombine.Add(new TNamed("TGRAPHERRORS_MEAN_CHAMBERGAIN_MEDIUM_BEAM_ALL","1"));
  graphsToCombine.Add(new TNamed("TGRAPHERRORS_MEAN_CHAMBERGAIN_LONG_BEAM_ALL","1"));
  graphsToCombine.Add(new TNamed("TGRAPHERRORS_QMAX_DIPANGLE_SHORT_BEAM_ALL","1"));
  graphsToCombine.Add(new TNamed("TGRAPHERRORS_QTOT_DIPANGLE_SHORT_BEAM_ALL","1"));
  graphsToCombine.Add(new TNamed("TGRAPHERRORS_QMAX_DIPANGLE_MEDIUM_BEAM_ALL","1"));
  graphsToCombine.Add(new TNamed("TGRAPHERRORS_QTOT_DIPANGLE_MEDIUM_BEAM_ALL","1"));
  graphsToCombine.Add(new TNamed("TGRAPHERRORS_QMAX_DIPANGLE_LONG_BEAM_ALL","1"));
  graphsToCombine.Add(new TNamed("TGRAPHERRORS_QTOT_DIPANGLE_LONG_BEAM_ALL","1"));
  graphsToCombine.Add(new TNamed("TGRAPHERRORS_QMAX_DIPANGLE_ABSOLUTE_BEAM_ALL","1"));
  graphsToCombine.Add(new TNamed("TGRAPHERRORS_QTOT_DIPANGLE_ABSOLUTE_BEAM_ALL","1"));

  for (Int_t igraph=0; igraph<graphsToCombine.GetEntriesFast(); ++igraph) {
    TObject *graphObj=graphsToCombine.At(igraph);
    GetGraphs(graphObj->GetName(), grOCDB, grThis);

    // --- check graphs
    if (!grOCDB || !grThis) {
      AliError(TString::Format("%s cannot be processed missing object(s) OCDB: %d, This: %d", graphObj->GetName(), grOCDB!=0x0, grThis!=0x0));
      if (TString(graphObj->GetTitle()).Atoi()) {
        error=kTRUE;
      }
      continue;
    }

    // --- combine graphs
    grCombined = CombineGraphs(grOCDB, grThis);
    fGainArrayCombined->AddLast(grCombined);

    // --- normalize pad region gain
    if ( TString(graphObj->GetName()).Contains("PADREGIONGAIN") ) {
      // first normalize to the weighted mean using the number of
      // rows in the different pad regions
      NormaliseYToWeightedMeandEdx(grCombined);
      // now correct for the effect the truncated mean in the dE/dx
      // calculation would have
      NormaliseYToTruncateddEdx(grCombined);
    }

    // --- for dip angle graphs normalize and do fit
    if ( TString(graphObj->GetName()).Contains("DIPANGLE") ) {
      NormaliseYToMean(grCombined);
      TF1 * fun= new TF1("","1++abs(x)++abs(x*x)");
      grCombined->Fit(fun,"w","rob=0.9",-0.8,0.8);
      TString funName(graphObj->GetName());
      funName.ReplaceAll("TGRAPHERRORS","TF1");
      fun->SetNameTitle(funName,funName);
      fGainArrayCombined->AddLast(fun);
    }
  }

  return !error;
}

//_____________________________________________________________________________
void AliTPCPreprocessorOffline::GetGraphs(const char* name, TGraphErrors* &grOCDB, TGraphErrors* &grThis)
{
  //
  // get graphs from OCDB gain array and local gain array
  //

  grOCDB=0x0;
  grThis=0x0;

  TObjArray *gainOCDB = AliTPCcalibDB::Instance()->GetTimeGainSplines();
  if (!gainOCDB) return;

  grOCDB = dynamic_cast<TGraphErrors*>(gainOCDB  ->FindObject(name));
  grThis = dynamic_cast<TGraphErrors*>(fGainArray->FindObject(name));

  if (!grOCDB) {
    AliError(TString::Format("Could not find graph '%s' in OCDB",name));
  }

  if (!grThis) {
    AliError(TString::Format("Could not find graph '%s' in present calibration",name));
  }
}

//_____________________________________________________________________________
TGraphErrors* AliTPCPreprocessorOffline::CombineGraphs(TGraphErrors *grOCDB, TGraphErrors *grThis, const Int_t type/*=0*/, const Bool_t multiply/*=kTRUE*/)
{
  //
  // Combine two graphs.
  // type
  //      0: Combine point by point, only do interpolation using Eval
  //      1: Combine using EvalConst
  //      2: Combine using Eval
  //
  // If mult is true, multiply the data points, otherwise add
  //

  Double_t x1,y1,x2,y2,x1Err,y1Err,y2Err,yNew,yNewErr;
  x1=y1=x2=y2=x1Err=y1Err=y2Err=yNew=yNewErr=0;

  const Int_t nOCDB=grOCDB->GetN();
  const Int_t nThis=grThis->GetN();

  // ===| output graph |========================================================
  //      copy from gThis to retain the attributes
  TGraphErrors *grCombined=new TGraphErrors(*grThis);
  grCombined->Set(0);

  // ===| sort graphs for better comparison |===================================
  grOCDB->Sort();
  grThis->Sort();

  const Double_t xMinOCDB = grOCDB->GetX()[0];
  const Double_t xMaxOCDB = grOCDB->GetX()[nOCDB-1];
  const Double_t xMinThis = grThis->GetX()[0];
  const Double_t xMaxThis = grThis->GetX()[nThis-1];

  // ===========================================================================
  // ===| treat combination types |=============================================
  //
  const Double_t kVerySmall=1e-10;
  Int_t ninterpol=0;

  const Bool_t interpol  = (type==1) || (type==2);
  const Bool_t evalConst = type!=2;

  Bool_t fallback=kFALSE;

  if (type==0) {
    if ( (xMinOCDB+kVerySmall < xMinThis) || (xMaxOCDB-kVerySmall > xMaxThis) ) {
      AliErrorF("Point by point combination of %s requested, but ranges are not compatible: [%.2f, %.2f] != [%.2f,%.2f]", grThis->GetName(), xMinOCDB, xMaxOCDB, xMinThis, xMaxThis);
      return 0x0;
    }
  }

  for (Int_t ipoint=0; ipoint<nThis; ++ipoint) {
    x2=y2=0;
    grThis->GetPoint(ipoint, x1, y1);
    x1Err=grThis->GetErrorX(ipoint);
    y1Err=grThis->GetErrorY(ipoint);

    if (type==0 && !fallback) {
      if (ipoint<nOCDB) {
        grOCDB->GetPoint(ipoint, x2, y2);
        y2Err=grOCDB->GetErrorY(ipoint);
        fallback=!TMath::AreEqualAbs(x1, x2, kVerySmall);
      }
      else {
        fallback=kTRUE;
      }
    }

    if (interpol || fallback) {
      GetPointWithError(grOCDB, x1, y2, y2Err, evalConst);
      ++ninterpol;
    }

    if (multiply) {
      yNew = y1*y2;
      yNewErr=0.;
      if (!TMath::AreEqualAbs(yNew, 0., kVerySmall)) {
        yNewErr = TMath::Sqrt( (y1Err*y1Err)/(y1*y1) + (y2Err*y2Err)/(y2*y2) ) * yNew;
      } else {
        AliErrorF("Cannot combined points from graph %s This y = %.4g +- %.4g, OCDB y = %.4g +- %.4g", grCombined->GetName(), y1, y1Err, y2, y2Err);
      }
    }
    else {
      yNew = y1+y2;
      yNewErr = TMath::Sqrt( (y1Err*y1Err) + (y2Err*y2Err) );
    }

    const Int_t point=grCombined->GetN();
    grCombined->SetPoint     (point, x1   , yNew   );
    grCombined->SetPointError(point, x1Err, yNewErr);
  }

  if ( (type==0) && ninterpol ) AliWarningF("%d number of points were interpolated for %s, although point-by-point combination was requested",ninterpol, grCombined->GetName());

  return grCombined;
}

//_____________________________________________________________________________
Bool_t AliTPCPreprocessorOffline::GetPointWithError(const TGraphErrors *gr, const Double_t xPos, Double_t &y, Double_t &ey, Bool_t evalConst/*=kTRUE*/)
{
  //
  // The function assumes the points to be sorted in x
  //
  // Evaluate the graph at xPos, do linear interpolation between neighbouring points
  // if 'evalConst' is true, the first/last point will be returned in case xPos is outside the graph range
  //
  // Errors in y are also calculated by linear interpolation between neighbouring points
  // Errors in x are not treated
  //
  // Return value indicates if xPos was inside the range

  // ===| reset input values |==================================================
  y=ey=0.;

  // ===| treat case of 0 point |===============================================
  if (!gr || gr->GetN()==0) return kFALSE;

  // ===| init variables |======================================================
  Double_t x1,x2,y1,y2, ey1,ey2;
  x1=x2=y1=y2=ey1=ey2=0.;

  const Int_t npoints=gr->GetN();
  const Double_t xmin=gr->GetX()[0];
  const Double_t xmax=gr->GetX()[npoints-1];

  const Bool_t belowBound = xPos<xmin;
  const Bool_t aboveBound = xPos>xmax;
  const Bool_t returnValue = !(belowBound || aboveBound);

  // ===| treat case of 1 point |===============================================
  if (npoints==1) {
    y  = gr->GetY()[0];
    ey = gr->GetErrorY(0);
    return returnValue;
  }

  // ===| treat eval const case |===============================================
  if (evalConst) {
    if (belowBound) {
      y  = gr->GetY()[0];
      ey = gr->GetErrorY(0);
      return returnValue;
    }

    if (aboveBound) {
      y  = gr->GetY()[npoints-1];
      ey = gr->GetErrorY(npoints-1);
      return returnValue;
    }
  }

  // ===| 2 and more points |===================================================
  Int_t point = TMath::BinarySearch(npoints, gr->GetX(), xPos);
  Printf("n, i: %d, %d", npoints, point);
  if (point==-1)        point=0;
  if (point==npoints-1) --point;

  gr->GetPoint(point, x1, y1);
  ey1 = gr->GetErrorY(point);

  gr->GetPoint(point+1, x2, y2);
  ey2 = gr->GetErrorY(point+1);

  Printf("%d, (%.2f, %.2f), (%.2f, %.2f)", point, x1, y1, x2, y2);

  if ( !(x2>x1) ) {
     AliTPCPreprocessorOffline p;
     p.Error("GetPointWithError","Graph not sorted, or error in extracting points");
     return kFALSE;
  }

  y = y1 + (y2-y1)  /(x2-x1) * (xPos-x1);
  ey=ey1 + (ey2-ey1)/(x2-x1) * (xPos-x1);

  // ===| linear error increase outside bounds |================================
  if (belowBound) {
    ey=ey1 + (ey1)/(x2-x1) * TMath::Abs(xPos-x1);
  }
  else if (aboveBound) {
    ey=ey2 + (ey1)/(x2-x1) * TMath::Abs(xPos-x2);
  }

  return returnValue;
}

//_____________________________________________________________________________
Bool_t AliTPCPreprocessorOffline::ValidateTimeDrift()
{
  //
  // Validate time drift velocity corrections 
  //
  AliInfo("ValidateTimeDrift..." );

  Float_t maxVDriftCorr = fMaxVdriftCorr;

  TGraphErrors* gr = (TGraphErrors*)fVdriftArray->FindObject("ALIGN_ITSB_TPC_DRIFTVD");
  AliInfo(Form("ALIGN_ITSB_TPC_DRIFTVD graph = %p",gr));
  if (!gr)
  {
    gr = (TGraphErrors*)fVdriftArray->FindObject("ALIGN_TOFB_TPC_DRIFTVD");
    AliInfo(Form("ALIGN_TOFB_TPC_DRIFTVD graph = %p",gr));
  }

  if(!gr) 
  {
    fCalibrationStatus|=kCalibFailedTimeDrift;
    return kFALSE;
  }
  
  // for now we validate even with low statistics
  ////check if we have enough statistics
  //if (fNtracksVdrift<fMinTracksVdrift) 
  //{
  //  fCalibrationStatus|=kCalibFailedTimeDrift;
  //  return kFALSE;
  //}

  if(gr->GetN()<1)  { 
    AliInfo(Form("ALIGN_ITSB_TPC_DRIFTVD number of points = %d",gr->GetN()));
    {
      fCalibrationStatus|=kCalibFailedTimeDrift;
      return kFALSE;
    }
  }

  // check whether drift velocity corrections in the range
  for(Int_t iPoint = 0; iPoint<gr->GetN(); iPoint++) 
  {
    //AliInfo(Form("Y value from the graph: %f",TMath::Abs(gr->GetY()[iPoint])));
    if(TMath::Abs(gr->GetY()[iPoint]) > maxVDriftCorr)  
    {
      fCalibrationStatus|=kCalibFailedTimeDrift;
      return kFALSE;
    }
  }

return kTRUE;
}

//_____________________________________________________________________________
void AliTPCPreprocessorOffline::UpdateDriftParam(AliTPCParam *param, TObjArray *const arr, Int_t lstartRun){
  //
  //  update the OCDB entry for the nominal time0
  //
  //
  //  AliTPCParam * param = AliTPCcalibDB::Instance()->GetParameters();
  AliTPCParam *paramNew = (AliTPCParam *)param->Clone();
  TGraphErrors *grT =  (TGraphErrors *)arr->FindObject("ALIGN_ITSM_TPC_T0");
  Double_t deltaTcm = TMath::Median(grT->GetN(),grT->GetY());
  Double_t deltaT   = deltaTcm/param->GetDriftV();
  paramNew->SetL1Delay(param->GetL1Delay()-deltaT);
  paramNew->Update();

  AliCDBMetaData *metaData= new AliCDBMetaData();
  metaData->SetObjectClassName("TObjArray");
  metaData->SetResponsible("Marian Ivanov");
  metaData->SetBeamPeriod(1);
  metaData->SetAliRootVersion("05-25-02"); //root version
  metaData->SetComment("Updated calibration of nominal time 0");
  AliCDBId* id1=NULL;
  id1=new AliCDBId("TPC/Calib/Parameters", lstartRun, AliCDBRunRange::Infinity());
  Bool_t status = fOCDBstorage->Put(param, (*id1), metaData);
  if (status==kFALSE) fCalibrationStatus|=kCalibFailedExport ;
}

//_____________________________________________________________________________
void AliTPCPreprocessorOffline::PrintArray(TObjArray *array){
  //
  // Print the names of the entries in array
  //
  Int_t entries = array->GetEntries();
  for (Int_t i=0; i<entries; i++){
    if (!array->At(i)) continue;
    Printf("%d\t %s", i,  array->At(i)->GetName());
  }
}

//_____________________________________________________________________________
TGraphErrors* AliTPCPreprocessorOffline::FilterGraphDrift(TGraphErrors * graph, Float_t errSigmaCut, Float_t medianCutAbs){
  // 2 filters:
  //    1. filter graph - error cut errSigmaCut
  //    2. filter graph - medianCutAbs around median
  //
  // errSigmaCut   - cut on error
  // medianCutAbs  - cut on value around median
  Double_t dummy=0;               //   
  //
  // 1. filter graph - error cut errSigmaCut
  //              
  TGraphErrors *graphF; 
  graphF = AliTPCcalibDButil::FilterGraphMedianErr(graph,errSigmaCut,dummy);
  delete graph;
  if (!graphF) return 0;
  graph = AliTPCcalibDButil::FilterGraphMedianErr(graphF,errSigmaCut,dummy);
  delete graphF;
  if (!graph) return 0;
  //
  // filter graph - kMedianCutAbs around median
  // 
  graphF=FilterGraphMedianAbs(graph, medianCutAbs,dummy);
  delete graph;
  if (!graphF) return 0;
  graph=FilterGraphMedianAbs(graphF, medianCutAbs,dummy);
  delete graphF;
  if (!graph) return 0;
  return graph;
}

//_____________________________________________________________________________
TGraphErrors* AliTPCPreprocessorOffline::FilterGraphMedianAbs(TGraphErrors * graph, Float_t cut,Double_t &medianY){
  //
  // filter outlyer measurement
  // Only points around median +- cut filtered 
  //
  if (!graph) return  0;
  Int_t kMinPoints=2;
  Int_t npoints0 = graph->GetN();
  Int_t npoints=0;
  Float_t  rmsY=0;
  Double_t *outx=new Double_t[npoints0];
  Double_t *outy=new Double_t[npoints0];
  Double_t *errx=new Double_t[npoints0];
  Double_t *erry=new Double_t[npoints0];
  //
  //
  if (npoints0<kMinPoints) {
    delete []outx;
    delete []outy;
    delete []errx;
    delete []erry;
    return 0;
  }
  for (Int_t iter=0; iter<3; iter++){
    npoints=0;
    for (Int_t ipoint=0; ipoint<npoints0; ipoint++){
      if (graph->GetY()[ipoint]==0) continue;
      if (iter>0 &&TMath::Abs(graph->GetY()[ipoint]-medianY)>cut) continue;  
      outx[npoints]  = graph->GetX()[ipoint];
      outy[npoints]  = graph->GetY()[ipoint];
      errx[npoints]  = graph->GetErrorX(ipoint);
      erry[npoints]  = graph->GetErrorY(ipoint);
      npoints++;
    }
    if (npoints<=1) break;
    medianY  =TMath::Median(npoints,outy);
    rmsY   =TMath::RMS(npoints,outy);
  }
  TGraphErrors *graphOut=0;
  if (npoints>1) graphOut= new TGraphErrors(npoints,outx,outy,errx,erry); 
  delete []outx;
  delete []outy;
  delete []errx;
  delete []erry;
  return graphOut;
}

//_____________________________________________________________________________
void AliTPCPreprocessorOffline::AddHistoGraphs(  TObjArray * vdriftArray, AliTPCcalibTime * const timeDrift, Int_t minEntries){
  //
  // Add graphs corresponding to the alignment
  //
  const Double_t kErrSigmaCut=5;      // error sigma cut - for filtering
  const Double_t kMedianCutAbs=0.03;  // error sigma cut - for filtering
  //
  TObjArray * array=timeDrift->GetHistoDrift();
  if (array){
    THnSparse* hist=NULL;
    // 2.a) cosmics with different triggers
    for (Int_t i=0; i<array->GetEntriesFast();i++){
      hist=(THnSparseF*)array->UncheckedAt(i);
      if(!hist) continue;
      if (hist->GetEntries()<minEntries) continue;
      //hist->Print();
      TString name=hist->GetName();
      Int_t dim[4]={0,1,2,3};
      THnSparse* newHist=hist->Projection(4,dim);
      newHist->SetName(name);
      TGraphErrors* graph=AliTPCcalibBase::FitSlices(newHist,2,0,400,100,0.05,0.95, kTRUE);
      delete newHist;
      if (!graph) {
	AliInfo(Form("Graph =%s filtered out", name.Data()));
	continue;
      }
      AliInfo(Form("name=%s graph=%i, N=%i", name.Data(), graph==0, graph->GetN()));
      Int_t pos=name.Index("_");
      name=name(pos,name.Capacity()-pos);
      TString graphName=graph->ClassName();
      graphName+=name;
      graphName.ToUpper();
      //
      graph = FilterGraphDrift(graph, kErrSigmaCut, kMedianCutAbs);
      //
      if (graph){
        graph->SetMarkerStyle(i%8+20);
        graph->SetMarkerColor(i%7);
        graph->GetXaxis()->SetTitle("Time");
        graph->GetYaxis()->SetTitle("v_{dcor}");
        graph->SetName(graphName);
        graph->SetTitle(graphName);
        AliInfo(Form("Graph %d\t=\t%s", i, graphName.Data()));
        vdriftArray->Add(graph);
      }
    }
  }
}

//_____________________________________________________________________________
void AliTPCPreprocessorOffline::AddAlignmentGraphs(  TObjArray * vdriftArray, AliTPCcalibTime *const timeDrift){
  //
  // Add graphs corresponding to alignment to the object array
  //
  TObjArray *arrayITS=0;
  TObjArray *arrayTOF=0;
  TObjArray *arrayTRD=0;
  TMatrixD *mstatITS=0;
  TMatrixD *mstatTOF=0;
  TMatrixD *mstatTRD=0;
  //
  arrayITS=timeDrift->GetAlignITSTPC();
  arrayTRD=timeDrift->GetAlignTRDTPC();
  arrayTOF=timeDrift->GetAlignTOFTPC();

  if (arrayITS->GetEntries()>0) mstatITS= AliTPCcalibDButil::MakeStatRelKalman(arrayITS,0.7,50,fMaxVdriftCorr);
  if (arrayTOF->GetEntries()>0) mstatTOF= AliTPCcalibDButil::MakeStatRelKalman(arrayTOF,0.7,1000,fMaxVdriftCorr);
  if (arrayTRD->GetEntries()>0) mstatTRD= AliTPCcalibDButil::MakeStatRelKalman(arrayTRD,0.7,50,fMaxVdriftCorr);
  //
  TObjArray * arrayITSP= AliTPCcalibDButil::SmoothRelKalman(arrayITS,mstatITS, 0, 5.);
  TObjArray * arrayITSM= AliTPCcalibDButil::SmoothRelKalman(arrayITS,mstatITS, 1, 5.);
  TObjArray * arrayITSB= AliTPCcalibDButil::SmoothRelKalman(arrayITSP,arrayITSM);
  TObjArray * arrayTOFP= AliTPCcalibDButil::SmoothRelKalman(arrayTOF,mstatTOF, 0, 5.);
  TObjArray * arrayTOFM= AliTPCcalibDButil::SmoothRelKalman(arrayTOF,mstatTOF, 1, 5.);
  TObjArray * arrayTOFB= AliTPCcalibDButil::SmoothRelKalman(arrayTOFP,arrayTOFM);

  TObjArray * arrayTRDP= 0x0;
  TObjArray * arrayTRDM= 0x0;
  TObjArray * arrayTRDB= 0x0;
  arrayTRDP= AliTPCcalibDButil::SmoothRelKalman(arrayTRD,mstatTRD, 0, 5.);
  arrayTRDM= AliTPCcalibDButil::SmoothRelKalman(arrayTRD,mstatTRD, 1, 5.);
  arrayTRDB= AliTPCcalibDButil::SmoothRelKalman(arrayTRDP,arrayTRDM);
  //
  //
  Int_t entries=TMath::Max(arrayITS->GetEntriesFast(),arrayTOF->GetEntriesFast());
  TObjArray *arrays[12]={arrayITS, arrayITSP, arrayITSM, arrayITSB,
			 arrayTRD, arrayTRDP, arrayTRDM, arrayTRDB,
			 arrayTOF, arrayTOFP, arrayTOFM, arrayTOFB};
  TString   grnames[12]={"ALIGN_ITS", "ALIGN_ITSP", "ALIGN_ITSM", "ALIGN_ITSB",
			 "ALIGN_TRD", "ALIGN_TRDP", "ALIGN_TRDM","ALIGN_TRDB",
			 "ALIGN_TOF", "ALIGN_TOFP", "ALIGN_TOFM","ALIGN_TOFB"};
  TString   grpar[9]={"DELTAPSI", "DELTATHETA", "DELTAPHI",
		      "DELTAX", "DELTAY", "DELTAZ",
		      "DRIFTVD", "T0", "VDGY"};

  
  TVectorD vX(entries);
  TVectorD vY(entries);
  TVectorD vEx(entries);
  TVectorD vEy(entries);
  TObjArray *arr=0;
  for (Int_t iarray=0; iarray<12; iarray++){
    arr = arrays[iarray];
    if (arr==0) continue;
    for (Int_t ipar=0; ipar<9; ipar++){      
      Int_t counter=0;
      for (Int_t itime=0; itime<arr->GetEntriesFast(); itime++){
	AliRelAlignerKalman * kalman = (AliRelAlignerKalman *) arr->UncheckedAt(itime);
	if (!kalman) continue;
	vX[counter]=kalman->GetTimeStamp();
	vY[counter]=(*(kalman->GetState()))[ipar];
	if (ipar==6) vY[counter]=1./(*(kalman->GetState()))[ipar]-1;
	vEx[counter]=0;
	vEy[counter]=TMath::Sqrt((*(kalman->GetStateCov()))(ipar,ipar));
	counter++;
      }
    
      TGraphErrors * graph=new TGraphErrors(counter, vX.GetMatrixArray(),
					  vY.GetMatrixArray(),
					  vEx.GetMatrixArray(),
					  vEy.GetMatrixArray());
      TString grName=grnames[iarray];
      grName+="_TPC_";
      grName+=grpar[ipar];
      graph->SetName(grName.Data());
      vdriftArray->AddLast(graph);
    }
    delete arrays[iarray];
  }  
}

//_____________________________________________________________________________
void AliTPCPreprocessorOffline::AddLaserGraphs(  TObjArray * vdriftArray, AliTPCcalibTime *timeDrift){
  //
  // add graphs for laser
  //
  const Double_t delayL0L1 = 0.071;  //this is hack for 1/2 weeks
  //THnSparse *hisN=0;
  TGraphErrors *grLaser[6]={0,0,0,0,0,0};
  //hisN = timeDrift->GetHistVdriftLaserA(0);
  if (timeDrift->GetHistVdriftLaserA(0)){
    grLaser[0]=MakeGraphFilter0(timeDrift->GetHistVdriftLaserA(0),0,2,5,delayL0L1);
    grLaser[0]->SetName("GRAPH_MEAN_DELAY_LASER_ALL_A");
    vdriftArray->AddLast(grLaser[0]);
  }    
  if (timeDrift->GetHistVdriftLaserA(1)){
    grLaser[1]=MakeGraphFilter0(timeDrift->GetHistVdriftLaserA(1),0,2,5);
    grLaser[1]->SetName("GRAPH_MEAN_DRIFT_LASER_ALL_A");
    vdriftArray->AddLast(grLaser[1]);
  }    
  if (timeDrift->GetHistVdriftLaserA(2)){
    grLaser[2]=MakeGraphFilter0(timeDrift->GetHistVdriftLaserA(2),0,2,5);
    grLaser[2]->SetName("GRAPH_MEAN_GLOBALYGRADIENT_LASER_ALL_A");
    vdriftArray->AddLast(grLaser[2]);
  }    
  if (timeDrift->GetHistVdriftLaserC(0)){
    grLaser[3]=MakeGraphFilter0(timeDrift->GetHistVdriftLaserC(0),0,2,5,delayL0L1);
    grLaser[3]->SetName("GRAPH_MEAN_DELAY_LASER_ALL_C");
    vdriftArray->AddLast(grLaser[3]);
  }    
  if (timeDrift->GetHistVdriftLaserC(1)){
    grLaser[4]=MakeGraphFilter0(timeDrift->GetHistVdriftLaserC(1),0,2,5);
    grLaser[4]->SetName("GRAPH_MEAN_DRIFT_LASER_ALL_C");
    vdriftArray->AddLast(grLaser[4]);
  }    
  if (timeDrift->GetHistVdriftLaserC(2)){
    grLaser[5]=MakeGraphFilter0(timeDrift->GetHistVdriftLaserC(2),0,2,5);
    grLaser[5]->SetName("GRAPH_MEAN_GLOBALYGRADIENT_LASER_ALL_C");    
    vdriftArray->AddLast(grLaser[5]);
  }    
  for (Int_t i=0; i<6;i++){
    if (grLaser[i]) {
      SetDefaultGraphDrift(grLaser[i], 1,(i+20));
      grLaser[i]->GetYaxis()->SetTitle("Laser Correction");
    }
  }
}
 
//_____________________________________________________________________________
TGraphErrors * AliTPCPreprocessorOffline::MakeGraphFilter0(THnSparse *hisN, Int_t itime, Int_t ival, Int_t minEntries, Double_t offset){
  //
  // Make graph with mean values and rms
  //
  hisN->GetAxis(itime)->SetRange(0,100000000);
  hisN->GetAxis(ival)->SetRange(0,100000000);
  TH1 * hisT      = hisN->Projection(itime);
  TH1 * hisV      = hisN->Projection(ival);
  //
  Int_t firstBinA = hisT->FindFirstBinAbove(2);
  Int_t lastBinA  = hisT->FindLastBinAbove(2);    
  Int_t firstBinV = hisV->FindFirstBinAbove(0);
  Int_t lastBinV  = hisV->FindLastBinAbove(0);    
  hisN->GetAxis(itime)->SetRange(firstBinA,lastBinA);
  hisN->GetAxis(ival)->SetRange(firstBinV,lastBinV);
  Int_t entries=0;
  for (Int_t ibin=firstBinA; ibin<=lastBinA; ibin++){
    Double_t cont = hisT->GetBinContent(ibin);
    if (cont<minEntries) continue;
    entries++;
  }
  TVectorD vecTime(entries);
  TVectorD vecMean0(entries);
  TVectorD vecRMS0(entries);
  TVectorD vecMean1(entries);
  TVectorD vecRMS1(entries);
  entries=0;
  for (Int_t ibin=firstBinA; ibin<=lastBinA; ibin++){
      Double_t cont = hisT->GetBinContent(ibin);
      if (cont<minEntries) continue;
      //hisN->GetAxis(itime)->SetRange(ibin-1,ibin+1);
      Int_t minBin = ibin-1;
      Int_t maxBin = ibin+1;
      if(minBin <= 0) minBin = 1;
      if(maxBin >= hisN->GetAxis(itime)->GetNbins()) maxBin = hisN->GetAxis(itime)->GetNbins()-1;
      hisN->GetAxis(itime)->SetRange(minBin,maxBin);

      Double_t time = hisT->GetBinCenter(ibin);
      TH1 * his = hisN->Projection(ival);
      Double_t nentries0= his->GetBinContent(his->FindBin(0));
      if (cont-nentries0<minEntries) continue;
      //
      his->SetBinContent(his->FindBin(0),0);
      vecTime[entries]=time;
      vecMean0[entries]=his->GetMean()+offset;
      vecMean1[entries]=his->GetMeanError();
      vecRMS0[entries] =his->GetRMS();
      vecRMS1[entries] =his->GetRMSError();
      delete his;  
      entries++;
  }
  delete hisT;
  delete hisV;
  TGraphErrors * graph =  new TGraphErrors(entries,vecTime.GetMatrixArray(), vecMean0.GetMatrixArray(),					   0, vecMean1.GetMatrixArray());

  return graph;
}

//_____________________________________________________________________________
void AliTPCPreprocessorOffline::SetDefaultGraphDrift(TGraph *graph, Int_t color, Int_t style){
  //
  // Set default style for QA views
  //
  graph->GetXaxis()->SetTimeDisplay(kTRUE);
  graph->GetXaxis()->SetTimeFormat("#splitline{%d/%m}{%H:%M}");
  graph->SetMaximum( 0.025);
  graph->SetMinimum(-0.025);
  graph->GetXaxis()->SetTitle("Time");
  graph->GetYaxis()->SetTitle("v_{dcorr}");
  //
  graph->GetYaxis()->SetLabelSize(0.03);
  graph->GetXaxis()->SetLabelSize(0.03);
  //
  graph->GetXaxis()->SetNdivisions(10,5,0);
  graph->GetYaxis()->SetNdivisions(10,5,0);
  //
  graph->GetXaxis()->SetLabelOffset(0.02);
  graph->GetYaxis()->SetLabelOffset(0.005);
  //
  graph->GetXaxis()->SetTitleOffset(1.3);
  graph->GetYaxis()->SetTitleOffset(1.2);
  //
  graph->SetMarkerColor(color);
  graph->SetLineColor(color);
  graph->SetMarkerStyle(style);
}

//_____________________________________________________________________________
void AliTPCPreprocessorOffline::SetPadStyle(TPad *pad, Float_t mx0, Float_t mx1, Float_t my0, Float_t my1){
  // 
  // Set default pad style for QA
  // 
  pad->SetTicks(1,1);
  pad->SetMargin(mx0,mx1,my0,my1);
}

//_____________________________________________________________________________
void AliTPCPreprocessorOffline::MakeDefaultPlots(TObjArray * const arr, TObjArray * /*picArray*/){
  //
  // 0. make a default QA plots
  // 1. Store them in the array
  //
  //
  Float_t mx0=0.12, mx1=0.1, my0=0.15, my1=0.1;
  //
  TGraphErrors* laserA       =(TGraphErrors*)arr->FindObject("GRAPH_MEAN_DRIFT_LASER_ALL_A");
  TGraphErrors* laserC       =(TGraphErrors*)arr->FindObject("GRAPH_MEAN_DRIFT_LASER_ALL_C");
  TGraphErrors* cosmic       =(TGraphErrors*)arr->FindObject("TGRAPHERRORS_MEAN_VDRIFT_COSMICS_ALL");
  TGraphErrors* cross        =(TGraphErrors*)arr->FindObject("TGRAPHERRORS_VDRIFT_CROSS_ALL");
  TGraphErrors* itstpcP       =(TGraphErrors*)arr->FindObject("ALIGN_ITSP_TPC_DRIFTVD");
  TGraphErrors* itstpcM       =(TGraphErrors*)arr->FindObject("ALIGN_ITSM_TPC_DRIFTVD");
  TGraphErrors* itstpcB       =(TGraphErrors*)arr->FindObject("ALIGN_ITSB_TPC_DRIFTVD");
  //
  if (laserA)  SetDefaultGraphDrift(laserA,2,25);
  if (laserC)  SetDefaultGraphDrift(laserC,4,26);
  if (cosmic)  SetDefaultGraphDrift(cosmic,3,27);
  if (cross)   SetDefaultGraphDrift(cross,4,28);
  if (itstpcP) SetDefaultGraphDrift(itstpcP,2,29);
  if (itstpcM) SetDefaultGraphDrift(itstpcM,4,30);
  if (itstpcB) SetDefaultGraphDrift(itstpcB,1,31);
  //
  //
  TPad *pad=0;
  //
  // Laser-Laser
  //
  if (laserA&&laserC){
    pad = new TCanvas("TPCLaserVDrift","TPCLaserVDrift");
    laserA->Draw("alp");
    SetPadStyle(pad,mx0,mx1,my0,my1);
    laserA->Draw("apl");
    laserC->Draw("p");
    TLegend *legend = new TLegend(mx0+0.01,1-my1-0.2, 0.5, 1-my1-0.01, "Drift velocity correction");
    legend->AddEntry(laserA,"Laser A side");
    legend->AddEntry(laserC,"Laser C side");
    legend->Draw();    
    //picArray->AddLast(pad);
  }

  if (itstpcP&&itstpcM&&itstpcB){
    pad = new TCanvas("ITSTPC","ITSTPC");
    itstpcP->Draw("alp");
    SetPadStyle(pad,mx0,mx1,my0,my1);    
    itstpcP->Draw("alp");
    gPad->Clear();
    itstpcM->Draw("apl");
    itstpcP->Draw("p");
    itstpcB->Draw("p");
    TLegend *legend = new TLegend(mx0+0.01,1-my1-0.2, 0.5, 1-my1-0.01, "Drift velocity correction");
    legend->AddEntry(itstpcP,"ITS-TPC smooth plus");
    legend->AddEntry(itstpcM,"ITS-TPC smooth minus");
    legend->AddEntry(itstpcB,"ITS-TPC smooth ");
    legend->Draw();    
    //picArray->AddLast(pad);
  }

  if (itstpcB&&laserA&&itstpcP&&itstpcM){
    pad = new TCanvas("ITSTPC_LASER","ITSTPC_LASER");
    SetPadStyle(pad,mx0,mx1,my0,my1);    
    laserA->Draw("alp");
    itstpcP->Draw("p");
    itstpcM->Draw("p");
    itstpcB->Draw("p");
    TLegend *legend = new TLegend(mx0+0.01,1-my1-0.2, 0.5, 1-my1-0.01, "Drift velocity correction");
    legend->AddEntry(laserA,"TPC laser");
    legend->AddEntry(itstpcP,"ITS-TPC smooth plus");   
    legend->AddEntry(itstpcM,"ITS-TPC smooth minus");   
    legend->AddEntry(itstpcB,"ITS-TPC smooth ");
    legend->Draw();
    //picArray->AddLast(pad);
  }

  if (itstpcP&&cross){ 
    pad = new TCanvas("ITSTPC_CROSS","ITSTPC_CROSS");
    SetPadStyle(pad,mx0,mx1,my0,my1);    
    itstpcP->Draw("alp");
    pad->Clear();
    cross->Draw("ap");
    itstpcP->Draw("p");
    //
    TLegend *legend = new TLegend(mx0+0.01,1-my1-0.2, 0.5, 1-my1-0.01, "Drift velocity correction");

    legend->AddEntry(cross,"TPC cross tracks");
    legend->AddEntry(itstpcB,"ITS-TPC smooth");
    legend->Draw();        
    //picArray->AddLast(pad);
  }
  if (itstpcP&&cosmic){ 
    pad = new TCanvas("ITSTPC_COSMIC","ITSTPC_COSMIC");
    SetPadStyle(pad,mx0,mx1,my0,my1);    
    itstpcP->Draw("alp");
    pad->Clear();
    cosmic->Draw("ap");
    itstpcP->Draw("p");
    //
    TLegend *legend = new TLegend(mx0+0.01,1-my1-0.2, 0.5, 1-my1-0.01, "Drift velocity correction");

    legend->AddEntry(cosmic,"TPC cross tracks0 up-down");
    legend->AddEntry(itstpcB,"ITS-TPC smooth");
    legend->Draw();        
    //picArray->AddLast(pad);
  }
}

//_____________________________________________________________________________
void AliTPCPreprocessorOffline::CalibTimeGain(const Char_t* fileName, Int_t startRunNumber, Int_t endRunNumber,  AliCDBStorage* fullStorage, AliCDBStorage* residualStorage){
  //
  // Update OCDB gain
  // fullStorage is where the full calibration object should go
  // residualStorage is where the residual calibration for QA purposes should go
  //
  if (fullStorage==0) {
    const TString localStorage = "local://"+gSystem->GetFromPipe("pwd")+"/OCDB";
    fullStorage = AliCDBManager::Instance()->GetStorage(localStorage.Data());
  }
  if ( (fGainCalibrationType==kResidualGainQA || fGainCalibrationType==kCombinedGainCalib) && residualStorage==0x0) {
    const TString localStorage = "local://"+gSystem->GetFromPipe("pwd")+"/residualOCDB";
    residualStorage = AliCDBManager::Instance()->GetStorage(localStorage.Data());
  }

  // dump info
  AliInfoF("Analyse gain from file %s", fileName);

  //
  // 1. Read gain values
  //
  ReadGainGlobal(fileName);

  //
  // 2. Extract calibration values
  //
  AnalyzeGain(startRunNumber,endRunNumber, 1000,1.43);
  AnalyzeAttachment(startRunNumber,endRunNumber);
  AnalyzePadRegionGain();
  AnalyzeGainMultiplicity();
  AnalyzeGainChamberByChamber();
  //
  AnalyzeGainDipAngle(0); // short pads
  AnalyzeGainDipAngle(1); // medium pads
  AnalyzeGainDipAngle(2); // long pads
  AnalyzeGainDipAngle(3); // absolute calibration on full track

  //
  // 2.a produce combined calibration if requested
  //

  Bool_t combinedGainSuccessful=kTRUE;
  if (fGainCalibrationType==kCombinedGainCalib) {
    combinedGainSuccessful=ProduceCombinedGainCalibration();
    AliInfoF("Result of combined gain calibration: %d", combinedGainSuccessful);
  }

  //
  // 3. Make control plots
  //
  MakeQAPlot(1.43);  

  //
  // 4. validate OCDB entries
  //
  if(fSwitchOnValidation==kTRUE &&
     (fGainCalibrationType==kFullGainCalib || fGainCalibrationType==kCombinedGainCalib) &&
     (!combinedGainSuccessful || !ValidateTimeGain()) ) {
    AliWarning("TPC time gain OCDB parameters out of range!");
    return;
  }

  //
  // 5. Update OCDB
  //
  if (fGainCalibrationType != kNoGainCalib ) {
    UpdateOCDBGain( startRunNumber, endRunNumber, fullStorage, residualStorage);
  }
}

//_____________________________________________________________________________
void AliTPCPreprocessorOffline::ReadGainGlobal(const Char_t* fileName){
  //
  // read calibration entries from file
  // 
  TFile *fcalib=TFile::Open(fileName);
  gROOT->cd();

  TObject* obj = dynamic_cast<TObject*>(fcalib->Get("TPCCalib"));
  TObjArray * array = dynamic_cast<TObjArray*>(obj);
  TDirectory * dir = dynamic_cast<TDirectory*>(obj);
  if (dir) {
    fGainMIP    = dynamic_cast<AliTPCcalibTimeGain *>(dir->Get("calibTimeGain"));
    fGainCosmic = dynamic_cast<AliTPCcalibTimeGain *>(dir->Get("calibTimeGainCosmic"));
    fGainMult   = dynamic_cast<AliTPCcalibGainMult *>(dir->Get("calibGainMult"));
  }
  else if (array){
    fGainMIP    = ( AliTPCcalibTimeGain *)array->FindObject("calibTimeGain");
    fGainCosmic = ( AliTPCcalibTimeGain *)array->FindObject("calibTimeGainCosmic");
    fGainMult   = ( AliTPCcalibGainMult *)array->FindObject("calibGainMult");
  }else{
    fGainMIP    = ( AliTPCcalibTimeGain *)fcalib->Get("calibTimeGain");
    fGainCosmic = ( AliTPCcalibTimeGain *)fcalib->Get("calibTimeGainCosmic");
    fGainMult   = ( AliTPCcalibGainMult *)fcalib->Get("calibGainMult");
  }
  if (!fGainMult){
    TFile calibMultFile("TPCMultObjects.root");
    fGainMult   = ( AliTPCcalibGainMult *)calibMultFile.Get("calibGainMult");
  }
  TH1 * hisT=0;
  Int_t firstBinA =0, lastBinA=0;

  if (fGainCosmic){ 
    hisT= fGainCosmic->GetHistGainTime()->Projection(1);
    firstBinA = hisT->FindFirstBinAbove(2);
    lastBinA  = hisT->FindLastBinAbove(2);    
    fGainCosmic->GetHistGainTime()->GetAxis(1)->SetRange(firstBinA,lastBinA);
    delete hisT;
  }

  if (fGainMIP){ 
    hisT= fGainMIP->GetHistGainTime()->Projection(1);
    firstBinA = hisT->FindFirstBinAbove(2);
    lastBinA  = hisT->FindLastBinAbove(2);    
    fGainMIP->GetHistGainTime()->GetAxis(1)->SetRange(firstBinA,lastBinA);
    delete hisT;
  }

  delete fcalib;
}

//_____________________________________________________________________________
Bool_t AliTPCPreprocessorOffline::AnalyzeGain(Int_t startRunNumber, Int_t endRunNumber, Int_t minEntriesGaussFit,  Float_t FPtoMIPratio){
  //
  // Analyze gain - produce the calibration graphs
  // Uses the MIP pions integrated over all chambers, eta and multiplicity
  // This provides the absolute normalization to channel 50 in the end
  //
  // TODO: o What happens in case the MIPs are not flat in eta, and/or multiplicity?
  //         In CPass0 it will not be flat.
  //         Can we fill the histograms flat in multiplicity?
  //       o Does the phi distribution play a role (inactive sector parts)?
  //       o In principle one could normalise over eta, by first projecting in slices
  //         of eta, normalize each slice and then sum
  //

  // 1.) try to create MIP spline
  if (fGainMIP) 
  {
    fGainMIP->GetHistGainTime()->GetAxis(5)->SetRangeUser(startRunNumber, endRunNumber);
    fGainMIP->GetHistGainTime()->GetAxis(2)->SetRangeUser(1.51,2.49); // only beam data
    fGainMIP->GetHistGainTime()->GetAxis(4)->SetRangeUser(0.39,0.51); // only MIP pions

    // ---| QA histogram |---
    if (fArrQAhist) {
      TH2D *hQA = fGainMIP->GetHistGainTime()->Projection(0,1);
      hQA->SetName("TGRAPHERRORS_MEAN_GAIN_BEAM_ALL_QA");
      hQA->SetTitle("MIP calibration collisions; time;d#it{E}/d#it{x} (arb. unit)");
      hQA->GetXaxis()->SetRange(hQA->FindFirstBinAbove(1), hQA->FindLastBinAbove(1));
      fArrQAhist->Add(hQA);
    }

    //
    fGraphMIP = AliTPCcalibBase::FitSlices(fGainMIP->GetHistGainTime(),0,1,minEntriesGaussFit,10,0.1,0.7);
    if (fGraphMIP->GetN()==0) fGraphMIP = 0x0;
    if (fGraphMIP) fFitMIP = AliTPCcalibTimeGain::MakeSplineFit(fGraphMIP);
    if (fGraphMIP) fGraphMIP->SetName("TGRAPHERRORS_MEAN_GAIN_BEAM_ALL");// set proper names according to naming convention
    fGainArray->AddAt(fFitMIP,0);
  } 

  // 2.) try to create Cosmic spline
  if (fGainCosmic)
  {
    fGainCosmic->GetHistGainTime()->GetAxis(2)->SetRangeUser(0.51,1.49); // only cosmics
    fGainCosmic->GetHistGainTime()->GetAxis(4)->SetRangeUser(20,100);    // only Fermi-Plateau muons

    // ---| QA histogram |---
    if (fArrQAhist) {
      TH2D *hQA = fGainMIP->GetHistGainTime()->Projection(0,1);
      hQA->SetName("TGRAPHERRORS_MEAN_GAIN_COSMIC_ALL_QA");
      hQA->SetTitle("MIP calibration cosmics; time;d#it{E}/d#it{x} (arb. unit)");
      hQA->GetXaxis()->SetRange(hQA->FindFirstBinAbove(1), hQA->FindLastBinAbove(1));
      fArrQAhist->Add(hQA);
    }

    //
    fGraphCosmic = AliTPCcalibBase::FitSlices(fGainCosmic->GetHistGainTime(),0,1,minEntriesGaussFit,10);
    if (fGraphCosmic->GetN()==0) fGraphCosmic = 0x0;
    //
    if (fGraphCosmic) {
      for(Int_t i=0; i < fGraphCosmic->GetN(); i++) {
	fGraphCosmic->GetY()[i] /= FPtoMIPratio;
	fGraphCosmic->GetEY()[i] /= FPtoMIPratio;
      }
    }
    //
    if (fGraphCosmic) fFitCosmic = AliTPCcalibTimeGain::MakeSplineFit(fGraphCosmic);
    if (fGraphCosmic) fGraphCosmic->SetName("TGRAPHERRORS_MEAN_GAIN_COSMIC_ALL"); // set proper names according to naming convention
    fGainArray->AddAt(fFitCosmic,1);
  }
  // with naming convention and backward compatibility
  fGainArray->AddAt(fGraphMIP,2);
  fGainArray->AddAt(fGraphCosmic,3);
  //
  // 3.) Add HV and PT correction parameterization which was used
  //
  AliTPCParam *param= AliTPCcalibDB::Instance()->GetParameters();
  if (param->GetGainSlopesHV())  fGainArray->AddLast(param->GetGainSlopesHV());
  if (param->GetGainSlopesPT())  fGainArray->AddLast(param->GetGainSlopesPT());
  //
  // 4.) Add the RecoParam and ClusterParam - for compatibility checks -deffrent sets of paramters can invalidate calibration 
  //
  AliTPCClusterParam *clParam =   AliTPCcalibDB::Instance()->GetClusterParam();
  TObjArray *recoParams = new TObjArray(4) ;
  for (Int_t i=0;i<4;i++) recoParams->AddAt(AliTPCcalibDB::Instance()->GetRecoParam(i),i);
  fGainArray->AddLast(clParam);
  fGainArray->AddLast(recoParams);
  //
  cout << "fGraphCosmic: " << fGraphCosmic << " fGraphMIP " << fGraphMIP << endl;
  return kTRUE;

}

//_____________________________________________________________________________
Bool_t AliTPCPreprocessorOffline::AnalyzeAttachment(Int_t startRunNumber, Int_t endRunNumber, Int_t minEntriesFit) {
  //
  // determine slope as a function of mean driftlength
  //
  if(!fGainMIP) return kFALSE;

  fGainMIP->GetHistGainTime()->GetAxis(5)->SetRangeUser(startRunNumber, endRunNumber);
  //
  fGainMIP->GetHistGainTime()->GetAxis(2)->SetRangeUser(1.51,2.49); // only beam data
  fGainMIP->GetHistGainTime()->GetAxis(4)->SetRangeUser(0.39,0.51); // only MIP pions
  //
  fGainMIP->GetHistGainTime()->GetAxis(3)->SetRangeUser(125,250);// only full tracking region (driftlength)
  fGainMIP->GetHistGainTime()->GetAxis(0)->SetRangeUser(1.5,3.5);// only full tracking region (driftlength)
  //
  TH3D * hist = fGainMIP->GetHistGainTime()->Projection(1, 0, 3);
  //
  Double_t *xvec = new Double_t[hist->GetNbinsX()];
  Double_t *yvec = new Double_t[hist->GetNbinsX()];
  Double_t *xerr = new Double_t[hist->GetNbinsX()];
  Double_t *yerr = new Double_t[hist->GetNbinsX()];
  Int_t counter  = 0;
  //
  for(Int_t i=1; i < hist->GetNbinsX(); i++) {
    Int_t nsum=0;
    Int_t imin   =  i;
    Int_t imax   =  i;    
    for (Int_t idelta=0; idelta<5; idelta++){
      //
      imin   =  TMath::Max(i-idelta,1);
      imax   =  TMath::Min(i+idelta,hist->GetNbinsX());
      nsum = TMath::Nint(hist->Integral(imin,imax,1,hist->GetNbinsY()-1,1,hist->GetNbinsZ()-1));
      //if (nsum==0) break;
      if (nsum>minEntriesFit) break;
    }
    if (nsum<minEntriesFit) continue;
    //
    fGainMIP->GetHistGainTime()->GetAxis(1)->SetRangeUser(hist->GetXaxis()->GetBinCenter(imin-1),hist->GetXaxis()->GetBinCenter(imax+1)); // define time range
    TH2D * histZdep = fGainMIP->GetHistGainTime()->Projection(0,3);
    TObjArray arr;
    histZdep->FitSlicesY(0,0,-1,0,"QNR",&arr);
    TH1D * driftDep = (TH1D*)arr.At(1);
    delete histZdep;
    //TGraphErrors * driftDep = AliTPCcalibBase::FitSlices(fGainMIP->GetHistGainTime(),0,3,100,1,0.,1);
    /*if (driftDep->GetN() < 4) {
      delete driftDep;
      continue;
      }*/
    //
    //TObjArray arr;
    //
    TF1 pol1("polynom1","pol1",125,240);
    //driftDep->Fit(&pol1,"QNRROB=0.8");
    driftDep->Fit(&pol1,"QNR");
    xvec[counter] = 0.5*(hist->GetXaxis()->GetBinCenter(imin-1)+hist->GetXaxis()->GetBinCenter(imax+1));
    yvec[counter] = pol1.GetParameter(1)/pol1.GetParameter(0);
    xerr[counter] = hist->GetXaxis()->GetBinCenter(imax+1)-hist->GetXaxis()->GetBinCenter(imin-1);
    yerr[counter] = pol1.GetParError(1)/pol1.GetParameter(0);
    counter++;
    //
    //delete driftDep;
  }
  //
  fGraphAttachmentMIP = new TGraphErrors(counter, xvec, yvec, xerr, yerr);
  if (fGraphAttachmentMIP) fGraphAttachmentMIP->SetName("TGRAPHERRORS_MEAN_ATTACHMENT_BEAM_ALL");// set proper names according to naming convention
  fGainArray->AddLast(fGraphAttachmentMIP);
  //
  delete [] xvec;
  delete [] yvec;
  delete [] xerr;
  delete [] yerr;
  delete hist;
  //
  if (counter < 1) return kFALSE;
  return kTRUE;

}

//_____________________________________________________________________________
Bool_t AliTPCPreprocessorOffline::AnalyzePadRegionGain(){
  //
  // Analyze gain for different pad regions - produce the calibration graphs 0,1,2
  // Uses the MIP pions integrated over all chambers, eta and multiplicity
  // TODO: o What happens in case the MIPs are not equally distributed over
  //         eta and multiplicity? In CPass0 neither will be flat.
  //         Can we fill the histograms flat in eta (and multiplicity)?
  //       o In principle one could normalise over eta, by first projecting in slices
  //         of eta, normalize each slice and then sum
  //
  if (!fGainMult) return kFALSE;

  THnSparseF *histPadEqual=fGainMult->GetHistPadEqual();
//   histPadEqual->GetAxis(4)->SetRangeUser(.35,.65);
//   TH2D * histQmaxTmp = (TH2D*) histPadEqual->Projection(0,2);
//   TH2D * histQtotTmp = (TH2D*) histPadEqual->Projection(1,2);
//   histQmaxTmp->SetDirectory(0);
//   histQtotTmp->SetDirectory(0);
//   histPadEqual->GetAxis(4)->SetRangeUser(-.65,-.35);
//   TH2D * histQmax = (TH2D*) histPadEqual->Projection(0,2);
//   TH2D * histQtot = (TH2D*) histPadEqual->Projection(1,2);
//   histQmax->Add(histQmaxTmp);
//   histQtot->Add(histQtotTmp);
//   delete histQmaxTmp;
//   delete histQtotTmp;
  TH2 * histQmax = AliTPCcalibBase::NormalizedProjection(histPadEqual,0,2,4);
  TH2 * histQtot = AliTPCcalibBase::NormalizedProjection(histPadEqual,1,2,4);

// === old part
//   TObjArray arr;
//   histQmax->FitSlicesY(0,0,-1,0,"QNR",&arr);
//   Double_t xMax[3] = {0,1,2};
//   Double_t yMax[3]    = {((TH1D*)arr.At(1))->GetBinContent(1),
//                          ((TH1D*)arr.At(1))->GetBinContent(2),
//                          ((TH1D*)arr.At(1))->GetBinContent(3)};
//   Double_t yMaxErr[3] = {((TH1D*)arr.At(1))->GetBinError(1),
//                          ((TH1D*)arr.At(1))->GetBinError(2),
//                          ((TH1D*)arr.At(1))->GetBinError(3)};
//   //
//   histQtot->FitSlicesY(0,0,-1,0,"QNR",&arr);
//   Double_t xTot[3] = {0,1,2};
//   Double_t yTot[3]    = {((TH1D*)arr.At(1))->GetBinContent(1),
//                          ((TH1D*)arr.At(1))->GetBinContent(2),
//                          ((TH1D*)arr.At(1))->GetBinContent(3)};
//   Double_t yTotErr[3] = {((TH1D*)arr.At(1))->GetBinError(1),
//                          ((TH1D*)arr.At(1))->GetBinError(2),
//                          ((TH1D*)arr.At(1))->GetBinError(3)};
//   Double_t truncationFactor=AliTPCcalibGainMult::GetTruncatedMeanPosition(yTot[0],yTot[1],yTot[2],1000);
//   for (Int_t i=0;i<3; i++){
//     yMax[i]*=truncationFactor;
//     yTot[i]*=truncationFactor;
//   }
//
//   TGraphErrors * fitPadRegionQmax = new TGraphErrors(3, xMax, yMax, 0, yMaxErr);
//   TGraphErrors * fitPadRegionQtot = new TGraphErrors(3, xTot, yTot, 0, yTotErr);
// ^^^ End old part

//   histQmax->GetXaxis()->SetRange(1,3);
//   histQtot->GetXaxis()->SetRange(1,3);
  TGraphErrors *fitPadRegionQmax = AliTPCcalibBase::FitSlices(histQmax,200,1,.15,.85);
  TGraphErrors *fitPadRegionQtot = AliTPCcalibBase::FitSlices(histQtot,200,1,.15,.85);
  histQmax->GetXaxis()->SetRange(0,-1);
  histQtot->GetXaxis()->SetRange(0,-1);
  fitPadRegionQmax->RemovePoint(3);
  fitPadRegionQtot->RemovePoint(3);

  Double_t *yMax =fitPadRegionQmax->GetY();
  Double_t *yTot =fitPadRegionQtot->GetY();
  Double_t *eyMax=fitPadRegionQmax->GetEY();
  Double_t *eyTot=fitPadRegionQtot->GetEY();
  Double_t truncationFactorMax=AliTPCcalibGainMult::GetTruncatedMeanPosition(yMax[0],yMax[1],yMax[2],1000);
  Double_t truncationFactorTot=AliTPCcalibGainMult::GetTruncatedMeanPosition(yTot[0],yTot[1],yTot[2],1000);

  if (truncationFactorMax>0) {
    for (Int_t i=0;i<3; i++){
      yMax[i] /=truncationFactorMax;
      eyMax[i]/=truncationFactorMax;
    }
  }
  else {
    AliError("truncationFactorMax<=0, could not normalize");
  }

  if (truncationFactorTot>0) {
    for (Int_t i=0;i<3; i++){
      yTot[i] /=truncationFactorTot;
      eyTot[i]/=truncationFactorTot;
    }
  }
  else {
    AliError("truncationFactorTot<=0, could not normalize");
  }

  //
  fitPadRegionQtot->SetName("TGRAPHERRORS_MEANQTOT_PADREGIONGAIN_BEAM_ALL");// set proper names according to naming convention
  fitPadRegionQmax->SetName("TGRAPHERRORS_MEANQMAX_PADREGIONGAIN_BEAM_ALL");// set proper names according to naming convention
  //
  fGainArray->AddLast(fitPadRegionQtot);
  fGainArray->AddLast(fitPadRegionQmax);

  // ---| QA histograms |---
  if (fArrQAhist) {
    // --- Qmax ---
    histQmax->SetName("TGRAPHERRORS_MEANQMAX_PADREGIONGAIN_BEAM_ALL_QA");
    histQmax->SetTitle("Pad region calibration Q_{max}; pad region;d#it{E}/d#it{x}_{Qmax} (arb. unit)");
    fArrQAhist->Add(histQmax);

    // --- Qtot ---
    histQtot->SetName("TGRAPHERRORS_MEANQTOT_PADREGIONGAIN_BEAM_ALL_QA");
    histQtot->SetTitle("Pad region calibration Q_{tot}; pad region;d#it{E}/d#it{x}_{Qtot} (arb. unit)");
    fArrQAhist->Add(histQtot);

  } else {
    delete histQmax;
    delete histQtot;
  }

  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliTPCPreprocessorOffline::AnalyzeGainDipAngle(Int_t padRegion)  {
  //
  // Analyze gain as a function of multiplicity and produce calibration graphs
  // padRegion -- 0: short, 1: medium, 2: long, 3: absolute calibration of full track
  // Uses the MIP pions integrated over all chambers and multiplicity
  // Noramlisation is done to the mean of the slices fit (i.e. tan(lambda)=0.5)
  //
  // TODO: What happens in case the MIPs are not equally distributed over
  //       multiplicity? In CPass0 it will not be flat.
  //       Can we fill the histograms flat in multiplicity?
  //
  Int_t kMarkers[10]={25,24,20,21,22};
  Int_t kColors[10]={1,2,4,3,6};
  if (!fGainMult) return kFALSE;
  if (!(fGainMult->GetHistTopology())) return kFALSE;
  const Double_t kMinStat=100;
  //
  // "dEdxRatioMax","dEdxRatioTot","padType","mult","driftlength"
  TObjArray arrMax;
  TObjArray arrTot;
  //
  TH2D * histQmax = 0x0;
  TH2D * histQtot = 0x0;
  fGainMult->GetHistPadEqual()->GetAxis(4)->SetRangeUser(-0.85,0.85);
  fGainMult->GetHistTopology()->GetAxis(2)->SetRangeUser(-0.85,0.85);
  if (padRegion < 3) {
    fGainMult->GetHistPadEqual()->GetAxis(2)->SetRangeUser(padRegion,padRegion); // short,medium,long
    histQmax = (TH2D*) fGainMult->GetHistPadEqual()->Projection(0,4);
    histQtot = (TH2D*) fGainMult->GetHistPadEqual()->Projection(1,4);
  } else {
    fGainMult->GetHistTopology()->GetAxis(1)->SetRangeUser(1,1); //Qmax
    histQmax = (TH2D*) fGainMult->GetHistTopology()->Projection(0,2);
    histQmax->SetName("fGainMult_GetHistPadEqual_11");
    fGainMult->GetHistTopology()->GetAxis(1)->SetRangeUser(0,0); //Qtot
    histQtot = (TH2D*) fGainMult->GetHistTopology()->Projection(0,2);
    histQtot->SetName("fGainMult_GetHistPadEqual_00");
  }
  //  
  
  if (histQmax->GetEntries()<=kMinStat || histQtot->GetEntries()<=kMinStat) {
    AliError(Form("hisQtot.GetEntries()=%f",histQtot->GetEntries()));
    AliError(Form("hisQmax.GetEntries()=%f",histQmax->GetEntries()));
    return kFALSE;
  }

//   TGraphErrors * graphMax = TStatToolkit::MakeStat1D( histQmax,0,0.8,4,kMarkers[padRegion],kColors[padRegion]);
//   TGraphErrors * graphTot = TStatToolkit::MakeStat1D( histQtot,0,0.8,4,kMarkers[padRegion],kColors[padRegion]);
  TGraphErrors * graphMax = TStatToolkit::MakeStat1D( histQmax,0,0.9,6,kMarkers[padRegion],kColors[padRegion]);
  TGraphErrors * graphTot = TStatToolkit::MakeStat1D( histQtot,0,0.9,6,kMarkers[padRegion],kColors[padRegion]);

  //
  const char* names[4]={"SHORT","MEDIUM","LONG","ABSOLUTE"};
  //
  const Double_t meanMax = TMath::Mean(graphMax->GetN(), graphMax->GetY());
  const Double_t meanTot = TMath::Mean(graphTot->GetN(), graphTot->GetY());

  // ---| QA histograms |---
  if (fArrQAhist) {
    // --- Qmax ---
    histQmax->SetName(Form("TGRAPHERRORS_QMAX_DIPANGLE_%s_BEAM_ALL_QA",names[padRegion]));
    histQmax->SetTitle(Form("tan(#lambda) calibration Q_{max} %s; tan(#lambda);d#it{E}/d#it{x}_{Qmax} (arb. unit)",names[padRegion]));
    fArrQAhist->Add(histQmax);

    // --- Qtot ---
    histQtot->SetName(Form("TGRAPHERRORS_QTOT_DIPANGLE_%s_BEAM_ALL_QA",names[padRegion]));
    histQtot->SetTitle(Form("tan(#lambda) calibration Q_{tot} %s; tan(#lambda);d#it{E}/d#it{x}_{Qtot} (arb. unit)",names[padRegion]));
    fArrQAhist->Add(histQtot);

    // ---| scale to mean multiplicity |---
    if(fNormaliseQA && meanMax>0) {
      TAxis *a=histQmax->GetYaxis();
      a->SetLimits(a->GetXmin()/meanMax, a->GetXmax()/meanMax);
    }
    if(fNormaliseQA && meanTot>0) {
      TAxis *a=histQtot->GetYaxis();
      a->SetLimits(a->GetXmin()/meanTot, a->GetXmax()/meanTot);
    }

  } else {
    delete histQmax;
    delete histQtot;
  }

  //
  if (meanMax<=0 || meanTot<=0){
    AliError(Form("meanMax=%f",meanMax));
    AliError(Form("meanTot=%f",meanTot));
    return kFALSE;
  }
  //
  for (Int_t ipoint=0; ipoint<graphMax->GetN(); ipoint++) {
    graphMax->GetY()[ipoint]/=meanMax;
    graphMax->GetEY()[ipoint]/=meanMax;
  }
  for (Int_t ipoint=0; ipoint<graphTot->GetN(); ipoint++) {
    graphTot->GetY()[ipoint]/=meanTot;
    graphTot->GetEY()[ipoint]/=meanTot;
  }
  //
  graphMax->SetNameTitle(Form("TGRAPHERRORS_QMAX_DIPANGLE_%s_BEAM_ALL",names[padRegion]),
			Form("TGRAPHERRORS_QMAX_DIPANGLE_%s_BEAM_ALL",names[padRegion]));
  graphTot->SetNameTitle(Form("TGRAPHERRORS_QTOT_DIPANGLE_%s_BEAM_ALL",names[padRegion]),
			Form("TGRAPHERRORS_QTOT_DIPANGLE_%s_BEAM_ALL",names[padRegion]));
  //
  fGainArray->AddLast(graphMax);
  fGainArray->AddLast(graphTot);
  //
  // Normalization to 1 (mean of the graph.fY --> 1)
  //
  TF1 * funMax= new TF1("","1++abs(x)++abs(x*x)");
  TF1 * funTot= new TF1("","1++abs(x)++abs(x*x)");
  graphMax->Fit(funMax,"w","rob=0.9",-0.8,0.8);
  graphTot->Fit(funTot,"w","rob=0.9",-0.8,0.8);
  funMax->SetNameTitle(Form("TF1_QMAX_DIPANGLE_%s_BEAM_ALL",names[padRegion]),
			Form("TF1_QMAX_DIPANGLE_%s_BEAM_ALL",names[padRegion]));
  funTot->SetNameTitle(Form("TF1_QTOT_DIPANGLE_%s_BEAM_ALL",names[padRegion]),
			Form("TF1_QTOT_DIPANGLE_%s_BEAM_ALL",names[padRegion]));

  //
  fGainArray->AddLast(funMax);
  fGainArray->AddLast(funTot);
  //
  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliTPCPreprocessorOffline::AnalyzeGainMultiplicity() {
  //
  // Analyze gain as a function of multiplicity and produce calibration graphs
  // Uses the MIP pions integrated over all chambers and eta
  // Noramlisation is done to the mean of the slices fit (i.e. sector average)
  //
  // TODO: o What happens in case the MIPs are not flat in eta?
  //         In CPass0 it will not be flat.
  //         Can we fill the histograms flat in multiplicity?
  //       o Is it better to use the median in case one chamber fit fails?
  //       o The gain mult does not have tan(lambda) as dim
  //
  if (!fGainMult) return kFALSE;
  fGainMult->GetHistGainMult()->GetAxis(3)->SetRangeUser(3,3);
  TH2D * histMultMax = fGainMult->GetHistGainMult()->Projection(0,4);
  TH2D * histMultTot = fGainMult->GetHistGainMult()->Projection(1,4);
  histMultMax->RebinX(4);
  histMultTot->RebinX(4);
  //
  TObjArray arrMax;
  TObjArray arrTot;
  TF1 fitGaus("fitGaus","gaus(0)",histMultMax->GetYaxis()->GetXmin(),histMultMax->GetYaxis()->GetXmax());
  fitGaus.SetParameters(histMultMax->GetEntries()/10., histMultMax->GetMean(2), TMath::Sqrt(TMath::Abs(histMultMax->GetMean(2))));
  histMultMax->FitSlicesY(&fitGaus,0,-1,1,"QNRB",&arrMax);
  fitGaus.SetParameters(histMultTot->GetEntries()/10., histMultTot->GetMean(2), TMath::Sqrt(TMath::Abs(histMultTot->GetMean(2))));
  histMultTot->FitSlicesY(&fitGaus,0,-1,1,"QNRB",&arrTot);
  //
  TH1D * meanMax = (TH1D*)arrMax.At(1);
  TH1D * meanTot = (TH1D*)arrTot.At(1);
  Float_t meanMult = histMultMax->GetMean();
  const Double_t qMaxCont=meanMax->GetBinContent(meanMax->FindBin(meanMult));
  const Double_t qTotCont=meanTot->GetBinContent(meanTot->FindBin(meanMult));

  // ---| QA histograms |---
  if (fArrQAhist) {
    // --- Qmax ---
    histMultMax->SetName("TGRAPHERRORS_MEANQMAX_MULTIPLICITYDEPENDENCE_BEAM_ALL_QA");
    histMultMax->SetTitle("Multiplicity correction Q_{max};#ESD tracks;d#it{E}/d#it{x}_{Qmax} (arb. unit)");
    fArrQAhist->Add(histMultMax);

    // --- Qtot ---
    histMultTot->SetName("TGRAPHERRORS_MEANQTOT_MULTIPLICITYDEPENDENCE_BEAM_ALL_QA");
    histMultTot->SetTitle("Multiplicity correction Q_{tot};#ESD tracks;d#it{E}/d#it{x}_{Qtot} (arb. unit)");
    fArrQAhist->Add(histMultTot);

    // ---| scale to mean multiplicity |---
    if(fNormaliseQA && qMaxCont>0) {
      TAxis *a=histMultMax->GetYaxis();
      a->SetLimits(a->GetXmin()/qMaxCont, a->GetXmax()/qMaxCont);
    }
    if(fNormaliseQA && qTotCont>0) {
      TAxis *a=histMultTot->GetYaxis();
      a->SetLimits(a->GetXmin()/qTotCont, a->GetXmax()/qTotCont);
    }

  } else {
    delete histMultMax;
    delete histMultTot;
  }

  // ---| scale to mean multiplicity |---
  if(qMaxCont) {
    meanMax->Scale(1./qMaxCont);
  }
  else {
   return kFALSE;
  }
  if(qTotCont) {
    meanTot->Scale(1./qTotCont);
  }
  else {
   return kFALSE;
  }
  Float_t xMultMax[50];
  Float_t yMultMax[50];
  Float_t yMultErrMax[50];
  Float_t xMultTot[50];
  Float_t yMultTot[50];
  Float_t yMultErrTot[50];
  //
  Int_t nCountMax = 0;
  for(Int_t iBin = 1; iBin < meanMax->GetXaxis()->GetNbins(); iBin++) {
    Float_t yValMax = meanMax->GetBinContent(iBin);
    if (yValMax < 0.7) continue;
    if (yValMax > 1.3) continue;
    if (meanMax->GetBinError(iBin)/yValMax > 0.01) continue;
    xMultMax[nCountMax] = meanMax->GetXaxis()->GetBinCenter(iBin);
    yMultMax[nCountMax] = yValMax;
    yMultErrMax[nCountMax] = meanMax->GetBinError(iBin);
    nCountMax++;
  }
  //
  if (nCountMax < 10) return kFALSE;
  TGraphErrors * fitMultMax = new TGraphErrors(nCountMax, xMultMax, yMultMax, 0, yMultErrMax);
  fitMultMax->SetName("TGRAPHERRORS_MEANQMAX_MULTIPLICITYDEPENDENCE_BEAM_ALL");
  //
  Int_t nCountTot = 0;
  for(Int_t iBin = 1; iBin < meanTot->GetXaxis()->GetNbins(); iBin++) {
    Float_t yValTot = meanTot->GetBinContent(iBin);
    if (yValTot < 0.7) continue;
    if (yValTot > 1.3) continue;
    if (meanTot->GetBinError(iBin)/yValTot > 0.1) continue;
    xMultTot[nCountTot] = meanTot->GetXaxis()->GetBinCenter(iBin);
    yMultTot[nCountTot] = yValTot;
    yMultErrTot[nCountTot] = meanTot->GetBinError(iBin);
    nCountTot++;
  }
  //
  if (nCountTot < 10) return kFALSE;
  TGraphErrors *  fitMultTot = new TGraphErrors(nCountTot, xMultTot, yMultTot, 0, yMultErrTot);
  fitMultTot->SetName("TGRAPHERRORS_MEANQTOT_MULTIPLICITYDEPENDENCE_BEAM_ALL");
  //
  fGainArray->AddLast(fitMultMax);
  fGainArray->AddLast(fitMultTot);
  //
  return kTRUE;

}

//_____________________________________________________________________________
Bool_t AliTPCPreprocessorOffline::AnalyzeGainChamberByChamber(){
  //
  // get chamber by chamber gain
  //
  if (!fGainMult) return kFALSE;
  TGraphErrors *grShort  = fGainMult->GetGainPerChamberRobust(0, kFALSE, fArrQAhist, fNormaliseQA);
  TGraphErrors *grMedium = fGainMult->GetGainPerChamberRobust(1, kFALSE, fArrQAhist, fNormaliseQA);
  TGraphErrors *grLong   = fGainMult->GetGainPerChamberRobust(2, kFALSE, fArrQAhist, fNormaliseQA);
  if (grShort==0x0 || grMedium==0x0 || grLong==0x0) {
    delete grShort;
    delete grMedium;
    delete grLong;
    return kFALSE;
  }

  fGainArray->AddLast(grShort);
  fGainArray->AddLast(grMedium);
  fGainArray->AddLast(grLong);

  return kTRUE;
}

//_____________________________________________________________________________
void AliTPCPreprocessorOffline::UpdateOCDBGain(Int_t startRunNumber, Int_t endRunNumber, AliCDBStorage *fullStorage, AliCDBStorage* residualStorage/*=0x0*/)
{
  //
  // Update OCDB entry
  //
  AliCDBMetaData *metaData= new AliCDBMetaData();
  metaData->SetObjectClassName("TObjArray");
  metaData->SetResponsible("Jens Wiechula");
  metaData->SetBeamPeriod(1);
  metaData->SetAliRootVersion("05-24-00"); //root version
  metaData->SetComment("Calibration of the time dependence of the gain due to pressure and temperature changes.");
  AliCDBId id1("TPC/Calib/TimeGain", startRunNumber, endRunNumber);
  if (fGainCalibrationType==kFullGainCalib) {
    AliInfoF("Writing gain calibration object to the full storage: %s", fullStorage->GetURI().Data());
    fullStorage->Put(fGainArray, id1, metaData);
  }
  else if (fGainCalibrationType==kResidualGainQA) {
    AliInfoF("Writing gain calibration object to the residual storage: %s", residualStorage->GetURI().Data());
    residualStorage->Put(fGainArray, id1, metaData);
  }
  else if (fGainCalibrationType==kCombinedGainCalib) {
    if (residualStorage){
      AliInfoF("Writing residual gain calibration object to the residual storage: %s", residualStorage->GetURI().Data());
      residualStorage->Put(fGainArray, id1, metaData);
    }
    else {
      AliError("No residual storage set, but calibration type is combined + residual QA");
    }

    // write the combined calibration after the residual calibration to make sure this is the
    //   latest object in case fullStorage and residualStorage are identical
    AliInfoF("Writing combined gain calibration object to the full storage: %s", fullStorage->GetURI().Data());
    fullStorage->Put(fGainArrayCombined, id1, metaData);
  }
  else {
    AliFatalF("Unsupported gain calibration type: %d", Int_t(fGainCalibrationType));
  }
}

//_____________________________________________________________________________
void AliTPCPreprocessorOffline::MakeQAPlot(Float_t  FPtoMIPratio) {
  //
  // Make QA plot to visualize results
  //
  //
  //
  if (fGraphCosmic) {
    TCanvas * canvasCosmic = new TCanvas("gain Cosmic", "time dependent gain QA histogram cosmic");
    canvasCosmic->cd();
    TH2D * gainHistoCosmic = fGainCosmic->GetHistGainTime()->Projection(0,1);
    gainHistoCosmic->SetDirectory(0);
    gainHistoCosmic->SetName("GainHistoCosmic");
    gainHistoCosmic->GetXaxis()->SetTimeDisplay(kTRUE);
    gainHistoCosmic->GetXaxis()->SetTimeFormat("#splitline{%d/%m}{%H:%M}");
    gainHistoCosmic->Draw("colz");
    fGraphCosmic->SetMarkerStyle(25);
    fGraphCosmic->Draw("lp");
    fGraphCosmic->SetMarkerStyle(25);
    TGraph * grfFitCosmic = fFitCosmic->MakeGraph(fGraphCosmic->GetX()[0],fGraphCosmic->GetX()[fGraphCosmic->GetN()-1],50000,0);
    if (grfFitCosmic) {
      for(Int_t i=0; i < grfFitCosmic->GetN(); i++) {
 	grfFitCosmic->GetY()[i] *= FPtoMIPratio;	
      }
      for(Int_t i=0; i < fGraphCosmic->GetN(); i++) {
 	fGraphCosmic->GetY()[i] *= FPtoMIPratio;	
      }
    }
    fGraphCosmic->Draw("lp"); 
    if (grfFitCosmic) {
      grfFitCosmic->SetLineColor(2);
      grfFitCosmic->Draw("lu");
    }
    fGainArray->AddLast(gainHistoCosmic);
    //fGainArray->AddLast(canvasCosmic->Clone());
    delete canvasCosmic;    
  }
  if (fFitMIP) {
    TCanvas * canvasMIP = new TCanvas("gain MIP", "time dependent gain QA histogram MIP");
    canvasMIP->cd();
    TH2D * gainHistoMIP    = fGainMIP->GetHistGainTime()->Projection(0,1);
    gainHistoMIP->SetName("GainHistoCosmic");
    gainHistoMIP->SetDirectory(0);
    gainHistoMIP->GetXaxis()->SetTimeDisplay(kTRUE);
    gainHistoMIP->GetXaxis()->SetTimeFormat("#splitline{%d/%m}{%H:%M}");
    gainHistoMIP->Draw("colz");
    fGraphMIP->SetMarkerStyle(25);
    fGraphMIP->Draw("lp");
    TGraph * grfFitMIP = fFitMIP->MakeGraph(fGraphMIP->GetX()[0],fGraphMIP->GetX()[fGraphMIP->GetN()-1],50000,0);
    grfFitMIP->SetLineColor(2);
    grfFitMIP->Draw("lu");    
    fGainArray->AddLast(gainHistoMIP);
    //fGainArray->AddLast(canvasMIP->Clone());
    delete canvasMIP;  
    delete grfFitMIP;
  }  
}

//_____________________________________________________________________________
void AliTPCPreprocessorOffline::MakeFitTime(){
  //
  // make aligment fit - store results in the file
  //
  const Int_t kMinEntries=1000;
  MakeChainTime();
  MakePrimitivesTime();
  if (!fAlignTree) return;
  if (fAlignTree->GetEntries()<kMinEntries) return;
  fAlignTree->SetAlias("ptype","type");
  fAlignTree->SetAlias("hasITS","(1+0)");
  fAlignTree->SetAlias("dITS","1-2*(refX<40)");
  fAlignTree->SetAlias("isITS","refX>10");
  fAlignTree->SetAlias("isVertex","refX<10");
  // 
  Int_t  npointsMax=30000000;
  Double_t chi2=0;
  Int_t    npoints=0;
  TVectorD param;
  TMatrixD covar;

  TString fstringFast="";
  fstringFast+="FExBTwistX++";
  fstringFast+="FExBTwistY++";
  fstringFast+="FAlignRot0D++";
  fstringFast+="FAlignTrans0D++";
  fstringFast+="FAlignTrans1D++";
  //
  fstringFast+="hasITS*FAlignTrans0++";
  fstringFast+="hasITS*FAlignTrans1++";
  fstringFast+="hasITS*FAlignRot0++";
  fstringFast+="hasITS*FAlignRot1++";
  fstringFast+="hasITS*FAlignRot2++";
  //
  fstringFast+="dITS*FAlignTrans0++";
  fstringFast+="dITS*FAlignTrans1++";
  fstringFast+="dITS*FAlignRot0++";
  fstringFast+="dITS*FAlignRot1++";
  fstringFast+="dITS*FAlignRot2++";
  
  TCut cutFit="entries>10&&abs(mean)>0.00001&&rms>0";
  fAlignTree->SetAlias("err","rms");

  TString *strDeltaITS = TStatToolkit::FitPlaneConstrain(fAlignTree,"mean:err", fstringFast.Data(),cutFit, chi2,npoints,param,covar,-1,0, npointsMax, 1);
  TObjArray* tokArr = strDeltaITS->Tokenize("++");
  static bool verboseOutput = !(getenv("HLT_ONLINE_MODE") && strcmp(getenv("HLT_ONLINE_MODE"), "on") == 0);
  if (verboseOutput) tokArr->Print();
  delete tokArr;
  fAlignTree->SetAlias("fitYFast",strDeltaITS->Data());
  delete strDeltaITS;

  // 
  TVectorD paramC= param;
  TMatrixD covarC= covar;
  TStatToolkit::Constrain1D(fstringFast,"Trans0D",paramC,covarC,0, 0.1);
  TStatToolkit::Constrain1D(fstringFast,"Trans1D",paramC,covarC,0, 0.1);
  TStatToolkit::Constrain1D(fstringFast,"TwistX",paramC,covarC,0, 0.1);
  TStatToolkit::Constrain1D(fstringFast,"TwistY",paramC,covarC,0, 0.1);
  TString strFitConst=TStatToolkit::MakeFitString(fstringFast, paramC,covar);
  fAlignTree->SetAlias("fitYFastC",strFitConst.Data());
  CreateAlignTime(fstringFast,paramC);

}

//_____________________________________________________________________________
void AliTPCPreprocessorOffline::MakeChainTime(){
  //
  //
  //
  TFile f("CalibObjects.root");
  
  //  const char *cdtype[7]={"ITS","TRD","Vertex","TOF","TPC","TPC0","TPC1"};
  //const char *cptype[5]={"dy","dz","dsnp","dtheta","d1pt"}; 
  const char * hname[5]={"dy","dz","dsnp","dtheta","d1pt"};
  Int_t run=0;
  AliTPCcalibTime  *calibTime = 0;
  TObject* obj = dynamic_cast<TObject*>(f.Get("TPCCalib"));
  TObjArray * array = dynamic_cast<TObjArray*>(obj);
  TDirectory * dir = dynamic_cast<TDirectory*>(obj);
  if (dir) {
    calibTime = dynamic_cast<AliTPCcalibTime*>(dir->Get("calibTime"));
  }
  else if (array){
    calibTime = (AliTPCcalibTime *)array->FindObject("calibTime");
  } else {
    calibTime = (AliTPCcalibTime*)f.Get("calibTime");
  }
  if (!calibTime) return;
  AliTPCCorrectionFit::CreateAlignMaps(AliTracker::GetBz(), run);
  TTreeSRedirector *pcstream = new TTreeSRedirector("meanITSVertex.root");
  //
  Int_t ihis=0;
  THnSparse *his = calibTime->GetResHistoTPCITS(ihis);
  if (his){
    his->GetAxis(1)->SetRangeUser(-1.1,1.1);
    his->GetAxis(2)->SetRange(0,1000000);
    his->GetAxis(3)->SetRangeUser(-0.35,0.35);
    AliTPCCorrection::MakeDistortionMap(his,pcstream, Form("ITS%s",hname[ihis]),run,85.,ihis,3);
  }
  ihis=1;
  his = calibTime->GetResHistoTPCITS(ihis);
  if (his){
    his->GetAxis(1)->SetRangeUser(-1.1,1.1);
    his->GetAxis(2)->SetRange(0,1000000);
    his->GetAxis(3)->SetRangeUser(-0.35,0.35);
    AliTPCCorrection::MakeDistortionMap(his,pcstream, Form("ITS%s",hname[ihis]),run,85.,ihis,3);
  }
  ihis=2;
  his = calibTime->GetResHistoTPCITS(ihis);
  if (his){
    his->GetAxis(1)->SetRangeUser(-1.1,1.1);
    his->GetAxis(2)->SetRange(0,1000000);
    his->GetAxis(3)->SetRangeUser(-0.35,0.35);
    AliTPCCorrection::MakeDistortionMap(his,pcstream, Form("ITS%s",hname[ihis]),run,85.,ihis,3);
  }
  ihis=0;
  his = calibTime->GetResHistoTPCvertex(ihis);
  if (his){
    his->GetAxis(1)->SetRangeUser(-1.1,1.1);
    his->GetAxis(2)->SetRange(0,1000000);
    his->GetAxis(3)->SetRangeUser(-0.35,0.35);
    AliTPCCorrection::MakeDistortionMap(his,pcstream, Form("Vertex%s",hname[ihis]),run,0.,ihis,3);
  }
  ihis=2;
  his = calibTime->GetResHistoTPCvertex(ihis);
  if (his){
    his->GetAxis(1)->SetRangeUser(-1.1,1.1);
    his->GetAxis(2)->SetRange(0,1000000);
    his->GetAxis(3)->SetRangeUser(-0.35,0.35);
    AliTPCCorrection::MakeDistortionMap(his,pcstream, Form("Vertex%s",hname[ihis]),run,0.,ihis,3);

  }
  ihis=1;
  his = calibTime->GetResHistoTPCvertex(ihis);
  if (his){
    his->GetAxis(1)->SetRangeUser(-1.1,1.1);
    his->GetAxis(2)->SetRange(0,1000000);
    his->GetAxis(3)->SetRangeUser(-0.35,0.35);
    AliTPCCorrection::MakeDistortionMap(his,pcstream, Form("Vertex%s",hname[ihis]),run,0.,ihis,3);

  }
  ihis=0;
  his = calibTime->GetResHistoTPCTOF(ihis);
  if (his){
    his->GetAxis(1)->SetRangeUser(-1.1,1.1);
    his->GetAxis(2)->SetRange(0,1000000);
    his->GetAxis(3)->SetRangeUser(-0.35,0.35);
    AliTPCCorrection::MakeDistortionMap(his,pcstream, Form("TOF%s",hname[ihis]),run,0.,ihis,3);

  }
  ihis=0;
  his = calibTime->GetResHistoTPCTRD(ihis);
  if (his){
    his->GetAxis(1)->SetRangeUser(-1.1,1.1);
    his->GetAxis(2)->SetRange(0,1000000);
    his->GetAxis(3)->SetRangeUser(-0.35,0.35);
    AliTPCCorrection::MakeDistortionMap(his,pcstream, Form("TRD%s",hname[ihis]),run,0.,ihis,3);

  }
  if (dir || (!dir && !array)) { // the object is taken from a directory or a file
    delete calibTime;
  }
  delete pcstream;
}

//_____________________________________________________________________________
Double_t AliTPCPreprocessorOffline::EvalAt(Double_t phi, Double_t refX, Double_t theta, Int_t corr, Int_t ptype){
  //
  //
  //
  Double_t sector = 9*phi/TMath::Pi();
  if (sector<0) sector+=18;
  Double_t y85=AliTPCCorrection::GetCorrSector(sector,85,theta,1,corr);
  Double_t y245=AliTPCCorrection::GetCorrSector(sector,245,theta,1,corr);
  if (ptype==0) return y85+(y245-y85)*(refX-85.)/(245.-85.);
  if (ptype==2) return (y245-y85)/(245.-85.);
  return 0;
}

//_____________________________________________________________________________
Double_t AliTPCPreprocessorOffline::EvalAtPar(Double_t phi0, Double_t snp, Double_t refX, Double_t theta, Int_t corr, Int_t ptype, Int_t nsteps){
  //
  // Fit the distortion along the line with the parabolic model
  // Parameters:
  //  phi0 - phi at the entrance of the TPC
  //  snp  - local inclination angle at the entrance of the TPC
  //  refX - ref X where the distortion is evanluated
  //  theta
  //  
  static TLinearFitter fitter(3,"pol2"); 
  fitter.ClearPoints();
  if (nsteps<3) nsteps=3;
  Double_t deltaX=(245-85)/(nsteps);
  for (Int_t istep=0; istep<(nsteps+1); istep++){
    //
    Double_t localX =85.+deltaX*istep;
    Double_t localPhi=phi0+deltaX*snp*istep;
    Double_t sector = 9*localPhi/TMath::Pi();
    if (sector<0) sector+=18;
    Double_t y=AliTPCCorrection::GetCorrSector(sector,localX,theta,1,corr);
    Double_t dlocalX=AliTPCCorrection::GetCorrSector(sector,localX,theta,0,corr);
    Double_t x[1]={localX-dlocalX};
    fitter.AddPoint(x,y);
  }
  fitter.Eval();
  Double_t par[3];
  par[0]=fitter.GetParameter(0);
  par[1]=fitter.GetParameter(1);
  par[2]=fitter.GetParameter(2);

  if (ptype==0) return par[0]+par[1]*refX+par[2]*refX*refX;
  if (ptype==2) return par[1]+2*par[2]*refX;
  if (ptype==4) return par[2];
  return 0;
}

//_____________________________________________________________________________
void AliTPCPreprocessorOffline::MakePrimitivesTime(){
  //
  // Create primitive transformation to fit
  //
  fAlignTree=new TChain("fit","fit");
  fAlignTree->AddFile("meanITSVertex.root",10000000,"ITSdy");
  fAlignTree->AddFile("meanITSVertex.root",10000000,"ITSdsnp");
  fAlignTree->AddFile("meanITSVertex.root",10000000,"Vertexdy");
  fAlignTree->AddFile("meanITSVertex.root",10000000,"Vertexdsnp");
  // 
  AliTPCParam *param= AliTPCcalibDB::Instance()->GetParameters();
  Double_t bzField=AliTrackerBase::GetBz(); 
  Double_t vdrift = param->GetDriftV()/1000000.; // [cm/us]   // From dataBase: to be updated: per second (ideally)
  Double_t ezField = 400; // [V/cm]   // to be updated: never (hopefully)
  Double_t wtP = -10.0 * (bzField) * vdrift /  ezField ; 
  AliTPCExBTwist *fitExBTwistX= new  AliTPCExBTwist;
  AliTPCExBTwist *fitExBTwistY= new  AliTPCExBTwist;
  AliTPCCalibGlobalMisalignment *trans0   =new  AliTPCCalibGlobalMisalignment;
  AliTPCCalibGlobalMisalignment *trans1   =new  AliTPCCalibGlobalMisalignment;
  AliTPCCalibGlobalMisalignment *trans0D  =new  AliTPCCalibGlobalMisalignment;
  AliTPCCalibGlobalMisalignment *trans1D  =new  AliTPCCalibGlobalMisalignment;
  AliTPCCalibGlobalMisalignment *rot0     =new  AliTPCCalibGlobalMisalignment;
  AliTPCCalibGlobalMisalignment *rot1     =new  AliTPCCalibGlobalMisalignment;
  AliTPCCalibGlobalMisalignment *rot2     =new  AliTPCCalibGlobalMisalignment;
  AliTPCCalibGlobalMisalignment *rot3     =new  AliTPCCalibGlobalMisalignment;
  //
  //
  fitExBTwistX->SetXTwist(0.001);
  fitExBTwistX->SetOmegaTauT1T2(wtP,1,1);  
  //
  fitExBTwistY->SetYTwist(0.001);
  fitExBTwistY->SetOmegaTauT1T2(wtP,1,1);  
  //
  TGeoHMatrix *matrixRot = new TGeoHMatrix; 
  TGeoHMatrix *matrixX = new TGeoHMatrix; 
  TGeoHMatrix *matrixY = new TGeoHMatrix; 
  matrixX->SetDx(0.1);
  matrixY->SetDy(0.1);
  Double_t rotAngles0[9]={0};
  Double_t rotAngles1[9]={0};
  Double_t rotAngles2[9]={0};
  //
  Double_t rotAngles3[9]={0};

  rotAngles0[0]=1; rotAngles0[4]=1; rotAngles0[8]=1;
  rotAngles1[0]=1; rotAngles1[4]=1; rotAngles1[8]=1;
  rotAngles2[0]=1; rotAngles2[4]=1; rotAngles2[8]=1;
  rotAngles3[0]=1; rotAngles3[4]=1; rotAngles3[8]=1;

  rotAngles0[1]=-0.001;rotAngles0[3]=0.001;
  rotAngles1[5]=-0.001;rotAngles1[7]=0.001;
  rotAngles2[2]=0.001;rotAngles2[6]=-0.001;
  rotAngles3[1]=0.001;rotAngles3[3]=-0.001;
  matrixRot->SetRotation(rotAngles0);
  rot0->SetAlignGlobal(matrixRot);
  matrixRot->SetRotation(rotAngles1);
  rot1->SetAlignGlobal(matrixRot);
  matrixRot->SetRotation(rotAngles2);
  rot2->SetAlignGlobal(matrixRot); 
  matrixRot->SetRotation(rotAngles3);
  rot3->SetAlignGlobalDelta(matrixRot); 
  //
  trans0->SetAlignGlobal(matrixX);
  trans1->SetAlignGlobal(matrixY);
  trans0D->SetAlignGlobalDelta(matrixX);
  trans1D->SetAlignGlobalDelta(matrixY);
  fitExBTwistX->Init();
  fitExBTwistY->Init();
  //
  fitExBTwistX->AddVisualCorrection((AliTPCExBTwist*)(fitExBTwistX->Clone()),100);
  fitExBTwistY->AddVisualCorrection((AliTPCExBTwist*)(fitExBTwistY->Clone()),101);
  //
  fitExBTwistX->AddVisualCorrection((AliTPCExBTwist*)(rot0->Clone()),102);
  fitExBTwistY->AddVisualCorrection((AliTPCExBTwist*)(rot1->Clone()),103);
  fitExBTwistX->AddVisualCorrection((AliTPCExBTwist*)(rot2->Clone()),104);
  fitExBTwistX->AddVisualCorrection((AliTPCExBTwist*)(rot3->Clone()),105);

  fitExBTwistX->AddVisualCorrection((AliTPCExBTwist*)(trans0->Clone()),106);
  fitExBTwistX->AddVisualCorrection((AliTPCExBTwist*)(trans1->Clone()),107);
  fitExBTwistX->AddVisualCorrection((AliTPCExBTwist*)(trans0D->Clone()),108);
  fitExBTwistX->AddVisualCorrection((AliTPCExBTwist*)(trans1D->Clone()),109);
  //
  fAlignTree->SetAlias("FExBTwistX", "AliTPCPreprocessorOffline::EvalAt(phi,refX,theta,100,ptype)+0");
  fAlignTree->SetAlias("FExBTwistY","AliTPCPreprocessorOffline::EvalAt(phi,refX,theta,101,ptype)+0");
  fAlignTree->SetAlias("FAlignRot0","AliTPCPreprocessorOffline::EvalAt(phi,refX,theta,102,ptype)+0");
  fAlignTree->SetAlias("FAlignRot0D","AliTPCPreprocessorOffline::EvalAt(phi,refX,theta,105,ptype)+0");
  fAlignTree->SetAlias("FAlignRot1","AliTPCPreprocessorOffline::EvalAt(phi,refX,theta,103,ptype)+0");
  fAlignTree->SetAlias("FAlignRot2","AliTPCPreprocessorOffline::EvalAt(phi,refX,theta,104,ptype)+0");
  fAlignTree->SetAlias("FAlignTrans0","AliTPCPreprocessorOffline::EvalAt(phi,refX,theta,106,ptype)+0");
  fAlignTree->SetAlias("FAlignTrans1","AliTPCPreprocessorOffline::EvalAt(phi,refX,theta,107,ptype)+0");
  fAlignTree->SetAlias("FAlignTrans0D","AliTPCPreprocessorOffline::EvalAt(phi,refX,theta,108,ptype)+0");
  fAlignTree->SetAlias("FAlignTrans1D","AliTPCPreprocessorOffline::EvalAt(phi,refX,theta,109,ptype)+0");
  //
  // test fast function
  //
//   fAlignTree->Draw("FExBTwistX:ExBTwistX","isITS&&ptype==0&&abs(snp)<0.05","");
//   fAlignTree->Draw("FExBTwistY:ExBTwistY","isITS&&ptype==0&&abs(snp)<0.05","");
//   fAlignTree->Draw("FAlignRot0:alignRot0","isITS&&ptype==0&&abs(snp)<0.05","");
//   fAlignTree->Draw("FAlignRot1:alignRot1","isITS&&ptype==0&&abs(snp)<0.05","");
//   fAlignTree->Draw("FAlignRot2:alignRot2","isITS&&ptype==0&&abs(snp)<0.05","");
//   //
//   fAlignTree->Draw("FAlignTrans0:alignTrans0","isITS&&ptype==0&&abs(snp)<0.05","");
//   fAlignTree->Draw("FAlignTrans1:alignTrans1","isITS&&ptype==0&&abs(snp)<0.05","");

} 

//_____________________________________________________________________________
void AliTPCPreprocessorOffline::CreateAlignTime(TString fstring, TVectorD paramC){
  //
  //
  //
  //
  TGeoHMatrix *matrixDelta     = new TGeoHMatrix; 
  TGeoHMatrix *matrixGlobal    = new TGeoHMatrix; 
  Double_t rAngles[9];
  Int_t index=0;
  //
  index=TStatToolkit::GetFitIndex(fstring,"FAlignTrans0D");
  if (index>=0) matrixDelta->SetDx(paramC[index+1]*0.1);
  index=TStatToolkit::GetFitIndex(fstring,"FAlignTrans1D");
  if (index>=0) matrixDelta->SetDy(paramC[index+1]*0.1);
  rAngles[0]=1; rAngles[4]=1; rAngles[8]=1;
  index=TStatToolkit::GetFitIndex(fstring,"FAlignRot0D");
  rAngles[1]=-paramC[index+1]*0.001; rAngles[3]=paramC[index+1]*0.001;
  rAngles[5]=0; rAngles[7] =0;
  rAngles[2]=0; rAngles[6] =0;
  matrixDelta->SetRotation(rAngles);
  //
  //
  //
  index=TStatToolkit::GetFitIndex(fstring,"hasITS*FAlignTrans0");
  if (index>=0) matrixGlobal->SetDx(paramC[index+1]*0.1);
  index=TStatToolkit::GetFitIndex(fstring,"hasITS*FAlignTrans1");
  if (index>=0) matrixGlobal->SetDy(paramC[index+1]*0.1);
  rAngles[0]=1; rAngles[4]=1; rAngles[8]=1;
  index=TStatToolkit::GetFitIndex(fstring,"hasITS*FAlignRot0");
  rAngles[1]=-paramC[index+1]*0.001; rAngles[3]=paramC[index+1]*0.001;
  index=TStatToolkit::GetFitIndex(fstring,"hasITS*FAlignRot1");  
  rAngles[5]=-paramC[index+1]*0.001; rAngles[7]=paramC[index+1]*0.001;
  index=TStatToolkit::GetFitIndex(fstring,"hasITS*FAlignRot2");  
  rAngles[2]=paramC[index+1]*0.001; rAngles[6] =-paramC[index+1]*0.001;
  matrixGlobal->SetRotation(rAngles);
  //
  AliTPCCalibGlobalMisalignment *fitAlignTime  =0;
  fitAlignTime  =new  AliTPCCalibGlobalMisalignment;
  fitAlignTime->SetName("FitAlignTime");
  fitAlignTime->SetTitle("FitAlignTime");
  fitAlignTime->SetAlignGlobalDelta(matrixDelta);
  fitAlignTime->SetAlignGlobal(matrixGlobal);
  //
  AliTPCExBTwist * fitExBTwist= new  AliTPCExBTwist;
  Int_t indexX=TStatToolkit::GetFitIndex(fstring,"ExBTwistX");
  Int_t indexY=TStatToolkit::GetFitIndex(fstring,"ExBTwistY");  
  fitExBTwist->SetXTwist(0.001*paramC[indexX+1]);  // 1 mrad twist in x
  fitExBTwist->SetYTwist(0.001*paramC[indexY+1]);  // 1 mrad twist in x
  fitExBTwist->SetName("FitExBTwistTime");
  fitExBTwist->SetTitle("FitExBTwistTime"); 
  AliTPCParam *param= AliTPCcalibDB::Instance()->GetParameters();
  Double_t bzField=AliTrackerBase::GetBz();
  Double_t vdrift = param->GetDriftV()/1000000.; // [cm/us]   // From dataBase: to be updated: per second (ideally)

  Double_t ezField = 400; // [V/cm]   // to be updated: never (hopefully)
  Double_t wt = -10.0 * (bzField) * vdrift /  ezField ; 
  //
  fitExBTwist->SetOmegaTauT1T2(wt,1,1);  
  fitExBTwist->Init();  

  AliTPCComposedCorrection *corrTime =  new AliTPCComposedCorrection;
  TObjArray *arr = new TObjArray;
  corrTime->SetCorrections(arr);
  
  corrTime->GetCorrections()->Add(fitExBTwist);
  corrTime->GetCorrections()->Add(fitAlignTime);
  corrTime->SetName("FitCorrectionTime");
  corrTime->SetTitle("FitCorrectionTime");

  fitExBTwist->AddVisualCorrection((AliTPCExBTwist*)(fitExBTwist->Clone()),1001);
  fitAlignTime->AddVisualCorrection((AliTPCExBTwist*)(fitAlignTime->Clone()),1002);
  fitAlignTime->AddVisualCorrection((AliTPCExBTwist*)(corrTime->Clone()),1003);
  
  
  fAlignTree->SetAlias("ExBTwistTime","AliTPCPreprocessorOffline::EvalAt(phi,refX,theta,1001,ptype)+0");
  fAlignTree->SetAlias("AlignTime","AliTPCPreprocessorOffline::EvalAt(phi,refX,theta,1002,ptype)+0");
  fAlignTree->SetAlias("FitCorrectionTime","AliTPCPreprocessorOffline::EvalAt(phi,refX,theta,1003,ptype)+0");


  TFile *f = new TFile("fitITSVertex.root","update");
  corrTime->Write("FitCorrectionTime");
  f->Close();
}

//_____________________________________________________________________________
Int_t AliTPCPreprocessorOffline::GetStatus()
{
  //get the calibration status
  // 0 means OK
  // positive numbers invalidate for unknown reasons.
  // negative numbers invalidate with a known reason (e.g. low statistics).
  // the returned integer has one bit set for every component that failed.

  Bool_t enoughStatistics = (fNtracksVdrift>fMinTracksVdrift && fNeventsVdrift>fMinEventsVdrift);
  
  if (!enoughStatistics) 
  {
    fCalibrationStatus=-TMath::Abs(fCalibrationStatus);
  }

  return fCalibrationStatus;
}

//_____________________________________________________________________________
void AliTPCPreprocessorOffline::FillQA(Bool_t qa, Bool_t norm/*=kTRUE*/)
{
  // setup QA histogram array
  if (qa && !fArrQAhist) {
    fArrQAhist = new TObjArray;
    fArrQAhist->SetOwner();
  } else {
    delete fArrQAhist;
    fArrQAhist = 0x0;
  }

  fNormaliseQA=norm;
}

//_____________________________________________________________________________
void AliTPCPreprocessorOffline::MakeQAPlotsGain(TString outputDirectory/*=""*/, TString fileTypes/*="png"*/)
{
  // Draw QA histograms, one per file
  // if outputDirectory is non empty, the QA canvases will be written to the output directory
  // fileTypes can contain output formats, comma separated. Default is one png per canvas.
  //           also recognized jpg, gif, root.
  //           in case of root, all canvases are written to one root file
  if (!fArrQAhist) return;

  TDirectory *dir=gDirectory;
  TFile *f=0x0;
  if (fileTypes.Contains("root")) {
    f=new TFile(TString::Format("%s/GainQA.root", outputDirectory.Data()),"recreate");
  }

  const Int_t ntypes=5;
  const TString ftypes[ntypes]={"png","jpg","gif","pdf","eps"};

  for (Int_t ihist=0; ihist<fArrQAhist->GetEntriesFast(); ++ihist) {
    dir->cd();
    TH2 *h = static_cast<TH2*>(fArrQAhist->UncheckedAt(ihist));
    if (!h) continue;
    TString histName = h->GetName();
    TString canvName=histName;
    canvName.Prepend("c_");
    TCanvas *c = static_cast<TCanvas*>(gROOT->GetListOfCanvases()->FindObject(canvName));
    if (!c) {
      c=new TCanvas(canvName, canvName, 700,500);
    }
    c->Clear();
    c->SetLogz();

    h->Draw("colz");

    // ---| check for derived histogram and draw it on top |---
    histName.ReplaceAll("_QA","");
    TObject *o = fGainArray->FindObject(histName);
    if (o) {
      o->Draw("same");
    }

    // ---| save output |---
    if (!outputDirectory.IsNull()) {
      // --- loop over file types ---
      for (Int_t itype=0; itype<ntypes; ++itype) {
        if (fileTypes.Contains(ftypes[itype])) {
          c->SaveAs(TString::Format("%s/%s.%s", outputDirectory.Data(), c->GetName(), ftypes[itype].Data()));
        }
      }

      // --- save to file if requested ---
      if (f) {
        f->cd();
        c->Write();
        dir->cd();
      }
    }
  }

  delete f;
}


void AliTPCPreprocessorOffline::ScaleY(TGraphErrors *graph, Double_t normval)
{
  for (Int_t ipoint=0; ipoint<graph->GetN(); ipoint++) {
    graph->GetY()[ipoint]*=normval;
    graph->GetEY()[ipoint]*=normval;
  }
}

Bool_t AliTPCPreprocessorOffline::NormaliseYToMean(TGraphErrors *graph)
{
  const Double_t mean = TMath::Mean(graph->GetN(), graph->GetY());

  if (mean<=0){
    AliError(Form("mean=%f",mean));
    return kFALSE;
  }

  ScaleY(graph, 1./mean);
}

Bool_t AliTPCPreprocessorOffline::NormaliseYToWeightedMeandEdx(TGraphErrors *graph)
{
  const Double_t *y=graph->GetY();
  Double_t scaleFactor=(63.*y[0] + 64.*y[1] + 32*y[2])/159.;

  if (scaleFactor<=0){
    AliError(Form("scaleFactor=%f",scaleFactor));
    return kFALSE;
  }

  ScaleY(graph, 1./scaleFactor);
}

Bool_t AliTPCPreprocessorOffline::NormaliseYToTruncateddEdx(TGraphErrors *graph)
{
  const Double_t *y=graph->GetY();
  Double_t truncationFactor=AliTPCcalibGainMult::GetTruncatedMeanPosition(y[0],y[1],y[2],1000);

  if (truncationFactor<=0){
    AliError(Form("truncationFactor=%f",truncationFactor));
    return kFALSE;
  }

  ScaleY(graph, 1./truncationFactor);
}

/*
  Short sequence to acces the calbration entry:
  TFile *f = TFile::Open("CalibObjects.root");
  AliTPCcalibGainMult      * fGainMult = (AliTPCcalibGainMult      *)f->Get("TPCCalib/calibGainMult");
  
 
*/
