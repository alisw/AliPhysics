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
 *                                                                        *
 **************************************************************************/

/********************************************************************************
*                                                                               *
*   AliLHCData: summary of the LHC related information from LHC DIP.            *
*   Created from the TMap provided by the AliLHCReader with optional beginning  *
*                                                                               *
*   The data are (wrapped in the AliLHCDipValT):                                *
*   made of TimeStamp (double) and array of values                              *
*                                                                               *
*   Multiple entries for each type of data are possible. To obtaine number of   *
*   records (with distinct timestamp) for give type od records use:             *
*   int GetNBunchConfigMeasured(int beam) (with beam=0,1) etc.                  *
*                                                                               *
*   To get i-th entry, use brec= AliLHCDipValI* GetBunchConfigMeasured(bm,i);   *
*   Note: exact type of templated AliLHCDipValT pointer depends on the record   *
*   type, concult getters to know it.                                           *
*                                                                               *
*   Then, once the pointer is obtained, details can be accessed:                *
*   int nBunches = brec->GetSize();                                             *
*   for (int i=0;i<nBunches;i++) printf("Bunch#%d: %d\n",i,(*brec)[i]);         *
*                                                                               *
*                                                                               *
*   Author: ruben.shahoyan@cern.ch                                              *
*                                                                               *
********************************************************************************/

#include "AliLHCData.h"
#include "TMap.h"
#include "AliDCSArray.h"
#include "AliLHCReader.h"
#include "AliTriggerBCMask.h"
#include <TString.h>
#include <TObjArray.h>
#include <TGraph.h>

ClassImp(AliLHCData)

const Char_t* AliLHCData::fgkDCSNames[] = {
  "LHC_IntensityBeam%1d_totalIntensity"
  ,"LHC_BeamIntensityPerBunchBeam%1d_averageBeamIntensity"
  ,"LHC_BeamIntensityPerBunchBeam%1d_Average_BunchIntensities"
  //
  ,"LHC_LumAverageBRANB_4%c2_acqMode"
  ,"LHC_LumAverageBRANB_4%c2_meanLuminosity"
  ,"LHC_LumAverageBRANB_4%c2_meanLuminosityError"
  ,"LHC_BeamLuminosityPerBunchBRANB_4%c2_Average_BunchLuminosity"
  ,"LHC_BeamLuminosityPerBunchBRANB_4%c2_BunchLuminosityError"
  ,"LHC_LumAverageBRANB_4%c2_meanCrossingAngle"
  ,"LHC_LumAverageBRANB_4%c2_meanCrossingAngleError"
  ,"LHC_CirculatingBunchConfig_Beam%d"
  ,"LHC_FillNumber"
  //
  ,"LHC_BunchLengthBeam%1d_nBunches"
  ,"LHC_BunchLengthBeam%1d_bunchesLenghts"
  ,"LHC_BunchLengthBeam%1d_filledBuckets"
  //
  ,"LHC_RunControl_ActiveInjectionScheme"
  ,"LHC_RunControl_BetaStar"
  ,"LHC_RunControl_IP2_Xing_Murad"
  ,"LHC_RunControl_IP2_ALICE_Murad"

  ,"LHC_BeamSizeBeam%1d_acqMode"
  ,"LHC_BeamSizeBeam%1d_sigmaH"
  ,"LHC_BeamSizeBeam%1d_sigmaV"
  ,"LHC_BeamSizeBeam%1d_emittanceH"
  ,"LHC_BeamSizeBeam%1d_emittanceV"
  ,"LHC_BeamSizeBeam%1d_errorSigmaH"
  ,"LHC_BeamSizeBeam%1d_errorSigmaV"
  //
  ,"LHC_CollimatorPos_%s_lvdt_%s"
  //
  ,"BPTX_deltaT_B1_B2"
  ,"BPTX_deltaTRMS_B1_B2"
  ,"BPTX_Phase_B%1d"
  ,"BPTX_PhaseRMS_B%1d"
  ,"BPTX_Phase_Shift_B%1d"
  //
  ,"ALI_Lumi_Total_Inst"
  ,"ALI_Lumi_Total_Delivered_StabBeam"
  ,"ALI_Lumi_Bunch_Inst"
  ,"ALI_Background%1d"
};

const Char_t* AliLHCData::fgkDCSColNames[] = {
  "TCTVB_4L2",
  "TCTVB_4R2",
  "TCLIA_4R2"
};

const Char_t* AliLHCData::fgkDCSColJaws[] = {
  "gap_downstream","gap_upstream","left_downstream",
  "left_upstream","right_downstream","right_upstream"};

//___________________________________________________________________
AliLHCData::AliLHCData(const TMap* dcsMap, double tmin, double tmax)
  : fTMin(0),fTMax(0),fFillNumber(0),fData(0),fkFile2Process(0),fkMap2Process(0)
{
  FillData(dcsMap,tmin,tmax);
}

//___________________________________________________________________
AliLHCData::AliLHCData(const Char_t* dcsFile, double tmin, double tmax)
  : fTMin(0),fTMax(0),fFillNumber(0),fData(0),fkFile2Process(dcsFile),fkMap2Process(0)
{
  FillData(dcsFile,tmin,tmax);
}

//___________________________________________________________________
Bool_t AliLHCData::FillData(const TMap*   dcsMap,double tmin, double tmax)
{
  // process DCS map and fill all fields. 
  Clear();
  fkMap2Process = dcsMap;
  FillData(tmin,tmax);
  return kTRUE;
}

//___________________________________________________________________
Bool_t AliLHCData::FillData(const Char_t* dcsFile,double tmin, double tmax)
{
  // process DCS file and fill all fields. 
  Clear();
  fkFile2Process = dcsFile;
  FillData(tmin,tmax);
  return kTRUE;
}

//___________________________________________________________________
Bool_t AliLHCData::FillData(double tmin, double tmax)
{
  // process DCS map and fill all fields. 
  // Accept only entries with timestamp between tmin and tmax
  //
  char buff[100],buff1[100];
  //
  SetTMin(tmin);
  SetTMax(tmax);
  //
  // -------------------------- extract Fill Number
  int iFirst=0,iLast=0;
  TObjArray* arr = GetDCSEntry(fgkDCSNames[kFillNum],iFirst,iLast,fTMin,fTMax);
  if (arr) SetFillNumber( ExtractInt( (AliDCSArray*)arr->At(iFirst), 0) );
  if (fkFile2Process) delete arr; // array was created on demand
  //
  for (int ibm=0;ibm<2;ibm++) {
    //
    snprintf(buff,99,fgkDCSNames[kBunchConf],ibm+1);         // ----- declared bunch configuration
    FillBunchConfig(fBunchConfDecl[ibm], buff);
    //
    snprintf(buff,99,fgkDCSNames[kBunchLgtFillB],ibm+1);     // ----- measured bunch configuration
    FillBunchConfig(fBunchConfMeas[ibm], buff);
    //
    snprintf(buff,99,fgkDCSNames[kBunchLgt],ibm+1);          // ----- measured bunch lenghts
    FillBunchInfo(fBunchLengths[ibm],buff,ibm,kFALSE);  
    //
    snprintf(buff,99,fgkDCSNames[kIntBunchAv],ibm+1);        // ----- B-by-B intensities
    FillBunchInfo(fIntensPerBunch[ibm],buff,ibm,kTRUE);
    //
    //
    snprintf(buff,99,fgkDCSNames[kIntTot],ibm+1);            // ----- total intensities for beam 1 and 2
    FillScalarRecord(fIntensTotal[ibm], buff);
    //
    snprintf(buff,99,fgkDCSNames[kIntTotAv],ibm+1);          // ----- total intensities for beam 1 and 2 from B-by-B average
    FillScalarRecord(fIntensTotalAv[ibm], buff);
    //
    snprintf(buff,99,fgkDCSNames[kBeamSzEmittH],ibm+1);      // ----- H emittance for beam 1 and 2 
    FillScalarRecord(fEmittanceH[ibm], buff);
    //
    snprintf(buff,99,fgkDCSNames[kBeamSzEmittV],ibm+1);      // ----- V emittance for beam 1 and 2 
    FillScalarRecord(fEmittanceV[ibm], buff);
    //
    snprintf(buff,99 ,fgkDCSNames[kBeamSzSigH],   ibm+1);    // ----- H sigmas and errors for beam 1 and 2 
    snprintf(buff1,99,fgkDCSNames[kBeamSzSigHErr],ibm+1);
    FillScalarRecord(fBeamSigmaH[ibm], buff, buff1);
    //
    snprintf(buff,99 ,fgkDCSNames[kBeamSzSigV],   ibm+1);    // ----- V sigmas and errors for beam 1 and 2 
    snprintf(buff1,99,fgkDCSNames[kBeamSzSigVErr],ibm+1);
    FillScalarRecord(fBeamSigmaV[ibm], buff, buff1);
    //
    //
    snprintf(buff,99,fgkDCSNames[kBPTXPhase],ibm+1);         // ----- BPTXPhase beam 1 and 2
    FillScalarRecord(fBPTXPhase[ibm], buff);
    //
    snprintf(buff,99,fgkDCSNames[kBPTXPhaseRMS],ibm+1);      // ----- BPTXPhaseRMS beam 1 and 2
    FillScalarRecord(fBPTXPhaseRMS[ibm], buff);
    //
    snprintf(buff,99,fgkDCSNames[kBPTXPhaseShift],ibm+1);    // ----- BPTXPhaseShift beam 1 and 2
    FillScalarRecord(fBPTXPhaseShift[ibm], buff);
    //
  }
  //
  for (int ibg=0;ibg<kNBGs;ibg++) {
    snprintf(buff,99,fgkDCSNames[kALIBackground],ibg+1);    // ----- Alice backgrounds 1,2,3
    FillScalarRecord(fBckgAlice[ibg], buff);
  }
  //
  FlagInteractingBunches(fBunchConfMeas[0],fBunchConfMeas[1]);
  FlagInteractingBunches(fBunchConfDecl[0],fBunchConfDecl[1]);
  //
  for (int ilr=0;ilr<2;ilr++) {
    //
    snprintf(buff,99 ,fgkDCSNames[kLumBunch], ilr ? 'R':'L');       // ---- BC-by-BC luminosity at IP2 and its error
    snprintf(buff1,99,fgkDCSNames[kLumBunchErr], ilr ? 'R':'L');
    FillBCLuminosities(fLuminPerBC[ilr], buff, buff1, 0); // BRAN L uses beam2 as a reference, BRAN R - beam1
    //
    snprintf(buff,99 ,fgkDCSNames[kLumTot]   , ilr ? 'R':'L');       // ---- total luminosity at IP2 and its error
    snprintf(buff1,99,fgkDCSNames[kLumTotErr], ilr ? 'R':'L');
    FillScalarRecord(fLuminTotal[ilr], buff, buff1);
    //
    snprintf(buff,99 ,fgkDCSNames[kLumAcqMode], ilr ? 'R':'L');      // ---- luminosity acquisition mode
    FillAcqMode(fLuminAcqMode[ilr], buff);
    //
    snprintf(buff,99, fgkDCSNames[kLumCrossAng]   , ilr ? 'R':'L');  //----- crossing angle at IP2 and its error
    snprintf(buff1,99,fgkDCSNames[kLumCrossAngErr], ilr ? 'R':'L');
    FillScalarRecord(fCrossAngle[ilr], buff, buff1);
    //    
  }
  //
  for (int icl=0;icl<kNCollimators;icl++) {             // ----- collimators positions
    for (int jaw=0;jaw<kNJaws;jaw++) {
      snprintf(buff,99,fgkDCSNames[kCollPos], fgkDCSColNames[icl],fgkDCSColJaws[jaw]);        
      FillScalarRecord(fCollimators[icl][jaw], buff);
    } // jaws
  } // collimators
  //
  //
  // RunControl info
  FillStringRecord(fRCInjScheme, fgkDCSNames[kRCInjSch]);   // ---- active injection scheme
  FillScalarRecord(fRCBeta, fgkDCSNames[kRCBeta]);          // ---- target beta 
  FillScalarRecord(fRCAngH, fgkDCSNames[kRCCrossAng]);      // ---- horisontal angle
  FillScalarRecord(fRCAngV,fgkDCSNames[kRCVang] );          // ---- vertical angle
  //
  FillScalarRecord(fBPTXdTB1B2,    fgkDCSNames[kBPTXdeltaTB1B2]);
  FillScalarRecord(fBPTXdTRMSB1B2, fgkDCSNames[kBPTXdeltaTRMSB1B2]);
  FillScalarRecord(fLumiAlice,     fgkDCSNames[kALILumiTotalInst]);
  FillScalarRecord(fLumiAliceStB,  fgkDCSNames[kALILumiTotalDeliveredStabBeam]);  
  FillBCLuminosities(fLumiAliceBbB,fgkDCSNames[kALILumiBunchInst],0,0);  
  //
  return kTRUE;
}

//___________________________________________________________________
TObjArray* AliLHCData::GetDCSEntry(const char* key,int &entry,int &last,double tmin,double tmax) const
{
  // extract array from the DCS map or file and find the first entry within the time limits
  entry = -1;
  last = -2;
  TObjArray* arr;
  if (fkMap2Process) arr = (TObjArray*)fkMap2Process->GetValue(key);
  else if (fkFile2Process) {
    AliLHCReader rd;
    arr = rd.ReadSingleLHCDP(fkFile2Process,key);
  }
  else {
    AliError("Neither DCS map nor DCS filename are set");
    return 0;  
  }
  //
  if (!arr || !arr->GetEntriesFast()) { 
    AliWarning(Form("No data for %s",key)); 
    if (fkMap2Process) delete arr; // created on demand
    return 0;
  }
  int ntot = arr->GetEntriesFast();
  //
  // search 1st entry before or at tmin
  AliDCSArray* ent = 0;
  Bool_t found = kFALSE;
  for (entry=0;entry<ntot;entry++) {
    ent = (AliDCSArray*)arr->At(entry);
    if (ent->GetTimeStamp()>=tmin-kMarginSOR && ent->GetTimeStamp()<=tmax+kMarginEOR) {
      found = kTRUE;
      if (ent->GetTimeStamp()>tmin) break;
    }
  }
  if (!found) {
    entry = -1;
    TString str;
    str += AliLHCDipValF::TimeAsString(tmin);
    str += " : ";
    str += AliLHCDipValF::TimeAsString(tmax);
    AliWarning(Form("All entries for %s are outside the requested range:\n%s",key,str.Data()));
    if (fkMap2Process) delete arr; // created on demand
    return 0;
  }
  if (entry>0) entry--;
  //
  // search last entry at or after tmin
  ent = 0;
  for (last=entry;last<ntot;last++) {
    ent = (AliDCSArray*)arr->At(last);
    if (ent->GetTimeStamp()>tmax) break;
  }
  if (last == ntot) last--;
  else if (ent->GetTimeStamp()>tmax+kMarginEOR) last--;
  //
  return arr;
}

//___________________________________________________________________
Int_t AliLHCData::TimeDifference(double v1,double v2,double tol) const
{
  // return 0 if the times are the same within the tolerance
  //        1 if v1>v2
  //       -1 if v1<v2
  v1-=v2;
  if (v1>tol)  return  1;
  if (v1<-tol) return -1;
  return 0;
}

//___________________________________________________________________
Bool_t AliLHCData::GoodPairID(int beam) const
{
  // check for correct beam identifier 
  if (beam>kBeam2||beam<0) {AliError(Form("BeamID can be 0 or 1, %d requested",beam)); return kFALSE;}
  return kTRUE;
}

//___________________________________________________________________
AliLHCDipValI* AliLHCData::GetBunchConfigMeasured(int beam,double tstamp) const
{
  // find measured bunch configuration valid for given tstamp
  if (!GoodPairID(beam)) return 0;
  return (AliLHCDipValI*)FindRecValidFor(fBunchConfMeas[beam][kStart],fBunchConfMeas[beam][kNStor],tstamp);
}

//___________________________________________________________________
AliLHCDipValI* AliLHCData::GetBunchConfigDeclared(int beam,double tstamp) const
{
  // find declared bunch configuration valid for given tstamp
  if (!GoodPairID(beam)) return 0;
  return (AliLHCDipValI*)FindRecValidFor(fBunchConfDecl[beam][kStart],fBunchConfDecl[beam][kNStor],tstamp);
}

//___________________________________________________________________
TObject* AliLHCData::FindRecValidFor(int start,int nrec, double tstamp) const
{
  // find record within this limits valid for given tstamp (i.e. the last one before or equal to tstamp)
  AliLHCDipValI *prevObj = 0;
  for (int i=0;i<nrec;i++) {
    AliLHCDipValI* curObj = (AliLHCDipValI*)fData[start+i];
    if (TimeDifference(tstamp,curObj->GetTimeStamp())<=0) break;
    prevObj = curObj;
  }
  if (!prevObj && nrec>0) prevObj = (AliLHCDipValI*)fData[start]; // if no exact match, return the 1st one
  return prevObj;
}

//___________________________________________________________________
Int_t AliLHCData::FillScalarRecord(int refs[2], const char* rec, const char* recErr, Double_t maxAbsVal)
{
  // fill record for scalar value, optionally accompanied by measurement error 
  //
  AliInfo(Form("Acquiring record: %s",rec));
  //
  TObjArray *arr=0,*arrE=0;
  Int_t iLast=0,iLastE=0,iFirst=0,iFirstE=0;
  //
  refs[kStart] = fData.GetEntriesFast();
  refs[kNStor] = 0;
  //
  if ( !(arr=GetDCSEntry(rec,iFirst,iLast,fTMin,fTMax)) ) return -1;
  //
  int dim = 1;
  if (recErr) {
    arrE = GetDCSEntry(recErr,iFirstE,iLastE,fTMin,fTMax);
    dim += 1;
  }
  //
  // Bool_t last = kFALSE;
  while (iFirst<=iLast) {
    AliDCSArray *dcsVal = (AliDCSArray*) arr->At(iFirst++);
    double tstamp = dcsVal->GetTimeStamp();
    //
    AliLHCDipValF* curValD = new AliLHCDipValF(dim,tstamp);  // start new period
    double vcheck = ExtractDouble(dcsVal,0);     // value
    if (TMath::Abs(vcheck) > maxAbsVal) {
      AliError(Form("ANOMALOUS VALUE %e for slot %d of %s: exceeds %e",vcheck, 0, rec, maxAbsVal));
      vcheck = 0.;
    }    
    (*curValD)[0] = vcheck;
    //
    if (recErr) {
      double errVal = -1;
      while (iFirstE<=iLastE) {       // try to find corresponding error
	AliDCSArray *dcsValE = (AliDCSArray*) arrE->At(iFirstE);
        double tstampE = dcsValE->GetTimeStamp();
        int tdif = TimeDifference(tstamp,tstampE);
        if (!tdif) { // error matches to value
	  errVal = ExtractDouble(dcsVal,0);     // value
	  if (TMath::Abs(errVal) > maxAbsVal) {
	    AliError(Form("ANOMALOUS VALUE %e for slot %d of %s: exceeds %e",errVal, 0, recErr, maxAbsVal));
	    errVal = 0.;
	  }
	  iFirstE++; 
	  break;
	}
        else if (tdif>0) iFirstE++; // error time lags behind, read the next one
        else break;                 // error time is ahead of value, no error associated
      }
      (*curValD)[dim-1] = errVal;   // error
      curValD->SetLastSpecial();    // lable the last entry as an error
    }
    //
    fData.Add(curValD);
    refs[kNStor]++;
    // if (last) break;
  }
  //
  if (fkFile2Process) {
    delete arr;
    delete arrE;
  }
  return refs[kNStor];
}

//___________________________________________________________________
Int_t AliLHCData::FillBunchConfig(int refs[2],const char* rec)
{
  // fill record for bunch configuration
  //
  AliInfo(Form("Acquiring record: %s",rec));
  TObjArray *arr;
  Int_t iLast,iFirst;
  //
  refs[kStart] = fData.GetEntriesFast();
  refs[kNStor] = 0;
  //
  if ( !(arr=GetDCSEntry(rec,iFirst,iLast,fTMin,fTMax)) ) return -1;
  //
  AliLHCDipValI* prevRecI=0;
  //  
  while (iFirst<=iLast) {
    AliDCSArray *dcsVal = (AliDCSArray*) arr->At(iFirst++);
    double tstamp = dcsVal->GetTimeStamp();
    //
    int bucket=0, nbunch=0, ndiff=0;
    int nSlots = dcsVal->GetNEntries();     // count number of actual bunches (non-zeros)
    int* dcsArr = dcsVal->GetInt();
    while(nbunch<nSlots && (bucket=dcsArr[nbunch])) {
      if (prevRecI && prevRecI->GetSize()>nbunch && bucket!=prevRecI->GetValue(nbunch)) ndiff++;
      nbunch++;
    }
    if (!nbunch) AliWarning(Form("%s record is present but empty: no beam?",rec));
    if (prevRecI && !ndiff && nbunch==prevRecI->GetSize()) continue; // record similar to previous one
    AliLHCDipValI* curValI = new AliLHCDipValI(nbunch,tstamp);      
    for (int i=nbunch;i--;) (*curValI)[i] = dcsArr[i];
    fData.Add(curValI);
    refs[kNStor]++;
    prevRecI = curValI;
  }
  //
  if (fkFile2Process) delete arr;
  return refs[kNStor];
}
 
//___________________________________________________________________
Int_t AliLHCData::FillAcqMode(int refs[2],const char* rec)
{
  // fill acquisition mode
  //
  AliInfo(Form("Acquiring record: %s",rec));
  TObjArray *arr;
  Int_t iLast,iFirst;
  //
  refs[kStart] = fData.GetEntriesFast();
  refs[kNStor] = 0;
  //
  if ( !(arr=GetDCSEntry(rec,iFirst,iLast,fTMin,fTMax)) ) return -1;
  //
  AliLHCDipValI* prevRecI=0;
  while (iFirst<=iLast) {
    AliDCSArray *dcsVal = (AliDCSArray*) arr->At(iFirst++);
    double tstamp = dcsVal->GetTimeStamp();
    //
    int nSlots = dcsVal->GetNEntries();
    if (nSlots<1) continue;
    int acqMode = dcsVal->GetInt()[0];
    if (prevRecI && (*prevRecI)[0] == acqMode) continue; // record similar to previous one
    AliLHCDipValI* curValI = new AliLHCDipValI(1,tstamp);      
    (*curValI)[0] = acqMode;
    fData.Add(curValI);
    refs[kNStor]++;
    prevRecI = curValI;
  }
  //
  if (fkFile2Process) delete arr;
  return refs[kNStor];
}
 
//___________________________________________________________________
Int_t AliLHCData::FillStringRecord(int refs[2],const char* rec)
{
  // fill record with string value
  //
  AliInfo(Form("Acquiring record: %s",rec));
  TString prevRec;
  TObjArray *arr;
  Int_t iLast,iFirst;
  //
  refs[kStart] = fData.GetEntriesFast();
  refs[kNStor] = 0;
  //
  if ( !(arr=GetDCSEntry(rec,iFirst,iLast,fTMin,fTMax)) ) return -1;
  //
  while (iFirst<=iLast) {
    AliDCSArray *dcsVal = (AliDCSArray*) arr->At(iFirst++);
    double tstamp = dcsVal->GetTimeStamp();
    //
    TString &str = ExtractString(dcsVal);
    if (!prevRec.IsNull()) {if (str == prevRec) continue;} // skip similar record
    else prevRec = str;
    //
    AliLHCDipValC* curValS = new AliLHCDipValC(1,tstamp);      
    curValS->SetValues(str.Data(),str.Length()+1);
    //
    fData.Add(curValS);
    refs[kNStor]++;
  }
  if (fkFile2Process) delete arr;
  return refs[kNStor];
}

//___________________________________________________________________
Int_t AliLHCData::FillBunchInfo(int refs[2],const char* rec, int ibm, Bool_t inRealSlots, Double_t maxAbsVal)
{
  // fill bunch properties for beam ibm
  // if inRealSlots = true, then the value is taken from bunchRFbucket/10, otherwise, the value 
  // for the i-th bunch is taken from the i-th element
  //
  AliInfo(Form("Acquiring record: %s",rec));
  TObjArray *arr;
  Int_t iLast,iFirst;
  //
  refs[kStart] = fData.GetEntriesFast();
  refs[kNStor] = 0;
  //
  if ( !(arr=GetDCSEntry(rec,iFirst,iLast,fTMin,fTMax)) ) return -1;
  //
  while (iFirst<=iLast) {
    AliDCSArray *dcsVal = (AliDCSArray*) arr->At(iFirst++);
    double tstamp = dcsVal->GetTimeStamp();
    //
    AliLHCDipValI *bconf = GetBunchConfigMeasured(ibm,tstamp);
    if (!bconf) {
      AliWarning(Form("Mearured bunch configuration for beam %d at t=%.1f is not available, trying declared one",ibm+1,tstamp));
      bconf = GetBunchConfigDeclared(ibm,tstamp);
    }
    if (!bconf) {
      AliWarning(Form("Declared bunch configuration for beam %d at t=%.1f is not available, skip this record",ibm+1,tstamp));
      return -1;
    }
    int nSlots = dcsVal->GetNEntries();     // count number of actual bunches (non-zeros)
    int nbunch = bconf->GetSize();
    if (nbunch>nSlots) {
      AliWarning(Form("More N bunches than slots in %s at time %.1f",rec,tstamp));
      continue;
    }
    double* dcsArr = dcsVal->GetDouble();
    AliLHCDipValF* curValD = new AliLHCDipValF(nbunch,tstamp);
    for (int i=nbunch;i--;) {
      int ind = inRealSlots ? (*bconf)[i]/10 : i;
      if (ind>nSlots) {
	AliError(Form("Bunch %d refers to wrong slot %d, set to -1",i,(*bconf)[i]));
	(*curValD)[i] = -1;
      }
      else {
	double vcheck = dcsArr[ind];
	if (TMath::Abs(vcheck) > maxAbsVal) {
	  AliError(Form("ANOMALOUS VALUE %e for slot %d of %s: exceeds %e",vcheck, ind, rec, maxAbsVal));
	  vcheck = 0.;
	}
	(*curValD)[i] = vcheck;
      }
    }
    fData.Add(curValD);
    refs[kNStor]++;
  }
  if (fkFile2Process) delete arr;
  return refs[kNStor];
  //
}
 
//___________________________________________________________________
Int_t AliLHCData::FillBCLuminosities(int refs[2],const char* rec, const char* recErr, int useBeam, Double_t maxAbsVal)
{
  // fill luminosities per bunch crossing
  //
  AliInfo(Form("Acquiring record: %s",rec));
  TObjArray *arr,*arrE=0;
  Int_t iLast=0,iLastE=0,iFirst=0,iFirstE=0;
  //
  refs[kStart] = fData.GetEntriesFast();
  refs[kNStor] = 0;
  //
  if ( !(arr=GetDCSEntry(rec,iFirst,iLast,fTMin,fTMax)) ) return -1;
  //
  while (iFirst<=iLast) {
    AliDCSArray *dcsVal = (AliDCSArray*) arr->At(iFirst++);
    double tstamp = dcsVal->GetTimeStamp();
    //
    AliLHCDipValI *bconf;
    bconf = GetBunchConfigMeasured(useBeam,tstamp);  // luminosities are stored according to beam bunches
    if (!bconf) {
      AliWarning(Form("Mearured bunch configuration for beam%d at t=%.1f is not available, trying declared one",useBeam,tstamp));
      bconf = GetBunchConfigDeclared(useBeam,tstamp);
    }
    if (!bconf) {
      AliWarning(Form("Declared bunch configuration for beam%i at t=%.1f is not available, skip this record",useBeam,tstamp));
      return -1;
    }
    int nSlots = dcsVal->GetNEntries();     // count number of actual bunches (non-zeros)
    int nbunch = bconf->GetSize();
    double* dcsArr = dcsVal->GetDouble();
    //
    if (nbunch>nSlots) {
      AliWarning(Form("More N bunches than slots in %s at time %.1f",rec,tstamp));
      continue;
    }
    int dim = 0;
    if (!bconf->IsProcessed1()) {
      AliWarningF("Bunch conf. for beam%d has no marked interacting bunches, store all luminosity for all %d filled bunches",useBeam,nbunch);
      dim = nbunch;
    }
    else { // count number of interacting bunches
      for (int i=nbunch;i--;) if ((*bconf)[i]<0) dim++;
    }
    //
    if (recErr) {
      arrE=GetDCSEntry(recErr,iFirstE,iLastE,fTMin,fTMax);
      dim += 1;
    }
    AliLHCDipValF* curValD = new AliLHCDipValF(dim,tstamp);
    int cnt = 0;
    for (int i=0;i<nbunch;i++) {
      int slot = (*bconf)[i];
      if (bconf->IsProcessed1() && slot>0) continue;
      //
      int ind = TMath::Abs(slot)/10;
      if (ind>nSlots) {
	AliError(Form("Bunch %d refers to wrong slot %d, set to -1",cnt,slot));
	(*curValD)[cnt] = -1;
      }
      else {
	double vcheck = dcsArr[ind];
	if (TMath::Abs(vcheck) > maxAbsVal) {
	  AliError(Form("ANOMALOUS VALUE %e for slot %d of %s: exceeds %e",vcheck, ind, rec, maxAbsVal));
	  vcheck = 0.;
	}
	(*curValD)[cnt] = vcheck;
      }
      cnt++;
    }
    //
    if (recErr) {
      double errVal = -1;
      while (iFirstE<=iLastE) {       // try to find corresponding error
	AliDCSArray *dcsValE = (AliDCSArray*) arrE->At(iFirstE);
	double tstamp1 = dcsValE->GetTimeStamp();
	int tdif = TimeDifference(tstamp,tstamp1);
	if (!tdif) { // error matches to value
	  errVal = dcsValE->GetDouble()[0];
	  if (TMath::Abs(errVal) > maxAbsVal) {
	    AliError(Form("ANOMALOUS VALUE %e for slot %d of %s: exceeds %e",errVal,0, rec, maxAbsVal));
	    errVal = 0.;
	  }	  
	  iFirstE++; 
	  break;
	}
	else if (tdif>0) iFirstE++; // error time lags behind, read the next one
	else break;                 // error time is ahead of value, no error associated
      }
      (*curValD)[dim-1] = errVal;   // error
      curValD->SetLastSpecial();    // lable the last entry as an error
    }
    //
    fData.Add(curValD);
    refs[kNStor]++;
  }
  if (fkFile2Process) {
    delete arr;
    delete arrE;
  }
  return refs[kNStor];
  //
}

//___________________________________________________________________
Int_t AliLHCData::ExtractInt(AliDCSArray* dcsArray,Int_t el) const
{
  // extract integer from the dcsArray
  int val = -1;
  //
  int sz = dcsArray->GetNEntries();
  if (sz<=el) return val;
  //
  if (dcsArray->GetType()==AliDCSArray::kInt)  val = dcsArray->GetInt(el);
  else if (dcsArray->GetType()==AliDCSArray::kString) {
    TObjString *stro = dcsArray->GetStringArray(el);
    if (stro) val = stro->GetString().Atoi();
    else AliError(Form("DCSArray TObjString for element %d is missing",el));
  }
  else if (dcsArray->GetType()==AliDCSArray::kUInt) val = dcsArray->GetUInt(el);
  else AliError(Form("Integer requested from DCSArray of type %d",dcsArray->GetType()));
  return val;
}

//___________________________________________________________________
Double_t AliLHCData::ExtractDouble(AliDCSArray* dcsArray,Int_t el) const
{
  // extract double from the dcsArray
  double val = 0;
  //
  int sz = dcsArray->GetNEntries();
  if (sz<=el) return val;
  //
  if      (dcsArray->GetType()==AliDCSArray::kDouble) val = dcsArray->GetDouble(el);
  else if (dcsArray->GetType()==AliDCSArray::kFloat)  val = dcsArray->GetFloat(el);
  else if (dcsArray->GetType()==AliDCSArray::kString) {
    TObjString *stro = dcsArray->GetStringArray(el);
    if (stro) val = stro->GetString().Atof();
    else AliError(Form("DCSArray has TObjString for element %d is missing",el));
  }
  else if (dcsArray->GetType()==AliDCSArray::kChar)   val = dcsArray->GetChar(el);
  else if (dcsArray->GetType()==AliDCSArray::kInt)    val = dcsArray->GetInt(el);
  else if (dcsArray->GetType()==AliDCSArray::kUInt)   val = dcsArray->GetUInt(el);
  else     AliError(Form("Double requested from DCSArray of type %d",dcsArray->GetType()));
  return val;
}

//___________________________________________________________________
TString& AliLHCData::ExtractString(AliDCSArray* dcsArray) const
{
  // extract string from the dcsArray
  static TString str;
  str = "";
  //
  int sz = dcsArray->GetNEntries();
  if (dcsArray->GetType()!=AliDCSArray::kString)  {
    AliError(Form("String requested from DCSArray of type %d",dcsArray->GetType()));
    return str;
  }
  //
  for (int i=0;i<sz;i++) {
    str += dcsArray->GetStringArray(i)->GetString();
    if (i<sz-1) str += " ";
  }
  return str;
}

//___________________________________________________________________
void AliLHCData::Print(const Option_t* opt) const
{
  // print full info
  TString opts = opt;
  opts.ToLower();
  Bool_t utcTime = opts.Contains("loc") ? kFALSE:kTRUE;
  //
  Bool_t includeMissing = opts.Contains("m");
  if (!includeMissing) printf("Missing records are skept, use Print(\"m\") to print missing record names\n");
  else                 printf("Missing records are printed, remove \"m\" from Print options to skip them\n");
  Bool_t full = kTRUE;
  if (!opts.Contains("f")) {
    printf("Use Print(\"f\") to print full info\n");
    printf("Printing short info:\n<RecordType>(number of records): <TimeStamp, value> for 1st record only\n");
    full = kFALSE;
  }
  TString sdtmn = AliLHCDipValI::TimeAsString(fTMin,utcTime);
  TString sdtmx = AliLHCDipValI::TimeAsString(fTMax,utcTime);
  printf("Fill#%6d Validity: %s - %s (%s)\n",fFillNumber,sdtmn.Data(),sdtmx.Data(),utcTime ? "UTC":"LOC");
  //
  printf("********** SETTINGS FROM RUN CONTROL **********\n");
  //
  if (fRCInjScheme[kNStor] || includeMissing) {
    printf("* %-38s","Injection Scheme");
    PrintAux(full,fRCInjScheme,opts);
  }
  //
  if (fRCBeta[kNStor] || includeMissing) {
    printf("* %-38s","Beta Star");
    PrintAux(full,fRCBeta,opts);
  }
  //
  if (fRCAngH[kNStor] || includeMissing) {
    printf("* %-38s","Horisontal Crossing Angle");
    PrintAux(full,fRCAngH,opts);
  }
  //
  if (fRCAngV[kNStor] || includeMissing) {
    printf("* %-38s","Vertical   Crossing Angle");
    PrintAux(full,fRCAngV,opts);
  }
  //
  for (int ib=0;ib<2;ib++) {
    if (fBunchConfDecl[ib][kNStor] || includeMissing) {
      printf("* Beam%d filling  [- interacts at IR2!]  ",ib+1);
      PrintAux(full,fBunchConfDecl[ib],opts);
    }
  }
  //
  printf("\n**********       MEASURED DATA       **********\n");
  //
  for (int ib=0;ib<2;ib++) {
    if (fBunchConfMeas[ib][kNStor] || includeMissing) {
      printf("* Beam%d filling  [- interacts at IR2!]  ",ib+1);
      PrintAux(full,fBunchConfMeas[ib],opts);
    }
  } 
  //
  for (int ib=0;ib<2;ib++) {
    if (fIntensTotal[ib][kNStor] || includeMissing) {
      printf("* Beam%d total intensity                 ",ib+1);
      PrintAux(full,fIntensTotal[ib],opts);
    }
  } 
  //
  for (int ib=0;ib<2;ib++) {
    if (fIntensTotalAv[ib][kNStor] || includeMissing) {
      printf("* Beam%d total intensity (bunch average) ",ib+1);
      PrintAux(full,fIntensTotalAv[ib],opts);
    }
  } 
  //
  for (int ib=0;ib<2;ib++) {
    if (fIntensPerBunch[ib][kNStor] || includeMissing) {
      printf("* Beam%d intensity per bunch             ",ib+1);
      PrintAux(full,fIntensPerBunch[ib],opts);
    }
  }
  //
  for (int ib=0;ib<2;ib++) {
    if (fBunchLengths[ib][kNStor] || includeMissing) {
      printf("* Beam%d bunch lengths                   ",ib+1);
      PrintAux(full,fBunchLengths[ib],opts);
    }
  } 
  //
  for (int ib=0;ib<2;ib++) {
    if (fEmittanceH[ib][kNStor] || includeMissing) {
      printf("* Beam%d Horisontal emittance            ",ib+1);
      PrintAux(full,fEmittanceH[ib],opts);
    }
  }
  //
  for (int ib=0;ib<2;ib++) {
    if (fEmittanceV[ib][kNStor] || includeMissing) {
      printf("* Beam%d Vertical emittance              ",ib+1);
      PrintAux(full,fEmittanceV[ib],opts);
    }
  }
  //
  for (int ib=0;ib<2;ib++) {
    if (fBeamSigmaH[ib][kNStor] || includeMissing) {
      printf("* Beam%d Horisontal sigma                ",ib+1);
      PrintAux(full,fBeamSigmaH[ib],opts);
    }
  }
  //
  for (int ib=0;ib<2;ib++) {
    if (fBeamSigmaV[ib][kNStor] || includeMissing) {
      printf("* Beam%d Vertical sigma                  ",ib+1);
      PrintAux(full,fBeamSigmaV[ib],opts);
    }
  }
  //
  for (int lr=0;lr<2;lr++) {
    if (fLuminTotal[lr][kNStor] || includeMissing) {
      printf("* Total luminosity from BRANB_4%c2       ",lr ? 'R':'L');
      PrintAux(full,fLuminTotal[lr],opts);
    }
  } 
  //
  for (int lr=0;lr<2;lr++) {
    if (fLuminAcqMode[lr][kNStor] || includeMissing) {
      printf("* Luminosity acq.mode, BRANB_4%c2        ",lr ? 'R':'L');
      PrintAux(full,fLuminAcqMode[lr],opts+"bit");
    }
  } 
  //
  for (int lr=0;lr<2;lr++) {
    if (fLuminPerBC[lr][kNStor] || includeMissing) {
      printf("* Luminosity per BC from BRANB_4%c2      ",lr ? 'R':'L');
      PrintAux(full,fLuminPerBC[lr],opts);
    }
  }
  //
  for (int lr=0;lr<2;lr++) {
    if (fCrossAngle[lr][kNStor] || includeMissing) {
      printf("* Crossing angle, side %c                ",lr ? 'R':'L');
      PrintAux(full,fCrossAngle[lr],opts);
    }
  }
  //
  for (int coll=0;coll<kNCollimators;coll++)
    for (int jaw=0;jaw<kNJaws;jaw++) {
      if (fCollimators[coll][jaw][kNStor] || includeMissing) {
	printf("* Collimator %10s:%16s",fgkDCSColNames[coll],fgkDCSColJaws[jaw]);
	PrintAux(full,fCollimators[coll][jaw],opts);
      }
    }
  //
  printf("\n**********    ALICE MEASURED DATA    **********\n");
  //
  if (fLumiAlice[kNStor] || includeMissing) {
    printf("* %-38s","Alice luminosity total");
    PrintAux(full,fLumiAlice,opts);
  }
  //
  if (fLumiAliceStB[kNStor] || includeMissing) {
    printf("* %-38s","Alice luminosity delivered stable beam");
    PrintAux(full,fLumiAliceStB,opts);
  }
  //
  if (fLumiAliceBbB[kNStor] || includeMissing) {
    printf("* %-38s","Alice luminosity B-by-B, stable beam");
    PrintAux(full,fLumiAliceBbB,opts);
  }
  //
  for (int ib=0;ib<3;ib++) {
    if (fBckgAlice[ib][kNStor] || includeMissing) {
      printf("* Alice background%d                     ",ib+1);
      PrintAux(full,fBckgAlice[ib],opts);
    }
  }
  //
  if (fBPTXdTB1B2[kNStor] || includeMissing) {
    printf("* %-38s","BPTX DeltaT Beam1 Beam2");
    PrintAux(full,fBPTXdTB1B2,opts);
  }
  //
  if (fBPTXdTRMSB1B2[kNStor] || includeMissing) {
    printf("* %-38s","BPTX DeltaT RMS Beam1 Beam2");
    PrintAux(full,fBPTXdTRMSB1B2,opts);
  }
  //
  for (int ib=0;ib<2;ib++) {
    if (fBPTXPhase[ib][kNStor] || includeMissing) {
      printf("* BPTX Phase Beam%d                      ",ib+1);
      PrintAux(full,fBPTXPhase[ib],opts);
    }
  }
  //
  for (int ib=0;ib<2;ib++) {
    if (fBPTXPhaseRMS[ib][kNStor] || includeMissing) {
      printf("* BPTX Phase RMS Beam%d                  ",ib+1);
      PrintAux(full,fBPTXPhaseRMS[ib],opts);
    }
  }
  //
  for (int ib=0;ib<2;ib++) {
    if (fBPTXPhaseShift[ib][kNStor] || includeMissing) {
      printf("* BPTX Phase Shift Beam%d                ",ib+1);
      PrintAux(full,fBPTXPhaseShift[ib],opts);
    }
  }
  //
}

//___________________________________________________________________
void AliLHCData::PrintAux(Bool_t full, const Int_t refs[2], const Option_t *opt) const
{
  // aux method to print the reocrds of the same type
  int nrec = refs[kNStor];
  if (nrec<1) {
    printf(": N/A\n"); 
    return;
  }
  printf(": (%3d):\t",nrec); // number of records
  if (!full) nrec = 1;
  int sz = ((AliLHCDipValI*)fData[refs[kStart]])->GetSizeTotal(); // dimension of the record
  Bool_t isStr = ((AliLHCDipValI*)fData[refs[kStart]])->IsTypeC();
  if ((!isStr && sz>2) || nrec>1) printf("\n"); // long record, open new line
  for (int i=0;i<nrec;i++) fData[refs[kStart]+i]->Print(opt);
  //
}

//___________________________________________________________________
void AliLHCData::Clear(const Option_t *)
{
  // clear all info
  fData.Delete();
  fFillNumber = 0;
  fTMin = 0;
  fTMax = 1e10;
  fkFile2Process = 0;
  fkMap2Process  = 0;
  //
  for (int i=2;i--;) {
    fRCInjScheme[i] = 0;
    fRCBeta[i] = 0;
    fRCAngH[i] = 0;
    fRCAngV[i] = 0;
    fLumiAlice[i] = 0;
    fLumiAliceStB[i] = 0;
    fLumiAliceBbB[i] = 0;
    for (int ibg=kNBGs;ibg--;) fBckgAlice[ibg][i] = 0;
    fBPTXdTB1B2[i] = 0;
    fBPTXdTRMSB1B2[i] = 0;
    //
    for (int icl=kNCollimators;icl--;) for (int jaw=kNJaws;jaw--;) fCollimators[icl][jaw][i]=0;
    //
    for (int j=2;j--;) {
      fBPTXPhase[j][i] = 0;
      fBPTXPhaseRMS[j][i] = 0;
      fBPTXPhaseShift[j][i] = 0;
      //
      fBunchConfDecl[j][i] = 0;
      fBunchConfMeas[j][i] = 0;
      fBunchLengths[j][i] = 0;
      fIntensTotal[j][i] = 0;
      fIntensTotalAv[j][i] = 0;
      fIntensPerBunch[j][i] = 0;      
      fCrossAngle[j][i] = 0;
      fEmittanceH[j][i] = 0;
      fEmittanceV[j][i] = 0;
      fBeamSigmaH[j][i] = 0;
      fBeamSigmaV[j][i] = 0;
      fLuminTotal[j][i] = 0;
      fLuminPerBC[j][i] = 0;
      fLuminAcqMode[j][i] = 0;
    }
  }
}

//___________________________________________________________________
Int_t AliLHCData::GetNInteractingBunchesMeasured(int i) const
{
  // get number of interacting bunches at IR2
  AliLHCDipValI* rec = GetBunchConfigMeasured(kBeam1,i);
  if (!rec) {AliInfo(Form("No record %d found",i)); return -1;}
  if (!rec->IsProcessed1()) { AliInfo("Interacting bunches were not marked"); return -1;}
  int n = 0;
  for (int j=rec->GetSize();j--;) if ( (*rec)[j]<0 ) n++;
  return n;
}

//___________________________________________________________________
Int_t AliLHCData::GetNInteractingBunchesDeclared(int i) const
{
  // get number of interacting bunches at IR2
  AliLHCDipValI* rec = GetBunchConfigMeasured(kBeam1,i);
  if (!rec) {AliInfo(Form("No record %d found",i)); return -1;}
  if (!rec->IsProcessed1()) { AliInfo("Interacting bunches were not marked"); return -1; }
  int n = 0;
  for (int j=rec->GetSize();j--;) if ( (*rec)[j]<0 ) n++;
  return n;
}

//___________________________________________________________________
Int_t AliLHCData::IsPilotPresent(int i) const
{
  // check in the filling scheme is the pilot bunch is present
  AliLHCDipValC* rec = GetInjectionScheme();
  if (!rec) {AliInfo(Form("No record %d found",i)); return -1;}
  TString scheme = rec->GetValues();
  return scheme.Contains("wp",TString::kIgnoreCase);
}

//___________________________________________________________________
void AliLHCData::FlagInteractingBunches(const Int_t beam1[2],const Int_t beam2[2])
{
  // assign - sign to interacting bunches
  //
  for (int ib1=0;ib1<beam1[kNStor];ib1++) {
    AliLHCDipValI *bm1 = (AliLHCDipValI*)fData[ beam1[kStart] + ib1];
    if (!bm1) continue;
    AliLHCDipValI *bm2 = (AliLHCDipValI*)FindRecValidFor(beam2[kStart],beam2[kNStor], bm1->GetTimeStamp());
    if (!bm2) continue;
    //
    int nb1 = bm1->GetSize();
    int nb2 = bm2->GetSize();
    int i1,i2;
    for (i1=0;i1<nb1;i1++) {
      int bunch2=-1, bunch1 = TMath::Abs((*bm1)[i1]);
      int slot2 =-1, slot1  = GetBCId(bunch1,0);
      for (i2=0;i2<nb2;i2++) {
	bunch2 = TMath::Abs((*bm2)[i2]);
	slot2 = GetBCId(bunch2,1);
	if (slot1==slot2) break;
      }
      if (slot1!=slot2) continue;
      (*bm1)[i1] = -bunch1;
      (*bm2)[i2] = -bunch2;
      bm1->SetProcessed1();
      bm2->SetProcessed1();
    }
  }
}

//___________________________________________________________________
Int_t AliLHCData::GetMeanIntensity(int beamID, Double_t &colliding, Double_t &noncolliding, const TObjArray* bcmasks) const
{
  // get average intensity for all, colliding and non-colliding bunches
  // on success returns number of intensity records used (1 per ~10 min)
  // If triggered BC masks are provided, calculation is done for Triggered BC only
  colliding = noncolliding = -1.;
  if (beamID<0||beamID>1) {
    AliError(Form("BeamID must be either 0 or 1, %d requested",beamID));
    return -10;
  }
  //
  AliTriggerBCMask *bcMaskBoth=0,*bcMaskSingle=0;
  int nbcm = 0;
  if (bcmasks && (nbcm=bcmasks->GetEntries())) {
    if (nbcm>1) bcMaskBoth = (AliTriggerBCMask*)bcmasks->At(1);
    if      (nbcm>0 && beamID==kBeam1) bcMaskSingle = (AliTriggerBCMask*)bcmasks->At(0);
    else if (nbcm>2 && beamID==kBeam2) bcMaskSingle = (AliTriggerBCMask*)bcmasks->At(2);
    //
    if (!bcMaskSingle) AliError(Form("Only triggered BSs are requested but %c mask is not provided",beamID ? 'C':'A'));
    if (!bcMaskBoth)   AliError("Only triggered BSs are requested but B mask is not provided");
  }
  else {
    AliWarning("No BC masks are provided");
  }
  //
  int nrec = GetNIntensityPerBunch(beamID);
  if (nrec<1) return -1;
  AliLHCDipValI *conf = GetBunchConfigMeasured(beamID);
  if (!conf) conf = GetBunchConfigDeclared(beamID);
  if (!conf) return -2;
  int nb = conf->GetSize();
  //
  for (int irec=0;irec<nrec;irec++) {
    //
    AliLHCDipValD* rIntD = 0;
    AliLHCDipValF* rIntF = GetIntensityPerBunch(beamID,irec);
    // for BWD compatibility of some periods
    if (rIntF->IsA() == AliLHCDipValD::Class()) {rIntD=(AliLHCDipValD*)rIntF; rIntF=0;}
    if (!rIntF && !rIntD) {
      AliError(Form("Failed to get GetIntensityPerBunch(%d,%d)",beamID,irec));
      continue;
    }
    for (int ib=0;ib<nb;ib++) {
      double val = 0;
      if (rIntF) val = rIntF->GetValue(ib);
      else if (rIntD) val = rIntD->GetValue(ib);
      if (val<0) continue;
      int bID = conf->GetValue(ib);
      // check if this is a triggered bunch
      int bcID = GetBCId(bID, beamID);
      if (bID<0) { // interacting
	if (bcMaskBoth && bcMaskBoth->GetMask(bcID)) continue; // masked
	colliding += val;
      }
      else {
	if (bcMaskSingle && bcMaskSingle->GetMask(bcID)) continue; // masked	
	noncolliding += val;
      }
    }
  }
  colliding /= nrec;
  noncolliding /= nrec;
  return nrec;
}

/***************************************************************************************
// this is for the obsolete retrofitting

//___________________________________________________________________
void AliLHCData::FillLumiAliceOFL(Int_t nrec, Int_t* timeArr, Double_t* valArr)
{
  // Create a record for lumi integrated from the beginning of fill
  // We need dedicated method since this info comes from Alice (not LHCDip) as instantaneous values
  // and is retrofitted for past runs
  fLumiAlice[kStart] = fData.GetEntriesFast();
  fLumiAlice[kNStor] = 0;
  if (nrec<2 || !timeArr || !valArr) return;
  double tprv,period,currTime;
  if ((currTime=double(UInt_t(timeArr[0])))>fTMin) {
    AliError(Form("TimeStamp of 1st record: %s > TimeStamp of SOR: %s, STOP",
		  AliLHCDipValI::TimeAsString(currTime),AliLHCDipValI::TimeAsString(fTMin)));
    return;
  }
  //
  if ((tprv=double(UInt_t(timeArr[nrec-1])))<fTMax) {
    AliWarning(Form("TimeStamp of last (%d) record: %s < TimeStamp of EOR: %s, Data will be truncated",nrec-1,
		    AliLHCDipValI::TimeAsString(tprv),AliLHCDipValI::TimeAsString(fTMax)));
  }
  //
  // init the average time step
  period = (tprv - currTime)/(nrec-1);
  double lumiInt  = 0;
  //
  for (int i=0;i<nrec;i++) {
    tprv = currTime;
    currTime = double(UInt_t(timeArr[i]));
    if (i>0) period = currTime - tprv;
    if (currTime-period>fTMax) continue;    
    lumiInt += valArr[i]*period;    
    //    printf("%d %.2f V:%f Int:%f\n",i,period,valArr[i],lumiInt);
    if (currTime+period<fTMin) continue;
    AliLHCDipValF* curValF = new AliLHCDipValF(1,currTime);
    (*curValF)[0] = lumiInt;
    fData.Add(curValF);
    fLumiAlice[kNStor]++;
  }
  //
  printf("Stored %d Alice Integrated luminosity records out of %d provided\n",fLumiAlice[kNStor],nrec);
}

//___________________________________________________________________
void AliLHCData::FillBckgAliceOFL(Int_t nrec, Int_t* timeArr, Double_t* valArr)
{
  // Create a record for lumi integrated from the beginning of fill
  // We need dedicated method since this info comes from Alice (not LHCDip) as instantaneous values
  // and is retrofitted for past runs
  fBckgAlice[kStart] = fData.GetEntriesFast();
  fBckgAlice[kNStor] = 0;
  if (nrec<2 || !timeArr || !valArr) return;
  double tprv,period,currTime;
  if ((currTime=double(UInt_t(timeArr[0])))>fTMin) {
    AliError(Form("TimeStamp of 1st record: %s > TimeStamp of SOR: %s, STOP",
		  AliLHCDipValI::TimeAsString(currTime),AliLHCDipValI::TimeAsString(fTMin)));
    return;
  }
  //
  if ((tprv=double(UInt_t(timeArr[nrec-1])))<fTMax) {
    AliWarning(Form("TimeStamp of last (%d) record: %s < TimeStamp of EOR: %s, Data will be truncated",nrec-1,
		    AliLHCDipValI::TimeAsString(tprv),AliLHCDipValI::TimeAsString(fTMax)));
  }
  //
  // init the average time step
  period = (tprv - currTime)/(nrec-1);
  double bckgInt  = 0;
  //
  for (int i=0;i<nrec;i++) {
    tprv = currTime;
    currTime = double(UInt_t(timeArr[i]));
    if (i>0) period = currTime - tprv;
    if (currTime-period>fTMax) continue;    
    bckgInt += valArr[i]*period;    
    if (currTime+period<fTMin) continue;
    AliLHCDipValF* curValF = new AliLHCDipValF(1,currTime);
    (*curValF)[0] = bckgInt;
    fData.Add(curValF);
    fBckgAlice[kNStor]++;
  }
  //
  printf("Stored %d Alice Integrated Background records out of %d provided\n",fBckgAlice[kNStor],nrec);
}

***************************************************************************************/

/*
//_____________________________________________________________________________
Float_t AliLHCData::GetLumiInstAlice(Double_t tStamp) const 
{
  // get closest in time value on inst luminosity
  int idx = FindEntryValidFor(fLumiAlice[kStart],fLumiAlice[kNStor],tStamp);
  if (idx<0) return -1;
  AliLHCDipValF *rec=GetLumiAlice(idx),*rec1=GetLumiAlice(idx>0 ? idx-1:idx+1);
  if (!rec || !rec1) return -1;
  double dt = rec->GetTimeStamp() - rec1->GetTimeStamp();
  return TMath::Abs(dt)>1e-6 ? (rec->GetValue()-rec1->GetValue())/dt : -1;
}

//_____________________________________________________________________________
Float_t AliLHCData::GetBckgInstAlice(Double_t tStamp) const 
{
  // get closest in time value on inst luminosity
  int idx = FindEntryValidFor(fBckgAlice[kStart],fBckgAlice[kNStor],tStamp);
  if (idx<0) return -1;
  AliLHCDipValF *rec=GetBckgAliceRecord(idx),*rec1=GetBckgAliceRecord(idx>0 ? idx-1:idx+1);
  if (!rec || !rec1) return -1;
  double dt = rec->GetTimeStamp() - rec1->GetTimeStamp();
  return TMath::Abs(dt)>1e-6 ? (rec->GetValue()-rec1->GetValue())/dt : -1;
}
*/

//_____________________________________________________________________________
TGraph* AliLHCData::ExportGraph(Int_t *coord, Int_t elID) const
{
  // export time/values to graph:
  // coord: int[2] array with 1st entry and number of entries stored, obtained via GetOffs... method 
  // elID - element of the AliLHCDipValT array to extract
  if (!coord || coord[1]<1) return 0;
  TGraph* gr = new TGraph(coord[1]);
  for (int i=0;i<coord[1];i++) {
    TObject* obj = fData.At(coord[0]+i);
    if (!obj) {
      AliError(Form("Entry %d does not exist",i));
      continue;
    }
    if (obj->IsA()==AliLHCDipValD::Class()) {
      AliLHCDipValD* objD =  (AliLHCDipValD*)obj;
      gr->SetPoint(i,objD->GetTimeStamp(),objD->GetValue(elID));
    }
    else if (obj->IsA()==AliLHCDipValF::Class()) {
      AliLHCDipValF* objF =  (AliLHCDipValF*)obj;
      gr->SetPoint(i,objF->GetTimeStamp(),objF->GetValue(elID));
    }
    else if (obj->IsA()==AliLHCDipValI::Class()) {
      AliLHCDipValI* objI =  (AliLHCDipValI*)obj;
      gr->SetPoint(i,objI->GetTimeStamp(),objI->GetValue(elID));
    }
    else if (obj->IsA()==AliLHCDipValC::Class()) {
      AliLHCDipValC* objC =  (AliLHCDipValC*)obj;
      gr->SetPoint(i,objC->GetTimeStamp(),objC->GetValue(elID));
    }
    else {
      AliError(Form("Graph cannot be exported for records of type %s",obj->IsA()->GetName()));
    }
  }
  return gr;
}
