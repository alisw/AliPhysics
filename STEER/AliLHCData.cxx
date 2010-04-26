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
#include <TString.h>
#include <TObjArray.h>

ClassImp(AliLHCData)

const Char_t* AliLHCData::fgkDCSNames[] = {
  "LHC_IntensityBeam%d_totalIntensity",
  "LHC_BeamIntensityPerBunchBeam%d_averageBeamIntensity",
  "LHC_BeamIntensityPerBunchBeam%d_Average_BunchIntensities",
  //
  "LHC_LumAverageBRANB_4%c2_acqMode",
  "LHC_LumAverageBRANB_4%c2_meanLuminosity",
  "LHC_LumAverageBRANB_4%c2_meanLuminosityError",
  "LHC_BeamLuminosityPerBunchBRANB_4%c2_Average_BunchLuminosity",
  "LHC_BeamLuminosityPerBunchBRANB_4%c2_BunchLuminosityError",
  "LHC_LumAverageBRANB_4%c2_meanCrossingAngle",
  "LHC_LumAverageBRANB_4%c2_meanCrossingAngleError",
  "LHC_CirculatingBunchConfig_Beam%d",
  "LHC_FillNumber",
  //
  "LHC_BunchLengthBeam%d_nBunches",
  "LHC_BunchLengthBeam%d_bunchesLenghts",
  "LHC_BunchLengthBeam%d_filledBuckets",
  //
  "LHC_RunControl_ActiveInjectionScheme",
  "LHC_RunControl_BetaStar",
  "LHC_RunControl_IP2_Xing_Murad",
  "LHC_RunControl_IP2_ALICE_Murad",

  "LHC_BeamSizeBeam%d_acqMode",
  "LHC_BeamSizeBeam%d_sigmaH",
  "LHC_BeamSizeBeam%d_sigmaV",
  "LHC_BeamSizeBeam%d_emittanceH",
  "LHC_BeamSizeBeam%d_emittanceV",
  "LHC_BeamSizeBeam%d_errorSigmaH",
  "LHC_BeamSizeBeam%d_errorSigmaV",
  //
  "LHC_CollimatorPos_%s_lvdt_%s"
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
AliLHCData::AliLHCData(const TMap* dcsMap, double tmin, double tmax, Bool_t ignoreSOR,Bool_t ignoreEOR)
  : fTMin(0),fTMax(0),fFillNumber(0),fData(0),fFile2Process(0),fMap2Process(0)
{
  SetIgnoreSOR(ignoreSOR);
  SetIgnoreEOR(ignoreEOR);
  FillData(dcsMap,tmin,tmax);
}

//___________________________________________________________________
AliLHCData::AliLHCData(const Char_t* dcsFile, double tmin, double tmax, Bool_t ignoreSOR,Bool_t ignoreEOR)
  : fTMin(0),fTMax(0),fFillNumber(0),fData(0),fFile2Process(dcsFile),fMap2Process(0)
{
  SetIgnoreSOR(ignoreSOR);
  SetIgnoreEOR(ignoreEOR);
  FillData(dcsFile,tmin,tmax);
}

//___________________________________________________________________
Bool_t AliLHCData::FillData(const TMap*   dcsMap,double tmin, double tmax)
{
  // process DCS map and fill all fields. 
  Clear();
  fMap2Process = dcsMap;
  FillData(tmin,tmax);
  return kTRUE;
}

//___________________________________________________________________
Bool_t AliLHCData::FillData(const Char_t* dcsFile,double tmin, double tmax)
{
  // process DCS file and fill all fields. 
  Clear();
  fFile2Process = dcsFile;
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
  int iEntry;
  if (IsSORIgnored()) tmin -= kTimeMargin;
  if (IsSORIgnored()) tmax += kTimeMargin;
  TObjArray* arr = GetDCSEntry(fgkDCSNames[kFillNum],iEntry,tmin,tmax);
  if (arr) SetFillNumber( ExtractInt( (AliDCSArray*)arr->At(iEntry), 0) );
  if (fFile2Process) delete arr; // array was created on demand
  //
  for (int ibm=0;ibm<2;ibm++) {
    //
    sprintf(buff,fgkDCSNames[kBunchConf],ibm+1);         // ----- declared bunch configuration
    FillBunchConfig(fBunchConfDecl[ibm], buff);
    //
    sprintf(buff,fgkDCSNames[kBunchLgtFillB],ibm+1);     // ----- measured bunch configuration
    FillBunchConfig(fBunchConfMeas[ibm], buff);
    //
    sprintf(buff,fgkDCSNames[kBunchLgt],ibm+1);          // ----- measured bunch lenghts
    FillBunchInfo(fBunchLengths[ibm],buff,ibm,kFALSE);  
    //
    sprintf(buff,fgkDCSNames[kIntBunchAv],ibm+1);        // ----- B-by-B intensities
    FillBunchInfo(fIntensPerBunch[ibm],buff,ibm,kTRUE);
    //
    //
    sprintf(buff,fgkDCSNames[kIntTot],ibm+1);            // ----- total intensities for beam 1 and 2
    FillScalarRecord(fIntensTotal[ibm], buff);
    //
    sprintf(buff,fgkDCSNames[kIntTotAv],ibm+1);          // ----- total intensities for beam 1 and 2 from B-by-B average
    FillScalarRecord(fIntensTotalAv[ibm], buff);
    //
    sprintf(buff,fgkDCSNames[kBeamSzEmittH],ibm+1);      // ----- H emittance for beam 1 and 2 
    FillScalarRecord(fEmittanceH[ibm], buff);
    //
    sprintf(buff,fgkDCSNames[kBeamSzEmittV],ibm+1);      // ----- V emittance for beam 1 and 2 
    FillScalarRecord(fEmittanceV[ibm], buff);
    //
    sprintf(buff ,fgkDCSNames[kBeamSzSigH],   ibm+1);    // ----- H sigmas and errors for beam 1 and 2 
    sprintf(buff1,fgkDCSNames[kBeamSzSigHErr],ibm+1);
    FillScalarRecord(fBeamSigmaH[ibm], buff, buff1);
    //
    sprintf(buff ,fgkDCSNames[kBeamSzSigV],   ibm+1);    // ----- V sigmas and errors for beam 1 and 2 
    sprintf(buff1,fgkDCSNames[kBeamSzSigVErr],ibm+1);
    FillScalarRecord(fBeamSigmaV[ibm], buff, buff1);
    //
  }
  //
  FlagInteractingBunches(fBunchConfMeas[0],fBunchConfMeas[1]);
  FlagInteractingBunches(fBunchConfDecl[0],fBunchConfDecl[1]);
  //
  for (int ilr=0;ilr<2;ilr++) {
    //
    sprintf(buff ,fgkDCSNames[kLumBunch], ilr ? 'R':'L');       // ---- BC-by-BC luminosity at IP2 and its error
    sprintf(buff1,fgkDCSNames[kLumBunchErr], ilr ? 'R':'L');
    FillBCLuminosities(fLuminPerBC[ilr], buff, buff1, 0); // BRAN L uses beam2 as a reference, BRAN R - beam1
    //
    sprintf(buff ,fgkDCSNames[kLumTot]   , ilr ? 'R':'L');       // ---- total luminosity at IP2 and its error
    sprintf(buff1,fgkDCSNames[kLumTotErr], ilr ? 'R':'L');
    FillScalarRecord(fLuminTotal[ilr], buff, buff1);
    //
    sprintf(buff ,fgkDCSNames[kLumAcqMode], ilr ? 'R':'L');      // ---- luminosity acquisition mode
    FillAcqMode(fLuminAcqMode[ilr], buff);
    //
    sprintf(buff, fgkDCSNames[kLumCrossAng]   , ilr ? 'R':'L');  //----- crossing angle at IP2 and its error
    sprintf(buff1,fgkDCSNames[kLumCrossAngErr], ilr ? 'R':'L');
    FillScalarRecord(fCrossAngle[ilr], buff, buff1);
    //    
  }
  //
  for (int icl=0;icl<kNCollimators;icl++) {             // ----- collimators positions
    for (int jaw=0;jaw<kNJaws;jaw++) {
      sprintf(buff,fgkDCSNames[kCollPos], fgkDCSColNames[icl],fgkDCSColJaws[jaw]);        
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
  return kTRUE;
}

//___________________________________________________________________
TObjArray* AliLHCData::GetDCSEntry(const char* key,int &entry,double tmin,double tmax) const
{
  // extract array from the DCS map or file and find the first entry within the time limits
  entry = -1;
  TObjArray* arr;
  if (fMap2Process) arr = (TObjArray*)fMap2Process->GetValue(key);
  else if (fFile2Process) {
    AliLHCReader rd;
    arr = rd.ReadSingleLHCDP(fFile2Process,key);
  }
  else {
    AliError("Neither DCS map nor DCS filename are set");
    return 0;  
  }
  //
  if (!arr || !arr->GetEntriesFast()) { 
    AliWarning(Form("No data for %s",key)); 
    if (fMap2Process) delete arr; // created on demand
    return 0;
  }
  int ntot = arr->GetEntriesFast();
  for (entry=0;entry<ntot;entry++) {
    AliDCSArray* ent = (AliDCSArray*)arr->At(entry);
    if (ent->GetTimeStamp()>=tmin && ent->GetTimeStamp()<=tmax) break;
  }
  if (entry==ntot) {
    entry = -1;
    TString str;
    str += AliLHCDipValD::TimeAsString(tmin);
    str += " : ";
    str += AliLHCDipValD::TimeAsString(tmax);
    AliWarning(Form("All entries for %s are outside the requested range:\n%s",key,str.Data()));
    if (fMap2Process) delete arr; // created on demand
    return 0;
  }
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
    if (TimeDifference(tstamp,curObj->GetTimeStamp())<0) break;
    prevObj = curObj;
  }
  if (!prevObj && nrec>0) prevObj = (AliLHCDipValI*)fData[start]; // if no exact match, return the 1st one
  return prevObj;
}

//___________________________________________________________________
Int_t AliLHCData::FillScalarRecord(int refs[2], const char* rec, const char* recErr)
{
  // fill record for scalar value, optionally accompanied by measurement error 
  //
  AliInfo(Form("Acquiring record: %s",rec));
  //
  TObjArray *arr=0,*arrE=0;
  Int_t nEntries=0,nEntriesE=0,iEntry=0,iEntryE=0;
  //
  refs[kStart] = fData.GetEntriesFast();
  refs[kNStor] = 0;
  //
  double tmin = fTMin;
  double tmax = fTMax;
  if (IsSORIgnored()) tmin -= kTimeMargin;
  if (IsSORIgnored()) tmax += kTimeMargin;
  //
  if ( !(arr=GetDCSEntry(rec,iEntry,tmin,tmax)) ) return -1;
  nEntries = arr->GetEntriesFast();
  //
  int dim = 1;
  if (recErr) {
    if ( !(arrE=GetDCSEntry(recErr,iEntryE,tmin,tmax)) ) nEntriesE = -999;
    else nEntriesE = arrE->GetEntriesFast();
    dim += 1;
  }
  //
  while (iEntry<nEntries) {
    AliDCSArray *dcsVal = (AliDCSArray*) arr->At(iEntry++);
    double tstamp = dcsVal->GetTimeStamp();
    if (tstamp>fTMax) break;
    //
    AliLHCDipValF* curValF = new AliLHCDipValF(dim,tstamp);  // start new period
    (*curValF)[0] = ExtractDouble(dcsVal,0);     // value
    //
    if (recErr) {
      double errVal = -1;
      while (iEntryE<nEntriesE) {       // try to find corresponding error
	AliDCSArray *dcsValE = (AliDCSArray*) arrE->At(iEntryE);
        double tstampE = dcsValE->GetTimeStamp();
        if (tstampE>fTMax) break;
        int tdif = TimeDifference(tstamp,tstampE);
        if (!tdif) { // error matches to value
          errVal = ExtractDouble(dcsValE,0);
	  iEntryE++; 
	  break;
	}
        else if (tdif>0) iEntryE++; // error time lags behind, read the next one
        else break;                 // error time is ahead of value, no error associated
      }
      (*curValF)[dim-1] = errVal;   // error
      curValF->SetLastSpecial();    // lable the last entry as an error
    }
    //
    fData.Add(curValF);
    refs[kNStor]++;
  }
  //
  if (fFile2Process) {
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
  Int_t nEntries,iEntry;
  //
  refs[kStart] = fData.GetEntriesFast();
  refs[kNStor] = 0;
  //
  double tmin = fTMin;
  double tmax = fTMax;
  if (IsSORIgnored()) tmin -= kTimeMargin;
  if (IsSORIgnored()) tmax += kTimeMargin;
  //
  if ( !(arr=GetDCSEntry(rec,iEntry,tmin,tmax)) ) return -1;
  nEntries = arr->GetEntriesFast();
  //
  AliLHCDipValI* prevRecI=0;
  while (iEntry<nEntries) {
    AliDCSArray *dcsVal = (AliDCSArray*) arr->At(iEntry++);
    double tstamp = dcsVal->GetTimeStamp();
    if (tstamp>fTMax) break;
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
  if (fFile2Process) delete arr;
  return refs[kNStor];
}
 
//___________________________________________________________________
Int_t AliLHCData::FillAcqMode(int refs[2],const char* rec)
{
  // fill acquisition mode
  //
  AliInfo(Form("Acquiring record: %s",rec));
  TObjArray *arr;
  Int_t nEntries,iEntry;
  //
  refs[kStart] = fData.GetEntriesFast();
  refs[kNStor] = 0;
  //
  double tmin = fTMin;
  double tmax = fTMax;
  if (IsSORIgnored()) tmin -= kTimeMargin;
  if (IsSORIgnored()) tmax += kTimeMargin;
  //
  if ( !(arr=GetDCSEntry(rec,iEntry,tmin,tmax)) ) return -1;
  nEntries = arr->GetEntriesFast();
  //
  AliLHCDipValI* prevRecI=0;
  while (iEntry<nEntries) {
    AliDCSArray *dcsVal = (AliDCSArray*) arr->At(iEntry++);
    double tstamp = dcsVal->GetTimeStamp();
    if (tstamp>fTMax) break;
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
  if (fFile2Process) delete arr;
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
  Int_t nEntries,iEntry;
  //
  refs[kStart] = fData.GetEntriesFast();
  refs[kNStor] = 0;
  //
  double tmin = fTMin;
  double tmax = fTMax;
  if (IsSORIgnored()) tmin -= kTimeMargin;
  if (IsSORIgnored()) tmax += kTimeMargin;
  //
  if ( !(arr=GetDCSEntry(rec,iEntry,tmin,tmax)) ) return -1;
  nEntries = arr->GetEntriesFast();
  //
  while (iEntry<nEntries) {
    AliDCSArray *dcsVal = (AliDCSArray*) arr->At(iEntry++);
    double tstamp = dcsVal->GetTimeStamp();
    if (tstamp>fTMax) break;
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
  if (fFile2Process) delete arr;
  return refs[kNStor];
}

//___________________________________________________________________
Int_t AliLHCData::FillBunchInfo(int refs[2],const char* rec, int ibm, Bool_t inRealSlots)
{
  // fill bunch properties for beam ibm
  // if inRealSlots = true, then the value is taken from bunchRFbucket/10, otherwise, the value 
  // for the i-th bunch is taken from the i-th element
  //
  AliInfo(Form("Acquiring record: %s",rec));
  TObjArray *arr;
  Int_t nEntries,iEntry;
  //
  refs[kStart] = fData.GetEntriesFast();
  refs[kNStor] = 0;
  //
  double tmin = fTMin;
  double tmax = fTMax;
  if (IsSORIgnored()) tmin -= kTimeMargin;
  if (IsSORIgnored()) tmax += kTimeMargin;
  //
  if ( !(arr=GetDCSEntry(rec,iEntry,tmin,tmax)) ) return -1;
  nEntries = arr->GetEntriesFast();
  //
  while (iEntry<nEntries) {
    AliDCSArray *dcsVal = (AliDCSArray*) arr->At(iEntry++);
    double tstamp = dcsVal->GetTimeStamp();
    if (tstamp>fTMax) break;
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
    AliLHCDipValF* curValF = new AliLHCDipValF(nbunch,tstamp);
    for (int i=nbunch;i--;) {
      int ind = inRealSlots ? (*bconf)[i]/10 : i;
      if (ind>nSlots) {
	AliError(Form("Bunch %d refers to wrong slot %d, set to -1",i,(*bconf)[i]));
	(*curValF)[i] = -1;
      }
      else (*curValF)[i] = dcsArr[ind];
    }
    fData.Add(curValF);
    refs[kNStor]++;
  }
  if (fFile2Process) delete arr;
  return refs[kNStor];
  //
}
 
//___________________________________________________________________
Int_t AliLHCData::FillBCLuminosities(int refs[2],const char* rec, const char* recErr, int useBeam)
{
  // fill luminosities per bunch crossing
  //
  AliInfo(Form("Acquiring record: %s",rec));
  TObjArray *arr,*arrE=0;
  Int_t nEntries=0,nEntriesE=0,iEntry=0,iEntryE=0;
  //
  refs[kStart] = fData.GetEntriesFast();
  refs[kNStor] = 0;
  //
  double tmin = fTMin;
  double tmax = fTMax;
  if (IsSORIgnored()) tmin -= kTimeMargin;
  if (IsSORIgnored()) tmax += kTimeMargin;
  //
  if ( !(arr=GetDCSEntry(rec,iEntry,tmin,tmax)) ) return -1;
  nEntries = arr->GetEntriesFast();
  //
  while (iEntry<nEntries) {
    AliDCSArray *dcsVal = (AliDCSArray*) arr->At(iEntry++);
    double tstamp = dcsVal->GetTimeStamp();
    if (tstamp>fTMax) break;
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
      AliWarning(Form("Bunch conf. for beam%d has no marked interacting bunches, store all luminosity for all filled bunches",useBeam));
      dim = nbunch;
    }
    else { // count number of interacting bunches
      for (int i=nbunch;i--;) if ((*bconf)[i]<0) dim++;
    }
    //
    if (recErr) {
      if ( !(arrE=GetDCSEntry(recErr,iEntryE,tmin,tmax)) || iEntryE<0 ) nEntriesE = -999;
      else nEntriesE = arrE->GetEntriesFast();
      dim += 1;
    }
    AliLHCDipValF* curValF = new AliLHCDipValF(dim,tstamp);
    int cnt = 0;
    for (int i=0;i<nbunch;i++) {
      int slot = (*bconf)[i];
      if (bconf->IsProcessed1() && slot>0) continue;
      //
      int ind = TMath::Abs(slot)/10;
      if (ind>nSlots) {
	AliError(Form("Bunch %d refers to wrong slot %d, set to -1",cnt,slot));
	(*curValF)[cnt] = -1;
      }
      else (*curValF)[cnt] = dcsArr[ind];
      cnt++;
    }
    //
    if (recErr) {
      double errVal = -1;
      while (iEntryE<nEntriesE) {       // try to find corresponding error
	AliDCSArray *dcsValE = (AliDCSArray*) arrE->At(iEntryE);
	double tstamp1 = dcsValE->GetTimeStamp();
	if (tstamp1>fTMax) break;
	int tdif = TimeDifference(tstamp,tstamp1);
	if (!tdif) { // error matches to value
	  errVal = dcsValE->GetDouble()[0];
	  iEntryE++; 
	  break;
	}
	else if (tdif>0) iEntryE++; // error time lags behind, read the next one
	else break;                 // error time is ahead of value, no error associated
      }
      (*curValF)[dim-1] = errVal;   // error
      curValF->SetLastSpecial();    // lable the last entry as an error
    }
    //
    fData.Add(curValF);
    refs[kNStor]++;
  }
  if (fFile2Process) {
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
  Bool_t full = kTRUE;
  if (!opts.Contains("f")) {
    printf("Use Print(\"f\") to print full info\n");
    printf("Printing short info:\n<RecordType>(number of records): <TimeStamp, value> for 1st record only\n");
    printf("Ignoring strict time of SOR: %s, EOR: %s\n",IsSORIgnored()?"ON":"OFF",IsEORIgnored()?"ON":"OFF");
    full = kFALSE;
  }
  TString sdtmn = AliLHCDipValI::TimeAsString(fTMin,utcTime);
  TString sdtmx = AliLHCDipValI::TimeAsString(fTMax,utcTime);
  printf("Fill#%6d Validity: %s - %s (%s)\n",fFillNumber,sdtmn.Data(),sdtmx.Data(),utcTime ? "UTC":"LOC");
  //
  printf("********** SETTINGS FROM RUN CONTROL **********\n");
  //
  printf("* %-38s","Injection Scheme");
  PrintAux(full,fRCInjScheme,opts);
  //
  printf("* %-38s","Beta Star");
  PrintAux(full,fRCBeta,opts);
  //
  printf("* %-38s","Horisontal Crossing Angle");
  PrintAux(full,fRCAngH,opts);
  //
  printf("* %-38s","Vertical   Crossing Angle");
  PrintAux(full,fRCAngV,opts);
  //
  for (int ib=0;ib<2;ib++) {
    printf("* Beam%d filling  [- interacts at IR2!]  ",ib+1);
    PrintAux(full,fBunchConfDecl[ib],opts);
  }
  //
  printf("\n**********       MEASURED DATA       **********\n");
  //
  for (int ib=0;ib<2;ib++) {
    printf("* Beam%d filling  [- interacts at IR2!]  ",ib+1);
    PrintAux(full,fBunchConfMeas[ib],opts);
  } 
  //
  for (int ib=0;ib<2;ib++) {
    printf("* Beam%d total intensity                 ",ib+1);
    PrintAux(full,fIntensTotal[ib],opts);
  } 
  //
  for (int ib=0;ib<2;ib++) {
    printf("* Beam%d total intensity (bunch average) ",ib+1);
    PrintAux(full,fIntensTotalAv[ib],opts);
  } 
  //
  for (int ib=0;ib<2;ib++) {
    printf("* Beam%d intensity per bunch             ",ib+1);
    PrintAux(full,fIntensPerBunch[ib],opts);
  }
  //
  for (int ib=0;ib<2;ib++) {
    printf("* Beam%d bunch lengths                   ",ib+1);
    PrintAux(full,fBunchLengths[ib],opts);
  } 
  //
  for (int ib=0;ib<2;ib++) {
    printf("* Beam%d Horisontal emittance            ",ib+1);
    PrintAux(full,fEmittanceH[ib],opts);
  }
  //
  for (int ib=0;ib<2;ib++) {
    printf("* Beam%d Vertical emittance              ",ib+1);
    PrintAux(full,fEmittanceV[ib],opts);
  }
  //
  for (int ib=0;ib<2;ib++) {
    printf("* Beam%d Horisontal sigma                ",ib+1);
    PrintAux(full,fBeamSigmaH[ib],opts);
  }
  //
  for (int ib=0;ib<2;ib++) {
    printf("* Beam%d Vertical sigma                  ",ib+1);
    PrintAux(full,fBeamSigmaV[ib],opts);
  }
  //
  for (int lr=0;lr<2;lr++) {
    printf("* Total luminosity from BRANB_4%c2       ",lr ? 'R':'L');
    PrintAux(full,fLuminTotal[lr],opts);
  } 
  //
  for (int lr=0;lr<2;lr++) {
    printf("* Luminosity acq.mode, BRANB_4%c2        ",lr ? 'R':'L');
    PrintAux(full,fLuminAcqMode[lr],opts+"bit");
  } 
  //
  for (int lr=0;lr<2;lr++) {
    printf("* Luminosity per BC from BRANB_4%c2      ",lr ? 'R':'L');
    PrintAux(full,fLuminPerBC[lr],opts);
  }
  //
  for (int lr=0;lr<2;lr++) {
    printf("* Crossing angle, side %c                ",lr ? 'R':'L');
    PrintAux(full,fCrossAngle[lr],opts);
  }
  //
  for (int coll=0;coll<kNCollimators;coll++)
    for (int jaw=0;jaw<kNJaws;jaw++) {
      printf("* Collimator %10s:%16s",fgkDCSColNames[coll],fgkDCSColJaws[jaw]);
      PrintAux(full,fCollimators[coll][jaw],opts);
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
  fFile2Process = 0;
  fMap2Process  = 0;
  //
  for (int i=2;i--;) {
    fRCInjScheme[i] = 0;
    fRCBeta[i] = 0;
    fRCAngH[i] = 0;
    fRCAngV[i] = 0;
    //
    for (int icl=kNCollimators;icl--;) for (int jaw=kNJaws;jaw--;) fCollimators[icl][jaw][i]=0;
    //
    for (int j=2;j--;) {
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
  const int kMaxSlots  = 3564;
  const int kOffsBeam1 = 346;
  const int kOffsBeam2 = 3019;
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
      int slot2 =-1, slot1  = (bunch1/10 + kOffsBeam1)%kMaxSlots;
      for (i2=0;i2<nb2;i2++) {
	bunch2 = TMath::Abs((*bm2)[i2]);
	slot2 = (bunch2/10 + kOffsBeam2)%kMaxSlots;
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
