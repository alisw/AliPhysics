// @(#) $Id$

/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 *                                                                        *
 * Primary Authors: Matthias Richter <Matthias.Richter@ift.uib.no>        *
 *                  Kenneth Aamodt   <Kenneth.aamodt@ift.uib.no>          *
 *                  for The ALICE HLT Project.                            *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/** @file   AliHLTTPCPad.cxx
    @author Matthias Richter, Kenneth Aamodt
    @date   
    @brief  Container Class for TPC Pads.
*/

// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#if __GNUC__>= 3
using namespace std;
#endif

#include <cerrno>
#include "AliHLTTPCPad.h"
#include "AliHLTStdIncludes.h"


//added by kenneth
#include "AliHLTTPCTransform.h"
#include "AliHLTTPCClusters.h"
#include <sys/time.h>
#include "TMath.h"
#include "TFile.h"
//------------------------------

/** margin for the base line be re-avaluated */
#define ALIHLTPAD_BASELINE_MARGIN (2*fAverage)

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTTPCPad)

AliHLTTPCPad::AliHLTTPCPad()
  :
  fClusterCandidates(),
  fUsedClusterCandidates(),
  fRowNo(-1),
  fPadNo(-1),
  fThreshold(0),
  fAverage(-1),
  fNofEvents(0),
  fSum(0),
  fCount(0),
  fTotal(0),
  fBLMax(-1),
  fBLMaxBin(-1),
  fBLMin(-1),
  fBLMinBin(-1),
  fFirstBLBin(0),
  fNofBins(0),
  fReadPos(0),
  fpRawData(NULL),
  fDataSignals(NULL),
  fSignalPositionArray(NULL),
  fSizeOfSignalPositionArray(0),
  fNSigmaThreshold(0),
  fSignalThreshold(0),
  fModeSwitch(0),
  fNGoodSignalsSent(0),
  fDebugHistoBeforeZS(NULL),
  fDebugHistoAfterZS(NULL)
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
  //  HLTInfo("Entering default constructor");
  fDataSignals= new AliHLTTPCSignal_t[AliHLTTPCTransform::GetNTimeBins()];
  memset( fDataSignals, 0xFF, sizeof(Int_t)*(AliHLTTPCTransform::GetNTimeBins()));
  
  fSignalPositionArray= new AliHLTTPCSignal_t[AliHLTTPCTransform::GetNTimeBins()];
  memset( fSignalPositionArray, 0xFF, sizeof(Int_t)*(AliHLTTPCTransform::GetNTimeBins()));
  fSizeOfSignalPositionArray=0;

}

AliHLTTPCPad::AliHLTTPCPad(Int_t mode)
  :
  fClusterCandidates(0),
  fUsedClusterCandidates(0),
  fRowNo(-1),
  fPadNo(-1),
  fThreshold(0),
  fAverage(-1),
  fNofEvents(0),
  fSum(0),
  fCount(0),
  fTotal(0),
  fBLMax(-1),
  fBLMaxBin(-1),
  fBLMin(-1),
  fBLMinBin(-1),
  fFirstBLBin(0),
  fNofBins(0),
  fReadPos(0),
  fpRawData(NULL),
  fDataSignals(NULL),
  fSignalPositionArray(NULL),
  fSizeOfSignalPositionArray(0),
  fNSigmaThreshold(0),
  fSignalThreshold(0),
  fModeSwitch(mode),
  fNGoodSignalsSent(0),
  fDebugHistoBeforeZS(NULL),
  fDebugHistoAfterZS(NULL)
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

AliHLTTPCPad::AliHLTTPCPad(Int_t offset, Int_t nofBins)
  :
  fClusterCandidates(),
  fUsedClusterCandidates(),
  fRowNo(-1),
  fPadNo(-1),
  fThreshold(0),
  fAverage(-1),
  fNofEvents(0),
  fSum(0),
  fCount(0),
  fTotal(0),
  fBLMax(-1),
  fBLMaxBin(-1),
  fBLMin(-1),
  fBLMinBin(-1),
  fFirstBLBin(offset),
  fNofBins(nofBins),
  fReadPos(0),
  fpRawData(NULL),
  fDataSignals(NULL),
  fSignalPositionArray(NULL),
  fSizeOfSignalPositionArray(0),
  fNSigmaThreshold(0),
  fSignalThreshold(0),
  fModeSwitch(0),
  fNGoodSignalsSent(0),
  fDebugHistoBeforeZS(NULL),
  fDebugHistoAfterZS(NULL)
{
  // see header file for class documentation
}

AliHLTTPCPad::~AliHLTTPCPad()
{
  // see header file for class documentation
  if (fpRawData) {
    HLTWarning("event data acquisition not stopped");
    StopEvent();
  }
  if (fDataSignals) {
    delete [] fDataSignals;
    fDataSignals=NULL;
  }
  if (fSignalPositionArray) {
    delete [] fSignalPositionArray;
    fSignalPositionArray=NULL;
  }
  if(fDebugHistoBeforeZS){
    delete fDebugHistoBeforeZS;
    fDebugHistoBeforeZS=NULL;
  }
  if(fDebugHistoAfterZS){
    delete fDebugHistoAfterZS;
    fDebugHistoAfterZS=NULL;
  }
}

Int_t AliHLTTPCPad::SetID(Int_t rowno, Int_t padno)
{
  // see header file for class documentation
  fRowNo=rowno;
  fPadNo=padno;

#if DebugHisto
  char *nameBefore;
  sprintf(nameBefore,"beforeRow%dPad%d",fRowNo,fPadNo);
  char *nameAfter;
  sprintf(nameAfter,"afterRow%dPad%d",fRowNo,fPadNo);
  fDebugHistoBeforeZS = new TH1F(nameBefore,nameBefore,1024,0,1024);
  fDebugHistoAfterZS = new TH1F(nameAfter,nameAfter,1024,0,1024);
#endif

  return 0;
}

Int_t AliHLTTPCPad::StartEvent()
{
  // see header file for class documentation
  Int_t iResult=0;
  if (fpRawData==NULL) {
    fBLMax=-1;
    fBLMaxBin=-1;
    fBLMin=-1;
    fBLMinBin=-1;
    fSum=0;
    fCount=0;
    fTotal=0;
    if (fNofBins>0) {
      fpRawData=new AliHLTTPCSignal_t[fNofBins];
      if (fpRawData) {
	for (int i=0; i<fNofBins; i++) fpRawData[i]=-1;
      } else {
	HLTError("memory allocation failed");
	iResult=-ENOMEM;
      }
    }
  } else {
    HLTWarning("event data acquisition already started");
    iResult=-EALREADY;
  }
  return iResult;
}

Int_t AliHLTTPCPad::CalculateBaseLine(Int_t reqMinCount)
{
  // see header file for class documentation
  Int_t iResult=0;
  AliHLTTPCSignal_t avBackup=fAverage;
  //HLTDebug("reqMinCount=%d fCount=%d fTotal=%d fSum=%d fBLMax=%d fBLMin=%d", reqMinCount, fCount, fTotal, fSum, fBLMax, fBLMin);
  if (fCount>=reqMinCount && fCount>=fTotal/2) {
    fAverage=fCount>0?fSum/fCount:0;
    if (fAverage>0) {
      //HLTDebug("average for current event %d (%d - %d)", fAverage, fBLMax, fBLMin);
      fCount=0;fSum=-1;
      if (fBLMax>ALIHLTPAD_BASELINE_MARGIN) {
	// calculate again
	//HLTDebug("maximum value %d exceeds margin for base line (%d) "
	//	 "-> re-evaluate base line", fBLMax, ALIHLTPAD_BASELINE_MARGIN);
	if (fpRawData) {
	  for (Int_t i=fFirstBLBin; i<fNofBins; i++)
	    if (fpRawData[i]>=0) AddBaseLineValue(i, fpRawData[i]);
	  if (fCount>0 && fCount>=reqMinCount && fCount>=fTotal/2) {
	    fAverage=fSum/fCount;
	    //HLTDebug("new average %d", fAverage);
	  } else {
	    // 	    HLTDebug("baseline re-eveluation skipped because of to few "
	    // 		       "contributing bins: total=%d, contributing=%d, req=%d"
	    // 		       "\ndata might be already zero suppressed"
	    // 		       , fTotal, fCount, reqMinCount);
	    iResult=-ENODATA;
	  }
	  fCount=0;fSum=-1;
	} else {
	  HLTError("missing raw data for base line calculation");
	  iResult=-ENOBUFS;
	}
      }
      if (iResult>=0) {
	// calculate average for all events
	fAverage=((avBackup*fNofEvents)+fAverage)/(fNofEvents+1);
	//HLTDebug("base line average for %d event(s): %d", fNofEvents+1, fAverage);
      } else {
	fAverage=avBackup;      
      }
    } else {
      fAverage=avBackup;
    }
  } else {
    //     HLTDebug("baseline calculation skipped because of to few contributing "
    // 	       "bins: total=%d, contributing=%d, required=%d \ndata might be "
    // 	       "already zero suppressed", fTotal, fCount, reqMinCount);
  }

  return iResult;
}

Int_t AliHLTTPCPad::StopEvent()
{
  // see header file for class documentation
  Int_t iResult=0;
  if (fpRawData) {
    AliHLTTPCSignal_t* pData=fpRawData;
    fpRawData=NULL;
    delete [] pData;
    fTotal=0;
    fNofEvents++;
    Rewind();
  } else if (fNofBins>0) {
    HLTError("event data acquisition not started");
    iResult=-EBADF;
  }
  return iResult;
}

Int_t AliHLTTPCPad::ResetHistory()
{
  // see header file for class documentation
  Int_t iResult=0;
  fAverage=-1;
  fNofEvents=0;
  return iResult;
}

Int_t AliHLTTPCPad::SetThreshold(AliHLTTPCSignal_t thresh)
{
  // see header file for class documentation
  Int_t iResult=0;
  fThreshold=thresh;
  return iResult;
}

Int_t AliHLTTPCPad::AddBaseLineValue(Int_t bin, AliHLTTPCSignal_t value)
{
  // see header file for class documentation
  Int_t iResult=0;
  if (bin>=fFirstBLBin) {
    if (fAverage<0 || value<ALIHLTPAD_BASELINE_MARGIN) {
      // add to the current sum and count
      fSum+=value;
      fCount++;
      if (fBLMax<value) {
	// keep the maximum value for later quality control of the base 
	// line calculation
	fBLMax=value;
	fBLMaxBin=bin;
      }
      if (fBLMin<0 || fBLMin>value) {
	// keep the minimum value for later quality control of the base 
	// line calculation
	fBLMin=value;
	fBLMinBin=bin;
      }
    } else {
      //       HLTDebug("ignoring value %d (bin %d) for base line calculation "
      // 	       "(current average is %d)",
      // 	       value, bin, fAverage);
    }
  }
  return iResult;
}

Int_t AliHLTTPCPad::SetRawData(Int_t bin, AliHLTTPCSignal_t value)
{
  // see header file for class documentation
  //  printf("Row: %d    Pad: %d  Time: %d Charge %d", fRowNo, fPadNo, bin, value);
  Int_t iResult=0;
  if (fpRawData) {
    if (bin<fNofBins) {
      if (value>=0) {
	if (fpRawData[bin]<0) {
	  AddBaseLineValue(bin, value);
	  fTotal++;
	} else {
	  // ignore value for average calculation
	  HLTWarning("overriding content of bin %d (%d)", bin, fpRawData[bin]);
	}
	fpRawData[bin]=value;
      } else {
	HLTWarning("ignoring neg. raw data");
      }
    } else {
      HLTWarning("bin %d out of range (%d)", bin, fNofBins);
      iResult=-ERANGE;
    }
  } else if (fNofBins>0) {
    HLTError("event cycle not started");
    iResult=-EBADF;
  }
  return iResult;
}

Int_t AliHLTTPCPad::Next(Int_t bZeroSuppression) 
{
  // see header file for class documentation
  if (fpRawData==NULL) return 0;
  Int_t iResult=fReadPos<fNofBins;
  if (iResult>0 && (iResult=(++fReadPos<fNofBins))>0) {
    if (bZeroSuppression) {
      while ((iResult=(fReadPos<fNofBins))>0 &&
	     GetCorrectedData(fReadPos)<=0)
	fReadPos++;
    }
  }
  return iResult;
}

Int_t AliHLTTPCPad::Rewind(Int_t bZeroSuppression)
{
  // see header file for class documentation
  fReadPos=(bZeroSuppression>0?0:fFirstBLBin)-1;
  return Next(bZeroSuppression);
}

AliHLTTPCSignal_t AliHLTTPCPad::GetRawData(Int_t bin) const
{
  // see header file for class documentation
  AliHLTTPCSignal_t data=0;
  if (fpRawData) {
    if (bin<fNofBins) {
      data=fpRawData[bin];
    } else {
      HLTWarning("requested bin %d out of range (%d)", bin, fNofBins);
    }
  } else if (fNofBins>0) {
    HLTWarning("data only available within event cycle");
  }
  return data;
}

AliHLTTPCSignal_t AliHLTTPCPad::GetCorrectedData(Int_t bin) const
{
  // see header file for class documentation
  AliHLTTPCSignal_t data=GetRawData(bin)-GetBaseLine(bin);
  AliHLTTPCSignal_t prev=0;
  if (bin>1) prev=GetRawData(bin-1)-GetBaseLine(bin-1);
  AliHLTTPCSignal_t succ=0;
  if (bin+1<GetSize()) succ=GetRawData(bin+1)-GetBaseLine(bin+1);
  if (fThreshold>0) {
    data-=fThreshold;
    prev-=fThreshold;
    succ-=fThreshold;
  }
 
  // case 1:
  // the signal is below the base-line and threshold
  if (data<0) data=0;

  //case 2:
  // the neighboring bins are both below base-line/threshold
  // a real signal is always more than one bin wide because of the shaper 
  if (prev<=0 && succ<=0) data=0;
 
  // case 3:
  // the bin is inside the range of ignored bins
  if (bin<fFirstBLBin) data=0;
  //HLTDebug("fReadPos=%d data=%d threshold=%d raw data=%d base line=%d", fReadPos, data, fThreshold, GetRawData(bin), GetBaseLine(bin));
  return data;
}

AliHLTTPCSignal_t AliHLTTPCPad::GetBaseLine(Int_t /*bin*/) const //TODO: Why is bin being ignored?
{
  // see header file for class documentation
  AliHLTTPCSignal_t val=0;
  if (fAverage>0) {
    // we take the minumum value as the base line if it doesn't differ from
    // the average to much
    val=fAverage;
#ifdef KEEP_NOISE
    const AliHLTTPCSignal_t kMaxDifference=15;
    if ((fAverage-fBLMin)<=kMaxDifference) val=fBLMin;
    else val>kMaxDifference?val-=kMaxDifference:0;
#endif
  }
  if (val<0) {
    // here we should never get
    val=0;
    HLTFatal("wrong base line value");
  }
  return val;
}

AliHLTTPCSignal_t AliHLTTPCPad::GetAverage() const
{
  // see header file for class documentation
  return fAverage>0?fAverage:0;
}

Float_t AliHLTTPCPad::GetOccupancy() const
{
  // see header file for class documentation
  Float_t occupancy=0;
  if (fpRawData && fNofBins>0) {
    for (Int_t i=fFirstBLBin; i<fNofBins; i++) {
      if (GetCorrectedData(i)>0) occupancy+=1;
    }
    if (fNofBins-fFirstBLBin>0)
      occupancy/=fNofBins-fFirstBLBin;
  }
  return occupancy;
}

Float_t AliHLTTPCPad::GetAveragedOccupancy() const
{
  // see header file for class documentation

  // history is not yet implemented
  return GetOccupancy();
}
void AliHLTTPCPad::PrintRawData()
{
  // see header file for class documentation
  for(Int_t bin=0;bin<AliHLTTPCTransform::GetNTimeBins();bin++){
    if(GetDataSignal(bin)>0)
      cout<<fRowNo<<"\t"<<fPadNo<<"\t"<<bin<<"\t"<<GetDataSignal(bin)<<endl;;
  }
  //  cout<<"bins: "<<AliHLTTPCTransform::GetNTimeBins()<<endl;
}

void AliHLTTPCPad::ClearCandidates(){
  fClusterCandidates.clear();
  fUsedClusterCandidates.clear();
}

void AliHLTTPCPad::SetDataToDefault()
{
  // see header file for class documentation
  if(fDataSignals && fSignalPositionArray){
    for(Int_t i =0;i<fSizeOfSignalPositionArray;i++){
      fDataSignals[fSignalPositionArray[i]]=-1;
    }
    fSizeOfSignalPositionArray=0;
  }
}

void AliHLTTPCPad::SetDataSignal(Int_t bin,Int_t signal)
{
  // see header file for class documentation
  fDataSignals[bin]=signal;
  fSignalPositionArray[fSizeOfSignalPositionArray]=bin;
  fSizeOfSignalPositionArray++;
#if DebugHisto 
  fDebugHistoBeforeZS->Fill(bin,signal);
#endif
}

Bool_t AliHLTTPCPad::GetNextGoodSignal(Int_t &time, Int_t &signal ){
  /*  for(Int_t i=70;i<900;i++){
    if(fDataSignals[i]>0){
      printf("Signals which are good: Bin: %d Signal: %d\n",i,fDataSignals[i]);
    }
    }*/
  if(fNGoodSignalsSent<fSizeOfSignalPositionArray&&fSizeOfSignalPositionArray>0){
    time = fSignalPositionArray[fNGoodSignalsSent];
    signal = GetDataSignal(time);
    //    printf("GoodSignal: Row: %d Pad: %d time %d  signal %d  signalsSent: %d\n",fRowNo,fPadNo,fSignalPositionArray[fNGoodSignalsSent],GetDataSignal(time), fNGoodSignalsSent);
    fNGoodSignalsSent++;
    return kTRUE;
  }
  return kFALSE;
}

Int_t AliHLTTPCPad::GetDataSignal(Int_t bin) const
{
  // see header file for class documentation
  return fDataSignals[bin];
}

void AliHLTTPCPad::ZeroSuppress(Double_t nRMS, Int_t threshold, Int_t reqMinPoint, Int_t beginTime, Int_t endTime, Int_t timebinsLeft, Int_t timebinsRight, Int_t valueUnderAverage){
  //see headerfile for documentation
 
  //HLTDebug("In Pad: nRMS=%d, threshold=%d, reqMinPoint=%d, beginTime=%d, endTime=%d, timebinsLeft=%d timebinsRight=%d valueUnderAverage=%d \n",nRMS,threshold,reqMinPoint,beginTime,endTime,timebinsLeft,timebinsRight,valueUnderAverage);

  Bool_t useRMS= kFALSE;
  if(nRMS>0){
    useRMS=kTRUE;
    if(threshold>0){
      HLTInfo("Both RMSThreshold and SignalThreshold defined, using RMSThreshold");
    }
  }
  if(threshold<1 && nRMS<=0){
    //setting the data to -1 for this pad
    HLTInfo("Neither of RMSThreshold and SignalThreshold set, zerosuppression aborted");
    return;
  }
 
  Int_t fThresholdUsed=threshold;
 
  Int_t nAdded=0;
  Int_t sumNAdded=0;
  fSizeOfSignalPositionArray=0;
  if(useRMS){
    for(Int_t i=beginTime;i<endTime+1;i++){
      if(fDataSignals[i]>0){
	nAdded++;
	sumNAdded+=fDataSignals[i]*fDataSignals[i];
      }
    }
  }
  else if(threshold>0){
    for(Int_t i=beginTime;i<endTime+1;i++){
      if(fDataSignals[i]>0){
	nAdded++;
	sumNAdded+=fDataSignals[i];
      }
    }
  }
  else{
    HLTFatal("This should never happen because this is tested earlier in the code.(nRMSThreshold<1&&signal-threshold<1)");
  }
  if(nAdded<reqMinPoint){
    HLTInfo("Number of signals is less than required, zero suppression aborted");
    return;
  }
 
  if(nAdded==0){
    HLTInfo("No signals added for this pad, zerosuppression aborted: pad %d row %d",fPadNo,fRowNo);
    return;
  }
  // HLTInfo("sumNAdded=%d    nAdded=%d pad %d ",sumNAdded,nAdded,fPadNo);
  Double_t averageValue=(Double_t)sumNAdded/nAdded;//true average for threshold approach, average of signals squared for rms approach
 
  //  Double_t rms=0;
  if(useRMS){
    //Calculate the RMS
    if(averageValue>0){
      fThresholdUsed =(Int_t)(TMath::Sqrt(averageValue)*nRMS);
    }
    else{
      HLTFatal("average value in ZeroSuppression less than 0, investigation needed. This should never happen");
    }
  }
  else{
    fThresholdUsed = (Int_t)(averageValue + threshold); 
  }

  averageValue = 55.3;
  // Do zero suppression on the adc values within [beginTime,endTime]
  for(Int_t i=beginTime;i<endTime;i++){
    if(fDataSignals[i]>fThresholdUsed){
      //  HLTInfo("Signal Larger in pad %d time %d signal %d  ,   threshold: %d  averageValue %e",fPadNo,i,fDataSignals[i],fThresholdUsed, averageValue);
      Int_t firstSignalTime=i;
      for(Int_t left=1;left<timebinsLeft;left++){//looking 5 to the left of the signal to add tail
	if(fDataSignals[i-left]-averageValue+valueUnderAverage>0&&i-left>=beginTime){
	  firstSignalTime--;
	}
	else{
	  break;
	}
      }
      Int_t lastSignalTime=i;
      for(Int_t right=1;right<timebinsRight;right++){//looking 5 to the left of the signal to add tail
	if(fDataSignals[i+right]-averageValue+valueUnderAverage>0&&i+right<endTime){
	  lastSignalTime++;
	}
	else{
	  break;
	}	
      }
      for(Int_t t=firstSignalTime;t<lastSignalTime;t++){
	//	cout<<"Row: "<<fRowNo<<" Pad: "<<fPadNo<<"   Adding to tmebin: "<<t<<" signal: "<<(AliHLTTPCSignal_t)(fDataSignals[t]-averageValue + valueUnderAverage)<<endl;
	fDataSignals[t]=(AliHLTTPCSignal_t)(fDataSignals[t]-averageValue + valueUnderAverage);
	//	cout<<"Adding to signalPosition array bin number: "<<fSizeOfSignalPositionArray<<"    timebin number: "<<t<<endl;
	fSignalPositionArray[fSizeOfSignalPositionArray]=t;
	fSizeOfSignalPositionArray++;
	//	cout<<"Number of signals added so far: "<<fSizeOfSignalPositionArray<<"     firstSignalTimeBin: "<<firstSignalTime<<"     lastSignalTimeBin: "<<lastSignalTime<<endl;
	/*	if(fRowNo==29&&fPadNo==58){
	  cout<<"Signal added: Row: "<<fRowNo<<" Pad: "<<fPadNo<<"  Time: "<<t<<" signal: "<<fDataSignals[t]<<"  #signals: "<<fSizeOfSignalPositionArray<<endl;
	}
	*/
#if DebugHisto
      fDebugHistoAfterZS->Fill(t,fDataSignals[t]);
#endif
      }
      i+=lastSignalTime;
    }
    else{
      fDataSignals[i]=-1;
    }
  }
}

void AliHLTTPCPad::AddClusterCandidate(AliHLTTPCClusters candidate){
  fClusterCandidates.push_back(candidate);
  fUsedClusterCandidates.push_back(0);
}

void AliHLTTPCPad::SaveHistograms(){
#if DebugHisto
  if(fSizeOfSignalPositionArray==0){
    return;
  }
  char* filename;
  sprintf(filename,"/afsuser/kenneth/SimpleComponentWrapper/histos/HistogramsRow%dPad%d.root",fRowNo,fPadNo);
  TFile file(filename,"RECREATE");
  fDebugHistoBeforeZS->Write();
  fDebugHistoAfterZS->Write();
  file.Close();
#endif
}

void AliHLTTPCPad::FindClusterCandidates()
{
  // see header file for class documentation
  /*
  if(fSizeOfSignalPositionArray<2){
    return;
  }

  if(fNSigmaThreshold>0){
    ZeroSuppress(fNSigmaThreshold);
  }
  else if(fSignalThreshold>0){
    ZeroSuppress((Double_t)0,(Int_t)fSignalThreshold);
  }
  UInt_t seqcharge=0;
  UInt_t seqaverage=0;
  UInt_t seqerror=0;
  vector<Int_t> tmpPos;
  vector<Int_t> tmpSig;
  UInt_t isFalling=0;

  for(Int_t pos=fSizeOfSignalPositionArray-2;pos>=0;pos--){
    if(fSignalPositionArray[pos]==fSignalPositionArray[pos+1]+1){
      seqcharge+=fDataSignals[fSignalPositionArray[pos+1]];	
      seqaverage += fSignalPositionArray[pos+1]*fDataSignals[fSignalPositionArray[pos+1]];
      seqerror += fSignalPositionArray[pos+1]*fSignalPositionArray[pos+1]*fDataSignals[fSignalPositionArray[pos+1]];
	  
      tmpPos.push_back(fSignalPositionArray[pos+1]);
      tmpSig.push_back(fDataSignals[fSignalPositionArray[pos+1]]);

      if(fDataSignals[fSignalPositionArray[pos+1]]>fDataSignals[fSignalPositionArray[pos]]){
	isFalling=1;
      }
      if(fDataSignals[fSignalPositionArray[pos+1]]<fDataSignals[fSignalPositionArray[pos]]&&isFalling){
	Int_t seqmean=0;
	seqmean = seqaverage/seqcharge;
	
	//Calculate mean in pad direction:
	Int_t padmean = seqcharge*fPadNo;
	Int_t paderror = fPadNo*padmean;
	AliHLTTPCClusters candidate;
	candidate.fTotalCharge   = seqcharge;
	candidate.fPad       = padmean;
	candidate.fPad2      = paderror;
	candidate.fTime      = seqaverage;
	candidate.fTime2     = seqerror;
	candidate.fMean          = seqmean;
	candidate.fLastMergedPad = fPadNo;
	fClusterCandidates.push_back(candidate);
	fUsedClusterCandidates.push_back(0);
	isFalling=0;
	seqcharge=0;
	seqaverage=0;
	seqerror=0;

	tmpPos.clear();
	tmpSig.clear();

	continue;
      }
	 
      if(pos<1){
	seqcharge+=fDataSignals[fSignalPositionArray[0]];	
	seqaverage += fSignalPositionArray[0]*fDataSignals[fSignalPositionArray[0]];
	seqerror += fSignalPositionArray[0]*fSignalPositionArray[0]*fDataSignals[fSignalPositionArray[0]];
	tmpPos.push_back(fSignalPositionArray[0]);
	tmpSig.push_back(fDataSignals[fSignalPositionArray[0]]);
	  
	//Calculate mean of sequence:
	Int_t seqmean=0;
	seqmean = seqaverage/seqcharge;
	  
	//Calculate mean in pad direction:
	Int_t padmean = seqcharge*fPadNo;
	Int_t paderror = fPadNo*padmean;
	AliHLTTPCClusters candidate;
	candidate.fTotalCharge   = seqcharge;
	candidate.fPad       = padmean;
	candidate.fPad2      = paderror;
	candidate.fTime      = seqaverage;
	candidate.fTime2     = seqerror;
	candidate.fMean          = seqmean;
	candidate.fLastMergedPad = fPadNo;
	fClusterCandidates.push_back(candidate);
	fUsedClusterCandidates.push_back(0);
	isFalling=0;
	seqcharge=0;
	seqaverage=0;
	seqerror=0;

	tmpPos.clear();
	tmpSig.clear();
      }
    }
    else if(seqcharge>0){
      seqcharge+=fDataSignals[fSignalPositionArray[pos+1]];	
      seqaverage += fSignalPositionArray[pos+1]*fDataSignals[fSignalPositionArray[pos+1]];
      seqerror += fSignalPositionArray[pos+1]*fSignalPositionArray[pos+1]*fDataSignals[fSignalPositionArray[pos+1]];
      tmpPos.push_back(fSignalPositionArray[pos+1]);
      tmpSig.push_back(fDataSignals[fSignalPositionArray[pos+1]]);

      //Calculate mean of sequence:
      Int_t seqmean=0;
      seqmean = seqaverage/seqcharge;
	
      //Calculate mean in pad direction:
      Int_t padmean = seqcharge*fPadNo;
      Int_t paderror = fPadNo*padmean;
      AliHLTTPCClusters candidate;
      candidate.fTotalCharge   = seqcharge;
      candidate.fPad       = padmean;
      candidate.fPad2      = paderror;
      candidate.fTime      = seqaverage;
      candidate.fTime2     = seqerror;
      candidate.fMean          = seqmean;
      candidate.fLastMergedPad = fPadNo;
      fClusterCandidates.push_back(candidate);
      fUsedClusterCandidates.push_back(0);
      isFalling=0;
      seqcharge=0;
      seqaverage=0;
      seqerror=0;

      tmpPos.clear();
      tmpSig.clear();
    }
  }
  */
}

