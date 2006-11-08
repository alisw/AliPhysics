/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Authors: Matthias Richter <Matthias.Richter@ift.uib.no>                *
 *          for The ALICE Off-line Project.                               *
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
    @author Matthias Richter
    @date   
    @brief  Container Class for TPC Pads.
*/

#if __GNUC__>= 3
using namespace std;
#endif

#include <cerrno>
#include "AliHLTTPCPad.h"
#include "AliHLTStdIncludes.h"

/** margin for the base line be re-avaluated */
#define ALIHLTPAD_BASELINE_MARGIN (2*fAverage)

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTTPCPad)

AliHLTTPCPad::AliHLTTPCPad()
  :
  fRowNo(-1),
  fPadNo(-1),
  fThreshold(0),
  fAverage(-1),
  fNofEvents(0),
  fSum(0),
  fBLMax(-1),
  fBLMaxBin(-1),
  fBLMin(-1),
  fBLMinBin(-1),
  fCount(0),
  fTotal(0),
  fpRawData(NULL),
  fFirstBLBin(0),
  fNofBins(0),
  fReadPos(0)
{
}

AliHLTTPCPad::AliHLTTPCPad(Int_t offset, Int_t nofBins)
  :
  fRowNo(-1),
  fPadNo(-1),
  fThreshold(0),
  fAverage(-1),
  fNofEvents(0),
  fSum(0),
  fBLMax(-1),
  fBLMaxBin(-1),
  fBLMin(-1),
  fBLMinBin(-1),
  fCount(0),
  fTotal(0),
  fpRawData(NULL),
  fFirstBLBin(offset),
  fNofBins(nofBins),
  fReadPos(0)
{
}

AliHLTTPCPad::AliHLTTPCPad(const AliHLTTPCPad& srcPad)
  :
  fRowNo(srcPad.fRowNo),
  fPadNo(srcPad.fPadNo),
  fThreshold(0),
  fAverage(-1),
  fNofEvents(0),
  fSum(0),
  fBLMax(-1),
  fBLMaxBin(-1),
  fBLMin(-1),
  fBLMinBin(-1),
  fCount(0),
  fTotal(0),
  fpRawData(NULL),
  fFirstBLBin(0),
  fNofBins(0),
  fReadPos(0)
{
  HLTFatal("copy constructor not implemented");
}

AliHLTTPCPad& AliHLTTPCPad::operator=(const AliHLTTPCPad&)
{
  HLTFatal("assignment operator not implemented");
  return (*this);
}

AliHLTTPCPad::~AliHLTTPCPad()
{
  if (fpRawData) {
    HLTWarning("event data acquisition not stopped");
    StopEvent();
  }
}

Int_t AliHLTTPCPad::SetID(Int_t rowno, Int_t padno)
{
  fRowNo=rowno;
  fPadNo=padno;
}

Int_t AliHLTTPCPad::StartEvent()
{
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
	HLTDebug("maximum value %d exceeds margin for base line (%d) "
		 "-> re-evaluate base line", fBLMax, ALIHLTPAD_BASELINE_MARGIN);
	if (fpRawData) {
	  for (Int_t i=fFirstBLBin; i<fNofBins; i++)
	    if (fpRawData[i]>=0) AddBaseLineValue(i, fpRawData[i]);
	  if (fCount>0 && fCount>=reqMinCount && fCount>=fTotal/2) {
	    fAverage=fSum/fCount;
	    HLTDebug("new average %d", fAverage);
	  } else {
	    HLTDebug("baseline re-eveluation skipped because of to few "
		       "contributing bins: total=%d, contributing=%d, req=%d"
		       "\ndata might be already zero suppressed"
		       , fTotal, fCount, reqMinCount);
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
    HLTDebug("baseline calculation skipped because of to few contributing "
	       "bins: total=%d, contributing=%d, required=%d \ndata might be "
	       "already zero suppressed", fTotal, fCount, reqMinCount);
  }

  return iResult;
}

Int_t AliHLTTPCPad::StopEvent()
{
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
  Int_t iResult=0;
  fAverage=-1;
  fNofEvents=0;
  return iResult;
}

Int_t AliHLTTPCPad::SetThreshold(AliHLTTPCSignal_t thresh)
{
  Int_t iResult=0;
  fThreshold=thresh;
  return iResult;
}

Int_t AliHLTTPCPad::AddBaseLineValue(Int_t bin, AliHLTTPCSignal_t value)
{
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
      HLTDebug("ignoring value %d (bin %d) for base line calculation "
	       "(current average is %d)",
	       value, bin, fAverage);
    }
  }
  return iResult;
}

Int_t AliHLTTPCPad::SetRawData(Int_t bin, AliHLTTPCSignal_t value)
{
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
  fReadPos=(bZeroSuppression>0?0:fFirstBLBin)-1;
  return Next(bZeroSuppression);
}

AliHLTTPCSignal_t AliHLTTPCPad::GetRawData(Int_t bin) const
{
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

AliHLTTPCSignal_t AliHLTTPCPad::GetBaseLine(Int_t bin) const
{
  AliHLTTPCSignal_t val=0;
  if (fAverage>0) {
    // we take the minumum value as the base line if it doesn't differ from
    // the average to much
    const AliHLTTPCSignal_t kMaxDifference=15;
    val=fAverage;
#ifdef KEEP_NOISE
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
  return fAverage>0?fAverage:0;
}

Float_t AliHLTTPCPad::GetOccupancy() const
{
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
  // history is not yet implemented
  return GetOccupancy();
}
