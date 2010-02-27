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
$Log: AliTOFLvHvDataPoints.cxx,v $
*/  

// AliTOFLvHvDataPoints class
// main aim to introduce the aliases for the TOF LV and HV DCS
// data points to be then
// stored in the OCDB, and to process them. 
// Process() method called by TOF preprocessor

#include "TString.h"
#include "TTimeStamp.h"
#include "TMap.h"
#include "TMath.h"

#include "AliDCSValue.h"
#include "AliLog.h"
#include "AliBitPacking.h"

#include "AliTOFLvHvDataPoints.h"
#include "AliTOFGeometry.h"

class AliCDBMetaData;
class TDatime;

ClassImp(AliTOFLvHvDataPoints)

//---------------------------------------------------------------
AliTOFLvHvDataPoints::AliTOFLvHvDataPoints():
  TObject(),
  fRun(0),
  fStartTime(0),
  fEndTime(0),
  fStartTimeDCSQuery(0),
  fEndTimeDCSQuery(0),
  fIsProcessed(kFALSE),
  fFDR(kFALSE),
  fNumberOfLVdataPoints(0),
  fNumberOfHVdataPoints(0),
  fNumberOfHVandLVmaps(0),
  fHisto(0x0)
{
  // main constructor 

}

//---------------------------------------------------------------
AliTOFLvHvDataPoints::AliTOFLvHvDataPoints(Int_t nRun, UInt_t startTime, UInt_t endTime, UInt_t startTimeDCSQuery, UInt_t endTimeDCSQuery):
  TObject(),
  fRun(nRun),
  fStartTime(startTime),
  fEndTime(endTime),
  fStartTimeDCSQuery(startTimeDCSQuery),
  fEndTimeDCSQuery(endTimeDCSQuery),
  fIsProcessed(kFALSE),
  fFDR(kFALSE),
  fNumberOfLVdataPoints(0),
  fNumberOfHVdataPoints(0),
  fNumberOfHVandLVmaps(0),
  fHisto(new TH1C("histo","",kNpads,-0.5,kNpads-0.5))
{

  // constructor with arguments

  AliInfo(Form("\n\tRun %d \n\tStartTime %s \n\tEndTime %s \n\tStartTime DCS Query %s \n\tEndTime DCS Query %s", nRun,
	       TTimeStamp(startTime).AsString(),
	       TTimeStamp(endTime).AsString(), 
	       TTimeStamp(startTimeDCSQuery).AsString(), 
	       TTimeStamp(endTimeDCSQuery).AsString()));

  Init();

}

//---------------------------------------------------------------

AliTOFLvHvDataPoints::AliTOFLvHvDataPoints(const AliTOFLvHvDataPoints & data):
  TObject(data), 
  fRun(data.fRun),
  fStartTime(data.fStartTime),
  fEndTime(data.fEndTime),
  fStartTimeDCSQuery(data.fStartTimeDCSQuery),
  fEndTimeDCSQuery(data.fEndTimeDCSQuery),
  fIsProcessed(data.fIsProcessed),
  fFDR(data.fFDR),
  fNumberOfLVdataPoints(data.fNumberOfLVdataPoints),
  fNumberOfHVdataPoints(data.fNumberOfHVdataPoints),
  fNumberOfHVandLVmaps(data.fNumberOfHVandLVmaps),
  fHisto(data.fHisto)
{

// copy constructor

  for(int i=0;i<kNddl;i++)
    fAliasNamesXLVmap[i]=data.fAliasNamesXLVmap[i];
    
  for(int i=0;i<kNsectors;i++)
    for(int j=0;j<kNplates;j++)
      fAliasNamesXHVmap[i][j]=data.fAliasNamesXHVmap[i][j];
 
}
//---------------------------------------------------------------

AliTOFLvHvDataPoints& AliTOFLvHvDataPoints:: operator=(const AliTOFLvHvDataPoints & data) { 

// assignment operator

  if (this == &data)
    return *this;

  TObject::operator=(data);
  fRun=data.GetRun();
  fStartTime=data.GetStartTime();
  fEndTime=data.GetEndTime();
  fStartTimeDCSQuery=data.GetStartTimeDCSQuery();
  fEndTimeDCSQuery=data.GetEndTimeDCSQuery();

  fNumberOfLVdataPoints=data.fNumberOfLVdataPoints;
  fNumberOfHVdataPoints=data.fNumberOfHVdataPoints;
  fNumberOfHVandLVmaps=data.fNumberOfHVandLVmaps;

  for(int i=0;i<kNddl;i++)
    fAliasNamesXLVmap[i]=data.fAliasNamesXLVmap[i];

  for(int i=0;i<kNsectors;i++)
    for(int j=0;j<kNplates;j++)
      fAliasNamesXHVmap[i][j]=data.fAliasNamesXHVmap[i][j];

  fHisto=data.fHisto;

  return *this;
}
//---------------------------------------------------------------
AliTOFLvHvDataPoints::~AliTOFLvHvDataPoints() {

  // destructor

}

//---------------------------------------------------------------
Bool_t AliTOFLvHvDataPoints::ProcessData(TMap& aliasMap) {
  //
  // method to process the data
  //

  if(!(fAliasNamesXHVmap[0][0]) || !(fAliasNamesXLVmap[0])) Init();

  AliInfo(Form(" Start Time = %i",fStartTime));
  AliInfo(Form(" End Time = %i",fEndTime));
  AliInfo(Form(" Start Time DCS Query= %i",fStartTimeDCSQuery));
  AliInfo(Form(" End Time DCS Query= %i",fEndTimeDCSQuery));

  if (fEndTime==fStartTime){
    AliError(Form(" Run with null time length: start time = %d = end time = %d",fStartTime,fEndTime));
    return kFALSE;
  }


  if (!ReadHVDataPoints(aliasMap)) return kFALSE;
  AliDebug(1,Form(" Number of HV dps = %d",fNumberOfHVdataPoints));
  if (!ReadLVDataPoints(aliasMap)) return kFALSE;
  AliDebug(1,Form(" Number of LV dps = %d",fNumberOfLVdataPoints));

  if (!MergeHVmap()) return kFALSE;
  if (!MergeLVmap()) return kFALSE;

  if (!MergeMaps()) return kFALSE;

  fIsProcessed=kTRUE;

  return kTRUE;

}

//---------------------------------------------------------------
Bool_t AliTOFLvHvDataPoints::MergeMaps() {
  //
  // Merge together LV and HV maps
  //

  Int_t counter=0;


  for (Int_t ii=1; ii<fNumberOfHVdataPoints-1; ii++)
    for (Int_t jj=1; jj<fNumberOfLVdataPoints-1; jj++) {
      if (fHVDataPoints[ii]->GetTime()==fLVDataPoints[jj]->GetTime())
	counter++;
    }


  fNumberOfHVandLVmaps = fNumberOfHVdataPoints-2+fNumberOfLVdataPoints-2-counter+2;
  AliInfo(Form(" Number of TOF HVandLV dps in the current run = %d",fNumberOfHVandLVmaps));

  counter=-1;
  Int_t timeMaps[kNmaxDataPoints];
  for (Int_t ii=0; ii<kNmaxDataPoints; ii++) timeMaps[ii]=0;

  for (Int_t ii=0; ii<fNumberOfHVdataPoints; ii++)
    timeMaps[ii]=fHVDataPoints[ii]->GetTime();

  Bool_t check = kTRUE;
  for (Int_t jj=0; jj<fNumberOfLVdataPoints; jj++) {
    check = kTRUE;
    for (Int_t ii=0; ii<fNumberOfHVdataPoints; ii++)
      check=check&&(fLVDataPoints[jj]->GetTime()!=timeMaps[ii]);

    if (check) {
      counter++;
      timeMaps[fNumberOfHVdataPoints+counter]=fLVDataPoints[jj]->GetTime();
    }
  }
  if (fNumberOfHVdataPoints+counter+1!=fNumberOfHVandLVmaps) AliWarning("Something is wrong!");

  Int_t controller[kNmaxDataPoints];
  for (Int_t ii=0; ii<kNmaxDataPoints; ii++) controller[ii]=-1;
  TMath::Sort(fNumberOfHVandLVmaps,timeMaps,controller,kFALSE); // increasing order


  Short_t array[kNpads];
  for (Int_t iPad=0; iPad<kNpads; iPad++) array[iPad]=-1;
  Int_t time = 0;

  // HVandLV status map @ SOR
  for (Int_t iPad=0; iPad<kNpads; iPad++)
    array[iPad] = fHVDataPoints[0]->GetCellValue(iPad)*fLVDataPoints[0]->GetCellValue(iPad);
  time = timeMaps[controller[0]];
  AliTOFDCSmaps *object0 = new AliTOFDCSmaps(time,array);
  fMap[0]= object0;

  // HVandLV status map during run
  for (Int_t index=1; index<fNumberOfHVandLVmaps-1; index++) {
    time = timeMaps[controller[index]];
    AliDebug(2,Form(" controller=%d time=%d ", controller[index],time));

    for (Int_t ii=0; ii<fNumberOfHVdataPoints-1; ii++)
      for (Int_t jj=0; jj<fNumberOfLVdataPoints-1; jj++) {

	if ( ( ( fHVDataPoints[ii]->GetTime()==time && fHVDataPoints[ii+1]->GetTime()>time ) &&
	       ( fLVDataPoints[jj]->GetTime()<=time && fLVDataPoints[jj+1]->GetTime()>time ) ) ||
	     ( ( fLVDataPoints[jj]->GetTime()==time && fLVDataPoints[jj+1]->GetTime()>time ) &&
	       ( fHVDataPoints[ii]->GetTime()<=time && fHVDataPoints[ii+1]->GetTime()>time ) ) ) {

	  AliDebug(2,Form(" HVdp_time=%d, LVdp_time=%d ",fHVDataPoints[ii]->GetTime(), fLVDataPoints[jj]->GetTime()));

	  for (Int_t iPad=0; iPad<kNpads; iPad++)
	    array[iPad] = fHVDataPoints[ii]->GetCellValue(iPad)*fLVDataPoints[jj]->GetCellValue(iPad);
	  AliTOFDCSmaps *object = new AliTOFDCSmaps(time,array);
	  fMap[index]= object;
	  break;
	}
      }

  }

  // HVandLV status map @ EOR
  for (Int_t iPad=0; iPad<kNpads; iPad++)
    array[iPad] = fHVDataPoints[fNumberOfHVdataPoints-1]->GetCellValue(iPad)*fLVDataPoints[fNumberOfLVdataPoints-1]->GetCellValue(iPad);
  time = timeMaps[controller[fNumberOfHVandLVmaps-1]];
  AliTOFDCSmaps *object1 = new AliTOFDCSmaps(time,array);
  fMap[fNumberOfHVandLVmaps-1]= object1;


  return kTRUE;

}

//---------------------------------------------------------------
Bool_t AliTOFLvHvDataPoints::MergeHVmap() {
  //
  // Create HV maps from HV dps
  //

  Bool_t check= kFALSE;

  for (Int_t iPad=0; iPad<kNpads; iPad+=96*91)
    AliDebug(2,Form(" HVdp0: channel=%6d -> %1d",iPad,fHVDataPoints[0]->GetCellValue(iPad)));

  if (fNumberOfHVdataPoints==2)
    check=kTRUE;
  else {
    for (Int_t ii=1; ii<fNumberOfHVdataPoints-1; ii++) {

      for (Int_t iPad=0; iPad<kNpads; iPad++) {
	if (fHVDataPoints[ii]->GetCellValue(iPad)==-1)
	  fHVDataPoints[ii]->SetCellValue(iPad,fHVDataPoints[ii-1]->GetCellValue(iPad));
      }
      for (Int_t iPad=0; iPad<kNpads; iPad+=96*91)
	AliDebug(2,Form("HVdp%d: channel=%6d -> %1d",ii,iPad,fHVDataPoints[ii]->GetCellValue(iPad)));
    }

    check=kTRUE;
  }

  return kTRUE;

}

//---------------------------------------------------------------
Bool_t AliTOFLvHvDataPoints::MergeLVmap() {
  //
  // Create LV maps from LV dps
  //

  Bool_t check= kFALSE;

  for (Int_t iPad=0; iPad<kNpads; iPad+=96*91)
    AliDebug(2,Form(" LVdp0: channel=%6d -> %1d",iPad,fLVDataPoints[0]->GetCellValue(iPad)));

  if (fNumberOfLVdataPoints==2)
    check=kTRUE;
  else {
    for (Int_t ii=1; ii<fNumberOfLVdataPoints-1; ii++) {

      for (Int_t iPad=0; iPad<kNpads; iPad++) {
	if (fLVDataPoints[ii]->GetCellValue(iPad)==-1)
	  fLVDataPoints[ii]->SetCellValue(iPad,fLVDataPoints[ii-1]->GetCellValue(iPad));
      }
      for (Int_t iPad=0; iPad<kNpads; iPad+=96*91)
	AliDebug(2,Form("LVpd%d: channel=%6d -> %1d",ii,iPad,fLVDataPoints[ii]->GetCellValue(iPad)));
    }

    check=kTRUE;
  }

  return check;

}

//---------------------------------------------------------------
Bool_t AliTOFLvHvDataPoints::ReadHVDataPoints(TMap& aliasMap) {
  //
  // Read HV dps
  //

  TObjArray *aliasArr;
  AliDCSValue* aValue;
  AliDCSValue* aValuePrev;
  UInt_t val = 0;
  Int_t time = 0;
  Int_t nEntries = 0;

  Short_t dummy[kNpads];

  for (Int_t iBin=0; iBin<kNpads; iBin++) dummy[iBin]=-1;
  // starting loop on aliases
  for (int i=0; i<kNsectors; i++)
    for (int j=0; j<kNplates; j++) {
      aliasArr = (TObjArray*) aliasMap.GetValue(fAliasNamesXHVmap[i][j].Data());
      if (!aliasArr) {
	AliError(Form("Alias %s not found!", fAliasNamesXHVmap[i][j].Data()));
	if (!fFDR)
	  return kFALSE;    // returning only in case we are not in a FDR run
	else
	  continue;
      }

      nEntries = aliasArr->GetEntries();
    
      if (nEntries==0) {
	AliError(Form("Alias %s has no entries! Nothing will be stored",
		      fAliasNamesXHVmap[i][j].Data()));
	continue;
      }
      else {

	for (Int_t iEntry=0; iEntry<nEntries; iEntry++) {
	  for (Int_t iBin=0; iBin<kNpads; iBin++) dummy[iBin]=-1;
	  aValue = (AliDCSValue*) aliasArr->At(iEntry);
	  val = aValue->GetUInt();
	  time = aValue->GetTimeStamp();
	  if (iEntry==0 || iEntry==nEntries-1) {
	    FillHVarrayPerDataPoint(i,j,val,dummy);
	    AliTOFDCSmaps *object = new AliTOFDCSmaps(time,dummy);
	    InsertHVDataPoint(object);
	  }
	  else {
	    aValuePrev = (AliDCSValue*) aliasArr->At(iEntry-1);
	    if (aValuePrev->GetUInt()!=val) {
	      FillHVarrayPerDataPoint(i,j,val,dummy);
	      AliTOFDCSmaps *object = new AliTOFDCSmaps(time,dummy);
	      InsertHVDataPoint(object);
	    }
	  }

	}
      }
    }

  return kTRUE;

}

//---------------------------------------------------------------
Bool_t AliTOFLvHvDataPoints::ReadLVDataPoints(TMap& aliasMap) {
  //
  // Read LV dps
  //

  TObjArray *aliasArr;
  AliDCSValue* aValue;
  AliDCSValue* aValuePrev;
  UInt_t val = 0;
  Int_t time = 0;
  Int_t nEntries = 0;

  Short_t dummy[kNpads];

  for (Int_t iBin=0; iBin<kNpads; iBin++) dummy[iBin]=-1;
  // starting loop on aliases
  for (int i=0; i<kNddl; i++) {
    aliasArr = (TObjArray*) aliasMap.GetValue(fAliasNamesXLVmap[i].Data());
    if (!aliasArr) {
      AliError(Form("Alias %s not found!", fAliasNamesXLVmap[i].Data()));
      if (!fFDR)
	return kFALSE;    // returning only in case we are not in a FDR run
      else
	continue;
    }

    nEntries = aliasArr->GetEntries();
    
    if (nEntries==0) {
      AliError(Form("Alias %s has no entries! Nothing will be stored",
		    fAliasNamesXLVmap[i].Data()));
      continue;
    }
    else {

      for (Int_t iEntry=0; iEntry<nEntries; iEntry++) {
	for (Int_t iBin=0; iBin<kNpads; iBin++) dummy[iBin]=-1;
	aValue = (AliDCSValue*) aliasArr->At(iEntry);
	val = aValue->GetUInt();
	time = aValue->GetTimeStamp();
	if (iEntry==0 || iEntry==nEntries-1) {
	  FillLVarrayPerDataPoint(i,val,dummy);
	  AliTOFDCSmaps *object = new AliTOFDCSmaps(time,dummy);
	  InsertLVDataPoint(object);
	}
	else {
	  aValuePrev = (AliDCSValue*) aliasArr->At(iEntry-1);
	  if (aValuePrev->GetUInt()!=val) {
	    FillLVarrayPerDataPoint(i,val,dummy);
	    AliTOFDCSmaps *object = new AliTOFDCSmaps(time,dummy);
	    InsertLVDataPoint(object);
	  }
	}
      }
    }
  }

  return kTRUE;

}

//---------------------------------------------------------------
Int_t AliTOFLvHvDataPoints::InsertHVDataPoint(AliTOFDCSmaps *object)
{
  //
  // Insert HV dp in the HV dps array.
  // The HV dps array is sorted according to increasing dp timeStamp value
  //

  if (fNumberOfHVdataPoints==kNmaxDataPoints) {
    AliError("Too many HV data points!");
    return 1;
  }

  if (fNumberOfHVdataPoints==0) {
    fHVDataPoints[fNumberOfHVdataPoints++] = object;
    return 0;
  }

  for (Int_t index=0; index<fNumberOfHVdataPoints; index++) {
    if (object->GetTime()==fHVDataPoints[index]->GetTime()) {
      fHVDataPoints[index]->Update(object);
      return 0;
    }
  }

  Int_t ii = FindHVdpIndex(object->GetTime());
  memmove(fHVDataPoints+ii+1 ,fHVDataPoints+ii,(fNumberOfHVdataPoints-ii)*sizeof(AliTOFDCSmaps*));
  fHVDataPoints[ii] = object;
  fNumberOfHVdataPoints++;
  
  return 0;

}

//_________________________________________________________________________
Int_t AliTOFLvHvDataPoints::FindHVdpIndex(Int_t z) const {
  //
  // This function returns the index of the nearest HV DP in time
  //

  if (fNumberOfHVdataPoints==0) return 0;
  if (z <= fHVDataPoints[0]->GetTime()) return 0;
  if (z > fHVDataPoints[fNumberOfHVdataPoints-1]->GetTime()) return fNumberOfHVdataPoints;
  Int_t b = 0, e = fNumberOfHVdataPoints-1, m = (b+e)/2;
  for (; b<e; m=(b+e)/2) {
    if (z > fHVDataPoints[m]->GetTime()) b=m+1;
    else e=m;
  }

  return m;

}

//---------------------------------------------------------------
Int_t AliTOFLvHvDataPoints::InsertLVDataPoint(AliTOFDCSmaps *object)
{
  //
  // Insert LV dp in the LV dps array.
  // The LV dps array is sorted according to increasing dp timeStamp value
  //

  if (fNumberOfLVdataPoints==kNmaxDataPoints) {
    AliError("Too many LV data points!");
    return 1;
  }

  if (fNumberOfLVdataPoints==0) {
    fLVDataPoints[fNumberOfLVdataPoints++] = object;
    return 0;
  }

  for (Int_t index=0; index<fNumberOfLVdataPoints; index++) {
    if (object->GetTime()==fLVDataPoints[index]->GetTime()) {
      fLVDataPoints[index]->Update(object);
      return 0;
    }
  }

  Int_t ii = FindLVdpIndex(object->GetTime());
  memmove(fLVDataPoints+ii+1 ,fLVDataPoints+ii,(fNumberOfLVdataPoints-ii)*sizeof(AliTOFDCSmaps*));
  fLVDataPoints[ii] = object;
  fNumberOfLVdataPoints++;
  
  return 0;

}

//_________________________________________________________________________
Int_t AliTOFLvHvDataPoints::FindLVdpIndex(Int_t z) const {
  //
  // This function returns the index of the nearest LV DP in time
  //

  if (fNumberOfLVdataPoints==0) return 0;
  if (z <= fLVDataPoints[0]->GetTime()) return 0;
  if (z > fLVDataPoints[fNumberOfLVdataPoints-1]->GetTime()) return fNumberOfLVdataPoints;
  Int_t b = 0, e = fNumberOfLVdataPoints-1, m = (b+e)/2;
  for (; b<e; m=(b+e)/2) {
    if (z > fLVDataPoints[m]->GetTime()) b=m+1;
    else e=m;
  }

  return m;

}

//---------------------------------------------------------------
void AliTOFLvHvDataPoints::Init(){
  //
  // Initialize aliases and DCS data
  //

  TString sindex;
  for(int i=0;i<kNsectors;i++)
    for(int j=0;j<kNplates;j++) {
      fAliasNamesXHVmap[i][j] = "TOF_HVSTATUS_";
      sindex.Form("SM%02dMOD%1d",i,j);
      fAliasNamesXHVmap[i][j] += sindex;
    }


  for(int i=0;i<kNddl;i++) {
    fAliasNamesXLVmap[i] = "TOF_FEACSTATUS_";
    sindex.Form("%02d",i);
    fAliasNamesXLVmap[i] += sindex;
  }

}

//---------------------------------------------------------------
void AliTOFLvHvDataPoints::FillHVarrayPerDataPoint(Int_t sector, Int_t plate, UInt_t baseWord, Short_t *array)
{
  //
  // Set the status of the TOF pads connected to the HV dp
  // labelled by sector and plate numbers
  //

  Int_t det[5] = {sector, plate, -1, -1, -1};
  UInt_t checkBit = 0;

  Int_t channel = -1;

  for (Int_t iStrip=0; iStrip<AliTOFGeometry::NStrip(plate); iStrip++) {
    checkBit = AliBitPacking::UnpackWord(baseWord,iStrip,iStrip);
    for (Int_t iPadZ=0; iPadZ<AliTOFGeometry::NpadZ(); iPadZ++)
      for (Int_t iPadX=0; iPadX<AliTOFGeometry::NpadX(); iPadX++) {
	det[2] = iStrip;
	det[3] = iPadZ;
	det[4] = iPadX;
	channel = AliTOFGeometry::GetIndex(det);
	array[channel]=checkBit;
    }
  }


}

//---------------------------------------------------------------
void AliTOFLvHvDataPoints::FillLVarrayPerDataPoint(Int_t nDDL, UInt_t baseWord, Short_t *array)
{
  //
  // Set the status of the TOF pads connected to the LV dp
  // labelled by TOF crate number
  //

  Int_t det[5] = {nDDL/4, -1, -1, -1, -1};
  UInt_t checkBit = 0;

  Int_t iStripXsm[6] = {-1,-1,-1,-1,-1,-1};
  Int_t firstPadX = -1;
  Int_t lastPadX = -1;
  Int_t plate = -1;
  Int_t strip = -1;

  Int_t channel = -1;

  for (Int_t nFEAC=0; nFEAC<8; nFEAC++) {
    checkBit = AliBitPacking::UnpackWord(baseWord,nFEAC,nFEAC);
    firstPadX = -1;
    lastPadX = -1;
    GetStripsConnectedToFEAC(nDDL, nFEAC, iStripXsm, firstPadX,lastPadX);
    for (Int_t index=0; index<6; index++) {
      if (iStripXsm[index]==-1) continue;

      for (Int_t iPadZ=0; iPadZ<AliTOFGeometry::NpadZ(); iPadZ++)
	for (Int_t iPadX=firstPadX; iPadX<=lastPadX; iPadX++) {
	  AliTOFGeometry::GetStripAndModule(iStripXsm[index],plate,strip);
	  det[1] = plate;
	  det[2] = strip;
	  det[3] = iPadZ;
	  det[4] = iPadX;
	  channel = AliTOFGeometry::GetIndex(det);
	  array[channel]=checkBit;
	}
    }
  }


}

//---------------------------------------------------------------
void AliTOFLvHvDataPoints::GetStripsConnectedToFEAC(Int_t nDDL, Int_t nFEAC, Int_t *iStrip, Int_t &firstPadX, Int_t &lastPadX)
{
  //
  // FEAC-strip mapping:
  // return the strips and first PadX numbers
  // connected to the FEAC number nFEAC in the crate number nDDL
  //

  for (Int_t ii=0; ii<6; ii++) iStrip[ii]=-1;

  if (nDDL<0 || nDDL>=kNddl || nFEAC<0 || nFEAC>=8) return;

  switch (nDDL%4) {
  case 0:
    firstPadX =  0;
    lastPadX  = AliTOFGeometry::NpadX()/2-1;

    if (nFEAC<=2)
      for (Int_t ii=0; ii<6; ii++) iStrip[ii]=ii+6*nFEAC;
    else if (nFEAC==3)
      for (Int_t ii=0; ii<5; ii++) iStrip[ii]=ii+6*nFEAC;
    else if (nFEAC==4)
      for (Int_t ii=0; ii<6; ii++) iStrip[ii]=ii+6*nFEAC-1;
    else if (nFEAC==5)
      for (Int_t ii=0; ii<5; ii++) iStrip[ii]=ii+6*nFEAC-1;
    else if (nFEAC==6)
      for (Int_t ii=0; ii<6; ii++) iStrip[ii]=ii+6*nFEAC-2;
    else if (nFEAC==7)
      for (Int_t ii=0; ii<5; ii++) iStrip[ii]=ii+6*nFEAC-2;

    break; 
  case 1:
    firstPadX = AliTOFGeometry::NpadX()/2;
    lastPadX  = AliTOFGeometry::NpadX()-1;

    if (nFEAC<=2)
      for (Int_t ii=0; ii<6; ii++) iStrip[ii]=ii+6*nFEAC;
    else if (nFEAC==3)
      for (Int_t ii=0; ii<5; ii++) iStrip[ii]=ii+6*nFEAC;
    else if (nFEAC==4)
      for (Int_t ii=0; ii<6; ii++) iStrip[ii]=ii+6*nFEAC-1;
    else if (nFEAC==5)
      for (Int_t ii=0; ii<6; ii++) iStrip[ii]=ii+6*nFEAC-1;
    else if (nFEAC==6)
      for (Int_t ii=0; ii<5; ii++) iStrip[ii]=ii+6*nFEAC-1;
    else if (nFEAC==7)
      for (Int_t ii=0; ii<6; ii++) iStrip[ii]=ii+6*nFEAC-2;

    break;
  case 2:
    firstPadX = AliTOFGeometry::NpadX()/2;
    lastPadX  = AliTOFGeometry::NpadX()-1;

    if (nFEAC<=2)
      for (Int_t ii=0; ii<6; ii++) iStrip[ii]=90-(ii+6*nFEAC);
    else if (nFEAC==3)
      for (Int_t ii=0; ii<5; ii++) iStrip[ii]=90-(ii+6*nFEAC);
    else if (nFEAC==4)
      for (Int_t ii=0; ii<6; ii++) iStrip[ii]=90-(ii+6*nFEAC-1);
    else if (nFEAC==5)
      for (Int_t ii=0; ii<5; ii++) iStrip[ii]=90-(ii+6*nFEAC-1);
    else if (nFEAC==6)
      for (Int_t ii=0; ii<6; ii++) iStrip[ii]=90-(ii+6*nFEAC-2);
    else if (nFEAC==7)
      for (Int_t ii=0; ii<5; ii++) iStrip[ii]=90-(ii+6*nFEAC-2);

    break;
  case 3:
    firstPadX =  0;
    lastPadX  = AliTOFGeometry::NpadX()/2-1;

    if (nFEAC<=2)
      for (Int_t ii=0; ii<6; ii++) iStrip[ii]=90-(ii+6*nFEAC);
    else if (nFEAC==3)
      for (Int_t ii=0; ii<5; ii++) iStrip[ii]=90-(ii+6*nFEAC);
    else if (nFEAC==4)
      for (Int_t ii=0; ii<6; ii++) iStrip[ii]=90-(ii+6*nFEAC-1);
    else if (nFEAC==5)
      for (Int_t ii=0; ii<6; ii++) iStrip[ii]=90-(ii+6*nFEAC-1);
    else if (nFEAC==6)
      for (Int_t ii=0; ii<5; ii++) iStrip[ii]=90-(ii+6*nFEAC-1);
    else if (nFEAC==7)
      for (Int_t ii=0; ii<6; ii++) iStrip[ii]=90-(ii+6*nFEAC-2);

    break;
 }

}

//---------------------------------------------------------------
void AliTOFLvHvDataPoints::Draw(const Option_t* /*option*/)
{
  //
  // Draw all histos and graphs
  //

  if(!fIsProcessed) return;
  /*
  TCanvas *ch;
  TString canvasHistoName="Histos";
  ch = new TCanvas(canvasHistoName,canvasHistoName,20,20,600,600);
  ch->cd();
  */
  // to be implemented

}


//---------------------------------------------------------------
void AliTOFLvHvDataPoints::DrawHVandLVMap(Int_t index)
{
  //
  // Draw HV+LV map labelled as index
  //

  if(!fIsProcessed) return;

  if (index>=fNumberOfHVandLVmaps) return;

  AliTOFDCSmaps *mappa=(AliTOFDCSmaps*)GetHVandLVmap(index);

  char title[100];
  if (index==0) sprintf(title,"HVandLV map at time %d (%dst map)",mappa->GetTime(),index+1);
  else if (index==1) sprintf(title,"HVandLV map at time %d (%dnd map)",mappa->GetTime(),index+1);
  else if (index==2) sprintf(title,"HVandLV map at time %d (%drd map)",mappa->GetTime(),index+1);
  else if (index>=3) sprintf(title,"HVandLV map at time %d (%dth map)",mappa->GetTime(),index+1);
  fHisto->Delete();
  fHisto = new TH1C("histo","",kNpads,-0.5,kNpads-0.5);
  //fHisto->Clear();
  fHisto->SetTitle(title);

  for (Int_t ii=0; ii<kNpads; ii++)
    fHisto->SetBinContent(ii+1,mappa->GetCellValue(ii));

  fHisto->Draw();

}

//---------------------------------------------------------------
void AliTOFLvHvDataPoints::DrawLVMap(Int_t index)
{
  //
  // Draw LV map labelled as index
  //

  if(!fIsProcessed) return;

  if (index>=fNumberOfLVdataPoints) return;

  AliTOFDCSmaps *mappa=(AliTOFDCSmaps*)GetLVmap(index);

  char title[100];
  if (index==0) sprintf(title,"LV map at time %d (%dst map)",mappa->GetTime(),index+1);
  else if (index==1) sprintf(title,"LV map at time %d (%dnd map)",mappa->GetTime(),index+1);
  else if (index==2) sprintf(title,"LV map at time %d (%drd map)",mappa->GetTime(),index+1);
  else if (index>=3) sprintf(title,"LV map at time %d (%dth map)",mappa->GetTime(),index+1);
  fHisto->Delete();
  fHisto = new TH1C("histo","",kNpads,-0.5,kNpads-0.5);
  //fHisto->Clear();
  fHisto->SetTitle(title);

  for (Int_t ii=0; ii<kNpads; ii++)
    fHisto->SetBinContent(ii+1,mappa->GetCellValue(ii));

  fHisto->Draw();

}

//---------------------------------------------------------------
void AliTOFLvHvDataPoints::DrawHVMap(Int_t index)
{
  //
  // Draw HV map labelled as index
  //

  if(!fIsProcessed) return;

  if (index>=fNumberOfHVdataPoints) return;

  AliTOFDCSmaps *mappa=(AliTOFDCSmaps*)GetHVmap(index);

  char title[100];
  if (index==0) sprintf(title,"HV map at time %d (%dst map)",mappa->GetTime(),index+1);
  else if (index==1) sprintf(title,"HV map at time %d (%dnd map)",mappa->GetTime(),index+1);
  else if (index==2) sprintf(title,"HV map at time %d (%drd map)",mappa->GetTime(),index+1);
  else if (index>=3) sprintf(title,"HV map at time %d (%dth map)",mappa->GetTime(),index+1);
  fHisto->Delete();
  fHisto = new TH1C("histo","",kNpads,-0.5,kNpads-0.5);
  //fHisto->Clear();
  fHisto->SetTitle(title);

  for (Int_t ii=0; ii<kNpads; ii++)
    fHisto->SetBinContent(ii+1,mappa->GetCellValue(ii));

  fHisto->Draw();

}

//---------------------------------------------------------------
AliTOFDCSmaps *AliTOFLvHvDataPoints::GetHVandLVmapAtEOR()
{
  //
  // Returns HVandLV status map at EOR.
  // If the end-of-run has been caused by TOF self,
  // the last but two value of HVandLV status map
  // will be taken into account.
  // This last condition is true
  // if the time interval between the second-last DP and the last one
  // is less than 60s.
  //

  AliTOFDCSmaps * lvANDhvMap = 0;

  if (fNumberOfHVandLVmaps==2) { // nothing changed during the run

    lvANDhvMap = fMap[fNumberOfHVandLVmaps-1];

  }
  else {

    if (fMap[fNumberOfHVandLVmaps-1]->GetTime()-fMap[fNumberOfHVandLVmaps-2]->GetTime()<=60)
      lvANDhvMap = (AliTOFDCSmaps*)fMap[fNumberOfHVandLVmaps-3];
    else
      lvANDhvMap = (AliTOFDCSmaps*)fMap[fNumberOfHVandLVmaps-1];

  }

  return lvANDhvMap;

}
