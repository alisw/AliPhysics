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


///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Calibration class for DCS sensors                                        //
//  Authors: Marian Ivanov and Haavard Helstrup                              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliDCSSensorArray.h"
#include "AliLog.h"
#include <TMath.h>

ClassImp(AliDCSSensorArray)

const Double_t kSecInHour = 3600.; // seconds in one hour
const UInt_t   kMinMapTime = 60;   // don't fit maps shorter than one minute

//_____________________________________________________________________________
AliDCSSensorArray::AliDCSSensorArray():TNamed(), 
  fMinGraph(10),
  fMinPoints(10),
  fIter(10),
  fMaxDelta(0.0),
  fFitReq(2),
  fValCut(-1),
  fDiffCut(-1),
  fStartTime (2000,1,1,0,0,0),
  fEndTime   (2000,1,1,0,0,0),
  fSensors(0)
{
  //
  // AliDCSSensorArray default constructor
  //

}
//_____________________________________________________________________________
AliDCSSensorArray::AliDCSSensorArray(TClonesArray *arr):TNamed(),
  fMinGraph(10),
  fMinPoints(10),
  fIter(10),
  fMaxDelta(0.0),
  fFitReq(2),
  fValCut(-1),
  fDiffCut(-1),
  fStartTime (2000,1,1,0,0,0),
  fEndTime   (2000,1,1,0,0,0),
  fSensors(arr)
{
  //
  // AliDCSSensorArray special constructor taking TClonesArray from ReadList
  //

}
//_____________________________________________________________________________
AliDCSSensorArray::AliDCSSensorArray(Int_t run, const char* dbEntry) :
  TNamed(),
  fMinGraph(10),
  fMinPoints(10),
  fIter(10),
  fMaxDelta(0.0),
  fFitReq(2),
  fValCut(-1),
  fDiffCut(-1),
  fStartTime (2000,1,1,0,0,0),
  fEndTime   (2000,1,1,0,0,0),
  fSensors(0)
{
  //
  // Read configuration from OCDB
  //

  AliCDBEntry *entry = AliCDBManager::Instance()->Get(dbEntry,run);
  if (entry) {
    TTree *tree = (TTree*) entry->GetObject();
    fSensors = AliDCSSensor::ReadTree(tree);
  } else {
    AliError("Unable to load configuration from CDB!");
  }
}
//_____________________________________________________________________________
AliDCSSensorArray::AliDCSSensorArray(UInt_t startTime, UInt_t endTime,
                       TTree* confTree) :
  TNamed(),
  fMinGraph(10),
  fMinPoints(10),
  fIter(10),
  fMaxDelta(0.0),
  fFitReq(2),
  fValCut(-1),
  fDiffCut(-1),
  fStartTime (2000,1,1,0,0,0),
  fEndTime   (2000,1,1,0,0,0),
  fSensors(0)

{
  //
  // AliDCSSensorArray constructor for Shuttle preprocessor
  //  (confTree read from OCDB)
  //
  fSensors = AliDCSSensor::ReadTree(confTree);
  fSensors->BypassStreamer(kFALSE);
  fStartTime = TTimeStamp((time_t)startTime,0);
  fEndTime   = TTimeStamp((time_t)endTime,0);
}


//_____________________________________________________________________________
AliDCSSensorArray::AliDCSSensorArray(UInt_t startTime, UInt_t endTime,
                       TClonesArray *sensors) :
  TNamed(),
  fMinGraph(10),
  fMinPoints(10),
  fIter(10),
  fMaxDelta(0.0),
  fFitReq(2),
  fValCut(-1),
  fDiffCut(-1),
  fStartTime (2000,1,1,0,0,0),
  fEndTime   (2000,1,1,0,0,0),
  fSensors(sensors)

{
  //
  // AliDCSSensorArray constructor for Shuttle preprocessor
  //  (TClonesArray of AliDCSSensor objects)
  //
  fStartTime = TTimeStamp((time_t)startTime,0);
  fEndTime   = TTimeStamp((time_t)endTime,0);
}

//_____________________________________________________________________________
AliDCSSensorArray::AliDCSSensorArray(const AliDCSSensorArray &c):TNamed(c),
  fMinGraph(c.fMinGraph),
  fMinPoints(c.fMinPoints),
  fIter(c.fIter),
  fMaxDelta(c.fMaxDelta),
  fFitReq(c.fFitReq),
  fValCut(c.fValCut),
  fDiffCut(c.fDiffCut),
  fStartTime (c.fStartTime),
  fEndTime   (c.fEndTime),
  fSensors(0)

{
  //
  // AliDCSSensorArray copy constructor
  //

  fSensors = (TClonesArray*)c.fSensors->Clone();
}

///_____________________________________________________________________________
AliDCSSensorArray::~AliDCSSensorArray()
{
  //
  // AliDCSSensorArray destructor
  //
  fSensors->Delete();
  delete fSensors;

}

//_____________________________________________________________________________
AliDCSSensorArray &AliDCSSensorArray::operator=(const AliDCSSensorArray &c)
{
  //
  // Assignment operator
  //
  if (this != &c) {
     fSensors->Delete();
     new (this) AliDCSSensorArray(c);
     fSensors = (TClonesArray*)c.fSensors->Clone();
  }
  return *this;
}

//____________________________________________________________________________

void AliDCSSensorArray::SetGraph(TMap *map)
{
  //
  // Read graphs from DCS maps
  //
  Int_t nsensors = fSensors->GetEntries();
  for ( Int_t isensor=0; isensor<nsensors; isensor++) {
    AliDCSSensor *entry = (AliDCSSensor*)fSensors->At(isensor);
    TString stringID = entry->GetStringID();
    TGraph *gr = (TGraph*)map->GetValue(stringID.Data());
    if ( gr !=0 ) {
       entry->SetGraph((TGraph*)gr->Clone());
    } else {
       entry->SetGraph(0);
    }
  }
}
//_____________________________________________________________________________
void AliDCSSensorArray::MakeSplineFit(TMap *map, Bool_t keepMap)
{
  //
  // Make spline fits from DCS maps
  //
  Int_t nsensors = fSensors->GetEntries();
  for ( Int_t isensor=0; isensor<nsensors; isensor++) {
    AliDCSSensor *entry = (AliDCSSensor*)fSensors->At(isensor);
    TString stringID = entry->GetStringID();
    TGraph *gr = (TGraph*)map->GetValue(stringID.Data());
    if (!gr ) {
      entry->SetFit(0);
      entry->SetGraph(0);
      AliWarning(Form("sensor %s: no input graph",stringID.Data()));
      continue;
    }
    UInt_t timeDiff = entry->GetEndTime() - entry->GetStartTime();
    if ( timeDiff < kMinMapTime ) {
      AliWarning(Form("sensor %s: map length < 60 s, DCS graph kept.",stringID.Data()));
      entry->SetGraph((TGraph*)gr->Clone());
    } else {
      AliSplineFit *fit = new AliSplineFit();
      fit->SetMinPoints(fMinGraph);
      fit->InitKnots(gr,fMinPoints,fIter,fMaxDelta);
      fit->SplineFit(fFitReq);
      fit->Cleanup();
      if (fit->GetKnots()>0) {
        entry->SetFit(fit);
      } else {
        AliWarning(Form("sensor %s: no fit performed, DCS graph kept.",stringID.Data()));
        entry->SetGraph((TGraph*)gr->Clone());
      }
    }
    if (keepMap) entry->SetGraph((TGraph*)gr->Clone());
  }
}
//_____________________________________________________________________________
void AliDCSSensorArray::MakeSplineFitAddPoints(TMap *map)
{
  //
  // Make spline fits from DCS maps
  //
  Int_t nsensors = fSensors->GetEntries();
  for ( Int_t isensor=0; isensor<nsensors; isensor++) {
    AliDCSSensor *entry = (AliDCSSensor*)fSensors->At(isensor);

  // fetch old points from existing graph

    TGraph *gr = entry->GetGraph();
    if (!gr) {
      gr = new TGraph();
      entry->SetGraph(gr);
    } 
    TString stringID = entry->GetStringID();

  // fetch new points from DCS map
  
    TGraph *grAdd = (TGraph*)map->GetValue(stringID.Data());
    if (!grAdd ) return;

  // add new points to end of graph
  
    Int_t nPointsOld=gr->GetN();
    Int_t nPointsAdd=grAdd->GetN();
    gr->Expand(nPointsOld+nPointsAdd);
    gr->Set(nPointsOld+nPointsAdd);
    Double_t *addX=grAdd->GetX();
    Double_t *addY=grAdd->GetY();
    for (Int_t i=0;i<nPointsAdd;i++) {
      gr->SetPoint(nPointsOld+i,addX[i],addY[i]);
    }
 
   // make fit to complete graph
   
    AliSplineFit *fit = new AliSplineFit();
    fit->SetMinPoints(fMinGraph);
    fit->InitKnots(gr,fMinPoints,fIter,fMaxDelta);
    fit->SplineFit(fFitReq);
    fit->Cleanup();
    if (fit->GetKnots()>0) {
      AliSplineFit *oldFit = entry->GetFit();
      if (oldFit) delete oldFit;
      entry->SetFit(fit);
    } else {
      AliWarning(Form("sensor %s: no new fit performed. If available, old fit kept.",stringID.Data()));
    }
  }
}

//_____________________________________________________________________________
Int_t AliDCSSensorArray::NumFits() const 
{
 //
 // Return number of sensors where a succesful fit has been made
 //
  Int_t nfit=0;
  Int_t nsensors = fSensors->GetEntries();
  for ( Int_t isensor=0; isensor<nsensors; isensor++) {
    AliDCSSensor *entry = (AliDCSSensor*)fSensors->At(isensor);
    if (entry->GetFit()) nfit++;
  }    
  return nfit;
}
//_____________________________________________________________________________
Double_t AliDCSSensorArray::GetValue(UInt_t timeSec, Int_t sensor) 
{
  //
  // Return sensor value at time timeSec (obtained from fitted function)
  //  timeSec = time in seconds from start of run
  //

  AliDCSSensor *entry = (AliDCSSensor*)fSensors->At(sensor);
  return entry->GetValue(TTimeStamp((time_t)fStartTime.GetSec()+timeSec,0));
}


//_____________________________________________________________________________
TMap* AliDCSSensorArray::ExtractDCS(TMap *dcsMap, Bool_t keepStart)
{
 //
 // Extract temperature graphs from DCS maps
 //
 TMap *values = new TMap;
 TObjArray * valueSet;
 //
 // Keep global start/end times
 //    to avoid extrapolations, the fits will only be valid from first 
 //    measured point to last measured point. This is consistent with hardware,
 //    as there would be a new measured point if the value changed.
 
 TTimeStamp startTime=fStartTime;
 TTimeStamp endTime=fEndTime;
 
 Int_t nsensors = fSensors->GetEntries();
 for ( Int_t isensor=0; isensor<nsensors; isensor++) {
   AliDCSSensor *entry = (AliDCSSensor*)fSensors->At(isensor);
   TString stringID = entry->GetStringID();
   TPair *pair = (TPair*)dcsMap->FindObject(stringID.Data());
   if ( pair ) {                            // only try to read values
                                            // if DCS object available
     valueSet = (TObjArray*)pair->Value();
     TGraph *graph = MakeGraph(valueSet,keepStart);   // MakeGraph sets start/end time
                                            // per sensor
     values->Add(new TObjString(stringID.Data()),graph);
     entry->SetStartTime(fStartTime);
     entry->SetEndTime(fEndTime);
   }
 }
 // Reset global start/end time 
 //    ..... yes, I know this won't get a prize for structured programming..:-)

 fStartTime=startTime;
 fEndTime=endTime;
 return values;
}


//_____________________________________________________________________________
TGraph* AliDCSSensorArray::MakeGraph(TObjArray* valueSet, Bool_t keepStart){
  //
  // Make graph of temperature values read from DCS map
  //   (spline fit parameters will subsequently be obtained from this graph) 
  //
  Int_t nentries = valueSet->GetEntriesFast(); 
  if ( nentries == 0 ) return 0;
  
  Float_t *x = new Float_t[nentries];
  Float_t *y = new Float_t[nentries];
  Int_t time0=0, previousTime=0;
  TTimeStamp firstTime(0);
  TTimeStamp lastTime(0);
  if (keepStart) { 
     firstTime = fStartTime;
     time0 = firstTime.GetSec();
  }
  Int_t out=0;
  Int_t skipped=0;
  AliDCSValue *val = (AliDCSValue *)valueSet->At(0);
  AliDCSValue::Type type = val->GetType();
  if ( type == AliDCSValue::kInvalid || type == AliDCSValue::kBool ) {
     delete [] x;
     delete [] y;
     return 0;
  }
  Float_t value;
  for (Int_t i=0; i<nentries; i++){
    val = (AliDCSValue *)valueSet->At(i);
    if (!val) continue;
    if (time0==0){
      time0=val->GetTimeStamp();
      firstTime= TTimeStamp((time_t)val->GetTimeStamp(),0);
      lastTime=TTimeStamp((time_t)val->GetTimeStamp(),0);
     }
    switch ( type )
    { 
      case AliDCSValue::kFloat:
        value = val->GetFloat();
        break;
      case AliDCSValue::kChar:
        value = static_cast<Float_t>(val->GetChar());
	break;
      case AliDCSValue::kInt:
        value = static_cast<Float_t>(val->GetInt());
	break;
      case AliDCSValue::kUInt:
        value = static_cast<Float_t>(val->GetUInt());
	break;
      default:
        continue;
    }
    if (fValCut>0 && TMath::Abs(value)>fValCut) continue;   // refuse values greater than cut
    if (fDiffCut>0 ) {
      if ( out>0 && skipped<10 && TMath::Abs(value-y[out-1])>fDiffCut) {
        skipped++;                               // refuse values changing 
        continue;                                // by > cut  in one time step
      }                                          
      skipped=0;
    }					      
    if (val->GetTimeStamp()-time0>1000000) continue;
    if (val->GetTimeStamp()-previousTime < 1 ) continue;   // refuse duplicate recordings
    previousTime=val->GetTimeStamp();
    lastTime=TTimeStamp((time_t)val->GetTimeStamp(),0);
    x[out] = (val->GetTimeStamp()-time0)/kSecInHour; // give times in fractions of hours 
    y[out] = value;
    out++;    
  }
  if (!keepStart) fStartTime=firstTime;
  fEndTime=lastTime;
  TGraph * graph = new TGraph(out,x,y);
  delete [] x;
  delete [] y;
  return graph;
}

//_____________________________________________________________________________
void AliDCSSensorArray::RemoveGraphDuplicates(Double_t tolerance){
//
//   Remove points with same y value as the previous measured point
//   (to save space for non-fitted graphs -- i.e. last measured point used)
//
  Int_t nsensors = fSensors->GetEntries();
  for ( Int_t isensor=0; isensor<nsensors; isensor++) {
    AliDCSSensor *entry = (AliDCSSensor*)fSensors->At(isensor);
    TGraph *graph = entry->GetGraph();
    Double_t x=-999.,y=-999., x0=-999.,y0=-999.;
    if (graph) {
      Int_t npoints=graph->GetN();
      if (npoints>1) {
        for (Int_t i=npoints-1;i>0;i--) {
	   graph->GetPoint(i,x,y);
	   graph->GetPoint(i-1,x0,y0);
	   if ( TMath::Abs(y-y0) < TMath::Abs(tolerance*y0) ) graph->RemovePoint(i);
	 }
      }
    }
   }
}    

  
//_____________________________________________________________________________
AliDCSSensor* AliDCSSensorArray::GetSensor(Int_t IdDCS) 
{
 //
 //  Return sensor information for sensor specified by IdDCS
 //
 Int_t nsensors = fSensors->GetEntries();
 for (Int_t isensor=0; isensor<nsensors; isensor++) {
   AliDCSSensor *entry = (AliDCSSensor*)fSensors->At(isensor);
   if (entry->GetIdDCS() == IdDCS) return entry;
 }
 return 0;
}
//_____________________________________________________________________________
AliDCSSensor* AliDCSSensorArray::GetSensor(const TString& stringID)
{
 //
 //  Return sensor information for sensor specified by stringID
 //
 Int_t nsensors = fSensors->GetEntries();
 for (Int_t isensor=0; isensor<nsensors; isensor++) {
   AliDCSSensor *entry = (AliDCSSensor*)fSensors->At(isensor);
   if (entry->GetStringID() == stringID) return entry;
 }
 return 0;
}
//_____________________________________________________________________________
AliDCSSensor* AliDCSSensorArray::GetSensor(Double_t x, Double_t y, Double_t z)
{
 //
 //  Return sensor closest to given position
 //
 Int_t nsensors = fSensors->GetEntries();
 Double_t dist2min=1e99;
 Double_t xs,ys,zs,dist2;
 Int_t ind=-1;
 for (Int_t isensor=0; isensor<nsensors; isensor++) {
   AliDCSSensor *entry = (AliDCSSensor*)fSensors->At(isensor);
   xs = entry->GetX();
   ys = entry->GetY();
   zs = entry->GetZ();
   dist2 = (x-xs)*(x-xs) + (y-ys)*(y-ys) + (z-zs)*(z-zs);
   if (dist2 < dist2min) {
      ind=isensor;
      dist2min = dist2;
   }
 }
 if ( ind >= 0 ) {
    return (AliDCSSensor*)fSensors->At(ind);
 } else {
    return 0;
 }
}
//_____________________________________________________________________________
AliDCSSensor* AliDCSSensorArray::GetSensorNum(Int_t ind)
{
 //
 //  Return sensor given by array index
 //
 return (AliDCSSensor*)fSensors->At(ind);
}

//_____________________________________________________________________________
Int_t AliDCSSensorArray::SetSensor(const TString& stringID,
                          const  AliDCSSensor& sensor)
{
 //
 //  Update sensor information for sensor specified by stringID
 //
 Int_t nsensors = fSensors->GetEntries();
 for (Int_t isensor=0; isensor<nsensors; isensor++) {
   AliDCSSensor *entry = (AliDCSSensor*)fSensors->At(isensor);
   if (entry->GetStringID() == stringID) 
     {
      new ((*fSensors)[isensor])AliDCSSensor(sensor);
      return isensor;
     }
 }
 return -1;
}
//_____________________________________________________________________________
void AliDCSSensorArray::SetSensorNum(const Int_t ind, const AliDCSSensor& sensor)
{
 //
 //  Update sensor information for sensor at index ind
 //
   new ((*fSensors)[ind])AliDCSSensor(sensor);
   return;
}
//_____________________________________________________________________________
void AliDCSSensorArray::RemoveSensorNum(Int_t ind)
{
 //
 //  Return sensor given by array index
 //

  delete fSensors->RemoveAt(ind);
  fSensors->Compress();
}
//_____________________________________________________________________________
void AliDCSSensorArray::RemoveSensor(Int_t IdDCS)
{
 //
 //  Deletes Sensor by given IdDCS
 //

  Int_t nsensors = fSensors->GetEntries();
  for (Int_t isensor=0; isensor<nsensors; isensor++) { // loop over sensors
    AliDCSSensor *entry = (AliDCSSensor*)fSensors->At(isensor);
    if (entry->GetIdDCS()==IdDCS) {
      delete fSensors->RemoveAt(isensor);
      break;
    }
  }
  fSensors->Compress();
}
//_____________________________________________________________________________
TArrayI AliDCSSensorArray::OutsideThreshold(Double_t threshold, UInt_t timeSec, Bool_t below) const
{
 //
 // Return sensors with values outside threshold at time defined by second
 // parameter
 // By default sensors with values below threshold are listed, if third
 // parameter is set to kFALSE sensors with values above threshold are listed
 //
  Int_t nsensors = fSensors->GetEntries();
  TArrayI array(nsensors);
  Int_t outside=0;
  for (Int_t isensor=0; isensor<nsensors; isensor++) { // loop over sensors
    AliDCSSensor *entry = (AliDCSSensor*)fSensors->At(isensor);
    Double_t val=entry->GetValue(timeSec);
    if (below) {
      if (val<threshold) array[outside++] = entry->GetIdDCS();
    } else {
      if (val>threshold) array[outside++] = entry->GetIdDCS();
    }    
  }
  array.Set(outside);
  return array;
}
 
//_____________________________________________________________________________
Int_t AliDCSSensorArray::GetFirstIdDCS() const
{
 //
 //  Return DCS Id of first sensor
 //
 if ( fSensors != 0 ) {
    return ((AliDCSSensor*)fSensors->At(0))->GetIdDCS();
 } else {
    return 0;
 }
}

//_____________________________________________________________________________
Int_t AliDCSSensorArray::GetLastIdDCS() const 
{
 //
 //  Return DCS Id of last sensor
 //
 if ( fSensors != 0 ) {
    Int_t last = fSensors->GetEntries();
    return ((AliDCSSensor*)fSensors->At(last-1))->GetIdDCS();
 } else {
    return 0;
 }
}
//_____________________________________________________________________________
void AliDCSSensorArray::ClearGraph()
{
  //
  // Delete DCS graphs from all sensors in array
  //
   
   Int_t nsensors = fSensors->GetEntries();
   for ( Int_t isensor=0; isensor<nsensors; isensor++) {
     AliDCSSensor *sensor = (AliDCSSensor*)fSensors->At(isensor);
     TGraph *gr = sensor->GetGraph();
     if ( gr != 0 ) {
       delete gr;
       gr = 0;
     }
     sensor->SetGraph(0);
   }
}
//_____________________________________________________________________________
void AliDCSSensorArray::ClearFit()
{
  //
  // Delete spline fits from all sensors in array
  //

   Int_t nsensors = fSensors->GetEntries();
   for ( Int_t isensor=0; isensor<nsensors; isensor++) {
     AliDCSSensor *sensor = (AliDCSSensor*)fSensors->At(isensor);
     AliSplineFit *fit = sensor->GetFit();
     if ( fit != 0 ) {
       delete fit;
       fit = 0;
     }
     sensor->SetFit(0);
   }
}
//_____________________________________________________________________________
void AliDCSSensorArray::AddSensors(AliDCSSensorArray *newSensors)
{
  //
  // add sensors from two sensor arrays
  //
  
  Int_t numNew = newSensors->NumSensors();
  Int_t numOld = fSensors->GetEntries();
  fSensors->Expand(numOld+numNew);
  for (Int_t i=0;i<numNew;i++) {
    AliDCSSensor *sens = newSensors->GetSensorNum(i);
    new ((*fSensors)[numOld+i]) AliDCSSensor(*sens);
  }
}  
  
