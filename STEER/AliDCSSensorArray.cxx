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

ClassImp(AliDCSSensorArray)

const Double_t kSecInHour = 3600.; // seconds in one hour

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
  TTree *tree = (TTree*) entry->GetObject();
  fSensors = AliDCSSensor::ReadTree(tree);
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
  fStartTime = TTimeStamp(startTime);
  fEndTime   = TTimeStamp(endTime);
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
    AliSplineFit *fit = new AliSplineFit();
    fit->SetMinPoints(fMinGraph);
    fit->InitKnots(gr,fMinPoints,fIter,fMaxDelta);
    fit->SplineFit(fFitReq);
    entry->SetStartTime(fStartTime);
    entry->SetEndTime(fEndTime);
    fit->Cleanup();
    if (fit) {
      entry->SetFit(fit);
    } else {
      AliWarning(Form("sensor %s: no fit performed, DCS graph kept.",stringID.Data()));
      entry->SetGraph(gr);
    }
    if (keepMap) entry->SetGraph(gr);
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
  return entry->GetValue(timeSec);
}


//_____________________________________________________________________________
TMap* AliDCSSensorArray::ExtractDCS(TMap *dcsMap)
{
 //
 // Extract temperature graphs from DCS maps
 //
 TMap *values = new TMap;
 TObjArray * valueSet;
 Int_t nsensors = fSensors->GetEntries();
 for ( Int_t isensor=0; isensor<nsensors; isensor++) {
   AliDCSSensor *entry = (AliDCSSensor*)fSensors->At(isensor);
   TString stringID = entry->GetStringID();
   TPair *pair = (TPair*)dcsMap->FindObject(stringID.Data());
   if ( pair ) {                            // only try to read values
                                            // if DCS object available
     valueSet = (TObjArray*)pair->Value();
     TGraph *graph = MakeGraph(valueSet);
     values->Add(new TObjString(stringID.Data()),graph);
   }
 }
 return values;
}

//_____________________________________________________________________________
TGraph* AliDCSSensorArray::MakeGraph(TObjArray* valueSet){
  //
  // Make graph of temperature values read from DCS map
  //   (spline fit parameters will subsequently be obtained from this graph) 
  //
  Int_t nentries = valueSet->GetEntriesFast(); 
  if ( nentries == 0 ) return 0;
  
  Float_t *x = new Float_t[nentries];
  Float_t *y = new Float_t[nentries];
  Int_t time0=fStartTime.GetSec();
  Int_t out=0;
  Int_t skipped=0;
  AliDCSValue *val = (AliDCSValue *)valueSet->At(0);
  AliDCSValue::Type type = val->GetType();
  if ( type == AliDCSValue::kInvalid || type == AliDCSValue::kBool ) return 0;
  Float_t value;
  for (Int_t i=0; i<nentries; i++){
    val = (AliDCSValue *)valueSet->At(i);
    if (!val) continue;
    if (time0==0){
      time0=val->GetTimeStamp();
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
    x[out] = (val->GetTimeStamp()-time0)/kSecInHour; // give times in fractions of hours 
    y[out] = val->GetFloat();
    out++;
    
  }
  TGraph * graph = new TGraph(out,x,y);
  delete [] x;
  delete [] y;
  return graph;
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
 //  Return sensor information for sensor specified by IdDCS
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

AliDCSSensor* AliDCSSensorArray::GetSensorNum(Int_t ind)
{
 //
 //  Return sensor given by array index
 //
 return (AliDCSSensor*)fSensors->At(ind);
}

void AliDCSSensorArray::RemoveSensorNum(Int_t ind)
{
 //
 //  Return sensor given by array index
 //

  delete fSensors->RemoveAt(ind);
  fSensors->Compress();
}
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
