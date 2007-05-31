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

/* $Id$ */

////////////////////////////////////////////////////////////////////////////
//                                                                        //
// This class collects the DCS information and does an AliSplineFit       //
//                                                                        //
// Author:                                                                //
//   W. Monange   (wilfried.monange@free.fr)                              //
//                                                                        //
////////////////////////////////////////////////////////////////////////////


#include <TGraph.h>
#include <TObjArray.h>
#include <AliCDBMetaData.h>
#include <TMap.h>

#include "AliDCSValue.h"
#include "AliDCSValue.h"
#include "AliLog.h"
#include "AliTRDDataDCS.h"
#include "AliSplineFit.h"

ClassImp(AliTRDDataDCS)

/*
Typical use :

AliTRDDataDCS dataDcs;
dataDcs.ExtractDCS (dcs);
->call a storeRef function
dataDcs.PerformFit ();
dataDcs.ClearGraphs ();
->call a store function
dataDcs.ClearFits ();

*/

//_____________________________________________________________________________
AliTRDDataDCS::AliTRDDataDCS () 
  :  chamberByteStatus 	(0),
     preTrigger		(1),
     goofyHv 		(2),
     goofyPeakPos 	(3),
     goofyPeakArea 	(4),
     goofyTemp 		(5),
     goofyPressure 	(6),
     goofyVelocity 	(7),
     goofyGain 		(8),
     goofyCO2 		(9),
     goofyN2 		(10),
     gasO2 		(11),
     gasOverpressure 	(12),
     envTemp  		(13),
     hvAnodeImon 	(14),
     hvDriftImon 	(15),
     hvAnodeUmon 	(16),
     hvDriftUmon 	(17),
     adcClkPhase 	(18),
     atmPressure 	(19),
     luminosity 	(20),
     magneticField 	(21),
     fNAlias		(22),
     graphsAreIni 	(false),
     fitsAreIni 	(false)
{
  Init ();
}

//_____________________________________________________________________________
AliTRDDataDCS::~AliTRDDataDCS ()
{
  ClearFits ();
  ClearGraphs ();
}

//_____________________________________________________________________________
void AliTRDDataDCS::Init ()
{
  SetConf (chamberByteStatus, "trd_chamberByteStatus%03d", 'c', 540, true, false, 10,10,0,2); 
  SetConf (preTrigger, 		"trd_preTrigger", 		'c', 1, true, true, 	10,10,0,2);
  SetConf (goofyHv, 			"trd_goofyHv", 			'f', 1, true, true, 	10,10,0,2);
  SetConf (goofyPeakPos, 		"trd_goofyPeakPos%02d", 'f', 2, true, true, 	10,10,0,2);
  SetConf (goofyPeakArea, 	"trd_goofyPeakArea%02d",'f', 2, true, true, 	10,10,0,2);
  SetConf (goofyTemp, 		"trd_goofyTemp%02d", 	'f', 2, true, true, 	10,10,0,2);
  SetConf (goofyPressure, 	"trd_goofyPressure", 	'f', 1, true, true, 	10,10,0,2);
  SetConf (goofyVelocity, 	"trd_goofyVelocity", 	'f', 1, true, true, 	10,10,0,2);
  SetConf (goofyGain, 		"trd_goofyGain%02d", 	'f', 2, true, true, 	10,10,0,2);
  SetConf (goofyCO2, 			"trd_goofyCO2", 		'f', 1, true, true, 	10,10,0,2);
  SetConf (goofyN2, 			"trd_goofyN2", 			'f', 1, true, true, 	10,10,0,2);
  SetConf (gasO2, 			"trd_gasO2", 			'f', 1, true, true, 	10,10,0,2);
  SetConf (gasOverpressure, 	"trd_gasOverpressure", 	'f', 1, true, true, 	10,10,0,2);
  SetConf (envTemp, 			"trd_envTemp%03d", 		'f', 540, true, true, 	10,10,0,2);
  SetConf (hvAnodeImon, 		"trd_hvAnodeImon%03d", 	'f', 540, true, true, 	10,10,0,2);
  SetConf (hvDriftImon, 		"trd_hvDriftImon%03d", 	'f', 540, true, true, 	10,10,0,2);
  SetConf (hvAnodeUmon, 		"trd_hvAnodeUmon%03d", 	'f', 540, true, true, 	10,10,0,2);
  SetConf (hvDriftUmon, 		"trd_hvDriftUmon%03d", 	'f', 540, true, true, 	10,10,0,2);
  SetConf (adcClkPhase, 		"trd_adcClkPhase", 		'f', 1, true, true, 	10,10,0,2);
  SetConf (atmPressure, 		"trd_atmPressure", 		'f', 1, true, true, 	10,10,0,2);
  SetConf (luminosity, 		"trd_luminosity", 		'f', 1, true, true, 	10,10,0,2);
  SetConf (magneticField, 	"trd_magneticField", 	'f', 1, true, true, 	10,10,0,2);		
  
  for (UInt_t iAlias = 0; iAlias < fNAlias; iAlias++) {
    fDatas[iAlias].graph.SetOwner (1);
    fDatas[iAlias].fit.SetOwner (1);
  }
}

//_____________________________________________________________________________
void AliTRDDataDCS::SetConf (UInt_t iAlias, const char * amanda, 
			     char dataType, UInt_t nChannel, 
			     Bool_t enableGraph, Bool_t enableFit, Int_t kMinPoints, 
			     Int_t kIter, Double_t kMaxDelta, Int_t kFitReq)
{
  if (iAlias >= fNAlias) {
    AliWarning (Form("Alias %d is not correct", iAlias));
    return;
  }
  
  fConfs[iAlias].amanda = amanda;
  fConfs[iAlias].dataType = dataType;
  fConfs[iAlias].nChannel = nChannel;
  fConfs[iAlias].enableGraph = enableGraph;
  fConfs[iAlias].enableFit = enableFit;
  fConfs[iAlias].kMinPoints = kMinPoints;
  fConfs[iAlias].kIter = kIter;
  fConfs[iAlias].kMaxDelta = kMaxDelta;
  fConfs[iAlias].kFitReq = kFitReq;
  
}

//_____________________________________________________________________________
void AliTRDDataDCS::InitFits ()
{
  if (fitsAreIni)
    return;
  
  UInt_t nChannel;
  for (UInt_t iAlias = 0; iAlias < fNAlias; iAlias++) {
    nChannel = fConfs[iAlias].enableFit ? fConfs[iAlias].nChannel : 0;
    fDatas[iAlias].fit.Expand (nChannel);	
  }
  
  fitsAreIni = true;
}

//_____________________________________________________________________________
void AliTRDDataDCS::InitGraphs ()
{
  if (graphsAreIni)
    return;
  
  UInt_t nChannel;	
  for (UInt_t iAlias = 0; iAlias < fNAlias; iAlias++) {
    nChannel = fConfs[iAlias].enableGraph ? fConfs[iAlias].nChannel : 0;	 
    fDatas[iAlias].graph.Expand (nChannel);
  }
  
  graphsAreIni = true;
}

//_____________________________________________________________________________
void AliTRDDataDCS::ClearFits ()
{
  for (UInt_t iAlias = 0; iAlias < fNAlias; iAlias++) {
    fDatas[iAlias].fit.Clear ();
    fDatas[iAlias].fit.Expand (0);
  }
  
  fitsAreIni = false;
}

//_____________________________________________________________________________
void AliTRDDataDCS::ClearGraphs ()
{	
  for (UInt_t iAlias = 0; iAlias < fNAlias; iAlias++) {
    fDatas[iAlias].graph.Clear ();
    fDatas[iAlias].graph.Expand (0);
  }
  
  graphsAreIni = false;
}

//_____________________________________________________________________________
void AliTRDDataDCS::Streamer(TBuffer &R__b) {
  
  if (R__b.IsReading()) {
    
    R__b.ReadBool (graphsAreIni);
    R__b.ReadBool (fitsAreIni);
    
    
    for (UInt_t iAlias=0; iAlias<fNAlias; iAlias++) {
      TString::ReadString (R__b, TString::Class());
      R__b.ReadChar (fConfs[iAlias].dataType);
      R__b.ReadUInt (fConfs[iAlias].nChannel);
      R__b.ReadBool (fConfs[iAlias].enableGraph);
      R__b.ReadBool (fConfs[iAlias].enableFit);
      R__b.ReadInt  (fConfs[iAlias].kMinPoints);
      R__b.ReadInt  (fConfs[iAlias].kIter);
      R__b.ReadDouble  (fConfs[iAlias].kMaxDelta);
      R__b.ReadInt  (fConfs[iAlias].kFitReq);
    }
    
    
    if (graphsAreIni) {
      for (UInt_t iAlias=0; iAlias<fNAlias; iAlias++)
	fDatas[iAlias].graph = *(TObjArray*)R__b.ReadObject (TObjArray::Class());
    }
    
    if (fitsAreIni) {
      for (UInt_t iAlias=0; iAlias<fNAlias; iAlias++)
	fDatas[iAlias].fit = *(TObjArray*)R__b.ReadObject (TObjArray::Class());
    }
    
    AliInfo (Form("Read %d octets to the stream (%d Ko)", R__b.Length(), R__b.Length()/1024));
    
  } else {
    R__b.WriteBool (graphsAreIni);
    R__b.WriteBool (fitsAreIni);
    
    for (UInt_t iAlias = 0; iAlias < fNAlias; iAlias++) {
      TString::WriteString (R__b, &fConfs[iAlias].amanda);
      R__b.WriteChar (fConfs[iAlias].dataType);
      R__b.WriteUInt (fConfs[iAlias].nChannel);
      R__b.WriteBool (fConfs[iAlias].enableGraph);
      R__b.WriteBool (fConfs[iAlias].enableFit);
      R__b.WriteInt  (fConfs[iAlias].kMinPoints);
      R__b.WriteInt  (fConfs[iAlias].kIter);
      R__b.WriteDouble  (fConfs[iAlias].kMaxDelta);
      R__b.WriteInt  (fConfs[iAlias].kFitReq);
    }
    
    if (graphsAreIni) {
      for (UInt_t iAlias = 0; iAlias < fNAlias; iAlias++) 
	R__b.WriteObject ((TObject*)&fDatas[iAlias].graph);
    }
    
    if (fitsAreIni) {
      for (UInt_t iAlias = 0; iAlias < fNAlias; iAlias++) 
	R__b.WriteObject ((TObject*)&fDatas[iAlias].fit);
    }
    
    AliInfo (Form("Write %d octets to the stream (%d Ko)", R__b.Length(), R__b.Length()/1024));
  }
}

//_____________________________________________________________________________
Bool_t AliTRDDataDCS::ExtractDCS (TMap * dcsAlias)
{
  if (dcsAlias == 0x0) {
    AliWarning ("No DCS Map");
    return kFALSE;
  }
  
  ClearGraphs ();
  InitGraphs ();
  
  TGraph * graphTemp;
  UInt_t nChannel;
  
  for (UInt_t iAlias = 0; iAlias < fNAlias; iAlias++){
    
    // extract dcs only when it is needed
    nChannel = fDatas[iAlias].graph.GetSize ();
    
    for (UInt_t iSensor=0; iSensor < nChannel; iSensor++) {
      
      graphTemp = FindAndMakeGraph (dcsAlias, 
				    Form (fConfs[iAlias].amanda.Data(), iSensor), 
				    fConfs[iAlias].dataType); 
      
      fDatas[iAlias].graph.AddAt (graphTemp, iSensor);
    }
  }
  
  return kTRUE;
}

//_____________________________________________________________________________
TGraph * AliTRDDataDCS::FindAndMakeGraph (TMap * dcsMap, const char * amandaStr, 
					  char dataType)
{
  TGraph * graph;
  
  TPair * pair = (TPair*)dcsMap->FindObject(amandaStr);
  if (pair == 0x0) {
    AliWarning (Form("Can't find %s in dcsMap", amandaStr));
    return 0x0;
  }
  
  TObjArray * valueSet = (TObjArray*)pair->Value();
  
  //
  // Make graph of values read from DCS map
  //   (spline fit parameters will subsequently be obtained from this graph) 
  //
  
  Int_t nEntries = valueSet->GetEntriesFast();
  if (nEntries == 0) {
    AliWarning (Form("Entry %s in dcsMap contain no datas", amandaStr));
    return 0x0;
  }
  
  Float_t * x = new Float_t [nEntries];
  Float_t * y = new Float_t [nEntries];
  
  Bool_t ok = true;
  Int_t time0 = 0;
  Int_t iEntries;
  for (iEntries = 0;  iEntries< nEntries; iEntries++) {
    
    AliDCSValue * val = (AliDCSValue *)valueSet->At(iEntries);
    
    if (val == 0x0) { 
      ok = false;
      AliError (Form("Entry %s at %d contain no datas", amandaStr, iEntries));
      break;
    }
    
    if (time0 == 0) 
      time0 = val->GetTimeStamp();
    
    
    x[iEntries] = (val->GetTimeStamp() - time0)/3600.0; // give times in fractions of hours 
    
    switch (dataType) {
    case 'f' :
      y[iEntries] = val->GetFloat();
      break;
      
    case 'c' :
      y[iEntries] = (Float_t) val->GetChar();
      break;
      
    default :
      ok = false;
      AliError (Form("Bad type for entry %s", amandaStr));
      break;
    }
  }
  
  if (ok)
    graph = new TGraph (iEntries, x, y);
  else
    graph = 0x0;
  
  delete [] x;
  delete [] y;
  
  return graph;
}


//_____________________________________________________________________________
Bool_t AliTRDDataDCS::PerformFit ()
{
  AliSplineFit * fitTemp;
  
  ClearFits ();
  InitFits ();
  
  UInt_t nChannel;
  
  for (UInt_t iAlias = 0; iAlias < fNAlias; iAlias++){
    
    // perform fit only when it is needed
    nChannel = fDatas[iAlias].fit.GetSize ();
    
    for (UInt_t iSensor=0; iSensor < nChannel; iSensor++) {
      
      fitTemp = Fit ((TGraph*)fDatas[iAlias].graph[iSensor],
		     fConfs[iAlias].kMinPoints,
		     fConfs[iAlias].kIter,
		     fConfs[iAlias].kMaxDelta,
		     fConfs[iAlias].kFitReq); 
      
      if (fitTemp == 0x0)
	AliInfo (Form ("Can't fit %s", Form (fConfs[iAlias].amanda.Data(), iSensor)));
      
      fDatas[iAlias].fit.AddAt (fitTemp, iSensor);
    }
  }
  
  return kTRUE;
  
}

//_____________________________________________________________________________
AliSplineFit * AliTRDDataDCS::Fit (TGraph * graph, 
				   Int_t  kMinPoints, Int_t  kIter, 
				   Double_t  kMaxDelta, Int_t  kFitReq)
{
  if (graph == 0x0) {
    AliError ("No graph for fit");
    return 0x0;
  }
  
  AliSplineFit * spline = new AliSplineFit ();
  spline->InitKnots (new TGraph (*graph), kMinPoints, kIter, kMaxDelta);
  spline->SplineFit (kFitReq);
  spline->Cleanup();   // delete also new TGraph (*graph)
  
  return spline;
}


//_____________________________________________________________________________
void AliTRDDataDCS::Print (Option_t* option) const
{
  if (option[0]=='g' || option[0]=='\0'){
    
    if (graphsAreIni){
      
      for (UInt_t iAlias = 0; iAlias < fNAlias; iAlias++){
	
	for (Int_t iSensor = 0; iSensor < fDatas[iAlias].graph.GetSize(); iSensor++) {
	  
	  if (fDatas[iAlias].graph[iSensor] != 0x0)
	    AliInfo (Form("Graph %s contain %d point(s)", 
			  Form(fConfs[iAlias].amanda, iSensor), 
			  ((TGraph*)(fDatas[iAlias].graph[iSensor]))->GetN ()));
	}
      }
    }else{
      AliInfo ("Graphs don't exist");
    }
  }
  
  
  if (option[0] == 'f' || option[0]=='\0'){
    
    AliInfo ("no print for fit");
    
  }
  
}

//_____________________________________________________________________________
TGraph * AliTRDDataDCS::GetGraph (UInt_t iAlias, UInt_t iChannel) const
{
  if (iAlias >= fNAlias) {
    AliWarning (Form("Alias %d is not correct", iAlias));
    return 0x0;
  }
  
  if (iChannel >= (UInt_t)fDatas[iAlias].graph.GetSize()) {
    AliWarning (Form("Alias %s is not correct", 
		     Form(fConfs[iAlias].amanda.Data(), iChannel)));
    return 0x0;
  }
  
  return (TGraph *)fDatas[iAlias].graph[iChannel];
}

//_____________________________________________________________________________
AliSplineFit * AliTRDDataDCS::GetFit (UInt_t iAlias, UInt_t iChannel) const
{
  if (iAlias >= fNAlias) {
    AliWarning (Form("Alias %d is not correct", iAlias));
    return 0x0;
  }
  
  if (iChannel >= (UInt_t)fDatas[iAlias].fit.GetSize()) {
    AliWarning (Form("Alias %s is not correct", 
		     Form(fConfs[iAlias].amanda.Data(), iChannel)));
    return 0x0;
  }
  
  return (AliSplineFit *)fDatas[iAlias].fit[iChannel];
}

//_____________________________________________________________________________
TString AliTRDDataDCS::GetAmandaStr (UInt_t iAlias) const 
{
  if (iAlias < fNAlias)
    return fConfs[iAlias].amanda;
  else return TString ();
}

//_____________________________________________________________________________
UInt_t AliTRDDataDCS::GetNChannel (UInt_t iAlias) const 
{
  if (iAlias < fNAlias)
    return fConfs[iAlias].nChannel;
  else return 0;
}
