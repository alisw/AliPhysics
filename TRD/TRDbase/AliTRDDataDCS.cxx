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
// Typical use :                                                          //
//                                                                        //
//   AliTRDDataDCS dataDcs;                                               //
//   dataDcs.ExtractDCS(dcs);                                             //
//   ->call a storeRef function                                           //
//   dataDcs.PerformFit();                                                //                        
//   dataDcs.ClearGraphs();                                               //
//   ->call a store function                                              //
//   dataDcs.ClearFits();                                                 //
//                                                                        //
// Author:                                                                //
//   W. Monange   (wilfried.monange@free.fr)                              //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#include <TGraph.h>
#include <TObjArray.h>
#include <TMap.h>

#include "AliDCSValue.h"
#include "AliDCSValue.h"
#include "AliLog.h"
#include "AliSplineFit.h"

#include "AliTRDDataDCS.h"

ClassImp(AliTRDDataDCS)

//_____________________________________________________________________________
AliTRDDataDCS::AliTRDDataDCS() 
  :TNamed()       
  ,fGraphsAreIni(kFALSE)
  ,fFitsAreIni(kFALSE)
  ,fNAlias(22)
{
  //
  // Default constructor
  //

  Init ();

}

//_____________________________________________________________________________
AliTRDDataDCS::~AliTRDDataDCS()
{
  //
  // Destructor
  //

  ClearFits();
  ClearGraphs();

}

//_____________________________________________________________________________
void AliTRDDataDCS::Init()
{
  //
  // Initialization
  //

  SetConf(kChamberByteStatus, "trd_chamberByteStatus%03d",'c', 540, kTRUE,kFALSE,       10,10,0,2); 
  SetConf(kPreTrigger, 	      "trd_preTrigger", 	  'c',   1, kTRUE, kTRUE, 	10,10,0,2);
  SetConf(kGoofyHv, 	      "trd_goofyHv", 		  'f',   1, kTRUE, kTRUE, 	10,10,0,2);
  SetConf(kGoofyPeakPos,      "trd_goofyPeakPos%02d",     'f',   2, kTRUE, kTRUE, 	10,10,0,2);
  SetConf(kGoofyPeakArea,     "trd_goofyPeakArea%02d",    'f',   2, kTRUE, kTRUE, 	10,10,0,2);
  SetConf(kGoofyTemp, 	      "trd_goofyTemp%02d", 	  'f',   2, kTRUE, kTRUE, 	10,10,0,2);
  SetConf(kGoofyPressure,     "trd_goofyPressure", 	  'f',   1, kTRUE, kTRUE, 	10,10,0,2);
  SetConf(kGoofyVelocity,     "trd_goofyVelocity", 	  'f',   1, kTRUE, kTRUE, 	10,10,0,2);
  SetConf(kGoofyGain, 	      "trd_goofyGain%02d", 	  'f',   2, kTRUE, kTRUE, 	10,10,0,2);
  SetConf(kGoofyCO2, 	      "trd_goofyCO2", 		  'f',   1, kTRUE, kTRUE, 	10,10,0,2);
  SetConf(kGoofyN2, 	      "trd_goofyN2", 		  'f',   1, kTRUE, kTRUE, 	10,10,0,2);
  SetConf(kGasO2, 	      "trd_gasO2", 		  'f',   1, kTRUE, kTRUE, 	10,10,0,2);
  SetConf(kGasOverpressure,   "trd_gasOverpressure", 	  'f',   1, kTRUE, kTRUE, 	10,10,0,2);
  SetConf(kEnvTemp, 	      "trd_envTemp%03d", 	  'f', 540, kTRUE, kTRUE, 	10,10,0,2);
  SetConf(kHvAnodeImon,       "trd_hvAnodeImon%03d",   	  'f', 540, kTRUE, kTRUE, 	10,10,0,2);
  SetConf(kHvDriftImon,       "trd_hvDriftImon%03d", 	  'f', 540, kTRUE, kTRUE, 	10,10,0,2);
  SetConf(kHvAnodeUmon,       "trd_hvAnodeUmon%03d", 	  'f', 540, kTRUE, kTRUE, 	10,10,0,2);
  SetConf(kHvDriftUmon,       "trd_hvDriftUmon%03d", 	  'f', 540, kTRUE, kTRUE, 	10,10,0,2);
  SetConf(kAdcClkPhase,       "trd_adcClkPhase", 	  'f',   1, kTRUE, kTRUE, 	10,10,0,2);
  SetConf(kAtmPressure,       "trd_atmPressure", 	  'f',   1, kTRUE, kTRUE, 	10,10,0,2);
  SetConf(kLuminosity, 	      "trd_luminosity", 	  'f',   1, kTRUE, kTRUE, 	10,10,0,2);
  SetConf(kMagneticField,     "trd_magneticField", 	  'f',   1, kTRUE, kTRUE, 	10,10,0,2);		
  
  for (UInt_t iAlias = 0; iAlias < fNAlias; iAlias++) {
    fDatas[iAlias].GetGraph().SetOwner(1);
    fDatas[iAlias].GetFit().SetOwner(1);
  }

}

//_____________________________________________________________________________
void AliTRDDataDCS::SetConf(UInt_t iAlias, const char *amanda, 
			     char dataType, UInt_t nChannel, 
			     Bool_t enableGraph, Bool_t enableFit, Int_t minPoints, 
			     Int_t iter, Double_t maxDelta, Int_t fitReq)
{
  //
  // Configure a DCS alias
  //

  if (iAlias >= fNAlias) {
    AliWarning (Form("Alias %d is not correct", iAlias));
    return;
  }
  
  fConfs[iAlias].SetAmanda(amanda);
  fConfs[iAlias].SetDataType(dataType);
  fConfs[iAlias].SetNChannel(nChannel);
  fConfs[iAlias].SetEnableGraph(enableGraph);
  fConfs[iAlias].SetEnableFit(enableFit);
  fConfs[iAlias].SetMinPoints(minPoints);
  fConfs[iAlias].SetIter(iter);
  fConfs[iAlias].SetMaxDelta(maxDelta);
  fConfs[iAlias].SetFitReq(fitReq);
  
}

//_____________________________________________________________________________
void AliTRDDataDCS::InitFits()
{
  //
  // Initialize the fits
  //

  if (fFitsAreIni)
    return;
  
  UInt_t nChannel;
  for (UInt_t iAlias = 0; iAlias < fNAlias; iAlias++) {
    nChannel = fConfs[iAlias].GetEnableFit() ? fConfs[iAlias].GetNChannel() : 0;
    fDatas[iAlias].GetFit().Expand(nChannel);	
  }
  
  fFitsAreIni = kTRUE;
}

//_____________________________________________________________________________
void AliTRDDataDCS::InitGraphs()
{
  //
  // Initialize the graphs
  //

  if (fGraphsAreIni)
    return;
  
  UInt_t nChannel;	
  for (UInt_t iAlias = 0; iAlias < fNAlias; iAlias++) {
    nChannel = fConfs[iAlias].GetEnableGraph() ? fConfs[iAlias].GetNChannel() : 0;	 
    fDatas[iAlias].GetGraph().Expand(nChannel);
  }
  
  fGraphsAreIni = kTRUE;
}

//_____________________________________________________________________________
void AliTRDDataDCS::ClearFits()
{
  //
  // Clear the fits
  //

  for (UInt_t iAlias = 0; iAlias < fNAlias; iAlias++) {
    fDatas[iAlias].GetFit().Clear();
    fDatas[iAlias].GetFit().Expand(0);
  }
  
  fFitsAreIni = kFALSE;
}

//_____________________________________________________________________________
void AliTRDDataDCS::ClearGraphs()
{
  //
  // Clear the grpahs
  //

  for (UInt_t iAlias = 0; iAlias < fNAlias; iAlias++) {
    fDatas[iAlias].GetGraph().Clear();
    fDatas[iAlias].GetGraph().Expand(0);
  }
  
  fGraphsAreIni = kFALSE;
}

//_____________________________________________________________________________
Bool_t AliTRDDataDCS::ExtractDCS (TMap *dcsAlias)
{
  //
  // Extract the DCS information
  //

  if (dcsAlias == 0x0) {
    AliWarning ("No DCS Map");
    return kFALSE;
  }
  
  ClearGraphs();
  InitGraphs();
  
  TGraph *graphTemp;
  UInt_t nChannel;
  
  for (UInt_t iAlias = 0; iAlias < fNAlias; iAlias++){
    
    // extract dcs only when it is needed
    nChannel = fDatas[iAlias].GetGraph().GetSize();
    
    for (UInt_t iSensor=0; iSensor < nChannel; iSensor++) {
      
      graphTemp = FindAndMakeGraph(dcsAlias, 
				   Form(fConfs[iAlias].GetAmanda().Data(), iSensor), 
				   fConfs[iAlias].GetDataType()); 
      
      fDatas[iAlias].GetGraph().AddAt(graphTemp, iSensor);
    }
  }
  
  return kTRUE;
}

//_____________________________________________________________________________
TGraph *AliTRDDataDCS::FindAndMakeGraph (TMap * const dcsMap
                                       , const char *amandaStr
				       , char dataType)
{
  //
  // Create the graphs
  //

  TGraph *graph;
  
  TPair *pair = (TPair *) dcsMap->FindObject(amandaStr);
  if (pair == 0x0) {
    AliWarning (Form("Can't find %s in dcsMap", amandaStr));
    return 0x0;
  }
  
  TObjArray *valueSet = (TObjArray *) pair->Value();
  
  //
  // Make graph of values read from DCS map
  //   (spline fit parameters will subsequently be obtained from this graph) 
  //
  
  Int_t nEntries = valueSet->GetEntriesFast();
  if (nEntries == 0) {
    AliWarning (Form("Entry %s in dcsMap contain no datas", amandaStr));
    return 0x0;
  }
  
  Float_t * x = new Float_t[nEntries];
  Float_t * y = new Float_t[nEntries];
  
  Bool_t ok = kTRUE;
  Int_t time0 = 0;
  Int_t iEntries;
  for (iEntries = 0;  iEntries< nEntries; iEntries++) {
    
    AliDCSValue *val = (AliDCSValue *) valueSet->At(iEntries);
    
    if (val == 0x0) { 
      ok = false;
      AliError(Form("Entry %s at %d contain no datas", amandaStr, iEntries));
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
      AliError(Form("Bad type for entry %s", amandaStr));
      break;
    }
  }
  
  if (ok)
    graph = new TGraph(iEntries, x, y);
  else
    graph = 0x0;
  
  delete [] x;
  delete [] y;
  
  return graph;
}

//_____________________________________________________________________________
Bool_t AliTRDDataDCS::PerformFit()
{
  //
  // Do the fit
  //

  AliSplineFit *fitTemp;
  
  ClearFits();
  InitFits();
  
  UInt_t nChannel;
  
  for (UInt_t iAlias = 0; iAlias < fNAlias; iAlias++){
    
    // perform fit only when it is needed
    nChannel = fDatas[iAlias].GetFit().GetSize();
    
    for (UInt_t iSensor=0; iSensor < nChannel; iSensor++) {
      
      fitTemp = Fit((TGraph*)fDatas[iAlias].GetGraph(iSensor),
		     fConfs[iAlias].GetMinPoints(),
		     fConfs[iAlias].GetIter(),
		     fConfs[iAlias].GetMaxDelta(),
		     fConfs[iAlias].GetFitReq()); 
      
      if (fitTemp == 0x0)
	AliInfo(Form("Can't fit %s", Form(fConfs[iAlias].GetAmanda().Data(), iSensor)));
      
      fDatas[iAlias].GetFit().AddAt(fitTemp, iSensor);
    }
  }
  
  return kTRUE;
  
}

//_____________________________________________________________________________
AliSplineFit *AliTRDDataDCS::Fit(const TGraph * const graph, 
			         Int_t  minPoints, Int_t  iter, 
			         Double_t  maxDelta, Int_t  fitReq)
{
  //
  // Do the spline fit
  //

  if (graph == 0x0) {
    AliError("No graph for fit");
    return 0x0;
  }
  
  AliSplineFit *spline = new AliSplineFit();
  spline->InitKnots(new TGraph (*graph), minPoints, iter, maxDelta);
  spline->SplineFit(fitReq);
  spline->Cleanup();   // delete also new TGraph (*graph)
  
  return spline;
}

//_____________________________________________________________________________
void AliTRDDataDCS::Print(const Option_t * const option) const
{
  //
  // Print function
  //

  if (option[0]=='g' || option[0]=='\0'){
    
    if (fGraphsAreIni){
      
      for (UInt_t iAlias = 0; iAlias < fNAlias; iAlias++){
	
	for (Int_t iSensor = 0; iSensor < fDatas[iAlias].GetGraph().GetSize(); iSensor++) {
	  
	  if (fDatas[iAlias].GetGraph(iSensor) != 0x0)
	    AliInfo(Form("Graph %s contain %d point(s)", 
			  Form(fConfs[iAlias].GetAmanda(), iSensor), 
			  ((TGraph*)(fDatas[iAlias].GetGraph(iSensor)))->GetN()));
	}
      }
    }
    else{
      AliInfo("Graphs don't exist");
    }
  }
  
  
  if (option[0] == 'f' || option[0]=='\0'){
    
    AliInfo("no print for fit");
    
  }
  
}

//_____________________________________________________________________________
TGraph *AliTRDDataDCS::GetGraph(UInt_t iAlias, UInt_t iChannel) const
{
  //
  // Get a graph
  //

  if (iAlias >= fNAlias) {
    AliWarning(Form("Alias %d is not correct", iAlias));
    return 0x0;
  }
  
  if (iChannel >= (UInt_t) fDatas[iAlias].GetGraph().GetSize()) {
    AliWarning(Form("Alias %s is not correct", 
		     Form(fConfs[iAlias].GetAmanda().Data(), iChannel)));
    return 0x0;
  }
  
  return (TGraph *)fDatas[iAlias].GetGraph(iChannel);
}

//_____________________________________________________________________________
AliSplineFit *AliTRDDataDCS::GetFit(UInt_t iAlias, UInt_t iChannel) const
{
  //
  // Get the spline fit
  //

  if (iAlias >= fNAlias) {
    AliWarning (Form("Alias %d is not correct", iAlias));
    return 0x0;
  }
  
  if (iChannel >= (UInt_t) fDatas[iAlias].GetFit().GetSize()) {
    AliWarning(Form("Alias %s is not correct", 
		     Form(fConfs[iAlias].GetAmanda().Data(), iChannel)));
    return 0x0;
  }
  
  return (AliSplineFit *) fDatas[iAlias].GetFit(iChannel);
}

//_____________________________________________________________________________
TString AliTRDDataDCS::GetAmandaStr (UInt_t iAlias) const 
{
  //
  // Return the AMANDA string
  //

  if (iAlias < fNAlias)
    return fConfs[iAlias].GetAmanda();
  else return TString ();
}

//_____________________________________________________________________________
UInt_t AliTRDDataDCS::GetNChannel (UInt_t iAlias) const 
{
  //
  // Get the channel number
  //

  if (iAlias < fNAlias)
    return fConfs[iAlias].GetNChannel();
  else return 0;
}
