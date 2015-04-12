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
// This class perform operation on DCS Sensor                             //
// The configuration of each sensor is included inside this class         //
// Use the methode GetList to get all the configuration                   //
//                                                                        //
// Author:                                                                //
//   W. Monange   (w.monange@gsi.de)                                      //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#include <TObjArray.h>
#include <TClonesArray.h>

#include "AliTRDSensorArray.h"
#include "AliTRDSensor.h"

ClassImp(AliTRDSensorArray)

//_____________________________________________________________________________
AliTRDSensorArray::AliTRDSensorArray() 
  :AliDCSSensorArray()
  ,fAmanda("")
  ,fStoreName("")					
{
  //
  // Default constructor
  //
	
}

//_____________________________________________________________________________
AliTRDSensorArray::AliTRDSensorArray(const char *amanda 
				   , const char *storeName
				   , Float_t /*diffCut*/
				   , TClonesArray * const trdSensor) 
  :AliDCSSensorArray()
  ,fAmanda(amanda)
  ,fStoreName(storeName)				
{
  //
  // Constructor set fMinGraph to 0, fValCut to 0, fDiffCut to 0
  //

	fSensors 	= trdSensor;
	fMinGraph 	= 0;
	fValCut		= -1;
	fDiffCut	= -1;
	Int_t entries = fSensors->GetEntriesFast();
	if(entries > 1){
	  for(Int_t k = 0; k < entries; k++){
	    TString name (Form(amanda, k));
	    //printf("name is %s of %d\n",(const char*)name,k);
	    ((AliDCSSensor *) fSensors->UncheckedAt(k))->SetStringID(name);
	  }
	}
	else{
	  TString name (amanda);
	  //printf("name is %s\n",(const char*)name);
	  ((AliDCSSensor *) fSensors->UncheckedAt(0))->SetStringID(name);
	}

}

//_____________________________________________________________________________
AliTRDSensorArray::AliTRDSensorArray(const AliTRDSensorArray & source) 
  :AliDCSSensorArray(source)
  ,fAmanda(source.fAmanda)
  ,fStoreName(source.fStoreName)				
{
  //
  // Copy constructor
  //

	fSensors 	= source.fSensors;
	fMinGraph	= 0;
	fValCut		= -1;
	fDiffCut	= -1;

}
	
//_____________________________________________________________________________
AliTRDSensorArray::~AliTRDSensorArray()
{
  //
  // Destructor
  //
	
}

//_____________________________________________________________________________		
AliTRDSensorArray &AliTRDSensorArray::operator=(const AliTRDSensorArray &source)
{
  //
  // Assignment operator
  //

	if (&source == this) return *this;
	new (this) AliTRDSensorArray (source);
  
	return *this;  
}

//_____________________________________________________________________________
TObjArray *AliTRDSensorArray::GetList() 
{
  //
  // Return TObjArray with a list of AliTRDSensorArray corresponding to each
  // group of sensor 
  //
	
	TObjArray * list = new TObjArray (23);
	list->SetOwner (kTRUE);
	AliTRDSensorArray * aH = 0x0;
	
	// generic list of sensors
	TClonesArray listSensor540 ("AliTRDSensor", 540);	
	TClonesArray listSensor2   ("AliTRDSensor", 2); 
	TClonesArray listSensor1   ("AliTRDSensor", 1); 
	listSensor540.SetOwner (kTRUE);
	listSensor2.SetOwner (kTRUE);
	listSensor1.SetOwner (kTRUE);
	for (Int_t i = 0; i < 540; i++) {
		new(listSensor540[i]) AliTRDSensor (i, 0, 0, 0);
	}
	for (Int_t i = 0; i < 2; i++)
		new(listSensor2[i]) AliTRDSensor (i, 0, 0, 0);
	new(listSensor1[0]) AliTRDSensor (0, 0, 0, 0);
	
	
	// now create and populate   
	aH = new AliTRDSensorArray ("trd_chamberStatus%03d", "trd_chamberStatus", 
				    0.5, (TClonesArray*)listSensor540.Clone ());
	list->Add (aH);
	
	aH = new AliTRDSensorArray ("trd_goofieHv", 			"trd_goofieHv",  
				    -1, (TClonesArray*)listSensor1.Clone ());
	list->Add (aH);
	aH = new AliTRDSensorArray ("trd_goofiePeakPos%02d",	"trd_goofiePeakPos",   
				    -1, (TClonesArray*)listSensor2.Clone ());
	list->Add (aH);
	aH = new AliTRDSensorArray ("trd_goofiePeakArea%02d","trd_goofiePeakArea",   
				    -1, (TClonesArray*)listSensor2.Clone ());
	list->Add (aH);
	aH = new AliTRDSensorArray ("trd_goofieTemp%02d", 	"trd_goofieTemp",  
				    -1, (TClonesArray*)listSensor2.Clone ());
	list->Add (aH);
	aH = new AliTRDSensorArray ("trd_goofiePressure", 	"trd_goofiePressure",  
				    -1, (TClonesArray*)listSensor1.Clone ());
	list->Add (aH);
	aH = new AliTRDSensorArray ("trd_goofieVelocity", 	"trd_goofieVelocity",  
				    -1, (TClonesArray*)listSensor1.Clone ());
	list->Add (aH);
	aH = new AliTRDSensorArray ("trd_goofieGain%02d", 	"trd_goofieGain",  
				    -1, (TClonesArray*)listSensor2.Clone ());
	list->Add (aH);
	aH = new AliTRDSensorArray ("trd_goofieCO2", 		"trd_goofieCO2",  
				    -1, (TClonesArray*)listSensor1.Clone ());
	list->Add (aH);
	aH = new AliTRDSensorArray ("trd_goofieN2", 			"trd_goofieN2", 
				    -1, (TClonesArray*)listSensor1.Clone ());
	list->Add (aH);
	aH = new AliTRDSensorArray ("trd_gasO2", 			"trd_gasO2",  
				    -1, (TClonesArray*)listSensor1.Clone ());
	list->Add (aH);
	aH = new AliTRDSensorArray ("trd_gasH2O", 			"trd_gasH2O",  
				    -1, (TClonesArray*)listSensor1.Clone ());
	list->Add (aH);
	aH = new AliTRDSensorArray ("trd_gasCO2", 			"trd_gasCO2",  
				    -1, (TClonesArray*)listSensor1.Clone ());
	list->Add (aH);
	aH = new AliTRDSensorArray ("trd_gasOverpressure", 	"trd_gasOverpressure",  
				    -1, (TClonesArray*)listSensor1.Clone ());
	list->Add (aH);
	aH = new AliTRDSensorArray ("trd_envTemp%03d", 		"trd_envTemp",  
				    -1, (TClonesArray*)listSensor540.Clone ());
	list->Add (aH);
	aH = new AliTRDSensorArray ("trd_hvAnodeImon%03d", 	"trd_hvAnodeImon",	 
				    -1, (TClonesArray*)listSensor540.Clone ());
	list->Add (aH);
	aH = new AliTRDSensorArray ("trd_hvDriftImon%03d", 	"trd_hvDriftImon",  
				    -1, (TClonesArray*)listSensor540.Clone ());
	list->Add (aH);
	aH = new AliTRDSensorArray ("trd_hvAnodeUmon%03d",	"trd_hvAnodeUmon",  
				    -1, (TClonesArray*)listSensor540.Clone ());
	list->Add (aH);
	aH = new AliTRDSensorArray ("trd_hvDriftUmon%03d", 	"trd_hvDriftUmon",  
				    -1, (TClonesArray*)listSensor540.Clone ());
	list->Add (aH);
	aH = new AliTRDSensorArray ("trd_gaschromatographXe", 	"trd_gaschromatographXe",  
				    -1, (TClonesArray*)listSensor1.Clone ());
	list->Add (aH);
	aH = new AliTRDSensorArray ("trd_gaschromatographCO2", 	"trd_gaschromatographCO2",  
				    -1, (TClonesArray*)listSensor1.Clone ());
	list->Add (aH);
	aH = new AliTRDSensorArray ("trd_gaschromatographN2", 	"trd_gaschromatographN2",  
				    -1, (TClonesArray*)listSensor1.Clone ());
	list->Add (aH);
	                          
	return list;
}

//_____________________________________________________________________________
TMap* AliTRDSensorArray::ExtractDCS(TMap *dcsMap)
{
  //
  // Return Tmap with TGraph inside corresponding to values in dcsMap
  //

	return AliDCSSensorArray::ExtractDCS(dcsMap);
}

//_____________________________________________________________________________
void AliTRDSensorArray::SetGraph(TMap *map)
{
  //
  // Assign list of TGraph to the current instance
  //

	AliDCSSensorArray::SetGraph(map);
}

//_____________________________________________________________________________
Int_t AliTRDSensorArray::GetNGraph() const
{
  //
  // Return the number of TGraph
  //

	Int_t nGraph = 0;
	Int_t nsensors = fSensors->GetEntries();
	
	for (Int_t isensor = 0; isensor < nsensors; isensor++) {
		AliDCSSensor *entry = (AliDCSSensor*)fSensors->At(isensor);
		if (entry->GetGraph () != 0x0)
			nGraph ++;
	} 
	return nGraph;
}


