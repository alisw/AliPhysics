//-*- Mode: C++ -*-
// $Id: AliHLTGlobalCompareFlatComponent.cxx $
/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 *                                                                        *
 * Primary Authors: Steffen Weber                                         *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/** @file    AliHLTGlobalCompareFlatComponent.cxx
    @author  Steffen Weber
    @brief   Compare flat events from different inputs
*/

#include "TMap.h"
#include "TSystem.h"
#include "TTimeStamp.h"
#include "TObjString.h"
#include "TList.h"
#include "THnSparse.h"
#include "AliESDEvent.h"
#include "AliFlatESDEvent.h"
#include "AliFlatESDTrigger.h"
#include "AliFlatESDV0.h"
#include "AliHLTErrorGuard.h"
#include "AliHLTDataTypes.h"
#include "AliHLTGlobalCompareFlatComponent.h"
#include "AliHLTITSClusterDataFormat.h"
#include "AliHLTTPCDefinitions.h"
#include "TTree.h"
#include "AliCDBEntry.h"
#include "AliCDBManager.h"

using namespace std;

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTGlobalCompareFlatComponent)

void AliHLTGlobalCompareFlatComponent::printDiff( string name, double val1, double val2){
	double relDiff = ( val1 != 0 || val2!=0 ) ? (val1-val2)/(fabs(val1) + fabs(val2)): 0;
	int diff = 0;
	if (relDiff > 1e-6) diff = 1;
	else if(relDiff < -1e-6) diff = -1;
	outFile<<name<<"\t" << val1 << "\t" << val2 <<"\t" << diff << "\n";
}
void AliHLTGlobalCompareFlatComponent::printDiff( string name, TString val1, TString val2){
	outFile << name << "\t" << "\t\"" << val1 <<"\"\t\"" << val2 <<"\"\t" << (val1.EqualTo(val2) ?0:1)<<"\n";
}

/*
 * ---------------------------------------------------------------------------------
 *                            Constructor / Destructor
 * ---------------------------------------------------------------------------------
 */

// #################################################################################
AliHLTGlobalCompareFlatComponent::AliHLTGlobalCompareFlatComponent() :
  AliHLTProcessor()
  {
  // an example component which implements the ALICE HLT processor
  // interface and does some analysis on the input raw data
  //
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
  //
  // NOTE: all helper classes should be instantiated in DoInit()
}

// #################################################################################
AliHLTGlobalCompareFlatComponent::~AliHLTGlobalCompareFlatComponent() {
  // see header file for class documentation

	
	
}

/*
 * ---------------------------------------------------------------------------------
 * Public functions to implement AliHLTComponent's interface.
 * These functions are required for the registration process
 * ---------------------------------------------------------------------------------
 */

// #################################################################################
const Char_t* AliHLTGlobalCompareFlatComponent::GetComponentID() { 
  // see header file for class documentation
  return "GlobalCompareFlat";
}

// #################################################################################
void AliHLTGlobalCompareFlatComponent::GetInputDataTypes( vector<AliHLTComponentDataType>& list) {
  // see header file for class documentation
	list.push_back(kAliHLTDataTypeFlatESD|kAliHLTDataOriginOut );
}

// #################################################################################
AliHLTComponentDataType AliHLTGlobalCompareFlatComponent::GetOutputDataType() {
  // see header file for class documentation
  return kAliHLTDataTypeTObject|kAliHLTDataOriginOut;
}

// #################################################################################
void AliHLTGlobalCompareFlatComponent::GetOutputDataSize( ULong_t& constBase, Double_t& inputMultiplier ) {
  // see header file for class documentation
  constBase = 100000;
  inputMultiplier = 10.0;
}


// #################################################################################
AliHLTComponent* AliHLTGlobalCompareFlatComponent::Spawn() {
  // see header file for class documentation
  return new AliHLTGlobalCompareFlatComponent;
}

/*
 * ---------------------------------------------------------------------------------
 * Protected functions to implement AliHLTComponent's interface.
 * These functions provide initialization as well as the actual processing
 * capabilities of the component. 
 * ---------------------------------------------------------------------------------
 */

// #################################################################################
Int_t AliHLTGlobalCompareFlatComponent::DoInit( Int_t argc, const Char_t** argv ) {
  // see header file for class documentation
  printf("AliHLTGlobalCompareFlatComponent::DoInit\n");
  // see header file for class documentation
  int iResult=0;

/*
	
		Int_t bins[fDim] = {3};
		Double_t mins[fDim] = {0};
		Double_t maxs[fDim] = {2};
		fhDiff = new THnSparseD("Differences","Differences",fDim,bins,mins,maxs);
		
		
		Int_t tmp = 0;
		
		fhDiff->GetAxis(tmp)->SetName("Overall");
		fhDiff->GetAxis(tmp)->SetBinLabel(1,"no differences");
		fhDiff->GetAxis(tmp)->SetBinLabel(2,"sizes differ");
		fhDiff->GetAxis(tmp)->SetBinLabel(3,"other differences");
		
		fhDiff->GetAxis(++tmp)->SetName("GetSize");
		fhDiff->GetAxis(++tmp)->SetName("GetMagneticField");
		fhDiff->GetAxis(++tmp)->SetName("GetPeriodNumber");
		fhDiff->GetAxis(++tmp)->SetName("GetRunNumber");
		fhDiff->GetAxis(++tmp)->SetName("GetOrbitNumber");
		fhDiff->GetAxis(++tmp)->SetName("GetBunchCrossNumber");
		fhDiff->GetAxis(++tmp)->SetName("GetTriggerMask");
		fhDiff->GetAxis(++tmp)->SetName("GetTriggerMaskNext50");
		fhDiff->GetAxis(++tmp)->SetName("GetFiredTriggerClasses");
		fhDiff->GetAxis(++tmp)->SetName("GetNumberOfTracks");
		fhDiff->GetAxis(++tmp)->SetName("GetNumberOfV0s");
		fhDiff->GetAxis(++tmp)->SetName("GetTimeStamp");
		fhDiff->GetAxis(++tmp)->SetName("GetEventSpecie");
	*/



		
		
  return iResult;
}



// #################################################################################
Int_t AliHLTGlobalCompareFlatComponent::DoDeinit() {
  // see header file for class documentation
  printf("AliHLTGlobalCompareFlatComponent::DoDeInit\n");

	Int_t iResult = 0;
	/*
	TFile * f = TFile::Open("histograms.root","RECREATE");
	f->Add(fhDiff);
	f->Write();
	f->Close();
	*/
  //delete fhDiff;
	
  return iResult;
}

// #################################################################################
Int_t AliHLTGlobalCompareFlatComponent::DoEvent(const AliHLTComponentEventData& /*evtData*/,
						    const AliHLTComponentBlockData* /*blocks*/, 
						    AliHLTComponentTriggerData& /*trigData*/,
						    AliHLTUInt8_t* /*outputPtr*/, 
						    AliHLTUInt32_t& /*size*/,
						    AliHLTComponentBlockDataList& /*outputBlocks*/) {
  // see header file for class documentation

  printf("AliHLTGlobalCompareFlatComponent::DoEvent\n");
  Int_t iResult=0;
	
	
  // -- Only use data event
 if (!IsDataEvent()) 
   return 0;

 AliFlatESDEvent *flatEsd[2] ;
	
  printf("search for input onbjects\n");
	{
	 int i=0;
	for ( const AliHLTComponentBlockData* pBlock=GetFirstInputBlock(kAliHLTDataTypeFlatESD|kAliHLTDataOriginOut);
    pBlock!=NULL && i<2; pBlock = GetNextInputBlock(),i++ ) {
			flatEsd[i] = reinterpret_cast<AliFlatESDEvent*>( pBlock->fPtr );
  }
	}
 cout<<"size event 1: "<<flatEsd[0]->GetSize()<<endl;
 cout<<"size event 2: "<<flatEsd[1]->GetSize()<<endl;


 cout<<"nTracks event 1: "<<flatEsd[0]->GetNumberOfTracks()<<endl;
 cout<<"nTracks event 2: "<<flatEsd[1]->GetNumberOfTracks()<<endl;
 
 outFile.open("comparison.txt",ios::app);
 outFile<<"\n\n------------------\nnew event\n------------------\n";
 printDiff( "AliFlatESDEvent::GetSize" ,flatEsd[0]->GetSize(), flatEsd[1]->GetSize() ) ;
	printDiff( "AliFlatESDEvent::GetMagneticField",flatEsd[0]->GetMagneticField(),flatEsd[1]->GetMagneticField() );
	printDiff( "AliFlatESDEvent::GetPeriodNumber",flatEsd[0]->GetPeriodNumber(),flatEsd[1]->GetPeriodNumber() );
	printDiff( "AliFlatESDEvent::GetRunNumber",flatEsd[0]->GetRunNumber(),	flatEsd[1]->GetRunNumber() );
	printDiff( "AliFlatESDEvent::GetOrbitNumber",flatEsd[0]->GetOrbitNumber(),flatEsd[1]->GetOrbitNumber() );
	printDiff( "AliFlatESDEvent::GetBunchCrossNumber",flatEsd[0]->GetBunchCrossNumber(),flatEsd[1]->GetBunchCrossNumber() );
	printDiff( "AliFlatESDEvent::GetTriggerMask",flatEsd[0]->GetTriggerMask(),flatEsd[1]->GetTriggerMask() );
	printDiff( "AliFlatESDEvent::GetTriggerMaskNext50",flatEsd[0]->GetTriggerMaskNext50(),flatEsd[1]->GetTriggerMaskNext50() );
	printDiff( "AliFlatESDEvent::GetFiredTriggerClasses",flatEsd[0]->GetFiredTriggerClasses() ,flatEsd[1]->GetFiredTriggerClasses() );
	printDiff( "AliFlatESDEvent::GetNumberOfTracks",flatEsd[0]->GetNumberOfTracks(),	flatEsd[1]->GetNumberOfTracks() );
	printDiff( "AliFlatESDEvent::GetNumberOfV0s",flatEsd[0]->GetNumberOfV0s(),flatEsd[1]->GetNumberOfV0s() );
	printDiff( "AliFlatESDEvent::GetTimeStamp",flatEsd[0]->GetTimeStamp(),flatEsd[1]->GetTimeStamp() );
	printDiff( "AliFlatESDEvent::GetEventSpecie",flatEsd[0]->GetEventSpecie(),flatEsd[1]->GetEventSpecie() ); 
	printDiff( "AliFlatESDEvent::GetNumberOfTriggerClasses",flatEsd[0]->GetNumberOfTriggerClasses(),flatEsd[1]->GetNumberOfTriggerClasses() ); 
	
	if(flatEsd[0]->GetNumberOfTriggerClasses()  && flatEsd[1]->GetNumberOfTriggerClasses() ){
		outFile<<"------------------\ntriggers\n------------------\n";
    AliFlatESDTrigger * trigger[2] = { const_cast<AliFlatESDTrigger*>(flatEsd[0]->GetTriggerClasses() ) , const_cast<AliFlatESDTrigger*>(flatEsd[1]->GetTriggerClasses() ) };
    for( Int_t i = 0; i < flatEsd[0]->GetNumberOfTriggerClasses()  && i < flatEsd[1]->GetNumberOfTriggerClasses()  ; i++ ){
      outFile<<"\nnew trigger\n";
			printDiff( "AliFlatESDTrigger::GetSize",trigger[0]->GetSize(),trigger[1]->GetSize() ); 
			printDiff( "AliFlatESDTrigger::GetTriggerIndex",trigger[0]->GetTriggerIndex(),trigger[1]->GetTriggerIndex() ); 
			printDiff( "AliFlatESDTrigger::GetTriggerClassName",trigger[0]->GetTriggerClassName(),trigger[1]->GetTriggerClassName() ); 
			
      trigger[0] = trigger[0]->GetNextTriggerNonConst();
			trigger[1] = trigger[1]->GetNextTriggerNonConst();
    }
	}
	if(flatEsd[0]->GetNumberOfV0s()  && flatEsd[1]->GetNumberOfV0s() ){
		outFile<<"------------------\nv0s\n------------------\n";
		
    AliFlatESDV0 * v0[2] = { const_cast<AliFlatESDV0*>(flatEsd[0]->GetV0s() ) , const_cast<AliFlatESDV0*>(flatEsd[1]->GetV0s() ) };
    for( Int_t i = 0; i < flatEsd[0]->GetNumberOfV0s()  && i < flatEsd[1]->GetNumberOfV0s()  ; i++ ){
      outFile<<"\nnew v0\n";
			printDiff( "AliFlatESDV0::GetSize",v0[0]->GetSize(),v0[1]->GetSize() ); 
			printDiff( "AliFlatESDV0::GetNegTrackID",v0[0]->GetNegTrackID(),v0[1]->GetNegTrackID() ); 
			printDiff( "AliFlatESDV0::GetPosTrackID",v0[0]->GetPosTrackID(),v0[1]->GetPosTrackID() ); 
			
      v0[0] = v0[0]->GetNextV0NonConst();
			v0[1] = v0[1]->GetNextV0NonConst();
    }
	}
	
	
	
	
 outFile.close();
 
 
  return iResult;
}


// #################################################################################
Int_t AliHLTGlobalCompareFlatComponent::ReadPreprocessorValues(const Char_t* /*modules*/) {
  // see header file for class documentation
  ALIHLTERRORGUARD(5, "ReadPreProcessorValues not implemented for this component");
  return 0;
}


int AliHLTGlobalCompareFlatComponent::Configure(const char*/* arguments*/)
{
  // see header file for class documentation
  int iResult=0;

  return iResult;
}

int AliHLTGlobalCompareFlatComponent::Reconfigure(const char* /*cdbEntry*/, const char* /*chainId*/)
{
  // see header file for class documentation
  int iResult=0;
  
  return iResult;
}


