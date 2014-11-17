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
#include "AliFlatESDFriend.h"
#include "AliFlatESDTrigger.h"
#include "AliFlatESDV0.h"
#include "AliFlatESDVertex.h"
#include "AliHLTErrorGuard.h"
#include "AliHLTDataTypes.h"
#include "AliHLTGlobalCompareFlatComponent.h"
#include "AliHLTITSClusterDataFormat.h"
#include "AliHLTTPCDefinitions.h"
#include "AliFlatTPCseed.h"
#include "AliExternalTrackParam.h"
#include "TTree.h"
#include "AliCDBEntry.h"
#include "AliCDBManager.h"

using namespace std;

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTGlobalCompareFlatComponent)

void AliHLTGlobalCompareFlatComponent::printDiff( string name, double val1, double val2){
	double sum = fabs(val1) + fabs(val2);
	double relDiff = ( val1 != 0 || val2!=0 ) ? (val1-val2)/sum: 0;
	
	int diff = 0;
	if (relDiff > 1e-3 && sum > 1e-6) diff = 1;
	else if(relDiff < -1e-3 && sum > 1e-6) diff = -1;
	outFile<<name<<"\t" << val1 << "\t" << val2 <<"\t" << diff << "\n";
}



void AliHLTGlobalCompareFlatComponent::printDiff( string name, int n , Float_t* vals1, Float_t* vals2 ){
	Float_t sum = 0;
	Float_t relDiff = 0;
	int diff = 0;
	for(int i=0; i<n && diff == 0; i++){
		sum = fabs(vals1[i]) + fabs(vals2[i]) ; 
		relDiff = ( vals1[i] != 0 || vals2[i] !=0 ) ? (vals1[i]-vals2[i])/sum : 0;
		if (relDiff > 1e-3 && sum > 1e-6) diff = 1;
		else if(relDiff < -1e-3 && sum > 1e-6) diff = -1;
	}
		
	outFile<<name<<"\t";
	for(int i=0; i<n;i++){
			outFile<<vals1[i]<<" ";
	}
	outFile<<"\t";
	for(int i=0; i<n;i++){
			outFile<<vals2[i]<<" ";
	}
	outFile<<"\t" << diff << "\n";
}

void AliHLTGlobalCompareFlatComponent::printDiff( string name, int n , Double_t* vals1, Double_t* vals2 ){
	Double_t sum = 0;
	Double_t relDiff = 0;
	int diff = 0;
	
	for(int i=0; i<n && diff == 0; i++){
		sum = fabs(vals1[i]) + fabs(vals2[i]) ; 
		relDiff = ( vals1[i] != 0 || vals2[i] !=0 ) ? (vals1[i]-vals2[i])/sum : 0;
		if (relDiff > 1e-3 && sum > 1e-6) diff = 1;
		else if(relDiff < -1e-3 && sum > 1e-6) diff = -1;
	}
		
	outFile<<name<<"\t";
	for(int i=0; i<n;i++){
			outFile<<vals1[i]<<" ";
	}
	outFile<<"\t";
	for(int i=0; i<n;i++){
			outFile<<vals2[i]<<" ";
	}
	outFile<<"\t" << diff << "\n";
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

 AliFlatESDEvent *flatEsd[2] ={0,0};
 AliFlatESDFriend *flatFriend[2] ={0,0};
	
  printf("search for input onbjects\n");
	{
	 int i=0;
	for ( const AliHLTComponentBlockData* pBlock=GetFirstInputBlock(kAliHLTDataTypeFlatESD|kAliHLTDataOriginOut);
    pBlock!=NULL && i<2; pBlock = GetNextInputBlock(),i++ ) {
			flatEsd[i] = reinterpret_cast<AliFlatESDEvent*>( pBlock->fPtr );
  }
  i=0;
	for ( const AliHLTComponentBlockData* pBlock=GetFirstInputBlock(kAliHLTDataTypeFlatESDFriend|kAliHLTDataOriginOut);
    pBlock!=NULL && i<2; pBlock = GetNextInputBlock(),i++ ) {
			flatFriend[i] = reinterpret_cast<AliFlatESDFriend*>( pBlock->fPtr );
  }
  
	}
	
	if( flatEsd[0] == NULL ){
			cout<<"flatEsd[0] not set!"<<endl;
			return 0;
	}
	if( flatEsd[1] == NULL ){
			cout<<"flatEsd[1] not set!"<<endl;
			return 0;
	}
	
	if( flatFriend[0] == NULL  ){
			cout<<"flatFriend[0] not set!"<<endl;
			return 0;
	}
	
	if( flatFriend[1] == NULL  ){
			cout<<"flatFriend[1] not set!"<<endl;
			return 0;
	}
	
	
 cout<<"size event : "<<flatEsd[0]->GetSize() << " "<<flatEsd[1]->GetSize()<<endl;
 cout<<"nTracks : "<<flatEsd[0]->GetNumberOfTracks()<<" "<<flatEsd[1]->GetNumberOfTracks()<<endl;
 cout<<"nV0s : "<<flatEsd[0]->GetNumberOfV0s()<<" "<<flatEsd[1]->GetNumberOfV0s()<<endl;
 
 cout<<"size friend : "<<flatFriend[0]->GetSize() << " "<<flatFriend[1]->GetSize()<<endl;
 cout<<"nFriendTracks : "<<flatFriend[0]->GetNumberOfTracks()<<" "<<flatFriend[1]->GetNumberOfTracks()<<endl;
 
 outFile.open("comparison.txt",ios::app);
 
 // Compare Event variables
 
 outFile<<"\n\n------------------\nnew AliFlatESDEvent\n------------------\n";
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
	
 const AliFlatESDVertex * vertexTracks[2] = {flatEsd[0]->GetFlatPrimaryVertexTracks(), flatEsd[1]->GetFlatPrimaryVertexTracks()};
	printDiff("AliFlatESDEvent::GetFlatPrimaryVertexTracks", (vertexTracks[0] ? 1:0), (vertexTracks[1] ? 1:0) );
	
 const AliFlatESDVertex * vertexSPD[2] = {flatEsd[0]->GetFlatPrimaryVertexSPD(), flatEsd[1]->GetFlatPrimaryVertexSPD()};
	printDiff("AliFlatESDEvent::GetFlatPrimaryVertexSPD", (vertexSPD[0] ? 1:0), (vertexSPD[1] ? 1:0) );
	
 // Compare primary vertices
	
	if(vertexTracks[0] && vertexTracks[1]){
      outFile<<"\nnew AliFlatESDVertexTracks\n";
			printDiff( "AliFlatESDVertexTracks::GetSize",vertexTracks[0]->GetSize(),vertexTracks[1]->GetSize() ); 
			printDiff( "AliFlatESDVertexTracks::GetX",vertexTracks[0]->GetX(),vertexTracks[1]->GetX() ); 
			printDiff( "AliFlatESDVertexTracks::GetY",vertexTracks[0]->GetY(),vertexTracks[1]->GetY() ); 
			printDiff( "AliFlatESDVertexTracks::GetZ",vertexTracks[0]->GetZ(),vertexTracks[1]->GetZ() ); 
	}
 
	if(vertexSPD[0] && vertexSPD[1]){
      outFile<<"\nnew AliFlatESDVertexSPD\n";
			outFile<<"AliFlatESDVertexSPD::memcmp "<<memcmp(vertexSPD[0], vertexSPD[1], max(vertexSPD[0]->GetSize(), vertexSPD[1]->GetSize()))<<"\n";
			printDiff( "AliFlatESDVertexSPD::GetSize",vertexSPD[0]->GetSize(),vertexSPD[1]->GetSize() ); 
			printDiff( "AliFlatESDVertexSPD::GetX",vertexSPD[0]->GetX(),vertexSPD[1]->GetX() ); 
			printDiff( "AliFlatESDVertexSPD::GetY",vertexSPD[0]->GetY(),vertexSPD[1]->GetY() ); 
			printDiff( "AliFlatESDVertexSPD::GetZ",vertexSPD[0]->GetZ(),vertexSPD[1]->GetZ() ); 
	}
 
 
 // Compare triggers
	
	if(flatEsd[0]->GetNumberOfTriggerClasses()  && flatEsd[1]->GetNumberOfTriggerClasses() ){
		outFile<<"------------------\ntriggers\n------------------\n";
    AliFlatESDTrigger * trigger[2] = { const_cast<AliFlatESDTrigger*>(flatEsd[0]->GetTriggerClasses() ) , const_cast<AliFlatESDTrigger*>(flatEsd[1]->GetTriggerClasses() ) };
    for( Int_t i = 0; i < flatEsd[0]->GetNumberOfTriggerClasses()  && i < flatEsd[1]->GetNumberOfTriggerClasses()  ; i++ ){
      outFile<<"\nnew AliFlatESDTrigger\n";
			printDiff( "AliFlatESDTrigger::GetSize",trigger[0]->GetSize(),trigger[1]->GetSize() ); 
			printDiff( "AliFlatESDTrigger::GetTriggerIndex",trigger[0]->GetTriggerIndex(),trigger[1]->GetTriggerIndex() ); 
			printDiff( "AliFlatESDTrigger::GetTriggerClassName",trigger[0]->GetTriggerClassName(),trigger[1]->GetTriggerClassName() ); 
			
      trigger[0] = trigger[0]->GetNextTriggerNonConst();
			trigger[1] = trigger[1]->GetNextTriggerNonConst();
    }
	}
	
 // Compare v0s
	
	if(flatEsd[0]->GetNumberOfV0s()  && flatEsd[1]->GetNumberOfV0s() ){
		outFile<<"------------------\nv0s\n------------------\n";
		
    AliFlatESDV0 * v0[2] = { const_cast<AliFlatESDV0*>(flatEsd[0]->GetV0s() ) , const_cast<AliFlatESDV0*>(flatEsd[1]->GetV0s() ) };
    for( Int_t i = 0; i < flatEsd[0]->GetNumberOfV0s()  && i < flatEsd[1]->GetNumberOfV0s()  ; i++ ){
      outFile<<"\nnew AliFlatESDV0\n";
			printDiff( "AliFlatESDV0::GetSize",v0[0]->GetSize(),v0[1]->GetSize() ); 
			printDiff( "AliFlatESDV0::GetNegTrackID",v0[0]->GetNegTrackID(),v0[1]->GetNegTrackID() ); 
			printDiff( "AliFlatESDV0::GetPosTrackID",v0[0]->GetPosTrackID(),v0[1]->GetPosTrackID() ); 
			
      v0[0] = v0[0]->GetNextV0NonConst();
			v0[1] = v0[1]->GetNextV0NonConst();
    }
	}
	
 // Compare tracks
	
	if(flatEsd[0]->GetNumberOfTracks()  && flatEsd[1]->GetNumberOfTracks() ){
		outFile<<"------------------\ntracks\n------------------\n";
		
    AliFlatESDTrack * track[2] = { const_cast<AliFlatESDTrack*>(flatEsd[0]->GetTracks() ) , const_cast<AliFlatESDTrack*>(flatEsd[1]->GetTracks() ) };
    for( Int_t t = 0; t < flatEsd[0]->GetNumberOfTracks()  && t < flatEsd[1]->GetNumberOfTracks()  ; t++ ){
      outFile<<"\nnew AliFlatESDTrack\n";
			printDiff( "AliFlatESDTrack::GetSize",track[0]->GetSize(),track[1]->GetSize() ); 
			printDiff( "AliFlatESDTrack::GetNumberOfTPCClusters",track[0]->GetNumberOfTPCClusters(),track[1]->GetNumberOfTPCClusters() ); 
			printDiff( "AliFlatESDTrack::GetNumberOfITSClusters",track[0]->GetNumberOfITSClusters(),track[1]->GetNumberOfITSClusters() ); 
			
			const char* pNames[7] = {"", "Refitted", "Ip", "TPCInner", "Op", "Cp", "ITSOUT"};
			
			const AliFlatExternalTrackParam * p[7][2] = {
				{track[0]->GetFlatTrackParam(), track[1]->GetFlatTrackParam()},
				{track[0]->GetFlatTrackParamRefitted(), track[1]->GetFlatTrackParamRefitted()},
				{track[0]->GetFlatTrackParamIp(), track[1]->GetFlatTrackParamIp()},
				{track[0]->GetFlatTrackParamTPCInner(), track[1]->GetFlatTrackParamTPCInner()},
				{track[0]->GetFlatTrackParamOp(), track[1]->GetFlatTrackParamOp()},
				{track[0]->GetFlatTrackParamCp(), track[1]->GetFlatTrackParamCp()},
				{track[0]->GetFlatTrackParamITSOut(), track[1]->GetFlatTrackParamITSOut()}
			};
			
			for(int i = 0 ; i<7; i++){
				printDiff( Form("AliFlatESDTrack::GetFlatTrackParam%s",pNames[i]), (p[i][0] ? 1:0), (p[i][1] ? 1:0) );
			}

			for(int i = 0 ; i<7 && p[i][0] && p[i][1]; i++){
				outFile<<"\nnew AliFlatExternalTrackParam" << pNames[i] << "\n";
				printDiff( Form("AliFlatExternalTrackParam%s::GetAlpha",pNames[i]),p[i][0]->GetAlpha(),p[i][1]->GetAlpha() ); 
				printDiff( Form("AliFlatExternalTrackParam%s::GetX",pNames[i]),p[i][0]->GetX(),p[i][1]->GetX() ); 
				printDiff( Form("AliFlatExternalTrackParam%s::GetY",pNames[i]),p[i][0]->GetY(),p[i][1]->GetY() ); 
				printDiff( Form("AliFlatExternalTrackParam%s::GetZ",pNames[i]),p[i][0]->GetZ(),p[i][1]->GetZ() ); 
				printDiff( Form("AliFlatExternalTrackParam%s::GetSnp",pNames[i]),p[i][0]->GetSnp(),p[i][1]->GetSnp() ); 
				printDiff( Form("AliFlatExternalTrackParam%s::GetTgl",pNames[i]),p[i][0]->GetTgl(),p[i][1]->GetTgl() ); 
				printDiff( Form("AliFlatExternalTrackParam%s::GetSigned1Pt",pNames[i]),p[i][0]->GetSigned1Pt(),p[i][1]->GetSigned1Pt() ); 
				
				
				Float_t* cov[2] = { p[i][0]->GetCov() , p[i][1]->GetCov() };
				printDiff( Form("AliFlatExternalTrackParam%s::GetCov",pNames[i]) , 15, cov[0], cov[1]); 
			}
			
			
			
			
      track[0] = track[0]->GetNextTrackNonConst();
			track[1] = track[1]->GetNextTrackNonConst();			
			
    }
	}
	
	
 // Compare Friend variables
 
	outFile<<"\n\n------------------\nnew AliFlatESDFriend\n------------------\n";
	printDiff( "AliFlatESDFriend::GetSize" ,flatFriend[0]->GetSize(), flatFriend[1]->GetSize() ) ;
	printDiff( "AliFlatESDFriend::GetNumberOfTracks" ,flatFriend[0]->GetNumberOfTracks(), flatFriend[1]->GetNumberOfTracks());
	printDiff( "AliFlatESDFriend::GetEntriesInTracks" ,flatFriend[0]->GetEntriesInTracks(), flatFriend[1]->GetEntriesInTracks());
	printDiff( "AliFlatESDFriend::TestSkipBit" ,flatFriend[0]->TestSkipBit(), flatFriend[1]->TestSkipBit() ) ;
 
	Double_t clTpc[2][72];
	Double_t clTpcU[2][72];
	for(int i = 0 ; i<72; i++){
		clTpc[0][i] = flatFriend[0]->GetNclustersTPC(i);
		clTpc[1][i] = flatFriend[1]->GetNclustersTPC(i);
		
		clTpcU[0][i] = flatFriend[0]->GetNclustersTPCused(i);
		clTpcU[1][i] = flatFriend[1]->GetNclustersTPCused(i);
	}
	
	
	printDiff( "AliFlatESDFriend::GetNclustersTPC" ,72, clTpc[0], clTpc[1] ) ;
	printDiff( "AliFlatESDFriend::GetNclustersTPCused" ,72, clTpcU[0], clTpcU[1] ) ;
	
	
	// compare friend tracks
	
	if(flatFriend[0]->GetEntriesInTracks()  && flatFriend[1]->GetEntriesInTracks() ){
		outFile<<"------------------\nfriend tracks\n------------------\n";
		
    AliFlatESDFriendTrack * track[2] = { flatFriend[0]->GetFlatTrackEntryNonConst(0), flatFriend[1]->GetFlatTrackEntryNonConst(0)};
		
		
    for( Int_t t = 0; t < flatFriend[0]->GetEntriesInTracks()  && t < flatFriend[1]->GetEntriesInTracks()  ; t++ ){
			
			
		//	cout<<"track0.size"<<track[0]->GetSize()<<endl;
	//		cout<<"track1.size"<<track[1]->GetSize()<<endl;
			
			
			if(!track[0] || !track[1]) continue;
      outFile<<"\nnew AliFlatESDFriendTrack\n";
			printDiff( "AliFlatESDFriendTrack::GetSize",track[0]->GetSize(),track[1]->GetSize() ); 
			
			AliExternalTrackParam p[3][2]; 
			int pp[3][2]= {{0,0},{0,0},{0,0}}; 
			
			const char* pNames[3] = {"TPCOut", "ITSOut", "TRDIn"};
			
			pp[0][0] = track[0]->GetTrackParamTPCOut(p[0][0] ) >-1 ? 1: 0;
			pp[0][1] = track[1]->GetTrackParamTPCOut(p[0][1] ) >-1 ? 1: 0;
			printDiff( "AliFlatESDFriendTrack::GetTrackParamTPCOut",pp[0][0], pp[0][1] ); 
			
			pp[1][0] = track[0]->GetTrackParamITSOut(p[1][0] ) >-1 ? 1: 0;
			pp[1][1] = track[1]->GetTrackParamITSOut(p[1][1] ) >-1 ? 1: 0;
			printDiff( "AliFlatESDFriendTrack::GetTrackParamITSOut",pp[1][0], pp[1][1] ); 
			
			
			pp[2][0] = track[0]->GetTrackParamTRDIn(p[2][0] ) >-1 ? 1: 0;
			pp[2][1] = track[1]->GetTrackParamTRDIn(p[2][1] ) >-1 ? 1: 0;
			printDiff( "AliFlatESDFriendTrack::GetTrackParamTRDIn",pp[2][0], pp[2][1] ); 
			
			for(int i = 0 ; pp[0][i] && pp[1][i] && i<3; i++){
				outFile<<"\nnew AliExternalTrackParam" << pNames[i] << "\n";
				printDiff( Form("AliExternalTrackParam%s::GetAlpha",pNames[i]),p[i][0].GetAlpha(),p[i][1].GetAlpha() ); 
				printDiff( Form("AliExternalTrackParam%s::GetX",pNames[i]),p[i][0].GetX(),p[i][1].GetX() ); 
				printDiff( Form("AliExternalTrackParam%s::GetY",pNames[i]),p[i][0].GetY(),p[i][1].GetY() ); 
				printDiff( Form("AliExternalTrackParam%s::GetZ",pNames[i]),p[i][0].GetZ(),p[i][1].GetZ() ); 
				printDiff( Form("AliExternalTrackParam%s::GetSnp",pNames[i]),p[i][0].GetSnp(),p[i][1].GetSnp() ); 
				printDiff( Form("AliExternalTrackParam%s::GetTgl",pNames[i]),p[i][0].GetTgl(),p[i][1].GetTgl() ); 
				printDiff( Form("AliExternalTrackParam%s::GetSigned1Pt",pNames[i]),p[i][0].GetSigned1Pt(),p[i][1].GetSigned1Pt() ); 
					
					
				Double_t* cov[2] = { const_cast<Double_t*>( p[i][0].GetCovariance()) , const_cast<Double_t*>( p[i][1].GetCovariance() ) };
				printDiff( Form("AliExternalTrackParam%s::GetCovariance",pNames[i]) , 15, cov[0], cov[1]); 
			
			}
			
			
			
			
			
			const AliFlatTPCseed* s[2]={ track[0]->GetFlatTPCseed(), track[1]->GetFlatTPCseed()};
			printDiff( "AliFlatTPCseed::GetSize",s[0]->GetSize(),s[1]->GetSize() ); 
			printDiff( "AliFlatTPCseed::GetLabel",s[0]->GetLabel(),s[1]->GetLabel() ); 
			printDiff( "AliFlatTPCseed::GetNClusters",s[0]->GetNClusters(),s[1]->GetNClusters() ); 
			
			
			//printf("track0: %p next: %p", track[0], track[0]->GetNextTrackNonConst() );
			//printf("track1: %p next: %p", track[1], track[1]->GetNextTrackNonConst() );
			
      track[0] = track[0]->GetNextTrackNonConst();
			track[1] = track[1]->GetNextTrackNonConst();	
			
			
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


