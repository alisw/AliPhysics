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
#include "AliFlatTPCCluster.h"
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
	if(diff!=0){
		conflictsFile<<fCurrentClass<<"\t"<<name<<"\t" << val1 << "\t" << val2  << "\n";
	}
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
	if(diff!=0){
		conflictsFile<<fCurrentClass<<"\t"<<name<< "\n";
	}
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
	if(diff!=0){
		conflictsFile<<fCurrentClass<<"\t"<<name<< "\n";
	}
}




void AliHLTGlobalCompareFlatComponent::printDiff( string name, TString val1, TString val2){
	outFile << name << "\t" << "\t\"" << val1 <<"\"\t\"" << val2 <<"\"\t" << (val1.EqualTo(val2) ?0:1)<<"\n";
	if(! val1.EqualTo(val2) ){
		conflictsFile<<fCurrentClass<<"\t"<<name<<"\t" << val1 << "\t" << val2 << "\n";
	}
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
  return iResult;
}



// #################################################################################
Int_t AliHLTGlobalCompareFlatComponent::DoDeinit() {
  // see header file for class documentation
	Int_t iResult = 0;
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

  Int_t iResult=0;
	
	
  // -- Only use data event
 if (!IsDataEvent()) 
   return 0;

 AliFlatESDEvent *flatEsd[2] ={0,0};
 AliFlatESDFriend *flatFriend[2] ={0,0};
	
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
			HLTError("flatEsd[0] not set!");
			return 0;
	}
	if( flatEsd[1] == NULL ){
			HLTError("flatEsd[1] not set!");
			return 0;
	}
	
	if( flatFriend[0] == NULL  ){
			HLTError("flatFriend[0] not set!");
			return 0;
	}
	
	if( flatFriend[1] == NULL  ){
			HLTError("flatFriend[1] not set!");
			return 0;
	}
	
 outFile.open("comparison.txt",ios::app);
 
 conflictsFile.open("conflicts.txt",ios::app);
 
 // Compare Event variables
 
 outFile<<"_FlatESDEvent\n";
 fCurrentClass = "FlatESDEvent";
 printDiff( "GetSize" ,flatEsd[0]->GetSize(), flatEsd[1]->GetSize() ) ;
	printDiff( "GetMagneticField",flatEsd[0]->GetMagneticField(),flatEsd[1]->GetMagneticField() );
	printDiff( "GetPeriodNumber",flatEsd[0]->GetPeriodNumber(),flatEsd[1]->GetPeriodNumber() );
	printDiff( "GetRunNumber",flatEsd[0]->GetRunNumber(),	flatEsd[1]->GetRunNumber() );
	printDiff( "GetOrbitNumber",flatEsd[0]->GetOrbitNumber(),flatEsd[1]->GetOrbitNumber() );
	printDiff( "GetBunchCrossNumber",flatEsd[0]->GetBunchCrossNumber(),flatEsd[1]->GetBunchCrossNumber() );
	printDiff( "GetTriggerMask",flatEsd[0]->GetTriggerMask(),flatEsd[1]->GetTriggerMask() );
	printDiff( "GetTriggerMaskNext50",flatEsd[0]->GetTriggerMaskNext50(),flatEsd[1]->GetTriggerMaskNext50() );
	printDiff( "GetFiredTriggerClasses",flatEsd[0]->GetFiredTriggerClasses() ,flatEsd[1]->GetFiredTriggerClasses() );
	printDiff( "GetNumberOfTracks",flatEsd[0]->GetNumberOfTracks(),	flatEsd[1]->GetNumberOfTracks() );
	printDiff( "GetNumberOfV0s",flatEsd[0]->GetNumberOfV0s(),flatEsd[1]->GetNumberOfV0s() );
	printDiff( "GetTimeStamp",flatEsd[0]->GetTimeStamp(),flatEsd[1]->GetTimeStamp() );
	printDiff( "GetEventSpecie",flatEsd[0]->GetEventSpecie(),flatEsd[1]->GetEventSpecie() ); 
	printDiff( "GetNumberOfTriggerClasses",flatEsd[0]->GetNumberOfTriggerClasses(),flatEsd[1]->GetNumberOfTriggerClasses() ); 
	
 const AliFlatESDVertex * vertexTracks[2] = {flatEsd[0]->GetFlatPrimaryVertexTracks(), flatEsd[1]->GetFlatPrimaryVertexTracks()};
	printDiff("GetFlatPrimaryVertexTracks", (vertexTracks[0] ? 1:0), (vertexTracks[1] ? 1:0) );
	
 const AliFlatESDVertex * vertexSPD[2] = {flatEsd[0]->GetFlatPrimaryVertexSPD(), flatEsd[1]->GetFlatPrimaryVertexSPD()};
	printDiff("GetFlatPrimaryVertexSPD", (vertexSPD[0] ? 1:0), (vertexSPD[1] ? 1:0) );
	
 // Compare primary vertices
	
	if(vertexTracks[0] && vertexTracks[1]){
      outFile<<"_FlatESDVertexTracks\n";
			fCurrentClass = "FlatESDVertexTracks";
			printDiff( "GetSize",vertexTracks[0]->GetSize(),vertexTracks[1]->GetSize() ); 
			printDiff( "GetX",vertexTracks[0]->GetX(),vertexTracks[1]->GetX() ); 
			printDiff( "GetY",vertexTracks[0]->GetY(),vertexTracks[1]->GetY() ); 
			printDiff( "GetZ",vertexTracks[0]->GetZ(),vertexTracks[1]->GetZ() ); 
	}
 
	if(vertexSPD[0] && vertexSPD[1]){
      outFile<<"_FlatESDVertexSPD\n";
			fCurrentClass = "FlatESDVertexSPD";
			printDiff( "GetSize",vertexSPD[0]->GetSize(),vertexSPD[1]->GetSize() ); 
			printDiff( "GetX",vertexSPD[0]->GetX(),vertexSPD[1]->GetX() ); 
			printDiff( "GetY",vertexSPD[0]->GetY(),vertexSPD[1]->GetY() ); 
			printDiff( "GetZ",vertexSPD[0]->GetZ(),vertexSPD[1]->GetZ() ); 
	}
 
 
 // Compare triggers
	
	if(flatEsd[0]->GetNumberOfTriggerClasses()  && flatEsd[1]->GetNumberOfTriggerClasses() ){
    AliFlatESDTrigger * trigger[2] = { const_cast<AliFlatESDTrigger*>(flatEsd[0]->GetTriggerClasses() ) , const_cast<AliFlatESDTrigger*>(flatEsd[1]->GetTriggerClasses() ) };
    for( Int_t i = 0; i < flatEsd[0]->GetNumberOfTriggerClasses()  && i < flatEsd[1]->GetNumberOfTriggerClasses()  ; i++ ){
      outFile<<"_FlatESDTrigger\n";
			fCurrentClass = "FlatESDTrigger";
			printDiff( "GetSize",trigger[0]->GetSize(),trigger[1]->GetSize() ); 
			printDiff( "GetTriggerIndex",trigger[0]->GetTriggerIndex(),trigger[1]->GetTriggerIndex() ); 
			printDiff( "GetTriggerClassName",trigger[0]->GetTriggerClassName(),trigger[1]->GetTriggerClassName() ); 
			
      trigger[0] = trigger[0]->GetNextTriggerNonConst();
			trigger[1] = trigger[1]->GetNextTriggerNonConst();
    }
	}
	
 // Compare v0s
	
	if(flatEsd[0]->GetNumberOfV0s()  && flatEsd[1]->GetNumberOfV0s() ){		
    AliFlatESDV0 * v0[2] = { const_cast<AliFlatESDV0*>(flatEsd[0]->GetV0s() ) , const_cast<AliFlatESDV0*>(flatEsd[1]->GetV0s() ) };
    for( Int_t i = 0; i < flatEsd[0]->GetNumberOfV0s()  && i < flatEsd[1]->GetNumberOfV0s()  ; i++ ){
      outFile<<"_FlatESDV0\n";
			fCurrentClass = "FlatESDV0";
			printDiff( "GetSize",v0[0]->GetSize(),v0[1]->GetSize() ); 
			printDiff( "GetNegTrackID",v0[0]->GetNegTrackID(),v0[1]->GetNegTrackID() ); 
			printDiff( "GetPosTrackID",v0[0]->GetPosTrackID(),v0[1]->GetPosTrackID() ); 
			
      v0[0] = v0[0]->GetNextV0NonConst();
			v0[1] = v0[1]->GetNextV0NonConst();
    }
	}
	
 // Compare tracks
	
	if(flatEsd[0]->GetNumberOfTracks()  && flatEsd[1]->GetNumberOfTracks() ){		
    AliFlatESDTrack * track[2] = { const_cast<AliFlatESDTrack*>(flatEsd[0]->GetTracks() ) , const_cast<AliFlatESDTrack*>(flatEsd[1]->GetTracks() ) };
    for( Int_t t = 0; t < flatEsd[0]->GetNumberOfTracks()  && t < flatEsd[1]->GetNumberOfTracks()  ; t++ ){
      outFile<<"_FlatESDTrack\n";
			fCurrentClass = "FlatESDTrack";
			printDiff( "GetSize",track[0]->GetSize(),track[1]->GetSize() ); 
			printDiff( "GetNumberOfTPCClusters",track[0]->GetNumberOfTPCClusters(),track[1]->GetNumberOfTPCClusters() ); 
			printDiff( "GetNumberOfITSClusters",track[0]->GetNumberOfITSClusters(),track[1]->GetNumberOfITSClusters() ); 
			
			const char* pNames[6] =  {"", "Refitted", "Ip", "TPCInner", "Op", "Cp"};
			
			const AliFlatExternalTrackParam * p[7][2] = {
				{track[0]->GetFlatTrackParam(), track[1]->GetFlatTrackParam()},
				{track[0]->GetFlatTrackParamRefitted(), track[1]->GetFlatTrackParamRefitted()},
				{track[0]->GetFlatTrackParamIp(), track[1]->GetFlatTrackParamIp()},
				{track[0]->GetFlatTrackParamTPCInner(), track[1]->GetFlatTrackParamTPCInner()},
				{track[0]->GetFlatTrackParamOp(), track[1]->GetFlatTrackParamOp()},
				{track[0]->GetFlatTrackParamCp(), track[1]->GetFlatTrackParamCp()}
			};
			
			for(int i = 0 ; i<6; i++){
				printDiff( Form("GetFlatTrackParam%s",pNames[i]) ,(p[i][0] ? 1:0), (p[i][1] ? 1:0) );
			}

			for(int i = 0 ; i<6 ; i++){
				if(p[i][0] && p[i][1]){
				outFile<<"_FlatExternalTrackParam" << pNames[i] << "\n";
				fCurrentClass = "FlatExternalTrackParam";
				printDiff( "GetAlpha",p[i][0]->GetAlpha(),p[i][1]->GetAlpha() ); 
				printDiff( "GetX",p[i][0]->GetX(),p[i][1]->GetX() ); 
				printDiff( "GetY",p[i][0]->GetY(),p[i][1]->GetY() ); 
				printDiff( "GetZ",p[i][0]->GetZ(),p[i][1]->GetZ() ); 
				printDiff( "GetSnp",p[i][0]->GetSnp(),p[i][1]->GetSnp() ); 
				printDiff( "GetTgl",p[i][0]->GetTgl(),p[i][1]->GetTgl() ); 
				printDiff( "GetSigned1Pt",p[i][0]->GetSigned1Pt(),p[i][1]->GetSigned1Pt() ); 
				Float_t* cov[2] = { p[i][0]->GetCov() , p[i][1]->GetCov() };
				printDiff( "GetCov", 15, cov[0], cov[1]); 
				}
			}
      track[0] = track[0]->GetNextTrackNonConst();
			track[1] = track[1]->GetNextTrackNonConst();
    }
	}
	
	
 // Compare Friend variables
 
	outFile<<"_FlatESDFriend\n";
	fCurrentClass = "FlatESDFriend";
	printDiff( "GetSize" ,flatFriend[0]->GetSize(), flatFriend[1]->GetSize() ) ;
	printDiff( "GetNumberOfTracks" ,flatFriend[0]->GetNumberOfTracks(), flatFriend[1]->GetNumberOfTracks());
	printDiff( "GetEntriesInTracks" ,flatFriend[0]->GetEntriesInTracks(), flatFriend[1]->GetEntriesInTracks());
	printDiff( "TestSkipBit" ,flatFriend[0]->TestSkipBit(), flatFriend[1]->TestSkipBit() ) ;
 
	Double_t clTpc[2][72];
	Double_t clTpcU[2][72];
	for(int i = 0 ; i<72; i++){
		clTpc[0][i] = flatFriend[0]->GetNclustersTPC(i);
		clTpc[1][i] = flatFriend[1]->GetNclustersTPC(i);
		
		clTpcU[0][i] = flatFriend[0]->GetNclustersTPCused(i);
		clTpcU[1][i] = flatFriend[1]->GetNclustersTPCused(i);
	}
	
	
	printDiff( "GetNclustersTPC" ,72, clTpc[0], clTpc[1] ) ;
	printDiff( "GetNclustersTPCused" ,72, clTpcU[0], clTpcU[1] ) ;
	
	
	// compare friend tracks
	
	if(flatFriend[0]->GetEntriesInTracks()  && flatFriend[1]->GetEntriesInTracks() ){
		
	  AliFlatESDFriendTrack * track[2] = { flatFriend[0]->GetFirstTrackEntryNonConst(), flatFriend[1]->GetFirstTrackEntryNonConst()};
		
		
    for( Int_t t = 0; t < flatFriend[0]->GetEntriesInTracks()  && t < flatFriend[1]->GetEntriesInTracks()  ; t++ ){
			if(!track[0] || !track[1]) continue;
      outFile<<"_FlatESDFriendTrack\n";
			fCurrentClass = "FlatESDFriendTrack";
			printDiff( "GetSize",track[0]->GetSize(),track[1]->GetSize() ); 
			
			AliExternalTrackParam p[3][2]; 
			int pp[3][2]= {{0,0},{0,0},{0,0}}; 
			
			const char* pNames[3] = {"TPCOut", "ITSOut", "TRDIn"};
			
			pp[0][0] = track[0]->GetTrackParamTPCOut(p[0][0] ) >-1 ? 1: 0;
			pp[0][1] = track[1]->GetTrackParamTPCOut(p[0][1] ) >-1 ? 1: 0;
			printDiff( "GetTrackParamTPCOut",pp[0][0], pp[0][1] ); 
			
			pp[1][0] = track[0]->GetTrackParamITSOut(p[1][0] ) >-1 ? 1: 0;
			pp[1][1] = track[1]->GetTrackParamITSOut(p[1][1] ) >-1 ? 1: 0;
			printDiff( "GetTrackParamITSOut",pp[1][0], pp[1][1] ); 
			
			
			pp[2][0] = track[0]->GetTrackParamTRDIn(p[2][0] ) >-1 ? 1: 0;
			pp[2][1] = track[1]->GetTrackParamTRDIn(p[2][1] ) >-1 ? 1: 0;
			printDiff( "GetTrackParamTRDIn",pp[2][0], pp[2][1] ); 
			
 			for(int i = 0 ; i<3; i++){
				
				if(pp[i][0] && pp[i][1]){
					
				outFile<<"_ExternalTrackParam" << pNames[i] << "\n";
				fCurrentClass = "ExternalTrackParam";
				printDiff( "GetAlpha" ,p[i][0].GetAlpha(),p[i][1].GetAlpha() ); 
				printDiff( "GetX",p[i][0].GetX(),p[i][1].GetX() ); 
				printDiff( "GetY",p[i][0].GetY(),p[i][1].GetY() ); 
				printDiff( "GetZ",p[i][0].GetZ(),p[i][1].GetZ() ); 
				printDiff( "GetSnp",p[i][0].GetSnp(),p[i][1].GetSnp() ); 
				printDiff( "GetTgl",p[i][0].GetTgl(),p[i][1].GetTgl() ); 
				printDiff( "GetSigned1Pt",p[i][0].GetSigned1Pt(),p[i][1].GetSigned1Pt() ); 
					
					
				Double_t* cov[2] = { const_cast<Double_t*>( p[i][0].GetCovariance()) , const_cast<Double_t*>( p[i][1].GetCovariance() ) };
				printDiff( Form("ExternalTrackParam%s::GetCovariance",pNames[i]) , 15, cov[0], cov[1]); 
				}
			}
			
			
			
			const AliFlatTPCseed* s[2]={ track[0]->GetFlatTPCseed(), track[1]->GetFlatTPCseed()};
			if(s[0] && s[1]){
			
			
      outFile<<"_FlatTPCseed\n";
			fCurrentClass = "FlatTPCseed";
			
			printDiff( "GetSize",s[0]->GetSize(),s[1]->GetSize() ); 
			printDiff( "GetLabel",s[0]->GetLabel(),s[1]->GetLabel() ); 
			printDiff( "GetNClusters",s[0]->GetNClusters(),s[1]->GetNClusters() ); 
			
			printDiff( "GetSize",s[0]->GetSize(),s[1]->GetSize() ); 
			printDiff( "GetLabel",s[0]->GetLabel(),s[1]->GetLabel() ); 
			printDiff( "GetNClusters",s[0]->GetNClusters(),s[1]->GetNClusters() ); 
			
			
				
			// loop over clusters 
			if(s[0]->GetNClusters() == s[1]->GetNClusters()){
			int ncl = s[0]->GetNClusters();
			//	cout<<"number of clusters: "<<ncl<<endl;
			AliFlatTPCCluster* cl[72][160][2];
					
			for(int j=0; j<72;j++){
				for( int i=0; i<160; i++ ) {
					cl[j][i][0]=0;
					cl[j][i][1]=0;
				}
			}
			for(int icl=0; icl < ncl; icl++){
				//		if(cl[ (s[0]->GetClusters()[icl]).GetPadRow() ][0])  cout<<"ERROR: cluster " << icl <<" [0] already set!!!"<<endl;
					//	if(cl[ (s[0]->GetClusters()[icl]).GetPadRow() ][1])  cout<<"ERROR: cluster " << icl <<" [1] already set!!!"<<endl;
						cl[ (s[0]->GetClusters()[icl]).GetSector() ][ (s[0]->GetClusters()[icl]).GetPadRow() ][0] = const_cast<AliFlatTPCCluster*>( &(s[0]->GetClusters()[icl])  );
						cl[ (s[1]->GetClusters()[icl]).GetSector() ][ (s[1]->GetClusters()[icl]).GetPadRow() ][1] = const_cast<AliFlatTPCCluster*>( &(s[1]->GetClusters()[icl])  );
					}
					for(int iSector=0; iSector<72; iSector++){
					
						for(int irow= 0; irow<160;irow++){
							if( cl[iSector][irow][0] && cl[iSector][irow][1] ){
								outFile<<"_FlatTPCCluster\n";
								fCurrentClass = "FlatTPCCluster";
								printDiff( "GetX",cl[iSector][irow][0]->GetX(),cl[iSector][irow][1]->GetX() ); 
								printDiff( "GetY",cl[iSector][irow][0]->GetY(),cl[iSector][irow][1]->GetY() ); 
								printDiff( "GetZ",cl[iSector][irow][0]->GetZ(),cl[iSector][irow][1]->GetZ() ); 
								printDiff( "GetSector",cl[iSector][irow][0]->GetSector(),cl[iSector][irow][1]->GetSector() ); 
								printDiff( "GetPadRow",cl[iSector][irow][0]->GetPadRow(),cl[iSector][irow][1]->GetPadRow() ); 
								printDiff( "GetSigmaY2",cl[iSector][irow][0]->GetSigmaY2(),cl[iSector][irow][1]->GetSigmaY2() ); 
								printDiff( "GetSigmaZ2",cl[iSector][irow][0]->GetSigmaZ2(),cl[iSector][irow][1]->GetSigmaZ2() ); 
								printDiff( "GetCharge",cl[iSector][irow][0]->GetCharge(),cl[iSector][irow][1]->GetCharge() ); 
								printDiff( "GetQMax",cl[iSector][irow][0]->GetQMax(),cl[iSector][irow][1]->GetQMax() ); 
								printDiff( "GetTrackAngleY",cl[iSector][irow][0]->GetTrackAngleY(),cl[iSector][irow][1]->GetTrackAngleY() ); 
								printDiff( "GetTrackAngleZ",cl[iSector][irow][0]->GetTrackAngleZ(),cl[iSector][irow][1]->GetTrackAngleZ() ); 
							}
							else if( cl[iSector][irow][0] || cl[iSector][irow][1] ){
								printDiff( "GetClusters(i)",  cl[iSector][irow][0] ?1:0 ,cl[iSector][irow][1] ?1:0 ); 
							}
						}
					}
					
				}
			
			}
      track[0] = track[0]->GetNextTrackNonConst();
			track[1] = track[1]->GetNextTrackNonConst();	
		}
	
	}
	
	outFile.close();
	conflictsFile.close();
 
 
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


