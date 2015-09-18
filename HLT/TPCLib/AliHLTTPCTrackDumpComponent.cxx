// $Id$

//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Gaute Ovrebekk <ovrebekk@ift.uib.no>                  *
//*                  for The ALICE HLT Project.                            *
//*                                                                        *
//* Permission to use, copy, modify and distribute this software and its   *
//* documentation strictly for non-commercial purposes is hereby granted   *
//* without fee, provided that the above copyright notice appears in all   *
//* copies and that both the copyright notice and this permission notice   *
//* appear in the supporting documentation. The authors make no claims     *
//* about the suitability of this software for any purpose. It is          *
//* provided "as is" without express or implied warranty.                  *
//**************************************************************************

/** @file   AliHLTTPCTrackDumpComponent.cxx
    @author Gaute Ovrebekk
    @date   
    @brief  Special file writer converting TPC tracks input to ASCII. */

#include <cassert>
#include "AliHLTTPCTrackDumpComponent.h"
#include "AliHLTTPCTrackletDataFormat.h"
#include "AliHLTTPCDefinitions.h"
#include <TSystem.h>

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTTPCTrackDumpComponent)

AliHLTTPCTrackDumpComponent::AliHLTTPCTrackDumpComponent()
  :
AliHLTFileWriter()
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

AliHLTTPCTrackDumpComponent::~AliHLTTPCTrackDumpComponent()
{
  // see header file for class documentation
}

const char* AliHLTTPCTrackDumpComponent::GetComponentID()
{
  // see header file for class documentation
  return "TPCTrackDump";
}

void AliHLTTPCTrackDumpComponent::GetInputDataTypes( vector<AliHLTComponentDataType>& list)
{
  // see header file for class documentation
  list.clear();
  list.push_back(AliHLTTPCDefinitions::fgkTrackSegmentsDataType);
  list.push_back(AliHLTTPCDefinitions::fgkTracksDataType);

}

AliHLTComponent* AliHLTTPCTrackDumpComponent::Spawn()
{
  // see header file for class documentation
  return new AliHLTTPCTrackDumpComponent;
}

int AliHLTTPCTrackDumpComponent::InitWriter()
{
  // see header file for class documentation
  return 0;
}

int AliHLTTPCTrackDumpComponent::ScanArgument(int /*argc*/, const char** /*argv*/)
{
  // see header file for class documentation
  int iResult=0;
  TString argument="";
  bool bMissingParam=0;
  int i=0;
  
  if (bMissingParam) iResult=-EPROTO;
  else if (iResult>=0) iResult=i;

  return iResult;
}

int AliHLTTPCTrackDumpComponent::CloseWriter()
{
  // see header file for class documentation
  return 0;
}

int AliHLTTPCTrackDumpComponent::DumpEvent( const AliHLTComponentEventData& evtData,
					    AliHLTComponentTriggerData& /*trigData*/ )
{
  // see header file for class documentation
  int iResult=0;
  int blockno=0;
  const AliHLTComponentBlockData* pDesc=NULL;

  Int_t TotalTracks=0;

  if ( GetFirstInputBlock( kAliHLTDataTypeSOR ) || GetFirstInputBlock( kAliHLTDataTypeEOR ) )
    return 0;
  
  for (pDesc=GetFirstInputBlock(AliHLTTPCDefinitions::fgkTrackSegmentsDataType); pDesc!=NULL; pDesc=GetNextInputBlock(), blockno++) {
    HLTDebug("event %Lu block %d: %s 0x%08x size %d", evtData.fEventID, blockno, DataType2Text(pDesc->fDataType).c_str(), pDesc->fSpecification, pDesc->fSize);
    
    if(pDesc->fDataType!=AliHLTTPCDefinitions::fgkTrackSegmentsDataType){continue;}
    
    iResult=PrintTrack(evtData,pDesc,TotalTracks);
  }
  
  for (pDesc=GetFirstInputBlock(AliHLTTPCDefinitions::fgkTracksDataType); pDesc!=NULL; pDesc=GetNextInputBlock(), blockno++) {
    HLTDebug("event %Lu block %d: %s 0x%08x size %d", evtData.fEventID, blockno, DataType2Text(pDesc->fDataType).c_str(), pDesc->fSpecification, pDesc->fSize);
    
    if(pDesc->fDataType!=AliHLTTPCDefinitions::fgkTracksDataType){continue;}
    
    iResult=PrintTrack(evtData,pDesc,TotalTracks);
  }
  HLTInfo("TrackDump found %d Tracks", TotalTracks);
  
  return iResult;
}

int AliHLTTPCTrackDumpComponent::PrintTrack(const AliHLTComponentEventData& evtData,const AliHLTComponentBlockData* bl,Int_t &nT){
  
  AliHLTUInt8_t slice = AliHLTTPCDefinitions::GetMinSliceNr( *bl );
  AliHLTUInt8_t patch = AliHLTTPCDefinitions::GetMinPatchNr( *bl );
  
  TString filename;
  int iResult=BuildFileName(evtData.fEventID, 0, bl->fDataType, 0, filename);
  ios::openmode filemode=(ios::openmode)0;
  if (fCurrentFileName.CompareTo(filename)==0) {
    filemode=ios::app;
  } else {
    fCurrentFileName=filename;
  }
  
  if (iResult>=0) {
    ofstream dump(filename.Data(), filemode);
    if (dump.good()) {
      
      HLTDebug ( "Input Data - TPC Tracks - Slice/Patch: %d/%d.", slice, patch );
      const AliHLTTPCTrackletData* trackData = (const AliHLTTPCTrackletData*) bl->fPtr;
      AliHLTUInt32_t nTracks = trackData->fTrackletCnt;
      AliHLTTPCTrackSegmentData *tracks = (AliHLTTPCTrackSegmentData*) trackData->fTracklets;
      
      for(AliHLTUInt32_t i=0;i<nTracks;i++){
	dump << "====================================================================" << endl;
	dump << "TrackNumber:   " << nT+i  << endl;
	dump << "Slice:         " << (unsigned int)slice << "     Partition:     "<<(unsigned int)patch <<endl;
	dump << "[X,Y,Z]:       [" << tracks->fX<<" , "<<tracks->fY<<" , "<<tracks->fZ <<"]"<< endl;
	dump << "[X,Y,Z](Last): [" << tracks->fLastX<<" , "<<tracks->fLastY<<" , "<<tracks->fLastZ <<"]"<< endl;
	dump << "pT:            " << tracks->fPt <<"\t\tpT Error:   " << tracks->fPterr <<endl;
	dump << "Psi:           " << tracks->fPsi <<"\t\tPsi Error:  " << tracks->fPsierr <<endl;
	dump << "Tgl:           " << tracks->fPt <<"\t\tTgl Error:  " << tracks->fPterr <<endl;
	dump << "Charge:        " << tracks->fCharge << "\t\tnClusters:  " << tracks->fNPoints << endl; 
	UChar_t *tmpP = (UChar_t*)tracks;
	tmpP += sizeof(AliHLTTPCTrackSegmentData)+tracks->fNPoints*sizeof(UInt_t);
	tracks = (AliHLTTPCTrackSegmentData*)tmpP;
      }
      nT+=nTracks;
    } 
    else {
      HLTError("can not open file %s for writing", fCurrentFileName.Data());
      iResult=-EBADF;
    }
    dump.close();
  } 
  return iResult;
}
