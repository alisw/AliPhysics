#ifdef __CINT__
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
 
#pragma link C++ enum   AliESDEvent::ESDListIndex;


#pragma link C++ class  AliESD+;
#pragma link C++ class  AliESDEvent+;
#pragma link C++ class  AliESDInputHandler+;
#pragma link C++ class  AliESDInputHandlerRP+;
#pragma link C++ class  AliESDRun+;
#pragma link C++ class  AliESDHeader+;
#pragma link C++ class  AliESDHLTDecision+;
#pragma link C++ class  AliESDZDC+;
#pragma link C++ class  AliESDCaloTrigger+;

#pragma read sourceClass="AliESDCaloTrigger" targetClass="AliESDCaloTrigger" source="Char_t fTriggerBits[48][64]" version="[2]" \
  target="fNEntries, fColumn, fRow, fTriggerBits" targetType="Int, Int_t*, Int_t*, Int_t*" code="{fTriggerBits = new Int_t[fNEntries]; for (Int_t i=0; i<fNEntries; ++i) fTriggerBits[i]=onfile.fTriggerBits[fColumn[i]][fRow[i]];}"

#pragma link C++ class  AliESDfriend+;                                                                                                           
#pragma read sourceClass="AliESDtrack" targetClass="AliESDtrack" source="UChar_t fTRDpidQuality"  version="[-49]" target="fTRDntracklets" targetType="UChar_t" code="{fTRDntracklets = onfile.fTRDpidQuality;}"
// see http://root.cern.ch/svn/root/trunk/io/doc/DataModelEvolution.txt
#pragma link C++ class  AliESDtrack+;
#pragma read sourceClass="AliESDfriendTrack" targetClass="AliESDfriendTrack" source="Int_t fITSindex" version="[-3]" \
        target="fnMaxITScluster, fITSindex" targetType="Int_t, Int_t*" code="{fnMaxITScluster = 12; fITSindex= new Int_t[fnMaxITScluster]; memcpy(fITSindex, &(onfile.fITSindex), fnMaxITScluster*sizeof(Int_t));}"
#pragma read sourceClass="AliESDfriendTrack" targetClass="AliESDfriendTrack" source="Int_t fTPCindex" version="[-3]" \
        target="fnMaxTPCcluster, fTPCindex" targetType="Int_t, Int_t*" code="{fnMaxTPCcluster = 160; fTPCindex= new Int_t[fnMaxTPCcluster]; memcpy(fTPCindex, &(onfile.fTPCindex), fnMaxTPCcluster*sizeof(Int_t));}"
#pragma read sourceClass="AliESDfriendTrack" targetClass="AliESDfriendTrack" source="Int_t fTRDindex" version="[-3]" \
        target="fnMaxTRDcluster, fTRDindex" targetType="Int_t, Int_t*" code="{fnMaxTRDcluster = 180; fTRDindex= new Int_t[fnMaxTRDcluster]; memcpy(fTRDindex, &(onfile.fTRDindex), fnMaxTRDcluster*sizeof(Int_t));}"

#pragma link C++ class  AliESDfriendTrack+;
#pragma link C++ class  AliESDMuonTrack+;
#pragma link C++ class  AliESDPmdTrack+;
#pragma link C++ class  AliESDTrdTrigger+;
#pragma link C++ class  AliESDTrdTrack+;
#pragma link C++ class  AliESDTrdTracklet+;
#pragma link C++ class  AliESDHLTtrack+;
#pragma link C++ class  AliESDv0+;
#pragma link C++ class  AliESDcascade+;
#pragma link C++ class  AliVertex+;
#pragma link C++ class  AliESDVertex+;
#pragma link C++ class  AliESDpid+;
#pragma link C++ class  AliESDkink+;
#pragma link C++ class  AliESDV0Params+;
#pragma link C++ class  AliESDCaloCluster+;
#pragma link C++ class  AliESDMuonCluster+;
#pragma link C++ class  AliESDMuonPad+;

#pragma link C++ class  AliKFParticleBase+;
#pragma link C++ class  AliKFParticle+;
#pragma link C++ class  AliKFVertex+;

#pragma link C++ class  AliKalmanTrack+;
#pragma link C++ class  AliVertexerTracks+;
#pragma link C++ class  AliStrLine+;
#pragma link C++ class  AliTrackPointArray+;
#pragma link C++ class  AliTrackPoint+;

#pragma link C++ class AliTrackPointArray+;
#pragma link C++ class AliTrackPoint+;

#pragma link C++ class  AliESDFMD+;
#pragma link C++ class  AliFMDMap+;
#pragma link C++ class  AliFMDFloatMap+;

#pragma link C++ class  AliESDVZERO+;
#pragma link C++ class  AliESDTZERO+;
#pragma link C++ class  AliESDACORDE+;
#ifdef MFT_UPGRADE
//#pragma link C++ class  AliESDMFT+;
#endif

#pragma link C++ class  AliESDMultITS+;
#pragma link C++ class  AliMultiplicity+;

#pragma link C++ class  AliSelector+;

#pragma link C++ class  AliRawDataErrorLog+;

#pragma link C++ class  AliMeanVertex+;
#pragma link C++ class  AliESDCaloCells+;

#pragma link C++ class  AliESDVZEROfriend+;
#pragma link C++ class  AliESDTZEROfriend+;

#pragma link C++ class  AliESDHandler+;
#pragma link C++ class  AliTrackerBase+;
#pragma link C++ class  AliTOFHeader+;

#pragma link C++ namespace AliESDUtils;

#pragma link C++ class  AliTriggerIR+;
#pragma link C++ class  AliTriggerScalersESD+;
#pragma link C++ class  AliTriggerScalersRecordESD+;
#pragma link C++ class AliTriggerCluster+;
#pragma link C++ class AliTriggerDescriptor+;
#pragma link C++ class AliTriggerInput+;
#pragma link C++ class AliTriggerInteraction+;
#pragma link C++ class AliTriggerPFProtection+;
#pragma link C++ class AliTriggerBCMask+;
#pragma link C++ class AliTriggerClass+;
#pragma link C++ class AliTriggerConfiguration+;
#pragma link C++ class AliExpression+;
#pragma link C++ class AliVariableExpression+;
#pragma link C++ class AliTPCdEdxInfo+;

#pragma link C++ function AliESDUtils::GetCorrV0(const AliESDEvent*,Float_t &);
#pragma link C++ function AliESDUtils::GetCorrSPD2(Float_t,Float_t);
#pragma link C++ function operator*(const AliFMDMap&,const AliFMDMap&);
#pragma link C++ function operator/(const AliFMDMap&,const AliFMDMap&);
#pragma link C++ function operator+(const AliFMDMap&,const AliFMDMap&);
#pragma link C++ function operator-(const AliFMDMap&,const AliFMDMap&);
  
#endif
