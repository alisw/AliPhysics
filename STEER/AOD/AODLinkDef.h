#ifdef __CINT__
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link off all methods;
 
#pragma link C++ enum   AliAODVertex::AODVtx_t;
#pragma link C++ enum   AliAODTrack::AODTrk_t;
#pragma link C++ enum   AliAODTrack::AODTrkPID_t;
//#pragma link C++ enum   AliAODCluster::AODClu_t;
//#pragma link C++ enum   AliAODCluster::AODCluPID_t;

#pragma link C++ class AliAODEvent+;
#pragma link C++ class AliAODHeader+;

#pragma read                                              \
    sourceClass="AliAODPid"                               \
    targetClass="AliAODPid"                               \
    source="UShort_t fTPCsignalN; Double_t fTPCmomentum; Float_t fTRDmomentum[6]"  \
    version="[-10]"                                       \
    target="fTPCsignalN, fTPCmomentum, fTRDmomentum"                                          \
    code="{newObj->SetTPCsignalN((UChar_t)onfile.fTPCsignalN); newObj->SetTPCmomentum(onfile.fTPCmomentum); for (Int_t i=0;i<6;++i) newObj->SetTRDmomentum(i,onfile.fTRDmomentum[i]);}" 

#pragma link C++ class AliAODPid+;


#pragma read sourceClass="AliAODTrack" targetClass="AliAODTrack" source="Double32_t fPID[10]"  version="[-22]" \
 target="fPID" targetType="Double32_t*" \
  code="{if (!fPID) fPID = new Double32_t[10];for(Int_t isp=10;isp--;) fPID[isp]=onfile.fPID[isp];}"

#pragma link C++ class AliAODTrack+;

#pragma link C++ class AliAODVertex+;
#pragma link C++ class AliAODCluster+;
#pragma link C++ class AliAODCaloCluster+;
#pragma link C++ class AliAODPmdCluster+;
#pragma link C++ class AliAODFmdCluster+;
#pragma link C++ class AliAODJet+;
#pragma link C++ class AliAODJetEventBackground+;
#pragma link C++ class AliAODPhoton+;
#pragma link C++ class AliAODRedCov<3>+;
#pragma link C++ class AliAODRedCov<4>+;
#pragma link C++ class AliAODRedCov<6>+;
#pragma link C++ class AliAODRecoDecay+;
#pragma link C++ class AliAODv0+;
#pragma link C++ class AliAODcascade+;
#pragma link C++ class AliAODHandler+;
#pragma link C++ class AliAODExtension+;
#pragma link C++ class AliAODBranchReplicator+;
#pragma link C++ class AliAODInputHandler+;
#pragma link C++ class AliAODTracklets+;
#pragma link C++ class AliAODTagCreator+;
#pragma link C++ class AliAODCaloCells+;
#pragma link C++ class AliAODCaloTrigger+;
#pragma link C++ class AliAODDiJet+;
#pragma link C++ class AliAODMCParticle+;
#pragma link C++ class AliAODMCHeader+;
#pragma link C++ class AliAODPWG4Particle+;
#pragma link C++ class AliAODPWG4ParticleCorrelation+;
#pragma link C++ class AliAODDimuon+;
#pragma link C++ class AliAODpidUtil+;
#pragma link C++ class AliAODTZERO+;
#pragma link C++ class AliAODVZERO+;
#pragma link C++ class AliAODZDC+;
#pragma link C++ class AliAODAD+;
#pragma link C++ class AliAODHMPIDrings+;
#pragma link C++ class AliAODTrdTrack+;
#pragma link C++ class AliAODTrdTracklet+;
#pragma link C++ class AliNanoAODTrackMapping+;
#pragma link C++ class AliNanoAODStorage+;
#pragma link C++ class AliNanoAODHeader+;

#pragma link C++ method AliAODTrack::SetPosition<double>(double const*, bool);

#endif
