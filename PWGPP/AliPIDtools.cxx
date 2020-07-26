#include "AliESDEvent.h"
#include "AliITSPIDResponse.h"
#include "AliTPCPIDResponse.h"
#include "AliTOFPIDResponse.h"
#include "AliPIDResponse.h"
#include "AliESDtrack.h"
#include "AliPIDtools.h"
#include "TLeaf.h"

std::map<Int_t, AliTPCPIDResponse *> AliPIDtools::pidTPC;     /// we should use better hash map
std::map<Int_t, AliPIDResponse *> AliPIDtools::pidAll;        /// we should use better hash map
AliESDtrack  AliPIDtools::dummyTrack;/// dummy value to save CPU - unfortunately PID object use AliVtrack - for the moment create global variable t avoid object constructions
TTree *       AliPIDtools::fFilteredTree = NULL;
TTree *       AliPIDtools::fFilteredTreeV0 = NULL;

AliPIDResponse* AliPIDtools::GetPID(Int_t hash ) {return pidAll[hash];}
AliTPCPIDResponse& AliPIDtools::GetTPCPID(Int_t hash ) {return pidAll[hash]->GetTPCResponse();}
AliITSPIDResponse& AliPIDtools::GetITSPID(Int_t hash ) {return pidAll[hash]->GetITSResponse();}
AliTOFPIDResponse& AliPIDtools::GetTOFPID(Int_t hash ) {return pidAll[hash]->GetTOFResponse();}

Int_t AliPIDtools::GetHash(Int_t run, Int_t passNumber, TString recoPass,Bool_t isMC){
  recoPass+=run;
  recoPass+=passNumber;
  recoPass+=isMC;
  return recoPass.Hash();
}

Double_t AliPIDtools::BetheBlochAleph(Int_t hash, Double_t bg){
  AliTPCPIDResponse *tpcPID=pidTPC[hash];
  if (tpcPID) return tpcPID->Bethe(bg);
  return 0;
}
Double_t AliPIDtools::BetheBlochAleph(Int_t hash, Double_t p,Int_t type){
  AliTPCPIDResponse *tpcPID=pidTPC[hash];
  Float_t bg = p/AliPID::ParticleMass(type);
  if (tpcPID) return tpcPID->Bethe(bg);
  return 0;
}

///  AliPIDtools::BetheBlochITS
/// \param hash   - hash value
/// \param p      - momentum (where?)
/// \param mass   - mass
/// \return
Double_t AliPIDtools::BetheBlochITS(Int_t hash, Double_t p, Double_t mass){
  if (pidAll[hash]== nullptr) return 0;
  AliITSPIDResponse &itsPID=GetITSPID(hash);
  return itsPID.Bethe(p, mass);
}

/// AliPIDtools::GetExpectedITSSignal(
/// \param hash   - hash value
/// \param p      - momentum (where?)
/// \return
Double_t AliPIDtools::GetExpectedITSSignal(Int_t hash, Double_t p, Int_t  particle){
  if (pidAll[hash]== nullptr) return 0;
  AliITSPIDResponse &itsPID=GetITSPID(hash);
  return itsPID.Bethe(p, (AliPID::EParticleType)particle);
}


/// GetExpectedTPCSignal
/// \param hash       - hash value of the PID version
/// \param p          - momenta
/// \param particle   - particle type
/// \return           - mean TPCdedx
Double_t AliPIDtools::GetExpectedTPCSignal(Int_t hash, Double_t p, Int_t  particle) {
  Double_t xyz[3] = {0., 0., 0.};
  Double_t pxyz[3] = {0, 0., 0.};
  Double_t cv[21] = {0.}; // dummy parameters for dummy tracks
  AliTPCPIDResponse *tpcPID=pidTPC[hash];
  if (tpcPID==0) return 0;
  pxyz[0]=p;
  dummyTrack.Set(xyz, pxyz, cv, 1);
  Double_t dEdx = tpcPID->GetExpectedSignal(&dummyTrack, (AliPID::EParticleType)particle, AliTPCPIDResponse::kdEdxDefault, kFALSE, kTRUE);
  return dEdx;
}
/// Load and reguster PID objects in hash maps
/// \param run
/// \param passNumber
/// \param recoPass
/// \param isMC
/// \return  - hash value of PID
Int_t AliPIDtools::LoadPID(Int_t run, Int_t passNumber, TString recoPass, Bool_t isMC){
  // Int_t run=246751, passNumber=1; TString recoPass("pass1"); Bool_t isMC=0;
  AliESDEvent ev;
  AliPIDResponse *pid = new AliPIDResponse(isMC);
  pid->SetUseTPCMultiplicityCorrection(kTRUE);
  pid->SetUseTPCEtaCorrection(kTRUE);
  pid->SetUseTPCPileupCorrection(kTRUE);
  pid->SetOADBPath("$ALICE_PHYSICS/OADB");
  pid->InitialiseEvent(&ev,passNumber, recoPass, run);
  AliTPCPIDResponse &tpcpid=pid->GetTPCResponse();
  // pid.InitFromOADB(246751,1,"pass1");
  Int_t  hash=GetHash(run,passNumber, recoPass,isMC);
  pidAll[hash]=pid;     /// we should clone them
  pidTPC[hash]=&tpcpid;  ///
  return hash;
}

Double_t AliPIDtools::GetExpectedTOFSigma(Int_t hash, Float_t mom, Int_t  type){
  Double_t dummyTime=0;
  if (pidAll[hash]== nullptr) return 0;
  AliTOFPIDResponse &tofPID=GetTOFPID(hash);
  return tofPID.GetExpectedSigma(mom,dummyTime,(AliPID::EParticleType)type);

}
Double_t AliPIDtools::GetExpectedTOFSignal(Int_t hash, const AliVTrack *track, Int_t type){
  if (pidAll[hash]== nullptr) return 0;
  AliTOFPIDResponse &tofPID=GetTOFPID(hash);
  return tofPID.GetExpectedSignal(track, (AliPID::EParticleType)type);
}

///  SetFiltered tree
/// \param filteredTree   - pointer to filtered tree
/// \return
Bool_t AliPIDtools::SetFilteredTree(TTree * filteredTree){
  if (filteredTree==NULL) return kFALSE;
  TBranch * branch = filteredTree->GetBranch("esdTrack.");
  if (branch==NULL) {
    ::Error("AliPIDtools::SetFilteredTreeV0","Invalid tree. Branch esdTrack does not exist");
    return kFALSE;
  }
  fFilteredTree=filteredTree;
  return kTRUE;
}

///  SetFiltered tree
Bool_t AliPIDtools::SetFilteredTreeV0(TTree * filteredTreeV0){
  if (filteredTreeV0==NULL) return kFALSE;
  TBranch * branch = filteredTreeV0->GetBranch("v0.");
  if (branch==NULL) {
    ::Error("AliPIDtools::SetFilteredTreeV0","Invalid tree. Branch v0 does not exist");
    return kFALSE;
  }
  fFilteredTreeV0=filteredTreeV0;
  return kTRUE;
}




/// GetExpected TPC signal for current track
/// \param hash                 - PID hash
/// \param particleType         - assumed particle type
/// \param corrMask             - corr mask
///                                 0x1 - eta correction
///                                 0x2 - multiplicity correction
///                                 0x4 - pile-up correction kPileUp
///                                 0x8  - return pileup correction
/// \param returnType            0
///                                 0 - expected signal
///                                 1 - corrected signal
/// \return                     - expected dEdx signal
Double_t AliPIDtools::GetExpectedTPCSignal(Int_t hash, Int_t particleType, Int_t corrMask, Int_t returnType){
  //
  AliTPCPIDResponse *tpcPID=pidTPC[hash];
  Double_t dEdx=0;
  AliESDtrack **pptrack=0;
  TVectorF   **pptpcVertexInfo=0;
  TVectorF   **ppitsClustersPerLayer=0;
  //Float_t primMult=0;
  Bool_t corrPileUp=corrMask&kPileUpCorr;
  if (fFilteredTree){  // data from filtered trees
    Int_t entry = fFilteredTree->GetReadEntry();
    static  TBranch * branch = NULL;
    static TBranch *branchVertex=0;
    static TBranch *branchITS=0;
    static TLeaf * leafPrim=0;
    static  Int_t treeNumber=-1;
    static TLeaf * leaftpcTrackBeforeClean=0;
    fFilteredTree->GetEntry(entry);   // load full tree - branch GetEntry is loading only for fit file in TChain  //TODO fix
    if (treeNumber!=fFilteredTree->GetTreeNumber()){
      branch=fFilteredTree->GetTree()->GetBranch("esdTrack.");
      if (fFilteredTree->GetFriend("E")) {
        branchVertex = fFilteredTreeV0->GetFriend("E")->GetBranch("tpcVertexInfoESD.");
        branchITS = fFilteredTreeV0->GetFriend("E")->GetBranch("itsClustersPerLayer.");
        leafPrim = fFilteredTreeV0->GetFriend("E")->GetLeaf("primMult");
        leaftpcTrackBeforeClean = fFilteredTree->GetFriend("E")->GetLeaf("tpcTrackBeforeClean");
      }
      treeNumber=fFilteredTree->GetTreeNumber();
    }
    pptrack = (AliESDtrack **)(branch->GetAddress());
    if (corrPileUp && branchVertex!=NULL) {
      pptpcVertexInfo = (branchVertex != NULL) ? (TVectorF **) (branchVertex->GetAddress()) : NULL;
      ppitsClustersPerLayer = (branchITS != NULL) ? (TVectorF **) (branchITS->GetAddress()) : NULL;
      SetPileUpProperties(**pptpcVertexInfo, **ppitsClustersPerLayer, leafPrim->GetValue(), tpcPID);
    }
    if (corrMask==-1) return entry;
    if (corrMask==-2) return (*pptrack)->Pt();
    if (corrMask==-3) return treeNumber;
    if (corrMask==-4 && pptpcVertexInfo) return (*(*pptpcVertexInfo))[0];
    //if (corrMask==-100) return tpcPID->GetP ileUpProperties(0);   - wait until new AliRoot distributed
    //if (corrMask==-101) return tpcPID->GetPileUpProperties(1);
    //if (corrMask==-102) return tpcPID->GetPileUpProperties(2);
    tpcPID->SetCurrentEventMultiplicity(leaftpcTrackBeforeClean->GetValue());
  }
  if (pptrack==0) return 0;
  if (corrMask==0x8) return tpcPID->GetPileupCorrectionValue(*pptrack);
  if (returnType==0) {
    dEdx = tpcPID->GetExpectedSignal(*pptrack, (AliPID::EParticleType) particleType, AliTPCPIDResponse::kdEdxDefault, corrMask & kEtaCorr, corrMask & kMultCorr, corrMask & kPileUpCorr);
    return dEdx;
  }
  if (returnType==1) {
    dEdx = tpcPID->GetCorrectedTrackdEdx(*pptrack, (AliPID::EParticleType) particleType, corrMask & kEtaCorr, corrMask & kMultCorr, corrMask & kPileUpCorr, AliTPCPIDResponse::kdEdxDefault);
    return dEdx;
  }
  return 0;
}

/// GetExpected TPC signal for current V0 track
/// \param hash                 - PID hash
/// \param particleType         - assumed particle type
/// \param corrMask             - corr mask (see TPCCorrFlag)
///                                 0x1 - eta correction
///                                 0x2 - multiplicity correction
///                                 0x4 - pile-up correction
///                                 0x8  - return pileup correction
/// \param index                - track index
/// \return                     - expected dEdx signal
Double_t AliPIDtools::GetExpectedTPCSignalV0(Int_t hash, Int_t particleType, Int_t corrMask, Int_t index, Int_t returnType){
  //
  AliTPCPIDResponse *tpcPID=pidTPC[hash];
  Double_t dEdx=0;
  AliESDtrack **pptrack=0;
  TVectorF   **pptpcVertexInfo=0;
  TVectorF   **ppitsClustersPerLayer=0;
  Bool_t corrPileUp=corrMask&kPileUpCorr;
  if (fFilteredTreeV0){  // data from filtered trees
    Int_t entry = fFilteredTreeV0->GetReadEntry();
    static TBranch * branch0, *branch1 = NULL;
    static TBranch *branchVertex=0;
    static TBranch *branchITS=0;
    static Int_t treeNumber=-1;
    static TLeaf * leafPrim=0;
    static TLeaf * leaftpcTrackBeforeClean=0;
    fFilteredTreeV0->GetEntry(entry);   // load full tree - branch GetEntry is loading only for fit file in TChain  //TODO fix
    if (treeNumber!=fFilteredTreeV0->GetTreeNumber()){
      branch0=fFilteredTreeV0->GetTree()->GetBranch("track0.");
      branch1=fFilteredTreeV0->GetTree()->GetBranch("track1.");
      if (fFilteredTree->GetFriend("E")) {
          branchVertex = fFilteredTreeV0->GetFriend("E")->GetBranch("tpcVertexInfoESD.");
          branchITS = fFilteredTreeV0->GetFriend("E")->GetBranch("itsClustersPerLayer.");
          leafPrim = fFilteredTreeV0->GetFriend("E")->GetLeaf("primMult");
          leaftpcTrackBeforeClean = fFilteredTreeV0->GetFriend("E")->GetLeaf("tpcTrackBeforeClean");
      }
      treeNumber=fFilteredTreeV0->GetTreeNumber();
    }
    pptrack = (index==0) ? (AliESDtrack **)(branch0->GetAddress()):(AliESDtrack **)(branch1->GetAddress());
    if (corrPileUp&& branchVertex!=NULL) {
      pptpcVertexInfo = (branchVertex != NULL) ? (TVectorF **) (branchVertex->GetAddress()) : NULL;
      ppitsClustersPerLayer = (branchITS != NULL) ? (TVectorF **) (branchITS->GetAddress()) : NULL;
      SetPileUpProperties(**pptpcVertexInfo, **ppitsClustersPerLayer, leafPrim->GetValue(), tpcPID);
    }
    if (corrMask==-1) return entry;
    if (corrMask==-2) return (*pptrack)->Pt();
    if (corrMask==-3) return treeNumber;
    if (corrMask==-4 && pptpcVertexInfo) return (*(*pptpcVertexInfo))[0];
    tpcPID->SetCurrentEventMultiplicity(leaftpcTrackBeforeClean->GetValue());
  }

  if (pptrack==0) return 0;
  if (corrMask==0x8) return tpcPID->GetPileupCorrectionValue(*pptrack);
  if (returnType==0) {
    dEdx = tpcPID->GetExpectedSignal(*pptrack, (AliPID::EParticleType) particleType, AliTPCPIDResponse::kdEdxDefault, corrMask & kEtaCorr, corrMask & kMultCorr, corrMask & kPileUpCorr);
    return dEdx;
  }
  if (returnType==1) {
    dEdx = tpcPID->GetCorrectedTrackdEdx(*pptrack, (AliPID::EParticleType) particleType, corrMask & kEtaCorr, corrMask & kMultCorr, corrMask & kPileUpCorr, AliTPCPIDResponse::kdEdxDefault);
    return dEdx;
  }
  return 0;
}


Bool_t AliPIDtools::SetPileUpProperties(const TVectorF & tpcVertexInfo, const TVectorF &itsClustersPerLayer, Int_t primMult, AliTPCPIDResponse *pidTPC){
  const Float_t  itsToTPC=2.38;        // conversion from SDD+SSD to TPC multiplicity
  const Float_t  multFraction=0.05;   //
  // ===| calculate derived variables |=========================================
  const Double_t shiftM = 0.5 * (tpcVertexInfo[1] + tpcVertexInfo[0]) - 25.;
  const Double_t multSSD = itsClustersPerLayer[4] + itsClustersPerLayer[5];
  const Double_t multSDD = itsClustersPerLayer[2] + itsClustersPerLayer[3];
  const Double_t multITSTPC = (multSSD + multSDD) / itsToTPC;
  const Double_t nPileUpSumCorr = (tpcVertexInfo[3] + tpcVertexInfo[4]) - multFraction * multITSTPC;
  const Double_t nPileUpPrim = nPileUpSumCorr / (1. - TMath::Abs(shiftM / 210.));
  // ===| set pileup event properties |=========================================
  pidTPC->SetEventPileupProperties(shiftM,nPileUpPrim,primMult);
  return kTRUE;
}

AliESDtrack* AliPIDtools::GetCurrentTrack() {
  AliESDtrack **pptrack = 0;
  if (fFilteredTree) {  // data from filtered trees
    Int_t entry = fFilteredTree->GetReadEntry();
    static Int_t lastEntry=-1;
    static TBranch *branch = NULL;
    static Int_t treeNumber = -1;
    if (lastEntry!=entry) {
      fFilteredTree->GetEntry(entry);   // load full tree - branch GetEntry is loading only for fit file in TChain
      if (treeNumber != fFilteredTree->GetTreeNumber()) {
        branch = fFilteredTree->GetTree()->GetBranch("esdTrack.");
        treeNumber = fFilteredTree->GetTreeNumber();
      }
      lastEntry=entry;
    }
    pptrack = (AliESDtrack **) (branch->GetAddress());
  }
  return *pptrack;
}

AliESDtrack* AliPIDtools::GetCurrentTrackV0(Int_t index) {
  AliESDtrack **pptrack = 0;
  if (fFilteredTreeV0) {  // data from filtered trees
    Int_t entry = fFilteredTreeV0->GetReadEntry();
    static TBranch *branch0, *branch1 = NULL;
    static Int_t treeNumber = -1;
    fFilteredTreeV0->GetEntry(entry);   // load full tree - branch GetEntry is loading only for fit file in TChain  //TODO fix
    if (treeNumber != fFilteredTreeV0->GetTreeNumber()) {
      branch0 = fFilteredTreeV0->GetTree()->GetBranch("track0.");
      branch1 = fFilteredTreeV0->GetTree()->GetBranch("track1.");
      treeNumber = fFilteredTreeV0->GetTreeNumber();
    }
    pptrack = (index == 0) ? (AliESDtrack **) (branch0->GetAddress()) : (AliESDtrack **) (branch1->GetAddress());
    return *pptrack;
  }
  return 0;
}

/// return TOF info for the V0 track
/// \param source         - 0 track, 1 track
/// \param infoType       - 0 cluster info, 1 - nsigma info
/// \return
TVectorD*     AliPIDtools::GetTOFInfo(Int_t infoType){
  if (fFilteredTree==0)  return 0;
  TVectorD ** tofInfo=0;
  Int_t entry = fFilteredTree->GetReadEntry();
  static TBranch *branchCl=NULL, *branchSigma = NULL;
  static Int_t treeNumber = -1;
  fFilteredTree->GetEntry(entry);   // load full tree - branch GetEntry is loading only for fit file in TChain  //TODO fix
  if (treeNumber != fFilteredTree->GetTreeNumber()) {
      branchCl = fFilteredTree->GetTree()->GetBranch("tofClInfo.");
      branchSigma = fFilteredTree->GetTree()->GetBranch("tofNsigma.");
      treeNumber = fFilteredTree->GetTreeNumber();
  }
  tofInfo =(infoType==0) ? (TVectorD**)branchCl->GetAddress(): (TVectorD**)branchSigma->GetAddress();
  if (tofInfo==NULL) return NULL;
  return *tofInfo;
}

/// return TOF info for the V0 track
/// \param source         - 0 track, 1 track
/// \param infoType       - 0 cluster info, 1 - nsigma info
/// \return
TVectorD*     AliPIDtools::GetTOFInfoV0(Int_t source, Int_t infoType){
  if (fFilteredTreeV0==0)  return 0;
  TVectorD ** tofInfo=0;
  Int_t entry = fFilteredTreeV0->GetReadEntry();
  static TBranch *branchCl0, *branchCl1 = NULL;
  static TBranch *branchSigma0, *branchSigma1 = NULL;
  static Int_t treeNumber = -1;
  fFilteredTreeV0->GetEntry(entry);   // load full tree - branch GetEntry is loading only for fit file in TChain  //TODO fix
  if (treeNumber != fFilteredTreeV0->GetTreeNumber()) {
      branchCl0 = fFilteredTreeV0->GetTree()->GetBranch("tofClInfo0.");
      branchCl1 = fFilteredTreeV0->GetTree()->GetBranch("tofClInfo1.");
      branchSigma0 = fFilteredTreeV0->GetTree()->GetBranch("tofNsigma0.");
      branchSigma1 = fFilteredTreeV0->GetTree()->GetBranch("tofNsigma1.");
      treeNumber = fFilteredTreeV0->GetTreeNumber();
  }
  if (source==0) tofInfo =(infoType==0) ? (TVectorD**)branchCl0->GetAddress(): (TVectorD**)branchSigma0->GetAddress();
  if (source==1) tofInfo =(infoType==0) ? (TVectorD**)branchCl1->GetAddress(): (TVectorD**)branchSigma1->GetAddress();
  if (tofInfo==NULL) return NULL;
  return *tofInfo;
}


/// Set TPCPIDResponse event information
/// \param pidHash      - pidhash  index of the response
/// \param corrPileUp   - switch -
/// \return
Bool_t       AliPIDtools::SetTPCEventInfo(Int_t pidHash,Int_t corrMaskTPC){
  if (fFilteredTree==NULL) return kFALSE;
  AliTPCPIDResponse *tpcPID=pidTPC[pidHash];
  if (tpcPID == NULL) return kFALSE;
  TVectorF   **pptpcVertexInfo=0;
  TVectorF   **ppitsClustersPerLayer=0;
  Int_t entry = fFilteredTree->GetReadEntry();
  static TBranch *branchVertex=0;
  static TBranch *branchITS=0;
  static Int_t treeNumber=-1;
  static TLeaf * leafPrim=0;
  //static TLeaf * leaftpcClusterMult=0;
  static TLeaf * leaftpcTrackBeforeClean=0;
    static TLeaf *leafGID=0;
  static TLeaf *leafGIDEv=0;
  fFilteredTree->GetEntry(entry);   // load full tree - branch GetEntry is loading only for fit file in TChain  //TODO fix
  Bool_t reset= (branchVertex)? (branchVertex->GetTree() != fFilteredTree->GetTree()):kFALSE;
  if (reset||treeNumber!=fFilteredTree->GetTreeNumber()) {
    if (fFilteredTree->GetFriend("E")) {
      branchVertex = fFilteredTree->GetFriend("E")->GetBranch("tpcVertexInfoESD.");
      branchITS = fFilteredTree->GetFriend("E")->GetBranch("itsClustersPerLayer.");
      leafPrim = fFilteredTree->GetFriend("E")->GetLeaf("primMult");
      //leaftpcClusterMult = fFilteredTree->GetFriend("E")->GetLeaf("tpcClusterMult");
      leaftpcTrackBeforeClean = fFilteredTree->GetFriend("E")->GetLeaf("tpcTrackBeforeClean");
      leafGIDEv=fFilteredTree->GetFriend("E")->GetLeaf("gid");
      leafGID=fFilteredTree->GetLeaf("gid");
    }else{
      return kFALSE;
    }
    treeNumber = fFilteredTree->GetTreeNumber();
  }
  if (leafGID->GetValueLong64()!=leafGIDEv->GetValueLong64()){
    Long64_t gid0=leafGID->GetValueLong64();
    Long64_t gidEv=leafGIDEv->GetValueLong64();
    ::Error("AliPIDtools::SetTPCEventInfo", "invalid gid number Ev=%d, Tree= %d, gid0 =%llu gidEv=%llu",entry,treeNumber, gid0,gidEv);
    return kFALSE;
  }
  if ((corrMaskTPC&kPileUpCorr)&& branchVertex!=NULL) {
    pptpcVertexInfo = (branchVertex != NULL) ? (TVectorF **) (branchVertex->GetAddress()) : NULL;
    ppitsClustersPerLayer = (branchITS != NULL) ? (TVectorF **) (branchITS->GetAddress()) : NULL;
    SetPileUpProperties(**pptpcVertexInfo, **ppitsClustersPerLayer, leafPrim->GetValue(), tpcPID);
  }
  tpcPID->SetCurrentEventMultiplicity(leaftpcTrackBeforeClean->GetValue());
  return kTRUE;
}


/// Set TPCPIDResponse event information
/// \param pidHash      - pidhash  index of the response
/// \param corrPileUp   - switch -
/// \return
Bool_t       AliPIDtools::SetTPCEventInfoV0(Int_t pidHash,Int_t corrMaskTPC){
  if (fFilteredTreeV0==NULL) return kFALSE;
  AliTPCPIDResponse *tpcPID=pidTPC[pidHash];
  if (tpcPID == NULL) return kFALSE;
  TVectorF   **pptpcVertexInfo=0;
  TVectorF   **ppitsClustersPerLayer=0;
  Int_t entry = fFilteredTreeV0->GetReadEntry();
  static TBranch *branchVertex=0;
  static TBranch *branchITS=0;
  static Int_t treeNumber=-1;
  static TLeaf * leafPrim=0;
  //static TLeaf * leaftpcClusterMult=0;
  static TLeaf * leaftpcTrackBeforeClean=0;
  static TLeaf *leafGID=0;
  static TLeaf *leafGIDEv=0;
  fFilteredTreeV0->GetEntry(entry);   // load full tree - branch GetEntry is loading only for fit file in TChain  //TODO fix
  Bool_t reset= (branchVertex)? (branchVertex->GetTree() != fFilteredTreeV0->GetTree()):kFALSE;
  Int_t entryEv=0;
  TTree *treeEv=0;
  if (reset||treeNumber!=fFilteredTreeV0->GetTreeNumber()) {
    leafGID=fFilteredTreeV0->GetLeaf("gid");
    treeEv=fFilteredTreeV0->GetFriend("E");
    if (treeEv) {
      entryEv=treeEv->GetReadEntry();
      branchVertex = treeEv->GetBranch("tpcVertexInfoESD.");
      branchITS = treeEv->GetBranch("itsClustersPerLayer.");
      leafPrim = treeEv->GetLeaf("primMult");
      //leaftpcClusterMult = treeEv->GetLeaf("tpcClusterMult");
      leaftpcTrackBeforeClean = treeEv->GetLeaf("tpcTrackBeforeClean");
      leafGIDEv=treeEv->GetLeaf("gid");
    }else{
      return kFALSE;
    }
    treeNumber = fFilteredTreeV0->GetTreeNumber();
  }
  if (leafGID->GetValueLong64()!=leafGIDEv->GetValueLong64()){
    Long64_t gid0=leafGID->GetValueLong64();
    Long64_t gidEv=leafGIDEv->GetValueLong64();
    ::Error("AliPIDtools::SetTPCEventInfoV0", "invalid gid number Ev=%d, Ev=%d,  Tree= %d, gid0 =%llu gidEv=%llu",entry,entryEv, treeNumber, gid0,gidEv);
    return kFALSE;
  }else{
      //::Info("AliPIDtools::SetTPCEventInfoV0", "invalid gid number Ev=%d, Tree= %d",entry,treeNumber);
  }
  if ((corrMaskTPC&kPileUpCorr)&& branchVertex!=NULL) {
    pptpcVertexInfo = (branchVertex != NULL) ? (TVectorF **) (branchVertex->GetAddress()) : NULL;
    ppitsClustersPerLayer = (branchITS != NULL) ? (TVectorF **) (branchITS->GetAddress()) : NULL;
    SetPileUpProperties(**pptpcVertexInfo, **ppitsClustersPerLayer, leafPrim->GetValue(), tpcPID);
  }
  Int_t tpcTrackBeforeClean=leaftpcTrackBeforeClean->GetValue();
  tpcPID->SetCurrentEventMultiplicity(tpcTrackBeforeClean);
  return kTRUE;
}




/// GetITSPID for given particle and valu type
/// TODO - interface other PID options
/// \param hash              - hash value of PID
/// \param particleType      - particle type
/// \param valueType         - return type  of value
///                            0 - delta
///                            1 - n sigma
///                            >1 - custom likelihood
/// \param resol
/// \return  value
Double_t AliPIDtools::GetITSPID(Int_t hash, Int_t particleType, Int_t valueType, Float_t resol){
  if (particleType>AliPID::kSPECIESC) return 0;
  AliITSPIDResponse &itsPID=pidAll[hash]->GetITSResponse();
  AliESDtrack *track=GetCurrentTrack();
  if (valueType==0) return itsPID.GetSignalDelta(track,(AliPID::EParticleType)particleType);
  if (valueType==1) return itsPID.GetNumberOfSigmas(track,(AliPID::EParticleType)particleType);
  Double_t prob[AliPID::kSPECIESC]={};
  if (resol>0){
    for (Int_t i=0;i<AliPID::kSPECIESC;i++){
      Float_t delta=itsPID.GetSignalDelta(track,(AliPID::EParticleType)i);
      prob[i]=TMath::Gaus(delta,resol);
    }
  }else{
    for (Int_t i=0;i<AliPID::kSPECIESC;i++){
      Float_t delta=itsPID.GetNumberOfSigmas(track,(AliPID::EParticleType)i);
      prob[i]=TMath::Gaus(delta,1);
    }
  }
  Double_t sumProb=0;
  for (Int_t i=0;i<AliPID::kSPECIESC;i++){sumProb+=prob[i];}
  if (sumProb==0) return 0;
  return prob[particleType]/sumProb;
}


/// Return PIDnsigma
/// \param hash           - hash value of PID correction
/// \param detCode        - detector code (0-ITS, 1-TPC, 2-TRD, 3-TOF)  AliPIDResponse::enum EDetector
/// \param particleType   - see enum
/// \param source         - track index
/// \param corrMask       - correction bitMask - AliPIDTools:: enum TPCCorrFlag
/// \return
Float_t AliPIDtools::NumberOfSigmas(Int_t hash, Int_t detCode, Int_t particleType, Int_t source, Int_t corrMask){
  if (pidAll[hash]==NULL) return 0;
  AliPIDResponse *pid = pidAll[hash];
  //
  Int_t maskBackup=0;                     // make backup of PID state
  if (pid->UseTPCEtaCorrection()) maskBackup+=kEtaCorr;
  if (pid->UseTPCMultiplicityCorrection()) maskBackup+=kMultCorr;
  if (pid->UseTPCPileupCorrection()) maskBackup+=kPileUpCorr;
  /// TODO - backup also strategy
  //
  if (corrMask<0) {
    corrMask=maskBackup;
  }else{
    pid->SetUseTPCEtaCorrection(corrMask&kEtaCorr);
    pid->SetUseTPCMultiplicityCorrection(corrMask&kMultCorr);
    pid->SetUseTPCPileupCorrection(corrMask&kPileUpCorr);
    if (corrMask&kPileUpCorr) pid->GetTPCResponse().SetPileupCorrectionStrategy(AliTPCPIDResponse::kPileupCorrectionInExpectedSignal);
  }
  AliESDtrack *track=NULL;
  Bool_t status=kTRUE;
  if (source<0){
    track=GetCurrentTrack();
    status=SetTPCEventInfo(hash,corrMask);
  }
  if (source>=0){
    track=GetCurrentTrackV0(source%2);
    status=SetTPCEventInfoV0(hash,corrMask);
  }
  if (status==kFALSE || track==NULL) return 0;
  if (detCode==3){ /// for TOF use tree values
    if (source<0) return AliPIDtools::GetTOFInfoAt(1,particleType);
    if (source>=0) return AliPIDtools::GetTOFInfoV0At(source,1,particleType);
  }
  Double_t value=pidAll[hash]->NumberOfSigmas((AliPIDResponse::EDetector) detCode, track, (AliPID::EParticleType)particleType);
  // restore flags
  pid->SetUseTPCEtaCorrection(kEtaCorr&maskBackup);
  pid->SetUseTPCMultiplicityCorrection(maskBackup&kMultCorr);
  pid->SetUseTPCPileupCorrection(maskBackup&kPileUpCorr);
  return value;
}

/// Return GetSignalDelta
/// \param hash           - hash value of PID correction
/// \param detCode        - detector code (0-ITS, 1-TPC, 2-TRD, 3-TOF)  AliPIDResponse::enum EDetector
/// \param particleType   - see enum
/// \param source         - track index
/// \param corrMask       - correction bitMask - AliPIDTools:: enum TPCCorrFlag
/// \return
Float_t AliPIDtools::GetSignalDelta(Int_t hash, Int_t detCode, Int_t particleType, Int_t source, Int_t corrMask){
  if (pidAll[hash]==NULL) return 0;
  AliPIDResponse *pid = pidAll[hash];
  //
  Int_t maskBackup=0;                     // make backup of PID state
  if (pid->UseTPCEtaCorrection()) maskBackup+=kEtaCorr;
  if (pid->UseTPCMultiplicityCorrection()) maskBackup+=kMultCorr;
  if (pid->UseTPCPileupCorrection()) maskBackup+=kPileUpCorr;
  //
  if (corrMask<0) {
    corrMask=maskBackup;
  }else{
    pid->SetUseTPCEtaCorrection(corrMask&kEtaCorr);
    pid->SetUseTPCMultiplicityCorrection(corrMask&kMultCorr);
    pid->SetUseTPCPileupCorrection(corrMask&kPileUpCorr);
  }
  AliESDtrack *track=NULL;
  if (source<0){
    track=GetCurrentTrack();
    SetTPCEventInfo(hash,corrMask);
  }
  if (source>=0){
    track=GetCurrentTrackV0(source%2);
    SetTPCEventInfoV0(hash,corrMask);
  }
  Double_t value=pidAll[hash]->GetSignalDelta((AliPIDResponse::EDetector) detCode, track, (AliPID::EParticleType)particleType);
  // restore flags
  pid->SetUseTPCEtaCorrection(kEtaCorr&maskBackup);
  pid->SetUseTPCMultiplicityCorrection(maskBackup&kMultCorr);
  pid->SetUseTPCPileupCorrection(maskBackup&kPileUpCorr);
  return value;
}


/// Return Compute probability
/// TOF information ofr ALIPIDResponse is not avalaible, tocNSigma used instead
/// \param hash           - hash value of PID correction
/// \param detCode        - detector code (0-ITS, 1-TPC, 2-TRD, 3-TOF)  AliPIDResponse::enum EDetector
/// \param particleType   - see enum
/// \param source         - track index
/// \param corrMask       - correction bitMask - AliPIDTools:: enum TPCCorrFlag
/// \param norm           - include normalization to all species
/// \param fakeProb       -  user defined fake probability (normaly scales with mult*(1+1/pt))  - detector dependent
/// \return
Float_t AliPIDtools::ComputePIDProbability(Int_t hash, Int_t detCode, Int_t particleType, Int_t source, Int_t corrMask,Int_t norm,Float_t fakeProb,Float_t *pidVector){
  if (pidAll[hash]==NULL) return 0;
  const Double_t kMaxSigma=4;
  AliPIDResponse *pid = pidAll[hash];
  //
  Int_t maskBackup=0;                     // make backup of PID state
  if (pid->UseTPCEtaCorrection()) maskBackup+=kEtaCorr;
  if (pid->UseTPCMultiplicityCorrection()) maskBackup+=kMultCorr;
  if (pid->UseTPCPileupCorrection()) maskBackup+=kPileUpCorr;
  //
  if (corrMask<0) {
    corrMask=maskBackup;
  }else{
    pid->SetUseTPCEtaCorrection(corrMask&kEtaCorr);
    pid->SetUseTPCMultiplicityCorrection(corrMask&kMultCorr);
    pid->SetUseTPCPileupCorrection(corrMask&kPileUpCorr);
    if (corrMask&kPileUpCorr) pid->GetTPCResponse().SetPileupCorrectionStrategy(AliTPCPIDResponse::kPileupCorrectionInExpectedSignal);
  }
  AliESDtrack *track=NULL;
  if (source<0){
    track=GetCurrentTrack();
    SetTPCEventInfo(hash,corrMask);
  }
  if (source>=0){
    track=GetCurrentTrackV0(source%2);
    SetTPCEventInfoV0(hash,corrMask);
  }
  //Double_t value=pidAll[hash]->GetSignalDelta((AliPIDResponse::EDetector) detCode, track, (AliPID::EParticleType)particleType);
  Double_t prob[AliPID::kSPECIESCN]={0};
  Bool_t status  =kTRUE;
  if (detCode!=3) status = pidAll[hash]->ComputePIDProbability( (AliPIDResponse::EDetector) detCode, track, AliPID::kSPECIESC, prob);
  else{ //special treatment for TOF
    TVectorD *tofSigma=(source==-1) ? GetTOFInfo(1):GetTOFInfoV0(source,1);
    //TVectorD *tofInfo=(source==-1) ? GetTOFInfo(0):GetTOFInfoV0(source,0);
    status=kFALSE;                    // time assigned to TOF cluster
    if (tofSigma) for (Int_t i=0; i<tofSigma->GetNrows();i++){
      Float_t nsigma=(*tofSigma)[i];
      if (TMath::Abs((*tofSigma)[0])<kMaxSigma)   {
        prob[i]=TMath::Exp(-0.5*nsigma*nsigma);
        status=kTRUE;      // assing status if measurement
      }
    }
  }
  Double_t value = (status==kTRUE) ? prob[particleType%AliPID::kSPECIESC]:0;
  if (norm>0){
    Double_t sumP=0;
    for (Int_t i=0; i<AliPID::kSPECIESC; i++) sumP+=prob[i];
    sumP+=fakeProb;
    value/=sumP;
  }
  //
  if (pidVector){
    for (Int_t i=0; i<AliPID::kSPECIESC; i++) pidVector[i]=prob[i];
    pidVector[AliPID::kSPECIESC]=status;
  }
  // restore flags
  pid->SetUseTPCEtaCorrection(kEtaCorr&maskBackup);
  pid->SetUseTPCMultiplicityCorrection(maskBackup&kMultCorr);
  pid->SetUseTPCPileupCorrection(maskBackup&kPileUpCorr);
  return value;
}

///
/// \param hash               - index of PID
/// \param detMask            - det mask as in AliPID
/// \param particleType       - particle type
/// \param source             -
/// \param corrMask           - corr mask for TPC only
/// \param norm               - flag - normalization
/// \param fakeProb           - fake probability  ( used for  normalization)
/// \return
Float_t AliPIDtools::ComputePIDProbabilityCombined(Int_t hash, Int_t detMask, Int_t particleType, Int_t source, Int_t corrMask,Int_t norm, Float_t fakeProb){
  Float_t pidVector[AliPID::kSPECIESC+1]={};
  for (Int_t i=0; i<AliPID::kSPECIESC; i++) pidVector[i]=1;
  Int_t nDetectors=0;
  for (Int_t iDet=0; iDet<AliPIDResponse::kNdetectors; ++iDet) {
    if ((detMask & (1 << iDet))) {
      Float_t pidVectorDet[AliPID::kSPECIESC + 1] = {1};
      ComputePIDProbability(hash,iDet,particleType%AliPID::kSPECIESC,source,corrMask,0,fakeProb, pidVectorDet);
      Float_t status =  pidVectorDet[AliPID::kSPECIESC];
      if (status>0){
        nDetectors++;
        for (Int_t i=0; i<AliPID::kSPECIESC; i++) pidVector[i]*=pidVectorDet[i]+fakeProb/AliPID::kSPECIESC;
      }
    }
  }
  if (nDetectors==0) return 1; // return default value 1 standard species

  if (norm&0x2){
    for (Int_t i=0; i<AliPID::kSPECIESC; i++) {
      pidVector[i]=TMath::Power(pidVector[i],1./nDetectors);
    }
  }
  if (norm&0x1){
    Float_t sum=0;
    for (Int_t i=0; i<AliPID::kSPECIESC; i++) { sum+=pidVector[i];}
    sum+=fakeProb*AliPID::kSPECIESC;
    for (Int_t i=0; i<AliPID::kSPECIESC; i++) { pidVector[i]/=sum;}
  }
  if (particleType==AliPID::kSPECIESC) return pidVector[1]+pidVector[2]; //return muon+pion
  return pidVector[particleType];
}


///
/// \param hash               - index of PID
/// \param detMask            - det mask as in AliPID
/// \param particleMask       - particle bitmask used to calculate sum .eg to join electron,muon and electron probabitity for ITS
/// \param source             -
/// \param corrMask           - corr mask for TPC only
/// \param norm               - flag - normalization
/// \param fakeProb           - fake probability  ( used for  normalization)
/// \return
Float_t AliPIDtools::ComputePIDProbabilityCombinedMask(Int_t hash, Int_t detMask, Int_t particleMask, Int_t source, Int_t corrMask,Int_t norm, Float_t fakeProb){
  Float_t pidVector[AliPID::kSPECIESC+1]={};
  for (Int_t i=0; i<AliPID::kSPECIESC; i++) pidVector[i]=1;
  Int_t nDetectors=0;
  for (Int_t iDet=0; iDet<AliPIDResponse::kNdetectors; ++iDet) {
    if ((detMask & (1 << iDet))) {
      Float_t pidVectorDet[AliPID::kSPECIESC + 1] = {1};
      ComputePIDProbability(hash,iDet,0,source,corrMask,0,fakeProb, pidVectorDet);
      Float_t status =  pidVectorDet[AliPID::kSPECIESC];
      if (status>0){
        nDetectors++;
        for (Int_t i=0; i<AliPID::kSPECIESC; i++) pidVector[i]*=pidVectorDet[i]+fakeProb/AliPID::kSPECIESC;
      }
    }
  }
  if (nDetectors==0) return 0.2; // return default value 1/5 standard species

  if (norm&0x2){
    for (Int_t i=0; i<AliPID::kSPECIESC; i++) {
      pidVector[i]=TMath::Power(pidVector[i],1./nDetectors);
    }
  }
  if (norm&0x1){
    Float_t sum=0;
    for (Int_t i=0; i<AliPID::kSPECIESC; i++) { sum+=pidVector[i];}
    sum+=fakeProb*AliPID::kSPECIESC;
    for (Int_t i=0; i<AliPID::kSPECIESC; i++) { pidVector[i]/=sum;}
  }
  Float_t probSum=0;
  for (Int_t i=0; i<AliPID::kSPECIESC; i++) {
     if ((particleMask & (1 << i))) {
       probSum+=pidVector[i];
     }
  }
  return probSum;
}



/// Unit test of invariants - check internal consistency of wrappers
void AliPIDtools::UnitTest() {
  Bool_t status=0;
  const Float_t kEpsilon=0.00001;
  Int_t entries=0;
  // Test TOF info interface
  entries=fFilteredTree->Draw("AliPIDtools::GetTOFInfoAt(1,2)-tofNsigma.fElements[2]","1","goff",100);
  status=TMath::RMS(entries, fFilteredTree->GetV1())<kEpsilon;
  ::Info("UnitTest","AliPIDtools::GetTOFInfoAt(1,2)-tofNsigma.fElements[2]\tStatus=%d",status);
  //   TOF Compute PID
  entries=fFilteredTree->Draw("AliPIDtools::ComputePIDProbability(pidHash,3,5,-1,3,0,0.01)-exp(-0.5*AliPIDtools::NumberOfSigmas(pidHash,3,5,-1,3)**2)","ITSRefit","goff",1000);
  status=TMath::RMS(entries, fFilteredTree->GetV1())<kEpsilon;
  ::Info("UnitTest","AliPIDtools::ComputePIDProbability(pidHash,3,5,-1,3,0,0.01)-exp(-0.5*AliPIDtools::NumberOfSigmas(pidHash,3,5,-1,3)**2)\tStatus=%d",status);
  // Test combined PID consistency
  //    ITS check
  entries=fFilteredTree->Draw("AliPIDtools::ComputePIDProbabilityCombined(pidHash,1,5,-1,3,0,0.0)-AliPIDtools::ComputePIDProbability(pidHash,0,5,-1,3,0,0.0)","ITSRefit","goff",1000);
  status=TMath::RMS(entries, fFilteredTree->GetV1())<kEpsilon;
  ::Info("UnitTest","AliPIDtools::ComputePIDProbabilityCombined(pidHash,1,5,-1,3,0,0.0)-AliPIDtools::ComputePIDProbability(pidHash,0,5,-1,3,0,0.0)\tStatus=%d",status);
  //   TOF combined check
  entries=fFilteredTree->Draw("AliPIDtools::ComputePIDProbabilityCombined(pidHash,8,5,-1,3,0,0.0)-exp(-0.5*AliPIDtools::NumberOfSigmas(pidHash,3,5,-1,3)**2)","ITSRefit","goff",1000);
  status=TMath::RMS(entries, fFilteredTree->GetV1())<kEpsilon;
  ::Info("UnitTest","AliPIDtools::ComputePIDProbabilityCombined(pidHash,8,5,-1,3,0,0.0)-exp(-0.5*AliPIDtools::NumberOfSigmas(pidHash,3,5,-1,3)**2)\tStatus=%d",status);
  //
  entries=fFilteredTree->Draw("AliPIDtools::ComputePIDProbabilityCombined(pidHash,10,2,-1,3+0,0,0.0)-AliPIDtools::ComputePIDProbability(pidHash,1,2,-1,3+0,0,0.0)*AliPIDtools::ComputePIDProbability(pidHash,3,2,-1,3+0,0,0.0)",
          "TOFOn&&abs(nSigma1_2)<5&&abs(nSigma3_2)<5","goff",1000);
  status=TMath::RMS(entries, fFilteredTree->GetV1())<kEpsilon;
  ::Info("UnitTest","AliPIDtools::ComputePIDProbabilityCombined(pidHash,10,2,-1,3+0,0,0.0)-AliPIDtools::ComputePIDProbability(pidHash,1,2,-1,3+0,0,0.0)*AliPIDtools::ComputePIDProbability(pidHash,3,2,-1,3+0,0,0.0)\tStatus=%d",status);

}

///
/// \param pidHash          - pid hash used to evaluate (more than one responce can be used)
/// \param fakeRate         - fake rate for particles - default is 10 % "0.1" - can be any expression used in trees
/// \param suffix           - suffix to add
/// \return
Bool_t AliPIDtools::RegisterPIDAliases(Int_t pidHash, TString fakeRate, Int_t suffix){
  if (pidAll[pidHash]==NULL){
    ::Error("AliPIDtools::RegisterPIDAliases","Invalid PID hash %d",pidHash);
    return kFALSE;
  }
  if (fFilteredTree==NULL){
    ::Error("AliPIDtools::RegisterPIDAliases","Non initialized tree");
    return kFALSE;
  }
  TString sSufix="";
  if (suffix>=0) sSufix=TString::Format("_%d",suffix);
  //
  for (Int_t jPID=0; jPID<AliPID::kSPECIESC+1; jPID++) {
    Int_t iPID = (jPID < AliPID::kSPECIESC) ? jPID : AliPID::kPion;
    fFilteredTree->SetAlias(Form("ldEdxITS%d%s", iPID, sSufix.Data()),
                            Form("log(esdTrack.fITSsignal/AliPIDtools::GetExpectedITSSignal(pidHash,esdTrack.fIp.P(),%d+0))", iPID));
    fFilteredTree->SetAlias(Form("ldEdxTRD%d%s", iPID, sSufix.Data()),
                            Form("log(50*esdTrack.fTRDsignal/AliPIDtools::GetExpectedITSSignal(pidHash,esdTrack.fOp.P(),%d+0))", iPID));
    fFilteredTree->SetAlias(Form("ldEdxTPC%d%s", iPID, sSufix.Data()),
                            Form("log(esdTrack.fTPCsignal/AliPIDtools::GetExpectedTPCSignal(pidHash,esdTrack.fIp.P(),%d+0))", iPID));
    fFilteredTreeV0->SetAlias(Form("ldEdx0ITS%d%s", iPID, sSufix.Data()),
                              Form("log(track0.fITSsignal/AliPIDtools::GetExpectedITSSignal(pidHash,esdTrack.fIp.P(),%d+0))", iPID));
    fFilteredTreeV0->SetAlias(Form("ldEdx0TPC%d%s", iPID, sSufix.Data()),
                              Form("log(track0.fTPCsignal/AliPIDtools::GetExpectedTPCSignal(pidHash,esdTrack.fIp.P(),%d+0))", iPID));
    fFilteredTreeV0->SetAlias(Form("ldEdx0TRD%d%s", iPID, sSufix.Data()),
                              Form("log(50*track0.fTRDsignal/AliPIDtools::GetExpectedTPCSignal(pidHash,esdTrack.fOp.P(),%d+0))", iPID));
    fFilteredTreeV0->SetAlias(Form("ldEdx1ITS%d%s", iPID, sSufix.Data()),
                              Form("log(track1.fITSsignal/AliPIDtools::GetExpectedITSSignal(pidHash,esdTrack.fIp.P(),%d+0))", iPID));
    fFilteredTreeV0->SetAlias(Form("ldEdx1TPC%d%s", iPID, sSufix.Data()),
                              Form("log(track1.fTPCsignal/AliPIDtools::GetExpectedTPCSignal(pidHash,esdTrack.fIp.P(),%d+0))", iPID));
    fFilteredTreeV0->SetAlias(Form("ldEdx1TRD%d%s", iPID, sSufix.Data()),
                              Form("log(track1.fTRDsignal/AliPIDtools::GetExpectedTPCSignal(pidHash,esdTrack.fOp.P(),%d+0))", iPID));

    for (Int_t iR = 0; iR <= 3; iR++) {
      fFilteredTree->SetAlias(Form("ldEdxMax%d_%d%s", iR, iPID, sSufix.Data()),
                              Form("log(fTPCdEdxInfo.GetSignalMax(%d)/AliPIDtools::GetExpectedTPCSignal(pidHash,esdTrack.fIp.P(),%d+0))", iR, iPID));
      fFilteredTree->SetAlias(Form("ldEdxTot%d_%d%s", iR, iPID, sSufix.Data()),
                              Form("log(fTPCdEdxInfo.GetSignalTot(%d)/AliPIDtools::GetExpectedTPCSignal(pidHash,esdTrack.fIp.P(),%d+0))", iR, iPID));
      fFilteredTree->SetAlias(Form("ldEdxMaxTot%d_%d%s", iR, iPID, sSufix.Data()),
                              Form("log(fTPCdEdxInfo.GetSignalMax(%d)/fTPCdEdxInfo.GetSignalTot(%d+0))", iR, iR));
      fFilteredTree->SetAlias(Form("ldEdxMax%d%d_%d%s", iR, (iR + 1) % 3, iPID, sSufix.Data()),
                              Form("log(fTPCdEdxInfo.GetSignalMax(%d)/fTPCdEdxInfo.GetSignalMax(%d+0))", iR, (iR + 1) % 3));
      fFilteredTree->SetAlias(Form("ldEdxTot%d%d_%d%s", iR, (iR + 1) % 3, iPID, sSufix.Data()),
                              Form("log(fTPCdEdxInfo.GetSignalTot(%d)/fTPCdEdxInfo.GetSignalTot(%d+0))", iR, (iR + 1) % 3));
      //
      //
      fFilteredTreeV0->SetAlias(Form("ldEdx0Max%d+%d%s", iR, iPID, sSufix.Data()),
                                Form("log(track0.fTPCdEdxInfo.GetSignalMax(%d)/AliPIDtools::GetExpectedTPCSignal(pidHash,track0.fIp.P(),%d+0))", iR, iPID));
      fFilteredTreeV0->SetAlias(Form("ldEdx0Tot%d_%d%s", iR, iPID, sSufix.Data()),
                                Form("log(track0.fTPCdEdxInfo.GetSignalTot(%d)/AliPIDtools::GetExpectedTPCSignal(pidHash,track0.fIp.P(),%d+0))", iR, iPID));
      fFilteredTreeV0->SetAlias(Form("ldEdx0MaxTot%d_%d%s", iR, iPID, sSufix.Data()),
                                Form("log(track0.fTPCdEdxInfo.GetSignalMax(%d)/track0.fTPCdEdxInfo.GetSignalTot(%d+0))", iR, iR));
      fFilteredTreeV0->SetAlias(Form("ldEdx1Max%d_%d%s", iR, iPID, sSufix.Data()),
                                Form("log(track1.fTPCdEdxInfo.GetSignalMax(%d)/AliPIDtools::GetExpectedTPCSignal(pidHash,track1.fIp.P(),%d+0))", iR, iPID));
      fFilteredTreeV0->SetAlias(Form("ldEdx1Tot%d_%d%s", iR, iPID, sSufix.Data()),
                                Form("log(track1.fTPCdEdxInfo.GetSignalTot(%d)/AliPIDtools::GetExpectedTPCSignal(pidHash,track1.fIp.P(),%d+0))", iR, iPID));
      fFilteredTreeV0->SetAlias(Form("ldEdx1MaxTot%d_%d%s", iR,iPID, sSufix.Data()),
                                Form("log(track1.fTPCdEdxInfo.GetSignalMax(%d)/track1.fTPCdEdxInfo.GetSignalTot(%d+0))", iR, iR));
      //
      fFilteredTreeV0->SetAlias(Form("ldEdx0Max%d%d%s", iR, (iR + 1) % 3, sSufix.Data()),
                                Form("log(track0.fTPCdEdxInfo.GetSignalMax(%d)/track0.fTPCdEdxInfo.GetSignalMax(%d+0))", iR, (iR + 1) % 3));
      fFilteredTreeV0->SetAlias(Form("ldEdx0Tot%d%d%s", iR, (iR + 1) % 3, sSufix.Data()),
                                Form("log(track0.fTPCdEdxInfo.GetSignalTot(%d)/track0.fTPCdEdxInfo.GetSignalTot(%d+0))", iR, (iR + 1) % 3));
      fFilteredTreeV0->SetAlias(Form("ldEdx1Max%d%d%s", iR, (iR + 1) % 3, sSufix.Data()),
                                Form("log(track1.fTPCdEdxInfo.GetSignalMax(%d)/track1.fTPCdEdxInfo.GetSignalMax(%d+0))", iR, (iR + 1) % 3));
      fFilteredTreeV0->SetAlias(Form("ldEdx1Tot%d%d%s", iR, (iR + 1) % 3, sSufix.Data()),
                                Form("log(track1.fTPCdEdxInfo.GetSignalTot(%d)/track1.fTPCdEdxInfo.GetSignalTot(%d+0))", iR, (iR + 1) % 3));
    }
  }

  for (Int_t iPID=0; iPID<AliPID::kSPECIESC; iPID++){
    for (Int_t iDet=0; iDet<5; iDet++){
      Float_t mass=AliPID::ParticleMass(iPID);
      Float_t charge=AliPID::ParticleCharge(iPID);
      fFilteredTree->SetAlias(Form("nSigma%d_%d%s",iDet,iPID,sSufix.Data()),Form("AliPIDtools::NumberOfSigmas(%d,%d,%d,-1,3+0)",pidHash,iDet,iPID));
      fFilteredTree->SetAlias(Form("prob%d_%d%s",iDet,iPID,sSufix.Data()),Form("AliPIDtools::ComputePIDProbability(%d,%d,%d,-1,3+0,0)",pidHash,iDet,iPID));
      fFilteredTree->SetAlias(Form("probN%d_%d%s",iDet,iPID,sSufix.Data()),Form("AliPIDtools::ComputePIDProbability(%d,%d,%d,-1,3+0,1,%s)",pidHash,iDet,iPID,fakeRate.Data()));
      fFilteredTree->SetAlias(Form("dEdxRatioITS75_%d",iPID), Form("log(AliPIDtools::BetheBlochITS(%d,AliPIDResponse::interpolateP(esdTrack.P(),esdTrack.fIp.P(), %f,0.75,%f),%f)/"
                                   "AliPIDtools::BetheBlochITS(%d,esdTrack.P(),%f))",pidHash,mass,charge,mass,pidHash,mass));
      fFilteredTree->SetAlias(Form("dEdxRatioITS50_%d",iPID), Form("log(AliPIDtools::BetheBlochITS(%d,AliPIDResponse::interpolateP(esdTrack.P(),esdTrack.fIp.P(), %f,0.5,%f),%f)/"
                                   "AliPIDtools::BetheBlochITS(%d,esdTrack.P(),%f))",pidHash,mass,charge,mass,pidHash,mass));
      //
      fFilteredTreeV0->SetAlias(Form("nSigma0_%d_%d%s",iDet,iPID,sSufix.Data()),Form("AliPIDtools::NumberOfSigmas(%d,%d,%d,0,3+0)",pidHash,iDet,iPID));
      fFilteredTreeV0->SetAlias(Form("prob0_%d_%d%s",iDet,iPID,sSufix.Data()),Form("AliPIDtools::ComputePIDProbability(%d,%d,%d,0,3+0,0)",pidHash,iDet,iPID));
      fFilteredTreeV0->SetAlias(Form("probN0_%d_%d%s",iDet,iPID,sSufix.Data()),Form("AliPIDtools::ComputePIDProbability(%d,%d,%d,0,3+0,1,%s)",pidHash,iDet,iPID,fakeRate.Data()));
      //
      fFilteredTreeV0->SetAlias(Form("nSigma1_%d_%d%s",iDet,iPID,sSufix.Data()),Form("AliPIDtools::NumberOfSigmas(%d,%d,%d,1,3+0)",pidHash,iDet,iPID));
      fFilteredTreeV0->SetAlias(Form("prob1_%d_%d%s",iDet,iPID,sSufix.Data()),Form("AliPIDtools::ComputePIDProbability(%d,%d,%d,1,3+0,0)",pidHash,iDet,iPID));
      fFilteredTreeV0->SetAlias(Form("probN1_%d_%d%s",iDet,iPID,sSufix.Data()),Form("AliPIDtools::ComputePIDProbability(%d,%d,%d,1,3+0,1,%s)",pidHash,iDet,iPID,fakeRate.Data()));
    }
  }
  /// combined PID prob
  for (Int_t iPID=0; iPID<8; iPID++){
    fFilteredTree->SetAlias(Form("probCN01_%d%s",iPID,sSufix.Data()),Form("AliPIDtools::ComputePIDProbabilityCombined(%d,3,%d,-1,3+0,0,%s)",pidHash,iPID,fakeRate.Data()));
    fFilteredTree->SetAlias(Form("probC01_%d%s", iPID,sSufix.Data()),Form("AliPIDtools::ComputePIDProbabilityCombined(%d,3,%d,-1,3+0,1,%s)",pidHash,iPID,fakeRate.Data()));
    fFilteredTree->SetAlias(Form("probCN03_%d%s",iPID,sSufix.Data()),Form("AliPIDtools::ComputePIDProbabilityCombined(%d,9,%d,-1,3+0,0,%s)",pidHash,iPID,fakeRate.Data()));
    fFilteredTree->SetAlias(Form("probC03_%d%s", iPID,sSufix.Data()),Form("AliPIDtools::ComputePIDProbabilityCombined(%d,9,%d,-1,3+0,1,%s)",pidHash,iPID,fakeRate.Data()));
    fFilteredTree->SetAlias(Form("probCN13_%d%s",iPID,sSufix.Data()),Form("AliPIDtools::ComputePIDProbabilityCombined(%d,10,%d,-1,3+0,0,%s)",pidHash,iPID,fakeRate.Data()));
    fFilteredTree->SetAlias(Form("probC13_%d%s", iPID,sSufix.Data()),Form("AliPIDtools::ComputePIDProbabilityCombined(%d,10,%d,-1,3+0,1,%s)",pidHash,iPID,fakeRate.Data()));
    //
    fFilteredTreeV0->SetAlias(Form("probCN01_0_%d%s",iPID,sSufix.Data()),Form("AliPIDtools::ComputePIDProbabilityCombined(%d,3,%d,0,3+0,0,%s)",pidHash,iPID,fakeRate.Data()));
    fFilteredTreeV0->SetAlias(Form("probC01_0_%d%s", iPID,sSufix.Data()),Form("AliPIDtools::ComputePIDProbabilityCombined(%d,3,%d,0,3+0,1,%s)",pidHash,iPID,fakeRate.Data()));
    fFilteredTreeV0->SetAlias(Form("probCN03_0_%d%s",iPID,sSufix.Data()),Form("AliPIDtools::ComputePIDProbabilityCombined(%d,9,%d,0,3+0,0,%s)",pidHash,iPID,fakeRate.Data()));
    fFilteredTreeV0->SetAlias(Form("probC03_0_%d%s", iPID,sSufix.Data()),Form("AliPIDtools::ComputePIDProbabilityCombined(%d,9,%d,0,3+0,1,%s)",pidHash,iPID,fakeRate.Data()));
    fFilteredTreeV0->SetAlias(Form("probCN13_0_%d%s",iPID,sSufix.Data()),Form("AliPIDtools::ComputePIDProbabilityCombined(%d,10,%d,0,3+0,0,%s)",pidHash,iPID,fakeRate.Data()));
    fFilteredTreeV0->SetAlias(Form("probC13_0_%d%s", iPID,sSufix.Data()),Form("AliPIDtools::ComputePIDProbabilityCombined(%d,10,%d,0,3+0,1,%s)",pidHash,iPID,fakeRate.Data()));
    //
    fFilteredTreeV0->SetAlias(Form("probCN01_1_%d%s",iPID,sSufix.Data()),Form("AliPIDtools::ComputePIDProbabilityCombined(%d,3,%d,1,3+0,0,%s)",pidHash,iPID,fakeRate.Data()));
    fFilteredTreeV0->SetAlias(Form("probC01_1_%d%s", iPID,sSufix.Data()),Form("AliPIDtools::ComputePIDProbabilityCombined(%d,3,%d,1,3+0,1,%s)",pidHash,iPID,fakeRate.Data()));
    fFilteredTreeV0->SetAlias(Form("probCN03_1_%d%s",iPID,sSufix.Data()),Form("AliPIDtools::ComputePIDProbabilityCombined(%d,9,%d,1,3+0,0,%s)",pidHash,iPID,fakeRate.Data()));
    fFilteredTreeV0->SetAlias(Form("probC03_1_%d%s", iPID,sSufix.Data()),Form("AliPIDtools::ComputePIDProbabilityCombined(%d,9,%d,1,3+0,1,%s)",pidHash,iPID,fakeRate.Data()));
    fFilteredTreeV0->SetAlias(Form("probCN13_1_%d%s",iPID,sSufix.Data()),Form("AliPIDtools::ComputePIDProbabilityCombined(%d,10,%d,1,3+0,0,%s)",pidHash,iPID,fakeRate.Data()));
    fFilteredTreeV0->SetAlias(Form("probC13_1_%d%s", iPID,sSufix.Data()),Form("AliPIDtools::ComputePIDProbabilityCombined(%d,10,%d,1,3+0,1,%s)",pidHash,iPID,fakeRate.Data()));
  }
  return kTRUE;
}

/// RegisterPIDAliasesV0 - to selec clean PID samples using V0 likelihood and PID selection on the leg
/// define "likelihood" for clean PID selection
///        * V0 likelihood
///        * defined in other place
///        * likelihod - combined PID for complementary track
/// * likelihood - combined PID for track except of detectro of interest
/// \param pidHash     - PID hash identifier
/// \param powerLike   - power for the V0 likelihood - defailt 0.6
/// \param powerLegN   - power leg for second track
/// \param powerLeg    - power leg for signal track
/// \param fakeRate    - fake rate query
/// \param suffix      - suffix to add
/// \return
 Bool_t    AliPIDtools::RegisterPIDAliasesV0(Int_t pidHash, Float_t powerLike, Float_t powerLegN, Float_t powerLeg,  const char * fakeR, const char * suffix){
  TTree * treeV0=fFilteredTreeV0;
  if (!treeV0) return kFALSE;
  /// Likelihood for TPC clean sample selection
  {
    // electron like
    treeV0->SetAlias(Form("likeTPCEl0%s",suffix), \
            Form("(ELike**%f)* (AliPIDtools::ComputePIDProbabilityCombined(%d,11,0,1,3,1,%s)**%f)*(AliPIDtools::ComputePIDProbabilityCombined(%d,9,0,0,3,1,%s)**%f)",powerLike,pidHash,fakeR,powerLegN,pidHash,fakeR,powerLeg));
    treeV0->SetAlias(Form("likeTPCEl1%s",suffix), \
            Form("(ELike**%f)* (AliPIDtools::ComputePIDProbabilityCombined(%d,11,0,0,3,1,%s)**%f)*(AliPIDtools::ComputePIDProbabilityCombined(%d,9,0,1,3,1,%s)**%f)",powerLike,pidHash,fakeR,powerLegN,pidHash,fakeR,powerLeg));
    // pion like
    treeV0->SetAlias(Form("likeTPCK0Pi0%s",suffix), \
            Form("(K0Like**%f)*(AliPIDtools::ComputePIDProbabilityCombinedMask(%d,10,6,1,3,1,%s)**%f)*(AliPIDtools::ComputePIDProbabilityCombinedMask(%d,9,7,0,3,1,%s)**%f)",powerLike,pidHash,fakeR,powerLegN,pidHash,fakeR,powerLeg));
    treeV0->SetAlias(Form("likeTPCK0Pi1%s",suffix), \
            Form("(K0Like**%f)*(AliPIDtools::ComputePIDProbabilityCombinedMask(%d,10,6,0,3,1,%s)**%f)*(AliPIDtools::ComputePIDProbabilityCombinedMask(%d,9,7,1,3,1,%s)**%f)",powerLike,pidHash,fakeR,powerLegN,pidHash,fakeR,powerLeg));
    treeV0->SetAlias(Form("likeTPCLPi0%s",suffix), \
            Form("(ALLike**%f)*(AliPIDtools::ComputePIDProbabilityCombinedMask(%d,10,16,1,3,1,%s)**%f)*(AliPIDtools::ComputePIDProbabilityCombinedMask(%d,9,7,0,3,1,%s)**%f)",powerLike,pidHash,fakeR,powerLegN,pidHash,fakeR,powerLeg));
    treeV0->SetAlias(Form("likeTPCLPi1%s",suffix), \
            Form("(LLike**%f)*(AliPIDtools::ComputePIDProbabilityCombinedMask(%d,10,16,0,3,0,%s)**%f)*(AliPIDtools::ComputePIDProbabilityCombinedMask(%d,9,7,1,3,1,%s)**%f)",powerLike,pidHash,fakeR,powerLegN,pidHash,fakeR,powerLeg));
    // proton Like
    treeV0->SetAlias(Form("likeTPCLPr0%s",suffix), \
            Form("(LLike**%f)*(AliPIDtools::ComputePIDProbabilityCombinedMask(%d,10,6,1,3,1,%s)**%f)*(AliPIDtools::ComputePIDProbabilityCombinedMask(%d,9,16,0,3,1,%s)**%f)",powerLike,pidHash,fakeR,powerLegN,pidHash,fakeR,powerLeg));
    treeV0->SetAlias(Form("likeTPCLPr1%s",suffix), \
            Form("(ALLike**%f)*(AliPIDtools::ComputePIDProbabilityCombinedMask(%d,10,6,0,3,0,%s)**%f)*(AliPIDtools::ComputePIDProbabilityCombinedMask(%d,9,16,1,3,1,%s)**%f)",powerLike,pidHash,fakeR,powerLegN,pidHash,fakeR,powerLeg));
  }
  /// Likelihood for TOF clean sample selection
  {
    // electron like
    treeV0->SetAlias(Form("likeTOFEl0%s",suffix), \
            Form("(ELike**%f)* (AliPIDtools::ComputePIDProbabilityCombined(%d,11,0,1,3,1,%s)**%f)*(AliPIDtools::ComputePIDProbabilityCombined(%d,2,0,0,3,1,%s)**%f)",powerLike,pidHash,fakeR,powerLegN,pidHash,fakeR,powerLeg));
    treeV0->SetAlias(Form("likeTOFEl1%s",suffix), \
            Form("(ELike**%f)* (AliPIDtools::ComputePIDProbabilityCombined(%d,11,0,0,3,1,%s)**%f)*(AliPIDtools::ComputePIDProbabilityCombined(%d,2,0,1,3,1,%s)**%f)",powerLike,pidHash,fakeR,powerLegN,pidHash,fakeR,powerLeg));
    // pion like
    treeV0->SetAlias(Form("likeTOFK0Pi0%s",suffix), \
            Form("(K0Like**%f)*(AliPIDtools::ComputePIDProbabilityCombinedMask(%d,10,6,1,3,1,%s)**%f)*(AliPIDtools::ComputePIDProbabilityCombinedMask(%d,3,7,0,3,1,%s)**%f)",powerLike,pidHash,fakeR,powerLegN,pidHash,fakeR,powerLeg));
    treeV0->SetAlias(Form("likeTOFK0Pi1%s",suffix), \
            Form("(K0Like**%f)*(AliPIDtools::ComputePIDProbabilityCombinedMask(%d,10,6,0,3,1,%s)**%f)*(AliPIDtools::ComputePIDProbabilityCombinedMask(%d,3,7,1,3,1,%s)**%f)",powerLike,pidHash,fakeR,powerLegN,pidHash,fakeR,powerLeg));
    treeV0->SetAlias(Form("likeTOFLPi0%s",suffix), \
            Form("(ALLike**%f)*(AliPIDtools::ComputePIDProbabilityCombinedMask(%d,10,16,1,3,1,%s)**%f)*(AliPIDtools::ComputePIDProbabilityCombinedMask(%d,3,7,0,3,1,%s)**%f)",powerLike,pidHash,fakeR,powerLegN,pidHash,fakeR,powerLeg));
    treeV0->SetAlias(Form("likeTOFLPi1%s",suffix), \
            Form("(LLike**%f)*(AliPIDtools::ComputePIDProbabilityCombinedMask(%d,10,16,0,3,0,%s)**%f)*(AliPIDtools::ComputePIDProbabilityCombinedMask(%d,3,7,1,3,1,%s)**%f)",powerLike,pidHash,fakeR,powerLegN,pidHash,fakeR,powerLeg));
    // proton Like
    treeV0->SetAlias(Form("likeTOFLPr0%s",suffix), \
            Form("(LLike**%f)*(AliPIDtools::ComputePIDProbabilityCombinedMask(%d,10,6,1,3,1,%s)**%f)*(AliPIDtools::ComputePIDProbabilityCombinedMask(%d,3,16,0,3,1,%s)**%f)",powerLike,pidHash,fakeR,powerLegN,pidHash,fakeR,powerLeg));
    treeV0->SetAlias(Form("likeTOFLPr1%s",suffix), \
            Form("(ALLike**%f)*(AliPIDtools::ComputePIDProbabilityCombinedMask(%d,10,6,0,3,0,%s)**%f)*(AliPIDtools::ComputePIDProbabilityCombinedMask(%d,3,16,1,3,1,%s)**%f)",powerLike,pidHash,fakeR,powerLegN,pidHash,fakeR,powerLeg));
  }
   /// Likelihood for ITS clean sample selection
  {
    // electron like
    treeV0->SetAlias(Form("likeITSEl0%s",suffix), \
            Form("(ELike**%f)* (AliPIDtools::ComputePIDProbabilityCombined(%d,11,0,1,3,1,%s)**%f)*(AliPIDtools::ComputePIDProbabilityCombined(%d,2,0,0,3,1,%s)**%f)",powerLike,pidHash,fakeR,powerLegN,pidHash,fakeR,powerLeg));
    treeV0->SetAlias(Form("likeITSEl1%s",suffix), \
            Form("(ELike**%f)* (AliPIDtools::ComputePIDProbabilityCombined(%d,11,0,0,3,1,%s)**%f)*(AliPIDtools::ComputePIDProbabilityCombined(%d,2,0,1,3,1,%s)**%f)",powerLike,pidHash,fakeR,powerLegN,pidHash,fakeR,powerLeg));
    // pion like
    treeV0->SetAlias(Form("likeITSK0Pi0%s",suffix), \
            Form("(K0Like**%f)*(AliPIDtools::ComputePIDProbabilityCombinedMask(%d,10,6,1,3,1,%s)**%f)*(AliPIDtools::ComputePIDProbabilityCombinedMask(%d,10,7,0,3,1,%s)**%f)",powerLike,pidHash,fakeR,powerLegN,pidHash,fakeR,powerLeg));
    treeV0->SetAlias(Form("likeITSK0Pi1%s",suffix), \
            Form("(K0Like**%f)*(AliPIDtools::ComputePIDProbabilityCombinedMask(%d,10,6,0,3,1,%s)**%f)*(AliPIDtools::ComputePIDProbabilityCombinedMask(%d,10,7,1,3,1,%s)**%f)",powerLike,pidHash,fakeR,powerLegN,pidHash,fakeR,powerLeg));
    treeV0->SetAlias(Form("likeITSLPi0%s",suffix), \
            Form("(ALLike**%f)*(AliPIDtools::ComputePIDProbabilityCombinedMask(%d,10,16,1,3,1,%s)**%f)*(AliPIDtools::ComputePIDProbabilityCombinedMask(%d,10,7,0,3,1,%s)**%f)",powerLike,pidHash,fakeR,powerLegN,pidHash,fakeR,powerLeg));
    treeV0->SetAlias(Form("likeITSLPi1%s",suffix), \
            Form("(LLike**%f)*(AliPIDtools::ComputePIDProbabilityCombinedMask(%d,10,16,0,3,0,%s)**%f)*(AliPIDtools::ComputePIDProbabilityCombinedMask(%d,10,7,1,3,1,%s)**%f)",powerLike,pidHash,fakeR,powerLegN,pidHash,fakeR,powerLeg));
    // proton Like
    treeV0->SetAlias(Form("likeITSLPr0%s",suffix), \
            Form("(LLike**%f)*(AliPIDtools::ComputePIDProbabilityCombinedMask(%d,10,6,1,3,1,%s)**%f)*(AliPIDtools::ComputePIDProbabilityCombinedMask(%d,10,16,0,3,1,%s)**%f)",powerLike,pidHash,fakeR,powerLegN,pidHash,fakeR,powerLeg));
    treeV0->SetAlias(Form("likeITSLPr1%s",suffix), \
            Form("(ALLike**%f)*(AliPIDtools::ComputePIDProbabilityCombinedMask(%d,10,6,0,3,0,%s)**%f)*(AliPIDtools::ComputePIDProbabilityCombinedMask(%d,10,16,1,3,1,%s)**%f)",powerLike,pidHash,fakeR,powerLegN,pidHash,fakeR,powerLeg));
  }
  /// TRD
    {
    // electron like
    treeV0->SetAlias(Form("likeTRDEl0%s",suffix), \
            Form("(ELike**%f)* (AliPIDtools::ComputePIDProbabilityCombined(%d,11,0,1,3,1,%s)**%f)*(AliPIDtools::ComputePIDProbabilityCombined(%d,2,0,0,3,1,%s)**%f)",powerLike,pidHash,fakeR,powerLegN,pidHash,fakeR,powerLeg));
    treeV0->SetAlias(Form("likeTRDEl1%s",suffix), \
            Form("(ELike**%f)* (AliPIDtools::ComputePIDProbabilityCombined(%d,11,0,0,3,1,%s)**%f)*(AliPIDtools::ComputePIDProbabilityCombined(%d,2,0,1,3,1,%s)**%f)",powerLike,pidHash,fakeR,powerLegN,pidHash,fakeR,powerLeg));
    // pion like
    treeV0->SetAlias(Form("likeTRDK0Pi0%s",suffix), \
            Form("(K0Like**%f)*(AliPIDtools::ComputePIDProbabilityCombinedMask(%d,10,6,1,3,1,%s)**%f)*(AliPIDtools::ComputePIDProbabilityCombinedMask(%d,11,7,0,3,1,%s)**%f)",powerLike,pidHash,fakeR,powerLegN,pidHash,fakeR,powerLeg));
    treeV0->SetAlias(Form("likeTRDK0Pi1%s",suffix), \
            Form("(K0Like**%f)*(AliPIDtools::ComputePIDProbabilityCombinedMask(%d,10,6,0,3,1,%s)**%f)*(AliPIDtools::ComputePIDProbabilityCombinedMask(%d,11,7,1,3,1,%s)**%f)",powerLike,pidHash,fakeR,powerLegN,pidHash,fakeR,powerLeg));
    treeV0->SetAlias(Form("likeTRDLPi0%s",suffix), \
            Form("(ALLike**%f)*(AliPIDtools::ComputePIDProbabilityCombinedMask(%d,10,16,1,3,1,%s)**%f)*(AliPIDtools::ComputePIDProbabilityCombinedMask(%d,11,7,0,3,1,%s)**%f)",powerLike,pidHash,fakeR,powerLegN,pidHash,fakeR,powerLeg));
    treeV0->SetAlias(Form("likeTRDLPi1%s",suffix), \
            Form("(LLike**%f)*(AliPIDtools::ComputePIDProbabilityCombinedMask(%d,10,16,0,3,0,%s)**%f)*(AliPIDtools::ComputePIDProbabilityCombinedMask(%d,11,7,1,3,1,%s)**%f)",powerLike,pidHash,fakeR,powerLegN,pidHash,fakeR,powerLeg));
    // proton Like
    treeV0->SetAlias(Form("likeTRDLPr0%s",suffix), \
            Form("(LLike**%f)*(AliPIDtools::ComputePIDProbabilityCombinedMask(%d,10,6,1,3,1,%s)**%f)*(AliPIDtools::ComputePIDProbabilityCombinedMask(%d,11,16,0,3,1,%s)**%f)",powerLike,pidHash,fakeR,powerLegN,pidHash,fakeR,powerLeg));
    treeV0->SetAlias(Form("likeTRDLPr1%s",suffix), \
            Form("(ALLike**%f)*(AliPIDtools::ComputePIDProbabilityCombinedMask(%d,10,6,0,3,0,%s)**%f)*(AliPIDtools::ComputePIDProbabilityCombinedMask(%d,11,16,1,3,1,%s)**%f)",powerLike,pidHash,fakeR,powerLegN,pidHash,fakeR,powerLeg));
  }
  return kTRUE;
}
