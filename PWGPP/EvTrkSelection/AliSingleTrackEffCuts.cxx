#include "AliVEvent.h"
#include "AliMCEvent.h"
#include "TParticle.h"
#include "AliVParticle.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliVVertex.h"
#include "AliLog.h"
#include "AliGenEventHeader.h"
#include "AliInputEventHandler.h"
#include "AliAnalysisManager.h"
#include "AliAODMCHeader.h"
#include "AliAODMCParticle.h"
#include "AliMCParticle.h"
#include "AliVParticle.h"
#include "AliMCEventHandler.h"
#include "AliPIDResponse.h"
#include "AliPIDCombined.h"
#include "AliSingleTrackEffCuts.h"

using std::cout;
using std::endl;

ClassImp(AliSingleTrackEffCuts)

//_____________________________________________________
AliSingleTrackEffCuts::AliSingleTrackEffCuts():
AliAnalysisCuts(),
  fisAOD(kTRUE),
  fIsPdgCode(kFALSE),
  fPdgCode(0),
  fEtaMin(-12),
  fEtaMax(12),
  fYMin(-12),
  fYMax(12),
  fPtMin(-15),
  fPtMax(15),
  fIsCharged(kTRUE),
  fTriggerMask(AliVEvent::kAny),
  fMinVtxType(0),
  fMinVtxContr(1),
  fMaxVtxZ(10.),
  fCutOnZVertexSPD(0),
  fnClusITS(0),
  fnClusTPC(0),
  fnClusTOF(0),
  fnClusMUON(0),
  fusePid(kFALSE),
  fParticlePid(AliPID::kPion),
  fuseTPCPid(kTRUE),
  fnPTPCBins(0),
  fnPTPCBinLimits(0),
  fPTPCBinLimits(0),
  fnSigmaTPC(0),
  fPmaxTPC(9999.),
  fuseTOFPid(kTRUE),
  fnPTOFBins(0),
  fnPTOFBinLimits(0),
  fPTOFBinLimits(0),
  fnSigmaTOF(0),
  fPmaxTOF(9999.),
  fuseCombinPid(0),
  fThreshold(0.3)
{
  //
  // Default constructor
  //
}

//________________________________________________________________________________
AliSingleTrackEffCuts::AliSingleTrackEffCuts(const char* name, const char* title):
AliAnalysisCuts(name,title),
  fisAOD(kTRUE),
  fIsPdgCode(kFALSE),
  fPdgCode(0),
  fEtaMin(-12),
  fEtaMax(12),
  fYMin(-12),
  fYMax(12),
  fPtMin(-15),
  fPtMax(15),
  fIsCharged(kTRUE),
  fTriggerMask(AliVEvent::kAny),
  fMinVtxType(0),
  fMinVtxContr(1),
  fMaxVtxZ(10.),
  fCutOnZVertexSPD(0),
  fnClusITS(0),
  fnClusTPC(0),
  fnClusTOF(0),
  fnClusMUON(0),
  fusePid(kFALSE),
  fParticlePid(AliPID::kPion),
  fuseTPCPid(kTRUE),
  fnPTPCBins(0),
  fnPTPCBinLimits(0),
  fPTPCBinLimits(0),
  fnSigmaTPC(0),
  fPmaxTPC(9999.),
  fuseTOFPid(kTRUE),
  fnPTOFBins(0),
  fnPTOFBinLimits(0),
  fPTOFBinLimits(0),
  fnSigmaTOF(0),
  fPmaxTOF(9999.),
  fuseCombinPid(0),
  fThreshold(0.3)
{
  //
  // Default constructor
  //
}

//_________________________________________________________________________________
AliSingleTrackEffCuts::AliSingleTrackEffCuts(const AliSingleTrackEffCuts &source):
  AliAnalysisCuts(source),
  fisAOD(source.fisAOD),
  fIsPdgCode(source.fIsPdgCode),
  fPdgCode(source.fPdgCode),
  fEtaMin(source.fEtaMin),
  fEtaMax(source.fEtaMax),
  fYMin(source.fYMin),
  fYMax(source.fYMax),
  fPtMin(source.fPtMin),
  fPtMax(source.fPtMax),
  fIsCharged(source.fIsCharged),
  fTriggerMask(source.fTriggerMask),
  fMinVtxType(source.fMinVtxType),
  fMinVtxContr(source.fMinVtxContr),
  fMaxVtxZ(source.fMaxVtxZ),
  fCutOnZVertexSPD(source.fCutOnZVertexSPD),
  fnClusITS(source.fnClusITS),
  fnClusTPC(source.fnClusTPC),
  fnClusTOF(source.fnClusTOF),
  fnClusMUON(source.fnClusMUON),
  fusePid(source.fusePid),
  fParticlePid(source.fParticlePid),
  fuseTPCPid(source.fuseTPCPid),
  fnPTPCBins(source.fnPTPCBins),
  fnPTPCBinLimits(source.fnPTPCBinLimits),
  fPTPCBinLimits(source.fPTPCBinLimits),
  fnSigmaTPC(source.fnSigmaTPC),
  fPmaxTPC(source.fPmaxTPC),
  fuseTOFPid(source.fuseTOFPid),
  fnPTOFBins(source.fnPTOFBins),
  fnPTOFBinLimits(source.fnPTOFBinLimits),
  fPTOFBinLimits(source.fPTOFBinLimits),
  fnSigmaTOF(source.fnSigmaTOF),
  fPmaxTOF(source.fPmaxTOF),
  fuseCombinPid(source.fuseCombinPid),
  fThreshold(source.fThreshold)
{
  //
  // Copy constructor
  //
}


//_________________________________________________________________________________________
AliSingleTrackEffCuts &AliSingleTrackEffCuts::operator=(const AliSingleTrackEffCuts &source)
{
  //
  // assignment operator
  //
  if(&source == this) return *this;

  fisAOD = source.fisAOD;

  fIsPdgCode = source.fIsPdgCode;
  fPdgCode = source.fPdgCode;
  fEtaMin = source.fEtaMin;
  fEtaMax = source.fEtaMax;
  fYMin = source.fYMin;
  fYMax = source.fYMax;
  fPtMin = source.fPtMin;
  fPtMax = source.fPtMax;
  fIsCharged = source.fIsCharged;

  fTriggerMask = source.fTriggerMask;
  fMinVtxType = source.fMinVtxType;
  fMinVtxContr = source.fMinVtxContr;
  fMaxVtxZ = source.fMaxVtxZ;
  fCutOnZVertexSPD = source.fCutOnZVertexSPD;

  fnClusITS = source.fnClusITS;
  fnClusTPC = source.fnClusTPC;
  fnClusTOF = source.fnClusTOF;
  fnClusMUON = source.fnClusMUON;

  fusePid = source.fusePid;
  fParticlePid = source.fParticlePid;
  fuseTPCPid = source.fuseTPCPid;
  fnPTPCBins = source.fnPTPCBins;
  fnPTPCBinLimits = source.fnPTPCBinLimits;
  if(source.fPTPCBinLimits && source.fnSigmaTPC) 
    SetTPCSigmaPtBins(source.fnPTPCBins,source.fPTPCBinLimits,source.fnSigmaTPC);
  fPmaxTPC = source.fPmaxTPC;
  fuseTOFPid = source.fuseTOFPid;
  fnPTOFBins = source.fnPTOFBins;
  fnPTOFBinLimits = source.fnPTOFBinLimits;
  if(source.fPTOFBinLimits && source.fnSigmaTOF) 
    SetTOFSigmaPtBins(source.fnPTOFBins,source.fPTOFBinLimits,source.fnSigmaTOF);
  fPmaxTOF = source.fPmaxTOF;
  fuseCombinPid = source.fuseCombinPid;
  fThreshold = source.fThreshold;
  return *this;
}

//______________________________________________
AliSingleTrackEffCuts::~AliSingleTrackEffCuts()
{
  //
  // Destructor
  //

  if(fPTPCBinLimits){
    delete [] fPTPCBinLimits;
    fPTPCBinLimits=NULL;
  }
  if(fnSigmaTPC){
    delete [] fnSigmaTPC;
    fnSigmaTPC=NULL;
  }
  if(fPTOFBinLimits){
    delete [] fPTOFBinLimits;
    fPTOFBinLimits=NULL;
  }
  if(fnSigmaTOF){
    delete [] fnSigmaTOF;
    fnSigmaTOF=NULL;
  }

}

//______________________________________________
Bool_t AliSingleTrackEffCuts::IsMCEventSelected(TObject* obj)
{
  //
  // Event selection at MC level 
  //  |zvtx|<=fMaxVtxZ
  //
  if(!obj) return  kFALSE;
  Bool_t isSelected = kTRUE;
  
  AliMCEvent *event=0;
  AliGenEventHeader *genHeader=0; //ESDs
  AliAODMCHeader *mcHeader=0; // AODs
  Bool_t isAOD = obj->IsA()->InheritsFrom("AliAODEvent");

  if(!isAOD) {
    event = dynamic_cast<AliMCEvent*>(obj);
    if (!event) return  kFALSE;
    genHeader = event->GenEventHeader();
    if (!genHeader) return  kFALSE;
  } else {
    AliAODEvent *aodEvent = dynamic_cast<AliAODEvent*> (obj);
    if (!aodEvent) return  kFALSE;
    mcHeader = dynamic_cast<AliAODMCHeader*>(aodEvent->GetList()->FindObject(AliAODMCHeader::StdBranchName()));
    if (!mcHeader) {
      AliError("Could not find MC Header in AOD");
      return kFALSE;
    }

    // Check for not null trigger mask to evict pPb MC buggy-events
    Int_t runnumber = aodEvent->GetRunNumber();
    if(aodEvent->GetTriggerMask()==0 && 
       (runnumber>=195344 && runnumber<=195677)){
      AliDebug(3,"Event rejected because of null trigger mask");
      return kFALSE;
    }
  }


  // Cut on Z Vertex ( |z-vtx| <= fMaxVtxZ )
  TArrayF vtxPos(3);
  Double_t zMCVertex=-1000;
  if(!isAOD) {
    genHeader->PrimaryVertex(vtxPos);
    zMCVertex = vtxPos[2];
  } else {
    zMCVertex = mcHeader->GetVtxZ();
  }  
  if( TMath::Abs(zMCVertex)>fMaxVtxZ ) isSelected = kFALSE;


  return isSelected;
}

//__________________________________________________________
Bool_t AliSingleTrackEffCuts::IsMCParticleGenerated(TObject* obj)
{
  //
  // Check generated particles (primary, charged, pdgcode) 
  //
  
  if(!obj) return kFALSE;
  if(!obj->InheritsFrom("AliVParticle")) {
    AliError("object must derived from AliVParticle !");
    return kFALSE;
  }
  AliVParticle* particle = dynamic_cast<AliVParticle *>(obj);
  if(!particle) return kFALSE;

  Bool_t isSelected = kTRUE;

  // Check particle Pdg Code
  if(fIsPdgCode && TMath::Abs( particle->PdgCode() )!= fPdgCode) isSelected = kFALSE;

  // Charge selection
  if(fIsCharged && (particle->Charge()==0)) isSelected = kFALSE;

  // Selection of Physical Primary particles
  if(!fisAOD) { // check on ESDs
    AliMCEventHandler* mcinfo = (AliMCEventHandler*) (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());  	     
    AliMCEvent* mcevent = mcinfo->MCEvent();
    AliMCParticle* mcPart = dynamic_cast<AliMCParticle *>(obj);	
    if (!mcPart) return kFALSE;
    if(!mcevent->IsPhysicalPrimary(mcPart->GetLabel())) {
      isSelected = kFALSE;
    }
  } else { // Check on AODs
    AliAODMCParticle* mcPart = dynamic_cast<AliAODMCParticle *>(obj);	    
    if (!mcPart) return kFALSE;
    if(!mcPart->IsPhysicalPrimary()) isSelected = kFALSE;
  }

  return isSelected;
}


//_____________________________________________________________________
Bool_t AliSingleTrackEffCuts::IsMCParticleInKineAcceptance(TObject *obj)
{
  //
  // Check generated particles (eta, y, pt)
  //

  if(!obj) return  kFALSE;
  if(!obj->InheritsFrom("AliVParticle")) AliError("object must derived from AliVParticle !");
  AliVParticle* particle = dynamic_cast<AliVParticle *>(obj);
  if(!particle) return kFALSE;

  Bool_t isSelected = kTRUE;

  // Cut on eta
  if(particle->Eta()<fEtaMin || particle->Eta()>fEtaMax) isSelected = kFALSE;

  // Cut on y
  Double_t energy = particle->E();
  Double_t pz2e = energy>0 ? particle->Pz()/energy : 1.0;
  Double_t particleY = (TMath::Abs(pz2e)<(1-1.e-12)) ? 0.5*TMath::Log( (1+pz2e)/(1-pz2e) ) : 1e6;
 
  if(particleY<fYMin || particleY>fYMax) isSelected = kFALSE;

  // Cut on pt
  if(particle->Pt()<fPtMin || particle->Pt()>fPtMax) isSelected = kFALSE;

  return isSelected;
}


//_______________________________________________________________________
Bool_t AliSingleTrackEffCuts::IsMCParticleInReconstructable(TObject *obj)
{
  //
  // Check if particle has left enough hits in the detectors (only at ESD level)
  //

  if(!obj) return kFALSE;
  TString className(obj->ClassName());
  if (className.CompareTo("AliMCParticle") != 0) {
    AliError("obj must point to an AliMCParticle !");
    return kTRUE; // <===================================== FIX ME !!
  }

  AliMCParticle * part = dynamic_cast<AliMCParticle*>(obj);
  if(!part) return kFALSE;

  Bool_t isSelected = kTRUE;

  Int_t nHitsITS=0, nHitsTPC=0, nHitsTRD=0, nHitsTOF=0, nHitsMUON=0;
  for(Int_t iTrackRef=0; iTrackRef<part->GetNumberOfTrackReferences(); iTrackRef++) {
    AliTrackReference * trackRef = part->GetTrackReference(iTrackRef);
    if(trackRef){
      Int_t detectorId = trackRef->DetectorId();
      switch(detectorId) {
      case AliTrackReference::kITS  : 
	nHitsITS++; 
	break;
      case AliTrackReference::kTPC  : 
	nHitsTPC++; 
	break;
      case AliTrackReference::kTRD  : 
	nHitsTRD++; 
	break;
      case AliTrackReference::kTOF  : 
	nHitsTOF++; 
	break;
      case AliTrackReference::kMUON : 
	nHitsMUON++; 
	break;
      default : break;
      }
    }
  }

  if(nHitsITS<fnClusITS) isSelected = kFALSE;
  if(nHitsTPC<fnClusTPC) isSelected = kFALSE;
  if(nHitsTOF<fnClusTOF) isSelected = kFALSE;
  if(nHitsMUON<fnClusMUON) isSelected = kFALSE;

  return isSelected;
}


//_____________________________________________________________
Bool_t AliSingleTrackEffCuts::IsRecoEventSelected(TObject* obj)
{
  //
  // Event selection at reconstructed level (trigger, zvtx)
  //
 
  AliVEvent *event = dynamic_cast<AliVEvent*>(obj);
  if(!event) return kFALSE;

  Bool_t isSelected = kTRUE;

  // Check if event is accepted by the Physics selection
  UInt_t trigFired = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
  Bool_t isEvtSelected = (trigFired & fTriggerMask);
  if(!isEvtSelected) isSelected = kFALSE;

  // Vertex selection
  Bool_t isVtxSelected = IsVertexSelected(event);
  if(!isVtxSelected) isSelected = kFALSE;

  return isSelected;
}


//______________________________________________________________________
Bool_t AliSingleTrackEffCuts::IsRecoParticleKineAcceptance(TObject *obj)
{
  //
  // Check if reconstructed particle is in the acceptance (eta, y, pt)
  //

  Bool_t isSelected = kTRUE;

  AliVParticle *track = dynamic_cast<AliVParticle*>(obj);
  if (!track) return kFALSE;
 
  // Cut on eta
  if(track->Eta()<fEtaMin || track->Eta()>fEtaMax) isSelected = kFALSE;

  // Cut on y
  if(track->Y()<fYMin || track->Y()>fYMax) isSelected = kFALSE;

  // Cut on pt
  if(track->Pt()<fPtMin || track->Pt() >fPtMax) isSelected = kFALSE;

  return isSelected;
}

//_______________________________________________________________
Bool_t AliSingleTrackEffCuts::IsVertexSelected(AliVEvent *event)
{
  //
  // Check if the reconstructed vertex is selected
  //

  Bool_t accept = kTRUE;
  Bool_t isAOD = event->IsA()->InheritsFrom("AliAODEvent");

  const AliVVertex *vertex = event->GetPrimaryVertex();
  if(!vertex){
    accept = kFALSE;
    AliInfo("no vtx");
    return accept;
  }

  // Cut on vertex type  
  TString title=vertex->GetTitle();
  if(title.Contains("Z") && fMinVtxType>1){
    accept=kFALSE;
  } else if(title.Contains("3D") && fMinVtxType>2){
    accept=kFALSE;
  }
  
  // cut on minimum number of contributors
  if(vertex->GetNContributors()<fMinVtxContr){
    AliInfo(Form("too few contributors %d",vertex->GetNContributors()));
    accept=kFALSE;
  }

  // cut on absolute |z| of the vertex
  if(TMath::Abs(vertex->GetZ())>fMaxVtxZ) {
    AliInfo("outside the Vtx range");
    accept=kFALSE;
  }

  // cut on distance of SPD and TRK vertexes
  if(fCutOnZVertexSPD==1) {
    const AliVVertex *vSPD = NULL;
    if(isAOD) {
      vSPD = ((AliAODEvent*)event)->GetPrimaryVertexSPD(); 
    }else {    
      vSPD = ((AliESDEvent*)event)->GetPrimaryVertexSPD();  
    }

    if(vSPD && vSPD->GetNContributors()>=fMinVtxContr) {      
      if(TMath::Abs(vSPD->GetZ()-vertex->GetZ())>0.5) accept = kFALSE;
    }
  }

  return accept;
}

//_________________________________________________________________________________________________
void AliSingleTrackEffCuts::SetTPCSigmaPtBins(Int_t nPtBins, Float_t *pBinLimits, Float_t *sigmaBin)
{
  //
  // Set TPC Pid number-of-Sigma in P bins
  //

  // Set the pt bins
  if(fPTPCBinLimits) {
    delete [] fPTPCBinLimits;
    fPTPCBinLimits = NULL;
    printf("Changing the TPC cut P bins\n");
  }
  if(fnSigmaTPC) {
    delete [] fnSigmaTPC;
    fnSigmaTPC = NULL;
    printf("Changing the TPC Sigma cut per p bin\n");
  }

  fnPTPCBins = nPtBins;
  fnPTPCBinLimits = nPtBins+1;
  fPTPCBinLimits = new Float_t[fnPTPCBinLimits];
  for(Int_t ib=0; ib<nPtBins+1; ib++) fPTPCBinLimits[ib]=pBinLimits[ib];
  fnSigmaTPC = new Float_t[fnPTPCBins];
  for(Int_t ib=0; ib<nPtBins; ib++) fnSigmaTPC[ib]=sigmaBin[ib];

  return;
}
//_____________________________________________________________________________________________________
void AliSingleTrackEffCuts::SetTOFSigmaPtBins(Int_t nPtBins, Float_t *pBinLimits, Float_t *sigmaBin)
{
  //
  // Set TOF Pid number-of-Sigma in P bins
  //

  // Set the pt bins
  if(fPTOFBinLimits) {
    delete [] fPTOFBinLimits;
    fPTOFBinLimits = NULL;
    printf("Changing the TOF cut P bins\n");
  }
  if(fnSigmaTOF) {
    delete [] fnSigmaTOF;
    fnSigmaTOF = NULL;
    printf("Changing the TOF Sigma cut per p bin\n");
  }

  fnPTOFBins = nPtBins;
  fnPTOFBinLimits = nPtBins+1;
  fPTOFBinLimits = new Float_t[fnPTOFBinLimits];
  for(Int_t ib=0; ib<nPtBins+1; ib++) fPTOFBinLimits[ib]=pBinLimits[ib];
  fnSigmaTOF = new Float_t[fnPTOFBins];
  for(Int_t ib=0; ib<nPtBins; ib++) fnSigmaTOF[ib]=sigmaBin[ib];

  return;
}

//___________________________________________________________________________
Bool_t AliSingleTrackEffCuts::CheckTPCPIDStatus(AliAODTrack *track) const{
  //
  // Check TPC PID status
  //

  if ((track->GetStatus() & AliESDtrack::kTPCin)==0) return kFALSE;
  UShort_t nTPCClus=track->GetTPCClusterMap().CountBits();
  if (nTPCClus<70) return kFALSE;
  return kTRUE;
}

//___________________________________________________________________________
Bool_t AliSingleTrackEffCuts::CheckTOFPIDStatus(AliAODTrack *track) const{
  //
  // Check TOC PID status
  //

  if ((track->GetStatus()&AliESDtrack::kTOFout)==0)   return kFALSE;
  if ((track->GetStatus()&AliESDtrack::kTIME)==0)     return kFALSE;
  if ((track->GetStatus()&AliESDtrack::kTOFpid)==0)   return kFALSE;
  if ((track->GetStatus()&AliESDtrack::kTOFmismatch)!=0) return kFALSE;
  return kTRUE;
}

//___________________________________________________________________________
Bool_t AliSingleTrackEffCuts::IsRecoParticlePID(TObject *obj)
{
  //
  // Check Particle PID (AOD only for now!)
  //

  Bool_t isSelected = kFALSE;
  Bool_t isAOD = obj->IsA()->InheritsFrom("AliAODTrack");

  if(!isAOD || !GetUsePid()) return isSelected;
  if(!fuseTPCPid && !fuseTOFPid) return isSelected;

  AliAODTrack *track = dynamic_cast<AliAODTrack*>(obj);
  if(!track) { cout<<"No track found"<<endl; return isSelected; }

  // AliAODPid *pid = track->GetDetPid();
  // if(!pid) { cout<<"No AliAODPid found"<<endl; return isSelected; }

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler *inputHandler=(AliInputEventHandler*)mgr->GetInputEventHandler();
  AliPIDResponse *pidResp=inputHandler->GetPIDResponse();
  if(!pidResp) { cout<<"No PidResponse found"<<endl; return isSelected;}

  Double_t pPart = track->P();
  
  // Check detector status
  Bool_t okTPC = CheckTPCPIDStatus(track);
  Bool_t okTOF = CheckTOFPIDStatus(track);

  //
  // Bayesian PID (fuseCombinPid>0)
  // maximum probability selection(fuseCombinPid==1)
  if(fuseCombinPid>0){
    isSelected = false;
    AliPIDCombined *combin = new AliPIDCombined();
    combin->SetDetectorMask(AliPIDResponse::kDetTPC|AliPIDResponse::kDetTOF);
    combin->SetDefaultTPCPriors();
    Double_t prob[AliPID::kSPECIES];
    combin->ComputeProbabilities(track,pidResp,prob);

    if(fuseCombinPid==AliSingleTrackEffCuts::kMaximumBayesianProb) {
    // maximum probability selection
      if( prob[(AliPID::EParticleType)fParticlePid]>0. &&
	  (TMath::MaxElement(AliPID::kSPECIES,prob) == prob[(AliPID::EParticleType)fParticlePid]) ){
	isSelected = true;
      }
    }
    else if(fuseCombinPid==AliSingleTrackEffCuts::kThresholdBayesianProb) {
      // probability threshold selection
      if( prob[(AliPID::EParticleType)fParticlePid] > GetPIDThreshold()){
	isSelected = true;
      }
    }
    else {
      AliWarning("Method for Bayesian PID not yet implemented");
    }

    delete combin;
    return isSelected;
  }


  //  Check Number of Sigmas
  Double_t nsigmaTPC=pidResp->NumberOfSigmasTPC((AliVParticle*)track,(AliPID::EParticleType)fParticlePid);
  Double_t nsigmaTOF=pidResp->NumberOfSigmasTOF((AliVParticle*)track,(AliPID::EParticleType)fParticlePid);
  Bool_t isTPCPid=false, isTOFPid=false;

  // If use TPC and TPC infos are ok, check whether the sigma is ok in the given p range
  if(fuseTPCPid && okTPC) {
    for(Int_t j=0; j<fnPTPCBins; j++) {
      //      cout<<" checking bin: ("<<fPTPCBinLimits[j]<<","<<fPTPCBinLimits[j+1]<<") should be nsigma < "<<fnSigmaTPC[j]<<endl;
      if ((pPart>fPTPCBinLimits[j]) && (pPart<fPTPCBinLimits[j+1]) && nsigmaTPC<fnSigmaTPC[j]) isTPCPid=true;
    }
    if(pPart>fPmaxTPC) isTPCPid=true;
  }

  // If use TPC and TPC infos are ok, check whether the sigma is ok in the given p range
  if(fuseTOFPid && okTOF) {
    for(Int_t j=0; j<fnPTOFBins; j++) {
      if ((pPart>fPTOFBinLimits[j]) && (pPart<fPTOFBinLimits[j+1]) && nsigmaTOF<fnSigmaTOF[j]) isTPCPid=true;
    }
    if(pPart>fPmaxTOF) isTOFPid=true;
  }


  isSelected = (isTPCPid || isTOFPid) ? true : false;

  return isSelected;
}
