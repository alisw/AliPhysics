/**************************************************************************
*                                                                        *
* Permission to use, copy, modify and distribute this software and its   *
* documentation strictly for non-commercial purposes is hereby granted   *
* without fee, provided that the above copyright notice appears in all   *
* copies and that both the copyright notice and this permission notice   *
* appear in the supporting documentation. The authors make no claims     *
* about the suitability of this software for any purpose. It is          *
* provided "as is" without express or implied warranty.                  *
**************************************************************************/

/*__|______________________________________________________________________________|
 |                              -----Info(i)-----                                  |
 |                                                                                 |
 |  Cuts Class for single track  efficiecy                                         |
 |                                                                                 |
 |   ESDs<-->AODs (ON/OFF)                                                         |
 |                                                      Authors:                   |
 |_____________________________________________________________________________|___*/


#include "AliVEvent.h"
#include "AliMCEvent.h"
#include "TParticle.h"
#include "AliVParticle.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliVVertex.h"
#include "AliLog.h"
#include "AliAODVertex.h"
#include "AliGenEventHeader.h"
#include "AliInputEventHandler.h"
#include "AliAnalysisManager.h"
#include "AliAODMCHeader.h"
#include "AliAODMCParticle.h"
#include "AliMCParticle.h"
#include "AliVParticle.h"
#include "AliStack.h"
#include "AliMCEventHandler.h"

#include "AliSingleTrackEffCuts.h"


ClassImp(AliSingleTrackEffCuts)

//____________________| Constructor |__________________________
AliSingleTrackEffCuts::AliSingleTrackEffCuts():
  TObject(),
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
  fRequireVtxCuts(kFALSE),
  fMCinfo(0),
  fTriggerMask(AliVEvent::kAny),
  fMinVtxType(0),
  fMinVtxContr(1),
  fMaxVtxZ(10.),
  fnClusITS(0),
  fnClusTPC(0),
  fnClusTOF(0),
  fnClusMUON(0),
  fCutOnZVertexSPD(0)
{
  // Default constructor
}

//____________________| Copy Constructor |__________________________
AliSingleTrackEffCuts::AliSingleTrackEffCuts(const AliSingleTrackEffCuts &source):
  TObject(source),
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
  fRequireVtxCuts(source.fRequireVtxCuts),
  fMCinfo(source.fMCinfo),
  fTriggerMask(source.fTriggerMask),
  fMinVtxType(source.fMinVtxType),
  fMinVtxContr(source.fMinVtxContr),
  fMaxVtxZ(source.fMaxVtxZ),
  fnClusITS(source.fnClusITS),
  fnClusTPC(source.fnClusTPC),
  fnClusTOF(source.fnClusTOF),
  fnClusMUON(source.fnClusMUON),
  fCutOnZVertexSPD(source.fCutOnZVertexSPD)
{
  // Copy constructor
}


//____________________| Source Operator |__________________________
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
  fRequireVtxCuts = source.fRequireVtxCuts;
  fMCinfo = source.fMCinfo;
  fTriggerMask = source.fTriggerMask;
  fMinVtxType = source.fMinVtxType;
  fMinVtxContr = source.fMinVtxContr;
  fMaxVtxZ = source.fMaxVtxZ;
  fnClusITS = source.fnClusITS;
  fnClusTPC = source.fnClusTPC;
  fnClusTOF = source.fnClusTOF;
  fnClusMUON = source.fnClusMUON;
  fCutOnZVertexSPD = source.fCutOnZVertexSPD;

  return *this;
}

//____________________| MC Gen Event Selection|__________________________
Bool_t AliSingleTrackEffCuts::IsMCEventSelected(TObject* obj)
{

        if (!obj) return  kFALSE;
	Bool_t isSelected = kTRUE;
    
	//a) Event interface from ESD/AOD
	AliMCEvent *event;
	AliGenEventHeader *genHeader; //ESDs
	AliAODMCHeader *mcHeader; // AODs
  
	Bool_t isAOD = obj->IsA()->InheritsFrom("AliAODEvent");
	
	if(!isAOD) {
	  event = dynamic_cast<AliMCEvent*>(obj);
	  if (!event) return  kFALSE;
	  genHeader = event->GenEventHeader();
	  if (!genHeader) return  kFALSE;
	} else {
	  AliAODEvent * aodEvent = dynamic_cast<AliAODEvent*> (obj);
	  if (!aodEvent) return  kFALSE;
	  mcHeader = dynamic_cast<AliAODMCHeader*>(aodEvent->GetList()->FindObject(AliAODMCHeader::StdBranchName()));
	  if (!mcHeader) {
	    AliError("Could not find MC Header in AOD");
	    return kFALSE;
	  }
	}
	
	
	//b) Cut on Z Vertex  ( -10cm <= Vz <= 10cm)
	TArrayF vtxPos(3);
	Double_t zMCVertex=-1000;
	if(!isAOD) {

	  genHeader->PrimaryVertex(vtxPos);
	  zMCVertex = vtxPos[2];

	} else {

	  zMCVertex = mcHeader->GetVtxZ();

	}
	
	if( TMath::Abs(zMCVertex)>fMaxVtxZ) isSelected = kFALSE;
	
	// Others cuts ?
	//c) Cut on Multiplicity ?
	//d) Type of MC process non||single||double Diffractive ?
	//e) Trigger Mask and ZDC Energy things?
	
	return isSelected;
}


//____________________| MC Particle Cuts|__________________
//____________________| Step 1. MC Particle Generation Cuts|__________________
Bool_t AliSingleTrackEffCuts::IsMCParticleGenerated(TObject* obj)
{
  
          //  AliMCParticle*  OR AliAODMCParticle*;
          if (!obj) return  kFALSE;
	  if (!obj->InheritsFrom("AliVParticle")) AliError("object must derived from AliVParticle !");
	  AliVParticle* particle = dynamic_cast<AliVParticle *>(obj);
	  if ( !particle ) return kFALSE;
	  
	  Bool_t isSelected = kTRUE;
	  
       	  // Pdg Code
	  if( fIsPdgCode && TMath::Abs( particle->PdgCode() )!= fPdgCode) isSelected = kFALSE;
	  
	  // Charge selection
	  if(fIsCharged && !(particle->Charge()!=0)) isSelected = kFALSE;
	 

	  // Physical Primary 
	  if(!fisAOD) { 
	    //cout <<" I am selecting Physics Primary only in ESDs"<< endl;
	    AliMCEventHandler* mcinfo = (AliMCEventHandler*) (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());  	     
	    fMCinfo = mcinfo->MCEvent();
	    AliStack* stack = ((AliMCEvent*)fMCinfo)->Stack();

	    AliMCParticle* mcPart = dynamic_cast<AliMCParticle *>(obj);	     
	    if (!stack->IsPhysicalPrimary(mcPart->GetLabel())) isSelected = kFALSE;
	    
	  } else {
	    // cout <<" I am selecting Physics Primary only in AODs"<< endl;
	    AliAODMCParticle* mcPart = dynamic_cast<AliAODMCParticle *>(obj);	    
	    if (!mcPart->IsPhysicalPrimary()) isSelected = kFALSE;

	  }
	    
	 
	  return isSelected;

}



//____________________| Step 2. MC Particle Kinematics Cuts|__________________
Bool_t AliSingleTrackEffCuts::IsMCParticleInKineAcceptance(TObject *obj)
{

          if (!obj) return  kFALSE;
	  if (!obj->InheritsFrom("AliVParticle")) AliError("object must derived from AliVParticle !");
	  AliVParticle* particle = dynamic_cast<AliVParticle *>(obj);
	  if (!particle) return kFALSE;
	  
	  Bool_t isSelected = kTRUE;
	  
	  // Cuts on eta
	  if( particle->Eta() < fEtaMin || particle->Eta() > fEtaMax ) isSelected = kFALSE;
	  
	  // Cuts on y
	  // if( particle->Y() < fYMin || particle->Y() > fYMax ) isSelected = kFALSE;
	  
	  // Cuts on pt
	  if( particle->Pt() < fPtMin || particle->Pt() > fPtMax ) isSelected = kFALSE;
	  
	  return isSelected;
}



//____________________| Step 3. MC Particle Track Ref Cuts "ONLY FOR ESDs"|__________________
Bool_t AliSingleTrackEffCuts::IsMCParticleInReconstructable(TObject *obj)
{
  
          if (!obj) return kFALSE;

	  TString className(obj->ClassName());
	  if (className.CompareTo("AliMCParticle") != 0) {
	    AliError("obj must point to an AliMCParticle !");
	    return kTRUE; // <===================================== FIX ME !!
	  }
	  
	  AliMCParticle * part = dynamic_cast<AliMCParticle*>(obj);
	  if(!part) return kFALSE;
	  
	  Bool_t isSelected = kTRUE;
	  
	  Int_t nHitsITS=0, nHitsTPC=0, nHitsTRD=0, nHitsTOF=0, nHitsMUON=0 ;
	  for (Int_t iTrackRef=0; iTrackRef<part->GetNumberOfTrackReferences(); iTrackRef++) {
	    AliTrackReference * trackRef = part->GetTrackReference(iTrackRef);
	    if(trackRef){
	      Int_t detectorId = trackRef->DetectorId();
	      switch(detectorId) {
	      case AliTrackReference::kITS  : nHitsITS++  ;
	      case AliTrackReference::kTPC  : nHitsTPC++  ;
	      case AliTrackReference::kTRD  : nHitsTRD++  ; 
	      case AliTrackReference::kTOF  : nHitsTOF++  ;
	      case AliTrackReference::kMUON : nHitsMUON++ ; 
		//      default : break ;
	      }
	    }
	  }
	  
	  if (nHitsITS<fnClusITS) isSelected = kFALSE;
	  if (nHitsTPC<fnClusTPC) isSelected = kFALSE;
	  if (nHitsTOF<fnClusTOF) isSelected = kFALSE;
	  if (nHitsMUON<fnClusMUON) isSelected = kFALSE;
	  
	  return isSelected;


}


//____________________| MC Reco Event Selection |__________________
Bool_t AliSingleTrackEffCuts::IsRecoEventSelected(TObject* obj)
{
 
  //  AliInfo("checking reco event");
  
      AliVEvent *event = dynamic_cast<AliVEvent*>(obj);
	  if (!event) return kFALSE;
	  
	  Bool_t isSelected = kTRUE;
	  
	  //a) Physics selection
	  UInt_t trigFired = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
	  Bool_t isEvtSelected = (trigFired & fTriggerMask);
	  if(!isEvtSelected) isSelected = kFALSE;
	  
	  //b) Vertex selection
	  Bool_t isVtxSelected = IsVertexSelected(event);
	  if(!isVtxSelected) isSelected = kFALSE;
	  
	  
	  //c) Cut on Multiplicity ?
	  //d) others cuts?
	  
	  return isSelected;


}


//____________________| MC Reco Track Cuts|__________________
//____________________| Step 4. MC Particle Generation Cuts|__________________
Bool_t AliSingleTrackEffCuts::IsRecoParticleKineAcceptance(TObject *obj)
{
  // AliInfo("checking particle reco info");
    Bool_t isSelected = kTRUE;
    
    AliVParticle *track = dynamic_cast<AliVParticle*>(obj);
    
    
    // Cuts on eta
    if( track->Eta() < fEtaMin || track->Eta() > fEtaMax ) isSelected = kFALSE;
    
    // Cuts on y
    // if( track->Y()   < fYMin   || track->Y()   > fYMax )   isSelected = kFALSE;
    
  // Cuts on pt
    if( track->Pt()  < fPtMin  || track->Pt()  > fPtMax ) isSelected = kFALSE;
    


    return isSelected;


}

//____________________| Function for Vertex finding at Reco |__________________
Bool_t AliSingleTrackEffCuts::IsVertexSelected(AliVEvent *event)
{


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
	  //  std::cout << " vtx tittle "<< title<< std::endl;
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
	  const AliVVertex *vSPD;
	  if(isAOD) {
	    vSPD = ((AliAODEvent*)event)->GetPrimaryVertexSPD(); 
	  }else {    
	    vSPD = ((AliESDEvent*)event)->GetPrimaryVertexSPD();  }

	  if(fCutOnZVertexSPD==1 && vSPD && vSPD->GetNContributors()>=fMinVtxContr) {
	    
	    if(TMath::Abs(vSPD->GetZ()-vertex->GetZ())>0.5) accept = kFALSE;
	  }
	  
	  return accept;



}
