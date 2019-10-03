/**
 * @file   AliSimpleHeaderTask.C
 * @author Christian Holm Christensen <cholm@nbi.dk>
 * @date   Wed Apr 27 16:47:16 2016
 * 
 * @brief  A task to make a simplified AOD header 
 * 
 * @ingroup pwglf_forward_tracklets 
 */
#ifndef ALISIMPLEHEADERTASK_C
#define ALISIMPLEHEADERTASK_C
#include <AliAnalysisTaskSE.h>
#ifndef __CINT__
# include "AliAODSimpleHeader.C"
# include <AliVVertex.h>
# include <AliVertex.h>
# include <AliAnalysisManager.h>
# include <AliVEventHandler.h>
# include <AliInputEventHandler.h>
# include <AliMultSelection.h>
# include <AliAODHandler.h>
# include <AliCollisionGeometry.h>
# include <AliGenEventHeader.h>
# include <AliGenCocktailEventHeader.h>
# include <AliMCEvent.h>
# include <AliLog.h>
# include <AliCentrality.h>
#else
class AliAODSimpleHeader;
class AliVEvent;
class AliMCEvent;
class AliMultSelection;  // Auto-load 
#endif

/**
 * A task to make a simple header in AOD 
 * 
 * @ingroup pwglf_forward_tracklets
 */
class AliSimpleHeaderTask : public AliAnalysisTaskSE
{
public:
  /** 
   * Default constructor - for ROOT I/O only
   */
  AliSimpleHeaderTask();
  /** 
   * Named (user) constructor 
   * 
   * @param name Name of the task 
   */
  AliSimpleHeaderTask(const char* name);
  /** 
   * Copy constructor 
   * 
   * @param o Object to copy from 
   */
  AliSimpleHeaderTask(const AliSimpleHeaderTask& o);
  /** 
   * Assignment operator 
   * 
   * @param o Object to assign from 
   * 
   * @return Reference to this object 
   */
  AliSimpleHeaderTask& operator=(const AliSimpleHeaderTask& o);
  /** 
   * Connect the task to the manager 
   * 
   * @return true on success 
   */
  Bool_t Connect();
  /** 
   * Create an object of this task class, add to manager and connect
   * the input.
   * 
   * @return Pointer to task on success 
   */
  static AliSimpleHeaderTask* Create();
  /** 
   * Set-up outputs 
   * 
   */
  void UserCreateOutputObjects();
  /** 
   * Get the reconstructed interaction point and fill into header 
   * 
   * @param event Event structure 
   * 
   * @return true on success
   */
  Bool_t GetRecIP(AliVEvent* event);
  /** 
   * Get simulation information and fill into header 
   * 
   * @param event Event structure 
   * 
   * @return true on success
   */
  Bool_t GetSim(AliMCEvent* event);
  /** 
   * Get new-style centrality and fill into header 
   * 
   * @param event Event structure 
   * 
   * @return true on success
   */
  Bool_t GetCent(AliVEvent* event);
  /** 
   * Get old-style centrality and fill into header 
   * 
   * @param event Event structure 
   * 
   * @return true on success
   */
  Bool_t GetCentOld(AliVEvent* event);
  /** 
   * Event processing 
   * 
   */
  void UserExec(Option_t* /*option*/);
  /** 
   * Print this task x
   * 
   * @param option 
   */
  void Print(Option_t* option="") const;
  /** the header to fill */
  AliAODSimpleHeader* fHeader; //!

  ClassDef(AliSimpleHeaderTask,1); // Task to make simple header 
};


//____________________________________________________________________
AliSimpleHeaderTask::AliSimpleHeaderTask()
  : AliAnalysisTaskSE(),
    fHeader(0)
{}
//____________________________________________________________________
AliSimpleHeaderTask::AliSimpleHeaderTask(const char* name)
  : AliAnalysisTaskSE(name),
    fHeader(0)
{}
//____________________________________________________________________
AliSimpleHeaderTask::AliSimpleHeaderTask(const AliSimpleHeaderTask& o)
  : AliAnalysisTaskSE(o),
    fHeader(0)
{}
//____________________________________________________________________
AliSimpleHeaderTask&
AliSimpleHeaderTask::operator=(const AliSimpleHeaderTask& o)
{
  return *this;
}
//____________________________________________________________________
Bool_t AliSimpleHeaderTask::Connect()
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    AliError("No analysis manager to connect to.");
    return false;
  }   
  AliAODHandler*      ah = 
    dynamic_cast<AliAODHandler*>(mgr->GetOutputEventHandler());
  if (!ah) {
    AliError("No AOD output handler!");
    return false;
  }
  
  // Add to the manager 
  mgr->AddTask(this);
  
  // Always connect input 
  mgr->ConnectInput(this, 0, mgr->GetCommonInputContainer());
}

//____________________________________________________________________
void AliSimpleHeaderTask::UserCreateOutputObjects()
{
  AliAnalysisManager* am  = AliAnalysisManager::GetAnalysisManager();
  AliAODHandler*      ah  = 
    dynamic_cast<AliAODHandler*>(am->GetOutputEventHandler());
  if (!ah) {
    AliWarning("No AOD output handler set in analysis manager");
    SetZombie();
    return;
  }
  
  fHeader      = new AliAODSimpleHeader;
  TObject* obj = fHeader;
  
  ah->AddBranch("AliAODSimpleHeader", &obj);
}

//____________________________________________________________________
Bool_t AliSimpleHeaderTask::GetRecIP(AliVEvent* event)
{
  const Double_t    maxDispersion = 0.04;
  const Double_t    maxZError     = 0.25;
  const AliVVertex* ip   = event->GetPrimaryVertex();
  if (ip->GetNContributors() <= 0) return false;
  // If this is from the Z vertexer, do some checks 
  if (ip->IsFromVertexerZ()) {
    // Get covariance matrix
    Double_t covar[6];
    ip->GetCovarianceMatrix(covar);
    Double_t sigmaZ = TMath::Sqrt(covar[5]);
    if (sigmaZ >= maxZError) {
      AliWarningF("IPz resolution = %f >= %f", sigmaZ, maxZError);
      return false;
    }
      
    // If this IP doesn not derive from AliVertex, don't check dispersion. 
    if (ip->IsA()->InheritsFrom(AliVertex::Class())) {
      const AliVertex* ipv = static_cast<const AliVertex*>(ip);
      // Dispersion is the parameter used by the vertexer for finding the IP. 
      if (ipv->GetDispersion() >= maxDispersion) {
	AliWarningF("IP dispersion = %f >= %f",
		    ipv->GetDispersion(), maxDispersion);
	return false;
      }
    }
  }
  fHeader->fRecIP.SetXYZ(ip->GetX(), ip->GetY(), ip->GetZ());
}
//____________________________________________________________________
Bool_t AliSimpleHeaderTask::GetSim(AliMCEvent* event)
{
  if (!event) return false;

  TArrayF            genIP(3);
  AliGenEventHeader* genHeader = event->GenEventHeader();
  if (!genHeader) {
    AliWarning("No generator header found in MC event");
    return false;
  }
  genHeader->PrimaryVertex(genIP);
  fHeader->fSimIP.SetXYZ(genIP[0], genIP[1], genIP[2]);
  AliCollisionGeometry* geomHeader =
    dynamic_cast<AliCollisionGeometry*>(genHeader);
  if (!geomHeader && 
      genHeader->IsA()->InheritsFrom(AliGenCocktailEventHeader::Class())) {
    // Cocktail, so we need to find the geometry in one of the
    // cocktail headers.
    AliGenCocktailEventHeader* ctHeader =
      static_cast<AliGenCocktailEventHeader*>(genHeader);     
    TIter next(ctHeader->GetHeaders());
    AliGenEventHeader* subHeader = 0;
    while ((subHeader = static_cast<AliGenEventHeader*>(next()))) {
      geomHeader = dynamic_cast<AliCollisionGeometry*>(subHeader);
      if (geomHeader) break;
    }
  } // end-cocktail
  if (!geomHeader) return true;

  fHeader->fImpactParameter = geomHeader->ImpactParameter();
  fHeader->fReactionPlane   = geomHeader->ReactionPlaneAngle();
  fHeader->fProjectileNpart = geomHeader->ProjectileParticipants();
  fHeader->fTargetNpart     = geomHeader->TargetParticipants();
  fHeader->fNcoll           = (geomHeader->NN()+
			       geomHeader->NNw()+
			       geomHeader->NwN()+
			       geomHeader->NwNw());
  geomHeader->GetNDiffractive(fHeader->fProjectileNsd,
			      fHeader->fTargetNsd,
			      fHeader->fNdd);
  return true;
}
//____________________________________________________________________
Bool_t AliSimpleHeaderTask::GetCent(AliVEvent* event)
{
  AliMultSelection* cent =
    static_cast<AliMultSelection*>(event->FindListObject("MultSelection"));
  if (!cent) {
    AliWarning("No centrality in event");
    return false;
  }
  const Double_t safety  = 1e-3;
  Double_t       centPer = cent->GetMultiplicityPercentile("V0M");
  if      (centPer < -safety)    return false;
  if      (centPer < +safety)    centPer = safety;
  else if (centPer > 100-safety) centPer = 100-safety;

  fHeader->fCent = centPer;
    
  AliMultEstimator* estTracklets = cent->GetEstimator("SPDTracklets");
  if (estTracklets)    
    fHeader->fNTracklets = estTracklets->GetValue();

  return true;
}
//____________________________________________________________________
Bool_t AliSimpleHeaderTask::GetCentOld(AliVEvent* event)
{
  AliCentrality* c = event->GetCentrality();
  if (!c) return false;

  fHeader->fCentOld = c->GetCentralityPercentileUnchecked("V0M");
  return true;
}
//____________________________________________________________________
void  AliSimpleHeaderTask::UserExec(Option_t* /*option*/)
{
  if (!fHeader) return;
    
  AliVEvent* event = InputEvent();
  fHeader->Clear();
  fHeader->fTriggers = fInputHandler->IsEventSelected();
  Bool_t centOK    = GetCent(event);
  Bool_t centOldOK = GetCentOld(event);
  Bool_t recOK     = GetRecIP(event);
  Bool_t simOK     = GetSim(MCEvent());

  if (!simOK || !recOK) return; // Do not insist that event is written

  AliAnalysisManager* am = AliAnalysisManager::GetAnalysisManager();
  AliAODHandler*      ah = 
    dynamic_cast<AliAODHandler*>(am->GetOutputEventHandler());
  if (!ah) return;

  // Mark event for store 
  ah->SetFillAOD(kTRUE);
}

//____________________________________________________________________
void  AliSimpleHeaderTask::Print(Option_t* /*option*/) const
{
  Printf("%s: %s", ClassName(), GetName());
  Printf("  Task to copy content of header to a simplified header");
}
//====================================================================
AliSimpleHeaderTask* AliSimpleHeaderTask::Create()
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("Create","No analysis manager to connect to.");
    return 0;
  }   
  AliAODHandler*      ah = 
    dynamic_cast<AliAODHandler*>(mgr->GetOutputEventHandler());
  if (!ah) {
    ::Error("Create","No AOD output handler!");
    return 0;
  }
  AliSimpleHeaderTask* ret = new AliSimpleHeaderTask("simpleHeader");
  if (!ret->Connect()) return 0;

  return ret;
}

#endif
//
// EOF
//

  
