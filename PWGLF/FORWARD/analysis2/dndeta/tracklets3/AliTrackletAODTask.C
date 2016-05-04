/**
 * @file   AliTrackletAODTask.C
 * @author Christian Holm Christensen <cholm@nbi.dk>
 * @date   Wed Apr 27 16:48:04 2016
 * 
 * @brief  Tasks to make tracklet AOD output
 * 
 * @ingroup pwglf_forward_tracklets
 * 
 */
#include <AliAnalysisTaskSE.h>
#ifndef __CINT__
#include "AliAODTracklet.C"
#include "AliVVertex.h"
#include <AliVEvent.h>
#include <AliAODEvent.h>
#include <AliVVertex.h>
#include <AliVertex.h>
#include <AliVMultiplicity.h>
#include <AliMultiplicity.h>
#include <AliAnalysisManager.h>
#include <AliVEventHandler.h>
#include <AliESDInputHandlerRP.h>
#include <AliITSMultRecBg.h>
#include <AliCDBManager.h>
#include <AliCDBEntry.h>
#include <AliCDBPath.h>
#include <AliCDBId.h>
#include <AliGeomManager.h>
#include <AliAODHandler.h>
#include <TGeoManager.h>
#include <TGeoGlobalMagField.h>
#else
class AliAODTracklet;
class AliVEvent;
class AliVVertex;
class AliMultiplicity;
class AliAnalysisDataContainer;
class AliITSMultRecBg;
class AliGeomManager;           // Auto-load
class AliCDBManager;            // Auto-load
class AliMultSelection;         // Auto-load
class TClonesArray;
class TTree;
class TGeoGlobalMagField;       // Auto-load
class TGeoManager;              // Auto-load
class TParticle;
#endif


//====================================================================
/**
 * Store tracklets on AOD 
 * 
 * @ingroup pwglf_forward_tracklets
 */
class AliTrackletAODTask : public AliAnalysisTaskSE
{
public:
  /** 
   * Default constructor - for ROOT I/O only 
   */
  AliTrackletAODTask();
  /** 
   * Named user constructor 
   * 
   * @param name Name of the task 
   */  
  AliTrackletAODTask(const char* name);
  /** 
   * Copy constructor 
   * 
   * @param other Object to copy from 
   */
  AliTrackletAODTask(const AliTrackletAODTask& other);
  /** 
   * Assignment operator 
   * 
   * @param other Object to assign from 
   * 
   * @return Reference to this object 
   */  
  AliTrackletAODTask& operator=(const AliTrackletAODTask& other);
  /** 
   * Print information to standard output 
   */
  void Print(Option_t*) const;
  /** 
   * Create our task and connect it 
   * 
   * @return Pointer to newly created task, or null
   */
  static AliTrackletAODTask* Create();
  /** 
   * @{ 
   * @name Task interface 
   */
  /**
   * Delegate worker initialization 
   * 
   */
  void UserCreateOutputObjects();
  /**
   * Event processing 
   */
  void UserExec(Option_t*); 
  /** 
   * Called at end of worker job. 
   * 
   */
  void FinishTaskOutput() { /*WorkerFinalize();*/ }
  /** 
   * Called at end of master job on merged results. 
   * 
   */
  void Terminate(Option_t*) {}
  /** 
   * Connect this task to the train
   * 
   * @return true on success 
   */
  Bool_t Connect();
  /* @} */
  /** 
   * @{ 
   * @name Set parameters of the reconstruction
   */
  /**
   * Set wether to scale @f$\Delta\theta@f$ by @f$\sin^2\theta@f$ 
   *
   * @param x If true, scale 
   */
  void SetScaleDTheta(Bool_t x=false) { fScaleDTheta = x; }
  /**
   * Set @f$\delta_{\phi}@f$
   *
   * @param x Shift of @f$\Delta\phi@f$ 
   */
  void SetDPhiShift(Double_t x=0.0045) { fDPhiShift = x; }
  /**
   * Set Maximum @f$ \Delta@f$ to consider 
   *
   * @param x Value 
   */
  void SetMaxDelta(Double_t x=25) { fMaxDelta = x; }
  /**
   * Set DThetaWindow
   *
   * @param x Value 
   */
  void SetDThetaWindow(Double_t x=0.025) { fDThetaWindow = x; }
  /**
   * Set @f$ d\phi@f$ window - used for reconstruction only
   *
   * @param x Value 
   */
  void SetDPhiWindow(Double_t x=0.06) { fDPhiWindow = x; }
  /**
   * Set PhiOverlapCut - used for reconstruction only
   *
   * @param x Value
   */
  void SetPhiOverlapCut(Double_t x=0.005) { fPhiOverlapCut = x; }
  /**
   * Set ZEtaOverlapCut - used for reconstruction only
   *
   * @param x Value
   */
  void SetZEtaOverlapCut(Double_t x=0.05) { fZEtaOverlapCut = x; }
  /* @} */
protected:
  /** 
   * @{ 
   * @name Worker initialization 
   */
  /** 
   * Initialize the worker 
   */
  Bool_t WorkerInit();
  /** 
   * Make sure CDB is initialized 
   * 
   * @return true on success
   */
  Bool_t InitCDB();
  /** 
   * Get the CDB reference run number 
   * 
   * @return A run from LHC10h
   */
  virtual Int_t GetCDBReferenceRun() const { return 137161; }
  /** 
   * Get the CDB reference URL 
   * 
   * @return A fixed string pointing to 2010 
   */
  virtual const char* GetCDBReferenceURL() const
  {
    return "alien://Folder=/alice/data/2010/OCDB";
  }
  /** 
   * Initialize our output branch
   * 
   * @return true on success 
   */
  Bool_t InitBranch();
  /** 
   * @{ 
   * @name Event processing 
   */
  /** 
   * Process a single event
   * 
   * 
   * @return 
   */
  virtual Bool_t ProcessEvent();
  /** 
   * Check that we have an initialized geometry if needed
   * 
   * @return true if all is good 
   */
  Bool_t HasGeometry();
  /** 
   * Check if we have the magnetic field 
   * 
   * @return true on success 
   */
  Bool_t HasField(AliVEvent* event);
  /** 
   * Find cluster (rec.point) tree
   * 
   * @return Pointer to tree or null
   */
  TTree* FindClusters();
  /** 
   * Get the event 
   * 
   * @return Pointer to event or null
   */
  AliVEvent* FindEvent();
  /** 
   * Find the interaction point location 
   * 
   * @param event         Event 
   * @param maxDispersion Max dispersion 
   * @param maxZError     Max error along Z
   * 
   * @return Pointer to vertex, or null in case of problems 
   */
  const AliVVertex* FindIP(AliVEvent* event,
			   Double_t maxDispersion=0.04,
			   Double_t maxZError=0.25);
  /** 
   * Reconstruct the tracklets from passed clusters 
   * 
   * @param clusters Clusters 
   * @param ip       Interaction point coordinates 
   */
  Bool_t Reconstruct(TTree* clusters, const AliVVertex* ip);
  /** 
   * Mark the event as one to store on AOD tree
   * 
   * @return true on success
   */
  Bool_t MarkForStore();
  /* @} */
  /**
   * @{ 
   * @name Tracklet creation 
   */
  /** 
   * The name of the tracklet class to use 
   * 
   * @return Class name as a string 
   */
  virtual const char* TrackletClassName() const { return "AliAODTracklet"; }
  /** 
   * Create a tracklet 
   *
   * @param normal If true, create for normal reconstruction 
   * 
   * @return Pointer to newly allocated tracklet 
   */
  virtual AliAODTracklet* MakeTracklet(Bool_t normal);
  /** 
   * Process tracklets create by the reconstruction 
   * 
   * @param normal Whether this this is normal reconstruction 
   * @param mult   The created tracklets 
   * 
   * @return true on success 
   */
  Bool_t ProcessTracklets(Bool_t normal, AliMultiplicity* mult);
  /** 
   * Process a single tracklet 
   * 
   * @param normal Whether this this is normal reconstruction 
   * @param mult   The created tracklets 
   * @param no     Tracklet number to investigate 
   * 
   * @return Newly allocated tracklet or null
   */
  virtual AliAODTracklet* ProcessTracklet(Bool_t           normal,
					  AliMultiplicity* mult,
					  Int_t            no);
  /* @} */
  
  /** Container of tracklets */
  TClonesArray* fTracklets;
  /* Output container */
  // Container* fContainer;
  /** Whether we should scale @f$ \Delta\theta@f$ by @f$\sin^{-2}(\theta)@f$ */ 
  Bool_t     fScaleDTheta;
  /** Maximum @f$\Delta@f$ to consider */
  Double_t   fMaxDelta;
  /** Shift in @f$\Delta\phi@f$ */
  Double_t   fDPhiShift;
  /** Window in @f$ \Delta\theta @f$ */
  Double_t   fDThetaWindow;  
  /** Window in @f$ \Delta\phi @f$ */
  Double_t   fDPhiWindow;
  /** Overlap cut in @f$\phi@f$ plane */
  Double_t   fPhiOverlapCut;
  /** Overlap cut in @f$ z,\eta@f$ plane */
  Double_t   fZEtaOverlapCut;
  /** Pointer to current reconstruction object */
  AliITSMultRecBg* fReco; //!
  
  ClassDef(AliTrackletAODTask,1); 
};
//____________________________________________________________________
AliTrackletAODTask::AliTrackletAODTask()
  : AliAnalysisTaskSE(),
    fTracklets(0),
    // fContainer(0),
    fScaleDTheta(true),
    fMaxDelta(25),
    fDPhiShift(0.0045),
    fDThetaWindow(0.025),
    fDPhiWindow(0.06),
    fPhiOverlapCut(0.005),
    fZEtaOverlapCut(0.05),
    fReco(0)
{}
//____________________________________________________________________
AliTrackletAODTask::AliTrackletAODTask(const char*)
  : AliAnalysisTaskSE("tracklets"),
    fTracklets(0),
    // fContainer(0),
    fScaleDTheta(true),
    fMaxDelta(25),
    fDPhiShift(0.0045),
    fDThetaWindow(0.025),
    fDPhiWindow(0.06),
    fPhiOverlapCut(0.005),
    fZEtaOverlapCut(0.05),
    fReco(0)
{
  // DefineOutput(1,TList::Class());
}
//____________________________________________________________________
AliTrackletAODTask::AliTrackletAODTask(const AliTrackletAODTask& other)
  : AliAnalysisTaskSE(other),
    fTracklets(0),
    // fContainer(0),
    fScaleDTheta(other.fScaleDTheta),
    fMaxDelta(other.fMaxDelta),
    fDPhiShift(other.fDPhiShift),
    fDThetaWindow(other.fDThetaWindow),
    fDPhiWindow(other.fDPhiWindow),
    fPhiOverlapCut(other.fPhiOverlapCut),
    fZEtaOverlapCut(other.fZEtaOverlapCut),
    fReco(0)
{}
//____________________________________________________________________
Bool_t AliTrackletAODTask::Connect()
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
  
  return true;
}
//____________________________________________________________________
void AliTrackletAODTask::Print(Option_t*) const
{
  Printf("%s: %s", ClassName(), GetName());
  Printf(" %22s: %d",   "Scale by sin^2(theta)",   fScaleDTheta);
  Printf(" %22s: %f",   "Delta phi shift",	   fDPhiShift);
  Printf(" %22s: %f",   "max Delta",	           fMaxDelta);
  Printf(" %22s: %f",   "Delta theta window",	   fDThetaWindow);
  Printf(" %22s: %f",   "Delta phi window",	   fDPhiWindow);
  Printf(" %22s: %f",   "phi overlap cut",	   fPhiOverlapCut);
  Printf(" %22s: %f",   "z-eta overlap cut",	   fZEtaOverlapCut);
}

//____________________________________________________________________
AliTrackletAODTask& AliTrackletAODTask::operator=(const AliTrackletAODTask&)
{
  return *this;
}
//____________________________________________________________________
void AliTrackletAODTask::UserCreateOutputObjects()
{
  // fContainer = new Container;
  // fContainer->SetOwner();
  // fContainer->SetName(GetName());
  // PostData(1, fContainer);
  
  // Create container of tracklets
  if (WorkerInit()) return;

  AliWarning("Failed to initialize task on worker, making zombie");
  SetZombie();
  
}

//____________________________________________________________________
Bool_t AliTrackletAODTask::InitCDB()
{
  Printf("Now initialising CDB");
  AliAnalysisManager* anaMgr = AliAnalysisManager::GetAnalysisManager();
  if (!anaMgr) {
    AliError("No manager defined!");
    return false;
  }
  // Check if we have the CDB connect task, and if so, do nothing as
  // we rely on that task set up things properly.
  const char*  cdbNames[] = {"CDBconnect", "cdb", 0 };
  const char** ptr        = cdbNames;
  while (*ptr) { 
    AliAnalysisTask* cdbConnect = anaMgr->GetTask(*ptr);
    if (cdbConnect && cdbConnect->IsA()->InheritsFrom("AliTaskCDBconnect")) {
      AliInfoF("CDB-connect task (%s: %s) present, do nothing",
	       cdbConnect->ClassName(), *ptr);
      return true;
    }
    ptr++;
  }
  // Otherwise, we need to do stuff ourselves
  Printf("Get the CDB manager");
  AliCDBManager* cdbMgr = AliCDBManager::Instance();
  if (!cdbMgr) {
    AliError("Failed to get instance of CDB manager");
    return false;
  }
  Int_t   refRun = GetCDBReferenceRun();
  TString refUrl = GetCDBReferenceURL();
  AliWarningF("Using reference CDB storage \"%s\" and run \"%d\"",
	      refUrl.Data(), refRun);
  // Set a very particular default storage. Perhaps we can do this
  // with specific storages instead!  Depends on whether the
  // reconstruction also uses CDB - probably does - in which case we
  // need specific storages for that too.
  cdbMgr->SetDefaultStorage(refUrl);
  // Now load our geometry - from a LHC10h run 
  Printf("Get Geometry entry");
  AliCDBEntry* cdbEnt = cdbMgr->Get("GRP/Geometry/Data", refRun);
  if (!cdbEnt) {
    AliErrorF("No geometry found from %d", refRun);
    return false;
  }
  // Initialize the geometry manager
  Printf("Set Geometry");
  AliGeomManager::SetGeometry(static_cast<TGeoManager*>(cdbEnt->GetObject()));
  // Now perform mis-alignment - again based on an LHC10h run!
  Printf("Misalign geometry");
  if (!AliGeomManager::ApplyAlignObjsToGeom("ITS",refRun,-1,-1)) {
    AliErrorF("Failed to misalign geometry from %d", refRun);
    return false;
  }
  return true;
}
//____________________________________________________________________
Bool_t AliTrackletAODTask::InitBranch()
{
  fTracklets              = new TClonesArray(TrackletClassName());
  // fTracklets->SetName(GetName());
  TObject*            obj = fTracklets;
  AliAnalysisManager* am  = AliAnalysisManager::GetAnalysisManager();
  AliAODHandler*      ah  = 
    dynamic_cast<AliAODHandler*>(am->GetOutputEventHandler());
  if (!ah) {
    AliWarning("No AOD output handler set in analysis manager");
    return false;
  }
  ah->AddBranch("TClonesArray", &obj);
  // AliAODEvent* aod = ah->GetAOD();
  // aod->Print();
  // ah->GetTree()->Print();
  return true;
}

//____________________________________________________________________
Bool_t AliTrackletAODTask::WorkerInit()
{
  if (!InitCDB())      return false;
  if (!InitBranch())   return false;
  
  // PostData(1, fContainer);

  return true;
}
//____________________________________________________________________
void AliTrackletAODTask::UserExec(Option_t*)
{
  // Clear container of tracklets
  if (!fTracklets) {
    AliWarning("Tracklets container not initialized, init must have failed!");
    return;
  }
  fTracklets->Clear();

  if (!ProcessEvent()) return;
  
  MarkForStore();
}

//____________________________________________________________________
Bool_t AliTrackletAODTask::ProcessEvent()
{
  AliVEvent*        event    = 0;
  if (!HasGeometry())         return false;
  if (!(event = FindEvent())) return false;
  if (!HasField(event))       return false;

  TTree*            clusters = 0;
  const AliVVertex* ip       = 0;
  if (!(clusters = FindClusters())) return false;
  if (!(ip       = FindIP(event)))  return false;

  return Reconstruct(clusters, ip);
}

//____________________________________________________________________
Bool_t AliTrackletAODTask::HasGeometry()
{
  if (!AliGeomManager::GetGeometry()) {
    AliError("No geometry loaded, needed for reconstruction");
    AliError("Add the AliTaskCDBconnect to the train");
    return false;
  }
  return true;
}
//____________________________________________________________________
Bool_t AliTrackletAODTask::HasField(AliVEvent* event)
{
  if (!TGeoGlobalMagField::Instance()->GetField() &&
      !event->InitMagneticField()) {
    AliWarning("Failed to initialize magnetic field");
    return false;
  }
  return true;
}

//____________________________________________________________________
AliVEvent* AliTrackletAODTask::FindEvent()
{
  AliVEvent* event = InputEvent();
  if (!event) AliWarning("No event");
  return event;
}
  
//____________________________________________________________________
TTree* AliTrackletAODTask::FindClusters()
{
  AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
  AliVEventHandler*   inh = mgr->GetInputEventHandler();
  if (!inh->IsA()->InheritsFrom(AliESDInputHandlerRP::Class())) {
    AliErrorF("Not the right kind of input handler: %s",
	      inh->ClassName());
    return 0;
  }
  AliESDInputHandlerRP* rph = static_cast<AliESDInputHandlerRP*>(inh);
  TTree* tree = rph->GetTreeR("ITS");
  if (!tree) {
    AliError("Tree of clusters (rec.points) not found");
    return 0;
  }
  return tree;
}

//____________________________________________________________________
const AliVVertex* AliTrackletAODTask::FindIP(AliVEvent* event,
					     Double_t   maxDispersion,
					     Double_t   maxZError)
{
  const AliVVertex* ip   = event->GetPrimaryVertex();
  if (ip->GetNContributors() <= 0) {
    AliWarning("Not enough contributors for IP");
    return 0;
  }   
  // If this is from the Z vertexer, do some checks 
  if (ip->IsFromVertexerZ()) {
    // Get covariance matrix
    Double_t covar[6];
    ip->GetCovarianceMatrix(covar);
    Double_t sigmaZ = TMath::Sqrt(covar[5]);
    if (sigmaZ >= maxZError) {
      AliWarningF("IPz resolution = %f >= %f", sigmaZ, maxZError);
      return 0;
    }
      
    // If this IP doesn not derive from AliVertex, don't check dispersion. 
    if (ip->IsA()->InheritsFrom(AliVertex::Class())) {
      const AliVertex* ipv = static_cast<const AliVertex*>(ip);
      // Dispersion is the parameter used by the vertexer for finding the IP. 
      if (ipv->GetDispersion() >= maxDispersion) {
	AliWarningF("IP dispersion = %f >= %f",
		    ipv->GetDispersion(), maxDispersion);
	return 0;
      }
    }
  }
    
#if 0
  // If we get here, we either have a full 3D vertex or track
  // vertex, and we should check if it is in range
  if (ip->GetZ() < fIPzAxis.GetXmin() || ip->GetZ() > fIPzAxis.GetXmax()) {
    AliWarningF("IPz = %fcm out of range [%f,%f]cm",
		ip->GetZ(), fIPzAxis.GetXmin(), fIPzAxis.GetXmax());
    return 0;
  }
#endif 
  // Good vertex, return it
  return ip;
}

//____________________________________________________________________
Bool_t AliTrackletAODTask::MarkForStore()
{
  
  AliAnalysisManager* am = AliAnalysisManager::GetAnalysisManager();
  AliAODHandler*      ah = 
    dynamic_cast<AliAODHandler*>(am->GetOutputEventHandler());
  if (!ah) return false;
  ah->SetFillAOD(kTRUE);
  return true;
}

//____________________________________________________________________
Bool_t AliTrackletAODTask::Reconstruct(TTree* clusters, const AliVVertex* ip)
{
  // Make an IP array for the reconstructor 
  Float_t ipv[3];
  ipv[0] = ip->GetX();
  ipv[1] = ip->GetY();
  ipv[2] = ip->GetZ();

  // Make the reconstructor
  fReco =  0;
  AliITSMultRecBg reco;
  reco.SetCreateClustersCopy        (true);
  reco.SetScaleDThetaBySin2T        (fScaleDTheta);
  reco.SetNStdDev                   (fMaxDelta);
  reco.SetPhiWindow                 (fDPhiWindow);
  reco.SetThetaWindow               (fDThetaWindow);
  reco.SetPhiShift                  (fDPhiShift);
  reco.SetRemoveClustersFromOverlaps(fPhiOverlapCut > 0 ||
				     fZEtaOverlapCut   > 0);
  reco.SetPhiOverlapCut             (fPhiOverlapCut);
  reco.SetZetaOverlapCut            (fZEtaOverlapCut);
  reco.SetHistOn                    (false);
  reco.SetRecType                   (AliITSMultRecBg::kData);

  // Run normal reconstruction 
  reco.Run(clusters, ipv);

  // Save pointer to reco object - for MC combinatorics
  fReco = &reco; 
  // And fill results into output branch 
  if (!ProcessTracklets(true, reco.GetMultiplicity())) return false;
  fReco = 0;
  
  // Run again, but in injection mode
  reco.SetRecType(AliITSMultRecBg::kBgInj);
  reco.Run(clusters, ipv);
  
  // And fill results into output branch 
  if(!ProcessTracklets(false, reco.GetMultiplicity())) return false;

  return true;
}

//____________________________________________________________________
Bool_t AliTrackletAODTask::ProcessTracklets(Bool_t normal,AliMultiplicity* mult)
{
  // Now loop over all tracklets.  Since this method is only called
  // once per "reconstruction" pass, we save a lot of computing time
  // since we loop over the tracklets of each reconstruction exactly
  // once. 
  Int_t nTracklets = mult->GetNumberOfTracklets();
  for (Int_t trackletNo = 0; trackletNo < nTracklets; trackletNo++) {
    if (!ProcessTracklet(normal,mult,trackletNo)) return false;
  }
  return true; 
}

//____________________________________________________________________
AliAODTracklet* AliTrackletAODTask::MakeTracklet(Bool_t)
{
  if (!fTracklets) return 0;
  Int_t n = fTracklets->GetEntries();
  return new((*fTracklets)[n]) AliAODTracklet;
}
//____________________________________________________________________
AliAODTracklet*
AliTrackletAODTask::ProcessTracklet(Bool_t            normal,
				    AliMultiplicity*  mult,
				    Int_t             no)
{
  Double_t theta   = mult->GetTheta     (no);
  Double_t phi     = mult->GetPhi       (no);
  Double_t dTheta  = mult->GetDeltaTheta(no);
  Double_t dPhi    = mult->GetDeltaPhi  (no);
  Double_t delta   = mult->CalcDist     (no);
  AliAODTracklet* tracklet = MakeTracklet(normal);
  tracklet->SetTheta (theta);
  tracklet->SetPhi   (phi);
  tracklet->SetDTheta(dTheta);
  tracklet->SetDPhi  (dPhi);
  tracklet->SetDelta (delta);
  if (!normal) tracklet->SetInjection();
  
  return tracklet;
}

/*********************************************************************
 *
 * Code for processing simulations 
 */
#ifndef __CINT__
#include <AliStack.h>
#include <AliMCEvent.h>
#include <AliGenEventHeader.h>
#include <TParticle.h>
#include <TParticlePDG.h>
#include <TDatabasePDG.h>
#else
class TParticlePDG; 
class TParticle;
class AliStack;
#endif 
//====================================================================
/**
 * Store tracklets on AOD 
 * 
 * @ingroup pwglf_forward_tracklets
 */
class AliTrackletAODMCTask : public AliTrackletAODTask
{
public:
  /** 
   * Default constructor - for ROOT I/O only 
   */
  AliTrackletAODMCTask() : AliTrackletAODTask() {}
  /** 
   * Named user constructor 
   * 
   * @param name Name of the task 
   */  
  AliTrackletAODMCTask(const char* name) : AliTrackletAODTask(name) {}
  /** 
   * Copy constructor 
   * 
   * @param other Object to copy from 
   */
  AliTrackletAODMCTask(const AliTrackletAODMCTask& other)
    : AliTrackletAODTask(other)
  {}
  /** 
   * Assignment operator 
   * 
   * @param other Object to assign from 
   * 
   * @return Reference to this object 
   */  
  AliTrackletAODMCTask& operator=(const AliTrackletAODMCTask& other)
  {
    return *this;
  }
private:
  /** 
   * @{ 
   * @name Event processing 
   */
  /* 
   * Initialize the worker 
   */
  // Bool_t WorkerInit();
  /** 
   * Process a single event
   * 
   * 
   * @return true on success
   */
  virtual Bool_t ProcessEvent();
  /** 
   * Process MC truth 

   * @return true on success
   */
  Bool_t ProcessGenerated();
  /* @} */
  /** 
   * @{ 
   * @name Worker initialization 
   */     
  /** 
   * Get the CDB reference URL 
   * 
   * @return A fixed string pointing to 2010 
   */
  virtual const char* GetCDBReferenceURL() const
  {
    return "alien://Folder=/alice/simulation/2008/v4-15-Release/Residual";
  }
  /* @} */
  /**
   * @{ 
   * @name Tracklet creation 
   */
  /** 
   * The name of the tracklet class to use 
   * 
   * @return Class name as a string 
   */
  virtual const char* TrackletClassName() const { return "AliAODMCTracklet"; }
  /** 
   * Create a tracklet 
   *
   * @param normal If true, create for normal reconstruction 
   * 
   * @return Pointer to newly allocated tracklet 
   */
  virtual AliAODTracklet* MakeTracklet(Bool_t normal);
  /** 
   * Process a single tracklet 
   * 
   * @param normal Whether this this is normal reconstruction 
   * @param mult   The created tracklets 
   * @param no     Tracklet number to investigate 
   * 
   * @return Newly allocated tracklet or null
   */
  virtual AliAODTracklet* ProcessTracklet(Bool_t           normal,
					  AliMultiplicity* mult,
					  Int_t            no);
  /* @} */
  /** 
   * Find first primary parent of particle identified by passed label
   * 
   * @param label Label particle to search for primary parent of 
   * 
   * @return Pointer to parent, or null if not found.  Note, if label
   * corresponds to a primary, then than particle is returned.
   */
  TParticle* FindPrimaryParent(Int_t label) const;
  /** 
   * Find parent of particle with label @a label.  If no parent is
   * found, or it goes beyond the stack, a negative value is returned.
   * 
   * @param label Label particle to get the parent for 
   * 
   * @return Label of parent to particle with label @a label. 
   */
  Int_t FindParent(Int_t label) const;
  /** 
   * Find list of parents corresponding to passed label.  
   *
   * We find the labels of all parent particles and store into cache
   * @a fill.  Elements are assigned starting from @a offset.  In this
   * way, we can invoke this member function multiple times and fill
   * into the same cache array.
   * 
   * @param label   Starting particle 
   * @param fill    Cache to fill into 
   * @param offset  Offset in @a fill to start assigning into 
   * 
   * @return Number of assigned elements
   */
  Int_t FindParents(Int_t label, TArrayI& fill,	Int_t offset) const;
  /** 
   * Find the the common parent of @a label and those given in @a
   * fill, if any.  If a common parent is found, return the index into
   * @a fill of that parent.
   * 
   * @param label Particles who's parents we're checking 
   * @param fill  List of known parents 
   * 
   * @return Index of common parent of @a label and those in @a fill 
   */
  Int_t CommonParent(Int_t label, const TArrayI& fill) const;
  ClassDef(AliTrackletAODMCTask,1); 
};

//____________________________________________________________________
Bool_t AliTrackletAODMCTask::ProcessEvent()
{
  if (!AliTrackletAODTask::ProcessEvent()) return false;

  // create generated "tracklets"
  if (!ProcessGenerated()) return false;

  return true;
}

//____________________________________________________________________
Bool_t AliTrackletAODMCTask::ProcessGenerated()
{
  AliStack* stack   = MCEvent()->Stack();
  Int_t     nTracks = stack->GetNtrack();
  for (Int_t trackNo = nTracks; trackNo--; ) {
    if (!stack->IsPhysicalPrimary(trackNo)) {
      // AliWarningF("Track # %6d not a primary", trackNo);
      // Not a primary, go on
      continue;
    }
    // Get the particle 
    TParticle* particle = stack->Particle(trackNo);
    if (!particle) {
      AliWarningF("No particle found for track # %d", trackNo);
      continue;
    }
    
    TParticlePDG* pdg = particle->GetPDG();
    // if (pdg->Charge() == 0) {
    // Uncharged - don't care
    // continue;
    // }
    
    // Get theta
    Double_t theta = particle->Theta();
    // Check for beam-like particle 
    if (theta < 1e-6 || TMath::Abs(theta-TMath::Pi()) < 1e-6) {
      AliWarningF("Track # %6d is beam-like (%f)", trackNo,
		  TMath::RadToDeg()*theta);    
      continue;
    }
    Double_t          phi      = particle->Phi();
    AliAODTracklet*   tracklet = MakeTracklet(false);
    AliAODMCTracklet* mc       = static_cast<AliAODMCTracklet*>(tracklet);
    mc->SetGenerated();
    mc->SetTheta(theta);
    mc->SetPhi(phi);
    mc->SetDTheta(0);
    mc->SetDPhi  (0);
    mc->SetDelta (0);
    mc->SetParentPdg(particle->GetPdgCode());
    mc->SetParentPt (particle->Pt());
    if (pdg->Charge() == 0) mc->SetNeutral();
  }
  return true;
}

//____________________________________________________________________
AliAODTracklet* AliTrackletAODMCTask::MakeTracklet(Bool_t)
{
  if (!fTracklets) return 0;
  Int_t n = fTracklets->GetEntries();
  return new((*fTracklets)[n]) AliAODMCTracklet;
}
//____________________________________________________________________
AliAODTracklet*
AliTrackletAODMCTask::ProcessTracklet(Bool_t            normal,
				      AliMultiplicity*  mult,
				      Int_t             no)
{
  AliAODTracklet* tracklet =
    AliTrackletAODTask::ProcessTracklet(normal,mult,no);
  if (!normal) return tracklet;

  AliAODMCTracklet* mc      = static_cast<AliAODMCTracklet*>(tracklet);
  Int_t             label0  = mult->GetLabel(no, 0);
  Int_t             label1  = mult->GetLabel(no, 1);
  TParticle*        parent0 = FindPrimaryParent(label0);
  if (parent0) {
    mc->SetParentPdg(parent0->GetPdgCode());
    mc->SetParentPt (parent0->Pt());
  }
  if (label0 != label1) {
    TParticle* parent1 = FindPrimaryParent(label1);
    if (parent1) { 
      mc->SetParentPdg(parent1->GetPdgCode(), true);
      mc->SetParentPt (parent1->Pt(),         true);
    }
    mc->SetCombinatorics();
    // Here, we could track back in the cluster labels to see if we
    // have the same ultimate mother of both clusters.  Since this
    // isn't really used, we do not do that
    if (fReco) {
      TArrayI parents(50); // At most 50 levels deep
      // Tracklet parameters
      Float_t* ftrack = fReco->GetTracklet(no); 
      // Cluster identifiers 
      Int_t    clus0  = Int_t(ftrack[AliITSMultReconstructor::kClID1]);
      Int_t    clus1  = Int_t(ftrack[AliITSMultReconstructor::kClID2]);
      // Cluster labelsx
      Float_t* fclus0 = (fReco->GetClusterOfLayer(0,clus0) +
			 AliITSMultReconstructor::kClMC0);
      Float_t* fclus1 = (fReco->GetClusterOfLayer(1,clus1) +
			 AliITSMultReconstructor::kClMC0);
      // Loop over three inner layer cluster labels
      Int_t offset = 0;
      for (Int_t i = 0; i < 3; i++)
	offset = FindParents(Int_t(fclus0[i]), parents, offset);
      // Loop over three outer layer cluster labels
      Bool_t distinct = true;
      for (Int_t i = 0; i < 3; i++) {
	Float_t flbl = fclus1[i];
	// Be careful not to get overflows on conversions
	if (flbl > nextafter(INT_MAX, 0) || flbl < nextafter(INT_MIN, 0))
	  continue; // Would overflow 
	if (CommonParent(Int_t(flbl), parents)) {
	  distinct = false;
	  // We break out as soon as we find a common parent. 
	  break;
	}
      } // loop over outer layer
      if (distinct) mc->SetDistinct();
    } // if (fReco)
  }
  else {
    if (!MCEvent()->Stack()->IsPhysicalPrimary(label0))
      mc->SetSecondary();
  }
  return tracklet;
}

//____________________________________________________________________
Int_t AliTrackletAODMCTask::FindParent(Int_t label) const
{
  AliStack*   stack     = MCEvent()->Stack();
  Int_t       nTracks   = stack->GetNtrack();
  TParticle*  particle  = stack->Particle(label);
  if (!particle) return -1;
  Int_t       ret       = particle->GetFirstMother();
  if (ret > nTracks) return -1;
  return ret;
}
    
//____________________________________________________________________
Int_t AliTrackletAODMCTask::FindParents(Int_t    label,
					TArrayI& fill,
					Int_t    offset) const
{
  // If offset is negative, then we don't really store the labels, but
  // only check against those already stored.
  Int_t       i         = offset;
  Int_t       lbl       = label;
  while (i < fill.GetSize()-1 && lbl >= 1) {
    fill[i]      = lbl;
    i++;
    lbl = FindParent(lbl);
  }
  // If we get here, and we're checking, that means we did not find a
  // common ancestor, so we return 0 (false).  If we're not checking,
  // then we return the number of elements assigned in the passed
  // cache.
  return i;
}

//____________________________________________________________________
Int_t AliTrackletAODMCTask::CommonParent(Int_t          label,
					 const TArrayI& fill) const
{
  Int_t       i         = 0;
  Int_t       lbl       = label;
  while (i < fill.GetSize()-1 && lbl >= 1) {
    // If we're checking, just see if we have the label in the list
    // already.  If we do, then return the index 
    for (Int_t j = 0; j < fill.GetSize(); j++) {
      if (fill[j] == lbl) return j;
    }
    i++;
    lbl = FindParent(lbl);
  }
  // If we get here, and we're checking, that means we did not find a
  // common ancestor, so we return 0 (false).  If we're not checking,
  // then we return the number of elements assigned in the passed
  // cache.
  return false;
}

  
//____________________________________________________________________
TParticle* AliTrackletAODMCTask::FindPrimaryParent(Int_t label) const
{
  AliStack*   stack     = MCEvent()->Stack();
  Int_t       nTracks   = stack->GetNtrack();
  TParticle*  particle  = stack->Particle(label);
  Int_t       trackNo   = label;
#if 1
  while (!stack->IsPhysicalPrimary(trackNo)) {
    trackNo  = particle->GetFirstMother();
    // If we have hit the top 
    if (trackNo < 0) return 0;
    // Partice first next iteration 
    particle = stack->Particle(trackNo);
  }
      
#else 
  Int_t       parentID  = particle->GetFirstMother();
  while (parentID >= 0 && parentID < nTracks) {
    particle  = stack->Particle(parentID);
    parentID  = particle->GetFirstMother();
  }
#endif 
  return particle;
}
//====================================================================
AliTrackletAODTask* AliTrackletAODTask::Create()
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
  Bool_t              mc  = mgr->GetMCtruthEventHandler() != 0;
  AliTrackletAODTask* ret = 0;
  if (mc)             ret = new AliTrackletAODMCTask("MidRapidityMC");
  else                ret = new AliTrackletAODTask("MidRapidity");
  if (ret)            ret->Connect();

  return ret;  
}


//====================================================================

//
// EOF
//


