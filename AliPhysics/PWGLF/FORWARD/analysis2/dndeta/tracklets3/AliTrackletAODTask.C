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
#include <TParameter.h>
#include "AliTrackletAODUtils.C"
#ifndef __CINT__
#include "AliAODTracklet.C"
#include "AliTrackletWeights.C"
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
#include <TH1.h>
#include <TProfile2D.h>
#include <TUrl.h>
#include <climits>
#else
class AliTrackletAODUtils;
class AliTrackletBaseWeights;
class AliAODTracklet;
class AliVEvent;
class AliVVertex;
class AliMultiplicity;
class AliVMultiplicity;
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
class TH1;
class TProfile2D;
#endif


//====================================================================
/**
 * Store tracklets on AOD 
 * 
 * @ingroup pwglf_forward_tracklets
 */
class AliTrackletAODTask : public AliAnalysisTaskSE, public AliTrackletAODUtils
{
public:
  /** 
   * Status bins 
   */
  enum EStatus {
    /** Event was seen */
    kSeen,
    /** Event has clusters */
    kClusters,
    /** Clusters was filtered */
    kFilter,
    /** Event has IP */
    kIP,
    /** Reconstruction was run */
    kReconstruct,
    /** Injection was run */
    kInject,
    /** Tracklets was processed */
    kTracklets,
    /** Event was stored */
    kStored,
    /** Counter of bins */
    kNStatus
  };
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
  static AliTrackletAODTask* Create(const char* weights="");
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
  /** 
   * Whether to filter clusters corresponding to strange particles
   *
   * @param mode Mode of filtering 
   *
   * - 0 No filtering 
   * - 1 Random filtering 
   * - 2 Track filtering  
   */
  virtual void SetFilterMode(Int_t mode) { fFilterMode = mode; }
  /** 
   * Set weights to use for filtering. 
   * 
   * @param w Weights 
   */
  virtual void SetFilterWeights(AliTrackletBaseWeights* w) {}
  /* @} */
protected:
  /** 
   * @{ 
   * @name Worker initialization 
   */
  /** 
   * Initialize the worker 
   */
  virtual Bool_t WorkerInit();
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
   * Pre-process clusters.  This can remove clusters, etc. 
   *
   * @param t Tree of clusters
   */
  virtual TTree* FilterClusters(TTree* t) { return t; }
  /** 
   * Clean up possible copy of tree of clusters 
   * 
   * @param t 
   */
  virtual void CleanClusters(TTree*& t) {};
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
  Bool_t ProcessTracklets(Bool_t normal, AliVMultiplicity* mult);
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
  TClonesArray*       fTracklets;
  /* Output container */
  Container*          fContainer;
  /** Status histogram */
  TH1*                fStatus; //!
  /** Correlation of clusters and tracklets */
  TProfile2D*         fNClustersVsNTracklets; //!
  /** Whether we should scale @f$ \Delta\theta@f$ by @f$\sin^{-2}(\theta)@f$ */ 
  Bool_t              fScaleDTheta;
  /** Maximum @f$\Delta@f$ to consider */
  Double_t            fMaxDelta;
  /** Shift in @f$\Delta\phi@f$ */
  Double_t            fDPhiShift;
  /** Window in @f$ \Delta\theta @f$ */
  Double_t            fDThetaWindow;  
  /** Window in @f$ \Delta\phi @f$ */
  Double_t            fDPhiWindow;
  /** Overlap cut in @f$\phi@f$ plane */
  Double_t            fPhiOverlapCut;
  /** Overlap cut in @f$ z,\eta@f$ plane */
  Double_t            fZEtaOverlapCut;
  /** Pointer to current reconstruction object */
  AliITSMultRecBg*    fReco; //!
  /** Whether to remove clusters corresponding to strange primary parents */
  Int_t fFilterMode;
  /** Fraction of K^0_S clusters removed */
  TParameter<double>* fStrangeLoss; //! 
  
  ClassDef(AliTrackletAODTask,1); 
};
//____________________________________________________________________
AliTrackletAODTask::AliTrackletAODTask()
  : AliAnalysisTaskSE(),
    fTracklets(0),
    fContainer(0),
    fStatus(0),
    fNClustersVsNTracklets(0),
    fScaleDTheta(true),
    fMaxDelta(25),
    fDPhiShift(0.0045),
    fDThetaWindow(0.025),
    fDPhiWindow(0.06),
    fPhiOverlapCut(0.005),
    fZEtaOverlapCut(0.05),
    fReco(0),
    fFilterMode(0),
    fStrangeLoss(0)
{}
//____________________________________________________________________
AliTrackletAODTask::AliTrackletAODTask(const char*)
  : AliAnalysisTaskSE("tracklets"),
    fTracklets(0),
    fContainer(0),
    fStatus(0),
    fNClustersVsNTracklets(0),
    fScaleDTheta(true),
    fMaxDelta(25),
    fDPhiShift(0.0045),
    fDThetaWindow(0.025),
    fDPhiWindow(0.06),
    fPhiOverlapCut(0.005),
    fZEtaOverlapCut(0.05),
    fReco(0),
    fFilterMode(0),
    fStrangeLoss(0)
{
  DefineOutput(1,TList::Class());
}
//____________________________________________________________________
AliTrackletAODTask::AliTrackletAODTask(const AliTrackletAODTask& other)
  : AliAnalysisTaskSE(other),
    fTracklets(0),
    fContainer(0),
    fStatus(0),
    fNClustersVsNTracklets(0),
    fScaleDTheta(other.fScaleDTheta),
    fMaxDelta(other.fMaxDelta),
    fDPhiShift(other.fDPhiShift),
    fDThetaWindow(other.fDThetaWindow),
    fDPhiWindow(other.fDPhiWindow),
    fPhiOverlapCut(other.fPhiOverlapCut),
    fZEtaOverlapCut(other.fZEtaOverlapCut),
    fReco(0),
    fFilterMode(other.fFilterMode),
    fStrangeLoss(0)
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

  AliAnalysisDataContainer* sumCon =
    mgr->CreateContainer(Form("%sSums", GetName()), TList::Class(),
			 AliAnalysisManager::kOutputContainer,
			 AliAnalysisManager::GetCommonFileName());
  
  // Always connect input 
  mgr->ConnectInput(this, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(this, 1, sumCon);
  
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
  Printf(" %22s: %s",   "Filter strange",
	 fFilterMode==0 ? "no" :
	 fFilterMode==1 ? "random" : "track");
}

//____________________________________________________________________
AliTrackletAODTask& AliTrackletAODTask::operator=(const AliTrackletAODTask&)
{
  return *this;
}
//____________________________________________________________________
void AliTrackletAODTask::UserCreateOutputObjects()
{
  fContainer = new Container;
  fContainer->SetOwner();
  fContainer->SetName(GetName());
  PostData(1, fContainer);
  
  // Create container of tracklets
  if (WorkerInit()) return;

  AliWarning("Failed to initialize task on worker, making zombie");
  SetZombie();
  
}

//____________________________________________________________________
Bool_t AliTrackletAODTask::InitCDB()
{
  AliInfo("Now initialising CDB");
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
  AliInfo("Get the CDB manager");
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
  AliInfo("Get Geometry entry");
  AliCDBEntry* cdbEnt = cdbMgr->Get("GRP/Geometry/Data", refRun);
  if (!cdbEnt) {
    AliErrorF("No geometry found from %d", refRun);
    return false;
  }
  // Initialize the geometry manager
  AliInfo("Set Geometry");
  AliGeomManager::SetGeometry(static_cast<TGeoManager*>(cdbEnt->GetObject()));
  // Now perform mis-alignment - again based on an LHC10h run!
  AliInfo("Misalign geometry");
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

  if (fFilterMode > 0) {
    fStrangeLoss = new TParameter<double>("strLoss", 0);
    ah->AddBranch("TParameter<double>", &fStrangeLoss);
  }
  return true;
}

//____________________________________________________________________
Bool_t AliTrackletAODTask::WorkerInit()
{
  fStatus = new TH1D("status","Processing status",kNStatus,0,kNStatus);
  fStatus->SetDirectory(0);
  fStatus->SetYTitle("Events");
  fStatus->GetXaxis()->SetBinLabel(kSeen       +1,"Seen");
  fStatus->GetXaxis()->SetBinLabel(kClusters   +1,"w/Clusters");
  fStatus->GetXaxis()->SetBinLabel(kFilter     +1,"Filtered");
  fStatus->GetXaxis()->SetBinLabel(kIP         +1,"w/IP");
  fStatus->GetXaxis()->SetBinLabel(kReconstruct+1,"Reconstructed");
  fStatus->GetXaxis()->SetBinLabel(kInject     +1,"Injection");
  fStatus->GetXaxis()->SetBinLabel(kTracklets  +1,"Processed");
  fStatus->GetXaxis()->SetBinLabel(kStored     +1,"Stored");
  fContainer->Add(fStatus);
  
  fNClustersVsNTracklets =
    new TProfile2D("clustersVsTracklets",
		   "Correlation of clusters and tracklets",
		   1000,0,10000,
		   1000,0,10000);
  fNClustersVsNTracklets->SetXTitle("#it{N}_{cluster,0}");
  fNClustersVsNTracklets->SetYTitle("#it{N}_{cluster,1}");
  fNClustersVsNTracklets->SetZTitle("#it{N}_{tracklets}");
  fNClustersVsNTracklets->SetDirectory(0);
  fContainer->Add(fNClustersVsNTracklets);

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
  PostData(1, fContainer);
}

//____________________________________________________________________
Bool_t AliTrackletAODTask::ProcessEvent()
{
  AliVEvent* event = FindEvent();
  if (!event) return false;
  fStatus->Fill(kSeen);

  TTree*            clusters = FindClusters();
  const AliVVertex* ip       = FindIP(event);
  Bool_t            ret      = true;
  if (fDebug > 0) AliInfoF("Cluster tree: %p,  ip: %p", clusters, ip);
  if (!ip)                        return false;
  if (clusters) {
    if (!HasGeometry())           ret = false;
    if (fDebug > 1) AliInfoF("After geo check: %d", ret);
    if (ret && !HasField(event))  ret = false;
    if (fDebug > 1) AliInfoF("After field check: %d", ret);
    // If we have clusters, then try to reconstruct the event (again)
    if (ret)                      ret = Reconstruct(clusters, ip);
    if (fDebug > 1) AliInfoF("After reconstruction: %d", ret);
  }
  else
    // Otherwise, use the stored tracklets (max Delta = 1.5)
    ret = ProcessTracklets(true, event->GetMultiplicity());

  CleanClusters(clusters);
  if (fDebug > 1) AliInfoF("Return value: %d", ret);
  return ret;
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
  if (fDebug > 0 && !event) AliWarning("No event");
  return event;
}
  
//____________________________________________________________________
TTree* AliTrackletAODTask::FindClusters()
{
  AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
  AliVEventHandler*   inh = mgr->GetInputEventHandler();
  if (!inh->IsA()->InheritsFrom(AliESDInputHandlerRP::Class())) {
    if (fDebug > 0)
      AliWarningF("Clusters not available via input handler of class: %s",
		  inh->ClassName());
    return 0;
  }
  AliESDInputHandlerRP* rph = static_cast<AliESDInputHandlerRP*>(inh);
  TTree* tree = rph->GetTreeR("ITS");
  if (!tree) {
    AliWarning("Tree of clusters (rec.points) not found");
    return 0;
  }
  fStatus->Fill(kClusters);
  return FilterClusters(tree);
}

//____________________________________________________________________
const AliVVertex* AliTrackletAODTask::FindIP(AliVEvent* event,
					     Double_t   maxDispersion,
					     Double_t   maxZError)
{
  const AliVVertex* ip   = event->GetPrimaryVertex();
  if (!ip) return 0;
  if (ip->GetNContributors() <= 0) {
    if (fDebug > 0)
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
      if (fDebug) 
	AliWarningF("IPz resolution = %f >= %f", sigmaZ, maxZError);
      return 0;
    }
      
    // If this IP doesn not derive from AliVertex, don't check dispersion. 
    if (ip->IsA()->InheritsFrom(AliVertex::Class())) {
      const AliVertex* ipv = static_cast<const AliVertex*>(ip);
      // Dispersion is the parameter used by the vertexer for finding the IP. 
      if (ipv->GetDispersion() >= maxDispersion) {
	if (fDebug > 0)
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
  fStatus->Fill(kIP);
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
  if (fDebug > 0) AliInfo("Storing event");
  fStatus->Fill(kStored);
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
  if (fDebug > 0)
    AliInfoF("Reconstructing from cluster tree %p with (%4.2f,%4.2f,%7.4f",
	     clusters, ipv[0], ipv[1], ipv[2]);
  
  // --- Make the reconstructor --------------------------------------
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

  // --- Run normal reconstruction -----------------------------------
  reco.Run(clusters, ipv);
  fStatus->Fill(kReconstruct);
  fNClustersVsNTracklets->Fill(reco.GetNClustersLayer1(),
			       reco.GetNClustersLayer2(),
			       reco.GetMultiplicity()->GetNumberOfTracklets());

  // Save pointer to reco object - for MC combinatorics
  fReco = &reco; 
  // And fill results into output branch 
  if (!ProcessTracklets(true, reco.GetMultiplicity())) {
    AliWarning("Process tracklets (normal) failed"); 
    return false;
  }
  fReco = 0;
  
  // --- Run again, but in injection mode ----------------------------
  reco.SetRecType(AliITSMultRecBg::kBgInj);
  reco.Run(clusters, ipv);
  fStatus->Fill(kInject);
  
  // And fill results into output branch 
  if(!ProcessTracklets(false, reco.GetMultiplicity())) {
    AliWarning("Process tracklets (injection) failed"); 
    return false;
  }
  if (fDebug > 1) AliInfo("Returning true");
  return true;
}

//____________________________________________________________________
Bool_t AliTrackletAODTask::ProcessTracklets(Bool_t normal,
					    AliVMultiplicity* vmult)
{
  // Now loop over all tracklets.  Since this method is only called
  // once per "reconstruction" pass, we save a lot of computing time
  // since we loop over the tracklets of each reconstruction exactly
  // once.
  // if (!mult->IsA()->InheritsFrom(AliMultiplicity::Class())) return false;
  AliMultiplicity* mult = static_cast<AliMultiplicity*>(vmult);
  Int_t nTracklets = mult->GetNumberOfTracklets();
  for (Int_t trackletNo = 0; trackletNo < nTracklets; trackletNo++) {
    if (!ProcessTracklet(normal,mult,trackletNo)) return false;
  }
  if (normal) {
    Int_t n0 = mult->GetNumberOfITSClusters(0);
    Int_t n1 = mult->GetNumberOfITSClusters(1); 
    if (fDebug > 0)
      Printf("%d x %d clusters -> %d", n0, n1, nTracklets);
    if (n0 > 0 || n1 >> 0) 
      fNClustersVsNTracklets->Fill(n0, n1, nTracklets);
    fStatus->Fill(kTracklets);
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
#include <TArrayF.h>
#include <TParticle.h>
#include <TParticlePDG.h>
#include <TDatabasePDG.h>
#include <AliITSgeomTGeo.h>
#include <TRandom.h>
#include <TROOT.h>
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
  AliTrackletAODMCTask()
    : AliTrackletAODTask(),
      fFilterWeights(0),
      fSeenTrackPDGs(0),
      fUsedTrackPDGs(0),
      fSeenClusterPDGs(0),
      fUsedClusterPDGs(0)
  {}
  /** 
   * Named user constructor 
   * 
   * @param name Name of the task 
   */  
  AliTrackletAODMCTask(const char* name)
    : AliTrackletAODTask(name),
      fFilterWeights(0),
      fSeenTrackPDGs(0),
      fUsedTrackPDGs(0),
      fSeenClusterPDGs(0),
      fUsedClusterPDGs(0)
  {}
  /** 
   * Copy constructor 
   * 
   * @param other Object to copy from 
   */
  AliTrackletAODMCTask(const AliTrackletAODMCTask& other)
    : AliTrackletAODTask(other),
      fFilterWeights(0),
      fSeenTrackPDGs(0),
      fUsedTrackPDGs(0),
      fSeenClusterPDGs(0),
      fUsedClusterPDGs(0)
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
  /** 
   * Set weights to use for filtering. 
   * 
   * @param w Weights 
   */
  virtual void SetFilterWeights(AliTrackletBaseWeights* w) { fFilterWeights=w; }
protected:
  /** 
   * @{ 
   * @name Event processing 
   */
  /* 
   * Initialize the worker 
   */
  Bool_t WorkerInit();
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
   * @name Filtering of tracks and clusters 
   */
  /** 
   * Pre-process clusters.  This can remove clusters, etc. 
   *
   * @param t Tree of clusters
   */
  virtual TTree* FilterClusters(TTree* t);
  /** 
   * Pre-process clusters.  This can remove clusters independently of
   * each other.  That is, for each cluster we look up the parent
   * primary particle(s), and find the corresponding weight @f$ w@f$.
   * We then throw a dice and if number is larger than @f$ 1-1/w@f$,
   * we remove the cluster.
   *
   * @param t    Tree of clusters
   * @param copy Tree to fill with clusters 
   * @param cent Centrality (mostly 0)
   * @param ipz  Collision point z-coordinate
   */
  virtual void FilterClustersRandom(TTree* t, TTree* copy,
				    Double_t cent, Double_t ipz);
  /** 
   * Pre-process clusters.  For all primary particles in the stack, we
   * compute a weight @f$ w@f$.  We then throw a dice, and if the
   * number is larger than @f$ 1-1/w@f$, we mark that particle for
   * removal.  We then loop over all clusters in the SPD, and find the
   * primary particle(s) that each cluster corresponds to.  If such a
   * primary particle is marked for removal, we remove the clusters.
   *
   * @param t    Tree of clusters
   * @param copy Tree to fill with clusters 
   * @param cent Centrality (mostly 0)
   * @param ipz  Collision point z-coordinate
   */
  virtual void FilterClustersTrack(TTree* t, TTree* copy,
				    Double_t cent, Double_t ipz);
  /** 
   * Clean up possible copy of tree of clusters 
   * 
   * @param t 
   */
  virtual void CleanClusters(TTree*& t);
  /** 
   * Find the particle weight 
   * 
   * @param particle The particle 
   * 
   * @return The weight
   */
  Double_t LookupWeight(TParticle* particle,
			Double_t   cent=0,
			Double_t   ipz=0) const;
  /** 
   * Pick random number to see if we should keep a track 
   * 
   * @param weight The weight 
   * 
   * @return true if we're to keep it 
   */
  Bool_t KeepIt(Double_t weight) const;
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
   * @{ 
   * @name Investigating primary parents 
   */
  /** 
   * Find identifier of the first primary parent of particle
   * identified by passed label
   * 
   * @param label Label particle to search for primary parent of 
   * 
   * @return Identifier of parent, or negative if not found.  Note, if label
   * corresponds to a primary, then that ID is returned.
   */
  Int_t FindPrimaryParentID(Int_t label) const;
  /** 
   * Find first primary parent of particle identified by passed label
   * 
   * @param label Label particle to search for primary parent of 
   * 
   * @return Pointer to parent, or null if not found.  Note, if label
   * corresponds to a primary, then that particle is returned.
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
  /* @} */
  /** Filter weights */
  AliTrackletBaseWeights* fFilterWeights;
  /** Number of tracks per PDG */
  TH1* fSeenTrackPDGs; //!
  /** Tracks that survive the filtering */
  TH1* fUsedTrackPDGs; //!
  /** Number of clusters per layer per PDG */
  TH2* fSeenClusterPDGs; //!
  /** Number of clusters that survive filter per layer per PDG */
  TH2* fUsedClusterPDGs; //! 
  
  ClassDef(AliTrackletAODMCTask,1); 
};

#define ADDP(P,M,L,Q) \
  db->AddParticle("p" #P, "p" #P,M,false,L,Q,"Baryon",P); \
  db->AddAntiParticle("p" #P "_bar", -P)

//____________________________________________________________________
Bool_t AliTrackletAODMCTask::WorkerInit()
{
  if (!AliTrackletAODTask::WorkerInit()) return false;

  const TAxis& pdgAxis = PdgAxis();
  TAxis layerAxis(2,0,2);
  layerAxis.SetBinLabel(1,"Layer 0");
  layerAxis.SetBinLabel(2,"Layer 1");
  FixAxis(layerAxis,"Layer");
  fSeenTrackPDGs = Make1D(fContainer, "seenTrackPdg",
			 "Seen track PDGs", kGreen+1,20, pdgAxis);  
  fUsedTrackPDGs = Make1D(fContainer, "usedTrackPdg",
			  "Used track PDGs", kBlue+1,24, pdgAxis);
  fSeenClusterPDGs = Make2D(fContainer, "seenClusterPdg",
			   "Seen cluster PDGs", kGreen+1,21,pdgAxis,layerAxis);
  fUsedClusterPDGs = Make2D(fContainer, "usedClusterPdg",
			   "Used cluster PDGs", kBlue+1,25,pdgAxis,layerAxis);
  
  
  // Register additional particles from EPOS-LHC 
  TDatabasePDG* db  = TDatabasePDG::Instance();
  db->GetParticle(2212);
  db->AddParticle("deut","deut",1.876560,false,0.000000,3,"Baryon",1000010020);
  db->AddAntiParticle("deut_bar",-1000010020);
  db->AddParticle("trit","trit",2.816700,false,0.000000,3,"Baryon",1000010030);
  db->AddAntiParticle("trit_bar",-1000010030);
  db->AddParticle("alph","alph",3.755000,false,0.000000,6,"Baryon",1000020040);
  db->AddAntiParticle("alph_bar",-1000020040);
  ADDP(32224,  1.524000,0.145000,6);
  ADDP(12224,  1.816000,0.350000,6);
  ADDP(12222,  2.108000,0.350000,6);
  ADDP(12212,  1.525720,0.145000,3);
  ADDP(2124,   1.819440,0.150000,3);
  ADDP(32214,  2.113160,0.150000,3);
  ADDP(32212,  2.406880,0.145000,3);
  ADDP(12214,  2.700600,0.300000,3);
  ADDP(22124,  2.994320,0.150000,3);
  ADDP(12122,  3.288040,0.300000,3);
  ADDP(13222,  1.575200,0.080000,3);
  ADDP(23222,  1.768100,0.100000,3);
  ADDP(13226,  1.961000,0.170000,3);
  ADDP(12112,  1.524430,0.350000,0);
  ADDP(1214,   1.816860,0.150000,0);
  ADDP(32114,  2.109290,0.150000,0);
  ADDP(32112,  2.401720,0.145000,0);
  ADDP(12114,  2.694150,0.300000,0);
  ADDP(21214,  2.986579,0.150000,0);
  ADDP(11212,  3.279010,0.300000,0);
  ADDP(13122,  1.761000,0.040000,0);
  ADDP(3124,   1.950500,0.016000,0);
  ADDP(23122,  2.140000,0.090000,0);
  ADDP(13212,  2.329500,0.080000,0);
  ADDP(23212,  2.519000,0.100000,0);
  ADDP(43122,  2.708500,0.145000,0);
  ADDP(13216,  2.898000,0.170000,0);
  ADDP(31114,  1.524000,0.145000,-3);
  ADDP(11114,  1.816000,0.350000,-3);
  ADDP(11112,  2.108000,0.350000,-3);
  ADDP(13112,  1.577600,0.080000,-3);
  ADDP(23112,  1.767700,0.100000,-3);
  ADDP(13116,  1.957800,0.170000,-3);
  ADDP(9900110,0.000000,0.000000,0);
  ADDP(9900210,0.000000,0.000000,0);
  ADDP(9900220,0.000000,0.000000,0);
  ADDP(9900330,0.000000,0.000000,0);
  ADDP(9900440,0.000000,0.000000,0);
  ADDP(9902210,0.000000,0.000000,0);
  ADDP(9902110,0.000000,0.000000,0);
  ADDP(88,     0.000000,0.000000,0);
  ADDP(90,     0.000000,0.000000,0);
  ADDP(990,    0.000000,0.000000,0);
  ADDP(99999,  0.000000,0.000000,0);

  return true;
}  

//____________________________________________________________________
TTree* AliTrackletAODMCTask::FilterClusters(TTree* t)
{
  if (fStrangeLoss) fStrangeLoss->SetVal(0);
  if (!t || fFilterMode <= 0) {
    if (fDebug > 0) AliInfo("Returning original cluster tree");
    return t;
  }

  TDirectory*    savDir = gDirectory;
  gDirectory            = gROOT;
  TTree*         copy   = new TTree("TreeR", "TreeR");
  copy->SetAutoFlush(0);// Keep in memory 

  Double_t cent = 0;
  Double_t ipz  = 0;
  if (MCEvent()->GenEventHeader()) {
    TArrayF v(3);
    MCEvent()->GenEventHeader()->PrimaryVertex(v);
    ipz = v[2];
  }
  
  if (fFilterMode == 1) FilterClustersRandom(t, copy, cent, ipz);
  else                  FilterClustersTrack(t, copy, cent, ipz);
    
  savDir->cd();
  if (fDebug > 0)
    AliInfoF("Returning reduced cluster tree %p (original %p)", copy, t);
  fStatus->Fill(kFilter);
  return copy;
}

  
//____________________________________________________________________
void AliTrackletAODMCTask::FilterClustersRandom(TTree*   t,
						TTree*   copy,
						Double_t cent,
						Double_t ipz)
{
  TClonesArray*  in     = new TClonesArray("AliITSRecPoint");
  TClonesArray*  out    = new TClonesArray("AliITSRecPoint");
  Int_t          outN   = 0;
  copy->Branch("ITSRecPoints", &out);
  t->SetBranchAddress("ITSRecPoints", &in);

  // Printf("Filtering clusters from strange particles");

  Int_t min1   = AliITSgeomTGeo::GetModuleIndex(1,1,1);
  Int_t max1   = AliITSgeomTGeo::GetModuleIndex(2,1,1);
  Int_t min2   = AliITSgeomTGeo::GetModuleIndex(2,1,1);
  Int_t max2   = AliITSgeomTGeo::GetModuleIndex(3,1,1);
  Int_t nTotal = 0;
  Int_t nKept  = 0;
  
  // Loop over the modules of the SPD 
  for (Int_t i = 0; i < max2; i++) {
    in->Clear();
    out->Clear();
    outN = 0;

    // Read in module data 
    t->GetEntry(i);

    // Loop over all clusters in the module 
    Int_t inN = in->GetEntries();
    for (Int_t j = 0; j < inN; j++) {
      AliITSRecPoint* inCl = static_cast<AliITSRecPoint*>(in->At(j));
      if (!inCl) continue;

      // Loop over labels in the module
      Double_t weight = 1;
      for (Int_t k = 0; k < 3; k++) {
	Int_t label = inCl->GetLabel(k);
	if (label <= 0) continue;

	// Check primary parent particle type 
	TParticle* parent = FindPrimaryParent(label);
	if (!parent) continue;

	// weight *= GetWeight(parent->GetPdgCode());
	weight *= LookupWeight(parent,cent,ipz);
	if (weight > 1)
	  Printf("Cluster from %+5d: %7.5f", parent->GetPdgCode(), weight);
	// Printf("Parent %d is a K^0_S: %d", k+1, parent->GetPdgCode());
	// Randomly remove clusters from the right kind of primary
	// parent particles
      }
      if (weight > 1) {
	nTotal++;
	if (!KeepIt(weight)) continue;
	nKept++;
      }
      new ((*out)[outN++]) AliITSRecPoint(*inCl);
    }
    // Printf("Kept %d out of %d clusters", outN, inN);
    copy->Fill();
  }
  if (nTotal > 0)  {
    Double_t loss = 1-Float_t(nKept)/nTotal;
    Printf("Kept %d out of %d clusters from strange primaries (%4.1f%%)",
	   nKept, nTotal, 100.*loss);
    if (fStrangeLoss) fStrangeLoss->SetVal(loss);
  }
  delete in;
  delete out;
}

//____________________________________________________________________
void AliTrackletAODMCTask::FilterClustersTrack(TTree*   t,
					       TTree*   copy,
					       Double_t cent,
					       Double_t ipz)
{
  TClonesArray*  in     = new TClonesArray("AliITSRecPoint");
  TClonesArray*  out    = new TClonesArray("AliITSRecPoint");
  Int_t          outN   = 0;
  copy->Branch("ITSRecPoints", &out);
  t->SetBranchAddress("ITSRecPoints", &in);

  // Loop over the stack
  Int_t     nTotal  = 0;
  Int_t     nKept   = 0;
  AliStack* stack   = MCEvent()->Stack();
  Int_t     nTracks = stack->GetNtrack();
  TBits     kept(nTracks);
  TBits     seen(nTracks);
  kept.ResetAllBits(true);
  seen.ResetAllBits(false);

  for (Int_t trackNo = nTracks; trackNo--; ) {
    // Get the primary parent identifier 
    Int_t      parent = FindPrimaryParentID(trackNo);
    if (parent < 0) continue;
    // Check if we've seen this parent already
    if (seen.TestBitNumber(parent)) continue;
    seen.SetBitNumber(parent,true);
    
    // Get the primary parent track and weight
    TParticle* par    = stack->Particle(parent);
    // Double_t   weight = GetWeight(par->GetPdgCode());
    Double_t   weight = LookupWeight(par,cent,ipz);
    Bool_t     keep   = true;
    if (weight > 1) {
      nTotal++;
      keep = KeepIt(weight);
      if (keep) nKept++;
    }
    Int_t bin = PdgBin(par->GetPdgCode());
    fSeenTrackPDGs->Fill(bin);
    if (keep) fUsedTrackPDGs->Fill(bin);
    if (fDebug > 1) 
      Printf("Primary parent %6d from a %6d %6s (%f)",
	     parent, par->GetPdgCode(), keep ? "kept" : "marked", weight);
    kept.SetBitNumber(parent, keep);
    par->SetWeight(keep ? 1 : weight);
  }
  // At this point, we have investigated all relevant primary
  // particles and thrown a dice to see if we should remove clusters
  // that correspond to these primary particles.
  if (fDebug > 0 && nTotal > 0)  {
    Double_t loss = 1-Float_t(nKept)/nTotal;
    Printf("Kept %d out of %d strange primaries (%4.1f%% loss)",
	   nKept, nTotal, 100.*loss);
    if (fStrangeLoss) fStrangeLoss->SetVal(loss);
  }

  // Next thing is to actually remove the clusters
  // Printf("Filtering clusters from strange particles");
  Int_t min1   = AliITSgeomTGeo::GetModuleIndex(1,1,1);
  Int_t max1   = AliITSgeomTGeo::GetModuleIndex(2,1,1);
  Int_t min2   = AliITSgeomTGeo::GetModuleIndex(2,1,1);
  Int_t max2   = AliITSgeomTGeo::GetModuleIndex(3,1,1);
  Int_t inT    = 0;
  Int_t outT   = 0;
  // Loop over the modules of the SPD 
  for (Int_t i = 0; i < max2; i++) {
    in->Clear();
    out->Clear();
    outN = 0;

    // Read in module data 
    t->GetEntry(i);

    // Loop over all clusters in the module 
    Int_t inN = in->GetEntries();
    inT       += inN;
    for (Int_t j = 0; j < inN; j++) {
      AliITSRecPoint* inCl = static_cast<AliITSRecPoint*>(in->At(j));
      if (!inCl) continue;
      
      // Loop over labels of the cluster
      Bool_t toRemove = false;
      Int_t  pdg      = 0;
      for (Int_t k = 0; k < 3; k++) {
	Int_t label = inCl->GetLabel(k);
	if (label <= 0) continue;

	// Check primary parent particle type 
	Int_t parent = FindPrimaryParentID(label);
	if (parent < 0) continue;

	if (pdg == 0) pdg = stack->Particle(parent)->GetPdgCode();
	if (!kept.TestBitNumber(parent)) toRemove = true;
	if (fDebug > 3) {
	  Printf("Cluster %6d from parent %6d (%7.5f) of type %6d %s",
		 j, parent,
		 stack->Particle(parent)->GetPdgCode(),
		 stack->Particle(parent)->GetWeight(),
		 toRemove ? "removed" : "kept");
	}
      }
      // We've now looked at all primaries for this cluster, and if
      // any of these primaries was flagged for removal, then we
      // should remove the cluster.
      if (fDebug > 3)
	Printf("Cluster %2d/%3d from %4d is %s",
	       i, j, pdg, toRemove ? "removed" : "kept");
      fSeenClusterPDGs->Fill(PdgBin(pdg), i < max1);
      if (toRemove) continue;
      fUsedClusterPDGs->Fill(PdgBin(pdg), i < max1);
      new ((*out)[outN++]) AliITSRecPoint(*inCl);
    }
    // Printf("Kept %d out of %d clusters", outN, inN);
    outT += outN;
    copy->Fill();
  }
  if (fDebug > 0) 
    Printf("Wrote out %6d out of %6d clusters (%4.1f%% loss)",
	   outT, inT, (inT > 0 ? 100*float(inT-outT)/inT : 100));
  delete in;
  delete out;
}

//____________________________________________________________________
Double_t AliTrackletAODMCTask::LookupWeight(TParticle* particle,
					    Double_t   cent,
					    Double_t   ipz) const
{
#if 1
  if (!fFilterWeights) return 1;
  return fFilterWeights->LookupWeight(particle,cent,ipz);
#else
  const Double_t k0s    = 1.52233299626516083e+00; // 310 - K^0_S weight
  const Double_t kpm    = (1.43744204476109627e+00*
			   9.82150320171071400e-01); // 321  - K^{+/-}
  const Double_t lam    = 2.75002089647900005e+00;   // 3122 - lambda
  const Double_t sig    = 2.75002089647899961e+00;   // 3212 - sigma
  const Double_t xi     = 3.24109605656453548e+00;   // 3322 - Xi

  switch (TMath::Abs(particle->GetPdgCode())) {
  case 310:  return k0s;
  case 321:  return kpm;
  case 3122: return lam;
    // case 3212: return sig; // Old, wrong code 
  case 3112:
  case 3222: return sig;
  case 3312: return xi;
    // case 3322: return xi; // Old, wrong code 
  }
  return 1;
  }
#endif
}

//____________________________________________________________________
Bool_t AliTrackletAODMCTask::KeepIt(Double_t weight) const
{
  if (weight <= 1) return true;
  Double_t chance = 1 - 1 / weight; // chance is 1 minus inverse 	
  return (gRandom->Uniform() >= chance);
}

//____________________________________________________________________
void AliTrackletAODMCTask::CleanClusters(TTree*& t)
{
  if (!t || fFilterMode <= 0) return;
 
  delete t;
  t = 0;
}

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
  if (fDebug > 0) AliInfo("Processing generated particles");
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
    if (!pdg) {
      AliWarningF("Unknown PDG code: %d", particle->GetPdgCode());
      // continue; // Unknown particle
    }
    // if (pdg->Charge() == 0) {
    // Uncharged - don't care
    // continue;
    // }
    
    // Get theta
    Double_t theta = particle->Theta();
    // Check for beam-like particle 
    if (theta < 1e-6 || TMath::Abs(theta-TMath::Pi()) < 1e-6) {
      if (fDebug > 0)
	AliWarningF("Track # %6d is beam-like (%f)", trackNo,
		    TMath::RadToDeg()*theta);    
      continue;
    }
    Double_t eta = -TMath::Log(TMath::ATan(theta/2));
    // If the pseudorapidity is way beyond the SPD acceptance, do not
    // write down this generated "tracklet" - to save space.
    if (TMath::Abs(eta) > 3) continue;
    
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
    if (particle->GetWeight() > 1) mc->SetSuppressed();
    if (pdg && pdg->Charge() == 0) mc->SetNeutral();
  }
  if (fDebug > 2) AliInfo("Returning true from generated");
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
  Int_t parent = FindPrimaryParentID(label);
  if (parent < 0) return 0;
  return MCEvent()->Stack()->Particle(parent);
}
//____________________________________________________________________
Int_t AliTrackletAODMCTask::FindPrimaryParentID(Int_t label) const
{
  AliStack*   stack     = MCEvent()->Stack();
  Int_t       nTracks   = stack->GetNtrack();
  Int_t       trackNo   = label;
  if (trackNo > nTracks || trackNo <= 0) return 0;
  TParticle*  particle  = stack->Particle(label);
  while (!stack->IsPhysicalPrimary(trackNo)) {
    trackNo  = particle->GetFirstMother();
    // If we have hit the top 
    if (trackNo < 0) return 0;
    // Partice first next iteration 
    particle = stack->Particle(trackNo);
  }
  return trackNo;
}
//====================================================================
AliTrackletAODTask* AliTrackletAODTask::Create(const char* weights)
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

  if (weights && weights[0] != '\0') {
    TUrl   wurl(weights);
    TFile* wfile = TFile::Open(wurl.GetFile());
    if (!wfile) {
      ::Warning("Create", "Failed to open weights file: %s",
		wurl.GetUrl());
      return 0;
    }
    TString wnam(wurl.GetAnchor());
    if (wnam.IsNull()) wnam = "weights";

    TObject* wobj = wfile->Get(wnam);
    if (!wobj) {
      ::Warning("Create", "Failed to get weights %s from file %s",
		wnam.Data(), wfile->GetName());
      return 0;
    }
    if (!wobj->IsA()->InheritsFrom(AliTrackletBaseWeights::Class())) {
      ::Warning("Create", "Object %s from file %s not an "
		"AliTrackletBaseWeights but a %s",
		wnam.Data(), wfile->GetName(), wobj->ClassName());
      return 0;
    }
    ret->SetFilterWeights(static_cast<AliTrackletBaseWeights*>(wobj));
  }
  return ret;  
}


//====================================================================

//
// EOF
//


