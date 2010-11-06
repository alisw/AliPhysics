#ifndef ALIEVEHLTEVENTMANAGER_H
#define ALIEVEHLTEVENTMANAGER_H


#include "TEveEventManager.h"
#include "AliHLTLoggingVariadicFree.h"

class AliHLTHOMERBlockDesc;
class TList;

class TEveManager;
class TEveScene;
class TEveProjectionManager;
class TTimer;
class TEveViewer;

class AliHLTEvePhos;
class AliHLTEveEmcal;
class AliHLTEveTPC;
class AliHLTEveHLT;
class AliHLTEveITS;
class AliHLTEveISSD;
class AliHLTEveISDD;
class AliHLTEveISPD;
class AliHLTEveTRD;
class AliHLTEveAny;
class AliHLTEveMuon;
class AliHLTEveMultCorr;
class AliEveEventBuffer;
class AliESDEvent;

class AliEveHLTEventManager : public TEveElementList { 

public:

  ///Constructor
  AliEveHLTEventManager();
  
  virtual ~AliEveHLTEventManager();



  /**Set and get the global instance of the Eve manager */
  void SetEveManager(TEveManager * manager) {fEveManager = manager;}
  TEveManager * GetEveManager() const {return fEveManager;}

  /**Set and get the global instance of TGeoManager */
  void SetGeoManager(TGeoManager * manager) {fGeoManager = manager;}
  TGeoManager * GetGeoManager() const {return fGeoManager;}

  /** Set the projection scenes and their managers */
  void SetRPhiManager (TEveProjectionManager * mgr) {fRPhiManager = mgr;}
  void SetRPhiEventScene (TEveScene * scene ) {fRPhiEventScene = scene;}
  void SetRPhiViewer(TEveViewer * viewer ) {fRPhiViewer = viewer;}
  void SetRhoZManager(TEveProjectionManager * mgr) {fRhoZManager = mgr;}
  void SetRhoZEventScene(TEveScene * scene ) {fRhoZEventScene = scene;}
  void SetRhoZViewer(TEveViewer * viewer ) {fRhoZViewer = viewer;}

  /** Start and stop the automatic event loop */
  void StartLoop();
  void StopLoop();

  //* Show muon arm*//
  void SetShowMuon(Bool_t showmuon) { fShowMuon = showmuon; }

  /** Print the screens to a file **/
  void PrintScreens();


  virtual void NavigateBack() = 0;
  virtual void NavigateFwd() = 0;



  //Event buffer stuff
  virtual void ConnectEventBuffer();  
  virtual void StartBufferMonitor();
  virtual void NextEvent() = 0;

  virtual void SaveEveryThing();
  
 protected:

  /** copy constructor prohibited */
  AliEveHLTEventManager(const AliEveHLTEventManager&);

  /** assignment operator prohibited */
  AliEveHLTEventManager& operator=(const AliEveHLTEventManager&);

  /** Process the event data */
  Int_t ProcessEvent(TList * blockList);
  Int_t ProcessEvent(AliESDEvent * event);

  /** Set flag for event loop */
  void SetEventLoopStarted (Bool_t started) {fEventLoopStarted = started;}

  void DestroyDetectorElements();
  
  virtual AliEveEventBuffer * GetEventBuffer() {return NULL;}
  
  /** Process block */
  void ProcessBlock(AliHLTHOMERBlockDesc * block);  //Process block
  /** Reset the elements in the display */
  void ResetDisplay();  
  /** Update the display  */
  void UpdateDisplay(); 


  Int_t GetRunNumber() const { return fRunNumber; }
  ULong64_t GetEventId() const { return fEventId; }
  void SetRunNumber(Int_t rn) { fRunNumber = rn; }
  void SetEventId(ULong64_t id) { fEventId = id; }


  void CreatePhosElement();
  void CreateMultCorrElement();
  void CreateEmcalElement();
  void CreateTPCElement();
  void CreateITSElement();
  void CreateISPDElement();
  void CreateISDDElement();
  void CreateISSDElement();
  void CreateTRDElement();
  void CreateHLTElement();


  TGeoManager * fGeoManager;              //The global TGeoManager instance
  TEveManager * fEveManager;              //The global TEveManager instance
  TEveProjectionManager * fRPhiManager;   //The R - Phi projection scene manager
  TEveProjectionManager * fRhoZManager;   //The Rho- Z projection sene manager
  TEveScene * fRPhiEventScene;            //The R - Phi projection scene
  TEveScene * fRhoZEventScene;            //The Rho - Z projection sene
  TEveViewer * fRhoZViewer;
  TEveViewer * fRPhiViewer;
  


  TTimer * fTimer;                   //Timer for event loop
  //TTimer * fSourceListTimer;       //Timer for source list loop
 
  AliHLTEvePhos  * fPhosElement;     //Phos eve processor
  AliHLTEveEmcal * fEmcalElement;    //Emcal eve processor
  AliHLTEveTPC   * fTPCElement;      //TPC eve processor
  AliHLTEveHLT   * fHLTElement;      //HLT
  AliHLTEveITS   * fITSElement;      //ITS
  AliHLTEveISPD  * fISPDElement;     //ISPD
  AliHLTEveISSD  * fISSDElement;     //ISSD
  AliHLTEveISDD  * fISDDElement;     //ISDD
  AliHLTEveTRD   * fTRDElement;      //TRD
  AliHLTEveMuon  * fMuonElement;     //MUON
  AliHLTEveMultCorr  * fMultCorrElement;     //Mult corr
  AliHLTEveAny   * fAnyElement;      //Catch all



  Bool_t fEventLoopStarted;                    // Flag indicating whether the loop is running
  Bool_t fCenterProjectionsAtPrimaryVertex;    // Flag indicating whether to center the projection scenes at primary vertex (as opposed to 0, 0, 0)
  Bool_t fShowBarrel;                               // Display barrel detectors ?
  Bool_t fShowMuon;                                 // Display Muon arm ?

  Int_t fRunNumber;
  ULong64_t fEventId;

  ClassDef(AliEveHLTEventManager, 0);

};

#endif
