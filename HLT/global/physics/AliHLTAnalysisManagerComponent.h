//-*- Mode: C++ -*-
// $Id: AliHLTAnalysisManagerComponent $

#ifndef ALIHLTANALYSISMANAGERCOMPONENT_H
#define ALIHLTANALYSISMANAGERCOMPONENT_H

/* This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/** @file    AliHLTAnalysisManagerComponent.cxx
    @author  David Rohr, Jens Wiechula, C. Zampolli, I. Vorobyev, M.Krzewicki
    @brief   Runs AliAnalysisTasks in a HLT component
*/

#include "AliHLTProcessor.h"
#include <map>
#include <string>
#include "AliHLTAsyncMemberProcessor.h"


class TH1F;
class TList;
class AliVEvent;
class AliVfriendEvent;

class AliHLTCTPData;
class AliHLTGlobalTriggerDecision;
class AliHLTAnalysisManager;
class AliHLTVEventInputHandler;

class AliHLTAnalysisManagerComponent : public AliHLTProcessor {
public:

  /*
   * ---------------------------------------------------------------------------------
   *                            Constructor / Destructor
   * ---------------------------------------------------------------------------------
   */

  typedef map<std::string,std::string> stringMap;

  /** constructor */
  AliHLTAnalysisManagerComponent();
  
  /** destructor */
  virtual ~AliHLTAnalysisManagerComponent();

  /*
   * ---------------------------------------------------------------------------------
   * Public functions to implement AliHLTComponent's interface.
   * These functions are required for the registration process
   * ---------------------------------------------------------------------------------
   */

  /** interface function, see @ref AliHLTComponent for description */
  const Char_t* GetComponentID();

  /** interface function, see @ref AliHLTComponent for description */
  void GetInputDataTypes( vector<AliHLTComponentDataType>& list);

  /** interface function, see @ref AliHLTComponent for description */
  AliHLTComponentDataType GetOutputDataType();

  /** interface function, see @ref AliHLTComponent for description */
  void GetOutputDataSize( ULong_t& constBase, Double_t& inputMultiplier );

  /** interface function, see @ref AliHLTComponent for description */
  void GetOCDBObjectDescription( TMap* const targetMap);

  /** interface function, see @ref AliHLTComponent for description */
  AliHLTComponent* Spawn();

  //new option parser
  static stringMap* TokenizeOptionString(const TString str);
  int ProcessOptionString(TString arguments);
  int ProcessOption(TString option, TString value);

 protected:

  /*
   * ---------------------------------------------------------------------------------
   * Protected functions to implement AliHLTComponent's interface.
   * These functions provide initialization as well as the actual processing
   * capabilities of the component. 
   * ---------------------------------------------------------------------------------
   */

  // AliHLTComponent interface functions

  /** interface function, see @ref AliHLTComponent for description */
  Int_t DoInit( Int_t /*argc*/, const Char_t** /*argv*/ );

  /** interface function, see @ref AliHLTComponent for description */
  Int_t DoDeinit();

  /** interface function, see @ref AliHLTComponent for description */
  Int_t DoEvent( const AliHLTComponentEventData& evtData, AliHLTComponentTriggerData& trigData);

  using AliHLTProcessor::DoEvent;


  /** interface function, see @ref AliHLTComponent for description */
  Int_t Reconfigure(const Char_t* cdbEntry, const Char_t* chainId);

  /** interface function, see @ref AliHLTComponent for description */
  Int_t ReadPreprocessorValues(const Char_t* modules);
 
  ///////////////////////////////////////////////////////////////////////////////////
  
private:

  /*
   * ---------------------------------------------------------------------------------
   * Private functions to implement AliHLTComponent's interface.
   * These functions provide initialization as well as the actual processing
   * capabilities of the component. 
   * ---------------------------------------------------------------------------------
   */

  /** copy constructor prohibited */
  AliHLTAnalysisManagerComponent(const AliHLTAnalysisManagerComponent&);

  /** assignment operator prohibited */
  AliHLTAnalysisManagerComponent& operator=(const AliHLTAnalysisManagerComponent&);


  /*
   * ---------------------------------------------------------------------------------
   *                              Helper
   * ---------------------------------------------------------------------------------
   */
  
  struct AnalysisManagerQueueData
  {
	  AnalysisManagerQueueData() : fEvent(NULL), fFriend(NULL), fRequestPush(false) {}
	  AliVEvent* fEvent;
	  AliVfriendEvent* fFriend;
	  bool fRequestPush;
  };
  
  Int_t ReadInput(AnalysisManagerQueueData* eventData);
  void CleanEventData(AnalysisManagerQueueData* eventData);
  void* AnalysisManagerInit(void*);
  void* AnalysisManagerExit(void*);
  void* AnalysisManagerDoEvent(void*);

  /*
   * ---------------------------------------------------------------------------------
   *                             Members - private
   * ---------------------------------------------------------------------------------
   */
  
  /** UID for merging */
  AliHLTUInt32_t fUID;                        // see above
  bool fQuickEndRun;				//After end of run we should not further process async tasks but stop asap
  bool fAnalysisInitialized;			//Is the AnalysisManager initialized correctly? If not, we don't process anything

  AliHLTAnalysisManager *fAnalysisManager;        // Manger

  AliHLTVEventInputHandler *fInputHandler;    // input handler

  //config stuff
  TString fAddTaskMacro;
  Bool_t fWriteAnalysisToFile;
  Bool_t fInitializeGeometry; //
  Bool_t fEnableDebug; //enable debug output - sysinfo,debug streamer, other files
  Bool_t fResetAfterPush; //reset the AnalysisManager after pushing for merging
  Int_t fPushEventModulo; //Push every n-th event
  Int_t fNEvents;		//Number of events processed
  Int_t fMinTracks;      //Min number of tracks to run AnalysisManager
  Int_t fNumEvents;      //Number of events used for analysis
  Int_t fQueueDepth;	//Depth of asynchronous Queue
  Int_t fAsyncProcess;	//Use an async process instead of an async thread
  Int_t fForceKillAsyncProcess; //Kill async process after n msec, -1 to disable
  Bool_t fPushRequestOngoing; //Have we already requested a push? If so, wait until push before request again
  AliHLTAsyncMemberProcessor<AliHLTAnalysisManagerComponent> fAsyncProcessor; //Processor for asynchronous processing

  ClassDef(AliHLTAnalysisManagerComponent, 2)
};
#endif
