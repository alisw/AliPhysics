//-*- Mode: C++ -*-
// $Id: AliHLTGlobalCompareFlatComponent $

#ifndef ALIHLTGLOBALCOMPAREFLATCOMPONENT_H
#define ALIHLTGLOBALCOMPAREFLATCOMPONENT_H

/* This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/** @file    AliHLTGlobalCompareFlatComponent.h
    @author  Steffen Weber <s.weber@gsi.de>
    @brief   Compare flat events from different inputs
*/

// see below for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include "AliHLTProcessor.h"

class THnSparse;

class AliESDVZERO;
class AliESDtrackCuts;
class AliHLTCTPData;
class AliHLTMultiplicityCorrelations;
class AliHLTGlobalTriggerDecision;
class AliHLTVEventInputHandler;


class AliHLTGlobalCompareFlatComponent : public AliHLTProcessor {
public:

  /*
   * ---------------------------------------------------------------------------------
   *                            Constructor / Destructor
   * ---------------------------------------------------------------------------------
   */

  /** constructor */
  AliHLTGlobalCompareFlatComponent();
  
  /** destructor */
  virtual ~AliHLTGlobalCompareFlatComponent();

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
 // void GetOCDBObjectDescription( TMap* const targetMap);

  /** interface function, see @ref AliHLTComponent for description */
  AliHLTComponent* Spawn();

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
  int DoEvent( const AliHLTComponentEventData& evtData,
	       const AliHLTComponentBlockData* blocks, 
	       AliHLTComponentTriggerData& trigData,
	       AliHLTUInt8_t* outputPtr, 
	       AliHLTUInt32_t& size,
	       AliHLTComponentBlockDataList& outputBlocks);

  using AliHLTProcessor::DoEvent;


  /**
   * Configure the component.
   * Parse a string for the configuration arguments and set the component
   * properties.
   */
  int Configure(const char* arguments);

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
  AliHLTGlobalCompareFlatComponent(const AliHLTGlobalCompareFlatComponent&);

  /** assignment operator prohibited */
  AliHLTGlobalCompareFlatComponent& operator=(const AliHLTGlobalCompareFlatComponent&);
void printDiff( string name, double val1, double val2);
void printDiff( string name, int n , Float_t* vals1, Float_t* vals2 );
void printDiff( string name, int n , Double_t* vals1, Double_t* vals2 );
void printDiff( string name, TString val1, TString val2);

  /*
   * ---------------------------------------------------------------------------------
   *                              Helper
   * ---------------------------------------------------------------------------------
   */

  /*
   * ---------------------------------------------------------------------------------
   *                             Members - private
   * ---------------------------------------------------------------------------------
   */
  
  /** UID for merging */

	/*
	THnSparse * fhDiff;
	
	static const	Int_t fDim = 14;
	*/
	ofstream outFile;
	ofstream conflictsFile;
	string fCurrentClass;
	
	
  ClassDef(AliHLTGlobalCompareFlatComponent, 0)
};
#endif
