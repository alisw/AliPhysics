//-*- Mode: C++ -*-

// $Id$

#ifndef ALIHLTTPCCALIBRATIONAGENT_H
#define ALIHLTTPCCALIBRATIONAGENT_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/** @file   AliHLTTPCCalibrationAgent.h
    @author Kalliopi Kanaki
    @date   
    @brief  Agent of the libAliHLTTPCCalibration library
*/

#include "AliHLTModuleAgent.h"

// raw data handler of HLTOUT data
#include "AliHLTOUTHandlerEquId.h"

/**
 * @class AliHLTTPCCalibrationAgent
 * This is the agent for the AliHLTTPCCalibration library.<br>
 *
 * @ingroup alihlt_tpc
 */
class AliHLTTPCCalibrationAgent : public AliHLTModuleAgent {
 public:
  
  /*
   * ---------------------------------------------------------------------------------
   *                            Constructor / Destructor
   * ---------------------------------------------------------------------------------
   */

  /**
   * standard constructor. The agent is automatically registered in the
   * global agent manager
   */
  AliHLTTPCCalibrationAgent();
  
  /** destructor */
  virtual ~AliHLTTPCCalibrationAgent();

  /*
   * ---------------------------------------------------------------------------------
   *                            
   * ---------------------------------------------------------------------------------
   */

  /**
   * Register all configurations belonging to the sample library with the
   * AliHLTConfigurationHandler. The agent can adapt the configurations
   * to be registered to the current AliRoot setup by checking the
   * runloader.
   * @param [in] handler   the configuration handler
   * @param [in] rawReader AliRoot RawReader instance 
   * @param [in] runloader AliRoot runloader
   * @return neg. error code if failed
   */
  Int_t CreateConfigurations(AliHLTConfigurationHandler* handler,
			     AliRawReader* rawReader=NULL,
			     AliRunLoader* runloader=NULL) const;

  /**
   * Get the top configurations for local event reconstruction.
   * A top configuration describes a processing chain. It can simply be
   * described by the last configuration(s) in the chain. 
   * The agent can adapt the configurations to be registered to the current
   * AliRoot setup by checking the runloader.
   * @param [in] rawReader AliRoot RawReader instance 
   * @param [in] runloader AliRoot runloader
   * @return string containing the top configurations separated by blanks
   */
  const Char_t* GetReconstructionChains(AliRawReader* rawReader=NULL,
					AliRunLoader* runloader=NULL) const;

  /**
   * Component libraries which the configurations of this agent depend on.
   * @return list of component libraries as a blank-separated string.
   */
  const Char_t* GetRequiredComponentLibraries() const;

  /**
   * Register components for the AliHLTSample library.
   * @param [in] pHandler  instance of the component handler          
   */
  Int_t RegisterComponents(AliHLTComponentHandler* pHandler) const;

  /*
   * ---------------------------------------------------------------------------------
   *                            
   * ---------------------------------------------------------------------------------
   */

  /**
   *
   */
  AliHLTModulePreprocessor* GetPreprocessor();
  
  /*
   * ---------------------------------------------------------------------------------
   *                            
   * ---------------------------------------------------------------------------------
   */

  /**
   *
   */
  Int_t GetHandlerDescription(AliHLTComponentDataType dt,
			      AliHLTUInt32_t spec,
			      AliHLTOUTHandlerDesc& desc) const;

  /**
   *
   */
  AliHLTOUTHandler* GetOutputHandler(AliHLTComponentDataType dt,
				     AliHLTUInt32_t spec);

  /**
   *
   */
  Int_t DeleteOutputHandler(AliHLTOUTHandler* pInstance);

  /*
     class AliHLTOUTSDDRawDataHandler: public AliHLTOUTHandlerEquId {
     public:
     AliHLTOUTSDDRawDataHandler() {}
     ~AliHLTOUTSDDRawDataHandler() {}
     int ProcessData(AliHLTOUT* pData);
     private:
     };
  */

 protected:

 private:

  /*
   * ---------------------------------------------------------------------------------
   * Private functions to implement AliHLTComponent's interface.
   * These functions provide initialization as well as the actual processing
   * capabilities of the component. 
   * ---------------------------------------------------------------------------------
   */
  
  /** Copy constructor prohibited */
  AliHLTTPCCalibrationAgent(const AliHLTTPCCalibrationAgent&);
  
  /** Assignment operator prohibited */
  AliHLTTPCCalibrationAgent& operator=(const AliHLTTPCCalibrationAgent&);

  /*
   * ---------------------------------------------------------------------------------
   *                                     Members
   * ---------------------------------------------------------------------------------
   */

  /** Handler for TPC calibration data in the HLTOUT stream */
  AliHLTOUTHandlerEquId* fRawDataHandler; //!transient

  /** ROOT specific member definition */
  ClassDef(AliHLTTPCCalibrationAgent, 0);
};

#endif
