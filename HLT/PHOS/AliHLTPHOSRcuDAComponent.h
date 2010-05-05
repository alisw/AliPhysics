//-*- Mode: C++ -*-
// $Id$

/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * All rights reserved.                                                   *
 *                                                                        *
 * Primary Authors: Oystein Djuvsland                                     *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          * 
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/


#ifndef ALIHLTPHOSRCUDACOMPONENT_H
#define ALIHLTPHOSRCUDACOMPONENT_H

/**
 * Module calibration processor component for PHOS HLT
 *
 * @file   AliHLTPHOSRcuDAComponent.h
 * @author Oystein Djuvsland
 * @date   
 * @brief  Module calibration processor component for PHOS HLT
*/

// see below for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include "AliHLTCalibrationProcessor.h"

// #include "AliHLTPHOSConstant.h"
//#include "AliHLTPHOSRcuProperties.h"

using namespace PhosHLTConst;

class AliPHOSRcuDA1;
class AliHLTPHOSSharedMemoryInterfacev2;
class TObjArray;

//class AliHLTPHOSRcuDAComponent: public AliHLTCalibrationProcessor,  public AliHLTPHOSRcuProcessor
//class AliHLTPHOSRcuDAComponent:  public AliHLTPHOSRcuProperties, public AliHLTCalibrationProcessor


class AliHLTPHOSRcuDAComponent:  public AliHLTCalibrationProcessor
{
public:
  /** constructor */
  AliHLTPHOSRcuDAComponent();
  /** destructor */
  virtual ~AliHLTPHOSRcuDAComponent();
 /*  virtual int Deinit(); */
/*   virtual int DoInit(int, const char**); */

  // Public functions to implement AliHLTComponent's interface.
  // These functions are required for the registration process

  const char* GetComponentID();
  void GetInputDataTypes( vector<AliHLTComponentDataType>& list);
  AliHLTComponentDataType GetOutputDataType();
  virtual void GetOutputDataSize( unsigned long& constBase, double& inputMultiplier );
  AliHLTComponent* Spawn();

protected:
  // Protected functions to implement AliHLTComponent's interface.
  // These functions provide initialization as well as the actual processing
  // capabilities of the component. 
      
  /** Initialize the calibration component. */
  Int_t InitCalibration();

  /** Scan commandline arguments of the calibration component. */
  Int_t ScanArgument( Int_t argc, const char** argv );

  /** DeInitialize the calibration component. */
  Int_t DeinitCalibration();

  /** Process the data in the calibration component. */
  using  AliHLTCalibrationProcessor::ProcessCalibration;
  Int_t ProcessCalibration( const AliHLTComponentEventData& evtData, AliHLTComponentTriggerData& trigData );

  /** Ship the data to the FXS at end of run or eventmodulo. */
  using AliHLTCalibrationProcessor::ShipDataToFXS; 
  Int_t ShipDataToFXS( const AliHLTComponentEventData& evtData, AliHLTComponentTriggerData& trigData );

  /** Reset the energy and timing arrays */
  //   void ResetArrays(Float_t e[NXCOLUMNSMOD][NZROWSMOD][NGAINS], Float_t t[NXCOLUMNSMOD][NZROWSMOD][NGAINS]);

  //  int fCnt; ///TODO, remove this

private:
    /** not a valid copy constructor, defined according to effective C++ style */
  AliHLTPHOSRcuDAComponent(const AliHLTPHOSRcuDAComponent&);
  /** not a valid assignment op, but defined according to effective C++ style */
  AliHLTPHOSRcuDAComponent& operator=(const AliHLTPHOSRcuDAComponent&);
  int fPhosEventCount;  ///TODO, remove this

  /** Pointer to the DA */
  AliPHOSRcuDA1* fPHOSDAPtr;   //! transient

  /** Interface to read altro channel data from shared memory */
  AliHLTPHOSSharedMemoryInterfacev2 *fShmPtr;   //! transient

};

#endif
  
