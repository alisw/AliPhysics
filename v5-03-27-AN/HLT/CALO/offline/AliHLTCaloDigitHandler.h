/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        *
 * ALICE Experiment at CERN, All rights reserved.                         *
 *                                                                        *
 * Primary Authors: Oystein Djuvsland <oysteind@ift.uib.no>               *
 *                  for The ALICE HLT Project.                            *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

#ifndef ALIHLTCALODIGITHANDLER_H
#define ALIHLTCALODIGITHANDLER_H

#include "Rtypes.h"
#include "AliHLTLogging.h"
#include "../AliHLTCaloConstantsHandler.h"

class AliDigitNew;
class AliHLTCaloGeometry;
class TTree;
class AliRunLoader;
struct AliHLTCaloDigitDataStruct;
class AliLoader;

class AliHLTCaloDigitHandler : public AliHLTLogging, public AliHLTCaloConstantsHandler
{

public:
    
    /** Constructor */
    AliHLTCaloDigitHandler(TString detName);
  
    /** Destructor */
    virtual ~AliHLTCaloDigitHandler();
    
    /** 
     * Initialise the digit handler 
     * @param runLoader is a pointer to the run loader we're going 
     * to use.
     * @return number of events on success, negative on error
     */
    virtual Int_t Init(AliRunLoader *runLoader);
    
    /** 
     * Get digits for the specified module 
     * @param module is the module number 
     * @param buffer is the data buffer to fill the digits
     * @param return number of digits 
     */
    Int_t GetDigits(Int_t module, AliHLTCaloDigitDataStruct *buffer);
    
    /** 
     * Process event
     * @param ev is the event number
     * @return 0 on success
     */
    virtual Int_t ProcessEvent(UInt_t ev);
    
    /** Return the data type produced */
    virtual AliHLTComponentDataType GetDataType() = 0;
    
protected:

    /** Convert an offline digit to a HLT digit */
    virtual Int_t ConvertDigit(AliDigitNew *digit) = 0;
  
    /** Initialise the digit array */
    void InitDigitArray();

    /** Reset the digit array */
    void ResetDigitArray();
  
    /** Run loader */
    AliRunLoader *fRunLoader;
    
    /** Detector loader */
    AliLoader *fDetLoader;
    
    /** Number of events */
    UInt_t fNumberOfEvents;
    
    /** Tree of digits */
    TTree *fDigitTree;

    /** Array of translated digits for each module */
    AliHLTCaloDigitDataStruct **fDigits;
    
    /** Number of digits in each module */
    Int_t *fDigitsInModule;
    
    /** Geometry class, must be initialised by child class */
    AliHLTCaloGeometry *fGeometry;
    
    /** The current event */
    UInt_t fCurrentEvent;

private:
  
  
    /** Constructor Prohibited*/
    AliHLTCaloDigitHandler();

    /** Prohibited */
    AliHLTCaloDigitHandler(const AliHLTCaloDigitHandler& );

    /** Prohibited */
    AliHLTCaloDigitHandler& operator=(const AliHLTCaloDigitHandler& );
    
    ClassDef(AliHLTCaloDigitHandler, 0);
};

#endif // ALIHLTCALODIGITHANDLER_H
