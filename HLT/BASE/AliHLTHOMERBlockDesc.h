//-*- Mode: C++ -*-
#ifndef ALIHLTHOMERBLOCKDESC_H
#define ALIHLTHOMERBLOCKDESC_H

/* This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/** @file   AliHLTHOMERBlockDesc.h
    @author Jochen Thaeder
    @date   
    @brief  Container for HOMER Blocks
*/

// see below for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

class AliHLTMessage;

#include "TString.h"
#include "TObject.h"

/**
 * @class AliHLTHOMERBlockDesc
 * This class contains the data which comes from 1 block, delivered via the 
 * HOMER interface. It used in order to fill these block in TLists. ( It has 
 * to inherit from TObject ). Further more it reads the specification of the 
 * block and classifies the data as TObject, raw data, or something else. Mainly 
 * used in the AliEVEHOMERManager, as data objects for AliEVE.
 * 
 * @ingroup alihlt_homer
 */
class AliHLTHOMERBlockDesc : public TObject {

public:

 /*
  * ---------------------------------------------------------------------------------
  *                            Constructor / Destructor 
  * --------------------------------------------------------------------------------- 
  */

  /** standard constructor */
  AliHLTHOMERBlockDesc();       

  /** Constructor, set block data 
   * @param data           Pointer to data
   * @param size           Size of data
   * @param origin         Detector
   * @param dataType       HLT data type  
   * @param specification  HLT specification
   */
  AliHLTHOMERBlockDesc( void * data, ULong_t size, TString origin, TString dataType, ULong_t specification );  

  /** destructor */
  virtual ~AliHLTHOMERBlockDesc();

  /*
   * ---------------------------------------------------------------------------------
   *                            Data Handling - Setter - public
   * --------------------------------------------------------------------------------- 
   */

  /** Set block data 
   * @param data           Pointer to data
   * @param size           Size of data
   * @param origin         Detector
   * @param dataType       HLT data type  
   * @param specification  HLT specification
   */
  void SetBlock( void * data, ULong_t size, TString origin, TString dataType, ULong_t specification );

  /*
   * ---------------------------------------------------------------------------------
   *                            TObject Handling - public
   * --------------------------------------------------------------------------------- 
   */

  /** Returns if block contains a TObject
   * @return           Returns kTRUE if block contains a TObject, kFALSE otherwise.
   */
  Bool_t IsTObject() { return fIsTObject;} 

  /** Returns Pointer to TObject */
  TObject* GetTObject();

  /*
   * ---------------------------------------------------------------------------------
   *                            Data Handling - Getter - public
   * --------------------------------------------------------------------------------- 
   */
  
  /** Returns if block contains a raw data, fClassName is not set in this case.
   * @return           Returns kTRUE if block contains raw data, kFALSE otherwise.
   */
  Bool_t IsRawData() { return fIsRawData;} 

  /** Get pointer to data 
   * @return           Pointer to data
   */
  void* GetData() { return fData; }

  /** Get size of data 
   * @return           Size of data
   */
  ULong_t GetSize() { return fSize; }

  /** Get detector of this block 
   * @return           Detector name
   */
  TString GetDetector() { return fDetector; }
  
  /** Get HLT data type of this block 
   * @return           HLT data type
   */
  TString GetDataType() { return fDataType; }
  
  /** Get HLT specification of this block 
   * @return           HLT specification
   */
  ULong_t GetSpecification() { return fSpecification; }

  /** Get class name of this block 
   * @return           class name
   */
  TString GetClassName() { return fClassName; }

  /** Get sub detector of this block 
   * @return           subdetector
   */
  TString GetSubDetector() { return fSubDetector; }

  /** Get sub sub detector of this block
   * @return           subsubdetector
   */
  TString GetSubSubDetector() { return fSubSubDetector; }

  /** Returns kTRUE if HLT specification indicates a subdetector range
   * @return           kTRUE if subdetector range
   */
  Bool_t HasSubDetectorRange() { return fHasSubDetectorRange; }

  /** Returns kTRUE if HLT specification indicates a subsubdetector range
   * @return           kTRUE if subsubdetector range
   */
  Bool_t HasSubSubDetectorRange() { return fHasSubSubDetectorRange; }


private:
  /** copy constructor prohibited */
  AliHLTHOMERBlockDesc(const AliHLTHOMERBlockDesc&);

  /** assignment operator prohibited */
  AliHLTHOMERBlockDesc& operator=(const AliHLTHOMERBlockDesc&);

  /** Set all additional members*/
  void SetBlockParameters();

  /** Checks if Block contains a TObject. 
   * If so, set fIsTObject to kTRUE, otherwise kFALSE
   * @return           fIsTObject
   */
  Bool_t CheckIfTObject();

  /** Checks if Block contains a TObject raw data.
   * If so, set fIsRawData to kTRUE, otherwise kFALSE
   * @return           fIsRawData
   */
  Bool_t CheckIfRawData();
  /*
   * ---------------------------------------------------------------------------------
   *                            Members - private
   * --------------------------------------------------------------------------------- 
   */

  /** Pointer to data of the block */
  void* fData;                 //! transient
	      
  /** Size of data */
  ULong_t fSize;               // see above

  /** States if block contains a TObject */
  Bool_t fIsTObject;           // see above

  /** States if block contains a raw data */
  Bool_t fIsRawData;           // see above

  /** AliHTMessage object containg a TObject */
  AliHLTMessage* fMessage;     //! transient

  /** Class Name of the block */
  TString fClassName;          // see above

  /** Detector Name, e.g. PHOS */
  TString fDetector;           // see above

  /** SubDetector Name e.g. MODULE */
  TString fSubDetector;        // see above

  /** SubSubDetector Name e.g. PARTITION */
  TString fSubSubDetector;     // see above

  /** HLT Specification */
  ULong_t fSpecification;      // see above

  /** HLT DataType */
  TString fDataType;           // see above 

  /** States id block has a subdetector range */
  Bool_t fHasSubDetectorRange;      // see above

  /** States id block has a subsubdetector range */
  Bool_t fHasSubSubDetectorRange;   // see above

  ClassDef( AliHLTHOMERBlockDesc, 0 )
};

#endif
