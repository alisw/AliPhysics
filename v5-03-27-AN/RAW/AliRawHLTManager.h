//-*- Mode: C++ -*-
// @(#) $Id: AliRawHLTManager.h 23318 2008-01-14 12:43:28Z hristov $

#ifndef ALIRAWHLTMANAGER_H
#define ALIRAWHLTMANAGER_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/** @file   AliRawHLTManager.h
    @author Matthias Richter
    @date   
    @brief  dynamic generation of HLT RAW readers and streams
*/

// see below for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
 
class AliRawReader;
#include "TObject.h"

extern "C" {
  typedef AliRawReader* (*AliRawReaderHLTCreateInstance_t)(AliRawReader* pParentReader, const char* options);
}
                        
/**
 * @class AliRawHLTManager
 * The class gives dynamic access to creater methods for HLT RAW readers and
 * streams without any library dependencies to HLT libraries.
 */
class AliRawHLTManager {
 public:
  AliRawHLTManager();
  virtual ~AliRawHLTManager();

  /**
   * Create an instance of the AliRawReaderHLT.
   * The AliRawReaderHLT instance needs the parent RAW reader and a list
   * of detectors for which it should access data from the HLT stream.
   */
  static AliRawReader* CreateRawReaderHLT(AliRawReader* pParent, const char* detectors);

  /**
   * Create an instance of a RAW stream.
   * There is no common base class for RAW streams due to the different nature of the
   * detectors and the data. The least common class is the TObject. The calling code
   * should check if the right class has been created by
   * <pre>
   * TObject pObject=AliRawHLTManager::CreateRawStream("MyClass");
   * MyClass* pClass=dynamic_cast<MyClass*>(pObject)
   * </pre>
   *
   * \b NOTE: The function redirects the request to the HLT framework, a handler
   * to actually create the RAW stream must be implemented in the corresponding
   * component library.
   */
  static TObject* CreateRawStream(const char* className);
 protected:
 private:
  enum {kUnloaded=0, kLoaded, kUnavailable};

  /**
   * Load the HLT interface library
   */
  static int LoadLibrary();

  /** status of the loading of the HOMER library */
  static int fLibraryStatus; //!transient

  /** entry in the HOMER library */
  static AliRawReaderHLTCreateInstance_t fFctCreateRawReaderHLT; //!transient

  /** entry in the HOMER library */
  static void* fFctCreateRawStream; //!transient

  ClassDef(AliRawHLTManager, 0)
};

// those definitions have been copied one to one from rec/AliRawReaderHLT.h
// to avoid including this header file
#define ALIHLTREC_LIBRARY                   "libHLTrec.so"
#define ALIHLTREC_LIBRARY_VERSION           0
#define ALIRAWREADERHLT_CREATE_INSTANCE     "AliRawReaderHLTCreateInstance"


#endif //ALIRAWHLTMANAGER_H
