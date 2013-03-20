
//-*- Mode: C++ -*-
// $Id$

#ifndef ALIHLTHOMERLIBMANAGER_H
#define ALIHLTHOMERLIBMANAGER_H
/* This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/// @file   AliHLTHOMERLibManager.h
/// @author Matthias Richter
/// @date   
/// @brief  dynamic HLT HOMER reader/writer generation and destruction

#include "AliHLTDataTypes.h"
#include "TObject.h" // for ClassDef/Imp

class AliHLTHOMERReader;
class AliHLTHOMERWriter;

/**
 * @class AliHLTHOMERLibManager
 * Dynamic manager of HOMER library.
 * The class allows to generate objects of HOMER readers and writers
 * dynamically and loads also the HOMER lib. In order to write HOMER library
 * independent code it is important to use the AliHLTMonitoringWriter/
 * AliHLTMonitoringReader classes when ever class methods are used. Those
 * classes just define a virtual interface. <br>
 *
 * Instead of creating a reader or writer by \em new and deleting it with
 * \em delete, one has to use the Open and Delete methods of this class.
 *
 * <pre>
 * AliHLTHOMERLibManager manager;
 *
 * // open a HOMER reader listening at port 23000 of the localhost
 * AliHLTMonitoringReader* pReader=manager.OpenReader(localhost, 23000);
 *
 * // read next event, timeout 5s
 * while (pReader && pReader->ReadNextEvent(5000000)==0) {
 *   unsigned long count=pReader->GetBlockCnt();
 *   gSystem->Sleep(5);
 *   ...
 * }
 *
 * // delete reader
 * manager.DeleteReader(pReader);
 * </pre>
 *
 * The manager does not provide methods to create a HOMER reader on
 * basis of shared memory. This is most likely a depricated functionality,
 * although kept for the sake of completeness. However, at some point it
 * might become useful. Please notify the developers if you need that
 * functionality.
 *
 * @ingroup alihlt_homer
 */
class AliHLTHOMERLibManager {
 public:
  /** standard constructor */
  AliHLTHOMERLibManager();
  /** destructor */
  virtual ~AliHLTHOMERLibManager();

  /**
   * Open a homer reader working on a TCP port.
   */
  AliHLTHOMERReader* OpenReader(const char* hostname, unsigned short port );
  
  /**
   * Open a homer reader working on multiple TCP ports.
   */
  AliHLTHOMERReader* OpenReader(unsigned int tcpCnt, const char** hostnames, unsigned short* ports);
	
  /**
   * Open a HOMER reader for reading from a System V shared memory segment.
  AliHLTHOMERReader* OpenReader(key_t shmKey, int shmSize );
   */
	
  /**
   * Open a HOMER reader for reading from multiple System V shared memory segments
  AliHLTHOMERReader* OpenReader(unsigned int shmCnt, key_t* shmKey, int* shmSize );
   */
	
  /**
   * Open a HOMER reader for reading from multiple TCP ports and multiple System V shared memory segments
  AliHLTHOMERReader* OpenReader(unsigned int tcpCnt, const char** hostnames, unsigned short* ports, 
				    unsigned int shmCnt, key_t* shmKey, int* shmSize );
   */

  /**
   * Open a HOMER reader.
   * Load HOMER library dynamically and create object working on the provided
   * buffer.
   */
  AliHLTHOMERReader* OpenReaderBuffer(const AliHLTUInt8_t* pBuffer, int size);

  /**
   * Delete a HOMER reader.
   * Clean-up of the object is done inside the HOMER library.
   */
  int DeleteReader(AliHLTHOMERReader* pReader);

  /**
   * Open a HOMER writer.
   * Load HOMER library dynamically and create object working on the provided
   * buffer.
   */
  AliHLTHOMERWriter* OpenWriter();

  /**
   * Delete a HOMER writer.
   * Clean-up of the object is done inside the HOMER library.
   */
  int DeleteWriter(AliHLTHOMERWriter* pWriter);

 protected:

 private:
  /** copy constructor prohibited */
  AliHLTHOMERLibManager(const AliHLTHOMERLibManager&);
  /** assignment operator prohibited */
  AliHLTHOMERLibManager& operator=(const AliHLTHOMERLibManager&);

  /**
   * Load the HOMER library.
   */
  int LoadHOMERLibrary();

  /**
   * Unloads the HOMER library.
   */
  int UnloadHOMERLibrary();

  /** status of the loading of the HOMER library */
static  int fgLibraryStatus; //!transient

  /** entry in the HOMER library */
  void (*fFctCreateReaderFromTCPPort)(); //!transient

  /** entry in the HOMER library */
  void (*fFctCreateReaderFromTCPPorts)(); //!transient

  /** entry in the HOMER library */
  void (*fFctCreateReaderFromBuffer)(); //!transient

  /** entry in the HOMER library */
  void (*fFctDeleteReader)(); //!transient

  /** entry in the HOMER library */
  void (*fFctCreateWriter)(); //!transient

  /** entry in the HOMER library */
  void (*fFctDeleteWriter)(); //!transient

  /** Indicates the library that was actually (and if) loaded in LoadHOMERLibrary(). */
  const char* fLoadedLib;  //!transient

  static const char* fgkLibraries[]; /// List of libraries to try and load.
  static int fgkLibRefCount[]; /// The library reference count to control when to unload the library.

  ClassDef(AliHLTHOMERLibManager, 0)
};
#endif
