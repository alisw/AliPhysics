//-*- Mode: C++ -*-
// @(#) $Id$

#ifndef ALIHLTMEMORYFILE_H
#define ALIHLTMEMORYFILE_H
/* This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/** @file   AliHLTMemoryFile.h
    @author Matthias Richter
    @date   
    @brief  Serialization of complete ROOT files.

// see below for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
                                                                          */
#include "TFile.h"
#include "AliHLTLogging.h"

/**
 * @class AliHLTMemoryFile
 * Serialization of ROOT files for transport in the Alice HLT analysis
 * chain.
 *
 * The file behaves like a normal ROOT file except that it is written to
 * a memory buffer instead of disk.
 */
class AliHLTMemoryFile : public TFile, public AliHLTLogging {
 public:
  /** standard constructor */
  AliHLTMemoryFile();

  /** constructor */
  AliHLTMemoryFile(void* pBuffer, int iSize);

  /** standard destructor */
  virtual ~AliHLTMemoryFile();

  /**
   * Write a header at the beginning of the file.
   * The header is not part of the ROOT file. It should be written before any
   * other object in order to avoid data to be moved.
   * @param pHeader     buffer to write
   * @param iSize       size of the buffer
   * @return neg. error code if failed
   *         - -ENOSPC    buffer size too small
   */
  int WriteHeaderBuffer(const char* pHeader, int iSize);

  /**
   * Write a header at the beginning of the file.
   * The trailer is not part of the ROOT file. It can only be written if the
   * file already has been closed.
   * @param pTrailer    buffer to write
   * @param iSize       size of the buffer
   * @return neg. error code if failed
   *         - -ENOSPC    buffer size too small
   */
  // not yet stable
  //int WriteTrailerBuffer(const char* pTrailer, int size);

  /**
   * Close file and flush output.
   * Forwards to @ref CloseMemoryFile and is introduced to avoid
   * compilation warnings about hidden virtual functions in the base class.
   */
  void Close(const Option_t*);

  /**
   * Close file and flush output.
   * @param bFlush       write remaining data
   * @return neg. error code if failed
   *         - -ENOSPC    buffer size too small
   */
  int CloseMemoryFile(int bFlush=1);

  /**
   * Check if file has been closed.
   * @return 1 if closed
   */
  int IsClosed() const {return fbClosed;}

  /**
   * Get the last error code.
   * @return error code
   */
  int GetErrno() const {return fErrno;}

  /**
   * Get header size.
   */
  int GetHeaderSize() const {return fHeaderSize;}

 protected:
  // Interface to basic system I/O routines
  Int_t    SysOpen(const char *pathname, Int_t flags, UInt_t mode);
  Int_t    SysClose(Int_t fd);
  Int_t    SysRead(Int_t fd, void *buf, Int_t len);
  Int_t    SysWrite(Int_t fd, const void *buf, Int_t len);
  Long64_t SysSeek(Int_t fd, Long64_t offset, Int_t whence);
  Int_t    SysStat(Int_t fd, Long_t *id, Long64_t *size, Long_t *flags, Long_t *modtime);
  Int_t    SysSync(Int_t fd);

 private:
  /** copy constructor prohibited */
  AliHLTMemoryFile(const AliHLTMemoryFile&);
  /** assignment operator prohibited */
  AliHLTMemoryFile& operator=(const AliHLTMemoryFile&);

  /** target buffer */
  char* fpBuffer;                                                 //! transient

  /** size of buffer */
  int fBufferSize;                                                // see above

  /** position */
  int fPosition;                                                  // see above

  /** filled posrtion of the buffer */
  int fSize;                                                      // see above

  /** result of last operation */
  int fErrno;                                                     // see above

  /** file closed */
  int fbClosed;                                                   // see above

  /** size of header */
  int fHeaderSize;                                                // see above

  /** size of trailer */
  int fTrailerSize;                                               // see above

  ClassDef(AliHLTMemoryFile, 1)
};
#endif // ALIHLTMEMORYFILE_H
