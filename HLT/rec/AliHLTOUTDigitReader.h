//-*- Mode: C++ -*-
// $Id$

#ifndef ALIHLTOUTDIGITREADER_H
#define ALIHLTOUTDIGITREADER_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               */

/** @file   AliHLTOUTDigitReader.h
    @author Matthias Richter
    @date   
    @brief  HLTOUT data wrapper for simulated AliRoot HLT digit data.
*/

#include "AliHLTOUTHomerCollection.h"
#include "TString.h"

class AliRawReader;
class AliHLTHOMERReader;
class TTree;
class TFile;
class TArrayC;

/**
 * @class AliHLTOUTDigitReader
 * Handler of HLTOUT data for simulated HLT digit input.
 */
class AliHLTOUTDigitReader : public AliHLTOUTHomerCollection {
 public:
  /** constructor */
  AliHLTOUTDigitReader(int event=-1, AliHLTEsdManager* pEsdManager=NULL, const char* digitFile="HLT.Digits.root");
  /** destructor */
  virtual ~AliHLTOUTDigitReader();

 protected:
  // interface functions of AliHLTOUTHomerCollection
  Bool_t ReadNextData(UChar_t*& data);
  int Reset();
  int GetDataSize();
  const AliRawDataHeader* GetDataHeader();
  void SelectEquipment(int equipmentType, int minEquipmentId = -1, int maxEquipmentId = -1);
  int GetEquipmentId();

 private:
  /** copy constructor prohibited */
  AliHLTOUTDigitReader(const AliHLTOUTDigitReader&);
  /** assignment operator prohibited */
  AliHLTOUTDigitReader& operator=(const AliHLTOUTDigitReader&);

  /**
   * Read the data from the root file and HLTOUT raw tree.
   * Retrieve the number of branches and allocate arrays acording
   * to that. After initialization of the arrays and variables, the
   * event fEnvent is read.
   */
  bool ReadArrays();

  /**
   * Cleanup tree and data arrays.
   */
  int CloseTree();

  /**
   * Set the RunLoader as parameter
   * The function is for internal use only in conjunction with the
   * AliHLTOUT::New() functions.
   */
  void SetParam(TTree* pDigitTree, int event=-1);

  /**
   * Set name of the digit file as parameter
   * Overloaded from AliHLTOUT
   */ 
  void SetParam(const char* filename, int event=-1);

  /** name of the digit file */
  TString fDigitFileName; //! transient

  /** the root file for the HLT 'digit' output */
  TFile* fpDigitFile; //!transient

  /** the tree for the HLT 'digit' output */
  TTree* fpDigitTree; //!transient

  /** min DDL id for equipment selection */
  int fMinDDL; //!transient

  /** max DDL id for equipment selection */
  int fMaxDDL; //!transient

  /** array of digit data read from tree */
  TArrayC** fppDigitArrays; //!transient

  /** array of equipment ids for the corresponding data blocks in fppDigitArrays */
  int* fpEquipments; //!transient

  /** number of DDL objects -> size of the arrays */
  int fNofDDLs; //!transient

  /** current position in the array */
  int fCurrent; //!transient

  ClassDef(AliHLTOUTDigitReader, 0)
};
#endif
