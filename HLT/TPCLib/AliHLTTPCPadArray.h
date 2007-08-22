// -*- Mode: C++ -*-
// @(#) $Id$

#ifndef ALIHLTPADARRAY_H
#define ALIHLTPADARRAY_H
/* This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/** @file   AliHLTTPCPadArray.h
    @author Kenneth Aamodt
    @date   
    @brief  Class containing arrays of TPC Pads.
*/

#include "AliHLTLogging.h"
#include "AliHLTTPCTransform.h"
#include "AliHLTTPCPad.h"
#include "AliHLTTPCClusters.h"
#include <vector>

typedef Int_t AliHLTTPCSignal_t;
class AliHLTTPCDigitReader;

/**
 * @class AliHLTTPCPadArray
 * TODO
 */

class AliHLTTPCPadArray : public AliHLTLogging {

public:

  /** standard constructor */
  AliHLTTPCPadArray();

  /** 
   * Constructor
   * @param patch   Patch number, either use this constructor or 
   * use the default constructor and the SetPatch method
   *
   */
  AliHLTTPCPadArray(Int_t patch);

  /** not a valid copy constructor, defined according to effective C++ style */
  AliHLTTPCPadArray(const AliHLTTPCPadArray&);
  /** not a valid assignment op, but defined according to effective C++ style */
  AliHLTTPCPadArray& operator=(const AliHLTTPCPadArray&);
  /** standard destructor */
  virtual ~AliHLTTPCPadArray();

  /**
   * Initialize the pad vector for the patch set.
   */
  Int_t InitializeVector();

  /**
   * Deinitialize the pad vector for the patch set.
   */
  Int_t DeInitializeVector();
  
  /**
   * Loop over all pads setting their data array to -1.
   */
  void DataToDefault(){
    for(Int_t i=0;i<fNumberOfRows;i++){
      for(Int_t j=0;j<fNumberOfPadsInRow[i];j++){
	fRowPadVector[i][j]->SetDataToDefault();
      }
    }
  }
  
  /**
   * Set the patch number.
   */
  void SetPatch(Int_t patch);

  /**
   * Set the digit reader.
   */
  void SetDigitReader(AliHLTTPCDigitReader* digitReader);

  /**
   * Reads the data, and set it in the Pad objects.
   */
  Int_t ReadData();

  /**
   * Retuns number of pads in this row.
   */
  Int_t GetNumberOfPads(Int_t row){return fNumberOfPadsInRow[row];}

  /**
   * Loop over all pads, checking for clustercandidates.
   */
  void FindClusterCandidates(){
    for(Int_t row=0;row<fNumberOfRows;row++){
      for(Int_t pad=0;pad<fNumberOfPadsInRow[row];pad++){
	fRowPadVector[row][pad]->FindClusterCandidates();
      }
    }
  }
  
  /**
   *
   * Loop over all pads looking for clusters, if cluster candidates on two neighbouring
   * pads have a mean time difference of <match it is said to be a cluster.
   *
   */
  void FindClusters(Int_t match);

  /**
   * Print the values of the cluster, used for debugging purposes.
   */
  void PrintClusters();

  typedef vector<AliHLTTPCPad*> fPadVector;

  vector<fPadVector> fRowPadVector;                                //! transient

  vector<AliHLTTPCClusters> fClusters;                             //! transient

private:

  /** The patch number */
  Int_t fPatch;

  Int_t fFirstRow;                                                 //! transient

  Int_t fLastRow;                                                  //! transient

  Int_t fThreshold;                                                //! transient

  Int_t* fNumberOfPadsInRow;                                       //! transient

  Int_t fNumberOfRows;                                             //! transient

  AliHLTTPCDigitReader* fDigitReader;                              //! transient

  ClassDef(AliHLTTPCPadArray, 0);
};
#endif
