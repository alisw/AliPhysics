// -*- Mode: C++ -*-
// @(#) $Id$

#ifndef ALIHLTTPCPADARRAY_H
#define ALIHLTTPCPADARRAY_H
/* This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/** @file   AliHLTTPCPadArray.h
    @author Kenneth Aamodt
    @date   
    @brief  Class containing arrays of TPC Pads.
*/

// see below for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include "AliHLTLogging.h"
#include "AliHLTTPCTransform.h"
#include "AliHLTTPCClusters.h"
#include <vector>

typedef Int_t AliHLTTPCSignal_t;
class AliHLTTPCDigitReader;
class AliHLTTPCPad;

/**
 * @class AliHLTTPCPadArray
 * TODO
 */

class AliHLTTPCPadArray : public AliHLTLogging {

public:

  struct AliHLTTPCActivePads
  {
    UInt_t fPad;           //pad value
    UInt_t fRow;           //row value
  };
  typedef struct AliHLTTPCActivePads AliHLTTPCActivePads; //!


  /** standard constructor */
  AliHLTTPCPadArray();

  /** 
   * Constructor
   * @param patch   Patch number, either use this constructor or 
   * use the default constructor and the SetPatch method
   *
   */
  AliHLTTPCPadArray(Int_t patch);

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
  void DataToDefault();
  
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
  Int_t GetNumberOfPads(Int_t row) const {return fNumberOfPadsInRow[row];}

  /**
   * Loop over all pads, checking for clustercandidates.
   */
  void FindClusterCandidates();
  
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

  /**
   * Set the Signal Threshold
   */
  void SetSignalThreshold(Int_t i){fSignalThreshold=i;}

  /**
   * Set the nSigma threshold
   */
  void SetNSigmaThreshold(Double_t i){fNSigmaThreshold=i;}


  /**
   * Loop over all pads adding pads with signal to active pads
   * Returns number of active pads
   *
   */
  Int_t GetActivePads(AliHLTTPCActivePads* activePads, Int_t maxActivePads);

  typedef vector<AliHLTTPCPad*> AliHLTTPCPadVector;

  vector<AliHLTTPCPadVector> fRowPadVector;                        //! transient

  vector<AliHLTTPCClusters> fClusters;                             //! transient

private:
  /** copy constructor prohibited */
  AliHLTTPCPadArray(const AliHLTTPCPadArray&);
  /** assignment operator prohibited */
  AliHLTTPCPadArray& operator=(const AliHLTTPCPadArray&);

  /** The patch number */
  Int_t fPatch;                                                    //! transient

  Int_t fFirstRow;                                                 //! transient

  Int_t fLastRow;                                                  //! transient

  //TODO: I suggest making the following UInt_t if it is never supposed to be negative. Will do!
  /* total charge of Cluster threshold*/
  Int_t fThreshold;                                                //! transient

  //TODO: I suggest making the following UInt_t if it is never supposed to be negative.
  Int_t fSignalThreshold;                                                //! transient

  //TODO: I suggest making the following UInt_t if it is never supposed to be negative.
  Double_t fNSigmaThreshold;                                                //! transient

  Int_t* fNumberOfPadsInRow;                                       //! transient

  Int_t fNumberOfRows;                                             //! transient

  AliHLTTPCDigitReader* fDigitReader;                              //! transient

  ClassDef(AliHLTTPCPadArray, 0);
};
#endif // ALIHLTTPCPADARRAY_H
