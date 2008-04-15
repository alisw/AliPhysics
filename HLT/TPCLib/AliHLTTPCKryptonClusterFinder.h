// $Id$

#ifndef AliHLTTPC_KRYPTONCLUSTERFINDER
#define AliHLTTPC_KRYPTONCLUSTERFINDER
/* This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/** @file   AliHLTTPCKryptonClusterFinder.h
    @author Kenneth Aamodt kenneth.aamodt@student.uib.no
    @date   
    @brief  Krypton Cluster Finder for the TPC
*/

// see below for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt


//#include "AliHLTLogging.h"
//#include "AliHLTTPCPad.h"
#include "AliHLTTPCClusterFinder.h"
#include "TString.h"
#include "TH1F.h"
#include "TObjArray.h"

class AliHLTTPCSpacePointData;
class AliHLTTPCDigitReader;

class AliHLTTPCKryptonClusterFinder : public AliHLTTPCClusterFinder {

 public:
  /** standard constructor */
  AliHLTTPCKryptonClusterFinder();
  /** destructor */
  virtual ~AliHLTTPCKryptonClusterFinder();

  /** Rebunches the data, use on real data which has "wrong" bunches due to keeping 0 data */
  void ReBunch(const UInt_t * bunchData,Int_t bunchSize);

  /** rads the data insorted */
  void ReadDataUnsorted(void* ptr,unsigned long size);

  /** compare one pads combining neighbouring clustercandidates to a cluster */
  Bool_t ComparePads(AliHLTTPCPad *nextPad,AliHLTTPCClusters* cluster,Int_t nextPadToRead);

  /** Find clusters on the rows */
  void FindRowClusters();

  /** combines the row clusters to a krypton cluster */
  void FindKryptonClusters();

  /** checks if there is a candidate on the previous row */
  void CheckForCandidateOnPreviousRow(AliHLTTPCClusters* tmpCluster);

  /** initialize the histograms */
  void InitializeHistograms();

  /** resets the histograms */
  void ResetHistograms(TString histoName);

  /*set the selection from minrow to maxrow, used to look at a certain interval of rows */
  void SetSelection(Int_t minRow, Int_t maxRow);

  /** puts the histograms in the objarray*/
  void GetHistogramObjectArray(TObjArray& histos);

  /** write histograms to a file, filename must also include the path. Used only for debugging by developers */
  void WriteHistograms(TString filename);

 private: 
  /** copy constructor prohibited */
  AliHLTTPCKryptonClusterFinder(const AliHLTTPCKryptonClusterFinder&);
  /** assignment operator prohibited */
  AliHLTTPCKryptonClusterFinder& operator=(const AliHLTTPCKryptonClusterFinder&);

  vector<Int_t> fTimebinsInBunch;                                  //! transient

  vector<Int_t> fIndexOfBunchStart;                                //! transient

  //histograms
  TH1F * fHKryptonSpectrumFullPatch;                                //! transient
  TH1F * fHKryptonSpectrumSelection;                                //! transient
  TH1F * fHNumberOfKryptonClusters;                                 //! transient
  TH1F * fHNumberOfKryptonClustersSelection;                        //! transient
  TH1F * fHMaxQofKryptonClusterLast1000;                            //! transient
  TH1F * fHMaxQofKryptonClusterSelection;                           //! transient
  
  Int_t fStartBinKryptonSpectrum;                                  //! transient
  Int_t fEndBinKryptonSpectrum;                                    //! transient
  Int_t fStartBinMaxQ;                                             //! transient
  Int_t fEndBinMaxQ;                                               //! transient
  Int_t fStartBinNumberOfKryptonClusters;                          //! transient
  Int_t fEndBinNumberOfKryptonClusters;                            //! transient

  Int_t fSelectionMinRowNumber;                                    //! transient
  Int_t fSelectionMaxRowNumber;                                    //! transient

  Int_t fMaxQOfCluster;                                            //! transient
  Int_t fMaxQOfClusterBin;                                         //! transient
  
  Int_t fNumberOfKryptonClusters;                                  //! transient
  Int_t fNumberOfKryptonClustersBin;                               //! transient

  Bool_t fHistogramsInitialized;                                   //! transient                 
  
  ClassDef(AliHLTTPCKryptonClusterFinder,1) //Fast cluster finder
};
#endif
