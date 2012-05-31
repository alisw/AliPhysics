#ifndef ALIHLTITSCLUSTERFINDERSSD_H
#define ALIHLTITSCLUSTERFINDERSSD_H
/* Copyright(c) 2007-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliHLTITSClusterFinderSSD.h 32604 2009-05-29 10:41:46Z masera $ */

//--------------------------------------------------------------
//                       ITS clusterer V2 for SSD
//
//   This can be a "wrapping" for the V1 cluster finding classes
//   if compiled with uncommented "#define V1" line 
//   in the AliITSclustererV2.cxx file.
//
//   Origin: Iouri Belikov, CERN, Jouri.Belikov@cern.ch 
//--------------------------------------------------------------
#include "AliITSClusterFinder.h"
#include <vector>

class TClonesArray;
class AliRawReader;
class AliITSRawStream;
class AliITSRawStreamSSD;
class AliITSCalibrationSSD;
class AliITSRecoParam;

class AliHLTITSClusterFinderSSD : public AliITSClusterFinder {
public:
  AliHLTITSClusterFinderSSD(AliITSDetTypeRec* dettyp, AliRawReader *reader);
  virtual ~AliHLTITSClusterFinderSSD();
  virtual void FindRawClusters(Int_t /*mod*/){;}
  virtual void RawdataToClusters(AliRawReader* /*rawReader*/,TClonesArray** /*clusters*/){;}

  virtual void RawdataToClusters( std::vector<AliITSRecPoint> &clusters);
 protected:

  void FindClustersSSD(Ali1Dcluster* neg, Int_t nn, 
		       Ali1Dcluster* pos, Int_t np,
		       std::vector<AliITSRecPoint> &);

  AliITSRecoParam *fRecoParam; //!
  AliRawReader *fRawReader; //!
  AliITSRawStreamSSD *fRawStream;//!

  Int_t fLastSSD1;        //index of the last SSD1 detector   
  static Short_t* fgPairs;       //array used to build positive-negative pairs
  static Int_t    fgPairsSize;    //actual size of pairs array
  static const Float_t fgkCosmic2008StripShifts[16][9]; // Shifts for 2007/2008 Cosmic data (timing problem)
  static const Float_t fgkThreshold; // threshold for the seed

 private:
  AliHLTITSClusterFinderSSD(const AliHLTITSClusterFinderSSD& );
  AliHLTITSClusterFinderSSD& operator=(const AliHLTITSClusterFinderSSD&  );

  ClassDef(AliHLTITSClusterFinderSSD,0)  // ITS cluster finder V2 for SDD
};

#endif
