#ifndef ALIHLTTRDCLUSTER_H
#define ALIHLTTRDCLUSTER_H

#include "AliTRDcluster.h"
#include "AliHLTDataTypes.h"

class AliHLTTRDCluster
{
 public:
  AliHLTTRDCluster();
  AliHLTTRDCluster(AliTRDcluster* inCluster);
  void ExportTRDCluster(AliTRDcluster* outCluster);

  AliHLTUInt8_t *GetEndPointer() // Returns pointer to the end of the cluster
  { return ((AliHLTUInt8_t *) this + sizeof(*this)); };
  AliHLTUInt32_t GetSize(){ return sizeof(*this); };
  void Print();      
  
 private:
  // From AliCluster
  Float_t  fX;        // X of the cluster in the tracking c.s.
  Float_t  fY;        // Y of the cluster in the tracking c.s.
  Float_t  fZ;        // Z of the cluster in the tracking c.s.
  Float_t  fQ;        //  Amplitude 

  Bool_t  fIsInChamber;
  Bool_t  fIsShared;
  Short_t fDetector;       //  TRD detector number
  Char_t  fLocalTimeBin;   //  T0-calibrated time bin number
  UChar_t fClusterMasking; //  Bit field containing cluster status information;

  // From AliTRDcluster
  UChar_t fPadCol;         //  Central pad number in column direction 
  UChar_t fPadRow;         //  Central pad number in row direction 
  UChar_t fPadTime;        //  Uncalibrated time bin number 
  //   Short_t fSignals[7];     //  Signals in the cluster 
  //   UChar_t fNPads;          //  Number of pads in cluster 
  //   Float_t fCenter;         //  Center of the cluster relative to the pad  
   
};

#endif
