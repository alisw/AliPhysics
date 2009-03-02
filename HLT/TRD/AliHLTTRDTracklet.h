#ifndef ALIHLTTRDTRACKLET_H
#define ALIHLTTRDTRACKLET_H

#include "AliTRDseedV1.h"
#include "AliHLTDataTypes.h"
#include "AliHLTLogging.h"
class AliHLTTRDCluster;

class AliHLTTRDTracklet
{
 public:
  AliHLTTRDTracklet();
  AliHLTTRDTracklet(AliTRDseedV1* inTracklet);
  
  void ExportTRDTracklet(AliTRDseedV1* outTracklet);
  void AddClusters();
  void CopyDataMembers();
  AliHLTUInt8_t *GetEndPointer() // Returns pointer to the end of the tracklet
    { return ((AliHLTUInt8_t *) this + fSize); };
  AliHLTUInt32_t GetSize(){ return fSize; };
  void Print(Bool_t printClusters = kTRUE);
  void ReadClustersFromMemory(void *input);
  
 private:
  AliHLTTRDTracklet(const AliHLTTRDTracklet&);
  AliHLTTRDTracklet& operator=(const AliHLTTRDTracklet&);
  void InitArrays();
  
  AliTRDseedV1* fTRDtracklet;
  AliHLTUInt32_t fSize; // Size of the tracklet with clusters in the memory
  
  /* Defenitely need */
  AliHLTTRDCluster *fClusters[AliTRDseedV1::kNTimeBins]; // Clusters
  Float_t        fYref[2];              //  Reference y
  Float_t        fZref[2];              //  Reference z
  //Float_t        fSigmaY;               //  "Robust" sigma in Y - constant fit
  Float_t        fSigmaY2;              //  "Robust" sigma in Y - line fit

  /* Probably need */
  Float_t        fTilt;                 //  Tilting angle
  Float_t        fPadLength;            //  Pad length
  Float_t        fX0;                   //  X0 position
//  Float_t        fX[knTimebins];        //! X position
//  Float_t        fY[knTimebins];        //! Y position
//  Float_t        fZ[knTimebins];        //! Z position
  Int_t          fIndexes[AliTRDseedV1::kNTimeBins];  //! Indexes
  Float_t        fYfit[2];              //  Y fit position +derivation
  //Float_t        fYfitR[2];             //  Y fit position +derivation
  Float_t        fZfit[2];              //  Z fit position
  //Float_t        fZfitR[2];             //  Z fit position
  //Float_t        fMeanz;                //  Mean vaue of z
  //Float_t        fZProb;                //  Max probbable z
  Int_t          fLabels[3];            //  Labels
  //Int_t          fN;                    //  Number of associated clusters
  Int_t          fN2;                   //  Number of not crossed
  Int_t          fNUsed;                //  Number of used clusters
  //Int_t          fFreq;                 //  Frequency
  //Int_t          fNChange;              //  Change z counter
  //Float_t        fMPads;                //  Mean number of pads per cluster

  Float_t        fC;                    //  Curvature
  //Float_t        fCC;                   //  Curvature with constrain
  Float_t        fChi2;                 //  Global chi2
  //Float_t        fChi2Z;                //  Global chi2

  /* ======= From AliTRDseedV1 ======== */

  /* Defenitely need */
  Int_t            fDet;                    //  TRD detector
  Float_t          fMom;                    //  Momentum estimate for  tracklet [GeV/c]
  Float_t          fdX;                     // length of time bin

  /* Probably need */

};

#endif
