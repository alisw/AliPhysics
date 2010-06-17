#ifndef ALIMULTIPLICITY_H
#define ALIMULTIPLICITY_H

#include<TObject.h>
#include <TBits.h>
#include<TMath.h>

////////////////////////////////////////////////////////
////   Class containing multiplicity information      //
////   to stored in the ESD                           //
////////////////////////////////////////////////////////

class AliMultiplicity : public TObject {

 public:

  AliMultiplicity();               // default constructor
  // standard constructor
  AliMultiplicity(Int_t ntr,Float_t *th, Float_t *ph, Float_t *dth, Float_t *dph, Int_t *labels,
         Int_t* labelsL2, Int_t ns, Float_t *ts, Float_t *ps, Short_t nfcL1, Short_t nfcL2, const TBits & fFastOrFiredChips);
  AliMultiplicity(Int_t ntr, Int_t ns, Short_t nfcL1, Short_t nfcL2, const TBits & fFastOr);
  AliMultiplicity(const AliMultiplicity& m);
  AliMultiplicity& operator=(const AliMultiplicity& m);
  virtual void Copy(TObject &obj) const;
  virtual void Clear(Option_t* opt="");
  virtual ~AliMultiplicity();
// methods to access tracklet information
  Int_t GetNumberOfTracklets() const {return fNtracks;}
  Double_t GetTheta(Int_t i) const { 
    if(i>=0 && i<fNtracks) return fTh[i];
    Error("GetTheta","Invalid track number %d",i); return -9999.;
  }
  Double_t GetEta(Int_t i) const { 
    if(i>=0 && i<fNtracks) return -TMath::Log(TMath::Tan(fTh[i]/2.));
    Error("GetEta","Invalid track number %d",i); return -9999.;
  }
  Double_t GetPhi(Int_t i) const { 
    if(i>=0 && i<fNtracks) return fPhi[i];
    Error("GetPhi","Invalid track number %d",i); return -9999.;
  }
  Double_t GetDeltaTheta(Int_t i) const {
    if(fDeltTh && i>=0 && i<fNtracks) return fDeltTh[i];
    Error("GetDeltaTheta","DeltaTheta not available in data or Invalid track number %d(max %d)",i, fNtracks); return -9999.;
  }
  Double_t GetDeltaPhi(Int_t i) const {
    if(i>=0 && i<fNtracks) return fDeltPhi[i];
    Error("GetDeltaPhi","Invalid track number %d",i); return -9999.;
  }

  Int_t GetLabel(Int_t i, Int_t layer) const;
  void  SetLabel(Int_t i, Int_t layer, Int_t label);

  Bool_t FreeClustersTracklet(Int_t i, Int_t mode) const {return (mode>=0&&mode<2&&i>=0&&i<fNtracks) ? !fUsedClusT[mode].TestBitNumber(i):kFALSE;}
  Bool_t FreeSingleCluster(Int_t i, Int_t mode)    const {return (mode>=0&&mode<2&&i>=0&&i<fNsingle) ? !fUsedClusS[mode].TestBitNumber(i):kFALSE;}

  
// methods to access single cluster information
  Int_t GetNumberOfSingleClusters() const {return fNsingle;}
  Double_t GetThetaSingle(Int_t i) const { 
    if(i>=0 && i<fNsingle) return fThsingle[i];
    Error("GetThetaSingle","Invalid cluster number %d",i); return -9999.;
  }

  Double_t GetPhiSingle(Int_t i) const { 
    if(i>=0 && i<fNsingle) return fPhisingle[i];
    Error("GetPhisingle","Invalid cluster number %d",i); return -9999.;
  }

  Short_t GetNumberOfFiredChips(Int_t layer) const { return fFiredChips[layer]; }
  void SetFiredChips(Int_t layer, Short_t firedChips) { fFiredChips[layer] = firedChips; }

  UInt_t GetNumberOfITSClusters(Int_t layer) const { return layer<6 ? fITSClusters[layer] : 0; }
  UInt_t GetNumberOfITSClusters(Int_t layMin, Int_t layMax) const ;
  void SetITSClusters(Int_t layer, UInt_t clusters) { fITSClusters[layer] = clusters; }

  void   SetFastOrFiredChips(UInt_t chipKey){fFastOrFiredChips.SetBitNumber(chipKey);}
  const TBits & GetFastOrFiredChips() const {return fFastOrFiredChips;}
  Bool_t TestFastOrFiredChips(UInt_t chipKey) const {return fFastOrFiredChips.TestBitNumber(chipKey);}

  void   SetFiredChipMap(TBits & firedChips){fClusterFiredChips = firedChips;}
  void   SetFiredChipMap(UInt_t chipKey){fClusterFiredChips.SetBitNumber(chipKey);}
  const TBits & GetFiredChipMap() const {return fClusterFiredChips;}
  Bool_t TestFiredChipMap(UInt_t chipKey) const {return fClusterFiredChips.TestBitNumber(chipKey);}



  // array getters
  Double_t* GetTheta()       const {return (Double_t*)fTh;}
  Double_t* GetPhi()         const {return (Double_t*)fPhi;}
  Double_t* GetDeltTheta()   const {return (Double_t*)fDeltTh;}
  Double_t* GetDeltPhi()     const {return (Double_t*)fDeltPhi;}
  Double_t* GetThetaSingle() const {return (Double_t*)fThsingle;}
  Double_t* GetPhiSingle()   const {return (Double_t*)fPhisingle;}
  Int_t*    GetLabels()      const {return (Int_t*)fLabels;}
  Int_t*    GetLabels2()     const {return (Int_t*)fLabelsL2;}
  
  void SetTrackletData(Int_t id, const Float_t* tlet, UInt_t bits);
  void SetSingleClusterData(Int_t id, const Float_t* scl,UInt_t bits);
  void CompactBits();

  virtual void Print(Option_t *opt="") const;

  protected:
  void Duplicate(const AliMultiplicity &m);  // used by copy ctr.

  Int_t fNtracks;            // Number of tracklets
  Int_t fNsingle;            // Number of clusters on SPD layer 1, not associated
                             // with a tracklet on SPD layer 2
  Int_t *fLabels;            //[fNtracks] array with labels of cluster in L1 used for tracklet
  Int_t *fLabelsL2;          //[fNtracks] array with labels of cluster in L2 used for tracklet
  Double32_t *fTh;           //[fNtracks] array with theta values
  Double32_t *fPhi;          //[fNtracks] array with phi values
  Double32_t *fDeltTh;       //[fNtracks] array with delta theta values
  Double32_t *fDeltPhi;      //[fNtracks] array with delta phi values
  Double32_t *fThsingle;     //[fNsingle] array with theta values of L1 clusters
  Double32_t *fPhisingle;    //[fNsingle] array with phi values of L1 clusters
  Short_t fFiredChips[2];    // Number of fired chips in the two SPD layers
  UInt_t fITSClusters[6];   // Number of ITS cluster per layer
  TBits fFastOrFiredChips;   // Map of FastOr fired chips
  TBits fClusterFiredChips;  // Map of fired chips (= at least one cluster)
  TBits fUsedClusT[2];       // flag that at least one of tracklets clusters is used by TPC/ITS+ITS_SA (0) or ITS_SA_PURE tracks (1)
  TBits fUsedClusS[2];       // flag that the single clusters is used by TPC/ITS+ITS_SA (0) or ITS_SA_PURE tracks

  ClassDef(AliMultiplicity,13);
};

inline Int_t AliMultiplicity::GetLabel(Int_t i, Int_t layer) const
{
    if(i>=0 && i<fNtracks) {
	if (layer == 0) {
	    return fLabels[i];
	} else if (layer == 1) {
	    if (fLabelsL2) {
		return fLabelsL2[i];
	    } else {
		Warning("GetLabel", "No information for layer 2 available !");
		return -9999;
	    }
	} else {
	    Error("GetLabel","Invalid layer number %d",layer); return -9999;
	}
    } else {
	Error("GetLabel","Invalid track number %d",i); return -9999;
    }
    return -9999;
}
#endif
