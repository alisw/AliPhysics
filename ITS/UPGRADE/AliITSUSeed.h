#ifndef ALIITSUSEED_H
#define ALIITSUSEED_H

#include "AliExternalTrackParam.h"
#include "AliITSUAux.h"
using namespace AliITSUAux;


class AliITSUSeed: public AliExternalTrackParam
{
 public:
  AliITSUSeed();
  AliITSUSeed(const AliITSUSeed& src);
  AliITSUSeed &operator=(const AliITSUSeed &src);
  virtual ~AliITSUSeed();
  virtual void    Print(Option_t* option = "") const;
  //
  void            SetLrClusterID(Int_t lr, Int_t cl);
  void            SetLr(Int_t lr)                        {SetLrClusterID(lr,-1);} // lr w/o cluster
  void            SetLrClusterID(UInt_t id)              {fClID = id;}
  void            SetParent(TObject* par)                {fParent = par;}
  void            SetChi2Cl(Double_t v)                  {fChi2Glo += fChi2Cl= v;}
  //
  UInt_t          GetLrClusterID()                 const {return fClID;}
  Int_t           GetLrCluster(Int_t &lr)          const {return UnpackCluster(fClID,lr);}
  Int_t           GetLayerID()                     const {return UnpackLayer(fClID);}
  Int_t           GetClusterID()                   const {return UnpackCluster(fClID);}
  Bool_t          HasClusterOnLayer(Int_t lr)      const {return fHitsPattern&(0x1<<lr);}
  Int_t           GetNLayersHit()                  const {return NumberOfBitsSet(fHitsPattern);}
  UShort_t        GetHitsPattern()                 const {return fHitsPattern;}
  Float_t         GetChi2Cl()                      const {return fChi2Cl;}
  Float_t         GetChi2Glo()                     const {return fChi2Glo;}
  //
  TObject*        GetParent()                      const {return fParent;}
  //
 protected:
  UShort_t              fHitsPattern;       // bit pattern of hits
  UInt_t                fClID;              // packed cluster info (see AliITSUAux::PackCluster)
  Float_t               fChi2Glo;           // current chi2 global
  Float_t               fChi2Cl;            // track-cluster chi2
  TObject*              fParent;            // parent track (in higher tree hierarchy)
  
  ClassDef(AliITSUSeed,1)
};

//_________________________________________________________________________
inline void AliITSUSeed::SetLrClusterID(Int_t lr, Int_t cl)
{
  // assign layer, cluster (if -1 - no hit on this layer)
  fClID = PackCluster(lr,cl);
  if (cl>=0) fHitsPattern |= 0x1<<lr;
}


#endif
