#ifndef ALIITSUTRACKCOND_H
#define ALIITSUTRACKCOND_H

#include <TObject.h>
#include <TArrayS.h>

//------------------------------------------------------------------------------
//
// This class defines a set of hit patterns (conditions) to consider the track reconstructable
// Each condidition (few consequitive elements from fConditions array) is set of bit patterns, 
// with each element of condition defining a group of layers which must be present in the tracks.
// For instance, if we require the track to have contributions from
// {lr0 or lr1 or lr2} AND {lr2 or lr 3 or lr 4} AND {lr5 or lr6} then the condition should
// be {BIT(0)|BIT(1)|BIT(2), BIT(2)|BIT(3)|BIT(4), BIT(5)|BIT(6)}.
// Additionally, each condition may request min number of hits to be present
//
// Each AliITSUTrackCond should correspond to single track finding pass and may contain multiple
// conditions. To consider the track reconstructable it is enough to satisfy 1 condition
//
//------------------------------------------------------------------------------


class AliITSUTrackCond : public TObject
{
 public:
  enum {kCondStart,kNGroups,kMinClus,kNAuxSz};
  //
  AliITSUTrackCond(Int_t nLayers=0);
  AliITSUTrackCond(const AliITSUTrackCond& src);
  AliITSUTrackCond &operator=(const AliITSUTrackCond& src);

  ~AliITSUTrackCond() {}
  
  void        SetNLayers(Int_t nl);
  void        SetClSharing(Int_t lr, Char_t v=0)   {fClSharing[lr] = v;}
  void        SetMaxBranches(Int_t lr, Int_t mb)   {fMaxBranches[lr] = mb;}
  void        SetMaxCandidates(Int_t lr, Int_t nc) {fMaxCandidates[lr] = nc;}
  void        SetID(Int_t id)                      {SetUniqueID(id);}
  void        AddNewCondition(Int_t minClusters);
  void        AddGroupPattern(UShort_t patt);

  Int_t       GetID()                                  const {return GetUniqueID();}
  Int_t       GetMaxBranches(Int_t lr)                 const {return fMaxBranches[lr];}
  Int_t       GetMaxCandidates(Int_t lr)               const {return fMaxCandidates[lr];}
  Int_t       GetNConditions()                         const {return fNConditions;}
  UShort_t    GetGroup(Int_t condID,Int_t grID)        const {return fConditions[fAuxData[condID*kNAuxSz+kCondStart]+grID];}
  Bool_t      CheckPattern(UShort_t patt)    const;
  //
  virtual void  Print(Option_t* option = "")           const;

  void        SetMaxITSTPCMatchChi2(Float_t v)               {fMaxITSTPCMatchChi2 = v;}
  void        SetMaxITSSAChi2(Float_t v)                     {fMaxITSSAChi2 = v;}
  void        SetMaxTr2ClChi2(Int_t lr, Float_t v)           {fMaxTr2ClChi2[lr] = v;}
  void        SetMaxChi2GloNrm(Int_t lr, Float_t v)           {fMaxChi2GloNrm[lr] = v;}
  void        SetMissPenalty(Int_t lr,  Float_t v)           {fMissPenalty[lr] = v;}
  void        SetNSigmaRoadY(Int_t lr,  Float_t v)           {fNSigmaRoadY[lr] = v;}
  void        SetNSigmaRoadZ(Int_t lr,  Float_t v)           {fNSigmaRoadZ[lr] = v;}
  //
  Float_t     GetMaxITSTPCMatchChi2()                  const {return fMaxITSTPCMatchChi2;}
  Float_t     GetMaxITSSAChi2()                        const {return fMaxITSSAChi2;}
  Char_t      GetClSharing(Int_t lr)                   const {return fClSharing[lr];}
  Float_t     GetMissPenalty(Int_t lr)                 const {return fMissPenalty[lr];}
  Float_t     GetMaxTr2ClChi2(Int_t lr)                const {return fMaxTr2ClChi2[lr];}
  Float_t     GetMaxChi2GloNrm(Int_t lr)               const {return fMaxChi2GloNrm[lr];}
  Float_t     GetNSigmaRoadY(Int_t lr)                 const {return fNSigmaRoadY[lr];}
  Float_t     GetNSigmaRoadZ(Int_t lr)                 const {return fNSigmaRoadZ[lr];}
  Bool_t      IsLayerExcluded(Int_t lr)                const {return GetMaxTr2ClChi2(lr)<=0;}
  //
  void        Init();
  Bool_t      IsInitDone()                             const {return fInitDone;}
  //
 protected: 
  //
  Bool_t      fInitDone;                 // initialization flag
  Int_t       fNLayers;                  // total number of layers
  Float_t     fMaxITSTPCMatchChi2;       // max chi2 for ITS/TPC matching
  Float_t     fMaxITSSAChi2;             // max chi2 for ITS standalone fit (backward)
  Char_t*     fClSharing;                // [fNLayers] is cluster sharing allowed
  Short_t*    fMaxBranches;              // [fNLayers] max allowed branches per seed on each layer
  Short_t*    fMaxCandidates;            // [fNLayers] max allowed candidates per TPC seed on each layer
  Float_t*    fMaxTr2ClChi2;             // [fNLayers] max track-to-cluster chi2
  Float_t*    fMaxChi2GloNrm;            // [fNLayers] max norm global chi2
  Float_t*    fMissPenalty;              // [fNLayers] chi2 penalty for missing hit on the layer
  Float_t*    fNSigmaRoadY;              // [fNLayers] number of sigmas in Y
  Float_t*    fNSigmaRoadZ;              // [fNLayers] number of sigmas in Z
  //
  Short_t     fNConditions;              // number of conditions defined
  TArrayS     fConditions;               // fNConditions  set of conditions
  TArrayS     fAuxData;                  // condition beginning (1st group), n groups, min clus
  //
  static Char_t  fgkClSharing;           // def cl.sharing allowed level
  static Int_t   fgkMaxBranches;          // def max number of branches per seed on current layer 
  static Int_t   fgkMaxCandidates;        // def max number of total candidates on current layer 
  static Float_t fgkMaxTr2ClChi2;         // def track-to-cluster chi2 cut
  static Float_t fgkMaxChi2GloNrm;        // def global norm chi2 cut
  static Float_t fgkMissPenalty;          // penalty for missing cluster
  static Float_t fgkMaxMatchChi2;         // max acceptable matching chi2
  static Float_t fgkMaxITSSAChi2;         // max acceptable standalone ITS backward fit chi2
  //
  ClassDef(AliITSUTrackCond,5)           // set of requirements on track hits pattern
};



#endif
