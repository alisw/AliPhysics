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

  AliITSUTrackCond(Int_t nLayers=0);
  AliITSUTrackCond(const AliITSUTrackCond& src);
  AliITSUTrackCond &operator=(const AliITSUTrackCond& src);

  ~AliITSUTrackCond() {}
  
  void        SetNLayers(Int_t nl)                 {fNLayers = nl;}
  void        SetID(Int_t id)                      {SetUniqueID(id);}
  void        AddNewCondition(Int_t minClusters);
  void        AddGroupPattern(UShort_t patt);

  Int_t       GetID()                                  const {return GetUniqueID();}
  Int_t       GetNConditions()                         const {return fNConditions;}
  UShort_t    GetGroup(Int_t condID,Int_t grID)        const {return fConditions[fAuxData[condID*kNAuxSz+kCondStart]+grID];}
  Bool_t      CheckPattern(UShort_t patt)    const;
  //
  virtual void  Print(Option_t* option = "")           const;

 protected:
  //
  Short_t     fNLayers;                  // total number of layers
  Short_t     fNConditions;              // number of conditions defined
  TArrayS     fConditions;               //[fNConditions] set of conditions
  TArrayS     fAuxData;                  // condition beginning (1st group), n groups, min clus
  //
  ClassDef(AliITSUTrackCond,1)           // set of requirements on track hits pattern
};



#endif
