#ifndef ALIMILLEPEDERECORD_H
#define ALIMILLEPEDERECORD_H

/**********************************************************************************************/
/* AliMillePedeRecords: class to store the data of single track processing                    */
/* Format: for each measured point the data is stored consequtively                           */
/* INDEX                                                      VALUE                           */
/* -1                                                         residual                        */
/* Local_param_id                                             dResidual/dLocal_param          */
/* ...                                                        ...                             */
/* -2                                                         weight of the measurement       */
/* Global_param_od                                            dResidual/dGlobal_param         */
/* ...                                                        ...                             */
/*                                                                                            */
/* The records for all processed tracks are stored in the temporary tree in orgder to be      */
/* reused for multiple iterations of MillePede                                                */
/*                                                                                            */
/* Author: ruben.shahoyan@cern.ch                                                             */
/*                                                                                            */
/**********************************************************************************************/
#include <TObject.h>

class AliMillePedeRecord : public TObject
{
 public:
  AliMillePedeRecord();
  AliMillePedeRecord(const AliMillePedeRecord& src);
  AliMillePedeRecord& operator=(const AliMillePedeRecord& rhs);
  //
  virtual    ~AliMillePedeRecord();
  void       Reset();
  void       Print(const Option_t *opt="")                   const;
  //
  Int_t      GetSize()                                       const {return fSize;}
  Int_t     *GetIndex()                                      const {return fIndex;}
  Int_t      GetIndex(int i)                                 const {return fIndex[i];}
  //
  void       GetIndexValue(Int_t i,Int_t &ind,Double_t &val) const {ind=fIndex[i]; val=fValue[i];}
  void       AddIndexValue(Int_t ind, Double_t val);
  void       AddResidual(Double_t val)                             {AddIndexValue(-1,val);}
  void       AddWeight(Double_t val)                               {AddIndexValue(-2,val);}
  Bool_t     IsResidual(Int_t i)                             const {return fIndex[i]==-1;}
  Bool_t     IsWeight(Int_t i)                               const {return fIndex[i]==-2;}
  //
  Double_t  *GetValue()                                      const {return fValue;}
  Double_t   GetValue(Int_t i)                               const {return fValue[i];}
  //
  void       MarkGroup(Int_t id);
  Int_t      GetNGroups()                                    const {return fNGroups;}
  Int_t      GetGroupID(Int_t i)                             const {return fGroupID[i];}
  Bool_t     IsGroupPresent(Int_t id)                        const;
  //
 protected:
  Int_t      GetDtBufferSize()                               const {return GetUniqueID()&0x0000ffff;}
  Int_t      GetGrBufferSize()                               const {return GetUniqueID()>>16;}
  void       SetDtBufferSize(Int_t sz)                             {SetUniqueID((GetGrBufferSize()<<16)+sz);}
  void       SetGrBufferSize(Int_t sz)                             {SetUniqueID(GetDtBufferSize()+(sz<<16));}
  void       ExpandDtBuffer(Int_t bfsize);
  void       ExpandGrBuffer(Int_t bfsize);
  //
 protected:
  Int_t      fSize;                             // size of the record
  Int_t      fNGroups;                          // number of groups (e.g. detectors) contributing
  Int_t   *  fGroupID;                          //[fNGroups] groups id's (in increasing order)
  Int_t   *  fIndex;                            //[fSize] index of variables
  Double_t*  fValue;                            //[fSize] array of values: derivs,residuals
  //
  ClassDef(AliMillePedeRecord,1)                // Record of track residuals and local/global deriavtives
};

inline void  AliMillePedeRecord::AddIndexValue(Int_t ind, Double_t val) 
{
  if (fSize>=GetDtBufferSize()) ExpandDtBuffer(2*(fSize+1));
  fIndex[fSize]=ind; 
  fValue[fSize++]=val;
}

#endif
