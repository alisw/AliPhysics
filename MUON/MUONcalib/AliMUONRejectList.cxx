/**************************************************************************
* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
*                                                                        *
* Author: The ALICE Off-line Project.                                    *
* Contributors are mentioned in the code where appropriate.              *
*                                                                        *
* Permission to use, copy, modify and distribute this software and its   *
* documentation strictly for non-commercial purposes is hereby granted   *
* without fee, provided that the above copyright notice appears in all   *
* copies and that both the copyright notice and this permission notice   *
* appear in the supporting documentation. The authors make no claims     *
* about the suitability of this software for any purpose. It is          *
* provided "as is" without express or implied warranty.                  *
**************************************************************************/

// $Id$

///
/// \class AliMUONRejectList
///
/// Object to hold the probability to reject elements during reconstruction.
///
/// Those elements are either channels, manus, 
/// bus patches, detection elements, pcbs (for slats) and HV channels,
/// or all of them.
///
/// (we do not consider the next level, chamber, because if a full
/// chamber is missing, we assume we'll remove that run from the
/// list of usable runs anyway).
/// 
/// The probability of rejection can be different from 0.0 or 1.0 only for
/// simulations, in which case it means that :
/// - this object was created from inspection of real data occupancy maps + ad-hoc rejections over
/// several runs
/// - the probability then reflects the chance a given element was useable
/// during data taking. For instance, if one DE has a probability of 0.8, it means
/// it was on (or correctly behaving during data taking) only for 20% of the
/// events.
///
/// \author Laurent Aphecetche, Subatech
///

#include "AliMUONRejectList.h"

#include "AliLog.h"
#include "AliMpArea.h"
#include "AliMpConstants.h"
#include "AliMpDCSNamer.h"
#include "AliMpDDLStore.h"
#include "AliMpDEStore.h"
#include "AliMpDetElement.h"
#include "AliMpMotifPosition.h"
#include "AliMpPCB.h"
#include "AliMpSegmentation.h"
#include "AliMpSlat.h"
#include "AliMpVSegmentation.h"
#include "AliMUON2DMap.h"
#include "AliMUONCalibParamNF.h"
#include "Riostream.h"
#include "TMath.h"

using std::cout;
using std::endl;
/// \cond CLASSIMP
ClassImp(AliMUONRejectList)
/// \endcond

namespace
{
  /// The functions below are here to help re-invent the wheel,
  /// i.e. code something which acts as std::map<int,float>
  /// to circumvent the fact that AliRoot does not allow STL to be used.
  
  void Dump(const char* str, UInt_t n, UInt_t* ids, Float_t* values, Bool_t debug)
  {
    /// Dump the values array
    
    TString s(str);
    s += " PROBA %e";

    for ( UInt_t i = 0; i < n; ++i )
    {
      UInt_t key = ids[i];
      Int_t a,b;
      AliMUONVCalibParam::DecodeUniqueID(key,a,b);
      if ( s.CountChar('%')==3 )
      {
        cout << Form(s.Data(),a,b,values[i]) << endl;
      }
      else
      {
        cout << Form(s.Data(),a,values[i]) << endl;
      }
    }
    
    if ( debug ) 
    {
      cout << "------" << endl;
      for ( UInt_t i = 0; i < n; ++i )
      {
        UInt_t key = ids[i];
        Int_t a,b;
        AliMUONVCalibParam::DecodeUniqueID(key,a,b);
        cout << Form("ids[%5d]=%d values[%5d]=%e (a,b)=(%5d,%5d)",
        i,ids[i],i,values[i],a,b) << endl;
      }
      
    }
  }
  
  Float_t GetValue(UInt_t n, UInt_t* ids, Float_t* values, UInt_t key)
  {
    /// Get the value corresponding to key, or zero if not found
    Long64_t index = TMath::BinarySearch(n,ids,key);

    Bool_t found = ( ( index >= 0 ) && ( ids[index] == key ) );

    if ( found )
    {
      return values[index];
    }
    else
    {
      return 0.0;
    }
  }
  
  void Insert(UInt_t n, UInt_t* ids, Float_t* values, UInt_t index, UInt_t key, Float_t value)
  {
    /// Insert (key,value) into arrays ids and values.
    
    for ( UInt_t i = n; i > index; --i )
    {
      ids[i] = ids[i-1];
      values[i] = values[i-1];
    }
    ids[index] = key;
    values[index] = value;
  }
  
  Bool_t SetValue(UInt_t n, UInt_t* ids, Float_t* values, UInt_t key, Float_t value)
  {
    /// Set the value for a given key
    Long64_t index = TMath::BinarySearch(n,ids, key);

    Bool_t alreadyThere = ( ( index >= 0 ) && ( ids[index] == key ) );

    if ( alreadyThere )
    {
      // replacement
      values[index] = value;
      return kFALSE;
    }
    Insert(n,ids,values,index+1,key,value);
    return kTRUE;
  }
  
  void Copy(UInt_t n, UInt_t* src, UInt_t*& dest)
  {
    /// Copy src into dest
    delete[] dest;
    dest = 0;
    if ( src && n )
    {
      dest = new UInt_t[n];
      memcpy(dest,src,n*sizeof(UInt_t));
    }
  }

  void Copy(UInt_t n, Float_t* src, Float_t*& dest)
  {
    /// Copy src into dest
    delete[] dest;
    dest = 0;
    if ( src && n )
    {
      dest = new Float_t[n];
      memcpy(dest,src,n*sizeof(Float_t));
    }
  }
  
}
    
//_____________________________________________________________________________
AliMUONRejectList::AliMUONRejectList()
: TObject(),
fIsBinary(kTRUE),
fMaxNofDEs(156), // not nice to put a constant here, but that way this object does not need the mapping at creation time...
fMaxNofBPs(888), // same remark as above
fMaxNofManus(16828), // same as above...
fNofDEs(),
fNofBPs(),
fNofManus(),
fDEIds(new UInt_t[fMaxNofDEs]),
fDEProbas(new Float_t[fMaxNofDEs]),
fBPIds(new UInt_t[fMaxNofBPs]),
fBPProbas(new Float_t[fMaxNofBPs]),
fManuIds(new UInt_t[fMaxNofManus]),
fManuProbas(new Float_t[fMaxNofManus]),
fChannels(new AliMUON2DMap(kTRUE))
{
  /// normal ctor
  memset(fDEIds,0,fMaxNofDEs*sizeof(UInt_t));
  memset(fDEProbas,0,fMaxNofDEs*sizeof(Float_t));
  memset(fBPIds,0,fMaxNofBPs*sizeof(UInt_t));
  memset(fBPProbas,0,fMaxNofBPs*sizeof(Float_t));
  memset(fManuIds,0,fMaxNofManus*sizeof(UInt_t));
  memset(fManuProbas,0,fMaxNofManus*sizeof(Float_t));
}

//_____________________________________________________________________________
AliMUONRejectList::AliMUONRejectList(TRootIOCtor* /*ioCtor*/)
: TObject(),
fIsBinary(kTRUE),
fMaxNofDEs(),
fMaxNofBPs(),
fMaxNofManus(),
fNofDEs(), 
fNofBPs(), 
fNofManus(0), 
fDEIds(0x0),
fDEProbas(0x0),
fBPIds(0x0),
fBPProbas(0x0),
fManuIds(0x0),
fManuProbas(0x0),
fChannels(0x0)
{
  /// ctor from root i/o
}

//_____________________________________________________________________________
AliMUONRejectList::AliMUONRejectList(const AliMUONRejectList& rl)
: TObject(rl),
fIsBinary(rl.fIsBinary),
fMaxNofDEs(rl.fMaxNofDEs),
fMaxNofBPs(rl.fMaxNofBPs),
fMaxNofManus(rl.fMaxNofManus),
fNofDEs(rl.fNofDEs), 
fNofBPs(rl.fNofBPs), 
fNofManus(rl.fNofManus), 
fDEIds(0x0),
fDEProbas(0x0),
fBPIds(0x0),
fBPProbas(0x0),
fManuIds(0x0),
fManuProbas(0x0),
fChannels(0x0)
{
  /// Copy ctor
  
  ::Copy(rl.fMaxNofDEs,rl.fDEIds,fDEIds);
  ::Copy(rl.fMaxNofDEs,rl.fDEProbas,fDEProbas);
  ::Copy(rl.fMaxNofBPs,rl.fBPIds,fBPIds);
  ::Copy(rl.fMaxNofBPs,rl.fBPProbas,fBPProbas);
  ::Copy(rl.fMaxNofManus,rl.fManuIds,fManuIds);
  ::Copy(rl.fMaxNofManus,rl.fManuProbas,fManuProbas);
  
  if ( rl.fChannels ) 
  {
    fChannels = static_cast<AliMUONVStore*>(rl.fChannels->Clone());
  }
}

//_____________________________________________________________________________
AliMUONRejectList& AliMUONRejectList::operator=(const AliMUONRejectList& rl)
{
  /// assignement operator
  if ( this != &rl ) 
  {
    static_cast<TObject&>(*this)=rl;
    
    fIsBinary = rl.fIsBinary;
    fMaxNofDEs = rl.fMaxNofDEs;
    fMaxNofBPs = rl.fMaxNofBPs;
    fMaxNofManus = rl.fMaxNofManus;
    fNofDEs = rl.fNofDEs;
    fNofBPs = rl.fNofBPs;
    fNofManus = rl.fNofManus;
    
    ::Copy(rl.fMaxNofDEs,rl.fDEIds,fDEIds);
    ::Copy(rl.fMaxNofDEs,rl.fDEProbas,fDEProbas);
    ::Copy(rl.fMaxNofBPs,rl.fBPIds,fBPIds);
    ::Copy(rl.fMaxNofBPs,rl.fBPProbas,fBPProbas);
    ::Copy(rl.fMaxNofManus,rl.fManuIds,fManuIds);
    ::Copy(rl.fMaxNofManus,rl.fManuProbas,fManuProbas);
    
    delete fChannels;
    fChannels = 0x0;
    
    if ( rl.fChannels ) 
    {
      fChannels = static_cast<AliMUONVStore*>(rl.fChannels->Clone());
    }
    
  }  
  return *this;
}

//_____________________________________________________________________________
AliMUONRejectList::~AliMUONRejectList()
{
  /// dtor
  delete fChannels;
  delete[] fDEIds;
  delete[] fDEProbas;
  delete[] fBPIds;
  delete[] fBPProbas;
  delete[] fManuIds;
  delete[] fManuProbas;
}

//_____________________________________________________________________________
Float_t AliMUONRejectList::DetectionElementProbability(Int_t detElemId) const
{
  /// Get the probability to reject a given detection element
  return ::GetValue(fNofDEs,fDEIds,fDEProbas,AliMUONVCalibParam::BuildUniqueID(detElemId,0));
}

//_____________________________________________________________________________
Float_t AliMUONRejectList::BusPatchProbability(Int_t busPatchId) const
{
  /// Get the probability to reject a given bus patch
  return ::GetValue(fNofBPs,fBPIds,fBPProbas,AliMUONVCalibParam::BuildUniqueID(busPatchId,0));
}

//_____________________________________________________________________________
Float_t AliMUONRejectList::ManuProbability(Int_t detElemId, Int_t manuId) const
{
  /// Get the probability to reject a given manu
  return ::GetValue(fNofManus,fManuIds,fManuProbas,AliMUONVCalibParam::BuildUniqueID(detElemId,manuId));
}

//_____________________________________________________________________________
Float_t AliMUONRejectList::ChannelProbability(Int_t detElemId, Int_t manuId, Int_t manuChannel) const
{
  /// Get the probability to reject a given channel
  AliMUONVCalibParam* param = static_cast<AliMUONVCalibParam*>(fChannels->FindObject(detElemId,manuId));
  if (!param) return 0.0;
  return param->ValueAsFloat(manuChannel);
}

//_____________________________________________________________________________
void AliMUONRejectList::SetDetectionElementProbability(Int_t detElemId, Float_t proba)
{
  /// Set the probability to reject a given detection element
  if ( ::SetValue(fNofDEs,fDEIds,fDEProbas,AliMUONVCalibParam::BuildUniqueID(detElemId,0),proba) ) 
  {
    ++fNofDEs;
  }
  
  ZeroOrOne(proba);
}

//_____________________________________________________________________________
void AliMUONRejectList::ZeroOrOne(Float_t proba)
{
  /// If proba is anything else than 0 or 1, we set fIsBinary to kFALSe
  
  Bool_t zeroorone = ( proba == 0.0 || proba == 1.0 );
  if (!zeroorone) fIsBinary = kFALSE;
}

//_____________________________________________________________________________
void AliMUONRejectList::SetBusPatchProbability(Int_t busPatchId, Float_t proba)
{
  /// Set the probability to reject a given bus patch
  if ( ::SetValue(fNofBPs,fBPIds,fBPProbas,AliMUONVCalibParam::BuildUniqueID(busPatchId,0),proba) )
  {
    ++fNofBPs;
  }
  ZeroOrOne(proba);
}

//_____________________________________________________________________________
void AliMUONRejectList::SetManuProbability(Int_t detElemId, Int_t manuId, Float_t proba)
{
  /// Set the probability to reject a given manu
  if ( ::SetValue(fNofManus,fManuIds,fManuProbas,AliMUONVCalibParam::BuildUniqueID(detElemId,manuId),proba) )
  {
    ++fNofManus;
  }
  ZeroOrOne(proba);
}

//_____________________________________________________________________________
void AliMUONRejectList::SetChannelProbability(Int_t detElemId, Int_t manuId, Int_t manuChannel, Float_t proba)
{
  /// Set the probability to reject a given channel
  AliMUONVCalibParam* param = static_cast<AliMUONVCalibParam*>(fChannels->FindObject(detElemId,manuId));
  if (!param) 
  {
    param = new AliMUONCalibParamNF(1,AliMpConstants::ManuNofChannels(),detElemId,manuId,0.0);    
    fChannels->Add(param);
  }
  param->SetValueAsFloat(manuChannel,0,proba);
  ZeroOrOne(proba);
}

//_____________________________________________________________________________
void AliMUONRejectList::SetPCBProbability(Int_t detElemId, Int_t pcbNumber, Float_t proba)
{
  /// Set the probability to reject all the manus of a given (slat) PCB
  AliMpSegmentation* seg = AliMpSegmentation::Instance();
  AliMp::CathodType ct[] = { AliMp::kCath0, AliMp::kCath1 };
  
  for ( Int_t i = 0; i < 2; ++i )
  {
    const AliMpVSegmentation* vseg = seg->GetMpSegmentation(detElemId,ct[i]);
    if (!vseg)
    {
      AliError(Form("Could not get segmentation of DE %d",detElemId));
      continue;
    }
    const AliMpSlat* slat = seg->GetSlat(vseg);
    if (!slat)
    {
      AliError(Form("Could not get slat from DE %d",detElemId));
      continue;      
    }
    AliMpPCB* pcb = slat->GetPCB(pcbNumber);
    for ( Int_t j = 0; j < pcb->GetSize(); ++j )
    {
      AliMpMotifPosition* mp = pcb->GetMotifPosition(j);
      SetManuProbability(detElemId,mp->GetID(),proba);
    }
  }
}

//_____________________________________________________________________________
void AliMUONRejectList::SetHVProbability(const char* dcsName, Float_t proba)
{
  /// Set the probability to reject all the manus of a given HV part
  /// Caution : the dcs string is a dcs NAME, _not_ an alias
  
  AliMpDCSNamer hv("TRACKER");

  TString alias = hv.DCSAliasFromName(dcsName);

  Int_t detElemId = hv.DetElemIdFromDCSAlias(alias.Data());
  Int_t index = hv.DCSIndexFromDCSAlias(alias.Data());

  AliMpDetElement* de = AliMpDDLStore::Instance()->GetDetElement(detElemId);

  const AliMpArrayI* manus = de->ManusForHV(index);

  for ( Int_t i = 0; i < manus->GetSize(); ++ i )
  {
    Int_t manuId = manus->GetValue(i);
    SetManuProbability(detElemId,manuId,proba);
  }
}

//_____________________________________________________________________________
void 
AliMUONRejectList::Print(Option_t* opt) const
{
  /// Printout

  TString sopt(opt);  
  sopt.ToUpper();
  Bool_t debug(kFALSE);
  
  if ( sopt.Contains("DEBUG") ) debug=kTRUE;
  
  cout << Form("We have probabilities for %d detection element(s), %d bus patch(es), %d manu(s)",
               fNofDEs,fNofBPs,fNofManus) << endl;
  
  ::Dump("DE %04d",fNofDEs,fDEIds,fDEProbas,debug);
  ::Dump("BusPatch %04d",fNofBPs,fBPIds,fBPProbas,debug);
  ::Dump("DE %04d MANU %4d",fNofManus,fManuIds,fManuProbas,debug);
}
