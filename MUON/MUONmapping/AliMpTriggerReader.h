/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$
// $MpId: AliMpTriggerReader.h,v 1.5 2006/05/24 13:58:27 ivana Exp $

/// \ingroup mptrigger
/// \class AliMpTriggerReader
/// \brief Read trigger slat ASCII files
///
//  Author: Laurent Aphecetche

#ifndef ALI_MP_TRIGGER_READER_H
#define ALI_MP_TRIGGER_READER_H

#ifndef ROOT_TObject
#  include "TObject.h"
#endif

#ifndef ROOT_TMap
#  include "TMap.h"
#endif

#ifndef ROOT_TString
#  include "TString.h"
#endif

#ifndef ALI_MP_PLANE_TYPE_H
#  include "AliMpPlaneType.h"
#endif

class AliMpSlatMotifMap;
class AliMpSlat;
class AliMpTrigger;
class AliMpPCB;
class AliMpDataStreams;

class TList;

class AliMpTriggerReader : public TObject
{
 public:
  AliMpTriggerReader(AliMpSlatMotifMap* motifMap);
  virtual ~AliMpTriggerReader();

  AliMpTrigger* ReadSlat(const AliMpDataStreams&  dataStreams,
                         const char* slatType, AliMp::PlaneType planeType);

  AliMpPCB* ReadPCB(const AliMpDataStreams&  dataStreams, const char* pcbType);
  
private:
    
  AliMpSlat* BuildSlat(const AliMpDataStreams&  dataStreams,
                              const char* slatName,
                              AliMp::PlaneType planeType,
                              const TList& descriptionLines,
                              Double_t scale=1.0);
  
  Int_t DecodeFlipLine(const TString& sline,
                              TString& slatType2,
                              Bool_t& flipX, Bool_t& flipY);
  
  Int_t DecodeScaleLine(const TString& sline, 
                               Double_t& scale, TString& slatType);
  
  void FlipLines(const AliMpDataStreams&  dataStreams,
                        TList& lines, Bool_t flipX, Bool_t flipY,
                        Int_t srcLine, Int_t destLine);
  
  TString GetBoardNameFromPCBLine(const TString& sline);
  
  Int_t GetLine(const TString& slatType);
  
  Int_t IsLayerLine(const TString& sline) const;
    
  int LocalBoardNumber(const AliMpDataStreams&  dataStreams,
                       const char* localBoardName);
  
  // AliMpPCB* PCB(const char* pcbType); 
  
  void ReadLines(const AliMpDataStreams&  dataStreams,
                        const char* slatType,
                        AliMp::PlaneType planeType,
                        TList& lines,
                        Double_t& scale, Bool_t& flipX, Bool_t& flipY,
                        Int_t& srcLine, Int_t& destLine);
  
  void ReadLocalBoardMapping(const AliMpDataStreams&  dataStreams);
  
private:
  /// Not implemented
  AliMpTriggerReader(const AliMpTriggerReader& rhs);
  /// Not implemented
  AliMpTriggerReader& operator=(const AliMpTriggerReader& rhs);
    
  // static methods
  static const TString& GetKeywordLayer();
  static const TString& GetKeywordScale();
  static const TString& GetKeywordPcb();
  static const TString& GetKeywordFlipX();
  static const TString& GetKeywordFlipY();
  
  // data members
  AliMpSlatMotifMap* fMotifMap; //!<! storage for motifTypes and motifs...
  
  TMap fLocalBoardMap; //!<! map of TObjString to TObjString

  ClassDef(AliMpTriggerReader,0) // Reader for trigger slats mapping files 
};

#endif
