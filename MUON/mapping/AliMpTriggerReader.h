/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$
// $MpId: AliMpTriggerReader.h,v 1.5 2006/05/24 13:58:27 ivana Exp $

/// \ingroup trigger
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

#ifndef ALI_MP_STATION_TYPE_H
#  include "AliMpStationType.h"
#endif

class AliMpSlatMotifMap;
class AliMpSlat;
class AliMpTrigger;
class AliMpPCB;
class TList;

class AliMpTriggerReader : public TObject
{
 public:
  AliMpTriggerReader(AliMpSlatMotifMap& motifMap);
  virtual ~AliMpTriggerReader();

  AliMpTrigger* ReadSlat(const char* slatType, AliMpPlaneType planeType);

  AliMpPCB* ReadPCB(const char* pcbType);
  
private:
    
  AliMpSlat* BuildSlat(const char* slatName, 
                              AliMpPlaneType planeType,
                              const TList& descriptionLines,
                              Double_t scale=1.0);
  
  Int_t DecodeFlipLine(const TString& sline,
                              TString& slatType2,
                              Bool_t& flipX, Bool_t& flipY);
  
  Int_t DecodeScaleLine(const TString& sline, 
                               Double_t& scale, TString& slatType);
  
  void FlipLines(TList& lines, Bool_t flipX, Bool_t flipY, 
                        Int_t srcLine, Int_t destLine);
  
  TString GetBoardNameFromPCBLine(const TString& sline);
  
  Int_t GetLine(const TString& slatType);
  
  Int_t IsLayerLine(const TString& sline);
    
  int LocalBoardNumber(const char* localBoardName);
  
  AliMpPCB* PCB(const char* pcbType); 
  
  void ReadLines(const char* slatType,
                        AliMpPlaneType planeType,
                        TList& lines,
                        Double_t& scale, Bool_t& flipX, Bool_t& flipY,
                        Int_t& srcLine, Int_t& destLine);
  
  void ReadLocalBoardMapping();
  
private:
    
  AliMpSlatMotifMap& fMotifMap; //!< storage for motifTypes and motifs...
  
  TMap fLocalBoardMap; //!< map of TObjString to TObjString

  static const TString fgkKeywordLayer; //!< Keyword: LAYER
  static const TString fgkKeywordScale; //!< Keyword: SCALE
  static const TString fgkKeywordPcb; //!< Keyword : PCB
  static const TString fgkKeywordFlipX; //!< Keyword : FLIPX
  static const TString fgkKeywordFlipY; //!< Keyword : FLIPY
  
  ClassDef(AliMpTriggerReader,0) // Reader for trigger slats mapping files 
};

#endif
