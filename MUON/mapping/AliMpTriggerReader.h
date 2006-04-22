/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$
// $MpId: AliMpTriggerReader.h,v 1.2 2006/03/02 16:36:26 ivana Exp $

/// \ingroup trigger
/// \class AliMpTriggerReader
/// \brief Read trigger slat ASCII files
/// \author Laurent Aphecetche

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

class AliMpSlat;
class AliMpTrigger;
class AliMpPCB;
class TList;

class AliMpTriggerReader : public TObject
{
 public:
  AliMpTriggerReader();
  virtual ~AliMpTriggerReader();

  static AliMpTrigger* ReadSlat(const char* slatType, AliMpPlaneType planeType);

  static AliMpPCB* ReadPCB(const char* pcbType);
  
  static void Reset();
  
private:
    
  static AliMpSlat* BuildSlat(const char* slatName, 
                              AliMpPlaneType planeType,
                              const TList& descriptionLines,
                              Double_t scale=1.0);
  
  static Int_t DecodeFlipLine(const TString& sline,
                              TString& slatType2,
                              Bool_t& flipX, Bool_t& flipY);
  
  static Int_t DecodeScaleLine(const TString& sline, 
                               Double_t& scale, TString& slatType);
  
  static void FlipLines(TList& lines, Bool_t flipX, Bool_t flipY, 
                        Int_t srcLine, Int_t destLine);
  
  static TString GetBoardNameFromPCBLine(const TString& sline);
  
  static Int_t GetLine(const TString& slatType);
  
  static Int_t IsLayerLine(const TString& sline);
    
  static int LocalBoardNumber(const char* localBoardName);
  
  static AliMpPCB* PCB(const char* pcbType); 
  
  static void ReadLines(const char* slatType,
                        AliMpPlaneType planeType,
                        TList& lines,
                        Double_t& scale, Bool_t& flipX, Bool_t& flipY,
                        Int_t& srcLine, Int_t& destLine);
  
  static void ReadLocalBoardMapping();
  
private:
    
  static TMap fgPCBMap; //! map of TObjString to AliMpPCB*
  
  static TMap fgLocalBoardMap; //! map of TObjString to TObjString

  static const TString fgkKeywordLayer; //! Keyword: LAYER
  static const TString fgkKeywordScale; //! Keyword: SCALE
  static const TString fgkKeywordPcb; //! Keyword : PCB
  static const TString fgkKeywordFlipX; //! Keyword : FLIPX
  static const TString fgkKeywordFlipY; //! Keyword : FLIPY
  
  ClassDef(AliMpTriggerReader,1) // Reader for trigger slats mapping files 
};

#endif
