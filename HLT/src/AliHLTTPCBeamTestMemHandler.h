// @(#) $Id$

#ifndef ALIL3TPCBEAMTESTMEMHANDLER_H
#define ALIL3TPCBEAMTESTMEMHANDLER_H

//_____________________________________________________________
// AliHLTTPCBeamTestMemHandler
//
// Class for converting the test beam data of May 2004 
// to the HLT file format using R. Bramms tables.
//
// Author: C. Loizides <loizides@ikf.uni-frankfurt.de>
// -- Copyright &copy ALICE HLT Group


#include "AliHLTMemHandler.h"

class AliHLTTPCBeamTestMemHandler : public AliHLTMemHandler{

 public:
  AliHLTTPCBeamTestMemHandler(Char_t *fPathToMappingFile="./MappingHWAdress.data");
  virtual ~AliHLTTPCBeamTestMemHandler();

  void SetInputPointer(Short_t *ptr, Int_t s) {fInputPtr=ptr;fInputSize=s;}
  void SetNTimeBins(Int_t i) {
   fNTimeBins=i;AliHLTTransform::SetNTimeBins(i);
  }
  void SetMinTimeBin(Int_t i) {
    fMinTimeBin=i;
  }

  Int_t MappingGetPadRow(Int_t channel) const {
    Int_t retval;
    if(channel < fNumOfChannels) retval = fMapping[channel][2]-fRowMin;
    else retval = -1;
    return retval;
  }

  Int_t MappingGetPad(Int_t channel) const {
    Int_t retval;
    if(channel < fNumOfChannels) retval = fMapping[channel][3];
    else retval = -1;
    return retval;
  }

  AliHLTDigitRowData* RawData2Memory(UInt_t &nrow,Int_t event=-1);
  Bool_t RawData2CompBinary(Int_t event=-1);

 private:

  struct AliRowStructure {
    Int_t fRow;       //row
    Int_t fNDigits;   //digits
    Int_t *fPadPos;   // pad position
  };
  
  Short_t **fMapping;//!          // mapping of channels to pads
  Short_t *fMappingEmptyRow;//!   // helper field
  Int_t fNumOfChannels;           // number of channels

  Short_t *fInputPtr;//!         // input pointer from shared memory
  Int_t fInputSize;//!           // size of input data
  AliRowStructure *fRows;//!     // rows
  Int_t fNTimeBins;              // number of timebins
  Int_t fMinTimeBin;             // min timebin (not zero because of altro)

  ClassDef(AliHLTTPCBeamTestMemHandler,1)   //RawData Filehandler class
};

typedef AliHLTTPCBeamTestMemHandler AliL3TPCBeamTestMemHandler; // for backward compatibility

#endif

