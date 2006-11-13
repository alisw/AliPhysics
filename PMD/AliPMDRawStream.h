#ifndef ALIPMDRAWSTREAM_H
#define ALIPMDRAWSTREAM_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
///
/// This class provides access to PMD digits in raw data.
///
///////////////////////////////////////////////////////////////////////////////

#include <TObject.h>

class AliRawReader;


class AliPMDRawStream: public TObject {
  public :
    AliPMDRawStream(AliRawReader* rawReader);
    virtual ~AliPMDRawStream();

    Bool_t DdlData(Int_t indexDDL, TObjArray *pmdddlcont);

    enum {kDDLOffset = 0xC00};      // offset for DDL numbers

  private :
    AliPMDRawStream(const AliPMDRawStream& stream);
    AliPMDRawStream& operator = (const AliPMDRawStream& stream);

    void             GetRowCol(Int_t ddlno, Int_t pbusid, 
			       UInt_t mcmno, UInt_t chno,
			       Int_t startRowBus[], Int_t endRowBus[],
			       Int_t startColBus[], Int_t endColBus[],
			       Int_t &row, Int_t &col) const;
    void             ConvertDDL2SMN(Int_t iddl, Int_t imodule,
				    Int_t &smn, Int_t &detector) const;
    void             TransformH2S(Int_t smn, Int_t &row, Int_t &col) const;

    AliRawReader*    fRawReader;    // object for reading the raw data


    ClassDef(AliPMDRawStream, 3)    // class for reading PMD raw digits
};

#endif
