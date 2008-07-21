#ifndef ALIEMCALSTURAWSTREAM_H
#define ALIEMCALSTURAWSTREAM_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

///////////////////////////////////////////////////////////////////////////////
///
/// This is a class for reading the EMCAL STU (trigger) DDL raw data
///
///  Two different formats foreseen (different by their length and Format Version):
///
///  - The regular one, containing jet patch indices (7x 32 bit words) and
///  gamma jet patch indices (~96 combinations per region x 32 regions / 32
///  bit per words = 96 words).
///
///  - The debug one containing the previous informations (103 words) plus
///  all the time integrated 2x2 sums used for generating the accepted
///  trigger (96 x 32 x 16 bit / 32 = 1536 words).
///
///////////////////////////////////////////////////////////////////////////////

#include <TObject.h>

class AliRawReader;

class AliEMCALSTURawStream: public TObject {

  public :

    enum { kNumJetPatchWords = 7, // max number of jet patch indices words
	   kNumGammaJetPatchWords = 96, // max number of gamma-jet patch indices words
	   kNum2x2Words = 1536,     //  max number of 2x2 indices words
	   kEMCALSTUDDL = 0x12c0 }; // 4800, in EMCAL DDL range but separated from regular FEE

    AliEMCALSTURawStream(AliRawReader* rawReader);
    virtual ~AliEMCALSTURawStream();
  
    virtual void             Reset();
    virtual Bool_t           Next();

    UInt_t    GetJetPatchWords(int i) const { return fJetPatchWords[i]; }
    UInt_t    GetGammaJetPatchWords(int i) const { return fGammaJetPatchWords[i]; }
    UInt_t    Get2x2Words(int i) const { return f2x2Words[i]; }
    UInt_t    GetNum2x2Words() const { return fNum2x2Words; }

  protected:
    AliEMCALSTURawStream(const AliEMCALSTURawStream& stream);
    AliEMCALSTURawStream& operator = (const AliEMCALSTURawStream& stream);

  private:

    AliRawReader*    fRawReader;   // object for reading the raw data
    UInt_t           fJetPatchWords[kNumJetPatchWords]; // jet patch indices
    UInt_t           fGammaJetPatchWords[kNumGammaJetPatchWords]; // gamma jet patch indices
    UInt_t           f2x2Words[kNum2x2Words]; // 2x2 sums
    UInt_t           fNum2x2Words; // how many 2x2 sums did we actually read?

    ClassDef(AliEMCALSTURawStream, 0)   // class for reading EMCAL STU DDL raw data
};

#endif
