#ifndef ALIJETPARTICLESREADERESD_H
#define ALIJETPARTICLESREADERESD_H

//___________________________________________________________________________
/////////////////////////////////////////////////////////////////////////////
//                                                                         //
// Multi file reader for ESD                                               //
//                                                                         //
// This reader reads tracks from Event Summary Data                        //
// do not read particles                                                   //
// Piotr.Skowronski@cern.ch                                                //
// more info: http://alisoft.cern.ch/people/skowron/analyzer/index.html    //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////

#include <TString.h>
#include <AliESDtrack.h>
#include "AliJetParticlesReader.h"

class TFile;
class AliESD;

class AliJetParticlesReaderESD: public AliJetParticlesReader
{
  public:
  AliJetParticlesReaderESD(const Char_t* esdfilename = "AliESDs.root") ;
  AliJetParticlesReaderESD(TObjArray* dirs,const Char_t* esdfilename = "AliESDs.root");

  void SetCompareFlag(ULong_t f){fPassFlag=f;}
  void SetCompareFlagTPC()
    {fPassFlag=AliESDtrack::kTPCrefit & AliESDtrack::kTPCrefit;}

#if 0
    kITSin=0x0001,kITSout=0x0002,kITSrefit=0x0004,kITSpid=0x0008,
    kTPCin=0x0010,kTPCout=0x0020,kTPCrefit=0x0040,kTPCpid=0x0080,
    kTRDin=0x0100,kTRDout=0x0200,kTRDrefit=0x0400,kTRDpid=0x0800,
    kTOFin=0x1000,kTOFout=0x2000,kTOFrefit=0x4000,kTOFpid=0x8000,
    kPHOSpid=0x10000, kRICHpid=0x20000,
    kTRDStop=0x20000000,
    kESDpid=0x40000000,
    kTIME=0x80000000
#endif

  virtual ~AliJetParticlesReaderESD();

  //Int_t Next();
  void Rewind();

  protected:
  Int_t    ReadESD(AliESD* esd); //read esd file/objects
  Int_t    ReadNext();           //read the next event
  TFile*   OpenFile(Int_t evno); //opens file to be read for given event
  Bool_t   IsAcceptedParticle(Float_t px, Float_t py, Float_t pz) const;
    
  TString fESDFileName; // name of the file with tracks
  TFile*  fFile;        //! pointer to current ESD file
  TIter*  fKeyIterator; //! key iterator through file
  ULong_t fPassFlag;    //flag to compare esd flag with 

  ClassDef(AliJetParticlesReaderESD,1) //
};

inline Bool_t AliJetParticlesReaderESD::IsAcceptedParticle(Float_t pt, Float_t phi, Float_t eta) const
{
  if((pt<fPtMin)||(pt>fPtMax)) return kFALSE;
  if((eta<fEtaMin)||(eta>fEtaMax)) return kFALSE;
  if((phi<fPhiMin)||(phi>fPhiMax)) return kFALSE;

  return kTRUE;
}
#endif
