#ifndef ALIJETPARTICLESREADERESD_H
#define ALIJETPARTICLESREADERESD_H

//___________________________________________________________________________
/////////////////////////////////////////////////////////////////////////////
//                                                                         //
// File reader for ESD                                                     //
//                                                                         //
// This reader reads tracks from Event Summary Data                        //
// taken from Piotr.Skowronski@cern.ch                                     //
// more info: http://alisoft.cern.ch/people/skowron/analyzer/index.html    //
//                                                                         //
// loizides@ikf.uni-frankfurt.de                                           //
/////////////////////////////////////////////////////////////////////////////

#include <TString.h>
#include <AliESDtrack.h>
#include "AliJetParticlesReader.h"

class TFile;
class TTree;
class AliESD;

class AliJetParticlesReaderESD: public AliJetParticlesReader
{
  public:
  AliJetParticlesReaderESD(const Char_t* esdfilename = "AliESDs.root") ;
  AliJetParticlesReaderESD(TObjArray* dirs,const Char_t* esdfilename = "AliESDs.root");

  void SetCompareFlag(ULong_t f){fPassFlag=f;}
  void SetCompareFlagTPC() {fPassFlag=AliESDtrack::kTPCrefit;}

  virtual ~AliJetParticlesReaderESD();

  //Int_t Next(); //in base class
  void Rewind();

  const AliESD* GetCurrentESD() const {return fESD;}

  protected:
  virtual Int_t ReadESD(AliESD* esd); //read esd file/objects
  Int_t    ReadNext();                //read the next event
  TFile*   OpenFile(Int_t evno);      //opens file to be read for given event
  Bool_t   IsAcceptedParticle(Float_t px, Float_t py, Float_t pz) const;
    
  TString fESDFileName; // name of the file with tracks
  AliESD *fESD;         //! pointer to current esd object
  TFile*  fFile;        //! pointer to current ESD file
  TTree*  fTree;        //! pointer to current tree with ESD objects
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
