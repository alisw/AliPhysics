#ifndef ALIDNDPTBACKGROUNDCUTS_H
#define ALIDNDPTBACKGROUNDCUTS_H

//------------------------------------------------------------------------------
// Class to keep selection cuts for 
// cosmic and kinks and splitted tracks. 
// 
// Author: J.Otwinowski 25/01/2010
//------------------------------------------------------------------------------

class TParticle;
class AliESDtrack;
class AliExternalTrackParam;

#include "AliAnalysisCuts.h"

class AlidNdPtBackgroundCuts : public AliAnalysisCuts
{
public:
  AlidNdPtBackgroundCuts(const Char_t* name ="AlidNdPtBackgroundCuts", const Char_t *title ="");
  virtual ~AlidNdPtBackgroundCuts(); 
 
  // setters 
  void SetEtaWindow(const Float_t min=-10., const Float_t max=10.)  { fMinEta=min; fMaxEta=max; }
  void SetPhiWindow(const Float_t min=0., const Float_t max=1e99)  { fMinPhi=min; fMaxPhi=max;}
  void SetPtWindow(const Float_t min=0., const Float_t max=1e99)   { fMinPt=min;  fMaxPt=max;}
  void SetMaxFracSharedClust(const Float_t max=1.)   {fMaxFracSharedClust=max;}

  // getters 
  Float_t GetMinEta() const {return fMinEta;}
  Float_t GetMaxEta() const {return fMaxEta;}
  Float_t GetMinPhi() const {return fMinPhi;}
  Float_t GetMaxPhi() const {return fMaxPhi;}
  Float_t GetMinPt() const  {return fMinPt;}
  Float_t GetMaxPt() const  {return fMaxPt;}

  Float_t GetMaxFracSharedClust() const {return fMaxFracSharedClust;}

  // Get control histo
  THnSparseF *GetControlHisto() const {return fControlHisto;} 

  // cuts init function
  void Init();

  // check MC tracks
  virtual Bool_t IsSelected(TObject *) {return kTRUE;}
  virtual Bool_t IsSelected(TList *) {return kTRUE;}

  //
  Bool_t IsBackgroundTrack(AliESDtrack *track1, AliESDtrack *track2);
  Bool_t IsCosmicTrack(AliESDtrack *track1, AliESDtrack *track2);
  Bool_t IsSplittedTrack(AliESDtrack *track1, AliESDtrack *track2);
  
  // Merge output objects (needed by PROOF) 
  virtual Long64_t Merge(TCollection* list);

  // fill control histograms
  void  SetHistogramsOn(Bool_t fill=kTRUE) {fFillControlHisto = fill; }
  Bool_t  IsHistogramsOn() const {return fFillControlHisto; }

private:
  Float_t fMinEta; // min pseudorapidity limit
  Float_t fMaxEta; // max pseudorapidity limit
  Float_t fMinPhi; // min azimuthal angle (rad) limit
  Float_t fMaxPhi; // max azimuthal angle (rad) limit
  Float_t fMinPt;  // min pt limit
  Float_t fMaxPt;  // max pt limit
  Float_t fMaxFracSharedClust; // max fraction of track shared clusters 

  Bool_t  fFillControlHisto;  // flag to fill control histograms 
  THnSparseF *fControlHisto;  //-> etasum:dphi:dpt:pt1:fracSharedClust1:qsum

  AlidNdPtBackgroundCuts(const AlidNdPtBackgroundCuts&); // not implemented
  AlidNdPtBackgroundCuts& operator=(const AlidNdPtBackgroundCuts&); // not implemented

  ClassDef(AlidNdPtBackgroundCuts, 1)
};

#endif // 
