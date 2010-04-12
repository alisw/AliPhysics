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
//
// Utility class for V0 PID
// Provides smaples of electrons, pions and protons
// More information can be found in the implementation file
//
#ifndef ALIHFEV0PID_H
#define ALIHFEV0PID_H

#ifndef ROOT_TObject
#include <TObject.h>
#endif

#ifndef ALIHFECOLLECTION_H
#include "AliHFEcollection.h"
#endif

class TObjArray;
class TList;

class AliESDv0;
class AliESDtrack;
class AliKFParticle;
class AliKFVertex;
class AliVEvent;
class AliVTrack;

class AliHFEV0pid : public TObject{
  public:
  enum{ // Reconstructed V0
    kUndef = 0,
      kRecoGamma = 1,
      kRecoK0s = 2,
      kRecoPhi = 3,
      kRecoLambda = 4
      
    };
    enum{ // Identified Daughter particles
      kRecoElectron = 0,
	    kRecoPionK0 = 1,
	    kRecoPionL = 2,
	    kRecoKaon = 3,
	    kRecoProton = 4
    };
    AliHFEV0pid();
    ~AliHFEV0pid();

    Double_t  GetEffMass(AliESDv0 *v0, UInt_t p1, UInt_t p2) const;

    void Process(AliVEvent * const inputEvent);
    Int_t ProcessV0(TObject *v0);
    void Flush();

    void InitQA();
    inline TList *GetListOfQAhistograms();

    TObjArray *GetListOfElectrons() const { return fElectrons; }
    TObjArray *GetListOfPionsK0() const { return fPionsK0; }
    TObjArray *GetListOfPionsL() const { return fPionsL; }
    TObjArray *GetListOfKaons() const { return fKaons; }
    TObjArray *GetListOfProtons() const { return fProtons; }

    Bool_t IsAODanalysis() const { return TestBit(kAODanalysis); }
    Bool_t IsESDanalysis() const { return !TestBit(kAODanalysis); }
    void SetAODanalysis(Bool_t isAOD = kTRUE) { SetBit(kAODanalysis, isAOD); };
    void SetESDanalysis(Bool_t isESD = kTRUE) { SetBit(kAODanalysis, !isESD); }; 
 private:
    Float_t PsiPair(AliESDv0 *esdv0);//angle between daughters in plane perpendicular to magnetic field (characteristically around zero for conversions)
    Float_t OpenAngle(AliESDv0 *esdv0) const;//opening angle between V0 daughters; close to zero for conversions
   

  protected:
    enum{
      kAODanalysis = BIT(14)
    };
    AliKFParticle *CreateMotherParticle(AliVTrack *pdaughter, AliVTrack *ndaughter, Int_t pspec, Int_t nspec);
    void AddTrackToKFVertex(AliVTrack *track, Int_t species);

    Bool_t IsGammaConv(TObject *v0);
    Bool_t IsK0s(TObject *v0);
    Bool_t IsPhi(TObject *v0);
    Bool_t IsLambda(TObject *v0);

    Bool_t CutESDtrack(AliESDtrack *track);
    Bool_t CutV0(AliESDv0 *v0, Int_t species);

    Bool_t LooseRejectK0(AliESDv0 * const v0) const;
    Bool_t LooseRejectLambda(AliESDv0 * const v0) const;
    Bool_t LooseRejectGamma(AliESDv0 * const v0) const;

  private:
    class AliHFEV0pidTrackIndex{
      public:
        AliHFEV0pidTrackIndex();
        ~AliHFEV0pidTrackIndex();
        void Init(Int_t capacity);
        void Add(Int_t index, Int_t species);
        Bool_t Find(Int_t index) const;
        Bool_t Find(Int_t index, Int_t species) const;
        Int_t GetNumberOfElectrons() const { return fNElectrons; };
        Int_t GetNumberOfPionsK0() const { return fNPionsK0; };
        Int_t GetNumberOfPionsL() const { return fNPionsL; };
	      Int_t GetNumberOfKaons() const { return fNKaons; };
        Int_t GetNumberOfProtons() const { return fNProtons; };
        void Flush();

      private:
        AliHFEV0pidTrackIndex(const AliHFEV0pidTrackIndex &ref);
        AliHFEV0pidTrackIndex &operator=(const AliHFEV0pidTrackIndex &ref);
        Int_t fNElectrons;        // Number of identified electrons
        Int_t fNPionsK0;          // Number of identified pions from K0s
        Int_t fNPionsL;           // Lumber of identified pions from Lambda
	      Int_t fNKaons;            // Number of identified kaons
        Int_t fNProtons;          // Number of identified protons
        Int_t *fIndexElectron;    // Indices of identified electrons
        Int_t *fIndexPionK0;      // Indices of identified pions from K0s
        Int_t *fIndexPionL;       // Indices of identified pions from Lambda
	      Int_t *fIndexKaon;        // Indices of identified kaons
        Int_t *fIndexProton;      // Indices of identified protons
    };

    class AliHFELambdaInf{
      public:
        AliHFELambdaInf(AliKFParticle *track, AliKFVertex * const primaryVertex);
        ~AliHFELambdaInf();

        Double_t GetInvariantMass() const { return fInvariantMass; };
        Double_t GetChi2NDF() const { return fChi2NDF; };
        Double_t GetDistanceFromPrimaryVertex() const { return fDistVtx; };
      private:
        Double_t fInvariantMass;    // invariant mass before constraint
        Double_t fChi2NDF;          // chi2/ndf after constraints
        Double_t fDistVtx;          // Distance to primary Vertex
    };

    AliHFEV0pid(const AliHFEV0pid &ref);
    AliHFEV0pid&operator=(const AliHFEV0pid &ref);

    AliVEvent   *fInputEvent;        // Input Event
    AliKFVertex *fPrimaryVertex;     // Primary Vertex
    TObjArray   *fElectrons;         // List of Electron tracks coming from Conversions
    TObjArray   *fPionsK0;           // List of Pion tracks coming from K0
    TObjArray   *fPionsL;            // List of Pion tracks coming from L
    TObjArray   *fKaons;             // List of Kaon tracks from Phi decay
    TObjArray   *fProtons;           // List of Proton Tracks coming from Lambdas
    AliHFEV0pidTrackIndex *fIndices; // Container for Track indices
    AliHFEcollection *fQA;           // Collection of QA histograms

    ClassDef(AliHFEV0pid, 1)          // V0 PID Class

};

//____________________________________________________________
TList *AliHFEV0pid::GetListOfQAhistograms(){
  //
  // Get QA histograms
  //
  if(fQA)
    return fQA->GetList();
  return NULL;
}
#endif
