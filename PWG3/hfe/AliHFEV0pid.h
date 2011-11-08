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

#ifndef ROOT_TNamed
#include <TNamed.h>
#endif

class TObjArray;
class TList;
class TString;

class AliESDv0;
class AliESDtrack;
class AliKFParticle;
class AliKFVertex;
class AliVEvent;
class AliVTrack;
class AliMCEvent;

class AliHFEV0cuts;
class AliHFEcollection;

class AliHFEV0pid : public TNamed{
  public:
    AliHFEV0pid();
    AliHFEV0pid(const char *name);
    ~AliHFEV0pid();

    void  Process(AliVEvent * const inputEvent);
    Int_t ProcessV0(TObject *v0);
    void  Flush();

    void  InitQA();
    TList *GetListOfQAhistograms();

    TObjArray *GetListOfElectrons() const { return fElectrons; }
    TObjArray *GetListOfPionsK0() const { return fPionsK0; }
    TObjArray *GetListOfPionsL() const { return fPionsL; }
    TObjArray *GetListOfKaons() const { return fKaons; }
    TObjArray *GetListOfProtons() const { return fProtons; }

    Bool_t   IsAODanalysis() const { return TestBit(kAODanalysis); }
    Bool_t   IsESDanalysis() const { return !TestBit(kAODanalysis); }
    void     SetAODanalysis(Bool_t isAOD = kTRUE) { SetBit(kAODanalysis, isAOD); };
    void     SetESDanalysis(Bool_t isESD = kTRUE) { SetBit(kAODanalysis, !isESD); }; 

    void     SetMCEvent(AliMCEvent* const ev) { fMCEvent = ev; };
    
 protected:
    enum{
      kAODanalysis = BIT(14)
	};

    Int_t PreselectV0(AliESDv0 * const v0, Int_t idMC);

    void   ArmenterosPlotMC(AliESDv0 * const v0, Int_t idMC);
    Bool_t IsGammaConv(TObject *v0);
    Bool_t IsK0s(TObject *v0);
    Bool_t IsPhi(const TObject *v0) const;
    Bool_t IsLambda(TObject *v0);        
    TList *GetV0pidQA(); 

    Int_t IdentifyV0(TObject *v0, Int_t d[2]);

    void  BenchmarkV0finder();

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
    AliHFEV0pid(const AliHFEV0pid &ref);
    AliHFEV0pid&operator=(const AliHFEV0pid &ref);
    
    AliVEvent   *fInputEvent;        // Input Event
    Int_t        fNtracks;           // number of tracks in current event
    AliMCEvent  *fMCEvent;           // MC evnet
    Bool_t       fMCon;              // availability of MC information
    AliKFVertex *fPrimaryVertex;     // Primary Vertex
    TObjArray   *fElectrons;         // List of Electron tracks coming from Conversions
    TObjArray   *fPionsK0;           // List of Pion tracks coming from K0
    TObjArray   *fPionsL;            // List of Pion tracks coming from L
    TObjArray   *fKaons;             // List of Kaon tracks from Phi decay
    TObjArray   *fProtons;           // List of Proton Tracks coming from Lambdas

    TObjArray   *fGammas;            // for MC purposes - list of found gammas
    TObjArray   *fK0s;               // for MC purposes - list of found K0s
    TObjArray   *fLambdas;           // for MC purposes - list of found lambdas
    TObjArray   *fAntiLambdas;       // for MC purposes - list of found anti lambdas

    AliHFEV0pidTrackIndex *fIndices; // Container for Track indices
    AliHFEcollection *fQA;           // Collection of QA histograms
    AliHFEV0cuts     *fV0cuts;       // separate class for studying and applying the V0 cuts
    TList       *fOutput;            // collection list

    UInt_t       fDestBits;              // logical bits for destructor

    ClassDef(AliHFEV0pid, 1)          // V0 PID Class

};

#endif
