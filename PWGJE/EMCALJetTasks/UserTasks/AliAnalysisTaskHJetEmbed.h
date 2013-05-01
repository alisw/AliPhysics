#ifndef ALIANALYSISTASKHJETEMBED_H
#define ALIANALYSISTASKHJETEMBED_H

#include <vector>

class TH2F;
class TH1F;
class TF1;
class TH3F;
class THnSparse;
class TClonesArray;
class TObject;
class TString;
class AliAODEvent;
class AliESDEvent;
class AliMCEvent;
class AliRhoParameter;
class TRandom3;
class AliEmcalJet;
class AliVTrack;

template<class T> class TParameter;

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskHJetEmbed : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskHJetEmbed();
  AliAnalysisTaskHJetEmbed(const char *name);
  virtual ~AliAnalysisTaskHJetEmbed();

  void UserCreateOutputObjects();
  void UserExec(Option_t *option);
  void Terminate(Option_t *);
  void PrintConfig();

  void SetMaxVtxZ(Double_t z)                     { fMaxVtxZ=z;                     }
  void SetCollisionSystem(char *s)                { fCollisionSystem=s;             }
  void SetRunPeriod(char *p)                      { fPeriod=p;                      }
  void SetMCParticleArrName(char *s)              { fMCParticleArrName=s;           }
  void SetTrackArrName(char *s)                   { fTrackArrName=s;                }
  void SetTrkPtRange(Double_t min, Double_t max)  { fMinTrkPt=min; fMaxTrkPt=max;   }
  void SetTrkPhiRange(Double_t min, Double_t max) { fMinTrkPhi=min; fMaxTrkPhi=max; }
  void SetTrkEtaRange(Double_t min, Double_t max) { fMinTrkEta=min; fMaxTrkEta=max; }
  void SetRadius(Double_t rad)                    { fRadius=rad;                    }
  void SetPLJetArrName(char *s)                   { fPLJetArrName=s;                }
  void SetDLJetArrName(char *s)                   { fDLJetArrName=s;                }
  void SetJetArrName(char *s)                     { fJetArrName=s;                  }
  void SetRhoName(char *s)                        { fRhoName=s;                     }
  void SetRunQA(Bool_t run)                       { fRunQA=run;                     }
  void SetRunHJet(Bool_t run)                     { fRunHJet=run;                   }
  void SetRunMatch(Bool_t run)                    { fRunMatch=run;                  }



protected:
  void         RunQA();
  void         RunHJet();
  void         RunMatch();
  void         FillHJetCor(const TClonesArray *tracks, const Int_t leadingIndex, const TClonesArray *jetArray, THnSparse *hn, Bool_t isBkg = kFALSE);
  Int_t        FindGeoMatchedJet(const AliEmcalJet* jet, const TClonesArray *jetArray, Double_t &dR);
  Int_t        FindEnergyMatchedJet(const AliEmcalJet* jet, const TClonesArray *jetArray, Double_t &dR);
  Bool_t       AcceptTrack(const AliVParticle *track);
  Bool_t       IsGoodJet(const AliEmcalJet* jet);
  Double_t     GetZ(const Double_t trkPx, const Double_t trkPy, const Double_t trkPz, const Double_t jetPx, const Double_t jetPy, const Double_t jetPz);
  Double_t     GetDPhi(const Double_t phi1, const Double_t phi2);
  Double_t     CalculateDPhi(const Double_t phi1, const Double_t phi2);

  enum { kNTrig = 3  };

 private:
  Int_t             fVerbosity;            //! Control output
  Int_t             fAnaType;              //  0-local; 1-grid
  TString           fPeriod;               //  Run period
  TString           fCollisionSystem;      //  Collision system
  AliVEvent         *fEvent;               //! Input event 
  Int_t             fTriggerType;          //! Trigger type of the event
  Double_t          fCentrality;           //! V0M for current event
  Double_t          fMaxVtxZ;              //  Maximum z of vertices
  TString           fMCParticleArrName;    //
  TClonesArray      *fMCParticleArray;     //!
  TString           fTrackArrName;         //  Name of the input track array
  TClonesArray      *fTrackArray;          //! Array of input tracks
  Int_t             fTriggerTrkIndex;      //! Index of the trigger track in the event 
  Double_t          fMinTrkPt;             //
  Double_t          fMaxTrkPt;             //
  Double_t          fMinTrkEta;            //
  Double_t          fMaxTrkEta;            //
  Double_t          fMinTrkPhi;            //
  Double_t          fMaxTrkPhi;            //
  Double_t          fRadius;               //  Jet radius
  TString           fJetArrName;           //  Name of the found jet array
  TString           fPLJetArrName;         //  Name of the embedded PYTHIA jet array on particle level
  TString           fDLJetArrName;         //  Name of the embedded PYTHIA jet array on detector level
  TClonesArray      *fJetArray;            //! Array of the found jets
  TClonesArray      *fPLJetArray;          //! Array of the embedded PYTHIA jet array on particle level
  TClonesArray      *fDLJetArray;          //! Array of the embedded PYTHIA jet array on detector level
  TString           fRhoName;              //  Name of the rho parameter
  AliRhoParameter   *fRho;                 //! Rho parameter
  Double_t          fRhoValue;             //! Value of the rho parameter
  TParameter<int>   *fPtHardBinParam;      //!Pt hard bin param
  Int_t             fPtHardBin;            //!            

  Bool_t            fRunQA;                //  Flag to run QA
  Bool_t            fRunHJet;              //  Flag to run h+jet 
  Bool_t            fRunMatch;             //  Flag to run matching

  TList             *fOutputList;                //! Output list
  TH1F              *fhEventStat;                //!
  TH1F              *fhPtHardBins;               //!

  TH1F              *fhVtxZ[kNTrig];             //!
  TH1F              *fhCentrality[kNTrig];       //!
  TH2F              *fhRhoVsCent[kNTrig];        //!

  // QA histograms
  THnSparse         *fhPLJetPtVsCent[kNTrig];    //!
  THnSparse         *fhDLJetPtVsCent[kNTrig];    //!

  // HJet
  THnSparse         *fhPLTT[kNTrig];             //!
  THnSparse         *fhDLTT[kNTrig];             //!          
  THnSparse         *fhPLHJet[kNTrig];           //!
  THnSparse         *fhDLHJet[kNTrig];           //!
  THnSparse         *fhTTPtQA[kNTrig];           //!
  THnSparse         *fhTTPt[kNTrig];             //!
  THnSparse         *fhHJet[kNTrig];             //!


  // Matching
  THnSparse         *fhJetPtGeoMatch[kNTrig];    //!
  THnSparse         *fhJetPtEnMatch[kNTrig];     //!

  THnSparse         *fhJetPhiGeoMatch[kNTrig];   //!
  THnSparse         *fhJetPhiEnMatch[kNTrig];    //!

  AliAnalysisTaskHJetEmbed(const AliAnalysisTaskHJetEmbed&);            // not implemented
  AliAnalysisTaskHJetEmbed &operator=(const AliAnalysisTaskHJetEmbed&); // not implemented

  ClassDef(AliAnalysisTaskHJetEmbed, 1);
};

#endif
