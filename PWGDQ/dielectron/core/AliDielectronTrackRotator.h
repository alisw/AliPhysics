#ifndef ALIDIELECTRONTRACKROTATOR_H
#define ALIDIELECTRONTRACKROTATOR_H

/* Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//#############################################################
//#         Class AliDielectronTrackRotator                   #
//#                                                           #
//#  Authors:                                                 #
//#   Anton     Andronic, GSI / A.Andronic@gsi.de             #
//#   Ionut C.  Arsene,   GSI / I.C.Arsene@gsi.de             #
//#   Julian    Book,     Uni Ffm / Julian.Book@cern.ch       #
//#   Frederick Kramer,   Uni Ffm, / Frederick.Kramer@cern.ch #
//#   Magnus    Mager,    CERN / Magnus.Mager@cern.ch         #
//#   WooJin J. Park,     GSI / W.J.Park@gsi.de               #
//#   Jens      Wiechula, Uni HD / Jens.Wiechula@cern.ch      #
//#                                                           #
//#############################################################

#include <TNamed.h>


#include <AliKFParticle.h>
#include <AliLog.h>

#include <vector>

class TObjArray;
class AliVTrack;
class AliVEvent;


class AliDielectronTrackRotator : public TNamed {
public:
  enum ERotationType {kRotatePositive, kRotateNegative, kRotateBothRandom};
  enum TypeOfModernRotation{fRotateULS, fRotateWithTGenPhaseSpace, fRotateLS};       // Use TGenPhaseSpace
  AliDielectronTrackRotator();
  AliDielectronTrackRotator(const char*name, const char* title);

  virtual ~AliDielectronTrackRotator();

  void SetTrackArrays(const TObjArray * const arrP, const TObjArray * const arrN) {fkArrTracksP=arrP;fkArrTracksN=arrN;}
  void SetTrackArraysRotation(TObjArray * const arrP, TObjArray * const arrN) {fkArrTracksPRotation=arrP;fkArrTracksNRotation=arrN;}
  void Reset();
  Bool_t NextCombination();

  //Setters
  void SetIterations(UInt_t niter)         { fIterations=niter;  }
  void SetRotationType(ERotationType type) { fRotationType=type; }
  void SetStartAnglePhi(Double_t phi)      { fStartAnglePhi=phi; }
  void SetConeAnglePhi(Double_t phi)       { fConeAnglePhi=phi;  }
  void SetKeepLocalY(Bool_t keep)          { fKeepLocalY=keep;  }
  void SetRotateAroundMother(Bool_t mother){ fRotateAroundMother=mother; }
  void SetUseTypeOfModernRotation(TypeOfModernRotation b)      { fTypeOfModernRotation = b; }
  void SetPtCut(Double_t ptCutMin, Double_t ptCutMax)     { fMinimalPtCut = ptCutMin; fMaximalPtCut = ptCutMax;}
  void SetEtaCut(Double_t etaCutMin, Double_t etaCutMax)     { fMinimalEtaCut = etaCutMin; fMaximalEtaCut = etaCutMax;}
  void SetRotatedTrackWeightMap(TString filename, TString histoname);
  void SetRotatedPairWeightMap(TString filename, TString histoname);
  void SetRotatedPairWeightMap2(TString filename, TString histoname);
  void SetUseAcceptanceMap(Bool_t use)     {fUseAccMap = use;}
  void SetRotWeightMinPtBin (Int_t pTbin)  {fRotWeight_minPtBin = pTbin;}
  void SetRotWeightMaxPtBin (Int_t pTbin)  {fRotWeight_maxPtBin = pTbin;}


  //Getters
  Int_t GetIterations() const           { return fIterations;    }
  ERotationType GetRotationType() const { return fRotationType;  }
  Double_t GetStartAnglePhi() const     { return fStartAnglePhi; }
  Double_t GetConeAnglePhi() const      { return fConeAnglePhi;  }
  Bool_t GetKeepLocalY() const          { return fKeepLocalY;  }
  Bool_t GetRotateAroundMother() const  { return fRotateAroundMother; }

  void SetEvent(AliVEvent * const ev)   { fEvent = ev;           }
  void SetPdgLegs(Int_t pdfLeg1, Int_t pdfLeg2) { fPdgLeg1=pdfLeg1; fPdgLeg2=pdfLeg2; }

  const AliKFParticle& GetKFTrack1() const {return fTrack1;}
  const AliKFParticle& GetKFTrack2() const {return fTrack2;}

  AliVTrack* GetVTrack1() const {return fVTrackP;}
  AliVTrack* GetVTrack2() const {return fVTrackN;}
  Int_t GetChargeTrack1() const {return fChargeTrack1;}
  Int_t GetChargeTrack2() const {return fChargeTrack2;}
  Bool_t SameTracks() const {return fSameTracks;}
  void ClearRotatedTrackPool() {  fRotatedTracksP.clear(); fRotatedTracksN.clear(); }
  void ClearRotatedPairPool()  {  fArrTrackPairs.clear();  fArrTrackPairsPM.clear();  fArrTrackPairsPP.clear();  fArrTrackPairsMM.clear();  }
  Int_t GetRotatedTrackPSize() { return fRotatedTracksP.size(); }
  Int_t GetRotatedTrackNSize() { return fRotatedTracksN.size(); }
  AliKFParticle* GetRotatedTrackP(Int_t track) { return &fRotatedTracksP[track]; }
  AliKFParticle* GetRotatedTrackN(Int_t track) { return &fRotatedTracksN[track]; }
  Double_t GetRotatedTrackWeightP(Int_t track) { return fRotatedTracksWeightP[track]; }
  Double_t GetRotatedTrackWeightN(Int_t track) { return fRotatedTracksWeightN[track]; }
  Double_t GetWeightFromRotation(AliKFParticle* part);
  Double_t GetWeightFromRotation2(Double_t rotAng);
  Double_t GetWeightFromOpeningAngle(AliKFParticle* KFpos, AliKFParticle* KFneg);
  Double_t GetWeightForPair(){return fWeight;};
  //void RotateKFParticle(AliKFParticle * kfParticle,Double_t angle, TVector3 *axis, const AliVEvent * const ev=0x0);
  void RotateKFParticle(AliKFParticle * kfParticle,Double_t angle, AliKFParticle * kfMother, const AliVEvent * const ev=0x0);
  AliESDtrack RotateVTrack(AliVTrack* Vtrack, Double_t angle, TLorentzVector* LvecMother, const AliVEvent * const ev);
  Int_t IJ( Int_t i, Int_t j );
  Int_t fRotWeight_minPtBin;
  Int_t fRotWeight_maxPtBin;
  

  class KFClass{
  public:
    KFClass (AliKFParticle kf1, AliKFParticle kf2, AliESDtrack vtrack1, AliESDtrack vtrack2, Short_t charged_tracks, Double_t weight) :
            kf1(kf1), kf2(kf2), vtrack1(vtrack1), vtrack2(vtrack2), charged_tracks(charged_tracks), weight(weight), rotAng(0.){}
    KFClass() : kf1(), kf2(), vtrack1(), vtrack2(), charged_tracks(0), weight(0.), rotAng(0.){}
    AliKFParticle kf1;
    AliKFParticle kf2;
    AliESDtrack vtrack1;
    AliESDtrack vtrack2;
    Short_t charged_tracks; //0: ULS, 1: LS_PP, 2: LS_MM
    Double_t weight;
    Double_t rotAng;
  };
  // std::vector<AliAODTrack> fRotatedTracks;
 
  KFClass* GetRotatedPair(Int_t pair) { return &fArrTrackPairs[pair]; } 
  KFClass* GetRotatedPairPM(Int_t pair) { return &fArrTrackPairsPM[pair]; }
  KFClass* GetRotatedPairPP(Int_t pair) { return &fArrTrackPairsPP[pair]; } 
  KFClass* GetRotatedPairMM(Int_t pair) { return &fArrTrackPairsMM[pair]; } 

private:
  UInt_t   fIterations;             // number of iterations

  ERotationType fRotationType;      // which track to rotate

  Double_t fStartAnglePhi;          // starting angle for rotation
  Double_t fConeAnglePhi;           // opening angle in phi for multiple rotation
  Bool_t fKeepLocalY;               // rotate on such an angle that the position wrt TPC chamber borders stays the same
  Bool_t fRotateAroundMother;       // Rotate around mother
  TypeOfModernRotation fTypeOfModernRotation;

  const TObjArray *fkArrTracksP;    //! array of positive tracks
  const TObjArray *fkArrTracksN;    //! array of negative tracks
  TObjArray *fkArrTracksPRotation;    //! array of positive rotated tracks
  TObjArray *fkArrTracksNRotation;    //! array of negative rotated tracks
  
  UInt_t   fCurrentIteration;       //! current iteration step
  Int_t    fCurrentPairPP;           //! current positive track in array
  Int_t    fCurrentPairMM;           //! current negative track in array
  Int_t    fCurrentTackP;           //! current positive track in array
  Int_t    fCurrentTackN;           //! current negative track in array

  Bool_t   fLastPairSent;

  AliVEvent *fEvent;                //! current event
  
  AliKFParticle fTrack1;            //! Positive track
  AliKFParticle fTrack2;            //! Negative track
  
  AliVTrack *fVTrackP;              //! Positive track
  AliVTrack *fVTrackN;              //! Negative track

  Int_t    fChargeTrack1;
  Int_t    fChargeTrack2;
  
  Int_t fPdgLeg1;                   //! pdg code leg1
  Int_t fPdgLeg2;                   //! pdg code leg2
  
  Bool_t fSameTracks;               //! tracks in both arrays at current position are the same
  Bool_t fUseAccMap;

  std::vector<AliKFParticle> fRotatedTracksP;
  std::vector<AliKFParticle> fRotatedTracksN;
  std::vector<Double_t> fRotatedTracksWeightP;
  std::vector<Double_t> fRotatedTracksWeightN;

  std::vector<KFClass> fArrTrackPairs;
  std::vector<KFClass> fArrTrackPairsPM;
  std::vector<KFClass> fArrTrackPairsPP;
  std::vector<KFClass> fArrTrackPairsMM;
  Double_t fMinimalPtCut;
  Double_t fMaximalPtCut;
  Double_t fMinimalEtaCut;
  Double_t fMaximalEtaCut;
  Double_t fWeight;

  TH3F fRotateTrackCorrectionMap;
  TH2F fRotatePairCorrectionMap;
  TH1F fRotatePairCorrectionMap2;

  Bool_t RotateTracks();
  void CalculatePairsFromRotationAroundMother();
  void CalculateLikeSignPairs();
  Double_t PhivPair(Double_t MagField, Int_t charge1, Int_t charge2, TVector3 dau1, TVector3 dau2); //const


  AliDielectronTrackRotator(const AliDielectronTrackRotator &c);
  AliDielectronTrackRotator &operator=(const AliDielectronTrackRotator &c);

  ClassDef(AliDielectronTrackRotator,2)         // Dielectron TrackRotator
};



#endif
