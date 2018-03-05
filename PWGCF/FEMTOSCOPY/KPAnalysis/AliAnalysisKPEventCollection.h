#ifndef AliAnalysisKPEventCollection_cxx
#define AliAnalysisKPEventCollection_cxx
#include <iostream>

#include <string>

#include "TH1.h"

#include "TH2.h"

// author: ramona.lea@cer.ch from 
// maria.nicassio@cern.ch (derived and adapted from D. Gangadharan PWGCF/FEMTOSCOPY/Chaoticity/AliChaoticityEventCollection
//                                             and J. Salzwedel PWGCF/FEMTOSCOPY/V0LamAnalysis/AliAnalysisV0LamEventCollection)

using namespace std;

class AliReconstructedFirst {
  
 public:
  AliReconstructedFirst();
  
  ~AliReconstructedFirst();
  
  /* AliReconstructedFirst(const AliReconstructedFirst&); */
  /* AliReconstructedFirst & operator=(const AliReconstructedFirst&); */
  
   enum MCFirstOrigin_t {kUnassigned, kFake, kFakeP, kPrimaryP, kPrimaryL, kOtherOriginP, kPrimaryAntiP, kPrimaryAntiL, kOtherOriginAntiP};
   double fMomentum[3]; // 3 reconstructed momentum
   double fMomentumTruth[3]; // 3 true momentum, used in momentum smearing analysis
   double fPt;
   double fEta;
   double fTheta;
   double fPhi;
   double fRap;
   short  fCharge;
   double fDCAxy;
   double fDCAz;
   double nSigmaFirstTPC[5];
   double nSigmaFirstTOF[5];
   bool   isTOFmismatch;
   bool   isMCptc;
   int    fMCcode;
   int    fPDGcode;
   int    fMCmumIdx;
   int    fMCmumPDG;
   int    fMCgrandmumIdx;
   int    fMCgrandmumPDG;
   int    index;
   MCFirstOrigin_t mcFirstOriginType;
   bool   doSkipOver;
   double iptoPV[2]; 

   Double_t fShiftedGlobalPosition[3];
   Double_t fEtaS;
   Double_t fPhiS;

   bool   isP;
   bool   isaP;

   ClassDef(AliReconstructedFirst, 1);   

};

class AliReconstructedSecond {


  public:
   AliReconstructedSecond();

   ~AliReconstructedSecond();

  /* AliReconstructedSecond(const AliReconstructedSecond&); */
  /* AliReconstructedSecond & operator=(const AliReconstructedSecond&); */

   enum MCSecondOrigin_t {kUnassigned, kFake, kFakeP, kPrimaryP, kPrimaryL, kOtherOriginP, kPrimaryAntiP, kPrimaryAntiL, kOtherOriginAntiP};
   double sMomentum[3]; // 3 reconstructed momentum
   double sMomentumTruth[3]; // 3 true momentum, used in momentum smearing analysis
   double sPt;
   double sEta;
   double sTheta;
   double sPhi;
   double sRap;
   short  sCharge;
   double sDCAxy;
   double sDCAz;
   double nSigmaSecondTPC[5];
   double nSigmaSecondTOF[5];
   bool   isTOFmismatch;
   bool   isMCptc;
   int    sMCcode;
   int    sPDGcode;
   int    sMCmumIdx;
   int    sMCmumPDG;
   int    sMCgrandmumIdx;
   int    sMCgrandmumPDG;
   int    index;
   MCSecondOrigin_t mcSecondOriginType;
   bool   doSkipOver;
   double iptoPV[2]; 

   Double_t sShiftedGlobalPosition[3];
   Double_t sEtaS;
   Double_t sPhiS;
   bool   isP;
   bool   isaP;
   
   ClassDef(AliReconstructedSecond, 1);   

};

class AliAnalysisKPEvent {

   public:

    int fNumberCandidateFirst;
    int fNumberCandidateSecond; 
    double fPrimaryVertex[3]; //Location of the primary vertex

    AliReconstructedFirst *fReconstructedFirst;
    AliReconstructedSecond *fReconstructedSecond;

    ClassDef(AliAnalysisKPEvent, 1);

};



class AliAnalysisKPEventCollection  {

   public:

    AliAnalysisKPEventCollection();

    AliAnalysisKPEventCollection(short eventBuffSize, int maxFirstMult, int maxSecondMult);
    /* AliAnalysisKPEventCollection(const AliAnalysisKPEventCollection&); */
    /* AliAnalysisKPEventCollection & operator=(const AliAnalysisKPEventCollection&); */
    
    ~AliAnalysisKPEventCollection();

    void FifoShift();

    AliAnalysisKPEvent *fEvt;

   private:

    short fifo; //Size of the Event Storage buffer

    void SetBuffSize(short eventBuffSize){fifo = eventBuffSize;}

   ClassDef(AliAnalysisKPEventCollection, 1);

};




#endif
