#ifndef AliAnalysisTaskTGReducedTree_cxx
#define AliAnalysisTaskTGReducedTree_cxx

// example of an analysis task creating a single electron Tree with V0Reader
// Authors: Taku Gunji

class TH1F;
class AliESDEvent;
class TTree;
class AliESDtrackCuts; 
class AliAnalysisFilter;
class AliDielectronPID;
class TList;
class AliV0ReaderV1;
//class AliDielectronTGReducedPair;
//class AliDielectronTGReducedTrack;
//class AliDielectronTGReducedInfo;
class AliESDv0KineCuts;

#include "AliAnalysisTaskSE.h"
#include <TObject.h>
#include <TClonesArray.h>


//_____________________________________________________________________

class AliDielectronTGReducedTrack : public TObject {

  friend class AliAnalysisTaskTGReducedTree; 

 public:
  AliDielectronTGReducedTrack();
  virtual ~AliDielectronTGReducedTrack();
  
  //getters
  Double_t  Px() const {return  fPx; }
  Double_t  Py() const {return  fPy; }
  Double_t  Pz() const {return  fPz; }
  Double_t  Pt() const {return  fPt; }
  Double_t  Xv() const {return  fXv; }
  Double_t  Yv() const {return  fYv; }
  Double_t  Zv() const {return  fZv; }
  Double_t  Phi() const {return  fPhi; }
  Double_t  Theta() const {return  fTheta; }
  Int_t  Charge() const {return  fCharge; }
  Double_t  ImpactParXY() const {return  fImpactParXY; }
  Double_t  ImpactParZ() const {return  fImpactParZ; }
  Double_t  ITSsignal() const {return  fITSsignal; }
  Double_t  ITSnSigmaEleRaw() const {return  fITSnSigmaEleRaw; }
  Double_t  ITSnSigmaEle() const {return  fITSnSigmaEle; }
  Double_t  NclsITS() const {return  fNclsITS; }
  Double_t  NclsSITS() const {return  fNclsSITS; }
  Double_t  NclsSFracITS() const {return  fNclsSFracITS; }
  Double_t  ITSchi2Cl() const {return  fITSchi2Cl; }
  Double_t  PIn() const {return  fPIn; }
  Double_t  POut() const {return  fPOut; }
  Double_t  TPCsignal() const {return  fTPCsignal; }
  Double_t  TOFsignal() const {return  fTOFsignal; }
  Double_t  TPCnSigmaEleRaw() const {return  fTPCnSigmaEleRaw; }
  Double_t  TPCnSigmaEle() const {return  fTPCnSigmaEle; }
  Double_t  TOFnSigmaEle() const {return  fTOFnSigmaEle; }
  Double_t  GoldenChi2() const {return  fGoldenChi2; }
  Double_t  TPCsignalN() const {return  fTPCsignalN; }
  Double_t  TPCsignalNfrac() const {return  fTPCsignalNfrac; }
  Double_t  TPCchi2Cl() const {return  fTPCchi2Cl; }
  Double_t  NclsTPC() const {return  fNclsTPC; }
  Double_t  NclsSTPC() const {return  fNclsSTPC; }
  Double_t  NFclsTPC() const {return  fNFclsTPC; }
  Double_t  NFclsTPCrFrac() const {return  fNFclsTPCrFrac; }
  Double_t  NFclsTPCfCross() const {return  fNFclsTPCfCross; }
  Double_t  TrackStatus() const {return  fTrackStatus; }
  Int_t  SelectInfo() const {return  fSelectInfo; }
  Int_t  ContVtx() const {return  fContVtx; }
  Int_t  PdgCode() const {return  fPdgCode; }
  Double_t  MCPx() const {return  fMCPx; }
  Double_t  MCPy() const {return  fMCPy; }
  Double_t  MCPz() const {return  fMCPz; }
  Double_t  MCXv() const {return  fMCXv; }
  Double_t  MCYv() const {return  fMCYv; }
  Double_t  MCZv() const {return  fMCZv; }
  Int_t  PdgMother() const {return  fPdgMother; }
  Int_t  MotherID() const {return  fMotherID; }
  Int_t  MotherLabel() const {return  fMotherLabel; }
  Int_t  GeneratorIndex() const {return  fGeneratorIndex; }
  Int_t  StatusAnotherLeg() const {return fStatusAnotherLeg;}
  Int_t  LabelAnotherLeg() const {return fLabelAnotherLeg;}
  Int_t  IDAnotherLeg() const {return fIDAnotherLeg;}
  Int_t  ESDV0Conv() const {return  fESDV0Conv; }
  Int_t  ESDV0ConvID2() const {return  fESDV0ConvID2; }
  Double_t  ESDV0Mass() const {return  fESDV0Mass; }
  Double_t  ESDV0R() const {return  fESDV0R; }
  Double_t  ESDV0Phi() const {return  fESDV0Phi; }
  Double_t  ESDV0PsiPair() const {return  fESDV0PsiPair; }
  Double_t  ESDV0PhivPair() const {return  fESDV0PhivPair; }
  Double_t  ESDV0CosP() const {return  fESDV0CosP; }                              
  Double_t  ESDV0ArmAlpha() const {return  fESDV0ArmAlpha; }                         
  Double_t  ESDV0ArmPt() const {return  fESDV0ArmPt; }                              
  Double_t  ESDV0OA() const {return  fESDV0OA; }                              
  Double_t  ESDV0cXv() const {return  fESDV0cXv; }                              
  Double_t  ESDV0cYv() const {return  fESDV0cYv; }                              
  Double_t  ESDV0cZv() const {return  fESDV0cZv; }                              
  Double_t  ESDV0Xv() const {return  fESDV0Xv; }                              
  Double_t  ESDV0Yv() const {return  fESDV0Yv; }                              
  Double_t  ESDV0Zv() const {return  fESDV0Zv; }                              
  Int_t  V0Conv() const {return  fV0Conv; }
  Int_t  V0ConvID2() const {return  fV0ConvID2; }
  Double_t  V0Mass() const {return  fV0Mass; }
  Double_t  V0Pt() const {return  fV0Pt; }
  Double_t  V0Phi() const {return  fV0Phi; }
  Double_t  V0Eta() const {return  fV0Eta; }
  Double_t  V0DCAr() const {return  fV0DCAr; }
  Double_t  V0DCAz() const {return  fV0DCAz; }
  Double_t  V0X() const {return  fV0X; }
  Double_t  V0Y() const {return  fV0Y; }
  Double_t  V0Z() const {return  fV0Z; }
  Double_t  V0R() const {return  fV0R; }
  Double_t  V0Alpha() const {return  fV0Alpha; }
  Double_t  V0Qt() const {return  fV0Qt; }
  Double_t  V0Psi() const {return  fV0Psi; }
  Int_t  ID() const {return  fID; }
  Int_t  Label() const {return  fLabel; }
  ULong_t  QualityFlags() const {return  fQualityFlags; }


  //setters
  void  Px(Double_t val){ fPx=val;}
  void  Py(Double_t val){ fPy=val;}
  void  Pz(Double_t val){ fPz=val;}
  void  Pt(Double_t val){ fPt=val;}
  void  Xv(Double_t val){ fXv=val;}
  void  Yv(Double_t val){ fYv=val;}
  void  Zv(Double_t val){ fZv=val;}
  void  Phi(Double_t val){ fPhi=val;}
  void  Theta(Double_t val){ fTheta=val;}
  void  Charge(Int_t val){ fCharge=val;}
  void  ImpactParXY(Double_t val){ fImpactParXY=val;}
  void  ImpactParZ(Double_t val){ fImpactParZ=val;}
  void  ITSsignal(Double_t val){ fITSsignal=val;}
  void  ITSnSigmaEleRaw(Double_t val){ fITSnSigmaEleRaw=val;}
  void  ITSnSigmaEle(Double_t val){ fITSnSigmaEle=val;}
  void  NclsITS(Double_t val){ fNclsITS=val;}
  void  NclsSITS(Double_t val){ fNclsSITS=val;}
  void  NclsSFracITS(Double_t val){ fNclsSFracITS=val;}
  void  ITSchi2Cl(Double_t val){ fITSchi2Cl=val;}
  void  PIn(Double_t val){ fPIn=val;}
  void  POut(Double_t val){ fPOut=val;}
  void  TPCsignal(Double_t val){ fTPCsignal=val;}
  void  TOFsignal(Double_t val){ fTOFsignal=val;}
  void  TPCnSigmaEleRaw(Double_t val){ fTPCnSigmaEleRaw=val;}
  void  TPCnSigmaEle(Double_t val){ fTPCnSigmaEle=val;}
  void  TOFnSigmaEle(Double_t val){ fTOFnSigmaEle=val;}
  void  GoldenChi2(Double_t val){ fGoldenChi2=val;}
  void  TPCsignalN(Double_t val){ fTPCsignalN=val;}
  void  TPCsignalNfrac(Double_t val){ fTPCsignalNfrac=val;}
  void  TPCchi2Cl(Double_t val){ fTPCchi2Cl=val;}
  void  NclsTPC(Double_t val){ fNclsTPC=val;}
  void  NclsSTPC(Double_t val){ fNclsSTPC=val;}
  void  NFclsTPC(Double_t val){ fNFclsTPC=val;}
  void  NFclsTPCrFrac(Double_t val){ fNFclsTPCrFrac=val;}
  void  NFclsTPCfCross(Double_t val){ fNFclsTPCfCross=val;}
  void  TrackStatus(Double_t val){ fTrackStatus=val;}
  void  SelectInfo(Int_t val){ fSelectInfo=val;}
  void  ContVtx(Int_t val){ fContVtx=val;}
  void  PdgCode(Int_t val){ fPdgCode=val;}
  void  MCPx(Double_t val){ fMCPx=val;}
  void  MCPy(Double_t val){ fMCPy=val;}
  void  MCPz(Double_t val){ fMCPz=val;}
  void  MCXv(Double_t val){ fMCXv=val;}
  void  MCYv(Double_t val){ fMCYv=val;}
  void  MCZv(Double_t val){ fMCZv=val;}
  void  PdgMother(Int_t val){ fPdgMother=val;}
  void  MotherID(Int_t val){ fMotherID=val;}
  void  MotherLabel(Int_t val){ fMotherLabel=val;}
  void  GeneratorIndex(Int_t val){ fGeneratorIndex=val;}
  void  StatusAnotherLeg(Int_t val){ fStatusAnotherLeg = val;}
  void  LabelAnotherLeg(Int_t val){ fLabelAnotherLeg = val;}
  void  IDAnotherLeg(Int_t val){ fIDAnotherLeg = val;}
  void  ESDV0Conv(Int_t val){ fESDV0Conv=val;}
  void  ESDV0ConvID2(Int_t val){ fESDV0ConvID2=val;}
  void  ESDV0Mass(Double_t val){ fESDV0Mass=val;}
  void  ESDV0R(Double_t val){ fESDV0R=val;}
  void  ESDV0Phi(Double_t val){ fESDV0Phi=val;}
  void  ESDV0PsiPair(Double_t val){ fESDV0PsiPair=val;}
  void  ESDV0PhivPair(Double_t val){ fESDV0PhivPair=val;}
  void  ESDV0CosP(Double_t val){ fESDV0CosP=val;}
  void  ESDV0ArmAlpha(Double_t val){ fESDV0ArmAlpha=val;}
  void  ESDV0ArmPt(Double_t val){ fESDV0ArmPt=val;}
  void  ESDV0OA(Double_t val){ fESDV0OA=val;}
  void  ESDV0cXv(Double_t val){ fESDV0cXv=val;}
  void  ESDV0cYv(Double_t val){ fESDV0cYv=val;}
  void  ESDV0cZv(Double_t val){ fESDV0cZv=val;}
  void  ESDV0Xv(Double_t val){ fESDV0Xv=val;}
  void  ESDV0Yv(Double_t val){ fESDV0Yv=val;}
  void  ESDV0Zv(Double_t val){ fESDV0Zv=val;}
  void  V0Conv(Int_t val){ fV0Conv=val;}
  void  V0ConvID2(Int_t val){ fV0ConvID2=val;}
  void  V0Mass(Double_t val){ fV0Mass=val;}
  void  V0Pt(Double_t val){ fV0Pt=val;}
  void  V0Phi(Double_t val){ fV0Phi=val;}
  void  V0Eta(Double_t val){ fV0Eta=val;}
  void  V0DCAr(Double_t val){ fV0DCAr=val;}
  void  V0DCAz(Double_t val){ fV0DCAz=val;}
  void  V0X(Double_t val){ fV0X=val;}
  void  V0Y(Double_t val){ fV0Y=val;}
  void  V0Z(Double_t val){ fV0Z=val;}
  void  V0R(Double_t val){ fV0R=val;}
  void  V0Alpha(Double_t val){ fV0Alpha=val;}
  void  V0Qt(Double_t val){ fV0Qt=val;}
  void  V0Psi(Double_t val){ fV0Psi=val;}
  void  ID(Int_t val){ fID=val;}
  void  Label(Int_t val){ fLabel=val;}
  void  QualityFlags(Int_t val){ fQualityFlags=val;}



 protected:
  Double_t fPx;              // px
  Double_t fPy;              // py
  Double_t fPz;              // pz
  Double_t fPt;                     // transverse momentum
  Double_t fXv;                     // vertex position in x
  Double_t fYv;                     // vertex position in y
  Double_t fZv;                     // vertex position in z
  Double_t fPhi;                    // phi angle
  Double_t fTheta;                  // theta angle
  Int_t fCharge;                 // charge
  Double_t fImpactParXY;            // Impact parameter in XY plane
  Double_t fImpactParZ;             // Impact parameter in Z
  Double_t fITSsignal;              // ITS dE/dx signal
  Double_t fITSnSigmaEleRaw;        // raw number of sigmas to the dE/dx electron line in the ITS
  Double_t fITSnSigmaEle;           // number of sigmas to the dE/dx electron line in the ITS
  Double_t fNclsITS;                // number of clusters assigned in the ITS
  Double_t fNclsSITS;               // number of shared clusters assigned in the ITS
  Double_t fNclsSFracITS;           // number of shared clusters assigned in the ITS
  Double_t fITSchi2Cl;              // chi2/cl in the ITS
  Double_t fPIn;                    // momentum at inner wall of TPC (if available), used for PID
  Double_t fPOut;                   // momentum at outer wall of TPC, used for TRD studies
  Double_t fTPCsignal;              // TPC dE/dx signal
  Double_t fTOFsignal;              // TOF signal
  Double_t fTPCnSigmaEleRaw;        // raw number of sigmas to the dE/dx electron line in the TPC
  Double_t fTPCnSigmaEle;           // number of sigmas to the dE/dx electron line in the TPC
  Double_t fTOFnSigmaEle;           // number of sigmas to the electron line in the TOF
  Double_t fGoldenChi2;             // Golden chi2
  Double_t fTPCsignalN;             // number of points used for dEdx
  Double_t fTPCsignalNfrac;         // fraction of points used for dEdx / cluster used for tracking
  Double_t fTPCchi2Cl;              // chi2/cl in TPC
  Double_t fNclsTPC;                // number of clusters assigned in the TPC
  Double_t fNclsSTPC;               // number of shared clusters assigned in the TPC
  Double_t fNFclsTPC;               // number of findable clusters in the TPC
  Double_t fNFclsTPCrFrac;          // number of found/findable clusters in the TPC with more robust definition
  Double_t fNFclsTPCfCross;         // fraction crossed rows/findable clusters in the TPC, as done in AliESDtrackCuts
  Double_t fTrackStatus;            // track status bits
  Int_t fSelectInfo; // selectInfo
  Int_t fContVtx; // track is used for vertex or not
  Int_t fPdgCode; // MC PDG code
  Double_t fMCPx;  // MC Px
  Double_t fMCPy;  // MC Py
  Double_t fMCPz; // MC Pz
  Double_t fMCXv;  // MC Xv
  Double_t fMCYv;  // MC Yv
  Double_t fMCZv; // MC Zv
  Int_t fPdgMother; // MC PDG mother
  Int_t fMotherID; // MC mother ID
  Int_t fMotherLabel; // MC mother label
  Int_t fGeneratorIndex; // generation index
  Int_t fStatusAnotherLeg; // to see whether another leg is reconstructed or not
  Int_t fLabelAnotherLeg; // label for another leg reconstructed
  Int_t fIDAnotherLeg; // ID for another leg reconstructed
  Int_t fESDV0Conv; // this track in V0 on the fly reconstructed 
  Int_t fESDV0ConvID2; // this track in V0 on the fly reconstructed 
  Double_t fESDV0Mass; // this track in V0 on the fly reconstructed 
  Double_t fESDV0R; // this track in V0 on the fly reconstructed 
  Double_t fESDV0Phi; // this track in V0 on the fly reconstructed 
  Double_t fESDV0PsiPair; // psi pair
  Double_t fESDV0PhivPair; // phiv pair
  Double_t fESDV0CosP; // CosPointingAngle
  Double_t fESDV0ArmAlpha; // ArmAlpha
  Double_t fESDV0ArmPt; // ArmPt
  Double_t fESDV0OA; // Opening Angle
  Double_t fESDV0cXv; // Xv from AliDielectronPair
  Double_t fESDV0cYv; // Yv
  Double_t fESDV0cZv; // Zv
  Double_t fESDV0Xv; // Xv from AliESDv0
  Double_t fESDV0Yv; // Yv
  Double_t fESDV0Zv; // Zv
  Int_t fV0Conv; // this track in V0 and passed PCM task
  Int_t fV0ConvID2; // counter part of fV0Conv
  Double_t fV0Mass; // V0 Mass
  Double_t fV0Pt; // V0 Pt
  Double_t fV0Phi; // V0 Phi of conversion
  Double_t fV0Eta; // V0 Eta
  Double_t fV0DCAr; // V0 DCAr
  Double_t fV0DCAz; // V0 DCAz
  Double_t fV0X; // V0 conversion X
  Double_t fV0Y; // V0 conversion Y
  Double_t fV0Z; // V0 conversion Z
  Double_t fV0R; // V0 conversion R
  Double_t fV0Alpha; // V0 Armenteros Alpha
  Double_t fV0Qt; // V0 Armenteros Qt
  Double_t fV0Psi; // V0 Psi Pair
  Int_t fID;  // ID  
  Int_t fLabel;  // Label
  Int_t fQualityFlags; // from ReducedTree
    

  //Double_t fgData[kNMaxValues];        // data


  AliDielectronTGReducedTrack(const AliDielectronTGReducedTrack &c);      
  AliDielectronTGReducedTrack& operator= (const AliDielectronTGReducedTrack &c);

  ClassDef(AliDielectronTGReducedTrack, 2);

};


//_____________________________________________________________________
class AliDielectronTGReducedPair : public TObject {
  
  friend class AliAnalysisTaskTGReducedTree; 

 public:
  AliDielectronTGReducedPair();
  virtual ~AliDielectronTGReducedPair();


  //getters
  Double_t  M() const {return  fM; }
  Double_t  Px() const {return  fPx; }
  Double_t  Py() const {return  fPy; }
  Double_t  Pz() const {return  fPz; }
  Double_t  Pt() const {return  fPt; }
  Double_t  Xv() const {return  fXv; }
  Double_t  Yv() const {return  fYv; }
  Double_t  Zv() const {return  fZv; }
  Double_t  Phiv() const {return  fPhiv; }
  Double_t  OpeningAngle() const {return  fOpeningAngle; }
  Double_t  LegDistXY() const {return  fLegDistXY; }
  Int_t     Label1() const{ return fLabel1;}
  Int_t     Label2() const{ return fLabel2;}
  Int_t     Index1() const{ return fIndex1;}
  Int_t     Index2() const{ return fIndex2;}
  Int_t     C1() const {return fC1;}
  Int_t     C2() const {return fC2;}

  //setters 
  void  M(Double_t val) { fM = val; }
  void  Px(Double_t val) {  fPx = val; }
  void  Py(Double_t val) {  fPy = val; }
  void  Pz(Double_t val) {  fPz = val; }
  void  Pt(Double_t val) {  fPt = val; }
  void  Xv(Double_t val) {  fXv = val; }
  void  Yv(Double_t val) {  fYv = val; }
  void  Zv(Double_t val) {  fZv = val; }
  void  Phiv(Double_t val) {  fPhiv = val; }
  void  OpeningAngle(Double_t val) {  fOpeningAngle = val; }
  void  LegDistXY(Double_t val) {  fLegDistXY = val; }
  void  Label1(Int_t val) { fLabel1 = val;}
  void  Label2(Int_t val) { fLabel2 = val;}
  void  Index1(Int_t val) { fIndex1 = val;}
  void  Index2(Int_t val) { fIndex2 = val;}
  void  C1(Int_t val) { fC1 = val;}
  void  C2(Int_t val) { fC2 = val;}

  
 private:
  Double_t fM;              // mass
  Double_t fPx;              // px
  Double_t fPy;              // py
  Double_t fPz;              // pz
  Double_t fPt;                     // transverse momentum
  Double_t fXv;                     // pair vertex position in x
  Double_t fYv;                     // pair vertex position in y
  Double_t fZv;                     // pair vertex position in z
  Double_t fPhiv; // phiv 
  Double_t fOpeningAngle ; // opening angle
  Double_t fLegDistXY; // leg dist XY
  Int_t fLabel1; // label1 in ESD
  Int_t fLabel2; // label2 in ESD
  Int_t fIndex1; // index1 in ReducedTrack
  Int_t fIndex2; // index2 in ReducedTrack
  Int_t fC1; // charge1
  Int_t fC2; // charge2

  AliDielectronTGReducedPair(const AliDielectronTGReducedPair &c) ;
  AliDielectronTGReducedPair& operator= (const AliDielectronTGReducedPair &c) ;


  ClassDef(AliDielectronTGReducedPair, 2)
  
};





//_____________________________________________________________________
class AliDielectronTGReducedInfo : public TObject {

  friend class AliAnalysisTaskTGReducedTree; 

 public:
  AliDielectronTGReducedInfo();
  virtual ~AliDielectronTGReducedInfo();
  
  
  AliDielectronTGReducedTrack* GetTrack(Int_t i) const {
    return (AliDielectronTGReducedTrack*)fTracks->At(i) ;}
  TClonesArray* GetTracks() const {return fTracks;}

  AliDielectronTGReducedPair* GetPair(Int_t i) const {
    return (AliDielectronTGReducedPair*)fPairs->At(i) ;}
  TClonesArray* GetPairs() const {return fPairs;}

  void ClearEvent();

  //getters
  Int_t     RunNo()  const {return fRun;}
  Int_t     EvtNo()  const {return fEvt;}
  Double_t  Vx()     const {return fVtx[0];}
  Double_t  Vy()     const {return fVtx[1];}
  Double_t  Vz()     const {return fVtx[2];}
  Int_t     VertexNContributors() const {return fVtxCont;}
  Int_t     NElectrons() const {return fNe;}
  Int_t     NV0Candidates() const {return fV0Cand;}

  //setters 
 protected:  
  ULong64_t fEventTag;              // Event tags to be used either during analysis or to filter events   
  Int_t fRun;  // Run number 
  Int_t fEvt;  // Event number 
  Double_t fVtx[3]; // vertex 
  Int_t  fVtxCont; // number of vertex contributors 
  Int_t fNe; // number of electrons 
  Int_t fV0Cand; // number of vertex candidates

  TClonesArray* fTracks;                // array containing particles
  TClonesArray* fPairs;                // array containing Pairs
  
  AliDielectronTGReducedInfo(const AliDielectronTGReducedInfo &c);      
  AliDielectronTGReducedInfo& operator= (const AliDielectronTGReducedInfo &c);

  ClassDef(AliDielectronTGReducedInfo, 3);

};




//_____________________________________________________________________
class AliAnalysisTaskTGReducedTree : public AliAnalysisTaskSE {

 public:
 AliAnalysisTaskTGReducedTree(): 
  AliAnalysisTaskSE(), fESD(0), fOutputList(0), fTree(0), 
    fEventStat(0), fESDtrackCuts(0), fTrackFilter(0), fMCEvent(0), 
    fPIDResponse(0),
    fPIDCuts(0), 
    fReducedInfo(0), fV0OpenCuts(0),
    fTriggerMask(0), fSelectPhysics(0),
    fFiredTrigger(""), fFiredExclude(kFALSE), fConvCut(117), 
    fCutArray(0), fMesonCutArray(0), fGammaCandidates(0),
    fV0Reader(0), fInputEvent(0), fReaderGammas(0), hasMC(kFALSE), 
    fEvalEfficiencyFlag(kFALSE), 
    fEvalEfficiencyMCset(-1), 
    fEvalEfficiencyIndex(-1), 
    fEvalEfficiencyParticle(-1),
    hGen(NULL), hMul(NULL)
    {
      hPtRap[0] = NULL; hPtRap[1] = NULL; hPtRap[2]=NULL;
      hPtRapConv[0] = NULL; hPtRapConv[1] = NULL; hPtRapConv[2]=NULL;
      hRPtConv[0] = NULL; hRPtConv[1] = NULL; hRPtConv[2]=NULL;
    }

  AliAnalysisTaskTGReducedTree(const char *name);
  virtual ~AliAnalysisTaskTGReducedTree(); 
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);

  void SetTriggerMask(ULong64_t mask) {fTriggerMask=mask;}
  void UsePhysicsSelection(Bool_t phy) {fSelectPhysics=phy;}
  void SetGammaConvCut(TString photoncut="00000069227300008250400000");
  void SetESDTrackCuts(AliESDtrackCuts *fCuts){ fESDtrackCuts = fCuts;}
  void SetPIDCuts(AliDielectronPID *fCuts){ fPIDCuts = fCuts;}
  void SetV0OpenCuts(AliESDv0KineCuts* const cuts) {fV0OpenCuts = cuts;}
  void SetEvalEfficiency(Bool_t flag, 
			 TString ProdName,
			 Int_t GeneratorIndex, 
			 Int_t ParticleIndex){
    fEvalEfficiencyFlag = flag;
    fEvalEfficiencyMCset = ProdName;
    fEvalEfficiencyIndex = GeneratorIndex;
    fEvalEfficiencyParticle = ParticleIndex;
  } 
  void SetFiredTriggerName(const char* select, Bool_t exclude=kFALSE){fFiredTrigger=select; fFiredExclude=exclude;}
  

 private:
  
  AliESDEvent *fESD;    // ESD object
  TList       *fOutputList; //! Output list
  TTree       *fTree; //! Tree object
  TH1D        *fEventStat; //! Event stat
  AliESDtrackCuts *fESDtrackCuts; // ESD track cuts (now defined in the task)
  AliAnalysisFilter *fTrackFilter; // track filter 
  AliMCEvent  *fMCEvent;  // MC event object
  AliPIDResponse *fPIDResponse;     // PID response object  
  AliDielectronPID *fPIDCuts;    // PID cuts

  AliDielectronTGReducedInfo *fReducedInfo; // reduced info 
  AliESDv0KineCuts *fV0OpenCuts ;// V0 open cuts


  UInt_t fTriggerMask;  // Event trigger mask
  Bool_t fSelectPhysics; // Whether to use physics selection
  TString fFiredTrigger; // online trigger class name
  Bool_t  fFiredExclude; // online trigger class name excluded
  Int_t  fConvCut; // Conversion cut 

  TList* fCutArray; // Photon cut array
  TList* fMesonCutArray; // Meson cut array
  TList* fGammaCandidates; // Photon candidates from on-the-fly V0

  AliV0ReaderV1*  fV0Reader;    // V0Reader from PCM
  AliVEvent*      fInputEvent;  // input event
  TClonesArray*   fReaderGammas; // Gamma info

  Bool_t hasMC;   //hasMC
  Bool_t fEvalEfficiencyFlag; // flag to look at pure MC
  TString fEvalEfficiencyMCset ; // MC name for the efficiency evaluation
  Int_t fEvalEfficiencyIndex; // GeneratorIndex to look at pure MC for this Generator index
  Int_t fEvalEfficiencyParticle; // Mother Particle to look at pure MC for this Generator index

  TH2F  *hPtRap[3]; //! All electrons/positrons generated/reco/accpeted
  TH2F  *hPtRapConv[3]; //! All conversion electrons/positrons generated/reco/accpeted
  TH2F  *hRPtConv[3]; //! All electrons/positrons from conversions generated/reco/accpeted (AliESDTrackCuts)
  TH1F   *hGen ; //! generator index
  TH2F   *hMul ; //! generator index vs. event activity
  

  AliAnalysisTaskTGReducedTree(const AliAnalysisTaskTGReducedTree&); // not implemented
  AliAnalysisTaskTGReducedTree& operator=(const AliAnalysisTaskTGReducedTree&); // not implemented
  
  ClassDef(AliAnalysisTaskTGReducedTree, 2); // example of analysis
};

#endif



