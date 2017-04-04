#ifndef ALIANALYSISPIDTRACK_H
#define ALIANALYSISPIDTRACK_H

#include "TObject.h"
#include "TMath.h"
#include "TLorentzVector.h"
//#include "AliTOFGeometry.h"
//#include "AliTOFcalibHisto.h"
#include "AliTPCPIDResponse.h"
#include "AliTOFPIDResponse.h"
#include "AliAnalysisPIDEvent.h"
#include "AliESDtrack.h"
#include "AliPID.h"
#include "TF1.h"
#include "AliPIDResponse.h"

class AliStack;
class AliMCEvent;
class TH2F;
class TH1;

class AliAnalysisPIDTrack :
public TObject
{

 public:

  AliAnalysisPIDTrack(); // default constructor
  AliAnalysisPIDTrack(const AliAnalysisPIDTrack &source); // copy constructor
  AliAnalysisPIDTrack &operator=(const AliAnalysisPIDTrack &source); // operator=
  virtual ~AliAnalysisPIDTrack(); // default destructor

  Float_t GetP() const {return fP;}; // get p
  Float_t GetPt() const {return fPt;}; // get pt
  Float_t GetEta() const {return fEta;}; // get eta
  Float_t GetPhi() const {return fPhi;}; // get phi
  Float_t GetY(Float_t mass) const; // get Y
  Double_t GetSign() const {return fSign;}; // get sign
  ULong_t GetStatus() const {return fStatus;}; // get status
  Int_t GetLabel() const {return fLabel;}; // get label
  Float_t GetImpactParameter(Int_t i) const {return fImpactParameter[i];}; // get impact parameter
  //Float_t GetImpactParameterCov(Int_t i) const {return fImpactParameterCov[i];}; // get impact parameter covariance
  Float_t GetTPCmomentum() const {return fTPCmomentum;}; // get TPC inner wall momentum
  Float_t GetTPCdEdx() const {return fTPCdEdx;}; // get TPC dEdx
  UShort_t GetTPCdEdxN() const {return fTPCdEdxN;}; // get TPC dEdx clusters
  UShort_t GetTPCNcls() const {return fTPCNcls;}; // get number of clusters TPC
  UShort_t GetTPCNclsF() const {return fTPCNclsF;}; // get number of findable clusters TPC
  Float_t GetTPCnc() const {return fTPCNcr;}; // get number of crossed rows TPC
  Int_t GetTOFIndex() const {return fTOFIndex;}; // get TOF index
  Float_t GetTOFTime() const {return fTOFTime;}; // get TOF time
  Float_t GetTOFExpTime(Int_t i) const {return fTOFExpTime[i];}; // get TOF integrated times
  Float_t GetTOFExpTimeCorrection(Int_t i, Int_t chargeCorr = 0) const; // get TOF integrated times correction
  Float_t GetTOFExpTimeSigma(Int_t i) const; // get TOF integrated times sigma
  Float_t GetTOFExpectedSigma(Int_t i) const; // get TOF expected sigma
  Float_t GetTOFLength() const {return fTOFLength;}; // get TOF length
  Int_t GetTOFLabel(Int_t i) const {return fTOFLabel[i];}; // get TOF label
  Float_t GetTOFDeltaX() const {return fTOFDeltaX;}; // get TOF deltaX
  Float_t GetTOFDeltaZ() const {return fTOFDeltaZ;}; // get TOF deltaZ
  Bool_t IsMCPrimary() const {return fMCPrimary;}; // is MC primary
  Bool_t IsMCSecondaryWeakDecay() const {return fMCSecondaryWeak;}; // is MC weak decay
  Bool_t IsMCSecondaryMaterial() const {return fMCSecondaryMaterial;}; // is MC material
  Bool_t IsMCPrimary(Int_t ipart) const {return (fMCPrimary && TMath::Abs(fMCPdgCode) == AliPID::ParticleCode(ipart));}; // is MC primary
  Int_t GetMCPdgCode() const {return fMCPdgCode;}; // get MC PDG code
  Int_t GetMCMotherPdgCode() const {return fMCMotherPdgCode;}; // get MC mother PDG code
  Int_t GetMCMotherPrimary() const {return fMCMotherPrimary;}; // get MC mother primary
  Int_t GetMCMotherLabel() const {return fMCMotherLabel; }; //get MC mother label
  Bool_t IsMCTOFMatchPrimary() const {return fMCTOFMatchPrimary;};
  Int_t GetMCTOFMatchPdgCode() const {return fMCTOFMatchPdgCode;};
  Int_t GetMCTOFMatchLevel() const {return fMCTOFMatchLevel;};
  Float_t GetMCTOFTime() const {return fMCTOFTime;};
  Float_t GetMCTOFLength() const {return fMCTOFLength;};

  Bool_t IsInTOFPad() const {return ((TMath::Abs(fTOFDeltaX) < 1.25) && (TMath::Abs(fTOFDeltaZ) < 1.75));}; // is in TOF pad

  //Int_t GetTOFVolumeIndex(Int_t i); // get TOF volume index
  Float_t GetTOFBeta() const {return HasTOFPID() ? fTOFLength / (2.99792457999999984e-02 * (fTOFTime - fgTOFResponse->GetStartTime(fP))) : 0.;}; // get TOF beta
  Float_t GetTOFBetaSigma() const; // get TOF beta sigma
  Float_t GetTOFExpBeta(Int_t ipart) const {return HasTOFPID() ? fTOFLength / (2.99792457999999984e-02 * fTOFExpTime[ipart]) : 0.;}; // get TOF beta
  Float_t GetTOFBetaTh(Int_t ipart) const {return GetTOFBetaTh(AliPID::ParticleMass(ipart));}; // get TOF beta th
  Float_t GetTOFBetaTPCin(Int_t ipart) const {return GetTOFBetaTPCin(AliPID::ParticleMass(ipart));}; // get TOF beta th
  Float_t GetTOFBetaTh(Float_t mass) const {return TMath::Sqrt(1. / (1. + mass * mass / (fP * fP)));}; // get TOF beta th
  Float_t GetTOFBetaTPCin(Float_t mass) const {return TMath::Sqrt(1. / (1. + mass * mass / (fTPCmomentum * fTPCmomentum)));}; // get TOF beta th
  Float_t GetTOFMass2() const {return fP * fP * (1. / (GetTOFBeta() * GetTOFBeta()) - 1.);}; // get TOF mass^2
  Float_t GetTOFMass() const {return GetTOFMass2() > 0. ? TMath::Sqrt(TMath::Abs(GetTOFMass2())) : -TMath::Sqrt(TMath::Abs(GetTOFMass2()));}; // get TOF mass
  Float_t GetTPCdEdxTh(Float_t betagamma) const; // get TPC dEdx th
  Float_t GetTOFExpTimeTh(Int_t ipart) const {return fTOFLength / 2.99792457999999984e-02 / GetTOFBetaTh(ipart);}; // get TOF exp time th
  Float_t GetTOFExpTimeTPCin(Int_t ipart) const {return fTOFLength / 2.99792457999999984e-02 / GetTOFBetaTPCin(ipart);}; // get TOF exp time th
  Float_t GetTuningExpTimeTh(Int_t ipart) const; // get tuning exp time th
  Float_t GetTPCBetaGamma(Int_t ipart) const {return fTPCmomentum / AliPID::ParticleMass(ipart);}; // get TPC beta-gamma

  void Reset(); // reset
  void Update(AliESDtrack *track, AliStack *stack, AliMCEvent *mcevent, AliPIDResponse *PIDRes, Int_t TrackCutFlag); // update
  Bool_t HasTOFMatch() const {return (fStatus & AliESDtrack::kTOFout);}; // has TOF match
  Bool_t HasIntegratedTimes() const{ return (fStatus & AliESDtrack::kTIME);}; // has integrated times
  Bool_t HasTPCPID() const; // has TPC PID
  Bool_t HasTOFPID(TH1 *henabled = NULL) const; // has TOF PID
  Bool_t MakeTPCPID(Float_t *nSigma) const; // make TPC PID
  Bool_t MakeTOFPID(Float_t *nSigma) const; // make TOF PID
  Bool_t IsMismatchMC() const {return (fMCTOFMatchLevel != 0);}; // is mismatch MC
  Int_t GetMCPID() const; // get MC PID
  Int_t GetMCCharge() const; // get MC charge
  Bool_t IsUncorrelatedMismatchMC() const {return (fMCTOFMatchLevel < 0);}; // is uncorreltaed mismatch MC
  Bool_t IsMismatch(Float_t cutTPC = 5., Float_t cutTOF = 5.) const; // is mismatch
  Bool_t IsMismatch(const Float_t *nSigmaTPC, const Float_t *nSigmaTOF, Float_t cutTPC = 5., Float_t cutTOF = 5.) const; // is mismatch
  Bool_t IsBetaGammaCompatible(Float_t cutTPC = 5., Float_t cutTOF = 5.) const; // is beta-gamma compatible
  Bool_t IsHeavyAndCompatible(const Float_t *nSigmaTOF, Float_t cutTPC = 5., Float_t cutTOF = 5.) const; // is heavy and compatible
  void RemoveTimeZero(const AliAnalysisPIDEvent *analysisEvent); // remove time-zero
  Bool_t IsTPCPID(Int_t ipart, Float_t cutTPC = 5.) const; // is TPC PID
  Bool_t IsTPCDeuton() const; // is TPC deuton
  Bool_t IsTPCTriton() const; // is TPC deuton
  Bool_t IsTPCHeavy() const; // is TPC heavy
  //Bool_t HasPrimaryDCA(Float_t nSigmaXY = 7., Float_t nSigmaZ = 5.) const; // has primary DCA
  Bool_t IsTRDin() const {return fStatus & AliESDtrack::kTRDin;}; // is TRD in
  Bool_t IsTRDout() const {return fStatus & AliESDtrack::kTRDout;}; // is TRD out
  Bool_t IsTRDrefit() const {return fStatus & AliESDtrack::kTRDrefit;}; // is TRD refit
  Bool_t IsTOFout() const {return fStatus & AliESDtrack::kTOFout;}; // is TOF out
  Float_t GetNSigmaPionTPC() {return nSigmaPionTPC;};
  Float_t GetNSigmaKaonTPC() {return nSigmaKaonTPC;};
  Float_t GetNSigmaProtonTPC() {return nSigmaProtonTPC;};
  Float_t GetNSigmaPionTOF() {return nSigmaPionTOF;};
  Float_t GetNSigmaKaonTOF() {return nSigmaKaonTOF;};
  Float_t GetNSigmaProtonTOF() {return nSigmaProtonTOF;};
  Int_t GetTrackCutFlag() {return fTrackCutFlag; };
  Bool_t HasEMCal() { return fHasEMCal; };
  Float_t GetEMCalE() { return fEMCalE; };
  Float_t GetEMCalP() { return fEMCalP; };
  Float_t GetChiSq() { return fTPCchi2; };
  void SetEMCalPars(Float_t lEMCalE, Float_t lEMCalP) { fHasEMCal = kTRUE; fEMCalE = lEMCalE; fEMCalP = lEMCalP;};
  static Bool_t LoadTuningExpTimeTh(const Char_t *filename); // load tuning exp time th
  Bool_t IsAcceptedByTrackCuts(Int_t CutFlag) { return GetTrackCutFlag()&CutFlag;};
  //  Int_t GetTOFCalibIndex(Int_t imap) {return (Int_t)fgTOFcalibHisto.GetCalibMap(imap, fTOFIndex);}; // get TOF calib index

  static AliTOFPIDResponse *GetTOFResponse() {return fgTOFResponse;}; // getter
  static AliTPCPIDResponse *GetTPCResponse() {return fgTPCResponse;}; // getter
  static void SetTOFResponse(AliTOFPIDResponse *value) {fgTOFResponse = value;}; // setter
  static void SetTPCResponse(AliTPCPIDResponse *value) {fgTPCResponse = value;}; // setter

  static void UpdateTOFResponse(AliAnalysisPIDEvent *analysisEvent); // update TOF response

  void ApplyTOFExpectedTimeCorrection(Int_t chargeCorr = 0); // apply expected time correction

  Bool_t AcceptTrack(Bool_t selPrimaries = kTRUE); // accept track
  static void SetAcceptTrackMinNClustersTPC(UShort_t value) {fgMinNClustersTPC = value;}; // setter
  static void SetAcceptTrackClusterCut(Int_t value) {fgAcceptTrackClusterCut = value;}; // setter
  static void SetAcceptTrackMaxDCAToVertexXYPtDepFormula(TFormula *formula) {fgMaxDCAToVertexXYPtDepFormula = formula;}; // setter
  static void SetAcceptTrackEtaCut(Float_t value) {fgEtaCut = value;}; // setter
  static void SetAcceptTrackEtaReject(Float_t value) {fgEtaReject = value;}; // setter
  static void SetAcceptTrackStatusCut(ULong_t value) {fgAcceptTrackStatusCut = value;}; // setter
  static void SetRejectTrackStatusCut(ULong_t value) {fgRejectTrackStatusCut = value;}; // setter
  static void SetAcceptTrackMaxDCAToVertexZCut(Float_t value) {fgMaxDCAToVertexZCut = value;}; // setter
  static void SetAcceptTrackMaxChi2PerClusterTPC(Float_t value) {fgMaxChi2PerClusterTPC = value;}; // setter
  static void SetRejectITSFakes(Bool_t value) {fgRejectITSFakes = value;}; // setter
  static void SetMatchTrackDeltaX(Float_t value) {fgMatchTrackDeltaX = value;}; // setter
  static void SetMatchTrackDeltaZ(Float_t value) {fgMatchTrackDeltaZ = value;}; // setter

 private:

  /*** global track info ***/
  Float_t fP; // p
  Float_t fPt; // pt
  Float_t fEta; // eta
  Float_t fPhi; // phi
  Double_t fSign; // sign
  ULong_t fStatus; // status
  Int_t fLabel; // label
  Float_t fImpactParameter[2]; // impact parameters 
  //Float_t fImpactParameterCov[3]; // impact parameters covariance
  /*** TPC PID info ***/
  Float_t fTPCmomentum; // TPC inner wall momentum
  Float_t fTPCdEdx; // dEdx
  UShort_t fTPCdEdxN; // dEdx clusters
  UShort_t fTPCNcls; // number of clusters TPC
  UShort_t fTPCNclsF; // number of findable clusters TPC
  Float_t fTPCNcr; // number of crossed rows TPC
  /*** TOF PID info ***/
  Int_t fTOFIndex; // index
  Float_t fTOFTime; // time
  Float_t fTOFExpTime[AliPID::kSPECIES]; // integrated time array
  Float_t fTOFLength; // track length
  Float_t fTOFDeltaX; // TOF deltaX
  Float_t fTOFDeltaZ; // TOF deltaZ
  Int_t fTOFLabel[3]; // TOF label
  /*** MC info ***/
  Bool_t fMCPrimary; // MC primary flag
  Int_t fMCPdgCode; // MC PDG code
  Bool_t fMCMotherPrimary; // MC mother primary flag
  Int_t fMCMotherPdgCode; // MC mother PDG code
  Int_t fMCMotherLabel; //MC mother label
  Bool_t fMCTOFMatchPrimary; // MC TOF match primary flag
  Int_t fMCTOFMatchPdgCode; // MC TOF match PDG code
  Short_t fMCTOFMatchLevel; // MC TOF match level
  Float_t fMCTOFTime; // MC TOF time
  Float_t fMCTOFLength; // MC TOF length
  Bool_t fMCSecondaryWeak; 
  Bool_t fMCSecondaryMaterial; 
  /*** HMPID PID info ***/
  //Float_t fHMPIDmomentum;
  //Float_t fHMPIDsignal;
  Float_t nSigmaPionTPC;
  Float_t nSigmaKaonTPC;
  Float_t nSigmaProtonTPC;
  Float_t nSigmaPionTOF;
  Float_t nSigmaKaonTOF;
  Float_t nSigmaProtonTOF;
  /*** extras ***/
  Float_t fTPCchi2; // TPC chi2
  Bool_t fITSFakeFlag; // ITS fake flag
  Int_t fTrackCutFlag;
  Float_t fEMCalE;
  Float_t fEMCalP;
  Bool_t fHasEMCal;
  /*** cut paramters */
  static Float_t fgEtaCut; // eta cut
  static Float_t fgEtaReject; // eta reject
  static TFormula *fgMaxDCAToVertexXYPtDepFormula; // DCA-xy cut formula
  static UShort_t fgMinNClustersTPC; // cut
  static UShort_t fgMinNCrossedRowsTPC; // cut
  static Float_t fgMinRatioCrossedRowsOverFindableClustersTPC; // cut
  static Int_t fgAcceptTrackClusterCut; // cluster cut
  static ULong_t fgAcceptTrackStatusCut; // accept track status cut
  static ULong_t fgRejectTrackStatusCut; // reject track status cut
  static Float_t fgMaxDCAToVertexZCut; // DCA-z cut
  static Float_t fgMaxChi2PerClusterTPC; // max chi2 per cluster TPC cut
  static Bool_t fgRejectITSFakes; // reject ITS fakes cut
  static Float_t fgMatchTrackDeltaX; // match track deltaX
  static Float_t fgMatchTrackDeltaZ; // match track deltaZ

  /*** tools ***/
  static TLorentzVector fgLorentzVector;
  //static AliTOFGeometry fgTOFGeometry;
  //  static AliTOFcalibHisto fgTOFcalibHisto;
  //  static Bool_t fgTOFcalibHistoFlag;
  static AliTPCPIDResponse *fgTPCResponse;
  static AliTOFPIDResponse *fgTOFResponse;
  static TH2F *hTOFtuned_th[AliPID::kSPECIES];

  Float_t fTimeZeroSigma; //!

  ClassDef(AliAnalysisPIDTrack, 3);
};

#endif /* ALIANALYSISPIDTRACK_H */
