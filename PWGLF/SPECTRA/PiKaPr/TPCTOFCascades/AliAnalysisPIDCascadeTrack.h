#ifndef ALIANALYSISPIDCASCADETRACK_H
#define ALIANALYSISPIDCASCADETRACK_H

#include "TObject.h"
#include "TMath.h"
#include "TLorentzVector.h"
//#include "AliTOFGeometry.h"
//#include "AliTOFcalibHisto.h"
#include "AliTPCPIDResponse.h"
#include "AliTOFPIDResponse.h"
#include "AliAnalysisPIDCascadeEvent.h"
#include "AliESDtrack.h"
#include "AliPID.h"
#include "TF1.h"
#include "AliPIDResponse.h"

//class AliStack;
class AliMCEvent;
class TH2F;
class TH1;

class AliAnalysisPIDCascadeTrack :
public TObject
{

 public:

  AliAnalysisPIDCascadeTrack();//default constructor
  AliAnalysisPIDCascadeTrack(const AliAnalysisPIDCascadeTrack &source); //copy constructor
  AliAnalysisPIDCascadeTrack &operator=(const AliAnalysisPIDCascadeTrack &source);
  virtual ~AliAnalysisPIDCascadeTrack(); //default destructor.

  
  //General Track Info
  Float_t GetP() const {return fP;}; // get p
  Float_t GetPt() const {return fPt;}; // get pt
  Float_t GetEta() const {return fEta;}; // get eta
  Float_t GetPhi() const {return fPhi;}; // get phi
  Float_t GetY(Float_t mass) const; // get Y
  Double_t GetSign() const {return fSign;}; // get sign
  ULong64_t GetStatus() const {return fStatus;}; // get status
  Int_t GetLabel() const {return fLabel;}; // get label
  Float_t GetImpactParameter(Int_t i) const {return fImpactParameter[i];}; // get impact parameter
  //TPC Info
  Float_t GetTPCdEdx() const {return fTPCdEdx;}; // get TPC dEdx
  UShort_t GetTPCdEdxN() const {return fTPCdEdxN;}; // get TPC dEdx clusters

  //TOF Info
  Int_t GetTOFIndex() const {return fTOFIndex;}; // get TOF index

 Float_t GetTOFLength() const {return fTOFLength;}; // get TOF length
  Int_t GetTOFLabel(Int_t i) const {return fTOFLabel[i];}; // get TOF label
  Float_t GetTOFDeltaX() const {return fTOFDeltaX;}; // get TOF deltaX
  Float_t GetTOFDeltaZ() const {return fTOFDeltaZ;}; // get TOF deltaZ
  //MC Info
  Bool_t IsMCPrimary() const {return fMCPrimary;}; // is MC primary
  Bool_t IsMCSecondaryWeakDecay() const {return fMCSecondaryWeak;}; // is MC weak decay
  Bool_t IsMCSecondaryMaterial() const {return fMCSecondaryMaterial;}; // is MC material
  Bool_t IsMCPrimary(Int_t ipart) const {return (fMCPrimary && TMath::Abs(fMCPdgCode) == AliPID::ParticleCode(ipart));}; // is MC primary
  Int_t GetMCPdgCode() const {return fMCPdgCode;}; // get MC PDG code
  Int_t GetMCMotherPdgCode() const {return fMCMotherPdgCode;}; // get MC mother PDG code
  Int_t GetMCMotherPrimary() const {return fMCMotherPrimary;}; // get MC mother primary
  Int_t GetMCMotherLabel() const {return fMCMotherLabel; }; //get MC mother label
  Int_t GetMCPrimaryPdgCode() const {return fMCPrimaryPdgCode;}; // get MC PDG code of primary (grand)mother
  Int_t GetMCPrimaryLabel() const {return fMCPrimaryLabel; }; //get MC label of primary (grand)mother
  Bool_t IsMCTOFMatchPrimary() const {return fMCTOFMatchPrimary;};
  Int_t GetMCTOFMatchPdgCode() const {return fMCTOFMatchPdgCode;};
  Int_t GetMCTOFMatchLevel() const {return fMCTOFMatchLevel;};
  Float_t GetMCTOFTime() const {return fMCTOFTime;};
  Float_t GetMCTOFLength() const {return fMCTOFLength;};


  void Reset(); // reset
  void Update(AliESDtrack *track, AliMCEvent *mcevent, AliPIDResponse *PIDRes, Int_t TrackCutFlag); // update
  Bool_t HasTOFMatch() const {return (fStatus & AliESDtrack::kTOFout);}; // has TOF match
  Bool_t HasTPCPID() const; // has TPC PID
  Bool_t HasTOFPID(TH1 *henabled = NULL) const; // has TOF PID
  Int_t GetMCPID() const; // get MC PID
  Int_t GetMCCharge() const; // get MC charge
  Bool_t IsTOFout() const {return fStatus & AliESDtrack::kTOFout;}; // is TOF out
  Float_t GetNSigmaPionTPC() {return nSigmaPionTPC;};
  Float_t GetNSigmaKaonTPC() {return nSigmaKaonTPC;};
  Float_t GetNSigmaProtonTPC() {return nSigmaProtonTPC;};
  Float_t GetNSigmaPionTOF() {return nSigmaPionTOF;};
  Float_t GetNSigmaKaonTOF() {return nSigmaKaonTOF;};
  Float_t GetNSigmaProtonTOF() {return nSigmaProtonTOF;};
  Int_t GetTrackCutFlag() {return fTrackCutFlag; };

  // Bool_t IsAcceptedByTrackCuts(Int_t CutFlag) { return GetTrackCutFlag()&CutFlag;};


  static AliTOFPIDResponse *GetTOFResponse() {return fgTOFResponse;}; // getter
  static AliTPCPIDResponse *GetTPCResponse() {return fgTPCResponse;}; // getter
  static void SetTOFResponse(AliTOFPIDResponse *value) {fgTOFResponse = value;}; // setter
  static void SetTPCResponse(AliTPCPIDResponse *value) {fgTPCResponse = value;}; // setter
  static void UpdateTOFResponse(AliAnalysisPIDCascadeEvent *analysisEvent); // update TOF response

  //Bool_t AcceptTrack(Bool_t selPrimaries = kTRUE); // accept track

  static void SetMatchTrackDeltaX(Float_t value) {fgMatchTrackDeltaX = value;}; // setter
  static void SetMatchTrackDeltaZ(Float_t value) {fgMatchTrackDeltaZ = value;}; // setter
  static AliTPCPIDResponse *fgTPCResponse;
  static AliTOFPIDResponse *fgTOFResponse;


 private:

  /*** global track info ***/
  Float_t fP; // p
  Float_t fPt; // pt
  Float_t fEta; // eta
  Float_t fPhi; // phi
  Double_t fSign; // sign
  ULong64_t fStatus; // status
  Int_t fLabel; // label
  Float_t fImpactParameter[2]; // impact parameters
  /*** TPC PID info ***/
  Float_t fTPCdEdx; // dEdx
  UShort_t fTPCdEdxN; // dEdx clusters
  /*** TOF PID info ***/
  Int_t fTOFIndex; // index
  Float_t fTOFLength; // track length
  Float_t fTOFDeltaX; // TOF deltaX
  Float_t fTOFDeltaZ; // TOF deltaZ
  Int_t fTOFLabel[3]; // TOF label
  /*** MC info ***/
  Bool_t fMCPrimary; // MC primary flag
  Int_t fMCPdgCode; // MC PDG code
  Bool_t fMCMotherPrimary; // MC mother primary flag
  Int_t fMCMotherPdgCode; // MC mother PDG code
  Int_t fMCMotherLabel; // MC mother label
  Int_t fMCPrimaryPdgCode; // MC PDG code of primary (grand)mother
  Int_t fMCPrimaryLabel; // MC label of primary (grand)mother
  Bool_t fMCTOFMatchPrimary; // MC TOF match primary flag
  Int_t fMCTOFMatchPdgCode; // MC TOF match PDG code
  Short_t fMCTOFMatchLevel; // MC TOF match level
  Float_t fMCTOFTime; // MC TOF time
  Float_t fMCTOFLength; // MC TOF length
  Bool_t fMCSecondaryWeak;
  Bool_t fMCSecondaryMaterial;
  /*** PID info ***/
  Float_t nSigmaPionTPC;
  Float_t nSigmaKaonTPC;
  Float_t nSigmaProtonTPC;
  Float_t nSigmaPionTOF;
  Float_t nSigmaKaonTOF;
  Float_t nSigmaProtonTOF;
  Int_t fTrackCutFlag;

  /*** cut paramters */
  static Float_t fgMatchTrackDeltaX; // match track deltaX
  static Float_t fgMatchTrackDeltaZ; // match track deltaZ

  /*** tools ***/
  static TLorentzVector fgLorentzVector;
  //static AliTOFGeometry fgTOFGeometry;
  //  static AliTOFcalibHisto fgTOFcalibHisto;
  //  static Bool_t fgTOFcalibHistoFlag;
  static TH2F *hTOFtuned_th[AliPID::kSPECIES];

  Float_t fTimeZeroSigma; //!

  ClassDef(AliAnalysisPIDCascadeTrack, 1);
};

#endif /* ALIANALYSISPIDCASCADETRACK_H */
