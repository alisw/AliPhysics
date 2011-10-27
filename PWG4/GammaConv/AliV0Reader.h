#ifndef ALIV0READER_H
#define ALIV0READER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

////////////////////////////////////////////////
//--------------------------------------------- 
// Class used to do analysis on conversion pairs
//---------------------------------------------
////////////////////////////////////////////////

// --- ROOT system ---
#include "TObject.h" 
#include "AliMCEvent.h"   // for CF
#include "AliESDv0.h"
#include "AliESDEvent.h"
#include "AliKFParticle.h"
#include "TParticle.h"
#include "AliGammaConversionHistograms.h"
#include <vector>
#include "AliCFManager.h"
#include "AliGammaConversionBGHandler.h"
#include "AliESDpid.h"
#include "TF1.h"
#include "TRandom3.h"

class TClonesArray; 
class TFormula;
class Riostream;
class TChain;

//--- AliRoot system ---

class AliStack;
class AliMCEvent;       // for CF
class AliESDEvent; 
class AliMCEventHandler;
class AliESDInputHandler;
class AliESDVertex;
class AliLog;
class TChain;
class TChain;
class AliCFManager;   // for CF
class AliCFContainer;  // for CF
//class AliESDpid; // for dEdx cut based on nSigma to particle lines 
class AliESDtrackCuts; 
class TF1;
class TRandom3;

class AliV0Reader : public TObject {
	
 public: 
	
	
  // for CF
  enum{
    kStepGenerated = 0,
    kStepReconstructable = 1, 
    kStepGetOnFly = 2,
    kStepLikeSign = 3,
    kStepTPCRefit = 4,
    kStepKinks = 5,
    kStepdEdxElectronselection = 6,
    kStepdEdxPionrejection = 7,
    kStepNContributors = 8,
    kStepTPCPID = 9,
    kStepR = 10,
    kStepLine = 11,
    kStepZ = 12,
    kStepMinClsTPC = 13,
    kStepSinglePt= 14,  
    kStepNDF = 15,
    kStepChi2 = 16,
    kStepEta = 17,
    kStepPt = 18,
    kStepTrueGamma = 19
  };
	
  AliV0Reader();                                        //constructor
  AliV0Reader(const AliV0Reader & g);                   //copy constructor
  AliV0Reader & operator = (const AliV0Reader & g);     //assignment operator
  //  virtual ~AliV0Reader() {;}                            //virtual destructor
  virtual ~AliV0Reader();                            //virtual destructor
  /*
   *Initialize the reader
   */
  void Initialize();
  void SetInputAndMCEvent(AliVEvent* esd, AliMCEvent* mc) ;


  virtual void SetInputEvent(AliVEvent* const input)  {fESDEvent  = dynamic_cast<AliESDEvent*>(input);}
  virtual void SetMC(AliMCEvent* const mc)            {fMCEvent = mc;}

	
  // for CF
  void SetCFManager(AliCFManager * const io){fCFManager = io;};
  AliCFManager *GetCFManager() const {return fCFManager;}
	
	
	
	
  /*
   * Returns  AliESDEvent
   */			
  AliESDEvent* GetESDEvent() const{return fESDEvent;}	
	
  /*
   *Returns the number of v0s in the event, no cuts applied.
   */
  Int_t GetNumberOfV0s() const{return fESDEvent->GetNumberOfV0s();}
	
  /*
   *Returns the number of contributors to the vertex
   */
  //  Int_t GetNumberOfContributorsVtx() const{return fESDEvent->GetPrimaryVertex()->GetNContributors();}
  Int_t GetNumberOfContributorsVtx();
  
  /*
   * Check if there are any more good v0s left in the v0 stack
   * if so, fCurrent v0 is set to this v0 and can be retrieved
   * by GetCurrentV0 function.
   * returns kFALSE if there is no more good v0s in the v0 stack
   */
  Bool_t NextV0();
	
  /*
   * Returns the v0 at the given index, no checks are done on the v0. 
   */
  AliESDv0* GetV0(Int_t index);
	
  /*
   * Returns the current v0
   */
  AliESDv0* GetCurrentV0() const{return fCurrentV0;}
	
  /*
   * Returns the negative ESD track which belongs to fCurrentV0
   */
  //  AliESDtrack* GetNegativeESDTrack(){return fESDEvent->GetTrack(fCurrentV0->GetNindex());}
  AliESDtrack* GetNegativeESDTrack(){return fCurrentNegativeESDTrack;}
	
  /*
   * Returns the positive ESD track which belongs to fCurrentV0
   */
  //  AliESDtrack* GetPositiveESDTrack(){return fESDEvent->GetTrack(fCurrentV0->GetPindex());}
  AliESDtrack* GetPositiveESDTrack(){return fCurrentPositiveESDTrack;}
	
  /*
   * Returns the negative KF particle which belongs to fCurrentV0
   */
  AliKFParticle* GetNegativeKFParticle() const{return fCurrentNegativeKFParticle;}
	
  /*
   * Returns the positive KF particle which belongs to fCurrentV0
   */
  AliKFParticle* GetPositiveKFParticle() const{return fCurrentPositiveKFParticle;}
	
  /*
   * Returns the KFParticle object of the 2 tracks.
   */
  AliKFParticle* GetMotherCandidateKFCombination() const{return fCurrentMotherKFCandidate;}
	
  /*
   * Checks the probablity that the PID of the particle is what we want it to be.
   */
  Bool_t CheckPIDProbability(Double_t negProbCut, Double_t posProbCut);
	
  /*
   * Checks if the PID of the two particles are within our cuts.
   */
  void GetPIDProbability(Double_t &negPIDProb, Double_t &posPIDProb);

  /*
   * Checks if the PID of the two particles are within our cuts.
   */
  void GetPIDProbabilityMuonPion(Double_t &negPIDProb, Double_t &posPIDProb);
	
  /*
   *Get the negative MC TParticle from the stack 
   */
  TParticle * GetNegativeMCParticle() const{return fNegativeMCParticle;}
	
  /*
   *Get the positive MC TParticle from the stack 
   */
  TParticle * GetPositiveMCParticle() const{return fPositiveMCParticle;}
	
  /*
   *Get the mother MC TParticle from the stack 
   */
  TParticle * GetMotherMCParticle() const{return fMotherMCParticle;}
	
  /*
   * Flag to see if the v0 particles share the same mother
   */
  Bool_t HasSameMCMother();
	
	
  /*
   *Get the PID of the MC mother particle
   */
  Int_t GetMotherMCParticlePDGCode() const{return fMotherMCParticle->GetPdgCode();}
	
  /*
   *Get the MC stack 
   */
  AliStack* GetMCStack() const{return fMCStack;}
	
	
  /*
   * Setup  AliMCEventHandler
   */			
  //  AliMCEventHandler* GetMCTruth() const{return fMCTruth;}	// for CF
	
	
  /*
   *Get the MC stack 
   */
  AliMCEvent* GetMCEvent() const{return fMCEvent;}   // for CF
	
	
  /*
   *Get the magnetic field from the ESD event 
   */
  Double_t GetMagneticField() const{return fESDEvent->GetMagneticField();}
	
  /*
   *Get the primary vertex from the esd event
   */
  const AliESDVertex *GetPrimaryVertex() const {return fESDEvent->GetPrimaryVertex();}
	
  /*
   * Set the PID of the negative track
   */
  void SetNegativeTrackPID(Int_t negTrackPID){fNegativeTrackPID=negTrackPID;}
	
  /*
   * Set the PID of the positive track
   */
  void SetPositiveTrackPID(Int_t posTrackPID){fPositiveTrackPID=posTrackPID;}
	
  /*
   * Set the flag to use the kfparticle class. Will also disable the use of esd tracks
   */
  void UseKFParticle(){fUseKFParticle = kTRUE; fUseESDTrack = kFALSE;}
	
  /*
   *  Set the flag to use the esd track class. Will also disable the use of kf particles
   */
  void UseESDTrack(){fUseESDTrack = kTRUE; fUseKFParticle = kFALSE;}
	
  /*
   *  Set the flag to use improved vertex or not
   */
  void SetUseImprovedVertex(Bool_t useImprovedVertex){fUseImprovedVertex=useImprovedVertex;}
	
  /*
   * Return the number in the species array belonging to the negative or positive track pid.
   */
  Int_t GetSpeciesIndex(Int_t chargeOfTrack);
	
  /*
   * Return the x coordinate of the v0
   */
  Double_t GetX() const{return fCurrentXValue;}
	
  /*
   * Return the y coordinate of the v0
   */
  Double_t GetY() const{return fCurrentYValue;}
	
  /*
   * Return the Z coordinate of the v0
   */
  Double_t GetZ() const{return fCurrentZValue;}
	
  /*
   * Return the radius of the v0
   */
  Double_t GetXYRadius() const{return sqrt((Double_t)(fCurrentXValue*fCurrentXValue + fCurrentYValue*fCurrentYValue));}
	
  /*
   * Get the opening angle between the two tracks
   */
  Double_t GetOpeningAngle(){return fNegativeTrackLorentzVector->Angle(fPositiveTrackLorentzVector->Vect());}

  /*
   * Get the Cos Pointing angle between the two tracks
   */
  Double_t GetCosPointingAngle(){return fCurrentV0->GetV0CosineOfPointingAngle();}

  /*
   * Get the DCA between the two tracks
   */
  Double_t GetDcaDaughters(){return fCurrentV0->GetDcaV0Daughters();}

  /*
   * Get the Normalized DCA between the two tracks
   */
  Double_t GetNormDcaDistDaughters(){return fCurrentV0->GetDcaV0Daughters()/fCurrentV0->GetDistSigma();}

  /*
   * Get the Likelihood for a Conversion
   */
  Double_t GetLikelihoodAP(){return fCurrentV0->GetLikelihoodAP(0,0);}
      
  /*
   * Gets the Energy of the negative track.
   */
  Double_t GetNegativeTrackEnergy() const{return fCurrentNegativeKFParticle->E();}
	
  /*
   * Gets the Energy of the positive track.
   */
  Double_t GetPositiveTrackEnergy() const{return fCurrentPositiveKFParticle->E();}
	
  /*
   * Gets the Energy of the mother candidate.
   */
  Double_t GetMotherCandidateEnergy() const{return fCurrentMotherKFCandidate->E();}
	
  /*
   * Gets the Pt of the negative track.
   */
  Double_t GetNegativeTrackPt() const{return fNegativeTrackLorentzVector->Pt();}
	
  /*
   * Gets the Pt of the positive track.
   */
  Double_t GetPositiveTrackPt() const{return fPositiveTrackLorentzVector->Pt();}
	

  /*
   * Gets the Pt of the mother candidate.
   */
  Double_t GetMotherCandidatePt() const{return fMotherCandidateLorentzVector->Pt();}


  /*
   * Gets the P of the mother candidate.
   */
  Double_t GetMotherCandidateP() const{return fMotherCandidateLorentzVector->P();}
	

  /*
   * Gets the Eta of the negative track.
   */
  Double_t GetNegativeTrackEta() const{return fNegativeTrackLorentzVector->Eta();}
  /*
   * Gets the Eta of the positive track.
   */
  Double_t GetPositiveTrackEta() const{return fPositiveTrackLorentzVector->Eta();}
  /*
   * Gets the Eta of the mother candidate.
   */
  Double_t GetMotherCandidateEta() const{return fMotherCandidateLorentzVector->Eta();}
	
  /*
   * Gets the NDF of the mother candidate.
   */
  Double_t GetMotherCandidateNDF() const{return fCurrentMotherKFCandidate->GetNDF();}
	
  /*
   * Gets the Chi2 of the mother candidate.
   */
  Double_t GetMotherCandidateChi2() const{return fCurrentMotherKFCandidate->GetChi2();}
	
  /*
   * Gets the Mass of the mother candidate.
   */
  Double_t GetMotherCandidateMass() const{return fMotherCandidateKFMass;}
	
  /*
   * Gets the Width of the mother candidate.
   */
  Double_t GetMotherCandidateWidth() const{return fMotherCandidateKFWidth;}
	
  /*
   * Gets the Phi of the negative track.
   */
  Double_t GetNegativeTrackPhi() const;
	
  /*
   * Gets the Phi of the positive track.
   */
  Double_t GetPositiveTrackPhi() const;
	
  /*
   * Gets the Phi of the mother candidate.
   */
  Double_t GetMotherCandidatePhi() const;
	
  /*
   * Gets the Rapidity of the mother candidate.
   */
  Double_t GetMotherCandidateRapidity() const;
	

  /*
   * Gets the P of the negative track.
   */
  Double_t GetNegativeTrackP() const{return fNegativeTrackLorentzVector->P();}
	
  /*
   * Gets the P of the positive track.
   */
  Double_t GetPositiveTrackP() const{return fPositiveTrackLorentzVector->P();}

  /*
   * Gets the dE/dx in the TPC of the negative track.
   */
  Double_t GetNegativeTrackTPCdEdx() const{return fCurrentNegativeESDTrack->GetTPCsignal();}
	
  /*
   * Gets the dE/dx in the TPC of the positive track.
   */
  Double_t GetPositiveTrackTPCdEdx() const{return fCurrentPositiveESDTrack->GetTPCsignal();}

  /*
   * Gets the Number of the TPC clusters of the negative track.
   */
  Int_t GetNegativeTracknTPCClusters() const{return fCurrentNegativeESDTrack->GetNcls(1);}

  /*
   * Gets the Number of the TPC clusters of the positive track.
   */
  Int_t GetPositiveTracknTPCClusters() const{return fCurrentPositiveESDTrack->GetNcls(1);}

  /*
   * Get the TOFsignal for negative/positive track. RRnewTOF
   */
  Double_t GetNegativeTrackTOFsignal() const{return fCurrentNegativeESDTrack->GetTOFsignal();}
  Double_t GetPositiveTrackTOFsignal() const{return fCurrentPositiveESDTrack->GetTOFsignal();}	

  /*
   * Gets the Number of the TPC findable clusters of the negative track.
   */
  Int_t GetNegativeTracknTPCFClusters() const{return fCurrentNegativeESDTrack->GetTPCNclsF();}

  /*
   * Gets the Number of the TPC findable clusters of the positive track.
   */
  Int_t GetPositiveTracknTPCFClusters() const{return fCurrentPositiveESDTrack->GetTPCNclsF();}

  /*
   * Gets the Number of the ITS clusters of the negative track.
   */
  Int_t GetNegativeTracknITSClusters() const{return fCurrentNegativeESDTrack->GetNcls(0);}

  /*
   * Gets the Number of the ITS clusters of the positive track.
   */
  Int_t GetPositiveTracknITSClusters() const{return fCurrentPositiveESDTrack->GetNcls(0);}

  /*
   * Gets the chi2 of the TPC  negative track.
   */
  Double_t GetNegativeTrackTPCchi2() const{return fCurrentNegativeESDTrack->GetTPCchi2();}

  /*
   * Gets the chi2 of the TPC  the positive track.
   */
  Double_t GetPositiveTrackTPCchi2() const{return fCurrentPositiveESDTrack->GetTPCchi2();}
	
  /*
   * Update data which need to be updated every event.
   */
  void UpdateEventByEventData();

  /*
   * Gets the MaxRCut value.
   */
  Double_t GetMaxVertexZ() const{return fMaxVertexZ;}
	
  /*
   * Gets the MaxRCut value.
   */
  Double_t GetMaxRCut() const{return fMaxR;}

   /*
    * Gets the MinRCut value.
    */
   Double_t GetMinRCut() const{return fMinR;}
	
  /*
   * Gets the Eta cut value.
   */
  Double_t GetEtaCut() const{return fEtaCut;}

  /*
   * Gets the Eta cut value.
   */
  Double_t GetEtaCutMin() const{return fEtaCutMin;}

  /*
   * Gets the Rapidity Meson cut value.
   */
  Double_t GetRapidityMesonCut() const{return fRapidityMesonCut;}
	
  /*
   * Gets the Pt cut value.
   */
  Double_t GetPtCut() const{return fPtCut;}
  Double_t GetSinglePtCut() const{return fSinglePtCut;}	
	
  /*
   * Gets the MaxZCut value.
   */
  Double_t GetMaxZCut() const{return fMaxZ;}
	
	
  /*
   * Gets the MinClsTPC value.
   */
  Double_t GetMinClsTPCCut() const{return fMinClsTPC;}
	
  /*
   * Gets the MinClsTPC value.
   */
  Double_t GetMinClsTPCCutToF() const{return fMinClsTPCToF;}
	


  /*
   * Gets the line cut values.
   */
  Double_t GetLineCutZRSlope() const{return fLineCutZRSlope;}
  Double_t GetLineCutZRSlopeMin() const{return fLineCutZRSlopeMin;}
  Double_t GetLineCutZValue() const{return fLineCutZValue;}
  Double_t GetLineCutZValueMin() const{return fLineCutZValueMin;}	
  /*
   * Gets the Chi2 cut value for the conversions.
   */
  Double_t GetChi2CutConversion() const{return fChi2CutConversion;}
	
  /*
   * Gets the Chi2 cut value for the mesons.
   */
  Double_t GetChi2CutMeson() const{return fChi2CutMeson;}
	
  /*
   * Gets the alpha cut value for the mesons.
   */
  Double_t GetAlphaCutMeson() const{return fAlphaCutMeson;}

  /*
   * Gets the Minimum alpha cut value for the mesons.
   */
  Double_t GetAlphaMinCutMeson() const{return fAlphaMinCutMeson;}

  Double_t GetPositiveTrackLength() const{return fCurrentPositiveESDTrack->GetIntegratedLength();}
  Double_t GetNegativeTrackLength() const{return fCurrentNegativeESDTrack->GetIntegratedLength();}
	
  Double_t GetPositiveNTPCClusters() const{return fCurrentPositiveESDTrack->GetTPCNcls();}
  Double_t GetNegativeNTPCClusters() const{return fCurrentNegativeESDTrack->GetTPCNcls();}
	
  /*
   * Sets the MaxVertexZ value.
   */
  void SetMaxVertexZ(Double_t maxVertexZ){fMaxVertexZ=maxVertexZ;}

  /*
   * Sets the MaxRCut value.
   */
  void SetMaxRCut(Double_t maxR){fMaxR=maxR;}
  /*	
   * Sets the MinRCut value.
   */
  void SetMinRCut(Double_t minR){fMinR=minR;}

  /*
   * Sets the EtaCut value.
   */
  void SetEtaCut(Double_t etaCut){fEtaCut=etaCut;}

  /*
   * Sets the EtaCutMin value.
   */
  void SetEtaCutMin(Double_t etaCutMin){fEtaCutMin=etaCutMin;}


  /*
   * Sets the Rapidity Meson Cut value.
   */
  void SetRapidityMesonCut(Double_t RapidityMesonCut){fRapidityMesonCut=RapidityMesonCut;}
	
  /*
   * Sets the PtCut value.
   */
  void SetPtCut(Double_t ptCut){fPtCut=ptCut;}
	
  /*
   * Sets the PtCut value.
   */
  void SetSinglePtCut(Double_t singleptCut){fSinglePtCut=singleptCut;}
	
    
  /*
   * Sets the MaxZCut value.
   */
  void SetMaxZCut(Double_t maxZ){fMaxZ=maxZ;}
	
  /*
   * Sets the MinClsTPC value.
   */
  void SetMinClsTPCCut(Double_t minClsTPC){fMinClsTPC=minClsTPC;}

  /*
   * Sets the MinClsTPC value.
   */
  void SetMinClsTPCCutToF(Double_t minClsTPCToF){fMinClsTPCToF=minClsTPCToF;}
	

  /*
   * Sets the LineCut values.
   */
  void SetLineCutZRSlope(Double_t LineCutZRSlope){fLineCutZRSlope=LineCutZRSlope;}
  void SetLineCutZValue(Double_t LineCutZValue){fLineCutZValue=LineCutZValue;}
	
  void SetLineCutZRSlopeMin(Double_t LineCutZRSlopeMin){fLineCutZRSlopeMin=LineCutZRSlopeMin;}
  void SetLineCutZValueMin(Double_t LineCutZValueMin){fLineCutZValueMin=LineCutZValueMin;}
		
  /*
   * Sets the Chi2Cut value for conversions.
   */
  void SetChi2CutConversion(Double_t chi2){fChi2CutConversion=chi2;}
	
  /*
   * Sets the Chi2Cut for the mesons.
   */
  void SetChi2CutMeson(Double_t chi2){fChi2CutMeson=chi2;}
	
  /*
   * Sets the AlphaCut for the mesons.
   */
  void SetAlphaCutMeson(Double_t alpha){fAlphaCutMeson=alpha;}
	

  /*
   * Sets the AlphaCut for the mesons.
   */
  void SetAlphaMinCutMeson(Double_t alpha){fAlphaMinCutMeson=alpha;}


  /*
   * Sets the XVertexCut value.
   */
  void SetXVertexCut(Double_t xVtx){fCurrentXValue=xVtx;}
	
  /*
   * Sets the YVertexCut value.
   */
  void SetYVertexCut(Double_t yVtx){fCurrentYValue=yVtx;}
	
  /*
   * Sets the ZVertexCut value.
   */
  void SetZVertexCut(Double_t zVtx){fCurrentZValue=zVtx;}
	
  /*
   * Sets the PIDProbabilityCut value for track particles.
   */
  void SetPIDProbability(Double_t pidProb){fPIDProbabilityCutPositiveParticle=pidProb; fPIDProbabilityCutNegativeParticle=pidProb;}
	
  /*
   * Sets the PIDProbability cut value for the negative track.
   */
  void SetPIDProbabilityNegativeParticle(Double_t pidProb){fPIDProbabilityCutNegativeParticle=pidProb;}
	
  /*
   * Sets the PIDProbability cut value for the positive track.
   */
  void SetPIDProbabilityPositiveParticle(Double_t pidProb){fPIDProbabilityCutPositiveParticle=pidProb;}

  /*
   * Sets the PIDnSigmaAboveElectron cut value for the tracks.
   */
  void SetPIDnSigmaAboveElectronLine(Double_t nSigmaAbove){fPIDnSigmaAboveElectronLine=nSigmaAbove;}
  void SetTofPIDnSigmaAboveElectronLine(Double_t nTofSigmaAbove){fTofPIDnSigmaAboveElectronLine=nTofSigmaAbove;} // RRnewTOF
	
  /*
   * Sets the PIDnSigmaBelowElectron cut value for the tracks.
   */
  void SetPIDnSigmaBelowElectronLine(Double_t nSigmaBelow){fPIDnSigmaBelowElectronLine=nSigmaBelow;}
  void SetTofPIDnSigmaBelowElectronLine(Double_t nTofSigmaBelow){fTofPIDnSigmaBelowElectronLine=nTofSigmaBelow;} // RRnewTOF
	
  /*
   * Sets the PIDnSigmaAbovePion cut value for the tracks.
   */
  void SetPIDnSigmaAbovePionLine(Double_t nSigmaAbovePion){fPIDnSigmaAbovePionLine=nSigmaAbovePion;}

  /*
   * Sets the PIDnSigmaAbovePion cut value for the tracks.
   */
  void SetPIDnSigmaAbovePionLineHighPt(Double_t nSigmaAbovePionHighPt){fPIDnSigmaAbovePionLineHighPt=nSigmaAbovePionHighPt;}

  /*
   * Sets the PIDMinPnSigmaAbovePion cut value for the tracks.
   */
  void SetPIDMinPnSigmaAbovePionLine(Double_t MinPnSigmaAbovePion){fPIDMinPnSigmaAbovePionLine=MinPnSigmaAbovePion;}

 /*
   * Sets the PIDMinPnSigmaAbovePion cut value for the tracks.
   */
  void SetPIDMaxPnSigmaAbovePionLine(Double_t MaxPnSigmaAbovePion){fPIDMaxPnSigmaAbovePionLine=MaxPnSigmaAbovePion;}

  /*
   * Sets the SigmaMassCut value.
   */
  void SetSigmaMass(Double_t sigmaMass){fNSigmaMass=sigmaMass;}
	
  /*
   * Sets the flag to enable/disable the usage of MC information. 
   */
  void SetDoMCTruth(Bool_t doMC){fDoMC = doMC;}

  /*
   * Sets the flag to enable/disable the usage of MC information. 
   */
  Bool_t GetDoMCTruth() const {return fDoMC;}
	
  /*
   * Sets the flag to enable/disable the cut dedx N sigma 
   */

  void SetDodEdxSigmaCut( Bool_t dodEdxSigmaCut){fDodEdxSigmaCut=dodEdxSigmaCut;}
  void SetDoTOFsigmaCut( Bool_t doTOFsigmaCut){fDoTOFsigmaCut=doTOFsigmaCut;} //RRnewTOF
  void SetDoPhotonAsymmetryCut( Bool_t doPhotonAsymmetryCut){fDoPhotonAsymmetryCut=doPhotonAsymmetryCut;}

  void SetMinPPhotonAsymmetryCut(Double_t minPPhotonAsymmetryCut){fMinPPhotonAsymmetryCut=minPPhotonAsymmetryCut;}
  void SetMinPhotonAsymmetry(Double_t minPhotonAsymmetry){fMinPhotonAsymmetry=minPhotonAsymmetry;}
  /*
   * Sets the flag to enable/disable the cut dedx N sigma for Kaon Rejection at low p 
   */
  void SetDoKaonRejectionLowP( Bool_t doKaonRejectionLowP){fDoKaonRejectionLowP=doKaonRejectionLowP;}
  /*
   * Sets the flag to enable/disable the cut dedx N sigma for Proton Rejection at low p 
   */
  void SetDoProtonRejectionLowP( Bool_t doProtonRejectionLowP){fDoProtonRejectionLowP=doProtonRejectionLowP;}

  /*
   * Sets the flag to enable/disable the cut dedx N sigma for Pion Rejection at low p 
   */
  void SetDoPionRejectionLowP( Bool_t doPionRejectionLowP){fDoPionRejectionLowP=doPionRejectionLowP;}

  /*
   * Sets the PIDMinPnSigmaAroundKaon cut value for the tracks.
   */
  void SetPIDnSigmaAtLowPAroundKaonLine(Double_t nSigmaAtLowPAroundKaon){fPIDnSigmaAtLowPAroundKaonLine =nSigmaAtLowPAroundKaon;}

  /*
   * Sets the PIDMinPnSigmaAroundProton cut value for the tracks.
   */
  void SetPIDnSigmaAtLowPAroundProtonLine(Double_t nSigmaAtLowPAroundProton){fPIDnSigmaAtLowPAroundProtonLine =nSigmaAtLowPAroundProton;}

  /*
   * Sets the PIDMinPnSigmaAroundPion cut value for the tracks.
   */
  void SetPIDnSigmaAtLowPAroundPionLine(Double_t nSigmaAtLowPAroundPion){fPIDnSigmaAtLowPAroundPionLine =nSigmaAtLowPAroundPion;}

  /*
   * Sets the PIDMinPnSigmaAbovePion cut value for the tracks.
   */
  void SetPIDMinPKaonRejectionLowP(Double_t PIDMinPKaonRejectionLowP ){fPIDMinPKaonRejectionLowP=PIDMinPKaonRejectionLowP;}

  /*
   * Sets the PIDMinPnSigmaAbovePion cut value for the tracks.
   */
  void SetPIDMinPProtonRejectionLowP(Double_t PIDMinPProtonRejectionLowP ){fPIDMinPProtonRejectionLowP=PIDMinPProtonRejectionLowP;}
  /*
   * Sets the PIDMinPnSigmaAbovePion cut value for the tracks.
   */
  void SetPIDMinPPionRejectionLowP(Double_t PIDMinPPionRejectionLowP ){fPIDMinPPionRejectionLowP=PIDMinPPionRejectionLowP;}

  /*
   *Set if we want to use Gamma Selection based on Qt from Armenteros
   */
  void SetDoQtGammaSelection(Bool_t doQtGammaSelection){fDoQtGammaSelection=doQtGammaSelection;}
  void SetDoHighPtQtGammaSelection(Bool_t doHighPtQtGammaSelection){fDoHighPtQtGammaSelection=doHighPtQtGammaSelection;} // RRnew
  /*
   * Sets the MaxQtCut value.
   */
  void SetQtMax(Double_t qtMax){fQtMax=qtMax;}
  void SetHighPtQtMax(Double_t qtMaxHighPt){fHighPtQtMax=qtMaxHighPt;} // RRnew
  void SetPtBorderForQt(Double_t ptBorderForQt){fPtBorderForQt=ptBorderForQt;} // RRnew

  /*
   * Updates the V0 information of the current V0.
   */
  Bool_t UpdateV0Information();
	
  /*
   * Resets the V0 index.
   */
  void ResetV0IndexNumber(){fCurrentV0IndexNumber=0;}
  

  /*
   * Returns number of good v0s in the event
   */
  Int_t GetNGoodV0s() const {return fNumberOfGoodV0s;}

  /*
   * Sets the histograms.
   */
  void SetHistograms(AliGammaConversionHistograms * const histograms){fHistograms=histograms;}
	
  /*
   * Check for primary vertex.
   */
  Bool_t CheckForPrimaryVertex();
	
  /*
   * Check for primary vertex Z.
   */
  Bool_t CheckForPrimaryVertexZ();

  /*
   * Gets a vector of good v0s.
   */
  TClonesArray* GetCurrentEventGoodV0s() const{return fCurrentEventGoodV0s;}
	
  /*
   * Gets the vector of previous events v0s (for bacground analysis)
   */
  AliGammaConversionKFVector* GetBGGoodV0s(Int_t event) const;
  //  vector<AliKFParticle> GetPreviousEventGoodV0s() const{return fPreviousEventGoodV0s;}

  void SetUseOwnXYZCalculation(Bool_t flag){fUseOwnXYZCalculation=flag;}

  void SetUseConstructGamma(Bool_t flag){fUseConstructGamma=flag;}

  //   Bool_t GetHelixCenter(const AliExternalTrackParam *track, Double_t b,Int_t charge, Double_t center[2]);
  Bool_t GetHelixCenter(AliESDtrack* track, Double_t b,Int_t charge, Double_t center[2]); 

  Bool_t GetConvPosXY(AliESDtrack* ptrack,AliESDtrack* ntrack, Double_t b, Double_t convpos[2]);
	
  Double_t GetConvPosZ(AliESDtrack* ptrack,AliESDtrack* ntrack, Double_t b);

	void GetConversionPoint(const AliExternalTrackParam *pparam,const AliExternalTrackParam *nparam,Double_t convpos[3]);

	Bool_t GetArmenterosQtAlfa(const AliKFParticle * posKFparticle, const AliKFParticle * negKFparticle, const AliKFParticle * gamKFparticle,Double_t armenterosQtAlfa[2]);
	Bool_t GetArmenterosQtAlfa(const TParticle * posKFparticle, const TParticle * negKFparticle, const AliKFParticle * gamKFparticle,Double_t armenterosQtAlfa[2]);
	Bool_t GetArmenterosQtAlfa(const TParticle * posKFparticle, const TParticle * negKFparticle, const TParticle * gamKFparticle,Double_t armenterosQtAlfa[2]);
  
  void SetDoCF(Bool_t flag){fDoCF = flag;}

  Bool_t CheckV0FinderStatus(Int_t index);

  void SetOnFlyFlag(Bool_t flag){fUseOnFlyV0Finder = flag;}
  
  Int_t GetNBGEvents(){return fBGEventHandler->GetNBGEvents();}

  void SetCalculateBackground(Bool_t flag){fCalculateBackground=flag;}

  AliGammaConversionBGHandler* GetBGHandler() const {return fBGEventHandler;}

  Double_t GetVertexZ(){return fESDEvent->GetPrimaryVertex()->GetZ();}

  Int_t GetMultiplicity(){return CountESDTracks();}

  void SetESDtrackCuts(AliESDtrackCuts * const trackCuts){fEsdTrackCuts = trackCuts;}

  void SetNEventsForBG(Int_t nev){fNEventsForBGCalculation = nev;}

  Int_t CountESDTracks();

  Int_t GetCurrentV0IndexNumber() const {return fCurrentV0IndexNumber;}

  Bool_t CheckIfPi0IsMother(Int_t label);
  Bool_t CheckIfEtaIsMother(Int_t label);

  static void InitESDpid(Int_t type=0);
  static void SetESDpid(AliESDpid * const pid) {fgESDpid=pid;}
  static AliESDpid* GetESDpid() {return fgESDpid;}

  void SetUseChargedTracksMultiplicityForBG(Bool_t flag){fUseChargedTrackMultiplicityForBG = flag;}
  void SetIsHeavyIon(Int_t isHeavyIon) {fIsHeavyIon=isHeavyIon;}
  Int_t GetIsHeavyIon() const { return fIsHeavyIon;}
 
  Int_t GetPindex(Int_t i) {return fV0Pindex.at(i);}
  Int_t GetNindex(Int_t i) {return fV0Nindex.at(i);}

  void ResetNGoodV0s(){fNumberOfGoodV0s=0;}
  Int_t GetFirstTPCRow(Double_t radius);

  void SetUseCorrectedTPCClsInfo(Bool_t flag){fUseCorrectedTPCClsInfo = flag;}
  Bool_t GetUseCorrectedTPCClsInfo() const {return fUseCorrectedTPCClsInfo;}

  void SetUseMCPSmearing(Int_t useMCPSmearing) {fUseMCPSmearing=useMCPSmearing;}
  void SetPBremSmearing(Double_t pBremSmearing){fPBremSmearing=pBremSmearing;}
  void SetPSigSmearing(Double_t pSigSmearing){fPSigSmearing=pSigSmearing;}
  void SetPSigSmearingCte(Double_t pSigSmearingCte){fPSigSmearingCte=pSigSmearingCte;}
  void SmearKFParticle(AliKFParticle * kfParticle);

 private:
  AliStack * fMCStack;           // pointer to MonteCarlo particle stack 
  //  AliMCEventHandler* fMCTruth;   // for CF    pointer to the MC object
  AliMCEvent *fMCEvent;			//  for CF      pointer to MC event
  TChain * fChain;               // pointer to the TChain
	
  //  AliESDInputHandler* fESDHandler;      //! pointer to esd object
  AliESDEvent *fESDEvent;               //! pointer to esd object
	
	
  // for CF
  AliCFManager *fCFManager; // pointer to the cf manager
  //  AliCFContainer *container;
	
  // for dEdx cut based on nSigma to a particle line
  //AliESDpid * fESDpid; // esd pid
	
  AliGammaConversionHistograms *fHistograms; // pointer to histogram handling class
	
  Int_t fCurrentV0IndexNumber;  // the current v0 index number
  AliESDv0 * fCurrentV0;                //! pointer to the current v0
  AliKFParticle * fCurrentNegativeKFParticle;  //! pointer to the negative KF particle
  AliKFParticle * fCurrentPositiveKFParticle;  //! pointer to the positive KF particle
  AliKFParticle * fCurrentMotherKFCandidate;   //! pointer to the positive KF particle
	
  AliESDtrack * fCurrentNegativeESDTrack;      //! pointer to the negative ESD track
  AliESDtrack * fCurrentPositiveESDTrack;      //! pointer to the positive ESD track
	
  TLorentzVector * fNegativeTrackLorentzVector; //! pointer to the negative Track Lorentz Vector
  TLorentzVector * fPositiveTrackLorentzVector; //! pointer to the positive Track Lorentz Vector
  TLorentzVector * fMotherCandidateLorentzVector;   //! pointer to the mother candidate Track Lorentz Vector
	
  Double_t fCurrentXValue;   // current x value
  Double_t fCurrentYValue;   // current y value
  Double_t fCurrentZValue;   // current z value
	
  Int_t fPositiveTrackPID;   // positive track pid
  Int_t fNegativeTrackPID;   // negative track pid
	
  TParticle *fNegativeMCParticle;      //!
  TParticle *fPositiveMCParticle;      //!
  TParticle *fMotherMCParticle;        //!
	
  Double_t fMotherCandidateKFMass;   // mass of mother candidate KF particle
  Double_t fMotherCandidateKFWidth;  // width of mother candidate KF particle
	
  Bool_t fUseKFParticle;   // flag 
  Bool_t fUseESDTrack;     // flag 
  Bool_t fDoMC;            // flag 

  //Event Cuts
  Double_t fMaxVertexZ; // max z vertex cut
  //cuts
  Double_t fMaxR; //r cut
  Double_t fMinR; //r cut
  Double_t fEtaCut; //eta cut
  Double_t fEtaCutMin; //eta cut
  Double_t fRapidityMesonCut; //rapidity for meson cut
  Double_t fPtCut; // pt cut
  Double_t fSinglePtCut; // pt cut for electron/positron
  Double_t fMaxZ; //z cut
  Double_t fMinClsTPC; // minimum clusters in the TPC
  Double_t fMinClsTPCToF; // minimum clusters to findable clusters
  Double_t fLineCutZRSlope; //linecut
  Double_t fLineCutZValue; //linecut
  Double_t fLineCutZRSlopeMin; //linecut
  Double_t fLineCutZValueMin; //linecut
  Double_t fChi2CutConversion; //chi2cut
  Double_t fChi2CutMeson;  //chi2cut
  Double_t fAlphaCutMeson;  //alphacut
  Double_t fAlphaMinCutMeson;  //alphacut
  Double_t fPIDProbabilityCutNegativeParticle; //pid cut
  Double_t fPIDProbabilityCutPositiveParticle; //pid cut
  Bool_t   fDodEdxSigmaCut; // flag to use the dEdxCut based on sigmas
  Bool_t   fDoTOFsigmaCut; // flag to use TOF pid cut RRnewTOF
  Double_t fPIDnSigmaAboveElectronLine; // sigma cut
  Double_t fPIDnSigmaBelowElectronLine; // sigma cut
  Double_t fTofPIDnSigmaAboveElectronLine; // sigma cut RRnewTOF
  Double_t fTofPIDnSigmaBelowElectronLine; // sigma cut RRnewTOF 
  Double_t fPIDnSigmaAbovePionLine;     // sigma cut
  Double_t fPIDnSigmaAbovePionLineHighPt;     // sigma cut
  Double_t fPIDMinPnSigmaAbovePionLine; // sigma cut
  Double_t fPIDMaxPnSigmaAbovePionLine; // sigma cut
  Double_t fDoKaonRejectionLowP;   // Kaon rejection at low p
  Double_t fDoProtonRejectionLowP; // Proton rejection at low p
  Double_t fDoPionRejectionLowP;   // Pion rejection at low p
  Double_t fPIDnSigmaAtLowPAroundKaonLine; // sigma cut
  Double_t fPIDnSigmaAtLowPAroundProtonLine; // sigma cut
  Double_t fPIDnSigmaAtLowPAroundPionLine; // sigma cut
  Double_t fPIDMinPKaonRejectionLowP; // Momentum limit to apply kaon rejection
  Double_t fPIDMinPProtonRejectionLowP; // Momentum limit to apply proton rejection
  Double_t fPIDMinPPionRejectionLowP; // Momentum limit to apply proton rejection
  Bool_t   fDoQtGammaSelection; // Select gammas using qtMax
  Bool_t   fDoHighPtQtGammaSelection; // RRnew Select gammas using qtMax for high pT
  Double_t fQtMax; // Maximum Qt from Armenteros to select Gammas
  Double_t fHighPtQtMax; // RRnew Maximum Qt for High pT from Armenteros to select Gammas
  Double_t fPtBorderForQt; // RRnew 
  Double_t fXVertexCut; //vertex cut
  Double_t fYVertexCut; //vertex cut
  Double_t fZVertexCut; // vertexcut
	
  Double_t fNSigmaMass; //nsigma cut
	
  Bool_t fUseImprovedVertex; //flag

  Bool_t fUseOwnXYZCalculation; //flag that determines if we use our own calculation of xyz (markus)

  Bool_t fUseConstructGamma; //flag that determines if we use ConstructGamma method from AliKF

  Bool_t fDoCF; //flag

  Bool_t fUseEtaMinCut; //flag

  Bool_t fUseOnFlyV0Finder; //flag

  Bool_t fUpdateV0AlreadyCalled; //flag
	
  TClonesArray* fCurrentEventGoodV0s; //vector of good v0s
 
  vector<Int_t> fV0Pindex; // index of positive track belonging to a V0
  vector<Int_t> fV0Nindex; // index of positive track belonging to a V0
  //  vector<AliKFParticle> fPreviousEventGoodV0s; // vector of good v0s from prevous events

  Bool_t fCalculateBackground; //flag
  AliGammaConversionBGHandler *fBGEventHandler; // background handler
  Bool_t fBGEventInitialized; //flag
	
  AliESDtrackCuts *fEsdTrackCuts; // track cuts
  Int_t fNumberOfESDTracks; //track counter

  static AliESDpid* fgESDpid;                 // ESD pid object

  Int_t fNEventsForBGCalculation; // Number of events used for background calculation
  
  Bool_t fUseChargedTrackMultiplicityForBG; // flag
  Int_t fNumberOfGoodV0s; // number of good V0s
  Int_t fIsHeavyIon; // flag
  Bool_t fUseCorrectedTPCClsInfo;
  Int_t fUseMCPSmearing;
  Double_t fPBremSmearing;
  Double_t fPSigSmearing;
  Double_t fPSigSmearingCte;
  TRandom3 fRandom;
  TF1 * fBrem;	
  Bool_t   fDoPhotonAsymmetryCut; // flag to use the PhotonAsymetryCut
  Double_t fMinPPhotonAsymmetryCut;
  Double_t fMinPhotonAsymmetry;
  ClassDef(AliV0Reader,22) // RRnewTOF
};

inline void AliV0Reader::InitESDpid(Int_t type)
{
  //
  // initialize PID parameters
  // type=0 is simulation
  // type=1 is data

  if (!fgESDpid) fgESDpid=new AliESDpid;
  Double_t alephParameters[5];
  // simulation
  alephParameters[0] = 2.15898e+00/50.;
  alephParameters[1] = 1.75295e+01;
  alephParameters[2] = 3.40030e-09;
  alephParameters[3] = 1.96178e+00;
  alephParameters[4] = 3.91720e+00;
  fgESDpid->GetTOFResponse().SetTimeResolution(80.);

  // data
  if (type==1){
    alephParameters[0] = 0.0283086/0.97;
    alephParameters[1] = 2.63394e+01;
    alephParameters[2] = 5.04114e-11;
    alephParameters[3] = 2.12543e+00;
    alephParameters[4] = 4.88663e+00;
    fgESDpid->GetTOFResponse().SetTimeResolution(130.);
    fgESDpid->GetTPCResponse().SetMip(50.);
  }

  fgESDpid->GetTPCResponse().SetBetheBlochParameters(
    alephParameters[0],alephParameters[1],alephParameters[2],
    alephParameters[3],alephParameters[4]);

  fgESDpid->GetTPCResponse().SetSigma(3.79301e-03, 2.21280e+04);

}

#endif





