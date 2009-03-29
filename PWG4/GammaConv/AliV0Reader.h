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
#include "AliESDv0.h"
#include "AliESDEvent.h"
#include "AliKFParticle.h"
#include "TParticle.h"
#include "AliGammaConversionHistograms.h"
#include <vector>

class TClonesArray; 
class TFormula;
class Riostream;
class TChain;
//--- AliRoot system ---

class AliStack;
class AliESDEvent; 
class AliMCEventHandler;
class AliESDInputHandler;
class AliESDVertex;
class AliLog;
class TChain;
class TChain;



class AliV0Reader : public TObject {

 public: 

  AliV0Reader();                                        //constructor
  AliV0Reader(const AliV0Reader & g);                   //copy constructor
  AliV0Reader & operator = (const AliV0Reader & g);     //assignment operator
  virtual ~AliV0Reader() {;}                            //virtual destructor
  /*
   *Initialize the reader
   */
  void Initialize();

  /*
   *Returns the number of v0s in the event, no cuts applied.
   */
  Int_t GetNumberOfV0s() const{return fESDEvent->GetNumberOfV0s();}

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
  AliESDtrack* GetNegativeESDTrack(){return fESDEvent->GetTrack(fCurrentV0->GetNindex());}

  /*
   * Returns the positive ESD track which belongs to fCurrentV0
   */
  AliESDtrack* GetPositiveESDTrack(){return fESDEvent->GetTrack(fCurrentV0->GetPindex());}

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
  Int_t GetMotherMCParticlePDGCode() const{if(fMotherMCParticle != NULL){ cout<<"MCParticle exists"<<endl;} return fMotherMCParticle->GetPdgCode();}

  /*
   *Get the MC stack 
   */
  AliStack* GetMCStack() const{return fMCStack;}

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
   * Update data which need to be updated every event.
   */
  void UpdateEventByEventData();
  
  /*
   * Gets the MaxRCut value.
   */
  Double_t GetMaxRCut() const{return fMaxR;}

  /*
   * Gets the Eta cut value.
   */
  Double_t GetEtaCut() const{return fEtaCut;}

  /*
   * Gets the Pt cut value.
   */
  Double_t GetPtCut() const{return fPtCut;}

  /*
   * Gets the line cut values.
   */
  Double_t GetLineCutZRSlope() const{return fLineCutZRSlope;}
  Double_t GetLineCutZValue() const{return fLineCutZValue;}
  
  /*
   * Gets the Chi2 cut value for the conversions.
   */
  Double_t GetChi2CutConversion() const{return fChi2CutConversion;}

  /*
   * Gets the Chi2 cut value for the mesons.
   */
  Double_t GetChi2CutMeson() const{return fChi2CutMeson;}

  Double_t GetPositiveTrackLength() const{return fCurrentPositiveESDTrack->GetIntegratedLength();}
  Double_t GetNegativeTrackLength() const{return fCurrentNegativeESDTrack->GetIntegratedLength();}

  Double_t GetPositiveNTPCClusters() const{return fCurrentPositiveESDTrack->GetTPCNcls();}
  Double_t GetNegativeNTPCClusters() const{return fCurrentNegativeESDTrack->GetTPCNcls();}

  /*
   * Sets the MaxRCut value.
   */
  void SetMaxRCut(Double_t maxR){fMaxR=maxR;}

  /*
   * Sets the EtaCut value.
   */
  void SetEtaCut(Double_t etaCut){fEtaCut=etaCut;}

  /*
   * Sets the PtCut value.
   */
  void SetPtCut(Double_t ptCut){fPtCut=ptCut;}

  /*
   * Sets the LineCut values.
   */
  void SetLineCutZRSlope(Double_t LineCutZRSlope){fLineCutZRSlope=LineCutZRSlope;}
  void SetLineCutZValue(Double_t LineCutZValue){fLineCutZValue=LineCutZValue;}

  /*
   * Sets the Chi2Cut value for conversions.
   */
  void SetChi2CutConversion(Double_t chi2){fChi2CutConversion=chi2;}

  /*
   * Sets the Chi2Cut for the mesons.
   */
  void SetChi2CutMeson(Double_t chi2){fChi2CutMeson=chi2;}
  
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
   * Sets the SigmaMassCut value.
   */
  void SetSigmaMass(Double_t sigmaMass){fNSigmaMass=sigmaMass;}

  /*
   * Sets the flag to enable/disable the usage of MC information. 
   */
  void SetDoMCTruth(Bool_t doMC){fDoMC = doMC;}

  /*
   * Updates the V0 information of the current V0.
   */
  Bool_t UpdateV0Information();

  /*
   * Resets the V0 index.
   */
  void ResetV0IndexNumber(){fCurrentV0IndexNumber=0;}

  /*
   * Sets the histograms.
   */
  void SetHistograms(AliGammaConversionHistograms *histograms){fHistograms=histograms;}

  /*
   * Check for primary vertex.
   */
  Bool_t CheckForPrimaryVertex();

  /*
   * Gets a vector of good v0s.
   */
  vector<AliKFParticle> GetCurrentEventGoodV0s() const{return fCurrentEventGoodV0s;}

  /*
   * Gets the vector of previous events v0s (for bacground analysis)
   */
  vector<AliKFParticle> GetPreviousEventGoodV0s() const{return fPreviousEventGoodV0s;}

 private:
  AliStack * fMCStack;           // pointer to MonteCarlo particle stack 
  AliMCEventHandler* fMCTruth;   // pointer to the MC event handler
  TChain * fChain;               // pointer to the TChain
  
  AliESDInputHandler* fESDHandler;      //! pointer to esd object
  AliESDEvent *fESDEvent;               //! pointer to esd object

  AliGammaConversionHistograms *fHistograms;
  
  Int_t fCurrentV0IndexNumber;
  AliESDv0 * fCurrentV0;                //! pointer to the current v0
  AliKFParticle * fCurrentNegativeKFParticle;  //! pointer to the negative KF particle
  AliKFParticle * fCurrentPositiveKFParticle;  //! pointer to the positive KF particle
  AliKFParticle * fCurrentMotherKFCandidate;   //! pointer to the positive KF particle

  AliESDtrack * fCurrentNegativeESDTrack;      //! pointer to the negative ESD track
  AliESDtrack * fCurrentPositiveESDTrack;      //! pointer to the positive ESD track
 
  TLorentzVector * fNegativeTrackLorentzVector; //! pointer to the negative Track Lorentz Vector
  TLorentzVector * fPositiveTrackLorentzVector; //! pointer to the positive Track Lorentz Vector
  TLorentzVector * fMotherCandidateLorentzVector;   //! pointer to the mother candidate Track Lorentz Vector

  Double_t fCurrentXValue;
  Double_t fCurrentYValue;
  Double_t fCurrentZValue;

  Int_t fPositiveTrackPID;
  Int_t fNegativeTrackPID;

  TParticle *fNegativeMCParticle;      //!
  TParticle *fPositiveMCParticle;      //!
  TParticle *fMotherMCParticle;        //!

  Double_t fMotherCandidateKFMass;
  Double_t fMotherCandidateKFWidth;

  Bool_t fUseKFParticle;
  Bool_t fUseESDTrack;
  Bool_t fDoMC;

  //cuts
  Double_t fMaxR;
  Double_t fEtaCut;
  Double_t fPtCut;
  Double_t fLineCutZRSlope;
  Double_t fLineCutZValue;
  Double_t fChi2CutConversion;
  Double_t fChi2CutMeson;
  Double_t fPIDProbabilityCutNegativeParticle;
  Double_t fPIDProbabilityCutPositiveParticle;
  Double_t fXVertexCut;
  Double_t fYVertexCut;
  Double_t fZVertexCut;
  
  Double_t fNSigmaMass;
  
  Bool_t fUseImprovedVertex;
  
  vector<AliKFParticle> fCurrentEventGoodV0s;
  vector<AliKFParticle> fPreviousEventGoodV0s;

  ClassDef(AliV0Reader,1)
};


#endif



