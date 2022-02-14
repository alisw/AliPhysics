#ifndef ALIANALYSISTASKEFFMATRIX_H
#define ALIANALYSISTASKEFFMATRIX_H

/*
 * Efficiency Matrix Analysis
 * author: Roberto Preghenella (R+)
 * email:  preghenella@bo.infn.it
 *
 */

#include "AliAnalysisTaskSE.h"
#include "THnSparse.h"

class AliVEvent;
class AliVVertex;
class AliVTrack;
class AliVParticle;
class TList;
class TClonesArray;
class AliESDtrackCuts;

class AliAnalysisTaskEffMatrix :
public AliAnalysisTaskSE
{
    
public:
    
    AliAnalysisTaskEffMatrix(const Char_t *name = "EffMatrix"); // default constructor
    virtual ~AliAnalysisTaskEffMatrix(); // default destructor
    
    void UserCreateOutputObjects();
    void UserExec(Option_t *option);
    
    /* setters */
    void SetAODfilterBit(Int_t value) {fAODfilterBit = value;}; // set AOD filter bit
        
    enum EData_t {
        kData_logpt,
        kData_eta,
	kData_y,
        kData_phi,
	kData_zv,
	kData_centrality,
        kNData
    };
    enum EStage_t {
        kStage_generated,
        kStage_reconstructed,
        kStage_matched,
        kNStages
    };
    
private:
    
    AliAnalysisTaskEffMatrix(const AliAnalysisTaskEffMatrix &); // not implemented
    AliAnalysisTaskEffMatrix &operator=(const AliAnalysisTaskEffMatrix &); // not implemented
    
    /*** private functions ***/
    const AliVVertex *GetPrimaryVertex(AliVEvent *event) const; // get primary vertex
    Bool_t AcceptEvent(AliVEvent *event) const; // accept event
    Bool_t AcceptTrack(AliVTrack *track) const; // accept track
    Bool_t AcceptParticle(AliVParticle *particle) const; // accept particle
    Bool_t AcceptResonance(AliVParticle *particle, Int_t dpdg1, Int_t dpdg2) const; // accept resonance
    Bool_t HasTPCPID(AliVTrack *track) const; // has TPC PID
    Bool_t HasTOFPID(AliVTrack *track) const; // has TOF PID
    
    /*** data members ***/
    Int_t fAODfilterBit; // AOD filter bit
    AliESDtrackCuts *fESDtrackCuts; // ESD track cuts
    Float_t fEtaMin; // eta min
    Float_t fEtaMax; // eta max
    Float_t fPtMin; // pt min
    Float_t fPtMax; // pt max
    Float_t fRapidityMin; // rapidity min
    Float_t fRapidityMax; // rapidity max
    TClonesArray *fMCArray; //! MC array
    TList *fOutputList; //! output list
    TList *fHistoList[kNStages]; //! histo list
    
    ClassDef(AliAnalysisTaskEffMatrix, 1);
};

#endif /* ALIANALYSISTASKEFFMATRIX_H */
