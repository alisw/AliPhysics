// This class is creatig TH2D histogramms for Nch - Nch , Nch - Pt , Pt - Pt
// dirtributions for given ETA windows and some supplementary data for Long Range Correlation (LRC) analysis .  
// Class is designid to work with AliAnalysisTaskLRC

// Author : Andrey Ivanov , St.Peterburg State University
// Email: Andrey.Ivanov@cern.ch

/*  See cxx source for full Copyright notice */


#ifndef ALILRCPROCESS_H
#define ALILRCPROCESS_H



// Include
#include "TString.h"
#include "AliLRCBase.h"
#include "THnSparse.h"

class TH1F;
class TH1D;
class TH2D;
class TProfile;
class TList;
//class TClonesArray;//TArrayF;
//class THnSparseD;


enum sparse_pid_list
{
    kSparsePIDany = 0,
    kSparsePIDdefined,
    kSparsePIDpion,
    kSparsePIDkaon,
    kSparsePIDproton,
    kSparsePIDtotal,
};

enum sparse_data_structure
{
    kSparseNf = 0,
    kSparseNb,
    kSparsePtF,
    kSparsePtB,
    /*en_sparse_Et_f,
    en_sparse_Et_b,*/
    en_sparse_N2_f,
    en_sparse_Nf_Nb,

    en_sparse_Ptb_Nf,
    en_sparse_Pt2_f,
    en_sparse_Ptf_Ptb,

    /*en_sparse_Nf_plus,
    en_sparse_Nf_minus,
    en_sparse_Nb_plus,
    en_sparse_Nb_minus,
    en_sparse_Nf_plus_Nb_minus,
    en_sparse_Nb_plus_Nf_minus,*/
    //en_sparse_N2_b,
    kSparseTotal,
};



enum array_labels
{
    en_arr_labels_NN_Nevents = 0,
    en_arr_labels_NN_Nf,
    en_arr_labels_NN_Nb,
    en_arr_labels_NN_N2_f,
    en_arr_labels_NN_Nf_Nb,

    en_arr_labels_PtN_Nevents,
    en_arr_labels_PtN_Nf,
    en_arr_labels_PtN_PtB,
    en_arr_labels_PtN_N2_f,
    en_arr_labels_PtN_Ptb_Nf,

    en_arr_labels_PtPt_Nevents,
    en_arr_labels_PtPt_PtF,
    en_arr_labels_PtPt_PtB,
    en_arr_labels_PtPt_Pt2_f,
    en_arr_labels_PtPt_Ptf_Ptb,

    en_arr_labels_total,
};



class AliLRCProcess: public AliLRCBase
{
public:
    // Constructors
    AliLRCProcess();

    AliLRCProcess(Double_t _StartForwardETA,Double_t _EndForwardETA,Double_t _StartBakwardETA,Double_t _EndBakwardETA ); // Constructor with window setup makes ready-to-run processor

    Bool_t InitDataMembers(); //Is to be called in CreateOutputObjects method

    // Destructor
    virtual ~AliLRCProcess();

    // Setters
    virtual void SetForwardWindow(Double_t StartETA,Double_t EndETA);  //Sets Forward ETA window
    virtual void SetBackwardWindow(Double_t StartETA,Double_t EndETA);  // Sets Backward window
    virtual void SetETAWindows(Double_t _StartForwardETA,Double_t _EndForwardETA,Double_t _StartBakwardETA,Double_t _EndBakwardETA); // Sets both forward and backward windows
    virtual void SetForwardWindowPhi(Double_t StartForwardPhi,Double_t EndForwardPhi){fStartForwardPhi=StartForwardPhi;fEndForwardPhi=EndForwardPhi;}
    virtual void SetBackwardWindowPhi(Double_t StartBackwardPhi,Double_t EndBackwardPhi){fStartBackwardPhi=StartBackwardPhi;fEndBackwardPhi=EndBackwardPhi;}

    virtual void SetDoubleSidedBackwardWindowPhi( Bool_t isDoubleSided ) { fDoubleSidedBackwardPhiWindow = isDoubleSided; }// Number of phi division, for rotating sectors

    virtual void SetEventCentrality(Double_t centrality){ fEventCentrality = centrality;}
    virtual void SetEventPlane(Double_t ) {}


    //virtual void SetNumberOfPhiSectors(Int_t nSectors){ fNumberOfSectors = nSectors; }//fNeedToRotateSector = kTRUE; }


    virtual  void SetHistPtRange(Double_t LoPt,Double_t HiPt,Int_t PtBins);  // Sets range for Pt histos axis
    virtual  void SetHistPtRangeForwardWindowRebinFactor( Int_t ptHistXaxisRebinFactor );
    virtual  void SetHistMultRange(Int_t whichWindow, Int_t LoMult,Int_t HiMult,Int_t MultBins=0); // Sets range for Nch histos axis
    virtual  void SetHistMultRangeHor(Int_t LoMult,Int_t HiMult,Int_t MultBins=0); // Sets range for Nch histos axis
    virtual  void SetHistMultRangeVert(Int_t LoMult,Int_t HiMult,Int_t MultBins=0); // Sets range for Nch histos axis
    virtual  void SetOutputSlotNumber(Int_t SlotNumber); // Sets number of output slot for LRCProcessor

    // Getters
    virtual void GetETAWindows(Double_t &_StartForwardETA,Double_t &_EndForwardETA,Double_t &_StartBakwardETA,Double_t &_EndBakwardETA);
    virtual void GetPhiWindows(Double_t &_StartForwardPhi,Double_t &_EndForwardPhi,Double_t &_StartBakwardPhi,Double_t &_EndBakwardPhi);
    virtual TList* CreateOutput() const ; // Creates output object
    virtual TString GetShortDef() const ; // Returns fShortDef value
    virtual Int_t GetOutputSlotNumber() const; // Returns number of output slot for LRCProcessor

    // Track by track event import
    virtual void StartEvent();  // Open new Event for track by track event import
    virtual void AddTrackPtEta(Double_t Pt, Double_t Eta , Double_t Phi, Short_t Charge = 100, Int_t particleType = -1 ); // Imports track to the event
    virtual void AddTrackForward(Double_t Pt, Double_t Eta , Double_t Phi, Short_t Charge, Int_t particleType  );  // Imports track to the event directly to Forward window
    virtual void AddTrackBackward(Double_t Pt, Double_t Eta , Double_t Phi, Short_t Charge, Int_t particleType  );  // Imports track to the event directly to Backward window
    virtual void AddTrackPtEtaMixing( Int_t winFB, Double_t Pt, Double_t Eta ,Double_t Phi, Short_t Charge, Int_t particleType = -1 );
    virtual void FinishEvent(Bool_t kDontCount = kFALSE); // Close opened event

    //virtual void SetPidToFill(  ) { fWhichParticleToProcess = particleType; }  // which particle type in this processor is considered
    //virtual void SetPidFillCondition( LRCparticleType particleType, LRCpidFillCondition fillCond )
    //	{ fPidFillCondition = fillCond;
    //		fWhichParticleToProcess = particleType; }  // how to fill bckwd and fwd hists by pid
    //virtual LRCparticleType GetPidToFill() { return fWhichParticleToProcess;  }
    //virtual void SetPidForward( LRCparticleType particleType );// { fPidForward = particleType; }  // which pid for fwd
    //virtual void SetPidBackward( LRCparticleType particleType );// { fPidBackward = particleType; }  // which pid for fwd
    virtual void SetParticleType( char* strForwardOrBackward, char* strPid );
    virtual void   Terminate() {} //was done for yield study

    virtual void SetUseSparse( Bool_t flag )            { fUseSparse             = flag;  }
    virtual void SetUseAccumulatingHist( Bool_t flag )  { fUseAccumulatingHist   = flag;  }
    inline Bool_t IsPhiInRange( Double_t phi, Double_t phiBoundMin, Double_t phiBoundMax );

private:

    virtual  void SetShortDef();  // Sets fShortDef according to window paramiters

    //Data Init and checks
    //Bool_t StateCheck();   // Check if data is OK
    Bool_t fIsEventOpend;  // Indicates opened event
    Bool_t fIsOnline;  // Are data objects created
    Bool_t fDisplayInitOnDemandWarning; // Switch warning when InitDataInDemand is called;
    //Bool_t InitDataOnDemand(); // Create data objects


    Bool_t fUseSparse; //flag to use THnSparse
    Bool_t fUseAccumulatingHist; //flag to accumulate some values

    // Statistics
    Int_t fEventCount; //Number of processed events

    // Windows paramiters -----------------------------------

    Double_t fStartForwardETA;  // Forward windos lover rapidity
    Double_t fEndForwardETA;    // Forward window higer rapidity
    Double_t fStartForwardPhi;  // Forward phi angle interval
    Double_t fEndForwardPhi;   // Forward phi angle interval


    Double_t fStartBackwardETA;  // Backward window lover rapidity
    Double_t fEndBackwardETA;    // Backward window higer rapidity

    Double_t fStartBackwardPhi;  // Backward window phi angle interval
    Double_t fEndBackwardPhi;    // Backward window phi angle interval

    Bool_t fDoubleSidedBackwardPhiWindow; // Flag, if true - we look at phi sector and it's opposite one


    Double_t fHiPt;		// Max Pt for histos
    Double_t fLoPt; 		// Min Pt for histos
    Int_t fHiMultHor;		// Min multiplicity for histos (horizontal)
    Int_t fLowMultHor;		// Max multiplicity for histos (horizontal)

    Int_t fHiMultVert;		// Min multiplicity for histos (vertical)
    Int_t fLowMultVert;	// Max multiplicity for histos (vertical)


    Int_t fMultBinsHor;		// N bins for multiplicity (horizontal)
    Int_t fMultBinsVert;	// N bins for multiplicity (vertical)
    Int_t fPtBins;		// N bins for Pt
    Int_t fPtHistXaxisRebinFactor;		// rebinning of Xaxis on PtN and PtPt histos


    // Track by track import values
    Double_t fSumPtFw; 		// Summary Pt forward
    Double_t fSumPtBw; 		// Summary Pt backward
    Double_t fSumPtBw2;		// Summary Pt^2 backward
    Int_t fNchFw; 			// Number of tracks Forward
    //Int_t fNchFwPtPt;		// Number of tracks Forward for PtPt accept conditions
    Int_t fNchBw; 			// Number of tracks Backward

    Int_t fNchFwPlus; 			// Number of plus charged tracks Forward
    Int_t fNchBwPlus; 			// Number of plus charged tracks Backward
    Int_t fNchFwMinus; 			// Number of minus charged tracks Forward
    Int_t fNchBwMinus; 			// Number of minus charged tracks Backward


    //PID data arrays
    Int_t fCorrespondanceWithAliROOTpid[kSparsePIDtotal];
    
    Double_t fSumPtFwPID[kSparsePIDtotal]; 		// Summary Pt forward
    Double_t fSumPtBwPID[kSparsePIDtotal]; 		// Summary Pt backward
    Double_t fSumEtFwPID[kSparsePIDtotal]; 		// Summary Pt forward
    Double_t fSumEtBwPID[kSparsePIDtotal]; 		// Summary Pt backward
    //Double_t fSumPtBw2;		// Summary Pt^2 backward
    Int_t fNchFwPID[kSparsePIDtotal]; 			// Number of tracks Forward
    Int_t fNchBwPID[kSparsePIDtotal]; 			// Number of tracks Backward

    Int_t fNchFwPlusPID[kSparsePIDtotal]; 			// Number of plus charged tracks Forward
    Int_t fNchBwPlusPID[kSparsePIDtotal]; 			// Number of plus charged tracks Backward
    Int_t fNchFwMinusPID[kSparsePIDtotal]; 			// Number of minus charged tracks Forward
    Int_t fNchBwMinusPID[kSparsePIDtotal]; 			// Number of minus charged tracks Backward
    
    
    //Output data
    TList* fOutList;    //! List of output data

    TString fShortDef; // Short desctiption
    Int_t fOutputSlot; // Output slot number for this Processor

    // Total spectras (debugging for TAG selection)
    TH1F        *fHistPt; //! Overal Pt spectrum
    TH1F        *fHistEta; //! Overal Eta spectrum


    // Output histogramms -----------------------------------
    THnSparseD *fHistClouds; //! all LRC Clouds

    TH2D* fHistNN;        //! N-N 2D Profile
    TH2D* fHistPtN;	//! Pt-N 2D Profile
    TH2D* fHistPtPt;	//! Pt-Pt 2D Profile
    TProfile* fProfNberr;	//! Nbackward error Profile PtN
    TProfile* fProfNberrPtPt;	//! Nbackward error Profile PtPt
    TProfile* fProfdPtB;  //! Used to store (in first bin) summary of PtB and its std diviation
    TProfile* fProfTestLRC; //! Diognostic LRC Pt - N correlation

    // Supp. info for windows
    TH1D *fHistSparseDimensionLabeling; //correspondance of THnSparse enums to what we are filling in each dimension
    TH1D *fHistSparsePIDblocksLabeling;
    //Forward
    TH1D* fHistPtForward;   //! Pt spectrum in Forward windows
    TH1D* fHistEtaForward;  //! Eta spectrum in Forward windows
    TH1D* fHistNchForward;  //! Nch spectrum in Forward windows
    TH1D* fHistNchForwardPtPt;  //! Nch spectrum in Forward windows for PtPt accept conditions
    TH1D* fHistPhiForward;  //! Phi spectrum in Forward windows
    TH1D *fHistTracksChargeForward;   //! Charge of accepted tracks in Forward windows

    //Backward
    TH1D* fHistPtBackward;   //! Pt spectrum in Backward windows
    TH1D* fHistEtaBackward;  //! Eta spectrum in Backward windows
    TH1D* fHistNchBackward;  //! Nch spectrum in Backward windows
    TH1D* fHistPhiBackward;  //! Phi spectrum in Backward windows
    TH1D *fHistTracksChargeBackward;   //! Charge of accepted tracks in Backward windows

    /*TClonesArray */
    TH1D *fArrAccumulatedValues;   //! accumulated values for observables
    //TH1D *fHistArrayLabeling; //correspondance of THnSparse enums to what we are filling in each dimension

    
    AliLRCProcess(const AliLRCProcess&); // not implemented
    AliLRCProcess& operator=(const AliLRCProcess&); // not implemented

    //Int_t fWhichParticleToFill; //! LRC processor sense only this type of particles
    //Int_t 	fNumberOfSectors; //! Number of phi division, for rotating sectors
    //Bool_t fNeedToRotateSector;	//! Flag for rotating sectors

    Double_t fEventCentrality;		// for filling hist nF-centrality (for PbPb)
    TH2D* fHistNfCentrality;	//! Nf-Centrality plot
    
    TH1D* fHistTestPIDForward;  //! PID in Forward windows
    TH1D* fHistTestPIDBackward;  //! PID in Backward windows

    TH2D* fHistDifferenceNf;	//! nB-nF as func of nF
    TH2D* fHistDifferenceNb;	//! nF-nB as func of nB

    int fPidForward; //! LRC processor sense only this type of particles (-1 by default)
    int fPidBackward; //! LRC processor sense only this type of particles (-1 by default)


    //LRCparticleType fWhichParticleToProcess; // ! LRC processor sense only this type of particles (-1 by default)
    //LRCpidFillCondition fPidFillCondition;

    ClassDef(AliLRCProcess, 3);
};

#endif

