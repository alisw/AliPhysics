//-------------------------------------------------------------------------
//    Description: 
//    This class is abstact (interface) base class for different classes 
//    used for track-by-track data analysis in AliAnalysisTaskLRC
//    Origin: Andrey Ivanov (SPbSU-CERN) anivanov@cern.ch, 
//    Igor Altsybeev (SPbSU-CERN) 
//-------------------------------------------------------------------------
/*  See cxx source for full Copyright notice */

#ifndef ALILRCBASE_H
#define ALILRCBASE_H

#include <TObject.h>



// Include
class TList;
class AliLRCBase: public TObject
{
public:	
	/*enum LRCparticleType
	{
		kLRCany = -1,
		kLRCother = 0,
		kLRCelectron,
		kLRCmuon,
		kLRCpion,
		kLRCkaon,
		kLRCproton,		
	};
		
	enum LRCpidFillCondition
	{
		kLRCpidIgnored,
		kLRCpidForBackwardOnly,
		kLRCpidForForwardOnly,
		kLRCpidForBoth,
	};*/


    // Constructors
    AliLRCBase(): TObject() {}; // Constructor

    // Destructor
    virtual ~AliLRCBase(){};

    virtual Bool_t InitDataMembers() = 0; //Is to be called in CreateOutputObjects method

    // Setters
    virtual  void SetOutputSlotNumber(Int_t SlotNumber) = 0; // Sets number of output slot for LRCProcessor
    //	virtual  void SetPrintInfo( Bool_t _flagPrintInfo )  = 0; // Print info flag

    // Getters
    virtual void GetETAWindows(Double_t &_StartForwardETA,Double_t &_EndForwardETA,Double_t &_StartBakwardETA,Double_t &_EndBakwardETA) = 0;
    virtual void GetPhiWindows(Double_t &_StartForwardPhi,Double_t &_EndForwardPhi,Double_t &_StartBakwardPhi,Double_t &_EndBakwardPhi) = 0;
    virtual TList* CreateOutput() const 	= 0 ; // Creates output object
    virtual TString GetShortDef() const 	= 0  ; // Returns fShortDef value
    virtual Int_t GetOutputSlotNumber() const = 0 ; // Returns number of output slot for LRCProcessor

    // Track by track event import
    virtual void StartEvent() 	= 0;  // Open new Event for track by track event import
    virtual void FinishEvent(Bool_t kDontCount = kFALSE) 	= 0; // Close opened event

    virtual void SetETAWindows( Double_t _StartForwardETA, Double_t _EndForwardETA, Double_t _StartBackwardETA, Double_t _EndBackwardETA ) = 0; // Sets both forward and backward

    virtual void AddTrackPtEta( Double_t Pt, Double_t Eta, Double_t Phi, Short_t Charge = 100, Int_t particleType = -1  )	= 0; // Imports track from the event
    virtual void AddTrackForward( Double_t Pt, Double_t Eta, Double_t Phi, Short_t Charge, Int_t particleType  ) 	= 0; // Imports track to the event directly to Forward window
    virtual void AddTrackBackward( Double_t Pt, Double_t Eta, Double_t Phi, Short_t Charge, Int_t particleType   ) 	= 0; // Imports track to the event directly to Backward window
    virtual void AddTrackPtEtaMixing( Int_t winFB, Double_t Pt, Double_t Eta ,Double_t Phi, Short_t Charge, Int_t particleType = -1 ) = 0;

    virtual void SetEventCentrality(Double_t centrality) = 0;
    virtual void SetEventPlane(Double_t eventPlane) = 0;

    virtual void SetParticleType( char* strForwardOrBackward, char* strPid ) = 0;

    virtual void SetForwardWindowPhi(Double_t StartForwardPhi,Double_t EndForwardPhi) = 0;
    virtual void SetBackwardWindowPhi(Double_t StartBackwardPhi,Double_t EndBackwardPhi) = 0;
 
    virtual  void SetHistPtRange(Double_t LoPt,Double_t HiPt,Int_t PtBins) = 0;  // Sets range for Pt histos axis

    virtual void   Terminate() = 0; //was done for yield study

	//virtual void SetPidFillCondition( LRCparticleType particleType, LRCpidFillCondition fillCond ) 
	//	{ fPidFillCondition = fillCond; 
	//		fWhichParticleToProcess = particleType; }  // how to fill bckwd and fwd hists by pid

    //void SetParticleTypeToFill( LRCparticleType particleType ) { fWhichParticleToProcess = particleType; }  // which particle type in this processor is considered

//protected:
	//LRCparticleType fWhichParticleToProcess; // ! LRC processor sense only this type of particles (-1 by default)
	//LRCpidFillCondition fPidFillCondition;
private:

    virtual  void SetShortDef()	=0;  // Sets fShortDef according to window paramiters

    AliLRCBase(const AliLRCBase&); // not implemented
    AliLRCBase& operator=(const AliLRCBase&); // not implemented



    ClassDef(AliLRCBase, 1);
};

#endif

