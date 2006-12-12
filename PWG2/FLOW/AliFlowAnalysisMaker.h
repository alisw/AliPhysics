//////////////////////////////////////////////////////////////////////
//
// $Id$
//
// Author: Emanuele Simili
//
//////////////////////////////////////////////////////////////////////
//
// Description: routines to study Flow from AliFlowEvent, adapted from STAR
// Original Authors:                 Raimond Snellings & Art Poskanzer
//
//////////////////////////////////////////////////////////////////////

#ifndef AliFlowAnalysisMaker_H
#define AliFlowAnalysisMaker_H

#include <stdlib.h>
#include <iostream>

#include "AliFlowTrack.h"
#include "AliFlowV0.h"
#include "AliFlowEvent.h"
#include "AliFlowSelection.h"
#include "AliFlowConstants.h"

#include "TTree.h"
#include "TFile.h"
#include "TString.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TOrdCollection.h"
#include "TMath.h"
#include "TText.h"
#include "TVector2.h"
#include "TString.h"
#include "TList.h" 
#include "TObjArray.h"
#include "TObjString.h"

class AliFlowTrack;
class AliFlowV0;
class AliFlowEvent;
class AliFlowSelection;
class Flow;
class TH1F;
class TH1D;
class TH2F;
class TH2D;
class TH3F;
class TProfile;
class TProfile2D;

class AliFlowAnalysisMaker {

  /*!
    \class AliFlowAnalysisMaker
    Makes histograms for the flow analysis. It reads event and particle quantities
    from AliFlowEvent. It removes autocorrelations of each particle with respect
    to the event plane. For each harmonic and each selection makes a 2D histogram
    of the anisotropic flow, v, vs. y and p_t.
  */
  
public:
  
  AliFlowAnalysisMaker(const Char_t* name="FlowAnalysis");			   // Default constructor
  AliFlowAnalysisMaker(const Char_t* name, const AliFlowSelection& pFlowSelect);   // Constructor with selection object
  AliFlowAnalysisMaker(const AliFlowAnalysisMaker &from) {};
  virtual  ~AliFlowAnalysisMaker();

 // Steps of the flow analysis
  void     InitDefault() ;					// Default sets
  Int_t    Init() ;						// Books histograms for flow analysis
  Int_t    Make() ;						// Open FlowEvents, Gets quantities, Fills histograms
  Int_t    Finish() ;						// Saves histograms, Closes stuff

 // Analysis of 1 event (can be called from outside)
  Bool_t   Analyze(AliFlowEvent* pFlowEvent = 0) ; 		// Fills the defaults histograms (init them first!) and performs the calculation for the given event         

 // Analysis options
  void	   SetEtaSub() ;					// Set eta subevents
  void     SetV1Ep1Ep2(Bool_t v1Ep1Ep2 = kTRUE);		// Switches the v_1{EP1,EP2} calculation on/off
  void	   SetShuffle(Bool_t sh = kTRUE) ;			// Set to re-shuffle evt tracks
  void     SetUsePhiWgt(Bool_t pw = kTRUE) ;			// Set to use phi weights (true by default...if there)
  void 	   SetUseBayWgt(Bool_t bw = kTRUE) ;			// Set to use bayesian weights for p.id. (false by default)

  void     SetUsePtWgt(Bool_t ptw = kTRUE);		        // uses pT as a weight for RP determination
  void     SetUseEtaWgt(Bool_t etw = kTRUE);		        // uses eta as a weight for RP determination
  void     SetUseOnePhiWgt();			        	// just one wgt histogram
  void     SetUseFirstLastPhiWgt();		        	// uses 3 wgt histograms

 // Other options
  void	   SetMakeAll(Bool_t mka=kTRUE) ;			// Set to calculate all events variables in one shoot
  void	   SetFlowForV0(Bool_t fV0 = kTRUE) ;			// Enables Flow study for v0
  void	   SetRedoWgt(Bool_t rd = kTRUE) ;			// Set recalculation of weights

 // Histograms
  void     SetHistoRanges(float etaMin = -1.,float etaMax = 1.,int EtaBins = 40) ; // Sets the histograms' limits
  void     SetPtRange_for_vEta(Float_t lo, Float_t hi);		// Sets the pt range for the v(eta) histograms.
  void     SetEtaRange_for_vPt(Float_t lo, Float_t hi);		// Sets the |eta| range for the v(pt) histograms.

 // Input/Output 
  void	   SetHistFileName(TString name) ;			// Sets output file name
  void	   SetPhiWgtFileName(TString name) ;			// Sets Wgt file name
  void	   SetInputFileName(TString name) ;			// Sets Input File name (the FlowEvents file) 
  void	   SetInputFileNames(TString name) ;			// Sets Input Files (more than 1)
  void	   SetOneInputFile(Bool_t in=kTRUE) ;			// Sets just one Input file 
  TString  GetHistFileName() const ;				// Output File Name
  TString  GetPhiWgtFileName() const ;				// Wgt File Name
  TString  GetInputFileName() const ;				// Input (Flow Events) file Name

 // MC simulation & debug
  void     SetDebugg(Int_t db = 1) ; 				// set the cout's for debug (default is 1)
  void     SetFillLabels(Bool_t lb=kTRUE) ;			// fills labels histogram (execution is heavier)
  void 	   SetMaxLabel(Int_t bin) ; 				// set the last bin in the MC-label histogram

 // Some Results
  Float_t  Res(Int_t eventN, Int_t harN) const ;		// Returns the calculated resolution for the RP
  Float_t  ResErr(Int_t eventN, Int_t harN) const ;		// Returns the estimated error on the resolution 
  void     GetRunBayesian(Double_t bayes[Flow::nPid],int sel=0) ;  // Normalized Particle abundance (all events up to here)
  void     PrintRunBayesian(int sel=0) ; 			// Prints the normalized Particle abundance (up to here)


private:

 // Weightening
  void     WgtChk() ; 						// Check for Read/Write weights, sets labels accordingly
  void     FillEvtPhiWgt() ; 					// Plugs the PhiWeights into the AliFlowEvent
  void     FillWgtArrays(TFile* wgtFile) ;			// Loads phi & bayesian weights from file (flowPhiWgt.hist.root) and fills the specific arrays
  void     FillBayesianWgt(int sel=0) ; 			// Plugs the Bayesian Weights into the AliFlowEvent (selection 0)
  void     Weightening() ;  					// Calculates weights and fills PhiWgt histograms

 // Open a new Imput File / Load a new Event
  AliFlowEvent* GetEvt(Int_t evt = -1) ;			// Gets flow event & V0s
  Bool_t   Open(const Char_t* filename="") ;			// Open FlowEvents file
  Bool_t   FillFromFlowEvent() ;				// Fills internal variables and array from Flow Events

 // Internal methods to fill the histogram
  void     FillEventHistograms() ;				// Fills Events' histograms (from AliFlowEvent)
  void     FillParticleHistograms() ;				// Fills Tracks' histograms (from AliFlowTrack)
  void     FillV0Histograms() ;                                 // Fills V0s' histograms
  int      HarmonicsLoop(float eta,float phi,float pt,int numPid=-1) ; 	// Harmonics & Selections relates histograms (from AliFlowTracks)
  void     FillLabels() ;					// fills an histogram of Labels (the ones from ESD) 

 // Resolution Calculation
  void     Resolution() ;  					// Calculates resolution and mean flow values
  Double_t chi(double res) ;  					// Calculates chi from the event plane resolution
  Double_t resEventPlane(double chi) ;	  			// Calculates the event plane resolution as a function of chi
  Double_t resEventPlaneK2(double chi) ;	  		// Calculates the event plane resolution as a function of chi for the case k=2.
  Double_t resEventPlaneK3(double chi) ;	  		// Calculates the event plane resolution as a function of chi for the case k=3.
  Double_t resEventPlaneK4(double chi) ;	  		// Calculates the event plane resolution as a function of chi for the case k=4.

 // Flags
  Bool_t   mV0 ;						//! correlation analysis is done also for neutral secundary vertex
  Bool_t   mShuffle ;						//! to randomly reshuffle tracks
  Bool_t   mV1Ep1Ep2;						//! Flag for v_1{EP1,EP2} calculation on/off
  Bool_t   mEtaSub;                                             //! eta subevents
  Bool_t   mReadPhiWgt ;			       		//! Flag for reading Wgt histograms from file (automatic if flowPhi.hist.root is there)
  Bool_t   mWritePhiWgt ;			   		//! Flag for writing Wgt histograms to file (automatic if flowPhi.hist.root is NOT there)
  Bool_t   mRedoWgt ;				   		//! Flag for recalculating Wgt histograms (even if flowPhi.hist.root exist)
  Bool_t   mPhiWgt ;						//! Phi Weights are applied to Phi distrib. (default is true)
  Bool_t   mBayWgt ;						//! Bayesian Weights are applied to P.Id. (default is false) 
  Bool_t   mRePid ;						//! Re-Calculates the P.Id. basing on the bayesian wgts (if plugged in)
  Bool_t   mOnlyConstrainable ;					//! to loop just over constrainable tracks
  Bool_t   mMakeAll ;						//! claculates all events vars in one shoot (should run faster)
  Bool_t   mLabellingFirst ;					//! creates label file (internal)
  Bool_t   mLabelling ;						//! takes note of Labels from ESD

  Bool_t   mPtWgt ;  	    					//! flag to use pT as a weight for RP determination
  Bool_t   mEtaWgt ; 	    					//! flag to use eta as a weight for RP determination
  Bool_t   mOnePhiWgt ;  		    			//! if kTRUE: just one phi-wgt histogram, if kFALSE: three phi-wgt histogram (TPC+,TPC-,cross)

 // ...for debugging
  Int_t    mDebug ;          					//! Debug level (0,1,2,3)
  Bool_t   Debug0 ;				        	//! Flag for Debug couts
  Bool_t   Debug1 ;				        	//! Flag for Debug couts lev.1
  Bool_t   Debug2 ;				        	//! Flag for Debug couts lev.2
  Bool_t   Debug3 ;				        	//! Flag for Debug couts lev.3
  Bool_t   Debug4 ;				        	//! Flag for Debug couts lev.3

 // File names
  TString	    mHistFileName ;				//! Output File Name (histograms from flow analysis)
  TString	    mPhiWgtFileName ;				//! Wgt File Name (histograms for weight)
  TString	    mFlowEvtFileName ;				//! Input file name (Flow Events)
  Bool_t 	    fOneInputFile ;				//! one/more imput file
  TObjArray         fFileNames ;

 // Couts
  TString	    mName ;  				        //! Name from the constructor
  TString	    FlowAnalysis ;				//! Default mark for cout's
  TString	    spaces ;					//! Default spaces for cout's

 // variables
  Int_t             nEvents ;					//! number of FlowEvents in file
  Int_t             nTracks ;					//! number of FlowTracks in FlowEvent
  Int_t             nV0s ;					//! number of FlowV0s in FlowEvent
  Int_t             evtN ;					//! current Event number (after selected ++)
  TString           evt_name ; 					//! current Event name (number) 
  Float_t           vertex[3] ;					//! Event's Vertex position 

 // For weights
  Flow::PhiWgt_t nPhiWgt ;					//! PhiWgt Array (all TPC)
  Flow::PhiWgt_t nPhiWgtPlus ;  				//! PhiWgt Array (TPC+)
  Flow::PhiWgt_t nPhiWgtMinus ; 				//! PhiWgt Array (TPC-)
  Flow::PhiWgt_t nPhiWgtCross ;  				//! PhiWgt Array (TPC/)     
 // For bayesian weights
  Double_t nBayWgt[Flow::nSels][Flow::nPid] ; 			//! Bayesian weights (expected particle abundance)

#ifndef __CINT__
  TVector2 mQ[Flow::nSels][Flow::nHars];			//! flow vector
  Float_t  mPsi[Flow::nSels][Flow::nHars];			//! event plane angle
  UInt_t   mMult[Flow::nSels][Flow::nHars];                  	//! multiplicity
  Float_t  m_q[Flow::nSels][Flow::nHars];                    	//! Q/::sqrt(Mult)
  TVector2 mQSub[Flow::nSubs][Flow::nSels][Flow::nHars];      	//! flow vector subs
  Float_t  mPsiSub[Flow::nSubs][Flow::nSels][Flow::nHars];    	//! plane angle of subevents
  UInt_t   mMultSub[Flow::nSubs][Flow::nSels][Flow::nHars];   	//! multiplicity subs
  Float_t  mRes[Flow::nSels][Flow::nHars];			//! event plane resolution
  Float_t  mResErr[Flow::nSels][Flow::nHars];			//! event plane resolution error
#endif /*__CINT__*/

 // Internal pointers
  AliFlowEvent*     pFlowEvent ;      				//! pointer to AliFlowEvent
  AliFlowTrack*     pFlowTrack ;      				//! pointer to AliFlowTrack
  AliFlowV0*        pFlowV0 ;      				//! pointer to AliFlowV0
  AliFlowSelection* pFlowSelect ;     				//! selection object
  TFile*            pFlowEventsFile ; 				//! pointer to FlowEvents file (input)
  TList*            pFlowEventsList ; 				//! FlowEvents in the file (input)
  TObjArray*        pFlowTracks ;     				//! pointer to the TrackCollection
  TObjArray*        pFlowV0s ;     				//! pointer to the V0Collection
  TFile*            histFile ; 				        //! histograms file (output)
  //TTree*            pFlowTree;             			//! flow events TTree (NEW input - not there)

 // for Histograms
  TString           xLabel ;             			//! label axis with rapidity or pseudorapidity
  Float_t 	    mEtaMin ;     				//! histo range (eta)
  Float_t 	    mEtaMax ;     				//! histo range (eta)
  Int_t   	    mNEtaBins ;     				//! histo bins  (eta)
  Float_t 	    mPtRange_for_vEta[2] ;			//! pt range for the v(eta) histograms.
  Float_t 	    mEtaRange_for_vPt[2] ;			//! |eta| range for the v(pt) histograms.
  Int_t 	    maxLabel ;             			//! for the MC labels histogram (max bin)

  TOrdCollection*   mPhiWgtHistList ;     			//! Weights:  histogram list
  TOrdCollection*   mVnResHistList ;     			//! Resolution and Vn:  histogram list
 
// for Single histograms

 // *****************
 // EVENTs HISTOGRAMS
 // *****************
  TH1F*     mHistTrigger;                 		   //!
  TH1F*     mHistMult;                    		   //!
  TH1F*     mHistV0Mult; 				   //!
  TH1F*     mHistOrigMult;                		   //!
  TH1F*     mHistMultOverOrig;            		   //!
  TH1F*     mHistMultEta;                 		   //!
  TH1F*     mHistCent;                    		   //!
  TH1F*     mHistVertexZ;                 		   //!
  TH2F*     mHistVertexXY2D;              		   //!
  TH2F*     mHistEnergyZDC;              		   //!
  TH1F*     mHistPartZDC;              		      	   //!
  TProfile* mHistPidMult;                 		   //!
  TH1F*     mHistBayPidMult;     		      	   //!
  TH1F*     mHistEtaSym;                  		   //!          // ...part also
  TH1F*     mHistEtaSymPart;                  		   //!
  TH2F*     mHistEtaSymVerZ2D;            		   //!          // ...part also
  TH2F*     mHistEtaSymVerZ2DPart;            		   //!
 // selected (TR & V0)
  TH1F*     mHistMultPart;                		   //!
  TH1F*     mHistV0MultPart;         			   //!
  TH1F*     mHistBayPidMultPart;     		      	   //!
  TH1F*     mHistMultPartUnit; 		  	    	   //!
  
 // *****************
 // TRACKs HISTOGRAMS (all tracks)
 // *****************
  TH1F*     mHistPtot ;                 		   //!
  TH1F*     mHistPt ;                 		   	   //!
  TH1F*     mHistCharge;                  		   //!
  TH1F*     mHistDcaGlobal;               		   //!
  TH1F*     mHistDca;                     		   //!
  TH1F*     mHistTransDca;                     		   //!
  TH1F*     mHistChi2;                  		   //!
  TH1F*     mHistLenght;           		           //!
  TH1F*     mHistInvMass ;				   //!
  TH1F*     mHistFitOverMax;           		           //!
  TH2D*     mHistPhiPtCon ;				   //!
  TH2D*     mHistPhiPtUnc ;				   //!
  TH2D*     mHistPtPhiPos ;                 		   //!
  TH2D*     mHistPtPhiNeg ;                 		   //!
  TH3F*     mHistAllEtaPtPhi3D;              		   //!
  TProfile* mHistCosPhi;                  		   //!
  TH2F*     mHistPidPt;              		   	   //!
  TH1F*     mHistPhi ;                 		   	   //!
  TH1F*     mHistPhiCons ;                 		   //!
  TH2D*     mHistYieldAll2D;              		   //!
  TH2D*     mHistYieldCon2D;              		   //!
  TH2D*     mHistYieldUnc2D;              		   //!
  TH3F*     mHistConsEtaPtPhi3D;              		   //!
  TH3F*     mHistGlobEtaPtPhi3D;             		   //! 
  TH3F*     mHistUncEtaPtPhi3D ;             		   //!
  // fit & dE/dX for each detector (all tracks) 
  TH1F*     mHistChi2ITS;                 		   //!
  TH1F*     mHistChi2normITS;                 		   //!
  TH1F*     mHistFitPtsITS;               		   //!
  TH1F*     mHistMaxPtsITS;               		   //!
  TH2F*     mHistMeanDedxPos2DITS;           		   //!
  TH2F*     mHistMeanDedxNeg2DITS;           		   //!
  // -
  TH1F*     mHistChi2TPC;                 		   //!
  TH1F*     mHistChi2normTPC;                 		   //!
  TH1F*     mHistFitPtsTPC;               		   //!
  TH1F*     mHistMaxPtsTPC;               		   //!
  TH1F*     mHistFitOverMaxTPC;           		   //!
  TH2F*     mHistMeanDedxPos2D;           		   //!
  TH2F*     mHistMeanDedxNeg2D;           		   //!
  // -
  TH1F*     mHistChi2TRD;                 		   //!
  TH1F*     mHistChi2normTRD;                 		   //!
  TH1F*     mHistFitPtsTRD;               		   //!
  TH1F*     mHistMaxPtsTRD;               		   //!
  TH2F*     mHistMeanDedxPos2DTRD;           		   //!
  TH2F*     mHistMeanDedxNeg2DTRD;           		   //!
  // -
  TH1F*     mHistChi2TOF;                 		   //!
  TH1F*     mHistChi2normTOF;                 		   //!
  TH1F*     mHistFitPtsTOF;               		   //!
  TH1F*     mHistMaxPtsTOF;               		   //!
  TH2F*     mHistMeanDedxPos2DTOF;           		   //!
  TH2F*     mHistMeanDedxNeg2DTOF;           		   //!
  // detector response for particle type (all tracks, based on Pid)
  TH2F*     mHistMeanTPCPiPlus ;        		   //!
  TH2F*     mHistMeanTPCPiMinus ;       		   //!
  TH2F*     mHistMeanTPCProton ;        		   //!
  TH2F*     mHistMeanTPCPbar ;          		   //!
  TH2F*     mHistMeanTPCKplus ;         		   //!
  TH2F*     mHistMeanTPCKminus ;        		   //!
  TH2F*     mHistMeanTPCDeuteron ;      		   //!
  TH2F*     mHistMeanTPCAntiDeuteron ;  		   //!
  TH2F*     mHistMeanTPCPositron ;      		   //!
  TH2F*     mHistMeanTPCElectron ;      		   //!
  TH2F*     mHistMeanTPCMuonPlus ;      		   //!
  TH2F*     mHistMeanTPCMuonMinus ;     		   //!
  // -
  TH2F*     mHistMeanITSPiPlus ;			   //!
  TH2F*     mHistMeanITSPiMinus ;			   //!
  TH2F*     mHistMeanITSProton ;			   //!
  TH2F*     mHistMeanITSPbar ;  			   //!
  TH2F*     mHistMeanITSKplus ; 			   //!
  TH2F*     mHistMeanITSKminus ;			   //!
  TH2F*     mHistMeanITSDeuteron ;			   //!
  TH2F*     mHistMeanITSAntiDeuteron ;  		   //!
  TH2F*     mHistMeanITSPositron ;			   //!
  TH2F*     mHistMeanITSElectron ;			   //!
  TH2F*     mHistMeanITSMuonPlus ;			   //!
  TH2F*     mHistMeanITSMuonMinus ;			   //!
  // -
  TH2F*     mHistMeanTOFPiPlus ;			   //!
  TH2F*     mHistMeanTOFPiMinus ;			   //!
  TH2F*     mHistMeanTOFProton ;			   //!
  TH2F*     mHistMeanTOFPbar ;  			   //!
  TH2F*     mHistMeanTOFKplus ; 			   //!
  TH2F*     mHistMeanTOFKminus ;			   //!
  TH2F*     mHistMeanTOFDeuteron ;			   //!
  TH2F*     mHistMeanTOFAntiDeuteron ;  		   //!
  TH2F*     mHistMeanTOFPositron ;			   //!
  TH2F*     mHistMeanTOFElectron ;			   //!
  TH2F*     mHistMeanTOFMuonPlus ;			   //!
  TH2F*     mHistMeanTOFMuonMinus ;			   //!
  // -
  TH2F*     mHistMeanTRDPiPlus ;			   //!
  TH2F*     mHistMeanTRDPiMinus ;			   //!
  TH2F*     mHistMeanTRDProton ;			   //!
  TH2F*     mHistMeanTRDPbar ;  			   //!
  TH2F*     mHistMeanTRDKplus ; 			   //!
  TH2F*     mHistMeanTRDKminus ;			   //!
  TH2F*     mHistMeanTRDDeuteron ;			   //!
  TH2F*     mHistMeanTRDAntiDeuteron ;  		   //!
  TH2F*     mHistMeanTRDPositron ;			   //!
  TH2F*     mHistMeanTRDElectron ;			   //!
  TH2F*     mHistMeanTRDMuonPlus ;			   //!
  TH2F*     mHistMeanTRDMuonMinus ;			   //!
  // pid probability for all particle (all tracks)
  TH1F*     mHistPidPiPlus;               		   //!
  TH1F*     mHistPidPiMinus;              		   //!
  TH1F*     mHistPidProton;               		   //!
  TH1F*     mHistPidAntiProton;           		   //!
  TH1F*     mHistPidKplus;                		   //!
  TH1F*     mHistPidKminus;               		   //!
  TH1F*     mHistPidDeuteron;             		   //!
  TH1F*     mHistPidAntiDeuteron;         		   //!
  TH1F*     mHistPidElectron;             		   //!
  TH1F*     mHistPidPositron;             		   //!
  TH1F*     mHistPidMuonMinus;            		   //!
  TH1F*     mHistPidMuonPlus;             		   //!
  // pid probability for particle type (all tracks, based on Pid)
  TH1F*     mHistPidPiPlusPart;           		   //!
  TH1F*     mHistPidPiMinusPart;          		   //!
  TH1F*     mHistPidProtonPart;           		   //!
  TH1F*     mHistPidAntiProtonPart;       		   //!
  TH1F*     mHistPidKplusPart;            		   //!
  TH1F*     mHistPidKminusPart;           		   //!
  TH1F*     mHistPidDeuteronPart;         		   //!
  TH1F*     mHistPidAntiDeuteronPart;     		   //!
  TH1F*     mHistPidElectronPart;         		   //!
  TH1F*     mHistPidPositronPart;         		   //!
  TH1F*     mHistPidMuonMinusPart;        		   //!
  TH1F*     mHistPidMuonPlusPart;         		   //!
  // MC labels from the simulation (all tracks)
  TH2F*     mLabHist;        	          		   //! 
 // *****************
 // selected TRACKS
 // *****************
  TProfile* mHistBinEta;                  		   //!
  TProfile* mHistBinPt;                   		   //!
  //
  TH3F*     mHistEtaPtPhi3DPart ;			   //!
  TH2D*     mHistYieldPart2D;             		   //!
  TH1F*     mHistDcaGlobalPart ;			   //!
  TH1F*     mHistInvMassPart ;			 	   //!
  TH3F*     mHistEtaPtPhi3DOut ;			   //!
  TH2D*     mHistYieldOut2D;             		   //!
  TH1F*     mHistDcaGlobalOut ;				   //!
  TH1F*     mHistInvMassOut ;				   //!
  TH3F*     mHistMeanDedxPos3DPart ;			   //!
  TH3F*     mHistMeanDedxNeg3DPart ;			   //!
  TH3F*     mHistMeanDedxPos3DPartITS ;			   //!
  TH3F*     mHistMeanDedxNeg3DPartITS ;			   //!
//

 // *****************
 // V0s HISTOGRAMS (all v0s)
 // *****************
  TH1F*     mHistV0Mass; 	  			   //!
  TH3F*     mHistV0EtaPtPhi3D;	  			   //!
  TH2D*     mHistV0YieldAll2D;    			   //!
  TH1F*     mHistV0Dca;	  	  			   //!
  TH1F*     mHistV0Chi2;	  			   //!
  TH1F*     mHistV0Lenght;	  			   //!
  TH1F*     mHistV0Sigma;	  			   //!
  TProfile* mHistV0CosPhi;	           		   //! 
  TH2D*     mHistV0MassPtSlices;    			   //!
 // *****************
 // selected V0s
 // *****************
  TProfile* mHistV0BinEta;	           		   //! 
  TProfile* mHistV0BinPt;	           		   //! 
  TProfile* mHistV0sbBinEta;	           		   //! 
  TProfile* mHistV0sbBinPt;	           		   //!
  //
  TH1F*     mHistV0MassWin ;     			   //!
  TH3F*     mHistV0EtaPtPhi3DPart ; 			   //!
  TH2D*     mHistV0YieldPart2D;    			   //!
  TH1F*     mHistV0DcaPart ;          			   //!
  TH1F*     mHistV0LenghtPart ;	  			   //!
  TH1F*     mHistV0sbMassSide ;    			   //!
  TH3F*     mHistV0sbEtaPtPhi3DPart ; 			   //!
  TH2D*     mHistV0sbYieldPart2D;    			   //!
  TH1F*     mHistV0sbDcaPart ;          		   //!
  TH1F*     mHistV0sbLenghtPart ;	  		   //!

// for each harmonic, each selection, and each sub-event

 // *****************
 // SUB-EVENTs HISTOGRAMS
 // *****************
  struct histSubHars {
   TH1F*     mHistPsiSubs;
  };
  struct histSubs;	
  friend struct histSubs;
  struct histSubs {
   struct histSubHars histSubHar[Flow::nHars];
  };
  struct histSubs histSub[Flow::nSels*Flow::nSubs];        //!

// for each harmonic and each selection

  struct histFullHars 
  {
   // weights
    TH1D*       mHistPhiPlus;
    TH1D*       mHistPhiMinus;
    TH1D*       mHistPhiAll;
    TH1D*       mHistPhiWgtPlus;
    TH1D*       mHistPhiWgtMinus;
    TH1D*       mHistPhiWgtAll;
    TH1D*       mHistPhiFlatPlus;
    TH1D*       mHistPhiFlatMinus;
    TH1D*       mHistPhiFlatAll;
    TH1D*       mHistPhi;
    TH1D*       mHistPhiWgt;
    TH1D*       mHistPhiFlat;
   // flow (events)
    TH1F*       mHistPsi;
    TH1F*       mHistPsiSubCorr;
    TH1F*       mHistPsiSubCorrDiff;
    TH1F*       mHistPsi_Diff;
    TH1F*       mHistMult;
    TH1F*       mHist_q;
   // flow (tracks)
    TH1F*       mHistPhiCorr;
    TProfile2D* mHist_vObs2D;
    TProfile*   mHist_vObsEta;
    TProfile*   mHist_vObsPt;
    TH2D*       mHist_v2D;
    TH1D*       mHist_vEta;
    TH1D*       mHist_vPt;
   // flow (v0s)
    TH1F*       mHistV0PhiCorr;   
    TProfile2D* mHistV0_vObs2D;   
    TProfile*   mHistV0_vObsEta;  
    TProfile*   mHistV0_vObsPt;   
    TH2D*	mHistV0_v2D;	  
    TH1D*	mHistV0_vEta;	  
    TH1D*	mHistV0_vPt;	  
   // flow (v0s sidebands)
    TProfile*   mHistV0sb_vObsEta_sx ;  
    TProfile*   mHistV0sb_vObsPt_sx ;   
    TProfile*   mHistV0sb_vObsEta_dx ;  
    TProfile*   mHistV0sb_vObsPt_dx ;   
    TH1F*       mHistV0sbPhiCorr ;   
    TProfile2D* mHistV0sb_vObs2D ;   
    TProfile*   mHistV0sb_vObsEta ;  
    TProfile*   mHistV0sb_vObsPt ;   
    TH2D*	mHistV0sb_v2D ;      
    TH1D*	mHistV0sb_vEta ;     
    TH1D*	mHistV0sb_vPt ;      
   // check (tracks used for R.P.)
    TH1F*       mHistYieldPt ;
    TH3F*       mHistEtaPtPhi3D ;
    TH2D*       mHistYield2D ;
    TH1F*       mHistDcaGlob ;
   // check (tracks excluded)
    TH1F*	mHistYieldPtout;
    TH3F*	mHistEtaPtPhi3Dout ;
    TH2D*	mHistYield2Dout ;
    TH1F*	mHistDcaGlobout ;
  };

// for each selection

  struct histFulls;	
  friend struct histFulls;
  struct histFulls 
  {
   TH1F*     mHistBayPidMult;
  // flow (events)
   TProfile* mHistCos;
   TH1F*     mHistRes;
   TProfile* mHist_vObs;
   TH1D*     mHist_v;
   TProfile* mHistV0_vObs;
   TProfile* mHistV0sb_vObs_sx;
   TProfile* mHistV0sb_vObs_dx;
   TH1D*     mHistV0_v;
  // wgt, evts, trks, v0s (as defined above)
   struct histFullHars histFullHar[Flow::nHars];
  };
  struct histFulls histFull[Flow::nSels];                     //!

  ClassDef(AliFlowAnalysisMaker,2)              // macro for rootcint
};

#endif



// lame = not productive; poorly designed; uncool ...
