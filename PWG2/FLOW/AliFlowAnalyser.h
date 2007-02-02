//////////////////////////////////////////////////////////////////////
//
// $Id$
//
// Author: Emanuele Simili
//
//////////////////////////////////////////////////////////////////////
//
// Description: flow analysis for AliFlowEvent(s), adapted from STAR .
// Original Authors:                 Raimond Snellings & Art Poskanzer
//
//////////////////////////////////////////////////////////////////////

#ifndef ALIFLOWANALYSER_H
#define ALIFLOWANALYSER_H

#include <TVector2.h>
#include <TFile.h>
#include "AliFlowConstants.h"

class TH1F;
class TH1D;
class TH2F;
class TH2D;
class TH3F;
class TProfile;
class TProfile2D;
class TOrdCollection;

class AliFlowTrack;
class AliFlowV0;
class AliFlowEvent;
class AliFlowSelection;
class Flow;

class AliFlowAnalyser {


public:
  
  AliFlowAnalyser(const AliFlowSelection* flowSelect = 0); 	// Constructor with selection object (default selection if no one given)
  virtual  ~AliFlowAnalyser(); 					// Default destructor (no actions)

 // Steps of the flow analysis
  Bool_t   Init() ;						// Books histograms for flow analysis
  Bool_t   Finish() ;						// Saves histograms, Closes stuff

 // Analysis of 1 event (can be called from outside)
  Bool_t   Analyse(AliFlowEvent* flowEvent = 0) ; 		// Fills the defaults histograms (init them first!) and performs the calculation for the given event         

 // Resolution corrections to v_n (call it at the end of the evts loop)
  Bool_t   Resolution() ;				 	// Calculates resolution and mean flow values
 
 // Weights calculation and saving (call it at the end of the evts loop)
  void     Weightening() ;  					// Calculates weights and fills PhiWgt histograms

 // Options
  void	   SetEtaSub(Bool_t es = kTRUE)          		{ fEtaSub = es ; }		// Set eta subevents
  void     SetV1Ep1Ep2(Bool_t v1Ep1Ep2 = kTRUE) 		{ fV1Ep1Ep2 = v1Ep1Ep2 ; }	// Switches the v_1{EP1,EP2} calculation on/off
  void	   SetShuffle(Bool_t sh = kTRUE)  			{ fShuffle = sh ; }		// Set to re-shuffle evt tracks
  void     SetUsePhiWgt(Bool_t pw = kTRUE)  			{ fReadPhiWgt = pw ; }		// Set to use phi weights (true by default...if there)
  void 	   SetUseBayWgt(Bool_t bw = kTRUE) 			{ fBayWgt = bw ; } 		// Set to use bayesian weights for p.id. (false by default)
  void     SetUsePtWgt(Bool_t ptw = kTRUE)		        { fPtWgt = ptw ; }	 	// uses pT as a weight for RP determination
  void     SetUseEtaWgt(Bool_t etw = kTRUE)		        { fEtaWgt = etw ; }	 	// uses eta as a weight for RP determination
  void     SetUseOnePhiWgt(Bool_t opw = kTRUE)			{ fOnePhiWgt = opw ; } 		// just one wgt histogram
  void     SetUseFirstLastPhiWgt(Bool_t flw = kTRUE)		{ fOnePhiWgt = !flw ; }		// uses 3 wgt histograms
  void	   SetFlowForV0(Bool_t v0 = kTRUE) 			{ fV0loop = v0 ; }		// Enables Flow study for v0
  void	   SetTrackLoop(Bool_t trkl = kTRUE) 			{ fTrackLoop = trkl ; }		// Enables Tracks loop (keep it kTRUE)
  //void     SetDebugg(Int_t db = 1) ; 				// set the cout's for debug (default is 1)

 // Histograms
  void     SetPtRangevEta(Float_t lo, Float_t hi)		{ fPtRangevEta[0] = lo ; fPtRangevEta[1] = hi ; }  // Sets the pt range for the v(eta) histograms.
  void     SetEtaRangevPt(Float_t lo, Float_t hi)		{ fEtaRangevPt[0] = lo ; fEtaRangevPt[1] = hi ; }  // Sets the |eta| range for the v(pt) histograms.
  //void     SetMaxLabel(Int_t lab = 100)    			{ fMaxLabel = lab ; }
  
 // Output 
  void	   SetHistFileName(TString name) 			{ fHistFileName = name ; }  	// Sets output file name
  TString  GetHistFileName() const				{ return fHistFileName ; }

 // Phi Weights 
  TString  GetWgtFileName() const 				{ return (TString)fPhiWgtFile->GetName() ; }
  void     FillWgtArrays(TFile* wgtFile) ;			// Loads phi & bayesian weights from file (flowPhiWgt.hist.root) and fills the arrays

 // Results
  Float_t  GetRunBayesian(Int_t nPid=2, Int_t selN=0) ;  	// Normalized Particle abundance (all events up to here)
  void     PrintRunBayesian(Int_t selN=0) ; 			// Prints the normalized Particle abundance (up to here)
  void     PrintEventQuantities() ; 				// Prints event by event calculated quantities
  Float_t  Res(Int_t eventN, Int_t harN) const 	     		{ return fRes[eventN][harN]; }	  // Returns the calculated resolution for the RP
  Float_t  ResErr(Int_t eventN, Int_t harN) const 	        { return fResErr[eventN][harN]; } // Returns the estimated error on the resolution 


 protected:
 
 // Internal methods to fill the histogram
  Bool_t   FillFromFlowEvent(AliFlowEvent* fFlowEvent) ;	// Fills internal variables and array from Flow Events
  void     FillEventHistograms(AliFlowEvent* fFlowEvent) ;	// Fills Events' histograms (from AliFlowEvent)
  void     FillParticleHistograms(TObjArray* fFlowTracks) ;  	// Fills Tracks' histograms (from AliFlowTrack)
  void     FillV0Histograms(TObjArray* fFlowV0s) ; 		// Fills V0s' histograms
  Int_t    HarmonicsLoop(AliFlowTrack* fFlowTrack) ; 		// Harmonics & Selections histograms (from AliFlowTracks)
  //void     FillLabels() ;					// fills an histogram of Labels (the ones from ESD) 

 // Weights plugged to the event
  void     FillBayesianWgt(AliFlowEvent* fFlowEvent) ; 		// Plugs the bayesian weights (fBayesianWgt[0]*) into the AliFlowEvent
  void 	   FillEvtPhiWgt(AliFlowEvent* fFlowEvent) ; 		// Plugs the PhiWeights (fPhiWgt*, etc.) into the AliFlowEvent
 
 // Resolution Calculation
  Double_t Chi(Double_t res) ;  			 	// Calculates chi from the event plane resolution
  Double_t ResEventPlane(Double_t chi) ;		 	// Calculates the event plane resolution as a function of chi
  Double_t ResEventPlaneK2(Double_t chi) ;		 	// Calculates the event plane resolution as a function of chi for the case k=2.
  Double_t ResEventPlaneK3(Double_t chi) ;		 	// Calculates the event plane resolution as a function of chi for the case k=3.
  Double_t ResEventPlaneK4(Double_t chi) ;		 	// Calculates the event plane resolution as a function of chi for the case k=4.


 private:

 // Flags
  Bool_t   	   fTrackLoop ;		     			//! tracks main loop
  Bool_t   	   fV0loop ;		     			//! correlation analysis is done also for neutral secundary vertex
  Bool_t   	   fShuffle ;		     			//! to randomly reshuffle tracks
  Bool_t   	   fV1Ep1Ep2;		     			//! Flag for v_1{EP1,EP2} calculation on/off
  Bool_t   	   fEtaSub;		     			//! eta subevents
  Bool_t   	   fReadPhiWgt ;		     		//! Phi Weights are applied to Phi distrib. (default is false)
  Bool_t   	   fBayWgt ;		     			//! Bayesian Weights are applied to P.Id. (default is false) 
  Bool_t   	   fRePid ;		     			//! Re-Calculates the P.Id. basing on the bayesian wgts (if plugged in)

  Bool_t   	   fPtWgt ;		     			//! flag to use pT as a weight for RP determination
  Bool_t   	   fEtaWgt ;		     			//! flag to use eta as a weight for RP determination
  Bool_t   	   fOnePhiWgt ; 	     			//! if kTRUE: just one phi-wgt histogram, if kFALSE: three phi-wgt histogram (TPC+,TPC-,cross)

 // Files
  TFile*           fHistFile ; 				        //! histograms file (output)
  TFile* 	   fPhiWgtFile ; 			 	//! phi weight file 
  TString	   fHistFileName ;			    	//! Output File Name (histograms from flow analysis)
  //TString	   fFlowEvtFileName ;			    	//! Input file name (Flow Events)

 // enumerators 			    
  Int_t            fEventNumber ;	  		    	//! progressive enumeration of AliFlowEvents
  Int_t            fTrackNumber ;	  		    	//! progressive enumeration of AliFlowTracks
  Int_t            fV0Number ;	  		            	//! progressive enumeration of AliFlowV0s
  //Int_t 	   fNumberOfEvents ;			    	//! total number of AliFlowEvents in file
  Int_t 	   fNumberOfTracks ;			    	//! total number of tracks in the current event
  Int_t 	   fNumberOfV0s ;			    	//! total number of v0s in the current event
  Int_t            fPidId ;	  		    		//! Particle Id hypothesys of the track (0..4 for e,mu,pi,k,p)

 // Internal pointers
  AliFlowEvent*     fFlowEvent ;      				//! pointer to AliFlowEvent
  AliFlowTrack*     fFlowTrack ;      				//! pointer to AliFlowTrack
  AliFlowV0*        fFlowV0 ;      				//! pointer to AliFlowV0
  AliFlowSelection* fFlowSelect ;     				//! selection object
  TObjArray*        fFlowTracks ;     				//! pointer to the TrackCollection
  TObjArray*        fFlowV0s ;     				//! pointer to the V0Collection

  Float_t           fVertex[3] ;				//! Event's Vertex position 

 // For weights
  Int_t      	    fPhiBins ;     				//! n. of phi bins     
  Float_t 	    fPhiMin ;    				//! wgt histo range (phi)
  Float_t  	    fPhiMax ;     				//! wgt histo range (phi) 

  Flow::PhiWgt_t    fPhiWgt ;	  			    	//! PhiWgt Array (all TPC)
  Flow::PhiWgt_t    fPhiWgtPlus ;  			    	//! PhiWgt Array (TPC+)
  Flow::PhiWgt_t    fPhiWgtMinus ;  			    	//! PhiWgt Array (TPC-)
  Flow::PhiWgt_t    fPhiWgtCross ;  			    	//! PhiWgt Array (TPC/)     

 // For bayesian weights
  Double_t 	    fBayesianWgt[Flow::nSels][Flow::nPid] ;  	//! Bayesian weights (expected particle abundance)

#ifndef __CINT__
  TVector2 fQ[Flow::nSels][Flow::nHars];			//! flow vector
  Float_t  fPsi[Flow::nSels][Flow::nHars];			//! event plane angle
  UInt_t   fMult[Flow::nSels][Flow::nHars];                  	//! multiplicity
  Float_t  fQnorm[Flow::nSels][Flow::nHars];                    //! Q/Sqrt(Mult)
  TVector2 fQSub[Flow::nSubs][Flow::nSels][Flow::nHars];      	//! flow vector subs
  Float_t  fPsiSub[Flow::nSubs][Flow::nSels][Flow::nHars];    	//! plane angle of subevents
  UInt_t   fMultSub[Flow::nSubs][Flow::nSels][Flow::nHars];   	//! multiplicity subs
  Float_t  fRes[Flow::nSels][Flow::nHars];			//! event plane resolution
  Float_t  fResErr[Flow::nSels][Flow::nHars];			//! event plane resolution error
#endif /*__CINT__*/

 // for Histograms
  TString           fLabel ;             			//! label axis : rapidity or pseudorapidity
  Float_t 	    fEtaMin ;     				//! histo range (eta)
  Float_t 	    fEtaMax ;     				//! histo range (eta) 
  Float_t 	    fPtMin ; 					//! histo range (pt)	   
  Float_t 	    fPtMax ; 					//! histo range (pt) 	   
  Float_t 	    fPtMaxPart ;   				//! max pt for _part histo
  Int_t   	    fEtaBins ;     				//! n. of eta bins
  Int_t   	    fPtBins ;      				//! n. of pt bins     
  Int_t 	    fPtBinsPart ;   				//! n. of pt bins for _part histo
  Float_t 	    fPtRangevEta[2] ;				//! pt range for the v(eta) histograms.
  Float_t 	    fEtaRangevPt[2] ;				//! |eta| range for the v(pt) histograms.
  Int_t 	    fMaxLabel ;             			//! for the MC labels histogram (max bin)

  TOrdCollection*   fPhiWgtHistList ;     			//! Weights:  histogram list
  TOrdCollection*   fVnResHistList ;     			//! Resolution and Vn:  histogram list
 
// for Single histograms

 // *****************
 // EVENTs HISTOGRAMS
 // *****************
  TH1F*     fHistTrigger;                 		   //!
  TH1F*     fHistMult;                    		   //!
  TH1F*     fHistV0Mult; 				   //!
  TH1F*     fHistOrigMult;                		   //!
  TH1F*     fHistMultOverOrig;            		   //!
  TH1F*     fHistMultEta;                 		   //!
  TH1F*     fHistCent;                    		   //!
  TH1F*     fHistVertexZ;                 		   //!
  TH2F*     fHistVertexXY2D;              		   //!
  TH2F*     fHistEnergyZDC;              		   //!
  TH1F*     fHistPartZDC;              		      	   //!
  TProfile* fHistPidMult;                 		   //!
  TH1F*     fHistBayPidMult;     		      	   //!
  TH1F*     fHistEtaSym;                  		   //!          // ...part also
  TH1F*     fHistEtaSymPart;                  		   //!
  TH2F*     fHistEtaSymVerZ2D;            		   //!          // ...part also
  TH2F*     fHistEtaSymVerZ2DPart;            		   //!
 // selected (TR & V0)
  TH1F*     fHistMultPart;                		   //!
  TH1F*     fHistV0MultPart;         			   //!
  TH1F*     fHistBayPidMultPart;     		      	   //!
  TH1F*     fHistMultPartUnit; 		  	    	   //!
  
 // *****************
 // TRACKs HISTOGRAMS (all tracks)
 // *****************
  TH1F*     fHistPtot ;                 		   //!
  TH1F*     fHistPt ;                 		   	   //!
  TH1F*     fHistCharge;                  		   //!
  TH1F*     fHistDcaGlobal;               		   //!
  TH1F*     fHistDca;                     		   //!
  TH1F*     fHistTransDca;                     		   //!
  TH1F*     fHistChi2;                  		   //!
  TH1F*     fHistLenght;           		           //!
  TH1F*     fHistInvMass ;				   //!
  TH1F*     fHistFitOverMax;           		           //!
  TH2D*     fHistPhiPtCon ;				   //!
  TH2D*     fHistPhiPtUnc ;				   //!
  TH2D*     fHistPtPhiPos ;                 		   //!
  TH2D*     fHistPtPhiNeg ;                 		   //!
  TH3F*     fHistAllEtaPtPhi3D;              		   //!
  TProfile* fHistCosPhi;                  		   //!
  TH2F*     fHistPidPt;              		   	   //!
  TH1F*     fHistPhi ;                 		   	   //!
  TH1F*     fHistPhiCons ;                 		   //!
  TH2D*     fHistYieldAll2D;              		   //!
  TH2D*     fHistYieldCon2D;              		   //!
  TH2D*     fHistYieldUnc2D;              		   //!
  TH3F*     fHistConsEtaPtPhi3D;              		   //!
  TH3F*     fHistGlobEtaPtPhi3D;             		   //! 
  TH3F*     fHistUncEtaPtPhi3D ;             		   //!
  // fit & dE/dX for each detector (all tracks) 
  TH1F*     fHistChi2ITS;                 		   //!
  TH1F*     fHistChi2normITS;                 		   //!
  TH1F*     fHistFitPtsITS;               		   //!
  TH1F*     fHistMaxPtsITS;               		   //!
  TH2F*     fHistMeanDedxPos2DITS;           		   //!
  TH2F*     fHistMeanDedxNeg2DITS;           		   //!
  // -
  TH1F*     fHistChi2TPC;                 		   //!
  TH1F*     fHistChi2normTPC;                 		   //!
  TH1F*     fHistFitPtsTPC;               		   //!
  TH1F*     fHistMaxPtsTPC;               		   //!
  TH1F*     fHistFitOverMaxTPC;           		   //!
  TH2F*     fHistMeanDedxPos2D;           		   //!
  TH2F*     fHistMeanDedxNeg2D;           		   //!
  // -
  TH1F*     fHistChi2TRD;                 		   //!
  TH1F*     fHistChi2normTRD;                 		   //!
  TH1F*     fHistFitPtsTRD;               		   //!
  TH1F*     fHistMaxPtsTRD;               		   //!
  TH2F*     fHistMeanDedxPos2DTRD;           		   //!
  TH2F*     fHistMeanDedxNeg2DTRD;           		   //!
  // -
  TH1F*     fHistChi2TOF;                 		   //!
  TH1F*     fHistChi2normTOF;                 		   //!
  TH1F*     fHistFitPtsTOF;               		   //!
  TH1F*     fHistMaxPtsTOF;               		   //!
  TH2F*     fHistMeanDedxPos2DTOF;           		   //!
  TH2F*     fHistMeanDedxNeg2DTOF;           		   //!
  // detector response for particle type (all tracks, based on Pid)
  TH2F*     fHistMeanTPCPiPlus ;        		   //!
  TH2F*     fHistMeanTPCPiMinus ;       		   //!
  TH2F*     fHistMeanTPCProton ;        		   //!
  TH2F*     fHistMeanTPCPbar ;          		   //!
  TH2F*     fHistMeanTPCKplus ;         		   //!
  TH2F*     fHistMeanTPCKminus ;        		   //!
  TH2F*     fHistMeanTPCDeuteron ;      		   //!
  TH2F*     fHistMeanTPCAntiDeuteron ;  		   //!
  TH2F*     fHistMeanTPCPositron ;      		   //!
  TH2F*     fHistMeanTPCElectron ;      		   //!
  TH2F*     fHistMeanTPCMuonPlus ;      		   //!
  TH2F*     fHistMeanTPCMuonMinus ;     		   //!
  // -
  TH2F*     fHistMeanITSPiPlus ;			   //!
  TH2F*     fHistMeanITSPiMinus ;			   //!
  TH2F*     fHistMeanITSProton ;			   //!
  TH2F*     fHistMeanITSPbar ;  			   //!
  TH2F*     fHistMeanITSKplus ; 			   //!
  TH2F*     fHistMeanITSKminus ;			   //!
  TH2F*     fHistMeanITSDeuteron ;			   //!
  TH2F*     fHistMeanITSAntiDeuteron ;  		   //!
  TH2F*     fHistMeanITSPositron ;			   //!
  TH2F*     fHistMeanITSElectron ;			   //!
  TH2F*     fHistMeanITSMuonPlus ;			   //!
  TH2F*     fHistMeanITSMuonMinus ;			   //!
  // -
  TH2F*     fHistMeanTOFPiPlus ;			   //!
  TH2F*     fHistMeanTOFPiMinus ;			   //!
  TH2F*     fHistMeanTOFProton ;			   //!
  TH2F*     fHistMeanTOFPbar ;  			   //!
  TH2F*     fHistMeanTOFKplus ; 			   //!
  TH2F*     fHistMeanTOFKminus ;			   //!
  TH2F*     fHistMeanTOFDeuteron ;			   //!
  TH2F*     fHistMeanTOFAntiDeuteron ;  		   //!
  TH2F*     fHistMeanTOFPositron ;			   //!
  TH2F*     fHistMeanTOFElectron ;			   //!
  TH2F*     fHistMeanTOFMuonPlus ;			   //!
  TH2F*     fHistMeanTOFMuonMinus ;			   //!
  // -
  TH2F*     fHistMeanTRDPiPlus ;			   //!
  TH2F*     fHistMeanTRDPiMinus ;			   //!
  TH2F*     fHistMeanTRDProton ;			   //!
  TH2F*     fHistMeanTRDPbar ;  			   //!
  TH2F*     fHistMeanTRDKplus ; 			   //!
  TH2F*     fHistMeanTRDKminus ;			   //!
  TH2F*     fHistMeanTRDDeuteron ;			   //!
  TH2F*     fHistMeanTRDAntiDeuteron ;  		   //!
  TH2F*     fHistMeanTRDPositron ;			   //!
  TH2F*     fHistMeanTRDElectron ;			   //!
  TH2F*     fHistMeanTRDMuonPlus ;			   //!
  TH2F*     fHistMeanTRDMuonMinus ;			   //!
  // pid probability for all particle (all tracks)
  TH1F*     fHistPidPiPlus;               		   //!
  TH1F*     fHistPidPiMinus;              		   //!
  TH1F*     fHistPidProton;               		   //!
  TH1F*     fHistPidAntiProton;           		   //!
  TH1F*     fHistPidKplus;                		   //!
  TH1F*     fHistPidKminus;               		   //!
  TH1F*     fHistPidDeuteron;             		   //!
  TH1F*     fHistPidAntiDeuteron;         		   //!
  TH1F*     fHistPidElectron;             		   //!
  TH1F*     fHistPidPositron;             		   //!
  TH1F*     fHistPidMuonMinus;            		   //!
  TH1F*     fHistPidMuonPlus;             		   //!
  // pid probability for particle type (all tracks, based on Pid)
  TH1F*     fHistPidPiPlusPart;           		   //!
  TH1F*     fHistPidPiMinusPart;          		   //!
  TH1F*     fHistPidProtonPart;           		   //!
  TH1F*     fHistPidAntiProtonPart;       		   //!
  TH1F*     fHistPidKplusPart;            		   //!
  TH1F*     fHistPidKminusPart;           		   //!
  TH1F*     fHistPidDeuteronPart;         		   //!
  TH1F*     fHistPidAntiDeuteronPart;     		   //!
  TH1F*     fHistPidElectronPart;         		   //!
  TH1F*     fHistPidPositronPart;         		   //!
  TH1F*     fHistPidMuonMinusPart;        		   //!
  TH1F*     fHistPidMuonPlusPart;         		   //!
  // MC labels from the simulation (all tracks)
  TH2F*     mLabHist;        	          		   //! 
 // *****************
 // selected TRACKS
 // *****************
  TProfile* fHistBinEta;                  		   //!
  TProfile* fHistBinPt;                   		   //!
  //
  TH3F*     fHistEtaPtPhi3DPart ;			   //!
  TH2D*     fHistYieldPart2D;             		   //!
  TH1F*     fHistDcaGlobalPart ;			   //!
  TH1F*     fHistInvMassPart ;			 	   //!
  TH3F*     fHistEtaPtPhi3DOut ;			   //!
  TH2D*     fHistYieldOut2D;             		   //!
  TH1F*     fHistDcaGlobalOut ;				   //!
  TH1F*     fHistInvMassOut ;				   //!
  TH3F*     fHistMeanDedxPos3DPart ;			   //!
  TH3F*     fHistMeanDedxNeg3DPart ;			   //!
  TH3F*     fHistMeanDedxPos3DPartITS ;			   //!
  TH3F*     fHistMeanDedxNeg3DPartITS ;			   //!
//

 // *****************
 // V0s HISTOGRAMS (all v0s)
 // *****************
  TH1F*     fHistV0Mass; 	  			   //!
  TH3F*     fHistV0EtaPtPhi3D;	  			   //!
  TH2D*     fHistV0YieldAll2D;    			   //!
  TH2D*     fHistV0PYall2D;    			   	   //!
  TH1F*     fHistV0Dca;	  	  			   //!
  TH1F*     fHistV0Chi2;	  			   //!
  TH1F*     fHistV0Lenght;	  			   //!
  TH1F*     fHistV0Sigma;	  			   //!
  TProfile* fHistV0CosPhi;	           		   //! 
  TH2D*     fHistV0MassPtSlices;    			   //!
 // *****************
 // selected V0s
 // *****************
  TProfile* fHistV0BinEta;	           		   //! 
  TProfile* fHistV0BinPt;	           		   //! 
  TProfile* fHistV0sbBinEta;	           		   //! 
  TProfile* fHistV0sbBinPt;	           		   //!
  //
  TH1F*     fHistV0MassWin ;     			   //!
  TH3F*     fHistV0EtaPtPhi3DPart ; 			   //!
  TH2D*     fHistV0YieldPart2D;    			   //!
  TH1F*     fHistV0DcaPart ;          			   //!
  TH1F*     fHistV0LenghtPart ;	  			   //!
  TH1F*     fHistV0sbMassSide ;    			   //!
  TH3F*     fHistV0sbEtaPtPhi3DPart ; 			   //!
  TH2D*     fHistV0sbYieldPart2D;    			   //!
  TH1F*     fHistV0sbDcaPart ;          		   //!
  TH1F*     fHistV0sbLenghtPart ;	  		   //!

// for each harmonic, each selection, and each sub-event

 // *****************
 // SUB-EVENTs HISTOGRAMS
 // *****************
  struct histSubHars {
   TH1F*     fHistPsiSubs;
  };
  struct histSubs;	
  friend struct histSubs;
  struct histSubs {
   struct histSubHars fHistSubHar[Flow::nHars];
  };
  struct histSubs fHistSub[Flow::nSels*Flow::nSubs];        //!

// for each harmonic and each selection

  struct histFullHars 
  {
   // weights
    TH1D*       fHistPhiPlus;
    TH1D*       fHistPhiMinus;
    TH1D*       fHistPhiAll;
    TH1D*       fHistPhiWgtPlus;
    TH1D*       fHistPhiWgtMinus;
    TH1D*       fHistPhiWgtAll;
    TH1D*       fHistPhiFlatPlus;
    TH1D*       fHistPhiFlatMinus;
    TH1D*       fHistPhiFlatAll;
    TH1D*       fHistPhi;
    TH1D*       fHistPhiWgt;
    TH1D*       fHistPhiFlat;
   // flow (events)
    TH1F*       fHistPsi;
    TH1F*       fHistPsiSubCorr;
    TH1F*       fHistPsiSubCorrDiff;
    TH1F*       fHistPsiDiff;
    TH1F*       fHistMult;
    TH1F*       fHistQnorm;
   // flow (tracks)
    TH1F*       fHistPhiCorr;
    TProfile2D* fHistvObs2D;
    TProfile*   fHistvObsEta;
    TProfile*   fHistvObsPt;
    TH2D*       fHistv2D;
    TH1D*       fHistvEta;
    TH1D*       fHistvPt;
   // flow (v0s)
    TH1F*       fHistV0PhiCorr;   
    TProfile2D* fHistV0vObs2D;   
    TProfile*   fHistV0vObsEta;  
    TProfile*   fHistV0vObsPt;   
    TH2D*	fHistV0v2D;	 
    TH1D*	fHistV0vEta;	 
    TH1D*	fHistV0vPt;	 
   // flow (v0s sidebands)
    TProfile*   fHistV0sbvObsEtaSx ;  
    TProfile*   fHistV0sbvObsPtSx ;   
    TProfile*   fHistV0sbvObsEtaDx ;  
    TProfile*   fHistV0sbvObsPtDx ;   
    TH1F*       fHistV0sbPhiCorr ;   
    TProfile2D* fHistV0sbvObs2D ;   
    TProfile*   fHistV0sbvObsEta ;  
    TProfile*   fHistV0sbvObsPt ;   
    TH2D*	fHistV0sbv2D ;      
    TH1D*	fHistV0sbvEta ;     
    TH1D*	fHistV0sbvPt ;      
   // check (tracks used for R.P.)
    TH1F*       fHistYieldPt ;
    TH3F*       fHistEtaPtPhi3D ;
    TH2D*       fHistYield2D ;
    TH1F*       fHistDcaGlob ;
   // check (tracks excluded)
    TH1F*	fHistYieldPtout;
    TH3F*	fHistEtaPtPhi3Dout ;
    TH2D*	fHistYield2Dout ;
    TH1F*	fHistDcaGlobout ;
  };

// for each selection

  struct histFulls;	
  friend struct histFulls;
  struct histFulls 
  {
   TH1F*     fHistBayPidMult;
  // flow (events)
   TProfile* fHistCos;
   TH1F*     fHistRes;
   TProfile* fHistvObs;
   TH1D*     fHistv;
   TProfile* fHistV0vObs;
   TProfile* fHistV0sbvObsSx;
   TProfile* fHistV0sbvObsDx;
   TH1D*     fHistV0v;
  // wgt, evts, trks, v0s (as defined above)
   struct histFullHars  fHistFullHar[Flow::nHars];
  };
  struct histFulls fHistFull[Flow::nSels];                     //!

  ClassDef(AliFlowAnalyser,0)              // macro for rootcint
};

#endif



// lame = not productive; poorly designed; uncool ...
