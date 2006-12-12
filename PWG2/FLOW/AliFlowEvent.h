//////////////////////////////////////////////////////////////////////
//
// $Id$
//
// Author: Emanuele Simili
//
//////////////////////////////////////////////////////////////////////
//
// Description: event format fitted to flow study, adapted from STAR 
// Original Authors:               Raimond Snellings & Art Poskanzer
//
//////////////////////////////////////////////////////////////////////

#ifndef AliFlowEvent_h
#define AliFlowEvent_h

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include "TVector.h"
#include "TVector2.h"
#include "TVector3.h"
#include "TObject.h"
#include "TRandom.h"
#include "TObjArray.h"
#include <TROOT.h>

#include "AliFlowTrack.h"
#include "AliFlowV0.h"
#include "AliFlowSelection.h"
#include "AliFlowConstants.h"

class AliFlowTrack ;
class AliFlowV0 ;
class AliFlowSelection ;
class Flow ;

class AliFlowEvent : public TNamed {

public: 									   
										   
  AliFlowEvent(Int_t length = 1000) ; 						   
  virtual        ~AliFlowEvent();						   
										   
 // Gets
  Int_t          EventID() 		const;					// Returns ID of the event
  Int_t          RunID() 		const;					// Returns ID of the run
  Long_t         L0TriggerWord()	const;					// Returns L0 trigger word
  Int_t          ZDCpart()	 	const;  			        // Returns estimated number of participants by the ZDC
  Float_t        ZDCenergy(Int_t npem)	const;	  			        // Returns reconstructed energy in the neutron(1), proton(2), em(3) ZDC
  void           VertexPos(Float_t vtx[3]) const;				// Returns primary vertex position
  TVector3       VertexPos() 		const;					// Returns primary vertex position as a TVector3
  UInt_t         Centrality() ;							// Returns centrality bin (based on MultEta() )
  UInt_t         OrigMult() 		const;					// Returns the original number of tracks (maybe some were trown away)
  Int_t    	 V0Mult() 		const;					// Returns the number of V0s 
  Int_t          FlowEventMult() 	const;					// Returns number of tracks stored in the event
  Int_t          UncorrNegMult(Float_t eta = Flow::fEtaGood) const ;	        // Returns number of - tracks in eta (-eta;eta)
  Int_t          UncorrPosMult(Float_t eta = Flow::fEtaGood) const ;	        // Returns number of + tracks in eta (-eta;eta)
  Int_t          MultEta() ;							// Returns multiplicity in |eta|<Flow::fEetaMid
  Double_t       CenterOfMassEnergy() 	const;					// Returns center of mass energy (5.5 TeV)
  Double_t       MagneticField() 	const;					// Returns magnetic field value
  Short_t        BeamMassNumberEast() 	const;					// Returns beam mass (Pb = 208)
  Short_t        BeamMassNumberWest() 	const;					// Returns beam mass (Pb = 208)

									  
 // Sets
  void  	 SetEventID(const Int_t&);
  void  	 SetRunID(const Int_t&);
  void  	 SetCenterOfMassEnergy(const Double_t&);
  void  	 SetMagneticField(const Double_t&);
  void           SetCentrality() ; 
  void           SetCentrality(Int_t cent) ; 
  void  	 SetBeamMassNumberEast(const Short_t&);
  void  	 SetBeamMassNumberWest(const Short_t&);
  void  	 SetOrigMult(const UInt_t&);
  void  	 SetL0TriggerWord(const Long_t&);
  void  	 SetVertexPos(Float_t v1=0.,Float_t v2=0.,Float_t v3=0.);
  void  	 SetZDCpart(Int_t zdcp);
  void  	 SetZDCenergy(Float_t n, Float_t p, Float_t em);

//  // new
//   void 	      	 SetExtPsi(Int_t harN, Float_t psi=0.) ;
//   void 	      	 SetExtRes(Int_t harN, Float_t res=0.) ;
//   Float_t 	 ExtPsi(Int_t harN) 	const ; 				// external RP angle
//   Float_t 	 ExtRes(Int_t harN)  	const ; 				// external RP resolution

 // Arrays
  TObjArray*     TrackCollection() 	const;					// Returns a pointer to the TObjArray of AliFlowTrack
  TObjArray*     V0Collection() 	const;					// Returns a pointer to the TObjArray of AliFlowV0

 // Calculations
  Int_t          Mult(AliFlowSelection*); 					// Returns Multiplicity of tracks selected for the event plane
  TVector2       Q(AliFlowSelection*);  					// Returns Event plane vector
  Float_t        q(AliFlowSelection*);					 	// Returns Magnitude of normalized Q vector without pt or eta weighting
  Float_t        MeanPt(AliFlowSelection*);					// Returns Mean pt of tracks selected for the event plane
  Float_t        Psi(AliFlowSelection*);  					// Returns Event plane angle
  Double_t       G_New(AliFlowSelection* pFlowSelect,Double_t Zx,Double_t Zy);	// Generating function for the new cumulant method (eq.3 in the Practical Guide)
  Double_t       G_Old(AliFlowSelection* pFlowSelect,Double_t Zx,Double_t Zy);	// Generating function for the old cumulant method (if expanded in Taylor series, one recovers G_New() in new new cumulant method)
  Double_t       SumWeightSquare(AliFlowSelection* pFlowSelect);  		// Returns Sum of weights^2
  Double_t       WgtMult_q4(AliFlowSelection* pFlowSelect);			// old comulants
  Double_t       WgtMult_q6(AliFlowSelection* pFlowSelect);			// old comulants
  TVector2       NormQ(AliFlowSelection* pFlowSelect);				// Returns normalized Q = Q/sqrt(weights^2++)
 // -
  void  	 SetSelections(AliFlowSelection* pFlowSelect) ;			// Sets the tracks selection for R.P. calculations (see static cuts in AliFlowSelection class)
  void  	 RandomShuffle() ;		  	 			// Randomly re-shuffles the ObjArray of tracks
  void  	 MakeSubEvents() ;						// Makes sub-events, eta based (if EtaSubs()) or Random (otherwise)
  void  	 MakeRndSubEvents() ;						// Makes random sub-events
  void  	 MakeEtaSubEvents() ;						// Makes eta sub-events
  void 	 	 SetPids() ;							// Re-sets the tracks P.id. (using the current fBayesianCs[] array)
 // -
  void  	 MakeAll() ;							// In just one loop, makes all the calculaton (Q, psi, mult) basing on the selection. 

 // Weights & settings
  Float_t        PtWgtSaturation() 	const;					// Returns saturation value for pt weighting
  Bool_t         PtWgt() 		const;					// Returns flag for pt weighting
  Bool_t         EtaWgt() 		const;					// Returns flag for eta weighting for odd harmonics
  Bool_t         FirstLastPhiWgt() 	const;					// Returns flag for using z of first and last points for phi weights (TPC +/-)
  Bool_t         OnePhiWgt() 		const;					// Returns flag for using just one phi weight
  Bool_t         NoWgt() 		const;					// returns kTRUE if weight are NOT used
  Double_t       PhiWeight(Int_t selN,Int_t harN,AliFlowTrack* pFlowTrack) const ; 	// Returns PhiWeightRaw()*Weight()  
  Double_t       PhiWeightRaw(Int_t selN,Int_t harN,AliFlowTrack* pFlowTrack) const ;	// Returns weights for making the R.P. isotropic in the lab
  Double_t       Weight(Int_t selN,Int_t harN,AliFlowTrack* pFlowTrack) const ; 	// Returns weights for enhancing the resolution (+/-Sign(eta) for odd harmonics)
  void  	 Bayesian(Double_t bayes[Flow::nPid]) ; 			// Returns the stored particles' abundances
  TVector   	 Bayesian() ; 	 						// Returns the stored particles' abundances as a TVector
 // -
  static void 	 SetPtWgt(Bool_t PtWgt = kTRUE);
  static void 	 SetEtaWgt(Bool_t EtaWgt = kTRUE);
  static void 	 SetOnePhiWgt();
  static void 	 SetFirstLastPhiWgt();
  static void 	 SetNoWgt(Bool_t nowgt = kTRUE) ;
 // -
  void  	 SetBayesian(Double_t bayes[Flow::nPid]) ; 			// Set the Bayesian vector of particles' abundances
#ifndef __CINT__		
  void     	 SetPhiWeight(const Flow::PhiWgt_t &pPhiWgt);		    	// Fills Weights from Arrays (from file: flowPhiWgt.hist.root)
  void  	 SetPhiWeightPlus(const Flow::PhiWgt_t &pPhiWgtPlus);
  void  	 SetPhiWeightMinus(const Flow::PhiWgt_t &pPhiWgtMinus);
  void  	 SetPhiWeightCross(const Flow::PhiWgt_t &pPhiWgtCross);
#endif

 // Analysis flags
  Bool_t         EtaSubs() const;						// Returns flag for eta sub-events
  static void 	 SetEtaSubs(Bool_t etasub = kTRUE) ;				// Sets the flag for eta sub-events
  void  	 PrintFlagList() const ;					// Prints a summary of the event's flag


private:

 // Data Members
  Int_t               fEventID;                                  // ID of the event
  Int_t               fRunID;                                    // ID of the run
  UInt_t              fOrigMult;                                 // Original number of tracks
  Long_t              fL0TriggerWord;                            // L0 trigger word
  Int_t               fZDCpart;                                  // ZDC estimated number of participants 
  Float_t             fZDCenergy[3];                             // ZDC reconstructed energy [neutron,proton,em]
  Float_t 	      fVertexPos[3];                             // primary vertex position
  Int_t               fCentrality;                               //! Centrality Class (calculated from mult.)

 // extension
  //Float_t 	      fExtPsi[Flow::nHars] ;			 // external RP angle (should be an input)
  //Float_t 	      fExtRes[Flow::nHars] ;			 // external RP resolution (should be an input as well)

 // Tracks & V0s
  TObjArray*	      fTrackCollection ;			 // collection of Flow Tracks
  TObjArray*	      fV0Collection ;				 // collection of Flow V0s

 // Weights
  Flow::PhiWgt_t      fPhiWgt;                                   //! flattening weights (single hist)
  Flow::PhiWgt_t      fPhiWgtPlus;                               //! flattening weights (3 hist) - plus Z
  Flow::PhiWgt_t      fPhiWgtMinus;                              //! flattening weights (3 hist) - minus Z
  Flow::PhiWgt_t      fPhiWgtCross;                              //! flattening weights (3 hist) - cross Z
  //Double_t     	fBayesianCs[Flow::nPid] ;        	 //! expected particles abundance (see Bayesian P.Id.)

 // Weighting & Settings
  static Bool_t       fPtWgt;                                    //! flag for pt weighting
  static Bool_t       fEtaWgt;                                   //! flag for eta weighting for odd harmonics
  static Bool_t       fOnePhiWgt;                                //! flag for phi weights (just one hist)
  static Bool_t       fNoWgt;                          		 //! No Weights (Wgt == 1)
  static Bool_t       fEtaSubs;                                  //! Flag for making Eta Subevents

#ifndef __CINT__
 // shortcuts (to speed up the execution)
  Bool_t   fDone ;						//! flag setted kTRUE when the loop is done
  TVector2 fQ[Flow::nSels][Flow::nHars];			//! flow vector
  UInt_t   fMult[Flow::nSels][Flow::nHars];                  	//! multiplicity
  Float_t  fSumOfWeightSqr[Flow::nSels][Flow::nHars];           //! Sqrt(Sum(wgt)) ~ Sqrt(Mult)
  TVector2 fQSub[Flow::nSubs][Flow::nSels][Flow::nHars];      	//! flow vector subs
  UInt_t   fMultSub[Flow::nSubs][Flow::nSels][Flow::nHars];   	//! multiplicity subs
#endif /*__CINT__*/

  ClassDef(AliFlowEvent,2) ;                    // macro for rootcint
};

#endif
