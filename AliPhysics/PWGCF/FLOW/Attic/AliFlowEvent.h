//////////////////////////////////////////////////////////////////////
//
// $Id: AliFlowEvent.h 18618 2007-05-16 15:38:22Z snelling $
//
// Author: Emanuele Simili
//
//////////////////////////////////////////////////////////////////////
//
// Description: event format fitted to flow study, adapted from STAR.
//  The AliFlowEvent object stores global event variables, tracks and  
//  v0s, and provides the method to calculate Flow observables (Q vector, 
//  R.P. angle, phi weights). It needs an AliFlowSelection object.
//
// Original Authors:               Raimond Snellings & Art Poskanzer
//
//////////////////////////////////////////////////////////////////////

#ifndef AliFlowEvent_h
#define AliFlowEvent_h

#include "TVector.h"
#include "TVector2.h"
#include "TVector3.h"
#include "TClonesArray.h"
#include "TNamed.h"

#include "AliFlowConstants.h"

class AliFlowTrack ;
class AliFlowV0 ;
class AliFlowSelection ;

class AliFlowEvent : public TNamed {

public: 									   
										   
  AliFlowEvent() ; 				   // default constructor (to read events)		   
  AliFlowEvent(Int_t length) ;		           // (new) constructor (to create events) 						   
  virtual        ~AliFlowEvent();						   
										   
 // Arrays
  TClonesArray*     TrackCollection() const        { return fTrackCollection; }  // Returns a pointer to the TClonesArray of AliFlowTrack
  TClonesArray*     V0Collection() const	   { return fV0Collection; }     // Returns a pointer to the TClonesArray of AliFlowV0

 // Calculations
  void              MakeAll() ; 						 // In just one loop, makes all the calculaton (Q, psi, mult) basing on the selection. 
  Int_t             Mult(AliFlowSelection* pFlowSelect) const ;			 // Returns Multiplicity of tracks selected for the event plane
  Float_t           MeanPt(AliFlowSelection* pFlowSelect) const ;  		 // Returns Mean pt of tracks selected for the event plane
  TVector2          Q(AliFlowSelection* pFlowSelect) const ;			 // Returns Event plane vector
  TVector2          NormQ(AliFlowSelection* pFlowSelect) const ;		 // Returns normalized Q = Q/sqrt(weights^2++)
  Float_t           OldQ(AliFlowSelection* pFlowSelect) const ;			 // Returns Magnitude of normalized Q vector without pt or eta weighting
  Float_t           Psi(AliFlowSelection* pFlowSelect) const ;			 // Returns Event plane angle
  Double_t          SumWeightSquare(AliFlowSelection* pFlowSelect) const ;	 // Returns Sum of weights^2
  Double_t          WgtMultQ4(AliFlowSelection* pFlowSelect) const ;		 // old comulants
  Double_t          WgtMultQ6(AliFlowSelection* pFlowSelect) const ;		 // old comulants
  void              SetSelections(AliFlowSelection* pFlowSelect) ;		 // Sets the tracks selection for R.P. calculations (see static cuts in AliFlowSelection class)
  void              MakeSubEvents() const ;					 // Makes sub-events, eta based (if EtaSubs()) or Random (otherwise)
  void              MakeRndSubEvents() const ;					 // Makes random sub-events
  void              MakeEtaSubEvents() const ;					 // Makes eta sub-events
  void              MakeChrSubEvents() const ;				 	 // Makes charged (+/-) sub-events
  void 	            SetPids() ; 						 // Re-sets the tracks P.id. (using the current fBayesianCs[] array)
  void              RandomShuffle() ;						 // Randomly re-shuffles the ObjArray of tracks
  Double_t          NewG(AliFlowSelection* pFlowSelect,Double_t Zx,Double_t Zy) const ;  // Generating function for the new cumulant method (eq.3 in the Practical Guide)
  Double_t          OldG(AliFlowSelection* pFlowSelect,Double_t Zx,Double_t Zy) const ;  // Generating function for the old cumulant method (if expanded in Taylor series, one recovers G_New() in new new cumulant method)
  Double_t          PhiWeight(Int_t selN,Int_t harN,AliFlowTrack* pFlowTrack) const ;	 // Returns PhiWeightRaw()*Weight()  
  Double_t          PhiWeightRaw(Int_t selN,Int_t harN,AliFlowTrack* pFlowTrack) const ; // Returns weights for making the R.P. isotropic in the lab
  Double_t          Weight(Int_t selN,Int_t harN,AliFlowTrack* pFlowTrack) const ;	 // Returns weights for enhancing the resolution (+/-Sign(eta) for odd harmonics)
  TVector           Bayesian() const ;						         // Returns the stored particle abundances as a TVector
  void 		    Bayesian(Double_t bayes[AliFlowConstants::kPid]) const { for(Int_t i=0;i<AliFlowConstants::kPid;i++) { bayes[i] = AliFlowConstants::fgBayesian[i] ; } } // Returns the stored particle abundances

 // Analysis flags
  void       PrintFlagList() const ;					 	      // Prints a summary of the event's flag
  Float_t    PtWgtSaturation() const  { return AliFlowConstants::fgPtWgtSaturation; } // Returns saturation value for pt weighting
  Bool_t     PtWgt() const	      { return fgPtWgt; }			      // Returns flag for pt weighting
  Bool_t     EtaWgt() const	      { return fgEtaWgt; }			      // Returns flag for eta weighting for odd harmonics
  Bool_t     FirstLastPhiWgt() const  { return !fgOnePhiWgt ; }			      // Returns flag for using z of first and last points for phi weights (TPC +/-)
  Bool_t     OnePhiWgt() const	      { return fgOnePhiWgt ; }  		      // Returns flag for using just one phi weight
  Bool_t     NoWgt() const	      { return fgNoWgt; }       		      // returns kTRUE if weight are NOT used
  Bool_t     CustomRespFunc() const   { return fgCustomRespFunc ; }
  Bool_t     EtaSubs() const	      { if(fgEtaSubs == 1)  { return kTRUE ; } return kFALSE ; }   // Returns flag for charged (+/-) sub-events
  Bool_t     RndSubs() const	      { if(fgEtaSubs == 0)  { return kTRUE ; } return kFALSE ; }   // Returns flag for random sub-events
  Bool_t     ChrSubs() const          { if(fgEtaSubs == -1) { return kTRUE ; } return kFALSE ; }   // Returns flag for eta sub-events
  Int_t      Subs() const             { return fgEtaSubs ; }   			      // Returns flag for sub-events type (0 = random , 1 = eta , -1 = charged)

 // Gets
  Int_t      EventID() const		    { return fEventID; }			  // Returns ID of the event
  Int_t      RunID() const		    { return fRunID; }  			  // Returns ID of the run
  UInt_t     OrigMult() const		    { return fOrigMult; }			  // Returns the original number of tracks (maybe some were trown away)
  Long_t     L0TriggerWord() const	    { return fL0Trigger; }			  // Returns L0 trigger word
  Int_t      ZDCpart() const		    { return fZDCpart; }			  // Returns estimated number of participants by the ZDC
  Float_t    ZDCenergy(Int_t npem) const    { return fZDCenergy[npem]; }		  // Returns reconstructed energy in the neutron(1), proton(2), em(3) ZDC
  Int_t      V0Mult() const		    { return fV0Collection->GetEntries() ; }	  // Returns the number of V0s 
  Int_t      TrackMult() const  	    { return fTrackCollection->GetEntries() ; }   // Returns number of tracks stored in the event

  Int_t      UncorrNegMult(Float_t eta = AliFlowConstants::fgEtaGood) const ; // Returns number of - tracks in eta (-eta;eta)
  Int_t      UncorrPosMult(Float_t eta = AliFlowConstants::fgEtaGood) const ; // Returns number of + tracks in eta (-eta;eta)
  Int_t      MultEta() const ;						      // Returns multiplicity in |eta|<AliFlowConstants::fgEetaMid
  UInt_t     Centrality() ;						      // Returns centrality class, if not there, it sets it with MultEta()
  TVector3   VertexPos() const; 					      // Returns primary vertex position as a TVector3
  void       VertexPos(Float_t vtx[3]) const { for(Int_t ii=0;ii<3;ii++) { vtx[ii] = fVertexPos[ii] ; }	}                         // Returns primary vertex position
									  
  Float_t    ExtPsi(Int_t harN = 0) const    { if(harN<AliFlowConstants::kHars) { return fExtPsi[harN] ; } else { return 0. ; } } // external RP angle
  Float_t    ExtRes(Int_t harN = 0) const    { if(harN<AliFlowConstants::kHars) { return fExtRes[harN] ; } else { return 0. ; } } // external RP resolution

  Double_t   CenterOfMassEnergy() const	     { return AliFlowConstants::fgCenterOfMassEnergy ; }   // Returns center of mass energy (5.5 TeV)
  Double_t   MagneticField() const	     { return AliFlowConstants::fgMagneticField ; }	   // Returns magnetic field value
  Short_t    BeamMassNumberEast() const	     { return AliFlowConstants::fgBeamMassNumberEast ; }   // Returns beam mass (Pb = 208)
  Short_t    BeamMassNumberWest() const	     { return AliFlowConstants::fgBeamMassNumberWest ; }   // Returns beam mass (Pb = 208)

 // Sets
  void     SetEventID(const Int_t& id)                    	    { fEventID = id ; }  // Sets Event ID and the Event name (name = evtNumber_runId)
  void     SetRunID(Int_t id)				  	    { fRunID = id; }
  void     SetOrigMult(UInt_t tracks)			  	    { fOrigMult = tracks; }
  void     SetL0TriggerWord(Long_t trigger)		  	    { fL0Trigger = trigger; }
  void     SetZDCpart(Int_t zdcp)			  	    { fZDCpart = zdcp ; }
  void     SetZDCenergy(Float_t n, Float_t p, Float_t em) 	    { fZDCenergy[0] = n ; fZDCenergy[1] = p ; fZDCenergy[2] = em ; }
  void     SetVertexPos(Float_t v1=0.,Float_t v2=0.,Float_t v3=0.)  { fVertexPos[0] = v1 ; fVertexPos[1] = v2 ; fVertexPos[2] = v3 ; }

  void 	   SetExtPsi(Int_t harN=0, Float_t psi=0.) 	  	    { if(harN<AliFlowConstants::kHars) { fExtPsi[harN] = psi ; } }
  void 	   SetExtRes(Int_t harN=0, Float_t res=0.)  	  	    { if(harN<AliFlowConstants::kHars) { fExtRes[harN] = res ; } }

  void     SetCentrality(Int_t cent) 		          	    { fCentrality = cent ; } // Set the Centrality Classes to "cent"
  void     SetCentrality() ; 				  	    // Sets the Centrality Classes basing on Multiplicity at mid rapidity

  static void SetChrSubs()		    			    { fgEtaSubs = -1 ; }
  static void SetRndSubs() 	      	    			    { fgEtaSubs = 0 ; }
  static void SetEtaSubs() 	      	    			    { fgEtaSubs = 1 ; }
  static void SetOnePhiWgt()				      	    { fgOnePhiWgt = kTRUE ; }
  static void SetFirstLastPhiWgt()			      	    { fgOnePhiWgt = kFALSE ; }
  static void SetPtWgt(Bool_t ptWgt = kTRUE)		      	    { fgPtWgt = ptWgt; }
  static void SetEtaWgt(Bool_t etaWgt = kTRUE)  	      	    { fgEtaWgt = etaWgt ; }
  static void SetNoWgt(Bool_t nowgt = kTRUE) 		      	    { fgNoWgt = nowgt ; }  // still for odd harmonics: Wgt = +1 (positive Eta) or -1 (negative Eta)
  static void SetCustomRespFunc(Bool_t crf = kTRUE) 		    { fgCustomRespFunc = crf ; }
		
  void    SetBayesian(Double_t bayes[AliFlowConstants::kPid]) const { for(Int_t i=0;i<AliFlowConstants::kPid;i++) { AliFlowConstants::fgBayesian[i] = bayes[i] ; } } // Set the Bayesian vector of particle abundances
  void    SetMagneticField(const Double_t& mf) const  	      	    { AliFlowConstants::fgMagneticField = mf; }
  void    SetCenterOfMassEnergy(const Double_t& cms) const	    { AliFlowConstants::fgCenterOfMassEnergy = cms; }
  void    SetBeamMassNumberEast(const Short_t& bme) const	    { AliFlowConstants::fgBeamMassNumberEast = bme; }
  void    SetBeamMassNumberWest(const Short_t& bmw) const	    { AliFlowConstants::fgBeamMassNumberWest = bmw; }

// Fills Weights from Arrays (from file: flowPhiWgt.hist.root)
  void    SetPhiWeight(const AliFlowConstants::PhiWgt_t& pPhiWgt)	    { memcpy (fPhiWgt, pPhiWgt, sizeof(AliFlowConstants::PhiWgt_t)); }
  void    SetPhiWeightPlus(const AliFlowConstants::PhiWgt_t& pPhiWgtPlus)   { memcpy (fPhiWgtPlus,  pPhiWgtPlus,  sizeof(AliFlowConstants::PhiWgt_t)); }
  void    SetPhiWeightMinus(const AliFlowConstants::PhiWgt_t& pPhiWgtMinus) { memcpy (fPhiWgtMinus, pPhiWgtMinus, sizeof(AliFlowConstants::PhiWgt_t)); }
  void    SetPhiWeightCross(const AliFlowConstants::PhiWgt_t& pPhiWgtCross) { memcpy (fPhiWgtCross, pPhiWgtCross, sizeof(AliFlowConstants::PhiWgt_t)); }


private:

 // to make the code checker happy
  AliFlowEvent(const AliFlowEvent &flowEvent) ; 		 // Copy Constructor (dummy)
  AliFlowEvent &operator=(const AliFlowEvent &flowEvent) ; 	 // Assignment Operator (dummy)

 // Data Members
  Int_t               fEventID;                                  // ID of the event
  Int_t               fRunID;                                    // ID of the run
  UInt_t              fOrigMult;                                 // Original number of tracks
  Long_t              fL0Trigger;                                // Level 0 trigger 
  Int_t               fZDCpart;                                  // ZDC estimated number of participants 
  Float_t             fZDCenergy[3];                             // ZDC reconstructed energy [neutron,proton,em]
  Float_t 	      fVertexPos[3];                             // primary vertex position
  Int_t               fCentrality;                               //! Centrality Class (calculated from mult.)

 // extension
  Float_t	      fExtPsi[AliFlowConstants::kHars] ;  	 // external RP angle (should be an input)
  Float_t	      fExtRes[AliFlowConstants::kHars] ;  	 // external RP resolution (should be an input as well)

 // Tracks & V0s
  TClonesArray*	      fTrackCollection ;			 // collection of Flow Tracks
  TClonesArray*	      fV0Collection ;				 // collection of Flow V0s

 // Weights
  AliFlowConstants::PhiWgt_t      fPhiWgt;                       //! flattening weights (single hist)
  AliFlowConstants::PhiWgt_t      fPhiWgtPlus;                   //! flattening weights (3 hist) - plus Z
  AliFlowConstants::PhiWgt_t      fPhiWgtMinus;                  //! flattening weights (3 hist) - minus Z
  AliFlowConstants::PhiWgt_t      fPhiWgtCross;                  //! flattening weights (3 hist) - cross Z
  //Double_t     	fBayesianCs[AliFlowConstants::kPid] ;    //! expected particles abundance (see Bayesian P.Id.)

 // Weighting & Settings
  static Bool_t       fgPtWgt;                                   //! flag for pt weighting
  static Bool_t       fgEtaWgt;                                  //! flag for eta weighting for odd harmonics
  static Bool_t       fgOnePhiWgt;                               //! flag for phi weights (just one hist)
  static Bool_t       fgNoWgt;                          	 //! No Weights (Wgt == 1)
  static Bool_t       fgCustomRespFunc ;  		 	 //! A custom "detector response function" is used for P.Id
  static Int_t        fgEtaSubs;                                 //! Flag type of Sub-Events (0 = random , 1 = eta , -1 = charged)


 // shortcuts (to speed up the execution)
  Bool_t   fDone ;										//! flag setted kTRUE when the loop is done
  TVector2 fQ[AliFlowConstants::kSels][AliFlowConstants::kHars];				//! flow vector
  UInt_t   fMult[AliFlowConstants::kSels][AliFlowConstants::kHars];                  		//! multiplicity
  Float_t  fSumOfWeightSqr[AliFlowConstants::kSels][AliFlowConstants::kHars];          		//! Sqrt(Sum(wgt)) ~ Sqrt(Mult)
  TVector2 fQSub[AliFlowConstants::kSubs][AliFlowConstants::kSels][AliFlowConstants::kHars];    //! flow vector subs
  UInt_t   fMultSub[AliFlowConstants::kSubs][AliFlowConstants::kSels][AliFlowConstants::kHars]; //! multiplicity subs

  ClassDef(AliFlowEvent,2) ;                    // macro for rootcint
};

#endif
