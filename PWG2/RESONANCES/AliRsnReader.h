/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               *
 **************************************************************************/

//-------------------------------------------------------------------------
//                      Class AliRsnEventReader
//             
//   Reader for conversion of ESD or Kinematics output into AliRsnEvent
//   .....
//   .....
//   .....
//   .....
// 
// author: A. Pulvirenti             (email: alberto.pulvirenti@ct.infn.it)
//-------------------------------------------------------------------------

#ifndef AliRSNREADER_H
#define AliRSNREADER_H

#include "AliPID.h"

class TH3D;
class TH1D;
class TOrdCollection;
class TTree;
class AliRunLoader;

class AliRsnReader : public TObject
{
public:

	enum EPIDMethod { 
		kNoPID = 0,                  // no PID performed
		kPerfectPID = 1,             // use Kinematics to simulate a 100% PID efficiency
		kESDPID = 2,                 // use experimental ESD weights
		kPerfectPIDwithRecSign = 3   // get particle type from Kine and charge sign from reconstruction
	};

	           AliRsnReader();
		  	   AliRsnReader(const AliRsnReader& copy);
			   AliRsnReader& operator=(const AliRsnReader& copy);
	virtual   ~AliRsnReader() {Clear("DELTREE");}
	
	void       Clear(Option_t *option = "");
	TTree*     GetEvents() const {return fEvents;}
	Double_t*  GetPIDprobabilities(AliRsnDaughter track) const;
	void       Identify(AliRsnDaughter &track);
	TTree*     ReadTracks(const char *path, Option_t *option="R");
	TTree*     ReadTracksAndParticles(const char *path, Option_t *option="R");
	void       SetMaxRadius(Double_t value) {fMaxRadius=value;}
	void       SetPIDMethod(AliRsnReader::EPIDMethod pm) {fPIDMethod=pm;}
	void       SetPriorProbabilities(Double_t *prior);
	void       SetPriorProbability(AliPID::EParticleType type, Double_t p);
	void       SetProbabilityThreshold(Double_t p) {fProbThreshold=p;}
	void       SetPtLimit4PID(Double_t value) {fPtLimit4PID=value;}
	void       SetUseKineInfo(Bool_t yesno = kTRUE) {fUseKineInfo=yesno;}
	
protected:

	AliPID::EParticleType FindType(Int_t pdg);
	AliRunLoader*         OpenRunLoader(const char *path);

	EPIDMethod fPIDMethod;                //  PID method
	Double_t   fPrior[AliPID::kSPECIES];  //  prior probabilities (in case of REAL pid)
	Double_t   fPtLimit4PID; 			  //  maximum transverse momentum to accept realistic PID
	Double_t   fProbThreshold;			  //  minimum required probability to accept realistic PID
	Double_t   fMaxRadius;                //  maximum allowed distance from primary vertex
	
	Bool_t     fUseKineInfo;              //  set to TRUE to fill the fields 'fTruePDG', 'fMother' and 'fMotherPDG'
	                                      //  of the AliRsnDaughter objects returned in the output
	
	TTree     *fEvents;                   //! tree of read events
	
	// Rsn event reader implementation
	ClassDef(AliRsnReader,1)
};

#endif

