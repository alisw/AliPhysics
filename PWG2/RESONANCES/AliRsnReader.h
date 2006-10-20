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
	Bool_t     EffSim(Int_t pdg, Double_t pt, Double_t eta, Double_t z);
	TTree*     GetEvents() const  						                  {return fEvents;}
	Double_t*  GetPIDprobabilities(AliRsnDaughter track) const;
	void       Identify(AliRsnDaughter &track);
	TTree*     ReadParticles(const char *path, Option_t *option="");
	TTree*     ReadTracks(const char *path, Option_t *option="R");
	TTree*     ReadTracksAndParticles(const char *path, Option_t *option="R");
	void       SetEfficiencyHistogram(AliPID::EParticleType type, TH3D *h, Bool_t pos=kTRUE);
	void       SetMaxRadius(Double_t value)				                  {fMaxRadius=value;}
	void       SetPIDMethod(AliRsnReader::EPIDMethod pm)                    {fPIDMethod=pm;}
	void       SetPriorProbabilities(Double_t *prior);
	void       SetPriorProbability(AliPID::EParticleType type, Double_t p);
	void       SetProbabilityThreshold(Double_t p)	                      {fProbThreshold=p;}
	void       SetPtLimit4PID(Double_t value) 			                  {fPtLimit4PID=value;}
	void       SetResolutionPt(TH1D *histo)                                 {fResPt=histo;}
	void       SetResolutionLambda(TH1D *histo)                             {fResLambda=histo;}
	void       SetResolutionPhi(TH1D *histo)                                {fResPhi=histo;}
	void       SmearMomentum(Double_t &px, Double_t &py, Double_t &pz);
	
protected:

	AliPID::EParticleType FindType(Int_t pdg);
	AliRunLoader*         OpenRunLoader(const char *path);

	Bool_t     fDoEffSim;                 //  set to true to introduce efficiency effects
	Bool_t     fDoSmearing;               //  when selecting only particles, set this to true 
	                                      //  to smear particle momentum
	
	EPIDMethod fPIDMethod;                //  PID method
	Double_t   fPrior[AliPID::kSPECIES];  //  prior probabilities (in case of REAL pid)
	Double_t   fPtLimit4PID; 			  //  maximum transverse momentum to accept realistic PID
	Double_t   fProbThreshold;			  //  minimum required probability to accept realistic PID
	Double_t   fMaxRadius;                //  maximum allowed distance from primary vertex
	
	TTree     *fEvents;                   //! tree of read events
	TH3D      *fEffPos[AliPID::kSPECIES]; //  tracking efficiency (in pt, eta and Zv) for positive particles
	TH3D      *fEffNeg[AliPID::kSPECIES]; //  tracking efficiency (in pt, eta and Zv) for negative particles
	TH1D      *fResPt;                    //  relative Pt resolution
	TH1D      *fResLambda;                //  relative Lambda resolution
	TH1D      *fResPhi;                   //  relative Phi resolution
	
	// Rsn event reader implementation
	ClassDef(AliRsnReader,1)
};

#endif

