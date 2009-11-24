#ifndef ALIANALYSISCENTRALEXTRAPOLATE_H
#define ALIANALYSISCENTRALEXTRAPOLATE_H
/*
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved.
 * See cxx source for full Copyright notice
 * $Id$
 */

//  *********************************************************
//  * pt spectra extrapolation in the 0. - 0.2 region using *
//  * Boltzmann-Gibbs Blast Wave model or Tsallis Blast     *
//  * Wave model for azimuthal isotropic  expansion in      *
//  * highly central collisions analysis                    *
//  * author: Cristian Andrei                               *
//  *         acristian@niham.nipne.ro                      *
//  * *******************************************************


#include "TObject.h"

class TList;
class TString;

class AliAnalysisCentralExtrapolate : public TObject {
	
public:
	AliAnalysisCentralExtrapolate(const char *name="AliAnalysisCentralExtrapolate");
	virtual ~AliAnalysisCentralExtrapolate();
	
	void SetInputList(TList* const list) {fInputList = list;} //get the input from the task
	
	void SetParticle(TString part) {fPartType = part;} //set the particle type to be processed
                                                       //suported: kPi*, kK* and kProton

	void ApplyEff(); //apply the eff correction
	
	void TsallisFit();  //fit the corrected spectrum with TBW model
	void BoltzmannFit(); //fit the corrected spectrum with BGBW model
	
	TList *GetOutputList() const {return fResultsList;} //return the results
	
private:
	AliAnalysisCentralExtrapolate(const AliAnalysisCentralExtrapolate& ref);
	AliAnalysisCentralExtrapolate& operator=(const AliAnalysisCentralExtrapolate& ref);
	
	TString fPartType; //can be kPi*, kK* or kProton
	
	TList *fInputList; //the input list = the output of the analysis Task
	
	TList *fResultsList; //List containing the results: corrected, normalized and extrapolated spectra
	
	ClassDef(AliAnalysisCentralExtrapolate, 1); 
};

#endif
