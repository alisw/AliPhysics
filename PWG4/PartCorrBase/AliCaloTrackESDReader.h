#ifndef ALICALOTRACKESDREADER_H
#define ALICALOTRACKESDREADER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
/* $Id:  $ */

//_________________________________________________________________________
// Class for reading data (ESDs) in order to do prompt gamma 
//  or other particle identification and correlations
// 

//*-- Author: Gustavo Conesa (INFN-LNF)

// --- ROOT system --- 

// --- AliRoot system ---
#include "AliCaloTrackReader.h" 

class AliCaloTrackESDReader : public AliCaloTrackReader {
	
public: 
	
	AliCaloTrackESDReader() ; // ctor
	AliCaloTrackESDReader(const AliCaloTrackESDReader & g) ; // cpy ctor
	//AliCaloTrackESDReader & operator = (const AliCaloTrackESDReader & g) ;//cpy assignment
	virtual ~AliCaloTrackESDReader() {;} //virtual dtor
	
	void FillInputCTS  () ;
	void FillInputEMCAL() ;
	void FillInputPHOS () ;  
	void FillInputEMCALCells() ;
	void FillInputPHOSCells() ;
	
	void GetVertex(Double_t v[3]) const ;
	
	AliVEvent*  GetInputEvent() const {return GetESD();}
	void SetInputEvent(TObject* esd, TObject* aod, TObject* mc) ; 
	
	ClassDef(AliCaloTrackESDReader,1)
} ;


#endif //ALICALOTRACKESDREADER_H



