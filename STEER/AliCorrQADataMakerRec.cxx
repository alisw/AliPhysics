/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/


/* $Id: AliCorrQADataMakerRec.cxx 27570 2008-07-24 21:49:27Z cvetan $ */

/*
  Produces the data needed to calculate the quality assurance. 
  All data must be mergeable objects.
  Y. Schutz CERN July 2007
*/

// --- ROOT system ---
#include <TClonesArray.h>
#include <TFile.h> 
#include <TH1F.h> 
#include <TH1I.h> 
#include <TH2F.h> 
#include <TNtupleD.h>
#include <TParameter.h>
#include <TMath.h> 

// --- Standard library ---

// --- AliRoot header files ---
#include "AliLog.h"
#include "AliCorrQADataMakerRec.h"
#include "AliQAChecker.h"

ClassImp(AliCorrQADataMakerRec)
           
//____________________________________________________________________________ 
AliCorrQADataMakerRec::AliCorrQADataMakerRec(AliQADataMaker ** qadm ) : 
AliQADataMakerRec(AliQAv1::GetDetName(AliQAv1::kCORR), "Corr Quality Assurance Data Maker"),
  fMaxRawVar(0),  
  fqadm(qadm)
{
  // ctor
  fCorrNt = new TNtupleD *[AliRecoParam::kNSpecies] ; 
  for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) 
    fCorrNt[specie] = NULL ; 
}

//____________________________________________________________________________ 
AliCorrQADataMakerRec::AliCorrQADataMakerRec(const AliCorrQADataMakerRec& qadm) :
  AliQADataMakerRec(),
  fMaxRawVar(qadm.fMaxRawVar), 
  fqadm(qadm.fqadm)
{
  //copy ctor 
  SetName((const char*)qadm.GetName()) ; 
  SetTitle((const char*)qadm.GetTitle()); 
}

//__________________________________________________________________
AliCorrQADataMakerRec& AliCorrQADataMakerRec::operator = (const AliCorrQADataMakerRec& qadm )
{
  // assign operator.
  this->~AliCorrQADataMakerRec();
  new(this) AliCorrQADataMakerRec(qadm);
  return *this;
}
   
//____________________________________________________________________________ 
AliCorrQADataMakerRec::~AliCorrQADataMakerRec()  
{
  // dtor only destroy the ntuple 
  if ( fCorrNt ) {
    for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) {
      if ( fCorrNt[specie] != NULL ) 
        delete fCorrNt[specie] ; 
    }
		delete[] fCorrNt ; 
    fCorrNt = 0x0;
	}  
}
  
//____________________________________________________________________________ 
void AliCorrQADataMakerRec::EndOfDetectorCycle(AliQAv1::TASKINDEX_t task, TObjArray ** /*list*/)
{
  //Detector specific actions at end of cycle
  // do the QA checking
  if (task == AliQAv1::kRAWS) {
     AliQAChecker::Instance()->Run(AliQAv1::kCORR, task, fCorrNt) ; 
  }
}

//____________________________________________________________________________ 
void AliCorrQADataMakerRec::InitESDs()
{
  //Create histograms to controll ESD

  AliInfo("TO BE IMPLEMENTED") ; 
}

//____________________________________________________________________________ 
void AliCorrQADataMakerRec::InitRecPoints()
{
  // create Reconstructed Points histograms in RecPoints subdir

  AliInfo("TO BE IMPLEMENTED") ; 
}

//____________________________________________________________________________ 
void AliCorrQADataMakerRec::InitRaws()
{
  // createa ntuple taking all the parameters declared by detectors
  if (fCorrNt) 
    return ; 
  delete fRawsQAList ; // not needed for the time being 
  fRawsQAList = NULL ; 
  TString varlist("") ;
  for ( Int_t detIndex = 0 ; detIndex < AliQAv1::kNDET ; detIndex++ ) {
    AliQADataMaker * qadm = fqadm[detIndex] ; 
    if ( ! qadm ) 
      continue ;
    TList * list = qadm->GetParameterList() ; 
    if (list) {
      TIter next(list) ; 
      TParameter<double> * p ; 
      while ( (p = static_cast<TParameter<double>*>(next()) ) ) {
        varlist.Append(p->GetName()) ; 
        varlist.Append(":") ; 
        fMaxRawVar++ ; 
      }
    }
  }
  varlist = varlist.Strip(TString::kTrailing, ':') ; 
  if (fMaxRawVar == 0) { 
    AliWarning("NTUPLE not created") ; 
  } else {
    char * name = Form("%s_%s", AliQAv1::GetQACorrName(), AliRecoParam::GetEventSpecieName(fEventSpecie)) ; 
    fCorrNt[AliRecoParam::AConvert(fEventSpecie)] = new TNtupleD(name, "Raws data correlation among detectors", varlist.Data()) ;  
  }  
}

//____________________________________________________________________________
void AliCorrQADataMakerRec::MakeESDs(AliESDEvent * /*esd*/)
{
  // make QA data from ESDs

  AliInfo("TO BE IMPLEMENTED") ; 

}

//____________________________________________________________________________
void AliCorrQADataMakerRec::MakeRaws(AliRawReader *)
{
  //Fill prepared histograms with Raw digit properties
  
  if ( ! fCorrNt[AliRecoParam::AConvert(fEventSpecie)])
      InitRaws() ; 
  
  if ( fMaxRawVar > 0 ) {
    const Int_t kSize = fMaxRawVar ; 
    Double_t  *varvalue = new Double_t[kSize] ;
    Int_t index = 0 ;
    for ( Int_t detIndex = 0 ; detIndex < AliQAv1::kNDET ; detIndex++ ) {
      AliQADataMaker * qadm = fqadm[detIndex] ; 
      if ( ! qadm ) 
        continue ;
      TList * list = qadm->GetParameterList() ; 
      TIter next(list) ; 
      TParameter<double> * p ; 
      while ( (p = static_cast<TParameter<double>*>(next()) ) ) {
        varvalue[index++] = p->GetVal() ; 
      }
    }
    static_cast<TNtupleD*>(fCorrNt[AliRecoParam::AConvert(fEventSpecie)])->Fill(varvalue);
    delete [] varvalue;
  }
}

//____________________________________________________________________________
void AliCorrQADataMakerRec::MakeRecPoints(TTree * /*clustersTree*/)
{
  AliInfo("TO BE IMPLEMENTED") ; 
}

//____________________________________________________________________________ 
void AliCorrQADataMakerRec::StartOfDetectorCycle()
{
  //Detector specific actions at start of cycle  

}
