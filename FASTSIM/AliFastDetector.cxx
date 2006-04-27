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

/* $Id$ */
//
// Base class for fast simulation of a detctor
// or a system of subdetectors.
// The detector response is described by resolution and efficiency.
// Author:
// Andreas Morsch
// andreas.morsch@cern.ch

#include "AliFastDetector.h"
#include "AliFastResponse.h"
#include "AliGeometry.h"

#include <TList.h>
#include <TIterator.h>
#include <TString.h>

ClassImp(AliFastDetector)
AliFastDetector::AliFastDetector()
{
// Default Constructor
    fName  = "FastDetector";
    fTitle = "Fast Detector Base Class";
    fLnkD  = 0;
    fLnkR  = 0;
    
    fResponses    = 0;
    fSubdetectors = 0;
}

AliFastDetector::AliFastDetector(char* Name, char* Title):
    TNamed(Name, Title)
{
// Constructor
    fSubdetectors = new TList();
    fResponses    = new TList();
}

AliFastDetector::AliFastDetector(const AliFastDetector & det)
    :TNamed(det)
{
// Copy constructor
    det.Copy(*this);
}

AliFastDetector::~AliFastDetector()
{
// Destructor
    delete fSubdetectors;
    delete fResponses;
}


void AliFastDetector::Init()
{
//
// Initialisation
//
    TIter nextRes(fResponses);
    AliFastResponse *res;
    //
    // Loop over responses  and initialize
    while((res = (AliFastResponse*)nextRes())) {
	res->Init();
    }  

    TIter nextDet(fSubdetectors);
    AliFastDetector *det;
    //
    // Loop over subdetectors  and initialize
    while((det = (AliFastDetector*)nextDet())) {
	det->Init();
    }  
    //
    TObject* obj;
    
    if ((obj = fResponses->FindObject("Efficiency")))
    {

	fEfficiency = (AliFastResponse*) obj;
	printf("Detector %s provides Efficiency: %s\n",
	       fName.Data(), fEfficiency->GetTitle());
    }
    
    if ((obj = fResponses->FindObject("Resolution"))) 
    {
	fResolution = (AliFastResponse*) obj;
	printf("Detector %s provides Resolution: %s\n",
	       fName.Data(), fResolution->GetTitle());
    }
}

Float_t AliFastDetector::EvaluateEfficiency(AliFastParticle* part)
{
//
//  Evaluate the efficiency for detecting particle part
//
    TIter nextDet(fSubdetectors);
    AliFastDetector *det;
    //
    // Loop over subdetectors  
    Float_t eff = 1;
    while((det = (AliFastDetector*)nextDet())) {
	eff *= det->EvaluateEfficiency(part);
    }  
    return eff;
}

Bool_t  AliFastDetector::EvaluateAcceptance(AliFastParticle* part)
{
    //
    // Loop over subdetectors 
    Bool_t acc = kFALSE;

    if (fSubdetectors) {
	TIter nextDet(fSubdetectors);
	AliFastDetector *det;
	while((det = (AliFastDetector*)nextDet()) && !acc) {
	    acc = (acc ||  det->EvaluateAcceptance(part));
	}  
    } else {
        if (fGeometry)
	acc = fGeometry->Impact((TParticle*) part);
    }
    
    return acc;
}

void    AliFastDetector::EvaluateResponse(AliFastParticle* /*part*/)
{
    ;
}

void AliFastDetector::
AddSubdetector(AliFastDetector *Detector, char* /*Name*/)
{
//
//  Add detector to list   
     fSubdetectors->Add(Detector);
}



void AliFastDetector::
AddResponse(AliFastResponse *Response)
{
//
//  Add detector to list   
     fResponses->Add(Response);
}


AliFastDetector*  AliFastDetector::FirstSubdetector()
{
// Iterator over generators: Initialisation
    fLnkD = fSubdetectors->FirstLink();
    if (fLnkD) {
	return (AliFastDetector*) (fLnkD->GetObject());
    } else {
	return 0;
    }
}

AliFastDetector*  AliFastDetector::NextSubdetector()
{
// Iterator over generators: Increment
    fLnkD = fLnkD->Next();
    if (fLnkD) {
	return (AliFastDetector*) (fLnkD->GetObject());
    } else {
	return 0;
    }
}


AliFastResponse*  AliFastDetector::FirstResponse()
{
// Iterator over generators: Initialisation
    fLnkR = fResponses->FirstLink();
    if (fLnkR) {
	return (AliFastResponse*) (fLnkR->GetObject());
    } else {
	return 0;
    }
}

AliFastResponse*  AliFastDetector::NextResponse()
{
// Iterator over generators: Increment
    fLnkR = fLnkR->Next();
    if (fLnkR) {
	return (AliFastResponse*) (fLnkR->GetObject());
    } else {
	return 0;
    }
}


AliFastDetector& AliFastDetector::operator=(const  AliFastDetector& rhs)
{
// Assignment operator
    rhs.Copy(*this);
    return *this;
}

void AliFastDetector::Copy(TObject&) const
{
    //
    // Copy 
    //
    Fatal("Copy","Not implemented!\n");
}

