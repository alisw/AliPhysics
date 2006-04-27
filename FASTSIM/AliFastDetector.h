#ifndef ALIFASTDETECTOR_H
#define ALIFASTDETECTOR_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//
// Base class for fast simulation of a detctor
// or a system of subdetectors.
// The detector response is described by resolution and efficiency.
// Author:
// Andreas Morsch
// andreas.morsch@cern.ch

#include <TNamed.h>
class TList;
class TObjLink;
class AliFastResponse;
class AliFastParticle;
class AliGeometry;

class AliFastDetector : public TNamed {
    
 public:
    AliFastDetector();
    AliFastDetector(char* Name, char* Title);
    AliFastDetector(const AliFastDetector& det);    
    virtual ~AliFastDetector();
    virtual void Init();
    virtual void SetGeometry(AliGeometry* geom) 
	{fGeometry = geom;}
    
    virtual AliGeometry* GetGeometry() const  
	{return fGeometry;}
    //
    // Add a new subdetector 
    virtual void AddSubdetector(AliFastDetector *Detector, char* Name);
    virtual TList* Subdetectors() {return fSubdetectors;}
    //
    // Add a new response
    virtual void AddResponse(AliFastResponse *Response);
    virtual TList* Responses() {return fResponses;}
    virtual Float_t EvaluateEfficiency(AliFastParticle* part);
    virtual Bool_t  EvaluateAcceptance(AliFastParticle* part);
    virtual void    EvaluateResponse(AliFastParticle* part);
    
    // Iterators
    AliFastDetector*  FirstSubdetector();
    AliFastDetector*  NextSubdetector();
    AliFastResponse*  FirstResponse();
    AliFastResponse*  NextResponse();
    // Copy
    AliFastDetector& operator=(const AliFastDetector & rhs);
    void Copy(TObject&) const;
 protected:
    TList            *fSubdetectors;      // List of Subdetectors
    TList            *fResponses;         // Responses
    TObjLink         *fLnkD;              // Pointer to detector in list 
    TObjLink         *fLnkR;              // Pointer to response in list
    AliFastResponse  *fEfficiency;        // Efficiency Simulation
    AliFastResponse  *fResolution;        // Resolution Simulation
    AliGeometry      *fGeometry;          // Geometry 
    ClassDef(AliFastDetector,1) // Base class for fast detector
};

#endif 



