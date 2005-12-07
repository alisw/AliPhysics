#ifndef ALIMUON_H
#define ALIMUON_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
// Revision of includes 07/05/2004

/// \ingroup base
/// \class AliMUON
/// \brief AliDetector Class for MUON subsystem

////////////////////////////////////////////////
//  AliDetector Class for MUON subsystem      //
////////////////////////////////////////////////

#include <TObjArray.h>

#include "AliDetector.h"
#include "AliMUONData.h"
#include "AliMUONChamber.h"

class TVector;
class TFile;
class TTree;

class AliLoader;
class AliMUONTriggerCircuit;
class AliMUONData;
class AliMUONResponse;
class AliMUONSegmentation;
class AliMUONHit;
class AliMUONRawCluster;
class AliMUONGeometry;
class AliMUONGeometryTransformer;
class AliMUONGeometryBuilder;
class AliMUONVGeometryBuilder;
class AliESD;

class AliMUON : public  AliDetector 
{
  public:
    AliMUON();
    AliMUON(const char *name, const char *title);
    virtual       ~AliMUON();
   
    // Geometry
    void           AddGeometryBuilder(AliMUONVGeometryBuilder* geomBuilder);
    virtual void   BuildGeometry();
    virtual Int_t  IsVersion()   const {return 0;}
    
    AliMUONGeometryBuilder*            GetGeometryBuilder() const {return fGeometryBuilder;}
    const AliMUONGeometryTransformer*  GetGeometryTransformer() const;
    AliMUONSegmentation*               GetSegmentation() const    { return fSegmentation; }

    // MUONData   
    AliMUONData*   GetMUONData() {return fMUONData;}

    // MUONLoader definition
    virtual AliLoader* MakeLoader(const char* topfoldername); //builds standard getter (AliLoader type)
    // Interface with AliMUONData
    virtual void       MakeBranch(Option_t *opt=" ") {GetMUONData()->MakeBranch(opt);}
    virtual void       SetTreeAddress();
    virtual void       ResetHits()                   {GetMUONData()->ResetHits();}
    virtual void       ResetDigits()                 {GetMUONData()->ResetDigits();}
    virtual void       ResetTrigger()                {GetMUONData()->ResetTrigger();}
    virtual void       ResetRawClusters()            {GetMUONData()->ResetRawClusters();}
    virtual void       SetSplitLevel(Int_t SplitLevel)     {fSplitLevel=SplitLevel;}

    // Digitisation 
    virtual AliDigitizer* CreateDigitizer(AliRunDigitizer* manager) const;
    virtual void   SDigits2Digits();      
    virtual void   Hits2SDigits();
    virtual void   Digits2Raw();


    // Configuration Methods (per station id)
    //
    // Set Signal Generation Parameters
    virtual void   SetSigmaIntegration(Int_t id, Float_t p1);
    virtual void   SetChargeSlope(Int_t id, Float_t p1);
    virtual void   SetChargeSpread(Int_t id, Float_t p1, Float_t p2);
    virtual void   SetMaxAdc(Int_t id, Int_t p1);
    // Set Response Model
    virtual void   SetResponseModel(Int_t id, AliMUONResponse *response);

    // Set Stepping Parameters
    virtual void   SetMaxStepGas(Float_t p1);
    virtual void   SetMaxStepAlu(Float_t p1);
    virtual void   SetMaxDestepGas(Float_t p1);
    virtual void   SetMaxDestepAlu(Float_t p1);
   
    // Get Stepping Parameters
    virtual Float_t  GetMaxStepGas() const;
    virtual Float_t  GetMaxStepAlu() const;
    virtual Float_t  GetMaxDestepGas() const;
    virtual Float_t  GetMaxDestepAlu() const;
    
    // Set alignement option
    virtual void  SetAlign(Bool_t align = true);
    virtual void  SetAlign(const TString& fileName, Bool_t align = true);
   
    // Return reference to Chamber #id
    virtual AliMUONChamber& Chamber(Int_t id)
      {return *((AliMUONChamber *) (*fChambers)[id]);}
    // Return reference to Circuit #id
    virtual AliMUONTriggerCircuit& TriggerCircuit(Int_t id)
      {return *((AliMUONTriggerCircuit *) (*fTriggerCircuits)[id]);}
    // Return pointers to digits
    AliMUONRawCluster    *RawCluster(Int_t ichamber, Int_t icathod,
				     Int_t icluster);
    // Inherited and overridden from AliModule:
    //PH    virtual void RemapTrackHitIDs(Int_t * map);

  protected:
    AliMUON(const AliMUON& rMUON);
    AliMUON& operator = (const AliMUON& rhs);

    const AliMUONGeometry* GetGeometry() const;

    Int_t                 fNCh;                // Number of chambers   
    Int_t                 fNTrackingCh;        // Number of tracking chambers*
    AliMUONData*          fMUONData;           // Data container for MUON subsystem  
    Int_t                 fSplitLevel;         // Splitlevel when making branches in outfiles.
    TObjArray*            fChambers;           // List of Tracking Chambers
    TObjArray*            fTriggerCircuits;    // List of Trigger Circuits
    AliMUONGeometryBuilder*  fGeometryBuilder; // Geometry builder 
    AliMUONSegmentation*  fSegmentation;       // New segmentation 
   
    //
    Bool_t   fAccCut;          //Transport acceptance cut
    Float_t  fAccMin;          //Minimum acceptance cut used during transport
    Float_t  fAccMax;          //Minimum acceptance cut used during transport
    //  
    //  Stepping Parameters
    Float_t fMaxStepGas;      // Maximum step size inside the chamber gas
    Float_t fMaxStepAlu;      // Maximum step size inside the chamber aluminum
    Float_t fMaxDestepGas;    // Maximum relative energy loss in gas
    Float_t fMaxDestepAlu;    // Maximum relative energy loss in aluminum
    
    // Pad Iterator
    Int_t fMaxIterPad;        // Maximum pad index
    Int_t fCurIterPad;        // Current pad index
    // Background eent for event mixing
    
    ClassDef(AliMUON,9)  // MUON Detector base class
};
#endif

