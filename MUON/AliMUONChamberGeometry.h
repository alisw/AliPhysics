/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
// Revision of includes 07/05/2004

//
// Class AliMUONChamberGeometry
// -----------------------------
// Class for definititon of the MUON chamber positions in ALIC.
//
// Author: Ivana Hrivnacova, IPN Orsay

#ifndef ALI_MUON_CHAMBER_GEOMETRY_H
#define ALI_MUON_CHAMBER_GEOMETRY_H

#include <TObject.h>

class TGeoTranslation;
class TGeoRotation;
class TGeoCombiTrans;
class TObjArray;
class TArrayI;

class AliMUONChamber;
class AliMUONGeometryEnvelope;

class AliMUONChamberGeometry : public TObject
{
  public:
    AliMUONChamberGeometry(Int_t chamberId);
    AliMUONChamberGeometry();
    virtual ~AliMUONChamberGeometry();

    // methods
        
          // adding virtual envelopes 	
          // (not placed in MC geometry, only logical assembly of volumes,
	  //  cannot have more copies)	
    void  AddEnvelope(const TString& name, Bool_t isVirtual, 
                      const char* only="ONLY"); 
    void  AddEnvelope(const TString& name, Bool_t isVirtual, 
                      const TGeoTranslation& translation, 
		      const char* only="ONLY"); 
    void  AddEnvelope(const TString& name, Bool_t isVirtual, 
                      const TGeoTranslation& translation, 
		      const TGeoRotation& rotation,
		      const char* only="ONLY");
		      
          // adding non-virtual envelopes 	
          // (placed in MC geometry with transformation composed
	  //  of transformation of chamber and their transformation, 
	  //  can have more copies )	
    void  AddEnvelope(const TString& name, Int_t copyNo, 
                      const char* only="ONLY"); 
    void  AddEnvelope(const TString& name, Int_t copyNo, 
                      const TGeoTranslation& translation,
		      const char* only="ONLY"); 
    void  AddEnvelope(const TString& name, Int_t copyNo, 
                      const TGeoTranslation& translation, 
		      const TGeoRotation& rotation,
		      const char* only="ONLY");

          // adding constituents to virtual envelopes 	
          // (placed in MC geometry with transformation composed
	  //  of transformation of chamber, envelope and their own
	  //  transformation )	
    void  AddEnvelopeConstituent(const TString& name, const TString& envName, 
                      Int_t copyNo); 
    void  AddEnvelopeConstituent(const TString& name, const TString& envName, 
                      Int_t copyNo, const TGeoTranslation& translation); 
    void  AddEnvelopeConstituent(const TString& name, const TString& envName, 
                      Int_t copyNo, const TGeoTranslation& translation, 
		      const TGeoRotation& rotation);
		      		      
          // adding constituents to virtual envelopes with specified shape
	  // parameters
          // (placed in MC geometry with transformation composed
	  //  of transformation of chamber, envelope and their own
	  //  transformation )	
    void  AddEnvelopeConstituentParam(const TString& name, const TString& envName, 
                      Int_t copyNo, Int_t npar, Double_t* param); 
    void  AddEnvelopeConstituentParam(const TString& name, const TString& envName, 
                      Int_t copyNo, const TGeoTranslation& translation,
		      Int_t npar, Double_t* param); 
    void  AddEnvelopeConstituentParam(const TString& name, const TString& envName, 
                      Int_t copyNo, const TGeoTranslation& translation, 
		      const TGeoRotation& rotation, Int_t npar, Double_t* param);
		      		      
    void  SetMotherVolume(const TString& motherVolumeName);
    void  SetTranslation(const TGeoTranslation& translation);
    void  SetRotation(const TGeoRotation& rotation);
    
    void  SetSensitiveVolume(Int_t volId);
    void  SetSensitiveVolume(const TString& name);
    void  SetDebug(Bool_t debug);

    // get methods
    TString                GetMotherVolume() const;
    const TGeoCombiTrans*  GetTransformation() const;
    const TObjArray*       GetEnvelopes() const;
    Bool_t IsSensitiveVolume(Int_t volId) const; 

  protected:
    AliMUONChamberGeometry(const AliMUONChamberGeometry& rhs);
    // operators  
    AliMUONChamberGeometry& operator = (const AliMUONChamberGeometry& rhs);

  private:
    // methods
    AliMUONGeometryEnvelope* FindEnvelope(const TString& name) const;
  
    // data members
    Int_t            fChamberId;     // the chamber Id
    TString          fMotherVolume;  // mother volume name
    TGeoCombiTrans*  fTransformation;// the chamber transformation wrt to mother
                                     // volume
    TObjArray*       fEnvelopes;     // the envelopes names and transformations
		                     // wrt to the chamber position in mother volume                                 
    Int_t            fNofSensVolumeIds; // Number of sensitive volumes IDs    
    TArrayI*         fSensVolumeIds; // Sensitive volumes IDs  
    Bool_t           fDebug;         // Switch for debugging  
 
  ClassDef(AliMUONChamberGeometry,1) // MUON chamber geometry base class
};

// inline functions

inline void  AliMUONChamberGeometry::SetMotherVolume(const TString& motherVolumeName)
{ fMotherVolume = motherVolumeName; }

inline void  AliMUONChamberGeometry::SetDebug(Bool_t debug)
{ fDebug = debug; }

inline TString  AliMUONChamberGeometry::GetMotherVolume() const
{ return fMotherVolume; }

inline const TGeoCombiTrans* AliMUONChamberGeometry::GetTransformation() const 
{ return fTransformation; }

inline const TObjArray* AliMUONChamberGeometry::GetEnvelopes() const
{ return fEnvelopes; }

#endif //ALI_MUON_V_CHAMBER_GEOMETRY_H
