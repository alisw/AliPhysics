#ifndef ALIGAMMACONVERSIONAODOBJECT_H
#define ALIGAMMACONVERSIONAODOBJECT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

////////////////////////////////////////////////
//--------------------------------------------- 
// Class containing the aod information
//---------------------------------------------
////////////////////////////////////////////////

// --- ROOT system ---
#include "TObject.h" 
class AliStack;
class AliESDEvent;

class AliGammaConversionAODObject : public TObject {

 public: 

  AliGammaConversionAODObject();                                        //constructor
  AliGammaConversionAODObject(const AliGammaConversionAODObject & g);                   //copy constructor
  AliGammaConversionAODObject & operator = (const AliGammaConversionAODObject & g);     //assignment operator
  virtual ~AliGammaConversionAODObject() {;}                            //virtual destructor

  /*
   * This function sets the Px
   */
  void SetPx(Float_t px){fPx = px;}

  /*
   * This function sets the Py
   */
  void SetPy(Float_t py){fPy = py;}

  /*
   * This function sets the Pz
   */
  void SetPz(Float_t pz){fPz = pz;}

  /*
   * This function sets the esd label of the first electron
   */
  void SetLabel1(Int_t label){fLabel1 = label;}

  /*
   * This function sets the esd label of the second electron
   */
  void SetLabel2(Int_t label){fLabel2 = label;}
  
  /*
   * This function returns the Px
   */
  Float_t GetGammaPx() const{return fPx;}

  /*
   * This function returns the Py
   */
  Float_t GetGammaPy() const{return fPy;}

  /*
   * This function returns the Pz
   */
  Float_t GetGammaPz() const{return fPz;}

  /*
   * This function returns the esd label of the first electron
   */
  Int_t GetElectronLabel1() const{return fLabel1;}

  /*
   * This function returns the esd label of the second electron
   */
  Int_t GetElectronLabel2()const {return fLabel2;}


  /*
   * This function sets the MC stack
   */
  void SetStack(AliStack* stack){fMCStack=stack;}

  /*
   * This function sets the ESD event
   */
  void SetESDEvent(AliESDEvent* esdEvent){fESDEvent = esdEvent;}

  /*
   * This function returns the Gamma MC label
   */
  Int_t GetGammaMCLabel() const;

  /*
   * This function returns the unique id  of the electrons (if they have the same mother and unique id)
   */
  Int_t GetElectronUniqueID() const;

  /*
   * This function returns the MC label of the first electron
   */
  Int_t GetElectronMCLabel1() const;

  /*
   * This function returns the MC label of the second electron
   */
  Int_t GetElectronMCLabel2() const;

 private:

  Float_t fPx;
  Float_t fPy;
  Float_t fPz;
  Int_t fLabel1;
  Int_t fLabel2;
  AliStack* fMCStack;
  AliESDEvent * fESDEvent;

  ClassDef(AliGammaConversionAODObject,0)
};


#endif



