#ifndef ALIITSDigitUpgrade_H
#define ALIITSDigitUpgrade_H
/* Copyright(c) 2004-2006, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
/////////////////////////////////////////////////////////////
//              Authors A.Mastroserio			   //		
//			C.Terrevoli			   // 
//			annalisa.mastroserio@cern.ch	   //
//			cristina.terrevoli@cern.ch         //
// 							   //
// 		Digit class for ITS upgrade                //
//							   //
/////////////////////////////////////////////////////////////

/* $Id$ */

#include <AliDigit.h>

//______________________________________________________________________
class AliITSDigitUpgrade: public AliDigit {

 public:
  AliITSDigitUpgrade();
  AliITSDigitUpgrade(Int_t *digits);
  AliITSDigitUpgrade(ULong_t pixid, Float_t eloss);
  AliITSDigitUpgrade(const AliITSDigitUpgrade &d); //copy ctor

  virtual ~AliITSDigitUpgrade(){/*destructor*/;}
  //____________________________________________________________________________________________________
    

  void AddTidOffset(Int_t offset) {for(Int_t i=0; i<3; i++) if (fTracks[i]>0) fTracks[i]+=offset;} //needed for merging
    
  // setters
       
  void SetSignal(Float_t sig) {fSignal = sig;}
  void SetLayer(Int_t layer) {fNLayer = layer;}
  void SetModule(Int_t module) {fModule = module ;}
  void SetNelectrons(Double_t nele) {fNelectrons = nele;}
  void SetTrackID(Int_t tid) {fTracks[0]=tid;}
  void SetPixId(ULong_t nx, ULong_t nz) {fPixId = 100000*nx + nz ;}
  void SetPixId(ULong_t pixid){fPixId = pixid;}
  void SetTids(Int_t tids[3]){ for(Int_t i=0; i<3; i++) fTracks[i]=tids[i];} // tracks participating to form the digit
  void SetSignalID(Float_t eloss[3]){for(Int_t i=0; i< 3; i++) fSignalID[i]=eloss[i];}

  // getters
    
  Float_t  GetSignal() const {return fSignal;}
  Int_t    GetLayer() const {return fNLayer;}
  Int_t    GetModule() const {return fModule;}
  Double_t GetNelectrons() const {return fNelectrons;}
  ULong_t  GetPixId(){return fPixId;}
  Int_t    GetxPixelNumber() const {return fPixId/100000;}
  Int_t    GetzPixelNumber() const {return fPixId%100000;}
    
  void GetPosition(Int_t ilayer, Int_t nx, Int_t nz, Double_t &xloc, Double_t &zloc);
  Float_t GetSignalID(Int_t ipart) {if(ipart<0 || ipart >2) return -1; else return fSignalID[ipart];}
    
  void PrintInfo(); 
  inline Int_t   Compare     (const TObject *pObj)const;
  Bool_t IsSortable() const {return kTRUE;}
   
 protected:
    
  ULong_t fPixId;
  Float_t fSignal;   // Signal as Eloss in the medium
  Int_t fNLayer;     
  Int_t fModule;
  Double_t fNelectrons; 
  Float_t fSignalID[3];

    
  ClassDef(AliITSDigitUpgrade,3)   // Simulated digit object for ITS upgrade

    };
#endif
Int_t AliITSDigitUpgrade::Compare(const TObject *pObj) const
{
  // Arguments: pObj - pointer to object to compare with
  //        

  Int_t result = -1;
   if (fModule>((AliITSDigitUpgrade*)pObj)->GetModule()) result=1;      

  else  if(fModule==((AliITSDigitUpgrade*)pObj)->GetModule()){
   if     (fPixId==((AliITSDigitUpgrade*)pObj)->GetPixId()) result=0;
   else if(fPixId >((AliITSDigitUpgrade*)pObj)->GetPixId()) result=1;
   }
  return result;
}


