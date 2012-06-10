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
#include <TArrayI.h>
//______________________________________________________________________
class AliITSDigitUpgrade: public AliDigit {

 public:
  AliITSDigitUpgrade();
  AliITSDigitUpgrade(Int_t *digits);
  AliITSDigitUpgrade(ULong_t pixid, Float_t eloss);
  AliITSDigitUpgrade(const AliITSDigitUpgrade &d); //copy ctor

  virtual ~AliITSDigitUpgrade() {/*destructor*/};
  //____________________________________________________________________________________________________
    
  enum {kMaxLab=12}; // maximum number of MC labels associated to the digit (4 times as much as can be stored in the "mother class")

  void AddTidOffset(Int_t offset) {for(Int_t i=0; i<kMaxLab; i++) if (fTrackIdMC[i]>-1) fTrackIdMC[i]+=offset;} //needed for merging
    
  // setters
       
  void SetSignal(Float_t sig) {fSignal = sig;}
  void SetLayer(Int_t layer) {fNLayer = layer;}
  void SetModule(Int_t module) {fModule = module ;}
  void SetNTracksIdMC(Int_t nLabels) {fNTracksIdMC = nLabels;}
  void SetNelectrons(Double_t nele) {fNelectrons = nele;}


  void SetTrackID(Int_t tid) {fTrackIdMC[0]=tid; }
  void SetTids(Int_t tids[kMaxLab]){ for(Int_t i=0; i<kMaxLab; i++) fTrackIdMC[i]=tids[i];} // tracks participating to form the digit
  void AddTrackID(Int_t tid); 
  void SetPixId(ULong_t nx, ULong_t nz) {fPixId = 100000*nx + nz;}
  void SetPixId(ULong_t pixid){fPixId = pixid;}

  void SetSignalID(Float_t eloss[kMaxLab]){for(Int_t i=0; i< kMaxLab; i++) fSignalID[i]=eloss[i];}

  // getters
    
  Float_t  GetSignal() const {return fSignal;}
  Int_t    GetLayer() const {return fNLayer;}
  Int_t    GetModule() const {return fModule;}
  Double_t GetNelectrons() const {return fNelectrons;}
  ULong_t  GetPixId(){return fPixId;}
  Int_t    GetxPixelNumber() const {return fPixId/100000;}
  Int_t    GetzPixelNumber() const {return fPixId%100000;}
  Int_t    GetNTracksIdMC() const {return fNTracksIdMC;}
  Int_t*   GetTracks() {return fTrackIdMC; }
  Int_t    GetTrackID(Int_t ipart) const  {if(ipart<0 || ipart>=kMaxLab) return -1; else return fTrackIdMC[ipart];}
  Float_t  GetSignalID(Int_t ipart) const {if(ipart<0 || ipart>=kMaxLab) return -1; else return fSignalID[ipart];}
    
  void GetPosition(Int_t ilayer, Int_t nx, Int_t nz, Double_t &xloc, Double_t &zloc);
    
  void PrintInfo(); 
  inline Int_t   Compare     (const TObject *pObj)const;
  Bool_t IsSortable() const {return kTRUE;}
   
 protected:
    
  ULong_t fPixId;    // ID number of the fired pixel in a module
  Float_t fSignal;   // Signal as Eloss in the medium
  Int_t fNLayer;     // Layer 
  Int_t fModule;     // Module in the layer
  Double_t fNelectrons; // released charge due to E loss
  Int_t fNTracksIdMC;  // Number of MC particles which produced the Digit
  Int_t fTrackIdMC[kMaxLab];  // MC track labels 
  Float_t fSignalID[kMaxLab]; // E loss
  

 private:
  AliITSDigitUpgrade &operator=(const AliITSDigitUpgrade &);    // not implemented


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


