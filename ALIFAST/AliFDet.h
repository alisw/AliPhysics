#ifndef AliFDet_H
#define AliFDet_H
////////////////////////////////////////////////////////////////////////////
//                                                                        //
// AliFast Detector Class                                                 //
//                                                                        //
// to provide information of effective material (X/Xo) of the detector    //
// needed for the multiple scattering formula used in AliFTrackMaker.     // 
//                                                                        // 
// the number and dimensions of cylindrical layers of material are        //
// initialised here for the TP status and are to be updated accordingly.  //
//                                                                        //
// this class is replacing the "init_geometry" routine of program "res.f  //
//                                                                        // 
////////////////////////////////////////////////////////////////////////////

//#ifndef ROOT_TObject
#include <TNamed.h>
//#endif

enum { kNMaxDet  = 100 };
enum { kNMaxDet2 = 200 };


class AliFDet : public TNamed {
  
private:
   //geometry parameters
   Double_t        fRDet[kNMaxDet];         // radius of detector material in cm
   Double_t        fRDetSQ[kNMaxDet];       // and the sq root of it
   Double_t        fThickDet[kNMaxDet];     // thickness divided by Xo
   //errors due to detector precision plus alignement given by groups. 
   //they are momentum dependent; for TPC are calculated properly.
   Double_t        fErrorRPhi[kNMaxDet];    // error in bending direction
   Double_t        fErrorZ[kNMaxDet];       // error in z direction
   Double_t        fErrorR[kNMaxDet];       // error in r direction,from alignement only 
   Int_t           fIFlagDet[kNMaxDet];     // 1: sensitive detector 
                                            // 2: errors will be calculated
   Int_t           fIFlagGas[kNMaxDet];     // for gas detectors  
   //vertex precision calculated for particle multiplicity mult_density
   //high multiplicity results in optimistic errors for vertex   
   Double_t        fErrorVertexX;           // vertex precision in x
   Double_t        fErrorVertexY;           // vertex precision in y
   Double_t        fErrorVertexZ;           // vertex precision in z
   Double_t        fBMag;                   // magnetic field in KGauss
   Double_t        fConstMag;               //
   Int_t           fNDetActive;             // n. of active detector layers
   Int_t           fNDet;                   // n. of detectors layers

public:
                  AliFDet() {}
                  AliFDet(const char *name, const char *title);
   virtual       ~AliFDet() {}

   // Initialise parameters for detector geometry
   void           InitDetParam();
   void           PrintDetInfo();


   // Getters
   Double_t        RDet(Int_t idDet) {return fRDet[idDet];}         
   Double_t        RDetSQ(Int_t idDet) {return fRDetSQ[idDet];}     
   Double_t        ThickDet(Int_t idDet) {return fThickDet[idDet];}    
   Double_t        ErrorRPhi(Int_t idDet) {return fErrorRPhi[idDet];}  
   Double_t        ErrorZ(Int_t idDet) {return fErrorZ[idDet];}     
   Double_t        ErrorR(Int_t idDet) {return fErrorR[idDet];} 
   Int_t           IFlagDet(Int_t idDet) {return fIFlagDet[idDet];}   
   Int_t           IFlagGas(Int_t idDet) {return fIFlagGas[idDet];}   
   Double_t        ErrorVertexX() {return fErrorVertexX;} 
   Double_t        ErrorVertexY() {return fErrorVertexY;} 
   Double_t        ErrorVertexZ() {return fErrorVertexZ;} 
   Double_t        BMag() {return fBMag;}                  
   Double_t        ConstMag() {return fConstMag;}         
   Int_t           NDetActive() {return fNDetActive;}      
   Int_t           NDet() {return fNDet;}                 


   ClassDef(AliFDet,1)   //AliFast Detector intialisation for AliFTrackMaker
};

#endif










