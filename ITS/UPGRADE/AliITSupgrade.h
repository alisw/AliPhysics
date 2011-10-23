#ifndef ALIITSUPGRADE_H
#define ALIITSUPGRADE_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//========================================================================
//
//      Geometry of the Inner Tracking System Upgrade based on TGeo
//
//
//   Authors: A. Mastroserio
//            C. Terrevoli
//            annalisa.mastroserio@cern.ch
//	      cristina.terrevoli@ba.infn.it
//
//========================================================================

//
//ITS class for new geometry based on TGeo
//

#include "AliITS.h"             //base class 
#include "AliITSupgradeDigitizer.h"
#include <TArrayS.h>
    
class TGeoVolume;
class TGeoVolumeAssembly;
class TArrayD;
class TClonesArray;
class TObjArray;
class AiITShit;
class AliITSsegmentationUpgrade;

class AliITSupgrade : public AliITS //TObject-TNamed-AliModule-AliDetector-AliITS-AliITSupgrade
{
 public:
  AliITSupgrade();                                                              //default ctor
  AliITSupgrade(const char *name, const char *title, Bool_t isBeamPipe=kFALSE);  //ctor for standard ITS 
  AliITSupgrade(const char *name, const char *title, TArrayD widths, TArrayD radii,TArrayD halfLengths, TArrayD radiiCu, TArrayD widthsCu, TArrayS copper,Bool_t bp, Double_t radiusBP, Double_t widthPB, Double_t halfLengthsBP);  //ctor
  virtual       ~AliITSupgrade();                                       //dtor
  Int_t IsVersion() const {return 110;}// same as ITSv11Hybrid to comply with the AliITSInitGeom cases
 
  
  //framework part+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
  void    AddAlignableVolumes() const;             //from AliModule invoked from AliMC           
  void    CreateMaterials  ();                     //from AliModule invoked from AliMC
  void    CreateGeometry   ();                     //from AliModule invoked from AliMC                     
  void    StepManager      ();                     //from AliModule invoked from AliMC::Stepping()
  void    Hits2SDigits();
  void    MakeBranch      (Option_t *opt="");        //from AliModule invokde from AliRun::Tree2Tree() to make requested ITSupgrade branch
  void    SetTreeAddress  (                );        //from AliModule invoked from AliRun::GetEvent(), AliLoader::SetTAddrInDet()
          
  //create containers++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  void     HitCreate()  	   {if(fHits)return; fHits=new TClonesArray("AliITShit"); fNhits=0;     }//create hits list 

  TObjArray*      SDigitsList()            const{return fSdigits;} //get digits list for all layers
  TClonesArray*   SDigitsList(Int_t layer) const{return fSdigits ? (TClonesArray *)fSdigits->At(layer):0;}//get sdigits list for a layer
  void            SDigitsCreate()            {if (fSdigits) return; 
    fSdigits=new  TObjArray(fNlayers);
    for(Int_t i=0;i<fNlayers;i++)fSdigits->AddAt(new TClonesArray("AliITSDigitUpgrade"),i);}//create digits list
  void           SDigitsReset ()              {if(fSdigits) for(int i=0;i<fNlayers;i++)fSdigits->At(i)->Clear();                     }//clean sdigits list  
     
     
     
  TObjArray*      DigitsList()            const{return fDigitArray;}//get digits list for all layers
  TClonesArray*   DigitsList(Int_t layer) const{return fDigitArray ? (TClonesArray *)fDigitArray->At(layer):0;}//get digits list for chamber
  void            DigitsCreate()            {if (fDigitArray) return; 
    fDigitArray=new TObjArray(fNlayers);
    for(Int_t i=0;i<fNlayers;i++)fDigitArray->AddAt(new TClonesArray("AliITSDigitUpgrade"),i);}//create digits list
  void          DigitsReset ()              {if(fDigitArray)for(int i=0;i<fNlayers;i++)fDigitArray->At(i)->Clear();                     }//clean digits list

  AliDigitizer*   CreateDigitizer  (AliDigitizationInput *m) const {return new AliITSupgradeDigitizer(m);}  //from AliModule invoked from AliSimulation::RunDigitization()
           
  void Hit2SumDig(TClonesArray *hits,const TObjArray *pSDig, Int_t *nSdigit);

  enum EMedia {kAir=0,kSi=1,kBe=2, kCu=3};                         //media ids used in CreateMaterials  
    

  // Geometry/segmentation creation part
  TGeoVolumeAssembly * CreateVol();
  TGeoVolumeAssembly * CreateVol(Bool_t cil);
  
  void SetFullSegmentation(TArrayD xsize, TArrayD zsize);

  void SetRadius(Double_t r, Int_t lay) {fRadii.AddAt(r,lay);}
  void SetWidth(Double_t w, Int_t lay) {fWidths.AddAt(w,lay);}
  void SetRadiusCu(Double_t rCu, Int_t lay) {fRadiiCu.AddAt(rCu,lay);}
  void SetWidthCu(Double_t wCu, Int_t lay) {fWidthsCu.AddAt(wCu,lay);}
  void SetSegmentationX(Double_t x, Int_t lay);
  void SetSegmentationZ(Double_t z, Int_t lay);
 void SetNSectors(Int_t nsect) {fNSectors=nsect;}
  void StepHistory();
  void PrintSummary();

 protected:
  TArrayD fWidths;   			      // thicknesses of silicon cylinders
  TArrayD fRadii;    			      // radii of silicon cylinders
  Int_t fNSectors;   			      // number of sectors
  TArrayD fRadiiCu;  			      // radii of copper cylinders 
  TArrayD fWidthsCu; 			      // thicknesses of copper cylinders
  TArrayS   fCopper; 			      // check to set the copper cylinder 
  Bool_t fBeampipe;                           // boolean to set beam pipe vol in geometry
  Double_t fRadiusBP;                         // beam pipe radius
  Double_t fWidthBP; 			      // thickness of the beam pipe		
  Double_t fHalfLengthBP;		      // 1/2 beam pipe cyclinder length
  Int_t   fNlayers;			      // # of layers
  TArrayD fHalfLength;			      // 1/2 silicon and copper cylinder lengths
  TObjArray            *fSdigits;             //![fNlayers] list of sdigits
  TObjArray            *fDigitArray;          //![fNlayers] list of digits per layer
  AliITSsegmentationUpgrade *fSegmentation;   //virtual segmentation class
  
 private:
  AliITSupgrade(const AliITSupgrade& u);              //dummy copy constructor
  AliITSupgrade &operator=(const AliITSupgrade& u);   //dummy assignment operator

  ClassDef(AliITSupgrade,2)                                                                //ITS Upgrade full version for simulation
  
    };

#endif

