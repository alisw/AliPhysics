#ifndef ALIITSMULTRECBG_H
#define ALIITSMULTRECBG_H

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Class for generating combinatorial backgroung                             //
// for the SPD tracklets                                                     //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

class TH2F;
class TTree;
#include <TArrayF.h>
#include "../ITS/AliITSMultReconstructor.h"
#include "../ITS/AliITSRecPoint.h"


class AliITSMultRecBg : public AliITSMultReconstructor 
{
 public:
  enum {kData,  // normal reconstruction
	kBgRot, // bg with rotation
	kBgInj, // bg with injection
	kBgMix  // bg with mixing
  };
  //
  AliITSMultRecBg();
  virtual ~AliITSMultRecBg();
  //
  void    Run(TTree* tree, Float_t* vtx, TTree* treeMix=0);
  //
  void    SetRecType(UInt_t m=kData)                                  {fRecType = m;}
  Int_t   GetRecType()                                          const {return fRecType;}
  //
  void    SetInjScale(Float_t s=1.)                  {fInjScale = s>0? s:1.;}
  Int_t   GetNInjTrials(Int_t lr=-1)           const {return lr<0 ? fInjNTrials[0]+fInjNTrials[1] : fInjNTrials[lr];}
  Int_t   GetNInjSuccsess(Int_t lr=-1)         const {return lr<0 ? fInjNSuccess[0]+fInjNSuccess[1] : fInjNSuccess[lr];}
  //
  void    Print(Option_t *opt=0)               const {AliITSMultReconstructor::Print(opt);}
  //
 protected:
  virtual void CreateMultiplicityObject();
  //
  AliITSMultRecBg(const AliITSMultRecBg &src);
  //     Injection stuff
  void    GenInjBgSample(TTree* treeRP, Float_t *vtx);
  Bool_t  PrepareInjBgGenerator(Float_t *vtx);
  void    InitInjBg();
  Bool_t  GenClusterToInject();
  void    PlaceInjCluster();
  Int_t   PickClusterPrototype(Int_t lr, Int_t ladInStave, Int_t stave2avoid);
  UInt_t* GenClusterPattern(Int_t &npix, Int_t &ny, Int_t &nz, Float_t cy,Float_t &cz);
  Int_t   SearchInjTracklet(const Float_t *vtx);
  Int_t   GetPixelBitC(int stave, int chip,  int col, int row) const {return row+((col+((chip+stave*20)<<5))<<8);} // row + (col+ (chip + stave*20)*32)*256;
  Int_t   GetPixelBitL(int stave, int ladder,int col, int row) const {return row+((col+(ladder+(stave<<2))*32*5)<<8);} // row + (col + (ladder + stave*4)*32*5)*256
  void    ExtractPixelL(int id, int &stave, int &ladder, int &col, int &row) const;


 private: 
  AliITSMultRecBg& operator=(const AliITSMultRecBg& src);

 protected:
  enum {kInjFakeLabel=-999};
  enum {kInjMaxNY=10, kInjMaxNZ=8};
  //
  Int_t fRecType;                            // what to generate: normal reco, injection etc.
  //
  // method specific members
  //
  // injection
  Int_t    fInjLr;                            // injection layer
  Int_t    fInjStave;                         // injection stave (counted from layer)
  Int_t    fInjModule;                        // injection module (counted from layer)
  Int_t    fInjModInStave;                    // injection module (counted from stave)
  Double_t fInjX;                             // injection local X
  Double_t fInjZ;                             // injection local Z
  Int_t    fInjHitsN;                         // number of injected pixels
  Int_t    fInjHits[256];                     // bit ids of injected pixels (GetPixelBitL/ExtractPixelL)
  Int_t*   fAssociations[2];                  //! associations [0]: from L1 to L2, [1]: L2 to L1
  Float_t  fInjScale;                         // scaling factor for generation
  Int_t    fInjCurrTrial;                     // current trial id
  Int_t    fInjNTrials[2];                    // number of injections per event for each layer
  Int_t    fInjNSuccess[2];                   // number of successful injected tracklet generation
  AliITSRecPoint fInjCluster;                 // fake cluster to inject
  TH2F*    fInjSPDOcc[2];                     // occupancy in 2 SPD layers
  TBits*   fInjSPDPatt[2];                    // occupancy pattern accounting for cluster types
  Int_t*   fInjModuleClStart[2][4];           //! start of clusters for modules in given Z belt of each stave
  Int_t*   fInjModuleClN[2][4];               //! number of clusters for modules in given Z belt of each stave
  Int_t    fInjModuleTotClN[2][4];            // total number of clusters for modules in given Z belt
  TArrayF  fInjBuffer;                        //! buffer for bg tracklets
  //
  ClassDef(AliITSMultRecBg,0)
};


inline void AliITSMultRecBg::ExtractPixelL(int id, int &stave, int &ladder, int &col, int &row) const
{
  // get sensor data
  enum {kNRow=256,kNCol=32*5,kNLadd=4};
  row = id%kNRow;
  id /= kNRow;
  col = id%kNCol;
  id /= kNCol;
  ladder = id%kNLadd;
  stave = id/kNLadd;
}


#endif
