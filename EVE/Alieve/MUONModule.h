/* HEAD11Jul06 */
#ifndef ALIEVE_MUONModule_H
#define ALIEVE_MUONModule_H

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// The main AliEVE drawing module for the MUON detector                 //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include <Reve/QuadSet.h>
#include <Reve/RenderElement.h>

#include <Alieve/MUONDigitsInfo.h>

namespace Alieve {

class MUONModule : public Reve::RenderElement,
                   public Reve::QuadSet
{

public:

  MUONModule(const Text_t* n="MUONModule", const Text_t* t=0, Color_t col=2);
  MUONModule(Int_t id, Int_t cath, MUONDigitsInfo* info, Bool_t dig, Bool_t clus, Color_t col=2);
  virtual ~MUONModule() {}

  MUONModule(const MUONModule&);           
  MUONModule& operator=(const MUONModule&);

  virtual void SetID(Int_t id, Int_t cath);

private:
    
  void LoadQuadsChambers(Int_t chamber1, Int_t chamber2, Int_t delElemId = -1, Int_t cat = -1);
  void LoadQuadsDigits();
  void LoadQuadsClusters();
  void LoadQuadsTracks(Int_t id);

  virtual void InitModule();

protected:
 
  MUONDigitsInfo* fInfo;            // pointer to the tree interface
  
  Int_t       fID;                  // detector element id for drawing
  Int_t       fCath;                // cathode number
 
  Bool_t      fShowDigits;          // digits drawing flag
  Bool_t      fShowClusters;        // clusters drawing flag
  Bool_t      fShowTracks;          // tracks drawing flag

  Color_t     fFrameCol;            // main color

  Int_t       fDetElemId;           // detector element id

  ClassDef(MUONModule,1);

};
}

#endif
