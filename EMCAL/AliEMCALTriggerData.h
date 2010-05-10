#ifndef ALIEMCALTRIGGERDATA_H
#define ALIEMCALTRIGGERDATA_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*
EMCal trigger data container
for persistency of produced data presently stored in TTreeD
Author: R. GUERNANE LPSC Grenoble CNRS/IN2P3
*/

#include <TObject.h>
#include <TVector2.h>

//#include <TObjArray.h>
//#include <Riostream.h>

#include <TClonesArray.h>

class AliEMCALTriggerData : public TObject 
{

public:
	             AliEMCALTriggerData();
	virtual     ~AliEMCALTriggerData();
	
	virtual void SetL0Patches(Int_t i, const TClonesArray& patches);
	virtual void SetL0RegionSize(       TVector2 size ) {    fL0RegionSize = size; }
	virtual void SetL0SubRegionSize(    TVector2 size ) { fL0SubRegionSize = size; }
	virtual void SetL0PatchSize(        TVector2 size ) {     fL0PatchSize = size; }
	
	virtual void SetL1GammaPatches(const TClonesArray& patches);
	virtual void SetL1JetPatches(const TClonesArray& patches);

	virtual void SetL1Region(Int_t**& region);
	virtual void SetL1V0(const Int_t*& arr);

	virtual void SetL1RegionSize(         TVector2 size ) {         fL1RegionSize = size; }
	virtual void SetL1GammaPatchSize(     TVector2 size ) {     fL1GammaPatchSize = size; }
	virtual void SetL1GammaSubRegionSize( TVector2 size ) { fL1GammaSubRegionSize = size; }
	virtual void SetL1JetPatchSize(       TVector2 size ) {       fL1JetPatchSize = size; }
	virtual void SetL1JetSubRegionSize(   TVector2 size ) {   fL1JetSubRegionSize = size; }

	virtual void          L0Patches(       TClonesArray& patches ) const { patches =  *fL0Patches;            }
	virtual TClonesArray* L0Patches(                             ) const {     return  fL0Patches;            }
 	virtual void          L0RegionSize(        TVector2 size     ) const { size =      fL0RegionSize;         }
 	virtual TVector2      L0RegionSize(                          ) const {     return  fL0RegionSize;         }
	virtual void          L0PatchSize(         TVector2 size     ) const { size =      fL0PatchSize;          }
	virtual TVector2      L0PatchSize(                           ) const {     return  fL0PatchSize;          }
	virtual void          L0SubRegionSize(     TVector2 size     ) const { size =      fL0SubRegionSize;      }
	virtual TVector2      L0SubRegionSize(                       ) const {     return  fL0SubRegionSize;      }
	virtual void          L0NPatches( Int_t arr[32]              ) const { for (Int_t i=0;i<32;i++) arr[i] = fL0NPatches[i]; }

	virtual void          L1GammaPatches(  TClonesArray& patches ) const { patches =  *fL1GammaPatches;       }
	virtual TClonesArray* L1GammaPatches(                        ) const {     return  fL1GammaPatches;       }
	virtual void          L1JetPatches(    TClonesArray& patches ) const { patches =  *fL1JetPatches;         }
	virtual TClonesArray* L1JetPatches(                          ) const {     return  fL1JetPatches;         }
	virtual void          L1Region( Int_t arr[][64]              ) const { for (Int_t i=0;i<48;i++) for (Int_t j=0;j<64;j++) { arr[i][j] = fL1Region[i][j]; } }
	virtual Int_t*        L1V0(                                  )       {     return &fL1V0[0];              }
	virtual void          L1RegionSize(         TVector2& size   ) const { size =      fL1RegionSize;         }
	virtual TVector2      L1RegionSize(                          ) const {     return  fL1RegionSize;         }
	virtual void          L1GammaPatchSize(     TVector2& size   ) const { size =      fL1GammaPatchSize;     }
	virtual TVector2      L1GammaPatchSize(                      ) const {     return  fL1GammaPatchSize;     }
	virtual void          L1GammaSubRegionSize( TVector2& size   ) const { size =      fL1GammaSubRegionSize; }
	virtual TVector2      L1GammaSubRegionSize(                  ) const {     return  fL1GammaSubRegionSize; }
	virtual void          L1JetPatchSize(       TVector2& size   ) const { size =      fL1JetPatchSize;       }
	virtual TVector2      L1JetPatchSize(                        ) const {     return  fL1JetPatchSize;       }
	virtual void          L1JetSubRegionSize(   TVector2& size   ) const { size =      fL1JetSubRegionSize;   }
	virtual TVector2      L1JetSubRegionSize(                    ) const {     return  fL1JetSubRegionSize;   }

	virtual void Scan() const;

	virtual void Reset();
	
private:

    AliEMCALTriggerData(const AliEMCALTriggerData& rhs);            // NOT implemented
	AliEMCALTriggerData& operator=(const AliEMCALTriggerData& rhs); // NOT implemented
	
	TClonesArray*  fL0Patches;            // array of patches  
	Int_t          fL0NPatches[32];
	TVector2       fL0RegionSize;         // region size in units of fast or
	TVector2       fL0SubRegionSize;      // subregion size in units of fast or
	TVector2       fL0PatchSize;          // patch size in units of subregion
	
	TClonesArray*  fL1GammaPatches;       // array of patches  
	TClonesArray*  fL1JetPatches;         // array of patches  
	Int_t          fL1Region[48][64];     // STU FastOR 48-by-124
	Int_t          fL1V0[2];              // V0A V0C multiplicity estimates	
	TVector2       fL1RegionSize;         // region size in units of fast or
	TVector2       fL1GammaPatchSize;     // patch size in units of subregion
	TVector2       fL1GammaSubRegionSize; // subregion size in units of fast or
	TVector2       fL1JetPatchSize;       // patch size in units of subregion
	TVector2       fL1JetSubRegionSize;   // subregion size in units of fast or 
	
	ClassDef(AliEMCALTriggerData,1)
};

#endif
