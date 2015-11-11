#ifndef ALIEMCALTRIGGERSTUDCSCONFIG_H
#define ALIEMCALTRIGGERSTUDCSCONFIG_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*


EMCAL STU DCS parameters to be stored in OCDB
Author: R. GUERNANE LPSC Grenoble CNRS/IN2P3
*/

#include "TObject.h"

class TVector2;

class AliEMCALTriggerSTUDCSConfig : public TObject 
{
 public:

  AliEMCALTriggerSTUDCSConfig();
  virtual ~AliEMCALTriggerSTUDCSConfig() {};
	  
  void    SetG(Int_t i, Int_t j, Int_t gv) { fG[i][j]    = gv; }
  void    SetJ(Int_t i, Int_t j, Int_t jv) { fJ[i][j]    = jv; }
  void    SetRawData(Int_t rd)             { fGetRawData = rd; }
  void    SetRegion(Int_t rg)              { fRegion     = rg; }
  void    SetFw(Int_t fv)                  { fFw         = fv; }
  void    SetPHOSScale(int iscale, int val) { fPHOSScale[iscale] = val; }
	
  Int_t   GetG(int i, int j) const { return fG[i][j];    }
  Int_t   GetJ(int i, int j) const { return fJ[i][j];    }
  Int_t   GetRawData()       const { return fGetRawData; }
  Int_t   GetRegion()        const { return fRegion;     }
  Int_t   GetFw()            const { return fFw;         }
  Int_t   GetPHOSScale(Int_t iscale) const { return fPHOSScale[iscale]; }

  void    GetSegmentation(TVector2& v1, TVector2& v2, TVector2& v3, TVector2& v4) const;
	
protected:

	AliEMCALTriggerSTUDCSConfig(const AliEMCALTriggerSTUDCSConfig &cd);
	AliEMCALTriggerSTUDCSConfig &operator=(const AliEMCALTriggerSTUDCSConfig &cd);

private:
	
  Int_t   fG[3][2];       // GA,B,C
  Int_t   fJ[3][2];       // JA,B,C
  Int_t   fGetRawData;    // GetRawData
  Int_t   fRegion;        // Region
  Int_t   fFw;            // Fw
  Int_t   fPHOSScale[4];  // PHOS scale factors
  
  ClassDef(AliEMCALTriggerSTUDCSConfig,3) //
};
#endif

