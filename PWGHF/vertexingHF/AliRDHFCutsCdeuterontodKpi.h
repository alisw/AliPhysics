#ifndef ALIRDHFCUTSCDEUTERONTODKPI_H
#define ALIRDHFCUTSCDEUTERONTODKPI_H
/* Copyright(c) 1998-2010, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */ 

//***********************************************************
/// \class Class AliRDHFCutsCdeuterontodKpi
/// \brief class for cuts on AOD reconstructed c-deuteron -> dKpi
/// \First implementation by copying the AliRDHFCutsXicTopKpi class
/// \inheriting from this class
/// \author Author: J. Norman (jaime.norman@cern.ch)
//***********************************************************

#include "AliRDHFCuts.h"
#include "AliAODPidHF.h"
#include "AliAODRecoDecayHF3Prong.h"
#include "AliRDHFCutsXictopKpi.h"
//#include "AliVertexerTracks.h"

class AliRDHFCutsCdeuterontodKpi : public AliRDHFCutsXictopKpi
{
 public:

  AliRDHFCutsCdeuterontodKpi(const char* name="CutsLctopKpi");
  
  virtual ~AliRDHFCutsCdeuterontodKpi();

  AliRDHFCutsCdeuterontodKpi(const AliRDHFCutsCdeuterontodKpi& source);
  AliRDHFCutsCdeuterontodKpi& operator=(const AliRDHFCutsCdeuterontodKpi& source); 

  // c deuteron specific
 
  using AliRDHFCuts::GetCutVarsForOpt;
  virtual void GetCutVarsForOpt(AliAODRecoDecayHF *d,Float_t *vars,Int_t nvars,Int_t *pdgdaughters){
    return GetCutVarsForOpt(d,vars,nvars,pdgdaughters,0x0);
  }
  virtual void GetCutVarsForOpt(AliAODRecoDecayHF *d,Float_t *vars,Int_t nvars,Int_t *pdgdaughters,AliAODEvent *aod);


  using AliRDHFCuts::IsSelected;
  virtual Int_t IsSelected(TObject* obj,Int_t selectionLevel)
                           {return IsSelected(obj,selectionLevel,0);}
  virtual Int_t IsSelected(TObject* obj,Int_t selectionLevel,AliAODEvent *aod);

  Int_t IsSelectedPID(AliAODRecoDecayHF* obj);

 protected:


private:

  /// \cond CLASSIMP    
  ClassDef(AliRDHFCutsCdeuterontodKpi,1);  /// class for cuts 
  /// \endcond
};

#endif
