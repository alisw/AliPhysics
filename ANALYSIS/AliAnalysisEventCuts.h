#ifndef ALIANALYSISEVENTCUTS_H
#define ALIANALYSISEVENTCUTS_H
/*  See cxx source for full Copyright notice */


/* $Id$ */

//-------------------------------------------------------------------------
//                       Class AliAnalysisEventCuts
//   This is the class for the cuts in event & track level
//
//    Origin: Panos Christakoglou, UOA-CERN, Panos.Christakoglou@cern.ch
//-------------------------------------------------------------------------

#include <TObject.h>

class TPaveText;
class AliESD;

class AliAnalysisEventCuts : public TObject
{
 public:
  AliAnalysisEventCuts();
  
  ~AliAnalysisEventCuts();

  void Reset();
  
  void SetPrimaryVertexXRange(Float_t r1, Float_t r2);
  void SetPrimaryVertexYRange(Float_t r1, Float_t r2);
  void SetPrimaryVertexZRange(Float_t r1, Float_t r2);
  void SetMultiplicityRange(Int_t n1, Int_t n2);
 
  Bool_t IsAccepted(AliESD *esd);

  TPaveText *GetEventCuts();
  void PrintEventCuts(); 
  void GetEventStats();
  void GetMultStats();
  void GetVxStats();
  void GetVyStats();
  void GetVzStats();
 
 private:
  Float_t fVxMin, fVxMax;  //Definition of the range of the Vx
  Float_t fVyMin, fVyMax;  //Definition of the range of the Vy
  Float_t fVzMin, fVzMax;  //Definition of the range of the Vz
  Int_t fMultMin, fMultMax;  //Definition of the range of the multiplicity

  Int_t fMult;  //Number of events rejected due to multiplicity cut
  Int_t fVx;  //Number of events rejected due to Vx cut
  Int_t fVy;  //Number of events rejected due to Vy cut
  Int_t fVz;  //Number of events rejected due to Vz cut
  Int_t fTotalEvents;  //Total number of events
  Int_t fAcceptedEvents;  //Total number of events accepted

  Int_t fFlagMult; //Flag that shows if the multiplicity cut was imposed
  Int_t fFlagVx; //Flag that shows if the Vx cut was imposed
  Int_t fFlagVy; //Flag that shows if the Vy cut was imposed
  Int_t fFlagVz; //Flag that shows ifthe Vz cut was imposed
 
  ClassDef(AliAnalysisEventCuts, 1)
} ;

#endif
