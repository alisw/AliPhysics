#ifndef Wa98Event_h
#define Wa98Event_h

// $Id$
 
#include "AliEvent.h"
#include "AliCalorimeter.h"
 
class Wa98Event : public AliEvent
{
 public:
  Wa98Event();                                              // Default constructor
  Wa98Event(Int_t n);                                       // Create an event to hold initially n tracks
  virtual ~Wa98Event();                                     // Default destructor
  Wa98Event(Wa98Event& evt);                                // Copy constructor
  void Reset();                                             // Reset the complete event
  void SetTrig(Int_t trig);                                 // Set the trigger class
  Int_t GetTrig();                                          // Provide the trigger class
  void SetWeight(Int_t w);                                  // Set the event weight factor
  Int_t GetWeight();                                        // Provide the event weight factor
  void SetZdc(Float_t zdc);                                 // Set ZDC signal
  Float_t GetZdc();                                         // Provide ZDC signal
  void SetMiracE(Float_t tot,Float_t em,Float_t had);       // Set the total, EM and hadr. MIRAC signals
  Float_t GetEmir();                                        // Provide the total MIRAC signal
  Float_t GetEmire();                                       // Provide the MIRAC EM signal
  Float_t GetEmirh();                                       // Provide the MIRAC hadronic signal
  void SetMiracEt(Float_t tot,Float_t em,Float_t had);      // Set the total, EM and hadr. MIRAC Et signals
  Float_t GetEtm();                                         // Provide the total MIRAC Et signal
  Float_t GetEtme();                                        // Provide the MIRAC EM Et signal
  Float_t GetEtmh();                                        // Provide the MIRAC hadronic Et signal
  void InitLeda(AliCalorimeter* cal);                       // Set positions and flag bad modules for LEDA part
  void ClusterLeda(AliCalorimeter* cal,Int_t n=2,Int_t mode=1); // Group (n rings) modules into clusters for LEDA part
  void VetoLeda(AliCalorimeter* cal,Float_t dtheta=0.2,Float_t dphi=1); // Associate veto hits with LEDA clusters

 protected:
  Int_t fTrig;        // The trigger class  1=nsc 3=cen 5=per 6=mbias 7=beam 8=inbeam ped
  Int_t fWeight;      // The event weight to account for the downscaling factor
  Float_t fZdc;       // The ZDC signal in GeV
  Float_t fEmir;      // The total MIRAC signal in GeV
  Float_t fEmire;     // The MIRAC EM signal in GeV
  Float_t fEmirh;     // The MIRAC hadronic signal in GeV
  Float_t fEtm;       // The total MIRAC Et signal in GeV
  Float_t fEtme;      // The MIRAC EM Et signal in GeV
  Float_t fEtmh;      // The MIRAC hadronic Et signal in GeV

 private:
  void SetPositionsLedaup(AliCalorimeter* cal);  // Set module positions for upper LEDA
  void SetPositionsLedalw(AliCalorimeter* cal);  // Set module positions for lower LEDA
  void SetBadModulesLedaup(AliCalorimeter* cal); // Flag bad modules for upper LEDA
  void SetBadModulesLedalw(AliCalorimeter* cal); // Flag bad modules for lower LEDA

 ClassDef(Wa98Event,4) // Creation and investigation of a Wa98 physics event.
};
#endif
