#ifndef ALIMESPPCOLTASK_H
#define ALIMESPPCOLTASK_H

////////////////////////////////////////////////////////////////////////////
//  PP collision task for Multiplicity and Event Shape group                      //
//  Authors:                                                              //
//    Madalina Tarzila <mtarzila@niham.nipne.ro>                          //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#ifndef ALIMESBASETASK_H
#include "AliMESbaseTask.h"
#include "AliVParticle.h"
#include "AliBasicParticle.h"
#endif
#define NMAXMULT 100

class AliMEStrackInfo;
class AliEventPoolManager;
class AliMESppColTask : public AliMESbaseTask
{
public:
  class AliMESppColTaskExchange: public TObject
  {
  public:
    AliMESppColTaskExchange();
    void        Add(Float_t deta, Float_t dphi);
  private:
    Int_t   fN;
    Float_t fDEta[NMAXMULT]; //[fN]
    Float_t fDPhi[NMAXMULT]; //[fN]
    ClassDef(AliMESppColTaskExchange, 1)
  };
  
  AliMESppColTask();
  AliMESppColTask(const char *name);
  virtual ~AliMESppColTask();
  
    //
  TObjArray*   FindLeadingObjects(TObjArray* obj);                                 // Looks for leading track 
  TObjArray*   CloneTracks(TObjArray* tracks);
  void         QSortTracks(TObjArray &a, Int_t first, Int_t last);               // Sort by pT an array of AliVParticles 
  Double_t     RangePhi(Double_t DPhi);
  Bool_t       DefineMixedEventPool(); // Definition of the Event pool parameters
  void FillCorrelationSE(Double_t MultipOrCent, TObjArray* selectedArray);
  void FillCorrelationMixing(Double_t MultipOrCentMix, Double_t Zvtx, Double_t poolmax,Double_t poolmin, TObjArray*selectedArray);

  Bool_t         BuildQAHistos();
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *opt);
  virtual Bool_t PostProcess();

private:
  AliMESppColTask(const AliMESppColTask&);
  AliMESppColTask& operator=(const AliMESppColTask&);

protected:
  AliEventPoolManager      * fPoolMgr;                  //  event pool manager for Event Mixing

  ClassDef(AliMESppColTask, 1)            // PP collision task for the Multi Event Shape
};
//---------------------------------------------------------------------------------------

#endif