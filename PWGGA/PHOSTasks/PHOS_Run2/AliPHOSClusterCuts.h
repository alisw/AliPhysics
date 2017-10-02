#ifndef AliPHOSClusterCuts_cxx
#define AliPHOSClusterCuts_cxx

//Author: Daiki Sekihata (Hiroshima University)
//this analsyis class provides cluster selection.

class AliCaloPhoton;

#include "AliAnalysisCuts.h"
#include "AliCaloPhoton.h"

class AliPHOSClusterCuts : public AliAnalysisCuts {

  public:
    AliPHOSClusterCuts(const char *name = "AliPHOSClusterCuts");
    virtual ~AliPHOSClusterCuts(); 

    virtual Bool_t IsSelected(TObject* obj) {return kTRUE;}
    virtual Bool_t IsSelected(TList* /*list*/) {return kTRUE;}

    void SetNsigmaCPV(Double_t nsigma)  {if(nsigma > 0) fUseCPV  = kTRUE; fNsigmaCPV  = nsigma;}
    void SetNsigmaDisp(Double_t nsigma) {if(nsigma > 0) fUseDisp = kTRUE; fNsigmaDisp = nsigma;}
    void SetUseCoreDispersion(Bool_t isCore) {fIsCore=isCore;}

    Bool_t AcceptPhoton(AliCaloPhoton *ph);
    Bool_t AcceptElectron(AliCaloPhoton *ph);//for E/p ratio of electron
    Bool_t AcceptChargedParticle(AliCaloPhoton *ph);//for track matching
    Double_t GetCPVParameter()  {return fNsigmaCPV;}
    Double_t GetDispParameter() {return fNsigmaDisp;}

  protected:
    Bool_t IsNeutral(AliCaloPhoton *ph);
    Bool_t AcceptDisp(AliCaloPhoton *ph);


  private:
    Bool_t fUseCPV;
    Bool_t fUseDisp;
    Double_t fNsigmaCPV;
    Double_t fNsigmaDisp;
    Bool_t fIsCore;

  private:
    AliPHOSClusterCuts(const AliPHOSClusterCuts&);
    AliPHOSClusterCuts& operator=(const AliPHOSClusterCuts&);

    ClassDef(AliPHOSClusterCuts, 11);

};

#endif

