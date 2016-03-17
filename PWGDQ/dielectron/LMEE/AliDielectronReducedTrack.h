#ifndef AliDielectronReducedTrack_H
#define AliDielectronReducedTrack_H

#include "AliLog.h"
#include "AliVParticle.h"


class AliDielectronReducedTrack : public AliVParticle
{
public:
    
    AliDielectronReducedTrack();
    AliDielectronReducedTrack(Double_t px, Double_t py, Double_t pz, Short_t charge, Bool_t IsTagged );
    ~AliDielectronReducedTrack();

    
    virtual Double_t Px()                              const { return fPx; }
    virtual Double_t Py()                              const { return fPy; }
    virtual Double_t Pz()                              const { return fPz; }
    virtual Short_t Charge()                           const { return fCharge; }
    virtual Bool_t IsConversionCandidate()             const { return fIsTagged; }
    
    
	//Kinematics
    virtual Double_t Pt()                              const { AliFatal("Not implemented"); return 0;}
    virtual Double_t Phi()                             const { AliFatal("Not implemented"); return 0;}
    virtual Double_t Eta()                             const { AliFatal("Not implemented"); return 0;}
    virtual Int_t   GetLabel()                         const { AliFatal("Not implemented"); return 0;}
    virtual Double_t P()                               const { AliFatal("Not implemented"); return 0; }
    virtual Bool_t   PxPyPz(Double_t[3])               const { AliFatal("Not implemented"); return 0; }
    virtual Double_t Xv()                              const { AliFatal("Not implemented"); return 0; }
    virtual Double_t Yv()                              const { AliFatal("Not implemented"); return 0; }
    virtual Double_t Zv()                              const { AliFatal("Not implemented"); return 0; }
    virtual Bool_t   XvYvZv(Double_t[3])               const { AliFatal("Not implemented"); return 0; }
    virtual Double_t OneOverPt()                       const { AliFatal("Not implemented"); return 0; }
    virtual Double_t Theta()                           const { AliFatal("Not implemented"); return 0; }
    virtual Double_t E()                               const { AliFatal("Not implemented"); return 0; }
    virtual Double_t M()                               const { AliFatal("Not implemented"); return 0; }
    virtual Double_t Y()                               const { AliFatal("Not implemented"); return 0; }
   
    //PID
    virtual Int_t   PdgCode()                          const { AliFatal("Not implemented"); return 0; }
    virtual const Double_t *PID()                      const { AliFatal("Not implemented"); return 0; }
	
    
private:
    Double_t fPx;//
    Double_t fPy;//
    Double_t fPz;//
    Short_t  fCharge;//
    Bool_t   fIsTagged;//
    
    ClassDef(AliDielectronReducedTrack, 1);
};
#endif
