#ifndef AliFParticle_H
#define AliFParticle_H

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// AliFParticle                                                         //
//                                                                      //
// Graphics interface to event generators particle                      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#ifndef ROOT_TNamed
#include <TNamed.h>
#endif

class TMCParticle;
class TPolyLine3D;
class TList;
class AliFDisplay;

class AliFParticle : public TNamed {

private:
   TList            *fParticles;            //List for particles
   AliFDisplay      *fDisplay;              //pointer to AliFDisplay object
   TMCParticle      *fMCParticle;           //pointer to selected particle
   TPolyLine3D      *fLine;                 //pointer to line3D
   
public:
                     AliFParticle() {;}
                     AliFParticle(const char *name);
   virtual          ~AliFParticle();
   virtual void      Clear(Option_t *option="");
   virtual void      Delete(Option_t *option="");
   virtual Int_t     DistancetoPrimitive(Int_t px, Int_t py);
   AliFDisplay      *Display() {return fDisplay;}
   virtual void      ExecuteEvent(Int_t event, Int_t px, Int_t py);
   virtual char     *GetObjectInfo(Int_t px, Int_t py);
   TPolyLine3D      *HelixCurve(Float_t field, Float_t pmom, Float_t *vin);
   virtual void      HelixStep(Float_t field, Float_t step, Float_t pmom, Float_t *vin, Float_t *vout);
   virtual void      Paint(Option_t *option="");
   virtual void      SetLineAttributes(); // *MENU*
   virtual void      SizeParticles() const;

   ClassDef(AliFParticle, 0)   //Graphics interface to event generators particle
};

#endif
