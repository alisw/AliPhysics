#ifndef AliFDisplay_H
#define AliFDisplay_H

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// AliFDisplay                                                          //
//                                                                      //
// Utility class to display ALICE outline, tracks, clusters, jets,..    //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#ifndef AliFVirtualDisplay_H
#include "AliFVirtualDisplay.h"
#endif

class TCanvas;
class TPad;
class TArc;
class TTUBE;
class TNode;
class TPolyLine3D;
class TList;
class AliFParticle;

class AliFDisplay : public AliFVirtualDisplay {

private:
   Bool_t            fDrawAllViews;         //Flag True if AllViews selected
   Bool_t            fDrawParticles;        //Flag True if particles to be drawn
   Float_t           fPTcut;                //PT cut to display objects
   Float_t           fPTcutEGMUNU;          //PT cut for Electrons, Gammas, MUons, Neutrinos
   Float_t           fRin;                  //Inner ALIAS radius
   Float_t           fRout;                 //Outer ALIAS radius
   Float_t           fZin;                  //Inner ALIAS length along Z
   Float_t           fZout;                 //Outer ALIAS length along Z
   Float_t           fTheta;                //Viewing angle theta
   Float_t           fPhi;                  //Viewing angle phi
   TCanvas          *fCanvas;               //Pointer to the display canvas
   TPad             *fTrigPad;              //Pointer to the trigger pad 
   TPad             *fButtons;              //Pointer to the buttons pad
   TPad             *fPad;                  //Pointer to the event display main pad
   TTUBE            *fTubin;                //Inner tube
   TTUBE            *fTubout;               //outer tube
   TNode            *fNodin;                //Node for detector outline
   TList            *fFruits;               //List for fruits
   AliFParticle     *fParticle;             //Pointer to Particle graphics manager
   
public:
                     AliFDisplay();
                     AliFDisplay(const char *title);
   virtual          ~AliFDisplay();
   virtual Bool_t    AllViews() {return fDrawAllViews;}
   virtual void      Clear(Option_t *option="");
   virtual void      DisplayButtons();
   virtual Int_t     DistancetoPrimitive(Int_t px, Int_t py);
   virtual void      Draw(Option_t *option="");
   virtual void      DrawAllViews();
   Bool_t            DrawParticles() {return fDrawParticles;}
   virtual void      DrawTitle(Option_t *option="");
   virtual void      DrawView(Float_t theta, Float_t phi);
   virtual void      DrawViewGL();
   virtual void      DrawViewX3D();
   virtual void      ExecuteEvent(Int_t event, Int_t px, Int_t py);
   virtual void      GetEvent(Int_t event); //*MENU*
   TNode            *Nodin() {return fNodin;}
   TTUBE            *Tubin() {return fTubin;}
   TPad             *Pad() {return fPad;}
   virtual void      Paint(Option_t *option="");
   virtual void      PaintFruit(TObject *obj, Float_t eta, Float_t phi, Float_t pt, Int_t type, Option_t *option="");
   virtual void      PaintParticles(Option_t *option="");
   Float_t           PTcut() {return fPTcut;}
   Float_t           PTcutEGMUNU() {return fPTcutEGMUNU;}
   Float_t           Rin() {return fRin;}
   Float_t           Rout() {return fRout;}
   virtual void      SetDrawParticles(Bool_t draw=kTRUE) {fDrawParticles=draw;} // *MENU*
   virtual void      SetPTcut(Float_t ptcut=0.4); // *MENU*
   virtual void      SetPTcutEGMUNU(Float_t ptcut=5); // *MENU*
   virtual void      SetGeometry(Float_t rin=115); // *MENU*
   virtual void      SetView(Float_t theta, Float_t phi);
   virtual void      ShowNextEvent(Int_t delta=1);
   virtual void      SizeFruit() const;
   virtual void      SizeParticles() const;
   Float_t           Zin() {return fZin;}
   Float_t           Zout() {return fZout;}

   ClassDef(AliFDisplay, 0)   //Utility class to display ALIAS outline, tracks, clusters, jets,..
};

#endif







