
//////////////////////////////////////////////////////////////////////////
//                                                                      //
// AliFParticle                                                         //
//                                                                      //
//  Graphics interface to event generators particle                     //
//////////////////////////////////////////////////////////////////////////

#include <TROOT.h>
#include "TMCParticle.h"
#include "TClonesArray.h"
#include "TPolyLine3D.h"
#include "TCanvas.h"
#include <TList.h>
#include <TMath.h>

#include "AliFParticle.h"
#include "AliFMCMaker.h"
#include "AliFast.h"
#include "AliFDisplay.h"

const Int_t kRECONS = BIT(16);

ClassImp(AliFParticle)


//_____________________________________________________________________________
AliFParticle::AliFParticle(const char * name) :TNamed(name,name)
{
   // Create list to support list of particles
   fParticles = new TList();
   fDisplay   = (AliFDisplay*)gAliFast->Display();
}

//_____________________________________________________________________________
AliFParticle::~AliFParticle()
{
   if (fParticles) fParticles->Delete();
   delete fParticles;
}

//_____________________________________________________________________________
void AliFParticle::Clear(Option_t *)
{
//    Delete graphics temporary objects

   fParticles->Delete();

}

//_____________________________________________________________________________
void AliFParticle::Delete(Option_t *)
{
//    Dummy

}

//_____________________________________________________________________________
Int_t AliFParticle::DistancetoPrimitive(Int_t px, Int_t py)
{
    // scan list of particles
   TClonesArray *particles = gAliFast->MCMaker()->Fruits();
   Int_t uid, dist;
   TIter next(fParticles);
   TPolyLine3D *line;
   while((line=(TPolyLine3D*)next())) {
      dist = line->DistancetoPrimitive(px, py);
      if (dist < 2) {
         uid = line->GetUniqueID();
         fMCParticle = (TMCParticle*)particles->UncheckedAt(uid);
         if (!fMCParticle) continue;
         fLine = line;
         SetName(fMCParticle->GetName());
         gPad->SetSelected(this);
         gPad->SetCursor(kCross);
         return 0;
      }
   }
   return 999;
}

//______________________________________________________________________________
void AliFParticle::ExecuteEvent(Int_t event, Int_t , Int_t )
{
   switch (event) {

   case kButton1Down:
      gGXW->SetLineColor(-1);
      gPad->AbsCoordinates(kTRUE);
      fLine->SetLineColor(6);
      fLine->SetLineWidth(6);
      fLine->Paint(); 
      break;

   case kMouseMotion:
      break;

   case kButton1Motion:
      break;

   case kButton1Up:
      gGXW->SetLineColor(-1);
      fLine->SetLineColor(kYellow);
      fLine->SetLineWidth(1);
      fLine->Paint(); 
      gPad->AbsCoordinates(kFALSE);
   }
}

//______________________________________________________________________________
char *AliFParticle::GetObjectInfo(Int_t , Int_t )
{
   static char info[100];
   sprintf(info,"px=%f, py=%f, pz=%f, E=%f",
            fMCParticle->GetPx(),
            fMCParticle->GetPy(),
            fMCParticle->GetPz(),
            fMCParticle->GetEnergy());


   return info;
}


//______________________________________________________________________________
TPolyLine3D *AliFParticle::HelixCurve(Float_t field, Float_t pmom, Float_t *vin)
{
//    Estimate step size in function of field.
//    Create a 3-D polyline with points computed with this step size

   Float_t step = 10;
   const Int_t kMAXSTEP = 2000;
   Float_t sx[kMAXSTEP], sy[kMAXSTEP], sz[kMAXSTEP];
   if (pmom > 0)   step = 10*pmom;
   if (step > 100) step = 100;
   if (step <  .5) step = .5;

   Float_t vout[6];
   Int_t i,j;
   sx[0] = vin[0];
   sy[0] = vin[1];
   sz[0] = vin[2];
   Int_t nsteps = 1;
   Float_t rin  = Display()->Rin();
   Float_t zout = Display()->Zout();
   for (i=1;i<kMAXSTEP;i++) {
      HelixStep(field,step,pmom,vin,vout);
      if (TMath::Abs(vout[2]) > zout) break;
      if (vout[0]*vout[0] + vout[1]*vout[1] > 1.1*rin*rin) break;
      sx[nsteps] = vout[0];
      sy[nsteps] = vout[1];
      sz[nsteps] = vout[2];
      nsteps++;
      for (j=0;j<6;j++) vin[j] = vout[j];
   }
   if (nsteps < 2) return 0;
   TPolyLine3D *line = new TPolyLine3D(nsteps,sx, sy,sz);
   line->SetBit(kCanDelete);
   return line;
}

//______________________________________________________________________________
void AliFParticle::HelixStep(Float_t field, Float_t step, Float_t pmom, Float_t *vin, Float_t *vout)
{
//     extrapolate track with parameters in vector vin in a constant field
//     oriented along Z axis (in tesla/meters).
//     Output in vector vout
//     vin[0-->6] = x,y,z,px,py,pz
//     translated to C++ from GEANT3 routine GHELX3

   Float_t sint, sintt, tsint, cos1t, sin2;
//      units are tesla,centimeters,gev/c
   const Float_t ec = 2.9979251e-3;
   Float_t h4  = field*ec;
   Float_t hp  = vin[2];
   Float_t tet = -h4*step/pmom;
   if (TMath::Abs(tet) > 0.15) {
      sint  = TMath::Sin(tet);
      sintt = sint/tet;
      tsint = (tet-sint)/tet;
      sin2  = TMath::Sin(0.5*tet);
      cos1t = 2*sin2*sin2/tet;
   } else {
      tsint = tet*tet/6;
      sintt = 1 - tsint;
      sint  = tet*sintt;
      cos1t = 0.5*tet;
   }
   Float_t f1 = step*sintt;
   Float_t f2 = step*cos1t;
   Float_t f3 = step*tsint*hp;
   Float_t f4 = -tet*cos1t;
   Float_t f5 = sint;
   Float_t f6 = tet*cos1t*hp;

   vout[0] = vin[0] + (f1*vin[3] - f2*vin[4]);
   vout[1] = vin[1] + (f1*vin[4] + f2*vin[3]);
   vout[2] = vin[2] + (f1*vin[5] + f3);

   vout[3] = vin[3] + (f4*vin[3] - f5*vin[4]);
   vout[4] = vin[4] + (f4*vin[4] + f5*vin[3]);
   vout[5] = vin[5] + (f4*vin[5] + f6);
}

//_____________________________________________________________________________
void AliFParticle::Paint(Option_t *option)
{
//    Paint particles generated by AliFMCMaker
//    Only particles above fPTcut are drawn
//    Particle trajectory is computed along an helix in a constant field


   // clean list of particles
   fParticles->Delete();

   TClonesArray *particles = gAliFast->MCMaker()->Fruits();
   Int_t nparticles = particles->GetEntriesFast();
   TMCParticle *part;
   Float_t pmom, vx, vy, vz, pt;
   Float_t vin[6];
   Float_t field = 2; // 2 tesla
   Int_t KF, aKF, charge;
   for (Int_t i=0;i<nparticles;i++) {
      part = (TMCParticle*)particles->UncheckedAt(i);     
      if (part->GetKS() != 1) continue;
      KF = part->GetKF();
      aKF = TMath::Abs(KF);
      charge = gAliFast->MCMaker()->Charge(KF);
      pt = TMath::Sqrt(part->GetPx()*part->GetPx() + part->GetPy()*part->GetPy());
      if (pt < gAliFast->Display()->PTcut()) continue;
      if (charge == 0 && pt < gAliFast->Display()->PTcutEGMUNU()) continue;
      pmom = TMath::Sqrt(part->GetPx()*part->GetPx() + part->GetPy()*part->GetPy() + part->GetPz()*part->GetPz());
      vx = part->GetVx();
      vy = part->GetVy();
      vz = part->GetVz();
      vin[0] = vx;
      vin[1] = vy;
      vin[2] = vz;
      vin[3] = part->GetPx()/pmom;
      vin[4] = part->GetPy()/pmom;
      vin[5] = part->GetPz()/pmom;
      TPolyLine3D *line = HelixCurve(charge*field/3, pmom, vin);
      if (line == 0) continue;
      fParticles->Add(line);
      if (part->GetLineColor() == 1) {
         if (charge == 0) part->SetLineStyle(2);
         part->SetLineColor(kYellow);
         if (aKF == 11) part->SetLineColor(kRed);
         if (aKF == 13) part->SetLineColor(kMagenta);
         if (aKF == 22) part->SetLineColor(kGreen);
      }
      line->SetUniqueID(i);
      if (!part->TestBit(kRECONS)) {
         part->SetLineColor(7);
         part->SetLineWidth(1);
      }
      line->SetLineColor(part->GetLineColor());
      line->SetLineStyle(part->GetLineStyle());
      line->SetLineWidth(part->GetLineWidth());
      line->Paint(option);
   }
}

//______________________________________________________________________________
void AliFParticle::SetLineAttributes()
{
//*-*-*-*-*-*-*-*-*Invoke the DialogCanvas Line attributes*-*-*-*-*-*-*
//*-*              =======================================

   gROOT->SetSelectedPrimitive(fMCParticle);
   fMCParticle->SetLineAttributes();
}


//______________________________________________________________________________
void AliFParticle::SizeParticles() const
{
   TIter next(fParticles);
   TPolyLine3D *line;
   while((line=(TPolyLine3D*)next())) {
      line->Sizeof3D();
   }
}








