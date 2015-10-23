/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id$ */

#include "AliGenHBTosl.h" 
#include "AliLog.h"
#include <iostream>
#include <fstream>

//__________________________________________________________
/////////////////////////////////////////////////////////////
//                                                         //
//  class AliGenHBTosl                                     //
//                                                         //
//  Genarator simulating particle correlations             //
//                                                         //
//  The main idea of the generator is to produce particles //
//  according to some distribution of two particle         //
//  property. In HBT they are qout,qsie and qlong.         //
//  In order to be able to generate signal that produces   //
//  given two particle correlation background must be      //
//  known before in order to produce the shape of signal   //
//  to randomize given distribution from.                  //
//                                                         //
//  The generator works as follows:                        //
//  1. Coarse Background (fQCoarseBackground) is generated //
//     ade  from the particles                             //
//     given by the external generator (variable           //
//     fGenerator) by the mixing technique.                //
//  2. Coarse signal is prduced by multiplying Coarse      //
//     background by a required function                   //
//     See method FillCoarseSignal                         //
//  3. Signal is randomized out of the coarse signal       //
//     histogram (two particle property). First particle   //
//     is taken from the external generator, and the       //
//     second one is CALCULATED on the basis of the first  //
//     one and the two particle property (qout,qside,qlong)//
//     Background is made by the mixing out of the         //
//     genereted signal events.                            //
//     This step is cotinued up to the moment signal       //
//     histogram has enough statistics (data member        //
//     fMinFill)                                           //
//     See method StartSignalPass1()                       //
//  4. chi is calculated for each bin (chiarray variqable) // 
//     (not the chi2 because sign is important)            //
//     Two particle prioperty                              //
//     (qout,qside,qlong) is chosen at the points that     //
//     chi is the smallest. First particle is taken from   //
//     the the external generator (fGenerator) and second's /
//     momenta are caclulated out of their momenta and     //
//     (qout,qside,qlong). Background is updated           //
//     continuesely for all the events. This step is       //
//     continued until stability conditions are fullfiled  //
//     or maximum number of iteration is reached.          //
//  5. The same as step 4 but events are stored.           //
//                                                         //
////////////////////////////////////////////////////////////

#include <TCanvas.h>


#include <TH3D.h>
#include <TList.h>
#include <TPDGCode.h>
#include <TParticle.h>
#include <AliStack.h>
#include <TMath.h>
#include <TVector3.h>
#include <TStopwatch.h>
#include <TFile.h>

#include "AliGenCocktailAfterBurner.h"
#include "AliGeVSimParticle.h"
#include "AliGenGeVSim.h"
#include "AliGenHIJINGpara.h"


/***********************************************************/
using std::cout;
using std::endl;
using std::ofstream;
using std::ios;
ClassImp(AliGenHBTosl)

AliGenHBTosl::AliGenHBTosl():
 AliGenerator(),
 fQCoarseBackground(0x0),
 fQCoarseSignal(0x0),
 fQSignal(0x0),
 fQBackground(0x0),
 fQSecondSignal(0x0),
 fQSecondBackground(0x0),
 fQRange(0.06),
 fQNBins(60),
 fGenerator(0x0),
 fStackBuffer(0),
 fBufferSize(5),
 fNBinsToScale(Int_t(fQNBins*0.1)),
 fDebug(0),
 fSignalShapeCreated(0),
 fMaxIterations(1000),
 fMaxChiSquereChange(0.01),
 fMaxChiSquerePerNDF(1.5), 
 fQRadius(8.0),
 fPID(kPiPlus),
 fSamplePhiMin(-0.01),
 fSamplePhiMax(TMath::TwoPi()+0.01),
 fSignalRegion(0.0),
 fMinFill(1000),
 fSwapped(0),
 fLogFile(0x0)
{
//default constructor
}
/***********************************************************/

AliGenHBTosl::AliGenHBTosl(Int_t n, Int_t pid):
 AliGenerator(n),
 fQCoarseBackground(0x0),
 fQCoarseSignal(0x0),
 fQSignal(0x0),
 fQBackground(0x0),
 fQSecondSignal(0x0),
 fQSecondBackground(0x0),
 fQRange(0.06),
 fQNBins(60),
 fGenerator(0x0),
 fStackBuffer(0),
 fBufferSize(5),
 fNBinsToScale(Int_t(fQNBins*0.1)),
 fDebug(0),
 fSignalShapeCreated(kFALSE),
 fMaxIterations(1000),
 fMaxChiSquereChange(0.01),
 fMaxChiSquerePerNDF(1.5),
 fQRadius(8.0),
 fPID(pid),
 fSamplePhiMin(-0.01),
 fSamplePhiMax(TMath::TwoPi()+0.01),
 fSignalRegion(0.0),
 fMinFill(1000),
 fSwapped(0),
 fLogFile(0x0)
{
//default constructor
}

AliGenHBTosl::AliGenHBTosl(const AliGenHBTosl & hbt):
 AliGenerator(-1),
 fQCoarseBackground(0x0),
 fQCoarseSignal(0x0),
 fQSignal(0x0),
 fQBackground(0x0),
 fQSecondSignal(0x0),
 fQSecondBackground(0x0),
 fQRange(0.06),
 fQNBins(60),
 fGenerator(0x0),
 fStackBuffer(0),
 fBufferSize(5),
 fNBinsToScale(Int_t(fQNBins*0.1)),
 fDebug(0),
 fSignalShapeCreated(kFALSE),
 fMaxIterations(1000),
 fMaxChiSquereChange(0.01),
 fMaxChiSquerePerNDF(1.5),
 fQRadius(8.0),
 fPID(kPiPlus),
 fSamplePhiMin(-0.01),
 fSamplePhiMax(TMath::TwoPi()+0.01),
 fSignalRegion(0.0),
 fMinFill(1000),
 fSwapped(0),
 fLogFile(0x0)
{
// Copy constructor
    hbt.Copy(*this);
}
/***********************************************************/

AliGenHBTosl::~AliGenHBTosl()
{
//destructor
 delete fQCoarseSignal;
 delete fQCoarseBackground;
 delete fQSignal;
 delete fQBackground;
 delete fGenerator; 
 delete fQSecondSignal;
 delete fQSecondBackground;
 delete fLogFile;
}
/***********************************************************/

void AliGenHBTosl::Init()
{
  //Initializes generator
  if (fGenerator == 0x0)
   { 
   
     AliGenHIJINGpara* bkggen = new AliGenHIJINGpara(fNpart*4);
     fGenerator = bkggen;

/*     
     AliGenGeVSim * gevsim = new AliGenGeVSim(0.0);
     AliGeVSimParticle* kplus = new AliGeVSimParticle(fPID,1,fNpart, 0.17, 0.9);
     gevsim->AddParticleType(kplus);
     
     fGenerator = gevsim;
*/

/*   
   
    AliMevSimConfig *c = new AliMevSimConfig(1);
    c->SetRectPlane(1);                                 // reaction plane control, model 4
    c->SetGrid(80,80);

    AliGenMevSim *mevsim = new AliGenMevSim(c);
    mevsim->SetPtRange(0.001, 3);
    mevsim->SetMomentumRange(0.1, 3);
    mevsim->SetTrackingFlag(0);
    mevsim->SetOrigin(0.0, 0.0, 0.0);
    mevsim->SetSigma(0.0, 0.0, 0.0);
    AliMevSimParticle *kplus = new AliMevSimParticle(kKPlus, fNpart, 0, 0.25, 0.0, 2, 0.15, 0.0, 0.0 );
    mevsim->AddParticleType(kplus);
    fGenerator = mevsim;
*/
    
    fGenerator->SetOrigin(fOrigin[0],fOrigin[1],fOrigin[2]);
    static const Double_t kDegToRadCF = 180./TMath::Pi();
    fGenerator->SetMomentumRange(fPtMin,fPtMax);
    fGenerator->SetPhiRange(kDegToRadCF*fPhiMin,kDegToRadCF*fPhiMax);
    fGenerator->SetYRange(fYMin,fYMax);
    fGenerator->SetThetaRange(kDegToRadCF*fThetaMin,kDegToRadCF*fThetaMax);
    fGenerator->Init();
    
   }

//  fQCoarseBackground =  new TH3D("fQCoarseBackground","",fQNBins,-fQRange,fQRange, fQNBins,-fQRange,fQRange, fQNBins,-fQRange,fQRange);
//  fQCoarseSignal =  new TH3D("fQCoarseSignal","fQCoarseSignal",fQNBins,-fQRange,fQRange, fQNBins,-fQRange,fQRange, fQNBins,-fQRange,fQRange);
//  fQCoarseBackground->Sumw2();
//  fQCoarseSignal->Sumw2();
   
  fQSignal =  new TH3D("fQSignal1","fQSignal",fQNBins,-fQRange,fQRange, fQNBins,-fQRange,fQRange, fQNBins,-fQRange,fQRange);
  fQBackground =  new TH3D("fQBackground1","fQBackground",fQNBins,-fQRange,fQRange, fQNBins,-fQRange,fQRange, fQNBins,-fQRange,fQRange);

  fQSecondSignal =  new TH3D("fQSignal2","fQSignal",fQNBins,-fQRange,fQRange, fQNBins,-fQRange,fQRange, fQNBins,-fQRange,fQRange);
  fQSecondBackground =  new TH3D("fQBackground2","fQBackground",fQNBins,-fQRange,fQRange, fQNBins,-fQRange,fQRange, fQNBins,-fQRange,fQRange);

  fQSignal->Sumw2();
  fQBackground->Sumw2();
  fQSecondSignal->Sumw2();
  fQSecondBackground->Sumw2();
  
  fLogFile = new ofstream("BadEvent",ios::out);
  
}
/***********************************************************/

void AliGenHBTosl::Generate()
{
 //the main method
  
 ofstream& logfile = *fLogFile;
 logfile<<"Generate"<<"Attempts to generate "<<fNpart<<" particles.";
 
 
 if (fStackBuffer == 0x0) fStackBuffer = new TList();
 //Here is initialization level
 if (fSignalShapeCreated == kFALSE)
  {
   TH3D *hs = 0x0, *hb = 0x0;
   TFile* file;

   file = TFile::Open("QTSignal.root");
   if (file)
    {
     hs = (TH3D*)file->Get("fQSignal1");
     if (hs) hs->SetDirectory(0x0);
    }
   delete file;
   
   file = TFile::Open("QTBackground.root");
   if (file)
    {
     hb = (TH3D*)file->Get("fQBackground1");
     if (hb) hb->SetDirectory(0x0);
    }
   delete file;
   
   if (hs && hb)
    {
      Info("Generate","**********************************");
      Info("Generate","Found starting histograms in files");
      Info("Generate","**********************************");
      delete fQSignal;
      delete fQBackground;
      fQSignal = hs;
      fQBackground = hb;
    }
   else
    { 
      TH3D *cs = 0x0, *cb = 0x0;
      file = TFile::Open("QTCoarseBackground.root");
      if (file)
       {
        cb = (TH3D*)file->Get("fQCoarseBackground");
        if (cb) cb->SetDirectory(0x0);
       }
      delete file;

      file = TFile::Open("QTCoarseSignal.root");
      if (file)
       {
        cs = (TH3D*)file->Get("fQCoarseSignal");
        if (cs) cs->SetDirectory(0x0);
       }
      delete file;

      if (cs && cb)
       {

         Info("Generate","Got Coarse signal and bkg from files");
         delete fQCoarseBackground;
         delete fQCoarseSignal;
         fQCoarseSignal = cs;
         fQCoarseBackground = cb;
       }
      else
       { 
         if (cb)
           {
             Info("Generate","Got Coarse bkg from file");
             delete fQCoarseBackground;
             fQCoarseBackground = cb;
           }
         else
           {
             fQCoarseBackground =  new TH3D("fQCoarseBackground","",fQNBins,-fQRange,fQRange, fQNBins,-fQRange,fQRange, fQNBins,-fQRange,fQRange);
             fQCoarseBackground->Sumw2();

             FillCoarse();      //create coarse background - just to know the spectrum
           }

         fQCoarseSignal =  new TH3D("fQCoarseSignal","fQCoarseSignal",fQNBins,-fQRange,fQRange, fQNBins,-fQRange,fQRange, fQNBins,-fQRange,fQRange);
         fQCoarseSignal->Sumw2();
         FillCoarseSignal();//create first coarse signal by brutal multplication coarse background and required function shape
       }
       
      StartSignal();     //Initilizes the stack that is used for generation
    }  
   fSignalShapeCreated = kTRUE;
  }

 AliStack* stack = RotateStack();

 AliStack* genstack = fGenerator->GetStack();
 if (genstack == 0x0)
  {
   genstack = new AliStack(fNpart);
   fGenerator->SetStack(genstack);  
  }
 else
  {
   genstack->Reset();
  }
 
 fGenerator->Generate();
 Int_t ntr = 0;
 if ( genstack->GetNtrack() < fNpart/2)
  {
    Warning("Generate","************************************************************");
    Warning("Generate","Generator generated (%d) less particles then expected (%d).",
             stack->GetNtrack(),fNpart/2);
    Warning("Generate","************************************************************");
  }

 TH3D* work = new TH3D("work","work",fQNBins,-fQRange,fQRange,fQNBins,-fQRange,fQRange,fQNBins,-fQRange,fQRange);
 work->SetDirectory(0x0);
 work->Sumw2();
 
 Double_t*** chiarray = new Double_t** [fQNBins+1];
 Double_t*** sigarray = new Double_t** [fQNBins+1];
 
 for (Int_t i = 1; i<=fQNBins; i++)
   {
      chiarray[i] =  new Double_t* [fQNBins+1];
      sigarray[i] =  new Double_t* [fQNBins+1];
      
      for (Int_t k = 1; k<=fQNBins; k++)
       {
         chiarray[i][k] =  new Double_t [fQNBins+1];
         sigarray[i][k] =  new Double_t [fQNBins+1];
       }
   }

      
 Double_t scale = Scale(fQSignal,fQBackground);
 work->Divide(fQSignal,fQBackground,scale);
 
 Double_t binwdh = work->GetBinWidth(1)/2.;

 for (Int_t k = 1; k<=fQNBins; k++)
   {
    Double_t z = work->GetZaxis()->GetBinCenter(k);
    for (Int_t j = 1; j<=fQNBins; j++)
      {  
       Double_t y = work->GetYaxis()->GetBinCenter(j);
       for (Int_t i = 1; i<=fQNBins; i++)
         {
              sigarray[i][j][k] = fQSignal->GetBinContent(i,j,k);//store current value of signal histogram
              Double_t x = work->GetXaxis()->GetBinCenter(i);//get center value of a bin (qinv)
              Double_t v = GetQOutQSideQLongCorrTheorValue(x,y,z);//get expected value of CF in that qinv
              Double_t diff = v - work->GetBinContent(i,j,k);//store difference betweeon current value, and desired value
              chiarray[i][j][k] = diff; // no-x x is a weight to get good distribution
         }
       }
    }
 
 char msg[1000];
 logfile<<endl;
 snprintf(msg,1000, "\n");
 Int_t middlebin = fQNBins/2;
 
 for (Int_t k = middlebin-5; k < middlebin+5; k++)
   {
     Double_t tx = work->GetXaxis()->GetBinCenter(30);
     Double_t ty = work->GetYaxis()->GetBinCenter(30);
     Double_t tz = work->GetZaxis()->GetBinCenter(k);
     snprintf(msg,1000, "% 6.5f ",GetQOutQSideQLongCorrTheorValue(tx,ty,tz));
     logfile<<msg;
   }
 logfile<<endl;
  
 for (Int_t k = middlebin-5; k < middlebin+5; k++)
  {
    snprintf(msg,1000, "% 6.5f ",work->GetBinContent(30,30,k));
    logfile<<msg;
  }
 logfile<<endl;

 for (Int_t k = middlebin-5; k < middlebin+5; k++)
  {
    snprintf(msg, 1000, "% 6.5f ",chiarray[30][30][k]);
    logfile<<msg;
  }
 logfile<<endl;
 
 TParticle particle(fPID,0,-1,-1,-1,-1,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0);
 TParticle* second = &particle;

 Bool_t shortloop = kTRUE;
 Int_t sc = 0;//security check against infinite loop

 while ( (ntr+1) < fNpart)
   {
     Int_t xmax = 1; 
     Int_t ymax = 1; 
     Int_t zmax = 1; 
     Double_t qout;
     Double_t qside;
     Double_t qlong;
     
     
     Int_t loopmin;
     Int_t loopmax;
     
     if (shortloop)
      {
        shortloop = kFALSE;
        loopmax = fQNBins;
        loopmin = 1;
      }
     else 
      {
        shortloop = kTRUE;
        loopmax = fQNBins/2+fQNBins/4;
        loopmin = fQNBins/2-fQNBins/4;
      }
     
     
     for (Int_t k = loopmin; k <=loopmax; k++ )
      {
       qlong = work->GetZaxis()->GetBinCenter(k);
       for (Int_t j = loopmin; j<=loopmax; j++)
         {
          qside = work->GetYaxis()->GetBinCenter(j);
          for (Int_t i = loopmin; i<=loopmax; i++)
            {
              qout  = work->GetXaxis()->GetBinCenter(i);
              if (chiarray[xmax][ymax][zmax] < chiarray[i][j][k]) 
               {
                 xmax = i;
                 ymax = j;
                 zmax = k;
               }  
              
//              Double_t qdist = TMath::Sqrt(qout*qout + qside*qside + qlong*qlong);
              
//              Double_t fact = chiarray[i][j][k];//chiarray is chi2
//              if (fact > work->GetBinError(i,j,k))//if differece between what we want and 
//               {                                  //what we have is bigger than stat. error
//                 xmax = i;                        //we force to fill that bin
//                 ymax = j;
//                 zmax = k;
//                 break;
//               }
            }
         }
      }
     Double_t qlongc = work->GetZaxis()->GetBinCenter(zmax);
     Double_t qsidec = work->GetYaxis()->GetBinCenter(ymax);
     Double_t qoutc  = work->GetXaxis()->GetBinCenter(xmax);
     

     snprintf(msg,1000, "Generate Fill bin chi2(%d,%d,%d)=%f",xmax,ymax,zmax,chiarray[xmax][ymax][zmax]);
     logfile<<msg;
     
     qout  = gRandom->Uniform(qoutc -binwdh, qoutc +binwdh);
     qside = gRandom->Uniform(qsidec-binwdh, qsidec+binwdh);
     qlong = gRandom->Uniform(qlongc-binwdh, qlongc+binwdh);

     TParticle* first = 0;
     Int_t jj = 0;
     
     while (jj < genstack->GetNtrack())
      {
        TParticle* tmpp = genstack->Particle(jj++);
        if (tmpp->GetPdgCode() == fPID)
         {
           if (CheckParticle(tmpp,0x0,stack) == kFALSE)
            {
              first = tmpp;
              break;
            } 
         }
      } 
     
     if (first == 0x0)
      {
        if ( fDebug > 2 ) Info("StartSignal","No more particles of that type");
        break;
      }
     
     //Here put the check if this particle do not fall into signal region with other partticle
     
      Int_t retval = GetThreeD(first,second,qout,qside,qlong);
      if (retval)
       {
         //Info("StartSignal","Can not find momenta for this OSL and particle");
         continue;
       }
    //in case this particle is falling into signal area with another
    //particle we take a another pair
    //it can intruduce artificial correlations 
     Bool_t checkresult = CheckParticle(second,first,stack);
     if ( checkresult  && (sc < 10) ) 
      { 
        sc++;
        continue;
      }  
     sc = 0;
     
     //Put on output stack
     SetTrack(first,ntr);
     SetTrack(second,ntr);

     //Put on internal stack   
     Int_t etmp;
     SetTrack(first,etmp,stack);
     SetTrack(second,etmp,stack);
     
     Double_t y = GetQOutQSideQLongCorrTheorValue(qoutc,qsidec,qlongc);
      
     sigarray[xmax][ymax][zmax] ++; 
     chiarray[xmax][ymax][zmax] = scale*sigarray[xmax][ymax][zmax]/fQBackground->GetBinContent(xmax,ymax,zmax);
     chiarray[xmax][ymax][zmax] = (y - chiarray[xmax][ymax][zmax]);
      
    }
  
 Mix(fStackBuffer,fQBackground,fQSecondSignal); //upgrate background
 Mix(stack,fQSignal,fQSecondBackground); //upgrate signal
 
 delete work;
 
 for (Int_t i = 1; i<=fQNBins; i++)
   {
     for (Int_t k = 1; k<=fQNBins; k++)
      {
        delete [] chiarray[i][k]; 
        delete [] sigarray[i][k];
      }
     delete [] chiarray[i];
     delete [] sigarray[i];
   }
 delete [] chiarray;
 delete [] sigarray;
}
/***********************************************************/

void AliGenHBTosl::GetOneD(TParticle* first, TParticle* second,Double_t qinv)
{
//deprecated method that caclulates momenta of the second particle 
// out of qinv and the first particle    
//first particle is rotated that only X is non-zero

    
  Double_t m = first->GetMass();
  Double_t msqrd = m*m;
  Double_t fourmassSquered = 4.*msqrd;
  
  //Condition that R must fullfill to be possible to have qinv less smaller then randomized
//  Double_t rRange = qinv*TMath::Sqrt(qinv*qinv + fourmassSquered)/fourmassSquered;
//  Double_t r = gRandom->Uniform(rRange);

  Double_t r = gRandom->Uniform(qinv);
  Double_t phi = gRandom->Uniform(TMath::TwoPi());
  
  Double_t firstPx = first->P();//first particle is rotated that only X is non-zero  thus P==Px
  Double_t px = 2.*msqrd*firstPx + firstPx*qinv*qinv;
  Double_t temp = qinv*qinv*qinv*qinv  + fourmassSquered * (qinv*qinv - r*r );
  if (temp < 0.0)
   {
     Error("GetOneD","temp is less then 0: %f",temp);
     return;
   }
  temp = temp*(msqrd+firstPx*firstPx);
  
  px = (px - TMath::Sqrt(temp))/(2.*msqrd);
  
  Double_t py = r*TMath::Sin(phi);
  Double_t pz = r*TMath::Cos(phi);
  
  TVector3 firstpvector(first->Px(),first->Py(),first->Pz());
  TVector3 vector(px,py,pz);
  Rotate(firstpvector,vector);
 
  Double_t e = TMath::Sqrt(msqrd + vector.X()*vector.X() + vector.Y()*vector.Y() + vector.Z()*vector.Z());
  second->SetMomentum(vector.X(),vector.Y(),vector.Z(),e);
//  TParticle* f = new TParticle(first->GetPdgCode(),0,-1,-1,-1,-1, firstPx,0,0,e=TMath::Sqrt(msqrd+firstPx*firstPx),0.0,0.0,0.0,0.0);
//        TParticle(pdg, is, parent, -1, kFirstDaughter, kLastDaughter,
//                px, py, pz, e, vx, vy, vz, tof);
 
  AliDebug(1,Form("Randomized qinv = %f, obtained = %f",qinv,GetQInv(first,second)));

}
/***********************************************************/

Int_t AliGenHBTosl::GetThreeD(TParticle* first,TParticle* second, Double_t qout, Double_t qside, Double_t qlong)
{
//deprecated method that caclulates momenta of the second particle 
//out of  qout qside and qlong and the first particle    
  Double_t m = first->GetMass();
  Double_t m2 = m*m;
  
  Double_t px = first->P();//first particle is rotated that only X is non-zero  thus P==Px
  Double_t px2 = px*px;

  
  Double_t qout2 = qout*qout;
  Double_t qside2 = qside*qside;
  Double_t qlong2 = qlong*qlong;
  
  
  Double_t util1 = 4.*px2 - qside2;//4*P1x^2 - Y^2
  if (util1 < 0) 
   {
     Info("GetThreeD","4.*px2* - qside2 is negative px: %f, qside: %f",px,qside);
     return 1;
   }  
  Double_t util2 = TMath::Sqrt(px2*qout2*util1);
  

  Double_t p2x,p2y,p2z;
      
//  if ( (qside >= 0) && (qout >= 0) && (qlong >= 0))
  if (qout >= 0)
   {
     //p2x
     Double_t tmp = px*(2.*px2 - qside2);
     tmp -=  util2;
     p2x = tmp/(2.*px2);

     //p2y
     tmp =  qout - TMath::Sqrt(util1);
     p2y = - (tmp*qside)/(2.*px);
      
     //p2z
     tmp = 4.*m2 + 2.*qout2+qlong2;
     tmp *= px;
     tmp -= 2.*util2;//!!!
     tmp += 4.*px*px2;
     tmp *= qlong2;

     Double_t m2px2 = m2+px2;
     Double_t tmp2 = m2px2*tmp;

     tmp  = 4.*(m2px2+qout2) + qlong2;
     tmp *= px;
     tmp -= 4.*util2;
     tmp *= 4.*(m2px2) + qlong2;
     tmp *= qlong2*qlong2;
     tmp *= m2px2*m2px2;
     tmp *= px;
     if (tmp < 0)
      {
        Error("","Argument of sqrt is negative");
        return 1;
      }

     tmp2 += TMath::Sqrt(tmp); 

     tmp = 8.0*px*m2px2*m2px2;
     p2z = -TMath::Sqrt(tmp2/tmp);
     if (qlong < 0) p2z = -p2z;
   }
  else
   {
     //p2x
     Double_t tmp = px*(2.*px2 - qside2);
     tmp +=  util2;
     p2x = tmp/(2.*px2);

     //p2y
     tmp =  qout - TMath::Sqrt(util1);
     p2y = - (tmp*qside)/(2.*px);
      
     //p2z
     tmp = 4.*m2 + 2.*qout2+qlong2;
     tmp *= px;
     tmp += 2.*util2;//!!!
     tmp += 4.*px*px2;
     tmp *= qlong2;

     Double_t m2px2 = m2+px2;
     Double_t tmp2 = m2px2*tmp;

     tmp  = 4.*(m2px2+qout2) + qlong2;
     tmp *= px;
     tmp += 4.*util2;
     tmp *= 4.*(m2px2) + qlong2;
     tmp *= qlong2*qlong2;
     tmp *= m2px2*m2px2;
     tmp *= px;
     if (tmp < 0)
      {
        Error("","Argument of sqrt is negative");
        return 1;
      }

     tmp2 += TMath::Sqrt(tmp); 

     tmp = 8.0*px*m2px2*m2px2;
     p2z = -TMath::Sqrt(tmp2/tmp);
     if (qlong < 0) p2z = -p2z;
    }
     
//  if ( (qside >= 0) && (qout >= 0) && (qlong >= 0))  p2z = -p2z;
  
  TVector3 firstpvector(first->Px(),first->Py(),first->Pz());
  TVector3 vector(p2x,p2y,p2z);
  Rotate(firstpvector,vector);
 
  Double_t e = TMath::Sqrt(m2 + vector.X()*vector.X() + vector.Y()*vector.Y() + vector.Z()*vector.Z());
  second->SetMomentum(vector.X(),vector.Y(),vector.Z(),e);
  
////////////  
  if ( AliDebugLevel() > 3 )
   {
     e=TMath::Sqrt(m2+px*px);  
     TParticle* f = new TParticle(first->GetPdgCode(),0,-1,-1,-1,-1, px , 0.0, 0.0, e,0.0,0.0,0.0,0.0);

     e = TMath::Sqrt(m2 + p2x*p2x + p2y*p2y + p2z*p2z);
     TParticle* s = new TParticle(first->GetPdgCode(),0,-1,-1,-1,-1, p2x, p2y, p2z, e, 0.0, 0.0, 0.0, 0.0);

     Double_t qo, qs, ql;
     GetQOutQSideQLong(f,s,qo, qs, ql);

     Info("GetThreeD","TEST");
     f->Print();
     s->Print();
     Info("GetThreeD","Required %f %f %f",qout,qside,qlong);
     Info("GetThreeD","Got      %f %f %f",qo,qs,ql);
     if ( qout == qo)
        if (qside == qs)
          if  (qlong == ql)
            Info("GetThreeD","!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
   }  
///////////// 
  return 0;
}
/***********************************************************/

void AliGenHBTosl::StartSignal()
{ 
//Starts the signal histograms  
 ofstream& logfile = *fLogFile;
 logfile<<"\n\n\n\n";
 logfile<<"************************************************"<<endl;
 logfile<<"************************************************"<<endl;
 logfile<<"               StartSignal                      "<<endl;
 logfile<<"************************************************"<<endl;
 logfile<<"************************************************"<<endl;
 
 AliStack* stack;
 
 fSwapped = kFALSE;
 
 TParticle particle(fPID,0,-1,-1,-1,-1,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0);
 TParticle* second = &particle;
 
 TIter next(fStackBuffer);
 while(( stack=(AliStack*)next() ))
  {
    stack->Reset();
  }

 AliStack* genstack = fGenerator->GetStack();
 if (genstack == 0x0)
  {
   genstack = new AliStack(fNpart);
   fGenerator->SetStack(genstack);  
  }
 else
  {
   genstack->Reset();
  }
 
 StartSignalPass1();
 //We alread have detailed histograms and we do not need Coarse anymore
 delete fQCoarseSignal;
 delete fQCoarseBackground;
 fQCoarseSignal = 0x0;
 fQCoarseBackground = 0x0;

  
 const Double_t kNDF = fQNBins*fQNBins*fQNBins;
 
 TH3D* work = new TH3D("work","work",fQNBins,-fQRange,fQRange,fQNBins,-fQRange,fQRange,fQNBins,-fQRange,fQRange);
 work->Sumw2();
 work->SetDirectory(0x0);
 
 
 Double_t binwdh = work->GetBinWidth(1)/2.;
 
 Double_t*** chiarray = new Double_t** [fQNBins+1];
 Double_t*** sigarray = new Double_t** [fQNBins+1];
 
 for (Int_t i = 1; i<=fQNBins; i++)
   {
      chiarray[i] =  new Double_t* [fQNBins+1];
      sigarray[i] =  new Double_t* [fQNBins+1];
      
      for (Int_t k = 1; k<=fQNBins; k++)
       {
         chiarray[i][k] =  new Double_t [fQNBins+1];
         sigarray[i][k] =  new Double_t [fQNBins+1];
       }
   }

  
 Float_t chisqrchange = fMaxChiSquereChange + 1.;
 Float_t chisqrPerDF = fMaxChiSquerePerNDF;
 Float_t chisqrold = 0.0;
 
 Int_t counter = 1;
 Int_t niterations = 1;
 Int_t rotaxisorder = 1;//defines order of looping over 3D histo (X,Y,Z or Y,Z,X or Z,X,Y)
 
 Bool_t flag = kTRUE;
 Bool_t shortloop = kTRUE;
 TCanvas* c1 = new TCanvas();


 printf("\n");
 Info("StartSignal","\n\n\n\nSecond Pass\n\n\n\n");

 while ( ( (chisqrPerDF > fMaxChiSquereChange) || flag) && (niterations++ < fMaxIterations)  )
  {
   
   logfile<<"StartSignal\n";
   logfile<<" Row 1 Theory, 2 current value, 3 Chi2 \n";

   Double_t chisqrnew = 0.0;
   flag = kFALSE;
     
   Double_t scale = Scale(fQSignal,fQBackground);
   work->Divide(fQSignal,fQBackground,scale);
   
   if ( (counter%100) == 0) 
    {
      c1->cd();
      char buff[50];
      snprintf(buff,50, "QTWorkPass2.%3d.root",counter);
      TFile* file = TFile::Open(buff,"update");
      work->Write();
      work->SetDirectory(0x0);
      file->Close();
      delete file;

      snprintf(buff,50, "QTBackgroundPass2.%3d.root",counter);
      file = TFile::Open(buff,"update");
      fQBackground->Write();
      fQBackground->SetDirectory(0x0);
      file->Close();
      delete file;

      snprintf(buff,50, "QTSignalPass2.%3d.root",counter);
      file = TFile::Open(buff,"update");
      fQSignal->Write();
      fQSignal->SetDirectory(0x0);
      file->Close();
      delete file;
    }

   counter++;
   Int_t novertresh = 0;
   for (Int_t k = 1; k<=fQNBins; k++)
    {
      Double_t z = work->GetZaxis()->GetBinCenter(k);
      for (Int_t j = 1; j<=fQNBins; j++)
        {  
          Double_t y = work->GetYaxis()->GetBinCenter(j);
          for (Int_t i = 1; i<=fQNBins; i++)
            {
              Double_t x = work->GetXaxis()->GetBinCenter(i);//get center value of a bin (qout)
              sigarray[i][j][k] = fQSignal->GetBinContent(i,j,k);//store current value of signal histogram
              Double_t v = GetQOutQSideQLongCorrTheorValue(x,y,z);//get expected value of CF in that qinv
              Double_t diff = v - work->GetBinContent(i,j,k);//store difference betweeon current value, and desired value 
              chiarray[i][j][k] = diff; // no-x x is a weight to get good distribution
              Double_t be = work->GetBinError(i,j,k);
              chisqrnew += diff*diff/(be*be);//add up chisq

              //even if algorithm is stable (chi sqr change less then threshold)
              //and any bin differs more then 5% from expected value we continue
              Double_t fact = diff;
              if (TMath::Abs(fact) > 0.1) 
               {
                 flag = kTRUE; 
                 novertresh++;
               } 
            }
         }   
     }
    

   char msg[1000];

   printf("\n");
  
   for (Int_t k = 25; k < 36; k++)
    {
      Double_t tx = work->GetXaxis()->GetBinCenter(30);
      Double_t ty = work->GetYaxis()->GetBinCenter(30);
      Double_t tz = work->GetZaxis()->GetBinCenter(k);
      snprintf(msg,1000, "% 6.5f ",GetQOutQSideQLongCorrTheorValue(tx,ty,tz));
      logfile<<msg;
    }
   logfile<<endl;

   for (Int_t k = 25; k < 36; k++)
    {
      snprintf(msg, 1000, "%6.5f ",work->GetBinContent(30,30,k));
      logfile<<msg;
    }
   logfile<<endl;

   for (Int_t k = 25; k < 36; k++)
    {
      snprintf(msg,1000, "% 6.5f ",chiarray[30][30][k]);
      logfile<<msg;
    }
   logfile<<endl;
    
   chisqrchange = TMath::Abs(chisqrnew - chisqrold)/chisqrnew;
   chisqrold = chisqrnew;

   chisqrPerDF = chisqrnew/kNDF;
   
   Info("StartSignal","Iteration %d Chi-sq change %f%%",niterations,chisqrchange*100.0);
   Info("StartSignal","ChiSq = %f, NDF = %f, ChiSq/NDF = %f",chisqrnew, kNDF, chisqrPerDF );
   Info("StartSignal","novertresh = %d",novertresh);
   
   
   stack = RotateStack();
   genstack->Reset();
   fGenerator->Generate();
   Int_t ninputparticle = 0, ntr = 0;
   if ( genstack->GetNtrack() < fNpart/2)
    {
      Warning("StartSignal","**********************************");
      Warning("StartSignal","Generator generated (%d) less particles then expected (%d).",
               genstack->GetNtrack(),fNpart/2);
      Warning("StartSignal","**********************************");
    }
   
   Int_t sc = 0; //security check against infinite loop
   while ( (ntr+1) < fNpart)//ntr is number of track generated up to now
    {
     Int_t xmax = 1; 
     Int_t ymax = 1; 
     Int_t zmax = 1; 
     Double_t qout;
     Double_t qside;
     Double_t qlong;
     
     Int_t i,j,k;
     
     Int_t* cx = 0x0;
     Int_t* cy = 0x0;
     Int_t* cz = 0x0;
     
     switch (rotaxisorder)
      {
        case 1:
          cx = &i;
          cy = &j;
          cz = &k;
          break;
        case 2:
          cx = &j;
          cy = &k;
          cz = &i;
          break;
        case 3:
          cx = &k;
          cy = &i;
          cz = &j;
          break;
      }
     rotaxisorder++;
     if (rotaxisorder > 3) rotaxisorder = 1;
     Int_t nrange;
     
     if (shortloop)
      {
        shortloop = kFALSE;
        nrange = fQNBins;
      }
     else 
      {
        shortloop = kTRUE;
        nrange = fQNBins/4;
      }
       
//     Bool_t force = kFALSE; 
     for ( k = 1; k <=nrange;k++ )
      {
       for ( j = 1; j<=nrange; j++)
         {
          for ( i = 1; i<=nrange; i++)
            {
              if ( (chiarray[*cx][*cy][*cz]) > (chiarray[xmax][ymax][zmax]) ) 
               {
                 xmax = *cx;
                 ymax = *cy;
                 zmax = *cz;
               }  
        
//              Double_t fact = chiarray[*cx][*cy][*cz];//chiarray is chi2*qinv
//              if (fact > work->GetBinError(*cx,*cy,*cz))//if differece between what we want and 
//               {                                  //what we have is bigger than stat. error
//                                                  //we force to fill that bin
//                 qout  = work->GetXaxis()->GetBinCenter(*cx);
//                 qside = work->GetYaxis()->GetBinCenter(*cy);
//                 qlong = work->GetZaxis()->GetBinCenter(*cz);

//                 Info("StartSignal"," bin: (%d,%d,%d) loop status (%d,%d,%d) \nUsing Force: chiarray: %f \nq(o,s,l): (%f,%f,%f)  signal: %d background: %d binerror: %f",
//	  *cx,*cy,*cz,i,j,k,fact,qout,qside,qlong,
//	  (Int_t)sigarray[*cx][*cy][*cz],(Int_t)fQBackground->GetBinContent(*cx,*cy,*cz),work->GetBinError(*cx,*cy,*cz));
//                 force = kTRUE;
//                 break;
//               }
               
            }
//           if (force) break; 
          }
//         if (force) break;
        }


      qout  = work->GetXaxis()->GetBinCenter(xmax);
      qside = work->GetYaxis()->GetBinCenter(ymax);
      qlong = work->GetZaxis()->GetBinCenter(zmax);

//      Info("StartSignal"," bin: (%d,%d,%d) chiarray: %f \nq(o,s,l): (%f,%f,%f)  signal: %d background: %d binerror: %f",
//            xmax,ymax,zmax,chiarray[xmax][ymax][zmax],qout,qside,qlong,
//            (Int_t)sigarray[xmax][ymax][zmax],
//            (Int_t)fQBackground->GetBinContent(xmax,ymax,zmax),
//            work->GetBinError(xmax,ymax,zmax));
          
      qout  = gRandom->Uniform(qout-binwdh,qout+binwdh);
      qside = gRandom->Uniform(qside-binwdh,qside+binwdh);
      qlong = gRandom->Uniform(qlong-binwdh,qlong+binwdh);

      TParticle* first = 0;
      while (ninputparticle < genstack->GetNtrack())
       {
         TParticle* tmpp = genstack->Particle(ninputparticle++);
         if (tmpp->GetPdgCode() == fPID)
          {
           if (CheckParticle(tmpp,0x0,stack) == kFALSE)
            {
              first = tmpp;
              break;
            } 
          }
       } 

      if (first == 0x0)
       {
         if ( fDebug > 2 ) Info("StartSignal","No more particles of that type");
         break;
       }
      
      Int_t retval = GetThreeD(first,second,qout,qside,qlong);
      if (retval)
       {
         Info("StartSignal","Can not find momenta for this OSL and particle OSL = %f %f %f",qout,qside,qlong);
         first->Print();
         second->Print();
         
         continue;
       }
     //in case this particle is falling into signal area with another
     //particle we take a another pair
     //it can intruduce artificial correlations 
      Bool_t checkresult = CheckParticle(second,first,stack);
      if ( checkresult  && (sc < 10) ) 
       { 
         sc++;
         continue;
       }  
      sc = 0;

      //Put on output stack
      SetTrack(first,ntr,stack);
      SetTrack(second,ntr,stack);

      Double_t y = GetQOutQSideQLongCorrTheorValue(qout,qside,qlong);
      
      sigarray[xmax][ymax][zmax] ++; 
      chiarray[xmax][ymax][zmax] = scale*sigarray[xmax][ymax][zmax]/fQBackground->GetBinContent(xmax,ymax,zmax);
      chiarray[xmax][ymax][zmax] = (y - chiarray[xmax][ymax][zmax]);
      
    }
   Info("StartSignal","Mixing background...");
   Mix(fStackBuffer,fQBackground,fQSecondBackground); //upgrate background
   Info("StartSignal","Mixing signal...");
   Mix(stack,fQSignal,fQSecondSignal); //upgrate background

      
//   if ( (chisqrPerDF < 2.0) && (fSwapped == kFALSE) )
//     {
//         SwapGeneratingHistograms();
//     }
   
 }
 TFile* filef = TFile::Open("QTBackground.root","recreate");
 fQBackground->Write();
 fQBackground->SetDirectory(0x0);
 filef->Close();
 delete filef;

 filef = TFile::Open("QTSignal.root","recreate");
 fQSignal->Write();
 fQSignal->SetDirectory(0x0);
 filef->Close();
 delete filef;

 
 delete c1;
 delete work; 

 for (Int_t i = 1; i<=fQNBins; i++)
   {
     for (Int_t k = 1; k<=fQNBins; k++)
      {
        delete [] chiarray[i][k]; 
        delete [] sigarray[i][k];
      }
     delete [] chiarray[i];
     delete [] sigarray[i];
   }
 delete [] chiarray;
 delete [] sigarray;
 
}
/***********************************************************/

void AliGenHBTosl::StartSignalPass1()
{
 //This method makes first part of the initialization of working histograms
 //It randomizes qout, qside and qlong from the coarse signal histogram
 
 Bool_t flag = kTRUE;
 TParticle particle(fPID,0,-1,-1,-1,-1,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0);
 TParticle* second = &particle;
 Double_t qout;
 Double_t qside;
 Double_t qlong;
 
 Info("StartSignalPass1","\n\nFirst Pass\n\n");
 
 while (flag)
  {
    Info("StartSignalPass1","NextEvent");
    AliStack* stack = RotateStack();
    AliStack* genstack = fGenerator->GetStack();
    genstack->Reset();
    fGenerator->Generate();
    Int_t j = 0, ntr = 0;
    if ( genstack->GetNtrack() < fNpart/2)
     {
       Warning("StartSignalPass1","**********************************");
       Warning("StartSignalPass1","Generator generated (%d) less particles then expected (%d).",
                genstack->GetNtrack(),fNpart/2);
       Warning("StartSignalPass1","**********************************");
     }
   
    Int_t sc = 0;//security check against infinite loop
    while ((ntr+1)<fNpart)
     {
  
//    Info("StartSignal","Number of track on output stack: = %d", ntr);
//    Info("StartSignal","Number of track on input stack: = %d\n", j);
      
      TParticle* first = 0;
      while (j < genstack->GetNtrack())
       {
         TParticle* tmpp = genstack->Particle(j++);
         if (tmpp->GetPdgCode() == fPID)
          {
           if (CheckParticle(tmpp,0x0,stack) == kFALSE)
            {
              first = tmpp;
              break;
            }
           else
            {
              Info("StartSignalPass1","Particle did not pass the safety check 1");
              tmpp->Print();
            } 
          }
       } 

      if (first == 0x0)
       {
         if ( fDebug > 2 ) Info("StartSignalPass1","No more particles of that type");
         
         break;
       }

      SetTrack(first,ntr,stack);
      
      fQCoarseSignal->GetRandom3(qout,qside,qlong);

      Int_t retval = GetThreeD(first,second,qout,qside,qlong);
      if (retval)
       {
         //Info("StartSignal","Can not find momenta for this OSL and particle");
         continue;
       }
      //in case this particle is falling into signal area with another
      //particle we take a another pair
      //it can intruduce artificial correlations 
       Bool_t checkresult = CheckParticle(second,first,stack);
       if ( checkresult  && (sc < 10) ) 
        { 
          sc++;
          Info("StartSignalPass1","Particle did not pass the safety check 2");
          second->Print();
          continue;
        }
 
       sc = 0;
      
      SetTrack(second,ntr,stack);
     }
    
    Mix(stack,fQSignal,fQSecondSignal);
    Mix(fStackBuffer,fQBackground,fQSecondBackground);
    
    flag = kFALSE;
    
    for (Int_t k = 1; k<=fQNBins; k++)
      {
       for (j = 1; j<=fQNBins; j++)
         {  
           for (Int_t i = 1; i<=fQNBins; i++)
             {
                if ( (fQBackground->GetBinContent(i,j,k) < fMinFill) )
                  {
                    //(fQSignal->GetBinContent(i,j,k) < fMinFill) || 
	Info("StartSignalPass1","bin (%d,%d,%d): signal=%f background=%f",i,j,k,
                          fQSignal->GetBinContent(i,j,k),fQBackground->GetBinContent(i,j,k));
	flag = kTRUE;//continue while
	break;//breakes for not while
                  }
             }
            if (flag == kTRUE) break;
          }
        if (flag == kTRUE) break;
      }
    
  }//while (flag)


}
/***********************************************************/

void AliGenHBTosl::FillCoarseSignal()
{
 //Makes coarse signal by multiplying the coarse background and required function
 Info("FillCoarseSignal","START");
 for (Int_t k = 1; k <=fQNBins ;k++ )
   {
    Double_t z = fQCoarseBackground->GetZaxis()->GetBinCenter(k);
    for (Int_t j = 1; j <=fQNBins; j++)
      {
       Double_t y = fQCoarseBackground->GetYaxis()->GetBinCenter(j);
       for (Int_t i = 1; i <=fQNBins; i++)
         {
           Double_t x  = fQCoarseBackground->GetXaxis()->GetBinCenter(i);
           Double_t v = GetQOutQSideQLongCorrTheorValue(x,y,z);
           Info("FillCoarseSignal","Bin (%d,%d,%d): osl(%f,%f,%f)=%f",i,j,k,x,y,z,v);
           fQCoarseSignal->SetBinContent(i,j,k,v*fQCoarseBackground->GetBinContent(i,j,k));
         }   
       }  
    }   
 
  //if (AliDebugLevel()) 
  TestCoarseSignal();
 
  Info("FillCoarseSignal","DONE"); 
}
/***********************************************************/

void AliGenHBTosl::FillCoarse()
{
  //creates the statistical background histogram on the base of input from
  //fGenerator
  Info("FillCoarse","START");
   
  AliStack* stack;
  Int_t niter = 0;
  
  Bool_t cont;
  TH3D tmph("tmph","tmph",2,0,1,2,0,1,2,0,1);
  printf("\n");

  do 
   { 
//     if (niter > 20) break;
     
     cout<<niter++<<"  bincont "<<fQCoarseBackground->GetBinContent(30,30,28)
                  <<" "<<fQCoarseBackground->GetBinContent(30,30,29)
                  <<" "<<fQCoarseBackground->GetBinContent(30,30,30)
                  <<" "<<fQCoarseBackground->GetBinContent(30,30,31)
                  <<" "<<fQCoarseBackground->GetBinContent(30,30,32)
                  <<"\n";
     fflush(0);

     stack = RotateStack();
     fGenerator->SetStack(stack);
     fGenerator->Init();
     fGenerator->Generate();

     Mix(fStackBuffer,fQCoarseBackground,&tmph);
     
     cont = kFALSE;
     
     Info("FillCoarse","fMinFill = %d",fMinFill);
     for (Int_t k = 1; k<=fQNBins; k++)
       {
        for (Int_t j = 1; j<=fQNBins; j++)
          {  
            for (Int_t i = 1; i<=fQNBins; i++)
              {
                if ( fQCoarseBackground->GetBinContent(i,j,k) < fMinFill )
                {
                  cont = kTRUE;
                  Info("FillCoarse","bin (%d,%d,%d)=%f",i,j,k,fQCoarseBackground->GetBinContent(i,j,k));
                  break;
                }

              }
            if (cont) break;
          }
         if (cont) break;
        }
   }while(cont);
   
  printf("\n");
  fGenerator->SetStack(0x0);
  Info("FillCoarse","DONE");
  
}
/***********************************************************/
 
void AliGenHBTosl::Mix(TList* eventbuffer,TH3D* denominator,TH3D* denominator2)
{
  //Fills denominators
  //Mixes events stored in the eventbuffer and fills the background histograms
  static TStopwatch stoper;
  
  if (eventbuffer == 0x0)
   {
    Error("Mix","Buffer List is null.");
    return;
   }

  if (denominator == 0x0)
   {
    Error("Mix","Denominator histogram is null.");
    return;
   }

  if (denominator2 == 0x0)
   {
    Error("Mix","Denominator2 histogram is null.");
    return;
   }

  Info("Mix","%s",denominator->GetName());
  stoper.Start();
  
  TIter next(eventbuffer);
  AliStack* firstevent;
  AliStack* secondevent = 0x0;
  
  while(( firstevent=(AliStack*)next() ))
   {
    if (secondevent == 0x0) 
     {
       secondevent = firstevent;
       continue;
     }
//    Info("Mix","Mixing %#x with %#x",firstevent,secondevent);
    for(Int_t j = 0; j < firstevent->GetNtrack(); j++ )
     { 
       TParticle* firstpart = firstevent->Particle(j);
       
       Float_t phi = firstpart->Phi();
       if ( (phi < fSamplePhiMin) || ( phi > fSamplePhiMax)) continue;
       
//       Info("Mix","Mixing %d phi %f min %f max %f",j,phi,fSamplePhiMin,fSamplePhiMax);

       for(Int_t i = 0; i < secondevent->GetNtrack(); i++ )
         {
           TParticle* secondpart = secondevent->Particle(i);
           phi = secondpart->Phi();
           if ( (phi < fSamplePhiMin) || ( phi > fSamplePhiMax)) continue;
           
           Double_t qout;
           Double_t qside;
           Double_t qlong;
           GetQOutQSideQLong(firstpart,secondpart,qout,qside,qlong);
           denominator->Fill(qout,qside,qlong);
           denominator2->Fill(qout,qside,qlong);
         }
     }

    secondevent = firstevent;
   }
  stoper.Stop();
  stoper.Print();
  
}
/***********************************************************/

void AliGenHBTosl::Mix(AliStack* stack, TH3D* numerator, TH3D* numerator2)
{
//fils numerator with particles from stack
  static TStopwatch stoper;
  if (stack == 0x0)
   {
    Error("Mix","Stack is null.");
    return;
   }

  if ( (numerator == 0x0) || (numerator2 == 0x0) )
   {
    Error("Mix","Numerator histogram is null.");
    return;
   }

  Info("Mix","%s",numerator->GetName());
  stoper.Start();

  for(Int_t j = 0; j < stack->GetNtrack(); j++ )
    { 
      TParticle* firstpart = stack->Particle(j);
      Float_t phi = firstpart->Phi();
      if ( (phi < fSamplePhiMin) || ( phi > fSamplePhiMax)) continue;
       
      for(Int_t i = j+1; i < stack->GetNtrack(); i++ )
       {
         TParticle* secondpart = stack->Particle(i);
         phi = secondpart->Phi();
         if ( (phi < fSamplePhiMin) || ( phi > fSamplePhiMax)) continue;
         Double_t qout;
         Double_t qside;
         Double_t qlong;
         GetQOutQSideQLong(firstpart,secondpart,qout,qside,qlong);
         numerator->Fill(qout,qside,qlong);
         numerator2->Fill(qout,qside,qlong);
       }
    }
  stoper.Stop();
  stoper.Print();
  
}
/***********************************************************/

Double_t AliGenHBTosl::GetQInv(TParticle* f, TParticle* s)
{
//calculates qinv 
// cout<<f->Px()<<"   "<<s->Px()<<endl;
 Double_t pxdiff = f->Px() - s->Px();
 Double_t pydiff = f->Py() - s->Py();
 Double_t pzdiff = f->Pz() - s->Pz();
 Double_t ediff  = f->Energy() - s->Energy();
 
 Double_t qinvl = ediff*ediff - ( pxdiff*pxdiff + pydiff*pydiff + pzdiff*pzdiff );
 Double_t qinv = TMath::Sqrt(TMath::Abs(qinvl)); 
 return qinv;
}
/***********************************************************/

void  AliGenHBTosl::GetQOutQSideQLong(TParticle* f, TParticle* s,Double_t& out, Double_t& side, Double_t& lon)
{
 //returns qout,qside and qlong of the pair of particles
 out = side = lon = 10e5;
 
 Double_t pxsum = f->Px() + s->Px();
 Double_t pysum = f->Py() + s->Py();
 Double_t pzsum = f->Pz() + s->Pz();
 Double_t esum  = f->Energy() + s->Energy();
 Double_t pxdiff = f->Px() - s->Px();
 Double_t pydiff = f->Py() - s->Py();
 Double_t pzdiff = f->Pz() - s->Pz();
 Double_t ediff  = f->Energy() - s->Energy();
 Double_t kt =  0.5*TMath::Hypot(pxsum,pysum);

 Double_t k2 = pxsum*pxdiff+pysum*pydiff;
 
 if (kt == 0.0)
  {
     f->Print();
     s->Print();
     kt = 10e5;
  }
 else
  { 
    out = 0.5*k2/kt;
    side = (f->Px()*s->Py()-s->Px()*f->Py())/kt;
  }

 Double_t beta = pzsum/esum;
 Double_t gamma = 1.0/TMath::Sqrt((1.-beta)*(1.+beta));

 lon = gamma * ( pzdiff - beta*ediff );

// out = TMath::Abs(out);
// side = TMath::Abs(side);
// lon = TMath::Abs(lon);
}

/***********************************************************/

Double_t AliGenHBTosl::Scale(TH3D* num, TH3D* den)
{
 //Calculates the factor that should be used to scale
 //quatience of num and den to 1 at tail

  AliDebug(1,"Entered");
  if(!num)
   {
     AliError("No numerator");
     return 0.0;
   }
  if(!den)
   {
     AliError("No denominator");
     return 0.0;
   }

  if(fNBinsToScale < 1)
   {
    AliError("Number of bins for scaling is smaller than 1");
    return 0.0;
   }
  Int_t fNBinsToScaleX = fNBinsToScale;
  Int_t fNBinsToScaleY = fNBinsToScale;
  Int_t fNBinsToScaleZ = fNBinsToScale;

  Int_t nbinsX = num->GetNbinsX();
  if (fNBinsToScaleX > nbinsX) 
   {
    AliError("Number of X bins for scaling is bigger thnan number of bins in histograms");
    return 0.0;
   }
   
  Int_t nbinsY = num->GetNbinsX();
  if (fNBinsToScaleY > nbinsY) 
   {
    AliError("Number of Y bins for scaling is bigger thnan number of bins in histograms");
    return 0.0;
   }

  Int_t nbinsZ = num->GetNbinsZ();
  if (fNBinsToScaleZ > nbinsZ) 
   {
    AliError("Number of Z bins for scaling is bigger thnan number of bins in histograms");
    return 0.0;
   }

  AliDebug(1,"No errors detected");

  Int_t offsetX = nbinsX - fNBinsToScaleX - 1; //bin that we start loop over bins in axis X
  Int_t offsetY = nbinsY - fNBinsToScaleY - 1; //bin that we start loop over bins in axis Y
  Int_t offsetZ = nbinsZ - fNBinsToScaleZ - 1; //bin that we start loop over bins in axis Z

  Double_t densum = 0.0;
  Double_t numsum = 0.0;
  
  for (Int_t k = offsetZ; k<nbinsZ; k++)
    for (Int_t j = offsetY; j<nbinsY; j++)
      for (Int_t i = offsetX; i<nbinsX; i++)
       {
        if ( num->GetBinContent(i,j,k) > 0.0 )
         {
           
           densum += den->GetBinContent(i,j,k);
           numsum += num->GetBinContent(i,j,k);
         }
       }
  
  AliDebug(1,Form("numsum=%f densum=%f fNBinsToScaleX=%d fNBinsToScaleY=%d fNBinsToScaleZ=%d",
          numsum,densum,fNBinsToScaleX,fNBinsToScaleY,fNBinsToScaleZ));
  
  if (numsum == 0) return 0.0;
  Double_t ret = densum/numsum;

  AliDebug(1,Form("returning %f",ret));
  return ret;
   
}
/***********************************************************/

void AliGenHBTosl::TestCoarseSignal()
{
//Tests how works filling from generated histogram shape
  TH3D* work = new TH3D("work","work",fQNBins,-fQRange,fQRange,fQNBins,-fQRange,fQRange,fQNBins,-fQRange,fQRange);
     
//  for (Int_t i = 0; i < fQCoarseBackground->GetEntries() ;i++)
//   {
//     Double_t x,y,z;
//     fQCoarseSignal->GetRandom3(x,y,z);
//     work->Fill(x,y,z);
//   }

  TCanvas* c1 = new TCanvas();
  c1->cd();
  work->Draw();
  c1->SaveAs("QTwork.root");
  TFile* file = TFile::Open("QTwork.root","update");
//  work->Write();
  work->SetDirectory(0x0);
  file->Close();
  
  fQCoarseSignal->Draw(); 
  c1->SaveAs("QTCoarseSignal.root");
  file = TFile::Open("QTCoarseSignal.root","update");
  fQCoarseSignal->Write();
  fQCoarseSignal->SetDirectory(0x0);
  file->Close();
  
  fQCoarseBackground->Draw(); 
  c1->SaveAs("QTCoarseBackground.root");
  file = TFile::Open("QTCoarseBackground.root","update");
  fQCoarseBackground->Write();
  fQCoarseBackground->SetDirectory(0x0);
  file->Close();
  
  TH1 *result = (TH1*)fQCoarseBackground->Clone("ratio");
  result->SetTitle("ratio");
  Float_t normfactor = Scale(work,fQCoarseBackground);
  result->Divide(work,fQCoarseBackground,normfactor);//work
  
  
  c1->cd();
  result->Draw();
  c1->SaveAs("QTresult.root");
  file = TFile::Open("QTresult.root","update");
  result->Write();
  result->SetDirectory(0x0);
  file->Close();
  
  delete work;
  delete c1;
}
/***********************************************************/

void AliGenHBTosl::SetTrack(TParticle* p, Int_t& ntr) 
{
//Shortcut to PushTrack(bla,bla,bla,bla.............)
   if (p->P() == 0.0)
    {
      Error("SetTrack(TParticle*,Int_t&)","Particle has zero momentum");
      return;
    }
   
   
   Int_t pdg = p->GetPdgCode();
   Double_t px = p->Px();
   Double_t py = p->Py();
   Double_t pz = p->Pz();
   Double_t e  = p->Energy();
   Double_t vx = p->Vx();
   Double_t vy = p->Vy();
   Double_t vz = p->Vz();
   Double_t tof = p->T();

   TVector3 pol;
   p->GetPolarisation(pol);

   Double_t polx = pol.X();
   Double_t poly = pol.Y();
   Double_t polz = pol.Z();
   TMCProcess mech = AliGenCocktailAfterBurner::IntToMCProcess(p->GetUniqueID());
   Float_t weight = p->GetWeight();

   AliGenerator::PushTrack(fTrackIt, -1, pdg, px, py, pz, e, vx, vy, vz, tof,polx, poly, polz, mech, ntr, weight);
}
/***********************************************************/

void AliGenHBTosl::SetTrack(TParticle* p, Int_t& ntr, AliStack* stack) const 
{
//Shortcut to SetTrack(bla,bla,bla,bla.............)
   if (p->P() == 0.0)
    {
      Error("SetTrack(TParticle*,Int_t&,AliStack*)","Particle has zero momentum");
      return;
    }
 
   Int_t pdg = p->GetPdgCode();
   Double_t px = p->Px();
   Double_t py = p->Py();
   Double_t pz = p->Pz();
   Double_t e  = p->Energy();

   stack->PushTrack(fTrackIt, -1, pdg, px, py, pz, e, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, kPPrimary, ntr,1,0);
}
/***********************************************************/

void AliGenHBTosl::Rotate(TVector3& relvector, TVector3& vector)
{
//This method rotates vector about the angeles that are needed to rotate 
//relvector from postion (firstPx,0,0) to its actual positon
//In other words: To make equations easier

  static TVector3 first;
  if (AliDebugLevel()>=1) 
   {
     first.SetXYZ(relvector.x(),relvector.y(),relvector.z());
   }
  
  Double_t firstPx = TMath::Sqrt( relvector.x()*relvector.x() + 
                                  relvector.y()*relvector.y() + 
                                  relvector.z()*relvector.z()   );
	              
  Double_t rotAngleZ = -TMath::ATan2(relvector.y(),relvector.x());//calculating rot angles
  relvector.RotateZ(rotAngleZ);
  rotAngleZ = -rotAngleZ;
  Double_t rotAngleY = -TMath::ATan2(relvector.z(),relvector.x());
  
  vector.RotateY(rotAngleY);
  vector.RotateZ(rotAngleZ);
  
  if (AliDebugLevel()>5)
   {
    TVector3 test(firstPx,0.0,0.0);
    test.RotateY(rotAngleY);
    test.RotateZ(rotAngleZ);
    AliInfo(Form("Rotation test: px %f %f",first.x(),test.x()));
    AliInfo(Form("Rotation test: py %f %f",first.y(),test.y()));
    AliInfo(Form("Rotation test: pz %f %f",first.z(),test.z()));
   }
}
/***********************************************************/

Double_t AliGenHBTosl::Rotate(Double_t x,Double_t y,Double_t z)
{
//Rotates vector to base where only x - coordinate is no-zero, and returns that 

  Double_t xylength = TMath::Hypot(x,y);
  Double_t sinphi = -y/xylength;
  Double_t cosphi = x/xylength;
  
  Double_t xprime = cosphi*x - sinphi*y;
  Double_t yprime = sinphi*x + cosphi*y;
  
  TVector3 v(x,y,z);
  Double_t a1 = -TMath::ATan2(v.Y(),v.X());
  
  if (AliDebugLevel()>5)
   {
     AliInfo(Form("Xpr = %f  Ypr = %f",xprime,yprime));
     AliInfo(Form("Calc sin = %f, and %f",sinphi,TMath::Sin(a1)));
     AliInfo(Form("Calc cos = %f, and %f",cosphi,TMath::Cos(a1)));
   }

  Double_t xprimezlength = TMath::Hypot(xprime,z);
  
  Double_t sintheta = z/xprimezlength;
  Double_t costheta = xprime/xprimezlength;
  
  
  Double_t xbis = sintheta*z + costheta*(cosphi*x - sinphi*y);
  
  AliInfo(Form("Calculated rot %f, modulus %f",xbis,TMath::Sqrt(x*x+y*y+z*z)));
  return xbis;
}
/***********************************************************/

AliStack* AliGenHBTosl::RotateStack()
{ 
//swaps to next stack last goes to first and is reseted

 AliStack* stack;
 if ( fStackBuffer->GetSize() >= fBufferSize )
  {
    stack = (AliStack*)fStackBuffer->Remove(fStackBuffer->Last());
  }
 else
  {
    stack = new AliStack(fNpart);
  }
     
 fStackBuffer->AddFirst(stack);
 stack->Reset();
 return stack;
}
/***********************************************************/

Double_t AliGenHBTosl::GetQInvCorrTheorValue(Double_t qinv) const
{
//Function (deprecated)
 static const Double_t kFactorsqrd = 0.197*0.197;//squared conversion factor SI<->eV
 
 return 1.0 + 0.5*TMath::Exp(-qinv*qinv*fQRadius*fQRadius/kFactorsqrd);
}
/***********************************************************/

Double_t AliGenHBTosl::GetQOutQSideQLongCorrTheorValue(Double_t& out, Double_t& side, Double_t& lon) const
{
 //Theoretical function. Wa want to get correlation of the shape of this function
 static const Double_t kFactorsqrd = 0.197*0.197;//squared conversion factor SI<->eV
 return 1.0 + 0.7*TMath::Exp(-fQRadius*fQRadius*(out*out+side*side+lon*lon)/kFactorsqrd);
}
/***********************************************************/

Bool_t AliGenHBTosl::CheckParticle(TParticle* p, TParticle* aupair ,AliStack* stack)
{
 //Checks if a given particle is falling into signal region with any other particle
 //already existing on stack
  //PH return kFALSE;
  
 if (fSignalRegion <=0) return kFALSE;
 
 for (Int_t i = 0; i < stack->GetNtrack(); i++)
  {
    TParticle* part = stack->Particle(i);
    if (part == aupair) continue;
    Double_t qout = 10e5;
    Double_t qside= 10e5;
    Double_t qlong= 10e5;
    GetQOutQSideQLong(p,part,qout,qside,qlong);
    
    if (TMath::Abs(qout)  < fSignalRegion)
      if (TMath::Abs(qside) < fSignalRegion)
       if (TMath::Abs(qlong) < fSignalRegion) 
         return kTRUE;
  }
  return kFALSE; 
}
/***********************************************************/

void AliGenHBTosl::SwapGeneratingHistograms()
{
  //Checks if it is time to swap signal and background histograms
  //if yes it swaps them
  Int_t threshold = fMinFill;
  for (Int_t k = 1; k<=fQNBins; k++)
   {
     for (Int_t j = 1; j<=fQNBins; j++)
       {  
         for (Int_t i = 1; i<=fQNBins; i++)
           {
             if ( fQSecondBackground->GetBinContent(i,j,k) < threshold) return;
           }
       }
            
   }
  
  
  Info("SwapGeneratingHistograms","*******************************************");
  Info("SwapGeneratingHistograms","*******************************************");
  Info("SwapGeneratingHistograms","*******************************************");
  Info("SwapGeneratingHistograms","****   SWAPPING HISTOGRAMS             ****");
  Info("SwapGeneratingHistograms","*******************************************");
  Info("SwapGeneratingHistograms","*******************************************");
  Info("SwapGeneratingHistograms","*******************************************");
  
  
  TH3D* h = fQSignal;
  fQSignal = fQSecondSignal;
  fQSecondSignal = h;
  fQSecondSignal->Reset();
  fQSecondSignal->SetDirectory(0x0);
  
  h = fQBackground;
  fQBackground = fQSecondBackground;
  fQSecondBackground = h;
  fQSecondBackground->Reset();
  fQSecondBackground->SetDirectory(0x0);
  
  fSwapped = kTRUE;
  
}

AliGenHBTosl& AliGenHBTosl::operator=(const  AliGenHBTosl& rhs)
{
// Assignment operator
    rhs.Copy(*this);
    return *this;
}

void AliGenHBTosl::Copy(TObject&) const
{
    //
    // Copy 
    //
    Fatal("Copy","Not implemented!\n");
}


