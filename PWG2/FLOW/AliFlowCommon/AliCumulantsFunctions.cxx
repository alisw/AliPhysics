/*************************************************************************
* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
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

/********************************** 
 * functions and equations needed * 
 * for calculation of cumulants   *
 * and final flow estimates       *
 *                                *   
 * author: Ante Bilandzic         * 
 *          (anteb@nikhef.nl)     *
 *********************************/ 

#define AliCumulantsFunctions_cxx

#include "Riostream.h"
#include "TChain.h"
#include "TFile.h"
#include "TList.h"
#include "TParticle.h"
#include <TMatrixD.h>
#include <TVectorD.h>

#include "TProfile.h"
#include "TProfile2D.h" 
#include "TProfile3D.h"
#include "TH1.h"
#include "TH1D.h"

#include "AliFlowEventSimple.h"
#include "AliFlowTrackSimple.h"
#include "AliFlowAnalysisWithCumulants.h"
#include "AliFlowCumuConstants.h"
#include "AliFlowCommonConstants.h"
#include "AliFlowCommonHist.h"
#include "AliFlowCommonHistResults.h"
#include "AliCumulantsFunctions.h"

ClassImp(AliCumulantsFunctions)

//================================================================================================================_

AliCumulantsFunctions::AliCumulantsFunctions():  
 fIntGenFun(NULL),
 fIntGenFun4(NULL),//only for other system of Eq.
 fIntGenFun6(NULL),//only for other system of Eq.
 fIntGenFun8(NULL),//only for other system of Eq.
 fIntGenFun16(NULL),//only for other system of Eq.
 fAvMult4(NULL),//only for other system of Eq.
 fAvMult6(NULL),//only for other system of Eq.
 fAvMult8(NULL),//only for other system of Eq.
 fAvMult16(NULL),//only for other system of Eq.
 fDiffPtRPGenFunRe(NULL),
 fDiffPtRPGenFunIm(NULL),
 fPtBinRPNoOfParticles(NULL),
 fDiffEtaRPGenFunRe(NULL),
 fDiffEtaRPGenFunIm(NULL),
 fEtaBinRPNoOfParticles(NULL), 
 fDiffPtPOIGenFunRe(NULL),
 fDiffPtPOIGenFunIm(NULL),
 fPtBinPOINoOfParticles(NULL),
 fDiffEtaPOIGenFunRe(NULL),
 fDiffEtaPOIGenFunIm(NULL),
 fEtaBinPOINoOfParticles(NULL),
 fifr(NULL),
 fdfr2(NULL), 
 fdfr4(NULL), 
 fdfr6(NULL),
 fdfr8(NULL),
 fAvMult(NULL),
 fQVector(NULL),
 fAverageOfSquaredWeight(NULL),
 fchr2nd(NULL),
 fchr4th(NULL),
 fchr6th(NULL),
 fchr8th(NULL),
 fch(NULL)//common control histograms  
 /*
 fdRe0(NULL),
 fdRe1(NULL),
 fdRe2(NULL),
 fdRe3(NULL),
 fdRe4(NULL),
 fdRe5(NULL),
 fdRe6(NULL),
 fdRe7(NULL),
 fdIm0(NULL),
 fdIm1(NULL),
 fdIm2(NULL),
 fdIm3(NULL),
 fdIm4(NULL),
 fdIm5(NULL),
 fdIm6(NULL),
 fdIm7(NULL)
 */
{
 //default constructor 
}

AliCumulantsFunctions::~AliCumulantsFunctions()
{
 //destructor
}

AliCumulantsFunctions::AliCumulantsFunctions(TProfile2D *intGenFun, TProfile2D *intGenFun4, TProfile2D *intGenFun6, TProfile2D *intGenFun8, TProfile2D *intGenFun16, TProfile *avMult4, TProfile *avMult6, TProfile *avMult8, TProfile *avMult16, TProfile3D *diffPtRPGenFunRe, TProfile3D *diffPtRPGenFunIm, TProfile *ptBinRPNoOfParticles, TProfile3D *diffEtaRPGenFunRe, TProfile3D *diffEtaRPGenFunIm, TProfile *etaBinRPNoOfParticles, TProfile3D *diffPtPOIGenFunRe, TProfile3D *diffPtPOIGenFunIm, TProfile *ptBinPOINoOfParticles, TProfile3D *diffEtaPOIGenFunRe, TProfile3D *diffEtaPOIGenFunIm, TProfile *etaBinPOINoOfParticles, TH1D *ifr, TH1D *dfr2, TH1D *dfr4, TH1D *dfr6, TH1D *dfr8, TProfile *avMult, TProfile *qVector, TProfile *averageOfSquaredWeight, AliFlowCommonHistResults *chr2nd, AliFlowCommonHistResults *chr4th, AliFlowCommonHistResults *chr6th, AliFlowCommonHistResults *chr8th, AliFlowCommonHist *ch):
 fIntGenFun(intGenFun),
 fIntGenFun4(intGenFun4),//only for other system of Eq.
 fIntGenFun6(intGenFun6),//only for other system of Eq.
 fIntGenFun8(intGenFun8),//only for other system of Eq.
 fIntGenFun16(intGenFun16),//only for other system of Eq.
 fAvMult4(avMult4),//only for other system of Eq.
 fAvMult6(avMult6),//only for other system of Eq.
 fAvMult8(avMult8),//only for other system of Eq.
 fAvMult16(avMult16),//only for other system of Eq.
 fDiffPtRPGenFunRe(diffPtRPGenFunRe),
 fDiffPtRPGenFunIm(diffPtRPGenFunIm),
 fPtBinRPNoOfParticles(ptBinRPNoOfParticles),
 fDiffEtaRPGenFunRe(diffEtaRPGenFunRe),
 fDiffEtaRPGenFunIm(diffEtaRPGenFunIm),
 fEtaBinRPNoOfParticles(etaBinRPNoOfParticles),
 fDiffPtPOIGenFunRe(diffPtPOIGenFunRe),
 fDiffPtPOIGenFunIm(diffPtPOIGenFunIm),
 fPtBinPOINoOfParticles(ptBinPOINoOfParticles),
 fDiffEtaPOIGenFunRe(diffEtaPOIGenFunRe),
 fDiffEtaPOIGenFunIm(diffEtaPOIGenFunIm),
 fEtaBinPOINoOfParticles(etaBinPOINoOfParticles),
 fifr(ifr),
 fdfr2(dfr2), 
 fdfr4(dfr4), 
 fdfr6(dfr6),
 fdfr8(dfr8),
 fAvMult(avMult),
 fQVector(qVector),
 fAverageOfSquaredWeight(averageOfSquaredWeight),
 fchr2nd(chr2nd),
 fchr4th(chr4th),
 fchr6th(chr6th),
 fchr8th(chr8th),
 fch(ch)//common control histograms  
 /*
 fdRe0(dRe0),
 fdRe1(dRe1),
 fdRe2(dRe2),
 fdRe3(dRe3),
 fdRe4(dRe4),
 fdRe5(dRe5),
 fdRe6(dRe6),
 fdRe7(dRe7),
 fdIm0(dIm0),
 fdIm1(dIm1),
 fdIm2(dIm2),
 fdIm3(dIm3),
 fdIm4(dIm4),
 fdIm5(dIm5),
 fdIm6(dIm6),
 fdIm7(dIm7)
 */
{
 //custom constructor 
}
  
//================================================================================================================

void AliCumulantsFunctions::Calculate()
{
 //calculate cumulants and final integrated and differential flow estimates and store the results into output histograms
 const Int_t cQmax=AliFlowCumuConstants::GetMaster()->GetQmax();     //needed for numerics
 const Int_t cPmax=AliFlowCumuConstants::GetMaster()->GetPmax();     //needed for numerics  
 const Int_t cQmax4=AliFlowCumuConstants::GetMaster()->GetQmax4();   //needed for numerics
 const Int_t cPmax4=AliFlowCumuConstants::GetMaster()->GetPmax4();   //needed for numerics
 const Int_t cQmax6=AliFlowCumuConstants::GetMaster()->GetQmax6();   //needed for numerics
 const Int_t cPmax6=AliFlowCumuConstants::GetMaster()->GetPmax6();   //needed for numerics
 const Int_t cQmax8=AliFlowCumuConstants::GetMaster()->GetQmax8();   //needed for numerics
 const Int_t cPmax8=AliFlowCumuConstants::GetMaster()->GetPmax8();   //needed for numerics
 const Int_t cQmax16=AliFlowCumuConstants::GetMaster()->GetQmax16(); //needed for numerics
 const Int_t cPmax16=AliFlowCumuConstants::GetMaster()->GetPmax16(); //needed for numerics
 
 const Int_t cFlow=AliFlowCumuConstants::GetMaster()->GetFlow();   //integrated flow coefficient to be calculated
 const Int_t cMltpl=AliFlowCumuConstants::GetMaster()->GetMltpl(); //the multiple in p=m*n (diff. flow) 
 const Int_t cnBinsPt=100;                        //number of pt bins //to be improved
 const Int_t cnBinsEta=80;                        //number of eta bins //to be improved
 
 Double_t fR0=AliFlowCumuConstants::GetMaster()->GetR0();              //needed for numerics
 //Double_t fPtMax=AliFlowCommonConstants::GetMaster()->GetPtMax(); //maximum pt
 //Double_t fPtMin=AliFlowCommonConstants::GetMaster()->GetPtMin(); //minimum pt
 //Double_t fBinWidthPt=(fPtMax-fPtMin)/cnBinsPt;    //width of pt bin (in GeV)   
 
 Bool_t fOtherEquations=AliFlowCumuConstants::GetMaster()->GetOtherEquations();     //numerical equations for cumulants solved up to different highest order 
 
 //avarage selected multiplicity
 Double_t dAvM=0.;
 if(fAvMult)
 {
  dAvM=fAvMult->GetBinContent(1);
 }

 //number of events
 Int_t nEvents=0;
 if(fAvMult)
 {
  nEvents=(Int_t)(fAvMult->GetBinEntries(1));
 }
 
 //<Q-vector stuff>
 Double_t dAvQx=0.,dAvQy=0.,dAvQ2x=0.,dAvQ2y=0.;
 if(fQVector)
 {
  dAvQx  = fQVector->GetBinContent(1); //<Q_x>
  dAvQy  = fQVector->GetBinContent(2); //<Q_y>
  dAvQ2x = fQVector->GetBinContent(3); //<(Q_x)^2>
  dAvQ2y = fQVector->GetBinContent(4); //<(Q_y)^2>
 }
 
 //<w^2>
 Double_t dAvw2 = 1.;
 if(nEvents>0)
 { 
  dAvw2 = fAverageOfSquaredWeight->GetBinContent(1); 
  if(dAvw2 == 0) 
  {
   cout<<endl;
   cout<<"WARNING: Average squared weight is 0 in GFC. Most probably one of the histograms in <weights.root> is empty. Nothing will be calculated !!!!"<<endl;
   cout<<endl;
  }
 }  
 
 //<G[p][q]>
 TMatrixD dAvG(cPmax, cQmax); 
 Bool_t someAvGEntryIsNegative=kFALSE;
 if(fIntGenFun)
 {   
  for(Int_t p=0;p<cPmax;p++)
  {
   for(Int_t q=0;q<cQmax;q++)
   {
    dAvG(p,q)=fIntGenFun->GetBinContent(p+1,q+1);
    if(dAvG(p,q)<0.)
    {
     someAvGEntryIsNegative=kTRUE;
    } 
   }
  }  
 } 
   
 /////////////////////////////////////////////////////////////////////////////      
 //////////////////gen. function for the cumulants////////////////////////////
 /////////////////////////////////////////////////////////////////////////////
    
 TMatrixD dC(cPmax, cQmax);//C[p][q]
 if(dAvM>0 && someAvGEntryIsNegative==kFALSE)
 {
  for(Int_t p=0;p<cPmax;p++)
  {
   for(Int_t q=0;q<cQmax;q++)
   {
    dC(p,q)=1.*dAvM*(pow(dAvG(p,q),(1./dAvM))-1.); 
   }
  }
 }
    
 /////////////////////////////////////////////////////////////////////////////
 ///////avaraging the gen. function for the cumulants over azimuth////////////
 /////////////////////////////////////////////////////////////////////////////
  
 TVectorD dAvC(cPmax);//<C[p][q]>
 Double_t tempHere=0.; 
 for(Int_t p=0;p<cPmax;p++)
 {
  tempHere=0.; 
  for(Int_t q=0;q<cQmax;q++)
  {
   tempHere+=1.*dC(p,q);
  } 
  dAvC[p]=1.*tempHere/cQmax;
 }
 
 /////////////////////////////////////////////////////////////////////////////
 //////////////////////////////////final results//////////////////////////////
 /////////////////////////////////////////////////////////////////////////////
 
 TVectorD cumulant(cPmax);//array to store various order cumulants
 
 //system of eq. for the cumulants  
 cumulant[0] = (-1./(60*fR0*fR0))*((-300.)*dAvC[0]+300.*dAvC[1]-200.*dAvC[2]+75.*dAvC[3]-12.*dAvC[4]);
 cumulant[1] = (-1./(6.*pow(fR0,4.)))*(154.*dAvC[0]-214.*dAvC[1]+156.*dAvC[2]-61.*dAvC[3]+10.*dAvC[4]);
 cumulant[2] = (3./(2.*pow(fR0,6.)))*(71.*dAvC[0]-118.*dAvC[1]+98.*dAvC[2]-41.*dAvC[3]+7.*dAvC[4]);
 cumulant[3] = (-24./pow(fR0,8.))*(14.*dAvC[0]-26.*dAvC[1]+24.*dAvC[2]-11.*dAvC[3]+2.*dAvC[4]);
 cumulant[4] = (120./pow(fR0,10.))*(5.*dAvC[0]-10.*dAvC[1]+10.*dAvC[2]-5.*dAvC[3]+1.*dAvC[4]);
    
 /*
 cout<<endl;
 cout<<"*********************************"<<endl;
 cout<<"cumulants:"<<endl; 
 cout<<" c_"<<cFlow<<"{2} = "<<cumulant[0]<<endl; 
 cout<<" c_"<<cFlow<<"{4} = "<<cumulant[1]<<endl;
 cout<<" c_"<<cFlow<<"{6} = "<<cumulant[2]<<endl;
 cout<<" c_"<<cFlow<<"{8} = "<<cumulant[3]<<endl; 
 cout<<"c_"<<cFlow<<"{10} = "<<cumulant[4]<<endl;  
 cout<<endl;
 */
 
 Double_t dV2=0.,dV4=0.,dV6=0.,dV8=0.,dV10=0.;//integrated flow estimates
 
 if(cumulant[0]>=0.)
 {
  dV2=pow(cumulant[0],(1./2.));
 }
 if(cumulant[1]<=0.)
 {
  dV4=pow(-cumulant[1],(1./4.));
 }
 if(cumulant[2]>=0.)
 {
  dV6=pow((1./4.)*cumulant[2],(1./6.));
 }
 if(cumulant[3]<=0.)
 {
  dV8=pow(-(1./33.)*cumulant[3],(1./8.));
 }
 if(cumulant[4]>=0.)
 {
  dV10=pow((1./456.)*cumulant[4],(1./10.));
 }

 cout<<endl;
 cout<<"***************************************"<<endl;
 cout<<"***************************************"<<endl;
 cout<<"flow estimates from GF-cumulants:      "<<endl;
 cout<<endl;          
                
 Double_t sdQ[4]={0.};
 Double_t chiQ[4]={0.};
                            
 //v_2{2}
 if(nEvents>0 && dAvM>0 && dAvw2>0 && cumulant[0]>=0.)
 { 
  if((dAvQ2x+dAvQ2y-pow(dAvQx,2.)-pow(dAvQy,2.)-pow(pow(cumulant[0],(1./2.))*(dAvM/dAvw2),2.)>0.))       
  {
   chiQ[0]=(dAvM*dV2)/pow(dAvQ2x+dAvQ2y-pow(dAvQx,2.)-pow(dAvQy,2.)-pow(dV2*dAvM/dAvw2,2.),0.5); // to be improved,  analogously for higher orders
  } 
  if(chiQ[0])
  {  
   sdQ[0]=pow(((1./(2.*dAvM*nEvents))*((1.+2.*pow(chiQ[0],2))/(2.*pow(chiQ[0],2)))),0.5);
  }
  cout<<"   v_"<<cFlow<<"{2} = "<<dV2<<" +/- "<<sdQ[0]<<endl;
  //cout<<" v_"<<cFlow<<"{2} = "<<dV2<<" +/- "<<sdQ[0]<<", chi{2} = "<<chiQ[0]<<endl;//printing also the chi
  fifr->SetBinContent(1,dV2);
  fifr->SetBinError(1,sdQ[0]);
  //filling common histograms:
  fchr2nd->FillIntegratedFlow(dV2,sdQ[0]);
  fchr2nd->FillChi(chiQ[0]);
  
  
  //abTempDeleteMe
  fchr2nd->FillIntegratedFlowRP(dV2,sdQ[0]);
  fchr2nd->FillChiRP(chiQ[0]);
  //abTempDeleteMe
  
 } 
 else 
 {
  cout<<"   v_"<<cFlow<<"{2} = Im"<<endl; 
 }
   
 //v_2{4}   
 if(nEvents>0 && dAvM>0 && dAvw2>0 && cumulant[1]<=0.)
 {
  if((dAvQ2x+dAvQ2y-pow(dAvQx,2.)-pow(dAvQy,2.)-pow(pow(-cumulant[1],(1./4.))*(dAvM/dAvw2),2.)>0.))
  {
   chiQ[1]=(dAvM*dV4)/pow(dAvQ2x+dAvQ2y-pow(dAvQx,2.)-pow(dAvQy,2.)-pow(dV4*dAvM/dAvw2,2.),0.5);
  } 
  if(chiQ[1])
  {
   sdQ[1]=(1./(pow(2.*dAvM*nEvents,0.5)))*pow((1.+4.*pow(chiQ[1],2)+1.*pow(chiQ[1],4.)+2.*pow(chiQ[1],6.))/(2.*pow(chiQ[1],6.)),0.5);
  }
  cout<<"   v_"<<cFlow<<"{4} = "<<dV4<<" +/- "<<sdQ[1]<<endl;
  //cout<<" v_"<<cFlow<<"{4} = "<<dV4<<" +/- "<<sdQ[1]<<", chi{4} = "<<chiQ[1]<<endl;//printing also the chi
  fifr->SetBinContent(2,dV4);
  fifr->SetBinError(2,sdQ[1]);
  //filling common histograms:
  fchr4th->FillIntegratedFlow(dV4,sdQ[1]);
  fchr4th->FillChi(chiQ[1]);
  
  
  //abTempDeleteMe
  fchr4th->FillIntegratedFlowRP(dV4,sdQ[1]);
  fchr4th->FillChiRP(chiQ[1]);
  //abTempDeleteMe
  
 } 
 else 
 {
  cout<<"   v_"<<cFlow<<"{4} = Im"<<endl;  
 } 
  
 //v_2{6}
 if(nEvents>0 && dAvM>0 && dAvw2>0 && cumulant[2]>=0.)
 {
  if((dAvQ2x+dAvQ2y-pow(dAvQx,2.)-pow(dAvQy,2.)-pow(pow((1./4.)*cumulant[2],(1./6.))*(dAvM/dAvw2),2.)>0.))
  {
   chiQ[2]=(dAvM*dV6)/pow(dAvQ2x+dAvQ2y-pow(dAvQx,2.)-pow(dAvQy,2.)-pow(dV6*dAvM/dAvw2,2.),0.5);
  } 
  if(chiQ[2])
  {
   sdQ[2]=(1./(pow(2.*dAvM*nEvents,0.5)))*pow((3.+18.*pow(chiQ[2],2)+9.*pow(chiQ[2],4.)+28.*pow(chiQ[2],6.)+12.*pow(chiQ[2],8.)+24.*pow(chiQ[2],10.))/(24.*pow(chiQ[2],10.)),0.5);
  } 
  cout<<"   v_"<<cFlow<<"{6} = "<<dV6<<" +/- "<<sdQ[2]<<endl;
  //cout<<" v_"<<cFlow<<"{6} = "<<dV6<<" +/- "<<sdQ[2]<<", chi{6} = "<<chiQ[2]<<endl;//printing also the chi
  fifr->SetBinContent(3,dV6);
  fifr->SetBinError(3,sdQ[2]);
  //filling common histograms:
  fchr6th->FillIntegratedFlow(dV6,sdQ[2]);
  fchr6th->FillChi(chiQ[2]);
  
  
  //abTempDeleteMe
  fchr6th->FillIntegratedFlowRP(dV6,sdQ[2]);
  fchr6th->FillChiRP(chiQ[2]);
  //abTempDeleteMe
  
  
 }
 else
 {
  cout<<"   v_"<<cFlow<<"{6} = Im"<<endl;  
 }
  
 //v_2{8}
 if(nEvents>0 && dAvM>0 && dAvw2>0 && cumulant[3]<=0.)
 { 
  if((dAvQ2x+dAvQ2y-pow(dAvQx,2.)-pow(dAvQy,2.)-pow(pow(-(1./33.)*cumulant[3],(1./8.))*(dAvM/dAvw2),2.)>0.))
  { 
   chiQ[3]=(dAvM*dV8)/pow(dAvQ2x+dAvQ2y-pow(dAvQx,2.)-pow(dAvQy,2.)-pow(dV8*dAvM/dAvw2,2.),0.5);
  }
  if(chiQ[3])
  {
   sdQ[3]=(1./(pow(2.*dAvM*nEvents,0.5)))*pow((12.+96.*pow(chiQ[3],2.)+72.*pow(chiQ[3],4.)+304.*pow(chiQ[3],6.)+257.*pow(chiQ[3],8.)+804.*pow(chiQ[3],10.)+363.*pow(chiQ[3],12.)+726.*pow(chiQ[3],14.))/(726.*pow(chiQ[3],14.)),0.5);
  } 
  cout<<"   v_"<<cFlow<<"{8} = "<<dV8<<" +/- "<<sdQ[3]<<endl;
  //cout<<" v_"<<cFlow<<"{8} = "<<dV8<<" +/- "<<sdQ[3]<<", chi{8} = "<<chiQ[3]<<endl;//printing also the chi
  fifr->SetBinContent(4,dV8);
  fifr->SetBinError(4,sdQ[3]);
  //filling common histograms:
  fchr8th->FillIntegratedFlow(dV8,sdQ[3]);
  fchr8th->FillChi(chiQ[3]);
  
  //abTempDeleteMe
  fchr8th->FillIntegratedFlowRP(dV8,sdQ[3]);
  fchr8th->FillChiRP(chiQ[3]);  
  //abTempDeleteMe  
   
 } 
 else 
 {
  cout<<"   v_"<<cFlow<<"{8} = Im"<<endl;     
 }
  
 /*
 //v_2{10}
 if (dAvM && cumulant[4]>=0.)
 {
  cout<<"v_"<<cFlow<<"{10} = "<<dV10<<endl;
 }
 else 
 {
  cout<<"v_"<<cFlow<<"{10} = Im"<<endl; 
 }
 */
 
 cout<<endl;
 cout<<"      nEvts = "<<nEvents<<", AvM = "<<dAvM<<endl; 
 cout<<"***************************************"<<endl;
 cout<<"***************************************"<<endl;
  
 
 //===========================================================================
 //                      RP = Reaction Plane
 //===========================================================================
 
 /////////////////////////////////////////////////////////////////////////////
 ///////////////////////DIFFERENTIAL FLOW CALCULATIONS////////////////////////
 /////////////////////////////////////////////////////////////////////////////
  
  TVectorD ptXRP(cnBinsPt * cPmax * cQmax);
  TVectorD ptYRP(cnBinsPt * cPmax * cQmax);
  TVectorD ptBinRPNoOfParticles(cnBinsPt);
  
  //3D profiles (for pt)
  for(Int_t b=0;b<cnBinsPt;b++)
  {
   ptBinRPNoOfParticles[b]=fPtBinRPNoOfParticles->GetBinEntries(b+1);
     
   for(Int_t p=0;p<cPmax;p++)
   {
    for(Int_t q=0;q<cQmax;q++)
    {
     if(dAvG(p,q))
     {   
      if(fDiffPtRPGenFunRe)
      {
       ptXRP[index3d(b,p,q,cnBinsPt,cPmax)]=fDiffPtRPGenFunRe->GetBinContent(b+1,p+1,q+1)/dAvG(p,q);
      }
      if(fDiffPtRPGenFunIm)
      {  
       ptYRP[index3d(b,p,q,cnBinsPt,cPmax)]=fDiffPtRPGenFunIm->GetBinContent(b+1,p+1,q+1)/dAvG(p,q);
      } 
     } 
    }
   }   
  }   
  
  TVectorD etaXRP(cnBinsEta*cPmax*cQmax);
  TVectorD etaYRP(cnBinsEta*cPmax*cQmax);
  TVectorD etaBinRPNoOfParticles(cnBinsEta);
           
  //3D profiles (for eta)
  for(Int_t b=0;b<cnBinsEta;b++)
  {
   etaBinRPNoOfParticles[b]=fEtaBinRPNoOfParticles->GetBinEntries(b);
   for(Int_t p=0;p<cPmax;p++)
   {
    for(Int_t q=0;q<cQmax;q++)
    {
     if(dAvG(p,q))
     {   
      if(fDiffEtaRPGenFunRe)
      {
       etaXRP[index3d(b,p,q,cnBinsEta,cPmax)]=fDiffEtaRPGenFunRe->GetBinContent(b+1,p+1,q+1)/dAvG(p,q);
      }
      if(fDiffEtaRPGenFunIm)
      {  
       etaYRP[index3d(b,p,q,cnBinsEta,cPmax)]=fDiffEtaRPGenFunIm->GetBinContent(b+1,p+1,q+1)/dAvG(p,q);
      } 
     } 
    }
   }   
  }   

  
  //-------------------------------------------------------------------------------------------------------------------------------
  //final results for differential flow (in pt):
  TMatrixD ptDRP(cnBinsPt, cPmax);
  Double_t tempSumForptDRP=0.;
  for (Int_t b=0;b<cnBinsPt;b++)
  {
   for (Int_t p=0;p<cPmax;p++)
   {
    tempSumForptDRP=0.; 
    for (Int_t q=0;q<cQmax;q++)
    {
     tempSumForptDRP+=cos(cMltpl*2.*q*TMath::Pi()/cQmax)*ptXRP[index3d(b,p,q,cnBinsPt,cPmax)] + sin(cMltpl*2.*q*TMath::Pi()/cQmax)*ptYRP[index3d(b,p,q,cnBinsPt,cPmax)];
    } 
     ptDRP(b,p)=1.*(pow(fR0*pow(p+1.0,0.5),cMltpl)/cQmax)*tempSumForptDRP;
   }
  } 
  
  Double_t ptRPDiffCumulant2[cnBinsPt]={0.};
  Double_t ptRPDiffCumulant4[cnBinsPt]={0.};
  Double_t ptRPDiffCumulant6[cnBinsPt]={0.};
  Double_t ptRPDiffCumulant8[cnBinsPt]={0.};
  Double_t ptRPDiffCumulant10[cnBinsPt]={0.};
  
  for (Int_t b=0;b<cnBinsPt;b++)
  {
   ptRPDiffCumulant2[b]=(1./(fR0*fR0))*(5.*ptDRP(b,0)-5.*ptDRP(b,1)+(10./3.)*ptDRP(b,2)-(5./4.)*ptDRP(b,3)+(1./5.)*ptDRP(b,4)); 
   ptRPDiffCumulant4[b]=(1./pow(fR0,4.))*((-77./6.)*ptDRP(b,0)+(107./6.)*ptDRP(b,1)-(13./1.)*ptDRP(b,2)+(61./12.)*ptDRP(b,3)-(5./6.)*ptDRP(b,4));
   ptRPDiffCumulant6[b]=(1./pow(fR0,6.))*((71./2.)*ptDRP(b,0)-59.*ptDRP(b,1)+49.*ptDRP(b,2)-(41./2.)*ptDRP(b,3)+(7./2.)*ptDRP(b,4));
   ptRPDiffCumulant8[b]=(1./pow(fR0,8.))*(-84.*ptDRP(b,0)+156.*ptDRP(b,1)-144.*ptDRP(b,2)+66.*ptDRP(b,3)-12.*ptDRP(b,4));
   ptRPDiffCumulant10[b]=(1./pow(fR0,10.))*(120.*ptDRP(b,0)-240.*ptDRP(b,1)+240.*ptDRP(b,2)-120.*ptDRP(b,3)+24.*ptDRP(b,4));
  }
  
  //diff. flow values per pt bin:
  Double_t v2ptRP[cnBinsPt]={0.};
  Double_t v4ptRP[cnBinsPt]={0.};
  Double_t v6ptRP[cnBinsPt]={0.};
  Double_t v8ptRP[cnBinsPt]={0.};
  
  //errrors:
  Double_t sdRPDiff2pt[cnBinsPt]={0.};
  Double_t sdRPDiff4pt[cnBinsPt]={0.};
  //Double_t sdDiff6pt[cnBinsPt]={0.};//to be improved (calculation needed)
  //Double_t sdDiff8pt[cnBinsPt]={0.};//to be improved (calculation needed)

  //cout<<"number of pt bins: "<<cnBinsPt<<endl;
  //cout<<"****************************************"<<endl;
  for (Int_t b=0;b<cnBinsPt;b++){ 
    //cout<<"pt bin: "<<b*fBinWidthPt<<"-"<<(b+1)*fBinWidthPt<<" GeV"<<endl;
    
    //v'_{2/2}{2}
    if(cumulant[0]>0)
    {
      v2ptRP[b]=ptRPDiffCumulant2[b]/pow(cumulant[0],0.5);
      if (dAvQ2x+dAvQ2y-pow(dAvQx,2.)-pow(dAvQy,2.)-pow(dV2*dAvM,2.)>0. && ptBinRPNoOfParticles[b]>0.) 
      {
       if(chiQ[0]>0)
       { 
        sdRPDiff2pt[b]=pow((1./(2.*ptBinRPNoOfParticles[b]))*((1.+pow(chiQ[0],2.))/pow(chiQ[0],2.)),0.5);
       }  
       //cout<<"v'_2/2{2} = "<<v2ptRP[b]<<"%, "<<" "<<"sd{2} = "<<100.*sdRPDiff2pt[b]<<"%"<<endl;
       fdfr2->SetBinContent(b+1,v2ptRP[b]);
       fdfr2->SetBinError(b+1,sdRPDiff2pt[b]);
       //common histogram (to be removed):
       fchr2nd->FillDifferentialFlow(b+1,v2ptRP[b],sdRPDiff2pt[b]);
       // Fill common result histogram:
       if(TMath::Abs(v2ptRP[b])>1.e-44) fchr2nd->FillDifferentialFlowPtRP(b+1,v2ptRP[b],sdRPDiff2pt[b]);      
       
      } else {
         //cout<<"v'_2/2{2} = Im"<<endl;
      }
    }else{
      //cout<<"v'_2/2{2} = Im"<<endl;
    } 
    
    //v'_{2/2}{4}
    if(cumulant[1]<0)
    {
      v4ptRP[b]=-ptRPDiffCumulant4[b]/pow(-cumulant[1],.75);
      if (dAvQ2x+dAvQ2y-pow(dAvQx,2.)-pow(dAvQy,2.)-pow(dV4*dAvM,2.)>0.&&ptBinRPNoOfParticles[b]>0.) // to be improved
      {
       if(chiQ[1]>0)
       {
        sdRPDiff4pt[b]=pow((1./(2.*ptBinRPNoOfParticles[b]))*((2.+6.*pow(chiQ[1],2.)+pow(chiQ[1],4.)+pow(chiQ[1],6.))/pow(chiQ[1],6.)),0.5);
       }
       //cout<<"v'_2/2{4} = "<<v4ptRP[b]<<"%, "<<" "<<"sd{4} = "<<100.*sdRPDiff4pt[b]<<"%"<<endl;
       fdfr4->SetBinContent(b+1,v4ptRP[b]);
       fdfr4->SetBinError(b+1,sdRPDiff4pt[b]);
       //common histogram (to be removed):
       fchr4th->FillDifferentialFlow(b+1,v4ptRP[b],sdRPDiff4pt[b]);
       // Fill common result histogram:
       if(TMath::Abs(v4ptRP[b])>1.e-44) fchr4th->FillDifferentialFlowPtRP(b+1,v4ptRP[b],sdRPDiff4pt[b]);
       
      } else {
         //cout<<"v'_2/2{4} = Im"<<endl;
      } 
    }else{
      //cout<<"v'_2/2{4} = Im"<<endl;
    }  
    
    //v'_{2/2}{6}
    if(cumulant[2]>0){
      //cout<<"v'_2/2{6} = "<<100.*ptRPDiffCumulant6[b]/(4.*pow((1./4.)*cumulant[2],(5./6.)))<<"%"<<endl;
      v6ptRP[b]=ptRPDiffCumulant6[b]/(4.*pow((1./4.)*cumulant[2],(5./6.)));
      fdfr6->SetBinContent(b+1,v6ptRP[b]);
      //common histogram (to be removed):
      fchr6th->FillDifferentialFlow(b+1,v6ptRP[b],0.);
      // Fill common result histogram:
      if(TMath::Abs(v6ptRP[b])>1.e-44) fchr6th->FillDifferentialFlowPtRP(b+1,v6ptRP[b],0.);
      
    }else{
      //cout<<"v'_2/2{6} = Im"<<endl;
    }     
    
    //v'_{2/2}{8}
    if(cumulant[3]<0){
      //cout<<"v'_2/2{8} = "<<-100.*ptRPDiffCumulant8[b]/(33.*pow(-(1./33.)*cumulant[3],(7./8.)))<<"%"<<endl;
      v8ptRP[b]=-ptRPDiffCumulant8[b]/(33.*pow(-(1./33.)*cumulant[3],(7./8.))); 
      fdfr8->SetBinContent(b+1,v8ptRP[b]);
      //common histogram (to be removed):
      fchr8th->FillDifferentialFlow(b+1,v8ptRP[b],0.);
      // Fill common result histogram:
      if(TMath::Abs(v8ptRP[b])>1.e-44) fchr8th->FillDifferentialFlowPtRP(b+1,v8ptRP[b],0.);
      
    }else{
      //cout<<"v'_2/2{8} = Im"<<endl;
    }       
    //cout<<"****************************************"<<endl;
  }    
 //-------------------------------------------------------
 
 
 
  ///////////////////////////////////////////////////////////////////////////////////
 //////////////////////// INTEGRATED FLOW CALCULATIONS (RP) /////////////////////////
 ///////////////////////////////////////////////////////////////////////////////////
 
 Double_t dV2RP=0., dV4RP=0., dV6RP=0., dV8RP=0.;
 Double_t dV2RPError=0., dV4RPError=0., dV6RPError=0., dV8RPError=0.;
 Double_t dSumOfYieldRP=0.;
 for (Int_t b=1;b<cnBinsPt+1;b++)
 { 
  if(fch->GetHistPtRP())
  {  
   dSumOfYieldRP+=(fch->GetHistPtRP())->GetBinContent(b);
   if(fchr2nd->GetHistDiffFlowPtRP())
   {
    dV2RP+=((fchr2nd->GetHistDiffFlowPtRP())->GetBinContent(b))*(fch->GetHistPtRP())->GetBinContent(b);
    dV2RPError+=pow((fch->GetHistPtRP())->GetBinContent(b),2.)*pow((fchr2nd->GetHistDiffFlowPtRP())->GetBinError(b),2.);  
   }
   if(fchr4th->GetHistDiffFlowPtRP())
   {
    dV4RP+=((fchr4th->GetHistDiffFlowPtRP())->GetBinContent(b))*(fch->GetHistPtRP())->GetBinContent(b);
    dV4RPError+=pow((fch->GetHistPtRP())->GetBinContent(b),2.)*pow((fchr4th->GetHistDiffFlowPtRP())->GetBinError(b),2.);
   }
   if(fchr6th->GetHistDiffFlowPtRP())
   {
    dV6RP+=((fchr6th->GetHistDiffFlowPtRP())->GetBinContent(b))*(fch->GetHistPtRP())->GetBinContent(b);
    dV6RPError+=pow((fch->GetHistPtRP())->GetBinContent(b),2.)*pow((fchr6th->GetHistDiffFlowPtRP())->GetBinError(b),2.);
   }
   if(fchr8th->GetHistDiffFlowPtRP())
   {
    dV8RP+=((fchr8th->GetHistDiffFlowPtRP())->GetBinContent(b))*(fch->GetHistPtRP())->GetBinContent(b);
    dV8RPError+=pow((fch->GetHistPtRP())->GetBinContent(b),2.)*pow((fchr8th->GetHistDiffFlowPtRP())->GetBinError(b),2.);
   }      
  } 
 }
 
 if(dSumOfYieldRP)
 {
  dV2RP/=dSumOfYieldRP;
  dV2RPError/=(dSumOfYieldRP*dSumOfYieldRP);
  dV4RP/=dSumOfYieldRP;
  dV4RPError/=(dSumOfYieldRP*dSumOfYieldRP);
  dV6RP/=dSumOfYieldRP;
  dV6RPError/=(dSumOfYieldRP*dSumOfYieldRP);
  dV8RP/=dSumOfYieldRP;
  dV8RPError/=(dSumOfYieldRP*dSumOfYieldRP);
  fchr2nd->FillIntegratedFlowRP(dV2RP,pow(dV2RPError,0.5)); 
  fchr4th->FillIntegratedFlowRP(dV4RP,pow(dV4RPError,0.5)); 
  fchr6th->FillIntegratedFlowRP(dV6RP,pow(dV6RPError,0.5));//to be improved (errors needed) 
  fchr8th->FillIntegratedFlowRP(dV8RP,pow(dV8RPError,0.5));//to be improved (errors needed)  
 }
 
 cout<<endl;
 cout<<"***************************************"<<endl;
 cout<<"***************************************"<<endl;
 cout<<"flow estimates from GF-cumulants (RP):"<<endl;
 cout<<endl;     

 cout<<"   v_"<<cFlow<<"{2} = "<<dV2RP<<" +/- "<<pow(dV2RPError,0.5)<<endl;
 cout<<"   v_"<<cFlow<<"{4} = "<<dV4RP<<" +/- "<<pow(dV4RPError,0.5)<<endl;
 cout<<"   v_"<<cFlow<<"{6} = "<<dV6RP<<" +/- "<<pow(dV6RPError,0.5)<<endl;
 cout<<"   v_"<<cFlow<<"{8} = "<<dV8RP<<" +/- "<<pow(dV8RPError,0.5)<<endl;

 cout<<endl;
 cout<<"      nEvts = "<<(fch->GetHistMultRP())->GetEntries()<<", AvM = "<<(fch->GetHistMultRP())->GetMean()<<endl; 
 cout<<"***************************************"<<endl;
 cout<<"***************************************"<<endl;

 

 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 //-------------------------------------------------------------------------------------------------------------------------------
  //final results for differential flow (in eta):
  TMatrixD etaDRP(cnBinsEta,cPmax);
  Double_t tempSumForEtaDRP=0.;
  for (Int_t b=0;b<cnBinsEta;b++)
  {
   for (Int_t p=0;p<cPmax;p++)
   {
    tempSumForEtaDRP=0.; 
    for (Int_t q=0;q<cQmax;q++)
    {
     tempSumForEtaDRP+=cos(cMltpl*2.*q*TMath::Pi()/cQmax)*etaXRP[index3d(b,p,q,cnBinsEta,cPmax)] + sin(cMltpl*2.*q*TMath::Pi()/cQmax)*etaYRP[index3d(b,p,q,cnBinsEta,cPmax)];
    } 
     etaDRP(b,p)=1.*(pow(fR0*pow(p+1.,.5),cMltpl)/cQmax)*tempSumForEtaDRP;
   }
  } 
  
  Double_t etaRPDiffCumulant2[cnBinsEta]={0.};
  Double_t etaRPDiffCumulant4[cnBinsEta]={0.};
  Double_t etaRPDiffCumulant6[cnBinsEta]={0.};
  Double_t etaRPDiffCumulant8[cnBinsEta]={0.};
  Double_t etaRPDiffCumulant10[cnBinsEta]={0.};
  
  for (Int_t b=0;b<cnBinsEta;b++)
  {
   etaRPDiffCumulant2[b]=(1./(fR0*fR0))*(5.*etaDRP(b,0)-5.*etaDRP(b,1)+(10./3.)*etaDRP(b,2)-(5./4.)*etaDRP(b,3)+(1./5.)*etaDRP(b,4));
   etaRPDiffCumulant4[b]=(1./pow(fR0,4.))*((-77./6.)*etaDRP(b,0)+(107./6.)*etaDRP(b,1)-(13./1.)*etaDRP(b,2)+(61./12.)*etaDRP(b,3)-(5./6.)*etaDRP(b,4));
   etaRPDiffCumulant6[b]=(1./pow(fR0,6.))*((71./2.)*etaDRP(b,0)-59.*etaDRP(b,1)+49.*etaDRP(b,2)-(41./2.)*etaDRP(b,3)+(7./2.)*etaDRP(b,4));
   etaRPDiffCumulant8[b]=(1./pow(fR0,8.))*(-84.*etaDRP(b,0)+156.*etaDRP(b,1)-144.*etaDRP(b,2)+66.*etaDRP(b,3)-12.*etaDRP(b,4));
   etaRPDiffCumulant10[b]=(1./pow(fR0,10.))*(120.*etaDRP(b,0)-240.*etaDRP(b,1)+240.*etaDRP(b,2)-120.*etaDRP(b,3)+24.*etaDRP(b,4));
  }
  
  //diff. flow values per eta bin:
  Double_t v2etaRP[cnBinsEta]={0.};
  Double_t v4etaRP[cnBinsEta]={0.};
  Double_t v6etaRP[cnBinsEta]={0.};
  Double_t v8etaRP[cnBinsEta]={0.};
  
  //errrors:
  Double_t sdRPDiff2eta[cnBinsEta]={0.};
  Double_t sdRPDiff4eta[cnBinsEta]={0.};
  //Double_t sdDiff6eta[cnBinsEta]={0.};//to be improved (calculation needed)
  //Double_t sdDiff8eta[cnBinsEta]={0.};//to be improved (calculation needed)

  //cout<<"number of eta bins: "<<cnBinsEta<<endl;
  //cout<<"****************************************"<<endl;
  for (Int_t b=0;b<cnBinsEta;b++){ 
    //cout<<"eta bin: "<<b*fBinWidthPt<<"-"<<(b+1)*fBinWidthPt<<" GeV"<<endl;
    
    //v'_{2/2}{2}
    if(cumulant[0]>0)
    {
      v2etaRP[b]=etaRPDiffCumulant2[b]/pow(cumulant[0],0.5);
      if (dAvQ2x+dAvQ2y-pow(dAvQx,2.)-pow(dAvQy,2.)-pow(dV2*dAvM,2.)>0. && etaBinRPNoOfParticles[b]>0.) // to be improved
      {
       if(chiQ[0]>0)
       {
        sdRPDiff2eta[b]=pow((1./(2.*etaBinRPNoOfParticles[b]))*((1.+pow(chiQ[0],2.))/pow(chiQ[0],2.)),0.5);
       } 
       //cout<<"v'_2/2{2} = "<<v2etaRP[b]<<"%, "<<" "<<"sd{2} = "<<100.*sdDiff2eta[b]<<"%"<<endl;
       //fdfr2->SetBinContent(b+1,v2etaRP[b]);
       //fdfr2->SetBinError(b+1,sdDiff2eta[b]);
       //common histogram:
       //fchr2nd->FillDifferentialFlow(b+1,v2etaRP[b],sdDiff2eta[b])
       // Fill common result histogram:
       if(TMath::Abs(v2etaRP[b])>1.e-44) fchr2nd->FillDifferentialFlowEtaRP(b+1,v2etaRP[b],sdRPDiff2eta[b]);       
       
      } else {
         //cout<<"v'_2/2{2} = Im"<<endl;
      }
    }else{
      //cout<<"v'_2/2{2} = Im"<<endl;
    } 
    
    //v'_{2/2}{4}
    if(cumulant[1]<0)
    {
      v4etaRP[b]=-etaRPDiffCumulant4[b]/pow(-cumulant[1],0.75);
      if (dAvQ2x+dAvQ2y-pow(dAvQx,2.)-pow(dAvQy,2.)-pow(dV4*dAvM,2.)>0.&&etaBinRPNoOfParticles[b]>0.) // to be improved
      {
       if(chiQ[1])
       { 
        sdRPDiff4eta[b]=pow((1./(2.*etaBinRPNoOfParticles[b]))*((2.+6.*pow(chiQ[1],2.)+pow(chiQ[1],4.)+pow(chiQ[1],6.))/pow(chiQ[1],6.)),0.5);
       } 
       //cout<<"v'_2/2{4} = "<<v4eta[b]<<"%, "<<" "<<"sd{4} = "<<100.*sdDiff4eta[b]<<"%"<<endl;
       //fdfr4->SetBinContent(b+1,v4eta[b]);
       //fdfr4->SetBinError(b+1,sdDiff4eta[b]);
       //common histogram:
       //fchr4th->FillDifferentialFlow(b+1,v4eta[b],sdDiff4eta[b]);
       // Fill common result histogram:
       if(TMath::Abs(v4etaRP[b])>1.e-44) fchr4th->FillDifferentialFlowEtaRP(b+1,v4etaRP[b],sdRPDiff4eta[b]);
       
      } else {
         //cout<<"v'_2/2{4} = Im"<<endl;
      } 
    }else{
      //cout<<"v'_2/2{4} = Im"<<endl;
    }  
    
    //v'_{2/2}{6}
    if(cumulant[2]>0){
      //cout<<"v'_2/2{6} = "<<100.*etaRPDiffCumulant6[b]/(4.*pow((1./4.)*cumulant[2],(5./6.)))<<"%"<<endl;
      v6etaRP[b]=etaRPDiffCumulant6[b]/(4.*pow((1./4.)*cumulant[2],(5./6.)));
      //fdfr6->SetBinContent(b+1,v6eta[b]);
      //common histogram:
      //fchr6th->FillDifferentialFlow(b+1,v6eta[b],0.);
      // Fill common result histogram:
      if(TMath::Abs(v6etaRP[b])>1.e-44) fchr6th->FillDifferentialFlowEtaRP(b+1,v6etaRP[b],0.);
    
    }else{
      //cout<<"v'_2/2{6} = Im"<<endl;
    }     
    
    //v'_{2/2}{8}
    if(cumulant[3]<0){
      //cout<<"v'_2/2{8} = "<<-100.*etaRPDiffCumulant8[b]/(33.*pow(-(1./33.)*cumulant[3],(7./8.)))<<"%"<<endl;
      v8etaRP[b]=-etaRPDiffCumulant8[b]/(33.*pow(-(1./33.)*cumulant[3],(7./8.))); 
      //fdfr8->SetBinContent(b+1,v8eta[b]);
      //common histogram:
      //fchr8th->FillDifferentialFlow(b+1,v8eta[b],0.);      
      // Fill common result histogram:
      if(TMath::Abs(v8etaRP[b])>1.e-44) fchr8th->FillDifferentialFlowEtaRP(b+1,v8etaRP[b],0.);
      
    }else{
      //cout<<"v'_2/2{8} = Im"<<endl;
    }       
    //cout<<"****************************************"<<endl;
  }  
 //-------------------------------------------------------------------------------------------------------------------------------
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 //===========================================================================
 //                      POI = Particle of Interest
 //===========================================================================
 
 /////////////////////////////////////////////////////////////////////////////
 ///////////////////////DIFFERENTIAL FLOW CALCULATIONS////////////////////////
 /////////////////////////////////////////////////////////////////////////////
  
  TVectorD ptX(cnBinsPt*cPmax*cQmax);
  TVectorD ptY(cnBinsPt*cPmax*cQmax);
  TVectorD ptBinPOINoOfParticles(cnBinsPt);
  
  //3D profiles (for pt)
  for(Int_t b=0;b<cnBinsPt;b++)
  {
   ptBinPOINoOfParticles[b]=fPtBinPOINoOfParticles->GetBinEntries(b+1);
   for(Int_t p=0;p<cPmax;p++)
   {
    for(Int_t q=0;q<cQmax;q++)
    {
     if(dAvG(p,q))
     {   
      if(fDiffPtPOIGenFunRe)
      {
       ptX[index3d(b,p,q,cnBinsPt,cPmax)]=fDiffPtPOIGenFunRe->GetBinContent(b+1,p+1,q+1)/dAvG(p,q);
      }
      if(fDiffPtPOIGenFunIm)
      {  
       ptY[index3d(b,p,q,cnBinsPt,cPmax)]=fDiffPtPOIGenFunIm->GetBinContent(b+1,p+1,q+1)/dAvG(p,q);
      } 
     } 
    }
   }   
  }   
  
  TVectorD etaX(cnBinsEta*cPmax*cQmax);
  TVectorD etaY(cnBinsEta*cPmax*cQmax);
  TVectorD etaBinPOINoOfParticles(cnBinsEta);
           
  //3D profiles (for eta)
  for(Int_t b=0;b<cnBinsEta;b++)
  {
   etaBinPOINoOfParticles[b]=fEtaBinPOINoOfParticles->GetBinEntries(b+1);
   for(Int_t p=0;p<cPmax;p++)
   {
    for(Int_t q=0;q<cQmax;q++)
    {
     if(dAvG(p,q))
     {   
      if(fDiffEtaPOIGenFunRe)
      {
       etaX[index3d(b,p,q,cnBinsEta,cPmax)]=fDiffEtaPOIGenFunRe->GetBinContent(b+1,p+1,q+1)/dAvG(p,q);
      }
      if(fDiffEtaPOIGenFunIm)
      {  
       etaY[index3d(b,p,q,cnBinsEta,cPmax)]=fDiffEtaPOIGenFunIm->GetBinContent(b+1,p+1,q+1)/dAvG(p,q);
      } 
     } 
    }
   }   
  }   

  /*
  if(dAvM){
  for(Int_t b=0;b<cnBinsPt;b++){cout<<"***************************************"<<endl;
 cout<<"***************************************"<<endl;
    //for(Int_t p=0;p<cPmax;p++){
      for(Int_t q=0;q<cQmax;q++){
	X[b][0][q]=fdRe0->GetBinContent(b+1,q+1)/AvG[0][q];
	Y[b][0][q]=fdIm0->GetBinContent(b+1,q+1)/AvG[0][q];
	//--------------------------------------------------
	X[b][1][q]=fdRe1->GetBinContent(b+1,q+1)/AvG[1][q];
	Y[b][1][q]=fdIm1->GetBinContent(b+1,q+1)/AvG[1][q];
	//--------------------------------------------------
	X[b][2][q]=fdRe2->GetBinContent(b+1,q+1)/AvG[2][q];
	Y[b][2][q]=fdIm2->GetBinContent(b+1,q+1)/AvG[2][q];
	//--------------------------------------------------
	X[b][3][q]=fdRe3->GetBinContent(b+1,q+1)/AvG[3][q];
	Y[b][3][q]=fdIm3->GetBinContent(b+1,q+1)/AvG[3][q];
	//--------------------------------------------------
	X[b][4][q]=fdRe4->GetBinContent(b+1,q+1)/AvG[4][q];
	Y[b][4][q]=fdIm4->GetBinContent(b+1,q+1)/AvG[4][q];
	//--------------------------------------------------
	X[b][5][q]=fdRe5->GetBinContent(b+1,q+1)/AvG[5][q];
	Y[b][5][q]=fdIm5->GetBinContent(b+1,q+1)/AvG[5][q];
	//--------------------------------------------------
	X[b][6][q]=fdRe6->GetBinContent(b+1,q+1)/AvG[6][q];
	Y[b][6][q]=fdIm6->GetBinContent(b+1,q+1)/AvG[6][q];
	//--------------------------------------------------
	X[b][7][q]=fdRe7->GetBinContent(b+1,q+1)/AvG[7][q];
	Y[b][7][q]=fdIm7->GetBinContent(b+1,q+1)/AvG[7][q];
      }
    //}   
  }
  }
  */
  
  //-------------------------------------------------------------------------------------------------------------------------------
  //final results for differential flow (in pt):
  TMatrixD ptD(cnBinsPt,cPmax);
  Double_t tempSumForPtD=0.;
  for (Int_t b=0;b<cnBinsPt;b++)
  {
   for (Int_t p=0;p<cPmax;p++)
   {
    tempSumForPtD=0.; 
    for (Int_t q=0;q<cQmax;q++)
    {
     tempSumForPtD+=cos(cMltpl*2.*q*TMath::Pi()/cQmax)*ptX[index3d(b,p,q,cnBinsPt,cPmax)] + sin(cMltpl*2.*q*TMath::Pi()/cQmax)*ptY[index3d(b,p,q,cnBinsPt,cPmax)];
    } 
     ptD(b,p)=1.*(pow(fR0*pow(p+1.0,0.5),cMltpl)/cQmax)*tempSumForPtD;
   }
  } 
  
  Double_t ptDiffCumulant2[cnBinsPt]={0.};
  Double_t ptDiffCumulant4[cnBinsPt]={0.};
  Double_t ptDiffCumulant6[cnBinsPt]={0.};
  Double_t ptDiffCumulant8[cnBinsPt]={0.};
  Double_t ptDiffCumulant10[cnBinsPt]={0.};
  
  for (Int_t b=0;b<cnBinsPt;b++)
  {
   ptDiffCumulant2[b]=(1./(fR0*fR0))*(5.*ptD(b,0)-5.*ptD(b,1)+(10./3.)*ptD(b,2)-(5./4.)*ptD(b,3)+(1./5.)*ptD(b,4));
   ptDiffCumulant4[b]=(1./pow(fR0,4.))*((-77./6.)*ptD(b,0)+(107./6.)*ptD(b,1)-(13./1.)*ptD(b,2)+(61./12.)*ptD(b,3)-(5./6.)*ptD(b,4));
   ptDiffCumulant6[b]=(1./pow(fR0,6.))*((71./2.)*ptD(b,0)-59.*ptD(b,1)+49.*ptD(b,2)-(41./2.)*ptD(b,3)+(7./2.)*ptD(b,4));
   ptDiffCumulant8[b]=(1./pow(fR0,8.))*(-84.*ptD(b,0)+156.*ptD(b,1)-144.*ptD(b,2)+66.*ptD(b,3)-12.*ptD(b,4));
   ptDiffCumulant10[b]=(1./pow(fR0,10.))*(120.*ptD(b,0)-240.*ptD(b,1)+240.*ptD(b,2)-120.*ptD(b,3)+24.*ptD(b,4));
  }
  
  //diff. flow values per pt bin:
  Double_t v2pt[cnBinsPt]={0.};
  Double_t v4pt[cnBinsPt]={0.};
  Double_t v6pt[cnBinsPt]={0.};
  Double_t v8pt[cnBinsPt]={0.};
  
  //errrors:
  Double_t sdDiff2pt[cnBinsPt]={0.};
  Double_t sdDiff4pt[cnBinsPt]={0.};
  //Double_t sdDiff6pt[cnBinsPt]={0.};//to be improved (calculation needed)
  //Double_t sdDiff8pt[cnBinsPt]={0.};//to be improved (calculation needed)

  //cout<<"number of pt bins: "<<cnBinsPt<<endl;
  //cout<<"****************************************"<<endl;
   
    
    for (Int_t b=0;b<cnBinsPt;b++){ 
    //cout<<"pt bin: "<<b*fBinWidthPt<<"-"<<(b+1)*fBinWidthPt<<" GeV"<<endl;
    
    
    //v'_{2/2}{2}
    if(cumulant[0]>0)
    {
      v2pt[b]=ptDiffCumulant2[b]/pow(cumulant[0],0.5);
      if (dAvQ2x+dAvQ2y-pow(dAvQx,2.)-pow(dAvQy,2.)-pow(dV2*dAvM,2.)>0.&&ptBinPOINoOfParticles[b]>0.)
      {
       if(chiQ[0]>0)
       { 
        sdDiff2pt[b]=pow((1./(2.*ptBinPOINoOfParticles[b]))*((1.+pow(chiQ[0],2.))/pow(chiQ[0],2.)),0.5);
       } 
       //cout<<"v'_2/2{2} = "<<v2pt[b]<<"%, "<<" "<<"sd{2} = "<<100.*sdDiff2pt[b]<<"%"<<endl;
       fdfr2->SetBinContent(b+1,v2pt[b]);
       fdfr2->SetBinError(b+1,sdDiff2pt[b]);
       //common histogram:
       fchr2nd->FillDifferentialFlow(b+1,v2pt[b],sdDiff2pt[b]);       
       // Fill common result histogram:
       if(TMath::Abs(v2pt[b])>1.e-44) fchr2nd->FillDifferentialFlowPtPOI(b+1,v2pt[b],sdDiff2pt[b]);        
      
      } else {
        //cout<<"v'_2/2{2} = Im"<<endl;
      }
    }else{
      //cout<<"v'_2/2{2} = Im"<<endl;
    } 
    
    //v'_{2/2}{4}
    if(cumulant[1]<0)
    {
      v4pt[b]=-ptDiffCumulant4[b]/pow(-cumulant[1],.75);
      if (dAvQ2x+dAvQ2y-pow(dAvQx,2.)-pow(dAvQy,2.)-pow(dV4*dAvM,2.)>0.&&ptBinPOINoOfParticles[b]>0.)
      {
       if(chiQ[1]>0)
       {
        sdDiff4pt[b]=pow((1./(2.*ptBinPOINoOfParticles[b]))*((2.+6.*pow(chiQ[1],2.)+pow(chiQ[1],4.)+pow(chiQ[1],6.))/pow(chiQ[1],6.)),0.5);
       } 
       //cout<<"v'_2/2{4} = "<<v4pt[b]<<"%, "<<" "<<"sd{4} = "<<100.*sdDiff4pt[b]<<"%"<<endl;
       fdfr4->SetBinContent(b+1,v4pt[b]);
       fdfr4->SetBinError(b+1,sdDiff4pt[b]);
       //common histogram:
       fchr4th->FillDifferentialFlow(b+1,v4pt[b],sdDiff4pt[b]);
       // Fill common result histogram:
       if(TMath::Abs(v4pt[b])>1.e-44) fchr4th->FillDifferentialFlowPtPOI(b+1,v4pt[b],sdDiff4pt[b]);
          
      } else {
         //cout<<"v'_2/2{4} = Im"<<endl;
      } 
    }else{
      //cout<<"v'_2/2{4} = Im"<<endl;
    }  
    
    //v'_{2/2}{6}
    if(cumulant[2]>0){
      //cout<<"v'_2/2{6} = "<<100.*ptDiffCumulant6[b]/(4.*pow((1./4.)*cumulant[2],(5./6.)))<<"%"<<endl;
      v6pt[b]=ptDiffCumulant6[b]/(4.*pow((1./4.)*cumulant[2],(5./6.)));
      fdfr6->SetBinContent(b+1,v6pt[b]);
      //common histogram:
      fchr6th->FillDifferentialFlow(b+1,v6pt[b],0.);
      // Fill common result histogram:
      if(TMath::Abs(v6pt[b])>1.e-44) fchr6th->FillDifferentialFlowPtPOI(b+1,v6pt[b],0.);
      
    }else{
      //cout<<"v'_2/2{6} = Im"<<endl;
    }     
    
    //v'_{2/2}{8}
    if(cumulant[3]<0){
      //cout<<"v'_2/2{8} = "<<-100.*ptDiffCumulant8[b]/(33.*pow(-(1./33.)*cumulant[3],(7./8.)))<<"%"<<endl;
      v8pt[b]=-ptDiffCumulant8[b]/(33.*pow(-(1./33.)*cumulant[3],(7./8.))); 
      fdfr8->SetBinContent(b+1,v8pt[b]);
      //common histogram:
      fchr8th->FillDifferentialFlow(b+1,v8pt[b],0.);
      // Fill common result histogram:
       if(TMath::Abs(v8pt[b])>1.e-44) fchr8th->FillDifferentialFlowPtPOI(b+1,v8pt[b],0.);
      
    }else{
      //cout<<"v'_2/2{8} = Im"<<endl;
    }       
    //cout<<"****************************************"<<endl;
  }    
 //-------------------------------------------------------------------------------------------------------------------------------
 
  
   
     
  
 ///////////////////////////////////////////////////////////////////////////////////
 //////////////////////// INTEGRATED FLOW CALCULATIONS (POI) ///////////////////////
 ///////////////////////////////////////////////////////////////////////////////////
 
 Double_t dV2POI=0., dV4POI=0., dV6POI=0., dV8POI=0.;
 Double_t dV2POIError=0., dV4POIError=0., dV6POIError=0., dV8POIError=0.;
 Double_t dSumOfYieldPOI=0.;
 for (Int_t b=1;b<cnBinsPt+1;b++)
 { 
  if(fch->GetHistPtPOI())
  {  
   dSumOfYieldPOI+=(fch->GetHistPtPOI())->GetBinContent(b);
   if(fchr2nd->GetHistDiffFlowPtPOI())
   {
    dV2POI+=((fchr2nd->GetHistDiffFlowPtPOI())->GetBinContent(b))*(fch->GetHistPtPOI())->GetBinContent(b);
    dV2POIError+=pow((fch->GetHistPtPOI())->GetBinContent(b),2.)*pow((fchr2nd->GetHistDiffFlowPtPOI())->GetBinError(b),2.);  
   }
   if(fchr4th->GetHistDiffFlowPtPOI())
   {
    dV4POI+=((fchr4th->GetHistDiffFlowPtPOI())->GetBinContent(b))*(fch->GetHistPtPOI())->GetBinContent(b);
    dV4POIError+=pow((fch->GetHistPtPOI())->GetBinContent(b),2.)*pow((fchr4th->GetHistDiffFlowPtPOI())->GetBinError(b),2.);
   }
   if(fchr6th->GetHistDiffFlowPtPOI())
   {
    dV6POI+=((fchr6th->GetHistDiffFlowPtPOI())->GetBinContent(b))*(fch->GetHistPtPOI())->GetBinContent(b);
    dV6POIError+=pow((fch->GetHistPtPOI())->GetBinContent(b),2.)*pow((fchr6th->GetHistDiffFlowPtPOI())->GetBinError(b),2.);
   }
   if(fchr8th->GetHistDiffFlowPtPOI())
   {
    dV8POI+=((fchr8th->GetHistDiffFlowPtPOI())->GetBinContent(b))*(fch->GetHistPtPOI())->GetBinContent(b);
    dV8POIError+=pow((fch->GetHistPtPOI())->GetBinContent(b),2.)*pow((fchr8th->GetHistDiffFlowPtPOI())->GetBinError(b),2.);
   }      
  } 
 }
 
 if(dSumOfYieldPOI)
 {
  dV2POI/=dSumOfYieldPOI;
  dV2POIError/=(dSumOfYieldPOI*dSumOfYieldPOI);
  dV4POI/=dSumOfYieldPOI;
  dV4POIError/=(dSumOfYieldPOI*dSumOfYieldPOI);
  dV6POI/=dSumOfYieldPOI;
  dV6POIError/=(dSumOfYieldPOI*dSumOfYieldPOI);
  dV8POI/=dSumOfYieldPOI;
  dV8POIError/=(dSumOfYieldPOI*dSumOfYieldPOI);
  fchr2nd->FillIntegratedFlowPOI(dV2POI,pow(dV2POIError,0.5)); 
  fchr4th->FillIntegratedFlowPOI(dV4POI,pow(dV4POIError,0.5)); 
  fchr6th->FillIntegratedFlowPOI(dV6POI,pow(dV6POIError,0.5));//to be improved (errors needed) 
  fchr8th->FillIntegratedFlowPOI(dV8POI,pow(dV8POIError,0.5));//to be improved (errors needed)  
 }
 
 cout<<endl;
 cout<<"***************************************"<<endl;
 cout<<"***************************************"<<endl;
 cout<<"flow estimates from GF-cumulants (POI):"<<endl;
 cout<<endl;     

 cout<<"   v_"<<cFlow<<"{2} = "<<dV2POI<<" +/- "<<pow(dV2POIError,0.5)<<endl;
 cout<<"   v_"<<cFlow<<"{4} = "<<dV4POI<<" +/- "<<pow(dV4POIError,0.5)<<endl;
 cout<<"   v_"<<cFlow<<"{6} = "<<dV6POI<<" +/- "<<pow(dV6POIError,0.5)<<endl;
 cout<<"   v_"<<cFlow<<"{8} = "<<dV8POI<<" +/- "<<pow(dV8POIError,0.5)<<endl;

 cout<<endl;
 cout<<"      nEvts = "<<(fch->GetHistMultPOI())->GetEntries()<<", AvM = "<<(fch->GetHistMultPOI())->GetMean()<<endl; 
 cout<<"***************************************"<<endl;
 cout<<"***************************************"<<endl;
 cout<<endl;
 
  //-------------------------------------------------------------------------------------------------------------------------------
  //final results for differential flow (in eta):
  TMatrixD etaD(cnBinsEta, cPmax);
  Double_t tempSumForEtaD=0.;
  for (Int_t b=0;b<cnBinsEta;b++)
  {
   for (Int_t p=0;p<cPmax;p++)
   {
    tempSumForEtaD=0.; 
    for (Int_t q=0;q<cQmax;q++)
    {
     tempSumForEtaD+=cos(cMltpl*2.*q*TMath::Pi()/cQmax)*etaX[index3d(b,p,q,cnBinsEta,cPmax)] + sin(cMltpl*2.*q*TMath::Pi()/cQmax)*etaY[index3d(b,p,q,cnBinsEta,cPmax)];
    } 
     etaD(b,p)=1.*(pow(fR0*pow(p+1.,.5),cMltpl)/cQmax)*tempSumForEtaD;
   }
  } 
  
  Double_t etaDiffCumulant2[cnBinsEta]={0.};
  Double_t etaDiffCumulant4[cnBinsEta]={0.};
  Double_t etaDiffCumulant6[cnBinsEta]={0.};
  Double_t etaDiffCumulant8[cnBinsEta]={0.};
  Double_t etaDiffCumulant10[cnBinsEta]={0.};
  
  for (Int_t b=0;b<cnBinsEta;b++)
  {
   etaDiffCumulant2[b]=(1./(fR0*fR0))*(5.*etaD(b,0)-5.*etaD(b,1)+(10./3.)*etaD(b,2)-(5./4.)*etaD(b,3)+(1./5.)*etaD(b,4));
   etaDiffCumulant4[b]=(1./pow(fR0,4.))*((-77./6.)*etaD(b,0)+(107./6.)*etaD(b,1)-(13./1.)*etaD(b,2)+(61./12.)*etaD(b,3)-(5./6.)*etaD(b,4));
   etaDiffCumulant6[b]=(1./pow(fR0,6.))*((71./2.)*etaD(b,0)-59.*etaD(b,1)+49.*etaD(b,2)-(41./2.)*etaD(b,3)+(7./2.)*etaD(b,4));
   etaDiffCumulant8[b]=(1./pow(fR0,8.))*(-84.*etaD(b,0)+156.*etaD(b,1)-144.*etaD(b,2)+66.*etaD(b,3)-12.*etaD(b,4));
   etaDiffCumulant10[b]=(1./pow(fR0,10.))*(120.*etaD(b,0)-240.*etaD(b,1)+240.*etaD(b,2)-120.*etaD(b,3)+24.*etaD(b,4));
  }
  
  //diff. flow values per eta bin:
  Double_t v2eta[cnBinsEta]={0.};
  Double_t v4eta[cnBinsEta]={0.};
  Double_t v6eta[cnBinsEta]={0.};
  Double_t v8eta[cnBinsEta]={0.};
  
  //errrors:
  Double_t sdDiff2eta[cnBinsEta]={0.};
  Double_t sdDiff4eta[cnBinsEta]={0.};
  //Double_t sdDiff6eta[cnBinsEta]={0.};//to be improved (calculation needed)
  //Double_t sdDiff8eta[cnBinsEta]={0.};//to be improved (calculation needed)

  //cout<<"number of eta bins: "<<cnBinsEta<<endl;
  //cout<<"****************************************"<<endl;
  for (Int_t b=0;b<cnBinsEta;b++){ 
    //cout<<"eta bin: "<<b*fBinWidthPt<<"-"<<(b+1)*fBinWidthPt<<" GeV"<<endl;
    
    //v'_{2/2}{2}
    if(cumulant[0]>0)
    {
      v2eta[b]=etaDiffCumulant2[b]/pow(cumulant[0],0.5);
      if (dAvQ2x+dAvQ2y-pow(dAvQx,2.)-pow(dAvQy,2.)-pow(dV2*dAvM,2.)>0.&&etaBinPOINoOfParticles[b]>0.) // to be improved
      {
       if(chiQ[0]>0)
       {
        sdDiff2eta[b]=pow((1./(2.*etaBinPOINoOfParticles[b]))*((1.+pow(chiQ[0],2.))/pow(chiQ[0],2.)),0.5);
       }
       //cout<<"v'_2/2{2} = "<<v2eta[b]<<"%, "<<" "<<"sd{2} = "<<100.*sdDiff2eta[b]<<"%"<<endl;
       fdfr2->SetBinContent(b+1,v2eta[b]);
       fdfr2->SetBinError(b+1,sdDiff2eta[b]);
       //common histogram:
       //fchr2nd->FillDifferentialFlow(b+1,v2eta[b],sdDiff2eta[b]);        
       // Fill common result histogram:
       if(TMath::Abs(v2eta[b])>1.e-44) fchr2nd->FillDifferentialFlowEtaPOI(b+1,v2eta[b],sdDiff2eta[b]);       
       
      } else {
         //cout<<"v'_2/2{2} = Im"<<endl;
      }
    }else{
      //cout<<"v'_2/2{2} = Im"<<endl;
    } 
    
    //v'_{2/2}{4}
    if(cumulant[1]<0)
    {
      v4eta[b]=-etaDiffCumulant4[b]/pow(-cumulant[1],0.75);
      if (dAvQ2x+dAvQ2y-pow(dAvQx,2.)-pow(dAvQy,2.)-pow(dV4*dAvM,2.)>0.&&etaBinPOINoOfParticles[b]>0.) // to be improved
      {
       if(chiQ[1]>0)
       {
        sdDiff4eta[b]=pow((1./(2.*etaBinPOINoOfParticles[b]))*((2.+6.*pow(chiQ[1],2.)+pow(chiQ[1],4.)+pow(chiQ[1],6.))/pow(chiQ[1],6.)),0.5);
       }
       //cout<<"v'_2/2{4} = "<<v4eta[b]<<"%, "<<" "<<"sd{4} = "<<100.*sdDiff4eta[b]<<"%"<<endl;
       fdfr4->SetBinContent(b+1,v4eta[b]);
       fdfr4->SetBinError(b+1,sdDiff4eta[b]);
       //common histogram:
       //fchr4th->FillDifferentialFlow(b+1,v4eta[b],sdDiff4eta[b]);
       // Fill common result histogram:
       if(TMath::Abs(v4eta[b])>1.e-44) fchr4th->FillDifferentialFlowEtaPOI(b+1,v4eta[b],sdDiff4eta[b]);
       
      } else {
         //cout<<"v'_2/2{4} = Im"<<endl;
      } 
    }else{
      //cout<<"v'_2/2{4} = Im"<<endl;
    }  
    
    //v'_{2/2}{6}
    if(cumulant[2]>0){
      //cout<<"v'_2/2{6} = "<<100.*etaDiffCumulant6[b]/(4.*pow((1./4.)*cumulant[2],(5./6.)))<<"%"<<endl;
      v6eta[b]=etaDiffCumulant6[b]/(4.*pow((1./4.)*cumulant[2],(5./6.)));
      //fdfr6->SetBinContent(b+1,v6eta[b]);
      //common histogram:
      //fchr6th->FillDifferentialFlow(b+1,v6eta[b],0.);
      // Fill common result histogram:
      if(TMath::Abs(v6eta[b])>1.e-44) fchr6th->FillDifferentialFlowEtaPOI(b+1,v6eta[b],0.);
      
    }else{
      //cout<<"v'_2/2{6} = Im"<<endl;
    }     
    
    //v'_{2/2}{8}
    if(cumulant[3]<0){
      //cout<<"v'_2/2{8} = "<<-100.*etaDiffCumulant8[b]/(33.*pow(-(1./33.)*cumulant[3],(7./8.)))<<"%"<<endl;
      v8eta[b]=-etaDiffCumulant8[b]/(33.*pow(-(1./33.)*cumulant[3],(7./8.))); 
      //fdfr8->SetBinContent(b+1,v8eta[b]);
      //common histogram:
      //fchr8th->FillDifferentialFlow(b+1,v8eta[b],0.);
      // Fill common result histogram:
      if(TMath::Abs(v8eta[b])>1.e-44) fchr8th->FillDifferentialFlowEtaPOI(b+1,v8eta[b],0.);
      
    }else{
      //cout<<"v'_2/2{8} = Im"<<endl;
    }       
    //cout<<"****************************************"<<endl;
  }  
 //-------------------------------------------------------------------------------------------------------------------------------
 
 
 
 
 
 
 
 
 //off the record: numerical equations for cumulants solved up to different highest order  
 if(fOtherEquations)
 {
 
 //==============================================================================================================================================
  
 //up to 4th order
 //avarage selected multiplicity
 Double_t dAvM4=0.;
 if(fAvMult4)
 {
  dAvM4=fAvMult4->GetBinContent(1);
 }

 //number of events
 Int_t nEvents4=0;
 if(fAvMult4)
 {
  nEvents4=(Int_t)(fAvMult4->GetBinEntries(1));
 }
 
 TMatrixD dAvG4(cPmax4,cQmax4); 
 Bool_t someAvGEntryIsNegative4=kFALSE;   
 for(Int_t p=0;p<cPmax4;p++)
 {
  for(Int_t q=0;q<cQmax4;q++)
  {
   dAvG4(p,q)=fIntGenFun4->GetBinContent(p+1,q+1);
   if(dAvG4(p,q)<0.)
   {
    someAvGEntryIsNegative4=kTRUE;
   } 
  }
 }  
 /////////////////////////////////////////////////////////////////////////////      
 //////////////////gen. function for the cumulants////////////////////////////
 /////////////////////////////////////////////////////////////////////////////
  
 TMatrixD dC4(cPmax4,cQmax4);//C[p][q]
 if(dAvM4>0 && someAvGEntryIsNegative4==kFALSE)
 {
  for (Int_t p=0;p<cPmax4;p++)
  {
   for (Int_t q=0;q<cQmax4;q++)
   {
    dC4(p,q)=1.*dAvM4*(pow(dAvG4(p,q),(1./dAvM4))-1.); 
   }
  }
 }
 /////////////////////////////////////////////////////////////////////////////
 ///////avaraging the gen. function for the cumulants over azimuth////////////
 /////////////////////////////////////////////////////////////////////////////
  
 TVectorD dAvC4(cPmax4);//<C[p][q]>
 for (Int_t p=0;p<cPmax4;p++)
 {
  Double_t tempHere4=0.; 
  for (Int_t q=0;q<cQmax4;q++)
  {
   tempHere4+=1.*dC4(p,q);
  } 
  dAvC4[p]=1.*tempHere4/cQmax4;
 }
 
 /////////////////////////////////////////////////////////////////////////////
 //////////////////////////////////final results//////////////////////////////
 /////////////////////////////////////////////////////////////////////////////
 
 TVectorD cumulant4(cPmax4);//array to store various order cumulants
  
 cumulant4[0]=(1./(fR0*fR0))*(2.*dAvC[0]-(1./2.)*dAvC[1]);
 cumulant4[1]=(2./pow(fR0,4.))*((-2.)*dAvC[0]+1.*dAvC[1]);
 
 /*      
 cout<<endl;
 cout<<"*********************************"<<endl;
 cout<<"cumulants:"<<endl;
 cout<<" c_"<<cFlow<<"{2} = "<<cumulant4[0]<<endl; 
 cout<<" c_"<<cFlow<<"{4} = "<<cumulant4[1]<<endl;
 cout<<endl;
 */ 
   
 Double_t dV2o4=0.,dV4o4=0.;
 
 if(cumulant4[0]>=0.)
 {
  dV2o4 = pow(cumulant4[0],(1./2.));
 }
 if(cumulant4[1]<=0.)
 {
  dV4o4 = pow(-cumulant4[1],(1./4.));
 }

 cout<<endl;
 cout<<"***********************************"<<endl;
 cout<<"***********************************"<<endl;
 cout<<"flow estimates from GF-cumulants:"<<endl;
 cout<<"  (calculated up to 4th order)   "<<endl;
 cout<<endl;
 
 Double_t sdQo4[2]={0.};
 Double_t chiQo4[2]={0.};
          
 //v_2{2}
 if(nEvents4 && dAvM4 && cumulant4[0]>=0. && (dAvQ2x+dAvQ2y-pow(dAvQx,2.)-pow(dAvQy,2.)-pow(pow(cumulant4[0],(1./2.))*dAvM4,2.)>0.))
 {        
  chiQo4[0]=dAvM4*dV2o4/pow(dAvQ2x+dAvQ2y-pow(dAvQx,2.)-pow(dAvQy,2.)-pow(dV2o4*dAvM4,2.),0.5);
  if(chiQo4[0])
  {
   sdQo4[0]=pow(((1./(2.*dAvM4*nEvents4))*((1.+2.*pow(chiQo4[0],2))/(2.*pow(chiQo4[0],2)))),0.5);
  }
  cout<<" v_"<<cFlow<<"{2} = "<<dV2o4<<" +/- "<<sdQo4[0]<<endl;
  //cout<<" v_"<<cFlow<<"{2} = "<<dV2o4<<" +/- "<<sdQo4[0]<<", chi{2} = "<<chiQo4[0]<<endl;//printing also the chi
 } 
 else 
 {
  cout<<" v_"<<cFlow<<"{2} = Im"<<endl; 
 }
   
 //v_2{4}   
 if(nEvents4 && dAvM4 && cumulant4[1]<=0. && (dAvQ2x+dAvQ2y-pow(dAvQx,2.)-pow(dAvQy,2.)-pow(pow(-cumulant4[1],(1./4.))*dAvM4,2.)>0.))
 {
  chiQo4[1]=dAvM4*dV4o4/pow(dAvQ2x+dAvQ2y-pow(dAvQx,2.)-pow(dAvQy,2.)-pow(dV4o4*dAvM4,2.),0.5);
  if(chiQo4[1])
  {
   sdQo4[1]=(1./(pow(2.*dAvM4*nEvents4,0.5)))*pow((1.+4.*pow(chiQo4[1],2)+1.*pow(chiQo4[1],4.)+2.*pow(chiQo4[1],6.))/(2.*pow(chiQo4[1],6.)),0.5);
  }
  cout<<" v_"<<cFlow<<"{4} = "<<dV4o4<<" +/- "<<sdQo4[1]<<endl;   
  //cout<<" v_"<<cFlow<<"{4} = "<<dV4o4<<" +/- "<<sdQo4[1]<<", chi{4} = "<<chiQo4[1]<<endl;//printing also the chi 
 } 
 else 
 {
  cout<<" v_"<<cFlow<<"{4} = Im"<<endl;  
 } 
  
 cout<<endl;
 cout<<"   nEvts = "<<nEvents4<<", AvM = "<<dAvM4<<endl; 
 cout<<"***********************************"<<endl;    
 cout<<"***********************************"<<endl;   
 
 //============================================================================================================================================== 
 
 //up to 6th order
 //avarage selected multiplicity
 Double_t dAvM6=0.;
 if(fAvMult6)
 {
  dAvM6=fAvMult6->GetBinContent(1);
 }

 //number of events
 Int_t nEvents6=0;
 if(fAvMult6)
 {
  nEvents6=(Int_t)(fAvMult6->GetBinEntries(1));
 }
 
 TMatrixD dAvG6(cPmax6,cQmax6);  
 Bool_t someAvGEntryIsNegative6=kFALSE;   
 for(Int_t p=0;p<cPmax6;p++)
 {
  for(Int_t q=0;q<cQmax6;q++)
  {
   dAvG6(p,q)=fIntGenFun6->GetBinContent(p+1,q+1);
   if(dAvG6(p,q)<0.)
   {
    someAvGEntryIsNegative6=kTRUE;
   } 
  }
 }  
 /////////////////////////////////////////////////////////////////////////////      
 //////////////////gen. function for the cumulants////////////////////////////
 /////////////////////////////////////////////////////////////////////////////
  
 TMatrixD dC6(cPmax6,cQmax6);//C[p][q]
 if(dAvM6>0 && someAvGEntryIsNegative6==kFALSE)
 {
  for (Int_t p=0;p<cPmax6;p++)
  {
   for (Int_t q=0;q<cQmax6;q++)
   {
    dC6(p,q)=1.*dAvM6*(pow(dAvG6(p,q),(1./dAvM6))-1.); 
   }
  }
 }
 
 /////////////////////////////////////////////////////////////////////////////
 ///////avaraging the gen. function for the cumulants over azimuth////////////
 /////////////////////////////////////////////////////////////////////////////
  
 TVectorD dAvC6(cPmax6);//<etBinContent(1)C[p][q]>
 Double_t tempHere6=0.;
 for (Int_t p=0;p<cPmax6;p++){
  tempHere6=0.; 
  for (Int_t q=0;q<cQmax6;q++){
   tempHere6+=1.*dC6(p,q);
  } 
  dAvC6[p]=1.*tempHere6/cQmax6;
 }
 
 /////////////////////////////////////////////////////////////////////////////
 //////////////////////////////////final results//////////////////////////////
 /////////////////////////////////////////////////////////////////////////////
 
 TVectorD cumulant6(cPmax6);//array to store various order cumulants
 cumulant6[0] = (1./(fR0*fR0))*(3.*dAvC[0]-(3./2.)*dAvC[1]+(1./3.)*dAvC[2]);
 cumulant6[1] = (2./pow(fR0,4.))*((-5.)*dAvC[0]+4.*dAvC[1]-1.*dAvC[2]);
 cumulant6[2] = (6./pow(fR0,6.))*(3.*dAvC[0]-3.*dAvC[1]+1.*dAvC[2]);
    
 /*
 cout<<endl;
 cout<<"*********************************"<<endl;
 cout<<"cumulants:"<<endl; 
 cout<<" c_"<<cFlow<<"{2} = "<<cumulant6[0]<<endl; 
 cout<<" c_"<<cFlow<<"{4} = "<<cumulant6[1]<<endl;
 cout<<" c_"<<cFlow<<"{6} = "<<cumulant6[2]<<endl;
 cout<<endl;
 */
 
 Double_t dV2o6=0.,dV4o6=0.,dV6o6=0.;
 
 if(cumulant6[0]>=0.)
 {
  dV2o6 = pow(cumulant6[0],(1./2.));
 }
 if(cumulant6[1]<=0.)
 {
  dV4o6 = pow(-cumulant6[1],(1./4.));
 }
 if(cumulant6[2]>=0.)
 {
  dV6o6 = pow((1./4.)*cumulant6[2],(1./6.));
 }
 
 cout<<endl;
 cout<<"***********************************"<<endl;
 cout<<"***********************************"<<endl;
 cout<<"flow estimates from GF-cumulants:"<<endl;
 cout<<"  (calculated up to 6th order)   "<<endl;
 cout<<endl;
 
 Double_t sdQo6[3]={0.};
 Double_t chiQo6[3]={0.};
          
 //v_2{2}
 if(nEvents6 && dAvM6 && cumulant6[0]>=0. && (dAvQ2x+dAvQ2y-pow(dAvQx,2.)-pow(dAvQy,2.)-pow(pow(cumulant6[0],(1./2.))*dAvM6,2.)>0.))
 {        
  chiQo6[0]=dAvM6*dV2o6/pow(dAvQ2x+dAvQ2y-pow(dAvQx,2.)-pow(dAvQy,2.)-pow(dV2o6*dAvM6,2.),0.5);
  if(chiQo6[0])
  {
   sdQo6[0]=pow(((1./(2.*dAvM6*nEvents6))*((1.+2.*pow(chiQo6[0],2))/(2.*pow(chiQo6[0],2)))),0.5);
  }
  cout<<" v_"<<cFlow<<"{2} = "<<dV2o6<<" +/- "<<sdQo6[0]<<endl;
  //cout<<" v_"<<cFlow<<"{2} = "<<dV2o6<<" +/- "<<sdQo6[0]<<", chi{2} = "<<chiQo6[0]<<endl;//printing also the chi
 } 
 else 
 {
  cout<<" v_"<<cFlow<<"{2} = Im"<<endl; 
 }
   
 //v_2{4}   
 if(nEvents6 && dAvM6 && cumulant6[1]<=0. && (dAvQ2x+dAvQ2y-pow(dAvQx,2.)-pow(dAvQy,2.)-pow(pow(-cumulant6[1],(1./4.))*dAvM6,2.)>0.))
 {
  chiQo6[1]=dAvM6*dV4o6/pow(dAvQ2x+dAvQ2y-pow(dAvQx,2.)-pow(dAvQy,2.)-pow(dV4o6*dAvM6,2.),0.5);
  if(chiQo6[1])
  {
   sdQo6[1]=(1./(pow(2.*dAvM6*nEvents6,0.5)))*pow((1.+4.*pow(chiQo6[1],2)+1.*pow(chiQo6[1],4.)+2.*pow(chiQo6[1],6.))/(2.*pow(chiQo6[1],6.)),0.5);
  }
  cout<<" v_"<<cFlow<<"{4} = "<<dV4o6<<" +/- "<<sdQo6[1]<<endl;    
  //cout<<" v_"<<cFlow<<"{4} = "<<dV4o6<<" +/- "<<sdQo6[1]<<", chi{4} = "<<chiQo6[1]<<endl;//printing also the chi 
 } 
 else 
 {
  cout<<" v_"<<cFlow<<"{4} = Im"<<endl;  
 }   
  
 //v_2{6}
 if(nEvents6 && dAvM6 && cumulant6[2]>=0. && (dAvQ2x+dAvQ2y-pow(dAvQx,2.)-pow(dAvQy,2.)-pow(pow((1./4.)*cumulant6[2],(1./6.))*dAvM6,2.)>0.))
 {
  chiQo6[2]=dAvM6*dV6/pow(dAvQ2x+dAvQ2y-pow(dAvQx,2.)-pow(dAvQy,2.)-pow(dV6o6*dAvM6,2.),0.5);
  if(chiQo6[2])
  {
   sdQo6[2]=(1./(pow(2.*dAvM6*nEvents6,0.5)))*pow((3.+18.*pow(chiQo6[2],2.)+9.*pow(chiQo6[2],4.)+28.*pow(chiQo6[2],6.)+12.*pow(chiQo6[2],8.)+24.*pow(chiQo6[2],10.))/(24.*pow(chiQo6[2],10.)),0.5);
  } 
   cout<<" v_"<<cFlow<<"{6} = "<<dV6o6<<" +/- "<<sdQo6[2]<<endl;   
   //cout<<" v_"<<cFlow<<"{6} = "<<dV6o6<<" +/- "<<sdQo6[2]<<", chi{6} = "<<chiQo6[2]<<endl;//printing also the chi
 }
 else
 {
  cout<<" v_"<<cFlow<<"{6} = Im"<<endl;  
 }
 
 cout<<endl;
 cout<<"   nEvts = "<<nEvents6<<", AvM = "<<dAvM6<<endl; 
 cout<<"***********************************"<<endl;    
 cout<<"***********************************"<<endl;   
 
 //============================================================================================================================================== 
    
 //up to 8th order
 //avarage selected multiplicity
 Double_t dAvM8=0.;
 if(fAvMult8)
 {
  dAvM8=fAvMult8->GetBinContent(1);
 }

 //number of events
 Int_t nEvents8=0;
 if(fAvMult8)
 {
  nEvents8=(Int_t)(fAvMult8->GetBinEntries(1));
 }
 
 TMatrixD dAvG8(cPmax8,cQmax8); 
 Bool_t someAvGEntryIsNegative8=kFALSE;   
 for(Int_t p=0;p<cPmax8;p++)
 {
  for(Int_t q=0;q<cQmax8;q++)
  {
   dAvG8(p,q)=fIntGenFun8->GetBinContent(p+1,q+1);
   if(dAvG8(p,q)<0.)
   {
    someAvGEntryIsNegative8=kTRUE;
   } 
  }
 }  
 /////////////////////////////////////////////////////////////////////////////      
 //////////////////gen. function for the cumulants////////////////////////////
 /////////////////////////////////////////////////////////////////////////////
  
 TMatrixD dC8(cPmax8,cQmax8);//C[p][q]
 if(dAvM8>0 && someAvGEntryIsNegative8==kFALSE)
 {
  for (Int_t p=0;p<cPmax8;p++)
  {
   for (Int_t q=0;q<cQmax8;q++)
   {
    dC8(p,q)=1.*dAvM8*(pow(dAvG8(p,q),(1./dAvM8))-1.); 
   }
  }
 }
 
 /////////////////////////////////////////////////////////////////////////////
 ///////avaraging the gen. function for the cumulants over azimuth////////////
 /////////////////////////////////////////////////////////////////////////////
  
 TVectorD dAvC8(cPmax8);//<C[p][q]>
 Double_t tempHere8=0.;
 for (Int_t p=0;p<cPmax8;p++)
 {
  tempHere8=0.; 
  for (Int_t q=0;q<cQmax8;q++)
  {
   tempHere8+=1.*dC8(p,q);
  } 
  dAvC8[p]=1.*tempHere8/cQmax8;
 }
 
 /////////////////////////////////////////////////////////////////////////////
 //////////////////////////////////final results//////////////////////////////
 /////////////////////////////////////////////////////////////////////////////
 
 TVectorD cumulant8(cPmax8);//array to store various order cumulants
 cumulant8[0] = (1./(fR0*fR0))*(4.*dAvC[0]-3.*dAvC[1]+(4./3.)*dAvC[2]-(1./4.)*dAvC[3]);
 cumulant8[1] = (1./pow(fR0,4.))*((-52./3.)*dAvC[0]+19.*dAvC[1]-(28./3.)*dAvC[2]+(11./6.)*dAvC[3]);
 cumulant8[2] = (3./pow(fR0,6.))*(18.*dAvC[0]-24.*dAvC[1]+14.*dAvC[2]-3.*dAvC[3]);
 cumulant8[3] = (24./pow(fR0,8.))*((-4.)*dAvC[0]+6.*dAvC[1]-4.*dAvC[2]+1.*dAvC[3]);
  
 /* x0,y0 ~ 1/sqrt(p) (setting another complex mesh, doesn't work correctly)
 cumulant8[0] = (-1./(6.*fR0*fR0)) * (1.*dAvC[0] - 48.*dAvC[1] + 243.*dAvC[2] - 256.*dAvC[3]);
 cumulant8[1] = (2./pow(fR0,4.)) * (3.*dAvC[0] - 128.*dAvC[1] + 567.*dAvC[2] - 512.*dAvC[3]);
 cumulant8[2] = (-12./pow(fR0,6.)) * (13.*dAvC[0] - 456.*dAvC[1] + 1701.*dAvC[2] - 1408.*dAvC[3]);
 cumulant8[3] = (2304./pow(fR0,8.)) * (1.*dAvC[0] - 24.*dAvC[1] + 81.*dAvC[2] - 64.*dAvC[3]);       
 */
 
 /* x0,y0 ~ p (setting another complex mesh, doesn't work correctly)
 cumulant8[0] = (-1./(5040.*fR0*fR0)) * ((-8064.)*dAvC[0] + 1008.*dAvC[1] - 128.*dAvC[2] + 9.*dAvC[3]);
 cumulant8[1] = (1./(720.*pow(fR0,4.))) * (1952.*dAvC[0] - 676.*dAvC[1] + 96.*dAvC[2] - 7.*dAvC[3]);
 cumulant8[2] = (-1./(40.*pow(fR0,6.))) * ((-116.)*dAvC[0] + 52.*dAvC[1] - 12.*dAvC[2] + 1.*dAvC[3]);
 cumulant8[3] = (1./(35.*pow(fR0,8.))) * (56.*dAvC[0] - 28.*dAvC[1] + 8.*dAvC[2] - 1.*dAvC[3]);       
 */ 
                                           
 /*                                          
 cout<<endl;
 cout<<"*********************************"<<endl;
 cout<<"cumulants8:"<<endl; 
 cout<<" c_"<<cFlow<<"{2} = "<<cumulant8[0]<<endl; 
 cout<<" c_"<<cFlow<<"{4} = "<<cumulant8[1]<<endl;
 cout<<" c_"<<cFlow<<"{6} = "<<cumulant8[2]<<endl;
 cout<<" c_"<<cFlow<<"{8} = "<<cumulant8[3]<<endl; 
 cout<<endl;
 */
 
 Double_t dV2o8=0.,dV4o8=0.,dV6o8=0.,dV8o8=0.;
 
 if(cumulant8[0]>=0.)
 {
  dV2o8 = pow(cumulant8[0],(1./2.));
 }
 if(cumulant8[1]<=0.)
 {
  dV4o8 = pow(-cumulant8[1],(1./4.));
 }
 if(cumulant8[2]>=0.)
 {
  dV6o8 = pow((1./4.)*cumulant8[2],(1./6.));
 }
 if(cumulant8[3]<=0.)
 {
  dV8o8 = pow(-(1./33.)*cumulant8[3],(1./8.));
 }
 
 cout<<endl;
 cout<<"***********************************"<<endl;
 cout<<"***********************************"<<endl;
 cout<<"flow estimates from GF-cumulants:"<<endl;
 cout<<"  (calculated up to 8th order)   "<<endl;
 cout<<endl;   
             
 Double_t sdQo8[4]={0.};
 Double_t chiQo8[4]={0.};
          
 //v_2{2}
 if(nEvents8 && dAvM8 && cumulant8[0]>=0. && (dAvQ2x+dAvQ2y-pow(dAvQx,2.)-pow(dAvQy,2.)-pow(pow(cumulant8[0],(1./2.))*dAvM8,2.)>0.))
 {        
  chiQo8[0]=dAvM8*dV2o8/pow(dAvQ2x+dAvQ2y-pow(dAvQx,2.)-pow(dAvQy,2.)-pow(dV2o8*dAvM8,2.),0.5);
  if(chiQo8[0])
  {
   sdQo8[0]=pow(((1./(2.*dAvM8*nEvents8))*((1.+2.*pow(chiQo8[0],2.))/(2.*pow(chiQo8[0],2)))),0.5);
  }
  cout<<" v_"<<cFlow<<"{2} = "<<dV2o8<<" +/- "<<sdQo8[0]<<endl;    
  //cout<<" v_"<<cFlow<<"{2} = "<<dV2o8<<" +/- "<<sdQo8[0]<<", chi{2} = "<<chiQo8[0]<<endl;//printing also the chi
 } 
 else 
 {
  cout<<" v_"<<cFlow<<"{2} = Im"<<endl; 
 }
   
 //v_2{4}   
 if(nEvents8 && dAvM8 && cumulant8[1]<=0. && (dAvQ2x+dAvQ2y-pow(dAvQx,2.)-pow(dAvQy,2.)-pow(pow(-cumulant8[1],(1./4.))*dAvM8,2.)>0.))
 {
  chiQo8[1]=dAvM8*dV4o8/pow(dAvQ2x+dAvQ2y-pow(dAvQx,2.)-pow(dAvQy,2.)-pow(dV4o8*dAvM8,2.),0.5);
  if(chiQo8[1])
  {
   sdQo8[1]=(1./(pow(2.*dAvM8*nEvents8,0.5)))*pow((1.+4.*pow(chiQo8[1],2)+1.*pow(chiQo8[1],4.)+2.*pow(chiQo8[1],6.))/(2.*pow(chiQo8[1],6.)),0.5);
  }
  cout<<" v_"<<cFlow<<"{4} = "<<dV4o8<<" +/- "<<sdQo8[1]<<endl;    
  //cout<<" v_"<<cFlow<<"{4} = "<<dV4o8<<" +/- "<<sdQo8[1]<<", chi{4} = "<<chiQo8[1]<<endl;//printing also the chi 
 } 
 else 
 {
  cout<<" v_"<<cFlow<<"{4} = Im"<<endl;  
 } 
  
 //v_2{6}
 if(nEvents8 && dAvM8 && cumulant8[2]>=0. && (dAvQ2x+dAvQ2y-pow(dAvQx,2.)-pow(dAvQy,2.)-pow(pow((1./4.)*cumulant8[2],(1./6.))*dAvM8,2.)>0.))
 {
  chiQo8[2]=dAvM8*dV6o8/pow(dAvQ2x+dAvQ2y-pow(dAvQx,2.)-pow(dAvQy,2.)-pow(dV6o8*dAvM8,2.),0.5);
  if(chiQo8[2])
  {
   sdQo8[2]=(1./(pow(2.*dAvM8*nEvents8,0.5)))*pow((3.+18.*pow(chiQo8[2],2)+9.*pow(chiQo8[2],4.)+28.*pow(chiQo8[2],6.)+12.*pow(chiQo8[2],8.)+24.*pow(chiQo8[2],10.))/(24.*pow(chiQo8[2],10.)),0.5);
  }
  cout<<" v_"<<cFlow<<"{6} = "<<dV6o8<<" +/- "<<sdQo8[2]<<endl;
  //cout<<" v_"<<cFlow<<"{6} = "<<dV6o8<<" +/- "<<sdQo8[2]<<", chi{6} = "<<chiQo8[2]<<endl;//printing also the chi 
 }
 else
 {
  cout<<" v_"<<cFlow<<"{6} = Im"<<endl;  
 }
  
 //v_2{8}
 if(nEvents8 && dAvM8 && cumulant8[3]<=0. && (dAvQ2x+dAvQ2y-pow(dAvQx,2.)-pow(dAvQy,2.)-pow(pow(-(1./33.)*cumulant8[3],(1./8.))*dAvM8,2.)>0.))
 {  
  chiQo8[3]=dAvM8*dV8o8/pow(dAvQ2x+dAvQ2y-pow(dAvQx,2.)-pow(dAvQy,2.)-pow(dV8o8*dAvM8,2.),0.5);
  if(chiQo8[3])
  {
   sdQo8[3]=(1./(pow(2.*dAvM8*nEvents8,0.5)))*pow((12.+96.*pow(chiQo8[3],2)+72.*pow(chiQo8[3],4.)+304.*pow(chiQo8[3],6.)+257.*pow(chiQo8[3],8.)+804.*pow(chiQo8[3],10.)+363.*pow(chiQo8[3],12.)+726.*pow(chiQo8[3],14.))/(726.*pow(chiQo8[3],14.)),0.5);
  } 
  cout<<" v_"<<cFlow<<"{8} = "<<dV8o8<<" +/- "<<sdQo8[3]<<endl;  
  //cout<<" v_"<<cFlow<<"{8} = "<<dV8o8<<" +/- "<<sdQo8[3]<<", chi{8} = "<<chiQo8[3]<<endl;//printing also the chi 
 } 
 else 
 {
  cout<<" v_"<<cFlow<<"{8} = Im"<<endl;     
 }
  
 cout<<endl;
 cout<<"   nEvts = "<<nEvents8<<", AvM = "<<dAvM8<<endl; 
 cout<<"*********************************"<<endl; 

 //============================================================================================================================================== 
   
 //up to 16-th order: 
 //avarage selected multiplicity
 Double_t dAvM16=0.;
 if(fAvMult16)
 {
  dAvM16=fAvMult16->GetBinContent(1);
 }

 //number of events
 Int_t nEvents16=0;
 if(fAvMult16)
 {
  nEvents16=(Int_t)(fAvMult16->GetBinEntries(1));
 }
 
 TMatrixD dAvG16(cPmax16,cQmax16);  
 Bool_t someAvGEntryIsNegative16=kFALSE; 
 for(Int_t p=0;p<cPmax16;p++)
 {
  for(Int_t q=0;q<cQmax16;q++)
  {
   dAvG16(p,q)=fIntGenFun16->GetBinContent(p+1,q+1);
   if(dAvG16(p,q)<0.)
   {
    someAvGEntryIsNegative16=kTRUE;
   } 
  }
 }  
 
 /////////////////////////////////////////////////////////////////////////////      
 //////////////////gen. function for the cumulants////////////////////////////
 /////////////////////////////////////////////////////////////////////////////
  
 TMatrixD dC16(cPmax16,cQmax16);//C16[p][q]
 if(dAvM16>0 && someAvGEntryIsNegative16==kFALSE)
 {
  for(Int_t p=0;p<cPmax16;p++)
  {
   for(Int_t q=0;q<cQmax16;q++)
   {
    dC16(p,q)=1.*dAvM16*(pow(dAvG16(p,q),(1./dAvM16))-1.); 
   }
  }
 }
 
 /////////////////////////////////////////////////////////////////////////////
 ///////avaraging the gen. function for the cumulants over azimuth////////////
 /////////////////////////////////////////////////////////////////////////////
  
 TVectorD dAvC16(cPmax16);//<C16[p][q]>
 Double_t tempHere16=0.; 
 for (Int_t p=0;p<cPmax16;p++)
 {
  tempHere16=0.; 
  for (Int_t q=0;q<cQmax16;q++)
  {
   tempHere16+=1.*dC16(p,q);
  } 
  dAvC16[p]=1.*tempHere16/cQmax16;
 }
 
 /////////////////////////////////////////////////////////////////////////////
 //////////////////////////////////final results//////////////////////////////
 /////////////////////////////////////////////////////////////////////////////
 
 TVectorD cumulant16(cPmax16);//array to store various order cumulants
  
 cumulant16[0] = (1./(fR0*fR0)) * (8.*dAvC16[0] - 14.*dAvC16[1] + (56./3.)*dAvC16[2] - (35./2.)*dAvC16[3] + 
			      (56./5.)*dAvC16[4] - (14./3.)*dAvC16[5] + (8./7.)*dAvC16[6] - (1./8.)*dAvC16[7]);

 cumulant16[1] = (1./pow(fR0,4.)) * ((-1924./35.)*dAvC16[0] + (621./5.)*dAvC16[1] - (8012./45.)*dAvC16[2] + 
				 (691./4.)*dAvC16[3] - (564./5.)*dAvC16[4] + (2143./45.)*dAvC16[5] - 
				 (412./35.)*dAvC16[6] + (363./280.)*dAvC16[7]);

 cumulant16[2] = (1./pow(fR0,6.)) * (349.*dAvC16[0] - (18353./20.)*dAvC16[1] + (7173./5.)*dAvC16[2] - 
				 1457.*dAvC16[3] + (4891./5.)*dAvC16[4] - (1683./4.)*dAvC16[5] + 
				 (527./5.)*dAvC16[6] - (469./40.)*dAvC16[7]);

 cumulant16[3] = (1./pow(fR0,8.)) * ((-10528./5.)*dAvC16[0] + (30578./5.)*dAvC16[1] - (51456./5.)*dAvC16[2] + 
				 10993.*dAvC16[3] - (38176./5.)*dAvC16[4] + (16818./5.)*dAvC16[5] - 
				 (4288./5.)*dAvC16[6] + (967./10.)*dAvC16[7]);

 cumulant16[4] = (1./pow(fR0,10.)) * (11500.*dAvC16[0] - 35800.*dAvC16[1] + 63900.*dAvC16[2] - 71600.*dAvC16[3] + 
				  51620.*dAvC16[4] - 23400.*dAvC16[5] + 6100.*dAvC16[6] - 700.*dAvC16[7]);

 cumulant16[5] = (1./pow(fR0,12.)) * (-52560.*dAvC16[0] + 172080.*dAvC16[1] - 321840.*dAvC16[2] + 376200.*dAvC16[3] - 
				  281520.*dAvC16[4] + 131760.*dAvC16[5] - 35280.*dAvC16[6] + 4140.*dAvC16[7]);

 cumulant16[6] = (1./pow(fR0,14.)) * (176400.*dAvC16[0] - 599760.*dAvC16[1] + 1164240.*dAvC16[2] - 1411200.*dAvC16[3] + 
				  1093680.*dAvC16[4] - 529200.*dAvC16[5] + 146160.*dAvC16[6] - 17640.*dAvC16[7]);

 cumulant16[7] = (1./pow(fR0,16.)) * (-322560*dAvC16[0] + 1128960.*dAvC16[1] - 2257920.*dAvC16[2] + 2822400.*dAvC16[3] - 
				  2257920.*dAvC16[4] + 1128960.*dAvC16[5] - 322560.*dAvC16[6] + 40320.*dAvC16[7]);
    
 /*
 cout<<endl;
 cout<<"*********************************"<<endl;
 cout<<"cumulants:"<<endl;
 cout<<" c_"<<cFlow<<"{2} = "<<cumulant16[0]<<endl; 
 cout<<" c_"<<cFlow<<"{4} = "<<cumulant16[1]<<endl;
 cout<<" c_"<<cFlow<<"{6} = "<<cumulant16[2]<<endl;
 cout<<" c_"<<cFlow<<"{8} = "<<cumulant16[3]<<endl; 
 cout<<"c_"<<cFlow<<"{10} = "<<cumulant16[4]<<endl; 
 cout<<"c_"<<cFlow<<"{12} = "<<cumulant16[5]<<endl;
 cout<<"c_"<<cFlow<<"{14} = "<<cumulant16[6]<<endl; 
 cout<<"c_"<<cFlow<<"{16} = "<<cumulant16[7]<<endl; 
 cout<<endl;
 */
 
 Double_t dV2o16=0.,dV4o16=0.,dV6o16=0.,dV8o16=0.,V10o16=0.,V12o16=0.,V14o16=0.,V16o16=0.;
 
 if(cumulant16[0]>=0.)
 {
  dV2o16 = pow(cumulant16[0],(1./2.));
 }
 if(cumulant16[1]<=0.)
 {
  dV4o16 = pow(-cumulant16[1],(1./4.));
 }
 if(cumulant16[2]>=0.)
 {
  dV6o16 = pow((1./4.)*cumulant16[2],(1./6.));
 }
 if(cumulant16[3]<=0.)
 {
  dV8o16 = pow(-(1./33.)*cumulant16[3],(1./8.));
 }
 if(cumulant16[4]>=0.)
 {
  V10o16 = pow((1./456.)*cumulant16[4],(1./10.));
 }
 if(cumulant16[5]<=0.)
 {
  V12o16 = pow(-(1./9460.)*cumulant16[5],(1./12.));
 }
 if(cumulant16[6]>=0.)
 {
  V14o16 = pow((1./274800.)*cumulant16[6],(1./14.));
 }
 if(cumulant16[7]<=0.)
 {
  V16o16 = pow(-(1./10643745.)*cumulant16[7],(1./16.));
 }
 
 cout<<endl;
 cout<<"***********************************"<<endl;
 cout<<"***********************************"<<endl;
 cout<<"flow estimates from GF-cumulants:"<<endl;
 cout<<"  (calculated up to 16th order)   "<<endl;
 cout<<endl;     
          
 Double_t sdQo16[8]={0.};
 Double_t chiQo16[8]={0.};
          
 //v_2{2}
 if(nEvents16 && dAvM16 && cumulant16[0]>=0. && (dAvQ2x+dAvQ2y-pow(dAvQx,2.)-pow(dAvQy,2.)-pow(pow(cumulant16[0],(1./2.))*dAvM16,2.)>0.))
 {        
  chiQo16[0]=dAvM16*dV2o16/pow(dAvQ2x+dAvQ2y-pow(dAvQx,2.)-pow(dAvQy,2.)-pow(dV2o16*dAvM16,2.),0.5);
  if(chiQo16[0])
  {
   sdQo16[0]=pow(((1./(2.*dAvM16*nEvents16))*((1.+2.*pow(chiQo16[0],2))/(2.*pow(chiQo16[0],2)))),0.5);
  }
  cout<<" v_"<<cFlow<<"{2} = "<<dV2o16<<" +/- "<<sdQo16[0]<<endl;
  //cout<<" v_"<<cFlow<<"{2} = "<<dV2o16<<" +/- "<<sdQo16[0]<<", chi{2} = "<<chiQo16[0]<<endl;//printing also the chi
 } 
 else 
 {
  cout<<" v_"<<cFlow<<"{2} = Im"<<endl; 
 }
   
 //v_2{4}   
 if(nEvents16 && dAvM16 && cumulant16[1]<=0. && (dAvQ2x+dAvQ2y-pow(dAvQx,2.)-pow(dAvQy,2.)-pow(pow(-cumulant16[1],(1./4.))*dAvM16,2.)>0.))
 {
  chiQo16[1]=dAvM16*dV4o16/pow(dAvQ2x+dAvQ2y-pow(dAvQx,2.)-pow(dAvQy,2.)-pow(dV4o16*dAvM16,2.),0.5);
  if(chiQo16[1])
  {
   sdQo16[1]=(1./(pow(2.*dAvM16*nEvents16,0.5)))*pow((1.+4.*pow(chiQo16[1],2.)+1.*pow(chiQo16[1],4.)+2.*pow(chiQo16[1],6.))/(2.*pow(chiQo16[1],6.)),0.5);
  }
  cout<<" v_"<<cFlow<<"{4} = "<<dV4o16<<" +/- "<<sdQo16[1]<<endl; 
  //cout<<" v_"<<cFlow<<"{4} = "<<dV4o16<<" +/- "<<sdQo16[1]<<", chi{4} = "<<chiQo16[1]<<endl;//printing also the chi
 } 
 else 
 {
  cout<<" v_"<<cFlow<<"{4} = Im"<<endl;  
 } 
  
 //v_2{6}
 if(nEvents16 && dAvM16 && cumulant16[2]>=0. && (dAvQ2x+dAvQ2y-pow(dAvQx,2.)-pow(dAvQy,2.)-pow(pow((1./4.)*cumulant16[2],(1./6.))*dAvM16,2.)>0.))
 {
  chiQo16[2]=dAvM16*dV6o16/pow(dAvQ2x+dAvQ2y-pow(dAvQx,2.)-pow(dAvQy,2.)-pow(dV6o16*dAvM16,2.),0.5);
  if(chiQo16[2])
  {
   sdQo16[2]=(1./(pow(2.*dAvM16*nEvents16,0.5)))*pow((3.+18.*pow(chiQo16[2],2)+9.*pow(chiQo16[2],4.)+28.*pow(chiQo16[2],6.)+12.*pow(chiQo16[2],8.)+24.*pow(chiQo16[2],10.))/(24.*pow(chiQo16[2],10.)),0.5);
  }
  cout<<" v_"<<cFlow<<"{6} = "<<dV6o16<<" +/- "<<sdQo16[2]<<endl;
  //cout<<" v_"<<cFlow<<"{6} = "<<dV6o16<<" +/- "<<sdQo16[2]<<", chi{6} = "<<chiQo16[2]<<endl;//printing also the chi
 }
 else
 {
  cout<<" v_"<<cFlow<<"{6} = Im"<<endl;  
 }
  
 //v_2{8}
 if(nEvents16 && dAvM16 && cumulant16[3]<=0. && (dAvQ2x+dAvQ2y-pow(dAvQx,2.)-pow(dAvQy,2.)-pow(pow(-(1./33.)*cumulant16[3],(1./8.))*dAvM16,2.)>0.))
 {  
  chiQo16[3]=dAvM16*dV8o16/pow(dAvQ2x+dAvQ2y-pow(dAvQx,2.)-pow(dAvQy,2.)-pow(dV8o16*dAvM16,2.),0.5);
  if(chiQo16[3])
  {
   sdQo16[3]=(1./(pow(2.*dAvM16*nEvents16,0.5)))*pow((12.+96.*pow(chiQo16[3],2)+72.*pow(chiQo16[3],4.)+304.*pow(chiQo16[3],6.)+257.*pow(chiQo16[3],8.)+804.*pow(chiQo16[3],10.)+363.*pow(chiQo16[3],12.)+726.*pow(chiQo16[3],14.))/(726.*pow(chiQo16[3],14.)),0.5);
  } 
  cout<<" v_"<<cFlow<<"{8} = "<<dV8o16<<" +/- "<<sdQo16[3]<<endl;
  //cout<<" v_"<<cFlow<<"{8} = "<<dV8o16<<" +/- "<<sdQo16[3]<<", chi{8} = "<<chiQo16[3]<<endl;//printing also the chi
 } 
 else 
 {
  cout<<" v_"<<cFlow<<"{8} = Im"<<endl;     
 }
  
 //v_2{10}
 if(nEvents16 && dAvM16 && cumulant16[4]>=0.)
 {
  cout<<"v_"<<cFlow<<"{10} = "<<V10o16<<endl;
 } 
 else 
 {
  cout<<"v_"<<cFlow<<"{10} = Im"<<endl; 
 }
  
 //v_2{12}
 if(nEvents16 && dAvM16 && cumulant16[5]<=0.)
 {
  cout<<"v_"<<cFlow<<"{12} = "<<V12o16<<endl;
 } 
 else 
 {
  cout<<"v_"<<cFlow<<"{12} = Im"<<endl; 
 }
  
 //v_2{14}
 if(nEvents16 && dAvM16 && cumulant16[6]>=0.)
 {
  cout<<"v_"<<cFlow<<"{14} = "<<V14o16<<endl;
 } 
 else 
 {
  cout<<"v_"<<cFlow<<"{14} = Im"<<endl;  
 }
  
 //v_2{16}
 if(nEvents16 && dAvM16 && cumulant16[7]<=0.)
 {
  cout<<"v_"<<cFlow<<"{16} = "<<V16o16<<endl;
 } 
 else 
 {
  cout<<"v_"<<cFlow<<"{16} = Im"<<endl;  
 }
  
 cout<<endl;
 cout<<"   nEvts = "<<nEvents16<<", AvM = "<<dAvM16<<endl; 
 cout<<"*********************************"<<endl;    
 
 //==============================================================================================================================================
 
 }//end of if(otherEquations) 
}

  

