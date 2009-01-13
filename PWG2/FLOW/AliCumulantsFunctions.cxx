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
 fDiffGenFunRe(NULL),
 fDiffGenFunIm(NULL),
 fBinNoOfParticles(NULL),
 fifr(NULL),
 fdfr2(NULL), 
 fdfr4(NULL), 
 fdfr6(NULL),
 fdfr8(NULL),
 fAvMult(NULL),
 fQVector(NULL),
 fchr2nd(NULL),
 fchr4th(NULL),
 fchr6th(NULL),
 fchr8th(NULL) 
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

AliCumulantsFunctions::AliCumulantsFunctions(TProfile2D *IntGenFun, TProfile2D *IntGenFun4, TProfile2D *IntGenFun6, TProfile2D *IntGenFun8, TProfile2D *IntGenFun16, TProfile *AvMult4, TProfile *AvMult6, TProfile *AvMult8, TProfile *AvMult16, TProfile3D *DiffGenFunRe, TProfile3D *DiffGenFunIm, TProfile *BinNoOfParticles, TH1D *ifr, TH1D *dfr2, TH1D *dfr4, TH1D *dfr6, TH1D *dfr8, TProfile *AvMult, TProfile *QVector, AliFlowCommonHistResults *chr2nd, AliFlowCommonHistResults *chr4th, AliFlowCommonHistResults *chr6th, AliFlowCommonHistResults *chr8th):
 fIntGenFun(IntGenFun),
 fIntGenFun4(IntGenFun4),//only for other system of Eq.
 fIntGenFun6(IntGenFun6),//only for other system of Eq.
 fIntGenFun8(IntGenFun8),//only for other system of Eq.
 fIntGenFun16(IntGenFun16),//only for other system of Eq.
 fAvMult4(AvMult4),//only for other system of Eq.
 fAvMult6(AvMult6),//only for other system of Eq.
 fAvMult8(AvMult8),//only for other system of Eq.
 fAvMult16(AvMult16),//only for other system of Eq.
 fDiffGenFunRe(DiffGenFunRe),
 fDiffGenFunIm(DiffGenFunIm),
 fBinNoOfParticles(BinNoOfParticles),
 fifr(ifr),
 fdfr2(dfr2), 
 fdfr4(dfr4), 
 fdfr6(dfr6),
 fdfr8(dfr8),
 fAvMult(AvMult),
 fQVector(QVector),
 fchr2nd(chr2nd),
 fchr4th(chr4th),
 fchr6th(chr6th),
 fchr8th(chr8th) 
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
 static const Int_t fgkQmax=AliFlowCumuConstants::kQmax;     //needed for numerics
 static const Int_t fgkPmax=AliFlowCumuConstants::kPmax;     //needed for numerics  
 static const Int_t fgkQmax4=AliFlowCumuConstants::kQmax4;   //needed for numerics
 static const Int_t fgkPmax4=AliFlowCumuConstants::kPmax4;   //needed for numerics
 static const Int_t fgkQmax6=AliFlowCumuConstants::kQmax6;   //needed for numerics
 static const Int_t fgkPmax6=AliFlowCumuConstants::kPmax6;   //needed for numerics
 static const Int_t fgkQmax8=AliFlowCumuConstants::kQmax8;   //needed for numerics
 static const Int_t fgkPmax8=AliFlowCumuConstants::kPmax8;   //needed for numerics
 static const Int_t fgkQmax16=AliFlowCumuConstants::kQmax16; //needed for numerics
 static const Int_t fgkPmax16=AliFlowCumuConstants::kPmax16; //needed for numerics
 
 static const Int_t fgkFlow=AliFlowCumuConstants::kFlow;   //integrated flow coefficient to be calculated
 static const Int_t fgkMltpl=AliFlowCumuConstants::kMltpl; //the multiple in p=m*n (diff. flow) 
 static const Int_t fgknBins=100;                          //number of pt bins //to be improved
 
 Double_t fR0=AliFlowCumuConstants::fgR0;            //needed for numerics
 //Double_t fPtMax=AliFlowCommonConstants::GetPtMax(); //maximum pt
 //Double_t fPtMin=AliFlowCommonConstants::GetPtMin(); //minimum pt
 //Double_t fBinWidth=(fPtMax-fPtMin)/fgknBins;        //width of pt bin (in GeV)   
 
 Bool_t fOtherEquations=AliFlowCumuConstants::fgOtherEquations;     //numerical equations for cumulants solved up to different highest order 
 
 //avarage selected multiplicity
 Double_t AvM=0.;
 if(fAvMult)
 {
  AvM=fAvMult->GetBinContent(1);
 }

 //number of events
 Int_t nEvents=0;
 if(fAvMult)
 {
  nEvents=(Int_t)(fAvMult->GetBinEntries(1));
 }
 
 //<Q-vector stuff>
 Double_t AvQx=0.,AvQy=0.,AvQ2x=0.,AvQ2y=0.;
 if(fQVector)
 {
  AvQx  = fQVector->GetBinContent(1); //<Q_x>
  AvQy  = fQVector->GetBinContent(2); //<Q_y>
  AvQ2x = fQVector->GetBinContent(3); //<(Q_x)^2>
  AvQ2y = fQVector->GetBinContent(4); //<(Q_y)^2>
 }
 
 //<G[p][q]>
 Double_t AvG[fgkPmax][fgkQmax]={{0.}}; 
 Bool_t someAvGEntryIsNegative=kFALSE;
 if(fIntGenFun)
 {   
  for(Int_t p=0;p<fgkPmax;p++)
  {
   for(Int_t q=0;q<fgkQmax;q++)
   {
    AvG[p][q]=fIntGenFun->GetBinContent(p+1,q+1);
    if(AvG[p][q]<0.)
    {
     someAvGEntryIsNegative=kTRUE;
    } 
   }
  }  
 } 
   
 /////////////////////////////////////////////////////////////////////////////      
 //////////////////gen. function for the cumulants////////////////////////////
 /////////////////////////////////////////////////////////////////////////////
    
 Double_t C[fgkPmax][fgkQmax]={{0.}};//C[p][q]
 if(AvM>0 && someAvGEntryIsNegative==kFALSE)
 {
  for(Int_t p=0;p<fgkPmax;p++)
  {
   for(Int_t q=0;q<fgkQmax;q++)
   {
    C[p][q]=1.*AvM*(pow(AvG[p][q],(1./AvM))-1.); 
   }
  }
 }
    
 /////////////////////////////////////////////////////////////////////////////
 ///////avaraging the gen. function for the cumulants over azimuth////////////
 /////////////////////////////////////////////////////////////////////////////
  
 Double_t AvC[fgkPmax]={0.};//<C[p][q]>
 for(Int_t p=0;p<fgkPmax;p++)
 {
  Double_t tempHere=0.; 
  for(Int_t q=0;q<fgkQmax;q++)
  {
   tempHere+=1.*C[p][q];
  } 
  AvC[p]=1.*tempHere/fgkQmax;
 }
 
 /////////////////////////////////////////////////////////////////////////////
 //////////////////////////////////final results//////////////////////////////
 /////////////////////////////////////////////////////////////////////////////
 
 Double_t cumulant[fgkPmax];//array to store various order cumulants
 
 //system of eq. for the cumulants  
 cumulant[0] = (-1./(60*fR0*fR0))*((-300.)*AvC[0]+300.*AvC[1]-200.*AvC[2]+75.*AvC[3]-12.*AvC[4]);
 cumulant[1] = (-1./(6.*pow(fR0,4.)))*(154.*AvC[0]-214.*AvC[1]+156.*AvC[2]-61.*AvC[3]+10.*AvC[4]);
 cumulant[2] = (3./(2.*pow(fR0,6.)))*(71.*AvC[0]-118.*AvC[1]+98.*AvC[2]-41.*AvC[3]+7.*AvC[4]);
 cumulant[3] = (-24./pow(fR0,8.))*(14.*AvC[0]-26.*AvC[1]+24.*AvC[2]-11.*AvC[3]+2.*AvC[4]);
 cumulant[4] = (120./pow(fR0,10.))*(5.*AvC[0]-10.*AvC[1]+10.*AvC[2]-5.*AvC[3]+1.*AvC[4]);
    
 /*
 cout<<endl;
 cout<<"*********************************"<<endl;
 cout<<"cumulants:"<<endl; 
 cout<<" c_"<<fgkFlow<<"{2} = "<<cumulant[0]<<endl; 
 cout<<" c_"<<fgkFlow<<"{4} = "<<cumulant[1]<<endl;
 cout<<" c_"<<fgkFlow<<"{6} = "<<cumulant[2]<<endl;
 cout<<" c_"<<fgkFlow<<"{8} = "<<cumulant[3]<<endl; 
 cout<<"c_"<<fgkFlow<<"{10} = "<<cumulant[4]<<endl;  
 cout<<endl;
 */
 
 Double_t V2=0.,V4=0.,V6=0.,V8=0.,V10=0.;//integrated flow estimates
 
 if(cumulant[0]>=0.)
 {
  V2=pow(cumulant[0],(1./2.));
 }
 if(cumulant[1]<=0.)
 {
  V4=pow(-cumulant[1],(1./4.));
 }
 if(cumulant[2]>=0.)
 {
  V6=pow((1./4.)*cumulant[2],(1./6.));
 }
 if(cumulant[3]<=0.)
 {
  V8=pow(-(1./33.)*cumulant[3],(1./8.));
 }
 if(cumulant[4]>=0.)
 {
  V10=pow((1./456.)*cumulant[4],(1./10.));
 }

 cout<<endl;
 cout<<"***********************************"<<endl;
 cout<<"***********************************"<<endl;
 cout<<"flow estimates from GF-cumulants:"<<endl;
 cout<<endl;          
                
 Double_t SdQ[4]={0.};
 Double_t ChiQ[4]={0.};
                            
 //v_2{2}
 if(nEvents>0 && AvM>0 && cumulant[0]>=0. && (AvQ2x+AvQ2y-pow(AvQx,2.)-pow(AvQy,2.)-pow(pow(cumulant[0],(1./2.))*AvM,2.)>0.))
 {        
  ChiQ[0]=AvM*V2/pow(AvQ2x+AvQ2y-pow(AvQx,2.)-pow(AvQy,2.)-pow(V2*AvM,2.),0.5);
  if(ChiQ[0])
  {  
   SdQ[0]=pow(((1./(2.*AvM*nEvents))*((1.+2.*pow(ChiQ[0],2))/(2.*pow(ChiQ[0],2)))),0.5);
  }
  cout<<" v_"<<fgkFlow<<"{2} = "<<V2<<" +/- "<<SdQ[0]<<endl;
  //cout<<" v_"<<fgkFlow<<"{2} = "<<V2<<" +/- "<<SdQ[0]<<", chi{2} = "<<ChiQ[0]<<endl;//printing also the chi
  fifr->SetBinContent(1,V2);
  fifr->SetBinError(1,SdQ[0]);
  //filling common histograms:
  fchr2nd->FillIntegratedFlow(V2,SdQ[0]);
  fchr2nd->FillChi(ChiQ[0]);
 } 
 else 
 {
  cout<<" v_"<<fgkFlow<<"{2} = Im"<<endl; 
 }
   
 //v_2{4}   
 if(nEvents>0 && AvM>0 && cumulant[1]<=0. && (AvQ2x+AvQ2y-pow(AvQx,2.)-pow(AvQy,2.)-pow(pow(-cumulant[1],(1./4.))*AvM,2.)>0.))
 {
  ChiQ[1]=AvM*V4/pow(AvQ2x+AvQ2y-pow(AvQx,2.)-pow(AvQy,2.)-pow(V4*AvM,2.),0.5);
  if(ChiQ[1])
  {
   SdQ[1]=(1./(pow(2.*AvM*nEvents,0.5)))*pow((1.+4.*pow(ChiQ[1],2)+1.*pow(ChiQ[1],4.)+2.*pow(ChiQ[1],6.))/(2.*pow(ChiQ[1],6.)),0.5);
  }
  cout<<" v_"<<fgkFlow<<"{4} = "<<V4<<" +/- "<<SdQ[1]<<endl;
  //cout<<" v_"<<fgkFlow<<"{4} = "<<V4<<" +/- "<<SdQ[1]<<", chi{4} = "<<ChiQ[1]<<endl;//printing also the chi
  fifr->SetBinContent(2,V4);
  fifr->SetBinError(2,SdQ[1]);
  //filling common histograms:
  fchr4th->FillIntegratedFlow(V4,SdQ[1]);
  fchr4th->FillChi(ChiQ[1]);
 } 
 else 
 {
  cout<<" v_"<<fgkFlow<<"{4} = Im"<<endl;  
 } 
  
 //v_2{6}
 if(nEvents>0 && AvM>0 && cumulant[2]>=0. && (AvQ2x+AvQ2y-pow(AvQx,2.)-pow(AvQy,2.)-pow(pow((1./4.)*cumulant[2],(1./6.))*AvM,2.)>0.))
 {
  ChiQ[2]=AvM*V6/pow(AvQ2x+AvQ2y-pow(AvQx,2.)-pow(AvQy,2.)-pow(V6*AvM,2.),0.5);
  if(ChiQ[2])
  {
   SdQ[2]=(1./(pow(2.*AvM*nEvents,0.5)))*pow((3.+18.*pow(ChiQ[2],2)+9.*pow(ChiQ[2],4.)+28.*pow(ChiQ[2],6.)+12.*pow(ChiQ[2],8.)+24.*pow(ChiQ[2],10.))/(24.*pow(ChiQ[2],10.)),0.5);
  } 
  cout<<" v_"<<fgkFlow<<"{6} = "<<V6<<" +/- "<<SdQ[2]<<endl;
  //cout<<" v_"<<fgkFlow<<"{6} = "<<V6<<" +/- "<<SdQ[2]<<", chi{6} = "<<ChiQ[2]<<endl;//printing also the chi
  fifr->SetBinContent(3,V6);
  fifr->SetBinError(3,SdQ[2]);
  //filling common histograms:
  fchr6th->FillIntegratedFlow(V6,SdQ[2]);
  fchr6th->FillChi(ChiQ[2]);
 }
 else
 {
  cout<<" v_"<<fgkFlow<<"{6} = Im"<<endl;  
 }
  
 //v_2{8}
 if(nEvents>0 && AvM>0 && cumulant[3]<=0. && (AvQ2x+AvQ2y-pow(AvQx,2.)-pow(AvQy,2.)-pow(pow(-(1./33.)*cumulant[3],(1./8.))*AvM,2.)>0.))
 {  
  ChiQ[3]=AvM*V8/pow(AvQ2x+AvQ2y-pow(AvQx,2.)-pow(AvQy,2.)-pow(V8*AvM,2.),0.5);
  if(ChiQ[3])
  {
   SdQ[3]=(1./(pow(2.*AvM*nEvents,0.5)))*pow((12.+96.*pow(ChiQ[3],2.)+72.*pow(ChiQ[3],4.)+304.*pow(ChiQ[3],6.)+257.*pow(ChiQ[3],8.)+804.*pow(ChiQ[3],10.)+363.*pow(ChiQ[3],12.)+726.*pow(ChiQ[3],14.))/(726.*pow(ChiQ[3],14.)),0.5);
  } 
  cout<<" v_"<<fgkFlow<<"{8} = "<<V8<<" +/- "<<SdQ[3]<<endl;
  //cout<<" v_"<<fgkFlow<<"{8} = "<<V8<<" +/- "<<SdQ[3]<<", chi{8} = "<<ChiQ[3]<<endl;//printing also the chi
  fifr->SetBinContent(4,V8);
  fifr->SetBinError(4,SdQ[3]);
  //filling common histograms:
  fchr8th->FillIntegratedFlow(V8,SdQ[3]);
  fchr8th->FillChi(ChiQ[3]);
 } 
 else 
 {
  cout<<" v_"<<fgkFlow<<"{8} = Im"<<endl;     
 }
  
 /*
 //v_2{10}
 if (AvM && cumulant[4]>=0.)
 {
  cout<<"v_"<<fgkFlow<<"{10} = "<<V10<<endl;
 }
 else 
 {
  cout<<"v_"<<fgkFlow<<"{10} = Im"<<endl; 
 }
 */
 
 cout<<endl;
 cout<<"   nEvts = "<<nEvents<<", AvM = "<<AvM<<endl; 
 cout<<"***********************************"<<endl;    
 cout<<"***********************************"<<endl;   
  
 
 
 
 
 /////////////////////////////////////////////////////////////////////////////
 ///////////////////////DIFFERENTIAL FLOW CALCULATIONS////////////////////////
 /////////////////////////////////////////////////////////////////////////////
  
  Double_t X[fgknBins][fgkPmax][fgkQmax]={{{0.}}};
  Double_t Y[fgknBins][fgkPmax][fgkQmax]={{{0.}}};
  Double_t BinNoOfParticles[fgknBins]={0.};
      
  //3D profiles
  for(Int_t b=0;b<fgknBins;b++)
  {
   BinNoOfParticles[b]=fBinNoOfParticles->GetBinEntries(b);
   for(Int_t p=0;p<fgkPmax;p++)
   {
    for(Int_t q=0;q<fgkQmax;q++)
    {
     if(AvG[p][q])
     {   
      X[b][p][q]=fDiffGenFunRe->GetBinContent(b+1,p+1,q+1)/AvG[p][q];
      Y[b][p][q]=fDiffGenFunIm->GetBinContent(b+1,p+1,q+1)/AvG[p][q];
     } 
    }
   }   
  } 
  
  /*
  if(AvM){
  for(Int_t b=0;b<fgknBins;b++){
    //for(Int_t p=0;p<fgkPmax;p++){
      for(Int_t q=0;q<fgkQmax;q++){
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
  
  Double_t D[fgknBins][fgkPmax]={{0.}};
  
  for (Int_t b=0;b<fgknBins;b++){
    for (Int_t p=0;p<fgkPmax;p++){
      Double_t fTempHere3=0.; 
      for (Int_t q=0;q<fgkQmax;q++){
	fTempHere3+=cos(fgkMltpl*2.*q*TMath::Pi()/fgkQmax)*X[b][p][q] + sin(fgkMltpl*2.*q*TMath::Pi()/fgkQmax)*Y[b][p][q];
      } 
      D[b][p]=1.*(pow(fR0*pow(p+1.,.5),fgkMltpl)/fgkQmax)*fTempHere3;
    }
  } 
  
  Double_t DiffCumulant2[fgknBins]={0.};
  Double_t DiffCumulant4[fgknBins]={0.};
  Double_t DiffCumulant6[fgknBins]={0.};
  Double_t DiffCumulant8[fgknBins]={0.};
  
  for (Int_t b=0;b<fgknBins;b++){
    DiffCumulant2[b]=(1./(fR0*fR0))*(4.*D[b][0]-3.*D[b][1]+(4./3.)*D[b][2]-(1./4.)*D[b][3]);
    DiffCumulant4[b]=(1./pow(fR0,4.))*((-26./3.)*D[b][0]+(19./2.)*D[b][1]-(14./3.)*D[b][2]+(11./12.)*D[b][3]);
    DiffCumulant6[b]=(1./pow(fR0,6.))*(18.*D[b][0]-24.*D[b][1]+14.*D[b][2]-3.*D[b][3]);
    DiffCumulant8[b]=(1./pow(fR0,8.))*(-24.*D[b][0]+36.*D[b][1]-24.*D[b][2]+6.*D[b][3]);
  }
  
  Double_t v2[fgknBins],v4[fgknBins],v6[fgknBins],v8[fgknBins];
  Double_t Sddiff2[fgknBins],Sddiff4[fgknBins];

  //cout<<"number of pt bins: "<<fgknBins<<endl;
  //cout<<"****************************************"<<endl;
  for (Int_t b=0;b<fgknBins;b++){ 
    //cout<<"pt bin: "<<b*fBinWidth<<"-"<<(b+1)*fBinWidth<<" GeV"<<endl;
    
    //v'_{2/2}{2}
    if(cumulant[0]>0)
    {
      v2[b]=DiffCumulant2[b]/pow(cumulant[0],.5);
      if (AvQ2x+AvQ2y-pow(AvQx,2.)-pow(AvQy,2.)-pow(V2*AvM,2.)>0.&&BinNoOfParticles[b]>0.)
      {
       Sddiff2[b]=pow((1./(2.*BinNoOfParticles[b]))*((1.+pow(ChiQ[0],2.))/pow(ChiQ[0],2.)),0.5);
       //cout<<"v'_2/2{2} = "<<v2[b]<<"%, "<<" "<<"sd{2} = "<<100.*Sddiff2[b]<<"%"<<endl;
       fdfr2->SetBinContent(b+1,v2[b]);
       fdfr2->SetBinError(b+1,Sddiff2[b]);
       //common histogram:
       fchr2nd->FillDifferentialFlow(b+1,v2[b],Sddiff2[b]);
      } else {
         //cout<<"v'_2/2{2} = Im"<<endl;
      }
    }else{
      //cout<<"v'_2/2{2} = Im"<<endl;
    } 
    
    //v'_{2/2}{4}
    if(cumulant[1]<0)
    {
      v4[b]=-DiffCumulant4[b]/pow(-cumulant[1],.75);
      if (AvQ2x+AvQ2y-pow(AvQx,2.)-pow(AvQy,2.)-pow(V4*AvM,2.)>0.&&BinNoOfParticles[b]>0.)
      {
       Sddiff4[b]=pow((1./(2.*BinNoOfParticles[b]))*((2.+6.*pow(ChiQ[1],2.)+pow(ChiQ[1],4.)+pow(ChiQ[1],6.))/pow(ChiQ[1],6.)),0.5);
       //cout<<"v'_2/2{4} = "<<v4[b]<<"%, "<<" "<<"sd{4} = "<<100.*Sddiff4[b]<<"%"<<endl;
       fdfr4->SetBinContent(b+1,v4[b]);
       fdfr4->SetBinError(b+1,Sddiff4[b]);
       //common histogram:
       fchr4th->FillDifferentialFlow(b+1,v4[b],Sddiff4[b]);
      } else {
         //cout<<"v'_2/2{4} = Im"<<endl;
      } 
    }else{
      //cout<<"v'_2/2{4} = Im"<<endl;
    }  
    
    //v'_{2/2}{6}
    if(cumulant[2]>0){
      //cout<<"v'_2/2{6} = "<<100.*DiffCumulant6[b]/(4.*pow((1./4.)*cumulant[2],(5./6.)))<<"%"<<endl;
      v6[b]=DiffCumulant6[b]/(4.*pow((1./4.)*cumulant[2],(5./6.)));
      fdfr6->SetBinContent(b+1,v6[b]);
      //common histogram:
      fchr6th->FillDifferentialFlow(b+1,v6[b],0.);
    }else{
      //cout<<"v'_2/2{6} = Im"<<endl;
    }     
    
    //v'_{2/2}{8}
    if(cumulant[3]<0){
      //cout<<"v'_2/2{8} = "<<-100.*DiffCumulant8[b]/(33.*pow(-(1./33.)*cumulant[3],(7./8.)))<<"%"<<endl;
      v8[b]=-DiffCumulant8[b]/(33.*pow(-(1./33.)*cumulant[3],(7./8.))); 
      fdfr8->SetBinContent(b+1,v8[b]);
      //common histogram:
      fchr8th->FillDifferentialFlow(b+1,v8[b],0.);
    }else{
      //cout<<"v'_2/2{8} = Im"<<endl;
    }       
    //cout<<"****************************************"<<endl;
  }  
  
  
 
 
 
 
 
 
 //off the record: numerical equations for cumulants solved up to different highest order  
 if(fOtherEquations)
 {
 
 //==============================================================================================================================================
  
 //up to 4th order
 //avarage selected multiplicity
 Double_t AvM4=0.;
 if(fAvMult4)
 {
  AvM4=fAvMult4->GetBinContent(1);
 }

 //number of events
 Int_t nEvents4=0;
 if(fAvMult4)
 {
  nEvents4=(Int_t)(fAvMult4->GetBinEntries(1));
 }
 
 Double_t AvG4[fgkPmax4][fgkQmax4]={{0.}}; 
 Bool_t someAvGEntryIsNegative4=kFALSE;   
 for(Int_t p=0;p<fgkPmax4;p++)
 {
  for(Int_t q=0;q<fgkQmax4;q++)
  {
   AvG4[p][q]=fIntGenFun4->GetBinContent(p+1,q+1);
   if(AvG4[p][q]<0.)
   {
    someAvGEntryIsNegative4=kTRUE;
   } 
  }
 }  
 /////////////////////////////////////////////////////////////////////////////      
 //////////////////gen. function for the cumulants////////////////////////////
 /////////////////////////////////////////////////////////////////////////////
  
 Double_t C4[fgkPmax4][fgkQmax4]={{0.}};//C[p][q]
 if(AvM4>0 && someAvGEntryIsNegative4==kFALSE)
 {
  for (Int_t p=0;p<fgkPmax4;p++)
  {
   for (Int_t q=0;q<fgkQmax4;q++)
   {
    C4[p][q]=1.*AvM4*(pow(AvG4[p][q],(1./AvM4))-1.); 
   }
  }
 }
 /////////////////////////////////////////////////////////////////////////////
 ///////avaraging the gen. function for the cumulants over azimuth////////////
 /////////////////////////////////////////////////////////////////////////////
  
 Double_t AvC4[fgkPmax4]={0.};//<C[p][q]>
 for (Int_t p=0;p<fgkPmax4;p++)
 {
  Double_t tempHere4=0.; 
  for (Int_t q=0;q<fgkQmax4;q++)
  {
   tempHere4+=1.*C4[p][q];
  } 
  AvC4[p]=1.*tempHere4/fgkQmax4;
 }
 
 /////////////////////////////////////////////////////////////////////////////
 //////////////////////////////////final results//////////////////////////////
 /////////////////////////////////////////////////////////////////////////////
 
 Double_t cumulant4[fgkPmax4];//array to store various order cumulants
  
 cumulant4[0]=(1./(fR0*fR0))*(2.*AvC[0]-(1./2.)*AvC[1]);
 cumulant4[1]=(2./pow(fR0,4.))*((-2.)*AvC[0]+1.*AvC[1]);
 
 /*      
 cout<<endl;
 cout<<"*********************************"<<endl;
 cout<<"cumulants:"<<endl;
 cout<<" c_"<<fgkFlow<<"{2} = "<<cumulant4[0]<<endl; 
 cout<<" c_"<<fgkFlow<<"{4} = "<<cumulant4[1]<<endl;
 cout<<endl;
 */ 
   
 Double_t V2o4=0.,V4o4=0.;
 
 if(cumulant4[0]>=0.)
 {
  V2o4 = pow(cumulant4[0],(1./2.));
 }
 if(cumulant4[1]<=0.)
 {
  V4o4 = pow(-cumulant4[1],(1./4.));
 }

 cout<<endl;
 cout<<"***********************************"<<endl;
 cout<<"***********************************"<<endl;
 cout<<"flow estimates from GF-cumulants:"<<endl;
 cout<<"  (calculated up to 4th order)   "<<endl;
 cout<<endl;
 
 Double_t SdQo4[2]={0.};
 Double_t ChiQo4[2]={0.};
          
 //v_2{2}
 if(nEvents4 && AvM4 && cumulant4[0]>=0. && (AvQ2x+AvQ2y-pow(AvQx,2.)-pow(AvQy,2.)-pow(pow(cumulant4[0],(1./2.))*AvM4,2.)>0.))
 {        
  ChiQo4[0]=AvM4*V2o4/pow(AvQ2x+AvQ2y-pow(AvQx,2.)-pow(AvQy,2.)-pow(V2o4*AvM4,2.),0.5);
  if(ChiQo4[0])
  {
   SdQo4[0]=pow(((1./(2.*AvM4*nEvents4))*((1.+2.*pow(ChiQo4[0],2))/(2.*pow(ChiQo4[0],2)))),0.5);
  }
  cout<<" v_"<<fgkFlow<<"{2} = "<<V2o4<<" +/- "<<SdQo4[0]<<endl;
  //cout<<" v_"<<fgkFlow<<"{2} = "<<V2o4<<" +/- "<<SdQo4[0]<<", chi{2} = "<<ChiQo4[0]<<endl;//printing also the chi
 } 
 else 
 {
  cout<<" v_"<<fgkFlow<<"{2} = Im"<<endl; 
 }
   
 //v_2{4}   
 if(nEvents4 && AvM4 && cumulant4[1]<=0. && (AvQ2x+AvQ2y-pow(AvQx,2.)-pow(AvQy,2.)-pow(pow(-cumulant4[1],(1./4.))*AvM4,2.)>0.))
 {
  ChiQo4[1]=AvM4*V4o4/pow(AvQ2x+AvQ2y-pow(AvQx,2.)-pow(AvQy,2.)-pow(V4o4*AvM4,2.),0.5);
  if(ChiQo4[1])
  {
   SdQo4[1]=(1./(pow(2.*AvM4*nEvents4,0.5)))*pow((1.+4.*pow(ChiQo4[1],2)+1.*pow(ChiQo4[1],4.)+2.*pow(ChiQo4[1],6.))/(2.*pow(ChiQo4[1],6.)),0.5);
  }
  cout<<" v_"<<fgkFlow<<"{4} = "<<V4o4<<" +/- "<<SdQo4[1]<<endl;   
  //cout<<" v_"<<fgkFlow<<"{4} = "<<V4o4<<" +/- "<<SdQo4[1]<<", chi{4} = "<<ChiQo4[1]<<endl;//printing also the chi 
 } 
 else 
 {
  cout<<" v_"<<fgkFlow<<"{4} = Im"<<endl;  
 } 
  
 cout<<endl;
 cout<<"   nEvts = "<<nEvents4<<", AvM = "<<AvM4<<endl; 
 cout<<"***********************************"<<endl;    
 cout<<"***********************************"<<endl;   
 
 //============================================================================================================================================== 
 
 //up to 6th order
 //avarage selected multiplicity
 Double_t AvM6=0.;
 if(fAvMult6)
 {
  AvM6=fAvMult6->GetBinContent(1);
 }

 //number of events
 Int_t nEvents6=0;
 if(fAvMult6)
 {
  nEvents6=(Int_t)(fAvMult6->GetBinEntries(1));
 }
 
 Double_t AvG6[fgkPmax6][fgkQmax6]={{0.}};  
 Bool_t someAvGEntryIsNegative6=kFALSE;   
 for(Int_t p=0;p<fgkPmax6;p++)
 {
  for(Int_t q=0;q<fgkQmax6;q++)
  {
   AvG6[p][q]=fIntGenFun6->GetBinContent(p+1,q+1);
   if(AvG6[p][q]<0.)
   {
    someAvGEntryIsNegative6=kTRUE;
   } 
  }
 }  
 /////////////////////////////////////////////////////////////////////////////      
 //////////////////gen. function for the cumulants////////////////////////////
 /////////////////////////////////////////////////////////////////////////////
  
 Double_t C6[fgkPmax6][fgkQmax6]={{0.}};//C[p][q]
 if(AvM6>0 && someAvGEntryIsNegative6==kFALSE)
 {
  for (Int_t p=0;p<fgkPmax6;p++)
  {
   for (Int_t q=0;q<fgkQmax6;q++)
   {
    C6[p][q]=1.*AvM6*(pow(AvG6[p][q],(1./AvM6))-1.); 
   }
  }
 }
 
 /////////////////////////////////////////////////////////////////////////////
 ///////avaraging the gen. function for the cumulants over azimuth////////////
 /////////////////////////////////////////////////////////////////////////////
  
 Double_t AvC6[fgkPmax6]={0.};//<C[p][q]>
 for (Int_t p=0;p<fgkPmax6;p++){
  Double_t tempHere6=0.; 
  for (Int_t q=0;q<fgkQmax6;q++){
   tempHere6+=1.*C6[p][q];
  } 
  AvC6[p]=1.*tempHere6/fgkQmax6;
 }
 
 /////////////////////////////////////////////////////////////////////////////
 //////////////////////////////////final results//////////////////////////////
 /////////////////////////////////////////////////////////////////////////////
 
 Double_t cumulant6[fgkPmax6];//array to store various order cumulants
 cumulant6[0] = (1./(fR0*fR0))*(3.*AvC[0]-(3./2.)*AvC[1]+(1./3.)*AvC[2]);
 cumulant6[1] = (2./pow(fR0,4.))*((-5.)*AvC[0]+4.*AvC[1]-1.*AvC[2]);
 cumulant6[2] = (6./pow(fR0,6.))*(3.*AvC[0]-3.*AvC[1]+1.*AvC[2]);
    
 /*
 cout<<endl;
 cout<<"*********************************"<<endl;
 cout<<"cumulants:"<<endl; 
 cout<<" c_"<<fgkFlow<<"{2} = "<<cumulant6[0]<<endl; 
 cout<<" c_"<<fgkFlow<<"{4} = "<<cumulant6[1]<<endl;
 cout<<" c_"<<fgkFlow<<"{6} = "<<cumulant6[2]<<endl;
 cout<<endl;
 */
 
 Double_t V2o6=0.,V4o6=0.,V6o6=0.;
 
 if(cumulant6[0]>=0.)
 {
  V2o6 = pow(cumulant6[0],(1./2.));
 }
 if(cumulant6[1]<=0.)
 {
  V4o6 = pow(-cumulant6[1],(1./4.));
 }
 if(cumulant6[2]>=0.)
 {
  V6o6 = pow((1./4.)*cumulant6[2],(1./6.));
 }
 
 cout<<endl;
 cout<<"***********************************"<<endl;
 cout<<"***********************************"<<endl;
 cout<<"flow estimates from GF-cumulants:"<<endl;
 cout<<"  (calculated up to 6th order)   "<<endl;
 cout<<endl;
 
 Double_t SdQo6[3]={0.};
 Double_t ChiQo6[3]={0.};
          
 //v_2{2}
 if(nEvents6 && AvM6 && cumulant6[0]>=0. && (AvQ2x+AvQ2y-pow(AvQx,2.)-pow(AvQy,2.)-pow(pow(cumulant6[0],(1./2.))*AvM6,2.)>0.))
 {        
  ChiQo6[0]=AvM6*V2o6/pow(AvQ2x+AvQ2y-pow(AvQx,2.)-pow(AvQy,2.)-pow(V2o6*AvM6,2.),0.5);
  if(ChiQo6[0])
  {
   SdQo6[0]=pow(((1./(2.*AvM6*nEvents6))*((1.+2.*pow(ChiQo6[0],2))/(2.*pow(ChiQo6[0],2)))),0.5);
  }
  cout<<" v_"<<fgkFlow<<"{2} = "<<V2o6<<" +/- "<<SdQo6[0]<<endl;
  //cout<<" v_"<<fgkFlow<<"{2} = "<<V2o6<<" +/- "<<SdQo6[0]<<", chi{2} = "<<ChiQo6[0]<<endl;//printing also the chi
 } 
 else 
 {
  cout<<" v_"<<fgkFlow<<"{2} = Im"<<endl; 
 }
   
 //v_2{4}   
 if(nEvents6 && AvM6 && cumulant6[1]<=0. && (AvQ2x+AvQ2y-pow(AvQx,2.)-pow(AvQy,2.)-pow(pow(-cumulant6[1],(1./4.))*AvM6,2.)>0.))
 {
  ChiQo6[1]=AvM6*V4o6/pow(AvQ2x+AvQ2y-pow(AvQx,2.)-pow(AvQy,2.)-pow(V4o6*AvM6,2.),0.5);
  if(ChiQo6[1])
  {
   SdQo6[1]=(1./(pow(2.*AvM6*nEvents6,0.5)))*pow((1.+4.*pow(ChiQo6[1],2)+1.*pow(ChiQo6[1],4.)+2.*pow(ChiQo6[1],6.))/(2.*pow(ChiQo6[1],6.)),0.5);
  }
  cout<<" v_"<<fgkFlow<<"{4} = "<<V4o6<<" +/- "<<SdQo6[1]<<endl;    
  //cout<<" v_"<<fgkFlow<<"{4} = "<<V4o6<<" +/- "<<SdQo6[1]<<", chi{4} = "<<ChiQo6[1]<<endl;//printing also the chi 
 } 
 else 
 {
  cout<<" v_"<<fgkFlow<<"{4} = Im"<<endl;  
 }   
  
 //v_2{6}
 if(nEvents6 && AvM6 && cumulant6[2]>=0. && (AvQ2x+AvQ2y-pow(AvQx,2.)-pow(AvQy,2.)-pow(pow((1./4.)*cumulant6[2],(1./6.))*AvM6,2.)>0.))
 {
  ChiQo6[2]=AvM6*V6/pow(AvQ2x+AvQ2y-pow(AvQx,2.)-pow(AvQy,2.)-pow(V6o6*AvM6,2.),0.5);
  if(ChiQo6[2])
  {
   SdQo6[2]=(1./(pow(2.*AvM6*nEvents6,0.5)))*pow((3.+18.*pow(ChiQo6[2],2.)+9.*pow(ChiQo6[2],4.)+28.*pow(ChiQo6[2],6.)+12.*pow(ChiQo6[2],8.)+24.*pow(ChiQo6[2],10.))/(24.*pow(ChiQo6[2],10.)),0.5);
  } 
   cout<<" v_"<<fgkFlow<<"{6} = "<<V6o6<<" +/- "<<SdQo6[2]<<endl;   
   //cout<<" v_"<<fgkFlow<<"{6} = "<<V6o6<<" +/- "<<SdQo6[2]<<", chi{6} = "<<ChiQo6[2]<<endl;//printing also the chi
 }
 else
 {
  cout<<" v_"<<fgkFlow<<"{6} = Im"<<endl;  
 }
 
 cout<<endl;
 cout<<"   nEvts = "<<nEvents6<<", AvM = "<<AvM6<<endl; 
 cout<<"***********************************"<<endl;    
 cout<<"***********************************"<<endl;   
 
 //============================================================================================================================================== 
    
 //up to 8th order
 //avarage selected multiplicity
 Double_t AvM8=0.;
 if(fAvMult8)
 {
  AvM8=fAvMult8->GetBinContent(1);
 }

 //number of events
 Int_t nEvents8=0;
 if(fAvMult8)
 {
  nEvents8=(Int_t)(fAvMult8->GetBinEntries(1));
 }
 
 Double_t AvG8[fgkPmax8][fgkQmax8]={{0.}}; 
 Bool_t someAvGEntryIsNegative8=kFALSE;   
 for(Int_t p=0;p<fgkPmax8;p++)
 {
  for(Int_t q=0;q<fgkQmax8;q++)
  {
   AvG8[p][q]=fIntGenFun8->GetBinContent(p+1,q+1);
   if(AvG8[p][q]<0.)
   {
    someAvGEntryIsNegative8=kTRUE;
   } 
  }
 }  
 /////////////////////////////////////////////////////////////////////////////      
 //////////////////gen. function for the cumulants////////////////////////////
 /////////////////////////////////////////////////////////////////////////////
  
 Double_t C8[fgkPmax8][fgkQmax8]={{0.}};//C[p][q]
 if(AvM8>0 && someAvGEntryIsNegative8==kFALSE)
 {
  for (Int_t p=0;p<fgkPmax8;p++)
  {
   for (Int_t q=0;q<fgkQmax8;q++)
   {
    C8[p][q]=1.*AvM8*(pow(AvG8[p][q],(1./AvM8))-1.); 
   }
  }
 }
 
 /////////////////////////////////////////////////////////////////////////////
 ///////avaraging the gen. function for the cumulants over azimuth////////////
 /////////////////////////////////////////////////////////////////////////////
  
 Double_t AvC8[fgkPmax8]={0.};//<C[p][q]>
 for (Int_t p=0;p<fgkPmax8;p++)
 {
  Double_t tempHere8=0.; 
  for (Int_t q=0;q<fgkQmax8;q++)
  {
   tempHere8+=1.*C8[p][q];
  } 
  AvC8[p]=1.*tempHere8/fgkQmax8;
 }
 
 /////////////////////////////////////////////////////////////////////////////
 //////////////////////////////////final results//////////////////////////////
 /////////////////////////////////////////////////////////////////////////////
 
 Double_t cumulant8[fgkPmax8];//array to store various order cumulants
 cumulant8[0] = (1./(fR0*fR0))*(4.*AvC[0]-3.*AvC[1]+(4./3.)*AvC[2]-(1./4.)*AvC[3]);
 cumulant8[1] = (1./pow(fR0,4.))*((-52./3.)*AvC[0]+19.*AvC[1]-(28./3.)*AvC[2]+(11./6.)*AvC[3]);
 cumulant8[2] = (3./pow(fR0,6.))*(18.*AvC[0]-24.*AvC[1]+14.*AvC[2]-3.*AvC[3]);
 cumulant8[3] = (24./pow(fR0,8.))*((-4.)*AvC[0]+6.*AvC[1]-4.*AvC[2]+1.*AvC[3]);
  
 /* x0,y0 ~ 1/sqrt(p) (setting another complex mesh, doesn't work correctly)
 cumulant8[0] = (-1./(6.*fR0*fR0)) * (1.*AvC[0] - 48.*AvC[1] + 243.*AvC[2] - 256.*AvC[3]);
 cumulant8[1] = (2./pow(fR0,4.)) * (3.*AvC[0] - 128.*AvC[1] + 567.*AvC[2] - 512.*AvC[3]);
 cumulant8[2] = (-12./pow(fR0,6.)) * (13.*AvC[0] - 456.*AvC[1] + 1701.*AvC[2] - 1408.*AvC[3]);
 cumulant8[3] = (2304./pow(fR0,8.)) * (1.*AvC[0] - 24.*AvC[1] + 81.*AvC[2] - 64.*AvC[3]);       
 */
 
 /* x0,y0 ~ p (setting another complex mesh, doesn't work correctly)
 cumulant8[0] = (-1./(5040.*fR0*fR0)) * ((-8064.)*AvC[0] + 1008.*AvC[1] - 128.*AvC[2] + 9.*AvC[3]);
 cumulant8[1] = (1./(720.*pow(fR0,4.))) * (1952.*AvC[0] - 676.*AvC[1] + 96.*AvC[2] - 7.*AvC[3]);
 cumulant8[2] = (-1./(40.*pow(fR0,6.))) * ((-116.)*AvC[0] + 52.*AvC[1] - 12.*AvC[2] + 1.*AvC[3]);
 cumulant8[3] = (1./(35.*pow(fR0,8.))) * (56.*AvC[0] - 28.*AvC[1] + 8.*AvC[2] - 1.*AvC[3]);       
 */ 
                                           
 /*                                          
 cout<<endl;
 cout<<"*********************************"<<endl;
 cout<<"cumulants8:"<<endl; 
 cout<<" c_"<<fgkFlow<<"{2} = "<<cumulant8[0]<<endl; 
 cout<<" c_"<<fgkFlow<<"{4} = "<<cumulant8[1]<<endl;
 cout<<" c_"<<fgkFlow<<"{6} = "<<cumulant8[2]<<endl;
 cout<<" c_"<<fgkFlow<<"{8} = "<<cumulant8[3]<<endl; 
 cout<<endl;
 */
 
 Double_t V2o8=0.,V4o8=0.,V6o8=0.,V8o8=0.;
 
 if(cumulant8[0]>=0.)
 {
  V2o8 = pow(cumulant8[0],(1./2.));
 }
 if(cumulant8[1]<=0.)
 {
  V4o8 = pow(-cumulant8[1],(1./4.));
 }
 if(cumulant8[2]>=0.)
 {
  V6o8 = pow((1./4.)*cumulant8[2],(1./6.));
 }
 if(cumulant8[3]<=0.)
 {
  V8o8 = pow(-(1./33.)*cumulant8[3],(1./8.));
 }
 
 cout<<endl;
 cout<<"***********************************"<<endl;
 cout<<"***********************************"<<endl;
 cout<<"flow estimates from GF-cumulants:"<<endl;
 cout<<"  (calculated up to 8th order)   "<<endl;
 cout<<endl;   
             
 Double_t SdQo8[4]={0.};
 Double_t ChiQo8[4]={0.};
          
 //v_2{2}
 if(nEvents8 && AvM8 && cumulant8[0]>=0. && (AvQ2x+AvQ2y-pow(AvQx,2.)-pow(AvQy,2.)-pow(pow(cumulant8[0],(1./2.))*AvM8,2.)>0.))
 {        
  ChiQo8[0]=AvM8*V2o8/pow(AvQ2x+AvQ2y-pow(AvQx,2.)-pow(AvQy,2.)-pow(V2o8*AvM8,2.),0.5);
  if(ChiQo8[0])
  {
   SdQo8[0]=pow(((1./(2.*AvM8*nEvents8))*((1.+2.*pow(ChiQo8[0],2.))/(2.*pow(ChiQo8[0],2)))),0.5);
  }
  cout<<" v_"<<fgkFlow<<"{2} = "<<V2o8<<" +/- "<<SdQo8[0]<<endl;    
  //cout<<" v_"<<fgkFlow<<"{2} = "<<V2o8<<" +/- "<<SdQo8[0]<<", chi{2} = "<<ChiQo8[0]<<endl;//printing also the chi
 } 
 else 
 {
  cout<<" v_"<<fgkFlow<<"{2} = Im"<<endl; 
 }
   
 //v_2{4}   
 if(nEvents8 && AvM8 && cumulant8[1]<=0. && (AvQ2x+AvQ2y-pow(AvQx,2.)-pow(AvQy,2.)-pow(pow(-cumulant8[1],(1./4.))*AvM8,2.)>0.))
 {
  ChiQo8[1]=AvM8*V4o8/pow(AvQ2x+AvQ2y-pow(AvQx,2.)-pow(AvQy,2.)-pow(V4o8*AvM8,2.),0.5);
  if(ChiQo8[1])
  {
   SdQo8[1]=(1./(pow(2.*AvM8*nEvents8,0.5)))*pow((1.+4.*pow(ChiQo8[1],2)+1.*pow(ChiQo8[1],4.)+2.*pow(ChiQo8[1],6.))/(2.*pow(ChiQo8[1],6.)),0.5);
  }
  cout<<" v_"<<fgkFlow<<"{4} = "<<V4o8<<" +/- "<<SdQo8[1]<<endl;    
  //cout<<" v_"<<fgkFlow<<"{4} = "<<V4o8<<" +/- "<<SdQo8[1]<<", chi{4} = "<<ChiQo8[1]<<endl;//printing also the chi 
 } 
 else 
 {
  cout<<" v_"<<fgkFlow<<"{4} = Im"<<endl;  
 } 
  
 //v_2{6}
 if(nEvents8 && AvM8 && cumulant8[2]>=0. && (AvQ2x+AvQ2y-pow(AvQx,2.)-pow(AvQy,2.)-pow(pow((1./4.)*cumulant8[2],(1./6.))*AvM8,2.)>0.))
 {
  ChiQo8[2]=AvM8*V6o8/pow(AvQ2x+AvQ2y-pow(AvQx,2.)-pow(AvQy,2.)-pow(V6o8*AvM8,2.),0.5);
  if(ChiQo8[2])
  {
   SdQo8[2]=(1./(pow(2.*AvM8*nEvents8,0.5)))*pow((3.+18.*pow(ChiQo8[2],2)+9.*pow(ChiQo8[2],4.)+28.*pow(ChiQo8[2],6.)+12.*pow(ChiQo8[2],8.)+24.*pow(ChiQo8[2],10.))/(24.*pow(ChiQo8[2],10.)),0.5);
  }
  cout<<" v_"<<fgkFlow<<"{6} = "<<V6o8<<" +/- "<<SdQo8[2]<<endl;
  //cout<<" v_"<<fgkFlow<<"{6} = "<<V6o8<<" +/- "<<SdQo8[2]<<", chi{6} = "<<ChiQo8[2]<<endl;//printing also the chi 
 }
 else
 {
  cout<<" v_"<<fgkFlow<<"{6} = Im"<<endl;  
 }
  
 //v_2{8}
 if(nEvents8 && AvM8 && cumulant8[3]<=0. && (AvQ2x+AvQ2y-pow(AvQx,2.)-pow(AvQy,2.)-pow(pow(-(1./33.)*cumulant8[3],(1./8.))*AvM8,2.)>0.))
 {  
  ChiQo8[3]=AvM8*V8o8/pow(AvQ2x+AvQ2y-pow(AvQx,2.)-pow(AvQy,2.)-pow(V8o8*AvM8,2.),0.5);
  if(ChiQo8[3])
  {
   SdQo8[3]=(1./(pow(2.*AvM8*nEvents8,0.5)))*pow((12.+96.*pow(ChiQo8[3],2)+72.*pow(ChiQo8[3],4.)+304.*pow(ChiQo8[3],6.)+257.*pow(ChiQo8[3],8.)+804.*pow(ChiQo8[3],10.)+363.*pow(ChiQo8[3],12.)+726.*pow(ChiQo8[3],14.))/(726.*pow(ChiQo8[3],14.)),0.5);
  } 
  cout<<" v_"<<fgkFlow<<"{8} = "<<V8o8<<" +/- "<<SdQo8[3]<<endl;  
  //cout<<" v_"<<fgkFlow<<"{8} = "<<V8o8<<" +/- "<<SdQo8[3]<<", chi{8} = "<<ChiQo8[3]<<endl;//printing also the chi 
 } 
 else 
 {
  cout<<" v_"<<fgkFlow<<"{8} = Im"<<endl;     
 }
  
 cout<<endl;
 cout<<"   nEvts = "<<nEvents8<<", AvM = "<<AvM8<<endl; 
 cout<<"*********************************"<<endl; 

 //============================================================================================================================================== 
   
 //up to 16-th order: 
 //avarage selected multiplicity
 Double_t AvM16=0.;
 if(fAvMult16)
 {
  AvM16=fAvMult16->GetBinContent(1);
 }

 //number of events
 Int_t nEvents16=0;
 if(fAvMult16)
 {
  nEvents16=(Int_t)(fAvMult16->GetBinEntries(1));
 }
 
 Double_t AvG16[fgkPmax16][fgkQmax16]={{0.}};  
 Bool_t someAvGEntryIsNegative16=kFALSE; 
 for(Int_t p=0;p<fgkPmax16;p++)
 {
  for(Int_t q=0;q<fgkQmax16;q++)
  {
   AvG16[p][q]=fIntGenFun16->GetBinContent(p+1,q+1);
   if(AvG16[p][q]<0.)
   {
    someAvGEntryIsNegative16=kTRUE;
   } 
  }
 }  
 
 /////////////////////////////////////////////////////////////////////////////      
 //////////////////gen. function for the cumulants////////////////////////////
 /////////////////////////////////////////////////////////////////////////////
  
 Double_t C16[fgkPmax16][fgkQmax16]={{0.}};//C16[p][q]
 if(AvM16>0 && someAvGEntryIsNegative16==kFALSE)
 {
  for(Int_t p=0;p<fgkPmax16;p++)
  {
   for(Int_t q=0;q<fgkQmax16;q++)
   {
    C16[p][q]=1.*AvM16*(pow(AvG16[p][q],(1./AvM16))-1.); 
   }
  }
 }
 
 /////////////////////////////////////////////////////////////////////////////
 ///////avaraging the gen. function for the cumulants over azimuth////////////
 /////////////////////////////////////////////////////////////////////////////
  
 Double_t AvC16[fgkPmax16]={0.};//<C16[p][q]>
 for (Int_t p=0;p<fgkPmax16;p++)
 {
  Double_t tempHere16=0.; 
  for (Int_t q=0;q<fgkQmax16;q++)
  {
   tempHere16+=1.*C16[p][q];
  } 
  AvC16[p]=1.*tempHere16/fgkQmax16;
 }
 
 /////////////////////////////////////////////////////////////////////////////
 //////////////////////////////////final results//////////////////////////////
 /////////////////////////////////////////////////////////////////////////////
 
 Double_t cumulant16[fgkPmax16];//array to store various order cumulants
  
 cumulant16[0] = (1./(fR0*fR0)) * (8.*AvC16[0] - 14.*AvC16[1] + (56./3.)*AvC16[2] - (35./2.)*AvC16[3] + 
			      (56./5.)*AvC16[4] - (14./3.)*AvC16[5] + (8./7.)*AvC16[6] - (1./8.)*AvC16[7]);

 cumulant16[1] = (1./pow(fR0,4.)) * ((-1924./35.)*AvC16[0] + (621./5.)*AvC16[1] - (8012./45.)*AvC16[2] + 
				 (691./4.)*AvC16[3] - (564./5.)*AvC16[4] + (2143./45.)*AvC16[5] - 
				 (412./35.)*AvC16[6] + (363./280.)*AvC16[7]);

 cumulant16[2] = (1./pow(fR0,6.)) * (349.*AvC16[0] - (18353./20.)*AvC16[1] + (7173./5.)*AvC16[2] - 
				 1457.*AvC16[3] + (4891./5.)*AvC16[4] - (1683./4.)*AvC16[5] + 
				 (527./5.)*AvC16[6] - (469./40.)*AvC16[7]);

 cumulant16[3] = (1./pow(fR0,8.)) * ((-10528./5.)*AvC16[0] + (30578./5.)*AvC16[1] - (51456./5.)*AvC16[2] + 
				 10993.*AvC16[3] - (38176./5.)*AvC16[4] + (16818./5.)*AvC16[5] - 
				 (4288./5.)*AvC16[6] + (967./10.)*AvC16[7]);

 cumulant16[4] = (1./pow(fR0,10.)) * (11500.*AvC16[0] - 35800.*AvC16[1] + 63900.*AvC16[2] - 71600.*AvC16[3] + 
				  51620.*AvC16[4] - 23400.*AvC16[5] + 6100.*AvC16[6] - 700.*AvC16[7]);

 cumulant16[5] = (1./pow(fR0,12.)) * (-52560.*AvC16[0] + 172080.*AvC16[1] - 321840.*AvC16[2] + 376200.*AvC16[3] - 
				  281520.*AvC16[4] + 131760.*AvC16[5] - 35280.*AvC16[6] + 4140.*AvC16[7]);

 cumulant16[6] = (1./pow(fR0,14.)) * (176400.*AvC16[0] - 599760.*AvC16[1] + 1164240.*AvC16[2] - 1411200.*AvC16[3] + 
				  1093680.*AvC16[4] - 529200.*AvC16[5] + 146160.*AvC16[6] - 17640.*AvC16[7]);

 cumulant16[7] = (1./pow(fR0,16.)) * (-322560*AvC16[0] + 1128960.*AvC16[1] - 2257920.*AvC16[2] + 2822400.*AvC16[3] - 
				  2257920.*AvC16[4] + 1128960.*AvC16[5] - 322560.*AvC16[6] + 40320.*AvC16[7]);
    
 /*
 cout<<endl;
 cout<<"*********************************"<<endl;
 cout<<"cumulants:"<<endl;
 cout<<" c_"<<fgkFlow<<"{2} = "<<cumulant16[0]<<endl; 
 cout<<" c_"<<fgkFlow<<"{4} = "<<cumulant16[1]<<endl;
 cout<<" c_"<<fgkFlow<<"{6} = "<<cumulant16[2]<<endl;
 cout<<" c_"<<fgkFlow<<"{8} = "<<cumulant16[3]<<endl; 
 cout<<"c_"<<fgkFlow<<"{10} = "<<cumulant16[4]<<endl; 
 cout<<"c_"<<fgkFlow<<"{12} = "<<cumulant16[5]<<endl;
 cout<<"c_"<<fgkFlow<<"{14} = "<<cumulant16[6]<<endl; 
 cout<<"c_"<<fgkFlow<<"{16} = "<<cumulant16[7]<<endl; 
 cout<<endl;
 */
 
 Double_t V2o16=0.,V4o16=0.,V6o16=0.,V8o16=0.,V10o16=0.,V12o16=0.,V14o16=0.,V16o16=0.;
 
 if(cumulant16[0]>=0.)
 {
  V2o16 = pow(cumulant16[0],(1./2.));
 }
 if(cumulant16[1]<=0.)
 {
  V4o16 = pow(-cumulant16[1],(1./4.));
 }
 if(cumulant16[2]>=0.)
 {
  V6o16 = pow((1./4.)*cumulant16[2],(1./6.));
 }
 if(cumulant16[3]<=0.)
 {
  V8o16 = pow(-(1./33.)*cumulant16[3],(1./8.));
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
          
 Double_t SdQo16[8]={0.};
 Double_t ChiQo16[8]={0.};
          
 //v_2{2}
 if(nEvents16 && AvM16 && cumulant16[0]>=0. && (AvQ2x+AvQ2y-pow(AvQx,2.)-pow(AvQy,2.)-pow(pow(cumulant16[0],(1./2.))*AvM16,2.)>0.))
 {        
  ChiQo16[0]=AvM16*V2o16/pow(AvQ2x+AvQ2y-pow(AvQx,2.)-pow(AvQy,2.)-pow(V2o16*AvM16,2.),0.5);
  if(ChiQo16[0])
  {
   SdQo16[0]=pow(((1./(2.*AvM16*nEvents16))*((1.+2.*pow(ChiQo16[0],2))/(2.*pow(ChiQo16[0],2)))),0.5);
  }
  cout<<" v_"<<fgkFlow<<"{2} = "<<V2o16<<" +/- "<<SdQo16[0]<<endl;
  //cout<<" v_"<<fgkFlow<<"{2} = "<<V2o16<<" +/- "<<SdQo16[0]<<", chi{2} = "<<ChiQo16[0]<<endl;//printing also the chi
 } 
 else 
 {
  cout<<" v_"<<fgkFlow<<"{2} = Im"<<endl; 
 }
   
 //v_2{4}   
 if(nEvents16 && AvM16 && cumulant16[1]<=0. && (AvQ2x+AvQ2y-pow(AvQx,2.)-pow(AvQy,2.)-pow(pow(-cumulant16[1],(1./4.))*AvM16,2.)>0.))
 {
  ChiQo16[1]=AvM16*V4o16/pow(AvQ2x+AvQ2y-pow(AvQx,2.)-pow(AvQy,2.)-pow(V4o16*AvM16,2.),0.5);
  if(ChiQo16[1])
  {
   SdQo16[1]=(1./(pow(2.*AvM16*nEvents16,0.5)))*pow((1.+4.*pow(ChiQo16[1],2.)+1.*pow(ChiQo16[1],4.)+2.*pow(ChiQo16[1],6.))/(2.*pow(ChiQo16[1],6.)),0.5);
  }
  cout<<" v_"<<fgkFlow<<"{4} = "<<V4o16<<" +/- "<<SdQo16[1]<<endl; 
  //cout<<" v_"<<fgkFlow<<"{4} = "<<V4o16<<" +/- "<<SdQo16[1]<<", chi{4} = "<<ChiQo16[1]<<endl;//printing also the chi
 } 
 else 
 {
  cout<<" v_"<<fgkFlow<<"{4} = Im"<<endl;  
 } 
  
 //v_2{6}
 if(nEvents16 && AvM16 && cumulant16[2]>=0. && (AvQ2x+AvQ2y-pow(AvQx,2.)-pow(AvQy,2.)-pow(pow((1./4.)*cumulant16[2],(1./6.))*AvM16,2.)>0.))
 {
  ChiQo16[2]=AvM16*V6o16/pow(AvQ2x+AvQ2y-pow(AvQx,2.)-pow(AvQy,2.)-pow(V6o16*AvM16,2.),0.5);
  if(ChiQo16[2])
  {
   SdQo16[2]=(1./(pow(2.*AvM16*nEvents16,0.5)))*pow((3.+18.*pow(ChiQo16[2],2)+9.*pow(ChiQo16[2],4.)+28.*pow(ChiQo16[2],6.)+12.*pow(ChiQo16[2],8.)+24.*pow(ChiQo16[2],10.))/(24.*pow(ChiQo16[2],10.)),0.5);
  }
  cout<<" v_"<<fgkFlow<<"{6} = "<<V6o16<<" +/- "<<SdQo16[2]<<endl;
  //cout<<" v_"<<fgkFlow<<"{6} = "<<V6o16<<" +/- "<<SdQo16[2]<<", chi{6} = "<<ChiQo16[2]<<endl;//printing also the chi
 }
 else
 {
  cout<<" v_"<<fgkFlow<<"{6} = Im"<<endl;  
 }
  
 //v_2{8}
 if(nEvents16 && AvM16 && cumulant16[3]<=0. && (AvQ2x+AvQ2y-pow(AvQx,2.)-pow(AvQy,2.)-pow(pow(-(1./33.)*cumulant16[3],(1./8.))*AvM16,2.)>0.))
 {  
  ChiQo16[3]=AvM16*V8o16/pow(AvQ2x+AvQ2y-pow(AvQx,2.)-pow(AvQy,2.)-pow(V8o16*AvM16,2.),0.5);
  if(ChiQo16[3])
  {
   SdQo16[3]=(1./(pow(2.*AvM16*nEvents16,0.5)))*pow((12.+96.*pow(ChiQo16[3],2)+72.*pow(ChiQo16[3],4.)+304.*pow(ChiQo16[3],6.)+257.*pow(ChiQo16[3],8.)+804.*pow(ChiQo16[3],10.)+363.*pow(ChiQo16[3],12.)+726.*pow(ChiQo16[3],14.))/(726.*pow(ChiQo16[3],14.)),0.5);
  } 
  cout<<" v_"<<fgkFlow<<"{8} = "<<V8o16<<" +/- "<<SdQo16[3]<<endl;
  //cout<<" v_"<<fgkFlow<<"{8} = "<<V8o16<<" +/- "<<SdQo16[3]<<", chi{8} = "<<ChiQo16[3]<<endl;//printing also the chi
 } 
 else 
 {
  cout<<" v_"<<fgkFlow<<"{8} = Im"<<endl;     
 }
  
 //v_2{10}
 if(nEvents16 && AvM16 && cumulant16[4]>=0.)
 {
  cout<<"v_"<<fgkFlow<<"{10} = "<<V10o16<<endl;
 } 
 else 
 {
  cout<<"v_"<<fgkFlow<<"{10} = Im"<<endl; 
 }
  
 //v_2{12}
 if(nEvents16 && AvM16 && AvM16 && cumulant16[5]<=0.)
 {
  cout<<"v_"<<fgkFlow<<"{12} = "<<V12o16<<endl;
 } 
 else 
 {
  cout<<"v_"<<fgkFlow<<"{12} = Im"<<endl; 
 }
  
 //v_2{14}
 if(nEvents16 && AvM16 && cumulant16[6]>=0.)
 {
  cout<<"v_"<<fgkFlow<<"{14} = "<<V14o16<<endl;
 } 
 else 
 {
  cout<<"v_"<<fgkFlow<<"{14} = Im"<<endl;  
 }
  
 //v_2{16}
 if(nEvents16 && AvM16 && cumulant16[7]<=0.)
 {
  cout<<"v_"<<fgkFlow<<"{16} = "<<V16o16<<endl;
 } 
 else 
 {
  cout<<"v_"<<fgkFlow<<"{16} = Im"<<endl;  
 }
  
 cout<<endl;
 cout<<"   nEvts = "<<nEvents16<<", AvM = "<<AvM16<<endl; 
 cout<<"*********************************"<<endl;    
 
 //==============================================================================================================================================
 
 }//end of if(otherEquations) 
}

  

