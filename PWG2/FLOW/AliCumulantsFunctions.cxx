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


//, TProfile2D *dRe0, TProfile2D *dRe1, TProfile2D *dRe2, TProfile2D *dRe3, TProfile2D *dRe4, TProfile2D *dRe5, TProfile2D *dRe6, TProfile2D *dRe7, TProfile2D *dIm0, TProfile2D *dIm1, TProfile2D *dIm2, TProfile2D *dIm3, TProfile2D *dIm4, TProfile2D *dIm5, TProfile2D *dIm6, TProfile2D *dIm7


AliCumulantsFunctions::AliCumulantsFunctions(TProfile2D *IntGenFun, TProfile3D *DiffGenFunRe, TProfile3D *DiffGenFunIm, TProfile *BinNoOfParticles, TH1D *ifr, TH1D *dfr2, TH1D *dfr4, TH1D *dfr6, TH1D *dfr8, TProfile *AvMult, TProfile *QVector, AliFlowCommonHistResults *chr2nd, AliFlowCommonHistResults *chr4th, AliFlowCommonHistResults *chr6th, AliFlowCommonHistResults *chr8th):
 fIntGenFun(IntGenFun),
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
 static const Int_t fgkQmax=AliFlowCumuConstants::kQmax;   //needed for numerics
 static const Int_t fgkPmax=AliFlowCumuConstants::kPmax;   //needed for numerics  
 static const Int_t fgkFlow=AliFlowCumuConstants::kFlow;   //integrated flow coefficient to be calculated
 static const Int_t fgkMltpl=AliFlowCumuConstants::kMltpl; //the multiple in p=m*n (diff. flow) 
 static const Int_t fgknBins=100;                          //number of pt bins //to be improved
 
 Double_t fR0=AliFlowCumuConstants::fgR0;                  //needed for numerics
 //Double_t fPtMax=AliFlowCommonConstants::GetPtMax();       //maximum pt
 //Double_t fPtMin=AliFlowCommonConstants::GetPtMin();       //minimum pt
 //Double_t fBinWidth=(fPtMax-fPtMin)/fgknBins;              //width of pt bin (in GeV)   
     
 //<G[p][q]>
 Double_t AvG[fgkPmax][fgkQmax]={{0.}};                                  
 for(Int_t p=0;p<fgkPmax;p++){
  for(Int_t q=0;q<fgkQmax;q++){ 
   AvG[p][q]=fIntGenFun->GetBinContent(p+1,q+1);
  }
 } 

 //avarage selected multiplicity
 Double_t AvM = fAvMult->GetBinContent(1);
 
 //number of events
 Int_t nEvents = (Int_t)(fAvMult->GetBinEntries(1));
 
 //<Q-vector stuff>
 Double_t AvQx  = fQVector->GetBinContent(1); //<Q_x>
 Double_t AvQy  = fQVector->GetBinContent(2); //<Q_y>
 Double_t AvQ2x = fQVector->GetBinContent(3); //<(Q_x)^2>
 Double_t AvQ2y = fQVector->GetBinContent(4); //<(Q_y)^2>
  
 /////////////////////////////////////////////////////////////////////////////      
 //////////////////gen. function for the cumulants////////////////////////////
 /////////////////////////////////////////////////////////////////////////////
  
 Double_t C[fgkPmax][fgkQmax]={{0.}};//C[p][q]
 if(AvM!=0){
  for (Int_t p=0;p<fgkPmax;p++){
   for (Int_t q=0;q<fgkQmax;q++){
    C[p][q]=1.*AvM*(pow(AvG[p][q],(1./AvM))-1.); 
   }
  }
 }
 /////////////////////////////////////////////////////////////////////////////
 ///////avaraging the gen. function for the cumulants over azimuth////////////
 /////////////////////////////////////////////////////////////////////////////
  
 Double_t AvC[fgkPmax]={0.};//<C[p][q]>
 for (Int_t p=0;p<fgkPmax;p++){
  Double_t tempHere=0.; 
  for (Int_t q=0;q<fgkQmax;q++){
   tempHere+=1.*C[p][q];
  } 
  AvC[p]=1.*tempHere/fgkQmax;
 }
 
 /////////////////////////////////////////////////////////////////////////////
 //////////////////////////////////final results//////////////////////////////
 /////////////////////////////////////////////////////////////////////////////
 
 Double_t cumulant[fgkPmax];//array to store various order cumulants
  
 cumulant[0] = (1./(fR0*fR0)) * (8.*AvC[0] - 14.*AvC[1] + (56./3.)*AvC[2] - (35./2.)*AvC[3] + 
			      (56./5.)*AvC[4] - (14./3.)*AvC[5] + (8./7.)*AvC[6] - (1./8.)*AvC[7]);

 cumulant[1] = (1./pow(fR0,4.)) * ((-1924./35.)*AvC[0] + (621./5.)*AvC[1] - (8012./45.)*AvC[2] + 
				 (691./4.)*AvC[3] - (564./5.)*AvC[4] + (2143./45.)*AvC[5] - 
				 (412./35.)*AvC[6] + (363./280.)*AvC[7]);

 cumulant[2] = (1./pow(fR0,6.)) * (349.*AvC[0] - (18353./20.)*AvC[1] + (7173./5.)*AvC[2] - 
				 1457.*AvC[3] + (4891./5.)*AvC[4] - (1683./4.)*AvC[5] + 
				 (527./5.)*AvC[6] - (469./40.)*AvC[7]);

 cumulant[3] = (1./pow(fR0,8.)) * ((-10528./5.)*AvC[0] + (30578./5.)*AvC[1] - (51456./5.)*AvC[2] + 
				 10993.*AvC[3] - (38176./5.)*AvC[4] + (16818./5.)*AvC[5] - 
				 (4288./5.)*AvC[6] + (967./10.)*AvC[7]);

 cumulant[4] = (1./pow(fR0,10.)) * (11500.*AvC[0] - 35800.*AvC[1] + 63900.*AvC[2] - 71600.*AvC[3] + 
				  51620.*AvC[4] - 23400.*AvC[5] + 6100.*AvC[6] - 700.*AvC[7]);

 cumulant[5] = (1./pow(fR0,12.)) * (-52560.*AvC[0] + 172080.*AvC[1] - 321840.*AvC[2] + 376200.*AvC[3] - 
				  281520.*AvC[4] + 131760.*AvC[5] - 35280.*AvC[6] + 4140.*AvC[7]);

 cumulant[6] = (1./pow(fR0,14.)) * (176400.*AvC[0] - 599760.*AvC[1] + 1164240.*AvC[2] - 1411200.*AvC[3] + 
				  1093680.*AvC[4] - 529200.*AvC[5] + 146160.*AvC[6] - 17640.*AvC[7]);

 cumulant[7] = (1./pow(fR0,16.)) * (-322560*AvC[0] + 1128960.*AvC[1] - 2257920.*AvC[2] + 2822400.*AvC[3] - 
				  2257920.*AvC[4] + 1128960.*AvC[5] - 322560.*AvC[6] + 40320.*AvC[7]);
    
 cout<<""<<endl;
 cout<<"*********************************"<<endl;
 cout<<"cumulants:"<<endl;
  
 cout<<" c_"<<fgkFlow<<"{2} = "<<cumulant[0]<<endl; 
 cout<<" c_"<<fgkFlow<<"{4} = "<<cumulant[1]<<endl;
 cout<<" c_"<<fgkFlow<<"{6} = "<<cumulant[2]<<endl;
 cout<<" c_"<<fgkFlow<<"{8} = "<<cumulant[3]<<endl; 
 cout<<"c_"<<fgkFlow<<"{10} = "<<cumulant[4]<<endl; 
 cout<<"c_"<<fgkFlow<<"{12} = "<<cumulant[5]<<endl;
 cout<<"c_"<<fgkFlow<<"{14} = "<<cumulant[6]<<endl; 
 cout<<"c_"<<fgkFlow<<"{16} = "<<cumulant[7]<<endl; 
  
 cout<<""<<endl;
 cout<<"integrated flow: "<<endl;
 
 Double_t V2=0.,V4=0.,V6=0.,V8=0.,V10=0.,V12=0.,V14=0.,V16=0.;
 
 if(cumulant[0]>=0.){
  V2 = pow(cumulant[0],(1./2.));
 }
 if(cumulant[1]<=0.){
  V4 = pow(-cumulant[1],(1./4.));
 }
 if(cumulant[2]>=0.){
  V6 = pow((1./4.)*cumulant[2],(1./6.));
 }
 if(cumulant[3]<=0.){
  V8 = pow(-(1./33.)*cumulant[3],(1./8.));
 }
 if(cumulant[4]>=0.){
  V10 = pow((1./456.)*cumulant[4],(1./10.));
 }
 if(cumulant[5]<=0.){
  V12 = pow(-(1./9460.)*cumulant[5],(1./12.));
 }
 if(cumulant[6]>=0.){
  V14 = pow((1./274800.)*cumulant[6],(1./14.));
 }
 if(cumulant[7]<=0.){
  V16 = pow(-(1./10643745.)*cumulant[7],(1./16.));
 }
    
 Double_t SdQ[4]={0.};
 Double_t ChiQ[4]={0.};
          
   //v_2{2}
   if(AvM!=0 && (cumulant[0]>=0.) && (AvQ2x+AvQ2y-pow(AvQx,2.)-pow(AvQy,2.)-pow(pow(cumulant[0],(1./2.))*AvM,2.)>0.))
   {        
    ChiQ[0]=AvM*V2/pow(AvQ2x+AvQ2y-pow(AvQx,2.)-pow(AvQy,2.)-pow(V2*AvM,2.),0.5);
    SdQ[0]=pow(((1./(2.*AvM*nEvents))*((1.+1.*pow(ChiQ[0],2))/(1.*pow(ChiQ[0],2)))),0.5);
    cout<<" v_"<<fgkFlow<<"{2} = "<<V2<<" +/- "<<SdQ[0]<<", chi{2} = "<<ChiQ[0]<<endl;
    fifr->SetBinContent(1,V2);
    fifr->SetBinError(1,SdQ[0]);
    //common histograms
    fchr2nd->FillIntegratedFlow(V2,SdQ[0]);
    fchr2nd->FillChi(V2*pow(AvM,0.5));
   } 
   else 
   {
    cout<<" v_"<<fgkFlow<<"{2} = Im"<<endl; 
   }
   
   //v_2{4}   
   if(AvM!=0 && (cumulant[1]<=0.) && (AvQ2x+AvQ2y-pow(AvQx,2.)-pow(AvQy,2.)-pow(pow(-cumulant[1],(1./4.))*AvM,2.)>0.))
   {
    ChiQ[1]=AvM*V4/pow(AvQ2x+AvQ2y-pow(AvQx,2.)-pow(AvQy,2.)-pow(V4*AvM,2.),0.5);
    SdQ[1]=(1./(pow(2.*AvM*nEvents,.5)))*pow((1.+2.*pow(ChiQ[1],2)+(1./4.)*pow(ChiQ[1],4.)+(1./4.)*pow(ChiQ[1],6.))/((1./4.)*pow(ChiQ[1],6.)),.5);
    cout<<" v_"<<fgkFlow<<"{4} = "<<V4<<" +/- "<<SdQ[1]<<", chi{4} = "<<ChiQ[1]<<endl;
    fifr->SetBinContent(2,V4);
    fifr->SetBinError(2,SdQ[1]);
    //common histograms
    fchr4th->FillIntegratedFlow(V4,SdQ[1]);
    fchr4th->FillChi(V4*pow(AvM,0.5));
   } 
   else 
   {
    cout<<" v_"<<fgkFlow<<"{4} = Im"<<endl;  
   } 
  
  //v_2{6}
  if(AvM!=0 && (cumulant[2]>=0.) && (AvQ2x+AvQ2y-pow(AvQx,2.)-pow(AvQy,2.)-pow(pow((1./4.)*cumulant[2],(1./6.))*AvM,2.)>0.))
  {
   ChiQ[2]=AvM*V6/pow(AvQ2x+AvQ2y-pow(AvQx,2.)-pow(AvQy,2.)-pow(V6*AvM,2.),0.5);
   SdQ[2]=(1./(pow(2.*AvM*nEvents,.5)))*pow((3.+18.*pow(ChiQ[2],2)+9.*pow(ChiQ[2],4.)+28.*pow(ChiQ[2],6.)+12.*pow(ChiQ[2],8.)+24.*pow(ChiQ[2],10.))/(24.*pow(ChiQ[2],10.)),.5);
   cout<<" v_"<<fgkFlow<<"{6} = "<<V6<<" +/- "<<SdQ[2]<<", chi{6} = "<<ChiQ[2]<<endl;
   fifr->SetBinContent(3,V6);
   fifr->SetBinError(3,SdQ[2]);
   //common histograms
   fchr6th->FillIntegratedFlow(V6,SdQ[2]);
   fchr6th->FillChi(V6*pow(AvM,0.5));
  }
  else
  {
   cout<<" v_"<<fgkFlow<<"{6} = Im"<<endl;  
  }
  
  //v_2{8}
  if(AvM!=0 && (cumulant[3]<=0.) && (AvQ2x+AvQ2y-pow(AvQx,2.)-pow(AvQy,2.)-pow(pow(-(1./33.)*cumulant[3],(1./8.))*AvM,2.)>0.))
  {  
   ChiQ[3]=AvM*V8/pow(AvQ2x+AvQ2y-pow(AvQx,2.)-pow(AvQy,2.)-pow(V8*AvM,2.),0.5);
   SdQ[3]=(1./(pow(2.*AvM*nEvents,.5)))*pow((12.+96.*pow(ChiQ[3],2)+72.*pow(ChiQ[3],4.)+304.*pow(ChiQ[3],6.)+257.*pow(ChiQ[3],8.)+804.*pow(ChiQ[3],10.)+363.*pow(ChiQ[3],12.)+726.*pow(ChiQ[3],14.))/(726.*pow(ChiQ[3],14.)),.5);
   cout<<" v_"<<fgkFlow<<"{8} = "<<V8<<" +/- "<<SdQ[3]<<", chi{8} = "<<ChiQ[3]<<endl;
   fifr->SetBinContent(4,V8);
   fifr->SetBinError(4,SdQ[3]);
   //common histograms
   fchr8th->FillIntegratedFlow(V8,SdQ[3]);
   fchr8th->FillChi(V8*pow(AvM,0.5));
  } 
  else 
  {
   cout<<" v_"<<fgkFlow<<"{8} = Im"<<endl;     
  }
  
  //v_2{10}
  if (AvM!=0 && cumulant[4]>=0.){
    cout<<"v_"<<fgkFlow<<"{10} = "<<V10<<endl;
    fifr->SetBinContent(5,pow((1./456.)*cumulant[4],(1./10.)));
  } else {
      cout<<"v_"<<fgkFlow<<"{10} = Im"<<endl; 
  }
  
  //v_2{12}
  if (AvM!=0 && AvM!=0 && cumulant[5]<=0.){
    cout<<"v_"<<fgkFlow<<"{12} = "<<V12<<endl;
    fifr->SetBinContent(6,pow(-(1./9460.)*cumulant[5],(1./12.)));
  } else {
    cout<<"v_"<<fgkFlow<<"{12} = Im"<<endl; 
  }
  
  //v_2{14}
  if (AvM!=0 && cumulant[6]>=0.){
    cout<<"v_"<<fgkFlow<<"{14} = "<<V14<<endl;
    fifr->SetBinContent(7,pow((1./274800.)*cumulant[6],(1./14.)));
  } else {
    cout<<"v_"<<fgkFlow<<"{14} = Im"<<endl;  
  }
  
  //v_2{16}
  if (AvM!=0 && cumulant[7]<=0.){
    cout<<"v_"<<fgkFlow<<"{16} = "<<V16<<endl;
    fifr->SetBinContent(8,pow(-(1./10643745.)*cumulant[7],(1./16.)));
  } else {
    cout<<"v_"<<fgkFlow<<"{16} = Im"<<endl;  
  }
  cout<<" "<<endl;
  cout<<"nEvts = "<<nEvents<<", AvM = "<<AvM<<endl; 
  cout<<"*********************************"<<endl;    
  
  
 /////////////////////////////////////////////////////////////////////////////
 ///////////////////////DIFFERENTIAL FLOW CALCULATIONS////////////////////////
 /////////////////////////////////////////////////////////////////////////////
  
  Double_t X[fgknBins][fgkPmax][fgkQmax]={{{0.}}};
  Double_t Y[fgknBins][fgkPmax][fgkQmax]={{{0.}}};
  Double_t BinNoOfParticles[fgknBins]={0.};
      
  //3D profiles
  for(Int_t b=0;b<fgknBins;b++){
    BinNoOfParticles[b]=fBinNoOfParticles->GetBinEntries(b);
    for(Int_t p=0;p<fgkPmax;p++){
      for(Int_t q=0;q<fgkQmax;q++){
	X[b][p][q]=fDiffGenFunRe->GetBinContent(b+1,p+1,q+1)/AvG[p][q];
	Y[b][p][q]=fDiffGenFunIm->GetBinContent(b+1,p+1,q+1)/AvG[p][q];
      }
    }   
  } 
  
  /*
  if(AvM!=0){
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
      //fCommonHistsRes6->FillDifferentialFlow(b+1,v6[b],0.);
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
  
}

  

