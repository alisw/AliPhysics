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
 * for calculation of Q-cumulants *
 * and final flow estimates       *
 *                                *   
 * author:  Ante Bilandzic        * 
 *           (anteb@nikhef.nl)    *
 *********************************/ 

#define AliQCumulantsFunctions_cxx

#include "Riostream.h"
#include "TChain.h"
#include "TFile.h"
#include "TList.h"
#include "TParticle.h"

#include "TProfile.h"
#include "TProfile2D.h" 
#include "TProfile3D.h"
#include "TF1.h"
#include "TAxis.h"
#include "TH1.h"
#include "TH1D.h"

#include "AliFlowEventSimple.h"
#include "AliFlowTrackSimple.h"
#include "AliFlowAnalysisWithCumulants.h"
#include "AliFlowCommonConstants.h"
#include "AliFlowCommonHistResults.h"
#include "AliQCumulantsFunctions.h"

ClassImp(AliQCumulantsFunctions)

//================================================================================================================_

AliQCumulantsFunctions::AliQCumulantsFunctions():  
 fIntRes(NULL),
 fDiffRes2nd(NULL),
 fDiffRes4th(NULL),
 fCovar(NULL),
 fAvMult(NULL),
 fQVector(NULL),
 fQCorr(NULL),
 fQProd(NULL),
 fDirect(NULL),
 fbin2_1n1n(NULL),
 fbin2_2n2n(NULL),
 fbin3_2n1n1n(NULL),
 fbin3_1n1n2n(NULL),
 fbin4_1n1n1n1n(NULL),
 fchr2nd(NULL),
 fchr4th(NULL),
 fchr6th(NULL),
 fchr8th(NULL) 
{
 //default constructor 
}

AliQCumulantsFunctions::~AliQCumulantsFunctions()
{
 //destructor
}

AliQCumulantsFunctions::AliQCumulantsFunctions(TH1D *intRes, TH1D *diffRes2nd, TH1D *diffRes4th, TH1D *covar, TProfile *AvMult, TProfile *QVector, TProfile *QCorr, TProfile *QProd, TProfile *Direct, TProfile *bin2_1n1n, TProfile *bin2_2n2n, TProfile *bin3_2n1n1n, TProfile *bin3_1n1n2n, TProfile *bin4_1n1n1n1n, AliFlowCommonHistResults *chr2nd, AliFlowCommonHistResults *chr4th, AliFlowCommonHistResults *chr6th, AliFlowCommonHistResults *chr8th):
 fIntRes(intRes),
 fDiffRes2nd(diffRes2nd),
 fDiffRes4th(diffRes4th),
 fCovar(covar),
 fAvMult(AvMult),
 fQVector(QVector),
 fQCorr(QCorr),
 fQProd(QProd),
 fDirect(Direct),
 fbin2_1n1n(bin2_1n1n),
 fbin2_2n2n(bin2_2n2n),
 fbin3_2n1n1n(bin3_2n1n1n),
 fbin3_1n1n2n(bin3_1n1n2n),
 fbin4_1n1n1n1n(bin4_1n1n1n1n),
 fchr2nd(chr2nd),
 fchr4th(chr4th),
 fchr6th(chr6th),
 fchr8th(chr8th) 
{
 //custom constructor 
}
      
//================================================================================================================

void AliQCumulantsFunctions::Calculate()
{
 //final results

 //harmonics
 Int_t n = 2; //to be improved 
 
 //--------------------------------------------------------------------------------------------------------- 
 //avarage multiplicity
 Double_t AvM = fAvMult->GetBinContent(1);
 
 //number of events
 Double_t nEvts = fAvMult->GetBinEntries(1);
 
 //2-, 4- and 6-particle azimuthal correlation:
 Double_t two   = fQCorr->GetBinContent(1);  //<<2>>_{n|n}
 Double_t four  = fQCorr->GetBinContent(11); //<<4>>_{n,n|n,n}
 Double_t six   = fQCorr->GetBinContent(21); //<<6>>_{n,n,n|n,n,n}
 //Double_t eight = fQCorr->GetBinContent(); //<<8>>_{n,n,n,n|n,n,n,n}
  
 //errors of 2-, 4- and 6-particle azimuthal correlation:
 Double_t twoErr   = fQCorr->GetBinError(1);  //sigma_{<<2>>_{n|n}} 
 Double_t fourErr  = fQCorr->GetBinError(11); //sigma_{<<4>>_{n,n|n,n}} 
 Double_t sixErr   = fQCorr->GetBinError(21); //sigma_{<<6>>_{n,n,n|n,n,n}}
 //Double_t eightErr = fQCorr->GetBinError(); //sigma_{<<8>>_{n,n,n,n|n,n,n,n}}
 
 //covariances of multi-particle correlations
 Double_t cov24=fQProd->GetBinContent(1)-two*four;     //cov24=<<2>*<4>>-<<2>>*<<4>>
 Double_t cov26=fQProd->GetBinContent(2)-two*six;      //cov26=<<2>*<6>>-<<2>>*<<6>>
 //Double_t cov28=fQProd->GetBinContent(3)-two*eight;  //cov28=<<2>*<8>>-<<2>>*<<8>>
 Double_t cov46=fQProd->GetBinContent(4)-four*six;     //cov46=<<4>*<6>>-<<4>>*<<6>>
 //Double_t cov48=fQProd->GetBinContent(5)-four*eight; //cov48=<<4>*<8>>-<<4>>*<<8>>
 //Double_t cov68=fQProd->GetBinContent(6)-six*eight;  //cov68=<<6>*<8>>-<<6>>*<<8>>
 
 fCovar->SetBinContent(1,cov24);
 fCovar->SetBinContent(2,cov26); 
 //fCovar->SetBinContent(3,cov28); 
 fCovar->SetBinContent(4,cov46);
 //fCovar->SetBinContent(5,cov48); 
 //fCovar->SetBinContent(6,cov68); 
 
 //2nd, 4th and 6th order Q-cumulant and theirs errors:
 Double_t secondOrderQCumulant = two;                //c_n{2} 
 Double_t secondOrderQCumulantError = twoErr;        //sigma_{c_n{2}}
  
 Double_t fourthOrderQCumulant = four-2.*pow(two,2.);//c_n{4}
 Double_t fourthOrderQCumulantError = 0.;            //sigma_{c_n{4}}
 
 if(16.*pow(two,2.)*pow(twoErr,2.)+pow(fourErr,2.)-8.*two*cov24>0.)
 {
  fourthOrderQCumulantError = pow(16.*pow(two,2.)*pow(twoErr,2.)+pow(fourErr,2.)-8.*two*cov24,0.5);//sigma_{c_n{4}}
 }
 
 Double_t sixthOrderQCumulant = six-9.*two*four+12.*pow(two,3.); //c_n{6}
 Double_t sixthOrderQCumulantError = 0.;  //sigma_{c_n{6}}

 if(81.*pow(4.*two*two-four,2.)*pow(twoErr,2.)+81.*pow(two,2.)*pow(fourErr,2.)+pow(sixErr,2.)-162.*(4.*two*two-four)*two*cov24+18.*(4.*two*two-four)*cov26-18.*two*cov46>0.)
 {
  sixthOrderQCumulantError = pow(81.*pow(4.*two*two-four,2.)*pow(twoErr,2.)+81.*pow(two,2.)*pow(fourErr,2.)+pow(sixErr,2.)-162.*(4.*two*two-four)*two*cov24+18.*(4.*two*two-four)*cov26-18.*two*cov46,0.5);//sigma_{c_n{6}}
 }
      
 //integrated flow estimates from Q-cumulants:
 cout<<" "<<endl;
 cout<<"*********************************"<<endl;
 cout<<"*********************************"<<endl;
 cout<<"flow estimates from Q-cumulants:"<<endl;
 cout<<" "<<endl;
 
 Double_t vn2=0.,vn4=0.,vn6=0.;
 Double_t sd2=0.,sd4=0.,sd6=0.; 
 if(secondOrderQCumulant>0.){
  vn2 = pow(secondOrderQCumulant,0.5);//v_n{2}
  sd2 = 0.5*pow(secondOrderQCumulant,-0.5)*secondOrderQCumulantError;//sigma_{v_n{2}}
  cout<<" v_"<<n<<"{2} = "<<vn2<<" +/- "<<sd2<<endl;
  fIntRes->SetBinContent(1,vn2);
  fIntRes->SetBinError(1,sd2);
  //common histograms:
  fchr2nd->FillIntegratedFlow(vn2,sd2);
  fchr2nd->FillChi(vn2*pow(AvM,0.5));
 }else{
  cout<<" v_"<<n<<"{2} = Im"<<endl;
 }
 if(four!=0. && fourthOrderQCumulant<0.){
  vn4 = pow(-fourthOrderQCumulant,1./4.); //v_n{4}
  sd4 = 0.25*pow(-fourthOrderQCumulant,-3./4.)*fourthOrderQCumulantError;//sigma_{v_n{4}}
  cout<<" v_"<<n<<"{4} = "<<vn4<<" +/- "<<sd4<<endl;
  fIntRes->SetBinContent(2,vn4);
  fIntRes->SetBinError(2,sd4);
  //common histograms:
  fchr4th->FillIntegratedFlow(vn4,sd4);
  fchr4th->FillChi(vn4*pow(AvM,0.5));
 }else{
  cout<<" v_"<<n<<"{4} = Im"<<endl;
 }
 if(six!=0. && sixthOrderQCumulant>0.){
  vn6 = pow((1./4.)*sixthOrderQCumulant,1./6.); //v_n{6}
  sd6 = (1./6.)*pow(2.,-1./3.)*pow(sixthOrderQCumulant,-5./6.)*sixthOrderQCumulantError;
  cout<<" v_"<<n<<"{6} = "<<vn6<<" +/- "<<sd6<<endl;
  fIntRes->SetBinContent(3,vn6);
  fIntRes->SetBinError(3,sd6);
  //common histograms:
  fchr6th->FillIntegratedFlow(vn6,sd6);
  fchr6th->FillChi(vn6*pow(AvM,0.5));
 }else{
  cout<<" v_"<<n<<"{6} = Im"<<endl;
 }
 cout<<" "<<endl;
 cout<<"   nEvts = "<<nEvts<<", AvM = "<<AvM<<endl;
 cout<<"*********************************"<<endl;
 cout<<"*********************************"<<endl;
 cout<<" "<<endl; 
//--------------------------------------------------------------------------------------------------------- 
 
//---------------------------------------------------------------------------------------------------------    
//differential flow
Double_t secondOrderQCumulantDiffFlow = 0.;
Double_t fourthOrderQCumulantDiffFlow = 0.;

Int_t nBins = fbin2_1n1n->GetNbinsX();

for(Int_t bb=1;bb<nBins+1;bb++)
{
 if(fbin2_1n1n->GetBinEntries(bb)>0.&&vn2!=0)
 {
  secondOrderQCumulantDiffFlow = fbin2_1n1n->GetBinContent(bb);
  fDiffRes2nd->SetBinContent(bb,secondOrderQCumulantDiffFlow/vn2);
  //common histogram:
  fchr2nd->FillDifferentialFlow(bb,secondOrderQCumulantDiffFlow/vn2,0.);//to be improved (errors)
 }
 if(fbin4_1n1n1n1n->GetBinEntries(bb)>0.&&vn4!=0.)
 {
  fourthOrderQCumulantDiffFlow = fbin4_1n1n1n1n->GetBinContent(bb)-2.*fbin2_1n1n->GetBinContent(bb)*pow(vn2,2.);
  fDiffRes4th->SetBinContent(bb,-1.*fourthOrderQCumulantDiffFlow/pow(vn4,3.));
  //common histogram:
  fchr4th->FillDifferentialFlow(bb,-1.*fourthOrderQCumulantDiffFlow/pow(vn4,3.),0.);//to be improved (errors)
 }
}      
//---------------------------------------------------------------------------------------------------------       
        
         
                                                                                    
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 /*
 Double_t first=16.*pow(two*twoErr,2.);
 Double_t second=pow(fourErr,2.);
 Double_t third=-8.*two*cov24;
 
 cout<<" "<<endl;
 cout<<" "<<endl;
 cout<<" "<<endl;
 cout<<"       1st = "<<(1./16.)*first/(pow(-fourthOrderQCumulant,1.5))<<endl;
 cout<<"       2nd = "<<(1./16.)*second/(pow(-fourthOrderQCumulant,1.5))<<endl;
 cout<<"       3rd = "<<(1./16.)*third/(pow(-fourthOrderQCumulant,1.5))<<endl;
 cout<<" "<<endl;
 cout<<"   nEvts = "<<nEvts<<", AvM = "<<AvM<<endl;
 cout<<" "<<endl;
 cout<<" "<<endl;
 cout<<"c-c = "<<pow((1./16.)*((first+second+third)/(pow(-fourthOrderQCumulant,1.5))),0.5)<<endl; 
 */
 
           
             
 /*
 
 //<Q-vector components>
 Double_t AvQx  = fQVector->GetBinContent(1); //<Q_x>
 Double_t AvQy  = fQVector->GetBinContent(2); //<Q_y>
 Double_t AvQ2x = fQVector->GetBinContent(3); //<(Q_x)^2>
 Double_t AvQ2y = fQVector->GetBinContent(4); //<(Q_y)^2>
 Double_t AvQ2  = fQVector->GetBinContent(5); //<|Q|^2>
 Double_t AvQ4  = fQVector->GetBinContent(6); //<|Q|^4>
 Double_t AvQ_2n2 = fQVector->GetBinContent(7); //<|Q_2n|^2>
 Double_t ReQ2nQnstarQnstar = fQVector->GetBinContent(8); //<Re[Q_n^2 * Q_2n^*]> 
 Double_t ImQ2nQnstarQnstar = fQVector->GetBinContent(9); //<Im[Q_n^2 * Q_2n^*]>
 Double_t AvQ_3n2 = fQVector->GetBinContent(10); //<|Q_3n|^2>
 Double_t AvQ_4n2 = fQVector->GetBinContent(11); //<|Q_4n|^2>
 Double_t AvQ_5n2 = fQVector->GetBinContent(12); //<|Q_5n|^2>
 Double_t ReQ3nQ2nstarQnstar = fQVector->GetBinContent(13); //Re<Q_{3n} Q_{2n}^* Q_{n}^*>
 Double_t ReQ3nQnstarQnstarQnstar = fQVector->GetBinContent(14); //Re<Q_{3n} Q_{n}^* Q_{n}^* Q_{n}^*>
 Double_t HereQ2nQnQ2nstarQnstar = fQVector->GetBinContent(15); //<|Q_{2n}|^2 |Q_{n}|^2>
 Double_t HereQ2nQnQnstarQnstarQnstar = fQVector->GetBinContent(16); //Re<Q_{2n} Q_{n} Q_{n}^* Q_{n}^* Q_{n}^*>
 Double_t AvQ6 = fQVector->GetBinContent(17); //<|Q|^6>
 
 
 
 
*/
 
 /*
 //direct particle correlations
 Double_t twoDirect_n_n = fDirect->GetBinContent(1);
 Double_t fourDirect = fDirect->GetBinContent(2);
 Double_t threeDirect_2n_n_n = fDirect->GetBinContent(3);
 Double_t threeDirect_3n_2n_n = fDirect->GetBinContent(4);
 Double_t fourDirect_3n_n_n_n = fDirect->GetBinContent(5);
 Double_t fourDirect_2n_n_2n_n = fDirect->GetBinContent(6);
 Double_t fiveDirect_2n_n_n_n_n = fDirect->GetBinContent(7);
 Double_t sixDirect_n_n_n_n_n_n = fDirect->GetBinContent(8);
 Double_t twoDirect_2n_2n = fDirect->GetBinContent(9);
 Double_t twoDirect_3n_3n = fDirect->GetBinContent(10);
 Double_t twoDirect_4n_4n = fDirect->GetBinContent(11);
 Double_t twoDirect_5n_5n = fDirect->GetBinContent(12);
 
*/ 
 
 
 
 
 
 
 
/* 
 
 
 //CORRELATIONS
 //two particle correlations
 Double_t two_n_n=(AvQ2-AvM)/(AvM*(AvM-1));
 Double_t two_2n_2n=(AvQ_2n2-AvM)/(AvM*(AvM-1));
 Double_t two_3n_3n=(AvQ_3n2-AvM)/(AvM*(AvM-1));
 Double_t two_4n_4n=(AvQ_4n2-AvM)/(AvM*(AvM-1));
 Double_t two_5n_5n=(AvQ_5n2-AvM)/(AvM*(AvM-1));
 



 


 //four particle correlations:
 //<4>_{n,n|n,n}
 Double_t four_n_n_n_n = (AvQ4+AvQ_2n2-4.*(AvM-2.)*AvQ2-2.*ReQ2nQnstarQnstar)/(AvM*(AvM-1)*(AvM-2)*(AvM-3))+2./((AvM-1)*(AvM-2));
 
 
 //<4>_{3n|n,n,n}
 //Double_t four_3n_n_n_n=(ReQ3nQnstarQnstarQnstar)/(AvM*(AvM-1)*(AvM-2)*(AvM-3))-3.*(three_3n_2n_n+three_2n_n_n)/(AvM-3.)-(3.*two+3.*two_2n+two_3n)/((AvM-2.)*(AvM-3.))-1./((AvM-1.)*(AvM-2.)*(AvM-3.)); OK!!!
 //Double_t four_3n_n_n_n_=(ReQ3nQnstarQnstarQnstar-3.*ReQ3nQ2nstarQnstar-3.*ReQ2nQnstarQnstar)/(AvM*(AvM-1.)*(AvM-2.)*(AvM-3.))+ (2.*AvQ_3n2 + 3.*AvQ_2n2 + 6.*AvQ2-6.*AvM ) /(AvM*(AvM-1.)*(AvM-2.)*(AvM-3.)); OK!!! (final version)


 //<4>_{2n,n|2n,n}
 //Double_t four_2n_n_2n_n = (HereQ2nQnQ2nstarQnstar)/(AvM*(AvM-1)*(AvM-2)*(AvM-3))-2.*(three_2n_n_n+three_3n_2n_n)/(AvM-3.)-((AvM+1.)*two+AvM*two_2n+two_3n)/((AvM-2.)*(AvM-3.))-(AvM)/((AvM-1)*(AvM-2)*(AvM-3)); OK!!!
 //Double_t four_2n_n_2n_n_Alternative = (HereQ2nQnQ2nstarQnstar-2.*ReQ3nQ2nstarQnstar-2.*ReQ2nQnstarQnstar)/(AvM*(AvM-1)*(AvM-2)*(AvM-3))-((AvM-5)*AvQ2+(AvM-4)*AvQ_2n2-AvQ_3n2)/(AvM*(AvM-1)*(AvM-2)*(AvM-3))+(AvM-6)/((AvM-1)*(AvM-2)*(AvM-3)); OK!!! (final version)
 
  
 //<6>_{n,n,n|n,n,n} 
 //Double_t six_n_n_n_n_n_n=AvQ6/(AvM*(AvM-1)*(AvM-2)*(AvM-3)*(AvM-4)*(AvM-5)) - 3.*(2.*five_2n_n_n_n_n)/(AvM-5) - (2.*four_3n_n_n_n+9.*(AvM-2.)*four+9.*four_2n_n_2n_n)/((AvM-4)*(AvM-5)) - 3.*(2.*(3.*AvM-5.)*three_2n_n_n+2.*three_3n_2n_n)/((AvM-3)*(AvM-4)*(AvM-5)) - ((18.*AvM*AvM-45.*AvM+33.)*two+3.*(3.*AvM-4.)*two_2n+two_3n)/((AvM-2)*(AvM-3)*(AvM-4)*(AvM-5)) - (6.*AvM*AvM-9.*AvM+4.)/((AvM-1)*(AvM-2)*(AvM-3)*(AvM-4)*(AvM-5)); OK!!!
 //Double_t six_n_n_n_n_n_n = (AvQ6+9.*HereQ2nQnQ2nstarQnstar-6.*HereQ2nQnQnstarQnstarQnstar)/(AvM*(AvM-1)*(AvM-2)*(AvM-3)*(AvM-4)*(AvM-5)) + 4.*(ReQ3nQnstarQnstarQnstar-3.*ReQ3nQ2nstarQnstar)/(AvM*(AvM-1)*(AvM-2)*(AvM-3)*(AvM-4)*(AvM-5)) + 2.*(9.*(AvM-4.)*ReQ2nQnstarQnstar+2.*AvQ_3n2)/(AvM*(AvM-1)*(AvM-2)*(AvM-3)*(AvM-4)*(AvM-5)) - 9.*(AvQ4+AvQ_2n2)/(AvM*(AvM-1)*(AvM-2)*(AvM-3)*(AvM-5)) + (18.*AvQ2 )/(AvM*(AvM-1)*(AvM-3)*(AvM-4)) - (6.)/((AvM-1)*(AvM-2)*(AvM-3)); OK!!! (final version)
      
       
        


*/

 /*
 
 
 cout<<"<4>_{3n,n,n,n} correlation from Q-vector = "<<four_3n_n_n_n<<endl;AvQ2
 cout<<"<4>_{3n,n,n,n} correlation directly      = "<<fourDirect_3n_n_n_n<<endl;
 cout<<" "<<endl;
 cout<<"<4>_{2n,n,2n,n} correlation from Q-vector = "<<four_2n_n_2n_n<<endl;
 cout<<"<4>_{2n,n,2n,n} correlation directly      = "<<fourDirect_2n_n_2n_n<<endl;
 cout<<" "<<endl;
 cout<<"<5>_{2n,n,n,n,n} correlation from Q-vector = "<<five_2n_n_n_n_n<<endl;
 cout<<"<5>_{2n,n,n,n,n} correlation directly      = "<<fiveDirect_2n_n_n_n_n<<endl;
 cout<<" "<<endl;
 cout<<"<6>_{n,n,n,n,n,n} correlation from Q-vector = "<<six_n_n_n_n_n_n<<endl;
 cout<<"<6>_{n,n,n,n,n,n} correlation directly      = "<<sixDirect_n_n_n_n_n_n<<endl;
 cout<<" "<<endl;
 */
 
 

 
  /*
   
//Q-CUMULANTS
Double_t fourthOrderQCumulant = (AvQ4+AvQ_2n2-2.*ReQ2nQnstarQnstar)/(AvM*(AvM-1.)*(AvM-2.)*(AvM-3.)) - (2.*AvQ2*AvQ2)/(AvM*AvM*(AvM-1.)*(AvM-1.)) - (8.*AvQ2)/(AvM*(AvM-1.)*(AvM-1.)*(AvM-3.)) + 2./((AvM-1.)*(AvM-1.)*(AvM-2.));    
     
Double_t sixthOrderQCumulant = (AvQ6+9.*HereQ2nQnQ2nstarQnstar-6.*HereQ2nQnQnstarQnstarQnstar)/(AvM*(AvM-1)*(AvM-2)*(AvM-3)*(AvM-4)*(AvM-5)) + 4.*(ReQ3nQnstarQnstarQnstar - 3.*ReQ3nQ2nstarQnstar+AvQ_3n2)/(AvM*(AvM-1)*(AvM-2)*(AvM-3)*(AvM-4)*(AvM-5)) + 18.*(4.*AvM+(AvM-5.)*AvQ2)*ReQ2nQnstarQnstar/(AvM*AvM*(AvM-1)*(AvM-1)*(AvM-2)*(AvM-3)*(AvM-5)) - 9.*(4.*AvM+(AvM-5.)*AvQ2)*(AvQ4+AvQ_2n2)/(AvM*AvM*(AvM-1)*(AvM-1)*(AvM-2)*(AvM-3)*(AvM-5)) + 36.*(5.*AvM-11.)*AvQ2/(AvM*(AvM-1)*(AvM-1)*(AvM-1)*(AvM-2)*(AvM-3)*(AvM-4)) + 12.*AvQ2*AvQ2*(6.*AvM-3.*AvQ2+AvM*AvQ2)/(AvM*AvM*AvM*(AvM-1)*(AvM-1)*(AvM-1)*(AvM-3)) - 24./((AvM-1)*(AvM-1)*(AvM-1)*(AvM-2)*(AvM-3));
          
           
                           
                
                 
                   
 
 
 //cout<<" "<<endl;
 //cout<<"should be the same? "<<threeDirect_3n_2n_n<<" "<<three_2n_n_nTemp<<endl;
 //cout<<" "<<endl;
 
 
 
 
 
 
 
 //cout<<" "<<endl;
 cout<<"***************************"<<endl;
 cout<<" "<<endl;
 cout<<"multiplicity = "<<AvM<<endl;
 cout<<" "<<endl;
 
 Double_t two=(AvQ2-AvM)/(AvM*(AvM-1));
 Double_t two_2n=(AvQ_2n2-AvM)/(AvM*(AvM-1));
 Double_t two_3n=(AvQ_3n2-AvM)/(AvM*(AvM-1));
 Double_t two_4n=(AvQ_4n2-AvM)/(AvM*(AvM-1));
 Double_t two_5n=(AvQ_5n2-AvM)/(AvM*(AvM-1));
 

 //cout<<"two's = "<<two<<" "<<two_2n<<" "<<two_3n<<" "<<two_4n<<" "<<two_5n<<endl;
 cout<<" "<<endl;
 Double_t four=(2.*AvM*(AvM-3.)+AvQ4-4.*(AvM-2.)*AvQ2-2.*ReQ2nQnstarQnstar+AvQ_2n2)/(AvM*(AvM-1)*(AvM-2)*(AvM-3));
 Double_t three_2n_n_n=(ReQ2nQnstarQnstar-AvM-2*AvM*(AvM-1)*two-AvM*(AvM-1)*two_2n)/(AvM*(AvM-1)*(AvM-2));
 
 //<3>_{3n,2n,n}  
 Double_t three_3n_2n_n=(ReQ3nQ2nstarQnstar-AvM*(AvM-1)*(two+two_2n+two_3n)-AvM)/(AvM*(AvM-1)*(AvM-2));
 
 //<4>_{3n,n,n,n}  
 //Double_t four_3n_n_n_n=(ReQ3nQnstarQnstarQnstar-3.*AvM*(AvM-1)*(AvM-2)*(three_2n_n_n+three_3n_2n_n)-AvM*(AvM-1)*(3.*two+3.*two_2n+two_3n)-AvM)/(AvM*(AvM-1)*(AvM-2)*(AvM-3));
 
 Double_t four_3n_n_n_n=(ReQ3nQnstarQnstarQnstar)/(AvM*(AvM-1)*(AvM-2)*(AvM-3))-3.*(three_3n_2n_n+three_2n_n_n)/(AvM-3.)-(3.*two+3.*two_2n+two_3n)/((AvM-2.)*(AvM-3.))-1./((AvM-1.)*(AvM-2.)*(AvM-3.));
 

 
 //<4>_{2n,n,2n,n}  
 Double_t four_2n_n_2n_n=(HereQ2nQnQ2nstarQnstar-2.*AvM*(AvM-1)*(AvM-2)*(three_2n_n_n+three_3n_2n_n)-AvM*(AvM-1)*((AvM+1)*two+AvM*two_2n+two_3n)-AvM*AvM)/(AvM*(AvM-1)*(AvM-2)*(AvM-3));
 
 
 
 
 
 
 
 
 //cout<<"Alternative "<<four_2n_n_2n_n<<" "<<four_2n_n_2n_n_Alternative<<" "<<fourDirect_2n_n_2n_n<<endl;
 
 
 //<5>_{2n,n,n,n,n}
 Double_t five_2n_n_n_n_n=(HereQ2nQnQnstarQnstarQnstar-AvM*(AvM-1)*(AvM-2)*(AvM-3)*(four_3n_n_n_n+3.*four+3.*four_2n_n_2n_n)-AvM*(AvM-1)*(AvM-2)*(3.*AvM*three_2n_n_n+4.*three_3n_2n_n)-AvM*(AvM-1)*((9.*AvM-11.)*two+(3.*AvM-2)*two_2n+two_3n)-AvM*(3.*AvM-2))/(AvM*(AvM-1)*(AvM-2)*(AvM-3)*(AvM-4));
 
 
*/
 
 
  /*
  
 //<6>_{n,n,n,n,n,n}
 Double_t six_n_n_n_n_n_n=(AvQ6-6.*AvM*(AvM-1)*(AvM-2+)*(AvM-3)*(AvM-4)*five_2n_n_n_n_n-AvM*(AvM-1)*(AvM-2)*(AvM-3)*(2.*four_3n_n_n_n+9.*(AvM-2.)*four+9.*four_2n_n_2n_n)-6.*AvM*(AvM-1)*(AvM-2)*((3.*AvM-7.)*three_2n_n_n+3.*three_3n_2n_n)-AvM*(AvM-1)*(3.*(2.*AvM*AvM+5.*AvM-13.)*two+3.*(3.*AvM-4.)*two_2n+two_3n)-AvM*(6.*AvM*AvM-9.*AvM+4.))/(AvM*(AvM-1)*(AvM-2)*(AvM-3)*(AvM-4)*(AvM-5));
  
  */
  
  /* 
//<6>_{n,n,n,n,n,n}
 Double_t six_n_n_n_n_n_n=(AvQ6-6.*AvM*(AvM-1)*(AvM-2)*(AvM-3)*(AvM-4)*five_2n_n_n_n_n-AvM*(AvM-1)*(AvM-2)*(AvM-3)*(2.*four_3n_n_n_n+9.*(AvM-2.)*four+9.*four_2n_n_2n_n)-6.*AvM*(AvM-1)*(AvM-2)*((3.*AvM-5.)*three_2n_n_n+three_3n_2n_n)-AvM*(AvM-1)*((18.*AvM*AvM-45.*AvM+33.)*two+3.*(3.*AvM-4.)*two_2n+two_3n)-AvM*(6.*AvM*AvM-9.*AvM+4.))/(AvM*(AvM-1)*(AvM-2)*(AvM-3)*(AvM-4)*(AvM-5));
 
 
 
 
 
 
  
//cout<<"Alternative "<<six_n_n_n_n_n_n<<" "<<six_n_n_n_n_n_n_Alternative<<" "<<sixDirect_n_n_n_n_n_n<<endl;
 
  

 
  
    
      
  
 cout<<"<2> correlation from Q-vector           = "<<two<<endl;
 cout<<"<2> correlation directly                = "<<twoDirect_n_n<<endl;
 cout<<" "<<endl;
 cout<<"<4> correlation from Q-vector           = "<<four<<endl;
 cout<<"<4> correlation directly                = "<<fourDirect<<endl;
 cout<<" "<<endl;
 cout<<"<3>_{2n,n,n} correlation from Q-vector  = "<<three_2n_n_n<<endl;
 cout<<"<3>_{2n,n,n} correlation directly       = "<<threeDirect_2n_n_n<<endl;
 cout<<" "<<endl;
 cout<<"<3>_{3n,2n,n} correlation from Q-vector = "<<three_3n_2n_n<<endl;
 cout<<"<3>_{3n,2n,n} correlation directly      = "<<threeDirect_3n_2n_n<<endl;
 cout<<" "<<endl;
 cout<<"<4>_{3n,n,n,n} correlation from Q-vector = "<<four_3n_n_n_n<<endl;
 cout<<"<4>_{3n,n,n,n} correlation directly      = "<<fourDirect_3n_n_n_n<<endl;
 cout<<" "<<endl;
 cout<<"<4>_{2n,n,2n,n} correlation from Q-vector = "<<four_2n_n_2n_n<<endl;
 cout<<"<4>_{2n,n,2n,n} correlation directly      = "<<fourDirect_2n_n_2n_n<<endl;
 cout<<" "<<endl;
 cout<<"<5>_{2n,n,n,n,n} correlation from Q-vector = "<<five_2n_n_n_n_n<<endl;
 cout<<"<5>_{2n,n,n,n,n} correlation directly      = "<<fiveDirect_2n_n_n_n_n<<endl;
 cout<<" "<<endl;
 cout<<"<6>_{n,n,n,n,n,n} correlation from Q-vector = "<<six_n_n_n_n_n_n<<endl;
 cout<<"<6>_{n,n,n,n,n,n} correlation directly      = "<<sixDirect_n_n_n_n_n_n<<endl;
 cout<<" "<<endl;




 //have so far: <2>, <2>_{2n}, <3>_{2n,n,n}, <4>


 //cout<<"temp re = "<<ReQ2nQnstarQnstar<<endl;
 //cout<<"temp im = "<<ImQ2nQnstarQnstar<<endl;
 //cout<<AvQ2<<" "<<sqrt(AvQ4)<<endl;
 cout<<" "<<endl;

 
 
 
 if(AvQ2>AvM){
  cout<<"v_2{2} = "<<100*sqrt(two)<<"%"<<endl;
 }else{
  cout<<"v_2{2} = Im"<<endl;
 }

 if(four>0){
  cout<<"v_2{4} = "<<100*pow(four,1./4.)<<"%"<<endl;
 }else{
  cout<<"v_2{4} = Im"<<endl;
 }




 cout<<"AvM = "<<AvM<<endl;
 cout<<" "<<endl; 
 cout<<"*************************************"<<endl;
 cout<<"*************************************"<<endl;
 cout<<"flow estimates from Q-cumulants:"<<endl;
 Double_t cumulant2Q=two;
 Double_t cumulant4Q=four-2.*two*two;
 Double_t cumulant6Q=12.*two*two*two-9.*two*four+six_n_n_n_n_n_n;

 //cout<<" "<<endl;
 if(cumulant2Q>0.){
  cout<<"v_2{2} = "<<100*pow(cumulant2Q,1./2.)<<"%"<<endl;
 }else{
  cout<<"v_2{2} = Im"<<endl;
 }
 if(cumulant4Q<0.){
  cout<<"v_2{4} = "<<100*pow(-cumulant4Q,1./4.)<<"%"<<endl;
 }else{
  cout<<"v_2{4} = Im"<<endl;
 }
 if(cumulant6Q>0.){
  cout<<"v_2{6} = "<<100*pow((1./4.)*cumulant6Q,1./6.)<<"%"<<endl;
 }else{
  cout<<"v_2{6} = Im"<<endl;
 }
 
 cout<<"*************************************"<<endl;
 cout<<"*************************************"<<endl;
 cout<<" "<<endl;
 
 
 cout<<" "<<endl; 
 cout<<"*************************************"<<endl;
 cout<<"*************************************"<<endl;
 cout<<"flow estimates from Q-cumulants"<<endl; 
 cout<<"without multiplicities fluctuations:"<<endl;
 Double_t cumulant2QnoFluct=fQCorr->GetBinContent(1);
 Double_t cumulant4QnoFluct=fQCorr->GetBinContent(11)-2.*fQCorr->GetBinContent(1)*fQCorr->GetBinContent(1);
 Double_t cumulant6QnoFluct=12.*fQCorr->GetBinContent(1)*fQCorr->GetBinContent(1)*fQCorr->GetBinContent(1) - 9.*fQCorr->GetBinContent(1)*fQCorr->GetBinContent(11)+fQCorr->GetBinContent(21);

 //cout<<" "<<endl;
 if(cumulant2QnoFluct>0.){
  cout<<"v_2{2} = "<<100*pow(cumulant2QnoFluct,1./2.)<<"%"<<endl;
  fIntRes->SetBinContent(1,100*pow(cumulant2QnoFluct,1./2.));
 }else{
  cout<<"v_2{2} = Im"<<endl;
 }
 if(cumulant4QnoFluct<0.){
  cout<<"v_2{4} = "<<100*pow(-cumulant4QnoFluct,1./4.)<<"%"<<endl;
  fIntRes->SetBinContent(2,100*pow(-cumulant4QnoFluct,1./4.));
 }else{
  cout<<"v_2{4} = Im"<<endl;
 }
 if(cumulant6QnoFluct>0.){
  cout<<"v_2{6} = "<<100*pow((1./4.)*cumulant6QnoFluct,1./6.)<<"%"<<endl;
  fIntRes->SetBinContent(3,100*pow((1./4.)*cumulant6QnoFluct,1./6.));
 }else{
  cout<<"v_2{6} = Im"<<endl;
 }
 
 cout<<"*************************************"<<endl;
 cout<<"*************************************"<<endl;
 cout<<" "<<endl;
 

 
  
   
    
 cout<<"*****************************************************************"<<endl;
 cout<<"*****************************************************************"<<endl;
 cout<<"direct (nested loops) correlations vs correlations from Q-vectors"<<endl;
 cout<<"for ARBITRARY multiplicity:"<<endl;
 cout<<" "<<endl;
 cout<<"avarage multiplicity = "<<AvM<<endl;
 cout<<" "<<endl;
 cout<<"<2>_{n|n} correlation from Q-vector         = "<<fQCorr->GetBinContent(1)<<endl;
 cout<<"<2>_{n|n} correlation nested loops          = "<<fDirect->GetBinContent(1)<<endl;
 cout<<" "<<endl;
 cout<<"<4>_{n,n|n,n} correlation from Q-vector     = "<<fQCorr->GetBinContent(11)<<endl;
 cout<<"<4>_{n,n|n,n} correlation nested loops      = "<<fDirect->GetBinContent(2)<<endl; 
 cout<<" "<<endl;
 cout<<"<6>_{n,n,n|n,n,n} correlation from Q-vector = "<<fQCorr->GetBinContent(21)<<endl;
 cout<<"<6>_{n,n,n|n,n,n} correlation nested loops  = "<<fDirect->GetBinContent(8)<<endl; 
 cout<<"*****************************************************************"<<endl;
 cout<<"*****************************************************************"<<endl;
     
      
       
        
         
           
 
 cout<<"****** BDO ******"<<endl;
 if(BDO4<0.){
  cout<<"v_2{4} = "<<100*pow(-BDO4,1./4.)<<"%"<<endl;
 }else{
  cout<<"v_2{4} = Im"<<endl;
 }
 if(BDO6>0.){
  cout<<"v_2{6} = "<<100*pow((1./4.)*BDO6,1./6.)<<"%"<<endl;
 }else{
  cout<<"v_2{6} = Im"<<endl;
 }

 
 
 cout<<" "<<endl;
 cout<<"fourth order Q-cumulant should be the same and BDO4? "<<fourthOrderQCumulant<<" "<<cumulant4Q<<" "<<BDO4<<endl;
 cout<<" "<<endl;
 cout<<" "<<endl;
 cout<<"sixth order Q-cumulant should be the same and BDO6? "<<sixthOrderQCumulant<<" "<<cumulant6Q<<" "<<BDO6<<endl;
 cout<<" "<<endl;
 
 cout<<"*************************************"<<endl;
 cout<<"*************************************"<<endl;
 cout<<" "<<endl;
 
 
 
  cout<<"*********************"<<endl;
 cout<<"multiplicity = "<<AvM<<endl;
 cout<<"correlations: "<<endl;
 cout<<"direct = "<<twoDirect_n_n<<endl;
 cout<<"approx = "<<two_n_n<<endl;
 cout<<"weight = "<<fQCorr->GetBinContent(1)<<endl;
 cout<<"*********************"<<endl;
 cout<<"test = "<<fQCorr->GetBinContent(11)<<" "<<four<<endl;
 
 
 
 
 
 
 
 
 
 
 

*/
















/*






//needed for direct correlations


 cout<<" "<<endl;
 cout<<" "<<endl;
 cout<<"   **** cross-checking the formulas ****"<<endl;
 cout<<"   ****     for integrated flow     ****"<<endl;
 cout<<"(selected only events for which 8 < M < 24):"<<endl;
 cout<<" "<<endl;
 cout<<"   nEvts = "<<nEvts<<", AvM = "<<AvM<<endl;
 cout<<" "<<endl;
 cout<<"<2>_{n|n} from Q-vectors              = "<<fQCorr->GetBinContent(1)<<endl;
 cout<<"<2>_{n|n} from nested loops           = "<<fDirect->GetBinContent(1)<<endl;
 cout<<" "<<endl;
 cout<<"<2>_{2n|2n} from Q-vectors            = "<<fQCorr->GetBinContent(2)<<endl;
 cout<<"<2>_{2n|2n} from nested loops         = "<<fDirect->GetBinContent(2)<<endl;
 cout<<" "<<endl;
 cout<<"<2>_{3n|3n} from Q-vectors            = "<<fQCorr->GetBinContent(3)<<endl;
 cout<<"<2>_{3n|3n} from nested loops         = "<<fDirect->GetBinContent(3)<<endl;
 cout<<" "<<endl;
 cout<<"<2>_{4n|4n} from Q-vectors            = "<<fQCorr->GetBinContent(4)<<endl;
 cout<<"<2>_{4n|4n} from nested loops         = "<<fDirect->GetBinContent(4)<<endl;
 cout<<" "<<endl;
 cout<<"<3>_{2n,n,n} from Q-vectors           = "<<fQCorr->GetBinContent(6)<<endl;
 cout<<"<3>_{2n,n,n} from nested loops        = "<<fDirect->GetBinContent(6)<<endl;
 cout<<" "<<endl;
 cout<<"<3>_{3n,2n,n} from Q-vectors          = "<<fQCorr->GetBinContent(7)<<endl;
 cout<<"<3>_{3n,2n,n} from nested loops       = "<<fDirect->GetBinContent(7)<<endl;
 cout<<" "<<endl; 
 cout<<"<3>_{4n,2n,2n} from Q-vectors         = "<<fQCorr->GetBinContent(8)<<endl;
 cout<<"<3>_{4n,2n,2n} from nested loops      = "<<fDirect->GetBinContent(8)<<endl;
 cout<<" "<<endl;
 cout<<"<3>_{4n,3n,n} from Q-vectors          = "<<fQCorr->GetBinContent(9)<<endl;
 cout<<"<3>_{4n,3n,n} from nested loops       = "<<fDirect->GetBinContent(9)<<endl;
 cout<<" "<<endl;
 cout<<"<4>_{n,n|n,n} from Q-vectors          = "<<fQCorr->GetBinContent(11)<<endl;
 cout<<"<4>_{n,n|n,n} from nested loops       = "<<fDirect->GetBinContent(11)<<endl;
 cout<<" "<<endl;
 cout<<"<4>_{2n,2n|2n,2n} from Q-vectors      = "<<fQCorr->GetBinContent(31)<<endl;//to be improved
 cout<<"<4>_{2n,2n|2n,2n} from nested loops   = "<<fDirect->GetBinContent(31)<<endl;//to be improved 
 cout<<" "<<endl;
 cout<<"<4>_{3n,n|3n,n} from Q-vectors        = "<<fQCorr->GetBinContent(32)<<endl;//to be improved
 cout<<"<4>_{3n,n|3n,n} from nested loops     = "<<fDirect->GetBinContent(32)<<endl;//to be improved 
 cout<<" "<<endl;
 cout<<"<4>_{2n,n|2n,n} from Q-vectors        = "<<fQCorr->GetBinContent(12)<<endl;
 cout<<"<4>_{2n,n|2n,n} from nested loops     = "<<fDirect->GetBinContent(12)<<endl;
 cout<<" "<<endl;
 cout<<"<4>_{3n|n,n,n} from Q-vectors         = "<<fQCorr->GetBinContent(13)<<endl;
 cout<<"<4>_{3n|n,n,n} from nested loops      = "<<fDirect->GetBinContent(13)<<endl;
 cout<<" "<<endl;
 cout<<"<4>_{4n|2n,n,n} from Q-vectors        = "<<fQCorr->GetBinContent(14)<<endl;
 cout<<"<4>_{4n|2n,n,n} from nested loops     = "<<fDirect->GetBinContent(14)<<endl;
 cout<<" "<<endl;
 cout<<"<4>_{3n,n|2n,2n} from Q-vectors       = "<<fQCorr->GetBinContent(15)<<endl;
 cout<<"<4>_{3n,n|2n,2n} from nested loops    = "<<fDirect->GetBinContent(15)<<endl;
 cout<<" "<<endl;
 cout<<"<5>_{2n,n|n,n,n} from Q-vectors       = "<<fQCorr->GetBinContent(16)<<endl;
 cout<<"<5>_{2n,n|n,n,n} from nested loops    = "<<fDirect->GetBinContent(16)<<endl;
 cout<<" "<<endl;
 cout<<"<5>_{2n,2n|2n,n,n} from Q-vectors     = "<<fQCorr->GetBinContent(17)<<endl;
 cout<<"<5>_{2n,2n|2n,n,n} from nested loops  = "<<fDirect->GetBinContent(17)<<endl;
 cout<<" "<<endl;
 cout<<"<5>_{3n,n|2n,n,n} from Q-vectors      = "<<fQCorr->GetBinContent(18)<<endl;
 cout<<"<5>_{3n,n|2n,n,n} from nested loops   = "<<fDirect->GetBinContent(18)<<endl;
 cout<<" "<<endl;
 cout<<"<5>_{4n|n,n,n,n} from Q-vectors       = "<<fQCorr->GetBinContent(19)<<endl;
 cout<<"<5>_{4n|n,n,n,n} from nested loops    = "<<fDirect->GetBinContent(19)<<endl;
 cout<<" "<<endl;
 */

 /*
 
 cout<<"<6>_{n,n,n|n,n,n} from Q-vectors    = "<<fQCorr->GetBinContent(21)<<endl;
 cout<<"<6>_{n,n,n|n,n,n} from nested loops = "<<fDirect->GetBinContent(21)<<endl;
 cout<<" "<<endl; 
 */




//DIFFERENTIAL FLOW


 //41st bin: <2'>_{n|n}
 //42nd bin: <2'>_{2n|2n}
 //46th bin: <3'>_{2n|n,n}
 //47th bin: <3'>_{n,n|2n}
 //51st bin: <4'>_{n,n|n,n}



/*
 
 cout<<" "<<endl;
 cout<<" "<<endl;
 cout<<"   **** cross-checking the formulas ****"<<endl;
 cout<<"   ****    for differential flow    ****"<<endl;
 cout<<"(selected only events for which 8 < M < 24):"<<endl;
 cout<<" "<<endl; 
 cout<<"nEvts = "<<nEvts<<", AvM = "<<AvM<<endl;
 cout<<"0.5 < Pt < 0.6 GeV"<<endl;                                
 cout<<" "<<endl;                                       
 cout<<"<2'>_{n|n} from Q-vectors        = "<<fbin2_1n1n->GetBinContent(6)<<endl;
 cout<<"<2'>_{n|n} from nested loops     = "<<fDirect->GetBinContent(41)<<endl;
 cout<<" "<<endl;                                       
 cout<<"<2'>_{2n|2n} from Q-vectors      = "<<fbin2_2n2n->GetBinContent(6)<<endl;
 cout<<"<2'>_{2n|2n} from nested loops   = "<<fDirect->GetBinContent(42)<<endl;                                        
 cout<<" "<<endl;  
 cout<<"<3'>_{2n|n,n} from Q-vectors     = "<<fbin3_2n1n1n->GetBinContent(6)<<endl;
 cout<<"<3'>_{2n|n,n} from nested loops  = "<<fDirect->GetBinContent(46)<<endl;                   
 cout<<" "<<endl;              
 cout<<"<3'>_{n,n|2n} from Q-vectors     = "<<fbin3_1n1n2n->GetBinContent(6)<<endl;
 cout<<"<3'>_{n,n|2n} from nested loops  = "<<fDirect->GetBinContent(47)<<endl;                                 
 cout<<" "<<endl;              
 cout<<"<4'>_{n,n|n,n} from Q-vectors    = "<<fbin4_1n1n1n1n->GetBinContent(6)<<endl;
 cout<<"<4'>_{n,n|n,n} from nested loops = "<<fDirect->GetBinContent(51)<<endl;                                                                                   
 cout<<" "<<endl;   
 
 */     
           
        
 /*
                     
                                                             
cout<<" "<<endl;      
cout<<" "<<endl;
cout<<"   nEvts = "<<nEvts<<", AvM = "<<AvM<<endl;
cout<<" "<<endl;                         
cout<<"sigma_2 = "<<twoErr<<endl;             
cout<<"sigma_4 = "<<fourErr<<endl;
cout<<"sigma_6 = "<<sixErr<<endl;


*/
                                              
                                                

}

  













