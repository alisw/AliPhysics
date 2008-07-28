//******************************* 
// flow analysis with cumulants *   
// author: Ante Bilandzic       * 
// email: anteb@nikhef.nl       *
//******************************* 

#define AliCumulantsFunctions_cxx

#include "Riostream.h"
#include "AliFlowCommonConstants.h"
#include "TChain.h"
#include "TFile.h"
#include "TList.h" //NEW
#include "TParticle.h"

#include "TProfile.h"
#include "TProfile2D.h" 
#include "TProfile3D.h"

#include "AliFlowEventSimple.h"
#include "AliFlowTrackSimple.h"
#include "AliFlowAnalysisWithCumulants.h"
#include "AliFlowCumuConstants.h"
#include "AliFlowCommonConstants.h"
#include "AliCumulantsFunctions.h"


ClassImp(AliCumulantsFunctions)

//________________________________________________________________________

AliCumulantsFunctions::AliCumulantsFunctions():  
  fIG(NULL),
  fDGRe(NULL),
  fDGIm(NULL),
  fifr(NULL),
  fdfr2(NULL), 
  fdfr4(NULL), 
  fdfr6(NULL),
  fdfr8(NULL),
  fAvMult(0)
 {
   //default constructor 
 }

AliCumulantsFunctions::~AliCumulantsFunctions(){
  //desctructor
}

AliCumulantsFunctions::AliCumulantsFunctions(TProfile2D *IntGen, TProfile3D *DiffGenRe, TProfile3D *DiffGenIm, TH1D *ifr, TH1D *dfr2, TH1D *dfr4, TH1D *dfr6, TH1D *dfr8, Double_t CvM):
  fIG(IntGen),
  fDGRe(DiffGenRe),
  fDGIm(DiffGenIm),
  fifr(ifr),
  fdfr2(dfr2), 
  fdfr4(dfr4), 
  fdfr6(dfr6),
  fdfr8(dfr8),
  fAvMult(CvM)
 {
   //custom constructor 
 }
  
//___________________________________________________________________________
void AliCumulantsFunctions::Calculate(){
 //calculate final flow estimates
 
     static const Int_t fgkQmax=AliFlowCumuConstants::kQmax;//needed for numerics
     static const Int_t fgkPmax=AliFlowCumuConstants::kPmax;//needed for numerics  
     static const Int_t fgkFlow=AliFlowCumuConstants::kFlow;//integrated flow coefficient to be calculated
     static const Int_t fgkMltpl=AliFlowCumuConstants::kMltpl;//the multiple in p=m*n (diff. flow) 
     static const Int_t fgknBins=100;//number of pt bins //to be improved
     Double_t fR0=AliFlowCumuConstants::fgR0;//needed for numerics
     Double_t fPtMax=AliFlowCommonConstants::GetPtMax();//maximum pt
     Double_t fPtMin=AliFlowCommonConstants::GetPtMin();//minimum pt
     Double_t fBinWidth=(fPtMax-fPtMin)/fgknBins;//width of pt bin (in GeV)   
     
  Double_t fBvG[fgkPmax][fgkQmax]={0.};
        
  for(Int_t p=0;p<fgkPmax;p++){
   for(Int_t q=0;q<fgkQmax;q++){ 
    fBvG[p][q]=fIG->GetBinContent(p+1,q+1);
   }
  } 
      
 
  /////////////////////////////////////////////////////////////////////////////      
  //////////////////gen. function for the cumulants////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  
  Double_t fC[fgkPmax][fgkQmax]={0.};
  for (Int_t p=0;p<fgkPmax;p++){
    for (Int_t q=0;q<fgkQmax;q++){
      fC[p][q]=1.*fAvMult*(pow(fBvG[p][q],(1./fAvMult))-1.); //to be improved
    }
  }
 
  /////////////////////////////////////////////////////////////////////////////
  ///////avaraging the gen. function for the cumulants over azimuth////////////
  /////////////////////////////////////////////////////////////////////////////
  
  Double_t fCAv[fgkPmax]={0.}; 
  for (Int_t p=0;p<fgkPmax;p++){
    Double_t fTempHere=0.; 
    for (Int_t q=0;q<fgkQmax;q++){
      fTempHere+=1.*fC[p][q];
    } 
    fCAv[p]=1.*fTempHere/fgkQmax;
  }
  
  /////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////final results//////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  
  Double_t fCumulant[fgkPmax];//c_{iFlow}\{2(p+1)\}
  
  fCumulant[0]=(1./(fR0*fR0)) * (8.*fCAv[0] - 14.*fCAv[1] + (56./3.)*fCAv[2] - (35./2.)*fCAv[3] + 
			      (56./5.)*fCAv[4] - (14./3.)*fCAv[5] + (8./7.)*fCAv[6] - (1./8.)*fCAv[7]);

  fCumulant[1]=(1./pow(fR0,4.)) * ((-1924./35.)*fCAv[0] + (621./5.)*fCAv[1] - (8012./45.)*fCAv[2] + 
				 (691./4.)*fCAv[3] - (564./5.)*fCAv[4] + (2143./45.)*fCAv[5] - 
				 (412./35.)*fCAv[6] + (363./280.)*fCAv[7]);

  fCumulant[2]=(1./pow(fR0,6.)) * (349.*fCAv[0] - (18353./20.)*fCAv[1] + (7173./5.)*fCAv[2] - 
				 1457.*fCAv[3] + (4891./5.)*fCAv[4] - (1683./4.)*fCAv[5] + 
				 (527./5.)*fCAv[6] - (469./40.)*fCAv[7]);

  fCumulant[3]=(1./pow(fR0,8.)) * ((-10528./5.)*fCAv[0] + (30578./5.)*fCAv[1] - (51456./5.)*fCAv[2] + 
				 10993.*fCAv[3] - (38176./5.)*fCAv[4] + (16818./5.)*fCAv[5] - 
				 (4288./5.)*fCAv[6] + (967./10.)*fCAv[7]);

  fCumulant[4]=(1./pow(fR0,10.)) * (11500.*fCAv[0] - 35800.*fCAv[1] + 63900.*fCAv[2] - 71600.*fCAv[3] + 
				  51620.*fCAv[4] - 23400.*fCAv[5] + 6100.*fCAv[6] - 700.*fCAv[7]);

  fCumulant[5]=(1./pow(fR0,12.)) * (-52560.*fCAv[0] + 172080.*fCAv[1] - 321840.*fCAv[2] + 376200.*fCAv[3] - 
				  281520.*fCAv[4] + 131760.*fCAv[5] - 35280.*fCAv[6] + 4140.*fCAv[7]);

  fCumulant[6]=(1./pow(fR0,14.)) * (176400.*fCAv[0] - 599760.*fCAv[1] + 1164240.*fCAv[2] - 1411200.*fCAv[3] + 
				  1093680.*fCAv[4] - 529200.*fCAv[5] + 146160.*fCAv[6] - 17640.*fCAv[7]);

  fCumulant[7]=(1./pow(fR0,16.)) * (-322560*fCAv[0] + 1128960.*fCAv[1] - 2257920.*fCAv[2] + 2822400.*fCAv[3] - 
				  2257920.*fCAv[4] + 1128960.*fCAv[5] - 322560.*fCAv[6] + 40320.*fCAv[7]);
  
  
  cout<<""<<endl;
  cout<<"***************************"<<endl;
  cout<<"cumulants:"<<endl;
  
  cout<<" c_"<<fgkFlow<<"{2} = "<<fCumulant[0]<<endl; 
  cout<<" c_"<<fgkFlow<<"{4} = "<<fCumulant[1]<<endl;
  cout<<" c_"<<fgkFlow<<"{6} = "<<fCumulant[2]<<endl;
  cout<<" c_"<<fgkFlow<<"{8} = "<<fCumulant[3]<<endl; 
  cout<<"c_"<<fgkFlow<<"{10} = "<<fCumulant[4]<<endl; 
  cout<<"c_"<<fgkFlow<<"{12} = "<<fCumulant[5]<<endl;
  cout<<"c_"<<fgkFlow<<"{14} = "<<fCumulant[6]<<endl; 
  cout<<"c_"<<fgkFlow<<"{16} = "<<fCumulant[7]<<endl; 
  
  cout<<""<<endl;
  cout<<"integrated flow: "<<endl;
  
  
  Double_t fV2=0.,fV4=0.,fV6=0.,fV8=0.;
  Double_t fSdQ[4]={0.};
  Double_t fChiQ[4]={0.};
   if (fCumulant[0]>=0.){ 
    fV2=sqrt(fCumulant[0]);    
    //if (fAvQ2x+fAvQ2y-pow(fAvQx,2.)-pow(fAvQy,2.)-pow(fV2*fAvM,2.)>0.){
     //fChiQ[0]=fAvM*fV2/pow(fAvQ2x+fAvQ2y-pow(fAvQx,2.)-pow(fAvQy,2.)-pow(fV2*fAvM,2.),0.5);
     //fSdQ[0]=pow(((1./(2.*fAvM*nEvents))*((1.+1.*pow(fChiQ[0],2))/(1.*pow(fChiQ[0],2)))),0.5);
     cout<<" v_"<<fgkFlow<<"{2} = "<<100.*fV2<<"%, chi{2} = "<<fChiQ[0]<<", sd{2} = "<<100.*fSdQ[0]<<"%"<<endl;//to be improved (2->fgkFlow)
     //fCommonHistsRes2->FillIntegratedFlow(100.*fV2,100.*fSdQ[0]);
     //fCommonHistsRes2->FillChi(fChiQ[0]);
     fifr->SetBinContent(1,100.*fV2);
    //}
   //} else {
    //cout<<" v_"<<fgkFlow<<"{2} = Im"<<endl;  
   } else {
    //cout<<" v_"<<fgkFlow<<"{8} = Im"<<endl; 
   }
  if (fCumulant[1]<=0.){
    fV4=pow(-fCumulant[1],(1./4.));
    //if (fAvQ2x+fAvQ2y-pow(fAvQx,2.)-pow(fAvQy,2.)-pow(fV4*fAvM,2.)>0.){
     //fChiQ[1]=fAvM*fV4/pow(fAvQ2x+fAvQ2y-pow(fAvQx,2.)-pow(fAvQy,2.)-pow(fV4*fAvM,2.),0.5);
     //fSdQ[1]=(1./(pow(2.*fAvM*nEvents,.5)))*pow((1.+2.*pow(fChiQ[1],2)+(1./4.)*pow(fChiQ[1],4.)+(1./4.)*pow(fChiQ[1],6.))/((1./4.)*pow(fChiQ[1],6.)),.5);
     cout<<" v_"<<fgkFlow<<"{4} = "<<100.*fV4<<"%, chi{4} = "<<fChiQ[1]<<", sd{4} = "<<100.*fSdQ[1]<<"%"<<endl;//to be improved (2->fgkFlow)
     //fCommonHistsRes4->FillInteg8ratedFlow(100.*fV4,100.*fSdQ[1]);
     //fCommonHistsRes4->FillChi(fChiQ[1]);
     fifr->SetBinContent(2,100.*fV4);
    //} else {
      //cout<<" v_"<<fgkFlow<<"{4} = Im"<<endl;
    //}
 // } else {
   // cout<<" v_"<<fgkFlow<<"{4} = Im"<<endl;  
  } 
  if (fCumulant[2]>=0.){
    fV6=pow((1./4.)*fCumulant[2],(1./6.));
    //if (fAvQ2x+fAvQ2y-pow(fAvQx,2.)-pow(fAvQy,2.)-pow(fV6*fAvM,2.)>0.){
     //fChiQ[2]=fAvM*fV6/pow(fAvQ2x+fAvQ2y-pow(fAvQx,2.)-pow(fAvQy,2.)-pow(fV6*fAvM,2.),0.5);
     //fSdQ[2]=(1./(pow(2.*fAvM*nEvents,.5)))*pow((3.+18.*pow(fChiQ[2],2)+9.*pow(fChiQ[2],4.)+28.*pow(fChiQ[2],6.)+12.*pow(fChiQ[2],8.)+24.*pow(fChiQ[2],10.))/(24.*pow(fChiQ[2],10.)),.5);
     cout<<" v_"<<fgkFlow<<"{6} = "<<100.*fV6<<"%, chi{6} = "<<fChiQ[2]<<", sd{6} = "<<100.*fSdQ[2]<<"%"<<endl;//to be improved (2->fgkFlow)
     //fCommonHistsRes6->FillIntegratedFlow(100.*fV6,100.*fSdQ[2]);
     //fCommonHistsRes6->FillChi(fChiQ[2]);
     fifr->SetBinContent(3,100.*fV6);
    //} else {
     // cout<<" v_"<<fgkFlow<<"{6} = Im"<<endl; 
   // }
  } else {
    //cout<<" v_"<<fgkFlow<<"{6} = Im"<<endl;  
  }
  if (fCumulant[3]<=0.){
    fV8=pow(-(1./33.)*fCumulant[3],(1./8.));
    //if (fAvQ2x+fAvQ2y-pow(fAvQx,2.)-pow(fAvQy,2.)-pow(fV8*fAvM,2.)>0.){
     //fChiQ[3]=fAvM*fV8/pow(fAvQ2x+fAvQ2y-pow(fAvQx,2.)-pow(fAvQy,2.)-pow(fV8*fAvM,2.),0.5);
     //fSdQ[3]=(1./(pow(2.*fAvM*nEvents,.5)))*pow((12.+96.*pow(fChiQ[3],2)+72.*pow(fChiQ[3],4.)+304.*pow(fChiQ[3],6.)+257.*pow(fChiQ[3],8.)+804.*pow(fChiQ[3],10.)+363.*pow(fChiQ[3],12.)+726.*pow(fChiQ[3],14.))/(726.*pow(fChiQ[3],14.)),.5);
      cout<<" v_"<<fgkFlow<<"{8} = "<<100.*fV8<<"%, chi{8} = "<<fChiQ[3]<<", sd{8} = "<<100.*fSdQ[3]<<"%"<<endl;//to be improved (2->fgkFlow)
      //fCommonHistsRes8->FillIntegratedFlow(100.*fV8,100.*fSdQ[3]);
      //fCommonHistsRes8->FillChi(fChiQ[3]);
      fifr->SetBinContent(4,100.*fV8);
     //} else {
       //cout<<" v_"<<fgkFlow<<"{8} = Im"<<endl;a
     //}
  } else {
    //cout<<" v_"<<fgkFlow<<"{8} = Im"<<endl;     
  }
  if (fCumulant[4]>=0.){//fHistProIntFlow
    cout<<"v_"<<fgkFlow<<"{10} = "<<100.*pow((1./456.)*fCumulant[4],(1./10.))<<"%"<<endl;//to be improved (2->fgkFlow)
    fifr->SetBinContent(5,100.*pow((1./456.)*fCumulant[4],(1./10.)));
  } else {
      cout<<"v_"<<fgkFlow<<"{10} = Im"<<endl;  //to be improved (2->fgkFlow)
  }
  if (fCumulant[5]<=0.){
    cout<<"v_"<<fgkFlow<<"{12} = "<<100.*pow(-(1./9460.)*fCumulant[5],(1./12.))<<"%"<<endl;//to be improved (2->fgkFlow)
    fifr->SetBinContent(6,100.*pow(-(1./9460.)*fCumulant[5],(1./12.)));
  } else {
    cout<<"v_"<<fgkFlow<<"{12} = Im"<<endl;  //to be improved (2->fgkFlow)
  }
  if (fCumulant[6]>=0.){
    cout<<"v_"<<fgkFlow<<"{14} = "<<100.*pow((1./274800.)*fCumulant[6],(1./14.))<<"%"<<endl;//to be improved (2->fgkFlow)
    fifr->SetBinContent(7,100.*pow((1./274800.)*fCumulant[6],(1./14.)));
  } else {
    cout<<"v_"<<fgkFlow<<"{14} = Im"<<endl;  //to be improved (2->fgkFlow)
  }
  if (fCumulant[7]<=0.){
    cout<<"v_"<<fgkFlow<<"{16} = "<<100.*pow(-(1./10643745.)*fCumulant[7],(1./16.))<<"%"<<endl;//to be improved (2->fgkFlow)
    fifr->SetBinContent(8,100.*pow(-(1./10643745.)*fCumulant[7],(1./16.)));
  } else {
    cout<<"v_"<<fgkFlow<<"{16} = Im"<<endl;  //to be improved (2->fgkFlow)
  }
  //cout<<"***************************"<<endl;
      
  //DIFFERENTIAL FLOW CALCULATIONS STARTS HERE!!!
  
  Double_t fX[fgknBins][fgkPmax][fgkQmax]={0.};//see the text bellow relation (11) in PG
  Double_t fY[fgknBins][fgkPmax][fgkQmax]={0.};
  
  for(Int_t b=0;b<fgknBins;b++){
    for(Int_t p=0;p<fgkPmax;p++){
      for(Int_t q=0;q<fgkQmax;q++){
	fX[b][p][q]=fDGRe->GetBinContent(b+1,p+1,q+1)/fBvG[p][q];
	fY[b][p][q]=fDGIm->GetBinContent(b+1,p+1,q+1)/fBvG[p][q];
      }
    }   
  } 
  
 
  Double_t fD[fgknBins][fgkPmax]={0.};//implementing relation (11) from PG
  
  for (Int_t b=0;b<fgknBins;b++){
    for (Int_t p=0;p<fgkPmax;p++){
      Double_t fTempHere3=0.; 
      for (Int_t q=0;q<fgkQmax;q++){
	fTempHere3+=cos(fgkMltpl*2.*q*TMath::Pi()/fgkQmax)*fX[b][p][q] + sin(fgkMltpl*2.*q*TMath::Pi()/fgkQmax)*fY[b][p][q];
      } 
      fD[b][p]=1.*(pow(fR0*pow(p+1.,.5),fgkMltpl)/fgkQmax)*fTempHere3;
      if(fD[b][p]){
       //cout<<"And this "<<b<<" "<<fD[b][p]<<endl;
      }  
    }
  } 
  
  Double_t fDiffCumulant2[fgknBins]={0.};//implementing relation (12) from PG
  Double_t fDiffCumulant4[fgknBins]={0.};
  Double_t fDiffCumulant6[fgknBins]={0.};
  Double_t fDiffCumulant8[fgknBins]={0.};
  
  for (Int_t b=0;b<fgknBins;b++){
    fDiffCumulant2[b]=(1./(fR0*fR0))*(4.*fD[b][0]-3.*fD[b][1]+(4./3.)*fD[b][2]-(1./4.)*fD[b][3]);
    fDiffCumulant4[b]=(1./pow(fR0,4.))*((-26./3.)*fD[b][0]+(19./2.)*fD[b][1]-(14./3.)*fD[b][2]+(11./12.)*fD[b][3]);
    fDiffCumulant6[b]=(1./pow(fR0,6.))*(18.*fD[b][0]-24.*fD[b][1]+14.*fD[b][2]-3.*fD[b][3]);
    fDiffCumulant8[b]=(1./pow(fR0,8.))*(-24.*fD[b][0]+36.*fD[b][1]-24.*fD[b][2]+6.*fD[b][3]);
  }
  
  Double_t fv2[fgknBins],fv4[fgknBins],fv6[fgknBins],fv8[fgknBins];
  //Double_t fAvPt[fgknBins];
  Double_t fSddiff2[fgknBins],fSddiff4[fgknBins];

  cout<<"number of pt bins: "<<fgknBins<<endl;
  cout<<"****************************************"<<endl;
  for (Int_t b=0;b<fgknBins;b++){ 
    //if(fBinNoOfParticles[b]==0)continue;
    //fAvPt[b]=fBinMeanPt[b]/fBinNoOfParticles[b];
    cout<<"pt bin: "<<b*fBinWidth<<"-"<<(b+1)*fBinWidth<<" GeV"<<endl;
    //cout<<"number of particles in this pt bin: "<<tempNo->GetBinContent(b+1,1,1)<<endl;
    //cout<<"mean pt in this bin: "<<fAvPt[b]<<" GeV"<<endl;
    if(fCumulant[0]>=0){
      fv2[b]=100.*fDiffCumulant2[b]/pow(fCumulant[0],.5);
      //if (fAvQ2x+fAvQ2y-pow(fAvQx,2.)-pow(fAvQy,2.)-pow(fV2*fAvM,2.)>0.){
       //fSddiff2[b]=pow((1./(2.*fBinNoOfParticles[b]))*((1.+pow(fChiQ[0],2.))/pow(fChiQ[0],2.)),0.5);
       cout<<"v'_2/2{2} = "<<fv2[b]<<"%, "<<" "<<"sd{2} = "<<100.*fSddiff2[b]<<"%"<<endl;
       //fDiffFlowResults2->SetBinContent(b+1,fv2[b],100.*fSddiff2[b]);
       fdfr2->SetBinContent(b+1,fv2[b]);
      //} else {
         //cout<<"v'_2/2{2} = Im"<<endl;
      //}
    }else{
      cout<<"v'_2/2{2} = Im"<<endl;
    } 
    if(fCumulant[1]<=0){
      fv4[b]=-100.*fDiffCumulant4[b]/pow(-fCumulant[1],.75);
      //if (fAvQ2x+fAvQ2y-pow(fAvQx,2.)-pow(fAvQy,2.)-pow(fV4*fAvM,2.)>0.){
       //fSddiff4[b]=pow((1./(2.*fBinNoOfParticles[b]))*((2.+6.*pow(fChiQ[1],2.)+pow(fChiQ[1],4.)+pow(fChiQ[1],6.))/pow(fChiQ[1],6.)),0.5);
       cout<<"v'_2/2{4} = "<<fv4[b]<<"%, "<<" "<<"sd{4} = "<<100.*fSddiff4[b]<<"%"<<endl;
       //fCommonHistsRes4->FillDifferentialFlow(b+1,fv4[b],100.*fSddiff4[b]);
       fdfr4->SetBinContent(b+1,fv4[b]);
      //} else {
       // cout<<"v'_2/2{4} = Im"<<endl;
      //} 
    }else{
      cout<<"v'_2/2{4} = Im"<<endl;
    }  
    if(fCumulant[2]>=0){
      cout<<"v'_2/2{6} = "<<100.*fDiffCumulant6[b]/(4.*pow((1./4.)*fCumulant[2],(5./6.)))<<"%"<<endl;
      fv6[b]=100.*fDiffCumulant6[b]/(4.*pow((1./4.)*fCumulant[2],(5./6.)));
      //fCommonHistsRes6->FillDifferentialFlow(b+1,fv6[b],0.);
      fdfr6->SetBinContent(b+1,fv6[b]);
    }else{
      cout<<"v'_2/2{6} = Im"<<endl;
    }     
    if(fCumulant[3]<=0){
      cout<<"v'_2/2{8} = "<<-100.*fDiffCumulant8[b]/(33.*pow(-(1./33.)*fCumulant[3],(7./8.)))<<"%"<<endl;
      fv8[b]=-100.*fDiffCumulant8[b]/(33.*pow(-(1./33.)*fCumulant[3],(7./8.))); 
      //fCommonHistsRes8->FillDifferentialFlow(b+1,fv8[b],0.);
      fdfr8->SetBinContent(b+1,fv8[b]);
    }else{
      cout<<"v'_2/2{8} = Im"<<endl;
    }       
    cout<<"****************************************"<<endl;
  }  
}

  













