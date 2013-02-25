// Class AliAnalysisTaskJetShape
// Author martin.poghosyan@cern.ch
//

#include "TChain.h"
#include "TTree.h"
#include "TH1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TCanvas.h"
#include "TParticlePDG.h"
#include "TParticle.h"
#include "TRandom.h"
#include "TDatabasePDG.h"
#include <TPDGCode.h>
#include "THnSparse.h"
#include "TGraph.h" 
#include "TArrayF.h"
#include "TArrayD.h"
#include <TVector3.h>
#include <TClonesArray.h>
#include "AliLog.h"

#include "AliCDBManager.h"
#include "AliGeomManager.h"
#include "AliCDBEntry.h"
#include "AliCDBPath.h"
#include "AliAlignObjParams.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliMultiplicity.h"
#include "AliESDFMD.h"
#include "AliESDVZERO.h"
#include "AliESDInputHandler.h"
#include "AliVEvent.h"
#include "AliCentrality.h"
#include "AliAnalysisHelperJetTasks.h"
#include "AliInputEventHandler.h"
#include "AliAODJetEventBackground.h"
#include "AliAODEvent.h"
#include "AliAODHandler.h"
#include "AliAODJet.h"
#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "AliAODMCParticle.h"
#include "AliGenEventHeader.h"
#include "AliHeader.h"
#include "AliStack.h"
#include "AliMCEvent.h"
#include "AliMCEventHandler.h"
#include "AliInputEventHandler.h"
#include "AliLHCData.h"
#include "AliGenPythiaEventHeader.h"
#include "AliGenDPMjetEventHeader.h"
#include "AliESDtrackCuts.h"
#include "AliTriggerAnalysis.h"



#include "AliAnalysisTaskJetShape.h"

/////////////////////////////////////

ClassImp(AliAnalysisTaskJetShape)


  AliAnalysisTaskJetShape::AliAnalysisTaskJetShape(const char *name) 
: AliAnalysisTaskSE(name), 
  //  fOutputTree(0),
  fOutputList(0),
  fESD(0x0),
  fAODIn(0x0),
  fAODOut(0x0),
  fAODExtension(0x0),
  fkMC(0),       
  fCMSE(0),      
  fRunNb(0),       
  fkIsPhysSel(0),  
  fNonStdFile(""),
  fFilterMask(0),
  fOfflineTrgMask(AliVEvent::kAny),
  fCentMin(0.),
  fCentMax(80.),
  fEvtClassMin(0),
  fEvtClassMax(4),
  fPtJcorrMin(20),
  fPtJBcorrMin(20),
  fJpPtmin(1),
  fJaPtmin(0.5),
  fVtxMinContrib(1),
  fVtxZMin(-10),
  fVtxZMax(10),
  fkIsPbPb(0),
  fhPtJL(0),
  fhAreaJL(0),
  fanalJetSubStrHM(0),
  fanalJetSubStrMCHM(0)
{


  fgkbinCent[0] = -0.099;
  fgkbinCent[1] = 5;
  fgkbinCent[2] = 10;
  fgkbinCent[3] = 20;
  fgkbinCent[4] = 40;
  fgkbinCent[5] = 60;
  fgkbinCent[6] = 80;
  fgkbinCent[7] =100;

  fBackgroundBranch[0] = "",
  fBackgroundBranch[1] = "",
  fJetBranchName[0] = "";
  fJetBranchName[1] = "";

  for(Int_t i=0; i<3; i++)
    {
      fhPtresVsPt[i]=0;
      fhPhiresVsPhi[i]=0;
      fhEtaresVsEta[i]=0;
      fhDCAxy[i]=0;
      fhDCAz[i]=0;
      fhTrackPtEtaPhi[i]=0;
    }
  //  fkIsPbPb = kFALSE;
  //  fDebug=0;



  fanalJetSubStrHM   = new AliAnalysisTaskJetShapeHM("rec", kFALSE);
  fanalJetSubStrMCHM = new AliAnalysisTaskJetShapeHM("truth", kTRUE);

  DefineOutput(1, TList::Class());

  //   DefineOutput(1, TTree::Class());
}


AliAnalysisTaskJetShape::~AliAnalysisTaskJetShape()
{
   if (fOutputList && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()) {
    printf("Deleteing output\n");

    if(fOutputList){
      delete fOutputList;
      fOutputList = 0;
    }

    //    if(fTriggerAnalysis) 
    //      delete fTriggerAnalysis;

   } 
}


AliAnalysisTaskJetShape::AliAnalysisTaskJetShapeTool::AliAnalysisTaskJetShapeTool():TObject(),
  fvecJ(0,0,0),
  fvecB1(0,0,0),
  fvecB2(0,0,0),
  fRmax(0),
  fPRecInRJ(0,0,0),
  fList(0),
  fEntries(0)
{
  fList = new TClonesArray("TVector3",10000);

  TVector3 v(0,0,0);
  SetVecJ(v);
  fRmax = -0.5;
  fR[0] =0.1;
  fR[1] =0.2;
  fR[2] =0.3;
  fEntries=0;

  for(Int_t i1=0; i1<fgkbinR; i1++) 
    {
      for(Int_t i2=0; i2<2; i2++)
	{
	  fListJ[i1][i2].Set(1000); 
	  fListB1[i1][i2].Set(1000);
	  fListB2[i1][i2].Set(1000);
	  fListJc[i1][i2]  = 0; 
	  fListB1c[i1][i2] = 0;
	  fListB2c[i1][i2] = 0;
	  //	  fListJ[i1][i2]("TVector3",10000); 
	  //	  fListB1[i1][i2] = TClonesArray("TVector3",10000);
	  //	  fListB2[i1][i2] = TClonesArray("TVector3",10000);

	 	  fListJ[i1][i2].Reset(-1); 
	  	  fListB1[i1][i2].Reset(-1);
	  	  fListB2[i1][i2].Reset(-1);

		  fPRecJ[i1][i2].SetXYZ(0,0,0);

	}
    }

  fPRecInRJ.SetXYZ(0,0,0);

  fPtMinTr[0] = 2;
  fPtMinTr[1] = 0.5;


}

AliAnalysisTaskJetShape::AliAnalysisTaskJetShapeTool::AliAnalysisTaskJetShapeTool(TClonesArray *list):TObject(),
  fvecJ(0,0,0),
  fvecB1(0,0,0),
  fvecB2(0,0,0),
  fRmax(0),
  fPRecInRJ(0,0,0),
  fList(list),
  fEntries(0)
{

  TVector3 v(0,0,0);
  SetVecJ(v);
  fRmax = 0.4;
  fR[0] =0.1;
  fR[1] =0.2;
  fR[2] =0.3;
  fEntries=0;

  for(Int_t i1=0; i1<fgkbinR; i1++) 
    {
      for(Int_t i2=0; i2<2; i2++)
	{
	  fListJ[i1][i2].Set(1000); 
	  fListB1[i1][i2].Set(1000);
	  fListB2[i1][i2].Set(1000);
	  fListJc[i1][i2]  = 0; 
	  fListB1c[i1][i2] = 0;
	  fListB2c[i1][i2] = 0;
	  //	  fListJ[i1][i2]("TVector3",10000); 
	  //	  fListB1[i1][i2] = TClonesArray("TVector3",10000);
	  //	  fListB2[i1][i2] = TClonesArray("TVector3",10000);
	 	  fListJ[i1][i2].Reset(-1); 
	  	  fListB1[i1][i2].Reset(-1);
	  	  fListB2[i1][i2].Reset(-1);
		  fPRecJ[i1][i2].SetXYZ(0,0,0);

	}
    }

  fPRecInRJ.SetXYZ(0,0,0);

  fPtMinTr[0] = 2;
  fPtMinTr[1] = 0.5;

}




AliAnalysisTaskJetShape::AliAnalysisTaskJetShapeTool::~AliAnalysisTaskJetShapeTool()
{
  fList->Delete();
  fEntries=0;

  //  if(fList)
  //    delete fList;

}


void AliAnalysisTaskJetShape::AliAnalysisTaskJetShapeTool::SetVecJ(TVector3 v)
{
  fvecJ.SetXYZ(v.X(), v.Y(), v.Z());
  fvecB1.SetXYZ(-v.Y(), v.X(), v.Z());
  fvecB2.SetXYZ(v.Y(),-v.X(), v.Z());
}


void AliAnalysisTaskJetShape::AliAnalysisTaskJetShapeTool::Clean()
{
  //  printf("AnalJetSubStrTool::Clean()\n");
  /*
  fList->Delete();
  for(Int_t i1=0; i1<fgkbinR; i1++) 
    {
      for(Int_t i2=0; i2<2; i2++)
	{
	  fListJ[i1][i2]->Delete(); 
	  fListB1[i1][i2]->Delete();
	  fListB2[i1][i2]->Delete();
	}
    }

  */

  //  fList.SetOwner();
  //  fList->Delete();
  //  fEntries=0;
  fPRecInRJ.SetXYZ(0,0,0);

  for(Int_t i1=0; i1<fgkbinR; i1++) 
    {
      for(Int_t i2=0; i2<2; i2++)
	{
	 	  fListJ[i1][i2].Reset(-1); 
	  	  fListB1[i1][i2].Reset(-1);
	  	  fListB2[i1][i2].Reset(-1);
	  //	  fListJ[i1][i2].Set(0); 
	  //	  fListB1[i1][i2].Set(0);
	  //	  fListB2[i1][i2].Set(0);
	  fListJc[i1][i2]  = 0; 
	  fListB1c[i1][i2] = 0;
	  fListB2c[i1][i2] = 0;

	  fPRecJ[i1][i2].SetXYZ(0,0,0);
	}
    }
  
}


//void AnalJetSubStrTool::Print(Option_t* /*option*/) const
void AliAnalysisTaskJetShape::AliAnalysisTaskJetShapeTool::PRINT() const
{

  //	  PRINT(fList, "all"); 

	  //  fList->Print();
  for(Int_t i1=0; i1<fgkbinR; i1++) 
    {
 
      if(i1!=0) continue;
     for(Int_t i2=0; i2<2; i2++)
	{

	  //	  printf("# %d   %d\n",i1, i2);
	  PRINT(fListJ[i1][i2], fListJc[i1][i2], Form("J_%d_%d",i1,i2)); 
		  //	  	  PRINT(fListB1[i1][i2], Form("B1_%d_%d",i1,i2)); 
		  //	  	  PRINT(fListB2[i1][i2], Form("B2_%d_%d",i1,i2)); 
	  //	  fListJ[i1][i2].Print("",1); 
	  //	  fListB1[i1][i2]->Print();
	  //	  fListB2[i1][i2]->Print();
	}
    }

}

void AliAnalysisTaskJetShape::AliAnalysisTaskJetShapeTool::PRINT(TArrayI a, Int_t n, const char *txt) const
{
  printf("%s :%d \n",txt, n);
  Int_t count = 0;
  for(Int_t i=0; i<n; i++){
    printf("%4d   ", a.At(i));

    if(count==20) 
      {
	printf("\n");
	count=0;
      }
    else count++;
  }
  if(count!=20) printf("\n");

}



 void AliAnalysisTaskJetShape::AliAnalysisTaskJetShapeTool::Init()
{
  Int_t Nev = fEntries;

  //    PRINT();

  //  printf("Nev = %d\n",Nev);

  for(Int_t iP=0; iP<Nev; iP++) 
    {
      TVector3 *v = (TVector3*) fList->At(iP);
      if(!v) {
	printf("ERROR : entry #%d not found\n",iP);
	continue;
      }

      //      printf("#%3d   ",iP); v->Print();

      Double_t R = CalcR(fvecJ, *v);
      //     printf("R = %f\n",R);
      if(R<fRmax)
	{
	  for(Int_t iT = 0; iT < fgkbinR; iT++)
	    {
	      Int_t type = 0;
	      if(R>fR[iT]) type = 1;

	      if(v->Pt() < fPtMinTr[type]) continue;
	      fPRecJ[iT][type]+=*v;
	      AddToJ(iT, type, iP);
	    }

	  if(v->Pt() < fPtMinTr[0]) continue;
	  fPRecInRJ+=*v;

	  continue;
	}

      R = CalcR(fvecB1, *v);
      if(R<fRmax)
	{
	  for(Int_t iT = 0; iT < fgkbinR; iT++)
	    {
	      Int_t type = 0;
	      if(R>fR[iT]) type = 1;

	      if(v->Pt() < fPtMinTr[type]) continue;

	      AddToB1(iT, type, iP);
	    }
	  continue;
	}

      R = CalcR(fvecB2, *v);
      if(R<fRmax)
	{
	  for(Int_t iT = 0; iT < fgkbinR; iT++)
	    {
	      Int_t type = 0;
	      if(R>fR[iT]) type = 1;

	      if(v->Pt() < fPtMinTr[type]) continue;

	      AddToB2(iT, type, iP);
	    }
	  continue;
	}
    }



  /*
  for(Int_t i1=0; i1<fgkbinR; i1++) 
    {
      for(Int_t i2=0; i2<2; i2++)
	{
	  //	  if(fListJc[i1][i2]<2) fListJc[i1][i2]=0;
	  //	  if(fListB1c[i1][i2]<2) fListB1c[i1][i2]=0;
	  //	  if(fListB2c[i1][i2]<2) fListB2c[i1][i2]=0;
	  fListJ[i1][i2].Set(fListJc[i1][i2]); 
	  fListB1[i1][i2].Set(fListB1c[i1][i2]);
	  fListB2[i1][i2].Set(fListB2c[i1][i2]);
	}
    }
  */

}



Double_t AliAnalysisTaskJetShape::AliAnalysisTaskJetShapeTool::CalcR(TVector3 v1, TVector3 v2)
{

  Double_t dphi = CalcdPhi(v1.Phi(), v2.Phi());
  //  dphi*=(0.9/TMath::Pi());
  Double_t deta = v1.Eta() - v2.Eta();
  Double_t RB = TMath::Sqrt(dphi*dphi+deta*deta);
  return RB;
}


Double_t AliAnalysisTaskJetShape::AliAnalysisTaskJetShapeTool::CalcdPhi(Double_t phi1, Double_t phi2)
{

  while(phi1<0) phi1+=TMath::TwoPi();
  while(phi1>TMath::TwoPi()) phi1-=TMath::TwoPi();

  while(phi2<0) phi2+=TMath::TwoPi();
  while(phi2>TMath::TwoPi()) phi2-=TMath::TwoPi();

  Double_t dphi = phi1- phi2;
  if(dphi>TMath::Pi())dphi = dphi - 2.*TMath::Pi();
  if(dphi<(-1.*TMath::Pi()))dphi = dphi + 2.*TMath::Pi();

  return dphi;
}


void AliAnalysisTaskJetShape::AliAnalysisTaskJetShapeTool::GetLocalMom(TVector3 vecJ, TVector3 *q)
{
  // theta and phi of a particle in loc.syst of fvecJ 


  Double_t PT = vecJ.Perp();
  Double_t P0 = vecJ.Mag();

  Double_t q0=q->X();
  Double_t q1=q->Y();
  Double_t q2=q->Z();

  Double_t qx1 = vecJ(2)/P0/PT*(vecJ(0)*q0+vecJ(1)*q1) - PT/P0*q2;
  Double_t qy1 = -vecJ(1)/PT*q0 + vecJ(0)/PT*q1;
  Double_t qz1 = 1/P0*(vecJ(0)*q0+vecJ(1)*q1+vecJ(2)*q2);

  //  Double_t q00 = TMath::Sqrt(qx1*qx1 + qy1*qy1 + qz1*qz1);
  //  phi = TMath::ATan2(qy1, qx1);
  //  if(phi<0) phi+=fTwoPi;
  //  theta = TMath::ACos(qz1/q00);

  q->SetXYZ(qx1, qy1, qz1);

  return;

}

/*

Bool_t AnalJetSubStrTool::FindCorrelationAxesAnd(TArrayI list, TVector3 vec, Double_t &Phi, Int_t scenario)
{
  Double_t TwoPi = TMath::TwoPi();

//
// 1st track loop to determine the sphericity matrix
//   
      // Initialize Shericity matrix
      Float_t mxx    = 0.;
      Float_t myy    = 0.;
      Float_t mxy    = 0.;
      Int_t   nc     = 0;
      Float_t sump2  = 0.;
      Float_t ptMax  = 0.;
      Float_t etaMax = 0.;
      Float_t phiMax = 0.;
      Int_t   iMax   = -1;
      Float_t ptsum  = 0.;
 

  Int_t Nev =  list.GetSize();
  if( Nev <2) return kFALSE;

  Int_t indexpmax = -1;
  Double_t pmax = -1;
  Double_t phipmax = 0;

  Int_t indexpzmax = -1;
  Double_t pzmax = -1;
  Double_t phipzmax = 0;


  Int_t indexthetamin = -1;
  Double_t thetamin = 1000;
  Double_t phithetamin = 0;

  for(Int_t ip=0; ip< Nev; ip++) {

    TVector3* part = (TVector3*) GetAt(list.At(ip));
    if(!part) continue;

    TVector3 vtmp(*part);
    GetLocalMom(vec, &vtmp);

    Float_t ppjX = vtmp.X();
    Float_t ppjY = vtmp.Y();
    Float_t ppjT = TMath::Sqrt(ppjX * ppjX + ppjY * ppjY);


    Float_t pt = ppjT;//part->Z();
    Float_t deta = part->Eta();
    Float_t dphi = part->Phi();

	      mxx += (ppjX * ppjX / ppjT);
	      myy += (ppjY * ppjY / ppjT);
	      mxy += (ppjX * ppjY / ppjT);
	      nc++;
	      sump2 += ppjT;


    if(vtmp.Mag()>pmax)
      {
	pmax=vtmp.Mag();
	indexpmax = ip;
	phipmax=vtmp.Phi();
      }

    if(vtmp.Z()>pzmax)
      {
	pzmax=vtmp.Z();
	indexpzmax = ip;
	phipzmax=vtmp.Phi();
      }
    if(vtmp.Theta()<thetamin)
      {
	thetamin=vtmp.Theta();
	indexthetamin = ip;
	phithetamin=vtmp.Phi();
      }


	  }
  


    
//
// At this point we have mxx, myy, mxy
  if (nc == 0) return kFALSE;      
// Shericity Matrix	
      const Double_t ele[4] = {mxx / sump2, mxy / sump2, mxy / sump2, myy / sump2};	
      TMatrixDSym m0(2, ele);
// Find eigenvectors
      TMatrixDSymEigen m(m0);
      TVectorD eval(2);
      TMatrixD evecm = m.GetEigenVectors();
      eval  = m.GetEigenValues();
// Largest eigenvector
      Int_t jev = 0;
      if (eval[0] < eval[1]) jev = 1;
      TVectorD evec0(2);
// Principle axis
      evec0 = TMatrixDColumn(evecm, jev);
      TVector2 evec(evec0[0], evec0[1]); 
// Principle axis from leading partice

//      Float_t phiM = TMath::ATan2(phiMax, etaMax);
//      TVector2 evecM(TMath::Cos(phiM), TMath::Sin(phiM)); 
//      Float_t phistM = evecM.DeltaPhi(evec);
//      if (TMath::Abs(phistM) > TMath::Pi()/2.) evec*=(-1.);


      Double_t phiTA = evec.Phi();


      Double_t phiDir = 0;
  if(scenario ==0)
    {
      phiDir = phipmax;
    }
  else if(scenario ==1)
    {
      phiDir = phipzmax;
    }

  else if(scenario ==2)
    {
      phiDir = phithetamin;
    }
  else 
    return kFALSE;


      Phi = phiTA;

      if( TMath::Abs(CalcdPhi(phiDir, phiTA)) <TMath::Pi()/2)
	return kTRUE;
      else {
	Phi+=TMath::Pi();
	if(Phi>TwoPi) Phi-=TwoPi;
	return kTRUE;
      }


 return kTRUE;
} 
*/

Bool_t AliAnalysisTaskJetShape::AliAnalysisTaskJetShapeTool::FindCorrelationAxes(TArrayI list, TVector3 vec, Double_t &Phi, Int_t scenario)
{
  // scenario 0- Pmax, 1-Pzmax_Local, 2->cosMin, 3-PTmac, 4->PTmax_Local

  Double_t xmean , ymean;

  Double_t TwoPi = TMath::TwoPi();

  xmean = 0;
  ymean = 0;
  Phi= -999;
  //find axes

  Double_t x=0;
  Double_t x2=0;
  Double_t y=0; 
  Double_t y2=0;
  Double_t xy=0;
  
  //  Double_t phiLoc, thetaLoc;

  Int_t N=0;

  Int_t Nev =  list.GetSize();
  if( Nev <2) return kFALSE;




  Int_t Index = -1;
  Double_t val = -1;
  if(scenario==2) val = 100;
  Double_t phiDir = -99;



  for(Int_t ip=0; ip< Nev; ip++) {
    if(list.At(ip)<0) break;

    TVector3* part = (TVector3*) GetAt(list.At(ip));
    if(!part) continue;


    Double_t xx, yy;
    /*
    TVector3 vtmp(*part);
    GetLocalMom(vec, &vtmp);
    Double_t xx=vtmp.X(); //part->Mag()*TMath::Sin(thetaLoc)*TMath::Cos(phiLoc);
    Double_t yy=vtmp.Y(); //part->Mag()*TMath::Sin(thetaLoc)*TMath::Sin(phiLoc);
    */

    xx = CalcdPhi(part->Phi(), vec.Phi());
    yy = part->Eta() - vec.Eta();

    x+=xx;
    y+=yy;
    x2+=(xx*xx);
    y2+=(yy*yy);
    xy+=(xx*yy);
    N++;

    /*
  if(scenario ==0)
    {
      if(vtmp.Mag()>val)
	{
	  val=vtmp.Mag();
	  Index = ip;
	  phiDir=vtmp.Phi();
	}
    }
  else if(scenario ==1)
    {
      if(vtmp.Z()>val)
	{
	  val=vtmp.Z();
	  Index = ip;
	  phiDir=vtmp.Phi();
	}
    }
  else if(scenario ==2)
    {
      if(vtmp.Theta()<val)
	{
	  val=vtmp.Theta();
	  Index = ip;
	  phiDir=vtmp.Phi();
	}
    }
  else 
*/
if(scenario ==3)
    {
      if(part->Pt()> val)
	{
	  val =part->Pt();
	  Index = ip;
	  //	  phiDir=vtmp.Phi();
	  phiDir=TMath::ATan2(yy, xx); //vtmp.Phi();
	}
    }


  }

  if(N<2) return kFALSE;

  //  return kFALSE;

  x/=N;
  y/=N;
  x2/=N;
  y2/=N;
  xy/=N;

  
  Double_t A = xy-x*y;
  Double_t B = x2-x*x+y*y-y2;
  Double_t D = TMath::Sqrt(B*B+4*A*A);

  //  printf("N = %d  A = %f\n",N, A);
  //  list.Print();

  Double_t a1 = (-B+D)/2/A;
  //  Double_t a2 = (-B-D)/2/A;
  //  Double_t b1 = y-a1*x;
  //  Double_t b2 = y-a2*x;
  //  Double_t phiDir = TMath::ATan2(y, x);
  // while(phiDir<0) phiDir+=TwoPi;
    Double_t phi = TMath::ATan(a1);
    
  while(phi>TwoPi) phi-=TwoPi;
  while(phi< 0 )   phi+=TwoPi;

  Phi=CalcdPhi0To2pi(phi, 0);

  if(scenario!=4 && Index<0) return kFALSE;

  if(Index>=0)
    {
      if( TMath::Abs(CalcdPhi(phiDir, phi)) <TMath::Pi()/2)
	return kTRUE;
      else
	{
	  phi+=TMath::Pi();
	  Phi = CalcdPhi0To2pi(phi);
	  return kTRUE;
	}
    }


  Double_t xmax = -100;
  Double_t xmin =  100;

  xmean = x;
  ymean = y;


  for(Int_t ip=0; ip< Nev; ip++) {
    if(list.At(ip)<0) break;

    TVector3* part = (TVector3*) GetAt(list.At(ip));
    if(!part) continue;

    TVector3 vtmp(*part);
    GetLocalMom(vec, &vtmp);

    Double_t xx=vtmp.X(); //part->Mag()*TMath::Sin(thetaLoc)*TMath::Cos(phiLoc);
    Double_t yy=vtmp.Y(); //part->Mag()*TMath::Sin(thetaLoc)*TMath::Sin(phiLoc);

    Double_t x1 = (xx-xmean)*TMath::Cos(Phi) +(yy-ymean)*TMath::Sin(Phi);
    if(x1 > xmax) xmax = x1;
    if(x1 < xmin) xmin = x1;
    //    printf("c3\n");
  }
   
  if(TMath::Abs(xmin) > xmax) Phi+=TMath::Pi();
 



   while(Phi>TwoPi) Phi-=TwoPi;

  return kTRUE;
}




Bool_t AliAnalysisTaskJetShape::AliAnalysisTaskJetShapeTool::FindCorrelationAxesCorr(TArrayI list, TVector3 vec, Double_t &Phi, Int_t scenario, Double_t R)
{
  // scenario 0- Pmax, 1-Pzmax_Local, 2->cosMin, 3-PTmac, 4->PTmax_Local

  Double_t xmean , ymean;

  Double_t TwoPi = TMath::TwoPi();

  xmean = 0;
  ymean = 0;
  Phi= -999;
  //find axes

  Double_t x=0;
  Double_t x2=0;
  Double_t y=0; 
  Double_t y2=0;
  Double_t xy=0;
  
  //  Double_t phiLoc, thetaLoc;

  Int_t N=0;

  Int_t Nev =  list.GetSize();
  if( Nev <2) return kFALSE;




  Int_t Index = -1;
  Double_t val = -1;
  if(scenario==2) val = 100;
  Double_t phiDir = -99;



  for(Int_t ip=0; ip< Nev; ip++) {
    if(list.At(ip)<0) break;

    TVector3* part = (TVector3*) GetAt(list.At(ip));
    if(!part) continue;

    Double_t xx, yy;
    /*
    TVector3 vtmp(*part);
    GetLocalMom(vec, &vtmp);
    Double_t xx=vtmp.X(); //part->Mag()*TMath::Sin(thetaLoc)*TMath::Cos(phiLoc);
    Double_t yy=vtmp.Y(); //part->Mag()*TMath::Sin(thetaLoc)*TMath::Sin(phiLoc);
    */

    xx = CalcdPhi(part->Phi(), vec.Phi());
    yy = part->Eta() - vec.Eta();

    x+=xx;
    y+=yy;
    x2+=(xx*xx);
    y2+=(yy*yy);
    xy+=(xx*yy);
    N++;

    /*
  if(scenario ==0)
    {
      if(vtmp.Mag()>val)
	{
	  val=vtmp.Mag();
	  Index = ip;
	  phiDir=vtmp.Phi();
	}
    }
  else if(scenario ==1)
    {
      if(vtmp.Z()>val)
	{
	  val=vtmp.Z();
	  Index = ip;
	  phiDir=vtmp.Phi();
	}
    }
  else if(scenario ==2)
    {
      if(vtmp.Theta()<val)
	{
	  val=vtmp.Theta();
	  Index = ip;
	  phiDir=vtmp.Phi();
	}
    }
  else 
*/
if(scenario ==3)
    {
      if(part->Pt()> val)
	{
	  val =part->Pt();
	  Index = ip;
	  //	  phiDir=vtmp.Phi();
	  phiDir=TMath::ATan2(yy, xx); //vtmp.Phi();
	}
    }


  }

  if(N<2) return kFALSE;
  if(scenario!=4 && Index<0) return kFALSE;

  //  return kFALSE;

  x/=N;
  y/=N;
  x2/=N;
  y2/=N;
  xy/=N;

  
  Double_t A = xy-x*y;
  Double_t B = x2-x*x+y*y-y2;
  Double_t D = TMath::Sqrt(B*B+4*A*A);

  //  printf("N = %d  A = %f\n",N, A);
  //  list.Print();

  Double_t a1 = (-B+D)/2/A;
  //  Double_t a2 = (-B-D)/2/A;
  //  Double_t b1 = y-a1*x;
  //  Double_t b2 = y-a2*x;
  //  Double_t phiDir = TMath::ATan2(y, x);
  // while(phiDir<0) phiDir+=TwoPi;
    Double_t phi = TMath::ATan(a1);
    
  Phi=CalcdPhi0To2pi(phi, 0);
  if( TMath::Abs(CalcdPhi(phiDir, phi))  > TMath::Pi()/2)
    Phi=CalcdPhi0To2pi(phi+TMath::Pi(), 0);




  if(Index>=0)
    {
      if(N==2) 
	return kTRUE;



      Double_t phiMinus = 0;
      Int_t Nminus = 0;
      for(Int_t ip=0; ip< Nev; ip++) 
	{
	  if(list.At(ip)<0) break;

	  TVector3* part = (TVector3*) GetAt(list.At(ip));
	  if(!part) continue;
	  Double_t xx0 = CalcdPhi(part->Phi(), vec.Phi());
	  Double_t yy0 = part->Eta() - vec.Eta();

	  Double_t xx  = (x*N - xx0)/(N-1);
	  Double_t yy  = (y*N - yy0)/(N-1);
	  Double_t xx2 = (x2*N - xx0*xx0)/(N-1);
	  Double_t yy2 = (y2*N - yy0*yy0)/(N-1);
	  Double_t xxyy= (xy*N - xx0*yy0)/(N-1);


	  Double_t AA = xxyy-xx*yy;
	  Double_t BB = xx2-xx*xx+yy*yy-yy2;
	  Double_t DD = TMath::Sqrt(BB*BB+4*AA*AA);
	  Double_t aa1 = (-BB+DD)/2/AA;
	  Double_t phi1 = TMath::ATan(aa1);
    
	  if( TMath::Abs(CalcdPhi(phiDir, phi1)) > TMath::Pi()/2)
	    phi1+=TMath::Pi();

	  phiMinus+=CalcdPhi(Phi, phi1);
	  Nminus++;

	}


      Double_t phiPlus = 0;
      Int_t Nplus = 0;
      for(Int_t ip=0; ip< -Nev; ip++) 
	{
	  if(list.At(ip)<0) break;

	  TVector3* part = (TVector3*) GetAt(list.At(ip));
	  if(!part) continue;
	  Double_t rtmp = gRandom->Uniform(0, R);
	  Double_t phitmp = gRandom->Uniform(0, TMath::TwoPi());
	  Double_t xx0 = rtmp*TMath::Cos(phitmp);
	  Double_t yy0 = rtmp*TMath::Sin(phitmp);

	  Double_t xx  = (x*N + xx0)/(N+1);
	  Double_t yy  = (y*N + yy0)/(N+1);
	  Double_t xx2 = (x2*N + xx0*xx0)/(N+1);
	  Double_t yy2 = (y2*N + yy0*yy0)/(N+1);
	  Double_t xxyy= (xy*N + xx0*yy0)/(N+1);


	  Double_t AA = xxyy-xx*yy;
	  Double_t BB = xx2-xx*xx+yy*yy-yy2;
	  Double_t DD = TMath::Sqrt(BB*BB+4*AA*AA);
	  Double_t aa1 = (-BB+DD)/2/AA;
	  Double_t phi1 = TMath::ATan(aa1);
    
	  if( TMath::Abs(CalcdPhi(phiDir, phi1)) > TMath::Pi()/2)
	    phi1+=TMath::Pi();

	  phiPlus+=CalcdPhi(Phi, phi1);
	  Nplus++;

	}


      Phi+=((phiMinus+phiPlus)/(Nminus+Nplus+1));

      if( TMath::Abs(CalcdPhi(phiDir, Phi)) > TMath::Pi()/2)
	Phi+=TMath::Pi();

      Phi=CalcdPhi0To2pi(Phi, 0);

      return kTRUE;
    }


  Double_t xmax = -100;
  Double_t xmin =  100;

  xmean = x;
  ymean = y;


  for(Int_t ip=0; ip< Nev; ip++) {
    if(list.At(ip)<0) break;

    TVector3* part = (TVector3*) GetAt(list.At(ip));
    if(!part) continue;

    TVector3 vtmp(*part);
    GetLocalMom(vec, &vtmp);

    Double_t xx=vtmp.X(); //part->Mag()*TMath::Sin(thetaLoc)*TMath::Cos(phiLoc);
    Double_t yy=vtmp.Y(); //part->Mag()*TMath::Sin(thetaLoc)*TMath::Sin(phiLoc);

    Double_t x1 = (xx-xmean)*TMath::Cos(Phi) +(yy-ymean)*TMath::Sin(Phi);
    if(x1 > xmax) xmax = x1;
    if(x1 < xmin) xmin = x1;
    //    printf("c3\n");
  }
   
  if(TMath::Abs(xmin) > xmax) Phi+=TMath::Pi();
 



   while(Phi>TwoPi) Phi-=TwoPi;

  return kTRUE;
}

/*
/////////////////////////////////////

 TH2F *fhPhiEtaTrack;//! 

 TH1F *fhTMA_JAA[3];//! 
 TH1F *fhTMA_JAp[3];//! 
 TH1F *fhTMA_B1AA[3];//! 
 TH1F *fhTMA_B1Ap[3];//! 
 TH1F *fhTMA_B2AA[3];//! 
 TH1F *fhTMA_B2Ap[3];//! 


 Int_t     fPtJetNbin;
 TArrayD   fPtJetArray;
 Int_t     fPsiVsRNbin;
 Double_t  fRmax;
 Int_t     fPsiVsPhiNbin;

 UInt_t   fFilterMask;
 Double_t fEtaTrackMax;
 Double_t fPtTrackMin;
 Double_t fPtTrackMax;
 Double_t fPtJmin;
 Double_t fPtJmax;
 Double_t fPtJ;

 TVector3 fJvec;
*/

AliAnalysisTaskJetShape::AliAnalysisTaskJetShapeHM::AliAnalysisTaskJetShapeHM(TString comment, Bool_t kMC):TObject(),
fComment(comment),
fkMC(kMC),
fkMCprod(0),
fhEvent(0),
fhMult(0),
fhPtJ(0),
fhPtJvsPtCorr(0),
fhPsiVsR(0),
fhPsiVsRPtJ(0),
fhPhiEtaTrack(0),
fhPsiVsR_MCtr(0),
fhPsiVsRPtJ_MCtr(0),
fhJetTrPtVsPartPt(0),
fPtJetNbin(0),
fPtJetArray(0),
fPsiVsRNbin(0),
fRmax(0),
fPsiVsPhiNbin(0),
fFilterMask(0),
fEtaTrackMax(0),
fPtTrackMin(0),
fPtTrackMax(0),
fPtJmin(0),
fPtJmax(0),
fPtJ(0),
fJvec(0,0,0)
{

  for(Int_t i=0; i<3; i++)
    {
      fhTMA_JAA[i]=0;
      fhTMA_JAp[i]=0;
      fhTMA_B1AA[i]=0;
      fhTMA_B1Ap[i]=0;
      fhTMA_B2AA[i]=0;
      fhTMA_B2Ap[i]=0;
    }

  for(Int_t i=0; i<3; i++)
    {
      for(Int_t j=0; j<2; j++)
	{
	  fhPtresVsPt[i][j]=0;
	  fhPhiresVsPhi[i][j]=0;
	  fhEtaresVsEta[i][j]=0;
	  fhRresVsPt[i][j]=0;
	  fhDCAxy[i][j]=0;
	  fhDCAz[i][j]=0;
	  fhTrackPtEtaPhi[i][j]=0;
	}
    }


  fPtJetNbin = 0;
  fPtJetArray = 0;

  SetRBins();
  SetPhiNbins();
  SetFilterMask();
  SetEtaTrackMax();
  SetPtTrackRange();
  SetPtJetRange();

  fJvec.SetXYZ(0,0,0);
  fPtJ = 0;
}


AliAnalysisTaskJetShape::AliAnalysisTaskJetShapeHM::~AliAnalysisTaskJetShapeHM()
{
  /*
  if(fPtJetArray)
    delete [] fPtJetArray;
  fPtJetArray=0;
  */
}


void AliAnalysisTaskJetShape::AliAnalysisTaskJetShapeHM::SetPtJetBins(Int_t Nbin, Double_t *array)
{

  /*
  if(fPtJetArray)
    delete [] fPtJetArray;
  fPtJetArray=0;

  fPtJetNbin = Nbin;
  fPtJetArray = new Double_t[fPtJetNbin+1];
  for(Int_t i=0; i<= fPtJetNbin; i++) fPtJetArray[i]= array[i]; 
  */


  fPtJetNbin = Nbin;
  fPtJetArray.Set(fPtJetNbin+1, array);
}


void AliAnalysisTaskJetShape::AliAnalysisTaskJetShapeHM::InitHistos()
{

  fhEvent    = Hist1D("hEvent" , 3  , -0.5,   2.5, "Event");
  fhMult     = Hist1D("hMult"  , 101, -0.5, 100.5, "N_{ch}");

  fhPhiEtaTrack = Hist2D("hPhiEtaTrack", 100, 0, TMath::TwoPi(), 20, -1, 1., "#phi", "#eta");

  if(fPtJetNbin<1) {
    Int_t Nbin = 13;
    Double_t array[]= {0., 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 130, 160, 200};
    SetPtJetBins(Nbin, array);
  }

  fhPtJ      = Hist1D("hPtJ"  , fPtJetNbin, fPtJetArray.GetArray(), "Pt_{J} GeV/c");
  //  fhPsiVsR   = Hist1D("hPsiVsR", fPsiVsRNbin, 0., fRmax, "R", 1, "#Psi(R)");
  fhPsiVsR   = Hist3D("hPsiVsR", fPsiVsRNbin, 0., fRmax, 100, 0., 1., fPtJetNbin, fPtJetArray.GetArray(), "R", "P_{Jt, frac}", "P_{J} GeV/c");
  fhPsiVsRPtJ   = Hist2D("hPsiVsRPtJ", fPsiVsRNbin, 0., fRmax, fPtJetNbin, fPtJetArray.GetArray(), "R", "P_{tJ} GeV/c", 1, "#Psi(R)");

  fhPtJvsPtCorr  = Hist2D("fhPtJvsPtCorr", fPtJetNbin, fPtJetArray.GetArray(), 100, -100, 100, "P_{tJ} GeV/c", "P_{tJ} - P_{tJ,corr} GeV/c" , 1);




  for(Int_t i=0; i<3; i++)
    {
      fhTMA_JAA[i]  = Hist1D(Form("fhTMA_JAA_%d",i), fPsiVsPhiNbin, 0, TMath::TwoPi(), "#phi_{R}", 1, "#Psi(#phi_{R})");
      fhTMA_JAp[i]  = Hist1D(Form("fhTMA_JAp_%d",i), fPsiVsPhiNbin, 0, TMath::TwoPi(), "#phi_{R}", 1, "#Psi(#phi_{R})");

      fhTMA_B1AA[i]  = Hist1D(Form("fhTMA_B1AA_%d",i), fPsiVsPhiNbin, 0, TMath::TwoPi(), "#phi_{R}", 1, "#Psi(#phi_{R})");
      fhTMA_B1Ap[i]  = Hist1D(Form("fhTMA_B1Ap_%d",i), fPsiVsPhiNbin, 0, TMath::TwoPi(), "#phi_{R}", 1, "#Psi(#phi_{R})");

      fhTMA_B2AA[i]  = Hist1D(Form("fhTMA_B2AA_%d",i), fPsiVsPhiNbin, 0, TMath::TwoPi(), "#phi_{R}", 1, "#Psi(#phi_{R})");
      fhTMA_B2Ap[i]  = Hist1D(Form("fhTMA_B2Ap_%d",i), fPsiVsPhiNbin, 0, TMath::TwoPi(), "#phi_{R}", 1, "#Psi(#phi_{R})");
    }

  if(fkMCprod)
    {
  fhPsiVsR_MCtr      = Hist3D("hPsiVsR_MCtr", fPsiVsRNbin, 0., fRmax, 100, 0., 1., fPtJetNbin, fPtJetArray.GetArray(), "R", "P_{Jt, frac}", "P_{J} GeV/c");
  //  fhPsiVsR_MCtr      = Hist1D("hPsiVsR_MCtr", fPsiVsRNbin, 0., fRmax, "R", 1, "#Psi(R)");
  fhPsiVsRPtJ_MCtr   = Hist2D("hPsiVsRPtJ_MCtr", fPsiVsRNbin, 0., fRmax, fPtJetNbin, fPtJetArray.GetArray(), "R", "P_{tJ} GeV/c", 1, "#Psi(R)");
  fhJetTrPtVsPartPt  = Hist2D("fhJetTrPtVsPartPt", fPtJetNbin, fPtJetArray.GetArray(), 100, -1, 1, "P_{tJ,part} GeV/c", "1-P_{tJ,tr}/P_{tJ,part} GeV/c" , 1);
  const char *ch[2]={"m","p"};
      for(Int_t i=0; i<3; i++)
	{
	  for(Int_t j=0; j<2; j++)
	    {
	  fhPtresVsPt[i][j]    = Hist2D(Form("hPtresVsPt%d%s",i,ch[j]), 100, 0, 100, 100, -0.5, 0.5, "P_{t}^{gen} GeV/c", "1-P_{t}^{rec}/P_{t}^{gen} GeV/c" );    
	  fhPhiresVsPhi[i][j]  = Hist2D(Form("hPhiresVsPhi%d%s",i,ch[j]), 600, 0, TMath::TwoPi(), 128, -0.2, 0.2, "#phi^{gen} rad", "#phi^{rec} - #phi^{gen} rad" );
	  fhEtaresVsEta[i][j]  = Hist2D(Form("hEtaresVsEta%d%s",i,ch[j]), 200, -1, 1, 40, -0.2, 0.2, "#eta^{gen}", "#eta^{rec} - #eta^{gen}" ); 
	  fhRresVsPt[i][j]     = Hist2D(Form("hRresVsPt%d%s",i,ch[j])  , 100, 0, 100, 500, -fRmax, fRmax, "P_{t, track}^{gen} GeV/c", "R^{gen}-R^{rec}" );    
	  fhDCAxy[i][j]        = Hist2D(Form("hDCAxy%d%s",i,ch[j]), 100, 0, 100, 300, -3, 3, Form("p_{t} of part. type %d%s [GeV/c]",i,ch[j]), "DCAxy [cm]");
	  fhDCAz[i][j]         = Hist2D(Form("hDCAz%d%s",i,ch[j]) , 100, 0, 100, 300, -3, 3, Form("p_{t} of part. type %d%s [GeV/c]",i,ch[j]), "DCAz  [cm]") ;
	  fhTrackPtEtaPhi[i][j]= Hist3D(Form("hTrackPtEtaPhi%d%s",i,ch[j]), 100, 0, 100, 100, -1, 1, 100, 0, TMath::TwoPi(),"P_{t} GeV/c ","#eta","#phi rad");
	    }
	}
    }

}





void AliAnalysisTaskJetShape::AliAnalysisTaskJetShapeHM::AddToList(TList *list)
{
  list->Add(fhEvent);
  list->Add(fhMult);
  list->Add(fhPhiEtaTrack);
  list->Add(fhPtJ);
  list->Add(fhPsiVsR);
  list->Add(fhPsiVsRPtJ);
  list->Add(fhPtJvsPtCorr);


  for(Int_t i=0; i<3; i++)
    {
      list->Add(fhTMA_JAA[i]);
      list->Add(fhTMA_JAp[i]);

      list->Add(fhTMA_B1AA[i]);
      list->Add(fhTMA_B1Ap[i]);

      list->Add(fhTMA_B2AA[i]);
      list->Add(fhTMA_B2Ap[i]);
    }

  if(fkMCprod)
    {
      list->Add(fhPsiVsR_MCtr);
      list->Add(fhPsiVsRPtJ_MCtr);
      list->Add(fhJetTrPtVsPartPt);

      for(Int_t i=0; i<3; i++)
	{
	  for(Int_t j=0; j<2; j++)
	    {
	      list->Add(fhPtresVsPt[i][j]);
	      list->Add(fhPhiresVsPhi[i][j]);
	      list->Add(fhEtaresVsEta[i][j]);
	      list->Add(fhRresVsPt[i][j]);
	      list->Add(fhDCAxy[i][j]);
	      list->Add(fhDCAz[i][j]);
	      list->Add(fhTrackPtEtaPhi[i][j]);
	    }
	}
    }
}




Bool_t AliAnalysisTaskJetShape::AliAnalysisTaskJetShapeHM::AddEvent(AliAODEvent* aodE,  AliAODJet *jet, Double_t DeltaPtJ)
{



  Int_t size = aodE->GetNumberOfTracks();

  TClonesArray *arrP = 0;

  if(fkMC || fkMCprod)
    {
      arrP = dynamic_cast<TClonesArray*>(aodE->FindListObject(AliAODMCParticle::StdBranchName()));
     if(!arrP) 
       {
	 printf("ERROR: no Info about particles in AOD\n");
	 return kFALSE;
       }
    }

  if(fkMC)
     size = arrP->GetEntriesFast();

  Int_t IndexArray[size];
  Int_t IndexArrayMC[size];

  TClonesArray farray("TVector3", size);

  Int_t counter = 0;
  Int_t counterMC = 0;
  Int_t counterAll = 0;
  Double_t pJ[3] = {0, 0, 0};
  Double_t pJmc[3] = {0, 0, 0};

  AliAODVertex* primVtx = aodE->GetPrimaryVertex();
  Double_t bfield = aodE->GetMagneticField();
  Double_t dca[2];
  Double_t cov[3];

  TVector3 vecJ(jet->Px(),jet->Py(),jet->Pz());




  for(int it = 0;it < size;++it){
    TVector3 vec;
    Int_t ch = -999;
    Int_t label = 0;

 
    if(fkMC)
      {
	AliAODMCParticle *part = dynamic_cast<AliAODMCParticle*>(arrP->At(it));
	if(!part)continue;
	if(!part->IsPhysicalPrimary())continue;
	if(part->Charge()==0)continue;
	vec.SetXYZ(part->Px(), part->Py(), part->Pz());
      }
    else 
      {
	AliAODTrack *tr = aodE->GetTrack(it);
	if(!tr) continue;
	if((fFilterMask>0)&&!(tr->TestFilterBit(fFilterMask)))continue;
	label = tr->GetLabel();
	AliAODTrack tmp(*tr);
	tmp.PropagateToDCA(primVtx, bfield, 5., dca, cov);
	vec.SetXYZ(tr->Px(), tr->Py(), tr->Pz());
	ch = tr->Charge();
      }


    if(TMath::Abs(vec.Eta())>fEtaTrackMax)continue;
    if(vec.Pt()< fPtTrackMin)continue;
    if(vec.Pt()> fPtTrackMax) {return kFALSE;}


    new(farray[counterAll]) TVector3(vec);
    counterAll++;


    Double_t R = CalcR(vecJ, vec);
    if(R> fRmax) continue;

    pJ[0]+=vec.Px();
    pJ[1]+=vec.Py();
    pJ[2]+=vec.Pz();

    IndexArray[counter] = it;
    counter++;


    
    // effics
    if(fkMCprod)
      {
	//	  printf("A1\n");
	AliAODMCParticle *part = dynamic_cast<AliAODMCParticle*>(arrP->At(TMath::Abs(label)));
	if(!part)continue;

	Int_t type = 0;
	if(!part->IsPhysicalPrimary()) type=1;
	if(label <0) type=2;

	if(!(ch==-1 || ch==1)) AliFatal("ch != +/- 1!!!\n\n");
	  //	  printf("A4\n");
	if(ch==-1) ch=0;

	Double_t phip = AliAnalysisTaskJetShapeTool::CalcdPhi0To2pi(part->Phi());
	Double_t dphi = AliAnalysisTaskJetShapeTool::CalcdPhiSigned(part->Phi(), vec.Phi());
	Double_t phiT = AliAnalysisTaskJetShapeTool::CalcdPhi0To2pi(vec.Phi());

	fhPtresVsPt[type][ch]->Fill(part->Pt(), 1-vec.Pt()/part->Pt());
	fhPhiresVsPhi[type][ch]->Fill(phip, dphi);
	fhEtaresVsEta[type][ch]->Fill(part->Eta(), vec.Eta() - part->Eta());
	fhDCAxy[type][ch]->Fill(part->Pt(), dca[0]);
	fhDCAz[type][ch]->Fill( part->Pt(), dca[1]);
	fhTrackPtEtaPhi[type][ch]->Fill(vec.Pt(), vec.Eta(), phiT);

	TVector3 vecP(part->Px(), part->Py(), part->Pz());
	Double_t Rgen = CalcR(vecJ, vecP);
	fhRresVsPt[type][ch]->Fill(part->Pt(), Rgen - R);
        
	pJmc[0]+=part->Px();
	pJmc[1]+=part->Py();
	pJmc[2]+=part->Pz();

	IndexArrayMC[counterMC] = label;
	counterMC++;
	  //	  printf("A5\n");
      }
  
  }

  fhMult->Fill(counter);
  if(counter<1) return kFALSE;

  fJvec.SetXYZ(pJ[0], pJ[1], pJ[2]);


  fPtJ =  TMath::Sqrt(pJ[0]*pJ[0] + pJ[1]*pJ[1]);
  if((fPtJ < fPtJmin) || (fPtJ > fPtJmax)) return kFALSE;
  fhPtJ->Fill(fPtJ);





  fhPtJvsPtCorr->Fill(fPtJ, jet->Pt() - DeltaPtJ);

  for(Int_t it = 0; it<counter; it++)
    {
      TVector3 vec;

      if(fkMC)
	{
	  AliAODMCParticle *part = dynamic_cast<AliAODMCParticle*>(arrP->At(IndexArray[it]));
	  if(!part) continue;
	  vec.SetXYZ(part->Px(), part->Py(), part->Pz());
	}
      else 
	{
	  AliAODTrack *tr = aodE->GetTrack(IndexArray[it]);
	  if(!tr) continue;
          AliAODTrack tmp(*tr);
          tmp.PropagateToDCA(primVtx, bfield, 5., dca, cov);
	  if(tr)vec.SetXYZ(tmp.Px(), tmp.Py(), tmp.Pz());
	}

	Double_t R = CalcR(fJvec, vec);
	//	Double_t pt = (tr->Px()*pJ[0] + tr->Py()*pJ[1])/PtJ;
	Double_t probL = fJvec.Dot(vec)/fJvec.Mag2();
	//	fhPsiVsR->Fill(R, probL);
	//	fhPsiVsRPtJ->Fill(R, fPtJ, probL);
	fhPsiVsR->Fill(R,probL, fJvec.Mag());
	fhPsiVsRPtJ->Fill(R, fPtJ);
 
	Double_t phi = AliAnalysisTaskJetShapeTool::CalcdPhi0To2pi(vec.Phi());
	fhPhiEtaTrack->Fill(phi, vec.Eta());

    }


  if(fkMCprod)
    {
      TVector3 fJvecMCtr(pJmc[0], pJmc[1], pJmc[2]);
      Double_t fPtJMCtr=  fJvecMCtr.Perp();

      fhJetTrPtVsPartPt->Fill(fPtJMCtr, 1-fPtJ/fPtJMCtr);
      for(Int_t it = 0; it<counterMC; it++)
	{
	  TVector3 vec;

	  AliAODMCParticle *part = dynamic_cast<AliAODMCParticle*>(arrP->At(TMath::Abs(IndexArrayMC[it])));
	  if(!part) continue;
	  vec.SetXYZ(part->Px(), part->Py(), part->Pz());

	  Double_t R = CalcR(fJvecMCtr, vec);
	  Double_t probL = fJvecMCtr.Dot(vec)/fJvecMCtr.Mag2();
	  //	  fhPsiVsR->Fill(R, probL);
	  //	  fhPsiVsRPtJ->Fill(R, fPtJ, probL);
	  //Double_t probL = fJvecMCtr.Dot(vec)/fJvecMCtr.Mag2();
	  fhPsiVsR_MCtr->Fill(R,probL, fJvecMCtr.Mag());
	  fhPsiVsRPtJ_MCtr->Fill(R, fPtJMCtr);
 	}
    }

  // ang. distr.
  //     fMyTool->Clean();
  AliAnalysisTaskJetShapeTool *fMyTool = new AliAnalysisTaskJetShapeTool(&farray);

      fMyTool->SetNEntries(counterAll);
      fMyTool->SetVecJ(vecJ);
      fMyTool->SetPtMinTr(fPtTrackMin,fPtTrackMin);
      fMyTool->Init();
      Int_t scenario = 3;

      /*
       to be added!!!!!!!!!
      fhTMA_B1AA[i]=0;
      fhTMA_B2AA[i]=0;
      */

      for(Int_t l=0; l<3; l++)
	{
	  Double_t phiA, phi1;

	  //	  Double_t ptRJ0= fMyTool->GetPRecInRJ().Pt();
	  //	  Double_t ptRJ  = fMyTool->GetPRecJ(l,0).Pt();
	  Int_t N1 = fMyTool->GetSizeJ(l,1);
	  Double_t dphi = -999;


	  if(!fMyTool->FindCorrelationAxes(fMyTool->GetListJ(l, 0)  , vecJ, phiA, scenario))
	    continue;

	  //	  f2JEvent->SetPhiJL(l,0, phiA);

	  if(fMyTool->FindCorrelationAxes(fMyTool->GetListJ(l, 1)  , vecJ, phi1, scenario))
	    {
	      fhTMA_JAA[l]->Fill(fMyTool->CalcdPhi0To2pi(phiA-phi1));
	    }


	  for(Int_t ip =0; ip<N1; ip++) 
	    {
	      phi1 = fMyTool->GetLocPhiJ(l,1,ip);
	      dphi = fMyTool->CalcdPhi0To2pi(phiA, phi1);
	      fhTMA_JAp[l]->Fill(dphi);
	    }


	  Int_t NB1 = fMyTool->GetSizeB1(l, 1);
	  for(Int_t ip =0; ip<NB1; ip++) 
	    {
	      phi1 = fMyTool->GetLocPhiB1(l,1,ip);
	      dphi = fMyTool->CalcdPhi0To2pi(phiA, phi1);
	      fhTMA_B1Ap[l]->Fill(dphi);
	    }

	  Int_t NB2 = fMyTool->GetSizeB2(l, 1);
	  for(Int_t ip =0; ip<NB2; ip++) 
	    {
	      phi1 = fMyTool->GetLocPhiB2(l,1,ip);
	      dphi = fMyTool->CalcdPhi0To2pi(phiA, phi1);
	      fhTMA_B2Ap[l]->Fill(dphi);
	    }

	}

    fhEvent->Fill(1);

    delete fMyTool;

  return kTRUE;
}


Double_t AliAnalysisTaskJetShape::AliAnalysisTaskJetShapeHM::CalcR(TVector3 v1, TVector3 v2)
{

  Double_t dphi = CalcdPhi(v1.Phi(), v2.Phi());
  Double_t deta = v1.Eta() - v2.Eta();
  Double_t RB = TMath::Sqrt(dphi*dphi+deta*deta);
  return RB;
}


Double_t AliAnalysisTaskJetShape::AliAnalysisTaskJetShapeHM::CalcdPhi(Double_t phi1, Double_t phi2)
{

  while(phi1<0) phi1+=TMath::TwoPi();
  while(phi1>TMath::TwoPi()) phi1-=TMath::TwoPi();

  while(phi2<0) phi2+=TMath::TwoPi();
  while(phi2>TMath::TwoPi()) phi2-=TMath::TwoPi();

  Double_t dphi = phi1- phi2;
  if(dphi>TMath::Pi())dphi = dphi - 2.*TMath::Pi();
  if(dphi<(-1.*TMath::Pi()))dphi = dphi + 2.*TMath::Pi();

  return dphi;
}







TH1F* AliAnalysisTaskJetShape::AliAnalysisTaskJetShapeHM::Hist1D(const char* name, Int_t nBins, Double_t xMin, Double_t xMax,  const char* xLabel, Int_t color, const char* yLabel)
{
// create a 1D histogram

  TH1F* res = new TH1F(Form("%s_%s",fComment.Data(), name), Form("%s_%s",fComment.Data(), name), nBins, xMin, xMax);
  if (xLabel) res->GetXaxis()->SetTitle(xLabel);
  if (yLabel) res->GetYaxis()->SetTitle(yLabel);
  res->SetLineColor(color);
  return res;
}


TH1F* AliAnalysisTaskJetShape::AliAnalysisTaskJetShapeHM::Hist1D(const char* name, Int_t nBins, Double_t *xArray,  const char* xLabel, Int_t color, const char* yLabel)
{
// create a 1D histogram

  TH1F* res = new TH1F(Form("%s_%s",fComment.Data(), name), Form("%s_%s",fComment.Data(), name), nBins, xArray);
  if (xLabel) res->GetXaxis()->SetTitle(xLabel);
  if (yLabel) res->GetYaxis()->SetTitle(yLabel);
  res->SetLineColor(color);
  return res;
}

TH2F *AliAnalysisTaskJetShape::AliAnalysisTaskJetShapeHM::Hist2D(const char* name, Int_t nBinsx, Double_t xMin, Double_t xMax, Int_t nBinsy, Double_t yMin, Double_t yMax, const char* xLabel, const char* yLabel, Int_t color)
{
// create a 2D histogram

  TH2F *res = new TH2F(Form("%s_%s",fComment.Data(), name), Form("%s_%s",fComment.Data(), name), nBinsx, xMin, xMax, nBinsy, yMin, yMax);
  if (xLabel) res->GetXaxis()->SetTitle(xLabel);
  if (yLabel) res->GetYaxis()->SetTitle(yLabel);
  //  res->SetMarkerStyle(kFullCircle);
  //  res->SetOption("E");
  res->SetLineColor(color);
  //  fOutputList->Add(res);
  return res;
}


TH2F *AliAnalysisTaskJetShape::AliAnalysisTaskJetShapeHM::Hist2D(const char* name, Int_t nBinsx, Double_t xMin, Double_t xMax, Int_t nBinsy, Double_t *yArray, const char* xLabel, const char* yLabel, Int_t color, const char* zLabel)
{
// create a 2D histogram

  TH2F *res = new TH2F(Form("%s_%s",fComment.Data(), name), Form("%s_%s",fComment.Data(), name), nBinsx, xMin, xMax, nBinsy, yArray);
  if (xLabel) res->GetXaxis()->SetTitle(xLabel);
  if (yLabel) res->GetYaxis()->SetTitle(yLabel);
  if (zLabel) res->GetZaxis()->SetTitle(zLabel);
  //  res->SetMarkerStyle(kFullCircle);
  //  res->SetOption("E");
  res->SetLineColor(color);
  //  fOutputList->Add(res);
  return res;
}

TH2F *AliAnalysisTaskJetShape::AliAnalysisTaskJetShapeHM::Hist2D(const char* name, Int_t nBinsx, Double_t *yArrax, Int_t nBinsy, Double_t yMin, Double_t yMax, const char* xLabel, const char* yLabel, Int_t color, const char* zLabel)
{
// create a 2D histogram

  TH2F *res = new TH2F(Form("%s_%s",fComment.Data(), name), Form("%s_%s",fComment.Data(), name), nBinsx, yArrax, nBinsy, yMin, yMax);
  if (xLabel) res->GetXaxis()->SetTitle(xLabel);
  if (yLabel) res->GetYaxis()->SetTitle(yLabel);
  if (zLabel) res->GetZaxis()->SetTitle(zLabel);
  //  res->SetMarkerStyle(kFullCircle);
  //  res->SetOption("E");
  res->SetLineColor(color);
  //  fOutputList->Add(res);
  return res;
}

TH3F *AliAnalysisTaskJetShape::AliAnalysisTaskJetShapeHM::Hist3D(const char* name, Int_t nBinsx, Double_t xMin, Double_t xMax, 
                                                  Int_t nBinsy, Double_t yMin, Double_t yMax, 
				                  Int_t nBinsz, Double_t zMin, Double_t zMax, const char* xLabel, const char* yLabel, const char* zLabel, Int_t color)
{
// create a 3D histogram

  TH3F *res = new TH3F(Form("%s_%s",fComment.Data(), name), Form("%s_%s",fComment.Data(), name), nBinsx, xMin, xMax, nBinsy, yMin, yMax, nBinsz, zMin, zMax);
  if (xLabel) res->GetXaxis()->SetTitle(xLabel);
  if (yLabel) res->GetYaxis()->SetTitle(yLabel);
  if (zLabel) res->GetZaxis()->SetTitle(zLabel);
  //  res->SetMarkerStyle(kFullCircle);
  //  res->SetOption("E");
  res->SetLineColor(color);
  //  fOutputList->Add(res);
  return res;
}

TH3F *AliAnalysisTaskJetShape::AliAnalysisTaskJetShapeHM::Hist3D(const char* name, Int_t nBinsx, Double_t xMin, Double_t xMax, 
                                                  Int_t nBinsy, Double_t yMin, Double_t yMax, 
				                  Int_t nBinsz, Double_t *zArr, const char* xLabel, const char* yLabel, const char* zLabel, Int_t color)
{
// create a 3D histogram

  Double_t xArr[nBinsx+1], yArr[nBinsy+1];

  for(Int_t i=0; i<=nBinsx; i++) xArr[i]=xMin+i*(xMax-xMin)/nBinsx;
  for(Int_t i=0; i<=nBinsy; i++) yArr[i]=yMin+i*(yMax-yMin)/nBinsy;

  TH3F *res = new TH3F(Form("%s_%s",fComment.Data(), name), Form("%s_%s",fComment.Data(), name), nBinsx, xArr, nBinsy, yArr, nBinsz, zArr);
  if (xLabel) res->GetXaxis()->SetTitle(xLabel);
  if (yLabel) res->GetYaxis()->SetTitle(yLabel);
  if (zLabel) res->GetZaxis()->SetTitle(zLabel);
  //  res->SetMarkerStyle(kFullCircle);
  //  res->SetOption("E");
  res->SetLineColor(color);
  //  fOutputList->Add(res);
  return res;
}


//////////////////////////////////////


/*
//________________________________________________________________________
void AnalysisJetMain::ConnectInputData(Option_t *)
{
  // Connect ESD
  // Called once

  fESD=dynamic_cast<AliESDEvent*>(InputEvent());
  //   if (!fESD) {
  //     AliError("ESD not available");
  //     fAODIn = dynamic_cast<AliAODEvent*>(InputEvent());}
  //     fAODOut = dynamic_cast<AliAODEvent*>(AODEvent());
 
  fAODIn = dynamic_cast<AliAODEvent*>(InputEvent());
  //  assume that the AOD is in the general output...
  fAODOut  = AODEvent();


}
*/

void AliAnalysisTaskJetShape::SetBranchNames(const TString &branch1, const TString &branch2)
{
   fJetBranchName[0] = branch1;
   fJetBranchName[1] = branch2;
}

//________________________________________________________________________
void AliAnalysisTaskJetShape::UserCreateOutputObjects()
{
  //printf("Open1\n");
  //       const char *nameF = OpenFile(1)->GetName();
  //printf("Open2 %s\n",nameF);

  //  fTriggerAnalysis = new AliTriggerAnalysis();
  //  if (fkMC) fTriggerAnalysis->SetAnalyzeMC(1);


  fOutputList = new TList();
  fOutputList->SetOwner();

  fhPtJL  = Hist1D("hPtJL", 100, 0, 200, "Pt_{JL}"); fOutputList->Add(fhPtJL);
  fhAreaJL = Hist1D("hAreaJL", 100, 0., 4, "AreaJL"); fOutputList->Add(fhAreaJL);



  printf("AliAnalysisTaskJetShape::UserCreateOutputObjects() \n");


  fanalJetSubStrHM->SetFilterMask(fFilterMask);
  if(fkMC)   fanalJetSubStrHM->MCprod();
  fanalJetSubStrHM->InitHistos();
  fanalJetSubStrHM->AddToList(fOutputList);


  if(fkMC)
    {

      printf("AliAnalysisTaskJetShape::UserCreateOutputObjects() MC\n");
      fanalJetSubStrMCHM->InitHistos();
      fanalJetSubStrMCHM->AddToList(fOutputList);

      for(Int_t i=0; i<3; i++)
	{
	  fhPtresVsPt[i]    = Hist2D(Form("hPtresVsPt%d",i), 100, 0, 100, 100, -0.5, 0.5, "P_{t}^{gen} GeV/c", "1-P_{t}^{rec}/P_{t}^{gen} GeV/c" );      fOutputList->Add(fhPtresVsPt[i]);
	  fhPhiresVsPhi[i]  = Hist2D(Form("hPhiresVsPhi%d",i), 600, 0, TMath::TwoPi(), 128, -0.2, 0.2, "#phi^{gen} rad", "#phi^{rec} - #phi^{gen} rad" );      fOutputList->Add(fhPhiresVsPhi[i]);
	  fhEtaresVsEta[i]  = Hist2D(Form("hEtaresVsEta%d",i), 200, -1, 1, 40, -0.2, 0.2, "#eta^{gen}", "#eta^{rec} - #eta^{gen}" );      fOutputList->Add(fhEtaresVsEta[i]);
	  fhDCAxy[i]        = Hist2D(Form("hDCAxy%d",i), 100, 0, 100, 300, -3, 3, "DCAxy of prim [cm]"); fOutputList->Add(fhDCAxy[i]);
	  fhDCAz[i]         = Hist2D(Form("hDCAz%d",i) , 100, 0, 100, 300, -3, 3, "DCAz of prim [cm]") ; fOutputList->Add(fhDCAz[i]);
	  fhTrackPtEtaPhi[i]= Hist3D(Form("hTrackPtEtaPhi%d",i), 100, 0, 100, 100, -1, 1, 100, 0, TMath::TwoPi(),"P_{t} GeV/c ","#eta","#phi rad");  fOutputList->Add(fhTrackPtEtaPhi[i]);
	}
    }
  



   printf(" JetBranchName[0]=%s\n JetBranchName[1]=%s\n", fJetBranchName[0].Data(), fJetBranchName[1].Data());



  PostData(1, fOutputList);

  /*
   OpenFile(1);

  fOutputTree = new TTree("tree","Tree with Ali2JEvent");
  fOutputTree->Branch("event",&f2JEvent);
  */

//    PostData(1, fOutputTree);
}

/*
Bool_t AliAnalysisTaskJetSpectrum2::Notify()
{

  // Fetch the aod also from the input in,
  // have todo it in notify
  
  
  fAODIn = dynamic_cast<AliAODEvent*>(InputEvent());
  //  assume that the AOD is in the general output...
  fAODOut  = AODEvent();
  
  if(fNonStdFile.Length()!=0){
    // case that we have an AOD extension we need can fetch the jets from the extended output
    AliAODHandler *aodH = dynamic_cast<AliAODHandler*>(AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler());
    fAODExtension = (aodH?aodH->GetExtension(fNonStdFile.Data()):0);    
    if(!fAODExtension){
      if(fDebug>1)Printf("AODExtension found for %s",fNonStdFile.Data());
    }
  }



  return kTRUE;
}
*/




//________________________________________________________________________
void AliAnalysisTaskJetShape::UserExec(Option_t *) 
{
  //   return;
//    if(f2JEvent)
//      delete f2JEvent;
//    f2JEvent = new Ali2JEvent();
//    return;
  
  //  return;
 if(fDebug)  Printf("\n\n\nEvent #%5d", (Int_t) fEntry);
  if(fDebug) printf("NEW EVENT 0\n");

  if(!IsGoodEvent()) return;


  //   printf("\n\n\n NEW EVENT\n");

   AliAODEvent* aodE = 0;
   if(!aodE){
     if(!fESD)aodE = fAODIn;
     else aodE = fAODOut;}

   /*
   AliAODHandler *aodH = dynamic_cast<AliAODHandler*>(AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler());
if(!aodH){
    Printf("%s:%d no output aodHandler found Jet",(char*)__FILE__,__LINE__);
    return;
  }
   */
   if(fDebug) {
     printf("NEW EVENT 1:  Number of tracks = %d\n",aodE->GetNumberOfTracks());
     //     aodE->GetList()->Print();
     //     printf("Look for %s\n",fJetBranchName[0].Data());
   }


  // centrality selection
   AliCentrality *cent = 0x0;
   Double_t centrality = 0.; 

   if(fESD) {cent = fESD->GetCentrality();
     if(cent) centrality = cent->GetCentralityPercentile("V0M");}
   else     centrality=aodE->GetHeader()->GetCentrality();


   if(!fkIsPbPb) {
     centrality=1.;
   }

//  if(fDebug) printf("NEW EVENT 2\n");

   Bool_t  fUseAOD049 = kFALSE;// kTRUE;// 
      if(fUseAOD049&&centrality>=0){
	Float_t v0=0;
	AliAODVZERO* aodV0 = aodE->GetVZEROData();
	v0+=aodV0->GetMTotV0A();
	v0+=aodV0->GetMTotV0C();
	if(centrality==0 && v0 < 19500) return ;//filtering issue
	Float_t tkl = (Float_t)(aodE->GetTracklets()->GetNumberOfTracklets());
	Float_t val= 1.30552 +  0.147931 * v0;
	Float_t tklSigma[101]={176.644, 156.401, 153.789, 153.015, 142.476, 137.951, 136.127, 129.852, 127.436, 124.86, 120.788, 115.611, 113.172, 110.496, 
            109.127, 104.421, 102.479, 99.9766, 97.5152, 94.0654, 92.4602, 89.3364, 87.1342, 83.3497, 82.6216, 81.1084, 78.0793, 76.1234, 72.9434, 72.1334, 
            68.0056, 68.2755, 66.0376, 62.9666, 62.4274, 59.65, 58.3776, 56.6361, 54.5184, 53.4224, 51.932, 50.8922, 48.2848, 47.912, 46.5717, 43.4114, 
            43.2083, 41.3065, 40.1863, 38.5255, 37.2851, 37.5396, 34.4949, 33.8366, 31.8043, 31.7412, 30.8392, 30.0274, 28.8793, 27.6398, 26.6488, 25.0183, 
            25.1489, 24.4185, 22.9107, 21.2002, 21.6977, 20.1242, 20.4963, 19.0235, 19.298, 17.4103, 16.868, 15.2939, 15.2939, 16.0295, 14.186, 14.186, 15.2173, 
            12.9504, 12.9504, 12.9504, 15.264, 12.3674, 12.3674, 12.3674, 12.3674, 12.3674, 18.3811, 13.7544, 13.7544, 13.7544, 13.7544, 13.7544, 13.7544, 13.7544, 
            13.7544, 13.7544, 13.7544, 13.7544, 13.7544};
	if ( TMath::Abs(tkl-val) > 6.*tklSigma[(Int_t)centrality] ) 
	  return; //outlier
      }
   

    if(fDebug) 
      {
	printf("centrality: %f\n", centrality);
	aodE->Print();
	aodE->GetList()->Print();
      }


    if (centrality < fCentMin || centrality > fCentMax){
    //    PostData(1, fOutputList);
    return;
    }

   //   fhNEvent->Fill(0); 

   // accepted events  
   // -- end event selection --
  
   
   // get background
   AliAODJetEventBackground* externalBackground = 0;
   AliAODJetEventBackground* externalBackgroundMC = 0;


   Float_t rho = 0;
   Float_t rho_MC=0;

   if(fkIsPbPb)
     {
       if(fAODOut&&!externalBackground&&fBackgroundBranch[0].Length()){
	 externalBackground =  (AliAODJetEventBackground*)(fAODOut->FindListObject(fBackgroundBranch[0].Data()));
	 if(!externalBackground)Printf("%s:%d Background branch not found %s",(char*)__FILE__,__LINE__,fBackgroundBranch[0].Data());
       }
       if(fAODExtension&&!externalBackground&&fBackgroundBranch[0].Length()){
	 externalBackground =  (AliAODJetEventBackground*)(fAODExtension->GetAOD()->FindListObject(fBackgroundBranch[0].Data()));
	 if(!externalBackground)Printf("%s:%d Background branch not found %s",(char*)__FILE__,__LINE__,fBackgroundBranch[0].Data());
       }
       if(fAODIn&&!externalBackground&&fBackgroundBranch[0].Length()){
	 externalBackground =  (AliAODJetEventBackground*)(fAODIn->FindListObject(fBackgroundBranch[0].Data()));
	 if(!externalBackground)Printf("%s:%d Background branch not found %s",(char*)__FILE__,__LINE__,fBackgroundBranch[0].Data());
       }
       if(fAODOut&&!externalBackgroundMC&&fBackgroundBranch[1].Length()){
	 externalBackgroundMC =  (AliAODJetEventBackground*)(fAODOut->FindListObject(fBackgroundBranch[1].Data()));
	 if(!externalBackgroundMC)Printf("%s:%d MC Background branch not found %s",(char*)__FILE__,__LINE__,fBackgroundBranch[1].Data());
       }
       if(fAODExtension&&!externalBackgroundMC&&fBackgroundBranch[1].Length()){
	 externalBackgroundMC =  (AliAODJetEventBackground*)(fAODExtension->GetAOD()->FindListObject(fBackgroundBranch[1].Data()));
	 if(!externalBackgroundMC)Printf("%s:%d MC Background branch not found %s",(char*)__FILE__,__LINE__,fBackgroundBranch[1].Data());
       }
       if(fAODIn&&!externalBackgroundMC&&fBackgroundBranch[1].Length()){
	 externalBackgroundMC =  (AliAODJetEventBackground*)(fAODIn->FindListObject(fBackgroundBranch[1].Data()));
	 if(!externalBackgroundMC)Printf("%s:%d MC Background branch not found %s",(char*)__FILE__,__LINE__,fBackgroundBranch[1].Data());
       }
       if(externalBackground)rho = externalBackground->GetBackground(0);
       if(externalBackgroundMC)rho_MC = externalBackgroundMC->GetBackground(0);
     }



  if(fkMC)
    {
      AliAODVertex* primVtx = aodE->GetPrimaryVertex();
      Double_t bfield = aodE->GetMagneticField();

      TClonesArray *arrP = dynamic_cast<TClonesArray*>(aodE->FindListObject(AliAODMCParticle::StdBranchName()));
      if(!arrP)
	{
	  printf("ERROR: no Info about particles in AOD\n");
	  return;
	}

      for(int it = 0;it < aodE->GetNumberOfTracks(); it++)
	{
	  AliAODTrack *tr = aodE->GetTrack(it);
	  if(!tr) continue;
	  if((fFilterMask>0)&&!(tr->TestFilterBit(fFilterMask))) continue;
	  if(TMath::Abs(tr->Eta())>1.) continue;
	  Int_t label = TMath::Abs(tr->GetLabel());

	  AliAODMCParticle *part = dynamic_cast<AliAODMCParticle*>(arrP->At(label));
	  if(!part)continue;
	  //	  if(!part->IsPhysicalPrimary())continue;
	  //	  if(part->Charge()==0)continue;
	  //	  vec.SetXYZ(part->Px(), part->Py(), part->Pz());



	  Double_t dca[2];
	  Double_t cov[3];
       
	  AliAODTrack tmp(*tr);
	  tmp.PropagateToDCA(primVtx, bfield, 5., dca, cov);
       
	  Int_t type = 0;
	  if(!part->IsPhysicalPrimary()) type=1;
	  if(label<0) type =2;

	  fhPtresVsPt[type]->Fill(part->Pt(), 1-tr->Pt()/part->Pt());
	  fhPhiresVsPhi[type]->Fill(part->Phi(), tr->Phi() - part->Phi());
	  fhEtaresVsEta[type]->Fill(part->Eta(), tr->Eta() - part->Eta());
	  fhDCAxy[type]->Fill(part->Pt(), dca[0]);
	  fhDCAz[type]->Fill( part->Pt(), dca[1]);
	  fhTrackPtEtaPhi[type]->Fill(tr->Pt(), tr->Eta(), tr->Phi());
	}


    }





   //   printf("rho = %f\n",rho);

    //   return;

   TClonesArray *Jets[2];
   Jets[0]=0;
   Jets[1]=0;
   if(fAODOut&&!Jets[0]){
   Jets[0] = dynamic_cast<TClonesArray*>(fAODOut->FindListObject(fJetBranchName[0].Data())); 
   Jets[1] = dynamic_cast<TClonesArray*>(fAODOut->FindListObject(fJetBranchName[1].Data()));  }
   if(fAODExtension && !Jets[0]){ 
   Jets[0] = dynamic_cast<TClonesArray*>(fAODExtension->GetAOD()->FindListObject(fJetBranchName[0].Data())); 
   Jets[1] = dynamic_cast<TClonesArray*>(fAODExtension->GetAOD()->FindListObject(fJetBranchName[1].Data()));  }
   if(fAODIn&&!Jets[0]){
   Jets[0] = dynamic_cast<TClonesArray*>(fAODIn->FindListObject(fJetBranchName[0].Data())); 
   Jets[1] = dynamic_cast<TClonesArray*>(fAODIn->FindListObject(fJetBranchName[1].Data())); }


   if(!Jets[0]) {
     Printf("no friend!!!\n");
     return;
   }
  Int_t nJ = Jets[0]->GetEntries();
  Float_t ptmax = 0.;
  Int_t   imax  = -1;

  if(fDebug) {
    printf("NEW EVENT 3:  Number of Rec. jets  %d\n",nJ);
  }

//
// Find highest pT jet with pt > 20 GeV
//
//  fPtJcorrMin=0;

  Float_t etaJmax = 0.4;

     for (Int_t i = 0; i < nJ; i++) {
       AliAODJet* jet = dynamic_cast<AliAODJet*> (Jets[0]->At(i));
       if (!jet) continue;
       //              jet->Print("0");
       Float_t ptJ  = jet->Pt();
       Float_t etaJ = TMath::Abs(jet->Eta());

       Float_t area = jet->EffectiveAreaCharged();
       Float_t ptJcorr=ptJ-rho*area;
       
       if ((ptJcorr > fPtJcorrMin) && (ptJcorr  > ptmax) && etaJ < etaJmax) {
	 ptmax = ptJcorr;
	 imax = i;
       }
       
     }


  if(fDebug) {
    printf("NEW EVENT 4:\n");
  }

     TVector3 vecJ;

     AliAODJet* jetL = 0;

     Float_t areaJL, ptJLcorr;

     if (imax != -1)  {

       jetL = dynamic_cast<AliAODJet*> (Jets[0]->At(imax));
       if(jetL){

	 areaJL   = jetL->EffectiveAreaCharged();
	 ptJLcorr  = jetL->Pt()-rho*areaJL;
	 
	 fhPtJL->Fill(ptJLcorr);
	 fhAreaJL->Fill(areaJL);
	 vecJ.SetXYZ(jetL->Px(), jetL->Py(), jetL->Pz());
	 fanalJetSubStrHM->AddEvent(aodE, jetL, rho*areaJL);

	 if(fDebug) {
	   Printf("Leading Jet");
	   jetL->Print("0");
	 }

       }
     }


     if(!fkMC)
       {
	 PostData(1, fOutputList);
	 if(fDebug)  Printf("End of event #%5d", (Int_t) fEntry);
	 return;
       }




     nJ = Jets[1]->GetEntries();
  if(fDebug) {
    printf("NEW EVENT 5:  Number of Rec. jets  %d\n",nJ);
  }

     ptmax = 0;
     imax = -1;
     Float_t areaJL_MC=0;

     for (Int_t i = 0; i < nJ; i++) {
       AliAODJet* jet = dynamic_cast<AliAODJet*> (Jets[1]->At(i));
       if (!jet) continue;
       Float_t ptJ1  = jet->Pt();
       Float_t etaJ1 = TMath::Abs(jet->Eta());

       areaJL_MC = jet->EffectiveAreaCharged();
       Float_t ptJcorr=ptJ1-rho*areaJL_MC;

       //    jet->Print();

       if ((ptJcorr > fPtJcorrMin) && (ptJcorr  > ptmax) && etaJ1 < etaJmax) {
	 ptmax = ptJcorr;
	 imax = i;
       }
       
     }

  if(fDebug) {
    printf("NEW EVENT 6:\n");
  }

 
     AliAODJet* jetMCL=0;

    if (imax != -1)  {
      jetMCL = dynamic_cast<AliAODJet*> (Jets[1]->At(imax));
      if(jetMCL){
	fanalJetSubStrMCHM->AddEvent(aodE, jetMCL, rho_MC*areaJL_MC);
      }
    }





     PostData(1, fOutputList);

    if(fDebug)  Printf("End of event #%5d", (Int_t) fEntry);

     return;




}      













Double_t AliAnalysisTaskJetShape::CalcdPhi(Double_t phi1, Double_t phi2)
{
  while(phi1<0) phi1+=TMath::TwoPi();
  while(phi1>TMath::TwoPi()) phi1-=TMath::TwoPi();

  while(phi2<0) phi2+=TMath::TwoPi();
  while(phi2>TMath::TwoPi()) phi2-=TMath::TwoPi();

  Double_t dphi = phi1- phi2;
  if(dphi>TMath::Pi())dphi = dphi - 2.*TMath::Pi();
  if(dphi<(-1.*TMath::Pi()))dphi = dphi + 2.*TMath::Pi();

  //  Double_t dphi = phi1- phi2;
  //  while(dphi<0) dphi+=TMath::TwoPi();
  //  while(dphi>TMath::TwoPi()) dphi-=TMath::TwoPi();
  return dphi;
}

Double_t AliAnalysisTaskJetShape::CalcdPhi0To2pi(Double_t phi1, Double_t phi2)
{

  Double_t dphi = CalcdPhi(phi1, phi2);
  while(dphi<0) dphi+=TMath::TwoPi();
  while(dphi>TMath::TwoPi()) dphi-=TMath::TwoPi();
  return dphi;
}


Bool_t AliAnalysisTaskJetShape::IsGoodEvent()
{
  
  fESD=dynamic_cast<AliESDEvent*>(InputEvent());
  //   if (!fESD) {
  //     AliError("ESD not available");
  //     fAODIn = dynamic_cast<AliAODEvent*>(InputEvent());}
  //     fAODOut = dynamic_cast<AliAODEvent*>(AODEvent());
 
  fAODIn = dynamic_cast<AliAODEvent*>(InputEvent());
  //  assume that the AOD is in the general output...
  fAODOut  = AODEvent();
  
    
   static AliAODEvent* aod = 0;
   // take all other information from the aod we take the tracks from
   if(!aod){
     if(!fESD)aod = fAODIn;
     else aod = fAODOut;}
     

 
   if(fNonStdFile.Length()!=0){
    // case that we have an AOD extension we need can fetch the jets from the extended output
     AliAODHandler *aodH = dynamic_cast<AliAODHandler*>(AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler());
     fAODExtension = (aodH?aodH->GetExtension(fNonStdFile.Data()):0);
   }

   if(fDebug){
    if(fAODIn) printf("fAODIn\n");
    if(fAODOut) printf("fAODOut\n");
    if(fAODExtension) printf("fAODExtension\n");
   }


   /*
   AliAODHandler *aodH = dynamic_cast<AliAODHandler*>(AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler());
if(!aodH){
    Printf("%s:%d no output aodHandler found Jet",(char*)__FILE__,__LINE__);
      return kFALSE;
  }
   */

   // -- event selection --

   // physics selection
   AliInputEventHandler* inputHandler = (AliInputEventHandler*)
   ((AliAnalysisManager::GetAnalysisManager())->GetInputEventHandler());
   //     cout<<inputHandler->IsEventSelected()<<" "<<fOfflineTrgMask<<endl;
   if(!(inputHandler->IsEventSelected() & fOfflineTrgMask)){
      if(fDebug) Printf(" Trigger Selection: event REJECTED ... ");
      return kFALSE;
   }

   // vertex selection
   //   if(!aod){
   //     if(fDebug) Printf("%s:%d No AOD",(char*)__FILE__,__LINE__);
   //     fhNEvent->Fill(3);
   //      PostData(1, fOutputList);
   //     return kFALSE;
   //   }

   AliAODVertex* primVtx = aod->GetPrimaryVertex();

   if(!primVtx){
     if(fDebug) Printf("%s:%d No primVtx",(char*)__FILE__,__LINE__);
     return kFALSE;
   }

   Int_t nTracksPrim = primVtx->GetNContributors();
   if ((nTracksPrim < fVtxMinContrib) ||
         (primVtx->GetZ() < fVtxZMin) ||
         (primVtx->GetZ() > fVtxZMax) ){
      if(fDebug) Printf("%s:%d primary vertex z = %f: event REJECTED...",(char*)__FILE__,__LINE__,primVtx->GetZ());
      return  kFALSE;
   }

   /*
   // event class selection (from jet helper task)
   Int_t eventClass = AliAnalysisHelperJetTasks::EventClass();
   if(fDebug) Printf("Event class %d", eventClass);
   if (eventClass < fEvtClassMin || eventClass > fEvtClassMax){
      return  kFALSE;
   }
   */
    return kTRUE;
}






//________________________________________________________________________
void AliAnalysisTaskJetShape::Terminate(Option_t *) 
{
  // Draw result to the screen
  // Called once at the end of the query
 //    fModuleMap = dynamic_cast<TH1F*> (fOutput->FindObject("fModuleMap"));

  printf("MyTaskTestTrigger::Terminate()\n\n\n");

  
  //  fOutputList = dynamic_cast<TList*> (GetOutputData(1));
  //  fOutputTree = dynamic_cast<TTree*> (GetOutputData(1));

}



/*
Bool_t AnalysisJetMain::GetVertex(const AliESDEvent* esd,  Double_t Vxyz[3], Int_t type)
{
  // type
  //     0 -> SPD vtx
  //     1 -> ESD vtx
  const AliESDVertex* vtx = 0;

  if(type==0) {
    vtx = esd->GetPrimaryVertexSPD();
    if(!vtx) return kFALSE;
    if(vtx->GetNContributors() < 1) return kFALSE;
    TString vtxTyp = vtx->GetTitle();
    if (vtxTyp.Contains("vertexer: Z")) {
      if (vtx->GetDispersion()>0.04) return kFALSE;
      if (vtx->GetZRes()>0.25) return kFALSE;
    }
 
  }
  if(type==1) {
    vtx = esd->GetPrimaryVertexTracks();
    if(!vtx) return kFALSE;
    if(vtx->GetNContributors() < 1) return kFALSE;
  }

 
   Vxyz[0] = vtx->GetXv();
   Vxyz[1] = vtx->GetYv();
   Vxyz[2] = vtx->GetZv();
 
   return kTRUE;
}
*/









TH1F* AliAnalysisTaskJetShape::Hist1D(const char* name, Int_t nBins, Double_t xMin, Double_t xMax,  const char* xLabel, Int_t color)
{
// create a 1D histogram

  TH1F *h = new TH1F(name, name, nBins, xMin, xMax);
  if (xLabel) h->GetXaxis()->SetTitle(xLabel);
  //  res->SetMarkerStyle(kFullCircle);
  //  res->SetOption("E");
  h->SetLineColor(color);
  //  fOutputList->Add(h);
  return h;
}


TH2F *AliAnalysisTaskJetShape::Hist2D(const char* name, Int_t nBinsx, Double_t xMin, Double_t xMax, Int_t nBinsy, Double_t yMin, Double_t yMax, const char* xLabel, const char* yLabel, Int_t color)
{
// create a 2D histogram

  TH2F *res = new TH2F(name, name, nBinsx, xMin, xMax, nBinsy, yMin, yMax);
  if (xLabel) res->GetXaxis()->SetTitle(xLabel);
  if (xLabel) res->GetYaxis()->SetTitle(yLabel);
  //  res->SetMarkerStyle(kFullCircle);
  //  res->SetOption("E");
  res->SetLineColor(color);
  //  fOutputList->Add(res);
  return res;
}


TH3F *AliAnalysisTaskJetShape::Hist3D(const char* name, Int_t nBinsx, Double_t xMin, Double_t xMax, 
                                                  Int_t nBinsy, Double_t yMin, Double_t yMax, 
				                  Int_t nBinsz, Double_t zMin, Double_t zMax, const char* xLabel, const char* yLabel, const char* zLabel, Int_t color)
{
// create a 3D histogram

  TH3F *res = new TH3F(name, name, nBinsx, xMin, xMax, nBinsy, yMin, yMax, nBinsz, zMin, zMax);
  if (xLabel) res->GetXaxis()->SetTitle(xLabel);
  if (yLabel) res->GetYaxis()->SetTitle(yLabel);
  if (zLabel) res->GetZaxis()->SetTitle(zLabel);
  //  res->SetMarkerStyle(kFullCircle);
  //  res->SetOption("E");
  res->SetLineColor(color);
  //  fOutputList->Add(res);
  return res;
}




//#endif








