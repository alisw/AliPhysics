// ---- ROOT system ---
#include "TDirectory.h"
#include "TBranch.h"
#include "TTree.h"
#include "TROOT.h"
#include "TFolder.h"

// --- AliRoot header files ---
#include "AliConfig.h"
#include "AliPHOSClusterizer.h"
#include "AliPHOSEvalRecPoint.h"
#include "AliRun.h"
#include "AliPHOSLoader.h"
#include "AliPHOSRecCpvManager.h"
#include "AliPHOSRecEmcManager.h"
#include "AliPHOSDigitizer.h"

// --- Standard library ---


ClassImp(AliPHOSEvalRecPoint)

  AliPHOSEvalRecPoint::AliPHOSEvalRecPoint(): fEventFolderName(AliConfig::fgkDefaultEventFolderName)
{
  fParent=-333;
  fChi2Dof=-111;
  fIsCpv = kTRUE;
  fIsEmc = kFALSE;
}

AliPHOSEvalRecPoint::AliPHOSEvalRecPoint(Bool_t cpv, AliPHOSEvalRecPoint* parent) 
  : fEventFolderName(AliConfig::fgkDefaultEventFolderName)
{

  fParent=-333;
  fChi2Dof=-111;

//    fParent=parent;
  TObjArray* wPool = (TObjArray*)GetWorkingPool();
  if(!wPool) {
    Error("AliPHOSEvalRecPoint", "Couldn't find working pool. Exit.") ; 
    exit(1);
  }

  fParent = wPool->IndexOf((TObject*)parent);
  fChi2Dof = parent->Chi2Dof();

  if(cpv) {
    fIsEmc = kFALSE;
    fIsCpv = kTRUE;
  }
  else {
    fIsEmc = kTRUE;
    fIsCpv = kFALSE;
  }

  // Add itself to working pool
  AddToWorkingPool((TObject*)this);
  
}

AliPHOSEvalRecPoint::AliPHOSEvalRecPoint(Int_t i, Bool_t cpv) : fEventFolderName(AliConfig::fgkDefaultEventFolderName)
{

  fChi2Dof=-111;
  fParent=-333;

  AliPHOSEmcRecPoint* rp=0;

  AliPHOSLoader* fLoader = AliPHOSLoader::GetPHOSLoader(fEventFolderName);

  if(cpv) {
    rp = (AliPHOSCpvRecPoint*)fLoader->CpvRecPoints()->At(i);
    fIsEmc = kFALSE;
    fIsCpv = kTRUE;
  }
  else {
    rp = (AliPHOSEmcRecPoint*)fLoader->EmcRecPoints()->At(i);
    fIsEmc = kTRUE;
    fIsCpv = kFALSE;
  }

  Int_t* Digits = rp->GetDigitsList();
  Float_t* Energies = rp->GetEnergiesList();
  Int_t nDigits = rp->GetMultiplicity();
  
  for(Int_t iDigit=0; iDigit<nDigits; iDigit++) {
    AliPHOSDigit* digit = (AliPHOSDigit*)fLoader->Digits()->At( Digits[iDigit] );
    Float_t eDigit = Energies[iDigit];
    this->AddDigit(*digit,eDigit);
  }

  TVector3 locpos;
  rp->GetLocalPosition(locpos);
  fLocPos = locpos; 

  // Add itself to working pool
  AddToWorkingPool((TObject*)this);
   
}

AliPHOSClusterizer* AliPHOSEvalRecPoint::GetClusterizer()
{
  TFolder* aliceF = (TFolder*)gROOT->FindObjectAny("YSAlice");
  TFolder* wPoolF = (TFolder*)aliceF->FindObject("WhiteBoard/RecPoints/PHOS/SmP");
  AliPHOSClusterizer* clu = (AliPHOSClusterizer*)wPoolF->FindObject("PHOS:clu-v1");
  if(!clu) { 
    Error("GetClusterizer", "Couldn't find Clusterizer. Exit.") ; 
    exit(1); 
  }

  return clu;
}

Bool_t AliPHOSEvalRecPoint::TooClose(AliPHOSRecPoint* pt) const 
{
  TVector3 her_pos;
  TVector3 my_pos;
  pt->GetLocalPosition(her_pos);
  this->GetLocalPosition(my_pos);
  Float_t dx = her_pos.X() - my_pos.X();
  Float_t dz = her_pos.Z() - my_pos.Z();
  Float_t dr = TMath::Sqrt(dx*dx + dz*dz);

  if(dr<GetReconstructionManager()->MergeGammasMinDistanceCut())
    return kTRUE;
  else
    return kFALSE;

}

Bool_t AliPHOSEvalRecPoint::NeedToSplit() const 
{
  return kFALSE;
}

void AliPHOSEvalRecPoint::DeleteParent()
{
  fParent=-333;
}

void AliPHOSEvalRecPoint::UpdateWorkingPool()
{

  Int_t i; //loop variable
  
  for(i=0; i<this->InWorkingPool(); i++) {
     AliPHOSEvalRecPoint* parent = (AliPHOSEvalRecPoint*)GetFromWorkingPool(i);
     TObjArray children;
     Int_t nChild = parent->HasChild(children);
     for(Int_t iChild=0; iChild<nChild; iChild++)
       {
        ((AliPHOSEvalRecPoint*)children.At(iChild))->DeleteParent();
       }

     if(nChild) {
       RemoveFromWorkingPool(parent);
       delete parent;
     }

  }

  for(i=0; i<this->InWorkingPool(); i++) {
    AliPHOSEvalRecPoint* weak = (AliPHOSEvalRecPoint*)GetFromWorkingPool(i);
    if (weak->KillWeakPoint()) delete weak;
  }

  for(i=0; i<this->InWorkingPool(); i++) {
    AliPHOSEvalRecPoint* close = (AliPHOSEvalRecPoint*)GetFromWorkingPool(i);
    close->MergeClosePoint();
  }

  for(i=0; i<this->InWorkingPool(); i++)
    ((AliPHOSEvalRecPoint*)AliPHOSEvalRecPoint::GetFromWorkingPool(i))->SetIndexInList(i);

}

Int_t AliPHOSEvalRecPoint::HasChild(TObjArray& children)
{
  for(Int_t iChild=0; iChild<InWorkingPool(); iChild++)
    {
      AliPHOSEvalRecPoint* child = (AliPHOSEvalRecPoint*)GetFromWorkingPool(iChild);
      TObject* parent = (TObject*)child->Parent();
      TObject* me = (TObject*)this;
      if(parent==me) children.Add(child);
    }

  return children.GetEntries();
}

void AliPHOSEvalRecPoint::Init()
{

  AliPHOSClusterizer* clusterizer = GetClusterizer();
  if(!clusterizer) {
    Error("Init", "Cannot get clusterizer. Exit.") ;
    exit(1);
  }

  TClonesArray* digits = AliPHOSLoader::GetPHOSLoader(fEventFolderName)->Digits();
  Float_t LogWeight=0;  

  if(this->IsEmc()) {
    LogWeight = clusterizer->GetEmcLogWeight();
  }
  else {
    LogWeight = clusterizer->GetCpvLogWeight();
  }
  
  EvalLocalPosition(LogWeight,digits); // evaluate initial position
}


void AliPHOSEvalRecPoint::MakeJob()
{
  // Reconstruction algoritm implementation.

  AliPHOSRecManager* recMng = GetReconstructionManager();

  Init();

  UnfoldLocalMaxima();
  
  TObjArray children;
  Int_t nChild = HasChild(children);
  
  if(!nChild) {
    EvaluatePosition();
    if(Chi2Dof()>recMng->OneGamChisqCut()) SplitMergedShowers();
  }

  for(Int_t iChild=0; iChild<nChild; iChild++) {
    AliPHOSEvalRecPoint* child = (AliPHOSEvalRecPoint*)children.At(iChild);
    child->EvaluatePosition();

    if(child->Chi2Dof()>recMng->OneGamChisqCut())
      child->SplitMergedShowers();
  }  

} 

void AliPHOSEvalRecPoint::InitTwoGam(Float_t* gamma1, Float_t* gamma2)
{
  //Compute start values for two gamma fit algorithm.
  // gamma1[0], gamma1[1], gamma1[2] are Energy,x,z of the reconstructed "gammas".


  TVector3 lpos; // start point choosing.
  GetLocalPosition(lpos);
  Float_t xx = lpos.Z();
  Float_t yy = lpos.X();
  Float_t E = GetEnergy();

  Info("InitTwoGam", "(x,z,E)[old] = (%f, %f, %f)", yy, xx, E) ;

//    xx = XY(xx/E);
//    yy = XY(yy/E);


  Float_t eDigit ;
  AliPHOSDigit * digit ;
  Int_t nDigits = GetMultiplicity();  
  Int_t * Digits = GetDigitsList();
  Float_t * Energies = GetEnergiesList();

  Float_t ix ;
  Float_t iy ;
  Int_t relid[4] ; 

  Float_t exx = 0;
  Float_t eyy = 0;
  Float_t exy = 0;
  Float_t sqr;
  Float_t cos2fi = 1.;

  AliPHOSLoader* fLoader = AliPHOSLoader::GetPHOSLoader(fEventFolderName);
  const AliPHOSGeometry* fGeom = fLoader->PHOSGeometry();

  Int_t iDigit; //loop variable

  for(iDigit = 0 ; iDigit < nDigits ; iDigit ++)
    {
      digit = (AliPHOSDigit*)fLoader->Digits()->At( Digits[iDigit] ); 
      eDigit = Energies[iDigit];
      fGeom->AbsToRelNumbering(digit->GetId(), relid) ;
      fGeom->RelPosInModule(relid, iy, ix);
    
      Float_t dx =  ix - xx;
      Float_t dy =  iy - yy;
      exx += eDigit*dx*dx;
      eyy += eDigit*dy*dy;
      exy += eDigit*dx*dy;
      
    }

  sqr = TMath::Sqrt(4.*exy*exy + (exx-eyy)*(exx-eyy));
  Float_t euu = (exx+eyy+sqr)/2.;
  if(sqr>1.e-10) cos2fi = (exx - eyy)/sqr;
  Float_t cosfi = TMath::Sqrt((1.+cos2fi)/2.);
  Float_t sinfi = TMath::Sqrt((1.-cos2fi)/2.);
  if(exy<0) sinfi = -sinfi;

  Float_t eu3 = 0;
  for(iDigit = 0 ; iDigit < nDigits ; iDigit ++)
    {
      digit = (AliPHOSDigit*)fLoader->Digits()->At( Digits[iDigit] ); 
      eDigit = Energies[iDigit];
      fGeom->AbsToRelNumbering(digit->GetId(), relid) ;
      fGeom->RelPosInModule(relid, iy, ix);
    
      Float_t dx =  ix - xx;
      Float_t dy =  iy - yy;
      Float_t du = dx*cosfi + dy*sinfi;
      eu3 += eDigit*du*du*du;
    }
  
  Float_t c = E*eu3*eu3/(euu*euu*euu)/2.;
  Float_t sign = 1.;
  if(eu3<0) sign = -1.;
  Float_t r = 1.+ c + sign*TMath::Sqrt(2.*c + c*c);
  Float_t u = 0;
  if(TMath::Abs(r-1.)>0.1) u = eu3/euu*(r+1.)/(r-1.);
  if(TMath::Abs(r-1.)<0.1) u = TMath::Sqrt(sqr/E/r)*(1.+r);
  
  Float_t e2c = E/(1.+r);
  Float_t e1c = E-e2c;
  Float_t u1 = -u/(1.+r);
  Float_t u2 = u+u1;
  
  Float_t x1c = xx+u1*cosfi;
  Float_t y1c = yy+u1*sinfi;
  Float_t x2c = xx+u2*cosfi;
  Float_t y2c = yy+u2*sinfi;

//    printf("e1c -> %f\n",e1c);
//    printf("x1c -> %f\n",x1c);
//    printf("y1c -> %f\n",y1c);
//    printf("e2c -> %f\n",e2c);
//    printf("x2c -> %f\n",x2c);
//    printf("y2c -> %f\n",y2c);

  gamma1[0] = e1c;
  gamma1[1] = y1c;
  gamma1[2] = x1c;

  gamma2[0] = e2c;
  gamma2[1] = y2c;
  gamma2[2] = x2c;
  
}

void AliPHOSEvalRecPoint::TwoGam(Float_t* gamma1, Float_t* gamma2)
{
  //Fitting algorithm for the two very closed gammas 
  //that merged into the cluster with one maximum.
  // Starting points gamma1 and gamma2 must be provided before ( see InitTwoGam)
  //Set chisquare of the fit.


  Float_t st0 = GetReconstructionManager()->TwoGamInitialStep();
  Float_t chmin = GetReconstructionManager()->TwoGamChisqMin();
  Float_t emin = GetReconstructionManager()->TwoGamEmin();
  Float_t stpmin = GetReconstructionManager()->TwoGamStepMin();
  Int_t Niter = GetReconstructionManager()->TwoGamNumOfIterations(); // Number of iterations.

  Float_t chisq = 100.; //Initial chisquare.

  Int_t nadc = GetMultiplicity();
  if(nadc<3)
      fChi2Dof= -111.;
  
  Int_t dof = nadc - 5;
  if(dof<1) dof=1;
  Float_t chstop = chmin*dof;

  Float_t ch = 1.e+20;
  Float_t st = st0;
  Float_t grx1 = 0.;
  Float_t gry1 = 0.;
  Float_t grx2 = 0.;
  Float_t gry2 = 0.;
  Float_t gre = 0.;
  Float_t gr = 1.;
  Float_t ee1=0,xx1=0,yy1=0,ee2=0,xx2=0,yy2=0,cosi;

  Float_t e1c = gamma1[0];
  Float_t y1c = gamma1[1];
  Float_t x1c = gamma1[2];

  Float_t e2c = gamma2[0];
  Float_t y2c = gamma2[1];
  Float_t x2c = gamma2[2];

  Float_t E = GetEnergy();

  Float_t eDigit ;
  AliPHOSDigit * digit ;
  Int_t nDigits = GetMultiplicity();  
  Int_t * Digits = GetDigitsList();
  Float_t * Energies = GetEnergiesList();

  Float_t ix ;
  Float_t iy ;
  Int_t relid[4] ; 

  AliPHOSLoader* fLoader = AliPHOSLoader::GetPHOSLoader(fEventFolderName);
  const AliPHOSGeometry* fGeom = fLoader->PHOSGeometry();

  for(Int_t Iter=0; Iter<Niter; Iter++)
    {
      Float_t chisqc = 0.;
      Float_t grx1c = 0.;
      Float_t gry1c = 0.;
      Float_t grx2c = 0.;
      Float_t gry2c = 0.;
      Float_t grec = 0.;

      for(Int_t iDigit = 0 ; iDigit < nDigits ; iDigit ++)
	{
	  digit = (AliPHOSDigit*)fLoader->Digits()->At( Digits[iDigit] ); 
	  eDigit = Energies[iDigit];
	  fGeom->AbsToRelNumbering(digit->GetId(), relid) ;
	  fGeom->RelPosInModule(relid, iy, ix);
	  
	  Float_t a1,gx1,gy1;
	  Float_t a2,gx2,gy2;
	  
	  Float_t dx1 =  x1c - ix;
	  Float_t dy1 =  y1c - iy;
//  	  Info("TwoGam", "Mult %d dx1 %f dy1 %f", nDigits, dx1, dy1) ;
//  	  AG(e1c,dx1,dy1,a1,gx1,gy1);
	  GetReconstructionManager()->AG(e1c,dx1,dy1,a1,gx1,gy1);

	  Float_t dx2 =  x2c - ix;
	  Float_t dy2 =  y2c - iy;
//  	  Info("TwoGam", "      dx2 %f dy2 %f", dx2, dy2) ;
//  	  AG(e2c,dx2,dy2,a2,gx2,gy2);
	  GetReconstructionManager()->AG(e2c,dx2,dy2,a2,gx2,gy2);
      
	  Float_t A = a1+a2;
//  	  Float_t D = Const*A*(1. - A/E);
//  	  if(D<0) D=0;

//  //    	  D = 9.; // ????
	  
//  	  Float_t da = A - eDigit;
//  	  chisqc += da*da/D;
//  	  Float_t dd = da/D;
//  	  dd = dd*(2.-dd*Const*(1.-2.*A/E));

	  Float_t dd;
	  chisqc += GetReconstructionManager()->TwoGamChi2(A,eDigit,E,dd);
	  
	  grx1c += gx1*dd;
	  gry1c += gy1*dd;
	  grx2c += gx2*dd;
	  gry2c += gy2*dd;
	  grec += (a1/e1c - a2/e2c)*E*dd;

	}

      Float_t grc = TMath::Sqrt(grx1c*grx1c + gry1c*gry1c + grx2c*grx2c + gry2c*gry2c);
      if(grc<1.e-10) grc=1.e-10;

      Float_t sc = 1. + chisqc/ch;
      st = st/sc; 

      if(chisqc>ch) 
	goto loop20;
      cosi = (grx1*grx1c + gry1*gry1c + grx2*grx2c + gry2*gry2c + gre*grec)/gr/grc;
      st = st*sc/(1.4 - cosi);

      ee1 = e1c;
      xx1 = x1c;
      yy1 = y1c;
      ee2 = e2c;
      xx2 = x2c;
      yy2 = y2c;

      ch = chisqc;

      if(ch<chstop)
	goto loop101;

      grx1 = grx1c;
      gry1 = gry1c;
      grx2 = grx2c;
      gry2 = gry2c;
      gre = grec;
      gr = grc;
 
    loop20: ;
      Float_t step = st*gr;

      Info("TwoGam", "Iteration %d dof %d chisq/dof %f chstop/dof %f step %d stpmin %d",
	   Iter, dof, ch/dof, chstop/dof, step, stpmin) ;

      
      if(step<stpmin) 
	goto loop101;

      Float_t de = st*gre*E;
      e1c = ee1 - de;
      e2c = ee2 + de;
      
      if( (e1c>emin) && (e2c>emin) )
	goto loop25;
      st = st/2.;
      goto loop20;
      
    loop25: ;
      x1c = xx1 - st*grx1;
      y1c = yy1 - st*gry1;
      x2c = xx2 - st*grx2;
      y2c = yy2 - st*gry2;


    }

 loop101: 

//    if(ch>chisq*(nadc-2)-delch)
//      return ch/dof;  
  
  chisq = ch/dof;
  gamma1[0] = ee1;
  gamma1[1] = yy1;
  gamma1[2] = xx1;

  gamma2[0] = ee2;
  gamma2[1] = yy2;
  gamma2[2] = xx2;

  Float_t x1_new = yy1;
  Float_t z1_new = xx1;
  Float_t e1_new = ee1;
  Float_t x2_new = yy2;
  Float_t z2_new = xx2;
  Float_t e2_new = ee2;

  TString message ; 
  message  = "     (x,z,E)[1 fit] = (%f, %f, %f)\n" ;  
  message  = "     (x,z,E)[2 fit] = (%f, %f, %f)\n" ;  
  Info("TwoGam", message.Data(), 
       x1_new, z1_new, e1_new, 
       x2_new, z2_new, e2_new) ;

  fChi2Dof = chisq;

}

void AliPHOSEvalRecPoint::UnfoldTwoMergedPoints(Float_t* gamma1, Float_t* gamma2)
{
  //Unfold TWO merged rec. points in the case when cluster has only one maximum, 
  //but it's fitting to the one gamma shower is too bad.
  // Use TwoGam() to estimate the positions and energies of merged points.
 
  
  Int_t Nmax = 2;
  Float_t* gamma;

  Int_t* Digits = GetDigitsList();
  Int_t nDigits = GetMultiplicity();
  Float_t* Energies = GetEnergiesList();
  Float_t* eFit = new Float_t[nDigits];

  AliPHOSLoader* fLoader = AliPHOSLoader::GetPHOSLoader(fEventFolderName);
  const AliPHOSGeometry* fGeom = fLoader->PHOSGeometry();

  for(Int_t iDigit=0; iDigit<nDigits; iDigit++)
    {
      AliPHOSDigit* digit = (AliPHOSDigit*)fLoader->Digits()->At( Digits[iDigit] ); 
      Int_t relid[4] ;
      fGeom->AbsToRelNumbering(digit->GetId(), relid) ;
      Float_t x,z;     
      fGeom->RelPosInModule(relid, x, z);

      Float_t gain = 0.;
      for(Int_t iMax=0; iMax<Nmax; iMax++)
	{
	  if(iMax==0)
	    gamma = gamma1;
	  else
	    gamma = gamma2;

	  Float_t eMax = gamma[0];
	  Float_t xMax = gamma[1];
	  Float_t zMax = gamma[2];

	  Float_t dx = xMax - x;
	  Float_t dz = zMax - z;
	  Float_t Amp,gx,gy;
	  GetReconstructionManager()->AG(eMax,dz,dx,Amp,gx,gy); 
	  gain += Amp;
	}
 
      eFit[iDigit] = gain;
    }


  for( Int_t iMax=0; iMax<Nmax; iMax++)
    {      
      if(iMax==0)
	gamma = gamma1;
      else
	gamma = gamma2;

      Float_t eMax = gamma[0];
      Float_t xMax = gamma[1];
      Float_t zMax = gamma[2];

      AliPHOSEvalRecPoint* newRP = new AliPHOSEvalRecPoint(IsCPV(),this);
      TVector3 newpos(xMax,0,zMax);
      newRP->SetLocalPosition(newpos);

      for( Int_t iDigit = 0 ; iDigit < nDigits ; iDigit ++)
	{
	  AliPHOSDigit* digit = (AliPHOSDigit*)fLoader->Digits()->At( Digits[iDigit] ); 
	  Float_t eDigit = Energies[iDigit];
	  Int_t relid[4] ;
	  fGeom->AbsToRelNumbering(digit->GetId(), relid) ;
	  Float_t ix,iz;
	  fGeom->RelPosInModule(relid, ix, iz);
	  
	  Float_t dx = xMax - ix;
	  Float_t dz = zMax - iz;
	  Float_t single_shower_gain,gxMax,gyMax;
	  GetReconstructionManager()->AG(eMax,dz,dx,single_shower_gain,gxMax,gyMax);
	  Float_t total_gain = eFit[iDigit];
	  Float_t ratio = single_shower_gain/total_gain; 
	  eDigit = eDigit*ratio;
	  newRP->AddDigit(*digit,eDigit);
	}
      Info("UnfoldTwoMergedPoints", "======= Split: daughter rec point %d =================", iMax) ;
      newRP->Print("");

    }

  delete[] eFit;
  

}


void AliPHOSEvalRecPoint::EvaluatePosition() 
{
  // One gamma fit algorithm.
  // Set chisq/dof of the fit.
  // gamma1[0], gamma1[1], gamma1[2] are Energy,x,z of the reconstructed gamma, respectively.

  Int_t nDigits = GetMultiplicity();
  if(nDigits<2)
    return;

  Int_t Niter = GetReconstructionManager()->OneGamNumOfIterations(); // number of iterations
  Float_t St0 = GetReconstructionManager()->OneGamInitialStep(); // initial step
//    const Float_t Stpmin = 0.005;
  Float_t Stpmin = GetReconstructionManager()->OneGamStepMin();
  Float_t chmin =  GetReconstructionManager()->OneGamChisqMin();
 
  TVector3 locpos;
  AliPHOSDigit* digit;
  Float_t eDigit;
  Int_t relid[4] ; 
  Float_t ix, iy;

  GetLocalPosition(locpos);
  Float_t E = GetEnergy();
  Float_t xc = locpos.Z();
  Float_t yc = locpos.X();
  Float_t dx,dy,gx,gy,grxc,gryc;
  Float_t st = St0;
  Float_t chisq = 1.e+20;
  Float_t gr = 1.;
  Float_t grx = 0.;
  Float_t gry = 0.;
  Int_t dof = GetMultiplicity() - 2;
  if(dof<1)
    dof = 1;
  Float_t chstop = chmin*dof;
  Float_t cosi,x1=0,y1=0;
  Float_t chisqc;

  AliPHOSLoader* fLoader = AliPHOSLoader::GetPHOSLoader(fEventFolderName);
  const AliPHOSGeometry* fGeom = fLoader->PHOSGeometry();

  for(Int_t Iter=0; Iter<Niter; Iter++)
    {
 
      chisqc = 0.;
      grxc = 0.;
      gryc = 0.;
      for(Int_t iDigit = 0 ; iDigit < nDigits ; iDigit ++)
	{

	  Float_t* Energies = GetEnergiesList();
	  Int_t* Digits = GetDigitsList();
	  digit = (AliPHOSDigit*)fLoader->Digits()->At( Digits[iDigit] );
	  eDigit = Energies[iDigit];
	  fGeom->AbsToRelNumbering(digit->GetId(), relid) ;
	  fGeom->RelPosInModule(relid, iy, ix);
      
	  dx =  xc - ix;
	  dy =  yc - iy;

	  if(!dx) dx=dy;
	  if(!dy) dy=dx;
	  	  
	  Float_t A;
	  GetReconstructionManager()->AG(E,dx,dy,A,gx,gy);
	  Float_t dd;
	  Info("EvaluatePosition", "  (ix iy  xc yc  dx dy)   %f %f %f %f %f %f", ix, iy, xc, yc, dx, dy) ;
	  Float_t chi2dg = GetReconstructionManager()->OneGamChi2(A,eDigit,E,dd);

  	  // Exclude digit with too large chisquare.
	  if(chi2dg > 10) {  continue; }

	  chisqc += chi2dg;
	  grxc += gx*dd;
	  gryc += gy*dd;

	}

      Float_t grc = TMath::Sqrt(grxc*grxc + gryc*gryc);
      if(grc<1.e-10) grc=1.e-10;
      Float_t sc = 1. + chisqc/chisq;
       Info("EvaluatePosition", " chisq: %f", chisq) ;
      st = st/sc;
      if(chisqc>chisq)
	goto loop20;
      cosi = (grx*grxc + gry*gryc)/gr/grc;
      st = st*sc/(1.4 - cosi);
      x1 = xc;
      y1 = yc;
      grx = grxc;
      gry = gryc;
      gr = grc;
      chisq = chisqc;

      if(chisq<chstop)
	goto loop101;
  
    loop20: ;
      Float_t step = st*gr;

      Info("EvaluatePosition", " Iteration %d dof %d chisq/dof %f chstop/dof %f step %d Stpmin %d",
	   Iter, dof, chisq/dof, chisq/dof, chstop/dof, step, Stpmin) ;
	

      if(step<Stpmin)
	goto loop101;

      if(step>1.)
	st = 1./gr;
      xc = x1 - st*grx;
      yc = y1 - st*gry;
    }


 loop101:
  chisq = chisq/dof;
//    if(chisq>Chsqcut)
//      {
//  //        TwoGam();
//      }

  Float_t gamma1[3];
  gamma1[0] = E;
  gamma1[1] = y1;
  gamma1[2] = x1;

  TVector3 newpos(gamma1[1],0,gamma1[2]);
  //SetLocalPosition(newpos);
  
  fLocPos=newpos;
  fAmp=E;

//    TVector3 pos;
//    RP->GetLocalPosition(pos);

  
  
  fChi2Dof = chisq;

} 


Bool_t AliPHOSEvalRecPoint::KillWeakPoint()
{
  // Remove this point from further procession
  // if it's energy is too small.
  
  Float_t Thr0 = GetReconstructionManager()->KillGamMinEnergy();
  
  if(GetEnergy()<Thr0) {
    Info("KillWeakPoint", "+++++++ Killing this rec point ++++++++++") ;
    RemoveFromWorkingPool(this);
    return kTRUE;
  }

  return kFALSE;
}


void AliPHOSEvalRecPoint::SplitMergedShowers() 
{
  // Split two very close rec. points.

  if(GetMultiplicity()<2) 
    return; 

  Float_t gamma1[3];
  Float_t gamma2[3];

  InitTwoGam(gamma1,gamma2); // start points for fit
  TwoGam(gamma1,gamma2); // evaluate the positions and energies
  UnfoldTwoMergedPoints(gamma1,gamma2); // make the job
  
}


void AliPHOSEvalRecPoint::MergeClosePoint() 
{

  AliPHOSLoader* fLoader = AliPHOSLoader::GetPHOSLoader(fEventFolderName);
//    AliPHOSDigitizer* digitizer = fLoader->Digitizer();
//    Float_t fPedestal = digitizer->GetPedestal();  // YVK 30.09.2001
//    Float_t fSlope = digitizer->GetSlope();
  
  for (Int_t i=0;i<InWorkingPool(); i++)
    {
      TObject* obj = GetFromWorkingPool(i);
      if(!((TObject*)this)->IsEqual(obj))
	{
	  AliPHOSRecPoint* rp = (AliPHOSRecPoint*)obj;
	  if(GetPHOSMod() == rp->GetPHOSMod())
	    {
	      if(TooClose(rp))
		{
		  Info("MergeClosePoint", "+++++++ Merging point 1: ++++++") ;
		  this->Print("");
		  Info("MergeClosePoint", "+++++++ and point 2: ++++++++++") ;
		  ((AliPHOSEvalRecPoint*)rp)->Print("");

		  //merge two rec. points
		  TVector3 lpos1;
		  TVector3 lpos2;
		  this->GetLocalPosition(lpos1);
		  rp->GetLocalPosition(lpos2);

		  Int_t* Digits = rp->GetDigitsList();
		  Float_t dE = rp->GetEnergy()/(rp->GetEnergy()+this->GetEnergy());
		  Float_t new_x = lpos1.X()*dE + lpos2.X()*(1.-dE); 
		  Float_t new_z = lpos1.Z()*dE + lpos2.Z()*(1.-dE);
		  Float_t new_E = rp->GetEnergy()+this->GetEnergy();
		  Float_t* Energies = ((AliPHOSEmcRecPoint*)rp)->GetEnergiesList();

		  for(Int_t iDigit=0; iDigit<rp->GetDigitsMultiplicity(); iDigit++)
		    {
		      AliPHOSDigit* digit = (AliPHOSDigit*)fLoader->Digits()->At(Digits[iDigit]);
		      Float_t eDigit = Energies[iDigit];
		      this->AddDigit(*digit,eDigit);
		    }

		  TVector3 newpos(new_x,0,new_z);
		  fLocPos = newpos;
		  fAmp = new_E;
		  RemoveFromWorkingPool(rp);
		  delete rp;
		  
		  Info("MergeClosePoint", "++++++ Resulting point: ++++++++") ;
		  this->Print("");

		  break;
		} 
	    }
	}
    }

}

Int_t AliPHOSEvalRecPoint::UnfoldLocalMaxima() 
{
  // Make unfolding in the reconstruction point with several local maxima.
  // Return the number of local maxima.

  // if multiplicity less then 2 - nothing to unfold
  if(GetMultiplicity()<2) return 1; 

  AliPHOSDigit * maxAt[1000];
  Float_t maxAtEnergy[1000];  
  Float_t LocMaxCut, LogWeight;
  Int_t relid[4] ; 
  Float_t xMax;
  Float_t zMax;

//    AliPHOSClusterizer* clusterizer = fLoader->Clusterizer("PHOSclu-v1");  
  AliPHOSClusterizer* clusterizer = GetClusterizer();
  if(!clusterizer) {
    Error("UnfoldLocalMaxima", "Cannot get clusterizer. Exit.") ;
    exit(1);
  }

  if(this->IsEmc()) {
    LocMaxCut = clusterizer->GetEmcLocalMaxCut();
    LogWeight = clusterizer->GetEmcLogWeight();
  }
  else {
    LocMaxCut = clusterizer->GetCpvLocalMaxCut();
    LogWeight = clusterizer->GetCpvLogWeight();
  }

  AliPHOSLoader* fLoader = AliPHOSLoader::GetPHOSLoader(fEventFolderName);
  const AliPHOSGeometry* fGeom = fLoader->PHOSGeometry();
  TClonesArray* digits = fLoader->Digits();

  // if number of local maxima less then 2 - nothing to unfold
  Int_t Nmax = GetNumberOfLocalMax(maxAt,maxAtEnergy,LocMaxCut,digits); 
  if(Nmax<2) return Nmax;

  Int_t* Digits = GetDigitsList();
  Int_t nDigits = GetMultiplicity();
  Float_t* Energies = GetEnergiesList();
  Float_t* eFit = new Float_t[nDigits];

  Int_t iDigit; //loop variable

  for(iDigit=0; iDigit<nDigits; iDigit++)
    {
      eFit[iDigit] = 0.;
    }

  for(iDigit=0; iDigit<nDigits; iDigit++)
    {
      
      AliPHOSDigit* digit = (AliPHOSDigit*)fLoader->Digits()->At( Digits[iDigit] );
      fGeom->AbsToRelNumbering(digit->GetId(), relid) ;
      Float_t x,z;
      fGeom->RelPosInModule(relid, x, z);

      for(Int_t iMax=0; iMax<Nmax; iMax++)
	{
	  AliPHOSDigit* digitMax = maxAt[iMax];
	  Float_t eMax = maxAtEnergy[iMax]; 
	  fGeom->AbsToRelNumbering(digitMax->GetId(), relid) ;
	  fGeom->RelPosInModule(relid, xMax, zMax);
	  Float_t dx = xMax - x;
	  Float_t dz = zMax - z;
	  Float_t Amp,gx,gy;
	  GetReconstructionManager()->AG(eMax,dz,dx,Amp,gx,gy); 
//        Amp = Amp + 0.5;
	  eFit[iDigit] += Amp;
	}
    }


  for(Int_t iMax=0; iMax<Nmax; iMax++) 
    {
      AliPHOSDigit* digitMax = maxAt[iMax];
      fGeom->AbsToRelNumbering(digitMax->GetId(), relid) ;
      fGeom->RelPosInModule(relid, xMax, zMax);
      Float_t eMax = maxAtEnergy[iMax];

      AliPHOSEvalRecPoint* newRP = new AliPHOSEvalRecPoint(IsCPV(),this);    
      newRP->AddDigit(*digitMax,maxAtEnergy[iMax]);

      //Neighbous ( matrix 3x3 around the local maximum)
      for(Int_t iDigit=0; iDigit<nDigits; iDigit++)
	{     
	  AliPHOSDigit* digit = (AliPHOSDigit*)fLoader->Digits()->At( Digits[iDigit] );
  	  Float_t eDigit = Energies[iDigit];
	  fGeom->AbsToRelNumbering(digit->GetId(), relid) ;
	  Float_t ix,iz;
	  fGeom->RelPosInModule(relid, ix, iz);

  	  if(AreNeighbours(digitMax,digit))
  	    {
	      Float_t dx = xMax - ix;
	      Float_t dz = zMax - iz;
	      Float_t single_shower_gain,gxMax,gyMax;
	      GetReconstructionManager()->AG(eMax,dz,dx,single_shower_gain,gxMax,gyMax);
  	      Float_t total_gain = eFit[iDigit];
    	      Float_t ratio = single_shower_gain/total_gain; 
	      Info("UnfoldLocalMaxima", " ratio -> %f", ratio) ;
  	      eDigit = eDigit*ratio;
	      newRP->AddDigit(*digit,eDigit);
  	    }
	}

      newRP->EvalLocalPosition(LogWeight,digits);
      Info("UnfoldLocalMaxima", "======= Unfold: daughter rec point %d =================", iMax) ;
      newRP->Print("");
    }

//    RemoveFromWorkingPool(this);

  delete[] eFit;

  return Nmax;

}

void AliPHOSEvalRecPoint::PrintPoint(Option_t* opt)
{
  AliPHOSCpvRecPoint::Print(opt);

  TVector3 lpos;
  GetLocalPosition(lpos);

  TString message ; 
  message  = "       Chi2/dof = %f" ;
  message += "       Local (x,z) = (%f, %f) in module %d" ; 
  Info("Print", message.Data(), Chi2Dof(), lpos.X(), lpos.Z(), GetPHOSMod()) ;

}

AliPHOSRecManager* AliPHOSEvalRecPoint::GetReconstructionManager() const
{

  TFolder* aliceF = (TFolder*)gROOT->FindObjectAny("YSAlice");
  TFolder* wPoolF = (TFolder*)aliceF->FindObject("WhiteBoard/RecPoints/PHOS/SmP");
  AliPHOSRecManager* recMng = (AliPHOSRecManager*)wPoolF->FindObject("AliPHOSRecManager");
  if(!recMng) { 
    Error("GetReconstructionManager", "Couldn't find Reconstruction Manager. Exit.") ; 
    exit(1); 
  }

  return recMng;
}


AliPHOSRecPoint* AliPHOSEvalRecPoint::Parent()
{
  if(fParent<0) return NULL;
  else
    return (AliPHOSRecPoint*)GetFromWorkingPool(fParent);
//    return fParent;
}

Float_t AliPHOSEvalRecPoint::Chi2Dof() const 
{
  return fChi2Dof;
}

const TObject* AliPHOSEvalRecPoint::GetWorkingPool()
{

  TFolder* aliceF = (TFolder*)gROOT->FindObjectAny("YSAlice");
  TFolder* wPoolF = (TFolder*)aliceF->FindObject("WhiteBoard/RecPoints/PHOS/SmP");
  TObject* wPool = wPoolF->FindObject("SmartPoints");
  if(!wPool) { 
    Error("GetWorkingPool", "Couldn't find Working Pool. Exit.") ; 
    exit(1); 
  }

  return wPool;
}

void AliPHOSEvalRecPoint::AddToWorkingPool(TObject* obj)
{
  ((TObjArray*)GetWorkingPool())->Add(obj);
}

TObject* AliPHOSEvalRecPoint::GetFromWorkingPool(Int_t i)
{
//    return fWorkingPool->At(i);
  return ((TObjArray*)GetWorkingPool())->At(i);
}

Int_t AliPHOSEvalRecPoint::InWorkingPool()
{
  return ((TObjArray*)GetWorkingPool())->GetEntries();
}

void AliPHOSEvalRecPoint::RemoveFromWorkingPool(TObject* obj)
{
  ((TObjArray*)GetWorkingPool())->Remove(obj);
  ((TObjArray*)GetWorkingPool())->Compress();
}

void AliPHOSEvalRecPoint::PrintWorkingPool()
{
  ((TObjArray*)GetWorkingPool())->Print("");
}

void AliPHOSEvalRecPoint::SetEventFolderName(const char* evfname)
{//Sets event folder name
  fEventFolderName = evfname;
}














