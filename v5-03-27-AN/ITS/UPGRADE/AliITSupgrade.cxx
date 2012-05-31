// **************************************************************************
// * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
// *                                                                        *
// * Author: The ALICE Off-line Project.                                    *
// * Contributors are mentioned in the code where appropriate.              *
// *                                                                        *
// * Permission to use, copy, modify and distribute this software and its   *
// * documentation strictly for non-commercial purposes is hereby granted   *
// * without fee, provided that the above copyright notice appears in all   *
// * copies and that both the copyright notice and this permission notice   *
// * appear in the supporting documentation. The authors make no claims     *
// * about the suitability of this software for any purpose. It is          *
// * provided "as is" without express or implied warranty.                  *
// **************************************************************************

/* $Id$ */

#include <TTree.h>
#include <TArrayD.h>            //new constructor
#include <TFile.h>  
#include <TGeoManager.h>    
#include <TGeoVolume.h>         //CreateGeometry()
#include <TGeoMatrix.h>
#include <TVirtualMC.h>         //->gMC in StepManager
#include <TPDGCode.h>           //StepHistory
#include <TClonesArray.h>
#include <TGeoGlobalMagField.h>
#include "AliRun.h"             //CreateMaterials()    
#include "AliMC.h"             //CreateMaterials()    
#include "AliMagF.h"            //CreateMaterials()
#include "AliITSupgrade.h"      //class header
#include "AliITShit.h"           // StepManager()
#include "AliITSDigitUpgrade.h"
#include "AliTrackReference.h"   // StepManager()
#include "AliITSsegmentationUpgrade.h"

ClassImp(AliITSupgrade) 
//__________________________________________________________________________________________________
  AliITSupgrade::AliITSupgrade():
    AliITS(),
    fWidths(0),
    fRadii(0),
    fNSectors(20),
    fRadiiCu(0),
    fWidthsCu(0),
    fCopper(0),
    fBeampipe(0), 
    fRadiusBP(0),
    fWidthBP(0),
    fHalfLengthBP(0),
    fNlayers(0),
    fHalfLength(0),
    fSdigits(0),
    fDigitArray(0),
    fSegmentation(0x0)
{
  //
  // default ctor
  //
}

//__________________________________________________________________________________________________
AliITSupgrade::AliITSupgrade(const char *name,const char *title, Bool_t isBeamPipe):
  AliITS(name,title),
  fWidths(0),
  fRadii(0),
  fNSectors(20),
  fRadiiCu(0),
  fWidthsCu(0),
  fCopper(0),
  fBeampipe(isBeamPipe),
  fRadiusBP(0),
  fWidthBP(0),
  fHalfLengthBP(0),
  fNlayers(0),
  fHalfLength(0),
  fSdigits(0),
  fDigitArray(0),
  fSegmentation(0x0)
{
  // Default values are used in order to simulate the standard ITS material budget (see The ALICE Collaboration et al 2008 JINST 3 S08002 - Figure 3.2).
  // The cell segmentation is chosen to achieve the tracking resolution as described in Table 3.2 of The ALICE Collaboration et al 2008 JINST 3 S08002,
  // apart from SPD0 where the 9 um value is considered : resolution*sqrt(12)
  // Cilinder lenghts are chosen according to Table 1 of the Alignment paper : http://iopscience.iop.org/1748-0221/5/03/P03003

  fNlayers = 6;
  fWidths.Set(fNlayers);
  fRadii.Set(fNlayers);
  fRadiiCu.Set(fNlayers);
  fWidthsCu.Set(fNlayers);
  fCopper.Set(fNlayers);
  fHalfLength.Set(fNlayers);


  Double_t xsz[6]={31.18*1e-04,41.6*1e-04,121.2e-04,121.4*1e-04,69.3*1e-04,69.3*1e-04};
  Double_t zsz[6]={416*1e-04,416*1e-04,97*1e-04,97*1e-04,2875*1e-04,2875*1e-04};
  Double_t halfL[6]={14.1,14.1,22.2,29.7,43.1,48.9};
  Double_t r[6]={4.,7.6,14.9,23.8,39.1,43.6};
  Double_t thick[6]={75.*1e-04,150.*1e-04,150.*1e-04,150.*1e-04,150.*1e-04,150.*1e-04};

  Int_t  nlayers =6;
  TArrayD xsizeSi(nlayers);
  TArrayD zsizeSi(nlayers);

  // adding beam pipe upgrade
  if(isBeamPipe){
    fRadiusBP = 2.0;
    fWidthBP = 0.08;
    fHalfLengthBP = 400;
  }

  // setting the geonetry parameters and the segmentation
  Int_t npixHalf[6], npixR[6];
  Double_t halfLmod[6], xszInt[6];
  Int_t c[6]={1,1,1,1,1,1};

  for(Int_t i=0; i<nlayers; i++){
    // recalc length 
    npixHalf[i]=(Int_t)(halfL[i]/zsz[i]);
    halfLmod[i]=npixHalf[i]*zsz[i];
    // recalc segmentation 
    npixR[i] = (Int_t)(2*TMath::Pi()*r[i]/xsz[i]);
    xszInt[i]= 2*TMath::Pi()*r[i]/npixR[i];
    xsizeSi.AddAt(xszInt[i],i);
    zsizeSi.AddAt(zsz[i],i);

    fHalfLength.AddAt(halfLmod[i],i);

    fRadii.AddAt(r[i],i);
    fWidths.AddAt(thick[i],i);

    fRadiiCu.AddAt(r[i]+thick[i],i);
    fWidthsCu.AddAt(0.015,i);//micron
    fCopper.AddAt(c[i],i);
  }

  SetFullSegmentation(xsizeSi,zsizeSi);
}
//__________________________________________________________________________________________________
AliITSupgrade::AliITSupgrade(const char *name,const char *title, TArrayD widths, TArrayD radii,TArrayD halfLengths, TArrayD radiiCu, TArrayD widthsCu, TArrayS copper,Bool_t bp,Double_t radiusBP, Double_t widthBP, Double_t halfLengthBP):
  AliITS(name,title),
  fWidths(0),
  fRadii(0),
  fNSectors(20),
  fRadiiCu(radiiCu),
  fWidthsCu(widthsCu),
  fCopper(copper),
  fBeampipe(bp),
  fRadiusBP(0),
  fWidthBP(0),
  fHalfLengthBP(0),
  fNlayers(widths.GetSize()),
  fHalfLength(halfLengths),
  fSdigits(0),
  fDigitArray(0),
  fSegmentation(0x0)
{
  //
  // constructor to set all the geometrical/segmentation specs from outside
  // 
  for(Int_t i=0;i<fNlayers;i++){
    fWidths.Set(i+1);
    fWidths.AddAt(widths.At(i),i);
    fRadii.Set(i+1);
    fRadii.AddAt(radii.At(i),i);
    AliDebug(1,"Creating Volume");
  }

  if(bp){
    fRadiusBP=radiusBP;
    fWidthBP=widthBP;
    fHalfLengthBP=halfLengthBP;
  }
  Init();

}

AliITSupgrade::~AliITSupgrade(){

  if(fSdigits) {fSdigits->Delete(); delete fSdigits;}
  if(fDigitArray)  {fDigitArray->Delete(); delete fDigitArray;}
}


 
//_________________________________________________________________________________________________
void AliITSupgrade::AddAlignableVolumes()const
{
  //not needed
  return;
}

//__________________________________________________________________________________________________

void AliITSupgrade::CreateMaterials()
{
  //
  // Definition of ITS materials  
  //
  
  AliInfo("Start ITS materials");
  //data from PDG booklet 2002     density [gr/cm^3] rad len [cm] abs len [cm]    
  Float_t   aAir[4]={12,14,16,36} ,   zAir[4]={6,7,8,18} ,   wAir[4]={0.000124,0.755267,0.231781,0.012827} , dAir=0.00120479; Int_t nAir=4;//mixture
  Float_t       aBe = 9.012 ,          zBe   = 4 ,            dBe   =  1.848    ,   radBe   =  65.19/dBe , absBe   = 77.8/dBe  ; // UPGRADE -> beryllium beampipe
  Float_t       aSi = 28.085 ,         zSi   = 14 ,           dSi   =  2.329    ,   radSi   =  21.82/dSi , absSi   = 108.4/dSi  ; // UPGRADE -> Si tube
  Float_t aCu = 63.54, zCu = 29 , dCu= 8.96 , radCu = 12.86/dCu, absCu= 137.3/dCu;//Upgrade -> Copper Tube

  Int_t   matId=0;                           //tmp material id number
  Int_t   unsens = 0, sens=1;                //sensitive or unsensitive medium
  Int_t   itgfld = 3;			     //type of field intergration 0 no field -1 user in guswim 1 Runge Kutta 2 helix 3 const field along z
  Float_t maxfld = 5.; 		             //max field value
  Float_t tmaxfd = -10.0;                    //max deflection angle due to magnetic field in one step
  Float_t deemax = - 0.2;                    //max fractional energy loss in one step   
  Float_t stemax = - 0.1;                    //max step allowed [cm]
  Float_t epsil  =  0.001;                  //abs tracking precision [cm]   
  Float_t stmin  = - 0.001;                  //min step size [cm] in continius process transport, negative value: choose it automatically
  Float_t tmaxfdSi = 0.1; // .10000E+01; // Degree
  Float_t stemaxSi = 0.0075; //  .10000E+01; // cm
  Float_t deemaxSi = 0.1; // 0.30000E-02; // Fraction of particle's energy 0<deemax<=1
  Float_t epsilSi  = 1.0E-4;// .10000E+01;
  Float_t stminSi  = 0.0; // cm "Default value used"

  Float_t epsilBe  = .001;    // Tracking precision,
  Float_t stemaxBe = -0.01;   // Maximum displacement for multiple scat
  Float_t tmaxfdBe = -20.;    // Maximum angle due to field deflection
  Float_t deemaxBe = -.3;     // Maximum fractional energy loss, DLS
  Float_t stminBe  = -.8;
  Int_t   isxfld = ((AliMagF*)TGeoGlobalMagField::Instance()->GetField())->Integ();//from CreateMaterials in STRUCT/AliPIPEv3.cxx
  Float_t sxmgmx = ((AliMagF*)TGeoGlobalMagField::Instance()->GetField())->Max();//from CreateMaterials in STRUCT/AliPIPEv3.cxx

      
  AliMixture(++matId,"UpgradeAir"  ,aAir  ,zAir  ,dAir  ,nAir  ,wAir  ); 
  AliMedium(kAir  ,"UpgradeAir"  ,matId, unsens, itgfld, maxfld, tmaxfd, stemax, deemax, epsil, stmin);
    
  AliMaterial(++matId,"UpgradeSi"  ,aSi  ,zSi  ,dSi  ,radSi  ,absSi  );  
  AliMedium(kSi  ,"UpgradeSi"  , matId, sens, isxfld, sxmgmx, tmaxfdSi, stemaxSi, deemaxSi, epsilSi, stminSi);
    
  AliMaterial(++matId,"UpgradeBe"  ,aBe  ,zBe  ,dBe  ,radBe  ,absBe  );  
  AliMedium(kBe  ,"UpgradeBe"  , matId, unsens, isxfld, sxmgmx, tmaxfdBe, stemaxBe, deemaxBe, epsilBe, stminBe);

  AliMaterial(++matId, "UpgradeCu", aCu, zCu, dCu, radCu, absCu);
  AliMedium(kCu, "UpgradeCu", matId, unsens, itgfld, maxfld, tmaxfd, stemax, deemax, epsil, stmin);
   
  AliInfo("End ITS materials");
          
}//void AliITS::CreateMaterials()
//__________________________________________________________________________________________________
void AliITSupgrade::CreateGeometry()
{
  //
  //Creates detailed geometry simulation (currently GEANT volumes tree)        
  //
  AliInfo("Start ITS upgrade preliminary version building");
  if(!gMC->IsRootGeometrySupported()) return;                
  TGeoVolumeAssembly *vol= CreateVol();
  gGeoManager->GetVolume("ALIC")->AddNode(vol,0);
  AliInfo("Stop ITS upgrade preliminary version building");
  PrintSummary();
} 
//__________________________________________________________________________________________________
void AliITSupgrade::StepManager()
{
  // Full Step Manager.
  // Arguments: none
  //   Returns: none           
  // StepHistory(); return; //uncomment to print tracks history
  //  StepCount(); return;   //uncomment to count photons

  if(!fSegmentation) AliFatal("No segmentation available");

  if(!(this->IsActive())) return;
  if(!(gMC->TrackCharge())) return;
  TString volumeName=gMC->CurrentVolName();
  if(volumeName.Contains("Be")) return;
  if(volumeName.Contains("Cu")) return;
  if(gMC->IsTrackExiting()) {
    AddTrackReference(gAlice->GetMCApp()->GetCurrentTrackNumber(), AliTrackReference::kITS);
  } // if Outer ITS mother Volume

  static TLorentzVector position, momentum; // Saves on calls to construtors
  static AliITShit hit;// Saves on calls to constructors

  
  Int_t  status = 0;
  
  //
  // Track status
  if(gMC->IsTrackInside())      status +=  1;
  if(gMC->IsTrackEntering())    status +=  2;
  if(gMC->IsTrackExiting())     status +=  4;
  if(gMC->IsTrackOut())         status +=  8;
  if(gMC->IsTrackDisappeared()) status += 16;
  if(gMC->IsTrackStop())        status += 32;
  if(gMC->IsTrackAlive())       status += 64;

  //
  // Fill hit structure.
  //
  Int_t copy=-1;
  gMC->CurrentVolID(copy);   

  volumeName.Remove(0,12);          // remove letters to get the layer number
  hit.SetModule(fSegmentation->GetIdIndex(volumeName.Atoi(),copy)); // layer and sector information are together in the IdIndex (if copy=0 the idIndex is the layer)); 
  hit.SetTrack(gAlice->GetMCApp()->GetCurrentTrackNumber());
    
  gMC->TrackPosition(position);
  gMC->TrackMomentum(momentum);
  hit.SetPosition(position);
    
  hit.SetTime(gMC->TrackTime());
  hit.SetMomentum(momentum);
  hit.SetStatus(status);
  hit.SetEdep(gMC->Edep());
  hit.SetShunt(GetIshunt());
  if(gMC->IsTrackEntering()){
    hit.SetStartPosition(position);
    hit.SetStartTime(gMC->TrackTime());
    hit.SetStartStatus(status);
    return; // don't save entering hit.
  } 
  // Fill hit structure with this new hit.
    
  new((*fHits)[fNhits++]) AliITShit(hit); // Use Copy Construtor.
  // Save old position... for next hit.
  hit.SetStartPosition(position);
  hit.SetStartTime(gMC->TrackTime());
  hit.SetStartStatus(status);
  return;

}
//__________________________________________________________________________________________________
TGeoVolumeAssembly * AliITSupgrade::CreateVol()
{
//
// method to create the Upgrade Geometry (silicon, copper cylinders and beampipe)
//  
  TGeoVolumeAssembly *vol = new TGeoVolumeAssembly("ITSupgrade");
  TGeoMedium *si   =gGeoManager->GetMedium("ITS_UpgradeSi");
  TGeoMedium *cu   =gGeoManager->GetMedium("ITS_UpgradeCu");
  TGeoMedium *be   =gGeoManager->GetMedium("ITS_UpgradeBe");
  for(Int_t ivol=0;ivol<fNlayers;ivol++){

    if(fNSectors<1){
      TGeoVolume *layer=gGeoManager->MakeTube(Form("LayerSilicon%i",ivol),si   ,    fRadii.At(ivol)   ,   fRadii.At(ivol)+fWidths.At(ivol)  ,   fHalfLength.At(ivol)); //upgraded situation

      TGeoVolume *layerCu=gGeoManager->MakeTube(Form("LayerCu%i",ivol),cu   ,    fRadiiCu.At(ivol)   ,   fRadiiCu.At(ivol)+fWidthsCu.At(ivol) ,  fHalfLength.At(ivol)  ); //upgraded situation

      vol ->AddNode(layer,ivol);
      if(fCopper.At(ivol)){
        vol->AddNode(layerCu,ivol);
      }
    }else{

      TGeoVolume *layer = gGeoManager->MakeTubs(Form("LayerSilicon%i",ivol),si,  fRadii.At(ivol),   fRadii.At(ivol)+fWidths.At(ivol)  ,fHalfLength.At(ivol),0,(360./fNSectors));
      TGeoVolume *layerCu = gGeoManager->MakeTubs(Form("LayerCu%i",ivol),cu   ,    fRadiiCu.At(ivol)   ,   fRadiiCu.At(ivol)+fWidthsCu.At(ivol) ,  fHalfLength.At(ivol) ,0,(360./fNSectors));


      for(Int_t i=0;i<fNSectors;i++){
        TGeoRotation *rot1 = new TGeoRotation(" ",0.0,0.0,360.*i/fNSectors);//sector rotation
        TGeoRotation *rot2 = new TGeoRotation(" ",0.0,0.0,360.*i/fNSectors);//
        vol->AddNode(layer,i,rot1);
        if(fCopper.At(ivol)){
          vol->AddNode(layerCu,i,rot2);
        }
      }
    }
  }
  
  
  if(fBeampipe) {
    TGeoVolume *beampipe=gGeoManager->MakeTube("BeamPipe", be   ,    fRadiusBP   ,  fRadiusBP+ fWidthBP ,  fHalfLengthBP  ); //upgraded situation
     vol->AddNode(beampipe,fNlayers);
   }
  return vol;
}
//_________________________________________________________________________________________________
void AliITSupgrade::SetFullSegmentation(TArrayD xsize,TArrayD zsize){
//
// Upgrade detector specs specified in used constructor are stored in a file .
// for further usage (loc<->master<->tracking reference system changes)

  Bool_t check=kFALSE;
  for(Int_t lay = 0; lay< xsize.GetSize(); lay++){
    Double_t arch = fRadii.At(lay)*(TMath::Pi()*2/fNSectors);
    Int_t nPixRPhi = (Int_t)(arch/xsize.At(lay));
    Int_t nPixZed = (Int_t)((2*fHalfLength.At(lay))/zsize.At(lay));
    if(nPixRPhi>9999)check=kTRUE;
    if(nPixZed>99999)check=kTRUE;
  }
  if(check) AliFatal(" Segmentation is too small!! ");
  if(fSegmentation) fSegmentation->SetNSectors(fNSectors);
  TArrayD nSect(1);
  nSect.AddAt(fNSectors,0);
  TFile *file= TFile::Open("Segmentation.root","recreate");
  file->WriteObjectAny(&xsize,"TArrayD","CellSizeX");
  file->WriteObjectAny(&zsize,"TArrayD","CellSizeZ");
  file->WriteObjectAny(&nSect,"TArrayD","nSectors");
  file->Close();
}
//_________________________________________________________________________________________________
void AliITSupgrade::StepHistory()
{ 
  // This methode is invoked from StepManager() in order to print out
  TString volumeName=gMC->CurrentVolName();
  if(!volumeName.Contains("Silicon")) return;
  static Int_t iStepN;
  const char *sParticle;
  switch(gMC->TrackPid()){
  case kProton:      sParticle="PROTON"    ;break;
  case kNeutron:     sParticle="neutron"   ;break;
  case kGamma:       sParticle="gamma"     ;break;
  case kPi0:         sParticle="Pi0"       ;break;
  case kPiPlus:      sParticle="Pi+"       ;break;
  case kPiMinus:     sParticle="Pi-"       ;break;
  case kElectron:    sParticle="electron"  ;break;
  default:           sParticle="not known" ;break;
  }

  TString flag="funny combination";
  if(gMC->IsTrackAlive()) {
    if(gMC->IsTrackEntering())      flag="enters to";
    else if(gMC->IsTrackExiting())  flag="exits from";
    else if(gMC->IsTrackInside())   flag="inside";
    else if(gMC->IsTrackStop())     flag="stopped in";
  }

  Int_t vid=0,copy=0;
  TString path=gMC->CurrentVolName(); path.Prepend("-");path.Prepend(gMC->CurrentVolOffName(1));//current volume and his mother are always there
  vid=gMC->CurrentVolOffID(2,copy);  if(vid) {path.Prepend("-");path.Prepend(gMC->VolName(vid));}
  vid=gMC->CurrentVolOffID(3,copy);  if(vid) {path.Prepend("-");path.Prepend(gMC->VolName(vid));}
  

  AliInfo(Form("\n Step %i: %s (%i) %s %s m=%.6f GeV q=%.1f dEdX=%.4f Etot=%.4f",iStepN,sParticle,gMC->TrackPid(),flag.Data(),path.Data(),gMC->TrackMass(),gMC->TrackCharge(),gMC->Edep()*1e9,gMC->Etot()));

  Double_t gMcTrackPos[3]; gMC->TrackPosition(gMcTrackPos[0],gMcTrackPos[1],gMcTrackPos[2]);
  Double_t  gMcTrackPosLoc[3]; gMC->Gmtod(gMcTrackPos,gMcTrackPosLoc,1);
  TString v(volumeName.Data());
  v.ReplaceAll("LayerSilicon","");
  Int_t ilayer = v.Atoi();
  Double_t rXY = TMath::Sqrt(gMcTrackPos[0]*gMcTrackPos[0]+gMcTrackPos[1]*gMcTrackPos[1]);
  AliInfo(Form("gMC Track Position (MARS) x: %5.3lf, y: %5.3lf, z: %5.3lf (r: %5.3lf) (deltaR %5.5lf - width %5.5f)",gMcTrackPos[0],gMcTrackPos[1],gMcTrackPos[2], rXY , rXY - (fRadii.At(ilayer)),fWidths.At(ilayer)));


  AliDebug(10,Form("Step %i: tid=%i flags alive=%i disap=%i enter=%i exit=%i inside=%i out=%i stop=%i new=%i",
		   iStepN, gAlice->GetMCApp()->GetCurrentTrackNumber(),
		   gMC->IsTrackAlive(), gMC->IsTrackDisappeared(),gMC->IsTrackEntering(), gMC->IsTrackExiting(),
		   gMC->IsTrackInside(),gMC->IsTrackOut(),        gMC->IsTrackStop(),     gMC->IsNewTrack()));

  Float_t a,z,den,rad,abs; a=z=den=rad=abs=-1;
  Int_t mid=gMC->CurrentMaterial(a,z,den,rad,abs);
  AliDebug(10, Form("Step %i: mid=%i a=%7.2f z=%7.2f den=%9.4f rad=%9.2f abs=%9.2f\n\n",iStepN,mid,a,z,den,rad,abs));
/*
  TArrayI proc;  gMC->StepProcesses(proc);

  AliInfo("Processes in this step:");
  for ( int i = 0 ; i < proc.GetSize(); i++)
    {
      printf(Form(" - %s - ",TMCProcessName[proc.At(i)]));
    }
  printf("\n");
  AliInfo("End process list");
*/
  iStepN++;
}//StepHistory()
//______________________________________________________
void AliITSupgrade::Hits2SDigits(){
  
  // Interface method invoked from AliSimulation to create a list of sdigits corresponding to list of hits. 
  // Every hit generates one or more sdigits.
  // Arguments: none
  //   Returns: none
  AliDebug(1,"Start Hits2SDigits.");
  
  if(!fSegmentation)fSegmentation=new AliITSsegmentationUpgrade();
 
  if(!GetLoader()->TreeH()) {GetLoader()->LoadHits();}
  if(!GetLoader()->TreeS())

    for(Int_t iEvt=0;iEvt < GetLoader()->GetRunLoader()->GetNumberOfEvents();iEvt++){//events loop
      GetLoader()->GetRunLoader()->GetEvent(iEvt);                          //get next event
      {
	GetLoader()->MakeTree("S");
	MakeBranch("S");
      }
      SetTreeAddress();
      Int_t nSdigit[10]; for(Int_t i=0;i<10;i++)  nSdigit[i] =0; 
      for(Int_t iEnt=0;iEnt<GetLoader()->TreeH()->GetEntries();iEnt++){ // prims loop
        GetLoader()->TreeH()->GetEntry(iEnt);    
        Hit2SumDig(Hits(),SDigitsList(),nSdigit);//convert this hit to list of sdigits  
      }//prims loop
      GetLoader()->TreeS()->Fill();
      GetLoader()->WriteSDigits("OVERWRITE");
      SDigitsReset();
    }//events loop
  GetLoader()->UnloadHits();
  GetLoader()->UnloadSDigits();
  AliDebug(1,"Stop Hits2SDigits.");
  
}
//____________________________
void AliITSupgrade::Hit2SumDig(TClonesArray *hits,const TObjArray *pSDig, Int_t *nSdigit)
{
  // Adds  sdigits of this hit to the list
  //   Returns: none
  
  AliDebug(1,"Start Hit2SumDig");
    
  if(!fSegmentation){    
    AliDebug(1,"Segmentation Not inizialized!!");
    return ;
  }

  TClonesArray *pSdigList[8]; // is is the max number of layers allowed in the tracking 
  for(Int_t il=0; il<8; il++) pSdigList[il]=0x0; 

  for(Int_t i=0;i<fNlayers;i++){ 
    pSdigList[i]=(TClonesArray*)(*pSDig)[i];
    if(pSdigList[i]->GetEntries()!=0 && nSdigit[i]==0) {
      AliDebug(1,Form("Entries of pSdigList %d;   layer: %d,",pSdigList[i]->GetEntries(),i));
      AliErrorClass(" -> Some of sdigits lists is not empty");         //in principle those lists should be empty 
    }
  }
  
  for(Int_t iHit=0;iHit<hits->GetEntries();iHit++){         //hits loop
    AliITShit *hit = (AliITShit*)hits->At(iHit);
    Double_t xz[2];

    Int_t module=99;
    if(!fSegmentation->GlobalToDet(fSegmentation->GetLayerFromIdIndex(hit->GetModule()),hit->GetXG(),hit->GetYG(),hit->GetZG(),xz[0],xz[1],module)) continue;
    AliITSDigitUpgrade digit;
    digit.SetSignal(hit->GetIonization());
    digit.SetNelectrons(hit->GetIonization()/(3.62*1e-09));
    digit.SetLayer(fSegmentation->GetLayerFromIdIndex(hit->GetModule()));
    digit.SetModule(fSegmentation->GetModuleFromIdIndex(hit->GetModule()));//set the module (=sector) of ITSupgrade 
    
    digit.SetTrackID(hit->GetTrack()); 
    
    Int_t xpix = 999;
    Int_t zpix = 999; // shift at the next line to have zpix always positive. Such numbers are used to build the Pixel Id in the layer (> 0!) 
    fSegmentation->DetToPixID(xz[0], xz[1],fSegmentation->GetLayerFromIdIndex(hit->GetModule()), xpix, zpix);
    digit.SetPixId(xpix,zpix);
    new((*pSdigList[fSegmentation->GetLayerFromIdIndex(hit->GetModule())])[nSdigit[fSegmentation->GetLayerFromIdIndex(hit->GetModule())]++]) AliITSDigitUpgrade(digit);

  }
  
  AliDebug(1,"Stop Hit2SumDig.");
}
//_______________________________________________________________________________________________
void AliITSupgrade::MakeBranch(Option_t *option){
  //Create Tree branches 
  AliDebug(1,Form("Start with option= %s.",option));
  
  const Int_t kBufSize = 4000;
  
  const char *cH = strstr(option,"H");
  const char *cD = strstr(option,"D");
  //const char *cR = strstr(option,"R");
  const char *cS = strstr(option,"S");

  if(cH&&fLoader->TreeH()){

    HitCreate();
    MakeBranchInTree(fLoader->TreeH(),"ITSupgrade",&fHits,kBufSize,0);   
  }

    
  if(cS&&fLoader->TreeS()){
    SDigitsCreate();
    for(Int_t i=0;i<fNlayers;i++) MakeBranchInTree(fLoader->TreeS(),Form("Layer%d",i),&((*fSdigits)[i]),kBufSize,0);
  }


  if(cD&&fLoader->TreeD()){
    DigitsCreate();
    for(Int_t i=0;i<fNlayers;i++) MakeBranchInTree(fLoader->TreeD(),Form("Layer%d",i),&((*fDigitArray)[i]),kBufSize,0);
  }
  AliDebug(1,"Stop.");
}
//____________________________________________________________________________________________________ 
void AliITSupgrade::SetTreeAddress()
{
  //Set branch address for the Hits and Digits Tree.
  AliDebug(1,"Start.");
  if(fLoader->TreeH() && fLoader->TreeH()->GetBranch("ITSupgrade" )){
    HitCreate();
    fLoader->TreeH()->SetBranchAddress("ITSupgrade",&fHits);

  }
    
  if(fLoader->TreeS() && fLoader->TreeS()->GetBranch("Layer0" )){
    SDigitsCreate();
    for(int i=0;i<fNlayers;i++){
      fLoader->TreeS()->SetBranchAddress(Form("Layer%d",i),&((*fSdigits)[i]));
    }
  }
    
  if(fLoader->TreeD() && fLoader->TreeD()->GetBranch("Layer0")){
    DigitsCreate(); 
    for(int i=0;i<fNlayers;i++){
      fLoader->TreeD()->SetBranchAddress(Form("Layer%d",i),&((*fDigitArray)[i]));
    }
  }
  AliDebug(1,"Stop.");
}
//____________________________________________________________________________________________________ 
void AliITSupgrade::PrintSummary()
{

  // pdg booklet or http://pdg.lbl.gov/2010/reviews/rpp2010-rev-atomic-nuclear-prop.pdf

  // X0 as X0[g/cm^2]/density 
  Double_t beX0= 65.19/1.848; 
  Double_t  cuX0= 12.86/8.96; 
  Double_t siX0= 21.82/2.329;

  if(!fSegmentation) fSegmentation= new AliITSsegmentationUpgrade();
  printf(" %10s %10s %10s %10s %10s %10s %10s %10s\n","Name","R [cm]","x/X0 [%]","Phi res[um]","Z res[um]","X segm[um]","Y segm[um]","widths [um]");
  
   if(fBeampipe) printf(" Beampipe   %10.3f %10f %10s %10s %10s %10s %10.1f\n", fRadiusBP, 100.*fWidthBP/beX0,"-","-","-","-",fWidthBP*1.e+4); 
   else printf(" * No Beam Pipe Upgrade  \n");
  for(Int_t i=0; i< fNlayers; i++){
    printf(" Si Layer %1i %10.3f %10f %10.2f %10.2f %10.2f %10.2f %10.1f\n",i, fRadii.At(i), 100.*fWidths.At(i)/siX0, (fSegmentation->GetCellSizeX(i) *1.e+4)/TMath::Sqrt(12.), (fSegmentation->GetCellSizeZ(i)*1.e+4)/TMath::Sqrt(12.),fSegmentation->GetCellSizeX(i) *1.e+4, fSegmentation->GetCellSizeZ(i)*1.e+4,fWidths.At(i)*1.e+4);
    printf(" Cu Layer   %10.3f %10f %10s %10s %10s %10s %10.1f\n", fRadiiCu.At(i), 100.*fWidthsCu.At(i)/cuX0,"-","-","-","-",fWidthsCu.At(i)*1.e+4);
  }
}
//____________________________________________________________________________________________________ 
void AliITSupgrade::SetSegmentationX(Double_t x, Int_t lay){
//
// Virtual dimension in X direction is written into the file
// (-> see SetFullSegmentation)

TFile *f = TFile::Open("Segmentation.root","UPDATE");
if(!f) { 
  AliError("\n\n\n Segmentation.root file does not exist. The segmentation is 0! \n\n\n");
  return;
 }
else {
  TArrayD *a = (TArrayD*)f->Get("CellSizeX");
  if(!a)  {
   AliError("CessSizeX array does not exist!");
   return; 
   }
  a->AddAt(x,lay);
  f->WriteObjectAny(a,"TArrayD","CellSizeX");
  f->Close();
  delete f;
 }
}
//____________________________________________________________________________________________________ 
void AliITSupgrade::SetSegmentationZ(Double_t z, Int_t lay){
//
// Virtual dimension in X direction is written into the file
// (-> see SetFullSegmentation)
//


TFile *f = TFile::Open("Segmentation.root","UPDATE");
if(!f) { 
  AliError("\n\n\n Segmentation.root file does not exist. The segmentation is 0! \n\n\n");
  return;
 }
else {
  TArrayD *a = (TArrayD*)f->Get("CellSizeZ");
  if(!a)  {
   AliError("CellSizeZ array does not exist!");
   return;
   }
  a->AddAt(z,lay);
  f->WriteObjectAny(a,"TArrayD","CellSizeZ");
  f->Close();
  delete f;
 }

}
