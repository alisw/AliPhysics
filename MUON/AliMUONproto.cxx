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
/*
$Log$
Revision 1.2  2000/06/15 07:58:49  morsch
Code from MUON-dev joined

Revision 1.1.2.1  2000/04/18 09:11:15  morsch
Implementation of MUON Chamber Prototype Class
Both read digits from raw data or use the Monte-Carlo.
Rachid GUERNANE, IPN Lyon guernane@ipnl.in2p3.fr

*/

/*
Implementation of MUON Chamber Prototype Class 
Both read digits from raw data or use the Monte-Carlo. 
1-Feb-2000 Rachid GUERNANE, IPN Lyon guernane@ipnl.in2p3.fr
*/


////////////////////////////////////////////////
//  Manager and hits classes for set:PROTO    //
////////////////////////////////////////////////

#include <TTUBE.h>
#include <TBRIK.h>
#include <TRotMatrix.h>
#include <TNode.h> 
#include <TTree.h> 
#include <TRandom.h> 
#include <TObject.h>
#include <TVector.h>
#include <TObjArray.h>
#include <TMinuit.h>
#include <TParticle.h>
#include <TROOT.h>
#include <TFile.h>
#include <TNtuple.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TDirectory.h>
#include <TObjectTable.h>
#include <AliPDG.h>
#include <TArrayI.h>
#include <AliDetector.h>

#include "AliMUONChamber.h"
#include "AliMUONproto.h"
#include "AliMUONHit.h"
#include "TTUBE.h"
#include "AliMUONClusterFinder.h"
#include "AliRun.h"
#include "AliMC.h"
#include "iostream.h"
#include "AliCallf77.h" 
#include "AliConst.h"
//#include "chainalice2.h"
#include "AliMUONSegmentationV0.h"
//#include "AliMUONSegResV11.h"

ClassImp(AliMUONproto)

//___________________________________________    
AliMUONproto::AliMUONproto()
        : AliMUON()
{
    cout << "\n Calling AliMUONproto constructor..." << endl;
    
        //
        //
        //
}
 
//___________________________________________
AliMUONproto::AliMUONproto(const char *name, const char *title) 
        : AliMUON(name,title)
//  AliMUON defines everything, but the chamber for NCH=1
{
    
//  z-Positions of Chambers
    const Float_t zch[1] = {975.};
//
//  inner diameter
    const Float_t dmi[1] = {0.};
//
//  outer diameter
    const Float_t dma[1] = {29.2};
//
//    
//    Default Parameters for ALICE2 prototype

//
    (*fChambers)[0] = new AliMUONChamber();
    AliMUONChamber* chamber = (AliMUONChamber*) (*fChambers)[0];
    chamber->SetGid(0);
    chamber->SetZ(zch[0]);
//
    chamber->InitGeo(zch[0]); 
    chamber->SetRInner(dmi[0]/2);
    chamber->SetROuter(dma[0]/2);

    for (int i = 0; i <= 99; i++) {
        fThreshold[i] = 0.;
    }
    
}

#ifdef WE_FORGET_FOR_THE_MOMENT
//___________________________________________ 
void AliMUONproto::GetRawDigits(Int_t evnb, Int_t *lptr, Int_t ilen) {

    Int_t ip = 0;
    Int_t equip = 0;
    Int_t nochnl;
    Int_t loop;
    Int_t val;    
    Int_t itype;
    Int_t id;
    Int_t serial;
    Int_t equiplen;

    Int_t digits[5]; 
    
    AliMUON* MUON  = (AliMUON*)gAlice->GetModule("MUON");
    AliMUONSegmentationV0* seg = (AliMUONSegmentationV0*) Chamber(0).SegmentationModel(1);
    
    Int_t tracks[10];
    Int_t charges[10];
        
    for (Int_t i = 0; i < 10; i++) {
        tracks[i] = 0;
        charges[i] = 0;
    }
        
    Int_t ich = 0;
    
  nwtype:	      

    itype = lptr[ip++];
    id = lptr[ip++];
    equiplen = lptr[ip++];
    
    if (equiplen < 0 ) return;

    if (itype != (int)(0XCACCA008)) {        
        ip += equiplen;
        if (ip < ilen-2)  goto nwtype;
    }    
    else {
        serial = id >> 16;
        equip = id & 0x1;
        if ((serial == 190) && (equip == 1)) {
            for (loop = 0; loop < equiplen; loop++) {
                nochnl = (lptr[ip] & 0x7ff000 ) >> 12;
                val = lptr[ip] & 0x3ff;
                // fill digits from raw data according to cathode connexions
                if (group[nochnl][2][1]!=0)
                    digits[0] = group[nochnl][2][1] - 12;
                else if  (group[nochnl][1][1]!=0)
                    digits[0] = group[nochnl][1][1] - 12;
                else
                    digits[0] = group[nochnl][0][1] - 12;
                digits[1] = group[nochnl][0][0];
                if (digits[0] != seg->Ix(digits[0], digits[1]))
                    printf("Pb pour ix=%d,iy=%d\n", digits[0], digits[1]);
                digits[2] = val;
                digits[3] = 0;
                digits[4] = 0;
                if (digits[2] >= fThreshold[nochnl]) 
                    MUON->AddDigits(ich, tracks, charges, digits);

                ip++;
            }   
        } 
        else
            ip += equiplen;
        
        if (ip < ilen-2) goto nwtype;
    }
    
    gAlice->TreeD()->Fill();
    MUON->ResetDigits();

    gAlice->TreeD()->Fill();
    MUON->ResetDigits();
    
    char hname[30];
    sprintf(hname, "TreeD%d", evnb);
    gAlice->TreeD()->Write(hname);
    // reset tree
    gAlice->TreeD()->Reset();

}
#endif

//___________________________________________
void AliMUONproto::SetPadSize(Int_t id, Int_t isec, Float_t p1, Float_t p2)
{
    Int_t i=2*(id-1);
    ((AliMUONChamber*)(*fChambers)[i])->SetPadSize(isec,p1,p2);
}

//___________________________________________
void AliMUONproto::SetChargeSlope(Int_t id, Float_t p1)
{
    Int_t i=2*(id-1);
    ((AliMUONChamber*) (*fChambers)[i])->SetChargeSlope(p1);
}

//___________________________________________
void AliMUONproto::SetChargeSpread(Int_t id, Float_t p1, Float_t p2)
{
    Int_t i=2*(id-1);
    ((AliMUONChamber*) (*fChambers)[i])->SetChargeSpread(p1,p2);
}

//___________________________________________
void AliMUONproto::SetMaxAdc(Int_t id, Float_t p1)
{
    Int_t i=2*(id-1);
    ((AliMUONChamber*) (*fChambers)[i])->SetMaxAdc(p1);
}

//___________________________________________
void AliMUONproto::CreateGeometry()
{
    Int_t *idtmed = fIdtmed->GetArray()-1099;
//
    Float_t zpos;
    Float_t bpar[3];
    Int_t idGas=idtmed[1100];
    
    AliMUONChamber *iChamber;
//********************************************************************
//                          Prototype ALICE2                        **
//********************************************************************
    iChamber=(AliMUONChamber*) (*fChambers)[0];
    
    zpos=iChamber->Z(); 
//
    const Float_t X_POS = 11*.975/2; //half x
    const Float_t Y_POS = 18*.55/2;  //half y
    const Float_t Z_POS = 0.325;
    
    bpar[0] = X_POS; 
    bpar[1] = Y_POS;  
    bpar[2] = Z_POS;
    
    gMC->Gsvolu("C01G", "BOX", idGas, bpar, 3);
    gMC->Gspos("C01G", 1, "ALIC", -bpar[0], bpar[1], zpos, 0, "ONLY");
}

//___________________________________________
void AliMUONproto::CreateMaterials()
{
  // *** DEFINITION OF AVAILABLE MUON MATERIALS *** 
  //
  //     Ar-CO2 gas 
    Float_t ag1[3]   = { 39.95,12.01,16. };
    Float_t zg1[3]   = { 18.,6.,8. };
    Float_t wg1[3]   = { .8,.0667,.13333 };
    Float_t dg1      = .001821;
    
    Float_t epsil  = .001; // Tracking precision, 
    Float_t tmaxfd = -20.; // Maximum angle due to field deflection
    Float_t stmin  = -.8;
    Int_t ISXFLD   = gAlice->Field()->Integ();
    Float_t SXMGMX = gAlice->Field()->Max();

    //
    // --- Define the various materials for GEANT ---
    AliMixture(22, "ArCO2 80%$", ag1, zg1, dg1, 3, wg1);

    //
    //    Ar/CO2
    AliMedium(1, "ARG_CO2   ", 22, 1, ISXFLD, SXMGMX, tmaxfd, fMaxStepGas,
              fMaxDestepAlu, epsil, stmin);
    //    Air 
        //AliMedium(1, "AIR_CH_US         ", 15, 1, ISXFLD, SXMGMX, tmaxfd, -1., -.3, epsil, stmin);
}

//___________________________________________
void AliMUONproto::Init()
{
   printf("\n\n\n Start Init for Prototype ALICE2 - CPC chamber type\n\n\n");

   // 
   // Initialize Tracking Chambers
   //
   
   ( (AliMUONChamber*) (*fChambers)[0])->Init();
   
   //
   // Set the chamber (sensitive region) GEANT identifier
   AliMC* gMC = AliMC::GetMC(); 
   ((AliMUONChamber*)(*fChambers)[0])->SetGid(gMC->VolId("C01G"));

   printf("\n\n\n Finished Init for Prototype ALICE2 - CPC chamber type\n\n\n");

}

//___________________________________________
void AliMUONproto::StepManager()
{
  Int_t          copy, id;
  static Int_t   idvol;
  static Int_t   vol[2];
  Int_t          ipart;
  TLorentzVector pos;
  TLorentzVector mom;
  Float_t        theta,phi;
  Float_t        destep, step;
  
  static Float_t eloss, eloss2, xhit, yhit, tlength;
  const  Float_t big=1.e10;
  
  //  modifs perso
  static Float_t hits[14];

  TClonesArray &lhits = *fHits;

  //
  // Set maximum step size for gas
  // numed=gMC->GetMedium();
  //
  // Only charged tracks
  if( !(gMC->TrackCharge()) ) return; 
  //
  // Only gas gap inside chamber
  // Tag chambers and record hits when track enters 
  idvol=-1;
  id=gMC->CurrentVolID(copy);
  
    for (Int_t i=1; i<=kNCH; i++) {
      if(id==((AliMUONChamber*)(*fChambers)[i-1])->GetGid()){ 
	  vol[0]=i; 
	  idvol=i-1;
      }
    }
    if (idvol == -1) return;
  //
  // Get current particle id (ipart), track position (pos)  and momentum (mom) 
  gMC->TrackPosition(pos);
  gMC->TrackMomentum(mom);

  ipart  = gMC->TrackPid();
  //Int_t ipart1 = gMC->IdFromPDG(ipart);
  //printf("ich, ipart %d %d \n",vol[0],ipart1);

  //
  // momentum loss and steplength in last step
  destep = gMC->Edep();
  step   = gMC->TrackStep();
  
  //
  // record hits when track enters ...
  if( gMC->IsTrackEntering()) {
      gMC->SetMaxStep(fMaxStepGas);
      Double_t tc = mom[0]*mom[0]+mom[1]*mom[1];
      Double_t rt = TMath::Sqrt(tc);
      Double_t pmom = TMath::Sqrt(tc+mom[2]*mom[2]);
      Double_t tx=mom[0]/pmom;
      Double_t ty=mom[1]/pmom;
      Double_t tz=mom[2]/pmom;
      Double_t s=((AliMUONChamber*)(*fChambers)[idvol])
	  ->ResponseModel()
	  ->Pitch()/tz;
      theta   = Float_t(TMath::ATan2(rt,Double_t(mom[2])))*kRaddeg;
      phi     = Float_t(TMath::ATan2(Double_t(mom[1]),Double_t(mom[0])))*kRaddeg;
      hits[0] = Float_t(ipart);         // Geant3 particle type
      hits[1] = pos[0]+s*tx;                 // X-position for hit
      hits[2] = pos[1]+s*ty;                 // Y-position for hit
      hits[3] = pos[2]+s*tz;                 // Z-position for hit
      hits[4] = theta;                  // theta angle of incidence
      hits[5] = phi;                    // phi angle of incidence 
      hits[8] = (Float_t) fNPadHits;   // first padhit
      hits[9] = -1;                     // last pad hit

      // modifs perso
      hits[10] = mom[3]; // hit momentum P
      hits[11] = mom[0]; // Px/P
      hits[12] = mom[1]; // Py/P
      hits[13] = mom[2]; // Pz/P
      // fin modifs perso

      // phi angle of incidence
      tlength = 0;
      eloss   = 0;
      eloss2  = 0;
      xhit    = pos[0];
      yhit    = pos[1];      
      // Only if not trigger chamber
      if(idvol<10) {
	  //
	  //  Initialize hit position (cursor) in the segmentation model 
	  ((AliMUONChamber*) (*fChambers)[idvol])
	      ->SigGenInit(pos[0], pos[1], pos[2]);
      } else {
	  //geant3->Gpcxyz();
	  //printf("In the Trigger Chamber #%d\n",idvol-9);
      }
  }
  eloss2+=destep;
  
  // 
  // Calculate the charge induced on a pad (disintegration) in case 
  //
  // Mip left chamber ...
  if( gMC->IsTrackExiting() || gMC->IsTrackStop() || gMC->IsTrackDisappeared()){
      gMC->SetMaxStep(big);
      eloss   += destep;
      tlength += step;
      
      // Only if not trigger chamber
      if(idvol<10) {
	  if (eloss > 0) 
	      MakePadHits(0.5*(xhit+pos[0]),0.5*(yhit+pos[1]),eloss,0.0,idvol);
      }
	  
      hits[6]=tlength;
      hits[7]=eloss2;
      if (fNPadHits > (Int_t)hits[8]) {
	  hits[8]= hits[8]+1;
	  hits[9]= (Float_t) fNPadHits;
      }
    
      new(lhits[fNhits++]) 
	  AliMUONHit(fIshunt,gAlice->CurrentTrack(),vol,hits);
      eloss = 0; 
      //
      // Check additional signal generation conditions 
      // defined by the segmentation
      // model (boundary crossing conditions) 
  } else if 
      (((AliMUONChamber*) (*fChambers)[idvol])
       ->SigGenCond(pos[0], pos[1], pos[2]))
  {
      ((AliMUONChamber*) (*fChambers)[idvol])
	  ->SigGenInit(pos[0], pos[1], pos[2]);
//      printf("\n-> MakePadHits, reason special %d",ipart);
      if (eloss > 0)
	  MakePadHits(0.5*(xhit+pos[0]),0.5*(yhit+pos[1]),eloss,0.0,idvol);
      xhit     = pos[0];
      yhit     = pos[1]; 
      eloss    = destep;
      tlength += step ;
      //
      // nothing special  happened, add up energy loss
  } else {        
      eloss   += destep;
      tlength += step ;
  }
}

//___________________________________________
void AliMUONproto::BuildGeometry()
{
    TNode *Node;
    TNode *Top;
    
    const int kColorMUON = kBlue;
    //
    Top=gAlice->GetGeometry()->GetNode("alice");
    //
    //
    //
    float dx, dy, dz, zpos;

    const Float_t cz[1]={975.};
    
    zpos=cz[0];

    dx=11*.975/2;
    dy=18*.55/2;
    dz=0.325;

    new TBRIK("C_MUON101","Mother volume for Proto.","void",dx*2,dy*2,dz*2);
    TBRIK* PROTO = new TBRIK("PROT", "Proto. Sens. Region","void",dx,dy,dz);
    Top->cd();
    Node = new TNode("MUON101","ChamberNode","C_MUON101",0,0,zpos,"");
    Node->SetLineColor(kColorMUON);
    Node->SetVisibility(0);    
    fNodes->Add(Node);
    Node->cd();
    Node = new TNode("MUON201", "Proto. Sens. Region Node", PROTO, -dx, dy, dz);
    Node->SetLineColor(kColorMUON);    

}

//_____________________________________________________________________________
void AliMUONproto::FindClusters(Int_t nev, Int_t last_entry)
{

//
// Loop on chambers and on cathode planes
//
  for (Int_t icat = 0; icat < 2; icat++) {
	    gAlice->ResetDigits();
 	    gAlice->TreeD()->GetEvent(last_entry+icat);
            
      for (Int_t ich = 0; ich < kNTrackingCh; ich++) {
	  AliMUONChamber* iChamber=(AliMUONChamber*) (*fChambers)[ich];
	  TClonesArray *MUONdigits  = this->DigitsAddress(ich);
	  if (MUONdigits == 0) continue;
          //
	  // Get ready the current chamber stuff
	  //

          AliMUONResponse* response = iChamber->ResponseModel();
	  AliMUONSegmentation* seg = iChamber->SegmentationModel(icat+1);
	  AliMUONClusterFinder* rec = iChamber->ReconstructionModel();

          if (seg) {	  
	      rec->SetSegmentation(seg);
	      rec->SetResponse(response);
	      rec->SetDigits(MUONdigits);
	      rec->SetChamber(ich);
	      rec->FindRawClusters();
	  }
          
          TClonesArray *fRch;
	  fRch=RawClustAddress(ich);
	  fRch->Sort();
          // it seems to work
          
      } // for ich
      // fill the tree
      gAlice->TreeR()->Fill();

      ResetRawClusters();
  } // for icat

  char hname[30];
  sprintf(hname,"TreeR%d",nev);
  gAlice->TreeR()->Write(hname);
  gAlice->TreeR()->Reset();
}

#ifdef WE_FORGRET_THIS_SHIT
//_____________________________________________________________________________
void AliMUONproto::SetThreshold()
{

    ifstream inputFile("/home/alice/guernane/aliroot/pro/MUON/crped190.dat", ios::in);

    if (inputFile.fail()) {
        cout << "Error opening file" << endl;
        exit(2);
    }

    char buff[32];
    Int_t Serial;
    Int_t Ntrigger;
    Int_t Nchannel;
    Int_t i1;
    Int_t i2;
        
    inputFile >> buff;

    inputFile >> Serial;
    inputFile >> Ntrigger;
    inputFile >> Nchannel;
    inputFile >> i1;
    inputFile >> i2;

    Float_t ped0[Nchannel];
    Float_t sig0[Nchannel];
    Float_t ped1[Nchannel];
    Float_t sig1[Nchannel];
    Int_t ichannel;

    for (Int_t i = 0; i < Nchannel-1; i++) {
    ped0[i] = 0;
    sig0[i] = 0;
    ped1[i] = 0;
    sig1[i] = 0;        
    }
    
    for (Int_t i = 0; i < Nchannel-1; i++) {
        inputFile >> ichannel;
        inputFile >> ped0[i];
        inputFile >> sig0[i];
        inputFile >> ped1[i];
        inputFile >> sig1[i];
        fThreshold[i] = fNsigma*sig1[i];        
    }
    
    inputFile.close();    
}
#endif