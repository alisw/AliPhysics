
//THIS Also includes summing ALL cells in the jetcone towards the jet energy NOT just those above threshold!!!!!


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
Revision 1.2  2006/06/02 23:03:03  pavlinov
ALICE numbering scheme

Revision 1.1  2006/02/28 21:56:34  jklay
moving jetfinder code to subdirectory

Revision 1.16  2006/02/15 15:11:05  pavlinov
update of Jenn and Marco

Revision 1.15  2004/11/22 19:52:05  mhorner
Make sure pi0 get through

Revision 1.14  2004/04/02 17:11:33  mhorner
Marco's bug - fixes implemented

Revision 1.13  2004/03/26 01:00:38  mhorner
Implementing Marco's optimisations

Revision 1.12  2004/03/15 19:59:37  mhorner
Pyhtia comparison extended

Revision 1.11  2004/03/09 17:06:38  mhorner
Made more robust

Revision 1.10  2004/03/06 01:56:09  mhorner
Pythai comp code

Revision 1.9  2004/02/23 20:37:32  mhorner
changed geometry call

Revision 1.8  2004/01/29 23:28:44  mhorner
Jet Finder - hard coded geom parameters removed

Revision 1.7  2004/01/21 22:27:47  mhorner
Cleaning up coding conventions

Revision 1.6  2003/10/28 13:54:30  schutz
Compilation warnings fixed

Revision 1.5  2003/09/23 13:31:41  mhorner
Changed coordinate system

Revision 1.4  2003/09/19 13:16:20  mhorner
Added additional jet energy info


Revision 1.3  2003/09/04 12:49:56  mhorner
Changed hadron correction and added saving EMCAL and track contributions

*/


//*--Author: Sarah Blyth (LBL)
//*--Based on UA1 jet algorithm from LUND JETSET called from EMC-erj

#include "TTask.h"
#include "TMath.h"
#include "TParticle.h"
#include "AliRun.h"
#include "AliRunLoader.h"
#include "AliEMCALJetFinderInput.h"
#include "AliEMCALJetFinderOutput.h"
#include "AliEMCALJetFinderAlgo.h"
#include "AliEMCALJetFinderAlgoOmni.h"
#include "AliEMCALJetFinderAlgoUA1Unit.h"
#include "AliEMCALGeometry.h"
#include "AliEMCALLoader.h"
#include "AliEMCAL.h"
#include "AliEMCALDigit.h"
#include "AliEMCALJet.h"


ClassImp(AliEMCALJetFinderAlgoOmni) 

AliEMCALJetFinderAlgoOmni::AliEMCALJetFinderAlgoOmni() :
  fUnit(0),fUnitNoCuts(0),fHadCorr(0),fBGType(kRatio),fNumIter(0),
  fNumUnits(0),fESeed(0.),fConeRad(0.),fJetEMin(0.),fEtMin(0.),
  fMinMove(0.),fMaxMove(0.),fBGMaxMove(0.),fPtCut(0.),fBGPar(0.),
  fEBGTotal(0.),fEBGTotalOld(0.),fEBGAve(0.),fEnergy(0.),fJetEta(0.),
  fJetPhi(0.),fEtaInit(0.),fPhiInit(0.),fEtaB(0.),fPhiB(0.),fJetESum(0.),
  fJetEtaSum(0.),fJetPhiSum(0.),fDEta(0.),fDPhi(0.),fDistP(0.),fDistI(0.),
  fTempE(0.),fRad(0.),fNumInCone(0),fNumJets(0),fArrayInitialised(0)
{
 //Default constructor
if (fDebug>0) Info("AliEMCALJetFinderAlgoOmni","Beginning Default Constructor");
  AliEMCALGeometry * geom = AliEMCALGeometry::GetInstance("EMCAL_55_25","EMCAL");
  fNumIter           = 0;
  fNumUnits          = geom->GetNTowers();     //Number of towers in EMCAL
  fESeed             = 5.0;       //Default value
  fConeRad           = 0.3;       //Default value
  fJetEMin           = 10.0;      //Default value
  fEtMin             = 0.0;      //Default value
  fMinMove           = 0.05;      //From original UA1 JetFinder
  fMaxMove           = 0.15;      //From original UA1 JetFinder
  fBGMaxMove         = 0.035;     //From original UA1 JetFinder
  fPtCut             = 0;         
  fHadCorr           = 0;       
  fEBGTotal          = 1.0;       //Set to 1 so that no div by zero in first FindJets() loop
  fEBGTotalOld       = 0.0;
  fEBGAve            = 0.0;
  fEnergy            = 0.0;
  fJetEta            = 0.0;
  fJetPhi            = 0.0;
  fEtaInit           = 0.0;
  fPhiInit           = 0.0;
  fEtaB              = 0.0;
  fPhiB              = 0.0;
  fJetESum           = 0.0;
  fJetEtaSum         = 0.0;
  fJetPhiSum         = 0.0;
  fDEta              = 0.0;
  fDPhi              = 0.0;
  fDistP             = 0.0;
  fDistI             = 0.0;
  fTempE             = 0.0;
  fRad               = 2.0;      //Set to 2 to start 
  fNumInCone         = 0;
  fNumJets           = 0;
  fArrayInitialised  = 0;        //Set to FALSE to start
  fBGType            = kRatio;   //Set Ratio method as default BG subtraction method 
  fBGPar             = -1.0;      //Set to 1 to start
}


AliEMCALJetFinderAlgoOmni::AliEMCALJetFinderAlgoOmni(const AliEMCALJetFinderAlgoOmni& jfao) 
  : AliEMCALJetFinderAlgo(jfao),
    fUnit(jfao.fUnit),fUnitNoCuts(jfao.fUnitNoCuts),fHadCorr(jfao.fHadCorr),
    fBGType(jfao.fBGType),fNumIter(jfao.fNumIter),fNumUnits(jfao.fNumUnits),
    fESeed(jfao.fESeed),fConeRad(jfao.fConeRad),fJetEMin(jfao.fJetEMin),
    fEtMin(jfao.fEtMin), fMinMove(jfao.fMinMove),fMaxMove(jfao.fMaxMove),
    fBGMaxMove(jfao.fBGMaxMove),fPtCut(jfao.fPtCut),fBGPar(jfao.fBGPar),
    fEBGTotal(jfao.fEBGTotal),fEBGTotalOld(jfao.fEBGTotalOld),fEBGAve(jfao.fEBGAve),
    fEnergy(jfao.fEnergy),fJetEta(jfao.fJetEta), fJetPhi(jfao.fJetPhi),fEtaInit(jfao.fEtaInit),
    fPhiInit(jfao.fPhiInit),fEtaB(jfao.fEtaB),fPhiB(jfao.fPhiB),fJetESum(jfao.fJetESum),
    fJetEtaSum(jfao.fJetEtaSum),fJetPhiSum(jfao.fJetPhiSum),fDEta(jfao.fDEta),fDPhi(jfao.fDPhi),
    fDistP(jfao.fDistP),fDistI(jfao.fDistI),fTempE(jfao.fTempE),fRad(jfao.fRad),
    fNumInCone(jfao.fNumInCone),fNumJets(jfao.fNumJets),fArrayInitialised(jfao.fArrayInitialised)
{
  //copy ctor
}

 AliEMCALJetFinderAlgoOmni::~AliEMCALJetFinderAlgoOmni()
   {
     //Destructor
     if (fDebug>0) Info("AliEMCALJetFinderAlgoOmni","Beginning Destructor");
     delete[] fUnit;
     delete[] fUnitNoCuts;
   }

 void AliEMCALJetFinderAlgoOmni::SetJetFindingParameters
                               (Int_t numUnits, Float_t eSeed, Float_t coneRad, Float_t jetEMin, Float_t etMin, 
                               Float_t minMove, Float_t maxMove, Float_t bgMaxMove)
   {
     //Sets parameters for the JetFinding algorithm
     if (fDebug>1) Info("SetJetFindingParameters","Setting parameters for JetFinding");

     SetNumUnits(numUnits);
     SetJetESeed(eSeed);
     SetConeRad(coneRad);
     SetJetEMin(jetEMin);
     SetEtMin(etMin);
     SetMinMove(minMove);
     SetMaxMove(maxMove);
     SetBGMaxMove(bgMaxMove);
   }

 void AliEMCALJetFinderAlgoOmni::SetJetFindingParameters
                               (Int_t numUnits, Float_t eSeed, Float_t coneRad, Float_t jetEMin, Float_t etMin)
   {
     //Sets fewer parameters for the JetFinding algorithm
     if (fDebug>1) Info("SetJetFindingParameters","Setting parameters for JetFinding");

     SetNumUnits(numUnits);
     SetJetESeed(eSeed);
     SetConeRad(coneRad);
     SetJetEMin(jetEMin);
     SetEtMin(etMin);
     SetMinMove(fMinMove);
     SetMaxMove(fMaxMove);
     SetBGMaxMove(fBGMaxMove);
   }

 void AliEMCALJetFinderAlgoOmni::InitUnitArray()
   {
     //Initialises unit arrays
     if(fArrayInitialised) delete [] fUnit;
     if(fArrayInitialised) delete [] fUnitNoCuts;
     fUnit = new AliEMCALJetFinderAlgoUA1Unit[fNumUnits];
     fUnitNoCuts = new AliEMCALJetFinderAlgoUA1Unit[fNumUnits];
     fArrayInitialised = 1;
   }

 void AliEMCALJetFinderAlgoOmni::FillUnitArray(AliEMCALJetFinderAlgoUA1FillUnitFlagType_t flag)
   {
	   // Fill the unit array 
     if (fDebug>1) Info("FillUnitArray","Beginning FillUnitArray");
         //    AliEMCAL* pEMCAL = (AliEMCAL*) gAlice->GetModule("EMCAL");
         //   if (pEMCAL){ 
         //	     AliEMCALGeometry* geom =  AliEMCALGeometry::GetInstance(pEMCAL->GetTitle(), "");
         //     }else
         //    {

     AliEMCALGeometry * geom = AliEMCALGeometry::GetInstance();
     if (geom == 0)
      geom = AliEMCALGeometry::GetInstance("EMCAL_55_25","EMCAL");

        //    }
         
     AliEMCALJetFinderAlgoUA1FillUnitFlagType_t option = flag;
     Int_t         numTracks, numDigits;
    
     //Loops over all elements in the AliEMCALJetFinderAlgoUA1Unit array and 
     //fills the objects with relevant values from the Data Input object
     if (fDebug>10) Info("FillUnitArray","Filling array with Unit objects");
     if (fDebug>15) Info("FillUnitArray","NTracks %i NDigits %i",fInputPointer->GetNTracks(),fInputPointer->GetNDigits());
	 numTracks = fInputPointer->GetNTracks();
	 numDigits = fInputPointer->GetNDigits();
	 TParticle         *myPart;
	 AliEMCALDigit     *myDigit;
     if (!fPythiaComparison)
     {
	 //Fill units with Track info if appropriate
	 if(option==kFillTracksOnly || option ==kFillAll) 
	   {          
	    for(Int_t j=0; j<numTracks; j++)
	    {
	     myPart = fInputPointer->GetTrack(j);
	     Float_t eta = myPart->Eta();
	     //   have to be tune for TRD1; May 31,06 
	     // Float_t  phi = myPart->Phi();
	     // Int_t towerID = geom->TowerIndexFromEtaPhi(eta,180.0/TMath::Pi()*phi);
	     Int_t towerID = 0; 

	     Float_t  pT = myPart->Pt();
	     Float_t unitEnergy = fUnit[towerID-1].GetUnitEnergy(); 
	     Float_t unitEnergyNoCuts = fUnitNoCuts[towerID-1].GetUnitEnergy();

	     
	     //OLD WAY:   //Do Hadron Correction
              if(fHadCorr != 0)
	       {
		 Double_t   fullP = myPart->P();
		 Double_t   hCEnergy = fHadCorr->GetEnergy(fullP, (Double_t)eta);
		 unitEnergy -= hCEnergy*TMath::Sin(myPart->Theta());
		 unitEnergyNoCuts -= hCEnergy*TMath::Sin(myPart->Theta());
		 fUnit[towerID-1].SetUnitEnergy(unitEnergy);
		 fUnitNoCuts[towerID-1].SetUnitEnergy(unitEnergyNoCuts);
	       } //end Hadron Correction loop 
	      
    
	     /*
	      //Do Hadron Correction with propagate phi for the track
	      if(fHadCorr != 0)
		{
		  Bool_t curl = 1;
		  Float_t deltaPhi;
		  TParticlePDG *pdg = myPart->GetPDG();
		  if(pdg->Charge() < 0)
		    {
		      deltaPhi = PropagatePhi(myPart->Pt(), -1.0, curl);
		    }
		  else{
		    deltaPhi = PropagatePhi(myPart->Pt(), 1.0, curl);
		  }
		  phi += deltaPhi;
		  //Get new tower id for cell that track would curve into
                  Int_t towerID2;
                  if(phi>(TMath::Pi()/180.0)*geom->GetArm1PhiMax() || phi<(TMath::Pi()/180.0)*geom->GetArm1PhiMin())
                    {
                      towerID2 = -1;
                    }
                  else{
                      towerID2 = geom->TowerIndexFromEtaPhi(eta,180.0/TMath::Pi()*phi);
                    }
                  
		  if(towerID2 != -1)
		    {
		      //Find unit energy of new tower
		      Float_t unitEnergy2 = fUnit[towerID2-1].GetUnitEnergy();
		      Float_t unitEnergy2NoCuts = fUnitNoCuts[towerID2-1].GetUnitEnergy();
		      Double_t   fullP = myPart->P();
		      Double_t   hCEnergy = fHadCorr->GetEnergy(fullP, (Double_t)eta);
		      unitEnergy2 -= hCEnergy*TMath::Sin(myPart->Theta());
		      unitEnergy2NoCuts -= hCEnergy*TMath::Sin(myPart->Theta());
		      fUnit[towerID2-1].SetUnitEnergy(unitEnergy2);
		      fUnitNoCuts[towerID2-1].SetUnitEnergy(unitEnergy2NoCuts);
		    }//end if for towerID2
		}//end Hadron Correction loop
	     */

	      fUnitNoCuts[towerID-1].SetUnitEnergy(unitEnergyNoCuts + pT);
	     //Do Pt cut on tracks
	     if(fPtCut != 0 && pT < fPtCut) continue;

	     fUnit[towerID-1].SetUnitEnergy(unitEnergy+pT);

             }//end tracks loop
	   }//end Tracks condition


	 //Fill units with Digit info if appropriate
	 if(option ==kFillDigitsOnly || option ==kFillAll)
	   {
            for(Int_t k=0; k<numDigits; k++)
	    {
	     myDigit = fInputPointer->GetDigit(k);
             if (fDebug>10) Info("FillUnitArray","getting digits %i %i numdigits",k,numDigits );
	     Int_t towerID = myDigit->GetId();
	     Int_t amplitude = myDigit->GetAmp();     //Gets the integer valued amplitude of the digit
	     Float_t amp = (Float_t)amplitude;        //Need to typecast to Float_t before doing real energy conversion
	     Float_t digitEnergy = amp/10000000.0;    //Factor of 10 million needed to convert!
	     Float_t unitEnergy = fUnit[towerID-1].GetUnitEnergy() + digitEnergy;
	     Float_t unitEnergyNoCuts = fUnitNoCuts[towerID-1].GetUnitEnergy() + digitEnergy;
	     fUnit[towerID-1].SetUnitEnergy(unitEnergy);
	     fUnitNoCuts[towerID-1].SetUnitEnergy(unitEnergyNoCuts);
	    }//end digits loop
	   }//end digits condition

	 //Set all unit flags, Eta, Phi
	 for(Int_t i=0; i<fNumUnits; i++)
	   {
             if (fDebug>10) Info("FillUnitArray","Setting all units outside jets");
	     //Set all units to be outside a jet initially
	     fUnit[i].SetUnitFlag(kOutJet);           
	     fUnit[i].SetUnitID(i+1);
	     Float_t eta;
	     Float_t phi;
	     geom->EtaPhiFromIndex(fUnit[i].GetUnitID(), eta, phi);
	     fUnit[i].SetUnitEta(eta);
	     fUnit[i].SetUnitPhi(phi*TMath::Pi()/180.0);
	     //Set all units to be outside a jet initially
	     fUnitNoCuts[i].SetUnitFlag(kOutJet);          
	     fUnitNoCuts[i].SetUnitID(i+1);
	     eta = 0.0;
	     phi = 0.0;
	     geom->EtaPhiFromIndex(fUnitNoCuts[i].GetUnitID(), eta, phi);
	     fUnitNoCuts[i].SetUnitEta(eta);
	     fUnitNoCuts[i].SetUnitPhi(phi*TMath::Pi()/180.0);
	     //	     if(i>13000) cout<<"!!!!!!!!!!!!!!!!!For unit0, eta="<<eta<<" and phi="<<phi*TMath::Pi()/180.0<<" and ID="<<fUnit[i].GetUnitID()<<endl;
	     //  if(fUnit[i].GetUnitEnergy()>0) cout<<"Unit ID "<<fUnit[i].GetUnitID() <<"with eta="<<eta<<" and phi="<<phi*TMath::Pi()/180.0<<" has energy="<<fUnit[i].GetUnitEnergy()<<endl;
	   }//end loop over all units in array (same as all towers in EMCAL)

     }
     if (fPythiaComparison)
     {
	     for(Int_t j=0; j<numTracks; j++)		                 
	     {
		     myPart = fInputPointer->GetTrack(j);
		     fUnit[j].SetUnitID(j);
		     fUnit[j].SetUnitFlag(kOutJet);
		     fUnit[j].SetUnitEta(myPart->Eta());
		     fUnit[j].SetUnitPhi(myPart->Phi());
		     if (myPart->Energy()*TMath::Sin(myPart->Theta()) > fPtCut || myPart->GetPDG()->Charge() == 0.0  )
		     {
			     fUnit[j].SetUnitEnergy(myPart->Energy()*TMath::Sin(myPart->Theta()));
		     }else
		     {
			     fUnit[j].SetUnitEnergy(0.00);
		     }
		     fUnitNoCuts[j].SetUnitID(j);
		     fUnitNoCuts[j].SetUnitFlag(kOutJet);
		     fUnitNoCuts[j].SetUnitEta(myPart->Eta());
		     fUnitNoCuts[j].SetUnitPhi(myPart->Phi());
		     fUnitNoCuts[j].SetUnitEnergy(myPart->Energy()*TMath::Sin(myPart->Theta()));
	     }
	     for(Int_t k=numTracks; k < fNumUnits; k++)
	     {
		     fUnit[k].SetUnitID(k);
		     fUnit[k].SetUnitFlag(kBelowMinEt);
		     fUnit[k].SetUnitEta(0.0);
		     fUnit[k].SetUnitPhi(0.0);
		     fUnit[k].SetUnitEnergy(0.0);		     
		     fUnitNoCuts[k].SetUnitID(k);
		     fUnitNoCuts[k].SetUnitFlag(kBelowMinEt);
		     fUnitNoCuts[k].SetUnitEta(0.0);
		     fUnitNoCuts[k].SetUnitPhi(0.0);
		     fUnitNoCuts[k].SetUnitEnergy(0.0);		          
	     }
     }



   }


 void AliEMCALJetFinderAlgoOmni::Sort(AliEMCALJetFinderAlgoUA1Unit *unit, Int_t length)
 {
   //Calls the recursive quicksort method to sort unit objects in decending order of Energy
   if (fDebug>1) Info("Sort","Sorting Unit objects");
   QS(unit, 0, length-1);
 }
  

 void AliEMCALJetFinderAlgoOmni::QS(AliEMCALJetFinderAlgoUA1Unit *unit, Int_t left, Int_t right)
 {
  //Sorts the AliEMCALJetFinderAlgoUA1Unit objects in decending order of Energy
   if (fDebug>111) Info("QS","QuickSorting Unit objects");   

   Int_t    i;
   Int_t    j;
   AliEMCALJetFinderAlgoUA1Unit  unitFirst;
   AliEMCALJetFinderAlgoUA1Unit  unitSecond;

   i = left;
   j = right;
   unitFirst = unit[(left+right)/2];

 do
  {
    while( (unit[i].GetUnitEnergy() > unitFirst.GetUnitEnergy()) && (i < right)) i++;
    while( (unitFirst.GetUnitEnergy() > unit[j].GetUnitEnergy()) && (j > left)) j--;

    if(i <= j)
      {
	unitSecond = unit[i];
	unit[i] = unit[j];
	unit[j] = unitSecond;
	i++;
	j--;
      }//end if
  }while(i <= j);

 if(left < j) QS(unit, left, j);
 if(i < right) QS(unit, i, right);
 }


 void AliEMCALJetFinderAlgoOmni::FindBG()
 {
   if(fBGType == kRatio) RatioBG();
   else if(fBGType == kCone) ConeBG();
   else if(fBGType == kConstant) ConstantBG();
 }

 void AliEMCALJetFinderAlgoOmni::RatioBG()
   {
     //Finds the background energy for the iteration
     //using the Ratio method
     if (fDebug>1) Info("FindBG","Finding Average Background"); 
     //Store BGEperCell from previous iteration!
     fEBGTotalOld = fEBGTotal;
     fEBGTotal          = 0.0;
     Int_t numCone      = 0;

     //If user has not set fBGPar, set it to the default
     //for TPC = 90% efficiency, PtCut = 2GeV/c, timecut = 30ns
     if(fBGPar == -1) fBGPar = 0.4685;

     //Loop over all unit objects in the Unit array and link to same
     //unit ID in NoCuts Unit array
     for(Int_t i=0; i<fNumUnits; i++)
       {
	 if(fUnit[i].GetUnitFlag() != kInJet)
	   {
	     Int_t id = fUnit[i].GetUnitID();
	     fEBGTotal += fUnitNoCuts[id-1].GetUnitEnergy();  
	   }
	 else numCone++;
       }//end for

     fEBGTotal *= fBGPar;
     fEBGAve = fEBGTotal / (fNumUnits - numCone);
     if (fDebug>5) Info("FindBG","Average BG is %f: ",fEBGAve); 

     for(Int_t count=0; count<fNumUnits;count++)
       {
	 fUnit[count].SetUnitFlag(kOutJet);
       }//end for
   }

 void AliEMCALJetFinderAlgoOmni::ConeBG()
   {
     //Finds the background energy for the iteration
     //using all energy not contained inside a jet
     if (fDebug>1) Info("FindBG","Finding Average Background"); 
     //Store old value of BGEperCell!
     fEBGTotalOld = fEBGTotal;
     fEBGTotal          = 0.0;
     Int_t numCone      = 0;

     //Loop over all unit objects in the array and sum the energy of those not in a jet
     for(Int_t i=0; i<fNumUnits; i++)
       {
	 if(fUnit[i].GetUnitFlag() != kInJet)
	   fEBGTotal += fUnit[i].GetUnitEnergy();
	 else numCone++;
       }//end for

     fEBGAve = fEBGTotal / (fNumUnits - numCone);
     if (fDebug>5) Info("FindBG","Average BG is %f: ",fEBGAve);      

     for(Int_t count=0; count<fNumUnits;count++)
       {
	 fUnit[count].SetUnitFlag(kOutJet);
       }//end for
   }

 void AliEMCALJetFinderAlgoOmni::ConstantBG()
   {     
     //Finds the background energy for the iteration
     //using all energy not contained inside a jet
     if (fDebug>1) Info("FindBG","Finding Average Background"); 

     //If user has not set fBGPar, set it to the default
     //for TPC = 90% efficiency, PtCut = 2GeV/c, timecut = 30ns
     if(fBGPar == -1) fBGPar = 0.03378;

     fEBGAve = fBGPar;
     if (fDebug>5) Info("FindBG","Average BG is %f: ",fEBGAve);      

     fEBGTotal          = 0.0;
     Int_t numCone      = 0;
     for(Int_t count=0; count<fNumUnits;count++)
       {
	 if(fUnit[count].GetUnitFlag() == kInJet)
	   {
	     numCone++;
	   }
	 fUnit[count].SetUnitFlag(kOutJet);
       }//end for
     fEBGTotal = fEBGAve * (fNumUnits-numCone);
     fEBGTotalOld = fEBGTotal;
   }

 void AliEMCALJetFinderAlgoOmni::FindJetEtaPhi(Int_t counter)
   {
     //Finds the eta and phi of the jet axis
     if (fDebug>10) Info("FindJetEtaPhi","Finding Jet Eta and Phi");

     fDEta = fUnit[counter].GetUnitEta() - fEtaInit;
     fDPhi = fUnit[counter].GetUnitPhi() - fPhiInit;

     fEnergy = fUnit[counter].GetUnitEnergy() - fEBGAve;
     fJetEtaSum += fEnergy * fDEta;
     fJetPhiSum += fEnergy * fDPhi;
     fJetESum += fEnergy;
     if (fJetESum >1.0e-7)
     {
	     fJetEta = fEtaInit + (fJetEtaSum / fJetESum);
	     fJetPhi = fPhiInit + (fJetPhiSum / fJetESum);
     }
   }


 void AliEMCALJetFinderAlgoOmni::FindJetEnergy()
   {
     //Finds the energy of the jet after the final axis has been found
     if (fDebug>1) Info("FindJetEnergy","Finding Jet Energy");

     for(Int_t i=0; i<fNumUnits; i++)
       {
	 //Loop over all unit objects in the array and find if within cone radius
	 Float_t dEta = fUnit[i].GetUnitEta() - fJetEta;
	 Float_t dPhi = fUnit[i].GetUnitPhi() - fJetPhi;
	 Float_t rad;
	 if ((dEta*dEta) + (dPhi*dPhi)>1.e-7)
	 {
 		 rad = TMath::Sqrt( (dEta*dEta) + (dPhi*dPhi) );
	 }else
	 {
		 rad=0.0;
	 }

	 if(fUnit[i].GetUnitFlag()==kOutJet && rad<= fConeRad)
	   {
	     fUnit[i].SetUnitFlag(kInCurrentJet);
	     Float_t energy = fUnit[i].GetUnitEnergy() - fEBGAve;
	     fJetESum += energy;                             
	     fJetEtaSum += energy * dEta;
	     fJetPhiSum += energy * dPhi;
	     fNumInCone++;                     //Increment the number of cells in the jet cone
	   }//end if
       }//end for
   }


 void AliEMCALJetFinderAlgoOmni::StoreJetInfo()
   {
     //Stores the resulting jet information in appropriate storage structure (TO BE DECIDED!!!!)
     if (fDebug>1) Info("StoreJetInfo","Storing Jet Information");
     AliEMCALGeometry * geom = AliEMCALGeometry::GetInstance();
     if (geom == 0)
       geom = AliEMCALGeometry::GetInstance("EMCAL_55_25","EMCAL");
     //Store:
     //fJetESum is the final jet energy (background has been subtracted)
     //fJetEta is the final jet Eta
     //fJetPhi is the final jet Phi
     //fNumInCone is the final number of cells included in the jet cone
     //fEtaInit is the eta of the initiator cell
     //fPhiInit is the phi of the initiator cell

     AliEMCALJet jet(fJetESum,fJetPhi,fJetEta);

      cout<<"For iteration "<<fNumIter <<" and Jet number " <<fNumJets <<endl;
      cout<<"The jet energy is: " <<fJetESum <<endl;
      cout<<"The jet eta is ---->" <<fJetEta <<endl;
      cout<<"The jet phi is ---->" <<fJetPhi <<endl;

     Int_t             numberTracks = fInputPointer->GetNTracks();
     Int_t  	       numberDigits = fInputPointer->GetNDigits();
     AliEMCALDigit     *myD;
     TParticle         *myP;
     Int_t             numTracksInCone = 0;
     Float_t           trackEnergy = 0.0;
     Float_t           trackEnergyPtCut =0.0; 
     Float_t           emcalEnergy = 0.0;
     Float_t           emcalEnergyBGSub = 0.0;

     for(Int_t counter=0; counter<numberTracks; counter++)
       {
	myP = fInputPointer->GetTrack(counter);
	Float_t eta = myP->Eta();
	Float_t phi = myP->Phi(); 
	Float_t deta = fJetEta-eta;
	Float_t dphi = fJetPhi -phi;
	Float_t rad = TMath::Sqrt( (deta*deta) + (dphi*dphi));
	if(rad<=fConeRad) numTracksInCone++;
       }//end for

     Float_t    *pTArray = new Float_t[numTracksInCone];
     Float_t    *etaArray = new Float_t[numTracksInCone];
     Float_t    *phiArray = new Float_t[numTracksInCone];
     Int_t      *pdgArray = new Int_t[numTracksInCone];
     Int_t             index = 0;

     for(Int_t counter2=0; counter2<numberTracks; counter2++)
       {
	 myP = fInputPointer->GetTrack(counter2);
 	 Float_t eta = myP->Eta();
         Float_t phi = myP->Phi(); 
	 Float_t deta = fJetEta-eta;
	 Float_t dphi = fJetPhi -phi;
	 Float_t rad ;
	 if ((deta*deta) + (dphi*dphi)>1.e-7)
	 {
 		 rad = TMath::Sqrt( (deta*deta) + (dphi*dphi) );
	 }else
	 {
		 rad=0.0;
	 }

	 if(rad<=fConeRad)
	   {
	     pTArray[index] = myP->Pt();
	     //Calculate track contribution within jetcone
	     trackEnergy += myP->Pt();
	     if(myP->Pt() >= fPtCut) trackEnergyPtCut += myP->Pt();
	     etaArray[index] = eta;
	     phiArray[index] = phi;
	     pdgArray[index] = myP->GetPdgCode();
	     index++;
	     if(fHadCorr != 0)
	     {
		     Double_t   fullP = myP->P();
		     Double_t   hCEnergy = fHadCorr->GetEnergy(fullP, (Double_t)eta);
		     emcalEnergy -= hCEnergy*TMath::Sin(myP->Theta());
		     emcalEnergyBGSub -= hCEnergy*TMath::Sin(myP->Theta());
	     } //end Hadron Correction loop 
			   
	   }//end if
       }//end for

     //Loop over digits to find EMCal contribution within jetcone
     for(Int_t counter3=0; counter3<numberDigits; counter3++)
       {
	 myD = fInputPointer->GetDigit(counter3);
	 //GET DIGIT ETA, PHI so that can check if inside R!
	 Float_t eta = 0.0;
	 Float_t phi = 0.0;
	 Int_t iID = myD->GetId();
	 geom->EtaPhiFromIndex(iID, eta, phi);
	 Float_t deta = fJetEta-eta;
	 Float_t dphi = fJetPhi -(TMath::Pi()/180.0)*phi;
	 //Float_t rad = TMath::Sqrt( (deta*deta) + (dphi*dphi));
	 Float_t rad ;
	 if ((deta*deta) + (dphi*dphi)>1.e-7)
	 {
 		 rad = TMath::Sqrt( (deta*deta) + (dphi*dphi) );
	 }else
	 {
		 rad=0.0;
	 }

	 if(rad<=fConeRad)
	   {
	 Int_t amplitude = myD->GetAmp();     //Gets the integer valued amplitude of the digit
	 Float_t amp = (Float_t)amplitude;        //Need to typecast to Float_t before doing real energy conversion
	 Float_t digitEnergy = amp/10000000.0;    //Factor of 10 million needed to convert!
	 emcalEnergy += digitEnergy;
	 emcalEnergyBGSub += (digitEnergy - fEBGAve);
	   }//end if
       }//end count3 for

     //Save in JET object
     jet.SetTrackList(numTracksInCone,pTArray, etaArray, phiArray, pdgArray);
     jet.SetEMCALEnergy(emcalEnergy);
     jet.SetEMCALEnergyBGSub(emcalEnergyBGSub);
     jet.SetTrackEnergy(trackEnergy);
     jet.SetTrackEnergyPtCut(trackEnergyPtCut);
     fOutputPointer->AddJet(&jet);
     delete[] pTArray;
     delete[] etaArray;
     delete[] phiArray;
     delete[] pdgArray;
   }


 void AliEMCALJetFinderAlgoOmni::FindJets()
   {
     //Runs the complete UA1 JetFinding algorithm to find jets!
     if (fDebug>1) Info("FindJets","Starting Jet Finding!!!");

     //If the array of JetFinderUnit objects has not been initialised then initialise with default settings
     if(!fArrayInitialised) 
      {
       InitUnitArray();
       FillUnitArray(kFillAll);
      }//end if
     if (fDebug>1) Info("FindJets","Unit array filled");

     //Step 1. Sort the array in descending order of Energy
     Sort(fUnit,fNumUnits);

     //Step 2. Set the number of iterations and Number of jets found to zero to start
     fNumIter = 0;
     fNumJets = 0;

     //Step 3. Begin the iteration loop to find jets
     //Need to iterate the algorithm while number of iterations<2 OR number of iterations<10 AND 
     //the value of the average background has changed more than specified amount
     //Min iterations = 2, Max iterations = 10
     //while(fNumIter<2 || (fNumIter <10 && ( (fEBGTotal-fEBGTotalOld)/fEBGTotal) > fBGMaxMove) )

     while(fNumIter<2 || (fNumIter <10 && ( fEBGTotal-fEBGTotalOld) > fEBGTotal*fBGMaxMove) )
       {
        if (fDebug>1) Info("FindJets","Starting BIG iteration ---> %i",fNumIter);

         //Step 4. Find the value of the average background energy
	 FindBG();
	 fOutputPointer->Reset(kResetJets); //Reset output object to store info for new iteration
	 fNumJets=0;

	 //Loop over the array of unit objects and flag those with energy below MinCellEt
         Int_t numbelow = 0;
	 for(Int_t j=0; j<fNumUnits; j++)
	   {
	     if( (fUnit[j].GetUnitEnergy()-fEBGAve) < fEtMin)
        	{       
		  //          fUnit[j].SetUnitFlag(kBelowMinEt);    TAKING OUT kBelow flag
                  numbelow++;
                }//end if
	   }//end for
	 //cout<<"THERE WERE "<<numbelow<<" units with E <EtMin!!!!!!!!!!!!!!!"<<endl;

	 //Do quick check if there are no jets upfront
	 // if(fUnit[0].GetUnitFlag() == kBelowMinEt)
	 if( (fUnit[0].GetUnitEnergy()-fEBGAve) < fEtMin)
	   {
            cout <<"There are no jets for this event!" <<endl;
	    break;
	   }//end if

         //Step 5. Begin with the first jet candidate cell (JET SEED LOOP)
         if (fDebug>5) Info("FindJets","Beginning JET SEED LOOP");
	 for(Int_t count=0; count<fNumUnits; count++)
	   {

//CHECK CONDITION HERE _ NOT SURE IF SHOULD MAYBE BE: GetUnitEnergy()-fEBGAve >fESeed?????????????????????????????
	     if(fUnit[count].GetUnitEnergy()>=fESeed && fUnit[count].GetUnitFlag()==kOutJet)
	       {
		 fEnergy = fUnit[count].GetUnitEnergy() - fEBGAve;
		 fJetEta = fUnit[count].GetUnitEta();
		 fJetPhi = fUnit[count].GetUnitPhi();
		 Int_t seedID = fUnit[count].GetUnitID();
                 if (fDebug>5) Info("FindJets","Inside first candidate jet seed loop for time : %i", count);
                 if (fDebug>5) Info("FindJets","Found candidate energy %f ",fEnergy);
                 if (fDebug>5) Info("FindJets","Found candidate eta %f ", fJetEta);
                 if (fDebug>5) Info("FindJets","Found candidate phi %f ", fJetPhi);
		 if (fDebug>5) Info("FindJets","Found candidate ID %i", seedID);

		 fEtaInit = fJetEta;
		 fPhiInit = fJetPhi;
		 fEtaB = fJetEta;
		 fPhiB = fJetPhi;
		 Int_t testflag = 1;
		 do
		   {
		 fJetESum = 0.0;
		 fJetEtaSum = 0.0;
		 fJetPhiSum = 0.0;
       
         //Step 6. Find Jet Eta and Phi
		 //Loop over all units in the array to find the ones in the jet cone and determine contrib to Jet eta, phi
		     for(Int_t count1=0; count1<fNumUnits; count1++)		   
		       {
			 if(fUnit[count1].GetUnitID() == seedID && testflag)
			 {
				 testflag=0;
				 continue;   //skip unit if the jetseed to avoid doublecounting
			 }
			 if(fUnit[count1].GetUnitFlag() == kOutJet)
			   {
			     fDEta = fUnit[count1].GetUnitEta() - fJetEta;
			     fDPhi = fUnit[count1].GetUnitPhi() - fJetPhi;
			     if ( (fDEta*fDEta) + (fDPhi*fDPhi) >1.e-7)
			     {
				     fRad = TMath::Sqrt( (fDEta*fDEta) + (fDPhi*fDPhi) );
			     }else
			     {
				     fRad=0.000;
			     }
			     if(fRad <= fConeRad)
			       {
				 FindJetEtaPhi(count1); 
			       }//end if
			   }//end if
		       }//end for (Jet Eta, Phi LOOP)
			     
		      //Find the distance cone centre moved from previous cone centre
		      if (fDebug>10) Info("FindJets","Checking if cone move small enough");
		      if (((fJetEta-fEtaB)*(fJetEta-fEtaB)) + ((fJetPhi-fPhiB)*(fJetPhi-fPhiB))  >1.e-7)
		      {
			      fDistP = TMath::Sqrt( ((fJetEta-fEtaB)*(fJetEta-fEtaB)) + ((fJetPhi-fPhiB)*(fJetPhi-fPhiB)) );
		      }else
		      {
			      fDistP = 0.00;
		      }
		      //     if(fDistP <= fMinMove) break;
			     

		      //Find the distance cone centre is from initiator cell
		      if (fDebug>10) Info("FindJets","Checking if cone move too large");
		      if ( ((fJetEtaSum)*(fJetEtaSum))+((fJetPhiSum)*(fJetPhiSum)) >1.e-7)
		      {
			      fDistI = TMath::Sqrt( ((fJetEtaSum/fJetESum)*(fJetEtaSum/fJetESum)) + ((fJetPhiSum/fJetESum)*
												    (fJetPhiSum/fJetESum)));
		      }else
		      {
			      fDistI = 0.00;
		      }

		      if(fDistP>fMinMove && fDistI<fMaxMove)
			{
			  fEtaB = fJetEta;
			  fPhiB = fJetPhi;
			}//end if
		 
		   }while(fDistP>fMinMove && fDistI<fMaxMove);
			  
		 fJetEta = fEtaB;
		 fJetPhi = fPhiB;


       //Step 7. Find the Jet Energy
                 if (fDebug>1) Info("FindJets","Looking for Jet energy");
		 fJetESum = 0.0;
		 fJetEtaSum = 0.0;
		 fJetPhiSum = 0.0;
		 fNumInCone = 0;
		 FindJetEnergy();

		 //cout<<"Number of cells in jet cone is: "<<fNumInCone<<endl;

       //Step 8. Check if the jet is a valid jet
		 //Check if cluster energy is above Min allowed to be a jet
//DID NOT DO THE COSH COMPARISON HERE -> NEED TO CHECK WHICH COMPARISON IS BEST!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                 if (fDebug>5) Info("FindJets","Checking cluster is valid jet");
		 if(fJetESum < fJetEMin)
		   {
		     for(Int_t count2=0; count2<fNumUnits; count2++)
		       {
			 if(fUnit[count2].GetUnitFlag()==kInCurrentJet || fUnit[count2].GetUnitFlag()==kOutJet)
			   fUnit[count2].SetUnitFlag(kOutJet);
		       }//end for
                   if (fDebug>10) Info("FindJets","NOT a valid jet cell");
		  }else
		    {
		     for(Int_t count2=0; count2<fNumUnits; count2++)
		       {
			 if(fUnit[count2].GetUnitFlag()==kInCurrentJet)
			   {
			     //	     cout<<"Setting unit #"<<count2 <<" to be officially in a jet!"<<endl;
                           fUnit[count2].SetUnitFlag(kInJet);
			   }
		       }//end for			

 //NEED TO CHECK FINAL WEIRD ITERATION OF ETA AND PHI CHANGES!!!!!!!!!
		     //	 fJetPhi += fJetPhiSum/fJetESum;        //CHECK!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		     //  fJetEta += fJetEtaSum/fJetESum;        //CHECK!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		     fNumJets++;              //Incrementing number of jets found
		     StoreJetInfo();          //Storing jet info

		 }//end if (check cluster above Min Jet Energy)
	       }//end if (Jet Seed condition)
	   }//end (JET SEED LOOP)

if (fDebug>5) Info("FindJets","End of BIG iteration number %i",fNumIter);
// this->Dump();
	 fNumIter++;
       }//end 10 iteration WHILE LOOP
 }











