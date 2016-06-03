/**************************************************************************
 * Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
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

// Function definitions needed for calculating the 3D near side acceptance correction

#include "AliJAcceptanceFunctions.h"
#include <TMath.h>
#include <TF1.h>
#include "Math/WrappedTF1.h"
#include "Math/BrentRootFinder.h"


/*
 * Default constructor
 */
AliJAcceptanceFunctions::AliJAcceptanceFunctions()
{
  // constructor
}

/*
 * Function for calculating two dimensional acceptance correction for 3D near side
 *
 *  double x[0] = deltaEta for trigger-associated pair
 *  double x[1] = deltaPhi for trigger-associated pair
 *  double par[0] = maximum eta range for the acceptance
 */
double AliJAcceptanceFunctions::AcceptanceCorrection3DNearSide(double *x, double *par){
  double deltaEta = x[0];
  double deltaPhi = par[0];
  double etaRange = par[1];
  
  AliJAcceptanceFunctions *functionWrapper = new AliJAcceptanceFunctions();
  TF1 distance("boundaryDistance",functionWrapper,&AliJAcceptanceFunctions::SideBoundaryDistance,-3*etaRange,3*etaRange,2,"AliJAcceptanceFunctions","SideBoundaryDistance");
  distance.SetParameter(0,deltaPhi);
  distance.SetParameter(1,deltaEta);
  ROOT::Math::WrappedTF1 wf1(distance);
  
  // Create numerical root finder
  ROOT::Math::BrentRootFinder brf;
  
  // Set parameters for the root finder
  brf.SetFunction( wf1, -3*etaRange, 3*etaRange );
  
  // Find the points where the distance between the constant deltaEta line and near-away side
  // boundary is zero.
  brf.Solve();
  double solution = brf.Root();
  
  // Calculate the lenght of the deltaEta curve inside and outside of the acceptance
  double diameterLegth = sqrt(pow(2*etaRange,2)+pow(2*etaRange,2));
  double outsideAcceptance = sqrt(pow(deltaEta,2)+pow(deltaEta,2));
  
  // Variable for the length of deltaEta curve inside the acceptance which will be on the away side
  double awaySideLength = 0;
  
  // Check if solution exists (if the root finder does not find a root, it just returns something)
  if(TMath::Abs(distance(solution)) < 1e-6){
    
    // Unless etaTrigger = -etaAssociated, there will be two solutions for the acceptance equation
    // First we need to find the second solution. This is symmetric with respect to etaTrigger = -etaAssociated line
    double firstEtaAssociated = solution;
    double firstEtaTrigger = deltaEta + firstEtaAssociated;
    
    double secondEtaTrigger = -firstEtaAssociated;
    double secondEtaAssociated = -firstEtaTrigger;
    
    // Put the solutions in correct order in eta
    double smallEtaTrigger = firstEtaTrigger;
    double bigEtaTrigger = secondEtaTrigger;
    double smallEtaAssoc = firstEtaAssociated;
    double bigEtaAssoc = secondEtaAssociated;
    
    if(firstEtaTrigger > secondEtaTrigger){
      smallEtaTrigger = secondEtaTrigger;
      bigEtaTrigger = firstEtaTrigger;
    }
    
    if(firstEtaAssociated > secondEtaAssociated){
      smallEtaAssoc = secondEtaAssociated;
      bigEtaAssoc = firstEtaAssociated;
    }
    
    // Now that we know the coordinates of the points where the away side is reached, we need to calculate which
    // portion of this is inside the acceptance and which is inside ideal acceptance.
    awaySideLength = sqrt(pow(bigEtaTrigger-smallEtaTrigger,2)+pow(bigEtaAssoc-smallEtaAssoc,2));
    
    // If the away side region starts outside of the acceptance, by symmetry we cannot ever reach the region inside acceptance.
    if(smallEtaAssoc < -etaRange || smallEtaTrigger < -etaRange){ // Region starts outside of acceptance
      return 0;
    }
    
  }
  
  // Delete the function wrapped created with new-keyword
  delete functionWrapper;
  
  // Remove the away side part from the diameter length
  diameterLegth = diameterLegth - awaySideLength;
  
  // Return the fraction of the length of the deltaEta curve that is inside of the acceptance
  return (diameterLegth - outsideAcceptance)/diameterLegth;
}

/*
 * Function for calculating distance between near-away side boundary and constant deltaEta lines in
 * (eta_trigger,eta_associated) plane as a function of eta_associated
 *
 *  double x[0] = associated particle pseudorapidity
 *  double par[0] = deltaPhi for trigger-associated pair
 *  double par[1] = deltaEta for trigger-associated pair
 */
double AliJAcceptanceFunctions::SideBoundaryDistance(double *x, double *par){
  double eta = x[0];
  double deltaPhi = par[0];
  double deltaEta = par[1];
  
  // First, calculate theta angle for trigger particle in near-away side boundary
  double angle = atan(-cos(2*atan(exp(-eta)))/(sin(2*atan(exp(-eta)))*cos(deltaPhi)));
  
  // Transform the angle to the interval [0,Pi]
  if(angle < 0) angle = TMath::Pi() + angle;
  
  // Transform the angle to pseudorapidity and check the distance from constant deltaEta line
  double edgeWithEta = -log(tan(angle/2))-deltaEta-eta;
  
  // Return the answer
  return edgeWithEta;
}
