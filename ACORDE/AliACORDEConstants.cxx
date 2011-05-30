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

/* $Id$ */

////////////////////////////////////////////////////////////////////////
//
// AliACORDEConstants class
//
// This class includes the constants needed by ACORDE detector in 
// easily accessible place. All constants are public const static data 
// members. The class is never instatiated.
// Authors: Arturo Fernandez, Enrique Gamez, Mario Rodriguez Cahuantzi, Eleazar Cuautle(ICN-UNAM) 
//         FCFM-UAP, Mexico.
// Last update: Nov. 24th 08
////////////////////////////////////////////////////////////////////////

#include "AliACORDEConstants.h"

AliACORDEConstants* AliACORDEConstants::fgInstance = 0;

const Float_t AliACORDEConstants::fgkModuleLength          = 300.0;
const Float_t AliACORDEConstants::fgkModuleWidth           = 24.0;
const Float_t AliACORDEConstants::fgkModuleHeight          =  10.0;
const Float_t AliACORDEConstants::fgkPlasticLength = 190.0;
const Float_t AliACORDEConstants::fgkPlasticWidth  =  19.5;
const Float_t AliACORDEConstants::fgkPlasticHeight =   1.0;
const Float_t AliACORDEConstants::fgkProfileWidth =    3.8;
const Float_t AliACORDEConstants::fgkProfileThickness = 0.3;
const Float_t AliACORDEConstants::fgkDepth               =4420; 

const Float_t AliACORDEConstants::fgkHitEnergyThreshold = 1.52; // MeVs
const Float_t AliACORDEConstants::fgkMaxHitTimeDifference = 40.0; // nanoseconds
const Int_t AliACORDEConstants::fgkMultiMuonThreshold = 2;
const Float_t AliACORDEConstants::fgkMultiMuonWindow = 25;
const Float_t AliACORDEConstants::fgkInsideModulePositionX[60] ={
149.831, 743.687, 744.367, 744.4, 744.535, 745.991, 745.41,0,0, 151.197,
529.449, 529.76, 529.911, 529.911, 530.172,529.709, 529.692, 529.597, 528.859, 528.131,
304.946, 304.472, 304.092, 303.734, 303.165, 303.301, 303.195, 303.422, 303.927, 304.091,
-3.974, -3.806, -2.984, -2.855, -3.042, -3.124, -3.395, -2.774, -3.072, -2.897,
-319.384, -318.925, 0,-318.133, -317.404, -317.365, -316.973, -317.222,-317.564,-317.913,
149.892, -537.461, -537.75, -537.327,-536.877, -502.826, -506.978,-531.402,-530.587,149.541};
const Float_t AliACORDEConstants::fgkInsideModulePositionY[60] ={
860.235, 486.767, 486.763, 487.156, 487.018, 485.638, 486.394, 0,0,859.869,
700.202, 700.11, 700.345, 700.746, 701.481, 701.662, 701.925, 701.51, 701.64, 702.098,
859.937, 859.712, 859.738, 859.788, 859.88, 860.278, 860.155, 860.131, 860.14, 859.731,
860.096, 860.035, 860.416, 860.451, 860.655, 860.445, 860.601, 860.275, 860.623, 860.665,
916.198, 916.005, 0, 915.731, 915.768, 914.931, 914.708, 914.794, 915.021, 915.084,
860.287, 692.384, 692.392, 693.071, 692.86, 725.954, 722.077, 698.292, 698.883, 860.37};
const Float_t AliACORDEConstants::fgkInsideModulePositionZ[60] ={
88.372, 348.682, 246.52, 147.039, 48.754, -51.643, -120.342, 0,0, 15.526,
447.027, 348.189, 249.102, 147.577, 47.405, -50.559, -150.334, -251.987, -348.106, -449.947,
448.725, 348.555, 248.541, 148.55, 48.717, -51.631, -151.254, -251.632, -351.217,-451.428,
453.195, 349.899, 249.957, 150.162, 50.603, -49.562, -149.784, -250.068, -349.753, -450.307,
449.871, 351.462, 0, 144.424, 48.172, -52.382, -153.346, -252.389, -353.167, -454.27,
-13.83, 350.436, 248.14, 107.763, 46.085, -85.097, -184.909, -258.298, -349.324, -113.94}; 
const Float_t AliACORDEConstants::fgkCenterModulePositionX[60] = {
-1.733	, 637.6	 , 638.1095, 637.888, 637.8125, 639.579 , 638.63  , 639.332, 639.28  , -0.869,
423.5795, 423.693, 423.795 , 423.452, 423.274,  422.9885, 422.8995, 423.166, 422.7265, 422.1595,
153.119 , 152.362, 152.065 , 151.976, 151.518,  155.316,  151.427,  151.642, 152.465 , 151.93,
-156.171, -152.082,-155.098, -155.141,-154.922, -155.124, -155.629, -154.709,-155.223, -154.841,
-423.037,  -422.772, -426,-422.229,-421.756, -422.053, -422.1545, -422.0375,-422.135,-422.311,
1.637, -643.0205,-643.1815,-642.6285, -642.5675, -610.356, -614.177, -637.256, -636.576, -2.157};
const Float_t AliACORDEConstants::fgkCenterModulePositionY[60] = {
859.72, 592.766, 592.428, 592.81, 592.68, 591.3185, 592.017, 590.053, 590.077, 859.516,
806.5215, 806.3125, 806.312, 806.4895, 806.6705, 807.0455, 807.335, 807.187, 807.615, 808.141,
859.493, 859.044, 859.285, 859.422, 859.396, 859.597, 859.624, 859.677, 859.482, 859.417, 
859.669, 859.494, 859.527, 859.774, 859.486, 859.499, 859.491, 859.505, 859.823, 859.747,
807.771, 807.671, 807.6, 807.5765, 807.9485, 807.2915, 807.82, 807.445, 807.366, 807.331,
859.525, 585.937, 585.616, 586.0805, 586.221, 618.107, 614.02, 591.9165, 592.588, 859.739};
const Float_t AliACORDEConstants::fgkCenterModulePositionZ[60] = {
87.187, 348.1785, 247.3485, 147.058, 48.413, -49.9585, -121.0015, -281.09, -349.005, 16.161,
447.538, 348.676, 250.728, 146.9505, 47.299, -50.535, -150.6745, -249.6215, -348.2345, -449.8365,
449.018, 349.157, 249.406, 149.052, 49.198, -50.944, -150.735, -250.661,-350.989,  -450.826, 
452.428, 349.194, 249.399, 149.286, 49.493, -51.392, -150.955, -251.476, -351.018, -451.487, 
446.97,  350.215, 250.215, 146.4975,44.1585, -50.8225, -147.4875, -254.989, -352.524, -448.606, 
-14.146, 349.752, 249.0625, 105.4995, 44.571, -81.677, -181.26, -258.498, -349.4315, -113.948};
const Float_t AliACORDEConstants::fgkOutsideModulePositionX[60] ={
-149.874, 531.513, 531.852, 531.376, 531.09, 533.167, 531.85, 531.587, 531.878, -148.895,
317.71, 317.626, 317.679, 316.993, 316.376, 316.268, 316.105, 316.735, 316.594, 316.188,
5.086, 4.271, 4.228, 3.875, 3.38, 3.425, 3.402, 3.534, 4.237, 4.199,
-303.888, -303.866, -303.157, -302.97, -302.994, -303.264, -303.36, -302.872, -303.247,-302.837,
-526.69, -526.619, -526.568, -526.325, -526.108, -526.741, -527.336, -526.853, -526.706, -526.709,
-150.248, -748.58, -748.613, -747.93, -748.258, 0, 0,-743.11, -742.565, -150.26};
const Float_t AliACORDEConstants::fgkOutsideModulePositionY[60] ={
860.564, 698.765, 698.093, 698.464, 698.342, 696.999, 697.64, 697.851, 697.969, 860.44,
912.841, 912.515, 912.279, 912.233, 911.86, 912.429, 912.745, 912.864, 913.59, 914.184,
860.42, 860.067, 860.278, 860.318, 860.285, 860.28, 860.466, 860.362, 860.222, 860.322,
860.563, 860.247, 859.991, 860.354, 860.06, 860.056, 859.686, 860.074, 860.463, 860.343,
699.344, 699.337, 698.834, 699.422, 700.129, 699.652, 700.932, 700.096, 699.711, 699.578,
860.159, 479.49, 478.84, 479.09, 479.582, 0,0,485.541, 486.293, 860.416};
const Float_t AliACORDEConstants::fgkOutsideModulePositionZ[60] ={
85.538, 347.675, 248.177, 147.077, 48.072, -48.274, -121.661, -281.827, -347.063, 16.81,
448.049, 349.163, 252.354, 146.324,47.193, -50.511, -151.015, -247.256, -348.363, -449.726,
449.235, 349.438, 249.614, 149.883, 49.413, -49.848, -149.932, -250.458, -350.314, -450.194,
451.868, 348.663, 248.448, 148.685, 48.844, -52.617, -152.185, -251.881, -352.541, -452.442, 
444.069, 348.968, 247.813, 148.571, 40.145, -49.263, -141.629, -257.589, -351.881, -442.942,
-14.142, 349.068, 249.985, 103.236, 43.057, 0, 0, -258.698, -349.539, -114.16}; 
const Float_t AliACORDEConstants::fgkOldModulePositionX[60] = {
  641, 641, 641, 641, 641, 641, 641, 641, 641, 641,
  426, 426, 426, 426, 426, 426, 426, 426, 426, 426,
  153, 153, 153, 153, 153, 153, 153, 153, 153, 153,
  -153, -153, -153, -153, -153, -153, -153, -153, -153,
  -153, -426, -426, -426, -426, -426, -426, -426, -426,
  -426, -426, -644, -644, -644, -644, -644, -619, -623,
  -641, -641, -641};
const Float_t AliACORDEConstants::fgkOldModulePositionY[60] = {
  582, 574, 574, 574, 574, 574, 574, 574, 574, 582,
  789, 789, 789, 789, 789, 789, 789, 789, 789, 789,
  850, 850, 850, 850, 850, 850, 850, 850, 850, 850,
  850, 850, 850, 850, 850, 850, 850, 850, 850, 850,
  789, 789, 789, 789, 789, 789, 789, 789, 789, 789,
  582, 574, 574, 574, 574, 601, 597, 574, 574, 582};
const Float_t AliACORDEConstants::fgkOldModulePositionZ[60] = {
  450, 350, 250, 150, 50, -50, -120, -280, -350, -450,
  450, 350, 250, 150, 50, -50, -150, -250, -350, -450,
  450, 350, 250, 150, 50, -50, -150, -250, -350, -450,
  450, 350, 250, 150, 50, -50, -150, -250, -350, -450,
  450, 350, 250, 150, 50, -50, -150, -250, -350, -450,
  450, 350, 250, 104, 50, -85, -184, -258, -350, -450};


const Float_t AliACORDEConstants::fgkSupportModulePositionX[60] = {
  641, 641, 641, 641, 641, 641, 641, 641, 641, 641,
  426, 426, 426, 426, 426, 426, 426, 426, 426, 426,
  153, 153, 153, 153, 153, 153, 153, 153, 153, 153,
  -153, -153, -153, -153, -153, -153, -153, -153, -153,
  -153, -426, -426, -426, -426, -426, -426, -426, -426,
  -426, -426, -644, -644, -644, -644, -644, -619, -623,
  -641, -641, -641};
const Float_t AliACORDEConstants::fgkSupportModulePositionY[60] = {
  582, 582, 582, 582, 582, 582, 582, 582, 582, 582,
  797, 797, 797, 797, 797, 797, 797, 797, 797, 797,
  850, 850, 850, 850, 850, 850, 850, 850, 850, 850,
  850, 850, 850, 850, 850, 850, 850, 850, 850, 850,
  797, 797, 797, 797, 797, 797, 797, 797, 797, 797,
  582, 582, 582, 582, 582, 609, 605, 582, 582, 582};
const Float_t AliACORDEConstants::fgkSupportModulePositionZ[60] = {
  450, 350, 250, 150, 50, -50, -120, -280, -350, -450,
  450, 350, 250, 150, 50, -50, -150, -250, -350, -450,
  450, 350, 250, 150, 50, -50, -150, -250, -350, -450,
  450, 350, 250, 150, 50, -50, -150, -250, -350, -450,
  450, 350, 250, 150, 50, -50, -150, -250, -350, -450,
  450, 350, 250, 104, 50, -85, -176, -250, -350, -450};
  

const Float_t AliACORDEConstants::fgkOldExtraModulePositionZ[4] = {93.0, 18., -18, -93};
const Float_t AliACORDEConstants::fgkOldExtraModulePositionX = 0.0;
const Float_t AliACORDEConstants::fgkOldExtraModulePositionY = 850.0;
const Int_t AliACORDEConstants::fgkOldModuleElectronicChannel[60] = {
// Old configuration of ACORDE channels in patch panel cards used ONLY in cosmic runs from 2008
/* DCS 0_0 ITS-1*/ 10,
/* DCS 0_1 */ 4,
/* DCS 0_2 */ 8,
/* DCS 0_3 */ 7,
/* DCS 0_4 */ 6,
/* DCS 0_5 */ 5,
/* DCS 0_6 */ 9,
/* DCS 0_7 */ 3,
/* DCS 0_8 */ 2,
/* DCS 0_9 ITS-2*/ 42,
/* DCS 1_0 */ 20,
/* DCS 1_1 */ 19,
/* DCS 1_2 */ 18,
/* DCS 1_3 */ 17,
/* DCS 1_4 */ 16,
/* DCS 1_5 */ 15,
/* DCS 1_6 */ 14,
/* DCS 1_7 */ 13,
/* DCS 1_8 */ 12,
/* DCS 1_9 */ 11,
/* DCS 2_0 */ 60,
/* DCS 2_1 */ 59,
/* DCS 2_2 */ 58,
/* DCS 2_3 */ 57,
/* DCS 2_4 */ 56,
/* DCS 2_5 */ 55,
/* DCS 2_6 */ 54,
/* DCS 2_7 */ 53,
/* DCS 2_8 */ 52,
/* DCS 2_9 */ 51,
/* DCS 3_0 */ 40,
/* DCS 3_1 */ 39,
/* DCS 3_2 */ 38,
/* DCS 3_3 */ 37,
/* DCS 3_4 */ 36,
/* DCS 3_5 */ 35,
/* DCS 3_6 */ 34,
/* DCS 3_7 */ 33,
/* DCS 3_8 */ 32,
/* DCS 3_9 */ 31,
/* DCS 4_0 */ 30,
/* DCS 4_1 */ 29,
/* DCS 4_2 */ 28,
/* DCS 4_3 */ 27,
/* DCS 4_4 */ 26,
/* DCS 4_5 */ 25,
/* DCS 4_6 */ 24,
/* DCS 4_7 */ 23,
/* DCS 4_8 */ 22,
/* DCS 4_9 */ 21,
/* DCS 5_0 ITS-3*/ 1,
/* DCS 5_1 */ 49,
/* DCS 5_2 */ 48,
/* DCS 5_3 */ 47,
/* DCS 5_4 */ 46,
/* DCS 5_5 */ 45,
/* DCS 5_6 */ 44,
/* DCS 5_7 */ 43,
/* DCS 5_8 */ 50,
/* DCS 5_9 ITS-4*/ 41
};



ClassImp(AliACORDEConstants)

//_____________________________________________________________________________
AliACORDEConstants::AliACORDEConstants()
  : TObject()
{
  // Default constructor
}


//_____________________________________________________________________________
AliACORDEConstants* AliACORDEConstants::Instance()
{
// 
// Instance implementacion
//

  if ( !fgInstance ) {
    fgInstance = new AliACORDEConstants;
  }
  return fgInstance;
}

//_____________________________________________________________________________
AliACORDEConstants::~AliACORDEConstants()
{
// 
// destructor for instance
//
  fgInstance = 0;
}

//_____________________________________________________________________________
Float_t AliACORDEConstants::InsideModulePositionX(Int_t i) const
{
//
// Returns the InsideModulePositionX
//
	return fgkInsideModulePositionX[i];
}
//_____________________________________________________________________________
Float_t AliACORDEConstants::InsideModulePositionY(Int_t i) const
{
//
// returns the InsideModulePositionY
//
	return fgkInsideModulePositionY[i];
}
//_____________________________________________________________________________
Float_t AliACORDEConstants::InsideModulePositionZ(Int_t i) const
{
//
// returns the InsideModulePositionZ
//	
	return fgkInsideModulePositionZ[i];
}
//_____________________________________________________________________________
Float_t AliACORDEConstants::CenterModulePositionX(Int_t i) const
{
//
// returns the center module position X
//
	return fgkCenterModulePositionX[i];
}
//_____________________________________________________________________________
Float_t AliACORDEConstants::CenterModulePositionY(Int_t i) const
{
//
// returns the center module position Y
//
	return fgkCenterModulePositionY[i];
}
//_____________________________________________________________________________
Float_t AliACORDEConstants::CenterModulePositionZ(Int_t i) const
{
//
// returns the center module position Z
//	
	return fgkCenterModulePositionZ[i];
}
//_____________________________________________________________________________
Float_t AliACORDEConstants::OutsideModulePositionX(Int_t i) const
{
//
// returns the outside module position x
//
	return fgkOutsideModulePositionX[i];
}
//_____________________________________________________________________________
Float_t AliACORDEConstants::OutsideModulePositionY(Int_t i) const
{
//
// returns the out side module position y
//
	return fgkOutsideModulePositionY[i];
}
//_____________________________________________________________________________
Float_t AliACORDEConstants::OutsideModulePositionZ(Int_t i) const
{
//
// returns the out side module position z
//	
	return fgkOutsideModulePositionZ[i];
}
//_____________________________________________________________________________
Float_t AliACORDEConstants::OldModulePositionX(Int_t i) const
{
  // Module lenght
  return fgkOldModulePositionX[i];
}

//_____________________________________________________________________________
Float_t AliACORDEConstants::OldModulePositionY(Int_t i) const
{
  // Module lenght
  return fgkOldModulePositionY[i];
}
//_____________________________________________________________________________
Float_t AliACORDEConstants::OldModulePositionZ(Int_t i) const
{
  // Module lenght
  return fgkOldModulePositionZ[i];
}


//_____________________________________________________________________________
Float_t AliACORDEConstants::SupportModulePositionX(Int_t i) const
{
  // Module lenght
  return fgkSupportModulePositionX[i];
}

//_____________________________________________________________________________
Float_t AliACORDEConstants::SupportModulePositionY(Int_t i) const
{
  // Module lenght
  return fgkSupportModulePositionY[i];
}
//_____________________________________________________________________________
Float_t AliACORDEConstants::SupportModulePositionZ(Int_t i) const
{
  // Module lenght
  return fgkSupportModulePositionZ[i];
}



Float_t AliACORDEConstants::OldExtraModulePositionX() const
{
  // Module lenght
  return fgkOldExtraModulePositionX;
}

//_____________________________________________________________________________
Float_t AliACORDEConstants::OldExtraModulePositionY() const
{
  // Module lenght
  return fgkOldExtraModulePositionY;
}
//_____________________________________________________________________________
Float_t AliACORDEConstants::OldExtraModulePositionZ(Int_t i) const
{
  // Module lenght
  return fgkOldExtraModulePositionZ[i];
}
//_____________________________________________________________________________
Int_t AliACORDEConstants::OldModuleElectronicChannel(Int_t i) const
{
	// return de ID (electronic channel in ACORDE) of each module
	// acording to the current match between DCS and Electronic nomenclature
	return fgkOldModuleElectronicChannel[i];
}
//_____________________________________________________________________________
Float_t AliACORDEConstants::ModuleLength() const
{
  // Module lenght
  return fgkModuleLength;
}

//_____________________________________________________________________________
Float_t AliACORDEConstants::ModuleWidth() const
{
  // Module width
  return fgkModuleWidth;
}

//_____________________________________________________________________________
Float_t AliACORDEConstants::ModuleHeight() const
{
  // Module height
  return fgkModuleHeight;
}

//_____________________________________________________________________________
Float_t AliACORDEConstants::PlasticLength() const
{
  // Length of the scintillator active zone for a single counter
  return fgkPlasticLength;
}

//_____________________________________________________________________________
Float_t AliACORDEConstants::PlasticWidth() const
{
  // Width of the scintillator active zone for a single counter
  return fgkPlasticWidth;
}

//_____________________________________________________________________________
Float_t AliACORDEConstants::PlasticHeight() const
{
  // Height of the scintillator active zone for a single counter
  return fgkPlasticHeight;
}

Float_t AliACORDEConstants::ProfileWidth() const
{
  // Width of the profile of the Al box
  return fgkProfileWidth;
}

Float_t AliACORDEConstants::ProfileThickness() const
{
  // Thickness of the profile of the Al box
  return fgkProfileThickness;
}


//_____________________________________________________________________________
Float_t AliACORDEConstants::Depth() const
{
  // Alice IP depth
  return fgkDepth;
}
