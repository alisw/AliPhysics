// $Id$
// Category: visualization
//
// According to visualization/management/include/MyVisManager.*
// John Allison 24th January 1998.
// I. Hrivnacova 12.5.98
//
// Renamed to TG4VisManager
// I. Hrivnacova, 3.9.99
//
// Added AliMC implementation
// A. Gheata, 22.2.00
//
// Added OpenGL*Win32, RayTracer (as replacement of RayX) drivers 
// based on G4 suggestions.
// I. Gonzalez, 4.4.2000
//
// Example Visualization Manager description:
//   Implements virtual function RegisterGraphicsSystems.  
//   Exploits C-pre-processor variables G4VIS_USE_DAWN, etc.,
//   which are set by the G4 makefiles if
//   environment variables of the same name are set.
//
// So all you have to do is set environment variables and compile and
//   instantiate this in your main().
//
// Alternatively, you can implement an empty function here and just
//   register the systems you want in your main(), e.g.:
//   G4VisManager* myVisManager = new MyVisManager;
//   myVisManager -> RegisterGraphicsSystem (new MyGraphicsSystem);


#include "TG4Globals.h"

#include <G4LogicalVolumeStore.hh>
#include <G4PhysicalVolumeStore.hh>
#include <G4PhysicalVolumeModel.hh>
#include <G4LogicalVolumeModel.hh>
#include <G4TransportationManager.hh>
#include <G4Material.hh>

#ifdef G4VIS_USE
#include "TG4VisManager.h"

#include <G4VVisManager.hh>

// Supported drivers...

#ifdef G4VIS_USE_DAWN
#include <G4FukuiRenderer.hh>
#endif

#ifdef G4VIS_USE_DAWNFILE
#include <G4DAWNFILE.hh>
#endif

#ifdef G4VIS_USE_OPACS
#include <G4Wo.hh>
#include <G4Xo.hh>
#endif

#ifdef G4VIS_USE_OPENGLX
#include <G4OpenGLImmediateX.hh>
#include <G4OpenGLStoredX.hh>
#endif

#ifdef G4VIS_USE_OPENGLXM
#include <G4OpenGLImmediateXm.hh>
#include <G4OpenGLStoredXm.hh>
#endif

#ifdef G4VIS_USE_OPENGLWIN32
#include "G4OpenGLImmediateWin32.hh"
#include "G4OpenGLStoredWin32.hh"
#endif

#ifdef G4VIS_USE_OIX
#include <G4OpenInventorX.hh>
#endif

#ifdef G4VIS_USE_OIWIN32
#include <G4OpenInventorWin32.hh>
#endif

#ifdef G4VIS_USE_RAYTRACER
#include "G4RayTracer.hh"
#endif

#ifdef G4VIS_USE_VRML
#include <G4VRML1.hh>
#include <G4VRML2.hh>
#endif

#ifdef G4VIS_USE_VRMLFILE
#include <G4VRML1File.hh>
#include <G4VRML2File.hh>
#endif

TG4VisManager::TG4VisManager(G4int verboseLevel) {
//  
  fVerbose = verboseLevel; 
  fColourFlag = true;
}

TG4VisManager::TG4VisManager(const TG4VisManager& right) {
// 
  TG4Globals::Exception(
    "Attempt to copy TG4VisManager singleton.");
}

TG4VisManager::~TG4VisManager() {
//
}  

// operators

TG4VisManager& TG4VisManager::operator=(const TG4VisManager& right) 
{
  // check assignement to self
  if (this == &right) return *this;

  TG4Globals::Exception(
    "Attempt to assign TG4VisManager singleton.");
    
  return *this;  
}    
          
// private methods

void TG4VisManager::RegisterGraphicsSystems() 
{
// Registers the graphics systems.
// ---

  // all created graphics system instances are
  // deleted in G4VisManager::~G4VisManager()
#ifdef G4VIS_USE_DAWN
  RegisterGraphicsSystem(new G4FukuiRenderer);
#endif

#ifdef G4VIS_USE_DAWNFILE
  RegisterGraphicsSystem(new G4DAWNFILE);
#endif

#ifdef G4VIS_USE_OPACS
  RegisterGraphicsSystem(new G4Wo);
  RegisterGraphicsSystem(new G4Xo);
#endif

#ifdef G4VIS_USE_OPENGLX
  RegisterGraphicsSystem(new G4OpenGLImmediateX);
  RegisterGraphicsSystem(new G4OpenGLStoredX);
#endif

#ifdef G4VIS_USE_OPENGLXM
  RegisterGraphicsSystem(new G4OpenGLImmediateXm);
  RegisterGraphicsSystem(new G4OpenGLStoredXm);
#endif

#ifdef G4VIS_USE_OPENGLWIN32
  RegisterGraphicsSystem (new G4OpenGLImmediateWin32);
  RegisterGraphicsSystem (new G4OpenGLStoredWin32);
#endif

#ifdef G4VIS_USE_OIX
  RegisterGraphicsSystem(new G4OpenInventorX);
#endif

#ifdef G4VIS_USE_OIWIN32
  RegisterGraphicsSystem(new G4OpenInventorWin32);
#endif

#ifdef G4VIS_USE_RAYTRACER
  RegisterGraphicsSystem (new G4RayTracer);
#endif

#ifdef G4VIS_USE_VRML
  RegisterGraphicsSystem(new G4VRML1);
  RegisterGraphicsSystem(new G4VRML2);
#endif

#ifdef G4VIS_USE_VRMLFILE
  RegisterGraphicsSystem(new G4VRML1File);
  RegisterGraphicsSystem(new G4VRML2File);
#endif

  if (fVerbose > 0) {
    cout <<
      "\nYou have successfully chosen to use the following graphics systems."
	 << endl;
    PrintAvailableGraphicsSystems();
  }
}

//---------------------------------------------------------------
// private methods
//---------------------------------------------------------------


G4RWTPtrOrderedVector<G4LogicalVolume> TG4VisManager::GetLVList(G4String name)
{
// Get function returning the list of logical volumes
// associated to NAME; G4 built clones of a G3 volume (identified 
// with NAME_NUMBER will be added to the list)  
//  NAME can be the name of a logical or physical volume
// ---

 G4RWTPtrOrderedVector <G4LogicalVolume> lvList;
 G4LogicalVolumeStore* pLVStore = G4LogicalVolumeStore::GetInstance();
 G4LogicalVolume* pLV = 0; 
 if (pLVStore)
 {
   for (G4int i=0; i<pLVStore->entries(); i++)
   {
     pLV = pLVStore->at(i);  
     if (CaseInsensitiveEqual(name,pLV->GetName())) 
     {
       if (!lvList.contains(pLV)) lvList.append(pLV);
     }
   }
 }  
 if (!lvList.isEmpty()) return lvList;
 G4PhysicalVolumeStore* pPVStore = G4PhysicalVolumeStore::GetInstance();
 G4VPhysicalVolume* pPV = 0;
 if (pPVStore) 
 {
   for (G4int i=0; i<pPVStore->entries(); i++)
   {
     pPV = pPVStore->at(i); 
     if (CaseInsensitiveEqual(name,pPV->GetName())) 
     {
       pLV = pPV->GetLogicalVolume();
       if (!lvList.contains(pLV)) lvList.append(pLV);
     }     
   }
 }  
 return lvList;
}


G4RWTPtrOrderedVector<G4VPhysicalVolume> TG4VisManager::GetPVList(G4String name)
{
// Get function returning the physical volume pointer for NAME
// ---

  G4RWTPtrOrderedVector <G4VPhysicalVolume> pvList;
  G4PhysicalVolumeStore* pPVStore = G4PhysicalVolumeStore::GetInstance();
  if (!pPVStore)
  {
    TG4Globals::Warning("TG4VisManager::Gdraw : No volume store !");
    return pvList;
  }
  G4VPhysicalVolume* pPV = 0;
  for (G4int i=0; i<pPVStore->entries(); i++)
  {
    pPV = pPVStore->at(i);
    if (CaseInsensitiveEqual(name,pPV->GetName()))
    {
      if (!pvList.contains(pPV)) pvList.append(pPV);
    }  
  }
  return pvList;
}


G4bool TG4VisManager::CaseInsensitiveEqual(const G4String string1,
					   const G4String string2)
{
// Case insensitive comparison of 2 strings
// ---

  G4int i;
  if (string1 == string2) return true;
//  if (string1.length() != string2.length()) return false;
  for (i=0; i<string1.length(); i++)
  {
    G4int diff = abs((G4int)string1(i)-(G4int)string2(i));
    if (diff && (diff!=32)) return false;
  }
  if (string2.length() > string1.length())
  {
    if (string2(i) == '_') return true;
    return false;
  }    
  return true;
}
 

void TG4VisManager::SetAtt4Daughters(G4LogicalVolume* const lv, 
				     const TG3Attribute att, const G4int val)
{
// Iterator for setting a visual attribute for all daughters
// ---

  SetG4Attribute(lv,att,val);  
  
  G4String lvName = lv->GetName();
  G4int nOfDaughters = lv->GetNoDaughters();
  if (nOfDaughters>0)
  {
    G4String previousName = "";
    for (G4int i=0; i<nOfDaughters; i++) 
    { 
      G4LogicalVolume* lvd = lv->GetDaughter(i)->GetLogicalVolume(); 
      G4String currentName = lvd->GetName();
      if (currentName != lvName && currentName != previousName)
      {
        SetAtt4Daughters(lvd, att, val);
        previousName = currentName;
      } 
    }
  }   
}


G4bool TG4VisManager::IsSharedVisAttributes(const G4LogicalVolume* pLV)
{
// Function seeking if the volume's visible attributes are shared with
//  other volumes
// ---

  G4LogicalVolumeStore* pLVStore = G4LogicalVolumeStore::GetInstance();
  G4LogicalVolume* pLVCurrent = 0;
  const G4VisAttributes* pVisAtt = pLV->GetVisAttributes();
  if (!pVisAtt) return false;
  for (G4int i=0; i<pLVStore->entries(); i++)
  {
    pLVCurrent = pLVStore->at(i);
    if (pLVCurrent != pLV)
    {
      if (pLVCurrent->GetVisAttributes() == pVisAtt) 
      {
        return true;
      }
    }
  }
  return false;
}


void TG4VisManager::SetG4Attribute(G4LogicalVolume* const lv,
				   const TG3Attribute att, const G4int val)
{
// Set the G4 attribute fo volume LV accordingly to the G3 description
//  of (att- val)    
// --

 if (!lv) return;
 // Dupplicating old vis. attributes    
 const G4VisAttributes* visAttributes = lv->GetVisAttributes();
 G4VisAttributes* newVisAttributes;
 if (!visAttributes)
   newVisAttributes    = new G4VisAttributes(false);
 else
   newVisAttributes    = new G4VisAttributes(visAttributes);

 const G4int absVal = abs(val);	// the functionality is given by the abs value

 // Default visible attributes
 G4double red(0),green(0),blue(0); // default is black
 G4bool isVisible(false);
 G4bool isDaughtersInvisible(false);
 G4VisAttributes::LineStyle lineStyle = G4VisAttributes::unbroken;
 G4double lineWidth = 1.0;
 G4bool isForceDrawingStyle(false);
 G4VisAttributes::ForcedDrawingStyle drawingStyle = G4VisAttributes::wireframe;
 
 // a 'hardcopy' of old vis attributes is needed because the copy constructor
 // resets to defaults some of the data members of G4VisAttributes class
 if (visAttributes)
 {
   isVisible	      = visAttributes->IsVisible();
   isDaughtersInvisible = visAttributes->IsDaughtersInvisible();
   red   = visAttributes->GetColour().GetRed();
   green = visAttributes->GetColour().GetGreen(); 
   blue	 = visAttributes->GetColour().GetBlue(); 	// old RGB components
   lineStyle = visAttributes->GetLineStyle();
   lineWidth = visAttributes->GetLineWidth();
   isForceDrawingStyle = visAttributes->IsForceDrawingStyle();
   if (isForceDrawingStyle) 
       drawingStyle = visAttributes->GetForcedDrawingStyle();
   
 }
 G4double luminosityBin(0.04),  // bin for luminosity 
          luminosity(0); 	// colour luminosity

 // Delete old vis. attributes if they are not shared
 if (visAttributes && !IsSharedVisAttributes(lv)) delete visAttributes;

 // Set the required attribute
 switch (att)
 {
  case kSEEN:
   switch (val)
   {
    case  0:
     isVisible = false; 
     break;
    case  1:
     isVisible = true;  
     break;
    case -1:
     isVisible = false; 
     break;
    case -2:
     isVisible = false; 
     break;
    default:
     isVisible = false; 
   }       
   break;
  case kLSTY:
   switch (absVal)
   {
    case 1:
     lineStyle = G4VisAttributes::unbroken; break;
    case 2:
     lineStyle = G4VisAttributes::dashed; break;
    case 3:
     lineStyle = G4VisAttributes::dotted; break;
    default:
     if (fVerbose > 0)
      G4cout << "TG4VisManager::Gsatt() Usage of LSTY :" << endl
    	     << "ATT = 1,2,3 means line unbroken, dashed or dotted" << endl
	     << "any other value resets to the default : unbroken" << endl;
     lineStyle = G4VisAttributes::unbroken;
   }       
   break;
  case kLWID:
   lineWidth = absVal;
   if (lineWidth > 7) lineWidth = 7;
   if (fVerbose > 0) 
       G4cout << "TG4VisManager::Gsatt() Usage for LWID :" << endl
              << "  The VAL you supply means the width of lines in pixels "
	      << "for the screen and in 0.1*mm for paper." << endl
	      << "  Negative values means the same, but for all daughters"
	      << endl;
   break;
  case kCOLO:
   if (absVal < 8)	// G3 base colours
   {
    switch (absVal)
    {
     case 1:
       red=0; green=0; blue=0;    	//black
       break;	   	        
     case 2:
       red=1; green=0; blue=0;    	//red
       break;	   	        
     case 3:
       red=0; green=1; blue=0;    	//green
       break;	   	        
     case 4:
       red=0; green=0; blue=1;    	//blue
       break;	   	        
     case 5:
       red=1; green=1; blue=0;    	//yellow
       break;	   	        
     case 6:
       red=1; green=0; blue=1;    	//violet
       break;	   	        
     case 7:
       red=0; green=1; blue=1;    	//lightblue (almost !)
    }
    luminosity = 0.;
   }
   if (absVal>=8 && absVal<=16)
   {
     red=0; green=0; blue=0;
     luminosity = (absVal-7)*luminosityBin;
   }
   if (absVal>=17 && absVal<=41)
   {
     red=1; green=0; blue=0;
     luminosity = (absVal-16)*luminosityBin;
   }
   if (absVal>=67 && absVal<=91)
   {
     red=0; green=1; blue=0;
     luminosity = (absVal-66)*luminosityBin;
   }
   if (absVal>=117 && absVal<=141)
   {
     red=0; green=0; blue=1;
     luminosity = (absVal-116)*luminosityBin;
   }
   if (absVal>=42 && absVal<=66)
   {
     red=1; green=1; blue=0;
     luminosity = (absVal-41)*luminosityBin;
   }
   if (absVal>=142 && absVal<=166)
   {
     red=1; green=0; blue=1;
     luminosity = (absVal-141)*luminosityBin;
   }
   if (absVal>=92 && absVal<=116)
   {
     red=0; green=1; blue=1;
     luminosity = (absVal-91)*luminosityBin;
   }
   if (red < luminosityBin) 	red += luminosity;
   if (green < luminosityBin) green += luminosity;
   if (blue < luminosityBin)   blue += luminosity;
   break;
  case kFILL:
   isForceDrawingStyle = true;
   switch (absVal)
   {
    case 0:
     drawingStyle = G4VisAttributes::wireframe;
     break;
    case 1:
     drawingStyle = G4VisAttributes::solid;
     break;
    default:
     if (fVerbose > 0)
         G4cout << "TG4VisManager::Gsatt() FILL usage :" << endl
	        << "  The FILL values you can supply are only :" << endl
		<< "+/- 1 : forces wireframe drawing (default)" << endl
		<< "+/- 2 : forces solid drawing" << endl
		<< "other values sets the drawing style to solid"
		<< endl;		
     drawingStyle = G4VisAttributes::solid;
   }      
 }
 // Register vis. attributes
 newVisAttributes->SetVisibility(isVisible);
 newVisAttributes->SetDaughtersInvisible(isDaughtersInvisible);
 newVisAttributes->SetColour(red,green,blue);
 newVisAttributes->SetLineStyle(lineStyle);
 newVisAttributes->SetLineWidth(lineWidth);
 if (drawingStyle == G4VisAttributes::wireframe) 
       newVisAttributes->SetForceWireframe(isForceDrawingStyle);
 if (drawingStyle == G4VisAttributes::solid) 
       newVisAttributes->SetForceSolid(isForceDrawingStyle);
 
 lv->SetVisAttributes(newVisAttributes);
} 
  
//-----------------------------------------------------------------
// functions for drawing
//-----------------------------------------------------------------

void TG4VisManager::DrawOneSpec(const char* name)
{
// Function called when one double-clicks on a volume name
// in a TPaveLabel drawn by Gdtree
// ---

 G4cout << "TG4VisManager::DrawOneSpec() Not yet implemented";
}


void TG4VisManager::SetColors()
{
// Function for setting default volume colours
// ---

  G4LogicalVolumeStore* pLVStore = G4LogicalVolumeStore::GetInstance();
  const G4LogicalVolume* pLV = 0; 
  if (!pLVStore)
  {
    TG4Globals::Warning("TG4VisManager::SetColors : Ignored, geometry not built.");
    return;
  }
  // parse the LV tree and set colours according to material density
  for (G4int i=0; i<pLVStore->entries(); i++)
  {
    pLV = pLVStore->at(i);
//    G4cout << "VOLUME : " << pLV->GetName() << endl;
    const G4Material* pMaterial = pLV->GetMaterial();
    const G4State state = pMaterial->GetState();
    G4double density = (pMaterial->GetDensity())*cm3/g;
    G4String nState = "Undefined";
    G4int colour = 1;			//black by default
    G4double luminosity = 0.;
    if (state == kStateUndefined) 
    {
      nState = "Undefined";
    }
    if (state == kStateSolid) 
    {
      nState = "Solid";
      if (density < 2)
      {
        colour = 17;	//red
	luminosity = 25 - 25*density/2;
      }
      else if (density < 3)
      {
        colour = 117;	//blue
	luminosity = 25 - 25*(density-2);
      }
      else if (density < 10)
      {
        colour = 67;	//green
	luminosity = 25 - 25*(density-5)/5;
      }
      else if (density < 15)
      {
        colour = 92;	//cyan
	luminosity = 25 - 25*(density-10)/5;
      }
      else if (density < 20)
      {
        colour = 8;	//black
	luminosity = 9 - 9*(density-15)/5;
      }
    }
    if (state == kStateLiquid) 
    {
      nState = "Liquid";
      colour = 142;	//violet
      luminosity = 25 - 25*density/2;
    }
    if (state == kStateGas) 
    {
      nState = "Gas";
      if (density < 0.001)  {colour = 42;}	//yellow
      else if (density < 0.002) {colour = 27;}	//light red
      else if (density < 0.003) {colour = 77;}	//light green
      else {colour = 102;}			//light cyan
      luminosity = 0;
    }
    if (luminosity < 0) luminosity=0;
    colour += (G4int)luminosity;
    //  Setting the corresponding colour
    Gsatt(pLV->GetName(),"COLO",colour);
  }
} 
  

void TG4VisManager::Gsatt(const char* name, const char* att, Int_t val)
{
// Geant3 description :
//  
//    NAME   Volume name
//    IOPT   Name of the attribute to be set
//    IVAL   Value to which the attribute is to be set
//  
//    name= "*" stands for all the volumes.
//    iopt can be chosen among the following :
//    kWORK, kSEEN, kLSTY, kLWID, kCOLO, kFILL, kSET, kDET, kDTYP
// ---

 G4int ival = val;
 G4LogicalVolume* lv = 0; 
 G4RWTPtrOrderedVector<G4LogicalVolume> lvList;
 G4String sname(name),
          satt(att);		

 // seek for known attributes
 TG3Attribute attribute = kUNKNOWN;
 if (CaseInsensitiveEqual(att,"WORK"))
 {
   G4String message = "TG4VisManager::Gsatt: G3Attribute ";
   message += satt + " not used in G4";
   TG4Globals::Warning(message);
   return;
 }
 if (CaseInsensitiveEqual(att,"SEEN"))	attribute = kSEEN;  
 if (CaseInsensitiveEqual(att,"LSTY"))	attribute = kLSTY;  
 if (CaseInsensitiveEqual(att,"LWID"))	attribute = kLWID;
 if (CaseInsensitiveEqual(att,"COLO"))	attribute = kCOLO;  
 if (CaseInsensitiveEqual(att,"FILL"))	attribute = kFILL;  
 if (CaseInsensitiveEqual(att,"SET"))
 {
   G4String message = "TG4VisManager::Gsatt: G3Attribute ";
   message += satt + " not used in G4";
   TG4Globals::Warning(message);
   return;
 }
 if (CaseInsensitiveEqual(att,"DET"))
 {
   G4String message = "TG4VisManager::Gsatt: G3Attribute ";
   message += satt + " not used in G4";
   TG4Globals::Warning(message);
   return;
 }
 if (CaseInsensitiveEqual(att,"DTYP"))
 {
   G4String message = "TG4VisManager::Gsatt: G3Attribute ";
   message += satt + " not used in G4";
   TG4Globals::Warning(message);
   return;
 }
 if (attribute == kUNKNOWN)
 {
   G4String message = "TG4VisManager::Gsatt: G3Attribute ";
   message += satt + " unknown";
   TG4Globals::Warning(message);
   return;
 }
 G4bool	 doForDaughters(false),		// tree iterator flag 
         doForAll(false),		// activated if NAME is "*" 
	 topVisible(false);		// activated for kSEEN/-2
 if (sname == "*")	doForAll = true;
 if (val < 0 && sname!="*") doForDaughters = true;
 if (attribute==kSEEN && val==-2) topVisible = true;
 
 // parse all the tree
 if (doForAll)
 {
     G4LogicalVolumeStore* pLVStore = G4LogicalVolumeStore::GetInstance();
     for (G4int i=0; i<pLVStore->entries(); i++)
     {
         lv = pLVStore->at(i);
	 SetG4Attribute(lv,attribute,ival);
     }
     return;
 }

 // get the logical volume pointer corresponding to NAME
 lvList = GetLVList(name);
 if (lvList.isEmpty())
 {
     G4String message = "TG4VisManager::Gsatt(): Ignored\n";
     message += "    Logical volume " + sname + " has not been found.\n";
     TG4Globals::Warning(message); 
     return;
 }
// set attribute for all descendents
 if (doForDaughters)
 {
   for (G4int i=0; i<lvList.entries(); i++)
   {
     lv = lvList[i];
     SetAtt4Daughters(lv,attribute,ival);
   }     
 }
 else
 {
   for (G4int i=0; i<lvList.entries(); i++)
   {
     lv = lvList[i];
     SetG4Attribute(lv,attribute,ival);
   }     
 }
 if (topVisible) 
 {
   for (G4int i=0; i<lvList.entries(); i++)
   {
     lv = lvList[i];
     SetG4Attribute(lv,attribute,1); 
   }
 }       
} 


void TG4VisManager::Gdraw(const char *name,Float_t theta, Float_t phi, Float_t psi,
		    Float_t u0,Float_t v0,Float_t ul,Float_t vl)
{ 
// Draw the physical volume NAME and all descendents;
// Mandatory : the graphics system, scene and view must be
//	initialized, e.g. "/vis~/create_view/new_graphics_system OGLSX";
//	Any call of Gdraw() will use the current graphics system and
//	the current window.
// The result will be a centered view drawing of the designated volume,
//	lights moving with camera, viewpoint direction given by theta/phi
//	and rotation on the screen given by psi;
// The u0, v0, ul, vl factors are ignored since the object will be 
// 	automatically centered and will be confortable in the window
//	at any viewing angle.
//
// check if G4 graphics is ready for drawing
// ---

  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
  if (!pVVisManager) { 
    TG4Globals::Warning(
       "TG4VisManager::Gdraw: Ignored - No graphics driver is built. ");
    return;
  }     
  if (NeedSetColours())
  {
    SetColors();
    SetColourFlag(false);
  }
  
  const G4double rad = M_PI/180.;
  G4RWTPtrOrderedVector<G4VPhysicalVolume> pvList;
  G4String sname(name);
  G4bool  successful 		= false;

  pvList = GetPVList(name);
  if (pvList.isEmpty())
  {
    G4String message = "TG4VisManager::Gdraw() :\n";
    message += "Volume " + sname + " not found. Bailing out";
    TG4Globals::Warning(message); 
    return;
  }

  G4VPhysicalVolume *pPV = 0;

  // clear the current scene if not empty
  if (!fpScene->IsEmpty()) fpScene->Clear();

  // create and add object's model list to the runtime-duration model 
  // list and draw it
  // (it is deleted in the VisManager destructor within 
  // all the RWTPtrOrderedVectors of the scene)
  for (G4int i=0; i<pvList.entries(); i++)
  { 
    pPV = pvList[i];
    G4LogicalVolume* pLV = pPV->GetLogicalVolume();
    G4VSolid* pSolid = pLV->GetSolid();

    successful = fpScene->AddRunDurationModel(new G4PhysicalVolumeModel(pPV));
    if (!successful) 
    {
      G4String message = "TG4VisManager::Gdraw() Could not add ";
      message += pPV->GetName() + " to the drawing list. Probably ";
      message += "it is already in the list.";
      TG4Globals::Warning(message);
    }  
  }
  // get the standard target point of the scene
  const G4Point3D targetPoint = fpScene->GetStandardTargetPoint();

  // set the viewpoint and the rotation on the screen
  G4Vector3D viewpointDirection(sin(theta*rad)*cos(phi*rad), 
                                sin(theta*rad)*sin(phi*rad), cos(theta*rad)); 
  G4Vector3D upVector(sin(psi*rad), cos(psi*rad),0);			       

  // set and register view parameters to the viewer
  
  fVP.SetLightsMoveWithCamera(true);
  fVP.SetViewGeom();
  fVP.UnsetViewHits();
  fVP.UnsetViewDigis();
  fVP.SetNoOfSides(48);
  fVP.SetCurrentTargetPoint(targetPoint);
  fVP.SetViewpointDirection(viewpointDirection);
  fVP.SetUpVector(upVector);
  fVP.SetDensityCulling(true);
  fpViewer->SetViewParameters(fVP);

  if (IsValidView())
  {
    fpSceneHandler->SetScene(fpScene);
    fpSceneHandler->SetCurrentViewer(fpViewer);
    fpViewer->DrawView();
    fpViewer->ShowView();
  }
  else TG4Globals::Warning(
         "TG4VisManager::Gdraw: Ignored - Failed to register volume"); 
}
#endif //G4VIS_USE
