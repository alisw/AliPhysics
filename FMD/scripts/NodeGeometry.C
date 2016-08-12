//____________________________________________________________________
//
//
// $Id$
//
// Script I used for rapid prototyping of the FMD3 geometry - in
// particular the support cone 
//
/** @defgroup FMD_node_geom Simple geometry
    @ingroup FMD_script
*/
#include <TGeometry.h>
#include <TNode.h>
#include <TXTRU.h>
#include <TTUBE.h>
#include <TTUBS.h>
#include <TPCON.h>
#include <TBRIK.h>
#include <TCanvas.h>
#include <vector>
#include <algorithm>
#include <cmath>
#include <iostream>

//____________________________________________________________________
/** @brief A 2D point
    @ingroup FMD_node_geom
 */
struct point_t
{
  point_t(double x=0, double y=0) : first(x), second(y) {}
  double first;
  double second;
};

//____________________________________________________________________
/** @brief Shape of a ring
    @ingroup FMD_node_geom
 */
struct Ring 
{
  // typedef std::pair<double,double> point_t;
  typedef std::vector<point_t>     points_t;
  /** Constructor 
      @param rL         Lower radius
      @param rH         Higer radius
      @param theta      Opening angle
      @param waferR     Wafer radius
      @param siThick    Silicon thickness 
      @param staggering Staggering of modules */
  Ring(double rL, double rH, double theta,  double waferR, 
       double siThick, double staggering)
    : fStaggering(staggering),
      fInnerRadius(rL), 
      fOuterRadius(rH), 
      fAngle(theta), 
      fRadius(waferR),
      fThickness(siThick), 
      fVerticies(6)
  {
    double tan_theta  = tan(fAngle * TMath::Pi() / 180.);
    double tan_theta2 = pow(tan_theta,2);
    double r2         = pow(fRadius,2);
    double ir2        = pow(fInnerRadius,2);
    double or2        = pow(fOuterRadius,2);
    double y_A        = tan_theta * fInnerRadius;
    double x_D        = fInnerRadius + sqrt(r2 - tan_theta2 * ir2);
    double x_D2       = pow(x_D,2);
    double y_B        = sqrt(r2 - or2 + 2 * fOuterRadius * x_D - x_D2);
    double x_C        = ((x_D + sqrt(-tan_theta2 * x_D2 
				     + r2 * (1 + tan_theta2))) 
			 / (1 + tan_theta2));
    double y_C        = tan_theta * x_C;
    
    fVerticies[0] = point_t(fInnerRadius,  y_A);
    fVerticies[1] = point_t(x_C,           y_C);
    fVerticies[2] = point_t(fOuterRadius,  y_B);
    fVerticies[3] = point_t(fOuterRadius, -y_B);
    fVerticies[4] = point_t(x_C,          -y_C);
    fVerticies[5] = point_t(fInnerRadius, -y_A);
  }
  /** Destructor */
  virtual ~Ring() 
  {
    fVerticies.clear();
  }
  /** Create a shape 
      @return pointer to new shape  */
  TShape* CreateShape()
  {
    std::cout << "Creating Module shape" << std::flush;
    TXTRU* moduleShape = new TXTRU("Module","Module", "", 6, 2);
    for (Int_t i = 0; i  < 6; i++) {
      std::cout << "." << std::flush;
      point_t& p = fVerticies[i];
      moduleShape->DefineVertex(i, p.first, p.second);
    }
    moduleShape->DefineSection(0, -fThickness/2, 1, 0, 0);
    moduleShape->DefineSection(1,  fThickness/2, 1, 0, 0);
    std::cout << std::endl;
    return (TShape*)moduleShape;
  }
  /** Create a node that represents a ring. 
      @return Node */
  TNode* CreateRing(const char* name, double z) 
  {
    std::cout << "Creating Ring node for " << name << std::flush;
    double bredth = fStaggering + fThickness;
    TShape* ringShape   = new TTUBE(Form("%sShape", name), "Ring Shape", 
				    "", fInnerRadius, 
				    fOuterRadius,bredth/2);
    TNode*  ringNode    = new TNode(Form("%sNode", name), "Ring Node", 
				    ringShape, 0, 0, z+bredth/2, 0);
    TShape* moduleShape = CreateShape();
    Int_t n = Int_t(360 / 2 / fAngle);
    for (Int_t i = 0; i < n; i++) {
      std::cout << "." << std::flush;
      ringNode->cd();
      Double_t theta  = 2  * fAngle * i;
      Double_t z      = -(bredth+fThickness)/2+(i%2?0:fStaggering);
      TRotMatrix* rot = new TRotMatrix(Form("%sRotation%02d", name, i), 
				       "Rotation", 90, theta, 90, 
				       fmod(90 + theta, 360), 0, 0);
      TNode* moduleNode = new TNode(Form("%sModule%02d", name, i), 
				    "Module", moduleShape, 0, 0, z,
				    rot);
      moduleNode->SetFillColor(2);
      moduleNode->SetLineColor(2);
      moduleNode->SetLineWidth(2);
    }
    std::cout << std::endl;
    ringNode->SetVisibility(0);
    return ringNode;
  }
  double fStaggering;
  /** Inner radius */
  double fInnerRadius;
  /** Outer radius */
  double fOuterRadius;
  /** Opening angle (in degrees) */
  double fAngle;
  /** Radius (in centimeters) */
  double fRadius;
  /** Thickness */
  double fThickness;
  /** List of verticies */
  points_t fVerticies;
};

//____________________________________________________________________
/** @brief Shape of a detector
    @ingroup FMD_node_geom
 */
struct Detector
{
  /** Constructor 
      @param id 
      @param inner 
      @param outer */
  Detector(Ring* inner, double iZ, Ring* outer=0, double oZ=0) 
    : fInner(inner), fInnerZ(iZ), fOuter(outer), fOuterZ(oZ)
  {}
  /** Destructor */
  virtual ~Detector() {}
  /** Create rings */
  virtual void CreateRings() 
  {
    if (fInner) fInner->CreateRing("inner", fInnerZ);
    if (fOuter) fOuter->CreateRing("outer", fOuterZ);
  }
  /** Create a node that represents the support */
  virtual void CreateSupport(double) { }
  /** Pointer to inner ring */
  Ring* fInner;
  /** Position in z of inner ring */ 
  double fInnerZ;
  /** Pointer to outer ring */
  Ring* fOuter;
  /** Position in z of inner ring */ 
  double fOuterZ;
};

//____________________________________________________________________
/** @brief FMD3 simple node geometry 
    @ingroup FMD_node_geom
 */
struct FMD3 : public Detector
{
  /** Constructor 
      @param inner Inner ring representation  
      @param outer Outer ring representation */
  FMD3(Ring* inner, Ring* outer)
    : Detector(inner, -62.8,outer, -75.2)
  {
    fNoseRl     = 5.5;
    fNoseRh     = 6.7;
    fNoseDz     = 2.8 / 2;
    fNoseZ      = -46;
    fConeL      = 30.9;
    fBackRl     = 61 / 2;
    fBackRh     = 66.8 /2;
    fBackDz     = 1.4 / 2;
    fBeamDz     = .5 / 2;
    fBeamW      = 6;
    fFlangeR    = 49.25;
  }
  virtual ~FMD3() {}
  void CreateRings()
  {
    double zdist      = fConeL - 2 * fBackDz - 2 * fNoseDz;
    double tdist      = fBackRh - fNoseRh;
    double alpha      = tdist / zdist;
    double x, rl, rh, z;
    z  = fNoseZ - fConeL / 2;
    TPCON* fmd3Shape = new TPCON("fmd3Shape", "FMD 3 Shape", "", 0, 360, 7);
    x  = fNoseZ;
    rl = fNoseRl;
    rh = fNoseRh;
    fmd3Shape->DefineSection(0, x - z, rl, rh);
    x  = fNoseZ-2*fNoseDz;
    fmd3Shape->DefineSection(1, x - z, rl, rh);
    x  = fInnerZ - fInner->fStaggering - fInner->fThickness; 
    rl = fInner->fInnerRadius;
    rh = fNoseRh + alpha * TMath::Abs(x-fNoseZ + 2 * fNoseDz);
    fmd3Shape->DefineSection(2, x - z, rl, rh);
    x  = fOuterZ;
    rl = fOuter->fInnerRadius;
    rh = fBackRh;
    fmd3Shape->DefineSection(3, x - z, rl, rh);
    x  = fNoseZ - zdist - 2 * fNoseDz;
    rl = fOuter->fInnerRadius;
    rh = fBackRh;
    fmd3Shape->DefineSection(4, x - z, rl, rh);
    x  = fNoseZ - zdist - 2 * fNoseDz;
    rl = fOuter->fInnerRadius;
    rh = fFlangeR;
    fmd3Shape->DefineSection(5, x - z, rl, rh);
    x  = fNoseZ - fConeL;
    rl = fOuter->fInnerRadius;
    rh = fFlangeR;
    fmd3Shape->DefineSection(6, x - z, rl, rh);

    TNode* fmd3Node = new TNode("fmd3Node", "FMD3 Node", fmd3Shape,0,0,z,0);
    fmd3Node->SetLineColor(11);
    fmd3Node->SetFillColor(11);
    fmd3Node->SetVisibility(1);
    fmd3Node->cd();
    if (fInner) fInner->CreateRing("inner", fInnerZ-z);
    fmd3Node->cd();
    if (fOuter) fOuter->CreateRing("outer", fOuterZ-z);
    fmd3Node->cd();
    CreateSupport(fNoseZ - z);
  }
  
  /** Create support volumes  */
  void CreateSupport(double noseZ) 
  {
    TShape* noseShape = new TTUBE("noseShape", "Nose Shape", "", 
				  fNoseRl, fNoseRh, fNoseDz);
    TNode*  noseNode  = new TNode("noseNode", "noseNode", noseShape, 
				  0, 0, noseZ - fNoseDz, 0);
    noseNode->SetLineColor(0);
    double  zdist = fConeL - 2 * fBackDz - 2 * fNoseDz;
    double  tdist = fBackRh - fNoseRh;
    double  beamL = TMath::Sqrt(zdist * zdist + tdist * tdist);
    double  theta = -TMath::ATan2(tdist, zdist);
    TShape* backShape = new TTUBE("backShape", "Back Shape", "", 
				  fBackRl, fBackRh, fBackDz);
    TNode*  backNode  = new TNode("backNode", "backNode", backShape, 
				  0, 0, noseZ-2*fNoseDz-zdist-fBackDz, 0);
    backNode->SetLineColor(0);
    TShape* beamShape = new TBRIK("beamShape", "beamShape", "", 
				  fBeamDz, fBeamW / 2 , beamL / 2);
    Int_t    n = 8;
    Double_t r = fNoseRl + tdist / 2;
    for (Int_t i = 0; i < n; i++) {
      Double_t phi   = 360. / n * i;
      Double_t t     = 180. * theta / TMath::Pi();
      TRotMatrix* beamRotation = new TRotMatrix(Form("beamRotation%d", i), 
						Form("beamRotation%d", i),
						180-t,phi,90,90+phi,t,phi);
      TNode* beamNode = new TNode(Form("beamNode%d", i), 
				  Form("beamNode%d", i), beamShape, 
				  r * TMath::Cos(phi / 180 * TMath::Pi()), 
				  r * TMath::Sin(phi / 180 * TMath::Pi()), 
				  noseZ-2*fNoseDz-zdist/2, beamRotation);
      beamNode->SetLineColor(0);
    }
    Double_t flangel    = (fFlangeR - fBackRh) / 2;
    TShape* flangeShape = new TBRIK("flangeShape", "FlangeShape", "", 
				    flangel, fBeamW / 2, fBackDz);
    n = 4;
    r = fBackRh + flangel;
    for (Int_t i = 0; i < n; i++) {
      Double_t phi = 360. / n * i + 180. / n;
      TRotMatrix* flangeRotation = new TRotMatrix(Form("flangeRotation%d", i),
						  Form("Flange Rotation %d",i),
						  90,phi,90,90+phi,0,0);
      TNode* flangeNode = new TNode(Form("flangeNode%d", i), 
				    Form("flangeNode%d", i), 
				    flangeShape,
				    r * TMath::Cos(phi / 180 * TMath::Pi()), 
				    r * TMath::Sin(phi / 180 * TMath::Pi()),
				    noseZ-2*fNoseDz-zdist-fBackDz, 
				    flangeRotation);
      flangeNode->SetLineColor(0);
				  
    }
  }
  /** Nose inner radius */
  double  fNoseRl;
  /** Nose outer radius */
  double  fNoseRh;
  /** Nose depth */
  double  fNoseDz;
  /** Nose start position */
  double  fNoseZ;
  /** Length of whole support structure */
  double  fConeL;
  /** Inner radius of back ring */
  double  fBackRl;
  /** Outer radius of back ring */
  double  fBackRh;
  /** Thickness of back ring */
  double  fBackDz;
  /** Thickness of beams */
  double  fBeamDz;
  /** Width of beams */
  double  fBeamW;
  /** Ending radius of flanges */
  double  fFlangeR;
};

//____________________________________________________________________
/** @brief Create a node geometry 
    @ingroup FMD_node_geom
    @code 
    .x NodeGeometry.C++
    @endcode 
 */
void
NodeGeometry() 
{
  TGeometry* geometry = new TGeometry("geometry","geometry");
  TShape* topShape = new TBRIK("topShape", "topShape", "", 40, 40, 150);
  TNode* topNode = new TNode("topNode", "topNode", topShape, 0, 0, 0, 0);
  topNode->SetVisibility(0);
  topNode->cd();
  
  Ring inner( 4.3, 17.2, 18, 13.4 / 2, .03, 1);
  Ring outer(15.6, 28.0,  9, 13.4 / 2, .03, 1);
  FMD3 fmd3(&inner, &outer);
  fmd3.CreateRings();
  
  TCanvas* c = new TCanvas("c", "c", 800, 800);
  c->SetFillColor(1);
  geometry->Draw();
  // c->x3d("ogl");
}
//____________________________________________________________________
//
// EOF
//
