//____________________________________________________________________
//
//
// $Id$
//
// Script I used for rapid prototyping of the FMD3 geometry - in
// particular the support cone 
//
//____________________________________________________________________
TObjArray*
waferParameters(double dl, double dh, double theta, double r, 
		Bool_t verbose=kFALSE)
{      
  double tan_theta  = tan(theta * TMath::Pi() / 180.);
  double tan_theta2 = pow(tan_theta,2);
  double r2         = pow(r,2);
  double y_A        = tan_theta * dl;
  double x_D        = dl + sqrt(r2 - tan_theta2 * pow(dl,2));
  double x_D_2      = dl - sqrt(r2 - tan_theta2 * pow(dl,2));
  
  double y_B       = sqrt(r2 - pow(dh,2) + 2 * dh * x_D - pow(x_D,2));
  double x_C       = (x_D + sqrt(-tan_theta2 * pow(x_D,2) + r2 
				 + r2 * tan_theta2)) / (1 + tan_theta2);
  double y_C       = tan_theta * x_C;

  if (verbose) {
    cout << "A: (" << dl << "," << y_A << ")" << endl;
    cout << "B: (" << dh << "," << y_B << ")" << endl;
    cout << "C: (" << x_C << "," << y_C << ")" << endl;
    cout << "D: (" << x_D << ",0)" << endl;
    
    cout << "Recentred at D:"  << endl;
    cout << "A: (" << dl - x_D  << "," << y_A << ")" << endl;
    cout << "B: (" << dh - x_D  << "," << y_B << ")" << endl;
    cout << "C: (" << x_C - x_D << "," << y_C << ")" << endl;
  }
  
  TObjArray* verticies = new TObjArray(6);
  verticies->AddAt(new TVector2(dl,   y_A), 0);
  verticies->AddAt(new TVector2(x_C,  y_C), 1);
  verticies->AddAt(new TVector2(dh,   y_B), 2);
  verticies->AddAt(new TVector2(dh,  -y_B), 3);
  verticies->AddAt(new TVector2(x_C, -y_C), 4);
  verticies->AddAt(new TVector2(dl,  -y_A), 5);
  
  return verticies;
}

//____________________________________________________________________
TShape* 
createModuleShape(const Char_t* name, double rl, double rh, double th, 
		  double r, double dz) 
{
  std::cout << "Creating Module shape for " << name << std::flush;
  // TShape* virtualShape = new TTUBS(Form("%sVirtual", name), 
  //                                  Form("Virtual %s", name), 
  //                                  "", rl, rh, 1, -th, th);
  
  TObjArray* v = waferParameters(rl, rh, th, r);
  TXTRU* moduleShape = new TXTRU(Form("%sModule", name),
				 Form("Module %s", name), 
				 "", 6, 2);
  for (Int_t i = 0; i  < 6; i++) {
    std::cout << "." << std::flush;
    TVector2* vv = static_cast<TVector2*>(v->At(i));
    moduleShape->DefineVertex(i, vv->X(), vv->Y());
  }
  moduleShape->DefineSection(0, -dz, 1, 0, 0);
  moduleShape->DefineSection(1,  dz, 1, 0, 0);
  std::cout << std::endl;

  return (TShape*)moduleShape;
}

//____________________________________________________________________
TNode* 
createRing(const char* name, double rl, double rh, double th, 
	   double siThick, double waferR, double staggering, double z) 
{
  std::cout << "Creating Ring node for " << name << std::flush;
  TShape* ringShape   = new TTUBE(Form("%sShape", name), "Ring Shape", 
				  "", rl, rh, staggering + siThick);
  TNode*  ringNode    = new TNode(Form("%sNode", name), "Ring Node", 
				  ringShape, 0, 0, z, 0);
  TShape* moduleShape = createModuleShape(name, rl, rh, th, waferR, siThick);
  Int_t n = 360 / 2 / th;
  for (Int_t i = 0; i < n; i++) {
    std::cout << "." << std::flush;
    ringNode->cd();
    Double_t theta = 2  * th * i;
    Double_t z     = (i % 2 ? 0 : staggering);
    TRotMatrix* rot = new TRotMatrix(Form("%sRotation%02d", name, i), 
				     "Rotation", 90, theta, 90, 
				     fmod(90 + theta, 360), 0, 0);
    TNode* moduleNode = new TNode(Form("%sModule%02d", name, i), 
				  "Module", moduleShape, 0, 0, z,
				  rot);
    moduleNode->SetFillColor(5);
    moduleNode->SetLineColor(5);
    moduleNode->SetLineWidth(2);
  }
  std::cout << std::endl;
  ringNode->SetVisibility(0);
  return ringNode;
}

//____________________________________________________________________
TNode*
createSupport(double noseRl, double noseRh, double noseDz, double noseZ, 
	      double backRl, double backRh, double backDz, double coneL, 
	      double beamW, double beamDz,  double flangeR) 
{
  TShape* noseShape = new TTUBE("noseShape", "Nose Shape", "", 
				noseRl, noseRh, noseDz);
  TNode*  noseNode  = new TNode("noseNode", "noseNode", noseShape, 
				0, 0, noseZ - noseDz, 0);
  noseNode->SetLineColor(0);
  
  Double_t zdist = coneL - 2 * backDz - 2 * noseDz;
  Double_t tdist = backRh - noseRh;
  Double_t beamL = TMath::Sqrt(zdist * zdist + tdist * tdist);
  Double_t theta = -TMath::ATan2(tdist, zdist);
  

  TShape* backShape = new TTUBE("backShape", "Back Shape", "", 
				 backRl, backRh, backDz);
  TNode*  backNode  = new TNode("backNode", "backNode", backShape, 
				0, 0, noseZ - 2 * noseDz - zdist - backDz, 0);
  backNode->SetLineColor(0);


  TShape* beamShape = new TBRIK("beamShape", "beamShape", "", 
				beamDz, beamW / 2 , beamL / 2);
  Int_t    n = 8;
  Double_t r = noseRl + tdist / 2;
  for (Int_t i = 0; i < n; i++) {
    Double_t phi   = 360. / n * i;
    Double_t t     = 180. * theta / TMath::Pi();
    TRotMatrix* beamRotation = new TRotMatrix(Form("beamRotation%d", i), 
						Form("beamRotation%d", i),
						180-t, phi, 90, 90+phi, 
						t, phi);
    TNode* beamNode = new TNode(Form("beamNode%d", i), 
				Form("beamNode%d", i), beamShape, 
				r * TMath::Cos(phi / 180 * TMath::Pi()), 
				r * TMath::Sin(phi / 180 * TMath::Pi()), 
				noseZ - 2 * noseDz - zdist / 2,  beamRotation);
    beamNode->SetLineColor(0);
  }
  
  Double_t flangel = (flangeR - backRh) / 2;
  TShape* flangeShape = new TBRIK("flangeShape", "FlangeShape", "", 
				  flangel, beamW / 2, backDz);
  n = 4;
  r = backRh + flangel;
  for (Int_t i = 0; i < n; i++) {
    Double_t phi = 360. / n * i + 180. / n;
    TRotMatrix* flangeRotation = new TRotMatrix(Form("flangeRotation%d", i),
						Form("Flange Rotation %d", i),
						90, phi, 90, 90+phi, 0, 0);
    TNode* flangeNode = new TNode(Form("flangeNode%d", i), 
				  Form("flangeNode%d", i), 
				  flangeShape,
				  r * TMath::Cos(phi / 180 * TMath::Pi()), 
				  r * TMath::Sin(phi / 180 * TMath::Pi()),
				  noseZ - 2 * noseDz - zdist - backDz, 
				  flangeRotation);
    flangeNode->SetLineColor(0);
				  
  }
  return 0;
}

				 


//____________________________________________________________________
void
SimpleGeometry() 
{
  TGeometry* geometry = new TGeometry("geometry","geometry");
  TShape* topShape = new TBRIK("topShape", "topShape", "", 40, 40, 150);
  TNode* topNode = new TNode("topNode", "topNode", topShape, 0, 0, 0, 0);
  topNode->SetVisibility(0);
  topNode->cd();
  
  Double_t waferR     = 13.4 / 2;
  Double_t siThick    = .03;
  Double_t staggering = 1;
  Double_t innerRh    = 17.2;
  Double_t innerRl    = 4.3;
  Double_t innerTh    = 18;
  Double_t innerZ     = -62.8;
  Double_t outerRh    = 28;
  Double_t outerRl    = 15.6;
  Double_t outerTh    = 9;
  Double_t outerZ     = -75.2;
  Double_t noseRl     = 5.5;
  Double_t noseRh     = 6.7;
  Double_t noseDz     = 2.8 / 2;
  Double_t noseZ      = -46;
  Double_t coneL      = 30.9;
  Double_t backRl     = 61 / 2;
  Double_t backRh     = 66.8 /2;
  Double_t backDz     = 1.4 / 2;
  Double_t beamDz     = .5 / 2;
  Double_t beamW      = 6;
  Double_t flangeR    = 49.25;

  Double_t zdist      = coneL - 2 * backDz - 2 * noseDz;
  Double_t tdist      = backRh - noseRh;
  Double_t alpha      = tdist / zdist;

  Double_t x, rl, rh, z;
  z  = noseZ - coneL / 2;
  TPCON* fmd3Shape = new TPCON("fmd3Shape", "FMD 3 Shape", "", 0, 360, 7);
  
  x  = noseZ;
  rl = noseRl;
  rh = noseRh;
  fmd3Shape->DefineSection(0, x - z, rl, rh);

  x  = noseZ-2*noseDz;
  fmd3Shape->DefineSection(1, x - z, rl, rh);

  x  = innerZ - staggering - siThick; 
  rl = innerRl;
  rh = noseRh + alpha * TMath::Abs(x-noseZ + 2 * noseDz);
  fmd3Shape->DefineSection(2, x - z, rl, rh);

  x  = outerZ;
  rl = outerRl;
  rh = backRh;
  fmd3Shape->DefineSection(3, x - z, rl, rh);

  x  = noseZ - zdist - 2 * noseDz;
  rl = outerRl;
  rh = backRh;
  fmd3Shape->DefineSection(4, x - z, rl, rh);

  x  = noseZ - zdist - 2 * noseDz;
  rl = outerRl;
  rh = flangeR;
  fmd3Shape->DefineSection(5, x - z, rl, rh);

  x  = noseZ - coneL;
  rl = outerRl;
  rh = flangeR;
  fmd3Shape->DefineSection(6, x - z, rl, rh);

  TNode* fmd3Node = new TNode("fmd3Node", "FMD3 Node", fmd3Shape, 
			      0, 0, z, 0);
  fmd3Node->SetLineColor(3);
  fmd3Node->SetVisibility(1);

  fmd3Node->cd();
  TNode* innerNode = createRing("inner", innerRl, innerRh, innerTh, 
				siThick, waferR, staggering, innerZ-z);


  fmd3Node->cd();
  TNode* outerNode = createRing("outer", outerRl, outerRh, outerTh, 
				siThick, waferR, staggering, outerZ-z);
  

  fmd3Node->cd();
  TNode* supportNode = createSupport(noseRl, noseRh, noseDz, noseZ-z, 
				     backRl, backRh, backDz, coneL,
				     beamW, beamDz, flangeR);
  
  TCanvas* c = new TCanvas("c", "c", 800, 800);
  c->SetFillColor(1);
  geometry->Draw();
  // c->x3d("ogl");
}
//____________________________________________________________________
//
// EOF
//
