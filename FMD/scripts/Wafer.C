//____________________________________________________________________
//
// $Id$
//
// Small script that I used to make some intial testing of the wafer
// layout and geometry. 
//
// Christian 
// 
/** Calculate wafer parameters
    @param c 
    @return  
    @ingroup FMD_simple_script
*/
TObjArray*
WaferParameters(const char c)
{      
  double dl;
  double dh;
  double r     = 134 / 2;
  double theta;
  switch (c) {
  case 'i': 
    dl     = 43;
    dh     = 172;
    theta  = 18;
    break;
  case 'o':
    dl     = 156;
    dh     = 280;
    theta  = 9;
    break;
  default:
    cerr << "Unknown wafer type: " << c << endl;
    return;
  }
  

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

  cout << "A: (" << dl << "," << y_A << ")" << endl;
  cout << "B: (" << dh << "," << y_B << ")" << endl;
  cout << "C: (" << x_C << "," << y_C << ")" << endl;
  cout << "D: (" << x_D << ",0)" << endl;
  
  cout << "Recentred at D:"  << endl;
  cout << "A: (" << dl - x_D  << "," << y_A << ")" << endl;
  cout << "B: (" << dh - x_D  << "," << y_B << ")" << endl;
  cout << "C: (" << x_C - x_D << "," << y_C << ")" << endl;

  TObjArray* verticies = new TObjArray(6);
  verticies->AddAt(new TVector2(dl,   y_A), 0);
  verticies->AddAt(new TVector2(x_C,  y_C), 1);
  verticies->AddAt(new TVector2(dh,   y_B), 2);
  verticies->AddAt(new TVector2(dh,  -y_B), 3);
  verticies->AddAt(new TVector2(x_C, -y_C), 4);
  verticies->AddAt(new TVector2(dl,  -y_A), 5);
  
  return verticies;
}

/** Draw wafers
    @ingroup FMD_simple_script
 */
void
Wafer()
{
  TCanvas* can = new TCanvas("can", "c", 400, 600);
  can->SetBorderMode(0);
  can->SetFillColor(0);
  
  TGeometry* g = new TGeometry("g", "G");
  TShape* topShape  = new TBRIK("top", "top", "", 100, 100, 100);
  TNode*  topNode = new TNode("top", "top", "top", 0, 0, 0);
  topNode->SetLineWidth(0);
  topNode->SetVisibility(0);

  TShape* virtualShape = new TTUBS("virtual", "Virtual", "", 
				   43, 172, 1, -18, 18);
  
  TObjArray* v = WaferParameters('i');
  TXTRU* moduleShape = new TXTRU("module", "module", "", 6, 2);
  for (Int_t i = 0; i  < 6; i++) {
    TVector2* vv = static_cast<TVector2*>(v->At(i));
    moduleShape->DefineVertex(i, vv->X(), vv->Y());
  }
  moduleShape->DefineSection(0, -1, 1, 0, 0);
  moduleShape->DefineSection(1,  1, 1, 0, 0);

  for (Int_t i = 0; i < 10; i++) {
    topNode->cd();
    Double_t theta   = 36 * i;
    Double_t z = (i % 2 ? +5 : -5);
    TRotMatrix* rot = new TRotMatrix(Form("rotation%02d", i), "Rotation", 
				     90, theta, 90, 
				     fmod(90 + theta, 360), 0, 0);
    TNode* moduleNode = new TNode(Form("module%02d", i), 
				  "Module", moduleShape, 0, 0, z,
				  rot);
    if (i == 0) {
      moduleNode->SetFillColor(2);
      moduleNode->SetLineColor(2);
      moduleNode->SetLineWidth(2);
    }
  }
  g->Draw();
  TView* view = can->GetView();
  view->SetPerspective();
  Int_t a;
  view->SetView(1.81208, 66.6725, 90, a);
  view->Zoom();
  view->Zoom();  
  view->Zoom();  
  can->Modified();
  can->cd();

  can->Print("fmd_module1.gif");
  // std::cout << "Waiting ..." << std::flush;
  // Char_t c = std::cin.get();
  
  topNode->cd();
  TNode* virtualNode = new TNode("virtual", "Virtual", 
				 virtualShape, 0, 0, -5);
  virtualNode->SetLineColor(3);
  virtualNode->SetLineWidth(2);
  virtualNode->SetLineStyle(2);
  g->Draw();
  view->SetPerspective();
  view->SetView(1.81208, 66.6725, 90, a);
  view->Zoom();
  view->Zoom();  
  view->Zoom();  
  can->Modified();
  can->cd();
  can->Print("fmd_module2.gif");
  
}
//____________________________________________________________________
//
// EOF
//





  
