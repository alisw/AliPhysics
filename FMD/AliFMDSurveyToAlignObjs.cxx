#include "AliFMDSurveyToAlignObjs.h"
#include "AliSurveyPoint.h"
#include <TGraph2DErrors.h>
#include <TF2.h>
#include <TVector3.h>
#include <iostream>
#include <iomanip>
#include <TMath.h>
#include <TRotation.h>
#include <TGeoMatrix.h>
#include <TGeoManager.h>
#include <TGeoPhysicalNode.h>
#include "AliFMDGeometry.h"

//____________________________________________________________________
Double_t
AliFMDSurveyToAlignObjs::GetUnitFactor() const
{
  // Returns the conversion factor from the measured values to
  // centimeters. 
  TString units(fSurveyObj->GetUnits());
  if      (units.CompareTo("mm", TString::kIgnoreCase) == 0) return .1;
  else if (units.CompareTo("cm", TString::kIgnoreCase) == 0) return 1.;
  else if (units.CompareTo("m",  TString::kIgnoreCase) == 0) return 100.;
  return 1;
}

//____________________________________________________________________
Bool_t
AliFMDSurveyToAlignObjs::GetPoint(const char* name, 
				  TVector3&   point, 
				  TVector3&   error) const
{
  // Get named point.   On return, point will contain the point
  // coordinates in centimeters, and error will contain the
  // meassurement errors in centimeters too.  If no point is found,
  // returns false, otherwise true. 
  Double_t unit = GetUnitFactor();
  TObject* obj  = fSurveyPoints->FindObject(name);
  if (!obj) return kFALSE;
  
   AliSurveyPoint* p = static_cast<AliSurveyPoint*>(obj);
   point.SetXYZ(unit * p->GetX(), 
		unit * p->GetY(),
		unit * p->GetZ());
   error.SetXYZ(unit * p->GetPrecisionX(),
		unit * p->GetPrecisionY(),
		unit * p->GetPrecisionZ());
   return kTRUE;
}

//____________________________________________________________________
Bool_t 
AliFMDSurveyToAlignObjs::CalculatePlane(const     TVector3& a, 
					const     TVector3& b,
					const     TVector3& c, 
					Double_t* trans,
					Double_t* rot) const
{
  // Vector a->b, b->c, and normal to plane defined by these two
  // vectors. 
  TVector3 ab(b-a), bc(c-b);
  
  // Normal vector to the plane of the fiducial marks obtained
  // as cross product of the two vectors on the plane d0^d1
  TVector3 nn(ab.Cross(bc));
  if (nn.Mag() < 1e-8) { 
    Info("DoIt", "Normal vector is null vector");
    return kFALSE;
  }
  
  // We express the plane in Hessian normal form.
  //
  //   n x = -p,
  //
  // where n is the normalised normal vector given by 
  // 
  //   n_x = a / l,   n_y = b / l,   n_z = c / l,  p = d / l
  //
  // with l = sqrt(a^2+b^2+c^2) and a, b, c, and d are from the
  // normal plane equation 
  //
  //  ax + by + cz + d = 0
  // 
  // Normalize
  TVector3 n(nn.Unit());
  Double_t p = - (n * a);
  
  // The center of the square with the fiducial marks as the
  // corners.  The mid-point of one diagonal - md.  Used to get the
  // center of the surveyd box. 
  TVector3 md(a + c);
  md *= 1/2.;
  
  // The center of the box. 
  TVector3 orig(md);
  trans[0] = orig[0];
  trans[1] = orig[1];
  trans[2] = orig[2];
  
  // Normalize the spanning vectors 
  TVector3 uab(ab.Unit());
  TVector3 ubc(bc.Unit());
  
  for (size_t i = 0; i < 3; i++) { 
    rot[i * 3 + 0] = uab[i];
    rot[i * 3 + 1] = ubc[i];
    rot[i * 3 + 2] = n[i];
  }
  return kTRUE;
}


//____________________________________________________________________
Bool_t
AliFMDSurveyToAlignObjs::GetFMD1Plane(Double_t* rot, Double_t* trans) const
{

  // The possile survey points 
  TVector3  icb, ict, ocb, oct, dummy;
  Int_t     missing = 0;
  if (!GetPoint("V0L_ICB", icb, dummy)) missing++;
  if (!GetPoint("V0L_ICT", icb, dummy)) missing++;
  if (!GetPoint("V0L_OCB", icb, dummy)) missing++;
  if (!GetPoint("V0L_OCT", icb, dummy)) missing++;

  // Check that we have enough points
  if (missing > 1) { 
    Error("GetFMD1Plane", "Only got %d survey points - no good for FMD1 plane",
	  4-missing);
    return kFALSE;
  }

  if (!CalculatePlane(icb, ict, oct, rot, trans)) return kFALSE;

  return kTRUE;
}

//____________________________________________________________________
Bool_t
AliFMDSurveyToAlignObjs::DoFMD1()
{
  return kTRUE;
}

//____________________________________________________________________
Bool_t
AliFMDSurveyToAlignObjs::GetFMD2Plane(Double_t* rot, Double_t* trans) const
{

  // The possile survey points 
  const char*    names[] = { "FMD2_ITOP",  "FMD2_OTOP", 
			     "FMD2_IBOTM", "FMD2_OBOTM", 
			     "FMD2_IBOT",  "FMD2_OBOT", 
			     0 };
  const char**   name    = names;
  Int_t          i       = 0;
  TGraph2DErrors g;
  // Loop and fill graph 
  while (*name) {
    TVector3 p, e;
    if (!GetPoint(*name++, p, e)) continue;

    g.SetPoint(i, p.X(), p.Y(), p.Z());
    g.SetPointError(i, e.X(), e.Y(), e.Z());
    i++;
  }
  
  // Check that we have enough points
  if (g.GetN() < 4) { 
    Error("GetFMD2Plane", "Only got %d survey points - no good for FMD2 plane",
	  g.GetN());
    return kFALSE;
  }

  // Next, declare fitting function and fit to graph. 
  // Fit to the plane equation: 
  // 
  //   ax + by + cz + d = 0
  //
  // or 
  // 
  //   z = - ax/c - by/c - d/c
  //
  TF2 f("fmd2Plane", "-[0]*x-[1]*y-[2]", 
	g.GetXmin(), g.GetXmax(), g.GetYmin(), g.GetYmax());
  g.Fit(&f, "Q");

  // Now, extract the normal and offset
  TVector3 nv(f.GetParameter(0), f.GetParameter(1), 1);
  TVector3 n(nv.Unit());
  Double_t p = -f.GetParameter(2);

  // Create two vectors spanning the plane 
  TVector3 a(1, 0, f.Eval(1, 0)-p);
  TVector3 b(0, -1, f.Eval(0, -1)-p);
  TVector3 ua(a.Unit());
  TVector3 ub(b.Unit());
  Double_t angAb = ua.Angle(ub);
  // PrintVector("ua: ", ua);
  // PrintVector("ub: ", ub);
  // std::cout << "Angle: " << angAb * 180 / TMath::Pi() << std::endl;
  
  for (size_t i = 0; i < 3; i++) { 
    rot[i * 3 + 0] = ua[i];
    rot[i * 3 + 1] = ub[i];
    rot[i * 3 + 2] = n[i];
  }
  
  // The intersection of the plane is given by (0, 0, -d/c)
  trans[0] = 0;
  trans[1] = 0;
  trans[2] = p;

  return kTRUE;
}

#define M(I,J) rot[(J-1) * 3 + (I-1)]
//____________________________________________________________________
Bool_t
AliFMDSurveyToAlignObjs::DoFMD2()
{
  Double_t rot[9], trans[3];
  if (!GetFMD2Plane(rot, trans)) return kFALSE;
  
  PrintRotation("Rotation: ", rot);
  PrintVector("Translation: ", trans);

#if 0
  Double_t theta = TMath::ATan2(M(3,1), M(3,2));
  Double_t phi   = TMath::ACos(M(3,3));
  Double_t psi   = -TMath::ATan2(M(1,3), M(2,3));
  
  std::cout << " Theta=" << theta * 180 / TMath::Pi() 
	    << " Phi="   << phi   * 180 / TMath::Pi() 
	    << " Psi="   << psi   * 180 / TMath::Pi() 
	    << std::endl;

  TRotation r;
  r.SetXEulerAngles(theta, phi, psi);
  r.Print();
  
  TGeoRotation*   geoR = new TGeoRotation("FMD2_survey_rotation", 
					  r.GetYTheta(),
					  r.GetYPhi(),
					  r.GetYPsi());
  TGeoCombiTrans* geoM = new TGeoCombiTrans(trans[0], trans[1], trans[2], geoR);
#else
  TGeoHMatrix* geoM = new TGeoHMatrix;
  geoM->SetRotation(rot);
  geoM->SetTranslation(trans);
#endif
  geoM->Print();

  const char* path = "/ALIC_1/F2MT_2/FMD2_support_0/FMD2_back_cover_2";
  gGeoManager->cd(path);
  TGeoMatrix* globalBack = gGeoManager->GetCurrentMatrix();
  globalBack->Print();
  PrintRotation("Back rotation:",    globalBack->GetRotationMatrix());
  PrintVector("Back translation:",   globalBack->GetTranslation());
  
  // TObjArray* pns = gGeoManager->GetListOfPhysicalNodes();
  // TObject*   oT  = pns->FindObject("/ALIC_1/F2MT_2");
  // TObject*   oB  = pns->FindObject("/ALIC_1/F2MB_2");
  // if (!oT) { 
  //   Warning("DoFMD2", Form("Physical node /ALIC_1/F2MT_2 not found"));
  //   return kFALSE;
  // }
  // if (!oB) { 
  //   Warning("DoFMD2", Form("Physical node /ALIC_1/F2MB_2 not found"));
  //   return kFALSE;
  // }
  // TGeoPhysicalNode* top = static_cast<TGeoPhysicalNode*>(oT);
  // TGeoPhysicalNode* bot = static_cast<TGeoPhysicalNode*>(oB);
  TGeoHMatrix tDelta(globalBack->Inverse());
  TGeoHMatrix bDelta(globalBack->Inverse());
  PrintRotation("Back^-1 rotation:",    tDelta.GetRotationMatrix());
  PrintVector("Back^-1 translation:",   tDelta.GetTranslation());

  std::cout << "tDelta = 1? " << tDelta.IsIdentity() << std::endl;
  
  tDelta.MultiplyLeft(geoM);
  bDelta.MultiplyLeft(geoM);
  
  PrintRotation("tDelta rotation:",  tDelta.GetRotationMatrix());
  PrintVector("tDelta translation:", tDelta.GetTranslation());
  PrintRotation("bDelta rotation:",  bDelta.GetRotationMatrix());
  PrintVector("bDelta translation:", bDelta.GetTranslation());
  
  return kTRUE;
}

//____________________________________________________________________
void
AliFMDSurveyToAlignObjs::Run()
{
  AliFMDGeometry* geom = AliFMDGeometry::Instance();
  geom->Init();
  geom->InitTransformations();
  
  DoFMD2();
}

//____________________________________________________________________
void 
AliFMDSurveyToAlignObjs::PrintVector(const char* text, const TVector3& v)
{
  Double_t va[] = { v.X(), v.Y(), v.Z() };
  PrintVector(text, va);
}
//____________________________________________________________________
void 
AliFMDSurveyToAlignObjs::PrintVector(const char* text, const Double_t* v)
{
  std::cout << text 
	    << std::setw(15) << v[0] 
	    << std::setw(15) << v[1]
	    << std::setw(15) << v[2]
	    << std::endl;
}


//____________________________________________________________________
void 
AliFMDSurveyToAlignObjs::PrintRotation(const char* text, const Double_t* rot)
{
  std::cout << text << std::endl;
  for (size_t i = 0; i < 3; i++) { 
    for (size_t j = 0; j < 3; j++) 
      std::cout << std::setw(15) << rot[i * 3 + j];
    std::cout << std::endl;
  }
}

//____________________________________________________________________
//
// EOF
//
