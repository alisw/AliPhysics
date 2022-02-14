#ifndef PLOTTING_H
#define PLOTTING_H

#include <array>

#include <TCanvas.h>
#include <TColor.h>
#include <TH1.h>
#include <TStyle.h>

namespace plotting {
  constexpr std::array<int,11> kHighContrastColors{
    kRed,
    kOrange,
    kGreen+1,
    kAzure+1,
    kBlue,
    kViolet+1,
    kBlack,
    kGray+1,
    kOrange-7,
    kOrange+4,
    kYellow -2
  };

  constexpr std::array<int,10> kListColors{
    kBlack,
    kRed,
    kOrange,
    kGreen+3,
    kAzure+1,
    kRed,
    kOrange,
    kRed,
    kOrange+3,
    kGreen+3
  };

  const std::array<int,11> kSpectraColors{
    TColor::GetColor("#ff3300"),
    TColor::GetColor("#ec6e0a"),
    TColor::GetColor("#daaa14"),
    TColor::GetColor("#c7e51e"),
    TColor::GetColor("#85dd69"),
    TColor::GetColor("#42d6b4"),
    TColor::GetColor("#00ceff"),
    TColor::GetColor("#009adf"),
    TColor::GetColor("#0067c0"),
    TColor::GetColor("#595959"),
    TColor::GetColor("#0033a1")
  };

  const Style_t kChargeMarker[2] = {20,24};

  void SetHistStyle(TH1* h, int color, int marker = 20, const char* opt = "", int linew = 1, int fillstyle = 0) {
    h->SetMarkerColor(color);
    h->SetLineColor(color);
    h->SetFillStyle(fillstyle);
    h->SetLineWidth(linew);
    h->SetMarkerStyle(marker);
    h->SetOption(opt);
  }


  void CanvasPartition(TCanvas *C,TPad* pad[], const Int_t Nx = 2,const Int_t Ny = 2,Float_t lMargin = 0.05, Float_t rMargin = 0.05,Float_t bMargin = 0.10, Float_t tMargin = 0.05) {
    if (!C) return;
     // Setup Pad layout:
    double vSpacing = 0.0;
    double vStep  = (1.- bMargin - tMargin - (Ny-1) * vSpacing) / Ny;
    double hSpacing = 0.0;
    double hStep  = (1.- lMargin - rMargin - (Nx-1) * hSpacing) / Nx;
    double vposd = 0.,vposu = 0.,vmard = 0.,vmaru = 0.,vfactor = 0.;
    double hposl = 0.,hposr = 0.,hmarl = 0.,hmarr = 0.,hfactor = 0.;
    for (Int_t i=0;i<Nx;i++) {
      if (i==0) {
        hposl = 0.0;
        hposr = lMargin + hStep;
        hfactor = hposr-hposl;
        hmarl = lMargin / hfactor;
        hmarr = 0.0;
      }
      else if (i == Nx-1) {
        hposl = hposr + hSpacing;
        hposr = hposl + hStep + rMargin;
        hfactor = hposr-hposl;
        hmarl = 0.0;
        hmarr = rMargin / (hposr-hposl);
      }
      else {
        hposl = hposr + hSpacing;
        hposr = hposl + hStep;
        hfactor = hposr-hposl;
        hmarl = 0.0;
        hmarr = 0.0;
      }
      for (Int_t j=0;j<Ny;j++) {
        if (j==0) {
          vposd = 0.0;
          vposu = bMargin + vStep;
          vfactor = vposu-vposd;
          vmard = bMargin / vfactor;
          vmaru = 0.0;
        } else if (j == Ny-1) {
          vposd = vposu + vSpacing;
          vposu = vposd + vStep + tMargin;
          vfactor = vposu-vposd;
          vmard = 0.0;
          vmaru = tMargin / (vposu-vposd);
        } else {
          vposd = vposu + vSpacing;
          vposu = vposd + vStep;
          vfactor = vposu-vposd;
          vmard = 0.0;
          vmaru = 0.0;
        }
        C->cd(0);
        pad[i*Nx+Ny-1-j] = new TPad(Form("pad_%d_%d",i,j),"",hposl,vposd,hposr,vposu);
        pad[i*Nx+Ny-1-j]->SetLeftMargin(hmarl);
        pad[i*Nx+Ny-1-j]->SetRightMargin(hmarr);
        pad[i*Nx+Ny-1-j]->SetBottomMargin(vmard);
        pad[i*Nx+Ny-1-j]->SetTopMargin(vmaru);
        pad[i*Nx+Ny-1-j]->SetFrameBorderMode(0);
        pad[i*Nx+Ny-1-j]->SetBorderMode(0);
        pad[i*Nx+Ny-1-j]->SetBorderSize(0);
        pad[i*Nx+Ny-1-j]->Draw();
      }
    }
  }

}

#endif
