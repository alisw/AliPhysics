void bright_colors(Float_t value=1)
{
  // Notes:
  // you can store the original colors by creating a clone of
  // (TObjArray*)gROOT->GetListOfColors() and restore the colors by
  // assigning the vector with original values to the list of colors
  // that gROOT handles.

  if (value > 5)
  {
    printf("Value %f too high - maximum is 5.\n", value);
    return;
  }
  value *= 0.1f;

  TObjArray *colors = (TObjArray*) gROOT->GetListOfColors();
  TColor    *color  = 0;
  Float_t    r, g, b;
  for (int i = 0; i < colors->GetSize(); ++i)
  {
    if ((color = dynamic_cast<TColor*>(colors->At(i))) != 0)
    {
      color->GetRGB(r, g, b);
      if (r < 0.01 && g < 0.01 && b < 0.01) continue; // skip black
      if (r > 0.95 && g > 0.95 && b > 0.95) continue; // skip white
      r = TMath::Min(r + value, 1.0f);
      g = TMath::Min(g + value, 1.0f);
      b = TMath::Min(b + value, 1.0f);
      color->SetRGB(r, g, b);
    }
  }

  if (gEve)
    gEve->FullRedraw3D();
}
