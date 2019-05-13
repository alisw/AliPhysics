//**************************************************************************\
//* This file is property of and copyright by the ALICE Project            *\
//* ALICE Experiment at CERN, All rights reserved.                         *\
//*                                                                        *\
//* Primary Authors: Matthias Richter <Matthias.Richter@ift.uib.no>        *\
//*                  for The ALICE HLT Project.                            *\
//*                                                                        *\
//* Permission to use, copy, modify and distribute this software and its   *\
//* documentation strictly for non-commercial purposes is hereby granted   *\
//* without fee, provided that the above copyright notice appears in all   *\
//* copies and that both the copyright notice and this permission notice   *\
//* appear in the supporting documentation. The authors make no claims     *\
//* about the suitability of this software for any purpose. It is          *\
//* provided "as is" without express or implied warranty.                  *\
//**************************************************************************

/// \file bitmapfile.h
/// \author David Rohr

struct BITMAPFILEHEADER {
  unsigned short bfType;
  unsigned int bfSize;
  unsigned int bfReserved;
  unsigned int bfOffBits;
} __attribute__((packed));

struct BITMAPINFOHEADER {
  unsigned int biSize;
  unsigned int biWidth;
  unsigned int biHeight;
  unsigned short biPlanes;
  unsigned short biBitCount;
  unsigned int biCompression;
  unsigned int biSizeImage;
  unsigned int biXPelsPerMeter;
  unsigned int biYPelsPerMeter;
  unsigned int biClrUsed;
  unsigned int biClrImportant;
} __attribute__((packed));

enum BI_Compression { BI_RGB = 0x0000,
                      BI_RLE8 = 0x0001,
                      BI_RLE4 = 0x0002,
                      BI_BITFIELDS = 0x0003,
                      BI_JPEG = 0x0004,
                      BI_PNG = 0x0005,
                      BI_CMYK = 0x000B,
                      BI_CMYKRLE8 = 0x000C,
                      BI_CMYKRLE4 = 0x000D };
