UShort_t Encode(UShort_t d, Char_t r)
{
  UShort_t q = (r == 'I' || r == 'i') ? 0 : 1;
  UShort_t c = 1 << (d-1);
  UShort_t t = 1 << (c+q-1);

  return t;
 // return (1 << q) | (1 << (d+1));

}
void Decode(UShort_t bits, UShort_t& d, Char_t& r)
{
  if      (bits & 0x1)  { d = 1; r = 'I'; }
  else if (bits & 0x6)  { d = 2; r = (bits & 0x2 ? 'I' : 'O'); }
  else if (bits & 0x18) { d = 3; r = (bits & 0x8 ? 'I' : 'O'); }
}

const Char_t* ShowBits(UShort_t bits, Char_t* buf)
{
  // Char_t buf[7];
  // for (Int_t i = 0; i < 6; i++) buf[i] = ' ';
  for (Int_t i = 0; i < 6; i++) {
    buf[5-i] = (bits & (1 << i)) ? '1' : '0';
  }
  buf[6] = '\0';
  return buf;
}

void TestOne(UShort_t d, Char_t r)
{
  UShort_t bits = Encode(d, r);
  UShort_t rd   = 0;
  Char_t   rr   = 0;
  Char_t   buf[7];
  Decode(bits, rd, rr);
  ShowBits(bits, buf);
  
  Printf("FMD%d%c -> 0x%02x (%s) -> FMD%d%c", d, r, bits, buf, rd, rr);
}

enum {/*
  kFMD1i = 0x05,
  kFMD1  = kFMD1i,
  kFMD2i = 0x09,
  kFMD2o = 0x0a,
  kFMD2  = kFMD2i|kFMD2o,
  kFMD3i = 0x11,
  kFMD3o = 0x12,
  kFMD3  = kFMD3i|kFMD3o*/
  kFMD1I=0x01,
  kFMD1 =kFMD1I,
  kFMD2I=0x02,
  kFMD2O=0x04,
  kFMD2 =kFMD2I|kFMD2O,
  kFMD3I=0x08,
  kFMD3O=0x10,
  kFMD3 =kFMD3I|kFMD3O
};

UShort_t T(UShort_t m, UShort_t t)
{
  return (m & t) == t;
}
void TestEnum(UShort_t e, const char* n)
{
  Printf(" %10s | %5x | %5x | %5x | %5x | %5x |", 
	 n, T(e,kFMD1I), T(e,kFMD2I), T(e,kFMD2O), T(e,kFMD3I), T(e,kFMD3O));
}
void TestEnums()
{
  Printf(" Enum       | FMD1i | FMD2i | FMD2o | FMD3i | FMD3o |");
  TestEnum(kFMD1I, "FMD1i");
  TestEnum(kFMD1,  "FMD1");
  TestEnum(kFMD2I, "FMD2i");
  TestEnum(kFMD2O, "FMD2o");
  TestEnum(kFMD2,  "FMD2");
  TestEnum(kFMD3I, "FMD3i");
  TestEnum(kFMD3O, "FMD3o");
  TestEnum(kFMD3,  "FMD3");
  TestEnum(kFMD3|kFMD2,  "FMD23");
  TestEnum(kFMD1|kFMD2,  "FMD12");
  TestEnum(kFMD3|kFMD2O,  "FMD32o");
  TestEnum(kFMD3I|kFMD2O,  "FMD3i2o");
  TestEnum(kFMD1I|kFMD2O,  "FMD1i2o");
  TestEnum(kFMD1I|kFMD2O|kFMD3I,  "FMD3i2o1i");
  TestEnum(kFMD1I|kFMD2I|kFMD3I,  "FMD3i2i1i");
  TestEnum(kFMD2O|kFMD3O,  "FMD3o2o");
  TestEnum(0xff,  "All");
}

void RingBits()
{
  TestOne(1, 'I');
  TestOne(2, 'I');
  TestOne(2, 'O');
  TestOne(3, 'I');
  TestOne(3, 'O');

  TestEnums();
}
