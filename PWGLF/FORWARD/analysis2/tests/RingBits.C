UShort_t Encode(UShort_t d, Char_t r)
{
  UShort_t q = (r == 'I' || r == 'i') ? 0 : 1;
  
  return (1 << q) | (1 << (d+1));
}
void Decode(UShort_t bits, UShort_t& d, Char_t& r)
{
  d          = (bits & (1<<2) ? 1 : 
		bits & (1<<3) ? 2 : 
		bits & (1<<4) ? 3 : 0);
  UShort_t q = (bits & 0x3);
  r          = q == 1 ? 'I' : 'O';
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

enum {
  kFMD1i = 0x05,
  kFMD1  = kFMD1i,
  kFMD2i = 0x09,
  kFMD2o = 0x0a,
  kFMD2  = kFMD2i|kFMD2o,
  kFMD3i = 0x11,
  kFMD3o = 0x12,
  kFMD3  = kFMD3i|kFMD3o
};

UShort_t T(UShort_t m, UShort_t t)
{
  return (m & t) == t;
}
void TestEnum(UShort_t e, const char* n)
{
  Printf(" %6s | %5x | %5x | %5x | %5x | %5x |", 
	 n, T(e,kFMD1i), T(e,kFMD2i), T(e,kFMD2o), T(e,kFMD3i), T(e,kFMD3o));
}
void TestEnums()
{
  Printf(" Enum   | FMD1i | FMD2i | FMD2o | FMD3i | FMD3o |");
  TestEnum(kFMD1i, "FMD1i");
  TestEnum(kFMD1,  "FMD1");
  TestEnum(kFMD2i, "FMD2i");
  TestEnum(kFMD2o, "FMD2o");
  TestEnum(kFMD2,  "FMD2");
  TestEnum(kFMD3i, "FMD3i");
  TestEnum(kFMD3o, "FMD3o");
  TestEnum(kFMD3,  "FMD3");
  TestEnum(kFMD3|kFMD2,  "FMD23");
  TestEnum(kFMD1|kFMD2,  "FMD12");
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
