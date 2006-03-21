//____________________________________________________________________
//
// $Id$
//
// Test of AliFMDAltro{Reader,Writer}
//
void
TestRawIO()
{
  std::ofstream ofile("foo.dat");
  AliFMDAltroWriter w(ofile);
  for (size_t i = 0; i < 16; i++) 
    w.AddSignal((i << 4)  + i);
  w.AddChannelTrailer(0xabe);
  w.Close();
  ofile.close();

  std::ifstream ifile("foo.dat");
  AliRawDataHeader h;
  ifile.read((char*)&h, sizeof(h));
  AliFMDAltroReader r(ifile);
  UShort_t hwaddr, last;
  UShort_t data[1024];
  int ret = r.ReadChannel(hwaddr, last, data);
  printf("Read returned %d, w/addr=0x%x, last=%d\n", ret, hwaddr, last);
  if (ret < 0) return;
  for (size_t i = 0; i < last; i++)
    std::cout << i << "\t0x" << std::hex << data[i] << std::endl;
}
//____________________________________________________________________
//
// EOF
//
