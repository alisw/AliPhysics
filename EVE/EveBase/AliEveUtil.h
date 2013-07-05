#ifndef ALIEVEUTIL_H
#define ALIEVEUTIL_H

class TGPicturePool;

class AliEveUtil
{
public:
// used by a TGComboBox for selecting CDB Path
enum CDB_PATH_TYPE
{
      CDB_LOCAL,
      CDB_RAW,
      CDB_MCIDEAL,
      CDB_MCRESIDUAL,
      CDB_MCFULL
};

    static void Init(); // init utils - run this once before any call to any method
    static TGPicturePool* GetPicturePool();

private:
    AliEveUtil();
    ~AliEveUtil();

    AliEveUtil(const AliEveUtil& other); // Not Implemented
    AliEveUtil operator=(const AliEveUtil& other);  // Not Implemented

    static AliEveUtil* fgAliEveUtil;


    ClassDef(AliEveUtil, 0);
};

R__EXTERN TGPicturePool* gAliEvePicturePool;
#endif // ALIEVEUTIL_H
