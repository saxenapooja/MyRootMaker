
#ifndef _WHANALYSIS_COMBINEDMVA_H_
#define _WHANALYSIS_COMBINEDMVA_H_

#include <TMVA/Reader.h>

class CombinedMVA
{
public:
  CombinedMVA();

  double mvaValue(float BDT8Output, float MLPOutput) const;

protected:
  TMVA::Reader* reader;

  mutable float BDT8Output;
  mutable float MLPOutput;
};

#endif // _WHANALYSIS_COMBINEDMVA_H_
