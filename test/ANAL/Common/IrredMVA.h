
#ifndef _WHANALYSIS_IRREDMVA_H_
#define _WHANALYSIS_IRREDMVA_H_

#include <TMVA/Reader.h>

class IrredMVA
{
public:
  IrredMVA();

  double mvaValue(float muPt, float tauOSPt, float tauSSPt, float met) const;

protected:
  TMVA::Reader* reader;

  mutable float muPt;
  mutable float tauOSPt; 
  mutable float tauSSPt;
  mutable float met;
};

#endif // _WHANALYSIS_IRREDMVA_H_
