
#ifndef _WHANALYSIS_METWEIGHT_H_
#define _WHANALYSIS_METWEIGHT_H_

#include <TFile.h>
#include <TF1.h>

class METWeight
{
public:
  METWeight(const char* data_era, const char* channel);

  double getWeight(double met) const;

private:
  TFile* fileMETWeight;
  TF1* fMETWeight;
};

#endif // _WHANALYSIS_METWEIGHT_H_
