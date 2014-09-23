
#ifndef _WHANALYSIS_PILEUPWEIGHT_H_
#define _WHANALYSIS_PILEUPWEIGHT_H_

#include <TFile.h>
#include <TH1D.h>

class PileupWeight
{
public:
  PileupWeight(const char* data_era);

  float getWeight(float trueNumInteractions) const;

private:
  TFile* pileupDataFile;
  TFile* pileupMCFile;
  TH1D* pileupData;
  TH1D* pileupMC;

  TH1D* pileupWeight;
};

#endif // _WHANALYSIS_PILEUPWEIGHT_H_
