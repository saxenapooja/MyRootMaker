
#ifndef _WHANALYSIS_SCALEFACTORS_H_
#define _WHANALYSIS_SCALEFACTORS_H_

class ScaleFactors
{
public:
  ScaleFactors(const char* data_era);

  float getDoubleMuTrigWeight(const char* trigName, float muPt1, float muEta1, float muPt2, float muEta2) const;
  float getMuIdIsoWeight(float muPt, float muEta) const;

private:
  enum Era {
    Era2011,
    Era2012,
    Data
  };

  Era era_;
};

#endif // _WHANALYSIS_SCALEFACTORS_H_
