
#ifndef _WHANALYSIS_MVAEVALUATOR_H_
#define _WHANALYSIS_MVAEVALUATOR_H_

#include "Common/MvaTree.h"

#include <TMVA/Reader.h>

class MvaEvaluator
{
public:
  MvaEvaluator(const char* method, const char* training, const char* variables);

  bool isAvailable() const;
  double mvaValue(const MvaTree::Variables::MvaVars& vars) const;

protected:
  const std::string m_method;
  TMVA::Reader* reader; 

  mutable MvaTree::Variables::MvaVars mvaVars;
  mutable float absDeltaPhiDiTau;
  mutable float absDeltaEtaDiTau;
};

#endif // _WHANALYSIS_MVAEVALUATOR_H_
