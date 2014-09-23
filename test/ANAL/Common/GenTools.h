#ifndef DESYHTAUTAU_TAUTAUANALYSIS_INTERFACE_GENTOOLS_H
#define DESYHTAUTAU_TAUTAUANALYSIS_INTERFACE_GENTOOLS_H

#include "AnalysisTool/DataForm.h"
#include <vector>
namespace gen
{
  GenParticle GetStableParticle(const GenParticle& gen);
  GenParticle GetTauVisibleDecayProduct(const GenParticle& tau);
  GenParticle GetTauPionDecayProduct(const GenParticle& tau);
/*   GenParticle GetPPionFromTauDecay(const GenParticle& tau); */
/*   GenParticle GetNPionFromTauDecay(const GenParticle& tau); */
  std::vector<GenParticle> GetGenTausFromAllHiggs(const GenParticle& higgs);
  GenParticle GetTauInvisibleDecayProduct(const GenParticle& tau);
}

#endif // DESYHTAUTAU_TAUTAUANALYSIS_INTERFACE_GENTOOLS_H
