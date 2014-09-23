
#ifndef DESYHTAUTAU_TAUTAUANALYSIS_INTERFACE_GENTOOLS_H
#define DESYHTAUTAU_TAUTAUANALYSIS_INTERFACE_GENTOOLS_H

#include "AnalysisTool/DataForm.h"

namespace gen
{

GenParticle GetStableParticle(const GenParticle& gen);
GenParticle GetTauVisibleDecayProduct(const GenParticle& tau);

}

#endif // DESYHTAUTAU_TAUTAUANALYSIS_INTERFACE_GENTOOLS_H
