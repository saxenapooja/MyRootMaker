#include "Common/GenTools.h"

#include <cassert>

namespace gen
{

GenParticle GetStableParticle(const GenParticle& gen)
{
  GenParticle Gen = gen;
  while(Gen.Status() != 1)
  {
    bool Found = false;
    for(unsigned int i = 0; i < Gen.NumDaughters(); ++i)
    {
      const GenParticle daughter = Gen.GetDaughter(i);
      if(daughter.PDGId() == Gen.PDGId())
      {
        Found = true;
        Gen = daughter;
        break;
      }
    }
    assert(Found);
  }
  return Gen;
}

GenParticle GetTauVisibleDecayProduct(const GenParticle& tau)
{
  assert(TMath::Abs(tau.PDGId()) == 15);

  // Skip internal nodes
  GenParticle Gen = tau;
  bool haveSubTau;
  do
  {
    haveSubTau = false;
    for(unsigned int i = 0; i < Gen.NumDaughters(); ++i)
    {
      const GenParticle daughter = Gen.GetDaughter(i);
      if(daughter.PDGId() == Gen.PDGId())
      {
        Gen = daughter;
        haveSubTau = true;
        break;
      }
    }
  } while(haveSubTau);

#if 1
  // Find non-neutrino decay product
  assert(TMath::Abs(Gen.PDGId()) == 15);
  GenParticle VisTau; bool haveVisTau = false;
  for(unsigned int i = 0; i < Gen.NumDaughters(); ++i)
  {
    const GenParticle daughter = Gen.GetDaughter(i);
    if(TMath::Abs(daughter.PDGId()) != 12 && TMath::Abs(daughter.PDGId()) != 14 && TMath::Abs(daughter.PDGId()) != 16)
    {
      // We can have a decay like tau -> e + gamma, so ignore all gammas. They are visible in principle,
      // but we are after the lepton.
      if(TMath::Abs(daughter.PDGId()) != 22)
      {
        assert(!haveVisTau);
        VisTau = daughter;
        haveVisTau = true;
      }
    }
  }

  assert(haveVisTau);
  Gen = VisTau;
#endif

  return Gen;
}

} // namespace gen
