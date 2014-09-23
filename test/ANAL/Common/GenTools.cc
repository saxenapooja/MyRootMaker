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
    do {
      haveSubTau = false;
      for (unsigned int i = 0; i < Gen.NumDaughters(); ++i) {
	const GenParticle daughter = Gen.GetDaughter(i);
	if (daughter.PDGId() == Gen.PDGId())
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
    GenParticle VisTau; 
    bool haveVisTau = false;
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
  

  
  GenParticle GetTauPionDecayProduct(const GenParticle& tau)
  {
    assert(TMath::Abs(tau.PDGId()) == 15);
    
    // Skip internal nodes
    GenParticle Gen = tau;
    bool haveSubTau;
    do {
      haveSubTau = false;
      for (unsigned int i = 0; i < Gen.NumDaughters(); ++i) {
	const GenParticle daughter = Gen.GetDaughter(i);
	if (daughter.PDGId() == Gen.PDGId())
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
    GenParticle VisTau; 
    bool haveVisTau = false;
    for(unsigned int i = 0; i < Gen.NumDaughters(); ++i)
      {
	const GenParticle daughter = Gen.GetDaughter(i);

	if(TMath::Abs(daughter.PDGId()) != 12 && TMath::Abs(daughter.PDGId()) != 14 && TMath::Abs(daughter.PDGId()) != 16)
	  {
	    // We can have a decay like tau -> e + gamma, so ignore all gammas. They are visible in principle,
	    // but we are after the lepton.
	    if(TMath::Abs(daughter.PDGId()) == 211 && daughter.Status() == 1)
	      {
		assert(!haveVisTau);
		VisTau = daughter;
		haveVisTau = true;
	      }
	  }
      }
    
    // assert(haveVisTau);
    if(haveVisTau) Gen = VisTau;
    
#endif
    
    return Gen;
  }


  GenParticle GetTauInvisibleDecayProduct(const GenParticle& tau)
  {
    assert(TMath::Abs(tau.PDGId()) == 15);
    
    // Skip internal nodes
    GenParticle Gen = tau;
    bool haveSubTau;
    do {
      haveSubTau = false;
      for (unsigned int i = 0; i < Gen.NumDaughters(); ++i) {
	const GenParticle daughter = Gen.GetDaughter(i);
	if (daughter.PDGId() == Gen.PDGId())
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
    GenParticle VisTau; 
    bool haveVisTau = false;
    for(unsigned int i = 0; i < Gen.NumDaughters(); ++i)
      {
	const GenParticle daughter = Gen.GetDaughter(i);

	    if(TMath::Abs(daughter.PDGId()) == 16 && daughter.Status() == 1)
	      {
		assert(!haveVisTau);
		VisTau = daughter;
		haveVisTau = true;
	      }
	    
      }
    
    // assert(haveVisTau);
    if(haveVisTau) Gen = VisTau;
    
#endif
    
    return Gen;
  }
  
    std::vector<GenParticle> GetGenTausFromAllHiggs(const GenParticle& higgs)
    {
    vector<GenParticle> allGenTaus;
    
    //make sure its higgs
    //    assert(TMath::Abs(higgs.PDGId()) == 23 || TMath::Abs(higgs.PDGId()) == 25  || TMath::Abs(higgs.PDGId()) == 35 || TMath::Abs(higgs.PDGId()) == 36 );
    
    // Skip internal nodes
    GenParticle Gen = higgs;
    bool haveTau;
    for(unsigned int i = 0; i < Gen.NumDaughters(); ++i)
      {
	const GenParticle daughter = Gen.GetDaughter(i);
	
	//cout<<"daughther charge :"<< daughter.Charge() <<" pdgid :"<< daughter.PDGId() << endl;
 	if(TMath::Abs(daughter.PDGId()) == 15 && 
	   ( daughter.Status() == 2 || daughter.Status() == 3))
	  {
	    allGenTaus.push_back(daughter);
	    //Gen = daughter;
	    haveTau= true;
	  }
      } // for(unsigned int i = 0; i < Gen.NumDaughters(); ++i)
    assert(haveTau);
    return allGenTaus;
    }
  
  
  
//   GenParticle GetPPionFromTauDecay(const GenParticle& tau)
//   {
//     assert(TMath::Abs(tau.PDGId()) == 15);
    
//     // Skip internal nodes
//     GenParticle Gen = tau;
//     bool havePionPlus;
//     for(unsigned int i = 0; i < Gen.NumDaughters(); ++i)
//       {
// 	const GenParticle daughter = Gen.GetDaughter(i);
// 	if(TMath::Abs(daughter.PDGId()) != 15 && daughter.PDGId() == 211)
// 	  {
// 	    assert(!havePionPlus);
// 	    Gen = daughter;
// 	    havePionPlus = true;
// 	  }
//       } // for(unsigned int i = 0; i < Gen.NumDaughters(); ++i)
//     return Gen;
//   }
  
  
//   GenParticle GetNPionFromTauDecay(const GenParticle& tau)
//   {
//     assert(TMath::Abs(tau.PDGId()) == 15);
    
//     // Skip internal nodes
//     GenParticle gen = tau;
//     bool havePionNeg;

//     for(unsigned int i = 0; i < gen.NumDaughters(); ++i)
//       {
// 	const GenParticle pion = gen.GetDaughter(i);
// 	if(pion.PDGId() == -211)
// 	  {
// 	    assert(!havePionNeg);
// 	    gen = pion;
// 	    havePionNeg = true;
// 	  }
//       } // for(unsigned int i = 0; i < Gen.NumDaughters(); ++i)
//     assert(havePionNeg);
//     return gen;
//   }
  
  

}// namespace gen
