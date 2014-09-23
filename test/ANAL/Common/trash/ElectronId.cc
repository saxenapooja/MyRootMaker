#include "Common/ElectronId.h"

namespace ElectronId
{

bool tight(const Electron& el)
{
  // Tight working point from here: https://twiki.cern.ch/twiki/bin/viewauth/CMS/HiggsToTauTauWorking2012#Object_ID_AN1
  return (TMath::Abs(el.Eta()) < 0.8 && el.MVAIdNonTrig() > 0.925) ||
         (TMath::Abs(el.Eta()) > 0.8 && TMath::Abs(el.Eta()) < 1.479 && el.MVAIdNonTrig() > 0.975) ||
         (TMath::Abs(el.Eta()) > 1.479 && el.MVAIdNonTrig() > 0.985);
}

bool loose(const Electron& el)
{
  // Loose working point from here: https://twiki.cern.ch/twiki/bin/viewauth/CMS/HiggsToTauTauWorking2012#Object_ID_AN1
  if(el.Pt() < 20)
  {
    return (TMath::Abs(el.Eta()) < 0.8 && el.MVAIdNonTrig() > 0.925) ||
           (TMath::Abs(el.Eta()) > 0.8 && TMath::Abs(el.Eta()) < 1.479 && el.MVAIdNonTrig() > 0.915) ||
           (TMath::Abs(el.Eta()) > 1.479 && el.MVAIdNonTrig() > 0.965);
  }
  else
  {
    return (TMath::Abs(el.Eta()) < 0.8 && el.MVAIdNonTrig() > 0.905) ||
           (TMath::Abs(el.Eta()) > 0.8 && TMath::Abs(el.Eta()) < 1.479 && el.MVAIdNonTrig() > 0.955) ||
           (TMath::Abs(el.Eta()) > 1.479 && el.MVAIdNonTrig() > 0.975);
  }
}

} // namespace ElectronId
