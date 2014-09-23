#include "AnalysisTool/DataForm.h"
#include "Common/diTauCand.h"

diTauCand::diTauCand(const diTau& tau):
  diTau_Leg1Index(tau.diTau_Leg1Index()),
  diTau_Leg2Index(tau.diTau_Leg2Index()),
  diTau_Vtx(tau.diTau_Vtx()),
  diTau_VtxErr(tau.diTau_VtxErr()),
  diTau_reFitVtx(tau.diTau_reFitVtx()),
  diTau_reFitVtxErr(tau.diTau_reFitVtxErr()),
  diTau_Leg1_pos(tau.diTau_Leg1_pos()),
  diTau_Leg1_pos_OPV(tau.diTau_Leg1_pos_OPV()),
  diTau_Leg1_dxy(tau.diTau_Leg1_dxy()),
  diTau_Leg1_dz(tau.diTau_Leg1_dz()),
  diTau_Leg1_dxyErr(tau.diTau_Leg1_dxyErr()),
  diTau_Leg1_dzErr(tau.diTau_Leg1_dzErr()),
  diTau_Leg2_pos(tau.diTau_Leg2_pos()),
  diTau_Leg2_pos_OPV(tau.diTau_Leg2_pos_OPV()),
  diTau_Leg2_dxy(tau.diTau_Leg2_dxy()),
  diTau_Leg2_dz(tau.diTau_Leg2_dz()),
  diTau_Leg2_dxyErr(tau.diTau_Leg2_dxyErr()),
  diTau_Leg2_dzErr(tau.diTau_Leg2_dzErr())
{
}
