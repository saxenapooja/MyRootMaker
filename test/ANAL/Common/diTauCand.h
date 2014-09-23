#ifndef _HTAUTAUANALYSIS_DITAUCAND_H_
#define _HTAUTAUANALYSIS_DITAUCAND_H_

#include "AnalysisTool/DataForm.h"
#include "Common/Config.h"

class diTauCand
{
public:
  diTauCand(const diTau& tau);

  //ditau vertex
  float ditau_Vtx_X() const  { return diTau_Vtx.X(); }
  float ditau_Vtx_Y() const  { return diTau_Vtx.Y(); }
  float ditau_Vtx_Z() const  { return diTau_Vtx.Z(); }

  //ditau vertexErr
  float ditau_Vtx_XErr() const  { return diTau_VtxErr.X(); }
  float ditau_Vtx_YErr() const  { return diTau_VtxErr.Y(); }
  float ditau_Vtx_ZErr() const  { return diTau_VtxErr.Z(); }

  //ditau Leg1 position
  float diTau_Leg1_X()   const { return diTau_Leg1_pos.X(); }
  float diTau_Leg1_Y()   const { return diTau_Leg1_pos.Y(); }
  float diTau_Leg1_Z()   const { return diTau_Leg1_pos.Z(); }

  //ditau Leg1 position wrt Original Primary Vertex
  float diTau_Leg1_X_OPV()   const { return diTau_Leg1_pos_OPV.X(); }
  float diTau_Leg1_Y_OPV()   const { return diTau_Leg1_pos_OPV.Y(); }
  float diTau_Leg1_Z_OPV()   const { return diTau_Leg1_pos_OPV.Z(); }

  //ditau Leg2 position
  float diTau_Leg2_X()   const { return diTau_Leg2_pos.X(); }
  float diTau_Leg2_Y()   const { return diTau_Leg2_pos.Y(); }
  float diTau_Leg2_Z()   const { return diTau_Leg2_pos.Z(); }

  //ditau Leg2 position wrt Original Primary Vertex
  float diTau_Leg2_X_OPV()   const { return diTau_Leg2_pos_OPV.X(); }
  float diTau_Leg2_Y_OPV()   const { return diTau_Leg2_pos_OPV.Y(); }
  float diTau_Leg2_Z_OPV()   const { return diTau_Leg2_pos_OPV.Z(); }

  int diTau_Leg1Index;
  int diTau_Leg2Index;
  TVector3 diTau_Vtx;
  TVector3 diTau_VtxErr;
  TVector3 diTau_reFitVtx;
  TVector3 diTau_reFitVtxErr;
  TVector3 diTau_Leg1_pos;
  TVector3 diTau_Leg1_pos_OPV;
  TVector3 diTau_Leg2_pos;
  TVector3 diTau_Leg2_pos_OPV;
  float diTau_Leg1_dxy;
  float diTau_Leg1_dz;
  float diTau_Leg1_dxyErr;
  float diTau_Leg1_dzErr;
  float diTau_Leg2_dxy;
  float diTau_Leg2_dz;
  float diTau_Leg2_dxyErr;
  float diTau_Leg2_dzErr;
};

#endif 
