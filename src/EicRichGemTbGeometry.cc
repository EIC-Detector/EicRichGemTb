#include "EicRichGemTbGeometry.hh"

EicRichGemTbGeometry::EicRichGemTbGeometry(){
  SetDefault();
  CalcDependentParameters();
}

EicRichGemTbGeometry::~EicRichGemTbGeometry() {}

void EicRichGemTbGeometry::SetDefault(){

  world_x = 1.0 * m;
  world_y = 1.0 * m;
  world_z = 3.0 * m;

  pressVess_r_inner = 40.0 * cm;
  pressVess_dr = 10.0 * cm;
  pressVess_dz = 100.0 * cm;

  mirror_f = 1.0 * m;

  mirror_sphere_r = mirror_f * 2.0;

  mirror_cyl_r = pressVess_r_inner;
  mirror_cyl_l = 10.0 * cm;

  mirror_pos_z = (0.5 * pressVess_dz - 0.5 * mirror_cyl_l);

  // simple assumption: curvature extends half way along cylinder, i.e. minimum mirror thickness = 0.5 * maximum mirror thickness
  mirror_cyl_sphere_dist_z = -( mirror_sphere_r );


  csilayer_x = 10.0 * cm;
  csilayer_y = 10.0 * cm;
  csilayer_z = 0.001 * mm;
  csilayer_pos_z = -0.5 * GetPressureVesselLength() + 2.0 * cm;

  return;

}

void EicRichGemTbGeometry::CalcDependentParameters() {

  return;

}
