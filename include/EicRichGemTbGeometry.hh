#ifndef EicRichGemTbDETECTORGEOMETRY_H_
#define EicRichGemTbDETECTORGEOMETRY_H_ 1

#include "globals.hh"

class EicRichGemTbGeometry
{
public:
  EicRichGemTbGeometry();
  ~EicRichGemTbGeometry();

  void SetDefault();
  void CalcDependentParameters();

  G4double GetWorldX() { return world_x; }
  G4double GetWorldY() { return world_y; }
  G4double GetWorldZ() { return world_z; }

  G4double GetPressureVesselInnerRadius() { return pressVess_r_inner; }
  G4double GetPressureVesselThickness() { return pressVess_dr; }
  G4double GetPressureVesselLength() { return pressVess_dz; }

  G4double GetMirrorCylinderRadius() {return mirror_cyl_r; }
  G4double GetMirrorCylinderLength() {return mirror_cyl_l; }
  G4double GetMirrorSphereRadius() {return mirror_sphere_r; }
  G4double GetMirrorPositionZ() {return mirror_pos_z; }
  G4double GetMirrorCylinderSphereDistanceZ() {return mirror_cyl_sphere_dist_z; }

private:

  G4double world_x;
  G4double world_y;
  G4double world_z;

  G4double pressVess_r_inner;
  G4double pressVess_dr;
  G4double pressVess_dz;

  G4double mirror_f;
  G4double mirror_pos_z;
  G4double mirror_cyl_r;
  G4double mirror_cyl_l;
  G4double mirror_sphere_r;
  G4double mirror_cyl_sphere_dist_z;

};

#endif
