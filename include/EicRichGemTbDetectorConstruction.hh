/*!
 * \file ${file_name}
 * \brief
 * \author Nils Feege <nils.feege@stonybrook.edu>
 * \version $$Revision: 1.0 $$
 * \date $$Date: 2013/11/15 18:11:08 $$
 */

#ifndef EicRichGemTbDETECTORCONSTRUCTION_H_
#define EicRichGemTbDETECTORCONSTRUCTION_H_ 1

#include <cmath>
#include <map>
#include <utility>

#include "G4PhysicalConstants.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4String.hh"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"

#include "EicRichGemTbMaterial.hh"
#include "EicRichGemTbGeometry.hh"

class G4VPhysicalVolume;
class G4LogicalVolume;
class G4PVPlacement;

class EicRichGemTbDetectorConstruction : public G4VUserDetectorConstruction
{
public:
  EicRichGemTbDetectorConstruction();
  ~EicRichGemTbDetectorConstruction();

  G4VPhysicalVolume* Construct();

  EicRichGemTbMaterial* getRichTbMaterial() {return  richTbMaterial; }
  EicRichGemTbGeometry* getRichTbGeometry() {return  richTbGeometry; }

protected:

private:

  G4Box* fExperimentalHall_box;
  G4LogicalVolume* fExperimentalHall_log;
  G4VPhysicalVolume* fExperimentalHall_phys;

  EicRichGemTbMaterial* richTbMaterial;
  EicRichGemTbGeometry* richTbGeometry;

  // Geometry
  G4double expHall_x;
  G4double expHall_y;
  G4double expHall_z;

  G4double tank_r;
  G4double tank_dr;
  G4double tank_dz;

  G4double mirror_f;
  G4double mirror_r;
  G4double mirror_cr;
  G4double mirror_dz;

};

#endif /* EicRichGemTbDETECTORCONSTRUCTION_H_ */
