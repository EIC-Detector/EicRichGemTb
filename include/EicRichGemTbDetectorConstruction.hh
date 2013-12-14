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

};

#endif /* EicRichGemTbDETECTORCONSTRUCTION_H_ */
