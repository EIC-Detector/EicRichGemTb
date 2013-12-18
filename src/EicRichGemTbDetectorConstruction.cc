#include "EicRichGemTbDetectorConstruction.hh"
#include "EicRichGemTbSD.hh"

#include "globals.hh"
#include "G4ProductionCuts.hh"
#include "G4Tubs.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4SDManager.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4OpticalSurface.hh"

#include "G4NistManager.hh"

#include "G4Ellipsoid.hh"
#include "G4Sphere.hh"
#include "G4Orb.hh"
#include "G4Cons.hh"
#include "G4DisplacedSolid.hh"
#include "G4ExtrudedSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"

#include <cassert>
#include <iostream>
#include <sstream>
using namespace std;

EicRichGemTbDetectorConstruction::EicRichGemTbDetectorConstruction()
{

}


EicRichGemTbDetectorConstruction::~EicRichGemTbDetectorConstruction()
{

}


G4VPhysicalVolume* EicRichGemTbDetectorConstruction::Construct()
{
  richTbMaterial = new EicRichGemTbMaterial();
  richTbGeometry = new EicRichGemTbGeometry();

  // The World
  //
  G4Box* world_box = new G4Box("World",
                               0.5 * getRichTbGeometry()->GetWorldX(),
                               0.5 * getRichTbGeometry()->GetWorldY(),
                               0.5 * getRichTbGeometry()->GetWorldZ()
                               );

  G4LogicalVolume* world_log = new G4LogicalVolume(world_box,
                                                   getRichTbMaterial()->getAir(),
                                                   "World",
                                                   0,
                                                   0,
                                                   0
                                                   );

  G4VPhysicalVolume* world_phys = new G4PVPlacement(0,
                                                    G4ThreeVector(),
                                                    world_log,
                                                    "World",
                                                    0,
                                                    false,
                                                    0
                                                    );

  // The RICH Pressure Vessel
  //
  G4Tubs* pressVess_cyl = new G4Tubs("PressureVessel",
                                     getRichTbGeometry()->GetPressureVesselInnerRadius(),
                                     getRichTbGeometry()->GetPressureVesselInnerRadius()
                                     + getRichTbGeometry()->GetPressureVesselThickness(),
                                     0.5 * getRichTbGeometry()->GetPressureVesselLength(),
                                     0,
                                     2 * pi
                                     );

  G4LogicalVolume* pressVess_log = new G4LogicalVolume(pressVess_cyl,
                                                       getRichTbMaterial()->getStainlessSteel(),
                                                       "PressureVessel",
                                                       0, 0, 0 );

  G4VPhysicalVolume* pressVess_phys = new G4PVPlacement(0,
                                                        G4ThreeVector(),
                                                        pressVess_log,
                                                        "PressureVessel",
                                                        world_log,
                                                        false,
                                                        0);
  pressVess_phys->GetName(); // avoid warning for unused variable during compilation

  // The RICH Gas Volume
  //
  G4Tubs* gasVol_cyl = new G4Tubs( "GasVolume",
                                   0,
                                   getRichTbGeometry()->GetPressureVesselInnerRadius(),
                                   0.5 * getRichTbGeometry()->GetPressureVesselLength(),
                                   0,
                                   2 * pi);

  G4LogicalVolume* gasVol_log = new G4LogicalVolume( gasVol_cyl,
                                                     getRichTbMaterial()->getCF4(),
                                                     "GasVolume", 0, 0, 0 );

  G4VPhysicalVolume* gasVol_phys = new G4PVPlacement(0,
                                                     G4ThreeVector(),
                                                     gasVol_log,
                                                     "GasVolume",
                                                     world_log,false,0);

  // The RICH Mirror
  //
  G4Tubs* richMirrorCylinder = new G4Tubs("RICHMirrorCylinder"
                                          ,
                                          0,
                                          getRichTbGeometry()->GetMirrorCylinderRadius(),
                                          0.5 * getRichTbGeometry()->GetMirrorCylinderLength(),
                                          0, 2 * pi);

  G4Sphere* richMirrorSphere = new G4Sphere("RICHMirrorSphere", //
                                            0, getRichTbGeometry()->GetMirrorSphereRadius(), //            G4double pRmin, G4double pRmax,
                                            0, 2 * pi, //            G4double pSPhi, G4double pDPhi,
                                            0, pi //            G4double pSTheta, G4double pDTheta
                                            );

  G4VSolid* richMirrorSphere_shift = new G4DisplacedSolid( "RICHMirrorSphere_shift",
                                                           richMirrorSphere,
                                                           0,
                                                           G4ThreeVector(0, 0, getRichTbGeometry()->GetMirrorCylinderSphereDistanceZ()) );
  richMirrorSphere_shift->GetName(); // avoid warning for unused variable during compilation

  G4VSolid *richMirror = new G4SubtractionSolid("RICHMirror", richMirrorCylinder, richMirrorSphere_shift);

  G4LogicalVolume *richMirror_log = new G4LogicalVolume(richMirror,
                                                        //getRichTbMaterial()->getMirrorQuartz(),
                                                        getRichTbMaterial()->getAluminum(),
                                                        "RICHMirror",
                                                        0, 0, 0);

  G4VPhysicalVolume* richMirror_phys = new G4PVPlacement(0,
                                                         G4ThreeVector(0, 0, getRichTbGeometry()->GetMirrorPositionZ()),
                                                         richMirror_log,
                                                         "RICHMirror",
                                                         gasVol_log,false, false);


  // Mirror surface
  //
  G4LogicalBorderSurface* MirrorSurface = new G4LogicalBorderSurface("RichTbMirrorSurface",
                                                                     gasVol_phys, richMirror_phys,
                                                                     getRichTbMaterial()->getOpticalMirrorSurface() );

  //  G4LogicalBorderSurface* MirrorSurface = new G4LogicalBorderSurface("RichTbMirrorSurface",
  //                                                                     gasVol_phys, world_phys,
  //                                                                     getRichTbMaterial()->getOpticalMirrorSurface() );

  //G4LogicalSkinSurface* MirrorSurface = new G4LogicalSkinSurface("RichTbMirrorSurface",
  //                                                             richMirror_log,
  //                                                             getRichTbMaterial()->getOpticalMirrorSurface() );

  MirrorSurface->GetName(); // avoid warning for unused variable during compilation


  // The RICH Photon Detection: CsI Cathode, GEM Stack, segemnted readout plane
  // For now: implemented as vacuum volume with absorbing surface
  //
  G4Box* GEMStack_box = new G4Box("GEMStack",
                                  0.5 * getRichTbGeometry()->GetGEMStackX(),
                                  0.5 * getRichTbGeometry()->GetGEMStackY(),
                                  0.5 * getRichTbGeometry()->GetGEMStackZ()
                                  );

  G4LogicalVolume* GEMStack_log = new G4LogicalVolume(GEMStack_box,
                                                      getRichTbMaterial()->getVacuum(),
                                                      "GEMStack",
                                                      0,
                                                      0,
                                                      0
                                                      );

  G4VPhysicalVolume* GEMStack_phys = new G4PVPlacement(0,
                                                       G4ThreeVector(0, 0, getRichTbGeometry()->GetGEMStackPositionZ()),
                                                       GEMStack_log,
                                                       "GEMStack",
                                                       gasVol_log,
                                                       false,
                                                       0 );

  // Give GEM stack an optical surface
  //
  G4LogicalSkinSurface* PhotocathodeSurface = new G4LogicalSkinSurface("PhotocathodeSurface",GEMStack_log,getRichTbMaterial()->getOpticalPhotocathodeSurface());

  // Make photocathode sensitive Detector
  //
  G4SDManager* SDman = G4SDManager::GetSDMpointer();
  G4String photoSDname = "/EicRichGemTbDet/photoSD";
  G4String photoHitColname="PhotoHitsCollection";

  EicRichGemTbSD* photoSD = new EicRichGemTbSD( photoSDname, photoHitColname );

  SDman->AddNewDetector( photoSD );
  GEMStack_log->SetSensitiveDetector( photoSD );


  // always return the physical World
  //
  return world_phys;

}
