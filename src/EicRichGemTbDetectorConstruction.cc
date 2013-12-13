#include "EicRichGemTbDetectorConstruction.hh"

#include "G4ProductionCuts.hh"
//#include "G4ElementTable.hh"
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
  expHall_x = 1.0*m;
  expHall_y = 1.0*m;
  expHall_z = 3.0*m;

  tank_r = 40.0*cm;
  tank_dr = 10.0*cm;
  tank_dz = 100.0*cm;

  // mirror relation focal length f vs radius curvature r: mirror_f = mirror_r / 2
  // mirror_cr = mirror cylinder radius
  // mirror_dz = mirror thickness along z axis
  mirror_f = 1.0*m;
  mirror_r = mirror_f * 2.0;
  mirror_cr = 40.0*cm;
  mirror_dz = 10*cm;
}


EicRichGemTbDetectorConstruction::~EicRichGemTbDetectorConstruction()
{

}


G4VPhysicalVolume* EicRichGemTbDetectorConstruction::Construct()
{
  richTbMaterial = new EicRichGemTbMaterial();
  richTbGeometry = new EicRichGemTbGeometry();

  // The experimental Hall
  //
  G4Box* expHall_box = new G4Box("World",0.5*expHall_x,0.5*expHall_y,0.5*expHall_z);

  G4LogicalVolume* expHall_log
    = new G4LogicalVolume(expHall_box,getRichTbMaterial()->getAir(),"World",0,0,0);

  G4VPhysicalVolume* expHall_phys
    = new G4PVPlacement(0,G4ThreeVector(),expHall_log,"World",0,false,0);

  // The RICH Gas Tank
  //
  G4Tubs* gasTank_cyl = new G4Tubs("Tank",tank_r,tank_r+tank_dr,0.5*tank_dz,0,2 * pi);

  G4LogicalVolume* gasTank_log
    = new G4LogicalVolume(gasTank_cyl,getRichTbMaterial()->getAluminum(),"Tank",0,0,0);

  G4VPhysicalVolume* gasTank_phys
    = new G4PVPlacement(0,G4ThreeVector(),gasTank_log,"Tank",
                        expHall_log,false,0);

  // The RICH Gas Volume
  //
  G4Tubs* gasVol_cyl = new G4Tubs("Gas",0,tank_r,0.5*tank_dz,0,2 * pi);

  G4LogicalVolume* gasVol_log
    = new G4LogicalVolume(gasVol_cyl,getRichTbMaterial()->getCF4(),"Gas",0,0,0);

  G4VPhysicalVolume* gasVol_phys
    = new G4PVPlacement(0,G4ThreeVector(),gasVol_log,"Gas",
                        expHall_log,false,0);

  // The RICH Mirror
  //
  // LHCb website: 3mm thick Be base + 0.3mm glass surface layer coated with Al.
  // LHCb TDR: The mirrors are made of polished 6mm-thick glass coated by vacuum
  //           deposition with 900 nm of aluminium and overcoated with 200 nm of quartz.

  //  G4Sphere* richMirror = new G4Sphere("RICHMirror", //
  //                                      mirror_r,
  //                                      mirror_r + mirror_dr, //            G4double pRmin, G4double pRmax,
  //                                      0, 2 * pi, //            G4double pSPhi, G4double pDPhi,
  //                                      0, pi //            G4double pSTheta, G4double pDTheta
  //                                      );

  G4Sphere* richMirrorSphere = new G4Sphere("RICHMirrorSphere", //
                                            0,
                                            mirror_r, //            G4double pRmin, G4double pRmax,
                                            0, 2 * pi, //            G4double pSPhi, G4double pDPhi,
                                            0, pi //            G4double pSTheta, G4double pDTheta
                                            );

  G4double mirror_d = sqrt(mirror_r*mirror_r - mirror_cr*mirror_cr); // distance center sphere to beginning cylinder
  G4double mirror_cdz = ( mirror_r - mirror_d ) + mirror_dz; // full mirror cylinder length

  G4Tubs* richMirrorCyl = new G4Tubs("RICHMirrorCyl", 0, mirror_cr, 0.5 * mirror_cdz, 0, 2 * pi);

  G4VSolid* richMirrorSphere_place = new G4DisplacedSolid( "RICHMirrorSphere_place", richMirrorSphere, 0,
                                                           G4ThreeVector(0, 0, -(mirror_d+0.5*mirror_cdz)) );

  //G4VSolid *richMirror = new G4SubtractionSolid("RICHMirror", richMirrorCyl, richMirrorSphere_place);
  //G4VSolid *richMirror = new G4UnionSolid("RICHMirror", richMirrorCyl, richMirrorSphere_place);
  //G4VSolid *richMirror = new G4UnionSolid("RICHMirror", richMirrorCyl, richMirrorCyl);

//  G4LogicalVolume *richMirror_log = new G4LogicalVolume(richMirror,
//                                                        getRichTbMaterial()->getMirrorQuartz(),
//                                                        "RICHMirror", 0, 0, 0);

  G4LogicalVolume *richMirror_log = new G4LogicalVolume(richMirrorCyl,
                                                        getRichTbMaterial()->getMirrorQuartz(),
                                                        "RICHMirror", 0, 0, 0);

  G4VPhysicalVolume* richMirror_phys
    = new G4PVPlacement(0,G4ThreeVector(0, 0, (0.5*tank_dz - 0.5*mirror_cdz-0.3*m)),richMirror_log,"RICHMirror",
                        expHall_log,false,0);


  // The RICH Photocathode
  // ...

  // The RICH Readout and GEM Stack
  // ...



  //##########################################################

  // ------------- Surface --------------
  // Mirror
  //

  // Now for the material properties of Surfaces
  //
  //Front (reflecting surface of RichTb Mirror)

  // First define wavelength in nm.
  // For now assume that all segments have the same reflectivity.
  // Hence the reflectivity is defined outside the loop of the
  // the number of segments.
  // Only the front surface is created.
  // The abosorption length is set to a small value just to
  // avoid photons exiting from the back of the mirror.
  // the efficiency is for the absorption process.

  G4int NumPhotonRichMirrorReflWaveLengthBins = 63;
  G4double PhotMomWaveConv=1243.125;

  //Mirror reflectivity
  // In the following, the bins at 100 nm, and 1000nm are
  // defined just for convenience of interpolation.
  // They are not measured points.

  static const G4double PhotonWavelengthRefl[]=
    {100.0, 200.0, 210.0, 220.0, 230.0, 240.0, 250.0, 260.0, 270.0, 280.0,
     290.0, 300.0, 310.0, 320.0, 330.0, 340.0, 350.0, 360.0, 370.0, 380.0,
     390.0, 400.0, 410.0, 420.0, 430.0, 440.0, 450.0, 460.0, 470.0, 480.0,
     490.0, 500.0, 510.0, 520.0, 530.0, 540.0, 550.0, 560.0, 570.0, 580.0,
     590.0, 600.0, 610.0, 620.0, 630.0, 640.0, 650.0, 660.0, 670.0, 680.0,
     690.0, 700.0, 710.0, 720.0, 730.0, 740.0, 750.0, 760.0, 770.0, 780.0,
     790.0, 800.0, 1000.0 };

  static const G4double RichTbMirrorReflectivity[]=
    {0.0, 0.9106, 0.9232, 0.9285, 0.9314, 0.9323, 0.9312, 0.9287, 0.9264,
     0.9234, 0.9195, 0.9156, 0.9109, 0.9066, 0.9022, 0.8981, 0.8925, 0.8883,
     0.8836, 0.8796, 0.8756, 0.8727, 0.8697, 0.8672, 0.8653, 0.8636, 0.8624,
     0.8612, 0.8608, 0.8601, 0.8601, 0.8601, 0.8600, 0.8603, 0.8603, 0.8604,
     0.8605, 0.8608, 0.8609, 0.8608, 0.8608, 0.8606, 0.8604, 0.8600, 0.8598,
     0.8591, 0.8581, 0.8573, 0.8563, 0.8549, 0.8535, 0.8517, 0.8497, 0.8475,
     0.8447, 0.8417, 0.8382, 0.8388, 0.8296, 0.8258, 0.8204, 0.8172, 0.8172 };

  static const G4double RichTbMirrorEfficiency[]=
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0 };

  G4double* PhotonMomentumRefl
    =new G4double[NumPhotonRichMirrorReflWaveLengthBins];
  G4double* PhotWaveRefl =
    new  G4double[NumPhotonRichMirrorReflWaveLengthBins];
  G4double* PhotReflEff =new  G4double[NumPhotonRichMirrorReflWaveLengthBins];
  G4double* MirrorQuRefIndex
    =new  G4double[NumPhotonRichMirrorReflWaveLengthBins];

  for (G4int ibin=0; ibin<NumPhotonRichMirrorReflWaveLengthBins; ibin++){
    PhotonMomentumRefl[ibin]=PhotMomWaveConv*eV/ PhotonWavelengthRefl[ibin];
    PhotWaveRefl[ibin]=  RichTbMirrorReflectivity[ibin];
    PhotReflEff[ibin]= RichTbMirrorEfficiency[ibin];
    //the following lines to avoid reflection at the mirror.

    MirrorQuRefIndex[ibin] = 1.40;
  }

  G4OpticalSurface * OpRichTbMirrorSurface = new G4OpticalSurface("RichTbMirrorSurface");

//  G4LogicalBorderSurface* MirrorSurface = new G4LogicalBorderSurface("RichTbMirrorSurface",
//                                                                     gasVol_phys, richMirror_phys,
//                                                                     OpRichTbMirrorSurface);

  G4LogicalSkinSurface* MirrorSurface = new G4LogicalSkinSurface("RichTbMirrorSurface",
                                                                     richMirror_log,
                                                                     OpRichTbMirrorSurface);

//  G4LogicalSkinSurface* MirrorSurface = new G4LogicalSkinSurface("RichTbMirrorSurface",
//                                                                     gasVol_log,
//                                                                     OpRichTbMirrorSurface);

  OpRichTbMirrorSurface->SetType(dielectric_metal);
  OpRichTbMirrorSurface->SetModel(glisur);
  //OpRichTbMirrorSurface->SetModel(unified);
  OpRichTbMirrorSurface->SetFinish(polished);

//  // test from G4 example
//  const G4int NUM = 2;
//
//  G4double pp[NUM] = {2.038*eV, 4.144*eV};
//  G4double specularlobe[NUM] = {0.3, 0.3};
//  G4double specularspike[NUM] = {0.2, 0.2};
//  G4double backscatter[NUM] = {0.1, 0.1};
//  G4double rindex[NUM] = {1.35, 1.40};
//  G4double reflectivity[NUM] = {1.0, 1.0};
//  G4double efficiency[NUM] = {0.0, 0.0};
//
//  G4MaterialPropertiesTable* SMPT = new G4MaterialPropertiesTable();
//
//  SMPT -> AddProperty("RINDEX",pp,rindex,NUM);
//  SMPT -> AddProperty("SPECULARLOBECONSTANT",pp,specularlobe,NUM);
//  SMPT -> AddProperty("SPECULARSPIKECONSTANT",pp,specularspike,NUM);
//  SMPT -> AddProperty("BACKSCATTERCONSTANT",pp,backscatter,NUM);
//  SMPT -> AddProperty("REFLECTIVITY",pp,reflectivity,NUM);
//  SMPT -> AddProperty("EFFICIENCY",pp,efficiency,NUM);

//  OpRichTbMirrorSurface->SetMaterialPropertiesTable(SMPT);
  // end test

  G4MaterialPropertiesTable* OpRichTbMirrorSurfaceMPT = new G4MaterialPropertiesTable();

  OpRichTbMirrorSurfaceMPT->AddProperty("RINDEX",
                                        PhotonMomentumRefl,
                                        MirrorQuRefIndex,
                                        NumPhotonRichMirrorReflWaveLengthBins);
  OpRichTbMirrorSurfaceMPT->AddProperty("REFLECTIVITY",
                                        PhotonMomentumRefl,
                                        PhotWaveRefl,
                                        NumPhotonRichMirrorReflWaveLengthBins);
  OpRichTbMirrorSurfaceMPT->AddProperty("EFFICIENCY",
                                        PhotonMomentumRefl,
                                        PhotReflEff,
                                        NumPhotonRichMirrorReflWaveLengthBins);

  OpRichTbMirrorSurface->SetMaterialPropertiesTable(OpRichTbMirrorSurfaceMPT);

  OpRichTbMirrorSurface->DumpInfo();

  //##########################################################

  // always return the physical World
  return expHall_phys;

}
