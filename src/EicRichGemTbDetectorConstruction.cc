// $$Id: EicRichGemTbDetectorConstruction.cc,v 1.3 2013/10/13 03:47:10 jinhuang Exp $$

/*!
 * \file ${file_name}
 * \brief
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $$Revision: 1.3 $$
 * \date $$Date: 2013/10/13 03:47:10 $$
 */

#include "EicRichGemTbDetectorConstruction.hh"

#include "G4ProductionCuts.hh"
#include "G4ElementTable.hh"
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


void EicRichGemTbDetectorConstruction::DefineMaterials()
{
  G4double a; //atomic mass
  G4double z;  // atomic number
  G4double density;
  G4double temperature;
  G4double pressure;

  G4int nelements;
  G4int natoms;

  // Elements
  //
  fH  = new G4Element("Hydrogen", "H", z=1 , a=1.01*g/mole);
  fN  = new G4Element("Nitrogen", "N", z=7 , a=14.01*g/mole);
  fO  = new G4Element("Oxygen"  , "O", z=8 , a=16.00*g/mole);
  fC  = new G4Element("Carbon", "C", z=6 , a=12.01*g/mole);
  fF  = new G4Element("Fluorine", "F", z=9 , a=19.00*g/mole);
  fSi = new G4Element("Silicon", "Si", z=14, a=28.09*g/mole);

  // Air
  //
  fAir = new G4Material("Air", density=1.29*mg/cm3, nelements=2);
  fAir->AddElement(fN, 70.*perCent);
  fAir->AddElement(fO, 30.*perCent);

  // CF4
  //

  // LHCb; no data available at room temp and pressure;
  fCF4 = new G4Material("CF4", density=0.003884*g/cm3, nelements=2, kStateGas, temperature=273.*kelvin, pressure=1.0*atmosphere);

  // http://www.slac.stanford.edu/pubs/icfa/summer98/paper3/paper3.pdf : density=0.00393*g/cm^3, T = 20degC

  // http://encyclopedia.airliquide.com/Encyclopedia.asp?GasID=61#GeneralData :
  // Molecular weight  : 88.01 g/mol
  // Gas density (1.013 bar and 15 °C (59 °F)) : 3.72 kg/m3

  // example
  //fCF4 = new G4Material("CF4", density=3.72*g/cm3, nelements=2);

  fCF4->AddElement(fC, 1);
  fCF4->AddElement(fF, 4);

  // Aluminum
  //
  fAl = new G4Material("Aluminum", z=13 , a=26.98*g/mole , density=2.7*g/cm3);

  G4cout << *(G4Material::GetMaterialTable()) << endl;

  //  LHCb
  //  //Aluminium
  //  density=2.7*g/cm3;
  //  G4Material* Aluminium =new G4Material(name="Aluminium",density,numel=1);
  //  Aluminium->AddElement(elAL,natoms=1);
  //
  //  G4double* AluminiumAbsorpLength=new G4double[NumPhotWaveLengthBins];
  //
  //  for (ibin=0; ibin<NumPhotWaveLengthBins; ibin++){
  //    AluminiumAbsorpLength[ibin]=0.0*mm;
  //  }
  //
  //  G4MaterialPropertiesTable* AluminiumMPT =
  //    new G4MaterialPropertiesTable();
  //
  //  AluminiumMPT->AddProperty("ABSLENGTH",PhotonMomentum,
  //                            AluminiumAbsorpLength,NumPhotWaveLengthBins);
  //
  //  Aluminium->SetMaterialPropertiesTable(AluminiumMPT);
  //  RichTbAluminium=Aluminium;

  // LHCb: There is a quartz for the mirror and
  // another quartz which is used in aerogel and
  // yet another quartz used for the quartz window.
  // Mirrorquartz

  fSiO2MirrorQuartz = new G4Material("MirrorQuartz", density=2.200*g/cm3, nelements=2);
  fSiO2MirrorQuartz->AddElement(fSi,natoms=1);
  fSiO2MirrorQuartz->AddElement(fO,natoms=2);

  //  G4double* MirrorQuartzRindex=new G4double[NumPhotWaveLengthBins];
  //  G4double* MirrorQuartzAbsorpLength=new G4double[NumPhotWaveLengthBins];
  //  for (ibin=0; ibin<NumPhotWaveLengthBins; ibin++){
  //    MirrorQuartzAbsorpLength[ibin]=0.01*mm;
  //
  //  }
  //  G4MaterialPropertiesTable* MirrorQuartzMPT =
  //    new G4MaterialPropertiesTable();
  //
  //
  //  MirrorQuartzMPT->AddProperty("ABSLENGTH",PhotonMomentum,
  //                               MirrorQuartzAbsorpLength,NumPhotWaveLengthBins);
  //
  //
  //  SiO2MirrorQuartz->SetMaterialPropertiesTable(MirrorQuartzMPT);


  // Generate & Add Material Properties Table
  //

  // CF4
  //

  // Need to check values. These values are from ExN06 for Water.
  //
  const G4int cf4_nEntries = 32;

  G4double cf4_PhotonEnergy[cf4_nEntries]    =
    { 2.034*eV, 2.068*eV, 2.103*eV, 2.139*eV,
      2.177*eV, 2.216*eV, 2.256*eV, 2.298*eV,
      2.341*eV, 2.386*eV, 2.433*eV, 2.481*eV,
      2.532*eV, 2.585*eV, 2.640*eV, 2.697*eV,
      2.757*eV, 2.820*eV, 2.885*eV, 2.954*eV,
      3.026*eV, 3.102*eV, 3.181*eV, 3.265*eV,
      3.353*eV, 3.446*eV, 3.545*eV, 3.649*eV,
      3.760*eV, 3.877*eV, 4.002*eV, 4.136*eV };

  //  G4double cf4_RefractiveIndex[cf4_nEntries]  =
  //    { 1.3435, 1.344,  1.3445, 1.345,  1.3455,
  //      1.346,  1.3465, 1.347,  1.3475, 1.348,
  //      1.3485, 1.3492, 1.35,   1.3505, 1.351,
  //      1.3518, 1.3522, 1.3530, 1.3535, 1.354,
  //      1.3545, 1.355,  1.3555, 1.356,  1.3568,
  //      1.3572, 1.358,  1.3585, 1.359,  1.3595,
  //      1.36,   1.3608};

  G4double cf4_Absorption[cf4_nEntries]  =
    {3.448*m,  4.082*m,  6.329*m,  9.174*m, 12.346*m, 13.889*m,
     15.152*m, 17.241*m, 18.868*m, 20.000*m, 26.316*m, 35.714*m,
     45.455*m, 47.619*m, 52.632*m, 52.632*m, 55.556*m, 52.632*m,
     52.632*m, 47.619*m, 45.455*m, 41.667*m, 37.037*m, 33.333*m,
     30.000*m, 28.500*m, 27.000*m, 24.500*m, 22.000*m, 19.500*m,
     17.500*m, 14.500*m };

  //  G4double cf4_ScintilFast[cf4_nEntries] =
  //    { 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
  //      1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
  //      1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
  //      1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
  //      1.00, 1.00, 1.00, 1.00 };
  //
  //  G4double cf4_ScintilSlow[cf4_nEntries] =
  //    { 0.01, 1.00, 2.00, 3.00, 4.00, 5.00, 6.00,
  //      7.00, 8.00, 9.00, 8.00, 7.00, 6.00, 4.00,
  //      3.00, 2.00, 1.00, 0.01, 1.00, 2.00, 3.00,
  //      4.00, 5.00, 6.00, 7.00, 8.00, 9.00, 8.00,
  //      7.00, 6.00, 5.00, 4.00 };

  // use const refractive index from Thom's presentation
  G4double cf4_RefractiveIndex[cf4_nEntries]  =
    { 1.00054, 1.00054, 1.00054, 1.00054, 1.00054,
      1.00054, 1.00054, 1.00054, 1.00054, 1.00054,
      1.00054, 1.00054, 1.00054, 1.00054, 1.00054,
      1.00054, 1.00054, 1.00054, 1.00054, 1.00054,
      1.00054, 1.00054, 1.00054, 1.00054, 1.00054,
      1.00054, 1.00054, 1.00054, 1.00054, 1.00054,
      1.00054, 1.00054};

  fMPT_cf4 = new G4MaterialPropertiesTable();

  fMPT_cf4->AddProperty("RINDEX",        cf4_PhotonEnergy, cf4_RefractiveIndex, cf4_nEntries )->SetSpline(true);
  fMPT_cf4->AddProperty("ABSLENGTH",     cf4_PhotonEnergy, cf4_Absorption,      cf4_nEntries )->SetSpline(true);
  //  fMPT_cf4->AddProperty("FASTCOMPONENT", cf4_PhotonEnergy, cf4_ScintilFast,     cf4_nEntries )->SetSpline(true);
  //  fMPT_cf4->AddProperty("SLOWCOMPONENT", cf4_PhotonEnergy, cf4_ScintilSlow,     cf4_nEntries )->SetSpline(true);

  fMPT_cf4->AddConstProperty("SCINTILLATIONYIELD",50./MeV);
  //  fMPT_cf4->AddConstProperty("RESOLUTIONSCALE",1.0);
  //  fMPT_cf4->AddConstProperty("FASTTIMECONSTANT", 1.*ns);
  //  fMPT_cf4->AddConstProperty("SLOWTIMECONSTANT",10.*ns);
  //  fMPT_cf4->AddConstProperty("YIELDRATIO",0.8);

  //   const G4int NUMENTRIES_cf4 = 60;
  //
  //   G4double ENERGY_cf4[NUMENTRIES_cf4] = {
  //     1.56962*eV, 1.58974*eV, 1.61039*eV, 1.63157*eV,
  //     1.65333*eV, 1.67567*eV, 1.69863*eV, 1.72222*eV,
  //     1.74647*eV, 1.77142*eV, 1.7971 *eV, 1.82352*eV,
  //     1.85074*eV, 1.87878*eV, 1.90769*eV, 1.93749*eV,
  //     1.96825*eV, 1.99999*eV, 2.03278*eV, 2.06666*eV,
  //     2.10169*eV, 2.13793*eV, 2.17543*eV, 2.21428*eV,
  //     2.25454*eV, 2.29629*eV, 2.33962*eV, 2.38461*eV,
  //     2.43137*eV, 2.47999*eV, 2.53061*eV, 2.58333*eV,
  //     2.63829*eV, 2.69565*eV, 2.75555*eV, 2.81817*eV,
  //     2.88371*eV, 2.95237*eV, 3.02438*eV, 3.09999*eV,
  //     3.17948*eV, 3.26315*eV, 3.35134*eV, 3.44444*eV,
  //     3.54285*eV, 3.64705*eV, 3.75757*eV, 3.87499*eV,
  //     3.99999*eV, 4.13332*eV, 4.27585*eV, 4.42856*eV,
  //     4.59258*eV, 4.76922*eV, 4.95999*eV, 5.16665*eV,
  //     5.39129*eV, 5.63635*eV, 5.90475*eV, 6.19998*eV
  //   };
  //
  //   //assume 100 times larger than the rayleigh scattering for now.
  //   G4double MIE_cf4[NUMENTRIES_cf4] = {
  //     167024.4*m, 158726.7*m, 150742  *m,
  //     143062.5*m, 135680.2*m, 128587.4*m,
  //     121776.3*m, 115239.5*m, 108969.5*m,
  //     102958.8*m, 97200.35*m, 91686.86*m,
  //     86411.33*m, 81366.79*m, 76546.42*m,
  //     71943.46*m, 67551.29*m, 63363.36*m,
  //     59373.25*m, 55574.61*m, 51961.24*m,
  //     48527.00*m, 45265.87*m, 42171.94*m,
  //     39239.39*m, 36462.50*m, 33835.68*m,
  //     31353.41*m, 29010.30*m, 26801.03*m,
  //     24720.42*m, 22763.36*m, 20924.88*m,
  //     19200.07*m, 17584.16*m, 16072.45*m,
  //     14660.38*m, 13343.46*m, 12117.33*m,
  //     10977.70*m, 9920.416*m, 8941.407*m,
  //     8036.711*m, 7202.470*m, 6434.927*m,
  //     5730.429*m, 5085.425*m, 4496.467*m,
  //     3960.210*m, 3473.413*m, 3032.937*m,
  //     2635.746*m, 2278.907*m, 1959.588*m,
  //     1675.064*m, 1422.710*m, 1200.004*m,
  //     1004.528*m, 833.9666*m, 686.1063*m
  //   };
  //
  //   // gforward, gbackward, forward backward ratio
  //   G4double MIE_cf4_const[3]={0.99,0.99,0.8};
  //
  //   fMPT_cf4->AddProperty("MIEHG",ENERGY_cf4,MIE_cf4,NUMENTRIES_cf4)
  //     ->SetSpline(true);
  //   fMPT_cf4->AddConstProperty("MIEHG_FORWARD",MIE_cf4_const[0]);
  //   fMPT_cf4->AddConstProperty("MIEHG_BACKWARD",MIE_cf4_const[1]);
  //   fMPT_cf4->AddConstProperty("MIEHG_FORWARD_RATIO",MIE_cf4_const[2]);

  fCF4->SetMaterialPropertiesTable(fMPT_cf4);

  // Set the Birks Constant for the Cf4 scintillator
  //fCF4->GetIonisation()->SetBirksConstant(0.126*mm/MeV);

  return;

}

G4VPhysicalVolume* EicRichGemTbDetectorConstruction::ConstructDetector()
{

  // The experimental Hall
  //
  G4Box* expHall_box = new G4Box("World",0.5*expHall_x,0.5*expHall_y,0.5*expHall_z);

  G4LogicalVolume* expHall_log
    = new G4LogicalVolume(expHall_box,fAir,"World",0,0,0);

  G4VPhysicalVolume* expHall_phys
    = new G4PVPlacement(0,G4ThreeVector(),expHall_log,"World",0,false,0);

  // The RICH Gas Tank
  //
  G4Tubs* gasTank_cyl = new G4Tubs("Tank",tank_r,tank_r+tank_dr,0.5*tank_dz,0,2 * pi);

  G4LogicalVolume* gasTank_log
    = new G4LogicalVolume(gasTank_cyl,fAl,"Tank",0,0,0);

  G4VPhysicalVolume* gasTank_phys
    = new G4PVPlacement(0,G4ThreeVector(),gasTank_log,"Tank",
                        expHall_log,false,0);

  // The RICH Gas Volume
  //
  G4Tubs* gasVol_cyl = new G4Tubs("Gas",0,tank_r,0.5*tank_dz,0,2 * pi);

  G4LogicalVolume* gasVol_log
    = new G4LogicalVolume(gasVol_cyl,fCF4,"Gas",0,0,0);

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
//                                                        fSiO2MirrorQuartz,
//                                                        "RICHMirror", 0, 0, 0);

  G4LogicalVolume *richMirror_log = new G4LogicalVolume(richMirrorCyl,
                                                        fSiO2MirrorQuartz,
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

G4VPhysicalVolume* EicRichGemTbDetectorConstruction::Construct()
{
  DefineMaterials();
  return ConstructDetector();
}

G4VPhysicalVolume* EicRichGemTbDetectorConstruction::Construct2()
{

  G4double a, z, density;
  G4int nelements;

  // Air
  //
  G4Element* N = new G4Element("Nitrogen", "N", z=7 , a=14.01*g/mole);
  G4Element* O = new G4Element("Oxygen"  , "O", z=8 , a=16.00*g/mole);

  G4Material* Air = new G4Material("Air", density=1.29*mg/cm3, nelements=2);
  Air->AddElement(N, 70.*perCent);
  Air->AddElement(O, 30.*perCent);

  //
  // ------------- Volumes --------------

  // The experimental Hall
  //
  G4Box* expHall_box = new G4Box("World",0.5*expHall_x,0.5*expHall_y,0.5*expHall_z);

  G4LogicalVolume* expHall_log
    = new G4LogicalVolume(expHall_box,Air,"World",0,0,0);

  G4VPhysicalVolume* expHall_phys
    = new G4PVPlacement(0,G4ThreeVector(),expHall_log,"World",0,false,0);


  //  // ePHENIX stuff
  //  // -- Logical volume:
  //  G4VSolid *RICHOutSphereBoundary = new G4Orb("RICHOutSphereBoundary",
  //                                              geom.get_R_max());
  //  G4VSolid *RICHOutSphereBoundary_place = new G4DisplacedSolid(
  //                                                               "RICHOutSphereBoundary_place", RICHOutSphereBoundary, 0,
  //                                                               G4ThreeVector(0, geom.get_R_shift(), geom.get_z_shift()));
  //
  //  G4VSolid *RICHInnerSphereBoundary = new G4Orb("RICHInnerSphereBoundary",
  //                                                geom.get_R_frontwindow());
  //  G4VSolid *RICHInnerSphereBoundary_place = new G4DisplacedSolid(
  //                                                                 "RICHInnerSphereBoundary_place",
  //                                                                 RICHInnerSphereBoundary,
  //                                                                 0,
  //                                                                 G4ThreeVector(0,
  //                                                                               geom.get_R_shift() * geom.get_frontwindow_DisplaceRatio(),
  //                                                                               geom.get_z_shift() * geom.get_frontwindow_DisplaceRatio()));
  //
  //  G4VSolid *RICHConeBoundary = new G4Cons(
  //                                          "RICHConeBoundary", //
  //                                          geom.get_R_beam_pipe(),
  //                                          geom.get_z_shift() / 2
  //                                          * std::tan(2 * std::atan(std::exp(-geom.get_min_eta()))), //            G4double pRmin1, G4double pRmax1,
  //                                          geom.get_R_beam_pipe(),
  //                                          geom.get_cone_size_z()
  //                                          * std::tan(2 * std::atan(std::exp(-geom.get_min_eta()))), //            G4double pRmin2, G4double pRmax2,
  //                                          (geom.get_cone_size_z() - (geom.get_z_shift() / 2)) / 2, //            G4double pDz,
  //                                          0, 2 * pi //            G4double pSPhi, G4double pDPhi
  //                                          );
  //  G4VSolid *RICHConeBoundary_place = new G4DisplacedSolid(
  //                                                          "RICHConeBoundary_place",
  //                                                          RICHConeBoundary,
  //                                                          0,
  //                                                          G4ThreeVector(
  //                                                                        0,
  //                                                                        0,
  //                                                                        (geom.get_cone_size_z() - (geom.get_z_shift() / 2)) / 2
  //                                                                        + (geom.get_z_shift() / 2)));
  //
  //  G4VSolid *RICHSecBoundary = new G4Tubs("RICHSecBoundary", //
  //                                         geom.get_R_beam_pipe(), //            G4double pRMin,
  //                                         geom.get_cone_size_z(), //            G4double pRMax,
  //                                         geom.get_cone_size_z(), //            G4double pDz,
  //                                         pi / 2 - pi / geom.get_N_RICH_Sector(), //            G4double pSPhi,
  //                                         2 * pi / geom.get_N_RICH_Sector() //            G4double pDPhi
  //                                         );
  //
  //  G4VSolid *RICHSphereBoundary = new G4SubtractionSolid("RICHSphereBoundary",
  //                                                        RICHOutSphereBoundary_place, RICHInnerSphereBoundary_place);
  //  //      G4VSolid *RICHSecBox_ConeSphere = new G4IntersectionSolid(
  //  //              "RICHSecBox_ConeSphere", RICHSphereBoundary,
  //  //              RICHConeBoundary_place);
  //  G4VSolid *RICHSecBox_ConeSec = new G4IntersectionSolid("RICHSecBox_ConeSec",
  //                                                         RICHSecBoundary, RICHConeBoundary_place);
  //
  //  // RICH sector
  //  G4VSolid *RICHSecBox = new G4IntersectionSolid("RICHSecBox",
  //                                                 RICHSecBox_ConeSec, RICHSphereBoundary);
  //
  //  G4LogicalVolume *RICHSecLog = new G4LogicalVolume(RICHSecBox,
  //                                                    G4Material::GetMaterial(geom.get_RICH_gas_mat()), "RICHSecLogical",
  //                                                    0, 0, 0);
  //  RegisterLogicalVolume(RICHSecLog);
  //
  //  for (G4int sec = 0; sec < geom.get_N_RICH_Sector(); sec++)
  //    {
  //      G4RotateZ3D sec_rot(2 * pi / geom.get_N_RICH_Sector() * sec);
  //
  //      RegisterPhysicalVolume(new G4PVPlacement(sec_rot, RICHSecLog, "RICHSecPhysical", expHall_log,
  //                                               false, sec));
  //    }
  //  //      CLHEP::HepRotationZ sec_rot(0);
  //  //      G4RotationMatrix g4_sec_rot(sec_rot);
  //  //      new G4PVPlacement(transform1, "RICHSecPhysical", RICHSecLog, WorldPhys,
  //  //              false, 0);
  //
  //  //RICH mirror
  //  // LHCb website: 3mm thick Be base + 0.3mm glass surface layer coated with Al.
  //  // LHCb TDR: The mirrors are made of polished 6mm-thick glass coated by vacuum
  //  //           deposition with 900 nm of aluminium and overcoated with 200 nm of quartz.
  //  G4VSolid *RICHMirrorSphereBoundary = new G4Sphere(
  //                                                    "RICHMirrorSphereBoundary", //
  //                                                    geom.get_R_mirror_ref(),
  //                                                    geom.get_R_mirror_ref() + geom.get_dR_mirror(), //            G4double pRmin, G4double pRmax,
  //                                                    0, 2 * pi, //            G4double pSPhi, G4double pDPhi,
  //                                                    0, pi //            G4double pSTheta, G4double pDTheta
  //                                                    );
  //  G4VSolid *RICHMirrorSphereBoundary_place = new G4DisplacedSolid(
  //                                                                  "RICHMirrorSphereBoundary_place", RICHMirrorSphereBoundary, 0,
  //                                                                  G4ThreeVector(0, geom.get_R_shift(), geom.get_z_shift()));
  //  G4VSolid *RICHMirror = new G4IntersectionSolid("RICHMirror",
  //                                                 RICHSecBox_ConeSec, RICHMirrorSphereBoundary_place);
  //  G4LogicalVolume *RICHMirrorLog = new G4LogicalVolume(RICHMirror,
  //                                                       G4Material::GetMaterial(geom.get_RICH_Mirror_mat()),
  //                                                       "RICHMirrorLog");
  //  RegisterPhysicalVolume(new G4PVPlacement(0, G4ThreeVector(), RICHMirrorLog, "RICHMirrorPhysical",
  //                                           RICHSecLog, false, 0));
  //  RegisterLogicalVolume(RICHMirrorLog);
  //
  //  //RICH gas vessel window - back
  //  //LHCb TDR: The frame will be sealed to contain the C4F10 gas radiator. The vacuum chamber acts as part of the boundary to the gas volume.
  //  //                    Kapton foils of 150um thickness and 400 mm diameter will be glued to flanges on the vacuum chamber
  //  G4VSolid *RICHBackWindowSphereBoundary = new G4Sphere(
  //                                                        "RICHBackWindowSphereBoundary", //
  //                                                        geom.get_R_mirror_ref() + geom.get_dR_mirror()
  //                                                        + geom.get_dR_mirror_spt(),
  //                                                        geom.get_R_mirror_ref() + geom.get_dR_mirror()
  //                                                        + geom.get_dR_mirror_spt() + geom.get_dR_backwindow(), //            G4double pRmin, G4double pRmax,
  //                                                        0, 2 * pi, //            G4double pSPhi, G4double pDPhi,
  //                                                        0, pi //            G4double pSTheta, G4double pDTheta
  //                                                        );
  //  G4VSolid *RICHBackWindowSphereBoundary_place = new G4DisplacedSolid(
  //                                                                      "RICHBackWindowSphereBoundary_place", RICHBackWindowSphereBoundary,
  //                                                                      0, G4ThreeVector(0, geom.get_R_shift(), geom.get_z_shift()));
  //  G4VSolid *RICHBackWindow = new G4IntersectionSolid("RICHBackWindow",
  //                                                     RICHSecBox_ConeSec, RICHBackWindowSphereBoundary_place);
  //  G4LogicalVolume *RICHBackWindowLog = new G4LogicalVolume(RICHBackWindow,
  //                                                           G4Material::GetMaterial(geom.get_RICH_Gas_Window_mat()),
  //                                                           "RICHBackWindowLog");
  //  RegisterLogicalVolume(RICHBackWindowLog);
  //  RegisterPhysicalVolume(new G4PVPlacement(0, G4ThreeVector(), RICHBackWindowLog,
  //                                           "RICHBackWindowPhysical", RICHSecLog, false, 0));
  //
  //  //RICH gas vessel window - front
  //  //LHCb TDR: The frame will be sealed to contain the C4F10 gas radiator. The vacuum chamber acts as part of the boundary to the gas volume.
  //  //                    Kapton foils of 150um thickness and 400 mm diameter will be glued to flanges on the vacuum chamber
  //  G4VSolid *RICHFrontWindowSphereBoundary = new G4Sphere(
  //                                                         "RICHFrontWindowSphereBoundary", //
  //                                                         geom.get_R_frontwindow(),
  //                                                         geom.get_R_frontwindow() + geom.get_dR_frontwindow(), //            G4double pRmin, G4double pRmax,
  //                                                         0, 2 * pi, //            G4double pSPhi, G4double pDPhi,
  //                                                         0, pi //            G4double pSTheta, G4double pDTheta
  //                                                         );
  //  G4VSolid *RICHFrontWindowSphereBoundary_place = new G4DisplacedSolid(
  //                                                                       "RICHFrontWindowSphereBoundary_place",
  //                                                                       RICHFrontWindowSphereBoundary,
  //                                                                       0,
  //                                                                       G4ThreeVector(0,
  //                                                                                     geom.get_R_shift() * geom.get_frontwindow_DisplaceRatio(),
  //                                                                                     geom.get_z_shift() * geom.get_frontwindow_DisplaceRatio()));
  //  G4VSolid *RICHFrontWindow = new G4IntersectionSolid("RICHFrontWindow",
  //                                                      RICHSecBox_ConeSec, RICHFrontWindowSphereBoundary_place);
  //  G4LogicalVolume *RICHFrontWindowLog = new G4LogicalVolume(RICHFrontWindow,
  //                                                            G4Material::GetMaterial(geom.get_RICH_Gas_Window_mat()),
  //                                                            "RICHFrontWindowLog");
  //  RegisterLogicalVolume(RICHFrontWindowLog);
  //  RegisterPhysicalVolume(new G4PVPlacement(0, G4ThreeVector(), RICHFrontWindowLog,
  //                                           "RICHFrontWindowPhysical", RICHSecLog, false, 0));
  //
  //  // photon detector - HBD
  //  G4LogicalVolume* RICHHBDLog = Construct_HBD(RICHSecLog);
  //
  //  //G4VisAttributes
  //  G4VisAttributes * RICHSecTubeVisAtt = new G4VisAttributes(
  //                                                            G4Colour::White());
  //  RICHSecTubeVisAtt->SetForceWireframe(true);
  //  RICHSecTubeVisAtt->SetForceLineSegmentsPerCircle(50);
  //  //      RICHSecTubeVisAtt->SetForceSolid(true);
  //  RICHSecLog->SetVisAttributes(RICHSecTubeVisAtt);
  //
  //  G4VisAttributes * RICHMirrorVisAtt = new G4VisAttributes(G4Colour::Green());
  //  RICHMirrorVisAtt->SetForceWireframe(true);
  //  RICHMirrorVisAtt->SetForceSolid(true);
  //  RICHMirrorVisAtt->SetForceLineSegmentsPerCircle(50);
  //  RICHMirrorLog->SetVisAttributes(RICHMirrorVisAtt);
  //
  //  G4VisAttributes * RICHWindowVisAtt = new G4VisAttributes(G4Colour::Yellow());
  //  RICHWindowVisAtt->SetForceWireframe(true);
  //  RICHWindowVisAtt->SetForceSolid(true);
  //  RICHWindowVisAtt->SetForceLineSegmentsPerCircle(50);
  //  RICHBackWindowLog->SetVisAttributes(RICHWindowVisAtt);
  //  RICHFrontWindowLog->SetVisAttributes(RICHWindowVisAtt);
  //
  //  G4VisAttributes * RICHHBDVisAtt = new G4VisAttributes(G4Colour::Red());
  //  RICHHBDVisAtt->SetForceWireframe(true);
  //  RICHHBDVisAtt->SetForceSolid(true);
  //  RICHHBDLog->SetVisAttributes(RICHHBDVisAtt);
  //
  //  G4cout << "EicRichGemTbDetectorConstruction::Construct_RICH - " << map_log_vol.size()
  //         << " logical volume constructed" << G4endl;
  //  G4cout << "EicRichGemTbDetectorConstruction::Construct_RICH - " << map_phy_vol.size()
  //         << " physical volume constructed" << G4endl;
  //
  //
  //  // endePHENIX stuff

  return expHall_phys;

}

G4LogicalVolume * EicRichGemTbDetectorConstruction::RegisterLogicalVolume(
                                                                    G4LogicalVolume * v)
{
  //  if (!v)
  //    {
  //      G4cout
  //        << "EicRichGemTbDetectorConstruction::RegisterVolume - Error - invalid volume!"
  //        << G4endl;
  //      return v;
  //    }
  //  if (map_log_vol.find(v->GetName()) != map_log_vol.end())
  //    {
  //      G4cout
  //        << "EicRichGemTbDetectorConstruction::RegisterVolume - Warning - replacing "
  //        << v->GetName() << G4endl;
  //    }
  //
  //  map_log_vol[v->GetName()] = v;

  return v;
}

G4PVPlacement * EicRichGemTbDetectorConstruction::RegisterPhysicalVolume(
                                                                   G4PVPlacement * v)
{
  //  if (!v)
  //    {
  //      G4cout
  //        << "EicRichGemTbDetectorConstruction::RegisterPhysicalVolume - Error - invalid volume!"
  //        << G4endl;
  //      return v;
  //    }
  //
  //  phy_vol_idx_t id(v->GetName(),  v->GetCopyNo());
  //
  //  if (map_phy_vol.find(id) != map_phy_vol.end())
  //    {
  //      G4cout
  //        << "EicRichGemTbDetectorConstruction::RegisterPhysicalVolume - Warning - replacing "
  //        << v->GetName()<<"["<< v->GetCopyNo() <<"]" << G4endl;
  //    }
  //
  //  map_phy_vol[id] = v;
  return v;
}

G4LogicalVolume*
EicRichGemTbDetectorConstruction::Construct_HBD(G4LogicalVolume* RICHSecLog)
{

  const double HBD_thickness = geom.get_HBD_thickness();
  const int n_GEM_layers = geom.get_n_GEM_layers();

  assert(
         HBD_thickness < geom.get_dR_frontwindow_shrink() - geom.get_dR_frontwindow());
  // depth check

  G4VSolid *RICHHBDBox = new G4Cons("RICHHBDBox", //
                                    0, geom.get_RZ_Seg1_HBD() + geom.get_RZ_Seg2_HBD(), //            G4double pRmin1, G4double pRmax1,
                                    0, geom.get_RZ_Seg1_HBD() + geom.get_RZ_Seg2_HBD(), //            G4double pRmin2, G4double pRmax2,
                                    HBD_thickness / 2, //            G4double pDz,
                                    -geom.get_half_angle_HBD() + pi / 2, 2 * geom.get_half_angle_HBD() //            G4double pSPhi, G4double pDPhi
                                    );

  G4LogicalVolume *RICHHBDLog = new G4LogicalVolume(RICHHBDBox,
                                                    G4Material::GetMaterial(geom.get_RICH_gas_mat()), "RICHHBDLog");
  RegisterLogicalVolume(RICHHBDLog);

  G4Transform3D transform1 = G4Translate3D(0, geom.get_R_Tip_HBD(),
                                           geom.get_Z_Tip_HBD()) * G4RotateX3D(-geom.get_Rotation_HBD())
    * G4RotateY3D(pi) * G4TranslateZ3D(HBD_thickness / 2);
  //  G4Transform3D transform1 = G4TranslateZ3D(HBD_thickness / 2);

  RegisterPhysicalVolume(new G4PVPlacement(transform1, RICHHBDLog, "RICHHBDPhysical", RICHSecLog,
                                           false, 0));

  // Internal HBD structure
  // From doi:10.1016/j.nima.2011.04.015
  // Component Material X0 (cm) Thickness (cm) Area (%) Rad. Length (%)
  double current_z = -HBD_thickness / 2;

  //  Mesh SS 1.67 0.003 11.5 0.021
  Construct_HBD_Layers(RICHHBDLog, "Mesh", "Steel", current_z,
                       0.003 * cm * 11.5e-2);
  current_z += 0.003 * cm;

  //  GEM frames FR4 17.1 0.15x4 6.5 0.228
  Construct_HBD_Layers(RICHHBDLog, "Frame0", "G10", current_z,
                       0.15 * cm * 6.5e-2);
  current_z += 0.15 * cm;

  for (int gem = 1; gem <= n_GEM_layers; gem++)
    {
      stringstream sid;
      sid << gem;

      //  GEM Copper 1.43 0.0005x6 64 0.134
      Construct_HBD_Layers(RICHHBDLog,
                           G4String("GEMFrontCu") + G4String(sid.str()), "G4_Cu",
                           current_z, 0.0005 * cm * 64e-2);
      current_z += 0.0005 * cm;

      //  GEM Kapton 28.6 0.005x3 64 0.034
      Construct_HBD_Layers(RICHHBDLog,
                           G4String("GEMKapton") + G4String(sid.str()), "G4_KAPTON",
                           current_z, 0.005 * cm * 64e-2);
      current_z += 0.005 * cm;

      //  GEM Copper 1.43 0.0005x6 64 0.134
      Construct_HBD_Layers(RICHHBDLog,
                           G4String("GEMBackCu") + G4String(sid.str()), "G4_Cu", current_z,
                           0.0005 * cm * 64e-2);
      current_z += 0.0005 * cm;

      //  GEM frames FR4 17.1 0.15x4 6.5 0.228
      Construct_HBD_Layers(RICHHBDLog,
                           G4String("Frame") + G4String(sid.str()), "G10", current_z,
                           0.15 * cm * 6.5e-2);
      current_z += 0.15 * cm;
    }

  //  PCB Kapton 28.6 0.005 100 0.017
  Construct_HBD_Layers(RICHHBDLog, G4String("PCBKapton"), "G4_KAPTON",
                       current_z, 0.005 * cm * 100e-2);
  current_z += 0.005 * cm;

  //  PCB Copper 1.43 0.0005 80 0.028
  Construct_HBD_Layers(RICHHBDLog, G4String("PCBCu"), "G4_Cu", current_z,
                       0.0005 * cm * 80e-2);
  current_z += 0.0005 * cm;

  //  Facesheet FR4 17.1 0.025x2 100 0.292
  Construct_HBD_Layers(RICHHBDLog, "Facesheet", "G10", current_z,
                       0.025 * 2 * cm * 100e-2);
  current_z += 0.025 * 2 * cm;

  //  Panel core Honeycomb 8170 1.905 100 0.023 <- very thin-X0 stuff, ignore

  //  Total vessel 0.82
  //  Readout
  //  Readout board FR4/copper 17.1/1.43 0.05/0.001 100 0.367
  Construct_HBD_Layers(RICHHBDLog, G4String("ReadoutFR4"), "G10", current_z,
                       0.05 * cm * 100e-2);
  current_z += 0.05 * cm;
  Construct_HBD_Layers(RICHHBDLog, G4String("ReadoutCu"), "G4_Cu", current_z,
                       0.001 * cm * 100e-2);
  current_z += 0.001 * cm;

  //  Preamps + sockets Copper 1.43 0.0005 100 0.66
  //  Total readout 1.03
  Construct_HBD_Layers(RICHHBDLog, G4String("SocketsCu"), "G4_Cu", current_z,
                       0.0005 * cm * 100e-2);
  current_z += 0.0005 * cm;

  assert(current_z<HBD_thickness / 2);
  //boundary check

  return RICHHBDLog;
}

G4LogicalVolume*
EicRichGemTbDetectorConstruction::Construct_HBD_Layers(G4LogicalVolume* RICHHBDLog,
                                                 const G4String name, const G4String material, const double start_z,
                                                 const double thickness)
{
  G4VSolid * box = new G4Cons(G4String("RICHHBD") + name + G4String("Box"), //
                              0, geom.get_RZ_Seg1_HBD() + geom.get_RZ_Seg2_HBD(), //            G4double pRmin1, G4double pRmax1,
                              0, geom.get_RZ_Seg1_HBD() + geom.get_RZ_Seg2_HBD(), //            G4double pRmin2, G4double pRmax2,
                              thickness / 2, //            G4double pDz,
                              -geom.get_half_angle_HBD() + pi / 2, 2 * geom.get_half_angle_HBD() //            G4double pSPhi, G4double pDPhi
                              );

  G4LogicalVolume *Log = new G4LogicalVolume(box,
                                             G4Material::GetMaterial(material) //
                                             , G4String("RICHHBD") + name + G4String("Log"));
  RegisterLogicalVolume(Log);

  RegisterPhysicalVolume(new G4PVPlacement(0, G4ThreeVector(0, 0, start_z + thickness / 2), Log,
                                           G4String("RICHHBD") + name + G4String("Physical"), RICHHBDLog,
                                           false, 0));

  return Log;
}

void EicRichGemTb_Geometry::SetDefault()
{

  N_RICH_Sector = 8;
  min_eta = 1;
  R_beam_pipe = 3 * cm;
  z_shift = 100 * cm;
  R_shift = 40 * cm;
  frontwindow_DisplaceRatio = .85; // Displace R,Z and radius simultainously
  dR_frontwindow_shrink = 2 * cm;
  R_mirror_ref = 200 * cm;
  dR_mirror = 0.6 * cm;
  dR_mirror_spt = 1.4 * cm;
  dR_backwindow = 150 * um;
  dR_frontwindow = 150 * um;

  HBD_thickness = 1.5 * cm;
  n_GEM_layers = 4;

  RICH_gas_mat = "CF4";
  RICH_Mirror_mat = "G4_Pyrex_Glass";
  RICH_Gas_Window_mat = "G4_KAPTON";

}

double EicRichGemTb_Geometry::get_R_frontwindow() const
{
  return (R_mirror_ref / 2
          + sqrt(z_shift * z_shift + R_shift * R_shift)
          * (1 - frontwindow_DisplaceRatio)) - dR_frontwindow_shrink;
}

double EicRichGemTb_Geometry::get_half_angle_HBD() const
{
  return atan(
              tan(pi / N_RICH_Sector)
              / sqrt(1 + R_shift * R_shift / z_shift / z_shift));
}

double EicRichGemTb_Geometry::get_RZ_Seg1_HBD() const
{
  return R_shift
    * (R_mirror_ref + 2 * sqrt(z_shift * z_shift + R_shift * R_shift))
    / 4 / z_shift;
}

double EicRichGemTb_Geometry::get_RZ_Seg2_HBD() const
{
  return ((R_mirror_ref + 2 * sqrt(pow(z_shift, 2) + pow(R_shift, 2)))
          * tan(2 * atan(exp(-min_eta)) - atan(R_shift / z_shift))) / 4.;
}

double EicRichGemTb_Geometry::get_R_Tip_HBD() const
{
  const double l = sqrt(z_shift * z_shift + R_shift * R_shift)
    + R_mirror_ref / 2;

  return l * sin(get_Rotation_HBD())
    - get_RZ_Seg1_HBD() * cos(get_Rotation_HBD());
}

double EicRichGemTb_Geometry::get_Z_Tip_HBD() const
{
  const double l = sqrt(z_shift * z_shift + R_shift * R_shift)
    + R_mirror_ref / 2;

  return l * cos(get_Rotation_HBD())
    + get_RZ_Seg1_HBD() * sin(get_Rotation_HBD());
}

double EicRichGemTb_Geometry::get_Rotation_HBD() const
{
  return atan(R_shift / z_shift);
}
