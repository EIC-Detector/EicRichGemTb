#include <iostream>
#include <cmath>

#include "globals.hh"

#include "EicRichGemTbMaterial.hh"

#include "G4Isotope.hh"
#include "G4Element.hh"
#include "G4ElementTable.hh"
#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4UnitsTable.hh"

#include "G4OpticalSurface.hh"
//#include "G4LogicalBorderSurface.hh"
//#include "G4LogicalSkinSurface.hh"
//#include "G4OpBoundaryProcess.hh"
#include "G4MaterialPropertyVector.hh"

//#include "EicRichGemTbMaterialParameters.hh"
//#include "EicRichGemTbGeometryParameters.hh"


EicRichGemTbMaterial::EicRichGemTbMaterial(){

  G4double a; //atomic mass
  G4double z;  //atomic number
  G4double density; //density
  G4double temperature; //temperature
  G4double pressure; //pressure

  G4int nelements;
  G4int natoms;

  // Elements
  //
  H  = new G4Element("Hydrogen", "H", z=1 , a=1.01*g/mole);
  N  = new G4Element("Nitrogen", "N", z=7 , a=14.01*g/mole);
  O  = new G4Element("Oxygen"  , "O", z=8 , a=16.00*g/mole);
  C  = new G4Element("Carbon", "C", z=6 , a=12.01*g/mole);
  F  = new G4Element("Fluorine", "F", z=9 , a=19.00*g/mole);
  Si = new G4Element("Silicon", "Si", z=14, a=28.09*g/mole);


  //Vacuum
  //
  density=universe_mean_density;
  a=1.01*g/mole;
  pressure=1.e-19*pascal;
  temperature=0.1*kelvin;

  Vacuum = new G4Material("Galactic",density,nelements=1,kStateGas,temperature,pressure);
  Vacuum->AddElement(H,natoms=1);

  // Air
  //
  AmbientAir = new G4Material("Air", density=1.29*mg/cm3, nelements=2);
  AmbientAir->AddElement(N, 70.*perCent);
  AmbientAir->AddElement(O, 30.*perCent);

  // CF4
  //

  // LHCb; no data available at room temp and pressure;
  CF4 = new G4Material("CF4", density=0.003884*g/cm3, nelements=2, kStateGas, temperature=273.*kelvin, pressure=1.0*atmosphere);

  // http://www.slac.stanford.edu/pubs/icfa/summer98/paper3/paper3.pdf : density=0.00393*g/cm^3, T = 20degC

  // http://encyclopedia.airliquide.com/Encyclopedia.asp?GasID=61#GeneralData :
  // Molecular weight  : 88.01 g/mol
  // Gas density (1.013 bar and 15 °C (59 °F)) : 3.72 kg/m3

  // example
  //fCF4 = new G4Material("CF4", density=3.72*g/cm3, nelements=2);

  CF4->AddElement(C, 1);
  CF4->AddElement(F, 4);

  // Aluminum
  //
  Aluminum = new G4Material("Aluminum", z=13 , a=26.98*g/mole , density=2.7*g/cm3);

  G4cout << *(G4Material::GetMaterialTable()) << G4endl;

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

  SiO2MirrorQuartz = new G4Material("MirrorQuartz", density=2.200*g/cm3, nelements=2);
  SiO2MirrorQuartz->AddElement(Si,natoms=1);
  SiO2MirrorQuartz->AddElement(O,natoms=2);

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

  G4MaterialPropertiesTable* fMPT_cf4 = new G4MaterialPropertiesTable();

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

  CF4->SetMaterialPropertiesTable(fMPT_cf4);

  // Set the Birks Constant for the Cf4 scintillator
  //fCF4->GetIonisation()->SetBirksConstant(0.126*mm/MeV);

  return;

}
//  G4double a,z,density;  //a=mass of a mole;
//  // z=mean number of protons;
//  G4String name,symbol;
//  //G4int isz,isn;    //isz= number of protons in an isotope;
//  //isn= number of nucleons in an isotope;
//
//
//  G4int numel,natoms;  //numel=Number of elements constituting a material.
//  G4double fractionmass;
//  G4double temperature, pressure;
//  // G4double FactorOne=1.0;
//  G4UnitDefinition::BuildUnitsTable();
//
//  //PhotonEnergy
//  G4int ibin=0;
//  G4double PhotonEnergyStep=(PhotonMaxEnergy-PhotonMinEnergy)/
//    NumPhotWaveLengthBins;
//  G4double* PhotonMomentum=new G4double[NumPhotWaveLengthBins];
//  for (ibin=0; ibin<NumPhotWaveLengthBins; ibin++){
//    PhotonMomentum[ibin]=PhotonMinEnergy+PhotonEnergyStep*ibin;
//  }
//
//  G4cout << "\nNow Define Elements ..\n" <<G4endl;
//
//  // Nitrogen
//
//  a=14.01*g/mole;
//  G4Element* elN = new G4Element(name="Nitrogen",
//                                 symbol="N", z=7., a);
//
//  //Oxygen
//
//  a=16.00*g/mole;
//  G4Element* elO = new G4Element(name="Oxygen",
//                                 symbol="O", z=8., a);
//
//  //Hydrogen
//
//  a=1.01*g/mole;
//  G4Element* elH = new G4Element(name="Hydrogen",
//                                 symbol="H",z=1.,a);
//
//  //Carbon
//
//  a=12.01*g/mole;
//  G4Element* elC = new G4Element(name="Carbon",
//                                 symbol="C",z=6.,a);
//
//  //Silicon
//
//  a=28.09*g/mole;
//  G4Element* elSi = new G4Element(name="Silicon",
//                                  symbol="Si",z=14.,a);
//  //Fluorine
//  a=18.998*g/mole;
//  G4Element* elF = new G4Element(name="Fluorine",
//                                 symbol="F",z=9.,a);
//  //Aluminum
//  a=26.98*g/mole;
//  G4Element* elAL =new G4Element(name="Aluminium",
//                                 symbol="Al",z=13.,a);
//
//  //Sodium
//  a=22.99*g/mole;
//  G4Element* elNa = new G4Element(name="Sodium",
//                                  symbol="Na",z=11.,a);
//
//  //Potassium
//  a=39.10*g/mole;
//  G4Element* elK = new G4Element(name="Potassium",
//                                 symbol="K",z=19.,a);
//
//  //Cesium
//
//  // a=132.91*g/mole;
//  // G4Element* elCs = new G4Element(name="Cesium",
//  //                                symbol="Cs",z=55.,a);
//
//  //Antimony
//
//  a=121.76*g/mole;
//  G4Element* elSb = new G4Element(name="Antimony",
//                                  symbol="Sb",z=51.,a);
//
//
//  //Define Materials
//  G4cout << "\nNow Define Materials ..\n" <<G4endl;
//  //
//
//  //Air at 20 degree C and 1 atm for the ambiet air.
//  // Also Air as  a radiator material for inside the tubes.
//  //--
//  density = 1.205e-03*g/cm3;
//  pressure=1.*atmosphere;
//  temperature=293.*kelvin;
//  G4Material* Air = new G4Material(name="Air ", density, numel=2,
//                                   kStateGas,temperature,pressure);
//  Air->AddElement(elN, fractionmass=0.7);
//  Air->AddElement(elO, fractionmass=0.3);
//
//  G4double* AirAbsorpLength=new G4double[NumPhotWaveLengthBins];
//  G4double* AirRindex=new G4double[NumPhotWaveLengthBins];
//
//  for (ibin=0; ibin<NumPhotWaveLengthBins; ibin++){
//    AirAbsorpLength[ibin]=1.E32*mm;
//    AirRindex[ibin]=1.000273;
//  }
//  G4MaterialPropertiesTable* AirMPT =
//    new G4MaterialPropertiesTable();
//
//  AirMPT->AddProperty("ABSLENGTH",PhotonMomentum,
//                      AirAbsorpLength,NumPhotWaveLengthBins);
//
//  Air->SetMaterialPropertiesTable(AirMPT);
//  EicRichGemTbAmbientAir = Air;
//
//
//  density = 1.205e-03*g/cm3;
//  pressure=1.*atmosphere;
//  temperature=293.*kelvin;
//  G4Material* TAir = new G4Material(name="TAir ", density, numel=2,
//                                    kStateGas,temperature,pressure);
//  TAir->AddElement(elN, fractionmass=0.7);
//  TAir->AddElement(elO, fractionmass=0.3);
//
//  G4double* TAirAbsorpLength=new G4double[NumPhotWaveLengthBins];
//  G4double* TAirRindex=new G4double[NumPhotWaveLengthBins];
//
//  for (ibin=0; ibin<NumPhotWaveLengthBins; ibin++){
//    TAirAbsorpLength[ibin]=1.E32*mm;
//    TAirRindex[ibin]=1.000273;
//  }
//  G4MaterialPropertiesTable* TAirMPT =
//    new G4MaterialPropertiesTable();
//
//  TAirMPT->AddProperty("ABSLENGTH",PhotonMomentum,
//                       TAirAbsorpLength,NumPhotWaveLengthBins);
//
//  TAirMPT->AddProperty("RINDEX", PhotonMomentum,
//                       AirRindex,NumPhotWaveLengthBins);
//
//  TAir->SetMaterialPropertiesTable(TAirMPT);
//  EicRichGemTbTubeAir = TAir;
//
//  //Nitrogen gas.
//
//  density = 0.8073e-03*g/cm3;
//  pressure = RConfig -> getPressureN2();
//  temperature = RConfig ->getTemperatureN2();
//
//  G4Material* NitrogenGas = new G4Material(name="NitrogenGas ",
//                                           density, numel=1,
//                                           kStateGas,temperature,pressure);
//  NitrogenGas->AddElement(elN, natoms=2);
//
//  G4double* NitrogenGasAbsorpLength=new G4double[NumPhotWaveLengthBins];
//  G4double* NitrogenGasRindex=new G4double[NumPhotWaveLengthBins];
//  G4double* NitrogenGasPhotW=new G4double[NumPhotWaveLengthBins];
//
//
//  std::vector<G4double>N2RefInd= InitN2RefIndex(pressure,temperature);
//  std::vector<G4double>N2RefPhotW=InitN2RefPhotW();
//
//  for (ibin=0; ibin<NumPhotWaveLengthBins; ibin++){
//    NitrogenGasAbsorpLength[ibin]=1.E32*mm;
//
//    NitrogenGasRindex[ibin]=N2RefInd[ibin];
//    NitrogenGasPhotW[ibin]=N2RefPhotW[ibin];
//
//  }
//  G4MaterialPropertiesTable* NitrogenGasMPT =
//    new G4MaterialPropertiesTable();
//
//  NitrogenGasMPT->AddProperty("ABSLENGTH",NitrogenGasPhotW,
//                              NitrogenGasAbsorpLength,NumPhotWaveLengthBins);
//
//  NitrogenGasMPT->AddProperty("RINDEX", NitrogenGasPhotW,
//                              NitrogenGasRindex,NumPhotWaveLengthBins);
//
//  NitrogenGas->SetMaterialPropertiesTable(NitrogenGasMPT);
//  EicRichGemTbNitrogenGas = NitrogenGas;
//
//  //Water
//  density=1.000*g/cm3;
//  G4Material* H2O = new G4Material(name="Water",density,numel=2);
//  H2O->AddElement(elH,natoms=2);
//  H2O->AddElement(elO,natoms=1);
//
//  G4double* H2OAbsorpLength=new G4double[NumPhotWaveLengthBins];
//  G4double* H2ORindex=new G4double[NumPhotWaveLengthBins];
//
//  for (ibin=0; ibin<NumPhotWaveLengthBins; ibin++){
//    H2OAbsorpLength[ibin]=1.E32*mm;
//    H2ORindex[ibin]=1.33;
//  }
//
//
//  G4MaterialPropertiesTable* H2OMPT =
//    new G4MaterialPropertiesTable();
//
//  H2OMPT->AddProperty("ABSLENGTH",PhotonMomentum,
//                      H2OAbsorpLength,NumPhotWaveLengthBins);
//
//  H2OMPT->AddProperty("RINDEX", PhotonMomentum,
//                      H2ORindex,NumPhotWaveLengthBins);
//
//  H2O->SetMaterialPropertiesTable(H2OMPT);
//
//
//  EicRichGemTbH2O=H2O;
//  //Sio2
//  //There is a quartz for the mirror and
//  //another quartz which is used in aerogel and
//  // yet another quartz used for the quartz window.
//  //Mirrorquartz
//
//  density=2.200*g/cm3;
//  G4Material* SiO2MirrorQuartz = new G4Material(name="MirrorQuartz",
//                                                density,numel=2);
//  SiO2MirrorQuartz->AddElement(elSi,natoms=1);
//  SiO2MirrorQuartz->AddElement(elO,natoms=2);
//
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
//  EicRichGemTbMirrorQuartz=SiO2MirrorQuartz;
//
//  density=2.200*g/cm3;
//  G4Material* SiO2AerogelQuartz = new G4Material(name="AerogelQuartz",
//                                                 density,numel=2);
//  SiO2AerogelQuartz->AddElement(elSi,natoms=1);
//  SiO2AerogelQuartz->AddElement(elO,natoms=2);
//
//  // QuartzWindow Quartz
//  density=2.200*g/cm3;
//  G4Material* WindowQuartz = new G4Material(name="WindowQuartz",
//                                            density,numel=2);
//  WindowQuartz->AddElement(elSi,natoms=1);
//  WindowQuartz->AddElement(elO,natoms=2);
//  G4double* WindowQuartzRindex=new G4double[NumPhotWaveLengthBins];
//  G4double* WindowQuartzAbsorpLength=new G4double[NumPhotWaveLengthBins];
//  for (ibin=0; ibin<NumPhotWaveLengthBins; ibin++){
//    WindowQuartzAbsorpLength[ibin]=1.E32*mm;
//    WindowQuartzRindex[ibin]=1.4;
//  }
//  G4MaterialPropertiesTable* WindowQuartzMPT =
//    new G4MaterialPropertiesTable();
//
//  WindowQuartzMPT->AddProperty("ABSLENGTH",PhotonMomentum,
//                               WindowQuartzAbsorpLength,NumPhotWaveLengthBins);
//
//  WindowQuartzMPT->AddProperty("RINDEX", PhotonMomentum,
//                               WindowQuartzRindex,NumPhotWaveLengthBins);
//  WindowQuartz->SetMaterialPropertiesTable(WindowQuartzMPT);
//
//  EicRichGemTbQuartzWindowMaterial=WindowQuartz;
//  //for now this is kept to be same as the hpdquartz window.
//  density=2.200*g/cm3;
//  G4Material* HpdWindowQuartz = new G4Material(name="HpdWindowQuartz",
//                                               density,numel=2);
//  HpdWindowQuartz->AddElement(elSi,natoms=1);
//  HpdWindowQuartz->AddElement(elO,natoms=2);
//  G4double* HpdWindowQuartzRindex=new G4double[NumPhotWaveLengthBins];
//  G4double* HpdWindowQuartzAbsorpLength=new G4double[NumPhotWaveLengthBins];
//  for (ibin=0; ibin<NumPhotWaveLengthBins; ibin++){
//    HpdWindowQuartzAbsorpLength[ibin]=1.E32*mm;
//    HpdWindowQuartzRindex[ibin]=1.40;
//  }
//  G4MaterialPropertiesTable* HpdWindowQuartzMPT =
//    new G4MaterialPropertiesTable();
//
//  HpdWindowQuartzMPT->AddProperty("ABSLENGTH",PhotonMomentum,
//                                  HpdWindowQuartzAbsorpLength,NumPhotWaveLengthBins);
//
//  HpdWindowQuartzMPT->AddProperty("RINDEX", PhotonMomentum,
//                                  HpdWindowQuartzRindex,NumPhotWaveLengthBins);
//  HpdWindowQuartz->SetMaterialPropertiesTable(HpdWindowQuartzMPT);
//
//  HpdQuartzWindowMaterial=HpdWindowQuartz;
//
//  // Borosilcate window of the Pad Hpd
//  // for now kept same as the other Hpd.
//  density=2.200*g/cm3;
//  G4Material* PadHpdWindowQuartz = new G4Material(name="PadHpdWindowQuartz",
//                                                  density,numel=2);
//  PadHpdWindowQuartz->AddElement(elSi,natoms=1);
//  PadHpdWindowQuartz->AddElement(elO,natoms=2);
//  G4double* PadHpdWindowQuartzRindex=new G4double[NumPhotWaveLengthBins];
//  G4double* PadHpdWindowQuartzAbsorpLength=new G4double[NumPhotWaveLengthBins];
//  for (ibin=0; ibin<NumPhotWaveLengthBins; ibin++){
//    PadHpdWindowQuartzAbsorpLength[ibin]=1.E32*mm;
//    PadHpdWindowQuartzRindex[ibin]=1.40;
//  }
//  G4MaterialPropertiesTable* PadHpdWindowQuartzMPT =
//    new G4MaterialPropertiesTable();
//
//  PadHpdWindowQuartzMPT->AddProperty("ABSLENGTH",PhotonMomentum,
//                                     PadHpdWindowQuartzAbsorpLength,NumPhotWaveLengthBins);
//
//  PadHpdWindowQuartzMPT->AddProperty("RINDEX", PhotonMomentum,
//                                     PadHpdWindowQuartzRindex,NumPhotWaveLengthBins);
//  PadHpdWindowQuartz->SetMaterialPropertiesTable(PadHpdWindowQuartzMPT);
//
//  PadHpdQuartzWindowMaterial=PadHpdWindowQuartz;
//  //
//  //
//
//  G4int filterNumberThisRun=RConfig->GetFilterTNumber();
//  // now for the filter material glass d263
//  density=2.200*g/cm3;
//  G4Material* GlassD263 = new G4Material(name= FilterTypeString[0],
//                                         density,numel=2);
//  GlassD263->AddElement(elSi,natoms=1);
//  GlassD263->AddElement(elO,natoms=2);
//
//  if(filterNumberThisRun >= 0 ) {
//
//    //in the following the +2 is to match the materialproperty bins
//    // for the various materials, to avoid the tons of printout from G4.
//    // Please see the explanation below for getting the abosorption
//    // length of aerogel. The same comments apply here as well.
//    // Essentially the measured transmission input here is a combination of
//    // the bulk absorption and the fresnel surface loss. One needs to
//    // decouple them. Here a partial attempt is made to avoid
//    // modifying the G4OpBoundary process. SE. 15-11-2002.
//    G4double* GlassD263Rindex=new G4double[NumPhotBinGlassD263Trans+2];
//    G4double* GlassD263AbsorpLength=new G4double[NumPhotBinGlassD263Trans+2];
//    G4double* GlassD263MomValue = new G4double[NumPhotBinGlassD263Trans+2];
//    G4double* currBulkTransFilter = new G4double[NumPhotBinGlassD263Trans+2];
//    FilterTrData* CurFil = RConfig->GetFilterTrData();
//    std::vector<G4double>GlassD263TransWL = CurFil-> GetTransWL();
//    std::vector<G4double>GlassD263Transmis = CurFil->GetTransTotValue();
//    G4double FilterHalfZ= CurFil->GetCurFilterThickness();
//
//    for (ibin=0; ibin<NumPhotBinGlassD263Trans+2; ibin++){
//
//
//      GlassD263Rindex[ibin]=RefIndexGlassD263;
//      if(ibin > 0 && ibin < NumPhotBinGlassD263Trans+1 ){
//        //now using the formula trans=std::exp(-thickness/absorplength).
//        G4int ibina=ibin-1;
//
//        if(GlassD263TransWL[ibina] > 0.0 ) {
//          GlassD263MomValue[ibin]= PhotMomWaveConv*eV/GlassD263TransWL[ibina];
//        }
//        if(GlassD263Transmis[ibina] >0.0 ) {
//          // G4double currentfilterRefIndex= GlassD263Rindex[ibin];
//          G4double currentAdjacentMediumRefIndex=NitrogenNominalRefIndex;
//          // the following needs to be improved in the future
//          // to have a binary search and
//          // interpolation between the adjacent
//          // array elements etc. SE 15-11-2002.
//          for(size_t ibinr=0; ibinr<N2RefPhotW.size()-1 ; ibinr++){
//            G4double currMomA=GlassD263MomValue[ibin];
//            if(currMomA >= N2RefPhotW[ibinr] && currMomA <= N2RefPhotW[ibinr+1]){
//              currentAdjacentMediumRefIndex=N2RefInd[ibinr];
//            }
//
//          }
//
//          if( GlassD263Transmis[ibina] > 0.01 ) {
//            currBulkTransFilter[ibin]=
//              GetCurrentBulkTrans(GlassD263Rindex[ibin],
//                                  currentAdjacentMediumRefIndex,
//                                  GlassD263Transmis[ibina]);
//          } else {
//            currBulkTransFilter[ibin]=GlassD263Transmis[ibina];
//
//          }
//          if(currBulkTransFilter[ibin] > 0.0 &&
//             currBulkTransFilter[ibin] < 0.9995 ) {
//            GlassD263AbsorpLength[ibin]=
//              -(2.0*FilterHalfZ)/(std::log(currBulkTransFilter[ibin]));
//          }else if (currBulkTransFilter[ibin]== 0.0 ) {
//            GlassD263AbsorpLength[ibin]=FilterHalfZ/1.0E32;
//          }else {
//            GlassD263AbsorpLength[ibin]=DBL_MAX;
//          }
//        }else {
//
//          GlassD263AbsorpLength[ibin]=FilterHalfZ/1.0E32;
//        }
//      }
//
//    }
//    GlassD263MomValue[0]=PhotonMaxEnergy;
//    GlassD263AbsorpLength[0]=GlassD263AbsorpLength[1];
//    currBulkTransFilter[0]=currBulkTransFilter[1];
//
//    G4int mbin=NumPhotBinGlassD263Trans+1;
//    GlassD263MomValue[mbin]=PhotonMinEnergy;
//    GlassD263AbsorpLength[mbin]=GlassD263AbsorpLength[mbin-1];
//    currBulkTransFilter[mbin]=currBulkTransFilter[mbin-1];
//
//    G4MaterialPropertiesTable* GlassD263MPT =
//      new G4MaterialPropertiesTable();
//
//    GlassD263MPT->AddProperty("ABSLENGTH",GlassD263MomValue,
//                              GlassD263AbsorpLength,NumPhotBinGlassD263Trans+2);
//
//    GlassD263MPT->AddProperty("RINDEX",GlassD263MomValue,
//                              GlassD263Rindex,NumPhotBinGlassD263Trans+2);
//
//    GlassD263->SetMaterialPropertiesTable(GlassD263MPT);
//  }
//
//  GlassD263FilterMaterial=GlassD263;
//  EicRichGemTbFilterMaterial[0]=GlassD263;
//  //for the G4Example only 1 filter type is used.
//  G4cout << " Now Define Aerogel .." <<G4endl;
//
//
//  //Aerogel upto five types considered so far.
//  // in the G4example the same type is repeated 5 times.
//  //Now for TypeA
//
//  density=0.200*g/cm3;
//
//  G4Material* AerogTypeA =
//    new G4Material(name=AerogelTypeString[0], density, numel=2);
//  AerogTypeA->AddMaterial(SiO2AerogelQuartz, fractionmass=97.0*perCent);
//  AerogTypeA->AddMaterial(H2O, fractionmass=3.0*perCent);
//
//
//  G4double* AerogTypeARindex=new G4double[NumPhotWaveLengthBins];
//  G4double* AerogTypeAAbsorpLength=new G4double[NumPhotWaveLengthBins];
//  G4double* AerogTypeARScatLength = new G4double[NumPhotWaveLengthBins];
//  G4double* currentAgelTrans = new G4double[NumPhotWaveLengthBins];
//
//  std::vector<G4double>AerogelTypeASLength = GetAerogelRScatLength(AerogelTypeA);
//  G4int AerogNumber=0;
//  G4double AerogelLength=GetCurAerogelLength(AerogNumber);
//  G4double MaxTotTransmission=AerogelTypeATotTrans;
//  // Unfortunately the transmission measurement values only give the
//  // total transmission which includes the loss within aerogel
//  // and the Fresnel loss at the surface. In order to
//  // partially decouple this, the approximate loss at the
//  // the surface is calculated using the ref index of the
//  // aerogel and its surroundings. Then this is added to the
//  // measured transmission to get the transmission in the bulk of
//  // aerogel. This is then converted to an absorption length.
//  // In a more accurate implementation the loss at the surface
//  // should be calculated using a more precise formula. It is
//  // difficult since we do not know the direction of the photons
//  // at this point.
//  // One possibility is to modify the G4opBoundaryProcess
//  // for this, since we do know the direction of the photons by then.
//  // This is not done for this G4example, but only in the LHCb implementation.
//  // SE 15-11-2002.
//  // The aerogel is inside a volume made of Nitrogen
//
//  for (ibin=0; ibin<NumPhotWaveLengthBins; ibin++){
//    AerogTypeARindex[ibin]= ConvertAgelRIndex(PhotonMomentum[ibin],0);
//    AerogTypeARScatLength[ibin]=AerogelTypeASLength[ibin];
//    // G4double photwl = PhotMomWaveConv/ (PhotonMomentum[ibin]/eV);
//
//    G4double currentAgelRefIndex=  AerogTypeARindex[ibin];
//    G4double currentNeighbourRefIndex= N2RefInd[ibin];
//    currentAgelTrans[ibin]=
//      GetCurrentBulkTrans( currentAgelRefIndex,
//                           currentNeighbourRefIndex,MaxTotTransmission);
//    //now using the formula trans=std::exp(-thickness/absorplength)
//    // to get the absorplength.
//
//    if( currentAgelTrans[ibin] > 0.0 && currentAgelTrans[ibin] < 0.9995) {
//      AerogTypeAAbsorpLength[ibin]=
//        -(AerogelLength)/(std::log( currentAgelTrans[ibin]));
//    }else if (currentAgelTrans[ibin] == 0.0) {
//
//      AerogTypeAAbsorpLength[ibin]=AerogelLength/1.0E32;
//    }else {
//
//      AerogTypeAAbsorpLength[ibin]=DBL_MAX;
//    }
//
//  }
//
//  G4MaterialPropertiesTable* AerogTypeAMPT =
//    new G4MaterialPropertiesTable();
//
//  AerogTypeAMPT->AddProperty("ABSLENGTH",PhotonMomentum,
//                             AerogTypeAAbsorpLength,NumPhotWaveLengthBins);
//
//
//  AerogTypeAMPT->AddProperty("RAYLEIGH",PhotonMomentum,
//                             AerogTypeARScatLength,NumPhotWaveLengthBins);
//
//  AerogTypeAMPT->AddProperty("RINDEX", PhotonMomentum,
//                             AerogTypeARindex,NumPhotWaveLengthBins);
//
//  AerogTypeA->SetMaterialPropertiesTable(AerogTypeAMPT);
//
//
//  EicRichGemTbAerogelTypeA = AerogTypeA;
//  EicRichGemTbAerogelMaterial[0] = AerogTypeA;
//  // In the G4example the same type is repeated 5 times.
//  // in the LHCb implementation 5 types of aerogel materials used.
//  //Now for Aerogel TypeB
//
//  EicRichGemTbAerogelTypeB = AerogTypeA;
//  EicRichGemTbAerogelMaterial[1] = AerogTypeA;
//
//  //Now for aerogel TypeC
//
//  EicRichGemTbAerogelTypeC = AerogTypeA;
//  EicRichGemTbAerogelMaterial[2] = AerogTypeA;
//
//  //Now for aerogel TypeD
//
//  EicRichGemTbAerogelTypeD = AerogTypeA;
//  EicRichGemTbAerogelMaterial[3] = AerogTypeA;
//
//  //Now for aerogel Type E
//
//
//  EicRichGemTbAerogelTypeE = AerogTypeA;
//  EicRichGemTbAerogelMaterial[4] = AerogTypeA;
//
//
//
//  //Bialkali Photocathode
//
//  //the following numbers on the property of the BiAlkali Photocathode
//  // may not be accurate.
//  //Some number is is jut put in for initial program test purposes.
//  density=0.100*g/cm3;
//  G4Material* BiAlkaliPhCathode = new G4Material(name="BiAlkaliPhCathode",
//                                                 density, numel=3);
//  BiAlkaliPhCathode->AddElement(elNa, fractionmass=37.5*perCent);
//  BiAlkaliPhCathode->AddElement(elK, fractionmass=37.5*perCent);
//  BiAlkaliPhCathode->AddElement(elSb, fractionmass=25.0*perCent);
//
//  //for now  properties for the ph cathode material.
//
//  G4double* BiAlkaliPhCathodeRindex=new G4double[NumPhotWaveLengthBins];
//  G4double* BiAlkaliPhCathodeAbsorpLength=new G4double[NumPhotWaveLengthBins];
//  G4double CathLen=PhotoCathodeThickness;
//  G4double CathTrans=PhCathodeNominalTransmission;
//  G4double CathAbsorpLen;
//  if(CathTrans > 0.0 && CathTrans < 0.9995 ) {
//    CathAbsorpLen =  -(CathLen)/(std::log(CathTrans));
//  }else if (CathTrans > 0.0) {
//    CathAbsorpLen  = CathLen/1.0E32;
//  }else {
//    CathAbsorpLen  = DBL_MAX;
//  }
//
//  for (ibin=0; ibin<NumPhotWaveLengthBins; ibin++){
//    BiAlkaliPhCathodeAbsorpLength[ibin]=CathAbsorpLen;
//    BiAlkaliPhCathodeRindex[ibin]=1.40;
//  }
//  G4MaterialPropertiesTable* BiAlkaliPhCathodeMPT =
//    new G4MaterialPropertiesTable();
//  BiAlkaliPhCathodeMPT->AddProperty("ABSLENGTH",PhotonMomentum,
//                                    BiAlkaliPhCathodeAbsorpLength,NumPhotWaveLengthBins);
//
//  BiAlkaliPhCathodeMPT->AddProperty("RINDEX", PhotonMomentum,
//                                    BiAlkaliPhCathodeRindex,NumPhotWaveLengthBins);
//  BiAlkaliPhCathode->SetMaterialPropertiesTable(BiAlkaliPhCathodeMPT);
//  PadHpdPhCathodeMaterial=BiAlkaliPhCathode;
//
//  //CF4
//  //no data available at room temp and pressure;
//  density=0.003884*g/cm3;
//  temperature=273.*kelvin;
//  pressure=1.0*atmosphere;
//  a=88.01*g/mole;
//
//  G4Material* CF4 =new G4Material(name="CF4",density,numel=2,
//                                  kStateGas,temperature,pressure);
//  CF4->AddElement(elC,natoms=1);
//  CF4->AddElement(elF,natoms=4);
//  // Sellmeir coef to be added.
//  EicRichGemTbCF4=CF4;
//
//  G4cout << "\nNowDefineVacuum ..\n" <<G4endl;
//
//
//  G4double* VacAbsorpLength=new G4double[NumPhotWaveLengthBins];
//  G4double* VacRindex=new G4double[NumPhotWaveLengthBins];
//
//  for (ibin=0; ibin<NumPhotWaveLengthBins; ibin++){
//    VacAbsorpLength[ibin]=1.E32*mm;
//    // the following ref index is just artifical, just to
//    // avoid the refraction between nitrogen gas and hpd master.
//    VacRindex[ibin]=1.000273;
//  }
//  G4MaterialPropertiesTable* VacMPT =
//    new G4MaterialPropertiesTable();
//
//  VacMPT->AddProperty("ABSLENGTH",PhotonMomentum,
//                      VacAbsorpLength,NumPhotWaveLengthBins);
//  VacMPT->AddProperty("RINDEX", PhotonMomentum,
//                      VacRindex,NumPhotWaveLengthBins);
//  vacuum->SetMaterialPropertiesTable(VacMPT);
//
//  EicRichGemTbVacuum=vacuum;
//
//  //beamgas
//  //
//  density=1.e-5*g/cm3;
//  pressure=2.e-2*bar;
//  temperature=STP_Temperature;
//  G4Material* beamgas = new G4Material(name="Beamgas",density,numel=1,
//                                       kStateGas,temperature,pressure);
//  beamgas->AddMaterial(Air,fractionmass=1.); // beware that air is at 20 deg;
//
//  //
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
//  EicRichGemTbAluminium=Aluminium;
//  //PlasticAg , this is used as a wrap of aerogel and as upstream holder
//  // for aerogel frame. For now use same properties as that of Aluminium.
//  // this is just an opaque material.
//
//  density=2.7*g/cm3;
//  G4Material* PlasticAg =new G4Material(name="PlasticAg",density,numel=1);
//  PlasticAg->AddElement(elAL,natoms=1);
//
//  G4double* PlasticAgAbsorpLength=new G4double[NumPhotWaveLengthBins];
//
//  for (ibin=0; ibin<NumPhotWaveLengthBins; ibin++){
//    PlasticAgAbsorpLength[ibin]=0.0*mm;
//  }
//
//  G4MaterialPropertiesTable* PlasticAgMPT =
//    new G4MaterialPropertiesTable();
//
//  PlasticAgMPT->AddProperty("ABSLENGTH",PhotonMomentum,
//                            PlasticAgAbsorpLength,NumPhotWaveLengthBins);
//
//  PlasticAg->SetMaterialPropertiesTable(PlasticAgMPT);
//  EicRichGemTbPlasticAg=PlasticAg;
//  // Kovar
//  density=2.7*g/cm3;
//  G4Material* Kovar =new G4Material(name="Kovar",density,numel=1);
//  Kovar->AddElement(elAL,natoms=1);
//
//  G4double* KovarAbsorpLength=new G4double[NumPhotWaveLengthBins];
//
//  for (ibin=0; ibin<NumPhotWaveLengthBins; ibin++){
//    KovarAbsorpLength[ibin]=0.0*mm;
//  }
//
//  G4MaterialPropertiesTable* KovarMPT =
//    new G4MaterialPropertiesTable();
//
//  KovarMPT->AddProperty("ABSLENGTH",PhotonMomentum,
//                        KovarAbsorpLength,NumPhotWaveLengthBins);
//
//  Kovar->SetMaterialPropertiesTable(KovarMPT);
//  HpdTubeMaterial=Kovar;
//
//  // Silicon
//
//  density=2.33*g/cm3;
//  G4Material* Silicon =new G4Material(name="Silicon",density,numel=1);
//  Silicon->AddElement(elSi,natoms=1);
//
//  G4double* SiliconAbsorpLength=new G4double[NumPhotWaveLengthBins];
//
//  for (ibin=0; ibin<NumPhotWaveLengthBins; ibin++){
//    SiliconAbsorpLength[ibin]=0.0*mm;
//  }
//
//  G4MaterialPropertiesTable* SiliconMPT =
//    new G4MaterialPropertiesTable();
//
//  SiliconMPT->AddProperty("ABSLENGTH",PhotonMomentum,
//                          SiliconAbsorpLength,NumPhotWaveLengthBins);
//
//  Silicon->SetMaterialPropertiesTable(SiliconMPT);
//  HpdSiDetMaterial=Silicon;
//
//  // Silicon coating made of Si02.
//
//  density=2.33*g/cm3;
//  G4Material* SiliconCoating =new G4Material(name="SilCoat",density,numel=2);
//  SiliconCoating->AddElement(elSi,natoms=1);
//  SiliconCoating->AddElement(elO,natoms=2);
//
//  G4double* SiliconCoatingAbsorpLength=new G4double[NumPhotWaveLengthBins];
//  // G4double* SiliconCoatingRindex=new G4double[NumPhotWaveLengthBins];
//
//  for (ibin=0; ibin<NumPhotWaveLengthBins; ibin++){
//    SiliconCoatingAbsorpLength[ibin]=0.0001*mm;
//
//  }
//
//  G4MaterialPropertiesTable* SiliconCoatingMPT =
//    new G4MaterialPropertiesTable();
//
//  SiliconCoatingMPT->AddProperty("ABSLENGTH",PhotonMomentum,
//                                 SiliconCoatingAbsorpLength,NumPhotWaveLengthBins);
//
//  SiliconCoating->SetMaterialPropertiesTable(SiliconCoatingMPT);
//  HpdSiCoatingMaterial=SiliconCoating;
//
//  //
//  // Now for the material properties of Surfaces
//  //
//  //
//  //
//  //Front (reflecting surface of EicRichGemTb Mirror)
//
//  // First define wavelength in nm.
//  //For now assume that all segments have the same reflectivity.
//  // Hence the reflectivity is defined outside the loop of the
//  // the number of segments.
//  //Only the front surface is created.
//  // The abosorption length is set to a small value just to
//  // avoid photons exiting from the back of the mirror.
//  // the efficiency is for the absorption process.
//
//
//  G4double* PhotonMomentumRefl
//    =new G4double[NumPhotonRichMirrorReflWaveLengthBins];
//  G4double* PhotWaveRefl =
//    new  G4double[NumPhotonRichMirrorReflWaveLengthBins];
//  G4double* PhotReflEff =new  G4double[NumPhotonRichMirrorReflWaveLengthBins];
//  G4double* MirrorQuRefIndex
//    =new  G4double[NumPhotonRichMirrorReflWaveLengthBins];
//
//  for (ibin=0; ibin<NumPhotonRichMirrorReflWaveLengthBins; ibin++){
//    PhotonMomentumRefl[ibin]=PhotMomWaveConv*eV/ PhotonWavelengthRefl[ibin];
//    PhotWaveRefl[ibin]=  EicRichGemTbMirrorReflectivity[ibin];
//    PhotReflEff[ibin]= EicRichGemTbMirrorEfficiency[ibin];
//    //the following lines to avoid reflection at the mirror.
//
//    MirrorQuRefIndex[ibin] = 1.40;
//  }
//
//  G4OpticalSurface * OpEicRichGemTbMirrorSurface =
//    new G4OpticalSurface("EicRichGemTbMirrorSurface");
//
//  OpEicRichGemTbMirrorSurface->SetType(dielectric_metal);
//  OpEicRichGemTbMirrorSurface->SetFinish(polished);
//  OpEicRichGemTbMirrorSurface->SetModel(glisur);
//  G4MaterialPropertiesTable* OpEicRichGemTbMirrorSurfaceMPT =
//    new G4MaterialPropertiesTable();
//
//  OpEicRichGemTbMirrorSurfaceMPT->AddProperty("REFLECTIVITY",
//                                        PhotonMomentumRefl,
//                                        PhotWaveRefl,
//                                        NumPhotonRichMirrorReflWaveLengthBins);
//  OpEicRichGemTbMirrorSurfaceMPT->AddProperty("EFFICIENCY",
//                                        PhotonMomentumRefl,
//                                        PhotReflEff,
//                                        NumPhotonRichMirrorReflWaveLengthBins);
//  OpEicRichGemTbMirrorSurfaceMPT->AddProperty("RINDEX",
//                                        PhotonMomentumRefl,
//                                        MirrorQuRefIndex,
//                                        NumPhotonRichMirrorReflWaveLengthBins);
//
//  OpEicRichGemTbMirrorSurface->SetMaterialPropertiesTable(OpEicRichGemTbMirrorSurfaceMPT);
//  EicRichGemTbOpticalMirrorSurface=OpEicRichGemTbMirrorSurface;
//
//
//  //  OpEicRichGemTbMirrorSurface->DumpInfo();
//
//  // Now for the Surface of the Vessel Enclosure.
//
//
//  G4OpticalSurface * OpEicRichGemTbEnclosureSurface =
//    new G4OpticalSurface("EicRichGemTbEnclosureSurface");
//  OpEicRichGemTbEnclosureSurface->SetType(dielectric_metal);
//  OpEicRichGemTbEnclosureSurface->SetFinish(polished);
//  OpEicRichGemTbEnclosureSurface->SetModel(glisur);
//
//  G4double NumPhotonRichEnclosureSurfaceWaveLengthBins=10;
//  G4double  EicRichGemTbEnclosureSurfaceReflectivity[]=
//    {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
//
//  G4double EicRichGemTbEnclosureSurfaceEfficiency[]=
//    {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
//  G4double RichEnclosureSurfacePhotMom[]=
//    {1.0*eV,2.0*eV, 3.0*eV,4.0*eV,5.0*eV,6.0*eV,7.0*eV,8.0*eV,
//     9.0*eV,10.0*eV};
//
//  G4MaterialPropertiesTable* OpEicRichGemTbEnclosureSurfaceMPT =
//    new G4MaterialPropertiesTable();
//
//  OpEicRichGemTbEnclosureSurfaceMPT->AddProperty("REFLECTIVITY",
//                                           RichEnclosureSurfacePhotMom,
//                                           EicRichGemTbEnclosureSurfaceReflectivity,
//                                           static_cast<int>(NumPhotonRichEnclosureSurfaceWaveLengthBins));
//  OpEicRichGemTbEnclosureSurfaceMPT->AddProperty("EFFICIENCY",
//                                           RichEnclosureSurfacePhotMom,
//                                           EicRichGemTbEnclosureSurfaceEfficiency,
//                                           static_cast<int>(NumPhotonRichEnclosureSurfaceWaveLengthBins));
//
//  OpEicRichGemTbEnclosureSurface->
//    SetMaterialPropertiesTable(OpEicRichGemTbEnclosureSurfaceMPT);
//
//  EicRichGemTbOpticalEnclosureSurface=OpEicRichGemTbEnclosureSurface;
//
//  //Now for the surface between the TAir and Quartz Window of the HPD
//
//  G4OpticalSurface * OpHpdQuartzWTSurface =
//    new G4OpticalSurface("HpdQuartzWTSurface");
//  OpHpdQuartzWTSurface->SetType(dielectric_dielectric);
//  OpHpdQuartzWTSurface->SetFinish(polished);
//  OpHpdQuartzWTSurface->SetModel(glisur);
//  //OpHpdQuartzWTSurface->SetModel(unified);
//
//  G4double NumPhotonHpdQuartzWTSurfaceWaveLengthBins=10;
//  G4double  HpdQuartzWTSurfaceReflectivity[]=
//    {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
//
//  G4double HpdQuartzWTSurfaceEfficiency[]=
//    {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
//
//  G4double HpdQuartzWTSurfacePhotMom[]=
//    {1.0*eV,2.0*eV, 3.0*eV,4.0*eV,5.0*eV,6.0*eV,7.0*eV,8.0*eV,
//     9.0*eV,10.0*eV};
//
//  G4MaterialPropertiesTable* OpHpdQuartzWTSurfaceMPT =
//    new G4MaterialPropertiesTable();
//
//
//  OpHpdQuartzWTSurfaceMPT->AddProperty("REFLECTIVITY",
//                                       HpdQuartzWTSurfacePhotMom,
//                                       HpdQuartzWTSurfaceReflectivity,
//                                       static_cast<int>(NumPhotonHpdQuartzWTSurfaceWaveLengthBins));
//  OpHpdQuartzWTSurfaceMPT->AddProperty("EFFICIENCY",
//                                       HpdQuartzWTSurfacePhotMom,
//                                       HpdQuartzWTSurfaceEfficiency,
//                                       static_cast<int>(NumPhotonHpdQuartzWTSurfaceWaveLengthBins));
//
//  OpHpdQuartzWTSurface->
//    SetMaterialPropertiesTable(OpHpdQuartzWTSurfaceMPT);
//
//  HpdTQuartzWSurface=OpHpdQuartzWTSurface;
//
//
//
//  //Now for the surface between the Quartz Window and Ph cathode of the HPD
//
//  G4OpticalSurface * OpHpdQuartzWPSurface =
//    new G4OpticalSurface("HpdQuartzWPSurface");
//  OpHpdQuartzWPSurface->SetType(dielectric_dielectric);
//  OpHpdQuartzWPSurface->SetFinish(polished);
//  OpHpdQuartzWPSurface->SetModel(glisur);
//
//  G4double NumPhotonHpdQuartzWPSurfaceWaveLengthBins=10;
//  G4double  HpdQuartzWPSurfaceReflectivity[]=
//    {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
//
//
//  G4double HpdQuartzWPSurfaceEfficiency[]=
//    {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
//
//  G4double HpdQuartzWPSurfacePhotMom[]=
//    {1.0*eV,2.0*eV, 3.0*eV,4.0*eV,5.0*eV,6.0*eV,7.0*eV,8.0*eV,
//     9.0*eV,10.0*eV};
//
//  G4MaterialPropertiesTable* OpHpdQuartzWPSurfaceMPT =
//    new G4MaterialPropertiesTable();
//
//
//  OpHpdQuartzWPSurfaceMPT->AddProperty("REFLECTIVITY",
//                                       HpdQuartzWPSurfacePhotMom,
//                                       HpdQuartzWPSurfaceReflectivity,
//                                       static_cast<int>(NumPhotonHpdQuartzWPSurfaceWaveLengthBins));
//  OpHpdQuartzWPSurfaceMPT->AddProperty("EFFICIENCY",
//                                       HpdQuartzWPSurfacePhotMom,
//                                       HpdQuartzWPSurfaceEfficiency,
//                                       static_cast<int>(NumPhotonHpdQuartzWPSurfaceWaveLengthBins));
//
//  OpHpdQuartzWPSurface->
//    SetMaterialPropertiesTable(OpHpdQuartzWPSurfaceMPT);
//
//  HpdQuartzWPhCathodeSurface=OpHpdQuartzWPSurface;
//
//
//
//  //Now for the skin surface of the PhCathode so that photons do
//  // not come out of the Photocathode.
//  // Changed to dielectric-dielectric so that photons DO come out
//  // of the photocathode. SE 26-9-01.
//
//  G4OpticalSurface * OpPhCathodeSurface =
//    new G4OpticalSurface("PhCathodeSurface");
//
//  OpPhCathodeSurface->SetType(dielectric_dielectric);
//  OpPhCathodeSurface->SetFinish(polished);
//  OpPhCathodeSurface->SetModel(glisur);
//
//
//  G4double NumPhotonPhCathodeSurfaceWaveLengthBins=10;
//  G4double  PhCathodeSurfaceReflectivity[]=
//    {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
//
//
//  G4double PhCathodeSurfaceEfficiency[]=
//    {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
//
//  G4double PhCathodeSurfacePhotMom[]=
//    {1.0*eV,2.0*eV, 3.0*eV,4.0*eV,5.0*eV,6.0*eV,7.0*eV,8.0*eV,
//     9.0*eV,10.0*eV};
//
//  G4MaterialPropertiesTable* OpPhCathodeSurfaceMPT =
//    new G4MaterialPropertiesTable();
//
//
//  OpPhCathodeSurfaceMPT->AddProperty("REFLECTIVITY",
//                                     PhCathodeSurfacePhotMom,
//                                     PhCathodeSurfaceReflectivity,
//                                     static_cast<int>(NumPhotonPhCathodeSurfaceWaveLengthBins));
//  OpPhCathodeSurfaceMPT->AddProperty("EFFICIENCY",
//                                     PhCathodeSurfacePhotMom,
//                                     PhCathodeSurfaceEfficiency,
//                                     static_cast<int>(NumPhotonPhCathodeSurfaceWaveLengthBins));
//
//  OpPhCathodeSurface->
//    SetMaterialPropertiesTable(OpPhCathodeSurfaceMPT);
//
//  PhCathodeSkinSurface=OpPhCathodeSurface;
//  PhCathodeBorderSurface=OpPhCathodeSurface;
//
//
//
//
//  //Now for the surface between Interior of HPD and Silicon Coating.
//
//  G4OpticalSurface * OpHpdSiCoatSurface =
//    new G4OpticalSurface("HpdSiCoatSurface");
//  OpHpdSiCoatSurface->SetType(dielectric_metal);
//  OpHpdSiCoatSurface->SetFinish(polished);
//  OpHpdSiCoatSurface->SetModel(glisur);
//
//
//  G4double NumPhotonHpdSiCoatSurfaceWaveLengthBins=10;
//
//  G4double  HpdSiCoatSurfaceReflectivity[]=
//    {0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9};
//
//  G4double HpdSiCoatSurfaceEfficiency[]=
//    {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
//
//  G4double HpdSiCoatSurfacePhotMom[]=
//    {1.0*eV,2.0*eV, 3.0*eV,4.0*eV,5.0*eV,6.0*eV,7.0*eV,8.0*eV,
//     9.0*eV,10.0*eV};
//  G4double HpdSiCoatSurfaceRefInd[]=
//    {1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.4};
//
//
//  G4MaterialPropertiesTable* OpHpdSiCoatSurfaceMPT =
//    new G4MaterialPropertiesTable();
//
//
//  OpHpdSiCoatSurfaceMPT->AddProperty("REFLECTIVITY",
//                                     HpdSiCoatSurfacePhotMom,
//                                     HpdSiCoatSurfaceReflectivity,
//                                     static_cast<int>(NumPhotonHpdSiCoatSurfaceWaveLengthBins));
//  OpHpdSiCoatSurfaceMPT->AddProperty("EFFICIENCY",
//                                     HpdSiCoatSurfacePhotMom,
//                                     HpdSiCoatSurfaceEfficiency,
//                                     static_cast<int>(NumPhotonHpdSiCoatSurfaceWaveLengthBins));
//  OpHpdSiCoatSurfaceMPT->AddProperty("RINDEX",
//                                     HpdSiCoatSurfacePhotMom,
//                                     HpdSiCoatSurfaceRefInd,
//                                     static_cast<int>(NumPhotonHpdSiCoatSurfaceWaveLengthBins));
//
//  OpHpdSiCoatSurface->
//    SetMaterialPropertiesTable(OpHpdSiCoatSurfaceMPT);
//
//  HpdSiCoatSurface=OpHpdSiCoatSurface;
//
//
//
//  // Now for the Surface of the MetalTube of HPD.
//
//
//  G4OpticalSurface * OpEicRichGemTbHpdMetalSurface =
//    new G4OpticalSurface("EicRichGemTbHpdMetalSurface");
//  OpEicRichGemTbHpdMetalSurface->SetType(dielectric_metal);
//  OpEicRichGemTbHpdMetalSurface->SetFinish(polished);
//  OpEicRichGemTbHpdMetalSurface->SetModel(glisur);
//
//  G4double NumPhotonHpdMetalSurfaceWaveLengthBins=10;
//  G4double  RichHpdMetalSurfaceReflectivity[]=
//    {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
//
//  G4double RichHpdMetalSurfaceEfficiency[]=
//    {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
//  G4double RichHpdMetalSurfacePhotMom[]=
//    {1.0*eV,2.0*eV, 3.0*eV,4.0*eV,5.0*eV,6.0*eV,7.0*eV,8.0*eV,
//     9.0*eV,10.0*eV};
//
//  G4MaterialPropertiesTable* OpEicRichGemTbHpdMetalSurfaceMPT =
//    new G4MaterialPropertiesTable();
//
//  OpEicRichGemTbHpdMetalSurfaceMPT->AddProperty("REFLECTIVITY",
//                                          RichHpdMetalSurfacePhotMom,
//                                          RichHpdMetalSurfaceReflectivity,
//                                          static_cast<int>(NumPhotonHpdMetalSurfaceWaveLengthBins));
//  OpEicRichGemTbHpdMetalSurfaceMPT->AddProperty("EFFICIENCY",
//                                          RichHpdMetalSurfacePhotMom,
//                                          RichHpdMetalSurfaceEfficiency,
//                                          static_cast<int>(NumPhotonHpdMetalSurfaceWaveLengthBins));
//
//  OpEicRichGemTbHpdMetalSurface->
//    SetMaterialPropertiesTable(OpEicRichGemTbHpdMetalSurfaceMPT);
//
//  EicRichGemTbOpticalHpdMetalSurface=OpEicRichGemTbHpdMetalSurface;
//
//
//  //Now for the surface of the Filter
//
//  G4OpticalSurface * OpEicRichGemTbFilterSurface =
//    new G4OpticalSurface("EicRichGemTbFilterSurface");
//  OpEicRichGemTbFilterSurface->SetType(dielectric_dielectric);
//  OpEicRichGemTbFilterSurface->SetFinish(polished);
//  OpEicRichGemTbFilterSurface->SetModel(glisur);
//
//
//
//  if(filterNumberThisRun >= 0 ) {
//
//    G4int FilterNumbins=NumPhotonEicRichGemTbFilterSurfaceWaveLengthBins;
//
//
//    G4double*  FilterReflectivity = new G4double(FilterNumbins);
//    G4double* FilterEff =new G4double(FilterNumbins);
//    G4double* FilterPhotMom =new G4double(FilterNumbins);
//
//
//    for(G4int ibinf =0 ; ibinf < FilterNumbins; ibinf++ ){
//      FilterReflectivity[ibinf]= EicRichGemTbFilterSurfaceReflectivity[ibinf];
//      FilterEff[ibinf]= EicRichGemTbFilterSurfaceEfficiency[ibinf];
//      FilterPhotMom[ibinf]= EicRichGemTbFilterSurfacePhotMom[ibinf];
//
//      //       G4MaterialPropertiesTable* OpEicRichGemTbFilterSurfaceMPT =
//      //                         new G4MaterialPropertiesTable();
//
//    }
//    EicRichGemTbOpticalFilterSurface=OpEicRichGemTbFilterSurface;
//
//  }
//
//  delete [] PhotonMomentum;
//  delete [] AirAbsorpLength;
//  delete [] AirRindex;
//  delete [] MirrorQuartzRindex;
//  delete [] MirrorQuartzAbsorpLength;
//  delete [] WindowQuartzRindex;
//  delete [] WindowQuartzAbsorpLength;
//  delete [] AluminiumAbsorpLength;
//  delete [] KovarAbsorpLength;
//  delete [] PhotonMomentumRefl;
//
//
//}

EicRichGemTbMaterial::~EicRichGemTbMaterial(){ ; }



