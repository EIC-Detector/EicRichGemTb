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



  // Now for the material properties of Surfaces
  //
  //Front (reflecting surface of RichTb Mirror)

  G4OpticalSurface * OpticalMirrorSurface = new G4OpticalSurface("OpticalMirrorSurface");

  //string mirrortype = "LHCb";
  G4String mirrortype = "G4example";

  if ( mirrortype == "LHCb" )
    { // LHCb mirror
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

      OpticalMirrorSurface->SetType(dielectric_metal);
      OpticalMirrorSurface->SetModel(glisur);
      OpticalMirrorSurface->SetFinish(polished);

      G4MaterialPropertiesTable* OpticalMirrorSurfaceMPT = new G4MaterialPropertiesTable();

      OpticalMirrorSurfaceMPT->AddProperty("RINDEX",
                                           PhotonMomentumRefl,
                                           MirrorQuRefIndex,
                                           NumPhotonRichMirrorReflWaveLengthBins);
      OpticalMirrorSurfaceMPT->AddProperty("REFLECTIVITY",
                                           PhotonMomentumRefl,
                                           PhotWaveRefl,
                                           NumPhotonRichMirrorReflWaveLengthBins);
      OpticalMirrorSurfaceMPT->AddProperty("EFFICIENCY",
                                           PhotonMomentumRefl,
                                           PhotReflEff,
                                           NumPhotonRichMirrorReflWaveLengthBins);

      OpticalMirrorSurface->SetMaterialPropertiesTable(OpticalMirrorSurfaceMPT);

      OpticalMirrorSurface->DumpInfo();
    } // end LHCb mirror

  else if ( mirrortype == "G4example" )
    { // test from G4 example
      const G4int NUM = 2;

      G4double pp[NUM] = {2.038*eV, 4.144*eV};
      G4double specularlobe[NUM] = {0.3, 0.3};
      G4double specularspike[NUM] = {0.2, 0.2};
      G4double backscatter[NUM] = {0.1, 0.1};
      G4double rindex[NUM] = {1.35, 1.40};
      G4double reflectivity[NUM] = {1.0, 1.0};
      G4double efficiency[NUM] = {0.0, 0.0};

      G4MaterialPropertiesTable* SMPT = new G4MaterialPropertiesTable();

      SMPT -> AddProperty("RINDEX",pp,rindex,NUM);
      SMPT -> AddProperty("SPECULARLOBECONSTANT",pp,specularlobe,NUM);
      SMPT -> AddProperty("SPECULARSPIKECONSTANT",pp,specularspike,NUM);
      SMPT -> AddProperty("BACKSCATTERCONSTANT",pp,backscatter,NUM);
      SMPT -> AddProperty("REFLECTIVITY",pp,reflectivity,NUM);
      SMPT -> AddProperty("EFFICIENCY",pp,efficiency,NUM);

      OpticalMirrorSurface->SetMaterialPropertiesTable(SMPT);

      OpticalMirrorSurface->DumpInfo();
    } // end test

  return;

}


EicRichGemTbMaterial::~EicRichGemTbMaterial(){ }
