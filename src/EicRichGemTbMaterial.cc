#include <iostream>
#include <cmath>

#include "globals.hh"

#include "EicRichGemTbMaterial.hh"

#include "G4Isotope.hh"
#include "G4NistManager.hh"
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

  // G4 database on G4Elements
  G4NistManager* manager = G4NistManager::Instance();
  manager->SetVerbose(0);

  // Define Elements
  H  = manager->FindOrBuildElement(1);
  N  = manager->FindOrBuildElement(7);
  O  = manager->FindOrBuildElement(8);
  C  = manager->FindOrBuildElement(6);
  F  = manager->FindOrBuildElement(9);
  Mg = manager->FindOrBuildElement(12);
  Si = manager->FindOrBuildElement(14);
  Ar = manager->FindOrBuildElement(18);
  I  = manager->FindOrBuildElement(53);
  Cs = manager->FindOrBuildElement(55);

  //Vacuum
  density=universe_mean_density;
  a=1.01*g/mole;
  pressure=1.e-19*pascal;
  temperature=0.1*kelvin;

  Vacuum = new G4Material("Galactic",density,nelements=1,kStateGas,temperature,pressure);
  Vacuum->AddElement(H,natoms=1);

  // Air
  CO2 = new G4Material("C02", density = 1.977*kg/m3, nelements = 2);
  CO2->AddElement(C, natoms = 1);
  CO2->AddElement(O, natoms = 2);

  AmbientAir = new G4Material("Air", density = 1.293*mg/cm3, nelements=4);
  AmbientAir->AddElement(N, 75.53*perCent);
  AmbientAir->AddElement(O, 23.14*perCent);
  AmbientAir->AddElement(Ar, 1.29*perCent);
  AmbientAir->AddMaterial(CO2, 0.04*perCent);

  // CF4
  //

  CF4 = new G4Material("CF4", density=0.003884*g/cm3, nelements = 2, kStateGas);
  CF4->AddElement(C, 1);
  CF4->AddElement(F, 4);

  // CsI
  //
  CsI = new G4Material("CsI", density= 4.534*g/cm3, nelements=2);
  CsI->AddElement(Cs, natoms=1);
  CsI->AddElement(I , natoms=1);
  //CsI->GetIonisation()->SetMeanExcitationEnergy(553.1*eV);


  // Stainless Steel
  //
  StainlessSteel = new G4Material("StainlessSteel",density=8.96*g/cm3,nelements=1);
  StainlessSteel->AddElement(O, 1.);


  // Aluminum
  //
  Aluminum = new G4Material("Aluminum", z=13 , a=26.98*g/mole , density=2.7*g/cm3);

  //MgF2 
  //
  MgF2 = new G4Material("MgF2", density = 3.148*g/cm3, nelements = 2);
  MgF2->AddElement(Mg, natoms = 1);
  MgF2->AddElement(F, natoms = 2);

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
  //
  SiO2MirrorQuartz = new G4Material("MirrorQuartz", density=2.200*g/cm3, nelements=2);
  SiO2MirrorQuartz->AddElement(Si,natoms=1);
  SiO2MirrorQuartz->AddElement(O,natoms=2);

  G4int NumPhotWaveLengthBins = 2;

  G4double MirrorQuartzRindex[2]={1.35,1.35};
  G4double MirrorQuartzAbsorpLength[2]={1,1};

  G4double PhotonMomentum[2]={0.1*eV,1000*eV};

  //  G4double* MirrorQuartzRindex=new G4double[NumPhotWaveLengthBins];
  //  G4double* MirrorQuartzAbsorpLength=new G4double[NumPhotWaveLengthBins];

  //  for (ibin=0; ibin<NumPhotWaveLengthBins; ibin++){
  //    MirrorQuartzAbsorpLength[ibin]=0.01*mm;
  //
  //  }
  G4MaterialPropertiesTable* MirrorQuartzMPT =
    new G4MaterialPropertiesTable();


  MirrorQuartzMPT->AddProperty("ABSLENGTH",PhotonMomentum,
                               MirrorQuartzAbsorpLength,NumPhotWaveLengthBins);


  SiO2MirrorQuartz->SetMaterialPropertiesTable(MirrorQuartzMPT);


  // Generate & Add Material Properties Table
  //

  // CF4
  //

  // Need to check values. These values are from ExN06 for Water.
  //
  const G4int cf4_nEntries = 32;

  // Photon wavelengths 120 - 200 nm
  //

  G4double cf4_PhotonWavelength[cf4_nEntries]  =
    {200.000000000000*nm, 197.419354838710*nm, 194.838709677419*nm, 192.258064516129*nm,
     189.677419354839*nm, 187.096774193548*nm, 184.516129032258*nm, 181.935483870968*nm,
     179.354838709677*nm, 176.774193548387*nm, 174.193548387097*nm, 171.612903225807*nm,
     169.032258064516*nm, 166.451612903226*nm, 163.870967741936*nm, 161.290322580645*nm,
     158.709677419355*nm, 156.129032258065*nm, 153.548387096774*nm, 150.967741935484*nm,
     148.387096774194*nm, 145.806451612903*nm, 143.225806451613*nm, 140.645161290323*nm,
     138.064516129032*nm, 135.483870967742*nm, 132.903225806452*nm, 130.322580645161*nm,
     127.741935483871*nm, 125.161290322581*nm, 122.580645161290*nm, 120.000000000000*nm};

  // Photon energies
  //

  G4double* cf4_PhotonEnergy = new G4double[cf4_nEntries];
  G4int i;
  for(i=0;i<cf4_nEntries;i++)
    {
      cf4_PhotonEnergy[i] = 1.23984193*eV*1000*nm/cf4_PhotonWavelength[i];
    }

  //Refractive Index of CF4
  //

 G4double cf4_RefractiveIndex[cf4_nEntries]  =
   {1.00052729574769, 1.00052876876669, 1.00053030945866, 1.00053192216804,
    1.00053361160381, 1.00053538287781, 1.00053724154787, 1.00053919366665,
    1.00054124583687, 1.00054340527413, 1.00054567987839, 1.00054807831564,
    1.00055061011147, 1.00055328575860, 1.00055611684072, 1.00055911617590,
    1.00056229798296, 1.00056567807526, 1.00056927408754, 1.00057310574219,
    1.00057719516344, 1.00058156724978, 1.00058625011743, 1.00059127563130,
    1.00059668004417, 1.00060250477064, 1.00060879733031, 1.00061561250462,
    1.00062301376631, 1.00063107505917, 1.00063988303241, 1.00064953987108};

  G4MaterialPropertiesTable* fMPT_cf4 = new G4MaterialPropertiesTable();

  fMPT_cf4->AddProperty("RINDEX", cf4_PhotonEnergy, cf4_RefractiveIndex, cf4_nEntries )->SetSpline(true);

  fMPT_cf4->AddConstProperty("SCINTILLATIONYIELD",50./MeV);

  CF4->SetMaterialPropertiesTable(fMPT_cf4);

  // Now for the material properties of Surfaces
  //

  //Front (reflecting surface of RichTb Mirror)
  //

  OpticalMirrorSurface = new G4OpticalSurface("OpticalMirrorSurface");

  //G4String mirrortype = "LHCb";
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

      //      OpticalMirrorSurface->SetType(dielectric_metal);
      OpticalMirrorSurface->SetType(dielectric_dielectric);
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
      const G4int NUM = 32;
      
      // Photon energies
      //
      G4double* pp = new G4double[NUM];
      G4int i;
      for(i = 0; i < NUM; i++)
	{
	  pp[i] = cf4_PhotonEnergy[i];
	}

      //G4double specularlobe[NUM] = {0.3, 0.3};
      //G4double specularspike[NUM] = {0.2, 0.2};
      //G4double backscatter[NUM] = {0.1, 0.1};
      //G4double rindex[NUM] = {1.35, 1.40};
      G4double rindex[NUM] = 
	{1.4, 1.4, 1.4, 1.4,
	 1.4, 1.4, 1.4, 1.4,
	 1.4, 1.4, 1.4, 1.4,
	 1.4, 1.4, 1.4, 1.4,
	 1.4, 1.4, 1.4, 1.4,
	 1.4, 1.4, 1.4, 1.4,
	 1.4, 1.4, 1.4, 1.4,
	 1.4, 1.4, 1.4, 1.4};

      G4double reflectivity[NUM] = 
	{0.801070, 0.820321, 0.826738, 0.824599,
	 0.817112, 0.807487, 0.794652, 0.783957,
	 0.776471, 0.768984, 0.767914, 0.768984,
	 0.766845, 0.762567, 0.759358, 0.757219,
	 0.757219, 0.758289, 0.762567, 0.766845,
	 0.775401, 0.782888, 0.787166, 0.789305,
	 0.789305, 0.793583, 0.800000, 0.805348,
	 0.809626, 0.814973, 0.823529, 0.831016};
	
      G4double efficiency[NUM] = 
	{0.0, 0.0, 0.0, 0.0,
	 0.0, 0.0, 0.0, 0.0,
	 0.0, 0.0, 0.0, 0.0,
	 0.0, 0.0, 0.0, 0.0,
	 0.0, 0.0, 0.0, 0.0,
	 0.0, 0.0, 0.0, 0.0,
	 0.0, 0.0, 0.0, 0.0,
	 0.0, 0.0, 0.0, 0.0};

      G4MaterialPropertiesTable* SMPT = new G4MaterialPropertiesTable();

      SMPT -> AddProperty("RINDEX",pp,rindex,NUM);
      //SMPT -> AddProperty("SPECULARLOBECONSTANT",pp,specularlobe,NUM);
      //SMPT -> AddProperty("SPECULARSPIKECONSTANT",pp,specularspike,NUM);
      //SMPT -> AddProperty("BACKSCATTERCONSTANT",pp,backscatter,NUM);
      SMPT -> AddProperty("REFLECTIVITY",pp,reflectivity,NUM);
      SMPT -> AddProperty("EFFICIENCY",pp,efficiency,NUM);

      //OpticalMirrorSurface->SetMaterialPropertiesTable(SMPT);

      OpticalMirrorSurface->SetType(dielectric_metal);
      OpticalMirrorSurface->SetModel(glisur);
      OpticalMirrorSurface->SetFinish(polished);

      OpticalMirrorSurface->DumpInfo();
    } // end test

  // GEM stack photocathode (CsI coating of top GEM)
  //
  OpticalPhotocathodeSurface = new G4OpticalSurface("OpticalPhotocathodeSurface", glisur, polished, dielectric_metal);

  const G4int n = 32;
  G4double* photocath_EPHOTON = new G4double[n];  
  for(i = 0; i < n; i++)
    {
      photocath_EPHOTON[i] = cf4_PhotonEnergy[i];
    }
  G4double photocath_EFF[n] =
    {0.654219, 0.630591, 0.609916, 0.587764, 
     0.567089, 0.546414, 0.524262, 0.503586, 
     0.482911, 0.462236, 0.438608, 0.420886, 
     0.398734, 0.376582, 0.355907, 0.335232, 
     0.314557, 0.293882, 0.271730, 0.251055, 
     0.230380, 0.208228, 0.187553, 0.168354, 
     0.144726, 0.125527, 0.103376, 0.082700, 
     0.062025, 0.039873, 0.017721, 0.000000};
  G4double photocath_REFL[n] = 
    {0.,0.,0.,0.,
     0.,0.,0.,0.,
     0.,0.,0.,0.,
     0.,0.,0.,0.,
     0.,0.,0.,0.,
     0.,0.,0.,0.,
     0.,0.,0.,0.,
     0.,0.,0.,0.};
  G4MaterialPropertiesTable* OpticalPhotocathodeSurface_MPT = new G4MaterialPropertiesTable();
  OpticalPhotocathodeSurface_MPT->AddProperty("EFFICIENCY",photocath_EPHOTON,photocath_EFF,32);
  OpticalPhotocathodeSurface_MPT->AddProperty("REFLECTIVITY",photocath_EPHOTON,photocath_REFL,32);
  OpticalPhotocathodeSurface->SetMaterialPropertiesTable(OpticalPhotocathodeSurface_MPT);

  // done
  //
  return;

}


EicRichGemTbMaterial::~EicRichGemTbMaterial(){ }
