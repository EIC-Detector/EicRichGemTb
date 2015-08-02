#ifndef EicRichGemTbMaterial_h
#define EicRichGemTbMaterial_h 1
#include "globals.hh"
#include <vector>
#include "G4String.hh"
#include "G4Material.hh"
//#include "G4MaterialTable.hh"
#include "G4Element.hh"
#include "G4ElementTable.hh"
#include "G4OpticalSurface.hh"

class EicRichGemTbMaterial {

public:
  EicRichGemTbMaterial();
  virtual ~EicRichGemTbMaterial();

  G4Material* getAir() {return AmbientAir;}
  G4Material* getMirrorQuartz() {return SiO2MirrorQuartz;}

  G4Material* getVacuum() {return Vacuum;}

  G4Material* getCF4() {return CF4;}
  G4Material* getAluminum() {return Aluminum;}

  G4Material* getStainlessSteel() {return StainlessSteel;}
  G4Material* getCsI() {return CsI;}

  G4OpticalSurface* getOpticalMirrorSurface() {return OpticalMirrorSurface;}
  G4OpticalSurface* getOpticalPhotocathodeSurface() {return OpticalPhotocathodeSurface;}

private:
  G4Element* H;
  G4Element* N;
  G4Element* O;
  G4Element* C;
  G4Element* F;
  G4Element* Mg;
  G4Element* Si;
  G4Element* Ar;
  G4Element* I;
  G4Element* Cs;


  G4Material* Vacuum;
  G4Material* CO2;
  G4Material* AmbientAir;
  G4Material* CF4;
  G4Material* CsI;
  G4Material* StainlessSteel;
  G4Material* Aluminum;
  G4Material* MgF2;
  G4Material* SiO2MirrorQuartz;

  G4OpticalSurface* OpticalMirrorSurface;
  G4OpticalSurface* OpticalPhotocathodeSurface;

};

#endif

