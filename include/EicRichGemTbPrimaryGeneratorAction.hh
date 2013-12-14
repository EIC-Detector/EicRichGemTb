#ifndef EicRichGemTbPrimaryGeneratorAction_h
#define EicRichGemTbPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"

class G4ParticleGun;
class G4Event;
class EicRichGemTbPrimaryGeneratorMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/** @TODO Review this class completely
 */
class EicRichGemTbPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
  EicRichGemTbPrimaryGeneratorAction();
  ~EicRichGemTbPrimaryGeneratorAction();

public:
  void GeneratePrimaries(G4Event*);

  void SetOptPhotonPolar();
  void SetOptPhotonPolar(G4double);

private:
  G4ParticleGun* particleGun;
  EicRichGemTbPrimaryGeneratorMessenger* gunMessenger;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif /*EicRichGemTbPrimaryGeneratorAction_h*/
