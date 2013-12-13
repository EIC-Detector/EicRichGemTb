#ifndef EicRichGemTbPhysicsListMessenger_h
#define EicRichGemTbPhysicsListMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class EicRichGemTbPhysicsList;
class G4UIdirectory;
class G4UIcmdWithAnInteger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class EicRichGemTbPhysicsListMessenger: public G4UImessenger
{
public:
  EicRichGemTbPhysicsListMessenger(EicRichGemTbPhysicsList* );
  ~EicRichGemTbPhysicsListMessenger();

  void SetNewValue(G4UIcommand*, G4String);

private:
  EicRichGemTbPhysicsList*     pPhysicsList;

  G4UIdirectory*        EicRichGemTbDir;
  G4UIdirectory*        physDir;
  G4UIcmdWithAnInteger* verboseCmd;
  G4UIcmdWithAnInteger* cerenkovCmd;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

