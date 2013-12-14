#ifndef EicRichGemTbPrimaryGeneratorMessenger_h
#define EicRichGemTbPrimaryGeneratorMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

class EicRichGemTbPrimaryGeneratorAction;
class G4UIdirectory;
class G4UIcmdWithADoubleAndUnit;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/** @TODO Review this class completely
 */
class EicRichGemTbPrimaryGeneratorMessenger: public G4UImessenger
{
public:
  EicRichGemTbPrimaryGeneratorMessenger(EicRichGemTbPrimaryGeneratorAction*);
  ~EicRichGemTbPrimaryGeneratorMessenger();

  void SetNewValue(G4UIcommand*, G4String);

private:
  EicRichGemTbPrimaryGeneratorAction* EicRichGemTbAction;
  G4UIdirectory*               gunDir;
  G4UIcmdWithADoubleAndUnit*   polarCmd;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

