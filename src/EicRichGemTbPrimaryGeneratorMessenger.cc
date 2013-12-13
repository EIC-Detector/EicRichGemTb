#include "EicRichGemTbPrimaryGeneratorMessenger.hh"

#include "EicRichGemTbPrimaryGeneratorAction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EicRichGemTbPrimaryGeneratorMessenger::EicRichGemTbPrimaryGeneratorMessenger(
                                                                 EicRichGemTbPrimaryGeneratorAction* EicRichGemTbGun)
 :EicRichGemTbAction(EicRichGemTbGun)
{
  gunDir = new G4UIdirectory("/EicRichGemTb/gun/");
  gunDir->SetGuidance("PrimaryGenerator control");

  polarCmd = new G4UIcmdWithADoubleAndUnit("/EicRichGemTb/gun/optPhotonPolar",this);
  polarCmd->SetGuidance("Set linear polarization");
  polarCmd->SetGuidance("  angle w.r.t. (k,n) plane");
  polarCmd->SetParameterName("angle",true);
  polarCmd->SetUnitCategory("Angle");
  polarCmd->SetDefaultValue(-360.0);
  polarCmd->SetDefaultUnit("deg");
  polarCmd->AvailableForStates(G4State_Idle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EicRichGemTbPrimaryGeneratorMessenger::~EicRichGemTbPrimaryGeneratorMessenger()
{
  delete polarCmd;
  delete gunDir;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EicRichGemTbPrimaryGeneratorMessenger::SetNewValue(
                                                  G4UIcommand* command, G4String newValue)
{
  if( command == polarCmd ) {
    G4double angle = polarCmd->GetNewDoubleValue(newValue);
    if ( angle == -360.0*deg ) {
      EicRichGemTbAction->SetOptPhotonPolar();
    } else {
      EicRichGemTbAction->SetOptPhotonPolar(angle);
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

