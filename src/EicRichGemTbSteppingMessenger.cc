#include "EicRichGemTbSteppingMessenger.hh"
#include "EicRichGemTbSteppingAction.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithABool.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EicRichGemTbSteppingMessenger::EicRichGemTbSteppingMessenger(EicRichGemTbSteppingAction* step)
  : fStepping(step)
{
  fOneStepPrimariesCmd = new G4UIcmdWithABool("/EicRichGemTb/oneStepPrimaries",this);
  fOneStepPrimariesCmd->SetGuidance("Only allows primaries to go one step in the gas volume before being killed.");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EicRichGemTbSteppingMessenger::~EicRichGemTbSteppingMessenger(){
  delete fOneStepPrimariesCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void
EicRichGemTbSteppingMessenger::SetNewValue(G4UIcommand* command,G4String newValue){
  if( command == fOneStepPrimariesCmd ){
    fStepping->SetOneStepPrimaries(fOneStepPrimariesCmd->GetNewBoolValue(newValue));
  }
}
