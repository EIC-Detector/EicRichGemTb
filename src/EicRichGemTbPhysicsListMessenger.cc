
#include "EicRichGemTbPhysicsListMessenger.hh"

#include "EicRichGemTbPhysicsList.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAnInteger.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EicRichGemTbPhysicsListMessenger::EicRichGemTbPhysicsListMessenger(EicRichGemTbPhysicsList* pPhys)
 :pPhysicsList(pPhys)
{
  EicRichGemTbDir = new G4UIdirectory("/EicRichGemTb/");
  EicRichGemTbDir->SetGuidance("UI commands of this example");

  physDir = new G4UIdirectory("/EicRichGemTb/phys/");
  physDir->SetGuidance("PhysicsList control");

  verboseCmd = new G4UIcmdWithAnInteger("/EicRichGemTb/phys/verbose",this);
  verboseCmd->SetGuidance("set verbose for physics processes");
  verboseCmd->SetParameterName("verbose",true);
  verboseCmd->SetDefaultValue(1);
  verboseCmd->SetRange("verbose>=0");
  verboseCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  cerenkovCmd = new G4UIcmdWithAnInteger("/EicRichGemTb/phys/cerenkovMaxPhotons",this);
  cerenkovCmd->SetGuidance("set max nb of photons per step");
  cerenkovCmd->SetParameterName("MaxNumber",false);
  cerenkovCmd->SetRange("MaxNumber>=0");
  cerenkovCmd->AvailableForStates(G4State_Idle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EicRichGemTbPhysicsListMessenger::~EicRichGemTbPhysicsListMessenger()
{
  delete verboseCmd;
  delete cerenkovCmd;
  delete physDir;
  delete EicRichGemTbDir;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EicRichGemTbPhysicsListMessenger::SetNewValue(G4UIcommand* command,
                                             G4String newValue)
{
//  if( command == verboseCmd )
//    { pPhysicsList->SetVerbose(verboseCmd->GetNewIntValue(newValue));}
//
//  if( command == cerenkovCmd )
//    {pPhysicsList->SetNbOfPhotonsCerenkov(cerenkovCmd->GetNewIntValue(newValue));}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
