#include "G4RunManager.hh"
#include "G4UImanager.hh"

#include "EicRichGemTbPhysicsList.hh"
#include "EicRichGemTbPrimaryGeneratorAction.hh"
#include "EicRichGemTbDetectorConstruction.hh"
#include "EicRichGemTbRunAction.hh"
#include "EicRichGemTbStackingAction.hh"
#include "EicRichGemTbSteppingVerbose.hh"
#include "EicRichGemTbSteppingAction.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

using namespace std;

int main(int argc,char** argv)
{
  // Seed the random number generator manually
  //
  G4long myseed = 345354;
  CLHEP::HepRandom::setTheSeed(myseed);

  // User Verbose output class
  //
  G4VSteppingVerbose* verbosity = new EicRichGemTbSteppingVerbose;
  G4VSteppingVerbose::SetInstance(verbosity);

  // Run manager
  //
  G4RunManager* runManager = new G4RunManager;

  // UserInitialization classes - mandatory
  //
  G4VUserPhysicsList* physics = new EicRichGemTbPhysicsList;
  runManager-> SetUserInitialization(physics);
  //
  //G4VUserDetectorConstruction* detector = new EicRichGemTbDetectorConstruction;
  G4VUserDetectorConstruction* detector = new EicRichGemTbDetectorConstruction;
  runManager-> SetUserInitialization(detector);
  //
  G4VUserPrimaryGeneratorAction* gen_action = new EicRichGemTbPrimaryGeneratorAction;
  runManager->SetUserAction(gen_action);
  //
  runManager->SetUserAction(new EicRichGemTbSteppingAction());

#ifdef G4VIS_USE
  // visualization manager
  //
  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();
#endif

  // UserAction classes
  //
  G4UserRunAction* run_action = new EicRichGemTbRunAction;
  runManager->SetUserAction(run_action);
  //
  G4UserStackingAction* stacking_action = new EicRichGemTbStackingAction;
  runManager->SetUserAction(stacking_action);

  // Initialize G4 kernel
  //
  runManager->Initialize();

  // Get the pointer to the User Interface manager
  //
  G4UImanager* UImanager = G4UImanager::GetUIpointer();

  if (argc==1)   // Define UI session for interactive mode
    {
#ifdef G4UI_USE
      G4UIExecutive * ui = new G4UIExecutive(argc,argv);
#ifdef G4VIS_USE
      UImanager->ApplyCommand("/control/execute vis.mac");
#endif
      ui->SessionStart();
      delete ui;
#endif
    }
  else         // Batch mode
    {
      G4String command = "/control/execute ";
      G4String fileName = argv[1];
      UImanager->ApplyCommand(command+fileName);
    }

  // Job termination
  // Free the store: user actions, physics_list and detector_description are
  //                 owned and deleted by the run manager, so they should not
  //                 be deleted in the main() program !

#ifdef G4VIS_USE
  delete visManager;
#endif
  delete runManager;
  delete verbosity;

  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
