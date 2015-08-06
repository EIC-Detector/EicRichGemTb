// Make this appear first!
#include "G4Timer.hh"

#include "EicRichGemTbRunAction.hh"

#include "G4Run.hh"

#include "EicRichGemTbAnalysis.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EicRichGemTbRunAction::EicRichGemTbRunAction()
{
  timer = new G4Timer;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EicRichGemTbRunAction::~EicRichGemTbRunAction()
{
  delete timer;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EicRichGemTbRunAction::BeginOfRunAction(const G4Run* aRun)
{
  G4cout << "### Run " << aRun->GetRunID() << " starts." << G4endl;
  timer->Start();

  //inform the runManager to save random number seed
  //G4RunManager::GetRunManager()->SetRandomNumberStore(true);

  // Book histograms, ntuple
  //

  // Create analysis manager
  // The choice of analysis technology is done via selectin of a namespace
  // in B4Analysis.hh
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  G4cout << "Using " << analysisManager->GetType()
         << " analysis manager" << G4endl;

  // Create directories
  //analysisManager->SetHistoDirectoryName("histograms");
  //analysisManager->SetNtupleDirectoryName("ntuple");

  // Open an output file
  //
  G4String fileName = "simData_EicRichGemTb";
  analysisManager->OpenFile(fileName);

  // Creating histograms
  //
  analysisManager->SetFirstHistoId(1);
  analysisManager->CreateH1("1","Number of photons generated", 301, -0.5, 300.5);
  analysisManager->CreateH1("2","Number of photons recorded", 301, -0.5, 300.5);

  // Creating tree
  //
  analysisManager->CreateNtuple("photoHits", "photon hits on GEM stack");
  analysisManager->CreateNtupleDColumn("event");
  analysisManager->CreateNtupleDColumn("nhit");
  analysisManager->CreateNtupleDColumn("x");
  analysisManager->CreateNtupleDColumn("y");
  analysisManager->CreateNtupleDColumn("z");
  analysisManager->FinishNtuple();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EicRichGemTbRunAction::EndOfRunAction(const G4Run* aRun)
{
  timer->Stop();
  G4int nofEvents = aRun->GetNumberOfEvent();
  G4cout << "Event number = " << nofEvents << " " << *timer << G4endl;

  //  if ( nofEvents == 0 ) return;

  // save analysis output
  //
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

  analysisManager->Write();
  analysisManager->CloseFile();

  // complete cleanup
  //
  delete G4AnalysisManager::Instance();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
