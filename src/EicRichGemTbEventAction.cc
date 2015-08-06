#include "EicRichGemTbEventAction.hh"

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"

#include "EicRichGemTbAnalysis.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EicRichGemTbEventAction::EicRichGemTbEventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EicRichGemTbEventAction::~EicRichGemTbEventAction()
{}


EicRichGemTbPhotoHitsCollection* EicRichGemTbEventAction::GetPhotoHitsCollection(const G4String& hcName, const G4Event* event)
{
  G4int hcID = G4SDManager::GetSDMpointer()->GetCollectionID(hcName);

  EicRichGemTbPhotoHitsCollection* hitsCollection
    = static_cast<EicRichGemTbPhotoHitsCollection*>(
                                           event->GetHCofThisEvent()->GetHC(hcID));

  if ( ! hitsCollection ) {
    G4cerr << "Cannot access hitsCollection " << hcName << G4endl;
    exit(1);
  }

  return hitsCollection;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EicRichGemTbEventAction::BeginOfEventAction(const G4Event*)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EicRichGemTbEventAction::EndOfEventAction(const G4Event* evt)
{
  // Get hits collections
  //
  EicRichGemTbPhotoHitsCollection* photoHits = GetPhotoHitsCollection("PhotoHitsCollection", evt);

  //photoHits->PrintAllHits();

  // Print per event (modulo n)
  //
  G4int eventID = evt->GetEventID();

  // Fill histograms, ntuple
  //

  // get analysis manager
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

  // print out total number of photons registered on CsI readout plane
  G4cout << "Number of optical photons registered on CsI in this event: " << photoHits->GetSize() << G4endl;

  // fill ntuple
  for ( G4int ihit = 0; ihit < photoHits->GetSize(); ihit++ )
    {
      EicRichGemTbPhotoHit* photoHit_i = static_cast<EicRichGemTbPhotoHit*>( photoHits->GetHit(ihit) );
      analysisManager->FillNtupleDColumn(0, eventID);
      analysisManager->FillNtupleDColumn(1, ihit+1 );
      analysisManager->FillNtupleDColumn(2, photoHit_i->GetHitPos().getX() );
      analysisManager->FillNtupleDColumn(3, photoHit_i->GetHitPos().getY() );
      analysisManager->FillNtupleDColumn(4, photoHit_i->GetHitPos().getZ() );
      analysisManager->AddNtupleRow();
    }

  // periodic printing
  //
  if (eventID < 100 || eventID%100 == 0) {
    G4cout << ">>> Event " << evt->GetEventID() << " finished." << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
