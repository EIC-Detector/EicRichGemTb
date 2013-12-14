#ifndef EicRichGemTbStackingAction_H
#define EicRichGemTbStackingAction_H 1

#include "globals.hh"
#include "G4UserStackingAction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/** @TODO Review this class completely
 */
class EicRichGemTbStackingAction : public G4UserStackingAction
{
public:
  EicRichGemTbStackingAction();
  ~EicRichGemTbStackingAction();

public:
  G4ClassificationOfNewTrack ClassifyNewTrack(const G4Track* aTrack);
  void NewStage();
  void PrepareNewEvent();

private:
  G4int gammaCounter;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

