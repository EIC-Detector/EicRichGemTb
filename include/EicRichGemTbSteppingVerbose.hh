class EicRichGemTbSteppingVerbose;

#ifndef EicRichGemTbSteppingVerbose_h
#define EicRichGemTbSteppingVerbose_h 1

#include "G4SteppingVerbose.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/** @TODO Review this class completely
 */
class EicRichGemTbSteppingVerbose : public G4SteppingVerbose
{
public:

  EicRichGemTbSteppingVerbose();
  ~EicRichGemTbSteppingVerbose();

  void StepInfo();
  void TrackingStarted();

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
