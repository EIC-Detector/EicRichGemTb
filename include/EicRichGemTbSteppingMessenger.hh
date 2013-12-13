#ifndef EicRichGemTbSteppingMessenger_h
#define EicRichGemTbSteppingMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

class EicRichGemTbSteppingAction;
class G4UIcmdWithABool;

class EicRichGemTbSteppingMessenger: public G4UImessenger
{
public:
  EicRichGemTbSteppingMessenger(EicRichGemTbSteppingAction*);
  virtual ~EicRichGemTbSteppingMessenger();

  virtual void SetNewValue(G4UIcommand*, G4String);

private:

  EicRichGemTbSteppingAction*        fStepping;
  G4UIcmdWithABool*  fOneStepPrimariesCmd;
};

#endif
