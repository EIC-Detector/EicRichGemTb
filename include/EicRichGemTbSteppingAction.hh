#ifndef EicRichGemTbSteppingAction_H
#define EicRichGemTbSteppingACtion_H 1

#include "globals.hh"
#include "G4UserSteppingAction.hh"

#include "G4OpBoundaryProcess.hh"

//class EicRichGemTbRecorderBase;
//class EicRichGemTbEventAction;
class EicRichGemTbTrackingAction;
class EicRichGemTbSteppingMessenger;

class EicRichGemTbSteppingAction : public G4UserSteppingAction
{
public:

  //  EicRichGemTbSteppingAction(EicRichGemTbRecorderBase*);
  EicRichGemTbSteppingAction();
  virtual ~EicRichGemTbSteppingAction();
  virtual void UserSteppingAction(const G4Step*);

  void SetOneStepPrimaries(G4bool b){fOneStepPrimaries=b;}
  G4bool GetOneStepPrimaries(){return fOneStepPrimaries;}

private:

  //  EicRichGemTbRecorderBase* fRecorder;
  G4bool fOneStepPrimaries;
  EicRichGemTbSteppingMessenger* fSteppingMessenger;

  G4OpBoundaryProcessStatus fExpectedNextStatus;
};

#endif
