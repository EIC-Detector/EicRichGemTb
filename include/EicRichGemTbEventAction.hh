#ifndef EicRichGemTbEventAction_h
#define EicRichGemTbEventAction_h 1

#include "G4UserEventAction.hh"

#include "EicRichGemTbPhotoHit.hh"

class G4Event;

/// Event action

class EicRichGemTbEventAction : public G4UserEventAction
{
public:
  EicRichGemTbEventAction();
  ~EicRichGemTbEventAction();

public:
  void BeginOfEventAction(const G4Event*);
  void EndOfEventAction(const G4Event*);

private:
  EicRichGemTbPhotoHitsCollection* GetPhotoHitsCollection(const G4String& hcName, const G4Event* event);

};

#endif


