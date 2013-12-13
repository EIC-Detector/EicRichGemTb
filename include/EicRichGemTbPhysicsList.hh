
#ifndef EicRichGemTbPhysicsList_h
#define EicRichGemTbPhysicsList_h 1

#include "globals.hh"
#include "G4VUserPhysicsList.hh"

class G4Cerenkov;
class G4Scintillation;
class G4OpAbsorption;
class G4OpRayleigh;
class G4OpMieHG;
class G4OpBoundaryProcess;

class EicRichGemTbPhysicsListMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class EicRichGemTbPhysicsList : public G4VUserPhysicsList
{
public:
  EicRichGemTbPhysicsList();
  ~EicRichGemTbPhysicsList();

public:
  void ConstructParticle();
  void ConstructProcess();

  void SetCuts();

  //these methods Construct particles
  void ConstructBosons();
  void ConstructLeptons();
  void ConstructMesons();
  void ConstructBaryons();

  //these methods Construct physics processes and register them
  void ConstructGeneral();
  void ConstructEM();
  void ConstructOp();

  //for the Messenger
  void SetVerbose(G4int);
  void SetNbOfPhotonsCerenkov(G4int);

private:
  G4Cerenkov*          theCerenkovProcess;
  G4Scintillation*     theScintillationProcess;
  G4OpAbsorption*      theAbsorptionProcess;
  G4OpRayleigh*        theRayleighScatteringProcess;
  G4OpMieHG*           theMieHGScatteringProcess;
  G4OpBoundaryProcess* theBoundaryProcess;

  EicRichGemTbPhysicsListMessenger* pMessenger;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif /* EicRichGemTbPhysicsList_h */
