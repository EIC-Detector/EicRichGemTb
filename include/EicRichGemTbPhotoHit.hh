#ifndef EicRichGemTbPhotoHit_h
#define EicRichGemTbPhotoHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4LogicalVolume.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"
#include "G4VPhysicalVolume.hh"

class G4VTouchable;

class EicRichGemTbPhotoHit : public G4VHit
{
public:

  EicRichGemTbPhotoHit();
  virtual ~EicRichGemTbPhotoHit();
  EicRichGemTbPhotoHit(const EicRichGemTbPhotoHit &right);

  const EicRichGemTbPhotoHit& operator=(const EicRichGemTbPhotoHit &right);
  G4int operator==(const EicRichGemTbPhotoHit &right) const;

  inline void *operator new(size_t);
  inline void operator delete(void *aHit);

  //    virtual void Draw();
  virtual void Print();

  //    inline void SetDrawit(G4bool b){fDrawit=b;}
  //    inline G4bool GetDrawit(){return fDrawit;}
  //
  //    inline void IncPhotonCount(){fPhotons++;}
  //    inline G4int GetPhotonCount(){return fPhotons;}
  //
  //    inline void SetPMTNumber(G4int n) { fPmtNumber = n; }
  //    inline G4int GetPMTNumber() { return fPmtNumber; }
  //
  //    inline void SetPMTPhysVol(G4VPhysicalVolume* physVol){this->fPhysVol=physVol;}
  //    inline G4VPhysicalVolume* GetPMTPhysVol(){return fPhysVol;}

  inline void SetHitPos(G4double x,G4double y,G4double z){
    fHitPos=G4ThreeVector(x,y,z);
  }

  inline G4ThreeVector GetHitPos(){return fHitPos;}

private:

  //    G4int fPmtNumber;
  //    G4int fPhotons;
  G4ThreeVector fHitPos;
  //    G4VPhysicalVolume* fPhysVol;
  //   G4bool fDrawit;

};

typedef G4THitsCollection<EicRichGemTbPhotoHit> EicRichGemTbPhotoHitsCollection;

extern G4Allocator<EicRichGemTbPhotoHit> EicRichGemTbPhotoHitAllocator;

inline void* EicRichGemTbPhotoHit::operator new(size_t){
  void *aHit;
  aHit = (void *) EicRichGemTbPhotoHitAllocator.MallocSingle();
  return aHit;
}

inline void EicRichGemTbPhotoHit::operator delete(void *aHit){
  EicRichGemTbPhotoHitAllocator.FreeSingle((EicRichGemTbPhotoHit*) aHit);
}

#endif
