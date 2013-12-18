#include "EicRichGemTbPhotoHit.hh"
#include "G4ios.hh"
#include "G4VVisManager.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"

G4Allocator<EicRichGemTbPhotoHit> EicRichGemTbPhotoHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//EicRichGemTbPhotoHit::EicRichGemTbPhotoHit()
//  : fPmtNumber(-1),fPhotons(0),fPhysVol(0),fDrawit(false) {}

EicRichGemTbPhotoHit::EicRichGemTbPhotoHit() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EicRichGemTbPhotoHit::~EicRichGemTbPhotoHit() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EicRichGemTbPhotoHit::EicRichGemTbPhotoHit(const EicRichGemTbPhotoHit &right) : G4VHit()
{
  //  fPmtNumber=right.fPmtNumber;
  //  fPhotons=right.fPhotons;
  //  fPhysVol=right.fPhysVol;
  //  fDrawit=right.fDrawit;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const EicRichGemTbPhotoHit& EicRichGemTbPhotoHit::operator=(const EicRichGemTbPhotoHit &right){
  //  fPmtNumber = right.fPmtNumber;
  //  fPhotons=right.fPhotons;
  //  fPhysVol=right.fPhysVol;
  //  fDrawit=right.fDrawit;
  return *this;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int EicRichGemTbPhotoHit::operator==(const EicRichGemTbPhotoHit &right) const{
  //return (fPmtNumber==right.fPmtNumber);
  return (0);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//void EicRichGemTbPhotoHit::Draw(){
//  if(fDrawit&&fPhysVol){ //ReDraw only the PMTs that have hit counts > 0
//    //Also need a physical volume to be able to draw anything
//    G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
//    if(pVVisManager){//Make sure that the VisManager exists
//      G4VisAttributes attribs(G4Colour(1.,0.,0.));
//      attribs.SetForceSolid(true);
//      G4RotationMatrix rot;
//      if(fPhysVol->GetRotation())//If a rotation is defined use it
//        rot=*(fPhysVol->GetRotation());
//      G4Transform3D trans(rot,fPhysVol->GetTranslation());//Create transform
//      pVVisManager->Draw(*fPhysVol,attribs,trans);//Draw it
//    }
//  }
//}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EicRichGemTbPhotoHit::Print()
{
  G4cout << "HIT at " << fHitPos.getX() << " " << fHitPos.getY() << " " << fHitPos.getZ() << G4endl;
}
