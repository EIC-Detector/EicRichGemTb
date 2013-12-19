#include "EicRichGemTbSD.hh"
#include "G4Step.hh"
#include "G4HCofThisEvent.hh"
#include "G4Track.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "G4VProcess.hh"
#include "G4OpticalPhoton.hh"

EicRichGemTbSD::EicRichGemTbSD(G4String DetName, G4String HitColName="EicRichGemTbPhotoHitCollection" )
  :G4VSensitiveDetector(DetName),HCID(-1)
{
  G4String HCname;
  collectionName.insert(HitColName);
}

EicRichGemTbSD::~EicRichGemTbSD(){

}

void EicRichGemTbSD::Initialize(G4HCofThisEvent*)
{

  PhotoHitCollection = new EicRichGemTbPhotoHitsCollection
    (SensitiveDetectorName,collectionName[0]);

}


G4bool EicRichGemTbSD::ProcessHits(G4Step*aStep,G4TouchableHistory*ROhist)
{

  //  G4Track* aTrack = aStep->GetTrack();
  //  //if(aTrack ->GetParticleDefinition()->GetPDGCharge() == 0.0 ) return false;
  //  G4double edep = aStep->GetTotalEnergyDeposit();
  //
  //  G4cout << "SD HIT by " << aTrack->GetParticleDefinition()->GetParticleType() << " - " << aTrack->GetParticleDefinition()->GetParticleSubType()  <<  G4endl;
  //
  //  if(edep < 0.0001) return false;
  //  G4StepPoint* PrePosition =  aStep->GetPreStepPoint();
  //
  //
  //  G4int CurrentHpdNumber =  PrePosition->GetTouchableHandle()
  //    -> GetReplicaNumber(1);
  //  G4int CurrentSectorNumber= PrePosition->GetTouchableHandle()
  //    -> GetReplicaNumber();
  //  G4VPhysicalVolume* ROphysVol = ROhist -> GetVolume();
  //  G4int CurrentPixelNumber =  ROphysVol->GetCopyNo();
  //
  //  G4ThreeVector HitAtPhotoCathode;
  //  if(aTrack->GetCreatorProcess() -> GetProcessName() == "PadHpdPhot") {
  //    HitAtPhotoCathode=aTrack->GetVertexPosition();
  //  }
  //
  //  G4int CopyId= NumberOfSensitiveSectorsInHpd*
  //    NumberOfSensitivePixelsInSector*CurrentHpdNumber+
  //    NumberOfSensitivePixelsInSector*CurrentSectorNumber+
  //    CurrentPixelNumber;
  //
  //
  //  if( HpdSDID[CopyId ] == -1 ) {
  //    EicRichGemTbHit* newHit = new EicRichGemTbHit();
  //    newHit->SetEdep( edep );
  //    newHit->SetPos( aStep->GetPreStepPoint()->GetPosition() );
  //    newHit->SetPosPC(HitAtPhotoCathode);
  //    newHit->SetCurHpdNum ( CurrentHpdNumber );
  //    newHit->SetCurSectNum ( CurrentSectorNumber );
  //    newHit->SetCurPixNum (  CurrentPixelNumber );
  //    //For now the same SD detector cell is allowed to have
  //    //multiple hits from the same event. This is ok for
  //    // analogue readout for photoelectron counting.
  //    //This is to be changed in the future.
  //    G4int NumHits = EicRichGemTbHitCollection->insert( newHit );
  //
  //    HpdSDID[CopyId]= NumHits -1 ;
  //
  //  } else {
  //    // the current pixel is already hit in this event.
  //    // here we can add extra energy (adc counts) to the
  //    // existing hit. But this is not relevant to the
  //    // current Tb analysis.
  //    (*EicRichGemTbHitCollection)[HpdSDID[CopyId]]->AddEdep( edep );
  //
  //    G4cout << " Multiple hits in Hpd sector pixel " << CurrentHpdNumber
  //           <<"   "<< CurrentSectorNumber<<"   "<< CurrentPixelNumber<< G4endl;
  //  }

  return true;
}

G4bool EicRichGemTbSD::ProcessHits_constStep(const G4Step* aStep, G4TouchableHistory* ROhist)
{
  //need to know if this is an optical photon
  if(aStep->GetTrack()->GetDefinition() != G4OpticalPhoton::OpticalPhotonDefinition())
    return false;

  //  G4cout << "HIT!" << G4endl;

  G4StepPoint* thePrePoint = aStep->GetPreStepPoint();
  G4TouchableHistory* theTouchable = (G4TouchableHistory*)(aStep->GetPreStepPoint()->GetTouchable());
  G4VPhysicalVolume* thePrePV = theTouchable->GetVolume();

  G4StepPoint* thePostPoint = aStep->GetPostStepPoint();

  //User replica number 1 since photocathode is a daughter volume
  //to the pmt which was replicated
  //G4int pmtNumber=
  //  aStep->GetPostStepPoint()->GetTouchable()->GetReplicaNumber(1);
  //G4VPhysicalVolume* physVol=
  //  aStep->GetPostStepPoint()->GetTouchable()->GetVolume(1);

  //Find the correct hit collection
  //G4int n=fPMTHitCollection->entries();
  EicRichGemTbPhotoHit* hit=NULL;
  //for(G4int i=0;i<n;i++){
  //  if((*fPMTHitCollection)[i]->GetPMTNumber()==pmtNumber){
  //    hit=(*fPMTHitCollection)[i];
  //    break;
  //  }
  //}

  //Get the average position of the hit
  //  G4ThreeVector hitpos = thePrePoint->GetPosition() + thePostPoint->GetPosition();
  //  hitpos/=2.;

  G4ThreeVector hitpos = thePostPoint->GetPosition();

  if(hit==NULL){//this pad wasnt previously hit in this event
    hit = new EicRichGemTbPhotoHit(); //so create new hit
    //  hit->SetPMTNumber(pmtNumber);
    //  hit->SetPMTPhysVol(physVol);
    //    hit->SetPMTPos((*fPMTPositionsX)[pmtNumber],(*fPMTPositionsY)[pmtNumber],(*fPMTPositionsZ)[pmtNumber]);
    hit->SetHitPos(hitpos.getX(), hitpos.getY(), hitpos.getZ());
    PhotoHitCollection->insert(hit);
  }

  //hit->IncPhotonCount(); //increment hit for the selected pmt

  //if(!LXeDetectorConstruction::GetSphereOn()){
  //  hit->SetDrawit(true);
  //  //If the sphere is disabled then this hit is automaticaly drawn
  //}
  //else{//sphere enabled
  //  LXeUserTrackInformation* trackInfo=
  //    (LXeUserTrackInformation*)aStep->GetTrack()->GetUserInformation();
  //  if(trackInfo->GetTrackStatus()&hitSphere)
  //    //only draw this hit if the photon has hit the sphere first
  //    hit->SetDrawit(true);
  //}

  return true;
}


void  EicRichGemTbSD::EndOfEvent(G4HCofThisEvent*HCE){
  if( HCID < 0 ){
    HCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]); }

  HCE->AddHitsCollection( HCID, PhotoHitCollection );

}

void  EicRichGemTbSD::clear(){}

void   EicRichGemTbSD::DrawAll(){ }

void   EicRichGemTbSD::PrintAll(){ }




