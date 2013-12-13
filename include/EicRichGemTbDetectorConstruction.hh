/*!
 * \file ${file_name}
 * \brief
 * \author Nils Feege <nils.feege@stonybrook.edu>
 * \version $$Revision: 1.0 $$
 * \date $$Date: 2013/11/15 18:11:08 $$
 */

#ifndef EicRichGemTbDETECTORCONSTRUCTION_H_
#define EicRichGemTbDETECTORCONSTRUCTION_H_ 1

#include <cmath>
#include <map>
#include <utility>

#include "G4PhysicalConstants.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4String.hh"

#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4Element.hh"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"

class G4VPhysicalVolume;
class G4LogicalVolume;
class G4PVPlacement;

//! RICH geometry data
//! use GEANT units! Use
class EicRichGemTb_Geometry
{
public:
  EicRichGemTb_Geometry()
  {
    SetDefault();
  }
  void
  SetDefault();

  //Unit
  static double Unit_cm() {return cm;}

  // derivative constants
  double
  get_R_max() const
  {
    return R_mirror_ref + dR_mirror + dR_mirror_spt + dR_backwindow;
  }

  double
  get_cone_size_z() const
  {
    return (z_shift + get_R_max()) * 2;
  }

  double
  get_R_frontwindow() const;

  double
  get_half_angle_HBD() const;

  double
  get_RZ_Seg1_HBD() const;

  double
  get_RZ_Seg2_HBD() const;

  double
  get_R_Tip_HBD() const;

  double
  get_Z_Tip_HBD() const;

  double
  get_Rotation_HBD() const;

private:

  int N_RICH_Sector;

  double min_eta;
  double R_beam_pipe;

  double z_shift;
  double R_shift;
  double frontwindow_DisplaceRatio; // Displace R,Z and radius simultainously
  double dR_frontwindow_shrink;
  double R_mirror_ref;

  double dR_mirror;
  double dR_mirror_spt;
  double dR_backwindow;
  double dR_frontwindow;

  int n_GEM_layers;
  double HBD_thickness;

  G4String RICH_gas_mat;
  G4String RICH_Mirror_mat;
  G4String RICH_Gas_Window_mat;

public:
  double
  get_R_beam_pipe() const
  {
    return R_beam_pipe;
  }

  double
  get_frontwindow_DisplaceRatio() const
  {
    return frontwindow_DisplaceRatio;
  }

  double
  get_min_eta() const
  {
    return min_eta;
  }

  double
  get_R_mirror_ref() const
  {
    return R_mirror_ref;
  }

  double
  get_dR_backwindow() const
  {
    return dR_backwindow;
  }

  double
  get_dR_frontwindow() const
  {
    return dR_frontwindow;
  }

  double
  get_dR_frontwindow_shrink() const
  {
    return dR_frontwindow_shrink;
  }

  double
  get_dR_mirror() const
  {
    return dR_mirror;
  }

  double
  get_dR_mirror_spt() const
  {
    return dR_mirror_spt;
  }

  G4String
  get_RICH_gas_mat() const
  {
    return RICH_gas_mat;
  }

  G4String
  get_RICH_Gas_Window_mat() const
  {
    return RICH_Gas_Window_mat;
  }

  G4String
  get_RICH_Mirror_mat() const
  {
    return RICH_Mirror_mat;
  }

  int
  get_N_RICH_Sector() const
  {
    return N_RICH_Sector;
  }

  double
  get_z_shift() const
  {
    return z_shift;
  }

  double
  get_R_shift() const
  {
    return R_shift;
  }

  void
  set_R_beam_pipe(double beamPipe)
  {
    R_beam_pipe = beamPipe;
  }

  void
  set_frontwindow_DisplaceRatio(double frontwindowDisplaceRatio)
  {
    frontwindow_DisplaceRatio = frontwindowDisplaceRatio;
  }

  void
  set_min_eta(double minEta)
  {
    min_eta = minEta;
  }

  void
  set_R_mirror_ref(double mirrorRef)
  {
    R_mirror_ref = mirrorRef;
  }

  void
  set_dR_backwindow(double rBackwindow)
  {
    dR_backwindow = rBackwindow;
  }

  void
  set_dR_frontwindow(double rFrontwindow)
  {
    dR_frontwindow = rFrontwindow;
  }

  void
  set_dR_frontwindow_shrink(double rFrontwindowShrink)
  {
    dR_frontwindow_shrink = rFrontwindowShrink;
  }

  void
  set_dR_mirror(double rMirror)
  {
    dR_mirror = rMirror;
  }

  void
  set_dR_mirror_spt(double rMirrorSpt)
  {
    dR_mirror_spt = rMirrorSpt;
  }

  void
  set_RICH_gas_mat(G4String richGasMat)
  {
    RICH_gas_mat = richGasMat;
  }

  void
  set_RICH_Gas_Window_mat(G4String richGasWindowMat)
  {
    RICH_Gas_Window_mat = richGasWindowMat;
  }

  void
  set_RICH_Mirror_mat(G4String richMirrorMat)
  {
    RICH_Mirror_mat = richMirrorMat;
  }

  void
  set_N_RICH_Sector(int richSector)
  {
    N_RICH_Sector = richSector;
  }

  void
  set_z_shift(double shift)
  {
    z_shift = shift;
  }

  void
  set_R_shift(double shift)
  {
    R_shift = shift;
  }

  int
  get_n_GEM_layers() const
  {
    return n_GEM_layers;
  }

  double
  get_HBD_thickness() const
  {
    return HBD_thickness;
  }

  void
  set_n_GEM_layers(int gemLayers)
  {
    n_GEM_layers = gemLayers;
  }

  void
  set_HBD_thickness(double hbdThickness)
  {
    HBD_thickness = hbdThickness;
  }
};

class EicRichGemTbDetectorConstruction : public G4VUserDetectorConstruction
{
public:
  EicRichGemTbDetectorConstruction();
  ~EicRichGemTbDetectorConstruction();

  G4VPhysicalVolume* Construct();

  G4VPhysicalVolume* Construct2();

  G4LogicalVolume*
  Construct_HBD(G4LogicalVolume* RICHSecLog);

  G4LogicalVolume*
  Construct_HBD_Layers(G4LogicalVolume* RICHHBDLog, const G4String name,
                       const G4String material, const double start_z, const double thickness);

  EicRichGemTb_Geometry geom;

protected:

  G4LogicalVolume *
  RegisterLogicalVolume(G4LogicalVolume*);

  typedef std::map<G4String, G4LogicalVolume*> map_log_vol_t;
  map_log_vol_t map_log_vol;

  G4PVPlacement *
  RegisterPhysicalVolume(G4PVPlacement*);

  typedef std::pair<G4String, G4int> phy_vol_idx_t;
  typedef std::map<phy_vol_idx_t, G4PVPlacement*> map_phy_vol_t;
  map_phy_vol_t map_phy_vol;

private:
  void DefineMaterials();
  G4VPhysicalVolume* ConstructDetector();

  G4Box* fExperimentalHall_box;
  G4LogicalVolume* fExperimentalHall_log;
  G4VPhysicalVolume* fExperimentalHall_phys;

  // Materials & Elements
  G4Element* fH;
  G4Element* fN;
  G4Element* fO;
  G4Material* fAir;

  G4Element* fC;
  G4Element* fF;
  G4Material* fCF4;

  G4Element* fSi;

  G4Material* fAl;
  G4Material* fSiO2MirrorQuartz;

  G4MaterialPropertiesTable* fMPT_cf4;

  // Geometry
  G4double expHall_x;
  G4double expHall_y;
  G4double expHall_z;

  G4double tank_r;
  G4double tank_dr;
  G4double tank_dz;

  G4double mirror_f;
  G4double mirror_r;
  G4double mirror_cr;
  G4double mirror_dz;

};

#endif /* EicRichGemTbDETECTORCONSTRUCTION_H_ */
