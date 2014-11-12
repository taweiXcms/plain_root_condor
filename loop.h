#include <TTree.h>

#define MAX_XB 8192
#define MAX_MUON 64
#define MAX_TRACK 4096 //default 2048
#define MAX_GEN 4096 //default 2048
#define MAX_BX 128
#define N_TRIGGER_BOOKINGS 5842

class bNtuple
{
 public:
  float mass;
  float pt;
  float eta;
  float y;
  float d0;
  float px;
  float py;
  float vx;
  float vy;
  float d0Err;
  float chi2ndf;
  float dtheta;

  float trk1Pt;
  float trk1Dxy;
  float trk1D0Err;
  float trk1PixelHit;
  float trk1StripHit;
  float trk1Chi2ndf;

  float trk2Pt;
  float trk2Dxy;
  float trk2D0Err;
  float trk2PixelHit;
  float trk2StripHit;
  float trk2Chi2ndf;
  float gen;

  void buildBranch(TTree* nt){
    nt->Branch("mass",&mass);
    nt->Branch("pt",&pt);
    nt->Branch("eta",&eta);
    nt->Branch("y",&y);
    nt->Branch("gen",&gen);
    nt->Branch("px",&px);
    nt->Branch("py",&py);
    nt->Branch("d0",&d0);
    nt->Branch("vx",&vx);
    nt->Branch("vy",&vy);
    nt->Branch("d0Err",&d0Err);
    nt->Branch("chi2ndf",&chi2ndf);
    nt->Branch("dtheta",&dtheta);

    nt->Branch("trk1Pt",&trk1Pt);
    nt->Branch("trk1Dxy",&trk1Dxy);
    nt->Branch("trk1D0Err",&trk1D0Err);
    nt->Branch("trk1PixelHit",&trk1PixelHit);
    nt->Branch("trk1StripHit",&trk1StripHit);
    nt->Branch("trk1Chi2ndf",&trk1Chi2ndf);

    nt->Branch("trk2Pt",&trk2Pt);
    nt->Branch("trk2Dxy",&trk2Dxy);
    nt->Branch("trk2D0Err",&trk2D0Err);
    nt->Branch("trk2PixelHit",&trk2PixelHit);
    nt->Branch("trk2StripHit",&trk2StripHit);
    nt->Branch("trk2Chi2ndf",&trk2Chi2ndf);
  }
};



//Declaration of leaves types
   Int_t           EvtInfo_RunNo;
   Int_t           EvtInfo_EvtNo;
   Int_t           EvtInfo_BxNo;
   Int_t           EvtInfo_LumiNo;
   Int_t           EvtInfo_Orbit;
   Bool_t          EvtInfo_McFlag;
   Int_t           EvtInfo_trgCount;
   Int_t           EvtInfo_nTrgBook;
   Char_t          EvtInfo_trgBook[5842];
   Int_t           EvtInfo_nHLT;
   Bool_t          EvtInfo_hltBits[164];
   Int_t           EvtInfo_nBX;
   Int_t           EvtInfo_BXPU[0];
   Int_t           EvtInfo_nPU[0];
   Float_t         EvtInfo_trueIT[0];
   Double_t        EvtInfo_PVx;
   Double_t        EvtInfo_PVy;
   Double_t        EvtInfo_PVz;
   Double_t        EvtInfo_PVxE;
   Double_t        EvtInfo_PVyE;
   Double_t        EvtInfo_PVzE;
   Double_t        EvtInfo_PVnchi2;
   Double_t        EvtInfo_PVchi2;
   Int_t           MuonInfo_size;
   Int_t           MuonInfo_index[MAX_MUON];
   Int_t           MuonInfo_handle_index[MAX_MUON];
   Int_t           MuonInfo_charge[MAX_MUON];
   Double_t        MuonInfo_pt[MAX_MUON];
   Double_t        MuonInfo_eta[MAX_MUON];
   Double_t        MuonInfo_phi[MAX_MUON];
   Int_t           MuonInfo_i_striphit[MAX_MUON];
   Int_t           MuonInfo_i_pixelhit[MAX_MUON];
   Int_t           MuonInfo_i_nStripLayer[MAX_MUON];
   Int_t           MuonInfo_i_nPixelLayer[MAX_MUON];
   Int_t           MuonInfo_g_striphit[MAX_MUON];
   Int_t           MuonInfo_g_pixelhit[MAX_MUON];
   Double_t        MuonInfo_i_chi2[MAX_MUON];
   Double_t        MuonInfo_i_ndf[MAX_MUON];
   Double_t        MuonInfo_g_chi2[MAX_MUON];
   Double_t        MuonInfo_g_ndf[MAX_MUON];
   Int_t           MuonInfo_nmuhit[MAX_MUON];
   Double_t        MuonInfo_d0[MAX_MUON];
   Double_t        MuonInfo_dz[MAX_MUON];
   Double_t        MuonInfo_dzPV[MAX_MUON];
   Double_t        MuonInfo_dxyPV[MAX_MUON];
   Int_t           MuonInfo_fpbarrelhit[MAX_MUON];
   Int_t           MuonInfo_fpendcaphit[MAX_MUON];
   Int_t           MuonInfo_muqual[MAX_MUON];
   Double_t        MuonInfo_iso_trk[MAX_MUON];
   Double_t        MuonInfo_iso_ecal[MAX_MUON];
   Double_t        MuonInfo_iso_hcal[MAX_MUON];
   Double_t        MuonInfo_n_matches[MAX_MUON];
   Int_t           MuonInfo_isGoodCand[MAX_MUON];
   Int_t           MuonInfo_geninfo_index[MAX_MUON];
   Int_t           TrackInfo_size;
   Int_t           TrackInfo_index[MAX_TRACK];
   Int_t           TrackInfo_handle_index[MAX_TRACK];
   Int_t           TrackInfo_charge[MAX_TRACK];
   Double_t        TrackInfo_pt[MAX_TRACK];
   Double_t        TrackInfo_eta[MAX_TRACK];
   Double_t        TrackInfo_phi[MAX_TRACK];
   Int_t           TrackInfo_striphit[MAX_TRACK];
   Int_t           TrackInfo_pixelhit[MAX_TRACK];
   Int_t           TrackInfo_fpbarrelhit[MAX_TRACK];
   Int_t           TrackInfo_fpendcaphit[MAX_TRACK];
   Double_t        TrackInfo_chi2[MAX_TRACK];
   Double_t        TrackInfo_ndf[MAX_TRACK];
   Double_t        TrackInfo_d0[MAX_TRACK];
   Double_t        TrackInfo_d0error[MAX_TRACK];
   Double_t        TrackInfo_dzPV[MAX_TRACK];
   Double_t        TrackInfo_dxyPV[MAX_TRACK];
   Int_t           TrackInfo_isGoodCand[MAX_TRACK];
   Int_t           TrackInfo_geninfo_index[MAX_TRACK];
   Int_t           BInfo_uj_size;
   Int_t           BInfo_uj_index[MAX_XB];
   Double_t        BInfo_uj_mass[MAX_XB];
   Double_t        BInfo_uj_px[MAX_XB];
   Double_t        BInfo_uj_py[MAX_XB];
   Double_t        BInfo_uj_pz[MAX_XB];
   Double_t        BInfo_uj_vtxX[MAX_XB];
   Double_t        BInfo_uj_vtxY[MAX_XB];
   Double_t        BInfo_uj_vtxZ[MAX_XB];
   Double_t        BInfo_uj_vtxXE[MAX_XB];
   Double_t        BInfo_uj_vtxYE[MAX_XB];
   Double_t        BInfo_uj_vtxZE[MAX_XB];
   Double_t        BInfo_uj_vtxdof[MAX_XB];
   Double_t        BInfo_uj_vtxchi2[MAX_XB];
   Int_t           BInfo_uj_rfmu1_index[MAX_XB];
   Int_t           BInfo_uj_rfmu2_index[MAX_XB];
   Int_t           BInfo_uj_isGoodCand[MAX_XB];
   Double_t        BInfo_uj_rfmu1_px[MAX_XB];
   Double_t        BInfo_uj_rfmu1_py[MAX_XB];
   Double_t        BInfo_uj_rfmu1_pz[MAX_XB];
   Double_t        BInfo_uj_rfmu2_px[MAX_XB];
   Double_t        BInfo_uj_rfmu2_py[MAX_XB];
   Double_t        BInfo_uj_rfmu2_pz[MAX_XB];
   Int_t           BInfo_size;
   Int_t           BInfo_index[MAX_XB];
   Double_t        BInfo_mass[MAX_XB];
   Double_t        BInfo_eta[MAX_XB];
   Double_t        BInfo_px[MAX_XB];
   Double_t        BInfo_py[MAX_XB];
   Double_t        BInfo_pz[MAX_XB];
   Double_t        BInfo_pxE[MAX_XB];
   Double_t        BInfo_pyE[MAX_XB];
   Double_t        BInfo_pzE[MAX_XB];
   Double_t        BInfo_vtxX[MAX_XB];
   Double_t        BInfo_vtxY[MAX_XB];
   Double_t        BInfo_vtxZ[MAX_XB];
   Double_t        BInfo_vtxXE[MAX_XB];
   Double_t        BInfo_vtxYE[MAX_XB];
   Double_t        BInfo_vtxZE[MAX_XB];
   Double_t        BInfo_vtxdof[MAX_XB];
   Double_t        BInfo_vtxchi2[MAX_XB];
   Int_t           BInfo_rfuj_index[MAX_XB];
   Int_t           BInfo_rftk1_index[MAX_XB];
   Int_t           BInfo_rftk2_index[MAX_XB];
   Int_t           BInfo_isGoodCand[MAX_XB];
   Int_t           BInfo_type[MAX_XB];
   Double_t        BInfo_rfmu1_px[MAX_XB];
   Double_t        BInfo_rfmu1_py[MAX_XB];
   Double_t        BInfo_rfmu1_pz[MAX_XB];
   Double_t        BInfo_rfmu2_px[MAX_XB];
   Double_t        BInfo_rfmu2_py[MAX_XB];
   Double_t        BInfo_rfmu2_pz[MAX_XB];
   Double_t        BInfo_rftk1_px[MAX_XB];
   Double_t        BInfo_rftk1_py[MAX_XB];
   Double_t        BInfo_rftk1_pz[MAX_XB];
   Double_t        BInfo_rftk2_px[MAX_XB];
   Double_t        BInfo_rftk2_py[MAX_XB];
   Double_t        BInfo_rftk2_pz[MAX_XB];
   Double_t        BInfo_tktk_mass[MAX_XB];
   Double_t        BInfo_tktk_px[MAX_XB];
   Double_t        BInfo_tktk_py[MAX_XB];
   Double_t        BInfo_tktk_pz[MAX_XB];
   Double_t        BInfo_tktk_vtxX[MAX_XB];
   Double_t        BInfo_tktk_vtxY[MAX_XB];
   Double_t        BInfo_tktk_vtxZ[MAX_XB];
   Double_t        BInfo_tktk_vtxXE[MAX_XB];
   Double_t        BInfo_tktk_vtxYE[MAX_XB];
   Double_t        BInfo_tktk_vtxZE[MAX_XB];
   Double_t        BInfo_tktk_vtxdof[MAX_XB];
   Double_t        BInfo_tktk_vtxchi2[MAX_XB];
   Int_t           GenInfo_size;
   Int_t           GenInfo_index[MAX_GEN];
   Int_t           GenInfo_handle_index[MAX_GEN];
   Double_t        GenInfo_pt[MAX_GEN];
   Double_t        GenInfo_eta[MAX_GEN];
   Double_t        GenInfo_phi[MAX_GEN];
   Double_t        GenInfo_mass[MAX_GEN];
   Int_t           GenInfo_pdgId[MAX_GEN];
   Int_t           GenInfo_status[MAX_GEN];
   Int_t           GenInfo_nMo[MAX_GEN];
   Int_t           GenInfo_nDa[MAX_GEN];
   Int_t           GenInfo_mo1[MAX_GEN];
   Int_t           GenInfo_mo2[MAX_GEN];
   Int_t           GenInfo_da1[MAX_GEN];
   Int_t           GenInfo_da2[MAX_GEN];
//   Int_t           GenInfo_mhmu1_index;
//   Int_t           GenInfo_mhmu2_index;
//   Int_t           GenInfo_mhtk1_index;
//   Int_t           GenInfo_mhtk2_index;
//   Double_t        GenInfo_mhujMass;
//   Double_t        GenInfo_mhxbMass;


void setBranch(TTree *root) {
   // Set branch addresses.
  /*
   root->SetBranchAddress("EvtInfo.RunNo",&EvtInfo_RunNo);
   root->SetBranchAddress("EvtInfo.EvtNo",&EvtInfo_EvtNo);
   root->SetBranchAddress("EvtInfo.BxNo",&EvtInfo_BxNo);
   root->SetBranchAddress("EvtInfo.LumiNo",&EvtInfo_LumiNo);
   root->SetBranchAddress("EvtInfo.Orbit",&EvtInfo_Orbit);
   root->SetBranchAddress("EvtInfo.McFlag",&EvtInfo_McFlag);
   root->SetBranchAddress("EvtInfo.trgCount",&EvtInfo_trgCount);
   root->SetBranchAddress("EvtInfo.nTrgBook",&EvtInfo_nTrgBook);
   root->SetBranchAddress("EvtInfo.trgBook",EvtInfo_trgBook);
   root->SetBranchAddress("EvtInfo.nHLT",&EvtInfo_nHLT);
   root->SetBranchAddress("EvtInfo.hltBits",EvtInfo_hltBits);
   root->SetBranchAddress("EvtInfo.nBX",&EvtInfo_nBX);
   root->SetBranchAddress("EvtInfo.BXPU",&EvtInfo_BXPU);
   root->SetBranchAddress("EvtInfo.nPU",&EvtInfo_nPU);
   root->SetBranchAddress("EvtInfo.trueIT",&EvtInfo_trueIT);
  */  
   root->SetBranchAddress("EvtInfo.PVx",&EvtInfo_PVx);
   root->SetBranchAddress("EvtInfo.PVy",&EvtInfo_PVy);
   root->SetBranchAddress("EvtInfo.PVz",&EvtInfo_PVz);
   root->SetBranchAddress("EvtInfo.PVxE",&EvtInfo_PVxE);
   root->SetBranchAddress("EvtInfo.PVyE",&EvtInfo_PVyE);
   root->SetBranchAddress("EvtInfo.PVzE",&EvtInfo_PVzE);
   root->SetBranchAddress("EvtInfo.PVnchi2",&EvtInfo_PVnchi2);
   root->SetBranchAddress("EvtInfo.PVchi2",&EvtInfo_PVchi2);
   

   root->SetBranchAddress("MuonInfo.size",&MuonInfo_size);
   root->SetBranchAddress("MuonInfo.index",MuonInfo_index);
   root->SetBranchAddress("MuonInfo.handle_index",MuonInfo_handle_index);
   root->SetBranchAddress("MuonInfo.charge",MuonInfo_charge);
   root->SetBranchAddress("MuonInfo.pt",MuonInfo_pt);
   root->SetBranchAddress("MuonInfo.eta",MuonInfo_eta);
   root->SetBranchAddress("MuonInfo.phi",MuonInfo_phi);
   root->SetBranchAddress("MuonInfo.i_striphit",MuonInfo_i_striphit);
   root->SetBranchAddress("MuonInfo.i_pixelhit",MuonInfo_i_pixelhit);
   root->SetBranchAddress("MuonInfo.i_nStripLayer",MuonInfo_i_nStripLayer);
   root->SetBranchAddress("MuonInfo.i_nPixelLayer",MuonInfo_i_nPixelLayer);
   root->SetBranchAddress("MuonInfo.g_striphit",MuonInfo_g_striphit);
   root->SetBranchAddress("MuonInfo.g_pixelhit",MuonInfo_g_pixelhit);
   root->SetBranchAddress("MuonInfo.i_chi2",MuonInfo_i_chi2);
   root->SetBranchAddress("MuonInfo.i_ndf",MuonInfo_i_ndf);
   root->SetBranchAddress("MuonInfo.g_chi2",MuonInfo_g_chi2);
   root->SetBranchAddress("MuonInfo.g_ndf",MuonInfo_g_ndf);
   root->SetBranchAddress("MuonInfo.nmuhit",MuonInfo_nmuhit);
   root->SetBranchAddress("MuonInfo.d0",MuonInfo_d0);
   root->SetBranchAddress("MuonInfo.dz",MuonInfo_dz);
   root->SetBranchAddress("MuonInfo.dzPV",MuonInfo_dzPV);
   root->SetBranchAddress("MuonInfo.dxyPV",MuonInfo_dxyPV);
   root->SetBranchAddress("MuonInfo.fpbarrelhit",MuonInfo_fpbarrelhit);
   root->SetBranchAddress("MuonInfo.fpendcaphit",MuonInfo_fpendcaphit);
   root->SetBranchAddress("MuonInfo.muqual",MuonInfo_muqual);
   root->SetBranchAddress("MuonInfo.iso_trk",MuonInfo_iso_trk);
   root->SetBranchAddress("MuonInfo.iso_ecal",MuonInfo_iso_ecal);
   root->SetBranchAddress("MuonInfo.iso_hcal",MuonInfo_iso_hcal);
   root->SetBranchAddress("MuonInfo.n_matches",MuonInfo_n_matches);
   root->SetBranchAddress("MuonInfo.isGoodCand",MuonInfo_isGoodCand);
   root->SetBranchAddress("MuonInfo.geninfo_index",MuonInfo_geninfo_index);
   

   root->SetBranchAddress("TrackInfo.size",&TrackInfo_size);
   root->SetBranchAddress("TrackInfo.index",TrackInfo_index);
   root->SetBranchAddress("TrackInfo.handle_index",TrackInfo_handle_index);
   root->SetBranchAddress("TrackInfo.charge",TrackInfo_charge);
   root->SetBranchAddress("TrackInfo.pt",TrackInfo_pt);
   root->SetBranchAddress("TrackInfo.eta",TrackInfo_eta);
   root->SetBranchAddress("TrackInfo.phi",TrackInfo_phi);
   root->SetBranchAddress("TrackInfo.striphit",TrackInfo_striphit);
   root->SetBranchAddress("TrackInfo.pixelhit",TrackInfo_pixelhit);
   root->SetBranchAddress("TrackInfo.fpbarrelhit",TrackInfo_fpbarrelhit);
   root->SetBranchAddress("TrackInfo.fpendcaphit",TrackInfo_fpendcaphit);
   root->SetBranchAddress("TrackInfo.chi2",TrackInfo_chi2);
   root->SetBranchAddress("TrackInfo.ndf",TrackInfo_ndf);
   root->SetBranchAddress("TrackInfo.d0",TrackInfo_d0);
   root->SetBranchAddress("TrackInfo.d0error",TrackInfo_d0error);
   root->SetBranchAddress("TrackInfo.dzPV",TrackInfo_dzPV);
   root->SetBranchAddress("TrackInfo.dxyPV",TrackInfo_dxyPV);
   root->SetBranchAddress("TrackInfo.isGoodCand",TrackInfo_isGoodCand);
   root->SetBranchAddress("TrackInfo.geninfo_index",TrackInfo_geninfo_index);
   
   root->SetBranchAddress("BInfo.uj_size",&BInfo_uj_size);
   root->SetBranchAddress("BInfo.uj_index",BInfo_uj_index);
   root->SetBranchAddress("BInfo.uj_mass",BInfo_uj_mass);
   root->SetBranchAddress("BInfo.uj_px",BInfo_uj_px);
   root->SetBranchAddress("BInfo.uj_py",BInfo_uj_py);
   root->SetBranchAddress("BInfo.uj_pz",BInfo_uj_pz);
   root->SetBranchAddress("BInfo.uj_vtxX",BInfo_uj_vtxX);
   root->SetBranchAddress("BInfo.uj_vtxY",BInfo_uj_vtxY);
   root->SetBranchAddress("BInfo.uj_vtxZ",BInfo_uj_vtxZ);
   root->SetBranchAddress("BInfo.uj_vtxXE",BInfo_uj_vtxXE);
   root->SetBranchAddress("BInfo.uj_vtxYE",BInfo_uj_vtxYE);
   root->SetBranchAddress("BInfo.uj_vtxZE",BInfo_uj_vtxZE);
   root->SetBranchAddress("BInfo.uj_vtxdof",BInfo_uj_vtxdof);
   root->SetBranchAddress("BInfo.uj_vtxchi2",BInfo_uj_vtxchi2);
   root->SetBranchAddress("BInfo.uj_rfmu1_index",BInfo_uj_rfmu1_index);
   root->SetBranchAddress("BInfo.uj_rfmu2_index",BInfo_uj_rfmu2_index);
   root->SetBranchAddress("BInfo.uj_isGoodCand",BInfo_uj_isGoodCand);
   root->SetBranchAddress("BInfo.uj_rfmu1_px",BInfo_uj_rfmu1_px);
   root->SetBranchAddress("BInfo.uj_rfmu1_py",BInfo_uj_rfmu1_py);
   root->SetBranchAddress("BInfo.uj_rfmu1_pz",BInfo_uj_rfmu1_pz);
   root->SetBranchAddress("BInfo.uj_rfmu2_px",BInfo_uj_rfmu2_px);
   root->SetBranchAddress("BInfo.uj_rfmu2_py",BInfo_uj_rfmu2_py);
   root->SetBranchAddress("BInfo.uj_rfmu2_pz",BInfo_uj_rfmu2_pz);
   

   root->SetBranchAddress("BInfo.size",&BInfo_size);
   root->SetBranchAddress("BInfo.index",BInfo_index);
   root->SetBranchAddress("BInfo.mass",BInfo_mass);
   root->SetBranchAddress("BInfo.eta",BInfo_eta);
   root->SetBranchAddress("BInfo.px",BInfo_px);
   root->SetBranchAddress("BInfo.py",BInfo_py);
   root->SetBranchAddress("BInfo.pz",BInfo_pz);
   root->SetBranchAddress("BInfo.pxE",BInfo_pxE);
   root->SetBranchAddress("BInfo.pyE",BInfo_pyE);
   root->SetBranchAddress("BInfo.pzE",BInfo_pzE);
   root->SetBranchAddress("BInfo.vtxX",BInfo_vtxX);
   root->SetBranchAddress("BInfo.vtxY",BInfo_vtxY);
   root->SetBranchAddress("BInfo.vtxZ",BInfo_vtxZ);
   root->SetBranchAddress("BInfo.vtxXE",BInfo_vtxXE);
   root->SetBranchAddress("BInfo.vtxYE",BInfo_vtxYE);
   root->SetBranchAddress("BInfo.vtxZE",BInfo_vtxZE);
   root->SetBranchAddress("BInfo.vtxdof",BInfo_vtxdof);
   root->SetBranchAddress("BInfo.vtxchi2",BInfo_vtxchi2);

   root->SetBranchAddress("BInfo.rfuj_index",BInfo_rfuj_index);
   root->SetBranchAddress("BInfo.rftk1_index",BInfo_rftk1_index);
   root->SetBranchAddress("BInfo.rftk2_index",BInfo_rftk2_index);
   root->SetBranchAddress("BInfo.isGoodCand",BInfo_isGoodCand);
   root->SetBranchAddress("BInfo.type",BInfo_type);
   
   root->SetBranchAddress("BInfo.rfmu1_px",BInfo_rfmu1_px);
   root->SetBranchAddress("BInfo.rfmu1_py",BInfo_rfmu1_py);
   root->SetBranchAddress("BInfo.rfmu1_pz",BInfo_rfmu1_pz);
   root->SetBranchAddress("BInfo.rfmu2_px",BInfo_rfmu2_px);
   root->SetBranchAddress("BInfo.rfmu2_py",BInfo_rfmu2_py);
   root->SetBranchAddress("BInfo.rfmu2_pz",BInfo_rfmu2_pz);

   root->SetBranchAddress("BInfo.rftk1_px",BInfo_rftk1_px);
   root->SetBranchAddress("BInfo.rftk1_py",BInfo_rftk1_py);
   root->SetBranchAddress("BInfo.rftk1_pz",BInfo_rftk1_pz);
   root->SetBranchAddress("BInfo.rftk2_px",BInfo_rftk2_px);
   root->SetBranchAddress("BInfo.rftk2_py",BInfo_rftk2_py);
   root->SetBranchAddress("BInfo.rftk2_pz",BInfo_rftk2_pz);
      
   root->SetBranchAddress("BInfo.tktk_mass",BInfo_tktk_mass);
   root->SetBranchAddress("BInfo.tktk_px",BInfo_tktk_px);
   root->SetBranchAddress("BInfo.tktk_py",BInfo_tktk_py);
   root->SetBranchAddress("BInfo.tktk_pz",BInfo_tktk_pz);
   root->SetBranchAddress("BInfo.tktk_vtxX",BInfo_tktk_vtxX);
   root->SetBranchAddress("BInfo.tktk_vtxY",BInfo_tktk_vtxY);
   root->SetBranchAddress("BInfo.tktk_vtxZ",BInfo_tktk_vtxZ);
   root->SetBranchAddress("BInfo.tktk_vtxXE",BInfo_tktk_vtxXE);
   root->SetBranchAddress("BInfo.tktk_vtxYE",BInfo_tktk_vtxYE);
   root->SetBranchAddress("BInfo.tktk_vtxZE",BInfo_tktk_vtxZE);
   root->SetBranchAddress("BInfo.tktk_vtxdof",BInfo_tktk_vtxdof);
   root->SetBranchAddress("BInfo.tktk_vtxchi2",BInfo_tktk_vtxchi2);
   
   root->SetBranchAddress("GenInfo.size",&GenInfo_size);
   root->SetBranchAddress("GenInfo.index",&GenInfo_index);
   root->SetBranchAddress("GenInfo.handle_index",&GenInfo_handle_index);
   root->SetBranchAddress("GenInfo.pt",&GenInfo_pt);
   root->SetBranchAddress("GenInfo.eta",&GenInfo_eta);
   root->SetBranchAddress("GenInfo.phi",&GenInfo_phi);
   root->SetBranchAddress("GenInfo.mass",&GenInfo_mass);
   root->SetBranchAddress("GenInfo.pdgId",&GenInfo_pdgId);
   root->SetBranchAddress("GenInfo.status",&GenInfo_status);
   root->SetBranchAddress("GenInfo.nMo",&GenInfo_nMo);
   root->SetBranchAddress("GenInfo.nDa",&GenInfo_nDa);
   root->SetBranchAddress("GenInfo.mo1",&GenInfo_mo1);
   root->SetBranchAddress("GenInfo.mo2",&GenInfo_mo2);
   root->SetBranchAddress("GenInfo.da1",&GenInfo_da1);
   root->SetBranchAddress("GenInfo.da2",&GenInfo_da2);
   //root->SetBranchAddress("Geninfo.mhmu1_index",&GenInfo_mhmu1_index);
   //root->SetBranchAddress("Geninfo.mhmu2_index",&GenInfo_mhmu2_index);
   //root->SetBranchAddress("Geninfo.mhtk1_index",&GenInfo_mhtk1_index);
   //root->SetBranchAddress("Geninfo.mhtk2_index",&GenInfo_mhtk2_index);
   //root->SetBranchAddress("GenInfo.mhujMass",&GenInfo_mhujMass);
   //root->SetBranchAddress("GenInfo.mhxbMass",&GenInfo_mhxbMass);

//     This is the loop skeleton
//       To read only selected branches, Insert statements like:
// root->SetBranchStatus("*",0);  // disable all branches
// TTreePlayer->SetBranchStatus("branchname",1);  // activate branchname


}
