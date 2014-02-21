#include <TTree.h>
#include <TFile.h>
#include <iostream>
#include <TNtuple.h>
#include <TVector3.h>
#include <TLorentzVector.h>

#include "loop.h"

void fillTree(bNtuple* b, TVector3* bP, TVector3* bVtx, int j)
{
  bP->SetXYZ(BInfo_px[j],BInfo_py[j],BInfo_pz[j]*0);
  bVtx->SetXYZ(BInfo_vtxX[j]-EvtInfo_PVx,
	       BInfo_vtxY[j]-EvtInfo_PVy,
	       BInfo_vtxZ[j]*0-EvtInfo_PVz*0);
  TLorentzVector b4P;
  b4P.SetXYZM(BInfo_px[j],BInfo_py[j],BInfo_pz[j],BInfo_mass[j]);
  b->y = b4P.Rapidity();
  b->dtheta = bP->Angle(*bVtx);
  b->pt = sqrt(BInfo_px[j]*BInfo_px[j]+BInfo_py[j]*BInfo_py[j]);
  b->eta = BInfo_eta[j];
  b->px = BInfo_px[j];
  b->py = BInfo_py[j];
  b->d0 = sqrt((BInfo_vtxX[j]-EvtInfo_PVx)*(BInfo_vtxX[j]-EvtInfo_PVx)+(BInfo_vtxY[j]-EvtInfo_PVy)*(BInfo_vtxY[j]-EvtInfo_PVy));
  b->vx = BInfo_vtxX[j] - EvtInfo_PVx;
  b->vy = BInfo_vtxY[j] - EvtInfo_PVy;
  b->d0Err = sqrt(BInfo_vtxXE[j]*BInfo_vtxXE[j]+BInfo_vtxYE[j]*BInfo_vtxYE[j]);
  b->mass = BInfo_mass[j];
  b->chi2ndf = BInfo_vtxchi2[j]/BInfo_vtxdof[j];
  
  b->trk1Dxy = TrackInfo_dxyPV[BInfo_rftk1_index[j]];
  b->trk1D0Err = TrackInfo_d0error[BInfo_rftk1_index[j]];
  b->trk1PixelHit = TrackInfo_pixelhit[BInfo_rftk1_index[j]];
  b->trk1StripHit = TrackInfo_striphit[BInfo_rftk1_index[j]];
  b->trk1Pt = TrackInfo_pt[BInfo_rftk1_index[j]];
  b->trk1Chi2ndf = TrackInfo_chi2[BInfo_rftk1_index[j]]/TrackInfo_ndf[BInfo_rftk1_index[j]];
  
  if(BInfo_type[j]==1 || BInfo_type[j]==2)
    {
      b->trk2Dxy = 0;
      b->trk2D0Err = 0;
      b->trk2PixelHit = 0;
      b->trk2StripHit = 0;
      b->trk2Pt = 0;
      b->trk2Chi2ndf = 0;
    }  
  else
    {
      b->trk2Dxy = TrackInfo_dxyPV[BInfo_rftk2_index[j]];
      b->trk2D0Err = TrackInfo_d0error[BInfo_rftk2_index[j]];
      b->trk2PixelHit = TrackInfo_pixelhit[BInfo_rftk2_index[j]];
      b->trk2StripHit = TrackInfo_striphit[BInfo_rftk2_index[j]];
      b->trk2Pt = TrackInfo_pt[BInfo_rftk2_index[j]];
      b->trk2Chi2ndf = TrackInfo_chi2[BInfo_rftk2_index[j]]/TrackInfo_ndf[BInfo_rftk2_index[j]];
    }

  
  //gen info judgement

  b->gen=0;//gen init
  int mGenIdxTk1=-1;
  int mGenIdxTk2=-1;
  int bGenIdxTk1=-1;
  int bGenIdxTk2=-1;
  int bGenIdxMu1=-1;
  int bGenIdxMu2=-1;


  float BId,MId,tk1Id,tk2Id;
  //tk1:positive, tk2:negtive
  if(BInfo_type[j]==1)
    {
      BId = 521;//B+-
      MId = -1;
      tk1Id = 321;//K+-
      tk2Id = -1;
    }
  if(BInfo_type[j]==2)
    {
      BId = 521;//B+-
      MId = -1;
      tk1Id = 211;//pi+-
      tk2Id = -1;
    }
  if(BInfo_type[j]==3)
    {
      BId = 511;//B0
      MId = 310;//Ks
      tk1Id = 211;//pi+
      tk2Id = 211;//pi-
    }
  if(BInfo_type[j]==4)
    {
      BId = 511;//B0
      MId = 313;//K*0
      tk1Id = 321;//K+
      tk2Id = 211;//pi-
    }
  if(BInfo_type[j]==5)
    {
      BId = 511;//B0
      MId = 313;//K*0
      tk1Id = 211;//pi+
      tk2Id = 321;//K-
    }
  if(BInfo_type[j]==6)
    {
      BId = 531;//Bs
      MId = 333;//phi
      tk1Id = 321;//K+
      tk2Id = 321;//K-
    }

  int twoTks;
  if(BInfo_type[j]==1 || BInfo_type[j]==2) twoTks=0;
  else twoTks=1;

  // tk1
  if(TrackInfo_geninfo_index[BInfo_rftk1_index[j]]>-1)
    {
      int level =0;
      if(fabs(GenInfo_pdgId[TrackInfo_geninfo_index[BInfo_rftk1_index[j]]])==tk1Id)
	{
	  level = 1;
	  if(GenInfo_mo1[TrackInfo_geninfo_index[BInfo_rftk1_index[j]]]>-1)
	    {
	      if(!twoTks)//one trk channel
		{
		  mGenIdxTk1=0;
		  if(fabs(GenInfo_pdgId[GenInfo_mo1[TrackInfo_geninfo_index[BInfo_rftk1_index[j]]]])==BId)
		    {
		      level = 3;
		      bGenIdxTk1=GenInfo_mo1[TrackInfo_geninfo_index[BInfo_rftk1_index[j]]];
		    }		  
		}
	      else//two trk channel
		{
		  if(fabs(GenInfo_pdgId[GenInfo_mo1[TrackInfo_geninfo_index[BInfo_rftk1_index[j]]]])==MId)
		    {
		      level = 2;
		      if(GenInfo_mo1[GenInfo_mo1[TrackInfo_geninfo_index[BInfo_rftk1_index[j]]]]>-1)
			{
			  if(fabs(GenInfo_pdgId[GenInfo_mo1[GenInfo_mo1[TrackInfo_geninfo_index[BInfo_rftk1_index[j]]]]])==BId)
			    {
			      level = 3;
			      bGenIdxTk1=GenInfo_mo1[GenInfo_mo1[TrackInfo_geninfo_index[BInfo_rftk1_index[j]]]];
			    }
			}
		      mGenIdxTk1=GenInfo_mo1[TrackInfo_geninfo_index[BInfo_rftk1_index[j]]];
		    }
		}
	    }
	}
      b->gen=level;
    }

  //tk2
  if(!twoTks)//one trk channel
    {
      b->gen+=30;
      mGenIdxTk2=0;
      bGenIdxTk2=0;
    }
  else//two trk channel
    {
      if(TrackInfo_geninfo_index[BInfo_rftk2_index[j]]>-1)
	{
	  int level =0;
	  if(fabs(GenInfo_pdgId[TrackInfo_geninfo_index[BInfo_rftk2_index[j]]])==tk2Id)
	    {
	      level = 1;
	      if(GenInfo_mo1[TrackInfo_geninfo_index[BInfo_rftk2_index[j]]]>-1)
		{
		  if(fabs(GenInfo_pdgId[GenInfo_mo1[TrackInfo_geninfo_index[BInfo_rftk2_index[j]]]])==MId)
		    {
		      level = 2;
		      if(GenInfo_mo1[GenInfo_mo1[TrackInfo_geninfo_index[BInfo_rftk2_index[j]]]]>-1)
			{
			  if(fabs(GenInfo_pdgId[GenInfo_mo1[GenInfo_mo1[TrackInfo_geninfo_index[BInfo_rftk2_index[j]]]]])==BId)
			    {
			      level = 3;
			      bGenIdxTk2 = GenInfo_mo1[GenInfo_mo1[TrackInfo_geninfo_index[BInfo_rftk2_index[j]]]];
			    }
			}
		      mGenIdxTk2 = GenInfo_mo1[TrackInfo_geninfo_index[BInfo_rftk2_index[j]]];
		    }
		}
	    }
	  b->gen+=(level*10);
	}
    }

  //mu1
  if(MuonInfo_geninfo_index[BInfo_uj_rfmu1_index[BInfo_rfuj_index[j]]]>-1)
    {  
      int level =0;
      if(fabs(GenInfo_pdgId[MuonInfo_geninfo_index[BInfo_uj_rfmu1_index[BInfo_rfuj_index[j]]]])==13) level=1;
      if(GenInfo_mo1[MuonInfo_geninfo_index[BInfo_uj_rfmu1_index[BInfo_rfuj_index[j]]]]>-1)
	{
	  if(GenInfo_mo1[GenInfo_mo1[MuonInfo_geninfo_index[BInfo_uj_rfmu1_index[BInfo_rfuj_index[j]]]]]>-1)
	    {
	      if(fabs(GenInfo_pdgId[GenInfo_mo1[GenInfo_mo1[MuonInfo_geninfo_index[BInfo_uj_rfmu1_index[BInfo_rfuj_index[j]]]]]])==BId)
		{
		  level = 2;
		  bGenIdxMu1=GenInfo_mo1[GenInfo_mo1[MuonInfo_geninfo_index[BInfo_uj_rfmu1_index[BInfo_rfuj_index[j]]]]];
		}
	    } 
	}
      b->gen+=(level*100);
    }

  //mu2
  if(MuonInfo_geninfo_index[BInfo_uj_rfmu2_index[BInfo_rfuj_index[j]]]>-1)
    {  
      int level =0;
      if(fabs(GenInfo_pdgId[MuonInfo_geninfo_index[BInfo_uj_rfmu2_index[BInfo_rfuj_index[j]]]])==13) level = 1;
      if(GenInfo_mo1[MuonInfo_geninfo_index[BInfo_uj_rfmu2_index[BInfo_rfuj_index[j]]]]>-1)
	{
	  if(GenInfo_mo1[GenInfo_mo1[MuonInfo_geninfo_index[BInfo_uj_rfmu2_index[BInfo_rfuj_index[j]]]]]>-1)
	    {
	      if(fabs(GenInfo_pdgId[GenInfo_mo1[GenInfo_mo1[MuonInfo_geninfo_index[BInfo_uj_rfmu2_index[BInfo_rfuj_index[j]]]]]])==BId)
		{
		  level = 2;
		  bGenIdxMu2=GenInfo_mo1[GenInfo_mo1[MuonInfo_geninfo_index[BInfo_uj_rfmu2_index[BInfo_rfuj_index[j]]]]];
		}
	    }
	}
      b->gen+=(level*1000);
    }
  
  int level=0;
  if(mGenIdxTk1!=-1 && mGenIdxTk2!=-1)
    {
      if(!twoTks) level=1;
      else
	{
	  if(mGenIdxTk1==mGenIdxTk2) level=1;
	}
    }
  if(bGenIdxMu1!=-1 && bGenIdxMu1==bGenIdxMu2 && bGenIdxMu1==bGenIdxTk1)
    {
      if(!twoTks) level=2;
      else if(bGenIdxMu1==bGenIdxTk2) level=2;
    }
  b->gen+=(level*10000);
  
}

int signalGen(int Btype, int j)
{
  float BId,MId,tk1Id,tk2Id;
  int twoTks;
  //tk1:positive, tk2:negtive
  if(Btype==1)
    {
      BId = 521;//B+-
      MId = -1;
      tk1Id = 321;//K+-
      tk2Id = -1;
      twoTks = 0;
    }
  if(Btype==2)
    {
      BId = 521;//B+-
      MId = -1;
      tk1Id = 211;//pi+-
      tk2Id = -1;
      twoTks = 0;
    }
  if(Btype==3)
    {
      BId = 511;//B0
      MId = 310;//Ks
      tk1Id = 211;//pi+
      tk2Id = 211;//pi-
      twoTks = 1;
    }
  if(Btype==4)
    {
      BId = 511;//B0
      MId = 313;//K*0
      tk1Id = 321;//K+
      tk2Id = 211;//pi-
      twoTks = 1;
    }
  if(Btype==5)
    {
      BId = 511;//B0
      MId = 313;//K*0
      tk1Id = 211;//pi+
      tk2Id = 321;//K-
      twoTks = 1;
    }
  if(Btype==6)
    {
      BId = 531;//Bs
      MId = 333;//phi
      tk1Id = 321;//K+
      tk2Id = 321;//K-
      twoTks = 1;
    }

  int flag=0;
  if (fabs(GenInfo_pdgId[j])==BId&&GenInfo_nDa[j]==2&&GenInfo_da1[j]!=-1&&GenInfo_da2[j]!=-1)
    {
      if (fabs(GenInfo_pdgId[GenInfo_da1[j]]==443))//jpsi
	{
	  if(!twoTks)
	    {
	      if(fabs(GenInfo_pdgId[GenInfo_da2[j]])==tk1Id) flag++;
	    }
	  else
	    {
	      if (fabs(GenInfo_pdgId[GenInfo_da2[j]])==MId) 
		{
		  if(GenInfo_da1[GenInfo_da2[j]]!=-1 && GenInfo_da2[GenInfo_da2[j]]!=-1)
		    {
		      if(fabs(GenInfo_pdgId[GenInfo_da1[GenInfo_da2[j]]])==tk1Id && fabs(GenInfo_pdgId[GenInfo_da2[GenInfo_da2[j]]])==tk2Id) flag++;
		    }
		}
	    }
	}
    }
  return flag;
}

//////////////////Taking in/out file name as input
void loop(string infile, string outfile, bool REAL=0){
//void loop(bool REAL=0){

const   char* infname;
const   char* outfname;
/////////////////

   if(REAL)
     {
      cout<<"--- REAL DATA ---"<<endl;
      infname = "/net/hidsk0001/d00/scratch/yjlee/bmeson/merged_pPbData_20131114.root";
      outfname = "nt_data.root";
     }
   else
     {
      cout<<"--- MC ---"<<endl;
	  //infname = "/mnt/hadoop/cms/store/user/wangj/HI_Btuple/20140218_PAMuon_HIRun2013_PromptReco_v1/Bfinder_all_100_1_Jrd.root";
      //outfname = "nt_mc.root";

//////////////////
	  infname = infile.c_str();
	  outfname = outfile.c_str();
/////////////////

     }

   //File type
   TFile *f = new TFile(infname);
   TTree *root = (TTree*)f->Get("demo/root");

   //Chain type
   //TChain* root = new TChain("demo/root");
   //root->Add("/mnt/hadoop/cms/store/user/twang/HI_Btuple/20131202_PPMuon_Run2013A-PromptReco-v1_RECO/Bfinder_all_*");

   setBranch(root);
   TFile *outf = new TFile(outfname,"recreate");

   int ifchannel[7];
   ifchannel[0] = 1; //jpsi+Kp
   ifchannel[1] = 1; //jpsi+pi
   ifchannel[2] = 1; //jpsi+Ks(pi+,pi-)
   ifchannel[3] = 1; //jpsi+K*(K+,pi-)
   ifchannel[4] = 1; //jpsi+K*(K-,pi+)
   ifchannel[5] = 1; //jpsi+phi
   ifchannel[6] = 1; //jpsi+pi pi <= psi', X(3872), Bs->J/psi f0
   bNtuple* b0 = new bNtuple;
   bNtuple* b1 = new bNtuple;
   bNtuple* b2 = new bNtuple;
   bNtuple* b3 = new bNtuple;
   bNtuple* b4 = new bNtuple;
   bNtuple* b5 = new bNtuple;
   bNtuple* b6 = new bNtuple;
      
   TTree* nt0 = new TTree("ntKp","");
   b0->buildBranch(nt0);
   TTree* nt1 = new TTree("ntpi","");
   b1->buildBranch(nt1);
   TTree* nt2 = new TTree("ntKs","");
   b2->buildBranch(nt2);
   TTree* nt3 = new TTree("ntKstar1","");
   b3->buildBranch(nt3);
   TTree* nt4 = new TTree("ntKstar2","");
   b4->buildBranch(nt4);
   TTree* nt5 = new TTree("ntphi","");
   b5->buildBranch(nt5);
   TTree* nt6 = new TTree("ntmix","");
   b6->buildBranch(nt6);

   TNtuple* ntGen = new TNtuple("ntGen","","y:eta:phi:pt:pdgId");

   Long64_t nentries = root->GetEntries();
   Long64_t nbytes = 0;
   TVector3* bP = new TVector3;
   TVector3* bVtx = new TVector3;
   TLorentzVector bGen;
   int type;

   for (Long64_t i=0; i<100;i++) {
//   for (Long64_t i=0; i<nentries;i++) {
      nbytes += root->GetEntry(i);
      if (i%10000==0) cout <<i<<" / "<<nentries<<endl;
      for (int j=0;j<BInfo_size;j++) {
	if(BInfo_type[j]>7) continue;
	if (ifchannel[BInfo_type[j]-1]!=1) continue;
	if(BInfo_type[j]==1)
	  {
	    fillTree(b0,bP,bVtx,j);
	    nt0->Fill();
	  }
	if(BInfo_type[j]==2)
	  {
	    fillTree(b1,bP,bVtx,j);
	    nt1->Fill();
	  }
	if(BInfo_type[j]==3)
	  {
	    fillTree(b2,bP,bVtx,j);
	    nt2->Fill();
	  }
	if(BInfo_type[j]==4)
	  {
	    fillTree(b3,bP,bVtx,j);
	    nt3->Fill();
	  }
	if(BInfo_type[j]==5)
	  {
	    fillTree(b4,bP,bVtx,j);
	    nt4->Fill();
	  }
	if(BInfo_type[j]==6)
	  {
	    fillTree(b5,bP,bVtx,j);
	    nt5->Fill();
	  }
	if(BInfo_type[j]==7)
	  {
	    fillTree(b6,bP,bVtx,j);
	    nt6->Fill();
	  }
      }

      for (int j=0;j<GenInfo_size;j++)
	{
	  for(type=1;type<8;type++)
	    {
	      if(signalGen(type,j))
		{
		  bGen.SetPtEtaPhiM(GenInfo_pt[j],GenInfo_eta[j],GenInfo_phi[j],GenInfo_mass[j]);
		  ntGen->Fill(bGen.Rapidity(),bGen.Eta(),bGen.Phi(),bGen.Pt(),GenInfo_pdgId[j]);
		  break;
		}
	    }
	}
   }

  outf->Write();
  outf->Close();
}
