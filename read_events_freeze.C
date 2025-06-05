//************************************************************************************************************************************************
/** author: Susanne Glaessel (Universitaet Frankfurt)
 ** Example macro for reading EventFreeze made from PHQMD-output with convert_phqmd_unigen_freezeout.C.
 **/

#include "TROOT.h"
#include "TFile.h"
#include "TMath.h"
#include "TString.h"
#include <iostream>

class ParticleFreeze {
public:
  Int_t    fPdgId;        // PDG code
  TVector3 fP;            // 3-momentum at final time (pz, py, pz) (GeV/c)
  Float_t  fEnergy;       // Energy (GeV/c)
  Float_t  fTimeFreeze;   // Freezeout-time (fm/c)
  TVector3 fXFreeze;      // Position at freezeout-time (fm)
  TVector3 fPFreeze;      // 3-momentum at freezeout-time (pz, py, pz) (GeV/c)
  Int_t    fOrigin;       // Information about origin of deuterons: = 0: kinetic deuteron (from phsd.dat), = 1: potential/MST deuteron (from fort.891)
  ParticleFreeze() : fPdgId(0), fEnergy(0.), fTimeFreeze(0.), fOrigin(0) {fP.SetXYZ(0.,0.,0.); fXFreeze.SetXYZ(0.,0.,0.); fPFreeze.SetXYZ(0.,0.,0.); };
  ParticleFreeze(Int_t pdgId, TVector3 P, Float_t energy, TLorentzVector XFreeze, TVector3 PFreeze, Int_t Origin) : fPdgId(pdgId), fP(P), fEnergy(energy), fPFreeze(PFreeze), fOrigin(Origin) {fTimeFreeze = XFreeze.T(); fXFreeze.SetXYZ(XFreeze.X(), XFreeze.Y(), XFreeze.Z()); };
  ParticleFreeze(Int_t pdgId, Float_t Px, Float_t Py, Float_t Pz, Float_t energy, Float_t TimeFreeze, TVector3 XFreeze, TVector3 PFreeze, Int_t Origin) : fPdgId(pdgId), fEnergy(energy), fTimeFreeze(TimeFreeze), fXFreeze(XFreeze), fPFreeze(PFreeze), fOrigin(Origin) {fP.SetXYZ(Px, Py, Pz); };
  ParticleFreeze(Int_t pdgId, Float_t Px, Float_t Py, Float_t Pz, Float_t energy, TLorentzVector XFreeze, TVector3 PFreeze, Int_t Origin) : fPdgId(pdgId), fEnergy(energy), fPFreeze(PFreeze), fOrigin(Origin) {fP.SetXYZ(Px, Py, Pz); fTimeFreeze = XFreeze.T(); fXFreeze.SetXYZ(XFreeze.X(), XFreeze.Y(), XFreeze.Z()); };
};

class EventFreeze : public TObject  {
public:
  Int_t   fEventId;       // Event Number
  Float_t fB;             // Impact parameter (fm)
  Int_t   fNParticipants; // Number of participants
  Float_t fTime;          // Final time at which event is written out (fm/c)
  Float_t fPhi;           // Reaction plane angle
  Int_t   fNpa;           // Number of particles
  vector<ParticleFreeze> fparticles;
  void SetParameters(Int_t eventId, Float_t b, Float_t time, Float_t phi) { fEventId = eventId; fB = b; fTime = time; fPhi = phi; };
  void AddParticle(Int_t pdgId, TVector3 P, Float_t energy, TLorentzVector XFreeze, TVector3 PFreeze, Int_t Origin) { fparticles.push_back(ParticleFreeze(pdgId, P, energy, XFreeze,  PFreeze, Origin)); fNpa += 1; }; 
  void AddParticle(Int_t pdgId, Float_t Px, Float_t Py, Float_t Pz, Float_t energy, Float_t TimeFreezeCluster, TVector3 posfo_cluster, TVector3 pfo_cluster, Int_t Origin) { fparticles.push_back(ParticleFreeze(pdgId, Px, Py, Pz, energy, TimeFreezeCluster, posfo_cluster, pfo_cluster, Origin)); fNpa += 1; }; 
  void AddParticle(Int_t pdgId, Float_t Px, Float_t Py, Float_t Pz, Float_t energy, TLorentzVector XFreeze, TVector3 PFreeze, Int_t Origin) {
    fparticles.push_back(ParticleFreeze(pdgId, Px, Py, Pz, energy, XFreeze, PFreeze, Origin)); fNpa += 1; };
  void Clear() { fparticles.clear(); fNpa = 0; };
  ClassDef(EventFreeze, 1);
};

void read_events_freeze() {

  TString indir = "";
  TString outdir = indir;
  TString prefix = "";
  TString inFilename = Form("%s/%sphqmd_freeze.root",indir.Data(),prefix.Data());
  TString outFilename = Form("%s/%sout_freeze.root",outdir.Data(),prefix.Data());
  
  std::cout << "inFilename = " <<inFilename<< std::endl;
  std::cout << "outFilename = " <<outFilename<< std::endl;
 
  TFile *outFile;
  outFile = new TFile(outFilename,"RECREATE");

  TString histname = "timefreeze";
  TH1D *htimefreeze = new TH1D(histname,histname,150,0,150);
  
  TH1D *hp[3]; TH1D *hpfreeze[3]; TH1D *hxfreeze[3];
  TString xyz [] = {"x", "y", "z"};

  for (int i = 0; i < 3; i++) {

    histname=Form("p%s", xyz[i].Data());
    hp[i] = new TH1D(histname,histname,1401,-7.005,7.005);
    hp[i] ->Sumw2();
    
    histname=Form("p%s_freeze", xyz[i].Data());
    hpfreeze[i] = new TH1D(histname,histname,1401,-7.005,7.005);
    hpfreeze[i] ->Sumw2();

    histname=Form("%s_freeze", xyz[i].Data());
    hxfreeze[i] = new TH1D(histname,histname,4001,-200.05,200.05);
    hxfreeze[i] ->Sumw2();
  }

  Float_t binwidth_p = hp[0] -> GetXaxis() -> GetBinWidth(1);
  Float_t binwidth_x = hxfreeze[0] -> GetXaxis() -> GetBinWidth(1);

  TFile *inFile;
  inFile = new TFile(inFilename);

  TTree *tree;
  tree = (TTree *)inFile->Get("events");
  EventFreeze *event = new EventFreeze();
  tree->SetBranchAddress("event", &event);
  
  Int_t nevents = tree->GetEntries() ;
  cout<<<<nevents<<" nevents"endl;
   
  for (int ievent=0; ievent < tree->GetEntries(); ievent++) {
  
    inFile->cd();
    tree->GetEntry(ievent);

    // Get variables for events
    Int_t EventId = event->fEventId;
    Float_t b = event->fB;
    Float_t finaltime = event->fTime;
    Float_t phi = event->fPhi;
    Int_t nparticles = event->fNpa;
    
    for (auto particle : event->fparticles) {

      // Get variables for particles
      Int_t pdgId = particle.fPdgId;
      Float_t Px = particle.fP.X();
      Float_t Py = particle.fP.Y();
      Float_t Pz = particle.fP.Z();
      TVector3 mom = particle.fP;
      Float_t energy = particle.fEnergy;
      
      // For deuterons
      Int_t origin = particle.fOrigin; // = 0: kinetic deuteron; = 1: potential MST deuteron
      
      // Freezeout coordinates
      Float_t timeFreeze = particle.fTimeFreeze;
      TVector3 posFreeze = particle.fXFreeze;
      Float_t XFreeze = particle.fXFreeze.X();
      Float_t YFreeze = particle.fXFreeze.Y();
      Float_t ZFreeze = particle.fXFreeze.Z();
      TVector3 momFreeze = particle.fPFreeze;
      Float_t pXFreeze = particle.fPFreeze.X();
      Float_t pYFreeze = particle.fPFreeze.Y();
      Float_t pZFreeze = particle.fPFreeze.Z();

      // Make example histograms
      outFile->cd();

      Float_t rapidity =0.5*TMath::Log((energy+Pz)/(energy-Pz));

      htimefreeze->Fill(particle.fTimeFreeze,1.0/nevents);
      for (int i = 0; i < 3; i++) {
	hp[i]->Fill(mom(i),1.0/nevents/binwidth_p);
	hpfreeze[i]->Fill(momFreeze(i),1.0/nevents/binwidth_p);
	hxfreeze[i]->Fill(posFreeze(i),1.0/nevents/binwidth_x);
      }
    }
  }	

  // Plot example histograms
  outFile->cd();

  TCanvas* Canvas;
  TLegend *legend; 
  TLegend *legend2; 
  TString name[] = {"pX_final_freeze", "pY_final_freeze", "pZ_final_freeze"};
  TString nameX[] = {"#it{p}_{X} (GeV/c)", "#it{p}_{Y} (GeV/c)", "#it{p}_{Z} (GeV/c)"};
  TString nameY[] = {"#it{dN/dp}_{X} (GeV/c)^{-1}", "#it{dN/dp}_{Y} (GeV/c)^{-1}", "#it{dN/dp}_{Z} (GeV/c)^{-1}"};
  Float_t xmax[] = {2.0,2.0,6.5};
  Float_t ymax[] = {1200.,1200.0,400.};

  // Momentum
  for (int ip = 0; ip < 3; ip++) {
    
    Canvas = new TCanvas(name[ip],name[ip],500,500);
    legend = new TLegend(0.65,0.8,0.85,0.95); 
    legend2 = new TLegend(0.2,0.83,0.4,0.95);

    Canvas->SetBottomMargin(0.09);
    Canvas->SetLeftMargin(0.12);
    Canvas->SetTopMargin(0.04);
    Canvas->SetRightMargin(0.04);
   
    Canvas->cd(1);
    hp[ip]->SetStats(0);
    hp[ip]->GetYaxis()->SetTitleOffset(1.7);
    hp[ip]->GetXaxis()->SetTitle(nameX[ip]);
    hp[ip]->GetYaxis()->SetTitle(nameY[ip]);
  
    hp[ip]->SetTitle(0);   
    hp[ip]->GetXaxis()->SetRangeUser(-xmax[ip],xmax[ip]);
    hp[ip]->GetYaxis()->SetRangeUser(0,ymax[ip]);

    hp[ip]->SetLineColor(kSpring-6);
    hp[ip]->SetLineWidth(4);
    hp[ip]->Draw("hist c");
    legend->AddEntry(hp[0],"at final time","l");
  
    hpfreeze[ip]->SetLineColor(kBlue-7);
    hpfreeze[ip]->SetLineWidth(4);
    hpfreeze[ip]->SetLineStyle(2);
    hpfreeze[ip]->Draw("hist c same");
    legend->AddEntry(hpfreeze[ip],"at freezeout time","l");

    legend->SetTextSize(0.03);
    legend->SetFillColor(0);
    legend->SetBorderSize(0);
    legend->Draw();
    
    Canvas->SaveAs(Form("%s/%s.pdf",outdir.Data(),name[ip].Data()));
  }

  // Freeze-out position

  name[0] = "X_freeze"; name[1] = "Y_freeze"; name[2] =  "Z_freeze";
  nameX[0] = "#it{X} (fm)"; nameX[1] = "#it{Y} (fm)"; nameX[2] = "#it{Z} (fm)";
  nameY[0] = "#it{dN/dX} (fm)^{-1}"; nameY[1] = "#it{dN/dY} (fm)^{-1}"; nameY[2] = "#it{dN/dZ} (fm)^{-1}";
  xmax[0] = 30.0; xmax[1] = 30.0; xmax[2] = 50.0;
  ymax[0] = 50.; ymax[1] = 70.; ymax[2] = 50.; 
  
  for (int ix = 0; ix < 3; ix++) {
    
    Canvas = new TCanvas(name[ix],name[ix],500,500);
    legend = new TLegend(0.65,0.8,0.85,0.95); 
    legend2 = new TLegend(0.2,0.83,0.4,0.95);

    Canvas->SetBottomMargin(0.09);
    Canvas->SetLeftMargin(0.12);
    Canvas->SetTopMargin(0.04);
    Canvas->SetRightMargin(0.04);
   
    Canvas->cd(1);
    hxfreeze[ix]->SetStats(0);
    hxfreeze[ix]->GetYaxis()->SetTitleOffset(1.7);
    hxfreeze[ix]->GetXaxis()->SetTitle(nameX[ix]);
    hxfreeze[ix]->GetYaxis()->SetTitle(nameY[ix]);
  
    hxfreeze[ix]->SetTitle(0);   
    hxfreeze[ix]->GetXaxis()->SetRangeUser(-xmax[ix],xmax[ix]);
    hxfreeze[ix]->GetYaxis()->SetRangeUser(0,ymax[ix]);

    hxfreeze[ix]->SetLineColor(kSpring-6);
    hxfreeze[ix]->SetLineWidth(4);
    hxfreeze[ix]->Draw("hist c");
    legend->AddEntry(hxfreeze[ix],"at freezeout time","l");

    legend->SetTextSize(0.03);
    legend->SetFillColor(0);
    legend->SetBorderSize(0);
    legend->Draw();
    
    Canvas->SaveAs(Form("%s/%s.pdf",outdir.Data(),name[ix].Data()));
  }

  // Freeze-out time
  Canvas = new TCanvas("timefreeze","timefreeze",500,500); 
  Canvas->SetBottomMargin(0.09);
  Canvas->SetLeftMargin(0.12);
  Canvas->SetTopMargin(0.04);
  Canvas->SetRightMargin(0.04);
  legend2 = new TLegend(0.2,0.8,0.4,0.95);
   
  Canvas->cd(1);
  htimefreeze->SetStats(0);
  htimefreeze->GetYaxis()->SetTitleOffset(1.7);
  htimefreeze->GetXaxis()->SetTitle("#it{freezeout-time} (fm/c)" );
  htimefreeze->GetYaxis()->SetTitle("particles/event");  
  htimefreeze->SetTitle(0);   

  htimefreeze->SetLineColor(kSpring-6);
  htimefreeze->SetLineWidth(4);
  htimefreeze->Draw("hist");

  Canvas->SaveAs(Form("%s/timefreeze.pdf",outdir.Data()));
    
  inFile->Close();
  
  outFile->Write();
  outFile->Save(); 
  
}
