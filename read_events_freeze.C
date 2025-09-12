//************************************************************************************************************************************************
/** author: Susanne Glaessel (Universitaet Frankfurt)
 ** Example macro for reading EventFreeze made from PHQMD-output with convert_phqmd_unigen_freezeout.C.
 **/

#include "TROOT.h"
#include "TFile.h"
#include "TMath.h"
#include "TString.h"
#include <iostream>

class ParticleFreeze : public TObject {
private:
  Int_t fIndex;
  Int_t fPdgId;               // PDG code
  Int_t fParent;              // Index of parent (only available for phi and K* with current PHQMD version)
  Int_t fDecay;               // decay index = type of process / reaction from which particle comes
  TVector3 fP;                // 3-momentum at final time (pz, py, pz) (GeV/c)
  Float_t fEnergy;            // Energy (GeV/c)
  Float_t fTimeFreeze;        // Freezeout-time (fm/c)
  TVector3 fXFreeze;          // Position at freezeout-time (fm)
  TVector3 fPFreeze;          // 3-momentum at freezeout-time (pz, py, pz) (GeV/c)
  Float_t  fEnergyFreeze;     // energy at freezeout-time (GeV)
  Int_t fOrigin;              // Information about origin of deuterons: = 0: kinetic deuteron (from phsd.dat), = 1: potential/MST deuteron (from fort.891)
  Int_t fWeight;              // weight

public:
  inline Int_t   GetIndex()       const {return fIndex;}
  inline Int_t   GetPdg()         const {return fPdgId;}
  inline Int_t   GetParent()      const {return fParent;}
  inline Int_t   GetDecay()       const {return fDecay;}
  inline Float_t Px()             const {return fP.X();}
  inline Float_t Py()             const {return fP.Y();}
  inline Float_t Pz()             const {return fP.Z();}
  inline Float_t E()              const {return fEnergy;}
  inline TLorentzVector GetMomentum() const {return TLorentzVector(fP.X(),fP.Y(),fP.Z(),fEnergy);}
  inline Float_t XFreeze()        const {return fXFreeze.X();}
  inline Float_t YFreeze()        const {return fXFreeze.Y();}
  inline Float_t ZFreeze()        const {return fXFreeze.Z();}
  inline Float_t TFreeze()        const {return fTimeFreeze;}
  Float_t EFreeze();       
  inline TLorentzVector GetPositionFreeze() const {return TLorentzVector(fXFreeze.X(),fXFreeze.Y(),fXFreeze.Z(),fTimeFreeze);}
  inline Float_t PxFreeze()       const {return fPFreeze.X();}
  inline Float_t PyFreeze()       const {return fPFreeze.Y();}
  inline Float_t PzFreeze()       const {return fPFreeze.Z();}
  inline Float_t EnergyFreeze()   const {return fEnergyFreeze;}   
  inline TLorentzVector GetMomentumFreeze() const {return TLorentzVector(fPFreeze.X(),fPFreeze.Y(),fPFreeze.Z(),fEnergyFreeze);}
  inline Int_t   GetOrigin()      const {return fOrigin;}
  inline Int_t   GetWeight()      const {return fWeight;}
  
  ParticleFreeze() : fIndex(-1), fPdgId(0), fParent(-1), fDecay(-1), fEnergy(0.), fTimeFreeze(0.), fOrigin(0), fWeight(0) {fP.SetXYZ(0.,0.,0.); fXFreeze.SetXYZ(0.,0.,0.); fPFreeze.SetXYZ(0.,0.,0.); };

  virtual ~ParticleFreeze() = default;

  ClassDef(ParticleFreeze, 1);
};

class EventFreeze : public TObject  {
public:
  Int_t fEventId;
  Float_t fB;
  Int_t fNParticipants;
  Float_t fTime;
  Float_t fPhi;
  Int_t fNpa;
  vector<ParticleFreeze> fParticles;

  inline Int_t    GetEventId()       const {return fEventId;}
  inline Float_t  GetB()             const {return fB;}
  inline Int_t    GetNParticipants() const {return fNParticipants;}
  inline Int_t    GetTime()          const {return fTime;}
  inline Float_t  GetPhi()           const {return fPhi;}
  inline Int_t    GetNpa()           const {return fNpa;}
  inline std::vector<ParticleFreeze> GetParticleList() const {return fParticles;}
  ParticleFreeze GetParticle(Int_t index) const;

  void Clear() { fParticles.clear(); fNpa = 0; };
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
  cout<<nevents<<" nevents"<<endl;

  for (int ievent=0; ievent < tree->GetEntries(); ievent++) {
  
    inFile->cd();
    tree->GetEntry(ievent);

    // Get variables for events
    Int_t EventId = event->GetEventId();
    Float_t b = event->GetB();
    Float_t finaltime = event->GetTime();
    Float_t phi = event->GetPhi();
    Int_t nparticles = event->GetNpa();
    
    for (auto particle : event->GetParticleList()) {

      // Get variables for particles
      
      Int_t pdgId = particle.GetPdg();
      Float_t Px = particle.Px();
      Float_t Py = particle.Py();
      Float_t Pz = particle.Pz();
      TLorentzVector mom = particle.GetMomentum();
      Float_t energy = particle.E();

      // For deuterons
      Int_t origin = particle.GetOrigin(); // = 0: kinetic deuteron; = 1: potential MST deuteron
      
      // Freezeout coordinates
      Float_t timeFreeze = particle.TFreeze();
      TLorentzVector posFreeze = particle.GetPositionFreeze();
      Float_t XFreeze = particle.XFreeze();
      Float_t YFreeze = particle.YFreeze();
      Float_t ZFreeze = particle.ZFreeze();
      TLorentzVector momFreeze = particle.GetMomentumFreeze();
      Float_t pXFreeze = particle.PxFreeze();
      Float_t pYFreeze = particle.PyFreeze();
      Float_t pZFreeze = particle.PzFreeze();

      // Make example histograms
      outFile->cd();

      Float_t rapidity =0.5*TMath::Log((energy+Pz)/(energy-Pz));

      htimefreeze->Fill(timeFreeze,1.0/nevents);
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
