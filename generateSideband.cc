//
//
// 
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <dirent.h>
#include <sys/stat.h>
#include "Riostream.h"
#include <map>
#include <string>
//#include <boost/algorithm/string.hpp>
//#include <boost/algorithm/string/trim.hpp>
#include <vector>
#include <math.h>
//#include <TCint.h>
//#include <TGenericClassInfo.h> 
#include <TMath.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TF1.h>
#include <TF2.h>
#include <TF3.h>
#include <TROOT.h>
#include <TEnv.h>
#include <TSystem.h>
#include <TTree.h>
#include "TBranch.h"
#include <TApplication.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TStyle.h> 
#include <TF1.h>  
#include <TF2.h> 
#include <TLorentzVector.h>
#include <TVector3.h>
#include "TChain.h"
#include <time.h> 
#include <TSystemDirectory.h>
#include <TEfficiency.h>
#include <TLegend.h>
#include <TMinuit.h>
#include "Math/WrappedMultiTF1.h"
#include "TRandom.h" 
#include "TRandom3.h" 
#include  <TStopwatch.h>
#include "TH1F.h"
#include "TH2F.h"			// unused?
#include "TStyle.h"
#include "TCanvas.h"
#include <TGraphAsymmErrors.h>
#include <TFrame.h>
#include <TFitResult.h>
#include <TFitter.h>
#include "Fit/Fitter.h"
#include <TMatrixDSym.h>
#include <TMatrixD.h>
#include <TBinomialEfficiencyFitter.h>
#include <TKDTreeBinning.h>
#include <TH2Poly.h>
//#if !defined(__CINT__) || defined(__MAKECINT__)
#include <RooRandom.h>
#include <RooFit.h>
#include <RooMinuit.h>
#include <RooRealVar.h>
#include <RooConstVar.h>
#include <RooPlot.h>
#include <RooDataSet.h>
#include <RooDataHist.h>
#include <RooHistPdf.h>
#include <RooRealVar.h>
#include <RooAbsReal.h>
#include <RooArgList.h>
#include <RooRealProxy.h>
#include <RooCategoryProxy.h>
#include <RooAbsCategory.h>
#include <RooProdPdf.h>
#include <RooAddPdf.h>
#include <RooProduct.h>
#include <RooGaussian.h>
#include <RooCBShape.h>
#include <RooPolynomial.h>
#include <RooChebychev.h>
#include <RooWorkspace.h>
#include <RooExponential.h>
#include <RooErrorVar.h>
#include <RooFitResult.h>
#include <RooRangeBinning.h>
#include <RooBinning.h>
#include <RooNumGenConfig.h>
#include <RooBernstein.h>
#include <RooPolynomial.h>
#include <RooExtendPdf.h>
#include <RooMultiVarGaussian.h>
#include "RooBernsteinSideband.h"
#include "RooMCStudy.h"
#include "TRatioPlot.h"

#include <sys/time.h>
#include <sys/times.h>
#include <iostream>
// System stuff
#include <fstream> 
#include <sys/time.h>
#include <sys/times.h>

timeval startTime, stopTime, totalTime;
clock_t startCPU, stopCPU; 
tms startProc, stopProc; 

using namespace std; 
using namespace ROOT;
using namespace RooFit;
void ReadSBModel();
std::string OutSaveFileName =  "savesb_2016_b0.root";
std::string OutSaveDir =  "GenBin";
char *ReadOutSaveFileName;

//std::string PNGNameReadSB3D = "generateSideband";
RooDataSet *geneDataNew=0;
RooDataSet *geneMassNew=0;
RooBernsteinSideband *bkg_ang_pdf=0;
RooAbsPdf *bkg_mass_sb=0;
//RooExponential *bkg_mass_sb=0;
RooRealVar *max_sbl=0;
RooRealVar *min_sbr=0;
TMatrixD *covMatrix =0;

int    xCosLHBin =  25;
int    xCosKHBin =  25;
int    xPhiHBin  =  25;
double XMinCosThetaL	     = -1.;
double XMaxCosThetaL	     =  1.;
double XMinCosThetaK	     = -1.;
double XMaxCosThetaK	     =  1.;
double XMinPhi  	     =-TMath::Pi();
double XMaxPhi  	     = TMath::Pi();
double XMinSign = 5.0;
double XMaxSign = 5.6;
double XStepSign = 0.025;
float  xMassHBin = (XMaxSign -XMinSign)/XStepSign;
//
int MaxGenerate = 100;

int Q2Bin  =0 ;
int RunEra =0 ;

int GenEntries=1000;
//
int num_threads=20;
//
//=================================================
int main (int argc, char** argv) {
  TStopwatch TimeWatch;
  TimeWatch.Start();
  
  if (argc>1 ){
     ReadOutSaveFileName  = (char *) malloc(strlen(argv[1])+1);
     strcpy(ReadOutSaveFileName,argv[1]);
     OutSaveFileName=Form("%s",ReadOutSaveFileName);
     std::cout<<"========================================================================="<<endl;
     std::cout<<Form("Setting the File: %s",argv[1])<<std::endl;
     std::cout<<"========================================================================="<<endl;
  }else{
     std::cout<<"========================================================================="<<endl;
     std::cout<<Form("Setting DEFAULT File: %s",OutSaveFileName.c_str())<<std::endl;
     std::cout<<"========================================================================="<<endl;
  }
//  PNGNameReadSB3D=Form("%s-%s.png",PNGNameReadSB3D.c_str(),OutSaveFileName);
  

  ReadSBModel(); 
//  app.Run() ;
  cout<<"esco..." <<endl;
  TimeWatch.Stop();
  TimeWatch.Print();
  
  return 0 ;
}
void  ReadSBModel(){
    ROOT::EnableImplicitMT(num_threads);
    RooRealVar* ctL = new RooRealVar("ctL", "ctL",  XMinCosThetaK,XMaxCosThetaK);
    RooRealVar* ctK = new RooRealVar("ctK", "ctK",  XMinCosThetaL,XMaxCosThetaL);
    RooRealVar* phi = new RooRealVar("phi", "phi",  XMinPhi,XMaxPhi);
//    RooRealVar* tagged_mass = new RooRealVar("tagged_mass" , "#mu^{+}#mu^{-}K#pi mass", XMinSign,XMaxSign, "GeV/c^2");
    
    TFile* OutSaveFile = new TFile( OutSaveFileName.c_str(), "READ" );
    if ( !OutSaveFile || !OutSaveFile ->IsOpen() ) {
      cout<<Form("File not found: %s\n",OutSaveFileName.c_str())<<endl;
      exit(1);
    }
    RooWorkspace* w = (RooWorkspace*)OutSaveFile->Get("wsb");
    if ( !w || w->IsZombie() ) {
     cout<<Form("Workspace not found in file:%s\n",OutSaveFileName.c_str())<<endl;
     exit(1);
    } else {
     cout<<Form("Workspace Found!!! In file : %s\n",OutSaveFileName.c_str())<<endl;
    }
//    w->Print();
    
//
//  Cov Matrix
//   
//     TMatrixD *covMatrix = (TMatrixD*)w->obj(Form("covMatrix_bin%d_%d",Q2Bin,RunEra));
//     if(covMatrix==0){
//         cout<<Form("covMatrix_bin%d_%d not found in file:%s\n",OutSaveFileName.c_str())<<endl;
//     }else{  
//      cout<<Form("covMatrix found in file:%s  DUMP:\n",OutSaveFileName.c_str())<<endl;
//      covMatrix->Print("f=  %10.3e  ");
//     } 
    bool bfound=false;
    for(int i=0;i<8;i++){
     if (bfound) break;
     for(int j=6;j<9;j++){
      if (bfound) break;
      covMatrix = (TMatrixD*)w->obj(Form("covMatrix_bin%d_201%d",i,j));
      if ( !covMatrix) {
       if(i==7&&j==8){
//         cout<<Form("covMatrix_bin%d_201%d not found in file:%s\n",i,j,OutSaveFileName.c_str())<<endl;
//         exit(1);
       }	 
      } else {
       cout<<Form("covMatrix_bin%d_201%d Found!! in file:%s DUMP:\n",i,j,OutSaveFileName.c_str())<<endl;
       bfound=true;
      }
     }  
    }
    if(bfound) {
     covMatrix->Print("f=  %10.3e  ");
    }else{
     cout<<Form("NO covMatrix found in file:%s !!!\n",OutSaveFileName.c_str())<<endl;
    } 
//  
    bfound=false;
    for(int i=0;i<8;i++){
     if (bfound) break;
     for(int j=6;j<9;j++){
      if (bfound) break;
      bkg_ang_pdf= (RooBernsteinSideband *)w->pdf(Form("BernSideBand_bin%d_201%d",i,j));
      if ( !bkg_ang_pdf ) {
//       cout<<Form("BernSideBand_201%d not found\n",i)<<endl;
       if(i==7&&j==8){
         cout<<Form("Not found any BernSideBand function in  in file:%s\n",OutSaveFileName.c_str())<<endl;
//         cout<<Form("BernSideBand_bin%d_201%d not found in file:%s\n",i,j,OutSaveFileName.c_str())<<endl;
         exit(1);
       }	 
      } else {
       cout<<Form("BernSideBand_bin%d_201%d Found!! in file:%s\n",i,j,OutSaveFileName.c_str())<<endl;
       Q2Bin=i;
       RunEra=2010+j;
       bfound=true;
      }
     }  
    }
    bfound=false;
    for(int i=0;i<=7;i++){
     if (bfound) break;
     for(int j=6;j<=8;j++){
      if (bfound) break;
      bkg_mass_sb = (RooAbsPdf *)w->pdf(Form("bkg_mass_sb_bin%d_201%d",i,j));
//      bkg_mass_sb = (RooExponential *)w->pdf(Form("bkg_mass_sb_bin%d_201%d",i,j));
      if ( ! bkg_mass_sb) {
       if(i==7&&j==8){
         cout<<Form("bkg_mass_sb_bin*_201* not found in file:%s\n",OutSaveFileName.c_str())<<endl;
         exit(1);
       }
      } else {
       cout<<Form("bkg_mass_sb_bin%d_201%d Found!! in file:%s\n",i,j,OutSaveFileName.c_str())<<endl;
       bfound=true;
      } 
     }
    }
    bfound=false;
    for(int i=0;i<=7;i++){
     if (bfound) break;
     for(int j=6;j<=8;j++){
      if (bfound) break;
      max_sbl = (RooRealVar *)w->var(Form("max_sbl_bin%d_201%d",i,j));
      if ( ! max_sbl) {
       if(i==7&&j==8){
         cout<<Form("max_sbl_bin*_201* not found in file:%s\n",OutSaveFileName.c_str())<<endl;
         exit(1);
       }
      } else {
       cout<<Form("max_sbl_bin%d_201%d=%f\n",i,j,max_sbl->getVal())<<endl;
       bfound=true;
      } 
     }
    }
    bfound=false;
    for(int i=0;i<=7;i++){
     if (bfound) break;
     for(int j=6;j<=8;j++){
      if (bfound) break;
      min_sbr = (RooRealVar *)w->var(Form("min_sbr_bin%d_201%d",i,j));
      if ( ! min_sbr ) {
       if(i==7&&j==8){
         cout<<Form("min_sbr_bin*_201* not found in file:%s\n",OutSaveFileName.c_str())<<endl;
         exit(1);
       }
      } else {
       cout<<Form("min_sbr_bin%d_201%d=%f\n",i,j,min_sbr->getVal())<<endl;
       bfound=true;
      } 
     }
    }
//  if doesn't exite the output directory, create it.

    OutSaveDir=Form("%s%d_%d",OutSaveDir.c_str(),Q2Bin,RunEra);
    DIR* dir = opendir(OutSaveDir.c_str());
 
    if(!dir) {
     gSystem->Exec(Form("mkdir -p %s",OutSaveDir.c_str()));
     std::cout<<Form("Create Dir  %s",OutSaveDir.c_str())<<std::endl;
    }

    int  maxDegree1 =   bkg_ang_pdf->maxDegree1Get();
    int  maxDegree2 =   bkg_ang_pdf->maxDegree2Get();
    int  maxDegree3 =   bkg_ang_pdf->maxDegree3Get();
    int  numParameters = (maxDegree1+1)*(maxDegree2+1)*(maxDegree3+1);

    cout<<Form("maxDegree1=%d\n",maxDegree1); 
    cout<<Form("maxDegree2=%d\n",maxDegree2); 
    cout<<Form("maxDegree3=%d\n",maxDegree3); 
    cout<<Form("numParameters=%d\n",numParameters); 

    RooArgList observables (*ctK, *ctL, *phi);
    std::vector<TH1D*> HPar;
    std::vector<TH3D*> HBern;
    std::vector<TH1D*> HBernX;
    std::vector<TH1D*> HBernY;
    std::vector<TH1D*> HBernZ;
    std::vector<RooBernsteinSideband*> BernSideBand; 
    RooArgList* bkg_ang_params = (RooArgList*)bkg_ang_pdf->getParameters(observables);
//    auto iter = bkg_ang_params->createIterator();
//    RooRealVar* ivar =  (RooRealVar*)iter->Next();
    RooArgList* bkg_ang_params_rnd = new RooArgList("bkg_ang_params_rnd");
    RooArgList* bkg_ang_params_redu_tmp = new RooArgList("bkg_ang_params_redu_tmp");
    RooArgList* bkg_ang_params_rnd_redu_tmp = new RooArgList("bkg_ang_params_rnd_redu_tmp");
    RooRealVar * rndPar;
    RooRealVar * newPar;
    int i=0;
    cout<<"=============== Original parameters ==============="<<endl;
//    while (ivar) {
    for (int i=0;i<numParameters;++i){
      RooRealVar* ivar =  (RooRealVar* )bkg_ang_params->find(Form("p%03d_%d",i,RunEra)) ;
      if( !ivar) {
       cout<<"Error: "<<Form("p%03d",i)<<" not found !!!"<<endl;
       exit(1);
      }
//      bkg_ang_params_rnd->add( *ivar ) ;
      newPar = new RooRealVar(Form("%s_%d",ivar->GetName(),RunEra),ivar->GetTitle(),ivar->getVal(),ivar->getMin(),ivar->getMax());
      bkg_ang_params_rnd->add( *newPar ); 
//      ivar->setConstant(true);
      if( ivar->getVal()>0.0&&ivar->getVal()!=1.) {
       cout <<Form("%s=%f this can be randomized",ivar->GetName(),ivar->getVal())<<endl;
       bkg_ang_params_redu_tmp->add( *ivar );
       rndPar =new RooRealVar(Form("rnd_%s",ivar->GetName()),ivar->GetName(),0.,ivar->getMin(),ivar->getMax());
//       rndPar =(RooRealVar*)ivar->clone(Form("rnd_%s",ivar->GetName()));
       bkg_ang_params_rnd_redu_tmp->add( *rndPar );
       HPar.push_back(new TH1D(Form("HPar_%s",ivar->GetName()),ivar->GetName(),200,-5.,5.));
      }else{
       cout <<Form("%s=%f",ivar->GetName(),ivar->getVal())<<endl;
      } 
//      ivar = (RooRealVar*) iter->Next();
    }
//
    if(!covMatrix) {
     cout<<Form("covMatrix not found: Exit!!!")<<endl;
     exit(0);
    }else{
     cout<<Form("covMatrix Nrows=%d ",covMatrix->GetNrows())<<endl;
    }
    TMatrixDSym *covMatrixSym=(TMatrixDSym *)covMatrix;     
    covMatrixSym->Print("f=  %10.3e  ");
// 
    TFile* OutSaveGenerate = new TFile( Form("%s/generate-savesb_%d_b%d.root",OutSaveDir.c_str(),RunEra,Q2Bin), "RECREATE" );
// 	 covMatrixSym->Use(covMatrix->GetNrows(),covMatrix->GetMatrixArray());
    RooMultiVarGaussian SideBandStatGauss("SideBandStatGauss", "SideBandStatGauss", *bkg_ang_params_rnd_redu_tmp, *bkg_ang_params_redu_tmp, *covMatrixSym);
    RooDataSet* genParSB = SideBandStatGauss.generate(*bkg_ang_params_rnd_redu_tmp,MaxGenerate);
//    RooAbsArg *aptr = 0; 
    RooRealVar *aptr = 0;
//    RooBernsteinSideband *BernSideBand=0; 
    TH3D* HBern_ref  = (TH3D*) bkg_ang_pdf->createHistogram(Form("%s_ref",bkg_ang_pdf->GetName()), *ctL, Binning(50),YVar(*ctK,Binning(50)),ZVar(*phi,Binning(50)));
     
    TH1D*  HBernX_ref= (TH1D*) HBern_ref->ProjectionX("HBernX_ref",1,HBern_ref->GetNbinsY(),1,HBern_ref->GetNbinsZ());
    TH1D*  HBernY_ref= (TH1D*) HBern_ref->ProjectionY("HBernY_ref",1,HBern_ref->GetNbinsX(),1,HBern_ref->GetNbinsZ());
    TH1D*  HBernZ_ref= (TH1D*) HBern_ref->ProjectionZ("HBernZ_ref",1,HBern_ref->GetNbinsX(),1,HBern_ref->GetNbinsY());
//
    float MaxHBernX=HBernX_ref->GetMaximum()*1.2;
    float MaxHBernY=HBernY_ref->GetMaximum()*1.2;
    float MaxHBernZ=HBernZ_ref->GetMaximum()*1.2;
//    bkg_ang_pdf->Delete();
//    RooArgList* bkg_ang_params_rnd_era = new RooArgList("bkg_ang_params_rnd_era");
    for (int i=0 ; i<genParSB->sumEntries() ; i++) {
     const RooArgList* bkg_ang_params_rnd_redu = new RooArgList("bkg_ang_params_rnd_redu");
     bkg_ang_params_rnd_redu =  (const RooArgList*)genParSB->get(i);
     auto iter_rnd_redu  = bkg_ang_params_rnd_redu->createIterator();
     RooRealVar* irnd =  (RooRealVar*)iter_rnd_redu->Next();
     cout <<Form("============ Gen event N. %d",i)<<endl;
     int j=0;
     while (irnd) {
      cout <<Form("%s=%f [%s]",irnd->GetName(),irnd->getVal(),irnd->GetTitle())<<endl;
      aptr =  (RooRealVar* )bkg_ang_params_rnd->find(Form("%s_%d",irnd->GetTitle(),RunEra)) ;
      if(!aptr){
       cout<<Form("Par => %s  not found !!!",irnd->GetTitle())<<endl;
       exit(1);
      } 
      aptr->setVal(irnd->getVal());
      HPar[j]->Fill(irnd->getVal());
      irnd =  (RooRealVar*)iter_rnd_redu->Next();
      j++;
     }
     auto iter_rnd  = bkg_ang_params_rnd->createIterator();
     RooRealVar* ivar_rnd =  (RooRealVar*)iter_rnd->Next();
     cout <<Form("============ Rewrite event N. %d",i)<<endl;
     while (ivar_rnd) {
//       bkg_ang_params_rnd_era->add( new RooRealVar(Form("%s_%d",ivar_rnd->GetName(),RunEra),ivar_rnd->GetTitle(),ivar_rnd->getVal(),ivar_rnd->getMinimum(),ivar_rnd->getMaximum()); 
//       bkg_ang_params_rnd_era->add( *(RooRealVar*)ivar_rnd->Clone(Form("%s_%d",ivar_rnd->GetName(),RunEra)));
       cout <<Form("%s=%f ",ivar_rnd->GetName(),ivar_rnd ->getVal())<<endl;
       ivar_rnd =  (RooRealVar*)iter_rnd->Next();
     } 
     const char* BernName = Form("%s",bkg_ang_pdf->GetName());
//     const char* BernName = Form("%s_%d",bkg_ang_pdf->GetName(),i);
     const char* BernTitle= Form("%s_%d",bkg_ang_pdf->GetTitle(),i);
     BernSideBand.push_back(new RooBernsteinSideband(BernName,BernTitle,*ctL,*ctK,*phi,*bkg_ang_params_rnd,maxDegree1,maxDegree2,maxDegree3));
     HBern.push_back( (TH3D*) BernSideBand[i]->createHistogram(Form("%s_%d",BernName,i), *ctL, Binning(50),YVar(*ctK,Binning(50)),ZVar(*phi,Binning(50))) );
     RooWorkspace wsb("wsb","workspace sideband");
     wsb.import(*covMatrix,Form("covMatrix_bin%d_%d",Q2Bin,RunEra));
     wsb.import(*BernSideBand[i]);
     wsb.import(*bkg_mass_sb);
     wsb.import(*max_sbl);
     wsb.import(*min_sbr);
     wsb.writeToFile(Form("%s/savesb_%d_b%d.root.Gen%d",OutSaveDir.c_str(),RunEra,Q2Bin,i));
     cout<<Form("================= SAVE %d =================",i)<<endl; 
     OutSaveGenerate->cd();
     
     HBernX.push_back((TH1D*) HBern[i]->ProjectionX(Form("HBernX_%d",i),1,HBern[i]->GetNbinsY(),1,HBern[i]->GetNbinsZ()));
     HBernY.push_back((TH1D*) HBern[i]->ProjectionY(Form("HBernY_%d",i),1,HBern[i]->GetNbinsX(),1,HBern[i]->GetNbinsZ()));
     HBernZ.push_back((TH1D*) HBern[i]->ProjectionZ(Form("HBernZ_%d",i),1,HBern[i]->GetNbinsX(),1,HBern[i]->GetNbinsY()));
     if(MaxHBernX<HBernX[i]->GetMaximum()) MaxHBernX = HBernX[i]->GetMaximum();
     if(MaxHBernY<HBernY[i]->GetMaximum()) MaxHBernY = HBernY[i]->GetMaximum();
     if(MaxHBernZ<HBernZ[i]->GetMaximum()) MaxHBernZ = HBernZ[i]->GetMaximum();
    }

    TCanvas* cProj = new TCanvas("cProj","Projections",200,200,1600,800);
    cProj->Divide(2,2);
    i=0;
    for (std::vector<TH1D*>::iterator it = HPar.begin() ; it != HPar.end(); ++it){
     
     OutSaveGenerate->cd();
     HPar[i]->Write();
     i++;
    }
    i=0;
    OutSaveGenerate->cd();
    HBernX_ref->Write();
    HBernY_ref->Write();
    HBernZ_ref->Write();
//    for (std::vector<TH1D*>::iterator it = HBernX.begin() ; it != HBernX.end(); ++it){
    for (int i=0 ; i<genParSB->sumEntries() ; i++) {
     
     HBernX[i]->Write();
     HBernY[i]->Write();
     HBernZ[i]->Write();
     cProj->cd(1);
     HBernX[i]->SetLineColor(kRed);
     HBernX[i]->SetMinimum(0.);
     HBernX[i]->SetMaximum(MaxHBernX);
     HBernX[i]->Draw("same,HIST C");
     cProj->cd(2);
     HBernY[i]->SetLineColor(kRed);
     HBernY[i]->SetMinimum(0.);
     HBernY[i]->SetMaximum(MaxHBernY);
     HBernY[i]->Draw("same,HIST C");
     cProj->cd(3);
     HBernZ[i]->SetLineColor(kRed);
     HBernZ[i]->SetMinimum(0.);
     HBernZ[i]->SetMaximum(MaxHBernZ);
     HBernZ[i]->Draw("same,HIST C");
    }
    cProj->cd(1);
    HBernX_ref->SetLineColor(kBlue);
    HBernX_ref->SetMaximum(MaxHBernX);
    HBernX_ref->SetMinimum(0.);
    HBernX_ref->Draw("same,HIST C");
    cProj->cd(2);
    HBernY_ref->SetLineColor(kBlue);
    HBernY_ref->SetMaximum(MaxHBernY);
    HBernY_ref->SetMinimum(0.);
    HBernY_ref->Draw("same,HIST C");
    cProj->cd(3);
    HBernZ_ref->SetLineColor(kBlue);
    HBernZ_ref->SetMaximum(MaxHBernZ);
    HBernZ_ref->SetMinimum(0.);
    HBernZ_ref->Draw("same,HIST C");
     
    OutSaveGenerate->cd();
    cProj->Write();
    OutSaveGenerate->Close();
    
//     
    
        
  } 
