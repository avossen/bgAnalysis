TColor* glColorTable[30];

TLegend* getLegend(TH1D** histos,bool moveRight=false)
{
  char* legendNames[]={"Coulomb_LER","Coulomb_HER","Touschek_LER","Touschek_HER","Two photon","RBB","BHWide","BHWide LargeAngle"};

  TLegend* legend=0;
  if(moveRight)
    legend=new TLegend(0.30,0.7,0.9,0.99);
  else
    legend=new TLegend(0.10,0.8,0.9,0.99);
  legend->SetNColumns(4);
  if(moveRight)
  legend->SetNColumns(3);
  for(int i=0;i<8;i++)
    {
      legend->AddEntry(histos[i],legendNames[i],"f");
    }
  return legend;

}




//phaseIdx 0: 2
//phaseIdx 1: 3
THStack* getStack(TH2D* histo, TH1D** histos, float lower, float upper,char* name, float scaleFactor, int& nbinsX,bool flip=false, bool rebin=true)
{
  int numBgSources=8;
  int bgSourceIds[8];
  bgSourceIds[0]=2;  //coulomb_ler
  bgSourceIds[1]=3;   // coulob_her
  bgSourceIds[2]=6;    //Touschek_LER
  bgSourceIds[3]=7 ;    //Touschek_HER
  bgSourceIds[4]=8;    //two gamma
  bgSourceIds[5]=18;    //RBB
  bgSourceIds[6]=19;    //BHWide
  bgSourceIds[7]=20;    //BHWideLargeAngle

  char buffer[300];

  int iFip=0;
  if(flip==false)
    iFip=1;
  sprintf(buffer,"stack_%s_%d_%f",name,iFip,lower);
  THStack* stack=new THStack(buffer,buffer);
  int numBinsX=histo->GetNbinsX();
  int numBinsY=histo->GetNbinsY();
  if(flip==true)
    {
      cout <<"nbinsX: " << numBinsX <<" nbinsY: "<< numBinsY <<endl;
      cout <<"lower: "<< lower <<" upper: "<< upper <<endl;
    }
  int tmp=numBinsX;
  if(flip)
    {
      numBinsX=numBinsY;
      numBinsY=tmp;
      cout <<" now nbinsX: " << numBinsX <<" nbinsY: "<< numBinsY <<endl;
    }
  else
    {
      cout <<" don't flip.. " << endl;
    }



  for(int i=0;i<8;i++)
    {
      sprintf(buffer,"histo%d_%s_%d_%f",i,name,iFip,lower);
      if(flip)
	histos[i]=new TH1D(buffer,buffer,3*numBinsX,lower, upper);
      else
	histos[i]=new TH1D(buffer,buffer,numBinsX,lower, upper);

      for(int j=1;j<=numBinsX;j++)
	{
	  if(flip)
	    {
	      //subtract 0.5 so layer number is centered on integers
	      histos[i]->SetBinContent(3*j-3,histo->GetBinContent(bgSourceIds[i],j)*scaleFactor);
	      histos[i]->SetBinContent(3*j-2,histo->GetBinContent(bgSourceIds[i],j)*scaleFactor);
	    }
	  else
	    {
	      histos[i]->SetBinContent(j,histo->GetBinContent(j,bgSourceIds[i])*scaleFactor);
	    }
	}
      if(rebin)
	histos[i]->Rebin(2);
      nbinsX=histos[i]->GetNbinsX();
      histos[i]->SetFillStyle(1001);
      histos[i]->SetFillColor(glColorTable[i]->GetNumber());
      stack->Add(histos[i]);
    }

  return stack;
}

 
void drawStacks(TH2D* phase2Histo, TH2D* phase3Histo, bool isPhi,float scaleFactorPhase2, float scaleFactorPhase3)
{
  char buffer[300];
  TH1D** histosPhase2=new TH1D*[8];
  TH1D** histosPhase3=new TH1D*[8];

  int numBinsX2=0;
  int numBinsX3=0;

  float lowerBound=0.0;
  if(isPhi)
    lowerBound=-TMath::Pi();

  THStack* phase2Stack=getStack(phase2Histo,histosPhase2,lowerBound,TMath::Pi(),"Phase2",scaleFactorPhase2,numBinsX2,false);
  THStack* phase3Stack=getStack(phase3Histo,histosPhase3,lowerBound,TMath::Pi(),"Phase3",scaleFactorPhase3,numBinsX3,false);

  TCanvas c;
  phase2Stack->Draw();
  TLegend* legendPhase2=getLegend(histosPhase2);
  //  int numBinsX2=phase2Stack->GetNbinsX();
  //  int numBinsX3=phase2Stack->GetNbinsX();
  cout <<"phase 2 has " << numBinsX2 <<" phase2: "<< numBinsX3 <<endl;

  float rate2=1000.0/(50*1.8*scaleFactorPhase2);  //HZ 
  float rate3=1000.0/(50*1.8*scaleFactorPhase3);  //HZ 

  float rad2=(TMath::Pi()-lowerBound)/numBinsX2;
  float rad3=(TMath::Pi()-lowerBound)/numBinsX3;

  sprintf(buffer,"%.3f HZ/%.2f rad",rate2,rad2);
  phase2Stack->GetYaxis()->SetRangeUser(0,4000);
  phase2Stack->GetYaxis()->SetTitle(buffer);
  if(isPhi)
    phase2Stack->GetXaxis()->SetTitle("#phi [rad]");
  else
    phase2Stack->GetXaxis()->SetTitle("#theta [rad]");

  c.Update();
  phase2Stack->Draw();
  c.Update();
  legendPhase2->Draw();
  if(isPhi)
    c.SaveAs("phase2Stack_phi.png");
  else
    c.SaveAs("phase2Stack_theta.png");

  phase3Stack->Draw();
  phase3Stack->GetYaxis()->SetRangeUser(0,4000);
  sprintf(buffer,"%.3f HZ/%.2f rad",rate3,rad3);
  phase3Stack->GetXaxis()->SetTitle("#theta [rad]");
  if(isPhi)
    phase3Stack->GetXaxis()->SetTitle("#phi [rad]");

  phase3Stack->GetYaxis()->SetTitle(buffer);
  c.Update();
  TLegend* legendPhase3=getLegend(histosPhase3);
  legendPhase3->Draw();
  if(isPhi)
    c.SaveAs("phase3Stack_phi.png");
  else
    c.SaveAs("phase3Stack_theta.png");


}

void drawStacksVsLayer(TH2D* phase2Histo, TH2D* phase3Histo, float scaleFactorPhase2, float scaleFactorPhase3, bool is2D=false)
{
  char buffer[300];
  TH1D** histosPhase2=new TH1D*[8];
  TH1D** histosPhase3=new TH1D*[8];

  int numBinsX2=0;
  int numBinsX3=0;

  float lowerBound=0.0;
  float upperBound=16.0;


  THStack* phase2Stack=getStack(phase2Histo,histosPhase2,lowerBound,upperBound,buffer,scaleFactorPhase2,numBinsX2,true,false);

  if(is2D)
    sprintf(buffer,"Phase3_2DHits");
  else
    sprintf(buffer,"Phase3_1DHits");
  THStack* phase3Stack=getStack(phase3Histo,histosPhase3,lowerBound,upperBound,buffer,scaleFactorPhase3,numBinsX3,true,false);


  TCanvas c;
  phase2Stack->Draw();
  TLegend* legendPhase2=getLegend(histosPhase2,true);
  //  int numBinsX2=phase2Stack->GetNbinsX();
  //  int numBinsX3=phase2Stack->GetNbinsX();
  cout <<"phase 2 has " << numBinsX2 <<" phase2: "<< numBinsX3 <<endl;

  float rate2=1000.0/(50*1.8*scaleFactorPhase2);  //HZ 
  float rate3=1000.0/(50*1.8*scaleFactorPhase3);  //HZ 

  sprintf(buffer,"%.3f HZ/layer",rate2);
  phase2Stack->GetYaxis()->SetRangeUser(0,4000);
  phase2Stack->GetYaxis()->SetTitle(buffer);
  phase2Stack->GetXaxis()->SetTitle("layer");

  c.Update();
  phase2Stack->Draw();
  c.Update();
  legendPhase2->Draw();
  if(is2D)
    c.SaveAs("phase2Stack_layer_2DHits.png");
  else
    c.SaveAs("phase2Stack_layer_1DHits.png");

  phase3Stack->Draw();
  phase3Stack->GetYaxis()->SetRangeUser(0,4000);
  sprintf(buffer,"%.3f HZ/layer",rate3);
  phase3Stack->GetXaxis()->SetTitle("layer");
  phase3Stack->GetYaxis()->SetTitle(buffer);
  c.Update();
  TLegend* legendPhase3=getLegend(histosPhase3,true);
  legendPhase3->Draw();
  if(is2D)
    c.SaveAs("phase3Stack_layer_2DHits.png");
  else
    c.SaveAs("phase3Stack_layer_1DHits.png");
}


//changed this, so that only one phase at a time but we can give the filename and time constant
//we assume 50k events
//@time time in us
void doAnalysis(char* filename, float time)
{

  //phase 2: scale all to 1ms: (assuming we threw 50k events, this means a factor of 556/50k
  float scaleFactorPhase2=556.0/50000.0;
  ///campaing 2.1 has varying number of generated events, each file is 50us/one exception of 25, but number of files is either
  //100 or 200


  //phase3, most data is 600us
  float scaleFactorPhase3=333.0/50000.0;
  glColorTable[0]=gROOT->GetColor(kBlue);
  glColorTable[1]=gROOT->GetColor(kRed);
  glColorTable[2]=gROOT->GetColor(kYellow);
  glColorTable[3]=gROOT->GetColor(kBlack);
  glColorTable[4]=gROOT->GetColor(kMagenta);
  glColorTable[5]=gROOT->GetColor(kOrange);

  glColorTable[6]=gROOT->GetColor(kSpring);
  glColorTable[7]=gROOT->GetColor(kGray);
  glColorTable[8]=gROOT->GetColor(kCyan);
  //
  // compare phase 2/3
  // get overall rates
  // bg source vs theta/phi
 //

  //    TFile phase2("../campaign15_phase2/out15_phase2_sum_nakayama.root");
    TFile phase2("../campaign15_phase2/sumCampaign15Nakayama.root");
    TFile phase3("../campaign15_phase2/sumCampaign15Nakayama.root");
  //out15_phase2_sum_nakayama.root
    //    TFile phase3("../campaign15_phase3/Out_sum.root");

  TH2D* bgSourceVsPhi_phase2=(TH2D*)phase2.Get("bgSourceVsPhi");
  TH2D* bgSourceVsPhi_phase3=(TH2D*)phase3.Get("bgSourceVsPhi");

  TH2D* bgSourceVsTheta_phase2=(TH2D*)phase2.Get("bgSourceVsTheta");
  TH2D* bgSourceVsTheta_phase3=(TH2D*)phase3.Get("bgSourceVsTheta");

  drawStacks(bgSourceVsPhi_phase2,bgSourceVsPhi_phase3,true,scaleFactorPhase2,scaleFactorPhase3);
  drawStacks(bgSourceVsTheta_phase2,bgSourceVsTheta_phase3,false,scaleFactorPhase2,scaleFactorPhase3);

//

  //yes.... there is a misspelling in the root file
  TH2D* bgSourceVsLayer_phase2=(TH2D*)phase2.Get("bgSourcPerLayer");
  TH2D* bgSourceVsLayer_phase3=(TH2D*)phase3.Get("bgSourcPerLayer");

  drawStacksVsLayer(bgSourceVsLayer_phase2,bgSourceVsLayer_phase3,scaleFactorPhase2,scaleFactorPhase3);

  bgSourceVsLayer_phase2=(TH2D*)phase2.Get("bgSourcPerLayer2D");
  bgSourceVsLayer_phase3=(TH2D*)phase3.Get("bgSourcPerLayer2D");

  drawStacksVsLayer(bgSourceVsLayer_phase2,bgSourceVsLayer_phase3,scaleFactorPhase2,scaleFactorPhase3,true);


  ///normalized and draw the hits per channel and layer:
  TH2D* posPerChannelLayerPhase2=(TH2D*)phase2.Get("posPerChannelLayer");
  TH2D* posPerChannelLayerPhase3=(TH2D*)phase3.Get("posPerChannelLayer");
  posPerChannelLayerPhase2->GetXaxis()->SetRangeUser(0,3500);
  posPerChannelLayerPhase3->GetXaxis()->SetRangeUser(0,3500);
  posPerChannelLayerPhase2->SetStats(0);
  posPerChannelLayerPhase3->SetStats(0);

  //let's scale this so the plot is in HZ: 50k events * 1.8 us = 0.09s
  posPerChannelLayerPhase2->Scale(1.0/0.09);
  posPerChannelLayerPhase3->Scale(1.0/0.09);
  TCanvas c;
    c.SetLogz();

  posPerChannelLayerPhase2->Draw("colz");
  c.SaveAs("PosPerChannelLayerPhase2.png");
  posPerChannelLayerPhase3->Draw("colz");
  c.SaveAs("PosPerChannelLayerPhase3.png");


};
