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

 
void drawStacks(TH2D* Histo, bool isPhi,float scaleFactor, char* filenameAdd)
{
  char buffer[300];
  TH1D** histos=new TH1D*[8];

  int numBinsX2=0;
  int numBinsX3=0;

  float lowerBound=0.0;
  if(isPhi)
    lowerBound=-TMath::Pi();

  THStack* Stack=getStack(Histo,histos,lowerBound,TMath::Pi(),filenameAdd,scaleFactor,numBinsX2,false);


  TCanvas c;
  Stack->Draw();
  TLegend* legend=getLegend(histos);
    int numBinsX=Histo->GetNbinsX();
  //  int numBinsX3=phase2Stack->GetNbinsX();
    cout <<"phase 2 has " << numBinsX2 <<" phase2: "<< numBinsX3 <<endl;

  float rate=1000.0/(50*1.8*scaleFactor);  //HZ 

  float rad=(TMath::Pi()-lowerBound)/numBinsX;


  sprintf(buffer,"%.3f HZ/%.2f rad",rate,rad);
  Stack->GetYaxis()->SetRangeUser(0,4000);
  Stack->GetYaxis()->SetTitle(buffer);
  if(isPhi)
    Stack->GetXaxis()->SetTitle("#phi [rad]");
  else
    Stack->GetXaxis()->SetTitle("#theta [rad]");

  c.Update();
  Stack->Draw();
  c.Update();
  legend->Draw();
  if(isPhi)
    sprintf(buffer,"Stack_phi%s.png",filenameAdd);
  else
    sprintf(buffer,"Stack_theta%s.png",filenameAdd);
  c.SaveAs(buffer);

}

void drawStacksVsLayer(TH2D* Histo, float scaleFactor,char* filenameAdd,  bool is2D=false)
{
  char buffer[300];
  TH1D** histos=new TH1D*[8];
  TH1D** histosPhase3=new TH1D*[8];

  int numBinsX2=0;
  int numBinsX3=0;

  float lowerBound=0.0;
  float upperBound=16.0;


  THStack* Stack=getStack(Histo,histos,lowerBound,upperBound,buffer,scaleFactor,numBinsX2,true,false);

  if(is2D)
    sprintf(buffer,"Phase3_2DHits");
  else
    sprintf(buffer,"Phase3_1DHits");
  //  THStack* phase3Stack=getStack(phase3Histo,histosPhase3,lowerBound,upperBound,buffer,scaleFactorPhase3,numBinsX3,true,false);


  TCanvas c;
  Stack->Draw();
  TLegend* legend=getLegend(histos,true);
  //  int numBinsX2=Stack->GetNbinsX();
  //  int numBinsX3=Stack->GetNbinsX();
  cout <<"phase 2 has " << numBinsX2 <<" : "<< numBinsX3 <<endl;

  float rate=1000.0/(50*1.8*scaleFactor);  //HZ 

  sprintf(buffer,"%.3f HZ/layer",rate);
  Stack->GetYaxis()->SetRangeUser(0,4000);
  Stack->GetYaxis()->SetTitle(buffer);
  Stack->GetXaxis()->SetTitle("layer");

  c.Update();
  Stack->Draw();
  c.Update();
  legend->Draw();

  if(is2D)
    sprintf(buffer,"Stack_layer_2DHits%s.png",filenameAdd);
  else
    sprintf(buffer,"Stack_layer_1DHits%s.png",filenameAdd);
  c.SaveAs(buffer);


}


//changed this, so that only one phase at a time but we can give the filename and time constant
//we assume 50k events
//@time time in us
void doAnalysis(char* filename, float time, char* filenameAdd)
{
  char buffer[300];
  //phase 2: scale all to 1ms: (assuming we threw 50k events, this means a factor of 556/50k
  float scaleFactor=time/(1.8*50000.0);
  ///campaing 2.1 has varying number of generated events, each file is 50us/one exception of 25, but number of files is either
  //100 or 200

  //phase3, most data is 600us
  //  float scaleFactorPhase3=333.0/50000.0;
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

  //    TFile ("../campaign15_/out15__sum_nakayama.root");
  TFile phase(filename);
    //    TFile phase3("../campaign15_/sumCampaign15Nakayama.root");
  //out15__sum_nakayama.root
    //    TFile phase3("../campaign15_phase3/Out_sum.root");

  TH2D* bgSourceVsPhi=(TH2D*)phase.Get("bgSourceVsPhi");


  TH2D* bgSourceVsTheta=(TH2D*)phase.Get("bgSourceVsTheta");


  drawStacks(bgSourceVsPhi,true,scaleFactor,filenameAdd);
  drawStacks(bgSourceVsTheta,false,scaleFactor,filenameAdd);

//

  //yes.... there is a misspelling in the root file
  TH2D* bgSourceVsLayer=(TH2D*)phase.Get("bgSourcPerLayer");
  //  TH2D* bgSourceVsLayer_phase3=(TH2D*)phase3.Get("bgSourcPerLayer");

  drawStacksVsLayer(bgSourceVsLayer,scaleFactor, filenameAdd);

  bgSourceVsLayer=(TH2D*)phase.Get("bgSourcPerLayer2D");


  drawStacksVsLayer(bgSourceVsLayer,scaleFactor,filenameAdd,true);


  ///normalized and draw the hits per channel and layer:
  TH2D* posPerChannelLayer=(TH2D*)phase.Get("posPerChannelLayer");

  posPerChannelLayer->GetXaxis()->SetRangeUser(0,3500);
  posPerChannelLayer->SetStats(0);


  //let's scale this so the plot is in HZ: 50k events * 1.8 us = 0.09s
  posPerChannelLayer->Scale(1.0/0.09);
  TCanvas c;
  c.SetLogz();

  sprintf(buffer,"PosPerChannelLayer%s.png",filenameAdd);

  posPerChannelLayer->Draw("colz");
  c.SaveAs(buffer);



};
