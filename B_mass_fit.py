from ROOT import *
import datetime
#####The fitting of the B+ from the decay: B+ -> K+ J/psi
#####The best parameters of the fit are in the picture "Parameters_of_fit.jpg" in this directory

PDG_B_MASS = 5.279

time = datetime.datetime.now().timetuple()
time_st = str(time[3])+'_'+str(time[4])+'_'+str(time[5])+'_'

ch = TChain('mytree')
ch.Add('/eos/home-d/dshmygol/SimpleFile_Bplus_Kplus_Jpsi_2018_A.root')

nEvt = ch.GetEntries()
print("Number of events:" , nEvt)

x = RooRealVar('B_mass','B_mass',5.17, 5.40)
varset = RooArgSet(x)
dh = 0
dh = RooDataSet ('dh','Dataset', varset)

for evt in range(nEvt):
	if ch.GetEntry(evt) <= 0: break
	if evt % 100000 == 0:
	       print ("Events proceeded:" , evt)
	B_mass = ch.B_mass
	if B_mass < 5.17: continue
	if B_mass > 5.40: continue
	x.setVal(B_mass)
	dh.add(varset)
	##
S = RooRealVar('S','Signal',6 * 10**5,0,9999999)
S_mean  = RooRealVar ( "S_mean" , "mean "   , PDG_B_MASS, 5.17  , 5.40       )
S1_sigma = RooRealVar ("S1_sigma", "sigma"   , 2.7594 * 10**(-2), 0.001, 0.4   )
S2_sigma = RooRealVar ("S2_sigma", "sigma"   , 1.2833*10**(-2), 0.0001, 0.1    )
S_f     = RooRealVar ("S_f"     , "f"       , 0.54  , 0.000 , 1.0       )

pdfS1 = RooGaussian( "pdfS1"  , "gaus"    , x   , S_mean, S1_sigma)
pdfS2 = RooGaussian( "pdfS2"  , "gaus"    , x   , S_mean, S2_sigma)
pdfSig =RooAddPdf  ("pdfSig", "pdfSig", RooArgList(pdfS1, pdfS2), RooArgList(S_f))

B       = RooRealVar ( "B"      , "B"       , 6*10**5 , 1     , 9999999)
B_al    = RooRealVar ( "B_al"   , "B_{#alpha}", -0.326, -10, 1.0   )
pdfG    = RooGenericPdf ("pdfG"    , "exp(@1*@0)" , RooArgList(x, B_al))

alist1  = RooArgList (pdfSig, pdfG)
alist2 = RooArgList (S, B)
pdfSum  = RooAddPdf  ("model", "model", alist1, alist2)

S_mean.setConstant(True)
S_f.setConstant(True)
S2_sigma.setConstant(True)
rrr = pdfSum.fitTo( dh, RooFit.NumCPU(7), RooFit.PrintLevel(-1), RooFit.Save())
rrr = pdfSum.fitTo( dh, RooFit.NumCPU(7), RooFit.PrintLevel(-1), RooFit.Save())
S_mean.setConstant(False)
rrr = pdfSum.fitTo( dh, RooFit.NumCPU(7), RooFit.PrintLevel(-1), RooFit.Save())
S_f.setConstant(False)
rrr = pdfSum.fitTo( dh, RooFit.NumCPU(7), RooFit.PrintLevel(-1), RooFit.Save())
S2_sigma.setConstant(False)
rrr = pdfSum.fitTo( dh, RooFit.NumCPU(7), RooFit.Save() )
rrr = pdfSum.fitTo( dh, RooFit.NumCPU(7), RooFit.Save() )
rrr.Print()

xframe = x.frame(125)
dh.plotOn(xframe, RooFit.MarkerSize(0.6), RooFit.Name('data'))
pdfSum.plotOn(xframe, RooFit.Name('fitt'), RooFit.Range(x.getMin(), x.getMax()) )
pdfSum.plotOn(xframe,RooFit.Components('pdfG'), RooFit.LineColor(kRed), RooFit.LineWidth(3))
pdfSum.plotOn(xframe,RooFit.Components('pdfSig'), RooFit.LineColor(kGreen), RooFit.LineWidth(3))
dh.plotOn(xframe, RooFit.MarkerSize(0.6), RooFit.Name('data'))

canvas = TCanvas("canvas")
canvas.cd()
xframe.Draw()
#add fit to the plo
latex = TLatex()
latex.SetNDC()
#latex.SetTextSize(0.04)
latex.DrawText(0.5,0.75,'B = ' + str(pdfSum.getParameters(dh)[0].getValV()))
latex.DrawText(0.5,0.75,'S = ' + str(pdfSum.getParameters(dh)[2].getValV()))
latex.DrawText(0.5,0.70,'S_mean = ' + str(pdfSum.getParameters(dh)[6].getValV()))
latex.DrawText(0.5,0.65,'S1_sigma = ' + str(pdfSum.getParameters(dh)[3].getValV()))
latex.DrawText(0.5,0.60,'S2_sigma = ' + str(pdfSum.getParameters(dh)[4].getValV()))
#latex.DrawText(0.5,0.60,'Background = ' + str(pdfSum.getParameters(dh)[7].getValV()))

canvas.SaveAs('B_mass_fit.png')
