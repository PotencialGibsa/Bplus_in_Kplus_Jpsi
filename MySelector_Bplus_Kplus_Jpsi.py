from ROOT import *
#from variables import *
import glob
import numpy as np
import math

def sqrt(var):
      return math.sqrt(abs(var))
#############################################################
PDG_MUON_MASS    =   0.1056583745
PDG_PION_MASS    =   0.13957061
PDG_PIOZ_MASS    =   0.1349770
PDG_KAON_MASS    =   0.493677
PDG_PROTON_MASS  =   0.9382720813
PDG_KSHORT_MASS  =   0.497611
PDG_KSHORT_DM    =   0.000013
PDG_KSHORT_TIME  =   0.8954 * 0.0000000001
PDG_KS_MASS      =   PDG_KSHORT_MASS
PDG_LAMBDA_MASS  =   1.115683
PDG_LAMBDA_DM    =   0.000006
PDG_LAMBDA_TIME  =   2.632 * 0.0000000001
PDG_SIGMA0_MASS  =   1.192642
PDG_XImunus_MASS =   1.32171
PDG_XImunus_DM   =   0.00007
PDG_XImunus_TIME =   1.639 * 0.0000000001
PDG_OMmunus_MASS =   1.67245
PDG_OMmunus_DM   =   0.00029
PDG_OMmunus_TIME =   0.821 * 0.0000000001
PDG_DPM_MASS     =   1.86965
PDG_DPM_DM       =   0.00005
PDG_DPM_TIME     =   1.040 * 0.000000000001
PDG_DZ_MASS      =   1.86483
PDG_DZ_DM        =   0.00005
PDG_DZ_TIME      =   0.4101 * 0.000000000001
PDG_DS_MASS      =   1.96834
PDG_DS_DM        =   0.00007
PDG_DS_TIME      =   0.504 * 0.000000000001
PDG_LAMCP_MASS   =   2.28646
PDG_LAMCP_DM     =   0.00031
PDG_LAMCP_TIME   =   2.00 * 0.0000000000001
PDG_XICZ_MASS    =   2.47087
PDG_XICZ_DM      =   0.00031
PDG_XICZ_TIME    =   1.12 * 0.0000000000001
PDG_XICP_MASS    =   2.46787
PDG_XICP_DM      =   0.00030
PDG_XICP_TIME    =   4.42 * 0.0000000000001
PDG_KSTARZ_MASS  =   0.89555
PDG_KSTARZ_GAMMA =   0.0473
PDG_KSTARP_MASS  =   0.89176
PDG_KSTARP_GAMMA =   0.0503
PDG_PHI_MASS     =   1.019461
PDG_PHI_GAMMA    =   0.004249
PDG_JPSI_MASS    =   3.096900
PDG_PSI2S_MASS   =   3.686097
PDG_X3872_MASS   =   3.87169
PDG_BU_MASS      =   5.27932
PDG_BU_TIME      =   1.638 * 0.000000000001
PDG_B0_MASS      =   5.27963
PDG_B0_TIME      =   1.520 * 0.000000000001
PDG_BS_MASS      =   5.36689
PDG_BS_TIME      =   1.509 * 0.000000000001
PDG_BC_MASS      =   6.2749
PDG_BC_TIME      =   0.507 * 0.000000000001
PDG_LB_MASS      =   5.61960
PDG_LB_TIME      =   1.470 * 0.000000000001
PDG_XIBZ_MASS    =   5.7919
PDG_XIBZ_TIME    =   1.479 * 0.000000000001
PDG_XIBM_MASS    =   5.7970
PDG_XIBM_TIME    =   1.571 * 0.000000000001
PDG_OMBM_MASS    =   6.0461
PDG_OMBM_TIME    =   1.64 * 0.000000000001
PDG_C            =   29979245800. ### in cm/c
PDG_DSTR         =   2.01026
###  }}}
'''
md      = RooRealVar ( "md"     ,"M(D) [GeV]"               , PDG_DZ_MASS-0.03 , PDG_DZ_MASS+0.03  )
mds     = RooRealVar ( "mds"    ,"M(Dstar) [GeV]"           , 2.004, 2.019 )
#
dspt    = RooRealVar ( "dspt"   ,"dspt"                     , 3.0   , 33.0 )
dzpt    = RooRealVar ( "dzpt"   ,"dzpt"                     , 2.5   , 22.5  )
kspt    = RooRealVar ( "kspt"   ,"kspt"                     , 0.0   , 10.0 )
pipt    = RooRealVar ( "pipt"   ,"pipt"                     , 0.1   , 1.5 )
#
dset    = RooRealVar ( "dset"   ,"dset"                     , -2.6  , 2.6   )
dzet    = RooRealVar ( "dzet"   ,"dzet"                     , -2.6  , 2.6   )
kset    = RooRealVar ( "kset"   ,"kset"                     , -2.6  , 2.6   )
piet    = RooRealVar ( "piet"   ,"piet"                     , -2.6  , 2.6   )
#
dzds2   = RooRealVar ( "dzds2"  ,"DetSign"                  , 0.0   , 50    )
dzds3   = RooRealVar ( "dzds3"  ,"DetSign"                  , 0.0   , 80    )
#
dsvtxp  = RooRealVar ( "dsvtxp" ,"dsvtxp"                   , -0.1  , 1.1   )
dzvtxp  = RooRealVar ( "dzvtxp" ,"dzvtxp"                   , -0.1  , 1.1   )
#
ipsmin  = RooRealVar ( "ipsmin" ,"ipsmin"                   , 0.0   , 1000.0)
vari    = RooRealVar ( "vari"   ,"vari"                     , 0.0   , 5)
vdr     = RooRealVar ( "vdr"    ,"vdr"                      , 0.0   , 2)
vdz     = RooRealVar ( "vdz"    ,"vdz"                      , 0.0   , 5)
'''
def DetachSignificance2(vtx, vtxE1, vtxE2):
    return sqrt( vtx.X()**2 / (vtxE1.X()**2 + vtxE2.X()**2) + vtx.Y()**2 / (vtxE1.Y()**2 + vtxE2.Y()**2))

def DetachSignificance3(vtx, vtxE1, vtxE2):
    return sqrt( vtx.X()**2 / (vtxE1.X()**2 + vtxE2.X()**2) + vtx.Y()**2 / (vtxE1.Y()**2 + vtxE2.Y()**2) + vtx.Z()**2 / (vtxE1.Z()**2 + vtxE2.Z()**2))

def DirectionCos2 (v1, v2):
    r1 = sqrt(v1.X()**2 + v1.Y()**2)
    r2 = sqrt(v2.X()**2 + v2.Y()**2)
    return ( v1.X() * v2.X() + v1.Y() * v2.Y() ) / (r1*r2 + 0.0000001)

def DirectionCos3 (v1, v2):
    r1 = sqrt(v1.X()**2 + v1.Y()**2 + v1.Z()**2)
    r2 = sqrt(v2.X()**2 + v2.Y()**2 + v2.Z()**2)
    return ( v1.X() * v2.X() + v1.Y() * v2.Y() + v1.Z() * v2.Z() ) / (r1*r2 + 0.0000001)

def DirectionChi22 (vtx0, vtx0E, vtx1, vtx1E, P, PE):
    dvtx    = vtx1 - vtx0
    dvtxE   = TVector3( sqrt(vtx0E.X()**2 + vtx1E.X()**2), sqrt(vtx0E.Y()**2 + vtx1E.Y()**2), sqrt(vtx0E.Z()**2 + vtx1E.Z()**2) )
    Pscaled = P * (dvtx.Mag() / P.Mag())
    PscaledE= PE * (dvtx.Mag() / P.Mag())
    return DetachSignificance2 (Pscaled - dvtx, PscaledE, dvtxE)

def DirectionChi23 (vtx0, vtx0E, vtx1, vtx1E, P, PE):
    dvtx    = vtx1 - vtx0 ## vertex difference
    dvtxE   = TVector3( sqrt(vtx0E.X()**2 + vtx1E.X()**2), sqrt(vtx0E.Y()**2 + vtx1E.Y()**2), sqrt(vtx0E.Z()**2 + vtx1E.Z()**2) ) ## its error
    Pscaled = P * (dvtx.Mag() / P.Mag()) ## scaled momentum to be the same length as vertex difference
    PscaledE= PE * (dvtx.Mag() / P.Mag()) ## its error
    return DetachSignificance3 (Pscaled - dvtx, PscaledE, dvtxE)


def TH1NormBinWidth(hh):
    _htemp = hh.Clone()
    _int = _htemp.Integral()
    for _i in range(1, hh.GetNbinsX()+1):
        _htemp[_i] = _htemp[_i] / _htemp.GetBinWidth(_i)
    #
    _nam = hh.GetName() + '_no'
    _htemp.SetName('_nam')
    _htemp.Scale(_int / _htemp.Integral())
    return _htemp

def stri(n):
    if n>=0 and n<10:
        return '0'+str(n)
    else:
        return str(n)

MyFileNames = glob.glob('/eos/home-d/dshmygol/Bfinder/2018_UL_sanua/2018/Charmonium/*/*/000?/BFinder*.root')
print('The number of files is %i for %s' % (len(MyFileNames), iter))
print('\nAdding the files to the chain')
chain = TChain('wztree')
for name in MyFileNames:
    chain.Add(name)
_fileOUT = '/eos/home-d/dshmygol/SimpleFile_Bplus_Kplus_Jpsi_2018_A.root'
fileOUT = TFile (_fileOUT, "recreate")
mytree = TTree("mytree","mytree")

nevt = chain.GetEntries()
print('Number of events: %s' % (nevt))
my_vars = [
    'SAMEEVENT',
    #Mass B+
    'B_mass',  'B_masc0',
    'B_pt', 'B_eta',
    'B_pvdistsignif2', 'B_pvdistsignif3',
    'B_pvdist', 'B_vtxprob',
    'B_pvcos2', 'B_pvcos3',
    #J/psi
    'J_psi_mass', 'J_psi_prob',
    'J_psi_pt',
    #Kaon
    'K1_pt', 'K1_ips', 'K1_chrg',
    #mu1
    'mu1_pt',
    #mu2
    'mu2_pt',
]

for _var in my_vars:
    exec(_var+'=np.zeros(1, dtype=float)')
    exec('mytree.Branch("' + _var + '"' + ' '*(25-len(_var)) + ',' + _var + ' '*(25-len(_var)) + ', "'+ _var + '/D")')

### finished with deploying vars for SimpleFile }}}

# P4 ROOT Lorentz vector
K1P4,K1P4_CV, MU1P4, MU2P4, JPSIP4, BP4 = [TLorentzVector() for i in range(6)]

BBB = 0

# looping all over the events
for evt in range(nevt):
    if chain.GetEntry(evt) <= 0:
        break
    # securing that _nCand is set to 0 for every event before assign the real chain value
    _nCand = 0
    _nCand = chain.nCand
    if _nCand < 1:
        continue
    for cand in range(_nCand):
        #it only works when we set variable with [0] by its side
        # assigning PV info
        PV = TVector3(chain.PV_becos_XX[cand], chain.PV_becos_YY[cand], chain.PV_becos_ZZ[cand])
        PVE = TVector3(sqrt(chain.PV_becos_EX[cand]), sqrt(chain.PV_becos_EY[cand]), sqrt(chain.PV_becos_EZ[cand]))

        # assigning Kaon info
        K1P4.SetXYZM (chain.K1_px[cand], chain.K1_py[cand], chain.K1_pz[cand], PDG_PION_MASS)
        K1P4_CV.SetXYZM (chain.K1_px_CV[cand], chain.K1_py_CV[cand], chain.K1_pz_CV[cand], PDG_KAON_MASS)
        K1_pt[0] = K1P4.Pt()
        K1_ips[0] = chain.K1_ips[cand]

        # assigning muon+ info
        MU1P4.SetXYZM(chain.B_mu_px1[cand], chain.B_mu_px1[cand], chain.B_mu_px1[cand], PDG_MUON_MASS)
        mu1_pt[0] = MU1P4.Pt()

        # assigning muon- info
        MU2P4.SetXYZM(chain.B_mu_px2[cand], chain.B_mu_px2[cand], chain.B_mu_px2[cand], PDG_MUON_MASS)
        mu2_pt[0] = MU2P4.Pt()

        #    # assigning J/psi info
        J_psi_mass[0] = chain.B_J_mass[cand]
        J_psi_prob[0] = chain.B_J_Prob[cand]
        JPSIP4.SetXYZM(chain.B_J_px[cand], chain.B_J_py[cand], chain.B_J_pz[cand], chain.B_J_mass[cand])
        J_psi_pt[0] = JPSIP4.Pt()
        # assigning info about B+
        BV = TVector3(chain.B_DecayVtxX[cand], chain.B_DecayVtxY[cand], chain.B_DecayVtxZ[cand])
        BVE = TVector3(sqrt(chain.B_DecayVtxXE[cand]), sqrt(chain.B_DecayVtxYE[cand]), sqrt(chain.B_DecayVtxZE[cand]))
        BP4.SetXYZM( chain.B_px[cand], chain.B_py[cand], chain.B_pz[cand], chain.B_mass[cand])
        BP3 = BP4.Vect()
        B_pt[0]= BP4.Pt()
        B_eta[0]= BP4.Eta()
        B_vtxprob[0] = chain.B_Prob[cand]
        B_mass[0]= chain.B_mass[cand]
        B_masc0[0]= chain.B_mass_c0[cand]
        B_pvdistsignif2[0] = DetachSignificance2( BV - PV, BVE, PVE)
        B_pvdistsignif3[0] = DetachSignificance3( BV - PV, BVE, PVE)
        B_pvdist[0] = (BV - PV).Mag()
        B_pvcos2[0] = DirectionCos2 ( BV - PV, BP3)
        B_pvcos3[0] = DirectionCos3 ( BV - PV, BP3)
        if B_vtxprob[0] < 0.05 : continue

        SAMEEVENT[0] = 0;
        if (BBB > -1):
            SAMEEVENT[0] = 1
        BBB = 1; ## the next candidate is not the 1st in event

        # fill _var branches in the tree if passes preselection cuts
        mytree.Fill()
                #
    BBB = -1 ## when loop over candidates in event is finished, set this to -1, so the next candidate has SAMEEVENT=0
    if (evt % 10001 == 0):
        print('Running... Now in event %i / %i ' %(evt, nevt))


print('Total entries stored in MyTree: %s' %(mytree.GetEntries()))

fileOUT.Write()
