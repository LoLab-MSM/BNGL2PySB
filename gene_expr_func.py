from pysb import *
Model()
Parameter('a',     4       )
Parameter('b',     10      )
Parameter('gam',   10      )
Parameter('kap0',  0.6     )
Parameter('kap1',  0.2     )
Parameter('d1',   0.0005      )
Parameter('d0',   gam*d1      )
Parameter('k0',   kap0*d1     )
Parameter('k1',   kap1*d1     )
Parameter('v0',   a*d1        )
Parameter('v1',   b*d0        )
Monomer('DNA', ['promoter'], {'promoter': '1)'})
Monomer('mRNA', [')'])
Monomer('Protein', [')'])
Monomer('Src', [')'])
Monomer('Null', [')'])
Observable("obs_DNA_Active", DNA(promoter=1))
Observable("obs_DNA_Total", DNA())
Observable("obs_mRNA_Total", mRNA())
Observable("obs_Protein_Total", Protein())
Initial(DNA(promoter=0), Parameter("DNA(promoter=0_0",1))
Initial(mRNA(), Parameter("mRNA_0",0))
Initial(Protein(), Parameter("Protein_0",0))
Initial(Src(), Parameter("Src_0",1))
Initial(Null(), Parameter("Null_0",0))
Expression("cn_mRNA", v0*DNA_Activ )
Expression("cn_Prot", v1*mRNA_Tota )
Rule("rule_0", DNA(promoter=0)   >>   DNA(promoter=1)    , k0,k1)
Rule("rule_1", Src()   >>   Src() + mRNA()       , fcn_mRNA)
Rule("rule_2", Src()   >>   Src() + Protein()    , fcn_Prot)
Rule("rule_3", mRNA()      >>   Null()   , d0)
Rule("rule_4", Protein()   >>   Null()   , d1)