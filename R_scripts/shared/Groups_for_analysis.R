purinergic_P2rx = c("P2rx1", "P2rx2", "P2rx3", "P2rx4", "P2rx5", "P2rx6", "P2rx7")

purinergic_P2ry = c("P2ry1", "P2ry2", "P2ry4", "P2ry6", "P2ry10", "P2ry12", "P2ry13", "P2ry14")

ccr_family = c("Ccr1", "Ccr2", "Ccr3", "Ccr4", "Ccr5", "Ccr6", "Ccr7", "Ccr8", "Ccr9", "Ccr10", "Cx3cr1")

cxc_receptors = c("Cxcr1", "Cxcr2", "Cxcr3", "Cxcr4", "Cxcr5", "Cxcr6")

fc_receptors = c("Fcgr1", "Fcgr2b", "Fcgr3", "Fcgr4", "Fcgrt", "Fcer1a", "Fcer1g", "Fcrl1", "Fcrl6")

ifitms = c("Ifitm1", "Ifitm2", "Ifitm3", "Ifitm5", "Ifitm6", "Ifitm7")

TLRs = c("Tlr1", "Tlr2", "Tlr3", "Tlr4", "Tlr5", "Tlr6", "Tlr7", "Tlr8", "Tlr9", "Tlr11", "Tlr12", "Tlr13")

siglecs = c("Siglece", "Siglecf", "Siglecg", "Siglech", "Siglec1", "Cd33", "Siglec15")

MMPs = c("Mmp1a", "Mmp1b", "Mmp2", "Mmp3", "Mmp7", "Mmp8", "Mmp9", "Mmp10", "Mmp11", "Mmp12", "Mmp13", "Mmp14", "Mmp15", 
         "Mmp16", "Mmp17", "Mmp19", "Mmp20", "Mmp21", "Mmp23", "Mmp24", "Mmp25", "Mmp27", "Mmp28")

IHC_genes = c("Aif1", "Cd74", "Arg1", "Itgam")

mhc_co_receptors = c("Cd80", "Cd86", "Cd274", "Pdcd1lg2")

# GSEA subgroups

GSEA_mouse_IFN_gamma = c( "Adar", "Apol6", "Arid5b", "Arl4a", "Auts2", "B2m", "Bank1", "Batf2", "Bpgm", "Btg1", "C1rb", "C1ra", "C1s1", "C1s2", "Casp1", "Casp3", "Casp4",
                          "Casp7", "Casp8", "Ccl2", "Ccl5", "Ccl7", "Cd274", "Cd38", "Cd40", "Cd69", "Cd74", "Cd86", "Cdkn1a", "Cfb", "Cfh", "Ciita", "Cmklr1", "Cmpk2", "Csf2rb2", 
                          "Csf2rb", "Cxcl10", "Ddx58", "Ddx60", "Dhx58", "Eif2ak2", "Eif4e3", "Epsti1", "Fas", "Fcgr1", "Fgl2", "Fpr1", "Cmtr1", "Gbp3", "Gbp9", "Gbp10", "Gbp11", 
                          "Gbp4", "Gbp6", "Gch1", "Gpr18", "Gzma", "Herc6", "Hif1a", "H2-Q4", "H2-D1", "H2-Q6", "H2-Q7", "H2-Q10", "H2-Q2", "H2-Q1", "H2-Bl", "H2-K1", "Gm8909", 
                          "Gm10499", "H2-DMa", "H2-Aa", "H2-M3", "Icam1", "Ido1", "Ifi30", "Ifi35", "Ifi44", "Ifi44l", "Ifih1", "Ifit1bl1", "Ifit1", "Ifit2", "Ifit3", "Ifit3", 
                          "Ifitm2", "Ifitm3", "Ifnar2", "Il10ra", "Il15", "Il15ra", "Il18bp", "Il2rb", "Il4ra", "Il6", "Il7", "Irf1", "Irf2", "Irf4", "Irf5", "Irf7", "Irf8", "Irf9", 
                          "Isg15", "Isg20", "Isoc1", "Itgb7", "Jak2", "Klrk1", "Lap3", "Lats2", "Lcp2", "Lgals3bp", "Ly6e", "Lysmd2", "March1", "Mettl7b", "Mt2", "Mthfd2", "Mvp", 
                          "Myd88", "Nampt", "Ncoa3", "Nfkb1", "Nfkbia", "Nlrc5", "Nmi", "Nod1", "Nup93", "Oas2", "Oas3", "Oasl1", "Ogfr", "P2ry14", "Parp12", "Parp14", "Pde4b", 
                          "Peli1", "Pfkp", "Pim1", "Pla2g4a", "Plscr1", "Pml", "Pnp", "Pnp2", "Pnpt1", "Helz2", "Psma2", "Psma3", "Psmb10", "Psmb2", "Psmb8", "Psmb9", "Psme1", 
                          "Psme2", "Ptgs2", "Ptpn1", "Ptpn2", "Ptpn6", "Rapgef6", "Rbck1", "Ripk1", "Ripk2", "Rnf213", "Rnf31", "Rsad2", "Rtp4", "Samd9l", "Samhd1", "Sectm1a", 
                          "Selp", "Serping1", "Slamf7", "Slc25a28", "Socs1", "Socs3", "Sod2", "Sp110", "Sppl2a", "Sri", "Sspn", "St3gal5", "St8sia4", "Stat1", "Stat2", 
                          "Stat3", "Stat4", "Tap1", "Tapbp", "Tdrd7", "Tnfaip2", "Tnfaip3", "Tnfaip6", "Tnfsf10", "Tor1b", "Trafd1", "Trim14", "Trim21", "Trim25", "Trim26", "Txnip", 
                          "Ube2l6", "Upp1", "Usp18", "Vamp5", "Vamp8", "Vcam1", "Wars", "Xaf1", "Xcl1", "Zbp1", "Znfx1")

GSEA_mouse_IL6_STAT3 = c( "A2m", "Acvr1b", "Acvrl1", "Bak1", "Cbl", "Ccl7", "Ccr1", "Cd14", "Cd36", "Cd38", "Cd44", "Cd9", "Cntfr", "Crlf2", "Csf1", "Csf2", "Csf2ra", "Csf2rb", "Csf2rb2",
                          "Csf3r", "Cxcl10", "Cxcl13", "Cxcl2", "Dntt", "Ebi3", "Fas", "Grb2", "Hax1", "Hmox1", "Ifnar1", "Ifngr1", "Ifngr2", "Il10rb", "Il12rb1", "Il13ra1", "Il15ra", 
                          "Il17ra", "Il17rb", "Il18r1", "Il1b", "Il1r1", "Il1r2", "Il2ra", "Il2rg", "Il3ra", "Il4ra", "Il6", "Il6st", "Il7", "Il9r", "Inhbe", "Irf1", "Irf9", "Itga4", 
                          "Itgb3", "Jun", "Lepr", "Ltb", "Ltbr", "Map3k8", "Myd88", "Osmr", "Pdgfc", "Pf4", "Pik3r5", "Pim1", "Ptpn1", "Ptpn11", "Ptpn2", "Reg1", "Socs1", "Socs3", 
                          "Stam2", "Stat1", "Stat2", "Stat3", "Tgfb1", "Tlr2", "Tnf", "Tnfrsf12a", "Tnfrsf1a", "Tnfrsf1b", "Tnfrsf21", "Tyk2", "Cxcl1", "Cxcl11", "Cxcl9", "Pla2g2a")

GSEA_mouse_TGFbeta1 = c( "Acvr1", "Apc", "Arid4b", "Bcar3", "Bmp2", "Bmpr1a", "Bmpr2", "Cdh1", "Cdk9", "Cdkn1c", "Ctnnb1", "Eng", "Fkbp1a", "Fnta", "Furin", "Hdac1", "Hipk2", "Id1", 
                         "Id2", "Id3", "Ifngr2", "Junb", "Klf10", "Lefty2", "Ltbp2", "Ncor2", "Nog", "Pmepa1", "Ppm1a", "Ppp1ca", "Ppp1r15a", "Rab31", "Rhoa", "Serpine1", "Ski", "Skil", 
                         "Slc20a1", "Smad1", "Smad3", "Smad6", "Smad7", "Smurf1", "Smurf2", "Sptbn1", "Tgfb1", "Tgfbr1", "Tgif1", "Thbs1", "Tjp1", "Trim33", "Ube2d3", "Wwtr1", "Xiap",
                         "Map3k7")
  
GSEA_mouse_IL4_IL13 = c("Clec7a", "Chil3", "Spp1", "Lgals3", "Cd302", "F13a1", "Fgl2", "Cxcr4", "Mrc1", "Lta4h", "Ccl8", "Cxcl2", "Il2ra", "Il10", "Tgfbi", "Tlr1", "Msr1", "Cd163", 
                        "Car2", "Tgfb1", "Egr2", "Mmp12", "Marco", "Retnla", "Tlr8", "Arg1", "Gas7", "Il1rn", "Alox15", "Adk", "Igf1", "Cxcl1", "Cxcl3", "Ccl26", "Ccl17", "Ccl20")  

## XUE Sets converted to mouse genes

# Xue IL4 set to mouse
Xue_IL4_mouse = c("A2m", "Acot7", "Acta2", "Adam15", "Adam19", "Alox15", "Ankrd37", "Arhgap10", "Asap1", "Ascl2", "Auh", "Avpi1", "Batf3", "Bcar3", "1810010H24Rik", 
                  "BC028528", "C1qb", "Card9", "Cbr3", "Ccl2", "Ccl17", "Ccl22", "Ccl5", "Ccnd1", "Ccnd2", "Cd209e", "Cdkn1a", "Cdr2l", "Chn2", "Mgl2", "Clec10a", 
                  "Clec4g", "Cmtm8", "Col9a2", "Crip1", "Cst7", "Ctnnal1", "Dennd4c", "Dhrs11", "Dtna", "Egln3", "Espnl", "Evl", "F13a1", "Fabp5", "Fcer1a", "Fcer2a", 
                  "Fcrlb", "Fgl2", "Foxq1", "Fscn1", "Fzd2", "Gadd45a", "Galm", "Galnt12", "Gfra2", "Gpd1l", "Hlx", "Hmg20b", "Homer2", "Hsph1", "Igf2bp2", "Il1r2", 
                  "Il1rn", "Ints3", "Itpripl2", "Krt17", "Lamp3", "Lipa", "Lpl", "Lrrfip1", "Maoa", "Mat2a", "Matk", "Mmp12", "Myc", "Nagpa", "Ndp", "Nfil3", "Palld", 
                  "Pdxk", "Pfkp", "Pon2", "Ppm1f", "Ppp1r14a", "Ptrf", "Qprt", "Rab33a", "Rab7b", "Ramp1", "Rap1gap", "Rasgrp3", "Rasl10a", "Rassf7", "Rftn1", "Rhbdf1", 
                  "Rnase1", "Rrp1b", "Rtkn", "Sema4a", "Shpk", "Siglecg", "Sla", "Slc27a3", "Slc47a1", "Slc7a8", "Socs1", "Sox8", "Spint2", "Sptan1", "St6gal1", 
                  "Sucnr1", "Syt17", "Tacstd2", "Tgm2", "Tmem71", "Tnfrsf4", "Tsku", "Tspan32", "Vash1", "Vcl", "Wnt5a", "Zfp366",
                  "Ifn2", "Scimp", "Dennd1b", "Mcur1", "Ccl26", "Cd1d1", "Cd1d2", "Cd52", "Galnt18", "2900026A02Rik", "Nmnat3", "Plpp3", "Tmem194", "Tpm2")

Xue_IL4_mouse = unique(Xue_IL4_mouse)

# Xue IFNg set to mouse
Xue_IFNg_mouse = c("Acsl1", "Adamdec1", "Aim2", "Aldh1a2", "Ankrd22", "Aqp9", "Arid5b", "Ascl2", "Asns", "Atf3", "AA467197", "C2", "Car12", "Car2", "Camp", "Ccl22", 
                   "Ccl5", "Ccnd1", "Ccnd2", "Cd274", "Cd38", "Cd69", "Cd83", "Cd86", "Cfb", "Crabp2", "Crispld2", "Cxcl10", "Cyp27b1", "Dhrs9", "Dusp5", "Edn1", "Eno2", 
                   "Epsti1", "Exog", "F3", "Fam20c", "Fam26f", "Fbp1", "Fcgr1", "Fem1c", "Gbp2", "Gbp3", "Gbp5", "Gch1", "Gk", "Gpc4", "Grem1", "Hapln3", "Gm7030", 
                   "Gm11127", "Hspa1a", "Ifih1", "Ifit2", "Ifit3", "Ifit3b", "Il15", "Il1rn", "Il6", "Irf1", "Irf7", "Irf8", "Isg20", "Itgb7", "Jak2", "Lap3", "Lgals3bp", 
                   "Lilr4b", "Lilrb4a", "Lrrc32", "Matk", "Mthfd2", "Mthfsl", "Mthfs", "Ncf1", "Nfs1", "Ninj1", "Nrip3", "P2rx7", "Parp14", "Parp9", "Pde4b", "Ppargc1b", 
                   "Psat1", "Psmb9", "Psme2", "Ptpn1", "Rhbdf1", "Rhbdf2", "Ripk2", "Rnf19b", "Rsad2", "Samd9l", "Scg5", "Serpine1", "Serping1", "Slc1a5", "Slc6a12", 
                   "Slc7a5", "Slco4a1", "Snx10", "Sod2", "Sphk1", "Stat1", "Stx11", "Tap1", "Tap2", "Tgm2", "Tm4sf19", "Tmem140", "Tmem38b", "Tnfaip6", "Tnfsf10", 
                   "Traf3ip2", "Tsc22d1", "Tsku", "Ube2l6", "Vamp5", "Wars", "Xaf1", "Zc3h12a", "Apol7a", "Bhlhe41", "Rgcc", "Ocstamp", "Map3k7cl", "Ccl8", "Cxcl9", 
                   "Myof", "Gbp2b", "Hfe", "Ifi27", "Ifitm1", "Vwa8", "Mt1", "Mt2", "Mx1", "Pcnx")

Xue_IFNg_mouse = unique(Xue_IFNg_mouse)


# Xue IL10 set to mouse
# Changed:  
Xue_IL10_mouse = c("Acot7", "Acsl1", "Acta2", "Adamtsl4", "Adcy3", "Agrp", "Akr1c18", "Aldh1a2", "Ang", "Ankdd1a", "Ankrd33", "Apoe", "Aqp9", "Arrdc4", "Atp6v0e2", 
                   "Avpi1", "Birc3", "Bok", "Bscl2", "Camp", "Ccdc106", "Ccl2", "Ccl24", "Ccl3", "Cd83", "Cfd", "Ch25h", "Chst13", "Cldn23", "Clec12a", "Coro2a", "Crabp2", 
                   "Ctnnal1", "Ctsg", "Cxcl16", "Cxcr4", "Cxxc5", "Dhrs9", "Dpep2", "Ehd1", "Evl", "Fabp5", "Fbp1", "Mfsd7c", "Folr2", "Fuca1", "Gadd45a", "Gadd45g", "Gchfr", 
                   "Gfra2", "Gstm2", "Gstm7", "H2-Ab1", "Hmg20b", "Hspa1a", "Il17rb", "Il1r2", "Il1rn", "Il3ra", "Irf1", "Kat2a", "Klf11", "Klhdc8b", "Lipa", "Lpl", "Maoa", 
                   "Map2k3", "Matk", "Mmp7", "Ms4a4a", "Myc", "Nenf", "Ngfrap1", "Ninj1", "Nr1h3", "Nrp2", "Nupr1", "Pcolce2", "Pcsk5", "Pld3", "Plin2", "Pnpla7", "Ralgds", 
                   "Rbp1", "Rilpl1", "Rnase1", "Rnase4", "Sema6b", "Sepp1", "Sil1", "Slc11a1", "Slc12a8", "Slc16a10", "Slc44a2", "Slco2b1", "Sod2", "Sparc", "Spocd1", "Sulf2", 
                   "Tfrc", "Tgm2", "Tmem158", "Tnfrsf14", "Traf3ip2", "Tsku", "Tspan4", "Vat1", "Vegfb", "Vmo1", "Zfand2a", "Zmynd15",
                   "Acp5", "Apoc1", "Rgcc", "Ifn2", "Ccl26", "Ccl3", "Ccl4", "Cd52", "Ctsl", "Cxcl1", "Cxcl5", "Tmem236", "Mrc1", "Naprt", "S100a8", "Ppbp", "Dcstamp")

Xue_IL10_mouse = unique(Xue_IL10_mouse)

# MHC_class_II_genes
MHC_class_II_genes = c("Ciita", "Cd74", "H2-Aa", "H2-DMb1", "H2-Eb1")

# Sensome updated
microglia_sensome_updated = c("P2ry12", "P2ry13", "P2ry6", "Gpr34", "Adora3", "Entpd1", "Tmem173", "A630033H20Rik", "Csf1r", "Csf3r", "Tgfbr1", "Tgfbr2", "Ifngr1", "Il10ra", "Il6ra", "Il21r", "Tnfrsf17", "Tnfrsf1b",
                              "Cx3cr1", "Ccr5", "C3ar1", "Ptafr", "C5ar2", "Cmklr1", "Cysltr1", "Ccrl2", "Cmtm6", "C5ar1", "Fcgr3", "Fcer1g", "Fcgr2b", "Fcgr1", "Cmtm7", "Fcrl1", "Fcgr4", "Selplg", "Ly86",
                              "Cd68", "Trem2", "Cd180", "Tlr2", "Cd37", "Tlr7", "Cd14", "Clec4a3", "Tlr4", "Tlr13", "Clec5a", "Havcr2", "Clec7a", "Cxcl16", "Cd48", "Ltf", "Cd74", "Upk1b", "Tlr12", "Tlr1", "Pilra", "Tlr6", "Ifitm6",
                              "Itgam", "Itgb2", "Itgb5", "Adgre1", "Ecscr", "Lair1", "Siglech", "Slco2b1", "Slc2a5", "Lgals9", "Gpr183", "Tmem37", "Cd33", "Gpr84", "Slc7a7", "Cd52", "Siglece", "Cd79b", "Slc16a3", 
                              "Icam1", "Icam4", "Cd84", "Lag3", "Cd86", "Ptprc", "Tyrobp", "Tnfrsf13b", "Tnfrsf17", "Cd22", "Tmem119", "Cd53", "Vsir", "Slamf9", "I830077J02Rik", "Clec4a2", "Clec4b1", "Lilra5", "Tmem8c", "Gpr160", "Cd101")

## Sensome subsets

# Sensome metabolic
sensome_metabolic = c("Gpr183", "Adora3", "Il6ra", "Cx3cr1", "P2ry12", "P2ry13", "Csf1r", "Csf3r" )

# Sensome misc
sensome_misc = c("Tmem8c", "C5ar1", "Tmem37", "Gpr84", "Ly86", "Cd180", "Slc2a5", "Tmem119", "Slco2b1", "Cmtm7", "Tlr13")

# Microglia Debris
microglia_debris = c("Cd93", "Msr1", "Cd36", "Olr1", "Megf10", "Clec7a", "Scarf1")