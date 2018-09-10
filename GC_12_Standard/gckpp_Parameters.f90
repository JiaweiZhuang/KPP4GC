! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! 
! Parameter Module File
! 
! Generated by KPP-2.2.4_gc symbolic chemistry Kinetics PreProcessor
!       (http://www.cs.vt.edu/~asandu/Software/KPP)
! KPP is distributed under GPL, the general public licence
!       (http://www.gnu.org/copyleft/gpl.html)
! (C) 1995-1997, V. Damian & A. Sandu, CGRER, Univ. Iowa
! (C) 1997-2005, A. Sandu, Michigan Tech, Virginia Tech
!     With important contributions from:
!        M. Damian, Villanova University, USA
!        R. Sander, Max-Planck Institute for Chemistry, Mainz, Germany
! 
! File                 : gckpp_Parameters.f90
! Time                 : Mon Sep 10 15:27:56 2018
! Working directory    : /home/zhuangjw/KPP/GC_12_Standard
! Equation file        : gckpp.kpp
! Output root filename : gckpp
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



MODULE gckpp_Parameters

  USE gckpp_Precision
  PUBLIC
  SAVE


! NSPEC - Number of chemical species
  INTEGER, PARAMETER :: NSPEC = 240 
! NVAR - Number of Variable species
  INTEGER, PARAMETER :: NVAR = 235 
! NFLUX - Number of Reaction Flux species
  INTEGER, PARAMETER :: NFLUX = 1 
! NFAM - Number of Prod/Loss Families
  INTEGER, PARAMETER :: NFAM = 7 
! NVARACT - Number of Active species
  INTEGER, PARAMETER :: NVARACT = 205 
! NFIX - Number of Fixed species
  INTEGER, PARAMETER :: NFIX = 5 
! NREACT - Number of reactions
  INTEGER, PARAMETER :: NREACT = 725 
! NVARST - Starting of variables in conc. vect.
  INTEGER, PARAMETER :: NVARST = 1 
! NFIXST - Starting of fixed in conc. vect.
  INTEGER, PARAMETER :: NFIXST = 236 
! NONZERO - Number of nonzero entries in Jacobian
  INTEGER, PARAMETER :: NONZERO = 2768 
! LU_NONZERO - Number of nonzero entries in LU factoriz. of Jacobian
  INTEGER, PARAMETER :: LU_NONZERO = 3441 
! CNVAR - (NVAR+1) Number of elements in compressed row format
  INTEGER, PARAMETER :: CNVAR = 236 
! NLOOKAT - Number of species to look at
  INTEGER, PARAMETER :: NLOOKAT = 0 
! NMONITOR - Number of species to monitor
  INTEGER, PARAMETER :: NMONITOR = 0 
! NMASS - Number of atoms to check mass balance
  INTEGER, PARAMETER :: NMASS = 1 

! Index declaration for variable species in C and VAR
!   VAR(ind_spc) = C(ind_spc)

  INTEGER, PARAMETER :: ind_CH2I2 = 1 
  INTEGER, PARAMETER :: ind_CH2ICl = 2 
  INTEGER, PARAMETER :: ind_CH2IBr = 3 
  INTEGER, PARAMETER :: ind_AERI = 4 
  INTEGER, PARAMETER :: ind_CO2 = 5 
  INTEGER, PARAMETER :: ind_INDIOL = 6 
  INTEGER, PARAMETER :: ind_ISALA = 7 
  INTEGER, PARAMETER :: ind_ISALC = 8 
  INTEGER, PARAMETER :: ind_ISN1OA = 9 
  INTEGER, PARAMETER :: ind_ISN1OG = 10 
  INTEGER, PARAMETER :: ind_LBRO2H = 11 
  INTEGER, PARAMETER :: ind_LBRO2N = 12 
  INTEGER, PARAMETER :: ind_LISOPOH = 13 
  INTEGER, PARAMETER :: ind_LISOPNO3 = 14 
  INTEGER, PARAMETER :: ind_LTRO2H = 15 
  INTEGER, PARAMETER :: ind_LTRO2N = 16 
  INTEGER, PARAMETER :: ind_LVOCOA = 17 
  INTEGER, PARAMETER :: ind_LVOC = 18 
  INTEGER, PARAMETER :: ind_LXRO2H = 19 
  INTEGER, PARAMETER :: ind_LXRO2N = 20 
  INTEGER, PARAMETER :: ind_MSA = 21 
  INTEGER, PARAMETER :: ind_PYAC = 22 
  INTEGER, PARAMETER :: ind_SO4H1 = 23 
  INTEGER, PARAMETER :: ind_SO4H2 = 24 
  INTEGER, PARAMETER :: ind_SOAGX = 25 
  INTEGER, PARAMETER :: ind_SOAIE = 26 
  INTEGER, PARAMETER :: ind_SOAME = 27 
  INTEGER, PARAMETER :: ind_IMAE = 28 
  INTEGER, PARAMETER :: ind_SOAMG = 29 
  INTEGER, PARAMETER :: ind_POx = 30 
  INTEGER, PARAMETER :: ind_LOx = 31 
  INTEGER, PARAMETER :: ind_PCO = 32 
  INTEGER, PARAMETER :: ind_LCO = 33 
  INTEGER, PARAMETER :: ind_PSO4 = 34 
  INTEGER, PARAMETER :: ind_LCH4 = 35 
  INTEGER, PARAMETER :: ind_PH2O2 = 36 
  INTEGER, PARAMETER :: ind_I2O4 = 37 
  INTEGER, PARAMETER :: ind_DHDN = 38 
  INTEGER, PARAMETER :: ind_DHDC = 39 
  INTEGER, PARAMETER :: ind_I2O2 = 40 
  INTEGER, PARAMETER :: ind_MONITA = 41 
  INTEGER, PARAMETER :: ind_BENZ = 42 
  INTEGER, PARAMETER :: ind_CH3CCl3 = 43 
  INTEGER, PARAMETER :: ind_CH3I = 44 
  INTEGER, PARAMETER :: ind_H1301 = 45 
  INTEGER, PARAMETER :: ind_H2402 = 46 
  INTEGER, PARAMETER :: ind_I2O3 = 47 
  INTEGER, PARAMETER :: ind_PMNN = 48 
  INTEGER, PARAMETER :: ind_PPN = 49 
  INTEGER, PARAMETER :: ind_TOLU = 50 
  INTEGER, PARAMETER :: ind_BrNO2 = 51 
  INTEGER, PARAMETER :: ind_CCl4 = 52 
  INTEGER, PARAMETER :: ind_CFC11 = 53 
  INTEGER, PARAMETER :: ind_CFC12 = 54 
  INTEGER, PARAMETER :: ind_CFC113 = 55 
  INTEGER, PARAMETER :: ind_CFC114 = 56 
  INTEGER, PARAMETER :: ind_CFC115 = 57 
  INTEGER, PARAMETER :: ind_H1211 = 58 
  INTEGER, PARAMETER :: ind_IBr = 59 
  INTEGER, PARAMETER :: ind_IEPOXD = 60 
  INTEGER, PARAMETER :: ind_INO = 61 
  INTEGER, PARAMETER :: ind_N2O = 62 
  INTEGER, PARAMETER :: ind_TRO2 = 63 
  INTEGER, PARAMETER :: ind_BRO2 = 64 
  INTEGER, PARAMETER :: ind_IEPOXA = 65 
  INTEGER, PARAMETER :: ind_IEPOXB = 66 
  INTEGER, PARAMETER :: ind_IONITA = 67 
  INTEGER, PARAMETER :: ind_N = 68 
  INTEGER, PARAMETER :: ind_OCS = 69 
  INTEGER, PARAMETER :: ind_XRO2 = 70 
  INTEGER, PARAMETER :: ind_HI = 71 
  INTEGER, PARAMETER :: ind_MAP = 72 
  INTEGER, PARAMETER :: ind_CHBr3 = 73 
  INTEGER, PARAMETER :: ind_ICl = 74 
  INTEGER, PARAMETER :: ind_CH2Cl2 = 75 
  INTEGER, PARAMETER :: ind_IMAO3 = 76 
  INTEGER, PARAMETER :: ind_CHCl3 = 77 
  INTEGER, PARAMETER :: ind_MPN = 78 
  INTEGER, PARAMETER :: ind_Cl2O2 = 79 
  INTEGER, PARAMETER :: ind_CH2Br2 = 80 
  INTEGER, PARAMETER :: ind_ETP = 81 
  INTEGER, PARAMETER :: ind_HCFC123 = 82 
  INTEGER, PARAMETER :: ind_ClNO2 = 83 
  INTEGER, PARAMETER :: ind_HCFC141b = 84 
  INTEGER, PARAMETER :: ind_HCFC142b = 85 
  INTEGER, PARAMETER :: ind_IONO = 86 
  INTEGER, PARAMETER :: ind_HCFC22 = 87 
  INTEGER, PARAMETER :: ind_OIO = 88 
  INTEGER, PARAMETER :: ind_RA3P = 89 
  INTEGER, PARAMETER :: ind_RB3P = 90 
  INTEGER, PARAMETER :: ind_XYLE = 91 
  INTEGER, PARAMETER :: ind_DMS = 92 
  INTEGER, PARAMETER :: ind_CH3Cl = 93 
  INTEGER, PARAMETER :: ind_CH3Br = 94 
  INTEGER, PARAMETER :: ind_HNO4 = 95 
  INTEGER, PARAMETER :: ind_ClOO = 96 
  INTEGER, PARAMETER :: ind_HNO2 = 97 
  INTEGER, PARAMETER :: ind_OClO = 98 
  INTEGER, PARAMETER :: ind_PAN = 99 
  INTEGER, PARAMETER :: ind_RP = 100 
  INTEGER, PARAMETER :: ind_PP = 101 
  INTEGER, PARAMETER :: ind_PRPN = 102 
  INTEGER, PARAMETER :: ind_SO4 = 103 
  INTEGER, PARAMETER :: ind_ALK4 = 104 
  INTEGER, PARAMETER :: ind_PIP = 105 
  INTEGER, PARAMETER :: ind_R4P = 106 
  INTEGER, PARAMETER :: ind_HPALD = 107 
  INTEGER, PARAMETER :: ind_BrCl = 108 
  INTEGER, PARAMETER :: ind_C3H8 = 109 
  INTEGER, PARAMETER :: ind_DHPCARP = 110 
  INTEGER, PARAMETER :: ind_HOI = 111 
  INTEGER, PARAMETER :: ind_IAP = 112 
  INTEGER, PARAMETER :: ind_HPC52O2 = 113 
  INTEGER, PARAMETER :: ind_VRP = 114 
  INTEGER, PARAMETER :: ind_ATOOH = 115 
  INTEGER, PARAMETER :: ind_Br2 = 116 
  INTEGER, PARAMETER :: ind_HC187 = 117 
  INTEGER, PARAMETER :: ind_MOBA = 118 
  INTEGER, PARAMETER :: ind_HONIT = 119 
  INTEGER, PARAMETER :: ind_DHMOB = 120 
  INTEGER, PARAMETER :: ind_RIPB = 121 
  INTEGER, PARAMETER :: ind_BrSALC = 122 
  INTEGER, PARAMETER :: ind_ISNP = 123 
  INTEGER, PARAMETER :: ind_MP = 124 
  INTEGER, PARAMETER :: ind_BrSALA = 125 
  INTEGER, PARAMETER :: ind_MAOP = 126 
  INTEGER, PARAMETER :: ind_MRP = 127 
  INTEGER, PARAMETER :: ind_RIPA = 128 
  INTEGER, PARAMETER :: ind_RIPD = 129 
  INTEGER, PARAMETER :: ind_EOH = 130 
  INTEGER, PARAMETER :: ind_ETHLN = 131 
  INTEGER, PARAMETER :: ind_N2O5 = 132 
  INTEGER, PARAMETER :: ind_INPN = 133 
  INTEGER, PARAMETER :: ind_MTPA = 134 
  INTEGER, PARAMETER :: ind_MTPO = 135 
  INTEGER, PARAMETER :: ind_NPMN = 136 
  INTEGER, PARAMETER :: ind_C2H6 = 137 
  INTEGER, PARAMETER :: ind_IONO2 = 138 
  INTEGER, PARAMETER :: ind_MOBAOO = 139 
  INTEGER, PARAMETER :: ind_DIBOO = 140 
  INTEGER, PARAMETER :: ind_IPMN = 141 
  INTEGER, PARAMETER :: ind_LIMO = 142 
  INTEGER, PARAMETER :: ind_H = 143 
  INTEGER, PARAMETER :: ind_BrNO3 = 144 
  INTEGER, PARAMETER :: ind_MACRNO2 = 145 
  INTEGER, PARAMETER :: ind_ROH = 146 
  INTEGER, PARAMETER :: ind_I2 = 147 
  INTEGER, PARAMETER :: ind_MONITS = 148 
  INTEGER, PARAMETER :: ind_Cl2 = 149 
  INTEGER, PARAMETER :: ind_ISOPNB = 150 
  INTEGER, PARAMETER :: ind_CH4 = 151 
  INTEGER, PARAMETER :: ind_ISNOHOO = 152 
  INTEGER, PARAMETER :: ind_MVKOO = 153 
  INTEGER, PARAMETER :: ind_ISNOOB = 154 
  INTEGER, PARAMETER :: ind_GAOO = 155 
  INTEGER, PARAMETER :: ind_CH3CHOO = 156 
  INTEGER, PARAMETER :: ind_IEPOXOO = 157 
  INTEGER, PARAMETER :: ind_GLYX = 158 
  INTEGER, PARAMETER :: ind_MVKN = 159 
  INTEGER, PARAMETER :: ind_MGLYOO = 160 
  INTEGER, PARAMETER :: ind_PRN1 = 161 
  INTEGER, PARAMETER :: ind_MONITU = 162 
  INTEGER, PARAMETER :: ind_MGLOO = 163 
  INTEGER, PARAMETER :: ind_A3O2 = 164 
  INTEGER, PARAMETER :: ind_PROPNN = 165 
  INTEGER, PARAMETER :: ind_MAN2 = 166 
  INTEGER, PARAMETER :: ind_ISNOOA = 167 
  INTEGER, PARAMETER :: ind_PO2 = 168 
  INTEGER, PARAMETER :: ind_ISOPNDO2 = 169 
  INTEGER, PARAMETER :: ind_HCOOH = 170 
  INTEGER, PARAMETER :: ind_B3O2 = 171 
  INTEGER, PARAMETER :: ind_MACROO = 172 
  INTEGER, PARAMETER :: ind_R4N1 = 173 
  INTEGER, PARAMETER :: ind_ISOP = 174 
  INTEGER, PARAMETER :: ind_MAOPO2 = 175 
  INTEGER, PARAMETER :: ind_H2O2 = 176 
  INTEGER, PARAMETER :: ind_ATO2 = 177 
  INTEGER, PARAMETER :: ind_I = 178 
  INTEGER, PARAMETER :: ind_RCO3 = 179 
  INTEGER, PARAMETER :: ind_LIMO2 = 180 
  INTEGER, PARAMETER :: ind_MACRN = 181 
  INTEGER, PARAMETER :: ind_OLND = 182 
  INTEGER, PARAMETER :: ind_OLNN = 183 
  INTEGER, PARAMETER :: ind_IO = 184 
  INTEGER, PARAMETER :: ind_KO2 = 185 
  INTEGER, PARAMETER :: ind_HOBr = 186 
  INTEGER, PARAMETER :: ind_ISOPNBO2 = 187 
  INTEGER, PARAMETER :: ind_PIO2 = 188 
  INTEGER, PARAMETER :: ind_HC5OO = 189 
  INTEGER, PARAMETER :: ind_HNO3 = 190 
  INTEGER, PARAMETER :: ind_ISOPND = 191 
  INTEGER, PARAMETER :: ind_GLYC = 192 
  INTEGER, PARAMETER :: ind_NMAO3 = 193 
  INTEGER, PARAMETER :: ind_ACTA = 194 
  INTEGER, PARAMETER :: ind_VRO2 = 195 
  INTEGER, PARAMETER :: ind_HOCl = 196 
  INTEGER, PARAMETER :: ind_CH2OO = 197 
  INTEGER, PARAMETER :: ind_ISN1 = 198 
  INTEGER, PARAMETER :: ind_ClNO3 = 199 
  INTEGER, PARAMETER :: ind_MGLY = 200 
  INTEGER, PARAMETER :: ind_ACET = 201 
  INTEGER, PARAMETER :: ind_HC5 = 202 
  INTEGER, PARAMETER :: ind_RIO2 = 203 
  INTEGER, PARAMETER :: ind_ETO2 = 204 
  INTEGER, PARAMETER :: ind_INO2 = 205 
  INTEGER, PARAMETER :: ind_R4O2 = 206 
  INTEGER, PARAMETER :: ind_R4N2 = 207 
  INTEGER, PARAMETER :: ind_HAC = 208 
  INTEGER, PARAMETER :: ind_MRO2 = 209 
  INTEGER, PARAMETER :: ind_BrO = 210 
  INTEGER, PARAMETER :: ind_PRPE = 211 
  INTEGER, PARAMETER :: ind_RCHO = 212 
  INTEGER, PARAMETER :: ind_MEK = 213 
  INTEGER, PARAMETER :: ind_CH2O = 214 
  INTEGER, PARAMETER :: ind_MACR = 215 
  INTEGER, PARAMETER :: ind_ALD2 = 216 
  INTEGER, PARAMETER :: ind_MVK = 217 
  INTEGER, PARAMETER :: ind_MCO3 = 218 
  INTEGER, PARAMETER :: ind_SO2 = 219 
  INTEGER, PARAMETER :: ind_CO = 220 
  INTEGER, PARAMETER :: ind_MO2 = 221 
  INTEGER, PARAMETER :: ind_Br = 222 
  INTEGER, PARAMETER :: ind_NO = 223 
  INTEGER, PARAMETER :: ind_HBr = 224 
  INTEGER, PARAMETER :: ind_HCl = 225 
  INTEGER, PARAMETER :: ind_O1D = 226 
  INTEGER, PARAMETER :: ind_Cl = 227 
  INTEGER, PARAMETER :: ind_O = 228 
  INTEGER, PARAMETER :: ind_NO3 = 229 
  INTEGER, PARAMETER :: ind_NO2 = 230 
  INTEGER, PARAMETER :: ind_O3 = 231 
  INTEGER, PARAMETER :: ind_HO2 = 232 
  INTEGER, PARAMETER :: ind_ClO = 233 
  INTEGER, PARAMETER :: ind_OH = 234 
  INTEGER, PARAMETER :: ind_H2O = 235 

! Index declaration for fixed species in C
!   C(ind_spc)

  INTEGER, PARAMETER :: ind_H2 = 236 
  INTEGER, PARAMETER :: ind_MOH = 237 
  INTEGER, PARAMETER :: ind_N2 = 238 
  INTEGER, PARAMETER :: ind_O2 = 239 
  INTEGER, PARAMETER :: ind_RCOOH = 240 

! Index declaration for fixed species in FIX
!    FIX(indf_spc) = C(ind_spc) = C(NVAR+indf_spc)

  INTEGER, PARAMETER :: indf_H2 = 1 
  INTEGER, PARAMETER :: indf_MOH = 2 
  INTEGER, PARAMETER :: indf_N2 = 3 
  INTEGER, PARAMETER :: indf_O2 = 4 
  INTEGER, PARAMETER :: indf_RCOOH = 5 

END MODULE gckpp_Parameters
