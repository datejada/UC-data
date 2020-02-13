$Title Stochastic hydrothermal security constraint UC and simulation model

$OnText

Developed by

   Diego Tejada (based on model by Andrés Ramos)
   Instituto de Investigación Tecnológica
   Escuela Técnica Superior de Ingeniería - ICAI
   UNIVERSIDAD PONTIFICIA COMILLAS
   Alberto Aguilera 23
   28015 Madrid, Spain

$OffText

*-------------------------------------------------------------------------------
*                              MODEL OPTIONS
*-------------------------------------------------------------------------------

* To use a algebraic formulation (formulation without defining the sets)
$OnEmpty OnMulti OffListing

* definition of symbol for comments at the end of the line
$EOLCOM //

* optimizer definition
option   lp = cplex ;
option  mip = cplex ;
option rmip = cplex ;

* general options
option optcr    =   1e-3 ;   // tolerance to solve MIP until IntGap < OptcR
option reslim   =   1800 ;   // maximum run time [sec]
option threads  =     -1 ;   // number of cores
option solprint =    off ;   // print the final solution in the .lst file
option limrow   =      0 ;   // maximum number of equations in the .lst file
option limcol   =      0 ;   // maximum number of variables in the .lst file
option bratio   =      0 ;   // bratio = 0 forces GAMS to always try to construct a basis
option savepoint=      1 ;   // save into a gdx file solution (0=no save, 1=only the last one, 2=for each solve)

* profile options
option profile=1, profileTol = 0.01 ;

*-------------------------------------------------------------------------------
*                     MODEL SETS, PARAMETERS, VARIABLES
*-------------------------------------------------------------------------------

* definitions
sets
   n      "hour                 "
   n1(n)  "first hour of the day"
   sc     "scenario             "

   g      "generating unit "
   t (g)  "thermal    unit "
   h (g)  "hydro      plant"

   l      "type of startup"

*Network sets
   i             "bus                                              "
   c             "circuit ID                                       "
   iic (i,i,c)   "lines in service in steady state                 "
   ii2 (i,i  )   "pair of busses with at least one line connection "
   jjc (i,i,c)   "lines with enable contingency   flag activated   "
   kkc (i,i,c)   "lines in service in steady state alias           "
   is  (i    )   "busses without slack bus                         "
   gi  (g  ,i)   "generator connected to a bus                     "

   cont(sc,n,i,i,c,i,i,c) "contingency set for each scenario and hour"

   alias (n,nn), (sc,scc), (i,ii,j,jj), (c,cc)

parameters
   p2ndResUpCost              "cost factor of the up   secondary reserve" /  0.4 /
   p2ndResDwCost              "cost factor of the down secondary reserve" /  0.2 /
   pDuration                  "duration of load level         [h]       " /    1 /
   pDemand         (n)        "hourly load                    [GW]      "
   pOperReserve    (n)        "hourly operating reserve       [GW]      "
   pDemand2ndResUp (n)        "hourly operating reserve up    [GW]      "
   pDemand2ndResDw (n)        "hourly operating reserve down  [GW]      "
   pIntermGen      (n,i,sc)   "stochastic IG generation       [GW]      "
   pScenProb       (  sc  )   "probability of scenarios       [p.u.]    "
   pCommitt        (   g,n)   "commitment  of the unit        [0-1]     "

   pMaxProd           (g)     "maximum output                 [GW]        "
   pMinProd           (g)     "minimum output                 [GW]        "
   pMaxCons           (g)     "maximum consumption            [GW]        "
   pIniOut            (g)     "initial output > min load      [GW]        "
   pIniUC             (g)     "initial commitment             [0-1]       "
   pRampUp            (g)     "ramp up                        [GW per   h]"
   pRampDw            (g)     "ramp down                      [GW per   h]"
   pMinTU             (g)     "minimun time down              [h]         "
   pMinTD             (g)     "minimun time up                [h]         "
   pSlopeVarCost      (g)     "slope     variable cost        [M€ per GWh]"
   pInterVarCost      (g)     "intercept variable cost        [M€ per   h]"
   pShutdownCost      (g)     "shutdown           cost        [M€]        "
   pENSCost                   "energy not served  cost        [M€ per GWh]"

   pDurationSD        (g    ) "shutdown duration time         [h] "
   pDurationSU        (g,  l) "startup  duration time         [h] "
   pTMinSU            (g,  l) "min duration time for SU type  [h] "
   pTypeSUCost        (g,  l) "startup cost for SU type       [M€]"
   pProductRampDw     (g,n  ) "ramp down output               [GW]"
   pProductRampUp     (g,n,l) "ramp up   output               [GW]"
   pNumSUTypes        (g    ) "number of active SU types          "

   pTotalEnergy       (g)     "total energy available per hydro unit [MWh]"

   pMax2ndResUP       (g)     "maximum 2nd reserve up         [GW]"
   pMin2ndResUP       (g)     "minimum 2nd reserve up         [GW]"
   pTot2ndResUP       (g)     "total   2nd reserve up         [GW]"
   pMax2ndResDW       (g)     "maximum 2nd reserve down       [GW]"
   pMin2ndResDW       (g)     "minimum 2nd reserve down       [GW]"
   pTot2ndResDW       (g)     "total   2nd reserve down       [GW]"

*Network Parameters
   pSlack                           "slack bus. Bus 1 as default                                           " / 1 /
   pSBase                           "S Base                                                          [GW]  "
   pNetworkConst                    "binary parameter to consider or not the network constraints     [0-1] "
   pEnableIterAlg                   "binary parameter to enable iterative algorithm                  [0-1] "
   pBusProfile    (i             )  "load profile per bus                                            [p.u.]"
   pYBUS          (i,ii          )  "Susceptance matrix                                              [p.u.]"
   pYBUSInv       (i,ii          )  "Susceptance matrix inverse Z = B^-1                             [p.u.]"
   pPmax          (i,ii,c        )  "maximum Power Transfer of line                                  [GW]  "
   pPctmax        (i,ii,c        )  "maximum Power Transfer of line during overload emergency state  [GW]  "
   pISF           (i,ii,c,i      )  "injection Shift Factors (for node i to line i->ii circuitID c)  [p.u.]"
   pPTDF          (i,ii,c,j,jj   )  "Power Transfer Distribution Factor                              [p.u.]"
   pLODF          (i,ii,c,j,jj,cc)  "Line  Outage   Distribution Factor                              [p.u.]"
   pXline         (i,ii,c        )  "Reactance  X of line                                            [p.u.]"
   pRline         (i,ii,c        )  "Resistance R of line                                            [p.u.]"
   pZline         (i,ii,c        )  "Impedance  Z of line                                            [p.u.]"

*N-1 Iterative Algorithm Parameters
   pCirPFContAux1 (sc,n,i,ii,c,  j,jj,cc) "Auxiliary calculation of contingency power flow [GW] "
   pCirPFContAux2 (sc,n,i,ii,c,  j,jj,cc) "Auxiliary calculation of contingency power flow [0-1]"
   pCPUTime                               "CPU Time for each iteration                          "
   pFCL                                   "Additional factor for Contingency Limit              "   /1.05/
   aux                                    "Auxiliar parameter to control iterations             "
   iter                                   "Number of iteration needed                           "

*State-Transition Parameters
   pModelForm                             "Formulation: 0-TC, 1-ST, 2-2bin                 [0-2]"

*Output Parameters
   pSummary       (*                    ) "Results summary                                                      "
   pProduct       (sc,g,n               ) "output      of the unit                                   [MW]       "
   pConsump       (sc,g,n               ) "consumption of the unit                                   [MW]       "
   pIG            (sc,n,i               ) "intermittent generation                                   [MW]       "
   pENS           (sc,n,i               ) "energy not served per node                                [MW]       "
   pSRMC          (sc,n,i               ) "short run marginal cost                                   [€ per MWh]"
   p2ndResUP      (   g,n               ) "up   secondary reserve                                    [MW]       "
   p2ndResDW      (   g,n               ) "down secondary reserve                                    [MW]       "
   p2ResUpMC      (     n               ) "marginal cost of up   secondary reserve                   [€ per MW] "
   p2ResDwMC      (     n               ) "marginal cost of down secondary reserve                   [€ per MW] "
   pPrintCirPF    (sc,  i,ii,c          ) "Lines to be printed in steady state                                  "
   pCirPF         (sc,n,i,ii,c          ) "Power flow in steady state                                [MW]       "
   pPrintCirPFCont(sc,  i,ii,c,  j,jj,cc) "Lines to be printed in contingency state                             "
   pCirPFCont     (sc,  i,ii,c,n,j,jj,cc) "Power flow in line iic under N-1 of line jjc              [MW]       "

variables
   vTotalVCost                     "total system variable cost      [M€]"
*Network variables
   vCirPF    (sc,n,i,ii,c        ) "power flow through a line       [GW]"
   vCirPFCont(sc,n,i,ii,c,j,jj,cc) "power flow through a line i->ii CircuitID c due to N-1 in line j->jj Circuit ID cc [GW]"

binary   variables
   vCommitt  (   n,g)  "commitment of the unit         [0-1]"
   vStartUp  (   n,g)  "startup    of the unit         [0-1]"
   vShutDown (   n,g)  "shutdown   of the unit         [0-1]"
   vTypeSU   ( n,g,l)  "type of start up ramp          [0-1]"

positive variables
   vProduct  (sc,n,g) "output of the unit              [GW]"
   vProduct1 (sc,n,g) "output of the unit > minload    [GW]"
   vConsump  (sc,n,g) "consumption of the unit         [GW]"
   v2ndResUP (   n,g) "2nd res. up   allocation        [GW]"
   v2ndResDW (   n,g) "2nd res. down allocation        [GW]"
   vENS      (sc,n,i) "energy not served per node      [GW]"
   vIG       (sc,n,i) "intermittent generation per bus [GW]"

*-------------------------------------------------------------------------------
*                                EQUATIONS
*-------------------------------------------------------------------------------

equations
   eTotalVCost          "total system variable cost     [M€]"
*Reserve requirements Constraints
   e2ReserveUp(   n   ) "operating reserve upwards      [GW]"
   e2ReserveDw(   n   ) "operating reserve downwards    [GW]"
*Tight and Compact UC Formulation
   eUCStrShut (   n, g) "relation among committment startup and shutdown"
   eMaxOutput (sc,n, g) "max output of a committed unit [GW]            "
   eMinOutput (sc,n, g) "min output of a committed unit [GW]            "
   eTotOutput (sc,n, g) "tot output of a committed unit [GW]            "
   eMinTUp    (sc,n, g) "minimun up time   (    committed)              "
   eMinTDown  (sc,n, g) "minimun down time (not committed)              "
   eRampUp    (sc,n, g) "bound on ramp up               [GW]            "
   eRampDw    (sc,n, g) "bound on ramp down             [GW]            "
   eTypeSU1   (n ,g, l) "Activation of start up type                    "
   eTypeSU2   (n ,g   ) "Just one start up type is selected             "
*State Transition Formulation
   eSTF_UCStrShut (   n,g) "state-transition constraints for UC                             "
   eSTF_MaxOutput (sc,n,g) "max output of a committed unit in ST formulation [GW]           "
   eSTF_MinOutput (sc,n,g) "min output of a committed unit in ST formulation [GW]           "
   eSTF_TotOutput (sc,n,g) "tot output of a committed unit in ST formulation [GW]           "
   eSTF_MinTUp    (sc,n,g) "minimun up time                in ST formulation (    committed)"
   eSTF_MinTDown  (sc,n,g) "minimun down time              in ST formulation (not committed)"
   eSTF_RampUp    (sc,n,g) "bound on ramp up               in ST formulation [GW]           "
   eSTF_RampDw    (sc,n,g) "bound on ramp down             in ST formulation [GW]           "
   eSTF_TypeSU1   (n ,g,l) "Activation of start up type                                     "
   eSTF_TypeSU2   (n ,g  ) "Just one start up type is selected                              "
*Projected two-binary-variable Formulation
   e2bin_UCStrShut(   n, g) "relation among committment startup and shutdown"
   e2bin_MaxOutput(sc,n, g) "max output of a committed unit [GW]            "
   e2bin_MinOutput(sc,n, g) "min output of a committed unit [GW]            "
   e2bin_TotOutput(sc,n, g) "tot output of a committed unit [GW]            "
   e2bin_MinTUp   (sc,n, g) "minimun up time   (    committed)              "
   e2bin_MinTDown (sc,n, g) "minimun down time (not committed)              "
   e2bin_RampUp   (sc,n, g) "bound on ramp up               [GW]            "
   e2bin_RampDw   (sc,n, g) "bound on ramp down             [GW]            "
   e2bin_TypeSU1  (n ,g, l) "Activation of start up type                    "
   e2bin_TypeSU2  (n ,g   ) "Just one start up type is selected             "
*Hydro units constraints
   eTotEnergy (sc,g   ) "Energy available per hydro unit[GWh]"
   eTot2ResUP (   g   ) "upward   secondary reserve          "
   eTot2ResDW (   g   ) "downward secondary reserve          "
*Network constraints
   eBalanceBus(sc,n, i)             "load generation balance per bus                    [GW]"
   eBalance   (sc,n   )             "load generation balance                            [GW]"
   eCirPF     (sc,n,i,ii,c        ) "power flow through a line                          [GW]"
   eCirPFCont (sc,n,i,ii,c,j,jj,cc) "power flow through line iic due to N-1 in line jjc [GW]"
;

*--------------------------- Objective Function --------------------------------

eTotalVCost.. vTotalVCost =e=
                sum[(sc,n,i   ) $[pScenProb(sc)*pDemand(n)], pENSCost                      *pDuration*vENS     (sc,n,i  )] +
                sum[(sc,n,g   ) $[pScenProb(sc)*pDemand(n)], pSlopeVarCost(g)              *pDuration*vProduct (sc,n,g  )] +
                sum[(   n,g   ) $               pDemand(n) , pInterVarCost(g)              *pDuration*vCommitt (   n,g  )]   $[pModelForm <> 1] +
                sum[(   n,g   ) $               pDemand(n) , pInterVarCost(g)              *pDuration*[vStartUp(   n,g  )  +
                                                                                                       vCommitt(   n,g  )]]  $[pModelForm =  1] +
                sum[(   n,g,l ) $               pDemand(n) , pTypeSUCost  (g,l)                      *vTypeSU  (   n,g,l)] +
                sum[(   n,g   ) $               pDemand(n) , pShutdownCost(g)                        *vShutDown(   n,g  )]   $[pModelForm <  2] +
                sum[(   n,g   ) $               pDemand(n) , pShutdownCost(g)                        *[vCommitt( n-1,g  )  -
                                                                                                       vCommitt(   n,g  )  +
                                                                                                       vStartUp(   n,g  )  +
                                                                                                       pIniUC  (g)$n1(n)]]   $[pModelForm =  2] +
                sum[(   n,g   ) $               pDemand(n) , pInterVarCost(g)*p2ndResUpCost*pDuration*v2ndResUP(   n,g  )]   $[pModelForm <  2] +
                sum[(   n,g   ) $               pDemand(n) , pInterVarCost(g)*p2ndResDwCost*pDuration*v2ndResDW(   n,g  )]   $[pModelForm <  2] +
                sum[(   n,g   ) $               pDemand(n) , pInterVarCost(g)*p2ndResUpCost*pDuration*v2ndResUP(   n,g  )  *
                                                                                                (pMaxProd(g)-pMinProd(g))]   $[pModelForm =  2] +
                sum[(   n,g   ) $               pDemand(n) , pInterVarCost(g)*p2ndResDwCost*pDuration*v2ndResDW(   n,g  )  *
                                                                                                (pMaxProd(g)-pMinProd(g))]   $[pModelForm =  2] ;

*------------------- Tight and Compact UC Formulation --------------------------

eUCStrShut (   n,t)               $[pDemand(n) and pModelForm=0]  ..
   vCommitt(n,t) - vCommitt(n-1,t) - pIniUC(t) $n1(n) =e= vStartUp(n,t) - vShutDown(n,t) ;

eMaxOutput (sc,n,t) $[pScenProb(sc)*pDemand(n) and pModelForm=0] ..
   vProduct1(sc,n,t) + v2ndResUP(n,t) =l= vCommitt(n,t)*(pMaxProd(t)-pMinProd(t)) ;

eMinOutput (sc,n,t) $[pScenProb(sc)*pDemand(n) and pModelForm=0] ..
   vProduct1(sc,n,t) - v2ndResDW(n,t) =g= 0 ;

eTotOutput (sc,n,t) $[pScenProb(sc)*pDemand(n) and pModelForm=0] ..
   vProduct(sc,n,t) =e= pMinProd(t)*vCommitt(n,t) + vProduct1(sc,n,t) +
   sum[(nn  ) $(ord(nn)<=pDurationSD(t  )), [vShutDown(n-[ord(nn)-1                 ],t  )*pProductRampDw(t,nn  )] $pDemand(n-[ord(nn)-1                 ])] +
   sum[(nn,l) $(ord(nn)<=pDurationSU(t,l)), [vTypeSU  (n-[ord(nn)-pDurationSU(t,l)-1],t,l)*pProductRampUp(t,nn,l)] $pDemand(n-[ord(nn)-pDurationSU(t,l)-1])] ;

eMinTUp    (sc,n,t) $[pScenProb(sc)*pDemand(n) and pModelForm=0] ..
   sum[nn $(ord(nn)>=ord(n)+1-pMinTU(t) and ord(nn)<=ord(n) and pDemand(nn)), vStartUp (nn,t)] =l=     vCommitt(n,t) ;

eMinTDown  (sc,n,t) $[pScenProb(sc)*pDemand(n) and pModelForm=0] ..
   sum[nn $(ord(nn)>=ord(n)+1-pMinTD(t) and ord(nn)<=ord(n) and pDemand(nn)), vShutDown(nn,t)] =l= 1 - vCommitt(n,t) ;

eRampUp    (sc,n,t) $[pScenProb(sc)*pDemand(n) and pModelForm=0] ..
   vProduct1(sc,n,t) - vProduct1(sc,n-1,t) - (pIniOut(t)-pMinProd(t)*pIniUC(t)) $(n1(n) and pIniUC(t)) + v2ndResUP(n,t) =l=   pRampUp(t)*pDuration*[vCommitt(n,t)-vStartUp(n,t)$n1(n)] ;

eRampDw    (sc,n,t) $[pScenProb(sc)*pDemand(n) and pModelForm=0] ..
   vProduct1(sc,n,t) - vProduct1(sc,n-1,t) - (pIniOut(t)-pMinProd(t)*pIniUC(t)) $(n1(n) and pIniUC(t)) - v2ndResDW(n,t) =g= - pRampDw(t)*pDuration*[vCommitt(n,t)+vShutDown(n,t)];

eTypeSU1   ( n,t,l) $[ord(l)<pNumSUTypes(t) and pDemand(n) and pModelForm=0] ..
   vTypeSU (n,t,l) =l= sum[nn $(ord(nn)>=pTMinSU(t,l) and ord(nn)<=pTMinSU(t,l+1)-1 and pDemand(nn)), vShutDown(n-ord(nn),t)$pDemand(n-ord(nn))] ;

eTypeSU2   ( n,t  )                    $[pDemand(n) and pModelForm=0]  ..
   vStartUp(n,t  ) =e= sum[l, vTypeSU(n,t,l)] ;

*-------------------- State-Transition UC Formulation --------------------------

eSTF_UCStrShut(   n, t) $[              pDemand(n) and pModelForm=1] ..
   vStartUp(n-1,t) + vCommitt(n-1,t) + pIniUC(t) $n1(n) =e=  vShutDown(n,t) + vCommitt(n,t) ;

eSTF_MaxOutput (sc,n,t) $[pScenProb(sc)*pDemand(n) and pModelForm=1] ..
   vProduct1(sc,n,t) + v2ndResUP(n,t) =l= [vStartUp(n,t) + vCommitt(n,t)] *(pMaxProd(t)-pMinProd(t)) ;

eSTF_MinOutput (sc,n,t) $[pScenProb(sc)*pDemand(n) and pModelForm=1] ..
   vProduct1(sc,n,t) - v2ndResDW(n,t) =g= 0 ;

eSTF_TotOutput (sc,n,t) $[pScenProb(sc)*pDemand(n) and pModelForm=1] ..
   vProduct(sc,n,t) =e= pMinProd(t)*[vStartUp(n,t) + vCommitt(n,t)] + vProduct1(sc,n,t) +
   sum[(nn  ) $(ord(nn)<=pDurationSD(t  )), [vShutDown(n-[ord(nn)-1                 ],t  )*pProductRampDw(t,nn  )] $pDemand(n-[ord(nn)-1                 ])] +
   sum[(nn,l) $(ord(nn)<=pDurationSU(t,l)), [vTypeSU  (n-[ord(nn)-pDurationSU(t,l)-1],t,l)*pProductRampUp(t,nn,l)] $pDemand(n-[ord(nn)-pDurationSU(t,l)-1])] ;

eSTF_MinTUp    (sc,n,t) $[pScenProb(sc)*pDemand(n) and pModelForm=1] ..
   sum[nn $(ord(nn)>=ord(n)+1-pMinTU(t) and [ord(nn)<=ord(n)-1] and pDemand(nn)), vStartUp (nn,t)] =l=     vCommitt(n,t) ;

eSTF_MinTDown  (sc,n,t) $[pScenProb(sc)*pDemand(n) and pModelForm=1] ..
   sum[nn $(ord(nn)>=ord(n)+1-pMinTD(t) and ord(nn)<=ord(n) and pDemand(nn)), vShutDown(nn,t)] =l= 1 - [vStartUp(n,t) + vCommitt(n,t)] ;

eSTF_RampUp    (sc,n,t) $[pScenProb(sc)*pDemand(n) and pModelForm=1] ..
   vProduct1(sc,n,t) - vProduct1(sc,n-1,t) - (pIniOut(t)-pMinProd(t)*pIniUC(t)) $(n1(n) and pIniUC(t)) + v2ndResUP(n,t) =l=   pRampUp(t)*pDuration * [vStartUp(n,t) + vCommitt(n,t)-vStartUp(n,t)$n1(n)];

eSTF_RampDw    (sc,n,t) $[pScenProb(sc)*pDemand(n) and pModelForm=1] ..
   vProduct1(sc,n,t) - vProduct1(sc,n-1,t) - (pIniOut(t)-pMinProd(t)*pIniUC(t)) $(n1(n) and pIniUC(t)) - v2ndResDW(n,t) =g= - pRampDw(t)*pDuration *[vStartUp(n,t) + vCommitt(n,t)+vShutDown(n,t)];

eSTF_TypeSU1   ( n,t,l) $[ord(l)<pNumSUTypes(t) and pDemand(n) and pModelForm=1] ..
   vTypeSU (n,t,l) =l= sum[nn $(ord(nn)>=pTMinSU(t,l) and ord(nn)<=pTMinSU(t,l+1)-1 and pDemand(nn)), vShutDown(n-ord(nn),t)$pDemand(n-ord(nn))] ;

eSTF_TypeSU2   ( n,t  )                    $[pDemand(n) and pModelForm=1]  ..
   vStartUp(n,t  ) =e= sum[l, vTypeSU(n,t,l)] ;

*------------------- Projected two-binary-variable Formulation --------------------------

e2bin_UCStrShut (   n,t)               $[pDemand(n) and pModelForm=2]  ..
   vCommitt(n,t) - vCommitt(n-1,t) - pIniUC(t) $n1(n) =l= vStartUp(n,t) ;

e2bin_MaxOutput (sc,n,t) $[pScenProb(sc)*pDemand(n) and pModelForm=2] ..
   vProduct1(sc,n,t) + v2ndResUP(n,t) =l= vCommitt(n,t) ;

e2bin_MinOutput (sc,n,t) $[pScenProb(sc)*pDemand(n) and pModelForm=2] ..
   vProduct1(sc,n,t) - v2ndResDW(n,t) =g= 0 ;

e2bin_TotOutput (sc,n,t) $[pScenProb(sc)*pDemand(n) and pModelForm=2] ..
   vProduct(sc,n,t) =e= pMinProd(t)*vCommitt(n,t) + vProduct1(sc,n,t)*(pMaxProd(t)-pMinProd(t)) +
   sum[(nn  ) $(ord(nn)<=pDurationSD(t  )), [[vCommitt(n-[ord(nn)                   ],t)  -
                                              vCommitt(n-[ord(nn)-1                 ],t)  +
                                              vStartUp(n-[ord(nn)-1                 ],t)  +
                                         pIniUC(t) $n1(n-[ord(nn)-1                 ])  ] * pProductRampDw(t,nn  )] $pDemand(n-[ord(nn)-1                 ])] +
   sum[(nn,l) $(ord(nn)<=pDurationSU(t,l)), [vTypeSU  (n-[ord(nn)-pDurationSU(t,l)-1],t,l)* pProductRampUp(t,nn,l)] $pDemand(n-[ord(nn)-pDurationSU(t,l)-1])] ;

e2bin_MinTUp    (sc,n,t) $[pScenProb(sc)*pDemand(n) and pModelForm=2] ..
   sum[nn $(ord(nn)>=ord(n)+1-pMinTU(t) and ord(nn)<=ord(n) and pDemand(nn)), vStartUp (nn,t)] =l=     vCommitt(n,t) ;

e2bin_MinTDown  (sc,n,t) $[pScenProb(sc)*pDemand(n) and pModelForm=2] ..
   sum[nn $(ord(nn)>=ord(n)+1-pMinTD(t) and ord(nn)<=ord(n) and pDemand(nn)), vStartUp (nn,t)] =l= 1 - vCommitt(n-pMinTD(t),t) $[ord(n)-pMinTD(t)> 0]
                                                                                                     - pIniUC  (            t) $[ord(n)-pMinTD(t)<=0];

e2bin_RampUp    (sc,n,t) $[pScenProb(sc)*pDemand(n) and pModelForm=2] ..
   vProduct1(sc,n,t) - vProduct1(sc,n-1,t) - [(pIniOut(t)-pMinProd(t)*pIniUC(t))/(pMaxProd(t)-pMinProd(t))]$(n1(n) and pIniUC(t)) + v2ndResUP(n,t) =l=   [pRampUp(t)/(pMaxProd(t)-pMinProd(t))]*pDuration*[vCommitt(n  ,t)-vStartUp(n,t)            $n1(n)] ;

e2bin_RampDw    (sc,n,t) $[pScenProb(sc)*pDemand(n) and pModelForm=2] ..
   vProduct1(sc,n,t) - vProduct1(sc,n-1,t) - [(pIniOut(t)-pMinProd(t)*pIniUC(t))/(pMaxProd(t)-pMinProd(t))]$(n1(n) and pIniUC(t)) - v2ndResDW(n,t) =g= - [pRampDw(t)/(pMaxProd(t)-pMinProd(t))]*pDuration*[vCommitt(n-1,t)+vStartUp(n,t)+ pIniUC(t) $n1(n)] ;

e2bin_TypeSU1   ( n,t,l) $[ord(l)<pNumSUTypes(t) and pDemand(n) and pModelForm=2] ..
   vTypeSU (n,t,l) =l= sum[nn $(ord(nn)>=pTMinSU(t,l) and ord(nn)<=pTMinSU(t,l+1)-1 and pDemand(nn)), [vCommitt(n-[ord(nn)+1],t) -
                                                                                                       vCommitt(n- ord(nn)   ,t) +
                                                                                                       vStartUp(n- ord(nn)   ,t) +
                                                                                                  pIniUC(t) $n1(n- ord(nn)     )]$pDemand(n-ord(nn))] ;

e2bin_TypeSU2   ( n,t  ) $[pDemand(n) and pModelForm=2]  ..
   vStartUp(n,t  ) =e= sum[l, vTypeSU(n,t,l)] ;

*------------------------ Hydro Units Constraints ------------------------------

eTotEnergy (sc,h  ) $ pScenProb(sc) ..
   sum[n $pDemand(n), vProduct (sc,n,h)] =l= pTotalEnergy(h)/pDuration ;

eTot2ResUP (   h  )                 ..
   sum[n $pDemand(n), v2ndResUP(   n,h)] =l= pTot2ndResUP(h)/pDuration ;

eTot2ResDW (   h  )                 ..
   sum[n $pDemand(n), v2ndResDW(   n,h)] =l= pTot2ndResDW(h)/pDuration ;

*------------------- Reserve requirements Constraints --------------------------

e2ReserveUp(   n  )               $[pDemand(n)] ..
   + sum[g, v2ndResUP(n,g)                            ] $[pModelForm < 2]
   + sum[g, v2ndResUP(n,g) * (pMaxProd(g)-pMinProd(g))] $[pModelForm = 2]
   =g=
   pDemand2ndResUp(n) ;

e2ReserveDw(   n  )               $[pDemand(n)] ..
   + sum[g, v2ndResDW(n,g)                            ] $[pModelForm < 2]
   + sum[g, v2ndResDW(n,g) * (pMaxProd(g)-pMinProd(g))] $[pModelForm = 2]
   =g=
   pDemand2ndResDw(n) ;

*------------------------- Network Constraints ---------------------------------

eBalanceBus(sc,n,i) $[pScenProb(sc)*pDemand(n) and pNetworkConst=1] ..
   sum[gi(t,i), vProduct(sc,n,t)] + sum[gi(h,i), vProduct(sc,n,h)] - sum[gi(h,i), vConsump(sc,n,h)] +        vIG(sc,n,i) + vENS(sc,n,i) + sum[iic(ii,i,c),vCirPF(sc,n,ii,i,c)]-sum[iic(i,ii,c),vCirPF(sc,n,i,ii,c)] =g=        pDemand(n)*pBusProfile(i)  ;

eBalance   (sc,n  ) $[pScenProb(sc)*pDemand(n) and pNetworkConst=0] ..
   sum[   t   , vProduct(sc,n,t)] + sum[   h   , vProduct(sc,n,h)] - sum[   h   , vConsump(sc,n,h)] + sum[i, vIG(sc,n,i) + vENS(sc,n,i)]                                                                            =g= sum[i, pDemand(n)*pBusProfile(i)] ;

eCirPF     (sc,n,iic(i,ii,c)) $[pScenProb(sc)*pDemand(n) and pNetworkConst=1] ..
   vCirPF(sc,n,iic) =e= sum[is, pISF(iic,is)              *
                                       [  + sum[gi(t,is), vProduct(sc,n,t)]
                                          + sum[gi(h,is), vProduct(sc,n,h)]
                                          - sum[gi(h,is), vConsump(sc,n,h)]
                                          + vIG(sc,n,is)
                                          - pDemand(n)*pBusProfile(is)
                                       ]
                                   ];

eCirPFCont (sc,n,iic(i,ii,c),jjc(j,jj,cc))
                     $[cont(sc,n,iic,jjc) and pScenProb(sc)*pDemand(n) and
                       not(sameas(i,j) and sameas(ii,jj) and sameas(c,cc))] ..
                       vCirPFCont(sc,n,iic,jjc) =e= vCirPF(sc,n,iic) + pLODF(iic,jjc)*vCirPF(sc,n,jjc) ;

*-------------------------------------------------------------------------------
*                                MODELS
*-------------------------------------------------------------------------------

model SDUC / all / ;
SDUC.holdfixed = 1 ; SDUC.optfile = 1 ;

*-------------------------------------------------------------------------------
*                             DATA FROM EXCEL FILE
*-------------------------------------------------------------------------------

* read input data from Excel and include into the model
file TMP / tmp.txt /
$onecho  > tmp.txt
   i="%gams.user1%.xlsm"
   r1=indices
   o1=indices
   r2=param
   o2=param
   r3=demand
   o3=demand
   r4=oprres
   o4=oprres
   r5=oprresup
   o5=oprresup
   r6=oprresdw
   o6=oprresdw
   r7=IGgen
   o7=IGgen
   r8=thermalgen
   o8=thermalgen
   r9=hydrogen
   o9=hydrogen
   r10=demandbus
   o10=demandbus
   r11=network
   o11=network
$offecho
$call xls2gms m @"tmp.txt"

sets
$include  indices
;
$include  param
;
parameter pDemand(n)         hourly load              [MW] /
$include  demand
                                                           /
parameter pOperReserve(n)    hourly operating reserve [MW] /
$include  oprres
                                                           /
parameter pDemand2ndResUp(n) hourly operating reserve [MW] /
$include  oprresup
                                                           /
parameter pDemand2ndResDw(n) hourly operating reserve [MW] /
$include  oprresdw
                                                           /

parameter pBusProfile(i)     load profile per bus   [p.u.] /
$include  demandbus
                                                           /

table     pIntermGen(n,i,sc) stochastic IG generation [MW]
$include  IGgen
table     pThermalGen(g,*)   thermal generation data
$include  thermalgen
table     pHydroGen  (g,*)   hydro   generation data
$include  hydrogen

table     pNetwork(i,ii,c,*) Network parameters
$INCLUDE network

execute 'del tmp.txt indices param demand oprres oprresup oprresdw demandbus' ;
execute 'del tmp.txt IGgen thermalgen hydrogen network' ;

*-------------------------------------------------------------------------------
*                   OPTIONS DEFINITION FOR SOLVERS
*-------------------------------------------------------------------------------

FILE     GOPT / gurobi.opt /
PUT      GOPT / 'IIS 1'    /  'writeprob gurobi.lp'
IF(pEnableIterAlg=1,
   PUT   GOPT / 'mipstart 1' / 'lazyconstraints 1'/ 'eCirPFCont.lazy 3' /
);
PUTCLOSE GOPT ;

FILE     COPT / cplex.opt  /
PUT      COPT / 'IIS yes'  / 'WritePre presolved.lp'
IF(pEnableIterAlg=1,
   PUT   COPT / 'mipstart 1' /
);
PUTCLOSE COPT ;

*-------------------------------------------------------------------------------
*                      Injection Sensitivity Factors
*-------------------------------------------------------------------------------

* Line conecttions
iic (i,ii,c) $ pNetwork(i,ii,c,'InService') = YES                ;
ii2 (i,ii  )                                = SUM[c,iic(i,ii,c)] ;

*Busses without the slack bus
is(i) $[ORD(i) <> pSlack] = YES;

* Obtaining the parameters
pXline (i,ii,c) = pNetwork(i,ii,c,'X');

* Obtaining susceptance matrix with multiple circuits defined as i->j and j->i
pYBUS(ii2(i,ii)) = -SUM[c $ (iic(i,ii,c) AND pXline (i,ii,c)),1/pXline(i,ii,c)]
                   -SUM[c $ (iic(ii,i,c) AND pXline (ii,i,c)),1/pXline(ii,i,c)];

* Creating the symmetric matrix
pYBUS(ii,i) $ii2(i,ii) = pYBUS(i,ii);

* Creating the diagonal for the B matriz
pYBUS(ii,ii) = SUM(i,-pYBUS(i,ii));

*Eliminaing the Slack/reference node
pYBUS(i,ii) $[NOT(is(i) )] = 0;
pYBUS(i,ii) $[NOT(is(ii))] = 0;

* Creating gdx files for inverse calculation
EXECUTE_UNLOAD 'gdxfrominverse.gdx' is ;
EXECUTE_UNLOAD 'gdxforinverse.gdx ' is ;

*Obtaining the inverse of pYBUS and saving it in pYBUSInv
EXECUTE_UNLOAD 'gdxforinverse.gdx' is,pYBUS;
EXECUTE 'invert gdxforinverse.gdx is pYBUS gdxfrominverse.gdx pYBUSInv';
EXECUTE_LOAD 'gdxfrominverse.gdx' , pYBUSInv;

* Obtaining the Injection Shift Factors
pISF(i,ii,c,is) $ [pXline(i,ii,c) AND iic(i,ii,c)]
                = (pYBUSInv(i,is) - pYBUSInv(ii,is))/pXline(i,ii,c);
* eliminating the pISF lower than 1E-9 (contribution lower than 1W)
pISF(i,ii,c,is) $ (ABS(pISF(i,ii,c,is)) < 1E-9) = 0;

*-------------------------------------------------------------------------------
*                        SET VALUE VARIABLES AND DYNAMIC SETS
*-------------------------------------------------------------------------------

* determine the first hour of the day
n1(n) $[ord(n) = 1] = yes ;

* assignment of thermal units, storage hydro and pumped storage hydro plants
t (g) $[    pThermalGen(g,'FuelCost') *   pThermalGen(g,'MaxProd')] = yes ;
h (g) $[not pThermalGen(g,'FuelCost') and pHydroGen  (g,'MaxProd')] = yes ;

*determine the connection between the generator and the bus
gi (t,i) $[ pThermalGen(t,'Bus') and [pThermalGen(t,'Bus')=ORD(i)]] = yes ;
gi (h,i) $[ pHydroGen  (h,'Bus') and [pHydroGen  (h,'Bus')=ORD(i)]] = yes ;


* alias for dynamic iic
kkc(iic) = yes ;

IF ( [pNetworkConst = 1],
*     lines with "enable contingency" flag activated
      jjc(iic) $[ pNetwork(iic,'EnableCont') ] = yes ;
      cont(sc,n,iic,jjc) = yes ;
ELSE
      iic (i,ii,c )      = no ;
      kkc (i,ii,c )      = no ;
      jjc (j,jj,cc)      = no ;
      cont(sc,n,iic,jjc) = no ;
);

* scaling of parameters to GW and M€
pDemand        (n     )  = pDemand        (n     )         * 1e-3 ;
pDemand2ndResUp(n     )  = pDemand2ndResUp(n     )         * 1e-3 ;
pDemand2ndResDw(n     )  = pDemand2ndResDw(n     )         * 1e-3 ;
pIntermGen     (n,i,sc)  = pIntermGen     (n,i,sc)         * 1e-3 ;

pENSCost               = pENSCost                          * 1e-3 ;
pMaxProd     (t)       = pThermalGen (t,'MaxProd'        ) * 1e-3 ;
pMinProd     (t)       = pThermalGen (t,'MinProd'        ) * 1e-3 ;
pIniOut      (t)       = pThermalGen (t,'IniProd'        ) * 1e-3 ;
pRampUp      (t)       = pThermalGen (t,'RampUp'         ) * 1e-3 ;
pRampDw      (t)       = pThermalGen (t,'RampDw'         ) * 1e-3 ;
pMinTU       (t)       = pThermalGen (t,'MinTU'          ) ;
pMinTD       (t)       = pThermalGen (t,'MinTD'          ) ;
pSlopeVarCost(t)       = pThermalGen (t,'OMVarCost'      ) * 1e-3 +
                         pThermalGen (t,'SlopeVarCost'   ) * 1e-3 * pThermalGen(t,'FuelCost') ;
pDurationSD  (t)       = pThermalGen (t,'SDduration'     ) ;
pInterVarCost(t)       = pThermalGen (t,'InterVarCost'   ) * 1e-6 * pThermalGen(t,'FuelCost') ;
pShutdownCost(t)       = pThermalGen (t,'ShutdownCost'   ) * 1e-6 * pThermalGen(t,'FuelCost') ;
pDurationSU  (t,'su1') = pThermalGen (t,'SUduration1'    ) ;
pDurationSU  (t,'su2') = pThermalGen (t,'SUduration2'    ) ;
pDurationSU  (t,'su3') = pThermalGen (t,'SUduration3'    ) ;
pDurationSU  (t,'su4') = pThermalGen (t,'SUduration4'    ) ;
pDurationSU  (t,'su5') = pThermalGen (t,'SUduration5'    ) ;
pTMinSU      (t,'su1') = pThermalGen (t,'DownTtimeforSU1') ;
pTMinSU      (t,'su2') = pThermalGen (t,'DownTtimeforSU2') ;
pTMinSU      (t,'su3') = pThermalGen (t,'DownTtimeforSU3') ;
pTMinSU      (t,'su4') = pThermalGen (t,'DownTtimeforSU4') ;
pTMinSU      (t,'su5') = pThermalGen (t,'DownTtimeforSU5') ;
pTypeSUCost  (t,'su1') = pThermalGen (t,'SUcost1'        ) * 1e-6 * pThermalGen(t,'FuelCost') ;
pTypeSUCost  (t,'su2') = pThermalGen (t,'SUcost2'        ) * 1e-6 * pThermalGen(t,'FuelCost') ;
pTypeSUCost  (t,'su3') = pThermalGen (t,'SUcost3'        ) * 1e-6 * pThermalGen(t,'FuelCost') ;
pTypeSUCost  (t,'su4') = pThermalGen (t,'SUcost4'        ) * 1e-6 * pThermalGen(t,'FuelCost') ;
pTypeSUCost  (t,'su5') = pThermalGen (t,'SUcost5'        ) * 1e-6 * pThermalGen(t,'FuelCost') ;
pNumSUTypes  (t      ) = sum[l $pTypeSUCost(t,l), 1]       ;

pMaxProd     (h)       = pHydroGen   (h,'MaxProd'        ) * 1e-3 ;
pMinProd     (h)       = pHydroGen   (h,'MinProd'        ) * 1e-3 ;
pMaxCons     (h)       = pHydroGen   (h,'MaxCons'        ) * 1e-3 ;
pTotalEnergy (h)       = pHydroGen   (h,'TotEnergy'      ) * 1e-3 ;
pMax2ndResUP (h)       = pHydroGen   (h,'Max2ndResUP'    ) * 1e-3 ;
pMin2ndResUP (h)       = pHydroGen   (h,'Min2ndResUP'    ) * 1e-3 ;
pTot2ndResUP (h)       = pHydroGen   (h,'Tot2ndResUP'    ) * 1e-3 ;
pMax2ndResDW (h)       = pHydroGen   (h,'Max2ndResDW'    ) * 1e-3 ;
pMin2ndResDW (h)       = pHydroGen   (h,'Min2ndResDW'    ) * 1e-3 ;
pTot2ndResDW (h)       = pHydroGen   (h,'Tot2ndResDW'    ) * 1e-3 ;

*Network parameters
pPmax   (iic)  = pNetwork    (iic,'Pmax'      ) * 1e-3  ;
pPctmax (iic)  = pNetwork    (iic,'Pctmax'    ) * 1e-3  ;
pRline  (iic)  = pNetwork    (iic,'R'         )         ;
pXline  (iic)  = pNetwork    (iic,'X'         )         ;
pZline  (iic)  = SQRT(SQR(pRline(iic))+SQR(pXline(iic)));
pSBase         = pSBase                         * 1e-3  ;

pPTDF   (iic(i,ii,c),j,jj) = pISF(iic,j)-pISF(iic,jj);

pLODF   (iic(i,ii,c),kkc(j,jj,cc)) $[ not pPTDF(kkc,j,jj) = 1 ]
                           = round(pPTDF(iic,j,jj)/[1-pPTDF(kkc,j,jj)],6);

* if the initial output of the unit is above its minimum load then the unit is committed, otherwise it is not committed
pIniUC(g)   = 1 $[pIniOut(g) >= pMinProd(g)] ;

* if the the minimun up or down times are 0, they are changed to 1
pMinTU(t) $[pMinTU(t) = 0] = 1 ;
pMinTD(t) $[pMinTD(t) = 0] = 1 ;

* determine the ramp output
pProductRampDw(t,n  ) $[pDurationSD(t  )*pDemand(n)] = (2*ord(n)-1)/(2*pDurationSD(t  ))*pMinProd(t) ;
pProductRampUp(t,n,l) $[pDurationSU(t,l)*pDemand(n)] = (2*ord(n)-1)/(2*pDurationSU(t,l))*pMinProd(t) ;

*-------------------------------------------------------------------------------
*                        SET BOUNDS FOR VARIABLES
*-------------------------------------------------------------------------------

* bounds on variables
vProduct.up (sc,n,g) $[pScenProb(sc)*pDemand(n)] = pMaxProd    (g     ) ;
vProduct1.up(sc,n,g) $[pScenProb(sc)*pDemand(n)] = pMaxProd    (g     ) -
                                                   PMinProd    (g     ) ;
vProduct1.up(sc,n,g) $[pScenProb(sc)*pDemand(n) and pModelForm = 2] = 1 ;
vIG.up      (sc,n,i) $[pScenProb(sc)*pDemand(n)] = pIntermGen  (n,i,sc) ;
vConsump.up (sc,n,g) $[pScenProb(sc)*pDemand(n)] = pMaxCons    (g     ) ;
v2ndResUP.up(   n,h) $               pDemand(n)  = pMax2ndResUP(h     ) ;
v2ndResUP.lo(   n,h) $               pDemand(n)  = pMin2ndResUP(h     ) ;
v2ndResDW.up(   n,h) $               pDemand(n)  = pMax2ndResDW(h     ) ;
v2ndResDW.lo(   n,h) $               pDemand(n)  = pMin2ndResDW(h     ) ;
v2ndResUP.up(   n,t) $               pDemand(n)  = pMaxProd    (t     ) -
                                                   PMinProd    (t     ) ;
v2ndResUP.up(   n,t) $[              pDemand(n) and pModelForm = 2] = 1 ;
v2ndResDW.up(   n,t) $               pDemand(n)  = pMaxProd    (t     ) -
                                                   PMinProd    (t     ) ;
v2ndResDW.up(   n,t) $[              pDemand(n) and pModelForm = 2] = 1 ;
vENS.up     (sc,n,i) $[pScenProb(sc)*pDemand(n)] = pDemand     (n     ) *
                                                   pBusProfile (i     ) ;
vTypeSU.fx  ( n,t,l) $[not pTypeSUCost(t,l)    ] = 0                    ;

* bounds for circuit power flow in steady state
vCirPF.up    (sc,n,iic    )  $[pScenProb(sc)*pDemand(n)] =  pPmax  (iic);
vCirPF.lo    (sc,n,iic    )  $[pScenProb(sc)*pDemand(n)] = -pPmax  (iic);
* bounds for circuit power flow in contingency state
vCirPFCont.up(sc,n,iic,kkc)  $[pScenProb(sc)*pDemand(n)] =  pPctmax(iic);
vCirPFCont.lo(sc,n,iic,kkc)  $[pScenProb(sc)*pDemand(n)] = -pPctmax(iic);
* if no network constraints then it's only necessary one vENS
vENS.fx      (sc,n,is(i)  )  $[pNetworkConst=0         ] =  0           ;

*-------------------------------------------------------------------------------
*                       ITERATIVE METHOD FOR SCUC
*-------------------------------------------------------------------------------

* initial parameters for the first iteration
aux      = 1 ;
iter     = 0 ;
pCPUTime = 0 ;

WHILE(aux > 0,
         // 1. Solve the model
         solve SDUC using MIP minimizing vTotalVCost ;

         // 2. Estimation of Contingency Power Flow using LODF
         pCirPFContAux1(sc,n,iic(i,ii,c),kkc(j,jj,cc)) $[pScenProb(sc)*pDemand(n) and
                         not(sameas(i,j) and sameas(ii,jj) and sameas(c,cc))
                         and (ABS(pPTDF(iic,j,jj)) gt 1e-6)]
                         = ABS[vCirPF.l(sc,n,iic) + pLODF(iic,kkc)*vCirPF.l(sc,n,kkc)] ;

         // 3. Calculation lines out of limits.
         pCirPFContAux2(sc,n,iic,kkc) $[ pCirPFContAux1(sc,n,iic,kkc) gt pFCL*pPctmax(iic)] = 1;

         // 4. Update of subsets of contingencies
         jjc (j,jj,cc) $[SUM[(sc,n,iic), pCirPFContAux2(sc,n,iic,j,jj,cc)]]= yes ;
         cont(sc,n,iic,jjc)            $ pCirPFContAux2(sc,n,iic,jjc)      = yes ;

         // 5. Calculating the number of contingencies added at this iteration.
         //    The iterations will stop when no more lines are included (aux=0)
         aux = SUM[cont(sc,n,iic,jjc), pCirPFContAux2(sc,n,iic,jjc)];
         aux $[not pEnableIterAlg] = 0;

         // 6. Scalars to determine the number of the iteration and CPUtime
         iter = iter+1;
         pCPUTime = pCPUTime + SDUC.resUsd;

         pSummary('Final Obj. Function') = vTotalVCost.l              ;
         pSummary('Total Iterations   ') $[    pEnableIterAlg] = iter ;
         pSummary('Total Iterations   ') $[not pEnableIterAlg] =    0 ;
         pSummary('Total CPUTime   [s]') = pCPUTime                   ;
         pSummary('Total N-1 Added    ') = aux      + pSummary('Total N-1 Added    ') ;
         pSummary('UC Formulation Used') = pModelForm                 ;

         // 7. Reset of auxiliary parameters for the next iteration
         pCirPFContAux1(sc,n,iic,kkc)=0;
         pCirPFContAux2(sc,n,iic,kkc)=0;
);

*-------------------------------------------------------------------------------
*                            CALCULATE RESULTS
*-------------------------------------------------------------------------------

* scaling of results
if(pModelForm = 0 or pModelForm = 2,
    pCommitt (t,n)$ pDemand(n)  = vCommitt.l(n,t)                   + eps ;
elseif pModelForm = 1,
    pCommitt (t,n)$ pDemand(n)  = vCommitt.l(n,t) + vStartUp.l(n,t) + eps ;
);

pProduct (sc,g,n) $[pScenProb(sc) *  pDemand(n)] = vProduct.l   (sc,n,g)*1e3 + eps ;
pConsump (sc,h,n) $[pScenProb(sc) *  pDemand(n)] = vConsump.l   (sc,n,h)*1e3 + eps ;
pENS     (sc,n,i) $[pScenProb(sc) *  pDemand(n)] = vENS.l       (sc,n,i)*1e3 + eps ;
pSRMC    (sc,n,i) $[pScenProb(sc) *  pDemand(n)
                       and pNetworkConst = 1   ] = eBalanceBus.m(sc,n,i)*1e3 + eps ;
pSRMC    (sc,n,i) $[pScenProb(sc) *  pDemand(n)
                       and pNetworkConst = 0   ] = eBalance.m   (sc,n  )*1e3 + eps ;
p2ResUpMC(     n) $                  pDemand(n)  = e2ReserveUp.m(   n  )*1e3 + eps ;
p2ResDwMC(     n) $                  pDemand(n)  = e2ReserveDw.m(   n  )*1e3 + eps ;
p2ndResUP(   g,n) $[pModelForm<2 and pDemand(n)] = v2ndResUP.l  (   n,g)*1e3 + eps ;
p2ndResDW(   g,n) $[pModelForm<2 and pDemand(n)] = v2ndResDW.l  (   n,g)*1e3 + eps ;
p2ndResUP(   g,n) $[pModelForm=2 and pDemand(n)] = v2ndResUP.l  (   n,g)*1e3 * (pMaxProd(g)-pMinProd(g)) + eps ;
p2ndResDW(   g,n) $[pModelForm=2 and pDemand(n)] = v2ndResDW.l  (   n,g)*1e3 * (pMaxProd(g)-pMinProd(g)) + eps ;

pIG       (sc,n,i  ) $[pScenProb(sc)*pDemand(n) and pIntermGen(n,i,sc)]
                                         = vIG.l       (sc,n,i      )*1e3 + eps ;

*Parameters for only print the lines that are congested in contingency or steady state
pPrintCirPF(sc,iic) $[sum(n,sum(jjc,abs(vCirPFCont.m(sc,n,iic,jjc))>0))
                      or    sum(n  ,abs(vCirPF.m    (sc,n,iic    ))>0)] = 1 ;
pPrintCirPF(sc,iic) $[not pEnableIterAlg] = 1 ;
pPrintCirPFCont(sc,iic,jjc) $[sum(n,abs(vCirPFCont.m(sc,n,iic,jjc))>0)] = 1 ;

pCirPF    (sc,n,iic    ) $[pScenProb(sc)*pDemand(n) and
                                 pPrintCirPF(sc,iic)]
                                         = vCirPF.l    (sc,n,iic    )*1e3 + eps ;

pCirPFCont(sc,iic,n,jjc) $[pScenProb(sc)*pDemand(n) and
                                 pPrintCirPFCont(sc,iic,jjc)]
                                         = vCirPFCont.l(sc,n,iic,jjc)*1e3 + eps ;

* data output to xls file
put TMP put 'par=pCommitt rdim=1 rng=UC!a1' / 'par=pProduct rdim=2 rng=Output!a1' / 'par=pConsump rdim=2 rng=Consumption!a1' /'par=pIG rdim=2 rng=IG!a1' /
put TMP put 'par=pSRMC rdim=2 rng=SRMC!a1' / 'par=p2ndResUP rdim=1 rng=SecResUP!a1' / 'par=p2ndResDW rdim=1 rng=SecResDW!a1' / 'par=p2ResUpMC rdim=1 rng=SecResMC!a1' /
put TMP put 'par=p2ResDwMC rdim=1 rng=SecResMC!c1' / 'par=pCirPF rdim=2 rng=CirPF!a1' / 'par=pCirPFCont rdim=5 rng=CirPFCont!a1' / 'par=pSummary rdim=1 rng=Summary!a1' /
put TMP put 'par=pENS      rdim=2 rng=ENS!a1' /
putclose
execute_unload   'tmp.gdx' pProduct pCommitt pConsump pIG pSRMC p2ndResUP p2ndResDW p2ResUpMC p2ResDwMC pCirPF pCirPFCont pSummary pENS
execute          'gdxxrw.exe tmp.gdx SQ=n EpsOut=0 O="tmp.xlsx" @tmp.txt'
execute          'del        tmp.gdx                             tmp.txt' ;

$onlisting
