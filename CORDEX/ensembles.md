ensembles
================
Rasmus Benestad
April 27, 2016

R Markdown
----------

This is an R Markdown document. Markdown is a simple formatting syntax for authoring HTML, PDF, and MS Word documents. For more details on using R Markdown see <http://rmarkdown.rstudio.com>.

When you click the **Knit** button a document will be generated that includes both content as well as the output of any embedded R code chunks within the document. You can embed an R code chunk like this:

**Simulated global mean temperartue**
-------------------------------------

Combine the simulated global mean temperature with information about the model configurations taken from Table 9-a-1 in the IPCC WG 1 2013 report ("AR5):

``` r
library('esd')
data("IPCC.AR5.Table.9.A.1")
data("global.t2m.gcm")
gcmnm <- unlist(lapply(attr(global.t2m.gcm$global.t2m.cmip5.rcp45,'meta'),function(x) x$GCM))
n <- length(gcmnm); m <- 12
X <- matrix(rep(NA,n*m),n,m)

IPCC.AR5.Table.9.A.1$Atmosphere.Component.Name <- as.character(IPCC.AR5.Table.9.A.1$Atmosphere.Component.Name)
IPCC.AR5.Table.9.A.1$Ocean.Component.Name <- as.character(IPCC.AR5.Table.9.A.1$Ocean.Component.Name)
```

Index some of the entries which do not contain useful information.

``` r
iam <- is.element(substr(IPCC.AR5.Table.9.A.1$Atmosphere.Component.Name,1,8),'Included') |
       is.element(IPCC.AR5.Table.9.A.1$Atmosphere.Component.Name,'<NA>')
IPCC.AR5.Table.9.A.1$Atmosphere.Component.Name[iam] <- IPCC.AR5.Table.9.A.1$Model.Name[iam]
iom <- is.element(IPCC.AR5.Table.9.A.1$Ocean.Component.Name,'Included')
IPCC.AR5.Table.9.A.1$Ocean.Component.Name[iom] <- IPCC.AR5.Table.9.A.1$Model.Name[iom]
```

There are many different similar settings with different names. Try to organise and sort the different types under fewer and simpler categories. Here is a mix of conventions and ways to refer to the resolution.

``` r
fixagcmnm <- grep('Included',as.character(IPCC.AR5.Table.9.A.1$Atmosphere.Component.Name))
IPCC.AR5.Table.9.A.1$Atmosphere.Component.Name[fixagcmnm] <- IPCC.AR5.Table.9.A.1$Model.Name[fixagcmnm]
```

``` r
IPCC.AR5.Table.9.A.1$Horizontal.Grid <- as.character(IPCC.AR5.Table.9.A.1$Horizontal.Grid)
IPCC.AR5.Table.9.A.1$Horizontal.Grid <- gsub(' ','',IPCC.AR5.Table.9.A.1$Horizontal.Grid)
IPCC.AR5.Table.9.A.1$Horizontal.Grid <- gsub('-','',IPCC.AR5.Table.9.A.1$Horizontal.Grid)
IPCC.AR5.Table.9.A.1$Horizontal.Grid <- gsub(',','',IPCC.AR5.Table.9.A.1$Horizontal.Grid)
IPCC.AR5.Table.9.A.1$Horizontal.Grid <- gsub('longitude','',IPCC.AR5.Table.9.A.1$Horizontal.Grid)
IPCC.AR5.Table.9.A.1$Horizontal.Grid <- gsub('latitude','',IPCC.AR5.Table.9.A.1$Horizontal.Grid)
IPCC.AR5.Table.9.A.1$Horizontal.Grid <- gsub('Nominally','',IPCC.AR5.Table.9.A.1$Horizontal.Grid)
IPCC.AR5.Table.9.A.1$Horizontal.Grid <- gsub('equivalentto','',IPCC.AR5.Table.9.A.1$Horizontal.Grid)
IPCC.AR5.Table.9.A.1$Horizontal.Grid <- gsub('FiniteVolume','',IPCC.AR5.Table.9.A.1$Horizontal.Grid)
IPCC.AR5.Table.9.A.1$Horizontal.Grid <- gsub('Averagedcellsize','',IPCC.AR5.Table.9.A.1$Horizontal.Grid)
IPCC.AR5.Table.9.A.1$Horizontal.Grid <- gsub('in','',IPCC.AR5.Table.9.A.1$Horizontal.Grid)
IPCC.AR5.Table.9.A.1$Horizontal.Grid <- gsub('and','',IPCC.AR5.Table.9.A.1$Horizontal.Grid)
IPCC.AR5.Table.9.A.1$Horizontal.Grid <- gsub('TL','T',IPCC.AR5.Table.9.A.1$Horizontal.Grid)
Tres<- c('T42','T31','T159','T63','T159','T959','T213','T319','N48','N49','C180','C360','R42','N96')
for (ires in Tres) {
  ii <- grep(ires,IPCC.AR5.Table.9.A.1$Horizontal.Grid)
  IPCC.AR5.Table.9.A.1$Horizontal.Grid[ii] <- ires
}
IPCC.AR5.Table.9.A.1$Horizontal.Grid <- as.factor(IPCC.AR5.Table.9.A.1$Horizontal.Grid)
```

Fix the GCM names:

``` r
IPCC.AR5.Table.9.A.1$Model.Name <- as.character(IPCC.AR5.Table.9.A.1$Model.Name)
IPCC.AR5.Table.9.A.1$Model.Name <- gsub('-','.',IPCC.AR5.Table.9.A.1$Model.Name)
IPCC.AR5.Table.9.A.1$Model.Name <- gsub('_','.',IPCC.AR5.Table.9.A.1$Model.Name)
IPCC.AR5.Table.9.A.1$Model.Name <- gsub(' v1.0','',IPCC.AR5.Table.9.A.1$Model.Name)
IPCC.AR5.Table.9.A.1$Model.Name <- gsub('(m)','.m',IPCC.AR5.Table.9.A.1$Model.Name,perl=TRUE,fixed=TRUE)
```

    ## Warning in gsub("(m)", ".m", IPCC.AR5.Table.9.A.1$Model.Name, perl =
    ## TRUE, : argument 'perl = TRUE' will be ignored

``` r
cesm1 <- grep('CESM1',IPCC.AR5.Table.9.A.1$Model.Name)
IPCC.AR5.Table.9.A.1$Model.Name[cesm1] <- 'CESM1'
IPCC.AR5.Table.9.A.1$Model.Name <- as.factor(IPCC.AR5.Table.9.A.1$Model.Name)
gcmnm <- gsub('-','.',gcmnm)
gcmnm <- gsub('_','.',gcmnm)
cesm1 <- grep('CESM1',gcmnm)
gcmnm[cesm1] <- 'CESM1'
cnrm5 <- grep('CNRM.CM5',gcmnm)
gcmnm[cnrm5] <- 'CNRM.CM51'
inmcm4 <- grep('inmcm4',gcmnm)
gcmnm[inmcm4] <- 'INM.CM4'
```

Also need to tidy up some of the description concerning the ocean models

``` r
IPCC.AR5.Table.9.A.1$z.Co.ord <- as.character(IPCC.AR5.Table.9.A.1$z.Co.ord)
IPCC.AR5.Table.9.A.1$z.Co.ord <- gsub('-coordinate','',IPCC.AR5.Table.9.A.1$z.Co.ord)
IPCC.AR5.Table.9.A.1$z.Co.ord <- gsub('hybrid','',IPCC.AR5.Table.9.A.1$z.Co.ord)
IPCC.AR5.Table.9.A.1$z.Co.ord <- gsub('sigma-z','sigma',IPCC.AR5.Table.9.A.1$z.Co.ord)
IPCC.AR5.Table.9.A.1$z.Co.ord <- gsub(' ','',IPCC.AR5.Table.9.A.1$z.Co.ord)
IPCC.AR5.Table.9.A.1$z.Co.ord <- gsub('Isopycnic','isopycnic',IPCC.AR5.Table.9.A.1$z.Co.ord)
IPCC.AR5.Table.9.A.1$z.Co.ord <- gsub('zisopycnic','isopycnic',IPCC.AR5.Table.9.A.1$z.Co.ord)
IPCC.AR5.Table.9.A.1$z.Co.ord <- as.factor(IPCC.AR5.Table.9.A.1$z.Co.ord)
```

Make a data.frame that organises the model information for each simulation.

    ##  [1] "ACCESS1.0"          "N96"                "38"                
    ##  [4] "13"                 "CLASSIC"            "Not implemented"   
    ##  [7] "MOSES2.2"           "ACCESS-OM (MOM4p1)" "50"                
    ## [10] "z*"                 "Not implemented"    "CICE4.1"           
    ##  [1] "ACCESS1.3"          "N96"                "38"                
    ##  [4] "13"                 "CLASSIC"            "Not implemented"   
    ##  [7] "CABLE"              "ACCESS-OM (MOM4p1)" "50"                
    ## [10] "z*"                 "Not implemented"    "CICE4.1"           
    ##  [1] "BCC_AGCM2.1"                  "T42"                         
    ##  [3] "26"                           "2.917"                       
    ##  [5] "Prescribed"                   "Not implemented"             
    ##  [7] "BCC-AVIM1.0"                  "MOM4-L40"                    
    ##  [9] "40"                           "z"                           
    ## [11] "Included"                     "GFDL Sea Ice Simulator (SIS)"
    ##  [1] "BCC_AGCM2.1"                  "T106"                        
    ##  [3] "26"                           "2.917"                       
    ##  [5] "Prescribed"                   "Not implemented"             
    ##  [7] "BCC-AVIM1.0"                  "MOM4-L40"                    
    ##  [9] "40"                           "z"                           
    ## [11] "Included"                     "GFDL Sea Ice Simulator (SIS)"
    ##  [1] "CAM3.5"            "T42"               "26"               
    ##  [4] "2.194"             "Semi-interactive"  "Not implemented"  
    ##  [7] "CoLM+BNUDGVM(C/N)" "MOM4p1"            "50"               
    ## [10] NA                  "IBGC"              "CICE4.1"          
    ##  [1] "CanESM2"         "T63"             "35"             
    ##  [4] "0.5"             "Interactive"     "Included"       
    ##  [7] "CLASS 2.7, CTEM" "CanESM2"         "40"             
    ## [10] "z"               "CMOC"            "Included"       
    ##  [1] "CanESM2"         "T63"             "35"             
    ##  [4] "0.5"             "Interactive"     "Included"       
    ##  [7] "CLASS 2.7, CTEM" "CanESM2"         "40"             
    ## [10] "z"               "CMOC"            "Included"       
    ##  [1] "CanESM2"         "T63"             "35"             
    ##  [4] "0.5"             "Interactive"     "Included"       
    ##  [7] "CLASS 2.7, CTEM" "CanESM2"         "40"             
    ## [10] "z"               "CMOC"            "Included"       
    ##  [1] "CanESM2"         "T63"             "35"             
    ##  [4] "0.5"             "Interactive"     "Included"       
    ##  [7] "CLASS 2.7, CTEM" "CanESM2"         "40"             
    ## [10] "z"               "CMOC"            "Included"       
    ##  [1] "CanESM2"         "T63"             "35"             
    ##  [4] "0.5"             "Interactive"     "Included"       
    ##  [7] "CLASS 2.7, CTEM" "CanESM2"         "40"             
    ## [10] "z"               "CMOC"            "Included"       
    ##  [1] "CAM4"                          "0.9x1.25deg"                  
    ##  [3] "27"                            "2.194067"                     
    ##  [5] "Interactive"                   "Not implemented"              
    ##  [7] "Community Land Model 4 (CLM4)" "POP2 with modifications"      
    ##  [9] "60"                            "z"                            
    ## [11] "Not implemented"               "CICE4 with modifications"     
    ##  [1] "CAM4"                          "0.9x1.25deg"                  
    ##  [3] "27"                            "2.194067"                     
    ##  [5] "Interactive"                   "Not implemented"              
    ##  [7] "Community Land Model 4 (CLM4)" "POP2 with modifications"      
    ##  [9] "60"                            "z"                            
    ## [11] "Not implemented"               "CICE4 with modifications"     
    ##  [1] "CAM4"                          "0.9x1.25deg"                  
    ##  [3] "27"                            "2.194067"                     
    ##  [5] "Interactive"                   "Not implemented"              
    ##  [7] "Community Land Model 4 (CLM4)" "POP2 with modifications"      
    ##  [9] "60"                            "z"                            
    ## [11] "Not implemented"               "CICE4 with modifications"     
    ##  [1] "CAM4"                          "0.9x1.25deg"                  
    ##  [3] "27"                            "2.194067"                     
    ##  [5] "Interactive"                   "Not implemented"              
    ##  [7] "Community Land Model 4 (CLM4)" "POP2 with modifications"      
    ##  [9] "60"                            "z"                            
    ## [11] "Not implemented"               "CICE4 with modifications"     
    ##  [1] "CAM4"                          "0.9x1.25deg"                  
    ##  [3] "27"                            "2.194067"                     
    ##  [5] "Interactive"                   "Not implemented"              
    ##  [7] "Community Land Model 4 (CLM4)" "POP2 with modifications"      
    ##  [9] "60"                            "z"                            
    ## [11] "Not implemented"               "CICE4 with modifications"     
    ##  [1] "CAM4"                          "0.9x1.25deg"                  
    ##  [3] "27"                            "2.194067"                     
    ##  [5] "Interactive"                   "Not implemented"              
    ##  [7] "Community Land Model 4 (CLM4)" "POP2 with modifications"      
    ##  [9] "60"                            "z"                            
    ## [11] "Not implemented"               "CICE4 with modifications"     
    ##  [1] "CAM4"                                  
    ##  [2] "0.9x1.25deg"                           
    ##  [3] "27"                                    
    ##  [4] "2.194067"                              
    ##  [5] "Semi-interactive"                      
    ##  [6] "Not implemented"                       
    ##  [7] "CLM4"                                  
    ##  [8] "POP2 with modifications"               
    ##  [9] "60"                                    
    ## [10] "z"                                     
    ## [11] "Biogeochemical Elemental Cycling (BEC)"
    ## [12] "CICE4 with modifications"              
    ##  [1] "CAM4"                                  
    ##  [2] "0.9x1.25deg"                           
    ##  [3] "27"                                    
    ##  [4] "2.194067"                              
    ##  [5] "Semi-interactive"                      
    ##  [6] "Not implemented"                       
    ##  [7] "CLM4"                                  
    ##  [8] "POP2 with modifications"               
    ##  [9] "60"                                    
    ## [10] "z"                                     
    ## [11] "Biogeochemical Elemental Cycling (BEC)"
    ## [12] "CICE4 with modifications"              
    ##  [1] "CAM4"                                  
    ##  [2] "0.9x1.25deg"                           
    ##  [3] "27"                                    
    ##  [4] "2.194067"                              
    ##  [5] "Semi-interactive"                      
    ##  [6] "Not implemented"                       
    ##  [7] "CLM4"                                  
    ##  [8] "POP2 with modifications"               
    ##  [9] "60"                                    
    ## [10] "z"                                     
    ## [11] "Biogeochemical Elemental Cycling (BEC)"
    ## [12] "CICE4 with modifications"              
    ##  [1] "CAM4"                                  
    ##  [2] "0.9x1.25deg"                           
    ##  [3] "27"                                    
    ##  [4] "2.194067"                              
    ##  [5] "Semi-interactive"                      
    ##  [6] "Not implemented"                       
    ##  [7] "CLM4"                                  
    ##  [8] "POP2 with modifications"               
    ##  [9] "60"                                    
    ## [10] "z"                                     
    ## [11] "Biogeochemical Elemental Cycling (BEC)"
    ## [12] "CICE4 with modifications"              
    ##  [1] "ECHAM5"           "T159"             "31"              
    ##  [4] "10"               "Semi-interactive" "Not implemented" 
    ##  [7] "Not implemented"  "OPA8.2"           "31"              
    ## [10] "z"                "Not implemented"  "LIM2"            
    ##  [1] "ECHAM5"           "T63"              "95"              
    ##  [4] "0.01"             "Semi-interactive" "Not implemented" 
    ##  [7] "Not implemented"  "OPA8.2"           "31"              
    ## [10] "z"                "Not implemented"  "LIM2"            
    ##  [1] "ARPEGE-Climat"                     
    ##  [2] "T127"                              
    ##  [3] "31"                                
    ##  [4] "10"                                
    ##  [5] "Prescribed"                        
    ##  [6] "(3-D linear ozone chemistry model)"
    ##  [7] "SURFEX (Land and Ocean Surface)"   
    ##  [8] "NEMO"                              
    ##  [9] "42"                                
    ## [10] "z"                                 
    ## [11] "PISCES"                            
    ## [12] "Gelato5 (Sea Ice)"                 
    ##  [1] "CSIRO-Mk3.6.0"   "T63"             "18"             
    ##  [4] "4.5"             "Interactive"     "Not implemented"
    ##  [7] "Included"        "Modified MOM2.2" "31"             
    ## [10] "z"               "Not implemented" "Included"       
    ##  [1] "CSIRO-Mk3.6.0"   "T63"             "18"             
    ##  [4] "4.5"             "Interactive"     "Not implemented"
    ##  [7] "Included"        "Modified MOM2.2" "31"             
    ## [10] "z"               "Not implemented" "Included"       
    ##  [1] "CSIRO-Mk3.6.0"   "T63"             "18"             
    ##  [4] "4.5"             "Interactive"     "Not implemented"
    ##  [7] "Included"        "Modified MOM2.2" "31"             
    ## [10] "z"               "Not implemented" "Included"       
    ##  [1] "CSIRO-Mk3.6.0"   "T63"             "18"             
    ##  [4] "4.5"             "Interactive"     "Not implemented"
    ##  [7] "Included"        "Modified MOM2.2" "31"             
    ## [10] "z"               "Not implemented" "Included"       
    ##  [1] "CSIRO-Mk3.6.0"   "T63"             "18"             
    ##  [4] "4.5"             "Interactive"     "Not implemented"
    ##  [7] "Included"        "Modified MOM2.2" "31"             
    ## [10] "z"               "Not implemented" "Included"       
    ##  [1] "CSIRO-Mk3.6.0"   "T63"             "18"             
    ##  [4] "4.5"             "Interactive"     "Not implemented"
    ##  [7] "Included"        "Modified MOM2.2" "31"             
    ## [10] "z"               "Not implemented" "Included"       
    ##  [1] "CSIRO-Mk3.6.0"   "T63"             "18"             
    ##  [4] "4.5"             "Interactive"     "Not implemented"
    ##  [7] "Included"        "Modified MOM2.2" "31"             
    ## [10] "z"               "Not implemented" "Included"       
    ##  [1] "CSIRO-Mk3.6.0"   "T63"             "18"             
    ##  [4] "4.5"             "Interactive"     "Not implemented"
    ##  [7] "Included"        "Modified MOM2.2" "31"             
    ## [10] "z"               "Not implemented" "Included"       
    ##  [1] "CSIRO-Mk3.6.0"   "T63"             "18"             
    ##  [4] "4.5"             "Interactive"     "Not implemented"
    ##  [7] "Included"        "Modified MOM2.2" "31"             
    ## [10] "z"               "Not implemented" "Included"       
    ##  [1] "CSIRO-Mk3.6.0"   "T63"             "18"             
    ##  [4] "4.5"             "Interactive"     "Not implemented"
    ##  [7] "Included"        "Modified MOM2.2" "31"             
    ## [10] "z"               "Not implemented" "Included"       
    ##  [1] "IFS c31r1"       "T159"            "62"             
    ##  [4] "1"               "Prescribed"      "Not implemented"
    ##  [7] "HTESSEL"         "NEMO_ecmwf"      "31"             
    ## [10] "z"               "Not implemented" "LIM2"           
    ##  [1] "IFS c31r1"       "T159"            "62"             
    ##  [4] "1"               "Prescribed"      "Not implemented"
    ##  [7] "HTESSEL"         "NEMO_ecmwf"      "31"             
    ## [10] "z"               "Not implemented" "LIM2"           
    ##  [1] "IFS c31r1"       "T159"            "62"             
    ##  [4] "1"               "Prescribed"      "Not implemented"
    ##  [7] "HTESSEL"         "NEMO_ecmwf"      "31"             
    ## [10] "z"               "Not implemented" "LIM2"           
    ##  [1] "IFS c31r1"       "T159"            "62"             
    ##  [4] "1"               "Prescribed"      "Not implemented"
    ##  [7] "HTESSEL"         "NEMO_ecmwf"      "31"             
    ## [10] "z"               "Not implemented" "LIM2"           
    ##  [1] "IFS c31r1"       "T159"            "62"             
    ##  [4] "1"               "Prescribed"      "Not implemented"
    ##  [7] "HTESSEL"         "NEMO_ecmwf"      "31"             
    ## [10] "z"               "Not implemented" "LIM2"           
    ##  [1] "IFS c31r1"       "T159"            "62"             
    ##  [4] "1"               "Prescribed"      "Not implemented"
    ##  [7] "HTESSEL"         "NEMO_ecmwf"      "31"             
    ## [10] "z"               "Not implemented" "LIM2"           
    ##  [1] "IFS c31r1"       "T159"            "62"             
    ##  [4] "1"               "Prescribed"      "Not implemented"
    ##  [7] "HTESSEL"         "NEMO_ecmwf"      "31"             
    ## [10] "z"               "Not implemented" "LIM2"           
    ##  [1] "GAMIL2"           "2.8125x2.8125deg" "26"              
    ##  [4] "2.194"            "Semi-interactive" "Not implemented" 
    ##  [7] "CLM3"             "LICOM2"           "30"              
    ## [10] "etaco-ordinate"   "Not implemented"  "CICE4-LASG"      
    ##  [1] "CAM3.0"                                                                            
    ##  [2] "T42"                                                                               
    ##  [3] "26"                                                                                
    ##  [4] "3.545"                                                                             
    ##  [5] "Prescribed"                                                                        
    ##  [6] "Not implemented"                                                                   
    ##  [7] "CLM3.5"                                                                            
    ##  [8] "Modified POP2.0 through incorporating the non-breaking surface wave-induced mixing"
    ##  [9] "40"                                                                                
    ## [10] "z"                                                                                 
    ## [11] "Improved OCMIP-2 biogeochemical model"                                             
    ## [12] "CICE4.0"                                                                           
    ##  [1] "CAM3.0"                                                                            
    ##  [2] "T42"                                                                               
    ##  [3] "26"                                                                                
    ##  [4] "3.545"                                                                             
    ##  [5] "Prescribed"                                                                        
    ##  [6] "Not implemented"                                                                   
    ##  [7] "CLM3.5"                                                                            
    ##  [8] "Modified POP2.0 through incorporating the non-breaking surface wave-induced mixing"
    ##  [9] "40"                                                                                
    ## [10] "z"                                                                                 
    ## [11] "Improved OCMIP-2 biogeochemical model"                                             
    ## [12] "CICE4.0"                                                                           
    ##  [1] "CAM3.0"                                                                            
    ##  [2] "T42"                                                                               
    ##  [3] "26"                                                                                
    ##  [4] "3.545"                                                                             
    ##  [5] "Prescribed"                                                                        
    ##  [6] "Not implemented"                                                                   
    ##  [7] "CLM3.5"                                                                            
    ##  [8] "Modified POP2.0 through incorporating the non-breaking surface wave-induced mixing"
    ##  [9] "40"                                                                                
    ## [10] "z"                                                                                 
    ## [11] "Improved OCMIP-2 biogeochemical model"                                             
    ## [12] "CICE4.0"                                                                           
    ##  [1] "GFDL-CM3"              "200kmC48L48"          
    ##  [3] "48"                    "0.01"                 
    ##  [5] "Interactive"           "Atmospheric Chemistry"
    ##  [7] "Included"              "MOM4.1"               
    ##  [9] "50"                    "z*"                   
    ## [11] "Not implemented"       "SIS"                  
    ##  [1] "GFDL-ESM2G"       "2.5deg2degM45L24" "24"              
    ##  [4] "3.65"             "Semi-interactive" "Not implemented" 
    ##  [7] "Included"         "GOLD"             "63"              
    ## [10] "isopycnic"        "TOPAZ"            "SIS"             
    ##  [1] "GFDL-ESM2M"       "2.5deg2degM45L24" "24"              
    ##  [4] "3.65"             "Semi-interactive" "Not implemented" 
    ##  [7] "Included"         "MOM4.1"           "50"              
    ## [10] "z*"               "TOPAZ"            "SIS"             
    ##  [1] "GISS-E2-H"       "2degx2.5deg"     "40"             
    ##  [4] "0.1"             "Interactive"     "G-PUCCINI"      
    ##  [7] "Included"        "HYCOM Ocean"     "26"             
    ## [10] "isopycnic"       "Not implemented" "Included"       
    ##  [1] "GISS-E2-H"       "2degx2.5deg"     "40"             
    ##  [4] "0.1"             "Interactive"     "G-PUCCINI"      
    ##  [7] "Included"        "HYCOM Ocean"     "26"             
    ## [10] "isopycnic"       "Not implemented" "Included"       
    ##  [1] "GISS-E2-H"       "2degx2.5deg"     "40"             
    ##  [4] "0.1"             "Interactive"     "G-PUCCINI"      
    ##  [7] "Included"        "HYCOM Ocean"     "26"             
    ## [10] "isopycnic"       "Not implemented" "Included"       
    ##  [1] "GISS-E2-H"       "2degx2.5deg"     "40"             
    ##  [4] "0.1"             "Interactive"     "G-PUCCINI"      
    ##  [7] "Included"        "HYCOM Ocean"     "26"             
    ## [10] "isopycnic"       "Not implemented" "Included"       
    ##  [1] "GISS-E2-H"       "2degx2.5deg"     "40"             
    ##  [4] "0.1"             "Interactive"     "G-PUCCINI"      
    ##  [7] "Included"        "HYCOM Ocean"     "26"             
    ## [10] "isopycnic"       "Not implemented" "Included"       
    ##  [1] "GISS-E2-H"       "2degx2.5deg"     "40"             
    ##  [4] "0.1"             "Interactive"     "G-PUCCINI"      
    ##  [7] "Included"        "HYCOM Ocean"     "26"             
    ## [10] "isopycnic"       "Not implemented" "Included"       
    ##  [1] "GISS-E2-H"       "2degx2.5deg"     "40"             
    ##  [4] "0.1"             "Interactive"     "G-PUCCINI"      
    ##  [7] "Included"        "HYCOM Ocean"     "26"             
    ## [10] "isopycnic"       "Not implemented" "Included"       
    ##  [1] "GISS-E2-H"       "2degx2.5deg"     "40"             
    ##  [4] "0.1"             "Interactive"     "G-PUCCINI"      
    ##  [7] "Included"        "HYCOM Ocean"     "26"             
    ## [10] "isopycnic"       "Not implemented" "Included"       
    ##  [1] "GISS-E2-H"       "2degx2.5deg"     "40"             
    ##  [4] "0.1"             "Interactive"     "G-PUCCINI"      
    ##  [7] "Included"        "HYCOM Ocean"     "26"             
    ## [10] "isopycnic"       "Not implemented" "Included"       
    ##  [1] "GISS-E2-H"       "2degx2.5deg"     "40"             
    ##  [4] "0.1"             "Interactive"     "G-PUCCINI"      
    ##  [7] "Included"        "HYCOM Ocean"     "26"             
    ## [10] "isopycnic"       "Not implemented" "Included"       
    ##  [1] "GISS-E2-H"       "2degx2.5deg"     "40"             
    ##  [4] "0.1"             "Interactive"     "G-PUCCINI"      
    ##  [7] "Included"        "HYCOM Ocean"     "26"             
    ## [10] "isopycnic"       "Not implemented" "Included"       
    ##  [1] "GISS-E2-H"       "2degx2.5deg"     "40"             
    ##  [4] "0.1"             "Interactive"     "G-PUCCINI"      
    ##  [7] "Included"        "HYCOM Ocean"     "26"             
    ## [10] "isopycnic"       "Not implemented" "Included"       
    ##  [1] "GISS-E2-H"       "2degx2.5deg"     "40"             
    ##  [4] "0.1"             "Interactive"     "G-PUCCINI"      
    ##  [7] "Included"        "HYCOM Ocean"     "26"             
    ## [10] "isopycnic"       "Not implemented" "Included"       
    ##  [1] "GISS-E2-H"       "2degx2.5deg"     "40"             
    ##  [4] "0.1"             "Interactive"     "G-PUCCINI"      
    ##  [7] "Included"        "HYCOM Ocean"     "26"             
    ## [10] "isopycnic"       "Not implemented" "Included"       
    ##  [1] "GISS-E2-H"       "2degx2.5deg"     "40"             
    ##  [4] "0.1"             "Interactive"     "G-PUCCINI"      
    ##  [7] "Included"        "HYCOM Ocean"     "26"             
    ## [10] "isopycnic"       "Not implemented" "Included"       
    ##  [1] "GISS-E2-H-CC"          "1deg"                 
    ##  [3] "40"                    "0.1"                  
    ##  [5] "Interactive (p1 only)" "G-PUCCINI"            
    ##  [7] "Included"              "HYCOM Ocean"          
    ##  [9] "26"                    "isopycnic"            
    ## [11] "Included"              "Included"             
    ##  [1] "GISS-E2-R"       "2degx2.5deg"     "40"             
    ##  [4] "0.1"             "Interactive"     "G-PUCCINI"      
    ##  [7] "Included"        "Russell Ocean"   "32"             
    ## [10] "z*"              "Not implemented" "Included"       
    ##  [1] "GISS-E2-R"       "2degx2.5deg"     "40"             
    ##  [4] "0.1"             "Interactive"     "G-PUCCINI"      
    ##  [7] "Included"        "Russell Ocean"   "32"             
    ## [10] "z*"              "Not implemented" "Included"       
    ##  [1] "GISS-E2-R"       "2degx2.5deg"     "40"             
    ##  [4] "0.1"             "Interactive"     "G-PUCCINI"      
    ##  [7] "Included"        "Russell Ocean"   "32"             
    ## [10] "z*"              "Not implemented" "Included"       
    ##  [1] "GISS-E2-R"       "2degx2.5deg"     "40"             
    ##  [4] "0.1"             "Interactive"     "G-PUCCINI"      
    ##  [7] "Included"        "Russell Ocean"   "32"             
    ## [10] "z*"              "Not implemented" "Included"       
    ##  [1] "GISS-E2-R"       "2degx2.5deg"     "40"             
    ##  [4] "0.1"             "Interactive"     "G-PUCCINI"      
    ##  [7] "Included"        "Russell Ocean"   "32"             
    ## [10] "z*"              "Not implemented" "Included"       
    ##  [1] "GISS-E2-R"       "2degx2.5deg"     "40"             
    ##  [4] "0.1"             "Interactive"     "G-PUCCINI"      
    ##  [7] "Included"        "Russell Ocean"   "32"             
    ## [10] "z*"              "Not implemented" "Included"       
    ##  [1] "GISS-E2-R"       "2degx2.5deg"     "40"             
    ##  [4] "0.1"             "Interactive"     "G-PUCCINI"      
    ##  [7] "Included"        "Russell Ocean"   "32"             
    ## [10] "z*"              "Not implemented" "Included"       
    ##  [1] "GISS-E2-R"       "2degx2.5deg"     "40"             
    ##  [4] "0.1"             "Interactive"     "G-PUCCINI"      
    ##  [7] "Included"        "Russell Ocean"   "32"             
    ## [10] "z*"              "Not implemented" "Included"       
    ##  [1] "GISS-E2-R"       "2degx2.5deg"     "40"             
    ##  [4] "0.1"             "Interactive"     "G-PUCCINI"      
    ##  [7] "Included"        "Russell Ocean"   "32"             
    ## [10] "z*"              "Not implemented" "Included"       
    ##  [1] "GISS-E2-R"       "2degx2.5deg"     "40"             
    ##  [4] "0.1"             "Interactive"     "G-PUCCINI"      
    ##  [7] "Included"        "Russell Ocean"   "32"             
    ## [10] "z*"              "Not implemented" "Included"       
    ##  [1] "GISS-E2-R"       "2degx2.5deg"     "40"             
    ##  [4] "0.1"             "Interactive"     "G-PUCCINI"      
    ##  [7] "Included"        "Russell Ocean"   "32"             
    ## [10] "z*"              "Not implemented" "Included"       
    ##  [1] "GISS-E2-R"       "2degx2.5deg"     "40"             
    ##  [4] "0.1"             "Interactive"     "G-PUCCINI"      
    ##  [7] "Included"        "Russell Ocean"   "32"             
    ## [10] "z*"              "Not implemented" "Included"       
    ##  [1] "GISS-E2-R"       "2degx2.5deg"     "40"             
    ##  [4] "0.1"             "Interactive"     "G-PUCCINI"      
    ##  [7] "Included"        "Russell Ocean"   "32"             
    ## [10] "z*"              "Not implemented" "Included"       
    ##  [1] "GISS-E2-R"       "2degx2.5deg"     "40"             
    ##  [4] "0.1"             "Interactive"     "G-PUCCINI"      
    ##  [7] "Included"        "Russell Ocean"   "32"             
    ## [10] "z*"              "Not implemented" "Included"       
    ##  [1] "GISS-E2-R"       "2degx2.5deg"     "40"             
    ##  [4] "0.1"             "Interactive"     "G-PUCCINI"      
    ##  [7] "Included"        "Russell Ocean"   "32"             
    ## [10] "z*"              "Not implemented" "Included"       
    ##  [1] "GISS-E2-R"       "2degx2.5deg"     "40"             
    ##  [4] "0.1"             "Interactive"     "G-PUCCINI"      
    ##  [7] "Included"        "Russell Ocean"   "32"             
    ## [10] "z*"              "Not implemented" "Included"       
    ##  [1] "GISS-E2-R"       "2degx2.5deg"     "40"             
    ##  [4] "0.1"             "Interactive"     "G-PUCCINI"      
    ##  [7] "Included"        "Russell Ocean"   "32"             
    ## [10] "z*"              "Not implemented" "Included"       
    ##  [1] "GISS-E2-R-CC"          "1deg"                 
    ##  [3] "40"                    "0.1"                  
    ##  [5] "Interactive (p1 only)" "G-PUCCINI"            
    ##  [7] "Included"              "Russell Ocean"        
    ##  [9] "32"                    "z*"                   
    ## [11] "Included"              "Included"             
    ##  [1] "HadGAM2"         "N96"             "60"             
    ##  [4] "0.091"           "Interactive"     "Not implemented"
    ##  [7] "Included"        "HadGEM2-AO"      NA               
    ## [10] "z"               "Not implemented" "Included"       
    ##  [1] "HadGAM2"               "N96"                  
    ##  [3] "60"                    "0.091"                
    ##  [5] "Interactive"           "Atmospheric Chemistry"
    ##  [7] "Included"              "HadGEM2-CC"           
    ##  [9] NA                      "z"                    
    ## [11] "Included"              "Included"             
    ##  [1] "HadGAM2"               "N96"                  
    ##  [3] "38"                    "13"                   
    ##  [5] "Interactive"           "Atmospheric Chemistry"
    ##  [7] "Included"              "HadGEM2-ES"           
    ##  [9] "40"                    "z"                    
    ## [11] "Included"              "Included"             
    ##  [1] "HadGAM2"               "N96"                  
    ##  [3] "38"                    "13"                   
    ##  [5] "Interactive"           "Atmospheric Chemistry"
    ##  [7] "Included"              "HadGEM2-ES"           
    ##  [9] "40"                    "z"                    
    ## [11] "Included"              "Included"             
    ##  [1] "HadGAM2"               "N96"                  
    ##  [3] "38"                    "13"                   
    ##  [5] "Interactive"           "Atmospheric Chemistry"
    ##  [7] "Included"              "HadGEM2-ES"           
    ##  [9] "40"                    "z"                    
    ## [11] "Included"              "Included"             
    ##  [1] "HadGAM2"               "N96"                  
    ##  [3] "38"                    "13"                   
    ##  [5] "Interactive"           "Atmospheric Chemistry"
    ##  [7] "Included"              "HadGEM2-ES"           
    ##  [9] "40"                    "z"                    
    ## [11] "Included"              "Included"             
    ##  [1] "INM-CM4"         "2x1.5deg"        "21"             
    ##  [4] "10"              "Prescribed"      "Not implemented"
    ##  [7] "Included"        "INM-CM4"         "40"             
    ## [10] "sigma"           "Included"        "Included"       
    ##  [1] "LMDZ5"                     "96x951.9x3.75degLMDZ96x95"
    ##  [3] "39"                        "0.04"                     
    ##  [5] "Semi-interactive"          "Not implemented"          
    ##  [7] "Included"                  "IPSL-CM5A-LR"             
    ##  [9] "31"                        "z"                        
    ## [11] "PISCES"                    "LIM2"                     
    ##  [1] "LMDZ5"                     "96x951.9x3.75degLMDZ96x95"
    ##  [3] "39"                        "0.04"                     
    ##  [5] "Semi-interactive"          "Not implemented"          
    ##  [7] "Included"                  "IPSL-CM5A-LR"             
    ##  [9] "31"                        "z"                        
    ## [11] "PISCES"                    "LIM2"                     
    ##  [1] "LMDZ5"                     "96x951.9x3.75degLMDZ96x95"
    ##  [3] "39"                        "0.04"                     
    ##  [5] "Semi-interactive"          "Not implemented"          
    ##  [7] "Included"                  "IPSL-CM5A-LR"             
    ##  [9] "31"                        "z"                        
    ## [11] "PISCES"                    "LIM2"                     
    ##  [1] "LMDZ5"                     "96x951.9x3.75degLMDZ96x95"
    ##  [3] "39"                        "0.04"                     
    ##  [5] "Semi-interactive"          "Not implemented"          
    ##  [7] "Included"                  "IPSL-CM5A-LR"             
    ##  [9] "31"                        "z"                        
    ## [11] "PISCES"                    "LIM2"                     
    ##  [1] "LMDZ5"                         "144x1431.25x2.5degLMDZ144x143"
    ##  [3] "39"                            "0.04"                         
    ##  [5] "Semi-interactive"              "Not implemented"              
    ##  [7] "Included"                      "IPSL-CM5A-MR"                 
    ##  [9] "31"                            "z"                            
    ## [11] "PISCES"                        "Included"                     
    ##  [1] "LMDZ5"                     "96x951.9x3.75degLMDZ96x95"
    ##  [3] "39"                        "0.04"                     
    ##  [5] "Semi-interactive"          "Not implemented"          
    ##  [7] "Included"                  "IPSL-CM5B-LR"             
    ##  [9] "31"                        "z"                        
    ## [11] "PISCES"                    "Included"                 
    ##  [1] "CCSR/NIES/ FRCGC AGCM6" "1.40625x1.40625degT85" 
    ##  [3] "40"                     "2.9"                   
    ##  [5] "SPRINTARS"              "Not implemented"       
    ##  [7] "MATSIRO"                "COCO4.5"               
    ##  [9] "50"                     "z-s"                   
    ## [11] "Not implemented"        "Included"              
    ##  [1] "CCSR/NIES/ FRCGC AGCM6" "1.40625x1.40625degT85" 
    ##  [3] "40"                     "2.9"                   
    ##  [5] "SPRINTARS"              "Not implemented"       
    ##  [7] "MATSIRO"                "COCO4.5"               
    ##  [9] "50"                     "z-s"                   
    ## [11] "Not implemented"        "Included"              
    ##  [1] "CCSR/NIES/ FRCGC AGCM6" "1.40625x1.40625degT85" 
    ##  [3] "40"                     "2.9"                   
    ##  [5] "SPRINTARS"              "Not implemented"       
    ##  [7] "MATSIRO"                "COCO4.5"               
    ##  [9] "50"                     "z-s"                   
    ## [11] "Not implemented"        "Included"              
    ##  [1] "MIROC-AGCM"      "T42"             "80"             
    ##  [4] "0.003"           "SPRINTARS"       "Not implemented"
    ##  [7] "MATSIRO"         "COCO3.4"         "44"             
    ## [10] "z-s"             "NPZD-type"       "Included"       
    ##  [1] "MIROC-AGCM" "T42"        "80"         "0.003"      "SPRINTARS" 
    ##  [6] "CHASER"     "MATSIRO"    "COCO3.4"    "44"         "z-s"       
    ## [11] "NPZD-type"  "Included"  
    ##  [1] "ECHAM6"          "T63"             "47"             
    ##  [4] "0.01"            "Prescribed"      "Not implemented"
    ##  [7] "JSBACH"          "MPIOM"           "40"             
    ## [10] "z"               "HAMOCC"          "Included"       
    ##  [1] "ECHAM6"          "T63"             "47"             
    ##  [4] "0.01"            "Prescribed"      "Not implemented"
    ##  [7] "JSBACH"          "MPIOM"           "40"             
    ## [10] "z"               "HAMOCC"          "Included"       
    ##  [1] "ECHAM6"          "T63"             "47"             
    ##  [4] "0.01"            "Prescribed"      "Not implemented"
    ##  [7] "JSBACH"          "MPIOM"           "40"             
    ## [10] "z"               "HAMOCC"          "Included"       
    ##  [1] "ECHAM6"          "T63"             "95"             
    ##  [4] "0.01"            "Prescribed"      "Not implemented"
    ##  [7] "JSBACH"          "MPIOM"           "40"             
    ## [10] "z"               "HAMOCC"          "Included"       
    ##  [1] "ECHAM6"          "T63"             "95"             
    ##  [4] "0.01"            "Prescribed"      "Not implemented"
    ##  [7] "JSBACH"          "MPIOM"           "40"             
    ## [10] "z"               "HAMOCC"          "Included"       
    ##  [1] "ECHAM6"          "T63"             "95"             
    ##  [4] "0.01"            "Prescribed"      "Not implemented"
    ##  [7] "JSBACH"          "MPIOM"           "40"             
    ## [10] "z"               "HAMOCC"          "Included"       
    ##  [1] "MRI-AGCM3.3"         "T159"                "48"                 
    ##  [4] "0.01"                "MASINGAR mk-2"       "Not implemented"    
    ##  [7] "HAL"                 "MRI.COM3"            "51"                 
    ## [10] "sigma"               "Not implemented"     "Included (MRI.COM3)"
    ##  [1] "CAM4-Oslo"       "1.9deg2.5deg"    "26"             
    ##  [4] "2.194067"        "CAM4-Oslo"       "CAM4-Oslo"      
    ##  [7] "CLM4"            "NorESM-Ocean"    "53"             
    ## [10] "isopycnic"       "Not implemented" "CICE4"          
    ##  [1] "CAM4-Oslo"    "1.9deg2.5deg" "26"           "2.194067"    
    ##  [5] "CAM4-Oslo"    "CAM4-Oslo"    "CLM4"         "NorESM-Ocean"
    ##  [9] "53"           "isopycnic"    "HAMOCC5"      "CICE4"

Set the column names

``` r
X <- as.data.frame(X)
colnames(X) <- c('atm.model','atm.grid',
                 'atm.lev','atm.top',
                 'atm.type','Atm.chem',
                 'Land','oce.model',
                 'oce.lev','oce.z',
                 'oce.bio','Seaice')
X$atm.lev <- as.integer(X$atm.lev)
X$atm.top <- as.numeric(X$atm.top)
X$oce.lev <- as.integer(X$oce.lev)
```

Combine the simulated temperature change with corresponding information about the model configuration

``` r
dT <- colMeans(window(zoo(global.t2m.gcm[[1]]),start=2070,end=2099))-
      colMeans(window(zoo(global.t2m.gcm[[1]]),start=1970,end=1999))
X <- as.data.frame(cbind(round(dT,2),X))
```

**Results**
-----------

Is the simulated temperature change sensitive to how the atmospheric model is configured?

``` r
fit.atm <- lm(dT ~ atm.grid + atm.lev + atm.top + oce.lev + oce.z,data=X)
print(summary(fit.atm))
```

    ## 
    ## Call:
    ## lm(formula = dT ~ atm.grid + atm.lev + atm.top + oce.lev + oce.z, 
    ##     data = X)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.52714 -0.13088 -0.00847  0.05716  0.41453 
    ## 
    ## Coefficients: (1 not defined because of singularities)
    ##                                       Estimate Std. Error t value Pr(>|t|)
    ## (Intercept)                            4.99886    0.68574   7.290 1.85e-10
    ## atm.grid1.40625x1.40625degT85         -1.47991    0.29691  -4.984 3.47e-06
    ## atm.grid144x1431.25x2.5degLMDZ144x143 -0.11496    0.45692  -0.252 0.801994
    ## atm.grid1.9deg2.5deg                  -0.26562    0.26824  -0.990 0.325007
    ## atm.grid1deg                          -0.87050    0.52287  -1.665 0.099804
    ## atm.grid200kmC48L48                    0.97442    0.42606   2.287 0.024802
    ## atm.grid2.5deg2degM45L24              -0.57825    0.26874  -2.152 0.034396
    ## atm.grid2.8125x2.8125deg              -2.68479    0.50132  -5.355 7.80e-07
    ## atm.grid2degx2.5deg                   -0.63709    0.50105  -1.272 0.207188
    ## atm.grid2x1.5deg                      -1.04733    0.58120  -1.802 0.075264
    ## atm.grid96x951.9x3.75degLMDZ96x95     -0.28997    0.41437  -0.700 0.486071
    ## atm.gridN96                            0.38798    0.28482   1.362 0.176915
    ## atm.gridT106                          -0.57904    0.32513  -1.781 0.078672
    ## atm.gridT127                          -0.27618    0.32501  -0.850 0.397970
    ## atm.gridT159                          -0.09278    0.35873  -0.259 0.796582
    ## atm.gridT42                           -0.91035    0.26141  -3.482 0.000804
    ## atm.gridT63                           -0.09327    0.31026  -0.301 0.764467
    ## atm.lev                               -0.06627    0.01987  -3.336 0.001287
    ## atm.top                               -0.04503    0.02062  -2.184 0.031866
    ## oce.lev                               -0.05633    0.04129  -1.364 0.176239
    ## oce.zisopycnic                        -1.59526    0.28861  -5.527 3.85e-07
    ## oce.zsigma                            -1.86663    0.36669  -5.090 2.28e-06
    ## oce.zz                                -1.53341    0.20845  -7.356 1.37e-10
    ## oce.zz*                               -1.69053    0.30159  -5.605 2.79e-07
    ## oce.zz-s                                    NA         NA      NA       NA
    ##                                          
    ## (Intercept)                           ***
    ## atm.grid1.40625x1.40625degT85         ***
    ## atm.grid144x1431.25x2.5degLMDZ144x143    
    ## atm.grid1.9deg2.5deg                     
    ## atm.grid1deg                          .  
    ## atm.grid200kmC48L48                   *  
    ## atm.grid2.5deg2degM45L24              *  
    ## atm.grid2.8125x2.8125deg              ***
    ## atm.grid2degx2.5deg                      
    ## atm.grid2x1.5deg                      .  
    ## atm.grid96x951.9x3.75degLMDZ96x95        
    ## atm.gridN96                              
    ## atm.gridT106                          .  
    ## atm.gridT127                             
    ## atm.gridT159                             
    ## atm.gridT42                           ***
    ## atm.gridT63                              
    ## atm.lev                               ** 
    ## atm.top                               *  
    ## oce.lev                                  
    ## oce.zisopycnic                        ***
    ## oce.zsigma                            ***
    ## oce.zz                                ***
    ## oce.zz*                               ***
    ## oce.zz-s                                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.2153 on 81 degrees of freedom
    ##   (3 observations deleted due to missingness)
    ## Multiple R-squared:  0.7967, Adjusted R-squared:  0.739 
    ## F-statistic:  13.8 on 23 and 81 DF,  p-value: < 2.2e-16
